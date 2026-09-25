"""HTTP API 测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_balance_success():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2", "O2"], "products": ["H2O"]}
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["balanced"] is True
    vector = data["coefficients"]["vector"]
    assert vector == [2, 1, 2]
    assert data["normalized"]["equation"] == "2H2 + O2 -> 2H2O"
    assert data["charge_balance"]["balanced"] is True
    assert data["cached"] is False


def test_idempotent_same_normalized_reaction():
    payload_a = {"reactants": ["Ca(OH)2", "HCl"], "products": ["CaCl2", "H2O"]}
    # 等价写法：方括号替代圆括号、内部空白不同
    payload_b = {
        "reactants": ["Ca [OH]2", " HCl"],
        "products": ["CaCl2", "H2 O"],
    }

    r1 = client.post("/api/v1/balance", json=payload_a)
    r2 = client.post("/api/v1/balance", json=payload_b)
    assert r1.status_code == 200 and r2.status_code == 200
    d1, d2 = r1.json(), r2.json()
    assert d1["reaction_id"] == d2["reaction_id"]
    assert d2["cached"] is True
    assert d1["coefficients"]["vector"] == d2["coefficients"]["vector"]


def test_idempotent_order_independent():
    r1 = client.post(
        "/api/v1/balance",
        json={"reactants": ["O2", "H2"], "products": ["H2O"]},
    )
    r2 = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    assert r1.json()["reaction_id"] == r2.json()["reaction_id"]


def test_get_reaction_by_id():
    created = client.post(
        "/api/v1/balance", json={"reactants": ["Fe", "O2"], "products": ["Fe2O3"]}
    ).json()
    rid = created["reaction_id"]
    resp = client.get(f"/api/v1/reactions/{rid}")
    assert resp.status_code == 200
    assert resp.json()["reaction_id"] == rid
    assert resp.json()["cached"] is True


def test_get_reaction_not_found():
    resp = client.get("/api/v1/reactions/rxn_doesnotexist")
    assert resp.status_code == 404
    assert resp.json()["error"]["code"] == "NOT_FOUND"


def test_balance_formula_parse_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["Xx"], "products": ["H2O"]}
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert "reaction_id" in body


def test_balance_non_unique_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]}
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "NON_UNIQUE_BALANCE"
    assert body["error"]["details"]["nullity"] == 2


def test_balance_impossible_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2O"], "products": ["Au"]}
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_duplicate_species_conflict():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2", "H2"], "products": ["H2O"]}
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_invalid_request_body():
    resp = client.post(
        "/api/v1/balance", json={"reactants": [], "products": ["H2O"]}
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    resp = client.post("/api/v1/balance", json={"reactants": ["H2"]})
    assert resp.status_code == 422


def test_analyze_shows_raw_and_balanced():
    resp = client.post(
        "/api/v1/analyze", json={"reactants": ["H2", "O2"], "products": ["H2O"]}
    )
    assert resp.status_code == 200
    data = resp.json()
    raw_rows = {row["element"]: row for row in data["raw"]["elements"]}
    assert raw_rows["H"]["reactant_total"] == 2
    assert raw_rows["H"]["product_total"] == 2
    assert raw_rows["H"]["difference"] == 0
    assert raw_rows["O"]["reactant_total"] == 2
    assert raw_rows["O"]["product_total"] == 1
    assert raw_rows["O"]["difference"] == 1
    assert data["raw"]["coefficients_assumed"] == 1
    assert data["balanced_ok"] is True
    balanced_rows = {row["element"]: row for row in data["balanced"]["elements"]}
    assert balanced_rows["O"]["reactant_total"] == 2
    assert balanced_rows["O"]["product_total"] == 2
    assert all(row["balanced"] for row in balanced_rows.values())


def test_analyze_non_unique_returns_error_in_body():
    resp = client.post(
        "/api/v1/analyze",
        json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["balanced_ok"] is False
    assert data["balanced"] is None
    assert data["error"]["code"] == "NON_UNIQUE_BALANCE"
    # 原始视图仍然给出差异
    assert {row["element"] for row in data["raw"]["elements"]} == {"C", "O"}


def test_analyze_charge_rows():
    resp = client.post(
        "/api/v1/analyze",
        json={"reactants": ["MnO4-", "H+", "e-"], "products": ["Mn2+", "H2O"]},
    )
    data = resp.json()
    assert data["balanced"]["charge"] == {
        "reactant_total": 2,
        "product_total": 2,
        "difference": 0,
        "balanced": True,
    }


def test_redox_full_reaction():
    resp = client.post(
        "/api/v1/balance",
        json={
            "reactants": ["Cu", "NO3-", "H+"],
            "products": ["Cu2+", "NO", "H2O"],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    coef = {item["normalized"]: item["coefficient"]
            for item in data["coefficients"]["reactants"]
            + data["coefficients"]["products"]}
    assert coef == {"Cu": 3, "NO3^-": 2, "H^+": 8, "Cu^2+": 3, "NO": 2, "H2O": 4}


# ------------------------------------------------------------ /balance/variants
def test_variants_unique_reaction_single_solution():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "max_coefficient": 30,
            "limit": 100,
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["nullity"] == 1
    assert data["total"] == data["count"] == 1
    assert data["truncated"] is False
    item = data["variants"][0]
    assert item["coefficients"]["vector"] == [2, 1, 2]
    assert item["sum"] == 5
    assert item["equation"] == "2H2 + O2 -> 2H2O"
    # 归一化块不携带系数，保持物种清单语义
    assert data["normalized"]["equation"] == "H2 + O2 -> H2O"


def test_variants_underdetermined_sorted_and_truncated():
    body = {
        "reactants": ["C", "O2"],
        "products": ["CO", "CO2"],
        "max_coefficient": 10,
        "limit": 5,
    }
    resp = client.post("/api/v1/balance/variants", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["nullity"] == 2
    assert data["total"] == 15
    assert data["count"] == 5
    assert data["truncated"] is True
    vectors = [v["coefficients"]["vector"] for v in data["variants"]]
    assert vectors == [
        [3, 2, 2, 1],
        [4, 3, 2, 2],
        [5, 3, 4, 1],
        [5, 4, 2, 3],
        [6, 5, 2, 4],
    ]
    sums = [v["sum"] for v in data["variants"]]
    assert sums == sorted(sums)


def test_variants_limit_not_exceeded_sets_truncated_false():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["C", "O2"],
            "products": ["CO", "CO2"],
            "max_coefficient": 4,
            "limit": 100,
        },
    )
    data = resp.json()
    assert data["truncated"] is False
    assert data["count"] == data["total"] == 2
    assert [v["coefficients"]["vector"] for v in data["variants"]] == [
        [3, 2, 2, 1],
        [4, 3, 2, 2],
    ]


def test_variants_no_positive_solution_empty_list():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["H2", "O2", "N2"],
            "products": ["H2O"],
            "max_coefficient": 10,
            "limit": 10,
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["variants"] == []
    assert data["count"] == data["total"] == 0
    assert data["truncated"] is False


def test_variants_space_too_large():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["H2", "O2", "H2O", "H2O2", "O3"],
            "products": ["HO2"],
            "max_coefficient": 10,
            "limit": 10,
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "VARIANT_SPACE_TOO_LARGE"
    assert body["error"]["details"]["nullity"] == 4
    assert "reaction_id" in body


def test_variants_impossible_balance():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["H2O"],
            "products": ["Au"],
            "max_coefficient": 10,
            "limit": 10,
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_variants_parse_error():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["Xx"],
            "products": ["H2O"],
            "max_coefficient": 5,
            "limit": 10,
        },
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


def test_variants_duplicate_species():
    resp = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["H2", "H2"],
            "products": ["H2O"],
            "max_coefficient": 5,
            "limit": 10,
        },
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


@pytest.mark.parametrize("payload", [
    {"reactants": ["H2"], "products": ["H2O"], "max_coefficient": 1, "limit": 10},
    {"reactants": ["H2"], "products": ["H2O"], "max_coefficient": 31, "limit": 10},
    {"reactants": ["H2"], "products": ["H2O"], "max_coefficient": 10, "limit": 0},
    {"reactants": ["H2"], "products": ["H2O"], "max_coefficient": 10, "limit": 101},
    {"reactants": ["H2"], "products": ["H2O"], "max_coefficient": "x", "limit": 10},
    {"reactants": [], "products": ["H2O"], "max_coefficient": 10, "limit": 10},
])
def test_variants_invalid_params(payload):
    resp = client.post("/api/v1/balance/variants", json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_variants_reaction_id_matches_balance_endpoint():
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["Ca(OH)2", "HCl"], "products": ["CaCl2", "H2O"]},
    ).json()
    variants = client.post(
        "/api/v1/balance/variants",
        json={
            "reactants": ["Ca [OH]2", " HCl"],
            "products": ["CaCl2", "H2 O"],
            "max_coefficient": 20,
            "limit": 10,
        },
    ).json()
    assert balance["reaction_id"] == variants["reaction_id"]


def test_balance_endpoint_unchanged_payload_shape():
    # 回归保护：/balance 的字段与缓存语义不受新接口影响
    resp = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    data = resp.json()
    assert set(data) == {
        "reaction_id", "cached", "normalized", "coefficients",
        "balanced", "elements", "charge_balance",
    }
    assert data["normalized"]["equation"] == "2H2 + O2 -> 2H2O"
