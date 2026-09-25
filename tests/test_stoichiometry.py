"""POST /api/v1/stoichiometry/batch 的 HTTP API 测试。

测试反应均选用其它测试文件未提交过的方程式（CH4 燃烧、Mg(OH)2 中和、
Al 氧化等），避免幂等缓存造成跨文件的顺序依赖。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/stoichiometry/batch"


# ------------------------------------------------------------ 成功路径
def test_basic_limiting_reactant_and_yields():
    # CH4 + 2O2 -> CO2 + 2H2O；1 mol CH4 + 3 mol O2 -> ξ=1，CH4 限量
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "O2", "moles": "3"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"] == {"value": "1", "unit": "mol"}
    assert data["limiting_reactants"] == ["CH4"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["CH4"]["fed"] == "1"
    assert by_name["CH4"]["consumed"] == "1"
    assert by_name["CH4"]["remaining"] == "0"
    assert by_name["CH4"]["is_limiting"] is True
    assert by_name["O2"]["fed"] == "3"
    assert by_name["O2"]["consumed"] == "2"
    assert by_name["O2"]["remaining"] == "1"
    assert by_name["O2"]["is_limiting"] is False

    products = {p["normalized"]: p for p in data["products"]}
    assert products["CO2"]["coefficient"] == 1
    assert products["CO2"]["theoretical_yield"] == "1"
    assert products["H2O"]["coefficient"] == 2
    assert products["H2O"]["theoretical_yield"] == "2"
    # 所有十进制结果均以字符串返回
    assert isinstance(data["reaction_extent"]["value"], str)


def test_tied_limiting_reactants_when_fed_in_stoichiometric_ratio():
    # CH4 + 2O2 -> ...；1 mol CH4 与 2 mol O2 恰好按计量比 -> 两者并列限量
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "2.5"},
                {"species": "O2", "moles": "5"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["limiting_reactants"] == ["CH4", "O2"]
    assert data["reaction_extent"]["value"] == "2.5"
    assert all(r["remaining"] == "0" for r in data["reactants"])
    products = {p["normalized"]: p for p in data["products"]}
    assert products["CO2"]["theoretical_yield"] == "2.5"
    assert products["H2O"]["theoretical_yield"] == "5"


def test_repeating_decimal_is_string_without_float_artifacts():
    # 4Al + 3O2 -> 2Al2O3；各投 1 mol -> Al 限量，ξ=1/4（有限小数）
    resp = client.post(
        PATH,
        json={
            "reactants": ["Al", "O2"],
            "products": ["Al2O3"],
            "amounts": [
                {"species": "Al", "moles": "1"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.25"
    o2 = next(r for r in data["reactants"] if r["normalized"] == "O2")
    assert o2["consumed"] == "0.75"
    assert o2["remaining"] == "0.25"
    assert data["products"][0]["theoretical_yield"] == "0.5"

    # 无限小数：2C2H6 + 7O2 -> 4CO2 + 6H2O；各 1 mol -> O2 限量，ξ=1/7
    resp2 = client.post(
        PATH,
        json={
            "reactants": ["C2H6", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "C2H6", "moles": "1"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    data2 = resp2.json()
    assert data2["reaction_extent"]["value"] == (
        "0.1428571428571428571428571429"  # 1/7，28 位有效数字
    )
    products = {p["normalized"]: p for p in data2["products"]}
    assert products["CO2"]["theoretical_yield"] == "0.5714285714285714285714285714"
    assert products["H2O"]["theoretical_yield"] == "0.8571428571428571428571428571"


def test_decimal_inputs_avoid_binary_float_error():
    # Mg(OH)2 + 2HCl -> MgCl2 + 2H2O；各 0.1 mol -> HCl 限量，ξ=0.05
    # 若经 float，0.1/2 会出现 0.05000…001 一类二进制尾巴
    resp = client.post(
        PATH,
        json={
            "reactants": ["Mg(OH)2", "HCl"],
            "products": ["MgCl2", "H2O"],
            "amounts": [
                {"species": "Mg(OH)2", "moles": "0.1"},
                {"species": "HCl", "moles": "0.1"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.05"
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["HCl"]["consumed"] == "0.1"
    assert by_name["HCl"]["remaining"] == "0"
    assert by_name["Mg(OH)2"]["consumed"] == "0.05"
    assert by_name["Mg(OH)2"]["remaining"] == "0.05"
    water = next(p for p in data["products"] if p["normalized"] == "H2O")
    assert water["theoretical_yield"] == "0.1"


def test_charged_species_three_way_tie_with_alias_notation():
    # MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
    # 投料 0.1 / 0.8 / 0.5 mol 恰好按系数比，三者并列限量；物种用等价写法给出
    resp = client.post(
        PATH,
        json={
            "reactants": ["MnO4-", "H+", "e-"],
            "products": ["Mn2+", "H2O"],
            "amounts": [
                {"species": "MnO4^-", "moles": "0.1"},
                {"species": "H^+", "moles": "0.8"},
                {"species": "e", "moles": "0.5"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.1"
    assert data["limiting_reactants"] == ["e^-", "H^+", "MnO4^-"]
    products = {p["normalized"]: p["theoretical_yield"] for p in data["products"]}
    assert products == {"Mn^2+": "0.1", "H2O": "0.4"}


def test_all_zero_feed_marks_every_reactant_limiting():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "0"},
                {"species": "O2", "moles": "0"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0"
    assert data["limiting_reactants"] == ["CH4", "O2"]
    assert all(p["theoretical_yield"] == "0" for p in data["products"])


def test_scientific_notation_and_underscores_accepted():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "2e-1"},
                {"species": "O2", "moles": "1_000"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["CH4"]["fed"] == "0.2"
    assert by_name["O2"]["fed"] == "1000"
    # ξ = 0.2/1 = 0.2（CH4 限量）
    assert data["reaction_extent"]["value"] == "0.2"


# ------------------------------------------------------------ 顺序/写法无关性
def test_reaction_id_and_results_stable_across_order_and_notation():
    payload_a = {
        "reactants": ["Mg(OH)2", "HCl"],
        "products": ["MgCl2", "H2O"],
        "amounts": [
            {"species": "Mg(OH)2", "moles": "0.3"},
            {"species": "HCl", "moles": "0.4"},
        ],
    }
    payload_b = {
        "reactants": ["HCl", "Mg [OH]2"],
        "products": ["H2 O", "MgCl2"],
        "amounts": [
            {"species": " HCl", "moles": "0.40"},
            {"species": "Mg(OH)2", "moles": "0.30"},
        ],
    }
    a = client.post(PATH, json=payload_a).json()
    b = client.post(PATH, json=payload_b).json()
    assert a["reaction_id"] == b["reaction_id"]
    assert a["reaction_extent"] == b["reaction_extent"]
    assert a["limiting_reactants"] == b["limiting_reactants"]

    def computed(blob):
        return (
            [
                (r["normalized"], r["coefficient"], r["fed"], r["consumed"],
                 r["remaining"], r["is_limiting"])
                for r in blob["reactants"]
            ],
            [
                (p["normalized"], p["coefficient"], p["theoretical_yield"])
                for p in blob["products"]
            ],
        )

    assert computed(a) == computed(b)


def test_reaction_id_matches_balance_endpoint():
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["CH4", "O2"], "products": ["CO2", "H2O"]},
    ).json()
    batch = client.post(
        PATH,
        json={
            "reactants": ["O2", "CH4"],
            "products": ["H2O", "CO2"],
            "amounts": [
                {"species": "O2", "moles": "5"},
                {"species": "CH4", "moles": "5"},
            ],
        },
    ).json()
    assert balance["reaction_id"] == batch["reaction_id"]


def test_cached_flag_reflects_balance_cache():
    # 使用本文件内其它用例未提交过的反应，避免幂等缓存造成顺序依赖
    body = {
        "reactants": ["P4", "O2"],
        "products": ["P2O5"],
        "amounts": [
            {"species": "P4", "moles": "1"},
            {"species": "O2", "moles": "5"},
        ],
    }
    first = client.post(PATH, json=body).json()
    second = client.post(PATH, json=body).json()
    assert first["cached"] is False
    assert second["cached"] is True


# ------------------------------------------------------------ 结构化错误
def test_missing_amount_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [{"species": "CH4", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"
    assert body["error"]["details"]["missing"] == ["O2"]
    assert "reaction_id" in body


def test_duplicate_amount_entry_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "CH4", "moles": "2"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate"] == ["CH4"]


def test_unknown_species_in_amounts_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "N2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["unknown"] == ["N2"]
    assert details["missing"] == ["O2"]


@pytest.mark.parametrize(
    "bad",
    ["-1", "-0.0", "abc", "1.2.3", "NaN", "Infinity", "", "1e10001", "1e-10001"],
)
def test_invalid_moles_values(bad):
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": bad},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    # 空字符串在 Pydantic 层即被拦截；其余在业务层聚合报告
    assert body["error"]["code"] == "INVALID_REQUEST"


def test_moles_must_be_string_not_json_number():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4"],
            "products": ["CO2"],
            "amounts": [{"species": "CH4", "moles": 0.1}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unparseable_species_in_amounts_is_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "Xx", "moles": "1"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["input"] == "Xx"


def test_non_unique_balance_reuses_structured_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["C", "O2"],
            "products": ["CO", "CO2"],
            "amounts": [
                {"species": "C", "moles": "1"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NON_UNIQUE_BALANCE"


def test_impossible_balance_reuses_structured_error():
    # 选用 test_api 未用过的不可能反应，避免顺序依赖
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4"],
            "products": ["Au"],
            "amounts": [{"species": "CH4", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_reaction_side_duplicate_reuses_conflict():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "CH4"],
            "products": ["CO2"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "CH4", "moles": "2"},
            ],
        },
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_reaction_formula_parse_error_reuses_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["Zz", "O2"],
            "products": ["H2O"],
            "amounts": [
                {"species": "Zz", "moles": "1"},
                {"species": "O2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


@pytest.mark.parametrize(
    "payload",
    [
        {"reactants": ["CH4"], "products": ["CO2"]},  # 缺 amounts
        {"reactants": [], "products": ["CO2"],
         "amounts": [{"species": "CH4", "moles": "1"}]},
        {"reactants": ["CH4"], "products": ["CO2"], "amounts": []},
    ],
)
def test_request_shape_validation(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"
