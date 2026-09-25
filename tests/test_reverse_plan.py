"""POST /api/v1/stoichiometry/plan 的 HTTP API 测试（反向投料规划）。

测试反应均选用其它测试文件未提交过的方程式（NH3 合成、NH3 催化氧化、
Fe/I 氧化还原、Fe2O3 还原、KClO3 分解等），避免幂等缓存造成跨文件的
顺序依赖。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/stoichiometry/plan"


# ------------------------------------------------------------ 成功路径
def test_single_target_reverse_planning():
    # 3H2 + N2 -> 2NH3；目标 NH3=4 mol -> ξ=2，需 H2 6 mol、N2 2 mol
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "4"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"] == {"value": "2", "unit": "mol"}
    assert data["driving_targets"] == ["NH3"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["N2"]["coefficient"] == 1
    assert by_name["N2"]["required"] == "2"
    assert by_name["H2"]["coefficient"] == 3
    assert by_name["H2"]["required"] == "6"

    assert len(data["targets"]) == 1
    target = data["targets"][0]
    assert target["normalized"] == "NH3"
    assert target["target"] == "4"
    assert target["achieved"] == "4"
    assert target["excess"] == "0"
    assert target["is_driving"] is True

    products = {p["normalized"]: p for p in data["products"]}
    assert products["NH3"]["coefficient"] == 2
    assert products["NH3"]["achieved"] == "4"
    # 所有十进制结果均以字符串返回
    assert isinstance(by_name["H2"]["required"], str)


def test_multiple_targets_pick_max_ratio_with_excess():
    # 4NH3 + 5O2 -> 4NO + 6H2O
    # 目标 NO=1（比值 0.25）、H2O=6（比值 1）-> ξ=1，由 H2O 驱动；NO 超额 3
    resp = client.post(
        PATH,
        json={
            "reactants": ["NH3", "O2"],
            "products": ["NO", "H2O"],
            "targets": [
                {"species": "NO", "moles": "1"},
                {"species": "H2O", "moles": "6"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["H2O"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["NH3"]["required"] == "4"
    assert by_name["O2"]["required"] == "5"

    targets = {t["normalized"]: t for t in data["targets"]}
    assert targets["H2O"] == {
        "input": "H2O",
        "normalized": "H2O",
        "coefficient": 6,
        "target": "6",
        "achieved": "6",
        "excess": "0",
        "is_driving": True,
    }
    assert targets["NO"]["target"] == "1"
    assert targets["NO"]["achieved"] == "4"
    assert targets["NO"]["excess"] == "3"
    assert targets["NO"]["is_driving"] is False

    # 全部生成物都给达成量（不只是提交目标的那些）
    products = {p["normalized"]: p["achieved"] for p in data["products"]}
    assert products == {"NO": "4", "H2O": "6"}


def test_tied_driving_targets_when_targets_match_coefficient_ratio():
    # 4NH3 + 5O2 -> 4NO + 6H2O；目标恰好按系数比 4:6 -> 两者并列驱动，无超额
    resp = client.post(
        PATH,
        json={
            "reactants": ["NH3", "O2"],
            "products": ["NO", "H2O"],
            "targets": [
                {"species": "NO", "moles": "4"},
                {"species": "H2O", "moles": "6"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["H2O", "NO"]
    assert all(t["excess"] == "0" for t in data["targets"])
    assert all(t["is_driving"] for t in data["targets"])


def test_subset_of_products_may_be_targeted():
    # Fe2O3 + 3CO -> 2Fe + 3CO2；只把 Fe 列为目标（CO2 未提交）
    resp = client.post(
        PATH,
        json={
            "reactants": ["Fe2O3", "CO"],
            "products": ["Fe", "CO2"],
            "targets": [{"species": "Fe", "moles": "2"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["Fe"]
    assert [t["normalized"] for t in data["targets"]] == ["Fe"]
    products = {p["normalized"]: p["achieved"] for p in data["products"]}
    # 未作为目标的 CO2 仍给出在该进度下的达成量
    assert products == {"CO2": "3", "Fe": "2"}
    by_name = {r["normalized"]: r["required"] for r in data["reactants"]}
    assert by_name == {"CO": "3", "Fe2O3": "1"}


def test_non_terminating_decimal_is_rounded_to_28_digits():
    # Fe2O3 + 3CO -> 2Fe + 3CO2；目标 CO2=1 -> ξ=1/3（无限小数）
    resp = client.post(
        PATH,
        json={
            "reactants": ["Fe2O3", "CO"],
            "products": ["Fe", "CO2"],
            "targets": [{"species": "CO2", "moles": "1"}],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == (
        "0.3333333333333333333333333333"  # 1/3，28 位有效数字
    )
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["Fe2O3"]["required"] == "0.3333333333333333333333333333"
    assert by_name["CO"]["required"] == "1"
    target = data["targets"][0]
    assert target["achieved"] == "1"
    assert target["excess"] == "0"
    products = {p["normalized"]: p["achieved"] for p in data["products"]}
    assert products["Fe"] == "0.6666666666666666666666666667"
    assert products["CO2"] == "1"


def test_decimal_inputs_avoid_binary_float_error():
    # C2H4 + 3O2 -> 2CO2 + 2H2O；目标 CO2=0.1 -> ξ=0.05，O2 需 0.15
    resp = client.post(
        PATH,
        json={
            "reactants": ["C2H4", "O2"],
            "products": ["CO2", "H2O"],
            "targets": [{"species": "CO2", "moles": "0.1"}],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.05"
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["O2"]["required"] == "0.15"
    assert by_name["C2H4"]["required"] == "0.05"
    products = {p["normalized"]: p["achieved"] for p in data["products"]}
    assert products == {"CO2": "0.1", "H2O": "0.1"}


def test_zero_target_gives_zero_extent_and_requirements():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [
                {"species": "NH3", "moles": "0"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0"
    assert data["driving_targets"] == ["NH3"]
    assert all(r["required"] == "0" for r in data["reactants"])
    assert data["targets"][0]["achieved"] == "0"
    assert data["targets"][0]["excess"] == "0"
    assert all(p["achieved"] == "0" for p in data["products"])


def test_scientific_notation_and_underscores_accepted():
    # C2H4 + 3O2 -> 2CO2 + 2H2O；目标 CO2=2e-1=0.2 -> ξ=0.1
    resp = client.post(
        PATH,
        json={
            "reactants": ["C2H4", "O2"],
            "products": ["CO2", "H2O"],
            "targets": [
                {"species": "CO2", "moles": "2e-1"},
                {"species": "H2O", "moles": "1_000"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    # ξ = max(0.2/2, 1000/2) = 500，由 H2O 驱动
    assert data["reaction_extent"]["value"] == "500"
    targets = {t["normalized"]: t for t in data["targets"]}
    assert targets["CO2"]["target"] == "0.2"
    assert targets["CO2"]["achieved"] == "1000"
    assert targets["CO2"]["excess"] == "999.8"


def test_charged_species_alias_notation_recognized():
    # 2Fe^3+ + 2I^- -> 2Fe^2+ + I2；目标用等价写法 Fe2+ 给出
    resp = client.post(
        PATH,
        json={
            "reactants": ["Fe3+", "I-"],
            "products": ["Fe2+", "I2"],
            "targets": [{"species": "Fe^2+", "moles": "0.2"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["normalized"]["equation"] == "2Fe^3+ + 2I^- -> 2Fe^2+ + I2"
    assert data["reaction_extent"]["value"] == "0.1"
    target = data["targets"][0]
    assert target["input"] == "Fe^2+"
    assert target["normalized"] == "Fe^2+"
    assert target["achieved"] == "0.2"
    by_name = {r["normalized"]: r["required"] for r in data["reactants"]}
    assert by_name == {"Fe^3+": "0.2", "I^-": "0.2"}
    products = {p["normalized"]: p["achieved"] for p in data["products"]}
    assert products == {"Fe^2+": "0.2", "I2": "0.1"}


# ------------------------------------------------------------ 顺序/写法无关性
def test_reaction_id_and_results_stable_across_order_and_notation():
    payload_a = {
        "reactants": ["N2", "H2"],
        "products": ["NH3"],
        "targets": [{"species": "NH3", "moles": "4"}],
    }
    payload_b = {
        "reactants": ["H2", "N2"],
        "products": ["N H3"],
        "targets": [{"species": "N H3", "moles": "4.0"}],
    }
    a = client.post(PATH, json=payload_a).json()
    b = client.post(PATH, json=payload_b).json()
    assert a["reaction_id"] == b["reaction_id"]
    assert a["reaction_extent"] == b["reaction_extent"]
    assert a["driving_targets"] == b["driving_targets"]

    def computed(blob):
        return (
            [(r["normalized"], r["coefficient"], r["required"])
             for r in blob["reactants"]],
            [(t["normalized"], t["coefficient"], t["target"], t["achieved"],
              t["excess"], t["is_driving"]) for t in blob["targets"]],
            [(p["normalized"], p["coefficient"], p["achieved"])
             for p in blob["products"]],
        )

    assert computed(a) == computed(b)
    # input 字段保留各自原始写法（仅去首尾空白，内部空白照原样回显）
    assert b["targets"][0]["input"] == "N H3"
    assert b["products"][0]["input"] == "N H3"


def test_reaction_id_matches_balance_endpoint():
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["N2", "H2"], "products": ["NH3"]},
    ).json()
    plan = client.post(
        PATH,
        json={
            "reactants": ["H2", "N2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "4"}],
        },
    ).json()
    assert balance["reaction_id"] == plan["reaction_id"]


def test_cached_flag_reflects_balance_cache():
    # 使用本文件内其它用例未提交过的反应，避免幂等缓存造成顺序依赖
    body = {
        "reactants": ["KClO3"],
        "products": ["KCl", "O2"],
        "targets": [{"species": "O2", "moles": "3"}],
    }
    first = client.post(PATH, json=body).json()
    second = client.post(PATH, json=body).json()
    assert first["cached"] is False
    assert second["cached"] is True
    assert first["reaction_id"] == second["reaction_id"]
    # 2KClO3 -> 2KCl + 3O2；目标 O2=3 -> ξ=1
    assert first["reaction_extent"]["value"] == "1"
    assert first["reactants"][0]["required"] == "2"


# ------------------------------------------------------------ 结构化错误
def test_unknown_target_species_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "CH4", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"
    assert body["error"]["details"]["unknown"] == ["CH4"]
    assert "reaction_id" in body


def test_reactant_given_as_target_is_unknown():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "N2", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["unknown"] == ["N2"]


def test_duplicate_target_entry_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [
                {"species": "NH3", "moles": "1"},
                {"species": "N H3", "moles": "4"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate"] == ["NH3"]


def test_invalid_moles_and_duplicate_aggregated_in_one_response():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [
                {"species": "NH3", "moles": "1"},
                {"species": "NH3", "moles": "-3"},
                {"species": "CH4", "moles": "2"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate"] == ["NH3"]
    assert details["unknown"] == ["CH4"]
    assert details["invalid_moles"] == [{"species": "NH3", "moles": "-3"}]


@pytest.mark.parametrize(
    "bad",
    ["-1", "-0.0", "abc", "1.2.3", "NaN", "Infinity", "", "1e10001", "1e-10001"],
)
def test_invalid_moles_values(bad):
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": bad}],
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
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": 0.1}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unparseable_species_in_targets_is_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "Xx", "moles": "1"}],
        },
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"] == {"field": "targets.species", "input": "Xx"}


def test_non_unique_balance_reuses_structured_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["C", "O2"],
            "products": ["CO", "CO2"],
            "targets": [{"species": "CO", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NON_UNIQUE_BALANCE"


def test_impossible_balance_reuses_structured_error():
    # 选用其它测试未用过的不可能反应，避免顺序依赖
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2"],
            "products": ["Au"],
            "targets": [{"species": "Au", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_reaction_side_duplicate_reuses_conflict():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "N2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "1"}],
        },
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_reaction_formula_parse_error_reuses_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["Zz", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "1"}],
        },
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


@pytest.mark.parametrize(
    "payload",
    [
        {"reactants": ["N2", "H2"], "products": ["NH3"]},  # 缺 targets
        {"reactants": ["N2", "H2"], "products": ["NH3"], "targets": []},
        {"reactants": [], "products": ["NH3"],
         "targets": [{"species": "NH3", "moles": "1"}]},
    ],
)
def test_request_shape_validation(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


# ------------------------------------------------------------ 响应形状与 OpenAPI
def test_success_payload_shape():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "4"}],
        },
    )
    data = resp.json()
    assert set(data) == {
        "reaction_id", "cached", "normalized", "reaction_extent",
        "driving_targets", "reactants", "targets", "products",
    }
    assert set(data["reactants"][0]) == {
        "input", "normalized", "coefficient", "required",
    }
    assert set(data["targets"][0]) == {
        "input", "normalized", "coefficient",
        "target", "achieved", "excess", "is_driving",
    }
    assert set(data["products"][0]) == {
        "input", "normalized", "coefficient", "achieved",
    }


def test_openapi_documents_plan_endpoint():
    schema = client.get("/openapi.json").json()
    path = schema["paths"]["/api/v1/stoichiometry/plan"]
    assert "post" in path
    assert path["post"]["tags"] == ["计量"]
    assert "ReversePlanRequest" in path["post"]["requestBody"]["content"][
        "application/json"
    ]["schema"]["$ref"]
    assert "ReversePlanResponse" in path["post"]["responses"]["200"][
        "content"
    ]["application/json"]["schema"]["$ref"]
    component_names = schema["components"]["schemas"]
    for name in (
        "ReversePlanRequest",
        "ReversePlanResponse",
        "ReversePlanTarget",
        "PlanReactantResult",
        "PlanTargetResult",
        "PlanProductResult",
    ):
        assert name in component_names


def test_batch_endpoint_payload_unchanged():
    # 回归保护：/stoichiometry/batch 的字段不受新接口影响
    resp = client.post(
        "/api/v1/stoichiometry/batch",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "amounts": [
                {"species": "H2", "moles": "2"},
                {"species": "O2", "moles": "2"},
            ],
        },
    )
    data = resp.json()
    assert set(data) == {
        "reaction_id", "cached", "normalized", "reaction_extent",
        "limiting_reactants", "reactants", "products",
    }
    assert data["reaction_extent"]["value"] == "1"
