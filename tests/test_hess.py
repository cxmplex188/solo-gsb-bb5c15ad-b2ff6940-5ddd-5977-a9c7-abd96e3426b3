"""赫斯定律反应焓推导（``POST /api/v1/thermochemistry/hess``）测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


PATH = "/api/v1/thermochemistry/hess"


# ------------------------------------------------------------- 经典正向/反向
def test_basic_hess_reverse_usage():
    """C + 2H2 -> CH4：第 3 条燃烧反应需反向使用（倍数 -1）。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "H2"], "products": ["CH4"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
                {"reactants": ["CH4", "O2"], "products": ["CO2", "H2O"],
                 "enthalpy_kj": "-890.3"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["target"]["equation"] == "C + 2H2 -> CH4"
    multipliers = [k["multiplier"] for k in data["known_reactions"]]
    assert multipliers == ["1", "1", "-1"]
    # -393.5 + (-571.6) + 890.3 = -74.8
    assert data["derived_enthalpy_kj"] == "-74.8"
    # 中间物种 CO2 / H2O / O2 完全抵消
    net = {item["normalized"]: item for item in data["net_coefficients"]}
    assert net["C"]["coefficient"] == "-1" and net["C"]["side"] == "reactant"
    assert net["H2"]["coefficient"] == "-2" and net["H2"]["side"] == "reactant"
    assert net["CH4"]["coefficient"] == "1" and net["CH4"]["side"] == "product"
    for name in ("CO2", "H2O", "O2"):
        assert net[name]["coefficient"] == "0"
        assert net[name]["side"] == "cancelled"


def test_fractional_multipliers():
    """C + O2 -> CO2 由两条「双倍」反应各取 1/2 组合而成。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO"],
                 "enthalpy_kj": "-221"},
                {"reactants": ["CO", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-566"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert [k["multiplier"] for k in data["known_reactions"]] == ["1/2", "1/2"]
    assert [k["multiplier_used"] for k in data["known_reactions"]] == [
        "-110.5",
        "-283",
    ]
    assert data["derived_enthalpy_kj"] == "-393.5"
    net = {item["normalized"]: item["coefficient"] for item in data["net_coefficients"]}
    assert net == {"C": "-1", "CO": "0", "CO2": "1", "O2": "-1"}


def test_three_halves_ammonia_synthesis():
    """N2 + 3H2 -> 2NH3：倍数 (-1/2, 1, 3/2)，ΔH = -92.5。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["N2", "H2"], "products": ["NH3"]},
            "known_reactions": [
                {"reactants": ["NH3", "O2"], "products": ["NO", "H2O"],
                 "enthalpy_kj": "-905"},
                {"reactants": ["N2", "O2"], "products": ["NO"],
                 "enthalpy_kj": "181"},
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-484"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["target"]["equation"] == "3H2 + N2 -> 2NH3"
    assert [k["multiplier"] for k in data["known_reactions"]] == [
        "-1/2",
        "1",
        "3/2",
    ]
    assert data["derived_enthalpy_kj"] == "-92.5"


def test_zero_enthalpy_and_integer_contribution_strings():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["known_reactions"][0]["multiplier"] == "1"
    assert data["known_reactions"][0]["multiplier_used"] == "-393.5"


def test_equivalent_formulas_recognized_in_combination():
    """已知反应用等价写法（空白、电荷风格）不影响物种匹配。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "H2"], "products": ["CH4"]},
            "known_reactions": [
                {"reactants": [" C ", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
                {"reactants": ["[CH4]", "O2"], "products": ["C O2", "H2O"],
                 "enthalpy_kj": "-890.3"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["derived_enthalpy_kj"] == "-74.8"


# ------------------------------------------------------------- 不可推导
def test_underivable_target_only_species():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["N2", "H2"], "products": ["NH3"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-285.8"},
            ],
        },
    )
    assert resp.status_code == 422
    error = resp.json()["error"]
    assert error["code"] == "HESS_UNDERIVABLE"
    assert sorted(error["details"]["target_only_species"]) == ["N2", "NH3"]


def test_underivable_inconsistent_stoichiometry():
    """目标物种都出现，但 H2 只存在于第 2 条反应 -> 系数必为 0，无法复现目标。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["CO", "H2O"], "products": ["CO2", "H2"],
                 "enthalpy_kj": "-41.2"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "HESS_UNDERIVABLE"


def test_no_known_reactions_is_422():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [],
        },
    )
    assert resp.status_code == 422


def test_more_than_six_known_reactions_is_422():
    one = {"reactants": ["C", "O2"], "products": ["CO2"], "enthalpy_kj": "-1"}
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [dict(one) for _ in range(7)],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


# ------------------------------------------------------------- 歧义
def test_ambiguous_different_enthalpies():
    """两条同计量已知反应给出不同反应焓 -> 组合无穷多、焓不唯一 -> 歧义。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-400"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "HESS_AMBIGUOUS"


def test_redundant_combination_same_enthalpy_is_ok():
    """组合方式不唯一但各组合导出的焓相同 -> 正常返回（取自由变量为 0 的特解）。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert [k["multiplier"] for k in data["known_reactions"]] == ["1", "0"]
    assert data["derived_enthalpy_kj"] == "-393.5"


def test_cycle_closure_consistent_enthalpies_is_ok():
    """三条反应构成闭合循环且 ΔH 之和为 0：组合多解但焓一致，正常返回。"""
    # v1: C+O2->CO2(-393.5); v2: 2C+O2->2CO(-221); v3: 2CO+O2->2CO2(-566)
    # 齐次关系 v1 = v2/2 + v3/2，焓：-110.5 + -283 = -393.5 恰好自洽
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["C", "O2"], "products": ["CO"],
                 "enthalpy_kj": "-221"},
                {"reactants": ["CO", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-566"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["derived_enthalpy_kj"] == "-393.5"


# ------------------------------------------------------------- 输入校验
def test_invalid_enthalpy_string():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "abc"},
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-1"},
            ],
        },
    )
    assert resp.status_code == 422
    error = resp.json()["error"]
    assert error["code"] == "INVALID_REQUEST"
    assert error["details"]["invalid_enthalpies"] == [
        {"index": 0, "enthalpy_kj": "abc"}
    ]


def test_numeric_enthalpy_rejected():
    """JSON 数字而非十进制字符串 -> 422（统一由请求体校验拦截）。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": -393.5},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_nan_infinity_enthalpy_rejected():
    for bad in ("NaN", "Infinity", "-Infinity"):
        resp = client.post(
            PATH,
            json={
                "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
                "known_reactions": [
                    {"reactants": ["C", "O2"], "products": ["CO2"],
                     "enthalpy_kj": bad},
                ],
            },
        )
        assert resp.status_code == 422, bad


def test_scientific_notation_and_positive_enthalpy():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "1e2"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["derived_enthalpy_kj"] == "100"


def test_missing_target_field():
    resp = client.post(
        PATH,
        json={
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-1"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_missing_enthalpy_field():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"]},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


# ------------------------------------------------------------- 复用配平错误
def test_target_formula_parse_error_is_400():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["Zz2"], "products": ["H2O"]},
            "known_reactions": [
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
            ],
        },
    )
    assert resp.status_code == 400
    error = resp.json()["error"]
    assert error["code"] == "FORMULA_PARSE_ERROR"
    assert error["details"]["where"] == "target"


def test_known_formula_parse_error_is_400_with_index():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["H2", "O2"], "products": ["H2O"]},
            "known_reactions": [
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
                {"reactants": ["Qx"], "products": ["H2O"],
                 "enthalpy_kj": "10"},
            ],
        },
    )
    assert resp.status_code == 400
    error = resp.json()["error"]
    assert error["code"] == "FORMULA_PARSE_ERROR"
    assert error["details"]["where"] == "known_reactions"
    assert error["details"]["index"] == 1


def test_known_duplicate_species_is_409():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["H2", "O2"], "products": ["H2O"]},
            "known_reactions": [
                {"reactants": ["H2", "H2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
            ],
        },
    )
    assert resp.status_code == 409
    error = resp.json()["error"]
    assert error["code"] == "DUPLICATE_SPECIES"
    assert error["details"]["where"] == "known_reactions"
    assert error["details"]["index"] == 0


def test_known_non_unique_balance_is_422():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["H2", "O2"], "products": ["H2O"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO", "CO2"],
                 "enthalpy_kj": "10"},
            ],
        },
    )
    assert resp.status_code == 422
    error = resp.json()["error"]
    assert error["code"] == "NON_UNIQUE_BALANCE"
    assert error["details"]["where"] == "known_reactions"
    assert error["details"]["index"] == 0


def test_target_impossible_balance_is_422():
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["H2"], "products": ["CO2"]},
            "known_reactions": [
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
            ],
        },
    )
    assert resp.status_code == 422
    error = resp.json()["error"]
    assert error["code"] == "IMPOSSIBLE_BALANCE"
    assert error["details"]["where"] == "target"


# ------------------------------------------------------------- 离子反应
def test_ionic_half_reactions():
    """两个还原半反应组合：Zn + Cu2+ -> Zn2+ + Cu，ΔH = E2 - E1（设值）。"""
    # Zn2+ + 2e- -> Zn，ΔH = -153
    # Cu2+ + 2e- -> Cu，ΔH = -64
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["Zn", "Cu2+"], "products": ["Zn2+", "Cu"]},
            "known_reactions": [
                {"reactants": ["Zn2+", "e-"], "products": ["Zn"],
                 "enthalpy_kj": "-153"},
                {"reactants": ["Cu2+", "e-"], "products": ["Cu"],
                 "enthalpy_kj": "-64"},
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["target"]["equation"] == "Cu^2+ + Zn -> Cu + Zn^2+"
    assert [k["multiplier"] for k in data["known_reactions"]] == ["-1", "1"]
    assert data["derived_enthalpy_kj"] == "89"
    net = {item["normalized"]: item["coefficient"] for item in data["net_coefficients"]}
    assert net["e^-"] == "0"
    assert net["Zn"] == "-1"
    assert net["Cu^2+"] == "-1"
    assert net["Zn^2+"] == "1"
    assert net["Cu"] == "1"


def test_response_matches_pydantic_model():
    """成功响应字段齐备且类型正确。"""
    resp = client.post(
        PATH,
        json={
            "target": {"reactants": ["C", "H2"], "products": ["CH4"]},
            "known_reactions": [
                {"reactants": ["C", "O2"], "products": ["CO2"],
                 "enthalpy_kj": "-393.5"},
                {"reactants": ["H2", "O2"], "products": ["H2O"],
                 "enthalpy_kj": "-571.6"},
                {"reactants": ["CH4", "O2"], "products": ["CO2", "H2O"],
                 "enthalpy_kj": "-890.3"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert set(data) == {
        "target",
        "known_reactions",
        "derived_enthalpy_kj",
        "net_coefficients",
    }
    known = data["known_reactions"][0]
    assert set(known) == {
        "index",
        "reactants",
        "products",
        "equation",
        "enthalpy_kj",
        "multiplier",
        "multiplier_used",
    }
    assert known["index"] == 0
    item = data["net_coefficients"][0]
    assert set(item) == {"normalized", "coefficient", "side"}
    assert item["side"] in {"product", "reactant", "cancelled"}


def test_hess_is_stateless():
    """重复请求结果完全一致，且不返回 reaction_id / cached 字段。"""
    body = {
        "target": {"reactants": ["C", "O2"], "products": ["CO2"]},
        "known_reactions": [
            {"reactants": ["C", "O2"], "products": ["CO2"],
             "enthalpy_kj": "-393.5"},
        ],
    }
    r1 = client.post(PATH, json=body)
    r2 = client.post(PATH, json=body)
    assert r1.json() == r2.json()
    assert "reaction_id" not in r1.json()
