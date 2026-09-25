"""POST /api/v1/stoichiometry/yield 的 HTTP API 测试。

覆盖：逐种理论/实际产量与收率百分比、规范化物种匹配、>100% 收率保留、
理论产量为 0 时收率 null 与原因说明、精确十进制（无限小数 28 位有效
数字），以及请求校验（缺失/重复/未知物种、摩尔数非法、JSON 数字、
配平错误透传）。
"""

from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.stoichiometry import fraction_to_decimal_string

client = TestClient(app)

PATH = "/api/v1/stoichiometry/yield"


# ------------------------------------------------------------ 成功路径
def test_basic_theoretical_actual_and_percent_yield():
    # N2 + 3H2 -> 2NH3；N2=1、H2=5 -> ξ=1（N2 限量），理论 NH3=2
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "5"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "1.5"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"] == {"value": "1", "unit": "mol"}
    assert data["limiting_reactants"] == ["N2"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["N2"]["fed"] == "1"
    assert by_name["N2"]["consumed"] == "1"
    assert by_name["N2"]["remaining"] == "0"
    assert by_name["H2"]["fed"] == "5"
    assert by_name["H2"]["consumed"] == "3"
    assert by_name["H2"]["remaining"] == "2"

    assert len(data["products"]) == 1
    product = data["products"][0]
    assert product["normalized"] == "NH3"
    assert product["coefficient"] == 2
    assert product["theoretical_yield"] == "2"
    assert product["actual_yield"] == "1.5"
    assert product["yield_percent"] == "75"
    assert product["yield_note"] is None
    # 数值字段均以十进制字符串返回
    assert isinstance(product["theoretical_yield"], str)
    assert isinstance(product["actual_yield"], str)
    assert isinstance(product["yield_percent"], str)


def test_yield_uses_normalized_species_matching():
    # 2HCl + Ca(OH)2 -> CaCl2 + 2H2O；等价写法（方括号/内部空白）应命中同一物种
    resp = client.post(
        PATH,
        json={
            "reactants": ["H Cl", "Ca [OH]2"],
            "products": ["Ca Cl2", "H2 O"],
            "amounts": [
                {"species": "HCl", "moles": "2"},
                {"species": "Ca(OH)2", "moles": "1"},
            ],
            "actual_yields": [
                {"species": "CaCl2", "moles": "0.8"},
                {"species": "H2O", "moles": "2"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["normalized"]["equation"] == "Ca(OH)2 + 2HCl -> CaCl2 + 2H2O"
    assert data["reaction_extent"]["value"] == "1"
    products = {p["normalized"]: p for p in data["products"]}
    assert products["CaCl2"]["theoretical_yield"] == "1"
    assert products["CaCl2"]["actual_yield"] == "0.8"
    assert products["CaCl2"]["yield_percent"] == "80"
    assert products["CaCl2"]["input"] == "CaCl2"
    assert products["H2O"]["theoretical_yield"] == "2"
    assert products["H2O"]["yield_percent"] == "100"


def test_actual_above_theoretical_keeps_over_100_percent():
    # 实际 2.5 mol / 理论 2 mol = 125%，不得封顶为 100
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "2.5"}],
        },
    )
    assert resp.status_code == 200
    product = resp.json()["products"][0]
    assert product["theoretical_yield"] == "2"
    assert product["actual_yield"] == "2.5"
    assert product["yield_percent"] == "125"
    assert product["yield_note"] is None


def test_zero_actual_yield_is_zero_percent():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "0"}],
        },
    )
    assert resp.status_code == 200
    product = resp.json()["products"][0]
    assert product["theoretical_yield"] == "2"
    assert product["actual_yield"] == "0"
    assert product["yield_percent"] == "0"
    assert product["yield_note"] is None


def test_zero_theoretical_yield_returns_null_percent_with_reason():
    # N2 投料为 0 -> ξ = 0，所有生成物理论产量为 0，收率无定义
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "0"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "0"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0"
    assert data["limiting_reactants"] == ["N2"]
    product = data["products"][0]
    assert product["theoretical_yield"] == "0"
    assert product["actual_yield"] == "0"
    assert product["yield_percent"] is None
    assert product["yield_note"] is not None
    assert "理论产量为 0" in product["yield_note"]
    assert "N2" in product["yield_note"]


def test_all_feeds_zero_lists_every_limiting_reactant_in_note():
    # 全部投料为 0 时全部反应物并列限量，原因说明须逐一列出
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "0"},
                {"species": "H2", "moles": "0"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "0"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["limiting_reactants"] == ["H2", "N2"]
    note = data["products"][0]["yield_note"]
    assert "N2" in note and "H2" in note
    assert data["products"][0]["yield_percent"] is None


def test_repeating_decimal_percent_uses_exact_arithmetic():
    # 2H2 + O2 -> 2H2O；H2=3、O2=3 -> ξ=1.5，理论 H2O=3，实际 2
    # 收率 = 2/3 × 100 = 200/3（无限循环小数，28 位有效数字银行家舍入）
    resp = client.post(
        PATH,
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "amounts": [
                {"species": "H2", "moles": "3"},
                {"species": "O2", "moles": "3"},
            ],
            "actual_yields": [{"species": "H2O", "moles": "2"}],
        },
    )
    assert resp.status_code == 200
    product = resp.json()["products"][0]
    assert product["theoretical_yield"] == "3"
    expected = fraction_to_decimal_string(Fraction(200, 3))
    assert product["yield_percent"] == expected
    # 28 位有效数字（整数 2 位 + 小数 26 位），末位 7 为银行家舍入
    assert product["yield_percent"] == "66." + "6" * 25 + "7"


def test_tied_limiting_reactants_still_computes_yield():
    # 投料恰按系数比时两种反应物并列限量，收率照常计算
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "2"},
                {"species": "H2", "moles": "6"},
            ],
            "actual_yields": [{"species": "NH3", "moles": "1"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["limiting_reactants"] == ["H2", "N2"]
    assert data["products"][0]["theoretical_yield"] == "4"
    assert data["products"][0]["yield_percent"] == "25"


# ------------------------------------------------------------ 结构化错误
def test_missing_actual_yield_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"


def test_missing_actual_yield_entry_for_product():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "O2", "moles": "3"},
            ],
            "actual_yields": [{"species": "CO2", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["actual_yields.missing"] == ["H2O"]


def test_duplicate_and_unknown_actual_yield_entries():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [
                {"species": "CH4", "moles": "1"},
                {"species": "O2", "moles": "3"},
            ],
            "actual_yields": [
                {"species": "CO2", "moles": "1"},
                {"species": "CO2", "moles": "0.5"},
                {"species": "N2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["actual_yields.duplicate"] == ["CO2"]
    assert details["actual_yields.unknown"] == ["N2"]
    assert details["actual_yields.missing"] == ["H2O"]


def test_feed_side_missing_and_actual_side_problem_aggregated_together():
    resp = client.post(
        PATH,
        json={
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "amounts": [{"species": "CH4", "moles": "1"}],  # 缺 O2
            "actual_yields": [
                {"species": "CO2", "moles": "1"},
                {"species": "H2O", "moles": "-2"},  # 摩尔数非法
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["missing"] == ["O2"]
    assert details["actual_yields.invalid_moles"] == [
        {"species": "H2O", "moles": "-2"}
    ]


@pytest.mark.parametrize("bad", ["-1", "abc", "1.2.3", "NaN", "Infinity", "1e10001"])
def test_invalid_actual_moles_values(bad):
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "NH3", "moles": bad}],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"
    assert body["error"]["details"]["actual_yields.invalid_moles"] == [
        {"species": "NH3", "moles": bad}
    ]


def test_moles_must_be_string_not_json_number():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "NH3", "moles": 0.5}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unparseable_actual_species_is_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
            "actual_yields": [{"species": "Xx", "moles": "1"}],
        },
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"] == {
        "field": "actual_yields.species",
        "input": "Xx",
    }


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
            "actual_yields": [
                {"species": "CO", "moles": "1"},
                {"species": "CO2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NON_UNIQUE_BALANCE"


def test_reaction_side_duplicate_reuses_conflict():
    resp = client.post(
        PATH,
        json={
            "reactants": ["H2", "H2"],
            "products": ["H2O"],
            "amounts": [
                {"species": "H2", "moles": "1"},
                {"species": "H2", "moles": "2"},
            ],
            "actual_yields": [{"species": "H2O", "moles": "1"}],
        },
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


@pytest.mark.parametrize(
    "payload",
    [
        # 缺 actual_yields
        {
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "3"},
            ],
        },
        # 缺 amounts
        {
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "actual_yields": [{"species": "NH3", "moles": "1"}],
        },
        # 反应式为空
        {
            "reactants": [],
            "products": ["NH3"],
            "amounts": [{"species": "N2", "moles": "1"}],
            "actual_yields": [{"species": "NH3", "moles": "1"}],
        },
        # actual_yields 为空列表
        {
            "reactants": ["N2"],
            "products": ["NH3"],
            "amounts": [{"species": "N2", "moles": "1"}],
            "actual_yields": [],
        },
    ],
)
def test_request_shape_validation(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"
