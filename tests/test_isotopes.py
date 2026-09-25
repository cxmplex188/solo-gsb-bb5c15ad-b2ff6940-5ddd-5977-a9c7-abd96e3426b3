"""同位素名义质量峰分布的核心算法与 HTTP API 测试。

核心（多项分布枚举、同质量合并、跨元素卷积、零概率峰保留、并列取
小质量、规模上限内的性能与概率守恒）直接测试 ``app.isotopes``；
HTTP 层测试 ``POST /api/v1/isotopes/envelope`` 的成功路径、最简分数
输出以及 400/422 输入校验。
"""

import math
import time
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.isotopes import (
    ElementIsotopes,
    Isotope,
    element_mass_distribution,
    isotope_envelope,
    most_probable_peak,
)
from app.main import app

client = TestClient(app)

PATH = "/api/v1/isotopes/envelope"


# ------------------------------------------------------------ 核心算法辅助
def iso_spec(element, n, pairs):
    """``pairs`` 为 (质量数, 分数字符串) 序列。"""
    return ElementIsotopes(
        element,
        n,
        tuple(Isotope(mass, Fraction(text)) for mass, text in pairs),
    )


# ------------------------------------------------------------ 核心算法：单元素多项分布
def test_single_atom_uniform_distribution():
    peaks = element_mass_distribution(iso_spec("X", 1, [(5, "1/2"), (6, "1/2")]))
    assert [(p.mass, p.probability) for p in peaks] == [
        (5, Fraction(1, 2)),
        (6, Fraction(1, 2)),
    ]


def test_binomial_two_isotopes_merges_same_mass_vectors():
    # Cl2：35/37 各一个的两种计数向量 (1,1) 内部即同一质量 72，
    # 多项系数 2 使概率为 2·3/4·1/4 = 3/8
    peaks = element_mass_distribution(iso_spec("Cl", 2, [(35, "3/4"), (37, "1/4")]))
    assert [(p.mass, p.probability) for p in peaks] == [
        (70, Fraction(9, 16)),
        (72, Fraction(3, 8)),
        (74, Fraction(1, 16)),
    ]
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


def test_trinomial_distribution_sums_to_one():
    peaks = element_mass_distribution(
        iso_spec("O", 4, [(16, "1/2"), (17, "1/3"), (18, "1/6")])
    )
    # 质量种数 = n·(k-1)+1 = 4·2+1 = 9
    assert [p.mass for p in peaks] == list(range(64, 73))
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


def test_coefficients_are_multinomial():
    # 两种等概率同位素、3 个原子：1/8, 3/8, 3/8, 1/8
    peaks = element_mass_distribution(iso_spec("X", 3, [(1, "1/2"), (2, "1/2")]))
    assert [str(p.probability) for p in peaks] == ["1/8", "3/8", "3/8", "1/8"]


def test_probabilities_are_reduced_fractions():
    # 0.75 = 3/4、0.25 = 1/4：两原子 35Cl2 的概率 9/16 已最简
    peaks = element_mass_distribution(iso_spec("Cl", 2, [(35, "3/4"), (37, "1/4")]))
    for peak in peaks:
        assert peak.probability == Fraction(peak.probability.numerator,
                                            peak.probability.denominator)
        assert Fraction(peak.probability.numerator,
                        peak.probability.denominator).denominator == \
            peak.probability.denominator


def test_zero_abundance_isotope_keeps_zero_probability_masses():
    # ¹⁴C 丰度 0：质量 41（一个 14C）与 42（两个 14C）可达但概率恰为 0，
    # 不得在合并前截断
    peaks = element_mass_distribution(iso_spec("C", 2, [(12, "1"), (14, "0")]))
    assert [(p.mass, p.probability) for p in peaks] == [
        (24, Fraction(1)),
        (26, Fraction(0)),
        (28, Fraction(0)),
    ]


def test_all_zero_abundance_except_one_gives_single_peak():
    peaks = element_mass_distribution(
        iso_spec("C", 3, [(12, "1"), (13, "0"), (14, "0")])
    )
    nonzero = [(p.mass, p.probability) for p in peaks if p.probability != 0]
    assert nonzero == [(36, Fraction(1))]
    # 零概率质量同样保留
    assert any(p.probability == 0 for p in peaks)
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


def test_zero_atoms_element_is_defensive_unit_peak():
    peaks = element_mass_distribution(iso_spec("X", 0, [(1, "1")]))
    assert [(p.mass, p.probability) for p in peaks] == [(0, Fraction(1))]


# ------------------------------------------------------------ 核心算法：跨元素卷积
def test_convolution_of_two_elements():
    # H: 1H/2H 各 1/2；Cl: 35Cl 3/4、37Cl 1/4 -> 四峰 3/8,3/8,1/8,1/8
    peaks = isotope_envelope([
        iso_spec("H", 1, [(1, "1/2"), (2, "1/2")]),
        iso_spec("Cl", 1, [(35, "3/4"), (37, "1/4")]),
    ])
    assert [(p.mass, str(p.probability)) for p in peaks] == [
        (36, "3/8"), (37, "3/8"), (38, "1/8"), (39, "1/8"),
    ]
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


def test_convolution_merges_cross_element_same_mass():
    # 元素 A（1 原子）：质量 1 -> 1/2，质量 3 -> 1/2
    # 元素 B（1 原子）：质量 2 -> 1/2，质量 4 -> 1/2
    # 质量 5 来自 (1+4) 与 (3+2) 两条路径，合并为 1/2
    peaks = isotope_envelope([
        iso_spec("A", 1, [(1, "1/2"), (3, "1/2")]),
        iso_spec("B", 1, [(2, "1/2"), (4, "1/2")]),
    ])
    by_mass = {p.mass: p.probability for p in peaks}
    assert by_mass == {3: Fraction(1, 4), 5: Fraction(1, 2), 7: Fraction(1, 4)}
    assert [p.mass for p in peaks] == sorted(by_mass)


def test_three_elements_total_probability_one():
    peaks = isotope_envelope([
        iso_spec("C", 10, [(12, "0.9893"), (13, "0.0107")]),
        iso_spec("H", 10, [(1, "0.9998"), (2, "0.0002")]),
        iso_spec("O", 10, [(16, "0.9976"), (17, "0.0004"), (18, "0.002")]),
    ])
    assert [p.mass for p in peaks] == sorted(p.mass for p in peaks)
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


def test_peaks_always_sorted_by_mass_ascending():
    peaks = isotope_envelope([
        iso_spec("H", 2, [(1, "1/3"), (2, "2/3")]),
        iso_spec("O", 2, [(16, "1/4"), (17, "1/2"), (18, "1/4")]),
    ])
    masses = [p.mass for p in peaks]
    assert masses == sorted(masses)
    assert sum((p.probability for p in peaks), Fraction(0)) == 1


# ------------------------------------------------------------ 核心算法：最高概率峰
def test_most_probable_picks_highest_probability():
    peaks = isotope_envelope([
        iso_spec("Cl", 2, [(35, "3/4"), (37, "1/4")]),
    ])
    top = most_probable_peak(peaks)
    assert (top.mass, top.probability) == (70, Fraction(9, 16))


def test_most_probable_tie_prefers_smaller_mass():
    # HCl 例中质量 36 与 37 概率并列 3/8：取质量较小者 36
    peaks = isotope_envelope([
        iso_spec("H", 1, [(1, "1/2"), (2, "1/2")]),
        iso_spec("Cl", 1, [(35, "3/4"), (37, "1/4")]),
    ])
    assert most_probable_peak(peaks).mass == 36


def test_most_probable_empty_raises():
    with pytest.raises(ValueError):
        most_probable_peak([])


# ------------------------------------------------------------ 核心算法：规模上限与性能
def test_max_scale_enumeration_is_fast_and_exact():
    # 30 个原子、3 种同位素：单元素质量种数 61，状态数 465
    start = time.perf_counter()
    peaks = element_mass_distribution(
        iso_spec("C", 30, [(12, "0.9893"), (13, "0.0107"), (14, "0")])
    )
    elapsed = time.perf_counter() - start
    assert len(peaks) == 61
    assert elapsed < 2.0
    assert sum((p.probability for p in peaks), Fraction(0)) == 1
    # 最高概率峰即全部为 12C
    assert most_probable_peak(peaks).mass == 360


# ------------------------------------------------------------ HTTP：请求构造
def element_block(symbol, isotopes):
    return {
        "element": symbol,
        "isotopes": [
            {"mass_number": mass, "abundance": abundance}
            for mass, abundance in isotopes
        ],
    }


def body(formula, blocks):
    return {"formula": formula, "elements": blocks}


# ------------------------------------------------------------ HTTP：成功路径
def test_hcl_envelope_fraction_strings():
    resp = client.post(PATH, json=body("HCl", [
        element_block("H", [(1, "0.5"), (2, "0.5")]),
        element_block("Cl", [(35, "0.75"), (37, "0.25")]),
    ]))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["input"] == "HCl"
    assert data["normalized"] == "HCl"
    assert data["total_atoms"] == 2
    assert data["elements"] == [
        {"element": "H", "atom_count": 1,
         "isotopes": [{"mass_number": 1, "abundance": "0.5"},
                      {"mass_number": 2, "abundance": "0.5"}]},
        {"element": "Cl", "atom_count": 1,
         "isotopes": [{"mass_number": 35, "abundance": "0.75"},
                      {"mass_number": 37, "abundance": "0.25"}]},
    ]
    assert data["peaks"] == [
        {"mass": 36, "probability": "3/8"},
        {"mass": 37, "probability": "3/8"},
        {"mass": 38, "probability": "1/8"},
        {"mass": 39, "probability": "1/8"},
    ]
    # 并列取质量较小者
    assert data["most_probable"] == {"mass": 36, "probability": "3/8"}


def test_water_three_elements_uses_parser_atom_counts():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "0.5"), (18, "0.5")]),
    ]))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["total_atoms"] == 3
    counts = {e["element"]: e["atom_count"] for e in data["elements"]}
    assert counts == {"H": 2, "O": 1}
    assert data["peaks"] == [
        {"mass": 18, "probability": "1/2"},
        {"mass": 20, "probability": "1/2"},
    ]
    assert data["most_probable"] == {"mass": 18, "probability": "1/2"}


def test_parenthesized_formula_atom_counts():
    resp = client.post(PATH, json=body("Ca(OH)2", [
        element_block("Ca", [(40, "1")]),
        element_block("O", [(16, "1")]),
        element_block("H", [(1, "1")]),
    ]))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["normalized"] == "Ca(OH)2"
    assert data["total_atoms"] == 5
    assert data["peaks"] == [{"mass": 74, "probability": "1"}]
    assert data["most_probable"]["probability"] == "1"


def test_charged_species_accepted():
    # 离子（非电子）是合法的，电荷不影响名义质量峰
    resp = client.post(PATH, json=body("SO4^2-", [
        element_block("S", [(32, "1")]),
        element_block("O", [(16, "1")]),
    ]))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["normalized"] == "SO4^2-"
    assert data["total_atoms"] == 5
    assert data["peaks"] == [{"mass": 96, "probability": "1"}]


def test_probability_one_rendered_as_plain_one():
    resp = client.post(PATH, json=body("H2", [
        element_block("H", [(1, "1")]),
    ]))
    data = resp.json()
    assert data["peaks"] == [{"mass": 2, "probability": "1"}]


def test_decimal_abundances_converted_exactly_and_reduced():
    # Cl2：0.75/0.25 -> 9/16、3/8、1/16（最简分数，非 75/100）
    resp = client.post(PATH, json=body("Cl2", [
        element_block("Cl", [(35, "0.75"), (37, "0.25")]),
    ]))
    assert resp.status_code == 200, resp.text
    assert resp.json()["peaks"] == [
        {"mass": 70, "probability": "9/16"},
        {"mass": 72, "probability": "3/8"},
        {"mass": 74, "probability": "1/16"},
    ]


def test_scientific_notation_and_zero_abundance_accepted():
    resp = client.post(PATH, json=body("H2", [
        element_block("H", [(1, "1e0"), (2, "0")]),
    ]))
    assert resp.status_code == 200, resp.text
    peaks = resp.json()["peaks"]
    # 质量 2（两个 1H）概率 1；质量 3、4 概率 0 仍保留
    assert peaks == [
        {"mass": 2, "probability": "1"},
        {"mass": 3, "probability": "0"},
        {"mass": 4, "probability": "0"},
    ]


def test_three_isotopes_fraction_probabilities():
    # 1 个原子、三种同位素丰度 0.2/0.3/0.5 -> 1/5、3/10、1/2（最简分数）
    resp = client.post(PATH, json=body("O", [
        element_block("O", [(16, "0.2"), (17, "0.3"), (18, "0.5")]),
    ]))
    assert resp.status_code == 200, resp.text
    assert resp.json()["peaks"] == [
        {"mass": 16, "probability": "1/5"},
        {"mass": 17, "probability": "3/10"},
        {"mass": 18, "probability": "1/2"},
    ]


def test_thirty_atoms_boundary_accepted():
    resp = client.post(PATH, json=body("H30", [
        element_block("H", [(1, "0.5"), (2, "0.5")]),
    ]))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["total_atoms"] == 30
    assert len(data["peaks"]) == 31
    # 概率恰为二项分布：中间两峰（44、45）并列，取较小质量
    assert data["most_probable"]["mass"] == 45
    assert data["most_probable"]["probability"] == str(
        Fraction(math.comb(30, 15), 2 ** 30)
    )


def test_response_shape():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "0.5"), (2, "0.5")]),
        element_block("O", [(16, "1")]),
    ]))
    data = resp.json()
    assert set(data) == {
        "input", "normalized", "total_atoms", "elements", "peaks", "most_probable",
    }
    assert set(data["peaks"][0]) == {"mass", "probability"}
    assert set(data["most_probable"]) == {"mass", "probability"}


def test_stateless_same_request_same_result():
    payload = body("HCl", [
        element_block("H", [(1, "0.5"), (2, "0.5")]),
        element_block("Cl", [(35, "0.75"), (37, "0.25")]),
    ])
    first = client.post(PATH, json=payload)
    second = client.post(PATH, json=payload)
    assert first.status_code == second.status_code == 200
    assert first.json() == second.json()


def test_input_order_of_elements_does_not_change_peaks():
    def post(blocks):
        return client.post(PATH, json=body("HCl", blocks)).json()["peaks"]
    peaks_ab = post([
        element_block("H", [(1, "0.5"), (2, "0.5")]),
        element_block("Cl", [(35, "0.75"), (37, "0.25")]),
    ])
    peaks_ba = post([
        element_block("Cl", [(35, "0.75"), (37, "0.25")]),
        element_block("H", [(1, "0.5"), (2, "0.5")]),
    ])
    assert peaks_ab == peaks_ba


# ------------------------------------------------------------ HTTP：化学式问题（400）
def test_unparseable_formula_returns_400():
    resp = client.post(PATH, json=body("H2!", [
        element_block("H", [(1, "1")]),
    ]))
    assert resp.status_code == 400
    err = resp.json()["error"]
    assert err["code"] == "FORMULA_PARSE_ERROR"
    assert err["details"]["field"] == "formula"


def test_electron_species_rejected_with_400():
    for text in ["e-", "e", "e+"]:
        resp = client.post(PATH, json=body(text, [
            element_block("H", [(1, "1")]),  # 即使给了无关元素也先拒绝电子
        ]))
        assert resp.status_code == 400, text
        assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


def test_unknown_element_symbol_in_block_is_422():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "1")]),
        element_block("Xx", [(10, "1")]),
    ]))
    # 3 个组但元素本身未知（不是结构问题，是语义问题）
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["invalid_elements"] == [{"element": "Xx"}]


# ------------------------------------------------------------ HTTP：元素组语义校验（422 INVALID_ISOTOPE）
def test_missing_element_group_rejected():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "1")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["missing_elements"] == ["O"]


def test_unexpected_element_group_rejected():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "1")]),
        element_block("N", [(14, "1")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["unexpected_elements"] == [
        {"element": "N", "formula": "H2O"}
    ]


def test_duplicate_element_group_rejected():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "0.5"), (2, "0.5")]),
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "1")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["duplicate_elements"] == ["H"]


def test_duplicate_mass_numbers_rejected():
    resp = client.post(PATH, json=body("H2", [
        element_block("H", [(1, "0.5"), (1, "0.5")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["duplicate_mass_numbers"] == [
        {"element": "H", "mass_numbers": [1]}
    ]


@pytest.mark.parametrize("bad", ["-1", "-0.0", "NaN", "Infinity", "abc", "1e20000"])
def test_invalid_abundance_rejected(bad):
    resp = client.post(PATH, json=body("H", [
        element_block("H", [(1, bad)]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["invalid_abundances"] == [
        {"element": "H", "mass_number": 1, "abundance": bad}
    ]


def test_whitespace_padded_abundance_is_trimmed_then_validated():
    # 两侧空白由模型剥除：合法值照常接受，回显为剥除后的文本
    resp = client.post(PATH, json=body("H", [
        element_block("H", [(1, "  0.5"), (2, " 0.5 ")]),
    ]))
    assert resp.status_code == 200, resp.text
    assert [iso["abundance"] for iso in resp.json()["elements"][0]["isotopes"]] == [
        "0.5", "0.5"
    ]


def test_abundance_sum_not_one_rejected_with_exact_fraction():
    resp = client.post(PATH, json=body("H", [
        element_block("H", [(1, "0.4"), (2, "0.5")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    # 0.4 + 0.5 = 9/10（精确分数，非 0.9 浮点）
    assert err["details"]["abundance_sums"] == [
        {"element": "H", "sum": "9/10"}
    ]


def test_abundance_sum_must_be_exact_one():
    # 三个 0.3333333333 之和为 0.9999999999，不等于精确的 1
    resp = client.post(PATH, json=body("O", [
        element_block("O", [(16, "0.3333333333"),
                            (17, "0.3333333333"),
                            (18, "0.3333333333")]),
    ]))
    assert resp.status_code == 422
    sums = resp.json()["error"]["details"]["abundance_sums"]
    assert sums[0]["element"] == "O"
    assert sums[0]["sum"] == "9999999999/10000000000"


def test_abundance_sum_one_accepted_with_long_decimals():
    resp = client.post(PATH, json=body("O", [
        element_block("O", [(16, "0.99757"), (17, "0.00038"), (18, "0.00205")]),
    ]))
    assert resp.status_code == 200, resp.text


def test_too_many_total_atoms_rejected():
    resp = client.post(PATH, json=body("H31", [
        element_block("H", [(1, "1")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["total_atoms"] == 31
    assert err["details"]["max_total_atoms"] == 30


def test_formula_with_four_elements_rejected():
    # 请求模型允许 1～3 个元素组；4 元素化学式无法满足一一对应
    resp = client.post(PATH, json=body("CHON", [
        element_block("C", [(12, "1")]),
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "1")]),
    ]))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_ISOTOPE"
    assert err["details"]["missing_elements"] == ["N"]


def test_multiple_semantic_problems_aggregated():
    resp = client.post(PATH, json=body("H2O", [
        element_block("H", [(1, "0.4"), (1, "0.5")]),  # 重复质量数 + 和不为 1
        element_block("N", [(14, "1")]),               # 不存在的元素
    ]))
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["missing_elements"] == ["O"]
    assert details["unexpected_elements"] == [{"element": "N", "formula": "H2O"}]
    assert details["duplicate_mass_numbers"] == [
        {"element": "H", "mass_numbers": [1]}
    ]
    # 和不为 1 仍然报告（丰度本身都能解析）
    assert details["abundance_sums"] == [{"element": "H", "sum": "9/10"}]


# ------------------------------------------------------------ HTTP：请求体结构校验（422 INVALID_REQUEST）
@pytest.mark.parametrize("payload", [
    {},
    {"formula": "H2O"},
    {"elements": []},
    {"formula": "H2O", "elements": []},
    {"formula": "", "elements": []},
    {"formula": "   ", "elements": [element_block("H", [(1, "1")])]},
    {"formula": "H2", "elements": [{"element": "H", "isotopes": []}]},
    {"formula": "H2", "elements": [{"element": " ", "isotopes": [{"mass_number": 1, "abundance": "1"}]}]},
    {"formula": "H2", "elements": [{"element": "H", "isotopes": [{"mass_number": 0, "abundance": "1"}]}]},
    {"formula": "H2", "elements": [{"element": "H", "isotopes": [{"mass_number": -1, "abundance": "1"}]}]},
    {"formula": "H2", "elements": [{"element": "H", "isotopes": [{"mass_number": 1.5, "abundance": "1"}]}]},
    {"formula": "H2", "elements": [{"element": "H", "isotopes": [{"mass_number": 1, "abundance": 1}]}]},
    {"formula": "H2", "elements": [
        element_block("H", [(1, "1")]),
        element_block("H", [(1, "1")]),
        element_block("H", [(1, "1")]),
        element_block("H", [(1, "1")]),
    ]},
    {"formula": "H2", "elements": [
        {"element": "H", "isotopes": [
            {"mass_number": 1, "abundance": "1/3"},
            {"mass_number": 2, "abundance": "1/3"},
            {"mass_number": 3, "abundance": "1/3"},
            {"mass_number": 4, "abundance": "0"},
        ]},
    ]},
])
def test_structural_validation_errors(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422, resp.text
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_whitespace_abundance_string_is_structural_error():
    resp = client.post(PATH, json=body("H", [
        element_block("H", [(1, "   ")]),
    ]))
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_formula_whitespace_is_trimmed_before_parse():
    resp = client.post(PATH, json=body("  H2O  ", [
        element_block("H", [(1, "1")]),
        element_block("O", [(16, "1")]),
    ]))
    assert resp.status_code == 200
    # 请求模型剥除两侧空白后交给解析器
    assert resp.json()["input"] == "H2O"
    assert resp.json()["normalized"] == "H2O"
