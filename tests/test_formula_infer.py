"""元素分析质量百分比反推分子式的核心算法与 HTTP API 测试。

核心（容差闭区间、整数计数枚举、总偏差排序与在线截取、与朴素四重循环
穷举逐案对比）直接测试 ``app.formula_infer``；HTTP 层测试
``POST /api/v1/formulas/infer`` 的成功路径、无解空列表、截断标志与
422 输入校验。
"""

import time
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.formula_infer import (
    ATOMIC_MASS,
    ELEMENTS,
    MAX_COUNT,
    MAX_MOLAR_MASS,
    RESULT_LIMIT,
    ElementMeasurement,
    infer_formulas,
    render_formula,
)
from app.main import app

client = TestClient(app)

PATH = "/api/v1/formulas/infer"


# ------------------------------------------------------------ 辅助
def measurements(pairs):
    """``pairs`` 为按 C/H/N/O 顺序的 (percent, tolerance) 十进制字符串对。"""
    return {
        element: ElementMeasurement(element, Fraction(percent), Fraction(tolerance))
        for element, (percent, tolerance) in zip(ELEMENTS, pairs)
    }


def brute_force(meas, mass_min, mass_max):
    """朴素四重循环穷举，作为核心枚举的参照实现。

    容差判定全部用预取的整数系数交叉相乘，仅对命中的候选构造总偏差
    Fraction，避免参照实现本身在全质量区间用例上过慢。
    """
    band_coef = {}
    for element in ELEMENTS:
        lower = meas[element].percent - meas[element].tolerance
        upper = meas[element].percent + meas[element].tolerance
        band_coef[element] = (
            lower.numerator, lower.denominator,
            upper.numerator, upper.denominator,
        )

    results = []
    for c in range(MAX_COUNT + 1):
        for h in range(MAX_COUNT + 1):
            for n in range(MAX_COUNT + 1):
                for o in range(MAX_COUNT + 1):
                    counts = (c, h, n, o)
                    if sum(counts) == 0:
                        continue
                    molar_mass = sum(
                        ATOMIC_MASS[e] * count
                        for e, count in zip(ELEMENTS, counts)
                    )
                    if not mass_min <= molar_mass <= mass_max:
                        continue
                    deviation = Fraction(0)
                    ok = True
                    for element, count in zip(ELEMENTS, counts):
                        q100 = 100 * ATOMIC_MASS[element] * count
                        ln, ld, un, ud = band_coef[element]
                        if q100 * ld < ln * molar_mass:
                            ok = False
                            break
                        if q100 * ud > un * molar_mass:
                            ok = False
                            break
                        deviation += abs(
                            Fraction(q100, molar_mass)
                            - meas[element].percent
                        )
                    if ok:
                        results.append((deviation, molar_mass, c, h, n, o))
    results.sort()
    return results


# ------------------------------------------------------------ render_formula
@pytest.mark.parametrize(
    "counts, expected",
    [
        ((1, 4, 0, 0), "CH4"),
        ((6, 12, 0, 6), "C6H12O6"),
        ((1, 4, 2, 1), "CH4N2O"),
        ((0, 0, 2, 0), "N2"),
        ((0, 1, 0, 0), "H"),
        ((30, 30, 30, 30), "C30H30N30O30"),
    ],
)
def test_render_formula(counts, expected):
    assert render_formula(*counts) == expected


# ------------------------------------------------------------ 核心：已知分子
def test_glucose_exact_mass_percent_single_match():
    meas = measurements([
        ("40", "0"), ("20/3", "0"), ("0", "0"), ("160/3", "0"),
    ])
    candidates, total = infer_formulas(meas, 180, 180)
    assert total == 1
    top = candidates[0]
    assert (top.c, top.h, top.n, top.o) == (6, 12, 0, 6)
    assert top.molar_mass == 180
    assert top.total_abs_deviation == 0
    assert top.mass_percent == {
        "C": Fraction(40),
        "H": Fraction(20, 3),
        "N": Fraction(0),
        "O": Fraction(160, 3),
    }


def test_urea_decimal_inputs_tolerance_only_match_in_range():
    # 尿素 CH4N2O，M=60；测量值为 2 位小数，存在非零总偏差
    meas = measurements([
        ("20", "0.5"), ("6.67", "0.01"),
        ("46.67", "0.01"), ("26.67", "0.01"),
    ])
    candidates, total = infer_formulas(meas, 58, 62)
    assert total == 1
    top = candidates[0]
    assert render_formula(top.c, top.h, top.n, top.o) == "CH4N2O"
    assert top.molar_mass == 60
    assert top.total_abs_deviation == Fraction(1, 100)


def test_methane_zero_tolerance_closed_boundary():
    # CH4，M=16：C=75%、H=25%，零容差闭区间恰命中
    meas = measurements([("75", "0"), ("25", "0"), ("0", "0"), ("0", "0")])
    candidates, total = infer_formulas(meas, 16, 16)
    assert total == 1
    assert (candidates[0].c, candidates[0].h) == (1, 4)


def test_methane_tolerance_boundary_inclusive():
    # 容差区间端点必须计入：测量值偏离恰好等于容差
    meas = measurements([("75.5", "0.5"), ("24.5", "0.5"),
                         ("0", "0"), ("0", "0")])
    candidates, total = infer_formulas(meas, 16, 16)
    assert total == 1
    assert candidates[0].total_abs_deviation == 1  # 0.5+0.5


def test_zero_count_formula_allowed():
    # 仅 H：容差允许 N/O 略大于 0 也能匹配 H 单质
    meas = measurements([("0", "1"), ("100", "0"), ("0", "1"), ("0", "1")])
    candidates, total = infer_formulas(meas, 1, 1)
    assert total == 1
    assert (candidates[0].c, candidates[0].h, candidates[0].n,
            candidates[0].o) == (0, 1, 0, 0)


def test_zero_zero_zero_zero_excluded():
    # 全 0 计数（M=0）不在正整数质量区间内，也不可能被枚举
    meas = measurements([("0", "100"), ("0", "100"),
                         ("0", "100"), ("0", "100")])
    _candidates, total = infer_formulas(meas, 1, MAX_MOLAR_MASS)
    assert total == 31 ** 4 - 1


# ------------------------------------------------------------ 核心：排序
def test_sorting_by_deviation_then_mass_then_counts():
    # 测量值全 0、容差全 100：所有非零计数都可行；总偏差恒为 100
    # （四个计算百分比之和恒为 100、测量值均为 0），故第一排序键失效，
    # 排序退化为先按摩尔质量、再按 C/H/N/O 计数
    meas = measurements([("0", "100"), ("0", "100"),
                         ("0", "100"), ("0", "100")])
    candidates, total = infer_formulas(meas, 1, 20)
    assert total == 41
    keys = [
        (candidate.total_abs_deviation, candidate.molar_mass,
         candidate.c, candidate.h, candidate.n, candidate.o)
        for candidate in candidates
    ]
    assert keys == sorted(keys)
    # 摩尔质量 1 只有 H（1,0,0,0 质量为 12），故最小键为 H、M=1
    assert (candidates[0].c, candidates[0].h, candidates[0].n,
            candidates[0].o) == (0, 1, 0, 0)
    assert candidates[0].molar_mass == 1
    assert all(candidate.total_abs_deviation == 100 for candidate in candidates)


def test_mass_tiebreak_before_counts():
    # M=16 下 O 与 H16 总偏差并列 100、按摩尔质量无法区分；计数次序
    # C,H,N,O：(0,0,0,1) O 中 H=0 < (0,16,0,0) H16 中 H=16，故 O 在前。
    # CH4（偏差 150）、H2N（175）随后。
    meas = measurements([
        ("0", "75"), ("50", "50"), ("0", "100"), ("50", "50"),
    ])
    candidates, total = infer_formulas(meas, 16, 16)
    formulas = [
        render_formula(c.c, c.h, c.n, c.o) for c in candidates
    ]
    assert formulas == ["O", "H16", "CH4", "H2N"]
    assert total == 4
    assert candidates[0].total_abs_deviation == \
        candidates[1].total_abs_deviation == 100


def test_truncation_returns_twenty_and_total():
    meas = measurements([("50", "100"), ("50", "100"),
                         ("50", "100"), ("50", "100")])
    candidates, total = infer_formulas(meas, 1, MAX_MOLAR_MASS)
    assert total == 31 ** 4 - 1
    assert len(candidates) == RESULT_LIMIT == 20
    keys = [
        (c.total_abs_deviation, c.molar_mass, c.c, c.h, c.n, c.o)
        for c in candidates
    ]
    assert keys == sorted(keys)


# ------------------------------------------------------------ 核心：无解
def test_no_match_returns_empty():
    # C=100% 但质量上限 11：唯一 100% C 的质量是 12 的倍数，区间内不存在
    meas = measurements([("100", "0"), ("0", "0"), ("0", "0"), ("0", "0")])
    candidates, total = infer_formulas(meas, 1, 11)
    assert candidates == []
    assert total == 0


def test_mass_range_above_max_molar_mass_is_empty():
    meas = measurements([("50", "100"), ("50", "100"),
                         ("50", "100"), ("50", "100")])
    candidates, total = infer_formulas(
        meas, MAX_MOLAR_MASS + 1, MAX_MOLAR_MASS + 10
    )
    assert candidates == [] and total == 0


# ------------------------------------------------------------ 核心：穷举交叉验证
@pytest.mark.parametrize("seed", range(6))
def test_matches_naive_brute_force(seed):
    # 简单确定性伪随机（不依赖固定的全局随机种子）
    state = 2463534242 ^ (seed * 2654435761 & 0xFFFFFFFF)

    def rnd():
        nonlocal state
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        state &= 0xFFFFFFFF
        return state

    pairs = [
        (Fraction(rnd() % 1001, 10), Fraction(rnd() % 301, 10))
        for _ in ELEMENTS
    ]
    meas = {
        element: ElementMeasurement(element, percent, tolerance)
        for element, (percent, tolerance) in zip(ELEMENTS, pairs)
    }
    mass_min = 1 + rnd() % 400
    mass_max = min(MAX_MOLAR_MASS, mass_min + rnd() % 400)

    candidates, total = infer_formulas(meas, mass_min, mass_max)
    expected = brute_force(meas, mass_min, mass_max)
    assert total == len(expected)
    assert [
        (c.total_abs_deviation, c.molar_mass, c.c, c.h, c.n, c.o)
        for c in candidates
    ] == expected[:RESULT_LIMIT]


def test_matches_brute_force_full_mass_range():
    # 全质量区间的穷举交叉验证（另由截断测试覆盖最宽容差场景）
    meas = measurements([
        ("12.3", "4.0"), ("7.7", "3.0"),
        ("30.1", "8.0"), ("49.9", "6.0"),
    ])
    candidates, total = infer_formulas(meas, 1, MAX_MOLAR_MASS)
    expected = brute_force(meas, 1, MAX_MOLAR_MASS)
    assert total == len(expected)
    assert [
        (c.total_abs_deviation, c.molar_mass, c.c, c.h, c.n, c.o)
        for c in candidates
    ] == expected[:RESULT_LIMIT]


def test_worst_case_performance():
    # 最宽容差 + 全质量区间：923,520 个非零向量，核心枚举须在合理时限内
    meas = measurements([("50", "100"), ("50", "100"),
                         ("50", "100"), ("50", "100")])
    start = time.perf_counter()
    _candidates, total = infer_formulas(meas, 1, MAX_MOLAR_MASS)
    elapsed = time.perf_counter() - start
    assert total == 31 ** 4 - 1
    assert elapsed < 8.0


# ------------------------------------------------------------ HTTP：请求构造
def element_block(element, percent, tolerance):
    return {"element": element, "percent": percent, "tolerance": tolerance}


def body(pairs, mass_min, mass_max):
    return {
        "elements": [
            element_block(element, percent, tolerance)
            for element, (percent, tolerance) in zip(ELEMENTS, pairs)
        ],
        "molar_mass": {"minimum": mass_min, "maximum": mass_max},
    }


# ------------------------------------------------------------ HTTP：成功路径
def test_http_glucose_match():
    resp = client.post(PATH, json=body([
        ("40.00", "0.05"), ("6.71", "0.05"),
        ("0", "0.05"), ("53.29", "0.05"),
    ], 170, 190))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["count"] == data["total"] == 1
    assert data["truncated"] is False
    match = data["matches"][0]
    assert match["formula"] == "C6H12O6"
    assert match["molar_mass"] == 180
    assert match["counts"] == {"C": 6, "H": 12, "N": 0, "O": 6}
    assert match["mass_percent"]["C"] == "40"
    assert match["mass_percent"]["N"] == "0"
    assert set(match["mass_percent"]) == {"C", "H", "N", "O"}
    assert match["total_abs_deviation"] == "0.08666666666666666666666666667"


def test_http_response_shape_and_echo():
    payload = body([
        ("20", "0.5"), ("6.67", "0.01"),
        ("46.67", "0.01"), ("26.67", "0.01"),
    ], 58, 62)
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert set(data) == {
        "elements", "molar_mass", "count", "total", "truncated", "matches",
    }
    assert data["elements"] == payload["elements"]
    assert data["molar_mass"] == {"minimum": 58, "maximum": 62}
    assert data["matches"][0]["formula"] == "CH4N2O"
    assert set(data["matches"][0]) == {
        "formula", "molar_mass", "counts", "mass_percent",
        "total_abs_deviation",
    }


def test_http_measurements_need_not_sum_to_one_hundred():
    # 合计恰为 99.95（典型实验误差），葡萄糖各计算百分比仍落入容差：
    # C 40∈[39.95,40.05]、H 6.6667∈[6.62,6.72]、O 53.3333∈[53.24,53.34]
    resp = client.post(PATH, json=body([
        ("39.99", "0.06"), ("6.67", "0.05"),
        ("0", "0.02"), ("53.29", "0.05"),
    ], 179, 181))
    assert resp.status_code == 200, resp.text
    assert resp.json()["matches"][0]["formula"] == "C6H12O6"


def test_http_scientific_notation_decimals():
    resp = client.post(PATH, json=body([
        ("7.5e1", "0"), ("2.5e1", "0"), ("0", "0"), ("0", "0"),
    ], 16, 16))
    assert resp.status_code == 200, resp.text
    assert resp.json()["matches"][0]["formula"] == "CH4"


def test_http_no_match_empty_list():
    resp = client.post(PATH, json=body([
        ("100", "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 1, 11))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["matches"] == []
    assert data["count"] == 0 and data["total"] == 0
    assert data["truncated"] is False


def test_http_truncated_flag():
    resp = client.post(PATH, json=body([
        ("50", "100"), ("50", "100"), ("50", "100"), ("50", "100"),
    ], 1, MAX_MOLAR_MASS))
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["total"] == 31 ** 4 - 1
    assert data["count"] == 20
    assert data["truncated"] is True


def test_http_wide_range_completes_quickly():
    start = time.perf_counter()
    resp = client.post(PATH, json=body([
        ("50", "100"), ("50", "100"), ("50", "100"), ("50", "100"),
    ], 1, MAX_MOLAR_MASS))
    assert time.perf_counter() - start < 10
    assert resp.status_code == 200


def test_http_stateless_same_request_same_result():
    payload = body([
        ("20", "0.5"), ("6.67", "0.01"),
        ("46.67", "0.01"), ("26.67", "0.01"),
    ], 58, 62)
    first = client.post(PATH, json=payload)
    second = client.post(PATH, json=payload)
    assert first.status_code == second.status_code == 200
    assert first.json() == second.json()


def test_http_sorting_order_in_response():
    resp = client.post(PATH, json=body([
        ("50", "100"), ("50", "100"), ("50", "100"), ("50", "100"),
    ], 1, 30))
    data = resp.json()
    triples = [
        (Fraction(m["total_abs_deviation"]), m["molar_mass"],
         (m["counts"]["C"], m["counts"]["H"], m["counts"]["N"],
          m["counts"]["O"]))
        for m in data["matches"]
    ]
    assert triples == sorted(triples)


# ------------------------------------------------------------ HTTP：语义校验（422 INVALID_FORMULA_INFER）
def test_http_percent_above_one_hundred_rejected():
    resp = client.post(PATH, json=body([
        ("100.1", "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 1, 100))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_FORMULA_INFER"
    assert err["details"]["invalid_percent"] == [
        {"element": "C", "percent": "100.1"}
    ]


@pytest.mark.parametrize("bad", ["-1", "-0", "-0.0", "NaN", "Infinity", "abc"])
def test_http_invalid_percent_rejected(bad):
    resp = client.post(PATH, json=body([
        (bad, "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 1, 100))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_FORMULA_INFER"
    assert err["details"]["invalid_percent"] == [
        {"element": "C", "percent": bad}
    ]


@pytest.mark.parametrize("bad", ["-0.5", "100.01", "NaN"])
def test_http_invalid_tolerance_rejected(bad):
    resp = client.post(PATH, json=body([
        ("0", bad), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 1, 100))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_FORMULA_INFER"
    assert err["details"]["invalid_tolerance"] == [
        {"element": "C", "tolerance": bad}
    ]


def test_http_missing_and_duplicate_elements_reported():
    payload = body([
        ("0", "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 1, 100)
    payload["elements"][2]["element"] = "C"  # C 重复、N 缺失
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["missing"] == ["N"]
    assert details["duplicate"] == ["C"]


def test_http_inverted_mass_range_rejected():
    resp = client.post(PATH, json=body([
        ("0", "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 62, 58))
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_FORMULA_INFER"
    assert err["details"]["invalid_molar_mass_range"] == {
        "minimum": 62, "maximum": 58
    }


def test_http_multiple_problems_aggregated():
    payload = body([
        ("0", "0"), ("0", "0"), ("0", "0"), ("0", "0"),
    ], 62, 58)
    payload["elements"][0]["percent"] = "120"
    payload["elements"][1]["tolerance"] = "-1"
    payload["elements"][3]["element"] = "H"  # H 重复、O 缺失
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["missing"] == ["O"]
    assert details["duplicate"] == ["H"]
    assert details["invalid_percent"] == [{"element": "C", "percent": "120"}]
    assert details["invalid_tolerance"] == [
        {"element": "H", "tolerance": "-1"}
    ]
    assert details["invalid_molar_mass_range"] == {
        "minimum": 62, "maximum": 58
    }


def test_http_whitespace_trimmed_before_validation():
    resp = client.post(PATH, json=body([
        (" 75 ", " 0 "), (" 25", "0"), ("0", "0"), ("0", "0"),
    ], 16, 16))
    assert resp.status_code == 200, resp.text
    assert resp.json()["matches"][0]["formula"] == "CH4"
    assert resp.json()["elements"][0]["percent"] == "75"


# ------------------------------------------------------------ HTTP：结构校验（422 INVALID_REQUEST）
@pytest.mark.parametrize("payload", [
    {},
    {"molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
        element_block("C", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        element_block("X", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        {"element": "C", "percent": 0, "tolerance": "0"},
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        {"element": "C", "percent": "  ", "tolerance": "0"},
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": 100}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": 0, "maximum": 100}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": 1, "maximum": -1}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ], "molar_mass": {"minimum": True, "maximum": 100}},
    {"elements": [
        element_block("C", "0", "0"),
        element_block("H", "0", "0"),
        element_block("N", "0", "0"),
        element_block("O", "0", "0"),
    ]},
])
def test_http_structural_validation_errors(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422, resp.text
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"
