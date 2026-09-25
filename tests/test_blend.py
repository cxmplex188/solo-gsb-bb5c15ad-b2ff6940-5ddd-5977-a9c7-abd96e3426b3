"""POST /api/v1/solutions/blend 的 HTTP API 测试。

覆盖：简单稀释、两母液两溶质精确求解、母液总量最小化、输入顺序字典序
tie-break、可用体积封顶、总量预算（余量补水）、零目标排斥母液、非目标
溶质报告、无限小数精确输出、无解（供应不足 / 总量预算 / 目标溶质缺失）
以及 4xx 输入校验（母液数量、标识/溶质重复、非法浓度/体积、零终体积、
JSON 数字、空白字符串等）。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/solutions/blend"


def stock(id_, solutes, available):
    return {
        "id": id_,
        "solutes": [{"solute": s, "concentration": c} for s, c in solutes],
        "available_volume": available,
    }


def target(solute, concentration):
    return {"solute": solute, "concentration": concentration}


def body(stocks, targets, target_volume="100"):
    return {"stocks": stocks, "target_volume": target_volume, "targets": targets}


# ------------------------------------------------------------ 成功路径
def test_simple_dilution():
    # 1.0 mol/L 母液稀释到 0.1 mol/L、终体积 250 mL：取 25 mL，补水 225 mL
    resp = client.post(
        PATH,
        json=body(
            [stock("stock-A", [("NaCl", "1.0")], "100")],
            [target("NaCl", "0.1")],
            target_volume="250",
        ),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["target_volume"] == "250"
    assert data["total_stock_volume"] == "25"
    assert data["water_added"] == "225"
    assert data["stocks"] == [
        {
            "id": "stock-A",
            "available_volume": "100",
            "volume_used": "25",
            "remaining_volume": "75",
        }
    ]
    assert data["final_concentrations"] == [
        {
            "solute": "NaCl",
            "target_concentration": "0.1",
            "final_concentration": "0.1",
            "target_met": True,
        }
    ]


def test_two_stocks_two_solutes_exact():
    # x_A = 280/9、x_B = 400/9、补水 220/9（无限小数按 28 位有效数字输出）
    resp = client.post(
        PATH,
        json=body(
            [
                stock("A", [("NaCl", "1.0"), ("KCl", "0.5")], "100"),
                stock("B", [("NaCl", "0.2"), ("KCl", "1.0")], "100"),
            ],
            [target("NaCl", "0.4"), target("KCl", "0.6")],
        ),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    by_id = {s["id"]: s for s in data["stocks"]}
    assert by_id["A"]["volume_used"] == "31.11111111111111111111111111"
    assert by_id["B"]["volume_used"] == "44.44444444444444444444444444"
    assert data["water_added"] == "24.44444444444444444444444444"
    final = {s["solute"]: s for s in data["final_concentrations"]}
    assert final["NaCl"]["final_concentration"] == "0.4"
    assert final["KCl"]["final_concentration"] == "0.6"
    assert all(s["target_met"] is True for s in final.values())


def test_min_total_volume_prefers_concentrated_stock():
    # 2.0 mol/L 与 1.0 mol/L 都能单独满足：优先用浓母液使总量最小
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "2.0")], "50"), stock("B", [("NaCl", "1.0")], "100")],
            [target("NaCl", "0.5")],
        ),
    )
    data = resp.json()
    by_id = {s["id"]: s for s in data["stocks"]}
    assert by_id["A"]["volume_used"] == "25"
    assert by_id["B"]["volume_used"] == "0"
    assert data["total_stock_volume"] == "25"
    assert data["water_added"] == "75"


def test_lexicographic_tie_break_by_input_order():
    # 两种完全相同的母液：总量相同的多解中按输入顺序字典序最小 -> 第一种取 0
    resp = client.post(
        PATH,
        json=body(
            [stock("first", [("NaCl", "1.0")], "100"),
             stock("second", [("NaCl", "1.0")], "100")],
            [target("NaCl", "0.5")],
        ),
    )
    data = resp.json()
    assert [s["volume_used"] for s in data["stocks"]] == ["0", "50"]
    assert [s["id"] for s in data["stocks"]] == ["first", "second"]


def test_available_volume_cap_is_respected():
    # 浓母液可用体积只有 20 mL：封顶后由稀母液补足
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "2.0")], "20"), stock("B", [("NaCl", "1.0")], "100")],
            [target("NaCl", "0.5")],
        ),
    )
    data = resp.json()
    by_id = {s["id"]: s for s in data["stocks"]}
    assert by_id["A"]["volume_used"] == "20"
    assert by_id["A"]["remaining_volume"] == "0"
    assert by_id["B"]["volume_used"] == "10"
    assert data["water_added"] == "70"


def test_zero_target_excludes_stock_carrying_that_solute():
    # KCl 目标为 0：含 KCl 的母液 A 不能取用，只能由纯 NaCl 母液 B 提供
    resp = client.post(
        PATH,
        json=body(
            [
                stock("A", [("NaCl", "1.0"), ("KCl", "0.5")], "100"),
                stock("B", [("NaCl", "1.0")], "100"),
            ],
            [target("NaCl", "0.5"), target("KCl", "0")],
        ),
    )
    data = resp.json()
    by_id = {s["id"]: s for s in data["stocks"]}
    assert by_id["A"]["volume_used"] == "0"
    assert by_id["B"]["volume_used"] == "50"
    final = {s["solute"]: s for s in data["final_concentrations"]}
    assert final["KCl"]["final_concentration"] == "0"
    assert final["KCl"]["target_met"] is True


def test_non_target_solute_is_reported_with_null_target():
    # 母液带入的 KCl 不在目标中：不约束，但最终浓度如实报告
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1.0"), ("KCl", "0.25")], "100")],
            [target("NaCl", "0.5")],
        ),
    )
    data = resp.json()
    final = {s["solute"]: s for s in data["final_concentrations"]}
    assert final["KCl"] == {
        "solute": "KCl",
        "target_concentration": None,
        "final_concentration": "0.125",
        "target_met": None,
    }
    assert final["NaCl"]["final_concentration"] == "0.5"


def test_all_zero_targets_use_no_stock():
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "1.0")], "100")], [target("NaCl", "0")]),
    )
    data = resp.json()
    assert data["stocks"][0]["volume_used"] == "0"
    assert data["total_stock_volume"] == "0"
    assert data["water_added"] == "100"


def test_repeating_decimal_exact_rational_output():
    # 0.3 mol/L 母液配 0.1 mol/L、终体积 1 mL：取用量恰为 1/3 mL
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "0.3")], "100")],
            [target("NaCl", "0.1")],
            target_volume="1",
        ),
    )
    data = resp.json()
    assert data["stocks"][0]["volume_used"] == "0.3333333333333333333333333333"
    assert data["water_added"] == "0.6666666666666666666666666667"
    assert data["final_concentrations"][0]["final_concentration"] == "0.1"


def test_four_stocks_maximum():
    resp = client.post(
        PATH,
        json=body(
            [
                stock("A", [("NaCl", "4")], "10"),
                stock("B", [("NaCl", "3")], "10"),
                stock("C", [("NaCl", "2")], "10"),
                stock("D", [("NaCl", "1")], "10"),
            ],
            [target("NaCl", "0.5")],
        ),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    # 需 50 mmol：最浓的 A 提供 40，余 10 由 B 提供（10/3 mL）
    assert [s["volume_used"] for s in data["stocks"]] == [
        "10",
        "3.333333333333333333333333333",
        "0",
        "0",
    ]


def test_stock_without_solutes_is_never_used():
    # 纯溶剂母液等价于水，不会优于直接补水
    resp = client.post(
        PATH,
        json=body(
            [stock("water-stock", [], "500"), stock("A", [("NaCl", "1.0")], "100")],
            [target("NaCl", "0.5")],
        ),
    )
    data = resp.json()
    by_id = {s["id"]: s for s in data["stocks"]}
    assert by_id["water-stock"]["volume_used"] == "0"
    assert by_id["A"]["volume_used"] == "50"


def test_stateless_same_request_same_result():
    payload = body(
        [stock("A", [("NaCl", "1.0")], "100")], [target("NaCl", "0.5")]
    )
    first = client.post(PATH, json=payload)
    second = client.post(PATH, json=payload)
    assert first.status_code == second.status_code == 200
    assert first.json() == second.json()


def test_response_shape():
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "1.0")], "100")], [target("NaCl", "0.5")]),
    )
    data = resp.json()
    assert set(data) == {
        "target_volume",
        "total_stock_volume",
        "water_added",
        "stocks",
        "final_concentrations",
    }
    assert set(data["stocks"][0]) == {
        "id",
        "available_volume",
        "volume_used",
        "remaining_volume",
    }
    assert set(data["final_concentrations"][0]) == {
        "solute",
        "target_concentration",
        "final_concentration",
        "target_met",
    }


# ------------------------------------------------------------ 无解（422 NO_FEASIBLE_BLEND）
def test_infeasible_target_exceeds_total_supply():
    # 目标 1.0 mol/L × 100 mL = 100 mmol，母液全部可用量仅 0.5 × 100 = 50 mmol
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "0.5")], "100")], [target("NaCl", "1.0")]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "NO_FEASIBLE_BLEND"
    insufficient = body_["error"]["details"]["insufficient_solutes"]
    assert insufficient == [
        {"solute": "NaCl", "required_mmol": "100", "available_mmol": "50"}
    ]


def test_infeasible_volume_budget_only():
    # 单溶质供应充足（0.4 × 500 = 200 mmol ≥ 50），但所需 125 mL 超过终体积
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "0.4")], "500")], [target("NaCl", "0.5")]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "NO_FEASIBLE_BLEND"
    assert "insufficient_solutes" not in body_["error"]["details"]


def test_infeasible_target_solute_missing_from_all_stocks():
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "1.0")], "100")], [target("KCl", "0.1")]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "NO_FEASIBLE_BLEND"
    insufficient = body_["error"]["details"]["insufficient_solutes"]
    assert insufficient[0]["solute"] == "KCl"
    assert insufficient[0]["available_mmol"] == "0"


def test_infeasible_zero_target_conflicts_with_required_solute():
    # 只有含 KCl 的母液能提供 NaCl，但 KCl 目标为 0 -> 无解
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1.0"), ("KCl", "0.5")], "100")],
            [target("NaCl", "0.5"), target("KCl", "0")],
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NO_FEASIBLE_BLEND"


# ------------------------------------------------------------ 输入校验（422）
def test_too_many_stocks_rejected():
    payload = body(
        [stock(f"S{i}", [("NaCl", "1")], "10") for i in range(5)],
        [target("NaCl", "0.5")],
    )
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


@pytest.mark.parametrize("payload", [
    body([], [target("NaCl", "0.5")]),                      # 母液为空
    body([stock("A", [("NaCl", "1")], "10")], []),          # 目标为空
    body([stock(" ", [("NaCl", "1")], "10")], [target("NaCl", "0.5")]),  # 空白标识
    body([stock("A", [(" ", "1")], "10")], [target("NaCl", "0.5")]),     # 空白溶质
    body([stock("A", [("NaCl", "1")], " ")], [target("NaCl", "0.5")]),   # 空白体积
    body([stock("A", [("NaCl", 1.0)], "10")], [target("NaCl", "0.5")]),  # JSON 数字浓度
    body([stock("A", [("NaCl", "1")], 100)], [target("NaCl", "0.5")]),   # JSON 数字体积
    body([stock("A", [("NaCl", "1")], "10")], [target("NaCl", 0.5)]),    # JSON 数字目标
    body([stock("A", [("NaCl", "1")], "10")], [target("NaCl", "0.5")],
         target_volume=250),                                # JSON 数字终体积
])
def test_structural_validation_errors(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_duplicate_stock_id_rejected():
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1")], "10"), stock("A", [("KCl", "1")], "10")],
            [target("NaCl", "0.5")],
        ),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["duplicate_stock_ids"] == ["A"]


def test_duplicate_solute_within_stock_rejected():
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1"), ("NaCl", "2")], "10")],
            [target("NaCl", "0.5")],
        ),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["duplicate_stock_solutes"] == [
        {"stock_id": "A", "solute": "NaCl"}
    ]


def test_duplicate_target_solute_rejected():
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1")], "10")],
            [target("NaCl", "0.5"), target("NaCl", "0.6")],
        ),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["duplicate_target_solutes"] == ["NaCl"]


@pytest.mark.parametrize("bad", ["-1", "-0.0", "NaN", "Infinity", "abc", "1e20000"])
def test_invalid_concentration_values_rejected(bad):
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", bad)], "10")], [target("NaCl", "0.5")]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["invalid_concentrations"] == [
        {"stock_id": "A", "solute": "NaCl", "concentration": bad}
    ]


@pytest.mark.parametrize("bad", ["-1", "NaN", "abc"])
def test_invalid_available_volume_rejected(bad):
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "1")], bad)], [target("NaCl", "0.5")]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["invalid_available_volumes"] == [
        {"stock_id": "A", "available_volume": bad}
    ]


@pytest.mark.parametrize("bad", ["-0.5", "NaN", "abc"])
def test_invalid_target_concentration_rejected(bad):
    resp = client.post(
        PATH,
        json=body([stock("A", [("NaCl", "1")], "10")], [target("NaCl", bad)]),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["invalid_target_concentrations"] == [
        {"solute": "NaCl", "concentration": bad}
    ]


@pytest.mark.parametrize("bad", ["0", "0.0", "-1", "abc"])
def test_non_positive_target_volume_rejected(bad):
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1")], "10")],
            [target("NaCl", "0.5")],
            target_volume=bad,
        ),
    )
    assert resp.status_code == 422
    body_ = resp.json()
    assert body_["error"]["code"] == "INVALID_BLEND"
    assert body_["error"]["details"]["invalid_target_volume"] == bad


def test_multiple_input_problems_aggregated():
    resp = client.post(
        PATH,
        json=body(
            [
                stock("A", [("NaCl", "-1")], "10"),
                stock("A", [("KCl", "1")], "xyz"),
            ],
            [target("NaCl", "0.5"), target("NaCl", "0.6")],
            target_volume="0",
        ),
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate_stock_ids"] == ["A"]
    assert details["invalid_concentrations"] == [
        {"stock_id": "A", "solute": "NaCl", "concentration": "-1"}
    ]
    assert details["invalid_available_volumes"] == [
        {"stock_id": "A", "available_volume": "xyz"}
    ]
    assert details["duplicate_target_solutes"] == ["NaCl"]
    assert details["invalid_target_volume"] == "0"


def test_scientific_notation_and_zero_concentration_accepted():
    # 科学计数法与 0 浓度均为合法的非负十进制字符串
    resp = client.post(
        PATH,
        json=body(
            [stock("A", [("NaCl", "1e-3"), ("KCl", "0")], "1e3")],
            [target("NaCl", "5e-4")],
        ),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["stocks"][0]["volume_used"] == "50"
    final = {s["solute"]: s for s in data["final_concentrations"]}
    assert final["KCl"]["final_concentration"] == "0"
