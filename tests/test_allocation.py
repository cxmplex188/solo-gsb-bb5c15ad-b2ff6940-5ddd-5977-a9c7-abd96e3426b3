"""POST /api/v1/inventory/allocate 的 HTTP API 测试。

覆盖：FEFO 排序（失效日期优先、同日批号升序）、名义量×纯度换算、
跨批次拆分、过期/隔离/零可用量跳过、未满足需求与汇总状态、未匹配物种
批次、核算日当天不算过期、等价化学式合并、精确十进制（无限小数）输出，
以及 4xx 输入校验（重复批号/需求、非法纯度/摩尔数、未知状态、坏日期、
坏化学式、空列表、JSON 数字等）。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/inventory/allocate"
AS_OF = "2026-09-23"


def batch(no, species, purity, nominal, expires, status="available"):
    return {
        "batch_no": no,
        "species": species,
        "purity": purity,
        "nominal_moles": nominal,
        "expires_on": expires,
        "status": status,
    }


def body(demands, batches_, as_of=AS_OF):
    return {"demands": demands, "batches": batches_, "accounting_date": as_of}


def demand(species, moles):
    return {"species": species, "effective_moles": moles}


def lines_by_no(data, species="HCl"):
    alloc = next(a for a in data["allocations"] if a["normalized"] == species)
    return {line["batch_no"]: line for line in alloc["batches"]}, alloc


# ------------------------------------------------------------ 成功路径
def test_fefo_order_and_full_allocation():
    # 同日两批按批号升序：B-001 先耗尽再用 B-002
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "2")],
            [
                batch("B-002", "HCl", "1", "2", "2026-10-01"),
                batch("B-001", "HCl", "1", "1", "2026-10-01"),
            ],
        ),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    alloc = data["allocations"][0]
    assert [line["batch_no"] for line in alloc["batches"]] == ["B-001", "B-002"]
    by_no = {line["batch_no"]: line for line in alloc["batches"]}
    assert by_no["B-001"]["used"] is True
    assert by_no["B-001"]["allocated_effective_moles"] == "1"
    assert by_no["B-001"]["remaining_effective_moles"] == "0"
    assert by_no["B-002"]["allocated_effective_moles"] == "1"
    assert by_no["B-002"]["remaining_effective_moles"] == "1"
    assert alloc["status"] == "fulfilled"
    assert alloc["unmet_effective_moles"] == "0"
    assert data["summary"]["status"] == "fulfilled"
    assert data["summary"]["batches_used"] == 2
    assert data["as_of"] == AS_OF


def test_expiry_date_takes_priority_over_request_order():
    # 请求顺序故意打乱：2027 年的批在前，也必须后用
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "2")],
            [
                batch("Z-LATE", "HCl", "1", "5", "2027-01-01"),
                batch("A-EARLY", "HCl", "1", "2", "2026-10-01"),
            ],
        ),
    )
    data = resp.json()
    assert [line["batch_no"] for line in data["allocations"][0]["batches"]] == [
        "A-EARLY",
        "Z-LATE",
    ]
    by_no, _ = lines_by_no(data)
    assert by_no["A-EARLY"]["allocated_effective_moles"] == "2"
    assert by_no["Z-LATE"]["allocated_effective_moles"] == "0"
    assert by_no["Z-LATE"]["used"] is False
    assert by_no["Z-LATE"]["skipped"] is False


def test_expired_batch_is_skipped_and_not_deducted():
    # expires_on 早于核算日一天 -> 过期跳过，仅靠另一批满足
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("OLD", "HCl", "1", "9", "2026-09-22"),
                batch("NEW", "HCl", "1", "1", "2026-10-01"),
            ],
        ),
    )
    data = resp.json()
    by_no, alloc = lines_by_no(data)
    old = by_no["OLD"]
    assert old["skipped"] is True
    assert old["skip_reason"] == "expired"
    assert old["eligible"] is False
    assert old["used"] is False
    assert old["allocated_nominal_moles"] == "0"
    assert old["remaining_nominal_moles"] == "9"  # 整批保留
    assert by_no["NEW"]["allocated_effective_moles"] == "1"
    assert data["summary"]["batches_skipped"] == 1


def test_expires_on_accounting_date_is_not_expired():
    # 失效日期恰为核算日：仍可领用
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [batch("EDGE", "HCl", "1", "1", AS_OF)],
        ),
    )
    data = resp.json()
    line = data["allocations"][0]["batches"][0]
    assert line["skip_reason"] is None
    assert line["used"] is True
    assert line["allocated_effective_moles"] == "1"


def test_quarantined_batch_is_skipped_even_if_earliest():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("Q", "HCl", "1", "5", "2026-09-30", "quarantined"),
                batch("OK", "HCl", "1", "1", "2027-01-01"),
            ],
        ),
    )
    data = resp.json()
    by_no, _ = lines_by_no(data)
    q = by_no["Q"]
    assert q["skipped"] is True
    assert q["skip_reason"] == "quarantined"
    assert q["allocated_effective_moles"] == "0"
    assert q["remaining_effective_moles"] == "5"
    assert by_no["OK"]["allocated_effective_moles"] == "1"


def test_quarantine_takes_precedence_over_expiry():
    # 既隔离又过期：报 quarantined（状态先判定）
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("BOTH", "HCl", "1", "1", "2020-01-01", "quarantined"),
                batch("OK", "HCl", "1", "1", "2027-01-01"),
            ],
        ),
    )
    by_no, _ = lines_by_no(resp.json())
    assert by_no["BOTH"]["skip_reason"] == "quarantined"


def test_purity_conversion_nominal_to_effective():
    # 纯度 0.98、名义 50 -> 有效 49；需求 25 -> 名义领用 25/0.98
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "25")],
            [batch("P98", "HCl", "0.98", "50", "2026-12-01")],
        ),
    )
    line = resp.json()["allocations"][0]["batches"][0]
    assert line["effective_moles"] == "49"
    assert line["allocated_effective_moles"] == "25"
    # 名义领用 = 25 / 0.98 = 25.5102040816326530612244898（28 位有效数字）
    assert line["allocated_nominal_moles"] == "25.5102040816326530612244898"
    assert line["remaining_nominal_moles"] == "24.4897959183673469387755102"
    assert line["remaining_effective_moles"] == "24"


def test_cross_batch_split_with_mixed_purities():
    # 需求有效 2.5：批1 有效 1（耗尽），批2 纯度0.5 再领有效 1.5 -> 名义 3
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "2.5")],
            [
                batch("B1", "HCl", "1", "1", "2026-10-01"),
                batch("B2", "HCl", "0.5", "6", "2026-11-01"),
            ],
        ),
    )
    by_no, alloc = lines_by_no(resp.json())
    assert by_no["B1"]["allocated_nominal_moles"] == "1"
    assert by_no["B1"]["remaining_nominal_moles"] == "0"
    assert by_no["B2"]["allocated_effective_moles"] == "1.5"
    assert by_no["B2"]["allocated_nominal_moles"] == "3"
    assert by_no["B2"]["remaining_nominal_moles"] == "3"
    assert by_no["B2"]["remaining_effective_moles"] == "1.5"
    assert alloc["allocated_effective_moles"] == "2.5"
    assert alloc["status"] == "fulfilled"


def test_shortage_when_effective_capacity_insufficient():
    # 有效合计 1.5（1 + 0.5）< 需求 2 -> 全部领完，缺 0.5
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "2")],
            [
                batch("B1", "HCl", "1", "1", "2026-10-01"),
                batch("B2", "HCl", "0.5", "1", "2026-11-01"),
            ],
        ),
    )
    data = resp.json()
    alloc = data["allocations"][0]
    assert alloc["status"] == "shortage"
    assert alloc["allocated_effective_moles"] == "1.5"
    assert alloc["unmet_effective_moles"] == "0.5"
    summary = data["summary"]
    assert summary["status"] == "shortage"
    assert summary["fulfilled_demands"] == 0
    assert summary["shortage_demands"] == 1
    assert summary["total_unmet_effective_moles"] == "0.5"


def test_zero_nominal_batch_is_skipped_as_exhausted():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("EMPTY", "HCl", "1", "0", "2026-10-01"),
                batch("FULL", "HCl", "1", "1", "2026-11-01"),
            ],
        ),
    )
    data = resp.json()
    by_no, _ = lines_by_no(data)
    empty = by_no["EMPTY"]
    assert empty["skipped"] is True
    assert empty["skip_reason"] == "exhausted"
    assert empty["used"] is False
    assert data["summary"]["batches_skipped"] == 1
    assert data["summary"]["batches_used"] == 1


def test_unmatched_species_batches_are_returned_separately():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("H1", "HCl", "1", "1", "2026-10-01"),
                batch("W1", "H2O", "1", "3", "2026-10-01"),
                batch("N1", "HNO3", "0.7", "2", "2026-12-01"),
            ],
        ),
    )
    data = resp.json()
    unmatched = data["unmatched_batches"]
    assert [line["batch_no"] for line in unmatched] == ["W1", "N1"]  # 按规范化物种排序
    assert all(line["skip_reason"] == "unmatched_species" for line in unmatched)
    assert all(line["allocated_nominal_moles"] == "0" for line in unmatched)
    w1 = next(line for line in unmatched if line["batch_no"] == "W1")
    assert w1["remaining_nominal_moles"] == "3"
    assert data["summary"]["unmatched_batches"] == 2
    # 未匹配批次不计入 batches_skipped（单独计数）
    assert data["summary"]["batches_skipped"] == 0
    assert data["summary"]["batches_total"] == 3


def test_equivalent_formula_writings_match_same_species():
    # 需求 "Ca(OH)2"，批次写作 "Ca [OH]2" -> 同一物种，可分配
    resp = client.post(
        PATH,
        json=body(
            [demand("Ca(OH)2", "1")],
            [batch("C1", "Ca [OH]2", "1", "1", "2026-10-01")],
        ),
    )
    assert resp.status_code == 200, resp.text
    alloc = resp.json()["allocations"][0]
    assert alloc["normalized"] == "Ca(OH)2"
    assert alloc["allocated_effective_moles"] == "1"


def test_multiple_species_independent_allocation_and_sorting():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1"), demand("NaOH", "5")],
            [
                batch("H1", "HCl", "1", "1", "2026-10-01"),
                batch("N1", "NaOH", "1", "2", "2026-10-01"),  # NaOH 只有 2，缺 3
            ],
        ),
    )
    data = resp.json()
    # 物种按规范化化学式排序
    assert [a["normalized"] for a in data["allocations"]] == ["HCl", "NaOH"]
    statuses = {a["normalized"]: a["status"] for a in data["allocations"]}
    assert statuses == {"HCl": "fulfilled", "NaOH": "shortage"}
    summary = data["summary"]
    assert summary["demands_total"] == 2
    assert summary["fulfilled_demands"] == 1
    assert summary["shortage_demands"] == 1
    assert summary["total_requested_effective_moles"] == "6"
    assert summary["total_allocated_effective_moles"] == "3"
    assert summary["total_unmet_effective_moles"] == "3"
    assert summary["status"] == "shortage"


def test_zero_demand_allocates_nothing_and_is_fulfilled():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "0")],
            [batch("B1", "HCl", "1", "5", "2026-10-01")],
        ),
    )
    data = resp.json()
    alloc = data["allocations"][0]
    assert alloc["status"] == "fulfilled"
    assert alloc["allocated_effective_moles"] == "0"
    line = alloc["batches"][0]
    # 无需求：批次合格但不扣减（不算 skipped、不算 used）
    assert line["eligible"] is True
    assert line["skipped"] is False
    assert line["used"] is False
    assert line["remaining_nominal_moles"] == "5"
    assert data["summary"]["batches_used"] == 0


def test_empty_batches_yields_shortage():
    resp = client.post(
        PATH, json=body([demand("HCl", "1")], [])
    )
    assert resp.status_code == 200
    data = resp.json()
    alloc = data["allocations"][0]
    assert alloc["batches"] == []
    assert alloc["status"] == "shortage"
    assert alloc["unmet_effective_moles"] == "1"
    assert data["summary"]["status"] == "shortage"


def test_exact_repeating_decimal_rendered_as_fixed_point_string():
    # 需求 1，纯度 1/3 无法直接表达；用名义 3、纯度 0.333...? 用可构造的
    # 无限小数：批纯度 0.7（=7/10），需求 1 -> 名义领用 1/0.7 = 10/7
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [batch("B1", "HCl", "0.7", "10", "2026-10-01")],
        ),
    )
    line = resp.json()["allocations"][0]["batches"][0]
    assert line["allocated_nominal_moles"] == "1.428571428571428571428571429"
    assert isinstance(line["allocated_nominal_moles"], str)


def test_each_batch_deducted_only_once():
    # 同一批号在输入中只允许出现一次 -> 合法请求中每批最多扣减一次；
    # 两个需求指向同一物种是非法的（下面的重复用例覆盖），这里用两个不同物种
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "3"), demand("H2O", "1")],
            [
                batch("B1", "HCl", "1", "2", "2026-10-01"),
                batch("B2", "H2O", "1", "1", "2026-10-01"),
            ],
        ),
    )
    data = resp.json()
    hcl = next(a for a in data["allocations"] if a["normalized"] == "HCl")
    # B1 只在 HCl 下出现一次，扣 2，缺 1
    assert len(hcl["batches"]) == 1
    assert hcl["batches"][0]["allocated_effective_moles"] == "2"
    assert hcl["unmet_effective_moles"] == "1"


# ------------------------------------------------------------ 4xx 校验
def test_duplicate_batch_number_rejected():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("DUP", "HCl", "1", "1", "2026-10-01"),
                batch("DUP", "HCl", "1", "1", "2026-11-01"),
            ],
        ),
    )
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_INVENTORY"
    assert err["details"]["duplicate_batch_nos"] == ["DUP"]


def test_duplicate_demand_species_rejected():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1"), demand(" HCl ", "1")],
            [batch("B1", "HCl", "1", "5", "2026-10-01")],
        ),
    )
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["code"] == "INVALID_INVENTORY"
    assert err["details"]["duplicate_demand_species"] == ["HCl"]


def test_invalid_purity_values_rejected():
    for bad in ["0", "1.01", "-0.5", "abc", "NaN", "Infinity"]:
        resp = client.post(
            PATH,
            json=body(
                [demand("HCl", "1")],
                [batch("B1", "HCl", bad, "5", "2026-10-01")],
            ),
        )
        assert resp.status_code == 422, bad
        err = resp.json()["error"]
        assert err["code"] == "INVALID_INVENTORY"
        assert err["details"]["invalid_purity"][0]["purity"] == bad


def test_invalid_moles_values_rejected():
    # 需求有效摩尔非法
    resp = client.post(
        PATH,
        json=body([demand("HCl", "-1")], [batch("B1", "HCl", "1", "5", "2026-10-01")]),
    )
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["details"]["invalid_effective_moles"][0]["effective_moles"] == "-1"

    # 批次名义摩尔非法
    resp = client.post(
        PATH,
        json=body([demand("HCl", "1")], [batch("B1", "HCl", "1", "abc", "2026-10-01")]),
    )
    assert resp.status_code == 422
    err = resp.json()["error"]
    assert err["details"]["invalid_nominal_moles"][0]["nominal_moles"] == "abc"


def test_json_number_instead_of_decimal_string_rejected():
    # Pydantic 层：数字不是字符串 -> 422 INVALID_REQUEST
    payload = body([demand("HCl", 1)], [])  # effective_moles 传成数字
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    payload = body(
        [demand("HCl", "1")],
        [{**batch("B1", "HCl", "1", "5", "2026-10-01"), "purity": 0.9}],
    )
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unknown_status_rejected():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                {
                    **batch("B1", "HCl", "1", "5", "2026-10-01"),
                    "status": "released",
                }
            ],
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_invalid_date_rejected():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [batch("B1", "HCl", "1", "5", "10/01/2026")],
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    resp = client.post(PATH, json=body([demand("HCl", "1")], [], "2026-13-01"))
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unparseable_demand_species_is_400():
    resp = client.post(
        PATH,
        json=body([demand("Xx", "1")], [batch("B1", "HCl", "1", "5", "2026-10-01")]),
    )
    assert resp.status_code == 400
    err = resp.json()["error"]
    assert err["code"] == "FORMULA_PARSE_ERROR"
    assert err["details"]["field"] == "demands.species"


def test_unparseable_batch_species_is_400():
    resp = client.post(
        PATH,
        json=body([demand("HCl", "1")], [batch("B1", "Yy", "1", "5", "2026-10-01")]),
    )
    assert resp.status_code == 400
    err = resp.json()["error"]
    assert err["code"] == "FORMULA_PARSE_ERROR"
    assert err["details"]["field"] == "batches.species"
    assert err["details"]["batch_no"] == "B1"


@pytest.mark.parametrize("payload", [
    {"batches": [], "accounting_date": AS_OF},                      # 缺 demands
    {"demands": [], "batches": [], "accounting_date": AS_OF},       # demands 为空
    {"demands": [demand("HCl", "1")], "accounting_date": AS_OF},    # 缺 batches
    {"demands": [demand("HCl", "1")], "batches": []},               # 缺 accounting_date
    {},                                                             # 空请求体
])
def test_missing_or_empty_fields_rejected(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_blank_string_fields_rejected():
    resp = client.post(
        PATH,
        json=body(
            [{"species": "  ", "effective_moles": "1"}],
            [],
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [batch("  ", "HCl", "1", "5", "2026-10-01")],
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_non_json_body_rejected():
    resp = client.post(PATH, content="not json", headers={"Content-Type": "application/json"})
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_multiple_validation_issues_aggregated():
    # 同时：批号重复 + 纯度非法 + 名义量非法，一次响应全部列出
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "1")],
            [
                batch("DUP", "HCl", "2", "-3", "2026-10-01"),
                batch("DUP", "HCl", "0", "1", "2026-11-01"),
            ],
        ),
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate_batch_nos"] == ["DUP"]
    assert {item["batch_no"] for item in details["invalid_purity"]} == {"DUP"}
    assert {item["batch_no"] for item in details["invalid_nominal_moles"]} == {"DUP"}


def test_response_matches_pydantic_model_shape():
    resp = client.post(
        PATH,
        json=body(
            [demand("HCl", "2")],
            [
                batch("B1", "HCl", "0.98", "1", "2026-10-01"),
                batch("OLD", "HCl", "1", "9", "2020-01-01"),
            ],
        ),
    )
    assert resp.status_code == 200
    data = resp.json()
    assert set(data) == {"as_of", "allocations", "unmatched_batches", "summary"}
    alloc = data["allocations"][0]
    assert set(alloc) == {
        "species", "normalized", "requested_effective_moles",
        "allocated_effective_moles", "unmet_effective_moles",
        "status", "batches",
    }
    line = alloc["batches"][0]
    assert set(line) == {
        "batch_no", "species", "normalized", "expires_on", "status", "purity",
        "nominal_moles", "effective_moles", "eligible", "skipped", "skip_reason",
        "used", "allocated_nominal_moles", "allocated_effective_moles",
        "remaining_nominal_moles", "remaining_effective_moles",
    }
    assert set(data["summary"]) == {
        "total_requested_effective_moles", "total_allocated_effective_moles",
        "total_unmet_effective_moles", "demands_total", "fulfilled_demands",
        "shortage_demands", "batches_total", "batches_used", "batches_skipped",
        "unmatched_batches", "status",
    }
