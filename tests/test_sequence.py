"""POST /api/v1/stoichiometry/sequence 的 HTTP API 测试（连续反应推演）。

主链路选用其它测试文件未提交过的反应（FeS2 焙烧制硫酸、BaSO4 沉淀、
ZnS 焙烧等），避免幂等缓存造成跨文件的顺序依赖。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/stoichiometry/sequence"


def _post(steps, inventory):
    return client.post(
        PATH,
        json={
            "steps": [
                {"reactants": r, "products": p} for r, p in steps
            ],
            "initial_inventory": [
                {"species": s, "moles": m} for s, m in inventory
            ],
        },
    )


# ------------------------------------------------------------ 成功路径
def test_three_step_chain_propagates_inventory():
    # 4FeS2 + 11O2 -> 2Fe2O3 + 8SO2
    # O2 + 2SO2 -> 2SO3
    # H2O + SO3 -> H2SO4
    # 初始 FeS2=4、O2=20、H2O=5：
    #   步0 ξ=1（FeS2 限量）-> O2 剩 9、SO2=8
    #   步1 ξ=4（SO2 限量）-> O2 剩 5、SO3=8
    #   步2 ξ=5（H2O 限量）-> H2O 耗尽、SO3 剩 3、H2SO4=5
    resp = _post(
        [
            (["FeS2", "O2"], ["Fe2O3", "SO2"]),
            (["SO2", "O2"], ["SO3"]),
            (["H2O", "SO3"], ["H2SO4"]),
        ],
        [("FeS2", "4"), ("O2", "20"), ("H2O", "5")],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert [s["index"] for s in data["steps"]] == [0, 1, 2]

    s0, s1, s2 = data["steps"]
    assert s0["normalized"]["equation"] == "4FeS2 + 11O2 -> 2Fe2O3 + 8SO2"
    assert s0["reaction_extent"] == {"value": "1", "unit": "mol"}
    assert s0["blocked"] is False
    assert s0["blocking_reactants"] == []
    assert s0["cached"] is False  # 本测试进程首次提交该反应

    r0 = {x["normalized"]: x for x in s0["reactants"]}
    assert (r0["FeS2"]["before"], r0["FeS2"]["consumed"],
            r0["FeS2"]["after"], r0["FeS2"]["is_blocking"]) == ("4", "4", "0", False)
    assert (r0["O2"]["before"], r0["O2"]["consumed"], r0["O2"]["after"]) == (
        "20", "11", "9",
    )
    p0 = {x["normalized"]: x for x in s0["products"]}
    assert p0["SO2"]["before"] == "0"
    assert (p0["SO2"]["produced"], p0["SO2"]["after"]) == ("8", "8")

    assert s1["normalized"]["equation"] == "O2 + 2SO2 -> 2SO3"
    assert s1["reaction_extent"]["value"] == "4"
    r1 = {x["normalized"]: x for x in s1["reactants"]}
    assert (r1["SO2"]["before"], r1["SO2"]["consumed"], r1["SO2"]["after"]) == (
        "8", "8", "0",
    )
    assert (r1["O2"]["before"], r1["O2"]["consumed"], r1["O2"]["after"]) == (
        "9", "4", "5",
    )

    assert s2["reaction_extent"]["value"] == "5"
    r2 = {x["normalized"]: x for x in s2["reactants"]}
    assert (r2["H2O"]["before"], r2["H2O"]["after"]) == ("5", "0")
    assert (r2["SO3"]["before"], r2["SO3"]["consumed"], r2["SO3"]["after"]) == (
        "8", "5", "3",
    )

    def inv(block):
        return {item["normalized"]: item["moles"] for item in block}

    assert inv(s0["inventory"])["O2"] == "9"
    assert inv(s1["inventory"])["SO3"] == "8"
    assert inv(s2["inventory"]) == inv(data["final_inventory"])
    assert inv(data["final_inventory"]) == {
        "Fe2O3": "2",
        "FeS2": "0",
        "H2O": "0",
        "H2SO4": "5",
        "O2": "5",
        "SO2": "0",
        "SO3": "3",
    }
    # 库存按规范化化学式排序
    assert [i["normalized"] for i in data["final_inventory"]] == sorted(
        i["normalized"] for i in data["final_inventory"]
    )
    # 所有数量字段都是字符串
    assert isinstance(data["final_inventory"][0]["moles"], str)


def test_intermediate_feeds_next_step_with_alias_notation():
    # 等价写法不影响配平 ID 与库存传递
    resp = _post(
        [
            (["ZnS", "O2"], ["ZnO", "SO2"]),   # 2ZnS + 3O2 -> 2ZnO + 2SO2
            (["O2", "S O2"], ["[SO3]"]),       # 等价写法的 SO2/SO3
        ],
        [("Zn S", "2"), ("O2", "5")],
    )
    assert resp.status_code == 200
    s0, s1 = resp.json()["steps"]
    assert s0["reaction_extent"]["value"] == "1"
    # 步0 后 O2 = 5-3 = 2，SO2 = 2
    r1 = {x["normalized"]: x for x in s1["reactants"]}
    assert (r1["O2"]["before"], r1["SO2"]["before"]) == ("2", "2")
    # 步1：O2 + 2SO2 -> 2SO3，ξ = min(2, 2/2) = 1（SO2 限量）
    assert s1["reaction_extent"]["value"] == "1"
    assert (r1["SO2"]["after"], r1["O2"]["after"]) == ("0", "1")


def test_missing_reactant_means_zero_and_blocks_step():
    # 步0 库存里没有 SO2（O2 有 3 mol）-> ξ=0，SO2 阻塞，库存原样传递；
    # 步1 与本步无关，仍可正常进行，证明阻塞后推演继续。
    resp = _post(
        [
            (["SO2", "O2"], ["SO3"]),
            (["H2", "O2"], ["H2O"]),  # 2H2 + O2 -> 2H2O
        ],
        [("H2", "2"), ("O2", "3")],
    )
    assert resp.status_code == 200
    s0, s1 = resp.json()["steps"]
    assert s0["blocked"] is True
    assert s0["reaction_extent"]["value"] == "0"
    assert s0["blocking_reactants"] == ["SO2"]
    assert all(x["consumed"] == "0" for x in s0["reactants"])
    assert all(x["produced"] == "0" for x in s0["products"])
    blocking = {x["normalized"]: x["is_blocking"] for x in s0["reactants"]}
    assert blocking == {"SO2": True, "O2": False}

    # 阻塞步不消耗库存：步1 看到的仍是初始 O2=3
    r1 = {x["normalized"]: x for x in s1["reactants"]}
    assert r1["O2"]["before"] == "3"
    assert s1["blocked"] is False
    assert s1["reaction_extent"]["value"] == "1"  # H2 限量
    final = {i["normalized"]: i["moles"] for i in resp.json()["final_inventory"]}
    assert final["H2"] == "0"
    assert final["O2"] == "2"
    assert final["H2O"] == "2"


def test_partial_shortage_blocks_on_zero_stock_only():
    # O2 有库存但 SO2 完全没有：只有 SO2 是阻塞物，且库存不被改动
    resp = _post(
        [(["SO2", "O2"], ["SO3"])],
        [("O2", "7")],
    )
    data = resp.json()
    s0 = data["steps"][0]
    assert s0["blocked"] is True
    assert s0["blocking_reactants"] == ["SO2"]
    r0 = {x["normalized"]: x for x in s0["reactants"]}
    assert r0["SO2"]["is_blocking"] is True
    assert r0["O2"]["is_blocking"] is False
    assert r0["O2"]["before"] == r0["O2"]["after"] == "7"
    assert r0["O2"]["consumed"] == "0"


def test_zero_initial_inventory_blocks_everything_but_runs():
    resp = _post(
        [
            (["H2", "O2"], ["H2O"]),
            (["H2O"], ["H2", "O2"]),  # 2H2O -> 2H2 + O2；同样无库存可分解
        ],
        [("H2", "0"), ("O2", "0")],
    )
    assert resp.status_code == 200
    steps = resp.json()["steps"]
    assert all(s["blocked"] for s in steps)
    assert [s["reaction_extent"]["value"] for s in steps] == ["0", "0"]
    assert steps[0]["blocking_reactants"] == ["H2", "O2"]
    assert steps[1]["blocking_reactants"] == ["H2O"]


def test_equivalent_species_in_inventory_are_merged():
    # Ba2+ + SO4^2- -> BaSO4；SO4^2- 以两种等价写法给出 0.1 + 0.2
    resp = _post(
        [(["Ba2+", "SO4^2-"], ["BaSO4"])],
        [("Ba2+", "0.3"), ("SO4^2-", "0.1"), ("[SO4]2-", "0.2")],
    )
    assert resp.status_code == 200
    step = resp.json()["steps"][0]
    sulfate = next(x for x in step["reactants"] if x["normalized"] == "SO4^2-")
    assert sulfate["before"] == "0.3"  # 精确合并，无 0.30000…004 浮点尾巴
    assert sulfate["consumed"] == "0.3"
    assert sulfate["after"] == "0"
    barium = next(x for x in step["reactants"] if x["normalized"] == "Ba^2+")
    assert barium["after"] == "0"
    product = step["products"][0]
    assert product["normalized"] == "BaSO4"
    assert product["produced"] == "0.3"


def test_exact_rational_progress_propagates_without_rounding():
    # 2C2H6 + 7O2 -> 4CO2 + 6H2O；各 1 mol -> ξ=1/7
    # CO2 + C -> 2CO；库存 C=1 -> ξ=min(4/7, 1)=4/7
    # 全部以精确有理数在步骤间传递，输出按 28 位有效数字舍入
    resp = _post(
        [
            (["C2H6", "O2"], ["CO2", "H2O"]),
            (["CO2", "C"], ["CO"]),
        ],
        [("C2H6", "1"), ("O2", "1"), ("C", "1")],
    )
    assert resp.status_code == 200
    s0, s1 = resp.json()["steps"]
    assert s0["reaction_extent"]["value"] == (
        "0.1428571428571428571428571429"
    )
    p0 = {x["normalized"]: x for x in s0["products"]}
    assert p0["CO2"]["after"] == "0.5714285714285714285714285714"
    # 步1 ξ=4/7；消耗 CO2 4/7，剩余 0；C 剩 3/7；CO 产出 8/7
    assert s1["reaction_extent"]["value"] == (
        "0.5714285714285714285714285714"
    )
    final = {i["normalized"]: i["moles"] for i in resp.json()["final_inventory"]}
    assert final["CO2"] == "0"
    assert final["C"] == "0.4285714285714285714285714286"  # 3/7
    assert final["CO"] == "1.142857142857142857142857143"  # 8/7


def test_step_reaction_ids_match_balance_endpoint():
    chain = _post(
        [
            (["FeS2", "O2"], ["Fe2O3", "SO2"]),
            (["SO2", "O2"], ["SO3"]),
        ],
        [("FeS2", "4"), ("O2", "20")],
    ).json()
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["SO2", "O2"], "products": ["SO3"]},
    ).json()
    assert chain["steps"][1]["reaction_id"] == balance["reaction_id"]


def test_second_submission_marks_steps_cached():
    body = {
        "steps": [
            {"reactants": ["FeS2", "O2"], "products": ["Fe2O3", "SO2"]},
            {"reactants": ["SO2", "O2"], "products": ["SO3"]},
        ],
        "initial_inventory": [
            {"species": "FeS2", "moles": "4"},
            {"species": "O2", "moles": "20"},
        ],
    }
    client.post(PATH, json=body)
    second = client.post(PATH, json=body).json()
    assert all(s["cached"] is True for s in second["steps"])


def test_inventory_independent_between_requests():
    # 同一反应连续提交两次，结果必须完全一致，不残留上一请求的库存
    body = {
        "steps": [{"reactants": ["H2", "O2"], "products": ["H2O"]}],
        "initial_inventory": [
            {"species": "H2", "moles": "2"},
            {"species": "O2", "moles": "1"},
        ],
    }
    a = client.post(PATH, json=body).json()
    b = client.post(PATH, json=body).json()
    assert a["final_inventory"] == b["final_inventory"]


def test_product_of_earlier_step_is_reactant_with_prior_stock():
    # 产物在初始库存中已有存量时，步前库存 = 初始 + 此前各步累计
    resp = _post(
        [
            (["H2", "O2"], ["H2O"]),  # 2H2 + O2 -> 2H2O
            (["H2O", "Na"], ["NaOH", "H2"]),  # 2H2O + 2Na -> 2NaOH + H2
        ],
        [("H2", "2"), ("O2", "10"), ("H2O", "4"), ("Na", "10")],
    )
    # 步0 ξ=1 -> 再产 2 mol H2O，步1 步前 H2O = 6；ξ=min(6/2,10/2)=3
    s1 = resp.json()["steps"][1]
    water = next(x for x in s1["reactants"] if x["normalized"] == "H2O")
    assert water["before"] == "6"
    assert s1["reaction_extent"]["value"] == "3"


# ------------------------------------------------------------ 结构化错误
def test_step_formula_parse_error_reports_step_index():
    resp = _post(
        [
            (["H2", "O2"], ["H2O"]),
            (["Zz"], ["H2"]),
        ],
        [],
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["step"] == 1


def test_first_step_parse_error_step_index_zero():
    # 初始库存物种均合法，错误来自第 0 步反应式
    resp = _post([(["Zz"], ["H2"])], [("H2", "1")])
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["step"] == 0


def test_step_duplicate_species_conflict():
    resp = _post([(["H2", "H2"], ["H2O"])], [])
    assert resp.status_code == 409
    body = resp.json()
    assert body["error"]["code"] == "DUPLICATE_SPECIES"
    assert body["error"]["details"]["step"] == 0


def test_step_non_unique_balance_reuses_error_with_step():
    resp = _post(
        [
            (["H2", "O2"], ["H2O"]),
            (["C", "O2"], ["CO", "CO2"]),
        ],
        [("H2", "1"), ("O2", "1")],
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "NON_UNIQUE_BALANCE"
    assert body["error"]["details"]["step"] == 1
    assert body["error"]["details"]["nullity"] == 2


def test_inventory_species_parse_error_is_400():
    resp = _post(
        [(["H2", "O2"], ["H2O"])],
        [("H2", "1"), ("Xy", "1")],
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["field"] == "initial_inventory.species"
    assert body["error"]["details"]["input"] == "Xy"


@pytest.mark.parametrize("bad", ["-1", "-0.0", "abc", "NaN", "Infinity", "1e10001"])
def test_invalid_inventory_moles_is_invalid_request(bad):
    resp = _post(
        [(["H2", "O2"], ["H2O"])],
        [("H2", bad), ("O2", "1")],
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"
    assert body["error"]["details"]["invalid_moles"] == [
        {"species": "H2", "moles": bad}
    ]


def test_inventory_moles_must_be_string_not_json_number():
    resp = client.post(
        PATH,
        json={
            "steps": [{"reactants": ["H2", "O2"], "products": ["H2O"]}],
            "initial_inventory": [{"species": "H2", "moles": 1}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


@pytest.mark.parametrize(
    "payload",
    [
        # 缺 steps / initial_inventory
        {"initial_inventory": []},
        {"steps": []},
        # steps 必须非空
        {"steps": [], "initial_inventory": []},
        # 每步仍需完整的 reactants / products
        {"steps": [{"reactants": ["H2"]}], "initial_inventory": []},
        {"steps": [{"products": ["H2O"]}], "initial_inventory": []},
        {"steps": [{"reactants": [], "products": ["H2O"]}], "initial_inventory": []},
        # 库存条目不能为空白
        {"steps": [{"reactants": ["H2"], "products": ["H2"]}],
         "initial_inventory": [{"species": "  ", "moles": "1"}]},
    ],
)
def test_request_shape_validation(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_response_shape_regression_guard():
    resp = _post(
        [
            (["H2", "O2"], ["H2O"]),
            (["H2O"], ["H2", "O2"]),
        ],
        [("H2", "4"), ("O2", "2")],
    )
    data = resp.json()
    assert set(data) == {"steps", "final_inventory"}
    step = data["steps"][0]
    assert set(step) == {
        "index", "reaction_id", "cached", "normalized", "reaction_extent",
        "blocked", "blocking_reactants", "reactants", "products", "inventory",
    }
    assert set(step["reactants"][0]) == {
        "input", "normalized", "coefficient",
        "before", "consumed", "after", "is_blocking",
    }
    assert set(step["products"][0]) == {
        "input", "normalized", "coefficient",
        "before", "produced", "after",
    }
    assert set(data["final_inventory"][0]) == {"normalized", "moles"}
