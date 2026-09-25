"""POST /api/v1/equilibrium/solve 的 HTTP API 测试。

测试反应（PCl5 解离、2O3 ⇌ 3O2、NO2 二聚等）均选用其它测试文件未提交过的
方程式，避免幂等缓存造成跨文件的顺序依赖。数值核验以精确
:class:`fractions.Fraction` 在返回进度 ±5e-10 mol 处重新构造 Qc，
确认 Kc 落在该邻域内（等价于进度绝对误差 ≤ 5e-10 < 1e-9 mol）。
"""

from decimal import Decimal, localcontext
from fractions import Fraction

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/api/v1/equilibrium/solve"


def _body(reactants, products, moles, volume="1", kc="1"):
    return {
        "reactants": reactants,
        "products": products,
        "amounts": [{"species": s, "moles": m} for s, m in moles.items()],
        "volume_l": volume,
        "kc": kc,
    }


def _exact_qc(d, xi, n0, coef, reactants, products, volume):
    """在进度 xi（Fraction）处用精确有理数构造 Qc。"""
    num = Fraction(1)
    for s in products:
        num *= (n0[s] + coef[s] * xi) ** coef[s]
    den = Fraction(1)
    for s in reactants:
        den *= (n0[s] - coef[s] * xi) ** coef[s]
    delta_n = sum(coef[s] for s in products) - sum(coef[s] for s in reactants)
    return num / (den * volume ** delta_n)


def assert_solution_consistent(
    d, *, tol=Fraction(5, 10**10), qc_bracket=True
):
    """核验：进度落在物料区间、±tol 邻域夹住 Kc、各物种摩尔数非负且自洽。

    ``qc_bracket=False`` 用于近完全反应（进度贴近物料端点，Qc 对 ξ 极
    陡峭，固定绝对 tol 的邻域无法夹住 Kc），此时只核验进度容差、区间
    与摩尔数自洽。
    """
    species = d["reactants"] + d["products"]
    reactants = [x["normalized"] for x in d["reactants"]]
    products = [x["normalized"] for x in d["products"]]
    n0 = {x["normalized"]: Fraction(x["initial_moles"]) for x in species}
    coef = {x["normalized"]: x["coefficient"] for x in species}
    xi = Fraction(Decimal(d["reaction_extent"]["value"]))
    volume = Fraction(d["volume_l"])
    kc = Fraction(d["kc"])

    lo_bound = Fraction(Decimal(d["extent_bounds"]["minimum"]))
    hi_bound = Fraction(Decimal(d["extent_bounds"]["maximum"]))
    assert lo_bound <= xi <= hi_bound

    if qc_bracket:
        q_lo = _exact_qc(d, xi - tol, n0, coef, reactants, products, volume)
        q_hi = _exact_qc(d, xi + tol, n0, coef, reactants, products, volume)
        assert q_lo <= kc <= q_hi or q_hi <= kc <= q_lo

    for x in species:
        s = x["normalized"]
        sign = -1 if s in reactants else 1
        expected = n0[s] + sign * coef[s] * xi
        reported = Fraction(Decimal(x["equilibrium_moles"]))
        assert reported >= 0
        assert abs(reported - expected) <= tol


# ------------------------------------------------------------ 成功路径
def test_dissociation_forward_analytic_root():
    # PCl5 ⇌ PCl3 + Cl2，V=1，纯 PCl5 1 mol（产物给极小正初值 ε=1e-12），Kc=0.8
    # 含 ε 的解析根 x = -(2ε+0.8)/2 + sqrt(((2ε+0.8)/2)^2 + 0.8(1-ε^2))
    # ε→0 时趋近 0.4(√33 - 1) ≈ 0.8970
    resp = client.post(
        PATH,
        json=_body(
            ["PCl5"],
            ["PCl3", "Cl2"],
            {"PCl5": "1", "PCl3": "1e-12", "Cl2": "1e-12"},
            kc="0.8",
        ),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "forward"
    xi = Decimal(d["reaction_extent"]["value"])
    assert xi > 0
    with localcontext() as ctx:
        ctx.prec = 40
        # (x+ε)^2 = 0.8(1-x)（ε=1e-12）：
        # x = -(2ε+0.8)/2 + sqrt(((2ε+0.8)/2)^2 + 0.8(1-ε^2))
        eps = Decimal("1e-12")
        half = (2 * eps + Decimal("0.8")) / 2
        expected = -half + (half * half + Decimal("0.8") * (1 - eps * eps)).sqrt()
    assert abs(xi - expected) < Decimal("1e-9")
    assert d["extent_bounds"] == {"minimum": "-0.000000000001", "maximum": "1"}
    assert_solution_consistent(d)
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    assert blocks["PCl5"]["coefficient"] == 1
    assert blocks["PCl3"]["coefficient"] == blocks["Cl2"]["coefficient"] == 1
    assert Decimal(blocks["PCl5"]["equilibrium_moles"]) > 0


def test_exact_rational_root_forward():
    # 2O3 ⇌ 3O2：初始 (2, 3)，V=1，Kc=729/8=91.125 -> ξ = 1/2 为精确有理根
    # （O3=1、O2=4.5，Qc=4.5^3/1^2=91.125）
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="91.125"),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "forward"
    assert d["reaction_extent"]["value"] == "0.5"
    assert d["extent_bounds"] == {"minimum": "-1", "maximum": "1"}
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    assert blocks["O3"]["equilibrium_moles"] == "1"
    assert blocks["O2"]["equilibrium_moles"] == "4.5"
    assert d["initial_qc"] == "6.75"
    assert d["equilibrium_qc"] == "91.125"
    assert_solution_consistent(d)


def test_initial_qc_greater_than_kc_negative_extent():
    # 2O3 ⇌ 3O2：初始 (2, 3)，Kc=3/8=0.375 < Qc0=27/4=6.75
    # -> 逆向，精确根 ξ = -1/2（O3=3、O2=1.5，Qc=1.5^3/3^2=0.375）
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="0.375"),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "reverse"
    assert d["reaction_extent"]["value"] == "-0.5"
    assert d["extent_bounds"] == {"minimum": "-1", "maximum": "1"}
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    assert blocks["O3"]["equilibrium_moles"] == "3"
    assert blocks["O2"]["equilibrium_moles"] == "1.5"
    assert d["equilibrium_qc"] == "0.375"
    assert_solution_consistent(d)


def test_initial_state_already_at_equilibrium():
    # 2NO2 ⇌ N2O4，V=1，初始 (2, 1) -> Qc = 1/4 = Kc
    resp = client.post(
        PATH,
        json=_body(["NO2"], ["N2O4"], {"NO2": "2", "N2O4": "1"}, kc="0.25"),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "at_equilibrium"
    assert d["reaction_extent"] == {"value": "0", "unit": "mol"}
    assert d["initial_qc"] == d["equilibrium_qc"] == "0.25"
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    assert blocks["NO2"]["equilibrium_moles"] == "2"
    assert blocks["N2O4"]["equilibrium_moles"] == "1"


def test_volume_enters_concentration_quotient():
    # 2NO2 ⇌ N2O4（Δn = -1），V=2，初始 n_NO2=2、n_N2O4=1，Kc=1：
    # c_NO2 = 1-ξ，c_N2O4 = (1+ξ)/2；[(1+ξ)/2] / (1-ξ)^2 = 1
    # -> 2ξ^2-5ξ+1=0 -> ξ = (5-√17)/4 ≈ 0.219224
    resp = client.post(
        PATH,
        json=_body(
            ["NO2"], ["N2O4"], {"NO2": "2", "N2O4": "1"},
            volume="2", kc="1",
        ),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "forward"
    xi = Decimal(d["reaction_extent"]["value"])
    with localcontext() as ctx:
        ctx.prec = 40
        expected = (Decimal(5) - Decimal(17).sqrt()) / 4
    assert abs(xi - expected) < Decimal("1e-9")
    assert_solution_consistent(d)


def test_very_large_kc_drives_near_completion():
    # N2 + 3H2 ⇌ 2NH3（选用未在其它文件出现的初值），Kc=1e20，V=2
    resp = client.post(
        PATH,
        json=_body(
            ["N2", "H2"], ["NH3"],
            {"N2": "1", "H2": "3", "NH3": "1e-9"},
            volume="2", kc="1e20",
        ),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "forward"
    xi = Decimal(d["reaction_extent"]["value"])
    # 进度应非常接近上限 1（H2 限量）但严格小于
    assert Decimal("0.9999") < xi < 1
    # 近完全反应 Qc 极陡峭，用更紧的邻域核验内部精度
    assert_solution_consistent(d, tol=Fraction(1, 10**12))


def test_very_small_kc_drives_reverse_near_completion():
    resp = client.post(
        PATH,
        json=_body(
            ["N2", "H2"], ["NH3"],
            {"N2": "0.001", "H2": "0.001", "NH3": "1"},
            volume="0.5", kc="1e-20",
        ),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "reverse"
    xi = Decimal(d["reaction_extent"]["value"])
    # ξ_min = -0.5（NH3 系数 2），逆向应非常接近 -0.5
    assert Decimal("-0.5") < xi < Decimal("-0.4999999")
    # 近完全反应 Qc 极陡峭，用更紧的邻域核验内部精度（绝对误差仍远 < 1e-9）
    assert_solution_consistent(
        d, tol=Fraction(1, 10**12), qc_bracket=True
    )


# ------------------------------------------------------------ 极端平衡常数
def _reference_root(d, iterations=600):
    """独立精确二分 Qc(ξ)=Kc 的参考根（与求解器终止准则无关）。"""
    species = d["reactants"] + d["products"]
    reactants = [x["normalized"] for x in d["reactants"]]
    products = [x["normalized"] for x in d["products"]]
    n0 = {x["normalized"]: Fraction(x["initial_moles"]) for x in species}
    coef = {x["normalized"]: x["coefficient"] for x in species}
    volume = Fraction(d["volume_l"])
    kc = Fraction(d["kc"])
    lo = Fraction(Decimal(d["extent_bounds"]["minimum"]))
    hi = Fraction(Decimal(d["extent_bounds"]["maximum"]))
    if _exact_qc(d, Fraction(0), n0, coef, reactants, products, volume) > kc:
        hi = Fraction(0)
    else:
        lo = Fraction(0)
    for _ in range(iterations):
        mid = (lo + hi) / 2
        if _exact_qc(d, mid, n0, coef, reactants, products, volume) > kc:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def assert_extreme_kc_consistent(d):
    """极端 Kc 核验：内部 Qc、按响应进度/摩尔数复算的 Qc 相对 Kc 误差
    均 ≤ 1e-6，且响应进度距独立参考根 < 1e-9 mol、摩尔数非负。"""
    species = d["reactants"] + d["products"]
    reactants = [x["normalized"] for x in d["reactants"]]
    products = [x["normalized"] for x in d["products"]]
    n0 = {x["normalized"]: Fraction(x["initial_moles"]) for x in species}
    coef = {x["normalized"]: x["coefficient"] for x in species}
    volume = Fraction(d["volume_l"])
    kc = Fraction(d["kc"])
    rel_tol = Fraction(1, 10**6)

    # 内部平衡反应商相对 Kc 的误差 ≤ 1e-6
    eq_qc = Fraction(Decimal(d["equilibrium_qc"]))
    assert abs(eq_qc / kc - 1) <= rel_tol

    # 响应进度落在物料区间且距独立参考根 < 1e-9 mol
    xi = Fraction(Decimal(d["reaction_extent"]["value"]))
    lo_bound = Fraction(Decimal(d["extent_bounds"]["minimum"]))
    hi_bound = Fraction(Decimal(d["extent_bounds"]["maximum"]))
    assert lo_bound <= xi <= hi_bound
    assert abs(xi - _reference_root(d)) < Fraction(1, 10**9)

    # 按响应进度复算的反应商相对 Kc 的误差 ≤ 1e-6
    q_from_extent = _exact_qc(d, xi, n0, coef, reactants, products, volume)
    assert abs(q_from_extent / kc - 1) <= rel_tol

    # 按响应平衡摩尔数复算的反应商相对 Kc 的误差 ≤ 1e-6，摩尔数非负
    eq_moles = {
        x["normalized"]: Fraction(Decimal(x["equilibrium_moles"]))
        for x in species
    }
    num = Fraction(1)
    for s in products:
        num *= eq_moles[s] ** coef[s]
    den = Fraction(1)
    for s in reactants:
        den *= eq_moles[s] ** coef[s]
    delta_n = sum(coef[s] for s in products) - sum(coef[s] for s in reactants)
    q_from_moles = num / (den * volume ** delta_n)
    assert abs(q_from_moles / kc - 1) <= rel_tol
    for value in eq_moles.values():
        assert value >= 0


def test_extreme_tiny_kc_reverse_near_endpoint():
    # 2O3 ⇌ 3O2：初始 (2, 3)，V=1，Kc=1e-100 -> 逆向贴近 ξ_min = -1
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="1e-100"),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "reverse"
    assert_extreme_kc_consistent(d)
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    # 近端点的微量剩余必须保留：O2 剩余 ~1e-33 mol，不得被抹成 0
    residual = Decimal(blocks["O2"]["equilibrium_moles"])
    assert Decimal(0) < residual < Decimal("1e-30")


def test_extreme_huge_kc_forward_near_endpoint():
    # 2O3 ⇌ 3O2：初始 (2, 3)，V=1，Kc=1e100 -> 正向贴近 ξ_max = 1
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="1e100"),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["direction"] == "forward"
    assert_extreme_kc_consistent(d)
    blocks = {x["normalized"]: x for x in d["reactants"] + d["products"]}
    # 近端点的微量剩余必须保留：O3 剩余 ~1e-49 mol，不得被抹成 0
    residual = Decimal(blocks["O3"]["equilibrium_moles"])
    assert Decimal(0) < residual < Decimal("1e-45")


def test_ionic_reaction_charge_balanced():
    # CH3COOH ⇌ CH3COO- + H+
    resp = client.post(
        PATH,
        json=_body(
            ["CH3COOH"], ["CH3COO-", "H+"],
            {"CH3COOH": "1", "CH3COO-": "1e-9", "H+": "1e-9"},
            kc="1.8e-5",
        ),
    )
    assert resp.status_code == 200, resp.text
    d = resp.json()
    assert d["normalized"]["equation"] == "CH3COOH -> CH3COO^- + H^+"
    assert d["direction"] == "forward"
    assert_solution_consistent(d)


def test_equivalent_formula_writings_match():
    # Ca(OH)2 + 2HCl ⇌ CaCl2 + 2H2O：括号风格与内部空白的等价写法共享签名
    def payload(reaction_writings):
        reactants_w, products_w = reaction_writings
        return {
            "reactants": reactants_w,
            "products": products_w,
            "amounts": [
                {"species": "Ca(OH)2", "moles": "1"},
                {"species": "HCl", "moles": "3"},
                {"species": "CaCl2", "moles": "0.001"},
                {"species": "H2O", "moles": "1"},
            ],
            "volume_l": "2",
            "kc": "7",
        }

    first = client.post(PATH, json=payload((["Ca(OH)2", "HCl"], ["CaCl2", "H2O"])))
    second = client.post(
        PATH,
        json=payload((["Ca [OH]2", " HCl"], ["CaCl2", "H2 O"])),
    )
    assert first.status_code == second.status_code == 200
    a, b = first.json(), second.json()
    assert a["reaction_id"] == b["reaction_id"]
    assert a["reaction_extent"]["value"] == b["reaction_extent"]["value"]
    assert a["normalized"]["equation"] == b["normalized"]["equation"]


def test_input_echo_preserves_request_writing():
    # 初始量条目的原始写法（去首尾空白后）通过 input 字段回显
    resp = client.post(
        PATH,
        json={
            "reactants": ["Ca(OH)2", "HCl"],
            "products": ["CaCl2", "H2O"],
            "amounts": [
                {"species": "Ca [OH]2", "moles": "1"},
                {"species": "HCl", "moles": "3"},
                {"species": "CaCl2", "moles": "0.001"},
                {"species": "H2O", "moles": "1"},
            ],
            "volume_l": "2",
            "kc": "7",
        },
    )
    assert resp.status_code == 200, resp.text
    inputs = {x["input"] for x in resp.json()["reactants"]}
    assert "Ca [OH]2" in inputs


def test_reaction_id_shared_with_balance():
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["PCl5"], "products": ["PCl3", "Cl2"]},
    ).json()
    solved = client.post(
        PATH,
        json=_body(
            ["PCl5"], ["PCl3", "Cl2"],
            {"PCl5": "1", "PCl3": "1e-12", "Cl2": "1e-12"},
            kc="0.8",
        ),
    ).json()
    assert solved["reaction_id"] == balance["reaction_id"]


def test_decimal_strings_used_throughout():
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="91.125"),
    )
    d = resp.json()
    assert isinstance(d["reaction_extent"]["value"], str)
    assert isinstance(d["initial_qc"], str)
    assert isinstance(d["equilibrium_qc"], str)
    for x in d["reactants"] + d["products"]:
        assert isinstance(x["initial_moles"], str)
        assert isinstance(x["equilibrium_moles"], str)


# ------------------------------------------------------------ 错误路径
def test_zero_initial_moles_rejected():
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "0", "O2": "3"}, kc="1"),
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "INVALID_REQUEST"
    assert body["error"]["details"]["invalid_moles"][0]["moles"] == "0"


def test_negative_initial_moles_rejected():
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "-1", "O2": "3"}, kc="1"),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_non_positive_volume_and_kc_rejected():
    body = _body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, volume="0", kc="1")
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["invalid_volume"] == "0"

    body = _body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, volume="2", kc="-3")
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["invalid_kc"] == "-3"

    body = _body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, volume="x", kc="1")
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert "invalid_volume" in resp.json()["error"]["details"]


def test_json_number_rejected():
    body = _body(["O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="1")
    body["volume_l"] = 1.5
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_amounts_missing_duplicate_unknown():
    # 缺失一个物种
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"], "products": ["NH3"],
            "amounts": [
                {"species": "N2", "moles": "1"},
                {"species": "H2", "moles": "1"},
            ],
            "volume_l": "1", "kc": "1",
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["missing"] == ["NH3"]

    # 重复初始量记录
    body = {
        "reactants": ["N2", "H2"], "products": ["NH3"],
        "amounts": [
            {"species": "N2", "moles": "1"},
            {"species": "H2", "moles": "1"},
            {"species": "NH3", "moles": "1"},
            {"species": "H2", "moles": "2"},
        ],
        "volume_l": "1", "kc": "1",
    }
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["duplicate"] == ["H2"]

    # 未知物种（在重复记录基础上替换第二条）
    body["amounts"][1] = {"species": "He", "moles": "1"}
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["unknown"] == ["He"]


def test_same_side_duplicate_species_conflict():
    resp = client.post(
        PATH,
        json=_body(["O3", "O3"], ["O2"], {"O3": "2", "O2": "3"}, kc="1"),
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_cross_side_duplicate_species_rejected():
    # N2 同时出现在两侧：可逆反应不允许，配平前直接报 422 INVALID_REQUEST
    body = {
        "reactants": ["N2", "H2"], "products": ["N2", "NH3"],
        "amounts": [
            {"species": "N2", "moles": "1"},
            {"species": "H2", "moles": "1"},
            {"species": "NH3", "moles": "1"},
        ],
        "volume_l": "1", "kc": "1",
    }
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    result = resp.json()
    assert result["error"]["code"] == "INVALID_REQUEST"
    assert result["error"]["details"]["cross_side_duplicate"] == ["N2"]


def test_non_unique_balance_rejected():
    body = {
        "reactants": ["C", "O2"], "products": ["CO", "CO2"],
        "amounts": [
            {"species": s, "moles": "1"} for s in ("C", "O2", "CO", "CO2")
        ],
        "volume_l": "1", "kc": "1",
    }
    resp = client.post(PATH, json=body)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NON_UNIQUE_BALANCE"


def test_formula_parse_errors():
    # 反应式物种无法解析
    resp = client.post(
        PATH,
        json=_body(["Xx"], ["Cl2"], {"Xx": "1", "Cl2": "2"}, kc="1"),
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"

    # 初始量物种无法解析（定位到 amounts.species）
    resp = client.post(
        PATH,
        json=_body(["O3"], ["O2"], {"O3": "2", "Zz": "3"}, kc="1"),
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["field"] == "amounts.species"


def test_missing_fields_422():
    resp = client.post(PATH, json={"reactants": ["N2"]})
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"
