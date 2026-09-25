"""配平核心测试。"""

from functools import reduce
from math import gcd

import pytest

from app.balancer import (
    BalanceError,
    VARIANT_SPACE_TOO_LARGE,
    build_balanced_view,
    build_raw_view,
    enumerate_variants,
    parse_side,
    solve_reaction,
)


def coefficient_map(solved) -> dict[str, int]:
    return {
        sp.normalized: c
        for sp, c in zip(
            solved.all_species, solved.coefficients
        )
    }


@pytest.mark.parametrize(
    "reactants,products,expected",
    [
        (["H2", "O2"], ["H2O"], {"H2": 2, "O2": 1, "H2O": 2}),
        (["Fe", "O2"], ["Fe2O3"], {"Fe": 4, "O2": 3, "Fe2O3": 2}),
        (["C3H8", "O2"], ["CO2", "H2O"],
         {"C3H8": 1, "O2": 5, "CO2": 3, "H2O": 4}),
        (["KMnO4", "HCl"], ["KCl", "MnCl2", "Cl2", "H2O"],
         {"KMnO4": 2, "HCl": 16, "KCl": 2, "MnCl2": 2, "Cl2": 5, "H2O": 8}),
        (["K4Fe(CN)6", "KMnO4", "H2SO4"],
         ["KHSO4", "Fe2(SO4)3", "MnSO4", "HNO3", "CO2", "H2O"],
         {"K4Fe(CN)6": 10, "KMnO4": 122, "H2SO4": 299, "KHSO4": 162,
          "Fe2(SO4)3": 5, "MnSO4": 122, "HNO3": 60, "CO2": 60, "H2O": 188}),
        (["Fe2+"], ["Fe3+", "e-"], {"Fe^2+": 1, "Fe^3+": 1, "e^-": 1}),
        (["MnO4-", "H+", "e-"], ["Mn2+", "H2O"],
         {"MnO4^-": 1, "H^+": 8, "e^-": 5, "Mn^2+": 1, "H2O": 4}),
        (["Cu", "NO3-", "H+"], ["Cu2+", "NO", "H2O"],
         {"Cu": 3, "NO3^-": 2, "H^+": 8, "Cu^2+": 3, "NO": 2, "H2O": 4}),
        (["H+", "OH-"], ["H2O"], {"H^+": 1, "OH^-": 1, "H2O": 1}),
        (["Ca(OH)2", "H3PO4"], ["Ca3(PO4)2", "H2O"],
         {"Ca(OH)2": 3, "H3PO4": 2, "Ca3(PO4)2": 1, "H2O": 6}),
        (["(NH4)2Cr2O7"], ["Cr2O3", "N2", "H2O"],
         {"(NH4)2Cr2O7": 1, "Cr2O3": 1, "N2": 1, "H2O": 4}),
        (["NH4Cl"], ["NH3", "HCl"], {"NH4Cl": 1, "NH3": 1, "HCl": 1}),
    ],
)
def test_balance(reactants, products, expected):
    solved = solve_reaction(reactants, products)
    assert coefficient_map(solved) == expected
    view = build_balanced_view(solved)
    assert view["all_balanced"] is True
    assert all(c > 0 for c in solved.coefficients)


def test_coefficients_are_coprime():
    solved = solve_reaction(["H2", "O2"], ["H2O"])
    assert reduce(gcd, solved.coefficients) == 1


def test_non_unique():
    with pytest.raises(BalanceError) as exc:
        solve_reaction(["C", "O2"], ["CO", "CO2"])
    assert exc.value.code == "NON_UNIQUE_BALANCE"
    assert exc.value.details["nullity"] == 2


def test_impossible():
    with pytest.raises(BalanceError) as exc:
        solve_reaction(["H2O"], ["Au"])
    assert exc.value.code == "IMPOSSIBLE_BALANCE"


def test_charge_not_conserved():
    with pytest.raises(BalanceError) as exc:
        solve_reaction(["Fe2+"], ["Fe3+"])
    assert exc.value.code == "IMPOSSIBLE_BALANCE"


def test_zero_coefficient():
    with pytest.raises(BalanceError) as exc:
        solve_reaction(["H2", "O2", "N2"], ["H2O"])
    assert exc.value.code == "ZERO_COEFFICIENT"
    assert "N2" in exc.value.details["zero_coefficient_species"]


def test_duplicate_species():
    with pytest.raises(BalanceError) as exc:
        solve_reaction(["H2", "H2"], ["H2O"])
    assert exc.value.code == "DUPLICATE_SPECIES"


def test_empty_side():
    with pytest.raises(BalanceError):
        solve_reaction([], ["H2O"])


def test_raw_view_difference():
    reactants = parse_side(["H2", "O2"], "反应物")
    products = parse_side(["H2O"], "生成物")
    view = build_raw_view(reactants, products)
    by_element = {row["element"]: row for row in view["elements"]}
    assert by_element["H"]["difference"] == 0
    assert by_element["O"]["difference"] == 1
    assert view["charge"]["difference"] == 0


# ------------------------------------------------------------ 变体枚举
def variant_vectors(reactants, products, max_coefficient):
    result = enumerate_variants(reactants, products, max_coefficient)
    return list(result.coefficients)


def test_variants_unique_reaction_returns_single_primitive():
    result = enumerate_variants(["H2", "O2"], ["H2O"], 30)
    assert result.nullity == 1
    assert result.coefficients == ((2, 1, 2),)


def test_variants_respects_max_coefficient_and_gcd():
    vectors = variant_vectors(["H2", "O2"], ["H2O"], 30)
    # (4,2,4)、(6,3,6)…… 都是 (2,1,2) 的倍数，必须被 gcd 过滤掉
    assert vectors == [(2, 1, 2)]


def test_variants_underdetermined_enumeration():
    # C + O2 -> CO + CO2，零空间维数 2
    result = enumerate_variants(["C", "O2"], ["CO", "CO2"], 10)
    assert result.nullity == 2
    vectors = [list(v) for v in result.coefficients]
    assert vectors[0] == [3, 2, 2, 1]
    assert len(vectors) == 15
    # 全部为 1..max_coefficient 的正整数且互质
    for vec in vectors:
        assert all(1 <= v <= 10 for v in vec)
        assert reduce(gcd, vec) == 1


def test_variants_sorting_by_sum_then_lexicographic():
    result = enumerate_variants(["C", "O2"], ["CO", "CO2"], 12)
    keys = [(sum(v), tuple(v)) for v in result.coefficients]
    assert keys == sorted(keys)


def test_variants_fractional_basis_still_yields_integer_solutions():
    # O2 -> O3：RREF 行带 1/2，本原解 (3,2) 必须能枚举到
    result = enumerate_variants(["O2"], ["O3"], 30)
    assert result.coefficients == ((3, 2),)
    # (6,4)、(9,6)…… 都是倍数解，只保留本原的一份
    assert len(result.coefficients) == 1


def test_variants_no_positive_solution_returns_empty():
    # 唯一解要求 N2 系数为 0
    result = enumerate_variants(["H2", "O2", "N2"], ["H2O"], 10)
    assert result.nullity == 1
    assert result.coefficients == ()


def test_variants_impossible():
    with pytest.raises(BalanceError) as exc:
        enumerate_variants(["H2O"], ["Au"], 10)
    assert exc.value.code == "IMPOSSIBLE_BALANCE"


def test_variants_nullity_3_allowed():
    result = enumerate_variants(["H2", "O2"], ["H2O", "H2O2", "O3"], 5)
    assert result.nullity == 3
    for vec in result.coefficients:
        assert all(1 <= v <= 5 for v in vec)
        assert reduce(gcd, vec) == 1


def test_variants_space_too_large():
    with pytest.raises(BalanceError) as exc:
        enumerate_variants(
            ["H2", "O2", "H2O", "H2O2", "O3"], ["HO2"], 10
        )
    assert exc.value.code == VARIANT_SPACE_TOO_LARGE
    assert exc.value.details["nullity"] == 4
    assert exc.value.details["max_free_variables"] == 3


def test_variants_charge_conserved_redox():
    result = enumerate_variants(
        ["MnO4-", "H+", "e-"], ["Mn2+", "H2O"], 30
    )
    # 物种按规范键排序：e^-, H^+, MnO4^-, H2O, Mn^2+
    assert result.coefficients == ((5, 8, 1, 4, 1),)
    view = build_balanced_view(
        solve_reaction(["MnO4-", "H+", "e-"], ["Mn2+", "H2O"])
    )
    assert view["all_balanced"] is True
