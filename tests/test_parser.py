"""化学式解析器测试。"""

import pytest

from app.parser import FormulaError, parse_formula


@pytest.mark.parametrize(
    "text,normalized,charge",
    [
        ("H2O", "H2O", 0),
        ("Ca(OH)2", "Ca(OH)2", 0),
        ("K4[Fe(CN)6]", "K4(Fe(CN)6)", 0),
        ("Mg3[Fe(CN)6]2", "Mg3(Fe(CN)6)2", 0),
        ("Al2(SO4)3", "Al2(SO4)3", 0),
        ("SO4^2-", "SO4^2-", -2),
        ("[Cu(NH3)4]^2+", "(Cu(NH3)4)^2+", 2),
        ("Fe3+", "Fe^3+", 3),           # 单原子：数字是电荷
        ("O2-", "O^2-", -2),
        ("Ba2+", "Ba^2+", 2),
        ("NH4+", "NH4^+", 1),           # 多原子：数字是下标
        ("MnO4-", "MnO4^-", -1),
        ("OH-", "OH^-", -1),
        ("H+", "H^+", 1),
        ("SO4 2-", "SO4^2-", -2),       # 电荷前空格
        ("O^{2-}", "O^2-", -2),
        ("Fe^{3+}", "Fe^3+", 3),
        ("e-", "e^-", -1),
        ("e", "e^-", -1),
        ("e+", "e^+", 1),
        ("  H2O  ", "H2O", 0),
        ("C2H5OH", "C2H5OH", 0),
        ("(NH4)2Cr2O7", "(NH4)2Cr2O7", 0),
    ],
)
def test_parse_valid(text, normalized, charge):
    parsed = parse_formula(text)
    assert parsed.normalized == normalized
    assert parsed.charge == charge


def test_composition():
    parsed = parse_formula("Mg3[Fe(CN)6]2")
    assert parsed.composition == {"Mg": 3, "Fe": 2, "C": 12, "N": 12}


def test_electron_key():
    assert parse_formula("e-").is_electron is True
    assert parse_formula("e").key == (((), -1))


@pytest.mark.parametrize(
    "text",
    [
        "",
        "   ",
        "H2O)",
        "(H2O",
        "Xx",
        "h2o",
        "2H",
        "H 2 O!",
        "()",
        "H0",
        "H^0+",
        "(e-)2",
        "H2^3",
        "Fe^",
        "e2-",
        "[Fe",
        "Fe++",
        "100",
        "H2O 3",
        "^-",
        "H^{}",
    ],
)
def test_parse_invalid(text):
    with pytest.raises(FormulaError):
        parse_formula(text)
