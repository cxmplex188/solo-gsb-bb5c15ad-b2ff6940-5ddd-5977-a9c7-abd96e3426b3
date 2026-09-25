"""化学反应配平核心。

方法：把各元素原子守恒与总电荷守恒写成齐次线性方程组 ``A x = 0``
（反应物列取正、生成物列取负），用 :class:`fractions.Fraction` 做精确
高斯-若尔当消元（RREF）求得零空间：

* 零空间维数为 1 且解向量所有分量同号非零时，把解化为最小互质正整数，
  得到唯一配平；
* 维数为 0 → 不可能守恒；维数 > 1 → 配平不唯一（自由度 ``d-1``）；
* 维数为 1 但存在 0 分量或符号不一致 → 该反应无法成立。

另外 :func:`enumerate_variants` 在零空间维数 1～3 时枚举所有
``1..max_coefficient`` 范围内的互质正整数解，供 ``/balance/variants`` 使用。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm
from typing import Iterable

from .parser import ParsedFormula, parse_formula


class BalanceError(Exception):
    """无法配平时抛出，携带机器可读的错误码。"""

    def __init__(self, code: str, message: str, **details: object) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


@dataclass(frozen=True)
class Species:
    """一侧中的一个化学物种。"""

    raw: str              # 用户输入（去首尾空白）
    normalized: str       # 规范化化学式
    formula: ParsedFormula


@dataclass(frozen=True)
class SolvedReaction:
    """配平结果。系数按 reactants + products 的顺序给出。"""

    reactants: tuple[Species, ...]
    products: tuple[Species, ...]
    coefficients: tuple[int, ...]
    elements: tuple[str, ...]

    @property
    def all_species(self) -> tuple[Species, ...]:
        return self.reactants + self.products

    @property
    def reactant_coefficients(self) -> tuple[int, ...]:
        return self.coefficients[: len(self.reactants)]

    @property
    def product_coefficients(self) -> tuple[int, ...]:
        return self.coefficients[len(self.reactants) :]


@dataclass(frozen=True)
class VariantSet:
    """零空间枚举结果：所有 1..max_coefficient 内的互质正整数配平。"""

    reactants: tuple[Species, ...]
    products: tuple[Species, ...]
    coefficients: tuple[tuple[int, ...], ...]
    elements: tuple[str, ...]
    nullity: int


# ---------------------------------------------------------------- 解析一侧
def parse_side(texts: Iterable[str], side_name: str) -> tuple[Species, ...]:
    """解析一侧的化学式，去重并按规范键排序。"""
    formulas: list[Species] = []
    seen_keys: set[tuple] = set()
    for raw in texts:
        formula = parse_formula(raw)
        if formula.key in seen_keys:
            raise BalanceError(
                "DUPLICATE_SPECIES",
                f"{side_name}中存在重复物种：{formula.normalized}",
                species=formula.normalized,
                side="reactants" if side_name == "反应物" else "products",
            )
        seen_keys.add(formula.key)
        formulas.append(Species(raw=raw.strip(), normalized=formula.normalized, formula=formula))
    formulas.sort(key=lambda sp: sp.formula.key)
    return tuple(formulas)


# ---------------------------------------------------------------- 配平入口
def solve_reaction(
    reactant_texts: Iterable[str],
    product_texts: Iterable[str],
) -> SolvedReaction:
    reactants = parse_side(reactant_texts, "反应物")
    products = parse_side(product_texts, "生成物")

    if not reactants:
        raise BalanceError("INVALID_REQUEST", "反应物列表不能为空")
    if not products:
        raise BalanceError("INVALID_REQUEST", "生成物列表不能为空")
    if len(reactants) + len(products) > 100:
        raise BalanceError("INVALID_REQUEST", "物种总数不能超过 100")

    species = reactants + products
    n_cols = len(species)

    element_set = sorted(
        {el for sp in species for el in sp.formula.composition}
    )
    matrix = _build_conservation_matrix(reactants, products, element_set)

    basis = _null_space(matrix)

    if len(basis) == 0:
        raise BalanceError(
            "IMPOSSIBLE_BALANCE",
            "元素或电荷守恒方程无解：该反应在给定反应物与生成物下无法配平",
            elements=element_set,
        )
    if len(basis) > 1:
        degrees = len(basis) - 1
        raise BalanceError(
            "NON_UNIQUE_BALANCE",
            f"配平不唯一：零空间维数为 {len(basis)}，存在 {degrees} 个自由度。"
            "通常意味着缺少物种（如电子 e-、H2O、H+ / OH-）或两侧物种不独立。",
            nullity=len(basis),
            degrees_of_freedom=degrees,
            elements=element_set,
        )

    vec = basis[0]
    # 合法解要求所有物种非零且整体同号（矩阵列已带方向，见上）
    if any(x == 0 for x in vec):
        zero_species = [
            species[i].normalized for i, x in enumerate(vec) if x == 0
        ]
        raise BalanceError(
            "ZERO_COEFFICIENT",
            "守恒解要求部分物种的系数为 0，该反应无法作为一个整体配平；"
            "请把无关物种从反应式中去掉。",
            zero_coefficient_species=zero_species,
            elements=element_set,
        )
    signs = {1 if x > 0 else -1 for x in vec}
    if len(signs) != 1:
        raise BalanceError(
            "INVALID_STOICHIOMETRY",
            "守恒解中各物种的系数符号不一致：部分物种必须移到另一侧，"
            "或该反应在化学上不成立。",
            elements=element_set,
        )
    if -1 in signs:
        vec = [-x for x in vec]

    ints = _to_positive_integers(vec)

    # 最后做一次整数复核，防止任何意外
    for element in element_set:
        left = sum(
            ints[i] * species[i].formula.composition.get(element, 0)
            for i in range(len(reactants))
        )
        right = sum(
            ints[i] * species[i].formula.composition.get(element, 0)
            for i in range(len(reactants), n_cols)
        )
        if left != right:
            raise BalanceError(  # 理论上不可达
                "IMPOSSIBLE_BALANCE",
                f"内部校验失败：元素 {element} 不守恒",
                elements=element_set,
            )
    charge_left = sum(
        ints[i] * species[i].formula.charge for i in range(len(reactants))
    )
    charge_right = sum(
        ints[i] * species[i].formula.charge
        for i in range(len(reactants), n_cols)
    )
    if charge_left != charge_right:
        raise BalanceError(
            "CHARGE_NOT_CONSERVED",
            f"总电荷不守恒：反应物侧 {charge_left:+d}，生成物侧 {charge_right:+d}。"
            "氧化还原半反应请记得加入电子 e-。",
            reactant_charge=charge_left,
            product_charge=charge_right,
        )

    return SolvedReaction(
        reactants=reactants,
        products=products,
        coefficients=tuple(ints),
        elements=tuple(element_set),
    )


# ------------------------------------------------------------ 多解枚举
MAX_VARIANT_FREE_VARIABLES = 3
VARIANT_SPACE_TOO_LARGE = "VARIANT_SPACE_TOO_LARGE"


def enumerate_variants(
    reactant_texts: Iterable[str],
    product_texts: Iterable[str],
    max_coefficient: int,
) -> VariantSet:
    """枚举零空间中所有分量落在 ``1..max_coefficient`` 的互质正整数解。

    零空间维数（自由变量数）为 0 时抛出 ``IMPOSSIBLE_BALANCE``；超过
    :data:`MAX_VARIANT_FREE_VARIABLES` 时抛出 ``VARIANT_SPACE_TOO_LARGE``。
    维数 1～3 时，在 ``1..max_coefficient`` 内对自由变量做精确有理数枚举
    （自由分量本身也必须是正系数），过滤掉非整数、超界、最大公约数大于 1
    的向量，结果按系数总和与向量字典序排序。不存在正整数解时返回空集合。
    """
    reactants = parse_side(reactant_texts, "反应物")
    products = parse_side(product_texts, "生成物")
    return enumerate_variants_parsed(reactants, products, max_coefficient)


def enumerate_variants_parsed(
    reactants: tuple[Species, ...],
    products: tuple[Species, ...],
    max_coefficient: int,
) -> VariantSet:
    """与 :func:`enumerate_variants` 相同，但入参为已解析、去重排序的物种。"""
    if not reactants:
        raise BalanceError("INVALID_REQUEST", "反应物列表不能为空")
    if not products:
        raise BalanceError("INVALID_REQUEST", "生成物列表不能为空")
    if len(reactants) + len(products) > 100:
        raise BalanceError("INVALID_REQUEST", "物种总数不能超过 100")

    species = reactants + products
    element_set = sorted(
        {el for sp in species for el in sp.formula.composition}
    )
    matrix = _build_conservation_matrix(reactants, products, element_set)
    basis = _null_space(matrix)
    nullity = len(basis)

    if nullity == 0:
        raise BalanceError(
            "IMPOSSIBLE_BALANCE",
            "元素或电荷守恒方程无解：该反应在给定反应物与生成物下无法配平",
            elements=element_set,
        )
    if nullity > MAX_VARIANT_FREE_VARIABLES:
        raise BalanceError(
            VARIANT_SPACE_TOO_LARGE,
            f"零空间维数为 {nullity}（{nullity - 1} 个自由度），"
            f"超过变体枚举上限 {MAX_VARIANT_FREE_VARIABLES}，"
            "请补充缺失物种（如电子 e-、H2O、H+ / OH-）后重试。",
            nullity=nullity,
            max_free_variables=MAX_VARIANT_FREE_VARIABLES,
            elements=element_set,
        )

    vectors = _enumerate_positive_primitive(basis, nullity, max_coefficient)
    vectors.sort(key=lambda vec: (sum(vec), vec))
    return VariantSet(
        reactants=reactants,
        products=products,
        coefficients=tuple(tuple(vec) for vec in vectors),
        elements=tuple(element_set),
        nullity=nullity,
    )


def _enumerate_positive_primitive(
    basis: list[list[Fraction]], nullity: int, max_coefficient: int
) -> list[list[int]]:
    """对 1～3 个自由变量枚举范围内的互质正整数解。

    RREF 基的第 k 个向量在第 k 个自由列上为 1、其余自由列为 0，因此令自由
    分量直接取 ``t_k in 1..max_coefficient``（正解的自由分量本身就是正系数），
    主元分量按 ``x_i = Σ b_k[i]·t_k`` 精确计算：全部为 1..max_coefficient 的
    整数且整体 gcd 为 1 才接受。每个本原正整数解恰有一个本原表示，gcd 过滤
    同时完成去重。自由值至多 ``30³`` 组。

    每个主元行的各基系数预先通分（分母 LCM），枚举热循环里只有整数乘加，
    并用“分子是否被分母整除”判定整数，逐行提前剪枝。
    """
    n = len(basis[0])
    # RREF 构造中第 k 个基向量在第 k 个自由列上为 1、其余自由列为 0
    free_positions = [
        next(
            i
            for i in range(n)
            if basis[k][i] == 1
            and all(j == k or basis[j][i] == 0 for j in range(nullity))
        )
        for k in range(nullity)
    ]
    pivot_positions = [i for i in range(n) if i not in free_positions]

    # (主元列, 分母 LCM, 通分后的整数权重)
    pivot_rows: list[tuple[int, int, tuple[int, ...]]] = []
    for i in pivot_positions:
        denom = lcm(*(basis[k][i].denominator for k in range(nullity)))
        weights = tuple(int(basis[k][i] * denom) for k in range(nullity))
        pivot_rows.append((i, denom, weights))

    results: list[list[int]] = []
    values = range(1, max_coefficient + 1)

    def accept(t: tuple[int, ...]) -> None:
        vec = [0] * n
        for pos, tk in zip(free_positions, t):
            vec[pos] = tk
        for i, denom, weights in pivot_rows:
            numerator = sum(w * tk for w, tk in zip(weights, t))
            if numerator % denom:
                return
            value = numerator // denom
            if value <= 0 or value > max_coefficient:
                return
            vec[i] = value
        common = 0
        for v in vec:
            common = gcd(common, v)
        if common != 1:
            return
        results.append(vec)

    if nullity == 1:
        for t0 in values:
            accept((t0,))
    elif nullity == 2:
        for t0 in values:
            for t1 in values:
                accept((t0, t1))
    else:
        for t0 in values:
            for t1 in values:
                for t2 in values:
                    accept((t0, t1, t2))
    return results


# ---------------------------------------------------------------- 线性代数
def _build_conservation_matrix(
    reactants: tuple[Species, ...],
    products: tuple[Species, ...],
    element_set: list[str],
) -> list[list[Fraction]]:
    """组装元素守恒 + 总电荷守恒的齐次方程组 A x = 0。

    守恒式（左侧原子数/电荷 = 右侧）写成齐次式时反应物列取正、生成物列取负；
    电子没有任何元素，电荷行保证矩阵至少有一行。
    """
    species = reactants + products
    n_cols = len(species)
    matrix: list[list[Fraction]] = []
    for element in element_set:
        row = [Fraction(sp.formula.composition.get(element, 0)) for sp in species]
        for j in range(len(reactants), n_cols):
            row[j] = -row[j]
        matrix.append(row)
    charge_row = [Fraction(sp.formula.charge) for sp in species]
    for j in range(len(reactants), n_cols):
        charge_row[j] = -charge_row[j]
    matrix.append(charge_row)
    return matrix


def _rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    """对有理数矩阵做高斯-若尔当消元，返回 (RREF, 主元列)。"""
    m = [row[:] for row in matrix]
    if not m:
        return m, []
    rows, cols = len(m), len(m[0])
    pivot_cols: list[int] = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        pivot = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [x - factor * y for x, y in zip(m[i], m[r])]
        pivot_cols.append(c)
        r += 1
    return m, pivot_cols


def _null_space(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """返回矩阵零空间的一组基（每列对应一个基向量）。"""
    if not matrix:
        n = 0
    else:
        n = len(matrix[0])
    rref, pivots = _rref(matrix)
    free_cols = [c for c in range(n) if c not in pivots]
    basis: list[list[Fraction]] = []
    for free in free_cols:
        vec = [Fraction(0) for _ in range(n)]
        vec[free] = Fraction(1)
        for row_idx, col_idx in enumerate(pivots):
            vec[col_idx] = -rref[row_idx][free]
        basis.append(vec)
    return basis


def _to_positive_integers(vec: list[Fraction]) -> list[int]:
    """把有理向量按最小公倍数放大，再除以最大公约数，得到互质正整数。"""
    denominators = [x.denominator for x in vec]
    scale = lcm(*denominators) if denominators else 1
    scaled = [int(x * scale) for x in vec]
    common = 0
    for value in scaled:
        common = gcd(common, abs(value))
    if common == 0:
        raise BalanceError("IMPOSSIBLE_BALANCE", "得到全零系数，无法配平")
    return [v // common for v in scaled]


# ---------------------------------------------------------------- 守恒视图
def side_totals(
    species: Iterable[Species],
    coefficients: Iterable[int],
    elements: Iterable[str],
) -> dict[str, int]:
    """计算一侧每种元素的原子总数与总电荷（键 ``__charge__``）。"""
    species = list(species)
    coefficients = list(coefficients)
    totals = {el: 0 for el in elements}
    charge = 0
    for sp, coef in zip(species, coefficients):
        for el in elements:
            totals[el] += coef * sp.formula.composition.get(el, 0)
        charge += coef * sp.formula.charge
    totals["__charge__"] = charge
    return totals


def build_balanced_view(solved: SolvedReaction) -> dict:
    """配平后的逐元素守恒表。"""
    r_tot = side_totals(
        solved.reactants, solved.reactant_coefficients, solved.elements
    )
    p_tot = side_totals(
        solved.products, solved.product_coefficients, solved.elements
    )
    rows = [
        {
            "element": el,
            "reactant_total": r_tot[el],
            "product_total": p_tot[el],
            "difference": r_tot[el] - p_tot[el],
            "balanced": r_tot[el] == p_tot[el],
        }
        for el in solved.elements
    ]
    return {
        "elements": rows,
        "charge": {
            "reactant_total": r_tot["__charge__"],
            "product_total": p_tot["__charge__"],
            "difference": r_tot["__charge__"] - p_tot["__charge__"],
            "balanced": r_tot["__charge__"] == p_tot["__charge__"],
        },
        "all_balanced": all(row["balanced"] for row in rows)
        and r_tot["__charge__"] == p_tot["__charge__"],
    }


def build_raw_view(
    reactants: Iterable[Species],
    products: Iterable[Species],
    elements: Iterable[str] | None = None,
) -> dict:
    """按单位系数（用户原始写法，不做配平）计算逐元素差异。"""
    reactants = tuple(reactants)
    products = tuple(products)
    if elements is None:
        elements = sorted(
            {el for sp in (*reactants, *products) for el in sp.formula.composition}
        )
    elements = tuple(elements)
    ones_r = [1] * len(reactants)
    ones_p = [1] * len(products)
    r_tot = side_totals(reactants, ones_r, elements)
    p_tot = side_totals(products, ones_p, elements)
    return {
        "coefficients_assumed": 1,
        "elements": [
            {
                "element": el,
                "reactant_total": r_tot[el],
                "product_total": p_tot[el],
                "difference": r_tot[el] - p_tot[el],
                "balanced": r_tot[el] == p_tot[el],
            }
            for el in elements
        ],
        "charge": {
            "reactant_total": r_tot["__charge__"],
            "product_total": p_tot["__charge__"],
            "difference": r_tot["__charge__"] - p_tot["__charge__"],
            "balanced": r_tot["__charge__"] == p_tot["__charge__"],
        },
    }
