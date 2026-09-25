"""溶液配制（母液混合）的精确线性规划求解。

给定不超过 4 种母液（各自携带若干溶质的浓度与可用体积）、目标终体积与
各溶质目标浓度，求满足全部目标浓度的非负母液取用量：

* 每个目标溶质 j 的摩尔守恒：Σ_i c_ij · x_i = t_j · V
  （浓度 mol/L × 体积 mL，等式两侧同为毫摩尔量级，单位自洽）；
* 每种母液取用量不超过可用体积：0 ≤ x_i ≤ a_i；
* 母液总取用量不超过终体积：Σ_i x_i ≤ V，余量用水补足
  （w = V − Σx ≥ 0）；
* 多解时先取母液总量最小（Σx 最小），再按输入顺序取用量字典序最小。

可行域是闭有界多面体（盒约束 0 ≤ x_i ≤ a_i），非空时必有顶点，且
"(Σx, x_1, …, x_n) 字典序最小" 的最优解必为某一顶点（逐目标极小化
得到的面套序列最终收缩为单点，而面的面仍是原多面体的面，零维面即顶点）。
因此枚举全部顶点——等式组之外再激活 n − rank(等式) 个不等式，解
n×n 精确有理数方程组并校验全部约束——在可行顶点中取字典序最小者。
n ≤ 4 时候选不超过 C(2n+1, n) = C(9, 4) = 126 个，全部使用
:class:`fractions.Fraction` 精确运算，不经过 float。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


@dataclass(frozen=True)
class Stock:
    """一条母液：标识、各溶质浓度（mol/L）与可用体积（mL）。"""

    id: str
    concentrations: dict[str, Fraction]  # 溶质 -> 非负浓度（mol/L）
    available_volume: Fraction           # 可用体积（mL，非负）


@dataclass(frozen=True)
class BlendSolution:
    """配制解：各母液取用量（mL，与输入顺序一致）与补水量（mL）。"""

    volumes: tuple[Fraction, ...]
    water: Fraction


def _solve_square(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> list[Fraction] | None:
    """精确高斯-若尔当消元解 n×n 线性方程组；奇异（秩亏）时返回 None。"""
    n = len(matrix)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col] != 0), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        inv = 1 / aug[col][col]
        aug[col] = [value * inv for value in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [v - factor * w for v, w in zip(aug[r], aug[col])]
    return [aug[i][n] for i in range(n)]


def _reduce_equalities(
    rows: list[list[Fraction]], rhs: list[Fraction]
) -> tuple[list[list[Fraction]], list[Fraction]] | None:
    """把等式系统约化为线性无关行；出现矛盾（零行右端非零）时返回 None。"""
    if not rows:
        return [], []
    n = len(rows[0])
    aug = [list(row) + [rhs[i]] for i, row in enumerate(rows)]
    independent: list[list[Fraction]] = []
    start = 0
    for col in range(n):
        pivot = next((r for r in range(start, len(aug)) if aug[r][col] != 0), None)
        if pivot is None:
            continue
        aug[start], aug[pivot] = aug[pivot], aug[start]
        inv = 1 / aug[start][col]
        aug[start] = [value * inv for value in aug[start]]
        for r in range(len(aug)):
            if r != start and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [v - factor * w for v, w in zip(aug[r], aug[start])]
        independent.append(aug[start])
        start += 1
        if start == len(aug):
            break
    # 剩余行系数已全零：右端非零即矛盾（如目标溶质在任何母液中都不存在）
    for r in range(start, len(aug)):
        if aug[r][n] != 0:
            return None
    return [row[:n] for row in independent], [row[n] for row in independent]


def solve_blend(
    stocks: list[Stock],
    targets: dict[str, Fraction],
    final_volume: Fraction,
) -> BlendSolution | None:
    """求满足全部目标浓度的母液取用量；无解返回 None。

    ``targets`` 为 溶质 -> 非负目标浓度（mol/L）；``final_volume`` 为正
    （mL）。返回的解满足：全部目标溶质摩尔守恒、每种母液取用量不超过其
    可用体积、母液总量不超过终体积（余量即补水量），并在全部可行解中使
    （母液总量, 各母液取用量按输入顺序）字典序最小。
    """
    n = len(stocks)

    # ---- 等式约束：每个目标溶质一行 Σ_i c_ij·x_i = t_j·V（毫摩尔守恒）
    eq_rows = [
        [stock.concentrations.get(solute, Fraction(0)) for stock in stocks]
        for solute in targets
    ]
    eq_rhs = [target * final_volume for target in targets.values()]
    reduced = _reduce_equalities(eq_rows, eq_rhs)
    if reduced is None:
        return None
    eq_rows, eq_rhs = reduced

    # ---- 不等式约束（统一写成 row·x ≤ b）：可用体积上限、非负性、总量预算
    bounds: list[tuple[list[Fraction], Fraction]] = []
    for i, stock in enumerate(stocks):
        row = [Fraction(0)] * n
        row[i] = Fraction(1)
        bounds.append((row, stock.available_volume))  # x_i ≤ a_i
    for i in range(n):
        row = [Fraction(0)] * n
        row[i] = Fraction(-1)
        bounds.append((row, Fraction(0)))             # -x_i ≤ 0（x_i ≥ 0）
    bounds.append(([Fraction(1)] * n, final_volume))  # Σx ≤ V

    def feasible(point: list[Fraction]) -> bool:
        return all(
            sum(coef * value for coef, value in zip(row, point)) <= bound
            for row, bound in bounds
        )

    # ---- 顶点枚举：等式（秩 r）之外再激活 n − r 个不等式，解 n×n 方程组
    best_key: tuple[Fraction, ...] | None = None
    best_point: list[Fraction] | None = None
    need = n - len(eq_rows)
    for active in combinations(range(len(bounds)), need):
        matrix = [list(row) for row in eq_rows]
        rhs = list(eq_rhs)
        for index in active:
            row, bound = bounds[index]
            matrix.append(list(row))
            rhs.append(bound)
        point = _solve_square(matrix, rhs)
        if point is None or not feasible(point):
            continue
        key = (sum(point), *point)
        if best_key is None or key < best_key:
            best_key, best_point = key, point

    if best_point is None:
        return None
    return BlendSolution(
        volumes=tuple(best_point),
        water=final_volume - sum(best_point),
    )
