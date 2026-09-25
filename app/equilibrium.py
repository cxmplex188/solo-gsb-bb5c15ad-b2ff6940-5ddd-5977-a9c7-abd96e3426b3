"""单个可逆反应的理想溶液平衡组成求解。

对唯一配平后的反应（最小互质正整数系数）

    a1 R1 + a2 R2 + …  ⇌  b1 P1 + b2 P2 + …

给定各参与物种的**正**初始摩尔数 n0_i、溶液体积 V（L）与浓度平衡常数 Kc，
以反应进度 ξ（mol，按配平式书写方向为正）写物料守恒：

    n_R(ξ) = n0_R − a·ξ        （反应物，正向进度被消耗）
    n_P(ξ) = n0_P + b·ξ        （生成物，正向进度被产出）

非负物料区间由各物种初始量与配平系数确定：

    ξ_min = −min_P(n0_P / b) ≤ 0 ≤ ξ_max = min_R(n0_R / a)

因各初始摩尔数均严格为正，必有 ξ_min < 0 < ξ_max，区间为开区间，
端点处至少一个物种浓度为 0。理想溶液中以各物种浓度的计量次方构造

    Qc(ξ) = Π_P (n_P/V)^b  /  Π_R (n_R/V)^a

Qc 在 (ξ_min, ξ_max) 上关于 ξ 严格单调递增（端点处取 0 / +∞），因此
``Qc = Kc``（Kc > 0）在区间内存在**唯一内点**：

* Qc(0) < Kc：反应正向进行，根在 (0, ξ_max)；
* Qc(0) > Kc：反应逆向进行（允许负进度），根在 (ξ_min, 0)；
* Qc(0) = Kc：初始即平衡，ξ = 0。

求解全程采用 :class:`fractions.Fraction` 二分（符号判定逐位精确，绝无
浮点抵消或误判）。每次迭代取中点并精确计算 f = num − Kc·den（与
Qc − Kc 同号且严格递增），以下两个终止准则**同时**满足时返回该中点：

* 进度精度：中点成为新区间端点后区间宽度 < 5e-10 mol，故返回点距真根
  不超过 5e-10 mol（规范要求绝对误差 < 1e-9 mol）；
* 反应商精度：|Qc(ξ)/Kc − 1| = |f| / (Kc·den) ≤ 1e-9。贴近物料区间
  端点时 Qc 对 ξ 极其平缓（dlnQc/dξ ~ 计量系数 ÷ 端点残余量），任何
  固定的进度宽度准则都会在 Qc 仍偏离 Kc 数个数量级时提前停下（如
  Kc = 1e∓100 时偏差达数百倍乃至 1e±33）；相对误差准则迫使二分在
  非负物料区间内继续逼近，保证任意数量级 Kc 下内部反应商都贴合 Kc，
  且远严于规范要求的 1e-6，为输出渲染留足余量。

极端输入（如 Kc = 1e±100）下对分次数约 350 次、分数位数数百位，仍在
毫秒级；最终进度与各物种平衡摩尔数均为精确有理数，仅在 API 渲染层转
十进制。返回解同时给出解处对数斜率 S = dlnQc/dξ（精确有理数），供
渲染层按 1/S 量级自适应选择绝对精度，避免把贴近端点的微量剩余抹掉。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

from .balancer import Species
from .cache import species_key

# 进度宽度容差：返回点距真根不超过当前区间宽度，宽度 < 5e-10 mol 时
# 进度绝对误差 < 5e-10 mol（规范要求 < 1e-9 mol）
WIDTH_TOL = Fraction(5, 10**10)

# 反应商相对容差：|Qc(ξ)/Kc − 1| ≤ 1e-9 才允许停止，
# 远严于规范要求的 1e-6，为输出渲染的舍入留足余量
QC_REL_TOL = Fraction(1, 10**9)

# 对分次数保护上限：常规数量级约 60 次，Kc = 1e±100 约 350 次，
# 输入允许的最大数量级（1e±10000）约 3.5 万次
MAX_BISECTION_ITERS = 100_000


@dataclass(frozen=True)
class EquilibriumSolution:
    """平衡求解结果（进度与各物种平衡摩尔数均为精确 :class:`Fraction`）。"""

    extent: Fraction
    extent_min: Fraction
    extent_max: Fraction
    direction: str  # "forward" / "reverse" / "at_equilibrium"
    initial_qc: Fraction
    equilibrium_qc: Fraction
    equilibrium_moles: dict[str, Fraction]
    iterations: int
    # 解处对数斜率 S = dlnQc/dξ = Σ系数²/摩尔数（恒正），
    # 供渲染层按 1/S 量级选择绝对精度，保留贴近端点的微量剩余
    qc_log_slope: Fraction


def _qc_parts(
    xi: Fraction,
    r_keys: list[str],
    p_keys: list[str],
    coef: dict[str, int],
    n0: dict[str, Fraction],
    volume: Fraction,
) -> tuple[Fraction, Fraction]:
    """Qc = num / den：num = Π n_P^b，den = Π n_R^a · V^Δn。

    正向进度产出生成物（n_P = n0_P + bξ）、消耗反应物（n_R = n0_R − aξ）。
    """
    num = Fraction(1)
    for k in p_keys:
        num *= (n0[k] + coef[k] * xi) ** coef[k]
    den = Fraction(1)
    for k in r_keys:
        den *= (n0[k] - coef[k] * xi) ** coef[k]
    delta_n = sum(coef[k] for k in p_keys) - sum(coef[k] for k in r_keys)
    den *= volume ** delta_n
    return num, den


def solve_equilibrium(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    initial: dict[str, Fraction],
    volume: Fraction,
    kc: Fraction,
) -> EquilibriumSolution:
    """求 Qc(ξ) = Kc 的唯一内点进度与各物种平衡摩尔数。

    ``coefficients`` / ``initial`` 均以 :func:`app.cache.species_key` 为键；
    ``initial`` 必须覆盖反应两侧全部物种且值严格为正，``volume`` / ``kc``
    严格为正（合法性由业务层保证）。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)
    r_keys = [species_key(sp) for sp in reactants]
    p_keys = [species_key(sp) for sp in products]
    coef = coefficients
    n0 = initial

    # 非负物料区间：ξ_min = -min_生成物(n0/b)，ξ_max = min_反应物(n0/a)
    extent_min = -min(n0[k] / coef[k] for k in p_keys)
    extent_max = min(n0[k] / coef[k] for k in r_keys)

    num0, den0 = _qc_parts(Fraction(0), r_keys, p_keys, coef, n0, volume)
    initial_qc = num0 / den0

    if initial_qc == kc:
        return _build_solution(
            Fraction(0), extent_min=extent_min, extent_max=extent_max,
            direction="at_equilibrium", initial_qc=initial_qc,
            r_keys=r_keys, p_keys=p_keys, coef=coef, n0=n0, volume=volume,
            iterations=0,
        )

    if initial_qc > kc:
        lo, hi, direction = extent_min, Fraction(0), "reverse"
    else:
        lo, hi, direction = Fraction(0), extent_max, "forward"

    # ---- 精确有理数二分：f = num - Kc·den，与 Qc-Kc 同号且严格递增
    iterations = 0
    extent: Fraction | None = None
    for iterations in range(1, MAX_BISECTION_ITERS + 1):
        mid = (lo + hi) / 2
        num, den = _qc_parts(mid, r_keys, p_keys, coef, n0, volume)
        diff = num - kc * den
        if diff == 0:
            extent = mid  # 精确命中根
            break
        if diff > 0:
            hi = mid
        else:
            lo = mid
        # mid 已成为新区间端点：|mid − 真根| < hi − lo；
        # Qc(mid)/Kc − 1 = diff / (Kc·den)，两个终止准则同时满足才返回
        if (hi - lo) < WIDTH_TOL and abs(diff) <= QC_REL_TOL * kc * den:
            extent = mid
            break
    if extent is None:  # pragma: no cover - 保护上限远超任何合法输入所需次数
        extent = (lo + hi) / 2

    return _build_solution(
        extent, extent_min=extent_min, extent_max=extent_max,
        direction=direction, initial_qc=initial_qc,
        r_keys=r_keys, p_keys=p_keys, coef=coef, n0=n0, volume=volume,
        iterations=iterations,
    )


def _build_solution(
    extent: Fraction,
    *,
    extent_min: Fraction,
    extent_max: Fraction,
    direction: str,
    initial_qc: Fraction,
    r_keys: list[str],
    p_keys: list[str],
    coef: dict[str, int],
    n0: dict[str, Fraction],
    volume: Fraction,
    iterations: int,
) -> EquilibriumSolution:
    """按精确有理进度计算各物种平衡摩尔数与平衡 Qc 并组装结果。"""
    moles: dict[str, Fraction] = {}
    for k in r_keys:
        moles[k] = n0[k] - coef[k] * extent
    for k in p_keys:
        moles[k] = n0[k] + coef[k] * extent
    num = Fraction(1)
    for k in p_keys:
        num *= moles[k] ** coef[k]
    den = Fraction(1)
    for k in r_keys:
        den *= moles[k] ** coef[k]
    delta_n = sum(coef[k] for k in p_keys) - sum(coef[k] for k in r_keys)
    den *= volume ** delta_n
    # 解处对数斜率 S = dlnQc/dξ = Σ_生成物 b²/n_P + Σ_反应物 a²/n_R（恒正）
    log_slope = sum(coef[k] ** 2 / moles[k] for k in (*r_keys, *p_keys))
    return EquilibriumSolution(
        extent=extent,
        extent_min=extent_min,
        extent_max=extent_max,
        direction=direction,
        initial_qc=initial_qc,
        equilibrium_qc=num / den,
        equilibrium_moles=moles,
        iterations=iterations,
        qc_log_slope=log_slope,
    )
