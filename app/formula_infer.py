"""元素分析质量百分比反推分子式（C/H/N/O）的精确枚举。

给定 C、H、N、O 四种元素的实测质量百分比与绝对容差（百分数，闭区间
0～100；各元素测量值不要求合计恰为 100），以及名义摩尔质量的正整数
闭区间 ``[mass_min, mass_max]``，枚举全部满足

    n_C, n_H, n_N, n_O ∈ [0, 30]，n_C + n_H + n_N + n_O > 0，
    mass_min ≤ 12·n_C + n_H + 14·n_N + 16·n_O ≤ mass_max，
    |w_e(n) − p_e| ≤ t_e，  w_e = 100·a_e·n_e / M，  e ∈ {C,H,N,O}

的整数计数向量（原子量取 C=12、H=1、N=14、O=16；百分比容差与质量
区间均为**闭区间**，边界值计入匹配）。输入十进制先无损转为
:class:`fractions.Fraction`，热循环只做整数乘除与比较（容差边界
``L_e·M ≤ 100·a_e·n_e ≤ U_e·M`` 全部交叉相乘为整数式），不经过
``float``，因此恰落在容差边界上的候选不会因浮点误差被误判。

枚举规模有界：31⁴ − 1 = 923,520 个非零计数向量。实现上先按质量区间
与各元素容差带把每种元素的全局计数范围收紧到闭区间，再以
``(n_C, n_N, n_O)`` 三重循环枚举；固定三者后质量 P = 12n_C + 14n_N +
16n_O 为常数，H 计数与四个容差带都表现为对总质量 M = P + n_H 的闭
区间约束，取交集后区间内的每个整数 M 即一个可行候选，无任何重复
枚举。

匹配结果排序键依次为：总百分比绝对偏差
``Σ_e |w_e − p_e|``（精确分数）、摩尔质量 M、n_C、n_H、n_N、n_O，
全部升序；仅保留前 :data:`RESULT_LIMIT` 个（以大小顶堆在线维护，
不保留全部候选对象）。
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass
from fractions import Fraction
from math import lcm
from typing import Mapping

#: 参与反推的四种元素（同时决定响应/回显的固定顺序）
ELEMENTS: tuple[str, ...] = ("C", "H", "N", "O")

#: 名义原子量：C=12、H=1、N=14、O=16
ATOMIC_MASS: dict[str, int] = {"C": 12, "H": 1, "N": 14, "O": 16}

#: 每种元素原子数上限（含端点）
MAX_COUNT = 30

#: 可达最大名义摩尔质量 = 30·(12+1+14+16)
MAX_MOLAR_MASS = MAX_COUNT * sum(ATOMIC_MASS.values())

#: 返回结果条数上限
RESULT_LIMIT = 20


@dataclass(frozen=True)
class ElementMeasurement:
    """一种元素的实测质量百分比与绝对容差（均为 0～100 的精确百分数）。"""

    element: str
    percent: Fraction
    tolerance: Fraction


@dataclass(frozen=True)
class FormulaCandidate:
    """一个可行分子式：C/H/N/O 计数、摩尔质量、各项计算百分比与总偏差。"""

    c: int
    h: int
    n: int
    o: int
    molar_mass: int
    mass_percent: dict[str, Fraction]
    total_abs_deviation: Fraction


# ------------------------------------------------------------ 整数区间工具
def _floor_div(numerator: int, denominator: int) -> int:
    """floor(numerator / denominator)，分母为正、分子允许为负。"""
    return numerator // denominator


def _ceil_div(numerator: int, denominator: int) -> int:
    """ceil(numerator / denominator)，分母为正、分子允许为负。"""
    return -((-numerator) // denominator)


def render_formula(c: int, h: int, n: int, o: int) -> str:
    """按 C/H/N/O 顺序渲染分子式；计数为 0 的元素省略、为 1 时省略下标。"""
    counts = {"C": c, "H": h, "N": n, "O": o}
    return "".join(
        element if count == 1 else element + str(count)
        for element, count in counts.items()
        if count > 0
    )


# ------------------------------------------------------------ 全局计数范围
def _global_count_range(
    element: str,
    mass_min: int,
    mass_max: int,
    lower: Fraction,
    upper: Fraction,
) -> tuple[int, int]:
    """由质量闭区间与容差带求该元素计数的必要全局范围（只剪枝、不漏解）。

    可行计数 n 必须同时满足：

    * 质量可实现性：``a·n ≤ mass_max``（其余元素质量非负）且
      ``a·n + 其余元素最大质量 ≥ mass_min``；
    * 容差带在某个 M ∈ [mass_min, mass_max] 内可达：
      ``lower·M ≤ 100·a·n ≤ upper·M``，下界随 M 增大故用 mass_min
      得必要下界，上界随 M 增大故用 mass_max 得必要上界。
    """
    mass = ATOMIC_MASS[element]
    rest_max = MAX_MOLAR_MASS - MAX_COUNT * mass

    low = max(0, _ceil_div(max(0, mass_min - rest_max), mass))
    high = min(MAX_COUNT, mass_max // mass)

    if lower > 0:
        low = max(
            low,
            _ceil_div(
                lower.numerator * mass_min,
                lower.denominator * 100 * mass,
            ),
        )
    if upper == 0:
        # 上界为 0 时该元素质量贡献必须为 0，即计数为 0
        high = min(high, 0)
    else:
        high = min(
            high,
            _floor_div(
                upper.numerator * mass_max,
                upper.denominator * 100 * mass,
            ),
        )
    return low, high


class _ReverseOrder:
    """反转比较顺序的包装：以最小堆模拟大小顶堆（堆顶为已保留的最大键）。

    排序键第一项是总偏差的未约分有理数 ``dev_num/(molar_mass·common_den)``
    （同一请求内 ``common_den`` 为常数），比较时交叉相乘、全程整数运算；
    其后依次为摩尔质量与 C/H/N/O 计数。
    """

    __slots__ = ("dev_num", "molar_mass", "counts", "common_den")

    def __init__(self, dev_num: int, molar_mass: int, counts: tuple[int, ...],
                 common_den: int):
        self.dev_num = dev_num
        self.molar_mass = molar_mass
        self.counts = counts
        self.common_den = common_den

    def _rank(self, other: "_ReverseOrder") -> int:
        """正常（未反转）排序键比较：-1 小于 / 0 相等 / 1 大于。"""
        left = self.dev_num * other.molar_mass
        right = other.dev_num * self.molar_mass
        if left != right:
            return -1 if left < right else 1
        if self.molar_mass != other.molar_mass:
            return -1 if self.molar_mass < other.molar_mass else 1
        if self.counts != other.counts:
            return -1 if self.counts < other.counts else 1
        return 0

    def __lt__(self, other: "_ReverseOrder") -> bool:
        # 反转：排序键更大的元素在堆中“更小”，使堆顶恒为已保留键的最大值
        return self._rank(other) > 0


# ------------------------------------------------------------ 主枚举
def infer_formulas(
    measurements: Mapping[str, ElementMeasurement],
    mass_min: int,
    mass_max: int,
) -> tuple[list[FormulaCandidate], int]:
    """枚举全部可行分子式，返回 ``(按排序键升序的前 RESULT_LIMIT 个, 总数)``。

    ``measurements`` 必须恰好提供 :data:`ELEMENTS` 中四种元素的测量值
    （百分比/容差均为 0～100 的百分数）；``mass_min ≤ mass_max`` 均为
    正整数。无可行解时返回 ``([], 0)``。
    """
    bands_frac: dict[str, tuple[Fraction, Fraction]] = {
        element: (
            measurements[element].percent - measurements[element].tolerance,
            measurements[element].percent + measurements[element].tolerance,
        )
        for element in ELEMENTS
    }
    # 容差带预取为整数系数（None 表示该方向不约束），热循环不再触碰 Fraction
    bands: dict[str, tuple[tuple[int, int] | None, tuple[int, int] | None]] = {}
    for element in ELEMENTS:
        lower, upper = bands_frac[element]
        # None 表示该方向不约束：下界 L>0 才有 M ≤ q/L；上界 U=0 退化为
        # 必须 q=0（在枚举内单独判定），否则给出 M ≥ q/U
        lower_coef = (
            (lower.numerator, lower.denominator) if lower > 0 else None
        )
        upper_coef = (
            (upper.numerator, upper.denominator) if upper > 0 else None
        )
        bands[element] = (lower_coef, upper_coef)

    ranges: dict[str, tuple[int, int]] = {
        element: _global_count_range(
            element, mass_min, mass_max, bands_frac[element][0],
            bands_frac[element][1],
        )
        for element in ELEMENTS
    }
    if any(low > high for low, high in ranges.values()):
        return [], 0

    # 总偏差 Σ|100·a·n/M − p| = Σ|100·a·n·pd − pn·M|·(D/pd) / (M·D)：
    # 热循环只算整数分子 dev_num，分数比较交叉相乘；D 为各 p 分母的最小公倍数
    measured = {element: measurements[element].percent for element in ELEMENTS}
    common_den = lcm(*(measured[element].denominator for element in ELEMENTS))
    deviation_parts = tuple(
        (
            ATOMIC_MASS[element],
            measured[element].numerator,
            measured[element].denominator,
            common_den // measured[element].denominator,
        )
        for element in ELEMENTS
    )

    c_low, c_high = ranges["C"]
    n_low, n_high = ranges["N"]
    o_low, o_high = ranges["O"]

    # H 容差带改写为对 M 的整数系数（q_H = 100(M−P)）：
    # (100−U_H)·M ≤ 100·P -> M ≤ floor(100·P·u_den/u_diff)
    # (100−L_H)·M ≥ 100·P -> M ≥ ceil (100·P·l_den/l_diff)
    c_band_lower, c_band_upper = bands["C"]
    n_band_lower, n_band_upper = bands["N"]
    o_band_lower, o_band_upper = bands["O"]
    h_lower, h_upper = bands_frac["H"]
    if h_upper >= 100:
        h_upper_bound: tuple[int, int] | None = None  # U≥100：上界不约束
    else:
        h_upper_bound = (
            h_upper.denominator,
            100 * h_upper.denominator - h_upper.numerator,
        )
    if h_lower >= 100:
        # w_H = 100 要求 P = 0（全部质量来自 H）
        h_lower_full = True
        h_lower_bound: tuple[int, int] | None = None
    else:
        h_lower_full = False
        h_lower_bound = (
            h_lower.denominator,
            100 * h_lower.denominator - h_lower.numerator,
        )

    heap: list[_ReverseOrder] = []
    total = 0

    def consider(dev_num: int, molar_mass: int, counts: tuple[int, ...]) -> None:
        """在线保留全局最小的 RESULT_LIMIT 个排序键。"""
        nonlocal total
        total += 1
        wrapped = _ReverseOrder(dev_num, molar_mass, counts, common_den)
        if len(heap) < RESULT_LIMIT:
            heapq.heappush(heap, wrapped)
        elif heap[0] < wrapped:
            # 堆顶是已保留键中的最大值；新键比它更小（正常顺序）时替换
            heapq.heapreplace(heap, wrapped)

    for c in range(c_low, c_high + 1):
        for n in range(n_low, n_high + 1):
            for o in range(o_low, o_high + 1):
                partial_mass = 12 * c + 14 * n + 16 * o

                # n_H ∈ [0,30] 与质量闭区间给出的总质量 M 基础窗口
                low_m = max(mass_min, partial_mass)
                high_m = min(mass_max, partial_mass + MAX_COUNT)
                if low_m > high_m:
                    continue

                feasible = True
                # 已固定的 C/N/O：容差带 L·M ≤ 100·a·count ≤ U·M 收 M 区间
                fixed = (
                    ("C", c, c_band_lower, c_band_upper),
                    ("N", n, n_band_lower, n_band_upper),
                    ("O", o, o_band_lower, o_band_upper),
                )
                for _element, count, lower_coef, upper_coef in fixed:
                    q = 100 * ATOMIC_MASS[_element] * count
                    if upper_coef is None:
                        if q != 0:
                            feasible = False
                            break
                    else:
                        unum, uden = upper_coef
                        # M ≥ q / U（U=0 时无下界）
                        bound = _ceil_div(q * uden, unum)
                        if bound > low_m:
                            low_m = bound
                    if lower_coef is not None:
                        lnum, lden = lower_coef
                        # M ≤ q / L
                        bound = _floor_div(q * lden, lnum)
                        if bound < high_m:
                            high_m = bound
                    if low_m > high_m:
                        feasible = False
                        break
                if not feasible:
                    continue

                # H 的百分比随 M 变化，容差带改写为对 M 的闭区间
                if h_upper_bound is not None:
                    u_den, u_diff = h_upper_bound
                    bound = _floor_div(100 * partial_mass * u_den, u_diff)
                    if bound < high_m:
                        high_m = bound
                if h_lower_full:
                    if partial_mass != 0:
                        continue
                else:
                    l_den, l_diff = h_lower_bound
                    bound = _ceil_div(100 * partial_mass * l_den, l_diff)
                    if bound > low_m:
                        low_m = bound
                if low_m > high_m:
                    continue

                for molar_mass in range(low_m, high_m + 1):
                    h = molar_mass - partial_mass
                    counts = (c, h, n, o)
                    dev_num = 0
                    for index, (mass, pn, pd, scale) in enumerate(deviation_parts):
                        term = 100 * mass * counts[index] * pd - pn * molar_mass
                        dev_num += abs(term) * scale
                    consider(dev_num, molar_mass, counts)

    kept = sorted(
        heap,
        key=lambda wrapped: (
            Fraction(wrapped.dev_num, wrapped.molar_mass * wrapped.common_den),
            wrapped.molar_mass,
            wrapped.counts,
        ),
    )
    candidates: list[FormulaCandidate] = []
    for wrapped in kept:
        c, h, n, o = wrapped.counts
        candidates.append(
            FormulaCandidate(
                c=c,
                h=h,
                n=n,
                o=o,
                molar_mass=wrapped.molar_mass,
                mass_percent={
                    element: Fraction(
                        100 * ATOMIC_MASS[element] * count, wrapped.molar_mass
                    )
                    for element, count in zip(ELEMENTS, wrapped.counts)
                },
                total_abs_deviation=Fraction(
                    wrapped.dev_num, wrapped.molar_mass * wrapped.common_den
                ),
            )
        )
    return candidates, total
