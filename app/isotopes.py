"""同位素名义质量峰分布（isotopic envelope）的精确枚举。

给定化学式（每种元素的原子数由既有解析器
:func:`app.parser.parse_formula` 给出）与每种元素的同位素清单
（质量数 + 非负丰度，丰度之和恰为 1），把分子的每个原子看作一次
**独立**的同位素抽样：

* 对含 n 个原子、k 种同位素（质量数 m_1..m_k、丰度 p_1..p_k）的元素，
  名义质量 M 与对应概率通过多项分布枚举全部计数向量
  ``(c_1, …, c_k)``（c_i ≥ 0、Σc_i = n）得到：

      P(c_1,…,c_k) = n! / (c_1!…c_k!) · p_1^{c_1}…p_k^{c_k}
      M(c_1,…,c_k) = Σ c_i·m_i

  同一元素内不同计数向量给出相同名义质量时先做一次合并（如 Cl 的
  35/37 各两原子时，``35+37`` 与 ``37+35`` 计入同一质量峰）；
* 再把各元素的（质量 -> 概率）分布做卷积，跨元素同质量组合同样合并。

枚举**不设概率截断**：所有可达的名义质量（含恰好为 0 的概率峰，由
零丰度同位素贡献）都会保留，合并在全部组合枚举完成之后进行。

约束（由业务层在调用前校验）：元素种类 ≤ 3、总原子数 ≤ 30、每种元素
至多 3 种同位素，因此组合规模至多 3^30 量级但经多项分布合并后实际
状态数为 C(n+k-1, k-1)（k ≤ 3、n ≤ 30 时单元素至多 465 个计数向量，
跨元素卷积后的峰数至多数百），全程使用
:class:`fractions.Fraction` 精确运算，不经过 float。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial
from typing import Sequence

#: 业务层允许的元素种类上限
MAX_ELEMENTS = 3
#: 总原子数上限
MAX_TOTAL_ATOMS = 30
#: 每种元素的同位素种数上限
MAX_ISOTOPES_PER_ELEMENT = 3


@dataclass(frozen=True)
class Isotope:
    """一种同位素：正整数质量数与非负丰度（精确有理数）。"""

    mass_number: int
    abundance: Fraction


@dataclass(frozen=True)
class ElementIsotopes:
    """一种元素的同位素清单及其在化学式中的原子数。"""

    element: str
    atom_count: int
    isotopes: tuple[Isotope, ...]


@dataclass(frozen=True)
class MassPeak:
    """一个名义质量峰：名义质量（正整数）与精确概率（最简有理数）。"""

    mass: int
    probability: Fraction


def element_mass_distribution(spec: ElementIsotopes) -> list[MassPeak]:
    """按多项分布枚举一种元素全部原子的（名义质量, 概率）分布。

    枚举该元素所有满足 Σc_i = n 的同位素计数向量（包括零丰度同位素以
    正计数出现、概率恰为 0 的向量——不做任何低概率截断），同质量向量在
    枚举完成后合并；结果按质量升序。
    """
    n = spec.atom_count
    isotopes = spec.isotopes
    masses = tuple(iso.mass_number for iso in isotopes)
    probabilities = tuple(iso.abundance for iso in isotopes)
    multinomial_factor = factorial(n)

    # 记忆化递归枚举计数向量；suffix(i, remaining) = 该后缀（同位素 i..k-1
    # 分配 remaining 个原子）的 {质量增量: 概率}，概率含组合系数但不含
    # 前面已选部分。零概率向量同样保留，仅在最后按质量合并。
    @lru_cache(maxsize=None)
    def suffix(i: int, remaining: int) -> dict[int, Fraction]:
        if i == len(isotopes) - 1:
            # 最后一种同位素吃掉全部剩余计数，只有一种向量
            c = remaining
            weight = Fraction(multinomial_factor, factorial(c))
            return {c * masses[i]: weight * (probabilities[i] ** c)}

        merged: dict[int, Fraction] = {}
        for c in range(remaining + 1):
            base_prob = Fraction(1, factorial(c)) * (probabilities[i] ** c)
            for tail_mass, tail_prob in suffix(i + 1, remaining - c).items():
                mass = c * masses[i] + tail_mass
                merged[mass] = merged.get(mass, Fraction(0)) + base_prob * tail_prob
        return merged

    distribution = suffix(0, n) if n > 0 else {0: Fraction(1)}
    return [MassPeak(mass, prob) for mass, prob in sorted(distribution.items())]


def _convolve(
    left: Sequence[MassPeak], right: Sequence[MassPeak]
) -> list[MassPeak]:
    """两个元素分布的卷积：质量相加、概率相乘，同质量组合合并。"""
    merged: dict[int, Fraction] = {}
    for a in left:
        for b in right:
            mass = a.mass + b.mass
            merged[mass] = merged.get(mass, Fraction(0)) + a.probability * b.probability
    return [MassPeak(mass, merged[mass]) for mass in sorted(merged)]


def isotope_envelope(specs: Sequence[ElementIsotopes]) -> list[MassPeak]:
    """计算整个化学式的名义质量峰分布。

    逐元素求多项分布后做卷积，跨元素同质量组合合并；返回按质量升序的
    全部峰，不截断低/零概率组合。空清单（理论上不会出现，化学式至少含
    一种元素）返回单个质量 0、概率 1 的峰。
    """
    peaks: list[MassPeak] = [MassPeak(0, Fraction(1))]
    for spec in specs:
        peaks = _convolve(peaks, element_mass_distribution(spec))
    return peaks


def most_probable_peak(peaks: Sequence[MassPeak]) -> MassPeak:
    """返回概率最高的峰；概率并列时取质量较小者。

    调用方须保证 ``peaks`` 按质量升序，因此直接以概率为键取最大值即可
    （``max`` 在并列时保留先出现者，即质量较小者）。
    """
    if not peaks:
        raise ValueError("峰列表不能为空")
    return max(peaks, key=lambda peak: peak.probability)
