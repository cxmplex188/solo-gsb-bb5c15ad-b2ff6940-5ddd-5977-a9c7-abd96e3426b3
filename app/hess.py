"""赫斯定律（Hess's law）反应焓推导核心。

把每个已唯一配平的反应表示为规范物种上的**有符号计量向量**：生成物
分量取正、反应物分量取负（向量即 ``B·ν = 0`` 中带方向的列，反应物在
箭头左侧故取负）。设已知反应矩阵 ``A``（行为规范物种、列为已知反应）、
目标向量 ``t``，问题是求有符号系数 ``x``（允许反向与分数倍使用）使

    A x = t

全程以 :class:`fractions.Fraction` 精确运算：

* 对增广矩阵 ``[A | t]`` 做有理数高斯-若尔当消元（RREF）；
* 非主元行最后一列非零 -> 目标向量不在 ``A`` 的列空间内，不可推导；
* 否则取一个特解（自由变量置 0）与齐次解基，逐组检查所有已知反应焓的
  线性组合：若沿某齐次方向组合焓会变化 -> 结果歧义；否则反应焓被唯一
  确定（即使组合方式不唯一），正常返回。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .balancer import SolvedReaction, _rref
from .cache import species_key


class HessError(Exception):
    """无法推导目标反应焓时抛出，携带机器可读的错误码。"""

    def __init__(self, code: str, message: str, **details: object) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


@dataclass(frozen=True)
class KnownReaction:
    """一条已知反应：唯一配平结果与以 kJ 表示的反应焓（精确有理数）。"""

    solved: SolvedReaction
    enthalpy: Fraction

    @property
    def vector(self) -> dict[str, Fraction]:
        """规范物种键 -> 有符号计量系数（生成物正、反应物负）。"""
        return reaction_vector(self.solved)


def reaction_vector(solved: SolvedReaction) -> dict[str, Fraction]:
    """配平反应 -> 规范物种键上的有符号计量向量（生成物正、反应物负）。"""
    vector: dict[str, Fraction] = {}
    for sp, coef in zip(solved.all_species, solved.coefficients):
        sign = 1 if sp in solved.products else -1
        vector[species_key(sp)] = Fraction(sign * coef)
    return vector


@dataclass(frozen=True)
class HessSolution:
    """赫斯定律求解结果。

    * ``multipliers``：各已知反应的一个可用倍数（与输入顺序一致），
      自由变量取 0 的特解；
    * ``enthalpy``：推导出的目标反应焓（kJ，精确有理数）；
    * ``net_coefficients``：合并抵消后的逐物种净系数（生成物正、反应物负，
      含恰好抵消为 0 的中间物种），键为规范物种键；
    * ``species_order``：``net_coefficients`` 的输出顺序（规范键，按规范
      物种排序）。
    """

    multipliers: tuple[Fraction, ...]
    enthalpy: Fraction
    net_coefficients: dict[str, Fraction]
    species_order: tuple[str, ...]


def solve_hess(
    target: SolvedReaction,
    known: list[KnownReaction],
) -> HessSolution:
    """求已知反应的有符号线性组合，使其计量向量等于目标向量。

    抛出 :class:`HessError`：

    * ``HESS_UNDERIVABLE``：目标向量不在已知向量张成空间内；
    * ``HESS_AMBIGUOUS``：存在多组组合且导出不同反应焓。
    """
    known_keys = [key for kr in known for key in kr.vector]
    all_keys = sorted(set(known_keys) | set(reaction_vector(target)))
    n = len(known)

    target_vec = reaction_vector(target)
    rows: list[list[Fraction]] = []
    unmatched: list[str] = []
    for key in all_keys:
        if key in target_vec and not any(key in kr.vector for kr in known):
            unmatched.append(key)
        row = [kr.vector.get(key, Fraction(0)) for kr in known]
        row.append(target_vec.get(key, Fraction(0)))
        rows.append(row)

    if unmatched:
        raise HessError(
            "HESS_UNDERIVABLE",
            "目标反应无法由已知反应推导：目标中的物种未出现在任何已知反应中",
            target_only_species=unmatched,
        )

    rref, pivot_cols = _rref(rows)

    # 非主元行（恒等约束 0 = 非零值）说明右端不在列空间
    for row in rref:
        if all(row[c] == 0 for c in range(n)) and row[n] != 0:
            raise HessError(
                "HESS_UNDERIVABLE",
                "目标反应不在已知反应向量的张成空间内："
                "不存在任何已知反应的有符号（含反向、分数倍）线性组合"
                "能同时在全部物种上复现目标反应的计量关系",
            )

    # 特解（自由变量置 0）：主元变量 x_p = RREF 对应行最后一列
    free_cols = [c for c in range(n) if c not in pivot_cols]
    particular = [Fraction(0)] * n
    for row_idx, col in enumerate(pivot_cols):
        particular[col] = rref[row_idx][n]

    # 齐次解基：每个自由变量取 1、其余取 0，主元变量取 -RREF 系数
    null_basis: list[list[Fraction]] = []
    for free in free_cols:
        vec = [Fraction(0)] * n
        vec[free] = Fraction(1)
        for row_idx, col in enumerate(pivot_cols):
            vec[col] = -rref[row_idx][free]
        null_basis.append(vec)

    enthalpies = [kr.enthalpy for kr in known]
    particular_enthalpy = sum(
        (m * h for m, h in zip(particular, enthalpies)), Fraction(0)
    )
    # 沿任一齐次方向改变组合若会改变总焓，则反应焓不唯一 -> 歧义
    for vec in null_basis:
        delta = sum((m * h for m, h in zip(vec, enthalpies)), Fraction(0))
        if delta != 0:
            raise HessError(
                "HESS_AMBIGUOUS",
                "已知反应存在多组线性组合可复现目标反应的计量向量，"
                "但按各组已知反应焓导出的目标反应焓不一致，"
                "无法唯一确定（已知反应数据可能不自洽）",
            )

    multipliers = tuple(particular)
    net: dict[str, Fraction] = {}
    for key in all_keys:
        net[key] = sum(
            (m * kr.vector.get(key, Fraction(0)) for m, kr in zip(multipliers, known)),
            Fraction(0),
        )

    return HessSolution(
        multipliers=multipliers,
        enthalpy=particular_enthalpy,
        net_coefficients=net,
        species_order=tuple(all_keys),
    )


def normalized_by_key(solved: SolvedReaction) -> dict[str, str]:
    """配平反应中规范物种键 -> 规范化化学式。"""
    return {species_key(sp): sp.normalized for sp in solved.all_species}
