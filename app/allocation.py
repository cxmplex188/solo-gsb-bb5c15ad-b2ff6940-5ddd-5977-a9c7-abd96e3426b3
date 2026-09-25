"""试剂批次领用分配（FEFO：先到期先出）。

给定待领用物种的**有效摩尔需求**与一组库存批次（批号、纯度、可用名义
摩尔数、失效日期、状态）以及核算日期，按以下规则分配：

* 批次按 **失效日期升序、同日按批号升序** 排序（FEFO）后依次领用；
* 每批的名义摩尔数乘纯度换算为**有效物质**，需求与扣减都以有效摩尔计；
  领用的名义摩尔数 = 有效领用 ÷ 纯度，保证不超过该批有效余量；
* **已过期**（失效日期早于核算日期）、**隔离**（``quarantined``）与
  **可用名义摩尔数为 0** 的批次跳过（前两类不参与扣减）；
* 需求可跨多个批次拆分；同一批号在请求中只允许出现一次（去重由业务层
  校验），因此同一批次在一次分配中至多被扣减一次。

为避免浮点误差，纯度、名义量、需求量在进入算法前由 :class:`Decimal`
无损转成 :class:`fractions.Fraction`，比较与扣减全部使用精确有理数，
仅在输出时经 :func:`fraction_to_decimal_string` 渲染为十进制字符串。
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from decimal import Decimal
from fractions import Fraction
from typing import Any

from .stoichiometry import (
    MAX_ABS_EXPONENT,
    decimal_to_fraction,
    fraction_to_decimal_string,
)

# 批次状态（与 Pydantic 模型中的字面量保持一致）
STATUS_AVAILABLE = "available"
STATUS_QUARANTINED = "quarantined"
BATCH_STATUSES = (STATUS_AVAILABLE, STATUS_QUARANTINED)

# 跳过原因：过期 / 隔离 / 名义量为零 / 批次物种不在领用清单中
REASON_EXPIRED = "expired"
REASON_QUARANTINED = "quarantined"
REASON_EXHAUSTED = "exhausted"
REASON_UNMATCHED = "unmatched_species"


@dataclass(frozen=True)
class Demand:
    """一条待领用物种记录（有效摩尔需求）。"""

    raw: str           # 客户端提交的原始化学式
    normalized: str    # 规范化化学式
    key: str           # 物种规范键（组成 + 电荷）
    moles: Fraction    # 有效摩尔需求（非负）


@dataclass(frozen=True)
class Batch:
    """一条库存批次记录。"""

    batch_no: str
    species_raw: str
    normalized: str
    key: str
    purity: Fraction           # 0 < 纯度 <= 1
    nominal_moles: Fraction    # 可用名义摩尔数（非负）
    expires_on: date
    status: str                # available / quarantined

    @property
    def effective_moles(self) -> Fraction:
        """该批可提供的有效物质 = 名义摩尔数 × 纯度。"""
        return self.nominal_moles * self.purity


def parse_purity(text: str) -> Decimal:
    """把客户端的纯度十进制字符串解析为 ``0 < p <= 1`` 的 :class:`Decimal`。

    拒绝非法输入、NaN / Infinity、零、负数、大于 1 与过大数量级；
    允许科学计数法（``Decimal`` 原生语法）。
    """
    value = Decimal(text.strip())
    if not value.is_finite():
        raise ValueError("纯度不能是 NaN 或无穷大")
    if value <= 0 or value > 1:
        raise ValueError("纯度必须在 (0, 1] 区间内（如 0.98、1.0）")
    if value.adjusted() > MAX_ABS_EXPONENT or value.as_tuple().exponent < -MAX_ABS_EXPONENT:
        raise ValueError(f"纯度的数量级超出允许范围（±1e{MAX_ABS_EXPONENT}）")
    return value


def _fefo_key(batch: Batch) -> tuple[date, str]:
    return batch.expires_on, batch.batch_no


def _batch_line(batch: Batch, *, as_of: date, remaining: Fraction) -> tuple[dict, Fraction, bool]:
    """处理 FEFO 序列中的一个批次，返回 (行明细, 处理后剩余需求, 是否实际扣减)。

    ``remaining`` 为进入该批时尚未满足的有效摩尔需求；过期、隔离、零名义量
    的批次原样保留余量并标注跳过原因，其余批次按 ``min(剩余需求, 有效余量)``
    扣减（同一批次只处理一次，不存在重复扣减）。
    """
    available = batch.effective_moles

    if batch.status == STATUS_QUARANTINED:
        skip_reason = REASON_QUARANTINED
    elif batch.expires_on < as_of:
        skip_reason = REASON_EXPIRED
    elif available == 0:
        skip_reason = REASON_EXHAUSTED
    else:
        skip_reason = None

    take_effective = Fraction(0)
    take_nominal = Fraction(0)
    used = False
    if skip_reason is None:
        take_effective = min(remaining, available)
        if take_effective > 0:
            used = True
            take_nominal = take_effective / batch.purity
            remaining -= take_effective

    line = {
        "batch_no": batch.batch_no,
        "species": batch.species_raw,
        "normalized": batch.normalized,
        "expires_on": batch.expires_on,
        "status": batch.status,
        "purity": batch.purity,
        "nominal_moles": batch.nominal_moles,
        "effective_moles": available,
        "eligible": skip_reason is None,
        "skipped": skip_reason is not None,
        "skip_reason": skip_reason,
        "used": used,
        "allocated_nominal_moles": take_nominal,
        "allocated_effective_moles": take_effective,
        "remaining_nominal_moles": batch.nominal_moles - take_nominal,
        "remaining_effective_moles": available - take_effective,
    }
    return line, remaining, used


def allocate(
    demands: list[Demand],
    batches: list[Batch],
    as_of: date,
) -> dict[str, Any]:
    """按 FEFO 规则把库存批次的有效物质分配给各待领用物种。

    返回原始（未经十进制渲染）的分配结果：物种按规范化化学式排序，物种下
    的批次按（失效日期, 批号）排序；物种之外的批次记入 ``unmatched_batches``
    （按 规范化物种、失效日期、批号 排序）。数值一律为 :class:`Fraction`，
    日期为 :class:`datetime.date`，由 :func:`render_allocation` 统一渲染。
    """
    batches_by_species: dict[str, list[Batch]] = {}
    for batch in batches:
        batches_by_species.setdefault(batch.key, []).append(batch)

    demand_keys = {demand.key for demand in demands}

    allocations: list[dict[str, Any]] = []
    total_requested = Fraction(0)
    total_allocated = Fraction(0)
    total_unmet = Fraction(0)
    fulfilled_demands = 0
    batches_used = 0
    batches_skipped = 0

    for demand in sorted(demands, key=lambda item: item.normalized):
        lines: list[dict[str, Any]] = []
        remaining = demand.moles
        allocated = Fraction(0)

        for batch in sorted(batches_by_species.get(demand.key, []), key=_fefo_key):
            line, remaining, used = _batch_line(batch, as_of=as_of, remaining=remaining)
            allocated += line["allocated_effective_moles"]
            if used:
                batches_used += 1
            if line["skipped"]:
                batches_skipped += 1
            lines.append(line)

        unmet = remaining
        fulfilled = unmet == 0
        if fulfilled:
            fulfilled_demands += 1

        total_requested += demand.moles
        total_allocated += allocated
        total_unmet += unmet

        allocations.append(
            {
                "species": demand.raw,
                "normalized": demand.normalized,
                "requested_effective_moles": demand.moles,
                "allocated_effective_moles": allocated,
                "unmet_effective_moles": unmet,
                "status": "fulfilled" if fulfilled else "shortage",
                "batches": lines,
            }
        )

    unmatched: list[dict[str, Any]] = []
    unmatched_batches = [b for b in batches if b.key not in demand_keys]
    for batch in sorted(
        unmatched_batches, key=lambda b: (b.normalized, b.expires_on, b.batch_no)
    ):
        line, _remaining, _used = _batch_line(
            batch, as_of=as_of, remaining=Fraction(0)
        )
        # 物种未被领用：跳过原因固定为 unmatched_species，整批保留
        line["eligible"] = False
        line["skipped"] = True
        line["skip_reason"] = REASON_UNMATCHED
        unmatched.append(line)

    return {
        "as_of": as_of,
        "allocations": allocations,
        "unmatched_batches": unmatched,
        "summary": {
            "total_requested_effective_moles": total_requested,
            "total_allocated_effective_moles": total_allocated,
            "total_unmet_effective_moles": total_unmet,
            "demands_total": len(demands),
            "fulfilled_demands": fulfilled_demands,
            "shortage_demands": len(demands) - fulfilled_demands,
            "batches_total": len(batches),
            "batches_used": batches_used,
            "batches_skipped": batches_skipped,
            "unmatched_batches": len(unmatched),
            "status": "fulfilled" if total_unmet == 0 else "shortage",
        },
    }


# 需要渲染成十进制字符串的数值字段
_LINE_DECIMAL_FIELDS = (
    "purity",
    "nominal_moles",
    "effective_moles",
    "allocated_nominal_moles",
    "allocated_effective_moles",
    "remaining_nominal_moles",
    "remaining_effective_moles",
)
_ALLOCATION_DECIMAL_FIELDS = (
    "requested_effective_moles",
    "allocated_effective_moles",
    "unmet_effective_moles",
)
_SUMMARY_DECIMAL_FIELDS = (
    "total_requested_effective_moles",
    "total_allocated_effective_moles",
    "total_unmet_effective_moles",
)


def _render_line(line: dict[str, Any]) -> dict[str, Any]:
    rendered = dict(line)
    for field_name in _LINE_DECIMAL_FIELDS:
        rendered[field_name] = fraction_to_decimal_string(rendered[field_name])
    rendered["expires_on"] = rendered["expires_on"].isoformat()
    return rendered


def render_allocation(result: dict[str, Any]) -> dict[str, Any]:
    """把 :func:`allocate` 的 Fraction/date 结果渲染为可 JSON 序列化的响应。"""
    allocations: list[dict[str, Any]] = []
    for block in result["allocations"]:
        rendered = dict(block)
        for field_name in _ALLOCATION_DECIMAL_FIELDS:
            rendered[field_name] = fraction_to_decimal_string(rendered[field_name])
        rendered["batches"] = [_render_line(line) for line in block["batches"]]
        allocations.append(rendered)

    summary = dict(result["summary"])
    for field_name in _SUMMARY_DECIMAL_FIELDS:
        summary[field_name] = fraction_to_decimal_string(summary[field_name])

    return {
        "as_of": result["as_of"].isoformat(),
        "allocations": allocations,
        "unmatched_batches": [_render_line(line) for line in result["unmatched_batches"]],
        "summary": summary,
    }
