"""限量试剂与反应进度的精确计量计算。

提供两个方向的计量：

* 正向（投料 -> 进度）：:func:`compute_stoichiometry`
  给定唯一配平后的反应（最小互质正整数系数）与每种反应物的投料摩尔数
  （:class:`decimal.Decimal`），计算：
  * 反应进度 ξ = min(投料量 ÷ 配平系数)，系数为正整数，结果为非负小数；
  * 所有并列取得最小值的限量试剂（可多于一种，如投料恰好按计量比）；
  * 每种反应物的消耗量（ξ·系数）与剩余量（投料量 − 消耗量）；
  * 每种生成物的理论产量（ξ·系数）。

* 多步连续反应推演：:func:`apply_sequence_step`
  给定唯一配平后的一步反应与当前库存（物种规范键 -> 非负 :class:`Fraction`
  摩尔数），以该步全部反应物的当前库存（缺省视为 0）计算最大反应进度
  ξ = min(库存 ÷ 系数)，返回扣减消耗、加入产物后的新库存与逐步明细。
  库存不足以让任何反应发生（ξ = 0）时不改动库存，并标出全部阻塞反应物。
  库存始终以精确 :class:`fractions.Fraction` 在步骤间传递，仅在输出时
  经 :func:`fraction_to_decimal_string` 渲染为十进制字符串。

* 反向（目标 -> 规划投料）：:func:`compute_reverse_plan`
  给定一个或多个生成物的目标摩尔数，计算：
  * 满足全部目标所需的最小反应进度 ξ = max(目标量 ÷ 配平系数)；
  * 达成该进度所需的每种反应物理论需求量（ξ·系数）；
  * 每种生成物在该进度下的实际达成量；目标产物同时给出目标量与
    超额量（达成 − 目标，非目标产物超额量为 null）。

为避免浮点误差，输入十进制先无损转成 :class:`fractions.Fraction`，
比较与极值判定全部使用精确有理数；最终结果再按
:class:`decimal.Decimal`（28 位有效数字）格式化为字符串返回。
有限小数（分母仅含因子 2、5）给出精确值，无限小数按 28 位有效数字
银行家舍入。
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Iterable

from .balancer import Species
from .cache import species_key

# 无限小数输出时保留的有效数字位数
SIGNIFICANT_DIGITS = 28

# 允许的数量级范围（Decimal 的 10 进制指数），防止 1e99999… 撑爆有理数内存
MAX_ABS_EXPONENT = 10_000


def parse_moles(text: str) -> Decimal:
    """把客户端的十进制字符串解析为非负 :class:`Decimal`。

    拒绝非法输入、NaN / Infinity、负数与过大数量级；允许千位下划线
    与科学计数法（``Decimal`` 原生语法）。
    """
    value = Decimal(text.strip())
    if not value.is_finite():
        raise ValueError("摩尔数不能是 NaN 或无穷大")
    # is_signed() 同时覆盖负数与 "-0" / "-0.0" 这类负零写法
    if value.is_signed() and value <= 0:
        raise ValueError("摩尔数不能为负数")
    if value.adjusted() > MAX_ABS_EXPONENT or value.as_tuple().exponent < -MAX_ABS_EXPONENT:
        raise ValueError(f"摩尔数的数量级超出允许范围（±1e{MAX_ABS_EXPONENT}）")
    return value


def decimal_to_fraction(value: Decimal) -> Fraction:
    """Decimal -> Fraction 的无损精确转换（不经过 float）。"""
    sign, digits, exponent = value.as_tuple()
    numerator = int("".join(map(str, digits)))
    if sign:
        numerator = -numerator
    if exponent >= 0:
        return Fraction(numerator * 10**exponent, 1)
    return Fraction(numerator, 10**-exponent)


def parse_positive_moles(text: str) -> Decimal:
    """把十进制字符串解析为**严格为正**的有限 :class:`Decimal`。

    与 :func:`parse_moles` 的区别仅在于拒绝 0 与负零写法（平衡组成要求
    各参与物种初始摩尔数、体积与 Kc 严格为正，以保证 Qc=Kc 存在唯一内点）。
    """
    value = parse_moles(text)
    if value == 0:
        raise ValueError("数值必须严格为正，不能为 0")
    return value


def parse_enthalpy(text: str) -> Decimal:
    """把反应焓十进制字符串解析为**可带符号**的有限 :class:`Decimal`。

    与 :func:`parse_moles` 的区别仅在于允许负值（放热反应 ΔH<0）；
    同样拒绝 NaN / Infinity、非数字与过大数量级，允许科学计数法。
    """
    value = Decimal(text.strip())
    if not value.is_finite():
        raise ValueError("反应焓不能是 NaN 或无穷大")
    if value.adjusted() > MAX_ABS_EXPONENT or value.as_tuple().exponent < -MAX_ABS_EXPONENT:
        raise ValueError(f"反应焓的数量级超出允许范围（±1e{MAX_ABS_EXPONENT}）")
    return value


def _factor_count(n: int, p: int) -> int:
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def _terminating_places(denominator: int) -> int:
    """分母只含 2、5 时，精确十进制表示所需的小数位数 max(v2, v5)。"""
    twos = _factor_count(denominator, 2)
    rest = denominator // 2**twos
    fives = _factor_count(rest, 5)
    return max(twos, fives)


def _is_terminating(denominator: int) -> bool:
    twos = _factor_count(denominator, 2)
    rest = denominator // 2**twos
    rest //= 5 ** _factor_count(rest, 5)
    return rest == 1


def _fixed_point(decimal_value: Decimal) -> str:
    """把 Decimal 渲染为不带指数、不带尾随零的定点字符串。"""
    sign, digits, exponent = decimal_value.as_tuple()
    digit_text = "".join(map(str, digits))
    if exponent >= 0:
        body = digit_text + "0" * exponent
    else:
        places = -exponent
        if places < len(digit_text):
            body = digit_text[:-places] + "." + digit_text[-places:]
        else:
            body = "0." + "0" * (places - len(digit_text)) + digit_text
        if "." in body:
            body = body.rstrip("0").rstrip(".")
    return ("-" if sign and digit_text != "0" else "") + (body or "0")


def fraction_to_decimal_string(
    value: Fraction, *, min_absolute_accuracy: Decimal | None = None
) -> str:
    """把精确有理数量渲染为十进制字符串。

    分母约分后只含 2、5 因子时输出精确有限小数；否则按
    :data:`SIGNIFICANT_DIGITS` 位有效数字舍入（ROUND_HALF_EVEN）。
    始终使用定点记法（不出现科学计数法），并去掉多余的尾随零。

    给定 ``min_absolute_accuracy`` 时，对大数量级值追加有效位数，使
    渲染引入的绝对舍入误差严格小于该阈值（如平衡进度的 1e-10 mol）。
    """
    if value == 0:
        return "0"
    denominator = value.denominator
    terminating = _is_terminating(denominator)
    places = _terminating_places(denominator) if terminating else 0

    extra_prec = 0
    if min_absolute_accuracy is not None and not terminating:
        # 总有效位 = 整数位数 + 阈值要求的小数位数（1e-20 -> 20 位小数），
        # 使大数量级值的绝对舍入误差小于阈值
        integer_digits = max(len(str(value.numerator)) - len(str(denominator)) + 1, 1)
        required_places = max(0, -min_absolute_accuracy.as_tuple().exponent)
        extra_prec = max(0, integer_digits + required_places - SIGNIFICANT_DIGITS)

    with localcontext() as ctx:
        if terminating:
            # 精确除法：精度覆盖整数部分与所需小数位
            integer_digits = max(len(str(value.numerator)) - len(str(denominator)) + 1, 1)
            ctx.prec = integer_digits + places + 2
            decimal_value = Decimal(value.numerator) / Decimal(denominator)
            if places:
                decimal_value = decimal_value.quantize(Decimal(1).scaleb(-places))
        else:
            ctx.prec = SIGNIFICANT_DIGITS + extra_prec
            decimal_value = Decimal(value.numerator) / Decimal(denominator)
    return _fixed_point(decimal_value)


def decimal_to_string(
    value: Decimal, *, min_absolute_accuracy: Decimal | None = None
) -> str:
    """把高精度 :class:`Decimal` 渲染为定点十进制字符串。

    默认按 28 位有效数字舍入；给定 ``min_absolute_accuracy``（如
    平衡进度要求的 1e-10 mol）时，对大数量级值自动追加精度，使渲染
    引入的绝对舍入误差严格小于该阈值。用于以 Decimal 二分得到的平衡
    结果；输出规则（定点记法、去尾随零）与
    :func:`fraction_to_decimal_string` 保持一致。
    """
    with localcontext() as ctx:
        ctx.prec = SIGNIFICANT_DIGITS
        rounded = +value
        if min_absolute_accuracy is not None and value != 0:
            # 总有效位 = 整数位数 + 阈值小数位（大数量级时自动追加）
            integer_digits = value.adjusted() + 1
            required_places = max(
                0, -min_absolute_accuracy.as_tuple().exponent
            )
            required_prec = max(
                SIGNIFICANT_DIGITS, integer_digits + required_places
            )
            if required_prec > SIGNIFICANT_DIGITS:
                ctx.prec = required_prec
                rounded = +value
    return _fixed_point(rounded)


def _extent_and_limiting(
    reactants: list[Species],
    coefficients: dict[str, int],
    amounts_fractions: dict[str, Fraction],
) -> tuple[Fraction, set[str]]:
    """反应进度 ξ = min(投料量 ÷ 系数) 与全部并列限量试剂的规范键集合。"""
    ratios: dict[str, Fraction] = {
        species_key(sp): amounts_fractions[species_key(sp)]
        / coefficients[species_key(sp)]
        for sp in reactants
    }
    extent = min(ratios.values())
    limiting = {key for key, ratio in ratios.items() if ratio == extent}
    return extent, limiting


def compute_stoichiometry(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    feed: dict[str, Decimal],
) -> dict:
    """计算反应进度、限量试剂、各反应物消耗/剩余与各生成物理论产量。

    ``coefficients`` 与 ``feed`` 均以 :func:`app.cache.species_key` 为键；
    ``feed`` 必须覆盖全部反应物且值非负。

    物种按规范键排序输出，因此同一反应调换书写顺序得到的结果完全一致。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)

    feed_fractions = {key: decimal_to_fraction(amount) for key, amount in feed.items()}

    # 反应进度 ξ = min(投料量 / 系数)
    extent, limiting = _extent_and_limiting(reactants, coefficients, feed_fractions)

    reactant_blocks: list[dict] = []
    for sp in reactants:
        key = species_key(sp)
        coef = coefficients[key]
        fed = feed_fractions[key]
        consumed = extent * coef
        remaining = fed - consumed
        reactant_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "fed": fraction_to_decimal_string(fed),
                "consumed": fraction_to_decimal_string(consumed),
                "remaining": fraction_to_decimal_string(remaining),
                "is_limiting": key in limiting,
            }
        )

    product_blocks: list[dict] = []
    for sp in products:
        key = species_key(sp)
        coef = coefficients[key]
        produced = extent * coef
        product_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "theoretical_yield": fraction_to_decimal_string(produced),
            }
        )

    return {
        "reaction_extent": {
            "value": fraction_to_decimal_string(extent),
            "unit": "mol",
        },
        "limiting_reactants": [block["normalized"] for block in reactant_blocks
                               if block["is_limiting"]],
        "reactants": reactant_blocks,
        "products": product_blocks,
    }


def compute_reverse_plan(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    targets: Iterable[tuple[Species, Decimal]],
) -> dict:
    """反向投料规划：由生成物目标摩尔数反推最小反应进度与理论投料。

    ``coefficients`` 以 :func:`app.cache.species_key` 为键；``targets`` 为
    ``(目标产物物种, 非负目标摩尔数)`` 序列，目标物种必须都在生成物中，
    且规范键唯一（去重由业务层负责）。

    最小反应进度 ξ = max(目标量 ÷ 配平系数)——取最大值是为了让每个目标
    产物的达成量（ξ·系数）都不低于其目标量；全部并列取得最大值的目标
    产物恰好在达成量上无超额，其余目标产物存在正超额。

    物种按规范键排序输出，因此同一反应调换书写顺序得到的结果完全一致。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)

    target_fractions: dict[str, Fraction] = {
        species_key(sp): decimal_to_fraction(amount) for sp, amount in targets
    }
    target_raw: dict[str, str] = {
        species_key(sp): sp.raw for sp, _amount in targets
    }

    # 最小反应进度 ξ = max(目标量 / 系数)
    ratios: dict[str, Fraction] = {
        key: fraction / coefficients[key]
        for key, fraction in target_fractions.items()
    }
    extent = max(ratios.values())
    driving = {key for key, ratio in ratios.items() if ratio == extent}

    reactant_blocks: list[dict] = []
    for sp in reactants:
        key = species_key(sp)
        coef = coefficients[key]
        required = extent * coef
        reactant_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "required": fraction_to_decimal_string(required),
            }
        )

    product_blocks: list[dict] = []
    target_blocks: list[dict] = []
    for sp in products:
        key = species_key(sp)
        coef = coefficients[key]
        achieved = extent * coef
        block = {
            "input": sp.raw,
            "normalized": sp.normalized,
            "coefficient": coef,
            "achieved": fraction_to_decimal_string(achieved),
        }
        product_blocks.append(block)
        if key in target_fractions:
            target = target_fractions[key]
            excess = achieved - target
            target_blocks.append(
                {
                    "input": target_raw[key],
                    "normalized": sp.normalized,
                    "coefficient": coef,
                    "target": fraction_to_decimal_string(target),
                    "achieved": fraction_to_decimal_string(achieved),
                    "excess": fraction_to_decimal_string(excess),
                    "is_driving": key in driving,
                }
            )

    return {
        "reaction_extent": {
            "value": fraction_to_decimal_string(extent),
            "unit": "mol",
        },
        "driving_targets": [block["normalized"] for block in target_blocks
                           if block["is_driving"]],
        "reactants": reactant_blocks,
        "targets": target_blocks,
        "products": product_blocks,
    }


# ------------------------------------------------------------ 实际收率
def compute_yields(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    feed: dict[str, Decimal],
    actual: dict[str, Decimal],
    actual_raw: dict[str, str],
) -> dict:
    """计算反应进度与逐种生成物的理论产量、实际产量和收率百分比。

    先完全复用 :func:`compute_stoichiometry` 的反应进度逻辑
    ξ = min(投料量 ÷ 系数)（同一精确有理数，理论产量与批量计量接口逐位
    一致）；每种生成物的理论产量 = ξ·系数，实际产量取自 ``actual``，
    收率 = 实际产量 ÷ 理论产量 × 100（单位：%）。

    ``feed`` 必须覆盖全部反应物、``actual`` 必须覆盖全部生成物，键均为
    :func:`app.cache.species_key`，值非负；``actual_raw`` 以规范键记录
    客户端提交的原始生成物写法，用于 ``input`` 回显。

    理论产量为 0（即 ξ = 0，至少一种反应物投料为 0）时收率百分比无定义，
    ``yield_percent`` 返回 ``None`` 并在 ``yield_note`` 中说明全部零投料
    的反应物；实际产量高于理论值时不做封顶，保留大于 100% 的计算结果。
    物种按规范键排序输出。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)

    feed_fractions = {key: decimal_to_fraction(amount) for key, amount in feed.items()}
    actual_fractions = {
        key: decimal_to_fraction(amount) for key, amount in actual.items()
    }

    extent, limiting = _extent_and_limiting(reactants, coefficients, feed_fractions)

    reactant_blocks: list[dict] = []
    for sp in reactants:
        key = species_key(sp)
        coef = coefficients[key]
        fed = feed_fractions[key]
        consumed = extent * coef
        remaining = fed - consumed
        reactant_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "fed": fraction_to_decimal_string(fed),
                "consumed": fraction_to_decimal_string(consumed),
                "remaining": fraction_to_decimal_string(remaining),
                "is_limiting": key in limiting,
            }
        )

    zero_fed_limiting = [
        block["normalized"]
        for block in reactant_blocks
        if block["is_limiting"] and block["fed"] == "0"
    ]
    note = None
    if extent == 0:
        note = (
            "理论产量为 0：反应进度 ξ = min(投料量 ÷ 系数) = 0，"
            f"投料为 0 的限量反应物为 {zero_fed_limiting}，收率百分比无定义"
        )

    product_blocks: list[dict] = []
    for sp in products:
        key = species_key(sp)
        coef = coefficients[key]
        theoretical = extent * coef
        actual_value = actual_fractions[key]
        if theoretical == 0:
            percent: Fraction | None = None
            product_note = note
        else:
            percent = actual_value * 100 / theoretical
            product_note = None
        product_blocks.append(
            {
                "input": actual_raw.get(key, sp.raw),
                "normalized": sp.normalized,
                "coefficient": coef,
                "theoretical_yield": fraction_to_decimal_string(theoretical),
                "actual_yield": fraction_to_decimal_string(actual_value),
                "yield_percent": (
                    None if percent is None else fraction_to_decimal_string(percent)
                ),
                "yield_note": product_note,
            }
        )

    return {
        "reaction_extent": {
            "value": fraction_to_decimal_string(extent),
            "unit": "mol",
        },
        "limiting_reactants": [block["normalized"] for block in reactant_blocks
                               if block["is_limiting"]],
        "reactants": reactant_blocks,
        "products": product_blocks,
    }


# ------------------------------------------------------------ 多步连续反应推演
def apply_sequence_step(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    inventory: dict[str, Fraction],
) -> dict:
    """以当前库存推演一步唯一配平反应，返回逐步明细与更新后的库存。

    ``coefficients`` 以 :func:`app.cache.species_key` 为键；``inventory``
    为物种规范键 -> 非负 :class:`Fraction` 摩尔数，未列入的物种库存视为 0。

    最大反应进度 ξ = min(当前库存 ÷ 反应物系数)；ξ = 0 时该步以零进度继续，
    库存不发生任何变化，库存为 0（或缺省）的反应物全部记入 ``blocking``。
    否则扣减各反应物 ξ·系数、向各生成物加入 ξ·系数，得到新库存。

    库存全程以精确有理数传递与比较，仅在输出时渲染为十进制字符串。
    返回的 ``inventory_after`` 是新字典，不修改入参 ``inventory``。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)

    # 最大反应进度 ξ = min(库存 / 系数)；未提供库存的反应物按 0 处理
    ratios: dict[str, Fraction] = {
        species_key(sp): inventory.get(species_key(sp), Fraction(0))
        / coefficients[species_key(sp)]
        for sp in reactants
    }
    extent = min(ratios.values())
    blocking = {key for key, ratio in ratios.items() if ratio == 0}

    after = dict(inventory)

    reactant_blocks: list[dict] = []
    for sp in reactants:
        key = species_key(sp)
        coef = coefficients[key]
        before = inventory.get(key, Fraction(0))
        consumed = extent * coef
        after[key] = before - consumed
        reactant_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "before": fraction_to_decimal_string(before),
                "consumed": fraction_to_decimal_string(consumed),
                "after": fraction_to_decimal_string(after[key]),
                "is_blocking": key in blocking,
            }
        )

    product_blocks: list[dict] = []
    for sp in products:
        key = species_key(sp)
        coef = coefficients[key]
        before = inventory.get(key, Fraction(0))
        produced = extent * coef
        after[key] = before + produced
        product_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "before": fraction_to_decimal_string(before),
                "produced": fraction_to_decimal_string(produced),
                "after": fraction_to_decimal_string(after[key]),
            }
        )

    return {
        "reaction_extent": extent,
        "blocking": [
            block["normalized"] for block in reactant_blocks
            if block["is_blocking"]
        ],
        "reactants": reactant_blocks,
        "products": product_blocks,
        "inventory_after": after,
    }
