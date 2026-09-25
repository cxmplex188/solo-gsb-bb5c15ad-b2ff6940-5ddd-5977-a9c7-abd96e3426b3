"""业务编排层：解析、配平、缓存与响应组装。"""

from __future__ import annotations

from datetime import date
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
from typing import Any

from .allocation import (
    BATCH_STATUSES,
    Batch,
    Demand,
    allocate,
    parse_purity,
    render_allocation,
)
from .balancer import (
    VARIANT_SPACE_TOO_LARGE,
    BalanceError,
    SolvedReaction,
    Species,
    build_balanced_view,
    build_raw_view,
    enumerate_variants_parsed,
    parse_side,
    solve_reaction,
)
from .blending import Stock, solve_blend
from .cache import ReactionRecord, build_signature, cache, reaction_id_for, species_key
from .equilibrium import solve_equilibrium
from .formula_infer import (
    ELEMENTS as FORMULA_ELEMENTS,
    ElementMeasurement,
    infer_formulas,
    render_formula,
)
from .hess import HessError, KnownReaction, normalized_by_key, solve_hess
from .isotopes import (
    MAX_ELEMENTS,
    MAX_TOTAL_ATOMS,
    ElementIsotopes,
    Isotope,
    isotope_envelope,
    most_probable_peak,
)
from .parser import ELEMENTS as ELEMENT_SYMBOLS, FormulaError, parse_formula
from .stoichiometry import (
    apply_sequence_step,
    compute_reverse_plan,
    compute_stoichiometry,
    compute_yields,
    decimal_to_fraction,
    fraction_to_decimal_string,
    parse_enthalpy,
    parse_moles,
    parse_positive_moles,
)


# 错误码 -> HTTP 状态码
STATUS_CODES = {
    "FORMULA_PARSE_ERROR": 400,
    "DUPLICATE_SPECIES": 409,
    "INVALID_REQUEST": 422,
    "IMPOSSIBLE_BALANCE": 422,
    "NON_UNIQUE_BALANCE": 422,
    "ZERO_COEFFICIENT": 422,
    "INVALID_STOICHIOMETRY": 422,
    "CHARGE_NOT_CONSERVED": 422,
    VARIANT_SPACE_TOO_LARGE: 422,
    "INVALID_INVENTORY": 422,
    "INVALID_BLEND": 422,
    "NO_FEASIBLE_BLEND": 422,
    "INVALID_ISOTOPE": 422,
    "INVALID_FORMULA_INFER": 422,
    "HESS_UNDERIVABLE": 422,
    "HESS_AMBIGUOUS": 422,
}


def error_payload(code: str, message: str, details: dict[str, Any]) -> dict[str, Any]:
    return {"code": code, "message": message, "details": details}


def _parse_all(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[list[Species], list[Species]]:
    try:
        return list(parse_side(reactant_texts, "反应物")), list(
            parse_side(product_texts, "生成物")
        )
    except FormulaError as exc:
        raise BalanceError("FORMULA_PARSE_ERROR", str(exc)) from exc


def _equation(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> str:
    def term(sp: Species) -> str:
        coef = coefficient_map.get(species_key(sp), 1) if coefficient_map else 1
        return ("" if coef == 1 else str(coef)) + sp.normalized

    return " + ".join(map(term, reactants)) + " -> " + " + ".join(map(term, products))


def _normalized_block(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> dict[str, Any]:
    return {
        "reactants": [sp.normalized for sp in reactants],
        "products": [sp.normalized for sp in products],
        "equation": _equation(reactants, products, coefficient_map),
    }


def _species_coefficients(
    side: list[Species], coefficient_map: dict[str, int]
) -> list[dict[str, Any]]:
    return [
        {
            "input": sp.raw,
            "normalized": sp.normalized,
            "coefficient": coefficient_map[species_key(sp)],
        }
        for sp in side
    ]


def _coef_map(solved) -> dict[str, int]:
    return {
        species_key(sp): coef
        for sp, coef in zip(solved.all_species, solved.coefficients)
    }


# ------------------------------------------------------------ 解析失败处理
def _text_signature(reactant_texts: list[str], product_texts: list[str]) -> str:
    """解析失败时的退化幂等键：压缩空白后排序。"""
    def norm(texts: list[str]) -> list[str]:
        return sorted("".join(t.split()) for t in texts)

    return f"raw:r={norm(reactant_texts)!r};p={norm(product_texts)!r}"


def _record_parse_error(
    reactant_texts: list[str], product_texts: list[str], exc: BalanceError
) -> tuple[ReactionRecord, bool]:
    signature = _text_signature(reactant_texts, product_texts)
    existing = cache.get_by_signature(signature)
    if existing is not None:
        return existing, True
    record = cache.put(
        signature, [], [], None, [],
        error_payload(exc.code, exc.message, dict(exc.details)),
    )
    return record, False


# ------------------------------------------------------------ 配平主流程
def process_balance(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """配平入口，返回 (http_status, payload)。成功与错误结果都按规范签名缓存。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    if record.error is not None:
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }
    return 200, build_balance_payload(record, reactants, products, cached)


def build_balance_payload(
    record: ReactionRecord,
    reactants: list[Species],
    products: list[Species],
    cached: bool,
) -> dict[str, Any]:
    coefficient_map = record.coefficients
    charge_left = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in reactants
    )
    charge_right = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in products
    )
    return {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(reactants, products, coefficient_map),
        "coefficients": {
            "reactants": _species_coefficients(reactants, coefficient_map),
            "products": _species_coefficients(products, coefficient_map),
            "vector": [
                coefficient_map[species_key(sp)] for sp in (*reactants, *products)
            ],
        },
        "balanced": True,
        "elements": list(record.elements),
        "charge_balance": {
            "reactant_total": charge_left,
            "product_total": charge_right,
            "difference": charge_left - charge_right,
            "balanced": charge_left == charge_right,
        },
    }


# ------------------------------------------------------------ 逐元素分析
def process_analyze(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """逐元素差异分析。可解析即返回 200（含原始/平衡两个视图）；无法解析返回 400。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "reaction_id": record.reaction_id,
            "cached": cached,
            "error": record.error,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None
    elements = sorted(
        {el for sp in (*reactants, *products) for el in sp.formula.composition}
    )

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    raw = build_raw_view(reactants, products, elements)
    payload: dict[str, Any] = {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(
            reactants, products, None if record.error else record.coefficients
        ),
        "raw": _view(reactants, products, raw, coefficient_map=None, assumed=1),
        "balanced_ok": record.error is None,
        "error": record.error,
    }
    if record.error is None:
        solved = solve_reaction(reactant_texts, product_texts)
        balanced_view = build_balanced_view(solved)
        payload["balanced"] = _view(
            list(solved.reactants),
            list(solved.products),
            balanced_view,
            coefficient_map=record.coefficients,
            assumed=None,
        )
    else:
        payload["balanced"] = None
    return 200, payload


def _view(
    reactants: list[Species],
    products: list[Species],
    view: dict[str, Any],
    *,
    coefficient_map: dict[str, int] | None,
    assumed: int | None,
) -> dict[str, Any]:
    def side_block(side: list[Species]) -> list[dict[str, Any]]:
        block = []
        for sp in side:
            coef = (
                coefficient_map[species_key(sp)] if coefficient_map else (assumed or 1)
            )
            block.append(
                {"input": sp.raw, "normalized": sp.normalized, "coefficient": coef}
            )
        return block

    result: dict[str, Any] = {
        "species": {
            "reactants": side_block(reactants),
            "products": side_block(products),
        },
        "elements": view["elements"],
        "charge": view["charge"],
    }
    if assumed is not None:
        result["coefficients_assumed"] = assumed
    return result


def get_reaction(reaction_id: str) -> dict[str, Any] | None:
    """按 reaction_id 取回缓存结果（幂等查询）。"""
    record = cache.get(reaction_id)
    if record is None:
        return None
    reactants = _record_species(record, "canonical_reactants")
    products = _record_species(record, "canonical_products")
    if record.error is not None:
        return {
            "reaction_id": record.reaction_id,
            "cached": True,
            "error": record.error,
        }
    return build_balance_payload(record, reactants, products, True)


def _record_species(record: ReactionRecord, attr: str) -> list[Species]:
    """缓存里只存了 normalized 文本；查询回显时构造轻量 Species。"""
    from .parser import parse_formula

    return [
        Species(
            raw=item["normalized"],
            normalized=item["normalized"],
            formula=parse_formula(item["normalized"]),
        )
        for item in getattr(record, attr)
    ]


# ------------------------------------------------------------ 多变体枚举
def process_variants(
    reactant_texts: list[str],
    product_texts: list[str],
    max_coefficient: int,
    limit: int,
) -> tuple[int, dict[str, Any]]:
    """枚举配平变体，返回 (http_status, payload)。结果不写幂等缓存。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        # 解析 / 去重阶段失败：沿用与 /balance 一致的退化签名生成 ID
        reaction_id = reaction_id_for(
            _text_signature(reactant_texts, product_texts)
        )
        return STATUS_CODES[exc.code], {
            "error": error_payload(exc.code, exc.message, dict(exc.details)),
            "reaction_id": reaction_id,
        }

    try:
        variant_set = enumerate_variants_parsed(reactants, products, max_coefficient)
    except BalanceError as exc:
        reaction_id = reaction_id_for(build_signature(reactants, products))
        return STATUS_CODES[exc.code], {
            "error": error_payload(exc.code, exc.message, dict(exc.details)),
            "reaction_id": reaction_id,
        }

    signature = build_signature(reactants, products)
    reaction_id = reaction_id_for(signature)

    total = len(variant_set.coefficients)
    page = variant_set.coefficients[:limit]

    items: list[dict[str, Any]] = []
    for vec in page:
        coefficient_map = {
            species_key(sp): coef
            for sp, coef in zip((*reactants, *products), vec)
        }
        items.append(
            {
                "coefficients": {
                    "reactants": _species_coefficients(reactants, coefficient_map),
                    "products": _species_coefficients(products, coefficient_map),
                    "vector": list(vec),
                },
                "sum": sum(vec),
                "equation": _equation(reactants, products, coefficient_map),
            }
        )

    return 200, {
        "reaction_id": reaction_id,
        "normalized": _normalized_block(reactants, products, None),
        "max_coefficient": max_coefficient,
        "limit": limit,
        "nullity": variant_set.nullity,
        "count": len(page),
        "total": total,
        "truncated": total > limit,
        "variants": items,
    }


# ------------------------------------------------------------ 限量试剂批量计量
def _parse_species_moles(
    items: list[dict[str, str]], field_name: str
) -> tuple[
    list[tuple[str, Species]],
    list[tuple[str, Species]],
    list[dict[str, str]],
    list[str],
]:
    """把 ``[{species, moles}]`` 条目解析为精确 Decimal，聚合可定位的问题。

    返回 ``(有效条目, 可解析物种条目, 非法摩尔条目, 重复规范化物种)``：

    * 有效条目 ``(原始摩尔文本, 物种, Decimal 摩尔数)``：物种可解析且摩尔数合法；
    * 可解析物种条目 ``(原始摩尔文本, 物种)`` 还包含摩尔数非法的记录，
      供匹配未知物种/缺失物种时不遗漏；
    * 物种本身无法解析时抛出 :class:`BalanceError`（FORMULA_PARSE_ERROR），
      由调用方按 ``field_name`` 组装 400 响应；
    * 摩尔数非法 / 物种重复不在此报错，由调用方聚合成 422 响应。
    """
    parsed_entries: list[tuple[str, Species, Decimal]] = []
    parseable_entries: list[tuple[str, Species]] = []
    invalid_moles: list[dict[str, str]] = []
    duplicate: list[str] = []
    seen_keys: set[str] = set()

    for item in items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            raise BalanceError(
                "FORMULA_PARSE_ERROR",
                str(exc),
                field=field_name,
                input=species_text,
            ) from exc
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        parseable_entries.append((moles_text, sp))

        key = species_key(sp)
        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            parsed_entries.append((moles_text, sp, amount))

    return parsed_entries, parseable_entries, invalid_moles, duplicate


def process_stoichiometry(
    reactant_texts: list[str],
    product_texts: list[str],
    amount_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """限量试剂批量计量，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的唯一配平、规范化与幂等缓存逻辑；
    配平失败（解析错误、重复物种、不可唯一配平等）沿用其结构化错误响应。
    随后校验投料与规范化反应物一一对应，并以精确十进制计算反应进度、
    并列限量试剂、各反应物消耗/剩余与各生成物理论产量。
    """
    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # 规范化反应物（parse_side 已按规范键排序，与配平向量顺序一致）
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    reactant_keys = {species_key(sp) for sp in reactants}

    # ---- 投料条目：逐条解析物种与摩尔数，聚合所有可定位的问题
    parsed_entries: list[tuple[str, Species]] = []
    invalid_moles: list[dict[str, str]] = []
    feed: dict[str, Decimal] = {}
    duplicate: list[str] = []
    seen_keys: set[str] = set()

    for item in amount_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            # 投料物种本身无法解析：沿用化学式解析错误（400）
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "reaction_id": reaction_id,
                "cached": cached,
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "amounts.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)

        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            parsed_entries.append((moles_text, sp))
            feed.setdefault(key, amount)

    unknown = sorted(
        {sp.normalized for _, sp in parsed_entries if species_key(sp) not in reactant_keys}
    )
    provided_keys = {species_key(sp) for _, sp in parsed_entries}
    normalized_by_key = {species_key(sp): sp.normalized for sp in reactants}
    missing = [
        normalized_by_key[key] for key in reactant_keys if key not in provided_keys
    ]

    if duplicate or unknown or missing or invalid_moles:
        details: dict[str, Any] = {}
        if missing:
            details["missing"] = sorted(missing)
        if duplicate:
            details["duplicate"] = sorted(set(duplicate))
        if unknown:
            details["unknown"] = unknown
        if invalid_moles:
            details["invalid_moles"] = invalid_moles
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "投料记录必须与规范化反应物一一对应：不能有缺失、重复或未知物种，"
                "且每条摩尔数必须是非负十进制字符串（不能为 JSON 数字）",
                details,
            ),
        }

    result = compute_stoichiometry(reactants, products, coefficient_map, feed)
    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        **result,
    }


# ------------------------------------------------------------ 实际收率
def process_yield(
    reactant_texts: list[str],
    product_texts: list[str],
    amount_items: list[dict[str, str]],
    actual_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """实际收率计算，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的唯一配平、规范化与幂等缓存逻辑；
    配平失败沿用其结构化错误响应。随后校验投料与规范化反应物一一对应、
    实际产量与规范化生成物一一对应（语义同
    :func:`process_stoichiometry`，错误明细分别以 ``amounts.`` /
    ``actual_yields.`` 为前缀），再复用 :func:`compute_stoichiometry`
    的反应进度精确计算，逐种给出理论产量、实际产量与收率百分比；
    理论产量为 0 时收率返回 null 并说明零投料的限量反应物，实际产量
    高于理论值时保留大于 100% 的结果。
    """
    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # parse_side 已按规范键排序，与配平向量顺序一致
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    reactant_keys = {species_key(sp) for sp in reactants}
    product_keys = {species_key(sp) for sp in products}

    try:
        feed_entries, feed_parseable, feed_invalid, feed_duplicate = (
            _parse_species_moles(amount_items, "amounts.species")
        )
        actual_entries, actual_parseable, actual_invalid, actual_duplicate = (
            _parse_species_moles(actual_items, "actual_yields.species")
        )
    except BalanceError as exc:
        return STATUS_CODES[exc.code], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(exc.code, exc.message, dict(exc.details)),
        }

    # ---- 投料侧（与 /stoichiometry/batch 相同的一一对应语义）
    feed_unknown = sorted(
        {sp.normalized for _, sp in feed_parseable
         if species_key(sp) not in reactant_keys}
    )
    feed_provided = {species_key(sp) for _, sp, _ in feed_entries}
    reactant_names = {species_key(sp): sp.normalized for sp in reactants}
    feed_missing = [
        reactant_names[key] for key in reactant_keys if key not in feed_provided
    ]

    # ---- 实际产量侧
    actual_unknown = sorted(
        {sp.normalized for _, sp in actual_parseable
         if species_key(sp) not in product_keys}
    )
    actual_provided = {species_key(sp) for _, sp, _ in actual_entries}
    product_names = {species_key(sp): sp.normalized for sp in products}
    actual_missing = [
        product_names[key] for key in product_keys if key not in actual_provided
    ]

    if (
        feed_duplicate or feed_unknown or feed_missing or feed_invalid
        or actual_duplicate or actual_unknown or actual_missing or actual_invalid
    ):
        details: dict[str, Any] = {}
        if feed_missing:
            details["missing"] = sorted(feed_missing)
        if feed_duplicate:
            details["duplicate"] = sorted(set(feed_duplicate))
        if feed_unknown:
            details["unknown"] = feed_unknown
        if feed_invalid:
            details["invalid_moles"] = feed_invalid
        if actual_missing:
            details["actual_yields.missing"] = sorted(actual_missing)
        if actual_duplicate:
            details["actual_yields.duplicate"] = sorted(set(actual_duplicate))
        if actual_unknown:
            details["actual_yields.unknown"] = actual_unknown
        if actual_invalid:
            details["actual_yields.invalid_moles"] = actual_invalid
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "投料记录必须与规范化反应物一一对应，实际产量记录必须与规范化"
                "生成物一一对应：均不能有缺失、重复或未知物种，且每条摩尔数必须"
                "是非负十进制字符串（不能为 JSON 数字）",
                details,
            ),
        }

    feed: dict[str, Decimal] = {}
    for _moles_text, sp, amount in feed_entries:
        feed[species_key(sp)] = amount
    actual: dict[str, Decimal] = {}
    actual_raw: dict[str, str] = {}
    for _moles_text, sp, amount in actual_entries:
        actual[species_key(sp)] = amount
        actual_raw[species_key(sp)] = sp.raw

    result = compute_yields(
        reactants, products, coefficient_map, feed, actual, actual_raw
    )
    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        **result,
    }


# ------------------------------------------------------------ 反向投料规划
def process_reverse_plan(
    reactant_texts: list[str],
    product_texts: list[str],
    target_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """反向投料规划，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的唯一配平、规范化与幂等缓存逻辑；
    配平失败沿用其结构化错误响应。随后校验目标条目都能映射到规范化
    生成物（允许只提交部分生成物，但不能重复或出现未知物种），并以精确
    十进制计算最小反应进度 ξ = max(目标量 ÷ 系数)、各反应物理论需求量、
    各目标产物达成量/超额量与全部生成物达成量。
    """
    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # parse_side 已按规范键排序，与配平向量顺序一致
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    product_keys = {species_key(sp) for sp in products}

    # ---- 目标条目：逐条解析物种与摩尔数，聚合所有可定位的问题
    parsed_entries: list[tuple[Species, Decimal]] = []
    invalid_moles: list[dict[str, str]] = []
    duplicate: list[str] = []
    seen_keys: set[str] = set()

    for item in target_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            # 目标物种本身无法解析：沿用化学式解析错误（400）
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "reaction_id": reaction_id,
                "cached": cached,
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "targets.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)

        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            parsed_entries.append((sp, amount))

    unknown = sorted(
        {sp.normalized for sp, _ in parsed_entries if species_key(sp) not in product_keys}
    )

    if duplicate or unknown or invalid_moles:
        details: dict[str, Any] = {}
        if duplicate:
            details["duplicate"] = sorted(set(duplicate))
        if unknown:
            details["unknown"] = unknown
        if invalid_moles:
            details["invalid_moles"] = invalid_moles
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "目标产物必须是规范化生成物之一：不能重复或包含未知物种"
                "（允许只提交部分生成物），且每条目标摩尔数必须是非负十进制"
                "字符串（不能为 JSON 数字）",
                details,
            ),
        }

    result = compute_reverse_plan(
        reactants, products, coefficient_map, parsed_entries
    )
    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        **result,
    }


# ------------------------------------------------------------ 理想溶液平衡组成
def _fraction_to_decimal(value: Fraction) -> Decimal:
    """把正有理数近似转换为 :class:`Decimal`（30 位有效数字）。

    仅用于向 :func:`fraction_to_decimal_string` 传递绝对精度阈值
    （该函数只取阈值的数量级决定追加多少有效位），不参与计算结果本身。
    """
    with localcontext() as ctx:
        ctx.prec = 30
        return Decimal(value.numerator) / Decimal(value.denominator)


def process_equilibrium(
    reactant_texts: list[str],
    product_texts: list[str],
    amount_items: list[dict[str, str]],
    volume_text: str,
    kc_text: str,
) -> tuple[int, dict[str, Any]]:
    """理想溶液平衡组成求解，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的解析、唯一配平与幂等缓存（解析失败
    400、同侧重复 409、不可唯一配平 422）；随后校验同一规范化物种不得
    同时出现在反应两侧，初始摩尔数记录必须与规范化的反应物 ∪ 生成物
    一一对应且为**正**十进制字符串，体积 V（L）与 Kc 同样必须为正十进制
    字符串，问题聚合为 422 INVALID_REQUEST。校验通过后以
    :func:`app.equilibrium.solve_equilibrium` 二分求 Qc=Kc 的唯一内点
    进度（绝对误差上界 5e-10 mol），返回进度、非负物料区间、初始/平衡
    反应商与各物种初始/平衡摩尔数。
    """
    # ---- 配平前预检：两侧重复物种对可逆反应永远非法（否则会被报成零系数，
    # 错误语义不直观）。任一侧无法解析时跳过预检，交由 process_balance 统一报 400。
    pre_reactants: list[Species] | None = None
    pre_products: list[Species] | None = None
    try:
        pre_reactants = list(parse_side(reactant_texts, "反应物"))
        pre_products = list(parse_side(product_texts, "生成物"))
    except (BalanceError, FormulaError):
        pass
    if pre_reactants is not None and pre_products is not None:
        product_key_set = {species_key(sp) for sp in pre_products}
        cross_duplicate = sorted(
            sp.normalized
            for sp in pre_reactants
            if species_key(sp) in product_key_set
        )
        if cross_duplicate:
            return STATUS_CODES["INVALID_REQUEST"], {
                "error": error_payload(
                    "INVALID_REQUEST",
                    "同一规范化物种不得同时出现在反应两侧：可逆反应要求每个"
                    "参与物种只属于反应物或生成物一侧",
                    {"cross_side_duplicate": cross_duplicate},
                ),
            }

    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # parse_side 已按规范键排序，与配平向量顺序一致
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    reactant_keys = {species_key(sp) for sp in reactants}
    product_keys = {species_key(sp) for sp in products}
    species_keys_all = reactant_keys | product_keys

    # ---- 初始摩尔条目：逐条解析物种与摩尔数，聚合全部可定位的问题
    invalid_moles: list[dict[str, str]] = []
    initial_fractions: dict[str, Fraction] = {}
    initial_raw: dict[str, str] = {}
    duplicate: list[str] = []
    unknown_species: list[str] = []
    seen_keys: set[str] = set()

    for item in amount_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "reaction_id": reaction_id,
                "cached": cached,
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "amounts.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)

        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)
        if key not in species_keys_all:
            unknown_species.append(formula.normalized)

        try:
            amount = parse_positive_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            initial_fractions[key] = decimal_to_fraction(amount)
            initial_raw.setdefault(key, species_text)

    all_species = [*reactants, *products]
    normalized_by_key = {species_key(sp): sp.normalized for sp in all_species}
    missing = [
        normalized_by_key[key] for key in species_keys_all if key not in seen_keys
    ]
    unknown = sorted(set(unknown_species))

    # ---- 体积与 Kc：严格为正的有限十进制字符串
    invalid_volume: str | None = None
    volume_fraction: Fraction | None = None
    try:
        volume_fraction = decimal_to_fraction(parse_positive_moles(volume_text))
    except (ValueError, InvalidOperation):
        invalid_volume = volume_text

    invalid_kc: str | None = None
    kc_fraction: Fraction | None = None
    try:
        kc_fraction = decimal_to_fraction(parse_positive_moles(kc_text))
    except (ValueError, InvalidOperation):
        invalid_kc = kc_text

    if (
        cross_duplicate or duplicate or unknown or missing or invalid_moles
        or invalid_volume is not None or invalid_kc is not None
    ):
        details: dict[str, Any] = {}
        if cross_duplicate:
            details["cross_side_duplicate"] = cross_duplicate
        if missing:
            details["missing"] = sorted(missing)
        if duplicate:
            details["duplicate"] = sorted(set(duplicate))
        if unknown:
            details["unknown"] = unknown
        if invalid_moles:
            details["invalid_moles"] = invalid_moles
        if invalid_volume is not None:
            details["invalid_volume"] = invalid_volume
        if invalid_kc is not None:
            details["invalid_kc"] = invalid_kc
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "初始摩尔数记录必须与规范化的反应物与生成物一一对应（不能缺失、"
                "重复或出现未知物种，同一物种不得同时出现在两侧），各初始摩尔数、"
                "体积 volume_l 与平衡常数 Kc 均必须是严格为正的十进制字符串"
                "（不能为 JSON 数字、0 或负数）",
                details,
            ),
        }

    solution = solve_equilibrium(
        reactants,
        products,
        coefficient_map,
        initial_fractions,
        volume_fraction,
        kc_fraction,
    )

    # 进度与摩尔数渲染精度：贴近物料区间端点时 Qc 对 ξ 的对数斜率
    # S = dlnQc/dξ 极大（~ 计量系数 ÷ 端点残余量），固定绝对精度会把
    # 微量剩余抹掉，使按响应数值复算的反应商严重偏离 Kc。取
    # min(1e-20, 1e-12/S)：渲染引入的 Qc 相对漂移不超过 ~1e-12，
    # 与求解器 1e-9 的内部误差合并后仍远严于规范要求的 1e-6；
    # 常规情形（S < 1e8）退化为既有的 1e-20 绝对精度，返回进度的
    # 总绝对误差上界 < 5e-10 + 1e-20 mol（规范要求 < 1e-9 mol）
    base_accuracy = Fraction(1, 10**20)
    slope_accuracy = Fraction(1, 10**12) / solution.qc_log_slope
    mole_accuracy = _fraction_to_decimal(min(base_accuracy, slope_accuracy))

    def render(value: Fraction) -> str:
        return fraction_to_decimal_string(
            value, min_absolute_accuracy=mole_accuracy
        )

    def species_block(side: list[Species]) -> list[dict[str, Any]]:
        blocks = []
        for sp in side:
            key = species_key(sp)
            blocks.append(
                {
                    "input": initial_raw.get(key, sp.raw),
                    "normalized": sp.normalized,
                    "coefficient": coefficient_map[key],
                    "initial_moles": fraction_to_decimal_string(initial_fractions[key]),
                    "equilibrium_moles": render(solution.equilibrium_moles[key]),
                }
            )
        return blocks

    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        "volume_l": volume_text,
        "kc": kc_text,
        "direction": solution.direction,
        "reaction_extent": {
            "value": render(solution.extent),
            "unit": "mol",
        },
        "extent_bounds": {
            "minimum": fraction_to_decimal_string(solution.extent_min),
            "maximum": fraction_to_decimal_string(solution.extent_max),
        },
        "initial_qc": fraction_to_decimal_string(solution.initial_qc),
        "equilibrium_qc": fraction_to_decimal_string(solution.equilibrium_qc),
        "reactants": species_block(reactants),
        "products": species_block(products),
    }


# ------------------------------------------------------------ 多步连续反应推演
def _step_error(
    step_index: int, payload: dict[str, Any], status: int
) -> tuple[int, dict[str, Any]]:
    """给单步配平的结构化错误补充步骤序号（不改动缓存内的错误对象）。"""
    body = dict(payload)
    error = dict(body["error"])
    details = dict(error.get("details") or {})
    details["step"] = step_index
    error["details"] = details
    body["error"] = error
    return status, body


def _render_inventory(
    inventory: dict[str, Fraction], normalized_by_key: dict[str, str]
) -> list[dict[str, str]]:
    """把规范键 -> Fraction 的库存渲染为按规范化化学式排序的记录列表。"""
    return [
        {
            "normalized": normalized_by_key[key],
            "moles": fraction_to_decimal_string(inventory[key]),
        }
        for key in sorted(inventory, key=lambda k: normalized_by_key[k])
    ]


def process_sequence(
    step_items: list[dict[str, Any]],
    inventory_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """多步连续反应推演，返回 (http_status, payload)。

    先把初始库存条目逐条解析为精确有理数（等价化学式合并为同一物种、
    摩尔数相加）；随后按顺序对每一步复用 :func:`process_balance` 的解析、
    唯一配平与幂等缓存，以该步全部反应物的当前库存计算最大反应进度，
    扣减消耗、加入产物，并把更新后的库存传给下一步。任何一步配平失败
    都沿用其结构化错误响应并在 details 中补充 ``step`` 序号；库存物种
    无法解析或摩尔数非法时在推演开始前聚合报 400 / 422。
    """
    # ---- 初始库存：逐条解析物种与摩尔数，聚合所有可定位的问题
    invalid_moles: list[dict[str, str]] = []
    inventory: dict[str, Fraction] = {}
    normalized_by_key: dict[str, str] = {}

    for item in inventory_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "initial_inventory.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)
        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append(
                {"species": formula.normalized, "moles": moles_text}
            )
            continue
        # 等价化学式合并为同一库存物种，摩尔数精确相加
        inventory[key] = inventory.get(key, Fraction(0)) + decimal_to_fraction(amount)
        normalized_by_key.setdefault(key, formula.normalized)

    if invalid_moles:
        return STATUS_CODES["INVALID_REQUEST"], {
            "error": error_payload(
                "INVALID_REQUEST",
                "初始库存摩尔数必须是非负十进制字符串（不能为 JSON 数字）",
                {"invalid_moles": invalid_moles},
            ),
        }

    # ---- 逐步解析、配平与库存推演
    step_results: list[dict[str, Any]] = []
    for index, step in enumerate(step_items):
        reactant_texts = step["reactants"]
        product_texts = step["products"]
        status, balance_payload = process_balance(reactant_texts, product_texts)
        if status != 200:
            return _step_error(index, balance_payload, status)

        reaction_id = balance_payload["reaction_id"]
        cached = balance_payload["cached"]

        # parse_side 已按规范键排序，与配平向量顺序一致
        reactants = list(parse_side(reactant_texts, "反应物"))
        products = list(parse_side(product_texts, "生成物"))
        coefficient_map = {
            species_key(sp): coef
            for sp, coef in zip(
                (*reactants, *products),
                balance_payload["coefficients"]["vector"],
            )
        }
        for sp in (*reactants, *products):
            normalized_by_key.setdefault(species_key(sp), sp.normalized)

        outcome = apply_sequence_step(
            reactants, products, coefficient_map, inventory
        )
        inventory = outcome.pop("inventory_after")
        extent = outcome.pop("reaction_extent")

        step_results.append(
            {
                "index": index,
                "reaction_id": reaction_id,
                "cached": cached,
                "normalized": balance_payload["normalized"],
                "reaction_extent": {
                    "value": fraction_to_decimal_string(extent),
                    "unit": "mol",
                },
                "blocked": extent == 0,
                "blocking_reactants": outcome["blocking"],
                "reactants": outcome["reactants"],
                "products": outcome["products"],
                "inventory": _render_inventory(inventory, normalized_by_key),
            }
        )

    return 200, {
        "steps": step_results,
        "final_inventory": _render_inventory(inventory, normalized_by_key),
    }


# ------------------------------------------------------------ 试剂批次领用分配
def process_allocate(
    demand_items: list[dict[str, Any]],
    batch_items: list[dict[str, Any]],
    accounting_date: date,
) -> tuple[int, dict[str, Any]]:
    """试剂批次 FEFO 领用分配，返回 (http_status, payload)。

    逐条解析领用需求（物种 + 非负有效摩尔）与库存批次（物种、0<纯度<=1、
    非负名义摩尔、失效日期、状态），聚合全部可定位的输入问题：物种无法
    解析沿用 400 FORMULA_PARSE_ERROR；需求物种重复、批号重复、纯度或摩尔
    数非法、状态未知则汇总为 422 INVALID_INVENTORY。校验通过后调用
    :func:`app.allocation.allocate` 按失效日期升序、同日批号升序分配，
    过期/隔离/零可用量批次跳过，结果以十进制字符串渲染。
    """
    # ---- 领用需求
    demands: list[Demand] = []
    invalid_demand_moles: list[dict[str, str]] = []
    duplicate_demands: list[str] = []
    seen_demand_keys: set[str] = set()

    for item in demand_items:
        species_text = item["species"]
        moles_text = item["effective_moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "demands.species", "input": species_text},
                ),
            }
        key = formula.key
        normalized = formula.normalized

        if key in seen_demand_keys:
            duplicate_demands.append(normalized)
        seen_demand_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_demand_moles.append(
                {"species": normalized, "effective_moles": moles_text}
            )
        else:
            demands.append(
                Demand(
                    raw=species_text,
                    normalized=normalized,
                    key=key,
                    moles=decimal_to_fraction(amount),
                )
            )

    # ---- 库存批次
    batches: list[Batch] = []
    invalid_purity: list[dict[str, str]] = []
    invalid_nominal: list[dict[str, str]] = []
    invalid_status: list[dict[str, str]] = []
    duplicate_batch_nos: list[str] = []
    seen_batch_nos: set[str] = set()

    for item in batch_items:
        species_text = item["species"]
        batch_no = item["batch_no"]
        purity_text = item["purity"]
        nominal_text = item["nominal_moles"]
        status = item["status"]

        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "batches.species", "batch_no": batch_no,
                     "input": species_text},
                ),
            }

        if batch_no in seen_batch_nos:
            duplicate_batch_nos.append(batch_no)
        seen_batch_nos.add(batch_no)

        if status not in BATCH_STATUSES:
            invalid_status.append({"batch_no": batch_no, "status": str(status)})

        try:
            purity = parse_purity(purity_text)
        except (ValueError, InvalidOperation):
            invalid_purity.append({"batch_no": batch_no, "purity": purity_text})
            purity = None

        try:
            nominal = parse_moles(nominal_text)
        except (ValueError, InvalidOperation):
            invalid_nominal.append(
                {"batch_no": batch_no, "nominal_moles": nominal_text}
            )
            nominal = None

        if purity is not None and nominal is not None:
            batches.append(
                Batch(
                    batch_no=batch_no,
                    species_raw=species_text,
                    normalized=formula.normalized,
                    key=formula.key,
                    purity=decimal_to_fraction(purity),
                    nominal_moles=decimal_to_fraction(nominal),
                    expires_on=item["expires_on"],
                    status=status,
                )
            )

    # ---- 聚合所有输入问题
    if (
        duplicate_demands
        or invalid_demand_moles
        or duplicate_batch_nos
        or invalid_purity
        or invalid_nominal
        or invalid_status
    ):
        details: dict[str, Any] = {}
        if duplicate_demands:
            details["duplicate_demand_species"] = sorted(set(duplicate_demands))
        if invalid_demand_moles:
            details["invalid_effective_moles"] = invalid_demand_moles
        if duplicate_batch_nos:
            details["duplicate_batch_nos"] = sorted(set(duplicate_batch_nos))
        if invalid_purity:
            details["invalid_purity"] = invalid_purity
        if invalid_nominal:
            details["invalid_nominal_moles"] = invalid_nominal
        if invalid_status:
            details["invalid_status"] = invalid_status
        return STATUS_CODES["INVALID_INVENTORY"], {
            "error": error_payload(
                "INVALID_INVENTORY",
                "领用需求与库存批次必须合法：需求物种与批号均不得重复，"
                "有效/名义摩尔数必须是非负十进制字符串（不能为 JSON 数字），"
                "纯度必须是 (0, 1] 区间内的十进制字符串，状态只能为 "
                "available 或 quarantined",
                details,
            ),
        }

    result = allocate(demands, batches, accounting_date)
    return 200, render_allocation(result)


# ------------------------------------------------------------ 溶液配制（母液混合）
def process_blend(
    stock_items: list[dict[str, Any]],
    target_volume_text: str,
    target_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """溶液配制（母液混合），返回 (http_status, payload)。

    逐条解析母液（标识、各溶质浓度、可用体积）与目标（终体积、各溶质目标
    浓度），聚合全部可定位的输入问题为 422 INVALID_BLEND（标识/溶质重复、
    浓度或体积非法、终体积非正）；随后调用
    :func:`app.blending.solve_blend` 以精确有理数求满足全部目标浓度、
    可用体积与总量预算的母液取用量（先母液总量最小、再按输入顺序字典序
    最小），余量用水补足；无解时返回 422 NO_FEASIBLE_BLEND 并附逐溶质
    供应诊断。该接口为无状态计算，不读写幂等缓存。
    """
    # ---- 母液：标识唯一性、溶质浓度与可用体积
    stocks: list[Stock] = []
    invalid_concentrations: list[dict[str, str]] = []
    invalid_volumes: list[dict[str, str]] = []
    duplicate_stock_ids: list[str] = []
    duplicate_stock_solutes: list[dict[str, str]] = []
    seen_ids: set[str] = set()

    for item in stock_items:
        stock_id = item["id"]
        if stock_id in seen_ids:
            duplicate_stock_ids.append(stock_id)
        seen_ids.add(stock_id)

        concentrations: dict[str, Fraction] = {}
        seen_solutes: set[str] = set()
        solutes_valid = True
        for entry in item["solutes"]:
            solute = entry["solute"]
            text = entry["concentration"]
            if solute in seen_solutes:
                duplicate_stock_solutes.append(
                    {"stock_id": stock_id, "solute": solute}
                )
            seen_solutes.add(solute)
            try:
                value = parse_moles(text)
            except (ValueError, InvalidOperation):
                invalid_concentrations.append(
                    {"stock_id": stock_id, "solute": solute, "concentration": text}
                )
                solutes_valid = False
            else:
                concentrations.setdefault(solute, decimal_to_fraction(value))

        try:
            available = parse_moles(item["available_volume"])
        except (ValueError, InvalidOperation):
            invalid_volumes.append(
                {"stock_id": stock_id, "available_volume": item["available_volume"]}
            )
            available = None

        if solutes_valid and available is not None:
            stocks.append(
                Stock(
                    id=stock_id,
                    concentrations=concentrations,
                    available_volume=decimal_to_fraction(available),
                )
            )

    # ---- 目标：终体积与各溶质目标浓度
    targets: dict[str, Fraction] = {}
    invalid_target_concentrations: list[dict[str, str]] = []
    duplicate_target_solutes: list[str] = []
    seen_target_solutes: set[str] = set()

    for entry in target_items:
        solute = entry["solute"]
        text = entry["concentration"]
        if solute in seen_target_solutes:
            duplicate_target_solutes.append(solute)
        seen_target_solutes.add(solute)
        try:
            value = parse_moles(text)
        except (ValueError, InvalidOperation):
            invalid_target_concentrations.append(
                {"solute": solute, "concentration": text}
            )
        else:
            targets.setdefault(solute, decimal_to_fraction(value))

    try:
        final_volume = parse_moles(target_volume_text)
    except (ValueError, InvalidOperation):
        final_volume = None
    # 终体积为 0 时最终浓度无定义，要求严格为正
    invalid_target_volume = final_volume is None or final_volume == 0

    # ---- 聚合所有输入问题
    if (
        duplicate_stock_ids
        or duplicate_stock_solutes
        or invalid_concentrations
        or invalid_volumes
        or duplicate_target_solutes
        or invalid_target_concentrations
        or invalid_target_volume
    ):
        details: dict[str, Any] = {}
        if duplicate_stock_ids:
            details["duplicate_stock_ids"] = sorted(set(duplicate_stock_ids))
        if duplicate_stock_solutes:
            details["duplicate_stock_solutes"] = duplicate_stock_solutes
        if invalid_concentrations:
            details["invalid_concentrations"] = invalid_concentrations
        if invalid_volumes:
            details["invalid_available_volumes"] = invalid_volumes
        if duplicate_target_solutes:
            details["duplicate_target_solutes"] = sorted(set(duplicate_target_solutes))
        if invalid_target_concentrations:
            details["invalid_target_concentrations"] = invalid_target_concentrations
        if invalid_target_volume:
            details["invalid_target_volume"] = target_volume_text
        return STATUS_CODES["INVALID_BLEND"], {
            "error": error_payload(
                "INVALID_BLEND",
                "母液与目标必须合法：母液标识不得重复，同一母液内溶质不得重复，"
                "目标溶质不得重复；浓度与体积必须是非负十进制字符串"
                "（不能为 JSON 数字），目标终体积必须为正数",
                details,
            ),
        }

    volume = decimal_to_fraction(final_volume)
    solution = solve_blend(stocks, targets, volume)
    if solution is None:
        # 逐溶质供应诊断：目标需求量（mmol）是否超出全部母液可用量合计
        insufficient = []
        for solute in sorted(targets):
            required = targets[solute] * volume
            supply = sum(
                (
                    stock.concentrations.get(solute, Fraction(0))
                    * stock.available_volume
                    for stock in stocks
                ),
                Fraction(0),
            )
            if required > supply:
                insufficient.append(
                    {
                        "solute": solute,
                        "required_mmol": fraction_to_decimal_string(required),
                        "available_mmol": fraction_to_decimal_string(supply),
                    }
                )
        details = {"target_volume": fraction_to_decimal_string(volume)}
        if insufficient:
            details["insufficient_solutes"] = insufficient
        return STATUS_CODES["NO_FEASIBLE_BLEND"], {
            "error": error_payload(
                "NO_FEASIBLE_BLEND",
                "无解：在可用体积与终体积约束下，不存在同时满足全部目标浓度的"
                "非负母液取用量（余量用水补足）；请降低目标浓度、增大终体积"
                "或补充母液",
                details,
            ),
        }

    # ---- 渲染：各母液取用量（输入顺序）、补水量与全部涉及溶质的最终浓度
    all_solutes = sorted(
        set(targets) | {s for stock in stocks for s in stock.concentrations}
    )
    final_concentrations: list[dict[str, Any]] = []
    for solute in all_solutes:
        final_value = sum(
            (
                stock.concentrations.get(solute, Fraction(0)) * used
                for stock, used in zip(stocks, solution.volumes)
            ),
            Fraction(0),
        ) / volume
        final_concentrations.append(
            {
                "solute": solute,
                "target_concentration": (
                    fraction_to_decimal_string(targets[solute])
                    if solute in targets
                    else None
                ),
                "final_concentration": fraction_to_decimal_string(final_value),
                "target_met": True if solute in targets else None,
            }
        )

    return 200, {
        "target_volume": fraction_to_decimal_string(volume),
        "total_stock_volume": fraction_to_decimal_string(sum(solution.volumes)),
        "water_added": fraction_to_decimal_string(solution.water),
        "stocks": [
            {
                "id": stock.id,
                "available_volume": fraction_to_decimal_string(
                    stock.available_volume
                ),
                "volume_used": fraction_to_decimal_string(used),
                "remaining_volume": fraction_to_decimal_string(
                    stock.available_volume - used
                ),
            }
            for stock, used in zip(stocks, solution.volumes)
        ],
        "final_concentrations": final_concentrations,
    }


# ------------------------------------------------------------ 同位素名义质量峰分布
def process_isotope_envelope(
    formula_text: str, element_items: list[dict[str, Any]]
) -> tuple[int, dict[str, Any]]:
    """同位素名义质量峰分布，返回 (http_status, payload)。

    先用既有解析器 :func:`parse_formula` 取得各元素原子数（电子物种
    拒绝，解析失败沿用 400 FORMULA_PARSE_ERROR）；随后聚合校验同位素
    输入：每个出现元素恰好一组、不得重复或给出不存在的元素，每组
    1～3 种互异正整数质量数，丰度均为非负十进制字符串且精确之和为 1，
    元素种类至多 3、总原子数至多 30——全部语义问题聚合为 422
    INVALID_ISOTOPE。校验通过后调用 :func:`app.isotopes.isotope_envelope`
    按独立原子抽样（多项分布）枚举全部名义质量与概率，同质量组合在
    枚举完成后合并，概率以最简分数字符串返回。该接口为无状态计算，
    不读写幂等缓存。
    """
    try:
        formula = parse_formula(formula_text)
    except FormulaError as exc:
        return STATUS_CODES["FORMULA_PARSE_ERROR"], {
            "error": error_payload(
                "FORMULA_PARSE_ERROR",
                str(exc),
                {"field": "formula", "input": formula_text},
            ),
        }

    if formula.is_electron:
        return STATUS_CODES["FORMULA_PARSE_ERROR"], {
            "error": error_payload(
                "FORMULA_PARSE_ERROR",
                "同位素质量峰分布不支持电子物种（e- / e+）：请提交含元素的化学式",
                {"field": "formula", "input": formula_text},
            ),
        }

    composition = formula.composition
    total_atoms = sum(composition.values())

    # ---- 逐元素校验同位素组
    specs: list[ElementIsotopes] = []
    seen_elements: set[str] = set()
    duplicate_elements: list[str] = []
    invalid_elements: list[dict[str, str]] = []
    unexpected_elements: list[dict[str, str]] = []
    invalid_abundances: list[dict[str, Any]] = []
    duplicate_mass_numbers: list[dict[str, Any]] = []
    abundance_sums: list[dict[str, str]] = []

    for item in element_items:
        element = item["element"]
        raw_isotopes = item["isotopes"]

        if element in seen_elements:
            duplicate_elements.append(element)
        seen_elements.add(element)

        if element not in ELEMENT_SYMBOLS:
            invalid_elements.append({"element": element})
            continue
        if element not in composition:
            unexpected_elements.append(
                {"element": element, "formula": formula.normalized}
            )
            continue

        seen_mass: set[int] = set()
        duplicates_here: list[int] = []
        isotopes: list[Isotope] = []
        abundance_total = Fraction(0)
        group_abundance_ok = True
        for raw in raw_isotopes:
            mass_number = raw["mass_number"]
            abundance_text = raw["abundance"]
            if mass_number in seen_mass and mass_number not in duplicates_here:
                duplicates_here.append(mass_number)
            seen_mass.add(mass_number)
            try:
                abundance = parse_moles(abundance_text)
            except (ValueError, InvalidOperation):
                invalid_abundances.append(
                    {
                        "element": element,
                        "mass_number": mass_number,
                        "abundance": abundance_text,
                    }
                )
                group_abundance_ok = False
                continue
            fraction_abundance = decimal_to_fraction(abundance)
            abundance_total += fraction_abundance
            isotopes.append(Isotope(mass_number, fraction_abundance))

        if duplicates_here:
            duplicate_mass_numbers.append(
                {"element": element, "mass_numbers": sorted(duplicates_here)}
            )
        # 仅当该组丰度全部可解析时才报告“和不为 1”，避免与单条非法丰度
        # 的明细重复；和以最简分数渲染（0.5 + 0.5 -> "1"）
        if group_abundance_ok and abundance_total != 1:
            abundance_sums.append(
                {"element": element, "sum": str(abundance_total)}
            )

        specs.append(
            ElementIsotopes(
                element=element,
                atom_count=composition[element],
                isotopes=tuple(isotopes),
            )
        )

    # ---- 缺失元素：化学式中出现但未提交同位素组的元素
    missing_elements = sorted(set(composition) - seen_elements)

    too_many_elements = len(composition) > MAX_ELEMENTS
    too_many_atoms = total_atoms > MAX_TOTAL_ATOMS

    if (
        invalid_elements
        or duplicate_elements
        or unexpected_elements
        or missing_elements
        or invalid_abundances
        or duplicate_mass_numbers
        or abundance_sums
        or too_many_elements
        or too_many_atoms
    ):
        details: dict[str, Any] = {}
        if invalid_elements:
            details["invalid_elements"] = invalid_elements
        if duplicate_elements:
            details["duplicate_elements"] = sorted(set(duplicate_elements))
        if missing_elements:
            details["missing_elements"] = missing_elements
        if unexpected_elements:
            details["unexpected_elements"] = unexpected_elements
        if invalid_abundances:
            details["invalid_abundances"] = invalid_abundances
        if duplicate_mass_numbers:
            details["duplicate_mass_numbers"] = duplicate_mass_numbers
        if abundance_sums:
            details["abundance_sums"] = abundance_sums
        if too_many_elements:
            details["element_count"] = len(composition)
            details["max_elements"] = MAX_ELEMENTS
        if too_many_atoms:
            details["total_atoms"] = total_atoms
            details["max_total_atoms"] = MAX_TOTAL_ATOMS
        return STATUS_CODES["INVALID_ISOTOPE"], {
            "error": error_payload(
                "INVALID_ISOTOPE",
                "同位素输入必须合法：化学式中的每种元素恰好给出一组同位素定义"
                "（不得缺失、重复或给出不存在的元素，元素种类至多 3 种、总原子"
                "数至多 30）；每组 1～3 种互异的正整数质量数，丰度均为非负"
                "十进制字符串（不能为 JSON 数字）且精确之和为 1",
                details,
            ),
        }

    # 与请求顺序一致（前面已保证元素恰好一一对应）
    ordered_specs = [
        next(spec for spec in specs if spec.element == item["element"])
        for item in element_items
    ]

    peaks = isotope_envelope(ordered_specs)
    top = most_probable_peak(peaks)

    def render_peak(peak) -> dict[str, Any]:
        return {"mass": peak.mass, "probability": str(peak.probability)}

    return 200, {
        "input": formula_text,
        "normalized": formula.normalized,
        "total_atoms": total_atoms,
        "elements": [
            {
                "element": item["element"],
                "atom_count": composition[item["element"]],
                "isotopes": [
                    {"mass_number": iso["mass_number"], "abundance": iso["abundance"]}
                    for iso in item["isotopes"]
                ],
            }
            for item in element_items
        ],
        "peaks": [render_peak(peak) for peak in peaks],
        "most_probable": render_peak(top),
    }


# ------------------------------------------------------------ 赫斯定律反应焓推导
def process_hess(
    target_item: dict[str, Any],
    known_items: list[dict[str, Any]],
) -> tuple[int, dict[str, Any]]:
    """赫斯定律反应焓推导，返回 (http_status, payload)。

    目标反应与每条已知反应都复用 :func:`process_balance` 的化学式解析与
    唯一配平（解析失败 400、同侧重复 409、不可唯一配平 422，错误 details
    中以 ``where``（``target`` / ``known_reactions``）与 ``index`` 定位）；
    各反应焓必须是可带符号的有限十进制字符串（kJ），非法时聚合为 422
    INVALID_REQUEST。随后调用 :func:`app.hess.solve_hess` 以精确有理数求
    已知反应计量向量的有符号线性组合（允许反向与分数倍），目标不在张成
    空间内返回 422 HESS_UNDERIVABLE，多组组合导出不同反应焓返回 422
    HESS_AMBIGUOUS；焓相同时正常返回各反应倍数、推导焓与逐物种净系数。
    配平复用 :func:`process_balance` 的共享幂等缓存（配平类错误因此可能
    带回 ``reaction_id``），但赫斯推导结果本身不单独缓存，成功响应不含
    ``reaction_id`` / ``cached``。
    """
    # ---- 目标反应：复用唯一配平（含解析、去重、幂等缓存）
    target_status, target_payload = process_balance(
        target_item["reactants"], target_item["products"]
    )
    if target_status != 200:
        error = target_payload.get("error")
        if error is not None:
            body = _hess_balance_error_body("target", None, error)
            rid = target_payload.get("reaction_id")
            if rid is not None:
                body["reaction_id"] = rid
            return target_status, body
        # 理论上不会走到：非 200 必为结构化错误
        return target_status, target_payload

    target_reactants = list(parse_side(target_item["reactants"], "反应物"))
    target_products = list(parse_side(target_item["products"], "生成物"))
    coef_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*target_reactants, *target_products),
            target_payload["coefficients"]["vector"],
        )
    }
    target_solved = _solved_from_parts(
        target_reactants, target_products, coef_map
    )

    # ---- 已知反应：逐条唯一配平 + 反应焓解析，聚合全部可定位的问题
    known_solved: list = []
    known_echos: list[dict[str, Any]] = []
    enthalpy_fractions: list[Fraction] = []
    invalid_enthalpies: list[dict[str, Any]] = []

    for index, item in enumerate(known_items):
        status, payload = process_balance(item["reactants"], item["products"])
        if status != 200:
            error = payload["error"]
            body = _hess_balance_error_body(
                "known_reactions", index, error
            )
            rid = payload.get("reaction_id")
            if rid is not None:
                body["reaction_id"] = rid
            return status, body

        reactants = list(parse_side(item["reactants"], "反应物"))
        products = list(parse_side(item["products"], "生成物"))
        cmap = {
            species_key(sp): coef
            for sp, coef in zip(
                (*reactants, *products),
                payload["coefficients"]["vector"],
            )
        }
        solved = _solved_from_parts(reactants, products, cmap)
        known_solved.append(solved)

        enthalpy_text = item["enthalpy_kj"]
        try:
            enthalpy = parse_enthalpy(enthalpy_text)
        except (ValueError, InvalidOperation):
            invalid_enthalpies.append(
                {"index": index, "enthalpy_kj": enthalpy_text}
            )
            enthalpy_fractions.append(Fraction(0))
        else:
            enthalpy_fractions.append(decimal_to_fraction(enthalpy))

        known_echos.append(
            {
                "index": index,
                "reactants": [sp.normalized for sp in reactants],
                "products": [sp.normalized for sp in products],
                "equation": payload["normalized"]["equation"],
                "enthalpy_kj": enthalpy_text,
            }
        )

    if invalid_enthalpies:
        return STATUS_CODES["INVALID_REQUEST"], {
            "error": error_payload(
                "INVALID_REQUEST",
                "已知反应焓必须是以 kJ 表示的有限十进制字符串"
                "（可带符号，如 \"-285.8\"，不能为 JSON 数字、NaN 或无穷大）",
                {"invalid_enthalpies": invalid_enthalpies},
            ),
        }

    known = [
        KnownReaction(solved=solved, enthalpy=h)
        for solved, h in zip(known_solved, enthalpy_fractions)
    ]

    try:
        solution = solve_hess(target_solved, known)
    except HessError as exc:
        details = dict(exc.details)
        # 规范键对客户端不友好：附带规范化物种名
        name_by_key = dict(normalized_by_key(target_solved))
        for kr in known:
            name_by_key.update(normalized_by_key(kr.solved))
        if "target_only_species" in details:
            details["target_only_species"] = [
                name_by_key.get(key, key) for key in details["target_only_species"]
            ]
        return STATUS_CODES[exc.code], {
            "error": error_payload(exc.code, exc.message, details),
        }

    # ---- 逐物种净系数：合并全部反应涉及的规范物种，按规范化名称排序
    name_by_key = dict(normalized_by_key(target_solved))
    for kr in known:
        name_by_key.update(normalized_by_key(kr.solved))

    ordered_keys = sorted(
        solution.species_order, key=lambda key: name_by_key.get(key, key)
    )
    net_coefficients = []
    for key in ordered_keys:
        value = solution.net_coefficients[key]
        if value > 0:
            side = "product"
        elif value < 0:
            side = "reactant"
        else:
            side = "cancelled"
        net_coefficients.append(
            {
                "normalized": name_by_key.get(key, key),
                "coefficient": str(value),
                "side": side,
            }
        )

    return 200, {
        "target": target_payload["normalized"],
        "known_reactions": [
            {
                **echo,
                "multiplier": str(multiplier),
                "multiplier_used": fraction_to_decimal_string(
                    multiplier * enthalpy_fractions[echo["index"]]
                ),
            }
            for echo, multiplier in zip(known_echos, solution.multipliers)
        ],
        "derived_enthalpy_kj": fraction_to_decimal_string(solution.enthalpy),
        "net_coefficients": net_coefficients,
    }


def _hess_balance_error_body(
    where: str, index: int | None, error: dict[str, Any]
) -> dict[str, Any]:
    """在既有结构化错误的 details 上补充 hess 定位信息（不修改原对象）。"""
    new_error = dict(error)
    details = dict(new_error.get("details") or {})
    details["where"] = where
    if index is not None:
        details["index"] = index
    new_error["details"] = details
    return {"error": new_error}


def _solved_from_parts(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int],
) -> SolvedReaction:
    """由已解析物种与配平系数映射构造 SolvedReaction（复用既有配平结果）。"""
    species = [*reactants, *products]
    coefficients = tuple(coefficient_map[species_key(sp)] for sp in species)
    elements = tuple(
        sorted({el for sp in species for el in sp.formula.composition})
    )
    return SolvedReaction(
        reactants=tuple(reactants),
        products=tuple(products),
        coefficients=coefficients,
        elements=elements,
    )


# ------------------------------------------------------------ 元素分析反推分子式
def process_formula_infer(
    element_items: list[dict[str, str]],
    mass_min: int,
    mass_max: int,
) -> tuple[int, dict[str, Any]]:
    """从 C/H/N/O 实测质量百分比反推分子式，返回 (http_status, payload)。

    逐条解析四种元素的测量值与绝对容差（均为 0～100 的非负有限十进制
    字符串，不要求合计为 100），并校验 C/H/N/O 各恰好一条、摩尔质量为
    正整数闭区间；问题聚合为 422 INVALID_FORMULA_INFER（缺元素/重复
    元素在 details.missing/duplicate 中给出）。校验通过后调用
    :func:`app.formula_infer.infer_formulas` 以精确有理数同时筛选落入
    全部百分比容差（闭区间）与质量区间的整数计数（逐项 0～30、总数
    大于零），按总百分比绝对偏差、摩尔质量、C/H/N/O 计数依次升序返回
    前 20 个。该接口为无状态计算，不读写幂等缓存。
    """
    # ---- 四种元素：一一对应（缺失/重复）+ 百分比与容差取值范围
    seen: list[str] = []
    duplicate_elements: list[str] = []
    invalid_percent: list[dict[str, Any]] = []
    invalid_tolerance: list[dict[str, Any]] = []
    measurements: dict[str, ElementMeasurement] = {}

    for item in element_items:
        element = item["element"]
        if element in seen:
            duplicate_elements.append(element)
        else:
            seen.append(element)

        percent_fraction = _parse_percent(item["percent"])
        if percent_fraction is None:
            invalid_percent.append(
                {"element": element, "percent": item["percent"]}
            )
        tolerance_fraction = _parse_percent(item["tolerance"])
        if tolerance_fraction is None:
            invalid_tolerance.append(
                {"element": element, "tolerance": item["tolerance"]}
            )
        if percent_fraction is not None and tolerance_fraction is not None:
            measurements[element] = ElementMeasurement(
                element=element,
                percent=percent_fraction,
                tolerance=tolerance_fraction,
            )

    missing_elements = sorted(set(FORMULA_ELEMENTS) - set(seen))

    # ---- 摩尔质量闭区间
    invalid_range = mass_min > mass_max

    if (
        duplicate_elements
        or missing_elements
        or invalid_percent
        or invalid_tolerance
        or invalid_range
    ):
        details: dict[str, Any] = {}
        if missing_elements:
            details["missing"] = missing_elements
        if duplicate_elements:
            details["duplicate"] = sorted(set(duplicate_elements))
        if invalid_percent:
            details["invalid_percent"] = invalid_percent
        if invalid_tolerance:
            details["invalid_tolerance"] = invalid_tolerance
        if invalid_range:
            details["invalid_molar_mass_range"] = {
                "minimum": mass_min,
                "maximum": mass_max,
            }
        return STATUS_CODES["INVALID_FORMULA_INFER"], {
            "error": error_payload(
                "INVALID_FORMULA_INFER",
                "分子式反推请求必须合法：C、H、N、O 四种元素各恰好一条测量"
                "记录（不能缺失或重复），实测百分比与绝对容差均必须是 0～100 "
                "范围内的非负十进制字符串（不能为 JSON 数字、负数、NaN 或"
                "无穷大），摩尔质量区间 minimum 不得大于 maximum",
                details,
            ),
        }

    candidates, total = infer_formulas(measurements, mass_min, mass_max)
    matches = []
    for candidate in candidates:
        matches.append(
            {
                "formula": render_formula(
                    candidate.c, candidate.h, candidate.n, candidate.o
                ),
                "molar_mass": candidate.molar_mass,
                "counts": {
                    "C": candidate.c,
                    "H": candidate.h,
                    "N": candidate.n,
                    "O": candidate.o,
                },
                "mass_percent": {
                    element: fraction_to_decimal_string(
                        candidate.mass_percent[element]
                    )
                    for element in FORMULA_ELEMENTS
                },
                "total_abs_deviation": fraction_to_decimal_string(
                    candidate.total_abs_deviation
                ),
            }
        )

    return 200, {
        "elements": element_items,
        "molar_mass": {"minimum": mass_min, "maximum": mass_max},
        "count": len(matches),
        "total": total,
        "truncated": total > len(matches),
        "matches": matches,
    }


def _parse_percent(text: str) -> Fraction | None:
    """解析 0～100 的非负有限十进制字符串为精确 Fraction；非法返回 None。"""
    try:
        value = parse_moles(text)
    except (ValueError, InvalidOperation):
        return None
    if value > 100:
        return None
    return decimal_to_fraction(value)
