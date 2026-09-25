"""FastAPI 应用入口。

路由（统一前缀 /api/v1）：

* ``POST /balance``          配平并返回最小正整数系数；
* ``POST /balance/variants`` 枚举零空间内全部互质正整数配平方案；
* ``POST /analyze``          逐元素差异说明（配平前 + 配平后视图）；
* ``POST /stoichiometry/batch`` 投料摩尔数 -> 反应进度、并列限量试剂、
  各反应物消耗/剩余与各生成物理论产量；
* ``POST /stoichiometry/plan``  生成物目标摩尔数 -> 最小反应进度、
  各反应物理论需求量、各目标产物达成量与超额量；
* ``POST /stoichiometry/yield`` 反应物投料摩尔数 + 各生成物实际摩尔产量
  -> 反应进度、逐种理论产量/实际产量/收率百分比（理论产量为 0 时收率
  为 null 并说明原因，实际高于理论时保留 >100%）；
* ``POST /stoichiometry/sequence`` 多步连续反应推演：按顺序执行多步
  反应，基于初始物种库存逐步计算最大反应进度、扣减消耗、加入产物，
  返回各步明细与最终库存；
* ``POST /equilibrium/solve`` 理想溶液平衡组成：单个可逆反应 + 各参与
  物种正初始摩尔数 + 体积（L）与 Kc，按配平系数确定双向非负物料区间，
  在区间内对 Qc(ξ)=Kc 二分求唯一内点反应进度（初始 Qc>Kc 时取负进度），
  返回进度、物料区间与各物种平衡摩尔数（绝对误差上界 1e-9 mol）；
* ``POST /inventory/allocate`` 试剂批次 FEFO 领用分配：按失效日期优先、
  同日批号升序跨批次扣减（名义量×纯度换算有效物质），跳过过期/隔离批次，
  返回每批领用/余量、未满足需求与汇总状态；
* ``POST /solutions/blend``   溶液配制：不超过 4 种母液（各溶质浓度 +
  可用体积）+ 目标终体积与各溶质目标浓度，以精确有理数求满足全部目标
  浓度的非负母液取用量（母液总量最小、再按输入顺序字典序最小），
  余量用水补足，返回各母液取用量、补水量与最终浓度；
* ``POST /isotopes/envelope`` 同位素名义质量峰分布：化学式 + 每种元素
  的同位素质量数与非负丰度（和为 1），按独立原子抽样（多项分布）枚举
  全部名义质量与概率，同质量组合合并，以最简分数返回按质量升序的峰
  列表与最高概率峰；
* ``POST /thermochemistry/hess`` 赫斯定律反应焓推导：目标反应 + 至多
  6 条已知反应（各含反应物、生成物、以 kJ 表示的反应焓），以精确有理
  数求已知反应计量向量的有符号（允许反向、分数倍）线性组合复现目标
  反应，返回各反应倍数、推导反应焓与合并抵消后的逐物种净系数；
* ``POST /formulas/infer`` 元素分析反推分子式：C/H/N/O 四种元素的实测
  质量百分比与绝对容差（0～100 十进制字符串，不要求合计 100）+ 摩尔
  质量正整数闭区间，以精确有理数枚举逐项 0～30、总数大于零且同时落入
  全部容差与质量区间的整数计数（C=12、H=1、N=14、O=16），按总百分比
  绝对偏差、摩尔质量、C/H/N/O 计数依次升序返回前 20 个分子式；
* ``GET  /reactions/{reaction_id}``  按幂等 ID 取回结果；
* ``GET  /health``   健康检查。
"""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from . import __version__
from .models import (
    AllocateRequest,
    AllocateResponse,
    AnalyzeResponse,
    BalanceResponse,
    BlendRequest,
    BlendResponse,
    EquilibriumRequest,
    EquilibriumResponse,
    FormulaInferRequest,
    FormulaInferResponse,
    HessRequest,
    HessResponse,
    IsotopeEnvelopeRequest,
    IsotopeEnvelopeResponse,
    ReactionRequest,
    ReversePlanRequest,
    ReversePlanResponse,
    SequenceRequest,
    SequenceResponse,
    StoichiometryResponse,
    StoichiometryRequest,
    VariantRequest,
    VariantsResponse,
    YieldRequest,
    YieldResponse,
)
from .service import (
    STATUS_CODES,
    get_reaction,
    process_allocate,
    process_analyze,
    process_balance,
    process_blend,
    process_equilibrium,
    process_formula_infer,
    process_hess,
    process_isotope_envelope,
    process_reverse_plan,
    process_sequence,
    process_stoichiometry,
    process_variants,
    process_yield,
)

app = FastAPI(
    title="化学反应计量平衡服务",
    description=(
        "解析含元素、括号与电荷的化学式，建立元素守恒与电荷守恒矩阵，"
        "通过有理数高斯-若尔当消元求解最小互质正整数系数；"
        "并提供限量试剂计量、连续反应推演与试剂批次 FEFO 领用分配。"
    ),
    version=__version__,
)


def _jsonable_errors(errors: list[dict]) -> list[dict]:
    """剔除 Pydantic v2 校验错误里不可 JSON 序列化的原始异常（ctx.error）。"""
    clean: list[dict] = []
    for err in errors:
        item = dict(err)
        ctx = item.get("ctx")
        if isinstance(ctx, dict) and "error" in ctx:
            item["ctx"] = {**ctx, "error": str(ctx["error"])}
        clean.append(item)
    return clean


@app.exception_handler(RequestValidationError)
async def handle_validation_error(request: Request, exc: RequestValidationError):
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "code": "INVALID_REQUEST",
                "message": "请求体校验失败：请检查必填字段、字段类型与取值范围"
                "（摩尔数/纯度/百分比须为十进制字符串、日期须为 YYYY-MM-DD、"
                "状态只能为 available 或 quarantined）",
                "details": {"validation": _jsonable_errors(exc.errors())},
            }
        },
    )


@app.get("/health", tags=["系统"])
async def health():
    return {"status": "ok", "version": __version__}


@app.post(
    "/api/v1/balance",
    response_model=BalanceResponse,
    tags=["配平"],
    summary="配平化学反应",
    responses={400: {"description": "化学式解析失败"}, 409: {"description": "同侧重复物种"},
               422: {"description": "无法配平 / 请求非法"}},
)
async def balance(request: ReactionRequest):
    status, payload = process_balance(request.reactants, request.products)
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/balance/variants",
    response_model=VariantsResponse,
    tags=["配平"],
    summary="枚举所有互质正整数配平方案（自由变量不超过 3 个）",
    responses={400: {"description": "化学式解析失败"}, 409: {"description": "同侧重复物种"},
               422: {"description": "无法配平 / 自由变量超过 3 个 / 请求非法"}},
)
async def balance_variants(request: VariantRequest):
    status, payload = process_variants(
        request.reactants,
        request.products,
        request.max_coefficient,
        request.limit,
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/analyze",
    response_model=AnalyzeResponse,
    tags=["配平"],
    summary="逐元素守恒差异分析",
    responses={400: {"description": "化学式解析失败"}},
)
async def analyze(request: ReactionRequest):
    status, payload = process_analyze(request.reactants, request.products)
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/batch",
    response_model=StoichiometryResponse,
    tags=["计量"],
    summary="批量限量试剂计量：反应进度、消耗/剩余与理论产量",
    responses={
        400: {"description": "化学式解析失败"},
        409: {"description": "同侧重复物种 / 投料重复"},
        422: {"description": "无法唯一配平 / 投料与反应物不对应 / 摩尔数非法"},
    },
)
async def stoichiometry_batch(request: StoichiometryRequest):
    status, payload = process_stoichiometry(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.amounts],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/plan",
    response_model=ReversePlanResponse,
    tags=["计量"],
    summary="反向投料规划：按生成物目标摩尔数反推最小反应进度与理论投料需求",
    responses={
        400: {"description": "化学式解析失败"},
        409: {"description": "同侧重复物种"},
        422: {"description": "无法唯一配平 / 目标与生成物不对应 / 摩尔数非法"},
    },
)
async def stoichiometry_plan(request: ReversePlanRequest):
    status, payload = process_reverse_plan(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.targets],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/yield",
    response_model=YieldResponse,
    tags=["计量"],
    summary="实际收率：投料摩尔数与各生成物实际产量 -> 理论产量、实际产量与收率百分比",
    responses={
        400: {"description": "化学式解析失败"},
        409: {"description": "同侧重复物种 / 投料或实际产量记录重复"},
        422: {"description": "无法唯一配平 / 投料或实际产量与反应物种不对应 / 摩尔数非法"},
    },
)
async def stoichiometry_yield(request: YieldRequest):
    status, payload = process_yield(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.amounts],
        [item.model_dump() for item in request.actual_yields],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/equilibrium/solve",
    response_model=EquilibriumResponse,
    tags=["平衡"],
    summary="理想溶液平衡组成：Qc=Kc 的唯一内点反应进度与各物种平衡摩尔数",
    responses={
        400: {"description": "化学式或初始量物种无法解析"},
        409: {"description": "同侧重复物种"},
        422: {"description": "无法唯一配平 / 初始量与物种不一一对应 / 两侧重复物种 / "
                            "初始摩尔数、体积或 Kc 非正十进制字符串"},
    },
)
async def equilibrium_solve(request: EquilibriumRequest):
    status, payload = process_equilibrium(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.amounts],
        request.volume_l,
        request.kc,
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/sequence",
    response_model=SequenceResponse,
    tags=["计量"],
    summary="连续反应推演：多步反应按顺序消耗/产出，逐步传递物种库存",
    responses={
        400: {"description": "化学式解析失败（反应式或初始库存物种）"},
        409: {"description": "某步同侧重复物种"},
        422: {"description": "某步无法唯一配平 / 步骤为空 / 库存摩尔数非法"},
    },
)
async def stoichiometry_sequence(request: SequenceRequest):
    status, payload = process_sequence(
        [item.model_dump() for item in request.steps],
        [item.model_dump() for item in request.initial_inventory],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/inventory/allocate",
    response_model=AllocateResponse,
    tags=["库存"],
    summary="试剂批次 FEFO 领用分配：按失效日期优先、同日批号升序跨批次扣减",
    responses={
        400: {"description": "领用物种或批次物种化学式无法解析"},
        422: {"description": "需求/批号重复、摩尔数或纯度非法、状态未知或请求体结构非法"},
    },
)
async def inventory_allocate(request: AllocateRequest):
    status, payload = process_allocate(
        [item.model_dump() for item in request.demands],
        [
            {**item.model_dump(), "expires_on": item.expires_on}
            for item in request.batches
        ],
        request.accounting_date,
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/solutions/blend",
    response_model=BlendResponse,
    tags=["溶液"],
    summary="溶液配制：满足全部目标浓度的最小母液取用量、补水量与最终浓度",
    responses={
        422: {"description": "输入非法（标识/溶质重复、浓度或体积非法、终体积非正）"
              " / 无可行解"},
    },
)
async def solutions_blend(request: BlendRequest):
    status, payload = process_blend(
        [item.model_dump() for item in request.stocks],
        request.target_volume,
        [item.model_dump() for item in request.targets],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/isotopes/envelope",
    response_model=IsotopeEnvelopeResponse,
    tags=["同位素"],
    summary="同位素名义质量峰分布：多项分布枚举全部名义质量与精确概率",
    responses={
        400: {"description": "化学式无法解析 / 电子物种"},
        422: {"description": "同位素输入非法（元素缺失/重复/多余、质量数重复"
              "或非正、丰度非法或和不为 1、超过规模上限）/ 请求体结构非法"},
    },
)
async def isotopes_envelope(request: IsotopeEnvelopeRequest):
    status, payload = process_isotope_envelope(
        request.formula,
        [item.model_dump() for item in request.elements],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/thermochemistry/hess",
    response_model=HessResponse,
    tags=["热化学"],
    summary="赫斯定律反应焓推导：已知反应的有符号线性组合 -> 目标反应焓",
    responses={
        400: {"description": "目标或某条已知反应的化学式无法解析"},
        409: {"description": "目标或某条已知反应同侧重复物种"},
        422: {"description": "某反应无法唯一配平 / 反应焓非法 / "
              "目标不在张成空间（HESS_UNDERIVABLE）/ 组合歧义（HESS_AMBIGUOUS）"},
    },
)
async def thermochemistry_hess(request: HessRequest):
    status, payload = process_hess(
        request.target.model_dump(),
        [item.model_dump() for item in request.known_reactions],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/formulas/infer",
    response_model=FormulaInferResponse,
    tags=["元素分析"],
    summary="由 C/H/N/O 实测质量百分比与摩尔质量区间反推分子式",
    responses={
        422: {"description": "元素缺失/重复、百分比或容差超出 0～100、"
              "摩尔质量区间倒置或请求体结构非法"},
    },
)
async def formulas_infer(request: FormulaInferRequest):
    status, payload = process_formula_infer(
        [item.model_dump() for item in request.elements],
        request.molar_mass.minimum,
        request.molar_mass.maximum,
    )
    return JSONResponse(status_code=status, content=payload)


@app.get(
    "/api/v1/reactions/{reaction_id}",
    tags=["配平"],
    summary="按反应 ID 查询缓存结果（幂等）",
    responses={404: {"description": "未找到该反应"}},
)
async def get_one(reaction_id: str):
    payload = get_reaction(reaction_id)
    if payload is None:
        return JSONResponse(
            status_code=404,
            content={
                "error": {
                    "code": "NOT_FOUND",
                    "message": f"未找到反应 {reaction_id}（缓存可能已淘汰）",
                    "details": {"reaction_id": reaction_id},
                }
            },
        )
    status = 200 if "error" not in payload else STATUS_CODES[payload["error"]["code"]]
    return JSONResponse(status_code=status, content=payload)
