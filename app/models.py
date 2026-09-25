"""Pydantic 请求/响应模型（同时驱动 OpenAPI 文档）。"""

from __future__ import annotations

from datetime import date
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator


class ReactionRequest(BaseModel):
    reactants: list[str] = Field(
        ...,
        min_length=1,
        description="反应物化学式，如 [\"KMnO4\", \"HCl\"]；电荷写作 SO4^2-、Fe3+、e-",
        examples=[["H2", "O2"]],
    )
    products: list[str] = Field(
        ...,
        min_length=1,
        description="生成物化学式，如 [\"KCl\", \"MnCl2\", \"Cl2\", \"H2O\"]",
        examples=[["H2O"]],
    )

    @field_validator("reactants", "products")
    @classmethod
    def _strip(cls, values: list[str]) -> list[str]:
        cleaned: list[str] = []
        for value in values:
            if not isinstance(value, str):
                raise ValueError("化学式必须是字符串")
            stripped = value.strip()
            if not stripped:
                raise ValueError("化学式不能为空字符串")
            cleaned.append(stripped)
        return cleaned


class SpeciesCoefficient(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数系数")


class CoefficientBlock(BaseModel):
    reactants: list[SpeciesCoefficient]
    products: list[SpeciesCoefficient]
    vector: list[int] = Field(description="按 reactants + products 顺序的系数向量")


class NormalizedReaction(BaseModel):
    reactants: list[str]
    products: list[str]
    equation: str = Field(description="规范化的化学方程式文本")


class BalanceResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="本次结果是否命中幂等缓存")
    normalized: NormalizedReaction
    coefficients: CoefficientBlock
    balanced: bool
    elements: list[str]
    charge_balance: dict[str, Any]


class ElementRow(BaseModel):
    element: str
    reactant_total: int
    product_total: int
    difference: int = Field(description="反应物合计 - 生成物合计；0 表示守恒")
    balanced: bool


class ChargeRow(BaseModel):
    reactant_total: int
    product_total: int
    difference: int
    balanced: bool


class SideSpecies(BaseModel):
    input: str
    normalized: str
    coefficient: int


class StoichView(BaseModel):
    coefficients_assumed: int | None = Field(
        description="该视图使用的统一假设系数；平衡视图省略此字段"
    )
    species: dict[str, list[SideSpecies]]
    elements: list[ElementRow]
    charge: ChargeRow


class AnalyzeResponse(BaseModel):
    reaction_id: str
    cached: bool
    normalized: NormalizedReaction
    raw: StoichView = Field(description="按单位系数统计的原始差异（配平前）")
    balanced: StoichView | None = Field(
        default=None, description="配平后的守恒视图；无法配平时为 null"
    )
    balanced_ok: bool
    error: dict[str, Any] | None = None


class ErrorResponse(BaseModel):
    error: dict[str, Any]


class VariantRequest(ReactionRequest):
    """``/balance/variants`` 请求：在 ReactionRequest 之上加枚举参数。"""

    max_coefficient: int = Field(
        ...,
        ge=2,
        le=30,
        description="允许的最大系数（含端点，范围 2～30）",
    )
    limit: int = Field(
        ...,
        ge=1,
        le=100,
        description="最多返回的配平方案数（含端点，范围 1～100）",
    )


class VariantItem(BaseModel):
    coefficients: CoefficientBlock = Field(
        description="按 reactants + products 顺序的一个互质正整数解"
    )
    sum: int = Field(description="该方案的系数总和（主排序键）")
    equation: str = Field(description="该方案对应的规范化化学方程式")


class VariantsResponse(BaseModel):
    reaction_id: str
    normalized: NormalizedReaction
    max_coefficient: int
    limit: int
    nullity: int = Field(description="守恒矩阵零空间维数（1～3）")
    count: int = Field(description="本次实际返回的方案数")
    total: int = Field(description="枚举到的全部方案数（截取前）")
    truncated: bool = Field(description="方案总数超过 limit、已截取时为 true")
    variants: list[VariantItem]


# ------------------------------------------------------------ 限量试剂计量
class StoichAmount(BaseModel):
    species: str = Field(
        ...,
        min_length=1,
        description="反应物化学式；必须与规范化后的反应物一一对应（等价写法可识别）",
        examples=["H2"],
    )
    moles: str = Field(
        ...,
        min_length=1,
        description="投料摩尔数的非负十进制字符串（如 \"1.5\"、\"0\"、\"1e-3\"），"
        "不要传 JSON 数字以避免浮点误差",
        examples=["2.0"],
    )

    @field_validator("species", "moles")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class StoichiometryRequest(ReactionRequest):
    """``/stoichiometry/batch`` 请求：反应式 + 每种反应物的投料摩尔数。"""

    amounts: list[StoichAmount] = Field(
        ...,
        min_length=1,
        description="每种反应物一条投料记录；数量必须恰好等于规范化反应物数量，"
        "且 species 一一对应（无缺失、重复、未知物种）",
    )


class ReactionExtent(BaseModel):
    value: str = Field(description="反应进度 ξ = min(投料量 ÷ 系数)，十进制字符串")
    unit: str = Field(default="mol", description="反应进度以 mol 计（按本反应配平式）")


class StoichReactantResult(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    fed: str = Field(description="投料摩尔数（回显）")
    consumed: str = Field(description="按反应进度计算的消耗量（mol）")
    remaining: str = Field(description="剩余量 = 投料 - 消耗（mol）")
    is_limiting: bool = Field(description="是否为并列限量试剂之一")


class StoichProductResult(BaseModel):
    input: str
    normalized: str
    coefficient: int
    theoretical_yield: str = Field(description="理论产量（mol）")


class StoichiometryResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="配平结果是否命中幂等缓存")
    normalized: NormalizedReaction
    reaction_extent: ReactionExtent
    limiting_reactants: list[str] = Field(
        description="全部并列限量试剂的规范化化学式（投料比恰好等于系数比时多于一个）"
    )
    reactants: list[StoichReactantResult]
    products: list[StoichProductResult]


# ------------------------------------------------------------ 理想溶液平衡组成
class EquilibriumInitialAmount(StoichAmount):
    """平衡计算的一条初始摩尔数记录；摩尔数为**正**十进制字符串。"""

    species: str = Field(
        ...,
        min_length=1,
        description="参与物种化学式；必须与规范化后的反应物+生成物一一对应"
        "（等价写法可识别），两侧均须给出",
        examples=["H2"],
    )
    moles: str = Field(
        ...,
        min_length=1,
        description="初始摩尔数的**正**十进制字符串（如 \"1.5\"、\"1e-3\"，"
        "不接受 \"0\" 或负数），不要传 JSON 数字以避免浮点误差",
        examples=["2.0"],
    )


class EquilibriumRequest(ReactionRequest):
    """``/equilibrium/solve`` 请求：可逆反应式 + 全部参与物种初始摩尔数 + 体积 + Kc。"""

    amounts: list[EquilibriumInitialAmount] = Field(
        ...,
        min_length=1,
        description="每个参与物种（反应物 ∪ 生成物）一条初始摩尔数记录；"
        "数量必须恰好等于规范化物种总数，且 species 一一对应"
        "（无缺失、重复、未知物种，同一物种不得同时出现在两侧），"
        "摩尔数必须严格为正",
    )
    volume_l: str = Field(
        ...,
        min_length=1,
        description="溶液体积，单位 L，**正**十进制字符串（如 \"2.0\"、\"0.5\"），"
        "不要传 JSON 数字",
        examples=["1.0"],
    )
    kc: str = Field(
        ...,
        min_length=1,
        description="浓度平衡常数 Kc（无量纲或以 mol/L 约定），**正**十进制"
        "字符串（如 \"4.0\"、\"1e-3\"），不要传 JSON 数字",
        examples=["1.0"],
    )

    @field_validator("volume_l", "kc")
    @classmethod
    def _strip_positive_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class EquilibriumSpeciesResult(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    initial_moles: str = Field(description="客户端提交的初始摩尔数（回显）")
    equilibrium_moles: str = Field(
        description="平衡摩尔数 = 初始摩尔数 − 生成物计量变化 / + 反应物计量"
        "变化（mol，非负）"
    )


class EquilibriumExtentBounds(BaseModel):
    minimum: str = Field(
        description="非负物料区间下界 ξ_min = −min_生成物(初始量 ÷ 系数) ≤ 0（mol）"
    )
    maximum: str = Field(
        description="非负物料区间上界 ξ_max = min_反应物(初始量 ÷ 系数) ≥ 0（mol）"
    )


class EquilibriumResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="配平结果是否命中幂等缓存")
    normalized: NormalizedReaction
    volume_l: str = Field(description="溶液体积（L，回显）")
    kc: str = Field(description="浓度平衡常数（回显）")
    direction: Literal["forward", "reverse", "at_equilibrium"] = Field(
        description="平衡移动方向：forward 正向（ξ>0）/ reverse 逆向（ξ<0，"
        "初始 Qc>Kc）/ at_equilibrium 初始即平衡（ξ=0）"
    )
    reaction_extent: ReactionExtent
    extent_bounds: EquilibriumExtentBounds = Field(
        description="由配平系数与初始量确定的非负物料区间 [ξ_min, ξ_max]"
    )
    initial_qc: str = Field(
        description="按初始浓度计量次方构造的反应商 Qc(0)（十进制字符串）"
    )
    equilibrium_qc: str = Field(
        description="按返回平衡组成构造的反应商，应与 Kc 一致（十进制字符串）"
    )
    reactants: list[EquilibriumSpeciesResult]
    products: list[EquilibriumSpeciesResult]


# ------------------------------------------------------------ 实际收率
class ActualYield(StoichAmount):
    """一条生成物实际摩尔产量记录；复用非负十进制字符串摩尔数校验。"""

    species: str = Field(
        ...,
        min_length=1,
        description="生成物化学式；必须与规范化后的生成物一一对应（等价写法可识别）",
        examples=["H2O"],
    )
    moles: str = Field(
        ...,
        min_length=1,
        description="该生成物的**实际**摩尔产量（非负十进制字符串，如 \"1.5\"、"
        "\"0\"、\"1e-3\"），允许高于理论值；不要传 JSON 数字以避免浮点误差",
        examples=["1.8"],
    )


class YieldRequest(ReactionRequest):
    """``/stoichiometry/yield`` 请求：反应式 + 各反应物投料 + 各生成物实际产量。"""

    amounts: list[StoichAmount] = Field(
        ...,
        min_length=1,
        description="每种反应物一条投料记录；数量必须恰好等于规范化反应物数量，"
        "且 species 一一对应（无缺失、重复、未知物种）",
    )
    actual_yields: list[ActualYield] = Field(
        ...,
        min_length=1,
        description="每种生成物一条实际产量记录；数量必须恰好等于规范化生成物"
        "数量，且 species 一一对应（无缺失、重复、未知物种）；实际产量允许为 0 "
        "或高于理论值",
    )


class YieldProductResult(BaseModel):
    input: str = Field(description="客户端提交的原始生成物化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    theoretical_yield: str = Field(
        description="理论产量 = 反应进度 × 系数（mol）"
    )
    actual_yield: str = Field(description="客户端提交的实际产量（mol，回显）")
    yield_percent: str | None = Field(
        description="收率百分比 = 实际产量 ÷ 理论产量 × 100（十进制字符串，"
        "不封顶，实际高于理论时大于 100）；理论产量为 0 时为 null"
    )
    yield_note: str | None = Field(
        default=None,
        description="收率无法计算时的说明（理论产量为 0 时给出零投料限量反应物）；"
        "正常计算时为 null",
    )


class YieldResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="配平结果是否命中幂等缓存")
    normalized: NormalizedReaction
    reaction_extent: ReactionExtent
    limiting_reactants: list[str] = Field(
        description="全部并列限量试剂的规范化化学式（投料比恰好等于系数比时多于一个）"
    )
    reactants: list[StoichReactantResult]
    products: list[YieldProductResult] = Field(
        description="逐种生成物的理论产量、实际产量与收率百分比"
    )


# ------------------------------------------------------------ 反向投料规划
class ReversePlanTarget(StoichAmount):
    """反向规划的一条目标产物记录；复用非负十进制字符串摩尔数校验。"""

    species: str = Field(
        ...,
        min_length=1,
        description="目标生成物化学式；必须是规范化生成物之一（等价写法可识别），可只提交部分生成物",
        examples=["H2O"],
    )


class ReversePlanRequest(ReactionRequest):
    """``/stoichiometry/plan`` 请求：反应式 + 一个或多个生成物的目标摩尔数。"""

    targets: list[ReversePlanTarget] = Field(
        ...,
        min_length=1,
        description="目标生成物记录，至少一条；species 必须是规范化生成物"
        "（无重复、未知物种），不必覆盖全部生成物",
    )


class PlanReactantResult(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    required: str = Field(
        description="满足全部目标所需的理论投料量 = 反应进度 × 系数（mol）"
    )


class PlanTargetResult(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    target: str = Field(description="客户端请求的目标摩尔数（回显）")
    achieved: str = Field(
        description="按最小反应进度该产物的实际达成量 = 反应进度 × 系数（mol）"
    )
    excess: str = Field(description="超额量 = 达成量 − 目标量（mol），驱动产物为 0")
    is_driving: bool = Field(
        description="是否为决定最小反应进度的目标产物之一"
        "（目标量 ÷ 系数取得最大值，可并列多个）"
    )


class PlanProductResult(BaseModel):
    input: str
    normalized: str
    coefficient: int
    achieved: str = Field(
        description="按最小反应进度该生成物的实际达成量 = 反应进度 × 系数（mol）"
    )


class ReversePlanResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="配平结果是否命中幂等缓存")
    normalized: NormalizedReaction
    reaction_extent: ReactionExtent
    driving_targets: list[str] = Field(
        description="全部并列决定最小反应进度的目标产物规范化化学式"
        "（目标比恰好等于系数比时多于一个）"
    )
    reactants: list[PlanReactantResult] = Field(description="每种反应物的理论需求量")
    targets: list[PlanTargetResult] = Field(
        description="客户端提交的各目标产物：目标量、达成量与超额量"
    )
    products: list[PlanProductResult] = Field(description="全部生成物的达成量")


# ------------------------------------------------------------ 连续反应推演
class SequenceStepRequest(ReactionRequest):
    """连续推演中的一步反应；本身就是一条完整的反应物/生成物反应式。"""


class InventoryEntry(StoichAmount):
    """连续推演的一条初始库存记录；摩尔数为非负十进制字符串。"""

    species: str = Field(
        ...,
        min_length=1,
        description="库存物种化学式；等价写法（括号风格、空白、电荷写法不同）"
        "会合并为同一库存物种，摩尔数相加",
        examples=["H2"],
    )
    moles: str = Field(
        ...,
        min_length=1,
        description="初始库存摩尔数的非负十进制字符串（如 \"1.5\"、\"0\"、\"1e-3\"），"
        "不要传 JSON 数字以避免浮点误差",
        examples=["2.0"],
    )


class SequenceRequest(BaseModel):
    """``/stoichiometry/sequence`` 请求：按顺序执行的多步反应与初始物种库存。"""

    steps: list[SequenceStepRequest] = Field(
        ...,
        min_length=1,
        description="按顺序执行的反应步骤；每步都必须能解析并得到唯一配平",
    )
    initial_inventory: list[InventoryEntry] = Field(
        ...,
        description="初始物种摩尔库存；等价化学式合并为同一物种（摩尔数相加），"
        "未列出的物种库存视为 0",
    )


class SequenceInventoryItem(BaseModel):
    normalized: str = Field(description="规范化化学式")
    moles: str = Field(description="库存摩尔数（非负十进制字符串）")


class SequenceReactantResult(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数配平系数")
    before: str = Field(description="该步开始前的库存摩尔数")
    consumed: str = Field(description="本步消耗量 = 反应进度 × 系数（mol）")
    after: str = Field(description="该步结束后的库存摩尔数")
    is_blocking: bool = Field(
        description="是否为本步的阻塞反应物之一（零进度时库存为 0 的全部反应物）"
    )


class SequenceProductResult(BaseModel):
    input: str
    normalized: str
    coefficient: int
    before: str = Field(description="该步开始前的库存摩尔数（含上一步产物）")
    produced: str = Field(description="本步产出量 = 反应进度 × 系数（mol）")
    after: str = Field(description="该步结束后的库存摩尔数")


class SequenceStepResult(BaseModel):
    index: int = Field(description="从 0 开始的步骤序号（与请求 steps 顺序一致）")
    reaction_id: str = Field(description="该步反应的幂等 ID（与 /balance 共享签名）")
    cached: bool = Field(description="该步配平结果是否命中幂等缓存")
    normalized: NormalizedReaction
    reaction_extent: ReactionExtent
    blocked: bool = Field(description="该步是否因库存不足以零进度继续")
    blocking_reactants: list[str] = Field(
        description="零进度时库存为 0 的全部反应物规范化化学式；正常进度时为空"
    )
    reactants: list[SequenceReactantResult] = Field(
        description="各反应物该步的步前库存、消耗与步后库存"
    )
    products: list[SequenceProductResult] = Field(
        description="各生成物该步的步前库存、产出与步后库存"
    )
    inventory: list[SequenceInventoryItem] = Field(
        description="该步结束后的全部物种库存快照（最后一步即 final_inventory）"
    )


class SequenceResponse(BaseModel):
    steps: list[SequenceStepResult] = Field(description="按顺序的各步推演结果")
    final_inventory: list[SequenceInventoryItem] = Field(
        description="全部步骤执行完毕后的最终物种库存"
    )


# ------------------------------------------------------------ 试剂批次领用分配
class AllocationDemand(BaseModel):
    """一条待领用物种及其有效摩尔需求。"""

    species: str = Field(
        ...,
        min_length=1,
        description="待领用物种化学式；等价写法（括号风格、空白、电荷写法不同）"
        "会识别为同一物种",
        examples=["HCl"],
    )
    effective_moles: str = Field(
        ...,
        min_length=1,
        description="需要领用的**有效摩尔数**十进制字符串（如 \"1.5\"、\"0\"、"
        "\"1e-3\"），不要传 JSON 数字以避免浮点误差；服务以 名义量 × 纯度 换算有效物质",
        examples=["2.0"],
    )

    @field_validator("species", "effective_moles")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class InventoryBatch(BaseModel):
    """一条带批号、纯度、可用名义摩尔数、失效日期与状态的库存批次。"""

    batch_no: str = Field(
        ...,
        min_length=1,
        description="库存批号；在单次请求内必须全局唯一，同一批号不得重复提交"
        "（同一批次不得重复扣减）",
        examples=["B-2025-001"],
    )
    species: str = Field(
        ...,
        min_length=1,
        description="批次内试剂化学式；等价写法（括号风格、空白、电荷写法不同）"
        "会识别为同一物种",
        examples=["HCl"],
    )
    purity: str = Field(
        ...,
        min_length=1,
        description="纯度（质量/摩尔分数）十进制字符串，必须在 (0, 1] 区间内"
        "（如 \"0.98\"、\"1.0\"）；有效摩尔 = 名义摩尔 × 纯度",
        examples=["0.98"],
    )
    nominal_moles: str = Field(
        ...,
        min_length=1,
        description="该批可用**名义摩尔数**的非负十进制字符串"
        "（如 \"5.0\"、\"0\"、\"1e-3\"），不要传 JSON 数字",
        examples=["5.0"],
    )
    expires_on: date = Field(
        ...,
        description="批次失效日期，ISO 8601（YYYY-MM-DD）；早于 accounting_date"
        "即视为已过期并跳过",
        examples=["2026-12-31"],
    )
    status: Literal["available", "quarantined"] = Field(
        ...,
        description="批次状态：available 可领用；quarantined 隔离中，跳过且不扣减",
        examples=["available"],
    )

    @field_validator("batch_no", "species", "purity", "nominal_moles")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class AllocateRequest(BaseModel):
    """``POST /inventory/allocate`` 请求：领用需求、库存批次与核算日期。"""

    demands: list[AllocationDemand] = Field(
        ...,
        min_length=1,
        description="待领用物种及其有效摩尔需求，至少一条；同一物种不得重复",
    )
    batches: list[InventoryBatch] = Field(
        ...,
        description="库存批次；批号必须全局唯一。按失效日期升序、同日批号升序"
        "（FEFO）分配，过期或隔离批次跳过，需求可跨批次拆分",
    )
    accounting_date: date = Field(
        ...,
        description="核算日期，ISO 8601（YYYY-MM-DD）；批次 expires_on 早于此日"
        "即视为已过期（恰为同一天不算过期）",
        examples=["2026-09-23"],
    )


class AllocationBatchLine(BaseModel):
    """一个库存批次在某次分配中的扣减与余量明细。"""

    batch_no: str = Field(description="库存批号（回显）")
    species: str = Field(description="客户端提交的批次物种原始化学式")
    normalized: str = Field(description="规范化后的物种化学式")
    expires_on: str = Field(description="批次失效日期（YYYY-MM-DD）")
    status: Literal["available", "quarantined"] = Field(description="批次状态（回显）")
    purity: str = Field(description="纯度（十进制字符串）")
    nominal_moles: str = Field(description="可用名义摩尔数（回显）")
    effective_moles: str = Field(
        description="该批可提供的有效摩尔 = 名义摩尔 × 纯度"
    )
    eligible: bool = Field(
        description="核算日是否具备领用资格（未过期、未隔离且有效量大于 0）"
    )
    skipped: bool = Field(description="该批是否被跳过（未参与扣减）")
    skip_reason: Literal["expired", "quarantined", "exhausted"] | None = Field(
        description="跳过原因：expired 已过期 / quarantined 隔离 / exhausted "
        "可用名义量为 0；实际领用的批次为 null"
    )
    used: bool = Field(description="该批是否实际承担了正的领用扣减")
    allocated_nominal_moles: str = Field(
        description="从该批领用的**名义**摩尔数 = 有效领用 ÷ 纯度；未领用为 \"0\""
    )
    allocated_effective_moles: str = Field(
        description="从该批领用的**有效**摩尔数；未领用为 \"0\""
    )
    remaining_nominal_moles: str = Field(
        description="领用后该批剩余名义摩尔数 = 名义量 − 名义领用"
    )
    remaining_effective_moles: str = Field(
        description="领用后该批剩余有效摩尔 = 剩余名义量 × 纯度"
    )


class SpeciesAllocation(BaseModel):
    """单个待领用物种的跨批次分配明细。"""

    species: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    requested_effective_moles: str = Field(description="请求的有效摩尔需求（回显）")
    allocated_effective_moles: str = Field(
        description="从全部批次实际领到的有效摩尔合计"
    )
    unmet_effective_moles: str = Field(
        description="未满足的有效摩尔需求 = 需求 − 已分配（为 \"0\" 即全部满足）"
    )
    status: Literal["fulfilled", "shortage"] = Field(
        description="fulfilled 需求全部满足；shortage 仍有未满足需求"
    )
    batches: list[AllocationBatchLine] = Field(
        description="该物种的全部提交批次，按失效日期升序、同日批号升序排列，"
        "逐批给出领用与余量（含被跳过的批次）"
    )


class UnmatchedBatch(AllocationBatchLine):
    """物种未出现在领用需求中的批次；``skip_reason`` 恒为 unmatched_species。"""

    skip_reason: Literal["unmatched_species"] = Field(
        description="批次物种未出现在 demands 中，整批保留、未参与分配"
    )


class AllocationSummary(BaseModel):
    total_requested_effective_moles: str = Field(
        description="全部物种的有效摩尔需求合计"
    )
    total_allocated_effective_moles: str = Field(
        description="全部物种实际领到的有效摩尔合计"
    )
    total_unmet_effective_moles: str = Field(description="未满足需求合计")
    demands_total: int = Field(description="领用需求条数")
    fulfilled_demands: int = Field(description="完全满足的需求条数")
    shortage_demands: int = Field(description="未完全满足的需求条数")
    batches_total: int = Field(description="请求提交的批次总数")
    batches_used: int = Field(description="实际承担领用扣减的批次数")
    batches_skipped: int = Field(
        description="被跳过的批次数（过期、隔离或名义量为 0；不含物种未匹配批次）"
    )
    unmatched_batches: int = Field(description="物种未出现在 demands 中的批次数")
    status: Literal["fulfilled", "shortage"] = Field(
        description="汇总状态：全部需求满足为 fulfilled，否则为 shortage"
    )


class AllocateResponse(BaseModel):
    as_of: str = Field(description="核算日期（YYYY-MM-DD，回显）")
    allocations: list[SpeciesAllocation] = Field(
        description="按规范化化学式排序的各物种分配明细；物种下批次按 FEFO 排序"
    )
    unmatched_batches: list[UnmatchedBatch] = Field(
        description="物种未出现在 demands 中的批次（整批保留、未参与分配），"
        "按规范化物种、失效日期、批号排序"
    )
    summary: AllocationSummary = Field(description="跨物种汇总与总体状态")


# ------------------------------------------------------------ 溶液配制（母液混合）
class BlendSolute(BaseModel):
    """母液中的一条溶质浓度记录。"""

    solute: str = Field(
        ...,
        min_length=1,
        description="溶质标识（任意非空字符串，如 \"NaCl\"）；同一母液内不得重复",
        examples=["NaCl"],
    )
    concentration: str = Field(
        ...,
        min_length=1,
        description="该溶质在母液中的浓度，非负十进制字符串，单位 mol/L"
        "（如 \"2.0\"、\"0\"、\"1e-3\"），不要传 JSON 数字以避免浮点误差",
        examples=["2.0"],
    )

    @field_validator("solute", "concentration")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class BlendStock(BaseModel):
    """一条母液：标识、各溶质浓度与可用体积。"""

    id: str = Field(
        ...,
        min_length=1,
        description="母液标识；单次请求内必须唯一",
        examples=["stock-A"],
    )
    solutes: list[BlendSolute] = Field(
        ...,
        description="母液中各溶质的浓度记录；同一溶质不得重复；空列表表示纯溶剂",
    )
    available_volume: str = Field(
        ...,
        min_length=1,
        description="该母液的可用体积，非负十进制字符串，单位 mL"
        "（如 \"100\"、\"0\"），不要传 JSON 数字",
        examples=["100"],
    )

    @field_validator("id", "available_volume")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class BlendTarget(BaseModel):
    """一条目标溶质浓度记录。"""

    solute: str = Field(
        ...,
        min_length=1,
        description="目标溶质标识；在 targets 内不得重复",
        examples=["NaCl"],
    )
    concentration: str = Field(
        ...,
        min_length=1,
        description="该溶质的目标终浓度，非负十进制字符串，单位 mol/L"
        "（如 \"0.1\"、\"0\"），不要传 JSON 数字",
        examples=["0.1"],
    )

    @field_validator("solute", "concentration")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class BlendRequest(BaseModel):
    """``POST /solutions/blend`` 请求：母液清单、目标终体积与各溶质目标浓度。"""

    stocks: list[BlendStock] = Field(
        ...,
        min_length=1,
        max_length=4,
        description="参与配制的母液，1～4 种；标识不得重复",
    )
    target_volume: str = Field(
        ...,
        min_length=1,
        description="目标终体积，**正**十进制字符串，单位 mL（如 \"250\"），"
        "不要传 JSON 数字",
        examples=["250"],
    )
    targets: list[BlendTarget] = Field(
        ...,
        min_length=1,
        description="各溶质目标终浓度，至少一条；同一溶质不得重复；"
        "母液带入但未列入目标的溶质不约束，仅在最终浓度中报告",
    )


class BlendStockResult(BaseModel):
    """一条母液的取用结果。"""

    id: str = Field(description="母液标识（回显，与请求顺序一致）")
    available_volume: str = Field(description="可用体积（mL，回显）")
    volume_used: str = Field(description="本次取用量（mL）")
    remaining_volume: str = Field(description="取用后剩余体积 = 可用体积 − 取用量（mL）")


class BlendSoluteResult(BaseModel):
    """一个溶质的最终浓度核对结果。"""

    solute: str = Field(description="溶质标识")
    target_concentration: str | None = Field(
        description="目标终浓度（mol/L，回显）；未列入目标的溶质为 null"
    )
    final_concentration: str = Field(
        description="配制后的最终浓度 = Σ(母液浓度 × 取用量) ÷ 终体积（mol/L）"
    )
    target_met: bool | None = Field(
        description="是否达到目标浓度；目标溶质恒为 true，非目标溶质为 null"
    )


class BlendResponse(BaseModel):
    target_volume: str = Field(description="目标终体积（mL，回显）")
    total_stock_volume: str = Field(description="母液总取用量（mL，已最小化）")
    water_added: str = Field(
        description="补水量 = 终体积 − 母液总取用量（mL，非负）"
    )
    stocks: list[BlendStockResult] = Field(
        description="各母液取用量，与请求 stocks 顺序一致"
    )
    final_concentrations: list[BlendSoluteResult] = Field(
        description="全部涉及溶质（目标 ∪ 母液带入）的最终浓度，按溶质标识排序"
    )


# ------------------------------------------------------------ 赫斯定律反应焓推导
class HessKnownReaction(ReactionRequest):
    """一条已知反应：反应物/生成物 + 以 kJ 表示的反应焓十进制字符串。"""

    enthalpy_kj: str = Field(
        ...,
        min_length=1,
        description="该已知反应的反应焓（ΔH），单位 kJ，**可带符号的十进制"
        "字符串**（如 \"-285.8\"、\"0\"、\"1e2\"）；不要传 JSON 数字以避免浮点误差。"
        "符号约定：放热为负、吸热为正，与反应书写方向一致（反向使用时自动变号）",
        examples=["-285.8"],
    )

    @field_validator("enthalpy_kj")
    @classmethod
    def _strip_enthalpy(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("反应焓必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("反应焓不能为空白字符串")
        return stripped


class HessRequest(BaseModel):
    """``POST /thermochemistry/hess`` 请求：目标反应 + 至多 6 条已知反应。"""

    target: ReactionRequest = Field(
        ...,
        description="待推导反应焓的目标反应（反应物/生成物，须能唯一配平）",
    )
    known_reactions: list[HessKnownReaction] = Field(
        ...,
        min_length=1,
        max_length=6,
        description="已知反应，1～6 条；每条含反应物、生成物与以 kJ 表示的"
        "反应焓十进制字符串，均须能唯一配平；允许以反向或分数倍组合使用",
    )


class HessKnownReactionResult(BaseModel):
    """一条已知反应的配平回显与推导所用倍数。"""

    index: int = Field(description="从 0 开始的已知反应序号（与请求顺序一致）")
    reactants: list[str] = Field(description="规范化反应物（按规范物种排序）")
    products: list[str] = Field(description="规范化生成物（按规范物种排序）")
    equation: str = Field(description="规范化并配平后的化学方程式")
    enthalpy_kj: str = Field(description="客户端提交的反应焓（kJ，回显）")
    multiplier: str = Field(
        description="该已知反应在组合中的倍数，最简分数字符串（如 \"2\"、\"-1/2\"）；"
        "负值表示反向使用，分数表示分数倍使用"
    )
    multiplier_used: str = Field(
        description="实际采用的反应焓贡献 = 倍数 × 反应焓（kJ，十进制字符串）"
    )


class HessNetCoefficient(BaseModel):
    """合并抵消后单个物种的净系数。"""

    normalized: str = Field(description="规范化化学式")
    coefficient: str = Field(
        description="净计量系数，最简分数字符串：生成物为正、反应物为负，"
        "中间物种完全抵消时为 \"0\""
    )
    side: Literal["product", "reactant", "cancelled"] = Field(
        description="净系数为正属生成物侧 / 为负属反应物侧 / 为 0 表示完全抵消"
    )


class HessResponse(BaseModel):
    target: NormalizedReaction = Field(description="规范化并配平后的目标反应")
    known_reactions: list[HessKnownReactionResult] = Field(
        description="各已知反应的配平回显、反应焓与推导所用倍数（与请求顺序一致）"
    )
    derived_enthalpy_kj: str = Field(
        description="推导出的目标反应焓 ΔH（kJ，十进制字符串）；等于各已知"
        "反应 倍数 × 反应焓 之和，负号表示放热"
    )
    net_coefficients: list[HessNetCoefficient] = Field(
        description="全部涉及物种（目标 ∪ 已知反应）合并抵消后的逐物种净系数，"
        "按规范化化学式排序；含完全抵消为 0 的中间物种"
    )


# ------------------------------------------------------------ 同位素名义质量峰分布
class IsotopeItem(BaseModel):
    """一种同位素：正整数质量数与非负丰度十进制字符串。"""

    mass_number: int = Field(
        ...,
        ge=1,
        description="同位素名义质量数（正整数，如 ¹²C 的 12、³⁷Cl 的 37）",
        examples=[12],
    )
    abundance: str = Field(
        ...,
        min_length=1,
        description="该同位素的天然（或指定）丰度，**非负十进制字符串**"
        "（如 \"0.9893\"、\"0\"、\"1e-3\"），不要传 JSON 数字；"
        "同一元素各同位素丰度必须精确之和为 1",
        examples=["0.9893"],
    )

    @field_validator("abundance")
    @classmethod
    def _strip_abundance(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("丰度必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("丰度不能为空白字符串")
        return stripped


class ElementIsotopesRequest(BaseModel):
    """一种元素的同位素清单。"""

    element: str = Field(
        ...,
        min_length=1,
        description="元素符号（如 \"H\"、\"Cl\"），必须是化学式中实际出现的元素",
        examples=["H"],
    )
    isotopes: list[IsotopeItem] = Field(
        ...,
        min_length=1,
        max_length=3,
        description="该元素的同位素，1～3 种；质量数必须两两互异且为正整数，"
        "丰度均非负且精确之和为 1",
    )

    @field_validator("element")
    @classmethod
    def _strip_element(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("元素符号必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("元素符号不能为空白字符串")
        return stripped


class IsotopeEnvelopeRequest(BaseModel):
    """``POST /isotopes/envelope`` 请求：化学式 + 每种元素的同位素质量数与丰度。"""

    formula: str = Field(
        ...,
        min_length=1,
        description="待计算的化学式（支持既有解析器的全部写法：括号、下标、电荷；"
        "不接受电子物种 e-/e+）",
        examples=["H2O"],
    )
    elements: list[ElementIsotopesRequest] = Field(
        ...,
        min_length=1,
        max_length=3,
        description="每种出现元素恰好一组同位素定义，1～3 种元素；"
        "不得缺失、重复或给出化学式中不存在的元素",
    )

    @field_validator("formula")
    @classmethod
    def _strip_formula(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("化学式必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("化学式不能为空字符串")
        return stripped


class IsotopePeak(BaseModel):
    """一个名义质量峰。"""

    mass: int = Field(description="名义质量数（各原子质量数之和）")
    probability: str = Field(
        description="该名义质量出现的精确概率，最简分数字符串（如 \"3/4\"、\"1\"）"
    )


class IsotopeElementEcho(BaseModel):
    """输入元素及其原子数、同位素清单的回显。"""

    element: str = Field(description="元素符号")
    atom_count: int = Field(description="该元素在化学式中的原子数（由解析器取得）")
    isotopes: list[IsotopeItem] = Field(description="客户端提交的同位素清单（回显）")


class IsotopeEnvelopeResponse(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    total_atoms: int = Field(description="化学式总原子数")
    elements: list[IsotopeElementEcho] = Field(
        description="各元素的原子数与同位素清单（与请求 elements 顺序一致）"
    )
    peaks: list[IsotopePeak] = Field(
        description="全部名义质量峰，按质量升序；同质量组合已合并，"
        "概率为最简分数字符串，未截断低/零概率峰"
    )
    most_probable: IsotopePeak = Field(
        description="概率最高的峰；概率并列时取质量较小者"
    )


# ------------------------------------------------------------ 元素分析反推分子式
class FormulaElementMeasurement(BaseModel):
    """一种元素（C/H/N/O）的实测质量百分比与绝对容差。"""

    element: Literal["C", "H", "N", "O"] = Field(
        ...,
        description="元素符号，仅支持 C、H、N、O；请求中每种元素恰好一条记录",
        examples=["C"],
    )
    percent: str = Field(
        ...,
        min_length=1,
        description="该元素的实测质量百分比，**非负十进制字符串**，取值在 "
        "0～100 之间（如 \"75.0\"、\"0\"、\"1e-2\"），不要传 JSON 数字；"
        "四种元素的测量值不要求恰好合计 100",
        examples=["75.0"],
    )
    tolerance: str = Field(
        ...,
        min_length=1,
        description="该元素质量百分比的**绝对容差**（百分数），非负十进制"
        "字符串，取值在 0～100 之间（如 \"0.5\"、\"0\"），匹配为闭区间"
        "[percent−tolerance, percent+tolerance]，不要传 JSON 数字",
        examples=["0.5"],
    )

    @field_validator("percent", "tolerance")
    @classmethod
    def _strip_text(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError("必须是字符串")
        stripped = value.strip()
        if not stripped:
            raise ValueError("不能为空白字符串")
        return stripped


class MolarMassRange(BaseModel):
    """名义摩尔质量的正整数闭区间（按 C=12/H=1/N=14/O=16 计算）。"""

    minimum: int = Field(
        ...,
        ge=1,
        description="摩尔质量闭区间下界（正整数，含端点）",
        examples=[16],
    )
    maximum: int = Field(
        ...,
        ge=1,
        description="摩尔质量闭区间上界（正整数，含端点，不得小于 minimum）",
        examples=[20],
    )

    @field_validator("minimum", "maximum", mode="before")
    @classmethod
    def _reject_bool(cls, value: Any) -> Any:
        # Pydantic v2 默认会把 JSON 布尔值当作 0/1 整数接受；必须在其整数
        # 强转之前（mode=before）拦截，否则 isinstance(bool) 已无从识别
        if isinstance(value, bool):
            raise ValueError("摩尔质量必须是整数，不能是布尔值")
        return value


class FormulaInferRequest(BaseModel):
    """``POST /formulas/infer`` 请求：C/H/N/O 实测质量百分比与摩尔质量区间。"""

    elements: list[FormulaElementMeasurement] = Field(
        ...,
        min_length=4,
        max_length=4,
        description="C、H、N、O 四种元素各恰好一条测量记录（每条含实测质量"
        "百分比 percent 与绝对容差 tolerance，均为 0～100 的十进制字符串），"
        "不得缺失或重复；各测量值不要求合计 100",
    )
    molar_mass: MolarMassRange = Field(
        ...,
        description="名义摩尔质量的正整数闭区间 [minimum, maximum]"
        "（均为含端点，按 C=12、H=1、N=14、O=16 计算）",
    )


class FormulaMassPercentBlock(BaseModel):
    """一个候选分子式中四种元素各自的计算质量百分比（%，十进制字符串）。"""

    C: str = Field(description="碳的计算质量百分比 100·12·n_C/M")
    H: str = Field(description="氢的计算质量百分比 100·n_H/M")
    N: str = Field(description="氮的计算质量百分比 100·14·n_N/M")
    O: str = Field(description="氧的计算质量百分比 100·16·n_O/M")


class FormulaMatch(BaseModel):
    """一个落入全部百分比容差与摩尔质量区间的分子式。"""

    formula: str = Field(
        description="分子式（按 C/H/N/O 顺序渲染，计数为 0 的元素省略、为 1 省略下标）"
    )
    molar_mass: int = Field(description="名义摩尔质量 = 12·n_C + n_H + 14·n_N + 16·n_O")
    counts: dict[str, int] = Field(
        description="四种元素的原子计数（C/H/N/O 键齐全；总数大于零、逐项 0～30）"
    )
    mass_percent: FormulaMassPercentBlock = Field(
        description="该分子式各项计算质量百分比（%）"
    )
    total_abs_deviation: str = Field(
        description="总百分比绝对偏差 = Σ|计算百分比 − 实测百分比|"
        "（精确有理数渲染的十进制字符串），第一排序键"
    )


class FormulaInferResponse(BaseModel):
    elements: list[FormulaElementMeasurement] = Field(
        description="客户端提交的四种元素测量记录（回显，与请求顺序一致）"
    )
    molar_mass: MolarMassRange = Field(description="摩尔质量闭区间（回显）")
    count: int = Field(description="本次实际返回的匹配分子式数（≤ 20）")
    total: int = Field(description="落入全部约束的匹配分子式总数（截断前）")
    truncated: bool = Field(description="匹配总数超过 20、已截取前 20 个时为 true")
    matches: list[FormulaMatch] = Field(
        description="按总百分比绝对偏差、摩尔质量、C/H/N/O 计数依次升序排列的"
        "前 20 个匹配分子式；无解时为空列表"
    )
