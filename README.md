# 化学反应计量平衡服务

基于 **Python + FastAPI** 的化学反应方程式配平微服务。客户端提交包含
**元素、括号与电荷**的反应物 / 生成物分子式，服务端：

1. 用递归下降解析器解析化学式（支持三种括号、嵌套、下标乘数、离子电荷、电子 `e-`）；
2. 建立 **元素守恒 + 总电荷守恒** 矩阵（反应物列取正、生成物列取负）；
3. 用 `fractions.Fraction` 做精确 **高斯-若尔当消元（RREF）** 求零空间，
   把唯一的有理数解放大、约分得到 **最小互质正整数系数**；
4. `/balance/variants` 在自由变量不超过 3 个时，枚举零空间内所有
   `1..max_coefficient` 的互质正整数配平方案；
5. `/stoichiometry/batch` 在唯一配平之上，用 `Decimal` 精确计算反应进度、
   并列限量试剂、各反应物消耗/剩余与各生成物理论产量；
6. `/stoichiometry/plan` 做反向投料规划：按一个或多个生成物的目标摩尔数，
   取各目标量/配平系数的最大值得到满足全部目标的最小反应进度，给出各反应物
   理论需求量、各目标产物达成量与超额量；
7. `/stoichiometry/yield` 计算实际收率：提交各反应物投料摩尔数与每种生成物
   的实际摩尔产量，复用唯一配平与精确计量逻辑计算反应进度与理论产量，逐种
   返回理论产量、实际产量与收率百分比（理论产量为 0 时收率返回 `null` 并说明
   原因，实际产量高于理论值时保留大于 100% 的结果）；
8. `/stoichiometry/sequence` 做连续反应推演：按顺序执行多步反应，从初始
   物种摩尔库存出发，每步以全部反应物的当前库存计算最大反应进度，扣减消耗、
   加入产物并把更新后的库存传给下一步（库存不足则该步零进度继续并标出全部
   阻塞反应物），全程精确有理数运算；
9. `/inventory/allocate` 做试剂批次 **FEFO 领用分配**：提交待领用物种的有效
   摩尔需求与带批号/纯度/名义摩尔数/失效日期/状态的库存批次及核算日期，
   服务按失效日期优先、同日批号升序跨批次扣减（名义量 × 纯度换算有效物质），
   跳过已过期与隔离批次，返回每批领用/余量、未满足需求与汇总状态；
10. `/solutions/blend` 做**溶液配制**：提交不超过 4 种母液（各溶质浓度 +
    可用体积）、目标终体积与各溶质目标浓度（均为非负十进制字符串，浓度
    mol/L、体积 mL），以精确有理数求满足全部目标浓度的非负母液取用量
    （受可用体积与终体积约束，余量用水补足；多解时先取母液总量最小、再按
    输入顺序取用量字典序最小），返回各母液取用量、补水量与最终浓度；
11. `/isotopes/envelope` 做**同位素名义质量峰分布**：提交化学式与每种元素
    的同位素质量数、丰度（非负十进制字符串，每组和恰为 1），按独立原子
    抽样（多项分布）枚举全部可能的名义质量与精确概率，同质量组合合并，
    以最简分数字符串返回按质量升序的全部峰与最高概率峰（并列取质量较小者）；
12. `/thermochemistry/hess` 做**赫斯定律反应焓推导**：提交目标反应与至多
    6 条已知反应（各含反应物、生成物与以 kJ 表示、可带符号的反应焓十进制
    字符串），复用唯一配平按规范物种建立生成物为正、反应物为负的计量向量，
    以精确有理数求已知向量的有符号线性组合（允许反向与分数倍使用）复现
    目标反应，返回各已知反应倍数、推导出的目标反应焓与合并抵消后的逐物种
    净系数；目标不在张成空间时报不可推导，多种组合导出不同反应焓时报歧义；
13. `/equilibrium/solve` 做**理想溶液平衡组成计算**：提交一个可逆反应、
    各参与物种（反应物 ∪ 生成物）的**正**初始摩尔数、溶液体积（L）与
    浓度平衡常数 Kc（均为正十进制字符串），复用唯一配平，以配平系数确定
    双向反应进度的非负物料区间，在区间内按各物种浓度的计量次方构造 Qc，
    精确有理二分求 Qc=Kc 的唯一内点（初始 Qc>Kc 时返回**负进度**），
    返回反应进度（绝对误差上界 1e-9 mol）、物料区间、初始/平衡反应商与
    各物种平衡摩尔数；
14. `/formulas/infer` 做**元素分析反推分子式**：提交 C、H、N、O 四种
    元素的实测质量百分比与绝对容差（均为 0～100 的十进制字符串，测量值
    不要求合计 100）以及摩尔质量的正整数闭区间，按 C=12、H=1、N=14、
    O=16 以精确有理数枚举逐元素原子数 0～30、总数大于零且**同时**落入
    全部百分比容差（闭区间）与质量区间的整数计数，按总百分比绝对偏差、
    摩尔质量、C/H/N/O 计数依次升序返回前 20 个分子式（含各项计算百分比、
    匹配总数与截断标志），无解返回空列表；
15. 对同一规范化反应做幂等缓存，重复请求返回相同的 `reaction_id` 与结果。

---

## 1. 一条命令启动

需要已安装 Docker 与 Docker Compose v2：

```bash
docker compose up --build -d
```

服务监听 **宿主机端口 `8000`**（映射容器 8000）。

| 入口 | 地址 |
| --- | --- |
| 健康检查 | <http://localhost:8000/health> |
| Swagger 交互文档 | <http://localhost:8000/docs> |
| ReDoc 文档 | <http://localhost:8000/redoc> |

> 想换宿主机端口：编辑 `docker-compose.yml` 的 `ports`，例如
> `"9000:8000"` 后访问 9000 端口。容器内固定监听 8000。

停止服务：

```bash
docker compose down
```

启动后可直接调用批量计量接口（宿主机端口 **8000**）：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/batch \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["CH4","O2"],"products":["CO2","H2O"],
       "amounts":[{"species":"CH4","moles":"1"},{"species":"O2","moles":"3"}]}'
```

反向投料规划——要求至少产出 4 mol NH₃，看需要多少 N₂ / H₂：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/plan \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["N2","H2"],"products":["NH3"],
       "targets":[{"species":"NH3","moles":"4"}]}'
```

实际收率——按投料与实际产出的 NH₃ 计算理论产量与收率百分比：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/yield \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["N2","H2"],"products":["NH3"],
       "amounts":[{"species":"N2","moles":"1"},{"species":"H2","moles":"3"}],
       "actual_yields":[{"species":"NH3","moles":"1.5"}]}'
```

试剂批次领用分配——按失效日期优先（FEFO）跨批领用 1.5 mol 有效 HCl：

```bash
curl -X POST http://localhost:8000/api/v1/inventory/allocate \
  -H 'Content-Type: application/json' \
  -d '{
    "demands": [{"species": "HCl", "effective_moles": "1.5"}],
    "batches": [
      {"batch_no": "A-1", "species": "HCl", "purity": "0.75",
       "nominal_moles": "1", "expires_on": "2026-10-01", "status": "available"},
      {"batch_no": "A-2", "species": "HCl", "purity": "1",
       "nominal_moles": "2", "expires_on": "2026-11-01", "status": "available"},
      {"batch_no": "OLD", "species": "HCl", "purity": "1",
       "nominal_moles": "8", "expires_on": "2020-01-01", "status": "available"}
    ],
    "accounting_date": "2026-09-23"
  }'
```

溶液配制——用 1.0 mol/L 的 NaCl 母液配 250 mL 0.1 mol/L 目标液（取 25 mL
母液、补水 225 mL）：

```bash
curl -X POST http://localhost:8000/api/v1/solutions/blend \
  -H 'Content-Type: application/json' \
  -d '{
    "stocks": [
      {"id": "stock-A", "available_volume": "100",
       "solutes": [{"solute": "NaCl", "concentration": "1.0"}]}
    ],
    "target_volume": "250",
    "targets": [{"solute": "NaCl", "concentration": "0.1"}]
  }'
```

同位素名义质量峰分布——HCl（氢 1H/2H 各 0.5，氯 35Cl 0.75、37Cl 0.25）：

```bash
curl -X POST http://localhost:8000/api/v1/isotopes/envelope \
  -H 'Content-Type: application/json' \
  -d '{
    "formula": "HCl",
    "elements": [
      {"element": "H", "isotopes": [
        {"mass_number": 1, "abundance": "0.5"},
        {"mass_number": 2, "abundance": "0.5"}]},
      {"element": "Cl", "isotopes": [
        {"mass_number": 35, "abundance": "0.75"},
        {"mass_number": 37, "abundance": "0.25"}]}
    ]
  }'
```

赫斯定律反应焓推导——由 C 与 H₂ 的燃烧焓、CH₄ 的燃烧焓求 CH₄ 标准生成焓
（第 3 条燃烧反应需反向使用，ΔH = −393.5 + (−571.6) + 890.3 = −74.8 kJ）：

```bash
curl -X POST http://localhost:8000/api/v1/thermochemistry/hess \
  -H 'Content-Type: application/json' \
  -d '{
    "target": {"reactants": ["C", "H2"], "products": ["CH4"]},
    "known_reactions": [
      {"reactants": ["C", "O2"], "products": ["CO2"], "enthalpy_kj": "-393.5"},
      {"reactants": ["H2", "O2"], "products": ["H2O"], "enthalpy_kj": "-571.6"},
      {"reactants": ["CH4", "O2"], "products": ["CO2", "H2O"], "enthalpy_kj": "-890.3"}
    ]
  }'
```

理想溶液平衡组成——`PCl5 ⇌ PCl3 + Cl2`，1 L 中初始 1 mol PCl5、各
0.01 mol 产物，Kc=0.8：

```bash
curl -X POST http://localhost:8000/api/v1/equilibrium/solve \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["PCl5"],
    "products": ["PCl3", "Cl2"],
    "amounts": [
      {"species": "PCl5", "moles": "1"},
      {"species": "PCl3", "moles": "0.01"},
      {"species": "Cl2",  "moles": "0.01"}
    ],
    "volume_l": "1",
    "kc": "0.8"
  }'
```

元素分析反推分子式——葡萄糖元素分析结果（C 40.00%、H 6.71%、O
53.29%，N 未检出，各自 ±0.05；测量值合计不必为 100），摩尔质量
170～190：

```bash
curl -X POST http://localhost:8000/api/v1/formulas/infer \
  -H 'Content-Type: application/json' \
  -d '{
    "elements": [
      {"element": "C", "percent": "40.00", "tolerance": "0.05"},
      {"element": "H", "percent": "6.71",  "tolerance": "0.05"},
      {"element": "N", "percent": "0",     "tolerance": "0.05"},
      {"element": "O", "percent": "53.29", "tolerance": "0.05"}
    ],
    "molar_mass": {"minimum": 170, "maximum": 190}
  }'
```

### 不用 Docker 的本地运行方式

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install -r requirements-dev.txt
pytest          # 全部测试（解析器/配平/计量/库存/溶液/同位素/元素分析/赫斯定律）
```

---

## 2. API 一览（前缀 `/api/v1`）

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/v1/balance` | 配平反应，返回最小正整数系数与电荷守恒结果 |
| POST | `/api/v1/balance/variants` | 枚举自由变量 ≤3 时所有互质正整数配平方案 |
| POST | `/api/v1/analyze` | 逐元素差异说明：配平前（单位系数）与配平后两个视图 |
| POST | `/api/v1/stoichiometry/batch` | 投料摩尔数批量计量：反应进度、并列限量试剂、消耗/剩余与理论产量 |
| POST | `/api/v1/stoichiometry/plan` | 反向投料规划：按生成物目标摩尔数求最小反应进度、理论投料需求与达成/超额量 |
| POST | `/api/v1/stoichiometry/yield` | 实际收率：投料摩尔数与各生成物实际产量 → 逐种理论产量、实际产量与收率百分比 |
| POST | `/api/v1/stoichiometry/sequence` | 连续反应推演：多步反应按顺序消耗/产出，逐步传递库存，返回各步明细与最终库存 |
| POST | `/api/v1/inventory/allocate` | 试剂批次 FEFO 领用分配：按失效日期优先、同日批号升序跨批次扣减，跳过过期/隔离批次，返回每批领用/余量、未满足需求与汇总状态 |
| POST | `/api/v1/solutions/blend` | 溶液配制：不超过 4 种母液 + 目标终体积与各溶质目标浓度 → 满足全部目标浓度的最小母液取用量、补水量与最终浓度 |
| POST | `/api/v1/isotopes/envelope` | 同位素名义质量峰分布：化学式 + 每元素同位素质量数/丰度 → 按质量升序的全部名义质量峰（最简分数概率）与最高概率峰 |
| POST | `/api/v1/thermochemistry/hess` | 赫斯定律反应焓推导：目标反应 + 至多 6 条已知反应（含 kJ 反应焓）→ 各反应倍数（可反向/分数倍）、推导反应焓与合并抵消后的逐物种净系数 |
| POST | `/api/v1/equilibrium/solve` | 理想溶液平衡组成：可逆反应 + 各参与物种正初始摩尔数 + 体积(L) 与 Kc → Qc=Kc 的唯一内点反应进度（初始 Qc>Kc 时为负进度）、物料区间与各物种平衡摩尔数（进度绝对误差上界 1e-9 mol） |
| POST | `/api/v1/formulas/infer` | 元素分析反推分子式：C/H/N/O 实测质量百分比与绝对容差（十进制字符串，不要求合计 100）+ 摩尔质量正整数闭区间 → 同时落入全部容差与质量区间的前 20 个分子式（按总偏差、摩尔质量、计数升序）及匹配总数/截断标志 |
| GET | `/api/v1/reactions/{reaction_id}` | 按幂等 ID 取回之前的结果 |
| GET | `/health` | 健康检查 |

### 2.1 配平 `POST /api/v1/balance`

请求体：

```json
{
  "reactants": ["KMnO4", "HCl"],
  "products": ["KCl", "MnCl2", "Cl2", "H2O"]
}
```

响应 `200`：

```json
{
  "reaction_id": "rxn_7ed134d0c6ec5699",
  "cached": false,
  "normalized": {
    "reactants": ["HCl", "KMnO4"],
    "products": ["KCl", "Cl2", "MnCl2", "H2O"],
    "equation": "16HCl + 2KMnO4 -> 2KCl + 5Cl2 + 2MnCl2 + 8H2O"
  },
  "coefficients": {
    "reactants": [
      {"input": "HCl", "normalized": "HCl", "coefficient": 16},
      {"input": "KMnO4", "normalized": "KMnO4", "coefficient": 2}
    ],
    "products": [
      {"input": "KCl", "normalized": "KCl", "coefficient": 2},
      {"input": "Cl2", "normalized": "Cl2", "coefficient": 5},
      {"input": "MnCl2", "normalized": "MnCl2", "coefficient": 2},
      {"input": "H2O", "normalized": "H2O", "coefficient": 8}
    ],
    "vector": [16, 2, 2, 5, 2, 8]
  },
  "balanced": true,
  "elements": ["Cl", "H", "K", "Mn", "O"],
  "charge_balance": {"reactant_total": 0, "product_total": 0, "difference": 0, "balanced": true}
}
```

`curl`：

```bash
curl -X POST http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["KMnO4","HCl"],"products":["KCl","MnCl2","Cl2","H2O"]}'
```

### 2.2 逐元素差异 `POST /api/v1/analyze`

返回两个视图：

* `raw`：所有物种按系数 1 统计（即“配平前”的差异，便于解释哪里不平衡）；
* `balanced`：用求得的系数统计，此时每个元素与电荷的 `difference` 均为 0。

以高锰酸根酸性还原半反应为例：

```bash
curl -X POST http://localhost:8000/api/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["MnO4-","H+","e-"],"products":["Mn2+","H2O"]}'
```

`raw.charge` 为 `{"reactant_total": -1, "product_total": 2, "difference": -3}`，
`balanced` 视图给出 `MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O`，电荷两侧均为 +2。

即使反应**无法唯一配平**，`/analyze` 也返回 `200`，同时带 `raw` 视图与
`error` 字段（`balanced_ok=false`、`balanced=null`），方便客户端展示诊断信息。

### 2.3 枚举所有配平方案 `POST /api/v1/balance/variants`

`/balance` 只返回唯一配平；当守恒方程存在自由变量（如不完整的燃烧反应）时，
本接口在零空间中枚举**所有**满足条件的正整数方案：

* 每个系数都在 `1 ~ max_coefficient`（**2～30**，必填）；
* 元素与总电荷守恒（由零空间解天然保证）；
* 系数最大公约数为 1（只保留本原解，`2H2+O2→2H2O` 与 `4H2+2O2→4H2O`
  视为同一方案）。

零空间自由变量（= 零空间维数 − 1）超过 **3** 个（即维数 > 3）时，枚举空间
过大，返回 `422 VARIANT_SPACE_TOO_LARGE`。结果先按**系数总和**升序、再按
**规范化物种顺序的系数向量字典序**排序，取前 `limit`（**1～100**，必填）条，
并用 `truncated` 标识是否发生截取。唯一配平的反应同样适用，返回只含单个
方案的列表。

```bash
curl -X POST http://localhost:8000/api/v1/balance/variants \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["C","O2"],"products":["CO","CO2"],"max_coefficient":10,"limit":5}'
```

`200` 响应：

```json
{
  "reaction_id": "rxn_aa766b24597bded4",
  "normalized": {
    "reactants": ["C", "O2"],
    "products": ["CO", "CO2"],
    "equation": "C + O2 -> CO + CO2"
  },
  "max_coefficient": 10,
  "limit": 5,
  "nullity": 2,
  "count": 5,
  "total": 15,
  "truncated": true,
  "variants": [
    {
      "coefficients": {
        "reactants": [
          {"input": "C", "normalized": "C", "coefficient": 3},
          {"input": "O2", "normalized": "O2", "coefficient": 2}
        ],
        "products": [
          {"input": "CO", "normalized": "CO", "coefficient": 2},
          {"input": "CO2", "normalized": "CO2", "coefficient": 1}
        ],
        "vector": [3, 2, 2, 1]
      },
      "sum": 8,
      "equation": "3C + 2O2 -> 2CO + CO2"
    }
  ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `nullity` | 守恒矩阵零空间维数（本接口只返回 1～3） |
| `total` | 在 `max_coefficient` 范围内枚举出的方案总数（截取前） |
| `count` | 本次实际返回的条数（`total` 与 `limit` 的较小值） |
| `truncated` | `total > limit` 时为 `true` |
| `variants[].sum` | 该方案系数总和，即第一排序键 |
| `variants[].coefficients.vector` | 按规范化后 reactants + products 顺序的系数向量 |

不存在正整数解时（如唯一解要求某物种系数为 0）返回 `200` 且
`variants=[]`、`total=0`；守恒本身无解返回 `422 IMPOSSIBLE_BALANCE`。
`reaction_id` 与 `/balance`、`/analyze` 采用同一规范化签名（变体结果不写
入幂等缓存）。

### 2.4 批量限量试剂计量 `POST /api/v1/stoichiometry/batch`

在唯一配平的基础上，给出每种**反应物的投料摩尔数**，服务端以精确十进制
（`Decimal`，全程不经过 `float`）计算：

* 反应进度 ξ = min(**投料量 ÷ 配平系数**)；
* **全部并列限量试剂**——当投料恰好按系数比时可能多于一种（全部投料为 0 时全部限量）；
* 每种反应物的投料量 `fed`、消耗量 `consumed = ξ·系数` 与剩余量
  `remaining = fed − consumed`；
* 每种生成物的理论产量 `theoretical_yield = ξ·系数`。

请求体在 `reactants` / `products` 之外增加 `amounts`，**每种规范化反应物
恰好一条记录**（等价写法可识别，但不能缺失、重复或出现未知物种）；
摩尔数必须是**非负十进制字符串**（不要传 JSON 数字，以避免浮点误差）：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/batch \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["CH4", "O2"],
    "products": ["CO2", "H2O"],
    "amounts": [
      {"species": "CH4", "moles": "1"},
      {"species": "O2",  "moles": "3"}
    ]
  }'
```

`200` 响应（十进制值均为字符串；物种按规范化顺序排列，调换书写顺序结果不变）：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["CH4", "O2"],
    "products": ["CO2", "H2O"],
    "equation": "CH4 + 2O2 -> CO2 + 2H2O"
  },
  "reaction_extent": {"value": "1", "unit": "mol"},
  "limiting_reactants": ["CH4"],
  "reactants": [
    {"input": "CH4", "normalized": "CH4", "coefficient": 1,
     "fed": "1", "consumed": "1", "remaining": "0", "is_limiting": true},
    {"input": "O2", "normalized": "O2", "coefficient": 2,
     "fed": "3", "consumed": "2", "remaining": "1", "is_limiting": false}
  ],
  "products": [
    {"input": "CO2", "normalized": "CO2", "coefficient": 1, "theoretical_yield": "1"},
    {"input": "H2O", "normalized": "H2O", "coefficient": 2, "theoretical_yield": "2"}
  ]
}
```

投料比恰好等于系数比（如 `CH4=2.5, O2=5`）时 `limiting_reactants` 同时返回
`["CH4", "O2"]`，二者 `remaining` 都是 `"0"`。

**无限小数**（分母含 2、5 以外的因子，如 ξ = 1/7）按 **28 位有效数字**
银行家舍入输出，例如 `"0.1428571428571428571428571429"`；有限小数给出精确值。

错误沿用统一结构化格式：

| 场景 | HTTP / code |
| --- | --- |
| 反应式或投料物种无法解析 | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 投料缺失、重复、含未知物种（详情在 `details.missing/duplicate/unknown`） | `422 INVALID_REQUEST` |
| 摩尔数非法：负数、`NaN`、`Infinity`、非数字字符串或 JSON 数字（详情 `details.invalid_moles`） | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

`reaction_id` 与 `/balance`、`/analyze` 共享同一规范化签名——同一反应调换
反应物 / 生成物书写顺序或使用等价化学式，得到相同的 `reaction_id` 与完全
一致的计算结果（仅 `input` 字段回显原始写法）。

### 2.5 反向投料规划 `POST /api/v1/stoichiometry/plan`

与 `/stoichiometry/batch` 方向相反：客户端为**一个或多个生成物**提交目标
摩尔数（允许只覆盖部分生成物），服务端以精确十进制（全程不经过 `float`）
计算：

* 满足全部目标所需的**最小反应进度** ξ = max(**目标量 ÷ 配平系数**)
  （取最大值是为了让每个目标产物的达成量都不低于其目标量）；
* 全部并列取得最大值的**驱动目标产物** `driving_targets`——目标比恰好
  等于系数比时可多于一种；
* 每种反应物的**理论需求量** `required = ξ·系数`；
* 每个提交目标的产物：目标量 `target`、达成量 `achieved = ξ·系数`、
  **超额量** `excess = achieved − target`（驱动产物为 `"0"`，其余为正）；
* 全部生成物（含未提交目标的）在该进度下的达成量。

请求体用 `targets` 代替 `amounts`；目标物种必须是规范化生成物之一
（等价写法可识别），不能重复或出现未知物种，摩尔数同样必须是非负十进制
字符串：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/plan \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["NH3", "O2"],
    "products": ["NO", "H2O"],
    "targets": [
      {"species": "NO",  "moles": "1"},
      {"species": "H2O", "moles": "6"}
    ]
  }'
```

`4NH3 + 5O2 -> 4NO + 6H2O` 中两个目标比值分别为 1/4 与 6/6，ξ 取最大值
`1`（由 `H2O` 驱动），因此 `NO` 实际达成 4 mol、超额 3 mol：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["NH3", "O2"],
    "products": ["NO", "H2O"],
    "equation": "4NH3 + 5O2 -> 4NO + 6H2O"
  },
  "reaction_extent": {"value": "1", "unit": "mol"},
  "driving_targets": ["H2O"],
  "reactants": [
    {"input": "NH3", "normalized": "NH3", "coefficient": 4, "required": "4"},
    {"input": "O2",  "normalized": "O2",  "coefficient": 5, "required": "5"}
  ],
  "targets": [
    {"input": "NO", "normalized": "NO", "coefficient": 4,
     "target": "1", "achieved": "4", "excess": "3", "is_driving": false},
    {"input": "H2O", "normalized": "H2O", "coefficient": 6,
     "target": "6", "achieved": "6", "excess": "0", "is_driving": true}
  ],
  "products": [
    {"input": "NO", "normalized": "NO", "coefficient": 4, "achieved": "4"},
    {"input": "H2O", "normalized": "H2O", "coefficient": 6, "achieved": "6"}
  ]
}
```

单目标示例：`3H2 + N2 -> 2NH3`，目标 `NH3 = 4` mol 时 ξ = 2，理论需求
`H2 = 6` mol、`N2 = 2` mol，达成量恰为 4 mol、超额为 0。目标摩尔数为
`"0"` 时 ξ 与所有理论需求均为 `"0"`。无限小数同样按 **28 位有效数字**
银行家舍入输出（如目标 `CO2 = 1` mol 于 `Fe2O3 + 3CO -> 2Fe + 3CO2`
时 ξ = 1/3 → `"0.3333333333333333333333333333"`）。

错误结构与 `/stoichiometry/batch` 一致：

| 场景 | HTTP / code |
| --- | --- |
| 反应式或目标物种无法解析 | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 目标重复或包含非生成物物种（详情 `details.duplicate/unknown`；允许只提交部分生成物） | `422 INVALID_REQUEST` |
| 摩尔数非法：负数、`NaN`、`Infinity`、非数字字符串或 JSON 数字（详情 `details.invalid_moles`） | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

### 2.6 实际收率 `POST /api/v1/stoichiometry/yield`

提交反应式、**每种反应物的投料摩尔数**（`amounts`，与
`/stoichiometry/batch` 完全一致的一一对应要求）以及**每种生成物的实际
摩尔产量**（`actual_yields`，同样与规范化生成物一一对应）。服务端复用
唯一配平与精确计量逻辑：

* 反应进度 ξ = min(**投料量 ÷ 配平系数**)、并列限量试剂与各反应物
  消耗/剩余，计算过程与结果同 `/stoichiometry/batch` 完全一致；
* 每种生成物的**理论产量** = ξ·系数、**实际产量**（回显请求值）与
  **收率百分比** = 实际产量 ÷ 理论产量 × 100；
* 物种按规范化签名匹配：`H2 O`、`[H]2` 等等价写法会识别为同一物种
  （`input` 回显原始写法，`normalized` 为规范化化学式）；
* 理论产量为 **0**（ξ = 0，即至少一种反应物投料为 0）时收率无定义，
  `yield_percent` 返回 `null`，并在 `yield_note` 中说明原因及投料为 0
  的全部限量反应物；
* 实际产量高于理论值时**不做封顶**，`yield_percent` 保留大于 100 的结果
  （如 125），请结合实验误差自行判断；
* 全程精确有理数运算（不经过 `float`），无限小数按 28 位有效数字
  银行家舍入输出，规则同其它计量接口。

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/yield \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["N2", "H2"],
    "products": ["NH3"],
    "amounts": [
      {"species": "N2", "moles": "1"},
      {"species": "H2", "moles": "3"}
    ],
    "actual_yields": [
      {"species": "NH3", "moles": "1.5"}
    ]
  }'
```

`N2 + 3H2 -> 2NH3`、ξ = 1 时理论 NH₃ = 2 mol，实际 1.5 mol，收率 75%：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["H2", "N2"],
    "products": ["NH3"],
    "equation": "N2 + 3H2 -> 2NH3"
  },
  "reaction_extent": {"value": "1", "unit": "mol"},
  "limiting_reactants": ["N2"],
  "reactants": [
    {"input": "H2", "normalized": "H2", "coefficient": 3,
     "fed": "3", "consumed": "3", "remaining": "0", "is_limiting": false},
    {"input": "N2", "normalized": "N2", "coefficient": 1,
     "fed": "1", "consumed": "1", "remaining": "0", "is_limiting": true}
  ],
  "products": [
    {"input": "NH3", "normalized": "NH3", "coefficient": 2,
     "theoretical_yield": "2", "actual_yield": "1.5",
     "yield_percent": "75", "yield_note": null}
  ]
}
```

零投料导致理论产量为 0 的示例（N2 投料 `"0"`）：

```json
"products": [
  {"input": "NH3", "normalized": "NH3", "coefficient": 2,
   "theoretical_yield": "0", "actual_yield": "0",
   "yield_percent": null,
   "yield_note": "理论产量为 0：反应进度 ξ = min(投料量 ÷ 系数) = 0，投料为 0 的限量反应物为 ['N2']，收率百分比无定义"}
]
```

错误结构与 `/stoichiometry/batch` 一致，投料侧问题沿用
`details.missing/duplicate/unknown/invalid_moles`，实际产量侧问题以
`actual_yields.` 前缀聚合：

| 场景 | HTTP / code |
| --- | --- |
| 反应式、投料或实际产量物种无法解析（`details.field` 指出 `amounts.species` / `actual_yields.species`） | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 投料缺失/重复/含未知物种（详情 `details.missing/duplicate/unknown`） | `422 INVALID_REQUEST` |
| 实际产量缺失/重复/含非生成物物种（详情 `details.actual_yields.missing/duplicate/unknown`） | `422 INVALID_REQUEST` |
| 任一侧摩尔数非法：负数、`NaN`、`Infinity`、非数字字符串或 JSON 数字 | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

`reaction_id` 与其它计量接口共享同一规范化签名。

### 2.7 连续反应推演 `POST /api/v1/stoichiometry/sequence`

给定**按顺序执行的多步反应**与**初始物种摩尔库存**，服务端逐步推演：

* 每一步先复用与 `/balance` 相同的解析与**唯一配平**（及幂等缓存）；
* 以该步**全部反应物的当前库存**（未列出的物种视为 0）计算最大反应进度
  ξ = min(**当前库存 ÷ 配平系数**)；
* 扣减各反应物 ξ·系数、向各生成物加入 ξ·系数，更新后的库存传给下一步
  （上一步的产物可直接成为下一步的反应物；产物在初始库存中已有时按
  初始 + 此前累计计算）；
* 库存不足以让反应发生（ξ = 0）时该步以**零进度继续**，库存不发生任何
  变化，并在 `blocking_reactants` 中列出库存为 0 的**全部**反应物；
* 等价化学式（括号风格、空白、电荷写法不同）在库存中**合并为同一物种**，
  摩尔数精确相加；库存全程以精确有理数传递，仅输出时按十进制字符串渲染
  （无限小数 28 位有效数字银行家舍入，规则同 `/stoichiometry/batch`）。

请求体：

```json
{
  "steps": [
    {"reactants": ["FeS2", "O2"], "products": ["Fe2O3", "SO2"]},
    {"reactants": ["SO2", "O2"], "products": ["SO3"]},
    {"reactants": ["SO3", "H2O"], "products": ["H2SO4"]}
  ],
  "initial_inventory": [
    {"species": "FeS2", "moles": "4"},
    {"species": "O2",   "moles": "20"},
    {"species": "H2O",  "moles": "5"}
  ]
}
```

`curl`：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/sequence \
  -H 'Content-Type: application/json' \
  -d '{"steps":[{"reactants":["FeS2","O2"],"products":["Fe2O3","SO2"]},
                {"reactants":["SO2","O2"],"products":["SO3"]},
                {"reactants":["SO3","H2O"],"products":["H2SO4"]}],
       "initial_inventory":[{"species":"FeS2","moles":"4"},
                            {"species":"O2","moles":"20"},
                            {"species":"H2O","moles":"5"}]}'
```

三步配平式分别为 `4FeS2 + 11O2 -> 2Fe2O3 + 8SO2`、`O2 + 2SO2 -> 2SO3`、
`H2O + SO3 -> H2SO4`；反应进度依次为 1、4、5，最终
`H2SO4=5、SO3=3、O2=5、Fe2O3=2`，耗尽的 `FeS2/SO2/H2O` 为 `"0"`。

`200` 响应（节选一步）：

```json
{
  "steps": [
    {
      "index": 0,
      "reaction_id": "rxn_a584b7435c0643ca",
      "cached": false,
      "normalized": {
        "reactants": ["FeS2", "O2"],
        "products": ["Fe2O3", "SO2"],
        "equation": "4FeS2 + 11O2 -> 2Fe2O3 + 8SO2"
      },
      "reaction_extent": {"value": "1", "unit": "mol"},
      "blocked": false,
      "blocking_reactants": [],
      "reactants": [
        {"input": "FeS2", "normalized": "FeS2", "coefficient": 4,
         "before": "4", "consumed": "4", "after": "0", "is_blocking": false},
        {"input": "O2", "normalized": "O2", "coefficient": 11,
         "before": "20", "consumed": "11", "after": "9", "is_blocking": false}
      ],
      "products": [
        {"input": "Fe2O3", "normalized": "Fe2O3", "coefficient": 2,
         "before": "0", "produced": "2", "after": "2"},
        {"input": "SO2", "normalized": "SO2", "coefficient": 8,
         "before": "0", "produced": "8", "after": "8"}
      ],
      "inventory": [
        {"normalized": "Fe2O3", "moles": "2"},
        {"normalized": "FeS2", "moles": "0"},
        {"normalized": "H2O", "moles": "5"},
        {"normalized": "O2", "moles": "9"},
        {"normalized": "SO2", "moles": "8"}
      ]
    }
  ],
  "final_inventory": [ /* 与最后一步 inventory 相同 */ ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `steps[].index` | 从 0 开始的步骤序号，与请求 `steps` 顺序一致 |
| `steps[].reaction_id` / `cached` | 该步反应的幂等 ID（与 `/balance` 共享签名）与缓存命中标记 |
| `steps[].reaction_extent` | 该步按当前库存可达的最大反应进度 ξ |
| `steps[].blocked` | ξ = 0（库存不足）时为 `true`，库存原样传递到下一步 |
| `steps[].blocking_reactants` | 零进度时库存为 0 的全部反应物规范化化学式；正常进度时为 `[]` |
| `steps[].reactants[]` | 各反应物该步 `before / consumed / after` 库存与 `is_blocking` |
| `steps[].products[]` | 各生成物该步 `before / produced / after` 库存 |
| `steps[].inventory` | 该步结束后的全部物种库存快照（含库存为 0 的物种） |
| `final_inventory` | 全部步骤执行完毕后的最终库存，等于最后一步的 `inventory` |

库存条目按规范化化学式名字排序；`steps[].inventory` / `final_inventory`
包含推演中出现过的全部物种（含恰好耗尽为 `"0"` 者），便于核对物料去向。

错误沿用统一结构化格式，并在 `details.step` 中标出从 0 开始的出错步骤：

| 场景 | HTTP / code |
| --- | --- |
| 某步反应式无法解析（`details.step` 指出步骤） | `400 FORMULA_PARSE_ERROR` |
| 初始库存物种无法解析（`details.field=initial_inventory.species`） | `400 FORMULA_PARSE_ERROR` |
| 某步同侧物种重复（`details.step`） | `409 DUPLICATE_SPECIES` |
| 某步无法唯一配平 / 无法配平（`details.step`） | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |
| `steps` 为空、缺字段、化学式为空白等请求体问题 | `422 INVALID_REQUEST` |
| 初始库存摩尔数非法（详情 `details.invalid_moles`） | `422 INVALID_REQUEST` |

### 2.8 试剂批次领用分配 `POST /api/v1/inventory/allocate`

提交**待领用物种的有效摩尔需求**、一组库存批次（批号、纯度、可用名义摩尔数、
失效日期、状态）与**核算日期**，服务按 FEFO（先到期先出）规则分配：

* 同一物种的可用批次按 **失效日期升序**、同日按 **批号字符串升序** 排序后依次领用；
* 每批的有效物质 = **名义摩尔数 × 纯度**；从一批领用的名义摩尔数 =
  有效领用 ÷ 纯度，因此不会超过该批有效余量；
* **已过期**（`expires_on < accounting_date`，恰为核算日当天不算过期）、
  **隔离**（`status="quarantined"`）与**可用名义量为 0** 的批次跳过，不参与扣减；
* 需求可**跨多个批次拆分**；批号在单次请求内必须全局唯一（同一批次不得重复扣减）；
* 物种不在 `demands` 中的批次整批保留，单独列入 `unmatched_batches`；
* 所有数值以精确有理数运算（全程不经过 `float`），输出为十进制字符串；
  无限小数按 28 位有效数字银行家舍入，规则同其它计量接口。

请求体：

```json
{
  "demands": [
    {"species": "HCl", "effective_moles": "1.5"}
  ],
  "batches": [
    {"batch_no": "A-1", "species": "HCl", "purity": "0.75",
     "nominal_moles": "1", "expires_on": "2026-10-01", "status": "available"},
    {"batch_no": "A-2", "species": "HCl", "purity": "1",
     "nominal_moles": "2", "expires_on": "2026-11-01", "status": "available"},
    {"batch_no": "OLD", "species": "HCl", "purity": "1",
     "nominal_moles": "8", "expires_on": "2020-01-01", "status": "available"}
  ],
  "accounting_date": "2026-09-23"
}
```

字段约束：

| 字段 | 约束 |
| --- | --- |
| `demands` | 至少 1 条；`species` 为可解析化学式且同物种不得重复；`effective_moles` 为**非负十进制字符串**（JSON 数字会被拒绝） |
| `batches[].batch_no` | 非空字符串，单次请求内全局唯一 |
| `batches[].species` | 可解析化学式；等价写法（括号/空白/电荷风格）识别为同一物种 |
| `batches[].purity` | 纯度十进制字符串，必须满足 **0 < purity ≤ 1**（如 `"0.98"`、`"1.0"`） |
| `batches[].nominal_moles` | 可用名义摩尔数，**非负十进制字符串** |
| `batches[].expires_on` | 失效日期 `YYYY-MM-DD` |
| `batches[].status` | `"available"` 或 `"quarantined"` |
| `accounting_date` | 核算日期 `YYYY-MM-DD`，必填 |

`200` 响应（节选；上例中过期批 `OLD` 跳过，`A-1` 有效 0.75 耗尽后再从
`A-2` 领有效 0.75 = 名义 0.75）：

```json
{
  "as_of": "2026-09-23",
  "allocations": [
    {
      "species": "HCl",
      "normalized": "HCl",
      "requested_effective_moles": "1.5",
      "allocated_effective_moles": "1.5",
      "unmet_effective_moles": "0",
      "status": "fulfilled",
      "batches": [
        {"batch_no": "OLD", "expires_on": "2020-01-01", "status": "available",
         "purity": "1", "nominal_moles": "8", "effective_moles": "8",
         "eligible": false, "skipped": true, "skip_reason": "expired", "used": false,
         "allocated_nominal_moles": "0", "allocated_effective_moles": "0",
         "remaining_nominal_moles": "8", "remaining_effective_moles": "8"},
        {"batch_no": "A-1", "expires_on": "2026-10-01", "status": "available",
         "purity": "0.75", "nominal_moles": "1", "effective_moles": "0.75",
         "eligible": true, "skipped": false, "skip_reason": null, "used": true,
         "allocated_nominal_moles": "1", "allocated_effective_moles": "0.75",
         "remaining_nominal_moles": "0", "remaining_effective_moles": "0"},
        {"batch_no": "A-2", "expires_on": "2026-11-01", "status": "available",
         "purity": "1", "nominal_moles": "2", "effective_moles": "2",
         "eligible": true, "skipped": false, "skip_reason": null, "used": true,
         "allocated_nominal_moles": "0.75", "allocated_effective_moles": "0.75",
         "remaining_nominal_moles": "1.25", "remaining_effective_moles": "1.25"}
      ]
    }
  ],
  "unmatched_batches": [],
  "summary": {
    "total_requested_effective_moles": "1.5",
    "total_allocated_effective_moles": "1.5",
    "total_unmet_effective_moles": "0",
    "demands_total": 1,
    "fulfilled_demands": 1,
    "shortage_demands": 0,
    "batches_total": 3,
    "batches_used": 2,
    "batches_skipped": 1,
    "unmatched_batches": 0,
    "status": "fulfilled"
  }
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `allocations[]` | 按规范化化学式排序的各物种分配；`batches` 为该物种的**全部提交批次**（含跳过批次），按（失效日期, 批号）排序 |
| `allocations[].status` | `fulfilled`（`unmet_effective_moles="0"`）或 `shortage` |
| `batches[].effective_moles` | 该批有效摩尔 = 名义摩尔 × 纯度 |
| `batches[].allocated_nominal_moles` / `allocated_effective_moles` | 从该批领用的名义 / 有效摩尔（名义 = 有效 ÷ 纯度） |
| `batches[].remaining_nominal_moles` / `remaining_effective_moles` | 领用后余量；跳过批次原样保留 |
| `batches[].skip_reason` | `expired`（已过期）/ `quarantined`（隔离）/ `exhausted`（名义量为 0）/ `null`（正常批次）；既隔离又过期时报 `quarantined` |
| `unmatched_batches[]` | 物种不在 `demands` 中的批次，`skip_reason` 恒为 `unmatched_species`，按（规范化物种, 失效日期, 批号）排序 |
| `summary` | 跨物种需求/满足/未满足合计、需求与批次计数及总体 `status`（任一物种短缺即为 `shortage`） |

错误响应：

| 场景 | HTTP / code |
| --- | --- |
| 需求物种或批次物种化学式无法解析（`details.field` 指出 `demands.species` / `batches.species`，后者附 `batch_no`） | `400 FORMULA_PARSE_ERROR` |
| 需求物种重复（`details.duplicate_demand_species`）、批号重复（`details.duplicate_batch_nos`） | `422 INVALID_INVENTORY` |
| 纯度非法（非 (0,1]、NaN、无穷、非数字；`details.invalid_purity`） | `422 INVALID_INVENTORY` |
| 有效/名义摩尔非法（负数、NaN、非数字；`details.invalid_effective_moles` / `details.invalid_nominal_moles`） | `422 INVALID_INVENTORY` |
| 摩尔数/纯度传成 JSON 数字、日期格式错误、状态非 available/quarantined、列表为空或缺字段、空白字符串 | `422 INVALID_REQUEST` |

该接口为无状态计算（不读写幂等缓存、不持久化扣减结果）；同一请求重复调用
得到完全一致的分配。

### 2.9 溶液配制 `POST /api/v1/solutions/blend`

提交 **1～4 种母液**（标识、各溶质浓度、可用体积）、**目标终体积**与
**各溶质目标浓度**，服务端以精确有理数（全程不经过 `float`）求一组
**非负母液取用量**，使得：

* 每个目标溶质 j 满足摩尔守恒 Σ(母液浓度 × 取用量) = 目标浓度 × 终体积
  （浓度 mol/L × 体积 mL，两侧同为毫摩尔量级）；
* 每种母液取用量不超过其可用体积（`0 ≤ x_i ≤ available_volume`）；
* 母液总取用量不超过终体积，**余量用水补足**（`water_added = 终体积 − 总量 ≥ 0`）；
* 多解时先取**母液总量最小**，再按**输入顺序取用量字典序最小**
  （如两种母液完全等价时，排在前面的取用量为 0）。

请求体（所有数值均为**非负十进制字符串**，不要传 JSON 数字；`target_volume`
必须为正）：

```bash
curl -X POST http://localhost:8000/api/v1/solutions/blend \
  -H 'Content-Type: application/json' \
  -d '{
    "stocks": [
      {"id": "conc-NaCl", "available_volume": "100",
       "solutes": [{"solute": "NaCl", "concentration": "1.0"},
                   {"solute": "KCl",  "concentration": "0.5"}]},
      {"id": "conc-KCl", "available_volume": "100",
       "solutes": [{"solute": "NaCl", "concentration": "0.2"},
                   {"solute": "KCl",  "concentration": "1.0"}]}
    ],
    "target_volume": "100",
    "targets": [
      {"solute": "NaCl", "concentration": "0.4"},
      {"solute": "KCl",  "concentration": "0.6"}
    ]
  }'
```

字段约束：

| 字段 | 约束 |
| --- | --- |
| `stocks` | 1～4 条；`id` 为非空字符串且请求内唯一 |
| `stocks[].solutes` | 该母液的溶质浓度列表（空列表表示纯溶剂）；同一溶质不得重复；`concentration` 为**非负十进制字符串**（mol/L） |
| `stocks[].available_volume` | 可用体积，**非负十进制字符串**（mL） |
| `target_volume` | 目标终体积，**正十进制字符串**（mL，必须 > 0） |
| `targets` | 至少 1 条；同一溶质不得重复；`concentration` 为**非负十进制字符串**（mol/L）；母液带入但未列入目标的溶质不约束，仅在 `final_concentrations` 中报告 |

`200` 响应（上例的解为 `conc-NaCl` 取 280/9 mL、`conc-KCl` 取 400/9 mL，
无限小数按 **28 位有效数字**银行家舍入输出，规则同其它计量接口）：

```json
{
  "target_volume": "100",
  "total_stock_volume": "75.55555555555555555555555556",
  "water_added": "24.44444444444444444444444444",
  "stocks": [
    {"id": "conc-NaCl", "available_volume": "100",
     "volume_used": "31.11111111111111111111111111",
     "remaining_volume": "68.88888888888888888888888889"},
    {"id": "conc-KCl", "available_volume": "100",
     "volume_used": "44.44444444444444444444444444",
     "remaining_volume": "55.55555555555555555555555556"}
  ],
  "final_concentrations": [
    {"solute": "KCl", "target_concentration": "0.6",
     "final_concentration": "0.6", "target_met": true},
    {"solute": "NaCl", "target_concentration": "0.4",
     "final_concentration": "0.4", "target_met": true}
  ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `total_stock_volume` | 母液总取用量（mL），已在全部可行解中最小化 |
| `water_added` | 补水量 = 终体积 − 母液总取用量（mL，非负） |
| `stocks[]` | 各母液取用量与取用后余量，**与请求 stocks 顺序一致** |
| `final_concentrations[]` | 全部涉及溶质（目标 ∪ 母液带入）的最终浓度，按溶质标识排序；非目标溶质的 `target_concentration` / `target_met` 为 `null` |

错误响应：

| 场景 | HTTP / code |
| --- | --- |
| 无可行解：目标超出供应、可用体积或终体积预算不足、零目标与必需母液冲突等（`details.insufficient_solutes` 逐溶质给出需求/可供 mmol，仅单溶质供应不足时出现） | `422 NO_FEASIBLE_BLEND` |
| 母液标识重复（`details.duplicate_stock_ids`）、母液内溶质重复（`details.duplicate_stock_solutes`）、目标溶质重复（`details.duplicate_target_solutes`） | `422 INVALID_BLEND` |
| 浓度/体积非法：负数、NaN、非数字字符串（`details.invalid_concentrations` / `invalid_available_volumes` / `invalid_target_concentrations`）；终体积非正（`details.invalid_target_volume`） | `422 INVALID_BLEND` |
| 母液超过 4 种、列表为空、缺字段、数值传成 JSON 数字、空白字符串 | `422 INVALID_REQUEST` |

该接口为无状态计算（不读写幂等缓存）；同一请求重复调用得到完全一致的结果。

### 2.10 同位素名义质量峰分布 `POST /api/v1/isotopes/envelope`

提交一个**化学式**与**化学式中每种元素的一组同位素定义**（质量数 + 丰度），
服务端把分子中的每个原子看作一次**独立**的同位素抽样：

* 元素原子数沿用既有的递归下降解析器取得（支持括号、嵌套、下标与离子
  电荷；**电子物种 `e-` / `e+` 拒绝**）；
* 对含 n 个原子、k 种同位素（质量数 m₁..mₖ、丰度 p₁..pₖ）的元素，枚举
  全部计数向量 (c₁,…,cₖ)（Σcᵢ = n），按**多项分布**给出质量 Σcᵢ·mᵢ 与
  概率 n!/(c₁!…cₖ!)·p₁ᶜ¹…pₖᶜᵏ；再对各元素分布做卷积，**同质量组合
  在全部枚举完成之后合并**——不截断任何低概率（含恰好为 0）的组合；
* 全程使用精确有理数（丰度十进制字符串先无损转为 `Fraction`，不经过
  `float`），峰概率以**最简分数字符串**返回（`"3/8"`、`"1"`、`"0"`）；
* 返回 `peaks` 按**名义质量升序**排列，`most_probable` 为概率最高的峰，
  **概率并列时取质量较小者**。

请求约束：

| 约束 | 上限 |
| --- | --- |
| 元素种类 | **至多 3 种**（化学式中出现的每种元素都必须恰好给出一组，不能缺失、重复或给出不存在的元素） |
| 总原子数 | **不超过 30** |
| 每组同位素种数 | **至多 3 种**，质量数必须是两两**互异的正整数** |
| 丰度 | 均为**非负十进制字符串**（不要传 JSON 数字），同一组丰度**精确之和为 1**（允许某个丰度为 `"0"`，对应质量峰概率为 `"0"` 但仍返回） |

请求示例（HCl：两个元素都只有 1 个原子，四个质量峰；质量 36 与 37
概率并列 3/8，最高概率峰取质量较小者 36）：

```bash
curl -X POST http://localhost:8000/api/v1/isotopes/envelope \
  -H 'Content-Type: application/json' \
  -d '{
    "formula": "HCl",
    "elements": [
      {"element": "H", "isotopes": [
        {"mass_number": 1, "abundance": "0.5"},
        {"mass_number": 2, "abundance": "0.5"}]},
      {"element": "Cl", "isotopes": [
        {"mass_number": 35, "abundance": "0.75"},
        {"mass_number": 37, "abundance": "0.25"}]}
    ]
  }'
```

`200` 响应：

```json
{
  "input": "HCl",
  "normalized": "HCl",
  "total_atoms": 2,
  "elements": [
    {"element": "H", "atom_count": 1, "isotopes": [
      {"mass_number": 1, "abundance": "0.5"},
      {"mass_number": 2, "abundance": "0.5"}]},
    {"element": "Cl", "atom_count": 1, "isotopes": [
      {"mass_number": 35, "abundance": "0.75"},
      {"mass_number": 37, "abundance": "0.25"}]}
  ],
  "peaks": [
    {"mass": 36, "probability": "3/8"},
    {"mass": 37, "probability": "3/8"},
    {"mass": 38, "probability": "1/8"},
    {"mass": 39, "probability": "1/8"}
  ],
  "most_probable": {"mass": 36, "probability": "3/8"}
}
```

`Cl2`（丰度 0.75/0.25）展示多项系数与最简分数：两个原子都是 35Cl
概率 9/16；一个 35Cl 一个 37Cl 有两条抽样路径，合并为质量 72、
概率 3/8；都是 37Cl 为 1/16。零丰度同位素（如
`[{"mass_number":12,"abundance":"1"},{"mass_number":14,"abundance":"0"}]`）
对应的质量峰以 `"0"` 概率保留返回。

字段说明：

| 字段 | 含义 |
| --- | --- |
| `normalized` | 规范化后的化学式（括号与电荷归一） |
| `total_atoms` | 化学式总原子数（≤ 30） |
| `elements[].atom_count` | 由解析器取得的该元素原子数 |
| `peaks[]` | 全部名义质量峰，按质量升序，概率为最简分数字符串；同质量组合已合并，未截断低/零概率峰 |
| `most_probable` | 概率最高的峰；概率并列时取质量较小者 |

错误响应：

| 场景 | HTTP / code |
| --- | --- |
| 化学式无法解析、**电子物种**（`details.field="formula"`） | `400 FORMULA_PARSE_ERROR` |
| 元素组缺失（`details.missing_elements`）、重复（`duplicate_elements`）、给出化学式中不存在的元素（`unexpected_elements`）、元素符号非法（`invalid_elements`） | `422 INVALID_ISOTOPE` |
| 质量数重复（`details.duplicate_mass_numbers`） | `422 INVALID_ISOTOPE` |
| 丰度非法：负数、NaN、非数字（`details.invalid_abundances`） | `422 INVALID_ISOTOPE` |
| 丰度之和不精确等于 1（`details.abundance_sums`，给出实际的最简分数和） | `422 INVALID_ISOTOPE` |
| 元素种类超过 3 种（`details.element_count`）或总原子数超过 30（`details.total_atoms`） | `422 INVALID_ISOTOPE` |
| 请求体缺字段、列表为空或超长、质量数非正整数、丰度传成 JSON 数字、空白字符串等结构问题 | `422 INVALID_REQUEST` |

该接口为无状态计算（不读写幂等缓存）；同一请求重复调用得到完全一致的结果。

### 2.11 赫斯定律反应焓推导 `POST /api/v1/thermochemistry/hess`

提交一个**目标反应**与 **1～6 条已知反应**（每条含反应物、生成物与
`enthalpy_kj`——以 kJ 表示、**可带符号的有限十进制字符串**，放热为负、
吸热为正，与反应书写方向一致）。服务端对每个反应复用既有解析器与
**唯一配平**，按规范物种建立**生成物为正、反应物为负**的有符号计量向量，
然后以精确有理数（全程不经过 `float`）求解线性组合

    A·x = t（A 的列为已知反应向量，t 为目标反应向量）

倍数 `x` 的每个分量都是有理数，**负值表示反向使用、分数表示分数倍使用**。
推导反应焓 = Σ 倍数 × 已知反应焓（反向自动变号）。

请求体：

```json
{
  "target": {"reactants": ["C", "H2"], "products": ["CH4"]},
  "known_reactions": [
    {"reactants": ["C", "O2"], "products": ["CO2"], "enthalpy_kj": "-393.5"},
    {"reactants": ["H2", "O2"], "products": ["H2O"], "enthalpy_kj": "-571.6"},
    {"reactants": ["CH4", "O2"], "products": ["CO2", "H2O"], "enthalpy_kj": "-890.3"}
  ]
}
```

上例即 CH₄ 标准生成焓：`C + O2 -> CO2`（×1）、`2H2 + O2 -> 2H2O`（×1）、
`CH4 + 2O2 -> CO2 + 2H2O`（**×−1 反向**）相加消去 O₂/CO₂/H₂O，得到
`C + 2H2 -> CH4`，ΔH = −393.5 + (−571.6) + 890.3 = **−74.8 kJ**。

`200` 响应：

```json
{
  "target": {
    "reactants": ["C", "H2"],
    "products": ["CH4"],
    "equation": "C + 2H2 -> CH4"
  },
  "known_reactions": [
    {"index": 0, "reactants": ["C", "O2"], "products": ["CO2"],
     "equation": "C + O2 -> CO2", "enthalpy_kj": "-393.5",
     "multiplier": "1", "multiplier_used": "-393.5"},
    {"index": 1, "reactants": ["H2", "O2"], "products": ["H2O"],
     "equation": "2H2 + O2 -> 2H2O", "enthalpy_kj": "-571.6",
     "multiplier": "1", "multiplier_used": "-571.6"},
    {"index": 2, "reactants": ["CH4", "O2"], "products": ["CO2", "H2O"],
     "equation": "CH4 + 2O2 -> CO2 + 2H2O", "enthalpy_kj": "-890.3",
     "multiplier": "-1", "multiplier_used": "890.3"}
  ],
  "derived_enthalpy_kj": "-74.8",
  "net_coefficients": [
    {"normalized": "C", "coefficient": "-1", "side": "reactant"},
    {"normalized": "CH4", "coefficient": "1", "side": "product"},
    {"normalized": "CO2", "coefficient": "0", "side": "cancelled"},
    {"normalized": "H2", "coefficient": "-2", "side": "reactant"},
    {"normalized": "H2O", "coefficient": "0", "side": "cancelled"},
    {"normalized": "O2", "coefficient": "0", "side": "cancelled"}
  ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `target` | 规范化并唯一配平后的目标反应 |
| `known_reactions[].index` | 从 0 开始的已知反应序号（与请求顺序一致） |
| `known_reactions[].multiplier` | 该反应在组合中的倍数，**最简分数字符串**（`"1"`、`"-1"`、`"1/2"`、`"3/2"`）；负为反向、分数为分数倍 |
| `known_reactions[].multiplier_used` | 该反应的实际焓贡献 = 倍数 × 反应焓（kJ 十进制字符串） |
| `derived_enthalpy_kj` | 推导出的目标反应焓 ΔH（kJ 十进制字符串），负号表示放热 |
| `net_coefficients[]` | 全部涉及物种（目标 ∪ 已知反应）**合并抵消后的逐物种净系数**，最简分数字符串，按规范化化学式排序 |
| `net_coefficients[].side` | `product`（净系数为正）/ `reactant`（为负）/ `cancelled`（为 0，中间物种完全抵消） |

倍数以最简分数给出（如由 `2C + O2 -> 2CO`（−221）与
`2CO + O2 -> 2CO2`（−566）各取 **1/2** 组合出 `C + O2 -> CO2` 时，
`multiplier` 为 `"1/2"`、`"1/2"`，推导 ΔH = −393.5 kJ）；反应焓与贡献
按十进制字符串输出，无限小数按 28 位有效数字银行家舍入，规则同其它
计量接口。

**可推导性与歧义判定**：对增广矩阵 `[A | t]` 做有理数 RREF——

* 出现恒等约束 `0 = 非零`（目标含已知反应从未涉及的物种，或计量关系
  互相矛盾）时返回 `422 HESS_UNDERIVABLE`；目标独有物种在
  `details.target_only_species` 中列出；
* 组合不唯一（已知反应线性相关）时，沿每个齐次方向检查反应焓的变化：
  存在使总焓改变的方向 -> `422 HESS_AMBIGUOUS`（已知数据不自洽）；
  所有方向的焓变化都为 0 -> 各组合导出的反应焓相同，**正常返回**
  （取自由变量为 0 的一个特解，如自洽的闭合反应循环）。

错误响应：

| 场景 | HTTP / code |
| --- | --- |
| 目标或某条已知反应化学式无法解析（`details.where` 为 `target` / `known_reactions`，后者附 `details.index`） | `400 FORMULA_PARSE_ERROR` |
| 目标或某条已知反应同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 目标或某条已知反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |
| 某条反应焓非法：NaN、Infinity、非数字、JSON 数字（`details.invalid_enthalpies` 附序号与原值）；已知反应为 0 条或超过 6 条、缺字段等结构问题 | `422 INVALID_REQUEST` |
| 目标反应不在已知反应向量的张成空间内 | `422 HESS_UNDERIVABLE` |
| 多组组合导出不同的目标反应焓 | `422 HESS_AMBIGUOUS` |

该接口为无状态计算（不读写幂等缓存，响应不含 `reaction_id`）；同一请求
重复调用得到完全一致的结果。离子反应同样适用（电子 `e-` 作为物种参与
配平与抵消，如两个还原半反应组合为置换反应时 `e-` 净系数为 `"0"`）。

### 2.12 理想溶液平衡组成 `POST /api/v1/equilibrium/solve`

对**单个可逆反应**（理想溶液）求平衡组成。服务端先复用与 `/balance` 相同
的解析与**唯一配平**（及幂等缓存），设反应进度为 ξ（mol，按配平式书写
方向为正），以配平系数确定所有物种摩尔数都非负的双向物料区间

    ξ_min = −min_生成物( n0_i ÷ b_i ) ≤ 0
    ξ_max =  min_反应物( n0_i ÷ a_i ) ≥ 0

再在区间内按各物种**浓度的计量次方**构造反应商

    Qc(ξ) = Π_生成物 (n_i/V)^b_i ÷ Π_反应物 (n_i/V)^a_i

Qc 在开区间 (ξ_min, ξ_max) 上关于 ξ 严格单调递增（端点处取 0 与 +∞），
因此 **Qc = Kc 有且仅有一个内点解**（Kc 严格为正）：

* Qc(0) < Kc：`direction="forward"`，ξ > 0，反应正向进行；
* Qc(0) > Kc：`direction="reverse"`，**ξ < 0（负进度）**，反应逆向进行；
* Qc(0) = Kc：`direction="at_equilibrium"`，ξ = 0。

求解全程使用精确 :class:`fractions.Fraction` 二分（符号判定无浮点抵消），
区间宽度同时满足 **< 1e-9 mol**（中点绝对误差 < 5e-10 mol）与相对宽度
< 区间尺度 × 1e-32（任意数量级输入都保留 32 位以上有效数字），故返回
**进度的绝对误差上界为 1e-9 mol**。

请求体在 `reactants` / `products` 之外增加 `amounts`、`volume_l` 与 `kc`：

| 字段 | 约束 |
| --- | --- |
| `amounts` | **每个参与物种（反应物 ∪ 生成物）恰好一条记录**，等价写法可识别；不能缺失、重复或出现未知物种；各 `moles` 必须是**严格为正**的十进制字符串（不接受 `"0"`、负数或 JSON 数字） |
| `volume_l` | 溶液体积（L），**正**十进制字符串 |
| `kc` | 浓度平衡常数 Kc，**正**十进制字符串（科学计数法如 `"1.8e-5"` 亦可） |
| 反应两侧 | 同一规范化物种**不得同时出现在两侧**（可逆反应要求每个物种只属于一侧） |

`PCl5 ⇌ PCl3 + Cl2`、1 L、初始 PCl5=1、PCl3=Cl2=0.01、Kc=0.8
（平衡 ξ ≈ 0.57387，PCl5 余约 0.42613、PCl3=Cl2≈0.58387）：

```bash
curl -X POST http://localhost:8000/api/v1/equilibrium/solve \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["PCl5"],
    "products": ["PCl3", "Cl2"],
    "amounts": [
      {"species": "PCl5", "moles": "1"},
      {"species": "PCl3", "moles": "0.01"},
      {"species": "Cl2",  "moles": "0.01"}
    ],
    "volume_l": "1",
    "kc": "0.8"
  }'
```

`200` 响应（十进制值均为字符串；物种按规范化顺序排列）：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["PCl5"],
    "products": ["Cl2", "PCl3"],
    "equation": "PCl5 -> Cl2 + PCl3"
  },
  "volume_l": "1",
  "kc": "0.8",
  "direction": "forward",
  "reaction_extent": {"value": "0.5738699100999075…", "unit": "mol"},
  "extent_bounds": {"minimum": "-0.01", "maximum": "1"},
  "initial_qc": "0.0001",
  "equilibrium_qc": "0.8",
  "reactants": [
    {"input": "PCl5", "normalized": "PCl5", "coefficient": 1,
     "initial_moles": "1", "equilibrium_moles": "0.4261300899000925…"}
  ],
  "products": [
    {"input": "Cl2", "normalized": "Cl2", "coefficient": 1,
     "initial_moles": "0.01", "equilibrium_moles": "0.5838699100999075…"},
    {"input": "PCl3", "normalized": "PCl3", "coefficient": 1,
     "initial_moles": "0.01", "equilibrium_moles": "0.5838699100999075…"}
  ]
}
```

**负进度示例**：`2O3 ⇌ 3O2`、初始 O3=2、O2=3、Kc=0.375（初始
Qc=6.75>Kc），平衡 ξ 为 **−0.5**（O3=3、O2=1.5，Qc=0.375）。
**初始即平衡示例**：`2NO2 ⇌ N2O4`，初始 NO2=2、N2O4=1、V=1、
Kc=0.25 时 ξ=0、`direction="at_equilibrium"`。体积通过浓度项
（Δn≠0 时带 V^Δn 因子）进入 Qc，不是仅按摩尔数计算。

字段说明：

| 字段 | 含义 |
| --- | --- |
| `direction` | `forward`（ξ>0）/ `reverse`（ξ<0，初始 Qc>Kc）/ `at_equilibrium`（ξ=0） |
| `reaction_extent.value` | 平衡反应进度 ξ（mol），绝对误差上界 1e-9 mol；负为逆向 |
| `extent_bounds.minimum` / `maximum` | 非负物料区间 ξ_min / ξ_max（mol），由初始量与配平系数确定 |
| `initial_qc` / `equilibrium_qc` | 初始浓度反应商 Qc(0) 与平衡反应商（后者应等于 Kc） |
| `reactants[]` / `products[]` | 各物种配平系数、初始摩尔数与平衡摩尔数（= 初始量 ∓ 系数·ξ，均非负） |

错误沿用统一结构化格式：

| 场景 | HTTP / code |
| --- | --- |
| 反应式或初始量物种化学式无法解析（`details.field=amounts.species`） | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 同一规范化物种同时出现在反应两侧（`details.cross_side_duplicate`） | `422 INVALID_REQUEST` |
| 初始量缺失、重复、含未知物种（`details.missing/duplicate/unknown`） | `422 INVALID_REQUEST` |
| 初始摩尔数为 `"0"`/负数/NaN/非数字、`volume_l` 或 `kc` 非正或非数字（`details.invalid_moles/invalid_volume/invalid_kc`） | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

`reaction_id` 与 `/balance`、`/analyze` 及其它计量接口共享同一规范化
签名；平衡求解本身不单独缓存，只复用配平缓存。

### 2.13 元素分析反推分子式 `POST /api/v1/formulas/infer`

提交 **C、H、N、O 四种元素各一条**元素分析测量记录（实测质量百分比
`percent` 与绝对容差 `tolerance`，均为 **0～100 的十进制字符串**）与
**名义摩尔质量的正整数闭区间** `molar_mass`，服务端按名义原子量
**C=12、H=1、N=14、O=16** 枚举整数原子计数

    n_C, n_H, n_N, n_O ∈ [0, 30]，  n_C + n_H + n_N + n_O > 0

要求每个候选**同时**满足：

* 摩尔质量（正整数，含端点）`minimum ≤ 12n_C + n_H + 14n_N + 16n_O ≤ maximum`；
* 每种元素的计算质量百分比落入**闭**容差带：
  `|100·a·n/M − percent| ≤ tolerance`（恰在边界上也算匹配）。

四种元素的测量值**不要求恰好合计 100**（允许实验归一化误差）；各元素
容差也是绝对值（百分数单位），不随其它元素缩放。全部输入先无损转为
`Fraction`，热循环只做整数乘除与交叉相乘（不经过 `float`），因此
恰落在容差边界上的候选不会被浮点误差误判。

匹配结果排序键依次为：**总百分比绝对偏差**
`Σ_e |计算百分比 − percent|`（精确有理数）、**摩尔质量**、**C/H/N/O
计数**，全部升序；只返回前 **20** 个，并以 `total` 给出截断前的匹配
总数、`truncated` 标识是否截取。**无任何匹配时返回 `200` 且
`matches=[]`、`total=0`**（无解放到正常响应中，不作为错误）。

请求体（`elements` 必须恰好 4 条，元素不得缺失或重复）：

```json
{
  "elements": [
    {"element": "C", "percent": "40.00", "tolerance": "0.05"},
    {"element": "H", "percent": "6.71",  "tolerance": "0.05"},
    {"element": "N", "percent": "0",     "tolerance": "0.05"},
    {"element": "O", "percent": "53.29", "tolerance": "0.05"}
  ],
  "molar_mass": {"minimum": 170, "maximum": 190}
}
```

上例唯一匹配葡萄糖 `C6H12O6`（M=180；计算百分比 C 40%、H 6.6̅%、
O 53.3̅%，均落入容差），`200` 响应：

```json
{
  "elements": [
    {"element": "C", "percent": "40.00", "tolerance": "0.05"},
    {"element": "H", "percent": "6.71",  "tolerance": "0.05"},
    {"element": "N", "percent": "0",     "tolerance": "0.05"},
    {"element": "O", "percent": "53.29", "tolerance": "0.05"}
  ],
  "molar_mass": {"minimum": 170, "maximum": 190},
  "count": 1,
  "total": 1,
  "truncated": false,
  "matches": [
    {
      "formula": "C6H12O6",
      "molar_mass": 180,
      "counts": {"C": 6, "H": 12, "N": 0, "O": 6},
      "mass_percent": {
        "C": "40",
        "H": "6.666666666666666666666666667",
        "N": "0",
        "O": "53.33333333333333333333333333"
      },
      "total_abs_deviation": "0.08666666666666666666666666667"
    }
  ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `elements` / `molar_mass` | 请求的测量记录与质量区间（原样回显，字符串保留客户端写法） |
| `matches[].formula` | 分子式（C/H/N/O 顺序；计数 0 的元素省略、计数 1 省略下标） |
| `matches[].counts` | 四种元素原子计数（四键齐全；逐项 0～30、合计大于零） |
| `matches[].molar_mass` | 名义摩尔质量 = 12n_C + n_H + 14n_N + 16n_O（正整数） |
| `matches[].mass_percent` | 各项**计算**质量百分比（%，十进制字符串；无限小数按 28 位有效数字银行家舍入） |
| `matches[].total_abs_deviation` | 总百分比绝对偏差，第一排序键 |
| `count` / `total` / `truncated` | 本次返回条数 / 截断前匹配总数 / 是否超过 20 已截取 |

误差说明：`mass_percent` 等十进制字符串仅在输出时按 28 位有效数字
舍入，匹配判定与排序都基于精确有理数；对同一请求重复调用得到完全
一致的结果。

错误响应：

| 场景 | HTTP / code |
| --- | --- |
| 元素缺失（`details.missing`）或重复（`details.duplicate`） | `422 INVALID_FORMULA_INFER` |
| `percent` / `tolerance` 超出 0～100、为负数、NaN、无穷或非数字字符串（`details.invalid_percent` / `invalid_tolerance` 附元素与原值） | `422 INVALID_FORMULA_INFER` |
| `molar_mass.minimum > maximum`（`details.invalid_molar_mass_range`） | `422 INVALID_FORMULA_INFER` |
| 请求体缺字段、`elements` 不是 4 条、元素符号非 C/H/N/O、百分比传成 JSON 数字、摩尔质量端点非正整数或为布尔值等结构问题 | `422 INVALID_REQUEST` |

该接口为无状态计算（不读写幂等缓存，响应不含 `reaction_id`）；同一
请求重复调用得到完全一致的结果。

### 2.14 幂等与查询

* 幂等键与书写顺序、括号风格、空白以及电荷的等价写法无关，只取决于每一侧
  物种的 **（元素组成多重集, 电荷）**。
* 因此下面两个请求得到 **同一个 `reaction_id`**，第二次 `cached=true`：

```bash
curl -X POST .../api/v1/balance -d '{"reactants":["Ca(OH)2","HCl"],"products":["CaCl2","H2O"]}'
curl -X POST .../api/v1/balance -d '{"reactants":["Ca [OH]2"," HCl"],"products":["CaCl2","H2 O"]}'
```

* 用 `GET /api/v1/reactions/{reaction_id}` 可随时取回缓存结果（LRU，默认 1 万条）。

---

## 3. 化学式书写约定

| 特性 | 写法示例 |
| --- | --- |
| 元素（118 种） | `H`、`Fe`、`O`、`Uuo`（必须大写开头） |
| 下标 | `H2O`、`Al2(SO4)3` |
| 三种括号（等价，可嵌套） | `Ca(OH)2`、`K4[Fe(CN)6]`、`Mg3{Fe(CN)6}2` |
| 显式电荷 | `SO4^2-`、`Fe^3+`、`NH4^+`、`O^{2-}`、`[Cu(NH3)4]^2+` |
| 纯文本电荷 | `Fe3+`、`O2-`、`NH4+`、`MnO4-`、`OH-` |
| 电荷前空格消歧 | `SO4 2-` 等价于 `SO4^2-` |
| 电子 / 正电子 | `e-`、`e`（默认 −1）、`e+` |
| 空白 | 化学式内部空白一律忽略 |

**末尾 `数字±`（无 `^`）消歧规则**（沿用常见化学软件约定）：

* 物种只含 **一个原子/基团** 时，数字是电荷：`Fe3+` → Fe³⁺、`O2-` → O²⁻；
* 含 **多个原子** 时，数字是下标、符号单独表示 ±1：`NH4+` → NH₄⁺、`MnO4-` → MnO₄⁻；
* 数字前有空格时一定是电荷：`SO4 2-` → SO₄²⁻；
* 任何不确定的场景都建议写显式 `^`，无歧义。

---

## 4. 配平原理

设所有物种系数向量为 `x`（未知）。对每种元素写一行守恒方程，再加一行
总电荷守恒；反应物一侧的列取正，生成物一侧的列取负，得到齐次方程组
`A x = 0`。

对整数矩阵的每个元素套上 `Fraction` 做 RREF：

* **零空间维数 = 1**：基向量所有分量非零且同号 → 整体取正，按分母最小公倍数
  放大、再除以最大公约数，得到最小互质正整数系数；
* **维数 = 0**：`IMPOSSIBLE_BALANCE`（元素/电荷守恒本身无解）；
* **维数 > 1**：`/balance` 返回 `NON_UNIQUE_BALANCE`（自由度 = 维数 − 1，
  常见于缺少电子、H₂O、H⁺/OH⁻ 等物种）；
* 维数为 1 但有零分量 → `ZERO_COEFFICIENT`；符号不一致 → `INVALID_STOICHIOMETRY`。

**多变体枚举（`/balance/variants`）**：RREF 后第 k 个基向量在第 k 个自由列
上为 1、其余自由列为 0，因此直接令自由分量取 `1..max_coefficient` 的整数，
主元分量按 `x_i = Σ b_k[i]·t_k` 用精确有理数计算，保留全部主元分量也是
`1..max_coefficient` 整数且整体 gcd 为 1 的向量。每个主元行预先通分（分母
LCM），热循环只剩整数乘加；自由变量最多 3 个、系数上限 30，枚举规模至多
`30³ = 27000` 组，通常在 0.2 秒内完成。

全程使用有理数，不会出现浮点误差，大系数反应（如
`10 K4Fe(CN)6 + 122 KMnO4 + 299 H2SO4 -> …`）也能精确求解。

---

## 5. 错误码

错误响应统一形如 `{"error": {"code": ..., "message": ..., "details": {...}}, "reaction_id": ..., "cached": ...}`。

| code | HTTP | 含义 |
| --- | --- | --- |
| `FORMULA_PARSE_ERROR` | 400 | 化学式无法解析（未知元素、括号不匹配、电荷非法等） |
| `DUPLICATE_SPECIES` | 409 | 同一侧出现重复物种 |
| `INVALID_REQUEST` | 422 | 请求体缺失/为空/物种过多（>100）；`/stoichiometry/batch` 投料缺失、重复、含未知物种或摩尔数非法（含 JSON 数字）；`/stoichiometry/plan` 目标重复、含非生成物物种或摩尔数非法；`/stoichiometry/yield` 投料或实际产量缺失/重复/含未知物种或任一侧摩尔数非法（实际产量侧问题在 `details.actual_yields.*`）；`/stoichiometry/sequence` 步骤为空、库存摩尔数非法（结构化错误的 `details.step` 另标出出错步骤）；`/equilibrium/solve` 同一物种同时出现在反应两侧、初始量缺失/重复/含未知物种、初始摩尔数或体积/Kc 非正十进制字符串；`/formulas/infer` 请求体结构问题（`elements` 不是 4 条、元素非 C/H/N/O、百分比/容差传成 JSON 数字或为空白、摩尔质量端点非正整数等）；`/inventory/allocate` 请求体结构问题（缺字段、空列表、日期格式错误、状态枚举非法、数字代替十进制字符串等） |
| `IMPOSSIBLE_BALANCE` | 422 | 守恒方程组无解（含电荷无法守恒） |
| `NON_UNIQUE_BALANCE` | 422 | 零空间维数 > 1，配平不唯一（返回 `nullity` 与自由度） |
| `ZERO_COEFFICIENT` | 422 | 唯一解要求某些物种系数为 0 |
| `INVALID_STOICHIOMETRY` | 422 | 唯一解符号不一致，物种可能需要换边 |
| `VARIANT_SPACE_TOO_LARGE` | 422 | `/balance/variants` 零空间维数（自由变量数）> 3，枚举空间过大（返回 `nullity`） |
| `INVALID_INVENTORY` | 422 | `/inventory/allocate` 需求物种或批号重复、纯度不在 (0,1]、有效/名义摩尔数非法（含 JSON 数字）或批次状态未知，问题明细在 `details` 中聚合返回 |
| `INVALID_BLEND` | 422 | `/solutions/blend` 母液标识或溶质重复、目标溶质重复、浓度/体积非法（含 JSON 数字）或目标终体积非正，问题明细在 `details` 中聚合返回 |
| `NO_FEASIBLE_BLEND` | 422 | `/solutions/blend` 无可行解：在可用体积与终体积约束下不存在满足全部目标浓度的非负取用量（`details.insufficient_solutes` 给出单溶质供应诊断） |
| `INVALID_ISOTOPE` | 422 | `/isotopes/envelope` 同位素输入非法：元素组缺失/重复/多余/非法、质量数重复、丰度非法或精确之和不为 1、超过元素/原子数上限，问题明细在 `details` 中聚合返回 |
| `HESS_UNDERIVABLE` | 422 | `/thermochemistry/hess` 目标反应不在已知反应计量向量的张成空间内（含目标独有物种，见 `details.target_only_species`），任何有符号/分数倍组合都无法复现目标反应 |
| `HESS_AMBIGUOUS` | 422 | `/thermochemistry/hess` 存在多组线性组合可复现目标反应的计量向量，但按各组已知反应焓导出的目标反应焓不一致（已知数据不自洽） |
| `INVALID_FORMULA_INFER` | 422 | `/formulas/infer` C/H/N/O 元素缺失/重复、实测百分比或绝对容差超出 0～100（含负数/NaN/非数字）或摩尔质量区间倒置，问题明细在 `details` 中聚合返回 |
| `NOT_FOUND` | 404 | `reaction_id` 不在缓存中 |

不唯一示例：`C + O2 -> CO + CO2` 缺少约束，返回

```json
{
  "error": {
    "code": "NON_UNIQUE_BALANCE",
    "message": "配平不唯一：零空间维数为 2，存在 1 个自由度……",
    "details": {"nullity": 2, "degrees_of_freedom": 1, "elements": ["C", "O"]}
  }
}
```

---

## 6. 项目结构

```
app/
  parser.py         递归下降化学式解析（元素/括号/电荷/电子）
  balancer.py       守恒矩阵、有理数 RREF、零空间、最小正整数、多变体枚举、守恒视图
  cache.py          规范化签名（SHA-256）与 LRU 幂等缓存
  stoichiometry.py  投料量 -> 反应进度（正向批量计量）、目标产量 -> 最小反应进度/理论投料（反向规划）、投料+实际产量 -> 理论/实际/收率百分比（收率计算）与多步连续反应库存推演，Decimal/Fraction 精确计算
  equilibrium.py    理想溶液平衡组成：配平系数确定双向非负物料区间，Qc(ξ) 精确有理二分求 Qc=Kc 唯一内点（支持负进度），Fraction 精确计算（进度绝对误差 < 1e-9 mol）
  allocation.py     试剂批次 FEFO 领用分配：纯度换算有效物质、跨批拆分、过期/隔离跳过、未满足需求与汇总，Fraction 精确计算
  blending.py       溶液配制：目标溶质摩尔守恒 + 可用体积/终体积约束的精确线性规划（顶点枚举），母液总量最小、再按输入顺序字典序最小，Fraction 精确计算
  isotopes.py       同位素名义质量峰分布：逐元素多项分布枚举全部计数向量（不截断零概率组合）、同质量合并、跨元素卷积，Fraction 精确计算
  formula_infer.py  元素分析反推分子式：C/H/N/O 计数 0～30 区间剪枝 + 容差闭区间整数交叉相乘枚举，按总偏差/摩尔质量/计数堆排前 20，Fraction 精确计算
  hess.py           赫斯定律反应焓推导：规范物种有符号计量向量、增广矩阵有理数 RREF 判定张成/歧义、有符号（反向/分数倍）组合求解，Fraction 精确计算
  service.py        业务编排与响应组装
  models.py         Pydantic 请求/响应模型
  main.py           FastAPI 路由与异常处理
tests/              解析器、配平核心、HTTP API、批量计量、反向投料规划、实际收率、连续反应推演、试剂批次领用分配、溶液配制、同位素名义质量峰分布、元素分析反推分子式、赫斯定律反应焓推导、理想溶液平衡组成测试
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```

## 7. 已验证的反应示例

```
2H2 + O2 -> 2H2O
4Fe + 3O2 -> 2Fe2O3
C3H8 + 5O2 -> 3CO2 + 4H2O
2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 5Cl2 + 8H2O
3Ca(OH)2 + 2H3PO4 -> Ca3(PO4)2 + 6H2O
(NH4)2Cr2O7 -> Cr2O3 + N2 + 4H2O
Fe^2+ -> Fe^3+ + e^-
MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
3Cu + 2NO3^- + 8H^+ -> 3Cu^2+ + 2NO + 4H2O
10K4Fe(CN)6 + 122KMnO4 + 299H2SO4
  -> 162KHSO4 + 5Fe2(SO4)3 + 122MnSO4 + 60HNO3 + 60CO2 + 188H2O
```
