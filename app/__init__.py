"""化学反应计量平衡服务。"""

import sys

# 输入校验允许数量级达 1e±10000 的十进制字符串（见
# app.stoichiometry.MAX_ABS_EXPONENT）；精确有理数在渲染、Decimal 转换与
# 极端 Kc 平衡进度（非负物料区间端点残余量可小至 ~1e-10010）下需要
# 万位以上的 int<->str 转换。默认 4300 位上限会让这类合法输入在输出
# 渲染时崩溃，这里放宽到 100 万位（仍保留拒绝恶意超长整数的防线）。
sys.set_int_max_str_digits(1_000_000)

__version__ = "1.5.0"
