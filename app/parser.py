"""化学式解析。

支持的文法（递归下降）::

    formula  := seq [charge]
    seq      := ( element [number] | "(" seq ")" [number] | "[" seq "]" [number]
                  | "{" seq "}" [number] )*
    charge   := ("^" ( sign [digits] | digits sign )
                     | sign [digits]
                     | digits sign )
    element  := 大写字母 [小写字母]

约定与扩展：

* 电荷只允许出现在最外层（如 ``[SO4]^2-``、``OH-``、``Fe3+``）；
* 电子写作 ``e-`` / ``e^-``，裸写 ``e`` 也视为电子（电荷 -1）；``e+`` 为正电子；
* 三种括号等价，可嵌套，可带下标乘数（``Ca(OH)2``、``Mg3[Fe(CN)6]2``）；
* 化学式内部的空白会被忽略；电荷前的空白可用来消歧（``SO4 2-``）。

末尾 ``数字+正负号``（无 ``^``）的消歧约定：

* 物种仅含一个原子/基团时视为电荷：``Fe3+`` → Fe³⁺、``O2-`` → O²⁻；
* 含多个原子时数字是下标、符号单独表示 ±1：``NH4+`` → NH₄⁺、``MnO4-`` → MnO₄⁻；
* 电荷前有空格时数字一定是电荷：``SO4 2-`` → SO₄²⁻；
* 任何无歧义场景都可显式写 ``^``：``SO4^2-``、``[Cu(NH3)4]^2+``、``O^{2-}``。
"""

from __future__ import annotations

from dataclasses import dataclass, field

# 118 种元素符号
ELEMENTS = frozenset(
    """H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni
    Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
    Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au
    Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf
    Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og""".split()
)

_OPEN = {"(": ")", "[": "]", "{": "}"}
_CLOSE = set(_OPEN.values())
_SIGN = {"+": 1, "-": -1}

# 节点表示（不可变元组）:
#   ("el", symbol, count)
#   ("grp", tuple[node, ...], count)
Node = tuple


class FormulaError(ValueError):
    """化学式无法解析时抛出。"""


@dataclass(frozen=True)
class ParsedFormula:
    """解析后的化学式。"""

    text: str                      # 去除空白后的原始文本
    normalized: str                # 规范化文本（圆括号 + ^q+ 电荷）
    composition: dict[str, int]    # 元素 -> 原子数（电子为空）
    charge: int                    # 总电荷
    is_electron: bool = False
    nodes: tuple[Node, ...] = field(default=())

    @property
    def key(self) -> tuple:
        """按化学组成判等的规范键：(原子数元组, 电荷)。"""
        return tuple(sorted(self.composition.items())), self.charge


def parse_formula(text: str) -> ParsedFormula:
    """解析单个化学式，失败时抛出 :class:`FormulaError`。"""
    if text is None:
        raise FormulaError("化学式不能为空")
    stripped = text.strip()
    if not stripped:
        raise FormulaError("化学式不能为空")
    if len(stripped) > 500:
        raise FormulaError("化学式过长（最多 500 个字符）")

    parser = _Parser(stripped)
    nodes, electron, charge = parser.parse_top()

    if not nodes and not electron:
        raise FormulaError("化学式不能为空")

    composition: dict[str, int] = {}
    for node in nodes:
        _accumulate(node, 1, composition)

    if electron:
        # 电子只能单独出现
        if nodes:
            raise FormulaError("电子 e- 只能单独出现，不能与其他物种组合")
        if abs(charge) != 1:
            raise FormulaError("电子的电荷必须为 ±1（写作 e- 或 e+）")
        return ParsedFormula(
            text=stripped,
            normalized=_render_charge("e", charge),
            composition={},
            charge=charge,
            is_electron=True,
        )

    normalized = _render_nodes(nodes) + _render_charge("", charge)
    return ParsedFormula(
        text=stripped,
        normalized=normalized,
        composition=composition,
        charge=charge,
        nodes=nodes,
    )


class _Parser:
    def __init__(self, text: str) -> None:
        self.s = text
        self.i = 0
        # 括号栈：栈项为 True 表示已看到该组内电荷
        self._stack: list[bool] = []

    def _skip_ws(self) -> None:
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    # ------------------------------------------------------------------ 顶层
    def parse_top(self) -> tuple[tuple[Node, ...], bool, int]:
        nodes, electron, charge_seen = self._parse_seq(None)
        self._skip_ws()
        if charge_seen:
            charge = self._parse_charge_after_prefix()
        elif self.i < len(self.s):
            ch = self.s[self.i]
            if ch in _CLOSE:
                raise FormulaError(f"多余的闭合括号 {ch!r}（位置 {self.i}）")
            if ch.isdigit():
                charge = self._parse_charge_with_leading_digit()
            elif ch in _SIGN:
                charge = self._parse_charge_after_prefix()
            else:
                raise FormulaError(f"无法识别的字符 {ch!r}（位置 {self.i}）")
        else:
            charge = 0
            if electron:
                # 裸 e 默认表示电子（电荷 -1）
                charge = -1
        return tuple(nodes), electron, charge

    # ------------------------------------------------------------------ 序列
    def _parse_seq(self, close: str | None) -> tuple[list[Node], bool, bool]:
        """返回 (节点, 是否电子, 是否已遇到电荷前缀 ^)。"""
        nodes: list[Node] = []
        electron = False
        charge_seen = False
        while True:
            self._skip_ws()
            if self.i >= len(self.s):
                break
            ch = self.s[self.i]
            if close is not None and ch == close:
                self.i += 1
                return nodes, electron, charge_seen
            if ch in _CLOSE:
                raise FormulaError(f"多余的闭合括号 {ch!r}（位置 {self.i}）")
            if ch in _OPEN:
                self._stack.append(False)
                self.i += 1
                inner, inner_electron, inner_charge = self._parse_seq(_OPEN[ch])
                if not inner and not inner_electron:
                    raise FormulaError(f"空括号（位置 {self.i - 1}）")
                if inner_electron:
                    raise FormulaError("电子 e- 不能放在括号内")
                if inner_charge:
                    raise FormulaError("电荷只能写在整个化学式的最外层")
                self._stack.pop()
                count = self._parse_count(single_token=len(nodes) == 0)
                nodes.append(("grp", tuple(inner), count))
            elif ch == "^":
                if self._stack:
                    raise FormulaError("电荷只能写在整个化学式的最外层")
                charge_seen = True
                break
            elif ch.isdigit():
                if self._stack:
                    raise FormulaError(
                        f"意外的数字（位置 {self.i}）：下标必须紧跟元素或右括号"
                    )
                # 顶层裸数字只可能是 "2-" / "2+" 形式的电荷开头
                break
            elif ch in _SIGN:
                if self._stack:
                    raise FormulaError("电荷只能写在整个化学式的最外层")
                break
            elif ch.isupper():
                symbol, start = self._parse_element()
                count = self._parse_count(single_token=len(nodes) == 0)
                nodes.append(("el", symbol, count))
            elif ch == "e":
                if electron or nodes:
                    raise FormulaError("电子 e- 只能单独出现，不能与其他物种组合")
                electron = True
                self.i += 1
                self._parse_count_for_electron()
            else:
                raise FormulaError(f"无法识别的字符 {ch!r}（位置 {self.i}）")

        if close is not None:
            self._skip_ws()
            if self.i >= len(self.s):
                raise FormulaError(f"缺少与开括号匹配的 {close!r}")
        return nodes, electron, charge_seen

    def _parse_element(self) -> tuple[str, int]:
        start = self.i
        symbol = self.s[self.i]
        self.i += 1
        if self.i < len(self.s) and self.s[self.i].islower():
            symbol += self.s[self.i]
            self.i += 1
        if symbol not in ELEMENTS:
            raise FormulaError(f"未知元素符号 {symbol!r}（位置 {start}）")
        return symbol, start

    def _parse_count(self, *, single_token: bool) -> int:
        """解析下标数字。

        顶层对 ``数字+正负号`` 结尾做电荷前瞻：

        * 数字前是空白（如 ``SO4 2-``）：数字一定属于电荷；
        * 目前只解析了一个原子/基团（如 ``Fe3+``、``(SO4)2-``）：数字属于电荷；
        * 多原子序列（如 ``NH4+``、``MnO4-``）：数字是下标，符号单独表示 ±1。

        括号内没有电荷歧义，数字始终是下标。
        """
        start = self.i
        self._skip_ws()
        space_before = self.i > start
        if space_before and self.i < len(self.s) and self.s[self.i].isdigit():
            # 下标不允许前置空白。若顶层的数字段后面跟正负号并结束，按电荷处理；
            # 其余情况直接报错。
            j = self.i
            while j < len(self.s) and self.s[j].isdigit():
                j += 1
            k = j
            while k < len(self.s) and self.s[k].isspace():
                k += 1
            if not self._stack and k < len(self.s) and self.s[k] in _SIGN and k + 1 == len(self.s):
                return 1
            raise FormulaError(
                f"意外的数字（位置 {self.i}）：下标前不能有空格，"
                "若想表示电荷请写成 \"SO4^2-\" 或 \"SO4 2-\""
            )
        digits = ""
        j = self.i
        while j < len(self.s) and self.s[j].isdigit():
            digits += self.s[j]
            j += 1
        if not digits:
            self.i = j
            return 1
        if not self._stack and j < len(self.s) and self.s[j] in _SIGN:
            # 跳过符号后的空白；后面没有字符才算结尾
            k = j + 1
            while k < len(self.s) and self.s[k].isspace():
                k += 1
            at_end = k == len(self.s)
            if at_end and (space_before or single_token):
                # 数字属于电荷，不消费
                return 1
        self.i = j
        value = int(digits)
        if value == 0:
            raise FormulaError(f"下标不能为 0（位置 {start}）")
        return value

    def _parse_count_for_electron(self) -> None:
        start = self.i
        n = self._parse_count(single_token=True)
        if n != 1:
            raise FormulaError(f"电子的个数必须为 1（位置 {start}）")

    # ------------------------------------------------------------------ 电荷
    def _parse_charge_with_leading_digit(self) -> int:
        """当前位置是数字，形如 ``2-`` / ``12+`` / ``SO4 2-`` 的纯文本电荷。"""
        start = self.i
        digits = ""
        while self.i < len(self.s) and self.s[self.i].isdigit():
            digits += self.s[self.i]
            self.i += 1
        self._skip_ws()
        if self.i >= len(self.s) or self.s[self.i] not in _SIGN:
            raise FormulaError(f"下标数字后缺少元素或括号（位置 {start}）")
        sign_char = self.s[self.i]
        self.i += 1
        self._skip_ws()
        if self.i < len(self.s):
            raise FormulaError(f"电荷符号 {sign_char!r} 后存在多余字符（位置 {self.i}）")
        magnitude = int(digits)
        if magnitude == 0:
            raise FormulaError("电荷不能为 0")
        return _SIGN[sign_char] * magnitude

    def _parse_charge_after_prefix(self) -> int:
        """解析电荷主体。

        进入本方法时，self.s[self.i]（跳过空白后）是以下之一：
        ``^``、``+``、``-``、数字。支持 ``^2-``、``^{-2}``、``^{2-}``、
        ``+``、``-3`` 等写法，结束后必须到达串尾。
        """
        self._skip_ws()
        if self.i < len(self.s) and self.s[self.i] == "^":
            self.i += 1
            self._skip_ws()
            if self.i >= len(self.s):
                raise FormulaError("缺少电荷值（^ 之后需要数字或 +/- 号）")
            if self.s[self.i] == "{":
                return self._parse_braced_charge()
        if self.i >= len(self.s):
            raise FormulaError("缺少电荷值")
        ch = self.s[self.i]
        if ch in _SIGN:
            sign = _SIGN[ch]
            self.i += 1
            self._skip_ws()
            magnitude = self._parse_tail_digits()
        elif ch.isdigit():
            start = self.i
            digits = ""
            while self.i < len(self.s) and self.s[self.i].isdigit():
                digits += self.s[self.i]
                self.i += 1
            self._skip_ws()
            if self.i >= len(self.s) or self.s[self.i] not in _SIGN:
                raise FormulaError(f"电荷数字之后必须是 + 或 - 号（位置 {start}）")
            sign = _SIGN[self.s[self.i]]
            self.i += 1
            self._skip_ws()
            magnitude = int(digits)
        else:
            raise FormulaError(f"电荷需要以 +/- 号或数字开头，得到 {ch!r}（位置 {self.i}）")
        if self.i < len(self.s):
            raise FormulaError(f"电荷之后存在多余字符（位置 {self.i}）")
        if magnitude == 0:
            raise FormulaError("电荷不能为 0")
        return sign * magnitude

    def _parse_braced_charge(self) -> int:
        """解析 ``^{...}`` 大括号包裹的电荷，内部允许 ``2-``、``-2``、``+`` 等。"""
        self.i += 1  # 消费 {
        depth_start = self.i
        body = ""
        while self.i < len(self.s) and self.s[self.i] != "}":
            if not self.s[self.i].isspace():
                body += self.s[self.i]
            self.i += 1
        if self.i >= len(self.s):
            raise FormulaError("电荷大括号缺少匹配的 }")
        self.i += 1  # 消费 }
        self._skip_ws()
        if self.i < len(self.s):
            raise FormulaError(f"电荷之后存在多余字符（位置 {self.i}）")
        if not body:
            raise FormulaError("大括号内缺少电荷值")
        sign = 0
        digits = ""
        for ch in body:
            if ch in _SIGN:
                if sign:
                    raise FormulaError(f"电荷中出现多个正负号（位置 {depth_start}）")
                sign = _SIGN[ch]
            elif ch.isdigit():
                digits += ch
            else:
                raise FormulaError(f"无法识别的电荷字符 {ch!r}（位置 {depth_start}）")
        if not sign:
            raise FormulaError("电荷必须包含 + 或 - 号")
        magnitude = int(digits) if digits else 1
        if magnitude == 0:
            raise FormulaError("电荷不能为 0")
        return sign * magnitude

    def _parse_tail_digits(self) -> int:
        digits = ""
        while self.i < len(self.s) and self.s[self.i].isdigit():
            digits += self.s[self.i]
            self.i += 1
        if not digits:
            return 1
        value = int(digits)
        if value == 0:
            raise FormulaError("电荷不能为 0")
        return value


# ---------------------------------------------------------------------- 辅助
def _accumulate(node: Node, multiplier: int, out: dict[str, int]) -> None:
    kind = node[0]
    if kind == "el":
        _, symbol, count = node
        out[symbol] = out.get(symbol, 0) + count * multiplier
    else:
        _, children, count = node
        for child in children:
            _accumulate(child, multiplier * count, out)


def _render_nodes(nodes: tuple[Node, ...] | list[Node]) -> str:
    parts: list[str] = []
    for node in nodes:
        if node[0] == "el":
            _, symbol, count = node
            parts.append(symbol if count == 1 else f"{symbol}{count}")
        else:
            _, children, count = node
            inner = _render_nodes(children)
            parts.append(f"({inner})" if count == 1 else f"({inner}){count}")
    return "".join(parts)


def _render_charge(prefix: str, charge: int) -> str:
    if charge == 0:
        return prefix
    magnitude = abs(charge)
    sign = "+" if charge > 0 else "-"
    suffix = sign if magnitude == 1 else f"{magnitude}{sign}"
    return f"{prefix}^{suffix}"
