"""基于规范化签名的幂等缓存。

缓存键与具体写法无关：只依赖每一侧“物种组成键 + 电荷”的多重集，
因此 ``SO4^2-`` 与 ``[SO4]2-``、不同括号风格、不同的物种书写顺序
都会命中同一条记录。条目数量有上限，按 LRU 淘汰。
"""

from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from dataclasses import dataclass, field

from .balancer import Species

# 内存中的一条“事实”：与书写顺序无关的配平系数映射与错误（若有）。
# 键: (侧别, 物种规范键)；值: 系数。错误记录携带 code/message/details。


@dataclass
class ReactionRecord:
    reaction_id: str
    signature: str
    canonical_reactants: list[dict]          # [{normalized, key}]
    canonical_products: list[dict]
    coefficients: dict[str, int] = field(default_factory=dict)
    error: dict | None = None
    elements: list[str] = field(default_factory=list)


def species_key(sp: Species) -> str:
    """物种规范键（组成 + 电荷）的可序列化字符串。"""
    comp = ",".join(f"{el}:{n}" for el, n in sp.formula.key[0])
    return f"[{comp}]{sp.formula.charge:+d}"


def build_signature(reactants: list[Species], products: list[Species]) -> str:
    payload = {
        "r": sorted((species_key(sp), sp.normalized) for sp in reactants),
        "p": sorted((species_key(sp), sp.normalized) for sp in products),
    }
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def reaction_id_for(signature: str) -> str:
    digest = hashlib.sha256(signature.encode("utf-8")).hexdigest()[:16]
    return f"rxn_{digest}"


class ReactionCache:
    def __init__(self, maxsize: int = 10_000) -> None:
        self._data: OrderedDict[str, ReactionRecord] = OrderedDict()
        self._maxsize = maxsize
        self._lock = threading.Lock()

    def get(self, reaction_id: str) -> ReactionRecord | None:
        with self._lock:
            record = self._data.get(reaction_id)
            if record is not None:
                self._data.move_to_end(reaction_id)
            return record

    def get_by_signature(self, signature: str) -> ReactionRecord | None:
        return self.get(reaction_id_for(signature))

    def put(
        self,
        signature: str,
        reactants: list[Species],
        products: list[Species],
        coefficients: dict[str, int] | None,
        elements: list[str],
        error: dict | None,
    ) -> ReactionRecord:
        reaction_id = reaction_id_for(signature)
        record = ReactionRecord(
            reaction_id=reaction_id,
            signature=signature,
            canonical_reactants=[
                {"normalized": sp.normalized, "key": species_key(sp)}
                for sp in sorted(reactants, key=lambda s: s.formula.key)
            ],
            canonical_products=[
                {"normalized": sp.normalized, "key": species_key(sp)}
                for sp in sorted(products, key=lambda s: s.formula.key)
            ],
            coefficients=coefficients or {},
            error=error,
            elements=list(elements),
        )
        with self._lock:
            self._data[reaction_id] = record
            self._data.move_to_end(reaction_id)
            while len(self._data) > self._maxsize:
                self._data.popitem(last=False)
        return record

    def __len__(self) -> int:  # pragma: no cover - 便于运维查看
        with self._lock:
            return len(self._data)


cache = ReactionCache()
