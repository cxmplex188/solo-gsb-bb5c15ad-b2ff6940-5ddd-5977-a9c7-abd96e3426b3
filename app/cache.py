"""基于规范化签名的幂等缓存（SQLite 持久化）。

缓存键与具体写法无关：只依赖每一侧“物种组成键 + 电荷”的多重集，
因此 ``SO4^2-`` 与 ``[SO4]2-``、不同括号风格、不同的物种书写顺序
都会命中同一条记录。条目数量有上限，按 LRU 淘汰。

与早期内存版不同，记录落盘在 SQLite 数据库中（默认
``data/reaction_cache.db``，可用环境变量 ``REACTION_CACHE_DB`` 覆盖，
传 ``:memory:`` 则退化为纯内存库），容器 / 进程重启后缓存仍可命中：

* 读取记录时在事务内把该记录的访问序号刷新为最新（LRU 的“最近访问”
  语义），访问顺序更新与记录读取不会在并发请求间交错出半新状态；
* 写入、访问序号更新以及超过 ``maxsize`` 条时淘汰最久未访问记录都在
  单个 SQLite 事务中完成——并发请求或进程中断不会留下半条记录；
* 建库时创建表并校验 ``meta.schema_version``，版本无法识别时抛出明确
  的 :class:`RuntimeError`，避免静默读写不兼容的存量数据库。
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import threading
from dataclasses import dataclass, field

from .balancer import Species

# 存储格式版本：结构变更时递增，旧库在启动时会被明确拒绝而不是误读。
SCHEMA_VERSION = "1"

DEFAULT_DB_PATH = os.path.join("data", "reaction_cache.db")

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


def _dump(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


class ReactionCache:
    """SQLite 持久化的反应记录缓存（LRU，默认上限 1 万条）。

    所有写路径（写入记录、刷新访问序号、淘汰最久未访问记录）都在
    ``with self._conn`` 开启的事务里提交；进程内用一把锁串行化，
    跨进程则依赖 SQLite 事务与 WAL 保证原子性。
    """

    def __init__(self, db_path: str | None = None, maxsize: int = 10_000) -> None:
        if db_path is None:
            db_path = os.environ.get("REACTION_CACHE_DB", DEFAULT_DB_PATH)
        self._db_path = db_path
        self._maxsize = maxsize
        self._lock = threading.Lock()
        if db_path != ":memory:":
            os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
        # check_same_thread=False：FastAPI 线程池 / 测试客户端都可能跨线程
        # 访问；并发安全由 self._lock 与 SQLite 事务共同保证。
        self._conn = sqlite3.connect(db_path, check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._initialize()

    # ------------------------------------------------------------ 启动初始化
    def _initialize(self) -> None:
        """建表并校验存储版本；版本无法识别时明确报错。"""
        with self._lock, self._conn:
            self._conn.execute(
                "CREATE TABLE IF NOT EXISTS meta ("
                "  key TEXT PRIMARY KEY,"
                "  value TEXT NOT NULL"
                ")"
            )
            row = self._conn.execute(
                "SELECT value FROM meta WHERE key = 'schema_version'"
            ).fetchone()
            if row is not None and row[0] != SCHEMA_VERSION:
                raise RuntimeError(
                    f"无法识别的配平缓存存储版本 {row[0]!r}（本服务仅支持 "
                    f"{SCHEMA_VERSION!r}）：请迁移或删除缓存数据库 {self._db_path}"
                )
            self._conn.execute(
                "CREATE TABLE IF NOT EXISTS reactions ("
                "  reaction_id TEXT PRIMARY KEY,"
                "  signature TEXT NOT NULL,"
                "  canonical_reactants TEXT NOT NULL,"
                "  canonical_products TEXT NOT NULL,"
                "  coefficients TEXT NOT NULL,"
                "  error TEXT,"
                "  elements TEXT NOT NULL,"
                "  last_accessed INTEGER NOT NULL"
                ")"
            )
            if row is None:
                self._conn.execute(
                    "INSERT INTO meta (key, value) VALUES ('schema_version', ?)",
                    (SCHEMA_VERSION,),
                )
                self._conn.execute(
                    "INSERT INTO meta (key, value) VALUES ('access_seq', '0')"
                )

    def _tick(self) -> int:
        """取下一个访问序号（调用方须已持有锁且处于事务中）。"""
        self._conn.execute(
            "UPDATE meta SET value = CAST(value AS INTEGER) + 1"
            " WHERE key = 'access_seq'"
        )
        row = self._conn.execute(
            "SELECT value FROM meta WHERE key = 'access_seq'"
        ).fetchone()
        return int(row[0])

    # ------------------------------------------------------------ 查询
    _COLUMNS = (
        "reaction_id, signature, canonical_reactants, canonical_products,"
        " coefficients, error, elements"
    )

    @staticmethod
    def _row_to_record(row: tuple) -> ReactionRecord:
        return ReactionRecord(
            reaction_id=row[0],
            signature=row[1],
            canonical_reactants=json.loads(row[2]),
            canonical_products=json.loads(row[3]),
            coefficients={key: int(value) for key, value in json.loads(row[4]).items()},
            error=json.loads(row[5]) if row[5] is not None else None,
            elements=json.loads(row[6]),
        )

    def get(self, reaction_id: str) -> ReactionRecord | None:
        with self._lock:
            row = self._conn.execute(
                f"SELECT {self._COLUMNS} FROM reactions WHERE reaction_id = ?",
                (reaction_id,),
            ).fetchone()
            if row is None:
                return None
            # 命中即刷新最近访问顺序；序号分配与更新在同一事务中提交。
            with self._conn:
                seq = self._tick()
                self._conn.execute(
                    "UPDATE reactions SET last_accessed = ? WHERE reaction_id = ?",
                    (seq, reaction_id),
                )
            return self._row_to_record(row)

    def get_by_signature(self, signature: str) -> ReactionRecord | None:
        return self.get(reaction_id_for(signature))

    # ------------------------------------------------------------ 写入
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
        # 写入、访问序号刷新与超上限淘汰在同一事务中提交。
        with self._lock, self._conn:
            seq = self._tick()
            self._conn.execute(
                "INSERT OR REPLACE INTO reactions ("
                "  reaction_id, signature, canonical_reactants,"
                "  canonical_products, coefficients, error, elements,"
                "  last_accessed"
                ") VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    record.reaction_id,
                    record.signature,
                    _dump(record.canonical_reactants),
                    _dump(record.canonical_products),
                    _dump(record.coefficients),
                    _dump(record.error) if record.error is not None else None,
                    _dump(record.elements),
                    seq,
                ),
            )
            excess = (
                self._conn.execute("SELECT COUNT(*) FROM reactions").fetchone()[0]
                - self._maxsize
            )
            if excess > 0:
                self._conn.execute(
                    "DELETE FROM reactions WHERE reaction_id IN ("
                    "  SELECT reaction_id FROM reactions"
                    "  ORDER BY last_accessed ASC, reaction_id ASC LIMIT ?"
                    ")",
                    (excess,),
                )
        return record

    def __len__(self) -> int:  # 便于运维查看
        with self._lock:
            return self._conn.execute("SELECT COUNT(*) FROM reactions").fetchone()[0]

    def close(self) -> None:  # 测试与优雅停机用
        with self._lock:
            self._conn.close()


cache = ReactionCache()
