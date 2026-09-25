"""基于规范化签名的幂等缓存（SQLite 持久化，LRU 淘汰）。

缓存键与具体写法无关：只依赖每一侧“物种组成键 + 电荷”的多重集，
因此 ``SO4^2-`` 与 ``[SO4]2-``、不同括号风格、不同的物种书写顺序
都会命中同一条记录。条目数量有上限（默认 10000 条），按最近访问
（读或写）时间 LRU 淘汰。

与早期纯内存实现不同，记录持久化在 SQLite 文件中，容器重启后仍可
取回。路径由环境变量 ``REACTION_CACHE_DB`` 指定（默认内存库 ``:memory:``，
便于测试）；未设置或为 ``:memory:`` 时不跨进程保留。

并发与完整性：

* 所有写操作（写入记录、访问时间推进、超额淘汰）都放在同一个
  ``BEGIN IMMEDIATE`` 事务里提交，任何一步失败整体回滚，磁盘上不会
  留下半条记录；
* 读命中时“读取 + 更新最近访问顺序”同样在单个事务内完成；
* 进程内以互斥锁串行化事务，并设置 ``busy_timeout`` 以便与外部只读
  访问共存；
* 启动时建立表结构并校验存储版本（``cache_meta`` 中的
  ``schema_version``），遇到无法识别的版本抛出
  :class:`CacheSchemaVersionError`，拒绝在未知格式的库上继续运行。
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

from .balancer import Species

# 当前存储格式版本。表结构或序列化语义发生不兼容变化时递增；
# 启动时若数据库里记录的版本不在 SCHEMA_VERSION 可识别集合内则报错。
SCHEMA_VERSION = 1
# 本代码可以安全打开的存储版本集合（将来兼容旧版本时可扩展）。
SUPPORTED_SCHEMA_VERSIONS = (1,)

DEFAULT_MAXSIZE = 10_000
DB_PATH_ENV = "REACTION_CACHE_DB"

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


class CacheSchemaVersionError(RuntimeError):
    """SQLite 缓存库的存储版本无法被当前代码识别。"""


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


# 读行时选取的列顺序，须与 _to_record 解包顺序一致。
_COLUMNS = (
    "reaction_id",
    "signature",
    "canonical_reactants",
    "canonical_products",
    "coefficients",
    "error",
    "elements",
)


class ReactionCache:
    """SQLite 支撑的 LRU 幂等缓存，跨进程/容器重启可恢复。"""

    def __init__(
        self,
        db_path: Union[str, os.PathLike[str], None] = None,
        maxsize: int = DEFAULT_MAXSIZE,
    ) -> None:
        if maxsize < 1:
            raise ValueError("maxsize 必须为正整数")
        if db_path is None:
            db_path = os.environ.get(DB_PATH_ENV) or ":memory:"
        self._db_path = str(db_path)
        self._persistent = self._db_path != ":memory:"
        self._maxsize = maxsize
        self._lock = threading.RLock()

        if self._persistent:
            Path(self._db_path).parent.mkdir(parents=True, exist_ok=True)
        # check_same_thread=False + 自管理互斥锁：FastAPI 的同步端点可能
        # 在线程池中并发访问；同一连接跨线程使用由 self._lock 串行化。
        self._conn = sqlite3.connect(
            self._db_path, check_same_thread=False, isolation_level=None
        )
        # 外部进程短暂持锁时等待而不是立刻报 SQLITE_BUSY。
        self._conn.execute(f"PRAGMA busy_timeout = 5000")
        # 以事务（DELETE 日志）+ FULL 同步保证断电/杀容器后不会读到
        # 半条写入；WAL 在断电窗口里可能留下 -shm/-wal 残片，这里取
        # 最保守的持久性设置。
        if self._persistent:
            self._conn.execute("PRAGMA journal_mode = DELETE")
            self._conn.execute("PRAGMA synchronous = FULL")
        self._initialize()

    # ------------------------------------------------------------ 建表/版本
    def _initialize(self) -> None:
        """建立表结构并校验存储版本（单个原子事务）。

        先读版本再建表：若库由不兼容的旧/新版本写入，直接报错，不在
        未知格式的库上创建或修改任何对象。
        """
        with self._lock:
            # ATTACH 等场景外，用户表与 sqlite 内部表不重名；先探测
            # cache_meta 是否已存在以区分“新库”与“已有库”。
            exists = self._conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='cache_meta'"
            ).fetchone()
            if exists is not None:
                row = self._conn.execute(
                    "SELECT schema_version FROM cache_meta WHERE key='schema_version'"
                ).fetchone()
                version = row[0] if row is not None else None
                if version not in SUPPORTED_SCHEMA_VERSIONS:
                    raise CacheSchemaVersionError(
                        f"无法识别的缓存存储版本：database={version!r}，"
                        f"当前代码支持 {SUPPORTED_SCHEMA_VERSIONS}（库文件："
                        f"{self._db_path}）。请清空或迁移该缓存库后再启动。"
                    )
                return

            self._conn.execute("BEGIN IMMEDIATE")
            try:
                self._conn.execute(
                    "CREATE TABLE cache_meta ("
                    " key TEXT PRIMARY KEY,"
                    " schema_version INTEGER NOT NULL,"
                    " access_clock INTEGER NOT NULL"
                    ")"
                )
                self._conn.execute(
                    "INSERT INTO cache_meta (key, schema_version, access_clock)"
                    " VALUES ('schema_version', ?, 0)",
                    (SCHEMA_VERSION,),
                )
                self._conn.execute(
                    "CREATE TABLE reactions ("
                    " reaction_id TEXT PRIMARY KEY,"
                    " signature TEXT NOT NULL,"
                    " canonical_reactants TEXT NOT NULL,"
                    " canonical_products TEXT NOT NULL,"
                    " coefficients TEXT NOT NULL,"
                    " error TEXT,"
                    " elements TEXT NOT NULL,"
                    " created_at TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')),"
                    " last_accessed INTEGER NOT NULL"
                    ")"
                )
                # 淘汰查询 ORDER BY last_accessed 的支撑索引
                self._conn.execute(
                    "CREATE INDEX idx_reactions_last_accessed"
                    " ON reactions(last_accessed)"
                )
            except BaseException:
                self._conn.execute("ROLLBACK")
                raise
            else:
                self._conn.execute("COMMIT")

    def _tick(self) -> int:
        """在当前事务内推进并返回单调访问时钟（须已开启事务）。"""
        self._conn.execute(
            "UPDATE cache_meta SET access_clock = access_clock + 1"
            " WHERE key='schema_version'"
        )
        return self._conn.execute(
            "SELECT access_clock FROM cache_meta WHERE key='schema_version'"
        ).fetchone()[0]

    # ------------------------------------------------------------ 读取
    def get(self, reaction_id: str) -> ReactionRecord | None:
        """按 ID 取回记录；命中时在同一事务内把它标记为最近访问。"""
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                row = self._conn.execute(
                    f"SELECT {', '.join(_COLUMNS)} FROM reactions"
                    " WHERE reaction_id = ?",
                    (reaction_id,),
                ).fetchone()
                if row is None:
                    self._conn.execute("COMMIT")
                    return None
                now = self._tick()
                self._conn.execute(
                    "UPDATE reactions SET last_accessed = ? WHERE reaction_id = ?",
                    (now, reaction_id),
                )
            except BaseException:
                self._conn.execute("ROLLBACK")
                raise
            else:
                self._conn.execute("COMMIT")
            return _to_record(row)

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
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                now = self._tick()
                # UPSERT：同一 reaction_id 重复写入时整体替换并刷新访问时间，
                # 记录数不因重写而增加。
                self._conn.execute(
                    "INSERT INTO reactions ("
                    " reaction_id, signature, canonical_reactants,"
                    " canonical_products, coefficients, error, elements,"
                    " last_accessed)"
                    " VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
                    " ON CONFLICT(reaction_id) DO UPDATE SET"
                    " signature = excluded.signature,"
                    " canonical_reactants = excluded.canonical_reactants,"
                    " canonical_products = excluded.canonical_products,"
                    " coefficients = excluded.coefficients,"
                    " error = excluded.error,"
                    " elements = excluded.elements,"
                    " last_accessed = excluded.last_accessed",
                    (
                        record.reaction_id,
                        record.signature,
                        json.dumps(record.canonical_reactants, ensure_ascii=False),
                        json.dumps(record.canonical_products, ensure_ascii=False),
                        json.dumps(record.coefficients, ensure_ascii=False),
                        None
                        if record.error is None
                        else json.dumps(record.error, ensure_ascii=False),
                        json.dumps(record.elements, ensure_ascii=False),
                        now,
                    ),
                )
                # 同一事务内淘汰最久未访问记录，直到不超过上限。
                self._conn.execute(
                    "DELETE FROM reactions WHERE reaction_id IN ("
                    " SELECT reaction_id FROM reactions"
                    " ORDER BY last_accessed ASC, reaction_id ASC"
                    " LIMIT MAX(0, (SELECT COUNT(*) FROM reactions) - ?)"
                    ")",
                    (self._maxsize,),
                )
            except BaseException:
                self._conn.execute("ROLLBACK")
                raise
            else:
                self._conn.execute("COMMIT")
        return record

    def close(self) -> None:
        with self._lock:
            self._conn.close()

    def __len__(self) -> int:  # pragma: no cover - 便于运维查看
        with self._lock:
            return self._conn.execute("SELECT COUNT(*) FROM reactions").fetchone()[0]


def _to_record(row: tuple) -> ReactionRecord:
    (
        reaction_id,
        signature,
        reactants_json,
        products_json,
        coefficients_json,
        error_json,
        elements_json,
    ) = row
    return ReactionRecord(
        reaction_id=reaction_id,
        signature=signature,
        canonical_reactants=json.loads(reactants_json),
        canonical_products=json.loads(products_json),
        coefficients=json.loads(coefficients_json),
        error=None if error_json is None else json.loads(error_json),
        elements=json.loads(elements_json),
    )


# 模块导入即建立/校验缓存库（启动时建表并校验版本）；路径由
# REACTION_CACHE_DB 环境变量决定，缺省为进程内内存库。
cache = ReactionCache()
