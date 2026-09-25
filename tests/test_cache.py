"""SQLite 持久幂等缓存测试。"""

import sqlite3

import pytest

from app.balancer import parse_side
from app.cache import (
    SCHEMA_VERSION,
    ReactionCache,
    build_signature,
    reaction_id_for,
    species_key,
)


def _put(cache: ReactionCache, reactants, products, error=None):
    """以真实解析结果写入一条记录，返回 (signature, record)。"""
    parsed_r = list(parse_side(reactants, "反应物"))
    parsed_p = list(parse_side(products, "生成物"))
    signature = build_signature(parsed_r, parsed_p)
    coefficients = (
        None
        if error is not None
        else {species_key(sp): 1 for sp in (*parsed_r, *parsed_p)}
    )
    record = cache.put(
        signature,
        parsed_r,
        parsed_p,
        coefficients,
        ["H", "O"],
        error,
    )
    return signature, record


def test_put_get_roundtrip(tmp_path):
    cache = ReactionCache(str(tmp_path / "cache.db"))
    signature, record = _put(cache, ["H2", "O2"], ["H2O"])

    by_id = cache.get(record.reaction_id)
    by_signature = cache.get_by_signature(signature)
    assert by_id is not None and by_signature is not None
    assert by_id.reaction_id == record.reaction_id == reaction_id_for(signature)
    assert by_signature.reaction_id == record.reaction_id
    assert by_id.signature == signature
    assert by_id.coefficients == record.coefficients
    assert by_id.error is None
    assert by_id.elements == ["H", "O"]
    assert [item["normalized"] for item in by_id.canonical_reactants] == [
        "H2",
        "O2",
    ]
    assert [item["normalized"] for item in by_id.canonical_products] == ["H2O"]
    cache.close()


def test_error_record_roundtrip(tmp_path):
    cache = ReactionCache(str(tmp_path / "cache.db"))
    error = {"code": "IMPOSSIBLE_BALANCE", "message": "守恒无解", "details": {}}
    _, record = _put(cache, ["H2O"], ["Au"], error=error)

    fetched = cache.get(record.reaction_id)
    assert fetched is not None
    assert fetched.error == error
    assert fetched.coefficients == {}
    cache.close()


def test_persistence_across_instances(tmp_path):
    """同一数据库文件重新建缓存实例后记录仍可命中（模拟容器重启）。"""
    path = str(tmp_path / "cache.db")
    first = ReactionCache(path)
    _, record = _put(first, ["Ca(OH)2", "HCl"], ["CaCl2", "H2O"])
    first.close()

    second = ReactionCache(path)
    fetched = second.get(record.reaction_id)
    assert fetched is not None
    assert fetched.reaction_id == record.reaction_id
    assert fetched.signature == record.signature
    assert fetched.coefficients == record.coefficients
    assert len(second) == 1
    second.close()


def test_get_refreshes_recency_and_lru_eviction(tmp_path):
    """读取会把记录刷成最近访问；超过上限时淘汰最久未访问的记录。"""
    cache = ReactionCache(str(tmp_path / "cache.db"), maxsize=3)
    _, rec_a = _put(cache, ["H2", "O2"], ["H2O"])
    _, rec_b = _put(cache, ["Fe", "O2"], ["Fe2O3"])
    _, rec_c = _put(cache, ["C", "O2"], ["CO2"])

    # 读取 A 使其成为最近访问；随后写入 D，应淘汰最久未访问的 B
    assert cache.get(rec_a.reaction_id) is not None
    _, rec_d = _put(cache, ["N2", "H2"], ["NH3"])

    assert len(cache) == 3
    assert cache.get(rec_a.reaction_id) is not None
    assert cache.get(rec_b.reaction_id) is None
    assert cache.get(rec_c.reaction_id) is not None
    assert cache.get(rec_d.reaction_id) is not None
    cache.close()


def test_put_existing_moves_to_most_recent(tmp_path):
    """重复写入同一签名会覆盖记录并刷新访问顺序，不产生重复条目。"""
    cache = ReactionCache(str(tmp_path / "cache.db"), maxsize=2)
    signature, rec_a = _put(cache, ["H2", "O2"], ["H2O"])
    _, rec_b = _put(cache, ["Fe", "O2"], ["Fe2O3"])
    # 再次写入 A（等价写法、同签名），随后写入 C：应淘汰 B 而不是 A
    _put(cache, ["O2", "H2"], ["H2O"])
    _, rec_c = _put(cache, ["C", "O2"], ["CO2"])

    assert len(cache) == 2
    assert cache.get(rec_a.reaction_id) is not None
    assert cache.get(rec_b.reaction_id) is None
    assert cache.get(rec_c.reaction_id) is not None
    cache.close()


def test_unknown_schema_version_rejected(tmp_path):
    """存量数据库版本无法识别时，启动建库必须明确报错。"""
    path = tmp_path / "cache.db"
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    conn.execute("INSERT INTO meta (key, value) VALUES ('schema_version', '0')")
    conn.commit()
    conn.close()

    with pytest.raises(RuntimeError, match="无法识别的配平缓存存储版本"):
        ReactionCache(str(path))


def test_schema_version_recorded(tmp_path):
    path = tmp_path / "cache.db"
    cache = ReactionCache(str(path))
    cache.close()
    conn = sqlite3.connect(path)
    row = conn.execute(
        "SELECT value FROM meta WHERE key = 'schema_version'"
    ).fetchone()
    conn.close()
    assert row is not None and row[0] == SCHEMA_VERSION


def test_missing_reaction_returns_none(tmp_path):
    cache = ReactionCache(str(tmp_path / "cache.db"))
    assert cache.get("rxn_doesnotexist") is None
    cache.close()
