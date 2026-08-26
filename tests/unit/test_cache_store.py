"""Tests for canonical content-addressed cache persistence."""

import json
from pathlib import Path

import pytest

from dnakit.cache import CacheKey, JSONCache
from dnakit.exceptions import CacheError, ConfigurationError


def test_cache_key_is_order_independent_and_values_roundtrip(tmp_path: Path) -> None:
    first = CacheKey.from_components("descriptor", {"sequence": "ACGT", "params": {"k": 2}})
    second = CacheKey.from_components("descriptor", {"params": {"k": 2}, "sequence": "ACGT"})
    cache = JSONCache(tmp_path)

    assert first == second
    artifact = cache.put(first, {"value": [1, 2], "method": "count"})
    assert cache.get(first) == {"method": "count", "value": [1, 2]}
    assert artifact.byte_size > 0
    assert cache.invalidate(first)
    assert cache.get(first) is None
    assert not cache.invalidate(first)


def test_cache_key_includes_installed_dnakit_version(monkeypatch: pytest.MonkeyPatch) -> None:
    import dnakit.cache.store as store

    monkeypatch.setattr(store, "_dnakit_version", lambda: "1.0.0")
    first = CacheKey.from_components("descriptor", {"sequence": "ACGT"})
    monkeypatch.setattr(store, "_dnakit_version", lambda: "2.0.0")
    second = CacheKey.from_components("descriptor", {"sequence": "ACGT"})

    assert first.digest != second.digest


def test_cache_detects_payload_tampering(tmp_path: Path) -> None:
    key = CacheKey.from_components("test", {"input": "A"})
    cache = JSONCache(tmp_path / "cache")
    artifact = cache.put(key, {"value": 1})
    path = Path(artifact.relative_path)
    envelope = json.loads(path.read_text(encoding="utf-8"))
    envelope["payload"]["value"] = 2
    path.write_text(json.dumps(envelope), encoding="utf-8")

    with pytest.raises(CacheError) as exc_info:
        cache.get(key)
    assert exc_info.value.code == "CACHE_INTEGRITY_ERROR"


def test_cache_clear_is_namespace_bounded(tmp_path: Path) -> None:
    cache = JSONCache(tmp_path / "cache")
    first = CacheKey.from_components("one", {"x": 1})
    second = CacheKey.from_components("two", {"x": 1})
    cache.put(first, 1)
    cache.put(second, 2)

    report = cache.clear("one")

    assert report.removed_count == 1
    assert cache.get(first) is None
    assert cache.get(second) == 2


@pytest.mark.parametrize("namespace", ["", "../escape", "two words", "/absolute"])
def test_cache_rejects_unsafe_namespace(namespace: str) -> None:
    with pytest.raises(ConfigurationError):
        CacheKey.from_components(namespace, {"x": 1})


def test_cache_refuses_existing_unmarked_directory_and_namespace_symlink(tmp_path: Path) -> None:
    occupied = tmp_path / "occupied"
    occupied.mkdir()
    (occupied / "keep.json").write_text("keep", encoding="utf-8")

    with pytest.raises(CacheError) as root_error:
        JSONCache(occupied)
    assert root_error.value.code == "UNSAFE_CACHE_ROOT"
    assert (occupied / "keep.json").read_text(encoding="utf-8") == "keep"

    cache_root = tmp_path / "cache"
    cache = JSONCache(cache_root)
    external = tmp_path / "external"
    external.mkdir()
    (cache_root / "linked").symlink_to(external, target_is_directory=True)
    key = CacheKey.from_components("linked", {"x": 1})

    with pytest.raises(CacheError) as namespace_error:
        cache.put(key, 1)
    assert namespace_error.value.code == "UNSAFE_CACHE_NAMESPACE"


def test_cache_rejects_invalid_public_types_and_entry_symlinks(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        CacheKey.from_components("test", [])  # type: ignore[arg-type]
    cache = JSONCache(tmp_path / "cache")
    with pytest.raises(ConfigurationError):
        cache.clear([])  # type: ignore[arg-type]

    key = CacheKey.from_components("test", {"x": 1})
    cache.put(key, {"x": 1})
    path = cache.root / key.namespace / f"{key.digest}.json"
    external = tmp_path / "external.json"
    external.write_text("{}", encoding="utf-8")
    path.unlink()
    path.symlink_to(external)

    with pytest.raises(CacheError) as error:
        cache.get(key)
    assert error.value.code == "UNSAFE_CACHE_ENTRY"


def test_cache_entry_byte_limit_applies_before_write_and_read(tmp_path: Path) -> None:
    cache = JSONCache(tmp_path / "cache", max_entry_bytes=300)
    key = CacheKey.from_components("bounded", {"x": 1})

    with pytest.raises(CacheError) as write_error:
        cache.put(key, {"large": "x" * 1_000})
    assert write_error.value.code == "CACHE_ENTRY_SIZE_LIMIT"
    assert cache.get(key) is None

    permissive = JSONCache(tmp_path / "cache", max_entry_bytes=2_000)
    artifact = permissive.put(key, {"small": True})
    path = Path(artifact.relative_path)
    path.write_bytes(path.read_bytes() + b" " * 2_000)

    with pytest.raises(CacheError) as read_error:
        cache.get(key)
    assert read_error.value.code == "CACHE_ENTRY_SIZE_LIMIT"


def test_cache_deep_json_failure_is_structured(tmp_path: Path) -> None:
    cache = JSONCache(tmp_path / "cache")
    key = CacheKey.from_components("bounded", {"x": 1})
    namespace = cache.root / key.namespace
    namespace.mkdir()
    (namespace / ".dnakit-cache-namespace-v1").touch()
    path = namespace / f"{key.digest}.json"
    path.write_text("[" * 1_200 + "]" * 1_200, encoding="utf-8")

    with pytest.raises(CacheError) as error:
        cache.get(key)
    assert error.value.code == "CACHE_DECODE_ERROR"


def test_cache_rejects_deep_values_and_duplicate_serialized_keys(tmp_path: Path) -> None:
    deep: object = "leaf"
    for _ in range(100):
        deep = [deep]
    with pytest.raises(ConfigurationError) as key_error:
        CacheKey.from_components("bounded", {"deep": deep})
    assert key_error.value.code == "INVALID_CACHE_VALUE"

    cache = JSONCache(tmp_path / "cache")
    key = CacheKey.from_components("bounded", {"x": 1})
    artifact = cache.put(key, {"value": 1})
    path = Path(artifact.relative_path)
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace('"payload":', '"payload":null,"payload":'), encoding="utf-8")
    with pytest.raises(CacheError) as duplicate_error:
        cache.get(key)
    assert duplicate_error.value.code == "CACHE_DECODE_ERROR"


def test_cache_accepts_shared_dag_but_rejects_a_true_cycle(tmp_path: Path) -> None:
    cache = JSONCache(tmp_path / "cache")
    key = CacheKey.from_components("bounded", {"x": 1})
    shared = {"x": 1}

    cache.put(key, {"left": shared, "right": shared})
    assert cache.get(key) == {"left": {"x": 1}, "right": {"x": 1}}

    cycle: list[object] = []
    cycle.append(cycle)
    with pytest.raises(ConfigurationError) as error:
        cache.put(key, cycle)
    assert error.value.code == "CACHE_STRUCTURE_LIMIT"


@pytest.mark.parametrize("limit", [True, 0, 1_000_000_001])
def test_cache_entry_byte_limit_is_strict(tmp_path: Path, limit: object) -> None:
    with pytest.raises(ConfigurationError) as error:
        JSONCache(tmp_path / "cache", max_entry_bytes=limit)  # type: ignore[arg-type]
    assert error.value.code == "INVALID_CACHE_ENTRY_LIMIT"
