"""Zarr engine: local directory stores through the optional ``zarr`` (>= 3) and ``numpy`` packages.

Every operation reads or writes a *directory* store.  Operations that produce a store (filter,
sort, merge, split, convert to zarr, write) require an explicit ``output`` / ``output_dir``.
"""

# ruff: noqa: RUF012, B905

from __future__ import annotations

import io
import math
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError
from dnakit.io_formats._impl_common import (
    DEFAULT_PAGE,
    Format,
    Opts,
    Payload,
    bad_input,
    bad_option,
    json_payload,
    looks_like_path,
    safe_label,
)


class OptionalBackendMissing(Exception):
    """Raised when zarr/numpy (or a recent enough zarr) is not importable."""

    def __init__(self, requires: str, message: str) -> None:
        super().__init__(message)
        self.requires = requires
        self.message = message


def _libs() -> tuple[Any, Any]:
    try:
        import numpy
        import zarr
    except ImportError as exc:
        raise OptionalBackendMissing(
            "pip install 'zarr>=3' numpy", f"zarr and numpy are required: {exc}"
        ) from exc
    major = int(re.match(r"\d+", getattr(zarr, "__version__", "0")).group(0))  # type: ignore[union-attr]
    if major < 3:
        raise OptionalBackendMissing(
            "pip install 'zarr>=3'",
            f"zarr {zarr.__version__} is installed but this adapter needs zarr >= 3.",
        )
    return zarr, numpy


_STORE_MARKERS = ("zarr.json", ".zgroup", ".zarray")


def _is_store(path: Path) -> bool:
    return path.is_dir() and any((path / marker).exists() for marker in _STORE_MARKERS)


def _open(source: object, mode: str = "r") -> tuple[Any, Path]:
    zarr, _ = _libs()
    if not (
        isinstance(source, os.PathLike) or (isinstance(source, str) and looks_like_path(source))
    ):
        raise bad_input(
            "A zarr input must be the path of a zarr directory store.",
            hint="Pass the directory path (a 'name.zarr' folder).",
        )
    path = Path(source).expanduser()
    if not _is_store(path):
        raise bad_input(
            f"{path} is not a zarr directory store (no zarr.json/.zgroup/.zarray).",
            code="IO_FORMATS_BAD_ZARR",
            hint="Pass the path of an existing zarr directory.",
        )
    try:
        return zarr.open(str(path), mode=mode), path
    except Exception as exc:
        raise bad_input(f"zarr cannot open {path}: {exc}", code="IO_FORMATS_BAD_ZARR") from exc


def _arrays(node: Any) -> dict[str, Any]:
    zarr, _ = _libs()
    if isinstance(node, zarr.Array):
        return {"/": node}
    return dict(
        sorted(
            (name, member)
            for name, member in node.members(max_depth=None)
            if isinstance(member, zarr.Array)
        )
    )


def _groups(node: Any) -> list[str]:
    zarr, _ = _libs()
    if isinstance(node, zarr.Array):
        return []
    return [name for name, member in node.members(max_depth=None) if isinstance(member, zarr.Group)]


def _info(array: Any) -> dict[str, object]:
    shape = [int(x) for x in array.shape]
    chunks = [int(x) for x in array.chunks]
    grid = [math.ceil(s / c) if c else 0 for s, c in zip(shape, chunks)]
    try:
        compressors = [type(c).__name__ for c in array.compressors]
    except Exception:
        compressors = []
    return {
        "shape": shape,
        "dtype": str(array.dtype),
        "chunks": chunks,
        "chunk_grid": grid,
        "nchunks": int(array.nchunks),
        "nchunks_initialized": int(array.nchunks_initialized),
        "compressors": compressors,
        "zarr_format": int(array.metadata.zarr_format),
        "nbytes": int(array.nbytes),
        "attrs": dict(array.attrs),
    }


def _to_list(values: Any) -> object:
    np = _libs()[1]
    arr = np.asarray(values)
    if arr.dtype.kind in "biuf":
        return arr.tolist()
    if arr.dtype.kind in "USmM":
        return arr.astype(str).tolist()
    raise bad_input(
        f"dtype {arr.dtype} cannot be converted to JSON values.", code="IO_FORMATS_BAD_ZARR"
    )


def _name(path: str) -> str:
    """Array name inside a new group: a root array ('/') is stored as 'data'."""

    return "data" if path == "/" else path.lstrip("/")


def _elements(array: Any) -> int:
    return int(math.prod(int(x) for x in array.shape)) if array.shape else 1


def _row_slice(spec: object, length: int) -> object:
    """``"start:stop[:step]"`` string, a mapping with start/stop/step, or a list of indices -> numpy-compatible index."""

    if isinstance(spec, str):
        match = re.fullmatch(r"\s*(-?\d*)\s*:\s*(-?\d*)\s*(?::\s*(-?\d*))?\s*", spec)
        if not match:
            raise bad_option(f"rows {spec!r} must look like 'start:stop[:step]'.")
        parts = [int(g) if g not in (None, "") else None for g in match.groups()]
        return slice(*parts)
    if isinstance(spec, Mapping):
        return slice(spec.get("start"), spec.get("stop"), spec.get("step"))
    if isinstance(spec, Sequence) and all(
        isinstance(i, int) and not isinstance(i, bool) for i in spec
    ):
        if any(not -length <= i < length for i in spec):
            raise bad_option(f"rows contains an index outside [0, {length}).")
        return list(spec)
    raise bad_option(
        "rows must be 'start:stop[:step]', {start, stop, step} or a list of integer indices."
    )


def _select(
    paths: dict[str, Any], keep: list[str] | None, drop: list[str] | None
) -> dict[str, Any]:
    def matches(path: str, patterns: list[str]) -> bool:
        return any(path == p.strip("/") or path.startswith(p.strip("/") + "/") for p in patterns)

    chosen = {
        p: a
        for p, a in paths.items()
        if (keep is None or matches(p.lstrip("/") or "/", keep))
        and (drop is None or not matches(p.lstrip("/") or "/", drop))
    }
    if keep is not None:
        missing = [k for k in keep if not any(matches(p.lstrip("/") or "/", [k]) for p in paths)]
        if missing:
            raise bad_option(
                f"Arrays not found in the store: {missing}", hint=f"Available: {list(paths)[:10]}"
            )
    return chosen


class ZarrFormat(Format):
    name = "zarr"
    title = "Zarr directory store"
    noun = "arrays"
    suffix = ".zarr"
    kind = "store"
    targets = ("zarr", "json", "tsv", "npz")
    native_ops = frozenset(
        {
            "read",
            "write",
            "validate",
            "filter",
            "sort",
            "merge",
            "split",
            "stream",
            "summarize",
            "index",
            "convert",
        }
    )
    _elements_opt = "int, default 5000000: refuse to materialise more array elements than this"
    _out = {
        "output": "directory path of the new zarr store (required)",
        "overwrite": "bool, default False: replace an existing zarr store at output",
    }
    options = {
        "read": {
            "array": "array path to read (default: all arrays)",
            "include_data": "bool, default True: return element values",
            "start": "int: first row along axis 0 (single array)",
            "stop": "int: end row (exclusive)",
            "max_elements": _elements_opt,
        },
        "write": {
            **_out,
            "array": "name of the array when the input is a single array-like (default 'data')",
            "chunks": "int or list of ints: chunk shape (default: chosen by zarr)",
            "attrs": "mapping stored as attributes of the root group",
            "zarr_format": "2 or 3 (default 3)",
        },
        "validate": {
            "max_chunk_bytes": "int, default 67108864: first chunks larger than this are not decoded"
        },
        "filter": {
            **_out,
            "arrays": "array paths (or group prefixes) to keep",
            "exclude": "array paths (or group prefixes) to drop",
            "rows": "rows to keep along axis 0: 'start:stop[:step]', {start, stop, step} or a list of indices "
            "(all selected arrays with at least one dimension must have the same length on axis 0)",
            "max_elements": _elements_opt,
        },
        "sort": {
            **_out,
            "by": "path of the 1-D key array (required); rows are ordered by its values (stable)",
            "reverse": "bool, default False: reverse the ascending order",
            "max_elements": _elements_opt,
        },
        "merge": {
            **_out,
            "axis": "int, default 0: axis to concatenate along",
            "max_elements": _elements_opt,
        },
        "split": {
            "output_dir": "directory that receives part_0001_<label>.zarr stores (required)",
            "overwrite": "bool, default False",
            "rows_per_part": "int (required): rows along axis 0 per part",
            "max_elements": _elements_opt,
        },
        "stream": {
            "array": "array path to page through (required when the store has several arrays)",
            "max_elements": "int, default 5000000: maximum elements per page",
        },
        "summarize": {},
        "index": {},
        "convert": {
            **_out,
            "array": "array path for the 'tsv' target (1-D or 2-D)",
            "max_elements": _elements_opt,
        },
    }
    # convert's `output` is optional for json/tsv/npz; requirement is checked per target
    options["convert"]["output"] = (
        "path to write (required for the 'zarr' target, optional otherwise)"
    )
    semantics = {
        "read": "Opens a zarr directory store read-only and returns array metadata (shape, dtype, chunks, codecs, attributes) "
        "and, unless include_data is false, the values as nested lists (bounded by max_elements).",
        "write": "Creates a zarr v3 (or v2) directory store from a {name: array-like} mapping, a list of {'name','data'} records or a "
        "single array-like; requires 'output'. Supported dtypes: bool, int, uint, float.",
        "validate": "Opens the store, checks every array's shape, chunk shape and chunk counts, and decodes the first chunk of each "
        "initialised array so that corrupt codecs or truncated chunk files are detected.",
        "filter": "Copies the selected arrays (by path / group prefix) and optionally only the chosen rows into a new store.",
        "sort": "Reorders axis 0 of every array whose first dimension equals the length of the key array (by) into a new store; "
        "other arrays are copied unchanged.",
        "merge": "Concatenates the arrays of several stores along an axis into a new store. Every store must contain the same "
        "array paths; shapes must agree except along the axis; dtypes are promoted with numpy.result_type.",
        "split": "Cuts axis 0 of all arrays into consecutive blocks of rows_per_part rows and writes one store per block.",
        "stream": "Pages through the rows along axis 0 of one array: records are the rows as nested lists; the cursor 'start' counts rows.",
        "summarize": "Counts groups/arrays and totals elements, bytes and initialised chunks; one entry per array.",
        "index": "Lists every array path with kind, shape, dtype, chunk shape, chunk grid and chunk counts (no chunk data is read).",
        "convert": "Targets: 'json' (metadata and values), 'tsv' (one 1-D/2-D array, tab separated), 'npz' (numpy archive of all "
        "arrays) and 'zarr' (a copy of the store at output, which is required).",
    }

    # The generic coerce/parse hooks are unused: every operation is native.
    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> Any:
        raise bad_input("zarr stores are directories; pass a path.")

    # -- helpers ----------------------------------------------------------------------------
    def _limit(self, opts: Opts) -> int:
        return opts.integer("max_elements", 5_000_000, 1) or 5_000_000

    def _target(self, opts: Opts, key: str = "output") -> Path:
        value = opts.text(key)
        if not value:
            raise ConfigurationError(
                f"This zarr operation writes a directory store; pass options={{'{key}': <path>}}.",
                code="IO_FORMATS_OUTPUT_REQUIRED",
            )
        return Path(value).expanduser()

    def _fresh(self, path: Path, overwrite: bool) -> None:
        if path.exists():
            if not overwrite:
                raise ConfigurationError(
                    f"Refusing to overwrite existing path: {path}",
                    code="IO_FORMATS_EXISTS",
                    hint="Set overwrite=True or choose another output.",
                )
            if not _is_store(path):
                raise ConfigurationError(
                    f"{path} exists and is not a zarr store; refusing to delete it.",
                    code="IO_FORMATS_EXISTS",
                )
            shutil.rmtree(path)
        path.parent.mkdir(parents=True, exist_ok=True)

    def _create(
        self,
        path: Path,
        arrays: Mapping[str, Any],
        attrs: Mapping[str, Any] | None,
        zarr_format: int = 3,
        chunks: dict[str, Any] | None = None,
        group_attrs: Mapping[str, Mapping[str, Any]] | None = None,
    ) -> None:
        zarr, np = _libs()
        group = zarr.open_group(str(path), mode="w-", zarr_format=zarr_format)
        if attrs:
            group.attrs.update(dict(attrs))
        for name, data in arrays.items():
            data = np.asarray(data)
            cs = (chunks or {}).get(name)
            kwargs: dict[str, Any] = {}
            if cs is not None and data.ndim:
                kwargs["chunks"] = tuple(
                    max(1, min(int(c), int(s))) if s else 1 for c, s in zip(cs, data.shape)
                )
            array = group.create_array(name, data=data, **kwargs)
            if (group_attrs or {}).get(name):
                array.attrs.update(dict(group_attrs[name]))  # type: ignore[index]

    def _snapshot(
        self, node: Any, opts: Opts
    ) -> tuple[dict[str, Any], dict[str, Any], dict[str, dict[str, Any]], dict[str, Any]]:
        """Return (numpy arrays, chunk shapes, array attrs, group attrs) read from a node, bounded by max_elements."""

        _, np = _libs()
        arrays = _arrays(node)
        limit = self._limit(opts)
        if sum(_elements(a) for a in arrays.values()) > limit:
            raise bad_option(
                f"The store has more than max_elements={limit} elements.",
                hint="Raise max_elements.",
            )
        data = {p: np.asarray(a[...]) for p, a in arrays.items()}
        return (
            data,
            {p: a.chunks for p, a in arrays.items()},
            {p: dict(a.attrs) for p, a in arrays.items()},
            dict(getattr(node, "attrs", {})),
        )

    # -- operations -------------------------------------------------------------------------
    def read_op(self, source: object, opts: Opts) -> dict[str, object]:
        node, path = _open(source)
        wanted = opts.text("array")
        include = opts.boolean("include_data", True)
        start, stop = opts.integer("start", None, 0), opts.integer("stop", None, 0)
        limit = self._limit(opts)
        arrays = _arrays(node)
        if wanted:
            arrays = _select(arrays, [wanted], None)
        if (start is not None or stop is not None) and len(arrays) != 1:
            raise bad_option("start/stop need exactly one array (use the 'array' option).")
        total = sum(_elements(a) for a in arrays.values())
        if include and total > limit and start is None and stop is None:
            raise bad_option(
                f"Reading {total} elements exceeds max_elements={limit}.",
                hint="Raise max_elements, select one array, use start/stop, or pass include_data=False.",
            )
        out: dict[str, object] = {}
        for name, array in arrays.items():
            entry = _info(array)
            if include:
                view = (
                    array[start:stop]
                    if (start is not None or stop is not None) and array.shape
                    else array[...]
                )
                if int(view.size) > limit:
                    raise bad_option(
                        f"Selection has {int(view.size)} elements, above max_elements={limit}."
                    )
                entry["data"] = _to_list(view)
            out[name] = entry
        return {
            "path": str(path),
            "arrays": out,
            "groups": _groups(node),
            "attrs": dict(getattr(node, "attrs", {})),
            "record_count": len(out),
            "records": [
                {"path": p, **{k: v for k, v in e.items() if k != "data"}} for p, e in out.items()
            ],
        }  # type: ignore[union-attr]

    def validate_op(self, source: object, opts: Opts) -> dict[str, object]:
        errors: list[str] = []
        warnings: list[str] = []
        max_chunk = opts.integer("max_chunk_bytes", 64 * 1024 * 1024, 1) or 64 * 1024 * 1024
        try:
            node, path = _open(source)
        except ConfigurationError:
            raise
        except Exception as exc:
            if getattr(exc, "code", "") == "IO_FORMATS_BAD_ZARR":
                return {
                    "valid": False,
                    "errors": [getattr(exc, "message", str(exc))],
                    "warnings": [],
                    "array_count": 0,
                    "group_count": 0,
                }
            raise
        arrays = _arrays(node)
        for name, array in arrays.items():
            info = _info(array)
            if any(s < 0 for s in info["shape"]) or any(c <= 0 for c in info["chunks"]):  # type: ignore[union-attr]
                errors.append(
                    f"{name}: invalid shape {info['shape']} or chunk shape {info['chunks']}"
                )
                continue
            if info["nchunks_initialized"] > info["nchunks"]:  # type: ignore[operator]
                errors.append(f"{name}: more initialised chunks than the chunk grid allows")
            if info["nchunks_initialized"] == 0 and _elements(array):
                warnings.append(
                    f"{name}: no chunk has been written (the array reads as fill values)"
                )
                continue
            chunk_bytes = math.prod(info["chunks"]) * array.dtype.itemsize  # type: ignore[arg-type]
            if _elements(array) and chunk_bytes <= max_chunk:
                try:
                    array[tuple(slice(0, min(c, s)) for c, s in zip(info["chunks"], info["shape"]))]  # type: ignore[arg-type]
                except Exception as exc:
                    errors.append(f"{name}: first chunk cannot be decoded: {exc}")
        return {
            "valid": not errors,
            "errors": errors,
            "warnings": warnings,
            "array_count": len(arrays),
            "group_count": len(_groups(node)),
            "path": str(path),
        }

    def summarize_op(self, source: object, opts: Opts) -> dict[str, object]:
        node, path = _open(source)
        arrays = _arrays(node)
        infos = {p: _info(a) for p, a in arrays.items()}
        return {
            "path": str(path),
            "kind": "array" if "/" in arrays else "group",
            "array_count": len(arrays),
            "group_count": len(_groups(node)),
            "total_elements": sum(_elements(a) for a in arrays.values()),
            "total_nbytes": sum(i["nbytes"] for i in infos.values()),  # type: ignore[misc]
            "initialized_chunks": sum(i["nchunks_initialized"] for i in infos.values()),  # type: ignore[misc]
            "total_chunks": sum(i["nchunks"] for i in infos.values()),  # type: ignore[misc]
            "attrs": dict(getattr(node, "attrs", {})),
            "arrays": infos,
        }

    def index_op(self, source: object, opts: Opts) -> dict[str, object]:
        node, path = _open(source)
        entries = {
            p: {
                k: v
                for k, v in _info(a).items()
                if k in ("shape", "dtype", "chunks", "chunk_grid", "nchunks", "nchunks_initialized")
            }
            for p, a in _arrays(node).items()
        }
        return {
            "key": "array path",
            "record_count": len(entries),
            "entries": entries,
            "groups": _groups(node),
            "path": str(path),
        }

    def stream_op(self, source: object, opts: Opts) -> dict[str, object]:
        node, _ = _open(source)
        arrays = _arrays(node)
        name = opts.text("array")
        if name:
            arrays = _select(arrays, [name], None)
        if len(arrays) != 1:
            raise bad_option(
                f"stream needs exactly one array; the store has {list(arrays)[:10]}.",
                hint="Pass the 'array' option.",
            )
        ((path, array),) = arrays.items()
        if not array.shape:
            raise bad_option("Cannot page through a 0-dimensional array.")
        start = opts.integer("start", 0, 0) or 0
        limit = opts.integer("limit", DEFAULT_PAGE, 1, 100_000) or DEFAULT_PAGE
        total = int(array.shape[0])
        row_size = max(1, math.prod(int(x) for x in array.shape[1:]))
        cap = self._limit(opts)
        rows = min(limit, max(0, total - start))
        if rows * row_size > cap:
            raise bad_option(
                f"The page would hold {rows * row_size} elements, above max_elements={cap}.",
                hint="Lower limit.",
            )
        values = _to_list(array[start : start + rows]) if rows else []
        end = start + rows
        return {
            "array": path,
            "start": start,
            "limit": limit,
            "returned": rows,
            "total": total,
            "next_start": end if end < total else None,
            "done": end >= total,
            "records": values,
        }

    # -- store-producing operations -----------------------------------------------------------
    def _result(self, path: Path, count: int, **info: object) -> Payload:
        return Payload("zarr", path=str(path), count=count, info=dict(info))

    def filter_op(self, source: object, opts: Opts) -> Payload:
        node, _ = _open(source)
        keep, drop = opts.names("arrays"), opts.names("exclude")
        rows = opts.get("rows")
        target = self._target(opts)
        overwrite = opts.boolean("overwrite", False)
        arrays = _select(_arrays(node), keep, drop)
        if not arrays:
            raise bad_option("No array is left after the arrays/exclude selection.")
        data, chunks, array_attrs, group_attrs = self._snapshot(node, opts)
        data = {p: data[p] for p in arrays}
        kept_rows = None
        if rows is not None:
            lengths = {d.shape[0] for d in data.values() if d.ndim}
            if len(lengths) != 1:
                raise bad_option(
                    "rows needs all selected arrays to share one length on axis 0.",
                    hint="Narrow the selection with 'arrays'.",
                )
            index = _row_slice(rows, lengths.pop())
            data = {p: (d[index] if d.ndim else d) for p, d in data.items()}
            kept_rows = int(next(d for d in data.values() if d.ndim).shape[0])
        self._fresh(target, overwrite)
        self._create(
            target,
            {_name(p): d for p, d in data.items()},
            group_attrs,
            chunks={_name(p): chunks[p] for p in data},
            group_attrs={_name(p): array_attrs[p] for p in data},
        )
        return self._result(target, len(data), kept_rows=kept_rows)

    def sort_op(self, source: object, opts: Opts) -> Payload:
        node, _ = _open(source)
        by = opts.text("by")
        reverse = opts.boolean("reverse", False)
        target = self._target(opts)
        overwrite = opts.boolean("overwrite", False)
        if not by:
            raise bad_option("sort_zarr needs the 'by' option: the path of a 1-D key array.")
        _, np = _libs()
        data, chunks, array_attrs, group_attrs = self._snapshot(node, opts)
        key_path = next((p for p in data if p.lstrip("/") == by.strip("/")), None)
        if key_path is None:
            raise bad_option(
                f"Key array {by!r} is not in the store.", hint=f"Available: {list(data)[:10]}"
            )
        key = data[key_path]
        if key.ndim != 1:
            raise bad_option("The key array must be 1-dimensional.")
        order = np.argsort(key, kind="stable")
        if reverse:
            order = order[::-1]
        permuted, unchanged = [], []
        for p, d in list(data.items()):
            if d.ndim and d.shape[0] == key.shape[0]:
                data[p] = d[order]
                permuted.append(p)
            else:
                unchanged.append(p)
        self._fresh(target, overwrite)
        self._create(
            target,
            {_name(p): d for p, d in data.items()},
            group_attrs,
            chunks={_name(p): c for p, c in chunks.items()},
            group_attrs={_name(p): a for p, a in array_attrs.items()},
        )
        return self._result(
            target, len(data), key=key_path, permuted_arrays=permuted, copied_unchanged=unchanged
        )

    def merge_op(self, sources: list[object], opts: Opts) -> Payload:
        _, np = _libs()
        axis = opts.integer("axis", 0, 0) or 0
        target = self._target(opts)
        overwrite = opts.boolean("overwrite", False)
        snapshots = []
        total = 0
        for source in sources:
            node, _ = _open(source)
            snapshots.append(self._snapshot(node, opts))
            total += sum(d.size for d in snapshots[-1][0].values())
            if total > self._limit(opts):
                raise bad_option(
                    f"The inputs hold more than max_elements={self._limit(opts)} elements.",
                    hint="Raise max_elements.",
                )
        paths = list(snapshots[0][0])
        for number, snap in enumerate(snapshots[1:], 2):
            if list(snap[0]) != paths and set(snap[0]) != set(paths):
                raise bad_input(
                    f"Store {number} has different array paths than store 1: {sorted(set(snap[0]) ^ set(paths))[:5]}",
                    code="IO_FORMATS_BAD_ZARR",
                )
        merged = {}
        for p in paths:
            parts = [snap[0][p] for snap in snapshots]
            if any(part.ndim <= axis for part in parts):
                raise bad_input(
                    f"Array {p!r} has too few dimensions to concatenate along axis {axis}."
                )
            try:
                merged[p] = np.concatenate(parts, axis=axis)
            except ValueError as exc:
                raise bad_input(
                    f"Array {p!r} cannot be concatenated along axis {axis}: {exc}",
                    code="IO_FORMATS_BAD_ZARR",
                ) from exc
        self._fresh(target, overwrite)
        self._create(
            target,
            {_name(p): d for p, d in merged.items()},
            snapshots[0][3],
            chunks={_name(p): c for p, c in snapshots[0][1].items()},
            group_attrs={_name(p): a for p, a in snapshots[0][2].items()},
        )
        return self._result(target, len(merged), input_count=len(sources), axis=axis)

    def split_op(self, source: object, opts: Opts) -> list[tuple[str, Payload]]:
        node, _ = _open(source)
        size = opts.integer("rows_per_part", None, 1)
        directory = self._target(opts, "output_dir")
        overwrite = opts.boolean("overwrite", False)
        if size is None:
            raise bad_option("split_zarr needs rows_per_part.")
        data, chunks, array_attrs, group_attrs = self._snapshot(node, opts)
        lengths = {d.shape[0] for d in data.values() if d.ndim}
        if len(lengths) != 1:
            raise bad_option(
                "split_zarr needs all arrays to share one length on axis 0.",
                hint="Use filter_zarr to select compatible arrays first.",
            )
        total = lengths.pop()
        parts = []
        for number, a in enumerate(range(0, total, size), 1):
            label = f"rows_{a}-{min(a + size, total)}"
            target = directory / f"part_{number:04d}_{safe_label(label)}.zarr"
            self._fresh(target, overwrite)
            self._create(
                target,
                {_name(p): (d[a : a + size] if d.ndim else d) for p, d in data.items()},
                group_attrs,
                chunks={_name(p): c for p, c in chunks.items()},
                group_attrs={_name(p): v for p, v in array_attrs.items()},
            )
            parts.append((label, self._result(target, len(data), row_count=min(size, total - a))))
        return parts

    def convert_op(self, source: object, opts: Opts) -> Payload:
        target_fmt = opts.text("to")
        if target_fmt not in self.targets:
            raise bad_option(f"convert_zarr: 'to' must be one of {list(self.targets)}.")
        node, path = _open(source)
        _, np = _libs()
        if target_fmt == "json":
            limit = self._limit(opts)
            arrays = _arrays(node)
            if sum(_elements(a) for a in arrays.values()) > limit:
                raise bad_option(
                    f"The store has more than max_elements={limit} elements.",
                    hint="Raise max_elements.",
                )
            body = {
                "attrs": dict(getattr(node, "attrs", {})),
                "arrays": {p: {**_info(a), "data": _to_list(a[...])} for p, a in arrays.items()},
            }
            return json_payload("zarr", body, len(arrays))
        if target_fmt == "tsv":
            arrays = _arrays(node)
            name = opts.text("array")
            if name:
                arrays = _select(arrays, [name], None)
            if len(arrays) != 1:
                raise bad_option(
                    "The 'tsv' target needs exactly one array (use the 'array' option)."
                )
            ((_, array),) = arrays.items()
            if array.ndim > 2:
                raise bad_option("The 'tsv' target supports 1-D and 2-D arrays only.")
            if _elements(array) > self._limit(opts):
                raise bad_option("The array is larger than max_elements.")
            values = np.asarray(array[...])
            rows = values.reshape(-1, 1) if values.ndim <= 1 else values
            lines = (
                ["\t".join(str(x) for x in row) for row in rows.astype(str).tolist()]
                if values.ndim
                else [str(values)]
            )
            return Payload("tsv", text="".join(line + "\n" for line in lines), count=len(lines))
        if target_fmt == "npz":
            data, _, _, _ = self._snapshot(node, opts)
            buffer = io.BytesIO()
            np.savez(
                buffer, **{p.strip("/").replace("/", "__") or "data": d for p, d in data.items()}
            )
            return Payload("npz", data=buffer.getvalue(), count=len(data))
        # 'zarr': copy into output
        data, chunks, array_attrs, group_attrs = self._snapshot(node, opts)
        target = self._target(opts)
        self._fresh(target, opts.boolean("overwrite", False))
        self._create(
            target,
            {_name(p): d for p, d in data.items()},
            group_attrs,
            chunks={_name(p): c for p, c in chunks.items()},
            group_attrs={_name(p): a for p, a in array_attrs.items()},
        )
        del path
        return self._result(target, len(data))

    def write_op(self, source: object, opts: Opts) -> Payload:
        _, np = _libs()
        target = self._target(opts)
        overwrite = opts.boolean("overwrite", False)
        single = opts.text("array", "data")
        chunks = opts.get("chunks")
        attrs = opts.mapping("attrs")
        fmt = opts.integer("zarr_format", 3, 2, 3) or 3
        pairs: list[tuple[str, Any]]
        if (
            isinstance(source, Mapping)
            and source
            and all(isinstance(k, str) for k in source)
            and "records" not in source
        ):
            pairs = [(str(k), v) for k, v in source.items()]
        elif isinstance(source, Mapping) and "records" in source:
            return self.write_op(source["records"], opts)
        elif (
            isinstance(source, Sequence)
            and not isinstance(source, (str, bytes))
            and source
            and all(isinstance(i, Mapping) and "name" in i and "data" in i for i in source)
        ):
            pairs = [(str(i["name"]), i["data"]) for i in source]
        elif isinstance(source, (Sequence, np.ndarray)) and not isinstance(source, (str, bytes)):
            pairs = [(single or "data", source)]
        else:
            raise bad_input(
                "write_zarr needs a {name: array} mapping, a list of {'name','data'} records or an array-like."
            )
        arrays = {}
        for name, value in pairs:
            if not name or name.startswith("/") or ".." in name.split("/"):
                raise bad_input(f"Invalid array name {name!r}.")
            try:
                arr = np.asarray(value)
            except ValueError as exc:
                raise bad_input(f"Array {name!r} is ragged or not array-like: {exc}") from exc
            if arr.dtype.kind not in "biuf":
                raise bad_input(
                    f"Array {name!r} has dtype {arr.dtype}; only bool, int, uint and float arrays are supported."
                )
            arrays[name] = arr
        per_chunks: dict[str, Any] = {}
        if chunks is not None:
            spec = [chunks] if isinstance(chunks, int) and not isinstance(chunks, bool) else chunks
            if not isinstance(spec, Sequence) or not all(
                isinstance(c, int) and not isinstance(c, bool) and c > 0 for c in spec
            ):
                raise bad_option(
                    "chunks must be a positive integer or a list of positive integers."
                )
            for name, arr in arrays.items():
                per_chunks[name] = (
                    tuple(spec) if len(spec) == arr.ndim else tuple([spec[0]] * arr.ndim)
                )
        self._fresh(target, overwrite)
        self._create(target, arrays, attrs, zarr_format=fmt, chunks=per_chunks)
        return self._result(target, len(arrays), zarr_format=fmt)
