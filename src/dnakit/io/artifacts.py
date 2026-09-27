"""Persistence for nested results, numeric arrays and file checksums."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit.core._json import to_json_compatible


@contextmanager
def _staged_path(target: str | os.PathLike[str], *, overwrite: bool) -> Iterator[Path]:
    """Stage a path beside its destination; replace only on successful exit."""
    path = Path(target)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(descriptor)
    temporary = Path(name)
    try:
        yield temporary
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def save_json(
    value: Any, target: str | os.PathLike[str], *, overwrite: bool = False, verbose: bool = True
) -> Path:
    """Save a nested JSON value or a result with ``to_dict()`` atomically.

    Unlike ``export_result``, this preserves the original object/array layout,
    without wrapping the result in a table schema. Non-finite numbers fail.
    """
    if hasattr(value, "to_dict"):
        value = value.to_dict()
    with (
        _staged_path(target, overwrite=overwrite) as temporary,
        temporary.open("w", encoding="utf-8", newline="") as stream,
    ):
        json.dump(to_json_compatible(value), stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return Path(target)


def load_json(source: str | os.PathLike[str], *, verbose: bool = True) -> Any:
    """Read a JSON object, array or scalar without changing its layout."""
    with Path(source).open(encoding="utf-8") as stream:
        return json.load(stream)


def save_arrays(
    arrays: Mapping[str, Any],
    target: str | os.PathLike[str],
    *,
    compressed: bool = True,
    overwrite: bool = False,
    verbose: bool = True,
) -> Path:
    """Save named arrays as NPZ; reject object arrays that would require pickle."""
    import numpy as np

    prepared = {name: np.asarray(value) for name, value in arrays.items()}
    if any(array.dtype.hasobject for array in prepared.values()):
        raise ValueError("Object arrays are unsupported; use numeric or string arrays.")
    with _staged_path(target, overwrite=overwrite) as temporary, temporary.open("wb") as stream:
        writer: Any = np.savez_compressed if compressed else np.savez
        writer(stream, **prepared)
    return Path(target)


def load_arrays(source: str | os.PathLike[str], *, verbose: bool = True) -> dict[str, Any]:
    """Load all NPZ arrays without pickle, closing the archive before returning."""
    import numpy as np

    with np.load(source, allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def file_digest(
    source: str | os.PathLike[str], algorithm: str = "sha256", *, verbose: bool = True
) -> str:
    """Compute a file digest with bounded memory."""
    digest = hashlib.new(algorithm)
    with Path(source).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
