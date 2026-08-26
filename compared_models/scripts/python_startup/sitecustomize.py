from __future__ import annotations

import importlib.abc
import importlib.machinery
import os
import sys


_ENV_FLAG = "VISUALDNA_SKIP_NUMPY_BLAS_FPE_CHECK"


class _NumpyBlasFpePatchLoader(importlib.machinery.SourceFileLoader):
    def get_code(self, fullname):
        source_path = self.get_filename(fullname)
        source_bytes = self.get_data(source_path)
        return self.source_to_code(source_bytes, source_path)

    def source_to_code(self, data, path, *, _optimize=-1):
        if os.environ.get(_ENV_FLAG) == "1":
            text = data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else str(data)
            old = "    blas_fpe_check()\n    del blas_fpe_check\n"
            new = (
                f"    if os.environ.get({_ENV_FLAG!r}) != '1':\n"
                "        blas_fpe_check()\n"
                "    del blas_fpe_check\n"
            )
            if old in text:
                data = text.replace(old, new, 1)
        return super().source_to_code(data, path, _optimize=_optimize)


class _NumpyBlasFpePatchFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname != "numpy" or os.environ.get(_ENV_FLAG) != "1":
            return None

        try:
            sys.meta_path.remove(self)
            spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        finally:
            sys.meta_path.insert(0, self)

        if spec is None or spec.loader is None:
            return None
        if not isinstance(spec.loader, importlib.machinery.SourceFileLoader):
            return None

        spec.loader = _NumpyBlasFpePatchLoader(spec.loader.name, spec.loader.path)
        return spec


if os.environ.get(_ENV_FLAG) == "1":
    sys.meta_path.insert(0, _NumpyBlasFpePatchFinder())
