"""Prokka annotation selecting its RNAmmer alternative rRNA predictor."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli


def annotate_prokka_rnammer(
    assembly: str | Path,
    output_dir: str | Path,
    *,
    prefix: str = "annotation",
    threads: int = 1,
    executable_path: str | None = None,
    timeout_seconds: float = 86400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run Prokka --rnammer, with an installed RNAmmer backend instead of Barrnap.

    Both Prokka and RNAmmer must be installed and configured by the caller. Prokka
    creates a new directory; existing output paths are refused. Returns GFF, GBK
    and FFN paths only when they really exist after a successful invocation.
    """
    source = cli.existing_file(assembly, "assembly")
    count = cli.integer(threads, "threads", minimum=1, maximum=256)
    name = cli.text(prefix, "prefix")
    if any(c in name for c in ("/", "\\")) or name in (".", ".."):
        from dnakit.exceptions import ConfigurationError

        raise ConfigurationError("prefix must be a simple file-name prefix.")
    target = cli.new_file(output_dir, "output_dir")
    return cli.run(
        "prokka",
        ["--rnammer", "--outdir", str(target), "--prefix", name, "--cpus", str(count), str(source)],
        cwd=target.parent,
        expected=[target / f"{name}.{ext}" for ext in ("gff", "gbk", "ffn")],
        explicit=executable_path,
        hint="Install Prokka and its optional RNAmmer backend; configure RNAmmer for Prokka.",
        timeout_seconds=timeout_seconds,
    )


__all__ = ["annotate_prokka_rnammer"]
_register_public_verbose(globals(), tuple(__all__))
