"""Actionable dependency diagnostics, printed before a caller raises its error."""

from __future__ import annotations

import sys


def dependency_hint(extras: str | None = None, *, guidance: str | None = None) -> str:
    """Print installation guidance to stderr and retain it in the exception hint.

    These are error diagnostics, so they remain visible with ``verbose=False``.
    No packages are installed here; compatible installed dependencies are left
    to pip's normal resolver rather than forcing upgrades or reinstallations.
    """
    parts = ["Required dependencies are missing or could not be imported."]
    if extras is not None:
        parts.extend(
            (
                f'From the DNAKit project root: python -m pip install -e ".[{extras}]"',
                f'For a packaged installation: python -m pip install "dnakit[{extras}]"',
            )
        )
    if guidance:
        parts.append(guidance)
    parts.append(
        "Compatible installed packages are reused; incompatible versions may need adjustment."
    )
    parts.append("Install in the Python environment running this script, then retry.")
    hint = "\n".join(parts)
    print(f"[DNAKit] {hint}", file=sys.stderr, flush=True)
    return hint


def model_dependency_hint(model: str, requirement: str | None = None) -> str:
    """Select the complete optional dependency group for a neural backend."""
    name = model.lower().replace(" ", "").replace("-", "")
    if name.startswith("enformer"):
        return dependency_hint("neural,neural-enformer")
    if name.startswith("evo2"):
        return dependency_hint(
            "neural,neural-evo2",
            guidance="The Evo 2 extra requires Python >=3.11,<3.13 "
            "and its official CUDA environment.",
        )
    if name.startswith("caduceus"):
        return dependency_hint("neural,neural-caduceus")
    if name.startswith("alphagenome"):
        return dependency_hint(
            guidance=(
                "Install the official AlphaGenome research dependencies with: "
                'python -m pip install "git+https://github.com/google-deepmind/'
                'alphagenome_research.git". Follow its JAX/device setup instructions, '
                "or pass model_source_path to a checkout with its dependencies installed."
            )
        )
    if name.startswith("lucaone"):
        return dependency_hint(
            guidance=(
                "Clone https://github.com/LucaOne/LucaOneTasks, install its official "
                "environment, and pass model_source_path to that checkout."
            )
        )
    if name.startswith("janusdna"):
        return dependency_hint(
            "neural",
            guidance="Clone https://github.com/Qihao-Duan/JanusDNA and install its official "
            "environment; pass model_source_path to that checkout.",
        )
    return dependency_hint(
        "neural", guidance=f"Required backend: {requirement}" if requirement else None
    )
