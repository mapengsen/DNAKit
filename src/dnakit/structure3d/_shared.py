"""Shared provenance for native coordinate analysis and 3DNA/DSSR adapters."""

from __future__ import annotations

from dnakit.core import (
    Citation,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
)


def native_structure_provenance() -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.NOVEL,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.DNAKIT,
        )
    )


def threedna_provenance() -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.HYBRID,
            origin_class=OriginClass.INTEGRATION,
            citations=(
                Citation(
                    "3dna2003",
                    title="3DNA: analysis, rebuilding and visualization of 3D nucleic acids",
                    doi="10.1093/nar/gkg680",
                ),
                Citation(
                    "dssr2015",
                    title="DSSR: an integrated tool for nucleic acid spatial structure",
                    doi="10.1093/nar/gkv716",
                ),
            ),
        )
    )


__all__ = ["native_structure_provenance", "threedna_provenance"]
