"""Import and signature checks for the coverage-gap modules added for Fig. 2 (see paper/write_paper/.../coverage_gap)."""

from __future__ import annotations

import importlib
import inspect

import pytest

MODULES = [
    "dnakit.bioinformatics.pysam_gap_functions",
    "dnakit.bioinformatics.pysam_gap2_functions",
    "dnakit.databases.bio_network_gap_functions",
    "dnakit.epigenome.cooler_gap_functions",
    "dnakit.epigenome.macs3_cli_functions",
    "dnakit.epigenome.deeptools_cli_facades",
    "dnakit.epigenome.deeptools_cli_facades2",
    "dnakit.epigenome.deeptools_cli_facades3",
    "dnakit.epigenome.deeptools_matrix_operations",
    "dnakit.epigenome.deeptools_gc_bias",
    "dnakit.metagenomics_sourmash_facades",
    "dnakit.molecular_lab.dnachisel_gap_functions",
    "dnakit.phylogeny.dendropy_cli_functions",
    "dnakit.phylogeny.toytree_gap_functions",
    "dnakit.population.dadi_gap_functions",
    "dnakit.population.tskit_cli_functions",
    "dnakit.population.tskit_gap_functions",
    "dnakit.sequence.bio_gap_functions",
]


@pytest.mark.parametrize("module_name", MODULES)
def test_module_imports_and_exports_callables(module_name: str) -> None:
    module = importlib.import_module(module_name)
    names = getattr(module, "__all__", None)
    assert names, f"{module_name} exports no __all__"
    for name in names:
        assert callable(getattr(module, name)), f"{module_name}.{name} is not callable"


@pytest.mark.parametrize("module_name", MODULES)
def test_public_functions_take_keyword_only_verbose_and_are_typed(module_name: str) -> None:
    module = importlib.import_module(module_name)
    for name in module.__all__:
        function = getattr(module, name)
        signature = inspect.signature(function)
        verbose = signature.parameters.get("verbose")
        assert verbose is not None, f"{module_name}.{name} lacks a verbose parameter"
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY, f"{module_name}.{name}: verbose must be keyword-only"
        assert function.__doc__, f"{module_name}.{name} has no docstring"
        for parameter in signature.parameters.values():
            if parameter.name in ("self", "cls"):
                continue
            assert parameter.annotation is not inspect.Parameter.empty, (
                f"{module_name}.{name}: parameter {parameter.name} is not annotated"
            )
