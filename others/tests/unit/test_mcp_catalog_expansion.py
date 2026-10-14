"""Catalog entries added on 2026-10-14: pinned sources, launch recipes and recorded checks."""

from __future__ import annotations

import re

import pytest

from dnakit.integrations import list_mcp_servers
from dnakit.integrations import mcp as bridge

NEW_IDS = {
    "seqeron-sequence": "MIT",
    "smartapi-biothings": "Apache-2.0",
    "just-prs": "MIT",
    "dnallm": "MIT",
    "just-dna-agents": "MIT",
    "genomi": "Apache-2.0",
    "galaxy-mcp": "MIT",
    "gnomad": "Apache-2.0",
    "addgene": "MIT",
    "alphagenome-jolab": "MIT",
    "locus": "MIT",
    "sequence-simulation": "BSD-2-Clause",
    "ncbi-sra": "MIT",
    "bioinformatics-mcp": "MIT",
    "rewire-genomics": "MIT",
    "bioinformatics-ncbi-uniprot": "MIT",
    "uniprot-provenance": "Apache-2.0",
}
# Excluded after reading the upstream code or licence files; they must not creep back in.
EXCLUDED_REPOSITORIES = {
    "bach-ai-tools/biotools-mcp-server",
    "dna-seq/ensembl-mcp",
    "longevity-genie/alphagenome-mcp",
    "longevity-genie/benchling-mcp",
    "sian-davies/mcp-ncbi-blast",
    "bio-mcp/bio-mcp-blast",
    "jvogan/motif",
    "cafferychen777/chatspatial",
    "fastmcp-me/uniprot-mcp",
}


@pytest.fixture(scope="module")
def catalog() -> dict[str, dict]:
    return {server["id"]: server for server in list_mcp_servers(verbose=False)}


def test_new_entries_are_present_with_the_licence_read_from_the_upstream_file(catalog) -> None:
    assert NEW_IDS.keys() <= catalog.keys()
    for identifier, licence in NEW_IDS.items():
        assert catalog[identifier]["license"] == licence, identifier


def test_new_entries_are_commit_pinned_and_use_installed_commands(catalog) -> None:
    for identifier in NEW_IDS:
        server = catalog[identifier]
        sha = server["source_commit"]
        assert re.fullmatch(r"[0-9a-f]{40}", sha), identifier
        assert all(sha in url for url in server["sources"]), identifier
        assert server["verified_date"] == "2026-10-14"
        assert server["command"] not in {"uv", "uvx", "npx", "npm", "bunx", "pip"}
        assert server["transports"] == ["stdio"]
        assert server["effects"] and server["installation"] and server["capabilities"]


def test_every_new_entry_records_tool_names_and_a_startup_status(catalog) -> None:
    for identifier in NEW_IDS:
        server = catalog[identifier]
        assert len(server["tool_names"]) >= 3 and len(set(server["tool_names"])) == len(
            server["tool_names"]
        )
        assert server["startup_check"].startswith("2026-10-14:")


def test_started_entries_name_the_tools_that_were_listed(catalog) -> None:
    started = {k: v for k, v in catalog.items() if k in NEW_IDS and "listed" in v["startup_check"]}
    assert started.keys() == NEW_IDS.keys() - {
        "dnallm",
        "just-dna-agents",
        "seqeron-sequence",
        "sequence-simulation",
        "alphagenome-jolab",
    }
    counts = {
        k: int(re.search(r"listed (\d+) tools", v["startup_check"]).group(1))
        for k, v in started.items()
    }
    for identifier, count in counts.items():
        assert count >= len(catalog[identifier]["tool_names"]), identifier


def test_secrets_are_declared_and_never_fixed(catalog) -> None:
    assert catalog["alphagenome-jolab"]["secret_environment"] == ["ALPHAGENOME_API_KEY"]
    assert "non-commercial" in catalog["alphagenome-jolab"]["effects"].lower()
    assert "GALAXY_API_KEY" in catalog["galaxy-mcp"]["secret_environment"]
    for identifier in NEW_IDS:
        server = catalog[identifier]
        assert not set(server.get("fixed_environment", {})) & set(server["secret_environment"])
        assert set(server["secret_environment"]) <= set(
            server["required_environment"]
            + server["optional_environment"]
            + [name for pair in server["environment_pairs"] for name in pair]
        )


def test_source_based_entries_require_the_files_they_launch(catalog) -> None:
    assert catalog["gnomad"]["required_files"] == ["server.py"]
    assert catalog["sequence-simulation"]["arguments"] == ["dist/server.js"]
    assert catalog["ncbi-sra"]["arguments"][0] == "-c"
    assert "sra_mcp_server.py" in catalog["ncbi-sra"]["required_files"]
    for identifier in ("seqeron-sequence", "dnallm", "gnomad", "sequence-simulation", "ncbi-sra"):
        assert catalog[identifier]["required_files"], identifier
    for identifier in ("gnomad", "sequence-simulation", "ncbi-sra"):
        with pytest.raises(bridge.ConfigurationError, match="workdir"):
            bridge._connection(bridge._definition(identifier), "stdio", None, None, None)


def test_known_upstream_problems_are_kept_in_the_catalog(catalog) -> None:
    assert "mcp<2" in catalog["bioinformatics-ncbi-uniprot"]["installation"]
    assert any("socksio" in note for note in catalog["smartapi-biothings"]["limitations"])
    assert any("CC-BY-NC" in note for note in catalog["locus"]["limitations"])


def test_excluded_repositories_are_not_registered(catalog) -> None:
    registered = {
        server["repository"].lower().split("github.com/")[-1] for server in catalog.values()
    }
    assert not (registered & EXCLUDED_REPOSITORIES)
