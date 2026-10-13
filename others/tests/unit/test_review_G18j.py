"""Offline wrappers for pyensembl and genomepy; empty gget and biothings_client review (G18j)."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.databases import biothings_client_review_functions as bc_api
from dnakit.databases import genomepy_review_functions as gp_api
from dnakit.databases import gget_review_functions as gg_api
from dnakit.databases import pyensembl_review_functions as pe_api
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

ALL_WRAPPERS = [
    pe_api.ensembl_dated_release_urls,
    pe_api.ensembl_dated_genome_fasta_mask,
    pe_api.ensembl_genome_fasta_source_options,
    pe_api.ensembl_species_record,
    gp_api.map_genome_contig_names,
]
ASSEMBLY_REPORT = (
    "Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
    "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    "1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t"
    "248956422\tchr1\n"
    "2\tassembled-molecule\t2\tChromosome\tCM000664.2\t=\tNC_000002.12\tPrimary Assembly\t"
    "242193529\tchr2\n"
    "MT\tassembled-molecule\tMT\tMitochondrion\tJ01415.2\t=\tNC_012920.1\tnon-nuclear\t"
    "16569\tchrM\n"
    "KI270706.1\tunplaced-scaffold\tna\tUnplaced-Scaffold\tKI270706.1\t=\tna\tUnlocalized\t"
    "175055\tna\n"
)


def _install_genome(root: Path, name: str, provider: str, report: str | None = ASSEMBLY_REPORT):
    """Create a genome folder the way genomepy leaves it: a README and, optionally, a report."""
    files = pytest.importorskip("genomepy.files")
    folder = root / name
    folder.mkdir(parents=True)
    metadata = {
        "name": name,
        "provider": provider,
        "assembly_accession": "GCF_000001405.40",
        "tax_id": "9606",
        "species": "Homo sapiens",
        "genome_file": f"{name}.fa",
    }
    files.write_readme(str(folder / "README.txt"), metadata)
    if report is not None:
        (folder / "assembly_report.txt").write_text(report, encoding="utf-8")
    return root


def test_every_wrapper_has_native_keyword_only_verbose() -> None:
    for function in ALL_WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_dormant_modules_export_nothing() -> None:
    assert list(gg_api.__all__) == []
    assert list(bc_api.__all__) == []


def test_dated_release_urls_match_hand_built_strings_and_upstream() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import ensembl_url_templates as upstream

    result = pe_api.ensembl_dated_release_urls(
        "GCA_000001405.29",
        "ensembl",
        "2026_04",
        include_alt=True,
        genome_fasta_mask="soft",
        verbose=False,
    )
    base = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04"
    assert result["directory"] == base
    assert result["gtf"] == f"{base}/geneset/genes-including_alt.gtf.gz"
    assert result["cdna"] == f"{base}/geneset/cdna.fa.bgz"
    assert result["pep"] == f"{base}/geneset/pep.fa.bgz"
    assert result["genome_fasta"] == f"{base}/genome/softmasked.fa.bgz"
    expected = upstream.make_dated_release_urls(
        "GCA_000001405.29", "ensembl", "2026_04", include_alt=True, genome_fasta_mask="soft"
    )
    assert result["gtf"] == expected.gtf and result["genome_fasta"] == expected.genome_fasta
    assert result["provenance"]["network"] == "none; URL string only"


def test_dated_release_urls_reject_bad_input() -> None:
    pytest.importorskip("pyensembl")
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_dated_release_urls("GCA_1.1", "ensembl", "2026_04", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_dated_release_urls("GCA_000001405.29", "ensembl", "2026-04", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_dated_release_urls(
            "GCA_000001405.29", "ensembl", "2026_04", genome_fasta_mask="dna", verbose=False
        )
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_dated_release_urls(
            "GCA_000001405.29", "ensembl", "2026_04", include_alt=1, verbose=False
        )


def test_dated_genome_fasta_mask_reads_each_mask_and_rejects_other_urls() -> None:
    pytest.importorskip("pyensembl")
    base = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome"
    expected = {"unmasked": "none", "softmasked": "soft", "hardmasked": "hard"}
    for file_stem, mask in expected.items():
        result = pe_api.ensembl_dated_genome_fasta_mask(f"{base}/{file_stem}.fa.bgz", verbose=False)
        assert result["genome_fasta_mask"] == mask
        assert result["matched"] is True
    other = pe_api.ensembl_dated_genome_fasta_mask(f"{base}/cdna.fa.bgz", verbose=False)
    assert other["genome_fasta_mask"] is None and other["matched"] is False


def test_source_options_read_type_and_mask_from_official_urls() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import genome_fasta_cache as upstream

    soft = (
        "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
        "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
    )
    plain = (
        "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    )
    result = pe_api.ensembl_genome_fasta_source_options(soft, verbose=False)
    assert (result["canonical"], result["genome_fasta_type"], result["genome_fasta_mask"]) == (
        True,
        "toplevel",
        "soft",
    )
    assert upstream.canonical_source_options(soft) == ("toplevel", "soft")
    plain_result = pe_api.ensembl_genome_fasta_source_options(plain, verbose=False)
    assert plain_result["genome_fasta_type"] == "primary_assembly"
    assert plain_result["genome_fasta_mask"] == "none"
    local = pe_api.ensembl_genome_fasta_source_options("/data/genome.fa.gz", verbose=False)
    assert local["canonical"] is False and local["genome_fasta_type"] is None


def test_species_record_matches_registry_by_name_and_reference() -> None:
    pytest.importorskip("pyensembl")
    by_name = pe_api.ensembl_species_record("human", verbose=False)
    assert by_name["latin_name"] == "homo_sapiens"
    assert by_name["synonyms"] == ["human"]
    assert {"reference": "GRCh38", "first_release": 76, "last_release": 116} in (
        by_name["reference_assemblies"]
    )
    assert by_name["dated_releases"] == {
        "assembly_accession": "GCA_000001405.29",
        "provider": "ensembl",
    }
    by_reference = pe_api.ensembl_species_record("grch37", by="reference", verbose=False)
    assert by_reference["latin_name"] == "homo_sapiens"
    assert by_reference["by"] == "reference"


def test_species_record_rejects_unknown_species_and_bad_mode() -> None:
    pytest.importorskip("pyensembl")
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_species_record("Ailuropoda melanoleuca", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_species_record("NoSuchReference", by="reference", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_api.ensembl_species_record("human", by="taxon", verbose=False)


def test_map_contig_names_matches_the_saved_report(tmp_path: Path) -> None:
    pytest.importorskip("genomepy")
    root = _install_genome(tmp_path, "demo", "NCBI")
    ucsc = gp_api.map_genome_contig_names("demo", "UCSC", genomes_dir=root, verbose=False)
    assert ucsc["available"] is True and ucsc["source_provider"] == "ncbi"
    assert [(r["source_name"], r["target_name"]) for r in ucsc["records"]] == [
        ("1", "chr1"),
        ("2", "chr2"),
        ("KI270706.1", None),
        ("MT", "chrM"),
    ]
    ensembl = gp_api.map_genome_contig_names(
        "demo", "ensembl", genomes_dir=str(root), verbose=False
    )
    assert ensembl["record_count"] == 4
    assert ensembl["records"][0] == {"source_name": "1", "target_name": "1"}


def test_map_contig_names_reports_no_mapping_for_the_same_provider(tmp_path: Path) -> None:
    pytest.importorskip("genomepy")
    root = _install_genome(tmp_path, "demo", "NCBI")
    result = gp_api.map_genome_contig_names("demo", "ncbi", genomes_dir=root, verbose=False)
    assert result["available"] is False and result["records"] == []


def test_map_contig_names_never_downloads_a_missing_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    providers = pytest.importorskip("genomepy.providers")

    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("network download attempted")

    monkeypatch.setattr(providers, "download_assembly_report", forbidden)
    root = _install_genome(tmp_path, "demo", "NCBI", report=None)
    with pytest.raises(BackendUnavailableError):
        gp_api.map_genome_contig_names("demo", "ucsc", genomes_dir=root, verbose=False)


def test_map_contig_names_rejects_bad_arguments(tmp_path: Path) -> None:
    pytest.importorskip("genomepy")
    root = _install_genome(tmp_path, "demo", "NCBI")
    with pytest.raises(ConfigurationError):
        gp_api.map_genome_contig_names("../demo", "ucsc", genomes_dir=root, verbose=False)
    with pytest.raises(ConfigurationError):
        gp_api.map_genome_contig_names("demo", "refseq", genomes_dir=root, verbose=False)
    with pytest.raises(ConfigurationError):
        gp_api.map_genome_contig_names("demo", "ucsc", genomes_dir=tmp_path / "missing")
    with pytest.raises(BackendUnavailableError):
        gp_api.map_genome_contig_names("absent", "ucsc", genomes_dir=root, verbose=False)
