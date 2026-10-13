"""Tests for the G8 database wrappers: pyensembl, genomepy, pybiomart and ncbi.datasets.

Each wrapped function runs the real upstream package on small fixtures created in tmp_path.
Network paths of ncbi.datasets run the real upstream generator with a fake API client. The
upstream packages are skipped with importorskip when they are not installed.
"""

from __future__ import annotations

import gzip
import inspect
import json
import random
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.databases import genomepy_functions as gpy
from dnakit.databases import ncbi_datasets_functions as ncbi
from dnakit.databases import pybiomart_functions as pbm
from dnakit.databases import pyensembl_functions as pye
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

MODULES = (pye, gpy, ncbi, pbm)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


FASTA_TEXT = ">chr1 test\nACGT\nGG\n>chr2\nTTTT\n"


# --- common contract -------------------------------------------------------------------------


@pytest.mark.parametrize("module", MODULES, ids=lambda module: module.__name__.rsplit(".", 1)[-1])
def test_public_functions_take_native_keyword_verbose(module: Any) -> None:
    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


# --- pyensembl ----------------------------------------------------------------------------------


def test_read_fasta_records_matches_upstream_and_hand_counts(tmp_path: Path) -> None:
    fasta = pytest.importorskip("pyensembl.fasta")
    path = _write(tmp_path / "toy.fa", FASTA_TEXT)
    result = pye.read_fasta_records(path, verbose=False)
    assert [(r["id"], r["length"], r["sequence"]) for r in result["records"]] == [
        ("chr1", 6, "ACGTGG"),
        ("chr2", 4, "TTTT"),
    ]
    assert result["record_count"] == 2
    assert result["records"] == [
        {"id": key, "length": len(value), "sequence": value}
        for key, value in sorted(fasta.parse_fasta_dictionary(str(path)).items())
    ]


def test_read_fasta_records_reads_gzip_and_bounds_size(tmp_path: Path) -> None:
    pytest.importorskip("pyensembl.fasta")
    path = tmp_path / "toy.fa.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(FASTA_TEXT)
    result = pye.read_fasta_records(path, verbose=False)
    assert [record["sequence"] for record in result["records"]] == ["ACGTGG", "TTTT"]
    with pytest.raises(ConfigurationError, match="max_bytes"):
        pye.read_fasta_records(path, max_bytes=1, verbose=False)


def test_read_first_fasta_id_and_empty_file(tmp_path: Path) -> None:
    pytest.importorskip("pyensembl.fasta")
    path = _write(tmp_path / "toy.fa", FASTA_TEXT)
    assert pye.read_first_fasta_id(path, verbose=False)["identifier"] == "chr1"
    empty = _write(tmp_path / "empty.fa", "")
    assert pye.read_first_fasta_id(empty, verbose=False)["identifier"] is None


def test_merge_genomic_intervals_hand_calculation(tmp_path: Path) -> None:
    common = pytest.importorskip("pyensembl.common")
    pairs = [(5, 9), (1, 3), (4, 4), (12, 15)]
    result = pye.merge_genomic_intervals(pairs, verbose=False)
    # After sorting, (1,3), (4,4) and (5,9) touch or overlap, so they form one block 1..9.
    assert result["intervals"] == [{"start": 1, "end": 9}, {"start": 12, "end": 15}]
    assert result["input_count"] == 4
    assert [list(item) for item in common.merge_intervals(pairs)] == [[1, 9], [12, 15]]
    with pytest.raises(ConfigurationError, match="end before start"):
        pye.merge_genomic_intervals([(9, 5)], verbose=False)


def test_find_nearest_genomic_locus_hand_distance_and_empty_list() -> None:
    pytest.importorskip("pyensembl.locus")
    loci = [
        {"contig": "1", "start": 100, "end": 200, "strand": "+"},
        {"contig": "2", "start": 500, "end": 600, "strand": "-"},
    ]
    result = pye.find_nearest_genomic_locus(300, 310, loci, verbose=False)
    # Interval 300-310 is 100 bases after locus 100-200 and 190 bases before 500-600.
    assert result["distance"] == 100
    assert result["nearest_index"] == 0
    assert result["nearest"] == {"contig": "1", "start": 100, "end": 200, "strand": "+"}
    empty = pye.find_nearest_genomic_locus(300, 310, [], verbose=False)
    assert empty["distance"] is None and empty["nearest"] is None
    with pytest.raises(ConfigurationError, match="strand"):
        pye.find_nearest_genomic_locus(
            1, 2, [{"contig": "1", "start": 1, "end": 2, "strand": "x"}], verbose=False
        )


def test_stable_ensembl_id_strips_version_only() -> None:
    gene_aliases = pytest.importorskip("pyensembl.gene_name_aliases")
    assert pye.stable_ensembl_id("ENSG00000139618.15", verbose=False)["stable_id"] == (
        "ENSG00000139618"
    )
    assert pye.stable_ensembl_id("AT1G01010.1", verbose=False)["stable_id"] == "AT1G01010.1"
    assert gene_aliases.stable_ensembl_gene_id("ENSG00000139618.15") == "ENSG00000139618"


def test_normalize_chromosome_name_rules() -> None:
    pytest.importorskip("pyensembl.normalization")
    assert pye.normalize_chromosome_name("chrx", verbose=False)["chromosome"] == "chrX"
    assert pye.normalize_chromosome_name("x", verbose=False)["chromosome"] == "X"
    assert pye.normalize_chromosome_name(7, verbose=False)["chromosome"] == "7"
    assert pye.normalize_chromosome_name("chrUn_gl000212", verbose=False)["chromosome"] == (
        "chrUn_gl000212"
    )
    with pytest.raises(ConfigurationError):
        pye.normalize_chromosome_name(0, verbose=False)
    with pytest.raises(ConfigurationError):
        pye.normalize_chromosome_name(True, verbose=False)


def test_normalize_strand_symbol_accepts_only_strand_notation() -> None:
    pytest.importorskip("pyensembl.normalization")
    assert pye.normalize_strand_symbol("1", verbose=False)["strand"] == "+"
    assert pye.normalize_strand_symbol(-1, verbose=False)["strand"] == "-"
    with pytest.raises(ConfigurationError):
        pye.normalize_strand_symbol(True, verbose=False)
    with pytest.raises(ConfigurationError):
        pye.normalize_strand_symbol("x", verbose=False)


def test_normalize_gene_aliases_sorts_names_and_ids() -> None:
    pytest.importorskip("pyensembl.gene_name_aliases")
    result = pye.normalize_gene_aliases({"BRCA2": "ENSG1.2", "X": ["b", "a", "a"]}, verbose=False)
    assert result["records"] == [
        {"name": "BRCA2", "gene_ids": ["ENSG1.2"]},
        {"name": "X", "gene_ids": ["a", "b"]},
    ]
    with pytest.raises(ConfigurationError, match="mapping"):
        pye.normalize_gene_aliases(["BRCA2"], verbose=False)


def test_match_identifier_version_cases() -> None:
    pytest.importorskip("pyensembl.versioned_ids")
    assert pye.match_identifier_version("ENSG1.2", {"ENSG1": 2}, verbose=False)["matched"] == (
        "ENSG1"
    )
    assert pye.match_identifier_version("ENSG1", {"ENSG1.2": 2}, verbose=False)["matched"] == (
        "ENSG1.2"
    )
    # A different stable ID is not matched even when it is the only installed form.
    absent = pye.match_identifier_version("ENSG9", {"ENSG1": 2}, verbose=False)
    assert absent["matched"] is None and absent["compared_count"] == 0
    with pytest.raises(ConfigurationError, match="several versions"):
        pye.match_identifier_version("ENSG1", {"ENSG1.1": 1, "ENSG1.2": 2}, verbose=False)
    with pytest.raises(ConfigurationError):
        pye.match_identifier_version("ENSG1.9", {"ENSG1": 2}, verbose=False)


def test_lookup_sequence_version_fallback_matches_upstream() -> None:
    sequence_data = pytest.importorskip("pyensembl.sequence_data")
    sequences = {"ENSP1.3": "MKV", "ENST2": "ACGT"}
    found = pye.lookup_sequence_version_fallback(sequences, "ENST2.5", verbose=False)
    assert found["found"] is True and found["sequence"] == "ACGT"
    assert sequence_data.lookup_sequence_with_version_fallback(sequences, "ENST2.5") == "ACGT"
    literal = pye.lookup_sequence_version_fallback(sequences, "ENSP1.3", verbose=False)
    assert literal["sequence"] == "MKV"
    missing = pye.lookup_sequence_version_fallback(sequences, "NOPE", verbose=False)
    assert missing["found"] is False and missing["sequence"] is None


def test_normalize_genome_reference_name_and_unknown_name() -> None:
    pytest.importorskip("pyensembl.reference_name")
    assert pye.normalize_genome_reference_name("grch38", verbose=False)["reference"] == "GRCh38"
    with pytest.raises(ConfigurationError, match="not found"):
        pye.normalize_genome_reference_name("NoSuchGenome", verbose=False)


def test_ensembl_reference_name_known_releases() -> None:
    reference = pytest.importorskip("pyensembl.reference_name")
    # Ensembl 75 is the last GRCh37 release; release 110 uses GRCh38.
    assert pye.ensembl_reference_name("human", 75, verbose=False)["reference"] == "GRCh37"
    assert pye.ensembl_reference_name("human", 110, verbose=False)["reference"] == "GRCh38"
    assert reference.which_reference("human", 110) == "GRCh38"
    with pytest.raises(ConfigurationError):
        pye.ensembl_reference_name("human", 0, verbose=False)


def test_ensembl_gtf_url_hand_written_path() -> None:
    templates = pytest.importorskip("pyensembl.ensembl_url_templates")
    expected = (
        "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/"
        "Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf.gz"
    )
    result = pye.ensembl_gtf_url(110, "human", verbose=False)
    assert result["url"] == expected
    assert templates.make_gtf_url(110, "human") == expected


def test_ensembl_genome_fasta_url_soft_primary_assembly() -> None:
    templates = pytest.importorskip("pyensembl.ensembl_url_templates")
    expected = (
        "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/"
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
    )
    result = pye.ensembl_genome_fasta_url(
        110, "human", fasta_type="primary_assembly", mask="soft", verbose=False
    )
    assert result["url"] == expected
    assert templates.make_genome_fasta_url(110, "human", "primary_assembly", "soft") == expected
    with pytest.raises(ConfigurationError, match="mask"):
        pye.ensembl_genome_fasta_url(110, "human", mask="weak", verbose=False)


def test_ensembl_sequence_fasta_url_protein_file() -> None:
    templates = pytest.importorskip("pyensembl.ensembl_url_templates")
    expected = (
        "https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/"
        "Mus_musculus.GRCm39.pep.all.fa.gz"
    )
    assert pye.ensembl_sequence_fasta_url(110, "mouse", "pep", verbose=False)["url"] == expected
    assert templates.make_fasta_url(110, "mouse", "pep") == expected
    with pytest.raises(ConfigurationError, match="sequence_type"):
        pye.ensembl_sequence_fasta_url(110, "mouse", "genome", verbose=False)


def test_pyensembl_missing_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path / "toy.fa", FASTA_TEXT)
    monkeypatch.setitem(sys.modules, "pyensembl", None)
    monkeypatch.setitem(sys.modules, "pyensembl.fasta", None)
    with pytest.raises(BackendUnavailableError, match="bio-pyensembl"):
        pye.read_fasta_records(path, verbose=False)


# --- genomepy -----------------------------------------------------------------------------------

GTF_TEXT = (
    'chr1\tsrc\tgene\t10\t100\t.\t+\t.\tgene_id "G1";\n'
    'chr2\tsrc\tgene\t5\t50\t.\t-\t.\tgene_id "G2";\n'
)
BED12_TEXT = "chr1\t0\t100\tname\t0\t+\t0\t100\t0\t1\t100\t0\n"


def test_read_annotation_records_gtf_and_bed12(tmp_path: Path) -> None:
    annotation = pytest.importorskip("genomepy.annotation.utils")
    gtf = _write(tmp_path / "a.gtf", GTF_TEXT)
    result = gpy.read_annotation_records(gtf, verbose=False)
    assert result["file_format"] == "gtf" and result["record_count"] == 2
    assert result["records"][0] == {
        "seqname": "chr1",
        "source": "src",
        "feature": "gene",
        "start": 10,
        "end": 100,
        "score": ".",
        "strand": "+",
        "frame": ".",
        "attribute": 'gene_id "G1";',
    }
    frame = annotation.read_annot(str(gtf))
    assert result["records"] == json.loads(frame.to_json(orient="records"))
    bed = _write(tmp_path / "b.bed", BED12_TEXT)
    bed_result = gpy.read_annotation_records(bed, verbose=False)
    assert bed_result["file_format"] == "bed12"
    assert bed_result["records"][0]["end"] == 100
    with pytest.raises(ConfigurationError, match="GTF or BED12"):
        gpy.read_annotation_records(_write(tmp_path / "x.txt", "a\tb\tc\n"), verbose=False)


def test_write_annotation_records_round_trip_and_no_overwrite(tmp_path: Path) -> None:
    pytest.importorskip("genomepy.annotation.utils")
    record = {
        "seqname": "chr1",
        "source": "s",
        "feature": "gene",
        "start": 1,
        "end": 5,
        "score": ".",
        "strand": "+",
        "frame": ".",
        "attribute": 'gene_id "X";',
    }
    target = tmp_path / "out.gtf"
    result = gpy.write_annotation_records([record], target, verbose=False)
    assert result["record_count"] == 1
    assert target.read_text(encoding="utf-8") == 'chr1\ts\tgene\t1\t5\t.\t+\t.\tgene_id "X";\n'
    assert gpy.read_annotation_records(target, verbose=False)["records"] == [record]
    with pytest.raises(ConfigurationError, match="already exists"):
        gpy.write_annotation_records([record], target, verbose=False)
    with pytest.raises(ConfigurationError, match="exactly the keys"):
        gpy.write_annotation_records([{"seqname": "chr1"}], tmp_path / "bad.gtf", verbose=False)


def test_filter_annotation_contigs_matches_upstream(tmp_path: Path) -> None:
    annotation = pytest.importorskip("genomepy.annotation")
    pandas = pytest.importorskip("pandas")
    records = gpy.read_annotation_records(_write(tmp_path / "a.gtf", GTF_TEXT), verbose=False)[
        "records"
    ]
    kept = gpy.filter_annotation_contigs(records, "chr2", verbose=False)
    assert [record["seqname"] for record in kept["records"]] == ["chr2"]
    assert kept["kept_count"] == 1 and kept["removed_count"] == 1
    inverted = gpy.filter_annotation_contigs(records, "chr2", invert_match=True, verbose=False)
    assert [record["seqname"] for record in inverted["records"]] == ["chr1"]
    frame = annotation.filter_regex(pandas.DataFrame(records), "chr2", invert_match=False)
    assert list(frame["seqname"]) == ["chr2"]
    with pytest.raises(ConfigurationError, match="regular expression"):
        gpy.filter_annotation_contigs(records, "chr(", verbose=False)
    numeric = [dict(records[0], seqname=1)]
    with pytest.raises(ConfigurationError, match="must be a string to filter"):
        gpy.filter_annotation_contigs(numeric, "1", verbose=False)


def test_write_gap_bed_hand_coordinates(tmp_path: Path) -> None:
    pytest.importorskip("genomepy.genome")
    fasta = _write(
        tmp_path / "n.fa",
        ">chr1\n" + "ACGT" * 10 + "NNNN" + "ACGT" * 5 + "\n>chr2 desc\nGGGGCCCCAAAATTTT\n"
        ">scaf_9\nNNNNACGTAC\n",
    )
    result = gpy.write_gap_bed(fasta, tmp_path / "gaps.bed", verbose=False)
    # The N run after 40 bases of sequence is 0-based [40, 44); scaf_9 starts with a run [0, 4).
    assert (tmp_path / "gaps.bed").read_text(encoding="utf-8") == "chr1\t40\t44\nscaf_9\t0\t4\n"
    assert result["gap_count"] == 2


def test_write_fasta_sizes_hand_lengths(tmp_path: Path) -> None:
    pytest.importorskip("genomepy.genome")
    fasta = _write(
        tmp_path / "n.fa",
        ">chr1\n" + "ACGT" * 10 + "NNNN" + "ACGT" * 5 + "\n>chr2 desc\nGGGGCCCCAAAATTTT\n"
        ">scaf_9\nNNNNACGTAC\n",
    )
    result = gpy.write_fasta_sizes(fasta, tmp_path / "n.sizes", verbose=False)
    assert (tmp_path / "n.sizes").read_text(encoding="utf-8") == "chr1\t64\nchr2\t16\nscaf_9\t10\n"
    assert result["sequence_count"] == 3


def test_filter_fasta_records_removes_matching_headers(tmp_path: Path) -> None:
    files = pytest.importorskip("genomepy.files")
    fasta = _write(tmp_path / "n.fa", ">chr1\nACGT\n>chr2 desc\nGGGG\n>scaf_9\nNNNN\n")
    result = gpy.filter_fasta_records(fasta, tmp_path / "f.fa", pattern="chr", verbose=False)
    assert result["removed"] == ["scaf_9"]
    assert result["kept_count"] == 2
    assert (tmp_path / "f.fa").read_text(encoding="utf-8") == ">chr1\nACGT\n>chr2 desc\nGGGG\n"
    upstream = files.filter_fasta(str(fasta), str(tmp_path / "g.fa"), regex="chr")
    assert upstream == ["scaf_9"]
    with pytest.raises(ConfigurationError, match="differ"):
        gpy.filter_fasta_records(fasta, fasta, verbose=False)


def _random_fasta(path: Path) -> Path:
    rng = random.Random(5)
    chr1 = "".join(rng.choice("ACGT") for _ in range(3000))
    chr2 = "".join(rng.choice("ACGT") for _ in range(2500))
    return _write(path, f">chr1\n{chr1}\n>chr2\n{chr2}\n>scaf_9\nNNNNACGTAC\n")


def test_sample_random_regions_seeded_and_matches_upstream(tmp_path: Path) -> None:
    genomepy = pytest.importorskip("genomepy")
    fasta = _random_fasta(tmp_path / "r.fa")
    before = random.getstate()
    first = gpy.sample_random_regions(fasta, n=3, length=20, seed=7, verbose=False)
    assert random.getstate() == before, "global random state must be restored"
    second = gpy.sample_random_regions(fasta, n=3, length=20, seed=7, verbose=False)
    assert first["regions"] == second["regions"]
    assert all(region["end"] - region["start"] == 20 for region in first["regions"])
    assert {region["chrom"] for region in first["regions"]} <= {"chr1", "chr2"}
    genome = genomepy.Genome(str(fasta))
    random.seed(7)
    direct = genome.get_random_sequences(n=3, length=20, outtype="list")
    assert first["regions"] == [
        {"chrom": chrom, "start": start, "end": end} for chrom, start, end in direct
    ]
    only_chr2 = gpy.sample_random_regions(
        fasta, n=2, length=20, chroms=["chr2"], seed=1, verbose=False
    )
    assert {region["chrom"] for region in only_chr2["regions"]} == {"chr2"}
    with pytest.raises(ConfigurationError, match="not in the FASTA"):
        gpy.sample_random_regions(fasta, chroms=["nope"], verbose=False)
    short = _write(tmp_path / "short.fa", ">chr1\nACGTACGT\n")
    with pytest.raises(ConfigurationError, match="cannot sample"):
        gpy.sample_random_regions(short, n=1, length=4, verbose=False)


def test_genomepy_missing_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fasta = _write(tmp_path / "n.fa", ">chr1\nACGT\n")
    monkeypatch.setitem(sys.modules, "genomepy.genome", None)
    with pytest.raises(BackendUnavailableError, match="bio-genomepy"):
        gpy.write_gap_bed(fasta, tmp_path / "gaps.bed", verbose=False)


# --- pybiomart ----------------------------------------------------------------------------------


def test_pybiomart_module_has_no_wrapped_functions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pybiomart = pytest.importorskip("pybiomart")
    monkeypatch.chdir(tmp_path)
    assert pbm.__all__ == ()
    public = [
        name
        for name, value in vars(pybiomart).items()
        if not name.startswith("_")
        and callable(value)
        and not inspect.isclass(value)
        and getattr(value, "__module__", "").startswith("pybiomart")
    ]
    assert public == []
    assert all(inspect.isclass(getattr(pybiomart, name)) for name in ("Server", "Mart", "Dataset"))


# --- ncbi.datasets ------------------------------------------------------------------------------


class _Reply:
    def __init__(self, assemblies: list[Any]) -> None:
        self.assemblies = assemblies
        self.next_page_token = None


class _FakeClient:
    def __enter__(self) -> _FakeClient:
        return self

    def __exit__(self, *exc: object) -> bool:
        return False


def _descriptor(accession: str) -> Any:
    model = pytest.importorskip("ncbi.datasets.openapi.model.v1_assembly_dataset_descriptor")
    return model.V1AssemblyDatasetDescriptor(accession=accession, current_accession=accession)


def _patch_api(monkeypatch: pytest.MonkeyPatch, api_class: type) -> None:
    genome = pytest.importorskip("ncbi.datasets.metadata.genome")
    monkeypatch.setattr(genome, "DatasetsApiClient", _FakeClient)
    monkeypatch.setattr(genome, "DatasetsGenomeApi", api_class)


def test_fetch_ncbi_assembly_metadata_runs_upstream_with_fake_client(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("ncbi.datasets.metadata.genome")
    seen: list[object] = []

    class FakeApi:
        def __init__(self, client: object) -> None:
            pass

        def assembly_descriptors_by_accessions(
            self, accessions: list[str], page_size: int, page_token: str
        ) -> _Reply:
            seen.append((tuple(accessions), page_size, page_token))
            return _Reply([_descriptor(item) for item in accessions])

    _patch_api(monkeypatch, FakeApi)
    result = ncbi.fetch_ncbi_assembly_metadata(
        ["GCF_000001405.40", "GCA_000001635.9"], verbose=False
    )
    assert result["record_count"] == 2
    assert [record["accession"] for record in result["records"]] == [
        "GCF_000001405.40",
        "GCA_000001635.9",
    ]
    assert result["messages"] == []
    assert [item[0] for item in seen] == [("GCF_000001405.40", "GCA_000001635.9")]


def test_fetch_ncbi_assembly_metadata_reports_empty_reply_and_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("ncbi.datasets.metadata.genome")

    class EmptyApi:
        def __init__(self, client: object) -> None:
            pass

        def assembly_descriptors_by_accessions(self, accessions: Any, **kwargs: Any) -> _Reply:
            return _Reply([])

    _patch_api(monkeypatch, EmptyApi)
    empty = ncbi.fetch_ncbi_assembly_metadata(["GCF_000001405.40"], verbose=False)
    assert empty["record_count"] == 0
    assert any(line.startswith("No assemblies found") for line in empty["messages"])

    class TwoApi(EmptyApi):
        def assembly_descriptors_by_accessions(self, accessions: Any, **kwargs: Any) -> _Reply:
            return _Reply([_descriptor("GCF_000001405.40"), _descriptor("GCF_000001405.39")])

    _patch_api(monkeypatch, TwoApi)
    with pytest.raises(ConfigurationError, match="at least 1"):
        ncbi.fetch_ncbi_assembly_metadata([], verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        ncbi.fetch_ncbi_assembly_metadata(["GCF_000001405.40"], max_records=1, verbose=False)
    with pytest.raises(ConfigurationError, match="unsupported format"):
        ncbi.fetch_ncbi_assembly_metadata(["bogus"], verbose=False)


def test_fetch_ncbi_assemblies_by_taxon_and_bioproject_use_mocked_client(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("ncbi.datasets.metadata.genome")
    calls: dict[str, object] = {}

    class FakeApi:
        def __init__(self, client: object) -> None:
            pass

        def assembly_descriptors_by_taxon(self, taxon: str, **kwargs: Any) -> _Reply:
            calls["taxon"] = taxon
            return _Reply([_descriptor("GCF_000001405.40")])

        def assembly_descriptors_by_bioproject(
            self, accessions: list[str], **kwargs: Any
        ) -> _Reply:
            calls["bioproject"] = list(accessions)
            return _Reply([_descriptor("GCF_000001405.40")])

    _patch_api(monkeypatch, FakeApi)
    by_taxon = ncbi.fetch_ncbi_assemblies_by_taxon(9606, verbose=False)
    assert calls["taxon"] == "9606"
    assert by_taxon["records"][0]["accession"] == "GCF_000001405.40"
    by_project = ncbi.fetch_ncbi_assemblies_by_bioproject(["PRJNA168"], verbose=False)
    assert calls["bioproject"] == ["PRJNA168"]
    assert by_project["record_count"] == 1
    with pytest.raises(ConfigurationError, match="taxon"):
        ncbi.fetch_ncbi_assemblies_by_taxon(True, verbose=False)
    with pytest.raises(ConfigurationError, match="bioproject"):
        ncbi.fetch_ncbi_assemblies_by_bioproject(["GCF_000001405.40"], verbose=False)


def test_select_ncbi_gene_fields_matches_upstream() -> None:
    gene_module = pytest.importorskip("ncbi.datasets.metadata.gene")
    model = pytest.importorskip("ncbi.datasets.openapi.model.v1_gene_descriptor")
    gene = {"gene_id": "7157", "symbol": "TP53"}
    assert ncbi.select_ncbi_gene_fields(gene, ["symbol"], verbose=False)["values"] == {
        "symbol": "TP53"
    }
    both = ncbi.select_ncbi_gene_fields(gene, ["gene_id", "symbol"], verbose=False)
    upstream = gene_module.gene_values_by_fields(
        model.V1GeneDescriptor(gene_id="7157", symbol="TP53"), ["gene_id", "symbol"]
    )
    assert both["values"] == upstream == {"gene_id": "7157", "symbol": "TP53"}
    with pytest.raises(ConfigurationError, match="subset"):
        ncbi.select_ncbi_gene_fields(gene, ["description"], verbose=False)
    with pytest.raises(ConfigurationError, match="gene_id"):
        ncbi.select_ncbi_gene_fields({"symbol": "TP53"}, ["symbol"], verbose=False)


def test_select_ncbi_assembly_fields_projection_and_full_record() -> None:
    pytest.importorskip("ncbi.datasets.metadata.genome")
    record = {"accession": "GCF_000001405.40", "organism": {"tax_id": 9606}}
    picked = ncbi.select_ncbi_assembly_fields(record, ["accession", "organism"], verbose=False)
    assert picked["values"] == {"accession": "GCF_000001405.40", "organism": {"tax_id": 9606}}
    only = ncbi.select_ncbi_assembly_fields(record, ["organism"], verbose=False)
    assert only["values"] == {"organism": {"tax_id": 9606}}
    assert only["missing_fields"] == []
    partial = ncbi.select_ncbi_assembly_fields(record, ["accession", "taxid"], verbose=False)
    assert partial["missing_fields"] == ["taxid"]
    full = ncbi.select_ncbi_assembly_fields(record, verbose=False)
    assert full["values"] == record


def test_ncbi_missing_package_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "ncbi.datasets.metadata.genome", None)
    with pytest.raises(BackendUnavailableError, match="bio-ncbi-datasets"):
        ncbi.fetch_ncbi_assembly_metadata(["GCF_000001405.40"], verbose=False)
