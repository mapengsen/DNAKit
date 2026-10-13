"""Contract tests for gene-level database clients over optional upstream packages.

Local operations (pyensembl, genomepy, biothings_client, gget.mutate) run the real upstream
code on tiny files written into tmp_path. Network operations replace only the upstream
client call and check the parsing and validation. One live test runs only when
DNAKIT_LIVE_TESTS=1 is set.
"""

from __future__ import annotations

import importlib
import io
import os
import sys
import types
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.databases import gene_backends as gb
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

# Periodic 400-base contig: S[p - 1] is the 1-based base p, so slices are easy to check by hand.
CONTIG = ("ACGTTGCA" * 50)[:400]


def _gtf_line(feature: str, start: int, end: int, strand: str, attributes: str) -> str:
    fields = ["chr1", "test", feature, str(start), str(end), ".", strand, "."]
    return "\t".join(fields) + "\t" + attributes


GENE_A = 'gene_id "ENSG00000000001"; gene_name "GENEA"'
GENE_B = 'gene_id "ENSG00000000002"; gene_name "GENEB"'
TX_A = f'{GENE_A}; transcript_id "ENST00000000011"; transcript_name "GENEA-201"'
TX_B = f'{GENE_B}; transcript_id "ENST00000000022"; transcript_name "GENEB-201"'
GTF_TEXT = (
    "##description: tiny Ensembl-style annotation for DNAKit examples\n"
    + _gtf_line("gene", 11, 200, "+", f'{GENE_A}; gene_biotype "protein_coding";')
    + "\n"
    + _gtf_line("transcript", 11, 200, "+", f'{TX_A}; transcript_biotype "protein_coding";')
    + "\n"
    + _gtf_line("exon", 11, 60, "+", f'{TX_A}; exon_number "1"; exon_id "ENSE00000000001";')
    + "\n"
    + _gtf_line("exon", 121, 200, "+", f'{TX_A}; exon_number "2"; exon_id "ENSE00000000002";')
    + "\n"
    + _gtf_line("gene", 301, 380, "-", f'{GENE_B}; gene_biotype "lncRNA";')
    + "\n"
    + _gtf_line("transcript", 301, 380, "-", f'{TX_B}; transcript_biotype "lncRNA";')
    + "\n"
    + _gtf_line("exon", 301, 380, "-", f'{TX_B}; exon_number "1"; exon_id "ENSE00000000003";')
    + "\n"
)

BED_TEXT = "chr1\t10\t20\tr1\t0\t+\nchr1\t300\t310\tr2\t0\t-\n"

VCF_TEXT = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\trs1\tA\tG,T\t.\tPASS\t.
chr1\t20\t.\tC\t.\t.\t.\t.
chrUn\t5\t.\tC\tT\t.\t.\t.
chr2\t30\t.\tGA\tG\t.\t.\t.
"""


def _import_or_skip(name: str) -> Any:
    """Skip when a package is missing or broken; importorskip re-raises broken installs."""

    try:
        return importlib.import_module(name)
    except ImportError as exc:
        pytest.skip(f"{name} could not be imported: {exc}")


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _fasta(text: str) -> str:
    lines = [">chr1 example contig"] + [text[i : i + 60] for i in range(0, len(text), 60)]
    return "\n".join(lines) + "\n"


@pytest.fixture()
def local_files(tmp_path: Path) -> dict[str, Path]:
    files = {
        "gtf": tmp_path / "genes.gtf",
        "fasta": tmp_path / "genome.fa",
        "bed": tmp_path / "regions.bed",
        "vcf": tmp_path / "sample.vcf",
        "cache": tmp_path / "ensembl_cache",
    }
    files["gtf"].write_text(GTF_TEXT, encoding="utf-8")
    files["fasta"].write_text(_fasta(CONTIG), encoding="utf-8")
    files["bed"].write_text(BED_TEXT, encoding="utf-8")
    files["vcf"].write_text(VCF_TEXT, encoding="utf-8")
    return files


# --- pyensembl: local GTF and genome FASTA -------------------------------------------------


def test_index_gene_models_counts_match_direct_pyensembl(local_files: dict[str, Path]) -> None:
    pyensembl = pytest.importorskip("pyensembl")
    result = gb.index_gene_models(local_files["gtf"], cache_dir=local_files["cache"], verbose=False)
    direct = pyensembl.Genome(
        reference_name="dnakit_direct",
        annotation_name="direct_check",
        gtf_path_or_url=str(local_files["gtf"]),
        cache_directory_path=str(local_files["cache"]),
    )
    direct.index()
    assert result["contigs"] == len(direct.contigs()) == 1
    assert result["genes"] == len(direct.gene_ids()) == 2
    assert result["transcripts"] == len(direct.transcript_ids()) == 2
    assert result["exons"] == len(direct.exon_ids()) == 3
    assert result["provenance"]["package"].startswith("pyensembl ")


def test_query_genes_by_name_returns_structure_and_empty_for_unknown(
    local_files: dict[str, Path],
) -> None:
    pytest.importorskip("pyensembl")
    result = gb.query_genes_by_name(
        local_files["gtf"], ["GENEA", "NOPE"], cache_dir=local_files["cache"], verbose=False
    )
    gene = result["results"]["GENEA"][0]
    assert gene == {
        "gene_id": "ENSG00000000001",
        "gene_name": "GENEA",
        "contig": "chr1",
        "start": 11,
        "end": 200,
        "strand": "+",
        "biotype": "protein_coding",
        "transcript_ids": ["ENST00000000011"],
    }
    assert result["results"]["NOPE"] == []
    assert result["matches"] == 1


def test_query_genes_at_locus_uses_inclusive_overlap_and_strand(
    local_files: dict[str, Path],
) -> None:
    pytest.importorskip("pyensembl")
    cache = local_files["cache"]
    inside = gb.query_genes_at_locus(local_files["gtf"], "chr1", 15, cache_dir=cache, verbose=False)
    assert [gene["gene_name"] for gene in inside["genes"]] == ["GENEA"]
    both = gb.query_genes_at_locus(
        local_files["gtf"], "chr1", 190, 310, cache_dir=cache, verbose=False
    )
    assert [gene["gene_name"] for gene in both["genes"]] == ["GENEA", "GENEB"]
    # GENEB spans 301-380 on the minus strand, so a plus-strand query finds nothing there.
    plus = gb.query_genes_at_locus(
        local_files["gtf"], "chr1", 305, strand="+", cache_dir=cache, verbose=False
    )
    assert plus["count"] == 0
    with pytest.raises(ConfigurationError, match="not in the annotation"):
        gb.query_genes_at_locus(local_files["gtf"], "chrZ", 1, cache_dir=cache, verbose=False)
    with pytest.raises(ConfigurationError, match="end"):
        gb.query_genes_at_locus(local_files["gtf"], "chr1", 50, 10, cache_dir=cache, verbose=False)


def test_fetch_gene_sequence_matches_hand_slices(local_files: dict[str, Path]) -> None:
    pytest.importorskip("pyensembl")
    cache = local_files["cache"]
    forward = gb.fetch_gene_sequence(
        local_files["gtf"], local_files["fasta"], "ENSG00000000001", cache_dir=cache, verbose=False
    )
    # GENEA is 11-200 on the plus strand: 190 bases, the contig's 1-based bases 11..200.
    assert forward["length"] == 190
    assert forward["sequence"] == CONTIG[10:200]
    reverse = gb.fetch_gene_sequence(
        local_files["gtf"], local_files["fasta"], "ENSG00000000002", cache_dir=cache, verbose=False
    )
    # GENEB is 301-380 on the minus strand: the reverse complement of the locus.
    assert reverse["length"] == 80
    assert reverse["sequence"] == _rc(CONTIG[300:380])
    with pytest.raises(ConfigurationError, match="not in the annotation"):
        gb.fetch_gene_sequence(
            local_files["gtf"], local_files["fasta"], "ENSG_MISSING", cache_dir=cache, verbose=False
        )


def test_fetch_gene_sequence_enforces_max_length(local_files: dict[str, Path]) -> None:
    pytest.importorskip("pyensembl")
    with pytest.raises(ConfigurationError, match="max_length=100"):
        gb.fetch_gene_sequence(
            local_files["gtf"],
            local_files["fasta"],
            "ENSG00000000001",
            cache_dir=local_files["cache"],
            max_length=100,
            verbose=False,
        )


# --- genomepy: local FASTA, spliced sequences and BED extraction ---------------------------


def test_read_genome_summary_reports_contigs_and_missing_metadata(
    local_files: dict[str, Path],
) -> None:
    pytest.importorskip("genomepy")
    result = gb.read_genome_summary(local_files["fasta"], verbose=False)
    assert result["name"] == "genome"
    assert result["contig_count"] == 1
    assert result["total_length"] == 400
    assert result["contigs"] == [{"name": "chr1", "length": 400, "n_bases": 0}]
    assert result["tax_id"] is None and result["assembly_accession"] is None


def test_fetch_spliced_sequence_joins_exons_like_direct_get_seq(
    local_files: dict[str, Path],
) -> None:
    genomepy = pytest.importorskip("genomepy")
    result = gb.fetch_spliced_sequence(
        local_files["fasta"], "chr1", [[11, 60], [121, 200]], verbose=False
    )
    # Exon 1 is 50 bases and exon 2 is 80 bases: 130 bases in total.
    assert result["length"] == 130
    assert result["sequence"] == CONTIG[10:60] + CONTIG[120:200]
    direct = genomepy.Genome(str(local_files["fasta"]))
    upstream = direct.get_seq("chr1", 11, 60).seq + direct.get_seq("chr1", 121, 200).seq
    assert result["sequence"] == upstream
    reverse = gb.fetch_spliced_sequence(
        local_files["fasta"], "chr1", [[301, 380]], reverse_complement=True, verbose=False
    )
    assert reverse["sequence"] == _rc(CONTIG[300:380])


def test_fetch_spliced_sequence_rejects_invalid_intervals(local_files: dict[str, Path]) -> None:
    pytest.importorskip("genomepy")
    fasta = local_files["fasta"]
    with pytest.raises(ConfigurationError, match="start"):
        gb.fetch_spliced_sequence(fasta, "chr1", [[0, 10]], verbose=False)
    with pytest.raises(ConfigurationError, match="start"):
        gb.fetch_spliced_sequence(fasta, "chr1", [[True, 10]], verbose=False)
    with pytest.raises(ConfigurationError, match="exceeds the contig"):
        gb.fetch_spliced_sequence(fasta, "chr1", [[1, 401]], verbose=False)
    with pytest.raises(ConfigurationError, match="pair"):
        gb.fetch_spliced_sequence(fasta, "chr1", [[1, 2, 3]], verbose=False)
    with pytest.raises(ConfigurationError, match="not in the FASTA"):
        gb.fetch_spliced_sequence(fasta, "chrX", [[1, 2]], verbose=False)
    with pytest.raises(ConfigurationError, match="max_length"):
        gb.fetch_spliced_sequence(fasta, "chr1", [[1, 100]], max_length=50, verbose=False)


def test_write_bed_sequences_writes_stranded_records_and_never_overwrites(
    local_files: dict[str, Path],
) -> None:
    pytest.importorskip("genomepy")
    output = local_files["fasta"].parent / "regions.fa"
    result = gb.write_bed_sequences(
        local_files["fasta"], local_files["bed"], output, stranded=True, verbose=False
    )
    assert result["regions"] == 2 and result["records"] == 2 and result["skipped"] == 0
    # BED is 0-based half-open: r1 covers 1-based 11..20; r2 is minus-strand 301..310.
    text = output.read_text(encoding="utf-8").splitlines()
    assert text[1] == CONTIG[10:20]
    assert text[3] == _rc(CONTIG[300:310])
    with pytest.raises(ConfigurationError, match="already exists"):
        gb.write_bed_sequences(
            local_files["fasta"], local_files["bed"], output, stranded=True, verbose=False
        )


def test_install_genome_assembly_calls_genomepy_with_explicit_arguments(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    genomepy = pytest.importorskip("genomepy")
    calls: dict[str, Any] = {}

    def fake_install(name: str, **kwargs: Any) -> SimpleNamespace:
        calls.update(name=name, **kwargs)
        return SimpleNamespace(genome_file=str(tmp_path / "hg38.fa"), annotation_gtf_file=None)

    monkeypatch.setattr(genomepy, "install_genome", fake_install)
    result = gb.install_genome_assembly(
        "hg38", provider="ucsc", genomes_dir=str(tmp_path), threads=2, verbose=False
    )
    assert calls["name"] == "hg38"
    assert calls["provider"] == "ucsc" and calls["threads"] == 2
    assert calls["annotation"] is False
    assert result["fasta"].endswith("hg38.fa") and result["annotation_gtf"] is None
    with pytest.raises(ConfigurationError, match="threads"):
        gb.install_genome_assembly("hg38", threads=True, verbose=False)


# --- biothings_client: HGVS identifiers (local) -----------------------------------------


def test_format_hgvs_id_matches_upstream_and_hand_values() -> None:
    biothings = pytest.importorskip("biothings_client")
    upstream = biothings.get_client("variant")
    assert gb.format_hgvs_id("chr9", 107620835, "G", "A", verbose=False)["hgvs_id"] == (
        "chr9:g.107620835G>A"
    )
    # Deletion of 9 bases after a shared first base: 8271-8279 (1-based).
    deletion = gb.format_hgvs_id("MT", 8270, "CACCCCCTCT", "C", verbose=False)
    assert deletion["hgvs_id"] == "chrMT:g.8271_8279del"
    assert deletion["hgvs_id"] == upstream.format_hgvs("MT", 8270, "CACCCCCTCT", "C")
    assert gb.format_hgvs_id("2", 17142, "G", "GA", verbose=False)["hgvs_id"] == (
        "chr2:g.17142_17143insA"
    )
    with pytest.raises(ConfigurationError, match="human chromosome"):
        gb.format_hgvs_id("chrUn_KI270", 5, "A", "G", verbose=False)
    with pytest.raises(ConfigurationError, match="must differ"):
        gb.format_hgvs_id("1", 5, "A", "A", verbose=False)
    with pytest.raises(ConfigurationError, match="A, C, G, T, N"):
        gb.format_hgvs_id("1", 5, "A", "R", verbose=False)


def test_hgvs_ids_from_vcf_skips_unsupported_records_and_matches_upstream(
    local_files: dict[str, Path],
) -> None:
    biothings = pytest.importorskip("biothings_client")
    result = gb.hgvs_ids_from_vcf(local_files["vcf"], verbose=False)
    assert result["records"] == 4
    assert result["hgvs_ids"] == ["chr1:g.10A>G", "chr1:g.10A>T", "chr2:g.31del"]
    assert result["skipped"] == {
        "unsupported_contig": 1,
        "non_sequence_allele": 1,
        "same_as_ref": 0,
    }
    # On records with valid alleles, the output equals the upstream VCF generator.
    valid = "#CHROM\tPOS\tID\tREF\tALT\n1\t10\trs1\tA\tG,T\nchr2\t30\t.\tGA\tG\n"
    upstream = list(biothings.get_client("variant").get_hgvs_from_vcf(io.StringIO(valid)))
    assert upstream == ["chr1:g.10A>G", "chr1:g.10A>T", "chr2:g.31del"]
    path = local_files["vcf"].parent / "valid.vcf"
    path.write_text(valid, encoding="utf-8")
    assert gb.hgvs_ids_from_vcf(path, verbose=False)["hgvs_ids"] == upstream
    with pytest.raises(ConfigurationError, match="max_records=2"):
        gb.hgvs_ids_from_vcf(local_files["vcf"], max_records=2, verbose=False)


# --- gget.mutate: local cDNA mutation windows ----------------------------------------------


@pytest.mark.filterwarnings("ignore:This pattern is interpreted as a regular expression")
def test_apply_cdna_mutations_applies_valid_and_drops_wrong_base() -> None:
    gget = pytest.importorskip("gget")
    direct = gget.mutate([CONTIG[:40]], ["c.2C>T"], k=3, merge_identical=False, verbose=False)
    result = gb.apply_cdna_mutations(
        {"s1": CONTIG[:40], "s2": CONTIG[:40]},
        {"s1": "c.10A>G", "s2": "c.2C>T"},
        flank=3,
        verbose=False,
    )
    by_id = {row["seq_id"]: row for row in result["mutations"]}
    # Base 2 of the contig is C, so c.2C>T is applied: window 1..5 with base 2 changed to T.
    assert direct == ["ATGTT"]
    assert by_id["s2"] == {
        "seq_id": "s2",
        "mutation": "c.2C>T",
        "status": "applied",
        "mutant_window": "ATGTT",
    }
    # Base 10 is C, not A, so the mutation is dropped.
    assert by_id["s1"]["status"] == "dropped" and by_id["s1"]["mutant_window"] is None
    assert result["applied"] == 1


def test_apply_cdna_mutations_validates_inputs() -> None:
    pytest.importorskip("gget")
    with pytest.raises(ConfigurationError, match="unknown sequence ID"):
        gb.apply_cdna_mutations({"s1": "ACGT"}, {"s9": "c.1A>G"}, verbose=False)
    with pytest.raises(ConfigurationError, match=r"c\. notation"):
        gb.apply_cdna_mutations({"s1": "ACGT"}, {"s1": "p.Met1Val"}, verbose=False)
    with pytest.raises(ConfigurationError, match="A, C, G, T, N"):
        gb.apply_cdna_mutations({"s1": "ACGR"}, {"s1": "c.1A>G"}, verbose=False)


# --- network clients: upstream call replaced, parsing and validation tested -----------------


def test_search_ensembl_genes_parses_results_and_passes_arguments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gget = pytest.importorskip("gget")
    seen: dict[str, Any] = {}

    def fake_search(words: list[str], **kwargs: Any) -> list[dict[str, Any]]:
        seen.update(words=words, **kwargs)
        return [{"ensembl_id": "ENSG00000123374", "gene_name": "CDK2", "biotype": "protein_coding"}]

    monkeypatch.setattr(gget, "search", fake_search)
    result = gb.search_ensembl_genes("CDK2", max_results=5, verbose=False)
    assert seen["words"] == ["CDK2"] and seen["limit"] == 5 and seen["json"] is True
    assert seen["species"] == "homo_sapiens"
    assert result["count"] == 1 and result["results"][0]["gene_name"] == "CDK2"
    with pytest.raises(ConfigurationError, match="max_results"):
        gb.search_ensembl_genes("CDK2", max_results=0, verbose=False)
    with pytest.raises(ConfigurationError, match="andor"):
        gb.search_ensembl_genes("CDK2", andor="xor", verbose=False)


def test_fetch_ensembl_gene_info_maps_requested_ids_and_missing_as_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gget = pytest.importorskip("gget")
    seen: dict[str, Any] = {}

    def fake_info(ids: list[str], **kwargs: Any) -> dict[str, Any]:
        seen.update(ids=ids, **kwargs)
        return {"ENSG00000123374": {"ensembl_gene_name": "CDK2", "start": 55966781}}

    monkeypatch.setattr(gget, "info", fake_info)
    result = gb.fetch_ensembl_gene_info(["ENSG00000123374", "ENSG00000000000"], verbose=False)
    assert seen["ncbi"] is False and seen["uniprot"] is False  # external sources off by default
    assert result["requested"] == 2 and result["found"] == 1
    assert result["records"]["ENSG00000000000"] is None
    assert result["records"]["ENSG00000123374"]["ensembl_gene_name"] == "CDK2"
    with pytest.raises(ConfigurationError, match="ensembl_ids"):
        gb.fetch_ensembl_gene_info(["bad id with spaces"], verbose=False)


def test_fetch_ensembl_sequences_pairs_headers_and_checks_layout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gget = pytest.importorskip("gget")
    raw = [
        ">ENSG00000123374 chromosome:GRCh38:12:55966781:55972789:1",
        "ATGGAGAACTTCCAAAAG",
    ]
    monkeypatch.setattr(gget, "seq", lambda ids, **kwargs: list(raw))
    result = gb.fetch_ensembl_sequences(["ENSG00000123374"], verbose=False)
    assert result["count"] == 1
    assert result["records"][0] == {
        "id": "ENSG00000123374",
        "description": "chromosome:GRCh38:12:55966781:55972789:1",
        "sequence": "ATGGAGAACTTCCAAAAG",
    }
    monkeypatch.setattr(gget, "seq", lambda ids, **kwargs: raw[:1])
    with pytest.raises(BackendExecutionError, match="without a sequence"):
        gb.fetch_ensembl_sequences(["ENSG00000123374"], verbose=False)
    monkeypatch.setattr(gget, "seq", lambda ids, **kwargs: ["ACGT", "ACGT"])
    with pytest.raises(BackendExecutionError, match="unexpected FASTA layout"):
        gb.fetch_ensembl_sequences(["ENSG00000123374"], verbose=False)


def test_upstream_failures_are_reported_as_backend_execution_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gget = pytest.importorskip("gget")

    def failing(*args: Any, **kwargs: Any) -> None:
        raise RuntimeError("service unavailable")

    monkeypatch.setattr(gget, "search", failing)
    with pytest.raises(BackendExecutionError, match="gget search failed: service unavailable"):
        gb.search_ensembl_genes("CDK2", verbose=False)


def test_mygene_wrappers_parse_hits_records_and_batches(monkeypatch: pytest.MonkeyPatch) -> None:
    mygene = pytest.importorskip("mygene")
    monkeypatch.setattr(
        mygene.MyGeneInfo,
        "query",
        lambda self, q, **kwargs: {"total": 1, "hits": [{"symbol": "CDK2", "_id": "1017"}]},
    )
    monkeypatch.setattr(
        mygene.MyGeneInfo,
        "getgene",
        lambda self, _id, **kwargs: {"symbol": "CDK2", "_id": _id} if _id == "1017" else None,
    )
    monkeypatch.setattr(
        mygene.MyGeneInfo,
        "querymany",
        lambda self, terms, **kwargs: [{"query": t, "symbol": "X"} for t in terms],
    )
    hits = gb.query_mygene_genes("symbol:CDK2", size=1, verbose=False)
    assert hits["total"] == 1 and hits["hits"][0]["symbol"] == "CDK2"
    found = gb.fetch_mygene_gene("1017", verbose=False)
    assert found["found"] is True and found["gene"]["symbol"] == "CDK2"
    missing = gb.fetch_mygene_gene("999999999", verbose=False)
    assert missing["found"] is False and missing["gene"] is None
    mapped = gb.map_gene_identifiers(["1017", "7157"], verbose=False)
    assert [row["query"] for row in mapped["results"]] == ["1017", "7157"]
    with pytest.raises(ConfigurationError, match="size"):
        gb.query_mygene_genes("symbol:CDK2", size=0, verbose=False)


def test_myvariant_wrappers_parse_results_and_check_hgvs_format(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    myvariant = pytest.importorskip("myvariant")
    seen: dict[str, Any] = {}

    def fake_getvariants(self: Any, ids: list[str], **kwargs: Any) -> list[dict[str, Any]]:
        seen.update(ids=ids, **kwargs)
        return [{"query": ids[0], "_id": ids[0], "cadd": {"phred": 19.22}}]

    monkeypatch.setattr(myvariant.MyVariantInfo, "getvariants", fake_getvariants)
    monkeypatch.setattr(
        myvariant.MyVariantInfo,
        "query",
        lambda self, q, **kwargs: {"total": 1, "hits": [{"_id": "chr1:g.1A>G"}]},
    )
    result = gb.fetch_variant_annotations(["chr9:g.107620835G>A"], assembly="hg19", verbose=False)
    assert seen["assembly"] == "hg19"
    assert result["results"][0]["cadd"]["phred"] == 19.22
    query = gb.query_variant_annotations("dbsnp.rsid:rs1", verbose=False)
    assert query["total"] == 1 and query["hits"][0]["_id"] == "chr1:g.1A>G"
    with pytest.raises(ConfigurationError, match="hgvs_ids"):
        gb.fetch_variant_annotations(["9:107620835"], verbose=False)
    with pytest.raises(ConfigurationError, match="assembly"):
        gb.fetch_variant_annotations(["chr9:g.107620835G>A"], assembly="hg17", verbose=False)


def test_biomart_wrappers_parse_frames_and_validate_names(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # pybiomart installs a requests cache in the working directory at import time.
    monkeypatch.chdir(tmp_path)
    pybiomart = _import_or_skip("pybiomart")
    from pybiomart.dataset import Attribute

    attributes = {
        "ensembl_gene_id": Attribute("ensembl_gene_id", "Gene stable ID", "Gene stable ID"),
        "external_gene_name": Attribute("external_gene_name", "Gene name", "Gene name"),
    }
    monkeypatch.setattr(
        pybiomart.Dataset, "attributes", property(lambda self: attributes), raising=True
    )
    monkeypatch.setattr(pybiomart.Dataset, "filters", property(lambda self: {}), raising=True)
    table = "Gene stable ID\tGene name\nENSG00000000003\tTSPAN6\nENSG00000000005\t\n"
    monkeypatch.setattr(
        pybiomart.base.ServerBase,
        "get",
        lambda self, **params: SimpleNamespace(text=table, content=b""),
    )
    rows = gb.query_biomart_dataset(
        "hsapiens_gene_ensembl",
        ["ensembl_gene_id", "external_gene_name"],
        max_rows=10,
        verbose=False,
    )
    assert rows["columns"] == ["Gene stable ID", "Gene name"]
    assert rows["rows"][0] == {"Gene stable ID": "ENSG00000000003", "Gene name": "TSPAN6"}
    assert rows["rows"][1]["Gene name"] is None  # empty cells are missing, not empty strings
    attr_list = gb.list_biomart_attributes("hsapiens_gene_ensembl", verbose=False)
    assert attr_list["count"] == 2
    assert {row["name"] for row in attr_list["attributes"]} == set(attributes)
    with pytest.raises(ConfigurationError, match="dataset"):
        gb.query_biomart_dataset("bad name", ["ensembl_gene_id"], verbose=False)
    with pytest.raises(ConfigurationError, match="max_rows=1"):
        gb.query_biomart_dataset(
            "hsapiens_gene_ensembl", ["ensembl_gene_id"], max_rows=1, verbose=False
        )


def test_biomart_list_datasets_uses_mart_listing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.chdir(tmp_path)
    pybiomart = _import_or_skip("pybiomart")
    import pandas as pd

    class FakeMart:
        def list_datasets(self) -> pd.DataFrame:
            return pd.DataFrame(
                [("hsapiens_gene_ensembl", "Human genes (GRCh38.p14)")],
                columns=["name", "display_name"],
            )

    monkeypatch.setattr(
        pybiomart.Server, "marts", property(lambda self: {"ENSEMBL_MART_ENSEMBL": FakeMart()})
    )
    result = gb.list_biomart_datasets(verbose=False)
    assert result["datasets"] == [
        {"name": "hsapiens_gene_ensembl", "display_name": "Human genes (GRCh38.p14)"}
    ]
    with pytest.raises(ConfigurationError, match="not listed"):
        gb.list_biomart_datasets("NOT_A_MART", verbose=False)


# --- ncbi-datasets-pylib: real signatures when installed, fake module otherwise ------------


def test_ncbi_gene_api_signatures_match_installed_package() -> None:
    # Skipped when ncbi-datasets-pylib is not importable; it has no published wheel for this venv.
    ncbi_datasets = pytest.importorskip("ncbi.datasets")
    api = ncbi_datasets.GeneApi
    for method in (
        "gene_metadata_by_id",
        "gene_metadata_by_tax_and_symbol",
        "gene_orthologs_by_id",
    ):
        assert callable(getattr(api, method))
    assert "limit" in (api.gene_metadata_by_id.__doc__ or "")


class _FakeNcbiResponse:
    def __init__(self, payload: dict[str, Any]) -> None:
        self._payload = payload

    def to_dict(self) -> dict[str, Any]:
        return self._payload


def _fake_ncbi_module(calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]]) -> Any:
    """A stand-in for ncbi.datasets whose GeneApi records calls and returns openapi-like models."""

    match = {
        "query": ["1017"],
        "gene": {"gene_id": "1017", "symbol": "CDK2", "tax_id": "9606"},
        "messages": [],
        "warnings": [],
        "errors": [],
    }

    class FakeGeneApi:
        def gene_metadata_by_id(self, gene_ids: Any, **kwargs: Any) -> _FakeNcbiResponse:
            calls.append(("gene_metadata_by_id", (gene_ids,), kwargs))
            return _FakeNcbiResponse({"genes": [match], "total_count": 1})

        def gene_metadata_by_tax_and_symbol(
            self, symbols: Any, taxon: Any, **kwargs: Any
        ) -> _FakeNcbiResponse:
            calls.append(("gene_metadata_by_tax_and_symbol", (symbols, taxon), kwargs))
            return _FakeNcbiResponse({"genes": [match], "total_count": 1})

        def gene_orthologs_by_id(self, gene_id: Any, **kwargs: Any) -> _FakeNcbiResponse:
            calls.append(("gene_orthologs_by_id", (gene_id,), kwargs))
            return _FakeNcbiResponse({"ortholog_set_id": 7, "method": "homology", "genes": [match]})

    module = types.ModuleType("ncbi.datasets")
    module.GeneApi = FakeGeneApi  # type: ignore[attr-defined]
    return module


def test_fetch_ncbi_gene_by_id_bounds_server_limit_and_returns_plain_genes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []
    monkeypatch.setitem(sys.modules, "ncbi.datasets", _fake_ncbi_module(calls))
    result = gb.fetch_ncbi_gene_by_id([1017], max_genes=5, verbose=False)
    assert calls[0] == ("gene_metadata_by_id", ([1017],), {"limit": "5"})
    assert result["total_count"] == 1
    assert result["genes"][0]["gene"]["symbol"] == "CDK2"
    with pytest.raises(ConfigurationError, match="gene_ids"):
        gb.fetch_ncbi_gene_by_id([True], verbose=False)
    with pytest.raises(ConfigurationError, match="gene_ids"):
        gb.fetch_ncbi_gene_by_id([0], verbose=False)
    with pytest.raises(ConfigurationError, match="max_genes"):
        gb.fetch_ncbi_gene_by_id([1017], max_genes=0, verbose=False)


def test_fetch_ncbi_gene_by_symbol_and_orthologs_pass_arguments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []
    monkeypatch.setitem(sys.modules, "ncbi.datasets", _fake_ncbi_module(calls))
    by_symbol = gb.fetch_ncbi_gene_by_symbol(["CDK2"], "human", max_genes=3, verbose=False)
    assert calls[0] == (
        "gene_metadata_by_tax_and_symbol",
        (["CDK2"], "human"),
        {"limit": "3"},
    )
    assert by_symbol["genes"][0]["gene"]["symbol"] == "CDK2"
    orthologs = gb.fetch_ncbi_gene_orthologs(1017, taxon_filter=["9606"], verbose=False)
    assert calls[1] == ("gene_orthologs_by_id", (1017,), {"taxon_filter": ["9606"]})
    assert orthologs["ortholog_set_id"] == 7 and orthologs["method"] == "homology"
    gb.fetch_ncbi_gene_orthologs(1017, verbose=False)
    assert calls[2] == ("gene_orthologs_by_id", (1017,), {})
    with pytest.raises(ConfigurationError, match="taxon"):
        gb.fetch_ncbi_gene_by_symbol(["CDK2"], "human;drop", verbose=False)


def test_ncbi_upstream_errors_are_backend_execution_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _fake_ncbi_module([])

    def failing(self: Any, gene_ids: Any, **kwargs: Any) -> None:
        raise OSError("connection reset")

    monkeypatch.setattr(module.GeneApi, "gene_metadata_by_id", failing)
    monkeypatch.setitem(sys.modules, "ncbi.datasets", module)
    with pytest.raises(BackendExecutionError, match="gene_metadata_by_id failed"):
        gb.fetch_ncbi_gene_by_id([1017], verbose=False)


# --- missing optional packages and validation ---------------------------------------------


@pytest.mark.parametrize(
    ("module", "extra", "call"),
    [
        ("gget", "bio-gget", lambda: gb.search_ensembl_genes("CDK2")),
        ("pyensembl", "bio-pyensembl", lambda: gb.index_gene_models("/nonexistent.gtf")),
        ("genomepy", "bio-genomepy", lambda: gb.read_genome_summary("/nonexistent.fa")),
        ("biothings_client", "bio-biothings", lambda: gb.format_hgvs_id("1", 1, "A", "G")),
        ("mygene", "bio-mygene", lambda: gb.query_mygene_genes("symbol:CDK2")),
        ("myvariant", "bio-myvariant", lambda: gb.fetch_variant_annotations(["chr1:g.1A>G"])),
        ("pybiomart", "bio-pybiomart", lambda: gb.list_biomart_datasets()),
        ("ncbi.datasets", "bio-ncbi-datasets", lambda: gb.fetch_ncbi_gene_by_id([1017])),
    ],
)
def test_missing_optional_package_raises_install_hint(
    monkeypatch: pytest.MonkeyPatch, module: str, extra: str, call: Any
) -> None:
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.raises(BackendUnavailableError, match=extra):
        call()


def test_local_functions_validate_inputs_before_running_upstream(
    local_files: dict[str, Path],
) -> None:
    pytest.importorskip("pyensembl")
    cache = local_files["cache"]
    with pytest.raises(ConfigurationError, match="max_results"):
        gb.query_genes_by_name(
            local_files["gtf"], ["GENEA"], cache_dir=cache, max_results=True, verbose=False
        )
    with pytest.raises(ConfigurationError, match="Input file does not exist"):
        gb.index_gene_models(local_files["gtf"].parent / "missing.gtf", verbose=False)
    with pytest.raises(ConfigurationError, match="gene_id"):
        gb.fetch_gene_sequence(
            local_files["gtf"], local_files["fasta"], "bad id!", cache_dir=cache, verbose=False
        )


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="live Ensembl and MyGene services; set DNAKIT_LIVE_TESTS=1 to run",
)
def test_live_gene_services_return_cdk2() -> None:
    info = gb.fetch_ensembl_gene_info(["ENSG00000123374"], verbose=False)
    assert info["records"]["ENSG00000123374"]["ensembl_gene_name"] == "CDK2"
    hits = gb.query_mygene_genes("symbol:CDK2", size=1, verbose=False)
    assert hits["hits"][0]["symbol"] == "CDK2"


def test_every_name_in_all_is_a_public_callable() -> None:
    assert len(gb.__all__) == len(set(gb.__all__)) == 25
    for name in gb.__all__:
        assert callable(getattr(gb, name)) and not name.startswith("_")
