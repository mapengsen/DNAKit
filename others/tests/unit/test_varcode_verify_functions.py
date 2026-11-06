"""Tests for dnakit.molecular_lab.varcode_verify_functions (varcode variant loading and CLI command equivalents)."""

from __future__ import annotations

import csv
import inspect
import json
import logging
import warnings
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("varcode")
pytest.importorskip("pyensembl")

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import varcode_verify_functions as vv

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENOME = DATA / "genome.fa"
VCF = DATA / "variants.vcf"
STOPS = {"TAA", "TAG", "TGA"}


def _read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def _orf(sequence: str) -> tuple[int, int]:
    """First ATG (0-based) at or after index 1000 whose frame reaches a stop; returns [start, end) with the stop."""
    for index in range(1000, len(sequence) - 3):
        if sequence[index : index + 3] != "ATG":
            continue
        end = index
        while end + 3 <= len(sequence) and sequence[end : end + 3] not in STOPS:
            end += 3
        if end + 3 <= len(sequence) and 60 <= end - index <= 600:
            return index, end + 3
    raise AssertionError("no ORF in the shared genome")


def _gtf_row(feature: str, start: int, end: int, attributes: str) -> str:
    score_or_frame = "0" if feature == "CDS" else "."
    return "\t".join(
        ["chr1", "test", feature, str(start), str(end), ".", "+", score_or_frame, attributes]
    )


@pytest.fixture(scope="module")
def reference(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    if not GENOME.exists() or not VCF.exists():
        pytest.skip("shared consistency data not present")
    chr1 = _read_fasta(GENOME)["chr1"].upper()
    start, end = _orf(chr1)
    span = (start + 1, end)
    tx = 'gene_id "G1"; transcript_id "T1"; transcript_biotype "protein_coding";'
    lines = [
        _gtf_row("gene", *span, 'gene_id "G1"; gene_name "G1"; gene_biotype "protein_coding";'),
        _gtf_row("transcript", *span, tx),
        _gtf_row("exon", *span, tx + ' exon_number "1";'),
        _gtf_row("CDS", *span, tx + ' exon_number "1"; protein_id "P1";'),
        _gtf_row("start_codon", start + 1, start + 3, tx + ' exon_number "1";'),
        _gtf_row("stop_codon", end - 2, end, tx + ' exon_number "1";'),
    ]
    directory = tmp_path_factory.mktemp("varcode_verify_reference")
    gtf = directory / "genes.gtf"
    gtf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "gtf": str(gtf),
        "fasta": str(GENOME),
        "cache": str(directory / "cache"),
        "chr1": chr1,
        "span": span,
    }


def _local_genome(reference: dict[str, Any]) -> Any:
    from dnakit.molecular_lab.varcode_gap_functions import _genome

    genome, _ = _genome(reference["gtf"], reference["fasta"], reference["cache"])
    return genome


def _quiet(function: Any, *args: Any, **kwargs: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        logging.disable(logging.WARNING)
        try:
            return function(*args, **kwargs)
        finally:
            logging.disable(logging.NOTSET)


def test_public_functions_take_keyword_only_verbose() -> None:
    for name in vv.__all__:
        parameter = inspect.signature(getattr(vv, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_no_variant_source_is_configuration_error() -> None:
    with pytest.raises(ConfigurationError, match="No variants were given"):
        vv.varcode_load_variant_sources(verbose=False)


def test_missing_vcf_is_configuration_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="does not exist"):
        vv.varcode_load_variant_sources([str(tmp_path / "absent.vcf")], verbose=False)


@pytest.mark.skipif(not VCF.exists(), reason="shared consistency variants.vcf not present")
def test_vcf_without_reference_needs_gtf_or_genome() -> None:
    with pytest.raises(ConfigurationError, match="could not load"):
        _quiet(vv.varcode_load_variant_sources, [str(VCF)], verbose=False)


def test_vcf_rows_match_upstream_load_with_the_same_reference(reference: dict[str, Any]) -> None:
    from varcode import load_vcf

    result = _quiet(
        vv.varcode_load_variant_sources,
        [str(VCF)],
        gtf=reference["gtf"],
        genome_fasta=reference["fasta"],
        cache_dir=reference["cache"],
        verbose=False,
    )
    expected = _quiet(load_vcf, str(VCF), genome=_local_genome(reference), only_passing=True)
    rows = sorted((v["contig"], v["start"], v["ref"], v["alt"]) for v in result["variants"])
    want = sorted((str(v.contig), int(v.start), str(v.ref), str(v.alt)) for v in expected)
    assert rows == want
    assert result["count"] == len(want) > 0
    assert result["sources"] == {"vcf": 1, "maf": 0, "json": 0, "variant": 0}
    assert result["reference"].startswith("local:")


def test_explicit_records_and_json_are_unioned(reference: dict[str, Any], tmp_path: Path) -> None:
    from varcode import Variant, VariantCollection

    start, _ = reference["span"]
    ref = reference["chr1"][
        start
    ]  # 1-based position start+1 is the first base after the 0-based index
    alt = "C" if ref != "C" else "G"
    genome = _local_genome(reference)
    collection = VariantCollection(
        [
            Variant(
                "chr1", start=start + 10, ref=reference["chr1"][start + 9], alt="A", genome=genome
            )
        ]
    )
    js = tmp_path / "variants.json"
    js.write_text(collection.to_json(), encoding="utf-8")
    records = [
        {"contig": "chr1", "start": start + 10, "ref": reference["chr1"][start + 9], "alt": "A"},
        {"contig": "chr1", "start": start + 20, "ref": reference["chr1"][start + 19], "alt": alt},
    ]
    result = _quiet(
        vv.varcode_load_variant_sources,
        json_paths=[str(js)],
        variants=records,
        gtf=reference["gtf"],
        genome_fasta=reference["fasta"],
        cache_dir=reference["cache"],
        verbose=False,
    )
    assert result["sources"] == {"vcf": 0, "maf": 0, "json": 1, "variant": 2}
    # Same union as varcode's own VariantCollection.union over the same inputs (the CLI behaviour).
    explicit = VariantCollection(
        [
            Variant(
                "chr1",
                start=start + 10,
                ref=records[0]["ref"],
                alt="A",
                genome=genome,
                convert_ucsc_contig_names=False,
            ),
            Variant(
                "chr1",
                start=start + 20,
                ref=records[1]["ref"],
                alt=alt,
                genome=genome,
                convert_ucsc_contig_names=False,
            ),
        ]
    )
    expected = _quiet(
        VariantCollection.union,
        VariantCollection.from_json(js.read_text(encoding="utf-8")),
        explicit,
    )
    assert result["count"] == len(expected) >= 2
    assert {v["start"] for v in result["variants"]} >= {start + 10, start + 20}


def test_explicit_record_without_reference_is_configuration_error() -> None:
    with pytest.raises(ConfigurationError, match="Explicit variants need"):
        vv.varcode_load_variant_sources(
            variants=[{"contig": "chr1", "start": 5, "ref": "A", "alt": "C"}], verbose=False
        )


def test_invalid_variant_record_is_configuration_error() -> None:
    with pytest.raises(ConfigurationError, match="contig, start, ref and alt"):
        vv.varcode_load_variant_sources(variants=[{"contig": "chr1"}], verbose=False)


def test_effects_table_matches_upstream_effects(reference: dict[str, Any], tmp_path: Path) -> None:
    from varcode import Variant, VariantCollection

    start, _ = reference["span"]
    position = start + 10  # 1-based position inside the CDS
    ref = reference["chr1"][position - 1]
    alt = "C" if ref != "C" else "G"
    out = tmp_path / "effects.csv"
    result = _quiet(
        vv.varcode_effects_table,
        str(out),
        gtf=reference["gtf"],
        genome_fasta=reference["fasta"],
        cache_dir=reference["cache"],
        variants=[{"contig": "chr1", "start": position, "ref": ref, "alt": alt}],
        verbose=False,
    )
    upstream = VariantCollection(
        [Variant("chr1", start=position, ref=ref, alt=alt, genome=_local_genome(reference))]
    ).effects(raise_on_error=True)
    with out.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert result["variant_count"] == 1
    assert result["effect_count"] == len(rows) == len(upstream) >= 1
    assert sorted(result["class_counts"].items()) == sorted(
        {
            name: sum(type(e).__name__ == name for e in upstream)
            for name in {type(e).__name__ for e in upstream}
        }.items()
    )
    assert set(rows[0]) == set(vv._EFFECT_FIELDS)
    assert rows[0]["transcript_id"] == "T1"


def test_effects_table_refuses_existing_output(reference: dict[str, Any], tmp_path: Path) -> None:
    existing = tmp_path / "exists.csv"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="must not exist"):
        vv.varcode_effects_table(
            str(existing),
            gtf=reference["gtf"],
            genome_fasta=reference["fasta"],
            variants=[{"contig": "chr1", "start": 1100, "ref": "A", "alt": "C"}],
            verbose=False,
        )


@pytest.mark.skipif(not VCF.exists(), reason="shared consistency variants.vcf not present")
def test_sample_identity_files_match_upstream_cli(tmp_path: Path) -> None:
    from varcode.cli import sample_checks

    upstream_json = tmp_path / "upstream.json"
    upstream_tsv = tmp_path / "upstream.tsv"
    code = _quiet(
        sample_checks.main,
        [
            "--germline",
            str(VCF),
            "--assembly",
            "GRCh38",
            "--json",
            str(upstream_json),
            "--tsv",
            str(upstream_tsv),
        ],
    )
    assert code in (None, 0)
    ours_json = tmp_path / "ours.json"
    ours_tsv = tmp_path / "ours.tsv"
    result = _quiet(
        vv.varcode_check_sample_identity_files,
        germline_vcfs=[str(VCF)],
        assembly="GRCh38",
        json_output=str(ours_json),
        tsv_output=str(ours_tsv),
        verbose=False,
    )
    assert json.loads(ours_json.read_text(encoding="utf-8")) == json.loads(
        upstream_json.read_text(encoding="utf-8")
    )
    assert ours_tsv.read_text(encoding="utf-8") == upstream_tsv.read_text(encoding="utf-8")
    assert result["sample_count"] == 24
    assert result["pair_count"] == len(result["report"]["pairs"]) > 0


@pytest.mark.skipif(not VCF.exists(), reason="shared consistency variants.vcf not present")
def test_sample_identity_outputs_must_not_overwrite_inputs(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="must differ"):
        _quiet(
            vv.varcode_check_sample_identity_files,
            germline_vcfs=[str(VCF)],
            assembly="GRCh38",
            json_output=str(VCF),
            verbose=False,
        )


def test_sample_identity_needs_samples(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="No samples"):
        vv.varcode_check_sample_identity_files(verbose=False)


def test_sv_comparison_matches_upstream_cli(tmp_path: Path) -> None:
    from varcode.sv_comparison import main as upstream_main

    header = ["call_id", "caller", "sample", "build", "chrom", "pos", "ref", "alt"]
    calls = [
        ["c1", "callerA", "S1", "GRCh38", "chr1", "1000", "N", "<DEL>"],
        ["c2", "callerB", "S1", "GRCh38", "chr1", "1040", "N", "<DEL>"],
        ["c3", "callerA", "S1", "GRCh38", "chr2", "500", "N", "<DUP>"],
    ]
    source = tmp_path / "calls.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(calls)
    upstream_dir = tmp_path / "upstream"
    _quiet(
        upstream_main,
        ["--calls", str(source), "--output", str(upstream_dir), "--max-distance", "100"],
    )
    ours_dir = tmp_path / "ours"
    result = _quiet(
        vv.varcode_compare_sv_calls_directory,
        str(source),
        str(ours_dir),
        max_distance=100,
        verbose=False,
    )
    for name in ("groups.csv", "members.csv"):
        assert _read_rows(ours_dir / name) == _read_rows(upstream_dir / name)
    ours = json.loads((ours_dir / "comparison.json").read_text(encoding="utf-8"))
    upstream = json.loads((upstream_dir / "comparison.json").read_text(encoding="utf-8"))
    assert ours == upstream
    assert result["call_count"] == len(calls)
    assert result["output_dir"] == str(ours_dir.resolve())


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_sv_output_directory_must_be_new(tmp_path: Path) -> None:
    source = tmp_path / "calls.csv"
    source.write_text("call_id,caller,sample,build,chrom,pos,ref,alt\n", encoding="utf-8")
    existing = tmp_path / "exists"
    existing.mkdir()
    with pytest.raises(ConfigurationError, match="must not exist"):
        vv.varcode_compare_sv_calls_directory(str(source), str(existing), verbose=False)


def test_sv_negative_distance_is_configuration_error(tmp_path: Path) -> None:
    source = tmp_path / "calls.csv"
    source.write_text("call_id,caller,sample,build,chrom,pos,ref,alt\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="max_distance"):
        vv.varcode_compare_sv_calls_directory(
            str(source), str(tmp_path / "out"), max_distance=-1, verbose=False
        )
