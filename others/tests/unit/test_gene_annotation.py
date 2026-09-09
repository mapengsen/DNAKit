"""Coordinate, resource-selection, cache and command-line gene annotation tests."""

from __future__ import annotations

import hashlib
import json
import os
import random
import shutil
import zipfile
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from typer.testing import CliRunner

import dnakit.annotation._gene_reference as reference_module
import dnakit.cli.gene_annotation as cli_module
from dnakit.annotation import GeneAnnotationConfig, GeneAnnotationResult, annotate_genes
from dnakit.annotation._gene_mapping import ReferenceGene, annotate_mappings, read_mappings
from dnakit.annotation._gene_reference import fasta_lengths, prepare_reference, read_genes
from dnakit.cli.app import app
from dnakit.core import Interval
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError, DownloadError, InputFormatError

ACCESSION = "GCF_000000001.1"


def _config(tmp_path: Path, **kwargs: Any) -> GeneAnnotationConfig:
    return GeneAnnotationConfig(cache_dir=str(tmp_path / "cache"), min_alignment_length=1, **kwargs)


def _paf(
    qlen: int,
    qs: int,
    qe: int,
    ts: int,
    te: int,
    cigar: str,
    *,
    target: str = "ref",
    target_length: int = 2000,
    strand: str = "+",
    matches: int | None = None,
    span: int | None = None,
    mapq: int = 60,
    tp: str = "P",
) -> str:
    return (
        "\t".join(
            map(
                str,
                (
                    "query",
                    qlen,
                    qs,
                    qe,
                    strand,
                    target,
                    target_length,
                    ts,
                    te,
                    qe - qs if matches is None else matches,
                    qe - qs if span is None else span,
                    mapq,
                    f"cg:Z:{cigar}",
                    f"tp:A:{tp}",
                ),
            )
        )
        + "\n"
    )


def _mappings(tmp_path: Path, text: str, length: int, **kwargs: Any) -> Any:
    path = tmp_path / "alignment.paf"
    path.write_text(text)
    return read_mappings(path, length, {"ref": 2000}, _config(tmp_path, **kwargs))


def _gene(name: str, start: int, end: int, strand: str = "+") -> ReferenceGene:
    return ReferenceGene(
        "ref", name, name, name, None, "protein_coding", strand, (Interval(start, end),)
    )


def test_requested_52000_base_example_preserves_overlapping_genes(tmp_path: Path) -> None:
    paf = tmp_path / "input.paf"
    paf.write_text(_paf(52000, 0, 52000, 0, 52000, "52000M", target_length=52000))
    mappings = read_mappings(paf, 52000, {"ref": 52000}, _config(tmp_path))
    records, segments = annotate_mappings(
        52000,
        mappings,
        (_gene("A", 0, 12000), _gene("B", 18000, 45000), _gene("C", 42999, 52000)),
        {"ref"},
    )
    assert [(r.gene_name, r.location) for r in records] == [
        ("A", Interval(0, 12000)),
        ("B", Interval(18000, 45000)),
        ("C", Interval(42999, 52000)),
    ]
    assert [(s.to_dict()["start"], s.to_dict()["end"], s.kind, s.gene_names) for s in segments] == [
        (1, 12000, "gene", ("A",)),
        (12001, 18000, "intergenic", ()),
        (18001, 42999, "gene", ("B",)),
        (43000, 45000, "gene", ("B", "C")),
        (45001, 52000, "gene", ("C",)),
    ]


@pytest.mark.parametrize(
    "strand,expected",
    [
        ("+", (Interval(2, 4), Interval(6, 10))),
        ("-", (Interval(2, 6), Interval(8, 10))),
    ],
)
def test_cigar_indels_and_reverse_coordinates(
    tmp_path: Path, strand: str, expected: tuple[Interval, ...]
) -> None:
    mappings = _mappings(
        tmp_path,
        _paf(12, 0, 12, 100, 112, "4M2I3M2D3M", strand=strand, matches=10, span=14),
        12,
        min_identity=0.5,
    )
    records, segments = annotate_mappings(
        12, mappings, (_gene("A", 102, 110, "-"), _gene("deleted", 107, 109)), {"ref"}
    )
    assert len(records) == 1
    assert records[0].parts == expected
    assert records[0].strand == ("-" if strand == "+" else "+")
    assert sum(len(s.location) for s in segments if s.kind == "unknown") == 2
    assert sum(len(s.location) for s in segments) == 12


def test_split_and_ambiguous_primary_mappings(tmp_path: Path) -> None:
    mappings = _mappings(
        tmp_path,
        _paf(100, 0, 70, 100, 170, "70M") + _paf(100, 50, 90, 200, 240, "40M"),
        100,
    )
    records, segments = annotate_mappings(100, mappings, (_gene("A", 100, 240),), {"ref"})
    assert [(s.location, s.kind, s.reason) for s in segments] == [
        (Interval(0, 50), "gene", None),
        (Interval(50, 70), "unknown", "ambiguous_mapping"),
        (Interval(70, 90), "gene", None),
        (Interval(90, 100), "unknown", "unmapped"),
    ]
    assert [r.parts for r in records] == [(Interval(0, 50),), (Interval(70, 90),)]


@pytest.mark.parametrize("mapq,tp", [(0, "P"), (255, "P"), (60, "S")])
def test_unreliable_or_secondary_mapping_never_implies_intergenic(
    tmp_path: Path, mapq: int, tp: str
) -> None:
    mappings = _mappings(tmp_path, _paf(20, 0, 20, 0, 20, "20M", mapq=mapq, tp=tp), 20)
    records, segments = annotate_mappings(20, mappings, (_gene("A", 100, 200),), {"ref"})
    assert not records
    assert len(segments) == 1 and segments[0].kind == "unknown"


def test_empty_mapping_and_unannotated_contig_are_unknown(tmp_path: Path) -> None:
    assert annotate_mappings(10, (), (), set())[1][0].reason == "unmapped"
    mappings = _mappings(tmp_path, _paf(10, 0, 10, 0, 10, "10M"), 10)
    segment = annotate_mappings(10, mappings, (), set())[1][0]
    assert segment.kind == "unknown" and segment.reason == "missing_contig_annotation"


@pytest.mark.parametrize("replacement", ["9M", "10S", "4M0I6M", "4M6Mgarbage", "10M10M"])
def test_malformed_cigar_is_rejected(tmp_path: Path, replacement: str) -> None:
    with pytest.raises(InputFormatError, match="INVALID_GENE_MAPPING"):
        _mappings(tmp_path, _paf(10, 0, 10, 0, 10, replacement), 10)


def _local_pair(tmp_path: Path, *, length: int = 2000) -> tuple[Path, Path]:
    fasta, gff = tmp_path / "reference.fna", tmp_path / "genes.gff3"
    fasta.write_text(">ref\n" + "A" * length + "\n")
    gff.write_text(
        "##gff-version 3\n"
        f"#!genome-build-accession NCBI_Assembly:{ACCESSION}\n"
        f"##sequence-region ref 1 {length}\n"
        "ref\ttest\tgene\t101\t200\t.\t+\t.\tID=g1;Name=alpha;Dbxref=GeneID:123;gene=alpha\n"
        "ref\ttest\tgene\t151\t300\t.\t-\t.\tID=g2;Name=LOC2;locus_tag=LOC2\n"
    )
    return fasta, gff


def test_gff_gene_name_id_fallback_and_overlap(tmp_path: Path) -> None:
    fasta, gff = _local_pair(tmp_path)
    lengths = fasta_lengths(fasta, _config(tmp_path))
    genes, represented, headers = read_genes(gff, lengths, _config(tmp_path))
    assert genes[0].gene_name == "alpha" and genes[0].gene_id == "GeneID:123"
    assert genes[1].gene_name is None and genes[1].locus_tag == "LOC2"
    assert represented == {"ref"}
    assert headers["genome-build-accession"].endswith(ACCESSION)


def test_gff_declared_circular_gene_can_cross_origin(tmp_path: Path) -> None:
    _, gff = _local_pair(tmp_path)
    gff.write_text(
        "##gff-version 3\n"
        "ref\ttest\tregion\t1\t2000\t.\t+\t.\tID=region;Is_circular=true\n"
        "ref\ttest\tgene\t1901\t2100\t.\t+\t.\tID=circular;gene=wrap\n"
    )
    genes, _, _ = read_genes(gff, {"ref": 2000}, _config(tmp_path))
    assert genes[0].parts == (Interval(1900, 2000), Interval(0, 100))


def test_gff_without_genes_and_mismatching_contigs_fail(tmp_path: Path) -> None:
    _, gff = _local_pair(tmp_path)
    with pytest.raises(InputFormatError):
        read_genes(gff, {"other_ref": 2000}, _config(tmp_path))
    gff.write_text("##gff-version 3\nref\ttest\tCDS\t1\t9\t.\t+\t0\tID=c1\n")
    with pytest.raises(InputFormatError, match="GENE_ANNOTATION_UNAVAILABLE"):
        read_genes(gff, {"ref": 2000}, _config(tmp_path))


def _report() -> dict[str, Any]:
    return {
        "accession": ACCESSION,
        "organism": {"tax_id": 9606, "organism_name": "Homo sapiens"},
        "assembly_info": {"assembly_name": "test-assembly"},
        "annotation_info": {"name": "test annotation", "release_date": "2026-01-01"},
    }


def _response(records: list[dict[str, Any]]) -> Any:
    return SimpleNamespace(records=records, total_count=len(records), next_page_token=None)


def test_optional_assembly_selects_reference_but_explicit_version_is_preserved(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    selected = []

    def query_assembly(query: str, **kwargs: Any) -> Any:
        selected.append((query, kwargs))
        return _response([_report()])

    monkeypatch.setattr(reference_module, "assembly", query_assembly)
    config = reference_module.SearchConfig()
    assert reference_module._resolve_report(9606, None, config)[0] == ACCESSION
    assert selected[0][0] == "9606"
    assert selected[0][1]["reference_only"] is True
    assert selected[0][1]["assembly_source"] == "refseq"
    selected.clear()
    assert reference_module._resolve_report(9606, ACCESSION, config)[0] == ACCESSION
    assert len(selected) == 1 and selected[0][0] == ACCESSION
    assert selected[0][1]["current_only"] is False


def test_species_and_annotation_availability_are_validated(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(reference_module, "assembly", lambda *a, **kw: _response([_report()]))
    monkeypatch.setattr(reference_module, "_taxon", lambda *a: {"tax_id": 9606, "parents": [9605]})
    with pytest.raises(ConfigurationError, match="ASSEMBLY_SPECIES_MISMATCH"):
        reference_module._resolve_report(10090, ACCESSION, reference_module.SearchConfig())
    report = _report()
    report.pop("annotation_info")
    monkeypatch.setattr(reference_module, "assembly", lambda *a, **kw: _response([report]))
    with pytest.raises(ConfigurationError, match="GENE_ANNOTATION_UNAVAILABLE"):
        reference_module._resolve_report(9606, ACCESSION, reference_module.SearchConfig())


def test_strain_reference_is_accepted_for_parent_species(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(reference_module, "_taxon", lambda *a: {"tax_id": 511145, "parents": [562]})
    reference_module._validate_species(562, {"tax_id": 511145}, reference_module.SearchConfig())
    with pytest.raises(ConfigurationError, match="ASSEMBLY_SPECIES_MISMATCH"):
        reference_module._validate_species(999, {"tax_id": 511145}, reference_module.SearchConfig())


def test_unavailable_or_ambiguous_defaults_are_not_replaced(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        reference_module,
        "assembly",
        lambda *a, **kw: _response([_report(), {**_report(), "accession": "GCF_000000002.1"}]),
    )
    with pytest.raises(ConfigurationError, match="AMBIGUOUS_GENOME_ASSEMBLY"):
        reference_module._resolve_report(9606, None, reference_module.SearchConfig())
    monkeypatch.setattr(reference_module, "assembly", lambda *a, **kw: _response([]))
    with pytest.raises(ConfigurationError, match="GENOME_ASSEMBLY_NOT_FOUND"):
        reference_module._resolve_report(9606, None, reference_module.SearchConfig())


def test_assembly_name_search_requires_a_unique_candidate(monkeypatch: pytest.MonkeyPatch) -> None:
    first = _report()
    second = {**first, "accession": "GCA_000000001.1"}
    monkeypatch.setattr(reference_module, "assembly", lambda *a, **kw: _response([first, second]))
    with pytest.raises(ConfigurationError, match="AMBIGUOUS_GENOME_ASSEMBLY"):
        reference_module._resolve_report(9606, "test-assembly", reference_module.SearchConfig())


def _archive(tmp_path: Path, *, bad_checksum: bool = False, camel_case: bool = True) -> Path:
    fasta, gff = _local_pair(tmp_path)
    files = {
        f"ncbi_dataset/data/{ACCESSION}/test_genomic.fna": fasta.read_bytes(),
        f"ncbi_dataset/data/{ACCESSION}/genomic.gff": gff.read_bytes(),
    }
    checksums = "".join(
        f"{('0' * 32) if bad_checksum else hashlib.md5(value).hexdigest()}  {name}\n"
        for name, value in files.items()
    )
    path = tmp_path / "genome.zip"
    report = _report()
    if camel_case:
        report = {
            "accession": ACCESSION,
            "organism": {"taxId": 9606, "organismName": "Homo sapiens"},
            "assemblyInfo": {"assemblyName": "test-assembly"},
            "annotationInfo": {"name": "test annotation", "releaseDate": "2026-01-01"},
        }
    with zipfile.ZipFile(path, "w") as archive:
        for name, data in files.items():
            archive.writestr(name, data)
        archive.writestr("md5sum.txt", checksums)
        archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl", json.dumps(report) + "\n")
    return path


def _mock_remote(monkeypatch: pytest.MonkeyPatch, archive: Path) -> list[str]:
    calls = []
    monkeypatch.setattr(reference_module, "_taxon", lambda *a: {"tax_id": 9606})
    monkeypatch.setattr(reference_module, "_resolve_report", lambda *a: (ACCESSION, _report()))

    def download(accession: str, *args: Any, **kwargs: Any) -> Any:
        calls.append(accession)
        return SimpleNamespace(files=(SimpleNamespace(path=str(archive), sha256="a" * 64),))

    monkeypatch.setattr(reference_module, "genome_package", download)
    return calls


@pytest.mark.parametrize("camel_case", [True, False])
def test_default_reference_and_annotation_are_pinned_and_reused_offline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, camel_case: bool
) -> None:
    archive = _archive(tmp_path, camel_case=camel_case)
    calls = _mock_remote(monkeypatch, archive)
    first = prepare_reference("human", None, _config(tmp_path))
    assert first[2]["selection"] == "automatic"
    assert first[2]["organism"]["tax_id"] == 9606
    assert first[2]["assembly_name"] == "test-assembly"
    assert first[2]["annotation"]["release_date"] == "2026-01-01"

    def no_network(*args: Any, **kwargs: Any) -> Any:
        pytest.fail("A pinned request must not query NCBI again")

    monkeypatch.setattr(reference_module, "_taxon", no_network)
    monkeypatch.setattr(reference_module, "_resolve_report", no_network)
    second = prepare_reference("human", None, _config(tmp_path, allow_download=False))
    assert second[:2] == first[:2]
    assert second[2]["selection"] == "cached"
    assert calls == [ACCESSION]
    second[1].write_text("corrupted\n")
    with pytest.raises(ConfigurationError, match="REFERENCE_CACHE_INVALID"):
        prepare_reference("human", None, _config(tmp_path, allow_download=False))


def test_bad_download_never_installs_reference_or_pin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = _archive(tmp_path, bad_checksum=True)
    _mock_remote(monkeypatch, archive)
    with pytest.raises(DownloadError, match="MD5"):
        prepare_reference("human", None, _config(tmp_path))
    assert not (tmp_path / "cache" / "assemblies" / ACCESSION).exists()
    assert not (tmp_path / "cache" / "requests").exists()


def test_no_download_without_cache_fails_before_any_network(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="REFERENCE_NOT_CACHED"):
        prepare_reference("human", None, _config(tmp_path, allow_download=False))


@pytest.mark.parametrize(
    "kwargs",
    [{"min_mapq": 0}, {"min_identity": float("nan")}, {"threads": True}, {"preset": "splice"}],
)
def test_invalid_mapping_configuration_is_rejected(kwargs: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        GeneAnnotationConfig(**kwargs)


def test_input_validation_precedes_reference_or_process_work(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="IUPAC"):
        annotate_genes("AT-GC", species="human")
    with pytest.raises(ConfigurationError, match="both"):
        annotate_genes("ATGC", species="human", reference_fasta=tmp_path / "not-present.fa")


def test_cli_optional_assembly_and_output_no_clobber(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">user-id\nACGTACGT\n")
    output = tmp_path / "result.json"
    calls = []

    def run(sequence: str, **kwargs: Any) -> GeneAnnotationResult:
        calls.append(kwargs)
        return GeneAnnotationResult(
            kwargs["sequence_id"],
            len(sequence),
            FrozenDict({"assembly": ACCESSION}),
            (),
            annotate_mappings(8, (), (), set())[1],
            (),
            FrozenDict(),
        )

    monkeypatch.setattr(cli_module, "annotate_genes", run)
    args = ["annotate-genes", "--input", str(query), "--species", "human", "--output", str(output)]
    runner = CliRunner()
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert calls[0]["assembly"] is None
    assert json.loads(output.read_text())["segments"][0] == {
        "start": 1,
        "end": 8,
        "type": "unknown",
        "gene_record_ids": [],
        "gene_ids": [],
        "gene_names": [],
        "mapping_id": None,
        "reason": "unmapped",
    }
    assert runner.invoke(app, args).exit_code == 2
    assert len(calls) == 1
    args += ["--overwrite", "--assembly", ACCESSION]
    assert runner.invoke(app, args).exit_code == 0
    assert calls[-1]["assembly"] == ACCESSION


def test_cli_rejects_multiple_records(tmp_path: Path) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">a\nACGT\n>b\nACGT\n")
    result = CliRunner().invoke(
        app, ["annotate-genes", "--input", str(query), "--species", "human"]
    )
    assert result.exit_code == 2 and "exactly one" in result.output


def test_real_minimap2_forward_reverse_and_index_reuse(tmp_path: Path) -> None:
    executable = os.environ.get("DNAKIT_TEST_MINIMAP2") or shutil.which("minimap2")
    if executable is None:
        pytest.skip("Set DNAKIT_TEST_MINIMAP2 or install minimap2 for the real-process test")
    rng = random.Random(483)
    reference = "".join(rng.choices("ACGT", k=6000))
    fasta, gff = tmp_path / "ref.fa", tmp_path / "genes.gff3"
    fasta.write_text(">ref\n" + reference + "\n")
    gff.write_text(
        "##gff-version 3\n##sequence-region ref 1 6000\n"
        "ref\ttest\tgene\t1001\t2000\t.\t+\t.\tID=A;gene=alpha\n"
        "ref\ttest\tgene\t3001\t4300\t.\t-\t.\tID=B;gene=beta\n"
        "ref\ttest\tgene\t4201\t4900\t.\t+\t.\tID=C;gene=gamma\n"
    )
    config = _config(tmp_path, minimap2_executable=executable, allow_download=False)
    sequence = reference[700:5000]
    forward = annotate_genes(
        sequence,
        species="synthetic test organism",
        reference_fasta=fasta,
        annotation_gff=gff,
        config=config,
    )
    assert [(r.gene_name, r.location) for r in forward.gene_records] == [
        ("alpha", Interval(300, 1300)),
        ("beta", Interval(2300, 3600)),
        ("gamma", Interval(3500, 4200)),
    ]
    assert any(segment.gene_names == ("beta", "gamma") for segment in forward.segments)
    reverse = annotate_genes(
        sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1],
        species="synthetic test organism",
        reference_fasta=fasta,
        annotation_gff=gff,
        config=config,
    )
    alpha = next(record for record in reverse.gene_records if record.gene_name == "alpha")
    assert alpha.location == Interval(3000, 4000) and alpha.strand == "-"
    assert len(list((tmp_path / "cache" / "indexes").glob("*/reference.mmi"))) == 1
    assert forward.to_dict()["coordinate_system"] == "1-based-closed"
