"""Reference consistency reports from genome aligner PAF, separate from syntax validity."""

from __future__ import annotations

import hashlib
import json
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.core import DNA, DNARecord, DNASequence, DNASet
from dnakit.evaluation import (
    ReferenceConsistencyConfig,
    evaluate_reference_consistency,
    evaluate_validity,
)
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _queries() -> DNASet:
    return DNASet(
        [
            DNARecord(DNASequence("ACGTACGT"), "copy"),
            DNARecord(DNASequence("ACGTACGT"), "partial"),
            DNARecord(DNASequence("ACGTACGT"), "unmapped"),
        ]
    )


def _paf(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_batch_dna_reports_per_sequence_and_aggregate_validity(
    tmp_path: Path,
) -> None:
    sequences = DNA(["ACGTACGT", "TGCATGCA", "ACGTACG"])
    paf = _paf(
        tmp_path / "batch.paf",
        [
            "sequence_1\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P\tcg:Z:8M",
            "sequence_2\t8\t0\t8\t-\tchr2\t100\t30\t38\t7\t8\t60\ttp:A:P\tcg:Z:8M",
        ],
    )

    syntax = evaluate_validity(sequences, sequence_length=8)
    reference = evaluate_reference_consistency(
        sequences, paf, reference_name="human", reference_version="GRCh38"
    )

    assert [entry.subject_id for entry in syntax.entries] == [
        "sequence_1",
        "sequence_2",
        "sequence_3",
    ]
    assert [entry.subject_id for entry in reference.entries] == [
        entry.subject_id for entry in syntax.entries
    ]
    assert syntax.metrics["record_count"] == len(reference.entries) == 3
    assert syntax.metrics["valid_fraction"] == pytest.approx(
        sum(entry.metrics["is_valid"] for entry in syntax.entries) / 3
    )
    assert syntax.metrics["valid_fraction"] == pytest.approx(2 / 3)
    assert set(reference.metrics) == {"consistency_rate"}
    assert reference.metrics["consistency_rate"] == pytest.approx((1 + 7 / 8 + 0) / 3)
    assert [entry.metrics["consistency_rate"] for entry in reference.entries] == pytest.approx(
        [1, 7 / 8, 0]
    )


def test_reference_consistency_distinguishes_exact_partial_and_unmapped(tmp_path: Path) -> None:
    queries = _queries()
    path = _paf(
        tmp_path / "alignments.paf",
        [
            "copy\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P\tcg:Z:8M",
            "partial\t8\t0\t4\t+\tchr1\t100\t30\t34\t4\t4\t60\ttp:A:P\tcg:Z:4M",
            "partial\t8\t0\t8\t-\tchr2\t100\t10\t18\t7\t8\t30\ttp:A:P\tcg:Z:8M",
            "unmapped\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:S\tcg:Z:8M",
        ],
    )
    report = evaluate_reference_consistency(
        queries,
        path,
        reference_name="human",
        reference_version="GRCh38",
        reference_sha256="a" * 64,
    )

    assert evaluate_validity(queries).metrics["valid_fraction"] == 1.0
    assert report.metrics["consistency_rate"] == pytest.approx((1 + 7 / 8 + 0) / 3)
    assert [entry.metrics["consistency_rate"] for entry in report.entries] == pytest.approx(
        [1, 7 / 8, 0]
    )
    assert all(set(entry.metrics) == {"consistency_rate"} for entry in report.entries)
    assert report.provenance.reference is not None
    assert report.provenance.reference.checksum == "a" * 64
    assert json.loads(json.dumps(report.to_dict()))["metrics"]["consistency_rate"] == pytest.approx(
        0.625
    )


def test_reference_consistency_counts_partial_and_unknown_mapq_hits(
    tmp_path: Path,
) -> None:
    path = _paf(
        tmp_path / "alignments.paf",
        [
            "short\t8\t0\t4\t+\tchr1\t100\t20\t24\t4\t4\t60\ttp:A:P\tcg:Z:4M",
            "poor_mapq\t8\t0\t8\t+\tchr1\t100\t20\t28\t7\t8\t255\ttp:A:P\tcg:Z:8M",
        ],
    )
    queries = DNASet(
        [
            DNARecord(DNASequence("ACGTACGT"), "short"),
            DNARecord(DNASequence("ACGTACGT"), "poor_mapq"),
        ]
    )
    report = evaluate_reference_consistency(
        queries,
        path,
        reference_name="human",
        reference_version="GRCh38",
    )
    assert report.metrics["consistency_rate"] == pytest.approx((4 / 8 + 7 / 8) / 2)
    assert [entry.metrics["consistency_rate"] for entry in report.entries] == pytest.approx(
        [0.5, 0.875]
    )


def test_reference_consistency_rejects_stale_or_ambiguous_paf(tmp_path: Path) -> None:
    path = _paf(
        tmp_path / "alignments.paf",
        ["wrong\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P\tcg:Z:8M"],
    )
    with pytest.raises(ConfigurationError, match="Query ID"):
        evaluate_reference_consistency(
            DNASequence("ACGTACGT"), path, reference_name="human", reference_version="GRCh38"
        )
    duplicate_queries = DNASet(
        [
            DNARecord(DNASequence("ACGT"), "same"),
            DNARecord(DNASequence("ACGT"), "same"),
        ]
    )
    with pytest.raises(ConfigurationError) as error:
        evaluate_reference_consistency(
            duplicate_queries, path, reference_name="human", reference_version="GRCh38"
        )
    assert error.value.code == "DUPLICATE_REFERENCE_CONSISTENCY_ID"


def test_reference_consistency_empty_query_set_has_undefined_rate(tmp_path: Path) -> None:
    path = _paf(tmp_path / "empty.paf", [])
    report = evaluate_reference_consistency(
        DNASet([]), path, reference_name="human", reference_version="GRCh38"
    )
    assert report.metrics["consistency_rate"] is None


def test_reference_consistency_all_unmapped_queries_score_zero(tmp_path: Path) -> None:
    path = _paf(tmp_path / "empty.paf", [])
    report = evaluate_reference_consistency(
        _queries(), path, reference_name="human", reference_version="GRCh38"
    )
    assert report.metrics["consistency_rate"] == 0.0
    assert [entry.metrics["consistency_rate"] for entry in report.entries] == [0.0] * 3


def test_reference_consistency_requires_base_level_alignment(tmp_path: Path) -> None:
    path = _paf(
        tmp_path / "approximate.paf",
        ["sequence\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P"],
    )
    with pytest.raises(ConfigurationError) as error:
        evaluate_reference_consistency(
            DNASequence("ACGTACGT"), path, reference_name="human", reference_version="GRCh38"
        )
    assert error.value.code == "INVALID_REFERENCE_CONSISTENCY_PAF"


def test_reference_consistency_rejects_invalid_config(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        ReferenceConsistencyConfig(max_paf_bytes=0)
    with pytest.raises(ConfigurationError, match="reference_sha256"):
        evaluate_reference_consistency(
            DNASet([]),
            tmp_path / "missing.paf",
            reference_name="human",
            reference_version="GRCh38",
            reference_sha256="not-a-digest",
        )


def test_reference_consistency_downloads_then_reuses_verified_hg38(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = import_module("dnakit.evaluation.reference_consistency")
    download_calls: list[str] = []
    accession = "GCF_000001405.40"
    fasta_name = f"{accession}_GRCh38.p14_genomic.fna"
    genome = b">chr1\nACGTACGT\n"
    digest = hashlib.sha256(genome).hexdigest()

    def fake_download(query: str, directory: Path, *, show_progress: bool) -> SimpleNamespace:
        assert query == "hg38"
        download_calls.append(query)
        directory.mkdir(parents=True, exist_ok=True)
        fasta = directory / fasta_name
        fasta.write_bytes(genome)
        checksum = hashlib.md5(genome, usedforsecurity=False).hexdigest()
        (directory / f"{accession}_md5sum.txt").write_text(
            f"{checksum}  ./ncbi_dataset/data/{accession}/{fasta_name}\n",
            encoding="utf-8",
        )
        return SimpleNamespace(fasta_path=str(fasta), fasta_sha256=digest, accession=accession)

    def fake_execute(executable: str, args: tuple[str, ...], **kwargs: object) -> SimpleNamespace:
        assert "-c" in args
        assert args[-2].endswith(fasta_name)
        Path(args[args.index("-o") + 1]).write_text(
            "sequence\t4\t0\t4\t+\tchr1\t8\t0\t4\t4\t4\t60\ttp:A:P\tcg:Z:4M\n",
            encoding="ascii",
        )
        return SimpleNamespace(return_code=0, output="")

    monkeypatch.setattr(module, "download_genome", fake_download)
    monkeypatch.setattr(module.shutil, "which", lambda _: "/usr/bin/true")
    monkeypatch.setattr(module, "execute_bounded_command", fake_execute)
    directory = tmp_path / "references" / "hg38"

    first = evaluate_reference_consistency(DNASequence("ACGT"), reference_dir=directory)
    second = evaluate_reference_consistency(DNASequence("ACGT"), reference_dir=directory)

    assert download_calls == ["hg38"]
    assert first.metrics["consistency_rate"] == second.metrics["consistency_rate"] == 1.0
    assert first.parameters["reference_version"] == accession
    assert first.parameters["reference_sha256"] == digest
    assert first.parameters["alignment_source"] == "minimap2-generated"
    assert first.parameters["paf_path"] is None


def test_reference_consistency_uses_supplied_fasta_without_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = import_module("dnakit.evaluation.reference_consistency")
    fasta = tmp_path / "custom.fasta"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="ascii")

    def fake_execute(executable: str, args: tuple[str, ...], **kwargs: object) -> SimpleNamespace:
        assert args[-2] == str(fasta)
        Path(args[args.index("-o") + 1]).write_text(
            "sequence\t4\t0\t4\t+\tchr1\t8\t0\t4\t4\t4\t60\ttp:A:P\tcg:Z:4M\n",
            encoding="ascii",
        )
        return SimpleNamespace(return_code=0, output="")

    monkeypatch.setattr(module, "download_genome", lambda *a, **k: pytest.fail("downloaded"))
    monkeypatch.setattr(module.shutil, "which", lambda _: "/usr/bin/true")
    monkeypatch.setattr(module, "execute_bounded_command", fake_execute)

    report = evaluate_reference_consistency(
        DNASequence("ACGT"), reference_fasta=fasta, reference_name="custom"
    )
    assert report.metrics["consistency_rate"] == 1.0
    assert report.parameters["reference_name"] == "custom"
    assert str(report.parameters["reference_version"]).startswith("local:")


def test_missing_minimap2_fails_before_any_download(monkeypatch: pytest.MonkeyPatch) -> None:
    module = import_module("dnakit.evaluation.reference_consistency")
    monkeypatch.setattr(module.shutil, "which", lambda _: None)
    monkeypatch.setattr(module, "download_genome", lambda *a, **k: pytest.fail("downloaded"))
    with pytest.raises(BackendUnavailableError) as error:
        evaluate_reference_consistency(DNASequence("ACGT"))
    assert error.value.code == "MINIMAP2_UNAVAILABLE"
