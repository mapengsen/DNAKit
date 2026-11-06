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
    assert set(reference.metrics) == {
        "consistency_rate",
        "query_count",
        "aligned_count",
        "aligned_fraction",
        "weighted_blast_identity",
        "weighted_gap_compressed_identity",
    }
    assert reference.metrics["aligned_fraction"] == pytest.approx(2 / 3)
    assert reference.metrics["weighted_blast_identity"] == pytest.approx(15 / 16)
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
    assert [entry.metrics["aligned"] for entry in report.entries] == [True, True, False]
    assert [entry.metrics["best_target"] for entry in report.entries] == ["chr1", "chr2", None]
    assert report.entries[2].metrics["blast_identity"] is None
    assert report.entries[2].metrics["query_coverage"] == 0.0
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


def test_paper_example_blast_and_gap_compressed_identity(tmp_path: Path) -> None:
    """Huang et al. Fig. 4: 10 columns, 7 matches, one 2-base gap -> BI 7/10, GI 7/9."""

    path = _paf(
        tmp_path / "paper.paf",
        ["query\t8\t0\t8\t+\tref\t100\t20\t30\t7\t10\t60\tcg:Z:3=1X2D4="],
    )
    report = evaluate_reference_consistency(
        DNASet([DNARecord(DNASequence("ACGTACGT"), "query")]),
        path,
        reference_name="example",
        reference_version="1",
    )
    entry = report.entries[0].metrics
    assert entry["blast_identity"] == pytest.approx(0.7)
    assert entry["gap_compressed_identity"] == pytest.approx(7 / 9)
    assert entry["consistency_rate"] == pytest.approx(7 / 8)
    assert report.metrics["weighted_blast_identity"] == pytest.approx(0.7)
    assert report.metrics["weighted_gap_compressed_identity"] == pytest.approx(7 / 9)
    assert report.metrics["aligned_fraction"] == 1.0


def test_minimap2_m_cigar_uses_alignment_columns_and_gap_runs(tmp_path: Path) -> None:
    # 6M2I2M: 8 aligned pairs, 2 inserted query bases in one run, 7 matches.
    path = _paf(
        tmp_path / "m.paf",
        ["query\t10\t0\t10\t+\tchr1\t100\t0\t8\t7\t10\t60\ttp:A:P\tcg:Z:6M2I2M"],
    )
    report = evaluate_reference_consistency(
        DNASet([DNARecord(DNASequence("ACGTACGTAC"), "query")]),
        path,
        reference_name="human",
        reference_version="GRCh38",
    )
    entry = report.entries[0].metrics
    assert entry["blast_identity"] == pytest.approx(7 / 10)
    assert entry["gap_compressed_identity"] == pytest.approx(7 / 9)
    assert entry["consistency_rate"] == pytest.approx(7 / 10)


def test_wfmash_fragments_pool_per_target_and_pick_best_target(tmp_path: Path) -> None:
    # wfmash splits one 10 bp query into overlapping fragments on chrA (with =/X
    # CIGARs, no tp:A tag); a weaker single hit on chrB must not be chosen.
    path = _paf(
        tmp_path / "wfmash.paf",
        [
            "q\t10\t0\t6\t+\tchrA\t100\t0\t6\t5\t6\t60\tgi:f:0.83\tcg:Z:2=1X3=",
            "q\t10\t4\t10\t+\tchrA\t100\t4\t10\t6\t6\t60\tgi:f:1\tcg:Z:6=",
            "q\t10\t0\t10\t-\tchrB\t100\t0\t10\t6\t10\t60\tcg:Z:6=4X",
        ],
    )
    report = evaluate_reference_consistency(
        DNASet([DNARecord(DNASequence("ACGTACGTAC"), "q")]),
        path,
        reference_name="training",
        reference_version="1",
    )
    entry = report.entries[0].metrics
    assert entry["best_target"] == "chrA"
    # Matched query positions: 0-2 and 3-6 from the first fragment, 4-10 from the second.
    assert entry["consistency_rate"] == pytest.approx(9 / 10)
    assert entry["query_coverage"] == 1.0
    assert entry["blast_identity"] == pytest.approx(11 / 12)
    assert entry["gap_compressed_identity"] == pytest.approx(11 / 12)


def test_reverse_strand_eq_cigar_maps_matches_from_query_end(tmp_path: Path) -> None:
    path = _paf(
        tmp_path / "minus.paf",
        ["q\t10\t2\t8\t-\tchr1\t100\t0\t6\t4\t6\t60\tcg:Z:4=2X"],
    )
    report = evaluate_reference_consistency(
        DNASet([DNARecord(DNASequence("ACGTACGTAC"), "q")]),
        path,
        reference_name="r",
        reference_version="1",
    )
    assert report.entries[0].metrics["consistency_rate"] == pytest.approx(0.4)
    assert report.entries[0].metrics["query_coverage"] == pytest.approx(0.6)


def test_min_query_coverage_controls_aligned_fraction(tmp_path: Path) -> None:
    path = _paf(
        tmp_path / "coverage.paf",
        [
            "full\t8\t0\t8\t+\tchr1\t100\t0\t8\t8\t8\t60\tcg:Z:8M",
            "half\t8\t0\t4\t+\tchr1\t100\t0\t4\t3\t4\t60\tcg:Z:4M",
        ],
    )
    queries = DNASet(
        [
            DNARecord(DNASequence("ACGTACGT"), "full"),
            DNARecord(DNASequence("ACGTACGT"), "half"),
        ]
    )
    loose = evaluate_reference_consistency(queries, path, reference_name="r", reference_version="1")
    strict = evaluate_reference_consistency(
        queries,
        path,
        reference_name="r",
        reference_version="1",
        config=ReferenceConsistencyConfig(min_query_coverage=0.9),
    )
    assert loose.metrics["aligned_fraction"] == 1.0
    assert loose.metrics["weighted_blast_identity"] == pytest.approx(11 / 12)
    assert strict.metrics["aligned_fraction"] == 0.5
    assert strict.metrics["weighted_blast_identity"] == 1.0
    assert [entry.metrics["aligned"] for entry in strict.entries] == [True, False]
    # The consistency rate does not depend on the threshold.
    assert loose.metrics["consistency_rate"] == strict.metrics["consistency_rate"]


def test_reference_consistency_rejects_cigar_inconsistent_with_paf(tmp_path: Path) -> None:
    path = _paf(
        tmp_path / "bad.paf",
        ["q\t8\t0\t8\t+\tchr1\t100\t0\t8\t8\t8\t60\tcg:Z:6M"],
    )
    with pytest.raises(ConfigurationError) as error:
        evaluate_reference_consistency(
            DNASet([DNARecord(DNASequence("ACGTACGT"), "q")]),
            path,
            reference_name="r",
            reference_version="1",
        )
    assert error.value.code == "INVALID_REFERENCE_CONSISTENCY_PAF"


def test_wfmash_config_validation() -> None:
    assert ReferenceConsistencyConfig(aligner="wfmash", wfmash_args=("-w", "1000")).aligner == (
        "wfmash"
    )
    with pytest.raises(ConfigurationError):
        ReferenceConsistencyConfig(aligner="bwa")
    with pytest.raises(ConfigurationError):
        ReferenceConsistencyConfig(aligner="wfmash", wfmash_args=("--sam",))
    with pytest.raises(ConfigurationError):
        ReferenceConsistencyConfig(wfmash_args=["-w", "1000"])  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        ReferenceConsistencyConfig(min_query_coverage=1.5)


def test_wfmash_runs_on_indexed_reference_and_reads_stdout_paf(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = import_module("dnakit.evaluation.reference_consistency")
    fasta = tmp_path / "training.fa"
    fasta.write_text(">ref1 desc\nACGTAC\nGTAC\n>ref2\nTTTT\n", encoding="ascii")
    calls: list[tuple[str, ...]] = []

    def fake_execute(executable: str, args: tuple[str, ...], **kwargs: object) -> SimpleNamespace:
        calls.append(args)
        query = Path(args[1])
        assert Path(str(query) + ".fai").read_text(encoding="ascii") == "sequence\t4\t10\t4\t5\n"
        assert kwargs["backend_id"] == "wfmash"
        Path(str(kwargs["stdout_path"])).write_text(
            "sequence\t4\t0\t4\t+\tref1\t10\t0\t4\t4\t4\t255\tcg:Z:4=\n",
            encoding="ascii",
        )
        return SimpleNamespace(return_code=0, output="")

    monkeypatch.setattr(module.shutil, "which", lambda _: "/usr/bin/true")
    monkeypatch.setattr(module, "execute_bounded_command", fake_execute)
    report = evaluate_reference_consistency(
        DNASequence("ACGT"),
        reference_fasta=fasta,
        config=ReferenceConsistencyConfig(aligner="wfmash", wfmash_args=("-S", "0")),
    )
    assert calls[0][0] == str(fasta)
    assert calls[0][2:] == ("-t", "1", "--quiet", "-S", "0")
    assert (tmp_path / "training.fa.fai").read_text(encoding="ascii") == (
        "ref1\t10\t11\t6\t7\nref2\t4\t29\t4\t5\n"
    )
    assert report.metrics["consistency_rate"] == 1.0
    assert report.parameters["aligner"] == "wfmash"
    assert report.parameters["alignment_source"] == "wfmash-generated"
    assert tuple(report.parameters["wfmash_args"]) == ("-S", "0")


def test_missing_wfmash_reports_its_own_error(monkeypatch: pytest.MonkeyPatch) -> None:
    module = import_module("dnakit.evaluation.reference_consistency")
    monkeypatch.setattr(module.shutil, "which", lambda _: None)
    with pytest.raises(BackendUnavailableError) as error:
        evaluate_reference_consistency(
            DNASequence("ACGT"), config=ReferenceConsistencyConfig(aligner="wfmash")
        )
    assert error.value.code == "WFMASH_UNAVAILABLE"
