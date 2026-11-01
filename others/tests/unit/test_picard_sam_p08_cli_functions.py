"""Tests for dnakit.alignment.picard_sam_p08_cli_functions (Picard SAM/BAM and Illumina-run tools).

The tests run the real Picard (3.5.0 tested) on a small synthetic dataset. They skip when Picard
(`picard` on PATH, or `PICARD_JAR` plus `java`) is not installed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

import _picard_p08_data as data

from dnakit.alignment import picard_sam_p08_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

requires_picard = pytest.mark.skipif(
    not data.picard_available(),
    reason="Picard (picard on PATH, or PICARD_JAR with java) is not installed",
)


@pytest.fixture(scope="module")
def files(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return data.write_dataset(tmp_path_factory.mktemp("picard_sam"))


def _ok(result: dict) -> dict:
    assert result["status"] == "completed", result["log_tail"][-1500:]
    assert result["backend"] == "picard"
    return result


def _reads(path: Path) -> list[dict]:
    """Records of a SAM/BAM as dictionaries (name, flag, contig, position)."""
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [
            {
                "name": read.query_name,
                "flag": read.flag,
                "contig": read.reference_name,
                "position": read.reference_start,
                "tags": dict(read.tags),
            }
            for read in handle.fetch(until_eof=True)
        ]


def _unmapped_bam(files: dict[str, Path], folder: Path) -> Path:
    """The example reads reverted to an unaligned, queryname-sorted file."""
    out = folder / "unmapped.sam"
    _ok(m.picard_revert_sam(str(files["sam"]), str(out), verbose=False))
    return out


# ---------------------------------------------------------------------------
# Behaviour that needs no Picard
# ---------------------------------------------------------------------------


def test_missing_picard_gives_install_hint(files: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        m.picard_revert_sam(
            str(files["sam"]),
            str(tmp_path / "x.sam"),
            picard_path=str(tmp_path / "no_picard"),
            verbose=False,
        )


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (
            lambda f, t: m.picard_filter_sam_reads(
                str(f["sam"]), str(t / "o.sam"), "includeJavascript", verbose=False
            ),
            "filter",
        ),
        (
            lambda f, t: m.picard_filter_sam_reads(
                str(f["sam"]), str(t / "o.sam"), "includeReadList", verbose=False
            ),
            "read_list_file",
        ),
        (
            lambda f, t: m.picard_filter_sam_reads(
                str(f["sam"]), str(t / "o.sam"), "includeTagValues", tag="RG", verbose=False
            ),
            "tag_values",
        ),
        (
            lambda f, t: m.picard_position_based_downsample_sam(
                str(f["sam"]), str(t / "o.sam"), 1.5, verbose=False
            ),
            "fraction",
        ),
        (
            lambda f, t: m.picard_split_sam_by_number_of_reads(
                str(f["sam"]), str(t), verbose=False
            ),
            "split_to_n",
        ),
        (
            lambda f, t: m.picard_split_sam_by_number_of_reads(
                str(f["sam"]), str(t), split_to_n_files=2, split_to_n_reads=5, verbose=False
            ),
            "exactly one",
        ),
        (
            lambda f, t: m.picard_merge_bam_alignment(
                str(f["sam"]), str(f["fasta"]), str(t / "o.bam"), verbose=False
            ),
            "aligned_bams",
        ),
        (
            lambda f, t: m.picard_mark_duplicates_with_mate_cigar(
                [str(f["sam"])], str(t / "o.sam"), str(t / "m.txt"), umi_aware=True, verbose=False
            ),
            "umi_metrics_file",
        ),
        (
            lambda f, t: m.picard_mark_duplicates_with_mate_cigar(
                [str(f["sam"])],
                str(t / "o.sam"),
                str(t / "m.txt"),
                umi_tag_name="RX",
                verbose=False,
            ),
            "UMI options",
        ),
        (
            lambda f, t: m.picard_illumina_basecalls_to_fastq(
                str(t), [1], "2T", "run1", verbose=False
            ),
            "output_prefix",
        ),
        (
            lambda f, t: m.picard_illumina_basecalls_to_sam(
                str(t), [1], "2T", "run1", "XX", verbose=False
            ),
            "output_path",
        ),
        (
            lambda f, t: m.picard_mark_illumina_adapters(
                str(f["sam"]), str(t / "m.txt"), adapters=["NOT_AN_ADAPTER"], verbose=False
            ),
            "adapters",
        ),
    ],
)
def test_arguments_are_validated(
    files: dict[str, Path], tmp_path: Path, call, message: str
) -> None:
    with pytest.raises(ConfigurationError) as error:
        call(files, tmp_path)
    assert message in str(error.value)


# ---------------------------------------------------------------------------
# Subsetting and downsampling
# ---------------------------------------------------------------------------


@requires_picard
def test_downsample_sam_keeps_templates_together_and_is_reproducible(
    files: dict[str, Path], tmp_path: Path
) -> None:
    everything = _ok(
        m.picard_downsample_sam(
            str(files["sam"]), str(tmp_path / "all.sam"), probability=1.0, verbose=False
        )
    )
    assert len(_reads(Path(everything["outputs"]["output"]))) == 44
    half = _ok(
        m.picard_downsample_sam(
            str(files["sam"]),
            str(tmp_path / "half.sam"),
            probability=0.5,
            random_seed=11,
            metrics_file=str(tmp_path / "half.txt"),
            verbose=False,
        )
    )
    kept = _reads(tmp_path / "half.sam")
    names = [read["name"] for read in kept]
    assert 0 < len(kept) < 44
    assert all(names.count(name) == 2 for name in set(names)), "pairs are kept or dropped together"
    assert half["metrics"][0]["TOTAL_READS"] == len(kept)
    _ok(
        m.picard_downsample_sam(
            str(files["sam"]),
            str(tmp_path / "half2.sam"),
            probability=0.5,
            random_seed=11,
            verbose=False,
        )
    )
    assert [read["name"] for read in _reads(tmp_path / "half2.sam")] == names
    with pytest.raises(ConfigurationError):
        m.picard_downsample_sam(
            str(files["sam"]), str(tmp_path / "bad.sam"), strategy="Quick", verbose=False
        )


@requires_picard
def test_position_based_downsample_sam_is_deterministic(
    files: dict[str, Path], tmp_path: Path
) -> None:
    for name in ("a", "b"):
        _ok(
            m.picard_position_based_downsample_sam(
                str(files["sam"]), str(tmp_path / f"{name}.sam"), 0.5, verbose=False
            )
        )
    first = [read["name"] for read in _reads(tmp_path / "a.sam")]
    assert first == [read["name"] for read in _reads(tmp_path / "b.sam")]
    assert len(first) < 44 and len(first) % 2 == 0


@requires_picard
def test_filter_sam_reads_by_tag_name_and_aligned_state(
    files: dict[str, Path], tmp_path: Path
) -> None:
    by_tag = _ok(
        m.picard_filter_sam_reads(
            str(files["sam"]),
            str(tmp_path / "rg1.sam"),
            "includeTagValues",
            tag="RG",
            tag_values=["rg1"],
            verbose=False,
        )
    )
    assert len(_reads(Path(by_tag["outputs"]["output"]))) == 44
    none_left = _ok(
        m.picard_filter_sam_reads(
            str(files["sam"]),
            str(tmp_path / "rg2.sam"),
            "includeTagValues",
            tag="RG",
            tag_values=["rg2"],
            verbose=False,
        )
    )
    assert _reads(Path(none_left["outputs"]["output"])) == []
    names = tmp_path / "names.txt"
    wanted = _reads(files["sam"])[0]["name"]
    names.write_text(wanted + "\n")
    listed = _ok(
        m.picard_filter_sam_reads(
            str(files["sam"]),
            str(tmp_path / "one.sam"),
            "includeReadList",
            read_list_file=str(names),
            verbose=False,
        )
    )
    assert {read["name"] for read in _reads(Path(listed["outputs"]["output"]))} == {wanted}
    rest = _ok(
        m.picard_filter_sam_reads(
            str(files["sam"]),
            str(tmp_path / "rest.sam"),
            "excludeReadList",
            read_list_file=str(names),
            verbose=False,
        )
    )
    assert len(_reads(Path(rest["outputs"]["output"]))) == 44 - len(_reads(tmp_path / "one.sam"))
    unmapped = _unmapped_bam(files, tmp_path)
    aligned = _ok(
        m.picard_filter_sam_reads(
            str(unmapped), str(tmp_path / "aln.sam"), "includeAligned", verbose=False
        )
    )
    assert _reads(Path(aligned["outputs"]["output"])) == []


@requires_picard
def test_split_sam_by_number_of_reads_keeps_every_read(
    files: dict[str, Path], tmp_path: Path
) -> None:
    unmapped = _unmapped_bam(files, tmp_path)
    folder = tmp_path / "shards"
    folder.mkdir()
    result = _ok(
        m.picard_split_sam_by_number_of_reads(
            str(unmapped), str(folder), split_to_n_files=3, verbose=False
        )
    )
    assert [Path(path).name for path in result["files"]] == [
        "shard_0001.sam",
        "shard_0002.sam",
        "shard_0003.sam",
    ]
    total = [read["name"] for path in result["files"] for read in _reads(Path(path))]
    assert sorted(total) == sorted(read["name"] for read in _reads(unmapped))
    with pytest.raises(ConfigurationError):
        m.picard_split_sam_by_number_of_reads(
            str(unmapped), str(folder), split_to_n_reads=10, verbose=False
        )
    other = tmp_path / "more"
    other.mkdir()
    by_reads = _ok(
        m.picard_split_sam_by_number_of_reads(
            str(unmapped), str(other), split_to_n_reads=20, out_prefix="part", verbose=False
        )
    )
    assert Path(by_reads["files"][0]).name == "part_0001.sam"
    assert len(by_reads["files"]) == 3  # 44 reads, about 20 per file


# ---------------------------------------------------------------------------
# Reshaping alignments
# ---------------------------------------------------------------------------


@requires_picard
def test_revert_sam_unmaps_reads_and_clears_alignment_tags(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(m.picard_revert_sam(str(files["sam"]), str(tmp_path / "u.sam"), verbose=False))
    reads = _reads(Path(result["outputs"]["output"]))
    assert len(reads) == 44
    assert all(read["flag"] & 4 for read in reads), "every read is unmapped"
    assert all("NM" not in read["tags"] and "MC" not in read["tags"] for read in reads)
    assert all(read["tags"].get("RG") == "rg1" for read in reads)
    keep = _ok(
        m.picard_revert_sam(
            str(files["sam"]),
            str(tmp_path / "kept.sam"),
            remove_alignment_information=False,
            restore_hardclips=False,
            sort_order="coordinate",
            verbose=False,
        )
    )
    assert any(not read["flag"] & 4 for read in _reads(Path(keep["outputs"]["output"])))
    folder = tmp_path / "per_rg"
    folder.mkdir()
    grouped = _ok(
        m.picard_revert_sam(str(files["sam"]), str(folder), output_by_readgroup=True, verbose=False)
    )
    assert [Path(path).name for path in grouped["files"]] == ["rg1.sam"]


@requires_picard
def test_reorder_sam_follows_the_new_contig_order(files: dict[str, Path], tmp_path: Path) -> None:
    dictionary = tmp_path / "swapped.dict"
    dictionary.write_text("@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:800\n@SQ\tSN:chr1\tLN:1200\n")
    result = _ok(
        m.picard_reorder_sam(
            str(files["sam"]), str(tmp_path / "re.sam"), str(dictionary), verbose=False
        )
    )
    header = [
        line.split("\t")[1]
        for line in Path(result["outputs"]["output"]).read_text().splitlines()
        if line.startswith("@SQ")
    ]
    assert header == ["SN:chr2", "SN:chr1"]
    contigs = [read["contig"] for read in _reads(Path(result["outputs"]["output"]))]
    assert contigs == sorted(contigs, key=["chr2", "chr1"].index)
    short = tmp_path / "chr1_only.dict"
    short.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1200\n")
    refused = m.picard_reorder_sam(
        str(files["sam"]), str(tmp_path / "re2.sam"), str(short), verbose=False
    )
    assert refused["status"] == "failed"
    allowed = _ok(
        m.picard_reorder_sam(
            str(files["sam"]),
            str(tmp_path / "re3.sam"),
            str(short),
            allow_incomplete_dict_concordance=True,
            verbose=False,
        )
    )
    assert any(read["flag"] & 4 for read in _reads(Path(allowed["outputs"]["output"])))


@requires_picard
def test_mark_illumina_adapters_writes_a_clipping_histogram(
    files: dict[str, Path], tmp_path: Path
) -> None:
    unmapped = _unmapped_bam(files, tmp_path)
    result = _ok(
        m.picard_mark_illumina_adapters(
            str(unmapped),
            str(tmp_path / "adapters.txt"),
            output_path=str(tmp_path / "marked.sam"),
            adapters=["PAIRED_END"],
            verbose=False,
        )
    )
    assert len(_reads(Path(result["outputs"]["output"]))) == 44
    only_metrics = _ok(
        m.picard_mark_illumina_adapters(str(unmapped), str(tmp_path / "only.txt"), verbose=False)
    )
    assert set(only_metrics["outputs"]) == {"metrics"}


@requires_picard
def test_merge_bam_alignment_restores_unmapped_metadata(
    files: dict[str, Path], tmp_path: Path
) -> None:
    unmapped = _unmapped_bam(files, tmp_path)
    result = _ok(
        m.picard_merge_bam_alignment(
            str(unmapped),
            str(files["fasta"]),
            str(tmp_path / "merged.sam"),
            aligned_bams=[str(files["sam"])],
            verbose=False,
        )
    )
    merged = _reads(Path(result["outputs"]["output"]))
    original = {(read["name"], read["flag"] & 192): read for read in _reads(files["sam"])}
    assert len(merged) == 44
    for read in merged:
        match = original[(read["name"], read["flag"] & 192)]
        assert (read["contig"], read["position"]) == (match["contig"], match["position"])
        assert read["tags"]["RG"] == "rg1"
    with pytest.raises(ConfigurationError):
        m.picard_merge_bam_alignment(
            str(unmapped),
            str(files["fasta"]),
            str(tmp_path / "x.sam"),
            aligned_bams=[str(files["sam"])],
            read1_aligned_bams=[str(files["sam"])],
            verbose=False,
        )


# ---------------------------------------------------------------------------
# Duplicate marking
# ---------------------------------------------------------------------------


@requires_picard
def test_mark_duplicates_with_mate_cigar_flags_the_two_duplicated_pairs(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_mark_duplicates_with_mate_cigar(
            [str(files["sam"])],
            str(tmp_path / "marked.sam"),
            str(tmp_path / "dup.txt"),
            verbose=False,
        )
    )
    row = result["metrics"][0]
    assert row["READ_PAIRS_EXAMINED"] == 22 and row["READ_PAIR_DUPLICATES"] == 2
    flagged = [read for read in _reads(Path(result["outputs"]["output"])) if read["flag"] & 1024]
    assert len(flagged) == 4  # both reads of the two duplicated pairs
    removed = _ok(
        m.picard_mark_duplicates_with_mate_cigar(
            [str(files["sam"])],
            str(tmp_path / "removed.sam"),
            str(tmp_path / "dup2.txt"),
            remove_duplicates=True,
            verbose=False,
        )
    )
    assert len(_reads(Path(removed["outputs"]["output"]))) == 40
    # The same duplicates are counted by CollectDuplicateMetrics on the marked file.
    from dnakit.sequencing_qc import picard_metrics_p08_cli_functions as metrics

    again = _ok(
        metrics.picard_duplicate_metrics(
            str(tmp_path / "marked.sam"), str(tmp_path / "again.txt"), verbose=False
        )
    )
    assert again["metrics"][0]["READ_PAIR_DUPLICATES"] == 2


@requires_picard
def test_umi_aware_duplicate_marking_reports_umi_metrics(
    files: dict[str, Path], tmp_path: Path
) -> None:
    lines = []
    for line in files["sam"].read_text().splitlines():
        if line.startswith("@"):
            lines.append(line)
        else:
            umi = "ACGT" if line.split("\t")[0].endswith(("1001:2001", "1004:2004")) else "TTTT"
            lines.append(f"{line}\tRX:Z:{umi}")
    umi_sam = tmp_path / "umi.sam"
    umi_sam.write_text("\n".join(lines) + "\n")
    result = _ok(
        m.picard_mark_duplicates_with_mate_cigar(
            [str(umi_sam)],
            str(tmp_path / "umi_marked.sam"),
            str(tmp_path / "umi_dup.txt"),
            umi_aware=True,
            umi_metrics_file=str(tmp_path / "umi_metrics.txt"),
            verbose=False,
        )
    )
    assert result["umi_metrics"][0]["LIBRARY"] == "lib1"
    assert result["umi_metrics"][0]["MEAN_UMI_LENGTH"] == 4
    assert result["metrics"][0]["READ_PAIRS_EXAMINED"] == 22


# ---------------------------------------------------------------------------
# Illumina run folders
# ---------------------------------------------------------------------------


@requires_picard
def test_illumina_basecalls_to_fastq_decodes_the_bcl_files(tmp_path: Path) -> None:
    run = data.write_illumina_run(tmp_path / "run")
    out = tmp_path / "fastq"
    out.mkdir()
    result = _ok(
        m.picard_illumina_basecalls_to_fastq(
            str(run),
            [1],
            "2T",
            "run1",
            output_prefix=str(out / "s"),
            machine_name="M1",
            flowcell_barcode="FC1",
            verbose=False,
        )
    )
    assert [Path(path).name for path in result["files"]] == ["s.1.fastq"]
    records = Path(result["files"][0]).read_text().splitlines()
    sequences = records[1::4]
    # Cycle 1 calls A, C, G, T and cycle 2 calls T, G, C, A (see write_illumina_run).
    assert sequences == ["AT", "CG", "GC", "TA"]
    assert records[0].startswith("@M1:run1:FC1:1:1101:")
    assert set("".join(records[3::4])) == {"?"}  # quality 30 + 33
    with pytest.raises(ConfigurationError):  # existing FASTQ files are never overwritten
        m.picard_illumina_basecalls_to_fastq(
            str(run),
            [1],
            "2T",
            "run1",
            output_prefix=str(out / "s"),
            machine_name="M1",
            flowcell_barcode="FC1",
            verbose=False,
        )
    plain = _ok(
        m.picard_illumina_basecalls_to_fastq(
            str(run),
            [1],
            "2T",
            "run1",
            output_prefix=str(out / "plain"),
            read_name_format="ILLUMINA",
            compress_outputs=True,
            verbose=False,
        )
    )
    assert [Path(path).name for path in plain["files"]] == ["plain.1.fastq.gz"]


@requires_picard
def test_illumina_basecalls_to_sam_writes_an_unmapped_bam(tmp_path: Path) -> None:
    run = data.write_illumina_run(tmp_path / "run")
    out = tmp_path / "bam"
    out.mkdir()
    result = _ok(
        m.picard_illumina_basecalls_to_sam(
            str(run),
            [1],
            "2T",
            "run1",
            "XX",
            output_path=str(out / "u.bam"),
            sample_alias="s1",
            library_name="l1",
            verbose=False,
        )
    )
    reads = _reads(Path(result["files"][0]))
    assert len(reads) == 4 and all(read["flag"] & 4 for read in reads)
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(result["files"][0], check_sq=False) as handle:
        group = handle.header.to_dict()["RG"][0]
        sequences = [read.query_sequence for read in handle.fetch(until_eof=True)]
    assert group["SM"] == "s1" and group["LB"] == "l1" and group["CN"] == "XX"
    assert sorted(sequences) == ["AT", "CG", "GC", "TA"]
