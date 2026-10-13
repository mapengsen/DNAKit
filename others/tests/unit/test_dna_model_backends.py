"""Contract tests for DNA-level adapters of DNA AI model packages.

Real upstream calls cover alphagenome, grelu, selene-sdk, enformer-pytorch and
genomic-benchmarks. evo2 and kipoi are exercised with fake modules because evo2
cannot be imported here and kipoi needs git clones or network access.
"""

from __future__ import annotations

import gzip
import inspect
import os
import random
import sys
import types
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.predictions import dna_model_backends as dmb

FASTA_TEXT = ">chr1\nACGTACGTNN\n>chr2\nGGGCCCAATT\n"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _gtf_line(chrom: str, feature: str, start: int, end: int, strand: str, frame: str, attrs: str):
    return "\t".join([chrom, "test", feature, str(start), str(end), ".", strand, frame, attrs])


def _gtf_text() -> str:
    gene = 'gene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";'
    tx = f'{gene} transcript_id "T1"; transcript_type "protein_coding";'
    lines = [
        _gtf_line("chr1", "gene", 101, 500, "+", ".", gene),
        _gtf_line("chr1", "transcript", 101, 500, "+", ".", tx),
        _gtf_line("chr1", "exon", 101, 200, "+", ".", tx),
        _gtf_line("chr1", "exon", 301, 500, "+", ".", tx),
        _gtf_line("chr1", "CDS", 151, 200, "+", "0", tx),
        _gtf_line("chr1", "CDS", 301, 450, "+", "2", tx),
        _gtf_line(
            "chr1",
            "gene",
            700,
            800,
            "-",
            ".",
            'gene_id "G2"; gene_name "GENE2"; gene_type "lincRNA";',
        ),
    ]
    return "\n".join(lines) + "\n"


def _fasta(tmp_path: Path) -> Path:
    return _write(tmp_path, "tiny.fa", FASTA_TEXT)


def _bed(tmp_path: Path, rows: list[str]) -> Path:
    return _write(tmp_path, "windows.bed", "\n".join(rows) + "\n")


def _dataset(root: Path) -> Path:
    """Write a two-class genomic-benchmarks dataset with 2 train and 2 test sequences per class."""

    folder = root / "demo_ds"
    (folder / "train").mkdir(parents=True)
    (folder / "test").mkdir(parents=True)
    (folder / "metadata.yaml").write_text("version: 0\nclasses:\n  - neg\n  - pos\n")
    for split in ("train", "test"):
        for label in ("neg", "pos"):
            payload = "id,region,start,end\n1,chr1,10,20\n2,chr1,30,45\n"
            with gzip.open(folder / split / f"{label}.csv.gz", "wt") as handle:
                handle.write(payload)
    return folder


# Module surface


EXPECTED_PUBLIC = {
    "check_genomic_benchmark_cache",
    "describe_genomic_benchmark",
    "describe_kipoi_model",
    "encode_one_hot_enformer",
    "encode_one_hot_selene",
    "extract_gencode_transcripts",
    "fetch_genome_sequences",
    "generate_evo2_sequences",
    "list_genomic_benchmark_datasets",
    "list_kipoi_models",
    "merge_alphagenome_intervals",
    "mutate_dna_sequence",
    "parse_alphagenome_interval",
    "parse_alphagenome_variant",
    "read_contig_lengths",
    "read_enformer_bed_windows",
    "read_gencode_gene_intervals",
    "resize_dna_sequences",
    "reverse_complement_with_grelu",
    "shuffle_dinucleotides",
}


def test_public_surface_is_exact_and_accepts_verbose() -> None:
    assert set(dmb.__all__) == EXPECTED_PUBLIC
    for name in dmb.__all__:
        assert "verbose" in inspect.signature(getattr(dmb, name)).parameters, name


@pytest.mark.parametrize(
    ("module", "extra", "call"),
    [
        (
            "alphagenome",
            "bio-alphagenome",
            lambda: dmb.parse_alphagenome_interval("chr1:1-2", verbose=False),
        ),
        ("grelu", "bio-grelu", lambda: dmb.reverse_complement_with_grelu("ACGT", verbose=False)),
        ("selene_sdk", "bio-selene", lambda: dmb.encode_one_hot_selene("ACGT", verbose=False)),
        (
            "enformer_pytorch",
            "neural-enformer",
            lambda: dmb.encode_one_hot_enformer("ACGT", verbose=False),
        ),
        (
            "genomic_benchmarks",
            "bio-genomic-benchmarks",
            lambda: dmb.list_genomic_benchmark_datasets(verbose=False),
        ),
        (
            "evo2",
            "neural-evo2",
            lambda: dmb.generate_evo2_sequences(
                ["ACGT"], weights_path=Path(__file__), verbose=False
            ),
        ),
        ("kipoi", "bio-kipoi", lambda: dmb.list_kipoi_models(verbose=False)),
    ],
)
def test_missing_package_raises_install_hint(
    monkeypatch: pytest.MonkeyPatch, module: str, extra: str, call
) -> None:
    # Block already-imported submodules too, otherwise require() can reuse them from the cache.
    for name in list(sys.modules):
        if name == module or name.startswith(module + "."):
            monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.raises(BackendUnavailableError, match=extra):
        call()


# Validation (no upstream calls needed)


def test_dna_inputs_are_checked_before_upstream_calls() -> None:
    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        dmb.reverse_complement_with_grelu("ACGX", verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty"):
        dmb.reverse_complement_with_grelu([], verbose=False)
    with pytest.raises(ConfigurationError, match="NUL"):
        dmb.reverse_complement_with_grelu("AC\x00GT", verbose=False)
    with pytest.raises(ConfigurationError, match="integer"):
        dmb.resize_dna_sequences(["ACGT"], True, verbose=False)
    with pytest.raises(ConfigurationError, match="n_shuffles"):
        dmb.shuffle_dinucleotides("AACCGGTT", n_shuffles=0, verbose=False)
    with pytest.raises(ConfigurationError, match="position"):
        dmb.mutate_dna_sequence("ACGT", "T", position=4, verbose=False)
    with pytest.raises(ConfigurationError, match="end must be one of"):
        dmb.resize_dna_sequences(["ACGT"], 5, end="middle", verbose=False)
    with pytest.raises(ConfigurationError, match="at most 5 characters"):
        dmb.encode_one_hot_selene("ACGT" * 10, max_length=5, verbose=False)


def test_alphagenome_gene_argument_rules(tmp_path: Path) -> None:
    gtf = _write(tmp_path, "genes.gtf", _gtf_text())
    with pytest.raises(ConfigurationError, match="exactly one"):
        dmb.read_gencode_gene_intervals(gtf, verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        dmb.read_gencode_gene_intervals(gtf, gene_symbols=["GENE1"], gene_ids=["G1"], verbose=False)


def test_evo2_validation_rejects_bad_settings_without_loading_weights(tmp_path: Path) -> None:
    weights = _write(tmp_path, "weights.pt", "placeholder")
    with pytest.raises(ConfigurationError, match="does not exist"):
        dmb.generate_evo2_sequences(["ACGT"], weights_path=tmp_path / "missing.pt", verbose=False)
    with pytest.raises(ConfigurationError, match="model_name"):
        dmb.generate_evo2_sequences(
            ["ACGT"], weights_path=weights, model_name="evo2_99b", verbose=False
        )
    with pytest.raises(ConfigurationError, match="temperature"):
        dmb.generate_evo2_sequences(["ACGT"], weights_path=weights, temperature=0.0, verbose=False)


# alphagenome: real upstream calls


def test_parse_alphagenome_interval_matches_direct_upstream(tmp_path: Path) -> None:
    genome = pytest.importorskip("alphagenome.data.genome")
    result = dmb.parse_alphagenome_interval("chr1:101-200:+", verbose=False)
    direct = genome.Interval.from_str("chr1:101-200:+")
    assert result["interval"] == {
        "chromosome": "chr1",
        "start": 101,
        "end": 200,
        "strand": "+",
        "width": 99,
    }
    assert result["interval"]["start"] == direct.start
    assert result["interval"]["end"] == direct.end
    with pytest.raises(ConfigurationError, match="rejected the input"):
        dmb.parse_alphagenome_interval("chr1:10-5", verbose=False)


def test_parse_alphagenome_variant_snv_and_deletion() -> None:
    genome = pytest.importorskip("alphagenome.data.genome")
    snv = dmb.parse_alphagenome_variant("chr1:10:A>G", verbose=False)["variant"]
    # VCF position 10 is 0-based start 9 and exclusive end 10.
    assert (snv["position"], snv["start"], snv["end"]) == (10, 9, 10)
    assert snv["is_snv"] is True and snv["is_deletion"] is False
    deletion = dmb.parse_alphagenome_variant("chr1:10:AT>A", verbose=False)["variant"]
    assert deletion["is_deletion"] is True and deletion["is_indel"] is True
    direct = genome.Variant.from_str("chr1:10:A>G")
    assert snv["start"] == direct.start and snv["end"] == direct.end
    with pytest.raises(ConfigurationError, match="variant_format"):
        dmb.parse_alphagenome_variant("chr1:10:A>G", variant_format="bogus", verbose=False)


def test_merge_alphagenome_intervals_matches_hand_merge() -> None:
    pytest.importorskip("alphagenome.data.genome")
    result = dmb.merge_alphagenome_intervals(
        ["chr1:0-10", "chr1:5-20", "chr1:30-40"], verbose=False
    )
    assert result["input_count"] == 3 and result["merged_count"] == 2
    assert [(item["start"], item["end"]) for item in result["intervals"]] == [(0, 20), (30, 40)]
    with pytest.raises(ConfigurationError, match="intervals must be"):
        dmb.merge_alphagenome_intervals([], verbose=False)


def test_read_gencode_gene_intervals_converts_to_zero_based(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome.data.gene_annotation")
    gtf = _write(tmp_path, "genes.gtf", _gtf_text())
    by_symbol = dmb.read_gencode_gene_intervals(gtf, gene_symbols=["GENE1"], verbose=False)
    # GTF 1-based 101..500 is 0-based start 100 with exclusive end 500.
    assert by_symbol["genes"] == [
        {
            "chromosome": "chr1",
            "start": 100,
            "end": 500,
            "strand": "+",
            "width": 400,
            "name": "GENE1",
        }
    ]
    by_id = dmb.read_gencode_gene_intervals(gtf, gene_ids=["G2"], verbose=False)
    assert (by_id["genes"][0]["start"], by_id["genes"][0]["end"]) == (699, 800)
    assert by_id["genes"][0]["strand"] == "-"


def test_read_gencode_gene_intervals_enforces_record_limit(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome.data.gene_annotation")
    gtf = _write(tmp_path, "genes.gtf", _gtf_text())
    with pytest.raises(ConfigurationError, match="max_records=1"):
        dmb.read_gencode_gene_intervals(gtf, gene_symbols=["GENE1"], max_records=1, verbose=False)


def test_extract_gencode_transcripts_matches_direct_upstream(tmp_path: Path) -> None:
    pd = pytest.importorskip("pandas")
    genome = pytest.importorskip("alphagenome.data.genome")
    transcript = pytest.importorskip("alphagenome.data.transcript")
    gtf = _write(tmp_path, "tx.gtf", _gtf_text())
    result = dmb.extract_gencode_transcripts(gtf, "chr1:0-1000", verbose=False)
    assert result["transcript_count"] == 1
    record = result["transcripts"][0]
    assert record["transcript_id"] == "T1" and record["gene_id"] == "G1"
    assert record["exons"] == [[100, 200], [300, 500]]
    assert record["strand"] == "+" and record["is_coding"] is True

    # Direct upstream call on the same rows; cds comes from Transcript.cds.
    rows = [
        {
            "Chromosome": "chr1",
            "Start": 150,
            "End": 200,
            "Strand": "+",
            "Feature": "CDS",
            "Frame": 0,
            "transcript_id": "T1",
            "gene_id": "G1",
            "transcript_type": "protein_coding",
            "gene_type": None,
            "protein_id": None,
            "gene_name": None,
        },
        {
            "Chromosome": "chr1",
            "Start": 300,
            "End": 450,
            "Strand": "+",
            "Feature": "CDS",
            "Frame": 2,
            "transcript_id": "T1",
            "gene_id": "G1",
            "transcript_type": "protein_coding",
            "gene_type": None,
            "protein_id": None,
            "gene_name": None,
        },
    ]
    base = dict(
        Chromosome="chr1",
        Strand="+",
        transcript_id="T1",
        gene_id="G1",
        transcript_type="protein_coding",
        gene_type=None,
        protein_id=None,
        gene_name=None,
        Frame=None,
    )
    frame = pd.DataFrame(
        [
            {**base, "Start": 100, "End": 500, "Feature": "transcript"},
            {**base, "Start": 100, "End": 200, "Feature": "exon"},
            {**base, "Start": 300, "End": 500, "Feature": "exon"},
            *rows,
        ]
    )
    extractor = transcript.TranscriptExtractor(frame)
    direct = extractor.extract(genome.Interval.from_str("chr1:0-1000"))
    assert len(direct) == 1
    assert record["cds"] == [[int(part.start), int(part.end)] for part in direct[0].cds]

    empty = dmb.extract_gencode_transcripts(gtf, "chr1:600-700", verbose=False)
    assert empty["transcript_count"] == 0


def test_gtf_reader_rejects_malformed_rows(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome.data.transcript")
    gtf = _write(tmp_path, "bad.gtf", "chr1\ttest\texon\t101\t200\n")
    with pytest.raises(ConfigurationError, match="fewer than 9 columns"):
        dmb.extract_gencode_transcripts(gtf, "chr1:0-1000", verbose=False)
    cds = _write(
        tmp_path,
        "badframe.gtf",
        _gtf_line("chr1", "CDS", 151, 200, "+", ".", 'transcript_id "T1";') + "\n",
    )
    with pytest.raises(ConfigurationError, match="numeric frame"):
        dmb.extract_gencode_transcripts(cds, "chr1:0-1000", verbose=False)


# grelu: real upstream calls


def test_reverse_complement_with_grelu_matches_hand_values() -> None:
    pytest.importorskip("grelu.sequence.utils")
    assert dmb.reverse_complement_with_grelu(["AAC", "GGT"], verbose=False)["sequences"] == [
        "GTT",
        "ACC",
    ]
    assert dmb.reverse_complement_with_grelu("ACGTN", verbose=False)["sequences"] == ["NACGT"]


def test_reverse_complement_is_an_involution() -> None:
    pytest.importorskip("grelu.sequence.utils")
    rng = random.Random(0)
    sequence = "".join(rng.choice("ACGT") for _ in range(200))
    once = dmb.reverse_complement_with_grelu(sequence, verbose=False)["sequences"][0]
    twice = dmb.reverse_complement_with_grelu(once, verbose=False)["sequences"][0]
    assert twice == sequence


def test_shuffle_dinucleotides_preserves_pairs_and_matches_upstream() -> None:
    utils = pytest.importorskip("grelu.sequence.utils")
    sequence = "AACCGGTTAGCTAGGATCCAGTAC"
    result = dmb.shuffle_dinucleotides(sequence, n_shuffles=3, seed=7, verbose=False)
    assert len(result["shuffles"]) == 3
    pairs = Counter(sequence[i : i + 2] for i in range(len(sequence) - 1))
    for shuffled in result["shuffles"]:
        assert len(shuffled) == len(sequence)
        assert Counter(shuffled[i : i + 2] for i in range(len(shuffled) - 1)) == pairs
        assert shuffled[0] == sequence[0] and shuffled[-1] == sequence[-1]
    direct = utils.dinuc_shuffle(sequence, n_shuffles=3, seed=7)
    assert result["shuffles"] == [str(item) for item in direct]


def test_shuffle_dinucleotides_rejects_low_diversity() -> None:
    # Upstream returns a single shuffle for one request, but raises when several
    # requested shuffles are all identical.
    pytest.importorskip("grelu.sequence.utils")
    assert dmb.shuffle_dinucleotides("ACGTAC", n_shuffles=1, verbose=False)["shuffles"] == [
        "ACGTAC"
    ]
    with pytest.raises(ConfigurationError, match="rejected the input"):
        dmb.shuffle_dinucleotides("ACGTAC", n_shuffles=2, seed=3, verbose=False)


def test_resize_dna_sequences_trims_and_pads() -> None:
    pytest.importorskip("grelu.sequence.utils")
    assert dmb.resize_dna_sequences(["ACGT"], 6, verbose=False)["sequences"] == ["NACGTN"]
    assert dmb.resize_dna_sequences(["ACGTAC"], 4, verbose=False)["sequences"] == ["CGTA"]
    assert dmb.resize_dna_sequences(["ACGT"], 5, end="left", verbose=False)["sequences"] == [
        "NACGT"
    ]
    assert dmb.resize_dna_sequences(["ACGT"], 5, end="right", verbose=False)["sequences"] == [
        "ACGTN"
    ]


def test_mutate_dna_sequence_replaces_one_base() -> None:
    pytest.importorskip("grelu.sequence.mutate")
    result = dmb.mutate_dna_sequence("ACGT", "T", position=1, verbose=False)
    assert result["sequence"] == "ATGT"
    assert result["reference_base"] == "C" and result["position"] == 1
    assert len(result["sequence"]) == len(result["reference"])


# selene-sdk: real upstream calls


def test_fetch_genome_sequences_matches_hand_values_and_direct_genome(tmp_path: Path) -> None:
    selene = pytest.importorskip("selene_sdk.sequences")
    fasta = _fasta(tmp_path)
    regions = [
        {"chrom": "chr1", "start": 0, "end": 4},
        {"chrom": "chr1", "start": 2, "end": 6, "strand": "-"},
        {"chrom": "chr1", "start": 8, "end": 12},
    ]
    result = dmb.fetch_genome_sequences(fasta, regions, pad=True, verbose=False)
    sequences = [item["sequence"] for item in result["regions"]]
    # chr1 is ACGTACGTNN; the last region is two N bases plus two N padding bases.
    assert sequences == ["ACGT", "GTAC", "NNNN"]
    direct = selene.Genome(str(fasta)).get_sequence_from_coords("chr1", 8, 12, pad=True)
    assert sequences[2] == direct


def test_fetch_genome_sequences_rejects_bad_regions(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk.sequences")
    fasta = _fasta(tmp_path)
    with pytest.raises(ConfigurationError, match="outside the FASTA"):
        dmb.fetch_genome_sequences(fasta, [{"chrom": "chr1", "start": 8, "end": 12}], verbose=False)
    with pytest.raises(ConfigurationError, match="unknown chromosome"):
        dmb.fetch_genome_sequences(fasta, [{"chrom": "chrZ", "start": 0, "end": 4}], verbose=False)
    with pytest.raises(ConfigurationError, match="greater than start"):
        dmb.fetch_genome_sequences(fasta, [{"chrom": "chr1", "start": 4, "end": 4}], verbose=False)
    with pytest.raises(ConfigurationError, match="max_region_length"):
        dmb.fetch_genome_sequences(
            fasta,
            [{"chrom": "chr1", "start": 0, "end": 9}],
            max_region_length=5,
            verbose=False,
        )


def test_read_contig_lengths_lists_each_contig(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk.sequences")
    result = dmb.read_contig_lengths(_fasta(tmp_path), verbose=False)
    assert result["contigs"] == [
        {"name": "chr1", "length": 10},
        {"name": "chr2", "length": 10},
    ]
    assert result["contig_count"] == 2


def test_encode_one_hot_selene_uses_uniform_unknown_base() -> None:
    selene = pytest.importorskip("selene_sdk.sequences")
    result = dmb.encode_one_hot_selene("ACGTN", verbose=False)
    assert result["one_hot"] == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.25, 0.25, 0.25, 0.25],
    ]
    direct = selene.Genome.sequence_to_encoding("ACGTN").tolist()
    assert result["one_hot"] == direct


# enformer-pytorch: real upstream calls


def test_encode_one_hot_enformer_uses_zero_unknown_base() -> None:
    data = pytest.importorskip("enformer_pytorch.data")
    result = dmb.encode_one_hot_enformer("ACGTN", verbose=False)
    assert result["one_hot"][4] == [0.0, 0.0, 0.0, 0.0]
    assert result["one_hot"][:4] == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    assert result["one_hot"] == data.str_to_one_hot("ACGTN").tolist()


def test_read_enformer_bed_windows_counts_padding_and_unknown(tmp_path: Path) -> None:
    pytest.importorskip("enformer_pytorch.data")
    fasta = _fasta(tmp_path)
    bed = _bed(tmp_path, ["chr1\t0\t4", "chr1\t6\t10"])
    result = dmb.read_enformer_bed_windows(bed, fasta, context_length=8, verbose=False)
    first, second = result["windows"]
    # Length-4 intervals are extended to 8 with 2 padding bases on each side.
    assert (first["length"], first["padding_positions"], first["unknown_positions"]) == (8, 2, 0)
    assert first["acgt_positions"] == 6
    assert (second["length"], second["padding_positions"], second["unknown_positions"]) == (8, 2, 2)
    assert second["acgt_positions"] == 4
    assert result["window_count"] == 2


def test_read_enformer_bed_windows_keeps_long_windows_at_bed_length(tmp_path: Path) -> None:
    pytest.importorskip("enformer_pytorch.data")
    fasta = _fasta(tmp_path)
    bed = _bed(tmp_path, ["chr2\t0\t10"])
    window = dmb.read_enformer_bed_windows(bed, fasta, context_length=8, verbose=False)["windows"][
        0
    ]
    assert window["length"] == 10 and window["padding_positions"] == 0


def test_read_enformer_bed_windows_enforces_window_limit(tmp_path: Path) -> None:
    pytest.importorskip("enformer_pytorch.data")
    fasta = _fasta(tmp_path)
    bed = _bed(tmp_path, ["chr1\t0\t4", "chr1\t6\t10"])
    with pytest.raises(ConfigurationError, match="between 1 and 1"):
        dmb.read_enformer_bed_windows(bed, fasta, context_length=8, max_windows=1, verbose=False)


# genomic-benchmarks: real offline calls


def test_list_genomic_benchmark_datasets_matches_upstream_index() -> None:
    gb = pytest.importorskip("genomic_benchmarks.data_check")
    result = dmb.list_genomic_benchmark_datasets(verbose=False)
    assert result["source"] == "package_index"
    assert result["datasets"] == sorted(gb.list_datasets())
    assert result["count"] == len(result["datasets"]) > 0


def test_list_genomic_benchmark_datasets_lists_local_directories(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks.data_check")
    (tmp_path / "beta").mkdir()
    (tmp_path / "alpha").mkdir()
    (tmp_path / "note.txt").write_text("not a dataset")
    result = dmb.list_genomic_benchmark_datasets(dataset_path=tmp_path, verbose=False)
    assert result["source"] == "local_directory"
    assert result["datasets"] == ["alpha", "beta"]


def test_describe_genomic_benchmark_counts_classes_and_splits(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks.data_check")
    folder = _dataset(tmp_path)
    result = dmb.describe_genomic_benchmark(folder, version=0, verbose=False)
    assert result["classes"] == {
        "neg": {"train": 2, "test": 2},
        "pos": {"train": 2, "test": 2},
    }
    assert result["total_sequences"] == 8
    assert (
        "Totally 8 sequences have been found, 4 for training and 4 for testing."
        in (result["messages"])
    )
    with pytest.raises(ConfigurationError, match="could not describe"):
        dmb.describe_genomic_benchmark(folder, version=1, verbose=False)


def test_check_genomic_benchmark_cache_reports_local_state(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks.data_check")
    cache = tmp_path / "cache"
    cache.mkdir()
    assert (
        dmb.check_genomic_benchmark_cache("demo_ds", cache_path=cache, verbose=False)["downloaded"]
        is False
    )
    (cache / "demo_ds").mkdir()
    assert (
        dmb.check_genomic_benchmark_cache("demo_ds", cache_path=cache, verbose=False)["downloaded"]
        is True
    )
    with pytest.raises(ConfigurationError, match="plain benchmark name"):
        dmb.check_genomic_benchmark_cache("../escape", verbose=False)


# evo2: fake-injected only (the real package cannot be imported in this environment)


_EVO2_CALLS: dict[str, dict] = {"init": {}, "generate": {}}


class _FakeEvo2:
    def __init__(self, model_name: str, local_path: str | None = None, use_kernels: bool = False):
        _EVO2_CALLS["init"] = {"model_name": model_name, "local_path": local_path}

    def generate(
        self,
        prompt_seqs,
        n_tokens=500,
        temperature=1.0,
        top_k=4,
        top_p=1.0,
        batched=True,
        cached_generation=True,
        verbose=1,
        force_prompt_threshold=None,
    ):
        _EVO2_CALLS["generate"] = {
            "prompt_seqs": list(prompt_seqs),
            "n_tokens": n_tokens,
            "temperature": temperature,
            "top_k": top_k,
            "top_p": top_p,
            "batched": batched,
            "verbose": verbose,
        }
        return ["ACGTACGT" for _ in prompt_seqs], [-1.25 for _ in prompt_seqs]


def _install_fake_evo2(monkeypatch: pytest.MonkeyPatch) -> None:
    module = types.ModuleType("evo2")
    module.Evo2 = _FakeEvo2
    monkeypatch.setitem(sys.modules, "evo2", module)


def test_generate_evo2_sequences_with_fake_module(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _install_fake_evo2(monkeypatch)
    weights = _write(tmp_path, "evo2_7b.pt", "placeholder")
    result = dmb.generate_evo2_sequences(
        ["ACGT", "GGCC"],
        weights_path=weights,
        n_tokens=8,
        temperature=0.7,
        top_k=2,
        top_p=0.9,
        verbose=False,
    )
    assert _EVO2_CALLS["init"] == {"model_name": "evo2_7b", "local_path": str(weights)}
    assert _EVO2_CALLS["generate"]["prompt_seqs"] == ["ACGT", "GGCC"]
    assert _EVO2_CALLS["generate"]["n_tokens"] == 8
    assert _EVO2_CALLS["generate"]["verbose"] == 0
    assert result["sequences"] == ["ACGTACGT", "ACGTACGT"]
    assert result["scores"] == [-1.25, -1.25]
    assert result["prompt_count"] == 2


def test_generate_evo2_sequences_maps_upstream_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    class _Failing(_FakeEvo2):
        def generate(self, *args, **kwargs):
            raise RuntimeError("CUDA out of memory")

    module = types.ModuleType("evo2")
    module.Evo2 = _Failing
    monkeypatch.setitem(sys.modules, "evo2", module)
    weights = _write(tmp_path, "evo2_7b.pt", "placeholder")
    with pytest.raises(BackendExecutionError, match="generation failed"):
        dmb.generate_evo2_sequences(["ACGT"], weights_path=weights, verbose=False)


def test_evo2_unavailable_in_this_environment_raises_install_hint(tmp_path: Path) -> None:
    """Runs only when the real evo2 package cannot be imported (no fakes installed here)."""

    try:
        import evo2  # noqa: F401
    except ImportError:
        pass
    else:
        pytest.skip("evo2 imports in this environment; the fake-module tests cover the wrapper")
    weights = _write(tmp_path, "evo2_7b.pt", "placeholder")
    with pytest.raises(BackendUnavailableError, match="neural-evo2"):
        dmb.generate_evo2_sequences(["ACGT"], weights_path=weights, verbose=False)


# kipoi: fake-injected only (kipoi reads git-backed sources and may clone or update them)


def _fake_kipoi(monkeypatch: pytest.MonkeyPatch) -> dict:
    pd = pytest.importorskip("pandas")
    calls: dict = {}
    module = types.ModuleType("kipoi")

    def list_models(sources=None):
        calls["sources"] = sources
        return pd.DataFrame(
            [
                {
                    "model": "m1",
                    "version": "0.1",
                    "type": "keras",
                    "license": "MIT",
                    "inputs": ["seq"],
                    "targets": ["y"],
                    "doc": "ignored",
                },
                {
                    "model": "m2",
                    "version": "0.2",
                    "type": "pytorch",
                    "license": None,
                    "inputs": ["seq"],
                    "targets": ["z"],
                    "doc": "ignored",
                },
            ]
        )

    def get_model_descr(model, source="kipoi"):
        calls["descr"] = (model, source)
        return SimpleNamespace(
            type="keras",
            info=SimpleNamespace(version="0.1", doc="demo", license="MIT", cite_as="doi:x"),
            schema=SimpleNamespace(
                inputs=[SimpleNamespace(name="seq")], targets=[SimpleNamespace(name="y")]
            ),
        )

    module.list_models = list_models
    module.get_model_descr = get_model_descr
    monkeypatch.setitem(sys.modules, "kipoi", module)
    return calls


def test_list_kipoi_models_with_fake_module(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _fake_kipoi(monkeypatch)
    result = dmb.list_kipoi_models(sources=["kipoi"], limit=1, verbose=False)
    assert calls["sources"] == ["kipoi"]
    assert result["total_models"] == 2 and result["returned"] == 1
    assert result["models"] == [
        {
            "model": "m1",
            "version": "0.1",
            "type": "keras",
            "license": "MIT",
            "inputs": ["seq"],
            "targets": ["y"],
        }
    ]
    with pytest.raises(ConfigurationError, match="sources must be"):
        dmb.list_kipoi_models(sources=[], verbose=False)


def test_describe_kipoi_model_with_fake_module(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _fake_kipoi(monkeypatch)
    result = dmb.describe_kipoi_model("deepsea", source="kipoi", verbose=False)
    assert calls["descr"] == ("deepsea", "kipoi")
    assert result["description"] == {
        "type": "keras",
        "version": "0.1",
        "doc": "demo",
        "license": "MIT",
        "cite_as": "doi:x",
        "inputs": ["seq"],
        "targets": ["y"],
    }


def test_describe_kipoi_model_accepts_single_array_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    # kipoi allows one array object in place of a list; the DeepSEA spec uses this form.
    _fake_kipoi(monkeypatch)
    module = sys.modules["kipoi"]
    descr = module.get_model_descr("DeepSEA/variantEffects")
    descr.schema = SimpleNamespace(
        inputs=SimpleNamespace(name="seq"), targets=SimpleNamespace(name="TFBS_DHS_probs")
    )
    module.get_model_descr = lambda model, source="kipoi": descr
    result = dmb.describe_kipoi_model("DeepSEA/variantEffects", verbose=False)
    assert result["description"]["inputs"] == ["seq"]
    assert result["description"]["targets"] == ["TFBS_DHS_probs"]


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="kipoi live test reads git-backed model sources; set DNAKIT_LIVE_TESTS=1",
)
def test_live_kipoi_model_description() -> None:
    pytest.importorskip("kipoi")
    result = dmb.describe_kipoi_model("DeepSEA/variantEffects", verbose=False)
    assert result["description"]["inputs"] == ["seq"]
