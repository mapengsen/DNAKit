"""Contract tests for the G11 epigenome wrappers: methylpy, MACS3 and deepTools."""

from __future__ import annotations

import inspect
import math
import random
import struct
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.epigenome import MACS3_functions as macs
from dnakit.epigenome import deeptools_functions as deep
from dnakit.epigenome import methylpy_functions as meth
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_CODES = {"T": 0, "C": 1, "A": 2, "G": 3}


def _write_twobit(path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a little-endian UCSC .2bit file; N is hard-masked, lowercase soft-masked."""

    def runs(flags: list[bool]) -> tuple[list[int], list[int]]:
        starts: list[int] = []
        sizes: list[int] = []
        index = 0
        while index < len(flags):
            if flags[index]:
                stop = index
                while stop < len(flags) and flags[stop]:
                    stop += 1
                starts.append(index)
                sizes.append(stop - index)
                index = stop
            else:
                index += 1
        return starts, sizes

    header = struct.pack("<IIII", 0x1A412743, 0, len(records), 0)
    index_bytes = b""
    body = b""
    offset = len(header) + sum(1 + len(name.encode()) + 4 for name, _ in records)
    for name, sequence in records:
        upper = sequence.upper()
        n_starts, n_sizes = runs([symbol == "N" for symbol in upper])
        mask_starts, mask_sizes = runs([symbol.islower() for symbol in sequence])
        packed = bytearray((len(sequence) + 3) // 4)
        for position, symbol in enumerate(upper):
            packed[position // 4] |= _CODES.get(symbol, 0) << (6 - 2 * (position % 4))
        record = struct.pack("<I", len(sequence))
        record += struct.pack("<I", len(n_starts)) + struct.pack(f"<{len(n_starts)}I", *n_starts)
        record += struct.pack(f"<{len(n_sizes)}I", *n_sizes)
        record += struct.pack("<I", len(mask_starts))
        record += struct.pack(f"<{len(mask_starts)}I", *mask_starts)
        record += struct.pack(f"<{len(mask_sizes)}I", *mask_sizes)
        record += struct.pack("<I", 0) + bytes(packed)
        raw = name.encode()
        index_bytes += bytes([len(raw)]) + raw + struct.pack("<I", offset)
        body += record
        offset += len(record)
    path.write_bytes(header + index_bytes + body)
    return path


def _write_bam(tmp_path: Path, name: str, sam_text: str) -> Path:
    """Sort, compress and index a SAM text as a BAM file."""

    pysam = pytest.importorskip("pysam")
    sam = tmp_path / f"{name}.sam"
    sam.write_text(sam_text, encoding="utf-8")
    unsorted = tmp_path / f"{name}.unsorted.bam"
    bam = tmp_path / f"{name}.bam"
    pysam.view("-b", "-o", str(unsorted), str(sam), catch_stdout=False)
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam


def _sam_header(length: int) -> str:
    return f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{length}\n"


def _signature_has_verbose(function: Any) -> bool:
    parameter = inspect.signature(function).parameters.get("verbose")
    return parameter is not None and parameter.kind is inspect.Parameter.KEYWORD_ONLY


# methylpy ------------------------------------------------------------------------------


def test_methylpy_all_public_functions_take_keyword_verbose() -> None:
    for module in (meth, macs, deep):
        for name in module.__all__:
            assert _signature_has_verbose(getattr(module, name)), name


def test_methylpy_expand_context_matches_upstream_and_is_sorted() -> None:
    utilities = pytest.importorskip("methylpy.utilities")
    result = meth.expand_methylation_context("CGN", verbose=False)
    assert result["contexts"] == sorted(utilities.expand_nucleotide_code("CGN"))
    assert "CGN" in result["contexts"]
    assert result["count"] == len(result["contexts"])
    with pytest.raises(ConfigurationError):
        meth.expand_methylation_context("CGZ", verbose=False)


def test_methylpy_encode_decode_round_trip() -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    encoded = meth.encode_read_c_positions("ACGTCCGA", verbose=False)
    assert encoded["encoding"] == call_mc.encode_c_positions("ACGTCCGA", False)
    decoded = meth.decode_read_c_positions("ATGTTTGA", encoded["encoding"], "+", verbose=False)
    assert decoded["sequence"] == "ACGTCCGA"
    with pytest.raises(ConfigurationError):
        meth.decode_read_c_positions("ACGT", "#", "x", verbose=False)


def test_methylpy_basecall_counts_hand_calculation() -> None:
    pytest.importorskip("methylpy.call_mc_se")
    # '.' is an unconverted C, 'T' a converted C, ',' a match on the reverse strand, 'A' a mismatch.
    result = meth.count_bisulfite_basecalls("C", ".,,TTA", verbose=False)
    assert result["consistent_calls"] == 3
    assert result["inconsistent_calls"] == 1
    assert result["unconverted_c"] == 1
    assert result["converted_c"] == 2
    tail = meth.count_bisulfite_basecalls("T", ",aA", verbose=False)
    assert (tail["consistent_calls"], tail["inconsistent_calls"]) == (1, 2)
    assert tail["unconverted_c"] is None


def test_methylpy_read_fasta_chromosome(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("methylpy.call_mc_se")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTNN\nACGT\n>chr2\nGGCC\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="FASTA index"):
        meth.read_fasta_chromosome(fasta, "chr1", verbose=False)
    pysam.faidx(str(fasta))
    result = meth.read_fasta_chromosome(fasta, "chr1", verbose=False)
    assert result["sequence"] == "ACGTACGTNNACGT"
    assert result["length"] == 14
    with pytest.raises(ConfigurationError, match="not listed"):
        meth.read_fasta_chromosome(fasta, "chrX", verbose=False)
    with pytest.raises(ConfigurationError, match="max_length"):
        meth.read_fasta_chromosome(fasta, "chr1", max_length=10, verbose=False)


def test_methylpy_convert_fastq_to_bisulfite(tmp_path: Path) -> None:
    pytest.importorskip("methylpy.call_mc_se")
    source = tmp_path / "reads.fq"
    source.write_text("@r1\nACGTCC\n+\nIIIIII\n@r2\nTTCG\n+\nIIII\n", encoding="utf-8")
    target = tmp_path / "converted.fq"
    result = meth.convert_fastq_to_bisulfite(source, target, verbose=False)
    assert result["reads"] == 2
    lines = target.read_text(encoding="utf-8").splitlines()
    # Every C becomes T: ACGTCC -> ATGTTT and TTCG -> TTTG.
    assert lines[1] == "ATGTTT"
    assert lines[5] == "TTTG"
    with pytest.raises(ConfigurationError, match="already exists"):
        meth.convert_fastq_to_bisulfite(source, target, verbose=False)


def test_methylpy_filter_bam_removes_methylated_reads(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("methylpy.call_mc_se")
    reference = "ACAACAACAACAACAACAAC" + "A" * 10
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr1\n{reference}\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    methylated = "ACAACAACAACAACAACAAC"
    unmethylated = "ATAATAATAATAATAATAAT"
    quality = "I" * 20
    sam = _sam_header(len(reference))
    sam += f"r1\t0\tchr1\t1\t60\t20M\t*\t0\t0\t{methylated}\t{quality}\n"
    sam += f"r2\t0\tchr1\t1\t60\t20M\t*\t0\t0\t{unmethylated}\t{quality}\n"
    bam = _write_bam(tmp_path, "mch", sam)
    output = tmp_path / "filtered.bam"
    result = meth.filter_bam_by_mch_level(bam, output, fasta, verbose=False)
    assert result["input_reads"] == 2
    assert result["kept_reads"] == 1
    assert result["removed_reads"] == 1
    with pysam.AlignmentFile(str(output), "rb") as handle:
        names = [read.query_name for read in handle]
    assert names == ["r2"]
    with pytest.raises(ConfigurationError, match=r"\.bam"):
        meth.filter_bam_by_mch_level(bam, tmp_path / "out.sam", fasta, verbose=False)


def test_methylpy_benjamini_hochberg_matches_upstream(tmp_path: Path) -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    allc = tmp_path / "allc.tsv"
    allc.write_text(
        "chr1\t10\t+\tCGN\t8\t10\nchr1\t20\t+\tCGN\t1\t10\nchr1\t30\t+\tCHH\t0\t9\n"
        "chr1\t40\t+\tCHH\t5\t9\n",
        encoding="utf-8",
    )
    binomial = tmp_path / "binomial.tsv"
    call_mc.allc_run_binom_tests(str(allc), str(binomial), 0.01, min_cov=1, sort_mem="64M")
    counts = {"CGN": 2, "CHH": 2}
    result = meth.benjamini_hochberg_methylation_cutoffs(
        [binomial], counts, sig_cutoff=0.5, verbose=False
    )
    expected = call_mc.benjamini_hochberg_correction_call_methylated_sites(
        [str(binomial)], counts, 0.5
    )
    assert result["cutoffs"] == pytest.approx(expected)
    unsorted = tmp_path / "unsorted.tsv"
    unsorted.write_text(
        "chr1\t10\t+\tCGN\t8\t10\t0.5\nchr1\t20\t+\tCGN\t1\t10\t0.1\n", encoding="utf-8"
    )
    with pytest.raises(ConfigurationError, match="sorted"):
        meth.benjamini_hochberg_methylation_cutoffs([unsorted], {"CGN": 2}, verbose=False)
    with pytest.raises(ConfigurationError, match="test count"):
        meth.benjamini_hochberg_methylation_cutoffs([binomial], {"CGN": 2}, verbose=False)


def test_methylpy_allc_index_offsets(tmp_path: Path) -> None:
    pytest.importorskip("methylpy.utilities")
    allc = tmp_path / "allc.tsv"
    chr1_lines = "chr1\t10\t+\tCGN\t8\t10\t1\nchr1\t20\t+\tCGN\t1\t10\t1\n"
    allc.write_text(chr1_lines + "chr2\t30\t+\tCHH\t0\t9\t1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="index_allc_table"):
        meth.read_allc_index(allc, verbose=False)
    built = meth.index_allc_table(allc, verbose=False)
    offsets = {item["chrom"]: item["offset"] for item in built["chromosomes"]}
    assert offsets == {"chr1": 0, "chr2": len(chr1_lines.encode())}
    assert meth.read_allc_index(allc, verbose=False)["chromosomes"] == built["chromosomes"]


def test_methylpy_sys_exit_becomes_configuration_error(tmp_path: Path, monkeypatch: Any) -> None:
    pysam = pytest.importorskip("pysam")
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    pysam.faidx(str(fasta))

    def refuse(*_args: Any, **_kwargs: Any) -> None:
        print("bad chromosome request")
        sys.exit(2)

    monkeypatch.setattr(call_mc, "get_chromosome_sequence", refuse)
    with pytest.raises(ConfigurationError, match="methylpy rejected"):
        meth.read_fasta_chromosome(fasta, "chr1", verbose=False)


def test_methylpy_missing_package_message(monkeypatch: Any, tmp_path: Path) -> None:
    for name in ("methylpy", "methylpy.call_mc_se", "methylpy.utilities"):
        monkeypatch.setitem(sys.modules, name, None)
    fastq = tmp_path / "r.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    with pytest.raises(BackendUnavailableError, match="bio-methylpy"):
        meth.convert_fastq_to_bisulfite(fastq, tmp_path / "out.fq", verbose=False)


# MACS3 ---------------------------------------------------------------------------------


def test_macs3_paired_fragment_pileup_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.PileupV2")
    result = macs.build_paired_fragment_pileup([10, 20, 25], [30, 40, 35], verbose=False)
    # Fragments [10,30), [20,40), [25,35): coverage 1, 2, 3, 2, 1 over the segments ending here.
    assert result["positions"] == [10, 20, 25, 30, 35, 40]
    assert result["values"] == [0.0, 1.0, 2.0, 3.0, 2.0, 1.0]
    with pytest.raises(ConfigurationError):
        macs.build_paired_fragment_pileup([30], [10], verbose=False)


def test_macs3_shifted_read_pileup_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.PileupV2")
    result = macs.build_shifted_read_pileup(
        [10, 12, 30],
        [40, 42],
        five_shift=0,
        three_shift=10,
        chrom_length=100,
        verbose=False,
    )
    # Plus reads cover [p, p+10); minus reads cover [p-10, p).
    assert result["positions"] == [10, 12, 20, 22, 30, 32, 40, 42]
    assert result["values"] == [0.0, 1.0, 2.0, 1.0, 0.0, 2.0, 3.0, 1.0]


def test_macs3_maxima_and_summit_filters_match_upstream() -> None:
    processing = pytest.importorskip("MACS3.Signal.SignalProcessing")
    numpy = pytest.importorskip("numpy")
    signal = [0, 1, 2, 3, 2, 1, 0, 1, 4, 6, 4, 1, 0, 0, 0]
    maxima = macs.find_signal_maxima(signal, window_size=3, verbose=False)["maxima"]
    expected = processing.maxima(numpy.asarray(signal, dtype=numpy.float32), 3)
    assert maxima == [int(item) for item in expected]
    assert maxima == [2, 3, 8, 9]
    valleys = macs.filter_summits_by_valleys(signal, [3, 9], verbose=False)
    assert valleys["summits"] == [3, 9]
    peaks = macs.filter_summits_by_peakiness(signal, [3, 9], verbose=False)
    assert peaks["summits"] == []
    with pytest.raises(ConfigurationError, match="odd"):
        macs.find_signal_maxima(signal, window_size=4, verbose=False)
    with pytest.raises(ConfigurationError, match="inside the signal"):
        macs.filter_summits_by_valleys(signal, [99], verbose=False)


def test_macs3_refine_peak_summit_matches_upstream() -> None:
    refine = pytest.importorskip("MACS3.Commands.refinepeak_cmd")
    numpy = pytest.importorskip("numpy")
    result = macs.refine_peak_summit(
        "chr1", [100, 110, 120], [105, 115], 90, 130, name="peak1", verbose=False
    )
    upstream = refine.find_summit(
        b"chr1",
        numpy.array([100, 110, 120], dtype=numpy.int32),
        numpy.array([105, 115], dtype=numpy.int32),
        90,
        130,
        b"peak1",
        100,
        5,
    )
    assert result["summit_start"] == upstream[1]
    assert result["summit_end"] == upstream[2]
    assert result["name"] == upstream[3].decode()
    assert result["score"] == pytest.approx(upstream[4])


def test_macs3_max_duplicate_tags_matches_upstream() -> None:
    callpeak = pytest.importorskip("MACS3.Commands.callpeak_cmd")
    result = macs.max_duplicate_tags(3_000_000_000, 10_000_000, verbose=False)
    assert result["max_duplicate_tags"] == int(
        callpeak.cal_max_dup_tags(3_000_000_000, 10_000_000, 1e-5)
    )
    assert result["max_duplicate_tags"] >= 1


def test_macs3_guess_tag_file_format(tmp_path: Path) -> None:
    pytest.importorskip("MACS3.IO.Parser")
    bed = tmp_path / "tags.bed"
    bed.write_text("chr1\t10\t30\tr1\t0\t+\nchr1\t20\t40\tr2\t0\t-\n", encoding="utf-8")
    result = macs.guess_tag_file_format(bed, verbose=False)
    assert result["parser"] == "BEDParser"
    assert result["format"] == "BED"
    broken = tmp_path / "broken.bin"
    broken.write_bytes(b"\x00\xff\x10")
    with pytest.raises(ConfigurationError, match="could not recognise"):
        macs.guess_tag_file_format(broken, verbose=False)


def test_macs3_genotype_quality_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.VariantStat")
    # GQ = -10 log10(exp(-20) / exp(-10)) = 100 / ln(10) = 43.43, truncated by MACS3 to 43.
    result = macs.calculate_genotype_quality(-10.0, -20.0, -30.0, verbose=False)
    assert result["genotype_quality"] == 43
    assert result["genotype_quality"] == pytest.approx(100 / math.log(10), abs=1.0)


def test_macs3_allele_model_matches_upstream() -> None:
    variant = pytest.importorskip("MACS3.Signal.VariantStat")
    numpy = pytest.importorskip("numpy")
    arrays = (
        [30, 25],
        [2, 3],
        [1, 1],
        [0, 2],
    )
    result = macs.calculate_allele_model_likelihood(*arrays, model="heterozygous", verbose=False)
    lnl, bic = variant.CalModel_Heter_noAS(
        *(numpy.array(item, dtype=numpy.int32) for item in arrays)
    )
    assert result["log_likelihood"] == pytest.approx(float(lnl))
    assert result["bic"] == pytest.approx(float(bic))
    homo = macs.calculate_allele_model_likelihood(*arrays, model="homozygous", verbose=False)
    assert homo["model"] == "homozygous"
    with pytest.raises(ConfigurationError):
        macs.calculate_allele_model_likelihood([1], [1], [1], [1, 2], verbose=False)


def test_macs3_merge_pileup_arrays_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.Pileup")
    first = {"positions": [10, 20], "values": [1.0, 2.0]}
    second = {"positions": [15, 20], "values": [3.0, 0.5]}
    result = macs.merge_pileup_arrays(first, second, reducer="max", verbose=False)
    # Max over the union of segment ends: 3 up to 10 and to 15, then 2 up to 20.
    assert result["positions"] == [10, 15, 20]
    assert result["values"] == [3.0, 3.0, 2.0]
    with pytest.raises(ConfigurationError):
        macs.merge_pileup_arrays(first, second, reducer="sum", verbose=False)


def test_macs3_hmmratac_weights_match_upstream() -> None:
    hmm = pytest.importorskip("MACS3.Signal.HMMR_Signal_Processing")
    means = [50.0, 150.0, 300.0, 450.0]
    stddevs = [20.0, 30.0, 40.0, 50.0]
    lengths = [50, 150, 300]
    result = macs.hmmratac_fragment_weights(lengths, means, stddevs, verbose=False)
    upstream = hmm.generate_weight_mapping(lengths, means, stddevs, 0.001)
    rows = {row["fragment_length"]: row for row in result["weights"]}
    for label, mapping in zip(result["classes"], upstream, strict=True):
        assert rows[50][label] == pytest.approx(mapping[50])
    assert rows[150]["mono"] > rows[150]["short"]
    with pytest.raises(ConfigurationError):
        macs.hmmratac_fragment_weights(lengths, means[:3], stddevs, verbose=False)


def test_macs3_missing_package_message(monkeypatch: Any) -> None:
    monkeypatch.setitem(sys.modules, "MACS3", None)
    monkeypatch.setitem(sys.modules, "MACS3.Signal.PileupV2", None)
    with pytest.raises(BackendUnavailableError, match="bio-macs3"):
        macs.build_shifted_read_pileup(
            [10],
            [20],
            five_shift=0,
            three_shift=10,
            chrom_length=100,
            verbose=False,
        )


# deepTools -----------------------------------------------------------------------------


def test_deeptools_gc_fraction_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    pytest.importorskip("deeptools.utilities")
    twobit = _write_twobit(tmp_path / "genome.2bit", [("chr1", "ACGTGGCCAANNGGTTACGT")])
    # ACGTGGCC has G and C at six of eight positions.
    fraction = deep.gc_fraction_from_twobit(twobit, "chr1", 0, 8, verbose=False)
    assert fraction["gc"] == pytest.approx(0.75)
    count = deep.gc_fraction_from_twobit(twobit, "chr1", 0, 8, fraction=False, verbose=False)
    assert count["gc"] == 6
    with pytest.raises(ConfigurationError, match="N"):
        deep.gc_fraction_from_twobit(twobit, "chr1", 0, 20, verbose=False)
    with pytest.raises(ConfigurationError, match="beyond"):
        deep.gc_fraction_from_twobit(twobit, "chr1", 0, 40, verbose=False)


def test_deeptools_parse_deeptools_region_matches_upstream_doctest() -> None:
    pytest.importorskip("deeptools.mapReduce")
    result = deep.parse_deeptools_region({"chr2": 1000}, "chr2:10:1001", verbose=False)
    assert (result["start"], result["end"], result["chunk_size"]) == (10, 1000, 990)
    assert result["regions"] == [{"chrom": "chr2", "length": 1000}]
    with pytest.raises(ConfigurationError, match="Unknown chromosome"):
        deep.parse_deeptools_region({"chr2": 1000}, "chr1:10:20", verbose=False)


def test_deeptools_scale_factors_identical_bams(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.estimateScaleFactor")
    random.seed(1)
    sequence = "".join(random.choice("ACGT") for _ in range(400))
    lines = []
    for index in range(60):
        position = 1 + (index * 6) % 300
        bases = sequence[position - 1 : position + 19]
        lines.append(f"r{index}\t0\tchr1\t{position}\t60\t20M\t*\t0\t0\t{bases}\t{'I' * 20}")
    sam = _sam_header(400) + "\n".join(lines) + "\n"
    first = _write_bam(tmp_path, "first", sam)
    second = _write_bam(tmp_path, "second", sam)
    result = deep.estimate_bam_scale_factors(
        [first, second],
        bin_length=50,
        number_of_samples=40,
        normalization_length=1000,
        verbose=False,
    )
    # Identical libraries get identical size factors of one.
    assert result["size_factors"] == pytest.approx([1.0, 1.0])
    assert result["sites_sampled"] > 0
    with pytest.raises(ConfigurationError, match="exactly two"):
        deep.estimate_bam_scale_factors([first], verbose=False)


def test_deeptools_bedgraph_to_bigwig_round_trip(tmp_path: Path) -> None:
    pyBigWig = pytest.importorskip("pyBigWig")
    pytest.importorskip("deeptools.writeBedGraph")
    bedgraph = tmp_path / "signal.bedgraph"
    bedgraph.write_text("chr1\t0\t10\t1.5\nchr1\t10\t20\t2\n", encoding="utf-8")
    target = tmp_path / "signal.bw"
    result = deep.bedgraph_to_bigwig({"chr1": 100}, [bedgraph], target, verbose=False)
    assert result["intervals"] == 2
    handle = pyBigWig.open(str(target))
    try:
        assert handle.chroms() == {"chr1": 100}
        assert list(handle.intervals("chr1")) == [(0, 10, 1.5), (10, 20, 2.0)]
    finally:
        handle.close()
    unsorted = tmp_path / "unsorted.bedgraph"
    unsorted.write_text("chr1\t10\t20\t1\nchr1\t0\t5\t1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="unsorted"):
        deep.bedgraph_to_bigwig({"chr1": 100}, [unsorted], tmp_path / "bad.bw", verbose=False)


def test_deeptools_chop_regions_hand_calculation() -> None:
    pytest.importorskip("deeptools.heatmapper")
    result = deep.chop_regions_for_profile(
        [(0, 200), (300, 400), (800, 900)], left=50, right=50, verbose=False
    )
    assert result["left_bins"] == [[0, 50]]
    assert result["body_bins"] == [[50, 200], [300, 400], [800, 850]]
    assert result["right_bins"] == [[850, 900]]
    with pytest.raises(ConfigurationError, match="non-overlapping"):
        deep.chop_regions_for_profile([(0, 200), (100, 300)], verbose=False)


def test_deeptools_chromosome_name_matching() -> None:
    pytest.importorskip("deeptools.utilities")
    assert deep.match_chromosome_name("chr1", ["1", "2"], verbose=False)["matched"] == "1"
    assert deep.match_chromosome_name("X", ["chr1"], verbose=False)["matched"] is None
    mapped = deep.map_twobit_to_bam_chromosomes(["chr1", "chr2"], ["1", "2"], verbose=False)
    assert mapped["mapping"] == {"chr1": "1", "chr2": "2"}
    assert mapped["complete"] is True
    partial = deep.map_twobit_to_bam_chromosomes(["chr1", "chrX"], ["chr1", "chr2"], verbose=False)
    assert partial["complete"] is False


def test_deeptools_fingerprint_jsd_matches_upstream() -> None:
    plot = pytest.importorskip("deeptools.plotFingerprint")
    chip = [0, 1, 2, 3, 10]
    control = [1, 1, 1, 1, 1]
    result = deep.fingerprint_jensen_shannon(chip, control, verbose=False)
    assert result["jsd"] == pytest.approx(float(plot.getJSDcommon(chip, control)))
    same = deep.fingerprint_jensen_shannon([1, 2, 3, 4], [1, 2, 3, 4], verbose=False)
    assert same["jsd"] == pytest.approx(0.0)


def test_deeptools_count_size_factors_match_documented_values() -> None:
    pytest.importorskip("deeptools.countReadsPerBin")
    matrix = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [0, 10, 0], [10, 5, 100]]
    result = deep.estimate_count_size_factors(matrix, verbose=False)
    # deepTools docstring values, returned as the inverse of the DESeq2 factors.
    assert result["size_factors"] == pytest.approx([1.305, 0.9932, 0.783], abs=1e-3)
    with pytest.raises(ConfigurationError):
        deep.estimate_count_size_factors([[0, 1], [0, 2]], verbose=False)


def test_deeptools_bigwig_chromosome_sizes(tmp_path: Path) -> None:
    pyBigWig = pytest.importorskip("pyBigWig")
    target = tmp_path / "track.bw"
    handle = pyBigWig.open(str(target), "w")
    handle.addHeader([("chr1", 100)])
    handle.addEntries(["chr1"], [0], ends=[10], values=[1.0])
    handle.close()
    result = deep.bigwig_chromosome_sizes([target], verbose=False)
    assert result["chromosomes"] == [{"chrom": "chr1", "length": 100}]
    assert result["non_matching"] == []


def test_deeptools_sys_exit_becomes_configuration_error(monkeypatch: Any) -> None:
    mapreduce = pytest.importorskip("deeptools.mapReduce")

    def refuse(*_args: Any, **_kwargs: Any) -> None:
        print("bad region")
        sys.exit(1)

    monkeypatch.setattr(mapreduce, "getUserRegion", refuse)
    with pytest.raises(ConfigurationError, match="deepTools rejected"):
        deep.parse_deeptools_region({"chr2": 1000}, "chr2:10:20", verbose=False)


def test_deeptools_missing_package_message(monkeypatch: Any) -> None:
    monkeypatch.setitem(sys.modules, "deeptools", None)
    monkeypatch.setitem(sys.modules, "deeptools.utilities", None)
    with pytest.raises(BackendUnavailableError, match="bio-deeptools"):
        deep.match_chromosome_name("chr1", ["1"], verbose=False)


def test_results_are_json_serializable(tmp_path: Path) -> None:
    import json

    pytest.importorskip("MACS3.Signal.PileupV2")
    pytest.importorskip("MACS3.Signal.VariantStat")
    pytest.importorskip("deeptools.heatmapper")
    pytest.importorskip("deeptools.countReadsPerBin")
    results = [
        macs.build_paired_fragment_pileup([10, 20], [30, 40], verbose=False),
        macs.calculate_allele_model_likelihood(
            [30, 25], [2, 3], [1, 1], [0, 2], model="homozygous", verbose=False
        ),
        macs.refine_peak_summit(
            "chr1", [100, 110, 120], [105, 115], 90, 130, name="peak1", verbose=False
        ),
        deep.chop_regions_for_profile([(0, 200), (300, 400)], left=50, verbose=False),
        deep.estimate_count_size_factors([[1, 2], [3, 4]], verbose=False),
        meth.expand_methylation_context("CHH", verbose=False),
    ]
    for result in results:
        json.dumps(result)
