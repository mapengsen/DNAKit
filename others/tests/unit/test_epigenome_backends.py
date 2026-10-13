"""Contract tests for optional cooler, cooltools, pairtools, MACS3 and deeptools adapters."""

from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest

from dnakit.epigenome.backends import (
    balance_contact_matrix,
    build_contact_matrix,
    build_extended_tag_pileup,
    build_fragment_pileup,
    calculate_coverage,
    calculate_expected_contacts,
    calculate_insulation_score,
    calculate_non_conversion_rate,
    calculate_site_binomial_pvalues,
    call_fragment_peaks,
    count_bam_reads_per_bin,
    deduplicate_pairs,
    estimate_fragment_lengths,
    fetch_contact_matrix,
    filter_allc_sites,
    read_pairs_header,
    summarize_bigwig_bins,
    summarize_pairs,
)
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

# Seven contacts on chr1 (1000 bp) and chr2 (600 bp); positions are 1-based.
# With binsize 100, chr1 owns bins 0-9 and chr2 owns bins 10-15.
TINY_PAIRS = """## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
r1\tchr1\t10\tchr1\t150\t+\t-
r2\tchr1\t20\tchr1\t160\t+\t-
r3\tchr1\t15\tchr1\t250\t+\t-
r4\tchr1\t510\tchr1\t520\t+\t-
r5\tchr1\t120\tchr2\t30\t+\t-
r6\tchr2\t40\tchr2\t90\t+\t-
r7\tchr2\t200\tchr2\t420\t-\t+
"""
TINY_CHROMSIZES = {"chr1": 1000, "chr2": 600}

# Pairs with a pair_type column, for pairtools statistics.
TYPED_PAIRS = """## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a\tchr1\t10\tchr1\t150\t+\t-\tUU
b\tchr1\t20\tchr1\t300\t+\t+\tUU
c\tchr1\t500\tchr1\t510\t-\t+\tUR
d\tchr1\t600\tchr2\t40\t+\t-\tUU
e\t!\t0\tchr1\t700\t+\t-\tNU
"""

# Six allc records: chromosome, position, strand, context, methylated reads, coverage, flag.
TINY_ALLC = """chrL\t1\t+\tCGN\t0\t10\t0
chrL\t2\t+\tCGG\t1\t10\t0
chrL\t3\t-\tCHH\t0\t4\t0
chr1\t5\t+\tCGN\t2\t2\t1
chr1\t9\t+\tCHG\t1\t2\t1
chr1\t12\t+\tCGN\t3\t10\t1
"""

# a-b-c form one duplicate cluster (each neighbour within 3 bp on both sides); d and e
# share positions but have opposite strands, so they are not duplicates.
DUPLICATE_PAIRS = """## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a\tchr1\t10\tchr1\t150\t+\t-\tUU
b\tchr1\t11\tchr1\t151\t+\t-\tUU
c\tchr1\t14\tchr1\t152\t+\t-\tUU
d\tchr1\t300\tchr1\t400\t+\t-\tUU
e\tchr1\t300\tchr1\t400\t-\t+\tUU
f\tchr1\t500\tchr1\t600\t+\t-\tUU
g\t!\t0\tchr1\t700\t+\t-\tNU
"""


def _hide_package(monkeypatch: pytest.MonkeyPatch, name: str) -> None:
    """Make a package and every already-imported submodule unimportable for one test."""

    for module in list(sys.modules):
        if module == name or module.startswith(name + "."):
            monkeypatch.setitem(sys.modules, module, None)
    monkeypatch.setitem(sys.modules, name, None)


def _write(tmp_path: Path, text: str, name: str = "input.txt") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _big_pairs(tmp_path: Path) -> Path:
    """4000 intra-chromosomal contacts on a 4 kb chromosome with a boundary near 2 kb."""

    rng = random.Random(7)
    lines = ["## pairs format v1.0", "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"]
    for index in range(4000):
        first = rng.randrange(1, 4001)
        same_side = first <= 2000
        while True:
            distance = int(rng.expovariate(1 / 300)) + 1
            second = first + distance if rng.random() < 0.5 else first - distance
            if 1 <= second <= 4000 and ((second <= 2000) == same_side or rng.random() < 0.1):
                break
        low, high = sorted((first, second))
        lines.append(f"r{index}\tchr1\t{low}\tchr1\t{high}\t+\t-")
    return _write(tmp_path, "\n".join(lines) + "\n", "big.pairs")


def _tiny_cool(tmp_path: Path, name: str = "tiny.cool") -> Path:
    pytest.importorskip("cooler")
    pairs = _write(tmp_path, TINY_PAIRS, "tiny.pairs")
    result = build_contact_matrix(
        pairs, tmp_path / name, chromsizes=TINY_CHROMSIZES, binsize=100, verbose=False
    )
    return Path(result["cool_path"])


def _bam(tmp_path: Path, name: str, records: list[tuple[str, int, int, int, int]]) -> Path:
    """Write a coordinate-sorted, indexed BAM from (name, flag, start, mapq, template_length)."""

    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    path = tmp_path / name
    ordered = sorted(records, key=lambda record: record[2])
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for query, flag, start, mapq, tlen in ordered:
            read = pysam.AlignedSegment()
            read.query_name = query
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigar = [(0, 50)]
            read.next_reference_id = 0
            read.template_length = tlen
            out.write(read)
    pysam.index(str(path))
    return path


def _pair_bam(tmp_path: Path) -> Path:
    """Five read pairs: 50 bp reads, template length 250, mate 2 starts 200 bp after mate 1."""

    records: list[tuple[str, int, int, int, int]] = []
    for index, start in enumerate([100, 110, 120, 400, 410]):
        records.append((f"q{index}", 99, start, 60, 250))
        records.append((f"q{index}", 147, start + 200, 60, -250))
    return _bam(tmp_path, "pairs.bam", records)


def _bigwig(tmp_path: Path, name: str, entries: list[tuple[int, int, float]]) -> Path:
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / name
    handle = pybigwig.open(str(path), "w")
    handle.addHeader([("chr1", 1000)])
    for start, end, value in entries:
        handle.addEntries(["chr1"], [start], ends=[end], values=[value])
    handle.close()
    return path


# cooler ----------------------------------------------------------------------


def test_build_contact_matrix_bins_pairs_by_hand_counts(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    result = build_contact_matrix(
        _write(tmp_path, TINY_PAIRS, "tiny.pairs"),
        tmp_path / "tiny.cool",
        chromsizes=TINY_CHROMSIZES,
        binsize=100,
        verbose=False,
    )

    # Hand count: (0,1)=2 from r1 and r2; (0,2)=1; (1,10)=1; (5,5)=1; (10,10)=1; (11,14)=1.
    assert result["nonzero_pixels"] == 6
    assert result["total_contacts"] == 7
    assert result["bins"] == 16
    pixels = cooler.Cooler(result["cool_path"]).pixels()[:]
    observed = {
        (int(a), int(b)): int(c)
        for a, b, c in zip(pixels["bin1_id"], pixels["bin2_id"], pixels["count"], strict=True)
    }
    assert observed == {(0, 1): 2, (0, 2): 1, (1, 10): 1, (5, 5): 1, (10, 10): 1, (11, 14): 1}


def test_build_contact_matrix_rejects_existing_output_bad_bins_and_empty_input(
    tmp_path: Path,
) -> None:
    pytest.importorskip("cooler")
    pairs = _write(tmp_path, TINY_PAIRS, "tiny.pairs")
    existing = tmp_path / "exists.cool"
    existing.write_text("occupied", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        build_contact_matrix(pairs, existing, chromsizes=TINY_CHROMSIZES, binsize=100)
    with pytest.raises(ConfigurationError, match="max_bins"):
        build_contact_matrix(
            pairs, tmp_path / "a.cool", chromsizes=TINY_CHROMSIZES, binsize=10, max_bins=20
        )
    with pytest.raises(ConfigurationError, match="binsize"):
        build_contact_matrix(pairs, tmp_path / "b.cool", chromsizes=TINY_CHROMSIZES, binsize=True)
    empty = _write(
        tmp_path, "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n", "e.pairs"
    )
    with pytest.raises(ConfigurationError, match="cooler rejected"):
        build_contact_matrix(empty, tmp_path / "c.cool", chromsizes=TINY_CHROMSIZES, binsize=100)
    assert not (tmp_path / "c.cool").exists()


def test_fetch_contact_matrix_matches_direct_cooler_fetch(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    np = pytest.importorskip("numpy")
    cool = _tiny_cool(tmp_path)

    result = fetch_contact_matrix(cool, "chr1", verbose=False)

    expected = cooler.Cooler(str(cool)).matrix(balance=False).fetch("chr1")
    assert np.array_equal(np.array(result["matrix"]), expected)
    assert result["matrix"][0][:3] == [0.0, 2.0, 1.0]  # hand count row for bin 0
    assert result["bins"][1] == {"chrom": "chr1", "start": 100, "end": 200}
    with pytest.raises(ConfigurationError, match="max_region_bins"):
        fetch_contact_matrix(cool, "chr1", max_region_bins=5, verbose=False)
    with pytest.raises(ConfigurationError, match="no balancing weights"):
        fetch_contact_matrix(cool, "chr1", balance=True, verbose=False)


def test_balance_contact_matrix_stores_weights_from_direct_call(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    big = build_contact_matrix(
        _big_pairs(tmp_path),
        tmp_path / "big.cool",
        chromsizes={"chr1": 4000},
        binsize=100,
        verbose=False,
    )
    cool = Path(big["cool_path"])

    result = balance_contact_matrix(cool, min_nnz=1, verbose=False)

    direct_weights, direct_stats = cooler.balance_cooler(
        cooler.Cooler(str(cool)), store=False, min_nnz=1
    )
    stored = cooler.Cooler(str(cool)).bins()[:]["weight"].to_numpy()
    assert result["bins"] == 40
    assert result["finite_weights"] == 40
    assert result["converged"] is True
    assert result["scale"] == pytest.approx(float(direct_stats["scale"]))
    assert list(stored) == pytest.approx(list(direct_weights))


def test_balance_contact_matrix_refuses_empty_filter_and_keeps_file(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cool = _tiny_cool(tmp_path)

    with pytest.raises(ConfigurationError, match="No bin passed"):
        balance_contact_matrix(cool, min_nnz=10, verbose=False)
    assert "weight" not in cooler.Cooler(str(cool)).bins()[:1].columns


def test_calculate_coverage_matches_hand_sums_and_direct_call(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    cool = _tiny_cool(tmp_path)

    result = calculate_coverage(cool, ignore_diags=0, verbose=False)

    rows = result["bins"]
    cis = {index: row["cis"] for index, row in enumerate(rows)}
    total = {index: row["total"] for index, row in enumerate(rows)}
    # Diagonal pixels count twice, and trans pixels only enter the total.
    assert cis[0] == 3.0 and total[0] == 3.0  # (0,1)=2 and (0,2)=1
    assert cis[1] == 2.0 and total[1] == 3.0  # (0,1)=2 is cis; (1,10)=1 is trans
    assert cis[5] == 2.0 and total[5] == 2.0  # (5,5)=1 on both sides
    assert cis[10] == 2.0 and total[10] == 3.0  # (10,10)=1 twice; (1,10)=1 in total only
    assert cis[11] == 1.0 and total[11] == 1.0
    assert cis[3] == 0.0 and total[3] == 0.0
    cis_direct, total_direct = cooltools.coverage(
        cooler.Cooler(str(cool)), ignore_diags=0, clr_weight_name=None
    )
    assert [row["cis"] for row in rows] == [float(value) for value in cis_direct]
    assert [row["total"] for row in rows] == [float(value) for value in total_direct]


def test_calculate_expected_contacts_matches_hand_means(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    cool = _tiny_cool(tmp_path)

    result = calculate_expected_contacts(cool, ignore_diags=0, verbose=False)

    by_key = {(row["chrom"], row["distance_bins"]): row for row in result["rows"]}
    # Distance 0 on chr1: 10 valid pixels, one contact (5,5).
    assert by_key[("chr1", 0)]["valid_pixels"] == 10
    assert by_key[("chr1", 0)]["sum"] == 1.0
    assert by_key[("chr1", 0)]["mean"] == pytest.approx(0.1)
    # Distance 1 on chr1: 9 valid pixels, contacts (0,1)=2 only.
    assert by_key[("chr1", 1)]["sum"] == 2.0
    assert by_key[("chr1", 1)]["mean"] == pytest.approx(2 / 9)
    # Distance 3 on chr2: 3 valid pixels, contact (11,14)=1.
    assert by_key[("chr2", 3)]["valid_pixels"] == 3
    assert by_key[("chr2", 3)]["mean"] == pytest.approx(1 / 3)
    direct = cooltools.expected_cis(
        cooler.Cooler(str(cool)),
        intra_only=True,
        smooth=False,
        aggregate_smoothed=False,
        clr_weight_name=None,
        ignore_diags=0,
    )
    assert len(result["rows"]) == len(direct)


def test_calculate_insulation_score_finds_synthetic_boundary(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    np = pytest.importorskip("numpy")
    big = build_contact_matrix(
        _big_pairs(tmp_path),
        tmp_path / "big.cool",
        chromsizes={"chr1": 4000},
        binsize=100,
        verbose=False,
    )
    cool = Path(big["cool_path"])
    balance_contact_matrix(cool, min_nnz=1, verbose=False)

    result = calculate_insulation_score(cool, 500, verbose=False)

    boundaries = [row["start"] for row in result["bins"] if row["is_boundary"]]
    assert any(1800 <= start <= 2100 for start in boundaries)
    direct = cooltools.insulation(cooler.Cooler(str(cool)), [500], verbose=False)
    observed = np.array(
        [
            np.nan if row["log2_insulation_score"] is None else row["log2_insulation_score"]
            for row in result["bins"]
        ],
        dtype=float,
    )
    assert np.allclose(observed, direct["log2_insulation_score_500"].to_numpy(), equal_nan=True)
    with pytest.raises(ConfigurationError, match="no balancing weights"):
        calculate_insulation_score(_tiny_cool(tmp_path, "unbalanced.cool"), 500, verbose=False)


def test_missing_cooler_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pairs = _write(tmp_path, TINY_PAIRS, "tiny.pairs")
    _hide_package(monkeypatch, "cooler")
    with pytest.raises(BackendUnavailableError, match="bio-cooler"):
        build_contact_matrix(pairs, tmp_path / "x.cool", chromsizes=TINY_CHROMSIZES, binsize=100)


def test_missing_cooltools_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("cooler")
    cool = tmp_path / "placeholder.cool"
    cool.write_text("not read before cooltools is imported", encoding="utf-8")
    _hide_package(monkeypatch, "cooltools")
    with pytest.raises(BackendUnavailableError, match="bio-cooltools"):
        calculate_coverage(cool, verbose=False)


# pairtools -------------------------------------------------------------------


def test_read_pairs_header_reports_columns_and_chromsizes(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    result = read_pairs_header(_write(tmp_path, DUPLICATE_PAIRS, "dup.pairs"), verbose=False)

    assert result["header_lines"] == 4
    assert result["columns"][:4] == ["readID", "chrom1", "pos1", "chrom2"]
    assert result["columns"][-1] == "pair_type"
    assert result["chromsizes"] == {"chr1": 1000}
    assert result["sorted"] == "chr1-pos1-chr2-pos2"


def test_summarize_pairs_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    result = summarize_pairs(_write(tmp_path, TYPED_PAIRS, "typed.pairs"), verbose=False)

    assert result["total"] == 5
    assert result["total_mapped"] == 4  # a, b, c, d
    assert result["total_single_sided_mapped"] == 1  # e
    assert result["total_unmapped"] == 0
    assert result["cis"] == 3  # a, b, c
    assert result["trans"] == 1  # d
    assert result["pair_types"] == {"UU": 3, "UR": 1, "NU": 1}


def test_summarize_pairs_needs_pair_type_column(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    untyped = _write(tmp_path, TINY_PAIRS, "untyped.pairs")
    with pytest.raises(ConfigurationError, match="pair_type"):
        summarize_pairs(untyped, verbose=False)


def test_deduplicate_pairs_counts_chained_duplicates_and_writes_outputs(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    source = _write(tmp_path, DUPLICATE_PAIRS, "dup.pairs")
    unique = tmp_path / "unique.pairs"
    unmapped = tmp_path / "unmapped.pairs"

    result = deduplicate_pairs(
        source, output_path=unique, unmapped_path=unmapped, max_mismatch=3, verbose=False
    )

    assert result["total"] == 7
    assert result["total_mapped"] == 6
    assert result["total_single_sided_mapped"] == 1
    assert result["total_dups"] == 2  # b and c join a's cluster
    assert result["total_nodups"] == 4  # a, d, e, f
    assert result["duplicate_fraction"] == pytest.approx(2 / 6)
    body = [
        line for line in unique.read_text(encoding="utf-8").splitlines() if not line.startswith("#")
    ]
    assert [line.split("\t")[0] for line in body] == ["a", "d", "e", "f"]
    assert "g\t!" in unmapped.read_text(encoding="utf-8")


def test_deduplicate_pairs_validates_method_and_output(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    source = _write(tmp_path, DUPLICATE_PAIRS, "dup.pairs")
    with pytest.raises(ConfigurationError, match="method"):
        deduplicate_pairs(source, method="mean", verbose=False)
    existing = _write(tmp_path, "occupied", "taken.pairs")
    with pytest.raises(ConfigurationError, match="already exists"):
        deduplicate_pairs(source, output_path=existing, verbose=False)


def test_missing_pairtools_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pairs = _write(tmp_path, DUPLICATE_PAIRS, "dup.pairs")
    _hide_package(monkeypatch, "pairtools")
    with pytest.raises(BackendUnavailableError, match="bio-pairtools"):
        read_pairs_header(pairs, verbose=False)


# MACS3 -----------------------------------------------------------------------


def test_build_fragment_pileup_matches_hand_segments() -> None:
    pytest.importorskip("MACS3")
    result = build_fragment_pileup([10, 50, 100], [40, 90, 140], verbose=False)

    assert [(row["start"], row["end"], row["value"]) for row in result["segments"]] == [
        (0, 10, 0.0),
        (10, 40, 1.0),
        (40, 50, 0.0),
        (50, 90, 1.0),
        (90, 100, 0.0),
        (100, 140, 1.0),
    ]
    scaled = build_fragment_pileup([10], [40], scale_factor=2.0, baseline_value=0.5, verbose=False)
    assert [row["value"] for row in scaled["segments"]] == [0.5, 2.0]


def test_build_extended_tag_pileup_matches_hand_segments() -> None:
    pytest.importorskip("MACS3")
    result = build_extended_tag_pileup([100, 200, 205], 50, verbose=False)

    # Tags cover [50,150), [150,250) and [155,255).
    assert [(row["start"], row["end"], row["value"]) for row in result["segments"]] == [
        (0, 50, 0.0),
        (50, 155, 1.0),
        (155, 250, 2.0),
        (250, 255, 1.0),
    ]


def test_call_fragment_peaks_reports_summits_above_threshold() -> None:
    pytest.importorskip("MACS3")
    starts = [100] * 5 + [400] * 2 + [1000] * 4
    ends = [300] * 5 + [500] * 2 + [1200] * 4

    strong = call_fragment_peaks(starts, ends, min_value=4.5, verbose=False)
    assert strong["peaks"] == [{"summit": 200, "height": 5.0}]

    default = call_fragment_peaks(starts, ends, min_value=1.0, verbose=False)
    assert default["peaks"] == [{"summit": 200, "height": 5.0}, {"summit": 1100, "height": 4.0}]

    short = call_fragment_peaks(starts, ends, min_value=1.5, min_length=100, verbose=False)
    assert short["peaks"] == [
        {"summit": 200, "height": 5.0},
        {"summit": 450, "height": 2.0},
        {"summit": 1100, "height": 4.0},
    ]


def test_fragment_inputs_are_validated() -> None:
    pytest.importorskip("MACS3")
    with pytest.raises(ConfigurationError, match="greater than its start"):
        build_fragment_pileup([50], [50], verbose=False)
    with pytest.raises(ConfigurationError, match="same length"):
        build_fragment_pileup([10, 20], [30], verbose=False)
    with pytest.raises(ConfigurationError, match="max_fragments"):
        build_fragment_pileup([10], [30], max_fragments=0, verbose=False)
    with pytest.raises(ConfigurationError, match="tag position"):
        build_extended_tag_pileup([True], 10, verbose=False)


def test_missing_macs3_package_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    _hide_package(monkeypatch, "MACS3")
    with pytest.raises(BackendUnavailableError, match="bio-macs3"):
        build_fragment_pileup([10], [40], verbose=False)


# methylpy --------------------------------------------------------------------


def test_filter_allc_sites_keeps_matching_context_chrom_and_coverage(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    allc = _write(tmp_path, TINY_ALLC, "tiny.allc")

    result = filter_allc_sites(
        allc,
        tmp_path / "filtered.allc",
        chromosomes=["chrL", "chr1"],
        mc_type="CGN",
        min_coverage=2,
        verbose=False,
    )

    # CGN expands to CGA/CGC/CGG/CGT/CGN, so CGG counts. CHH and CHG are dropped by context.
    kept = (tmp_path / "filtered.allc").read_text(encoding="utf-8").splitlines()
    assert [line.split("\t")[1] for line in kept] == ["1", "2", "5", "12"]
    assert result["input_sites"] == 6
    assert result["kept_sites"] == 4


def test_calculate_non_conversion_rate_matches_hand_ratio(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    allc = _write(tmp_path, TINY_ALLC, "tiny.allc")

    result = calculate_non_conversion_rate(allc, "chrL", verbose=False)
    single = calculate_non_conversion_rate(allc, "chrL:2-2", verbose=False)

    # chrL: methylated reads 0 + 1 + 0 = 1 out of coverage 10 + 10 + 4 = 24.
    assert result["non_conversion_rate"] == pytest.approx(1 / 24)
    assert single["non_conversion_rate"] == pytest.approx(0.1)
    assert allc.read_text(encoding="utf-8") == TINY_ALLC
    assert not list(tmp_path.glob("*.idx"))


def test_calculate_site_binomial_pvalues_matches_hand_and_scipy(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    stats = pytest.importorskip("scipy.stats")
    allc = _write(tmp_path, TINY_ALLC, "tiny.allc")

    result = calculate_site_binomial_pvalues(
        allc,
        tmp_path / "pvalues.tsv",
        non_conversion=0.01,
        min_coverage=2,
        max_report=3,
        verbose=False,
    )

    # Hand values: 2 of 2 reads gives 0.01**2; 1 of 2 reads gives 1 - 0.99**2.
    # The chr1:12 site is checked against scipy's binomial survival function.
    expected = [
        ("chr1", 5, 0.01**2),
        ("chr1", 12, float(stats.binom.sf(2, 10, 0.01))),
        ("chr1", 9, 1 - 0.99**2),
    ]
    lowest = [(row["chrom"], row["pos"], row["p_value"]) for row in result["lowest_p_values"]]
    assert [item[:2] for item in lowest] == [item[:2] for item in expected]
    for observed, wanted in zip(lowest, expected, strict=True):
        assert observed[2] == pytest.approx(wanted[2], rel=1e-9)
    assert result["tested_sites"] == 4  # chrL:2 (1/10), chr1:5, chr1:9 and chr1:12
    assert result["sites_written"] == 6
    assert result["mc_class_counts"] == {"CGG": 1, "CGN": 2, "CHG": 1}


def test_methylpy_inputs_are_validated(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    allc = _write(tmp_path, TINY_ALLC, "tiny.allc")
    with pytest.raises(ConfigurationError, match="control_region"):
        calculate_non_conversion_rate(allc, "chr L:1", verbose=False)
    with pytest.raises(ConfigurationError, match="not in the output allc file"):
        calculate_non_conversion_rate(allc, "chrZZ", verbose=False)
    spaced = tmp_path / "with space.tsv"
    with pytest.raises(ConfigurationError, match="may only contain"):
        calculate_site_binomial_pvalues(allc, spaced, non_conversion=0.01, verbose=False)
    with pytest.raises(ConfigurationError, match="sort_buffer"):
        calculate_site_binomial_pvalues(
            allc, tmp_path / "ok.tsv", non_conversion=0.01, sort_buffer="1G;rm", verbose=False
        )
    with pytest.raises(ConfigurationError, match="non_conversion"):
        calculate_site_binomial_pvalues(
            allc, tmp_path / "bad.tsv", non_conversion=1.5, verbose=False
        )


def test_missing_methylpy_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    allc = _write(tmp_path, TINY_ALLC, "tiny.allc")
    _hide_package(monkeypatch, "methylpy")
    with pytest.raises(BackendUnavailableError, match="bio-methylpy"):
        filter_allc_sites(allc, tmp_path / "x.allc", verbose=False)


# deeptools -------------------------------------------------------------------


def test_summarize_bigwig_bins_matches_hand_means(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    first = _bigwig(tmp_path, "a.bw", [(0, 500, 2.0), (500, 1000, 4.0)])
    second = _bigwig(tmp_path, "b.bw", [(0, 1000, 1.0)])

    result = summarize_bigwig_bins([first, second], bin_size=200, verbose=False)

    assert [row["values"] for row in result["bins"]] == [
        [2.0, 1.0],
        [2.0, 1.0],
        [3.0, 1.0],
        [4.0, 1.0],
        [4.0, 1.0],
    ]
    assert result["bins"][2] == {"chrom": "chr1", "start": 400, "end": 600, "values": [3.0, 1.0]}


def test_summarize_bigwig_bins_ignores_missing_bases_and_reports_null(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    gapped = _bigwig(tmp_path, "gap.bw", [(0, 100, 2.0), (300, 400, 4.0)])

    result = summarize_bigwig_bins([gapped], bin_size=200, verbose=False)

    assert [row["values"] for row in result["bins"]] == [[2.0], [4.0], [None], [None], [None]]


def test_count_bam_reads_per_bin_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    pairs = _pair_bam(tmp_path)

    result = count_bam_reads_per_bin([pairs], bin_size=100, verbose=False)

    counts = [row["counts"][0] for row in result["bins"]]
    # Starts 100-120 and 300-320 give three reads each, 400-410 two, 600-610 two.
    assert counts == [0, 3, 0, 3, 2, 0, 2, 0, 0, 0]
    assert result["bins"][1] == {"chrom": "chr1", "start": 100, "end": 200, "counts": [3]}


def test_count_bam_reads_per_bin_applies_mapping_quality(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    # The read at 170 has MAPQ 0 and is dropped by min_mapping_quality=10.
    reads = [
        ("s1", 0, 80, 60, 0),
        ("s2", 0, 150, 60, 0),
        ("dup", 1024, 160, 60, 0),
        ("lowq", 0, 170, 0, 0),
    ]
    bam = _bam(tmp_path, "single.bam", reads)

    kept = count_bam_reads_per_bin([bam], bin_size=100, min_mapping_quality=10, verbose=False)
    unfiltered = count_bam_reads_per_bin([bam], bin_size=100, verbose=False)

    assert [row["counts"][0] for row in kept["bins"]][:3] == [1, 3, 1]
    assert [row["counts"][0] for row in unfiltered["bins"]][:3] == [1, 4, 2]


def test_count_bam_reads_per_bin_rejects_unindexed_bam(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    indexed = _pair_bam(tmp_path)
    unindexed = tmp_path / "unindexed.bam"
    unindexed.write_bytes(indexed.read_bytes())
    with pytest.raises(ConfigurationError, match="index"):
        count_bam_reads_per_bin([unindexed], bin_size=100, verbose=False)


def test_estimate_fragment_lengths_reports_template_length(tmp_path: Path) -> None:
    pytest.importorskip("deeptools")
    pairs = _pair_bam(tmp_path)

    result = estimate_fragment_lengths(pairs, verbose=False)

    assert result["fragment_length"]["sample_size"] == 5
    assert result["fragment_length"]["median"] == 250.0
    assert result["fragment_length"]["std"] == 0.0
    assert result["read_length"]["median"] == 50.0


def test_missing_deeptools_package_names_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = _bigwig(tmp_path, "a.bw", [(0, 1000, 1.0)])
    _hide_package(monkeypatch, "deeptools")
    with pytest.raises(BackendUnavailableError, match="bio-deeptools"):
        summarize_bigwig_bins([bigwig], bin_size=200, verbose=False)


def test_public_functions_reject_bad_arguments(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="file path"):
        read_pairs_header(12, verbose=False)
    with pytest.raises(ConfigurationError, match="bigwig_files"):
        summarize_bigwig_bins([], verbose=False)
    placeholder = _write(tmp_path, "not a bam", "placeholder.bam")
    with pytest.raises(ConfigurationError, match="number_of_processors"):
        estimate_fragment_lengths(placeholder, number_of_processors=0, verbose=False)
