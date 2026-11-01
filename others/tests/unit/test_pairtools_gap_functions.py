"""Tests for dnakit.epigenome.pairtools_gap_functions (pairtools command-line facades)."""

from __future__ import annotations

import inspect
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("pairtools")  # the mark-split helper imports the pairtools package
from dnakit.epigenome import pairtools_gap_functions as facade
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_LOCAL_PAIRTOOLS = Path(sys.executable).parent / "pairtools"
_PAIRTOOLS = shutil.which("pairtools") or (str(_LOCAL_PAIRTOOLS) if _LOCAL_PAIRTOOLS.is_file() else None)
needs_pairtools = pytest.mark.skipif(_PAIRTOOLS is None, reason="pairtools executable is not installed")

_COLUMNS = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
_PAIRS = (
    "## pairs format v1.0.0\n"
    "#sorted: chr-pos\n"
    "#shape: upper triangle\n"
    "#genome_assembly: unknown\n"
    "#chromsize: chr1 1000\n"
    "#chromsize: chr2 1000\n"
    f"#columns: {_COLUMNS}\n"
    "r1\tchr1\t100\tchr1\t500\t+\t-\tUU\n"
    "r2\tchr1\t100\tchr1\t500\t+\t-\tUU\n"
    "r3\tchr1\t200\tchr2\t300\t+\t+\tUU\n"
    "r4\tchr1\t600\tchr2\t700\t+\t+\tUU\n"
)
_CHROMS = "chr1\t1000\nchr2\t1000\n"
_FRAGS = "chr1\t0\t300\nchr1\t300\t1000\nchr2\t0\t1000\n"
_SAM = (
    "@HD\tVN:1.0\tSO:queryname\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "@SQ\tSN:chr2\tLN:1000\n"
    "r1\t99\tchr1\t100\t60\t10M\t=\t400\t300\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
    "r1\t147\tchr1\t400\t60\t10M\t=\t100\t-300\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
)


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


def _data_rows(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]


def test_public_functions_take_keyword_only_verbose() -> None:
    assert facade.__all__
    for name in facade.__all__:
        verbose = inspect.signature(getattr(facade, name)).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True


def test_mark_split_pair_as_dup_sets_pair_type_and_keeps_newline() -> None:
    record = ["r1", "chr1", "100", "chr1", "500", "+", "-", "UU\n"]
    marked = facade.pairtools_mark_split_pair_as_dup(record)
    assert marked[7] == "DD\n"  # upstream appends the original newline to the last column
    assert marked[-1] == "DD\n"
    assert record[7] == "UU\n"  # the caller's list is not modified
    with pytest.raises(ConfigurationError):
        facade.pairtools_mark_split_pair_as_dup([])


def test_missing_executable_raises_with_install_hint(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    with pytest.raises(BackendUnavailableError) as info:
        facade.pairtools_sort(str(pairs), str(tmp_path / "out.pairs"), pairtools_path=str(tmp_path / "missing"))
    assert "conda install -c bioconda pairtools" in (info.value.hint or "")


def test_invalid_arguments_are_rejected(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    with pytest.raises(ConfigurationError):
        facade.pairtools_sample(str(pairs), 1.5, str(tmp_path / "s.pairs"), pairtools_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.pairtools_select("", str(pairs), str(tmp_path / "x.pairs"), pairtools_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.pairtools_stats([str(pairs), str(pairs)], str(tmp_path / "st.txt"), pairtools_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.pairtools_phase(str(pairs), str(tmp_path / "p.pairs"), phase_suffixes=["bad suffix"],
                               pairtools_path="/nope")


@needs_pairtools
def test_dedup_marks_identical_pairs_as_duplicates(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    out, dups = tmp_path / "out.pairs", tmp_path / "dups.pairs"
    result = facade.pairtools_dedup(str(pairs), str(out), output_dups=str(dups), pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed"
    assert len(_data_rows(dups)) == 1
    assert len(_data_rows(out)) == 3


@needs_pairtools
def test_filterbycov_and_sample(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    assert facade.pairtools_filterbycov(str(pairs), str(tmp_path / "f.pairs"), max_cov=2,
                                        pairtools_path=_PAIRTOOLS)["status"] == "completed"
    sampled = tmp_path / "s.pairs"
    assert facade.pairtools_sample(str(pairs), 1.0, str(sampled), seed=1,
                                   pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert len(_data_rows(sampled)) == 4


@needs_pairtools
def test_flip_sort_and_select(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    chroms = _write(tmp_path, "chroms.tsv", _CHROMS)
    flipped = tmp_path / "flip.pairs"
    assert facade.pairtools_flip(str(pairs), str(chroms), str(flipped), pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert len(_data_rows(flipped)) == 4
    sorted_out = tmp_path / "sorted.pairs"
    assert facade.pairtools_sort(str(pairs), str(sorted_out), pairtools_path=_PAIRTOOLS)["status"] == "completed"
    selected = tmp_path / "sel.pairs"
    result = facade.pairtools_select("chrom1 == 'chr1' and chrom2 == 'chr2'", str(pairs), str(selected),
                                     pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed"
    assert [row.split("\t")[0] for row in _data_rows(selected)] == ["r3", "r4"]


@needs_pairtools
def test_header_commands(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    chroms = _write(tmp_path, "chroms.tsv", _CHROMS)
    generated = tmp_path / "gen.pairs"
    assert facade.pairtools_header_generate(str(pairs), str(generated), chroms_path=str(chroms),
                                            pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert "#chromsize: chr1 1000" in generated.read_text(encoding="utf-8")
    transferred = tmp_path / "tr.pairs"
    assert facade.pairtools_header_transfer(str(pairs), str(pairs), str(transferred),
                                            pairtools_path=_PAIRTOOLS)["status"] == "completed"
    set_out = tmp_path / "cols.pairs"
    assert facade.pairtools_header_set_columns(str(pairs), _COLUMNS, str(set_out),
                                               pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert f"#columns: {_COLUMNS}" in set_out.read_text(encoding="utf-8")
    validated = tmp_path / "val.pairs"
    assert facade.pairtools_header_validate_columns(str(pairs), str(pairs), str(validated),
                                                    pairtools_path=_PAIRTOOLS)["status"] == "completed"


@needs_pairtools
def test_markasdup_and_merge(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    marked = tmp_path / "marked.pairs"
    assert facade.pairtools_markasdup(str(pairs), str(marked), pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert len(_data_rows(marked)) == 4
    merged = tmp_path / "merged.pairs"
    result = facade.pairtools_merge([str(pairs), str(pairs)], str(merged), pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed"
    assert len(_data_rows(merged)) == 8


@needs_pairtools
def test_parse_and_parse2_read_sam_alignments(tmp_path: Path) -> None:
    sam = _write(tmp_path, "in.sam", _SAM)
    chroms = _write(tmp_path, "chroms.tsv", _CHROMS)
    for function, name in ((facade.pairtools_parse, "parsed.pairs"), (facade.pairtools_parse2, "parsed2.pairs")):
        out = tmp_path / name
        result = function(str(sam), str(chroms), str(out), pairtools_path=_PAIRTOOLS)
        assert result["status"] == "completed", result["log_tail"]
        assert out.read_text(encoding="utf-8").startswith("## pairs format")


@needs_pairtools
def test_restrict_scaling_split_and_stats(tmp_path: Path) -> None:
    pairs = _write(tmp_path, "in.pairs", _PAIRS)
    frags = _write(tmp_path, "frags.bed", _FRAGS)
    restricted = tmp_path / "res.pairs"
    assert facade.pairtools_restrict(str(pairs), str(frags), str(restricted),
                                     pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert facade.pairtools_scaling(str(pairs), str(tmp_path / "scale.txt"),
                                    pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert facade.pairtools_split(str(pairs), str(tmp_path / "sp.pairs"), str(tmp_path / "sp.sam"),
                                  pairtools_path=_PAIRTOOLS)["status"] == "completed"
    stats_a = tmp_path / "a.stats"
    assert facade.pairtools_stats([str(pairs)], str(stats_a), pairtools_path=_PAIRTOOLS)["status"] == "completed"
    assert "total" in stats_a.read_text(encoding="utf-8")
    stats_merged = tmp_path / "merged.stats"
    result = facade.pairtools_stats([str(stats_a), str(stats_a)], str(stats_merged), merge=True,
                                    pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed"


_PHASE_CHROMS = "chr1-1\t1000\nchr1-2\t1000\nchr2-1\t1000\nchr2-2\t1000\n"
_PHASE_SAM = (
    "@HD\tVN:1.0\tSO:queryname\n"
    "@SQ\tSN:chr1-1\tLN:1000\n"
    "@SQ\tSN:chr1-2\tLN:1000\n"
    "@SQ\tSN:chr2-1\tLN:1000\n"
    "@SQ\tSN:chr2-2\tLN:1000\n"
    "r1\t99\tchr1-1\t100\t60\t10M\t=\t400\t300\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tXB:Z:A\tAS:i:10\tXS:i:0\n"
    "r1\t147\tchr1-1\t400\t60\t10M\t=\t100\t-300\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tXB:Z:A\tAS:i:10\tXS:i:0\n"
)


@needs_pairtools
def test_phase_after_parse_with_add_columns(tmp_path: Path) -> None:
    sam = _write(tmp_path, "phased.sam", _PHASE_SAM)
    chroms = _write(tmp_path, "phase_chroms.tsv", _PHASE_CHROMS)
    parsed = tmp_path / "parsed.pairs"
    result = facade.pairtools_parse(str(sam), str(chroms), str(parsed), add_columns="XB,AS,XS", min_mapq=0,
                                    pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed", result["log_tail"]
    phased = tmp_path / "phased.pairs"
    result = facade.pairtools_phase(str(parsed), str(phased), phase_suffixes=["-1", "-2"], pairtools_path=_PAIRTOOLS)
    assert result["status"] == "completed", result["log_tail"]
    assert "phase1" in phased.read_text(encoding="utf-8")
