"""Tests for dnakit.epigenome.cooler_gap2_functions (cooler command-line facades, headless plots).

Offline: the real ``cooler`` program runs on tiny coolers built in tmp_path. Tests that need an
external program that is not installed (bgzip, tabix, pairix, pigz, pypairix) are skipped with the
reason.
"""

from __future__ import annotations

import gzip
import importlib
import importlib.util
import inspect
import os
import shutil
import stat
import subprocess
import sys
import textwrap
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
pytest.importorskip("h5py")
cooler = pytest.importorskip("cooler")

from dnakit.epigenome import cooler_gap2_functions as gap2  # noqa: E402
from dnakit.epigenome.cooler_gap2_functions import (  # noqa: E402
    cooler_attrs,
    cooler_balance,
    cooler_balance_check,
    cooler_cload_hiclib,
    cooler_cload_pairix,
    cooler_cload_pairs,
    cooler_cload_tabix,
    cooler_coarsen,
    cooler_cp,
    cooler_csort,
    cooler_csort_flip_pairs,
    cooler_csort_index_pairs,
    cooler_csort_read_pairs,
    cooler_csort_sort_pairs,
    cooler_digest,
    cooler_dump,
    cooler_info,
    cooler_ln,
    cooler_load,
    cooler_ls,
    cooler_makebins,
    cooler_merge,
    cooler_mv,
    cooler_pairs_header,
    cooler_show,
    cooler_show_view,
    cooler_tree,
    cooler_zoomify,
)
from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402

DATA = (
    Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)


def _cooler_program() -> str | None:
    try:
        return gap2._find_program("cooler", None, "")
    except BackendUnavailableError:
        return None


requires_cooler = pytest.mark.skipif(
    _cooler_program() is None, reason="the cooler command-line program was not found"
)

SMALL_SIZES = {"chr1": 1000, "chr2": 600}
BIG_SIZES = {"chr1": 100_000, "chr2": 60_000}
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


# ---------------------------------------- data helpers


def _offsets(sizes: dict[str, int], binsize: int) -> dict[str, int]:
    offsets, total = {}, 0
    for name, length in sizes.items():
        offsets[name] = total
        total += -(-length // binsize)
    return offsets


def _random_pairs(sizes: dict[str, int], count: int, seed: int) -> list[tuple[str, int, str, int]]:
    rng = np.random.default_rng(seed)
    names = list(sizes)
    pairs = []
    for _ in range(count):
        c1, c2 = names[int(rng.integers(len(names)))], names[int(rng.integers(len(names)))]
        pairs.append(
            (c1, int(rng.integers(1, sizes[c1] + 1)), c2, int(rng.integers(1, sizes[c2] + 1)))
        )
    return pairs


def _expected_pixels(
    pairs: list[tuple[str, int, str, int]], sizes: dict[str, int], binsize: int
) -> Counter:
    offsets = _offsets(sizes, binsize)
    counter: Counter = Counter()
    for c1, p1, c2, p2 in pairs:
        b1, b2 = offsets[c1] + (p1 - 1) // binsize, offsets[c2] + (p2 - 1) // binsize
        counter[(min(b1, b2), max(b1, b2))] += 1
    return counter


def _write_sizes(path: Path, sizes: dict[str, int]) -> Path:
    path.write_text("".join(f"{name}\t{length}\n" for name, length in sizes.items()))
    return path


def _write_pairs(path: Path, pairs: list[tuple[str, int, str, int]]) -> Path:
    lines = [
        "## pairs format v1.0",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 score",
    ]
    lines += [
        f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\t1.5" for i, (c1, p1, c2, p2) in enumerate(pairs)
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def _make_cool(path: Path, sizes: dict[str, int], binsize: int, counter: Counter) -> Path:
    bins = cooler.util.binnify(pd.Series(sizes), binsize)
    pixels = pd.DataFrame(
        [(b1, b2, n) for (b1, b2), n in sorted(counter.items())],
        columns=["bin1_id", "bin2_id", "count"],
    )
    cooler.create_cooler(str(path), bins, pixels)
    return path


def _dataset(
    folder: Path, name: str, sizes: dict[str, int], binsize: int, count: int, seed: int
) -> SimpleNamespace:
    pairs = _random_pairs(sizes, count, seed)
    expected = _expected_pixels(pairs, sizes, binsize)
    return SimpleNamespace(
        sizes=sizes,
        binsize=binsize,
        pairs=pairs,
        expected=expected,
        sizes_file=_write_sizes(folder / f"{name}.sizes", sizes),
        pairs_file=_write_pairs(folder / f"{name}.pairs", pairs),
        cool=_make_cool(folder / f"{name}.cool", sizes, binsize, expected),
        bins_spec=f"{folder / f'{name}.sizes'}:{binsize}",
    )


@pytest.fixture(scope="module")
def workspace(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("cooler_gap2")


@pytest.fixture(scope="module")
def small(workspace: Path) -> SimpleNamespace:
    return _dataset(workspace, "small", SMALL_SIZES, 100, 600, 1)


@pytest.fixture(scope="module")
def big(workspace: Path) -> SimpleNamespace:
    return _dataset(workspace, "big", BIG_SIZES, 100, 5000, 2)


def _sum_and_nnz(path: Path | str) -> tuple[int, int]:
    clr = cooler.Cooler(str(path))
    return int(clr.info["sum"]), int(clr.info["nnz"])


def _pixel_counter(path: Path | str) -> Counter:
    pixels = cooler.Cooler(str(path)).pixels()[:]
    return Counter(
        {
            (int(a), int(b)): int(n)
            for a, b, n in zip(pixels.bin1_id, pixels.bin2_id, pixels["count"], strict=True)
        }
    )


def _copy(source: Path, destination: Path) -> Path:
    shutil.copyfile(source, destination)
    return destination


def _tsv(path: Path, **kwargs: Any) -> Any:
    return pd.read_csv(path, sep="\t", **kwargs)


# ---------------------------------------- contract


def test_every_public_function_follows_the_conventions() -> None:
    assert len(gap2.__all__) == 28
    for name in gap2.__all__:
        function = getattr(gap2, name)
        signature = inspect.signature(function)
        last = list(signature.parameters.values())[-1]
        assert (
            last.name == "verbose"
            and last.kind is inspect.Parameter.KEYWORD_ONLY
            and last.default is True
        ), name
        assert (function.__doc__ or "").strip(), name
        assert signature.return_annotation in (dict[str, Any], "dict[str, Any]"), name
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, (name, parameter.name)


def test_implemented_and_skipped_rows_exist_in_the_installed_cooler() -> None:
    qualnames = [
        "cooler.cli.balance.cli",
        "cooler.cli.cload.cli",
        "cooler.cli.cload.parse_bins",
        "cooler.cli.cload.parse_field_param",
        "cooler.cli.cload.parse_kv_list_param",
        "cooler.cli.cload.get_header",
        "cooler.cli.cload.register_subcommand",
        "cooler.cli.cload.add_arg_help",
        "cooler.cli.coarsen.cli",
        "cooler.cli.csort.make_read_command",
        "cooler.cli.csort.make_flip_command",
        "cooler.cli.csort.make_sort_command",
        "cooler.cli.csort.make_index_command",
        "cooler.cli.digest.cli",
        "cooler.cli.digest.exit_on_broken_pipe",
        "cooler.cli.dump.make_annotator",
        "cooler.cli.show.interactive",
        "cooler.cli.zoomify.legacy_zoomify",
        "cooler.cli.zoomify.preferred_sequence",
        "cooler.cli.zoomify.invoke_balance",
        "cooler.fileops.tree_widget",
        "cooler.fileops.pprint_attr_tree",
        "cooler.fileops.pprint_data_tree",
        "cooler.cli.cli",
    ]
    for qualname in qualnames:
        parts = qualname.split(".")
        for cut in range(len(parts), 0, -1):
            try:
                target: Any = importlib.import_module(".".join(parts[:cut]))
            except ImportError:
                continue
            for attribute in parts[cut:]:
                target = getattr(target, attribute)
            break
        else:  # pragma: no cover
            pytest.fail(f"cannot resolve {qualname}")
        assert target is not None


def test_group_rows_are_the_shared_root_group() -> None:
    root = importlib.import_module("cooler.cli").cli
    for module in (
        "balance",
        "cload",
        "coarsen",
        "csort",
        "digest",
        "dump",
        "fileops",
        "info",
        "load",
        "makebins",
        "merge",
        "show",
        "zoomify",
    ):
        assert importlib.import_module(f"cooler.cli.{module}").cli is root


# ---------------------------------------- program lookup


def test_program_lookup_explicit_path_then_path_then_interpreter_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    explicit = tmp_path / "my-cooler"
    explicit.write_text("#!/bin/sh\n")
    assert gap2._find_program("cooler", str(explicit), "hint") == str(explicit)
    with pytest.raises(BackendUnavailableError) as missing:
        gap2._find_program("cooler", str(tmp_path / "nope"), gap2._COOLER_HINT)
    assert "pip install cooler" in (missing.value.hint or "")
    folder = tmp_path / "bin"
    folder.mkdir()
    monkeypatch.setattr(sys, "executable", str(folder / "python"))
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    with pytest.raises(BackendUnavailableError) as absent:
        gap2._find_program("cooler", None, gap2._COOLER_HINT)
    assert "conda install" in (absent.value.hint or "")
    script = folder / "cooler"
    script.write_text("#!/bin/sh\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    assert gap2._find_program("cooler", None, "hint") == str(script)


def test_missing_program_raises_backend_unavailable_with_install_hint(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        cooler_info(str(small.cool), cooler_path=str(tmp_path / "no-such-cooler"))
    assert "cooler" in (info.value.hint or "").lower()


def test_timeout_must_be_positive(small: SimpleNamespace) -> None:
    with pytest.raises(ConfigurationError):
        gap2._timeout(0)
    with pytest.raises(ConfigurationError):
        gap2._timeout(True)


# ---------------------------------------- info, bins, digest


@requires_cooler
def test_info_returns_parsed_info_field_and_metadata(small: SimpleNamespace) -> None:
    result = cooler_info(str(small.cool))
    assert result["status"] == "completed" and result["backend_version"]
    assert result["info"]["bin-size"] == 100 and result["info"]["nbins"] == 16
    assert (result["info"]["sum"], result["info"]["nnz"]) == (600, len(small.expected))
    assert cooler_info(str(small.cool), field="nchroms")["value"] == 2
    assert cooler_info(str(small.cool), metadata=True)["metadata"] == {}
    missing = cooler_info(str(small.cool), field="no-such-field")
    assert missing["status"] == "failed" and "not found" in missing["log_tail"]
    with pytest.raises(ConfigurationError):
        cooler_info(str(small.cool), field="nbins", metadata=True)


@requires_cooler
def test_makebins_matches_binnify(small: SimpleNamespace, tmp_path: Path) -> None:
    result = cooler_makebins(
        str(small.sizes_file), 300, str(tmp_path / "bins.bed"), header=True, rel_ids=1
    )
    assert result["status"] == "completed", result["log_tail"]
    table = _tsv(Path(result["outputs"][0]))
    expected = cooler.util.binnify(pd.Series(SMALL_SIZES), 300)
    assert list(table.columns) == ["chrom", "start", "end", "id"]
    assert table[["chrom", "start", "end"]].values.tolist() == expected.values.tolist()
    assert table["id"].tolist() == [1, 2, 3, 4, 1, 2]
    assert result["line_count"] == 7 and result["preview"][0] == "chrom\tstart\tend\tid"
    with pytest.raises(ConfigurationError):
        cooler_makebins(str(small.sizes_file), 300, str(tmp_path / "bins.bed"))  # exists now
    with pytest.raises(ConfigurationError):
        cooler_makebins(str(small.sizes_file), 0, str(tmp_path / "other.bed"))
    with pytest.raises(ConfigurationError):
        cooler_makebins(str(small.sizes_file), 300, str(tmp_path / "other.bed"), rel_ids=2)  # type: ignore[arg-type]


@requires_cooler
def test_digest_matches_the_python_api_for_a_file_and_a_folder(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    pytest.importorskip("pyfaidx")
    genome = tmp_path / "genes" / "genome.fa"
    genome.parent.mkdir()
    shutil.copyfile(DATA / "genome.fa", genome)  # pyfaidx writes a .fai next to the FASTA
    sizes = _write_sizes(tmp_path / "genome.sizes", {"chr1": 120000, "chr2": 80000})
    result = cooler_digest(
        str(sizes), str(genome), "HindIII", str(tmp_path / "frags.bed"), header=True, rel_ids=0
    )
    assert result["status"] == "completed", result["log_tail"]
    table = _tsv(Path(result["outputs"][0]))
    records = cooler.util.load_fasta(["chr1", "chr2"], str(genome))
    expected = cooler.util.digest(records, "HindIII")
    assert len(table) == len(expected) and table["end"].iloc[-1] == 80000
    assert (
        table[["chrom", "start", "end"]].values.tolist()
        == expected[["chrom", "start", "end"]].values.tolist()
    )
    assert table.groupby("chrom")["id"].first().tolist() == [0, 0]
    folder = cooler_digest(str(sizes), str(genome.parent), "DpnII", str(tmp_path / "dpn.bed"))
    assert folder["status"] == "completed" and folder["line_count"] > result["line_count"]
    with pytest.raises(ConfigurationError):
        cooler_digest(str(sizes), str(genome), "Hind III", str(tmp_path / "bad.bed"))
    with pytest.raises(ConfigurationError):
        cooler_digest(
            str(sizes), str(tmp_path / "missing.fa"), "HindIII", str(tmp_path / "bad.bed")
        )


# ---------------------------------------- dump (make_annotator)


@requires_cooler
def test_dump_pixels_bins_chroms_and_regions(small: SimpleNamespace, tmp_path: Path) -> None:
    result = cooler_dump(str(small.cool), str(tmp_path / "px.tsv"), header=True)
    assert result["status"] == "completed", result["log_tail"]
    pixels = _tsv(Path(result["outputs"][0]))
    assert list(pixels.columns) == ["bin1_id", "bin2_id", "count"]
    assert pixels["count"].sum() == 600 and len(pixels) == len(small.expected)
    assert result["line_count"] == len(pixels) + 1
    bins = cooler_dump(
        str(small.cool),
        str(tmp_path / "bins.tsv.gz"),
        table="bins",
        columns=["chrom", "start"],
        header=True,
    )
    with gzip.open(bins["outputs"][0], "rt") as handle:
        assert handle.readline().strip() == "chrom\tstart" and len(handle.readlines()) == 16
    chroms = _tsv(
        Path(
            cooler_dump(str(small.cool), str(tmp_path / "chroms.tsv"), table="chroms", header=True)[
                "outputs"
            ][0]
        )
    )
    assert chroms.values.tolist() == [["chr1", 1000], ["chr2", 600]]
    region = _tsv(
        Path(
            cooler_dump(
                str(small.cool), str(tmp_path / "region.tsv"), region="chr1:0-500", header=True
            )["outputs"][0]
        )
    )
    assert region[["bin1_id", "bin2_id"]].max().max() < 5 and len(region) < len(pixels)
    lower = cooler_dump(
        str(small.cool),
        str(tmp_path / "lower.tsv"),
        region="chr1:0-300",
        fill_lower=True,
        header=True,
    )
    assert len(_tsv(Path(lower["outputs"][0]))) >= len(
        _tsv(
            Path(
                cooler_dump(
                    str(small.cool), str(tmp_path / "upper.tsv"), region="chr1:0-300", header=True
                )["outputs"][0]
            )
        )
    )
    with pytest.raises(ConfigurationError):
        cooler_dump(str(small.cool), str(tmp_path / "bad.tsv"), region2="chr1:0-300")
    with pytest.raises(ConfigurationError):
        cooler_dump(str(small.cool), str(tmp_path / "bad.tsv"), float_format="%g; rm")


@requires_cooler
def test_dump_join_balanced_annotate_and_one_based_options(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    balanced_cool = _copy(small.cool, tmp_path / "balanced.cool")
    assert (
        cooler_balance(str(balanced_cool), min_nnz=1, mad_max=0, nproc=1)["weights_stored"] is True
    )
    result = cooler_dump(
        str(balanced_cool),
        str(tmp_path / "joined.tsv"),
        header=True,
        join=True,
        balanced=True,
        annotate=["weight"],
        one_based_starts=True,
        one_based_ids=True,
    )
    assert result["status"] == "completed", result["log_tail"]
    table = _tsv(Path(result["outputs"][0]))
    assert list(table.columns)[:7] == [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "count",
    ]
    assert {"balanced", "weight1", "weight2"} <= set(table.columns)
    clr = cooler.Cooler(str(balanced_cool))
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]
    assert (
        table["start1"].tolist()
        == (bins["start"].to_numpy()[pixels["bin1_id"].to_numpy()] + 1).tolist()
    )
    expected = (
        bins["weight"].to_numpy()[pixels["bin1_id"]]
        * bins["weight"].to_numpy()[pixels["bin2_id"]]
        * pixels["count"]
    )
    both = np.isfinite(expected.to_numpy()) & table["balanced"].notna().to_numpy()
    assert both.any() and np.allclose(
        table["balanced"].to_numpy()[both], expected.to_numpy()[both], rtol=1e-4
    )
    ids = cooler_dump(
        str(balanced_cool),
        str(tmp_path / "ids.tsv"),
        header=True,
        annotate=["chrom"],
        one_based_ids=True,
    )
    assert _tsv(Path(ids["outputs"][0]))["bin1_id"].min() == 1
    plain = cooler_dump(str(small.cool), str(tmp_path / "nobalance.tsv"), balanced=True)
    assert plain["status"] == "failed" and "weights not found" in plain["log_tail"]


# ---------------------------------------- balance


@requires_cooler
def test_balance_writes_weights_and_needs_force_to_repeat(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    target = _copy(small.cool, tmp_path / "balance.cool")
    assert cooler_balance_check(str(target))["balanced"] is False
    result = cooler_balance(str(target), min_nnz=1, mad_max=0, nproc=1, tol=1e-6, max_iters=300)
    assert result["status"] == "completed" and result["weights_stored"] is True, result["log_tail"]
    assert cooler_balance_check(str(target))["balanced"] is True
    weights = cooler.Cooler(str(target)).bins()[:]["weight"].to_numpy()
    assert np.isfinite(weights).any()
    again = cooler_balance(str(target), min_nnz=1, mad_max=0, nproc=1)
    assert again["status"] == "failed" and "already exists" in again["log_tail"]
    forced = cooler_balance(str(target), min_nnz=1, mad_max=0, nproc=1, force=True)
    assert forced["status"] == "completed"
    named = cooler_balance(str(target), min_nnz=1, mad_max=0, nproc=1, name="other", cis_only=True)
    assert (
        named["weights_stored"] is True
        and cooler_balance_check(str(target), name="other")["balanced"] is True
    )
    with pytest.raises(ConfigurationError):
        cooler_balance(str(target), cis_only=True, trans_only=True)
    bad = cooler_balance_check(str(tmp_path / "balance.cool::/nope"))
    assert bad["balanced"] is None and bad["status"] == "failed"


@requires_cooler
def test_balance_blacklist_masks_bins_and_discard_policy(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    target = _copy(small.cool, tmp_path / "masked.cool")
    blacklist = tmp_path / "bad.bed"
    blacklist.write_text(
        "chr1\t0\t100\nchr2\t300\t400\n"
    )  # cooler reads a one-line file as a header and fails
    result = cooler_balance(
        str(target), min_nnz=1, mad_max=0, nproc=1, blacklist=str(blacklist), ignore_diags=1
    )
    assert result["status"] == "completed", result["log_tail"]
    weights = cooler.Cooler(str(target)).bins()[:]["weight"].to_numpy()
    assert np.isnan(weights[0]) and np.isnan(weights[13]) and np.isfinite(weights[1:13]).any()
    never = _copy(small.cool, tmp_path / "never.cool")
    discarded = cooler_balance(
        str(never),
        min_nnz=1,
        mad_max=0,
        nproc=1,
        max_iters=1,
        tol=1e-12,
        convergence_policy="discard",
    )
    assert discarded["status"] == "completed" and discarded["weights_stored"] is False


# ---------------------------------------- cload


@requires_cooler
def test_cload_pairs_bins_a_text_file(small: SimpleNamespace, tmp_path: Path) -> None:
    result = cooler_cload_pairs(
        small.bins_spec,
        str(small.pairs_file),
        str(tmp_path / "pairs.cool"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        fields=["score=8:dtype=float"],
        assembly="toy1",
        chunksize=200,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["subcommand"] == "cload pairs" and result["command"][:3] == [
        "cooler",
        "cload",
        "pairs",
    ]
    clr = cooler.Cooler(result["outputs"][0])
    assert _pixel_counter(result["outputs"][0]) == small.expected
    assert clr.info["genome-assembly"] == "toy1"
    assert clr.pixels()[:]["score"].sum() == pytest.approx(600 * 1.5)
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            small.bins_spec,
            str(small.pairs_file),
            result["outputs"][0],
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            "no-such.sizes:100",
            str(small.pairs_file),
            str(tmp_path / "x.cool"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            small.bins_spec,
            str(small.pairs_file),
            str(tmp_path / "x.cool"),
            chrom1=0,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            small.bins_spec,
            str(small.pairs_file),
            str(tmp_path / "x.cool"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            fields=["score; rm"],
        )


@requires_cooler
def test_cload_pairs_accepts_bed_bins_zero_based_and_appends_to_a_file(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    bed = Path(
        cooler_makebins(str(small.sizes_file), 100, str(tmp_path / "bins.bed"))["outputs"][0]
    )
    shifted = tmp_path / "zero.pairs"
    shifted.write_text(
        "".join(
            f"r{i}\t{c1}\t{p1 - 1}\t{c2}\t{p2 - 1}\n"
            for i, (c1, p1, c2, p2) in enumerate(small.pairs)
        )
    )
    target = tmp_path / "multi.cool"
    first = cooler_cload_pairs(
        str(bed),
        str(shifted),
        f"{target}::/first",
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        zero_based=True,
    )
    assert first["status"] == "completed", first["log_tail"]
    second = cooler_cload_pairs(
        str(bed),
        str(shifted),
        f"{target}::/second",
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        zero_based=True,
        append=True,
    )
    assert second["status"] == "completed", second["log_tail"]
    assert cooler_ls(str(target))["coolers"] == [f"{target}::/first", f"{target}::/second"]
    assert _pixel_counter(f"{target}::/first") == small.expected
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            str(bed), str(shifted), f"{target}::/third", chrom1=2, pos1=3, chrom2=4, pos2=5
        )  # file exists


def _sorted_flipped_pairs(
    small: SimpleNamespace,
) -> list[tuple[int, int, int, int, tuple[str, int, str, int]]]:
    rank = {name: index for index, name in enumerate(small.sizes)}
    rows = []
    for c1, p1, c2, p2 in small.pairs:
        if (rank[c1], p1) > (rank[c2], p2):
            c1, p1, c2, p2 = c2, p2, c1, p1
        rows.append((rank[c1], p1, rank[c2], p2, (c1, p1, c2, p2)))
    return sorted(rows)


def test_cload_tabix_bins_an_indexed_file(small: SimpleNamespace, tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("dill")
    if _cooler_program() is None:
        pytest.skip("the cooler command-line program was not found")
    plain = tmp_path / "sorted.pairs"
    plain.write_text(
        "".join(
            f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\n"
            for i, (*_, (c1, p1, c2, p2)) in enumerate(_sorted_flipped_pairs(small))
        )
    )
    pysam.tabix_compress(str(plain), str(tmp_path / "sorted.pairs.gz"))
    pysam.tabix_index(
        str(tmp_path / "sorted.pairs.gz"), seq_col=1, start_col=2, end_col=2, force=True
    )
    result = cooler_cload_tabix(
        small.bins_spec,
        str(tmp_path / "sorted.pairs.gz"),
        str(tmp_path / "tabix.cool"),
        chrom2=4,
        pos2=5,
        nproc=1,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _pixel_counter(result["outputs"][0]) == small.expected
    with pytest.raises(ConfigurationError):
        cooler_cload_tabix(small.bins_spec, str(small.pairs_file), str(tmp_path / "noindex.cool"))


@requires_cooler
def test_cload_hiclib_bins_a_frag_file(small: SimpleNamespace, tmp_path: Path) -> None:
    h5py = pytest.importorskip("h5py")
    rows = _sorted_flipped_pairs(small)
    with h5py.File(tmp_path / "frags.h5", "w") as handle:
        handle["chrms1"] = np.array([r[0] for r in rows])
        handle["cuts1"] = np.array([r[1] - 1 for r in rows])  # hiclib positions are zero-based
        handle["chrms2"] = np.array([r[2] for r in rows])
        handle["cuts2"] = np.array([r[3] - 1 for r in rows])
    result = cooler_cload_hiclib(
        small.bins_spec, str(tmp_path / "frags.h5"), str(tmp_path / "hiclib.cool")
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _pixel_counter(result["outputs"][0]) == small.expected


def test_cload_pairix_needs_the_index_and_reports_a_missing_pypairix(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    packed = tmp_path / "x.pairs.gz"
    packed.write_bytes(b"")
    with pytest.raises(ConfigurationError):
        cooler_cload_pairix(small.bins_spec, str(packed), str(tmp_path / "x.cool"))
    (tmp_path / "x.pairs.gz.px2").write_bytes(b"")
    if _cooler_program() is None or importlib.util.find_spec("pypairix") is not None:
        pytest.skip("needs the cooler program and no pypairix installed")
    result = cooler_cload_pairix(
        small.bins_spec,
        str(packed),
        str(tmp_path / "x.cool"),
        zero_based=True,
        nproc=1,
        max_split=2,
        block_char="|",
    )
    # every flag was accepted by cooler's option parser; the run stops at the missing Python package
    assert result["status"] == "failed" and "pypairix" in result["log_tail"]
    assert "No such option" not in result["log_tail"]


# ---------------------------------------- load, coarsen, merge


@requires_cooler
def test_load_reads_coo_and_bg2_text(small: SimpleNamespace, tmp_path: Path) -> None:
    coo = cooler_dump(str(small.cool), str(tmp_path / "pixels.coo"))
    loaded = cooler_load(
        small.bins_spec,
        coo["outputs"][0],
        str(tmp_path / "coo.cool"),
        input_format="coo",
        assembly="toy1",
    )
    assert loaded["status"] == "completed", loaded["log_tail"]
    assert _pixel_counter(loaded["outputs"][0]) == small.expected
    bg2 = cooler_dump(str(small.cool), str(tmp_path / "pixels.bg2.gz"), join=True)
    again = cooler_load(
        small.bins_spec,
        bg2["outputs"][0],
        str(tmp_path / "bg2.cool"),
        input_format="bg2",
        fields=["count=7:dtype=float"],
        chunksize=50,
    )
    assert again["status"] == "completed", again["log_tail"]
    clr = cooler.Cooler(again["outputs"][0])
    assert clr.pixels()[:]["count"].dtype.kind == "f" and _sum_and_nnz(again["outputs"][0]) == (
        600,
        len(small.expected),
    )
    with pytest.raises(ConfigurationError):
        cooler_load(
            small.bins_spec, coo["outputs"][0], str(tmp_path / "bad.cool"), input_format="tsv"
        )  # type: ignore[arg-type]


@requires_cooler
def test_coarsen_pools_pixels(small: SimpleNamespace, tmp_path: Path) -> None:
    result = cooler_coarsen(
        str(small.cool), str(tmp_path / "coarse.cool"), factor=2, nproc=1, fields=["count"]
    )
    assert result["status"] == "completed", result["log_tail"]
    clr = cooler.Cooler(result["outputs"][0])
    assert clr.binsize == 200 and int(clr.info["sum"]) == 600 and clr.info["nbins"] == 8
    inside = cooler_coarsen(
        str(small.cool), f"{tmp_path / 'coarse.cool'}::/factor4", factor=4, append=True
    )
    assert inside["status"] == "completed", inside["log_tail"]
    assert cooler.Cooler(f"{tmp_path / 'coarse.cool'}::/factor4").binsize == 400
    with pytest.raises(ConfigurationError):
        cooler_coarsen(
            str(small.cool), str(tmp_path / "coarse.cool"), factor=2
        )  # exists, no append
    with pytest.raises(ConfigurationError):
        cooler_coarsen(str(small.cool), str(tmp_path / "other.cool"), factor=1)


@requires_cooler
def test_merge_sums_identical_coolers(small: SimpleNamespace, tmp_path: Path) -> None:
    other = _make_cool(
        tmp_path / "other.cool",
        SMALL_SIZES,
        100,
        _expected_pixels(_random_pairs(SMALL_SIZES, 400, 9), SMALL_SIZES, 100),
    )
    result = cooler_merge(
        str(tmp_path / "merged.cool"), [str(small.cool), str(other)], chunksize=1000
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _sum_and_nnz(result["outputs"][0])[0] == 1000
    assert _pixel_counter(result["outputs"][0]) == small.expected + _pixel_counter(other)
    with pytest.raises(ConfigurationError):
        cooler_merge(str(tmp_path / "single.cool"), [str(small.cool)])
    with pytest.raises(ConfigurationError):
        cooler_merge(str(tmp_path / "single.cool"), str(small.cool))  # type: ignore[arg-type]


# ---------------------------------------- zoomify


@requires_cooler
def test_zoomify_explicit_binary_and_nice_progressions(
    big: SimpleNamespace, tmp_path: Path
) -> None:
    explicit = cooler_zoomify(
        str(big.cool), str(tmp_path / "explicit.mcool"), resolutions="100, 300,900"
    )
    assert explicit["status"] == "completed", explicit["log_tail"]
    levels = cooler_ls(explicit["outputs"][0], long=True)
    assert list(levels["bin_sizes"].values()) == [100, 300, 900]
    binary = cooler_zoomify(str(big.cool), str(tmp_path / "binary.mcool"))  # default progression B
    assert list(cooler_ls(binary["outputs"][0], long=True)["bin_sizes"].values()) == [100, 200, 400]
    nice = cooler_zoomify(str(big.cool), str(tmp_path / "nice.mcool"), resolutions="N", nproc=1)
    assert list(cooler_ls(nice["outputs"][0], long=True)["bin_sizes"].values()) == [100, 200, 500]
    for path in (explicit, binary, nice):
        assert int(cooler.Cooler(f"{path['outputs'][0]}::/resolutions/100").info["sum"]) == 5000
    with pytest.raises(ConfigurationError):
        cooler_zoomify(str(big.cool), str(tmp_path / "bad.mcool"), resolutions="100;rm")
    with pytest.raises(ConfigurationError):
        cooler_zoomify(str(big.cool), str(tmp_path / "bad.mcool"), balance_options={"min_nnz": 1})
    with pytest.raises(ConfigurationError):
        cooler_zoomify(
            str(big.cool), str(tmp_path / "bad.mcool"), balance=True, balance_options={"name": "x"}
        )


@requires_cooler
def test_zoomify_start_binary_progression_is_passed_through_and_reported(
    big: SimpleNamespace, tmp_path: Path
) -> None:
    # cooler 0.10.4 documents '<start>B' but cannot expand it ("int('100b')"); the facade passes
    # it on
    # unchanged and reports the outcome instead of raising.
    result = cooler_zoomify(str(big.cool), str(tmp_path / "startb.mcool"), resolutions="100B")
    assert result["command"][2:4] == ["--resolutions", "100B"]
    assert result["status"] in ("completed", "failed")
    if result["status"] == "failed":
        assert "100b" in result["log_tail"]


@requires_cooler
def test_zoomify_legacy_layout_and_balanced_levels(
    big: SimpleNamespace, small: SimpleNamespace, tmp_path: Path
) -> None:
    legacy = cooler_zoomify(str(big.cool), str(tmp_path / "legacy.mcool"), legacy=True)
    assert legacy["status"] == "completed", legacy["log_tail"]
    listed = cooler_ls(legacy["outputs"][0], long=True)["bin_sizes"]
    assert sorted(listed.values()) == [100, 200, 400, 800]
    assert all(uri.rsplit("::/", 1)[1] in {"0", "1", "2", "3"} for uri in listed)
    attrs = cooler_attrs(legacy["outputs"][0], level=0)["attributes"]
    assert attrs["@attrs"]["max-zoom"] == 3
    balanced = cooler_zoomify(
        str(small.cool),
        str(tmp_path / "balanced.mcool"),
        resolutions="100,200",
        balance=True,
        balance_options={"min_nnz": 1, "mad_max": 0, "nproc": 1},
    )
    assert balanced["status"] == "completed", balanced["log_tail"]
    for resolution in (100, 200):
        assert (
            cooler_balance_check(f"{balanced['outputs'][0]}::/resolutions/{resolution}")["balanced"]
            is True
        )


# ---------------------------------------- file operations


@requires_cooler
def test_ls_cp_ln_mv_tree_and_attrs(small: SimpleNamespace, tmp_path: Path) -> None:
    store = tmp_path / "store.cool"
    copied = cooler_cp(str(small.cool), str(store))
    assert copied["status"] == "completed", copied["log_tail"]
    assert cooler_ls(str(store))["coolers"] == [f"{store}::/"]
    assert cooler_ls(str(store), long=True)["bin_sizes"] == {f"{store}::/": 100}
    into = cooler_cp(f"{store}::/", f"{store}::/copy1")
    assert into["status"] == "completed", into["log_tail"]
    assert cooler_ln(f"{store}::/copy1", f"{store}::/hard")["status"] == "completed"
    assert cooler_ln(f"{store}::/copy1", f"{store}::/soft", soft=True)["status"] == "completed"
    assert cooler_mv(f"{store}::/hard", f"{store}::/moved")["status"] == "completed"
    names = sorted(uri.rsplit("::", 1)[1] for uri in cooler_ls(str(store))["coolers"])
    assert names == ["/", "/copy1", "/moved", "/soft"]
    clash = cooler_mv(f"{store}::/moved", f"{store}::/soft")
    assert clash["status"] == "failed" and "already exists" in clash["log_tail"]
    tree = cooler_tree(f"{store}::/copy1")
    assert tree["status"] == "completed" and "pixels" in tree["tree"] and "bin1_id" in tree["tree"]
    shallow = cooler_tree(str(store), level=1)["tree"]
    assert "bin1_id" not in shallow and "copy1" in shallow
    attrs = cooler_attrs(f"{store}::/copy1", level=0)
    assert (
        attrs["attributes"]["@attrs"]["bin-size"] == 100
        and "bin-size: 100" in attrs["attributes_yaml"]
    )
    other = tmp_path / "other.cool"
    assert cooler_cp(f"{store}::/copy1", f"{other}::/inner")["status"] == "completed"
    assert cooler_ls(str(other))["coolers"] == [f"{other}::/inner"]
    with pytest.raises(ConfigurationError):
        cooler_cp(str(tmp_path / "missing.cool"), str(tmp_path / "x.cool"))
    with pytest.raises(ConfigurationError):
        cooler_ls(str(tmp_path / "missing.cool"))
    with pytest.raises(ConfigurationError):
        cooler_tree(str(store), level=-1)


# ---------------------------------------- plots


@requires_cooler
def test_show_writes_an_image_with_the_values_behind_it(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    result = cooler_show(
        str(small.cool),
        "chr1:0-600",
        str(tmp_path / "show.png"),
        dpi=60,
        cmap="viridis",
        scale="linear",
        zmin=0.0,
        zmax=5.0,
    )
    assert result["status"] == "completed", result["log_tail"]
    image = Path(result["image_path"])
    assert (
        image.read_bytes()[:8] == PNG_SIGNATURE
        and result["image_bytes"] == image.stat().st_size > 500
    )
    clr = cooler.Cooler(str(small.cool))
    expected = clr.matrix(balance=False).fetch("chr1:0-600")
    assert result["data"]["summary"]["shape"] == [6, 6]
    assert result["data"]["matrix"] == expected.astype(float).tolist()
    assert result["data"]["summary"]["sum"] == pytest.approx(float(expected.sum()))
    trans = cooler_show(str(small.cool), "chr1", str(tmp_path / "trans.svg"), region2="chr2:0-600")
    assert trans["status"] == "completed" and trans["data"]["summary"]["shape"] == [10, 6]
    assert b"<svg" in Path(trans["image_path"]).read_bytes()[:400]
    unbalanced = cooler_show(str(small.cool), "chr1", str(tmp_path / "balanced.png"), balanced=True)
    assert unbalanced["status"] == "failed" and unbalanced["image_path"] is None
    with pytest.raises(ConfigurationError):
        cooler_show(str(small.cool), "chr1", str(tmp_path / "show.txt"))
    with pytest.raises(ConfigurationError):
        cooler_show(str(small.cool), "-chr1", str(tmp_path / "other.png"))


def test_show_view_applies_the_interactive_window_rounding(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    pytest.importorskip("matplotlib")
    whole = cooler_show_view(str(small.cool), "chr1", str(tmp_path / "whole.png"))
    assert whole["status"] == "completed" and not whole["placeholder"] and not whole["view_applied"]
    assert Path(whole["image_path"]).read_bytes()[:8] == PNG_SIGNATURE
    assert whole["row_region"] == ["chr1", 0, 1000] and whole["data"]["summary"]["shape"] == [
        10,
        10,
    ]
    zoom = cooler_show_view(
        str(small.cool), "chr1", str(tmp_path / "zoom.png"), view=[130, 480, 250, 910], dpi=50
    )
    assert (
        zoom["view_applied"]
        and zoom["column_region"] == ["chr1", 100, 500]
        and zoom["row_region"] == ["chr1", 200, 1000]
    )
    expected = (
        cooler.Cooler(str(small.cool))
        .matrix(balance=False)
        .fetch(("chr1", 200, 1000), ("chr1", 100, 500))
    )
    assert zoom["data"]["matrix"] == expected.astype(float).tolist()
    trimmed = cooler_show_view(
        str(small.cool), "chr2", str(tmp_path / "trim.png"), view=[-50, 5000, 0, 5000]
    )
    assert trimmed["column_region"] == ["chr2", 0, 600] and trimmed["row_region"] == [
        "chr2",
        0,
        600,
    ]
    last_bin = cooler_show_view(
        str(small.cool), "chr1", str(tmp_path / "last.png"), region2="chr2", view=[0, 599, 0, 1000]
    )
    assert last_bin["column_region"] == ["chr2", 0, 600]
    too_large = cooler_show_view(
        str(small.cool), "chr1", str(tmp_path / "large.png"), max_matrix_size=50
    )
    assert too_large["placeholder"] is True and too_large["data"]["matrix"] is None
    assert Path(too_large["image_path"]).stat().st_size > 500
    with pytest.raises(ConfigurationError):
        cooler_show_view(str(small.cool), "chr1", str(tmp_path / "bad.png"), view=[10, 10, 0, 100])
    with pytest.raises(ConfigurationError):
        cooler_show_view(str(small.cool), "chr1", str(tmp_path / "bad.png"), view=[1, 2, 3])
    with pytest.raises(ConfigurationError):
        cooler_show_view(str(small.cool), "chr1", str(tmp_path / "bad.png"), balanced=True)
    with pytest.raises(ConfigurationError):
        cooler_show_view(str(small.cool), "chr1", str(tmp_path / "bad.png"), scale="sqrt")  # type: ignore[arg-type]


# ---------------------------------------- pairs header


def test_pairs_header_reads_plain_and_gzipped_files(small: SimpleNamespace, tmp_path: Path) -> None:
    result = cooler_pairs_header(str(small.pairs_file))
    assert result["header"] == [
        "## pairs format v1.0",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 score",
    ]
    assert result["header_line_count"] == 2 and result["has_data_lines"] is True
    packed = tmp_path / "pairs.txt.gz"
    with gzip.open(packed, "wt") as handle:
        handle.write("@sq\tchr1  \n@sq\tchr2\nr1\tchr1\t1\n")
    from_gzip = cooler_pairs_header(str(packed), comment_char="@")
    assert from_gzip["header"] == ["@sq\tchr1", "@sq\tchr2"] and from_gzip["has_data_lines"] is True
    only_header = tmp_path / "only.txt"
    only_header.write_text("#a\n#b\n")
    assert cooler_pairs_header(str(only_header))["has_data_lines"] is False
    plain = tmp_path / "plain.txt"
    plain.write_text("r1\tchr1\n")
    assert cooler_pairs_header(str(plain))["header"] == []
    with pytest.raises(ConfigurationError):
        cooler_pairs_header(str(plain), comment_char="")


# ---------------------------------------- csort steps

CHROMS_TEXT = "chr1\t1000\nchr2\t600\nchr10\t500\n"
FLIP_TEXT = (
    "#header line\n"
    "r1\tchr2\t50\tchr1\t900\t+\t-\n"  # chromosome order reversed: flipped
    "r2\tchr1\t700\tchr1\t20\t+\t-\n"  # same chromosome, positions reversed: flipped
    "r3\tchr1\t10\tchr10\t30\t+\t-\n"  # already upper triangular
    "r4\tchrUn\t10\tchr1\t30\t+\t-\n"  # unknown chromosome: dropped
    "r5\tchr10\t5\tchr2\t8\t-\t-\n"  # chr10 is listed after chr2: flipped
    "r6\tchr2\t9\tchr2\t9\t+\t+\n"  # identical positions: kept as is
)


def test_flip_pairs_matches_the_flip_stage_of_cooler_csort(tmp_path: Path) -> None:
    csort_module = pytest.importorskip("cooler.cli.csort")
    chromosomes = tmp_path / "chroms.txt"
    chromosomes.write_text(CHROMS_TEXT)
    pairs = tmp_path / "pairs.txt"
    pairs.write_text(FLIP_TEXT)
    command = csort_module.make_flip_command(
        str(chromosomes), "\t", "#", {"C1": 2, "P1": 3, "C2": 4, "P2": 5}
    )
    expected = subprocess.run(
        command, input=FLIP_TEXT.encode(), capture_output=True, check=True
    ).stdout
    result = cooler_csort_flip_pairs(
        str(pairs),
        str(chromosomes),
        str(tmp_path / "flipped.txt"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
    )
    assert Path(result["outputs"][0]).read_bytes() == expected
    assert (
        result["lines_read"],
        result["comment_lines"],
        result["dropped"],
        result["flipped"],
        result["written"],
    ) == (7, 1, 1, 3, 5)
    lines = expected.decode().splitlines()
    assert lines[0] == "r1\tchr1\t900\tchr2\t50\t-\t+".replace("\t-\t+", "\t+\t-")
    packed = tmp_path / "pairs.txt.gz"
    with gzip.open(packed, "wb") as handle:
        handle.write(FLIP_TEXT.encode())
    from_gzip = cooler_csort_flip_pairs(
        str(packed),
        str(chromosomes),
        str(tmp_path / "flipped-gz.txt"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
    )
    assert Path(from_gzip["outputs"][0]).read_bytes() == expected


def test_flip_pairs_checks_columns_and_positions(tmp_path: Path) -> None:
    chromosomes = tmp_path / "chroms.txt"
    chromosomes.write_text(CHROMS_TEXT)
    short = tmp_path / "short.txt"
    short.write_text("r1\tchr1\t5\n")
    with pytest.raises(ConfigurationError):
        cooler_csort_flip_pairs(
            str(short),
            str(chromosomes),
            str(tmp_path / "out.txt"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    assert not (tmp_path / "out.txt").exists()
    text = tmp_path / "text.txt"
    text.write_text("r1\tchr1\tx\tchr1\t5\n")
    with pytest.raises(ConfigurationError):
        cooler_csort_flip_pairs(
            str(text),
            str(chromosomes),
            str(tmp_path / "out.txt"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    empty = tmp_path / "empty.chroms"
    empty.write_text("\n\n")
    with pytest.raises(ConfigurationError):
        cooler_csort_flip_pairs(
            str(text), str(empty), str(tmp_path / "out2.txt"), chrom1=2, pos1=3, chrom2=4, pos2=5
        )
    with pytest.raises(ConfigurationError):
        cooler_csort_flip_pairs(
            str(text),
            str(chromosomes),
            str(tmp_path / "out3.txt"),
            chrom1=2,
            pos1=2,
            chrom2=4,
            pos2=5,
        )


SORT_ROWS = [
    ("chr2", 40, "chr2", 90),
    ("chr10", 5, "chr10", 7),
    ("chr1", 300, "chr2", 10),
    ("chr1", 30, "chr1", 500),
    ("chr1", 4, "chr2", 700),
    ("chr1", 30, "chr1", 60),
    ("chr1", 30, "chr1", 60),
    ("chr2", 40, "chr2", 8),
    ("chr1", 100, "chr10", 1),
]


@pytest.mark.skipif(shutil.which("sort") is None, reason="the system sort program was not found")
@pytest.mark.parametrize("index", ["tabix", "pairix"])
def test_sort_pairs_uses_cooler_csort_keys_in_the_c_locale(tmp_path: Path, index: str) -> None:
    source = tmp_path / "flipped.txt"
    lines = [f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}" for i, (c1, p1, c2, p2) in enumerate(SORT_ROWS)]
    source.write_text("\n".join(lines) + "\n")
    result = cooler_csort_sort_pairs(
        str(source),
        str(tmp_path / "sorted.txt"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        index=index,
        nproc=2,
        buffer_size="50%",
        temp_dir=str(tmp_path),
    )  # type: ignore[arg-type]
    assert result["status"] == "completed", result["log_tail"]
    parsed = [
        (
            line.split("\t")[1].encode(),
            int(line.split("\t")[2]),
            line.split("\t")[3].encode(),
            int(line.split("\t")[4]),
            line.encode(),
        )
        for line in lines
    ]
    key = (
        (lambda r: (r[0], r[1], r[2], r[3], r[4]))
        if index == "tabix"
        else (lambda r: (r[0], r[2], r[1], r[3], r[4]))
    )
    expected = [row[4].decode() for row in sorted(parsed, key=key)]
    assert Path(result["outputs"][0]).read_text().splitlines() == expected
    assert (
        result["command"][0] == "sort"
        and "--parallel=2" in result["command"]
        and "--buffer-size=50%" in result["command"]
    )
    first_chromosomes = [row.split("\t")[1] for row in expected]
    assert first_chromosomes.index("chr10") < first_chromosomes.index(
        "chr2"
    )  # byte order, not natural order


def test_sort_and_index_commands_equal_the_ones_cooler_builds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    csort_module = pytest.importorskip("cooler.cli.csort")
    monkeypatch.setenv("LC_ALL", "C")  # make_sort_command sets this variable in the parent process
    fields = {"C1": 2, "P1": 3, "C2": 4, "P2": 5}
    for index in ("tabix", "pairix"):
        assert ["sort", *gap2._sort_keys(index, fields)] == csort_module.make_sort_command(
            index, fields, []
        )
        for zero_based in (False, True):
            ours = [index, *gap2._index_arguments(index, fields, zero_based, "out.gz")]
            assert ours == csort_module.make_index_command(index, fields, zero_based, "out.gz")
    shifted = {"C1": 5, "P1": 1, "C2": 3, "P2": 2}
    assert ["sort", *gap2._sort_keys("pairix", shifted)] == csort_module.make_sort_command(
        "pairix", shifted, []
    )


def test_sort_pairs_validates_arguments(tmp_path: Path) -> None:
    source = tmp_path / "in.txt"
    source.write_text("r\tchr1\t1\tchr1\t2\n")
    with pytest.raises(ConfigurationError):
        cooler_csort_sort_pairs(
            str(source),
            str(tmp_path / "o.txt"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            buffer_size="1G; ls",
        )
    with pytest.raises(ConfigurationError):
        cooler_csort_sort_pairs(
            str(source), str(tmp_path / "o.txt"), chrom1=2, pos1=3, chrom2=4, pos2=5, index="bed"
        )  # type: ignore[arg-type]
    with pytest.raises(BackendUnavailableError):
        cooler_csort_sort_pairs(
            str(source),
            str(tmp_path / "o.txt"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            sort_path=str(tmp_path / "no-sort"),
        )


def test_read_command_prefers_pigz_then_gzip() -> None:
    assert gap2._read_command("a.gz", "/bin/pigz", "/bin/gzip") == ["/bin/pigz", "-dc", "a.gz"]
    assert gap2._read_command("a.gz", None, "/bin/gzip") == ["/bin/gzip", "-dc", "a.gz"]
    with pytest.raises(BackendUnavailableError) as info:
        gap2._read_command("a.gz", None, None)
    assert "pigz" in (info.value.hint or "")


@pytest.mark.skipif(shutil.which("gzip") is None, reason="the gzip program was not found")
def test_read_pairs_decompresses_gz_and_copies_plain_text(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    text = "".join(f"r{i}\tchr1\t{i}\tchr1\t{i + 1}\n" for i in range(50))
    packed = tmp_path / "pairs.txt.gz"
    with gzip.open(packed, "wt") as handle:
        handle.write(text)
    monkeypatch.setattr(
        shutil, "which", lambda name, *a, **k: "/usr/bin/gzip" if name == "gzip" else None
    )
    result = cooler_csort_read_pairs(str(packed), str(tmp_path / "plain.txt"))
    assert result["status"] == "completed" and result["command"] == [
        "/usr/bin/gzip",
        "-dc",
        str(packed),
    ]
    assert Path(result["outputs"][0]).read_text() == text and result["bytes_written"] == len(text)
    copy = cooler_csort_read_pairs(str(tmp_path / "plain.txt"), str(tmp_path / "copy.txt"))
    assert copy["command"] is None and Path(copy["outputs"][0]).read_text() == text
    broken = tmp_path / "broken.txt.gz"
    broken.write_bytes(b"not gzip data")
    failed = cooler_csort_read_pairs(str(broken), str(tmp_path / "broken.txt"))
    assert (
        failed["status"] == "failed"
        and failed["outputs"] == []
        and not (tmp_path / "broken.txt").exists()
    )
    with pytest.raises(BackendUnavailableError):
        cooler_csort_read_pairs(
            str(packed), str(tmp_path / "x.txt"), gzip_path=str(tmp_path / "no-gzip")
        )


@pytest.mark.skipif(
    shutil.which("pigz") is None, reason="the pigz program was not found (gzip is used instead)"
)
def test_read_pairs_with_pigz(tmp_path: Path) -> None:
    packed = tmp_path / "pairs.txt.gz"
    with gzip.open(packed, "wt") as handle:
        handle.write("a\tb\n" * 20)
    result = cooler_csort_read_pairs(str(packed), str(tmp_path / "plain.txt"))
    assert (
        Path(result["command"][0]).name == "pigz"
        and Path(result["outputs"][0]).read_text() == "a\tb\n" * 20
    )


def test_index_pairs_runs_the_program_with_cooler_arguments(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    plain = tmp_path / "sorted.txt"
    plain.write_text("r\tchr1\t5\tchr1\t9\n")
    pysam.tabix_compress(str(plain), str(tmp_path / "sorted.txt.gz"))
    fake = tmp_path / "fake-index"
    fake.write_text(
        "#!/bin/sh\nfor last; do :; done\n"
        'printf "%s\\n" "$@" > "$last.args"\n'
        ': > "$last.px2"\n: > "$last.tbi"\n'
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    packed = tmp_path / "sorted.txt.gz"
    result = cooler_csort_index_pairs(
        str(packed),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        index="pairix",
        zero_based=True,
        index_path=str(fake),
    )
    assert result["status"] == "completed" and result["index_file"] == f"{packed}.px2"
    assert Path(f"{packed}.args").read_text().split() == [
        "-f",
        "-s2",
        "-d4",
        "-b3",
        "-e3",
        "-u5",
        "-v5",
        "-0",
        str(packed),
    ]
    assert result["command"][1:] == [
        "-f",
        "-s2",
        "-d4",
        "-b3",
        "-e3",
        "-u5",
        "-v5",
        "-0",
        str(packed),
    ]
    tabix = cooler_csort_index_pairs(
        str(packed), chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix", index_path=str(fake)
    )
    assert Path(f"{packed}.args").read_text().split() == ["-f", "-s2", "-b3", "-e3", str(packed)]
    assert tabix["index_file"] == f"{packed}.tbi"
    with pytest.raises(ConfigurationError):
        cooler_csort_index_pairs(
            str(plain), chrom1=2, pos1=3, chrom2=4, pos2=5, index_path=str(fake)
        )  # not bgzip data
    with pytest.raises(BackendUnavailableError):
        cooler_csort_index_pairs(
            str(packed), chrom1=2, pos1=3, chrom2=4, pos2=5, index_path=str(tmp_path / "no-tabix")
        )


@pytest.mark.skipif(
    shutil.which("tabix") is None, reason="the tabix program (htslib) was not found"
)
def test_index_pairs_with_the_real_tabix(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    plain = tmp_path / "sorted.txt"
    plain.write_text("r1\tchr1\t5\tchr1\t9\nr2\tchr1\t7\tchr2\t3\n")
    pysam.tabix_compress(str(plain), str(tmp_path / "sorted.txt.gz"))
    result = cooler_csort_index_pairs(
        str(tmp_path / "sorted.txt.gz"), chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix"
    )
    assert result["status"] == "completed", result["log_tail"]
    with pysam.TabixFile(str(tmp_path / "sorted.txt.gz")) as handle:
        assert len(list(handle.fetch("chr1", 0, 100))) == 2


def _stand_in_programs(folder: Path) -> Path:
    """bgzip and tabix stand-ins built on pysam (no htslib needed); they only drive cooler csort."""

    folder.mkdir()
    bgzip = folder / "bgzip"
    bgzip.write_text(
        f"#!{sys.executable}\n"
        + textwrap.dedent("""
        import os, sys, tempfile
        import pysam
        with tempfile.TemporaryDirectory() as scratch:
            plain, packed = os.path.join(scratch, "in.txt"), os.path.join(scratch, "out.gz")
            with open(plain, "wb") as handle:
                handle.write(sys.stdin.buffer.read())
            pysam.tabix_compress(plain, packed, force=True)
            with open(packed, "rb") as handle:
                sys.stdout.buffer.write(handle.read())
        """)
    )
    tabix = folder / "tabix"
    tabix.write_text(
        f"#!{sys.executable}\n"
        + textwrap.dedent("""
        import sys
        import pysam
        args = sys.argv[1:]
        numbers = {a[1]: int(a[2:]) for a in args if a[:2] in ("-s", "-b", "-e")}
        pysam.tabix_index(args[-1], force=True, seq_col=numbers["s"] - 1,
                          start_col=numbers["b"] - 1, end_col=numbers["e"] - 1,
                          zerobased="-0" in args)
        """)
    )
    for script in (bgzip, tabix):
        script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return folder


@requires_cooler
@pytest.mark.skipif(shutil.which("sort") is None, reason="the system sort program was not found")
def test_csort_end_to_end_then_cload_tabix(
    small: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("dill")
    if shutil.which("bgzip") is None or shutil.which("tabix") is None:
        stand_ins = _stand_in_programs(tmp_path / "standins")
        monkeypatch.setenv("PATH", f"{stand_ins}{os.pathsep}{os.environ.get('PATH', '')}")
    result = cooler_csort(
        str(small.pairs_file),
        str(small.sizes_file),
        str(tmp_path / "sorted.pairs.gz"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        index="tabix",
        nproc=1,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["outputs"][-1].endswith(".tbi") and result["command"][:2] == ["cooler", "csort"]
    with gzip.open(result["outputs"][0], "rt") as handle:
        rows = [line.split("\t") for line in handle.read().splitlines()]
    assert len(rows) == 600
    rank = {name: index for index, name in enumerate(small.sizes)}
    keys = [(rank[r[1]], int(r[2]), rank[r[3]], int(r[4])) for r in rows]
    assert keys == sorted(keys) and all((k[0], k[1]) <= (k[2], k[3]) for k in keys)
    with pysam.TabixFile(result["outputs"][0]) as handle:
        assert handle.contigs == ["chr1", "chr2"]
    loaded = cooler_cload_tabix(
        small.bins_spec,
        result["outputs"][0],
        str(tmp_path / "from-csort.cool"),
        chrom2=4,
        pos2=5,
        nproc=1,
    )
    assert loaded["status"] == "completed", loaded["log_tail"]
    assert _pixel_counter(loaded["outputs"][0]) == small.expected
    with pytest.raises(ConfigurationError):
        cooler_csort(
            str(small.pairs_file),
            str(small.sizes_file),
            result["outputs"][0],
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )


@requires_cooler
def test_csort_reports_missing_htslib_programs(small: SimpleNamespace, tmp_path: Path) -> None:
    if shutil.which("bgzip") is not None and shutil.which("pairix") is not None:
        pytest.skip("bgzip and pairix are installed")
    with pytest.raises(BackendUnavailableError) as info:
        cooler_csort(
            str(small.pairs_file),
            str(small.sizes_file),
            str(tmp_path / "sorted.gz"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
        )
    assert (info.value.hint or "").lower().startswith("install")
    with pytest.raises(ConfigurationError):
        cooler_csort(
            str(small.pairs_file),
            str(small.sizes_file),
            str(tmp_path / "sorted.gz"),
            chrom1=2,
            pos1=3,
            chrom2=3,
            pos2=5,
        )
    with pytest.raises(ConfigurationError):
        cooler_csort(
            str(small.pairs_file),
            str(small.sizes_file),
            str(tmp_path / "sorted.gz"),
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            sep="'",
        )


# ---------------------------------------- option coverage (each facade flag is run once)


def test_signatures_are_agent_registry_compatible() -> None:
    from dnakit.tools.adapters import analyze_callable

    for name in gap2.__all__:
        plan = analyze_callable(getattr(gap2, name), verbose=False)
        assert plan.compatible, (name, plan.incompatibilities)
        assert all(parameter.supported for parameter in plan.parameters), name


@requires_cooler
def test_cload_pairs_symmetry_metadata_storage_and_temporary_file_options(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    common = {"chrom1": 2, "pos1": 3, "chrom2": 4, "pos2": 5}
    square = cooler_cload_pairs(
        small.bins_spec,
        str(small.pairs_file),
        str(tmp_path / "square.cool"),
        comment_char="#",
        no_symmetric_upper=True,
        **common,
    )
    assert square["status"] == "completed", square["log_tail"]
    assert cooler.Cooler(square["outputs"][0]).info["storage-mode"] == "square"
    assert _sum_and_nnz(square["outputs"][0])[0] == 600
    duplex = cooler_cload_pairs(
        small.bins_spec,
        str(small.pairs_file),
        str(tmp_path / "duplex.cool"),
        input_copy_status="duplex",
        **common,
    )
    assert duplex["status"] == "completed", duplex["log_tail"]
    # records whose first mate lies after the second (by position, not by bin) are discarded
    rank = {name: index for index, name in enumerate(small.sizes)}
    kept = _expected_pixels(
        [pair for pair in small.pairs if (rank[pair[0]], pair[1]) <= (rank[pair[2]], pair[3])],
        small.sizes,
        small.binsize,
    )
    assert _pixel_counter(duplex["outputs"][0]) == kept
    metadata = tmp_path / "meta.json"
    metadata.write_text('{"project": "toy"}')
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    chunked = cooler_cload_pairs(
        small.bins_spec,
        str(small.pairs_file),
        str(tmp_path / "chunked.cool"),
        chunksize=100,
        mergebuf=50,
        max_merge=2,
        temp_dir=str(scratch),
        no_delete_temp=True,
        storage_options="compression=gzip,compression_opts=3",
        metadata_path=str(metadata),
        **common,
    )
    assert chunked["status"] == "completed", chunked["log_tail"]
    assert _pixel_counter(chunked["outputs"][0]) == small.expected
    assert cooler.Cooler(chunked["outputs"][0]).info["metadata"] == {"project": "toy"}
    assert any(scratch.iterdir())  # no_delete_temp keeps the partial coolers
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            small.bins_spec,
            str(small.pairs_file),
            str(tmp_path / "bad.cool"),
            storage_options="compression=gzip; ls",
            **common,
        )
    with pytest.raises(ConfigurationError):
        cooler_cload_pairs(
            small.bins_spec,
            str(small.pairs_file),
            str(tmp_path / "bad.cool"),
            temp_dir=str(tmp_path / "missing"),
            **common,
        )


def test_cload_tabix_and_hiclib_option_flags(small: SimpleNamespace, tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    h5py = pytest.importorskip("h5py")
    pytest.importorskip("dill")
    if _cooler_program() is None:
        pytest.skip("the cooler command-line program was not found")
    rows = _sorted_flipped_pairs(small)
    plain = tmp_path / "zero.pairs"
    plain.write_text(
        "".join(
            f"r{i}\t{c1}\t{p1 - 1}\t{c2}\t{p2 - 1}\n"
            for i, (*_, (c1, p1, c2, p2)) in enumerate(rows)
        )
    )
    pysam.tabix_compress(str(plain), str(tmp_path / "zero.pairs.gz"))
    pysam.tabix_index(
        str(tmp_path / "zero.pairs.gz"),
        seq_col=1,
        start_col=2,
        end_col=2,
        zerobased=True,
        force=True,
    )
    tabix = cooler_cload_tabix(
        small.bins_spec,
        str(tmp_path / "zero.pairs.gz"),
        str(tmp_path / "zero.cool"),
        chrom2=4,
        pos2=5,
        zero_based=True,
        nproc=2,
        max_split=3,
        assembly="toy1",
    )
    assert tabix["status"] == "completed", tabix["log_tail"]
    assert _pixel_counter(tabix["outputs"][0]) == small.expected
    with h5py.File(tmp_path / "frags.h5", "w") as handle:
        handle["chrms1"] = np.array([r[0] for r in rows])
        handle["cuts1"] = np.array([r[1] - 1 for r in rows])
        handle["chrms2"] = np.array([r[2] for r in rows])
        handle["cuts2"] = np.array([r[3] - 1 for r in rows])
    hiclib = cooler_cload_hiclib(
        small.bins_spec, str(tmp_path / "frags.h5"), str(tmp_path / "chunks.cool"), chunksize=100
    )
    assert hiclib["status"] == "completed", hiclib["log_tail"]
    assert _pixel_counter(hiclib["outputs"][0]) == small.expected


@requires_cooler
def test_load_options_one_based_comments_symmetry_and_storage(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    coo = tmp_path / "one-based.coo"
    coo.write_text(
        "% a header line\n"
        + "".join(f"{a + 1}\t{b + 1}\t{n}\n" for (a, b), n in sorted(small.expected.items()))
    )
    metadata = tmp_path / "meta.json"
    metadata.write_text('{"source": "coo"}')
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = cooler_load(
        small.bins_spec,
        str(coo),
        str(tmp_path / "options.cool"),
        input_format="coo",
        one_based=True,
        comment_char="%",
        chunksize=40,
        mergebuf=30,
        max_merge=2,
        temp_dir=str(scratch),
        no_delete_temp=True,
        storage_options="compression=gzip",
        metadata_path=str(metadata),
        input_copy_status="unique",
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _pixel_counter(result["outputs"][0]) == small.expected
    assert cooler.Cooler(result["outputs"][0]).info["metadata"] == {"source": "coo"}
    assert any(scratch.iterdir())
    square = cooler_load(
        small.bins_spec,
        str(tmp_path / "one-based.coo"),
        f"{tmp_path / 'options.cool'}::/square",
        input_format="coo",
        one_based=True,
        comment_char="%",
        no_symmetric_upper=True,
        append=True,
    )
    assert square["status"] == "completed", square["log_tail"]
    assert cooler.Cooler(f"{tmp_path / 'options.cool'}::/square").info["storage-mode"] == "square"


@requires_cooler
def test_balance_trans_only_distance_and_count_filters(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    target = _copy(small.cool, tmp_path / "trans.cool")
    result = cooler_balance(
        str(target),
        trans_only=True,
        ignore_dist=100,
        min_count=1,
        min_nnz=1,
        mad_max=0,
        chunksize=50,
        nproc=1,
        convergence_policy="store_nan",
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["weights_stored"] is True
    assert "--trans-only" in result["command"] and "--ignore-dist" in result["command"]


@requires_cooler
def test_dump_missing_value_format_second_range_and_chunk_options(
    small: SimpleNamespace, tmp_path: Path
) -> None:
    balanced_cool = _copy(small.cool, tmp_path / "weights.cool")
    blacklist = tmp_path / "bad.bed"
    blacklist.write_text("chr1\t0\t100\nchr2\t300\t400\n")
    cooler_balance(str(balanced_cool), min_nnz=1, mad_max=0, nproc=1, blacklist=str(blacklist))
    weights = cooler_dump(
        str(balanced_cool),
        str(tmp_path / "weights.tsv"),
        table="bins",
        columns=["weight"],
        header=True,
        na_rep="NA",
        float_format=".2f",
    )
    assert weights["status"] == "completed", weights["log_tail"]
    values = Path(weights["outputs"][0]).read_text().split()[1:]
    assert (
        values[0] == "NA"
        and values[13] == "NA"
        and all(len(v.split(".")[-1]) == 2 for v in values[1:13])
    )
    trans = cooler_dump(
        str(small.cool),
        str(tmp_path / "trans.tsv"),
        region="chr1",
        region2="chr2:0-600",
        chunksize=7,
        header=True,
    )
    table = _tsv(Path(trans["outputs"][0]))
    assert trans["status"] == "completed" and table["bin1_id"].max() < 10 <= table["bin2_id"].min()


@requires_cooler
def test_coarsen_merge_and_zoomify_option_flags(
    small: SimpleNamespace, big: SimpleNamespace, tmp_path: Path
) -> None:
    coarse = cooler_coarsen(
        str(small.cool),
        str(tmp_path / "coarse.cool"),
        factor=2,
        chunksize=40,
        nproc=2,
        fields=["count:dtype=float64,agg=sum"],
    )
    assert coarse["status"] == "completed", coarse["log_tail"]
    assert cooler.Cooler(coarse["outputs"][0]).pixels()[:]["count"].dtype.kind == "f"
    other = _make_cool(tmp_path / "other.cool", SMALL_SIZES, 100, small.expected)
    target = tmp_path / "merged.cool"
    first = cooler_merge(str(target), [str(small.cool), str(other)], fields=["count"])
    assert first["status"] == "completed", first["log_tail"]
    second = cooler_merge(
        f"{target}::/again",
        [str(small.cool), str(other)],
        fields=["count:dtype=int64"],
        append=True,
    )
    assert second["status"] == "completed", second["log_tail"]
    assert _sum_and_nnz(f"{target}::/again")[0] == 1200
    zoom = cooler_zoomify(
        str(small.cool),
        str(tmp_path / "options.mcool"),
        resolutions="100,200",
        chunksize=40,
        fields=["count:dtype=float64"],
        base_uris=[str(other)],
    )
    assert zoom["status"] == "completed", zoom["log_tail"]
    assert int(cooler.Cooler(f"{zoom['outputs'][0]}::/resolutions/100").info["sum"]) == 600
    assert str(big.cool)  # the larger dataset is only built once per module


@requires_cooler
def test_show_leaves_out_the_matrix_of_very_large_regions(
    small: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(gap2, "_DATA_CELL_LIMIT", 10)
    result = cooler_show(str(small.cool), "chr1", str(tmp_path / "large.png"))
    assert result["status"] == "completed" and Path(result["image_path"]).is_file()
    assert result["data"]["matrix"] is None and result["data"]["matrix_omitted"] is True
    assert "100 cells" in result["data"]["reason"]


@requires_cooler
def test_overwrite_flags_and_show_force_and_field(small: SimpleNamespace, tmp_path: Path) -> None:
    store = tmp_path / "store.cool"
    assert cooler_cp(str(small.cool), f"{store}::/a")["status"] == "completed"
    assert cooler_cp(str(small.cool), f"{store}::/b")["status"] == "completed"
    clash = cooler_cp(f"{store}::/a", f"{store}::/b")
    assert clash["status"] == "failed"
    # -w truncates the destination file: it replaces another file, and cooler 0.10.4 cannot do it
    # inside the file that holds the source ("unable to truncate a file which is already open")
    other = tmp_path / "other.cool"
    assert cooler_cp(str(small.cool), f"{other}::/old")["status"] == "completed"
    assert cooler_cp(str(small.cool), str(other), overwrite=True)["status"] == "completed"
    assert cooler_ls(str(other))["coolers"] == [f"{other}::/"]
    for result in (
        cooler_cp(f"{store}::/a", f"{store}::/b", overwrite=True),
        cooler_mv(f"{store}::/a", f"{store}::/c", overwrite=True),
        cooler_ln(f"{store}::/b", f"{store}::/d", overwrite=True, soft=True),
    ):
        assert result["status"] in ("completed", "failed")
        assert "No such option" not in result["log_tail"]
        if result["status"] == "failed":
            assert "already open" in result["log_tail"]
    shown = cooler_show(
        str(small.cool), "chr2", str(tmp_path / "forced.png"), force=True, field="count", zmax=3.0
    )
    assert shown["status"] == "completed", shown["log_tail"]
    assert shown["data"]["summary"]["shape"] == [6, 6]


@requires_cooler
@pytest.mark.skipif(shutil.which("sort") is None, reason="the system sort program was not found")
def test_csort_zero_based_separator_and_comment_options(
    small: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("pysam")
    if shutil.which("bgzip") is None or shutil.which("tabix") is None:
        stand_ins = _stand_in_programs(tmp_path / "standins")
        monkeypatch.setenv("PATH", f"{stand_ins}{os.pathsep}{os.environ.get('PATH', '')}")
    pairs = tmp_path / "semicolon.pairs"
    pairs.write_text(
        "!comment\n"
        + "".join(
            f"r{i};{c1};{p1 - 1};{c2};{p2 - 1}\n" for i, (c1, p1, c2, p2) in enumerate(small.pairs)
        )
    )
    result = cooler_csort(
        str(pairs),
        str(small.sizes_file),
        str(tmp_path / "semicolon.gz"),
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        index="tabix",
        zero_based=True,
        sep=";",
        comment_char="!",
        nproc=1,
    )
    assert result["status"] == "completed", result["log_tail"]
    with gzip.open(result["outputs"][0], "rt") as handle:
        assert len(handle.read().splitlines()) == 600
    assert "--zero-based" in result["command"] and "--sep" in result["command"]
