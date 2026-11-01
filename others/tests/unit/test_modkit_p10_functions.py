"""Tests for dnakit.epigenome.modkit_p10_functions (modkit subcommands as typed facades).

Most tests run a stand-in ``modkit`` executable written into tmp_path (it records its argument vector and writes the output
files the real program would write), so the argument mapping, the validation and the result shape are checked without the
real program. Tests marked ``requires_modkit`` run the real modkit and a small modBAM built with pysam; they skip when
modkit (https://github.com/nanoporetech/modkit/releases) is not installed.
"""

from __future__ import annotations

import json
import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import modkit_p10_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_STUB = r"""#!{python}
import json, os, sys
argv = sys.argv[1:]
here = os.path.dirname(os.path.abspath(__file__))
if argv == ["--version"]:
    print("modkit 0.6.4")
    sys.exit(0)
with open(os.path.join(here, "calls.jsonl"), "a") as handle:
    handle.write(json.dumps(argv) + "\n")
if os.path.exists(os.path.join(here, "FAIL")):
    print("Error! simulated failure")
    sys.exit(1)
nested = {{"motif", "bedmethyl", "modbam", "extract", "dmr", "open-chromatin"}}
path = argv[:2] if argv[0] in nested else argv[:1]
rest = argv[len(path):]
positional = []
for token in rest:
    if token.startswith("-"):
        break
    positional.append(token)
TABLE = "chrom\tstart\tend\tname\nchr1\t1\t2\tm\n"
def write(target, text=TABLE):
    if os.path.isdir(target):
        with open(os.path.join(target, "result.tsv"), "w") as out:
            out.write(text)
    else:
        with open(target, "w") as out:
            out.write(text)
if path in (["motif", "bed"], ["modbam", "summary"]):
    sys.stdout.write(TABLE)
    sys.exit(0)
if path in (["pileup"], ["extract", "full"], ["extract", "calls"], ["extract", "read-stats"], ["bedmethyl", "tobigwig"],
            ["bedmethyl", "map-to-genome"], ["modbam", "adjust-mods"], ["modbam", "update-tags"], ["modbam", "call-mods"],
            ["dmr", "isoform"]):
    write(positional[1])
for flag in ("--out-bed", "--out-table", "--out", "--output", "--out-file", "--out-filepath", "--out-path",
             "--output-bam", "--chart", "--plot", "--segment", "--out-dir"):
    if flag in rest:
        write(rest[rest.index(flag) + 1])
for index, token in enumerate(rest):
    if token == "-o" and index + 1 < len(rest):
        write(rest[index + 1])
"""


@pytest.fixture
def fake_modkit(tmp_path: Path) -> Path:
    directory = tmp_path / "bin"
    directory.mkdir()
    program = directory / "modkit"
    program.write_text(_STUB.format(python=sys.executable), encoding="utf-8")
    program.chmod(program.stat().st_mode | stat.S_IXUSR)
    return program


def _calls(program: Path) -> list[list[str]]:
    log = program.parent / "calls.jsonl"
    return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]


@pytest.fixture
def files(tmp_path: Path) -> dict[str, str]:
    """Small placeholder inputs; the stand-in program does not read them."""
    names = {
        "bam": "reads.bam",
        "bam2": "reads2.bam",
        "ref": "ref.fa",
        "bed": "regions.bed",
        "bm": "a.bed.gz",
        "bm2": "b.bed.gz",
        "sizes": "sizes.txt",
        "gtf": "genes.gtf",
        "model": None,
        "truth": "truth.bed",
    }
    out: dict[str, str] = {}
    for key, name in names.items():
        if name is None:
            continue
        path = tmp_path / name
        path.write_bytes(b"x")
        out[key] = str(path)
    model = tmp_path / "model"
    model.mkdir()
    out["model"] = str(model)
    return out


# ---------------------------------------------------------------------------
# Argument mapping: one case per public function
# ---------------------------------------------------------------------------


def _cases(f: dict[str, str], out: Path):
    return [
        (
            "modkit_pileup",
            (f["bam"], str(out / "p.bed")),
            dict(
                cpg=True,
                reference=f["ref"],
                modified_bases=["5mC"],
                threads=2,
                filter_percentile=20,
                combine_strands=True,
            ),
            [
                "pileup",
                f["bam"],
                str(out / "p.bed"),
                "--cpg",
                "--reference",
                f["ref"],
                "--threads",
                "2",
                "--filter-percentile",
                "20.0",
                "--combine-strands",
                "--modified-bases",
                "5mC",
            ],
        ),
        (
            "modkit_pileup_hemi",
            (f["bam"], f["ref"], str(out / "h.bed")),
            dict(cpg=True),
            ["pileup-hemi", f["bam"], "--ref", f["ref"], "--out-bed", str(out / "h.bed"), "--cpg"],
        ),
        (
            "modkit_validate",
            ([[f["bam"], f["truth"]]], str(out / "v.tsv")),
            dict(canonical_base="C", min_alignment_length=100),
            [
                "validate",
                "--bam-and-bed",
                f["bam"],
                f["truth"],
                "--canonical-base",
                "C",
                "--min-length",
                "100",
            ],
        ),
        (
            "modkit_motif_bed",
            (f["ref"], "CG", 0, str(out / "m.bed")),
            dict(mask=True),
            ["motif", "bed", f["ref"], "CG", "0", "--mask"],
        ),
        (
            "modkit_motif_search",
            (f["bm"], f["ref"], str(out / "s.tsv")),
            dict(min_coverage=5, known_motifs=[["CG", "0", "m"]]),
            [
                "motif",
                "search",
                "--in-bedmethyl",
                f["bm"],
                "--ref",
                f["ref"],
                "--known-motif",
                "CG",
                "0",
                "m",
                "--min-coverage",
                "5",
                "--out-table",
                str(out / "s.tsv"),
            ],
        ),
        (
            "modkit_motif_evaluate",
            (f["bm"], f["ref"], str(out / "e.tsv")),
            dict(known_motifs=[["GATC", "1", "a"]]),
            [
                "motif",
                "evaluate",
                "--in-bedmethyl",
                f["bm"],
                "--known-motif",
                "GATC",
                "1",
                "a",
                "--out",
                str(out / "e.tsv"),
            ],
        ),
        (
            "modkit_motif_refine",
            (f["bm"], f["ref"], str(out / "r.tsv")),
            dict(known_motifs=[["CG", "0", "m"]], min_refine_sites=3),
            ["motif", "refine", "--min-refine-sites", "3", "--out", str(out / "r.tsv")],
        ),
        (
            "modkit_entropy",
            ([f["bam"]], f["ref"], str(out / "en.bed")),
            dict(cpg=True, num_positions=3, window_size=40),
            [
                "entropy",
                "--in-bam",
                f["bam"],
                "--ref",
                f["ref"],
                "--cpg",
                "--num-positions",
                "3",
                "--window-size",
                "40",
                "--out-bed",
                str(out / "en.bed"),
            ],
        ),
        (
            "modkit_localize",
            (f["bm"], f["bed"], f["sizes"], str(out / "l.tsv")),
            dict(window=500, chart_path=str(out / "l.html")),
            [
                "localize",
                f["bm"],
                "--regions",
                f["bed"],
                "--genome-sizes",
                f["sizes"],
                "--window",
                "500",
                "--chart",
                str(out / "l.html"),
                "--out-file",
                str(out / "l.tsv"),
            ],
        ),
        (
            "modkit_stats",
            (f["bm"], f["bed"], str(out / "st.tsv")),
            dict(mod_codes="m,h", min_coverage=5),
            [
                "stats",
                f["bm"],
                "--regions",
                f["bed"],
                "--mod-codes",
                "m,h",
                "--min-coverage",
                "5",
                "--out-table",
                str(out / "st.tsv"),
            ],
        ),
        (
            "modkit_bedmethyl_merge",
            ([f["bm"], f["bm2"]], f["sizes"], str(out / "mg.bed")),
            dict(min_samples=2),
            [
                "bedmethyl",
                "merge",
                f["bm"],
                f["bm2"],
                "--genome-sizes",
                f["sizes"],
                "--min-samples",
                "2",
                "--out-bed",
                str(out / "mg.bed"),
            ],
        ),
        (
            "modkit_bedmethyl_to_bigwig",
            (f["bm"], str(out / "t.bw"), ["m"]),
            dict(chrom_sizes=f["sizes"]),
            [
                "bedmethyl",
                "tobigwig",
                f["bm"],
                str(out / "t.bw"),
                "--sizes",
                f["sizes"],
                "--mod-codes",
                "m",
            ],
        ),
        (
            "modkit_bedmethyl_map_to_genome",
            (f["bm"], str(out / "g.bed"), f["gtf"], "tx1"),
            dict(ignore_version=True),
            [
                "bedmethyl",
                "map-to-genome",
                f["bm"],
                str(out / "g.bed"),
                "--transcript-id",
                "tx1",
                "--ignore-version",
                "--gtf",
                f["gtf"],
            ],
        ),
        (
            "modkit_modbam_check_tags",
            (f["bam"], str(out / "chk")),
            dict(permissive=True, head=100),
            [
                "modbam",
                "check-tags",
                f["bam"],
                "--permissive",
                "--out-dir",
                str(out / "chk"),
                "--head",
                "100",
            ],
        ),
        (
            "modkit_modbam_adjust_mods",
            (f["bam"], str(out / "adj.bam")),
            dict(convert_from="h", convert_to="m"),
            ["modbam", "adjust-mods", f["bam"], str(out / "adj.bam"), "--convert", "h", "m"],
        ),
        (
            "modkit_modbam_update_tags",
            (f["bam"], str(out / "up.bam")),
            dict(mode="explicit"),
            ["modbam", "update-tags", f["bam"], str(out / "up.bam"), "--mode", "explicit"],
        ),
        (
            "modkit_modbam_sample_probs",
            (f["bam"], str(out / "sp")),
            dict(histogram=True, percentiles="10,50,90"),
            [
                "modbam",
                "sample-probs",
                f["bam"],
                "--percentiles",
                "10,50,90",
                "--out-dir",
                str(out / "sp"),
                "--hist",
            ],
        ),
        (
            "modkit_modbam_summary",
            (f["bam"], str(out / "sum.txt")),
            dict(tsv=True),
            ["modbam", "summary", f["bam"], "--tsv"],
        ),
        (
            "modkit_modbam_call_mods",
            (f["bam"], str(out / "call.bam")),
            dict(filter_threshold=0.8),
            ["modbam", "call-mods", f["bam"], str(out / "call.bam"), "--filter-threshold", "0.8"],
        ),
        (
            "modkit_modbam_repair",
            (f["bam"], f["bam2"], str(out / "rep.bam")),
            dict(),
            [
                "modbam",
                "repair",
                "--donor-bam",
                f["bam"],
                "--acceptor-bam",
                f["bam2"],
                "--output-bam",
                str(out / "rep.bam"),
            ],
        ),
        (
            "modkit_extract",
            (f["bam"], str(out / "x.tsv")),
            dict(table="calls", pass_only=True),
            ["extract", "calls", f["bam"], str(out / "x.tsv"), "--pass-only"],
        ),
        (
            "modkit_dmr",
            ({"a": f["bm"], "b": f["bm2"]}, f["bed"], f["ref"], str(out / "d.bed")),
            dict(),
            [
                "dmr",
                "pair",
                "-a",
                f["bm"],
                "-b",
                f["bm2"],
                "--out-path",
                str(out / "d.bed"),
                "--regions-bed",
                f["bed"],
                "--ref",
                f["ref"],
                "--base",
                "C",
            ],
        ),
        (
            "modkit_dmr_isoform",
            (f["bm"], f["gtf"], str(out / "i.tsv")),
            dict(gene_name="GENE1", full=True),
            [
                "dmr",
                "isoform",
                f["bm"],
                str(out / "i.tsv"),
                "--gtf",
                f["gtf"],
                "--full",
                "--gene-name",
                "GENE1",
            ],
        ),
        (
            "modkit_dmr_compare_tx_sites",
            (f["bm"], f["bm2"], f["gtf"], str(out / "c.tsv")),
            dict(top_k=10),
            [
                "dmr",
                "compare-tx-sites",
                "-a",
                f["bm"],
                "-b",
                f["bm2"],
                "--out",
                str(out / "c.tsv"),
                "--gtf",
                f["gtf"],
                "--top-k",
                "10",
            ],
        ),
        (
            "modkit_open_chromatin_predict",
            (f["bam"], f["model"], str(out / "o.bed")),
            dict(dry_run=True),
            [
                "open-chromatin",
                "predict",
                f["bam"],
                "--model",
                f["model"],
                "--output",
                str(out / "o.bed"),
                "--dryrun",
            ],
        ),
    ]


def _contiguous(argv: list[str], segment: list[str]) -> int | None:
    for index in range(len(argv) - len(segment) + 1):
        if argv[index : index + len(segment)] == segment:
            return index
    return None


def _assert_segments(name: str, argv: list[str], expected: list[str]) -> None:
    """The leading positional tokens must open the command; every ``--flag value...`` group must appear contiguously."""
    lead = []
    for token in expected:
        if token.startswith("-") and len(token) > 1 and not token[1:].isdigit():
            break
        lead.append(token)
    assert argv[: len(lead)] == lead, (name, lead, argv)
    segment: list[str] = []
    for token in [*expected[len(lead) :], "-"]:
        flag = token.startswith("-") and len(token) > 1 and not token[1:].replace(".", "").isdigit()
        if flag or token == "-":
            if segment:
                assert _contiguous(argv, segment) is not None, (name, segment, argv)
            segment = [token] if token != "-" else []
        else:
            segment.append(token)


def test_every_public_function_has_a_case(fake_modkit, files, tmp_path):
    cases = _cases(files, tmp_path)
    assert sorted(name for name, *_ in cases) == sorted(m.__all__)


def test_each_function_maps_its_arguments_to_modkit_flags(fake_modkit, files, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    for name, args, kwargs, expected in _cases(files, out):
        before = len(_calls(fake_modkit)) if (fake_modkit.parent / "calls.jsonl").exists() else 0
        result = getattr(m, name)(*args, modkit_path=str(fake_modkit), verbose=False, **kwargs)
        assert result["status"] == "completed", (name, result["log_tail"])
        assert result["operation"] == name
        assert result["backend"] == "modkit" and result["version"] == "0.6.4"
        assert result["provenance"]["tested_version"] == "0.6.4"
        json.dumps(result)
        argv = _calls(fake_modkit)[before]
        _assert_segments(name, argv, expected)
        # modkit's multi-value options must come after the positionals (they would swallow them otherwise)
        if "--modified-bases" in argv:
            assert argv.index("--modified-bases") > argv.index(args[1])


def test_tables_are_previewed(fake_modkit, files, tmp_path):
    result = m.modkit_stats(
        files["bm"],
        files["bed"],
        str(tmp_path / "stats.tsv"),
        modkit_path=str(fake_modkit),
        verbose=False,
    )
    assert result["columns"] == ["chrom", "start", "end", "name"]
    assert result["row_count"] == 1 and result["rows"] == [["chr1", "1", "2", "m"]]
    motifs = m.modkit_motif_bed(
        files["ref"], "CG", 0, str(tmp_path / "cg.bed"), modkit_path=str(fake_modkit), verbose=False
    )
    assert (
        motifs["record_count"] == 2
    )  # the stand-in prints a header line and one row, no comment prefix


def test_directory_outputs_are_listed(fake_modkit, files, tmp_path):
    result = m.modkit_dmr(
        {"a": files["bm"], "b": files["bm2"], "c": files["bm"]},
        files["bed"],
        files["ref"],
        str(tmp_path / "multi"),
        mode="multi",
        prefix="cmp",
        modkit_path=str(fake_modkit),
        verbose=False,
    )
    assert result["mode"] == "multi" and result["files"] == [str(tmp_path / "multi" / "result.tsv")]
    argv = _calls(fake_modkit)[-1]
    assert argv[:2] == ["dmr", "multi"] and argv.count("--sample") == 3 and "--prefix" in argv


def test_extract_read_stats_needs_mod_codes_and_rejects_other_options(fake_modkit, files, tmp_path):
    with pytest.raises(ConfigurationError, match="mod_codes"):
        m.modkit_extract(
            files["bam"],
            str(tmp_path / "rs.tsv"),
            table="read-stats",
            modkit_path=str(fake_modkit),
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="does not take"):
        m.modkit_extract(
            files["bam"],
            str(tmp_path / "rs.tsv"),
            table="read-stats",
            mod_codes=["C:m"],
            region="chr1:1-10",
            modkit_path=str(fake_modkit),
            verbose=False,
        )
    result = m.modkit_extract(
        files["bam"],
        str(tmp_path / "rs.csv"),
        table="read-stats",
        mod_codes=["C:m"],
        filter_threshold=0.7,
        modkit_path=str(fake_modkit),
        verbose=False,
    )
    argv = _calls(fake_modkit)[-1]
    assert argv[:2] == ["extract", "read-stats"] and "--filter-threshold" in argv
    assert argv.index("--mod-codes") > argv.index(str(tmp_path / "rs.csv"))
    assert result["table"] == "read-stats"
    with pytest.raises(ConfigurationError, match="threshold"):
        m.modkit_extract(
            files["bam"],
            str(tmp_path / "full.tsv"),
            table="full",
            filter_threshold=0.5,
            modkit_path=str(fake_modkit),
            verbose=False,
        )


# ---------------------------------------------------------------------------
# Validation, failure and missing program
# ---------------------------------------------------------------------------


def test_bad_inputs_are_rejected_before_modkit_runs(fake_modkit, files, tmp_path):
    kw = dict(modkit_path=str(fake_modkit), verbose=False)
    with pytest.raises(ConfigurationError, match="existing file"):
        m.modkit_pileup(str(tmp_path / "missing.bam"), str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="reference is required"):
        m.modkit_pileup(files["bam"], str(tmp_path / "x.bed"), cpg=True, **kw)
    with pytest.raises(ConfigurationError, match="Give exactly one"):
        m.modkit_pileup_hemi(files["bam"], files["ref"], str(tmp_path / "h.bed"), **kw)
    with pytest.raises(ConfigurationError, match="offset"):
        m.modkit_motif_bed(files["ref"], "CG", 5, str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="letters only"):
        m.modkit_motif_bed(files["ref"], "C-G", 0, str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="must be one of"):
        m.modkit_modbam_update_tags(files["bam"], str(tmp_path / "x.bam"), mode="maybe", **kw)
    with pytest.raises(ConfigurationError, match="at least two"):
        m.modkit_bedmethyl_merge([files["bm"]], files["sizes"], str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="exactly one of chrom_sizes"):
        m.modkit_bedmethyl_to_bigwig(files["bm"], str(tmp_path / "x.bw"), ["m"], **kw)
    with pytest.raises(ConfigurationError, match="together"):
        m.modkit_modbam_adjust_mods(files["bam"], str(tmp_path / "x.bam"), convert_from="h", **kw)
    with pytest.raises(ConfigurationError, match="exactly two"):
        m.modkit_dmr({"a": files["bm"]}, files["bed"], files["ref"], str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="model directory"):
        m.modkit_open_chromatin_predict(files["bam"], files["bed"], str(tmp_path / "x.bed"), **kw)
    with pytest.raises(ConfigurationError, match="must be an integer"):
        m.modkit_pileup(files["bam"], str(tmp_path / "x.bed"), threads="4", **kw)
    with pytest.raises(ConfigurationError, match=r"m:0\.8"):
        m.modkit_pileup(files["bam"], str(tmp_path / "x.bed"), mod_thresholds=["bad"], **kw)
    assert not (fake_modkit.parent / "calls.jsonl").exists()


def test_existing_outputs_are_never_overwritten(fake_modkit, files, tmp_path):
    existing = tmp_path / "exists.bed"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="must not exist"):
        m.modkit_pileup(files["bam"], str(existing), modkit_path=str(fake_modkit), verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep"
    (tmp_path / "dir").mkdir()
    with pytest.raises(ConfigurationError, match="must not exist"):
        m.modkit_modbam_check_tags(
            files["bam"], str(tmp_path / "dir"), modkit_path=str(fake_modkit), verbose=False
        )


def test_failed_run_is_reported_not_raised(fake_modkit, files, tmp_path):
    (fake_modkit.parent / "FAIL").write_text("", encoding="utf-8")
    result = m.modkit_stats(
        files["bm"],
        files["bed"],
        str(tmp_path / "s.tsv"),
        modkit_path=str(fake_modkit),
        verbose=False,
    )
    assert result["status"] == "failed" and result["return_code"] == 1
    assert "simulated failure" in result["log_tail"] and "row_count" not in result


def test_missing_program_raises_with_install_hint_and_creates_nothing(files, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setattr(sys, "executable", str(tmp_path / "nowhere" / "python"))
    target = tmp_path / "checks"
    with pytest.raises(BackendUnavailableError) as error:
        m.modkit_modbam_check_tags(files["bam"], str(target), verbose=False)
    assert "modkit" in str(error.value) and not target.exists()
    with pytest.raises(BackendUnavailableError):
        m.modkit_pileup(
            files["bam"],
            str(tmp_path / "x.bed"),
            modkit_path=str(tmp_path / "no_modkit"),
            verbose=False,
        )


def test_public_functions_carry_the_verbosity_parameter():
    import inspect

    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True


# ---------------------------------------------------------------------------
# Real modkit on a small modBAM
# ---------------------------------------------------------------------------


def _modkit_installed() -> bool:
    return shutil.which("modkit") is not None or (Path(sys.executable).parent / "modkit").is_file()


requires_modkit = pytest.mark.skipif(not _modkit_installed(), reason="modkit is not installed")


def _toy_modbam(directory: Path) -> tuple[Path, Path]:
    pysam = pytest.importorskip("pysam")
    import array

    rng = random.Random(7)
    pieces = []
    while sum(len(p) for p in pieces) < 600:
        pieces.append(
            rng.choice(["CG", "ACGT", "TCGA", "GCGC", "AT", "CCGG"])
            + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4)))
        )
    reference = "".join(pieces)[:600]
    ref = directory / "ref.fa"
    ref.write_text(f">chr1\n{reference}\n", encoding="utf-8")
    pysam.faidx(str(ref))
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": len(reference)}]}
    unsorted = directory / "u.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as handle:
        for number, start in enumerate(range(0, 300, 20)):
            read = reference[start : start + 200]
            segment = pysam.AlignedSegment()
            segment.query_name = f"r{number:03d}"
            segment.query_sequence = read
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigartuples = [(0, len(read))]
            segment.query_qualities = pysam.qualitystring_to_array("I" * len(read))
            cs = [i for i, base in enumerate(read) if base == "C"]
            cpg = [i for i in range(len(read) - 1) if read[i] == "C" and read[i + 1] == "G"]
            deltas, last, probabilities = [], -1, []
            for position in cpg:
                index = cs.index(position)
                deltas.append(index - (last + 1))
                last = index
                probabilities.append(
                    rng.randint(200, 255) if rng.random() < 0.8 else rng.randint(0, 40)
                )
            segment.set_tag("NM", 0)
            segment.set_tag("MM", "C+m?," + ",".join(map(str, deltas)) + ";", "Z")
            segment.set_tag("ML", array.array("B", probabilities))
            handle.write(segment)
    bam = directory / "reads.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam, ref


@requires_modkit
def test_real_pileup_stats_and_motif_bed(tmp_path):
    pysam = pytest.importorskip("pysam")
    bam, ref = _toy_modbam(tmp_path)
    pileup = m.modkit_pileup(
        str(bam),
        str(tmp_path / "cpg.bed.gz"),
        reference=str(ref),
        cpg=True,
        modified_bases=["5mC"],
        bgzf=True,
        verbose=False,
    )
    assert pileup["status"] == "completed" and pileup["record_count"] > 20
    pysam.tabix_index(str(tmp_path / "cpg.bed.gz"), preset="bed", force=True, keep_original=True)
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t0\t300\tfirst\nchr1\t300\t600\tsecond\n", encoding="utf-8")
    stats = m.modkit_stats(
        str(tmp_path / "cpg.bed.gz"), str(regions), str(tmp_path / "stats.tsv"), verbose=False
    )
    assert stats["status"] == "completed" and stats["row_count"] == 2
    assert stats["columns"][:3] == ["#chrom", "start", "end"]
    motifs = m.modkit_motif_bed(str(ref), "CG", 0, str(tmp_path / "cg.bed"), verbose=False)
    sites = [
        line.split("\t") for line in (tmp_path / "cg.bed").read_text(encoding="utf-8").splitlines()
    ]
    reference = "".join((ref.read_text(encoding="utf-8").splitlines())[1:])
    assert motifs["record_count"] == len(sites)
    assert [int(site[1]) for site in sites if site[5] == "+"] == [
        i for i in range(len(reference) - 1) if reference[i : i + 2] == "CG"
    ]
    summary = m.modkit_modbam_summary(
        str(bam), str(tmp_path / "summary.tsv"), tsv=True, verbose=False
    )
    assert summary["status"] == "completed" and summary["row_count"] >= 3


@requires_modkit
def test_real_validation_error_is_a_failed_status_not_an_exception(tmp_path):
    bam, _ = _toy_modbam(tmp_path)
    result = m.modkit_modbam_check_tags(str(bam), str(tmp_path / "chk"), verbose=False)
    assert result["status"] == "completed" and result["files"]
    bad = m.modkit_pileup(str(bam), str(tmp_path / "bad.bed"), region="chrNOPE:1-10", verbose=False)
    assert bad["status"] in ("completed", "failed")
