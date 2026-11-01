"""Tests for dnakit.variants.gemini_p10_functions (GEMINI sub-commands as typed facades).

A stand-in ``gemini`` executable written into tmp_path records its argument vector and prints a small table (or creates the
database for ``merge_chunks``), so the argument mapping, the validation and the result shape are checked without GEMINI.
A real GEMINI needs about 20 GB of annotation files to load a database, so no test starts the real program.
"""

from __future__ import annotations

import inspect
import json
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import gemini_p10_functions as g

_STUB = r'''#!{python}
import json, os, sys
argv = sys.argv[1:]
here = os.path.dirname(os.path.abspath(__file__))
if argv == ["--version"]:
    print("gemini 0.30.2")
    sys.exit(0)
with open(os.path.join(here, "calls.jsonl"), "a") as handle:
    handle.write(json.dumps(argv) + "\n")
if os.path.exists(os.path.join(here, "FAIL")):
    sys.stderr.write("simulated failure\n")
    sys.exit(1)
if "merge_chunks" in argv:
    open(argv[argv.index("--db") + 1], "w").write("merged")
elif not any(command in argv for command in ("amend", "annotate", "set_somatic")):
    sys.stdout.write("chrom\tstart\tend\nchr1\t10\t11\nchr1\t20\t21\n")
'''


@pytest.fixture
def fake(tmp_path: Path) -> Path:
    directory = tmp_path / "bin"
    directory.mkdir()
    program = directory / "gemini"
    program.write_text(_STUB.format(python=sys.executable), encoding="utf-8")
    program.chmod(program.stat().st_mode | stat.S_IXUSR)
    return program


def _calls(program: Path) -> list[list[str]]:
    return [json.loads(line) for line in (program.parent / "calls.jsonl").read_text().splitlines()]


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, str]:
    out = {}
    for name in ("variants.db", "chunk1.db", "chunk2.db", "samples.ped", "anno.bed.gz", "edges.txt", "calls.vcf"):
        path = tmp_path / name
        path.write_text("x", encoding="utf-8")
        out[name] = str(path)
    return out


def _cases(i: dict[str, str], out: Path):
    db = i["variants.db"]
    return [
        ("gemini_merge_chunks", ([i["chunk1.db"], i["chunk2.db"]], str(out / "m.db")), dict(index=True, anno_type="VEP"),
         ["merge_chunks", "--chunkdb", i["chunk1.db"], i["chunk2.db"], "--db", str(out / "m.db"), "-t", "VEP", "--index"]),
        ("gemini_amend", (db, i["samples.ped"]), dict(output_database=str(out / "a.db"), clear=True),
         ["amend", str(out / "a.db"), "--sample", i["samples.ped"], "--clear"]),
        ("gemini_annotate", (db, i["anno.bed.gz"], "score"), dict(output_database=str(out / "n.db"), annotation_type="count"),
         ["annotate", str(out / "n.db"), "-f", i["anno.bed.gz"], "-c", "score", "-a", "count"]),
        ("gemini_set_somatic", (db,), dict(output_database=str(out / "s.db"), min_depth=10, min_tumor_alt_count=3, dry_run=True),
         ["set_somatic", str(out / "s.db"), "--min-depth", "10.0", "--min-tumor-alt-count", "3", "--dry-run"]),
        ("gemini_carrier_summary_by_phenotype", (db, "phenotype", str(out / "c.tsv")), dict(),
         ["query", "--header", "--carrier-summary-by-phenotype", "phenotype", "-q", "SELECT chrom, start, end, ref, alt, gene FROM variants", db]),
        ("gemini_dump", (db, str(out / "d.tsv")), dict(table="genotypes", header=True),
         ["dump", "--genotypes", "--header", db]),
        ("gemini_region", (db, str(out / "r.tsv")), dict(gene="TP53", columns="chrom,start"),
         ["region", "--gene", "TP53", "--columns", "chrom,start", db]),
        ("gemini_stats", (db, str(out / "st.tsv")), dict(statistic="sfs"), ["stats", "--sfs", db]),
        ("gemini_windower", (db, str(out / "w.tsv")), dict(analysis="hwe", window_size=5000), ["windower", "-t", "hwe", "-w", "5000", "-s", "0", "-o", "mean", db]),
        ("gemini_inheritance_model", (db, str(out / "im.tsv")), dict(model="x_linked_recessive", x_chromosomes=["X"], min_gq=20),
         ["x_linked_recessive", "--min-gq", "20", "-X", "X", db]),
        ("gemini_pathways", (db, str(out / "p.tsv")), dict(lof=True), ["pathways", "-v", "68", "--lof", db]),
        ("gemini_lof_sieve", (db, str(out / "l.tsv")), dict(), ["lof_sieve", db]),
        ("gemini_burden", (db, str(out / "b.tsv")), dict(nonsynonymous=True, cases=["s1", "s2"], controls=["s3"]),
         ["burden", "--nonsynonymous", "--cases", "s1", "s2", "--controls", "s3", db]),
        ("gemini_calpha_test", (db, str(out / "ca.tsv"), ["s1"], ["s2"]), dict(permutations=100),
         ["burden", "--calpha", "--cases", "s1", "--controls", "s2", "--permutations", "100", db]),
        ("gemini_interactions", (db, str(out / "in.tsv")), dict(gene="TP53", radius=2),
         ["interactions", "-g", "TP53", "-r", "2", db]),
        ("gemini_fusions", (db, str(out / "f.tsv")), dict(evidence_type="PE,SR"), ["fusions", "--evidence_type", "PE,SR", db]),
        ("gemini_gene_wise", (db, str(out / "gw.tsv")), dict(genotype_filters=["gt_types.s1 == HET"], min_filters=1),
         ["gene_wise", "--min-filters", "1", "--gt-filter", "gt_types.s1 == HET", db]),
        ("gemini_sex_check", (db, str(out / "sx.tsv")), dict(chromosome="X"), ["qc", "--mode", "sex", "--chrom", "X", db]),
    ]


def _contains(argv: list[str], segment: list[str]) -> bool:
    return any(argv[k : k + len(segment)] == segment for k in range(len(argv) - len(segment) + 1))


def test_every_public_function_has_a_case(fake, inputs, tmp_path):
    assert sorted(name for name, *_ in _cases(inputs, tmp_path)) == sorted(g.__all__)


def test_each_function_maps_its_arguments_to_gemini_flags(fake, inputs, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    for index, (name, args, kwargs, expected) in enumerate(_cases(inputs, out)):
        result = getattr(g, name)(*args, gemini_path=str(fake), verbose=False, **kwargs)
        assert result["status"] == "completed", (name, result["log_tail"])
        assert result["operation"] == name and result["backend"] == "gemini"
        assert result["provenance"]["tested_version"] == "0.30.2"
        json.dumps(result)
        argv = _calls(fake)[index]
        assert argv[: len(expected[:1])] == expected[:1] or _contains(argv, expected[:1]), (name, argv)
        # the sub-command name and every flag group appear in order of the expected list
        cursor = 0
        for token in expected:
            assert token in argv[cursor:], (name, token, argv)
            cursor = argv.index(token, cursor)


def test_table_output_is_streamed_to_the_file_and_previewed(fake, inputs, tmp_path):
    result = g.gemini_region(inputs["variants.db"], str(tmp_path / "t.tsv"), region="chr1:1-50", header=True,
                             gemini_path=str(fake), verbose=False)
    assert Path(result["table"]).read_text().splitlines()[1] == "chr1\t10\t11"
    assert result["columns"] == ["chrom", "start", "end"] and result["row_count"] == 2


def test_mutating_commands_work_on_a_copy_unless_in_place(fake, inputs, tmp_path):
    original = Path(inputs["variants.db"])
    copy = tmp_path / "copy.db"
    result = g.gemini_amend(str(original), inputs["samples.ped"], output_database=str(copy), gemini_path=str(fake), verbose=False)
    assert result["changed_copy"] is True and copy.read_text() == "x" and _calls(fake)[-1][1] == str(copy)
    in_place = g.gemini_set_somatic(str(original), in_place=True, gemini_path=str(fake), verbose=False)
    assert in_place["changed_copy"] is False and _calls(fake)[-1][1] == str(original)
    with pytest.raises(ConfigurationError, match="exactly one"):
        g.gemini_amend(str(original), inputs["samples.ped"], gemini_path=str(fake), verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        g.gemini_amend(str(original), inputs["samples.ped"], output_database=str(tmp_path / "x.db"), in_place=True,
                       gemini_path=str(fake), verbose=False)


def test_bad_arguments_are_rejected_before_gemini_runs(fake, inputs, tmp_path):
    kw = dict(gemini_path=str(fake), verbose=False)
    db = inputs["variants.db"]
    with pytest.raises(ConfigurationError, match="exactly one of gene or region"):
        g.gemini_region(db, str(tmp_path / "a.tsv"), **kw)
    with pytest.raises(ConfigurationError, match="must be one of"):
        g.gemini_dump(db, str(tmp_path / "b.tsv"), table="everything", **kw)
    with pytest.raises(ConfigurationError, match="needs a query"):
        g.gemini_stats(db, str(tmp_path / "c.tsv"), statistic="summarize", **kw)
    with pytest.raises(ConfigurationError, match="autosomal"):
        g.gemini_inheritance_model(db, str(tmp_path / "d.tsv"), model="x_linked_dominant", lenient=True, **kw)
    with pytest.raises(ConfigurationError, match="X-linked"):
        g.gemini_inheritance_model(db, str(tmp_path / "e.tsv"), x_chromosomes=["X"], **kw)
    with pytest.raises(ConfigurationError, match="does not apply"):
        g.gemini_interactions(db, str(tmp_path / "f.tsv"), lof=True, gene="TP53", **kw)
    with pytest.raises(ConfigurationError, match="existing file"):
        g.gemini_lof_sieve(str(tmp_path / "missing.db"), str(tmp_path / "g.tsv"), **kw)
    with pytest.raises(ConfigurationError, match="NUL"):
        g.gemini_gene_wise(db, str(tmp_path / "h.tsv"), where="a\x00b", **kw)
    assert not (fake.parent / "calls.jsonl").exists()
    existing = tmp_path / "keep.tsv"
    existing.write_text("keep")
    with pytest.raises(ConfigurationError, match="must not exist"):
        g.gemini_lof_sieve(db, str(existing), **kw)


def test_failed_run_is_reported(fake, inputs, tmp_path):
    (fake.parent / "FAIL").write_text("")
    result = g.gemini_lof_sieve(inputs["variants.db"], str(tmp_path / "x.tsv"), gemini_path=str(fake), verbose=False)
    assert result["status"] == "failed" and "simulated failure" in result["log_tail"] and "row_count" not in result


def test_missing_program_raises_before_creating_output(inputs, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setattr(sys, "executable", str(tmp_path / "nowhere" / "python"))
    target = tmp_path / "x.tsv"
    with pytest.raises(BackendUnavailableError) as error:
        g.gemini_lof_sieve(inputs["variants.db"], str(target), verbose=False)
    assert "gemini" in str(error.value) and not target.exists()


def test_public_functions_carry_the_verbosity_parameter():
    for name in g.__all__:
        parameter = inspect.signature(getattr(g, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
