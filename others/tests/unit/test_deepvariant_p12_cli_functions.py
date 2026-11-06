"""Tests of the DeepVariant facades with stand-in programs (the real programs need a container image)."""

from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import deepvariant_p12_cli_functions as dv

STUB = """#!{python}
import json, pathlib, re, sys
args = sys.argv[1:]
open({record!r}, 'a').write(json.dumps(args) + '\\n')
if args == ['--version']:
    print('DeepVariant version 9.9.9'); sys.exit(0)
for arg in args:
    if arg.startswith('--outfile_base='):
        pathlib.Path(arg.split('=', 1)[1] + '.visual_report.html').write_text('x')
    m = re.match(r'--(output_vcf\\w*|output_gvcf\\w*|outfile|gvcf_outfile|output|candidate_positions|examples)=(.*)', arg)
    if not m: continue
    path = m.group(2)
    s = re.search(r'@(\\d+)', path)
    if s:
        path = path.replace(s.group(0), '-00000-of-%05d' % int(s.group(1)))
    pathlib.Path(path).write_text('stub')
    if m.group(1).startswith('output_vcf') and '--vcf_stats_report=true' in args:
        base = path[:-7] if path.endswith('.vcf.gz') else path[:-4]
        pathlib.Path(base + '.visual_report.html').write_text('x')
if '--metrics_dir=' in ' '.join(args):
    print('100\\t200.0\\t0.5')
"""


@pytest.fixture()
def env(tmp_path: Path) -> dict[str, Path]:
    for name in (
        "ref.fa",
        "reads.bam",
        "reads2.bam",
        "reads3.bam",
        "pan.gbz",
        "truth.vcf",
        "regions.bed",
        "ex.pbtxt",
        "m.json",
    ):
        (tmp_path / name).write_text("x")
    (tmp_path / "bin").mkdir()
    return {"dir": tmp_path, "bin": tmp_path / "bin"}


def _install(env: dict[str, Path], *names: str) -> None:
    record = env["dir"] / "record.jsonl"
    for name in names:
        path = env["bin"] / name
        path.parent.mkdir(exist_ok=True)
        path.write_text(STUB.format(python=sys.executable, record=str(record)))
        path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _call(env: dict[str, Path], index: int = 0) -> list[str]:
    lines = [json.loads(line) for line in (env["dir"] / "record.jsonl").read_text().splitlines()]
    return [line for line in lines if line != ["--version"]][index]


def _p(env: dict[str, Path], name: str) -> str:
    return str(env["dir"] / name)


def test_public_names() -> None:
    import inspect

    assert len(dv.__all__) == 15
    for name in dv.__all__:
        assert inspect.signature(getattr(dv, name)).parameters["verbose"].default is True


def test_missing_program(env: dict[str, Path]) -> None:
    with pytest.raises(BackendUnavailableError):
        dv.deepvariant_run(
            "WGS",
            _p(env, "ref.fa"),
            _p(env, "reads.bam"),
            _p(env, "o.vcf.gz"),
            deepvariant_dir=str(env["bin"]),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        dv.deepvariant_run(
            "WGS",
            _p(env, "ref.fa"),
            _p(env, "reads.bam"),
            _p(env, "o.vcf.gz"),
            deepvariant_dir="/nonexistent",
            verbose=False,
        )


def test_run_phase_and_report(env: dict[str, Path]) -> None:
    _install(env, "run_deepvariant")
    result = dv.deepvariant_run(
        "PACBIO", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "o.vcf.gz"), output_gvcf=_p(env, "o.g.vcf.gz"),
        phase_vcf=True, vcf_stats_report=True, report_title="T", num_shards=4, regions=["chr1:1-100", "r.bed"],
        make_examples_extra_args={"min_mapping_quality": 5, "keep_duplicates": False},
        deepvariant_dir=str(env["bin"]), verbose=False,
    )  # fmt: skip
    call = _call(env)
    for item in ("--model_type=PACBIO", "--phase_vcf=true", "--vcf_stats_report=true", "--num_shards=4",
                 "--regions=chr1:1-100 r.bed", "--report_title=T",
                 "--make_examples_extra_args=min_mapping_quality=5,keep_duplicates=false"):  # fmt: skip
        assert item in call, item
    assert (
        result["output_vcf"]
        and result["output_gvcf"]
        and result["stats_report"].endswith("o.visual_report.html")
    )
    assert result["version"] == "9.9.9"
    json.dumps(result)


def test_run_validation(env: dict[str, Path]) -> None:
    args = ("WGS", _p(env, "ref.fa"), _p(env, "reads.bam"))
    with pytest.raises(ConfigurationError):
        dv.deepvariant_run("XYZ", *args[1:], _p(env, "a.vcf"), verbose=False)
    with pytest.raises(ConfigurationError, match="not allowed"):
        dv.deepvariant_run(
            *args, _p(env, "a.vcf"), make_examples_extra_args={"rm": 1}, verbose=False
        )
    with pytest.raises(ConfigurationError, match="comma"):
        dv.deepvariant_run(
            *args, _p(env, "a.vcf"), make_examples_extra_args={"regions": "a,b"}, verbose=False
        )
    with pytest.raises(ConfigurationError, match="runtime_report"):
        dv.deepvariant_run(*args, _p(env, "a.vcf"), runtime_report=True, verbose=False)
    with pytest.raises(ConfigurationError):
        dv.deepvariant_run(*args, _p(env, "a.vcf"), regions=["a b"], verbose=False)
    (env["dir"] / "exists.vcf").write_text("x")
    with pytest.raises(ConfigurationError, match="must not exist"):
        dv.deepvariant_run(*args, _p(env, "exists.vcf"), verbose=False)


def test_deeptrio_duo_and_trio(env: dict[str, Path]) -> None:
    _install(env, "run_deeptrio")
    common = {"deepvariant_dir": str(env["bin"]), "verbose": False}
    with pytest.raises(ConfigurationError, match="parent"):
        dv.deepvariant_run_deeptrio(
            "WGS", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "c.vcf"), "child", **common
        )
    with pytest.raises(ConfigurationError, match="together"):
        dv.deepvariant_run_deeptrio("WGS", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "c.vcf"), "child",
                                    reads_parent1=_p(env, "reads2.bam"), **common)  # fmt: skip
    result = dv.deepvariant_run_deeptrio(
        "WGS", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "c.vcf"), "child",
        reads_parent1=_p(env, "reads2.bam"), output_vcf_parent1=_p(env, "p1.vcf"), sample_name_parent1="p1",
        reads_parent2=_p(env, "reads3.bam"), output_vcf_parent2=_p(env, "p2.vcf"), sample_name_parent2="p2",
        **common,
    )  # fmt: skip
    call = _call(env)
    assert (
        "--output_vcf_parent2=" + _p(env, "p2.vcf") in call and "--sample_name_child=child" in call
    )
    assert set(result["output_vcfs"]) == {"child", "parent1", "parent2"}


def test_pangenome_and_somatic(env: dict[str, Path]) -> None:
    _install(env, "run_pangenome_aware_deepvariant", "run_deepsomatic")
    common = {"deepvariant_dir": str(env["bin"]), "verbose": False}
    pan = dv.deepvariant_run_pangenome_aware(
        "WGS",
        _p(env, "ref.fa"),
        _p(env, "reads.bam"),
        _p(env, "pan.gbz"),
        _p(env, "pan.vcf"),
        sample_name_reads="s",
        **common,
    )
    assert "--pangenome=" + _p(env, "pan.gbz") in _call(env) and pan["output_vcf"]
    with pytest.raises(ConfigurationError, match="normal"):
        dv.deepvariant_run_deepsomatic(
            "WGS", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "s.vcf"), "tumor", **common
        )
    with pytest.raises(ConfigurationError, match="Tumor-only"):
        dv.deepvariant_run_deepsomatic("WGS_TUMOR_ONLY", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "s.vcf"), "tumor",
                                       reads_normal=_p(env, "reads2.bam"), sample_name_normal="n", **common)  # fmt: skip
    somatic = dv.deepvariant_run_deepsomatic(
        "WGS", _p(env, "ref.fa"), _p(env, "reads.bam"), _p(env, "s.vcf"), "tumor",
        reads_normal=_p(env, "reads2.bam"), sample_name_normal="normal", pon_filtering=_p(env, "truth.vcf"),
        process_somatic=True, **common,
    )  # fmt: skip
    call = _call(env, 1)
    assert (
        "--reads_normal=" + _p(env, "reads2.bam") in call
        and "--process_somatic=true" in call
        and somatic["output_vcf"]
    )


@pytest.mark.parametrize(
    "program,reads,expected",
    [
        ("make_examples", ["reads.bam"], ["--reads={0}"]),
        (
            "make_examples_somatic",
            ["reads.bam", "reads2.bam"],
            ["--reads_tumor={0}", "--reads_normal={1}"],
        ),
        ("multisample_make_examples", ["reads.bam", "reads2.bam"], ["--reads={0},{1}"]),
        ("make_examples_pangenome_aware_dv", ["reads.bam"], ["--reads={0}", "--pangenome={g}"]),
    ],
)
def test_make_examples_programs(
    env: dict[str, Path], program: str, reads: list[str], expected: list[str]
) -> None:
    _install(env, program)
    paths = [_p(env, name) for name in reads]
    result = dv.deepvariant_make_examples(
        program, _p(env, "ref.fa"), paths, examples_path=_p(env, "ex.tfrecord@2.gz"), task=1,
        pangenome_path=_p(env, "pan.gbz") if "pangenome" in program else None,
        extra_flags={"min_base_quality": 10}, deepvariant_dir=str(env["bin"]), verbose=False,
    )  # fmt: skip
    call = _call(env)
    for item in expected:
        assert item.format(*paths, g=_p(env, "pan.gbz")) in call
    assert "--mode=calling" in call and "--task=1" in call and "--min_base_quality=10" in call
    assert result["examples_files"] == [_p(env, "ex.tfrecord-00000-of-00002.gz")]


def test_make_examples_modes_and_trio(env: dict[str, Path]) -> None:
    _install(env, "make_examples")
    (env["bin"] / "deeptrio").mkdir()
    _install(env, "deeptrio/make_examples")
    ref, bam = _p(env, "ref.fa"), _p(env, "reads.bam")
    with pytest.raises(ConfigurationError, match="training"):
        dv.deepvariant_make_examples(
            "make_examples",
            ref,
            [bam],
            examples_path=_p(env, "t.gz"),
            mode="training",
            deepvariant_dir=str(env["bin"]),
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="examples_path"):
        dv.deepvariant_make_examples(
            "make_examples", ref, [bam], deepvariant_dir=str(env["bin"]), verbose=False
        )
    with pytest.raises(ConfigurationError, match="not allowed"):
        dv.deepvariant_make_examples("make_examples", ref, [bam], examples_path=_p(env, "t.gz"), extra_flags={"reads_tumor": "x"}, deepvariant_dir=str(env["bin"]), verbose=False)  # fmt: skip
    sweep = dv.deepvariant_make_examples(
        "make_examples", ref, [bam], mode="candidate_sweep", candidate_positions=_p(env, "cand.tsv"),
        deepvariant_dir=str(env["bin"]), verbose=False,
    )  # fmt: skip
    assert sweep["candidate_positions_files"] == [
        _p(env, "cand.tsv")
    ] and "--mode=candidate_sweep" in _call(env)
    training = dv.deepvariant_make_examples(
        "make_examples", ref, [bam], examples_path=_p(env, "tr.gz"), mode="training", truth_variants=_p(env, "truth.vcf"),
        confident_regions=_p(env, "regions.bed"), deepvariant_dir=str(env["bin"]), verbose=False,
    )  # fmt: skip
    assert training["examples_files"] == [_p(env, "tr.gz")]
    trio = dv.deepvariant_make_examples(
        "deeptrio_make_examples", ref, [bam, _p(env, "reads2.bam")], examples_path=_p(env, "trio.gz"),
        sample_names=["c", "p"], deepvariant_dir=str(env["bin"]), verbose=False,
    )  # fmt: skip
    assert "--reads_parent1=" + _p(env, "reads2.bam") in _call(env, 2) and trio["examples_files"]
    with pytest.raises(ConfigurationError, match="deepvariant_dir"):
        dv.deepvariant_make_examples(
            "deeptrio_make_examples", ref, [bam, bam], examples_path=_p(env, "x.gz"), verbose=False
        )


def test_call_postprocess_report_and_others(env: dict[str, Path]) -> None:
    _install(env, "call_variants", "postprocess_variants", "vcf_stats_report", "show_examples", "train",
             "labeled_examples_to_vcf", "compare_labelers", "preprocess_truth", "print_f1", "run_oracle_inference")  # fmt: skip
    kw = {"deepvariant_dir": str(env["bin"]), "verbose": False}
    (env["dir"] / "ex-00000-of-00002.gz").write_text("x")
    called = dv.deepvariant_call_variants(
        _p(env, "ex@2.gz"), _p(env, "cvo.gz"), "/model/ckpt", batch_size=64, **kw
    )
    assert called["input_example_files"] == [
        _p(env, "ex-00000-of-00002.gz")
    ] and "--batch_size=64" in _call(env)
    with pytest.raises(ConfigurationError, match="No example"):
        dv.deepvariant_call_variants(_p(env, "nothing@3.gz"), _p(env, "x.gz"), "/m", **kw)
    with pytest.raises(ConfigurationError, match="nonvariant"):
        dv.deepvariant_postprocess_variants(
            _p(env, "ref.fa"),
            _p(env, "cvo.gz"),
            _p(env, "p.vcf"),
            output_gvcf=_p(env, "p.g.vcf"),
            **kw,
        )
    post = dv.deepvariant_postprocess_variants(
        _p(env, "ref.fa"), _p(env, "cvo.gz"), _p(env, "p.vcf"), output_gvcf=_p(env, "p.g.vcf"),
        nonvariant_site_tfrecord_path=_p(env, "ex-00000-of-00002.gz"), cpus=2, **kw,
    )  # fmt: skip
    assert post["output_gvcf"] and "--gvcf_outfile=" + _p(env, "p.g.vcf") in _call(env, 1)
    (env["dir"] / "in.vcf").write_text("x")
    dv.deepvariant_vcf_stats_report(_p(env, "in.vcf"), _p(env, "rep"), title="X", **kw)
    assert "--outfile_base=" + _p(env, "rep") in _call(env, 2)
