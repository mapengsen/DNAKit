"""Tests for dnakit.alignment.minigraph_p11_cli_functions; real runs skip when minigraph / k8 / the scripts are missing."""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.alignment import minigraph_p11_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

have_minigraph = shutil.which("minigraph") is not None
have_scripts = all(shutil.which(n) for n in ("k8", "mgutils.js", "mgutils-es6.js", "gafcall.js"))
requires_minigraph = pytest.mark.skipif(not have_minigraph, reason="minigraph is not installed")
requires_scripts = pytest.mark.skipif(
    not (have_minigraph and have_scripts), reason="k8 / mgutils.js / gafcall.js missing"
)


def _mutate(seq: str, rng: random.Random) -> str:
    out = list(seq)
    for i in range(0, len(out), 40):
        out[i] = rng.choice("ACGT")
    out[3000:3000] = list("".join(rng.choice("ACGT") for _ in range(400)))  # 400 bp insertion
    return "".join(out)


@pytest.fixture(scope="module")
def assemblies(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    rng = random.Random(3)
    root = tmp_path_factory.mktemp("mg")
    ref = "".join(rng.choice("ACGT") for _ in range(8000))
    files = {"root": root}
    for name, seq in (("ref", ref), ("a", _mutate(ref, rng)), ("b", ref[:5000] + ref[5300:])):
        path = root / f"{name}.fa"
        path.write_text(f">{name}\n{seq}\n")
        files[name] = path
    return files


def test_missing_program(tmp_path: Path) -> None:
    fa = tmp_path / "x.fa"
    fa.write_text(">x\nACGT\n")
    with pytest.raises(BackendUnavailableError):
        m.map_sequences_to_graph(
            fa, [fa], tmp_path / "o", minigraph_path=str(tmp_path / "none"), verbose=False
        )


@requires_minigraph
def test_graph_build_map_cover_call(assemblies: dict[str, Path]) -> None:
    root = assemblies["root"]
    graph = m.build_sequence_graph(
        assemblies["ref"],
        [assemblies["a"], assemblies["b"]],
        root / "g",
        min_aln_length=1000,
        min_depth_aln_length=1000,
        threads=1,
        verbose=False,
    )
    assert graph["segments"] >= 2 and graph["links"] >= 1 and graph["version"].startswith("0.2")
    gfa = graph["outputs"]["output"]
    mapped = m.map_sequences_to_graph(
        gfa,
        [assemblies["a"]],
        root / "map",
        preset="asm",
        base_alignment=True,
        threads=1,
        verbose=False,
    )
    assert mapped["records"] >= 1 and mapped["mapped_queries"] == 1
    cov = m.graph_segment_coverage(
        gfa, [assemblies["a"]], root / "cov", preset="asm", threads=1, verbose=False
    )
    assert cov["segments"] == graph["segments"]
    call = m.call_graph_bubbles(
        gfa,
        assemblies["a"],
        root / "call",
        min_aln_length=1000,
        min_depth_aln_length=1000,
        threads=1,
        verbose=False,
    )
    assert call["bubbles"] >= 1


@requires_minigraph
def test_refuses_existing_directory(assemblies: dict[str, Path]) -> None:
    (assemblies["root"] / "exists").mkdir()
    with pytest.raises(ConfigurationError):
        m.map_sequences_to_graph(
            assemblies["ref"], [assemblies["a"]], assemblies["root"] / "exists", verbose=False
        )


@requires_scripts
def test_mgutils_and_gafcall(assemblies: dict[str, Path]) -> None:
    root = assemblies["root"]
    gfa = m.build_sequence_graph(
        assemblies["ref"],
        [assemblies["a"], assemblies["b"]],
        root / "g2",
        min_aln_length=1000,
        min_depth_aln_length=1000,
        threads=1,
        verbose=False,
    )["outputs"]["output"]
    calls = [
        m.call_graph_bubbles(
            gfa,
            assemblies[n],
            root / f"c_{n}",
            min_aln_length=1000,
            min_depth_aln_length=1000,
            threads=1,
            verbose=False,
        )["outputs"]["output"]
        for n in ("ref", "a", "b")
    ]
    merged = m.merge_bubble_calls(calls, root / "merged", samples=["ref", "a", "b"], verbose=False)
    assert merged["bubbles"] >= 1 and merged["sample_count"] == 3
    vcf = m.convert_calls_to_vcf(
        merged["outputs"]["output"], root / "vcf", samples=["ref", "a", "b"], verbose=False
    )
    assert vcf["variants"] == merged["bubbles"]
    named = m.add_sample_names_to_calls(
        merged["outputs"]["output"], ["x", "y", "z"], root / "named", verbose=False
    )
    assert "x\ty\tz" in Path(named["outputs"]["output"]).read_text()
    renamed = m.rename_fasta_sequences(assemblies["a"], "HS", root / "ren", verbose=False)
    assert Path(renamed["outputs"]["output"]).read_text().startswith(">HS#a")
    gaf = m.map_sequences_to_graph(
        gfa,
        [assemblies["a"]],
        root / "gaf",
        preset="asm",
        base_alignment=True,
        vertex_coordinates=True,
        threads=1,
        verbose=False,
    )
    stable = m.convert_to_stable_gaf(gfa, gaf["outputs"]["output"], root / "stable", verbose=False)
    assert stable["records"] >= 1
    extract = m.extract_long_indels(
        stable["outputs"]["output"],
        root / "ext",
        min_indel_length=50,
        sample_name="a",
        verbose=False,
    )
    assert extract["version"]
    sv = m.merge_sv_evidence(
        [extract["outputs"]["output"]],
        root / "sv",
        min_read_count=1,
        min_read_count_per_strand=0,
        verbose=False,
    )
    assert sv["calls"] >= 0


@requires_scripts
def test_small_text_tools(tmp_path: Path) -> None:
    paf = tmp_path / "ins.paf"
    paf.write_text(
        "q1\t800\t0\t800\t+\tchr1\t10000\t100\t900\t780\t800\t60\ttp:A:P\tde:f:0.05\n"
        "q2\t800\t0\t800\t+\tchr1\t10000\t2000\t2800\t800\t800\t60\ttp:A:P\tde:f:0.001\n"
    )
    assert m.blacklist_insert_alignments(paf, tmp_path / "bl", verbose=False)["regions"] == 1
    assert (
        m.blacklist_insert_alignments(paf, tmp_path / "bl2", identical_only=True, verbose=False)[
            "regions"
        ]
        == 1
    )
    bed = tmp_path / "ins.bed"
    bed.write_text("chr1\t100\t300\t200\nchr1\t500\t900\t400\n")
    anno = m.annotate_inserted_sequences(bed, tmp_path / "anno", verbose=False)
    assert anno["annotated"] == 2
    table = m.summarize_insert_annotation(
        anno["outputs"]["output"], tmp_path / "tbl", verbose=False
    )
    assert table["class_count"] == 1
    vcf_a = tmp_path / "a.vcf"
    vcf_a.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-500\n"
    )
    ev = m.evaluate_sv_calls([vcf_a, vcf_a], tmp_path / "ev", verbose=False)
    assert ev["recall"]["errors"] == 0
    view = m.view_sv_calls(vcf_a, tmp_path / "view", verbose=False)
    assert view["records"] == 1
