"""Tests for dnakit.sequence.gimmemotifs_precondition_functions (Fig. 2 coverage gap, batch F).

Offline by design: inputs are small FASTA, PFM, table and GTF files written by the tests, and the
gffread and orthofinder programs are replaced by small scripts placed on PATH. Nothing touches the
network. The tests import the module directly, so they run before the area package is wired.
"""

from __future__ import annotations

import importlib.util
import inspect
import os
import random
import stat
import sys
from pathlib import Path

import pytest

pytest.importorskip("gimmemotifs")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.sequence import gimmemotifs_precondition_functions as pre  # noqa: E402

PFM_TEXT = """>MA0001.1
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
>MA0002.1
0.1 0.1 0.7 0.1
0.1 0.7 0.1 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
0.4 0.2 0.2 0.2
0.2 0.4 0.2 0.2
"""

GTF_TEXT = (
    'chr1\tsrc\ttranscript\t100\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; gene_name "GATA1";\n'
    'chr1\tsrc\texon\t100\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; gene_name "GATA1";\n'
    'chr2\tsrc\ttranscript\t500\t1500\t.\t+\t.\tgene_id "ENSG2"; transcript_id "ENST2"; gene_name "CTCF";\n'
)

GFFREAD_SCRIPT = """#!__PYTHON__
import re, sys
args = sys.argv[1:]
pep = args[args.index("-y") + 1]
gtf = [a for a in args if a.endswith(".gtf")][0]
out = []
for line in open(gtf):
    fields = line.rstrip("\\n").split("\\t")
    if len(fields) < 9 or fields[2] != "transcript":
        continue
    name = re.search('gene_name "([^"]+)"', fields[8]).group(1)
    gene = re.search('gene_id "([^"]+)"', fields[8]).group(1)
    tid = re.search('transcript_id "([^"]+)"', fields[8]).group(1)
    out.append(">" + tid + " " + name + " " + gene + "\\nMKTAYIAKQRQISFVKSHFSRQLE\\n")
open(pep, "w").write("".join(out))
"""

ORTHOFINDER_SCRIPT = """#!__PYTHON__
import sys
from pathlib import Path
args = sys.argv[1:]
folder = Path(args[args.index("-f") + 1])
results = folder.parent / "Results_fake"
(results / "Orthogroups").mkdir(parents=True, exist_ok=True)
header = ["Orthogroup"]
cells = {}
for fasta in sorted(folder.glob("*.pep.fa")):
    asm = fasta.name[: -len(".pep.fa")]
    header.append(asm + ".pep")
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name, gene_id = line[1:].split("|")
            cells.setdefault(gene_id, {})[asm] = name + "|" + gene_id
lines = ["\\t".join(header)]
for index, (gene_id, per_asm) in enumerate(sorted(cells.items())):
    values = [per_asm.get(column[: -len(".pep")], "") for column in header[1:]]
    lines.append("\\t".join(["OG%07d" % index] + values))
(results / "Orthogroups" / "Orthogroups.tsv").write_text("\\n".join(lines) + "\\n")
(results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text("\\t".join(header) + "\\n")
print("Results directory:  " + str(results) + "/")
"""


def _seqs(seed: int, count: int, length: int) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(count)]


def _fasta(path: Path, sequences: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">seq{i}\n{seq}\n" for i, seq in enumerate(sequences)), encoding="utf-8")
    return path


def _genome(path: Path, length: int = 300_000) -> Path:
    rng = random.Random(3)
    body = "".join(rng.choice("ACGT") for _ in range(length))
    path.write_text(f">chr1\n{body}\n>chr2\n{body[::-1]}\n", encoding="utf-8")
    return path


def _program(bin_dir: Path, name: str, template: str) -> None:
    bin_dir.mkdir(parents=True, exist_ok=True)
    script = bin_dir / name
    script.write_text(template.replace("__PYTHON__", sys.executable), encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)


@pytest.fixture()
def motifs(tmp_path: Path) -> Path:
    path = tmp_path / "motifs.pfm"
    path.write_text(PFM_TEXT, encoding="utf-8")
    return path


@pytest.fixture()
def fasta_inputs(tmp_path: Path) -> dict[str, Path]:
    foreground = _fasta(tmp_path / "fg.fa", _seqs(7, 30, 60))
    background = _fasta(tmp_path / "bg.fa", _seqs(8, 30, 60))
    location = _fasta(tmp_path / "loc.fa", _seqs(9, 20, 200))
    return {"fg": foreground, "bg": background, "loc": location}


@pytest.fixture()
def fake_programs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bin_dir = tmp_path / "bin"
    _program(bin_dir, "gffread", GFFREAD_SCRIPT)
    _program(bin_dir, "orthofinder", ORTHOFINDER_SCRIPT)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return bin_dir


def _annotation(tmp_path: Path) -> tuple[Path, Path]:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF_TEXT, encoding="utf-8")
    return _genome(tmp_path / "genome.fa", 20_000), gtf


def _orthogroup_results(root: Path) -> Path:
    groups = root / "Results" / "Orthogroups"
    groups.mkdir(parents=True)
    (groups / "Orthogroups.tsv").write_text(
        "Orthogroup\tasmA.pep\tasmB.pep\n"
        "OG0000000\tGATA1|ENSG1\tGATA1|ENSG1\n"
        "OG0000001\tCTCF|ENSG2\t\n",
        encoding="utf-8",
    )
    (groups / "Orthogroups_UnassignedGenes.tsv").write_text(
        "Orthogroup\tasmA.pep\tasmB.pep\n", encoding="utf-8"
    )
    return root / "Results"


def test_public_functions_are_keyword_only_verbose() -> None:
    names = pre.__all__
    assert names, "the module must export its public functions"
    for name in names:
        function = getattr(pre, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_diff_plot_draws_png_and_returns_tables(tmp_path: Path, motifs: Path) -> None:
    pytest.importorskip("matplotlib")
    output = tmp_path / "diff.png"
    result = pre.gimmemotifs_diff_plot(
        ["MA0001.1", "MA0002.1"],
        str(motifs),
        ["first", "second"],
        [[0.5, 0.1], [0.4, 0.05]],
        [[10.0, 10.0], [8.0, 8.0]],
        [0.05, 0.05],
        [1.0, 1.0],
        str(output),
        verbose=False,
    )
    assert result["plot_written"] is True
    assert output.is_file() and output.stat().st_size > 0
    assert result["frequencies"] == [[0.5, 0.1], [0.4, 0.05]]


def test_diff_plot_rejects_wrong_table_shape(tmp_path: Path, motifs: Path) -> None:
    with pytest.raises(ConfigurationError):
        pre.gimmemotifs_diff_plot(
            ["MA0001.1", "MA0002.1"],
            str(motifs),
            ["sample"],
            [[0.5]],
            [[10.0]],
            [0.05, 0.05],
            [1.0, 1.0],
            str(tmp_path / "x.png"),
            verbose=False,
        )


def test_diff_motifs_counts_background_and_foreground(tmp_path: Path, motifs: Path, fasta_inputs) -> None:
    pytest.importorskip("matplotlib")
    second = _fasta(tmp_path / "fg2.fa", _seqs(10, 30, 60))
    result = pre.gimmemotifs_diff_motifs(
        [str(fasta_inputs["fg"]), str(second)],
        str(fasta_inputs["bg"]),
        str(motifs),
        str(tmp_path / "diff.png"),
        min_enrichment=0.0,
        min_frequency=0.0,
        verbose=False,
    )
    assert result["motif_ids"] == ["MA0001.1", "MA0002.1"]
    assert result["sample_names"] == ["fg", "fg2"]
    assert len(result["frequencies"]) == 2 and len(result["frequencies"][0]) == 2
    assert result["background_size"] == 30
    assert isinstance(result["plot_written"], bool)


def test_select_nonredundant_motifs_uses_report_and_tables(tmp_path: Path, motifs: Path) -> None:
    pytest.importorskip("sklearn")
    report = tmp_path / "roc.tsv"
    report.write_text(
        "\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
        "MA0001.1\t0.80\t0.60\t3.5\t0.40\n"
        "MA0002.1\t0.70\t0.50\t3.0\t0.30\n",
        encoding="utf-8",
    )
    rng = random.Random(1)
    for name in ("fg_table", "bg_table"):
        rows = ["seq\tMA0001.1\tMA0002.1"]
        rows += [f"s{i}\t{rng.random():.4f}\t{rng.random():.4f}" for i in range(30)]
        (tmp_path / f"{name}.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    result = pre.gimmemotifs_select_nonredundant_motifs(
        str(report),
        str(motifs),
        str(tmp_path / "fg_table.tsv"),
        str(tmp_path / "bg_table.tsv"),
        verbose=False,
    )
    assert result["selected_motifs"]
    assert set(result["selected_motifs"]) <= {"MA0001.1", "MA0002.1"}


def test_create_background_random_and_custom(tmp_path: Path, fasta_inputs) -> None:
    random_result = pre.gimmemotifs_denovo_create_background(
        "random",
        str(fasta_inputs["fg"]),
        str(tmp_path / "out" / "random.fa"),
        size=60,
        nr_times=2,
        verbose=False,
    )
    assert random_result["sequence_count"] == 60
    assert (tmp_path / "out" / "random.fa").read_text().count(">") == 60

    custom_result = pre.gimmemotifs_denovo_create_background(
        "custom",
        str(fasta_inputs["fg"]),
        str(tmp_path / "custom.fa"),
        custom_background=str(fasta_inputs["bg"]),
        size=60,
        verbose=False,
    )
    assert custom_result["sequence_count"] == 30


def test_create_background_needs_a_genome_for_genomic(tmp_path: Path, fasta_inputs) -> None:
    with pytest.raises(ConfigurationError):
        pre.gimmemotifs_denovo_create_background(
            "genomic", str(fasta_inputs["fg"]), str(tmp_path / "g.fa"), verbose=False
        )


def test_create_background_genomic_from_fasta_genome(tmp_path: Path, fasta_inputs) -> None:
    genome = _genome(tmp_path / "genome.fa")
    result = pre.gimmemotifs_denovo_create_background(
        "genomic",
        str(fasta_inputs["fg"]),
        str(tmp_path / "genomic.fa"),
        genome_fasta=str(genome),
        size=60,
        nr_times=2,
        verbose=False,
    )
    assert result["sequence_count"] > 0
    assert (tmp_path / "genomic.fa").is_file()


def test_create_backgrounds_random_layout(tmp_path: Path, fasta_inputs) -> None:
    result = pre.gimmemotifs_denovo_create_backgrounds(
        str(tmp_path / "bgs"),
        prediction_fasta=str(fasta_inputs["fg"]),
        validation_fasta=str(fasta_inputs["bg"]),
        background_types=["random"],
        size=60,
        verbose=False,
    )
    assert (tmp_path / "bgs" / "prediction.bg.fa").is_file()
    assert (tmp_path / "bgs" / "bg.random.fa").is_file()
    assert result["background_types"] == ["random"]


def test_filter_significant_motifs_applies_minimum_values(tmp_path: Path, motifs: Path) -> None:
    stats = {
        "MA0001.1": {"random": {"max_enrichment": 5.0, "roc_auc": 0.9, "enr_at_fpr": 0.8}},
        "MA0002.1": {"random": {"max_enrichment": 1.0, "roc_auc": 0.4, "enr_at_fpr": 0.2}},
    }
    result = pre.gimmemotifs_filter_significant_motifs(
        str(motifs),
        str(tmp_path / "significant.pfm"),
        background_type="random",
        stats=stats,
        verbose=False,
    )
    assert result["significant_motifs"] == ["MA0001.1"]
    assert "MA0001.1" in (tmp_path / "significant.pfm").read_text()


def test_best_motif_in_cluster_returns_one_motif_per_cluster(tmp_path: Path, motifs: Path, fasta_inputs) -> None:
    result = pre.gimmemotifs_best_motif_in_cluster(
        str(motifs),
        str(motifs),
        [
            {"representative": "MA0001.1", "members": ["MA0001.1", "MA0002.1"]},
            {"representative": "MA0002.1", "members": ["MA0002.1"]},
        ],
        str(fasta_inputs["fg"]),
        {"random": str(fasta_inputs["bg"])},
        str(_genome(tmp_path / "genome.fa", 20_000)),
        verbose=False,
    )
    assert len(result["best_motifs"]) == 2
    assert set(result["best_motifs"]) <= {"MA0001.1", "MA0002.1"}


def test_best_motif_in_cluster_rejects_unknown_motif(tmp_path: Path, motifs: Path, fasta_inputs) -> None:
    with pytest.raises(ConfigurationError):
        pre.gimmemotifs_best_motif_in_cluster(
            str(motifs),
            str(motifs),
            [{"representative": "NOPE", "members": ["MA0001.1"]}],
            str(fasta_inputs["fg"]),
            {"random": str(fasta_inputs["bg"])},
            str(_genome(tmp_path / "genome.fa", 20_000)),
            verbose=False,
        )


def test_denovo_motif_report_writes_html(tmp_path: Path, motifs: Path, fasta_inputs) -> None:
    result = pre.gimmemotifs_denovo_motif_report(
        str(fasta_inputs["fg"]),
        str(motifs),
        str(fasta_inputs["fg"]),
        {"random": str(fasta_inputs["bg"])},
        str(fasta_inputs["loc"]),
        str(tmp_path / "report"),
        genome_fasta=str(_genome(tmp_path / "genome.fa", 20_000)),
        verbose=False,
    )
    assert any(name.endswith(".html") for name in result["outputs"])


def test_roc_plots_writes_one_png_per_motif_and_background(tmp_path: Path, motifs: Path, fasta_inputs) -> None:
    pytest.importorskip("matplotlib")
    result = pre.gimmemotifs_roc_plots(
        str(motifs),
        str(fasta_inputs["fg"]),
        {"random": str(fasta_inputs["bg"])},
        str(tmp_path / "roc"),
        genome_fasta=str(_genome(tmp_path / "genome.fa", 20_000)),
        verbose=False,
    )
    assert result["motif_count"] == 2
    images = sorted((tmp_path / "roc" / "images").glob("*.png"))
    assert len(images) == 2


def test_roc_html_report_writes_named_file(tmp_path: Path, motifs: Path) -> None:
    table = tmp_path / "enrichment.txt"
    table.write_text(
        "\tP-value\t% matches input\t% matches background\t# matches\tROC AUC\tPR AUC\t"
        "Enr. at 1% FPR\tRecall at 10% FDR\n"
        "MA0001.1\t0.0001\t40.0\t10.0\t12\t0.80\t0.60\t3.5\t0.40\n"
        "MA0002.1\t0.02\t30.0\t15.0\t9\t0.70\t0.50\t2.5\t0.30\n",
        encoding="utf-8",
    )
    result = pre.gimmemotifs_roc_html_report(
        str(tmp_path / "html"), str(table), str(motifs), "report.html", verbose=False
    )
    assert (tmp_path / "html" / "report.html").is_file()
    assert result["outputs"]["report"].endswith("report.html")
    with pytest.raises(ConfigurationError):
        pre.gimmemotifs_roc_html_report(
            str(tmp_path / "html2"), str(table), str(motifs), "../escape.html", verbose=False
        )


def test_motif_tree_plot_needs_ete3(tmp_path: Path) -> None:
    if importlib.util.find_spec("ete3") is not None:
        pytest.skip("ete3 is installed; the tree rendering is not exercised by this offline test")
    with pytest.raises(BackendUnavailableError):
        pre.gimmemotifs_motif_tree_plot(
            str(tmp_path / "tree.png"),
            "((MA0001.1:1,MA0002.1:1):1,X:1);",
            {"MA0001.1": [1.0, 2.0], "MA0002.1": [0.5, 1.0], "X": [0.1, 0.2]},
            verbose=False,
        )


def test_annot2primpep_uses_gffread_output(tmp_path: Path, fake_programs: Path) -> None:
    genome, gtf = _annotation(tmp_path)
    result = pre.gimmemotifs_annot2primpep(
        str(genome), str(gtf), str(tmp_path / "annot"), name="asmA", verbose=False
    )
    peptides = Path(result["outputs"]["peptides"])
    assert peptides.name == "asmA.pep.fa"
    assert result["protein_count"] == 2
    assert ">GATA1|ENSG1" in peptides.read_text()


def test_annot2primpep_missing_gffread_is_backend_unavailable(tmp_path: Path, monkeypatch) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    genome, gtf = _annotation(tmp_path)
    with pytest.raises(BackendUnavailableError) as caught:
        pre.gimmemotifs_annot2primpep(str(genome), str(gtf), str(tmp_path / "annot"), verbose=False)
    assert "gffread" in str(caught.value)


def test_orthogroup_database_lookup_and_motif2factors(tmp_path: Path, motifs: Path) -> None:
    results = _orthogroup_results(tmp_path)
    database = tmp_path / "orthologs.sqlite"
    loaded = pre.gimmemotifs_load_orthogroups_in_db(
        str(database), ["asmA", "asmB"], str(results), verbose=False
    )
    assert loaded["orthogroup_count"] == 3
    assert loaded["gene_counts"] == {"asmA": 2, "asmB": 1}

    orthogroups = pre.gimmemotifs_factor2orthogroups(
        "GATA1", ["asmA"], str(database), verbose=False
    )
    assert orthogroups["orthogroups"] == ["OG0000000"]

    table = pre.gimmemotifs_make_motif2factors(
        str(tmp_path / "out" / "asmB"),
        motif_file=str(motifs),
        motif_factors={"MA0002.1": ["GATA1"], "MA0001.1": ["NOPE"]},
        database_path=str(database),
        new_reference="asmB",
        database_references=["asmA"],
        verbose=False,
    )
    text = Path(table["outputs"]["motif2factors"]).read_text(encoding="utf-8")
    assert "MA0002.1\tGATA1\tOrthologs\tN" in text
    assert "MA0001.1\tNO ORTHOLOGS FOUND" in text


def test_offline_strategy_only(tmp_path: Path) -> None:
    results = _orthogroup_results(tmp_path)
    database = tmp_path / "orthologs.sqlite"
    pre.gimmemotifs_load_orthogroups_in_db(str(database), ["asmA", "asmB"], str(results), verbose=False)
    with pytest.raises(ConfigurationError):
        pre.gimmemotifs_factor2orthogroups("GATA1", ["asmA"], str(database), strategy="lenient", verbose=False)


def test_motif2factor_from_orthologs_runs_pipeline_with_local_inputs(
    tmp_path: Path, motifs: Path, fake_programs: Path
) -> None:
    genome, gtf = _annotation(tmp_path)
    genomes = {
        "asmA": {"fasta": str(genome), "annotation": str(gtf)},
        "asmB": {"fasta": str(genome), "annotation": str(gtf)},
    }
    result = pre.gimmemotifs_motif2factor_from_orthologs(
        genomes,
        ["asmA"],
        ["asmB"],
        str(motifs),
        {"MA0002.1": ["GATA1"]},
        str(tmp_path / "orth"),
        verbose=False,
    )
    table = Path(result["outputs"]["asmB_motif2factors"])
    assert "MA0002.1\tGATA1\tOrthologs\tN" in table.read_text(encoding="utf-8")


def test_motif2factor_from_orthologs_missing_orthofinder(tmp_path: Path, motifs: Path, monkeypatch) -> None:
    genome, gtf = _annotation(tmp_path)
    bin_dir = tmp_path / "gffonly"
    _program(bin_dir, "gffread", GFFREAD_SCRIPT)
    monkeypatch.setenv("PATH", str(bin_dir))
    with pytest.raises(BackendUnavailableError) as caught:
        pre.gimmemotifs_motif2factor_from_orthologs(
            {"asmA": {"fasta": str(genome), "annotation": str(gtf)}},
            ["asmA"],
            ["asmA"],
            str(motifs),
            {},
            str(tmp_path / "orth"),
            verbose=False,
        )
    assert "orthofinder" in str(caught.value)
