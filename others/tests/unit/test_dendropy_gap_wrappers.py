"""Tests for dnakit.phylogeny.dendropy_gap_wrappers (dendropy formats, sumtrees, sumlabels, Entrez, MUSCLE, PAUP*)."""

from __future__ import annotations

import io
import os
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("dendropy")

from dendropy.interop import entrez  # noqa: E402

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.phylogeny import dendropy_gap_wrappers as wrappers  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
MSA = DATA / "msa.fa"
_VENV_BIN = Path(sys.executable).parent
_HAS_DENDROPY_FORMAT = (_VENV_BIN / "dendropy-format").is_file() and (_VENV_BIN / "sumtrees").is_file()
_HAS_MUSCLE = shutil.which("muscle") is not None
_HAS_PAUP = shutil.which("paup") is not None

_FASTA = ">A\nACGTACGT\n>B\nACGTACGA\n>C\nACGAACGT\n"
_PHYLIP = " 3 8\nA ACGTACGT\nB ACGTACGA\nC ACGAACGT\n"
_TREES = """#NEXUS
begin trees;
tree t1 = [&U] ((A:1,B:1)0.9:1,(C:1,D:1)0.8:1);
tree t2 = [&U] ((A:1,B:1)0.95:1,(C:1,D:1)0.7:1);
tree t3 = [&U] (((A:1,B:1):1,C:1)0.6:1,D:1);
end;
"""
_TARGET = "(((A,B),C),D);\n"


@pytest.fixture(autouse=True)
def _venv_programs_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_VENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "aln.fa"
    path.write_text(_FASTA)
    return path


def _write(path: Path, text: str) -> str:
    path.write_text(text)
    return str(path)


# --- format conversion -------------------------------------------------------------------------


def test_convert_fasta_to_nexus_and_phylip(tmp_path: Path, fasta: Path) -> None:
    nexus = wrappers.convert_format(str(fasta), str(tmp_path / "out.nex"), input_format="fasta", data_type="dna",
                                    output_format="nexus", verbose=False)
    assert nexus["status"] == "completed" and nexus["character_matrices"] == 1
    assert (tmp_path / "out.nex").read_text().startswith("#NEXUS")
    phylip = wrappers.format_to_phylip(str(fasta), str(tmp_path / "out.phy"), input_format="fasta", data_type="dna",
                                       verbose=False)
    assert phylip["output_format"] == "phylip"
    assert (tmp_path / "out.phy").read_text().split()[:2] == ["3", "8"]


def test_shortcuts_set_the_output_format(tmp_path: Path, fasta: Path) -> None:
    source = str(fasta)
    assert wrappers.format_to_fasta(source, str(tmp_path / "f.fa"), input_format="fasta", data_type="dna",
                                    verbose=False)["output_format"] == "fasta"
    assert wrappers.format_to_nexml(source, str(tmp_path / "f.xml"), input_format="fasta", data_type="dna",
                                    verbose=False)["output_format"] == "nexml"
    assert (tmp_path / "f.xml").read_text().lstrip().startswith("<?xml")
    newick_source = _write(tmp_path / "t.tre", _TREES)
    newick = wrappers.format_to_newick(newick_source, str(tmp_path / "t.nwk"), input_format="nexus",
                                       unquoted_underscores=True, verbose=False)
    assert newick["output_format"] == "newick"
    assert (tmp_path / "t.nwk").read_text().startswith("[&U] ((A:1.0")


def test_phylip_input_variants_are_read(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.phy", _PHYLIP)
    result = wrappers.convert_format(source, str(tmp_path / "out.fa"), input_format="phylip-relaxed-singlespace",
                                     data_type="dna", output_format="fasta", verbose=False)
    assert result["status"] == "completed"
    assert (tmp_path / "out.fa").read_text().startswith(">A\nACGTACGT")


def test_recode_uncertain_replaces_ambiguity_codes(tmp_path: Path) -> None:
    source = _write(tmp_path / "amb.fa", ">A\nACGTN\n>B\nACRTA\n")
    wrappers.convert_format(source, str(tmp_path / "gap.fa"), input_format="fasta", data_type="dna",
                            output_format="fasta", recode_uncertain="gap", verbose=False)
    wrappers.convert_format(source, str(tmp_path / "miss.fa"), input_format="fasta", data_type="dna",
                            output_format="fasta", recode_uncertain="missing", verbose=False)
    gap_lines = (tmp_path / "gap.fa").read_text().split()
    miss_lines = (tmp_path / "miss.fa").read_text().split()
    assert "ACGT-" in gap_lines and "AC-TA" in gap_lines
    assert "ACGT?" in miss_lines and "AC?TA" in miss_lines


def test_fasta_and_phylip_input_need_a_data_type(tmp_path: Path, fasta: Path) -> None:
    with pytest.raises(ConfigurationError):
        wrappers.convert_format(str(fasta), str(tmp_path / "x.nex"), input_format="fasta", verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.run_dendropy_format(str(fasta), str(tmp_path / "y.nex"), input_format="fasta", verbose=False)


def test_convert_does_not_close_stdout(tmp_path: Path, fasta: Path) -> None:
    """The console convert() closes sys.stdout; the DNAKit function must leave the caller's stream open."""
    wrappers.convert_format(str(fasta), str(tmp_path / "keep.nex"), input_format="fasta", data_type="dna",
                            output_format="nexus", verbose=False)
    print("stdout still usable", file=sys.stdout)
    assert not sys.stdout.closed


def test_converted_output_refuses_existing_path(tmp_path: Path, fasta: Path) -> None:
    existing = _write(tmp_path / "exists.nex", "x")
    with pytest.raises(ConfigurationError):
        wrappers.convert_format(str(fasta), existing, input_format="fasta", data_type="dna", verbose=False)


def test_run_dendropy_format_matches_in_process_convert(tmp_path: Path) -> None:
    if not _HAS_DENDROPY_FORMAT:
        pytest.skip("dendropy-format console script is not installed in the dendropy environment")
    if not MSA.is_file():
        pytest.skip("shared alignment msa.fa is missing")
    in_process = tmp_path / "in_process.nex"
    program = tmp_path / "program.nex"
    wrappers.convert_format(str(MSA), str(in_process), input_format="fasta", data_type="dna",
                            output_format="nexus", verbose=False)
    result = wrappers.run_dendropy_format(str(MSA), str(program), input_format="fasta", data_type="dna",
                                          output_format="nexus", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert result["return_code"] == 0
    assert in_process.read_text() == program.read_text()


def test_run_dendropy_format_recode_and_phylip_strict(tmp_path: Path) -> None:
    if not _HAS_DENDROPY_FORMAT:
        pytest.skip("dendropy-format console script is not installed in the dendropy environment")
    source = _write(tmp_path / "amb.fa", ">A\nACGTN\n>B\nACRTA\n")
    result = wrappers.run_dendropy_format(source, str(tmp_path / "rec.phy"), input_format="fasta", data_type="dna",
                                          output_format="phylip-strict", recode_uncertain="gap", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert (tmp_path / "rec.phy").read_text().split()[:2] == ["2", "5"]
    with pytest.raises(ConfigurationError):
        wrappers.run_dendropy_format(source, str(tmp_path / "rec2.fa"), input_format="fasta", data_type="dna",
                                     recode_uncertain="missing", verbose=False)


def test_run_dendropy_format_missing_program(tmp_path: Path, fasta: Path) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        wrappers.run_dendropy_format(str(fasta), str(tmp_path / "o.nex"), input_format="fasta", data_type="dna",
                                     dendropy_format_path=str(tmp_path / "missing-dendropy-format"), verbose=False)
    assert "dendropy" in info.value.hint


# --- sumtrees and sumlabels ------------------------------------------------------------------


def test_summarize_tree_samples_consensus(tmp_path: Path) -> None:
    if not _HAS_DENDROPY_FORMAT:
        pytest.skip("sumtrees console script is not installed in the dendropy environment")
    trees = _write(tmp_path / "samples.tre", _TREES)
    result = wrappers.summarize_tree_samples([trees], str(tmp_path / "cons.tre"), summary_target="consensus",
                                             min_clade_freq=0.5, edges="mean-length", burnin=0, verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    summary = (tmp_path / "cons.tre").read_text()
    assert "Total of 3 trees analyzed" in summary
    assert "consensus tree with minimum clade frequency threshold of 0.5" in summary


def test_summarize_tree_samples_target_support_and_burnin(tmp_path: Path) -> None:
    if not _HAS_DENDROPY_FORMAT:
        pytest.skip("sumtrees console script is not installed in the dendropy environment")
    trees = _write(tmp_path / "samples.tre", _TREES)
    target = _write(tmp_path / "target.tre", _TARGET)
    result = wrappers.summarize_tree_samples([trees], str(tmp_path / "sup.tre"), target_tree=target, burnin=1,
                                             percentages=True, decimals=0, output_format="newick", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "sup.tre").read_text()
    assert text.startswith("[&U] ((A[&support=100.0]")  # burnin=1 leaves t2 and t3, both containing clade (A,B)


def test_summarize_tree_samples_refuses_bad_arguments(tmp_path: Path) -> None:
    trees = _write(tmp_path / "samples.tre", _TREES)
    target = _write(tmp_path / "target.tre", _TARGET)
    with pytest.raises(ConfigurationError):
        wrappers.summarize_tree_samples([trees], str(tmp_path / "a.tre"), summary_target="consensus",
                                        target_tree=target, verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.summarize_tree_samples([trees], str(tmp_path / "b.tre"), min_clade_freq=1.5, verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.summarize_tree_samples([], str(tmp_path / "c.tre"), verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.summarize_tree_samples([trees], str(tmp_path / "d.tre"), edges="sideways", verbose=False)


def test_summarize_tree_samples_missing_program(tmp_path: Path) -> None:
    trees = _write(tmp_path / "samples.tre", _TREES)
    with pytest.raises(BackendUnavailableError) as info:
        wrappers.summarize_tree_samples([trees], str(tmp_path / "e.tre"),
                                        sumtrees_path=str(tmp_path / "no-sumtrees"), verbose=False)
    assert "dendropy" in info.value.hint


def test_transfer_tree_labels(tmp_path: Path) -> None:
    if not _HAS_DENDROPY_FORMAT:
        pytest.skip("sumlabels console script is not installed in the dendropy environment")
    trees = _write(tmp_path / "samples.tre", _TREES)
    target = _write(tmp_path / "target.tre", _TARGET)
    result = wrappers.transfer_tree_labels([trees], target, str(tmp_path / "labelled.tre"), to_newick=True,
                                           verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "labelled.tre").read_text()
    assert "0.9/0.95" in text


def test_transfer_tree_labels_missing_program(tmp_path: Path) -> None:
    trees = _write(tmp_path / "samples.tre", _TREES)
    target = _write(tmp_path / "target.tre", _TARGET)
    with pytest.raises(BackendUnavailableError):
        wrappers.transfer_tree_labels([trees], target, str(tmp_path / "x.tre"),
                                      sumlabels_path=str(tmp_path / "no-sumlabels"), verbose=False)


# --- NCBI Entrez efetch (network replaced by a fake urlopen) ---------------------------------


class _FakeResponse(io.BytesIO):
    def __init__(self, body: bytes, url: str) -> None:
        super().__init__(body)
        self._url = url

    def geturl(self) -> str:
        return self._url


def test_efetch_builds_the_request_and_returns_records(monkeypatch: pytest.MonkeyPatch) -> None:
    requested: list[str] = []

    def fake_urlopen(url: str, *args: object, **kwargs: object) -> _FakeResponse:
        requested.append(url)
        return _FakeResponse(b">NM_1 demo\nACGT\n>NM_2 demo\nTTGG\n", url)

    monkeypatch.setattr(entrez, "urlopen", fake_urlopen)
    result = wrappers.fetch_entrez_records("nuccore", [123, "456"], verbose=False)
    assert requested and "db=nuccore" in requested[0] and "id=123%2C456" in requested[0]
    assert "rettype=fasta" in requested[0] and "retmode=text" in requested[0]
    assert result["status"] == "completed"
    assert result["service"] == "efetch"
    assert result["record_count"] == 2
    assert result["identifier_count"] == 2
    assert result["records"].startswith(">NM_1")
    assert len(result["access_date"]) == 10 and result["access_date"][4] == "-"


def test_efetch_validates_arguments(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: object, **kwargs: object) -> None:  # pragma: no cover - must not be reached
        raise AssertionError("the network must not be used for invalid arguments")

    monkeypatch.setattr(entrez, "urlopen", refuse)
    with pytest.raises(ConfigurationError):
        wrappers.fetch_entrez_records("nuccore", [], verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.fetch_entrez_records("nuccore", ["NM_000546"], verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.fetch_entrez_records("nuccore", [True], verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.fetch_entrez_records("nuccore", ["1"] * 201, verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.fetch_entrez_records("nuccore", ["1"], email="not-an-address", verbose=False)


def test_efetch_xml_has_no_record_count(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(entrez, "urlopen", lambda url, *a, **k: _FakeResponse(b"<xml/>", url))
    result = wrappers.fetch_entrez_records("nuccore", [7], rettype="gb", retmode="xml", verbose=False)
    assert result["record_count"] is None and result["records"] == "<xml/>"


# --- MUSCLE -------------------------------------------------------------------------------


def test_muscle_missing_program_gives_install_hint(tmp_path: Path, fasta: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    empty = tmp_path / "empty-bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.setattr(sys, "executable", str(empty / "python"))
    with pytest.raises(BackendUnavailableError) as info:
        wrappers.align_with_muscle(str(fasta), str(tmp_path / "aligned.fa"), verbose=False)
    assert "drive5.org" in info.value.hint or "bioconda" in info.value.hint
    with pytest.raises(BackendUnavailableError):
        wrappers.align_with_muscle(str(fasta), str(tmp_path / "aligned2.fa"),
                                   muscle_path=str(tmp_path / "no-muscle"), verbose=False)


def test_muscle_plumbing_with_a_stand_in_program(tmp_path: Path, fasta: Path) -> None:
    """A stand-in 'muscle' that copies standard input to standard output checks the file handling only.

    It is not a MUSCLE result: the sequences come back unchanged, which is what dendropy's stdin/stdout
    protocol needs from any aligner.
    """
    stand_in = tmp_path / "muscle"
    stand_in.write_text("#!/bin/sh\ncat\n")
    stand_in.chmod(0o755)
    result = wrappers.align_with_muscle(str(fasta), str(tmp_path / "copy.fa"), muscle_path=str(stand_in), verbose=False)
    assert result["status"] == "completed"
    assert result["sequence_count"] == 3 and result["alignment_length"] == 8
    assert (tmp_path / "copy.fa").read_text().startswith(">A\nACGTACGT")


@pytest.mark.skipif(not _HAS_MUSCLE, reason="MUSCLE is not installed (drive5.org/muscle5 or conda install -c bioconda muscle)")
def test_muscle_real_alignment(tmp_path: Path) -> None:
    result = wrappers.align_with_muscle(str(MSA), str(tmp_path / "msa_muscle.fa"), muscle_args=None, verbose=False)
    assert result["status"] == "completed", result
    assert result["sequence_count"] == 12


# --- PAUP* --------------------------------------------------------------------------------


def test_paup_functions_need_the_program_first(tmp_path: Path, fasta: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    empty = tmp_path / "empty-bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.setattr(sys, "executable", str(empty / "python"))
    trees = _write(tmp_path / "samples.tre", _TREES)
    _write(tmp_path / "target.tre", _TARGET)
    alignment = str(fasta)
    calls = [
        lambda: wrappers.paup_estimate_tree(alignment, str(tmp_path / "t1.nex"), verbose=False),
        lambda: wrappers.paup_estimate_ultrametric_tree(alignment, str(tmp_path / "t2.nex"), verbose=False),
        lambda: wrappers.paup_estimate_model(alignment, str(tmp_path / "t3.nex"), verbose=False),
        lambda: wrappers.paup_pscore_trees(trees, alignment, verbose=False),
        lambda: wrappers.paup_prune_taxa(trees, ["A"], str(tmp_path / "t4.nex"), verbose=False),
    ]
    for call in calls:
        with pytest.raises(BackendUnavailableError) as info:
            call()
        assert "paup.sc.fsu.edu" in info.value.hint


def test_paup_explicit_missing_path_is_refused(tmp_path: Path, fasta: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "t.nex"), paup_path=str(tmp_path / "no-paup"),
                                    verbose=False)


def test_paup_path_with_spaces_is_refused(tmp_path: Path, fasta: Path) -> None:
    spaced = tmp_path / "paup dir" / "paup"
    spaced.parent.mkdir()
    spaced.write_text("#!/bin/sh\nexit 0\n")
    spaced.chmod(0o755)
    with pytest.raises(ConfigurationError):
        wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "t.nex"), paup_path=str(spaced), verbose=False)


def test_paup_arguments_are_checked_before_the_program(tmp_path: Path, fasta: Path) -> None:
    trees = _write(tmp_path / "samples.tre", _TREES)
    with pytest.raises(ConfigurationError):
        wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "t.nex"), criterion="bayes", verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "t.nex"), substitution_types=3, verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.paup_pscore_trees(trees, str(fasta), pset_options=["gapmode=missing; quit"], verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.paup_prune_taxa(trees, ["NotATaxon"], str(tmp_path / "p.nex"), verbose=False)


@pytest.mark.skipif(not _HAS_PAUP, reason="PAUP* is not installed (paup.sc.fsu.edu or conda install -c bioconda paup)")
def test_paup_estimate_tree_real(tmp_path: Path) -> None:
    result = wrappers.paup_estimate_tree(str(MSA), str(tmp_path / "paup_tree.nex"), criterion="nj", verbose=False)
    assert result["status"] == "completed", result
    assert (tmp_path / "paup_tree.nex").stat().st_size > 0


def _stub_paup(tmp_path: Path) -> str:
    """A PAUP* path that exists but is never started: the upstream call is replaced by a fake in these tests."""
    stub = tmp_path / "paup"
    stub.write_text("#!/bin/sh\nexit 1\n")
    stub.chmod(0o755)
    return str(stub)


def test_paup_wrappers_pass_arguments_through(tmp_path: Path, fasta: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """dendropy's PAUP* functions are replaced by fakes: this checks the arguments and outputs, not PAUP* itself."""
    from dendropy.interop import paup as upstream

    dendropy = pytest.importorskip("dendropy")
    stub = _stub_paup(tmp_path)
    trees = _write(tmp_path / "samples.tre", _TREES)
    seen: dict[str, object] = {}

    def fake_estimate_tree(matrix, **kwargs):  # type: ignore[no-untyped-def]
        seen["estimate_tree"] = kwargs
        return dendropy.Tree.get(data="((A,B),(C,D));", schema="newick", taxon_namespace=matrix.taxon_namespace)

    def fake_estimate_model(matrix, **kwargs):  # type: ignore[no-untyped-def]
        seen["estimate_model"] = kwargs
        tree = dendropy.Tree.get(data="((A,B),(C,D));", schema="newick", taxon_namespace=matrix.taxon_namespace)
        return tree, {"log_likelihood": -12.5, "kappa": 2.0, "pA": 0.25}

    def fake_pscore(tree_list, matrix, **kwargs):  # type: ignore[no-untyped-def]
        seen["pscore"] = kwargs
        return [5.0 for _ in tree_list]

    def fake_prune(tree_list, taxa, **kwargs):  # type: ignore[no-untyped-def]
        seen["prune"] = [t.label for t in taxa]
        return tree_list

    monkeypatch.setattr(upstream, "estimate_tree", fake_estimate_tree)
    monkeypatch.setattr(upstream, "estimate_model", fake_estimate_model)
    monkeypatch.setattr(upstream, "pscore_trees", fake_pscore)
    monkeypatch.setattr(upstream, "prune_taxa_from_trees", fake_prune)

    tree = wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "est.nwk"), criterion="parsimony",
                                       substitution_types=2, paup_path=stub, verbose=False)
    assert tree["status"] == "completed" and tree["taxon_count"] == 3
    assert seen["estimate_tree"]["tree_est_criterion"] == "parsimony"
    assert seen["estimate_tree"]["num_subst"] == 2
    assert seen["estimate_tree"]["paup_path"] == str(Path(stub).resolve())
    assert (tmp_path / "est.nwk").read_text().strip().endswith(";")

    model = wrappers.paup_estimate_model(str(fasta), str(tmp_path / "model.nwk"), paup_path=stub, verbose=False)
    assert model["estimates"] == {"log_likelihood": -12.5, "kappa": 2.0, "pA": 0.25}
    assert seen["estimate_model"]["tree_model"] is None

    scores = wrappers.paup_pscore_trees(trees, str(fasta), pset_options=["gapmode=missing"], paup_path=stub,
                                        verbose=False)
    assert scores["scores"] == [5, 5, 5] and scores["tree_count"] == 3
    assert seen["pscore"]["pset_option_list"] == ["gapmode=missing"]

    pruned = wrappers.paup_prune_taxa(trees, ["A"], str(tmp_path / "pruned.nex"), paup_path=stub, verbose=False)
    assert pruned["removed"] == ["A"] and pruned["tree_count"] == 3
    assert seen["prune"] == ["A"]
    assert (tmp_path / "pruned.nex").read_text().startswith("#NEXUS")


def test_paup_upstream_failure_is_reported(tmp_path: Path, fasta: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from dendropy.interop import paup as upstream

    def failing(*args: object, **kwargs: object) -> None:
        raise OSError("PAUP* exited with status 1")

    monkeypatch.setattr(upstream, "estimate_tree", failing)
    with pytest.raises(BackendExecutionError, match="exited with status 1"):
        wrappers.paup_estimate_tree(str(fasta), str(tmp_path / "bad.nwk"), paup_path=_stub_paup(tmp_path),
                                    verbose=False)


def test_paup_upstream_error_is_wrapped() -> None:
    """A failure inside dendropy's PAUP* code is reported as BackendExecutionError, not as a bare Exception."""

    def broken() -> None:
        raise RuntimeError("boom")

    with pytest.raises(BackendExecutionError):
        wrappers._paup_call("probe", broken)


# --- shared ---------------------------------------------------------------------------------


def test_public_names_are_exported() -> None:
    assert set(wrappers.__all__) >= {"convert_format", "fetch_entrez_records", "align_with_muscle",
                                     "paup_estimate_tree", "summarize_tree_samples", "transfer_tree_labels"}
    assert len(wrappers.__all__) == 16
