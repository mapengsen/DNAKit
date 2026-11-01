"""Tests for the Primer3 design, alignment and masking functions."""

from __future__ import annotations

import os
import random
import shutil

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import primer3_p14_functions as p3

_RNG = random.Random(2)
TEMPLATE = "".join(_RNG.choice("ACGT") for _ in range(600))
needs_primer3_py = pytest.mark.skipif(__import__("importlib").util.find_spec("primer3") is None, reason="primer3-py is not installed")


@needs_primer3_py
def test_pick_primer_list_returns_ranked_single_primers():
    result = p3.primer3_design_by_task(TEMPLATE, task="pick_primer_list", num_return=3, verbose=False)
    assert result["status"] == "completed" and len(result["left"]) == 3 and len(result["right"]) == 3 and result["pairs"] == []
    penalties = [row["penalty"] for row in result["left"]]
    assert penalties == sorted(penalties)
    first = result["left"][0]
    assert TEMPLATE[first["start"]:first["start"] + first["length"]] == first["sequence"]


@needs_primer3_py
def test_generic_and_check_primers():
    generic = p3.primer3_design_by_task(TEMPLATE, num_return=2, product_size_range=[[100, 300]], verbose=False)
    pair = generic["pairs"][0]
    assert 100 <= pair["product_size"] <= 300
    left = generic["left"][0]["sequence"]
    checked = p3.primer3_design_by_task(TEMPLATE, task="check_primers", left_primer=left, verbose=False)
    assert checked["left"][0]["sequence"] == left
    assert checked["left"][0]["start"] == generic["left"][0]["start"]
    assert abs(checked["left"][0]["tm"] - generic["left"][0]["tm"]) < 1e-6


@needs_primer3_py
def test_sequencing_cloning_and_discriminative_tasks_run():
    sequencing = p3.primer3_design_by_task(TEMPLATE, task="pick_sequencing_primers", target=[300, 20], verbose=False)
    assert sequencing["status"] == "completed" and sequencing["left"]
    discriminative = p3.primer3_design_by_task(TEMPLATE, task="pick_discriminative_primers", target=[300, 1], verbose=False)
    assert discriminative["status"] == "completed"
    cloning = p3.primer3_design_by_task(TEMPLATE, task="pick_cloning_primers", included_region=[50, 400], verbose=False)
    assert cloning["status"] == "completed" and "error" in cloning


@needs_primer3_py
def test_validation():
    with pytest.raises(ConfigurationError):
        p3.primer3_design_by_task(TEMPLATE, task="nonsense")
    with pytest.raises(ConfigurationError):
        p3.primer3_design_by_task(TEMPLATE, settings={"PRIMER_TASK": "generic"})
    with pytest.raises(ConfigurationError):
        p3.primer3_design_by_task(TEMPLATE, target=[590, 50])
    with pytest.raises(ConfigurationError):
        p3.primer3_design_by_task("ACGTX")


needs_ntdpal = pytest.mark.skipif(shutil.which("ntdpal") is None, reason="Primer3 programs are not installed")


@needs_ntdpal
@pytest.mark.parametrize("mode,score,length", [("local", 6.0, 6), ("global", 6.0, 6), ("global_end_anchored", 6.0, 6), ("local_end_anchored", 6.0, 6)])
def test_ntdpal_modes(mode, score, length):
    result = p3.ntdpal_align_oligos("GCTAGCTAGGATCC", "GGATCCAAAA", mode=mode, verbose=False)
    assert result["status"] == "completed" and (result["score"], result["length"]) == (score, length)
    assert result["pairs"][0] == [8, 0] and result["pairs"][-1] == [13, 5]


@needs_ntdpal
def test_ntdpal_local_reports_end_and_gc_scores():
    result = p3.ntdpal_align_oligos("AAAGGCCTTT", "GGCC", mode="local", gc_weighted_scores=True, verbose=False)
    assert result["score"] == 12.0 and result["end_1"] == 6 and result["end_2"] == 3


def test_ntdpal_validation_and_missing_program(tmp_path):
    with pytest.raises(ConfigurationError):
        p3.ntdpal_align_oligos("ACGT", "AC-T")
    with pytest.raises(ConfigurationError):
        p3.ntdpal_align_oligos("ACGT", "ACGT", mode="x")
    with pytest.raises(BackendUnavailableError):
        p3.ntdpal_align_oligos("ACGT", "ACGT", primer3_path=str(tmp_path))


KMER_LISTS = os.environ.get("PRIMER3_KMER_LISTS", "")


@pytest.mark.skipif(shutil.which("primer3_masker") is None or not os.path.isfile(os.path.join(KMER_LISTS, "homo_sapiens_16.list")), reason="primer3_masker or its k-mer lists (PRIMER3_KMER_LISTS) are not available")
def test_masker_probability_and_count_modes(tmp_path):
    fasta = tmp_path / "s.fa"
    fasta.write_text(">s\nACGTACGTAGCTAGCTAGCTAGCTAAAAAAAAAAAAAAAAAAAAAAAA\n")
    probability = p3.primer3_mask_repeats(str(fasta), str(tmp_path / "p.fa"), kmer_lists_path=KMER_LISTS, verbose=False)
    assert probability["status"] == "completed" and probability["mode"] == "probability"
    assert "N" in (tmp_path / "p.fa").read_text().split("\n")[1]
    counted = p3.primer3_mask_repeats(str(fasta), str(tmp_path / "c.fa"), kmer_lists_path=KMER_LISTS, absolute_cutoff=100, kmer_list=os.path.join(KMER_LISTS, "homo_sapiens_16.list"), soft_mask=True, verbose=False)
    assert counted["status"] == "completed" and counted["mode"] == "absolute_count"
    assert any(char.islower() for char in (tmp_path / "c.fa").read_text().split("\n")[1])


def test_masker_validation(tmp_path):
    fasta = tmp_path / "s.fa"
    fasta.write_text(">s\nACGT\n")
    with pytest.raises(ConfigurationError):
        p3.primer3_mask_repeats(str(fasta), str(tmp_path / "o"), kmer_lists_path=str(tmp_path), absolute_cutoff=5)
    with pytest.raises(ConfigurationError):
        p3.primer3_mask_repeats(str(fasta), str(tmp_path / "o"), kmer_lists_path=str(tmp_path), masking_direction="up")
