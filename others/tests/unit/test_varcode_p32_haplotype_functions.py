"""Compare P32 haplotype/model wrappers with direct varcode on a local reference."""

import pytest

pytest.importorskip("varcode")
pytest.importorskip("pyensembl")
pytest.importorskip("pysam")

import test_varcode_p32_functions as _shared
import varcode
from test_varcode_p32_functions import _variant as _record
from varcode.annotators.fast import FastEffectAnnotator

from dnakit.molecular_lab import varcode_p32_haplotype_functions as m
from dnakit.molecular_lab.varcode_gap_functions import (
    _effect_record,
    _genome,
    _mutant_record,
    _variant,
)

# Re-use the shared fixture; pytest collects it under the name `world`.
world = _shared.world


@pytest.mark.parametrize(
    "labels",
    [
        ("missense_a", "missense_b"),
        ("silent_c", "missense_a"),
        ("missense_a", "missense_b", "silent_c"),
    ],
)
def test_known_cis_haplotype_matches_upstream(world, labels):
    items = [_record(world, label) for label in labels]
    genome, _ = _genome(world["gtf"], world["fasta"], world["cache"])
    transcript = genome.transcript_by_id("T1")
    direct = FastEffectAnnotator().annotate_haplotype(
        [_variant(genome, item) for item in items], transcript
    )
    got = m.varcode_annotate_cis_haplotype(
        world["gtf"], world["fasta"], items, "T1", cache_dir=world["cache"], verbose=False
    )
    assert got["effect"] == _effect_record(direct)
    assert got["mutant_transcript"] == _mutant_record(direct.mutant_transcript)


@pytest.mark.parametrize("labels", [("missense_a",), ("missense_a", "missense_b"), ("silent_c",)])
def test_transcript_model_matches_upstream(world, labels):
    items = [_record(world, label) for label in labels]
    genome, _ = _genome(world["gtf"], world["fasta"], world["cache"])
    transcript = genome.transcript_by_id("T1")

    def provider(contig, start, end):
        return varcode.reference_range(genome, contig, start, end)

    direct = varcode.predict_transcript_model_effect(
        [_variant(genome, item) for item in items], transcript, sequence_provider=provider
    )
    got = m.varcode_predict_transcript_model(
        world["gtf"], world["fasta"], items, "T1", cache_dir=world["cache"], verbose=False
    )
    assert got["effect"] == _effect_record(direct)
    assert got["candidate_count"] == len(getattr(direct, "candidates", ()) or ())
