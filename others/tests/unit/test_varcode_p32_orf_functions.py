"""Reference-ORF annotations compared with varcode 10.16.3 and orfdata."""

import pytest

varcode = pytest.importorskip("varcode")
pytest.importorskip("pyensembl")
pytest.importorskip("pysam")
import test_varcode_p32_functions as _shared
from test_varcode_p32_functions import _variant as _record

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import varcode_p32_orf_functions as m
from dnakit.molecular_lab.varcode_gap_functions import _transcript, _variant

CONTEXT = {
    "taxon_id": 9606,
    "species": "homo_sapiens",
    "genome_build": "GRCh38",
    "annotation_name": "test",
    "annotation_version": "1",
}

# Re-use the shared fixture; pytest collects it under the name `world`.
world = _shared.world


@pytest.fixture()
def data(world):
    if not hasattr(varcode, "annotate_orf_variants"):
        pytest.skip("Installed varcode predates the reference-ORF API; not an annotation pass.")
    pytest.importorskip("orfdata")
    vc, od, genome, context = m._reference_genome(
        world["gtf"], world["fasta"], CONTEXT, world["cache"]
    )
    transcript = _transcript(genome, "T1")
    model = {
        "orf_id": "test-orf",
        "gene_id": "G1",
        "gene_name": "G1",
        "transcript_id": "T1",
        "contig": "chr1",
        "strand": "+",
        "coding_blocks": [[101, 280]],
        "protein_sequence": transcript.protein_sequence,
        "biotype": "CDS",
    }
    reference = m._model(od, genome, context, model).bind(transcript, context=context)
    return vc, od, genome, context, transcript, model, reference


def test_old_varcode_has_install_hint():
    if hasattr(varcode, "annotate_orf_variants"):
        pytest.skip("Old-version dependency path only.")
    with pytest.raises(BackendUnavailableError, match="reference-ORF"):
        m._api()


@pytest.mark.parametrize("label", ["missense_a", "silent_c", "missense_b"])
def test_explicit_orf_annotation_matches_upstream(world, data, label):
    vc, _od, genome, _context, _transcript, model, reference = data
    records = [_record(world, label)]
    expected = vc.annotate_orf_variants([_variant(genome, r) for r in records], reference).to_dict()
    got = m.varcode_annotate_reference_orf(
        world["gtf"],
        world["fasta"],
        records,
        model,
        reference_context=CONTEXT,
        cache_dir=world["cache"],
        verbose=False,
    )
    assert got["annotation"] == expected


@pytest.mark.parametrize("case", ["unchanged", "substitution", "start_missing"])
def test_observed_rna_orf_annotation_matches_upstream(world, data, case):
    vc, _od, _genome, _context, transcript, model, reference = data
    sequence = transcript.sequence
    origins = [["chr1", p, "+"] for p in range(101, 281)]
    if case == "substitution":
        sequence = sequence[:4] + ("C" if sequence[4] != "C" else "G") + sequence[5:]
    elif case == "start_missing":
        sequence, origins = sequence[3:], origins[3:]
    from varcode.effect_hypotheses import RealizedTranscriptProduct

    product = RealizedTranscriptProduct(
        transcript, sequence, None, evidence={"base_origins": [tuple(o) for o in origins]}
    )
    expected = vc.annotate_orf_rna_product(reference, product).to_dict()
    got = m.varcode_annotate_orf_rna_product(
        world["gtf"],
        world["fasta"],
        model,
        sequence,
        origins,
        reference_context=CONTEXT,
        cache_dir=world["cache"],
        verbose=False,
    )
    assert got["annotation"] == expected


@pytest.mark.parametrize("label", ["missense_a", "silent_c", "intergenic_e"])
def test_catalogue_prediction_matches_upstream(world, data, label):
    vc, od, genome, context, _transcript, model, reference = data
    records = [_record(world, label)]
    metadata = {"genome_build": "GRCh38"}
    catalog = od.ORFCatalog([reference.record], metadata)
    expected = vc.predict_orf_effects(
        [_variant(genome, r) for r in records], catalog=catalog, reference_context=context
    )
    got = m.varcode_predict_orf_catalogue(
        world["gtf"],
        world["fasta"],
        records,
        [model],
        reference_context=CONTEXT,
        cache_dir=world["cache"],
        verbose=False,
    )
    assert got["annotations"] == [a.to_dict() for a in expected]
