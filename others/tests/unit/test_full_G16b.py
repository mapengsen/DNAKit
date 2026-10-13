"""Upstream-checked wrappers for sbol3, snapgene_reader, DNA Features Viewer,
pyGenomeViz and primer3-py.

primer3-py is GPL-2.0 and optional; these tests import it only through pytest.importorskip.
"""

from __future__ import annotations

import ast
import inspect
import io
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.molecular_lab import dna_features_viewer_functions as dfv_api
from dnakit.molecular_lab import primer3_functions as p3_api
from dnakit.molecular_lab import pygenomeviz_functions as pgv_api
from dnakit.molecular_lab import sbol3_functions as sbol_api
from dnakit.molecular_lab import snapgene_reader_functions as snap_api

GENBANK_SEQUENCE = "ATGCGTACGATCGTTTAAGGCCATGCAAGTCGATCGATCGGGATCCTAG"
TEMPLATE = (
    "ATGCGTACGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTACGATCGTAGCTAGCTAGCTAGCTGATCGATCGATCG"
    "TAGCTAGCTAGCATGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATGCATGCATGCTAGCTAGCTAGCATCG"
    "ATCGATGCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTACGATCGATCGATCGTAGCTACGATCGA"
)
ALL_WRAPPERS = [
    sbol_api.copy_sbol3_design,
    sbol_api.to_sbol3_display_id,
    sbol_api.check_sbol3_display_id,
    sbol_api.extract_sbol3_display_id,
    snap_api.snapgene_seqrecord_summary,
    dfv_api.annotate_sequence_feature,
    dfv_api.graphical_translation,
    pgv_api.summarize_genbank_records,
    pgv_api.fetch_genbank_by_accession,
    p3_api.calculate_primer3_end_stability,
    p3_api.design_primer3_pairs,
    p3_api.sanitize_primer3_sequence,
]


def _reverse_complement(text: str) -> str:
    return text.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _sbol_design(sbol3, directory: Path) -> Path:
    sequence = sbol3.Sequence(
        "https://example.org/a/seq1",
        elements="ATGCGTACGATCG",
        encoding=sbol3.IUPAC_DNA_ENCODING,
    )
    component = sbol3.Component("https://example.org/a/part1", sbol3.SBO_DNA)
    component.sequences = [sequence]
    document = sbol3.Document()
    document.add([sequence, component])
    path = directory / "design.ttl"
    document.write(str(path), sbol3.TURTLE)
    return path


def _genbank_text(sequence: str, accession: str) -> str:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq(sequence),
        id=accession,
        name="TEST1",
        description="fixture",
        annotations={"molecule_type": "DNA"},
    )
    record.features.append(SeqFeature(FeatureLocation(0, 12, 1), type="CDS"))
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


# ---------------------------------------------------------------- signatures


@pytest.mark.parametrize("function", ALL_WRAPPERS, ids=lambda item: item.__name__)
def test_wrappers_expose_keyword_only_verbose_default_true(function) -> None:
    parameter = inspect.signature(function).parameters["verbose"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is True


def test_module_exports_match_the_wrapped_functions() -> None:
    assert sbol_api.__all__ == [
        "check_sbol3_display_id",
        "copy_sbol3_design",
        "extract_sbol3_display_id",
        "to_sbol3_display_id",
    ]
    assert snap_api.__all__ == ["snapgene_seqrecord_summary"]
    assert dfv_api.__all__ == ["annotate_sequence_feature", "graphical_translation"]
    assert pgv_api.__all__ == ["fetch_genbank_by_accession", "summarize_genbank_records"]
    assert p3_api.__all__ == [
        "calculate_primer3_end_stability",
        "design_primer3_pairs",
        "sanitize_primer3_sequence",
    ]


def test_primer3_module_imports_primer3_only_inside_functions() -> None:
    tree = ast.parse(Path(p3_api.__file__).read_text(encoding="utf-8"))
    top_level = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom))]
    imported = set()
    for node in top_level:
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif node.module:
            imported.add(node.module)
    assert not any(name.split(".")[0] == "primer3" for name in imported)


# -------------------------------------------------------------------- sbol3


def test_copy_sbol3_design_repairs_references_and_matches_upstream_order(tmp_path) -> None:
    sbol3 = pytest.importorskip("sbol3")
    source = _sbol_design(sbol3, tmp_path)
    output = tmp_path / "copied.ttl"

    result = sbol_api.copy_sbol3_design(
        source, into_namespace="https://example.org/b", output=output, verbose=False
    )

    assert result["copied_count"] == 2
    assert [row["identity"] for row in result["copied"]] == [
        "https://example.org/b/part1",
        "https://example.org/b/seq1",
    ]
    assert result["output_format"] == "TURTLE"
    copied = sbol3.Document()
    copied.read(str(output))
    component = next(item for item in copied.objects if isinstance(item, sbol3.Component))
    assert [str(ref) for ref in component.sequences] == ["https://example.org/b/seq1"]

    # Upstream sbol3.copy gives the same result when the Sequence comes first.
    sequence = sbol3.Sequence(
        "https://example.org/a/seq1", elements="ATGCGTACGATCG", encoding=sbol3.IUPAC_DNA_ENCODING
    )
    component = sbol3.Component("https://example.org/a/part1", sbol3.SBO_DNA)
    component.sequences = [sequence]
    upstream = sbol3.Document()
    sbol3.copy(
        [sequence, component], into_namespace="https://example.org/b", into_document=upstream
    )
    upstream_component = next(
        item for item in upstream.objects if isinstance(item, sbol3.Component)
    )
    assert [str(ref) for ref in upstream_component.sequences] == ["https://example.org/b/seq1"]


def test_copy_sbol3_design_rejects_bad_namespace_output_and_missing_file(tmp_path) -> None:
    sbol3 = pytest.importorskip("sbol3")
    source = _sbol_design(sbol3, tmp_path)
    with pytest.raises(ConfigurationError, match="absolute URI"):
        sbol_api.copy_sbol3_design(source, into_namespace="not a uri", verbose=False)
    with pytest.raises(ConfigurationError, match="must end with"):
        sbol_api.copy_sbol3_design(
            source, into_namespace="https://example.org/b", output=tmp_path / "x.txt", verbose=False
        )
    with pytest.raises(ConfigurationError, match="does not exist"):
        sbol_api.copy_sbol3_design(
            tmp_path / "absent.ttl", into_namespace="https://example.org/b", verbose=False
        )


def test_to_sbol3_display_id_matches_hand_rule_and_upstream() -> None:
    sbol3 = pytest.importorskip("sbol3")
    result = sbol_api.to_sbol3_display_id("pTac promoter 1", verbose=False)
    assert result["display_id"] == "pTac_promoter_1"
    assert result["valid"] is True
    assert sbol_api.to_sbol3_display_id("1abc", verbose=False)["display_id"] == "_1abc"
    assert result["display_id"] == sbol3.string_to_display_id("pTac promoter 1")
    with pytest.raises(ConfigurationError):
        sbol_api.to_sbol3_display_id("   ", verbose=False)


def test_check_sbol3_display_id_matches_upstream_validity() -> None:
    sbol3 = pytest.importorskip("sbol3")
    assert sbol_api.check_sbol3_display_id("pTac_1", verbose=False)["valid"] is True
    assert sbol_api.check_sbol3_display_id("1bad", verbose=False)["valid"] is False
    assert sbol_api.check_sbol3_display_id("a-b", verbose=False)["valid"] is False
    for text in ("pTac_1", "1bad", "a-b", "_x"):
        assert sbol_api.check_sbol3_display_id(text, verbose=False)["valid"] == bool(
            sbol3.is_valid_display_id(text)
        )


def test_extract_sbol3_display_id_reads_url_path_and_rejects_invalid_tail() -> None:
    sbol3 = pytest.importorskip("sbol3")
    url = "https://example.org/lab/pTac"
    result = sbol_api.extract_sbol3_display_id(url, verbose=False)
    assert result["display_id"] == "pTac"
    assert result["found"] is True
    assert result["display_id"] == sbol3.extract_display_id(url)
    plain = sbol_api.extract_sbol3_display_id("not a url", verbose=False)
    assert plain["display_id"] is None
    assert plain["found"] is False
    with pytest.raises(ConfigurationError, match="not a valid displayId"):
        sbol_api.extract_sbol3_display_id("https://example.org/lab/pTac/1.0", verbose=False)


# ----------------------------------------------------------- snapgene_reader


def test_snapgene_missing_file_raises_configuration_error(tmp_path) -> None:
    pytest.importorskip("snapgene_reader")
    with pytest.raises(ConfigurationError, match="does not exist"):
        snap_api.snapgene_seqrecord_summary(tmp_path / "absent.dna", verbose=False)


def test_snapgene_malformed_file_raises_configuration_error(tmp_path) -> None:
    snapgene = pytest.importorskip("snapgene_reader")
    bad = tmp_path / "bad.dna"
    bad.write_bytes(b"not a SnapGene file")
    with pytest.raises(ValueError, match="Wrong format"):
        snapgene.snapgene_file_to_seqrecord(filepath=str(bad))
    with pytest.raises(ConfigurationError, match=r"could not parse bad\.dna"):
        snap_api.snapgene_seqrecord_summary(bad, verbose=False)


@pytest.mark.skip(
    reason="no SnapGene .dna sample exists in the local environment; "
    "the positive path cannot run offline"
)
def test_snapgene_seqrecord_summary_positive_path(tmp_path) -> None:
    raise AssertionError("requires a SnapGene .dna sample file")


def test_snapgene_record_summary_helper_sorts_features_without_a_sample() -> None:
    pytest.importorskip("Bio.SeqRecord")
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ATGCATGCAT"), id="synthetic", name="synthetic")
    record.features.append(SeqFeature(FeatureLocation(6, 9, -1), type="terminator"))
    record.features.append(
        SeqFeature(FeatureLocation(0, 4, 1), type="promoter", qualifiers={"label": ["pTac"]})
    )
    summary = snap_api._record_summary(
        record, Path("synthetic.dna"), max_bases=100, max_features=10
    )
    assert summary["length"] == 10
    assert [row["type"] for row in summary["features"]] == ["promoter", "terminator"]
    assert summary["features"][0]["label"] == "pTac"
    assert summary["features"][1]["strand"] == -1
    with pytest.raises(ConfigurationError, match="max_features"):
        snap_api._record_summary(record, Path("synthetic.dna"), max_bases=100, max_features=1)


# -------------------------------------------------------- dna_features_viewer


def test_annotate_sequence_feature_matches_direct_upstream_annotation() -> None:
    dfv = pytest.importorskip("dna_features_viewer")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ATGCCCGGGTAA"), id="sequence")
    dfv.annotate_biopython_record(
        record, location=(3, 9, 1), feature_type="CDS", margin=0, label="core"
    )
    direct = record.features[-1]
    result = dfv_api.annotate_sequence_feature(
        "ATGCCCGGGTAA",
        start=3,
        end=9,
        feature_type="CDS",
        qualifiers={"label": "core"},
        verbose=False,
    )
    assert result["feature"] == {
        "type": "CDS",
        "start": 3,
        "end": 9,
        "strand": 1,
        "qualifiers": {"label": "core"},
    }
    assert result["feature"]["start"] == int(direct.location.start)
    assert result["feature"]["end"] == int(direct.location.end)


def test_annotate_sequence_feature_defaults_and_rejects_bad_input() -> None:
    pytest.importorskip("dna_features_viewer")
    full = dfv_api.annotate_sequence_feature("ATGCCCGGGTAA", strand=-1, verbose=False)
    assert full["feature"] == {
        "type": "misc_feature",
        "start": 0,
        "end": 12,
        "strand": -1,
        "qualifiers": {},
    }
    with pytest.raises(ConfigurationError, match="reserved"):
        dfv_api.annotate_sequence_feature("ATGC", qualifiers={"location": "x"}, verbose=False)
    with pytest.raises(ConfigurationError, match="smaller than end"):
        dfv_api.annotate_sequence_feature("ATGCCC", start=4, end=3, verbose=False)


@pytest.mark.filterwarnings("ignore:Partial codon")
def test_graphical_translation_hand_calculation_and_upstream_agreement() -> None:
    biotools = pytest.importorskip("dna_features_viewer.biotools")
    forward = dfv_api.graphical_translation("ATGGCCTAA", verbose=False)
    assert forward["translation"] == "MA*"
    assert forward["residue_count"] == 3

    reverse = dfv_api.graphical_translation("TTAGCCAT", end=8, strand=-1, verbose=False)
    assert reverse["translation"] == "AM"
    assert reverse["translation"] == biotools.extract_graphical_translation("TTAGCCAT", (0, 8, -1))

    long_form = dfv_api.graphical_translation("TTAGCCATG", strand=-1, long_form=True, verbose=False)
    assert long_form["translation"] == ["*", "Gly", "His"]
    assert long_form["translation"] == biotools.extract_graphical_translation(
        "TTAGCCATG", (0, 9, -1), long_form=True
    )
    with pytest.raises(ConfigurationError, match="smaller than end"):
        dfv_api.graphical_translation("ATGGCC", start=5, end=2, verbose=False)


# ------------------------------------------------------------- pygenomeviz


def test_summarize_genbank_records_matches_upstream_parser(tmp_path) -> None:
    pytest.importorskip("pygenomeviz")
    from pygenomeviz.parser.genbank import Genbank

    path = tmp_path / "fixture.gb"
    path.write_text(_genbank_text(GENBANK_SEQUENCE, "TEST1.1"), encoding="utf-8")
    result = pgv_api.summarize_genbank_records([path], verbose=False)
    genome = result["genomes"][0]
    direct = Genbank(str(path))
    assert genome["name"] == "fixture"
    assert genome["full_length"] == direct.full_genome_length == len(GENBANK_SEQUENCE)
    assert genome["range"] == [0, len(GENBANK_SEQUENCE)]
    assert genome["records"] == [
        {"seqid": "TEST1.1", "length": len(GENBANK_SEQUENCE), "feature_count": 1}
    ]
    gc_count = sum(base in "GC" for base in GENBANK_SEQUENCE)
    assert genome["gc_percent"] == pytest.approx(gc_count / len(GENBANK_SEQUENCE) * 100)
    assert genome["gc_percent"] == pytest.approx(direct.calc_genome_gc_content())


def test_summarize_genbank_records_range_and_reverse_specs(tmp_path) -> None:
    pytest.importorskip("pygenomeviz")
    path = tmp_path / "fixture.gb"
    path.write_text(_genbank_text(GENBANK_SEQUENCE, "TEST1.1"), encoding="utf-8")
    ranged = pgv_api.summarize_genbank_records([f"{path}:10-30"], verbose=False)["genomes"][0]
    assert ranged["range"] == [10, 30]
    assert ranged["range_size"] == 20
    reverse = pgv_api.summarize_genbank_records([f"{path}::-1"], verbose=False)["genomes"][0]
    assert reverse["reverse"] is True
    assert reverse["full_length"] == len(GENBANK_SEQUENCE)
    with pytest.raises(ConfigurationError, match="does not exist"):
        pgv_api.summarize_genbank_records([tmp_path / "absent.gb"], verbose=False)


def test_fetch_genbank_by_accession_uses_a_mocked_entrez_client(tmp_path, monkeypatch) -> None:
    pytest.importorskip("pygenomeviz")
    from Bio import Entrez

    calls: list[dict[str, object]] = []
    text = _genbank_text(GENBANK_SEQUENCE, "JX128258.1")

    def fake_efetch(**kwargs):
        calls.append(kwargs)
        return io.StringIO(text)

    monkeypatch.setattr(Entrez, "efetch", fake_efetch)
    output = tmp_path / "fetched.gb"
    result = pgv_api.fetch_genbank_by_accession(
        "JX128258.1", output=output, email="lab@example.org", verbose=False
    )
    assert len(calls) == 1
    assert calls[0]["id"] == "JX128258.1"
    assert calls[0]["db"] == "nucleotide"
    assert result["record_id"] == "JX128258.1"
    assert result["length"] == len(GENBANK_SEQUENCE)
    assert result["output"] == str(output.resolve())
    assert output.read_text(encoding="utf-8") == text
    assert result["provenance"]["network"] is True


def test_fetch_genbank_by_accession_maps_network_failure_and_bad_input(monkeypatch) -> None:
    pytest.importorskip("pygenomeviz")
    from Bio import Entrez

    def offline(**kwargs):
        raise OSError("network unreachable (mocked)")

    monkeypatch.setattr(Entrez, "efetch", offline)
    with pytest.raises(BackendExecutionError, match="network unreachable"):
        pgv_api.fetch_genbank_by_accession("JX128258.1", verbose=False)

    def must_not_call(**kwargs):
        raise AssertionError("validation must run before any network call")

    monkeypatch.setattr(Entrez, "efetch", must_not_call)
    with pytest.raises(ConfigurationError, match="accession"):
        pgv_api.fetch_genbank_by_accession("bad acc", verbose=False)
    with pytest.raises(ConfigurationError, match="email"):
        pgv_api.fetch_genbank_by_accession("JX128258.1", email="not-an-email", verbose=False)


# ---------------------------------------------------------------- primer3


def test_end_stability_matches_direct_upstream_call() -> None:
    primer3 = pytest.importorskip("primer3")
    primer, template = "GCTAGCTCAGTCCTAGG", "CCTAGGACTGAGCTAGC"
    result = p3_api.calculate_primer3_end_stability(primer, template, verbose=False)
    direct = primer3.calc_end_stability(primer, template)
    assert result["primer_length"] == 17
    assert result["structure_found"] is True
    assert result["dg_cal_per_mol"] == pytest.approx(direct.dg)
    assert result["dh_cal_per_mol"] == pytest.approx(direct.dh)
    assert result["tm_c"] == pytest.approx(direct.tm)
    assert result["provenance"]["license"] == "GPL-2.0"


def test_end_stability_rejects_invalid_bases_and_conditions() -> None:
    pytest.importorskip("primer3")
    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        p3_api.calculate_primer3_end_stability("ACGTR", "ACGT", verbose=False)
    with pytest.raises(ConfigurationError, match="mv_conc"):
        p3_api.calculate_primer3_end_stability("ACGT", "ACGT", mv_conc=-1.0, verbose=False)


def test_design_primer3_pairs_hand_checks_and_matches_upstream() -> None:
    primer3 = pytest.importorskip("primer3")
    result = p3_api.design_primer3_pairs(
        TEMPLATE, min_product=80, max_product=200, num_return=2, verbose=False
    )
    assert result["pair_count"] == 2
    first = result["pairs"][0]
    left, right = first["left"], first["right"]
    assert first["product_size"] == right["start"] + right["length"] - left["start"]
    assert left["sequence"] == TEMPLATE[left["start"] : left["start"] + left["length"]]
    assert right["sequence"] == _reverse_complement(
        TEMPLATE[right["start"] : right["start"] + right["length"]]
    )

    seq_args = {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": TEMPLATE}
    global_args = {
        "PRIMER_PRODUCT_SIZE_RANGE": [[80, 200]],
        "PRIMER_NUM_RETURN": 2,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 40.0,
        "PRIMER_MAX_GC": 60.0,
    }
    raw = primer3.design_primers(seq_args, global_args)
    assert left["sequence"] == raw["PRIMER_LEFT_0_SEQUENCE"]
    assert right["sequence"] == raw["PRIMER_RIGHT_0_SEQUENCE"]
    assert first["product_size"] == raw["PRIMER_PAIR_0_PRODUCT_SIZE"]


def test_design_primer3_pairs_passes_target_and_rejects_bad_regions() -> None:
    pytest.importorskip("primer3")
    result = p3_api.design_primer3_pairs(
        TEMPLATE, target=(100, 50), min_product=80, max_product=200, verbose=False
    )
    assert result["target"] == [100, 50]
    with pytest.raises(ConfigurationError, match="extends past"):
        p3_api.design_primer3_pairs(TEMPLATE, target=(len(TEMPLATE) - 10, 20), verbose=False)
    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        p3_api.design_primer3_pairs("ACGT" * 15 + "R", verbose=False)
    with pytest.raises(ConfigurationError, match="max_product"):
        p3_api.design_primer3_pairs(TEMPLATE, min_product=200, max_product=100, verbose=False)


def test_sanitize_primer3_sequence_replaces_ambiguity_codes() -> None:
    primer3 = pytest.importorskip("primer3")
    result = p3_api.sanitize_primer3_sequence("ACGTRYN", verbose=False)
    assert result["sanitized"] == "ACGTNNN"
    assert result["changed_count"] == 2
    assert result["sanitized"] == primer3.p3helpers.sanitize_sequence("ACGTRYN")
    assert p3_api.sanitize_primer3_sequence("acgtn", verbose=False)["sanitized"] == "acgtn"
    with pytest.raises(ConfigurationError, match="rejected the sequence"):
        p3_api.sanitize_primer3_sequence("ACGTX", verbose=False)
