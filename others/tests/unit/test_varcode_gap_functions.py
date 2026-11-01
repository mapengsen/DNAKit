"""Tests for dnakit.molecular_lab.varcode_gap_functions (varcode on a local GTF and genome FASTA).

Every test runs offline. The reference is built from the shared genome (consistency/data/genome.fa) and a small
GTF written into tmp_path: one transcript on the forward strand of chr1 and one on the reverse strand of chr2,
each spanning an open reading frame found in the genome. Expected CDS, protein and effect classes are derived from
the genome sequence in these tests, not copied from varcode output. Network access is never used; the one live
check (varcode_infer_genome_for_reference_name with a release name) runs only when DNAKIT_LIVE_TESTS=1.
"""

from __future__ import annotations

import inspect
import itertools
import os
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("varcode")
pytest.importorskip("pyensembl")
pysam = pytest.importorskip("pysam")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.molecular_lab import varcode_gap_functions as vg  # noqa: E402

DATA = (
    Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENOME = DATA / "genome.fa"
# Standard genetic code, codons in TCAG order (independent of the module's own table).
CODON = {
    "".join(bases): amino
    for bases, amino in zip(
        itertools.product("TCAG", repeat=3),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
}
STOPS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _revcomp(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


def _translate(coding: str) -> str:
    return "".join(CODON[coding[i : i + 3]] for i in range(0, len(coding) - 2, 3)).rstrip("*")


def _orf(sequence: str) -> tuple[int, int]:
    """First ATG (0-based index) whose frame reaches a stop codon; returns [start, end) including the stop."""
    for index in range(1000, len(sequence) - 3):
        if sequence[index : index + 3] != "ATG":
            continue
        end = index
        while end + 3 <= len(sequence) and sequence[end : end + 3] not in STOPS:
            end += 3
        if end + 3 <= len(sequence) and 60 <= end - index <= 600:
            return index, end + 3
    raise AssertionError("no ORF found in the shared genome")


def _gtf_row(contig: str, strand: str, feature: str, start: int, end: int, attributes: str) -> str:
    score_or_frame = "0" if feature == "CDS" else "."
    return "\t".join(
        [contig, "test", feature, str(start), str(end), ".", strand, score_or_frame, attributes]
    )


def _gtf_text(genes: list[dict[str, Any]]) -> str:
    lines: list[str] = []
    for gene in genes:
        contig, strand, gid, tx = gene["contig"], gene["strand"], gene["gene_id"], gene["tx"]
        gene_attrs = f'gene_id "{gid}"; gene_name "{gid}"; gene_biotype "protein_coding";'
        tx_attrs = f'gene_id "{gid}"; transcript_id "{tx}"; transcript_biotype "protein_coding";'
        lines.append(_gtf_row(contig, strand, "gene", *gene["span"], gene_attrs))
        lines.append(_gtf_row(contig, strand, "transcript", *gene["span"], tx_attrs))
        for number, (start, end) in enumerate(gene["exons"], start=1):
            lines.append(
                _gtf_row(contig, strand, "exon", start, end, f'{tx_attrs} exon_number "{number}";')
            )
        for number, (start, end) in enumerate(gene["cds"], start=1):
            cds_attrs = f'{tx_attrs} exon_number "{number}"; protein_id "{gene["protein"]}";'
            lines.append(_gtf_row(contig, strand, "CDS", start, end, cds_attrs))
        lines.append(
            _gtf_row(
                contig, strand, "start_codon", *gene["start_codon"], f'{tx_attrs} exon_number "1";'
            )
        )
        lines.append(
            _gtf_row(
                contig, strand, "stop_codon", *gene["stop_codon"], f'{tx_attrs} exon_number "1";'
            )
        )
    return "\n".join(lines) + "\n"


@pytest.fixture(scope="module")
def genome_text() -> str:
    if not GENOME.exists():
        pytest.skip("shared consistency genome.fa not present")
    with pysam.FastaFile(str(GENOME)) as fasta:
        return fasta.fetch("chr1")


@pytest.fixture(scope="module")
def layout(genome_text: str, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """A forward-strand transcript T1 on chr1 and a reverse-strand transcript T2 on chr2."""
    with pysam.FastaFile(str(GENOME)) as fasta:
        chr1 = fasta.fetch("chr1")
        chr2 = fasta.fetch("chr2")
    a, b = _orf(chr1)
    # Reverse strand: find the ORF in the reverse complement of chr2 and map it back to forward coordinates.
    rc = _revcomp(chr2)
    c, d = _orf(rc)
    length = len(chr2)
    forward = (length - d + 1, length - c)  # 1-based inclusive transcript span on chr2
    directory = tmp_path_factory.mktemp("varcode_reference")
    genes = [
        {
            "gene_id": "G1",
            "tx": "T1",
            "protein": "P1",
            "contig": "chr1",
            "strand": "+",
            "span": (a + 1, b),
            "exons": [(a + 1, b)],
            "cds": [(a + 1, b)],
            "start_codon": (a + 1, a + 3),
            "stop_codon": (b - 2, b),
        },
        {
            "gene_id": "G2",
            "tx": "T2",
            "protein": "P2",
            "contig": "chr2",
            "strand": "-",
            "span": forward,
            "exons": [forward],
            "cds": [forward],
            "start_codon": (forward[1] - 2, forward[1]),
            "stop_codon": (forward[0], forward[0] + 2),
        },
    ]
    gtf = directory / "genes.gtf"
    gtf.write_text(_gtf_text(genes), encoding="utf-8")
    cds_t1 = chr1[a:b]
    cds_t2 = _revcomp(chr2[forward[0] - 1 : forward[1]])
    return {
        "gtf": str(gtf),
        "cache": str(directory / "cache"),
        "chr1": chr1,
        "chr2": chr2,
        "t1_span": (a + 1, b),
        "t1_cds": cds_t1,
        "t1_protein": _translate(cds_t1),
        "t2_span": forward,
        "t2_cds": cds_t2,
        "t2_protein": _translate(cds_t2),
    }


def _call(function: Callable[..., Any], layout: dict[str, Any], *args: Any, **kwargs: Any) -> Any:
    return function(
        layout["gtf"], str(GENOME), *args, cache_dir=layout["cache"], verbose=False, **kwargs
    )


def _snv_in_cds(layout: dict[str, Any], wanted: str) -> tuple[int, str, str]:
    """Find a single-base change in T1's CDS whose protein effect is ``wanted`` (silent, missense or nonsense)."""
    cds = layout["t1_cds"]
    protein = layout["t1_protein"]
    start = layout["t1_span"][0]
    for codon in range(3, len(cds) // 3 - 2):
        offset = codon * 3 + 2  # third codon position, 0-based within the CDS
        ref = cds[offset]
        for alt in "ACGT":
            if alt == ref:
                continue
            mutated = cds[:offset] + alt + cds[offset + 1 :]
            new_protein = _translate(mutated)
            if wanted == "silent" and new_protein == protein:
                return start + offset, ref, alt
            if wanted == "nonsense" and "*" in new_protein[: len(protein)]:
                return start + offset, ref, alt
            if (
                wanted == "missense"
                and len(new_protein) == len(protein)
                and new_protein != protein
                and "*" not in new_protein
            ):
                return start + offset, ref, alt
    raise AssertionError(f"no {wanted} change found in the CDS")


# --- module contract ----------------------------------------------------------------------------


def test_public_functions_have_keyword_only_verbose_and_docstrings() -> None:
    for name in vg.__all__:
        function = getattr(vg, name)
        signature = inspect.signature(function)
        verbose = signature.parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert function.__doc__, name
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, f"{name}.{parameter.name}"


def test_module_exports_are_callable() -> None:
    assert len(vg.__all__) == len(set(vg.__all__))
    assert all(callable(getattr(vg, name)) for name in vg.__all__)


# --- local reference and reference sequence ---------------------------------------------------


def test_build_local_reference_indexes_both_contigs(layout: dict[str, Any], tmp_path: Path) -> None:
    result = vg.build_local_reference(
        layout["gtf"], str(GENOME), cache_dir=str(tmp_path / "fresh"), verbose=False
    )
    assert result["status"] == "completed"
    assert result["transcripts"] == 2
    assert result["genes"] == 2
    assert result["contigs"] == 2
    assert result["reused"] is False
    again = vg.build_local_reference(
        layout["gtf"], str(GENOME), cache_dir=str(tmp_path / "fresh"), verbose=False
    )
    assert again["reused"] is True
    assert again["annotation_name"] == result["annotation_name"]


def test_reference_base_and_range_match_the_fasta(layout: dict[str, Any]) -> None:
    start = layout["t1_span"][0]
    base = _call(vg.varcode_reference_base, layout, "chr1", start)["base"]
    assert base == layout["chr1"][start - 1]
    window = _call(vg.varcode_reference_range, layout, "chr1", start, start + 9)["sequence"]
    assert window == layout["chr1"][start - 1 : start + 9]


def test_reference_range_rejects_reversed_interval(layout: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        _call(vg.varcode_reference_range, layout, "chr1", 500, 400)


def test_missing_gtf_or_fasta_is_a_configuration_error(
    layout: dict[str, Any], tmp_path: Path
) -> None:
    with pytest.raises(ConfigurationError):
        vg.build_local_reference(
            str(tmp_path / "absent.gtf"), str(GENOME), cache_dir=str(tmp_path), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vg.build_local_reference(
            layout["gtf"], str(tmp_path / "absent.fa"), cache_dir=str(tmp_path), verbose=False
        )


def test_gtf_with_unknown_contig_is_rejected(tmp_path: Path) -> None:
    gtf = tmp_path / "bad.gtf"
    gtf.write_text(
        'chrZ\ttest\texon\t1\t50\t.\t+\t.\ttranscript_id "TX"; gene_id "GX";\n', encoding="utf-8"
    )
    with pytest.raises(ConfigurationError):
        vg.build_local_reference(
            str(gtf), str(GENOME), cache_dir=str(tmp_path / "c"), verbose=False
        )


# --- effect prediction on the local reference -------------------------------------------------


def test_silent_change_is_reported_as_silent(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "silent")
    record = _call(
        vg.varcode_predict_variant_effects,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
    )
    classes = [effect["class"] for effect in record["effects"]]
    assert classes == ["Silent"]
    assert record["effects"][0]["transcript_id"] == "T1"


def test_missense_and_nonsense_changes_are_not_silent(layout: dict[str, Any]) -> None:
    for wanted in ("missense", "nonsense"):
        start, ref, alt = _snv_in_cds(layout, wanted)
        record = _call(
            vg.varcode_predict_variant_effects,
            layout,
            {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        )
        assert record["effects"], wanted
        assert all(effect["class"] != "Silent" for effect in record["effects"]), wanted


def test_start_codon_change_is_start_loss(layout: dict[str, Any]) -> None:
    start = layout["t1_span"][0]
    first_codon = layout["t1_cds"][:3]
    start_codons = {"ATG", "CTG", "TTG"}  # standard table start codons, as varcode reports them
    alt = next(
        base
        for base in "ACGT"
        if base != first_codon[0] and (base + first_codon[1:]) not in start_codons
    )
    record = _call(
        vg.varcode_predict_variant_effect_on_transcript,
        layout,
        {"contig": "chr1", "start": start, "ref": first_codon[0], "alt": alt},
        "T1",
    )
    assert record["effect"]["class"] == "StartLoss"


def test_variant_on_reverse_strand_transcript_is_annotated(layout: dict[str, Any]) -> None:
    start, end = layout["t2_span"]
    position = (start + end) // 2
    ref = layout["chr2"][position - 1]
    alt = "A" if ref != "A" else "C"
    record = _call(
        vg.varcode_predict_variant_effects,
        layout,
        {"contig": "chr2", "start": position, "ref": ref, "alt": alt},
    )
    assert record["effects"]
    assert record["effects"][0]["transcript_id"] == "T2"


def test_variant_outside_transcripts_is_intergenic(layout: dict[str, Any]) -> None:
    record = _call(
        vg.varcode_predict_variant_effects,
        layout,
        {"contig": "chr1", "start": 5000, "ref": layout["chr1"][4999], "alt": "A"},
    )
    assert record["effects"]
    assert all("transcript_id" not in effect for effect in record["effects"])
    assert record["effects"][0]["class"] == "Intergenic"


def test_reference_mismatch_is_reported_as_failure_when_requested(layout: dict[str, Any]) -> None:
    start = layout["t1_span"][0] + 10
    wrong = "A" if layout["chr1"][start - 1] != "A" else "C"
    record = _call(
        vg.varcode_predict_variant_effect_on_transcript_or_failure,
        layout,
        {"contig": "chr1", "start": start, "ref": wrong, "alt": "G"},
        "T1",
    )
    assert record["status"] == "failed"
    assert record["effect"] is None
    assert record["failure"]


def test_mutant_transcript_protein_matches_translated_mutant_cds(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    record = _call(
        vg.varcode_apply_variant_to_transcript,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        "T1",
    )
    offset = start - layout["t1_span"][0]
    mutated = layout["t1_cds"][:offset] + alt + layout["t1_cds"][offset + 1 :]
    expected = _translate(mutated)
    observed = record["mutant_transcript"]["mutant_protein_sequence"]
    assert observed is not None
    assert observed.rstrip("*") == expected


def test_apply_variants_returns_one_mutant_per_result(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "silent")
    variants = [{"contig": "chr1", "start": start, "ref": ref, "alt": alt}]
    record = _call(vg.varcode_apply_variants_to_transcript, layout, variants, "T1")
    assert record["mutant_transcript"]["transcript_id"] == "T1"
    assert len(record["variants"]) == 1
    assert record["mutant_transcript"]["mutant_protein_sequence"] == layout["t1_protein"]


def test_codon_table_and_selenocysteine_for_annotated_transcript(layout: dict[str, Any]) -> None:
    table = _call(vg.varcode_codon_table_for_transcript, layout, "T1")["codon_table"]
    assert table["ncbi_id"] == 1
    assert "ATG" in table["start_codons"] or "CTG" in table["start_codons"]
    assert _call(vg.varcode_reference_selenocysteine, layout, "T1")["cdna_offsets"] == []


def test_classify_from_protein_diff_on_a_single_substitution(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    protein = layout["t1_protein"]
    index = (start - layout["t1_span"][0]) // 3
    replacement = "W" if protein[index] != "W" else "C"
    mutant = protein[:index] + replacement + protein[index + 1 :]
    record = _call(
        vg.varcode_classify_from_protein_diff,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        "T1",
        protein,
        mutant,
    )
    assert record["effect"]["class"] != "Silent"
    assert record["effect"]["transcript_id"] == "T1"


def test_unknown_transcript_is_a_configuration_error(layout: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        _call(vg.varcode_codon_table_for_transcript, layout, "ENST_NOT_HERE")


# --- random, left-alignment and realised effects ----------------------------------------------


def test_random_variants_are_reproducible_and_valid(layout: dict[str, Any]) -> None:
    first = _call(vg.varcode_random_variants, layout, 4, random_seed=11)
    second = _call(vg.varcode_random_variants, layout, 4, random_seed=11)
    assert first["variants"] == second["variants"]
    assert first["count"] == 4
    for variant in first["variants"]:
        assert variant["contig"] in {"chr1", "chr2"}
        assert variant["ref"] or variant["alt"]


def test_left_alignment_keeps_count_and_never_moves_right(layout: dict[str, Any]) -> None:
    # A one-base deletion in the T1 CDS; its left-aligned position is at or before the input position.
    start = layout["t1_span"][0] + 30
    ref = layout["chr1"][start - 1 : start + 1]
    record = _call(
        vg.varcode_left_align_indels,
        layout,
        [{"contig": "chr1", "start": start, "ref": ref, "alt": ref[0]}],
    )
    assert record["count"] == 1
    assert record["variants"][0]["start"] <= start


def test_realized_effect_reports_top_effect_and_candidates(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    record = _call(
        vg.varcode_predict_realized_effect,
        layout,
        [{"contig": "chr1", "start": start, "ref": ref, "alt": alt}],
        "T1",
    )
    assert record["top_effect"]["class"]
    assert isinstance(record["candidates"], list)


def test_germline_window_contains_the_somatic_position(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    window = _call(
        vg.varcode_default_germline_window,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        "T1",
    )
    assert window["contig"] == "chr1"
    assert window["start"] <= start <= window["end"]


def test_germline_aware_effect_with_a_nearby_germline_variant(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    germline_position = start + 1
    germline = {
        "contig": "chr1",
        "start": germline_position,
        "ref": layout["chr1"][germline_position - 1],
        "alt": "A" if layout["chr1"][germline_position - 1] != "A" else "C",
    }
    record = _call(
        vg.varcode_predict_germline_aware_effect,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        "T1",
        [germline],
    )
    assert record["germline_count"] == 1
    assert record["effect"]["class"]


def test_germline_applied_to_transcript_returns_a_mutant(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "silent")
    germline = {"contig": "chr1", "start": start, "ref": ref, "alt": alt}
    record = _call(vg.varcode_apply_germline_to_transcript, layout, "T1", [germline])
    assert record["mutant_transcript"]["transcript_id"] == "T1"


# --- structural variants, splice layouts and phase -----------------------------------------------


def test_structural_deletion_inside_the_transcript_is_classified(layout: dict[str, Any]) -> None:
    start, end = layout["t1_span"]
    deletion = {"contig": "chr1", "start": start + 30, "end": start + 60, "sv_type": "del"}
    record = _call(vg.varcode_predict_structural_variant_effect, layout, deletion, "T1")
    assert record["structural_variant"]["sv_type"] == "DEL"
    assert record["effect"]["class"] != "Intergenic"


def test_cryptic_exon_enumeration_returns_a_list(layout: dict[str, Any]) -> None:
    start = layout["t1_span"][0]
    record = _call(
        vg.varcode_enumerate_from_structural_variant,
        layout,
        {"contig": "chr1", "start": start + 30, "end": start + 60, "sv_type": "DEL"},
        flank=100,
    )
    assert isinstance(record["candidates"], list)
    assert record["count"] == len(record["candidates"])


def test_structural_variant_rejects_unknown_type(layout: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        _call(
            vg.varcode_predict_structural_variant_effect,
            layout,
            {"contig": "chr1", "start": 10, "end": 20, "sv_type": "XYZ"},
            "T1",
        )


def test_splice_layout_functions_run_on_a_transcript(layout: dict[str, Any]) -> None:
    runs = _call(vg.varcode_build_exon_runs, layout, "T1")
    assert runs["runs"], "the transcript has one exon run"
    plans = _call(vg.varcode_enumerate_layout_splice_plans, layout, "T1", max_plans=8)
    assert plans["count"] == len(plans["plans"]) >= 1
    kept = _call(vg.varcode_realize_exon_path, layout, "T1")
    assert kept["status"] == "completed"
    realised = _call(vg.varcode_realize_splice_plan, layout, "T1", plan_index=0)
    assert realised["plan_index"] == 0
    with pytest.raises(ConfigurationError):
        _call(vg.varcode_realize_splice_plan, layout, "T1", plan_index=50, max_plans=8)


def test_phase_resolver_reads_a_phased_vcf(layout: dict[str, Any], tmp_path: Path) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    vcf = tmp_path / "phased.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "##contig=<ID=chr1,length=120000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t{start + 2}\t.\t{layout['chr1'][start + 1]}\tA\t50\tPASS\t.\tGT:PS\t0|1:1\n",
        encoding="utf-8",
    )
    record = _call(
        vg.varcode_apply_phase_resolver_to_effects,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        str(vcf),
        "S1",
    )
    assert record["sample"] == "S1"
    assert record["count"] == len(record["effects"]) >= 1


def test_candidates_and_rna_outcome_wrap_the_predicted_effect(layout: dict[str, Any]) -> None:
    start, ref, alt = _snv_in_cds(layout, "missense")
    candidates = _call(
        vg.varcode_candidates_from_effects,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        source="unit",
    )
    assert candidates["source"] == "unit"
    assert candidates["count"] == len(candidates["candidates"]) >= 1
    outcome = _call(
        vg.varcode_make_rna_outcome,
        layout,
        {"contig": "chr1", "start": start, "ref": ref, "alt": alt},
        "T1",
        read_count=7,
    )
    assert outcome["candidate"]["source"] == "rna"
    assert outcome["candidate"]["evidence"] == {"read_count": 7}


# --- reference resolution (no network) --------------------------------------------------------


def test_infer_genome_for_an_integer_release_is_lazy_and_local() -> None:
    record = vg.varcode_infer_genome(116, verbose=False)
    assert record["release"] == 116
    assert record["reference_name"] == "GRCh38"
    assert record["species"] == "homo_sapiens"
    assert record["converted_from_ucsc"] is False


def test_infer_genome_functions_pass_through_monkeypatched_resolver(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    varcode = pytest.importorskip("varcode")

    class Species:
        latin_name = "homo_sapiens"

    class Release:
        reference_name = "GRCh37"
        release = 75
        species = Species()

    monkeypatch.setattr(varcode.reference, "infer_genome", lambda value: (Release(), True))
    monkeypatch.setattr(
        varcode.reference, "infer_genome_for_reference_name", lambda value: (Release(), False)
    )
    by_name = vg.varcode_infer_genome("hg19", verbose=False)
    assert by_name["converted_from_ucsc"] is True and by_name["reference_name"] == "GRCh37"
    by_reference = vg.varcode_infer_genome_for_reference_name("GRCh37.75", verbose=False)
    assert by_reference["release"] == 75 and by_reference["converted_from_ucsc"] is False


def test_infer_genome_rejects_booleans_and_empty_names() -> None:
    with pytest.raises(ConfigurationError):
        vg.varcode_infer_genome(True, verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        vg.varcode_infer_genome_for_reference_name("  ", verbose=False)


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live check: set DNAKIT_LIVE_TESTS=1"
)
def test_live_infer_genome_for_reference_name_uses_pyensembl_release_index() -> None:
    record = vg.varcode_infer_genome_for_reference_name("GRCh38", verbose=False)
    assert record["reference_name"] == "GRCh38"
    assert isinstance(record["release"], int)


# --- argument validation and dependency handling ------------------------------------------------


def test_variant_fields_and_bases_are_validated(layout: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        _call(
            vg.varcode_predict_variant_effects,
            layout,
            {"contig": "chr1", "start": 1000, "ref": "A"},
        )
    with pytest.raises(ConfigurationError):
        _call(
            vg.varcode_predict_variant_effects,
            layout,
            {"contig": "chr1", "start": 1000, "ref": "AXG", "alt": "A"},
        )
    with pytest.raises(ConfigurationError):
        _call(
            vg.varcode_predict_variant_effects,
            layout,
            {"contig": "chr1", "start": 0, "ref": "A", "alt": "C"},
        )


def test_annotator_name_must_be_registered(layout: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        _call(
            vg.varcode_predict_variant_effects,
            layout,
            {"contig": "chr1", "start": 5000, "ref": "A", "alt": "C"},
            annotator="not-an-annotator",
        )


def test_missing_varcode_reports_an_install_hint(
    layout: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    def unavailable(module: str, extra: str) -> Any:
        raise BackendUnavailableError(
            f"{module} is not installed", hint=f'python -m pip install "dnakit[{extra}]"'
        )

    monkeypatch.setattr(vg, "require", unavailable)
    with pytest.raises(BackendUnavailableError):
        vg.varcode_reference_base(
            layout["gtf"],
            str(GENOME),
            "chr1",
            10,
            cache_dir=str(layout["cache"] + "_missing"),
            verbose=False,
        )
