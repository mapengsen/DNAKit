"""Tests for dnakit.molecular_lab.varcode_precondition_functions (varcode on a local GTF and genome FASTA).

Every test runs offline on the shared genome (consistency/data/genome.fa). The GTF is written into a temporary
folder and has three genes: G1, a single-exon ORF on the forward strand of chr1; G2, a single-exon ORF on the
reverse strand of chr2; and G3, a two-exon non-coding transcript on chr1 used for the splice case. Observed RNA
rows are built from the genome sequence in these tests. The expected protein and effect classes are derived here,
not copied from varcode output. Nothing is downloaded.
"""

from __future__ import annotations

import csv
import inspect
import itertools
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("varcode")
pytest.importorskip("pyensembl")
pysam = pytest.importorskip("pysam")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import varcode_precondition_functions as vp

GENOME = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data/genome.fa"
)
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


def _spliced_orf(chromosome: str, minimum: int = 2000) -> tuple[tuple[int, int], tuple[int, int]]:
    """Two exons (1-based inclusive) on chr1 whose spliced sequence is an open reading frame.

    Exon 1 is 60 bases starting at an ATG with no in-frame stop; exon 2 starts after an intron of at least 100 bases
    and ends at its first in-frame stop codon, at least 50 codons later. The spliced coding sequence then has no
    internal stop, as a real transcript would.
    """
    for start in range(minimum, len(chromosome) - 2000):
        if chromosome[start : start + 3] != "ATG":
            continue
        exon1 = chromosome[start : start + 60]
        if any(exon1[i : i + 3] in STOPS for i in range(0, 60, 3)):
            continue
        for first in range(start + 160, start + 1500):
            stop = next(
                (i for i in range(0, 3000, 3) if chromosome[first + i : first + i + 3] in STOPS),
                None,
            )
            if stop is not None and stop >= 150:
                return (start + 1, start + 60), (first + 1, first + stop + 3)
    raise AssertionError("no spliced ORF found in the shared genome")


def _revcomp(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


def _translate(coding: str) -> str:
    return "".join(CODON[coding[i : i + 3]] for i in range(0, len(coding) - 2, 3)).rstrip("*")


def _orf(sequence: str) -> tuple[int, int]:
    """First in-frame ATG (0-based) that reaches a stop codon, between 60 and 600 bases; returns [start, end)."""
    for index in range(1000, len(sequence) - 3):
        if sequence[index : index + 3] != "ATG":
            continue
        end = index
        while end + 3 <= len(sequence) and sequence[end : end + 3] not in STOPS:
            end += 3
        if end + 3 <= len(sequence) and 60 <= end - index <= 600:
            return index, end + 3
    raise AssertionError("no ORF found in the shared genome")


def _row(contig: str, strand: str, feature: str, start: int, end: int, attributes: str) -> str:
    frame = "0" if feature == "CDS" else "."
    return "\t".join(
        [contig, "test", feature, str(start), str(end), ".", strand, frame, attributes]
    )


def _gtf(lines_by_gene: list[list[str]]) -> str:
    return "\n".join(line for group in lines_by_gene for line in group) + "\n"


def _gene_lines(
    contig: str,
    strand: str,
    gene_id: str,
    transcript_id: str,
    exons: list[tuple[int, int]],
    cds: list[tuple[int, int]],
) -> list[str]:
    span = (min(s for s, _ in exons), max(e for _, e in exons))
    gene_attr = f'gene_id "{gene_id}"; gene_name "{gene_id}"; gene_biotype "protein_coding";'
    tx_attr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; transcript_biotype "protein_coding";'
    lines = [
        _row(contig, strand, "gene", *span, gene_attr),
        _row(contig, strand, "transcript", *span, tx_attr),
    ]
    for number, (start, end) in enumerate(exons, start=1):
        lines.append(_row(contig, strand, "exon", start, end, f'{tx_attr} exon_number "{number}";'))
    for number, (start, end) in enumerate(cds, start=1):
        attr = f'{tx_attr} exon_number "{number}"; protein_id "P_{transcript_id}";'
        lines.append(_row(contig, strand, "CDS", start, end, attr))
    first, last = cds[0], cds[-1]
    start_codon = (first[0], first[0] + 2) if strand == "+" else (first[1] - 2, first[1])
    stop_codon = (last[1] - 2, last[1]) if strand == "+" else (last[0], last[0] + 2)
    lines.append(_row(contig, strand, "start_codon", *start_codon, f'{tx_attr} exon_number "1";'))
    lines.append(_row(contig, strand, "stop_codon", *stop_codon, f'{tx_attr} exon_number "1";'))
    return lines


@pytest.fixture(scope="module")
def layout(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    if not GENOME.exists():
        pytest.skip("shared consistency genome.fa not present")
    with pysam.FastaFile(str(GENOME)) as fasta:
        chr1 = fasta.fetch("chr1")
        chr2 = fasta.fetch("chr2")
    a, b = _orf(chr1)
    c, d = _orf(_revcomp(chr2))
    length = len(chr2)
    forward = (length - d + 1, length - c)
    splice = _spliced_orf(chr1)
    t1_cds = chr1[a:b]
    t2_cds = _revcomp(chr2[forward[0] - 1 : forward[1]])
    directory = tmp_path_factory.mktemp("varcode_precondition")
    gtf_lines = [
        _gene_lines("chr1", "+", "G1", "T1", [(a + 1, b)], [(a + 1, b)]),
        _gene_lines("chr2", "-", "G2", "T2", [forward], [forward]),
        _gene_lines("chr1", "+", "G3", "T3", list(splice), list(splice)),
    ]
    gtf = directory / "genes.gtf"
    gtf.write_text(_gtf(gtf_lines), encoding="utf-8")
    return {
        "gtf": str(gtf),
        "cache": str(directory / "cache"),
        "directory": directory,
        "chr1": chr1,
        "chr2": chr2,
        "a": a,
        "b": b,
        "forward": forward,
        "splice": splice,
        "t1_cds": t1_cds,
        "t2_cds": t2_cds,
        "t1_protein": _translate(t1_cds),
    }


def _call(function: Callable[..., Any], layout: dict[str, Any], *args: Any, **kwargs: Any) -> Any:
    kwargs.setdefault("verbose", False)
    return function(layout["gtf"], str(GENOME), *args, cache_dir=layout["cache"], **kwargs)


def _sv(layout: dict[str, Any]) -> dict[str, Any]:
    """A 21-base deletion inside the T1 exon, used as the DNA call of each fusion."""
    return {"contig": "chr1", "start": layout["a"] + 20, "end": layout["a"] + 40, "sv_type": "DEL"}


def _write_exacto(
    layout: dict[str, Any], *, five_offset: int = 10, structures: bool = True
) -> tuple[Path, Path]:
    """Write ExACTO structure and integration TSVs for one fusion from T1 (5') to T2 (3')."""
    a = layout["a"]
    chr1, chr2 = layout["chr1"], layout["chr2"]
    lo1 = a + 1 + five_offset
    hi1 = lo1 + 59
    seq1 = chr1[lo1 - 1 : hi1]
    hi2 = layout["forward"][1]
    lo2 = hi2 - 59
    seq2 = _revcomp(chr2[lo2 - 1 : hi2])
    fields = [
        "transcript_model_id", "reference_transcript_ids", "index", "read_start", "read_end",
        "sequence", "type", "kind", "context", "chromosome_1", "position_1", "strand_1",
        "chromosome_2", "position_2", "strand_2", "transcript_id_1", "transcript_id_2",
    ]  # fmt: skip
    base = {
        "transcript_model_id": "M1", "reference_transcript_ids": "T1,T2", "type": "base",
        "kind": "match", "context": "linear",
    }  # fmt: skip
    rows = [
        {**base, "index": 0, "read_start": 0, "read_end": 59, "sequence": seq1,
         "chromosome_1": "chr1", "position_1": lo1, "strand_1": "+", "chromosome_2": "chr1",
         "position_2": hi1, "strand_2": "+", "transcript_id_1": "T1", "transcript_id_2": "T1"},
        {**base, "index": 1, "read_start": 60, "read_end": 119, "sequence": seq2,
         "chromosome_1": "chr2", "position_1": hi2, "strand_1": "-", "chromosome_2": "chr2",
         "position_2": lo2, "strand_2": "-", "transcript_id_1": "T2", "transcript_id_2": "T2"},
    ]  # fmt: skip
    structures_path = layout["directory"] / f"structures_{five_offset}.tsv"
    integrated_path = layout["directory"] / f"integrated_{five_offset}.tsv"
    with structures_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        if structures:
            writer.writerows(rows)
    with integrated_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "transcript_model_id",
                "reference_transcript_ids",
                "rna_variant_call_id",
                "dna_variant_call_id",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"transcript_model_id": "M1", "reference_transcript_ids": "T1,T2",
             "rna_variant_call_id": "R1", "dna_variant_call_id": "D1"}
        )  # fmt: skip
    return structures_path, integrated_path


def test_public_functions_are_keyword_verbose() -> None:
    for name in (
        "varcode_make_fusion_outcome",
        "varcode_apply_rna_evidence_to_effects",
        "varcode_load_exacto_fusions",
        "varcode_reference_completion_hypotheses",
        "varcode_enumerate_splice_outcomes",
    ):
        parameter = inspect.signature(getattr(vp, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_missing_varcode_raises_backend_unavailable(
    layout: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "varcode", None)
    with pytest.raises(BackendUnavailableError, match="bio-varcode"):
        _call(vp.varcode_make_fusion_outcome, layout, _sv(layout), "T1", "ACGT", "M1")


def test_make_fusion_outcome_with_partner_is_gene_fusion(layout: dict[str, Any]) -> None:
    sequence = layout["t1_cds"][:60] + layout["t2_cds"][:60]
    result = _call(
        vp.varcode_make_fusion_outcome,
        layout,
        _sv(layout),
        "T1",
        sequence,
        "M1",
        partner_transcript_id="T2",
        read_count=7,
    )
    candidate = result["candidate"]
    assert candidate["effect"]["class"] == "GeneFusion"
    assert candidate["source"] == "rna"
    assert candidate["evidence"]["transcript_model_id"] == "M1"
    assert candidate["evidence"]["read_count"] == 7
    assert result["mutant_transcript"]["cdna_sequence"] == sequence
    assert result["mutant_transcript"]["mutant_protein_sequence"] is None


def test_make_fusion_outcome_without_partner_is_intergenic(layout: dict[str, Any]) -> None:
    result = _call(
        vp.varcode_make_fusion_outcome,
        layout,
        _sv(layout),
        "T1",
        layout["t1_cds"][:90],
        "M2",
    )
    assert result["candidate"]["effect"]["class"] == "TranslocationToIntergenic"


def test_make_fusion_outcome_translates_complete_orf(layout: dict[str, Any]) -> None:
    sequence = layout["t1_cds"] + layout["t2_cds"][:60]
    result = _call(
        vp.varcode_make_fusion_outcome,
        layout,
        _sv(layout),
        "T1",
        sequence,
        "M3",
        partner_transcript_id="T2",
        cds_start=0,
    )
    mutant = result["mutant_transcript"]
    assert mutant["mutant_protein_sequence"] == layout["t1_protein"]
    assert result["candidate"]["evidence"]["protein_status"] == "predicted_from_observed_rna"


def test_make_fusion_outcome_rejects_bad_input(layout: dict[str, Any]) -> None:
    sv = _sv(layout)
    with pytest.raises(ConfigurationError, match="read_count"):
        _call(vp.varcode_make_fusion_outcome, layout, sv, "T1", "ACGT", "M1", read_count=-1)
    with pytest.raises(ConfigurationError, match="start codon"):
        _call(vp.varcode_make_fusion_outcome, layout, sv, "T1", "ACGT" * 10, "M1", cds_start=0)
    with pytest.raises(ConfigurationError, match="sv_type"):
        _call(
            vp.varcode_make_fusion_outcome,
            layout,
            {"contig": "chr1", "start": 10, "end": 20, "sv_type": "TRA"},
            "T1",
            "ACGT",
            "M1",
        )
    with pytest.raises(ConfigurationError, match="only A, C, G, T or N"):
        _call(vp.varcode_make_fusion_outcome, layout, sv, "T1", "ACGTX", "M1")


def test_apply_rna_evidence_attaches_observed_candidate(layout: dict[str, Any]) -> None:
    sv = _sv(layout)
    fusion = {
        "transcript_id": "T1",
        "partner_transcript_id": "T2",
        "transcript_model_id": "M4",
        "sequence": layout["t1_cds"][:60] + layout["t2_cds"][:60],
        "read_count": 9,
    }
    result = _call(vp.varcode_apply_rna_evidence_to_effects, layout, sv, [fusion])
    assert result["fusion_count"] == 1
    assert result["observed_candidate_count"] >= 1
    observed = [c for row in result["effects"] for c in row["candidates"]]
    assert any(c["evidence"].get("transcript_model_id") == "M4" for c in observed)
    assert all(
        c["evidence"].get("read_count") == 9
        for c in observed
        if c["evidence"].get("transcript_model_id") == "M4"
    )


def test_apply_rna_evidence_rejects_empty_list(layout: dict[str, Any]) -> None:
    sv = _sv(layout)
    with pytest.raises(ConfigurationError, match="fusion_outcomes"):
        _call(vp.varcode_apply_rna_evidence_to_effects, layout, sv, [])


def test_load_exacto_fusions_selects_linked_call(layout: dict[str, Any]) -> None:
    structures, integrated = _write_exacto(layout, five_offset=0)
    sv = {"D1": _sv(layout)}
    result = _call(
        vp.varcode_load_exacto_fusions,
        layout,
        str(structures),
        str(integrated),
        sv,
    )
    assert result["count"] == 1
    candidate = result["candidates"][0]
    assert candidate["source"] == "exacto"
    assert candidate["effect"]["class"] == "GeneFusion"
    assert candidate["evidence"]["partner_status"] == "sense"
    assert candidate["evidence"]["dna_variant_call_id"] == "D1"


def test_load_exacto_fusions_requires_structures(layout: dict[str, Any]) -> None:
    structures, integrated = _write_exacto(layout, five_offset=0, structures=False)
    sv = {"D1": _sv(layout)}
    with pytest.raises(ConfigurationError, match="Missing structures"):
        _call(vp.varcode_load_exacto_fusions, layout, str(structures), str(integrated), sv)


def test_reference_completion_extends_five_prime_end(layout: dict[str, Any]) -> None:
    structures, integrated = _write_exacto(layout, five_offset=10)
    sv = {"D1": _sv(layout)}
    result = _call(
        vp.varcode_reference_completion_hypotheses,
        layout,
        str(structures),
        str(integrated),
        sv,
        ends=["five_prime"],
    )
    assert result["observed_count"] == 1
    assert result["hypothesis_count"] >= 1
    hypothesis = result["hypotheses"][0]
    assert hypothesis["source"] == "varcode_reference_completion"
    assert hypothesis["evidence"]["sequence_status"] == "reference_completion_hypothesis"
    assert hypothesis["effect"]["class"] == "GeneFusion"


def test_reference_completion_rejects_unknown_end(layout: dict[str, Any]) -> None:
    structures, integrated = _write_exacto(layout)
    sv = {"D1": _sv(layout)}
    with pytest.raises(ConfigurationError, match="ends"):
        _call(
            vp.varcode_reference_completion_hypotheses,
            layout,
            str(structures),
            str(integrated),
            sv,
            ends=["middle"],
        )


def test_enumerate_splice_outcomes_on_multi_exon_donor(layout: dict[str, Any]) -> None:
    donor_intron_base = layout["splice"][0][1] + 2  # second intron base of the donor (1-based)
    ref = layout["chr1"][donor_intron_base - 1]
    alt = "A" if ref != "A" else "C"
    variant = {"contig": "chr1", "start": donor_intron_base, "ref": ref, "alt": alt}
    result = _call(vp.varcode_enumerate_splice_outcomes, layout, variant, "T3")
    assert result["input_effect"]["class"] == "SpliceOutcomeSet"
    assert result["is_outcome_set"] is True
    assert result["genomic_sequence_used"] is True
    classes = {candidate["effect"]["class"] for candidate in result["candidates"]}
    assert {"NormalSplicing", "IntronRetention", "CrypticDonor"} <= classes
    assert result["count"] == len(result["candidates"])


def test_enumerate_splice_outcomes_passes_non_splice_effect_through(layout: dict[str, Any]) -> None:
    position = layout["a"] + 30
    ref = layout["chr1"][position - 1]
    alt = "A" if ref != "A" else "C"
    variant = {"contig": "chr1", "start": position, "ref": ref, "alt": alt}
    result = _call(
        vp.varcode_enumerate_splice_outcomes,
        layout,
        variant,
        "T1",
        use_genomic_sequence=False,
    )
    assert result["is_outcome_set"] is False
    assert result["count"] == 0
    assert result["genomic_sequence_used"] is False
