"""Fixed 240-field DNA descriptor schema with formulas, units, and sources."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from dnakit.descriptors._dinucleotide import DINUCLEOTIDE_PROPERTY_SPECS

DESCRIPTOR_SCHEMA_VERSION = "descriptor_schema_v1"


@dataclass(frozen=True)
class DescriptorField:
    """Metadata for one ordered field in the 240-value descriptor vector."""

    index: int
    name: str
    category: str
    unit: str
    formula: str
    source: str

    def to_dict(self) -> dict[str, int | str]:
        return {
            "index": self.index,
            "name": self.name,
            "category": self.category,
            "unit": self.unit,
            "formula": self.formula,
            "source": self.source,
        }


def _build_schema() -> tuple[DescriptorField, ...]:
    fields: list[DescriptorField] = []

    def add(name: str, category: str, unit: str, formula: str, source: str) -> None:
        fields.append(DescriptorField(len(fields) + 1, name, category, unit, formula, source))

    native = "DNAKit descriptor_schema_v1"
    for name, unit, formula in (
        ("symbol_length", "nt", "number of nucleotide symbols; explicit gaps excluded"),
        (
            "coordinate_span",
            "nt",
            "symbol_length + sum(known gap lengths); undefined with an unknown gap",
        ),
        ("canonical_base_count", "count", "count(A,C,G,T)"),
        ("ambiguity_symbol_count", "count", "symbol_length - canonical_base_count"),
        ("gap_object_count", "count", "number of explicit Gap objects"),
        ("known_gap_nt", "nt", "sum(length of known explicit gaps)"),
        ("unknown_gap_count", "count", "number of explicit gaps with unknown length"),
        ("canonical_symbol_fraction", "fraction", "canonical_base_count / symbol_length"),
        ("ambiguity_symbol_fraction", "fraction", "ambiguity_symbol_count / symbol_length"),
        ("known_gap_fraction", "fraction", "known_gap_nt / coordinate_span"),
        (
            "canonical_run_count",
            "count",
            "number of uninterrupted A/C/G/T runs split by ambiguity or gaps",
        ),
        (
            "longest_canonical_run_nt",
            "nt",
            "maximum uninterrupted A/C/G/T run length",
        ),
    ):
        add(name, "basic", unit, formula, native)

    composition_fields = (
        ("purine_count", "count", "count(A)+count(G)"),
        ("purine_fraction", "fraction", "purine_count / canonical_base_count"),
        ("pyrimidine_count", "count", "count(C)+count(T)"),
        (
            "pyrimidine_fraction",
            "fraction",
            "pyrimidine_count / canonical_base_count",
        ),
        ("amino_count", "count", "count(A)+count(C)"),
        ("amino_fraction", "fraction", "amino_count / canonical_base_count"),
        ("keto_count", "count", "count(G)+count(T)"),
        ("keto_fraction", "fraction", "keto_count / canonical_base_count"),
        ("weak_count", "count", "count(A)+count(T)"),
        ("weak_fraction", "fraction", "weak_count / canonical_base_count"),
        ("strong_count", "count", "count(C)+count(G)"),
        ("strong_fraction", "fraction", "strong_count / canonical_base_count"),
        (
            "purine_pyrimidine_skew",
            "ratio",
            "(purine_count-pyrimidine_count)/(purine_count+pyrimidine_count)",
        ),
        (
            "amino_keto_skew",
            "ratio",
            "(amino_count-keto_count)/(amino_count+keto_count)",
        ),
        (
            "weak_strong_skew",
            "ratio",
            "(weak_count-strong_count)/(weak_count+strong_count)",
        ),
        ("gc_at_ratio", "ratio", "strong_count / weak_count"),
    )
    for name, unit, formula in composition_fields:
        add(name, "composition", unit, formula, native)

    for k in range(1, 4):
        for symbols in product("ACGT", repeat=k):
            word = "".join(symbols)
            add(
                f"k{k}_{word}_frequency",
                "kmer",
                "fraction",
                f"overlapping count({word}) / valid canonical {k}-mer positions",
                "DNAKit exact overlapping k-mer definition",
            )

    skew_cpg_fields = (
        ("gc_skew", "ratio", "(count(G)-count(C))/(count(G)+count(C))"),
        ("at_skew", "ratio", "(count(A)-count(T))/(count(A)+count(T))"),
        ("cpg_count", "count", "overlapping count(CG)"),
        ("cpg_density", "fraction", "count(CG) / valid canonical dinucleotide positions"),
        (
            "cpg_observed_expected",
            "ratio",
            "count(CG)*canonical_base_count/(count(C)*count(G))",
        ),
        ("gpc_count", "count", "overlapping count(GC)"),
        ("gpc_density", "fraction", "count(GC) / valid canonical dinucleotide positions"),
        ("cpg_gpc_ratio", "ratio", "count(CG) / count(GC)"),
        (
            "cumulative_gc_skew_max",
            "count",
            "max prefix cumulative score where G=+1,C=-1,A/T=0",
        ),
        (
            "cumulative_gc_skew_min",
            "count",
            "min prefix cumulative score where G=+1,C=-1,A/T=0",
        ),
        (
            "cumulative_gc_skew_range",
            "count",
            "cumulative_gc_skew_max - cumulative_gc_skew_min",
        ),
        (
            "cumulative_at_skew_max",
            "count",
            "max prefix cumulative score where A=+1,T=-1,C/G=0",
        ),
        (
            "cumulative_at_skew_min",
            "count",
            "min prefix cumulative score where A=+1,T=-1,C/G=0",
        ),
        (
            "cumulative_at_skew_range",
            "count",
            "cumulative_at_skew_max - cumulative_at_skew_min",
        ),
        (
            "dinucleotide_rc_total_variation",
            "fraction",
            "0.5*sum_xy(abs(f_xy-f_reverse_complement(xy)))",
        ),
        (
            "mono_chargaff_l1_distance",
            "fraction",
            "abs(f_A-f_T)+abs(f_C-f_G)",
        ),
    )
    for name, unit, formula in skew_cpg_fields:
        add(name, "skew_cpg", unit, formula, native)

    shannon = "Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x"
    for k in range(1, 4):
        add(
            f"shannon_entropy_k{k}_bits",
            "complexity",
            "bits",
            f"-sum(p({k}-mer)*log2(p({k}-mer)))",
            shannon,
        )
    for k in range(1, 4):
        add(
            f"normalized_entropy_k{k}",
            "complexity",
            "fraction",
            f"shannon_entropy_k{k}_bits / log2(4**{k})",
            shannon,
        )
    for k in range(2, 7):
        add(
            f"linguistic_complexity_k{k}",
            "complexity",
            "fraction",
            f"unique {k}-mers / min(4**{k}, valid {k}-mer positions)",
            "Observed/possible k-word linguistic complexity",
        )
    add(
        "linguistic_complexity_product_k1_k6",
        "complexity",
        "fraction",
        "product of defined linguistic_complexity_k values for k=1..6",
        "Observed/possible k-word linguistic complexity",
    )
    add(
        "lz76_complexity",
        "complexity",
        "count",
        "number of phrases in exhaustive Lempel-Ziv 1976 parsing",
        "Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501",
    )
    add(
        "normalized_lz76_complexity",
        "complexity",
        "ratio",
        "lz76_complexity*log_base4(canonical_base_count)/canonical_base_count",
        "Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501",
    )
    add("longest_homopolymer_nt", "complexity", "nt", "max canonical homopolymer run", native)
    for base in "ACGT":
        add(
            f"longest_homopolymer_{base.lower()}_nt",
            "complexity",
            "nt",
            f"max homopolymer run of {base}",
            native,
        )
    add(
        "exact_tandem_repeat_coverage_fraction",
        "complexity",
        "fraction",
        "union bases covered by exact tandem repeats / canonical_base_count",
        "DNAKit exact tandem repeat scanner; units 1..20; minimum repeats 2",
    )

    coding_fields = (
        ("frame0_codon_count", "count", "number of valid forward frame-0 codons"),
        ("frame0_unique_codon_count", "count", "number of distinct forward frame-0 codons"),
        ("frame0_start_codon_count", "count", "count(ATG) in forward frame 0"),
        (
            "frame0_stop_codon_count",
            "count",
            "count(TAA,TAG,TGA) in forward frame 0",
        ),
        (
            "frame0_start_codon_fraction",
            "fraction",
            "frame0_start_codon_count / frame0_codon_count",
        ),
        (
            "frame0_stop_codon_fraction",
            "fraction",
            "frame0_stop_codon_count / frame0_codon_count",
        ),
        (
            "frame0_codon_entropy_bits",
            "bits",
            "-sum(frame0 codon frequency*log2(frequency))",
        ),
        (
            "frame0_effective_number_of_codons",
            "count",
            "2**frame0_codon_entropy_bits",
        ),
        ("frame0_gc1_fraction", "fraction", "GC bases at position 1 / valid frame-0 codons"),
        ("frame0_gc2_fraction", "fraction", "GC bases at position 2 / valid frame-0 codons"),
        ("frame0_gc3_fraction", "fraction", "GC bases at position 3 / valid frame-0 codons"),
        (
            "six_frame_complete_orf_count",
            "count",
            "complete start-to-next-stop ORFs across three frames on both strands",
        ),
        (
            "six_frame_forward_complete_orf_count",
            "count",
            "complete ORFs across three forward frames",
        ),
        (
            "six_frame_reverse_complete_orf_count",
            "count",
            "complete ORFs across three reverse-complement frames",
        ),
        (
            "six_frame_longest_complete_orf_nt",
            "nt",
            "maximum complete six-frame ORF length including terminal stop",
        ),
        (
            "six_frame_complete_orf_coverage_fraction",
            "fraction",
            "union of symbol positions covered by complete six-frame ORFs / canonical_base_count",
        ),
    )
    for name, unit, formula in coding_fields:
        add(name, "coding", unit, formula, "NCBI standard genetic code table 1; DNAKit ORF rules")

    physchem_fields = (
        ("mw_ss_oh_da", "Da", "anhydrous mass of one ssDNA strand with 5-prime OH"),
        (
            "mw_ss_5p_phosphate_da",
            "Da",
            "anhydrous mass of one ssDNA strand with 5-prime phosphate",
        ),
        (
            "mw_ds_oh_da",
            "Da",
            "anhydrous mass of sequence plus complete reverse complement with 5-prime OH",
        ),
        (
            "mw_ds_5p_phosphate_da",
            "Da",
            "anhydrous mass of sequence plus complete reverse complement; "
            "both 5-prime phosphorylated",
        ),
        (
            "epsilon260_ss_m_inverse_cm_inverse",
            "M^-1 cm^-1",
            "nearest-neighbor epsilon260 pair sum minus internal-base sum",
        ),
        (
            "nmol_per_a260_1ml_1cm",
            "nmol",
            "1e6 / epsilon260 for A260=1, volume=1 mL, path=1 cm",
        ),
        (
            "ug_per_a260_1ml_1cm",
            "ug",
            "1000*mw_ss_oh_da/epsilon260 for A260=1, volume=1 mL, path=1 cm",
        ),
        ("tm_wallace_c", "degree C", "2*(A+T)+4*(G+C)"),
        (
            "stacking_delta_h_kcal_per_mol",
            "kcal/mol",
            "sum SantaLucia nearest-neighbor stacking delta H",
        ),
        (
            "stacking_delta_s_cal_per_k_mol_k",
            "cal/(K mol)",
            "sum SantaLucia nearest-neighbor stacking delta S",
        ),
        (
            "stacking_delta_g37_kcal_per_mol",
            "kcal/mol",
            "stacking delta H - 310.15*stacking delta S/1000",
        ),
        (
            "nn_delta_h_kcal_per_mol",
            "kcal/mol",
            "SantaLucia complete-duplex delta H with initiation and symmetry",
        ),
        (
            "nn_delta_s_cal_per_mol_k",
            "cal/(K mol)",
            "SantaLucia complete-duplex delta S with initiation, symmetry, and salt",
        ),
        (
            "nn_delta_g37_kcal_per_mol",
            "kcal/mol",
            "complete-duplex delta H - 310.15*delta S/1000",
        ),
        ("nn_tm_c", "degree C", "SantaLucia concentration- and sodium-adjusted Tm"),
        ("self_complementary", "boolean", "sequence == reverse_complement(sequence)"),
    )
    for name, unit, formula in physchem_fields:
        if name.startswith("mw_"):
            source = "DNAKit anhydrous DNA residue mass table v1"
        elif "epsilon260" in name or "a260" in name:
            source = "Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970"
        elif name == "tm_wallace_c":
            source = "Wallace short-oligo 2AT+4GC rule"
        else:
            source = "SantaLucia 1998; DOI 10.1073/pnas.95.4.1460"
        add(name, "physicochemical", unit, formula, source)

    for item in DINUCLEOTIDE_PROPERTY_SPECS:
        for statistic in ("mean", "sd", "min", "max"):
            formula = (
                f"population {statistic} of {item.display_name} values over valid "
                "overlapping dinucleotides"
            )
            add(
                f"diprodb_{item.key}_{statistic}",
                "dinucleotide_property",
                item.unit,
                formula,
                "Caller-supplied table; DNAKit bundles no numerical values",
            )

    if len(fields) != 240:
        raise AssertionError(f"Descriptor schema v1 must contain 240 fields, got {len(fields)}.")
    if len({field.name for field in fields}) != len(fields):
        raise AssertionError("Descriptor schema v1 field names must be unique.")
    return tuple(fields)


DESCRIPTOR_SCHEMA_V1 = _build_schema()
DESCRIPTOR_NAMES_V1 = tuple(field.name for field in DESCRIPTOR_SCHEMA_V1)


def descriptor_schema_v1() -> tuple[DescriptorField, ...]:
    """Return the immutable, ordered 240-field schema."""

    return DESCRIPTOR_SCHEMA_V1


__all__ = [
    "DESCRIPTOR_NAMES_V1",
    "DESCRIPTOR_SCHEMA_V1",
    "DESCRIPTOR_SCHEMA_VERSION",
    "DescriptorField",
    "descriptor_schema_v1",
]
