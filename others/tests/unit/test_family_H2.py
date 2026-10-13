"""Contract tests for the H2 family check (pydna, DnaCauldron, DnaChisel and primer3).

New wrappers are compared with the upstream call or a hand calculation. Kept
``family-variant`` rows are checked by calling the claimed parameter of the existing
wrapper and comparing with the upstream variant. ``duplicate-of`` rows are checked by
object identity or by comparing the result of the same code. All tests run offline.
"""

from __future__ import annotations

import importlib
import math
import random
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import dnacauldron_family_functions as dcf
from dnakit.molecular_lab import dnachisel_family_functions as dchf
from dnakit.molecular_lab import dnachisel_functions as dch
from dnakit.molecular_lab import pydna_family_functions as pyf
from dnakit.molecular_lab import pydna_functions as pyd
from dnakit.molecular_lab import pydna_review_functions as pyr

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _rc(text: str) -> str:
    return text.translate(_COMPLEMENT)[::-1]


def _random_dna(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _resolve(qualname: str) -> object:
    parts = qualname.split(".")
    for cut in range(len(parts), 0, -1):
        try:
            obj: object = importlib.import_module(".".join(parts[:cut]))
        except ImportError:
            continue
        for name in parts[cut:]:
            obj = getattr(obj, name)
        return obj
    raise AssertionError(qualname)


TEMPLATE = _random_dna(3, 300)
FORWARD = TEMPLATE[50:75]
REVERSE = _rc(TEMPLATE[200:225])


# ---------------------------------------------------------------------------
# New pydna wrappers
# ---------------------------------------------------------------------------


def test_common_substrings_with_topology_wraps_the_origin_of_a_circular_molecule() -> None:
    pytest.importorskip("pydna")
    x, y = "ACGTTTGCAAC", "TTGCAACACGTT"
    circular = pyf.find_common_substrings_with_topology(x, y, min_length=6, circular_x=True)
    # Hand calculation: x read from index 4 is TTGCAAC + ACGT = TTGCAACACGT, which opens y.
    assert circular["matches"] == [{"x_start": 4, "y_start": 0, "length": 11}]
    linear = pyf.find_common_substrings_with_topology(x, y, min_length=6)
    assert linear["matches"] == [{"x_start": 4, "y_start": 0, "length": 7}]
    plain = pyd.find_common_substrings(x, y, min_length=6)
    assert linear["matches"] == plain["matches"]
    assert linear["match_count"] == plain["match_count"]
    with pytest.raises(ConfigurationError):
        pyf.find_common_substrings_with_topology(x, y, circular_x="yes")  # type: ignore[arg-type]


def test_common_substrings_with_topology_matches_pydna_docstring() -> None:
    pytest.importorskip("pydna")
    result = pyf.find_common_substrings_with_topology("TAAAAAAT", "CCaAaAaACC", min_length=5)
    triples = [(row["x_start"], row["y_start"], row["length"]) for row in result["matches"]]
    assert triples == [(1, 2, 6), (1, 3, 5), (2, 2, 5)]


def test_simulate_pcr_assembly_matches_pydna_and_hand_slice() -> None:
    pytest.importorskip("pydna")
    from pydna.assembly2 import pcr_assembly
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    result = pyf.simulate_pcr_assembly(TEMPLATE, FORWARD, REVERSE)
    assert result["product_count"] == 1
    assert result["products"][0]["sequence"] == TEMPLATE[50:225]
    upstream = pcr_assembly(Dseqrecord(TEMPLATE), Primer(FORWARD), Primer(REVERSE), limit=14)
    assert [row["sequence"] for row in result["products"]] == [str(p.seq).upper() for p in upstream]
    assert result["products"][0]["circular"] is False
    amplicon = pyd.simulate_pcr_amplicon(TEMPLATE, FORWARD, REVERSE)
    assert amplicon["product_sequence"] == result["products"][0]["sequence"]


def test_simulate_pcr_assembly_tolerates_primer_mismatches() -> None:
    pytest.importorskip("pydna")
    from pydna.assembly2 import pcr_assembly
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    changed = FORWARD[:-1] + ("A" if FORWARD[-1] != "A" else "C")
    strict = pyf.simulate_pcr_assembly(TEMPLATE, changed, REVERSE)
    tolerant = pyf.simulate_pcr_assembly(TEMPLATE, changed, REVERSE, mismatches=1)
    assert strict["product_count"] == 0
    upstream = pcr_assembly(
        Dseqrecord(TEMPLATE), Primer(changed), Primer(REVERSE), limit=14, mismatches=1
    )
    assert tolerant["product_count"] == len(upstream) >= 1
    assert [row["length"] for row in tolerant["products"]] == sorted(len(p) for p in upstream)
    with pytest.raises(ConfigurationError):
        pyf.simulate_pcr_assembly(TEMPLATE, FORWARD, REVERSE, mismatches=-1)
    with pytest.raises(ConfigurationError):
        pyf.simulate_pcr_assembly(TEMPLATE, "ACGN", REVERSE)


def test_simulate_pcr_assembly_circular_template_spans_the_origin() -> None:
    pytest.importorskip("pydna")
    circle = TEMPLATE
    forward = circle[-20:]
    reverse = _rc(circle[10:30])
    result = pyf.simulate_pcr_assembly(circle, forward, reverse, template_circular=True)
    assert [row["sequence"] for row in result["products"]] == [circle[-20:] + circle[:30]]
    assert pyf.simulate_pcr_assembly(circle, forward, reverse)["product_count"] == 0


def test_primer_melting_temperature_dbd_matches_pydna_and_differs_from_taq() -> None:
    pytest.importorskip("pydna")
    from pydna.tm import tm_dbd

    result = pyf.primer_melting_temperature_dbd(FORWARD)
    assert result["tm_celsius"] == round(float(tm_dbd(FORWARD)), 3)
    assert result["tm_celsius"] != pyd.primer_melting_temperature(FORWARD)["tm_celsius"]
    with pytest.raises(ConfigurationError):
        pyf.primer_melting_temperature_dbd("ACGT")
    with pytest.raises(ConfigurationError):
        pyf.primer_melting_temperature_dbd("ACGTNNACGT")


def test_primer_melting_temperature_breslauer_matches_pydna_and_concentration_effect() -> None:
    pytest.importorskip("pydna")
    from pydna.tm import tmbresluc

    result = pyf.primer_melting_temperature_breslauer(FORWARD)
    assert result["tm_celsius"] == round(float(tmbresluc(FORWARD)), 3)
    more_salt = pyf.primer_melting_temperature_breslauer(FORWARD, salt_mm=100.0)
    # Hand calculation: doubling the salt adds 16.6 * log10(2) = 5.0 degrees to the formula.
    assert more_salt["tm_celsius"] - result["tm_celsius"] == pytest.approx(
        16.6 * math.log10(2), abs=0.002
    )
    with pytest.raises(ConfigurationError):
        pyf.primer_melting_temperature_breslauer(FORWARD, primer_nm=0)
    with pytest.raises(ConfigurationError):
        pyf.primer_melting_temperature_breslauer(FORWARD, salt_mm=True)  # type: ignore[arg-type]


def test_amplicon_melting_temperature_hand_calculation() -> None:
    pytest.importorskip("pydna")
    product = TEMPLATE[50:225]
    gc = sum(base in "GC" for base in product) / len(product)
    expected = 81.5 + 0.41 * gc * 100 + 16.6 * math.log10(0.05) - 675 / len(product)
    result = pyf.amplicon_melting_temperature(product)
    assert result["tm_celsius"] == round(expected, 3)
    assert result["gc_fraction"] == round(gc, 6)
    higher = pyf.amplicon_melting_temperature(product, potassium_m=0.5)
    assert higher["tm_celsius"] - result["tm_celsius"] == pytest.approx(16.6, abs=0.002)
    with pytest.raises(ConfigurationError):
        pyf.amplicon_melting_temperature("ACGT")
    with pytest.raises(ConfigurationError):
        pyf.amplicon_melting_temperature(product, potassium_m=-1.0)


def test_pcr_annealing_temperature_taq_and_dbd_follow_the_pydna_rules() -> None:
    pytest.importorskip("pydna")
    from pydna.tm import ta_dbd, ta_default, tm_dbd, tm_default, tm_product

    product = TEMPLATE[50:225]
    taq = pyf.pcr_annealing_temperature(FORWARD, REVERSE, product)
    expected = (
        0.3 * min(tm_default(FORWARD), tm_default(REVERSE)) + 0.7 * tm_product(product) - 14.9
    )
    assert taq["annealing_celsius"] == round(expected, 3)
    assert taq["annealing_celsius"] == round(ta_default(FORWARD, REVERSE, product), 3)
    dbd = pyf.pcr_annealing_temperature(FORWARD, REVERSE, product, polymerase="dbd")
    assert dbd["annealing_celsius"] == round(min(tm_dbd(FORWARD), tm_dbd(REVERSE)) + 3, 3)
    assert dbd["annealing_celsius"] == round(ta_dbd(FORWARD, REVERSE, product), 3)
    assert dbd["forward_tm_celsius"] == round(tm_dbd(FORWARD), 3)
    with pytest.raises(ConfigurationError):
        pyf.pcr_annealing_temperature(FORWARD, REVERSE, product, polymerase="q5")


# ---------------------------------------------------------------------------
# New DnaCauldron wrappers
# ---------------------------------------------------------------------------


def test_read_dnacauldron_records_reads_multi_record_files_folders_and_zip(tmp_path: Path) -> None:
    pytest.importorskip("dnacauldron")
    from dnacauldron.biotools.sequence_io import load_records_from_file, load_records_from_files

    first = tmp_path / "a.fa"
    first.write_text(">partA\nATGCGTACGTTAGC\n>partB\nGGCATCGATCGGA\n", encoding="utf-8")
    second = tmp_path / "b.fasta"
    second.write_text(">\nACGTACGTAGCTAGCTAGC\n", encoding="utf-8")
    result = dcf.read_dnacauldron_records([first, str(second)])
    upstream = load_records_from_files(files=[str(first), str(second)])
    assert [row["id"] for row in result["records"]] == [rec.id for rec in upstream]
    assert [row["sequence"] for row in result["records"]] == [str(rec.seq) for rec in upstream]
    assert [row["id"] for row in result["records"]] == ["partA", "partB", "b"]
    assert [item["format"] for item in result["files"]] == [
        load_records_from_file(str(first))[1],
        "fasta",
    ]
    assert [item["record_count"] for item in result["files"]] == [2, 1]
    folder = dcf.read_dnacauldron_records(folder=tmp_path)
    assert [row["id"] for row in folder["records"]] == ["partA", "partB", "b"]
    named = dcf.read_dnacauldron_records([second], use_file_names_as_ids=True)
    assert named["records"][0]["id"] == "b"
    # read_dnacauldron_parts reads one record per file and so cannot read a two-record file.
    with pytest.raises(ConfigurationError):
        dcf.read_dnacauldron_records([])
    with pytest.raises(ConfigurationError):
        dcf.read_dnacauldron_records([tmp_path / "missing.fa"])
    with pytest.raises(ConfigurationError):
        dcf.read_dnacauldron_records([first], max_records=1)


def test_read_dnacauldron_records_opens_a_zip_archive(tmp_path: Path) -> None:
    pytest.importorskip("dnacauldron")
    import zipfile

    archive = tmp_path / "parts.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("one.fa", ">one\nACGTACGTACGT\n")
    result = dcf.read_dnacauldron_records([archive])
    assert result["files"][0]["format"] == "zip"
    assert [row["sequence"] for row in result["records"]] == ["ACGTACGTACGT"]


def test_parse_dnacauldron_sequence_text_detects_format_like_string_to_records() -> None:
    pytest.importorskip("dnacauldron")
    from dnacauldron.biotools.sequence_io import string_to_records

    bare = dcf.parse_dnacauldron_sequence_text("ACGTACGT")
    assert bare["format"] == "ATGC"
    assert bare["records"][0]["sequence"] == "ACGTACGT"
    assert string_to_records("ACGTACGT")[1] == bare["format"]
    fasta = dcf.parse_dnacauldron_sequence_text(">x\nACGT\n>y\nGGCC\n")
    assert fasta["format"] == "fasta"
    assert [row["id"] for row in fasta["records"]] == ["x", "y"]
    with pytest.raises(ConfigurationError):
        dcf.parse_dnacauldron_sequence_text("not a sequence")
    with pytest.raises(ConfigurationError):
        dcf.parse_dnacauldron_sequence_text("  ")


def test_parse_dnacauldron_sequence_text_accepts_what_parse_sequence_text_rejects() -> None:
    pytest.importorskip("dnacauldron")
    pytest.importorskip("pydna")
    with pytest.raises(ConfigurationError):
        pyr.parse_sequence_text("ACGTACGT")
    assert dcf.parse_dnacauldron_sequence_text("ACGTACGT")["record_count"] == 1


# ---------------------------------------------------------------------------
# New DnaChisel wrappers
# ---------------------------------------------------------------------------


def test_find_sequence_difference_positions_matches_dnachisel_and_hand_cases() -> None:
    pytest.importorskip("dnachisel")
    from dnachisel.biotools import sequences_differences_array

    result = dchf.find_sequence_difference_positions("ACGTACGTAC", "ACCTACGAAC")
    assert result["positions"] == [2, 7]
    assert result["difference_array"] == [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
    assert result["difference_array"] == [
        int(x) for x in sequences_differences_array("ACGTACGTAC", "ACCTACGAAC")
    ]
    segments = dch.find_sequence_difference_segments("ACGTACGTAC", "ACCTACGAAC")
    assert segments["differing_positions"] == result["difference_count"] == 2
    with pytest.raises(ConfigurationError):
        dchf.find_sequence_difference_positions("ACGT", "ACG")


def test_annotate_sequence_differences_labels_reference_bases() -> None:
    pytest.importorskip("dnachisel")
    result = dchf.annotate_sequence_differences("ACGTACGTAC", "ACCTACGAAC")
    assert result["features"] == [
        {"start": 2, "end": 3, "type": "misc_feature", "label": "#C"},
        {"start": 7, "end": 8, "type": "misc_feature", "label": "#A"},
    ]
    merged = dchf.annotate_sequence_differences("AATTGGCC", "AACCGGCC", prefix="!")
    assert merged["features"] == [{"start": 2, "end": 4, "type": "misc_feature", "label": "!CC"}]
    segments = dch.find_sequence_difference_segments("AATTGGCC", "AACCGGCC")["segments"]
    assert [[row["start"], row["end"]] for row in merged["features"]] == segments
    same = dchf.annotate_sequence_differences("ACGT", "ACGT")
    assert same["feature_count"] == 0
    with pytest.raises(ConfigurationError):
        dchf.annotate_sequence_differences("ACGT", "ACG")
    with pytest.raises(ConfigurationError):
        dchf.annotate_sequence_differences("ACGT", "ACGA", feature_type="bad type")


# ---------------------------------------------------------------------------
# Kept family-variant rows: the claimed parameter reaches the variant
# ---------------------------------------------------------------------------


def _sorted_sequences(products: object) -> list[str]:
    return sorted(str(item.seq).upper() for item in products)  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    ("method", "function"),
    [
        ("in_fusion", "in_fusion_assembly"),
        ("fusion_pcr", "fusion_pcr_assembly"),
        ("in_vivo", "in_vivo_assembly"),
    ],
)
def test_assemble_pydna_fragments_method_reaches_overlap_variants(
    method: str, function: str
) -> None:
    pytest.importorskip("pydna")
    import pydna.assembly2 as assembly
    from pydna.dseqrecord import Dseqrecord

    overlap = _random_dna(21, 30)
    first = _random_dna(22, 40) + overlap
    second = overlap + _random_dna(23, 40)
    result = pyd.assemble_pydna_fragments([first, second], method=method, min_overlap=25)
    upstream = getattr(assembly, function)([Dseqrecord(first), Dseqrecord(second)], limit=25)
    assert result["method"] == method
    assert result["provenance"]["method"] == f"pydna.assembly2.{function}"
    assert [row["sequence"] for row in result["products"]] == sorted(
        (str(p.seq).upper() for p in upstream), key=lambda s: (len(s), s)
    )
    assert result["product_count"] == len(upstream) >= 1


@pytest.mark.parametrize(
    ("method", "function"),
    [
        ("golden_gate", "golden_gate_assembly"),
        ("restriction_ligation", "restriction_ligation_assembly"),
    ],
)
def test_assemble_pydna_fragments_method_reaches_enzyme_variants(
    method: str, function: str
) -> None:
    pytest.importorskip("pydna")
    import pydna.assembly2 as assembly
    from Bio.Restriction import BsmBI
    from pydna.dseqrecord import Dseqrecord

    def part(left: str, right: str, insert: str) -> str:
        return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

    first = part("AAAA", "CCCC", _random_dna(31, 30))
    second = part("CCCC", "GGGG", _random_dna(32, 30))
    result = pyd.assemble_pydna_fragments(
        [first, second], method=method, enzymes=["BsmBI"], max_products=50
    )
    upstream = getattr(assembly, function)([Dseqrecord(first), Dseqrecord(second)], [BsmBI])
    assert result["provenance"]["method"] == f"pydna.assembly2.{function}"
    assert [row["sequence"] for row in result["products"]] == sorted(
        (str(p.seq).upper() for p in upstream), key=lambda s: (len(s), s)
    )
    assert result["product_count"] == len(upstream) >= 1
    with pytest.raises(ConfigurationError):
        pyd.assemble_pydna_fragments([first, second], method=method)


def test_fuse_by_pcr_gives_the_same_molecule_as_the_fusion_pcr_method() -> None:
    pytest.importorskip("pydna")
    from pydna.dseqrecord import Dseqrecord
    from pydna.fusionpcr import fuse_by_pcr

    overlap = _random_dna(41, 24)
    first = _random_dna(42, 40) + overlap
    second = overlap + _random_dna(43, 40)
    legacy = [str(p.seq).upper() for p in fuse_by_pcr([Dseqrecord(first), Dseqrecord(second)], 15)]
    method = pyd.assemble_pydna_fragments([first, second], method="fusion_pcr", min_overlap=15)
    assert len(legacy) == method["product_count"] == 1
    assert legacy[0] in {method["products"][0]["sequence"], _rc(method["products"][0]["sequence"])}
    assert legacy[0] == first + second[len(overlap) :]


def test_dbd_program_is_reached_by_polymerase_dbd() -> None:
    pytest.importorskip("pydna")
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from pydna.tm import dbd_program, pfu_sso7d_program

    amplicon = pcr(Primer(FORWARD), Primer(REVERSE), Dseqrecord(TEMPLATE), limit=14)
    dbd = pyr.suggest_pcr_program(TEMPLATE, FORWARD, REVERSE, polymerase="dbd")
    sso7d = pyr.suggest_pcr_program(TEMPLATE, FORWARD, REVERSE, polymerase="pfu_sso7d")
    assert dbd["program"] == dbd_program(amplicon) == pfu_sso7d_program(amplicon)
    assert sso7d["program"] == dbd["program"]
    assert dbd["provenance"]["method"] == "pydna.tm.dbd_program"


_GENOME = (
    "CCCAAATTTGGGCCCAAATTT"
    + "ACGTACGTAGGCTTAACCGGTTAGCCTAGGATCGATCGGA"
    + "GGGTTTAAACCCGGGTTTAAA"
    + "ACGTACGTAGGCTTAACCGGTTAGCCTAGGATCGATCGGA"
    + "TTTGGGAAACCCTTTGGGAAA"
)
_LOXP = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
_SITES = ("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
_RECOMBINASE_GENOME = "TTTTATGCCCTAAaaTTGGGGGGGGGGGGGAAaaTTTTTTTCCCTCCCC"
_CRISPR_GENOME = "aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac"


def _genome_cases() -> list[tuple[str, str, dict[str, object], str]]:
    """Return (mechanism, operation, keyword arguments, upstream function)."""

    return [
        (
            "homologous_recombination",
            "excision",
            {"genome": _GENOME.upper(), "limit": 40},
            "homologous_recombination_excision",
        ),
        (
            "homologous_recombination",
            "excision_or_inversion",
            {"genome": _GENOME.upper(), "limit": 40},
            "homologous_recombination_excision_or_inversion",
        ),
        (
            "cre_lox",
            "excision",
            {"genome": "CCCCCC" + _LOXP + "GGGGG" + _LOXP + "AAAAA"},
            "cre_lox_excision",
        ),
        (
            "cre_lox",
            "excision_or_inversion",
            {"genome": "CCCCCC" + _LOXP + "GGGGG" + _LOXP + "AAAAA"},
            "cre_lox_excision_or_inversion",
        ),
        (
            "cre_lox",
            "integration",
            {
                "genome": "CCCCCC" + _LOXP + "AAAAA",
                "inserts": [_LOXP + "GGGGG"],
                "inserts_circular": True,
            },
            "cre_lox_integration",
        ),
        (
            "recombinase",
            "excision",
            {"genome": _RECOMBINASE_GENOME.upper(), "recombinase_sites": _SITES},
            "recombinase_excision",
        ),
        (
            "recombinase",
            "excision_or_inversion",
            {"genome": _RECOMBINASE_GENOME.upper(), "recombinase_sites": _SITES},
            "recombinase_excision_or_inversion",
        ),
        (
            "recombinase",
            "integration",
            {
                "genome": "TTTTATGCCCTAAaaTTGGGG",
                "inserts": ["CCCCAAaaTTTTTTTCCCTGGGG"],
                "inserts_circular": True,
                "recombinase_sites": _SITES,
            },
            "recombinase_integration",
        ),
        (
            "crispr",
            "integration",
            {
                "genome": _CRISPR_GENOME,
                "inserts": ["aaccggttAAAAAAAAAttcaaagcac"],
                "guides": ["ttcaatgcaaacagtaatga"],
                "limit": 8,
            },
            "crispr_integration",
        ),
    ]


def _upstream_genome_call(function: str, arguments: dict[str, object]) -> list[str]:
    import pydna.assembly2 as assembly
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from pydna.recombinase import Recombinase

    genome = Dseqrecord(str(arguments["genome"]).upper())
    args: list[object] = [genome]
    if "inserts" in arguments:
        inserts = [
            Dseqrecord(str(item).upper(), circular=bool(arguments.get("inserts_circular")))
            for item in arguments["inserts"]  # type: ignore[attr-defined]
        ]
        args.append(inserts)
    if "recombinase_sites" in arguments:
        args.append(Recombinase(*arguments["recombinase_sites"]))  # type: ignore[misc]
    if "guides" in arguments:
        args.append([Primer(str(item).upper()) for item in arguments["guides"]])  # type: ignore[attr-defined]
    if function in {
        "homologous_recombination_excision",
        "homologous_recombination_excision_or_inversion",
        "crispr_integration",
    }:
        args.append(int(arguments["limit"]))  # type: ignore[call-overload]
    return _sorted_sequences(getattr(assembly, function)(*args))


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
@pytest.mark.parametrize(
    ("mechanism", "operation", "arguments", "function"),
    _genome_cases(),
    ids=lambda value: value if isinstance(value, str) else "",
)
def test_simulate_genome_edit_reaches_each_pydna_edit_function(
    mechanism: str, operation: str, arguments: dict[str, object], function: str
) -> None:
    pytest.importorskip("pydna")
    result = pyr.simulate_genome_edit(mechanism=mechanism, operation=operation, **arguments)  # type: ignore[arg-type]
    assert result["provenance"]["method"] == f"pydna.assembly2.{function}"
    expected = _upstream_genome_call(function, arguments)
    got = sorted(row["sequence"] for row in result["products"])
    assert got == sorted(item.upper() for item in expected)
    assert result["product_count"] == len(expected) >= 1


# ---------------------------------------------------------------------------
# duplicate-of rows: same object or same code
# ---------------------------------------------------------------------------

_SAME_OBJECT = [
    ("pydna.all.pcr", "pydna.amplify.pcr"),
    ("pydna.design.pcr", "pydna.amplify.pcr"),
    ("pydna.all.genbank", "pydna.genbank.genbank"),
    ("pydna.all.read", "pydna.readers.read"),
    ("pydna.genbank.read", "pydna.readers.read"),
    ("pydna.all.read_primer", "pydna.readers.read_primer"),
    ("pydna.all.parse_primers", "pydna.parsers.parse_primers"),
    ("pydna.all.eq", "pydna.utils.eq"),
    ("pydna.amplicon.dbd_program", "pydna.tm.dbd_program"),
    ("pydna.tm.pfu_sso7d_program", "pydna.tm.dbd_program"),
    ("pydna.tm.program", "pydna.amplicon.program"),
    ("pydna.tm.taq_program", "pydna.amplicon.program"),
    ("pydna.assembly.common_sub_strings", "pydna.common_sub_strings.common_sub_strings"),
    ("pydna.assembly2.common_sub_strings_str", "pydna.common_sub_strings.common_sub_strings"),
    ("pydna.dseq.common_sub_strings", "pydna.common_sub_strings.common_sub_strings"),
    ("pydna.dseqrecord.common_sub_strings", "pydna.common_sub_strings.common_sub_strings"),
    ("pydna.seqrecord.common_sub_strings", "pydna.common_sub_strings.common_sub_strings"),
    (
        "pydna.oligonucleotide_hybridization.common_sub_strings",
        "pydna.common_sub_strings.common_sub_strings",
    ),
    ("pydna.fusionpcr.terminal_overlap", "pydna.common_sub_strings.terminal_overlap"),
    ("pydna.design.tm_default", "pydna.tm.tm_default"),
    ("pydna.snapgene_history_parser.gibson_assembly", "pydna.assembly2.gibson_assembly"),
    ("pydna.snapgene_history_parser.in_fusion_assembly", "pydna.assembly2.in_fusion_assembly"),
    ("pydna.snapgene_history_parser.fusion_pcr_assembly", "pydna.assembly2.fusion_pcr_assembly"),
    (
        "pydna.snapgene_history_parser.restriction_ligation_assembly",
        "pydna.assembly2.restriction_ligation_assembly",
    ),
    ("pydna.snapgene_history_parser.gateway_assembly", "pydna.assembly2.gateway_assembly"),
    ("pydna.snapgene_history_parser.ligation_assembly", "pydna.assembly2.ligation_assembly"),
    ("pydna.snapgene_history_parser.pcr_assembly", "pydna.assembly2.pcr_assembly"),
    (
        "pydna.snapgene_history_parser.oligonucleotide_hybridization",
        "pydna.oligonucleotide_hybridization.oligonucleotide_hybridization",
    ),
    (
        "dnacauldron.AssemblyMix.generate_type2s_restriction_mix",
        "dnacauldron.generate_type2s_restriction_mix",
    ),
    (
        "dnacauldron.utils.utils.generate_type2s_restriction_mix",
        "dnacauldron.generate_type2s_restriction_mix",
    ),
    ("dnacauldron.biotools.autoselect_enzyme.autoselect_enzyme", "dnacauldron.autoselect_enzyme"),
    ("dnacauldron.utils.utils.autoselect_enzyme", "dnacauldron.autoselect_enzyme"),
    ("dnacauldron.utils.swap_donor_vector_part", "dnacauldron.swap_donor_vector_part"),
    (
        "dnacauldron.utils.insert_parts_on_backbones.insert_parts_on_backbones",
        "dnacauldron.insert_parts_on_backbones",
    ),
    (
        "dnacauldron.utils.insert_parts_on_backbones.record_contains_backbone",
        "dnacauldron.utils.record_contains_backbone",
    ),
    ("dnacauldron.utils.substitute_overhangs", "dnacauldron.utils.utils.substitute_overhangs"),
    (
        "dnacauldron.utils.utils.list_digestion_overhangs",
        "dnacauldron.utils.list_digestion_overhangs",
    ),
    (
        "dnacauldron.biotools.record_operations.crop_record_with_saddling_features",
        "dnacauldron.biotools.crop_record_with_saddling_features",
    ),
    ("dnachisel.biotools.list_common_enzymes", "dnachisel.list_common_enzymes"),
    ("dnachisel.biotools.enzymes_operations.list_common_enzymes", "dnachisel.list_common_enzymes"),
    (
        "dnachisel.biotools.sequences_differences_segments",
        "dnachisel.sequences_differences_segments",
    ),
    (
        "dnachisel.biotools.sequences_differences.sequences_differences_segments",
        "dnachisel.sequences_differences_segments",
    ),
    (
        "dnachisel.biotools.biotables.get_backtranslation_table",
        "dnachisel.biotools.get_backtranslation_table",
    ),
    (
        "dnachisel.biotools.indices_operations.group_nearby_indices",
        "dnachisel.biotools.group_nearby_indices",
    ),
    ("dnachisel.utils.random_compatible_dna_sequence", "dnachisel.random_compatible_dna_sequence"),
    (
        "dnachisel.utils.utils.random_compatible_dna_sequence",
        "dnachisel.random_compatible_dna_sequence",
    ),
    (
        "dnachisel.annotate_differences",
        "dnachisel.biotools.genbank_operations.annotate_differences",
    ),
    (
        "dnachisel.biotools.annotate_differences",
        "dnachisel.biotools.genbank_operations.annotate_differences",
    ),
    (
        "dnachisel.biotools.sequences_differences_array",
        "dnachisel.biotools.sequences_differences.sequences_differences_array",
    ),
    ("primer3.bindings.calc_hairpin_tm", "primer3.calc_hairpin_tm"),
    ("primer3.bindings.calc_homodimer_tm", "primer3.calc_homodimer_tm"),
    ("primer3.bindings.calc_heterodimer_tm", "primer3.calc_heterodimer_tm"),
]


@pytest.mark.parametrize(("alias", "target"), _SAME_OBJECT)
def test_duplicate_rows_are_the_same_object(alias: str, target: str) -> None:
    root = alias.split(".")[0]
    pytest.importorskip(root)
    assert _resolve(alias) is _resolve(target)


def test_terminal_overlap_of_records_equals_the_string_version() -> None:
    pytest.importorskip("pydna")
    from pydna.assembly2 import terminal_overlap
    from pydna.common_sub_strings import terminal_overlap as string_overlap
    from pydna.dseqrecord import Dseqrecord

    rng = random.Random(9)
    for _ in range(60):
        x = "".join(rng.choice("ACGT") for _ in range(rng.randint(10, 30)))
        y = x[rng.randint(0, 5) :] + "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 10)))
        assert sorted(terminal_overlap(Dseqrecord(x), Dseqrecord(y), 6)) == sorted(
            string_overlap(x, y, 6)
        )
    result = pyd.find_terminal_overlaps("TTACTAAAAAAA", "AAAAAAACGCACG", min_length=5)
    assert [(r["x_start"], r["y_start"], r["length"]) for r in result["matches"]] == sorted(
        terminal_overlap(Dseqrecord("TTACTAAAAAAA"), Dseqrecord("AAAAAAACGCACG"), 5)
    )


def test_read_and_read_primer_are_parse_with_a_single_record_check() -> None:
    pytest.importorskip("pydna")
    from pydna.readers import read, read_primer

    text = ">one\nACGTACGTAGCTAG\n"
    record = read(text)
    parsed = pyr.parse_sequence_text(text)
    assert parsed["record_count"] == 1
    assert parsed["records"][0]["sequence"] == str(record.seq).upper()
    assert parsed["records"][0]["name"] == record.name
    with pytest.raises(ValueError, match="More than one"):
        read(">a\nACGT\n>b\nGGCC\n")
    primer = read_primer(">p1\nACGTACGTAGCTAG\n")
    primers = pyd.parse_primer_sequences(">p1\nACGTACGTAGCTAG\n")
    assert primers["primer_count"] == 1
    assert primers["primers"][0]["sequence"] == str(primer.seq).upper()
    assert primers["primers"][0]["name"] == str(primer.name)


def test_hybridization_overhangs_are_returned_by_the_wrapper() -> None:
    pytest.importorskip("pydna")
    from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization_overhangs

    result = pyd.primer_pair_hybridization("ATGGC", "GCCAT", minimal_annealing=4)
    assert (
        result["overhangs"] == oligonucleotide_hybridization_overhangs("ATGGC", "GCCAT", 4) == [0]
    )


def test_annotate_sites_equal_the_find_recombination_sites_result() -> None:
    pytest.importorskip("pydna")
    from Bio.Data.IUPACData import ambiguous_dna_values
    from pydna.cre_lox import LOXP_SEQUENCE, annotate_loxP_sites
    from pydna.dseqrecord import Dseqrecord
    from pydna.gateway import annotate_gateway_sites, recombinase_dict

    from dnakit.molecular_lab import find_recombination_sites

    lox = "CCCCCC" + LOXP_SEQUENCE + "AAAAA" + LOXP_SEQUENCE + "GG"
    found = find_recombination_sites(lox, system="loxp", verbose=False)
    annotated = annotate_loxP_sites(Dseqrecord(lox))
    assert [
        (f.qualifiers["label"][0], int(f.location.start), int(f.location.end), f.location.strand)
        for f in annotated.features
    ] == [(s["name"], s["start"], s["end"], s["strand"]) for s in found["sites"]]
    assert found["site_count"] == 2

    site = recombinase_dict["BP"]["conservative"].recombinases[0]
    att_b = "".join(
        ambiguous_dna_values[c][0] if c in ambiguous_dna_values else c for c in site.site1
    )
    gateway = "TTTT" + att_b + "CCCC"
    gw_found = find_recombination_sites(gateway, system="gateway", verbose=False)
    gw_annotated = annotate_gateway_sites(Dseqrecord(gateway), False)
    assert gw_found["site_count"] == len(gw_annotated.features) >= 1
    assert sorted((int(f.location.start), int(f.location.end)) for f in gw_annotated.features) == [
        (s["start"], s["end"]) for s in gw_found["sites"]
    ]


def test_primer3_tm_only_functions_equal_the_wrapper_tm_field() -> None:
    primer3 = pytest.importorskip("primer3")
    from dnakit.molecular_lab import calculate_primer3_dimer, calculate_primer3_hairpin

    hairpin = "GCGCGCAAAAGCGCGC"
    partner = "ACGTACGTAGCTAGCTAGGC"
    assert calculate_primer3_hairpin(hairpin, verbose=False)["tm_c"] == primer3.calc_hairpin_tm(
        hairpin
    )
    assert calculate_primer3_dimer(hairpin, verbose=False)["tm_c"] == primer3.calc_homodimer_tm(
        hairpin
    )
    assert calculate_primer3_dimer(hairpin, partner=partner, verbose=False)[
        "tm_c"
    ] == primer3.calc_heterodimer_tm(hairpin, partner)
    assert primer3.calc_hairpin_tm(hairpin) == primer3.calc_hairpin(hairpin).tm


def test_primer3_bytes_helpers_equal_the_string_wrappers() -> None:
    pytest.importorskip("primer3")
    from primer3 import p3helpers

    from dnakit.molecular_lab import sanitize_primer3_sequence, validate_acgt_sequence

    checked = validate_acgt_sequence("acgtacgt", verbose=False)["sequence"]
    assert p3helpers.ensure_acgt_uppercase_b(b"acgtacgt").decode() == checked == "ACGTACGT"
    assert p3helpers.ensure_acgt_uppercase("acgtacgt") == checked
    cleaned = sanitize_primer3_sequence("acgtRYacgn", verbose=False)
    assert p3helpers.sanitize_sequence_b(b"acgtRYacgn").decode() == cleaned["sanitized"]
    with pytest.raises(ValueError):
        p3helpers.ensure_acgt_uppercase_b(b"acgn")
    with pytest.raises(ConfigurationError):
        validate_acgt_sequence("acgn", verbose=False)


@pytest.mark.parametrize(
    "function",
    ["homologous_recombination_excision", "cre_lox_excision", "recombinase_excision"],
)
def test_plain_excision_functions_are_deprecated_upstream(function: str) -> None:
    pytest.importorskip("pydna")

    cases = {item[3]: item[2] for item in _genome_cases()}
    with pytest.warns(DeprecationWarning, match="deprecated"):
        _upstream_genome_call(function, cases[function])


def test_tm_q5_is_not_implemented_and_tm_neb_needs_the_neb_service() -> None:
    pytest.importorskip("pydna")
    import inspect

    from pydna.tm import Q5, tm_neb

    with pytest.raises(NotImplementedError):
        Q5(FORWARD)
    assert "https://tmapi.neb.com/tm" in inspect.getsource(tm_neb)
