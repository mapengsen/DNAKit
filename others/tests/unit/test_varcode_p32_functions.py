"""Tests for dnakit.molecular_lab.varcode_p32_functions (varcode collections, effects, genotypes, phasing, germline).

A small synthetic genome (one contig of 900 bp, seeded random sequence) carries two transcripts: T1, a
single-exon ORF of 60 codons, and T3, a two-exon ORF of 40 codons with an intron of 140 bases. A multi-sample
VCF holds six SNVs (four coding, one silent, one intergenic, one at the T3 donor), and an RNA BAM is
written with pysam. Expected effects are derived here from the codon table, not copied from varcode output;
the varcode results are compared with direct varcode calls on the same reference. Nothing is downloaded.
"""

from __future__ import annotations

import inspect
import itertools
import random
import sys
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("varcode")
pytest.importorskip("pyensembl")
pysam = pytest.importorskip("pysam")

import varcode

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import varcode_p32_functions as vp
from dnakit.molecular_lab.varcode_gap_functions import _genome

CODON = {
    "".join(bases): amino
    for bases, amino in zip(
        itertools.product("TCAG", repeat=3),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
}
SENSE = sorted(codon for codon, amino in CODON.items() if amino != "*" and codon != "ATG")
PUBLIC = list(vp.__all__)


def _coding(rng: random.Random, codons: int, stop: bool) -> str:
    body = "ATG" + "".join(rng.choice(SENSE) for _ in range(codons - 1 - int(stop)))
    return body + ("TAA" if stop else "")


def _flank(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _gtf_rows(feature_lines: list[tuple[str, str, int, int, str]]) -> str:
    out = []
    for _contig_strand, feature, start, end, attributes in feature_lines:
        frame = "0" if feature == "CDS" else "."
        out.append(
            "\t".join(["chr1", "test", feature, str(start), str(end), ".", "+", frame, attributes])
        )
    return "\n".join(out) + "\n"


def _transcript_lines(
    gene: str, transcript: str, exons: list[tuple[int, int]]
) -> list[tuple[str, str, int, int, str]]:
    tx = f'gene_id "{gene}"; transcript_id "{transcript}"; transcript_name "{transcript}-001"; transcript_biotype "protein_coding";'
    span = (exons[0][0], exons[-1][1])
    lines = [
        (
            "+",
            "gene",
            *span,
            f'gene_id "{gene}"; gene_name "{gene}"; gene_biotype "protein_coding";',
        ),
        ("+", "transcript", *span, tx),
    ]
    for number, (start, end) in enumerate(exons, start=1):
        lines.append(("+", "exon", start, end, f'{tx} exon_number "{number}";'))
    for number, (start, end) in enumerate(exons, start=1):
        lines.append(
            ("+", "CDS", start, end, f'{tx} exon_number "{number}"; protein_id "P_{transcript}";')
        )
    lines.append(("+", "start_codon", exons[0][0], exons[0][0] + 2, f'{tx} exon_number "1";'))
    lines.append(
        ("+", "stop_codon", exons[-1][1] - 2, exons[-1][1], f'{tx} exon_number "{len(exons)}";')
    )
    return lines


def _alt_for(ref: str) -> str:
    return {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    rng = random.Random(20261020)
    t1 = _coding(rng, 60, True)  # 180 bases, 101..280
    t3a = _coding(rng, 20, False)  # 60 bases, 401..460
    t3b = "".join(rng.choice(SENSE) for _ in range(19)) + "TAA"  # 60 bases, 601..660
    intron = "GT" + _flank(rng, 136) + "AG"  # 461..600
    genome = _flank(rng, 100) + t1 + _flank(rng, 120) + t3a + intron + t3b + _flank(rng, 240)
    assert len(genome) == 900 and genome[100:280] == t1 and genome[400:460] == t3a
    directory = tmp_path_factory.mktemp("varcode_p32")
    fasta = directory / "genome.fa"
    fasta.write_text(
        ">chr1\n" + "\n".join(genome[i : i + 60] for i in range(0, 900, 60)) + "\n",
        encoding="utf-8",
    )
    lines = _transcript_lines("G1", "T1", [(101, 280)]) + _transcript_lines(
        "G3", "T3", [(401, 460), (601, 660)]
    )
    gtf = directory / "genes.gtf"
    gtf.write_text(_gtf_rows(lines), encoding="utf-8")

    def base(pos: int) -> str:
        return genome[pos - 1]

    def codon_at(index: int) -> tuple[int, str]:
        start = 101 + 3 * (index - 1)
        return start, genome[start - 1 : start + 2]

    # coding variants chosen from the codon table
    positions: dict[str, int] = {}
    for label, index in (("missense_a", 3), ("missense_b", 4), ("missense_d", 8)):
        start, codon = codon_at(index)
        for offset in range(3):
            changed = codon[:offset] + _alt_for(codon[offset]) + codon[offset + 1 :]
            if CODON[changed] not in (CODON[codon], "*"):
                positions[label] = start + offset
                break
        else:
            raise AssertionError("no missense position found")
    for index in range(5, 60):
        start, codon = codon_at(index)
        changed = codon[:2] + _alt_for(codon[2])
        if CODON[changed] == CODON[codon]:
            positions["silent_c"] = start + 2
            break
    else:
        raise AssertionError("no silent position found")
    positions["intergenic_e"] = 50
    positions["donor_f"] = 462
    variants = sorted(positions.items(), key=lambda item: item[1])
    genotypes = {
        "missense_a": ("0|1:1:5,7", "0/1:.:4,4"),
        "missense_b": ("0|1:1:6,6", "./.:.:."),
        "silent_c": ("1|1:.:0,9", "0/0:.:8,0"),
        "missense_d": ("0/1:.:5,5", "1/1:.:0,7"),
        "intergenic_e": ("0/1:.:5,5", "0/0:.:9,0"),
        "donor_f": ("0/1:.:5,5", "0/1:.:5,5"),
    }
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=900>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="p">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">\n'
    )
    body = []
    for label, pos in variants:
        ref = base(pos)
        s1, s2 = genotypes[label]
        body.append(f"chr1\t{pos}\t.\t{ref}\t{_alt_for(ref)}\t50\tPASS\t.\tGT:PS:AD\t{s1}\t{s2}")
    column = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    multi = directory / "multi.vcf"
    multi.write_text(header + column + "\tS1\tS2\n" + "\n".join(body) + "\n", encoding="utf-8")
    single = directory / "single.vcf"
    single.write_text(
        header
        + column
        + "\tS1\n"
        + "\n".join("\t".join(line.split("\t")[:10]) for line in body)
        + "\n",
        encoding="utf-8",
    )
    # structural variants: a breakend pair with inserted bases and a deletion
    sv = directory / "sv.vcf"
    sv.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=900>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t300\tb1\t{base(300)}\t{base(300)}ACGT[chr1:700[\t.\tPASS\tSVTYPE=BND\n"
        f"chr1\t800\td1\t{base(800)}\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=850\n",
        encoding="utf-8",
    )
    # RNA BAM: 6 fragments carry both missense_a and missense_b alt alleles, 4 are reference
    pos_a, pos_b = positions["missense_a"], positions["missense_b"]
    reads = []
    for index in range(10):
        seq = list(genome[95:145])
        if index < 6:
            seq[pos_a - 96] = _alt_for(seq[pos_a - 96])
            seq[pos_b - 96] = _alt_for(seq[pos_b - 96])
        read = pysam.AlignedSegment()
        read.query_name = f"r{index}"
        read.query_sequence = "".join(seq)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 95
        read.mapping_quality = 60
        read.cigartuples = [(0, 50)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        reads.append(read)
    bam = directory / "rna.bam"
    header_dict = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 900}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header_dict) as handle:
        for read in reads:
            handle.write(read)
    pysam.index(str(bam))
    return {
        "dir": directory,
        "fasta": str(fasta),
        "gtf": str(gtf),
        "cache": str(directory / "cache"),
        "multi": str(multi),
        "single": str(single),
        "sv": str(sv),
        "bam": str(bam),
        "genome": genome,
        "positions": positions,
        "base": base,
    }


def _call(function: Any, world: dict[str, Any], *args: Any, **kwargs: Any) -> Any:
    kwargs.setdefault("verbose", False)
    kwargs.setdefault("cache_dir", world["cache"])
    return function(world["gtf"], world["fasta"], *args, **kwargs)


def _variant(world: dict[str, Any], label: str) -> dict[str, Any]:
    pos = world["positions"][label]
    ref = world["base"](pos)
    return {"contig": "chr1", "start": pos, "ref": ref, "alt": _alt_for(ref)}


def _direct(world: dict[str, Any], path: str | None = None) -> tuple[Any, Any]:
    genome, _ = _genome(world["gtf"], world["fasta"], world["cache"])
    collection = varcode.load_vcf(
        path or world["multi"], genome=genome, convert_ucsc_contig_names=False
    )
    return genome, collection


def _keys(rows: list[dict[str, Any]]) -> list[tuple[Any, ...]]:
    return [(row["contig"], row["start"], row["ref"], row["alt"]) for row in rows]


def test_public_functions_are_keyword_verbose() -> None:
    assert len(PUBLIC) == 22
    for name in PUBLIC:
        parameter = inspect.signature(getattr(vp, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True


def test_missing_varcode_raises_backend_unavailable(
    world: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "varcode", None)
    with pytest.raises(BackendUnavailableError, match="bio-varcode"):
        vp.varcode_genotype_alt_summary("0/1", verbose=False)


# ---- collections -----------------------------------------------------------------------------


def test_variant_table_lists_genes_and_writes_csv(world: dict[str, Any]) -> None:
    result = _call(vp.varcode_variant_table, world, vcf_paths=[world["multi"]])
    assert result["count"] == 6
    by_start = {row["start"]: row for row in result["variants"]}
    positions = world["positions"]
    assert by_start[positions["missense_a"]]["gene_name"] == "G1"
    assert by_start[positions["donor_f"]]["gene_id"] == "G3"
    assert by_start[positions["intergenic_e"]]["gene_name"] == ""
    _, direct = _direct(world)
    assert result["variants"] == vp._table(direct.to_dataframe())
    target = world["dir"] / "variants_out.csv"
    again = _call(
        vp.varcode_variant_table, world, vcf_paths=[world["multi"]], output_csv=str(target)
    )
    assert again["output_csv"] == str(target)
    lines = [line for line in target.read_text().splitlines() if not line.startswith("#")]
    assert lines[0] == "chr,start,ref,alt,gene_name,gene_id" and len(lines) == 7
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_variant_table, world, vcf_paths=[world["multi"]], output_csv=str(target))


@pytest.mark.parametrize("ref,alt", [("C", "CAA"), ("ACT", "A"), ("", "AA")])
def test_variant_table_preserves_the_supplied_allele_representation(world, ref, alt):
    record = {"contig": "chr1", "start": 20, "ref": ref, "alt": alt}
    genome, _ = _genome(world["gtf"], world["fasta"], world["cache"])
    direct = varcode.VariantCollection([varcode.Variant("chr1", 20, ref, alt, genome=genome)])
    result = _call(vp.varcode_variant_table, world, variants=[record])
    assert result["variants"] == vp._table(direct.to_dataframe())


def test_read_variant_csv_round_trips(world: dict[str, Any]) -> None:
    target = world["dir"] / "roundtrip.csv"
    _call(vp.varcode_variant_table, world, vcf_paths=[world["multi"]], output_csv=str(target))
    result = _call(vp.varcode_read_variant_csv, world, str(target))
    original = _call(vp.varcode_variant_table, world, vcf_paths=[world["multi"]])
    assert result["count"] == 6
    assert _keys(result["variants"]) == [
        ("chr1", row["start"], row["ref"], row["alt"]) for row in original["variants"]
    ]
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_read_variant_csv, world, str(world["dir"] / "missing.csv"))


def test_variant_set_operations(world: dict[str, Any]) -> None:
    a = _variant(world, "missense_a")
    e = _variant(world, "intergenic_e")
    extra = {
        "contig": "chr1",
        "start": 700,
        "ref": world["base"](700),
        "alt": _alt_for(world["base"](700)),
    }
    first = {"vcf": [world["multi"]]}
    second = {"variants": [a, e, extra]}
    union = _call(vp.varcode_variant_set_operation, world, "union", first, second)
    inter = _call(vp.varcode_variant_set_operation, world, "intersection", first, second)
    diff = _call(vp.varcode_variant_set_operation, world, "difference", first, second)
    assert union["count"] == 7 and (inter["count"], diff["count"]) == (2, 4)
    assert {row["start"] for row in inter["variants"]} == {a["start"], e["start"]}
    assert {row["start"] for row in diff["variants"]} == set(world["positions"].values()) - {
        a["start"],
        e["start"],
    }
    swapped = _call(vp.varcode_variant_set_operation, world, "difference", second, first)
    assert [row["start"] for row in swapped["variants"]] == [700]
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_variant_set_operation, world, "xor", first, second)
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_variant_set_operation, world, "union", first, {"bogus": []})


def test_filter_by_expression_variants_and_effects(world: dict[str, Any]) -> None:
    positions = world["positions"]
    kept = _call(
        vp.varcode_filter_by_expression,
        world,
        {"T1": 5.0, "T3": 0.0},
        level="transcript",
        vcf_paths=[world["multi"]],
    )
    assert {row["start"] for row in kept["items"]} == {
        positions[k] for k in ("missense_a", "missense_b", "silent_c", "missense_d")
    }
    _, direct = _direct(world)
    assert len(direct.filter_by_transcript_expression({"T1": 5.0, "T3": 0.0})) == kept["count"]
    gene = _call(
        vp.varcode_filter_by_expression,
        world,
        {"G3": 2.5},
        level="gene",
        target="effects",
        vcf_paths=[world["multi"]],
    )
    assert gene["count"] == 1 and gene["items"][0]["transcript_id"] == "T3"
    none = _call(
        vp.varcode_filter_by_expression,
        world,
        {"T1": 5.0},
        level="transcript",
        min_expression_value=5.0,
        vcf_paths=[world["multi"]],
    )
    assert none["count"] == 0
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_filter_by_expression, world, {}, vcf_paths=[world["multi"]])
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_filter_by_expression, world, {"T1": -1}, vcf_paths=[world["multi"]])


# ---- genotypes -------------------------------------------------------------------------------


def test_variants_in_sample_by_zygosity(world: dict[str, Any]) -> None:
    positions = world["positions"]
    carried = _call(vp.varcode_variants_in_sample, world, [world["multi"]], "S1")
    het = _call(
        vp.varcode_variants_in_sample, world, [world["multi"]], "S1", genotype="heterozygous"
    )
    hom = _call(
        vp.varcode_variants_in_sample, world, [world["multi"]], "S1", genotype="homozygous_alt"
    )
    assert {row["start"] for row in hom["variants"]} == {positions["silent_c"]}
    assert {row["start"] for row in het["variants"]} == {
        positions[k] for k in ("missense_a", "missense_b", "missense_d", "intergenic_e", "donor_f")
    }
    _, direct = _direct(world)
    assert carried["count"] == len(direct.for_sample("S1"))
    s2_hom = _call(
        vp.varcode_variants_in_sample, world, [world["multi"]], "S2", genotype="homozygous_alt"
    )
    assert {row["start"] for row in s2_hom["variants"]} == {positions["missense_d"]}
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_variants_in_sample, world, [world["multi"]], "S9")
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_variants_in_sample, world, [world["multi"]], "S1", genotype="all")


def test_sample_genotype_matches_vcf_fields(world: dict[str, Any]) -> None:
    result = _call(
        vp.varcode_sample_genotype, world, [world["multi"]], "S1", _variant(world, "missense_a")
    )
    gt = result["genotype"]
    assert (gt["raw_gt"], gt["phased"], gt["phase_set"], gt["allele_depths"]) == (
        "0|1",
        True,
        1,
        [5, 7],
    )
    assert (gt["zygosity"], gt["carries_alt"], gt["copies_of_alt"], gt["depth_for_alt"]) == (
        "het",
        True,
        1,
        7,
    )
    assert gt["collection_zygosity"] == "het"
    hom = _call(
        vp.varcode_sample_genotype, world, [world["multi"]], "S1", _variant(world, "silent_c")
    )["genotype"]
    assert (hom["zygosity"], hom["copies_of_alt"], hom["depth_for_alt"]) == ("hom", 2, 9)
    ref = _call(
        vp.varcode_sample_genotype, world, [world["multi"]], "S2", _variant(world, "silent_c")
    )["genotype"]
    assert (ref["carries_alt"], ref["copies_of_alt"]) == (False, 0)
    missing = _call(
        vp.varcode_sample_genotype, world, [world["multi"]], "S2", _variant(world, "missense_b")
    )["genotype"]
    assert missing["alleles"] == [None, None] and missing["carries_alt"] is False
    with pytest.raises(ConfigurationError):
        _call(
            vp.varcode_sample_genotype,
            world,
            [world["multi"]],
            "S1",
            {"contig": "chr1", "start": 5, "ref": "A", "alt": "C"},
        )


def test_genotype_alt_summary_without_reference() -> None:
    het = vp.varcode_genotype_alt_summary("0/1", allele_depths=[3, 8], verbose=False)["genotype"]
    assert (
        het["carries_alt"],
        het["copies_of_alt"],
        het["depth_for_alt"],
        het["zygosity"],
        het["phased"],
    ) == (True, 1, 8, "het", False)
    hom = vp.varcode_genotype_alt_summary("1|1", allele_depths=[0, 12], phase_set=4, verbose=False)[
        "genotype"
    ]
    assert (hom["copies_of_alt"], hom["zygosity"], hom["phase_set"], hom["ploidy"]) == (
        2,
        "hom",
        4,
        2,
    )
    second = vp.varcode_genotype_alt_summary(
        "1/2", allele_depths=[0, 4, 6], alt_index=2, verbose=False
    )["genotype"]
    assert (second["copies_of_alt"], second["depth_for_alt"]) == (1, 6)
    none = vp.varcode_genotype_alt_summary("./.", verbose=False)["genotype"]
    assert (none["is_called"], none["carries_alt"], none["depth_for_alt"]) == (False, False, None)
    haploid = vp.varcode_genotype_alt_summary("1", verbose=False)["genotype"]
    assert (haploid["ploidy"], haploid["copies_of_alt"]) == (1, 1)
    with pytest.raises(ConfigurationError):
        vp.varcode_genotype_alt_summary("", verbose=False)
    with pytest.raises(ConfigurationError):
        vp.varcode_genotype_alt_summary("0/1", alt_index=0, verbose=False)


# ---- effects ---------------------------------------------------------------------------------


def _expected_classes(world: dict[str, Any]) -> dict[int, str]:
    positions = world["positions"]
    return {
        positions["missense_a"]: "Substitution",
        positions["missense_b"]: "Substitution",
        positions["missense_d"]: "Substitution",
        positions["silent_c"]: "Silent",
        positions["intergenic_e"]: "Intergenic",
        positions["donor_f"]: "SpliceOutcomeSet",
    }


def test_effect_table_classes_and_csv(world: dict[str, Any]) -> None:
    result = _call(vp.varcode_effect_table, world, vcf_paths=[world["multi"]])
    assert {row["start"]: row["effect_type"] for row in result["effects"]} == _expected_classes(
        world
    )
    _, direct = _direct(world)
    assert result["effects"] == vp._table(direct.effects().to_dataframe())
    target = world["dir"] / "effects_out.csv"
    _call(vp.varcode_effect_table, world, vcf_paths=[world["multi"]], output_csv=str(target))
    assert any(line.startswith("variant,contig,start") for line in target.read_text().splitlines())
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_effect_table, world, vcf_paths=[world["multi"]], output_csv=str(target))


def test_read_effect_csv_reannotates(world: dict[str, Any]) -> None:
    target = world["dir"] / "effects_in.csv"
    _call(vp.varcode_effect_table, world, vcf_paths=[world["multi"]], output_csv=str(target))
    result = _call(vp.varcode_read_effect_csv, world, str(target))
    assert result["count"] == 6
    assert sorted(row["class"] for row in result["effects"]) == sorted(
        _expected_classes(world).values()
    )


def test_effect_priority_ranks_are_ordered() -> None:
    result = vp.varcode_effect_priority(
        ["Intergenic", "Silent", "Substitution", "PrematureStop"], verbose=False
    )
    rank = {row["effect_class"]: row["priority"] for row in result["ranks"]}
    assert rank["PrematureStop"] > rank["Substitution"] > rank["Silent"] > rank["Intergenic"]
    everything = vp.varcode_effect_priority(verbose=False)
    priorities = [row["priority"] for row in everything["ranks"]]
    assert priorities == sorted(priorities, reverse=True) and everything["count"] > 20
    from varcode.effects import effect_ordering

    assert (
        rank["Substitution"]
        == effect_ordering.transcript_effect_priority_dict[varcode.effects.Substitution]
    )
    with pytest.raises(ConfigurationError):
        vp.varcode_effect_priority(["NotAnEffect"], verbose=False)


def test_top_priority_effects_by_scope(world: dict[str, Any]) -> None:
    kwargs = {"vcf_paths": [world["multi"]]}
    overall = _call(vp.varcode_top_priority_effects, world, scope="overall", **kwargs)
    per_variant = _call(vp.varcode_top_priority_effects, world, scope="variant", **kwargs)
    assert overall["count"] == 1
    assert overall["effects"][0]["priority"] == max(
        row["priority"] for row in per_variant["effects"]
    )
    assert per_variant["count"] == 6
    per_transcript = _call(vp.varcode_top_priority_effects, world, scope="transcript", **kwargs)
    assert {row["key"] for row in per_transcript["effects"]} == {"T1", "T3", None}
    per_gene = _call(vp.varcode_top_priority_effects, world, scope="gene", **kwargs)
    t1 = next(row for row in per_gene["effects"] if row["key"] == "G1")
    _, direct = _direct(world)
    assert (
        t1["description"]
        == direct.effects().top_priority_effect_per_gene_id()["G1"].short_description
    )
    single = _call(
        vp.varcode_top_priority_effects,
        world,
        scope="overall",
        variants=[_variant(world, "silent_c")],
    )
    assert single["effects"][0]["class"] == "Silent"
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_top_priority_effects, world, scope="all", **kwargs)


def test_filter_effects_drops_silent_and_applies_priority(world: dict[str, Any]) -> None:
    kwargs = {"vcf_paths": [world["multi"]]}
    coding = _call(vp.varcode_filter_effects, world, **kwargs)
    assert coding["input_count"] == 6
    assert sorted(row["class"] for row in coding["effects"]) == [
        "SpliceOutcomeSet",
        "Substitution",
        "Substitution",
        "Substitution",
    ]
    _, direct = _direct(world)
    assert coding["count"] == len(direct.effects().drop_silent_and_noncoding())
    nothing_dropped = _call(
        vp.varcode_filter_effects, world, drop_silent_and_noncoding=False, **kwargs
    )
    assert nothing_dropped["count"] == 6
    severe = _call(
        vp.varcode_filter_effects,
        world,
        drop_silent_and_noncoding=False,
        min_priority_class="Substitution",
        **kwargs,
    )
    assert [row["class"] for row in severe["effects"]] == [
        "Substitution"
    ] * 3  # SpliceDonor (19) ranks below Substitution (21)
    down_to_donor = _call(
        vp.varcode_filter_effects,
        world,
        drop_silent_and_noncoding=False,
        min_priority_class="SpliceDonor",
        **kwargs,
    )
    assert sorted(row["class"] for row in down_to_donor["effects"]) == [
        "SpliceOutcomeSet",
        "Substitution",
        "Substitution",
        "Substitution",
    ]
    assert down_to_donor["count"] == len(
        direct.effects().filter_by_effect_priority(varcode.effects.SpliceDonor)
    )
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_filter_effects, world, min_priority_class="NoSuchClass", **kwargs)


def test_expression_effects_match_varcode_where_hashable(world: dict[str, Any]) -> None:
    kwargs = {"vcf_paths": [world["multi"]]}
    levels = {"T1": 3.0, "T3": 1.0}
    annotate = _call(vp.varcode_expression_effects, world, levels, mode="annotate", **kwargs)
    assert annotate["count"] == 5  # intergenic has no transcript
    assert {row["transcript_id"]: row["expression"] for row in annotate["effects"]} == levels
    top = _call(vp.varcode_expression_effects, world, levels, mode="top", **kwargs)
    assert (
        top["count"] == 1
        and top["effects"][0]["transcript_id"] == "T1"
        and top["effects"][0]["expression"] == 3.0
    )
    only_t1 = _call(
        vp.varcode_expression_effects,
        world,
        {"T1": 3.0},
        mode="annotate",
        variants=[_variant(world, "missense_a"), _variant(world, "missense_d")],
    )
    _, direct = _direct(world)
    subset = varcode.VariantCollection(
        [
            v
            for v in direct
            if v.start in (world["positions"]["missense_a"], world["positions"]["missense_d"])
        ]
    )
    expected = {
        str(e.short_description): level
        for e, level in subset.effects().effect_expression({"T1": 3.0}).items()
    }
    assert {row["description"]: row["expression"] for row in only_t1["effects"]} == expected
    unlisted = _call(vp.varcode_expression_effects, world, {"T9": 4.0}, mode="annotate", **kwargs)
    assert {row["expression"] for row in unlisted["effects"]} == {0.0}
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_expression_effects, world, levels, mode="rank", **kwargs)


def test_reconcile_splice_outcomes_with_rna(world: dict[str, Any]) -> None:
    donor = _variant(world, "donor_f")
    result = _call(
        vp.varcode_reconcile_splice_outcomes_with_rna,
        world,
        donor,
        "T3",
        [{"mechanism": "IntronRetention", "read_count": 9, "source": "star"}],
    )
    predicted = {row["effect"]["class"] for row in result["predicted"]}
    assert {"NormalSplicing", "IntronRetention", "CrypticDonor"} <= predicted
    assert [row["effect"]["class"] for row in result["candidates"]] == ["IntronRetention"]
    evidence = result["candidates"][0]["rna_evidence"]
    assert (
        len(evidence) == 1
        and evidence[0]["source"] == "star"
        and evidence[0]["evidence"] == {"read_count": 9}
    )
    assert {row["effect"]["class"] for row in result["excluded"]} == predicted - {"IntronRetention"}
    with pytest.raises(ConfigurationError, match="not predicted"):
        _call(
            vp.varcode_reconcile_splice_outcomes_with_rna,
            world,
            donor,
            "T3",
            [{"mechanism": "ExonSkipping"}],
        )
    with pytest.raises(ConfigurationError, match="splice signal"):
        _call(
            vp.varcode_reconcile_splice_outcomes_with_rna,
            world,
            _variant(world, "missense_a"),
            "T1",
            [{"mechanism": "IntronRetention"}],
        )
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_reconcile_splice_outcomes_with_rna, world, donor, "T3", [])


# ---- phasing ---------------------------------------------------------------------------------


def test_query_in_cis_vcf_and_rna(world: dict[str, Any]) -> None:
    a, b, c = (_variant(world, k) for k in ("missense_a", "missense_b", "silent_c"))
    vcf = _call(vp.varcode_query_in_cis, world, a, b, vcf_path=world["multi"], sample="S1")
    assert (vcf["in_cis"], vcf["relation"], vcf["phase_source"]) == (True, "cis", "vcf_ps")
    _, direct = _direct(world)
    resolver = varcode.VCFPhaseResolver(direct, "S1")
    variants = {v.start: v for v in direct}
    assert resolver.in_cis(variants[a["start"]], variants[c["start"]]) is True
    assert (
        _call(vp.varcode_query_in_cis, world, a, c, vcf_path=world["multi"], sample="S1")["in_cis"]
        is True
    )
    unphased = _call(
        vp.varcode_query_in_cis,
        world,
        a,
        _variant(world, "missense_d"),
        vcf_path=world["multi"],
        sample="S1",
    )
    assert unphased["in_cis"] is None and unphased["relation"] == "unknown"
    rna = _call(vp.varcode_query_in_cis, world, a, b, method="rna_reads", bam_path=world["bam"])
    assert (rna["in_cis"], rna["phase_source"]) == (True, "rna_reads")
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_query_in_cis, world, a, b, method="vcf")
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_query_in_cis, world, a, b, method="rna_reads")
    with pytest.raises(ConfigurationError):
        _call(
            vp.varcode_query_in_cis, world, a, b, method="vcf", vcf_path=world["multi"], sample="S9"
        )


def test_phased_partners_vcf_and_rna(world: dict[str, Any]) -> None:
    a, b, c, d = (
        _variant(world, k) for k in ("missense_a", "missense_b", "silent_c", "missense_d")
    )
    vcf = _call(vp.varcode_phased_partners, world, a, vcf_path=world["multi"], sample="S1")
    assert {row["start"] for row in vcf["partners"]} == {b["start"], c["start"]}
    rna = _call(
        vp.varcode_phased_partners,
        world,
        a,
        method="rna_reads",
        bam_path=world["bam"],
        candidates=[b, d],
    )
    assert [row["start"] for row in rna["partners"]] == [b["start"]]
    alone = _call(
        vp.varcode_phased_partners,
        world,
        d,
        method="rna_reads",
        bam_path=world["bam"],
        candidates=[a, b],
    )
    assert alone["partners"] == [] and alone["count"] == 0


def test_rna_read_support_counts_fragments(world: dict[str, Any]) -> None:
    a, b, d = (_variant(world, k) for k in ("missense_a", "missense_b", "missense_d"))
    result = _call(vp.varcode_rna_read_support, world, world["bam"], [a, b, d])
    assert [(row["alt_reads"], row["has_evidence"]) for row in result["variants"]] == [
        (6, True),
        (6, True),
        (0, False),
    ]
    strict = _call(vp.varcode_rna_read_support, world, world["bam"], [a], min_alt_reads=7)
    assert (
        strict["variants"][0]["alt_reads"] == 6 and strict["variants"][0]["has_evidence"] is False
    )
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_rna_read_support, world, world["bam"], [])
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_rna_read_support, world, world["bam"], [a], phasing_error_rate=2)


def test_rna_haplotype_test_needs_whole_context(world: dict[str, Any]) -> None:
    a, b, d = (_variant(world, k) for k in ("missense_a", "missense_b", "missense_d"))
    result = _call(vp.varcode_test_rna_haplotype, world, world["bam"], [a, b])
    assert [row["alt_reads"] for row in result["support"]] == [6, 6]
    assert result["pairs"][0]["in_cis"] is True
    other = _call(vp.varcode_test_rna_haplotype, world, world["bam"], [a, d])
    assert [row["alt_reads"] for row in other["support"]] == [0, 0]
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_test_rna_haplotype, world, world["bam"], [a, a])
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_test_rna_haplotype, world, world["bam"], [a], flanking_bases=0)


# ---- germline --------------------------------------------------------------------------------


def test_load_germline_context_single_and_multi_sample(world: dict[str, Any]) -> None:
    single = _call(vp.varcode_load_germline_context, world, world["single"])
    assert (
        single["count"] == 6 and single["completeness"] == "complete" and single["sample"] is None
    )
    multi = _call(
        vp.varcode_load_germline_context, world, world["multi"], sample="S2", completeness="sparse"
    )
    positions = world["positions"]
    # S2 does not carry silent_c (0/0), intergenic_e (0/0) and its missense_b call is missing
    assert {row["start"] for row in multi["variants"]} == {
        positions["missense_a"],
        positions["missense_d"],
        positions["donor_f"],
    }
    assert multi["completeness"] == "sparse"
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_load_germline_context, world, world["multi"], sample="S9")
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_load_germline_context, world, world["single"], completeness="partial")


def test_germline_variants_window_and_homozygous(world: dict[str, Any]) -> None:
    positions = world["positions"]
    low, high = positions["missense_a"], positions["silent_c"]
    window = _call(
        vp.varcode_germline_variants, world, world["single"], contig="chr1", start=low, end=high
    )
    assert {row["start"] for row in window["variants"]} == {
        positions["missense_a"],
        positions["missense_b"],
        positions["silent_c"],
    }
    edge = _call(
        vp.varcode_germline_variants,
        world,
        world["single"],
        contig="chr1",
        start=high + 1,
        end=high + 2,
    )
    assert edge["count"] == 0
    homozygous = _call(vp.varcode_germline_variants, world, world["single"], homozygous_only=True)
    assert [row["start"] for row in homozygous["variants"]] == [positions["silent_c"]]
    both = _call(
        vp.varcode_germline_variants,
        world,
        world["single"],
        contig="chr1",
        start=1,
        end=100,
        homozygous_only=True,
    )
    assert both["count"] == 0
    with pytest.raises(ConfigurationError, match="together"):
        _call(vp.varcode_germline_variants, world, world["single"], contig="chr1", start=1)
    with pytest.raises(ConfigurationError, match="single-sample"):
        _call(
            vp.varcode_germline_variants, world, world["multi"], sample="S1", homozygous_only=True
        )


def test_validate_germline_context(world: dict[str, Any]) -> None:
    ok = _call(
        vp.varcode_validate_germline_context,
        world,
        world["single"],
        {"variants": [_variant(world, "missense_a")]},
    )
    assert ok["valid"] is True and ok["error"] is None
    assert ok["germline_count"] == 6 and ok["somatic_count"] == 1
    empty = world["dir"] / "empty.vcf"
    empty.write_text(
        Path(world["single"]).read_text().split("#CHROM")[0]
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
        encoding="utf-8",
    )
    warned = _call(
        vp.varcode_validate_germline_context,
        world,
        str(empty),
        {"variants": [_variant(world, "missense_a")]},
    )
    assert warned["valid"] is True and any("zero variants" in text for text in warned["warnings"])
    with pytest.raises(ConfigurationError):
        _call(vp.varcode_validate_germline_context, world, world["single"], {})


# ---- structural variants ---------------------------------------------------------------------


def test_sv_junction_inserted_sequence(world: dict[str, Any]) -> None:
    result = _call(vp.varcode_sv_junction_inserted_sequence, world, [world["sv"]])
    by_type = {row["sv_type"]: row for row in result["junctions"]}
    assert (
        by_type["BND"]["inserted_sequence"] == "ACGT"
    )  # the ALT is the anchor base followed by ACGT
    assert by_type["DEL"]["inserted_sequence"] == ""
    assert [end["position"] for end in by_type["BND"]["junction"]] == [300, 700]
    skipped = _call(vp.varcode_sv_junction_inserted_sequence, world, [world["multi"]])
    assert skipped["count"] == 0
