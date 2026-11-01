"""Tests for the minimap2 and paftools.js command-line facades (pack P09)."""

from __future__ import annotations

import inspect
import itertools
import random
import shutil
from pathlib import Path

import pytest

from dnakit.alignment import minimap2_p09_cli_functions as mm
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

needs_minimap2 = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 is not installed"
)
needs_paftools = pytest.mark.skipif(
    shutil.which("k8") is None or shutil.which("paftools.js") is None,
    reason="k8 and paftools.js are not installed",
)

_EXONS = ((500, 800), (1500, 1800), (2500, 2800))


def _random(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(sequence: str, positions: tuple[int, ...]) -> str:
    letters = list(sequence)
    for position in positions:
        letters[position] = "A" if letters[position] != "A" else "C"
    return "".join(letters)


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """A small genome with one spliced gene, cDNA, an assembly with SNPs and indels, simulated short reads."""
    root = tmp_path_factory.mktemp("minimap2data")
    rng = random.Random(7)
    chr_a = list(_random(rng, 12000))
    chr_b = _random(rng, 8000)
    for (_, end), (start_next, _) in itertools.pairwise(_EXONS):
        chr_a[end], chr_a[end + 1] = "G", "T"
        chr_a[start_next - 2], chr_a[start_next - 1] = "A", "G"
    reference = "".join(chr_a)
    paths = {
        name: root / name
        for name in (
            "ref.fa",
            "ref.fa.fai",
            "cdna.fa",
            "genes.gtf",
            "genes.gff3",
            "asm.fa",
            "sim.fq",
            "queries.bed",
            "truth.txt",
        )
    }
    paths["ref.fa"].write_text(f">chrA\n{reference}\n>chrB\n{chr_b}\n")
    paths["ref.fa.fai"].write_text(
        f"chrA\t{len(reference)}\t6\t{len(reference)}\t{len(reference) + 1}\n"
        f"chrB\t{len(chr_b)}\t{len(reference) + 13}\t{len(chr_b)}\t{len(chr_b) + 1}\n"
    )
    cdna = "".join(reference[start:end] for start, end in _EXONS)
    paths["cdna.fa"].write_text(f">tx1\n{cdna}\n>tx1mut\n{_mutate(cdna, (100,))}\n")
    paths["genes.gtf"].write_text(
        "".join(
            f'chrA\ttest\texon\t{start + 1}\t{end}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            for start, end in _EXONS
        )
    )
    paths["genes.gff3"].write_text(
        "##gff-version 3\nchrA\ttest\tgene\t501\t2800\t.\t+\t.\tID=g1\nchrA\ttest\tmRNA\t501\t2800\t.\t+\t.\tID=t1;Parent=g1\n"
        + "".join(
            f"chrA\ttest\tCDS\t{start + 1}\t{end}\t.\t+\t0\tID=c{index};Parent=t1\n"
            for index, (start, end) in enumerate(_EXONS)
        )
    )
    assembly = _mutate(reference, (1000, 3000, 5000, 7000))
    assembly = assembly[:9000] + assembly[9010:]
    paths["asm.fa"].write_text(
        f">ctgA\n{assembly}\n>ctgB\n{chr_b[:4000]}ACGTACGTAC{chr_b[4000:]}\n"
    )
    reads, truth = [], []
    for index in range(20):
        start = rng.randrange(0, len(chr_b) - 150)
        name = f"r{index}!chrB!{start}!{start + 150}!+"
        reads.append(f"@{name}\n{chr_b[start : start + 150]}\n+\n{'I' * 150}\n")
        truth.append(name)
    paths["sim.fq"].write_text("".join(reads))
    paths["truth.txt"].write_text("\n".join(truth) + "\n")
    paths["queries.bed"].write_text("ctgA\t100\t200\nctgA\t9500\t9600\n")
    return paths


@pytest.fixture(scope="module")
def alignments(data: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> dict[str, str]:
    """Alignments produced by minimap2_map_reads (so the mapper wrapper is exercised by every dependent test)."""
    if shutil.which("minimap2") is None:
        pytest.skip("minimap2 is not installed")
    root = tmp_path_factory.mktemp("minimap2aln")
    out: dict[str, str] = {}
    for name, preset, target, query, extra in (
        (
            "splice",
            "splice",
            "ref.fa",
            "cdna.fa",
            {"base_level_alignment": True, "cs_tag": "short"},
        ),
        ("splice_sam", "splice", "ref.fa", "cdna.fa", {"output_format": "sam"}),
        ("asm", "asm5", "ref.fa", "asm.fa", {"base_level_alignment": True, "cs_tag": "short"}),
        ("asm_long", "asm5", "ref.fa", "asm.fa", {"base_level_alignment": True, "cs_tag": "long"}),
        ("sim", "sr", "ref.fa", "sim.fq", {"base_level_alignment": True}),
        ("sim_sam", "sr", "ref.fa", "sim.fq", {"output_format": "sam"}),
        ("ava", "ava-ont", "sim.fq", "sim.fq", {}),
    ):
        result = mm.minimap2_map_reads(
            data[target], data[query], root / name, preset=preset, verbose=False, **extra
        )
        out[name] = result["outputs"]["alignments"]
    return out


# ---------------------------------------------------------------------------
# No program needed
# ---------------------------------------------------------------------------


def test_missing_programs_create_nothing(tmp_path: Path, data: dict[str, Path]) -> None:
    missing = str(tmp_path / "absent")
    with pytest.raises(BackendUnavailableError):
        mm.minimap2_map_reads(
            data["ref.fa"], data["cdna.fa"], tmp_path / "o1", minimap2_path=missing, verbose=False
        )
    with pytest.raises(BackendUnavailableError):
        mm.minimap2_vcf_stats(
            data["genes.gtf"],
            tmp_path / "o2",
            paftools_path=missing,
            k8_path=missing,
            verbose=False,
        )
    with pytest.raises(BackendUnavailableError):
        mm.minimap2_cluster_overlap_reads(
            data["genes.gtf"], tmp_path / "o3", pafcluster_path=missing, verbose=False
        )
    assert (
        not (tmp_path / "o1").exists()
        and not (tmp_path / "o2").exists()
        and not (tmp_path / "o3").exists()
    )


def test_validation_errors(tmp_path: Path, data: dict[str, Path]) -> None:
    with pytest.raises(ConfigurationError):
        mm.minimap2_map_reads(
            data["ref.fa"], data["cdna.fa"], tmp_path / "o1", preset="nonsense", verbose=False
        )  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        mm.minimap2_convert_simulated_reads(
            "pbsim", [data["sim.fq"]], tmp_path / "o2", verbose=False
        )
    with pytest.raises(ConfigurationError):
        mm.minimap2_pair_haplotype_vcf(data["genes.gtf"], tmp_path / "o3", male=True, verbose=False)
    with pytest.raises(ConfigurationError):
        mm.minimap2_detect_misjoins(
            data["genes.gtf"], tmp_path / "o4", max_gap="lots", verbose=False
        )
    with pytest.raises(ConfigurationError):
        mm.minimap2_annotation_to_bed12(
            data["genes.gtf"], tmp_path / "o5", output_gff=True, verbose=False
        )


def test_public_names_have_verbose() -> None:
    for name in mm.__all__:
        signature = inspect.signature(getattr(mm, name))
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name


# ---------------------------------------------------------------------------
# minimap2
# ---------------------------------------------------------------------------


@needs_minimap2
def test_spliced_alignment_and_overlaps(tmp_path: Path, data: dict[str, Path]) -> None:
    spliced = mm.minimap2_map_reads(
        data["ref.fa"],
        data["cdna.fa"],
        tmp_path / "splice",
        preset="splice:hq",
        base_level_alignment=True,
        verbose=False,
    )
    first = spliced["records"][0]
    assert first["query"] == "tx1" and first["target"] == "chrA"
    assert (first["target_start"], first["target_end"]) == (500, 2800)
    assert first["tags"]["cg"] == "300M700N300M700N300M"
    sam = mm.minimap2_map_reads(
        data["ref.fa"],
        data["cdna.fa"],
        tmp_path / "splice_sam",
        preset="splice",
        output_format="sam",
        splice_strand="f",
        verbose=False,
    )
    assert sam["alignment_count"] == 2 and "records" not in sam
    overlap = mm.minimap2_map_reads(
        data["sim.fq"], data["sim.fq"], tmp_path / "ava", preset="ava-ont", verbose=False
    )
    assert overlap["alignment_count"] > 0
    assert overlap["version"].startswith("2.")


@needs_minimap2
def test_paired_short_reads(tmp_path: Path, data: dict[str, Path]) -> None:
    result = mm.minimap2_map_reads(
        data["ref.fa"],
        data["sim.fq"],
        tmp_path / "sr",
        preset="sr",
        mate=data["sim.fq"],
        output_format="sam",
        verbose=False,
    )
    assert result["alignment_count"] >= 20


# ---------------------------------------------------------------------------
# paftools.js
# ---------------------------------------------------------------------------


@needs_minimap2
@needs_paftools
def test_format_conversions(
    tmp_path: Path, data: dict[str, Path], alignments: dict[str, str]
) -> None:
    view = mm.minimap2_view_paf(
        alignments["asm_long"], tmp_path / "view", line_length=60, verbose=False
    )
    text = Path(view["outputs"]["output"]).read_text()
    assert "Ref+:" in text and "Qry+:" in text
    converted = mm.minimap2_sam_to_paf(
        alignments["splice_sam"], tmp_path / "sam2paf", alignments="primary", verbose=False
    )
    assert converted["row_count"] == 2 and converted["rows"][0][5] == "chrA"
    bed = mm.minimap2_spliced_alignments_to_bed12(
        alignments["splice"], tmp_path / "bed12", verbose=False
    )
    assert (
        bed["records"][0]["block_sizes"] == "300,300,300,"
        and bed["records"][0]["block_starts"] == "0,1000,2000,"
    )
    from_sam = mm.minimap2_spliced_alignments_to_bed12(
        alignments["splice_sam"], tmp_path / "bed12sam", verbose=False
    )
    assert from_sam["record_count"] == 2
    annotation = mm.minimap2_annotation_to_bed12(data["genes.gtf"], tmp_path / "gtf", verbose=False)
    assert annotation["rows"][0][:3] == ["chrA", "500", "2800"] and annotation["rows"][0][9] == "3"
    junction_bed = mm.minimap2_annotation_to_bed12(
        data["genes.gtf"], tmp_path / "gtfj", junctions_only=True, verbose=False
    )
    assert junction_bed["record_count"] >= 1
    junctions = mm.minimap2_gff3_to_junctions(data["genes.gff3"], tmp_path / "junc", verbose=False)
    assert [(row["start"], row["end"]) for row in junctions["records"]] == [
        (800, 1500),
        (1800, 2500),
    ]
    sequences = mm.minimap2_longcs_to_sequences(
        alignments["asm_long"], tmp_path / "longcs", verbose=False
    )
    assert sequences["record_count"] == 2
    queries = mm.minimap2_longcs_to_sequences(
        alignments["asm_long"], tmp_path / "longcs_q", query=True, verbose=False
    )
    assert Path(queries["outputs"]["output"]).read_text().startswith(">ctgA_0_11990")


@needs_minimap2
@needs_paftools
def test_delta_and_gff_failure(tmp_path: Path, data: dict[str, Path]) -> None:
    delta = tmp_path / "toy.delta"
    delta.write_text("ref.fa qry.fa\nNUCMER\n>chrX qryX 1000 900\n1 900 1 900 0 0 0\n0\n")
    result = mm.minimap2_delta_to_paf(delta, tmp_path / "delta", verbose=False)
    row = result["rows"][0]
    assert row[:12] == [
        "qryX",
        "900",
        "0",
        "900",
        "+",
        "chrX",
        "1000",
        "0",
        "900",
        "900",
        "900",
        "0",
    ]
    with pytest.raises(BackendExecutionError):
        mm.minimap2_annotation_to_bed12(data["genes.gff3"], tmp_path / "gff3", verbose=False)


@needs_minimap2
@needs_paftools
def test_statistics(tmp_path: Path, data: dict[str, Path], alignments: dict[str, str]) -> None:
    stats = mm.minimap2_alignment_stats(alignments["asm"], tmp_path / "stat", verbose=False)
    assert stats["metrics"]["substitutions"] == 4 and stats["metrics"]["insertions_in_0_50"] == 1
    assert stats["metrics"]["deletions_in_0_50"] == 1
    assembly = mm.minimap2_assembly_stats(
        data["ref.fa.fai"], [alignments["asm"]], tmp_path / "asmstat", verbose=False
    )
    assert assembly["metrics"]["Length"] == [20000] and assembly["metrics"]["NG50"] == [11990]
    genes = mm.minimap2_assembly_gene_completeness(
        alignments["splice"], [alignments["splice"]], tmp_path / "asmgene", verbose=False
    )
    assert genes["metrics"]["full_sgl"] == [2, 2]
    misjoin = mm.minimap2_detect_misjoins(alignments["asm"], tmp_path / "misjoin", verbose=False)
    assert misjoin["metrics"]["inter_chromosomal_misjoins"] == [0, 0]
    vcf_stats = mm.minimap2_vcf_stats(
        _write_vcf(tmp_path / "toy.vcf"), tmp_path / "vcfstat", verbose=False
    )
    assert vcf_stats["metrics"]["substitutions"] == 2 and vcf_stats["metrics"]["insertions"] == 1
    assert vcf_stats["metrics"]["deletions"] == 1


def _write_vcf(path: Path) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\nchr1\t200\t.\tC\tA\t.\t.\t.\n"
        "chr1\t300\t.\tT\tTAAA\t.\t.\t.\nchr1\t400\t.\tGCCC\tG\t.\t.\t.\n"
    )
    return path


@needs_minimap2
@needs_paftools
def test_liftover_and_calling(
    tmp_path: Path, data: dict[str, Path], alignments: dict[str, str]
) -> None:
    lifted = mm.minimap2_liftover_coordinates(
        alignments["asm"],
        data["queries.bed"],
        tmp_path / "lift",
        min_alignment_length=1000,
        verbose=False,
    )
    assert lifted["rows"][0][:3] == ["chrA", "100", "200"]
    variants = mm.minimap2_call_assembly_variants(
        alignments["asm"],
        tmp_path / "call",
        reference_fasta=data["ref.fa"],
        min_coverage_length=1000,
        min_call_length=1000,
        verbose=False,
    )
    records = [
        line.split("\t")
        for line in Path(variants["outputs"]["output"]).read_text().splitlines()
        if not line.startswith("#")
    ]
    assert sum(1 for row in records if len(row[3]) == 1 and len(row[4]) == 1) == 4
    assert variants["summary"]["substitutions"] == 4
    text_only = mm.minimap2_call_assembly_variants(
        alignments["asm"],
        tmp_path / "call_txt",
        min_coverage_length=1000,
        min_call_length=1000,
        verbose=False,
    )
    assert Path(text_only["outputs"]["output"]).name == "variants.txt"


@needs_minimap2
@needs_paftools
def test_vcf_tools(tmp_path: Path) -> None:
    vcf = _write_vcf(tmp_path / "toy.vcf")
    selected = mm.minimap2_select_vcf_by_length(vcf, tmp_path / "sel", min_length=2, verbose=False)
    assert selected["record_count"] == 2
    pair = tmp_path / "pair.vcf"
    pair.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thap1\thap2\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,5\t0/0:5,0\n"
    )
    paired = mm.minimap2_pair_haplotype_vcf(
        pair, tmp_path / "pair", sample_name="demo", verbose=False
    )
    body = Path(paired["outputs"]["output"]).read_text()
    assert paired["record_count"] == 1 and "\t1|0:5,5" in body and "demo" in body
    base = tmp_path / "base.vcf"
    called = tmp_path / "called.vcf"
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    ref = "A" + "C" * 99
    base.write_text(
        header + f"chr1\t1000\t.\tA\t{'A' + 'G' * 99}\t.\t.\t.\nchr1\t5000\t.\t{ref}\tA\t.\t.\t.\n"
    )
    called.write_text(header + f"chr1\t1010\t.\tA\t{'A' + 'G' * 99}\t.\t.\t.\n")
    compared = mm.minimap2_compare_sv_callsets(base, called, tmp_path / "sveval", verbose=False)
    assert compared["row_count"] >= 1


@needs_minimap2
@needs_paftools
def test_simulation_and_evaluation(
    tmp_path: Path, data: dict[str, Path], alignments: dict[str, str]
) -> None:
    truth = mm.minimap2_simulated_names_to_bed(
        data["truth.txt"], tmp_path / "sim2bed", verbose=False
    )
    assert truth["record_count"] == 20 and truth["records"][0]["chrom"] == "chrB"
    evaluation = mm.minimap2_evaluate_mapping(
        alignments["sim"], tmp_path / "mapeval", verbose=False
    )
    assert (
        evaluation["records"][-1]["cumulative_reads"] == 20
        and evaluation["records"][-1]["wrong"] == 0
    )
    base, test = tmp_path / "base.paf", tmp_path / "test.paf"

    def line(name: str, start: int) -> str:
        return f"{name}\t6000\t0\t6000\t+\tchr1\t1000000\t{start}\t{start + 6000}\t5900\t6000\t60\n"

    base.write_text(line("q1", 100) + line("q2", 20000) + line("q3", 40000))
    test.write_text(line("q1", 150) + line("q3", 90000))
    compared = mm.minimap2_compare_paf(base, test, tmp_path / "pafcmp", verbose=False)
    assert compared["metrics"] == {
        "base_alignments_with_mapq_10": 3,
        "base_alignments_correctly_mapped_by_test": 1,
        "wrong_test_alignment": 1,
        "base_alignments_missing": 1,
        "additional_test_alignments_with_mapq_10": 0,
    }
    junctions = mm.minimap2_evaluate_junctions(
        data["genes.gtf"], alignments["splice_sam"], tmp_path / "junceval", verbose=False
    )
    assert (
        junctions["metrics"]["predicted_introns"] == 4
        and junctions["metrics"]["correct_introns_percent"] == 100.0
    )
    exons = mm.minimap2_evaluate_exons(
        data["genes.gtf"], alignments["splice_sam"], tmp_path / "exoneval", verbose=False
    )
    assert exons["metrics"]["predicted_exons"] == 6 and exons["metrics"]["correct_exons"] == 6
    overlap_eval = mm.minimap2_evaluate_overlaps(
        alignments["sim"],
        alignments["ava"],
        tmp_path / "oveval",
        min_overlap_length=50,
        verbose=False,
    )
    assert (
        overlap_eval["metrics"]["inferred_overlaps"] >= 1
        and "sensitivity_percent" in overlap_eval["metrics"]
    )


@needs_minimap2
@needs_paftools
def test_simulator_conversions(tmp_path: Path, data: dict[str, Path]) -> None:
    mason = tmp_path / "mason.sam"
    mason.write_text(
        "@HD\tVN:1.6\n"
        "simulated.r1\t0\tchrB\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXE:i:0\tXS:i:0\tXI:i:0\n"
        "simulated.r2\t16\tchrB\t201\t60\t10M\t*\t0\t0\tAAAACCCCGG\tABCDEFGHIJ\tXE:i:1\tXS:i:1\tXI:i:0\n"
    )
    converted = mm.minimap2_convert_simulated_reads(
        "mason2", [mason], tmp_path / "mason", verbose=False
    )
    lines = Path(converted["outputs"]["output"]).read_text().splitlines()
    assert converted["record_count"] == 2 and lines[0] == "@r1!chrB!100!110!+ 0:0:0"
    assert lines[5] == "CCGGGGTTTT" and lines[7] == "JIHGFEDCBA"
    maf = tmp_path / "pbsim.maf"
    maf.write_text("a\ns ref 100 8 + 5000 ACGTACGT\ns S1_1 0 8 - 8 ACGT-ACGT\n\n")
    pbsim = mm.minimap2_convert_simulated_reads(
        "pbsim", [maf], tmp_path / "pbsim", reference_fai=data["ref.fa.fai"], verbose=False
    )
    assert Path(pbsim["outputs"]["output"]).read_text().splitlines()[0] == ">S1_1!chrA!100!108!-"
    badread = tmp_path / "badread.fq"
    badread.write_text(
        "@id1 chrA,+strand,100-108 length=8 read_identity=91.5%\nACGTACGT\n+\nIIIIIIII\n"
    )
    converted_badread = mm.minimap2_convert_simulated_reads(
        "badread", [badread], tmp_path / "badread", reference_fai=data["ref.fa.fai"], verbose=False
    )
    assert (
        Path(converted_badread["outputs"]["output"])
        .read_text()
        .startswith(">S1!chrA!100!108!+\tri:f:91.5")
    )


@needs_minimap2
@needs_paftools
def test_protein_models_and_clusters(
    tmp_path: Path, data: dict[str, Path], alignments: dict[str, str]
) -> None:
    protein_paf = tmp_path / "miniprot.paf"
    protein_paf.write_text(
        "prot1\t100\t0\t100\t+\tchrA\t12000\t500\t800\t290\t300\t60\tAS:i:100\tms:i:100\tcg:Z:300M\n"
    )
    models = mm.minimap2_paf_to_gff3(protein_paf, tmp_path / "gff", verbose=False)
    lines = Path(models["outputs"]["output"]).read_text().splitlines()
    assert lines[0].split("\t")[2] == "transcript" and lines[1].split("\t")[2] == "CDS"
    assert lines[0].split("\t")[3:5] == ["501", "800"]
    finder = shutil.which("pafcluster.js")
    if finder is not None:
        clusters = mm.minimap2_cluster_overlap_reads(
            alignments["ava"], tmp_path / "cluster", verbose=False
        )
        assert clusters["status"] == "completed"
