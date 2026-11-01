"""Tests for the Ensembl VEP command-line facades (pack P09)."""

from __future__ import annotations

import inspect
import json
import os
import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.variants import ensembl_vep_p09_cli_functions as vep

needs_vep = pytest.mark.skipif(
    shutil.which("vep") is None or shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="vep, bgzip and tabix are not installed",
)


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """One two-exon gene on a 3 kb contig, a bgzipped GFF3, a VCF of four variants and a phased VCF."""
    root = tmp_path_factory.mktemp("vepdata")
    rng = random.Random(5)
    codons = [
        a + b + c
        for a in "ACGT"
        for b in "ACGT"
        for c in "ACGT"
        if a + b + c not in ("TAA", "TAG", "TGA")
    ]
    cds = "".join(
        ["ATG"]
        + [rng.choice(codons) for _ in range(49)]
        + [rng.choice(codons) for _ in range(49)]
        + ["TAA"]
    )
    genome = list("".join(rng.choice("ACGT") for _ in range(3000)))
    genome[500:650] = cds[:150]
    genome[1000:1150] = cds[150:]
    genome[650], genome[651], genome[998], genome[999] = "G", "T", "A", "G"
    sequence = "".join(genome)
    paths = {
        name: root / name
        for name in (
            "genome.fa",
            "genes.gff3",
            "genes.gff3.gz",
            "variants.vcf",
            "phased.vcf",
            "input.txt",
        )
    }
    paths["genome.fa"].write_text(
        ">chr1\n" + "\n".join(sequence[i : i + 60] for i in range(0, 3000, 60)) + "\n"
    )
    features = [
        ("gene", 501, 1150, "ID=gene1;Name=GENE1;biotype=protein_coding"),
        ("mRNA", 501, 1150, "ID=tx1;Parent=gene1;Name=GENE1-201;biotype=protein_coding"),
        ("exon", 501, 650, "ID=exon1;Parent=tx1"),
        ("exon", 1001, 1150, "ID=exon2;Parent=tx1"),
        ("CDS", 501, 650, "ID=cds1;Parent=tx1"),
        ("CDS", 1001, 1150, "ID=cds2;Parent=tx1"),
    ]
    lines = [
        f"chr1\ttest\t{kind}\t{start}\t{end}\t.\t+\t{'0' if kind == 'CDS' else '.'}\t{attributes}"
        for kind, start, end, attributes in features
    ]
    paths["genes.gff3"].write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
    header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=3000>\n"
    rows = []
    for position in (521, 530, 700, 2000):
        ref = sequence[position - 1]
        rows.append(f"chr1\t{position}\t.\t{ref}\t{'A' if ref != 'A' else 'C'}\t.\t.\t.")
    paths["variants.vcf"].write_text(
        header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(rows) + "\n"
    )
    phased = []
    for position, genotype in ((522, "0|1"), (531, "1|1"), (540, "0|1")):
        ref = sequence[position - 1]
        phased.append(
            f"chr1\t{position}\t.\t{ref}\t{'A' if ref != 'A' else 'C'}\t.\t.\t.\tGT\t{genotype}"
        )
    paths["phased.vcf"].write_text(
        header
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
        + "\n".join(phased)
        + "\n"
    )
    paths["input.txt"].write_text("5 140532 140532 T/C +\n5 140540 140540 A/G +\n")
    if shutil.which("bgzip") and shutil.which("tabix"):
        import subprocess

        body = [line for line in lines]
        body.sort(
            key=lambda line: (
                line.split("\t")[0],
                int(line.split("\t")[3]),
                int(line.split("\t")[4]),
            )
        )
        plain = root / "sorted.gff3"
        plain.write_text("##gff-version 3\n" + "\n".join(body) + "\n")
        subprocess.run(
            ["bgzip", "-c", str(plain)], stdout=paths["genes.gff3.gz"].open("wb"), check=True
        )
        subprocess.run(["tabix", "-p", "gff", str(paths["genes.gff3.gz"])], check=True)
    return paths


def _fake_vep(tmp_path: Path) -> Path:
    """A stand-in program that records its arguments and writes a minimal VCF, to test option mapping."""
    script = tmp_path / "vep"
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        "args = sys.argv[1:]\n"
        "if args == ['--help']:\n"
        "    print('  ensembl-vep          : 116.3')\n"
        "    sys.exit(0)\n"
        "out = args[args.index('--output_file') + 1]\n"
        "open(out, 'w').write('##INFO=<ID=CSQ,Number=.,Type=String,Description=\"x. Format: Allele|Consequence\">\\n"
        "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr1\\t5\\t.\\tA\\tC\\t.\\t.\\tCSQ=C|missense_variant\\n')\n"
        "open(out + '.args.json', 'w').write(json.dumps(args))\n"
        "if '--stats_file' in args:\n    open(args[args.index('--stats_file') + 1], 'w').write('stats')\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


# ---------------------------------------------------------------------------
# No VEP needed
# ---------------------------------------------------------------------------


def test_missing_programs_create_nothing(tmp_path: Path, data: dict[str, Path]) -> None:
    absent = str(tmp_path / "absent")
    for index, call in enumerate(
        (
            lambda target: vep.vep_annotate_variants(
                data["variants.vcf"], target, vep_path=absent, verbose=False
            ),
            lambda target: vep.vep_filter_results(
                data["variants.vcf"],
                target,
                filters=["x is y"],
                filter_vep_path=absent,
                verbose=False,
            ),
            lambda target: vep.vep_predict_haplotypes(
                data["phased.vcf"], target, haplo_path=absent, verbose=False
            ),
        )
    ):
        target = tmp_path / f"out{index}"
        with pytest.raises(BackendUnavailableError):
            call(target)
        assert not target.exists()
    with pytest.raises(BackendUnavailableError):
        vep.vep_validate_input_format(
            data["input.txt"], validator_path=tmp_path / "absent.py", verbose=False
        )
    with pytest.raises(BackendUnavailableError):
        vep.vep_convert_cache(
            tmp_path,
            species="homo_sapiens",
            version="all",
            convert_cache_path=absent,
            verbose=False,
        )


def test_validation_errors(tmp_path: Path, data: dict[str, Path]) -> None:
    fake = str(_fake_vep(tmp_path))
    with pytest.raises(ConfigurationError):
        vep.vep_annotate_variants(
            data["variants.vcf"], tmp_path / "o1", vep_path=fake, phylop=["100"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        vep.vep_annotate_variants(
            data["variants.vcf"],
            tmp_path / "o2",
            vep_path=fake,
            phylop=["3"],
            ucsc_assembly="hg38",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        vep.vep_annotate_variants(
            data["variants.vcf"], tmp_path / "o3", vep_path=fake, shift_3prime=2, verbose=False
        )
    with pytest.raises(ConfigurationError):
        vep.vep_annotate_variants(
            data["variants.vcf"],
            tmp_path / "o4",
            vep_path=fake,
            custom_annotations=[{"file": data["variants.vcf"]}],
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        vep.vep_filter_results(data["variants.vcf"], tmp_path / "o5", verbose=False)
    with pytest.raises(ConfigurationError):
        vep.vep_filter_results(
            data["variants.vcf"], tmp_path / "o6", filters=["x is y; rm -rf /"], verbose=False
        )
    for name in ("o1", "o2", "o3", "o4", "o5", "o6"):
        assert not (tmp_path / name).exists()


def test_public_names_have_verbose() -> None:
    for name in vep.__all__:
        assert (
            inspect.signature(getattr(vep, name)).parameters["verbose"].kind
            is inspect.Parameter.KEYWORD_ONLY
        ), name


def test_option_mapping_with_stand_in_program(tmp_path: Path, data: dict[str, Path]) -> None:
    """The flags that need a cache or a database are only checked for the command line they produce."""
    fake = _fake_vep(tmp_path)
    result = vep.vep_annotate_variants(
        data["variants.vcf"],
        tmp_path / "run",
        vep_path=str(fake),
        output_format="tab",
        offline=True,
        summary=True,
        shift_hgvs=0,
        shift_3prime=1,
        shift_genomic=1,
        check_ref=True,
        lookup_ref=True,
        individual=["s1", "s2"],
        individual_zyg="all",
        check_existing=True,
        af=True,
        af_1kg=True,
        af_gnomade=True,
        af_gnomadg=True,
        max_af=True,
        check_frequency=True,
        freq_pop="1KG_ALL",
        freq_freq=0.01,
        freq_gt_lt="gt",
        freq_filter="exclude",
        filter_common=True,
        clin_sig_allele=1,
        clinvar_somatic_classification=True,
        gene_phenotype=True,
        mirna=True,
        regulatory=True,
        cell_type="HepG2",
        phylop=["100"],
        phastcons=["46"],
        ucsc_assembly="hg38",
        domains=True,
        nearest="symbol",
        stats="html",
        custom_annotations=[
            {
                "file": data["variants.vcf"],
                "short_name": "mine",
                "format": "vcf",
                "type": "overlap",
                "coverage": True,
                "fields": ["AF", "AC"],
            }
        ],
        verbose=False,
    )
    arguments = json.loads(Path(result["outputs"]["output"] + ".args.json").read_text())
    pairs = {
        arguments[i]: arguments[i + 1]
        for i in range(len(arguments) - 1)
        if arguments[i].startswith("--")
    }
    assert (
        pairs["--shift_hgvs"] == "0"
        and pairs["--shift_3prime"] == "1"
        and pairs["--shift_genomic"] == "1"
    )
    assert pairs["--individual"] == "s1,s2" and pairs["--individual_zyg"] == "all"
    assert (
        pairs["--freq_pop"] == "1KG_ALL"
        and pairs["--freq_freq"] == "0.01"
        and pairs["--freq_gt_lt"] == "gt"
    )
    assert (
        pairs["--freq_filter"] == "exclude"
        and pairs["--clin_sig_allele"] == "1"
        and pairs["--cell_type"] == "HepG2"
    )
    assert (
        pairs["--phyloP"] == "100"
        and pairs["--phastCons"] == "46"
        and pairs["--ucsc_assembly"] == "hg38"
    )
    assert pairs["--nearest"] == "symbol"
    assert pairs["--custom"].endswith(",mine,vcf,overlap,1,AF%AC")
    for flag in (
        "--summary",
        "--check_ref",
        "--lookup_ref",
        "--check_existing",
        "--af",
        "--af_1kg",
        "--af_gnomade",
        "--af_gnomadg",
        "--max_af",
        "--check_frequency",
        "--filter_common",
        "--clinvar_somatic_classification",
        "--gene_phenotype",
        "--mirna",
        "--regulatory",
        "--domains",
        "--tab",
        "--offline",
        "--stats_html",
    ):
        assert flag in arguments, flag
    assert "--no_stats" not in arguments
    assert result["columns"] == ["Allele", "Consequence"] or result["output_format"] == "tab"
    assert result["version"] == "116.3"


# ---------------------------------------------------------------------------
# With VEP installed (offline, annotation from a bgzipped GFF3)
# ---------------------------------------------------------------------------


@needs_vep
def test_annotate_all_output_formats(tmp_path: Path, data: dict[str, Path]) -> None:
    common = {"gff": data["genes.gff3.gz"], "fasta": data["genome.fa"], "verbose": False}
    as_vcf = vep.vep_annotate_variants(data["variants.vcf"], tmp_path / "vcf", **common)
    assert as_vcf["record_count"] == 4 and "Consequence" in as_vcf["columns"]
    by_position = {row["pos"]: row["annotations"][0]["Consequence"] for row in as_vcf["records"]}
    assert by_position == {
        521: "synonymous_variant",
        530: "synonymous_variant",
        700: "intron_variant",
        2000: "downstream_gene_variant",
    }
    as_tab = vep.vep_annotate_variants(
        data["variants.vcf"],
        tmp_path / "tab",
        output_format="tab",
        variant_class=True,
        hgvs=True,
        hgvsg=True,
        spdi=True,
        numbers=True,
        nearest="symbol",
        shift_3prime=1,
        check_ref=True,
        **common,
    )
    first = as_tab["records"][0]
    assert (
        first["VARIANT_CLASS"] == "SNV"
        and first["EXON"] == "1/2"
        and first["HGVSg"] == "chr1:g.521A>C"
    )
    assert first["HGVSc"] == "tx1:c.21A>C" and first["SPDI"] == "chr1:520:A:C"
    as_json = vep.vep_annotate_variants(
        data["variants.vcf"], tmp_path / "json", output_format="json", stats="text", **common
    )
    assert as_json["records"][0]["transcript_consequences"][0]["transcript_id"] == "tx1"
    assert Path(as_json["outputs"]["stats"]).is_file() and as_json["version"] is not None
    summary = vep.vep_annotate_variants(
        data["variants.vcf"], tmp_path / "summary", output_format="tab", summary=True, **common
    )
    assert {row["Consequence"] for row in summary["records"]} == {
        "synonymous_variant",
        "intron_variant",
        "downstream_gene_variant",
    }
    with pytest.raises(BackendExecutionError):
        vep.vep_annotate_variants(
            data["variants.vcf"],
            tmp_path / "clash",
            output_format="tab",
            summary=True,
            numbers=True,
            **common,
        )


@needs_vep
def test_filter_results(tmp_path: Path, data: dict[str, Path]) -> None:
    annotated = vep.vep_annotate_variants(
        data["variants.vcf"],
        tmp_path / "annotated",
        gff=data["genes.gff3.gz"],
        fasta=data["genome.fa"],
        variant_class=True,
        verbose=False,
    )
    path = annotated["outputs"]["output"]
    columns = vep.vep_filter_results(path, tmp_path / "list", list_fields=True, verbose=False)
    assert "Consequence" in columns["fields"] and "IMPACT" in columns["fields"]
    kept = vep.vep_filter_results(
        path, tmp_path / "kept", filters=["Consequence is intron_variant"], verbose=False
    )
    assert kept["record_count"] == 1 and kept["records"][0][1] == "700"
    either = vep.vep_filter_results(
        path, tmp_path / "either", filters=["IMPACT is LOW or IMPACT is MODIFIER"], verbose=False
    )
    assert either["record_count"] == 4
    counted = vep.vep_filter_results(
        path, tmp_path / "count", filters=["IMPACT is LOW"], count_only=True, verbose=False
    )
    assert counted["count"] == 2


@needs_vep
def test_predict_haplotypes(tmp_path: Path, data: dict[str, Path]) -> None:
    result = vep.vep_predict_haplotypes(
        data["phased.vcf"],
        tmp_path / "haplo",
        gff=data["genes.gff3.gz"],
        fasta=data["genome.fa"],
        verbose=False,
    )
    assert result["record_count"] == 1
    haplotypes = result["records"][0]
    assert haplotypes["transcript_id"] == "tx1" and len(haplotypes["protein_haplotypes"]) >= 2
    names = {item["name"] for item in haplotypes["protein_haplotypes"]}
    assert any(name.startswith("tx1.p:") for name in names)


def test_validate_input_format(tmp_path: Path, data: dict[str, Path]) -> None:
    script = os.environ.get("VEP_INPUT_VALIDATOR") or shutil.which("vep_input_format_validator.py")
    if script is None:
        pytest.skip("vep_input_format_validator.py is not available (set VEP_INPUT_VALIDATOR)")
    good = vep.vep_validate_input_format(data["input.txt"], validator_path=script, verbose=False)
    assert good["valid"] is True and good["error_count"] == 0
    bad = tmp_path / "bad.txt"
    bad.write_text("5 140532 140532 T/C +\n5 140540 140540 A/G f\n3 100 100 T +\n")
    report = vep.vep_validate_input_format(bad, validator_path=script, verbose=False)
    assert report["valid"] is False and report["error_count"] == 2
    assert report["errors"][0]["line"] == 2 and "Strand" in report["errors"][0]["message"]


@pytest.mark.skipif(
    shutil.which("vep_convert_cache") is None and shutil.which("convert_cache.pl") is None,
    reason="convert_cache is not installed",
)
def test_convert_cache_on_an_empty_cache(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    (cache / "homo_sapiens" / "116_GRCh38" / "1").mkdir(parents=True)
    result = vep.vep_convert_cache(
        cache, species="homo_sapiens", version="116_GRCh38", verbose=False
    )
    assert result["status"] == "completed" and "All done" in result["log_tail"]
