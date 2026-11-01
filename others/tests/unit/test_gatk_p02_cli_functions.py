"""Tests of the GATK P02 facades: argument mapping with a fake launcher, real GATK runs when it is installed."""
# ruff: noqa: E501

from __future__ import annotations

import gzip
import json
import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import gatk_p02_cli_functions as gp

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">\n'
    '##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">\n'
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">\n'
    '##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "##contig=<ID=chr1,length=3000>\n"
    "##contig=<ID=chr2,length=2000>\n"
)

FAKE_LAUNCHER = """#!{python}
import gzip, json, os, sys

args = sys.argv[1:]
if args[:1] == ["--version"]:
    print("The Genome Analysis Toolkit (GATK) v4.7.0.0")
    sys.exit(0)
with open(os.environ["FAKE_GATK_LOG"], "a") as log:
    log.write(json.dumps(args) + "\\n")
if os.environ.get("FAKE_GATK_FAIL"):
    print("A USER ERROR has occurred: fake failure")
    sys.exit(2)
if args[0] == "--java-options":
    args = args[2:]
tool = args[0]
BODY = "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr1\\t5\\t.\\tA\\tC\\t9\\tPASS\\t.\\n"


def write(path, text=BODY):
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        with open(path, "w") as handle:
            handle.write(text)


outputs = {{"-O", "--OUTPUT", "--REJECT", "--SNP_OUTPUT", "--INDEL_OUTPUT", "--filtering-stats", "--bam-output", "--output-model"}}
if tool == "VariantRecalibrator":
    outputs.add("--tranches-file")
prefix_tools = {{
    "ExtractVariantAnnotations": ["{{p}}.annot.hdf5", "{{p}}.vcf.gz"],
    "TrainVariantAnnotationsModel": ["{{p}}.snp.scorer.pkl"],
    "ScoreVariantAnnotations": ["{{p}}.vcf.gz", "{{p}}.scores.hdf5"],
}}
for index, token in enumerate(args):
    if token in outputs:
        target = args[index + 1]
        if tool in prefix_tools:
            for pattern in prefix_tools[tool]:
                write(pattern.format(p=target)) if pattern.endswith(".gz") else open(pattern.format(p=target), "w").close()
        elif tool == "VariantEval":
            with open(target, "w") as handle:
                handle.write(
                    "#:GATKReport.v1.1:1\\n#:GATKTable:2:1:%s:%s:;\\n#:GATKTable:CountVariants:Counts\\n"
                    "CountVariants  nSNPs\\nCountVariants  3\\n\\n"
                )
        elif tool == "VcfToIntervalList":
            with open(target, "w") as handle:
                handle.write("@HD\\tVN:1.6\\nchr1\\t10\\t10\\t+\\tinterval-1\\n")
        else:
            write(target)
"""


def _vcf(path: Path, records: int = 3) -> Path:
    lines = [VCF_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"]
    for index in range(records):
        lines.append(f"chr1\t{10 + index * 7}\t.\tA\tC\t50\t.\tQD=12.0;DP=30\tGT\t0/1\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _file(path: Path, text: str = "x") -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _write_reference(path: Path, sequences: dict[str, str]) -> None:
    """Write a FASTA with 60-column lines and its .fai index (GATK needs both plus a .dict)."""
    offset = 0
    index_lines = []
    with path.open("w") as handle:
        for name, sequence in sequences.items():
            header = f">{name}\n"
            handle.write(header)
            offset += len(header)
            index_lines.append(f"{name}\t{len(sequence)}\t{offset}\t60\t61")
            for start in range(0, len(sequence), 60):
                line = sequence[start : start + 60] + "\n"
                handle.write(line)
                offset += len(line)
    Path(f"{path}.fai").write_text("\n".join(index_lines) + "\n")


@pytest.fixture
def fake(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    launcher = tmp_path / "bin" / "gatk"
    launcher.parent.mkdir()
    launcher.write_text(FAKE_LAUNCHER.format(python=sys.executable), encoding="utf-8")
    launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "argv.log"
    monkeypatch.setenv("FAKE_GATK_LOG", str(log))
    work = tmp_path / "work"
    work.mkdir()

    class Fake:
        path = str(launcher)
        folder = work

        @staticmethod
        def calls() -> list[list[str]]:
            return [json.loads(line) for line in log.read_text().splitlines()]

        @staticmethod
        def last() -> list[str]:
            return Fake.calls()[-1]

    return Fake


# ---------------------------------------------------------------------------
# Argument mapping with a fake launcher
# ---------------------------------------------------------------------------


def test_all_public_functions_have_the_verbosity_flag_and_gatk_path() -> None:
    import inspect

    assert len(gp.__all__) == 24
    for name in gp.__all__:
        parameters = inspect.signature(getattr(gp, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY
        assert parameters["gatk_path"].default is None


def test_missing_gatk_raises_backend_unavailable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PATH", str(tmp_path))
    vcf = _vcf(tmp_path / "in.vcf")
    with pytest.raises(BackendUnavailableError, match="gatk"):
        gp.fix_vcf_header_gatk(str(vcf), str(tmp_path / "out.vcf"), verbose=False)
    assert not (tmp_path / "out.vcf").exists()


def test_validate_vcf_maps_arguments_and_reports_verdict(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    ref = _file(fake.folder / "ref.fa")
    result = gp.validate_vcf_gatk(
        str(vcf),
        reference=str(ref),
        validate_gvcf=True,
        exclude_validation_types=["IDS"],
        intervals=["chr1:1-100"],
        gatk_path=fake.path,
        java_max_heap_gb=2,
        verbose=False,
    )
    argv = fake.last()
    assert argv[:3] == ["--java-options", "-Xmx2g", "ValidateVariants"]
    assert ["-V", str(vcf)] == argv[3:5]
    assert "--validate-GVCF" in argv and argv[argv.index("--validate-GVCF") + 1] == "true"
    assert argv[argv.index("--validation-type-to-exclude") + 1] == "IDS"
    assert argv[argv.index("-L") + 1] == "chr1:1-100"
    assert result["status"] == "completed" and result["valid"] is True
    assert result["version"] == "4.7.0.0" and result["tool"] == "ValidateVariants"


def test_validate_vcf_failure_is_reported_not_raised(fake, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("FAKE_GATK_FAIL", "1")
    vcf = _vcf(fake.folder / "in.vcf")
    result = gp.validate_vcf_gatk(str(vcf), gatk_path=fake.path, verbose=False)
    assert result["status"] == "failed" and result["valid"] is False
    assert (
        any("USER ERROR" in line for line in result["problems"])
        or "USER ERROR" in result["log_tail"]
    )


def test_evaluate_variants_parses_report_and_inverts_standard_module_switches(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    ref = _file(fake.folder / "ref.fa")
    out = fake.folder / "eval.report"
    result = gp.evaluate_variants_gatk(
        [str(vcf)],
        str(ref),
        str(out),
        eval_modules=["CountVariants"],
        use_all_standard_modules=False,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[argv.index("--do-not-use-all-standard-modules") + 1] == "true"
    assert "--do-not-use-all-standard-stratifications" not in argv
    assert argv[argv.index("-EV") + 1] == "CountVariants"
    assert result["tables"]["CountVariants"]["columns"] == ["CountVariants", "nSNPs"]
    assert result["tables"]["CountVariants"]["rows"] == [["CountVariants", 3]]


def test_vqsr_resources_become_tagged_feature_inputs(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    truth = _vcf(fake.folder / "truth.vcf")
    result = gp.recalibrate_variants_vqsr(
        [str(vcf)],
        [
            {
                "path": str(truth),
                "name": "hapmap",
                "tags": {"known": False, "training": True, "truth": True, "prior": 15.0},
            }
        ],
        ["QD", "FS"],
        str(fake.folder / "x.recal"),
        str(fake.folder / "x.tranches"),
        mode="BOTH",
        truth_sensitivity_tranches=[100.0, 99.0],
        max_gaussians=4,
        use_allele_specific_annotations=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert "--resource:hapmap,known=false,training=true,truth=true,prior=15.0" in argv
    assert argv[argv.index("--mode") + 1] == "BOTH"
    assert argv.count("-an") == 2 and argv.count("-tranche") == 2
    assert argv[argv.index("-AS") + 1] == "true"
    assert result["status"] == "completed"
    assert result["recal_file"].endswith("x.recal") and result["tranches_file"].endswith(
        "x.tranches"
    )


def test_apply_vqsr_and_tranche_filtering_arguments(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    recal = _file(fake.folder / "r.recal")
    tranches = _file(fake.folder / "r.tranches")
    truth = _vcf(fake.folder / "truth.vcf")
    result = gp.apply_vqsr(
        str(vcf),
        str(recal),
        str(fake.folder / "o.vcf"),
        tranches_file=str(tranches),
        truth_sensitivity_filter_level=99.5,
        exclude_filtered=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "ApplyVQSR" and "--exclude-filtered" in argv
    assert argv[argv.index("--truth-sensitivity-filter-level") + 1] == "99.5"
    assert result["summary"]["records"] == 1
    gp.filter_variant_tranches(
        str(vcf),
        [str(truth)],
        str(fake.folder / "t.vcf"),
        snp_tranches=[99.9],
        indel_tranches=[99.0],
        info_key="CNN_1D",
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "FilterVariantTranches"
    assert (
        argv[argv.index("--snp-tranche") + 1] == "99.9"
        and argv[argv.index("--info-key") + 1] == "CNN_1D"
    )


def test_annotation_model_chain_prefixes_and_backend_checks(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    truth = _vcf(fake.folder / "truth.vcf")
    result = gp.extract_variant_annotations(
        str(vcf),
        ["QD", "MQ"],
        str(fake.folder / "ext"),
        resources=[{"path": str(truth), "name": "truth", "tags": {"training": True}}],
        modes=["SNP"],
        maximum_number_of_unlabeled_variants=100,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert "--resource:truth,training=true" in argv and argv.count("-A") == 2
    assert result["annotations_hdf5"].endswith("ext.annot.hdf5")
    assert any(item.endswith("ext.vcf.gz") for item in result["prefix_outputs"])
    trained = gp.train_variant_annotations_model(
        result["annotations_hdf5"], str(fake.folder / "model"), gatk_path=fake.path, verbose=False
    )
    assert fake.last()[fake.last().index("--model-backend") + 1] == "PYTHON_IFOREST"
    assert trained["status"] == "completed" and trained["prefix_outputs"]
    scored = gp.score_variant_annotations(
        str(vcf),
        ["QD", "MQ"],
        str(fake.folder / "model"),
        str(fake.folder / "scored"),
        snp_calibration_sensitivity_threshold=0.99,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[argv.index("--model-prefix") + 1] == str(fake.folder / "model")
    assert scored["scored_vcf"].endswith("scored.vcf.gz") and scored["summary"]["records"] == 1
    with pytest.raises(ConfigurationError, match="hyperparameters_json"):
        gp.train_variant_annotations_model(
            result["annotations_hdf5"],
            str(fake.folder / "m2"),
            model_backend="JAVA_BGMM",
            gatk_path=fake.path,
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="python_script"):
        gp.score_variant_annotations(
            str(vcf),
            ["QD"],
            "m",
            str(fake.folder / "s2"),
            model_backend="PYTHON_SCRIPT",
            gatk_path=fake.path,
            verbose=False,
        )


def test_prefix_outputs_are_never_overwritten(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    _file(fake.folder / "ext.annot.hdf5")
    with pytest.raises(ConfigurationError, match="overwrite"):
        gp.extract_variant_annotations(
            str(vcf), ["QD"], str(fake.folder / "ext"), gatk_path=fake.path, verbose=False
        )


def test_cnn_scoring_arguments(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    ref = _file(fake.folder / "ref.fa")
    gp.score_variants_with_cnn(
        str(vcf),
        str(ref),
        str(fake.folder / "o.vcf.gz"),
        batch_size=32,
        accelerator="cpu",
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "NVScoreVariants" and argv[argv.index("--tensor-type") + 1] == "reference"
    assert (
        argv[argv.index("--batch-size") + 1] == "32"
        and argv[argv.index("--accelerator") + 1] == "cpu"
    )
    with pytest.raises(ConfigurationError, match="reads"):
        gp.score_variants_with_cnn(
            str(vcf),
            str(ref),
            str(fake.folder / "o2.vcf"),
            tensor_type="read_tensor",
            gatk_path=fake.path,
            verbose=False,
        )


def test_panel_of_normals_uses_a_gendb_url(fake) -> None:
    db = fake.folder / "pon_db"
    db.mkdir()
    result = gp.create_somatic_panel_of_normals(
        str(db),
        str(fake.folder / "pon.vcf"),
        min_sample_count=3,
        call_genotypes=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "CreateSomaticPanelOfNormals"
    assert argv[argv.index("-V") + 1] == f"gendb://{db}"
    assert (
        argv[argv.index("--min-sample-count") + 1] == "3"
        and argv[argv.index("--call-genotypes") + 1] == "true"
    )
    assert result["status"] == "completed"
    with pytest.raises(ConfigurationError, match="existing folder"):
        gp.create_somatic_panel_of_normals(
            str(db / "nope"), str(fake.folder / "p2.vcf"), gatk_path=fake.path, verbose=False
        )


def test_mutect_and_mitochondrial_filters_arguments(fake) -> None:
    vcf = _vcf(fake.folder / "m2.vcf")
    ref = _file(fake.folder / "ref.fa")
    table = _file(fake.folder / "contam.table")
    priors = _file(fake.folder / "priors.tar.gz")
    result = gp.filter_mutect_calls(
        str(vcf),
        str(ref),
        str(fake.folder / "f.vcf"),
        filtering_stats_path=str(fake.folder / "f.stats"),
        contamination_tables=[str(table)],
        orientation_bias_priors=[str(priors)],
        threshold_strategy="FALSE_DISCOVERY_RATE",
        false_discovery_rate=0.1,
        mitochondria_mode=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "FilterMutectCalls"
    assert argv[argv.index("--threshold-strategy") + 1] == "FALSE_DISCOVERY_RATE"
    assert argv[argv.index("--false-discovery-rate") + 1] == "0.1"
    assert argv[argv.index("--mitochondria-mode") + 1] == "true" and "--microbial-mode" not in argv
    assert result["filtering_stats"].endswith("f.stats")
    gp.filter_mt_low_heteroplasmy(
        str(vcf),
        str(fake.folder / "lh.vcf"),
        low_het_threshold=0.2,
        max_allowed_low_hets=5,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert (
        argv[0] == "MTLowHeteroplasmyFilterTool"
        and argv[argv.index("--low-het-threshold") + 1] == "0.2"
    )
    gp.filter_numt_variants(
        str(vcf),
        str(fake.folder / "nm.vcf"),
        autosomal_coverage=120.0,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "NuMTFilterTool" and argv[argv.index("--autosomal-coverage") + 1] == "120.0"
    bam = _file(fake.folder / "r.bam")
    image = _file(fake.folder / "ref.img")
    gp.filter_alignment_artifacts(
        str(vcf),
        str(bam),
        str(ref),
        str(image),
        str(fake.folder / "faa.vcf"),
        kmer_size=25,
        dont_skip_filtered_variants=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "FilterAlignmentArtifacts" and argv[argv.index("--kmer-size") + 1] == "25"
    assert argv[argv.index("--bwa-mem-index-image") + 1] == str(image)


def test_variant_filtration_pairs_expressions_with_names(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    result = gp.filter_vcf_variants_gatk(
        str(vcf),
        str(fake.folder / "vf.vcf"),
        filter_expressions=["QD < 2.0", "FS > 60.0"],
        filter_names=["lowQD", "highFS"],
        genotype_filter_expressions=["DP < 5"],
        genotype_filter_names=["lowDP"],
        set_filtered_genotype_to_no_call=True,
        cluster_size=3,
        cluster_window_size=10,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv.count("--filter-expression") == 2 and argv.count("--filter-name") == 2
    assert argv[argv.index("--genotype-filter-name") + 1] == "lowDP"
    assert argv[argv.index("--cluster-size") + 1] == "3"
    assert result["status"] == "completed"
    with pytest.raises(ConfigurationError, match="exactly one filter name"):
        gp.filter_vcf_variants_gatk(
            str(vcf),
            str(fake.folder / "bad.vcf"),
            filter_expressions=["QD < 2.0"],
            gatk_path=fake.path,
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="must not start with"):
        gp.filter_vcf_variants_gatk(
            str(vcf),
            str(fake.folder / "bad2.vcf"),
            filter_expressions=["-QD"],
            filter_names=["x"],
            gatk_path=fake.path,
            verbose=False,
        )


def test_annotator_header_and_dictionary_tools_arguments(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    ref = _file(fake.folder / "ref.fa")
    bam = _file(fake.folder / "r.bam")
    gp.annotate_vcf_gatk(
        str(vcf),
        str(fake.folder / "an.vcf"),
        annotations=["Coverage"],
        annotation_groups=["StandardAnnotation"],
        reads=[str(bam)],
        reference=str(ref),
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "VariantAnnotator" and argv[argv.index("-A") + 1] == "Coverage"
    assert argv[argv.index("-G") + 1] == "StandardAnnotation" and argv[argv.index("-I") + 1] == str(
        bam
    )
    header = _vcf(fake.folder / "header.vcf")
    gp.fix_vcf_header_gatk(
        str(vcf),
        str(fake.folder / "fx.vcf"),
        header=str(header),
        enforce_same_samples=False,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "FixVcfHeader" and argv[argv.index("--ENFORCE_SAME_SAMPLES") + 1] == "false"
    assert argv[argv.index("--HEADER") + 1] == str(header)
    gp.fix_vcf_header_gatk(
        str(vcf), str(fake.folder / "fx2.vcf"), gatk_path=fake.path, verbose=False
    )
    assert "--ENFORCE_SAME_SAMPLES" not in fake.last()
    gp.update_vcf_sequence_dictionary_gatk(
        str(vcf),
        str(fake.folder / "ud.vcf"),
        reference=str(ref),
        replace=True,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "UpdateVCFSequenceDictionary" and argv[argv.index("--replace") + 1] == "true"
    with pytest.raises(ConfigurationError, match="dictionary source"):
        gp.update_vcf_sequence_dictionary_gatk(
            str(vcf), str(fake.folder / "ud2.vcf"), gatk_path=fake.path, verbose=False
        )


def test_liftover_split_convert_interval_list_and_sites_only(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    chain = _file(fake.folder / "a.chain")
    ref = _file(fake.folder / "new.fa")
    result = gp.liftover_vcf_gatk(
        str(vcf),
        str(chain),
        str(ref),
        str(fake.folder / "lo.vcf"),
        str(fake.folder / "rej.vcf"),
        min_match=0.95,
        recover_swapped_ref_alt=True,
        tags_to_reverse=["AF"],
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "LiftoverVcf" and argv[argv.index("--CHAIN") + 1] == str(chain)
    assert argv[argv.index("--REFERENCE_SEQUENCE") + 1] == str(ref)
    assert argv[argv.index("--LIFTOVER_MIN_MATCH") + 1] == "0.95"
    assert result["reject"].endswith("rej.vcf") and result["reject_summary"]["records"] == 1
    with pytest.raises(ConfigurationError, match="must differ"):
        gp.liftover_vcf_gatk(
            str(vcf),
            str(chain),
            str(ref),
            str(fake.folder / "same.vcf"),
            str(fake.folder / "same.vcf"),
            gatk_path=fake.path,
            verbose=False,
        )
    split = gp.split_vcf_snps_indels_gatk(
        str(vcf),
        str(fake.folder / "s.vcf"),
        str(fake.folder / "i.vcf"),
        strict=False,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "SplitVcfs" and argv[argv.index("--STRICT") + 1] == "false"
    assert split["snp_summary"]["records"] == 1 and split["indel_summary"]["records"] == 1
    gp.convert_vcf_format_gatk(
        str(vcf),
        str(fake.folder / "c.bcf"),
        require_index=False,
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "VcfFormatConverter" and argv[argv.index("--REQUIRE_INDEX") + 1] == "false"
    with pytest.raises(ConfigurationError, match="must end with"):
        gp.convert_vcf_format_gatk(
            str(vcf), str(fake.folder / "c.txt"), gatk_path=fake.path, verbose=False
        )
    ivl = gp.vcf_to_interval_list_gatk(
        str(vcf),
        str(fake.folder / "v.interval_list"),
        include_filtered=True,
        variant_id_method="FIRST",
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert (
        argv[argv.index("--VARIANT_ID_METHOD") + 1] == "FIRST"
        and argv[argv.index("--INCLUDE_FILTERED") + 1] == "true"
    )
    assert ivl["intervals"] == [["chr1", 10, 10, "+", "interval-1"]]
    with pytest.raises(ConfigurationError, match="variant_id_method"):
        gp.vcf_to_interval_list_gatk(
            str(vcf),
            str(fake.folder / "v2.interval_list"),
            variant_id_method="BAD",
            gatk_path=fake.path,
            verbose=False,
        )
    gp.make_sites_only_vcf_gatk(
        str(vcf), str(fake.folder / "so.vcf"), samples=["S1"], gatk_path=fake.path, verbose=False
    )
    argv = fake.last()
    assert argv[0] == "MakeSitesOnlyVcf" and argv[argv.index("--SAMPLE") + 1] == "S1"


def test_print_variants_spark_arguments(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    gp.print_variants_spark(
        str(vcf),
        str(fake.folder / "sp.vcf"),
        variant_shard_size=500,
        spark_master="local[2]",
        gatk_path=fake.path,
        verbose=False,
    )
    argv = fake.last()
    assert argv[0] == "PrintVariantsSpark" and argv[argv.index("--variantShardSize") + 1] == "500"
    assert argv[argv.index("--") + 1 :] == ["--spark-runner", "LOCAL", "--spark-master", "local[2]"]


def test_existing_output_and_bad_values_are_refused(fake) -> None:
    vcf = _vcf(fake.folder / "in.vcf")
    existing = _file(fake.folder / "exists.vcf")
    with pytest.raises(ConfigurationError, match="must not exist"):
        gp.make_sites_only_vcf_gatk(str(vcf), str(existing), gatk_path=fake.path, verbose=False)
    with pytest.raises(ConfigurationError, match="mode"):
        gp.apply_vqsr(
            str(vcf),
            str(vcf),
            str(fake.folder / "o.vcf"),
            mode="BAD",
            gatk_path=fake.path,
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="at most"):
        gp.apply_vqsr(
            str(vcf),
            str(vcf),
            str(fake.folder / "o2.vcf"),
            truth_sensitivity_filter_level=101.0,
            gatk_path=fake.path,
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="True or False"):
        gp.apply_vqsr(
            str(vcf),
            str(vcf),
            str(fake.folder / "o3.vcf"),
            exclude_filtered="yes",
            gatk_path=fake.path,
            verbose=False,
        )  # type: ignore[arg-type]
    assert not (fake.folder / "o.vcf").exists()


def test_vcf_summary_counts_plain_and_gzipped(tmp_path: Path) -> None:
    plain = _vcf(tmp_path / "a.vcf", records=4)
    zipped = tmp_path / "a.vcf.gz"
    with gzip.open(zipped, "wt") as handle:
        handle.write(plain.read_text())
    for path in (plain, zipped):
        summary = gp._vcf_summary(path)
        assert summary == {"records": 4, "filters": {".": 4}, "samples": ["S1"]}
    assert gp._vcf_summary(tmp_path / "x.bcf") is None


def test_gatk_report_parser() -> None:
    text = (
        "#:GATKReport.v1.1:2\n#:GATKTable:2:2:%s:%.3f:;\n#:GATKTable:A:first table\nname value\nx 1\ny 2.5\n\n"
        "#:GATKTable:1:1:%s:;\n#:GATKTable:B:second\nk\nz\n\n"
    )
    tables = gp._parse_gatk_report(text)
    assert tables["A"] == {
        "description": "first table",
        "columns": ["name", "value"],
        "rows": [["x", 1], ["y", 2.5]],
    }
    assert tables["B"]["rows"] == [["z"]]


# ---------------------------------------------------------------------------
# Real GATK runs (skipped when GATK is not installed)
# ---------------------------------------------------------------------------


def _real_launcher() -> str | None:
    explicit = os.environ.get("DNAKIT_GATK")
    if explicit and Path(explicit).is_file():
        return explicit
    return shutil.which("gatk")


@pytest.fixture(scope="module")
def real(tmp_path_factory: pytest.TempPathFactory):
    launcher = _real_launcher()
    if launcher is None or shutil.which("java") is None:
        pytest.skip("GATK 4.x and Java are not installed.")
    folder = tmp_path_factory.mktemp("gatk_real")
    rng = random.Random(3)
    bases = "ACGT"
    fasta = folder / "ref.fa"
    lengths = {"chr1": 3000, "chr2": 2000}
    sequences = {
        name: "".join(rng.choice(bases) for _ in range(size)) for name, size in lengths.items()
    }
    _write_reference(fasta, sequences)
    subprocess.run(
        [launcher, "CreateSequenceDictionary", "-R", str(fasta)],
        check=True,
        capture_output=True,
        timeout=300,
    )
    rows = []
    for index in range(40):
        chrom = "chr1" if index < 25 else "chr2"
        position = 50 + 60 * (index if index < 25 else index - 25)
        ref = sequences[chrom][position - 1]
        if index % 10 == 9:  # a deletion
            ref = sequences[chrom][position - 1 : position + 2]
            alt = ref[0]
        else:
            alt = next(b for b in bases if b != ref)
        good = index % 3 != 0
        rows.append((chrom, position, ref, alt, 20.0 if good else 1.5, 60.0 if good else 30.0))
    call = folder / "call.vcf"
    with call.open("w") as handle:
        handle.write(VCF_HEADER.replace("chr1,length=3000", "chr1,length=3000"))
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for chrom, position, ref, alt, qd, mq in rows:
            handle.write(
                f"{chrom}\t{position}\trs{position}\t{ref}\t{alt}\t50\t.\tQD={qd};MQ={mq};FS=1.0;SOR=0.5;DP=30\tGT\t0/1\n"
            )

    class Real:
        gatk = launcher
        dir = folder
        reference = fasta
        vcf = call

    return Real


def test_real_validate_and_variant_filtration(real, tmp_path: Path) -> None:
    ok = gp.validate_vcf_gatk(
        str(real.vcf), reference=str(real.reference), gatk_path=real.gatk, verbose=False
    )
    assert ok["status"] == "completed" and ok["valid"] is True and ok["version"]
    bad_vcf = tmp_path / "bad.vcf"
    lines = real.vcf.read_text().splitlines()
    first = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    fields = lines[first].split("\t")
    fields[3] = "A" if fields[3] != "A" else "C"
    lines[first] = "\t".join(fields)
    bad_vcf.write_text("\n".join(lines) + "\n")
    bad = gp.validate_vcf_gatk(
        str(bad_vcf), reference=str(real.reference), gatk_path=real.gatk, verbose=False
    )
    assert bad["valid"] is False and bad["status"] == "failed"
    out = tmp_path / "vf.vcf"
    filtered = gp.filter_vcf_variants_gatk(
        str(real.vcf),
        str(out),
        filter_expressions=["QD < 10.0"],
        filter_names=["lowQD"],
        gatk_path=real.gatk,
        verbose=False,
    )
    assert filtered["status"] == "completed"
    assert filtered["summary"]["records"] == 40
    assert (
        filtered["summary"]["filters"]["lowQD"] == 14
        and filtered["summary"]["filters"]["PASS"] == 26
    )
    assert filtered["index_files"]


def test_real_vcf_utilities(real, tmp_path: Path) -> None:
    split = gp.split_vcf_snps_indels_gatk(
        str(real.vcf),
        str(tmp_path / "snp.vcf"),
        str(tmp_path / "indel.vcf"),
        gatk_path=real.gatk,
        verbose=False,
    )
    assert (split["snp_summary"]["records"], split["indel_summary"]["records"]) == (36, 4)
    sites = gp.make_sites_only_vcf_gatk(
        str(real.vcf), str(tmp_path / "sites.vcf"), gatk_path=real.gatk, verbose=False
    )
    assert sites["summary"]["records"] == 40 and sites["summary"]["samples"] == []
    to_bcf = gp.convert_vcf_format_gatk(
        str(real.vcf),
        str(tmp_path / "c.bcf"),
        require_index=False,
        gatk_path=real.gatk,
        verbose=False,
    )
    assert to_bcf["status"] == "completed" and to_bcf["summary"] is None
    back = gp.convert_vcf_format_gatk(
        str(tmp_path / "c.bcf"),
        str(tmp_path / "back.vcf"),
        require_index=False,
        gatk_path=real.gatk,
        verbose=False,
    )
    assert back["summary"]["records"] == 40
    ivl = gp.vcf_to_interval_list_gatk(
        str(real.vcf), str(tmp_path / "v.interval_list"), gatk_path=real.gatk, verbose=False
    )
    assert len(ivl["intervals"]) == 40 and ivl["intervals"][0][:3] == ["chr1", 50, 50]
    fixed = gp.fix_vcf_header_gatk(
        str(real.vcf), str(tmp_path / "fixed.vcf"), gatk_path=real.gatk, verbose=False
    )
    assert fixed["summary"]["records"] == 40
    updated = gp.update_vcf_sequence_dictionary_gatk(
        str(real.vcf),
        str(tmp_path / "upd.vcf"),
        reference=str(real.reference),
        replace=True,
        gatk_path=real.gatk,
        verbose=False,
    )
    assert updated["summary"]["records"] == 40


def test_real_evaluate_variants_counts_snps_and_indels(real, tmp_path: Path) -> None:
    report = gp.evaluate_variants_gatk(
        [str(real.vcf)],
        str(real.reference),
        str(tmp_path / "eval.report"),
        gatk_path=real.gatk,
        verbose=False,
    )
    assert report["status"] == "completed"
    table = report["tables"]["CountVariants"]
    row = table["rows"][0]
    named = dict(zip(table["columns"], row, strict=False))
    assert named["nSNPs"] == 36 and named["nDeletions"] == 4


@pytest.fixture(scope="module")
def real_scored(real, tmp_path_factory: pytest.TempPathFactory):
    """A larger annotated SNP callset (900 sites, a quarter of them poor) with a truth VCF, both indexed."""
    folder = tmp_path_factory.mktemp("gatk_scored")
    rng = random.Random(11)
    fasta = folder / "big.fa"
    sequence = "".join(rng.choice("ACGT") for _ in range(70_000))
    _write_reference(fasta, {"chr1": sequence})
    subprocess.run(
        [real.gatk, "CreateSequenceDictionary", "-R", str(fasta)],
        check=True,
        capture_output=True,
        timeout=300,
    )
    header = VCF_HEADER.replace(
        "##contig=<ID=chr1,length=3000>\n##contig=<ID=chr2,length=2000>\n",
        "##contig=<ID=chr1,length=70000>\n",
    )
    call, truth = folder / "call.vcf", folder / "truth.vcf"
    with call.open("w") as c_handle, truth.open("w") as t_handle:
        c_handle.write(header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        t_handle.write(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=70000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        for index in range(900):
            position = 100 + 75 * index
            ref = sequence[position - 1]
            alt = next(b for b in "ACGT" if b != ref)
            good = index % 4 != 0
            qd = rng.gauss(18, 3) if good else rng.gauss(4, 2)
            fs = abs(rng.gauss(2, 2)) if good else abs(rng.gauss(30, 10))
            mq = rng.gauss(60, 1) if good else rng.gauss(40, 5)
            sor = abs(rng.gauss(1, 0.5)) if good else abs(rng.gauss(4, 1))
            c_handle.write(
                f"chr1\t{position}\t.\t{ref}\t{alt}\t50\t.\tQD={qd:.3f};FS={fs:.3f};MQ={mq:.3f};SOR={sor:.3f};DP=40\tGT\t0/1\n"
            )
            if good and index % 3 == 0:
                t_handle.write(f"chr1\t{position}\trs{position}\t{ref}\t{alt}\t.\t.\t.\n")
    for path in (call, truth):
        subprocess.run(
            [real.gatk, "IndexFeatureFile", "-I", str(path)],
            check=True,
            capture_output=True,
            timeout=300,
        )

    class Scored:
        gatk = real.gatk
        dir = folder
        reference = fasta
        vcf = call
        truth_vcf = truth

    return Scored


def test_real_vqsr_chain_and_tranche_filtering(real_scored, tmp_path: Path) -> None:
    truth = {
        "path": str(real_scored.truth_vcf),
        "name": "truth",
        "tags": {"known": False, "training": True, "truth": True, "prior": 15.0},
    }
    recal = gp.recalibrate_variants_vqsr(
        [str(real_scored.vcf)],
        [truth],
        ["QD", "FS", "MQ", "SOR"],
        str(tmp_path / "s.recal"),
        str(tmp_path / "s.tranches"),
        truth_sensitivity_tranches=[100.0, 99.0, 90.0],
        max_gaussians=2,
        max_negative_gaussians=2,
        minimum_bad_variants=100,
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    assert recal["status"] == "completed", recal["log_tail"][-500:]
    assert recal["index_files"] and Path(recal["tranches_file"]).read_text().count("\n") > 3
    applied = gp.apply_vqsr(
        str(real_scored.vcf),
        recal["recal_file"],
        str(tmp_path / "v.vcf"),
        tranches_file=recal["tranches_file"],
        truth_sensitivity_filter_level=99.0,
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    assert applied["status"] == "completed" and applied["summary"]["records"] == 900
    assert "PASS" in applied["summary"]["filters"] and any(
        name.startswith("VQSRTranche") for name in applied["summary"]["filters"]
    )
    tranched = gp.filter_variant_tranches(
        applied["output"],
        [str(real_scored.truth_vcf)],
        str(tmp_path / "t.vcf"),
        snp_tranches=[99.9],
        info_key="VQSLOD",
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    assert tranched["status"] == "completed" and tranched["summary"]["records"] == 900


def test_real_annotation_model_chain_extract_train_score(
    real_scored, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("sklearn")
    pytest.importorskip("h5py")
    pytest.importorskip("dill")
    monkeypatch.setenv(
        "PATH", f"{Path(sys.executable).parent}{os.pathsep}{os.environ.get('PATH', '')}"
    )
    labels = [
        {
            "path": str(real_scored.truth_vcf),
            "name": "truth",
            "tags": {"training": True, "calibration": True},
        }
    ]
    names = ["QD", "FS", "MQ", "SOR"]
    extracted = gp.extract_variant_annotations(
        str(real_scored.vcf),
        names,
        str(tmp_path / "ext"),
        resources=labels,
        modes=["SNP"],
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    assert extracted["status"] == "completed" and extracted["annotations_hdf5"]
    trained = gp.train_variant_annotations_model(
        extracted["annotations_hdf5"],
        str(tmp_path / "trn"),
        modes=["SNP"],
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    if trained["status"] != "completed":
        pytest.skip(
            "the GATK Python environment (scikit-learn IsolationForest) is not usable here: "
            + trained["log_tail"][-300:]
        )
    assert any(name.endswith(".snp.scorer.pkl") for name in trained["prefix_outputs"])
    scored = gp.score_variant_annotations(
        str(real_scored.vcf),
        names,
        str(tmp_path / "trn"),
        str(tmp_path / "sc"),
        resources=labels,
        modes=["SNP"],
        snp_calibration_sensitivity_threshold=0.95,
        gatk_path=real_scored.gatk,
        verbose=False,
    )
    assert scored["status"] == "completed" and scored["summary"]["records"] == 900
    assert "LOW_SCORE" in scored["summary"]["filters"]


def test_real_liftover_with_shifted_reference(real, tmp_path: Path) -> None:
    sequences: dict[str, str] = {}
    name = ""
    for line in real.reference.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = ""
        else:
            sequences[name] += line
    shifted = {"chr1": sequences["chr1"][5:], "chr2": sequences["chr2"]}
    new_ref = tmp_path / "new.fa"
    _write_reference(new_ref, shifted)
    subprocess.run(
        [real.gatk, "CreateSequenceDictionary", "-R", str(new_ref)],
        check=True,
        capture_output=True,
        timeout=300,
    )
    chain = tmp_path / "a.chain"
    chain.write_text(
        "chain 1000 chr1 3000 + 5 3000 chr1 2995 + 0 2995 1\n2995\n\n"
        "chain 900 chr2 2000 + 0 2000 chr2 2000 + 0 2000 2\n2000\n\n"
    )
    result = gp.liftover_vcf_gatk(
        str(real.vcf),
        str(chain),
        str(new_ref),
        str(tmp_path / "lift.vcf"),
        str(tmp_path / "reject.vcf"),
        gatk_path=real.gatk,
        verbose=False,
    )
    assert result["status"] == "completed"

    def first_record(path: Path) -> list[str]:
        text = path.read_text().splitlines()
        return next(row for row in text if not row.startswith("#")).split("\t")

    first_in = first_record(real.vcf)
    first_out = first_record(tmp_path / "lift.vcf")
    assert int(first_out[1]) == int(first_in[1]) - 5
    assert result["summary"]["records"] + result["reject_summary"]["records"] == 40


def test_real_sites_only_spark_and_missing_python_environment_is_reported(
    real, tmp_path: Path
) -> None:
    sites = gp.make_sites_only_vcf_gatk(
        str(real.vcf), str(tmp_path / "sites.vcf"), gatk_path=real.gatk, verbose=False
    )
    spark = gp.print_variants_spark(
        sites["output"],
        str(tmp_path / "spark.vcf"),
        spark_master="local[1]",
        gatk_path=real.gatk,
        timeout_seconds=600,
        verbose=False,
    )
    assert spark["status"] == "completed" and spark["summary"]["records"] == 40
    cnn = gp.score_variants_with_cnn(
        str(real.vcf),
        str(real.reference),
        str(tmp_path / "cnn.vcf"),
        gatk_path=real.gatk,
        timeout_seconds=300,
        verbose=False,
    )
    assert cnn["status"] in {"completed", "failed"}
    if cnn["status"] == "failed":
        assert cnn["return_code"] != 0 and cnn["summary"] is None


def test_real_cnn_scoring_with_the_gatk_python_archive(
    real_scored, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = Path(real_scored.gatk).resolve().parent / "gatkPythonPackageArchive.zip"
    if not archive.is_file():
        pytest.skip("gatkPythonPackageArchive.zip is not next to the gatk launcher.")
    for module in ("torch", "pytorch_lightning", "Bio", "pysam"):
        pytest.importorskip(module)
    import zipfile

    packages = tmp_path / "gatk_python"
    with zipfile.ZipFile(archive) as handle:
        handle.extractall(packages)
    monkeypatch.setenv("PYTHONPATH", f"{packages}{os.pathsep}{os.environ.get('PYTHONPATH', '')}")
    monkeypatch.setenv("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "1")
    monkeypatch.setenv("PATH", f"{Path(sys.executable).parent}{os.pathsep}{os.environ.get('PATH', '')}")
    lines = real_scored.vcf.read_text().splitlines()
    header = [line for line in lines if line.startswith("##")]
    header += [
        '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="MQRankSum">',
        '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="ReadPosRankSum">',
    ]
    records = [line for line in lines if not line.startswith("#")][:30]
    source = tmp_path / "cnn_in.vcf"
    source.write_text(
        "\n".join([*header, next(line for line in lines if line.startswith("#CHROM")), *records]) + "\n"
    )
    result = gp.score_variants_with_cnn(
        str(source),
        str(real_scored.reference),
        str(tmp_path / "cnn.vcf"),
        accelerator="cpu",
        gatk_path=real_scored.gatk,
        timeout_seconds=900,
        verbose=False,
    )
    assert result["status"] == "completed", result["log_tail"][-500:]
    assert result["summary"]["records"] == 30
    scored = (tmp_path / "cnn.vcf").read_text()
    assert "CNN_1D=" in scored
