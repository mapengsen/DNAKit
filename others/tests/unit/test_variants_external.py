"""Contract tests for FreeBayes and VCFtools facades."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools import default_tool_registry
from dnakit.variants import (
    annotate_vcf_ensembl_vep,
    annotate_vcf_vcfanno,
    call_cutesv,
    call_deepvariant,
    call_delly,
    call_freebayes,
    call_gatk_haplotypecaller,
    call_sniffles2,
    filter_vcf_vcftools,
    load_gemini_database,
    merge_structural_variants_survivor,
    query_gemini_database,
)


def _fake_freebayes(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "print('##fileformat=VCFv4.2')\n"
        "print('#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO')\n"
        "print('chr1\\t1\\t.\\tA\\tG\\t30\\tPASS\\t.')\n"
        "print('freebayes diagnostic', file=__import__('sys').stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_vcftools(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "prefix = Path(sys.argv[sys.argv.index('--out') + 1])\n"
        "prefix.with_suffix('.recode.vcf').write_text("
        "'##fileformat=VCFv4.2\\n', encoding='ascii')\n"
        "print('vcftools diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_vcfanno(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "print('##fileformat=VCFv4.2')\n"
        "print('#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO')\n"
        "print('chr1\\t1\\t.\\tA\\tG\\t30\\tPASS\\tANNOT=1')\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_survivor(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[-1]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_gatk(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "target = Path(sys.argv[sys.argv.index('-O') + 1])\n"
        "target.write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n"
        "print('gatk diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_deepvariant(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "for argument in sys.argv:\n"
        "    if argument.startswith('--output_vcf='):\n"
        "        Path(argument.split('=', 1)[1]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n"
        "    if argument.startswith('--output_gvcf='):\n"
        "        Path(argument.split('=', 1)[1]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n"
        "print('deepvariant diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_cutesv(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[3]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_sniffles(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[sys.argv.index('--vcf') + 1]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_delly(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "print('##fileformat=VCFv4.2')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_gemini(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "args = sys.argv[1:]\n"
        "if args[0] == 'load':\n"
        "    Path(args[-1]).write_text('GEMINI DB\\n', encoding='ascii')\n"
        "else:\n"
        "    if '--header' in args:\n"
        "        print('chrom\\tstart\\tgene')\n"
        "    print('chr1\\t0\\tTEST')\n"
        "print('gemini diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_vep(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "output = Path(sys.argv[sys.argv.index('--output_file') + 1])\n"
        "output.write_text('##fileformat=VCFv4.2\\n##INFO=<ID=CSQ>\\n', encoding='ascii')\n"
        "print('vep diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_call_freebayes_writes_vcf_and_keeps_inputs(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    alignment = tmp_path / "sample.bam"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    executable = _fake_freebayes(tmp_path / "fake-freebayes")

    result = call_freebayes(
        reference,
        alignment,
        tmp_path / "calls",
        ploidy=2,
        min_alternate_fraction=0.1,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.vcf" in result["result"]["artifacts"]
    assert "##fileformat" in Path(result["result"]["artifacts"]["variants.vcf"]).read_text()


def test_gemini_load_and_query_facades_capture_outputs(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    executable = _fake_gemini(tmp_path / "fake-gemini")

    loaded = load_gemini_database(
        source,
        tmp_path / "gemini-load",
        annotation_type="VEP",
        cores=2,
        program=str(executable),
        verbose=False,
    )
    assert "variants.db" in loaded["result"]["artifacts"]
    database = Path(loaded["result"]["artifacts"]["variants.db"])

    queried = query_gemini_database(
        database,
        "select chrom, start, gene from variants",
        tmp_path / "gemini-query",
        result_format="sampledetail",
        program=str(executable),
        verbose=False,
    )
    assert "query.tsv" in queried["result"]["artifacts"]
    assert "TEST" in Path(queried["result"]["artifacts"]["query.tsv"]).read_text()


def test_ensembl_vep_facade_writes_annotated_vcf(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    cache = tmp_path / "vep-cache"
    cache.mkdir()
    executable = _fake_vep(tmp_path / "fake-vep")

    result = annotate_vcf_ensembl_vep(
        source,
        tmp_path / "vep-output",
        species="homo_sapiens",
        assembly="GRCh38",
        cache_directory=cache,
        offline=True,
        fork=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "annotated.vcf" in result["result"]["artifacts"]
    assert "CSQ" in Path(result["result"]["artifacts"]["annotated.vcf"]).read_text()


def test_filter_vcf_vcftools_writes_recoded_vcf(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    executable = _fake_vcftools(tmp_path / "fake-vcftools")

    result = filter_vcf_vcftools(
        source,
        tmp_path / "filtered",
        remove_indels=True,
        min_allele_frequency=0.05,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "filtered.recode.vcf" in result["result"]["artifacts"]


def test_annotate_vcf_vcfanno_captures_stdout(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    configuration = tmp_path / "annotations.toml"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    configuration.write_text("[[annotation]]\n", encoding="ascii")
    executable = _fake_vcfanno(tmp_path / "fake-vcfanno")

    result = annotate_vcf_vcfanno(
        source,
        configuration,
        tmp_path / "annotated",
        threads=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "annotated.vcf" in result["result"]["artifacts"]


def test_merge_structural_variants_survivor_writes_vcf(tmp_path: Path) -> None:
    first = tmp_path / "one.vcf"
    second = tmp_path / "two.vcf"
    first.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    second.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    executable = _fake_survivor(tmp_path / "fake-survivor")

    result = merge_structural_variants_survivor(
        (first, second),
        tmp_path / "merged",
        min_support=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "merged.vcf" in result["result"]["artifacts"]


def test_call_gatk_haplotypecaller_writes_vcf(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    alignment = tmp_path / "sample.bam"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    executable = _fake_gatk(tmp_path / "fake-gatk")

    result = call_gatk_haplotypecaller(
        reference,
        alignment,
        tmp_path / "gatk-calls",
        emit_gvcf=True,
        pair_hmm_threads=2,
        java_max_heap_gb=4,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.g.vcf.gz" in result["result"]["artifacts"]


def test_call_deepvariant_writes_vcf_and_gvcf(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    alignment = tmp_path / "sample.bam"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    executable = _fake_deepvariant(tmp_path / "fake-deepvariant")

    result = call_deepvariant(
        reference,
        alignment,
        tmp_path / "deepvariant-calls",
        model_type="WGS",
        num_shards=2,
        emit_gvcf=True,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "output.vcf.gz" in result["result"]["artifacts"]
    assert "output.g.vcf.gz" in result["result"]["artifacts"]


def test_call_cutesv_writes_vcf(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    alignment = tmp_path / "sample.bam"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    executable = _fake_cutesv(tmp_path / "fake-cutesv")

    result = call_cutesv(
        alignment,
        reference,
        tmp_path / "cutesv-calls",
        threads=2,
        sample="sample1",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.vcf" in result["result"]["artifacts"]


def test_call_sniffles2_writes_vcf(tmp_path: Path) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM")
    executable = _fake_sniffles(tmp_path / "fake-sniffles")

    result = call_sniffles2(
        alignment,
        tmp_path / "sniffles-calls",
        threads=2,
        min_sv_length=30,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.vcf" in result["result"]["artifacts"]


def test_call_delly_captures_vcf_stdout(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    alignment = tmp_path / "sample.bam"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    executable = _fake_delly(tmp_path / "fake-delly")

    result = call_delly(
        alignment,
        reference,
        tmp_path / "delly-calls",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.vcf" in result["result"]["artifacts"]


def test_variant_facades_validate_conflicting_filters(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    with pytest.raises(ConfigurationError, match="cannot both"):
        filter_vcf_vcftools(
            source,
            tmp_path / "filtered",
            remove_indels=True,
            remove_snps=True,
            verbose=False,
        )


def test_variant_facades_are_registered_as_external_tools() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.variants.annotate_vcf_ensembl_vep",
        "dnakit.variants.annotate_vcf_vcfanno",
        "dnakit.variants.call_freebayes",
        "dnakit.variants.call_gatk_haplotypecaller",
        "dnakit.variants.call_deepvariant",
        "dnakit.variants.call_cutesv",
        "dnakit.variants.call_delly",
        "dnakit.variants.call_sniffles2",
        "dnakit.variants.filter_vcf_vcftools",
        "dnakit.variants.merge_structural_variants_survivor",
        "dnakit.variants.load_gemini_database",
        "dnakit.variants.query_gemini_database",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True
