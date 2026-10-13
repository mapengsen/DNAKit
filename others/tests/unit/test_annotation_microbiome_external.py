"""Contract tests for DIAMOND, MMseqs2, annotation and Kraken2 facades."""

from __future__ import annotations

import sys
from pathlib import Path

from dnakit.alignment import (
    cluster_mmseqs_protein,
    search_diamond_blastp,
    search_diamond_blastx,
    search_mmseqs_protein,
)
from dnakit.annotation import annotate_bakta, annotate_prokka, predict_checkm2
from dnakit.microbiome import classify_kraken2
from dnakit.microbiome import (
    classify_centrifuge,
    classify_gtdbtk,
    classify_kaiju,
    estimate_bracken_abundance,
    profile_humann,
)
from dnakit.epigenome import extract_methylation_methyldackel
from dnakit.assembly import assemble_unicycler
from dnakit.tools import default_tool_registry


def _fake_program(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import os\n"
        "import sys\n"
        "args = sys.argv[1:]\n"
        "name = Path(sys.argv[0]).name\n"
        "def value(flag):\n"
        "    return Path(args[args.index(flag) + 1])\n"
        "if name == 'fake-diamond':\n"
        "    if args[0] == 'makedb':\n"
        "        Path(str(value('--db')) + '.dmnd').write_text('db')\n"
        "    else:\n"
        "        value('--out').write_text('query\\thit\\n')\n"
        "elif name == 'fake-mmseqs':\n"
        "    if args[0] == 'easy-search':\n"
        "        Path(args[3]).write_text('query\\thit\\n')\n"
        "    else:\n"
        "        prefix = Path(args[2])\n"
        "        for suffix in ('_cluster.tsv', '_rep_seq.fasta', '_all_seqs.fasta'):\n"
        "            Path(str(prefix) + suffix).write_text('cluster\\n')\n"
        "elif name == 'fake-prokka':\n"
        "    out = value('--outdir'); out.mkdir()\n"
        "    prefix = args[args.index('--prefix') + 1]\n"
        "    (out / f'{prefix}.gff').write_text('##gff-version 3\\n')\n"
        "    (out / f'{prefix}.faa').write_text('>protein\\nM\\n')\n"
        "elif name == 'fake-bakta':\n"
        "    out = value('--output'); out.mkdir()\n"
        "    prefix = args[args.index('--prefix') + 1]\n"
        "    for suffix in ('.gff3', '.tsv', '.json', '.fna'):\n"
        "        (out / f'{prefix}{suffix}').write_text('annotation\\n')\n"
        "elif name == 'fake-checkm2':\n"
        "    out = value('--output-directory'); out.mkdir()\n"
        "    (out / 'quality_report.tsv').write_text('Name\\tCompleteness\\n')\n"
        "elif name == 'fake-kraken2':\n"
        "    value('--report').write_text('0.0\\t0\\t0\\tU\\t0\\troot\\n')\n"
        "    value('--output').write_text('read\\tU\\t0\\t0\\n')\n"
        "elif name == 'fake-bracken':\n"
        "    value('-o').write_text('name\\t100\\n')\n"
        "elif name == 'fake-centrifuge':\n"
        "    value('-S').write_text('read\\t0\\n')\n"
        "    value('--report-file').write_text('name\\t100\\n')\n"
        "elif name == 'fake-kaiju':\n"
        "    value('-o').write_text('C\\tread\\t1\\n')\n"
        "elif name == 'fake-unicycler':\n"
        "    out = value('-o'); out.mkdir()\n"
        "    (out / 'assembly.fasta').write_text('>contig\\nATGC\\n')\n"
        "elif name == 'fake-gtdbtk':\n"
        "    assert Path(os.environ['GTDBTK_DATA_PATH']).is_dir()\n"
        "    out = value('--out_dir'); out.mkdir()\n"
        "    (out / 'gtdbtk.bac120.summary.tsv').write_text('user_genome\\ttaxonomy\\n')\n"
        "elif name == 'fake-humann':\n"
        "    out = value('--output'); out.mkdir(exist_ok=True)\n"
        "    stem = value('--input').name.split('.', 1)[0]\n"
        "    for suffix in ('_genefamilies.tsv', '_pathabundance.tsv', '_pathcoverage.tsv'):\n"
        "        (out / f'{stem}{suffix}').write_text('feature\\tvalue\\n')\n"
        "elif name == 'fake-MethylDackel':\n"
        "    prefix = Path(args[-1])\n"
        "    (Path(str(prefix) + '_CpG.bedGraph')).write_text('chr1\\t0\\t1\\n')\n"
        "    if '--CHG' in args:\n"
        "        (Path(str(prefix) + '_CHG.bedGraph')).write_text('chr1\\t0\\t1\\n')\n"
        "    if '--CHH' in args:\n"
        "        (Path(str(prefix) + '_CHH.bedGraph')).write_text('chr1\\t0\\t1\\n')\n"
        "print('fake diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def test_diamond_and_mmseqs_facades_capture_outputs(tmp_path: Path) -> None:
    query = tmp_path / "query.fa"
    reference = tmp_path / "reference.fa"
    query.write_text(">q\nM\n", encoding="ascii")
    reference.write_text(">r\nM\n", encoding="ascii")
    diamond = _fake_program(tmp_path / "fake-diamond")
    mmseqs = _fake_program(tmp_path / "fake-mmseqs")

    for function, output in (
        (search_diamond_blastp, "diamond-p"),
        (search_diamond_blastx, "diamond-x"),
    ):
        result = function(query, reference, tmp_path / output, program=str(diamond), verbose=False)
        assert "reference.dmnd" in result["result"]["artifacts"]
        assert "hits.tsv" in result["result"]["artifacts"]
    result = search_mmseqs_protein(query, reference, tmp_path / "mmseqs-search", program=str(mmseqs), verbose=False)
    assert "hits.tsv" in result["result"]["artifacts"]
    result = cluster_mmseqs_protein(reference, tmp_path / "mmseqs-cluster", program=str(mmseqs), verbose=False)
    assert "clusters_cluster.tsv" in result["result"]["artifacts"]


def test_annotation_and_kraken_facades_capture_outputs(tmp_path: Path) -> None:
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">contig\nATGC\n", encoding="ascii")
    database = tmp_path / "database"
    database.mkdir()
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "one.fna").write_text(">one\nATGC\n", encoding="ascii")
    checkm_database = tmp_path / "checkm.dmnd"
    checkm_database.write_text("db", encoding="ascii")
    kraken_database = tmp_path / "kraken-db"
    kraken_database.mkdir()
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nATGC\n+\nIIII\n", encoding="ascii")

    prokka = _fake_program(tmp_path / "fake-prokka")
    bakta = _fake_program(tmp_path / "fake-bakta")
    checkm2 = _fake_program(tmp_path / "fake-checkm2")
    kraken2 = _fake_program(tmp_path / "fake-kraken2")

    assert "annotation/annotation.gff" in annotate_prokka(
        assembly, tmp_path / "prokka", program=str(prokka), verbose=False
    )["result"]["artifacts"]
    assert "bakta/annotation.gff3" in annotate_bakta(
        assembly, database, tmp_path / "bakta", program=str(bakta), verbose=False
    )["result"]["artifacts"]
    assert "checkm2/quality_report.tsv" in predict_checkm2(
        bins, checkm_database, tmp_path / "checkm2", program=str(checkm2), verbose=False
    )["result"]["artifacts"]
    assert "report.tsv" in classify_kraken2(
        reads, kraken_database, tmp_path / "kraken", program=str(kraken2), verbose=False
    )["result"]["artifacts"]


def test_taxonomy_and_unicycler_facades_capture_outputs(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nATGC\n+\nIIII\n", encoding="ascii")
    report = tmp_path / "kraken.report"
    report.write_text("0.0\t0\t0\tU\t0\troot\n", encoding="ascii")
    database = tmp_path / "database"
    database.mkdir()
    index_prefix = tmp_path / "centrifuge-index"
    protein_database = tmp_path / "kaiju.fmi"
    taxonomy_nodes = tmp_path / "nodes.dmp"
    checkm_database = tmp_path / "unused.dmnd"
    protein_database.write_text("index", encoding="ascii")
    taxonomy_nodes.write_text("nodes", encoding="ascii")
    checkm_database.write_text("db", encoding="ascii")
    bracken = _fake_program(tmp_path / "fake-bracken")
    centrifuge = _fake_program(tmp_path / "fake-centrifuge")
    kaiju = _fake_program(tmp_path / "fake-kaiju")
    unicycler = _fake_program(tmp_path / "fake-unicycler")

    assert "bracken.tsv" in estimate_bracken_abundance(
        report, database, tmp_path / "bracken", program=str(bracken), verbose=False
    )["result"]["artifacts"]
    assert "assignments.tsv" in classify_centrifuge(
        index_prefix, reads, tmp_path / "centrifuge", program=str(centrifuge), verbose=False
    )["result"]["artifacts"]
    assert "assignments.tsv" in classify_kaiju(
        reads, protein_database, taxonomy_nodes, tmp_path / "kaiju", program=str(kaiju), verbose=False
    )["result"]["artifacts"]
    assert "unicycler/assembly.fasta" in assemble_unicycler(
        reads, tmp_path / "unicycler", program=str(unicycler), verbose=False
    )["result"]["artifacts"]


def test_gtdbtk_humann_and_methyldackel_facades_capture_outputs(tmp_path: Path) -> None:
    reads = tmp_path / "sample.fastq"
    reads.write_text("@r\nATGC\n+\nIIII\n", encoding="ascii")
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "one.fna").write_text(">one\nATGC\n", encoding="ascii")
    database = tmp_path / "database"
    database.mkdir()
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nATGC\n", encoding="ascii")
    alignment = tmp_path / "alignment.bam"
    alignment.write_text("BAM", encoding="ascii")
    gtdbtk = _fake_program(tmp_path / "fake-gtdbtk")
    humann = _fake_program(tmp_path / "fake-humann")
    methyldackel = _fake_program(tmp_path / "fake-MethylDackel")

    assert "gtdbtk/gtdbtk.bac120.summary.tsv" in classify_gtdbtk(
        genomes, database, tmp_path / "gtdbtk", program=str(gtdbtk), verbose=False
    )["result"]["artifacts"]
    assert "humann/sample_genefamilies.tsv" in profile_humann(
        reads, tmp_path / "humann", program=str(humann), verbose=False
    )["result"]["artifacts"]
    methylation = extract_methylation_methyldackel(
        reference, alignment, tmp_path / "methylation", program=str(methyldackel), verbose=False
    )
    assert "methylation_CpG.bedGraph" in methylation["result"]["artifacts"]
    assert "methylation_CHG.bedGraph" in methylation["result"]["artifacts"]


def test_new_facades_are_agent_compatible_external_tools() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.alignment.search_diamond_blastp",
        "dnakit.alignment.search_diamond_blastx",
        "dnakit.alignment.search_mmseqs_protein",
        "dnakit.alignment.cluster_mmseqs_protein",
        "dnakit.annotation.annotate_prokka",
        "dnakit.annotation.annotate_bakta",
        "dnakit.annotation.predict_checkm2",
        "dnakit.microbiome.classify_kraken2",
        "dnakit.microbiome.classify_centrifuge",
        "dnakit.microbiome.classify_kaiju",
        "dnakit.microbiome.estimate_bracken_abundance",
        "dnakit.assembly.assemble_unicycler",
        "dnakit.microbiome.classify_gtdbtk",
        "dnakit.microbiome.profile_humann",
        "dnakit.epigenome.extract_methylation_methyldackel",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True
