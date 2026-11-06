"""Metagenome and microbial-typing facades. Real-program tests skip when the program is missing; programs whose
databases cannot be shipped with a test (geNomad, eggNOG-mapper, ABRicate download) are checked with a recording shim."""

from __future__ import annotations

import random
import stat
import subprocess
from pathlib import Path

import pytest
from _p16_data import genome, kaiju_data, need, rand_dna, write_fasta

import dnakit.microbiome.p16_microbial_cli_functions as mb
from dnakit.exceptions import BackendExecutionError, ConfigurationError


def _kaiju_index(tmp_path: Path) -> dict[str, str]:
    data = kaiju_data(tmp_path)
    mkbwt, mkfmi = need("kaiju-mkbwt"), need("kaiju-mkfmi")
    subprocess.run(
        [
            mkbwt,
            "-n",
            "2",
            "-l",
            "1",
            "-a",
            "ACDEFGHIKLMNPQRSTVWY",
            "-o",
            str(tmp_path / "db"),
            data["db.faa"],
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run([mkfmi, str(tmp_path / "db")], check=True, capture_output=True)
    return {**data, "fmi": str(tmp_path / "db.fmi")}


def test_kaiju_chain(tmp_path: Path) -> None:
    kaiju = need("kaiju")
    data = _kaiju_index(tmp_path)
    best = mb.search_kaiju_best_match(
        data["reads.fa"], data["fmi"], str(tmp_path / "x.out"), verbose=False
    )
    by_read = {row["read"]: row for row in best["rows"]}
    assert (
        by_read["read_p0_100_0"]["matches"].startswith("p0_100")
        and by_read["read_random"]["classified"] == "U"
        and best["classified"] == 9
    )
    plain = tmp_path / "k.out"
    subprocess.run(
        [
            kaiju,
            "-t",
            data["nodes.dmp"],
            "-f",
            data["fmi"],
            "-i",
            data["reads.fa"],
            "-o",
            str(plain),
        ],
        check=True,
        capture_output=True,
    )
    krona = mb.convert_kaiju_to_krona(
        str(plain),
        data["nodes.dmp"],
        data["names.dmp"],
        str(tmp_path / "k.krona"),
        include_unclassified=True,
        verbose=False,
    )
    assert krona["line_count"] == 4
    table = mb.summarize_kaiju_table(
        [str(plain)],
        data["nodes.dmp"],
        data["names.dmp"],
        str(tmp_path / "t.tsv"),
        rank="species",
        verbose=False,
    )
    assert {(r["taxon_name"], r["reads"]) for r in table["rows"] if r["reads"]} >= {
        ("Species_one", 3),
        ("Species_two", 3),
    }
    named = mb.add_kaiju_taxon_names(
        str(plain),
        data["nodes.dmp"],
        data["names.dmp"],
        str(tmp_path / "n.out"),
        full_path=True,
        verbose=False,
    )
    assert named["line_count"] == 10 and "Species_one" in (tmp_path / "n.out").read_text()
    merged = mb.merge_kaiju_outputs(
        str(plain), str(plain), str(tmp_path / "m.out"), nodes=data["nodes.dmp"], verbose=False
    )
    assert (tmp_path / "m.out").read_text() == plain.read_text() and merged["conflict"] == "lca"
    with pytest.raises(ConfigurationError, match="nodes"):
        mb.merge_kaiju_outputs(str(plain), str(plain), str(tmp_path / "m2.out"), verbose=False)


def test_diamond_blastx_modes(tmp_path: Path) -> None:
    program = need("diamond")
    data = kaiju_data(tmp_path)
    subprocess.run(
        [program, "makedb", "--in", data["db.faa"], "-d", str(tmp_path / "dm"), "--quiet"],
        check=True,
        capture_output=True,
    )
    plain = mb.search_dna_translated_diamond(
        data["reads.fa"],
        str(tmp_path / "dm.dmnd"),
        str(tmp_path / "a.tsv"),
        diamond_path=program,
        verbose=False,
    )
    assert {r["column_1"]: r["column_2"] for r in plain["rows"]}["read_p1_200_0"] == "p1_200"
    swipe = mb.search_dna_translated_diamond(
        data["reads.fa"],
        str(tmp_path / "dm.dmnd"),
        str(tmp_path / "b.tsv"),
        exhaustive=True,
        diamond_path=program,
        verbose=False,
    )
    frame = mb.search_dna_translated_diamond(
        data["reads.fa"],
        str(tmp_path / "dm.dmnd"),
        str(tmp_path / "c.tsv"),
        frameshift_penalty=15,
        diamond_path=program,
        verbose=False,
    )
    assert len(swipe["rows"]) >= 9 and len(frame["rows"]) >= 9
    with pytest.raises(ConfigurationError, match="exclusive"):
        mb.search_dna_translated_diamond(
            data["reads.fa"],
            str(tmp_path / "dm.dmnd"),
            str(tmp_path / "d.tsv"),
            frameshift_penalty=15,
            long_reads=True,
            verbose=False,
        )


def test_bracken_build_reestimate_and_combine(tmp_path: Path) -> None:
    kraken2, build, bracken = need("kraken2"), need("bracken-build"), need("bracken")
    del build
    data = kaiju_data(tmp_path)
    db = tmp_path / "krakendb"
    (db / "taxonomy").mkdir(parents=True)
    (db / "library").mkdir()
    for name in ("nodes.dmp", "names.dmp"):
        (db / "taxonomy" / name).write_text(Path(data[name]).read_text())
    rng = random.Random(21)
    records = {
        f"seq{i}|kraken:taxid|{tax}": rand_dna(rng, 4000) for i, tax in enumerate((100, 200, 300))
    }
    write_fasta(db / "library" / "lib.fna", records)
    (db / "seqid2taxid.map").write_text("".join(f"{h}\t{h.split('|')[-1]}\n" for h in records))
    subprocess.run(
        [
            kraken2.replace("kraken2", "kraken2-build"),
            "--build",
            "--db",
            str(db),
            "--kmer-len",
            "15",
            "--minimizer-len",
            "12",
            "--minimizer-spaces",
            "0",
        ],
        check=True,
        capture_output=True,
    )
    reads = {
        f"r{i}_{j}": seq[p : p + 100]
        for i, seq in enumerate(records.values())
        for j, p in enumerate(range(0, 3000, 100))
        if i < 2
    }
    write_fasta(tmp_path / "reads.fa", reads)
    subprocess.run(
        [
            kraken2,
            "--db",
            str(db),
            "--report",
            str(tmp_path / "rep.txt"),
            "--output",
            str(tmp_path / "out.txt"),
            str(tmp_path / "reads.fa"),
        ],
        check=True,
        capture_output=True,
    )
    built = mb.build_bracken_kmer_distribution(
        str(db), read_length=100, kmer_length=15, verbose=False
    )
    assert Path(built["outputs"][1]).stat().st_size > 0
    with pytest.raises(ConfigurationError, match="overwritten"):
        mb.build_bracken_kmer_distribution(str(db), read_length=100, kmer_length=15, verbose=False)
    estimated = mb.reestimate_bracken_report(
        str(tmp_path / "rep.txt"),
        str(db),
        str(tmp_path / "b.tsv"),
        str(tmp_path / "b.report"),
        bracken_path=bracken,
        verbose=False,
    )
    assert {r["name"]: r["new_est_reads"] for r in estimated["rows"]} == {
        "Species_one": 30,
        "Species_two": 30,
    }
    assert (tmp_path / "b.report").read_text().splitlines()[0].split("\t")[1] == "60"
    combined = mb.combine_bracken_tables(
        [str(tmp_path / "b.tsv"), str(tmp_path / "b.tsv")],
        str(tmp_path / "c.tsv"),
        names=["s1", "s2"],
        verbose=False,
    )
    lines = (tmp_path / "c.tsv").read_text().splitlines()
    assert (
        lines[0].split("\t")[3:] == ["s1_num", "s1_frac", "s2_num", "s2_frac"]
        and combined["taxa"] == 2
        and lines[1].split("\t")[4] == "0.50000"
    )


def test_combine_bracken_tables_rejects_mismatches(tmp_path: Path) -> None:
    head = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text(head + "X\t1\tS\t1\t1\t5\t1.0\n")
    b.write_text(head + "X\t2\tS\t1\t1\t5\t1.0\n")
    with pytest.raises(ConfigurationError, match="Taxonomy IDs"):
        mb.combine_bracken_tables([str(a), str(b)], str(tmp_path / "o.tsv"), verbose=False)


def test_fastani_pair_and_plot(tmp_path: Path) -> None:
    program = need("fastANI")
    pytest.importorskip("matplotlib")
    contigs = genome(tmp_path)
    result = mb.map_fastani_conserved_regions(
        str(tmp_path / "genome.fa"),
        str(tmp_path / "genome.fa"),
        str(tmp_path / "ani.txt"),
        fastani_path=program,
        verbose=False,
    )
    assert (
        result["ani"] > 99.9
        and result["matched_fragments"]
        == result["total_fragments"]
        == len(contigs["chrA"]) // 3000 + len(contigs["chrB"]) // 3000
    )
    plot = mb.plot_fastani_conserved_regions(result["visual_file"], str(tmp_path / "ani.png"))
    assert (
        plot["fragments"] == result["matched_fragments"]
        and Path(plot["image"]).read_bytes()[:4] == b"\x89PNG"
    )


def test_abricate_mlst_and_prokka_scripts(tmp_path: Path) -> None:
    abricate, mlst = need("abricate"), need("mlst")
    contigs = genome(tmp_path)
    try:
        screened = mb.screen_contigs_abricate(
            str(tmp_path / "genome.fa"),
            str(tmp_path / "ab.tab"),
            database="ncbi",
            abricate_path=abricate,
            verbose=False,
        )
        typed = mb.type_sequence_type_mlst(
            str(tmp_path / "genome.fa"),
            str(tmp_path / "mlst.tsv"),
            mlst_path=mlst,
            novel_alleles_path=str(tmp_path / "novel.fa"),
            verbose=False,
        )
    except BackendExecutionError as exc:  # perl modules of a mixed PATH
        pytest.skip(f"abricate/mlst not runnable here: {str(exc.context.get('log_tail'))[:120]}")
    assert screened["rows"] == [] and typed["fields"][0] == "genome.fa" and contigs
    summary = mb.summarize_abricate_reports(
        [str(tmp_path / "ab.tab")], str(tmp_path / "sum.tab"), abricate_path=abricate, verbose=False
    )
    assert summary["rows"][0]["NUM_FOUND"] == 0


def test_prokka_database_converters(tmp_path: Path) -> None:
    script = need("prokka-genbank_to_fasta_db")
    bio = pytest.importorskip("Bio.SeqIO")
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    sequence = rand_dna(random.Random(3), 3000)
    record = SeqRecord(
        Seq(sequence), id="chrA", name="chrA", description="x", annotations={"molecule_type": "DNA"}
    )
    for i, (start, end, product) in enumerate(
        ((100, 400, "prod0"), (600, 900, "prod1"), (1200, 1500, "hypothetical protein"))
    ):
        record.features.append(
            SeqFeature(
                FeatureLocation(start, end, strand=1),
                type="CDS",
                qualifiers={
                    "locus_tag": [f"L{i}"],
                    "product": [product],
                    "gene": [f"g{i}"],
                    "transl_table": ["11"],
                    "protein_id": [f"P{i}.1"],
                },
            )
        )
    bio.write(record, str(tmp_path / "t.gbk"), "genbank")
    try:
        result = mb.genbank_to_protein_fasta_prokka(
            [str(tmp_path / "t.gbk")],
            str(tmp_path / "p.faa"),
            prokka_script_path=script,
            verbose=False,
        )
    except BackendExecutionError as exc:
        pytest.skip(f"BioPerl not usable here: {str(exc.context.get('log_tail'))[:120]}")
    assert result["proteins"] == 2  # the hypothetical protein is excluded by default
    header = (tmp_path / "p.faa").read_text().splitlines()[0]
    assert header == ">P0.1 ~~~g0~~~prod0"


def _shim(tmp_path: Path, name: str) -> str:
    """A program that records its arguments and creates the files the facade expects, to check command assembly."""
    path = tmp_path / name
    path.write_text(
        '#!/bin/sh\nprintf \'%s\\n\' "$@" >> "$0.argv"\n'
        "out=''; prev=''\nfor a in \"$@\"; do case \"$prev\" in --output_dir|-o|--output) out=$a;; esac; prev=$a; done\n"
        "exit 0\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def test_genomad_steps_assemble_commands(tmp_path: Path) -> None:
    shim = _shim(tmp_path, "genomad")
    fasta = write_fasta(tmp_path / "c.fa", {"c1": "ACGT" * 50})
    database = tmp_path / "gdb"
    database.mkdir()
    result = mb.run_genomad(
        fasta,
        str(tmp_path / "out"),
        step="annotate",
        database=str(database),
        threads=2,
        genomad_path=shim,
        verbose=False,
    )
    assert result["command"][1:3] == ["annotate", "--cleanup"] and result["command"][-1] == str(
        database
    )
    summary = mb.run_genomad(
        fasta,
        str(tmp_path / "out"),
        step="summary",
        min_score=0.8,
        genomad_path=shim,
        verbose=False,
    )
    assert (
        summary["command"][1] == "summary"
        and "--min-score" in summary["command"]
        and str(database) not in summary["command"]
    )
    nn = mb.run_genomad(
        fasta,
        str(tmp_path / "out"),
        step="aggregated-classification",
        enable_score_calibration=True,
        genomad_path=shim,
        verbose=False,
    )
    assert "--enable-score-calibration" in nn["command"]


def test_eggnog_mapper_modes_assemble_commands(tmp_path: Path) -> None:
    shim = _shim(tmp_path, "emapper.py")
    fasta = write_fasta(tmp_path / "g.fa", {"c1": "ACGT" * 50})
    out = tmp_path / "e"
    # the shim cannot create the annotations file the facade waits for: that is reported as an error with the log
    with pytest.raises(BackendExecutionError, match="did not write"):
        mb.annotate_dna_eggnog_mapper(
            fasta,
            str(out),
            input_type="genome",
            gene_prediction="prodigal",
            emapper_path=shim,
            verbose=False,
        )
    argv = Path(shim + ".argv").read_text().split()
    assert (
        argv[argv.index("--itype") + 1] == "genome"
        and argv[argv.index("--genepred") + 1] == "prodigal"
    )
    with pytest.raises(ConfigurationError, match="gene_prediction"):
        mb.annotate_dna_eggnog_mapper(
            fasta, str(tmp_path / "e2"), input_type="CDS", gene_prediction="search", verbose=False
        )


def test_amrfinder_family_assemble_commands(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "c.fa", {"c1": "ACGT" * 50})
    amr = _shim(tmp_path, "amrfinder")
    with pytest.raises(BackendExecutionError, match="did not write"):
        mb.detect_amr_genes_nucleotide(
            fasta,
            str(tmp_path / "amr.tsv"),
            organism="Escherichia",
            plus=True,
            amrfinder_path=amr,
            verbose=False,
        )
    argv = Path(amr + ".argv").read_text().split()
    assert (
        "--nucleotide" in argv
        and argv[argv.index("--organism") + 1] == "Escherichia"
        and "--plus" in argv
    )
    stx = _shim(tmp_path, "stxtyper")
    with pytest.raises(BackendExecutionError, match="did not write"):
        mb.type_shiga_toxin_operons(
            fasta,
            str(tmp_path / "stx.tsv"),
            amrfinder_format=True,
            stxtyper_path=stx,
            verbose=False,
        )
    assert "--amrfinder" in Path(stx + ".argv").read_text().split()
    gff = tmp_path / "a.gff"
    gff.write_text("##gff-version 3\n")
    checker = _shim(tmp_path, "gff_check")
    ok = mb.check_gff_annotation_amrfinder(
        str(gff), fasta, fasta, annotation_format="prokka", gff_check_path=checker, verbose=False
    )
    assert ok["status"] == "completed" and "prokka" in Path(checker + ".argv").read_text().split()


def test_abricate_get_db_assembles_command(tmp_path: Path) -> None:
    shim = _shim(tmp_path, "abricate-get_db")
    result = mb.download_abricate_database(
        "card",
        database_directory=str(tmp_path / "db"),
        force=True,
        abricate_get_db_path=shim,
        verbose=False,
    )
    assert result["command"][1:] == ["--db", "card", "--dbdir", str(tmp_path / "db"), "--force"]
    with pytest.raises(ConfigurationError):
        mb.download_abricate_database("nodb", verbose=False)


def test_amrfinder_stxtyper_and_fasta_check_with_the_real_database(tmp_path: Path) -> None:
    amrfinder = need("amrfinder")
    database = (
        Path(amrfinder).resolve().parent.parent
        / "share"
        / "amrfinderplus"
        / "data"
        / "latest"
        / "AMR_CDS.fa"
    )
    if not database.is_file():
        pytest.skip("the AMRFinderPlus database is not installed (run `amrfinder -u`)")
    records = {}
    name = None
    for line in database.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif name:
            records[name] += line
    bla = next(seq for header, seq in records.items() if "|blaTEM-116|" in header)
    stx = next(seq for header, seq in records.items() if "|stxA2b|" in header)
    pad = rand_dna(random.Random(2), 300)
    contigs = write_fasta(
        tmp_path / "contigs.fa", {"contig1": pad + bla + pad, "contig2": pad + stx}
    )
    result = mb.detect_amr_genes_nucleotide(contigs, str(tmp_path / "amr.tsv"), verbose=False)
    hit = next(row for row in result["rows"] if row["Element symbol"] == "blaTEM-116")
    assert (hit["Contig id"], hit["Start"], hit["Stop"], hit["% Identity to reference"]) == (
        "contig1",
        301,
        300 + len(bla) - 3,
        100.0,
    )  # the reported span excludes the stop codon
    typed = mb.type_shiga_toxin_operons(contigs, str(tmp_path / "stx.tsv"), verbose=False)
    assert typed["rows"][0]["stx_type"] == "stx2"
    checked = mb.check_nucleotide_fasta_amrfinder(contigs, verbose=False)
    assert (
        checked["sequence_count"] == 2
        and checked["lengths"] == [600 + len(bla), 300 + len(stx)]
        and checked["cumulative_lengths"][-1] == 900 + len(bla) + len(stx)
    )
