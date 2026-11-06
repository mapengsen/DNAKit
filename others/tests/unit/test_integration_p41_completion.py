"""Offline P41 comparisons against installed upstream implementations."""

import subprocess
import sys
from pathlib import Path

import pytest

from dnakit.alignment.edlib_p41_functions import edlib_prefix_alignment
from dnakit.annotation.ncls_p41_functions import ncls_limited_overlaps
from dnakit.annotation.pybigwig_p41_functions import inspect_bigwig_file
from dnakit.epigenome.epigenome_p41_cli_functions import (
    bam_region_count_matrix,
    bigwig_region_score_matrix,
    macs3_bedgraph_cutoff_report,
    macs3_duplicate_statistics,
    pairtools_tile_duplication_stats,
)
from dnakit.molecular_lab.dna_features_viewer_p41_functions import (
    export_feature_map_record,
    plot_feature_map,
)
from dnakit.phylogeny.pyvolve_p41_functions import (
    pyvolve_state_frequencies,
    simulate_pyvolve_sequences,
)
from dnakit.regulation.tangermeme_p41_functions import (
    tangermeme_extract_signal,
    tangermeme_hypothetical_attributions,
    tangermeme_position_entropy,
    tangermeme_unchunk,
)
from dnakit.sequence.py2bit_p41_functions import twobit_masked_blocks, twobit_region_composition
from dnakit.sequence.screed_p41_functions import (
    screed_is_dna,
    screed_record_at_index,
    screed_record_names,
    screed_write_record,
)
from dnakit.sequencing_qc.pyslow5_p41_functions import slow5_read_records
from dnakit.variants.cyvcf2_p41_functions import (
    cyvcf2_export_text,
    cyvcf2_ibd_bins,
    cyvcf2_informative_sites,
    vcf_call_is_filtered,
)

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)


def exe(name):
    path = Path(sys.executable).parent / name
    if not path.is_file():
        pytest.skip(name + " missing")
    return str(path)


@pytest.mark.parametrize("q,t", [("ACGT", "ACGTAAA"), ("A", "TTT"), ("", "AC")])
@pytest.mark.parametrize("task", ["distance", "locations", "path"])
def test_edlib(q, t, task):
    r = pytest.importorskip("edlib").align(q, t, mode="SHW", task=task)
    a = edlib_prefix_alignment(q, t, task=task, verbose=False)
    assert (
        a["edit_distance"] == r["editDistance"]
        and a["locations"] == [list(x) for x in r["locations"]]
        and a["cigar"] == r["cigar"]
    )


@pytest.mark.parametrize("mode", ["first", "last", "k"])
@pytest.mark.parametrize("qs,qe", [([1], [2]), ([8], [10]), ([], [])])
def test_ncls(mode, qs, qe):
    np = pytest.importorskip("numpy")
    m = pytest.importorskip("ncls")
    s = [0, 2, 3]
    e = [4, 5, 8]
    native = m.NCLS(
        np.array(s, dtype="int64"), np.array(e, dtype="int64"), np.arange(3, dtype="int64")
    )
    args = tuple(np.array(x, dtype="int64") for x in (qs, qe, list(range(len(qs)))))
    r = (
        native.k_overlaps_both(*args, 2)
        if mode == "k"
        else getattr(native, mode + "_overlap_both")(*args)
    )
    a = ncls_limited_overlaps(s, e, qs, qe, selection=mode, k=2, verbose=False)
    assert a["query_indices"] == r[0].tolist() and a["stored_indices"] == r[1].tolist()


@pytest.mark.parametrize("s", ["ACGTN", "acgt", ""])
def test_screed_alphabet(s):
    assert screed_is_dna(s, verbose=False)["is_dna"] == bool(
        pytest.importorskip("screed.dna").is_DNA(s)
    )


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
def test_screed_files(fmt, tmp_path):
    m = pytest.importorskip("screed")
    from dnakit.fasta_fastq.screed_gap_functions import create_screed_database

    files = []
    for i, s in enumerate(["ACGT", "N", "A"]):
        f = tmp_path / f"{i}.{fmt}"
        screed_write_record(
            dict(name=str(i), sequence=s, **({"quality": "I" * len(s)} if fmt == "fastq" else {})),
            str(f),
            file_format=fmt,
            verbose=False,
        )
        assert next(iter(m.open(str(f))))["sequence"] == s
        files.append(f)
    f = tmp_path / "all"
    f.write_text("".join(x.read_text() for x in files))
    db = create_screed_database(str(f), str(tmp_path / "db"), verbose=False)["database"]
    r = m.ScreedDB(db)
    try:
        assert screed_record_names(db, verbose=False)["names"] == list(r.keys())
        for i in range(3):
            assert screed_record_at_index(db, i, verbose=False)["record"] == {
                k: v if isinstance(v, (str, int, float, bool, type(None))) else str(v)
                for k, v in r.loadRecordByIndex(i).items()
            }
    finally:
        r.close()


@pytest.mark.parametrize("chrom,start,end", [("chr1", 0, 150), ("chr1", 24, 74), ("chr2", 0, 1)])
def test_twobit(chrom, start, end):
    m = pytest.importorskip("py2bit")
    f = Path(m.__file__).parent / "py2bitTest/foo.2bit"
    if not f.is_file():
        pytest.skip("upstream fixture absent")
    r = m.open(str(f), storeMasked=True)
    try:
        for frac in [True, False]:
            assert twobit_region_composition(
                str(f), chrom, start=start, end=end, fraction=frac, verbose=False
            )["bases"] == r.bases(chrom, start, end, fraction=frac)
        for kind in ["hard", "soft"]:
            assert twobit_masked_blocks(str(f), chrom, kind=kind, verbose=False)["blocks"] == [
                list(x) for x in getattr(r, kind + "MaskedBlocks")(chrom)
            ]
    finally:
        r.close()


@pytest.mark.parametrize("n", [1, 3, 7])
def test_tangermeme(n):
    t = pytest.importorskip("torch")
    u = pytest.importorskip("tangermeme.utils")
    d = pytest.importorskip("tangermeme.deep_lift_shap")
    pd = pytest.importorskip("pandas")
    x = t.arange(4 * n, dtype=t.float32).reshape(1, 4, n)
    p = t.full((1, 4, n), 0.25)
    assert (
        tangermeme_position_entropy(p.tolist(), verbose=False)["entropy"] == u.entropy(p).tolist()
    )
    assert (
        tangermeme_extract_signal(x.tolist(), [[0, 0, n]], verbose=False)["signal_sums"]
        == u.extract_signal(pd.DataFrame([[0, 0, n]]), x).tolist()
    )
    z = t.zeros_like(x)
    o = z.clone()
    o[:, 0] = 1
    assert (
        tangermeme_hypothetical_attributions(x.tolist(), o.tolist(), z.tolist(), verbose=False)[
            "attributions"
        ]
        == d.hypothetical_attributions((x,), (o,), (z,))[0].tolist()
    )
    chunks = t.arange(2 * n, dtype=t.float32).reshape(2, 1, n)
    assert tangermeme_unchunk(chunks.tolist(), [2 * n], verbose=False)["sequences"] == [
        v.tolist() for v in u.unchunk(chunks, [2 * n])
    ]


@pytest.mark.parametrize(
    "model,params",
    [("gy", {"omega": 0.5}), ("ecmrest", {}), ("mutsel", {"state_freqs": [0.25] * 4})],
)
def test_pyvolve(model, params):
    m = pytest.importorskip("pyvolve")
    tree = "(a:0.1,b:0.1);"
    a = simulate_pyvolve_sequences(
        tree, model=model, parameters=params, size=3, seed=17, verbose=False
    )
    r = m.Evolver(
        tree=m.read_tree(tree=tree),
        partitions=m.Partition(models=m.Model(model, parameters=params), size=3),
    )
    r(seqfile=None, ratefile=None, infofile=None, seed=17)
    assert a["sequences"] == r.get_sequences()


@pytest.mark.parametrize("s", ["ACGT", "AAAACG", "ATGAAATTT"])
def test_frequencies(s, tmp_path):
    m = pytest.importorskip("pyvolve")
    f = tmp_path / "x.fa"
    f.write_text(">s\n" + s + "\n")
    r = m.ReadFrequencies("nucleotide", file=str(f), format="fasta").compute_frequencies()
    assert (
        list(pyvolve_state_frequencies(str(f), verbose=False)["frequencies"].values()) == r.tolist()
    )


@pytest.mark.parametrize("s", ["ATG", "ATGAAATAG", "ATG" * 5])
def test_plots(s, tmp_path):
    pytest.importorskip("dna_features_viewer")
    bio = pytest.importorskip("Bio.SeqIO")
    features = [{"start": 0, "end": 3, "label": "gene"}]
    f = tmp_path / "x.gb"
    export_feature_map_record(s, features, str(f), verbose=False)
    r = bio.read(f, "genbank")
    assert str(r.seq) == s and len(r.features) == 1
    for layout in ["nucleotide", "translation", "multiline", "multipage"]:
        f = tmp_path / (layout + (".pdf" if layout == "multipage" else ".png"))
        plot_feature_map(s, features, str(f), layout=layout, nucleotides_per_line=3, verbose=False)
        assert f.stat().st_size > 0


@pytest.mark.parametrize("ft", [None, [], ["PASS"], ["q10"]])
def test_filter(ft):
    m = pytest.importorskip("vcfpy")
    d = {"GT": "0/1"}
    if ft is not None:
        d["FT"] = ft
    assert (
        vcf_call_is_filtered("0/1", filters=ft, verbose=False)["is_filtered"]
        == m.Call("s", d).is_filtered()
    )


@pytest.mark.parametrize("n", [1, 2, 3])
def test_slow5(n, tmp_path):
    m = pytest.importorskip("pyslow5")
    from dnakit.sequencing_qc.backends import write_slow5_signals

    f = tmp_path / "x.slow5"
    write_slow5_signals(
        [
            {
                "read_id": f"r{i}",
                "signal": [1, 2, 3][: i + 1],
                "digitisation": 8192,
                "offset": 0,
                "range": 1467,
                "sampling_rate": 4000,
            }
            for i in range(n)
        ],
        f,
        verbose=False,
    )
    r = m.Open(str(f), "r")
    try:
        rows = list(r.seq_reads(pA=False))
    finally:
        r.close()
    assert slow5_read_records(str(f), include_signal=False, verbose=False)["records"] == [
        r["read_id"] for r in rows
    ]
    assert [r["signal"] for r in slow5_read_records(str(f), verbose=False)["records"]] == [
        r["signal"].tolist() for r in rows
    ]


@pytest.mark.parametrize("policy", [1, 2, "all"])
def test_macs_stats(policy, tmp_path):
    f = tmp_path / "x.bed"
    f.write_text("chr1\t0\t10\ta\t0\t+\n" * 2 + "chr1\t20\t30\tb\t0\t+\n")
    a = macs3_duplicate_statistics(f, file_format="BED", keep_duplicates=policy, verbose=False)
    r = subprocess.run(
        [
            exe("macs3"),
            "filterdup",
            "-i",
            str(f),
            "-f",
            "BED",
            "-g",
            "hs",
            "--keep-dup",
            str(policy),
            "--dry-run",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (
        a["statistics"]["reads_before"] == 3
        and a["statistics"]["reads_after"] == (2 if policy == 1 else 3)
        and "total tags in alignment file: 3" in r.stderr
    )


@pytest.mark.parametrize("steps", [1, 3, 5])
def test_cutoff(steps, tmp_path):
    f = tmp_path / "x.bdg"
    f.write_text("chr1\t0\t100\t5\nchr1\t100\t200\t10\n")
    o = tmp_path / "o.tsv"
    d = tmp_path / "d.tsv"
    a = macs3_bedgraph_cutoff_report(f, o, min_length=10, cutoff_steps=steps, verbose=False)
    subprocess.run(
        [
            exe("macs3"),
            "bdgpeakcall",
            "-i",
            str(f),
            "--cutoff-analysis",
            "--outdir",
            str(tmp_path),
            "-o",
            d.name,
            "-l",
            "10",
            "--cutoff-analysis-steps",
            str(steps),
        ],
        check=True,
        capture_output=True,
    )
    assert a["status"] == "completed" and o.read_text() == d.read_text()


@pytest.mark.parametrize("value", [0.0, 1.5, -2.0])
def test_bigwig(value, tmp_path):
    m = pytest.importorskip("pyBigWig")
    np = pytest.importorskip("numpy")
    f = tmp_path / "x.bw"
    r = m.open(str(f), "w")
    r.addHeader([("chr1", 100)])
    r.addEntries(["chr1"], [0], ends=[10], values=[value])
    r.close()
    assert inspect_bigwig_file(str(f), verbose=False)["is_bigwig"]
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t0\t5\nchr1\t5\t10\n")
    o = tmp_path / "o.npz"
    d = tmp_path / "d.npz"
    a = bigwig_region_score_matrix([f], [bed], o, verbose=False)
    subprocess.run(
        [
            exe("multiBigwigSummary"),
            "BED-file",
            "--bwfiles",
            str(f),
            "--BED",
            str(bed),
            "-o",
            str(d),
            "--numberOfProcessors",
            "1",
        ],
        check=True,
        capture_output=True,
    )
    with np.load(d) as r:
        assert a["matrix_preview"] == r["matrix"].tolist()
    assert a["labels"] == ["x.bw"]


@pytest.mark.parametrize("end", [1, 100, 1000])
def test_bam(end, tmp_path):
    np = pytest.importorskip("numpy")
    bed = tmp_path / "x.bed"
    bed.write_text(f"chr1\t0\t{end}\n")
    f = DATA / "reads.bam"
    o = tmp_path / "o.npz"
    d = tmp_path / "d.npz"
    a = bam_region_count_matrix([f], [bed], o, verbose=False)
    subprocess.run(
        [
            exe("multiBamSummary"),
            "BED-file",
            "--bamfiles",
            str(f),
            "--BED",
            str(bed),
            "-o",
            str(d),
            "--numberOfProcessors",
            "1",
        ],
        check=True,
        capture_output=True,
    )
    with np.load(d) as r:
        assert a["matrix_preview"] == r["matrix"].tolist()


@pytest.mark.parametrize("count", [0, 1, 3])
def test_cyvcf2(count, tmp_path):
    m = pytest.importorskip("cyvcf2")
    f = tmp_path / "x.vcf"
    header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n'
    f.write_text(
        header
        + "".join(
            f"chr1\t{1 + 6000 * i}\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:20\t0/0:20\n"
            for i in range(count)
        )
    )
    r = m.VCF(str(f), gts012=True)
    try:
        expected = r.ibd(nmax=100)
    finally:
        r.close()
    assert cyvcf2_ibd_bins(str(f), max_variants=100, verbose=False)["pairs"] == [
        {"sample_a": a, "sample_b": b, "bins": x.tolist()} for (a, b), x in expected.items()
    ]
    r = m.VCF(str(f), gts012=True)
    try:
        expected = [(i, v.POS) for i, v in r.gen_variants(None)]
    finally:
        r.close()
    assert [
        (x["index"], x["position"])
        for x in cyvcf2_informative_sites(str(f), verbose=False)["records"]
    ] == expected
    out = tmp_path / "out.vcf"
    cyvcf2_export_text(str(f), str(out), verbose=False)
    r = m.VCF(str(f))
    try:
        expected = r.raw_header + "".join(str(v) for v in r)
    finally:
        r.close()
    assert out.read_text() == expected


@pytest.mark.parametrize("count", [0, 1, 3])
def test_gff(count, tmp_path):
    pytest.importorskip("BCBio.GFF")
    upstream = pytest.importorskip("dna_features_viewer.biotools")
    from dnakit.molecular_lab.dna_features_viewer_p41_functions import read_gff_feature_record

    f = tmp_path / "x.gff"
    f.write_text(
        "##gff-version 3\n##sequence-region chr1 1 100\n"
        + "".join(
            f"chr1\ttest\tgene\t{1 + i * 10}\t{5 + i * 10}\t.\t+\t.\tID=g{i}\n"
            for i in range(count)
        )
        + "##FASTA\n>chr1\n"
        + "A" * 100
        + "\n"
    )
    native = upstream.load_record(str(f), filetype="gff")
    a = read_gff_feature_record(str(f), verbose=False)
    assert a["id"] == native.id and a["length"] == len(native)
    assert [(r["start"], r["end"]) for r in a["features"]] == [
        (int(r.location.start), int(r.location.end)) for r in native.features
    ]


@pytest.mark.parametrize("duplicates", [0, 1, 2])
def test_tile_stats(duplicates, tmp_path):
    source = tmp_path / "x.pairs"
    source.write_text(
        "## pairs format v1.0.0\n#chromsize: chr1 1000\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type parent_readID\n"
        + "M:1:F:1:1101:1:1\tchr1\t10\tchr1\t100\t+\t-\tUU\t\n"
        + "".join(
            f"M:1:F:1:1102:{i}:1\tchr1\t10\tchr1\t100\t+\t-\tDD\tM:1:F:1:1101:1:1\n"
            for i in range(duplicates)
        )
    )
    out, tiles = tmp_path / "stats.tsv", tmp_path / "tiles.tsv"
    actual = pairtools_tile_duplication_stats(source, out, tile_output_path=tiles, verbose=False)
    direct, tmp_tiles = tmp_path / "direct.tsv", tmp_path / "directtiles.tsv"
    subprocess.run(
        [
            exe("pairtools"),
            "stats",
            "--bytile-dups",
            "--output",
            str(direct),
            "--output-bytile-stats",
            str(tmp_tiles),
            str(source),
        ],
        check=True,
        capture_output=True,
    )
    assert actual["status"] == "completed" and tiles.read_text() == tmp_tiles.read_text()


def test_network_contract(tmp_path, monkeypatch):
    import types

    import dnakit.databases.gget_p41_functions as functions
    import dnakit.download.files as downloader

    calls = []
    fake = types.SimpleNamespace(ref=lambda *a, **k: ["https://ftp.ensembl.org/pub/test.fa.gz"])
    monkeypatch.setattr(functions, "require", lambda *args: fake)

    def download(url, folder, *, filename, verbose):
        calls.append((url, Path(folder), filename))
        (Path(folder) / filename).write_bytes(b"toy")

    monkeypatch.setattr(downloader, "download_url", download)
    actual = functions.gget_download_reference("homo_sapiens", str(tmp_path), verbose=False)
    assert (
        actual["status"] == "completed" and calls[0][1] == tmp_path and calls[0][2] == "test.fa.gz"
    )
    from dnakit.exceptions import ConfigurationError

    with pytest.raises(ConfigurationError):
        functions.gget_download_reference("human", str(tmp_path), verbose=False)


def test_prokka_contract(tmp_path, monkeypatch):
    from dnakit.annotation import prokka_p41_cli_functions as functions

    source = tmp_path / "x.fa"
    source.write_text(">s\nACGT\n")
    calls = []

    def run(name, args, **kwargs):
        calls.append((name, args, kwargs))
        return {"status": "failed", "outputs": []}

    monkeypatch.setattr(functions.cli, "run", run)
    result = functions.annotate_prokka_rnammer(source, tmp_path / "out", verbose=False)
    assert result["status"] == "failed" and calls[0][0] == "prokka" and "--rnammer" in calls[0][1]
