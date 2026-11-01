"""Tests for dnakit.visualization.pygenomeviz_gap_functions (example data, figures and aligners).

The example-data loaders run against a fake of the upstream download, so no request is sent. The
plot and the alignment facades run on small GenBank and FASTA fixtures built in tmp_path.
MMseqs2, MUMmer and progressiveMauve are not installed in every environment: the functional tests of
the facades use stub executables that write the output formats those programs produce, the
missing-program path runs without any binary, and the real-program tests skip with a reason when
their program is absent. One live download runs only with DNAKIT_LIVE_TESTS=1.
"""

from __future__ import annotations

import inspect
import json
import os
import random
import shutil
import socket
import sys
from pathlib import Path

import pytest

# Importing pygenomeviz changes global matplotlib settings (savefig.bbox) for the whole test process.
# Restore them so that the plot-size tests of other modules are not affected by this file.
_mpl = pytest.importorskip("matplotlib")
_saved_rc = _mpl.rcParams.copy()
pgv = pytest.importorskip(
    "pygenomeviz", reason="pygenomeviz is not installed (pip install 'dnakit[bio-pygenomeviz]')"
)
_mpl.rcParams.update(_saved_rc)

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.visualization import pygenomeviz_gap_functions as pg  # noqa: E402

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(autouse=True)
def _no_network(monkeypatch, request):
    if os.environ.get("DNAKIT_LIVE_TESTS") == "1" and "live" in request.node.name:
        return

    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


def _sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _write_genbank(
    path: Path, name: str, sequence: str, cds: list[tuple[int, int, int, str]]
) -> Path:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq(sequence),
        id=name,
        name=name,
        description=f"{name} test genome",
        annotations={"molecule_type": "DNA"},
    )
    for index, (start, end, strand, product) in enumerate(cds, 1):
        record.features.append(
            SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type="CDS",
                qualifiers={
                    "locus_tag": [f"{name}_{index:03d}"],
                    "product": [product],
                    "protein_id": [f"{name}_P{index:03d}.1"],
                    "translation": [product.upper()],
                },
            )
        )
    SeqIO.write([record], str(path), "genbank")
    return path


@pytest.fixture
def genomes(tmp_path: Path) -> list[str]:
    """Two 420 bp genomes; genome A has one reverse-strand CDS."""
    genome_a = _write_genbank(
        tmp_path / "genomeA.gbk",
        "genomeA",
        _sequence(420, 1),
        [(0, 90, 1, "MKPGFK"), (120, 300, -1, "RDIASTRIL")],
    )
    genome_b = _write_genbank(
        tmp_path / "genomeB.gbk",
        "genomeB",
        _sequence(420, 2),
        [(10, 100, 1, "MKPGFK"), (150, 330, 1, "GSEQ")],
    )
    return [str(genome_a), str(genome_b)]


def _stub(path: Path, body: str, log: Path) -> Path:
    """An executable Python script that stands in for an external program and logs its arguments."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"#!{sys.executable}\nLOG = {str(log)!r}\n" + body
    path.write_text(text)
    path.chmod(0o755)
    return path


_RECORD_ARGS = """
import json
with open(LOG, "a") as handle:
    handle.write(json.dumps(sys.argv[1:]) + "\\n")
"""

_MMSEQS_BODY = (
    "import sys\nfrom pathlib import Path\n"
    + _RECORD_ARGS
    + """
args = sys.argv[1:]
assert args[0] == "easy-rbh", args

def first_id(path):
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            return line[1:].split()[0]
    raise SystemExit("empty fasta")

Path(args[3]).write_text(
    f"{first_id(args[1])}\\t{first_id(args[2])}"
    "\\t0.95215\\t100\\t0\\t0\\t1\\t90\\t1\\t90\\t1e-30\\t200\\n"
)
"""
)

_PROMER_BODY = (
    "import sys\nfrom pathlib import Path\n"
    + _RECORD_ARGS
    + """
args = sys.argv[1:]
assert args[0] == "--mum", args
prefix = [a for a in args if a.startswith("--prefix=")][0][len("--prefix="):]
Path(prefix + ".delta").write_text("/fake/a /fake/b\\nPROMER\\nNUCMER fake delta\\n")
"""
)

_DELTA_FILTER_BODY = (
    "import sys\nfrom pathlib import Path\n"
    + _RECORD_ARGS
    + """
args = sys.argv[1:]
assert args[0] in ("-m", "-1"), args
sys.stdout.write(Path(args[1]).read_text())
"""
)

_SHOW_COORDS_BODY = (
    "import sys\n"
    + _RECORD_ARGS
    + """
args = sys.argv[1:]
assert args[:4] == ["-H", "-T", "-r", "-k"], args
print("NOTICE: stub show-coords")
print("\\t".join(["1", "90", "1", "90", "90", "90", "99.00", "90", "90", "100.00", "100.00",
                  "genomeA", "genomeB"]))
"""
)

_MAUVE_BODY = (
    "import sys\nfrom pathlib import Path\n"
    + _RECORD_ARGS
    + """
args = sys.argv[1:]
bbone = [a for a in args if a.startswith("--backbone-output=")][0].split("=", 1)[1]
Path(bbone).write_text(
    "genomeA\\tgenomeA\\tgenomeB\\tgenomeB\\n"
    "1\\t300\\t1\\t300\\n"
    "301\\t450\\t-450\\t-301\\n"
    "1\\t5\\t1\\t5\\n"
)
"""
)


def _calls(log: Path) -> list[list[str]]:
    if not log.exists():
        return []
    return [json.loads(line) for line in log.read_text().splitlines()]


def test_public_functions_take_keyword_only_options_and_verbose_last():
    assert sorted(pg.__all__) == [
        "compare_genomes_mmseqs",
        "compare_genomes_mummer",
        "compare_genomes_progressive_mauve",
        "load_example_dataset",
        "load_example_gff",
        "plot_genome_features",
    ]
    for name in pg.__all__:
        function = getattr(pg, name)
        parameters = list(inspect.signature(function).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name
        assert function.__doc__ and function.__doc__.strip(), name


def test_load_example_dataset_wraps_upstream_download(monkeypatch, tmp_path):
    from pygenomeviz.utils import DatasetLink

    calls = []

    def fake_download(name, cache_dir=None, overwrite_cache=False):
        calls.append((name, cache_dir, overwrite_cache))
        gbk = Path(cache_dir) / "JX128258.gbk"
        gbk.write_text("LOCUS fake\n")
        return [gbk], [DatasetLink("JX128258", 1, 10, "MH051335", 2, 11, 98.5)]

    monkeypatch.setattr(pgv.utils, "load_example_dataset", fake_download)
    result = pg.load_example_dataset(
        "escherichia_phage", cache_dir=str(tmp_path), overwrite_cache=True
    )
    assert calls == [("escherichia_phage", str(tmp_path), True)]
    assert result["status"] == "completed"
    assert result["dataset"] == "escherichia_phage"
    assert result["source_url"].endswith("/escherichia_phage/")
    assert result["gbk_files"] == [str(tmp_path / "JX128258.gbk")]
    assert result["link_count"] == 1
    assert result["links"][0]["identity"] == 98.5
    assert result["retrieved_at_utc"].endswith("+00:00")
    json.dumps(result)


def test_load_example_dataset_rejects_unknown_name_before_downloading(monkeypatch):
    calls = []
    monkeypatch.setattr(pgv.utils, "load_example_dataset", lambda *a, **k: calls.append(a))
    with pytest.raises(ConfigurationError):
        pg.load_example_dataset("not_a_dataset")
    assert calls == []


def test_load_example_gff_wraps_upstream_download(monkeypatch, tmp_path):
    def fake_download(filename, cache_dir=None, overwrite_cache=False):
        calls.append((filename, cache_dir, overwrite_cache))
        path = Path(cache_dir) / filename
        path.write_text("##gff-version 3\n")
        return path

    calls = []
    monkeypatch.setattr(pgv.utils, "load_example_gff", fake_download)
    result = pg.load_example_gff("enterobacteria_phage.gff", cache_dir=str(tmp_path))
    assert calls == [("enterobacteria_phage.gff", str(tmp_path), False)]
    assert result["gff_path"] == str(tmp_path / "enterobacteria_phage.gff")
    assert result["source_url"].endswith("/gff/enterobacteria_phage.gff")
    assert result["retrieved_at_utc"].endswith("+00:00")


def test_load_example_gff_rejects_unknown_file(monkeypatch):
    calls = []
    monkeypatch.setattr(pgv.utils, "load_example_gff", lambda *a, **k: calls.append(a))
    with pytest.raises(ConfigurationError):
        pg.load_example_gff("unknown.gff")
    assert calls == []


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live download; set DNAKIT_LIVE_TESTS=1"
)
def test_live_example_gff_download(tmp_path):
    result = pg.load_example_gff("enterobacteria_phage.gff", cache_dir=str(tmp_path))
    assert Path(result["gff_path"]).stat().st_size > 0


def test_upstream_feature_extraction_defect_is_reproduced(genomes):
    gbk = pgv.Genbank(genomes[0])
    try:
        gbk.extract_features("CDS")
    except AttributeError as exc:
        assert "strand" in str(exc)
    else:
        pytest.skip("upstream Genbank.extract_features works with this Biopython")


def test_cds_spans_read_the_strand_from_the_location(genomes):
    spans = pg._cds_spans(pgv.Genbank(genomes[0]))
    assert [(span.start, span.end, span.strand) for span in spans] == [(0, 90, 1), (120, 300, -1)]


def test_cds_protein_fasta_uses_the_rbh_id_layout(genomes, tmp_path):
    out = tmp_path / "genomeA.faa"
    pg._write_cds_protein_fasta(pgv.Genbank(genomes[0]), out)
    ids = [line[1:].split()[0] for line in out.read_text().splitlines() if line.startswith(">")]
    assert ids == ["GENE000001_genomeA_P001.1|0_90_1|", "GENE000002_genomeA_P002.1|120_300_-1|"]


def test_parse_rbh_converts_identity_and_skips_duplicate_hits(tmp_path):
    rbh = tmp_path / "rbh.tsv"
    rbh.write_text(
        "GENE000001_genomeA_P001.1|0_90_1|\tGENE000001_genomeB_P001.1|10_100_1|\t0.95215\t100\n"
        "GENE000001_genomeA_P001.1|0_90_1|\tGENE000002_genomeB_P002.1|150_330_1|\t0.90\t100\n"
    )
    coords = pg._parse_rbh(rbh, "genomeA", "genomeB")
    assert len(coords) == 1
    assert coords[0].identity == 95.21
    assert (coords[0].ref_start, coords[0].ref_end) == (0, 90)
    assert (coords[0].query_start, coords[0].query_end) == (10, 100)


def test_parse_progressive_mauve_backbone_keeps_long_blocks_and_flips_reverse_blocks(tmp_path):
    bbone = tmp_path / "bbone.tsv"
    bbone.write_text("A\tA\tB\tB\n1\t300\t1\t300\n301\t450\t-450\t-301\n1\t5\t1\t5\n")
    coords = pg._parse_bbone(bbone, ["genomeA", "genomeB"], 0)
    assert len(coords) == 2
    assert (coords[0].ref_start, coords[0].ref_end, coords[0].query_length) == (1, 300, 300)
    assert (coords[1].query_start, coords[1].query_end) == (301, 450)


def test_plot_genome_features_writes_png_and_returns_records(genomes, tmp_path):
    image = tmp_path / "features.png"
    result = pg.plot_genome_features(genomes, str(image), dpi=80)
    assert result["status"] == "completed"
    assert Path(result["image_path"]) == image.resolve()
    assert image.read_bytes().startswith(PNG_MAGIC)
    summaries = {entry["name"]: entry for entry in result["data"]["genomes"]}
    assert summaries["genomeA"]["cds_count"] == 2
    assert summaries["genomeB"]["cds_count"] == 2
    assert summaries["genomeA"]["range_size"] == 420
    json.dumps(result)


def test_plot_refuses_existing_image_and_non_png_suffix(genomes, tmp_path):
    existing = tmp_path / "existing.png"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        pg.plot_genome_features(genomes, str(existing))
    with pytest.raises(ConfigurationError):
        pg.plot_genome_features(genomes, str(tmp_path / "figure.pdf"))


def test_genome_names_must_be_plain_identifiers(tmp_path):
    spaced = _write_genbank(
        tmp_path / "my genome.gbk", "spaced", _sequence(200, 3), [(0, 90, 1, "MKPGFK")]
    )
    other = _write_genbank(
        tmp_path / "other.gbk", "other", _sequence(200, 4), [(0, 90, 1, "MKPGFK")]
    )
    with pytest.raises(ConfigurationError):
        pg.plot_genome_features([str(spaced), str(other)], str(tmp_path / "x.png"))


def test_comparisons_need_two_genomes_and_a_new_output_folder(genomes, tmp_path):
    with pytest.raises(ConfigurationError):
        pg.compare_genomes_mmseqs(genomes[:1], str(tmp_path / "one"))
    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "old.txt").write_text("x")
    with pytest.raises(ConfigurationError):
        pg.compare_genomes_mummer(genomes, str(busy))


def test_mmseqs_program_missing_raises_with_install_hint(genomes, tmp_path, monkeypatch):
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    out = tmp_path / "mmseqs_out"
    with pytest.raises(BackendUnavailableError) as info:
        pg.compare_genomes_mmseqs(genomes, str(out))
    assert "mmseqs" in str(info.value)
    assert "MMseqs2" in info.value.hint
    assert not out.exists()


def test_mmseqs_explicit_path_that_does_not_exist_raises(genomes, tmp_path):
    with pytest.raises(BackendUnavailableError):
        pg.compare_genomes_mmseqs(
            genomes, str(tmp_path / "out"), mmseqs_path=str(tmp_path / "no_mmseqs")
        )


def test_mummer_programs_missing_raise_with_install_hint(genomes, tmp_path, monkeypatch):
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    with pytest.raises(BackendUnavailableError) as info:
        pg.compare_genomes_mummer(genomes, str(tmp_path / "mummer_out"), mummer_dir=str(empty))
    assert "MUMmer" in info.value.hint


def test_progressive_mauve_missing_raises_with_install_hint(genomes, tmp_path, monkeypatch):
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    with pytest.raises(BackendUnavailableError) as info:
        pg.compare_genomes_progressive_mauve(genomes, str(tmp_path / "mauve_out"))
    assert "progressiveMauve" in str(info.value)
    assert "darlinglab" in info.value.hint


def test_mmseqs_facade_runs_with_a_stub_program(genomes, tmp_path):
    log = tmp_path / "mmseqs_calls.jsonl"
    stub = _stub(tmp_path / "bin" / "mmseqs", _MMSEQS_BODY, log)
    out = tmp_path / "run mmseqs; not a shell"
    result = pg.compare_genomes_mmseqs(
        genomes,
        str(out),
        mmseqs_path=str(stub),
        threads=2,
        dpi=80,
    )
    assert result["status"] == "completed"
    assert result["backend"] == "mmseqs"
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)
    calls = _calls(log)
    assert len(calls) == 1
    assert calls[0][0] == "easy-rbh"
    assert calls[0][-6:] == ["--threads", "2", "-e", "0.001", "-v", "0"]
    assert result["data"]["alignment_count"] == 1
    assert result["data"]["alignments"][0]["identity"] == 95.21
    assert (out / "align_coords.tsv").is_file()


def test_mummer_facade_runs_with_stub_programs(genomes, tmp_path):
    log = tmp_path / "mummer_calls.jsonl"
    bin_dir = tmp_path / "bin"
    _stub(bin_dir / "promer", _PROMER_BODY, log)
    _stub(bin_dir / "delta-filter", _DELTA_FILTER_BODY, log)
    _stub(bin_dir / "show-coords", _SHOW_COORDS_BODY, log)
    result = pg.compare_genomes_mummer(
        genomes, str(tmp_path / "mummer_out"), mummer_dir=str(bin_dir), dpi=80
    )
    assert result["status"] == "completed"
    assert result["backend"] == "mummer"
    assert [call[0] for call in result["commands"]] == ["promer", "delta-filter", "show-coords"]
    assert _calls(log)[1] == ["-m", str(tmp_path / "mummer_out" / "out0.delta")]
    alignments = result["data"]["alignments"]
    assert len(alignments) == 1
    assert (alignments[0]["ref_start"], alignments[0]["ref_end"]) == (1, 90)
    assert alignments[0]["ref_name"] == "genomeA" and alignments[0]["query_name"] == "genomeB"
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)


def test_progressive_mauve_facade_runs_with_a_stub_program(genomes, tmp_path):
    log = tmp_path / "mauve_calls.jsonl"
    stub = _stub(tmp_path / "bin" / "progressiveMauve", _MAUVE_BODY, log)
    result = pg.compare_genomes_progressive_mauve(
        genomes, str(tmp_path / "mauve_out"), progressive_mauve_path=str(stub), dpi=80
    )
    assert result["status"] == "completed"
    assert result["backend"] == "progressiveMauve"
    assert result["data"]["alignment_count"] == 2
    assert [genome["name"] for genome in result["data"]["genomes"]] == ["genomeA", "genomeB"]
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)
    assert (tmp_path / "mauve_out" / "seqfiles" / "genomeA.gbk").is_file()


def test_failed_external_run_returns_status_failed_without_figure(genomes, tmp_path):
    body = "import sys\nsys.stdout.write('boom')\nsys.exit(3)\n"
    stub = _stub(tmp_path / "bin" / "mmseqs", body, tmp_path / "unused.jsonl")
    result = pg.compare_genomes_mmseqs(genomes, str(tmp_path / "failed_out"), mmseqs_path=str(stub))
    assert result["status"] == "failed"
    assert result["return_codes"] == [3]
    assert result["image_path"] is None
    assert "boom" in result["log_tail"]


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs (MMseqs2) is not installed")
def test_mmseqs_real_program(genomes, tmp_path):
    result = pg.compare_genomes_mmseqs(genomes, str(tmp_path / "real_mmseqs"), dpi=80)
    assert result["status"] == "completed"
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)


_MUMMER_PROGRAMS = ("nucmer", "promer", "delta-filter", "show-coords")


@pytest.mark.skipif(
    any(shutil.which(name) is None for name in _MUMMER_PROGRAMS),
    reason="MUMmer (nucmer, promer, delta-filter, show-coords) is not installed",
)
def test_mummer_real_programs(genomes, tmp_path):
    result = pg.compare_genomes_mummer(genomes, str(tmp_path / "real_mummer"), dpi=80)
    assert result["status"] == "completed"
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)


@pytest.mark.skipif(
    shutil.which("progressiveMauve") is None, reason="progressiveMauve is not installed"
)
def test_progressive_mauve_real_program(genomes, tmp_path):
    result = pg.compare_genomes_progressive_mauve(genomes, str(tmp_path / "real_mauve"), dpi=80)
    assert result["status"] == "completed"
    assert Path(result["image_path"]).read_bytes().startswith(PNG_MAGIC)
