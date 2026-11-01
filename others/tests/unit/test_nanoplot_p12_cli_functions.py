"""Tests of the NanoPlot facade: validation, argument vector (stand-in program) and a real NanoPlot run."""

from __future__ import annotations

import gzip
import json
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import nanoplot_p12_cli_functions as npl

needs_nanoplot = pytest.mark.skipif(
    shutil.which("NanoPlot") is None, reason="NanoPlot is not installed"
)


@pytest.fixture()
def fastq(tmp_path: Path) -> Path:
    path = tmp_path / "reads.fastq"
    lines = []
    for i in range(40):
        length = 100 + 7 * i
        lines += [f"@read{i}", "ACGT" * (length // 4), "+", "5" * (length // 4 * 4)]
    path.write_text("\n".join(lines) + "\n")
    return path


def _stand_in(root: Path) -> Path:
    script = root / "NanoPlot"
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, pathlib, sys\n"
        f"open({str(root / 'record.jsonl')!r}, 'a').write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:] == ['--version']:\n    print('NanoPlot 9.9'); sys.exit(0)\n"
        "out = pathlib.Path(sys.argv[sys.argv.index('-o') + 1]); out.mkdir(exist_ok=True)\n"
        "(out / 'NanoStats.txt').write_text('Metrics\\tdataset\\nnumber_of_reads\\t40\\nmean_qual\\t20.0\\n')\n"
        "(out / 'Yield_By_Length.html').write_text('x')\n"
        "(out / 'Non_weightedHistogramReadlength.png').write_text('x')\n"
        "(out / 'NanoPlot-data.tsv.gz').write_bytes(__import__('gzip').compress(b'readIDs\\tlengths\\nr1\\t10\\n'))\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def test_validation(fastq: Path, tmp_path: Path) -> None:
    out = str(tmp_path / "o")
    for kwargs in ({"input_type": "pdf"}, {"plot_names": ["nope"]}, {"bivariate_plots": ["pie"]},
                   {"image_formats": ["gif"]}, {"color": "red;rm"}, {"min_length": 0}, {"show_n50": "yes"}):  # fmt: skip
        with pytest.raises(ConfigurationError):
            npl.nanoplot_plots(
                [str(fastq)], out, nanoplot_path="/nonexistent", verbose=False, **kwargs
            )  # type: ignore[arg-type]
    with pytest.raises(BackendUnavailableError):
        npl.nanoplot_plots([str(fastq)], out, nanoplot_path="/nonexistent/NanoPlot", verbose=False)


def test_argument_vector_and_result(fastq: Path, tmp_path: Path) -> None:
    script = _stand_in(tmp_path)
    result = npl.nanoplot_plots(
        [str(fastq)], str(tmp_path / "out"), static_images=False, raw_table=True, log_length=True,
        min_length=50, bivariate_plots=["dot"], show_n50=False, nanoplot_path=str(script), verbose=False,
    )  # fmt: skip
    call = json.loads((tmp_path / "record.jsonl").read_text().splitlines()[0])
    assert (
        "--no_static" in call and "--raw" in call and "--loglength" in call and "--no-N50" in call
    )
    assert call[call.index("--plots") + 1] == "dot" and call[call.index("--minlength") + 1] == "50"
    assert call[-2:] == ["--fastq", str(fastq)]
    assert set(result["plots"]) == {"yield_by_minimal_length", "read_length_histogram"}
    assert "reads_over_time" in result["missing_plots"]
    assert result["statistics"] == {"number_of_reads": 40, "mean_qual": 20.0}
    assert result["raw_table_columns"] == ["readIDs", "lengths"] and result["raw_table_rows"] == 1
    json.dumps(result)


@needs_nanoplot
def test_real_nanoplot_on_fastq(fastq: Path, tmp_path: Path) -> None:
    result = npl.nanoplot_plots(
        [str(fastq)], str(tmp_path / "np"), static_images=False, raw_table=True,
        plot_names=["read_length_histogram", "yield_by_minimal_length", "reads_over_time"], verbose=False,
    )  # fmt: skip
    assert result["version"] and result["statistics"]["number_of_reads"] == 40
    assert len(result["plots"]["read_length_histogram"]) >= 2 and result["missing_plots"] == [
        "reads_over_time"
    ]
    with gzip.open(result["raw_table_file"], "rt") as handle:
        assert len(handle.read().splitlines()) == 41
    assert result["statistics"]["mean_qual"] == pytest.approx(20.0, abs=0.1)  # '5' = Q20


@needs_nanoplot
def test_real_nanoplot_on_summary_gives_time_plots(tmp_path: Path) -> None:
    header = "channel\tstart_time\tduration\tsequence_length_template\tmean_qscore_template\tpasses_filtering\n"
    rows = "".join(
        f"{1 + i % 50}\t{i * 30.0}\t{5 + i % 7}\t{500 + 10 * i}\t{8 + i % 5}\tTrue\n"
        for i in range(200)
    )
    summary = tmp_path / "sequencing_summary.txt"
    summary.write_text(header + rows)
    result = npl.nanoplot_plots(
        [str(summary)],
        str(tmp_path / "np"),
        input_type="summary",
        static_images=False,
        verbose=False,
    )
    for name in ("reads_over_time", "cumulative_yield_gigabases", "reads_per_channel_activity_map"):
        assert name in result["plots"], name
