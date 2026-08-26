"""Contract tests for the explicit, non-redistributed Dashing adapter."""

from __future__ import annotations

import json
import sys
import time
from collections.abc import Iterator
from pathlib import Path

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet, ExecutionMode, Gap, Topology
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.similarity import DashingAdapter


def _fake_dashing(tmp_path: Path) -> Path:
    executable = tmp_path / "fake-dashing"
    executable.write_text(
        f"""#!{sys.executable}
import json
import os
import subprocess
import sys
import time
from pathlib import Path

arguments = sys.argv[1:]
Path(__file__).with_suffix('.args.json').write_text(json.dumps(arguments), encoding='utf-8')
behavior = os.environ.get('FAKE_DASHING_BEHAVIOR', 'success')
if behavior == 'timeout':
    subprocess.Popen(['sleep', '3'])
    time.sleep(3)
if behavior == 'noisy':
    sys.stdout.write('x' * 10000)
    sys.stdout.flush()
    time.sleep(3)
if behavior != 'no-version':
    print('Dashing version: 1.2.3-test', file=sys.stderr, flush=True)
if behavior == 'nonzero':
    raise SystemExit(7)

def option(flag):
    index = arguments.index(flag)
    return Path(arguments[index + 1])

matrix_path = option('-O')
sizes_path = option('-o')
input_paths = option('-F').read_text(encoding='utf-8').splitlines()
if behavior == 'mutate':
    Path(input_paths[0]).write_text('>changed\\nAAAA\\n', encoding='ascii')
if behavior == 'huge-artifact':
    matrix_path.write_bytes(b'x' * 100000)
    sizes_path.write_text('sizes', encoding='utf-8')
    raise SystemExit(0)
if behavior == 'malformed':
    matrix_path.write_text('not-a-dashing-matrix\\n', encoding='utf-8')
    sizes_path.write_text('sizes', encoding='utf-8')
    raise SystemExit(0)

rows = ['##Names\\t' + '\\t'.join(input_paths)]
for row_index, input_path in enumerate(input_paths):
    fields = [input_path]
    for column_index in range(len(input_paths)):
        if column_index <= row_index:
            fields.append('-')
        else:
            fields.append(str(1.0 - abs(row_index - column_index) / 10.0))
    rows.append('\\t'.join(fields))
matrix_path.write_text('\\n'.join(rows) + '\\n', encoding='utf-8')
sizes_path.write_text('#Path\\tSize (est.)\\n', encoding='utf-8')
""",
        encoding="utf-8",
    )
    executable.chmod(0o700)
    return executable


def _arguments(executable: Path) -> list[str]:
    payload = executable.with_suffix(".args.json").read_text(encoding="utf-8")
    value = json.loads(payload)
    assert isinstance(value, list) and all(isinstance(item, str) for item in value)
    return value


def test_dashing_exact_matrix_uses_fixed_command_and_records_provenance(tmp_path: Path) -> None:
    executable = _fake_dashing(tmp_path)
    output_path = tmp_path / "raw-jaccard.tsv"

    result = DashingAdapter(executable).matrix(
        (DNASequence("AACCGG"), DNASequence("AACCTT")),
        k=2,
        mode="exact",
        canonical=False,
        threads=2,
        output_path=output_path,
    )

    assert result.method == "exact-kmer-set"
    assert result.labels == ("sequence_1", "sequence_2")
    assert result.values == ((1.0, 0.9), (0.9, 1.0))
    assert result.raw_output_path == str(output_path)
    assert output_path.is_file()
    assert result.provenance.backend is not None
    assert result.provenance.backend.version == "1.2.3-test"
    assert result.provenance.backend.license_expression == "GPL-3.0-only"
    assert result.provenance.backend.metadata["redistributed"] is False
    assert result.provenance.implementation.execution_mode is ExecutionMode.EXTERNAL
    assert len(result.input_sha256) == 2
    json.dumps(result.to_dict(), sort_keys=True)

    arguments = _arguments(executable)
    assert arguments[0] == "dist"
    assert "--avoid-sorting" in arguments
    assert "--use-full-khash-sets" in arguments
    assert "-C" in arguments
    assert "--nearest-neighbors" not in arguments
    assert "-S" not in arguments
    assert set(item for item in arguments if item.startswith("--")) <= {
        "--avoid-sorting",
        "--use-full-khash-sets",
    }
    temporary_matrix = Path(arguments[arguments.index("-O") + 1])
    assert not temporary_matrix.exists()


def test_dashing_sketch_top_k_is_deterministic_hybrid_postprocessing(tmp_path: Path) -> None:
    executable = _fake_dashing(tmp_path)
    records = DNASet(
        (
            DNARecord(DNASequence("AACCGG"), "a"),
            DNARecord(DNASequence("AACCTT"), "b"),
            DNARecord(DNASequence("AATTTT"), "c"),
        )
    )

    result = DashingAdapter(executable).top_k(
        records,
        top_k=1,
        k=2,
        sketch_size_log2=8,
    )

    assert tuple(row.query_label for row in result.rows) == ("a", "b", "c")
    assert result.rows[0].hits[0].index == 1
    assert result.rows[1].hits[0].index == 0
    assert result.rows[2].hits[0].index == 1
    assert result.parameters["top_k_selection"] == "dnakit-deterministic-score-index-order"
    assert result.provenance.implementation.execution_mode is ExecutionMode.HYBRID
    arguments = _arguments(executable)
    assert arguments[arguments.index("-S") + 1] == "8"
    assert "--nearest-neighbors" not in arguments


def test_dashing_accepts_explicit_paths_and_rejects_input_mutation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    executable = _fake_dashing(tmp_path)
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">first\nAACCGG\n", encoding="ascii")
    second.write_text(">second\nAACCTT\n", encoding="ascii")
    monkeypatch.setenv("FAKE_DASHING_BEHAVIOR", "mutate")

    with pytest.raises(BackendExecutionError) as error:
        DashingAdapter(executable).matrix((first, second), k=2)

    assert error.value.code == "DASHING_INPUT_CHANGED"


@pytest.mark.parametrize(
    ("behavior", "expected_code"),
    [
        ("nonzero", "DASHING_EXECUTION_FAILED"),
        ("malformed", "DASHING_OUTPUT_PROTOCOL_ERROR"),
        ("huge-artifact", "BACKEND_ARTIFACT_OUTPUT_LIMIT"),
        ("noisy", "BACKEND_COMMAND_OUTPUT_LIMIT"),
        ("no-version", "DASHING_VERSION_UNPARSEABLE"),
    ],
)
def test_dashing_external_failures_are_structured_and_bounded(
    behavior: str,
    expected_code: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    executable = _fake_dashing(tmp_path)
    monkeypatch.setenv("FAKE_DASHING_BEHAVIOR", behavior)

    with pytest.raises(BackendExecutionError) as error:
        DashingAdapter(executable).matrix(
            (DNASequence("AACCGG"), DNASequence("AACCTT")),
            k=2,
            max_capture_bytes=100 if behavior == "noisy" else 1_000_000,
            max_output_bytes=100 if behavior == "huge-artifact" else 1_000_000,
        )

    assert error.value.code == expected_code


def test_dashing_timeout_terminates_process_tree(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    executable = _fake_dashing(tmp_path)
    monkeypatch.setenv("FAKE_DASHING_BEHAVIOR", "timeout")
    started = time.monotonic()

    with pytest.raises(BackendTimeoutError):
        DashingAdapter(executable).matrix(
            (DNASequence("AACCGG"), DNASequence("AACCTT")),
            k=2,
            timeout_seconds=0.05,
        )

    assert time.monotonic() - started < 1.5


def test_dashing_input_and_configuration_boundaries_fail_before_execution(
    tmp_path: Path,
) -> None:
    executable = _fake_dashing(tmp_path)
    adapter = DashingAdapter(executable)

    def endless() -> Iterator[DNASequence]:
        while True:
            yield DNASequence("AACCGG")

    with pytest.raises(ConfigurationError) as count_error:
        adapter.matrix(endless(), k=2, max_items=2)
    assert count_error.value.code == "DASHING_ITEM_LIMIT_EXCEEDED"
    with pytest.raises(ConfigurationError) as memory_error:
        adapter.matrix(
            (DNASequence("AACCGG"), DNASequence("AACCTT")),
            k=2,
            sketch_size_log2=20,
            max_sketch_memory_bytes=1_000_000,
        )
    assert memory_error.value.code == "DASHING_SKETCH_MEMORY_LIMIT_EXCEEDED"
    with pytest.raises(ConfigurationError) as gap_error:
        adapter.matrix((DNASequence(["AA", Gap(2), "CC"]), DNASequence("AACCGG")), k=2)
    assert gap_error.value.code == "DASHING_GAPPED_INPUT_UNSUPPORTED"
    with pytest.raises(ConfigurationError) as circle_error:
        adapter.matrix(
            (
                DNASequence("AACCGG", topology=Topology.CIRCULAR),
                DNASequence("AACCTT"),
            ),
            k=2,
        )
    assert circle_error.value.code == "DASHING_CIRCULAR_INPUT_UNSUPPORTED"
    with pytest.raises(ConfigurationError):
        adapter.top_k(
            (DNASequence("AACCGG"), DNASequence("AACCTT")),
            top_k=2,
            k=2,
        )
    assert not executable.with_suffix(".args.json").exists()


def test_dashing_requires_an_explicit_executable_and_safe_output_path(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        DashingAdapter(tmp_path / "missing-dashing")

    executable = _fake_dashing(tmp_path)
    input_path = tmp_path / "input.fasta"
    other_path = tmp_path / "other.fasta"
    input_path.write_text(">input\nAACCGG\n", encoding="ascii")
    other_path.write_text(">other\nAACCTT\n", encoding="ascii")
    with pytest.raises(ConfigurationError) as conflict:
        DashingAdapter(executable).matrix(
            (input_path, other_path),
            k=2,
            output_path=input_path,
            overwrite=True,
        )
    assert conflict.value.code == "DASHING_OUTPUT_INPUT_CONFLICT"
