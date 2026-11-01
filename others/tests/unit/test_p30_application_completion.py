"""Operation routing, output validation and direct-upstream checks for P30.

Fake executables exercise the real Biotite parser; they do not validate ViennaRNA
or tantan algorithms. Optional direct executable comparisons skip when unavailable.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from dnakit.download import sra_p30_cli_functions as sra
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError
from dnakit.secondary_structure import viennarna_p30_functions as rna


def _program(tmp_path: Path, body: str) -> str:
    import sys
    executable = tmp_path / "fake_program"
    executable.write_text(f"#!{sys.executable}\n" + body)
    executable.chmod(0o755)
    return str(executable)


@pytest.mark.parametrize("value", ["", "--help", "SRR1;touch bad", "../../file"])
def test_invalid_accession_cannot_become_command(value, tmp_path):
    with pytest.raises(ConfigurationError):
        sra.prefetch_sra_run(value, tmp_path / "out", verbose=False)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("kind", ["fastq", "fasta"])
def test_sra_retains_both_paired_output_files_and_local_input(tmp_path, kind):
    local = tmp_path / "reads.sra"
    local.write_bytes(b"SRA")
    executable = _program(tmp_path, '''import sys
from pathlib import Path
if '--version' in sys.argv:
    print('fasterq-dump 3.2.1'); raise SystemExit
p = Path(sys.argv[sys.argv.index('-o') + 1])
assert Path(sys.argv[-1]).read_bytes() == b'SRA'
for i in (1, 2):
    p.with_name(p.stem + '_' + str(i) + p.suffix).write_text('>read\\nACGT\\n')
''')
    result = sra.extract_sra_reads("SRR123", tmp_path / "out", output_format=kind,
                                   local_sra=local, executable_path=executable, verbose=False)
    assert len(result["outputs"]) == 2
    assert result["network_operation"] is False
    assert all(Path(p).is_file() for p in result["outputs"])
    json.dumps(result, allow_nan=False)
    with pytest.raises(ConfigurationError):
        sra.extract_sra_reads("SRR123", tmp_path / "out", executable_path=executable, verbose=False)


def test_sra_zero_exit_without_outputs_is_not_success(tmp_path):
    executable = _program(tmp_path, "print('3.2.1')\n")
    with pytest.raises(DNAKitError):
        sra.extract_sra_reads("SRR123", tmp_path / "out", executable_path=executable, verbose=False)


@pytest.mark.parametrize("sequence", ["", "A\nC", "-AC", "AC;GT"])
def test_rna_rejects_invalid_sequences_before_starting_backend(sequence):
    with pytest.raises(ConfigurationError):
        rna.fold_rna(sequence, verbose=False)


@pytest.mark.parametrize("sequence", ["ACGU", "NNNNNN", "ACGTACGT"])
def test_real_biotite_parser_with_controlled_rnafold_output(tmp_path, sequence):
    pytest.importorskip("biotite")
    executable = _program(tmp_path, '''import sys
if '--version' in sys.argv:
    print('RNAfold 2.7.2'); raise SystemExit
s = sys.stdin.readline().strip()
print(s)
print('.' * len(s) + ' (0.00)')
''')
    result = rna.fold_rna(sequence, executable_path=executable, verbose=False)
    assert result["structures"][0]["dot_bracket"] == "." * len(sequence)
    assert result["structures"][0]["free_energy"] == 0
    assert result["structures"][0]["base_pairs"] == []
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("sequence", ["GGGAAACCC", "ACGUACGU", "AAAAAA"])
def test_installed_rnafold_matches_direct_upstream(sequence):
    pytest.importorskip("biotite")
    executable = shutil.which("RNAfold")
    if executable is None:
        pytest.skip("ViennaRNA executable not installed; algorithm not verified")
    import biotite.sequence as seq
    from biotite.application.viennarna import FoldApp

    expected = FoldApp(path=executable).run(seq.NucleotideSequence(sequence.replace("U", "T"))).result(timeout=30)
    actual = rna.fold_rna(sequence, executable_path=executable, verbose=False)["structures"][0]
    assert actual["dot_bracket"] == expected.dot_bracket
    assert actual["free_energy"] == pytest.approx(expected.free_energy)


def test_biotite_application_dispatch_contract(monkeypatch):
    """Every fixed adapter passes the appropriate input shape and upstream default."""
    calls = []

    class Structure:
        dot_bracket = "...."
        free_energy = 0
        def base_pairs(self):
            class Pairs:
                def tolist(self): return []
            return Pairs()

    class Future:
        def result(self, timeout): return Structure()

    class Application:
        def __init__(self, path): pass
        def run(self, *inputs, **options):
            calls.append((inputs, options))
            return Future()

    class Seq:
        NucleotideSequence = str

    class Align:
        class Alignment:
            @staticmethod
            def from_strings(rows, factory): return ("alignment", list(rows))

    class Apps:
        pass

    for program, name in rna._APPS.values():
        setattr(Apps, name, Application)
    monkeypatch.setattr(rna.cli, "program", lambda *a, **k: "/fixed/tool")
    monkeypatch.setattr(rna.cli, "_version_of", lambda p: "2.7.2")
    monkeypatch.setattr(rna, "require", lambda name, extra: Apps if "viennarna" in name else Align if name.endswith("align") else Seq)
    for name in ("fold_rna", "fold_rna_pseudoknots", "fold_rna_local"):
        getattr(rna, name)("ACGU", verbose=False)
        assert calls[-1][0] == ("ACGT",)
    for name in ("cofold_rna", "duplex_rna", "find_rna_interactions"):
        getattr(rna, name)("ACGU", "UGCA", verbose=False)
        assert calls[-1][0] == ("ACGT", "TGCA")
    rna.fold_rna_complex(["ACGU", "UGCA"], verbose=False)
    assert calls[-1][0] == (["ACGT", "TGCA"],)
    for name in ("fold_rna_alignment", "fold_rna_alignment_local"):
        getattr(rna, name)(["ACGU", "AC-U"], verbose=False)
        assert calls[-1][0] == (("alignment", ["ACGT", "AC-T"]),)
    rna.duplex_rna_alignments(["ACGU", "ACGU"], ["UGCA", "UGCA"], verbose=False)
    assert len(calls[-1][0]) == 2
    assert all(options["temperature"] == 37 for _, options in calls)


def test_missing_executable_has_install_hint():
    with pytest.raises(BackendUnavailableError):
        rna.fold_rna("ACGT", executable_path="/does/not/exist/RNAfold", verbose=False)


@pytest.mark.parametrize('coordinates, succeeds', [([[1.0, 2.0]], True), ([[float('nan'), 2.0]], False)])
def test_single_base_plot_coordinates_are_rows_and_json_finite(monkeypatch, coordinates, succeeds):
    np = pytest.importorskip('numpy')
    class Future:
        def result(self, timeout): return np.array(coordinates[0])
    class Plot:
        class Layout:
            NAVIEW = 1
            @classmethod
            def __class_getitem__(cls, key): return getattr(cls, key)
        def __init__(self, path): pass
        def run(self, *args, **kwargs): return Future()
    class Apps:
        PlotApp = Plot
    class Structure:
        @staticmethod
        def base_pairs_from_dot_bracket(value): return []
    monkeypatch.setattr(rna.cli, 'program', lambda *a, **k: '/fixed/tool')
    monkeypatch.setattr(rna.cli, '_version_of', lambda p: '2.7.2')
    monkeypatch.setattr(rna, 'require', lambda name, extra: Apps if 'viennarna' in name else Structure)
    if succeeds:
        result = rna.rna_structure_coordinates('.', verbose=False)
        assert result['coordinates'] == coordinates
        json.dumps(result, allow_nan=False)
    else:
        with pytest.raises(DNAKitError, match='non-finite'):
            rna.rna_structure_coordinates('.', verbose=False)
