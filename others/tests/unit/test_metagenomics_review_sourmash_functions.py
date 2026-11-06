"""Tests for dnakit.metagenomics_review_sourmash_functions (sourmash sketch fromfile, sig grep/inflate, SBT tools)."""

from __future__ import annotations

import inspect
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit import metagenomics_review_sourmash_functions as review
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_SOURMASH = Path(sys.executable).parent / "sourmash"
pytestmark = pytest.mark.skipif(
    not _SOURMASH.is_file() and shutil.which("sourmash") is None, reason="sourmash is not installed"
)


def _sm(*args: str) -> None:
    subprocess.run([str(_SOURMASH), *args], capture_output=True, text=True, timeout=300, check=True)


def _names(path: Path) -> list[str]:
    import sourmash

    return sorted(sig.name for sig in sourmash.load_file_as_signatures(str(path)))


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("sourmash_review_functions")
    for number in (1, 2, 3):
        rng = random.Random(number)
        seq = "".join(rng.choice("ACGT") for _ in range(30000))
        (root / f"g{number}.fa").write_text(f">g{number}\n{seq}\n")
    header = "name,genome_filename,protein_filename\n"
    rows = {number: f"g{number},{root / f'g{number}.fa'},\n" for number in (1, 2, 3)}
    (root / "all.csv").write_text(header + "".join(rows[i] for i in (1, 2, 3)))
    (root / "ab.csv").write_text(header + rows[1] + rows[2])
    (root / "c.csv").write_text(header + rows[3])
    _sm(
        "sketch",
        "fromfile",
        str(root / "all.csv"),
        "--output-signatures",
        str(root / "all.zip"),
        "-p",
        "dna,k=31,scaled=100",
    )
    _sm(
        "sketch",
        "fromfile",
        str(root / "ab.csv"),
        "--output-signatures",
        str(root / "ab.zip"),
        "-p",
        "dna,k=31,scaled=100",
    )
    _sm(
        "sketch",
        "fromfile",
        str(root / "c.csv"),
        "--output-signatures",
        str(root / "c.zip"),
        "-p",
        "dna,k=31,scaled=100",
    )
    # reference with abundances for sig inflate
    _sm(
        "sketch",
        "dna",
        str(root / "g1.fa"),
        "-p",
        "k=31,scaled=100,abund",
        "-o",
        str(root / "g1_abund.sig"),
    )
    _sm("index", str(root / "ab.sbt.zip"), str(root / "ab.zip"), "-k", "31", "--scaled", "100")
    _sm("index", str(root / "c.sbt.zip"), str(root / "c.zip"), "-k", "31", "--scaled", "100")
    _sm("index", str(root / "all.sbt.zip"), str(root / "all.zip"), "-k", "31", "--scaled", "100")
    return root


def test_sketch_fromfile_writes_every_listed_genome(data: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.zip"
    result = review.sketch_fromfile_signatures(
        [str(data / "all.csv")], str(out), param_strings=["dna,k=31,scaled=100"]
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _names(out) == ["g1", "g2", "g3"]


def test_sketch_fromfile_refuses_existing_output(data: Path, tmp_path: Path) -> None:
    existing = tmp_path / "exists.zip"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        review.sketch_fromfile_signatures([str(data / "all.csv")], str(existing))


def test_sketch_fromfile_needs_existing_csv(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        review.sketch_fromfile_signatures([str(tmp_path / "missing.csv")], str(tmp_path / "o.zip"))


def test_grep_keeps_matching_names_and_inverts(data: Path, tmp_path: Path) -> None:
    matched = tmp_path / "g2.zip"
    result = review.grep_signatures(
        "^g2$", [str(data / "all.zip")], str(matched), ksize=31, moltype="DNA"
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _names(matched) == ["g2"]
    inverted = tmp_path / "not_g2.zip"
    result = review.grep_signatures(
        "^g2$", [str(data / "all.zip")], str(inverted), invert_match=True
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _names(inverted) == ["g1", "g3"]


def test_grep_ignore_case_and_csv_manifest(data: Path, tmp_path: Path) -> None:
    out = tmp_path / "upper.zip"
    manifest = tmp_path / "manifest.csv"
    result = review.grep_signatures(
        "G3", [str(data / "all.zip")], str(out), ignore_case=True, csv_manifest=str(manifest)
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _names(out) == ["g3"]
    assert "g3" in manifest.read_text()


def test_grep_rejects_unknown_moltype(data: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        review.grep_signatures(
            "g", [str(data / "all.zip")], str(tmp_path / "x.zip"), moltype="codon"
        )


def test_inflate_gives_abundances_to_other_signatures(data: Path, tmp_path: Path) -> None:
    import sourmash

    out = tmp_path / "inflated.zip"
    result = review.inflate_signatures(
        str(data / "g1_abund.sig"), [str(data / "all.zip")], str(out), ksize=31
    )
    assert result["status"] == "completed", result["log_tail"]
    sigs = list(sourmash.load_file_as_signatures(str(out)))
    assert sorted(sig.name for sig in sigs) == ["g1", "g2", "g3"]
    assert all(sig.minhash.track_abundance for sig in sigs)


def test_inflate_fails_without_reference_abundances(data: Path, tmp_path: Path) -> None:
    result = review.inflate_signatures(
        str(data / "all.zip"), [str(data / "all.zip")], str(tmp_path / "bad.zip"), ksize=31
    )
    assert result["status"] == "failed"


def test_combine_sbt_indexes_contains_all_signatures(data: Path, tmp_path: Path) -> None:
    import sourmash

    out = tmp_path / "combined.sbt.zip"
    result = review.combine_sbt_indexes(
        str(out), [str(data / "ab.sbt.zip"), str(data / "c.sbt.zip")]
    )
    assert result["status"] == "completed", result["log_tail"]
    index = sourmash.load_file_as_index(str(out))
    assert sorted(sig.name for sig in index.signatures()) == ["g1", "g2", "g3"]


def test_migrate_keeps_original_and_signatures(data: Path, tmp_path: Path) -> None:
    import hashlib

    import sourmash

    source = data / "c.sbt.zip"
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    result = review.migrate_sbt_index(str(source), str(tmp_path / "migrated"))
    assert result["status"] == "completed", result["log_tail"]
    migrated = tmp_path / "migrated" / "c.sbt.zip"
    assert migrated.is_file()
    index = sourmash.load_file_as_index(str(migrated))
    assert [sig.name for sig in index.signatures()] == ["g3"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before


def test_convert_sbt_storage_to_zip_backend(data: Path, tmp_path: Path) -> None:
    import sourmash

    out = tmp_path / "converted.sbt.json"
    storage = tmp_path / "converted_data.zip"
    result = review.convert_sbt_storage(str(data / "ab.sbt.zip"), str(out), f"zip({storage})")
    assert result["status"] == "completed", result["log_tail"]
    assert storage.is_file()
    index = sourmash.load_file_as_index(str(out))
    assert sorted(sig.name for sig in index.signatures()) == ["g1", "g2"]


def test_convert_rejects_unknown_backend(data: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        review.convert_sbt_storage(str(data / "c.sbt.zip"), str(tmp_path / "x.json"), "redis(host)")


def test_missing_sourmash_path_is_reported(data: Path, tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        review.sketch_fromfile_signatures(
            [str(data / "all.csv")],
            str(tmp_path / "o.zip"),
            sourmash_path=str(tmp_path / "no_such_sourmash"),
        )


def test_public_functions_take_keyword_only_verbose() -> None:
    public = [
        review.sketch_fromfile_signatures,
        review.grep_signatures,
        review.inflate_signatures,
        review.combine_sbt_indexes,
        review.migrate_sbt_index,
        review.convert_sbt_storage,
    ]
    for function in public:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__
