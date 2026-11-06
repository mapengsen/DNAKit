"""Tests for dnakit.alignment.pyhmmer_p39_functions (HMMER pipelines on DNA, Easel alignment tools, SSI index)."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

pyhmmer = pytest.importorskip("pyhmmer")

from dnakit.alignment import pyhmmer_p39_functions as ph
from dnakit.exceptions import ConfigurationError

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("ph_p39")
    rng = random.Random(7)
    base = "".join(rng.choice("ACGT") for _ in range(120))
    rows = {}
    for index in range(6):
        row = list(base)
        for _ in range(4 + index):
            row[rng.randrange(len(row))] = rng.choice("ACGT")
        rows[f"s{index}"] = "".join(row)
    rows["s6"] = "".join(rng.choice("ACGT") for _ in range(120))
    alignment = root / "aln.afa"
    alignment.write_text("".join(f">{name}\n{row}\n" for name, row in rows.items()))
    sequences = root / "seqs.fa"
    sequences.write_text("".join(f">{name}\n{row}\n" for name, row in rows.items()))
    alphabet = pyhmmer.easel.Alphabet.dna()
    with pyhmmer.easel.MSAFile(str(alignment), digital=True, alphabet=alphabet) as handle:
        msa = handle.read()
    msa.name = b"fam"
    hmm, _profile, _optimized = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    profile = root / "fam.hmm"
    with profile.open("wb") as out:
        hmm.write(out)
    return {"root": root, "aln": alignment, "seqs": sequences, "hmm": profile}


def test_hmmsearch_and_hmmscan_find_family_members(data: dict[str, Path]) -> None:
    search = ph.pyhmmer_hmmsearch(str(data["hmm"]), str(data["seqs"]), verbose=False)
    hits = search["results"][0]["hits"]
    assert {hit["name"] for hit in hits} >= {"s0", "s1", "s2", "s3", "s4", "s5"}
    assert all(hit["included"] for hit in hits)
    assert hits[0]["domains"][0]["env_from"] == 1
    scan = ph.pyhmmer_hmmscan(str(data["seqs"]), str(data["hmm"]), verbose=False)
    by_query = {item["query"]: item["hit_count"] for item in scan["results"]}
    assert by_query["s0"] == 1 and set(by_query.values()) == {1}


def test_hmmsearch_matches_pyhmmer_direct(data: dict[str, Path]) -> None:
    alphabet = pyhmmer.easel.Alphabet.dna()
    with pyhmmer.plan7.HMMFile(str(data["hmm"])) as handle:
        hmm = handle.read()
    with pyhmmer.easel.SequenceFile(str(data["seqs"]), digital=True, alphabet=alphabet) as handle:
        block = handle.read_block()
    direct = next(iter(pyhmmer.hmmsearch([hmm], block)))
    mine = ph.pyhmmer_hmmsearch(str(data["hmm"]), str(data["seqs"]), verbose=False)["results"][0][
        "hits"
    ]
    assert [hit["name"] for hit in mine] == [str(hit.name) for hit in direct]
    assert [hit["score"] for hit in mine] == pytest.approx([hit.score for hit in direct])


def test_hmmsearch_rejects_missing_file(data: dict[str, Path]) -> None:
    with pytest.raises(ConfigurationError):
        ph.pyhmmer_hmmsearch(str(data["root"] / "none.hmm"), str(data["seqs"]), verbose=False)


def test_build_hmm_from_sequence_and_emit(data: dict[str, Path]) -> None:
    root = data["root"]
    built = ph.build_hmm_from_sequence(
        "ACGTACGTAGCTAGCTAGGATCGATCGATGCTAGC", str(root / "one.hmm"), name="one", verbose=False
    )
    assert built["model_length"] == 35
    emitted = ph.emit_hmm_sequences(
        str(root / "one.hmm"), 4, str(root / "emit.fa"), seed=3, verbose=False
    )
    assert emitted["sequences"] == 4 and (root / "emit.fa").read_text().count(">") == 4
    again = ph.emit_hmm_sequences(
        str(root / "one.hmm"), 4, str(root / "emit2.fa"), seed=3, verbose=False
    )
    assert (root / "emit.fa").read_text() == (root / "emit2.fa").read_text() and again[
        "lengths"
    ] == emitted["lengths"]
    profile_based = ph.emit_hmm_sequences(
        str(root / "one.hmm"), 2, str(root / "emit3.fa"), source="profile", verbose=False
    )
    assert profile_based["source"] == "profile"
    alignment = ph.emit_hmm_alignment(
        str(root / "one.hmm"), 3, str(root / "emit.afa"), verbose=False
    )
    assert alignment["sequences"] == 3
    with pytest.raises(ConfigurationError):
        ph.build_hmm_from_sequence("ACGT", str(root / "one.hmm"), verbose=False)


def test_reverse_complement_alignment(data: dict[str, Path]) -> None:
    result = ph.reverse_complement_alignment(
        str(data["aln"]), str(data["root"] / "rc.afa"), verbose=False
    )
    original = {
        line.split()[0][1:]: nxt
        for line, nxt in zip(
            data["aln"].read_text().splitlines()[::2],
            data["aln"].read_text().splitlines()[1::2],
            strict=True,
        )
    }
    for name, row in original.items():
        assert result["rows"][name] == row.translate(COMPLEMENT)[::-1]


def test_weights_and_identity_filter(data: dict[str, Path]) -> None:
    weights = ph.compute_alignment_weights(str(data["aln"]), method="gsc", verbose=False)
    assert len(weights["weights"]) == 7 and weights["weight_sum"] == pytest.approx(7.0)
    pb = ph.compute_alignment_weights(str(data["aln"]), verbose=False)
    assert pb["method"] == "pb"
    filtered = ph.filter_alignment_by_identity(
        str(data["aln"]), str(data["root"] / "f.afa"), max_identity=0.9, verbose=False
    )
    assert filtered["kept_sequences"] < filtered["input_sequences"]
    assert "s6" in filtered["kept_names"]
    with pytest.raises(ConfigurationError):
        ph.compute_alignment_weights(str(data["aln"]), method="nope", verbose=False)


def test_guess_alphabet(data: dict[str, Path]) -> None:
    assert ph.guess_file_alphabet(str(data["seqs"]), verbose=False)["alphabet"] == "DNA"
    assert (
        ph.guess_file_alphabet(str(data["aln"]), kind="alignment", verbose=False)["alphabet"]
        == "DNA"
    )
    protein = data["root"] / "p.fa"
    protein.write_text(">p\nMKVLAAGIVALLLAAGCSSEEQPEKKDEQ\n")
    assert ph.guess_file_alphabet(str(protein), verbose=False)["alphabet"] == "amino"


def test_ssi_index_round_trip(data: dict[str, Path]) -> None:
    copy = data["root"] / "ssi_copy.fa"
    copy.write_text(data["seqs"].read_text())
    built = ph.build_sequence_ssi_index(str(copy), aliases={"seven": "s6"}, verbose=False)
    assert built["keys"] == 7 and built["aliases"] == 1
    fetched = ph.fetch_sequence_by_ssi(str(copy), "s3", verbose=False)
    assert (
        fetched["sequence"]
        == dict(
            zip(
                [line[1:] for line in copy.read_text().splitlines()[::2]],
                copy.read_text().splitlines()[1::2],
                strict=True,
            )
        )["s3"]
    )
    assert ph.fetch_sequence_by_ssi(str(copy), "seven", verbose=False)["name"] == "s6"
    with pytest.raises(ConfigurationError):
        ph.fetch_sequence_by_ssi(str(copy), "absent", verbose=False)
    with pytest.raises(ConfigurationError):
        ph.build_sequence_ssi_index(str(copy), verbose=False)
