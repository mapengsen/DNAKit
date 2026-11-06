import random
import shutil

import pytest

from dnakit.alignment import mafft_p13_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")


@pytest.fixture()
def seqs(tmp_path):
    rng = random.Random(8)
    base = "".join(rng.choice("ACGT") for _ in range(160))
    records = []
    for i in range(8):
        s = list(base)
        for _ in range(8):
            s[rng.randrange(len(s))] = rng.choice("ACGT")
        if i % 2:
            del s[rng.randrange(20, 140)]
        records.append((f"s{i}", "".join(s)))
    (tmp_path / "all.fa").write_text("".join(f">{n}\n{s}\n" for n, s in records))
    (tmp_path / "first.fa").write_text("".join(f">{n}\n{s}\n" for n, s in records[:5]))
    (tmp_path / "rest.fa").write_text("".join(f">{n}\n{s}\n" for n, s in records[5:]))
    (tmp_path / "frag.fa").write_text(f">frag\n{records[6][1][30:90]}\n")
    return tmp_path


def read_fasta(path):
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        else:
            out[name] += line
    return out


def test_validation(seqs):
    with pytest.raises(ConfigurationError):
        m.mafft_align(str(seqs / "all.fa"), str(seqs / "o.fa"), strategy="bogus", verbose=False)
    with pytest.raises(BackendUnavailableError):
        m.mafft_align(
            str(seqs / "all.fa"), str(seqs / "o.fa"), mafft_path=str(seqs / "none"), verbose=False
        )


@needs
@pytest.mark.parametrize("strategy", ["fftnsi", "nwns", "nwnsi", "linsi", "ginsi", "einsi"])
def test_strategies_keep_residues(seqs, strategy):
    if shutil.which(f"mafft-{strategy}") is None:
        pytest.skip(f"mafft-{strategy} not installed")
    out = seqs / f"{strategy}.fa"
    result = m.mafft_align(str(seqs / "all.fa"), str(out), strategy=strategy, verbose=False)
    assert result["status"] == "completed" and result["sequence_count"] == 8
    aligned, original = read_fasta(out), read_fasta(seqs / "all.fa")
    assert list(aligned) == list(original)
    assert len({len(v) for v in aligned.values()}) == 1
    assert {k: v.replace("-", "").upper() for k, v in aligned.items()} == original


@needs
def test_add_seed_profile_partition_distance(seqs):
    base = seqs / "base.fa"
    m.mafft_align(str(seqs / "first.fa"), str(base), strategy="fftnsi", verbose=False)
    width = len(next(iter(read_fasta(base).values())))
    added = m.mafft_add_sequences(
        str(base), str(seqs / "rest.fa"), str(seqs / "add.fa"), verbose=False
    )
    assert added["sequence_count"] == 8
    frag = m.mafft_add_sequences(
        str(base),
        str(seqs / "frag.fa"),
        str(seqs / "fr.fa"),
        mode="addfragments",
        keep_length=True,
        verbose=False,
    )
    assert len(read_fasta(seqs / "fr.fa")["frag"]) == width and frag["mode"] == "addfragments"
    assert (
        m.mafft_seed_alignment(
            str(seqs / "rest.fa"), [str(base)], str(seqs / "seed.fa"), verbose=False
        )["sequence_count"]
        == 8
    )
    other = seqs / "other.fa"
    m.mafft_align(str(seqs / "rest.fa"), str(other), strategy="fftnsi", verbose=False)
    if shutil.which("mafft-profile"):
        assert (
            m.mafft_profile_alignment(str(base), str(other), str(seqs / "prof.fa"), verbose=False)[
                "sequence_count"
            ]
            == 8
        )
    pt = m.mafft_partition_tree_align(
        str(seqs / "all.fa"), str(seqs / "pt.fa"), partition_size=4, verbose=False
    )
    assert pt["sequence_count"] == 8
    if shutil.which("mafft-distance"):
        d = m.mafft_pairwise_distances(str(seqs / "all.fa"), str(seqs / "d.txt"), verbose=False)
        assert d["pair_count"] == 28 and d["distances"][0]["length_i"] > 100
    if shutil.which("ruby") and shutil.which("mafft-sparsecore.rb"):
        sc = m.mafft_sparsecore_align(
            str(seqs / "all.fa"),
            str(seqs / "sc.fa"),
            core_size=5,
            core_candidates="100%",
            verbose=False,
        )
        assert sc["sequence_count"] == 8
