import random
import shutil
import subprocess

import pytest

from dnakit.alignment import bwa_p13_cli_functions as b
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs = pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not installed")
COMP = str.maketrans("ACGT", "TGCA")


@pytest.fixture()
def data(tmp_path):
    rng = random.Random(2)
    genome = "".join(rng.choice("ACGT") for _ in range(30000))
    (tmp_path / "g.fa").write_text(">chr1\n" + genome + "\n")
    with (tmp_path / "p1.fq").open("w") as f1, (tmp_path / "p2.fq").open("w") as f2:
        for i in range(30):
            start = rng.randrange(0, 29000)
            frag = genome[start : start + rng.randrange(120, 200)]
            a, c = frag[:100], frag[-100:].translate(COMP)[::-1]
            f1.write(f"@p{i}/1\n{a}\n+\n{'I' * 100}\n")
            f2.write(f"@p{i}/2\n{c}\n+\n{'I' * 100}\n")
    (tmp_path / "q.fa").write_text(">q1\n" + genome[500:600] + "\n>q2\nACGTACGTA\n")
    if shutil.which("bwa"):
        subprocess.run(["bwa", "index", str(tmp_path / "g.fa")], capture_output=True, check=True)
    return tmp_path


def test_validation(tmp_path):
    (tmp_path / "x.fa").write_text(">a\nACGT\n")
    with pytest.raises(ConfigurationError):
        b.bwa_index_stage("nonsense", str(tmp_path / "x.fa"), str(tmp_path / "o"), verbose=False)
    with pytest.raises(BackendUnavailableError):
        b.bwa_index_stage(
            "fa2pac",
            str(tmp_path / "x.fa"),
            str(tmp_path / "o"),
            bwa_path=str(tmp_path / "none"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        b.bwa_sai_to_sam(
            str(tmp_path / "g"),
            [str(tmp_path / "x.fa")],
            [str(tmp_path / "x.fa")],
            str(tmp_path / "o"),
            verbose=False,
        )
    with pytest.raises(BackendUnavailableError):
        b.bwakit_generate_reference(
            "hs38", str(tmp_path / "ref"), bwakit_directory=str(tmp_path), verbose=False
        )


@needs
def test_index_stages_reproduce_bwa_index(data):
    t = data
    b.bwa_index_stage("fa2pac", str(t / "g.fa"), str(t / "f"), forward_only=True, verbose=False)
    for ext in ("pac", "ann", "amb"):
        assert (t / f"f.{ext}").read_bytes() == (t / f"g.fa.{ext}").read_bytes()
    b.bwa_index_stage(
        "fa2pac", str(t / "g.fa"), str(t / "p"), verbose=False
    )  # the BWT is built from the two-strand pac
    b.bwa_index_stage("pac2bwt", str(t / "p.pac"), str(t / "p.bwt"), verbose=False)
    b.bwa_index_stage("pac2bwtgen", str(t / "p.pac"), str(t / "q.bwt"), verbose=False)
    b.bwa_index_stage("bwtupdate", str(t / "p.bwt"), str(t / "pu.bwt"), verbose=False)
    b.bwa_index_stage("bwt2sa", str(t / "pu.bwt"), str(t / "p.sa"), verbose=False)
    assert (t / "pu.bwt").read_bytes() == (t / "g.fa.bwt").read_bytes()
    assert (t / "p.sa").read_bytes() == (t / "g.fa.sa").read_bytes()


@needs
def test_aln_sampe_samse_fastmap_pemerge_maxk(data):
    t = data
    prefix = str(t / "g.fa")
    for n in (1, 2):
        assert (
            b.bwa_aln(prefix, str(t / f"p{n}.fq"), str(t / f"{n}.sai"), verbose=False)["status"]
            == "completed"
        )
    pe = b.bwa_sai_to_sam(
        prefix,
        [str(t / "1.sai"), str(t / "2.sai")],
        [str(t / "p1.fq"), str(t / "p2.fq")],
        str(t / "pe.sam"),
        read_group="@RG\\tID:x\\tSM:y",
        verbose=False,
    )
    assert pe["paired"] and pe["alignment_count"] == 60
    se = b.bwa_sai_to_sam(
        prefix, [str(t / "1.sai")], [str(t / "p1.fq")], str(t / "se.sam"), verbose=False
    )
    assert not se["paired"] and se["alignment_count"] == 30
    with pytest.raises(ConfigurationError):
        b.bwa_sai_to_sam(
            prefix,
            [str(t / "1.sai")],
            [str(t / "p1.fq")],
            str(t / "bad.sam"),
            max_insert_size=300,
            verbose=False,
        )
    fm = b.bwa_fastmap(prefix, str(t / "q.fa"), str(t / "fm.txt"), verbose=False)
    assert fm["queries"][0]["matches"][0]["end"] == 100 and fm["queries"][1]["matches"] == []
    pm = b.bwa_pemerge(
        str(t / "p1.fq"),
        str(t / "pm.fq"),
        read2_path=str(t / "p2.fq"),
        output="merged",
        verbose=False,
    )
    assert pm["record_count"] > 0
    mk = b.bwa_maxk(prefix + ".bwt", str(t / "q.fa"), str(t / "mk.txt"), verbose=False)
    assert {"length": 100, "positions": 100} in mk["histogram"]


@pytest.mark.skipif(
    shutil.which("qualfa2fq.pl") is None or shutil.which("perl") is None,
    reason="qualfa2fq.pl not installed",
)
def test_qual_fasta_to_fastq(tmp_path):
    (tmp_path / "t.fa").write_text(">s1\nACGTAC\n>s2\nGGGTT\n")
    (tmp_path / "t.qual").write_text(">s1\n10 20 30 40 50 60\n>s2\n5 6 7 8 9\n")
    result = b.bwa_qual_fasta_to_fastq(
        str(tmp_path / "t.fa"), str(tmp_path / "t.qual"), str(tmp_path / "o.fq"), verbose=False
    )
    assert result["record_count"] == 2
    assert (tmp_path / "o.fq").read_text().splitlines()[:4] == ["@s1", "ACGTAC", "+", "+5?IS]"]
