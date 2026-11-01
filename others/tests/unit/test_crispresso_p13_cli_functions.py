import gzip
import random
import shutil

import pytest

from dnakit.design_crispr import crispresso_p13_cli_functions as c
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs = pytest.mark.skipif(shutil.which("CRISPResso") is None, reason="CRISPResso2 not installed")
GUIDE = "GGCACTGCGGCTGGAGGTGG"


@pytest.fixture()
def reads(tmp_path):
    rng = random.Random(11)
    amplicon = "".join(rng.choice("ACGT") for _ in range(60)) + GUIDE + "GGG" + "".join(rng.choice("ACGT") for _ in range(70))
    cut = amplicon.index(GUIDE) + 17
    with (tmp_path / "r1.fq").open("w") as handle:
        for i in range(200):
            seq = amplicon if rng.random() < 0.6 else amplicon[:cut] + amplicon[cut + 3 :]
            handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return tmp_path, amplicon


def test_validation(reads):
    path, amplicon = reads
    with pytest.raises(ConfigurationError):
        c.crispresso_quantify_editing(str(path / "r1.fq"), str(path / "o"), amplicon_sequences="ACGX", verbose=False)
    with pytest.raises(ConfigurationError):
        c.crispresso_quantify_editing(str(path / "r1.fq"), str(path / "o"), amplicon_sequences=amplicon, base_editor_output=True, verbose=False)
    with pytest.raises(BackendUnavailableError):
        c.crispresso_quantify_editing(str(path / "r1.fq"), str(path / "o"), amplicon_sequences=amplicon, crispresso_path=str(path / "none"), verbose=False)
    with pytest.raises(ConfigurationError):
        c.crispresso_pooled(str(path / "r1.fq"), str(path / "o"), verbose=False)


def test_filter_reads_matches_expected(reads, tmp_path):
    path, _ = reads
    result = c.crispresso_filter_reads_by_sequence(str(path / "r1.fq"), str(path / "kept.fq"), include_sequences=[GUIDE], verbose=False)
    total = sum(1 for line in (path / "r1.fq").read_text().splitlines()[1::4] if GUIDE in line)
    assert result["reads_written"] == total == len((path / "kept.fq").read_text().splitlines()) // 4
    gz = path / "kept.fq.gz"
    c.crispresso_filter_reads_by_sequence(str(path / "r1.fq"), str(gz), exclude_sequences=[GUIDE], verbose=False)
    with gzip.open(gz, "rt") as handle:
        assert len(handle.read().splitlines()) // 4 == 200 - total


def test_count_high_quality_bases(tmp_path):
    ref, read = "ACGTACGT", "ACGAACGT"
    with gzip.open(tmp_path / "a.fastq.gz", "wt") as handle:
        handle.write(f"@r1\n{read}\n+ ALN=Reference ALN_DETAILS=Reference,{read},{ref},100\n{'I' * 8}\n")
        handle.write(f"@r2\n{read}\n+ ALN=Other ALN_DETAILS=Other,{read},{ref},100\n{'I' * 8}\n")
    result = c.crispresso_count_high_quality_bases(str(tmp_path / "a.fastq.gz"), 3, "T", "A", verbose=False)
    assert result["counts"]["target_high"] == 1 and result["counts"]["unaligned"] == 1


@needs
def test_quantify_editing_and_run_based_helpers(reads):
    path, amplicon = reads
    result = c.crispresso_quantify_editing(
        str(path / "r1.fq"), str(path / "a"), amplicon_sequences=amplicon, guide_sequences=GUIDE, write_detailed_allele_table=True, suppress_plots=True, verbose=False
    )
    assert result["status"] == "completed" and result["mapping_statistics"]["READS IN INPUTS"] == 200
    row = result["quantification"][0]
    assert row["Reads_in_input"] == 200 and 0 < row["Modified%"] < 100
    base = c.crispresso_quantify_editing(
        str(path / "r1.fq"), str(path / "b"), amplicon_sequences=amplicon, guide_sequences=GUIDE, base_editor_output=True, conversion_nuc_from="C", conversion_nuc_to="T", suppress_plots=True, verbose=False
    )
    assert base["status"] == "completed"
    if shutil.which("python") and pytest.importorskip("CRISPResso2"):
        run = result["run_folders"][0]
        counts = c.crispresso_count_sgrna_specific_edits(run, verbose=False)
        assert counts["references"]["Reference"]["total_reads"] == 200
        rebuilt = c.crispresso_reconstitute_reads(run, str(path / "re.fq"), verbose=False)
        assert rebuilt["matches_original"]
