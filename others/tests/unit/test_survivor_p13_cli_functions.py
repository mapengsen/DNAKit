import random
import shutil

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import survivor_p13_cli_functions as s

needs = pytest.mark.skipif(shutil.which("SURVIVOR") is None, reason="SURVIVOR not installed")


@pytest.fixture()
def genome(tmp_path):
    rng = random.Random(4)
    path = tmp_path / "ref.fa"
    path.write_text("".join(f">c{i}\n" + "".join(rng.choice("ACGT") for _ in range(100000)) + "\n" for i in range(2)))
    return path


def test_validation_and_missing_program(tmp_path, genome):
    with pytest.raises(BackendUnavailableError):
        s.survivor_simulate_sv_genome(str(genome), str(tmp_path / "x"), survivor_path=str(tmp_path / "none"), verbose=False)
    with pytest.raises(ConfigurationError):
        s.survivor_simulate_sv_genome(str(genome), str(tmp_path / "x"), duplication=[10, 5, 1], verbose=False)
    with pytest.raises(ConfigurationError):
        s.survivor_convert_to_sv_vcf(str(genome), str(tmp_path / "o.vcf"), input_format="bed", verbose=False)
    with pytest.raises(ConfigurationError):
        s.survivor_convert_to_sv_vcf(str(genome), str(tmp_path / "o.vcf"), input_format="smap", sv_type="DEL", verbose=False)


@needs
def test_simulate_evaluate_filter_stats_convert(tmp_path, genome):
    params = {"duplication": [100, 500, 2], "indel": [50, 200, 2], "translocation": [500, 900, 1], "inversion": [300, 400, 2], "inversion_deletion": [300, 400, 1], "inversion_duplication": [300, 400, 1]}
    sim = s.survivor_simulate_sv_genome(str(genome), str(tmp_path / "sim"), snp_frequency=0.001, verbose=False, **params)
    assert sim["status"] == "completed" and sim["event_count"] >= 8
    ev = s.survivor_evaluate_sv_calls(str(tmp_path / "sim.vcf"), str(tmp_path / "sim.bed"), str(tmp_path / "ev"), max_distance=100, verbose=False)
    assert ev["evaluation"]["sensitivity"] == 1.0
    flt = s.survivor_filter_sv_vcf(str(tmp_path / "sim.vcf"), str(tmp_path / "f.vcf"), min_size=100, verbose=False)
    assert 0 < flt["records_out"] <= flt["records_in"]
    st = s.survivor_sv_vcf_stats(str(tmp_path / "sim.vcf"), str(tmp_path / "st.txt"), verbose=False)
    assert st["totals"]["Tot"] == sim["event_count"] and st["length_table"]
    bed = s.survivor_sv_vcf_to_bed(str(tmp_path / "sim.vcf"), str(tmp_path / "v.bed"), min_size=50, verbose=False)
    assert bed["record_count"] >= 1
    (tmp_path / "r.bed").write_text("chr1\t1000\t2000\n")
    conv = s.survivor_convert_to_sv_vcf(str(tmp_path / "r.bed"), str(tmp_path / "c.vcf"), sv_type="DEL", verbose=False)
    assert conv["record_count"] == 1
    (tmp_path / "cov.txt").write_text("chr1\t100\t3\nchr1\t101\t3\nchr1\t500\t4\n")
    bc = s.survivor_low_mapping_quality_regions(str(tmp_path / "cov.txt"), str(tmp_path / "bc.bed"), verbose=False)
    assert bc["region_count"] == 1
    (tmp_path / "files.txt").write_text("")


@needs
def test_error_profile(tmp_path, genome):
    rng = random.Random(1)
    seqs = {}
    name = None
    for line in genome.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    rows = []
    ref = seqs["c0"]
    for i in range(30):
        start = rng.randrange(0, 90000)
        read, md, run, mismatches = [], [], 0, 0
        for base in ref[start : start + 3000]:
            if rng.random() < 0.03:
                read.append(rng.choice([x for x in "ACGT" if x != base]))
                md += [str(run), base]
                run = 0
                mismatches += 1
            else:
                read.append(base)
                run += 1
        md.append(str(run))
        read = "".join(read)
        rows.append(f"r{i}\t0\tc0\t{start + 1}\t60\t3000M\t*\t0\t0\t{read}\t{'I' * 3000}\tNM:i:{mismatches}\tMD:Z:{''.join(md)}")
    sam = tmp_path / "long.sam"
    sam.write_text("\n".join(rows) + "\n")
    prof = s.survivor_scan_read_error_profile(str(sam), str(tmp_path / "err.txt"), verbose=False)
    assert prof["valid_reads"] == 30
    assert (tmp_path / "err.txt").read_text().startswith("Pos")
