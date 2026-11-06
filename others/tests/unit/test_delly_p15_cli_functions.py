from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import delly_p15_cli_functions as dl

pysam = pytest.importorskip("pysam")
needs_delly = pytest.mark.skipif(shutil.which("delly") is None, reason="delly is not installed")


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def build(folder: Path) -> dict[str, Path]:
    """Reference, short-read BAM, long-read BAM and assembly BAM that all carry one 800 bp deletion at 30000."""
    random.seed(11)
    length, start, size = 60000, 30000, 800
    ref = "".join(random.choice("ACGT") for _ in range(length))
    donor = ref[:start] + ref[start + size :]
    (folder / "ref.fa").write_text(
        ">chr1\n" + "\n".join(ref[i : i + 70] for i in range(0, length, 70)) + "\n"
    )
    pysam.faidx(str(folder / "ref.fa"))
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": length}],
        "RG": [{"ID": "s1", "SM": "sample1"}],
    }

    def seg(name, seq, pos, cigar, flag, tags=()):
        a = pysam.AlignedSegment()
        a.query_name, a.query_sequence, a.flag = name, seq, flag
        a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = 0, pos, 60, cigar
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.set_tag("RG", "s1")
        for tag, value in tags:
            a.set_tag(tag, value)
        return a

    def to_ref(pos):
        return pos if pos < start else pos + size

    def write(name, records):
        records.sort(key=lambda r: (r.reference_start, r.query_name))
        with pysam.AlignmentFile(str(folder / name), "wb", header=header) as out:
            for r in records:
                out.write(r)
        pysam.index(str(folder / name))

    short = []
    for i in range(int(len(donor) * 25 / 200)):
        insert = max(200, int(random.gauss(350, 30)))
        s = random.randrange(0, len(donor) - insert)
        frag = donor[s : s + insert]
        pair = []
        for seq, dpos, flag in (
            (frag[:100], s, 0x1 | 0x2 | 0x20 | 0x40),
            (revcomp(frag)[:100], s + insert - 100, 0x1 | 0x2 | 0x10 | 0x80),
        ):
            if dpos < start < dpos + 100:
                k = start - dpos
                pair.append(
                    seg(
                        f"p{i}",
                        seq,
                        to_ref(dpos),
                        f"{k}M{100 - k}S",
                        flag,
                        [("SA", f"chr1,{start + size + 1},+,{k}S{100 - k}M,60,0;")],
                    )
                )
            else:
                pair.append(seg(f"p{i}", seq, to_ref(dpos), "100M", flag))
        tlen = pair[1].reference_start + 100 - pair[0].reference_start
        for x, y, t in ((pair[0], pair[1], tlen), (pair[1], pair[0], -tlen)):
            x.next_reference_id, x.next_reference_start, x.template_length = 0, y.reference_start, t
            if tlen > 700:
                x.flag &= ~0x2
        short += pair
        for x in pair:
            if x.has_tag("SA"):
                k = x.cigartuples[0][1]
                sup = seg(
                    x.query_name,
                    x.query_sequence,
                    start + size,
                    f"{k}S{100 - k}M",
                    x.flag | 0x800,
                    [("SA", f"chr1,{x.reference_start + 1},+,{k}M{100 - k}S,60,0;")],
                )
                sup.next_reference_id, sup.next_reference_start, sup.template_length = (
                    0,
                    x.next_reference_start,
                    x.template_length,
                )
                short.append(sup)
    write("short.bam", short)
    long_reads = []
    for i in range(20):
        s = start - random.randrange(1500, 3500)
        seq = donor[s : s + 6000]
        k = start - s
        long_reads.append(
            seg(
                f"lr{i}",
                seq,
                s,
                f"{k}M{len(seq) - k}S",
                0,
                [("SA", f"chr1,{start + size + 1},+,{k}S{len(seq) - k}M,60,0;")],
            )
        )
        long_reads.append(
            seg(
                f"lr{i}",
                seq,
                start + size,
                f"{k}S{len(seq) - k}M",
                0x800,
                [("SA", f"chr1,{s + 1},+,{k}M{len(seq) - k}S,60,0;")],
            )
        )
    write("long.bam", long_reads)
    contig = donor[start - 10000 : start + 10000]
    write(
        "asm.bam",
        [
            seg(
                "contig1",
                contig,
                start - 10000,
                "10000M10000S",
                0,
                [("SA", f"chr1,{start + size + 1},+,10000S10000M,60,0;")],
            ),
            seg(
                "contig1",
                contig,
                start + size,
                "10000S10000M",
                0x800,
                [("SA", f"chr1,{start - 9999},+,10000M10000S,60,0;")],
            ),
        ],
    )
    return {name: folder / name for name in ("ref.fa", "short.bam", "long.bam", "asm.bam")}


def test_missing_program_and_validation(tmp_path: Path) -> None:
    f = tmp_path / "x.bcf"
    f.write_text("x")
    with pytest.raises(BackendUnavailableError):
        dl.delly_filter_call_set(
            str(f),
            str(tmp_path / "o.bcf"),
            filter_mode="germline",
            delly_path=str(tmp_path / "missing"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        dl.delly_filter_call_set(
            str(f), str(tmp_path / "o2.bcf"), verbose=False
        )  # somatic needs samples_file
    with pytest.raises(ConfigurationError):
        dl.delly_genotype_sv_sites(
            [str(f)], str(f), str(f), str(tmp_path / "o3.bcf"), mode="x", verbose=False
        )
    with pytest.raises(ConfigurationError):
        dl.delly_plot_read_depth(
            str(f), str(tmp_path / "plots"), str(tmp_path / "nonexistent.R"), verbose=False
        )


@needs_delly
def test_genotype_assembly_cnv_merge_filter(tmp_path: Path) -> None:
    data = build(tmp_path)
    import subprocess

    delly = shutil.which("delly")
    sr, lr = tmp_path / "sr.bcf", tmp_path / "lr.bcf"
    subprocess.run(
        [delly, "sr", "-g", str(data["ref.fa"]), "-o", str(sr), str(data["short.bam"])],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [delly, "lr", "-y", "pb", "-g", str(data["ref.fa"]), "-o", str(lr), str(data["long.bam"])],
        check=True,
        capture_output=True,
    )
    genotyped = dl.delly_genotype_sv_sites(
        [str(data["short.bam"])],
        str(data["ref.fa"]),
        str(sr),
        str(tmp_path / "g.bcf"),
        verbose=False,
    )
    assert genotyped["variant_count"] == 1 and genotyped["variants"][0]["type"] == "DEL"
    long_genotyped = dl.delly_genotype_sv_sites(
        [str(data["long.bam"])],
        str(data["ref.fa"]),
        str(lr),
        str(tmp_path / "gl.bcf"),
        mode="lr",
        technology="pb",
        verbose=False,
    )
    assert long_genotyped["variant_count"] == 1
    assembly = dl.delly_assembly_sv(
        [str(data["asm.bam"])], str(data["ref.fa"]), str(tmp_path / "a.bcf"), verbose=False
    )
    assert (
        assembly["variants"][0]["type"] == "DEL"
        and abs(assembly["variants"][0]["start"] - 30000) <= 5
    )
    cnv = dl.delly_copy_number_variants(
        str(data["short.bam"]),
        str(data["ref.fa"]),
        str(tmp_path / "c.bcf"),
        coverage_file=str(tmp_path / "c.cov.gz"),
        segmentation_file=str(tmp_path / "c.bed"),
        verbose=False,
    )
    assert cnv["status"] == "completed" and (tmp_path / "c.cov.gz").is_file()
    merged = dl.delly_merge_call_sets(
        [str(sr), str(lr)], str(tmp_path / "m.bcf"), min_quality=0, verbose=False
    )
    assert merged["variant_count"] == 1
    filtered = dl.delly_filter_call_set(
        str(sr), str(tmp_path / "f.bcf"), filter_mode="germline", min_quality=0, verbose=False
    )
    assert filtered["status"] == "completed"
