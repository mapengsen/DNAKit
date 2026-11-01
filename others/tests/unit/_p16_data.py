"""Deterministic toy data and tool lookup shared by the P16 command-line facade tests (no network)."""

from __future__ import annotations

import gzip
import random
import shutil
from pathlib import Path

import pytest


def need(program: str) -> str:
    """Path of ``program`` or skip the test."""
    found = shutil.which(program)
    if found is None:
        pytest.skip(f"{program} is not on PATH")
    return found


def rand_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path: Path, records: dict[str, str]) -> str:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in records.items()), encoding="utf-8")
    return str(path)


def write_fastq(path: Path, records: list[tuple[str, str]], quality: str = "I") -> str:
    path.write_text("".join(f"@{h}\n{s}\n+\n{quality * len(s)}\n" for h, s in records), encoding="utf-8")
    return str(path)


def genome(directory: Path) -> dict[str, str]:
    """Two random contigs (20 kb, 12 kb) with an A homopolymer; FASTA + .fai in ``directory``."""
    pysam = pytest.importorskip("pysam")
    rng = random.Random(11)
    contigs = {"chrA": rand_dna(rng, 20000), "chrB": rand_dna(rng, 12000)}
    contigs["chrA"] = contigs["chrA"][:1000] + "A" * 12 + contigs["chrA"][1012:]
    write_fasta(directory / "genome.fa", contigs)
    pysam.faidx(str(directory / "genome.fa"))
    return contigs


def read_bam(directory: Path, contigs: dict[str, str]) -> str:
    """Sorted, indexed BAM of 150 bp error-free reads every 25 bp; some reads carry a deletion in the homopolymer."""
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": k, "LN": len(v)} for k, v in contigs.items()]}
    path = directory / "reads.bam"
    n = 0
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for ref_id, (name, seq) in enumerate(contigs.items()):
            for pos in range(0, len(seq) - 150, 25):
                read, cigar = seq[pos : pos + 150], [(0, 150)]
                if name == "chrA" and 960 < pos < 1000:
                    read = seq[pos:1008] + seq[1009 : pos + 151]
                    cigar = [(0, 1008 - pos), (2, 1), (0, len(read) - (1008 - pos))]
                segment = pysam.AlignedSegment()
                segment.query_name, segment.query_sequence, segment.flag = f"r{n}", read, 0
                segment.reference_id, segment.reference_start, segment.mapping_quality = ref_id, pos, 60
                segment.cigar = cigar
                segment.query_qualities = pysam.qualitystring_to_array("I" * len(read))
                out.write(segment)
                n += 1
    pysam.index(str(path))
    return str(path)


def alignment_fasta(directory: Path) -> str:
    rng = random.Random(5)
    base = rand_dna(rng, 300)

    def mutate(seq: str, count: int) -> str:
        letters = list(seq)
        for _ in range(count):
            letters[rng.randrange(len(letters))] = rng.choice("ACGT")
        return "".join(letters)

    return write_fasta(directory / "aln.fa", {"t1": base, "t2": mutate(base, 5), "t3": mutate(base, 10), "t4": mutate(base, 20), "t5": mutate(base, 30)})


def long_reads(directory: Path, contigs: dict[str, str]) -> tuple[str, str]:
    """40 noisy long reads with ``start_time`` headers: plain and gzip FASTQ."""
    rng = random.Random(3)
    records = []
    for i in range(40):
        name = rng.choice(list(contigs))
        length = rng.randint(1500, 3000)
        start = rng.randrange(len(contigs[name]) - length)
        letters = list(contigs[name][start : start + length])
        for _ in range(length // 30):
            letters[rng.randrange(length)] = rng.choice("ACGT")
        records.append((f"lr{i} runid=x read={i} ch=1 start_time=2024-01-0{1 + i % 5}T10:00:00Z", "".join(letters)))
    plain = write_fastq(directory / "long.fastq", records, "5")
    zipped = directory / "long.fastq.gz"
    with gzip.open(zipped, "wt") as handle:
        handle.write(Path(plain).read_text())
    return plain, str(zipped)


def kaiju_data(directory: Path) -> dict[str, str]:
    """Three proteins, reads back-translated from them, and a small taxonomy (nodes.dmp / names.dmp)."""
    codons = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
        "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
    }
    rng = random.Random(5)
    proteins = {f"p{i}_{taxon}": "".join(rng.choice(list(codons)) for _ in range(300)) for i, taxon in enumerate((100, 200, 300))}
    write_fasta(directory / "db.faa", proteins)
    reads = {f"read_{name}_{j}": "".join(codons[a] for a in seq[j * 30 : j * 30 + 50]) for name, seq in proteins.items() for j in range(3)}
    reads["read_random"] = rand_dna(rng, 150)
    write_fasta(directory / "reads.fa", reads)
    (directory / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n10\t|\t1\t|\tsuperkingdom\t|\n20\t|\t10\t|\tphylum\t|\n30\t|\t20\t|\tclass\t|\n40\t|\t30\t|\torder\t|\n"
        "50\t|\t40\t|\tfamily\t|\n60\t|\t50\t|\tgenus\t|\n100\t|\t60\t|\tspecies\t|\n200\t|\t60\t|\tspecies\t|\n300\t|\t50\t|\tgenus\t|\n"
    )
    labels = ((1, "root"), (10, "Bacteria"), (20, "Phylum_A"), (30, "Class_A"), (40, "Order_A"), (50, "Family_A"), (60, "Genus_A"),
              (100, "Species_one"), (200, "Species_two"), (300, "Genus_B"))
    (directory / "names.dmp").write_text("".join(f"{i}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for i, n in labels))
    return {name: str(directory / name) for name in ("db.faa", "reads.fa", "nodes.dmp", "names.dmp")}


def noisy_reads(contigs: dict[str, str], coverage: int, path: Path, seed: int = 8) -> str:
    """Noisy long reads (3 % errors, both strands) at the given coverage."""
    rng = random.Random(seed)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    records, n = [], 0
    for seq in contigs.values():
        for _ in range(int(len(seq) * coverage / 3000)):
            length = min(rng.randint(2000, 4500), len(seq) - 1)
            start = rng.randrange(len(seq) - length)
            out: list[str] = []
            for ch in seq[start : start + length]:
                r = rng.random()
                if r < 0.01:
                    continue
                out.append(rng.choice("ACGT") if r < 0.02 else ch)
                if 0.02 <= r < 0.03:
                    out.append(rng.choice("ACGT"))
            if rng.random() < 0.5:
                out = [comp[c] for c in reversed(out)]
            records.append((f"lr{n}", "".join(out)))
            n += 1
    return write_fastq(path, records, "5")


def sv_bam(directory: Path, contigs: dict[str, str]) -> str:
    """Long reads of a sample with a 400 bp deletion (chrA:8000) and a 300 bp insertion (chrA:14000), mapped with minimap2."""
    import subprocess

    minimap2, samtools = need("minimap2"), need("samtools")
    rng = random.Random(9)
    a = contigs["chrA"]
    sample = {"chrA": a[:8000] + a[8400:14000] + rand_dna(rng, 300) + a[14000:], "chrB": contigs["chrB"]}
    reads = noisy_reads(sample, 30, directory / "sample_reads.fq")
    sam = subprocess.run([minimap2, "-ax", "map-ont", "-t", "4", str(directory / "genome.fa"), reads], capture_output=True, check=True).stdout
    bam = directory / "sample.bam"
    subprocess.run([samtools, "sort", "-o", str(bam), "-"], input=sam, check=True, capture_output=True)
    subprocess.run([samtools, "index", str(bam)], check=True)
    return str(bam)


def bisulfite_data(directory: Path) -> dict[str, str]:
    """A 1.5 kb reference rich in CpG, a directional bisulfite BAM (YD tag) and a bigWig of mappability 1."""
    pysam = pytest.importorskip("pysam")
    rng = random.Random(4)
    ref = "".join(rng.choice(["A", "C", "G", "T", "CG"]) for _ in range(1500))
    write_fasta(directory / "ref.fa", {"chrM": ref})
    pysam.faidx(str(directory / "ref.fa"))
    methylated = {i: rng.random() < 0.7 for i in range(len(ref) - 1) if ref[i : i + 2] == "CG"}
    reads = []
    for pos in range(0, len(ref) - 100, 3):
        for strand in (0, 16):
            letters = list(ref[pos : pos + 100])
            for j, ch in enumerate(letters):
                if strand == 0 and ch == "C" and not methylated.get(pos + j, False):
                    letters[j] = "T"
                if strand == 16 and ch == "G" and not methylated.get(pos + j - 1, False):
                    letters[j] = "A"
            reads.append((pos, strand, "".join(letters)))
    reads.sort()
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrM", "LN": len(ref)}]}
    with pysam.AlignmentFile(str(directory / "bs.bam"), "wb", header=header) as out:
        for n, (pos, strand, seq) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag, a.reference_id, a.reference_start = f"r{n}", seq, strand, 0, pos
            a.mapping_quality, a.cigar = 60, [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.set_tag("YD", "f" if strand == 0 else "r")
            out.write(a)
    pysam.index(str(directory / "bs.bam"))
    bigwig = pytest.importorskip("pyBigWig").open(str(directory / "map.bw"), "w")
    bigwig.addHeader([("chrM", len(ref))])
    bigwig.addEntries(["chrM"], [0], ends=[len(ref)], values=[1.0])
    bigwig.close()
    return {name: str(directory / name) for name in ("ref.fa", "bs.bam", "map.bw")}


def mageck_data(directory: Path) -> dict[str, str]:
    """40 genes x 4 sgRNAs, four FASTQ samples (G1 and G2 depleted in the treatment), a design matrix and gene sets."""
    rng = random.Random(2)
    genes = [f"G{i}" for i in range(40)]
    library = [(f"{g}_sg{j}", rand_dna(rng, 20), g) for g in genes for j in range(4)]
    (directory / "lib.txt").write_text("".join(f"{a}\t{b}\t{c}\n" for a, b, c in library))

    def sample(name: str, factor: float) -> None:
        records, n = [], 0
        for _sid, seq, gene in library:
            count = rng.randint(80, 150)
            count = int(count * factor) if gene in {"G1", "G2"} else count
            for _ in range(count):
                records.append((f"r{n}", "TT" + seq + "AA"))
                n += 1
        write_fastq(directory / name, records)

    for name, factor in (("c1.fq", 1.0), ("c2.fq", 1.0), ("t1.fq", 0.1), ("t2.fq", 0.1)):
        sample(name, factor)
    (directory / "design.txt").write_text("Samples\tbaseline\ttreat\nc1\t1\t0\nc2\t1\t0\nt1\t1\t1\nt2\t1\t1\n")
    (directory / "sets.gmt").write_text("setA\tna\t" + "\t".join(genes[:10]) + "\nsetB\tna\t" + "\t".join(genes[10:25]) + "\n")
    lines = ["GENE_CLONE\tGENE\tc1\tc2\td1\td2"]
    for sid, _seq, gene in library:
        base = rng.randint(100, 400)
        factor = 0.3 if gene == "G5" else 1.0
        lines.append(f"{gene}_{sid}\t{gene}\t{base + rng.randint(0, 20)}\t{base + rng.randint(0, 20)}\t{int(base * factor) + rng.randint(0, 20)}\t{int(base * factor) + rng.randint(0, 20)}")
    (directory / "drugz_counts.txt").write_text("\n".join(lines) + "\n")
    return {p.name: str(p) for p in directory.iterdir()}
