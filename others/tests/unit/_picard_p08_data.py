"""Small deterministic Picard test data for the P08 facade tests (reference, reads, maps, VCFs)."""

from __future__ import annotations

import os
import random
import shutil
from pathlib import Path

READ_LENGTH = 36
_BASES = "ACGT"


def picard_available() -> bool:
    """True when Picard and (for a jar) Java can be started."""
    if shutil.which("picard") is not None:
        return True
    jar = os.environ.get("PICARD_JAR")
    return bool(jar) and Path(jar).is_file() and shutil.which("java") is not None


def _rnd(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(_BASES) for _ in range(n))


def _revcomp(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def build_reference(rng: random.Random) -> dict[str, str]:
    return {"chr1": _rnd(rng, 1200), "chr2": _rnd(rng, 800)}


def build_pairs(rng: random.Random, reference: dict[str, str]) -> list[list[object]]:
    """22 read pairs of 36 bp (two exact duplicate pairs), coordinate sorted, read group rg1."""
    records: list[list[object]] = []
    index = 0
    cigar = f"{READ_LENGTH}M"
    for chrom, count in (("chr1", 14), ("chr2", 6)):
        for _ in range(count):
            position = rng.randint(10, len(reference[chrom]) - 300)
            insert = rng.randint(150, 230)
            start1, start2 = position, position + insert - READ_LENGTH
            seq1 = reference[chrom][start1 : start1 + READ_LENGTH]
            seq2 = _revcomp(reference[chrom][start2 : start2 + READ_LENGTH])
            qual1 = "".join(rng.choice("5?@ABCDEFGHI") for _ in range(READ_LENGTH))
            qual2 = "".join(rng.choice("5?@ABCDEFGHI") for _ in range(READ_LENGTH))
            name = f"M1:7:FC1:1:1101:{1000 + index * 37}:{2000 + index * 11}"
            records.append(
                [
                    name,
                    99,
                    chrom,
                    start1 + 1,
                    60,
                    cigar,
                    "=",
                    start2 + 1,
                    insert,
                    seq1,
                    qual1,
                    "RG:Z:rg1",
                    f"MC:Z:{cigar}",
                    "NM:i:0",
                ]
            )
            records.append(
                [
                    name,
                    147,
                    chrom,
                    start2 + 1,
                    60,
                    cigar,
                    "=",
                    start1 + 1,
                    -insert,
                    seq2,
                    qual2,
                    "RG:Z:rg1",
                    f"MC:Z:{cigar}",
                    "NM:i:0",
                ]
            )
            index += 1
    for source in (0, 3):
        name = f"M1:7:FC1:1:1101:{1001 + source}:{2001 + source}"
        for row in records[2 * source : 2 * source + 2]:
            clone = list(row)
            clone[0] = name
            records.append(clone)
    order = list(reference)
    records.sort(key=lambda row: (order.index(str(row[2])), int(row[3])))
    return records


def write_sam(path: Path, reference: dict[str, str], records: list[list[object]]) -> None:
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{name}\tLN:{len(sequence)}" for name, sequence in reference.items()]
    lines.append("@RG\tID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:FC1.1")
    lines += ["\t".join(str(value) for value in row) for row in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_reference(directory: Path, reference: dict[str, str]) -> None:
    fasta, index, dictionary = [], [], ["@HD\tVN:1.6\tSO:unsorted"]
    offset = 0
    for name, sequence in reference.items():
        fasta.append(f">{name}")
        fasta += [sequence[i : i + 60] for i in range(0, len(sequence), 60)]
        offset += len(name) + 2
        index.append(f"{name}\t{len(sequence)}\t{offset}\t60\t61")
        offset += len(sequence) + (len(sequence) + 59) // 60
        dictionary.append(f"@SQ\tSN:{name}\tLN:{len(sequence)}")
    (directory / "ref.fa").write_text("\n".join(fasta) + "\n", encoding="utf-8")
    (directory / "ref.fa.fai").write_text("\n".join(index) + "\n", encoding="utf-8")
    (directory / "ref.dict").write_text("\n".join(dictionary) + "\n", encoding="utf-8")


def fingerprint_sites(
    reference: dict[str, str], records: list[list[object]], count: int = 8
) -> list[tuple[str, int, str, str]]:
    """Chr1 positions covered by two or more reads: (contig, position, reference, alternative)."""
    depth: dict[int, int] = {}
    for row in records:
        if row[2] != "chr1":
            continue
        for offset in range(READ_LENGTH):
            position = int(row[3]) + offset
            depth[position] = depth.get(position, 0) + 1
    sites = []
    for position in sorted(depth):
        if depth[position] >= 2 and len(sites) < count and position % 7 == 0:
            ref = reference["chr1"][position - 1]
            sites.append(("chr1", position, ref, _BASES[(_BASES.index(ref) + 1) % 4]))
    return sites


def write_haplotype_map(
    path: Path, reference: dict[str, str], sites: list[tuple[str, int, str, str]]
) -> None:
    lines = ["@HD\tVN:1.6"] + [
        f"@SQ\tSN:{name}\tLN:{len(sequence)}" for name, sequence in reference.items()
    ]
    lines.append("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS")
    for contig, position, ref, alt in sites:
        lines.append(f"{contig}\t{position}\tsnp{position}\t{ref}\t{alt}\t0.3\tsnp{position}\t")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_vcf(
    path: Path,
    reference: dict[str, str],
    sites: list[tuple[str, int, str, str]],
    sample: str,
    genotype: str = "0/0",
    *,
    ids: bool = False,
) -> None:
    lines = ["##fileformat=VCFv4.2"]
    lines += [
        f"##contig=<ID={name},length={len(sequence)}>" for name, sequence in reference.items()
    ]
    lines += [
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">',
        '##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    for index, (contig, position, ref, alt) in enumerate(sites):
        snp_id = f"rs{position}" if ids else "."
        lines.append(
            f"{contig}\t{position}\t{snp_id}\t{ref}\t{alt}\t50\t.\t"
            f"DP=10;QD={5 + index};FS=1.0\tGT:DP:GQ\t{genotype}:{8 + index}:{30 + index}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_dataset(directory: Path, seed: int = 7) -> dict[str, Path]:
    """Write the shared example files into ``directory`` and return their paths by role."""
    rng = random.Random(seed)
    reference = build_reference(rng)
    records = build_pairs(rng, reference)
    write_reference(directory, reference)
    write_sam(directory / "example.sam", reference, records)
    sites = fingerprint_sites(reference, records)
    write_haplotype_map(directory / "hap.map", reference, sites)
    write_vcf(directory / "genotypes.vcf", reference, sites, "sample1")
    write_vcf(directory / "sample.vcf", reference, sites, "sample1", "0/1")
    write_vcf(directory / "dbsnp.vcf", reference, sites, "dbsnp", ids=True)
    for name in ("genotypes", "sample", "dbsnp"):
        index_vcf(directory / f"{name}.vcf")
    return {
        "sam": directory / "example.sam",
        "fasta": directory / "ref.fa",
        "dict": directory / "ref.dict",
        "haplotype_map": directory / "hap.map",
        "genotypes": directory / "genotypes.vcf.gz",
        "vcf": directory / "sample.vcf.gz",
        "dbsnp": directory / "dbsnp.vcf.gz",
    }


def index_vcf(path: Path) -> Path:
    """Compress a VCF with bgzip and write its tabix index; the plain file is kept."""
    pysam = __import__("pysam")
    pysam.tabix_index(str(path), preset="vcf", force=True, keep_original=True)
    return Path(f"{path}.gz")


def sam_to_bam(sam: Path, bam: Path) -> Path:
    """Write a BGZF BAM (with the end-of-file block) from a SAM with pysam."""
    pysam = __import__("pysam")
    with (
        pysam.AlignmentFile(str(sam), "r") as source,
        pysam.AlignmentFile(str(bam), "wb", template=source) as target,
    ):
        for read in source:
            target.write(read)
    return bam


CONTROL_CODES = (
    ("DNP(High)", "Staining"),
    ("DNP(Bgnd)", "Staining"),
    ("Biotin(High)", "Staining"),
    ("Biotin(Bgnd)", "Staining"),
    ("Extension(A)", "Extension"),
    ("Extension(T)", "Extension"),
    ("Extension(C)", "Extension"),
    ("Extension(G)", "Extension"),
    ("TargetRemoval", "TargetRemoval"),
    ("Hyb(High)", "Hybridization"),
    ("Hyb(Medium)", "Hybridization"),
    ("Hyb(Low)", "Hybridization"),
    ("String(PM)", "Stringency"),
    ("String(MM)", "Stringency"),
    ("NSB(Bgnd)Red", "Non-SpecificBinding"),
    ("NSB(Bgnd)Purple", "Non-SpecificBinding"),
    ("NSB(Bgnd)Blue", "Non-SpecificBinding"),
    ("NSB(Bgnd)Green", "Non-SpecificBinding"),
    ("NP(A)", "Non-Polymorphic"),
    ("NP(T)", "Non-Polymorphic"),
    ("NP(C)", "Non-Polymorphic"),
    ("NP(G)", "Non-Polymorphic"),
    ("Restore", "Restoration"),
)


def write_arrays_vcf(
    path: Path, reference: dict[str, str], sites: list[tuple[str, int, str, str]]
) -> Path:
    """A genotyping-array VCF with the Infinium header lines Picard requires."""
    lines = ["##fileformat=VCFv4.2"]
    lines += [
        f"##contig=<ID={name},length={len(sequence)}>" for name, sequence in reference.items()
    ]
    lines += [f"##{control}={control}|{category}|1|2" for control, category in CONTROL_CODES]
    meta = {
        "sampleAlias": "sample1",
        "arrayType": "TestChip",
        "autocallGender": "F",
        "autocallVersion": "1.0",
        "autocallDate": "09/21/2016 20:40",
        "imagingDate": "8/15/2015 7:28:52 AM",
        "extendedIlluminaManifestVersion": "1.4",
        "clusterFile": "test.egt",
        "p95Green": "1000",
        "p95Red": "900",
        "scannerName": "iScan1",
    }
    lines += [f"##{key}={value}" for key, value in meta.items()]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1")
    for index, (contig, position, ref, alt) in enumerate(sites):
        lines.append(
            f"{contig}\t{position}\trs{position}\t{ref}\t{alt}\t.\t.\t.\tGT\t"
            f"{'0/0' if index % 2 else '0/1'}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_illumina_run(directory: Path) -> Path:
    """A tiny Illumina run folder (4 clusters, 2 cycles, tile 1101, lane 1); returns BaseCalls."""
    import struct

    clusters = 4
    basecalls = directory / "Data" / "Intensities" / "BaseCalls"
    lane = basecalls / "L001"
    for cycle in (1, 2):
        (lane / f"C{cycle}.1").mkdir(parents=True)
    (directory / "InterOp").mkdir()
    (directory / "InterOp" / "TileMetricsOut.bin").write_bytes(
        bytes([2, 10]) + struct.pack("<HHHf", 1, 1101, 100, 1000.0)
    )
    (lane / "s_1_1101.filter").write_bytes(
        struct.pack("<III", 0, 3, clusters) + bytes([1, 1, 0, 1])
    )
    for cycle, bases in ((1, (0, 1, 2, 3)), (2, (3, 2, 1, 0))):
        (lane / f"C{cycle}.1" / "s_1_1101.bcl").write_bytes(
            struct.pack("<I", clusters) + bytes((30 << 2) | base for base in bases)
        )
    locs = directory / "Data" / "Intensities" / "L001"
    locs.mkdir()
    (locs / "s_1_1101.locs").write_bytes(
        struct.pack("<IfI", 1, 1.0, clusters)
        + b"".join(struct.pack("<ff", 100.0 + i * 10, 200.0 + i * 10) for i in range(clusters))
    )
    return basecalls


__all__ = [
    "CONTROL_CODES",
    "build_pairs",
    "build_reference",
    "fingerprint_sites",
    "index_vcf",
    "picard_available",
    "sam_to_bam",
    "write_arrays_vcf",
    "write_dataset",
    "write_haplotype_map",
    "write_illumina_run",
    "write_sam",
    "write_vcf",
]
