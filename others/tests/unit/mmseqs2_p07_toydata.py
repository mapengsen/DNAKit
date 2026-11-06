"""Deterministic toy inputs for the MMseqs2 P07 facade tests and documentation examples."""

from __future__ import annotations

import random
from pathlib import Path

AMINO = "ACDEFGHIKLMNPQRSTVWY"
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
    "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
    "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _mutate(rng: random.Random, protein: str, rate: float) -> str:
    return "".join(rng.choice(AMINO) if rng.random() < rate else residue for residue in protein)


def protein_families(
    *, families: int = 6, members: int = 4, length: int = 120, seed: int = 7, rate: float = 0.1
) -> dict[str, str]:
    """``families`` unrelated random proteins, each with ``members`` mutated copies (``rate`` per residue)."""
    rng = random.Random(seed)
    records: dict[str, str] = {}
    for family in range(families):
        root = "M" + "".join(rng.choice(AMINO) for _ in range(length - 1))
        for member in range(members):
            records[f"fam{family}_m{member}"] = root if member == 0 else _mutate(rng, root, rate)
    return records


def write_fasta(path: str | Path, records: dict[str, str]) -> Path:
    target = Path(path)
    target.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in records.items()), encoding="utf-8"
    )
    return target


def write_protein_fasta(path: str | Path, **kwargs: int) -> Path:
    return write_fasta(path, protein_families(**kwargs))


def contigs(*, seed: int = 11) -> tuple[dict[str, str], list[tuple[str, int, int, str]]]:
    """Four contigs, each with two genes (one on each strand) that encode proteins of ``protein_families()``.

    Returns the contig sequences and the gene features (contig, 1-based start, end, strand).
    """
    rng = random.Random(seed)
    reference = protein_families()
    names = ["fam0_m0", "fam3_m0", "fam1_m0", "fam4_m0", "fam2_m0", "fam5_m0", "fam0_m1", "fam3_m1"]
    proteins = [reference[name] for name in names]
    sequences: dict[str, str] = {}
    features: list[tuple[str, int, int, str]] = []
    for index in range(4):
        parts = ""
        for gene in range(2):
            spacer = "".join(rng.choice("ACGT") for _ in range(60))
            protein = proteins[index * 2 + gene]
            cds = "".join(CODON[residue] for residue in protein) + "TAA"
            strand = "+" if gene == 0 else "-"
            start = len(parts) + len(spacer) + 1
            parts += spacer + (cds if strand == "+" else cds.translate(COMPLEMENT)[::-1])
            features.append((f"contig{index}", start, start + len(cds) - 1, strand))
        parts += "".join(rng.choice("ACGT") for _ in range(60))
        sequences[f"contig{index}"] = parts
    return sequences, features


def write_contigs(path: str | Path) -> Path:
    return write_fasta(path, contigs()[0])


def write_cds_fasta(path: str | Path) -> Path:
    """Back-translated coding sequences of ``protein_families()`` with the same record names."""
    return write_fasta(
        path,
        {
            name: "".join(CODON[residue] for residue in protein)
            for name, protein in protein_families().items()
        },
    )


def write_gff3(path: str | Path) -> Path:
    lines = ["##gff-version 3"]
    for number, (contig, start, end, strand) in enumerate(contigs()[1]):
        lines.append(f"{contig}\ttoy\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds{number}")
    target = Path(path)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def write_taxdump(directory: str | Path) -> Path:
    """A five-node taxonomy: root, Bacteria, Eukaryota, Escherichia (genus), Homo (genus)."""
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    nodes = [
        (1, 1, "no rank"),
        (2, 1, "superkingdom"),
        (2759, 1, "superkingdom"),
        (561, 2, "genus"),
        (9605, 2759, "genus"),
    ]
    names = {1: "root", 2: "Bacteria", 2759: "Eukaryota", 561: "Escherichia", 9605: "Homo"}
    (target / "nodes.dmp").write_text(
        "".join(f"{node}\t|\t{parent}\t|\t{rank}\t|\t\t|\n" for node, parent, rank in nodes),
        encoding="utf-8",
    )
    (target / "names.dmp").write_text(
        "".join(f"{node}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for node, name in names.items()),
        encoding="utf-8",
    )
    (target / "merged.dmp").write_text("", encoding="utf-8")
    (target / "delnodes.dmp").write_text("", encoding="utf-8")
    return target


def write_tax_mapping(path: str | Path, *, prefix: str = "fam") -> Path:
    """Families 0-2 belong to Escherichia (561), families 3-5 to Homo (9605)."""
    lines = []
    for name in protein_families():
        family = int(name.split("_")[0].removeprefix(prefix))
        lines.append(f"{name}\t{561 if family < 3 else 9605}")
    target = Path(path)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


FAKE_MMSEQS = '''#!/usr/bin/env python3
"""A stand-in for ``mmseqs`` that records its command line and creates the outputs it is asked for."""
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
if args[:1] == ["version"]:
    print("fake-mmseqs-1")
    sys.exit(0)
with open(os.environ["FAKE_MMSEQS_LOG"], "a", encoding="utf-8") as handle:
    handle.write(json.dumps(args) + "\\n")
if os.environ.get("FAKE_MMSEQS_FAIL"):
    print("simulated failure")
    sys.exit(3)
for item in args[1:]:
    if not item.startswith("/"):
        continue
    path = Path(item)
    if path.exists() or Path(item + ".dbtype").exists() or not path.parent.is_dir():
        continue
    if path.name.startswith(".mmseqs_tmp_"):
        continue
    path.write_text("a\\tb\\n", encoding="utf-8")
    Path(item + ".dbtype").write_bytes(b"\\x00\\x00\\x00\\x00")
    Path(item + ".index").write_text("0\\t0\\t4\\n", encoding="utf-8")
for item in filter(None, os.environ.get("FAKE_MMSEQS_TOUCH", "").split(",")):
    Path(item).write_text("x\\ty\\n", encoding="utf-8")
'''


def make_fake_mmseqs(directory: str | Path) -> Path:
    """Write the stand-in executable into ``directory`` and return its path."""
    target = Path(directory) / "mmseqs"
    target.write_text(FAKE_MMSEQS, encoding="utf-8")
    target.chmod(0o755)
    return target


def make_database(prefix: str | Path) -> Path:
    """Create the three files that identify an MMseqs2 database (enough for the stand-in)."""
    base = Path(prefix)
    base.write_text(">a\nACGT\n", encoding="utf-8")
    Path(f"{base}.dbtype").write_bytes(b"\x00\x00\x00\x00")
    Path(f"{base}.index").write_text("0\t0\t5\n", encoding="utf-8")
    return base


def run_direct(*args: str | Path, threads: bool = True) -> str:
    """Run the real ``mmseqs`` directly (the reference side of the equivalence tests)."""
    import subprocess

    completed = subprocess.run(
        ["mmseqs", *map(str, args), *(["--threads", "1"] if threads else []), "-v", "1"],
        capture_output=True,
        text=True,
        check=True,
    )
    return completed.stdout + completed.stderr


def database_files(prefix: str | Path) -> dict[str, object]:
    """The content of a database: its type bytes and every entry (key to bytes).

    Entries are read through ``<prefix>.index``; the data is the file ``<prefix>`` or, when MMseqs2 wrote
    per-thread parts, the parts ``<prefix>.0``, ``<prefix>.1``, ... joined in numeric order. Which part
    holds an entry differs between runs, the entries do not.
    """
    base = Path(prefix)
    if base.is_file():
        data = base.read_bytes()
    else:
        parts = sorted(
            (item for item in base.parent.glob(base.name + ".*") if item.suffix[1:].isdigit()),
            key=lambda item: int(item.suffix[1:]),
        )
        data = b"".join(item.read_bytes() for item in parts)
    entries: dict[str, bytes] = {}
    index = Path(f"{base}.index")
    for line in index.read_text(encoding="utf-8").splitlines():
        key, offset, length = line.split("\t")[:3]
        entries[key] = data[int(offset) : int(offset) + int(length)]
    dbtype = Path(f"{base}.dbtype")
    return {"dbtype": dbtype.read_bytes() if dbtype.exists() else b"", "entries": entries}


def build_world(directory: str | Path) -> dict[str, Path]:
    """Build the toy inputs and the upstream databases with plain ``mmseqs`` calls (single thread)."""
    work = Path(directory)
    work.mkdir(parents=True, exist_ok=True)
    world: dict[str, Path] = {"dir": work}

    def path(name: str) -> Path:
        world[name] = work / name
        return world[name]

    write_protein_fasta(path("prot.fa"))
    write_contigs(path("contigs.fa"))
    write_cds_fasta(path("cds.fa"))
    write_gff3(path("features.gff3"))
    write_taxdump(path("taxdump"))
    write_tax_mapping(path("taxmap.tsv"))
    records = contigs()[0]
    names = list(records)
    write_fasta(path("set_a.fa"), {name: records[name] for name in names[:2]})
    write_fasta(path("set_b.fa"), {name: records[name] for name in names[2:]})
    run_direct("createdb", world["prot.fa"], path("protDB"))
    run_direct("createdb", world["contigs.fa"], path("ctgDB"))
    run_direct("createdb", world["cds.fa"], path("cdsDB"), "--shuffle", "0")
    run_direct("createdb", world["prot.fa"], path("prot0DB"), "--shuffle", "0")
    run_direct("prefilter", world["protDB"], world["protDB"], path("pf"))
    run_direct("rescorediagonal", world["protDB"], world["protDB"], world["pf"], path("aln"),
               "--rescore-mode", "2", "-a", "1")  # fmt: skip
    run_direct("clust", world["protDB"], world["aln"], path("clu"))
    run_direct("prefilter", world["prot0DB"], world["prot0DB"], path("pf0"))
    run_direct("rescorediagonal", world["prot0DB"], world["prot0DB"], world["pf0"], path("aln0"),
               "--rescore-mode", "2", "-a", "1")  # fmt: skip
    run_direct(
        "extractorfs", world["ctgDB"], path("orfaa"), "--translate", "1", "--min-length", "30"
    )
    run_direct("createdb", world["prot.fa"], path("refDB"))
    run_direct("createtaxdb", world["refDB"], work / "tmp_taxdb", "--ncbi-tax-dump", world["taxdump"],
               "--tax-mapping-file", world["taxmap.tsv"])  # fmt: skip
    run_direct("prefilter", world["orfaa"], world["refDB"], path("orfpf"))
    run_direct("rescorediagonal", world["orfaa"], world["refDB"], world["orfpf"], path("orfaln"),
               "--rescore-mode", "2", "-a", "1")  # fmt: skip
    run_direct("taxonomy", world["orfaa"], world["refDB"], path("taxDB"), work / "tmp_tax",
               "--tax-output-mode", "2")  # fmt: skip
    run_direct(
        "multihitdb",
        world["set_a.fa"],
        world["set_b.fa"],
        path("sets1"),
        work / "tmp_s1",
        "--min-length",
        "30",
    )
    run_direct(
        "multihitdb",
        world["set_b.fa"],
        world["set_a.fa"],
        path("sets2"),
        work / "tmp_s2",
        "--min-length",
        "30",
    )
    run_direct("search", world["sets1"], world["sets2"], path("setres"), work / "tmp_search")
    run_direct("multihitsearch", world["sets1"], world["sets2"], path("msr"), work / "tmp_mh")
    run_direct("createdb", world["contigs.fa"], path("ntDB"), "--dbtype", "2")
    run_direct("prefilter", world["ntDB"], world["ntDB"], path("ntpf"))
    run_direct("rescorediagonal", world["ntDB"], world["ntDB"], world["ntpf"], path("ntaln"),
               "--rescore-mode", "2", "-a", "1")  # fmt: skip
    return world
