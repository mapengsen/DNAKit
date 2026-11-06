"""Sequence-alignment format engines: PHYLIP, NEXUS (data + taxa + trees blocks) and aligned FASTA."""

# ruff: noqa: RUF012, B905

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from dnakit.io_formats._impl_common import (
    Format,
    Opts,
    Payload,
    bad_input,
    bad_option,
    decode_text,
    json_payload,
    read_source,
)
from dnakit.io_formats._impl_trees import (
    NewickSyntaxError,
    TreeEntry,
    _skip_comment,
    clone,
    nexus_taxa_block,
    nexus_token,
    nexus_trees_block,
    parse_one,
    preorder,
    prune,
    tip_names,
    to_newick,
    tree_stats,
)

_DATATYPES = ("dna", "rna", "protein", "any")
_ALPHABET = {
    "dna": set("ACGTUNRYKMSWBDHV-?.*"),
    "rna": set("ACGTUNRYKMSWBDHV-?.*"),
    "protein": set("ACDEFGHIKLMNPQRSTVWYBZXJUO*-?."),
}
_STATES = {"dna": set("ACGT"), "rna": set("ACGT"), "protein": set("ACDEFGHIKLMNPQRSTVWY")}


@dataclass
class Alignment:
    names: list[str] = field(default_factory=list)
    seqs: list[str] = field(default_factory=list)
    datatype: str = "dna"
    missing: str = "?"
    gap: str = "-"
    trees: list[TreeEntry] = field(default_factory=list)
    taxa: list[str] | None = None
    declared: tuple[int | None, int | None] = (None, None)
    layout: str = "sequential"

    @property
    def nchar(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0


def _missing_set(doc: Alignment) -> set[str]:
    chars = {doc.gap, doc.missing, "-", "?"}
    if doc.datatype in ("dna", "rna"):
        chars |= {"N", "n"}
    elif doc.datatype == "protein":
        chars |= {"X", "x"}
    return chars


def _squash(text: str) -> str:
    return "".join(text.split())


def parse_fasta(text: str, errors: list[str] | None, bad: Any) -> Alignment:
    names, seqs = [], []
    for number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            header = stripped[1:].split(None, 1)
            names.append(header[0] if header else "")
            seqs.append("")
        elif not names:
            bad(errors, f"line {number}", "sequence data before the first '>' header")
            return Alignment()
        else:
            seqs[-1] += _squash(stripped)
    if not names:
        bad(errors, "", "no FASTA records found")
    return Alignment(names, seqs)


def _ranges(spec: object, nchar: int) -> list[tuple[int, int]]:
    """``"1-100,201-300"`` or ``[[1,100],...]`` (1-based inclusive) -> 0-based half-open ranges."""

    pairs: list[tuple[int, int]] = []
    if isinstance(spec, str):
        for part in spec.split(","):
            match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", part)
            if not match:
                raise bad_option(f"Bad site range {part!r}; use 'start-end' (1-based, inclusive).")
            pairs.append((int(match.group(1)), int(match.group(2))))
    elif isinstance(spec, Sequence):
        for item in spec:
            if not (
                isinstance(item, Sequence)
                and len(item) == 2
                and all(isinstance(x, int) for x in item)
            ):
                raise bad_option(
                    "Site ranges must be [start, end] integer pairs (1-based, inclusive)."
                )
            pairs.append((item[0], item[1]))
    else:
        raise bad_option(
            "partition ranges must be a 'start-end[,start-end]' string or a list of [start, end] pairs."
        )
    for start, end in pairs:
        if not 1 <= start <= end <= nchar:
            raise bad_option(f"Site range {start}-{end} is outside the alignment (1-{nchar}).")
    return [(s - 1, e) for s, e in pairs]


class AlignmentFormat(Format):
    noun = "sequences"
    kind = "text"
    targets = ("phylip", "nexus", "fasta", "twobit", "json")
    _DEFAULT_TYPE = "dna"

    #: options shared by every operation that parses the alignment
    _parse_options: dict[str, str] = {}

    # -- construction from records --------------------------------------------------------
    def from_records(self, records: object, opts: Opts) -> Alignment:
        if isinstance(records, Mapping):
            for key in ("records", "sequences", "items"):
                if key in records:
                    return self.from_records(records[key], opts)
            pairs: list[Any] = [(str(k), v) for k, v in records.items()]
        elif isinstance(records, Sequence) and not isinstance(records, (str, bytes)):
            pairs = []
            for number, item in enumerate(records, 1):
                if isinstance(item, Mapping):
                    if "name" not in item or "sequence" not in item:
                        raise bad_input(f"sequence record {number} needs 'name' and 'sequence'.")
                    pairs.append((str(item["name"]), item["sequence"]))
                elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                    pairs.append((str(item[0]), item[1]))
                else:
                    raise bad_input(
                        f"sequence record {number} must be a {{'name','sequence'}} dict or a (name, sequence) pair."
                    )
        else:
            raise bad_input(
                "Expected FASTA text, a name->sequence mapping or a list of sequence records."
            )
        names, seqs = [], []
        for name, seq in pairs:
            if not isinstance(seq, str) or not name:
                raise bad_input(
                    f"Sequence {name!r} must be a non-empty name with a string sequence."
                )
            names.append(name)
            seqs.append(_squash(seq))
        return Alignment(names, seqs, datatype=self._datatype(opts))

    def _datatype(self, opts: Opts) -> str:
        return str(opts.text("datatype", self._DEFAULT_TYPE, _DATATYPES))

    def check_alphabet(self, doc: Alignment) -> None:
        """Raise InputFormatError when a sequence has characters outside the alphabet of ``doc.datatype``."""

        allowed = _ALPHABET.get(doc.datatype)
        if allowed is None:
            return
        allowed = allowed | {doc.missing, doc.gap}
        offending: dict[str, list[str]] = {}
        for name, seq in zip(doc.names, doc.seqs):
            for char in sorted({c for c in seq.upper() if c not in allowed}):
                offending.setdefault(char, []).append(name)
        if offending:
            listed = ", ".join(
                f"{char!r} (in {', '.join(map(repr, names[:3]))}{', ...' if len(names) > 3 else ''})"
                for char, names in sorted(offending.items())
            )
            raise bad_input(
                f"Characters outside the '{doc.datatype}' alphabet: {listed}.",
                code="IO_FORMATS_BAD_ALPHABET",
                hint="Fix the sequences or pass options={'datatype': 'protein'} (or 'any' to skip the check).",
            )

    def coerce_for_write(self, source: object, opts: Opts) -> Alignment:
        """Accept records, FASTA, PHYLIP or NEXUS content (sniffed) as the source of a conversion."""

        if isinstance(source, (Mapping, list, tuple)):
            return self.from_records(source, opts)
        data, _ = read_source(
            source, opts.integer("max_bytes", 256 * 1024 * 1024, 1) or 256 * 1024 * 1024
        )
        head = data.lstrip()[:16].lower()
        if head.startswith(b"#nexus"):
            return NexusFormat().parse(data, None, opts)
        if head.startswith(b">"):
            doc = parse_fasta(decode_text(data), None, self.bad)
            doc.datatype = self._datatype(opts)
            return doc
        return (self if isinstance(self, PhylipFormat) else PhylipFormat()).parse(data, None, opts)

    # -- records / validation ---------------------------------------------------------------
    def records(self, doc: Alignment) -> list[dict[str, object]]:
        return [{"name": n, "sequence": s} for n, s in zip(doc.names, doc.seqs)]

    def describe(self, doc: Alignment) -> dict[str, object]:
        return {"taxon_count": len(doc.names), "site_count": doc.nchar, "datatype": doc.datatype}

    def count(self, doc: Alignment) -> int:
        return len(doc.names)

    def validate_doc(
        self, doc: Alignment, opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        dups = [n for n, c in Counter(doc.names).items() if c > 1]
        if dups:
            errors.append(f"duplicate taxon names: {dups[:5]}")
        if any(not n for n in doc.names):
            errors.append("empty taxon name")
        lengths = {len(s) for s in doc.seqs}
        if len(lengths) > 1:
            errors.append(
                f"sequences have different lengths {sorted(lengths)}; an alignment needs equal lengths"
            )
        ntax, nchar = doc.declared
        if ntax is not None and ntax != len(doc.names):
            errors.append(f"header declares {ntax} taxa but {len(doc.names)} were found")
        if nchar is not None and doc.seqs and lengths != {nchar}:
            errors.append(
                f"header declares {nchar} characters but sequences have length {sorted(lengths)}"
            )
        allowed = _ALPHABET.get(doc.datatype)
        if allowed is not None:
            for name, seq in zip(doc.names, doc.seqs):
                bad_chars = sorted({c for c in seq.upper() if c not in allowed})
                if bad_chars:
                    errors.append(
                        f"{name}: characters {bad_chars[:8]} are not valid for datatype '{doc.datatype}'"
                        " (set the datatype option if this is a different alphabet)"
                    )
                    break
        if any(len(s) == 0 for s in doc.seqs):
            warnings.append("contains empty sequences")
        if doc.trees:
            known = set(doc.names)
            for i, entry in enumerate(doc.trees, 1):
                extra = {t for t in tip_names(entry.root) if t} - known
                if known and extra:
                    warnings.append(f"tree {i}: tips not in the matrix: {sorted(extra)[:5]}")
        return {
            "taxon_count": len(doc.names),
            "site_count": max(lengths, default=0),
            "tree_count": len(doc.trees),
        }

    # -- filter / sort ----------------------------------------------------------------------
    def _miss_frac(self, doc: Alignment, seq: str) -> float:
        miss = _missing_set(doc)
        return sum(1 for c in seq if c in miss) / len(seq) if seq else 0.0

    def filter(self, doc: Alignment, opts: Opts) -> Alignment:
        keep, drop = opts.name_set("keep"), opts.name_set("exclude")
        pattern = opts.text("name_pattern")
        max_taxon = opts.number("max_missing_fraction", None, 0.0, 1.0)
        max_site = opts.number("max_site_missing_fraction", None, 0.0, 1.0)
        rx = None
        if pattern:
            if len(pattern) > 256:
                raise bad_option("name_pattern is limited to 256 characters.")
            try:
                rx = re.compile(pattern)
            except re.error as exc:
                raise bad_option(f"name_pattern is not a valid regular expression: {exc}") from exc
        rows = []
        for name, seq in zip(doc.names, doc.seqs):
            if (keep is not None and name not in keep) or (drop is not None and name in drop):
                continue
            if rx is not None and not rx.search(name):
                continue
            if max_taxon is not None and self._miss_frac(doc, seq) > max_taxon:
                continue
            rows.append((name, seq))
        if max_site is not None and rows:
            miss = _missing_set(doc)
            width = len(rows[0][1])
            good = [
                j
                for j in range(width)
                if sum(1 for _, s in rows if s[j] in miss) / len(rows) <= max_site
            ]
            rows = [(n, "".join(s[j] for j in good)) for n, s in rows]
        return self._rebuild(doc, rows)

    def _rebuild(self, doc: Alignment, rows: list[tuple[str, str]]) -> Alignment:
        kept = {n for n, _ in rows}
        removed = set(doc.names) - kept
        trees = []
        for entry in doc.trees:
            tips = {t for t in tip_names(entry.root) if t}
            root = prune(clone(entry.root), tips - removed) if removed & tips else entry.root
            if root is not None:
                trees.append(TreeEntry(root, entry.name, entry.rooted))
        taxa = None if doc.taxa is None else [t for t in doc.taxa if t not in removed]
        return Alignment(
            [n for n, _ in rows],
            [s for _, s in rows],
            doc.datatype,
            doc.missing,
            doc.gap,
            trees,
            taxa,
            layout=doc.layout,
        )

    def sort(self, doc: Alignment, opts: Opts) -> Alignment:
        key = opts.text("key", "name", ["name", "missing"])
        reverse = opts.boolean("reverse", False)
        rows = list(zip(doc.names, doc.seqs))
        if key == "name":
            rows.sort(key=lambda r: r[0], reverse=reverse)
        else:
            rows.sort(key=lambda r: (self._miss_frac(doc, r[1]), r[0]), reverse=reverse)
        return self._rebuild(doc, rows)

    # -- merge / split ----------------------------------------------------------------------
    def merge(self, docs: list[Alignment], opts: Opts) -> Alignment:
        mode = opts.text("mode", "concatenate", ["concatenate", "append"])
        first = docs[0]
        if any(d.datatype != first.datatype for d in docs):
            raise bad_input("Cannot merge alignments with different datatypes.")
        trees: list[TreeEntry] = []
        used: set[str] = set()
        for d in docs:
            for entry in d.trees:
                name = entry.name or f"tree{len(trees) + 1}"
                base, k = name, 2
                while name in used:
                    name, k = f"{base}_{k}", k + 1
                used.add(name)
                trees.append(TreeEntry(entry.root, name, entry.rooted))
        if mode == "append":
            names: list[str] = []
            seqs: list[str] = []
            for d in docs:
                if d.seqs and first.seqs and d.nchar != first.nchar:
                    raise bad_input(
                        f"append needs equal site counts ({d.nchar} vs {first.nchar}); use mode='concatenate'."
                    )
                for n, s in zip(d.names, d.seqs):
                    if n in names:
                        raise bad_input(
                            f"Duplicate taxon {n!r} while appending alignments.",
                            code="IO_FORMATS_DUPLICATE_ID",
                            hint="Use mode='concatenate' to join sites of the same taxa.",
                        )
                    names.append(n)
                    seqs.append(s)
        else:
            names = list(dict.fromkeys(n for d in docs for n in d.names))
            seqs = []
            for n in names:
                parts = []
                for d in docs:
                    if n in d.names:
                        parts.append(d.seqs[d.names.index(n)])
                    else:
                        parts.append(first.missing * d.nchar)
                seqs.append("".join(parts))
        return Alignment(names, seqs, first.datatype, first.missing, first.gap, trees, None)

    def split(self, doc: Alignment, opts: Opts) -> list[tuple[str, Alignment]]:
        site_size = opts.integer("site_size", None, 1)
        taxa_size = opts.integer("taxa_size", None, 1)
        partitions = opts.mapping("partitions")
        trees_per = opts.integer("trees_per_part", None, 1)
        if trees_per is not None and self.name == "phylip":
            raise bad_option(
                "PHYLIP has no trees, so trees_per_part cannot be used.",
                hint="Give exactly one of site_size, taxa_size or partitions "
                "(trees_per_part applies to NEXUS files with a TREES block).",
            )
        chosen = [x for x in (site_size, taxa_size, partitions, trees_per) if x is not None]
        if len(chosen) != 1:
            raise bad_option(
                "Give exactly one of site_size, taxa_size, partitions or trees_per_part."
            )

        def part(names: list[str], seqs: list[str], trees: list[TreeEntry]) -> Alignment:
            return Alignment(names, seqs, doc.datatype, doc.missing, doc.gap, trees)

        if site_size is not None:
            return [
                (
                    f"sites_{i + 1}-{min(i + site_size, doc.nchar)}",
                    part(doc.names, [s[i : i + site_size] for s in doc.seqs], []),
                )
                for i in range(0, doc.nchar, site_size)
            ]
        if taxa_size is not None:
            return [
                (
                    f"taxa_{i // taxa_size + 1:04d}",
                    part(doc.names[i : i + taxa_size], doc.seqs[i : i + taxa_size], []),
                )
                for i in range(0, len(doc.names), taxa_size)
            ]
        if partitions is not None:
            out = []
            for label, spec in partitions.items():
                spans = _ranges(spec, doc.nchar)
                out.append(
                    (
                        label,
                        part(doc.names, ["".join(s[a:b] for a, b in spans) for s in doc.seqs], []),
                    )
                )
            return out
        assert trees_per is not None
        if not doc.trees:
            raise bad_option("trees_per_part needs a document with trees.")
        return [
            (f"trees_{i // trees_per + 1:04d}", part([], [], doc.trees[i : i + trees_per]))
            for i in range(0, len(doc.trees), trees_per)
        ]

    # -- summaries --------------------------------------------------------------------------
    def summarize(self, doc: Alignment, opts: Opts) -> dict[str, object]:
        miss = _missing_set(doc)
        states = _STATES.get(doc.datatype)
        upper = [s.upper() for s in doc.seqs]
        constant = variable = informative = all_missing = 0
        norm = (
            (lambda c: "T" if c == "U" else c) if doc.datatype in ("dna", "rna") else (lambda c: c)
        )
        for column in zip(*upper) if upper else ():
            kinds = Counter(
                norm(c) for c in column if c not in miss and (states is None or norm(c) in states)
            )
            if not kinds:
                all_missing += 1
            elif len(kinds) == 1:
                constant += 1
            else:
                variable += 1
                if sum(1 for v in kinds.values() if v >= 2) >= 2:
                    informative += 1
        composition = Counter(c for s in upper for c in s)
        total = sum(composition.values())
        fractions = [self._miss_frac(doc, s) for s in doc.seqs]
        out: dict[str, object] = {
            "taxon_count": len(doc.names),
            "site_count": doc.nchar,
            "datatype": doc.datatype,
            "equal_length": len({len(s) for s in doc.seqs}) <= 1,
            "constant_sites": constant,
            "variable_sites": variable,
            "parsimony_informative_sites": informative,
            "all_missing_sites": all_missing,
            "total_characters": total,
            "missing_fraction": sum(composition[c] for c in miss if c in composition) / total
            if total
            else None,
            "taxon_missing_fraction": {
                "min": min(fractions, default=None),
                "mean": sum(fractions) / len(fractions) if fractions else None,
                "max": max(fractions, default=None),
            },
            "duplicate_sequence_count": len(upper) - len(set(upper)),
            "composition": dict(sorted(composition.items())),
        }
        if doc.datatype in ("dna", "rna"):
            gc = composition["G"] + composition["C"]
            acgt = gc + composition["A"] + composition["T"] + composition["U"]
            out["gc_content"] = gc / acgt if acgt else None
        if doc.trees:
            out["tree_count"] = len(doc.trees)
            out["trees"] = [tree_stats(e.root) for e in doc.trees]
        return out

    def index(self, doc: Alignment, opts: Opts) -> dict[str, object]:
        miss = _missing_set(doc)
        entries = {
            n: {"ordinal": i, "length": len(s), "residues": sum(1 for c in s if c not in miss)}
            for i, (n, s) in enumerate(zip(doc.names, doc.seqs))
        }
        out: dict[str, object] = {
            "key": "taxon name",
            "record_count": len(entries),
            "entries": entries,
        }
        if doc.trees:
            out["trees"] = {(e.name or f"tree{i + 1}"): i for i, e in enumerate(doc.trees)}
        return out

    # -- conversion -------------------------------------------------------------------------
    def dump_fasta(self, doc: Alignment, opts: Opts) -> Payload:
        width = opts.integer("line_width", 60, 0) or 0
        lines: list[str] = []
        for name, seq in zip(doc.names, doc.seqs):
            lines.append(f">{name}")
            if width:
                lines.extend(seq[i : i + width] for i in range(0, len(seq), width))
            else:
                lines.append(seq)
        return Payload("fasta", text="".join(line + "\n" for line in lines), count=len(doc.names))

    def dump_phylip(self, doc: Alignment, opts: Opts) -> Payload:
        layout = str(opts.text("layout", "sequential", ["sequential", "interleaved"]))
        block = opts.integer("block_size", 60, 1) or 60
        strict = opts.boolean("strict_names", False)
        if len({len(s) for s in doc.seqs}) > 1:
            raise bad_input(
                "PHYLIP needs equal-length sequences.", hint="Align the sequences first."
            )
        if not doc.names:
            raise bad_input("PHYLIP needs at least one sequence.")
        if strict:
            if any(len(n) > 10 for n in doc.names):
                raise bad_input(
                    "strict_names: names longer than 10 characters cannot be written.",
                    hint="Shorten the names or use strict_names=False (relaxed PHYLIP).",
                )
            width = 10
        else:
            if any(any(c.isspace() for c in n) or not n for n in doc.names):
                raise bad_input(
                    "Relaxed PHYLIP names must be non-empty and contain no whitespace.",
                    hint="Rename the taxa, or pass strict_names=True (names of at most 10 characters may contain spaces).",
                )
            width = max(10, max(len(n) for n in doc.names)) + 1
        nchar = doc.nchar
        lines = [f" {len(doc.names)} {nchar}"]
        if layout == "sequential" or nchar <= block:
            lines += [f"{n:<{width}}{s}" for n, s in zip(doc.names, doc.seqs)]
        else:
            for start in range(0, nchar, block):
                for n, s in zip(doc.names, doc.seqs):
                    lines.append((f"{n:<{width}}" if start == 0 else "") + s[start : start + block])
                if start + block < nchar:
                    lines.append("")
        return Payload("phylip", text="".join(line + "\n" for line in lines), count=len(doc.names))

    def dump_nexus(self, doc: Alignment, opts: Opts) -> Payload:
        layout = str(opts.text("layout", "sequential", ["sequential", "interleaved"]))
        block = opts.integer("block_size", 60, 1) or 60
        out = ["#NEXUS", ""]
        if doc.names:
            if len({len(s) for s in doc.seqs}) > 1:
                raise bad_input(
                    "NEXUS DATA needs equal-length sequences.", hint="Align the sequences first."
                )
            nexus_type = {"dna": "DNA", "rna": "RNA", "protein": "PROTEIN"}.get(
                doc.datatype, "STANDARD"
            )
            tokens = [nexus_token(n) for n in doc.names]
            width = max(len(t) for t in tokens) + 2
            interleaved = layout == "interleaved" and doc.nchar > block
            out += [
                "BEGIN DATA;",
                f"    DIMENSIONS NTAX={len(doc.names)} NCHAR={doc.nchar};",
                f"    FORMAT DATATYPE={nexus_type} MISSING={doc.missing} GAP={doc.gap}"
                + (" INTERLEAVE=YES;" if interleaved else ";"),
                "    MATRIX",
            ]
            if interleaved:
                for start in range(0, doc.nchar, block):
                    out += [
                        f"    {t:<{width}}{s[start : start + block]}"
                        for t, s in zip(tokens, doc.seqs)
                    ]
                    if start + block < doc.nchar:
                        out.append("")
            else:
                out += [f"    {t:<{width}}{s}" for t, s in zip(tokens, doc.seqs)]
            out += ["    ;", "END;", ""]
        elif doc.trees:
            labels = (
                list(dict.fromkeys(t for e in doc.trees for t in tip_names(e.root) if t))
                if doc.taxa is None
                else doc.taxa
            )
            out += [nexus_taxa_block(labels)]
        if doc.trees:
            out.append(nexus_trees_block(doc.trees))
        return Payload("nexus", text="\n".join(out), count=len(doc.names) + len(doc.trees))

    def _convert_matrix(self, doc: Alignment, target: str, opts: Opts) -> Payload:
        if target == "phylip":
            return self.dump_phylip(doc, opts)
        if target == "fasta":
            return self.dump_fasta(doc, opts)
        from dnakit.io_formats._impl_twobit import TwoBitFormat

        engine = TwoBitFormat()
        return engine.dump(engine.from_records(AlignmentFormat.records(self, doc), opts), opts)

    def convert(self, doc: Alignment, target: str, opts: Opts) -> Payload:
        if target in ("phylip", "fasta", "twobit"):
            payload = self._convert_matrix(doc, target, opts)
            if doc.trees:
                count = len(doc.trees)
                payload.info["warnings"] = [
                    f"The {count} tree(s) of the TREES block were dropped: the {target} format holds only the matrix."
                ]
            return payload
        if target == "nexus":
            return self.dump_nexus(doc, opts)
        if target == "newick":
            if not doc.trees:
                raise bad_input("This document has no trees to convert to Newick.")
            return Payload(
                "newick",
                text="".join(to_newick(e.root) + "\n" for e in doc.trees),
                count=len(doc.trees),
            )
        info: dict[str, object] = {
            "taxon_count": len(doc.names),
            "site_count": doc.nchar,
            "datatype": doc.datatype,
        }
        if doc.trees:
            info["trees"] = [
                {"name": e.name, "rooted": e.rooted, "newick": to_newick(e.root)} for e in doc.trees
            ]
        return json_payload(self.name, {**info, "records": self.records(doc)}, len(doc.names))


# =============================================================================================
# PHYLIP
# =============================================================================================

_ALL_OPS = (
    "read",
    "write",
    "validate",
    "filter",
    "sort",
    "merge",
    "split",
    "stream",
    "summarize",
    "index",
    "convert",
)
_ALIGN_FILTER = {
    "keep": "list of taxon names to keep",
    "exclude": "list of taxon names to drop",
    "name_pattern": "regular expression (re.search, max 256 characters) a taxon name must match",
    "max_missing_fraction": "float in [0,1]: drop taxa whose gap/missing/ambiguous fraction exceeds it",
    "max_site_missing_fraction": "float in [0,1]: drop alignment columns whose missing fraction (over the kept taxa) exceeds it",
}
_ALIGN_CONVERT = {
    "layout": "'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output",
    "block_size": "int, default 60: sites per block of an interleaved output",
    "strict_names": "bool, default False: write strict 10-character PHYLIP names (error when longer)",
    "line_width": "int, default 60: FASTA line width (0 = one line per sequence)",
}
_ALIGN_MERGE = {
    "mode": "'concatenate' (default: join the sites of the same taxa into a supermatrix, filling taxa absent from "
    "an input with the missing character) or 'append' (add taxa; site counts must be equal)"
}
_ALIGN_SPLIT = {
    "site_size": "int: consecutive blocks of this many sites",
    "taxa_size": "int: consecutive groups of this many taxa",
    "partitions": "mapping label -> 'start-end[,start-end...]' (1-based, inclusive site ranges)",
}
_NEXUS_SPLIT = {
    **_ALIGN_SPLIT,
    "trees_per_part": "int: trees per part (documents without a matrix); give exactly one of the four split options",
}


class PhylipFormat(AlignmentFormat):
    name = "phylip"
    title = "PHYLIP alignment file"
    suffix = ".phy"
    targets = ("phylip", "nexus", "fasta", "twobit", "json")
    options = {
        **{
            op: {
                "datatype": "'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for "
                "counting missing characters (gap, ?, N for nucleotides, X for protein)"
            }
            for op in _ALL_OPS
        },
    }
    for _op in (
        "read",
        "validate",
        "filter",
        "sort",
        "merge",
        "split",
        "stream",
        "summarize",
        "index",
        "convert",
    ):
        options[_op] = {
            **options[_op],
            "input_layout": "'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input",
            "input_strict_names": "bool, default False: input names occupy the first 10 characters of each row "
            "(otherwise the first whitespace-delimited word)",
        }
    del _op
    options["filter"] = {**options["filter"], **_ALIGN_FILTER}
    options["sort"] = {
        **options["sort"],
        "key": "'name' (default) or 'missing' (missing fraction, then name)",
        "reverse": "bool, default False",
    }
    options["merge"] = {**options["merge"], **_ALIGN_MERGE}
    options["split"] = {**options["split"], **_ALIGN_SPLIT}
    rejected_options = {"split": ("trees_per_part",)}
    options["convert"] = {**options["convert"], **_ALIGN_CONVERT}
    options["write"] = {
        **options["write"],
        **_ALIGN_CONVERT,
        "input_layout": "'auto' (default), 'sequential' or 'interleaved' when the input is PHYLIP text",
        "input_strict_names": "bool, default False when the input is PHYLIP text",
    }
    semantics = {
        "read": "Reads relaxed PHYLIP (name = first whitespace-delimited word; blocks of residues may contain spaces) or strict "
        "PHYLIP; sequential and interleaved layouts are detected automatically. The header 'ntax nchar' must match.",
        "validate": "Checks the header counts, equal sequence lengths, unique taxon names and (by datatype) the alphabet.",
        "summarize": "Taxa, sites, constant/variable/parsimony-informative sites (a site is informative when at least two "
        "different unambiguous states each occur in at least two taxa; gaps and ambiguity codes are ignored), "
        "composition, GC content (nucleotides), per-taxon missing fractions and duplicate sequences.",
        "index": "Maps each taxon name to its row ordinal, sequence length and number of non-missing residues (PHYLIP rows are "
        "not stored one per line when interleaved, so no byte offsets are given).",
        "convert": "Targets: 'phylip' (re-layout), 'nexus', 'fasta', 'twobit' (nucleotides) and 'json'.",
        "write": "Input may be FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records; "
        "the alignment is written as relaxed PHYLIP (names padded to at least 10 characters) unless strict_names is true.",
    }

    def coerce_for_write(self, source: object, opts: Opts) -> Alignment:
        doc = super().coerce_for_write(source, opts)
        self.check_alphabet(doc)
        return doc

    def _layout_args(self, opts: Opts) -> tuple[str, bool]:
        return str(
            opts.text("input_layout", "auto", ["auto", "sequential", "interleaved"])
        ), opts.boolean("input_strict_names", False)

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> Alignment:
        text = decode_text(data)
        layout, strict = self._layout_args(opts)
        datatype = self._datatype(opts)
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            self.bad(errors, "", "empty input")
            return Alignment(datatype=datatype)
        header = lines[0].split()
        try:
            ntax, nchar = int(header[0]), int(header[1])
            if (
                ntax < 1
                or nchar < 0
                or (len(header) > 2 and not set("".join(header[2:]).upper()) <= set("ISC"))
            ):
                raise ValueError
        except (ValueError, IndexError):
            self.bad(
                errors,
                "line 1",
                f"header must be 'ntax nchar' (two integers), got {lines[0].strip()!r}",
            )
            return Alignment(datatype=datatype)
        flags = "".join(header[2:]).upper()
        if layout == "auto" and "I" in flags:
            layout = "interleaved"
        elif layout == "auto" and "S" in flags:
            layout = "sequential"
        body = lines[1:]

        def split_name(line: str) -> tuple[str, str]:
            if strict:
                return line[:10].strip(), _squash(line[10:])
            parts = line.strip().split(None, 1)
            return parts[0], _squash(parts[1]) if len(parts) > 1 else ""

        def sequential() -> Alignment | None:
            i, names, seqs = 0, [], []
            for _ in range(ntax):
                if i >= len(body):
                    return None
                name, seq = split_name(body[i])
                i += 1
                while len(seq) < nchar and i < len(body):
                    seq += _squash(body[i])
                    i += 1
                if len(seq) != nchar:
                    return None
                names.append(name)
                seqs.append(seq)
            return (
                Alignment(names, seqs, datatype, declared=(ntax, nchar), layout="sequential")
                if i == len(body)
                else None
            )

        def interleaved() -> Alignment | None:
            if len(body) < ntax or (len(body) - ntax) % ntax:
                return None
            pairs = [split_name(line) for line in body[:ntax]]
            names, seqs = [p[0] for p in pairs], [p[1] for p in pairs]
            for j, line in enumerate(body[ntax:]):
                seqs[j % ntax] += _squash(line)
            if any(len(s) != nchar for s in seqs):
                return None
            return Alignment(names, seqs, datatype, declared=(ntax, nchar), layout="interleaved")

        attempts = {
            "sequential": [sequential],
            "interleaved": [interleaved],
            "auto": [sequential, interleaved],
        }[layout]
        for attempt in attempts:
            result = attempt()
            if result is not None:
                return result
        self.bad(
            errors,
            "",
            f"the sequence rows do not form {ntax} taxa of exactly {nchar} characters in "
            f"{'a sequential or interleaved' if layout == 'auto' else layout} layout",
        )
        return Alignment(datatype=datatype, declared=(ntax, nchar))

    def dump(self, doc: Alignment, opts: Opts) -> Payload:
        return self.dump_phylip(doc, opts)


# =============================================================================================
# aligned FASTA (conversion source / target only)
# =============================================================================================


class FastaAlignmentFormat(AlignmentFormat):
    name = "fasta"
    title = "FASTA file"
    suffix = ".fasta"
    options = {op: {"datatype": "'dna' (default), 'rna', 'protein' or 'any'"} for op in _ALL_OPS}
    options["convert"] = {**options["convert"], **_ALIGN_CONVERT}
    targets = ("phylip", "nexus", "fasta", "twobit", "json")
    semantics: dict[str, str] = {}

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> Alignment:
        doc = parse_fasta(decode_text(data), errors, self.bad)
        doc.datatype = self._datatype(opts)
        return doc

    def dump(self, doc: Alignment, opts: Opts) -> Payload:
        return self.dump_fasta(doc, opts)


# =============================================================================================
# NEXUS
# =============================================================================================


def _scan_commands(text: str) -> list[str]:
    """Split NEXUS text into ';'-terminated commands; comments are removed (except [&R]/[&U] -> \\x01R/\\x01U)."""

    out: list[str] = []
    buf: list[str] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "'":
            j = i + 1
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        j += 2
                        continue
                    break
                j += 1
            if j >= n:
                raise NewickSyntaxError("unterminated quoted token")
            buf.append(text[i : j + 1])
            i = j + 1
        elif c == "[":
            j = _skip_comment(text, i)
            inner = text[i + 1 : j - 1].strip().upper()
            buf.append("\x01" + inner[1] if inner in ("&R", "&U") else " ")
            i = j
        elif c == ";":
            out.append("".join(buf).strip())
            buf = []
            i += 1
        else:
            buf.append(c)
            i += 1
    if "".join(buf).strip():
        raise NewickSyntaxError("the last command is not terminated by ';'")
    return out


_TOKEN = re.compile(r"'(?:[^']|'')*'|[^\s]+")


def _unquote(token: str) -> str:
    return (
        token[1:-1].replace("''", "'")
        if token.startswith("'") and token.endswith("'") and len(token) >= 2
        else token
    )


def _kv(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z_]+)\s*=\s*("[^"]*"|\'[^\']*\'|[^\s;]+)|([A-Za-z_]+)', text):
        if match.group(1):
            out[match.group(1).lower()] = match.group(2).strip("\"'")
        else:
            out.setdefault(match.group(3).lower(), "yes")
    return out


class NexusFormat(AlignmentFormat):
    name = "nexus"
    title = "NEXUS file"
    suffix = ".nex"
    targets = ("nexus", "phylip", "fasta", "twobit", "newick", "json")
    options = {op: {} for op in _ALL_OPS}
    options["filter"] = {**_ALIGN_FILTER}
    options["sort"] = {"key": "'name' (default) or 'missing'", "reverse": "bool, default False"}
    options["merge"] = {**_ALIGN_MERGE}
    options["split"] = {**_NEXUS_SPLIT}
    options["convert"] = {**_ALIGN_CONVERT}
    options["write"] = {
        **_ALIGN_CONVERT,
        "datatype": "'dna' (default), 'rna', 'protein' or 'any' for records without a declared type",
    }
    semantics = {
        "read": "Parses the TAXA, DATA/CHARACTERS and TREES blocks (other blocks are ignored): taxon labels, the matrix "
        "(sequential or interleaved, matchchar expanded), the declared datatype/missing/gap symbols and the trees "
        "(TRANSLATE tables applied, [&R]/[&U] kept as the 'rooted' flag). Comments in [] are dropped.",
        "validate": "Checks the #NEXUS header, block structure, declared NTAX/NCHAR against the matrix, equal row lengths, unique "
        "taxa, the alphabet of nucleotide/protein data and (as warnings) tree tips absent from the matrix.",
        "filter": "Taxon filters apply to the matrix, the TAXA labels and (by pruning tips) every tree.",
        "summarize": "Matrix statistics as for PHYLIP (constant/variable/parsimony-informative sites, composition, missing "
        "fractions) plus per-tree statistics when TREES blocks are present.",
        "index": "Maps taxon names to row ordinal, length and non-missing residues, and tree names to their position.",
        "convert": "Targets: 'nexus' (normalised: DATA block + TREES block), 'phylip', 'fasta', 'twobit' (nucleotides), "
        "'newick' (one tree per line; needs trees) and 'json'.",
        "write": "Input may be FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records; "
        "written as a DATA block (names quoted when needed).",
    }

    def _datatype(self, opts: Opts) -> str:
        return str(opts.text("datatype", "dna", _DATATYPES))

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> Alignment:
        text = decode_text(data).lstrip()
        if text[:6].lower() != "#nexus":
            self.bad(errors, "", "the file does not start with '#NEXUS'")
            return Alignment()
        try:
            commands = _scan_commands(text[6:])
        except NewickSyntaxError as exc:
            self.bad(errors, "", str(exc))
            return Alignment()
        doc = Alignment(datatype="any")
        block: str | None = None
        declared_ntax = declared_nchar = None
        taxlabels: list[str] | None = None
        translate: dict[str, str] = {}
        interleave = False
        matchchar = None
        matrix_text = None
        for cmd in commands:
            if not cmd:
                continue
            parts = cmd.split(None, 1)
            word, rest = parts[0].lower(), parts[1] if len(parts) > 1 else ""
            if word == "begin":
                block = rest.strip().lower()
                continue
            if word in ("end", "endblock"):
                block = None
                continue
            if block is None:
                continue
            if word == "dimensions":
                values = _kv(rest)
                if values.get("ntax", "").isdigit():
                    declared_ntax = int(values["ntax"])
                if values.get("nchar", "").isdigit():
                    declared_nchar = int(values["nchar"])
            elif block == "taxa" and word == "taxlabels":
                taxlabels = [_unquote(t) for t in _TOKEN.findall(rest)]
            elif block in ("data", "characters") and word == "format":
                values = _kv(rest)
                dt = values.get("datatype", "standard").lower()
                doc.datatype = {"nucleotide": "dna"}.get(dt, dt)
                doc.missing = values.get("missing", doc.missing)[:1] or "?"
                doc.gap = values.get("gap", doc.gap)[:1] or "-"
                interleave = values.get("interleave", "no").lower() in ("yes", "true")
                matchchar = values.get("matchchar")
            elif block in ("data", "characters") and word == "matrix":
                matrix_text = rest
            elif block == "trees" and word == "translate":
                for pair in rest.split(","):
                    tokens = _TOKEN.findall(pair)
                    if len(tokens) == 2:
                        translate[_unquote(tokens[0])] = _unquote(tokens[1])
            elif block == "trees" and word in ("tree", "utree"):
                match = re.match(r"\*?\s*('(?:[^']|'')*'|[^\s=]+)\s*=\s*(.*)\Z", rest, re.S)
                if not match:
                    self.bad(errors, "TREES block", f"cannot read tree command {cmd[:40]!r}")
                    continue
                body = match.group(2).strip()
                rooted: bool | None = None
                if body.startswith("\x01"):
                    rooted = body[1] == "R"
                    body = body[2:]
                body = body.replace("\x01R", "").replace("\x01U", "")
                try:
                    root, _ = parse_one(body + ";", 0)
                except NewickSyntaxError as exc:
                    self.bad(errors, f"tree {_unquote(match.group(1))}", str(exc))
                    continue
                if translate:
                    for node in preorder(root):
                        if node.name in translate:
                            node.name = translate[node.name]
                doc.trees.append(TreeEntry(root, _unquote(match.group(1)), rooted))
        if matrix_text is not None:
            names: list[str] = []
            rows: dict[str, str] = {}
            for line in matrix_text.splitlines():
                line = line.strip()
                if not line:
                    continue
                match = re.match(r"('(?:[^']|'')*'|\S+)\s*(.*)\Z", line)
                if not match:
                    continue
                name = _unquote(match.group(1))
                if name not in rows:
                    names.append(name)
                    rows[name] = ""
                rows[name] += _squash(match.group(2))
            doc.names = names
            doc.seqs = [rows[n] for n in names]
            if matchchar and doc.seqs:
                ref = doc.seqs[0]
                doc.seqs = [doc.seqs[0]] + [
                    "".join(
                        ref[j] if c == matchchar and j < len(ref) else c for j, c in enumerate(s)
                    )
                    for s in doc.seqs[1:]
                ]
            doc.layout = "interleaved" if interleave else "sequential"
        doc.taxa = taxlabels if taxlabels is not None else (list(doc.names) if doc.names else None)
        doc.declared = (declared_ntax if doc.names else None, declared_nchar if doc.names else None)
        if doc.datatype == "any" and doc.names:
            doc.datatype = "standard"
        if not doc.names and not doc.trees and taxlabels is None and errors is None:
            raise bad_input(
                "The NEXUS file contains no TAXA, DATA/CHARACTERS or TREES block.",
                code="IO_FORMATS_BAD_NEXUS",
            )
        return doc

    def records(self, doc: Alignment) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = [
            {"kind": "sequence", "name": n, "sequence": s} for n, s in zip(doc.names, doc.seqs)
        ]
        rows += [
            {"kind": "tree", "name": e.name, "rooted": e.rooted, "newick": to_newick(e.root)}
            for e in doc.trees
        ]
        return rows

    def count(self, doc: Alignment) -> int:
        return len(doc.names) + len(doc.trees)

    def describe(self, doc: Alignment) -> dict[str, object]:
        return {
            **super().describe(doc),
            "taxa": doc.taxa,
            "tree_count": len(doc.trees),
            "missing": doc.missing,
            "gap": doc.gap,
            "layout": doc.layout,
        }

    def validate_doc(
        self, doc: Alignment, opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        if doc.taxa is not None and doc.names and set(doc.taxa) != set(doc.names):
            errors.append("TAXLABELS and the matrix row names differ")
        return super().validate_doc(doc, opts, errors, warnings)

    def dump(self, doc: Alignment, opts: Opts) -> Payload:
        return self.dump_nexus(doc, opts)

    def from_records(self, records: object, opts: Opts) -> Alignment:
        return super().from_records(records, opts)
