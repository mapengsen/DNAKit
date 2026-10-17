"""Read-level and record-level single-cell operations: barcode/UMI extraction, UMI correction, matrix building."""

from __future__ import annotations

import bisect
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence

from dnakit.exceptions import ConfigurationError, InputFormatError

from ._impl_matrix import (
    CountMatrix,
    DEFAULT_MAX_ELEMENTS,
    matrix_result,
    opt_choice,
    opt_int,
    text_of,
)

Result = tuple[object, list[str]]
MAX_RECORDS = 5_000_000
MAX_LISTED = 1000
_DNA = frozenset("ACGTN")


# ---------------------------------------------------------------------------------------------
# reads
# ---------------------------------------------------------------------------------------------


def read_sequences(value: object) -> list[tuple[str, str]]:
    """Return ``(name, sequence)`` pairs from FASTQ/FASTA text or path, a list of sequences, or ``{name: sequence}``."""

    if value is None:
        raise InputFormatError("input is required.", code="EMPTY_INPUT", hint="Provide reads.")
    pairs: list[tuple[str, str]] = []
    if isinstance(value, Mapping):
        for name, seq in value.items():
            if isinstance(seq, Mapping):
                seq = seq.get("sequence")
            pairs.append((str(name), str(seq)))
    elif isinstance(value, (str, bytes)) or hasattr(value, "__fspath__"):
        text = text_of(value, "reads")
        lines = [line.rstrip("\r") for line in text.splitlines() if line.strip()]
        if lines and lines[0].startswith("@"):
            if len(lines) % 4:
                raise InputFormatError("FASTQ text must have 4 lines per record.", code="FASTQ_FORMAT", hint="Check for truncation.")
            for i in range(0, len(lines), 4):
                if not lines[i].startswith("@") or not lines[i + 2].startswith("+"):
                    raise InputFormatError(f"Malformed FASTQ record at line {i + 1}.", code="FASTQ_FORMAT", hint="Expected @name / sequence / + / quality.")
                pairs.append((lines[i][1:].split()[0] if len(lines[i]) > 1 else "", lines[i + 1].strip()))
        elif lines and lines[0].startswith(">"):
            name, parts = None, []
            for line in lines:
                if line.startswith(">"):
                    if name is not None:
                        pairs.append((name, "".join(parts)))
                    name, parts = line[1:].split()[0] if len(line) > 1 else "", []
                else:
                    parts.append(line.strip())
            if name is not None:
                pairs.append((name, "".join(parts)))
        else:
            pairs = [(f"read_{i}", line.strip()) for i, line in enumerate(lines, start=1)]
    elif isinstance(value, Sequence):
        for i, item in enumerate(value, start=1):
            if isinstance(item, Mapping):
                pairs.append((str(item.get("name", f"read_{i}")), str(item.get("sequence", ""))))
            else:
                pairs.append((f"read_{i}", str(item)))
    else:
        raise InputFormatError(f"Unsupported reads input of type {type(value).__name__}.", code="UNSUPPORTED_INPUT", hint="Use FASTQ/FASTA text, a path, a list of sequences or a mapping.")
    if not pairs:
        raise InputFormatError("No reads found.", code="EMPTY_INPUT", hint="Provide at least one read.")
    if len(pairs) > MAX_RECORDS:
        raise InputFormatError(f"More than {MAX_RECORDS} reads.", code="INPUT_TOO_LARGE", hint="Process the reads in chunks.")
    cleaned = []
    for name, seq in pairs:
        seq = seq.strip().upper()
        if seq and not set(seq) <= _DNA:
            raise InputFormatError(f"Read {name!r} contains non-ACGTN characters.", code="INVALID_SEQUENCE", hint="Reads must be DNA (A, C, G, T, N).")
        cleaned.append((name, seq))
    return cleaned


def _whitelist(raw: object) -> set[str] | None:
    if raw is None:
        return None
    if isinstance(raw, (str, bytes)) or hasattr(raw, "__fspath__"):
        items = [line.split()[0] for line in text_of(raw, "whitelist").splitlines() if line.strip() and not line.startswith("#")]
    elif isinstance(raw, Sequence) or isinstance(raw, (set, frozenset)):
        items = [str(x) for x in raw]
    else:
        raise ConfigurationError("options['whitelist'] must be a list, text or path.", code="INVALID_OPTION", hint="One barcode per line.")
    return {item.strip().upper().split("-")[0] for item in items if item.strip()}


def _correct_barcode(barcode: str, whitelist: set[str]) -> str | None:
    """Exact whitelist hit, or the unique whitelist barcode at Hamming distance 1 (else ``None``)."""

    if barcode in whitelist:
        return barcode
    hits = []
    for i, base in enumerate(barcode):
        for other in "ACGT":
            if other != base:
                candidate = barcode[:i] + other + barcode[i + 1:]
                if candidate in whitelist:
                    hits.append(candidate)
    return hits[0] if len(hits) == 1 else None


def _extract(value: object, options: dict[str, object], *, kind: str) -> Result:
    start_default, length_default = (0, 16) if kind == "barcode" else (16, 12)
    start = opt_int(options, f"{kind}_start", start_default, 0, 10_000)
    length = opt_int(options, f"{kind}_length", length_default, 1, 200)
    top_n = opt_int(options, "top_n", 100, 1, 10**6)
    reads = read_sequences(value)
    whitelist = _whitelist(options.get("whitelist")) if kind == "barcode" else None
    counts: Counter[str] = Counter()
    listed: list[dict[str, object]] = []
    too_short = with_n = exact = corrected = unmatched = 0
    cache: dict[str, str | None] = {}
    for name, seq in reads:
        piece = seq[start:start + length]
        if len(piece) < length:
            too_short += 1
            continue
        record: dict[str, object] = {"read": name, kind: piece}
        if "N" in piece:
            with_n += 1
        if whitelist is not None:
            if piece not in cache:
                cache[piece] = _correct_barcode(piece, whitelist) if "N" not in piece else None
            fixed = cache[piece]
            if fixed is None:
                unmatched += 1
                record["corrected"] = None
                if len(listed) < MAX_LISTED:
                    listed.append(record)
                continue
            if fixed == piece:
                exact += 1
            else:
                corrected += 1
            record["corrected"] = fixed
            piece = fixed
        counts[piece] += 1
        if len(listed) < MAX_LISTED:
            listed.append(record)
    result: dict[str, object] = {
        "definition": f"{kind} = read[{start}:{start + length}] (0-based, {length} bases)" + (" matched to the whitelist exactly or by a unique 1-mismatch correction" if whitelist is not None else ""),
        "n_reads": len(reads),
        "n_extracted": sum(counts.values()),
        **({"n_unmatched": unmatched} if kind == "barcode" else {}),  # whitelist given but no exact / unique 1-mismatch match (0 without a whitelist)
        "n_too_short": too_short,
        f"n_{kind}_with_N": with_n,
        f"n_unique_{kind}s": len(counts),
        f"{kind}_counts": dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:top_n]),
        "reads": listed,
        "reads_truncated": len(reads) - too_short > len(listed),
    }
    if whitelist is not None:
        result["whitelist"] = {"size": len(whitelist), "exact": exact, "corrected": corrected, "unmatched": unmatched}
    return result, []


def op_extract_barcodes(value: object, options: dict[str, object]) -> Result:
    return _extract(value, options, kind="barcode")


def op_extract_umis(value: object, options: dict[str, object]) -> Result:
    return _extract(value, options, kind="umi")


# ---------------------------------------------------------------------------------------------
# UMI correction (UMI-tools algorithms)
# ---------------------------------------------------------------------------------------------


def _umi_counts(value: object) -> Counter[str]:
    if value is None:
        raise InputFormatError("input is required.", code="EMPTY_INPUT", hint="Provide UMIs.")
    counts: Counter[str] = Counter()
    if isinstance(value, Mapping):
        for umi, n in value.items():
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                raise InputFormatError(f"Count of {umi!r} must be a positive integer.", code="INVALID_COUNT", hint="Use {umi: reads}.")
            counts[str(umi).strip().upper()] += n
    else:
        if isinstance(value, (str, bytes)) or hasattr(value, "__fspath__"):
            items: list[object] = [line for line in text_of(value, "UMIs").splitlines() if line.strip() and not line.startswith("#")]
        elif isinstance(value, Sequence):
            items = list(value)
        else:
            raise InputFormatError(f"Unsupported UMI input of type {type(value).__name__}.", code="UNSUPPORTED_INPUT", hint="Use a list, mapping or text.")
        for item in items:
            if isinstance(item, str) and ("\t" in item or "," in item or " " in item.strip()):
                parts = item.replace(",", "\t").split()
                try:
                    counts[parts[0].upper()] += int(parts[1])
                except (IndexError, ValueError):
                    raise InputFormatError(f"Bad UMI line {item!r}.", code="UMI_FORMAT", hint="Use 'UMI' or 'UMI<TAB>count'.") from None
            else:
                counts[str(item).strip().upper()] += 1
    if not counts:
        raise InputFormatError("No UMIs found.", code="EMPTY_INPUT", hint="Provide at least one UMI.")
    lengths = {len(u) for u in counts}
    if len(lengths) != 1:
        raise InputFormatError(f"UMIs have different lengths: {sorted(lengths)}.", code="UMI_LENGTH", hint="All UMIs must have the same length.")
    for umi in counts:
        if not umi or not set(umi) <= _DNA:
            raise InputFormatError(f"Invalid UMI {umi!r}.", code="INVALID_UMI", hint="UMIs use A, C, G, T or N.")
    return counts


def _neighbours(umis: Sequence[str], threshold: int) -> dict[str, list[str]]:
    present = set(umis)
    adjacency: dict[str, list[str]] = {u: [] for u in umis}
    if threshold == 1:
        for umi in umis:
            for i, base in enumerate(umi):
                for other in "ACGTN":
                    if other != base:
                        candidate = umi[:i] + other + umi[i + 1:]
                        if candidate in present:
                            adjacency[umi].append(candidate)
        return adjacency
    if len(umis) > 5000:
        raise InputFormatError("threshold > 1 supports up to 5000 distinct UMIs.", code="INPUT_TOO_LARGE", hint="Use threshold=1 or split by gene/cell first.")
    for i, a in enumerate(umis):
        for b in umis[i + 1:]:
            if sum(x != y for x, y in zip(a, b)) <= threshold:
                adjacency[a].append(b)
                adjacency[b].append(a)
    return adjacency


def group_umis(counts: Mapping[str, int], method: str, threshold: int) -> list[list[str]]:
    """Group UMIs; the first member of each group is its representative (most abundant, ties by sequence)."""

    order = sorted(counts, key=lambda u: (-counts[u], u))
    position = {u: i for i, u in enumerate(order)}
    if method == "unique":
        return [[u] for u in order]
    adjacency = {u: sorted(v, key=position.get) for u, v in _neighbours(order, threshold).items()}  # type: ignore[arg-type]
    seen: set[str] = set()
    groups: list[list[str]] = []
    for top in order:
        if top in seen:
            continue
        seen.add(top)
        group = [top]
        if method == "adjacency":
            for other in adjacency[top]:
                if other not in seen:
                    seen.add(other)
                    group.append(other)
        else:
            frontier = [top]
            while frontier:
                node = frontier.pop(0)
                for other in adjacency[node]:
                    if other in seen:
                        continue
                    if method == "cluster" or counts[node] >= 2 * counts[other] - 1:
                        seen.add(other)
                        group.append(other)
                        frontier.append(other)
        groups.append(group)
    return groups


def op_correct_umis(value: object, options: dict[str, object]) -> Result:
    counts = _umi_counts(value)
    method = opt_choice(options, "method", "directional", ("directional", "adjacency", "cluster", "unique"))
    threshold = opt_int(options, "threshold", 1, 1, 4)
    groups = group_umis(counts, method, threshold)  # type: ignore[arg-type]
    mapping = {member: group[0] for group in groups for member in group}
    listed = [
        {"representative": group[0], "members": group, "reads": sum(counts[m] for m in group)} for group in groups[:MAX_LISTED]
    ]
    definitions = {
        "directional": "UMI-tools directional: starting from the most abundant UMI, absorb neighbours within the Hamming threshold whose count n_b satisfies n_a >= 2*n_b - 1 (followed transitively)",
        "adjacency": "UMI-tools adjacency: most abundant UMI absorbs all unassigned neighbours within the Hamming threshold (one step)",
        "cluster": "UMI-tools cluster: connected components of the Hamming-threshold graph",
        "unique": "no correction: every distinct UMI is its own group",
    }
    return {
        "definition": definitions[method], "method": method, "threshold": threshold,
        "n_reads": sum(counts.values()), "n_unique_before": len(counts), "n_unique_after": len(groups),
        "corrected": {u: r for u, r in mapping.items() if u != r} if len(mapping) <= 50_000 else None,
        "groups": listed, "groups_truncated": len(groups) > MAX_LISTED,
    }, []


# ---------------------------------------------------------------------------------------------
# building matrices from records
# ---------------------------------------------------------------------------------------------

_CELL_KEYS = ("cell", "barcode", "cell_barcode", "cb", "cell_id")
_GENE_KEYS = ("gene", "feature", "gene_id", "gene_name", "gx", "gn", "peak")
_UMI_KEYS = ("umi", "ub", "umi_tag")


def _rows(value: object, what: str) -> tuple[list[str] | None, list[list[object]]]:
    """Return an optional header and the data rows of text, a path, a list of lists or a list of mappings."""

    if value is None:
        raise InputFormatError("input is required.", code="EMPTY_INPUT", hint=f"Provide {what}.")
    if isinstance(value, (str, bytes)) or hasattr(value, "__fspath__"):
        rows: list[list[object]] = []
        for line in text_of(value, what).splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else (line.split(",") if "," in line else line.split())
            rows.append([p.strip() for p in parts])
        header = None
        if rows and str(rows[0][0]).casefold() in _CELL_KEYS + ("chrom", "chr", "chromosome"):
            header = [str(h).casefold() for h in rows[0]]
            rows = rows[1:]
        return header, rows
    if isinstance(value, Sequence):
        items = list(value)
        if items and all(isinstance(i, Mapping) for i in items):
            keys = [str(k).casefold() for k in items[0]]
            return keys, [[{str(k).casefold(): v for k, v in i.items()}.get(key) for key in keys] for i in items]
        return None, [list(i) for i in items]  # type: ignore[call-overload]
    raise InputFormatError(f"Unsupported {what} input of type {type(value).__name__}.", code="UNSUPPORTED_INPUT", hint="Use delimited text, a path, a list of lists or a list of mappings.")


def _column(header: list[str] | None, names: Sequence[str], default: int | None) -> int | None:
    if header:
        for name in names:
            if name in header:
                return header.index(name)
    return default


def op_build_cell_gene_matrix(value: object, options: dict[str, object]) -> Result:
    header, rows = _rows(value, "cell-gene records")
    if len(rows) > MAX_RECORDS:
        raise InputFormatError(f"More than {MAX_RECORDS} records.", code="INPUT_TOO_LARGE", hint="Process in chunks.")
    c_i = _column(header, _CELL_KEYS, 0)
    g_i = _column(header, _GENE_KEYS, 1)
    u_i = _column(header, _UMI_KEYS, 2 if header is None else None)
    n_i = _column(header, ("count", "n", "reads"), None)
    if c_i is None or g_i is None:
        raise InputFormatError("Records need a cell and a gene column.", code="RECORD_COLUMNS", hint="Header names: cell/barcode and gene/feature.")
    counts: dict[tuple[str, str], int] = defaultdict(int)
    umis: dict[tuple[str, str], set[str]] = defaultdict(set)
    use_umi = False
    for number, row in enumerate(rows, start=1):
        if len(row) <= max(c_i, g_i):
            raise InputFormatError(f"Record {number} is too short.", code="RECORD_COLUMNS", hint="Each record needs a cell and a gene.")
        key = (str(row[c_i]), str(row[g_i]))
        if u_i is not None and len(row) > u_i and row[u_i] not in (None, ""):
            use_umi = True
            umis[key].add(str(row[u_i]).upper())
        else:
            try:
                counts[key] += int(row[n_i]) if n_i is not None and row[n_i] not in (None, "") else 1
            except (ValueError, TypeError):
                raise InputFormatError(f"Record {number} has a non-integer count.", code="INVALID_COUNT", hint="Counts are integers.") from None
    if use_umi:
        for key, members in umis.items():
            counts[key] += len(members)
    if not counts:
        raise InputFormatError("No records found.", code="EMPTY_INPUT", hint="Provide at least one record.")
    cells = sorted({k[0] for k in counts})
    genes = sorted({k[1] for k in counts})
    if len(cells) * len(genes) > opt_int(options, "max_matrix_elements", DEFAULT_MAX_ELEMENTS, 1, 10**9):  # type: ignore[operator]
        raise InputFormatError(f"The matrix would have {len(cells) * len(genes)} elements.", code="INPUT_TOO_LARGE", hint="Restrict the records to fewer genes/cells.")
    ci, gi = {c: i for i, c in enumerate(cells)}, {g: i for i, g in enumerate(genes)}
    values = [[0.0] * len(cells) for _ in genes]
    for (cell, gene), n in counts.items():
        values[gi[gene]][ci[cell]] = float(n)
    cm = CountMatrix(genes, cells, values)
    definition = "counts = number of distinct UMIs per (cell, gene)" if use_umi else "counts = number of records (or the sum of the count column) per (cell, gene)"
    return matrix_result(cm, options, {"definition": definition, "n_records": len(rows)}), []


def _parse_peak(item: object) -> tuple[str, int, int, str]:
    if isinstance(item, str):
        text = item.strip()
        if ":" in text and "-" in text.split(":", 1)[1]:
            chrom, span = text.rsplit(":", 1)
            start, end = span.split("-", 1)
        else:
            parts = text.replace(",", "\t").split()
            if len(parts) < 3:
                raise InputFormatError(f"Bad peak {item!r}.", code="PEAK_FORMAT", hint="Use 'chr1:100-200' or 'chr1<TAB>100<TAB>200'.")
            chrom, start, end = parts[:3]
    elif isinstance(item, Sequence) and len(item) >= 3:
        chrom, start, end = item[0], item[1], item[2]
    else:
        raise InputFormatError(f"Bad peak {item!r}.", code="PEAK_FORMAT", hint="Use 'chr1:100-200' or [chrom, start, end].")
    try:
        s, e = int(start), int(end)  # type: ignore[call-overload]
    except (TypeError, ValueError):
        raise InputFormatError(f"Bad peak coordinates in {item!r}.", code="PEAK_FORMAT", hint="Coordinates are integers.") from None
    if e <= s or s < 0:
        raise InputFormatError(f"Peak {item!r} must have 0 <= start < end.", code="PEAK_FORMAT", hint="BED-style 0-based half-open intervals.")
    return str(chrom), s, e, f"{chrom}:{s}-{e}"


def op_build_cell_peak_matrix(value: object, options: dict[str, object]) -> Result:
    raw_peaks = options.get("peaks")
    if raw_peaks is None:
        raise ConfigurationError("options['peaks'] (list of peaks) is required.", code="MISSING_OPTION", hint="e.g. ['chr1:100-600', ...] (0-based, half-open).")
    if isinstance(raw_peaks, (str, bytes)) or hasattr(raw_peaks, "__fspath__"):
        raw_peaks = [line for line in text_of(raw_peaks, "peaks").splitlines() if line.strip() and not line.startswith("#")]
    peaks = [_parse_peak(p) for p in raw_peaks]  # type: ignore[union-attr]
    if not peaks:
        raise InputFormatError("No peaks given.", code="EMPTY_INPUT", hint="Provide at least one peak.")
    mode = opt_choice(options, "count", "fragments", ("fragments", "cut_sites"))
    names = [p[3] for p in peaks]
    if len(set(names)) != len(names):
        raise InputFormatError("Duplicate peaks.", code="DUPLICATE_NAMES", hint="List each peak once.")
    by_chrom: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, (chrom, start, end, _name) in enumerate(peaks):
        by_chrom[chrom].append((start, end, index))
    starts: dict[str, list[int]] = {}
    longest: dict[str, int] = {}
    for chrom, items in by_chrom.items():
        items.sort()
        starts[chrom] = [s for s, _e, _i in items]
        longest[chrom] = max(e - s for s, e, _i in items)
    header, rows = _rows(value, "fragments")
    if len(rows) > MAX_RECORDS:
        raise InputFormatError(f"More than {MAX_RECORDS} fragments.", code="INPUT_TOO_LARGE", hint="Process in chunks.")
    counts: dict[tuple[int, str], int] = defaultdict(int)
    cells: set[str] = set()
    for number, row in enumerate(rows, start=1):
        if len(row) < 4:
            raise InputFormatError(f"Fragment {number} needs chrom, start, end, barcode.", code="FRAGMENT_FORMAT", hint="BED-like columns.")
        try:
            chrom, fs, fe = str(row[0]), int(row[1]), int(row[2])  # type: ignore[call-overload]
            weight = int(row[4]) if len(row) > 4 and row[4] not in (None, "") else 1  # type: ignore[call-overload]
        except (TypeError, ValueError):
            raise InputFormatError(f"Fragment {number} has a non-integer coordinate or count.", code="FRAGMENT_FORMAT", hint="Coordinates are integers.") from None
        if fe <= fs:
            raise InputFormatError(f"Fragment {number} has end <= start.", code="FRAGMENT_FORMAT", hint="BED-style 0-based half-open intervals.")
        barcode = str(row[3])
        cells.add(barcode)
        items = by_chrom.get(chrom)
        if not items:
            continue
        if mode == "fragments":
            spans = [(fs, fe)]
        else:
            spans = [(fs, fs + 1), (fe - 1, fe)]
        hits: list[int] = []
        for qs, qe in spans:
            high = bisect.bisect_left(starts[chrom], qe)  # peaks that start before the query ends
            low = bisect.bisect_left(starts[chrom], qs - longest[chrom])  # earlier peaks cannot reach qs
            hits.extend(index for s, e, index in items[low:high] if e > qs)
        for index in hits:
            counts[(index, barcode)] += weight
    cell_names = sorted(cells)
    if not cell_names:
        raise InputFormatError("No fragments found.", code="EMPTY_INPUT", hint="Provide at least one fragment.")
    if len(cell_names) * len(peaks) > opt_int(options, "max_matrix_elements", DEFAULT_MAX_ELEMENTS, 1, 10**9):  # type: ignore[operator]
        raise InputFormatError("The matrix would be too large.", code="INPUT_TOO_LARGE", hint="Use fewer peaks or cells.")
    ci = {c: i for i, c in enumerate(cell_names)}
    values = [[0.0] * len(cell_names) for _ in peaks]
    for (index, barcode), n in counts.items():
        values[index][ci[barcode]] = float(n)
    cm = CountMatrix(names, cell_names, values)
    definition = (
        "counts = fragments (weighted by the optional 5th column) overlapping each peak by >= 1 bp; each fragment counts once per peak"
        if mode == "fragments" else
        "counts = Tn5 cut sites (first and last base of each fragment) falling inside each peak; coordinates are 0-based half-open"
    )
    result = matrix_result(cm, options, {"definition": definition, "n_fragments": len(rows), "mode": mode})
    result = {("n_peaks" if key == "n_genes" else key): value for key, value in result.items()}  # rows are peaks, not genes
    return result, []
