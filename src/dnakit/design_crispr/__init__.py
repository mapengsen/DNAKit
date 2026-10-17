"""Local candidate adapters for the design crispr capability family."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.design_crispr import _impl_local as _impl
from dnakit.exceptions import ConfigurationError, InputFormatError

def assemble_dna_fragments(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assemble DNA fragments in the given order by exact terminal overlaps (in-silico Gibson/overlap joining).

    ``input``: ordered fragments as a list of DNA strings, a mapping ``{id: sequence}`` (insertion order is the
    assembly order) or FASTA text/path. ``options``: ``min_overlap`` (default 15), ``max_overlap`` (default 200):
    consecutive fragments are joined at the longest exact suffix/prefix overlap within these bounds;
    ``circularize`` (default false): also join the last fragment back to the first. Delegates to
    ``dnakit.molbio.simulate_assembly(method="gibson")``. Result: the ``product`` sequence, its length and the
    junction overlaps. No reaction kinetics, secondary structure or mis-annealing is modelled.
    """
    opts = dict(options or {})
    from dnakit.core import DNASequence
    from dnakit.molbio.assembly import AssemblyFragment, simulate_assembly

    items = _impl.parse_sequences(input, "fragments", max_items=1_000)
    if len(items) < 2:
        raise InputFormatError("Assembly needs at least two fragments.", code="CRISPR_INPUT_INVALID")
    ids = [name if name not in {n for n, _ in items[:i]} else f"{name}_{i}" for i, (name, _) in enumerate(items)]
    fragments = [AssemblyFragment(fid, DNASequence(seq)) for fid, (_, seq) in zip(ids, items)]
    result = simulate_assembly(
        fragments, method="gibson", min_overlap=_impl.opt_int(opts, "min_overlap", 15, 1, 100_000),
        max_overlap=_impl.opt_int(opts, "max_overlap", 200, 1, 100_000), circularize=_impl.opt_bool(opts, "circularize", False),
        verbose=False,
    )
    payload = {
        "product": result.product.symbols, "length": len(result.product.symbols), "circular": result.circularized,
        "junctions": [{"left": s.left_fragment_id, "right": s.right_fragment_id, "overlap": s.junction_sequence,
                       "overlap_length": len(s.junction_sequence)} for s in result.steps],
    }
    return candidate_ok("design_crispr", "assemble_dna_fragments", input, payload, options=opts)

def calculate_primer_tm(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Melting temperature (Tm) of one or more primers by the Wallace rule or the nearest-neighbour model.

    ``input``: a primer string, ``{"id": ..., "sequence": ...}``, a mapping ``{name: primer}`` or a list of
    strings / ``{"id", "sequence"}`` mappings (at most 1000 primers of at most 200 bases, 5'->3'). ``options``:
    ``method`` (default ``"auto"``): ``"wallace"`` is the rule of thumb Tm = 2 x (A+T) + 4 x (G+C), meant for
    oligos of about 14 bases or fewer; ``"nearest_neighbor"`` is the SantaLucia (1998) nearest-neighbour Tm of the
    perfectly matched primer/complement duplex from ``dnakit.thermodynamics.melting_temperature`` (2-60 nt);
    ``"auto"`` uses Wallace for length <= 13 and nearest-neighbour above. Nearest-neighbour conditions (mol/L):
    ``sodium_molar`` (default 0.05), ``potassium_molar`` (default 0), ``strand_concentration_molar`` (default
    250e-9); Mg2+, dNTP and cosolvent corrections are not modelled. ``rna`` (default false): read ``U`` as ``T``;
    otherwise ``U`` is rejected. Ambiguous bases (N, IUPAC codes) raise ``InputFormatError``: Tm is undefined for
    them. Result: ``primers``, one entry per primer with ``id``, ``sequence``, ``length``, ``gc_fraction``,
    ``tm_celsius`` (nearest-neighbour values rounded to 0.01) and the ``method`` actually used, plus ``conditions``
    when a nearest-neighbour Tm was computed.
    """
    opts = dict(options or {})
    method = _impl.opt_str(opts, "method", "auto", ("auto", "wallace", "nearest_neighbor"))
    rna = _impl.opt_bool(opts, "rna", False)
    conditions = _impl.tm_conditions(opts)
    rows, warnings = [], []
    for name, seq in _impl.parse_primers(input, rna):
        tm, used = _impl.primer_tm(seq, method, conditions)
        if used == "wallace" and len(seq) > 14:
            warnings.append(f"Primer {name!r}: the Wallace rule overestimates Tm for oligos longer than about 14 bases.")
        rows.append({"id": name, "sequence": seq, "length": len(seq), "gc_fraction": _impl.gc_fraction(seq),
                     "tm_celsius": tm, "method": used})
    result: dict[str, object] = {"primers": rows}
    if any(r["method"] == "nearest_neighbor" for r in rows):
        result["conditions"] = conditions
    return candidate_ok("design_crispr", "calculate_primer_tm", input, result, options=opts, warnings=warnings)

def check_primer_dimer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Screen primers for self- and cross-dimers by the longest perfectly complementary ungapped stretch.

    ``input``: one primer string (self-dimer only), a list of primers, or ``{"forward": ..., "reverse": ...}``.
    Every primer is checked against itself and every other primer. For primers A and B the longest stretch of A
    that is the reverse complement of a stretch of B is ``longest_complementary_run``; ``three_prime_run`` is the
    longest suffix of A (its 3' end) that is reverse-complementary to some part of B. ``options``:
    ``max_3prime_run`` (default 3): a pair is flagged ``three_prime_dimer`` when the run exceeds this value;
    ``max_run`` (default 6): flagged ``long_complementary_run`` when the longest stretch exceeds it.
    This counts contiguous Watson-Crick pairs only; it is not a free-energy (nearest-neighbour) calculation.
    """
    opts = dict(options or {})
    if isinstance(input, Mapping):
        raw = [(str(k), str(v)) for k, v in input.items()]
    elif isinstance(input, str):
        raw = [("primer1", input)]
    elif isinstance(input, Sequence):
        raw = [(f"primer{i + 1}", str(v)) for i, v in enumerate(input)]
    else:
        raise InputFormatError("Provide a primer string, a list or {'forward','reverse'}.", code="CRISPR_INPUT_INVALID")
    if not raw or len(raw) > 100:
        raise InputFormatError("Provide between 1 and 100 primers.", code="CRISPR_INPUT_INVALID")
    primers = [(n, _impl.clean_dna(s, f"primer {n!r}", max_length=200)) for n, s in raw]
    max_3 = _impl.opt_int(opts, "max_3prime_run", 3, 0, 100)
    max_run = _impl.opt_int(opts, "max_run", 6, 0, 100)
    pairs = []
    for i, (name_a, a) in enumerate(primers):
        for name_b, b in primers[i:]:
            rc_b = _impl.revcomp(b)
            length, start_a, _ = _impl.longest_common_substring(a, rc_b)
            run_a, run_b = _impl.three_prime_run(a, b), _impl.three_prime_run(b, a)
            pairs.append({
                "a": name_a, "b": name_b, "self": name_a == name_b, "longest_complementary_run": length,
                "run_start_in_a": start_a if length else None, "three_prime_run_a": run_a, "three_prime_run_b": run_b,
                "three_prime_dimer": max(run_a, run_b) > max_3, "long_complementary_run": length > max_run,
            })
    return candidate_ok("design_crispr", "check_primer_dimer", input, {"pairs": pairs}, options=opts)

def check_primer_secondary_structure(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Screen a primer for hairpins: maximal perfect Watson-Crick stems with a loop of >= ``min_loop`` bases.

    ``input``: one primer string, a list of primers or a mapping ``{name: primer}``. ``options``: ``min_stem``
    (default 4) minimum stem length in base pairs; ``min_loop`` (default 3) minimum unpaired loop length.
    Each primer reports its ``hairpins`` (``stem_length``, ``left_start``, ``right_start``, ``loop_length``, stem
    sequences, whether the stem reaches the 3' end), the ``longest_stem`` and ``three_prime_hairpin``. This is a
    stem-count screen without mismatches, bulges or free energy; it is not a thermodynamic fold (no NUPACK/mfold).
    """
    opts = dict(options or {})
    if isinstance(input, Mapping):
        raw = [(str(k), str(v)) for k, v in input.items()]
    elif isinstance(input, str):
        raw = [("primer1", input)]
    elif isinstance(input, Sequence):
        raw = [(f"primer{i + 1}", str(v)) for i, v in enumerate(input)]
    else:
        raise InputFormatError("Provide a primer string, a list or a mapping.", code="CRISPR_INPUT_INVALID")
    if not raw or len(raw) > 1_000:
        raise InputFormatError("Provide between 1 and 1000 primers.", code="CRISPR_INPUT_INVALID")
    min_stem = _impl.opt_int(opts, "min_stem", 4, 2, 50)
    min_loop = _impl.opt_int(opts, "min_loop", 3, 0, 50)
    out = []
    for name, text in raw:
        seq = _impl.clean_dna(text, f"primer {name!r}", max_length=200)
        hairpins = _impl.find_hairpins(seq, min_stem, min_loop)
        out.append({
            "primer": name, "sequence": seq, "hairpins": hairpins,
            "longest_stem": hairpins[0]["stem_length"] if hairpins else 0,
            "three_prime_hairpin": any(h["involves_3prime_end"] for h in hairpins),
        })
    return candidate_ok("design_crispr", "check_primer_secondary_structure", input, {"primers": out}, options=opts)

def check_primer_specificity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Count where a primer binds (ungapped, mismatches allowed) on each template sequence, both strands.

    ``input``: one primer string, a list of primers or ``{name: primer}``. ``options``: ``templates`` (required:
    mapping ``{id: sequence}``, list of sequences, or FASTA text/path to search); ``max_mismatches`` (default 0);
    ``strict_three_prime_bases`` (default 0: no mismatch allowed in the last N bases). Delegates to
    ``dnakit.molbio.match_primer``. Result per primer: all ``hits`` (template, strand, 0-based start/end, mismatch
    positions 5'->3') and ``total_hits``; ``specific`` is true when exactly one site is found in total.
    """
    opts = dict(options or {})
    from dnakit.core import DNASequence
    from dnakit.molbio.primers import match_primer

    if isinstance(input, Mapping):
        raw = [(str(k), str(v)) for k, v in input.items()]
    elif isinstance(input, str):
        raw = [("primer1", input)]
    elif isinstance(input, Sequence):
        raw = [(f"primer{i + 1}", str(v)) for i, v in enumerate(input)]
    else:
        raise InputFormatError("Provide a primer string, a list or a mapping.", code="CRISPR_INPUT_INVALID")
    if "templates" not in opts:
        raise ConfigurationError(
            "check_primer_specificity needs options['templates'].", code="CRISPR_OPTION_INVALID",
            hint="Pass {'templates': {'chr1': 'ACGT...'}}.",
        )
    templates = _impl.parse_sequences(opts["templates"], "templates", max_items=1_000, allow_iupac=False)
    if sum(len(s) for _, s in templates) > 20_000_000:
        raise ConfigurationError("Templates exceed 20 Mb.", code="CRISPR_SEQUENCE_TOO_LONG")
    max_mm = _impl.opt_int(opts, "max_mismatches", 0, 0, 10)
    strict = _impl.opt_int(opts, "strict_three_prime_bases", 0, 0, 50)
    out = []
    for name, text in raw:
        seq = _impl.clean_dna(text, f"primer {name!r}", max_length=200)
        hits: list[dict[str, object]] = []
        for tid, template in templates:
            result = match_primer(
                DNASequence(seq), DNASequence(template), primer_name=name, max_mismatches=min(max_mm, len(seq)),
                strict_three_prime_bases=min(strict, len(seq)), max_hits=10_000, verbose=False,
            )
            hits.extend({"template": tid, "strand": h.strand.value, "start": h.start, "end": h.end,
                         "mismatch_positions": list(h.mismatch_positions_5to3)} for h in result.hits)
        out.append({"primer": name, "sequence": seq, "total_hits": len(hits), "specific": len(hits) == 1, "hits": hits})
    return candidate_ok("design_crispr", "check_primer_specificity", input, {"primers": out}, options=opts)

def design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Design CRISPR guide RNAs for a target sequence, or an HDR donor template around a cut site.

    ``options["entity"]``:
    ``"grna"`` (default): ``input`` is the target DNA (string, FASTA text/path or ``{"sequence": ...}``). Guides are
    enumerated next to a PAM on both strands (``dnakit.patterns.scan_pam_candidates``), then filtered and ranked.
    Options: ``nuclease`` (``"SpCas9"`` NGG default, ``"SaCas9"``, ``"AsCas12a"``) or custom ``pam``/``pam_side``
    (``"3prime"``|``"5prime"``)/``guide_length``; ``strand`` (both|forward|reverse); ``min_gc`` 0.4, ``max_gc`` 0.7,
    ``max_homopolymer`` 4 (longer runs are rejected), ``avoid_poly_t`` true (rejects TTTT), ``exclude_motifs``
    (IUPAC list), ``avoid_enzymes`` (built-in restriction enzyme names), ``max_guides`` 100; optional
    ``references`` (mapping/FASTA) with ``max_mismatches`` (default 3) to count off-target sites per guide.
    Guides are ordered by ``proxy_score`` (a documented heuristic, see ``score``), then position. With
    ``references`` each guide gains ``exact_sites`` (perfect matches), ``near_sites`` (1..max_mismatches mismatches)
    and ``off_target_sites`` = near_sites + max(exact_sites - 1, 0). Coordinates are
    0-based half-open on the forward strand; ``cut_site`` is the boundary 3 bp from the PAM (SpCas9/SaCas9).
    ``"donor_template"``: ``input`` is the target DNA; options ``cut_site`` (0-based boundary, required),
    ``insert`` (default ""), ``delete_length`` (bases removed after the cut, default 0), ``left_arm``/``right_arm``
    (homology arm lengths, default 60). The donor is ``left arm + insert + right arm``.
    """
    opts = dict(options or {})
    which = _impl.opt_str(opts, "entity", "grna", ("grna", "donor_template"))
    items = _impl.parse_sequences(input, "target", max_items=1, allow_iupac=False)
    name, sequence = items[0]
    if which == "donor_template":
        cut = _impl.opt_int(opts, "cut_site", -1, 0, len(sequence)) if "cut_site" in opts else None
        if cut is None:
            raise ConfigurationError("The donor template needs options['cut_site'].", code="CRISPR_OPTION_INVALID")
        insert = _impl.clean_dna(opts["insert"], "insert") if opts.get("insert") else ""
        deleted = _impl.opt_int(opts, "delete_length", 0, 0, len(sequence))
        left, right = _impl.opt_int(opts, "left_arm", 60, 1, 5_000), _impl.opt_int(opts, "right_arm", 60, 1, 5_000)
        edit_end = cut + deleted
        if cut - left < 0 or edit_end + right > len(sequence):
            raise ConfigurationError(
                "Homology arms extend beyond the target sequence.", code="CRISPR_ARM_OUT_OF_RANGE",
                hint="Shorten left_arm/right_arm or provide a longer target.",
            )
        donor = sequence[cut - left : cut] + insert + sequence[edit_end : edit_end + right]
        result = {
            "target": name, "donor": donor, "donor_length": len(donor),
            "left_arm": {"start": cut - left, "end": cut, "sequence": sequence[cut - left : cut]},
            "right_arm": {"start": edit_end, "end": edit_end + right, "sequence": sequence[edit_end : edit_end + right]},
            "insert": insert, "deleted": sequence[cut:edit_end],
            "edited_sequence": sequence[:cut] + insert + sequence[edit_end:],
        }
        return candidate_ok("design_crispr", "design", input, result, options=opts)
    min_gc, max_gc = _impl.opt_float(opts, "min_gc", 0.4, 0, 1), _impl.opt_float(opts, "max_gc", 0.7, 0, 1)
    max_hp = _impl.opt_int(opts, "max_homopolymer", 4, 1, 100)
    poly_t = _impl.opt_bool(opts, "avoid_poly_t", True)
    motifs = _impl.resolve_motifs(opts)
    limit = _impl.opt_int(opts, "max_guides", 100, 1, 100_000)
    candidates = _impl.enumerate_guides(sequence, opts)
    kept = []
    for guide in candidates:
        checks = _impl.guide_checks(str(guide["sequence"]), min_gc=min_gc, max_gc=max_gc, max_homopolymer=max_hp,
                                    avoid_poly_t=poly_t, motifs=motifs)
        if checks["failed"]:
            continue
        guide.update({"gc": checks["gc"], "longest_homopolymer": checks["longest_homopolymer"],
                      "proxy_score": _impl.proxy_score(str(guide["sequence"]))})
        kept.append(guide)
    kept.sort(key=lambda g: (-g["proxy_score"], g["start"], g["strand"]))  # type: ignore[operator]
    kept = kept[:limit]
    for i, guide in enumerate(kept, 1):
        guide["id"] = f"{name}_g{i}"
    warnings = ["proxy_score is a rule-based heuristic, not a trained on-target efficiency model."]
    if "references" in opts:
        refs = _impl.parse_sequences(opts["references"], "references", max_items=1_000, allow_iupac=True)
        _impl.add_offtargets(kept, refs, opts)
    result = {"target": name, "nuclease": _impl.resolve_nuclease(opts)[0], "candidates_scanned": len(candidates),
              "guides": kept}
    return candidate_ok("design_crispr", "design", input, result, options=opts, warnings=warnings)

def encode_data_into_dna(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Encode bytes as a DNA string with a fixed 2-bits-per-base code (A=00, C=01, G=10, T=11, 4 bases per byte).

    ``input``: raw ``bytes``/``bytearray`` (encoded as they are), a ``pathlib.Path`` (the file's bytes), text
    (UTF-8 bytes) or a mapping/list (sorted-key JSON text, UTF-8). At most 1,000,000 bytes (larger input raises
    ``ConfigurationError``). ``options``: ``source`` (default ``"auto"``): ``"text"`` encodes a string as text,
    ``"path"`` reads a string as a file path, ``"auto"`` treats a ``str`` as text. The code is lossless for any byte
    sequence: ``dnakit.design_crispr._impl_local.dna_to_bytes`` restores the exact bytes. Result: the DNA string.
    This toy mapping has no index, checksum, error correction or synthesis constraints (use
    ``dnakit.genomics.encode_dna_storage`` for indexed segments with an integrity digest).
    """
    opts = dict(options or {})
    source = _impl.opt_str(opts, "source", "auto", ("auto", "text", "path"))
    if isinstance(input, (bytes, bytearray, memoryview)):
        data = bytes(input)
    elif isinstance(input, Path) or (source == "path" and isinstance(input, str)):
        path = Path(input)
        if not path.is_file():
            raise InputFormatError(f"{str(input)[:200]!r} is not an existing file.", code="CRISPR_INPUT_INVALID",
                                   hint="Pass an existing file path, or text with options['source'] = 'text'.")
        if path.stat().st_size > _impl.MAX_STORAGE_BYTES:
            raise ConfigurationError(f"File is larger than {_impl.MAX_STORAGE_BYTES} bytes.", code="CRISPR_INPUT_TOO_LARGE")
        data = path.read_bytes()
    elif isinstance(input, str):
        data = input.encode("utf-8", errors="surrogatepass")
    elif isinstance(input, (Mapping, list, tuple)):
        try:
            data = json.dumps(input, ensure_ascii=False, sort_keys=True).encode("utf-8")
        except (TypeError, ValueError) as exc:
            raise InputFormatError("Mapping/list input must be JSON-serialisable.", code="CRISPR_INPUT_INVALID") from exc
    else:
        raise InputFormatError("Provide bytes, a file path, text or a JSON-serialisable mapping/list.", code="CRISPR_INPUT_INVALID")
    if len(data) > _impl.MAX_STORAGE_BYTES:
        raise ConfigurationError(f"Input is larger than {_impl.MAX_STORAGE_BYTES} bytes.", code="CRISPR_INPUT_TOO_LARGE")
    return candidate_ok("design_crispr", "encode_data_into_dna", input, _impl.bytes_to_dna(data), options=opts)

def evaluate_primer_pair(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Evaluate a forward/reverse primer pair: Tm of each primer, delta Tm, GC, 3' GC clamp and complementarity flags.

    ``input``: ``{"forward": ..., "reverse": ...}`` or a list of two primers, both written 5'->3' as ordered.
    ``options``: ``method``, ``rna``, ``sodium_molar``, ``potassium_molar``, ``strand_concentration_molar`` as in
    ``calculate_primer_tm`` (same Tm for each primer); flag thresholds (defaults in brackets): ``max_delta_tm``
    [5.0 C], ``min_length`` [18], ``max_length`` [30], ``min_gc`` [0.40], ``max_gc`` [0.60], ``min_gc_clamp`` [1]
    and ``max_gc_clamp`` [3] (G/C count among the last 5 bases), ``max_self_complementary`` [4] and
    ``max_cross_complementary`` [4] (bases). Result per primer (``forward``, ``reverse``): ``sequence``, ``length``,
    ``gc_fraction``, ``tm_celsius``, ``method``, ``gc_clamp`` and ``longest_self_complementary`` (longest contiguous
    stretch whose reverse complement also occurs in the same primer, i.e. self-dimer type; hairpins are in
    ``check_primer_secondary_structure``). Pair level: ``delta_tm_celsius`` (absolute), ``longest_cross_complementary``
    (longest contiguous stretch of the forward primer that is the reverse complement of a stretch of the reverse
    primer), ``thresholds``, ``flags`` (each ``{"flag", "primer", "value", "threshold"}``) and ``pass`` (no flags).
    Complementarity is exact contiguous Watson-Crick matching, not a free-energy calculation.
    """
    opts = dict(options or {})
    method = _impl.opt_str(opts, "method", "auto", ("auto", "wallace", "nearest_neighbor"))
    rna = _impl.opt_bool(opts, "rna", False)
    conditions = _impl.tm_conditions(opts)
    limits = {
        "max_delta_tm": _impl.opt_float(opts, "max_delta_tm", 5.0, 0.0, 1000.0),
        "min_length": _impl.opt_int(opts, "min_length", 18, 1, 200), "max_length": _impl.opt_int(opts, "max_length", 30, 1, 200),
        "min_gc": _impl.opt_float(opts, "min_gc", 0.40, 0.0, 1.0), "max_gc": _impl.opt_float(opts, "max_gc", 0.60, 0.0, 1.0),
        "min_gc_clamp": _impl.opt_int(opts, "min_gc_clamp", 1, 0, 5), "max_gc_clamp": _impl.opt_int(opts, "max_gc_clamp", 3, 0, 5),
        "max_self_complementary": _impl.opt_int(opts, "max_self_complementary", 4, 0, 200),
        "max_cross_complementary": _impl.opt_int(opts, "max_cross_complementary", 4, 0, 200),
    }
    if isinstance(input, Mapping) and set(input) == {"forward", "reverse"}:
        raw = [("forward", input["forward"]), ("reverse", input["reverse"])]
    elif isinstance(input, Sequence) and not isinstance(input, (str, bytes, bytearray)) and len(input) == 2:
        raw = [("forward", input[0]), ("reverse", input[1])]
    else:
        raise InputFormatError("Provide {'forward': ..., 'reverse': ...} or a list of two primers.", code="CRISPR_INPUT_INVALID")
    primers: dict[str, dict[str, object]] = {}
    flags: list[dict[str, object]] = []

    def flag(name: str, primer: str | None, value: object, threshold: object) -> None:
        flags.append({"flag": name, "primer": primer, "value": value, "threshold": threshold})

    for role, text in raw:
        seq = _impl.clean_primer(text, f"{role} primer", rna)
        tm, used = _impl.primer_tm(seq, method, conditions)
        entry: dict[str, object] = {
            "sequence": seq, "length": len(seq), "gc_fraction": _impl.gc_fraction(seq), "tm_celsius": tm, "method": used,
            "gc_clamp": _impl.gc_clamp(seq), "longest_self_complementary": _impl.self_complementary_run(seq),
        }
        primers[role] = entry
        if not limits["min_length"] <= len(seq) <= limits["max_length"]:
            flag("length_out_of_range", role, len(seq), [limits["min_length"], limits["max_length"]])
        if not limits["min_gc"] <= entry["gc_fraction"] <= limits["max_gc"]:  # type: ignore[operator]
            flag("gc_out_of_range", role, entry["gc_fraction"], [limits["min_gc"], limits["max_gc"]])
        if not limits["min_gc_clamp"] <= entry["gc_clamp"] <= limits["max_gc_clamp"]:  # type: ignore[operator]
            flag("gc_clamp_out_of_range", role, entry["gc_clamp"], [limits["min_gc_clamp"], limits["max_gc_clamp"]])
        if entry["longest_self_complementary"] > limits["max_self_complementary"]:  # type: ignore[operator]
            flag("self_complementarity", role, entry["longest_self_complementary"], limits["max_self_complementary"])
    fwd, rev = str(primers["forward"]["sequence"]), str(primers["reverse"]["sequence"])
    delta = round(abs(float(primers["forward"]["tm_celsius"]) - float(primers["reverse"]["tm_celsius"])), 2)  # type: ignore[arg-type]
    cross = _impl.longest_common_substring(fwd, _impl.revcomp(rev))[0]
    if delta > limits["max_delta_tm"]:
        flag("tm_mismatch", None, delta, limits["max_delta_tm"])
    if cross > limits["max_cross_complementary"]:
        flag("cross_complementarity", None, cross, limits["max_cross_complementary"])
    warnings = []
    if primers["forward"]["method"] != primers["reverse"]["method"]:
        warnings.append("The two Tm values come from different methods; delta_tm_celsius is not meaningful. Set options['method'].")
    result = {"forward": primers["forward"], "reverse": primers["reverse"], "delta_tm_celsius": delta,
              "longest_cross_complementary": cross, "thresholds": limits, "flags": flags, "pass": not flags}
    if any(p["method"] == "nearest_neighbor" for p in primers.values()):
        result["conditions"] = conditions
    return candidate_ok("design_crispr", "evaluate_primer_pair", input, result, options=opts, warnings=warnings)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise guide RNAs as TSV, CSV, FASTA or JSON text, optionally writing a file.

    ``input``: guides as a list of strings, FASTA text/path, ``{id: guide}`` or a list/``{"guides": [...]}`` of
    records (as returned by ``design``; extra fields such as ``start``, ``strand``, ``proxy_score`` are kept).
    ``options``: ``format`` ``"tsv"`` (default), ``"csv"``, ``"fasta"`` or ``"json"``; ``output`` (optional file
    path, existing files are only replaced when ``overwrite`` is true). Result: ``text`` and ``output_path``.
    """
    opts = dict(options or {})
    fmt = _impl.opt_str(opts, "format", "tsv", ("tsv", "csv", "fasta", "json"))
    guides = _impl.parse_guides(input)
    if fmt == "json":
        text = json.dumps(guides, sort_keys=True)
    elif fmt == "fasta":
        text = "".join(f">{g['id']}\n{g['sequence']}\n" for g in guides)
    else:
        sep = "\t" if fmt == "tsv" else ","
        columns = ["id", "sequence"] + sorted({k for g in guides for k in g} - {"id", "sequence"})
        rows = [sep.join(columns)]
        for g in guides:
            rows.append(sep.join(str(g.get(c, "")) for c in columns))
        text = "\n".join(rows) + "\n"
    written = None
    if opts.get("output") is not None:
        path = Path(str(opts["output"]))
        if path.exists() and not _impl.opt_bool(opts, "overwrite", False):
            raise ConfigurationError(f"Output file {path} already exists.", code="CRISPR_OUTPUT_EXISTS",
                                     hint="Set options['overwrite']=True.")
        path.write_text(text, encoding="utf-8")
        written = str(path)
    return candidate_ok("design_crispr", "export", input, {"format": fmt, "count": len(guides), "text": text,
                                                           "output_path": written}, options=opts)

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter guide RNAs by GC content, homopolymer runs, poly-T, length and excluded motifs, reporting why each was dropped.

    ``input``: guides as strings, FASTA text/path, ``{id: guide}`` or records with a ``guide``/``sequence`` field.
    ``options``: ``min_gc`` (default 0.4), ``max_gc`` (0.7), ``max_homopolymer`` (4: runs longer than this fail),
    ``avoid_poly_t`` (true: any ``TTTT``, the Pol III terminator, fails), ``exclude_motifs`` (IUPAC list, both
    strands), ``avoid_enzymes`` (built-in restriction enzyme names), ``length`` (optional exact guide length).
    Result: ``kept`` and ``rejected`` (each with GC, longest homopolymer and the list of failed rules).
    """
    opts = dict(options or {})
    guides = _impl.parse_guides(input)
    min_gc, max_gc = _impl.opt_float(opts, "min_gc", 0.4, 0, 1), _impl.opt_float(opts, "max_gc", 0.7, 0, 1)
    if min_gc > max_gc:
        raise ConfigurationError("min_gc cannot exceed max_gc.", code="CRISPR_OPTION_INVALID")
    max_hp = _impl.opt_int(opts, "max_homopolymer", 4, 1, 100)
    poly_t = _impl.opt_bool(opts, "avoid_poly_t", True)
    motifs = _impl.resolve_motifs(opts)
    length = _impl.opt_int(opts, "length", 0, 0, 100) if "length" in opts else None
    kept, rejected = [], []
    for guide in guides:
        checks = _impl.guide_checks(str(guide["sequence"]), min_gc=min_gc, max_gc=max_gc, max_homopolymer=max_hp,
                                    avoid_poly_t=poly_t, motifs=motifs)
        failed = list(checks["failed"])  # type: ignore[call-overload]
        if length and len(str(guide["sequence"])) != length:
            failed.append("length")
        row = {**guide, "gc": checks["gc"], "longest_homopolymer": checks["longest_homopolymer"], "failed": failed}
        (rejected if failed else kept).append(row)
    return candidate_ok("design_crispr", "filter", input, {"kept": kept, "rejected": rejected, "kept_count": len(kept),
                                                           "rejected_count": len(rejected)}, options=opts)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw guide RNA positions along a target as an SVG map (forward-strand guides above the line, reverse below).

    ``input``: ``{"length": N, "guides": [...]}`` or ``{"sequence": ..., "guides": [...]}`` where each guide has
    ``start``, ``end`` (0-based half-open) and ``strand`` (``forward``/``reverse``), e.g. the ``result`` of
    ``design``. ``options``: ``output`` (optional ``.svg`` path). Each guide is a ``<rect class="guide">`` with
    ``data-start``, ``data-end`` and ``data-strand`` attributes. Result: ``{"svg": ...}``.
    """
    opts = dict(options or {})
    if not isinstance(input, Mapping) or "guides" not in input or not isinstance(input["guides"], Sequence):
        raise InputFormatError("plot expects {'length' or 'sequence', 'guides': [...]}.", code="CRISPR_INPUT_INVALID")
    guides = list(input["guides"])  # type: ignore[arg-type]
    for g in guides:
        if not isinstance(g, Mapping) or not {"start", "end"} <= set(g) or g["end"] <= g["start"] or g["start"] < 0:  # type: ignore[operator]
            raise InputFormatError("Every guide needs 0 <= start < end.", code="CRISPR_INPUT_INVALID")
    length = input.get("length") or len(str(input.get("sequence", "")))
    if not isinstance(length, int) or length <= 0:
        length = max((int(g["end"]) for g in guides), default=1)
    svg = _impl.svg_guide_map(length, guides)
    written = None
    if opts.get("output") is not None:
        Path(str(opts["output"])).write_text(svg, encoding="utf-8")
        written = str(opts["output"])
    return candidate_ok("design_crispr", "plot", input, {"svg": svg, "output_path": written}, options=opts)

def predict_editing_outcomes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict CRISPR repair outcomes (not implemented: needs a trained outcome model such as inDelphi or FORECasT)."""
    return candidate_not_implemented(
        "design_crispr", "predict_editing_outcomes", input, options=options,
        requires="a trained repair-outcome model (inDelphi, FORECasT, Lindel) and its weights",
        message="Indel/repair outcome prediction needs trained model weights that are not bundled.",
    )

def rank(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Rank guide RNAs by a numeric field; the default ``proxy_score`` is a rule-based heuristic (see ``score``).

    ``input``: guides as strings, FASTA text/path, ``{id: guide}`` or records. ``options``: ``by`` (default
    ``"proxy_score"``, computed when absent; any other value names a numeric field present in every record, for
    instance ``off_target_sites``); ``ascending`` (default false); ``top_n`` (optional). Ties keep input order.
    Result: ``ranked`` records with a 1-based ``rank``.
    """
    opts = dict(options or {})
    guides = _impl.parse_guides(input)
    by = str(opts.get("by", "proxy_score"))
    ascending = _impl.opt_bool(opts, "ascending", False)
    for g in guides:
        if by == "proxy_score" and "proxy_score" not in g:
            g["proxy_score"] = _impl.proxy_score(str(g["sequence"]))
        value = g.get(by)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise InputFormatError(f"Guide {g['id']!r} has no numeric field {by!r}.", code="CRISPR_INPUT_INVALID")
    ranked = sorted(guides, key=lambda g: g[by] if ascending else -g[by])  # type: ignore[operator]
    if "top_n" in opts:
        ranked = ranked[: _impl.opt_int(opts, "top_n", 1, 1, len(ranked))]
    ranked = [{**g, "rank": i} for i, g in enumerate(ranked, 1)]
    warnings = ["proxy_score is a rule-based heuristic, not a trained on-target efficiency model."] if by == "proxy_score" else []
    return candidate_ok("design_crispr", "rank", input, {"by": by, "ascending": ascending, "ranked": ranked},
                        options=opts, warnings=warnings)

def remove_restriction_sites(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Remove restriction-enzyme recognition sites from a DNA sequence by minimal base substitutions.

    ``input``: DNA string, FASTA text/path or ``{"sequence": ...}``. ``options``: ``enzymes`` (built-in names such as
    ``"EcoRI"``) and/or ``sites`` (IUPAC recognition sequences) - at least one is required; sites are removed on both
    strands; ``preserve_translation`` (default false): when true the sequence must be an in-frame coding sequence
    and only synonymous changes are made; ``seed`` (default 0). Delegates to ``dnakit.molbio.optimize_sequence_rules``
    (greedy, deterministic for a given seed). Result: ``sequence`` (edited), ``changes`` (position, before, after),
    ``sites_before``/``sites_after`` counts and ``removed_all``. Sites that cannot be removed under the constraints
    stay and ``removed_all`` is false.
    """
    opts = dict(options or {})
    from dnakit.core import DNASequence
    from dnakit.molbio.optimization import RuleOptimizationConfig, optimize_sequence_rules

    motifs = _impl.resolve_motifs({"exclude_motifs": opts.get("sites", ()), "avoid_enzymes": opts.get("enzymes", ())})
    if not motifs:
        raise ConfigurationError("Provide options['enzymes'] and/or options['sites'].", code="CRISPR_OPTION_INVALID")
    name, sequence = _impl.parse_sequences(input, "sequence", max_items=1)[0]

    def count(seq: str) -> int:
        return sum(len(_impl.motif_positions(seq, m)) for m in motifs.values())

    before = count(sequence)
    config = RuleOptimizationConfig(
        forbidden_motifs=tuple(motifs.values()), forbidden_motifs_both_strands=True,
        preserve_translation=_impl.opt_bool(opts, "preserve_translation", False), seed=_impl.opt_int(opts, "seed", 0, 0, 2**31),
    )
    outcome = optimize_sequence_rules(DNASequence(sequence), config, verbose=False)
    edited = outcome.optimized.symbols
    after = count(edited)
    result = {
        "sequence_id": name, "sequence": edited, "sites_before": before, "sites_after": after, "removed_all": after == 0,
        "changes": [{"position": c.position, "before": c.before, "after": c.after}
                    for c in sorted(outcome.changes, key=lambda c: c.position)],
    }
    return candidate_ok("design_crispr", "remove_restriction_sites", input, result, options=opts)

def run_crispr_screen_design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Design a pooled CRISPR guide library: the best filtered guides per target sequence.

    ``input``: targets as ``{target_id: DNA sequence}``, a list of sequences or FASTA text/path (e.g. exons).
    ``options``: ``guides_per_target`` (default 4); the guide enumeration and filter options of ``design``
    (``nuclease``/``pam``, ``min_gc`` 0.4, ``max_gc`` 0.7, ``max_homopolymer`` 4, ``avoid_poly_t``, ``exclude_motifs``,
    ``avoid_enzymes``); ``min_cut_spacing`` (default 0: selected guides of one target must have cut sites at least
    this many bp apart); optional ``references`` + ``max_mismatches`` (default 3) + ``max_off_target_sites``
    (default: no limit; guides with more non-self sites within the mismatch budget are dropped). Guides are chosen
    greedily by heuristic ``proxy_score``. Result: ``library`` rows, per-target counts and ``targets_short`` (targets
    that got fewer guides than requested). No non-targeting controls are generated.
    """
    opts = dict(options or {})
    targets = _impl.parse_sequences(input, "targets", max_items=10_000)
    per = _impl.opt_int(opts, "guides_per_target", 4, 1, 1_000)
    spacing = _impl.opt_int(opts, "min_cut_spacing", 0, 0, 100_000)
    min_gc, max_gc = _impl.opt_float(opts, "min_gc", 0.4, 0, 1), _impl.opt_float(opts, "max_gc", 0.7, 0, 1)
    max_hp = _impl.opt_int(opts, "max_homopolymer", 4, 1, 100)
    poly_t = _impl.opt_bool(opts, "avoid_poly_t", True)
    motifs = _impl.resolve_motifs(opts)
    refs = _impl.parse_sequences(opts["references"], "references", max_items=1_000, allow_iupac=True) if "references" in opts else None
    max_sites = _impl.opt_int(opts, "max_off_target_sites", -1, -1, 1_000_000)
    library: list[dict[str, object]] = []
    short = []
    for tid, sequence in targets:
        pool = []
        for guide in _impl.enumerate_guides(sequence, opts):
            checks = _impl.guide_checks(str(guide["sequence"]), min_gc=min_gc, max_gc=max_gc, max_homopolymer=max_hp,
                                        avoid_poly_t=poly_t, motifs=motifs)
            if not checks["failed"]:
                guide.update({"gc": checks["gc"], "proxy_score": _impl.proxy_score(str(guide["sequence"]))})
                pool.append(guide)
        pool.sort(key=lambda g: (-g["proxy_score"], g["start"], g["strand"]))  # type: ignore[operator]
        chosen: list[dict[str, object]] = []
        for guide in pool:
            if len(chosen) >= per:
                break
            cut = guide["cut_site"] if guide["cut_site"] is not None else guide["start"]
            if any(abs(int(cut) - int(c["cut_site"] if c["cut_site"] is not None else c["start"])) < spacing for c in chosen):  # type: ignore[call-overload]
                continue
            if refs is not None and max_sites >= 0:
                _impl.add_offtargets([guide], refs, opts)
                if guide["off_target_sites"] > max_sites:  # type: ignore[operator]
                    continue
            chosen.append(guide)
        if refs is not None:
            _impl.add_offtargets([g for g in chosen if "off_target_sites" not in g], refs, opts)
        for i, guide in enumerate(chosen, 1):
            library.append({"target": tid, "guide_id": f"{tid}_g{i}", **guide})
        if len(chosen) < per:
            short.append({"target": tid, "requested": per, "designed": len(chosen)})
    result = {"guides_per_target": per, "target_count": len(targets), "guide_count": len(library), "library": library,
              "targets_short": short}
    return candidate_ok("design_crispr", "run_crispr_screen_design", input, result, options=opts,
                        warnings=["proxy_score is a rule-based heuristic, not a trained on-target efficiency model."])

def scan_grna_offtargets(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find off-target sites of guide RNAs in supplied sequences by exact ungapped Hamming search on both strands.

    ``input``: guides as strings, FASTA text/path, ``{id: guide}`` or records. ``options``: ``references``
    (required: mapping ``{id: sequence}``, list of sequences, or FASTA text/path to search); ``max_mismatches``
    (default 3, 0..6); ``pam`` (IUPAC, default ``"NGG"``; ``null`` ignores the PAM) and ``pam_side``
    (``"3prime"`` default | ``"5prime"``). A site is reported when the protospacer differs from the guide at
    <= ``max_mismatches`` positions and (if ``pam``) the adjacent PAM matches. Hits give reference id, strand,
    0-based half-open forward coordinates, mismatch count and 1-based mismatch positions from the guide 5' end.
    The guide's own target is reported too (0 mismatches). No bulges/indels, no CFD/MIT scoring.
    """
    opts = dict(options or {})
    guides = _impl.parse_guides(input, max_items=1_000)
    if "references" not in opts:
        raise ConfigurationError("scan_grna_offtargets needs options['references'].", code="CRISPR_OPTION_INVALID",
                                 hint="Pass {'references': {'chr1': 'ACGT...'}}.")
    refs = _impl.parse_sequences(opts["references"], "references", max_items=10_000, allow_iupac=True)
    max_mm = _impl.opt_int(opts, "max_mismatches", 3, 0, 6)
    pam = "NGG" if "pam" not in opts else opts["pam"]
    pam = None if pam is None else _impl.clean_dna(str(pam), "pam", allow_iupac=True, max_length=20)
    side = _impl.opt_str(opts, "pam_side", "3prime", ("3prime", "5prime"))
    out = _impl.scan_offtargets(guides, refs, max_mm, pam, side)
    return candidate_ok("design_crispr", "scan_grna_offtargets", input, {"max_mismatches": max_mm, "pam": pam,
                                                                          "guides": out}, options=opts)

def score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score guide RNAs with a rule-based HEURISTIC on-target proxy in [0, 100] (not a trained efficiency model).

    ``input``: guides as strings, FASTA text/path, ``{id: guide}`` or records. The proxy starts at 100 and subtracts:
    15 when GC lies outside 40-70 % (40 instead when outside 20-80 %); 40 for a ``TTTT`` run (Pol III
    terminator); 15 for an A/C/G homopolymer of length >= 4; floor 0. It encodes only common design rules of
    thumb; it was not fitted or validated against editing data and must not be read as predicted efficiency
    (use Rule Set 2, DeepSpCas9 or similar for that). Result per guide: ``proxy_score``, ``gc``,
    ``longest_homopolymer``, ``poly_t``.
    """
    opts = dict(options or {})
    guides = _impl.parse_guides(input)
    rows = []
    for g in guides:
        seq = str(g["sequence"])
        if not 10 <= len(seq) <= 40:
            raise InputFormatError(f"Guide {g['id']!r} has length {len(seq)}; expected 10-40.", code="CRISPR_SEQUENCE_INVALID")
        rows.append({"id": g["id"], "sequence": seq, "proxy_score": _impl.proxy_score(seq), "gc": _impl.gc_fraction(seq),
                     "longest_homopolymer": _impl.longest_run(seq), "poly_t": "TTTT" in seq})
    return candidate_ok("design_crispr", "score", input, {"method": "heuristic_rule_proxy", "scores": rows}, options=opts,
                        warnings=["Heuristic rule-based proxy, not a trained on-target efficiency model."])

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate homology-directed repair in silico: swap the region between a donor's homology arms into the target.

    ``input``: target DNA (string, FASTA text/path or ``{"sequence": ...}``). ``options``: ``donor`` (required,
    donor DNA = ``left arm + payload + right arm``), ``left_arm`` and ``right_arm`` (required arm lengths).
    Both arms must occur exactly once on the forward strand of the target, left before right; the target
    sequence between them is replaced by the payload. Result: ``edited_sequence``, the replaced ``region`` (0-based
    half-open), ``deleted`` and ``inserted`` sequences and the net length change. No repair-pathway competition,
    cutting efficiency or NHEJ indels are modelled.
    """
    opts = dict(options or {})
    name, target = _impl.parse_sequences(input, "target", max_items=1)[0]
    if "donor" not in opts:
        raise ConfigurationError("simulate needs options['donor'].", code="CRISPR_OPTION_INVALID")
    donor = _impl.clean_dna(opts["donor"], "donor")
    left_len = _impl.opt_int(opts, "left_arm", 0, 1, len(donor)) if "left_arm" in opts else 0
    right_len = _impl.opt_int(opts, "right_arm", 0, 1, len(donor)) if "right_arm" in opts else 0
    if not left_len or not right_len or left_len + right_len > len(donor):
        raise ConfigurationError("simulate needs left_arm and right_arm lengths that fit in the donor.",
                                 code="CRISPR_OPTION_INVALID")
    left, right, payload = donor[:left_len], donor[len(donor) - right_len :], donor[left_len : len(donor) - right_len]
    if target.count(left) != 1 or target.count(right) != 1:
        raise InputFormatError(
            "Each homology arm must occur exactly once in the target.", code="CRISPR_ARM_NOT_UNIQUE",
            context={"left_arm_matches": target.count(left), "right_arm_matches": target.count(right)},
        )
    left_end = target.index(left) + left_len
    right_start = target.index(right)
    if right_start < left_end:
        raise InputFormatError("The right arm lies before the end of the left arm.", code="CRISPR_ARM_ORDER")
    edited = target[:left_end] + payload + target[right_start:]
    result = {
        "target": name, "edited_sequence": edited, "region": {"start": left_end, "end": right_start},
        "deleted": target[left_end:right_start], "inserted": payload, "length_change": len(edited) - len(target),
    }
    return candidate_ok("design_crispr", "simulate", input, result, options=opts)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate guide RNAs: alphabet, length, PAM adjacency in a target (optional) and design-rule warnings.

    ``input``: guides as strings, FASTA text/path, ``{id: guide}`` or records. ``options``: ``length`` (expected
    guide length; default the nuclease default, 20 for SpCas9); ``nuclease``/``pam``/``pam_side`` as in ``design``;
    ``target`` (optional DNA): each guide must then occur in the target on either strand next to a matching PAM
    (``target_sites`` lists the sites). ``errors`` make a guide invalid (bad characters, wrong length, no PAM site
    in the target); ``warnings`` (GC outside 40-70 %, homopolymer run > 4, ``TTTT``) do not.
    """
    opts = dict(options or {})
    records = [{"id": g["id"], "raw": g["sequence"]} for g in _impl.parse_guides(input, strict=False)]
    name, pam, side, default_len = _impl.resolve_nuclease(opts)
    length = _impl.opt_int(opts, "length", default_len, 10, 40)
    target = _impl.parse_sequences(opts["target"], "target", max_items=1)[0][1] if "target" in opts else None
    pam_re = re.compile("".join(_impl._IUPAC[c] for c in pam))
    rows = []
    for rec in records:
        text = re.sub(r"\s+", "", str(rec["raw"])).upper().replace("U", "T")
        errors, warnings, sites = [], [], []
        if not text or set(text) - set("ACGT"):
            errors.append("invalid_characters")
        else:
            if len(text) != length:
                errors.append(f"length_{len(text)}_expected_{length}")
            gc = _impl.gc_fraction(text)
            if not 0.4 <= gc <= 0.7:
                warnings.append("gc_outside_40_70")
            if _impl.longest_run(text) > 4:
                warnings.append("homopolymer_run_gt_4")
            if "TTTT" in text:
                warnings.append("poly_t")
            if target is not None and not errors:
                for strand, seq in (("forward", target), ("reverse", _impl.revcomp(target))):
                    start = seq.find(text)
                    while start != -1:
                        pam_text = seq[start + length : start + length + len(pam)] if side == "3prime" else seq[max(0, start - len(pam)) : start]
                        if len(pam_text) == len(pam) and pam_re.fullmatch(pam_text):
                            pos = start if strand == "forward" else len(target) - start - length
                            sites.append({"strand": strand, "start": pos, "end": pos + length, "pam": pam_text})
                        start = seq.find(text, start + 1)
                if not sites:
                    errors.append("no_pam_adjacent_site_in_target")
        rows.append({"id": rec["id"], "sequence": text, "valid": not errors, "errors": errors, "warnings": warnings,
                     "target_sites": sites})
    return candidate_ok("design_crispr", "validate", input,
                        {"nuclease": name, "pam": pam, "expected_length": length, "guides": rows,
                         "all_valid": all(r["valid"] for r in rows)}, options=opts)

__all__ = ['assemble_dna_fragments', 'calculate_primer_tm', 'check_primer_dimer', 'check_primer_secondary_structure', 'check_primer_specificity', 'design', 'encode_data_into_dna', 'evaluate_primer_pair', 'export', 'filter', 'plot', 'predict_editing_outcomes', 'rank', 'remove_restriction_sites', 'run_crispr_screen_design', 'scan_grna_offtargets', 'score', 'simulate', 'validate']

__candidate_adapter_operations__ = set(['assemble_dna_fragments', 'calculate_primer_tm', 'check_primer_dimer', 'check_primer_secondary_structure', 'check_primer_specificity', 'design', 'encode_data_into_dna', 'evaluate_primer_pair', 'export', 'filter', 'plot', 'predict_editing_outcomes', 'rank', 'remove_restriction_sites', 'run_crispr_screen_design', 'scan_grna_offtargets', 'score', 'simulate', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
