"""Local candidate adapters for the regulation capability family."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.regulation import _impl_local as _impl


def _link_regions(operation: str, input: object, opts: Mapping[str, object]) -> dict[str, object]:
    """Distance-based region -> gene assignment shared by ``link`` and ``annotate(entity="target_gene")``."""

    from dnakit.epigenome import _impl_local as epi

    if "genes" not in opts:
        raise ConfigurationError(f"{operation} needs options['genes'].", code="REGULATION_OPTION_INVALID",
                                 hint="Pass genes as rows/mappings with chrom, start, end, name and strand.")
    method = _impl.opt_str(opts, "method", "nearest", ("nearest", "window"))
    max_distance = _impl.opt_int(opts, "max_distance", 1_000_000, 0, 10**12)
    links = epi.link_regions_to_genes(epi.parse_intervals(input), epi.parse_genes(opts["genes"]), method, max_distance)
    return candidate_ok("regulation", operation, input, {"method": method, "max_distance": max_distance, "links": links},
                        options=opts, warnings=["Distance-based assignment is a proximity heuristic, not evidence of regulation."])

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate sequences with PWM motif hits, or regions with their nearest target gene.

    ``options["entity"]``:
    ``"motif"`` (default; alias ``"tf_binding_site"``): ``input`` is DNA (string, FASTA text/path, list or
    ``{id: sequence}``); ``options["motifs"]`` (required) is one PWM or several (list / ``{name: pwm}``), each given
    as counts or probabilities in the forms ``{"A": [...], "C": [...], "G": [...], "T": [...]}``, rows
    ``[[a, c, g, t], ...]`` or JASPAR text (see ``validate``). Score = sum over positions of log2(p / background)
    (``background`` default 0.25 each; ``pseudocount`` default 0 for counts; zero probabilities are floored at
    1e-6). Both strands are scanned (``strand``: both|forward|reverse); a window is a hit when its score is >=
    ``threshold`` (absolute bits) or, if no threshold is given, its relative score (score - min) / (max - min) is >=
    ``min_relative_score`` (default 0.8). Hits report 0-based half-open start/end on the forward strand, strand,
    score, relative score and the motif-oriented matched sequence. Windows containing N are skipped.
    ``"target_gene"``: ``input`` is regions; same as ``link`` with method ``nearest``.
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "motif", ("motif", "tf_binding_site", "target_gene"))
    if which == "target_gene":
        return _link_regions("annotate", input, {**opts, "method": "nearest"})
    pwms = _impl.parse_motifs(opts)
    threshold, min_rel = _impl.resolve_scan_thresholds(opts)
    strand = _impl.strand_option(opts)
    cap = _impl.opt_int(opts, "max_hits", 10_000, 1, 1_000_000)
    out = []
    for sid, seq in _impl.parse_sequences(input):
        hits = []
        for pwm in pwms:
            if len(seq) >= pwm.length:
                hits.extend({"motif": pwm.name, **h} for h in _impl.scan(seq, pwm, threshold, min_rel, strand))
        hits.sort(key=lambda h: (h["start"], h["motif"], h["strand"]))  # type: ignore[arg-type]
        out.append({"id": sid, "length": len(seq), "hit_count": len(hits), "hits": hits[:cap], "truncated": len(hits) > cap})
    return candidate_ok("regulation", "annotate", input, {"entity": which, "sequences": out}, options=opts)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate PWM-based regulatory effect deltas of variants, or the information content of motifs.

    ``options["entity"]``:
    ``"regulatory_effect_delta"`` (default): ``input`` is a variant record or list, either
    ``{"id", "ref_sequence", "alt_sequence"}`` or ``{"id", "sequence", "position", "ref", "alt"}`` (0-based position;
    ``ref`` must match). ``options["motifs"]`` (required) as in ``annotate``. For every motif the best log2-odds
    score over all windows and both strands is computed for the reference and the alternate sequence;
    ``delta = alt_score - ref_score`` (bits; negative means weaker binding). This is a position-weight-matrix
    proxy, not a trained model.
    ``"pwm_information_content"``: ``input`` is one PWM or several (the ``motifs`` forms of ``annotate``); reports
    per-position information content sum(p log2(p / bg)) and its total in bits (2 bits per fully conserved
    column with a uniform background).
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "regulatory_effect_delta", ("regulatory_effect_delta", "pwm_information_content"))
    if which == "pwm_information_content":
        pwms = _impl.parse_motifs({**opts, "motifs": input})
        result: dict[str, object] = {"entity": which, "motifs": [
            {"motif": p.name, "length": p.length, **_impl.information_content(p)} for p in pwms]}
    else:
        pwms = _impl.parse_motifs(opts)
        result = {"entity": which, "variants": _impl.variant_effects(_impl.parse_variants(input), pwms)}
    return candidate_ok("regulation", "calculate", input, result, options=opts,
                        warnings=["PWM log-odds is a simple binding-affinity proxy."] if which == "regulatory_effect_delta" else [])

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Classify regions as promoter, proximal or distal by position relative to gene TSSs (position rule only).

    ``input``: regions as BED text/path, rows or mappings (0-based half-open). ``options``: ``genes`` (required:
    ``chrom``, ``start``, ``end``, ``name``, optional ``strand``/``tss``; TSS = ``start`` for + or unstranded,
    ``end - 1`` for - genes); ``promoter_upstream`` (default 1000) and ``promoter_downstream`` (default 100): the
    promoter window is ``[tss - up, tss + down)`` for + genes and ``[tss - down + 1, tss + up + 1)`` for - genes.
    A region overlapping any promoter window is ``promoter``; otherwise ``proximal`` when its midpoint
    ``(start + end) // 2`` is within ``proximal_distance`` (default 10000) of the nearest TSS on the same
    chromosome, else ``distal`` (also when the chromosome has no genes). No chromatin or conservation evidence is used.
    """
    opts = dict(options or {})
    _impl.entity(opts, "promoter", ("promoter",))
    from dnakit.epigenome import _impl_local as epi

    if "genes" not in opts:
        raise ConfigurationError("classify needs options['genes'].", code="REGULATION_OPTION_INVALID")
    up = _impl.opt_int(opts, "promoter_upstream", 1000, 0, 10**7)
    down = _impl.opt_int(opts, "promoter_downstream", 100, 1, 10**7)
    proximal = _impl.opt_int(opts, "proximal_distance", 10_000, 0, 10**9)
    genes = epi.parse_genes(opts["genes"])
    regions = epi.parse_intervals(input)
    rows = []
    for r in regions:
        rs, re_ = int(r["start"]), int(r["end"])  # type: ignore[call-overload]
        mid = (rs + re_) // 2
        label, nearest = "distal", None
        for g in genes:
            if g["chrom"] != r["chrom"]:
                continue
            tss = epi.gene_tss(g)
            lo, hi = (tss - up, tss + down) if g.get("strand") != "-" else (tss - down + 1, tss + up + 1)
            dist = abs(tss - mid)
            if nearest is None or (dist, str(g["name"])) < (nearest[0], nearest[1]):
                nearest = (dist, str(g["name"]))
            if rs < hi and re_ > lo:
                label = "promoter"
            elif label != "promoter" and dist <= proximal:
                label = "proximal"
        rows.append({"chrom": r["chrom"], "start": rs, "end": re_, "class": label,
                     "nearest_gene": nearest[1] if nearest else None, "distance_to_tss": nearest[0] if nearest else None})
    return candidate_ok("regulation", "classify", input, {"regions": rows}, options=opts,
                        warnings=["Classification uses genomic position only."])

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare PWM motifs pairwise by best ungapped alignment (mean per-column Pearson correlation).

    ``input``: two or more motifs (list of matrices, ``{name: matrix}``; matrices as in ``annotate``, counts or
    probabilities). ``options``: ``entity`` (only ``"motif"``); ``min_overlap`` (default 3, reduced to the shortest
    motif length): minimum number of aligned columns. For motifs a and b every offset and both orientations of b
    (forward and reverse complement) are tried; the similarity of an alignment is the mean over aligned columns of the
    Pearson correlation between the A/C/G/T probability vectors (uniform columns count 0). Result: ``pairs`` with the
    best ``similarity`` (1.0 = identical), ``offset`` (position of b's first column relative to a's), ``orientation``
    and ``overlap``. This is a simple column-correlation score, not Tomtom's p-value/E-value.
    """
    opts = dict(options or {})
    _impl.entity(opts, "motif", ("motif",))
    pwms = _impl.parse_motifs({**opts, "motifs": input})
    if len(pwms) < 2:
        raise InputFormatError("compare needs at least two motifs.", code="REGULATION_INPUT_INVALID")
    wanted = _impl.opt_int(opts, "min_overlap", 3, 1, 1_000)
    pairs = []
    for i in range(len(pwms)):
        for j in range(i + 1, len(pwms)):
            overlap = min(wanted, pwms[i].length, pwms[j].length)
            pairs.append({"a": pwms[i].name, "b": pwms[j].name, **_impl.compare_pwms(pwms[i], pwms[j], overlap)})
    return candidate_ok("regulation", "compare", input, {"pairs": pairs}, options=opts)

def enrich(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test motif enrichment in foreground vs background sequences (hypergeometric or binomial upper tail).

    ``options["entity"]``: ``"motif"`` (default, alias ``"tf_binding_site"``). ``input`` is either counts
    ``{"foreground_hits", "foreground_total", "background_hits", "background_total"}`` (number of sequences with at least
    one motif hit) or sequences ``{"foreground": [...], "background": [...]}``; with sequences ``options["motifs"]``
    plus ``threshold``/``min_relative_score``/``strand`` (see ``annotate``) define a hit, and each motif is tested.
    ``options["test"]``: ``"hypergeometric"`` (default): with N = fg_total + bg_total, K = fg_hits + bg_hits,
    n = fg_total, p = P(X >= fg_hits) for X ~ Hypergeometric(N, K, n); ``"binomial"``: p = P(X >= fg_hits) for
    X ~ Binomial(fg_total, bg_hits / bg_total). Also reported: ``fold_enrichment = (fg_hits / fg_total) /
    (bg_hits / bg_total)`` (null when the background rate is 0) and Benjamini-Hochberg ``q_value`` over motifs.
    """
    opts = dict(options or {})
    _impl.entity(opts, "motif", ("motif", "tf_binding_site"))
    test = _impl.opt_str(opts, "test", "hypergeometric", ("hypergeometric", "binomial"))
    if not isinstance(input, Mapping):
        raise InputFormatError("enrich expects a mapping of counts or of foreground/background sequences.",
                               code="REGULATION_INPUT_INVALID")
    entries: list[tuple[str, int, int, int, int]] = []
    count_keys = ("foreground_hits", "foreground_total", "background_hits", "background_total")
    if all(k in input for k in count_keys):
        vals = [input[k] for k in count_keys]
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in vals):
            raise InputFormatError("Counts must be non-negative integers.", code="REGULATION_INPUT_INVALID")
        entries.append((str(opts.get("name", "motif")), *vals))  # type: ignore[arg-type]
    elif "foreground" in input and "background" in input:
        pwms = _impl.parse_motifs(opts)
        threshold, min_rel = _impl.resolve_scan_thresholds(opts)
        strand = _impl.strand_option(opts)
        fg, bg = _impl.parse_sequences(input["foreground"], "foreground"), _impl.parse_sequences(input["background"], "background")
        for pwm in pwms:
            def has_hit(seq: str, pwm: object = pwm) -> bool:
                return len(seq) >= pwm.length and bool(_impl.scan(seq, pwm, threshold, min_rel, strand))  # type: ignore[attr-defined]
            entries.append((pwm.name, sum(has_hit(s) for _, s in fg), len(fg), sum(has_hit(s) for _, s in bg), len(bg)))
    else:
        raise InputFormatError("Provide the four count keys or foreground/background sequences.", code="REGULATION_INPUT_INVALID")
    rows = []
    for name, k, n, bk, bn in entries:
        if k > n or bk > bn or n == 0 or bn == 0:
            raise InputFormatError(f"Invalid counts for {name!r}: hits cannot exceed totals and totals must be positive.",
                                   code="REGULATION_INPUT_INVALID")
        if test == "hypergeometric":
            p = _impl.hypergeom_sf(k, n + bn, k + bk, n)
        else:
            p = _impl.binom_sf(k, n, bk / bn)
        rows.append({"motif": name, "foreground_hits": k, "foreground_total": n, "background_hits": bk, "background_total": bn,
                     "fold_enrichment": (k / n) / (bk / bn) if bk else None, "p_value": p})
    for row, q in zip(rows, _impl.benjamini_hochberg([r["p_value"] for r in rows])):  # type: ignore[misc]
        row["q_value"] = q
    return candidate_ok("regulation", "enrich", input, {"test": test, "motifs": rows}, options=opts)

def link(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign regulatory regions (enhancers etc.) to genes by distance from the region midpoint to the gene TSS.

    ``input``: regions as BED text/path, rows or mappings (0-based half-open). ``options``: ``entity`` (one of
    ``enhancer`` (default), ``promoter``, ``silencer``, ``insulator``, ``target_gene``: only a label, the rule is the
    same); ``genes`` (required: ``chrom``, ``start``, ``end``, ``name``, optional ``strand``/``tss``); ``method``:
    ``"nearest"`` (default) or ``"window"``; ``max_distance`` (default 1,000,000 bp). TSS = ``start`` for + or
    unstranded genes, ``end - 1`` for - genes; distance from the region midpoint ``(start + end) // 2``; ties go to the
    alphabetically first gene. A proximity heuristic: no contact or expression evidence is used.
    """
    opts = dict(options or {})
    _impl.entity(opts, "enhancer", ("enhancer", "promoter", "silencer", "insulator", "target_gene"))
    return _link_regions("link", input, opts)

def perform_in_silico_mutagenesis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """In-silico saturation mutagenesis of a sequence scored by a PWM (every single-base substitution).

    ``input``: one DNA sequence (string, FASTA text/path or ``{"sequence": ...}``), at most 2000 bases, ACGT only.
    ``options["motifs"]`` (required, as in ``annotate``). For every position i and alternative base b the effect is
    ``best(mutated) - best(reference)`` where ``best`` is the maximum PWM log2-odds score over all windows and both
    strands. Result per motif: ``motif``, ``reference_score`` and ``positions``, a list with one entry per sequence
    position ``{"position", "ref", "effects"}`` where ``effects`` maps each of the three alternative bases to its
    ``delta_bits`` (null when no scorable window remains); the reference base itself is not listed. The model is
    the PWM itself, not a trained neural network.
    """
    opts = dict(options or {})
    pwms = _impl.parse_motifs(opts)
    items = _impl.parse_sequences(input, max_items=1)
    name, seq = items[0]
    if len(seq) > 2000 or "N" in seq:
        raise InputFormatError("In-silico mutagenesis needs an ACGT sequence of at most 2000 bases.", code="REGULATION_SEQUENCE_INVALID")
    out = []
    for pwm in pwms:
        if len(seq) < pwm.length:
            raise InputFormatError(f"Sequence is shorter than motif {pwm.name!r} ({pwm.length}).", code="REGULATION_SEQUENCE_INVALID")
        ref = _impl.best_score(seq, pwm)
        out.append({"motif": pwm.name, "reference_score": ref, "positions": _impl.ism_matrix(seq, pwm)})
    return candidate_ok("regulation", "perform_in_silico_mutagenesis", input, {"sequence": name, "motifs": out}, options=opts,
                        warnings=["The scoring model is a PWM, not a trained sequence model."])

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict regulatory activity from sequence (not implemented: needs a trained sequence-to-function model)."""
    return candidate_not_implemented(
        "regulation", "predict", input, options=options,
        requires="a trained regulatory sequence model and weights (Enformer, Borzoi, DeepSEA, Sei)",
        message="Regulatory-activity prediction needs trained model weights that are not bundled.",
    )

def prioritize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Rank variants by the largest absolute PWM score change across motifs (a proxy for regulatory impact).

    ``input`` and ``options`` follow ``calculate`` (``regulatory_effect_delta``): variants as ``ref_sequence`` +
    ``alt_sequence`` or ``sequence`` + ``position`` + ``ref`` + ``alt``, and ``options["motifs"]``. ``options``:
    ``entity`` (only ``"regulatory_variant"``); ``top_n`` (optional). Variants are sorted by ``max_abs_delta``
    (bits) descending; variants without any scorable motif window go last. This ranks by a PWM proxy and does not
    use conservation, chromatin or a trained model.
    """
    opts = dict(options or {})
    _impl.entity(opts, "regulatory_variant", ("regulatory_variant",))
    pwms = _impl.parse_motifs(opts)
    rows = _impl.variant_effects(_impl.parse_variants(input), pwms)
    ranked = sorted(rows, key=lambda r: (r["max_abs_delta"] is None, -(r["max_abs_delta"] or 0.0), str(r["id"])))  # type: ignore[operator]
    if "top_n" in opts:
        ranked = ranked[: _impl.opt_int(opts, "top_n", 1, 1, max(len(ranked), 1))]
    ranked = [{**r, "rank": i} for i, r in enumerate(ranked, 1)]
    return candidate_ok("regulation", "prioritize", input, {"ranked": ranked}, options=opts,
                        warnings=["Ranking uses a PWM score change proxy only."])

def score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score sequences against PWM motifs: best log2-odds window per sequence and motif.

    ``input``: DNA (string, FASTA text/path, list or ``{id: sequence}``). ``options``: ``motifs`` (required, as in
    ``annotate``: counts or probabilities, ``background``, ``pseudocount``, ``matrix_type``); ``entity``
    (``"motif"`` default or ``"tf_binding_site"``); ``strand`` (both|forward|reverse). Score = sum over positions of
    log2(p(base) / background(base)); relative score = (score - min) / (max - min) of the matrix. Result per
    sequence and motif: the best ``start``/``end`` (0-based half-open, forward coordinates), ``strand``, ``score``,
    ``relative_score``; null when the sequence is shorter than the motif or has only N windows.
    """
    opts = dict(options or {})
    _impl.entity(opts, "motif", ("motif", "tf_binding_site"))
    pwms = _impl.parse_motifs(opts)
    strand = _impl.strand_option(opts)
    rows = []
    for sid, seq in _impl.parse_sequences(input):
        rows.append({"id": sid, "motifs": [
            {"motif": p.name, "best": _impl.best_hit(seq, p, strand) if len(seq) >= p.length else None} for p in pwms]})
    return candidate_ok("regulation", "score", input, {"sequences": rows}, options=opts)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate PWM motifs: shape, finite non-negative entries, non-zero columns and whether columns sum to 1.

    ``input``: one PWM or several (list / ``{name: matrix}``) as counts or probabilities in the forms accepted by
    ``annotate``. ``options``: ``entity`` (only ``"motif"``); ``tolerance`` (default 1e-3) for the column-sum test.
    Result per motif: ``valid`` (false when the matrix is malformed), ``errors``, ``matrix_type`` (``probabilities``
    when every column sums to 1 within tolerance, else ``counts``), ``length``, and for valid matrices the total
    information content in bits.
    """
    opts = dict(options or {})
    _impl.entity(opts, "motif", ("motif",))
    tol = _impl.opt_float(opts, "tolerance", 1e-3, 0.0, 1.0)
    raw = input
    if isinstance(raw, Mapping) and not (all(b in raw for b in "ACGT") or any(k in raw for k in ("matrix", "pwm", "ppm", "pfm", "counts", "probabilities"))):
        named = [(str(k), v) for k, v in raw.items()]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)) and raw and not all(
        isinstance(r, Sequence) and not isinstance(r, str) and len(r) == 4 and all(isinstance(x, (int, float)) for x in r) for r in raw
    ):
        named = [(f"motif{i + 1}", v) for i, v in enumerate(raw)]
    else:
        named = [("motif1", raw)]
    rows = []
    for name, value in named:
        entry: dict[str, object] = {"motif": name, "valid": True, "errors": []}
        try:
            matrix = _impl._matrix_rows(value, name)
            entry["length"] = len(matrix)
            sums = [sum(r) for r in matrix]
            bad = [i for i, r in enumerate(matrix) if any(not math.isfinite(x) or x < 0 for x in r)]
            if bad:
                raise InputFormatError(f"{name}: negative or non-finite entries at positions {bad}.", code="REGULATION_PWM_INVALID")
            if any(s <= 0 for s in sums):
                raise InputFormatError(f"{name}: a column sums to zero.", code="REGULATION_PWM_INVALID")
            entry["matrix_type"] = "probabilities" if all(abs(s - 1.0) <= tol for s in sums) else "counts"
            entry["total_information_bits"] = _impl.information_content(_impl.build_pwm(value, name, opts))["total_bits"]
        except InputFormatError as exc:
            entry["valid"] = False
            entry["errors"] = [exc.message]
        rows.append(entry)
    return candidate_ok("regulation", "validate", input, {"motifs": rows, "all_valid": all(r["valid"] for r in rows)}, options=opts)

__all__ = ['annotate', 'calculate', 'classify', 'compare', 'enrich', 'link', 'perform_in_silico_mutagenesis', 'predict', 'prioritize', 'score', 'validate']

__candidate_adapter_operations__ = set(['annotate', 'calculate', 'classify', 'compare', 'enrich', 'link', 'perform_in_silico_mutagenesis', 'predict', 'prioritize', 'score', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    count_peak_fragments,
    count_tile_matrix,
    deep_lift_shap_attribution,
    marginalize_motif_effect,
    predict_sequences,
    read_fragment_qc,
    read_motif_file,
    scan_motif_hits,
    score_motif_kmers,
)

__all__ = [*__all__, "count_peak_fragments", "count_tile_matrix", "deep_lift_shap_attribution", "marginalize_motif_effect", "predict_sequences", "read_fragment_qc", "read_motif_file", "scan_motif_hits", "score_motif_kmers"]

# Full-pass exports (see docs/api/features)
from .gimmemotifs_functions import (
    cluster_motif_set,
    compare_motif_positions,
    equal_size_intervals,
    motif_enrichment_stats,
    motif_from_alignment,
    motif_from_consensus,
    motif_score_correlation,
    promoter_regions,
    random_genomic_regions,
    read_peak_summits,
    scan_best_motif_hits,
    summit_intervals_from_narrowpeak,
)

from .snapatac2_functions import (
    chromosome_fragment_counts,
    fragment_in_peaks_fraction,
    fragment_size_distribution,
    gene_activity_counts,
    merge_peak_sets,
    region_motif_enrichment,
    regulatory_regions_to_genes,
    tss_enrichment,
)

from .tangermeme_functions import (
    ablate_region,
    call_seqlets,
    deletion_effect,
    design_greedy_substitution,
    extract_loci_sequences,
    insertion_effect,
    kmer_counts,
    matched_background_loci,
    motif_pair_spacing,
    place_motif,
    pwm_information_content,
    saturation_mutagenesis,
    substitution_effect,
)

__all__ = [*__all__, "cluster_motif_set", "compare_motif_positions", "equal_size_intervals", "motif_enrichment_stats", "motif_from_alignment", "motif_from_consensus", "motif_score_correlation", "promoter_regions", "random_genomic_regions", "read_peak_summits", "scan_best_motif_hits", "summit_intervals_from_narrowpeak", "chromosome_fragment_counts", "fragment_in_peaks_fraction", "fragment_size_distribution", "gene_activity_counts", "merge_peak_sets", "region_motif_enrichment", "regulatory_regions_to_genes", "tss_enrichment", "ablate_region", "call_seqlets", "deletion_effect", "design_greedy_substitution", "extract_loci_sequences", "insertion_effect", "kmer_counts", "matched_background_loci", "motif_pair_spacing", "place_motif", "pwm_information_content", "saturation_mutagenesis", "substitution_effect"]

# Review-pass exports (see docs/api/features)
from .gimmemotifs_review_functions import (
    combine_summit_files,
    equal_length_motif_pair,
    fasta_record_count,
    motif_activity_scores,
    motif_position_ks,
    motif_score_metric,
    peak_coverage_table,
    scan_motifs_to_file,
    scan_region_scores,
)

from .snapatac2_review_functions import (
    add_tf_binding_edges,
    call_peaks_macs3,
    export_coverage_tracks,
    export_fragment_files,
    fragments_from_bam,
    import_base_pair_values,
    import_contact_pairs,
    tenx_atac_metrics,
)

from .tangermeme_review_functions import (
    ablate_annotations,
    annotate_seqlets,
    count_seqlet_annotations,
    delete_region,
    example_to_sequence_coords,
    gapped_kmer_features,
    greedy_marginalize_design,
    local_dinucleotide_shuffle,
    marginalize_annotations,
    one_hot_to_fasta_file,
    pairwise_annotation_counts,
    pisa_attribution,
    place_motif_set,
    pwm_consensus_sequence,
    random_one_hot_sequences,
    randomize_region,
    screen_sequence_design,
    shuffle_region,
    tfmodisco_seqlets,
)

__all__ = [*__all__, "combine_summit_files", "equal_length_motif_pair", "fasta_record_count", "motif_activity_scores", "motif_position_ks", "motif_score_metric", "peak_coverage_table", "scan_motifs_to_file", "scan_region_scores", "add_tf_binding_edges", "call_peaks_macs3", "export_coverage_tracks", "export_fragment_files", "fragments_from_bam", "import_base_pair_values", "import_contact_pairs", "tenx_atac_metrics", "ablate_annotations", "annotate_seqlets", "count_seqlet_annotations", "delete_region", "example_to_sequence_coords", "gapped_kmer_features", "greedy_marginalize_design", "local_dinucleotide_shuffle", "marginalize_annotations", "one_hot_to_fasta_file", "pairwise_annotation_counts", "pisa_attribution", "place_motif_set", "pwm_consensus_sequence", "random_one_hot_sequences", "randomize_region", "screen_sequence_design", "shuffle_region", "tfmodisco_seqlets"]

# Family-pass exports (see docs/api/features)
from .gimmemotifs_family_functions import (
    calculate_partial_motif_roc_auc,
    count_motif_matches_at_fpr,
    normalize_peak_widths,
)

from .tangermeme_family_functions import (
    design_beam_substitution,
    insert_motif,
)

__all__ = [*__all__, "calculate_partial_motif_roc_auc", "count_motif_matches_at_fpr", "normalize_peak_widths", "design_beam_substitution", "insert_motif"]


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .tangermeme_gap_functions import (
    plot_attribution_logo,
    plot_categorical_attribution_scatter,
    plot_pwm_logo,
    plot_sequence_attributions,
    render_interactive_logo_png,
)

__all__ = [*__all__, "plot_attribution_logo", "plot_categorical_attribution_scatter", "plot_pwm_logo", "plot_sequence_attributions", "render_interactive_logo_png"]

# Public integration of the reviewed P00–P45 operation packs.
from .biopython_p18_motif_data_functions import (
    cluster_expression_data,
    motif_gc_content,
    read_affymetrix_cel,
    read_cluster_expression_data,
    read_geo_soft,
    read_motif_collection,
    read_unigene,
)
from .biopython_p19_motif_functions import (
    jaspar_pseudocounts,
    motif_anticonsensus,
    motif_matrices,
    pssm_score_distribution,
    pssm_score_statistics,
    pssm_search,
    pssm_threshold,
    weblogo_service_logo,
)
from .fimo_p17_cli_functions import (
    scan_motifs_with_position_priors,
)
from .gimmemotifs_p28_data_functions import (
    gimmemotifs_annotate_motif_factors,
    gimmemotifs_count_motif_matches,
    gimmemotifs_download_motif_database,
    gimmemotifs_hardmask_fasta,
    gimmemotifs_motif_thresholds,
    gimmemotifs_narrowpeak_to_bed,
    gimmemotifs_prepare_denovo_input,
    gimmemotifs_random_fasta_subset,
    gimmemotifs_rank_motif_stats,
    gimmemotifs_sequence_dict,
    gimmemotifs_write_motif_stats,
)
from .gimmemotifs_p28_maelstrom_functions import (
    gimmemotifs_aggregate_motif_ranks,
    gimmemotifs_plot_maelstrom_heatmap,
    gimmemotifs_predict_motif_activity,
)
from .gimmemotifs_p28_motif_functions import (
    gimmemotifs_average_motifs,
    gimmemotifs_compare_motif_pair,
    gimmemotifs_consensus_scan,
    gimmemotifs_format_motif,
    gimmemotifs_motif_consensus,
    gimmemotifs_motif_ic_similarity,
    gimmemotifs_motif_pcc,
    gimmemotifs_motif_reverse_complement,
    gimmemotifs_motif_score_distribution,
    gimmemotifs_motif_score_range,
    gimmemotifs_motif_shuffle,
    gimmemotifs_pfm_to_ppm,
    gimmemotifs_ppm_pseudocount,
    gimmemotifs_pwm_scan,
    gimmemotifs_sample_motif_sequences,
)
from .snapatac2_p38_functions import (
    aggregate_atac_profiles,
    compute_atac_spectral_embedding,
    filter_atac_cells,
    find_atac_marker_regions,
    import_atac_fragments,
    link_tf_to_gene_network,
    read_meme_motifs,
    scan_dna_with_motif,
    score_atac_doublets,
    score_regulatory_network_edges,
    select_atac_features,
    test_atac_differential_accessibility,
    test_dna_motif_enrichment,
)
from .tangermeme_p41_functions import (
    tangermeme_extract_signal,
    tangermeme_hypothetical_attributions,
    tangermeme_position_entropy,
    tangermeme_unchunk,
)

__all__ = [
    *__all__,
    "aggregate_atac_profiles",
    "cluster_expression_data",
    "compute_atac_spectral_embedding",
    "filter_atac_cells",
    "find_atac_marker_regions",
    "gimmemotifs_aggregate_motif_ranks",
    "gimmemotifs_annotate_motif_factors",
    "gimmemotifs_average_motifs",
    "gimmemotifs_compare_motif_pair",
    "gimmemotifs_consensus_scan",
    "gimmemotifs_count_motif_matches",
    "gimmemotifs_download_motif_database",
    "gimmemotifs_format_motif",
    "gimmemotifs_hardmask_fasta",
    "gimmemotifs_motif_consensus",
    "gimmemotifs_motif_ic_similarity",
    "gimmemotifs_motif_pcc",
    "gimmemotifs_motif_reverse_complement",
    "gimmemotifs_motif_score_distribution",
    "gimmemotifs_motif_score_range",
    "gimmemotifs_motif_shuffle",
    "gimmemotifs_motif_thresholds",
    "gimmemotifs_narrowpeak_to_bed",
    "gimmemotifs_pfm_to_ppm",
    "gimmemotifs_plot_maelstrom_heatmap",
    "gimmemotifs_ppm_pseudocount",
    "gimmemotifs_predict_motif_activity",
    "gimmemotifs_prepare_denovo_input",
    "gimmemotifs_pwm_scan",
    "gimmemotifs_random_fasta_subset",
    "gimmemotifs_rank_motif_stats",
    "gimmemotifs_sample_motif_sequences",
    "gimmemotifs_sequence_dict",
    "gimmemotifs_write_motif_stats",
    "import_atac_fragments",
    "jaspar_pseudocounts",
    "link_tf_to_gene_network",
    "motif_anticonsensus",
    "motif_gc_content",
    "motif_matrices",
    "pssm_score_distribution",
    "pssm_score_statistics",
    "pssm_search",
    "pssm_threshold",
    "read_affymetrix_cel",
    "read_cluster_expression_data",
    "read_geo_soft",
    "read_meme_motifs",
    "read_motif_collection",
    "read_unigene",
    "scan_dna_with_motif",
    "scan_motifs_with_position_priors",
    "score_atac_doublets",
    "score_regulatory_network_edges",
    "select_atac_features",
    "tangermeme_extract_signal",
    "tangermeme_hypothetical_attributions",
    "tangermeme_position_entropy",
    "tangermeme_unchunk",
    "test_atac_differential_accessibility",
    "test_dna_motif_enrichment",
    "weblogo_service_logo",
]

_register_public_verbose(globals())
