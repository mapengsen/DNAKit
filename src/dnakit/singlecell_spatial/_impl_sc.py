"""Matrix-level single-cell and spatial operations (pure standard library unless noted).

Each ``op_*`` function takes the raw ``input`` and an options dict and returns ``(result, warnings)``.
Matrices are genes by cells; see ``_impl_matrix`` for the accepted input forms.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence

from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError

from . import _impl_stats as stats
from . import _impl_svg as svg
from ._impl_de import (
    de_settings,
    differential_expression,
    log1p_columns,
    rank_genes,
    scale_columns,
    two_groups,
)
from ._impl_matrix import (
    CountMatrix,
    labels_for,
    matrix_result,
    numpy_module,
    opt_bool,
    opt_choice,
    opt_float,
    opt_int,
    opt_str,
    option_labels,
    parse_coordinates,
    parse_matrix,
    require_option,
    scipy_optimize_module,
)

Result = tuple[object, list[str]]

DEFAULT_SC_TARGET = 1e4
_MAX_LISTED = 1000


def _mito_prefixes(options: Mapping[str, object]) -> tuple[str, ...]:
    raw = options.get("mito_prefix", ("MT-", "mt-"))
    if isinstance(raw, str):
        raw = (raw,)
    if (
        not isinstance(raw, Sequence)
        or not raw
        or any(not isinstance(p, str) or not p for p in raw)
    ):
        raise ConfigurationError(
            "options['mito_prefix'] must be a non-empty string or list of strings.",
            code="INVALID_OPTION",
            hint="e.g. 'MT-'.",
        )
    return tuple(raw)


def _spatial_summary(coords: Sequence[tuple[float, float]]) -> dict[str, object]:
    xs, ys = [c[0] for c in coords], [c[1] for c in coords]
    return {
        "n_spots_with_coordinates": len(coords),
        "x_range": [min(xs), max(xs)],
        "y_range": [min(ys), max(ys)],
    }


# ---------------------------------------------------------------------------------------------
# read / validate / qc / normalize
# ---------------------------------------------------------------------------------------------


def op_read(value: object, options: dict[str, object], *, spatial: bool = False) -> Result:
    cm = parse_matrix(value, options)
    totals = cm.column_totals()
    nnz = sum(1 for row in cm.values for x in row if x != 0)
    result: dict[str, object] = {
        "n_genes": cm.n_genes,
        "n_cells": cm.n_cells,
        "total_counts": math.fsum(totals),
        "nonzero_entries": nnz,
        "sparsity": 1.0 - nnz / (cm.n_genes * cm.n_cells),
        "genes_preview": cm.genes[:10],
        "cells_preview": cm.cells[:10],
        "orientation": "genes_by_cells",
    }
    warnings: list[str] = []
    if spatial and options.get("coordinates") is not None:
        coords = parse_coordinates(options["coordinates"], cm.cells, options)
        result["coordinates"] = {cell: list(xy) for cell, xy in zip(cm.cells, coords, strict=False)}
        result.update(_spatial_summary(coords))
    elif spatial:
        warnings.append("No options['coordinates'] given: only the count matrix was read.")
    limit = opt_int(options, "max_return_elements", 10_000, 0, 10**8)
    if cm.n_genes * cm.n_cells <= limit:
        result["matrix"] = cm.as_dict()
    else:
        warnings.append(
            f"Matrix ({cm.n_genes * cm.n_cells} elements) exceeds max_return_elements={limit}; values omitted."
        )
    return result, warnings


def op_validate(value: object, options: dict[str, object], *, spatial: bool = False) -> Result:
    if value is None:
        raise InputFormatError(
            "input is required.", code="EMPTY_INPUT", hint="Provide a count matrix."
        )
    errors: list[str] = []
    warnings: list[str] = []
    result: dict[str, object] = {}
    try:
        cm = parse_matrix(value, options, allow_negative=True, check_names=False)
    except InputFormatError as exc:
        return {"valid": False, "errors": [str(exc)], "warnings": []}, []
    require_integer = opt_bool(options, "require_integer", True)
    result.update(n_genes=cm.n_genes, n_cells=cm.n_cells)
    for what, names in (("gene", cm.genes), ("cell", cm.cells)):
        duplicates = (
            sorted({n for n in names if names.count(n) > 1})
            if len(set(names)) != len(names)
            else []
        )
        if duplicates:
            errors.append(f"duplicate {what} names: {duplicates[:5]}")
    negatives = sum(1 for row in cm.values for x in row if x < 0)
    if negatives:
        errors.append(f"{negatives} negative values")
    if require_integer:
        fractional = sum(1 for row in cm.values for x in row if x != math.floor(x))
        if fractional:
            errors.append(
                f"{fractional} non-integer values (raw counts must be integers; set require_integer=False for normalised data)"
            )
    totals = cm.column_totals()
    empty_cells = [c for c, t in zip(cm.cells, totals, strict=False) if t == 0]
    empty_genes = [g for g, row in zip(cm.genes, cm.values, strict=False) if not any(row)]
    if empty_cells:
        warnings.append(
            f"{len(empty_cells)} cells have zero total counts (first: {empty_cells[0]!r})"
        )
    if empty_genes:
        warnings.append(
            f"{len(empty_genes)} genes are zero in every cell (first: {empty_genes[0]!r})"
        )
    result.update(empty_cells=len(empty_cells), empty_genes=len(empty_genes))
    if spatial:
        if options.get("coordinates") is None:
            errors.append("options['coordinates'] is missing for a spatial dataset")
        else:
            try:
                coords = parse_coordinates(options["coordinates"], cm.cells, options)
                result.update(_spatial_summary(coords))
            except (DNAKitError, ValueError) as exc:
                errors.append(f"coordinates: {exc}")
    result.update(valid=not errors, errors=errors, warnings=warnings)
    return result, []


def op_qc(value: object, options: dict[str, object], *, spatial: bool = False) -> Result:
    cm = parse_matrix(value, options)
    prefixes = _mito_prefixes(options)
    min_genes = opt_int(options, "min_genes", None, 0)
    min_counts = opt_float(options, "min_counts", None, 0.0)
    max_counts = opt_float(options, "max_counts", None, 0.0)
    max_pct_mito = opt_float(options, "max_pct_mito", None, 0.0, 100.0)
    min_cells = opt_int(options, "min_cells", None, 0)
    mito_rows = [i for i, g in enumerate(cm.genes) if g.startswith(prefixes)]
    totals = cm.column_totals()
    n_detected = [sum(1 for row in cm.values if row[j] > 0) for j in range(cm.n_cells)]
    mito_totals = [math.fsum(cm.values[i][j] for i in mito_rows) for j in range(cm.n_cells)]
    pct_mito = [100.0 * m / t if t > 0 else 0.0 for m, t in zip(mito_totals, totals, strict=False)]
    per_cell: list[dict[str, object]] = []
    failed: list[str] = []
    for j, cell in enumerate(cm.cells):
        reasons = []
        if min_genes is not None and n_detected[j] < min_genes:
            reasons.append("min_genes")
        if min_counts is not None and totals[j] < min_counts:
            reasons.append("min_counts")
        if max_counts is not None and totals[j] > max_counts:
            reasons.append("max_counts")
        if max_pct_mito is not None and pct_mito[j] > max_pct_mito:
            reasons.append("max_pct_mito")
        if reasons:
            failed.append(cell)
        per_cell.append(
            {
                "cell": cell,
                "total_counts": totals[j],
                "n_genes_detected": n_detected[j],
                "pct_mito": pct_mito[j],
                "fails": reasons,
            }
        )
    cells_per_gene = [sum(1 for x in row if x > 0) for row in cm.values]
    low_genes = [
        g
        for g, n in zip(cm.genes, cells_per_gene, strict=False)
        if min_cells is not None and n < min_cells
    ]
    result: dict[str, object] = {
        "n_cells": cm.n_cells,
        "n_genes": cm.n_genes,
        "n_mito_genes": len(mito_rows),
        "median_total_counts": stats.median(totals),
        "median_genes_detected": stats.median([float(n) for n in n_detected]),
        "median_pct_mito": stats.median(pct_mito),
        "genes_detected_in_any_cell": sum(1 for n in cells_per_gene if n > 0),
        "thresholds": {
            "min_genes": min_genes,
            "min_counts": min_counts,
            "max_counts": max_counts,
            "max_pct_mito": max_pct_mito,
            "min_cells": min_cells,
        },
        "n_cells_failing": len(failed),
        "failing_cells": failed[:_MAX_LISTED],
        "n_genes_below_min_cells": len(low_genes),
        "genes_below_min_cells": low_genes[:_MAX_LISTED],
        "per_cell": per_cell[:_MAX_LISTED],
        "per_cell_truncated": len(per_cell) > _MAX_LISTED,
    }
    warnings = []
    if not mito_rows:
        warnings.append(f"No gene starts with {list(prefixes)}: pct_mito is 0 for every cell.")
    if spatial and options.get("coordinates") is not None:
        result.update(
            _spatial_summary(parse_coordinates(options["coordinates"], cm.cells, options))
        )
    return result, warnings


def op_normalize(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    method = opt_choice(options, "method", "log1p", ("log1p", "total", "cpm"))
    target = (
        1e6
        if method == "cpm"
        else opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True)
    )
    out = log1p_columns(cm, target) if method == "log1p" else scale_columns(cm, target)
    totals = cm.column_totals()
    warnings = []
    zero = sum(1 for t in totals if t == 0)
    if zero:
        warnings.append(f"{zero} cells have zero total counts and stay all-zero.")
    extra = {
        "method": method,
        "target_sum": target,
        "library_sizes": dict(zip(cm.cells, totals, strict=False)),
        "log_base": "e" if method == "log1p" else None,
    }
    return matrix_result(out, options, extra), warnings


# ---------------------------------------------------------------------------------------------
# differential features, markers, perturbations
# ---------------------------------------------------------------------------------------------


def op_find_differential(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    labels = option_labels(cm, options, "groups", "a label per cell, e.g. cluster or condition")
    group1, group2, first, second = two_groups(cm, labels, options)
    settings = de_settings(options, DEFAULT_SC_TARGET)
    rows = differential_expression(
        cm,
        first,
        second,
        method=settings["method"],
        target_sum=settings["target_sum"],
        pseudocount=settings["pseudocount"],
    )  # type: ignore[arg-type]
    result = {
        "group1": group1,
        "group2": group2,
        "n_group1": len(first),
        "n_group2": len(second),
        "method": settings["method"],
        "normalization": settings["normalization"],
        "target_sum": settings["target_sum"],
        "pseudocount": settings["pseudocount"],
        "n_genes_tested": len(rows),
        "genes": rank_genes(rows, top_n=settings["top_n"], min_pct=settings["min_pct"]),  # type: ignore[arg-type]
    }
    return result, []


def op_calculate_markers(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    labels = option_labels(cm, options, "groups", "a cluster/cell-type label per cell")
    settings = de_settings(options, DEFAULT_SC_TARGET)
    only_positive = opt_bool(options, "only_positive", True)
    warnings: list[str] = []
    markers: dict[str, object] = {}
    for group in sorted(set(labels)):
        first = [i for i, label in enumerate(labels) if label == group]
        rest = [i for i, label in enumerate(labels) if label != group]
        if not rest:
            raise InputFormatError(
                "calculate_markers needs at least two groups.",
                code="GROUPS_TOO_FEW",
                hint="Provide labels with two or more distinct values.",
            )
        if settings["method"] == "welch" and (len(first) < 2 or len(rest) < 2):
            warnings.append(
                f"Group {group!r} skipped: Welch's t-test needs >= 2 cells on each side."
            )
            continue
        rows = differential_expression(
            cm,
            first,
            rest,
            method=settings["method"],
            target_sum=settings["target_sum"],
            pseudocount=settings["pseudocount"],
        )  # type: ignore[arg-type]
        markers[group] = {
            "n_cells": len(first),
            "genes": rank_genes(
                rows,
                top_n=settings["top_n"],
                min_pct=settings["min_pct"],
                only_positive=only_positive,
            ),  # type: ignore[arg-type]
        }
    if not markers:
        raise InputFormatError(
            "No group had enough cells to test.",
            code="GROUPS_TOO_SMALL",
            hint="Use method='wilcoxon' or merge small clusters.",
        )
    return {
        "method": settings["method"],
        "only_positive": only_positive,
        "comparison": "each group versus all other cells",
        "markers": markers,
    }, warnings


def _perturbation_groups(
    cm: CountMatrix, options: Mapping[str, object]
) -> tuple[list[str], str, list[int]]:
    labels = option_labels(cm, options, "perturbations", "a perturbation/guide label per cell")
    control = opt_str(options, "control", "control")
    controls = [i for i, label in enumerate(labels) if label == control]
    if not controls:
        raise InputFormatError(
            f"No cell carries the control label {control!r}.",
            code="NO_CONTROL_CELLS",
            hint="Set options['control'] to the label of control cells.",
        )
    return labels, control, controls


def op_perturbation_effect(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    labels, control, controls = _perturbation_groups(cm, options)
    settings = de_settings(options, DEFAULT_SC_TARGET)
    warnings: list[str] = []
    effects: dict[str, object] = {}
    for label in sorted(set(labels) - {control}):
        members = [i for i, group in enumerate(labels) if group == label]
        if settings["method"] == "welch" and (len(members) < 2 or len(controls) < 2):
            warnings.append(
                f"Perturbation {label!r} skipped: Welch's t-test needs >= 2 cells in each group."
            )
            continue
        rows = differential_expression(
            cm,
            members,
            controls,
            method=settings["method"],
            target_sum=settings["target_sum"],
            pseudocount=settings["pseudocount"],
        )  # type: ignore[arg-type]
        effects[label] = {
            "n_cells": len(members),
            "genes": rank_genes(rows, top_n=settings["top_n"], min_pct=settings["min_pct"]),
        }  # type: ignore[arg-type]
    if not effects:
        raise InputFormatError(
            "No perturbation had enough cells to test.",
            code="GROUPS_TOO_SMALL",
            hint="Use method='wilcoxon' or add cells.",
        )
    return {
        "control": control,
        "n_control_cells": len(controls),
        "comparison": "each perturbation versus control cells",
        "method": settings["method"],
        "effects": effects,
    }, warnings


def op_analyze_crispr(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    labels, control, controls = _perturbation_groups(cm, options)
    targets_raw = options.get("targets") or {}
    if not isinstance(targets_raw, Mapping):
        raise ConfigurationError(
            "options['targets'] must map guide label -> target gene.",
            code="INVALID_OPTION",
            hint="e.g. {'sgTP53_1': 'TP53'}.",
        )
    target_sum = opt_float(
        options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True
    )
    scaled = scale_columns(cm, target_sum)
    index = cm.gene_index()
    guides: dict[str, object] = {}
    warnings: list[str] = []
    for label in sorted(set(labels) - {control}):
        members = [i for i, group in enumerate(labels) if group == label]
        target = str(targets_raw.get(label, label))
        entry: dict[str, object] = {"n_cells": len(members), "target_gene": target}
        if target in index:
            row = scaled.values[index[target]]
            mean_guide = stats.mean([row[j] for j in members])
            mean_ctrl = stats.mean([row[j] for j in controls])
            entry.update(
                target_in_matrix=True,
                mean_normalized_guide=mean_guide,
                mean_normalized_control=mean_ctrl,
                log2fc=math.log2(mean_guide + 1.0) - math.log2(mean_ctrl + 1.0),
                percent_knockdown=100.0 * (1.0 - mean_guide / mean_ctrl) if mean_ctrl > 0 else None,
            )
        else:
            entry["target_in_matrix"] = False
            warnings.append(f"Target gene {target!r} of {label!r} is not in the matrix.")
        guides[label] = entry
    assigned = defaultdict(int)
    for label in labels:
        assigned[label] += 1
    return {
        "control": control,
        "n_control_cells": len(controls),
        "cells_per_label": dict(sorted(assigned.items())),
        "target_sum": target_sum,
        "knockdown_definition": "100 * (1 - mean CP10k(guide cells) / mean CP10k(control cells)) of the target gene",
        "guides": guides,
    }, warnings


# ---------------------------------------------------------------------------------------------
# marker-based annotation, cell cycle, cell-cell interaction
# ---------------------------------------------------------------------------------------------


def _marker_sets(raw: object) -> dict[str, list[str]]:
    if not isinstance(raw, Mapping) or not raw:
        raise ConfigurationError(
            "options['markers'] must map cell type -> list of marker genes.",
            code="INVALID_OPTION",
            hint="e.g. {'T cell': ['CD3D', 'CD3E']}.",
        )
    result: dict[str, list[str]] = {}
    for name, genes in raw.items():
        if isinstance(genes, (str, bytes)) or not isinstance(genes, Sequence) or not genes:
            raise ConfigurationError(
                f"Marker list of {name!r} must be a non-empty list of genes.",
                code="INVALID_OPTION",
                hint="Give gene symbols as strings.",
            )
        result[str(name)] = [str(g) for g in genes]
    return result


def op_annotate_markers(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    markers = _marker_sets(options.get("markers"))
    min_score = opt_float(options, "min_score", 0.0)
    log = log1p_columns(
        cm, opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True)
    )
    z_rows: dict[str, list[float]] = {}
    index = cm.gene_index()
    warnings: list[str] = []
    for name, genes in markers.items():
        present = [g for g in genes if g in index]
        if len(present) < len(genes):
            warnings.append(
                f"{name}: markers missing from the matrix: {[g for g in genes if g not in index]}"
            )
        if not present:
            raise InputFormatError(
                f"None of the markers of {name!r} is in the matrix.",
                code="MARKERS_NOT_FOUND",
                hint="Check gene symbols and case.",
            )
        for gene in present:
            if gene not in z_rows:
                row = log.values[index[gene]]
                centre = stats.mean(row)
                sd = math.sqrt(stats.variance(row, ddof=0))
                z_rows[gene] = [(x - centre) / sd if sd > 0 else 0.0 for x in row]
    scores = {
        name: [stats.mean([z_rows[g][j] for g in genes if g in z_rows]) for j in range(cm.n_cells)]
        for name, genes in markers.items()
    }
    names = list(markers)
    assignments: dict[str, str] = {}
    confidences: dict[str, float] = {}
    for j, cell in enumerate(cm.cells):
        best = max(range(len(names)), key=lambda k: (scores[names[k]][j], -k))
        top = scores[names[best]][j]
        assignments[cell] = names[best] if top > min_score else "Unassigned"
        confidences[cell] = top
    counts: dict[str, int] = defaultdict(int)
    for label in assignments.values():
        counts[label] += 1
    return {
        "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
        "assignments": assignments,
        "best_score": confidences,
        "counts": dict(sorted(counts.items())),
        "scores": {
            name: dict(zip(cm.cells, values, strict=False)) for name, values in scores.items()
        },
    }, warnings


def op_cell_cycle(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    s_genes = [str(g) for g in require_option(options, "s_genes", "list of S-phase genes")]  # type: ignore[union-attr]
    g2m_genes = [str(g) for g in require_option(options, "g2m_genes", "list of G2/M-phase genes")]  # type: ignore[union-attr]
    n_bins = opt_int(options, "n_bins", 25, 1, 1000)
    ctrl_size = opt_int(options, "ctrl_size", 50, 1, 10_000)
    log = log1p_columns(
        cm, opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True)
    )
    index = cm.gene_index()
    gene_means = [stats.mean(row) for row in log.values]
    order = sorted(range(cm.n_genes), key=lambda i: (gene_means[i], i))
    bins: dict[int, list[int]] = defaultdict(list)
    for rank, _gene_i in enumerate(order):
        bins[rank * n_bins // cm.n_genes].append(rank)  # store the rank positions
    rank_of = {gene_i: rank for rank, gene_i in enumerate(order)}
    warnings: list[str] = []

    def score(genes: Sequence[str], label: str) -> list[float]:
        present = [index[g] for g in genes if g in index]
        if len(present) < len(genes):
            warnings.append(f"{label}: {len(genes) - len(present)} genes missing from the matrix.")
        if not present:
            raise InputFormatError(
                f"None of the {label} genes is in the matrix.",
                code="MARKERS_NOT_FOUND",
                hint="Check gene symbols and case.",
            )
        excluded = set(present)
        controls: set[int] = set()
        for gene_i in present:
            rank = rank_of[gene_i]
            same_bin = [
                order[r] for r in bins[rank * n_bins // cm.n_genes] if order[r] not in excluded
            ]
            same_bin.sort(key=lambda other: (abs(rank_of[other] - rank), rank_of[other]))
            controls.update(same_bin[:ctrl_size])
        if not controls:
            raise InputFormatError(
                f"No control genes available for {label}.",
                code="NO_CONTROL_GENES",
                hint="Use fewer bins or more genes.",
            )
        return [
            stats.mean([log.values[i][j] for i in present])
            - stats.mean([log.values[i][j] for i in controls])
            for j in range(cm.n_cells)
        ]

    s_scores = score(s_genes, "S")
    g2m_scores = score(g2m_genes, "G2M")
    phases: dict[str, str] = {}
    for j, cell in enumerate(cm.cells):
        if s_scores[j] < 0 and g2m_scores[j] < 0:
            phases[cell] = "G1"
        else:
            phases[cell] = "S" if s_scores[j] >= g2m_scores[j] else "G2M"
    return {
        "definition": "score = mean ln(1+CP10k) of the gene set minus mean of expression-matched control genes (genes in the same mean-expression rank bin, nearest rank first); phase = G1 if both scores < 0, else the larger of S and G2M",
        "s_score": dict(zip(cm.cells, s_scores, strict=False)),
        "g2m_score": dict(zip(cm.cells, g2m_scores, strict=False)),
        "phase": phases,
        "phase_counts": {p: list(phases.values()).count(p) for p in ("G1", "S", "G2M")},
    }, warnings


def _pairs(raw: object) -> list[tuple[str, str]]:
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or not raw:
        raise ConfigurationError(
            "options['pairs'] must be a list of [ligand, receptor] pairs.",
            code="INVALID_OPTION",
            hint="e.g. [['CXCL12', 'CXCR4']].",
        )
    pairs = []
    for item in raw:
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 2:
            raise ConfigurationError(
                f"Bad ligand-receptor pair {item!r}.",
                code="INVALID_OPTION",
                hint="Use [ligand, receptor].",
            )
        pairs.append((str(item[0]), str(item[1])))
    return pairs


def op_cell_interactions(
    value: object, options: dict[str, object], *, spatial: bool = False
) -> Result:
    cm = parse_matrix(value, options)
    labels = option_labels(cm, options, "cell_types", "a cell-type label per cell")
    pairs = _pairs(options.get("pairs"))
    min_pct = opt_float(options, "min_pct", 0.0, 0.0, 1.0)
    top_n = opt_int(options, "top_n", 100, 1, 10**6)
    log = log1p_columns(
        cm, opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True)
    )
    index = cm.gene_index()
    types = sorted(set(labels))
    members = {t: [j for j, group in enumerate(labels) if group == t] for t in types}
    warnings: list[str] = []
    interactions: list[dict[str, object]] = []
    radius = (
        opt_float(options, "radius", None, 0.0, math.inf, exclusive_min=True) if spatial else None
    )
    contacts: dict[tuple[str, str], list[tuple[int, int]]] = {}
    if spatial:
        if cm.n_cells > 5000:
            raise InputFormatError(
                "Spatial interactions support up to 5000 cells/spots.",
                code="INPUT_TOO_LARGE",
                hint="Subsample or split the tissue.",
            )
        if radius is None:
            raise ConfigurationError(
                "options['radius'] is required for spatial interactions.",
                code="MISSING_OPTION",
                hint="Maximum sender-receiver distance in coordinate units.",
            )
        coords = parse_coordinates(options.get("coordinates"), cm.cells, options)
        for a in range(cm.n_cells):
            for b in range(cm.n_cells):
                if a != b and math.dist(coords[a], coords[b]) <= radius:
                    contacts.setdefault((labels[a], labels[b]), []).append((a, b))
    for ligand, receptor in pairs:
        if ligand not in index or receptor not in index:
            warnings.append(f"Pair {ligand}-{receptor} skipped: gene missing from the matrix.")
            continue
        lrow, rrow = log.values[index[ligand]], log.values[index[receptor]]
        raw_l, raw_r = cm.values[index[ligand]], cm.values[index[receptor]]
        for sender in types:
            for receiver in types:
                if spatial:
                    cell_pairs = contacts.get((sender, receiver), [])
                    if not cell_pairs:
                        continue
                    ligand_mean = stats.mean([lrow[a] for a, _ in cell_pairs])
                    receptor_mean = stats.mean([rrow[b] for _, b in cell_pairs])
                    score = stats.mean([lrow[a] * rrow[b] for a, b in cell_pairs])
                    pct_l = sum(1 for a, _ in cell_pairs if raw_l[a] > 0) / len(cell_pairs)
                    pct_r = sum(1 for _, b in cell_pairs if raw_r[b] > 0) / len(cell_pairs)
                    extra: dict[str, object] = {"n_cell_pairs": len(cell_pairs)}
                else:
                    ligand_mean = stats.mean([lrow[j] for j in members[sender]])
                    receptor_mean = stats.mean([rrow[j] for j in members[receiver]])
                    score = ligand_mean * receptor_mean
                    pct_l = sum(1 for j in members[sender] if raw_l[j] > 0) / len(members[sender])
                    pct_r = sum(1 for j in members[receiver] if raw_r[j] > 0) / len(
                        members[receiver]
                    )
                    extra = {}
                if pct_l < min_pct or pct_r < min_pct:
                    continue
                interactions.append(
                    {
                        "ligand": ligand,
                        "receptor": receptor,
                        "sender": sender,
                        "receiver": receiver,
                        "ligand_mean": ligand_mean,
                        "receptor_mean": receptor_mean,
                        "score": score,
                        "pct_ligand": pct_l,
                        "pct_receptor": pct_r,
                        **extra,
                    }
                )
    interactions.sort(
        key=lambda r: (-r["score"], r["ligand"], r["receptor"], r["sender"], r["receiver"])
    )  # type: ignore[operator]
    definition = (
        "score = mean ln(1+CP10k) of the ligand in sender cells x mean of the receptor in receiver cells"
        if not spatial
        else "score = mean over sender-receiver cell pairs within `radius` of ln(1+CP10k)(ligand in sender) x ln(1+CP10k)(receptor in receiver)"
    )
    return {
        "definition": definition,
        "no_significance_test": True,
        "n_interactions": len(interactions),
        "interactions": interactions[:top_n],
    }, warnings


# ---------------------------------------------------------------------------------------------
# ambient RNA, hashtag demultiplexing, deconvolution, integration
# ---------------------------------------------------------------------------------------------


def op_remove_ambient(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    rho = opt_float(options, "contamination", None, 0.0, 1.0)
    if rho is None or rho >= 1.0:
        raise ConfigurationError(
            "options['contamination'] (fraction in [0, 1)) is required.",
            code="MISSING_OPTION",
            hint="The ambient fraction of every cell's counts, e.g. 0.05 (SoupX 'rho').",
        )
    totals = cm.column_totals()
    profile_raw = options.get("ambient_profile")
    empty_raw = options.get("empty_cells")
    empty_max = opt_float(options, "empty_max_counts", None, 0.0)
    warnings: list[str] = []
    keep = list(range(cm.n_cells))
    if profile_raw is not None:
        if not isinstance(profile_raw, Mapping):
            raise ConfigurationError(
                "options['ambient_profile'] must map gene -> weight.",
                code="INVALID_OPTION",
                hint="Weights are normalised to sum to 1.",
            )
        weights = [max(0.0, float(profile_raw.get(g, 0.0))) for g in cm.genes]
        total_w = math.fsum(weights)
        if total_w <= 0:
            raise InputFormatError(
                "ambient_profile has no weight on any gene of the matrix.",
                code="EMPTY_PROFILE",
                hint="Use gene names from the matrix.",
            )
        profile = [w / total_w for w in weights]
        source = "ambient_profile"
    else:
        if empty_raw is not None:
            names = [str(c) for c in empty_raw]  # type: ignore[union-attr]
            idx = cm.cell_index()
            unknown = [c for c in names if c not in idx]
            if unknown:
                raise InputFormatError(
                    f"empty_cells not in the matrix: {unknown[:3]}",
                    code="UNKNOWN_CELLS",
                    hint="Use column names of the matrix.",
                )
            empty = sorted(idx[c] for c in names)
        elif empty_max is not None:
            empty = [j for j, t in enumerate(totals) if t <= empty_max]
        else:
            raise ConfigurationError(
                "Provide options['ambient_profile'], options['empty_cells'] or options['empty_max_counts'].",
                code="MISSING_OPTION",
                hint="The ambient profile is estimated from empty droplets.",
            )
        if not empty:
            raise InputFormatError(
                "No empty droplets selected.",
                code="NO_EMPTY_CELLS",
                hint="Raise empty_max_counts or pass empty_cells.",
            )
        pooled = [math.fsum(cm.values[i][j] for j in empty) for i in range(cm.n_genes)]
        total_p = math.fsum(pooled)
        if total_p <= 0:
            raise InputFormatError(
                "The empty droplets contain no counts.",
                code="EMPTY_PROFILE",
                hint="Select other droplets.",
            )
        profile = [x / total_p for x in pooled]
        keep = [j for j in range(cm.n_cells) if j not in set(empty)]
        source = "empty droplets"
        if not keep:
            raise InputFormatError(
                "Every droplet was treated as empty.",
                code="NO_CELLS",
                hint="Lower empty_max_counts.",
            )
    rounding = opt_bool(options, "round", False)
    corrected = []
    for i in range(cm.n_genes):
        row = []
        for j in keep:
            value_ = max(0.0, cm.values[i][j] - rho * totals[j] * profile[i])
            row.append(float(round(value_)) if rounding else value_)
        corrected.append(row)
    out = CountMatrix(list(cm.genes), [cm.cells[j] for j in keep], corrected)
    removed = math.fsum(totals[j] for j in keep) - math.fsum(out.column_totals())
    extra = {
        "definition": "corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); simplified subtractive SoupX-style correction with a user-given contamination fraction",
        "contamination": rho,
        "ambient_source": source,
        "ambient_profile": dict(zip(cm.genes, profile, strict=False))
        if cm.n_genes <= 1000
        else None,
        "counts_removed": removed,
        "excluded_empty_droplets": [cm.cells[j] for j in range(cm.n_cells) if j not in set(keep)],
    }
    return matrix_result(out, options, extra), warnings


def op_demultiplex(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    min_counts = opt_float(options, "min_counts", 5.0, 0.0)
    ratio = opt_float(options, "ratio", 3.0, 1.0)
    if cm.n_genes < 2:
        raise InputFormatError(
            "Demultiplexing needs at least two hashtag rows.",
            code="TOO_FEW_HASHTAGS",
            hint="Rows are hashtags, columns are cells.",
        )
    assignments: dict[str, str] = {}
    for j, cell in enumerate(cm.cells):
        counts = sorted(
            ((cm.values[i][j], i) for i in range(cm.n_genes)), key=lambda item: (-item[0], item[1])
        )
        (c1, i1), c2 = counts[0], counts[1][0]
        if c1 < min_counts:
            assignments[cell] = "Negative"
        elif c2 >= min_counts and c1 < ratio * c2:
            assignments[cell] = "Doublet"
        else:
            assignments[cell] = cm.genes[i1]
    tally: dict[str, int] = defaultdict(int)
    for label in assignments.values():
        tally[label] += 1
    return {
        "definition": f"Negative if top hashtag count < {min_counts}; Doublet if second count >= {min_counts} and top < {ratio} x second; otherwise Singlet of the top hashtag. Simple count-threshold rule, not a mixture-model (HTODemux) fit",
        "min_counts": min_counts,
        "ratio": ratio,
        "assignments": assignments,
        "counts": dict(sorted(tally.items())),
    }, []


def _signature(options: Mapping[str, object], mixtures: CountMatrix) -> CountMatrix:
    if options.get("signature") is not None:
        return parse_matrix(
            options["signature"],
            {k: v for k, v in options.items() if k not in {"genes", "cells", "orientation"}}
            | {"orientation": "genes_by_cells"},
        )
    if options.get("reference") is None:
        raise ConfigurationError(
            "Provide options['signature'] (genes x cell types) or options['reference'] + options['reference_labels'].",
            code="MISSING_OPTION",
            hint="A signature matrix has one column per cell type.",
        )
    ref = parse_matrix(
        options["reference"],
        {k: v for k, v in options.items() if k in {"delimiter", "max_matrix_elements"}},
    )
    labels = labels_for(
        ref.cells,
        require_option(options, "reference_labels", "label of each reference cell"),
        "options['reference_labels']",
    )
    scaled = scale_columns(ref, 1.0)
    types = sorted(set(labels))
    columns = [[j for j, group in enumerate(labels) if group == t] for t in types]
    values = [
        [stats.mean([scaled.values[i][j] for j in cols]) for cols in columns]
        for i in range(ref.n_genes)
    ]
    return CountMatrix(ref.genes, types, values)


def op_deconvolve(value: object, options: dict[str, object]) -> Result:
    mixtures = parse_matrix(value, options)
    signature = _signature(options, mixtures)
    optimize = scipy_optimize_module()
    np = numpy_module()
    shared = [g for g in mixtures.genes if g in signature.gene_index()]
    if len(shared) < len(signature.cells):
        raise InputFormatError(
            f"Only {len(shared)} genes are shared between mixtures and signature.",
            code="TOO_FEW_SHARED_GENES",
            hint="Use the same gene identifiers.",
        )
    mi, si = mixtures.gene_index(), signature.gene_index()
    mix = CountMatrix(shared, mixtures.cells, [mixtures.values[mi[g]] for g in shared])
    sig = CountMatrix(shared, signature.cells, [signature.values[si[g]] for g in shared])
    normalize = opt_bool(options, "normalize", True)
    if normalize:
        mix, sig = scale_columns(mix, 1.0), scale_columns(sig, 1.0)
    design = np.array(sig.values, dtype=float)
    fractions: dict[str, dict[str, float]] = {}
    rmse: dict[str, float] = {}
    for j, cell in enumerate(mix.cells):
        target = np.array([row[j] for row in mix.values], dtype=float)
        coefficients, residual = optimize.nnls(design, target)
        total = float(coefficients.sum())
        fractions[cell] = {
            t: (float(c) / total if total > 0 else 0.0)
            for t, c in zip(sig.cells, coefficients, strict=False)
        }
        rmse[cell] = float(residual) / math.sqrt(len(shared))
    return {
        "definition": "non-negative least squares of each mixture column on the signature columns (scipy.optimize.nnls); fractions = coefficients / sum(coefficients); columns are library-size normalised to 1 first"
        if normalize
        else "non-negative least squares on the raw values",
        "cell_types": list(sig.cells),
        "n_shared_genes": len(shared),
        "fractions": fractions,
        "rmse": rmse,
    }, []


def op_integrate(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options, allow_negative=True)
    batches = option_labels(cm, options, "batches", "a batch label per cell")
    log_norm = opt_bool(options, "log_normalize", True)
    data = (
        log1p_columns(
            cm,
            opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True),
        )
        if log_norm
        else cm
    )
    names = sorted(set(batches))
    if len(names) < 2:
        raise InputFormatError(
            "Integration needs at least two batches.",
            code="BATCHES_TOO_FEW",
            hint="Provide labels with two or more distinct values.",
        )
    members = {b: [j for j, group in enumerate(batches) if group == b] for b in names}
    corrected = []
    for row in data.values:
        grand = stats.mean(row)
        out = list(row)
        for b in names:
            offset = stats.mean([row[j] for j in members[b]]) - grand
            for j in members[b]:
                out[j] = row[j] - offset
        corrected.append(out)
    result_cm = CountMatrix(list(cm.genes), list(cm.cells), corrected)
    extra = {
        "definition": "per-gene batch mean-centring: x' = x - mean(gene, batch) + mean(gene, all cells), on ln(1+CP10k) values"
        if log_norm
        else "per-gene batch mean-centring on the supplied values",
        "batches": {b: len(members[b]) for b in names},
        "scale": "ln(1+CP10k)" if log_norm else "as supplied",
        "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    }
    return matrix_result(result_cm, options, extra), []


# ---------------------------------------------------------------------------------------------
# RNA velocity (steady-state) and CNV inference
# ---------------------------------------------------------------------------------------------


def op_velocity(value: object, options: dict[str, object]) -> Result:
    spliced = parse_matrix(value, options)
    if options.get("unspliced") is None:
        raise ConfigurationError(
            "options['unspliced'] (unspliced count matrix) is required.",
            code="MISSING_OPTION",
            hint="Same genes and cells as the spliced input.",
        )
    unspliced = parse_matrix(
        options["unspliced"],
        {k: v for k, v in options.items() if k in {"orientation", "delimiter", "genes", "cells"}},
    )
    if unspliced.genes != spliced.genes or unspliced.cells != spliced.cells:
        raise InputFormatError(
            "The unspliced matrix must have the same genes and cells as the spliced matrix, in the same order.",
            code="LAYER_MISMATCH",
            hint="Reindex the layers first.",
        )
    quantile = opt_float(options, "quantile", 0.05, 0.0, 1.0, exclusive_min=True)
    if opt_bool(options, "normalize", False):
        spliced = scale_columns(spliced, DEFAULT_SC_TARGET)
        unspliced = scale_columns(unspliced, DEFAULT_SC_TARGET)
    gamma: dict[str, float | None] = {}
    velocity = []
    warnings: list[str] = []
    for g, (srow, urow) in enumerate(zip(spliced.values, unspliced.values, strict=False)):
        high_s = stats.quantile(srow, 1.0 - quantile)
        high_u = stats.quantile(urow, 1.0 - quantile)
        chosen = [j for j in range(len(srow)) if srow[j] >= high_s or urow[j] >= high_u]
        sss = math.fsum(srow[j] ** 2 for j in chosen)
        if sss == 0.0:
            gamma[spliced.genes[g]] = None
            velocity.append([0.0] * len(srow))
            warnings.append(
                f"{spliced.genes[g]}: no spliced signal in the fitted cells; velocity set to 0."
            )
            continue
        value_g = math.fsum(srow[j] * urow[j] for j in chosen) / sss
        gamma[spliced.genes[g]] = value_g
        velocity.append([urow[j] - value_g * srow[j] for j in range(len(srow))])
    out = CountMatrix(list(spliced.genes), list(spliced.cells), velocity)
    extra = {
        "definition": "steady-state model (La Manno et al. 2018): gamma_g = sum(u*s)/sum(s^2) over the cells whose s or u is in the top `quantile` of that gene (least squares through the origin); velocity = u - gamma*s",
        "quantile": quantile,
        "gamma": gamma,
    }
    return matrix_result(out, options, extra), warnings[:20]


def _chrom_key(chrom: str) -> tuple[int, int, str]:
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, 0, name)


def _gene_positions(raw: object) -> dict[str, tuple[str, float]]:
    if isinstance(raw, Mapping):
        out = {}
        for gene, item in raw.items():
            if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) < 2:
                raise ConfigurationError(
                    f"gene_positions[{gene!r}] must be [chromosome, start].",
                    code="INVALID_OPTION",
                    hint="e.g. {'TP53': ['chr17', 7661779]}.",
                )
            out[str(gene)] = (str(item[0]), float(item[1]))
        return out
    raise ConfigurationError(
        "options['gene_positions'] must map gene -> [chromosome, start].",
        code="MISSING_OPTION",
        hint="Genomic order is needed to smooth expression along chromosomes.",
    )


def op_infer_cnv(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    positions = _gene_positions(options.get("gene_positions"))
    labels = option_labels(cm, options, "groups", "a label per cell")
    reference = opt_str(options, "reference_group", None)
    if reference is None:
        raise ConfigurationError(
            "options['reference_group'] (label of normal reference cells) is required.",
            code="MISSING_OPTION",
            hint="Reference cells define the diploid baseline.",
        )
    ref_idx = [j for j, group in enumerate(labels) if group == reference]
    if not ref_idx:
        raise InputFormatError(
            f"No cell carries the reference label {reference!r}.",
            code="UNKNOWN_GROUP",
            hint=f"Labels present: {sorted(set(labels))}.",
        )
    window = opt_int(options, "window", 101, 1, 10_001)
    min_mean = opt_float(options, "min_mean", 0.1, 0.0)
    target_sum = opt_float(options, "target_sum", 1e6, 0.0, math.inf, exclusive_min=True)
    log2 = scale_columns(cm, target_sum)
    log2 = CountMatrix(
        log2.genes, log2.cells, [[math.log2(x + 1.0) for x in row] for row in log2.values]
    )
    kept = [
        i for i, g in enumerate(cm.genes) if g in positions and stats.mean(cm.values[i]) >= min_mean
    ]
    if not kept:
        raise InputFormatError(
            "No gene has a genomic position and mean count >= min_mean.",
            code="NO_GENES",
            hint="Check gene names or lower min_mean.",
        )
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for i in kept:
        by_chrom[positions[cm.genes[i]][0]].append(i)
    ordered_genes: list[str] = []
    chrom_of: list[str] = []
    rows: list[list[float]] = []
    for chrom in sorted(by_chrom, key=_chrom_key):
        members = sorted(by_chrom[chrom], key=lambda i: (positions[cm.genes[i]][1], cm.genes[i]))
        relative = []
        for i in members:
            baseline = stats.mean([log2.values[i][j] for j in ref_idx])
            relative.append([x - baseline for x in log2.values[i]])
        half = window // 2
        for position, i in enumerate(members):
            low, high = max(0, position - half), min(len(members), position + half + 1)
            rows.append(
                [stats.mean([relative[k][j] for k in range(low, high)]) for j in range(cm.n_cells)]
            )
            ordered_genes.append(cm.genes[i])
            chrom_of.append(chrom)
    # centre every cell on the median of its smoothed profile
    medians = [stats.median([row[j] for row in rows]) for j in range(cm.n_cells)]
    rows = [[row[j] - medians[j] for j in range(cm.n_cells)] for row in rows]
    out = CountMatrix(ordered_genes, list(cm.cells), rows)
    scores = {cell: stats.mean([abs(row[j]) for row in rows]) for j, cell in enumerate(cm.cells)}
    extra = {
        "definition": "inferCNV-style expression profile: log2(CPM+1) minus the mean of the reference cells per gene, moving average over `window` genes along each chromosome (window shrinks at chromosome ends), then median-centred per cell. No HMM / CNV calling.",
        "window": window,
        "reference_group": reference,
        "n_reference_cells": len(ref_idx),
        "chromosome_of_gene": dict(zip(ordered_genes, chrom_of, strict=False)),
        "cnv_score": scores,
        "cnv_score_definition": "mean absolute value of the modified expression of the cell",
    }
    return matrix_result(out, options, extra), []


# ---------------------------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------------------------


def _points(
    value: object, options: Mapping[str, object]
) -> tuple[list[str], list[float], list[float]]:
    """Point names and their first two coordinates; rows are points unless orientation says otherwise."""

    opts = {k: v for k, v in options.items() if k not in {"genes", "cells"}}
    opts["orientation"] = options.get("orientation", "cells_by_genes")
    table = parse_matrix(value, opts, allow_negative=True)
    # parse_matrix always returns dimensions as rows ("genes") and points as columns ("cells").
    if table.n_genes < 2:
        raise InputFormatError(
            "An embedding needs at least two coordinates (x and y) per point.",
            code="EMBEDDING_SHAPE",
            hint="Pass [x, y] per point.",
        )
    return table.cells, table.values[0], table.values[1]


def op_plot_embedding(
    value: object, options: dict[str, object], *, spatial: bool = False
) -> Result:
    names, xs, ys = _points(value, options)
    title = opt_str(options, "title", "Spatial map" if spatial else "Embedding") or ""
    labels = (
        labels_for(names, options["labels"], "options['labels']")
        if options.get("labels") is not None
        else None
    )
    values = None
    if options.get("values") is not None:
        raw = options["values"]
        if isinstance(raw, Mapping):
            values = [float(raw[n]) for n in names]
        else:
            values = [float(v) for v in raw]  # type: ignore[union-attr]
            if len(values) != len(names):
                raise InputFormatError(
                    "options['values'] must have one number per point.",
                    code="LABEL_COUNT_MISMATCH",
                    hint="Match the point order.",
                )
    drawing = svg.scatter_svg(
        xs, ys, names=names, labels=labels, values=values, title=title, invert_y=spatial
    )
    path = svg.write_svg(drawing, options)
    return {"format": "svg", "svg": drawing, "n_points": len(names), "output_path": path}, []


def op_plot_heatmap(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options, allow_negative=True)
    max_genes = opt_int(options, "max_genes", 50, 1, 500)
    max_cells = opt_int(options, "max_cells", 200, 1, 2000)
    log = opt_bool(options, "log1p", True)
    zscore = opt_bool(options, "row_zscore", False)
    order = sorted(
        range(cm.n_genes),
        key=lambda i: (-stats.variance(cm.values[i]) if cm.n_cells > 1 else 0.0, i),
    )[:max_genes]
    cols = list(range(min(cm.n_cells, max_cells)))
    data = [[cm.values[i][j] for j in cols] for i in order]
    if log:
        if any(x < 0 for row in data for x in row):
            raise InputFormatError(
                "log1p needs non-negative values.", code="NEGATIVE_VALUES", hint="Set log1p=False."
            )
        data = [[math.log1p(x) for x in row] for row in data]
    if zscore:
        scaled = []
        for row in data:
            centre, sd = stats.mean(row), math.sqrt(stats.variance(row, ddof=0))
            scaled.append([(x - centre) / sd if sd > 0 else 0.0 for x in row])
        data = scaled
    names = [cm.genes[i] for i in order]
    drawing = svg.heatmap_svg(
        data,
        names,
        [cm.cells[j] for j in cols],
        title=opt_str(options, "title", "Expression heatmap") or "",
    )
    path = svg.write_svg(drawing, options)
    return {
        "format": "svg",
        "svg": drawing,
        "genes_shown": names,
        "n_cells_shown": len(cols),
        "selection": "genes with the highest variance across the shown cells",
        "output_path": path,
    }, []


def op_plot_spatial_gene(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    gene = str(require_option(options, "gene", "gene symbol to colour spots by"))
    if gene not in cm.gene_index():
        raise InputFormatError(
            f"Gene {gene!r} is not in the matrix.", code="UNKNOWN_GENE", hint="Check the symbol."
        )
    coords = parse_coordinates(options.get("coordinates"), cm.cells, options)
    log = log1p_columns(
        cm, opt_float(options, "target_sum", DEFAULT_SC_TARGET, 0.0, math.inf, exclusive_min=True)
    )
    row = log.values[cm.gene_index()[gene]]
    drawing = svg.scatter_svg(
        [c[0] for c in coords],
        [c[1] for c in coords],
        names=cm.cells,
        values=row,
        title=opt_str(options, "title", f"{gene} (ln(1+CP10k))") or "",
        invert_y=True,
    )
    path = svg.write_svg(drawing, options)
    return {
        "format": "svg",
        "svg": drawing,
        "gene": gene,
        "n_spots": cm.n_cells,
        "min_value": min(row),
        "max_value": max(row),
        "output_path": path,
    }, []
