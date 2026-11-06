"""Bulk RNA-seq operations: normalisation, QC, differential expression/splicing, translation efficiency, annotation, tool runners.

Count matrices are genes x samples (see ``dnakit.singlecell_spatial._impl_matrix`` for the accepted forms).
"""

from __future__ import annotations

import gzip
import math
import os
import re
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
    InputFormatError,
)
from dnakit.singlecell_spatial import _impl_stats as stats
from dnakit.singlecell_spatial import _impl_svg as svg
from dnakit.singlecell_spatial._impl_de import (
    de_settings,
    differential_expression,
    rank_genes,
    scale_columns,
    two_groups,
)
from dnakit.singlecell_spatial._impl_matrix import (
    CountMatrix,
    matrix_result,
    opt_bool,
    opt_choice,
    opt_float,
    opt_int,
    opt_str,
    option_labels,
    parse_matrix,
    require_option,
    text_of,
)

Result = tuple[object, list[str]]
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
            "options['mito_prefix'] must be a string or list of strings.",
            code="INVALID_OPTION",
            hint="e.g. 'MT-'.",
        )
    return tuple(raw)


def _gene_lengths(cm: CountMatrix, raw: object) -> list[float]:
    if not isinstance(raw, Mapping):
        raise ConfigurationError(
            "options['gene_lengths'] must map gene -> length in bp.",
            code="MISSING_OPTION",
            hint="Effective or transcript lengths are required for TPM/FPKM.",
        )
    lengths = {str(k): v for k, v in raw.items()}
    missing = [g for g in cm.genes if g not in lengths]
    if missing:
        raise InputFormatError(
            f"gene_lengths lacks {len(missing)} genes (first: {missing[0]!r}).",
            code="MISSING_LENGTHS",
            hint="Provide a length for every gene.",
        )
    out = []
    for gene in cm.genes:
        value = lengths[gene]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise InputFormatError(
                f"Length of {gene!r} must be a positive number.",
                code="INVALID_LENGTH",
                hint="Lengths are in base pairs.",
            )
        out.append(float(value))
    return out


# ---------------------------------------------------------------------------------------------
# normalisation, QC, detection
# ---------------------------------------------------------------------------------------------


def op_normalize(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    method = opt_choice(
        options,
        "method",
        "cpm",
        ("cpm", "tpm", "fpkm", "rpkm", "median_of_ratios", "upper_quartile"),
    )
    totals = cm.column_totals()
    extra: dict[str, object] = {"method": method}
    warnings: list[str] = []
    if method in {"cpm", "fpkm", "rpkm"} and any(t == 0 for t in totals):
        raise InputFormatError(
            "A sample has zero total counts.",
            code="EMPTY_SAMPLE",
            hint="Remove empty samples first.",
        )
    if method == "cpm":
        values = [[x / t * 1e6 for x, t in zip(row, totals, strict=False)] for row in cm.values]
        definition = "CPM = count / library size * 1e6"
    elif method in {"fpkm", "rpkm"}:
        lengths = _gene_lengths(cm, options.get("gene_lengths"))
        values = [
            [x / (length / 1e3) / (t / 1e6) for x, t in zip(row, totals, strict=False)]
            for row, length in zip(cm.values, lengths, strict=False)
        ]
        definition = "FPKM = count / (gene length in kb * library size in millions)"
    elif method == "tpm":
        lengths = _gene_lengths(cm, options.get("gene_lengths"))
        rates = [
            [x / (length / 1e3) for x in row]
            for row, length in zip(cm.values, lengths, strict=False)
        ]
        sums = [math.fsum(col) for col in zip(*rates, strict=False)]
        if any(s == 0 for s in sums):
            raise InputFormatError(
                "A sample has zero total counts.",
                code="EMPTY_SAMPLE",
                hint="Remove empty samples first.",
            )
        values = [[x / s * 1e6 for x, s in zip(row, sums, strict=False)] for row in rates]
        definition = "TPM = (count / gene length in kb) / sum over genes of (count / gene length in kb) * 1e6"
    elif method == "median_of_ratios":
        usable = [i for i, row in enumerate(cm.values) if all(x > 0 for x in row)]
        if not usable:
            raise InputFormatError(
                "No gene has a positive count in every sample.",
                code="NO_COMMON_GENES",
                hint="Median-of-ratios needs genes expressed in all samples.",
            )
        logs = {i: stats.mean([math.log(x) for x in cm.values[i]]) for i in usable}
        factors = [
            math.exp(stats.median([math.log(cm.values[i][j]) - logs[i] for i in usable]))
            for j in range(cm.n_cells)
        ]
        values = [[x / f for x, f in zip(row, factors, strict=False)] for row in cm.values]
        extra["size_factors"] = dict(zip(cm.cells, factors, strict=False))
        definition = "DESeq2 median-of-ratios: size factor = median over genes with all counts > 0 of count / geometric mean across samples; value = count / size factor"
    else:
        q75 = []
        for j in range(cm.n_cells):
            positive = [row[j] for row in cm.values if row[j] > 0]
            if not positive:
                raise InputFormatError(
                    "A sample has no positive counts.",
                    code="EMPTY_SAMPLE",
                    hint="Remove empty samples first.",
                )
            q75.append(stats.quantile(positive, 0.75))
        reference = stats.mean(q75)
        factors = [q / reference for q in q75]
        values = [[x / f for x, f in zip(row, factors, strict=False)] for row in cm.values]
        extra["size_factors"] = dict(zip(cm.cells, factors, strict=False))
        definition = "upper quartile: size factor = 75th percentile of the sample's non-zero counts / mean of those percentiles over samples; value = count / size factor"
    out = CountMatrix(list(cm.genes), list(cm.cells), values)
    if opt_bool(options, "log2", False):
        pseudo = opt_float(options, "pseudocount", 1.0, 0.0, math.inf, exclusive_min=True)
        out = CountMatrix(
            out.genes, out.cells, [[math.log2(x + pseudo) for x in row] for row in out.values]
        )
        definition += f"; then log2(value + {pseudo})"
        extra["pseudocount"] = pseudo
    extra["definition"] = definition
    return matrix_result(out, options, extra), warnings


def op_qc(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    min_count = opt_float(options, "min_count", 1.0, 0.0)
    min_library = opt_float(options, "min_library_size", None, 0.0)
    min_genes = opt_int(options, "min_genes", None, 0)
    max_mito = opt_float(options, "max_pct_mito", None, 0.0, 100.0)
    prefixes = _mito_prefixes(options)
    totals = cm.column_totals()
    mito = [i for i, g in enumerate(cm.genes) if g.startswith(prefixes)]
    samples: dict[str, object] = {}
    failing: list[str] = []
    for j, sample in enumerate(cm.cells):
        column = [row[j] for row in cm.values]
        detected = sum(1 for x in column if x >= min_count)
        ordered = sorted(column, reverse=True)
        total = totals[j]
        top1 = ordered[0] / total if total else 0.0
        top10 = math.fsum(ordered[:10]) / total if total else 0.0
        pct_mito = 100.0 * math.fsum(column[i] for i in mito) / total if total else 0.0
        reasons = []
        if min_library is not None and total < min_library:
            reasons.append("min_library_size")
        if min_genes is not None and detected < min_genes:
            reasons.append("min_genes")
        if max_mito is not None and pct_mito > max_mito:
            reasons.append("max_pct_mito")
        if reasons:
            failing.append(sample)
        samples[sample] = {
            "library_size": total,
            "genes_detected": detected,
            "fraction_zero_genes": 1.0 - sum(1 for x in column if x > 0) / cm.n_genes,
            "top_gene_fraction": top1,
            "top10_genes_fraction": top10,
            "pct_mito": pct_mito,
            "fails": reasons,
        }
    return {
        "n_samples": cm.n_cells,
        "n_genes": cm.n_genes,
        "min_count": min_count,
        "n_mito_genes": len(mito),
        "median_library_size": stats.median(totals),
        "samples": samples,
        "failing_samples": failing,
        "thresholds": {
            "min_library_size": min_library,
            "min_genes": min_genes,
            "max_pct_mito": max_mito,
        },
    }, []


def op_detect(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    min_count = opt_float(options, "min_count", 1.0, 0.0)
    min_samples = opt_int(options, "min_samples", 1, 1, cm.n_cells)
    min_cpm = opt_float(options, "min_cpm", None, 0.0)
    totals = cm.column_totals()
    if min_cpm is not None:
        if any(t == 0 for t in totals):
            raise InputFormatError(
                "A sample has zero total counts.",
                code="EMPTY_SAMPLE",
                hint="Remove empty samples first.",
            )
        flags = [
            [x / t * 1e6 >= min_cpm for x, t in zip(row, totals, strict=False)] for row in cm.values
        ]
        rule = f"CPM >= {min_cpm}"
    else:
        flags = [[x >= min_count for x in row] for row in cm.values]
        rule = f"count >= {min_count}"
    detected = [g for g, row in zip(cm.genes, flags, strict=False) if sum(row) >= min_samples]
    return {
        "rule": f"a gene is detected in a sample when {rule}; it is expressed overall when detected in >= {min_samples} samples",
        "n_genes": cm.n_genes,
        "n_detected": len(detected),
        "detected_genes": detected[:5000],
        "detected_truncated": len(detected) > 5000,
        "per_sample_detected": {
            s: sum(1 for row in flags if row[j]) for j, s in enumerate(cm.cells)
        },
    }, []


# ---------------------------------------------------------------------------------------------
# differential expression, splicing, translation efficiency, comparison
# ---------------------------------------------------------------------------------------------


def op_differential_expression(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    labels = option_labels(cm, options, "groups", "a condition label per sample")
    group1, group2, first, second = two_groups(cm, labels, options)
    settings = de_settings(options, 1e6)
    rows = differential_expression(
        cm,
        first,
        second,
        method=settings["method"],
        target_sum=settings["target_sum"],
        pseudocount=settings["pseudocount"],
    )  # type: ignore[arg-type]
    return {
        "group1": group1,
        "group2": group2,
        "n_group1": len(first),
        "n_group2": len(second),
        "method": settings["method"],
        "normalization": settings["normalization"],
        "target_sum": settings["target_sum"],
        "pseudocount": settings["pseudocount"],
        "n_genes_tested": len(rows),
        "definition": "value = count / library size * target_sum; test on log2(value + pseudocount); log2fc = log2(mean(group1) + pc) - log2(mean(group2) + pc); q_value = Benjamini-Hochberg. Not a negative-binomial (DESeq2/edgeR) model.",
        "genes": rank_genes(rows, top_n=settings["top_n"], min_pct=settings["min_pct"]),  # type: ignore[arg-type]
    }, []


def op_differential_splicing(value: object, options: dict[str, object]) -> Result:
    inclusion = parse_matrix(value, options)
    if options.get("exclusion") is None:
        raise ConfigurationError(
            "options['exclusion'] (exclusion-junction count matrix) is required.",
            code="MISSING_OPTION",
            hint="Same events and samples as the inclusion input.",
        )
    exclusion = parse_matrix(
        options["exclusion"],
        {k: v for k, v in options.items() if k in {"orientation", "delimiter", "genes", "cells"}},
    )
    if exclusion.genes != inclusion.genes or exclusion.cells != inclusion.cells:
        raise InputFormatError(
            "Inclusion and exclusion matrices must have identical events and samples in the same order.",
            code="LAYER_MISMATCH",
            hint="Reindex the matrices.",
        )
    for matrix in (inclusion, exclusion):
        if any(x != math.floor(x) for row in matrix.values for x in row):
            raise InputFormatError(
                "Junction counts must be integers.",
                code="NON_INTEGER_COUNTS",
                hint="Use raw junction read counts.",
            )
    labels = option_labels(inclusion, options, "groups", "a condition label per sample")
    group1, group2, first, second = two_groups(inclusion, labels, options)
    len_inc = opt_float(options, "inclusion_length", 1.0, 0.0, math.inf, exclusive_min=True)
    len_exc = opt_float(options, "exclusion_length", 1.0, 0.0, math.inf, exclusive_min=True)
    min_reads = opt_int(options, "min_reads", 10, 0)

    def psi(inc: float, exc: float) -> float | None:
        a, b = inc / len_inc, exc / len_exc
        return a / (a + b) if a + b > 0 else None

    rows = []
    for event, inc_row, exc_row in zip(
        inclusion.genes, inclusion.values, exclusion.values, strict=False
    ):
        record: dict[str, object] = {"event": event}
        means = []
        pooled = []
        for _name, members in ((group1, first), (group2, second)):
            values = [p for p in (psi(inc_row[j], exc_row[j]) for j in members) if p is not None]
            inc_sum, exc_sum = (
                sum(int(inc_row[j]) for j in members),
                sum(int(exc_row[j]) for j in members),
            )
            pooled.append((inc_sum, exc_sum))
            means.append(stats.mean(values) if values else None)
        record.update(
            psi_group1=means[0],
            psi_group2=means[1],
            inclusion_group1=pooled[0][0],
            exclusion_group1=pooled[0][1],
            inclusion_group2=pooled[1][0],
            exclusion_group2=pooled[1][1],
        )
        enough = (
            all(sum(p) >= min_reads for p in pooled)
            and means[0] is not None
            and means[1] is not None
        )
        if enough:
            record["delta_psi"] = means[0] - means[1]  # type: ignore[operator]
            record["p_value"] = stats.fisher_exact_two_sided(
                pooled[0][0], pooled[0][1], pooled[1][0], pooled[1][1]
            )
        else:
            record["delta_psi"] = None
            record["p_value"] = None
        rows.append(record)
    tested = [r for r in rows if r["p_value"] is not None]
    for record, q in zip(
        tested, stats.benjamini_hochberg([r["p_value"] for r in tested]), strict=False
    ):  # type: ignore[misc]
        record["q_value"] = q
    for record in rows:
        record.setdefault("q_value", None)
    rows.sort(
        key=lambda r: (
            r["p_value"] is None,
            r["p_value"] if r["p_value"] is not None else 0.0,
            r["event"],
        )
    )  # type: ignore[operator]
    top_n = opt_int(options, "top_n", 100, 1, 10**7)
    return {
        "group1": group1,
        "group2": group2,
        "n_events": len(rows),
        "n_tested": len(tested),
        "definition": "PSI per sample = (inc/inclusion_length) / (inc/inclusion_length + exc/exclusion_length); psi_group = mean over samples with reads; delta_psi = psi_group1 - psi_group2; "
        "p_value = two-sided Fisher exact test on the pooled inclusion/exclusion counts of each group (ignores between-replicate variance, so it is anti-conservative); q_value = Benjamini-Hochberg",
        "min_reads": min_reads,
        "events": rows[:top_n],
    }, []


def op_translation_efficiency(value: object, options: dict[str, object]) -> Result:
    ribo = parse_matrix(value, options)
    if options.get("rna") is None:
        raise ConfigurationError(
            "options['rna'] (RNA-seq count matrix of the same genes) is required.",
            code="MISSING_OPTION",
            hint="Ribo-seq is the input; matched RNA-seq goes in options['rna'].",
        )
    rna = parse_matrix(
        options["rna"], {k: v for k, v in options.items() if k in {"orientation", "delimiter"}}
    )
    rna_index = rna.gene_index()
    genes = [g for g in ribo.genes if g in rna_index]
    if not genes:
        raise InputFormatError(
            "Ribo-seq and RNA-seq share no genes.",
            code="NO_SHARED_GENES",
            hint="Use the same gene identifiers.",
        )
    if set(ribo.cells) == set(rna.cells):
        pairs = [(ribo.cells.index(s), rna.cells.index(s), s) for s in ribo.cells]
    elif ribo.n_cells == rna.n_cells:
        pairs = [(j, j, ribo.cells[j]) for j in range(ribo.n_cells)]
    else:
        raise InputFormatError(
            "Ribo-seq and RNA-seq need the same sample names or the same number of samples.",
            code="SAMPLE_MISMATCH",
            hint="Name matching samples identically.",
        )
    ri = ribo.gene_index()
    lengths = (
        _gene_lengths(CountMatrix(genes, [], []), options["gene_lengths"])
        if options.get("gene_lengths") is not None
        else None
    )
    r_tot, n_tot = ribo.column_totals(), rna.column_totals()
    use_log2 = opt_bool(options, "log2", False)
    pseudo = opt_float(options, "pseudocount", 0.0, 0.0)
    per_gene: dict[str, dict[str, float | None]] = {g: {} for g in genes}
    for rj, nj, name in pairs:
        if r_tot[rj] == 0 or n_tot[nj] == 0:
            raise InputFormatError(
                f"Sample {name!r} has zero total counts.",
                code="EMPTY_SAMPLE",
                hint="Remove empty samples first.",
            )
        for k, gene in enumerate(genes):
            scale = (lengths[k] / 1e3) if lengths else 1.0
            ribo_abundance = ribo.values[ri[gene]][rj] / r_tot[rj] * 1e6 / scale
            rna_abundance = rna.values[rna_index[gene]][nj] / n_tot[nj] * 1e6 / scale
            if rna_abundance + pseudo <= 0 or (use_log2 and ribo_abundance + pseudo <= 0):
                per_gene[gene][name] = None
                continue
            te = (ribo_abundance + pseudo) / (rna_abundance + pseudo)
            per_gene[gene][name] = math.log2(te) if use_log2 else te
    summary = {}
    for gene, values in per_gene.items():
        valid = [v for v in values.values() if v is not None]
        summary[gene] = {"te": values, "mean_te": stats.mean(valid) if valid else None}
    return {
        "definition": (
            "TE = ribo abundance / RNA abundance, each as "
            + ("RPKM (counts per kb per million)" if lengths else "CPM (counts per million)")
            + " of its own library"
            + ("; reported as log2(TE)" if use_log2 else "")
            + "; null when the RNA abundance is zero"
        ),
        "n_genes": len(genes),
        "samples": [p[2] for p in pairs],
        "genes": summary,
    }, []


def op_compare(value: object, options: dict[str, object]) -> Result:
    left = parse_matrix(value, options)
    method = opt_choice(options, "method", "pearson", ("pearson", "spearman"))
    use_log2 = opt_bool(options, "log2", True)
    normalization = opt_choice(options, "normalization", "none", ("none", "cpm"))
    shared_opts = {k: v for k, v in options.items() if k in {"orientation", "delimiter"}}
    right = (
        parse_matrix(options["other"], shared_opts) if options.get("other") is not None else left
    )

    def prepare(cm: CountMatrix, genes: Sequence[str]) -> dict[str, list[float]]:
        scaled = scale_columns(cm, 1e6) if normalization == "cpm" else cm
        index = scaled.gene_index()
        columns = {}
        for j, name in enumerate(scaled.cells):
            column = [scaled.values[index[g]][j] for g in genes]
            columns[name] = [math.log2(x + 1.0) for x in column] if use_log2 else column
        return columns

    right_index = right.gene_index()
    genes = [g for g in left.genes if g in right_index]
    if len(genes) < 3:
        raise InputFormatError(
            f"Only {len(genes)} genes are shared; at least 3 are needed.",
            code="TOO_FEW_SHARED_GENES",
            hint="Use the same gene identifiers.",
        )
    a, b = prepare(left, genes), prepare(right, genes)
    func = stats.pearson if method == "pearson" else stats.spearman
    table = {x: {y: func(a[x], b[y]) for y in b} for x in a}
    return {
        "definition": f"{method} correlation across {len(genes)} shared genes of "
        + ("log2(value + 1)" if use_log2 else "the values")
        + (" after CPM scaling" if normalization == "cpm" else ""),
        "method": method,
        "n_shared_genes": len(genes),
        "mode": "input versus other"
        if options.get("other") is not None
        else "all pairs within the input",
        "correlation": table,
    }, []


def op_plot(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    kind = opt_choice(options, "kind", "heatmap", ("heatmap", "scatter"))
    use_log2 = opt_bool(options, "log2", True)
    title = opt_str(options, "title", None)
    if kind == "heatmap":
        max_genes = opt_int(options, "max_genes", 50, 1, 500)
        order = sorted(
            range(cm.n_genes),
            key=lambda i: (-stats.variance(cm.values[i]) if cm.n_cells > 1 else 0.0, i),
        )[:max_genes]
        data = [[math.log2(x + 1.0) if use_log2 else x for x in cm.values[i]] for i in order]
        names = [cm.genes[i] for i in order]
        drawing = svg.heatmap_svg(data, names, cm.cells, title=title or "Expression heatmap")
        extra: dict[str, object] = {"genes_shown": names}
    else:
        x_name, y_name = (
            str(require_option(options, "x", "sample on the x axis")),
            str(require_option(options, "y", "sample on the y axis")),
        )
        index = cm.cell_index()
        for name in (x_name, y_name):
            if name not in index:
                raise InputFormatError(
                    f"Sample {name!r} is not in the matrix.",
                    code="UNKNOWN_SAMPLE",
                    hint=f"Samples: {cm.cells[:10]}.",
                )
        xs = [row[index[x_name]] for row in cm.values]
        ys = [row[index[y_name]] for row in cm.values]
        if use_log2:
            xs, ys = [math.log2(v + 1.0) for v in xs], [math.log2(v + 1.0) for v in ys]
        drawing = svg.scatter_svg(
            xs,
            ys,
            names=cm.genes,
            title=title or f"{y_name} vs {x_name}" + (" (log2(x+1))" if use_log2 else ""),
        )
        extra = {"x": x_name, "y": y_name, "pearson_r": stats.pearson(xs, ys)}
    path = svg.write_svg(drawing, options)
    return {"format": "svg", "svg": drawing, "kind": kind, "output_path": path, **extra}, []


# ---------------------------------------------------------------------------------------------
# GTF annotation
# ---------------------------------------------------------------------------------------------

_ATTRIBUTE = re.compile(r'\s*([A-Za-z_][\w.]*)\s+(?:"([^"]*)"|([^;]*))\s*;?')


def _gene_names(value: object) -> list[str]:
    if value is None:
        raise InputFormatError(
            "input is required.", code="EMPTY_INPUT", hint="Provide gene identifiers."
        )
    if isinstance(value, Mapping):
        names = [str(k) for k in value]
    elif isinstance(value, (str, bytes)) or hasattr(value, "__fspath__"):
        names = [
            line.split()[0]
            for line in text_of(value, "genes").splitlines()
            if line.strip() and not line.startswith("#")
        ]
    elif isinstance(value, Sequence):
        names = [str(v) for v in value]
    else:
        raise InputFormatError(
            "Unsupported gene input.",
            code="UNSUPPORTED_INPUT",
            hint="Use a list, text lines or a mapping keyed by gene.",
        )
    if not names:
        raise InputFormatError(
            "No gene identifiers given.", code="EMPTY_INPUT", hint="Provide at least one gene."
        )
    return names


def op_annotate(value: object, options: dict[str, object]) -> Result:
    names = _gene_names(value)
    gtf = Path(str(require_option(options, "gtf", "path to a GTF annotation file"))).expanduser()
    if not gtf.is_file():
        raise InputFormatError(
            f"GTF file not found: {gtf}",
            code="FILE_NOT_FOUND",
            hint="Pass an existing .gtf or .gtf.gz path.",
        )
    if gtf.stat().st_size > 2 * 1024**3:
        raise InputFormatError(
            "The GTF file is larger than 2 GiB.",
            code="INPUT_TOO_LARGE",
            hint="Subset the annotation.",
        )
    wanted = set(names) | {re.sub(r"\.\d+$", "", n) for n in names}
    by_key: dict[str, dict[str, object]] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    n_genes = 0
    with opener(gtf, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            n_genes += 1
            attrs = {
                m.group(1): (m.group(2) if m.group(2) is not None else (m.group(3) or "").strip())
                for m in _ATTRIBUTE.finditer(fields[8])
            }
            gene_id, gene_name = attrs.get("gene_id", ""), attrs.get("gene_name", "")
            keys = {gene_id, re.sub(r"\.\d+$", "", gene_id), gene_name} - {""}
            if not keys & wanted:
                continue
            try:
                record = {
                    "gene_id": gene_id,
                    "gene_name": gene_name or None,
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "biotype": attrs.get("gene_biotype", attrs.get("gene_type")),
                    "source": fields[1],
                }
            except ValueError:
                raise InputFormatError(
                    f"Bad coordinates in GTF line: {line[:80]!r}",
                    code="GTF_FORMAT",
                    hint="Columns 4 and 5 must be integers.",
                ) from None
            for key in keys:
                by_key.setdefault(key, record)
    if n_genes == 0:
        raise InputFormatError(
            "The GTF file has no 'gene' features.",
            code="GTF_FORMAT",
            hint="Use a GENCODE/Ensembl GTF that contains gene lines.",
        )
    annotations: dict[str, object] = {}
    missing: list[str] = []
    for name in names:
        record = by_key.get(name) or by_key.get(re.sub(r"\.\d+$", "", name))
        if record is None:
            missing.append(name)
            annotations[name] = None
        else:
            annotations[name] = record
    return {
        "definition": "lookup of each identifier (gene_id with or without version suffix, or gene_name) among the 'gene' features of the GTF; coordinates are 1-based inclusive as in GTF",
        "gtf": str(gtf),
        "n_requested": len(names),
        "n_found": len(names) - len(missing),
        "not_found": missing[:_MAX_LISTED],
        "annotations": annotations,
    }, []


# ---------------------------------------------------------------------------------------------
# external programs
# ---------------------------------------------------------------------------------------------


def _paths(value: object, what: str, *, must_exist: bool = True) -> list[str]:
    if value is None:
        raise InputFormatError("input is required.", code="EMPTY_INPUT", hint=f"Provide {what}.")
    items = (
        [value]
        if isinstance(value, (str, Path))
        else list(value)
        if isinstance(value, Sequence)
        else None
    )  # type: ignore[arg-type]
    if not items:
        raise InputFormatError(
            f"Provide {what} as a path or list of paths.",
            code="UNSUPPORTED_INPUT",
            hint="Pass file paths.",
        )
    out = []
    for item in items:
        if not isinstance(item, (str, Path)) or not str(item).strip() or "\x00" in str(item):
            raise InputFormatError(
                f"Invalid path {item!r}.", code="INVALID_PATH", hint="Use non-empty path strings."
            )
        text = os.path.expanduser(str(item))
        if text.startswith("-"):
            raise InputFormatError(
                f"Path {item!r} must not start with '-'.",
                code="INVALID_PATH",
                hint="Use an absolute path.",
            )
        if must_exist and not os.path.isfile(text):
            raise InputFormatError(
                f"File not found: {item}", code="FILE_NOT_FOUND", hint="Pass existing files."
            )
        out.append(text)
    return out


def _new_dir_or_prefix(value: object, key: str, *, directory: bool) -> str:
    text = value
    if (
        not isinstance(text, (str, Path))
        or not str(text).strip()
        or "\x00" in str(text)
        or str(text).startswith("-")
    ):
        raise ConfigurationError(
            f"options[{key!r}] is required and must be a path.",
            code="MISSING_OPTION",
            hint=f"Pass options={{{key!r}: ...}}.",
        )
    path = Path(os.path.expanduser(str(text)))
    if directory:
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise ConfigurationError(
                f"Output directory already exists and is not empty: {path}",
                code="OUTPUT_EXISTS",
                hint="Choose a new directory.",
            )
        if not path.parent.is_dir():
            raise ConfigurationError(
                f"Parent directory does not exist: {path.parent}",
                code="OUTPUT_DIR",
                hint="Create the parent directory.",
            )
    else:
        if path.exists():
            raise ConfigurationError(
                f"Output file already exists: {path}",
                code="OUTPUT_EXISTS",
                hint="Choose a new output path.",
            )
        if not path.parent.is_dir():
            raise ConfigurationError(
                f"Output directory does not exist: {path.parent}",
                code="OUTPUT_DIR",
                hint="Create the directory first.",
            )
    return str(path)


def _existing(value: object, key: str, what: str, *, directory_ok: bool = False) -> str:
    if (
        not isinstance(value, (str, Path))
        or not str(value).strip()
        or "\x00" in str(value)
        or str(value).startswith("-")
    ):
        raise ConfigurationError(
            f"options[{key!r}] is required: {what}.",
            code="MISSING_OPTION",
            hint=f"Pass options={{{key!r}: ...}}.",
        )
    path = os.path.expanduser(str(value))
    if not (os.path.isfile(path) or (directory_ok and os.path.isdir(path))):
        raise ConfigurationError(
            f"options[{key!r}] does not exist: {path}", code="FILE_NOT_FOUND", hint=what
        )
    return path


def _execute(arguments: list[str], timeout: float) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(
            arguments, capture_output=True, text=True, timeout=timeout, check=False, shell=False
        )
    except subprocess.TimeoutExpired as exc:
        raise BackendTimeoutError(
            f"{Path(arguments[0]).name} exceeded the {timeout:.0f}s timeout.",
            hint="Raise options['timeout_seconds'].",
        ) from exc
    except OSError as exc:
        raise BackendExecutionError(
            f"Could not start {arguments[0]}: {exc}", hint="Check the installation."
        ) from exc
    if completed.returncode != 0:
        raise BackendExecutionError(
            f"{Path(arguments[0]).name} exited with status {completed.returncode}: {completed.stderr.strip()[-500:]}",
            context={"command": " ".join(arguments)[:500]},
            hint="Check the inputs and the tool's log.",
        )
    return completed


def _tool_missing(
    module: str, operation: str, value: object, options: dict[str, object], names: Sequence[str]
) -> dict[str, object]:
    return candidate_not_implemented(
        module,
        operation,
        value,
        options=options,
        requires=f"{' or '.join(names)} on PATH",
        message=f"{module}.{operation} runs the external program {' / '.join(names)}, which was not found on PATH.",
    )


def _run_options(options: Mapping[str, object]) -> tuple[int, float]:
    return opt_int(options, "threads", 1, 1, 256), opt_float(
        options, "timeout_seconds", 3600.0, 0.0, 86_400.0, exclusive_min=True
    )


def _read_table(path: Path, comment: str | None = None) -> list[dict[str, str]]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not (comment and line.startswith(comment))
    ]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=False)) for line in lines[1:]]


def _tail_listing(rows: Sequence[Mapping[str, object]]) -> list[Mapping[str, object]]:
    return list(rows[:_MAX_LISTED])


def build_featurecounts(
    input: object, options: dict[str, object]
) -> tuple[list[str], dict[str, object]]:
    bams = _paths(input, "alignment (BAM/SAM) files")
    annotation = _existing(options.get("annotation"), "annotation", "path to a GTF/GFF annotation")
    output = _new_dir_or_prefix(options.get("output"), "output", directory=False)
    threads, timeout = _run_options(options)
    strand = opt_int(options, "strand", 0, 0, 2)
    feature_type = opt_str(options, "feature_type", "exon")
    attribute = opt_str(options, "attribute", "gene_id")
    args = [
        "-a",
        annotation,
        "-o",
        output,
        "-T",
        str(threads),
        "-s",
        str(strand),
        "-t",
        feature_type,
        "-g",
        attribute,
    ]  # type: ignore[list-item]
    if opt_bool(options, "paired", False):
        args += ["-p", "--countReadPairs"]
    return args + bams, {"output": output, "timeout": timeout}


def op_run_featurecounts(value: object, options: dict[str, object]) -> Result | dict[str, object]:
    args, info = build_featurecounts(value, options)
    exe = shutil.which("featureCounts")
    if exe is None:
        return _tool_missing(
            "transcriptomics", "run_featurecounts", value, options, ["featureCounts (subread)"]
        )  # type: ignore[return-value]
    completed = _execute([exe, *args], info["timeout"])  # type: ignore[arg-type]
    output = Path(info["output"])  # type: ignore[arg-type]
    if not output.is_file():
        raise BackendExecutionError(
            "featureCounts finished but wrote no counts file.", hint="Inspect the tool's stderr."
        )
    rows = _read_table(output, "#")
    samples = list(rows[0])[6:] if rows else []
    totals = {s: sum(float(r[s]) for r in rows) for s in samples}
    return {
        "program": exe,
        "command": [exe, *args],
        "counts_file": str(output),
        "summary_file": str(output) + ".summary"
        if Path(str(output) + ".summary").exists()
        else None,
        "n_genes": len(rows),
        "samples": samples,
        "assigned_counts": totals,
        "stderr_tail": completed.stderr.strip()[-500:],
    }, []


def op_run_salmon(value: object, options: dict[str, object]) -> Result | dict[str, object]:
    reads = _paths(value, "FASTQ file(s) (one for single-end, two for paired-end)")
    if len(reads) > 2:
        raise InputFormatError(
            "Provide one (single-end) or two (paired-end) FASTQ files.",
            code="UNSUPPORTED_INPUT",
            hint="Run salmon once per sample.",
        )
    index = _existing(
        options.get("index"), "index", "path to a salmon index directory", directory_ok=True
    )
    out_dir = _new_dir_or_prefix(options.get("output_dir"), "output_dir", directory=True)
    threads, timeout = _run_options(options)
    library = opt_str(options, "library_type", "A")
    args = [
        "quant",
        "-i",
        index,
        "-l",
        library,
        "-o",
        out_dir,
        "-p",
        str(threads),
        "--validateMappings",
    ]  # type: ignore[list-item]
    args += ["-r", reads[0]] if len(reads) == 1 else ["-1", reads[0], "-2", reads[1]]
    exe = shutil.which("salmon")
    if exe is None:
        return _tool_missing(
            "transcriptomics", "run_salmon_quantification", value, options, ["salmon"]
        )  # type: ignore[return-value]
    completed = _execute([exe, *args], timeout)
    quant = Path(out_dir) / "quant.sf"
    if not quant.is_file():
        raise BackendExecutionError(
            "salmon finished but wrote no quant.sf.", hint="Inspect the tool's stderr."
        )
    rows = _read_table(quant)
    return {
        "program": exe,
        "command": [exe, *args],
        "quant_file": str(quant),
        "n_transcripts": len(rows),
        "total_reads_assigned": math.fsum(float(r["NumReads"]) for r in rows),
        "total_tpm": math.fsum(float(r["TPM"]) for r in rows),
        "top_transcripts": _tail_listing(sorted(rows, key=lambda r: -float(r["TPM"]))[:20]),
        "stderr_tail": completed.stderr.strip()[-500:],
    }, []


def op_run_kallisto(value: object, options: dict[str, object]) -> Result | dict[str, object]:
    reads = _paths(value, "FASTQ file(s) (one for single-end, two for paired-end)")
    if len(reads) > 2:
        raise InputFormatError(
            "Provide one (single-end) or two (paired-end) FASTQ files.",
            code="UNSUPPORTED_INPUT",
            hint="Run kallisto once per sample.",
        )
    index = _existing(options.get("index"), "index", "path to a kallisto index file")
    out_dir = _new_dir_or_prefix(options.get("output_dir"), "output_dir", directory=True)
    threads, timeout = _run_options(options)
    args = ["quant", "-i", index, "-o", out_dir, "-t", str(threads)]
    if len(reads) == 1:
        length = opt_float(options, "fragment_length", None, 0.0, exclusive_min=True)
        sd = opt_float(options, "fragment_sd", None, 0.0, exclusive_min=True)
        if length is None or sd is None:
            raise ConfigurationError(
                "Single-end data needs options['fragment_length'] and options['fragment_sd'].",
                code="MISSING_OPTION",
                hint="Typical values: 200 and 20.",
            )
        args += ["--single", "-l", str(length), "-s", str(sd)]
    args += reads
    exe = shutil.which("kallisto")
    if exe is None:
        return _tool_missing(
            "transcriptomics", "run_kallisto_quantification", value, options, ["kallisto"]
        )  # type: ignore[return-value]
    completed = _execute([exe, *args], timeout)
    table = Path(out_dir) / "abundance.tsv"
    if not table.is_file():
        raise BackendExecutionError(
            "kallisto finished but wrote no abundance.tsv.", hint="Inspect the tool's stderr."
        )
    rows = _read_table(table)
    return {
        "program": exe,
        "command": [exe, *args],
        "abundance_file": str(table),
        "n_transcripts": len(rows),
        "total_est_counts": math.fsum(float(r["est_counts"]) for r in rows),
        "total_tpm": math.fsum(float(r["tpm"]) for r in rows),
        "top_transcripts": _tail_listing(sorted(rows, key=lambda r: -float(r["tpm"]))[:20]),
        "stderr_tail": completed.stderr.strip()[-500:],
    }, []


def op_align(value: object, options: dict[str, object]) -> Result | dict[str, object]:
    reads = _paths(value, "FASTQ file(s) (one for single-end, two for paired-end)")
    if len(reads) > 2:
        raise InputFormatError(
            "Provide one (single-end) or two (paired-end) FASTQ files.",
            code="UNSUPPORTED_INPUT",
            hint="Align one sample at a time.",
        )
    aligner = opt_choice(options, "aligner", "auto", ("auto", "star", "hisat2"))
    index = options.get("index")
    if (
        not isinstance(index, (str, Path))
        or not str(index).strip()
        or "\x00" in str(index)
        or str(index).startswith("-")
    ):
        raise ConfigurationError(
            "options['index'] is required: STAR genome directory or HISAT2 index prefix.",
            code="MISSING_OPTION",
            hint="Pass options={'index': ...}.",
        )
    index = os.path.expanduser(str(index))
    prefix = _new_dir_or_prefix(options.get("output_prefix"), "output_prefix", directory=False)
    threads, timeout = _run_options(options)
    candidates = {"star": "STAR", "hisat2": "hisat2"}
    chosen = None
    for key, exe_name in candidates.items():
        if aligner in {"auto", key} and shutil.which(exe_name):
            chosen = key
            break
    if chosen is None:
        return _tool_missing(
            "transcriptomics",
            "align",
            value,
            options,
            ["STAR", "hisat2"] if aligner == "auto" else [candidates[aligner]],
        )  # type: ignore[return-value]
    if chosen == "star" and not os.path.isdir(index):
        raise ConfigurationError(
            f"STAR genome directory not found: {index}",
            code="FILE_NOT_FOUND",
            hint="Build it with STAR --runMode genomeGenerate.",
        )
    if chosen == "hisat2" and not any(
        os.path.isfile(f"{index}.1.{ext}") for ext in ("ht2", "ht2l")
    ):
        raise ConfigurationError(
            f"HISAT2 index files {index}.1.ht2 not found.",
            code="FILE_NOT_FOUND",
            hint="Build it with hisat2-build.",
        )
    exe = shutil.which(candidates[chosen])
    assert exe is not None
    if chosen == "star":
        args = [
            "--runThreadN",
            str(threads),
            "--genomeDir",
            index,
            "--readFilesIn",
            *reads,
            "--outFileNamePrefix",
            prefix,
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
        ]
        if reads[0].endswith(".gz"):
            args += ["--readFilesCommand", "zcat"]
        output = prefix + "Aligned.sortedByCoord.out.bam"
    else:
        args = ["-p", str(threads), "-x", index]
        args += ["-U", reads[0]] if len(reads) == 1 else ["-1", reads[0], "-2", reads[1]]
        output = prefix + ".sam"
        args += ["-S", output]
    completed = _execute([exe, *args], timeout)
    if not Path(output).is_file():
        raise BackendExecutionError(
            f"{candidates[chosen]} finished but wrote no alignment file.",
            hint="Inspect the tool's stderr.",
        )
    return {
        "program": exe,
        "aligner": chosen,
        "command": [exe, *args],
        "alignment_file": output,
        "alignment_bytes": Path(output).stat().st_size,
        "stderr_tail": completed.stderr.strip()[-1000:],
    }, []


def op_quantify(value: object, options: dict[str, object]) -> Result | dict[str, object]:
    method = opt_choice(options, "method", "auto", ("auto", "salmon", "kallisto", "featurecounts"))
    tools = {"salmon": "salmon", "kallisto": "kallisto", "featurecounts": "featureCounts"}
    first = (
        value[0]
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value
        else value
    )
    alignments = isinstance(first, (str, Path)) and str(first).lower().endswith((".bam", ".sam"))
    pool = ["featurecounts"] if alignments else ["salmon", "kallisto"]
    if method != "auto" and method not in pool:
        raise ConfigurationError(
            f"method={method!r} does not fit this input ({'alignment files' if alignments else 'FASTQ files'}).",
            code="INVALID_OPTION",
            hint=f"Use one of {pool}.",
        )
    chosen = next((m for m in pool if method in {"auto", m} and shutil.which(tools[m])), None)
    if chosen is None:
        return _tool_missing(
            "transcriptomics",
            "quantify",
            value,
            options,
            [tools[m] for m in pool if method in {"auto", m}],
        )  # type: ignore[return-value]
    runner = {
        "salmon": op_run_salmon,
        "kallisto": op_run_kallisto,
        "featurecounts": op_run_featurecounts,
    }[chosen]
    result = runner(value, options)
    if isinstance(result, dict):
        return result
    payload, warnings = result
    payload["method"] = chosen  # type: ignore[index]
    return payload, warnings
