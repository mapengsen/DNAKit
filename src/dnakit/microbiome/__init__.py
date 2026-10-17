"""Local candidate adapters for the microbiome capability family."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.microbiome import _impl_local as _impl
from dnakit.microbiome.external import (
    classify_centrifuge,
    classify_gtdbtk,
    classify_kaiju,
    classify_kraken2,
    estimate_bracken_abundance,
    profile_humann,
)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse taxonomy lineage strings into ranked fields (domain..species) and report the lowest assigned rank.

    ``input``: one lineage string, a list of lineage strings, or a mapping ``{feature_id: lineage}``.
    Lineages may carry rank prefixes (``k__``/``d__``, ``p__``, ``c__``, ``o__``, ``f__``, ``g__``, ``s__``;
    ``k__`` and ``d__`` both map to ``domain``) or be unprefixed ``;``/``|`` separated levels read positionally
    as domain, phylum, class, order, family, genus, species. Empty or placeholder names (``__``, ``Unassigned``,
    ``uncultured``...) become null. ``options``: ``entity`` (only ``"taxonomy"`` is supported). This does not
    classify sequences; it only parses lineage text that is already assigned.
    """
    opts = dict(options or {})
    _impl.entity(opts, "taxonomy", ("taxonomy",))
    if isinstance(input, Mapping):
        items = [(str(k), str(v)) for k, v in input.items()]
    elif isinstance(input, str):
        items = [("lineage1", input)]
    elif isinstance(input, Sequence) and all(isinstance(v, str) for v in input):
        items = [(f"lineage{i + 1}", v) for i, v in enumerate(input)]
    else:
        raise InputFormatError(
            "annotate expects a lineage string, a list of strings or a mapping id -> lineage.",
            code="MICROBIOME_TAXONOMY_INVALID",
            hint="Example: 'k__Bacteria;p__Firmicutes;c__Bacilli'.",
        )
    if not items or len(items) > 1_000_000:
        raise InputFormatError("Provide between 1 and 1,000,000 lineages.", code="MICROBIOME_TAXONOMY_INVALID")
    result = [{"id": name, "lineage": text, **_impl.parse_taxonomy(text)} for name, text in items]
    return candidate_ok("microbiome", "annotate", input, result, options=opts)

def assemble(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assemble metagenomic reads into contigs (not implemented: needs an assembler such as MEGAHIT or metaSPAdes)."""
    return candidate_not_implemented(
        "microbiome", "assemble", input, options=options,
        requires="a metagenome assembler (MEGAHIT, metaSPAdes) on PATH",
        message="Metagenome assembly is not implemented natively; no assembler is wired into this function.",
    )

def bin(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Bin metagenomic contigs into genome bins (not implemented: needs coverage data and a binner such as MetaBAT2)."""
    return candidate_not_implemented(
        "microbiome", "bin", input, options=options,
        requires="a contig binner (MetaBAT2, MaxBin2, CONCOCT) plus read-coverage data",
        message="Contig binning is not implemented natively.",
    )

def build_cooccurrence_network(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a taxon co-occurrence network from pairwise correlations of relative abundance across samples.

    ``input``: abundance table (``{sample: {taxon: count}}``, ``{"samples","taxa","counts"}``, rows of counts with
    ``options["samples"]``/``options["taxa"]``, or TSV text/path). ``options``: ``method`` ``"spearman"`` (default)
    or ``"pearson"``; ``min_abs_correlation`` (default 0.6); ``min_samples_present`` (default 3: a taxon needs
    non-zero abundance in at least this many samples); ``orientation`` for text tables (``"taxa_rows"`` default).
    Result: ``nodes``, ``edges`` (``source``, ``target``, ``correlation``, sorted by |r| descending) and
    ``n_samples``. Correlations are computed on per-sample relative abundance; no significance test is applied, and
    relative-abundance correlations are subject to compositionality bias.
    """
    opts = dict(options or {})
    table = _impl.parse_table(input, opts)
    method = _impl.opt_str(opts, "method", "spearman", ("spearman", "pearson"))
    cutoff = _impl.opt_float(opts, "min_abs_correlation", 0.6, 0.0, 1.0)
    min_present = _impl.opt_int(opts, "min_samples_present", 3, 1, 10_000)
    if len(table.samples) < 3:
        raise InputFormatError("A co-occurrence network needs at least 3 samples.", code="MICROBIOME_TABLE_INVALID")
    _impl.require_depth(table)
    n = len(table.samples)
    columns: dict[str, list[float]] = {}
    for t, taxon in enumerate(table.taxa):
        column = [table.counts[s][t] / sum(table.counts[s]) for s in range(n)]
        if sum(1 for v in column if v > 0) >= min_present:
            columns[taxon] = _impl.average_ranks(column) if method == "spearman" else column
    names = list(columns)
    if len(names) > 2_000:
        raise ConfigurationError(
            "Too many taxa for pairwise correlation.", code="MICROBIOME_TABLE_TOO_LARGE",
            hint="Filter to at most 2000 taxa first (raise min_samples_present).",
        )
    edges = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            r = _impl.pearson(columns[names[i]], columns[names[j]])
            if r is not None and abs(r) >= cutoff:
                edges.append({"source": names[i], "target": names[j], "correlation": r})
    edges.sort(key=lambda e: (-abs(e["correlation"]), e["source"], e["target"]))
    result = {"method": method, "min_abs_correlation": cutoff, "nodes": names, "edges": edges, "n_samples": n}
    return candidate_ok(
        "microbiome", "build_cooccurrence_network", input, result, options=opts,
        warnings=["Correlations of relative abundances are compositional; interpret with care."],
    )

def build_microbial_gene_catalog(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a non-redundant microbial gene catalog (not implemented: needs gene prediction and CD-HIT/MMseqs2 clustering)."""
    return candidate_not_implemented(
        "microbiome", "build_microbial_gene_catalog", input, options=options,
        requires="a gene predictor (Prodigal) and a sequence clusterer (CD-HIT, MMseqs2) on PATH",
        message="Gene-catalog construction needs gene calling and identity clustering, which are not implemented natively.",
    )

def build_pathogen_transmission_tree(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer a pathogen transmission tree (not implemented: needs an epidemiological/phylodynamic model such as outbreaker2)."""
    return candidate_not_implemented(
        "microbiome", "build_pathogen_transmission_tree", input, options=options,
        requires="a transmission inference model (outbreaker2, TransPhylo) with sampling dates",
        message="Transmission-tree inference requires a statistical model and dated samples; none is implemented.",
    )

def calculate_ani(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate average nucleotide identity between genomes from shared k-mers (Mash formula; an estimate, not alignment ANI).

    ``input``: at least two DNA sequences as a list of strings, a mapping ``{id: sequence}``, or FASTA text/path.
    ``options``: ``k`` (default 21, 2..32); ``canonical`` (default true: k-mer and reverse complement are the
    same word). For every pair the exact k-mer sets A and B give Jaccard ``j = |A&B| / |A|B|``, containment, the
    Mash distance ``D = -ln(2j / (1 + j)) / k`` and ``ani_estimate = 1 - D`` (null when no k-mer is shared). This is
    NOT FastANI/ANIb (no fragment alignment); it is only reliable for genomes of similar size with ANI above ~0.8.
    """
    opts = dict(options or {})
    k = _impl.opt_int(opts, "k", 21, 2, 32)
    canonical = _impl.opt_bool(opts, "canonical", True)
    seqs = _impl.parse_sequences(input, minimum=2, maximum=50, max_bases=20_000_000)
    sets = {name: _impl.kmer_set(seq, k, canonical) for name, seq in seqs}
    pairs = []
    names = [name for name, _ in seqs]
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pairs.append({"a": names[i], "b": names[j], "k": k, **_impl.kmer_ani(sets[names[i]], sets[names[j]], k)})
    return candidate_ok(
        "microbiome", "calculate_ani", input, {"method": "mash_kmer_jaccard", "pairs": pairs}, options=opts,
        warnings=["k-mer set (Mash-formula) ANI estimate, not an alignment-based ANI."],
    )

def calculate_differential_abundance(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test every taxon for differential relative abundance between two sample groups (Mann-Whitney U, BH-adjusted).

    ``input``: abundance table (see ``build_cooccurrence_network``). ``options``: ``groups`` (required: a list of
    group labels in sample order, or a mapping sample -> label); ``group_a``/``group_b`` (default: the two sorted
    labels, which must be exactly two); ``pseudocount`` (default 1e-6, added to mean relative abundance for the
    log2 fold change); ``min_prevalence`` (default 0.1: taxa present in fewer samples overall are skipped).
    Per taxon: group means of relative abundance, ``log2_fold_change = log2((mean_b + pc) / (mean_a + pc))``,
    the Mann-Whitney U statistic of group_a, a two-sided p-value from ``scipy.stats.mannwhitneyu`` (scipy is an
    optional dependency; without it the function returns ``not_implemented``) and Benjamini-Hochberg ``q_value``.
    """
    opts = dict(options or {})
    table = _impl.parse_table(input, opts)
    _impl.require_depth(table)
    groups = opts.get("groups")
    if isinstance(groups, Mapping):
        try:
            labels = [str(groups[s]) for s in table.samples]
        except KeyError as exc:
            raise ConfigurationError(f"Sample {exc.args[0]!r} has no group.", code="MICROBIOME_GROUPS_INVALID") from exc
    elif isinstance(groups, Sequence) and not isinstance(groups, str) and len(groups) == len(table.samples):
        labels = [str(g) for g in groups]
    else:
        raise ConfigurationError(
            "options['groups'] must list one group label per sample (or map sample -> label).",
            code="MICROBIOME_GROUPS_INVALID", hint="Example: {'groups': ['ctrl', 'ctrl', 'case', 'case']}.",
        )
    unique = sorted(set(labels))
    group_a = str(opts.get("group_a", unique[0]))
    group_b = str(opts.get("group_b", unique[-1]))
    if group_a == group_b or group_a not in unique or group_b not in unique or ("group_a" not in opts and len(unique) != 2):
        raise ConfigurationError(
            "Exactly two distinct groups are required.", code="MICROBIOME_GROUPS_INVALID",
            hint="Provide two labels in 'groups' or choose them with group_a/group_b.",
        )
    idx_a = [i for i, g in enumerate(labels) if g == group_a]
    idx_b = [i for i, g in enumerate(labels) if g == group_b]
    if len(idx_a) < 2 or len(idx_b) < 2:
        raise ConfigurationError("Each group needs at least two samples.", code="MICROBIOME_GROUPS_INVALID")
    try:
        from scipy import stats
    except ImportError:
        return candidate_not_implemented(
            "microbiome", "calculate_differential_abundance", input, options=opts, requires="pip install scipy",
            message="The Mann-Whitney p-values use scipy.stats, which is not installed.",
        )
    pseudocount = _impl.opt_float(opts, "pseudocount", 1e-6, 0.0)
    min_prev = _impl.opt_float(opts, "min_prevalence", 0.1, 0.0, 1.0)
    rel = [[x / sum(row) for x in row] for row in table.counts]
    rows = []
    for t, taxon in enumerate(table.taxa):
        column = [rel[s][t] for s in range(len(table.samples))]
        if sum(1 for v in column if v > 0) / len(column) < min_prev:
            continue
        a = [column[i] for i in idx_a]
        b = [column[i] for i in idx_b]
        mean_a, mean_b = sum(a) / len(a), sum(b) / len(b)
        if len(set(a + b)) == 1:
            u, p = len(a) * len(b) / 2.0, 1.0
        else:
            res = stats.mannwhitneyu(a, b, alternative="two-sided")
            u, p = float(res.statistic), float(res.pvalue)
        rows.append({
            "taxon": taxon, "mean_a": mean_a, "mean_b": mean_b,
            "log2_fold_change": math.log2((mean_b + pseudocount) / (mean_a + pseudocount)),
            "u_statistic": u, "p_value": p,
        })
    for row, q in zip(rows, _impl.benjamini_hochberg([r["p_value"] for r in rows])):
        row["q_value"] = q
    rows.sort(key=lambda r: (r["p_value"], r["taxon"]))
    result = {"group_a": group_a, "group_b": group_b, "n_a": len(idx_a), "n_b": len(idx_b), "taxa": rows, "test": "mann_whitney_u"}
    return candidate_ok("microbiome", "calculate_differential_abundance", input, result, options=opts)

def calculate_digital_ddh(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute digital DNA-DNA hybridisation (not implemented: needs GGDC formulas on BLAST/BLAT HSP tables)."""
    return candidate_not_implemented(
        "microbiome", "calculate_digital_ddh", input, options=options,
        requires="genome alignments (BLAST/BLAT/MUMmer HSPs) and the GGDC regression formulas",
        message="dDDH is not implemented; use calculate_ani for a k-mer identity estimate instead.",
    )

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign taxonomy to reads or contigs (not implemented: needs a classifier and reference database)."""
    return candidate_not_implemented(
        "microbiome", "classify", input, options=options,
        requires="a taxonomic classifier database (Kraken2, Kaiju, Centrifuge, GTDB-Tk; see classify_kraken2 etc.)",
        message="Native taxonomic classification is not implemented; use the classify_* wrappers with an installed database.",
    )

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute a between-sample beta-diversity distance matrix from an abundance table.

    ``input``: abundance table (``{sample: {taxon: count}}``, ``{"samples","taxa","counts"}``, rows of counts, or TSV
    text/path). ``options``: ``entity`` (only ``"beta_diversity"``); ``metric``: ``"braycurtis"`` (default;
    sum|a-b| / sum(a+b) on raw counts), ``"jaccard"`` (1 - shared/union of presence), ``"aitchison"`` (Euclidean
    distance of centred log-ratio profiles after adding ``pseudocount``, default 1.0) or ``"euclidean"``;
    ``relative`` (default false: convert each sample to relative abundance first).
    Result: ``samples``, ``metric`` and the symmetric ``distance_matrix``.
    """
    opts = dict(options or {})
    _impl.entity(opts, "beta_diversity", ("beta_diversity",))
    table = _impl.parse_table(input, opts)
    metric = _impl.opt_str(opts, "metric", "braycurtis", ("braycurtis", "jaccard", "aitchison", "euclidean"))
    pseudocount = _impl.opt_float(opts, "pseudocount", 1.0, 1e-12)
    if len(table.samples) < 2:
        raise InputFormatError("Beta diversity needs at least 2 samples.", code="MICROBIOME_TABLE_INVALID")
    if len(table.samples) > 2_000:
        raise ConfigurationError("At most 2000 samples are supported.", code="MICROBIOME_TABLE_TOO_LARGE")
    if _impl.opt_bool(opts, "relative", False):
        _impl.require_depth(table)
        table = _impl.Table(table.samples, table.taxa, [[x / sum(r) for x in r] for r in table.counts])
    result = {"metric": metric, "samples": table.samples, "distance_matrix": _impl.distance_matrix(table, metric, pseudocount)}
    return candidate_ok("microbiome", "compare", input, result, options=opts)

def dereplicate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Collapse redundant genome/sequence records into clusters with one representative each.

    ``input``: DNA sequences as a list of strings, a mapping ``{id: sequence}``, or FASTA text/path.
    ``options``: ``method``: ``"exact"`` (default; identical sequences, case-insensitive, optionally identical up
    to reverse complement via ``canonical`` default true) or ``"ani"`` (greedy: longest sequence first becomes a
    representative, later sequences join the first representative whose k-mer ANI estimate (see
    ``calculate_ani``; ``k`` default 16, 2..32) is >= ``min_ani`` (default 0.95), else start a new cluster).
    Result: ``clusters`` (``representative``, ``members``) and counts. The ``"ani"`` mode is a k-mer estimate, not
    an alignment-based dRep clustering.
    """
    opts = dict(options or {})
    method = _impl.opt_str(opts, "method", "exact", ("exact", "ani"))
    canonical = _impl.opt_bool(opts, "canonical", True)
    seqs = _impl.parse_sequences(input, minimum=1, maximum=2_000, max_bases=50_000_000)
    clusters: list[dict[str, object]] = []
    if method == "exact":
        index: dict[str, int] = {}
        for name, seq in seqs:
            key = min(seq, _impl.reverse_complement(seq)) if canonical else seq
            if key in index:
                clusters[index[key]]["members"].append(name)  # type: ignore[union-attr]
            else:
                index[key] = len(clusters)
                clusters.append({"representative": name, "members": [name]})
    else:
        k = _impl.opt_int(opts, "k", 16, 2, 32)
        min_ani = _impl.opt_float(opts, "min_ani", 0.95, 0.0, 1.0)
        ordered = sorted(seqs, key=lambda item: (-len(item[1]), item[0]))
        sets = {name: _impl.kmer_set(seq, k, canonical) for name, seq in ordered}
        for name, _seq in ordered:
            for cluster in clusters:
                ani = _impl.kmer_ani(sets[str(cluster["representative"])], sets[name], k)["ani_estimate"]
                if ani is not None and ani >= min_ani:
                    cluster["members"].append(name)  # type: ignore[union-attr]
                    break
            else:
                clusters.append({"representative": name, "members": [name]})
    result = {"method": method, "input_count": len(seqs), "cluster_count": len(clusters), "clusters": clusters}
    return candidate_ok("microbiome", "dereplicate", input, result, options=opts)

def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect taxon prevalence and the core microbiome (taxa present above a detection level in enough samples).

    ``input``: abundance table (see ``compare``). ``options``: ``entity``: ``"core_microbiome"`` (default; returns
    only the core taxa plus the full prevalence table) or ``"prevalence"``; ``detection`` (default 0.0): a taxon
    counts as present in a sample when its relative abundance is strictly greater than this value; ``prevalence``
    (default 0.5): core taxa have prevalence (fraction of samples where present) >= this value.
    Relative abundance of a taxon in a sample is count / sample total.
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "core_microbiome", ("core_microbiome", "prevalence"))
    detection = _impl.opt_float(opts, "detection", 0.0, 0.0, 1.0)
    threshold = _impl.opt_float(opts, "prevalence", 0.5, 0.0, 1.0)
    table = _impl.parse_table(input, opts)
    _impl.require_depth(table)
    n = len(table.samples)
    rel = [[x / sum(row) for x in row] for row in table.counts]
    rows = []
    for t, taxon in enumerate(table.taxa):
        present = sum(1 for s in range(n) if rel[s][t] > detection)
        rows.append({"taxon": taxon, "samples_present": present, "prevalence": present / n,
                     "mean_relative_abundance": sum(rel[s][t] for s in range(n)) / n})
    rows.sort(key=lambda r: (-r["prevalence"], -r["mean_relative_abundance"], r["taxon"]))
    result: dict[str, object] = {"detection": detection, "n_samples": n, "prevalence_table": rows}
    if which == "core_microbiome":
        result["prevalence_threshold"] = threshold
        result["core_taxa"] = [r["taxon"] for r in rows if r["prevalence"] >= threshold]
    return candidate_ok("microbiome", "detect", input, result, options=opts)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise an abundance table as QIIME-style TSV (taxa rows) or JSON, optionally writing it to a file.

    ``input``: abundance table (see ``compare``). ``options``: ``format`` ``"tsv"`` (default) or ``"json"``;
    ``output`` (optional file path: when given the text is written there, refusing to overwrite an existing
    file unless ``overwrite`` is true). The text is also returned in ``result["text"]``.
    """
    opts = dict(options or {})
    table = _impl.parse_table(input, opts)
    fmt = _impl.opt_str(opts, "format", "tsv", ("tsv", "json"))

    def fmt_num(x: float) -> str:
        return str(int(x)) if x == int(x) else repr(x)

    if fmt == "tsv":
        lines = ["#OTU ID\t" + "\t".join(table.samples)]
        for t, taxon in enumerate(table.taxa):
            lines.append(taxon + "\t" + "\t".join(fmt_num(table.counts[s][t]) for s in range(len(table.samples))))
        text = "\n".join(lines) + "\n"
    else:
        text = json.dumps(table.as_dict(), sort_keys=True)
    written = None
    output = opts.get("output")
    if output is not None:
        path = Path(str(output))
        if path.exists() and not _impl.opt_bool(opts, "overwrite", False):
            raise ConfigurationError(
                f"Output file {path} already exists.", code="MICROBIOME_OUTPUT_EXISTS", hint="Set options['overwrite']=True."
            )
        path.write_text(text, encoding="utf-8")
        written = str(path)
    return candidate_ok("microbiome", "export", input, {"format": fmt, "text": text, "output_path": written}, options=opts)

def infer_viral_quasispecies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Reconstruct viral quasispecies haplotypes (not implemented: needs a haplotype reconstruction program)."""
    return candidate_not_implemented(
        "microbiome", "infer_viral_quasispecies", input, options=options,
        requires="a viral haplotype reconstruction tool (QuRe, ShoRAH, PredictHaplo) and aligned reads",
        message="Quasispecies inference is not implemented natively.",
    )

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a stacked relative-abundance bar chart (one bar per sample) as an SVG string.

    ``input``: abundance table (see ``compare``). ``options``: ``top_n`` (default 8, 1..10): the most abundant taxa
    (by total count) get their own colour, the rest are pooled as ``Other``; ``output`` (optional ``.svg`` path
    to write). Each bar segment is a ``<rect class="bar">`` carrying ``data-sample``, ``data-taxon`` and
    ``data-fraction`` attributes. Result: ``{"svg": ...}``.
    """
    opts = dict(options or {})
    table = _impl.parse_table(input, opts)
    _impl.require_depth(table)
    top_n = _impl.opt_int(opts, "top_n", 8, 1, 10)
    svg = _impl.svg_stacked_bars(table, top_n)
    written = None
    if opts.get("output") is not None:
        Path(str(opts["output"])).write_text(svg, encoding="utf-8")
        written = str(opts["output"])
    return candidate_ok("microbiome", "plot", input, {"svg": svg, "output_path": written}, options=opts)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict microbiome phenotypes or functions (not implemented: needs a trained model or reference database)."""
    return candidate_not_implemented(
        "microbiome", "predict", input, options=options,
        requires="a trained model or reference database (e.g. PICRUSt2, BugBase)",
        message="Prediction from a microbiome profile needs trained reference data that is not bundled.",
    )

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quality-check an abundance table: per-sample depth, richness and zero fraction plus low-depth flags.

    ``input``: abundance table (see ``compare``). ``options``: ``entity`` (only ``"abundance_table"``);
    ``min_depth`` (default 1000): samples with a total count below this are flagged ``low_depth``.
    Result: per-sample ``total``, ``observed`` (taxa with count > 0), ``zero_fraction`` (share of taxa absent);
    table-level counts of samples/taxa, taxa never observed (``empty_taxa``) and ``flagged_samples``.
    """
    opts = dict(options or {})
    _impl.entity(opts, "abundance_table", ("abundance_table",))
    min_depth = _impl.opt_float(opts, "min_depth", 1000.0, 0.0)
    table = _impl.parse_table(input, opts)
    samples = []
    for name, row in zip(table.samples, table.counts):
        observed = sum(1 for x in row if x > 0)
        samples.append({"sample": name, "total": sum(row), "observed": observed,
                        "zero_fraction": 1 - observed / len(row), "low_depth": sum(row) < min_depth})
    empty = [taxon for t, taxon in enumerate(table.taxa) if all(row[t] == 0 for row in table.counts)]
    result = {
        "n_samples": len(table.samples), "n_taxa": len(table.taxa), "min_depth": min_depth, "samples": samples,
        "empty_taxa": empty, "flagged_samples": [s["sample"] for s in samples if s["low_depth"]],
    }
    return candidate_ok("microbiome", "qc", input, result, options=opts)

def quantify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quantify an abundance table: relative abundance, alpha diversity, rarefaction or centred log-ratio (CLR).

    ``input``: abundance table (see ``compare``); reads are not mapped here, counts must already exist.
    ``options``: ``entity``:
    ``"relative_abundance"`` (default; also ``"abundance"``) count / sample total;
    ``"alpha_diversity"`` per sample observed, Shannon (natural log unless ``log_base``), Simpson (1 - sum p^2),
    inverse Simpson, Pielou evenness (H / ln S_obs), bias-corrected Chao1 (S_obs + F1(F1-1)/(2(F2+1))) and Good's
    coverage (1 - F1/N), the latter two only for integer counts;
    ``"rarefaction"`` subsample each sample to ``depth`` reads without replacement using ``random.Random(seed)``
    (``seed`` default 0; samples below ``depth`` are dropped and listed unless ``drop_short`` is false);
    ``"clr"`` centred log-ratio of counts + ``pseudocount`` (default 1.0).
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "relative_abundance", ("relative_abundance", "abundance", "alpha_diversity", "rarefaction", "clr"))
    table = _impl.parse_table(input, opts)
    if which == "rarefaction":
        depth = _impl.opt_int(opts, "depth", 0, 1, _impl.MAX_READS_PER_SAMPLE) if "depth" in opts else None
        if depth is None:
            raise ConfigurationError("Rarefaction needs options['depth'].", code="MICROBIOME_OPTION_INVALID")
        rarefied, dropped = _impl.rarefy(table, depth, _impl.opt_int(opts, "seed", 0, 0, 2**63), _impl.opt_bool(opts, "drop_short", True))
        result: dict[str, object] = {"depth": depth, "dropped_samples": dropped, **rarefied.as_dict()}
        return candidate_ok("microbiome", "quantify", input, result, options=opts)
    _impl.require_depth(table)
    if which in ("relative_abundance", "abundance"):
        rel = [[x / sum(row) for x in row] for row in table.counts]
        result = {"samples": table.samples, "taxa": table.taxa, "relative_abundance": rel}
    elif which == "alpha_diversity":
        base = _impl.opt_float(opts, "log_base", math.e, 1.000001)
        result = {"samples": table.samples,
                  "metrics": {name: _impl.alpha_metrics(row, base) for name, row in zip(table.samples, table.counts)}}
    else:
        pc = _impl.opt_float(opts, "pseudocount", 1.0, 1e-12)
        result = {"samples": table.samples, "taxa": table.taxa, "clr": [_impl._clr(row, pc) for row in table.counts]}
    return candidate_ok("microbiome", "quantify", input, result, options=opts)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read an abundance table (QIIME/BIOM-TSV or CSV text, a file path, or in-memory structures) into samples x taxa counts.

    ``input``: table text, a path to a text/TSV/CSV file, ``{sample: {taxon: count}}``, ``{"samples","taxa","counts"}``
    or rows of counts. ``options``: ``orientation`` for text tables: ``"taxa_rows"`` (default; first column is the
    taxon id, remaining header cells are sample ids, lines starting with ``# `` are comments) or ``"samples_rows"``.
    Result: ``samples``, ``taxa``, ``counts`` (rows are samples). Values must be finite and non-negative.
    """
    opts = dict(options or {})
    table = _impl.parse_table(input, opts)
    return candidate_ok("microbiome", "read", input, table.as_dict(), options=opts)

def track(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Track strains or pathogens across samples (not implemented: needs strain-level markers or SNP profiles and a reference)."""
    return candidate_not_implemented(
        "microbiome", "track", input, options=options,
        requires="strain-level marker/SNP references (e.g. StrainPhlAn, inStrain) and epidemiological metadata",
        message="Strain/pathogen tracking is not implemented natively.",
    )

__all__ = ['annotate', 'assemble', 'bin', 'build_cooccurrence_network', 'build_microbial_gene_catalog', 'build_pathogen_transmission_tree', 'calculate_ani', 'calculate_differential_abundance', 'calculate_digital_ddh', 'classify', 'classify_centrifuge', 'classify_gtdbtk', 'classify_kaiju', 'classify_kraken2', 'compare', 'dereplicate', 'detect', 'estimate_bracken_abundance', 'export', 'infer_viral_quasispecies', 'plot', 'predict', 'profile_humann', 'qc', 'quantify', 'read', 'track']

__candidate_adapter_operations__ = set(['annotate', 'assemble', 'bin', 'build_cooccurrence_network', 'build_microbial_gene_catalog', 'build_pathogen_transmission_tree', 'calculate_ani', 'calculate_differential_abundance', 'calculate_digital_ddh', 'classify', 'compare', 'dereplicate', 'detect', 'export', 'infer_viral_quasispecies', 'plot', 'predict', 'qc', 'quantify', 'read', 'track'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    compare_minhash_sketches,
    estimate_unique_kmers,
    kmer_abundance_histogram,
    normalize_reads_by_median,
    query_kmer_counts,
    read_sourmash_signatures,
    sketch_fasta_minhash,
    write_sourmash_signatures,
)

__all__ = [*__all__, "compare_minhash_sketches", "estimate_unique_kmers", "kmer_abundance_histogram", "normalize_reads_by_median", "query_kmer_counts", "read_sourmash_signatures", "sketch_fasta_minhash", "write_sourmash_signatures"]

# Full-pass exports (see docs/api/features)
from .khmer_functions import (
    decode_kmer_2bit,
    encode_kmer_2bit,
    estimate_kmer_table_size,
    expected_kmer_table_collisions,
    hash_kmer_murmur3,
    inspect_khmer_graph,
    query_saved_countgraph,
    trim_low_abundance_reads,
)

from .sourmash_functions import (
    compare_minhash_matrix,
    downsample_minhash_sketch,
    estimate_mutation_rate,
    hash_canonical_kmer,
    intersect_minhash_sketches,
    lowest_common_lineage,
    search_minhash_database,
    set_size_error_probability,
)

__all__ = [*__all__, "decode_kmer_2bit", "encode_kmer_2bit", "estimate_kmer_table_size", "expected_kmer_table_collisions", "hash_kmer_murmur3", "inspect_khmer_graph", "query_saved_countgraph", "trim_low_abundance_reads", "compare_minhash_matrix", "downsample_minhash_sketch", "estimate_mutation_rate", "hash_canonical_kmer", "intersect_minhash_sketches", "lowest_common_lineage", "search_minhash_database", "set_size_error_probability"]

# Review-pass exports (see docs/api/features)
from .khmer_review_functions import (
    clean_khmer_read_sequences,
    read_paired_fastx,
)

from .sourmash_review_functions import (
    bp_threshold_to_hashes,
    count_kmer_lca_assignments,
    expected_mutated_kmers,
    kmer_mutation_probability,
    probability_no_shared_kmers,
    variance_mutated_kmers,
)

__all__ = [*__all__, "clean_khmer_read_sequences", "read_paired_fastx", "bp_threshold_to_hashes", "count_kmer_lca_assignments", "expected_mutated_kmers", "kmer_mutation_probability", "probability_no_shared_kmers", "variance_mutated_kmers"]

# Family-pass exports (see docs/api/features)
from .khmer_family_functions import (
    query_saved_khmer_table,
)

from .sourmash_family_functions import (
    compare_minhash_abundance_matrix,
    downsample_minhash_num_sketch,
    parse_sourmash_signature_json,
    search_minhash_database_ranked,
    search_sourmash_sbt_index,
    select_sourmash_signatures,
    write_sourmash_signature_archive,
)

__all__ = [*__all__, "query_saved_khmer_table", "compare_minhash_abundance_matrix", "downsample_minhash_num_sketch", "parse_sourmash_signature_json", "search_minhash_database_ranked", "search_sourmash_sbt_index", "select_sourmash_signatures", "write_sourmash_signature_archive"]

_register_public_verbose(globals())
