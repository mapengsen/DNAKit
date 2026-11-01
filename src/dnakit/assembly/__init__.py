"""Local assembly operations: native statistics and pangenome tools, plus honest ``not_implemented`` stubs.

Functions that need reads, alignments, databases or an assembler return
``status: "not_implemented"`` (see ``dnakit._candidate_runtime``).
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError, InputFormatError

from . import _impl_assembly as _impl
from .external import assemble_flye, assemble_hifiasm, assemble_megahit, assemble_shasta, assemble_unicycler, evaluate_merqury, run_busco
from .local import (
    calculate_l50 as _calculate_l50_native,
    calculate_n50 as _calculate_n50_native,
    calculate_ng50 as _calculate_ng50_native,
    calculate_nx as _calculate_nx_native,
    calculate_stats as _calculate_stats_native,
)

_MODULE = "assembly"


def _ok(operation: str, input: object, result: object, options: Mapping[str, object] | None, warnings: tuple[str, ...] = ()) -> object:
    return candidate_ok(_MODULE, operation, input, result, options=options, warnings=warnings)


def _not_implemented(operation: str, input: object, options: Mapping[str, object] | None, requires: str, message: str) -> object:
    return candidate_not_implemented(_MODULE, operation, input, options=options, requires=requires, message=message)


def _write_output(settings: Mapping[str, object], text: str) -> str | None:
    """Write ``text`` to ``options['output']`` only when the caller gave an explicit path."""

    target = settings.get("output")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or "\x00" in str(target):
        raise ConfigurationError("options.output must be a file path.", code="ASSEMBLY_OPTION")
    path = Path(target)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def add_genome_to_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Add one genome's gene set as a new column of a gene presence/absence matrix.

    input: the existing pangenome as {genome: [genes]}, {gene: {genome: 0/1}} or
    {"genomes": [...], "genes": [...], "matrix": [[...]]} (rows = genes, columns = genomes).
    options: genome (str, required, new unique genome name); genes (list of str, required,
    the genes present in the new genome). Genes not seen before become new all-absent rows
    for the old genomes. Returns the updated matrix, the number of new genes and the
    pangenome size before/after. No homology search is done: gene names are the clusters.
    """

    settings = _impl.options_dict(options)
    genomes, genes, rows = _impl.presence_matrix(input)
    name = settings.get("genome")
    new_genes = settings.get("genes")
    if not isinstance(name, str) or not name:
        raise ConfigurationError("options.genome must be the new genome's name.", code="ASSEMBLY_OPTION", hint="Pass options={'genome': 'G9', 'genes': [...]}.")
    if not isinstance(new_genes, (list, tuple, set, frozenset)) or not all(isinstance(g, str) for g in new_genes):
        raise ConfigurationError("options.genes must be a list of gene names.", code="ASSEMBLY_OPTION")
    if name in genomes:
        raise ConfigurationError(f"genome {name!r} is already in the pangenome.", code="ASSEMBLY_OPTION", hint="Use a different genome name.")
    present = set(new_genes)
    updated_rows = [[*row, 1 if gene in present else 0] for gene, row in zip(genes, rows)]
    added = sorted(present - set(genes))
    updated_genes = [*genes, *added]
    updated_rows.extend([[0] * len(genomes) + [1] for _ in added])
    result = {
        **_impl.matrix_dict([*genomes, name], updated_genes, updated_rows),
        "new_gene_count": len(added),
        "new_genes": added,
        "pangenome_size_before": len(genes),
        "pangenome_size_after": len(updated_genes),
    }
    return _ok("add_genome_to_pangenome", input, result, options)


def build_chromosome_level_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: chromosome-level assembly needs Hi-C/optical-map/synteny scaffolding evidence."""

    return _not_implemented(
        "build_chromosome_level_assembly", input, options, "scaffolding evidence (Hi-C contacts, optical maps or a reference) and a scaffolder such as YaHS/SALSA",
        "Ordering and orienting contigs into chromosomes needs linking evidence that is not part of the input; use scaffold() with an explicit layout instead.",
    )


_COMMAND_TOOLS = {
    "flye": "flye",
    "hifiasm": "hifiasm",
    "megahit": "megahit",
    "shasta": "shasta",
    "unicycler": "unicycler",
}


def _tool_settings(input: object, options: Mapping[str, object] | None) -> dict[str, object]:
    settings: dict[str, object] = dict(input) if isinstance(input, Mapping) else {}
    settings.update(_impl.options_dict(options))
    tool = settings.get("tool")
    if tool not in _COMMAND_TOOLS:
        raise ConfigurationError(f"tool must be one of {', '.join(sorted(_COMMAND_TOOLS))}; got {tool!r}.", code="ASSEMBLY_OPTION", hint="Pass {'tool': 'flye', 'reads': ..., 'output_directory': ...}.")
    for key in ("output_directory",):
        if not isinstance(settings.get(key), (str, Path)) or not str(settings[key]):
            raise ConfigurationError(f"{key} is required.", code="ASSEMBLY_OPTION")
    return settings


def _threads(settings: Mapping[str, object]) -> int:
    value = _impl.get_int(settings, "threads", 1, low=1, high=256)
    return 1 if value is None else value


def _required_path(settings: Mapping[str, object], key: str) -> str:
    value = settings.get(key)
    if not isinstance(value, (str, Path)) or not str(value) or "\x00" in str(value):
        raise ConfigurationError(f"{key} is required for tool {settings.get('tool')!r}.", code="ASSEMBLY_OPTION")
    return str(value)


def _assembler_argv(settings: Mapping[str, object]) -> list[str]:
    """Argument vector identical to the one used by the ``assemble_*`` wrappers (no file system access)."""

    tool = str(settings["tool"])
    program = str(settings.get("program", _COMMAND_TOOLS[tool]))
    out = Path(str(settings["output_directory"]))
    threads = _threads(settings)
    if tool == "flye":
        platform = str(settings.get("platform", "nano-hq"))
        if platform not in {"nano-raw", "nano-hq", "pacbio-raw", "pacbio-hifi", "pacbio-corr"}:
            raise ConfigurationError("platform must be a supported Flye read mode.", code="ASSEMBLY_OPTION")
        genome_size = settings.get("genome_size")
        if not isinstance(genome_size, str) or not genome_size:
            raise ConfigurationError("genome_size (for example '5m') is required for flye.", code="ASSEMBLY_OPTION")
        return [program, f"--{platform}", _required_path(settings, "reads"), "--genome-size", genome_size, "--out-dir", str(out), "--threads", str(threads)]
    if tool == "hifiasm":
        prefix = str(settings.get("prefix", "assembly"))
        return [program, "-o", str(out / prefix), "-t", str(threads), _required_path(settings, "reads")]
    if tool == "megahit":
        argv = [program, "-o", str(out / "megahit"), "-t", str(threads)]
        min_len = _impl.get_int(settings, "min_contig_length", None, low=1, high=10_000_000)
        if min_len is not None:
            argv += ["--min-contig-len", str(min_len)]
        if settings.get("read2") is None:
            return [*argv, "-r", _required_path(settings, "read1")]
        return [*argv, "-1", _required_path(settings, "read1"), "-2", _required_path(settings, "read2")]
    if tool == "shasta":
        argv = [program, "--input", _required_path(settings, "reads"), "--config", str(settings.get("config", "Nanopore-May2022")), "--assemblyDirectory", str(out)]
        if settings.get("threads") is not None:
            argv += ["--threads", str(threads)]
        return argv
    # unicycler
    mode = str(settings.get("mode", "normal"))
    if mode not in {"conservative", "normal", "bold"}:
        raise ConfigurationError("mode must be conservative, normal or bold.", code="ASSEMBLY_OPTION")
    argv = [program]
    if settings.get("read1") is not None:
        argv += ["-1", _required_path(settings, "read1")]
    if settings.get("read2") is not None:
        argv += ["-2", _required_path(settings, "read2")]
    if settings.get("long_reads") is not None:
        argv += ["-l", _required_path(settings, "long_reads")]
    if len(argv) == 1:
        raise ConfigurationError("unicycler needs read1 or long_reads.", code="ASSEMBLY_OPTION")
    return [*argv, "-o", str(out / "unicycler"), "-m", mode, "-t", str(threads)]


def build_command(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build (without running) the command line of one supported assembler as an argument list.

    Supported tools: flye, hifiasm, megahit, shasta, unicycler. The argument vectors are the same
    as those used by assemble_flye/assemble_hifiasm/assemble_megahit/assemble_shasta/assemble_unicycler.
    input: a mapping with the settings below (options are merged over it and win).
    keys: tool (required); output_directory (required); program (executable name, default = tool);
    reads (flye, hifiasm, shasta); read1/read2 (megahit, unicycler); long_reads (unicycler);
    genome_size (flye, e.g. '5m'); platform (flye, default 'nano-hq'); config (shasta); mode
    (unicycler); prefix (hifiasm, default 'assembly'); threads (1-256, default 1);
    min_contig_length (megahit). No file is read or created and nothing is executed.
    """

    settings = _tool_settings(input, options)
    argv = _assembler_argv(settings)
    result = {"tool": settings["tool"], "argv": argv, "command_line": " ".join(_quote(part) for part in argv)}
    return _ok("build_command", input, result, options)


def _quote(part: str) -> str:
    import shlex

    return shlex.quote(part)


def build_graph_genome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a compacted forward-strand de Bruijn graph (unitigs, GFA1) from input sequences.

    input: FASTA text/path, {name: sequence} or a list of sequences (each is one 'colour').
    options: k (int 3-127, default 21); output (optional path where the GFA1 text is written).
    Nodes are k-mers, edges are observed (k-1)-overlaps, unitigs are maximal non-branching
    paths; reverse-complement strands are NOT merged (it is not a bidirected graph) and no
    variation-graph/alignment step is performed. Each unitig lists the input sequences it occurs in.
    Input is limited to 2,000,000 bases.
    """

    settings = _impl.options_dict(options)
    k = _impl.get_int(settings, "k", 21, low=3, high=127) or 21
    contigs = _impl.parse_contigs(input, require_sequences=True)
    result = _impl.de_bruijn_unitigs(contigs, k)
    result["output_path"] = _write_output(settings, result["gfa"])
    return _ok("build_graph_genome", input, result, options)


def build_haplotype_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: haplotype-resolved assembly needs reads (HiFi/ONT, Hi-C or parental k-mers) and an assembler."""

    return _not_implemented(
        "build_haplotype_assembly", input, options, "an assembler with phasing support (for example hifiasm) and sequencing reads",
        "Haplotype-resolved assembly cannot be derived from contigs alone; run assemble_hifiasm on reads instead.",
    )


def build_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a gene presence/absence pangenome matrix with core/shell/cloud counts from per-genome gene sets.

    input: {genome: [gene cluster names]} (also accepts the gene->{genome: 0/1} and
    genomes/genes/matrix forms). Gene names are taken as already-clustered orthology groups; no
    sequence clustering is performed. options: core_threshold (fraction in (0, 1], default 0.99;
    a gene is core if present in >= ceil(threshold * genomes) genomes); cloud_threshold (default
    0.15; accessory genes below this fraction are 'cloud', the rest 'shell'). Result: matrix,
    pangenome size, partition counts, per-genome gene and unique-gene counts.
    """

    settings = _impl.options_dict(options)
    genomes, genes, rows = _impl.presence_matrix(input)
    partition = _impl.pangenome_partition(
        genomes, genes, rows, _impl.core_threshold_option(settings), _impl.get_float(settings, "cloud_threshold", 0.15, low=0.0, high=1.0)
    )
    unique = {g: sum(1 for row in rows if row[j] and sum(row) == 1) for j, g in enumerate(genomes)}
    result = {
        **_impl.matrix_dict(genomes, genes, rows),
        "pangenome_size": len(genes),
        "core_count": partition["core_count"],
        "shell_count": len(partition["shell_genes"]),
        "cloud_count": len(partition["cloud_genes"]),
        "core_threshold": partition["core_threshold"],
        "genes_per_genome": {g: sum(row[j] for row in rows) for j, g in enumerate(genomes)},
        "unique_genes_per_genome": unique,
    }
    return _ok("build_pangenome", input, result, options)


def calculate_assembly_completeness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate assembly completeness as length ratio and/or reference k-mer recovery.

    input: assembly (FASTA text/path, {name: sequence}, list of sequences/lengths).
    options (at least one of reference_length / reference): reference_length (positive number,
    expected genome size; gives length_ratio = assembly_length / reference_length and
    length_completeness = min(ratio, 1)); reference (assembly-like with sequences; gives
    kmer_completeness = fraction of the reference's distinct canonical k-mers present in the
    assembly); k (3-63, default 21). This is a length / k-mer measure, not gene-content
    completeness (see calculate_busco_completeness).
    """

    settings = _impl.options_dict(options)
    contigs = _impl.parse_contigs(input)
    return _ok("calculate_assembly_completeness", input, _impl.completeness(contigs, settings), options)


def calculate_busco_completeness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute BUSCO completeness percentages from BUSCO counts or a BUSCO short-summary line.

    input: a mapping with single (S), duplicated (D), fragmented (F), missing (M) integer counts,
    or the text/path of a BUSCO short summary containing 'C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614'.
    complete = S + D; percentages are relative to n = S + D + F + M (or the stated n). This only
    summarises existing BUSCO results; it does not run BUSCO (see run_busco). options: none.
    """

    return _ok("calculate_busco_completeness", input, _impl.busco_completeness(input), options)


def calculate_pangenome_openness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Fit Heaps' law to the number of new genes per added genome to classify a pangenome as open or closed.

    input: gene presence/absence data (see build_pangenome; at least 4 genomes). Model: the mean
    number of new genes found when the N-th genome is added is kappa * N^(-alpha); log-log OLS
    over N = 2..G; the pangenome is open if alpha < 1 and closed if alpha >= 1 (Tettelin 2008). The mean is taken over
    random genome orders. options: n_permutations (int 0-10000, default 100; 0 = use the given
    genome order once); seed (int, default 0). Deterministic for a fixed seed.
    """

    settings = _impl.options_dict(options)
    permutations = _impl.get_int(settings, "n_permutations", 100, low=0, high=10_000) or 0
    seed = _impl.get_int(settings, "seed", 0, low=0, high=2**32 - 1) or 0
    genomes, genes, rows = _impl.presence_matrix(input)
    return _ok("calculate_pangenome_openness", input, _impl.pangenome_openness(genomes, genes, rows, permutations, seed), options)


def check_input(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate assembly input (FASTA text/path, {name: sequence}, list of sequences or lengths) and report problems.

    Reports the detected format, contig count, total length, duplicate contig names and
    non-IUPAC symbols. An unusable input is reported as valid=False with the reasons rather than
    raising. options: none.
    """

    return _ok("check_input", input, _impl.check_assembly_input(input), options)


def _assembly_pair(input: object, settings: Mapping[str, object]) -> tuple[object, object]:
    if settings.get("other") is not None:
        return input, settings["other"]
    if isinstance(input, Mapping) and {"left", "right"} <= set(input):
        return input["left"], input["right"]
    if isinstance(input, Mapping) and {"a", "b"} <= set(input):
        return input["a"], input["b"]
    if isinstance(input, (list, tuple)) and len(input) == 2 and all(isinstance(x, (str, Path, Mapping, list, tuple)) for x in input):
        return input[0], input[1]
    raise InputFormatError("compare_assemblies needs two assemblies.", code="ASSEMBLY_COMPARE_INPUT", hint="Pass [assembly_a, assembly_b], {'left': .., 'right': ..} or options={'other': assembly_b}.")


def compare_assemblies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two assemblies by contiguity statistics (count, length, N50, auN, GC, gaps) and shared contig sequences.

    input: [assembly_a, assembly_b] or {'left': a, 'right': b} (or one assembly plus options.other);
    each assembly is FASTA text/path, {name: sequence}, or a list of lengths/sequences.
    options: other (second assembly); genome_size (optional, adds NG50). Shared contigs are
    contigs with identical sequence on either strand (only when both have sequences); no
    alignment is performed.
    """

    settings = _impl.options_dict(options)
    left_raw, right_raw = _assembly_pair(input, settings)
    left = _impl.parse_contigs(left_raw, what="left assembly")
    right = _impl.parse_contigs(right_raw, what="right assembly")
    genome_size = _impl.get_int(settings, "genome_size", None, low=1, high=10**13)
    ls, rs = _impl.assembly_stats(left, genome_size), _impl.assembly_stats(right, genome_size)
    result: dict[str, object] = {
        "left": ls,
        "right": rs,
        "difference_right_minus_left": {
            key: rs[key] - ls[key] for key in ("contig_count", "total_length", "n50", "l50", "max_length") if key in ls and key in rs
        },
        "n50_ratio_right_to_left": rs["n50"] / ls["n50"],
    }
    if all(seq is not None for _, seq, _ in left + right):
        lkeys = {_impl.canonical_sequence(seq or "") for _, seq, _ in left}
        rkeys = {_impl.canonical_sequence(seq or "") for _, seq, _ in right}
        shared = lkeys & rkeys
        result["shared_contig_sequences"] = len(shared)
        result["shared_bases"] = sum(len(s) for s in shared)
        result["left_only_contigs"] = len(lkeys - rkeys)
        result["right_only_contigs"] = len(rkeys - lkeys)
    return _ok("compare_assemblies", input, result, options)


def compare_assembly_graphs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two GFA1 assembly graphs by segment and link content.

    input: [gfa_a, gfa_b] (GFA1 text or paths) or one GFA plus options.other. Segments are
    matched by canonical sequence (either strand) when every segment has a sequence, otherwise by
    name; links are matched on matched segments, a link and its reverse-complement counterpart
    counting as the same. Reports counts, shared/unique segments and links and Jaccard indices.
    """

    settings = _impl.options_dict(options)
    if settings.get("other") is not None:
        left_raw, right_raw = input, settings["other"]
    elif isinstance(input, (list, tuple)) and len(input) == 2:
        left_raw, right_raw = input
    else:
        raise InputFormatError("compare_assembly_graphs needs two GFA graphs.", code="ASSEMBLY_COMPARE_INPUT", hint="Pass [gfa_a, gfa_b] or options={'other': gfa_b}.")
    result = _impl.compare_gfa(_impl.load_gfa(left_raw), _impl.load_gfa(right_raw))
    return _ok("compare_assembly_graphs", input, result, options)


def detect_assembly_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: contamination screening needs a taxonomy/vector/adapter database or read coverage."""

    return _not_implemented(
        "detect_assembly_contamination", input, options, "a taxonomic or contaminant reference database (for example BLAST nt, Kraken2, FCS-GX)",
        "Deciding whether a contig is foreign requires comparison against reference databases; a GC or length outlier is not a contamination call.",
    )


def detect_chimeric_contigs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: chimeric contig detection needs read or reference alignments."""

    return _not_implemented(
        "detect_chimeric_contigs", input, options, "read alignments (BAM/PAF) to the contigs or a reference assembly",
        "Chimeras are detected from split or discordant alignments, which are not part of the input.",
    )


def detect_misassemblies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: misassembly detection needs a reference alignment or read-mapping evidence."""

    return _not_implemented(
        "detect_misassemblies", input, options, "a reference genome alignment (QUAST/MUMmer) or read mappings (CRAQ/Inspector)",
        "Misassemblies are defined relative to a reference or to read evidence, neither of which is part of the input.",
    )


def estimate_heterozygosity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate the SNP heterozygosity rate of a diploid from its k-mer spectrum (simple two-peak model).

    input: k-mer histogram ({coverage: distinct k-mers}, [[coverage, count], ...] or 'coverage count'
    text as written by jellyfish histo), or sequences/reads (counted, limited to 5,000,000 bases).
    Heuristic model: k-mers with coverage between min_coverage and het_upper*peak are heterozygous
    (count H), those in (het_upper*peak, hom_upper*peak] are homozygous (count O); isolated SNPs
    give H = 2*k*r*L and O = L*(1-k*r), so r = R/(k*(2+R)) with R = H/O. Ignores repeats, errors
    above the cutoff and clustered SNPs. options: k (default 21, needed to interpret the
    histogram), min_coverage (error cutoff; default = first valley), peak_coverage (homozygous
    peak; default = highest bin after the valley), het_upper (default 0.75), hom_upper (default 1.5).
    """

    settings = _impl.options_dict(options)
    hist = _impl.kmer_histogram(input, settings)
    return _ok("estimate_heterozygosity", input, _impl.estimate_het(hist, settings), options, ("heuristic two-peak model; see docstring for assumptions",))


def estimate_repeats(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate repeat content from k-mer multiplicity: k-mers at more than repeat_factor times the unique peak are repeats.

    input: k-mer histogram (see estimate_heterozygosity) or sequences/reads. With T = sum_{x >=
    min_coverage} x*n(x) and R = sum_{x > repeat_factor*peak} x*n(x), repeat_fraction = R/T (share
    of genome positions in multi-copy k-mers) and the estimated copy number of a k-mer is round(x/peak).
    Heuristic: assumes uniform coverage. options: min_coverage (default = first valley),
    peak_coverage (default = highest bin after the valley), repeat_factor (default 1.5), k (for sequences).
    """

    settings = _impl.options_dict(options)
    hist = _impl.kmer_histogram(input, settings)
    return _ok("estimate_repeats", input, _impl.estimate_repeat_content(hist, settings), options, ("k-mer multiplicity estimate; assumes uniform coverage",))


def estimate_size(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate genome size from a k-mer histogram as total k-mer instances divided by the main coverage peak.

    Model (Liu et al. 2013): G = sum_{x >= min_coverage} x*n(x) / peak_coverage, where n(x) is the
    number of distinct k-mers seen x times. input: k-mer histogram ({coverage: count},
    [[coverage, count], ...], 'coverage count' text) or sequences/reads (counted, 5,000,000 bases max).
    options: min_coverage (error cutoff; default = first local valley of the histogram),
    peak_coverage (default = highest bin above the cutoff), max_coverage (optional upper cutoff), k
    (default 21, only for sequence input). Result is in bases.
    """

    settings = _impl.options_dict(options)
    hist = _impl.kmer_histogram(input, settings)
    return _ok("estimate_size", input, _impl.estimate_genome_size(hist, settings), options)


def extract_accessory_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """List accessory pangenome genes: those present in fewer than core_threshold of the genomes, split into shell and cloud.

    input: gene presence/absence data (see build_pangenome). options: core_threshold (default 0.99);
    cloud_threshold (default 0.15; accessory genes below this genome fraction are cloud, others shell).
    """

    settings = _impl.options_dict(options)
    genomes, genes, rows = _impl.presence_matrix(input)
    partition = _impl.pangenome_partition(
        genomes, genes, rows, _impl.core_threshold_option(settings), _impl.get_float(settings, "cloud_threshold", 0.15, low=0.0, high=1.0)
    )
    keep = ("genome_count", "total_genes", "core_threshold", "accessory_genes", "accessory_count", "shell_genes", "cloud_genes", "cloud_threshold")
    result = {key: partition[key] for key in keep}
    result["gene_fractions"] = {g: partition["gene_fractions"][g] for g in partition["accessory_genes"]}
    return _ok("extract_accessory_pangenome", input, result, options)


def extract_core_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """List core pangenome genes: those present in at least core_threshold of the genomes.

    input: gene presence/absence data (see build_pangenome). options: core_threshold (fraction in
    (0, 1], default 0.99); a gene is core if present in >= ceil(core_threshold * genome_count) genomes.
    """

    settings = _impl.options_dict(options)
    genomes, genes, rows = _impl.presence_matrix(input)
    partition = _impl.pangenome_partition(genomes, genes, rows, _impl.core_threshold_option(settings))
    keep = ("genome_count", "total_genes", "core_threshold", "core_min_genomes", "core_genes", "core_count")
    result = {key: partition[key] for key in keep}
    result["gene_fractions"] = {g: partition["gene_fractions"][g] for g in partition["core_genes"]}
    return _ok("extract_core_pangenome", input, result, options)


def identify_centromeres(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: centromere calling has no precise sequence-only rule and needs satellite/repeat annotation."""

    return _not_implemented(
        "identify_centromeres", input, options, "tandem-repeat/satellite annotation (TRF, SRF, RepeatMasker) and species-specific centromere repeat knowledge",
        "Centromere positions are species-specific (satellite monomers, CENH3 ChIP, Hi-C); no universal documented rule exists to apply to bare contigs.",
    )


def identify_repeats_in_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find repeated regions of an assembly as intervals covered by k-mers that occur at least min_copies times.

    input: assembly with sequences (FASTA text/path, {name: sequence}, list of sequences). Canonical
    k-mers (a k-mer and its reverse complement are the same) are counted over all contigs; every
    position covered by a k-mer with count >= min_copies is a repeat; overlapping/adjacent
    k-mers are merged into 0-based half-open intervals. No repeat library or classification is used.
    options: k (int 3-63, default 21); min_copies (int >= 2, default 2); min_length (default = k,
    shorter intervals are dropped). Limited to 5,000,000 bases.
    """

    settings = _impl.options_dict(options)
    k = _impl.get_int(settings, "k", 21, low=3, high=63) or 21
    copies = _impl.get_int(settings, "min_copies", 2, low=2, high=10**6) or 2
    min_length = _impl.get_int(settings, "min_length", k, low=1, high=10**9) or k
    contigs = _impl.parse_contigs(input, require_sequences=True)
    return _ok("identify_repeats_in_assembly", input, _impl.repeats_in_assembly(contigs, k, copies, min_length), options)


def identify_telomeres(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find telomeric repeat arrays at contig ends by perfect tandem repeats of a telomere motif.

    input: assembly with sequences. Rule: a telomere end is a perfect tandem run of at least
    min_repeats copies of the motif or its reverse complement (any rotation) that begins within
    max_distance bases of the contig's start (left) or ends within max_distance bases of its end (right),
    assigned to the nearer end (ties: left).
    options: motif (default 'TTAGGG', the vertebrate repeat; plants use 'TTTAGGG'); min_repeats (default 5);
    max_distance (default 100). Perfect repeats only: degenerate or short arrays are missed.
    """

    settings = _impl.options_dict(options)
    motif = str(settings.get("motif", "TTAGGG")).upper()
    min_repeats = _impl.get_int(settings, "min_repeats", 5, low=2, high=10_000) or 5
    max_distance = _impl.get_int(settings, "max_distance", 100, low=0, high=10_000_000) or 0
    contigs = _impl.parse_contigs(input, require_sequences=True)
    return _ok("identify_telomeres", input, _impl.telomeres(contigs, motif, min_repeats, max_distance), options)


def parse_log(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise an assembler log by severity level (counts of ERROR/WARNING/INFO/DEBUG lines and the first error and warning lines).

    input: log text or path. A line's level is the first of CRITICAL/FATAL, ERROR, WARNING/WARN, INFO, DEBUG
    found in it (case-insensitive, whole word). No tool-specific fields are parsed.
    options: max_lines (int 1-1000, default 20; number of error and of warning lines returned).
    """

    settings = _impl.options_dict(options)
    max_lines = _impl.get_int(settings, "max_lines", 20, low=1, high=1000) or 20
    text = _impl.read_text_input(input, what="log")
    if text is None:
        raise InputFormatError("log must be text or a path.", code="ASSEMBLY_INPUT_FORMAT")
    return _ok("parse_log", input, _impl.parse_log_text(text, max_lines), options)


def parse_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse assembler output (FASTA contigs or GFA1) into per-record length/GC/coverage fields and assembly statistics.

    input: FASTA or GFA1 text or path. SPAdes headers (NODE_1_length_100_cov_5.2) and MEGAHIT headers
    (flag=1 multi=2.0 len=300) are recognised and their coverage/multiplicity returned. For GFA the
    segment length (sequence or LN tag) and depth (DP/dp/KC tag) are returned. options: none.
    """

    return _ok("parse_output", input, _impl.parse_assembler_output(input), options)


def phase_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: phasing an assembly needs reads, Hi-C or parental data."""

    return _not_implemented(
        "phase_assembly", input, options, "long reads, Hi-C or trio k-mers and a phasing assembler (hifiasm, HapCUT2, WhatsHap)",
        "Assigning contigs or bases to haplotypes needs linking evidence that is not part of the input.",
    )


def plot_graph(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a GFA1 assembly graph as an SVG with segments on a circle and links as straight lines.

    input: GFA1 text or path (at most 2,000 segments). Node radius grows with log10 of segment
    length; the layout is deterministic (sorted segment names around the circle) and carries no
    biological meaning. options: output (optional path; the SVG text is written there).
    """

    settings = _impl.options_dict(options)
    graph = _impl.load_gfa(input)
    svg = _impl.graph_svg(graph)
    result = {"svg": svg, "node_count": len(graph["segments"]), "edge_count": len(graph["links"]), "output_path": _write_output(settings, svg)}
    return _ok("plot_graph", input, result, options)


def polish(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: polishing needs reads aligned to the assembly and a polisher (Racon, Medaka, Pilon)."""

    return _not_implemented(
        "polish", input, options, "sequencing reads, an aligner and a polishing tool (Racon, Medaka, Pilon, NextPolish)",
        "Consensus correction is computed from read pileups, which are not part of the input.",
    )


def purge_haplotigs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: purging haplotigs needs read-depth and self-alignment evidence (purge_dups, Purge_Haplotigs)."""

    return _not_implemented(
        "purge_haplotigs", input, options, "read coverage and contig self-alignments (purge_dups, Purge_Haplotigs)",
        "Haplotigs are identified from half-coverage and allelic alignment evidence that is not part of the input.",
    )


def remove_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Remove the contigs named in options.contaminants from an assembly (detection is not performed).

    input: assembly with sequences or lengths. options: contaminants (list of contig names to
    remove, required); names not found in the assembly are reported in `not_found`, not an error.
    Result: kept and removed contig names, statistics before and after, and FASTA text of the
    kept contigs when sequences are available. Call a contamination screen first to get the list.
    """

    settings = _impl.options_dict(options)
    names = settings.get("contaminants")
    if not isinstance(names, (list, tuple, set, frozenset)) or not all(isinstance(n, str) for n in names):
        raise ConfigurationError("options.contaminants must be a list of contig names.", code="ASSEMBLY_OPTION", hint="Pass options={'contaminants': ['contig_3']}.")
    contigs = _impl.parse_contigs(input)
    drop = set(names)
    kept = [c for c in contigs if c[0] not in drop]
    removed = [c[0] for c in contigs if c[0] in drop]
    if not kept:
        raise ConfigurationError("removing these contigs would leave an empty assembly.", code="ASSEMBLY_OPTION")
    result: dict[str, object] = {
        "kept": [c[0] for c in kept],
        "removed": removed,
        "not_found": sorted(drop - {c[0] for c in contigs}),
        "stats_before": _impl.assembly_stats(contigs),
        "stats_after": _impl.assembly_stats(kept),
    }
    if all(seq is not None for _, seq, _ in kept):
        result["fasta"] = _impl.fasta_text((name, seq or "") for name, seq, _ in kept)
    return _ok("remove_contamination", input, result, options)


def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run one caller-installed assembler (flye, hifiasm, megahit, shasta or unicycler) through the matching assemble_* wrapper.

    input: a mapping of settings (options are merged over it): tool, output_directory, and the tool's
    inputs/parameters exactly as in build_command (reads, read1, read2, long_reads, genome_size, platform,
    config, mode, prefix, threads, min_contig_length, program) plus optional timeout_seconds. If the
    executable is not on PATH the result is status 'not_implemented' naming the program. Output is
    written only under output_directory.
    """

    settings = _tool_settings(input, options)
    tool = str(settings["tool"])
    argv = _assembler_argv(settings)  # validates the settings and the tool-specific required keys
    program = argv[0]
    if shutil.which(program) is None:
        return _not_implemented("run", input, options, f"{program} executable on PATH", f"The assembler program {program!r} is not installed or not on PATH, so nothing was run.")
    common = {"program": program, "verbose": False}
    if settings.get("timeout_seconds") is not None:
        common["timeout_seconds"] = float(settings["timeout_seconds"])  # type: ignore[arg-type]
    out = str(settings["output_directory"])
    threads = _threads(settings)
    if tool == "flye":
        outcome = assemble_flye(_required_path(settings, "reads"), out, genome_size=str(settings["genome_size"]), platform=str(settings.get("platform", "nano-hq")), threads=threads, **common)  # type: ignore[arg-type]
    elif tool == "hifiasm":
        outcome = assemble_hifiasm(_required_path(settings, "reads"), out, threads=threads, prefix=str(settings.get("prefix", "assembly")), **common)  # type: ignore[arg-type]
    elif tool == "megahit":
        read2 = settings.get("read2")
        outcome = assemble_megahit(_required_path(settings, "read1"), out, read2=None if read2 is None else str(read2), threads=threads, min_contig_length=_impl.get_int(settings, "min_contig_length", None, low=1, high=10_000_000), **common)  # type: ignore[arg-type]
    elif tool == "shasta":
        outcome = assemble_shasta(_required_path(settings, "reads"), out, config=str(settings.get("config", "Nanopore-May2022")), threads=None if settings.get("threads") is None else threads, **common)  # type: ignore[arg-type]
    else:
        outcome = assemble_unicycler(
            None if settings.get("read1") is None else str(settings["read1"]), out,
            read2=None if settings.get("read2") is None else str(settings["read2"]),
            long_reads=None if settings.get("long_reads") is None else str(settings["long_reads"]),
            mode=str(settings.get("mode", "normal")), threads=threads, **common,  # type: ignore[arg-type]
        )
    return _ok("run", input, outcome, options)


def scaffold(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Join contigs into scaffolds following an explicit, caller-supplied layout (no linking evidence is inferred).

    input: assembly with sequences. options: layout (required): list of scaffolds
    {'name': str, 'contigs': [contig_name | {'contig': name, 'orientation': '+'|'-', 'gap': int}, ...]};
    'gap' is the number of N between this contig and the next (default `gap`). '-' reverse-complements
    the contig. gap (int >= 0, default 100). Result: scaffold sequences, an AGP-like table (0-based
    half-open coordinates) and the unplaced contigs. Ordering/orientation is NOT predicted.
    """

    settings = _impl.options_dict(options)
    gap = _impl.get_int(settings, "gap", 100, low=0, high=10_000_000) or 0
    contigs = _impl.parse_contigs(input, require_sequences=True)
    return _ok("scaffold", input, _impl.scaffold_layout(contigs, settings.get("layout"), gap), options)


# Promote deterministic contiguity statistics to the formal local namespace;
# all other assembly operations remain dependency-checked candidate adapters.
calculate_l50 = _calculate_l50_native
calculate_n50 = _calculate_n50_native
calculate_ng50 = _calculate_ng50_native
calculate_nx = _calculate_nx_native
calculate_stats = _calculate_stats_native

__all__ = ['add_genome_to_pangenome', 'assemble_flye', 'assemble_hifiasm', 'assemble_megahit', 'assemble_shasta', 'assemble_unicycler', 'build_chromosome_level_assembly', 'build_command', 'build_graph_genome', 'build_haplotype_assembly', 'build_pangenome', 'calculate_assembly_completeness', 'calculate_busco_completeness', 'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_pangenome_openness', 'calculate_stats', 'check_input', 'compare_assemblies', 'compare_assembly_graphs', 'detect_assembly_contamination', 'detect_chimeric_contigs', 'detect_misassemblies', 'estimate_heterozygosity', 'estimate_repeats', 'estimate_size', 'extract_accessory_pangenome', 'extract_core_pangenome', 'evaluate_merqury', 'identify_centromeres', 'identify_repeats_in_assembly', 'identify_telomeres', 'parse_log', 'parse_output', 'phase_assembly', 'plot_graph', 'polish', 'purge_haplotigs', 'remove_contamination', 'run', 'run_busco', 'scaffold']

__candidate_adapter_operations__ = set(['add_genome_to_pangenome', 'build_chromosome_level_assembly', 'build_command', 'build_graph_genome', 'build_haplotype_assembly', 'build_pangenome', 'calculate_assembly_completeness', 'calculate_busco_completeness', 'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_pangenome_openness', 'calculate_stats', 'check_input', 'compare_assemblies', 'compare_assembly_graphs', 'detect_assembly_contamination', 'detect_chimeric_contigs', 'detect_misassemblies', 'estimate_heterozygosity', 'estimate_repeats', 'estimate_size', 'extract_accessory_pangenome', 'extract_core_pangenome', 'identify_centromeres', 'identify_repeats_in_assembly', 'identify_telomeres', 'parse_log', 'parse_output', 'phase_assembly', 'plot_graph', 'polish', 'purge_haplotigs', 'remove_contamination', 'run', 'run_busco', 'scaffold']) - {'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_stats'}
# Public integration of the reviewed P00–P45 operation packs.
from .busco_p15_cli_functions import (
    busco_assess_completeness,
)
from .hifiasm_p13_cli_functions import (
    hifiasm_assemble,
)
from .megahit_p15_cli_functions import (
    megahit_contigs_to_fastg,
    megahit_local_assembly,
    megahit_read_statistics,
)
from .merqury_p10_functions import (
    merqury_assign_reads,
    merqury_bed_to_bigwig,
    merqury_best_k,
    merqury_block_n_stats,
    merqury_combine_kmer_dbs,
    merqury_exclude_reads,
    merqury_false_duplications,
    merqury_false_duplications_track,
    merqury_filter_solid_kmers,
    merqury_hamming_error,
    merqury_hapmer_counts,
    merqury_hapmer_spectra,
    merqury_hapmers,
    merqury_kmer_histogram_table,
    merqury_kmer_multiplicity_track,
    merqury_mash_qv,
    merqury_per_sequence_qv,
    merqury_phase_blocks,
    merqury_plot_blob,
    merqury_plot_block_n,
    merqury_plot_spectra_cn,
    merqury_split_10x_fastq,
)
from .quast_p15_cli_functions import (
    metaquast_assess_metagenome,
    quast_assess_assembly,
    quast_map_reads,
    quast_predict_genes,
)
from .shasta_p29_cli_functions import (
    shasta_assemble,
    shasta_read_assembly_output,
    shasta_read_run_table,
)
from .spades_p14_cli_functions import (
    binspreader_refine_bins,
    pathracer_align_hmm_to_sequences,
    pathracer_search_graph,
    spades_build_assembly_graph,
    spades_correct_reads,
    spades_estimate_kmer_count,
    spades_filter_reads_by_kmer_coverage,
    spades_map_reads_to_graph,
    spades_simplify_assembly_graph,
    spades_split_assembly_graph,
    spaligner_align_long_reads,
    splitter_resolve_repeats,
)
from .spades_p40_cli_functions import (
    spades_assemble,
)

__all__ = [
    *__all__,
    "binspreader_refine_bins",
    "busco_assess_completeness",
    "hifiasm_assemble",
    "megahit_contigs_to_fastg",
    "megahit_local_assembly",
    "megahit_read_statistics",
    "merqury_assign_reads",
    "merqury_bed_to_bigwig",
    "merqury_best_k",
    "merqury_block_n_stats",
    "merqury_combine_kmer_dbs",
    "merqury_exclude_reads",
    "merqury_false_duplications",
    "merqury_false_duplications_track",
    "merqury_filter_solid_kmers",
    "merqury_hamming_error",
    "merqury_hapmer_counts",
    "merqury_hapmer_spectra",
    "merqury_hapmers",
    "merqury_kmer_histogram_table",
    "merqury_kmer_multiplicity_track",
    "merqury_mash_qv",
    "merqury_per_sequence_qv",
    "merqury_phase_blocks",
    "merqury_plot_blob",
    "merqury_plot_block_n",
    "merqury_plot_spectra_cn",
    "merqury_split_10x_fastq",
    "metaquast_assess_metagenome",
    "pathracer_align_hmm_to_sequences",
    "pathracer_search_graph",
    "quast_assess_assembly",
    "quast_map_reads",
    "quast_predict_genes",
    "shasta_assemble",
    "shasta_read_assembly_output",
    "shasta_read_run_table",
    "spades_assemble",
    "spades_build_assembly_graph",
    "spades_correct_reads",
    "spades_estimate_kmer_count",
    "spades_filter_reads_by_kmer_coverage",
    "spades_map_reads_to_graph",
    "spades_simplify_assembly_graph",
    "spades_split_assembly_graph",
    "spaligner_align_long_reads",
    "splitter_resolve_repeats",
]

_register_public_verbose(globals())
