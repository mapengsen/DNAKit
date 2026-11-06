"""IQ-TREE 3 command-line facades (pack P11): fixed argument lists, no shell, bounded execution.

Each public function runs ``iqtree3`` (3.1.4 is the tested version) for one analysis. The user installs IQ-TREE
(``conda install -c bioconda iqtree``, or a release binary from https://github.com/iqtree/iqtree3/releases) and puts
``iqtree3`` on PATH or passes ``iqtree_path``. Every run writes into a fresh ``output_directory`` (never overwritten):
IQ-TREE's own files are ``<prefix>.iqtree`` (report), ``<prefix>.log`` and the analysis-specific files, and
``tool.log`` holds the standard output and error of the program. Results are JSON-serialisable dictionaries with the
status, version, command, ``outputs`` (every file written, keyed by its extension after the prefix), the main numbers
parsed from the report and a bounded log tail. Optional parameters default to ``None`` and are then not passed, so
IQ-TREE's own defaults apply (without ``model`` IQ-TREE runs ModelFinder first, as it does on the command line).
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p11_cli as _p
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_TESTED_VERSION = "3.1.4"
_CITATION = "https://doi.org/10.1093/molbev/msaa015"
_HINT = (
    "Install IQ-TREE 3 (for example `conda install -c bioconda iqtree`, or a release binary from "
    "https://github.com/iqtree/iqtree3/releases), then put iqtree3 on PATH or pass iqtree_path."
)
_SEQUENCE_TYPES = ("BIN", "DNA", "AA", "NT2AA", "CODON", "MORPH")
_SEARCHES = ("TESTONLY", "TEST", "MF", "MFP")
_MERGE_ALGORITHMS = ("greedy", "rcluster", "rclusterf")
_CRITERIA = ("AIC", "AICc", "BIC", "mAIC")
_PARTITION_OPTIONS = {"proportional": "-p", "equal": "-q", "unlinked": "-Q"}
_SUPPORT_METHODS = {
    "ufboot": "-B",
    "ufjack": "-J",
    "boot": "-b",
    "jack": "-j",
    "alrt": "--alrt",
    "abayes": "--abayes",
    "lbp": "--lbp",
    "sprta": "--sprta",
}
_NO_REPLICATES = ("abayes", "sprta")
_MODEL = re.compile(r"^[A-Za-z0-9_.+*,{}()\[\]/:=-]+$")
_FLOAT = r"(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _model(value: str | None, name: str = "model") -> str | None:
    if value is None:
        return None
    text = _cli.text(value, name)
    if _MODEL.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} contains characters that are not valid in an IQ-TREE model string: {text!r}"
        )
    return text


def _common(
    arguments: list[str],
    *,
    sequence_type: str | None,
    threads: int | str,
    seed: int | None,
    model: str | None = None,
    model_flag: bool = True,
) -> None:
    """Options shared by all analyses: sequence type, model, threads and random seed."""
    if sequence_type is not None:
        arguments.extend(
            ("--seqtype", _cli.choice(sequence_type, _SEQUENCE_TYPES, "sequence_type"))
        )
    if model_flag and model is not None:
        arguments.extend(("-m", model))
    if threads == "AUTO":
        arguments.extend(("-T", "AUTO"))
    else:
        arguments.extend(("-T", str(_cli.integer(threads, "threads", minimum=1, maximum=1024))))
    _p.opt_int(arguments, "--seed", seed, "seed", minimum=0)


def _iqtree(
    operation: str,
    arguments: Sequence[str],
    output_directory: str | Path,
    prefix: str,
    *,
    iqtree_path: str | None,
    timeout_seconds: float,
    required: Sequence[str] = ("iqtree",),
    **extra: object,
) -> dict[str, Any]:
    """Run ``iqtree3 <arguments> --prefix <output_directory>/<prefix>`` and collect the files it wrote."""
    executable = _cli.program("iqtree3", explicit=iqtree_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    stem = directory / prefix
    command = [*arguments, "--prefix", str(stem), "--quiet"]
    _p.run_checked(operation, executable, command, directory, timeout_seconds=timeout_seconds)
    outputs: dict[str, Path] = {}
    for path in sorted(directory.glob(f"{prefix}.*")):
        if path.is_file():
            outputs[path.name[len(prefix) + 1 :]] = path
    missing = [suffix for suffix in required if suffix not in outputs]
    if missing:
        raise BackendExecutionError(
            f"{operation} finished but IQ-TREE did not write: {', '.join(f'{prefix}.{item}' for item in missing)}.",
            code="P11_OUTPUT_MISSING",
            context={"operation": operation, "log_tail": _p.tail(directory / _p.LOG_NAME)},
        )
    report = _report(outputs["iqtree"]) if "iqtree" in outputs else {}
    return _p.result(
        operation,
        "iqtree3",
        _p.version_of(executable),
        ["iqtree3", *command],
        directory,
        outputs,
        tool="IQ-TREE",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        prefix=str(stem),
        **report,
        **extra,
    )


def _report(path: Path) -> dict[str, Any]:
    """The main numbers of an IQ-TREE ``.iqtree`` report (keys are present only when the report has them)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    found: dict[str, Any] = {}
    best = re.search(r"Best-fit model according to (\w+): (\S+)", text)
    if best:
        found["best_model"] = best.group(2)
        found["selection_criterion"] = best.group(1)
    used = re.search(r"^Model of substitution: (.+)$", text, re.MULTILINE)
    if used:
        found["model"] = used.group(1).strip()
    likelihood = re.search(r"Log-likelihood of the tree: " + _FLOAT, text)
    if likelihood:
        found["log_likelihood"] = float(likelihood.group(1))
    length = re.search(r"Total tree length \(sum of branch lengths\): " + _FLOAT, text)
    if length:
        found["tree_length"] = float(length.group(1))
    return found


def _alignment(alignment: str | Path) -> str:
    return _p.file_path(alignment, "alignment")


def _prefix(prefix: str) -> str:
    return _p.file_stem(prefix, "prefix")


def _replicates(value: int, name: str, minimum: int = 1) -> str:
    return str(_cli.integer(value, name, minimum=minimum, maximum=10_000_000))


# ---------------------------------------------------------------------------
# Likelihood and model selection
# ---------------------------------------------------------------------------


def evaluate_tree_likelihood(
    alignment: str | Path,
    tree: str | Path,
    output_directory: str | Path,
    *,
    model: str | None = None,
    optimise_branch_lengths: bool = True,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate a fixed tree topology on an alignment: its log-likelihood, with or without branch-length optimisation.

    With ``optimise_branch_lengths=True`` (IQ-TREE ``--tree-fix``) the topology is kept, and branch lengths and model
    parameters are optimised; the optimised tree is ``outputs['treefile']``. With ``False`` (``--show-lh``) the
    log-likelihood is computed on the branch lengths of ``tree`` without optimising them. ``model`` is an IQ-TREE model
    string such as ``GTR+F+G4``; without it IQ-TREE runs ModelFinder first. Result keys: ``log_likelihood``, ``model``,
    ``tree_length``.
    """
    arguments = ["-s", _alignment(alignment), "-t", _p.file_path(tree, "tree")]
    arguments.append(
        "--tree-fix" if _p.flag(optimise_branch_lengths, "optimise_branch_lengths") else "--show-lh"
    )
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    return _iqtree(
        "evaluate_tree_likelihood",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("iqtree", "treefile"),
    )


def select_substitution_model(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    search: str = "MFP",
    partitions: str | Path | None = None,
    partition_mode: str = "proportional",
    merge: bool = False,
    merge_algorithm: str | None = None,
    criterion: str | None = None,
    model_set: str | None = None,
    model_subset: str | None = None,
    rate_models: str | None = None,
    max_rate_categories: int | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the best-fit substitution model (ModelFinder) and optionally merge partitions (PartitionFinder-style).

    ``search`` is IQ-TREE's ``-m`` keyword: ``TESTONLY`` (standard model selection only), ``TEST`` (standard selection,
    then tree inference), ``MF`` (extended selection with FreeRate models, no tree search) or ``MFP`` (extended
    selection, then tree inference; IQ-TREE's default). With ``partitions`` (a NEXUS or RAxML partition file;
    ``partition_mode`` ``proportional`` ``-p``, ``equal`` ``-q`` or ``unlinked`` ``-Q``) a model is chosen per partition,
    and with ``merge=True`` the partitions are also merged to increase the model fit (``--merge``;
    ``merge_algorithm`` ``greedy``, ``rcluster`` or ``rclusterf``, IQ-TREE's default). ``criterion`` is ``AIC``,
    ``AICc``, ``BIC`` (default) or ``mAIC``; ``model_set`` (``--mset``), ``model_subset`` (``--msub``), ``rate_models``
    (``--mrate``, for example ``E,I,G,I+G,R``) and ``max_rate_categories`` (``--cmax``) restrict the search. Result keys:
    ``best_model``, ``selection_criterion``, ``log_likelihood`` and ``model_table`` (the first 20 rows of the ranked
    model list: ``model``, ``log_likelihood``, ``aic``, ``aicc``, ``bic``); with merging ``best_scheme`` holds the path
    of the merged partition scheme.
    """
    keyword = _cli.choice(search, _SEARCHES, "search")
    merging = _p.flag(merge, "merge")
    if merging and partitions is None:
        raise ConfigurationError("merge needs a partition file (partitions).")
    if not merging and merge_algorithm is not None:
        raise ConfigurationError("merge_algorithm is only used with merge=True.")
    model_string = (
        {"TESTONLY": "TESTMERGEONLY", "TEST": "TESTMERGE", "MF": "MF+MERGE", "MFP": "MFP+MERGE"}[
            keyword
        ]
        if merging
        else keyword
    )
    arguments = ["-s", _alignment(alignment)]
    if partitions is not None:
        mode = _cli.choice(partition_mode, tuple(_PARTITION_OPTIONS), "partition_mode")
        arguments.extend((_PARTITION_OPTIONS[mode], _p.file_path(partitions, "partitions")))
    arguments.extend(("-m", model_string))
    if merging:
        arguments.append("--merge")
        if merge_algorithm is not None:
            arguments.append(_cli.choice(merge_algorithm, _MERGE_ALGORITHMS, "merge_algorithm"))
    if criterion is not None:
        arguments.extend(("--merit", _cli.choice(criterion, _CRITERIA, "criterion")))
    if model_set is not None:
        arguments.extend(("--mset", _model(model_set, "model_set") or ""))
    if model_subset is not None:
        arguments.extend(
            (
                "--msub",
                _cli.choice(
                    model_subset,
                    ("nuclear", "mitochondrial", "chloroplast", "viral"),
                    "model_subset",
                ),
            )
        )
    if rate_models is not None:
        arguments.extend(("--mrate", _model(rate_models, "rate_models") or ""))
    _p.opt_int(
        arguments, "--cmax", max_rate_categories, "max_rate_categories", minimum=2, maximum=50
    )
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model_flag=False)
    result = _iqtree(
        "select_substitution_model",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
    )
    outputs = result["outputs"]
    rows: list[dict[str, Any]] = []
    if "iqtree" in outputs:
        text = Path(outputs["iqtree"]).read_text(encoding="utf-8", errors="replace")
        block = re.search(
            r"List of models sorted by \w+ scores:\s+Model\s+LogL[^\n]*\n(.*?)\n\n", text, re.DOTALL
        )
        if block:
            for line in block.group(1).splitlines()[:20]:
                fields = line.split()
                if len(fields) >= 11:
                    rows.append(
                        {
                            "model": fields[0],
                            "log_likelihood": float(fields[1]),
                            "aic": float(fields[2]),
                            "aicc": float(fields[5]),
                            "bic": float(fields[8]),
                        }
                    )
    result["model_table"] = rows
    if "best_scheme" in outputs:
        result["best_scheme"] = outputs["best_scheme"]
    return result


# ---------------------------------------------------------------------------
# Tree inference and branch support
# ---------------------------------------------------------------------------


def infer_pathogen_tree(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    model: str | None = None,
    force_cmaple: bool = False,
    alrt_replicates: int | None = None,
    sprta: bool = False,
    sprta_zero_length_branches: bool = False,
    sprta_alternative_placements: bool = False,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer a tree for low-divergence genomes (outbreak or pathogen data) with the CMAPLE algorithm (``--pathogen``).

    IQ-TREE applies CMAPLE when the sequence divergence is low and its standard search otherwise;
    ``force_cmaple=True`` (``--pathogen-force``) applies CMAPLE regardless of divergence. ``alrt_replicates`` adds SH-aLRT
    branch supports (``--alrt``), ``sprta`` adds SPRTA supports (``--sprta``), with ``sprta_zero_length_branches``
    (``--sprta-zero-branch``) and ``sprta_alternative_placements`` (``--sprta-other-places``). The tree is
    ``outputs['treefile']``.
    """
    arguments = [
        "-s",
        _alignment(alignment),
        "--pathogen-force" if _p.flag(force_cmaple, "force_cmaple") else "--pathogen",
    ]
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    if alrt_replicates is not None:
        arguments.extend(("--alrt", _replicates(alrt_replicates, "alrt_replicates", 0)))
    _p.switch(arguments, "--sprta", sprta, "sprta")
    _p.switch(
        arguments, "--sprta-zero-branch", sprta_zero_length_branches, "sprta_zero_length_branches"
    )
    _p.switch(
        arguments,
        "--sprta-other-places",
        sprta_alternative_placements,
        "sprta_alternative_placements",
    )
    return _iqtree(
        "infer_pathogen_tree",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("treefile",),
    )


def compute_branch_support(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    methods: Sequence[str] = ("ufboot",),
    replicates: int | Mapping[str, int] = 1000,
    tree: str | Path | None = None,
    model: str | None = None,
    jackknife_proportion: float | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer a maximum-likelihood tree and compute branch supports with one or several methods.

    ``methods`` may combine ``ufboot`` (ultrafast bootstrap, ``-B``), ``ufjack`` (ultrafast jackknife, ``-J``),
    ``boot`` (non-parametric bootstrap, ``-b``), ``jack`` (non-parametric jackknife, ``-j``), ``alrt`` (SH-aLRT,
    ``--alrt``), ``abayes`` (approximate Bayes test, ``--abayes``), ``lbp`` (fast local bootstrap probabilities,
    ``--lbp``) and ``sprta`` (SPRTA branch support, ``--sprta``; this selects IQ-TREE's CMAPLE tree search). IQ-TREE
    requires at least 1000 replicates for ``ufboot`` and ``ufjack``. ``replicates`` is one number for all methods or a
    mapping method to number. With ``tree`` the supports are computed on that topology instead of an inferred tree.
    ``jackknife_proportion`` is ``--jack-prop`` (default 0.5). Result keys: ``methods``, ``treefile`` (the ML tree with
    supports), ``consensus_tree`` (the ``contree`` file of the bootstrap/jackknife methods, or None), ``log_likelihood``,
    ``model``.
    """
    if isinstance(methods, str) or not methods:
        raise ConfigurationError("methods must be a non-empty list of method names.")
    chosen = [_cli.choice(method, tuple(_SUPPORT_METHODS), "methods") for method in methods]
    if len(set(chosen)) != len(chosen):
        raise ConfigurationError("methods must not repeat a method.")
    arguments = ["-s", _alignment(alignment)]
    if tree is not None:
        arguments.extend(("-t", _p.file_path(tree, "tree")))
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    for method in chosen:
        flag = _SUPPORT_METHODS[method]
        if method in _NO_REPLICATES:
            arguments.append(flag)
            continue
        count = replicates.get(method) if isinstance(replicates, Mapping) else replicates
        if count is None:
            raise ConfigurationError(f"replicates has no entry for {method}.")
        minimum = 1000 if method in ("ufboot", "ufjack") else 1
        arguments.extend((flag, _replicates(count, f"replicates[{method}]", minimum)))
    _p.opt_num(
        arguments,
        "--jack-prop",
        jackknife_proportion,
        "jackknife_proportion",
        minimum=0.0,
        maximum=1.0,
    )
    result = _iqtree(
        "compute_branch_support",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("treefile",),
        methods=chosen,
    )
    result["consensus_tree"] = result["outputs"].get("contree")
    return result


def compare_tree_topologies(
    alignment: str | Path,
    trees: str | Path,
    output_directory: str | Path,
    *,
    model: str | None = None,
    replicates: int = 10000,
    approximately_unbiased: bool = True,
    weighted: bool = False,
    site_likelihoods: bool = False,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare tree topologies on an alignment with the KH, SH and approximately unbiased (AU) tests.

    ``trees`` is a file with one tree per line. IQ-TREE evaluates the log-likelihood of every tree (``--trees``) and runs
    the tests with ``replicates`` RELL replicates (``--test``); ``approximately_unbiased`` adds the AU test
    (``--test-au``), ``weighted`` the weighted KH and SH tests (``--test-weight``) and ``site_likelihoods`` writes the
    per-site log-likelihoods (``--sitelh``). Result key: ``tests``, one row per tree with ``tree`` (1-based), ``logl``,
    ``delta_l``, and for each test (``bp_rell``, ``p_kh``, ``p_sh``, ``c_elw``, ``p_au``) ``value`` and ``passed``
    (IQ-TREE's ``+``/``-`` sign: the tree is not rejected at 5%).
    """
    arguments = ["-s", _alignment(alignment), "--trees", _p.file_path(trees, "trees")]
    arguments.extend(("--test", _replicates(replicates, "replicates")))
    _p.switch(arguments, "--test-au", approximately_unbiased, "approximately_unbiased")
    _p.switch(arguments, "--test-weight", weighted, "weighted")
    _p.switch(arguments, "--sitelh", site_likelihoods, "site_likelihoods")
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    result = _iqtree(
        "compare_tree_topologies",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
    )
    text = Path(result["outputs"]["iqtree"]).read_text(encoding="utf-8", errors="replace")
    block = re.search(r"USER TREES\n-+\n.*?\n(Tree\s+logL[^\n]*)\n-+\n(.*?)\n\n", text, re.DOTALL)
    rows: list[dict[str, Any]] = []
    if block:
        columns = [name for name in block.group(1).split() if name not in ("Tree",)]
        for line in block.group(2).splitlines():
            fields = line.split()
            row: dict[str, Any] = {"tree": int(fields[0])}
            index = 1
            for column in columns:
                key = column.lower().replace("-", "_")
                if key in ("logl", "deltal"):
                    row["logl" if key == "logl" else "delta_l"] = float(fields[index])
                    index += 1
                else:
                    row[key] = {"value": float(fields[index]), "passed": fields[index + 1] == "+"}
                    index += 2
            rows.append(row)
    result["tests"] = rows
    return result


def likelihood_mapping(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    quartets: int = 10000,
    clusters: str | Path | None = None,
    model: str | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assess the phylogenetic signal of an alignment by likelihood mapping of random quartets (``--lmap``).

    ``quartets`` is the number of quartets drawn (IQ-TREE uses all unique quartets when fewer exist); ``clusters`` is an
    optional NEXUS file defining sequence clusters (``--lmclust``). The maps are written as ``lmap.svg`` and
    ``lmap.eps`` in ``outputs``. Result keys: ``quartets`` and ``resolution`` (``fully_resolved``, ``partly_resolved`` and
    ``unresolved`` each with ``count`` and ``percent``).
    """
    arguments = ["-s", _alignment(alignment), "--lmap", _replicates(quartets, "quartets")]
    _p.opt_file(arguments, "--lmclust", clusters, "clusters")
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    result = _iqtree(
        "likelihood_mapping",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
    )
    text = Path(result["outputs"]["iqtree"]).read_text(encoding="utf-8", errors="replace")
    count = re.search(r"Number of quartets: (\d+)", text)
    result["quartets"] = int(count.group(1)) if count else None
    resolution: dict[str, Any] = {}
    for label, key in (
        ("fully resolved", "fully_resolved"),
        ("partly resolved", "partly_resolved"),
        ("unresolved", "unresolved"),
    ):
        match = re.search(rf"Number of {label}\s+quartets[^:]*: (\d+) \(=([\d.]+)%\)", text)
        if match:
            resolution[key] = {"count": int(match.group(1)), "percent": float(match.group(2))}
    result["resolution"] = resolution
    return result


def compute_site_concordance_factors(
    alignment: str | Path,
    tree: str | Path,
    output_directory: str | Path,
    *,
    quartets: int = 100,
    method: str = "likelihood",
    partitions: str | Path | None = None,
    model: str | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute site concordance factors (sCF) of the branches of a reference tree by sampling quartets of sequences.

    ``method`` is ``likelihood`` (``--scfl``, recommended by the IQ-TREE authors) or ``parsimony`` (``--scf``);
    ``quartets`` is the number of quartets sampled around each branch. ``partitions`` is an optional partition file for
    partitioned alignments. The annotated tree is ``outputs['cf.tree']``. Result key: ``branches``, one row per branch
    with ``id``, ``scf``, ``scf_n``, ``sdf1``, ``sdf1_n``, ``sdf2``, ``sdf2_n``, ``sn``, ``label`` and ``length``.
    """
    flag = {"likelihood": "--scfl", "parsimony": "--scf"}[
        _cli.choice(method, ("likelihood", "parsimony"), "method")
    ]
    arguments = [
        "-s",
        _alignment(alignment),
        "-t",
        _p.file_path(tree, "tree"),
        flag,
        _replicates(quartets, "quartets"),
    ]
    _p.opt_file(arguments, "-p", partitions, "partitions")
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    result = _iqtree(
        "compute_site_concordance_factors",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("cf.stat", "cf.tree"),
    )
    names = ("id", "scf", "scf_n", "sdf1", "sdf1_n", "sdf2", "sdf2_n", "sn", "label", "length")
    rows, _ = _p.table(Path(result["outputs"]["cf.stat"]))
    branches = []
    for row in rows[1:]:
        padded = row + [""] * (len(names) - len(row))
        branches.append(
            {
                name: _p.number(value) if value != "" else None
                for name, value in zip(names, padded, strict=False)
            }
        )
    result["branches"] = branches
    return result


# ---------------------------------------------------------------------------
# Ancestral states, symmetry, alignment statistics
# ---------------------------------------------------------------------------


def reconstruct_ancestral_sequences(
    alignment: str | Path,
    tree: str | Path,
    output_directory: str | Path,
    *,
    model: str | None = None,
    min_probability: float | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reconstruct the ancestral sequence of every internal node by empirical Bayes (``--ancestral``).

    The tree topology is fixed to ``tree``. ``min_probability`` is ``--asr-min`` (minimum posterior probability for a state
    to be called; below it IQ-TREE writes ``?``). The per-site posterior probabilities are in ``outputs['state']``.
    Result keys: ``nodes`` (internal node names), ``sites`` and ``sequences`` (node name to its most likely sequence; None
    when the table is larger than 2 million states).
    """
    arguments = ["-s", _alignment(alignment), "-t", _p.file_path(tree, "tree"), "--ancestral"]
    _p.opt_num(arguments, "--asr-min", min_probability, "min_probability", minimum=0.0, maximum=1.0)
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    arguments.append("--tree-fix")
    result = _iqtree(
        "reconstruct_ancestral_sequences",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("iqtree", "state"),
    )
    states: dict[str, list[str]] = {}
    total = 0
    sites = 0
    with Path(result["outputs"]["state"]).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("Node\t"):
                continue
            fields = line.split("\t", 3)
            if len(fields) < 3:
                continue
            total += 1
            sites = max(sites, int(fields[1]))
            if total <= 2_000_000:
                states.setdefault(fields[0], []).append(fields[2])
    result["nodes"] = list(states) if total <= 2_000_000 else None
    result["sites"] = sites
    result["sequences"] = (
        {node: "".join(values) for node, values in states.items()} if total <= 2_000_000 else None
    )
    return result


def run_symmetry_tests(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    only_test: bool = True,
    remove_partitions: str | None = None,
    test_type: str | None = None,
    p_value: float | None = None,
    keep_zero: bool = False,
    partitions: str | Path | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test the stationarity, reversibility and homogeneity assumptions with the matched-pair tests of symmetry.

    IQ-TREE runs three tests on every sequence pair (symmetry, marginal symmetry, internal symmetry; Naser-Khdour et al.
    2019). ``only_test=True`` runs them and stops (``--symtest-only``); with ``False`` the tests are followed by the
    tree analysis (``--symtest``). ``remove_partitions`` ``bad`` or ``good`` removes the partitions that fail or pass the
    tests (``--symtest-remove-bad`` / ``--symtest-remove-good``; needs ``partitions``), judged by ``test_type``
    ``MAR`` or ``INT`` and the cut-off ``p_value`` (default 0.05); ``keep_zero`` keeps NA results (``--symtest-keep-zero``).
    Result key: ``tests`` (one row per partition from ``outputs['symtest.csv']``: ``name``, ``sym_sig``, ``sym_non``,
    ``sym_pval``, ``mar_sig``, ``mar_non``, ``mar_pval``, ``int_sig``, ``int_non``, ``int_pval``).
    """
    arguments = ["-s", _alignment(alignment)]
    _p.opt_file(arguments, "-p", partitions, "partitions")
    if remove_partitions is not None:
        if partitions is None:
            raise ConfigurationError("remove_partitions needs a partition file (partitions).")
        arguments.append(
            "--symtest-remove-"
            + _cli.choice(remove_partitions, ("bad", "good"), "remove_partitions")
        )
    else:
        arguments.append("--symtest-only" if _p.flag(only_test, "only_test") else "--symtest")
    if test_type is not None:
        arguments.extend(("--symtest-type", _cli.choice(test_type, ("MAR", "INT"), "test_type")))
    _p.opt_num(arguments, "--symtest-pval", p_value, "p_value", minimum=0.0, maximum=1.0)
    _p.switch(arguments, "--symtest-keep-zero", keep_zero, "keep_zero")
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed)
    result = _iqtree(
        "run_symmetry_tests",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("symtest.csv",),
    )
    names = (
        "name",
        "sym_sig",
        "sym_non",
        "sym_pval",
        "mar_sig",
        "mar_non",
        "mar_pval",
        "int_sig",
        "int_non",
        "int_pval",
    )
    tests: list[dict[str, Any]] = []
    lines = [
        line
        for line in Path(result["outputs"]["symtest.csv"]).read_text(encoding="utf-8").splitlines()
        if line and line[0] != "#"
    ]
    for line in lines[1:]:
        fields = line.split(",")
        tests.append({name: _p.number(value) for name, value in zip(names, fields, strict=False)})
    result["tests"] = tests
    return result


def alignment_site_statistics(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    model: str = "JC",
    sequence_type: str | None = None,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify every alignment site as constant, informative or uninformative (``-alninfo``).

    IQ-TREE reads the alignment, writes ``<prefix>.alninfo`` and stops before any tree search (``-n 0`` with the cheap
    ``model``, default ``JC``). Site classes: ``I`` parsimony-informative, ``C`` constant, ``c`` constant with ambiguous
    characters, ``U`` uninformative but not constant, ``-`` all gaps. Result keys: ``sites`` and ``counts`` (class
    letter to number of sites).
    """
    arguments = ["-s", _alignment(alignment), "-alninfo", "-n", "0"]
    _common(arguments, sequence_type=sequence_type, threads=1, seed=seed, model=_model(model))
    result = _iqtree(
        "alignment_site_statistics",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("alninfo",),
    )
    counts: dict[str, int] = {}
    for line in Path(result["outputs"]["alninfo"]).read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0].isdigit():
            counts[fields[1]] = counts.get(fields[1], 0) + 1
    result["counts"] = counts
    result["sites"] = sum(counts.values())
    return result


# ---------------------------------------------------------------------------
# Simulation, dating, PoMo
# ---------------------------------------------------------------------------


def simulate_alignment(
    output_directory: str | Path,
    *,
    tree: str | Path | None = None,
    random_tree_taxa: int | None = None,
    random_tree_model: str = "yh",
    model: str = "JC",
    length: int = 1000,
    alignments: int = 1,
    output_format: str = "fasta",
    template_alignment: str | Path | None = None,
    branch_scale: float | None = None,
    partitions: str | Path | None = None,
    insertion_deletion_rates: Sequence[float] | None = None,
    sequence_type: str | None = None,
    seed: int | None = None,
    prefix: str = "simulated",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate sequence alignments along a tree with AliSim (IQ-TREE ``--alisim``).

    Give a ``tree`` file, or ``random_tree_taxa`` to generate a random tree (``-t RANDOM{<random_tree_model>/<taxa>}``;
    ``yh`` Yule-Harding, ``u`` uniform, ``bd`` birth-death, ...). ``model`` is an IQ-TREE model string (``HKY+F{0.3/0.2/0.2/0.3}+G4``);
    ``length`` is the length of the root sequence and ``alignments`` the number of datasets. ``output_format`` is
    ``fasta`` or ``phylip``. ``template_alignment`` (``-s``) makes AliSim copy the gap pattern and mimic the site
    frequencies and rates of a real alignment; ``partitions`` (``-p``) simulates a partitioned alignment;
    ``branch_scale`` multiplies all branch lengths; ``insertion_deletion_rates`` is the pair ``(insertion, deletion)`` of
    the indel model (``--indel``). Set ``seed`` for reproducible output. Result keys: ``alignment_files`` (paths) and
    ``sequences`` and ``sites`` of the first alignment.
    """
    if (tree is None) == (random_tree_taxa is None):
        raise ConfigurationError("Give exactly one of tree and random_tree_taxa.")
    stem = _prefix(prefix)
    arguments = ["--alisim", stem]
    if tree is not None:
        arguments.extend(("-t", _p.file_path(tree, "tree")))
    else:
        generator = _p.file_stem(random_tree_model, "random_tree_model")
        taxa = _cli.integer(random_tree_taxa, "random_tree_taxa", minimum=3, maximum=1_000_000)
        arguments.extend(("-t", f"RANDOM{{{generator}/{taxa}}}"))
    arguments.extend(("-m", _model(model) or "JC"))
    arguments.extend(
        ("--length", str(_cli.integer(length, "length", minimum=1, maximum=1_000_000_000)))
    )
    arguments.extend(
        (
            "--num-alignments",
            str(_cli.integer(alignments, "alignments", minimum=1, maximum=100_000)),
        )
    )
    arguments.extend(
        (
            "-af",
            {"fasta": "fasta", "phylip": "phy"}[
                _cli.choice(output_format, ("fasta", "phylip"), "output_format")
            ],
        )
    )
    _p.opt_file(arguments, "-s", template_alignment, "template_alignment")
    _p.opt_file(arguments, "-p", partitions, "partitions")
    _p.opt_num(arguments, "--branch-scale", branch_scale, "branch_scale", minimum=0.0)
    if insertion_deletion_rates is not None:
        if len(insertion_deletion_rates) != 2:
            raise ConfigurationError(
                "insertion_deletion_rates must be (insertion rate, deletion rate)."
            )
        rates = [
            format(_cli.number(item, "insertion_deletion_rates", minimum=0.0), ".12g")
            for item in insertion_deletion_rates
        ]
        arguments.extend(("--indel", ",".join(rates)))
    if sequence_type is not None:
        arguments.extend(
            ("--seqtype", _cli.choice(sequence_type, _SEQUENCE_TYPES, "sequence_type"))
        )
    _p.opt_int(arguments, "--seed", seed, "seed", minimum=0)
    executable = _cli.program("iqtree3", explicit=iqtree_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    arguments[1] = str(directory / stem)
    # AliSim writes <prefix>.treefile for a random tree and then refuses to continue unless --redo is given;
    # the output directory is fresh, so nothing of the user's can be overwritten.
    command = [*arguments, "--redo", "--quiet"]
    _p.run_checked(
        "simulate_alignment", executable, command, directory, timeout_seconds=timeout_seconds
    )
    files = sorted(
        path
        for path in directory.iterdir()
        if path.is_file()
        and path.name.startswith(stem)
        and path.suffix in (".fa", ".phy", ".fasta")
    )
    if not files:
        raise BackendExecutionError(
            "simulate_alignment finished but AliSim wrote no alignment.",
            code="P11_OUTPUT_MISSING",
            context={"log_tail": _p.tail(directory / _p.LOG_NAME)},
        )
    first = files[0].read_text(encoding="utf-8", errors="replace").splitlines()
    if output_format == "fasta":
        sequences = sum(1 for line in first if line.startswith(">"))
        sites = len(first[1]) if len(first) > 1 else 0
    else:
        header = first[0].split() if first else []
        sequences, sites = (int(header[0]), int(header[1])) if len(header) >= 2 else (0, 0)
    outputs = {
        path.name: path
        for path in sorted(directory.iterdir())
        if path.is_file() and path.name != _p.LOG_NAME
    }
    return _p.result(
        "simulate_alignment",
        "iqtree3",
        _p.version_of(executable),
        ["iqtree3", *command],
        directory,
        outputs,
        tool="IQ-TREE",
        tested_version=_TESTED_VERSION,
        citation="https://doi.org/10.1093/sysbio/syac027",
        alignment_files=[str(path) for path in files],
        sequences=sequences,
        sites=sites,
    )


def date_phylogeny(
    alignment: str | Path,
    dates: str | Path,
    output_directory: str | Path,
    *,
    tree: str | Path | None = None,
    model: str | None = None,
    tip_date: str | None = None,
    root_date: str | None = None,
    confidence_replicates: int | None = None,
    clock_sd: float | None = None,
    exclude_outgroup: bool = False,
    outlier_zscore: float | None = None,
    sequence_type: str | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate a time-scaled tree and divergence times from an alignment and sampling dates (least-squares dating, LSD2).

    ``dates`` is a file with one taxon name and date per line, or the word ``TAXNAME`` to read the date after the last
    ``|`` of each taxon name (``--date``). ``tip_date`` and ``root_date`` set the date of all tips or of the root
    (real number or ``YYYY-MM-DD``); ``confidence_replicates`` computes confidence intervals with that many resamples
    (``--date-ci``); ``clock_sd`` is the standard deviation of the log-normal relaxed clock (``--clock-sd``,
    default 0.2); ``exclude_outgroup`` is ``--date-no-outgroup`` and ``outlier_zscore`` removes tips whose z-score
    exceeds it (``--date-outlier``). With ``tree`` that topology is dated, otherwise IQ-TREE infers the ML tree first.
    The time tree is ``outputs['timetree.nwk']`` (also ``.nex`` and the LSD2 report ``.lsd``).
    """
    arguments = ["-s", _alignment(alignment)]
    if tree is not None:
        arguments.extend(("-t", _p.file_path(tree, "tree")))
    if isinstance(dates, str) and dates == "TAXNAME":
        arguments.extend(("--date", "TAXNAME"))
    else:
        arguments.extend(("--date", _p.file_path(dates, "dates")))
    _p.opt_text(arguments, "--date-tip", tip_date, "tip_date")
    _p.opt_text(arguments, "--date-root", root_date, "root_date")
    _p.opt_int(arguments, "--date-ci", confidence_replicates, "confidence_replicates", minimum=1)
    _p.opt_num(arguments, "--clock-sd", clock_sd, "clock_sd", minimum=0.0)
    _p.switch(arguments, "--date-no-outgroup", exclude_outgroup, "exclude_outgroup")
    _p.opt_num(arguments, "--date-outlier", outlier_zscore, "outlier_zscore", minimum=0.0)
    arguments.extend(("--dating", "LSD"))
    _common(arguments, sequence_type=sequence_type, threads=threads, seed=seed, model=_model(model))
    return _iqtree(
        "date_phylogeny",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("timetree.nwk",),
    )


def infer_polymorphism_aware_tree(
    counts: str | Path,
    output_directory: str | Path,
    *,
    substitution_model: str = "HKY",
    population_size: int | None = None,
    sampling: str | None = None,
    gamma_categories: int | None = None,
    tree: str | Path | None = None,
    threads: int | str = 1,
    seed: int | None = None,
    prefix: str = "result",
    iqtree_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer a population or species tree from allele counts with polymorphism-aware phylogenetic models (PoMo).

    ``counts`` is an IQ-TREE counts file (``COUNTSFILE NPOP <n> NSITES <m>``, then ``CHROM POS`` and one column of four
    comma-separated nucleotide counts per population). The model is ``<substitution_model>+P`` with optionally
    ``+N<population_size>`` (virtual population size, IQ-TREE default 9), a sampling method (``sampling`` ``WB``
    weighted binomial, ``WH`` weighted hypergeometric or ``S`` sampled) and ``+G<gamma_categories>``. ``tree`` fixes a
    starting tree (``-t``). The tree is ``outputs['treefile']``; ``model`` is the model IQ-TREE reports.
    """
    model_string = _model(substitution_model, "substitution_model") or "HKY"
    model_string += "+P"
    if population_size is not None:
        model_string += (
            f"+N{_cli.integer(population_size, 'population_size', minimum=2, maximum=1000)}"
        )
    if sampling is not None:
        model_string += "+" + _cli.choice(sampling, ("WB", "WH", "S"), "sampling")
    if gamma_categories is not None:
        model_string += (
            f"+G{_cli.integer(gamma_categories, 'gamma_categories', minimum=2, maximum=32)}"
        )
    arguments = ["-s", _p.file_path(counts, "counts"), "-m", model_string]
    _p.opt_file(arguments, "-t", tree, "tree")
    _common(arguments, sequence_type=None, threads=threads, seed=seed, model_flag=False)
    return _iqtree(
        "infer_polymorphism_aware_tree",
        arguments,
        output_directory,
        _prefix(prefix),
        iqtree_path=iqtree_path,
        timeout_seconds=timeout_seconds,
        required=("treefile",),
        pomo_model=model_string,
    )


__all__ = [
    "alignment_site_statistics",
    "compare_tree_topologies",
    "compute_branch_support",
    "compute_site_concordance_factors",
    "date_phylogeny",
    "evaluate_tree_likelihood",
    "infer_pathogen_tree",
    "infer_polymorphism_aware_tree",
    "likelihood_mapping",
    "reconstruct_ancestral_sequences",
    "run_symmetry_tests",
    "select_substitution_model",
    "simulate_alignment",
]
_register_public_verbose(globals())
