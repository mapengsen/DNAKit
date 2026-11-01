# IQ-TREE 3: model selection, branch support, tree tests, simulation and dating

This group of functions in `dnakit.phylogeny` is a set of command-line facades over [IQ-TREE 3](https://github.com/iqtree/iqtree3). They evaluate a fixed tree (`--tree-fix`, `--show-lh`), select substitution models and merge partitions (ModelFinder, `--merge`), infer trees for low-divergence genomes (CMAPLE, `--pathogen`), compute branch supports (UFBoot, ultrafast jackknife, bootstrap, jackknife, SH-aLRT, approximate Bayes, local bootstrap probabilities, SPRTA), run KH/SH/AU topology tests, likelihood mapping, site concordance factors, ancestral sequence reconstruction, tests of symmetry, alignment site statistics, AliSim simulation, LSD2 dating and polymorphism-aware (PoMo) inference. Related rows of the work list that differ only by a method are one function with a parameter (`compute_branch_support(methods=...)`, `select_substitution_model(search=..., merge=...)`, `evaluate_tree_likelihood(optimise_branch_lengths=...)`, `compare_tree_topologies(approximately_unbiased=...)`, `compute_site_concordance_factors(method=...)`). Ultrafast bootstrap alone is also available as `bootstrap` in the existing `infer_ml_tree`.

IQ-TREE 3.1.4 is the tested version (GPL-2.0; DNAKit starts it as a separate program). Install it with `conda install -c bioconda iqtree` or from the release page and put `iqtree3` on PATH, or pass `iqtree_path`. Every function writes into a new `output_directory` that must not exist, runs the program without a shell, and returns a JSON-serialisable dictionary; `outputs` maps each file extension after the prefix (`treefile`, `iqtree`, `contree`, ...) to its path. Options left at `None` or `False` are not passed, so IQ-TREE's defaults apply; without `model`, IQ-TREE runs ModelFinder first. The examples use `verbose=False`.

## Example data and examples

**Example file:** `tree.nwk`

```text
((A:0.1,B:0.12):0.05,((C:0.2,D:0.1):0.07,(E:0.15,F:0.11):0.06):0.04,(G:0.3,H:0.2):0.05);
```

**Example code:**

```python
from pathlib import Path
import tempfile
from dnakit.phylogeny import iqtree_p11_cli_functions as iq

work = Path(tempfile.mkdtemp())
sim = iq.simulate_alignment(work / "sim", tree="tree.nwk", model="HKY", length=400, seed=7, verbose=False)
print(sim["sequences"], sim["sites"], Path(sim["alignment_files"][0]).name)
stats = iq.alignment_site_statistics(sim["alignment_files"][0], work / "stats", verbose=False)
print(stats["sites"], sorted(stats["counts"]))
```

**Example output:**

```text
8 400 simulated.fa
400 ['C', 'I', 'U']
```

**Example code:**

```python
from pathlib import Path
import tempfile
from dnakit.phylogeny import iqtree_p11_cli_functions as iq

work = Path(tempfile.mkdtemp())
sim = iq.simulate_alignment(work / "sim", tree="tree.nwk", model="HKY", length=400, seed=7, verbose=False)
fixed = iq.evaluate_tree_likelihood(sim["alignment_files"][0], "tree.nwk", work / "fit", model="HKY", verbose=False)
print(fixed["version"], round(fixed["log_likelihood"]) < 0)
asr = iq.reconstruct_ancestral_sequences(sim["alignment_files"][0], "tree.nwk", work / "asr", model="HKY", verbose=False)
print(len(asr["nodes"]), asr["sites"])
```

**Example output:**

```text
3.1.4 True
6 400
```


## 1) `dnakit.phylogeny.iqtree_p11_cli_functions.alignment_site_statistics()`

- **Purpose:** Classify every alignment site as constant, informative or uninformative (``-alninfo``).
- **Method:** IQ-TREE reads the alignment, writes ``<prefix>.alninfo`` and stops before any tree search (``-n 0`` with the cheap ``model``, default ``JC``). Site classes: ``I`` parsimony-informative, ``C`` constant, ``c`` constant with ambiguous characters, ``U`` uninformative but not constant, ``-`` all gaps.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.alignment_site_statistics(alignment: 'str | Path', output_directory: 'str | Path', *, model: 'str' = 'JC', sequence_type: 'str | None' = None, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``sites`` and ``counts`` (class letter to number of sites).

## 2) `dnakit.phylogeny.iqtree_p11_cli_functions.compare_tree_topologies()`

- **Purpose:** Compare tree topologies on an alignment with the KH, SH and approximately unbiased (AU) tests.
- **Method:** ``trees`` is a file with one tree per line. IQ-TREE evaluates the log-likelihood of every tree (``--trees``) and runs the tests with ``replicates`` RELL replicates (``--test``); ``approximately_unbiased`` adds the AU test (``--test-au``), ``weighted`` the weighted KH and SH tests (``--test-weight``) and ``site_likelihoods`` writes the per-site log-likelihoods (``--sitelh``).
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.compare_tree_topologies(alignment: 'str | Path', trees: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, replicates: 'int' = 10000, approximately_unbiased: 'bool' = True, weighted: 'bool' = False, site_likelihoods: 'bool' = False, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `trees`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``tests``, one row per tree with ``tree`` (1-based), ``logl``, ``delta_l``, and for each test (``bp_rell``, ``p_kh``, ``p_sh``, ``c_elw``, ``p_au``) ``value`` and ``passed`` (IQ-TREE's ``+``/``-`` sign: the tree is not rejected at 5%).

## 3) `dnakit.phylogeny.iqtree_p11_cli_functions.compute_branch_support()`

- **Purpose:** Infer a maximum-likelihood tree and compute branch supports with one or several methods.
- **Method:** ``methods`` may combine ``ufboot`` (ultrafast bootstrap, ``-B``), ``ufjack`` (ultrafast jackknife, ``-J``), ``boot`` (non-parametric bootstrap, ``-b``), ``jack`` (non-parametric jackknife, ``-j``), ``alrt`` (SH-aLRT, ``--alrt``), ``abayes`` (approximate Bayes test, ``--abayes``), ``lbp`` (fast local bootstrap probabilities, ``--lbp``) and ``sprta`` (SPRTA branch support, ``--sprta``; this selects IQ-TREE's CMAPLE tree search). IQ-TREE requires at least 1000 replicates for ``ufboot`` and ``ufjack``. ``replicates`` is one number for all methods or a mapping method to number. With ``tree`` the supports are computed on that topology instead of an inferred tree. ``jackknife_proportion`` is ``--jack-prop`` (default 0.5).
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.compute_branch_support(alignment: 'str | Path', output_directory: 'str | Path', *, methods: 'Sequence[str]' = ('ufboot',), replicates: 'int | Mapping[str, int]' = 1000, tree: 'str | Path | None' = None, model: 'str | None' = None, jackknife_proportion: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``methods``, ``treefile`` (the ML tree with supports), ``consensus_tree`` (the ``contree`` file of the bootstrap/jackknife methods, or None), ``log_likelihood``, ``model``.

## 4) `dnakit.phylogeny.iqtree_p11_cli_functions.compute_site_concordance_factors()`

- **Purpose:** Compute site concordance factors (sCF) of the branches of a reference tree by sampling quartets of sequences.
- **Method:** ``method`` is ``likelihood`` (``--scfl``, recommended by the IQ-TREE authors) or ``parsimony`` (``--scf``); ``quartets`` is the number of quartets sampled around each branch. ``partitions`` is an optional partition file for partitioned alignments. The annotated tree is ``outputs['cf.tree']``.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.compute_site_concordance_factors(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, quartets: 'int' = 100, method: 'str' = 'likelihood', partitions: 'str | Path | None' = None, model: 'str | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `tree`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``branches``, one row per branch with ``id``, ``scf``, ``scf_n``, ``sdf1``, ``sdf1_n``, ``sdf2``, ``sdf2_n``, ``sn``, ``label`` and ``length``.

## 5) `dnakit.phylogeny.iqtree_p11_cli_functions.date_phylogeny()`

- **Purpose:** Estimate a time-scaled tree and divergence times from an alignment and sampling dates (least-squares dating, LSD2).
- **Method:** ``dates`` is a file with one taxon name and date per line, or the word ``TAXNAME`` to read the date after the last ``|`` of each taxon name (``--date``). ``tip_date`` and ``root_date`` set the date of all tips or of the root (real number or ``YYYY-MM-DD``); ``confidence_replicates`` computes confidence intervals with that many resamples (``--date-ci``); ``clock_sd`` is the standard deviation of the log-normal relaxed clock (``--clock-sd``, default 0.2); ``exclude_outgroup`` is ``--date-no-outgroup`` and ``outlier_zscore`` removes tips whose z-score exceeds it (``--date-outlier``). With ``tree`` that topology is dated, otherwise IQ-TREE infers the ML tree first. The time tree is ``outputs['timetree.nwk']`` (also ``.nex`` and the LSD2 report ``.lsd``).
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.date_phylogeny(alignment: 'str | Path', dates: 'str | Path', output_directory: 'str | Path', *, tree: 'str | Path | None' = None, model: 'str | None' = None, tip_date: 'str | None' = None, root_date: 'str | None' = None, confidence_replicates: 'int | None' = None, clock_sd: 'float | None' = None, exclude_outgroup: 'bool' = False, outlier_zscore: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `dates`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 6) `dnakit.phylogeny.iqtree_p11_cli_functions.evaluate_tree_likelihood()`

- **Purpose:** Evaluate a fixed tree topology on an alignment: its log-likelihood, with or without branch-length optimisation.
- **Method:** With ``optimise_branch_lengths=True`` (IQ-TREE ``--tree-fix``) the topology is kept, and branch lengths and model parameters are optimised; the optimised tree is ``outputs['treefile']``. With ``False`` (``--show-lh``) the log-likelihood is computed on the branch lengths of ``tree`` without optimising them. ``model`` is an IQ-TREE model string such as ``GTR+F+G4``; without it IQ-TREE runs ModelFinder first.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.evaluate_tree_likelihood(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, optimise_branch_lengths: 'bool' = True, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `tree`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``log_likelihood``, ``model``, ``tree_length``.

## 7) `dnakit.phylogeny.iqtree_p11_cli_functions.infer_pathogen_tree()`

- **Purpose:** Infer a tree for low-divergence genomes (outbreak or pathogen data) with the CMAPLE algorithm (``--pathogen``).
- **Method:** IQ-TREE applies CMAPLE when the sequence divergence is low and its standard search otherwise; ``force_cmaple=True`` (``--pathogen-force``) applies CMAPLE regardless of divergence. ``alrt_replicates`` adds SH-aLRT branch supports (``--alrt``), ``sprta`` adds SPRTA supports (``--sprta``), with ``sprta_zero_length_branches`` (``--sprta-zero-branch``) and ``sprta_alternative_placements`` (``--sprta-other-places``). The tree is ``outputs['treefile']``.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.infer_pathogen_tree(alignment: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, force_cmaple: 'bool' = False, alrt_replicates: 'int | None' = None, sprta: 'bool' = False, sprta_zero_length_branches: 'bool' = False, sprta_alternative_placements: 'bool' = False, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 8) `dnakit.phylogeny.iqtree_p11_cli_functions.infer_polymorphism_aware_tree()`

- **Purpose:** Infer a population or species tree from allele counts with polymorphism-aware phylogenetic models (PoMo).
- **Method:** ``counts`` is an IQ-TREE counts file (``COUNTSFILE NPOP <n> NSITES <m>``, then ``CHROM POS`` and one column of four comma-separated nucleotide counts per population). The model is ``<substitution_model>+P`` with optionally ``+N<population_size>`` (virtual population size, IQ-TREE default 9), a sampling method (``sampling`` ``WB`` weighted binomial, ``WH`` weighted hypergeometric or ``S`` sampled) and ``+G<gamma_categories>``. ``tree`` fixes a starting tree (``-t``). The tree is ``outputs['treefile']``; ``model`` is the model IQ-TREE reports.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.infer_polymorphism_aware_tree(counts: 'str | Path', output_directory: 'str | Path', *, substitution_model: 'str' = 'HKY', population_size: 'int | None' = None, sampling: 'str | None' = None, gamma_categories: 'int | None' = None, tree: 'str | Path | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `counts`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 9) `dnakit.phylogeny.iqtree_p11_cli_functions.likelihood_mapping()`

- **Purpose:** Assess the phylogenetic signal of an alignment by likelihood mapping of random quartets (``--lmap``).
- **Method:** ``quartets`` is the number of quartets drawn (IQ-TREE uses all unique quartets when fewer exist); ``clusters`` is an optional NEXUS file defining sequence clusters (``--lmclust``). The maps are written as ``lmap.svg`` and ``lmap.eps`` in ``outputs``.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.likelihood_mapping(alignment: 'str | Path', output_directory: 'str | Path', *, quartets: 'int' = 10000, clusters: 'str | Path | None' = None, model: 'str | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``quartets`` and ``resolution`` (``fully_resolved``, ``partly_resolved`` and ``unresolved`` each with ``count`` and ``percent``).

## 10) `dnakit.phylogeny.iqtree_p11_cli_functions.reconstruct_ancestral_sequences()`

- **Purpose:** Reconstruct the ancestral sequence of every internal node by empirical Bayes (``--ancestral``).
- **Method:** The tree topology is fixed to ``tree``. ``min_probability`` is ``--asr-min`` (minimum posterior probability for a state to be called; below it IQ-TREE writes ``?``). The per-site posterior probabilities are in ``outputs['state']``.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.reconstruct_ancestral_sequences(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, min_probability: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `tree`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``nodes`` (internal node names), ``sites`` and ``sequences`` (node name to its most likely sequence; None when the table is larger than 2 million states).

## 11) `dnakit.phylogeny.iqtree_p11_cli_functions.run_symmetry_tests()`

- **Purpose:** Test the stationarity, reversibility and homogeneity assumptions with the matched-pair tests of symmetry.
- **Method:** IQ-TREE runs three tests on every sequence pair (symmetry, marginal symmetry, internal symmetry; Naser-Khdour et al. 2019). ``only_test=True`` runs them and stops (``--symtest-only``); with ``False`` the tests are followed by the tree analysis (``--symtest``). ``remove_partitions`` ``bad`` or ``good`` removes the partitions that fail or pass the tests (``--symtest-remove-bad`` / ``--symtest-remove-good``; needs ``partitions``), judged by ``test_type`` ``MAR`` or ``INT`` and the cut-off ``p_value`` (default 0.05); ``keep_zero`` keeps NA results (``--symtest-keep-zero``).
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.run_symmetry_tests(alignment: 'str | Path', output_directory: 'str | Path', *, only_test: 'bool' = True, remove_partitions: 'str | None' = None, test_type: 'str | None' = None, p_value: 'float | None' = None, keep_zero: 'bool' = False, partitions: 'str | Path | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``tests`` (one row per partition from ``outputs['symtest.csv']``: ``name``, ``sym_sig``, ``sym_non``, ``sym_pval``, ``mar_sig``, ``mar_non``, ``mar_pval``, ``int_sig``, ``int_non``, ``int_pval``).

## 12) `dnakit.phylogeny.iqtree_p11_cli_functions.select_substitution_model()`

- **Purpose:** Select the best-fit substitution model (ModelFinder) and optionally merge partitions (PartitionFinder-style).
- **Method:** ``search`` is IQ-TREE's ``-m`` keyword: ``TESTONLY`` (standard model selection only), ``TEST`` (standard selection, then tree inference), ``MF`` (extended selection with FreeRate models, no tree search) or ``MFP`` (extended selection, then tree inference; IQ-TREE's default). With ``partitions`` (a NEXUS or RAxML partition file; ``partition_mode`` ``proportional`` ``-p``, ``equal`` ``-q`` or ``unlinked`` ``-Q``) a model is chosen per partition, and with ``merge=True`` the partitions are also merged to increase the model fit (``--merge``; ``merge_algorithm`` ``greedy``, ``rcluster`` or ``rclusterf``, IQ-TREE's default). ``criterion`` is ``AIC``, ``AICc``, ``BIC`` (default) or ``mAIC``; ``model_set`` (``--mset``), ``model_subset`` (``--msub``), ``rate_models`` (``--mrate``, for example ``E,I,G,I+G,R``) and ``max_rate_categories`` (``--cmax``) restrict the search.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.select_substitution_model(alignment: 'str | Path', output_directory: 'str | Path', *, search: 'str' = 'MFP', partitions: 'str | Path | None' = None, partition_mode: 'str' = 'proportional', merge: 'bool' = False, merge_algorithm: 'str | None' = None, criterion: 'str | None' = None, model_set: 'str | None' = None, model_subset: 'str | None' = None, rate_models: 'str | None' = None, max_rate_categories: 'int | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``best_model``, ``selection_criterion``, ``log_likelihood`` and ``model_table`` (the first 20 rows of the ranked model list: ``model``, ``log_likelihood``, ``aic``, ``aicc``, ``bic``); with merging ``best_scheme`` holds the path of the merged partition scheme.

## 13) `dnakit.phylogeny.iqtree_p11_cli_functions.simulate_alignment()`

- **Purpose:** Simulate sequence alignments along a tree with AliSim (IQ-TREE ``--alisim``).
- **Method:** Give a ``tree`` file, or ``random_tree_taxa`` to generate a random tree (``-t RANDOM{<random_tree_model>/<taxa>}``; ``yh`` Yule-Harding, ``u`` uniform, ``bd`` birth-death, ...). ``model`` is an IQ-TREE model string (``HKY+F{0.3/0.2/0.2/0.3}+G4``); ``length`` is the length of the root sequence and ``alignments`` the number of datasets. ``output_format`` is ``fasta`` or ``phylip``. ``template_alignment`` (``-s``) makes AliSim copy the gap pattern and mimic the site frequencies and rates of a real alignment; ``partitions`` (``-p``) simulates a partitioned alignment; ``branch_scale`` multiplies all branch lengths; ``insertion_deletion_rates`` is the pair ``(insertion, deletion)`` of the indel model (``--indel``). Set ``seed`` for reproducible output.
- **API:** `dnakit.phylogeny.iqtree_p11_cli_functions.simulate_alignment(output_directory: 'str | Path', *, tree: 'str | Path | None' = None, random_tree_taxa: 'int | None' = None, random_tree_model: 'str' = 'yh', model: 'str' = 'JC', length: 'int' = 1000, alignments: 'int' = 1, output_format: 'str' = 'fasta', template_alignment: 'str | Path | None' = None, branch_scale: 'float | None' = None, partitions: 'str | Path | None' = None, insertion_deletion_rates: 'Sequence[float] | None' = None, sequence_type: 'str | None' = None, seed: 'int | None' = None, prefix: 'str' = 'simulated', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``alignment_files`` (paths) and ``sequences`` and ``sites`` of the first alignment.

