# Bioinformatics Python APIs

`dnakit.bioinformatics` adds traditional bioinformatics, DNA representations, generative evaluation, statistics and molecular design. These development APIs are **Python only**. Coverage below describes actual behavior, not complete equivalence to Biopython, SeqKit, scikit-bio, UGENE, iLearnPlus, iFeatureOmega, seqme, pydna or DNA Chisel.

## Implementation and installation {#dependencies}

DNAKit implements domain logic itself. NumPy/SciPy supply arrays, linear algebra, distances and numerical solvers; scikit-learn supplies general machine-learning algorithms. There are no iFeatureOmega, iLearnPlus, seqme or DNA Chisel function adapters.

| Extra | Purpose | Packages |
| --- | --- | --- |
| `bio-analysis` | Native encodings, evaluation, statistics; general feature analysis | NumPy, SciPy, scikit-learn |
| `bio-formats` | Alignment, tree, motif, BLAST/HMMER, ABI, SnapGene/GenBank codecs | Biopython |
| `bio-reads` | FASTA indexing and random BAM access | pysam/htslib |
| `bio-plots` | Trees, dotplots, chromatograms, gels | Matplotlib |
| `bio-assembly` | Complex overlap assembly and permissive construct parsing | pydna |
| `bio` | All extras above | External executables are installed separately |
| `bio-validation` | Independent numerical references in tests, not runtime dependencies | scikit-bio, statsmodels |

Install the working tree with the complete Python dependencies:

```bash
python -m pip install -e ".[bio]"
```

Install only numerical and feature-analysis dependencies:

```bash
python -m pip install -e ".[bio-analysis]"
```

Import submodules explicitly. Dependencies load when needed; APIs do not install packages, download datasets or submit remote jobs automatically. Some codec APIs explicitly return Biopython objects. Native results generally use dataclasses, dictionaries and NumPy arrays.

## Alignments, BLAST, profile HMMs and trees

| Module | API | Behavior |
| --- | --- | --- |
| `alignment` | `MultipleAlignment.from_sequences`, `.distances` | Native equal-width, unique-ID alignment; p-distance with pairwise deletion or gaps as a state |
| `alignment` | `read_alignment`, `write_alignment` | Biopython single-alignment FASTA/Clustal/Stockholm/PHYLIP codecs |
| `alignment` | `align_sequences` | MAFFT `--auto`, MUSCLE **5**, Clustal Omega; restore IDs and verify sequence content |
| `alignment` | `build_blast_database`, `local_blast` | Local nucleotide database and BLAST XML output |
| `alignment` | `parse_search_results` | Streaming Bio.SearchIO QueryResult objects for BLAST and protein HMMER2/3 formats |
| `alignment` | `build_profile`, `search_profile`, `parse_nhmmer_table` | hmmbuild, nucleotide nhmmer / protein hmmsearch; native nhmmer table parser |
| `phylogeny` | `TreeNode`, `build_tree` | Native neighbor joining and size-weighted UPGMA |
| `phylogeny` | `read_tree`, `write_tree`, `plot_tree` | Newick/Nexus/phyloXML/NeXML codecs and tree plots |
| `phylogeny` | `robinson_foulds`, `consensus_tree`, `prune_tree` | Native rooted/unrooted RF, strict-majority/strict consensus, pruning |

Remote BLAST already has explicit submission/polling APIs in `dnakit.search.blast`; this layer does not duplicate them.

```python
from dnakit.bioinformatics.alignment import MultipleAlignment
from dnakit.bioinformatics.phylogeny import build_tree

msa = MultipleAlignment.from_sequences({"a": "ACGT", "b": "ACTT", "c": "TCGT"})
tree = build_tree(msa.distances(), msa.ids, method="nj")
print(tree.tips)
```

NJ retains negative branch estimates and has an arbitrary root; UPGMA assumes a molecular clock. RF requires identical taxa and optionally divides by the combined split count. Consensus assumes the same root interpretation: threshold `0.5` keeps strictly majority clades, `1` keeps shared clades. GME/BME, rearrangement search and maximum-likelihood phylogenetics are not provided.

## FASTX, intervals and NGS

| Module | API | Contract |
| --- | --- | --- |
| `reads` | `Read`, `read_fastx`, `write_fastx` | FASTA, multiline FASTQ, gzip, integer Phred scores; writer returns record count |
| `reads` | `pair_reads` | Two read iterators paired by IDs after removing `/1` or `/2`; ReadPair uses None for an orphan mate |
| `reads` | `repair_fastq` | Recover valid four-line records, report output/recovered count/discarded lines |
| `reads` | `arrange_reads`, `common_reads`, `rename_reads`, `sequence_digest` | Sorting/seeded shuffle, intersection, regex renaming, framed SHA256 |
| `reads` | `Region`, `read_regions` | BED or selected GTF feature rows; convert GTF coordinates exactly once |
| `reads` | `index_reference`, `IndexedReference.fetch` | Explicit `.fai` creation, context-managed indexed access, reverse complements |
| `reads` | `alignment_summary`, `coverage_consensus` | pysam decoding, native flag counts and A/C/G/T depth/majority consensus |
| `ngs` | `map_reads` | Private BWA MEM / Bowtie2 reference index and SAM mapping |
| `ngs` | `prepare_bam` | samtools name-sort/fixmate/coordinate-sort/markdup/index/stats |
| `ngs` | `call_variants` | bcftools mpileup/call/norm and indexed VCF; explicit ploidy |
| `ngs` | `assemble_reads` | SPAdes single/paired short-read or Nanopore hybrid assembly |
| `ngs` | `quality_control`, `dna_seq_workflow` | FastQC; QC→BWA→BAM→variant calling with all intermediates retained |
| `ngs` | `rna_seq_workflow` | HISAT2→BAM→StringTie→featureCounts, strandedness 0/1/2 |
| `ngs` | `chip_seq_workflow` | Prepared treatment/control BAMs → MACS3 peaks and coverage |

```python
from dnakit.bioinformatics.reads import IndexedReference, Region, read_fastx, write_fastx

# Explicit Phred+64 to Phred+33, not a Solexa log-odds conversion.
n = write_fastx(read_fastx("input.fq.gz", phred_offset=64), "phred33.fq", phred_offset=33)
with IndexedReference("reference.fa") as reference:  # Existing .fai required.
    sequence = reference.fetch(Region("chr1", 100, 180, "-"))
```

FASTA→FASTQ requires explicit `default_quality`. Repair never invents bases or qualities, handles four-line FASTQ only and cannot guarantee unambiguous resynchronization. Pairing bounds pending read payload with `max_pending`; seen-ID sets still grow with total records.

Python coordinates are **0-based half-open**. BAM consensus is reference-forward regardless of `Region.strand`. It excludes secondary, supplementary, duplicate, QC-failed and low-quality reads; ties and insufficient depth produce N. Indels are retained in alignment/VCF data, not inserted into consensus.

Install specialist programs separately on PATH. `ExecutionConfig(threads=4, timeout_seconds=3600)` controls execution without a shell. Output directories must not exist; intermediates and diagnostics are retained. `max_output_bytes` bounds process logs; `max_artifact_bytes` bounds stdout or explicitly monitored output files, not the entire third-party working directory. Inputs are not overwritten; the explicit indexing API creates sidecar indexes next to its input.

```python
from dnakit.bioinformatics import ExecutionConfig
from dnakit.bioinformatics.ngs import dna_seq_workflow

stages = dna_seq_workflow(
    "reference.fa", "R1.fq.gz", "new_analysis", second="R2.fq.gz",
    ploidy=2, config=ExecutionConfig(threads=4),
)
vcf = stages["variants"].artifacts["variants.vcf.gz"]
```

Variant calling requires an indexed reference and sorted indexed BAM. The RNA workflow processes one sample; inference across biological replicates is separate. ChIP begins with prepared BAMs and composes with mapping/BAM APIs. There is no assembly-browser GUI.

## Native DNA encodings and general feature analysis

`encoding.encode_dna({id: sequence}, descriptor, parameters=..., property_tables=..., labels=..., training_ids=...)` returns `EncodingResult(ids, names, values, descriptor, parameters, schema_version)`. The read-only matrix is samples×features. Schema `dnakit.dna-encoding.v1` accepts canonical A/C/G/T and rejects unknown parameters.

| Family | Parameters and convention |
| --- | --- |
| Kmer / RCKmer | `k=3, normalized=True`; overlapping frequencies / lexicographically canonical reverse-complement aggregation |
| Mismatch / Subsequence | `k=3, normalized=True`; `mismatches=1` Hamming-neighborhood sums / `delta=0.5` span-weighted subsequences |
| NAC / ANF / ENAC | Single-base frequency / current-base prefix frequency / `window=5` local composition |
| binary / PS2 / PS3 / PS4 | Positional k-mer one-hot in A,C,G,T order |
| CKSNAP / ASDC / LPDF / DBE | `gap=3, normalized=True` spaced pairs / all ordered position pairs / prefix pair frequency / four-bit dinucleotides |
| NCP / EIIP / PseEIIP | Three chemical classifications / positional EIIP / triplet frequency × EIIP sum |
| Z_curve_9bit/12bit/36bit/48bit/144bit | Phase/prefix-conditioned x/y/z contrasts divided by total overlapping windows |
| DPCP / TPCP | Positional standardized physicochemical profiles from explicit tables |
| DAC / DCC / DACC, TAC / TCC / TACC | `lag=3`, centered auto/cross covariance; pure cross requires at least two tables |
| NMBroto / Moran / Geary | `lag=3`, property products / normalized autocorrelation / squared-difference statistic |
| PseDNC / PseKNC, PCPseDNC/TNC, SCPseDNC/TNC | `lag=3, weight=0.05`; PseKNC also `k=3`; composition and lag terms share a denominator |
| MMI | Positional information contributions aggregated by nucleotide multiset: 10 pair-MI and 20 triplet co-information cells |
| PSTNPss / PSTNPds | Explicit binary-training-set positional triplet frequency differences; ds folds T→A and G→C |
| KNN | `neighbors=3`; positive fraction among Hamming-nearest training samples, ties ordered by ID |

All names are in `encoding.DNA_ENCODINGS`. Positional encodings require equal lengths. `property_tables` maps property names to all 16 dinucleotide or 64 trinucleotide values, standardized using the population SD over the complete alphabet. No unattributed property tables are bundled. NMBroto/Moran/Geary use dinucleotide tables. Short inputs, constant profiles and undefined denominators raise errors.

Pseudo type-I uses lagged mean squared differences across properties; SC uses per-property lagged product means. Subsequence weights are `delta ** span_length`, divided by total weight when normalized. MMI uses consistent marginals from same-length windows; pair/triplet cells sum to total MI/co-information. These explicit conventions **do not promise numerical identity with iFeatureOmega's similarly named vectors**. PSTNP/KNN read training labels only and exclude each training sample from its own reference set.

```python
from dnakit.bioinformatics.encoding import encode_dna
from dnakit.bioinformatics.features import analyze_features

encoded = encode_dna({"a": "AACGTG", "b": "CCGTAT", "c": "ATGCAA"}, "RCKmer", parameters={"k": 2})
fitted = analyze_features(encoded.values, feature_names=encoded.names, method="pca", components=2)
held_out = encode_dna({"test": "GCTAAC"}, "RCKmer", parameters={"k": 2})
values = fitted.transform(held_out.values, feature_names=held_out.names)
```

`features.analyze_features` supports standard/minmax/L1/L2/none scaling; PCA, t-SNE, supervised LDA, Latent Dirichlet Allocation, chi2, F-classification/regression, MI, Information Gain (training quantile bins, `parameters={"bins": 10}`), Pearson, variance and none. Supervised methods require labels; chi2/topic models require nonnegative data. t-SNE has no held-out transform. Reuse the fitted training pipeline without refitting on held-out samples.

`cluster_features` supports kmeans, minibatch_kmeans, gmm, spectral, affinity_propagation, mean_shift, dbscan, agglomerative and birch, returning fitted sklearn estimators. Native `mcl` applies stochastic-matrix expansion/inflation to symmetric nonnegative affinities and returns index groups, raising on nonconvergence. Clustering does not automatically scale inputs.

## Generative and multi-objective evaluation

| `metrics` API | Definition |
| --- | --- |
| `kernel_distance` | Gaussian MMD² or cubic polynomial KID on embeddings; default unbiased estimates may be negative |
| `manifold_metrics` | Conventional PRDC precision/recall/density/coverage, open kNN balls, self excluded |
| `clipped_manifold_metrics` | Salvy et al. median-clipped density with leave-one-out normalization and beta-binomial calibrated coverage; closed balls and equal sample counts |
| `authenticity` | Fraction where generated-to-nearest-training distance exceeds that training point's nearest-other-training distance, in [0,1] |
| `pareto_ranks` | Explicit per-column `maximize`; zero-based nondomination ranks, duplicates share ranks |
| `hypervolume` | Explicit directions and worse reference point; native exact recursion, at most five dimensions/1000 points with work limit |
| `property_hit_rate` | Inclusive `lower/upper` bounds, joint/per-property rates and Boolean mask |
| `property_kl` | KL(generated∥reference), shared histogram edges, explicit pseudocount, natural logarithms |
| `conformity_score` | Cross-fitted isotropic Gaussian KDE and validation-density empirical scores |
| `sequence_properties` | Named scalar Python predictors → stable IDs/property names/matrix for the metrics above |

Authenticity requires the **generator's training set**, not a held-out reference. Fix embedding model, scaling and k before comparing methods. Conventional PRDC and clipped metrics differ; KID is not FID. Histogram property KL differs from seqme's KDE estimate. Conformity returns the value, fold SE and fold scores; it is not a posterior class probability. Predictor callbacks may wrap existing DNAKit models; evaluation never fits them on generated samples.

## Statistics, compositions and count inference

| Module | API | Scope |
| --- | --- | --- |
| `statistics` | `distance_matrix`, `distance_test` | Distance validation; PERMANOVA, ANOSIM, centroid PERMDISP with seeded permutations |
| `statistics` | `mantel`, `bioenv` | Pearson/Spearman Mantel; bounded environmental-subset BIOENV |
| `statistics` | `pcoa`, `ordination` | PCoA/CA/CCA/RDA; report negative PCoA eigenvalues, constrained axes only for CCA/RDA |
| `statistics` | `composition_transform` | CLR/ILR/ALR and inverses, orthonormal Helmert ILR basis |
| `statistics` | `ancom` | Pairwise log ratios, per-taxon Holm correction and W, without an invented W rejection cutoff |
| `statistics` | `alpha_diversity`, `beta_diversity` | observed/Shannon/Simpson/inverse Simpson; Bray-Curtis/Jaccard/Euclidean/cityblock/cosine |
| `phylogeny` | `unifrac` | Native rooted weighted/unweighted UniFrac, exact tree-tip/taxa-column agreement |
| `counts` | `ancom_bc` | Original ANCOM-BC log-linear OLS, sandwich variance, three-component sampling-bias EM |
| `counts` | `differential_expression` | Unshrunk native NB2 GLMs, size factors, explicit contrasts, Wald/BH |
| `counts` | `adjust_pvalues` | One-dimensional Benjamini–Hochberg correction |

Permutation p-values use `(extreme+1)/(permutations+1)`; zero permutations return None. No restricted/repeated-measures permutations. PERMDISP currently requires Euclidean distances. Composition transforms and ANCOM-BC require explicit zero handling; ANCOM accepts an explicit pseudocount. UniFrac requires a biologically rooted tree and nonnegative branches.

ANCOM-BC takes samples×taxa and a full-rank numeric design including an explicit intercept; at least ten taxa are required. It returns feature×covariate **natural-log** effects, SE, Wald, p/q, sampling bias and EM iteration counts. BH correction is per covariate. This is not ANCOM-BC2 and excludes structural-zero detection, global tests and bias-variance adjustment. Degenerate variance/nonconvergence raise errors.

NB2 takes integer samples×genes counts, a numeric design and explicit contrast, with at least two residual degrees of freedom and actual biological replicates. Default size factors are median ratios over genes positive in every sample; dispersion defaults to within-design residual moments. Explicit size factors/dispersions are accepted. Results include log2 fold changes, log2 SE, Wald p/q, dispersions and convergence flags; all-zero genes and failed fits have NaN inference. **DESeq2/edgeR empirical-Bayes trends, effect shrinkage and outlier strategies are not implemented.**

## Specifications, motifs and cloning

`optimization.specification(kind, region=None, weight=1, **parameters)` creates a native rule. Custom `Specification(name, evaluator)` accepts `(sequence, circular) -> SpecEvaluation`, separating nonnegative violation, utility to maximize, failure locations and messages.

| kind | Parameters |
| --- | --- |
| `gc` | `minimum=0.4, maximum=0.6, window=whole-region length` |
| `avoid_pattern`, `pattern_count` | IUPAC `pattern`, `both_strands=False`; latter also `count=1` |
| `sequence`, `translation`, `length` | `pattern` / `protein` / `minimum=1, maximum=∞` |
| `avoid_changes` | Full original `sequence`, optional absolute `positions`, otherwise lock region/all positions |
| `unique_kmers`, `avoid_stop`, `avoid_codons` | `k=12, both_strands=True` / `allow_terminal=True` / `codons` |
| `avoid_rare_codons` | All-64 `usage`, relative-to-best-synonym `threshold=0.1` |
| `hairpin` | Exact complementary stem: `match_length=12, window=200, minimum_loop=3` |
| `heterodimer`, `avoid_matches` | `references, match_length=12`; latter also `both_strands=True` |
| `regions_compatibility` | At least two linear absolute `regions`, `match_length=12`, avoid inter-region complementarity |
| `tm` | `minimum, maximum`, optional existing thermodynamic `conditions` |
| `cai`, `codon_usage`, `harmonize` | All-64 `usage`; harmonize also per-codon `relative_adaptiveness` |

Structural rules are explicit exact-complement models, not full thermodynamic energy calculations. `blast_match_specification(database, minimum_identity=0.8, minimum_length=20)` uses specialized BLAST HSP detection with sequence caching. No HSPs does not prove no homology. `avoid_matches` is native exact k-mer matching; no Bowtie genomic-constraint wrapper is provided.

```python
from dnakit.bioinformatics.optimization import optimize_sequence, specification

result = optimize_sequence(
    "GCTGCT", preserve_translation=True,
    constraints=[specification("avoid_pattern", pattern="GCT")], seed=0,
)
print(result.sequence, result.constraints_satisfied, result.edits)
# result.write_report("new_report.json")  # Fresh file, failure spans and search counts.
```

Hard constraints precede weighted objectives. Search uses greedy single-base substitutions or whole synonymous-codon replacements with `preserve_translation=True`. Length is fixed; length constraints do not insert/delete bases. Circular windows/patterns cross the origin with end coordinates potentially beyond length. Search exhaustion is not proof of infeasibility; satisfied constraints do not establish a global optimum.

GenBank qualifiers `dnakit_constraint` and `dnakit_objective` accept JSON such as `{"kind":"gc","minimum":0.4,"maximum":0.6}` through `specifications_from_record`. Only exact forward/simple features convert automatically; reverse, fuzzy and compound features require explicit Python declarations.

| Module | API | Scope |
| --- | --- | --- |
| `motifs` | `MotifModel`, `build_motif`, `fpr_threshold`, `scan_motif_model` | Native counts/PWM/PSSM, background, degenerate consensus, conservative quantized FPR, both strands/circular scanning |
| `motifs` | `read_motifs`, `dotplot` | JASPAR/MEME/TRANSFAC codecs; native forward/reverse exact k-mer dotplot |
| `cloning` | `Duplex.from_sequence`, `.cut`, `.ligate`, `.reverse_complement`, `.ends` | Native explicit 5′→3′ strands, offsets, sticky/blunt ends and topology |
| `cloning` | `CloningPlan.add`, `.apply`, `.steps` | Executable history retaining every intermediate |
| `cloning` | `assemble_duplexes` | Specialist pydna Gibson/homologous overlap assembly |
| `cloning` | `read_sanger`, `plot_chromatogram`, `sanger_consensus` | ABI codecs, native Mott trim, signal plots, native bidirectional local mapping and reference-coordinate consensus |
| `cloning` | `read_constructs` | Strict SnapGene/GenBank codecs or explicit pydna permissive parsing |
| `cloning` | `simulate_gel` | Qualitative log-length migration image, returned as a Matplotlib Axes |

Duplex offset places the Crick 3′ end relative to the Watson 5′ end. Paired regions must be complementary; nicks, internal gaps and mismatches are unsupported. Gel rendering does not model dye/buffer dynamics or calibrated migration. Sanger consensus filters low-quality/weak mappings and reports N for ties/uncovered sites, without inserting indels. Trace-sample and base coordinates remain distinct.

## Python workflows and validation

`workflows.Step(name, function, dependencies={parameter: step_name}, parameters=...)` and `run_workflow` implement a native DAG with unique-name/dependency/cycle checks. Execution can fail fast or record failures and block dependents. `query_regions(left, right, operation='intersection')` supports union/intersection/difference per chromosome/strand. There is no visual designer.

Run the offline API example (requires only `bio-analysis`):

```bash
python others/examples/bioinformatics_api.py
```

Run the new tests; real-tool integration skips when programs are not on PATH:

```bash
python -m pytest others/tests/unit/test_bioinformatics_native.py others/tests/unit/test_bioinformatics_io_workflows.py others/tests/unit/test_bioinformatics_external.py -q
```

References: [Biopython motifs](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html), [PRDC](https://github.com/clovaai/generative-evaluation-prdc), [clipped metrics](https://arxiv.org/abs/2507.01761), [iLearnPlus encoding families](https://academic.oup.com/nar/article/49/10/e60/6154472), [ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7), [DNA Chisel specifications](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html). Formula conventions are defined here and in API docstrings.
