# Bioinformatics FAQ

Function usage and examples are included in their task tutorials. This page covers dependencies, parameter catalogs, implementation conventions and validation.

`dnakit.bioinformatics` adds traditional bioinformatics, DNA representations, generative evaluation, statistics and molecular design. These development APIs are **Python only**. Coverage below describes actual behavior, not complete equivalence to Biopython, SeqKit, scikit-bio, UGENE, iLearnPlus, iFeatureOmega, seqme, pydna, DNA Chisel, GeneBlocks or sourmash.

## Which dependencies are required? {#dependencies}

DNAKit implements domain logic itself. NumPy/SciPy supply arrays, linear algebra, distances and numerical solvers; scikit-learn supplies general machine-learning algorithms. scikit-bio, seqme, DNA Chisel, Biotite, Edlib, GeneBlocks, sourmash, Cutadapt, PyRanges, pyfastx, dnaio, mappy, pyBigWig, pybedtools, gffutils, pySBOL3, scadnano, polars-bio, DnaCauldron, DNA Features Viewer and gtars have explicit adapters; BWA, BWA-MEM2, Bowtie2, SAMtools, BCFtools, SPAdes, IQ-TREE2, FastTree, FreeBayes, VCFtools, vcfanno, SURVIVOR, GEMINI, Ensembl VEP, Picard, MACS3, NCBI BLAST+, minimap2, BEDTools, fastp, FastQC, MultiQC, NanoFilt, Filtlong, NanoPlot, Porechop, Trimmomatic, AfterQC, SeqMagick, Shasta, Merqury, BUSCO, DeepVariant, cuteSV, Sniffles2, DELLY, HTSlib, seqtk, WhatsHap, BamTools, VSEARCH, BEDOPS, Flye, hifiasm, MEGAHIT and GATK4 use caller-provided external executables, while iFeatureOmega and iLearnPlus remain external-backend boundaries.

| Extra | Purpose | Packages |
| --- | --- | --- |
| `bio-analysis` | Native encodings, evaluation, statistics; general feature selection, dimensionality reduction, and clustering | NumPy, SciPy, scikit-learn |
| `bio-formats` | Alignment, tree, motif, BLAST/HMMER, ABI, SnapGene/GenBank codecs | Biopython |
| `bio-reads` | FASTA indexing and random BAM access | pysam/htslib |
| `bio-plots` | Trees, dotplots, chromatograms, gels | Matplotlib |
| `bio-assembly` | Complex overlap assembly and permissive construct parsing | pydna |
| `bio-seqme` | seqme sequence-design metrics; model extras are separate | seqme (BSD-3-Clause) |
| `bio-dnachisel` | DNA Chisel constraint and objective optimization adapter | DNA Chisel (MIT) |
| `bio-skibio` | beta-diversity and PERMANOVA statistics adapter | scikit-bio (BSD-3-Clause) |
| `bio-biotite` | bounded optimal nucleotide pairwise alignment | Biotite (BSD-3-Clause) |
| `bio-edlib` | bounded edit-distance alignment | Edlib (MIT) |
| `bio-geneblocks` | bounded DNA sequence diff blocks | GeneBlocks (MIT) |
| `bio-sourmash` | bounded local DNA MinHash comparison | sourmash (BSD-3-Clause) |
| `bio-cutadapt` | Cutadapt FASTQ adapter trimming | Cutadapt (MIT) |
| `bio-pyranges` | BED interval summary and merge | PyRanges (MIT) |
| `bio-pyfastx` | bounded FASTA/FASTQ summary | pyfastx (MIT) |
| `bio-dnaio` | bounded read-only FASTA/FASTQ stream summary | dnaio (MIT) |
| `bio-mappy` | bounded local sequence mapping | mappy/minimap2 (MIT) |
| `bio-pybigwig` | bounded bigWig region statistics | pyBigWig (MIT) |
| `bio-pybedtools` | bounded BED intersection; bedtools is still separate | pybedtools (MIT) + BEDTools (MIT) |
| `bio-gffutils` | read-only queries against an existing GFF/GTF FeatureDB | gffutils (MIT) |
| `bio-pysbol3` | read-only inventory of a local SBOL3 RDF document | pySBOL3 (MIT) |
| `bio-scadnano` | read-only inventory of a local scadnano `.sc` design | scadnano (MIT) |
| `bio-polarsbio` | bounded local BED/CSV interval overlap | polars-bio (Apache-2.0) |
| `bio-dnacauldron` | Gibson/Type IIS assembly simulation | DnaCauldron (MIT) |
| `bio-dna-features-viewer` | read-only inventory of local GenBank/GFF features | DNA Features Viewer (MIT) |
| `bio-gtars` | read-only inventory of a local BED RegionSet | gtars (BSD-2-Clause) |
| `bio-pyfaidx` | read-only FASTA interval adapter | pyfaidx (BSD-3-Clause) |
| `bio-vcfpy` | bounded VCF summary parsing | vcfpy (MIT) |
| `bio-cyvcf2` | bounded VCF summary parsing | cyvcf2 (MIT) |
| `bio` | All extras above | External executables are installed separately |
| `bio-validation` | scikit-bio statistics adapter plus independent numerical references | scikit-bio, statsmodels |

BWA-MEM2, BLAST+, minimap2, BEDTools, MACS3, FastTree, NanoFilt, Filtlong, NanoPlot, Porechop, Trimmomatic, AfterQC, SeqMagick, Shasta, Merqury, BUSCO, DeepVariant, cuteSV, Sniffles2, DELLY, HTSlib, seqtk, WhatsHap, BamTools, VSEARCH, BEDOPS, Flye, hifiasm, MEGAHIT, GEMINI, Ensembl VEP and GATK4, DIAMOND, MMseqs2, Kraken2, Prokka, Bakta, CheckM2, Bracken, Centrifuge, Kaiju, Unicycler, GTDB-Tk, HUMAnN and MethylDackel also have formal Agent/MCP facades: `dnakit.ngs.map_reads(method="bwa-mem2")`, `dnakit.alignment.build_blast_database()` / `dnakit.alignment.local_blast()`, `dnakit.alignment.map_long_reads()`, `dnakit.alignment.intersect_bed()`, `dnakit.ngs.call_chip_peaks()`, `dnakit.phylogeny.infer_fast_tree()`, `dnakit.sequencing_qc.filter_nanopore_reads()`, `filter_long_reads()`, `summarize_long_reads()`, `trim_porechop()`, `trim_trimmomatic()`, `run_afterqc()`, `dnakit.fasta_fastq.convert_seqmagick()`, `dnakit.assembly.assemble_shasta()`, `dnakit.assembly.evaluate_merqury()`, `dnakit.assembly.run_busco()`, `dnakit.variants.call_deepvariant()`, `dnakit.variants.call_cutesv()`, `dnakit.variants.call_sniffles2()`, `dnakit.variants.call_delly()`, `dnakit.ngs.compress_index_vcf_htslib()`, `dnakit.ngs.sample_sequences_seqtk()`, `dnakit.ngs.phase_variants_whatshap()`, `dnakit.ngs.sort_bam_bamtools()`, `dnakit.fasta_fastq.cluster_vsearch()`, `dnakit.alignment.sort_bed_bedops()`, `dnakit.assembly.assemble_flye()`, `dnakit.assembly.assemble_hifiasm()`, `dnakit.assembly.assemble_megahit()` and `dnakit.variants.call_gatk_haplotypecaller()`, `dnakit.variants.load_gemini_database()`, `dnakit.variants.query_gemini_database()`, `dnakit.variants.annotate_vcf_ensembl_vep()`, `dnakit.alignment.search_diamond_blastp()`, `search_diamond_blastx()`, `search_mmseqs_protein()`, `cluster_mmseqs_protein()`, `dnakit.microbiome.classify_kraken2()`, `dnakit.annotation.annotate_prokka()`, `annotate_bakta()`, `predict_checkm2()`, `dnakit.microbiome.estimate_bracken_abundance()`, `classify_centrifuge()`, `classify_kaiju()`, `dnakit.assembly.assemble_unicycler()`, `dnakit.microbiome.classify_gtdbtk()`, `profile_humann()` and `dnakit.epigenome.extract_methylation_methyldackel()`. They run only caller-installed programs and expose bounded fixed subsets of each CLI. The cited records are Vasimuddin et al. (2019), DOI `10.1109/IPDPS.2019.00041`; Camacho et al. (2009), DOI `10.1186/1471-2105-10-421`; Li (2018), DOI `10.1093/bioinformatics/bty191`; Quinlan & Hall (2010), DOI `10.1093/bioinformatics/btq033`; Price, Dehal & Arkin (2010), DOI `10.1371/journal.pone.0009490`; De Coster et al. (2018), DOI `10.1093/bioinformatics/bty149`; Wick et al. (2017), DOI `10.1099/mgen.0.000132`; Bolger, Lohse & Usadel (2014), DOI `10.1093/bioinformatics/btu170`; and De Coster & Rademakers (2023), DOI `10.1093/bioinformatics/btad311`; Chen et al. (2017), DOI `10.1186/s12859-017-1469-3`; Shafin et al. (2020), DOI `10.1038/s41587-020-0503-6`. The survey identifies MACS3 and AfterQC as permissive upstream records and Filtlong/Porechop/Trimmomatic as GPL-3.0; the official Filtlong project supplies no DOI and Trimmomatic adapter sequences have separate Illumina terms. AfterQC's README recommends fastp for maintained C++ workflows, so this route preserves historical AfterQC reproducibility rather than claiming current upstream maintenance. SeqMagick has no DOI in its official project; GATK4 should be cited as McKenna et al. (2010), DOI `10.1101/gr.107524.110`, and the current GATK source carries Apache-2.0.

The VCF, duplicate-marking and consensus facades are `dnakit.variants.annotate_vcf_vcfanno()`, `dnakit.variants.merge_structural_variants_survivor()`, `dnakit.ngs.mark_duplicates_picard()` and `dnakit.alignment.consensus_spoa()`. They run caller-installed vcfanno, SURVIVOR, Picard and SPOA with fixed parameter subsets. Cite Pedersen et al. (2016), DOI `10.1186/s13059-016-0973-5`; Jeffares et al. (2017), DOI `10.1038/NCOMMS14061`; Broad Institute *Picard Toolkit* (2019); and SPOA, DOI `10.1093/bioinformatics/btg109`. The four official repositories record MIT licenses.

The GEMINI facades are `dnakit.variants.load_gemini_database()` and `dnakit.variants.query_gemini_database()` for local database construction and bounded SQL queries. Cite Paila et al. (2013), DOI `10.1371/journal.pcbi.1003153`; the official repository records MIT. Annotation databases, resource downloads and the complete GEMINI CLI remain caller-prepared.

The Ensembl VEP facade is `dnakit.variants.annotate_vcf_ensembl_vep()`. Cite McLaren et al. (2016), DOI `10.1186/s13059-016-0974-4`; the official repository records Apache-2.0. VEP caches, reference sequences, plugins and the complete CLI remain caller-prepared.

Install the working tree with the complete Python dependencies:

```bash
python -m pip install -e ".[bio]"
```

Install only numerical and feature-analysis dependencies:

```bash
python -m pip install -e ".[bio-analysis]"
```

Import submodules explicitly. Dependencies load when needed; APIs do not install packages, download datasets or submit remote jobs automatically. Some codec APIs explicitly return Biopython objects. Native results generally use dataclasses, dictionaries and NumPy arrays.

## What conventions apply to alignments, searches and trees? {#alignment}

| Module | API | Behavior |
| --- | --- | --- |
| `alignment` | `MultipleAlignment.from_sequences`, `.distances` | Native equal-width, unique-ID alignment; p-distance with pairwise deletion or gaps as a state |
| `alignment` | `read_alignment`, `write_alignment` | Biopython single-alignment FASTA/Clustal/Stockholm/PHYLIP codecs |
| `alignment` | `align_sequences` | Python: `dnakit.bioinformatics.alignment.align_sequences`; Agent/MCP: `dnakit.alignment.align_sequences`; MAFFT `--auto`, MUSCLE **5**, Clustal Omega; restore IDs and verify sequence content |
| `alignment` | `build_blast_database`, `local_blast` | Local nucleotide database and BLAST XML output |
| `alignment` | `parse_search_results` | Streaming Bio.SearchIO QueryResult objects for BLAST XML/tab and protein HMMER2/3 formats |
| `alignment` | `build_profile`, `search_profile`, `parse_nhmmer_table` | hmmbuild, nucleotide nhmmer / protein hmmsearch; native nhmmer table parser |
| `phylogeny` | `TreeNode`, `build_tree` | Native neighbor joining and size-weighted UPGMA |
| `phylogeny` | `read_tree`, `write_tree`, `plot_tree` | Newick/Nexus/phyloXML/NeXML codecs and tree plots |
| `phylogeny` | `robinson_foulds`, `consensus_tree`, `prune_tree` | Native rooted/unrooted RF, strict-majority/strict consensus, pruning |

Remote BLAST already has explicit submission/polling APIs in `dnakit.search.blast`; this layer does not duplicate them.

NJ retains negative branch estimates and has an arbitrary root; UPGMA assumes a molecular clock. RF requires identical taxa and optionally divides by the combined split count. Consensus assumes the same root interpretation: threshold `0.5` keeps strictly majority clades, `1` keeps shared clades. GME/BME, rearrangement search and maximum-likelihood phylogenetics are not provided.

## What conventions apply to FASTX, coordinates and sequencing tools? {#reads}

| Module | API | Contract |
| --- | --- | --- |
| `reads` | `Read`, `read_fastx`, `write_fastx` | FASTA, multiline FASTQ, gzip, integer Phred scores; writer returns record count |
| `reads` | `pair_reads` | Two read iterators paired by IDs after removing `/1` or `/2`; ReadPair uses None for an orphan mate |
| `reads` | `repair_fastq` | Recover valid four-line records, report output/recovered count/discarded lines |
| `reads` | `arrange_reads`, `common_reads`, `rename_reads`, `sequence_digest` | Sorting/seeded shuffle, intersection, regex renaming, framed SHA256 |
| `reads` | `Region`, `read_regions` | BED or selected GTF feature rows; convert GTF coordinates exactly once |
| `reads` | `index_reference`, `IndexedReference.fetch` | Explicit `.fai` creation, context-managed indexed access, reverse complements |
| `reads` | `alignment_summary`, `coverage_consensus` | pysam decoding, native flag counts and A/C/G/T depth/majority consensus |
| `ngs` | `map_reads` | `dnakit.ngs.map_reads`: private BWA MEM / BWA-MEM2 / Bowtie2 reference index and SAM mapping |
| `ngs` | `prepare_bam` | `dnakit.ngs.prepare_bam`: samtools name-sort/fixmate/coordinate-sort/markdup/index/stats |
| `ngs` | `call_variants` | `dnakit.ngs.call_variants`: bcftools mpileup/call/norm and indexed VCF; explicit ploidy |
| `ngs` | `assemble_reads` | `dnakit.ngs.assemble_reads`: SPAdes single/paired short-read or Nanopore hybrid assembly |
| `ngs` | `quality_control`, `dna_seq_workflow` | Existing FastQC workflow; new bounded `dnakit.sequencing_qc.run_fastqc()` and `run_multiqc()` report facades; QC→BWA→BAM→variant calling with all intermediates retained |
| `ngs` | `rna_seq_workflow` | HISAT2→BAM→StringTie→featureCounts, strandedness 0/1/2 |
| `ngs` | `chip_seq_workflow` | Prepared treatment/control BAMs → MACS3 peaks and coverage |

FASTA→FASTQ requires explicit `default_quality` and does not infer quality. Repair recovers only valid four-line records and cannot guarantee unambiguous resynchronization of damaged files; use strict reading for multiline files. Pairing bounds pending read payload with `max_pending`; seen-ID sets still grow with total records.

Python coordinates are **0-based half-open**. BAM consensus is reference-forward regardless of `Region.strand`. It excludes secondary, supplementary, duplicate, QC-failed and low-quality reads; ties and insufficient depth produce N. Indels are retained in alignment/VCF data, not inserted into consensus.

Install specialist programs separately on PATH. `ExecutionConfig(threads=4, timeout_seconds=3600)` controls execution without a shell. Output directories must not exist; intermediates and diagnostics are retained. `max_output_bytes` bounds process logs; `max_artifact_bytes` bounds stdout or explicitly monitored output files, not the entire third-party working directory. Inputs are not overwritten; the explicit indexing API creates sidecar indexes next to its input.

The reference FASTA must already have a `.fai`; variant calling requires sorted, indexed BAM. The RNA workflow processes one sample; inference across biological replicates is separate. ChIP begins with prepared BAMs and composes with mapping/BAM APIs. There is no assembly-browser GUI.

## Which encodings and feature-analysis methods are supported? {#encoding}

`encoding.encode_dna({id: sequence}, descriptor, parameters=..., property_tables=..., labels=..., training_ids=...)` returns `EncodingResult(ids, names, values, descriptor, parameters, schema_version)`. `values` is a read-only samples×features array. Schema `dnakit.dna-encoding.v1` accepts canonical A/C/G/T and rejects unknown parameters.

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

All names are in `encoding.DNA_ENCODINGS`. Positional encodings require equal lengths. Physicochemical tables have the form `{property_name: {word: finite_value}}`, must cover all 16 dinucleotides or 64 trinucleotides, and are standardized using population SD over the complete word alphabet. No unattributed property tables are bundled. `NMBroto`, Moran, and Geary use dinucleotide tables. Short inputs, constant profiles and undefined denominators raise errors.

Pseudo type-I uses lagged mean squared differences across properties; SC uses per-property lagged product means. Subsequence weights are `delta ** span_length`, divided by total weight when normalized. MMI uses consistent marginals from same-length windows; pair/triplet cells sum to total MI/co-information. These explicit conventions **do not promise numerical identity with iFeatureOmega's similarly named vectors**. PSTNP/KNN read only the labels in `training_ids`, exclude each sample itself, and never use test labels.

`features.analyze_features` supports standard/minmax/L1/L2/none scaling; PCA, t-SNE, supervised LDA, Latent Dirichlet Allocation, chi2, F-classification/regression, MI, Information Gain (training quantile bins, `parameters={"bins": 10}`), Pearson, variance and none. Supervised methods require `labels`; chi2/topic models require nonnegative data. t-SNE has no held-out transform. Reuse the fitted training pipeline without refitting on held-out samples.

`cluster_features` supports kmeans, minibatch_kmeans, gmm, spectral, affinity_propagation, mean_shift, dbscan, agglomerative and birch, returning fitted sklearn estimators. Native `mcl` applies stochastic-matrix expansion/inflation to symmetric nonnegative affinities and returns index groups, raising on nonconvergence. Clustering does not automatically scale inputs.

## How should generative metrics be interpreted? {#metrics}

| `metrics` API | Definition |
| --- | --- |
| `kernel_distance` | Gaussian MMD² or cubic polynomial KID on embeddings; default unbiased estimates may be negative |
| `manifold_metrics` | Conventional PRDC precision/recall/density/coverage, open kNN balls, self excluded |
| `clipped_manifold_metrics` | Salvy et al. median-clipped density with leave-one-out normalization and beta-binomial calibrated coverage; closed balls and equal sample counts |
| `authenticity` | Fraction where generated-to-nearest-training distance exceeds that training point's nearest-other-training distance, in [0,1] |
| `pareto_ranks` | Explicit per-column `maximize`; zero-based nondomination ranks, duplicates share ranks |
| `hypervolume` | Explicit directions and a worse `reference_point`; native exact recursion, at most five dimensions/1000 points with work limit |
| `property_hit_rate` | Inclusive `lower/upper` bounds, joint/per-property rates and Boolean mask |
| `property_kl` | KL(generated∥reference), shared histogram edges, explicit pseudocount, natural logarithms |
| `conformity_score` | Cross-fitted isotropic Gaussian KDE and validation-density empirical scores |
| `sequence_properties` | Named scalar Python predictors → stable IDs/property names/matrix for the metrics above |

`authenticity` requires the **generator’s training set**; other reference sets are usually independent reference/validation sets. Fix embedding model, scaling and k before comparing methods. Conventional PRDC and clipped metrics differ; KID is not FID. Histogram property KL differs from seqme's KDE estimate. Conformity returns `value`, `standard_error`, and fold scores. It measures distributional consistency, not posterior class probability.

## What assumptions apply to statistical and count inference? {#statistics}

`dnakit.integrations.beta_diversity_scikit_bio()` and
`dnakit.integrations.permanova_scikit_bio()` are bounded optional scikit-bio
routes. They accept a finite metric allowlist and JSON numeric matrices; the
full scikit-bio statistics namespace is not re-exported.

`dnakit.integrations.fetch_pyfaidx_region()` provides read-only FASTA extraction
using zero-based half-open coordinates. `summarize_vcf_vcfpy()` and
`summarize_vcf_cyvcf2()` provide bounded VCF summaries. Each requires its
explicit optional backend and preserves backend provenance.

`dnakit.integrations.align_biotite()` provides bounded optimal nucleotide
pairwise alignment; `align_edlib()` provides edit distance, locations or a CIGAR
path; PyRanges adapters read zero-based half-open BED and summarize or merge
intervals. Inputs over `max_rows` are rejected instead of silently truncated;
`merge_pyranges_bed(slack=0)` merges adjacent intervals, and `stranded=True`
requires `+` or `-` on every row.

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

Permutation p-values use `(extreme+1)/(permutations+1)`; `permutations=0` returns `None`. No restricted/repeated-measures permutations. PERMDISP currently requires Euclidean distances. Composition transforms and ANCOM-BC require explicit zero handling; ANCOM accepts an explicit pseudocount. UniFrac requires a biologically rooted tree and nonnegative branches.

ANCOM-BC takes samples×taxa and a full-rank numeric design including an explicit intercept; at least ten taxa are required. It returns feature×covariate **natural-log** effects, SE, Wald, p/q, sampling bias and EM iteration counts. BH correction is per covariate. This is not ANCOM-BC2 and excludes structural-zero detection, global tests and bias-variance adjustment. Degenerate variance/nonconvergence raise errors.

NB2 takes integer samples×genes counts, a numeric design and explicit contrast, with at least two residual degrees of freedom and actual biological replicates. Default size factors are median ratios over genes positive in every sample; dispersion defaults to within-design residual moments. Explicit `size_factors` and `dispersions` are also accepted. Results include log2 fold changes, log2 SE, Wald p/q, dispersions and convergence flags; all-zero genes and failed fits have NaN inference. **DESeq2/edgeR empirical-Bayes trends, effect shrinkage and outlier strategies are not implemented.**

## What rules apply to optimization, motifs and cloning? {#optimization}

`optimization.specification(kind, region=None, weight=1, **parameters)` creates a native rule. Custom `Specification(name, evaluator)` accepts `(sequence, circular) -> SpecEvaluation`. `SpecEvaluation` separately represents nonnegative `violation`, `utility` to maximize, failure locations, and messages.

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

Structural rules are explicit exact-complement models, not full thermodynamic energy calculations. `blast_match_specification(database, minimum_identity=0.8, minimum_length=20)` uses specialized BLAST HSP detection with sequence caching. An empty search does not prove that no homology exists. `avoid_matches` is native exact k-mer matching; no Bowtie genomic-constraint wrapper is provided.

Hard constraints precede weighted objectives. Search uses greedy single-base substitutions or whole synonymous-codon replacements with `preserve_translation=True`. Length is fixed; `length` only checks it and does not insert/delete bases. Circular windows/patterns cross the origin and allow `end > length`. `search_exhausted` is not proof of infeasibility; `constraints_satisfied` does not establish a global optimum.

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

Duplex `offset` places the Crick 3′ end relative to the Watson 5′ end. Paired regions must be complementary; nicks, internal gaps and mismatches are unsupported. Gel rendering does not model dye/buffer dynamics or calibrated migration. Sanger consensus counts only quality-passing alignments meeting identity and coverage thresholds and reports N for ties/uncovered sites, without inserting indels. Trace-sample and base coordinates remain distinct.

## How can workflows and validation examples be run? {#workflows}

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
