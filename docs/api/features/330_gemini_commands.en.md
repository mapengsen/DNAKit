# GEMINI: queries, exports, statistics, inheritance models and gene-level tools

This group of functions in `dnakit.variants` runs sub-commands of [GEMINI](https://github.com/arq5x/gemini) on a GEMINI variant database (a SQLite file made by `dnakit.variants.load_gemini_database`): merging chunked loads, amending sample information, annotating with a BED/VCF, tagging somatic variants, carrier summaries by phenotype, table exports, region and gene look-ups, statistics (SNP substitution counts, site frequency spectrum, MDS distances, variants per sample, summaries), windowed Hardy-Weinberg or nucleotide diversity, five inheritance models (autosomal and X-linked), KEGG pathways, LoF prioritisation, gene burden and the C-alpha test, interaction partners, gene fusions, gene-wise aggregation and the sex check. Each function runs one fixed sub-command with an explicit argument list (no shell). Sub-commands that print a table stream it into a new file; sub-commands that change a database work on a copy unless `in_place=True`. Results carry `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `operation` and `provenance`; tables also get `row_count` and a short `rows` preview.

GEMINI 0.30.2 is the tested version of the upstream interface (MIT licence); DNAKit starts it as a separate program and copies none of its code. Install it with `conda install -c bioconda gemini` and put it on PATH or pass `gemini_path`; the annotation files that loading needs (about 20 GB, `gemini update --dataonly`) are the caller's responsibility and are never downloaded. If the program is missing, `BackendUnavailableError` is raised before anything is created. **Verification status:** the argument mapping of every function was written from the 0.30.2 command-line definitions (`gemini_main.py`) and is tested with a stand-in executable; it has not been run against a real GEMINI database because the GEMINI annotation data could not be installed here, so this page has no runnable examples with captured output.

## Example data

There are no example files on this page.

## 1) `dnakit.variants.gemini_amend()`

- **Purpose:** Amend the sample or pedigree information of a GEMINI database from a PED file (``gemini amend --sample``).
- **Method:** ``sample_file`` is a PED(+) file with new columns or values. ``clear=True`` first sets every value of the columns that the file provides to NULL. The change is made on a copy at ``output_database`` or, with ``in_place=True``, on ``database`` itself. Result key: ``database``.
- **API:** `gemini_amend(database: 'str', sample_file: 'str | None' = None, *, output_database: 'str | None' = None, in_place: 'bool' = False, clear: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 2) `dnakit.variants.gemini_annotate()`

- **Purpose:** Annotate the variants of a GEMINI database with a tabix-indexed BED or VCF file (``gemini annotate``).
- **Method:** ``annotation_file`` is a bgzipped, tabix-indexed BED (or VCF) and ``column_names`` the name(s) of the new variant table column(s) (for a VCF: the INFO field(s) to pull). ``annotation_type`` is ``boolean`` (overlap flag), ``count`` (number of overlapping records) or ``extract`` (GEMINI's default; values taken from ``extract_columns`` and combined by ``operations`` such as ``sum``, ``mean``, ``list``, ``first``). ``column_types`` gives ``integer``, ``float`` or ``text`` per column. For a VCF annotation, ``region_only`` ignores REF/ALT. Works on a copy at ``output_database`` or, with ``in_place=True``, on ``database`` itself.
- **API:** `gemini_annotate(database: 'str', annotation_file: 'str', column_names: 'str', *, output_database: 'str | None' = None, in_place: 'bool' = False, annotation_type: 'str' = 'extract', extract_columns: 'str | None' = None, column_types: 'str | None' = None, operations: 'str | None' = None, region_only: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 3) `dnakit.variants.gemini_burden()`

- **Purpose:** Gene-level variant burden counts per sample (or case/control group) (``gemini burden``).
- **Method:** By default every high-impact variant in coding regions counts towards a gene's burden; ``nonsynonymous=True`` (``--nonsynonymous``) counts all nonsynonymous variants. ``cases`` and ``controls`` are sample names; with both, GEMINI adds the burden test of cases against controls. ``min_aaf``/``max_aaf`` bound the alternate allele frequency of the variants used (GEMINI's defaults 0.0 and 1.0). The C-alpha test is ``gemini_calpha_test``. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_burden(database: 'str', output_path: 'str', *, nonsynonymous: 'bool' = False, cases: 'Sequence[str] | None' = None, controls: 'Sequence[str] | None' = None, min_aaf: 'float | None' = None, max_aaf: 'float | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 4) `dnakit.variants.gemini_calpha_test()`

- **Purpose:** Gene-level C-alpha association test between cases and controls (``gemini burden --calpha``).
- **Method:** ``cases`` and ``controls`` are sample names of the database (both required). ``permutations`` is the number of permutations for the permutation p-value (0: asymptotic only; try 1000 to start); ``save_tscores`` keeps the permuted T-scores in a file next to the working directory of GEMINI. The other options are those of ``gemini_burden``. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_calpha_test(database: 'str', output_path: 'str', cases: 'Sequence[str]', controls: 'Sequence[str]', *, permutations: 'int' = 0, nonsynonymous: 'bool' = False, min_aaf: 'float | None' = None, max_aaf: 'float | None' = None, save_tscores: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 5) `dnakit.variants.gemini_carrier_summary_by_phenotype()`

- **Purpose:** Counts of carriers and non-carriers of each variant, stratified by a sample phenotype column (``gemini query --carrier-summary-by-phenotype``).
- **Method:** ``phenotype_column`` is a column of the samples table (for example ``phenotype``); the query's variants get extra columns with the number of carriers and non-carriers per phenotype value. ``query`` is the SQL that selects the variants (GEMINI requires one; the default lists the basic columns), ``gt_filter`` restricts genotypes and ``region`` the positions (``chr1:10-20``). The TSV goes to ``output_path``. Result keys: ``columns``, ``row_count``, ``rows``.
- **API:** `gemini_carrier_summary_by_phenotype(database: 'str', phenotype_column: 'str', output_path: 'str', *, query: 'str' = 'SELECT chrom, start, end, ref, alt, gene FROM variants', gt_filter: 'str | None' = None, region: 'str | None' = None, header: 'bool' = True, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 6) `dnakit.variants.gemini_dump()`

- **Purpose:** Export a table of a GEMINI database as delimited text (``gemini dump``).
- **Method:** ``table`` selects the export: ``"variants"`` (``--variants``, every column of the variants table), ``"genotypes"`` (``--genotypes``, one line per variant and sample), ``"samples"`` (``--samples``, the samples/pedigree table) or ``"tfam"`` (``--tfam``, the sample information in TFAM format). ``header=True`` adds the column names; ``separator`` is the output column separator (GEMINI's default is a tab). Result keys: ``columns``, ``row_count``, ``rows``.
- **API:** `gemini_dump(database: 'str', output_path: 'str', *, table: 'str' = 'variants', header: 'bool' = False, separator: 'str' = '\t', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 7) `dnakit.variants.gemini_fusions()`

- **Purpose:** Detect somatic gene fusions from structural-variant breakends in a GEMINI database (``gemini fusions``).
- **Method:** ``in_cosmic_census=True`` keeps fusions where one or both genes are in the COSMIC cancer census, ``min_quality`` is the minimum VCF QUAL and ``evidence_type`` the supporting evidence (``PE``, ``SR`` or ``PE,SR``). Result keys: ``row_count``, ``rows``.
- **API:** `gemini_fusions(database: 'str', output_path: 'str', *, in_cosmic_census: 'bool' = False, min_quality: 'float | None' = None, evidence_type: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 8) `dnakit.variants.gemini_gene_wise()`

- **Purpose:** Aggregate variants per gene with genotype filters, reporting genes that pass enough filters (``gemini gene_wise``).
- **Method:** ``genotype_filters`` (``--gt-filter``) are genotype expressions evaluated per variant; a gene is reported when at least ``min_filters`` of them (default 1) pass among its variants. ``required_genotype_filters`` (``--gt-filter-required``) must pass for a variant to count at all. ``where`` is the SQL WHERE clause that selects the variants (GEMINI's default is ``is_exonic = 1 AND impact_severity != 'LOW'``), ``filter`` an extra SQL restriction and ``columns`` the columns to return (default ``chrom,start,end,gene,impact,impact_severity,max_aaf_all``). Result keys: ``row_count``, ``rows``.
- **API:** `gemini_gene_wise(database: 'str', output_path: 'str', *, genotype_filters: 'Sequence[str] | None' = None, required_genotype_filters: 'Sequence[str] | None' = None, min_filters: 'int | None' = None, where: 'str | None' = None, filter: 'str | None' = None, columns: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 9) `dnakit.variants.gemini_inheritance_model()`

- **Purpose:** Find candidate variants that fit a Mendelian inheritance model in the families of a GEMINI database (``gemini autosomal_recessive`` and friends).
- **Method:** ``model`` is ``"autosomal_recessive"``, ``"autosomal_dominant"``, ``"x_linked_recessive"``, ``"x_linked_dominant"`` or ``"x_linked_de_novo"`` (the GEMINI sub-command of the same name). The pedigree comes from the database. ``columns`` selects the returned columns (default ``*``), ``filter`` restricts the variants (SQL), ``min_kindreds`` is the minimum number of families with a candidate in a gene (default 1), ``families`` a comma-separated list of family ids, ``allow_unaffected`` reports candidates that also affect unaffected samples, ``min_sample_depth`` (``-d``) and ``min_gq`` set per-sample depth and genotype-quality minima (default 0). Autosomal models also take ``lenient`` (looser family structure) and ``max_phred_likelihood`` (``--gt-pl-max``); X-linked models take ``x_chromosomes`` (``-X``, names of the X chromosome when not ``chrX``/``X``). Options of the other kind raise ConfigurationError. Result keys: ``model``, ``row_count``, ``rows``.
- **API:** `gemini_inheritance_model(database: 'str', output_path: 'str', *, model: 'str' = 'autosomal_recessive', columns: 'str | None' = None, filter: 'str | None' = None, min_kindreds: 'int | None' = None, families: 'str | None' = None, allow_unaffected: 'bool' = False, min_sample_depth: 'int | None' = None, min_gq: 'int | None' = None, lenient: 'bool' = False, max_phred_likelihood: 'int | None' = None, x_chromosomes: 'Sequence[str] | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 10) `dnakit.variants.gemini_interactions()`

- **Purpose:** Gene interaction partners of the genes that carry variants in the samples (``gemini interactions`` and ``gemini lof_interactions``).
- **Method:** ``lof=False`` uses ``gemini interactions``: with ``gene`` (``-g``) as the root and ``radius`` (``-r``) the breadth-first search depth (0 and up). ``lof=True`` uses ``gemini lof_interactions`` for genes with loss-of-function variants (``gene`` does not apply). ``edges_file`` replaces the default HPRD network with ``geneA|geneB`` lines. ``variant_mode`` (``--var``) also returns variant information (impact, biotype) for the interacting genes. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_interactions(database: 'str', output_path: 'str', *, lof: 'bool' = False, gene: 'str | None' = None, radius: 'int | None' = None, edges_file: 'str | None' = None, variant_mode: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 11) `dnakit.variants.gemini_lof_sieve()`

- **Purpose:** Prioritise loss-of-function variants by the likelihood that they are real LoF (``gemini lof_sieve``).
- **Method:** Reports each LoF variant with the transcript position, the fraction of the transcript affected and whether the transcript is canonical, which separate likely true LoF from variants near the 3' end. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_lof_sieve(database: 'str', output_path: 'str', *, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 12) `dnakit.variants.gemini_merge_chunks()`

- **Purpose:** Merge the partial databases written by a chunked parallel load into one database (``gemini merge_chunks``).
- **Method:** ``chunk_databases`` are the chunk databases (``gemini load_chunk`` output of one VCF). The merged database is the new file ``output_database``. ``vcf_path`` is the original VCF (GEMINI reads extra annotation fields from it) and ``anno_type`` one of ``snpEff``, ``VEP``, ``BCFT`` or ``all``. ``index=True`` builds all database indexes (do it only in the last of several merges); ``tempdir`` is a local, non-NFS directory for SQLite working files.
- **API:** `gemini_merge_chunks(chunk_databases: 'Sequence[str]', output_database: 'str', *, vcf_path: 'str | None' = None, anno_type: 'str | None' = None, index: 'bool' = False, skip_pls: 'bool' = False, tempdir: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 13) `dnakit.variants.gemini_pathways()`

- **Purpose:** Map the genes and variants of a GEMINI database to KEGG pathways (``gemini pathways``).
- **Method:** ``ensembl_version`` is the Ensembl gene version (GEMINI supports 66 to 71; default 68). ``lof=True`` (``--lof``) reports the pathways of individuals, genes and sites with loss-of-function variants. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_pathways(database: 'str', output_path: 'str', *, ensembl_version: 'str' = '68', lof: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 14) `dnakit.variants.gemini_region()`

- **Purpose:** List the variants in a gene (by symbol) or a genomic region (``gemini region --gene`` / ``--reg``).
- **Method:** Give exactly one of ``gene`` (a gene symbol such as ``TP53``) or ``region`` (``chr1:100-200``). ``columns`` is the comma-separated list of columns to return (GEMINI's default is ``*``), ``filter`` an SQL restriction on the variants, ``show_samples`` adds the names of the samples that carry each variant and ``output_format`` is ``JSON``, ``TPED`` or ``default``. Result keys: ``columns``, ``row_count``, ``rows``.
- **API:** `gemini_region(database: 'str', output_path: 'str', *, gene: 'str | None' = None, region: 'str | None' = None, columns: 'str | None' = None, filter: 'str | None' = None, header: 'bool' = False, show_samples: 'bool' = False, output_format: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 15) `dnakit.variants.gemini_set_somatic()`

- **Purpose:** Tag somatic variants (``is_somatic``) by comparing tumor and normal samples of a GEMINI database (``gemini set_somatic``).
- **Method:** The thresholds are GEMINI's options of the same name (all default to 0 in GEMINI): combined tumor+normal depth, VCF QUAL, somatic score, maximum alternate frequency/count in the normal, minimum depth in the normal and the minimum alternate frequency/count/depth in the tumor. ``chrom`` restricts the tagging to one chromosome; ``dry_run=True`` only reports what would be set (the database is not changed, but is still copied unless ``in_place``). The tumor/normal pairs come from the phenotype or PED information stored in the database.
- **API:** `gemini_set_somatic(database: 'str', *, output_database: 'str | None' = None, in_place: 'bool' = False, min_depth: 'float | None' = None, min_quality: 'float | None' = None, min_somatic_score: 'float | None' = None, max_normal_alt_freq: 'float | None' = None, max_normal_alt_count: 'int | None' = None, min_normal_depth: 'int | None' = None, min_tumor_alt_freq: 'float | None' = None, min_tumor_alt_count: 'int | None' = None, min_tumor_depth: 'int | None' = None, chrom: 'str | None' = None, dry_run: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 16) `dnakit.variants.gemini_sex_check()`

- **Purpose:** Check the recorded sex of each sample against its genotypes on the X chromosome (``gemini qc --mode sex``).
- **Method:** ``chromosome`` is the name of the X chromosome in the database (GEMINI's default ``chrX``; use ``X`` for Ensembl-style names). The table has one row per sample with the sex in the pedigree, the number of homozygous/ heterozygous sites and the inferred sex. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_sex_check(database: 'str', output_path: 'str', *, chromosome: 'str' = 'chrX', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 17) `dnakit.variants.gemini_stats()`

- **Purpose:** Basic statistics of a GEMINI database (``gemini stats``).
- **Method:** ``statistic`` selects the report: ``"snp_counts"`` (``--snp-counts``, count of each substitution type A->G, G->T, ...), ``"sfs"`` (``--sfs``, site frequency spectrum), ``"mds"`` (``--mds``, pairwise genetic distance between samples, i.e. multidimensional-scaling input), ``"vars_by_sample"`` (``--vars-by-sample``, variants per sample), ``"gts_by_sample"`` (``--gts-by-sample``, genotype-class counts per sample), ``"tstv"``, ``"tstv_coding"`` and ``"tstv_noncoding"`` (transition/transversion ratios) and ``"summarize"`` (``--summarize``, summarise the variants of your SQL ``query`` by sample). ``gt_filter`` restricts genotypes. Result keys: ``columns`` (None), ``row_count``, ``rows``.
- **API:** `gemini_stats(database: 'str', output_path: 'str', *, statistic: 'str' = 'snp_counts', query: 'str | None' = None, gt_filter: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).

## 18) `dnakit.variants.gemini_windower()`

- **Purpose:** Windowed population-genetic statistic along the genome: Hardy-Weinberg or nucleotide diversity (``gemini windower``).
- **Method:** ``analysis`` is ``"hwe"`` (``-t hwe``, windowed Hardy-Weinberg equilibrium test; GEMINI's default) or ``"nucl_div"`` (``-t nucl_div``, windowed nucleotide diversity). ``window_size`` is the window length in bp (default 1,000,000), ``step_size`` the step in bp (0 lets GEMINI use the window size) and ``operation`` how the per-variant values of a window are combined: ``mean`` (default), ``median``, ``min``, ``max`` or ``collapse``. Result keys: ``row_count``, ``rows``.
- **API:** `gemini_windower(database: 'str', output_path: 'str', *, analysis: 'str' = 'hwe', window_size: 'int' = 1000000, step_size: 'int' = 0, operation: 'str' = 'mean', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the GEMINI program (`conda install -c bioconda gemini`) and an existing GEMINI database; DNAKit looks for it at `gemini_path`, then on PATH.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
No runnable example: GEMINI and its annotation data are not installed in the documentation environment (see the verification status above).
