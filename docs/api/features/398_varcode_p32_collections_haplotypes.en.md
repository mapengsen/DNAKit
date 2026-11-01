# varcode collections, phasing and haplotypes

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: Known-cis haplotype annotation requires at least two variants. Transcript-model prediction retains unresolved outcomes. Local GTF/FASTA files and cache_dir are required for reference-dependent calls; no reference data are downloaded.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_read_variant_csv`

**Purpose:** Read a variant collection from a CSV file written by varcode (``VariantCollection.from_csv``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_read_variant_csv(gtf: str, genome_fasta: str, csv_path: str, *, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read a variant collection from a CSV file written by varcode (``VariantCollection.from_csv``).

The CSV has the columns ``chr,start,ref,alt`` (optionally gene columns) after ``#`` comment lines. The variants
are attached to the local reference built from ``gtf`` and ``genome_fasta``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_variant_table`

**Purpose:** Table of a variant collection with the genes each variant overlaps, optionally written to a CSV file.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_variant_table(gtf: str, genome_fasta: str, *, vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, output_csv: str | None=None, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Table of a variant collection with the genes each variant overlaps, optionally written to a CSV file.

Without ``output_csv`` the table comes from ``VariantCollection.to_dataframe`` (columns chr, start, ref, alt,
gene_name, gene_id), which is also the listing of the ``varcode-genes`` command. With ``output_csv`` the collection
is also written by ``VariantCollection.to_csv`` (an existing file is refused). Sources are VCF, MAF, varcode
CSV and JSON files and explicit records.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_variant_set_operation`

**Purpose:** Combine two variant call sets: ``operation`` is ``union``, ``intersection`` or ``difference`` (first minus second).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_variant_set_operation(gtf: str, genome_fasta: str, operation: str, first: Mapping[str, Any], second: Mapping[str, Any], *, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Combine two variant call sets: ``operation`` is ``union``, ``intersection`` or ``difference`` (first minus second).

``first`` and ``second`` are mappings with any of the keys ``vcf``, ``maf``, ``csv``, ``json`` (lists of file
paths) and ``variants`` (list of records). They map to ``VariantCollection.union``, ``.intersection`` and
``.difference``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_filter_by_expression`

**Purpose:** Keep the variants or effects whose transcript or gene is expressed above a threshold.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_filter_by_expression(gtf: str, genome_fasta: str, expression: Mapping[str, float], *, level: str='transcript', target: str='variants', min_expression_value: float=0.0, vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, skip_errors: bool=False, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Keep the variants or effects whose transcript or gene is expressed above a threshold.

``expression`` maps transcript ids (``level='transcript'``, ``filter_by_transcript_expression``) or gene ids
(``level='gene'``, ``filter_by_gene_expression``) to a level such as TPM. ``target='variants'`` filters the
``VariantCollection``; ``target='effects'`` annotates the variants first and filters the ``EffectCollection``.
Items with a level at or below ``min_expression_value`` are removed.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_variants_in_sample`

**Purpose:** Variants of a VCF sample, optionally restricted by zygosity.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_variants_in_sample(gtf: str, genome_fasta: str, vcf_paths: Sequence[str], sample: str, *, genotype: str='carried', only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Variants of a VCF sample, optionally restricted by zygosity.

``genotype='carried'`` returns the variants the sample carries (``for_sample``), ``'heterozygous'`` those called
heterozygous (``heterozygous_in``) and ``'homozygous_alt'`` those homozygous for the alternate allele
(``homozygous_alt_in``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_sample_genotype`

**Purpose:** Genotype and zygosity of one sample at one variant of a VCF (``VariantCollection.genotype`` and ``.zygosity``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_sample_genotype(gtf: str, genome_fasta: str, vcf_paths: Sequence[str], sample: str, variant: Mapping[str, Any], *, alt_index: int=1, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Genotype and zygosity of one sample at one variant of a VCF (``VariantCollection.genotype`` and ``.zygosity``).

The result holds the call (GT, alleles, phase, phase set, depths) and, for allele ``alt_index`` (1 = first ALT),
the zygosity (``Genotype.zygosity_for_alt``), whether the call carries it, its copy number and its read depth.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_genotype_alt_summary`

**Purpose:** Summarise one VCF genotype for an alternate allele without any reference (``varcode.Genotype``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_genotype_alt_summary(gt: str, *, allele_depths: Sequence[int] | None=None, alt_index: int=1, phase_set: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Summarise one VCF genotype for an alternate allele without any reference (``varcode.Genotype``).

``gt`` is a VCF GT string such as ``0/1`` or ``1|1``; ``allele_depths`` is the AD list. Returns
``carries_alt`` (``Genotype.carries_alt``), ``copies_of_alt`` (``copies_of_alt``), ``depth_for_alt``
(``depth_for_alt``, read depth of the alt allele from AD) and ``zygosity`` (``zygosity_for_alt``) for allele
``alt_index`` (1 = first ALT), plus ploidy and phase.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_effect_table`

**Purpose:** Annotate variants and return the effects as a table (``EffectCollection.to_dataframe``), optionally also as CSV.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_effect_table(gtf: str, genome_fasta: str, *, vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, skip_errors: bool=False, output_csv: str | None=None, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Annotate variants and return the effects as a table (``EffectCollection.to_dataframe``), optionally also as CSV.

Columns: variant, contig, start, ref, alt, gene and transcript ids and names, effect_type and effect description.
With ``output_csv`` the effects are also written by ``EffectCollection.to_csv`` (an existing file is refused).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_read_effect_csv`

**Purpose:** Read an effect CSV written by varcode and re-annotate its variants on the local reference (``EffectCollection.from_csv``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_read_effect_csv(gtf: str, genome_fasta: str, csv_path: str, *, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read an effect CSV written by varcode and re-annotate its variants on the local reference (``EffectCollection.from_csv``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_effect_priority`

**Purpose:** Consequence severity rank of varcode effect classes (``varcode.effect_priority``): a higher number is more severe.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_effect_priority(effect_classes: Sequence[str] | None=None, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Consequence severity rank of varcode effect classes (``varcode.effect_priority``): a higher number is more severe.

``effect_classes`` lists class names such as ``Substitution`` or ``PrematureStop``; without it every ranked class
is returned, most severe first.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_top_priority_effects`

**Purpose:** Highest-priority (most severe) effect of annotated variants.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_top_priority_effects(gtf: str, genome_fasta: str, *, scope: str='variant', vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, skip_errors: bool=False, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Highest-priority (most severe) effect of annotated variants.

``scope='overall'`` gives the single most severe effect of the whole collection (``top_priority_effect``; give one
variant for "the top effect of a variant"); ``'variant'`` one effect per variant
(``top_priority_effect_per_variant``); ``'transcript'`` one per transcript id
(``top_priority_effect_per_transcript_id``); ``'gene'`` one per gene id (``top_priority_effect_per_gene_id``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_filter_effects`

**Purpose:** Filter annotated effects by consequence.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_filter_effects(gtf: str, genome_fasta: str, *, drop_silent_and_noncoding: bool=True, keep_unresolved: bool=True, min_priority_class: str | None=None, vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, skip_errors: bool=False, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Filter annotated effects by consequence.

``drop_silent_and_noncoding`` removes silent and non-coding effects (``drop_silent_and_noncoding``, the
``varcode --only-coding`` option; ``keep_unresolved`` keeps effects whose coding status is unknown).
``min_priority_class`` keeps effects at least as severe as the named class, for example
``NonsynonymousCodingMutation`` or ``Substitution`` (``filter_by_effect_priority``). Both filters can be combined.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_expression_effects`

**Purpose:** Use transcript expression to rank the effects of annotated variants.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_expression_effects(gtf: str, genome_fasta: str, expression: Mapping[str, float], *, mode: str='annotate', vcf_paths: Sequence[str]=(), maf_paths: Sequence[str]=(), csv_paths: Sequence[str]=(), json_paths: Sequence[str]=(), variants: Sequence[Mapping[str, Any]]=(), only_passing: bool=True, skip_errors: bool=False, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Use transcript expression to rank the effects of annotated variants.

``expression`` maps transcript ids to a level such as FPKM or TPM. ``mode='annotate'`` attaches the level to every
effect that has a transcript (``EffectCollection.effect_expression``; effects without a transcript are left
out, transcripts not listed count as 0; splice-outcome sets are included, which varcode 10.11.10 itself cannot hash). ``mode='top'`` returns the one effect on the most highly expressed
transcript, ties broken by effect priority (``top_expression_effect``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_reconcile_splice_outcomes_with_rna`

**Purpose:** Reconcile the predicted splice mechanisms of a variant with RNA-observed mechanisms.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_reconcile_splice_outcomes_with_rna(gtf: str, genome_fasta: str, variant: Mapping[str, Any], transcript_id: str, observed: Sequence[Mapping[str, Any]], *, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Reconcile the predicted splice mechanisms of a variant with RNA-observed mechanisms.

The splice outcomes of the variant on the transcript are enumerated (``enumerate_splice_outcomes``). Each entry of
``observed`` names a predicted mechanism class (``mechanism``, for example ``IntronRetention``,
``ExonSkipping`` or ``CrypticDonor``) and may give ``read_count`` and ``source``. The set is then reconciled
with ``SpliceOutcomeSet.with_rna_evidence``: observed mechanisms replace their broad predictions, unobserved
predictions are listed as excluded, and ``SpliceOutcomeSet.rna_evidence_for`` gives the evidence kept for each
remaining candidate. Variants that do not disrupt a splice signal are refused.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_query_in_cis`

**Purpose:** Are two variants on the same haplotype (cis), on different ones (trans), or unknown (None)?

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_query_in_cis(gtf: str, genome_fasta: str, variant_a: Mapping[str, Any], variant_b: Mapping[str, Any], *, method: str='vcf', vcf_path: str | None=None, sample: str | None=None, bam_path: str | None=None, transcript_id: str | None=None, min_mapping_quality: int=20, min_base_quality: int=20, min_alt_reads: int=2, phasing_error_rate: float=0.01, max_p_value_for_phasing: float=0.05, require_proper_pair: bool=True, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Are two variants on the same haplotype (cis), on different ones (trans), or unknown (None)?

``method='vcf'`` reads the genotype and phase set of ``sample`` in ``vcf_path`` (``VCFPhaseResolver.in_cis``,
GT with ``|`` and PS). ``method='rna_reads'`` counts RNA fragments from ``bam_path`` that carry both alleles
(``RNAReadPhasingSource.in_cis`` through ``MolecularPhaseResolver``); the read filters are the varcode defaults.
Either resolver is queried through ``varcode.query_in_cis``. ``transcript_id`` limits the question to one transcript.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_phased_partners`

**Purpose:** Variants observed in cis with a variant.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_phased_partners(gtf: str, genome_fasta: str, variant: Mapping[str, Any], *, method: str='vcf', vcf_path: str | None=None, sample: str | None=None, bam_path: str | None=None, candidates: Sequence[Mapping[str, Any]]=(), transcript_id: str | None=None, min_mapping_quality: int=20, min_base_quality: int=20, min_alt_reads: int=2, phasing_error_rate: float=0.01, max_p_value_for_phasing: float=0.05, require_proper_pair: bool=True, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Variants observed in cis with a variant.

``method='vcf'`` lists the variants in the same VCF phase set and haplotype of ``sample``
(``VCFPhaseResolver.phased_partners``). ``method='rna_reads'`` lists the ``candidates`` that are co-observed with
the variant on the same RNA read or fragment of ``bam_path`` (``MolecularPhaseResolver.phased_partners`` and
``RNAReadPhasingSource.partners_in_cis``); an empty list means no direct evidence.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_rna_read_support`

**Purpose:** Count the quality-filtered RNA reads or fragments carrying each variant's alt allele.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_rna_read_support(gtf: str, genome_fasta: str, bam_path: str, variants: Sequence[Mapping[str, Any]], *, min_mapping_quality: int=20, min_base_quality: int=20, min_alt_reads: int=2, phasing_error_rate: float=0.01, max_p_value_for_phasing: float=0.05, require_proper_pair: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Count the quality-filtered RNA reads or fragments carrying each variant's alt allele.

For each variant returns ``alt_reads`` (``RNAReadPhasingSource.supports_variant``; None when the contig is not in the
BAM) and ``has_evidence`` (``has_evidence``: at least ``min_alt_reads`` alt fragments). The BAM must be
coordinate-sorted and indexed.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_test_rna_haplotype`

**Purpose:** Test whether RNA reads carry a known local combination of alleles on one alignment (``RNAReadPhasingSource.register_haplotype``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_test_rna_haplotype(gtf: str, genome_fasta: str, bam_path: str, variants: Sequence[Mapping[str, Any]], *, flanking_bases: int=5, min_mapping_quality: int=20, min_base_quality: int=20, min_alt_reads: int=2, phasing_error_rate: float=0.01, max_p_value_for_phasing: float=0.05, require_proper_pair: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Test whether RNA reads carry a known local combination of alleles on one alignment (``RNAReadPhasingSource.register_haplotype``).

``variants`` are non-overlapping substitutions or deletions on one contig that form one hypothesis. A read supports it
only when its whole sequence, including ``flanking_bases`` unchanged reference bases on each side, equals the
hypothesis. Returns per-variant support after registration and, for each pair, the cis/trans call.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_load_germline_context`

**Purpose:** Load a germline VCF into a varcode ``GermlineContext``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_load_germline_context(gtf: str, genome_fasta: str, vcf_path: str, *, sample: str | None=None, completeness: str='complete', only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Load a germline VCF into a varcode ``GermlineContext``.

Without ``sample`` the whole VCF is the germline call set (``GermlineContext.from_germline_vcf``). With ``sample`` that
sample's non-reference calls are extracted from a multi-sample VCF (``from_multi_sample_vcf``). ``completeness`` states
what a missing call means: ``complete`` (absent means reference), ``sparse``, ``hotspots_only`` or ``empty``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_germline_variants`

**Purpose:** Germline variants of a context inside a window, or the homozygous ones.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_germline_variants(gtf: str, genome_fasta: str, vcf_path: str, *, contig: str | None=None, start: int | None=None, end: int | None=None, homozygous_only: bool=False, sample: str | None=None, completeness: str='complete', only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Germline variants of a context inside a window, or the homozygous ones.

With ``contig``, ``start`` and ``end`` (1-based, inclusive) the variants overlapping that window are returned
(``GermlineContext.variants_in_window``). With ``homozygous_only`` the result is limited to the variants the sample carries
on both haplotypes (``GermlineContext.homozygous``); this needs a single-sample VCF, because varcode 10.11.10
keeps no genotypes when it extracts one sample from a multi-sample VCF (``sample`` is refused with
``homozygous_only``). Without a window the whole context is considered.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_validate_germline_context`

**Purpose:** Check a germline context against somatic variants before annotation (``GermlineContext.validate_against``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_validate_germline_context(gtf: str, genome_fasta: str, germline_vcf: str, somatic: Mapping[str, Any], *, sample: str | None=None, completeness: str='complete', validate_reference: bool=True, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Check a germline context against somatic variants before annotation (``GermlineContext.validate_against``).

``somatic`` is a source mapping (keys ``vcf``, ``maf``, ``csv``, ``json``, ``variants``). The check fails on a
reference-build mismatch (unless ``validate_reference`` is False) and warns about empty or suspicious call sets.
Returns ``valid`` with the error text when it fails, and the warnings varcode raised.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_functions.varcode_sv_junction_inserted_sequence`

**Purpose:** Bases inserted at each junction of the structural variants in VCF files (``StructuralVariant.junction_inserted_sequence``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_sv_junction_inserted_sequence(gtf: str, genome_fasta: str, vcf_paths: Sequence[str], *, only_passing: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Bases inserted at each junction of the structural variants in VCF files (``StructuralVariant.junction_inserted_sequence``).

The sequence is oriented from the retained first breakend towards the retained second one. An empty string means no
insertion; None means the allele is conflicting or unreadable. Small variants in the VCF are skipped.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_haplotype_functions.varcode_annotate_cis_haplotype`

**Purpose:** Annotate the joint effect of known-cis point edits with FastEffectAnnotator.annotate_haplotype.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_annotate_cis_haplotype(gtf: str, genome_fasta: str, variants: Sequence[Mapping[str, Any]], transcript_id: str, *, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Annotate the joint effect of known-cis point edits with FastEffectAnnotator.annotate_haplotype.

At least two variants are required; each has contig, 1-based start, ref and alt. Structural variants and
unsupported transcript models raise ConfigurationError, matching varcode's
explicit NotImplemented outcome. No germline or RNA evidence is inferred.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.varcode_p32_haplotype_functions.varcode_predict_transcript_model`

**Purpose:** Predict transcript-model effects with splice alternatives via predict_transcript_model_effect.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def varcode_predict_transcript_model(gtf: str, genome_fasta: str, variants: Sequence[Mapping[str, Any]], transcript_id: str, *, germline_variants: Sequence[Mapping[str, Any]]=(), max_hypotheses: int=64, use_genomic_sequence: bool=True, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Predict transcript-model effects with splice alternatives via predict_transcript_model_effect.

Known-cis variants are evaluated jointly; germline_variants retain upstream
phase uncertainty. Genomic sequence enables intron-retention/cryptic-splice
realisation. HypothesisLimit and Unresolved remain explicit effect classes.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.molecular_lab.varcode_p32_functions import varcode_genotype_alt_summary
r = varcode_genotype_alt_summary('0/1', allele_depths=[3, 7], verbose=False)
print(r['genotype']['carries_alt'], r['genotype']['copies_of_alt'], r['genotype']['depth_for_alt'])
```

**Example output:**

```text
True 1 7
```
