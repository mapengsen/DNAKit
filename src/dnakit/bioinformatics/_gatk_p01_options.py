"""Long options that ``gatk <Tool> --help`` lists for the tools of pack P01 (GATK 4.7.0.0).

Generated once from the help text of each tool. ``arguments_file``, ``gatk-config-file``, the help and version
switches, hidden-argument display and ``JAVASCRIPT_FILE`` (a script that GATK would execute) are left out on purpose.
"""

from __future__ import annotations

TOOL_OPTIONS: dict[str, frozenset[str]] = {
    "AddCommentsToBam": frozenset(
        (
            "COMMENT COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "AddOATag": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT INTERVAL_LIST MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "AddOrReplaceReadGroups": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE RGCN RGDS RGDT RGFO RGID RGKS RGLB RGPG RGPI RGPL RGPM RGPU RGSM SORT_ORDER "
            "TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "AddOriginalAlignmentTags": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "ApplyBQSR": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line allow-missing-read-group "
            "ambig-filter-bases ambig-filter-frac black-listed-lanes bqsr-recal-file cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends emit-original-quals exclude-intervals filter-too-short "
            "gcs-max-retries global-qscore-prior input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name preserve-qscores-less-than "
            "quantize-quals read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference round-down-quantized sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output static-quantized-quals tmp-dir use-jdk-deflater "
            "use-jdk-inflater use-original-qualities verbosity"
        ).split()
    ),
    "BQSRPipelineSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line allow-missing-read-group ambig-filter-bases ambig-filter-frac "
            "bam-partition-size binary-tag-name black-listed-lanes bqsr-baq-gap-open-penalty conf "
            "create-output-bam-index create-output-bam-splitting-index create-output-variant-index "
            "default-base-qualities deletions-default-quality disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends emit-original-quals exclude-intervals filter-too-short "
            "gcs-max-retries global-qscore-prior indels-context-size input insertions-default-quality "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only known-sites library "
            "low-quality-tail mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio maximum-cycle-value "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality "
            "mismatches-context-size mismatches-default-quality num-reducers output output-cram-version "
            "output-shard-tmp-dir platform-filter-name preserve-qscores-less-than program-name quantize-quals "
            "quantizing-levels read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference round-down-quantized sample sharded-output spark-master "
            "spark-verbosity splitting-index-granularity static-quantized-quals tmp-dir use-jdk-deflater "
            "use-jdk-inflater use-nio use-original-qualities verbosity"
        ).split()
    ),
    "BaitDesigner": frozenset(
        (
            "BAIT_OFFSET BAIT_SIZE COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE DESIGN_NAME "
            "DESIGN_ON_TARGET_STRAND DESIGN_STRATEGY FILL_POOLS LEFT_PRIMER MAX_RECORDS_IN_RAM "
            "MERGE_NEARBY_TARGETS MINIMUM_BAITS_PER_TARGET OUTPUT_AGILENT_FILES OUTPUT_DIRECTORY PADDING "
            "POOL_SIZE QUIET REFERENCE_SEQUENCE REPEAT_TOLERANCE RIGHT_PRIMER TARGETS TMP_DIR USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "BamToBfq": frozenset(
        (
            "ANALYSIS_DIR BASES_TO_WRITE CLIP_ADAPTERS COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "FLOWCELL_BARCODE INCLUDE_NON_PF_READS INPUT LANE MAX_RECORDS_IN_RAM OUTPUT_FILE_PREFIX PAIRED_RUN "
            "QUIET READS_TO_ALIGN READ_CHUNK_SIZE READ_NAME_PREFIX REFERENCE_SEQUENCE RUN_BARCODE TMP_DIR "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "BaseRecalibrator": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "binary-tag-name black-listed-lanes bqsr-baq-gap-open-penalty cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 default-base-qualities deletions-default-quality disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries "
            "indels-context-size input insertions-default-quality interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only known-sites lenient library low-quality-tail "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-cycle-value maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality mismatches-context-size mismatches-default-quality output "
            "output-cram-version platform-filter-name preserve-qscores-less-than quantizing-levels read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater use-original-qualities verbosity"
        ).split()
    ),
    "BwaAndMarkDuplicatesPipelineSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes bwa-mem-index-image conf create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index disable-read-filter "
            "disable-tool-default-read-filters do-not-mark-unmapped-mates dont-require-soft-clips-both-ends "
            "duplicate-scoring-strategy duplicate-tagging-policy exclude-intervals filter-too-short "
            "flow-end-pos-uncertainty flow-q-is-known-end flow-quality-sum-strategy flow-skip-start-homopolymers "
            "flowbased gcs-max-retries input interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only library mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality num-reducers "
            "output output-cram-version output-shard-tmp-dir platform-filter-name program-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference "
            "remove-all-duplicates remove-sequencing-duplicates sample sharded-output single-end-alignment "
            "single-end-reads-clipping-is-end single-end-reads-end-position-significant spark-master "
            "spark-verbosity splitting-index-granularity tmp-dir use-jdk-deflater use-jdk-inflater use-nio "
            "verbosity"
        ).split()
    ),
    "BwaMemIndexImageCreator": frozenset(
        (
            "QUIET gcs-max-retries input output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "BwaSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes bwa-mem-index-image conf create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index disable-read-filter "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only library mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality num-reducers "
            "output output-cram-version output-shard-tmp-dir platform-filter-name program-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "sharded-output single-end-alignment spark-master spark-verbosity splitting-index-granularity tmp-dir "
            "use-jdk-deflater use-jdk-inflater use-nio verbosity"
        ).split()
    ),
    "CRAMIssue8768Detector": frozenset(
        (
            "QUIET echo-to-stdout gcs-max-retries input mismatch-rate-threshold output output-tsv reference "
            "tmp-dir use-jdk-deflater use-jdk-inflater verbose verbosity"
        ).split()
    ),
    "CalculateGenotypePosteriors": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 de-novo-prior "
            "default-to-allele-count disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "discovered-allele-count-priors-off dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries global-prior-indel global-prior-snp ignore-input-samples input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality "
            "num-reference-samples-if-no-call output output-cram-version pedigree platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-index read-name read-validation-stringency reference "
            "sample seconds-between-progress-updates sequence-dictionary sites-only-vcf-output skip-family-priors "
            "skip-population-priors supporting-callsets tmp-dir use-flat-priors-for-indels use-jdk-deflater "
            "use-jdk-inflater variant variant-output-filtering verbosity"
        ).split()
    ),
    "CalibrateDragstrModel": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "api-mono-threshold api-values black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "debug-sites-output disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends down-sample-size exclude-intervals filter-too-short "
            "force-estimation gcs-max-retries gop-values gp-values het-to-hom-ratio input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-period max-read-length max-repeats "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-loci-count min-read-length minimum-depth "
            "minimum-mapping-quality output output-cram-version parallel pileup-padding platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-index read-name read-validation-stringency reference "
            "sample sampling-min-mq seconds-between-progress-updates sequence-dictionary shard-size "
            "sites-only-vcf-output str-table-path threads tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "CheckReferenceCompatibility": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality output "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference references-to-compare sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater variant "
            "verbosity"
        ).split()
    ),
    "CleanSam": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "ClipReads": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes clip-adapter clip-representation clip-sequence clip-sequences-file "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 cycles-to-trim disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length min-read-length-to-output "
            "minimum-mapping-quality output output-cram-version output-statistics platform-filter-name "
            "q-trimming-threshold read read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "CollectSVEvidence": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer compression-level "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "depth-evidence-file depth-evidence-intervals depth-evidence-min-mapq depth-summary-file "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output-cram-version pe-file platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample sample-name "
            "sd-file seconds-between-progress-updates sequence-dictionary site-depth-locs-vcf "
            "site-depth-min-baseq site-depth-min-mapq sites-only-vcf-output sr-file tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "CombineGVCFs": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line "
            "allow-old-rms-mapping-quality-annotation-data ambig-filter-bases ambig-filter-frac annotation "
            "annotation-group annotations-to-exclude assembly-complexity-reference-mode black-listed-lanes "
            "break-bands-at-multiples-of call-genotypes cloud-index-prefetch-buffer cloud-prefetch-buffer "
            "combine-variants-distance convert-to-base-pair-resolution create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 dbsnp "
            "denovo-depth-threshold denovo-parent-gq-threshold dirichlet-keep-prior-in-count "
            "dirichlet-prior-pseudo-count disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-annotations "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "drop-somatic-filtering-annotations enable-all-annotations exclude-intervals filter-too-short "
            "flow-order-for-annotations founder-id gcs-max-retries ignore-variants-starting-outside-interval "
            "input input-is-somatic interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases max-distance "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version pedigree platform-filter-name "
            "pseudo-count-weight-decay-rate read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency ref-padding reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater variant "
            "variant-output-filtering verbosity"
        ).split()
    ),
    "CompareReferences": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "base-comparison base-comparison-output black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "display-only-differing-sequences display-sequences-by-name dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality md5-calculation-mode "
            "min-fragment-length min-read-length minimum-mapping-quality output output-cram-version "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference references-to-compare sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "ComposeSTRTableFile": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 decimation "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-period max-read-length max-repeats max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "Concordance": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer confidence "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends evaluation exclude-intervals "
            "filter-analysis filter-too-short filtered-true-negatives-and-false-negatives gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output summary tmp-dir true-positives-and-false-negatives "
            "true-positives-and-false-positives truth use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "CondenseDepthEvidence": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer compression-level "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "depth-evidence disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-interval-size max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-interval-size min-read-length minimum-mapping-quality output output-cram-version "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "ConvertCountsToDepthFile": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer compression-level "
            "counts-filename create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality output "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample sample-name seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "ConvertHeaderlessHadoopBamShardToBam": frozenset(
        (
            "QUIET bam-shard bam-with-header gcs-max-retries output tmp-dir use-jdk-deflater use-jdk-inflater "
            "verbosity"
        ).split()
    ),
    "CountBasesInReference": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "CountFalsePositives": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater variant variant-output-filtering verbosity"
        ).split()
    ),
    "CountVariants": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater variant variant-output-filtering verbosity"
        ).split()
    ),
    "CreateHadoopBamSplittingIndex": frozenset(
        (
            "QUIET create-bai gcs-max-retries input output read-validation-stringency splitting-index-granularity "
            "tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "DownsampleByDuplicateSet": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short fraction-to-keep gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-per-strand-reads min-read-length min-reads minimum-mapping-quality output output-cram-version "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "DownsampleSam": frozenset(
        (
            "ACCURACY COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM METRICS_FILE OUTPUT "
            "PROBABILITY QUIET RANDOM_SEED REFERENCE_SEQUENCE STRATEGY TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "EvaluateInfoFieldConcordance": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer confidence "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends epsilon eval-info-key evaluation "
            "exclude-intervals filter-too-short gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length minimum-mapping-quality output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output summary tmp-dir truth "
            "truth-info-key use-jdk-deflater use-jdk-inflater verbosity warn-big-differences"
        ).split()
    ),
    "ExtractOriginalAlignmentRecordsByNameSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes conf create-output-bam-index create-output-bam-splitting-index "
            "create-output-variant-index disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals invert-match inverted-read-filter keep-intervals "
            "keep-read-group keep-reverse-strand-only library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality num-reducers "
            "output output-cram-version output-shard-tmp-dir platform-filter-name program-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-name-file read-validation-stringency reference "
            "sample sharded-output spark-master spark-verbosity splitting-index-granularity tmp-dir "
            "use-jdk-deflater use-jdk-inflater use-nio verbosity"
        ).split()
    ),
    "ExtractSequences": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INTERVAL_LIST LINE_LENGTH MAX_RECORDS_IN_RAM OUTPUT "
            "QUIET REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "FastaAlternateReferenceMaker": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library line-width mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output snp-mask "
            "snp-mask-priority tmp-dir use-iupac-sample use-jdk-deflater use-jdk-inflater variant verbosity"
        ).split()
    ),
    "FastaReferenceMaker": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library line-width mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "FilterFuncotations": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line allele-frequency-data-source "
            "ambig-filter-bases ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality output "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency ref-version reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater variant "
            "variant-output-filtering verbosity"
        ).split()
    ),
    "FilterSamReads": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FILTER INPUT INTERVAL_LIST MAX_RECORDS_IN_RAM OUTPUT "
            "QUIET READ_LIST_FILE REFERENCE_SEQUENCE SORT_ORDER TAG TAG_VALUE TMP_DIR USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY WRITE_READS_FILES"
        ).split()
    ),
    "FindBadGenomicKmersSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes conf create-output-bam-index create-output-bam-splitting-index "
            "create-output-variant-index disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries high-copy-fasta input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter k-size "
            "keep-intervals keep-read-group keep-reverse-strand-only kmer-max-dust-score library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality num-reducers output output-cram-version "
            "output-shard-tmp-dir platform-filter-name program-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample sharded-output spark-master "
            "spark-verbosity splitting-index-granularity tmp-dir use-jdk-deflater use-jdk-inflater use-nio "
            "verbosity"
        ).split()
    ),
    "FindMendelianViolations": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FEMALE_CHROMS INPUT MALE_CHROMS MAX_RECORDS_IN_RAM "
            "MIN_DP MIN_GQ MIN_HET_FRACTION OUTPUT PSEUDO_AUTOSOMAL_REGIONS QUIET REFERENCE_SEQUENCE SKIP_CHROMS "
            "TAB_MODE THREAD_COUNT TMP_DIR TRIOS USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VCF_DIR "
            "VERBOSITY"
        ).split()
    ),
    "FixMisencodedBaseQualityReads": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "FuncotateSegments": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line alias-to-key-mapping "
            "ambig-filter-bases ambig-filter-frac annotation-default annotation-override black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 custom-variant-classification-order "
            "data-sources-path disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-field exclude-intervals filter-too-short gcs-max-retries "
            "input interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "lookahead-cache-bp mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-num-bases-for-segment-funcotation min-read-length "
            "minimum-mapping-quality output output-cram-version output-file-format platform-filter-name "
            "prefer-mane-transcripts read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency ref-version reference sample seconds-between-progress-updates segments "
            "sequence-dictionary sites-only-vcf-output splice-site-window-size tmp-dir transcript-list "
            "transcript-selection-mode use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "Funcotator": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "annotation-default annotation-override black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 custom-variant-classification-order data-sources-path "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-field exclude-intervals "
            "filter-too-short five-prime-flank-size gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library lookahead-cache-bp "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-num-bases-for-segment-funcotation min-read-length "
            "minimum-mapping-quality output output-cram-version output-file-format platform-filter-name "
            "prefer-mane-transcripts read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reannotate-vcf ref-version reference remove-filtered-variants sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output splice-site-window-size "
            "three-prime-flank-size tmp-dir transcript-list transcript-selection-mode use-jdk-deflater "
            "use-jdk-inflater variant variant-output-filtering verbosity"
        ).split()
    ),
    "GenomicsDBImport": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "avoid-nio batch-size black-listed-lanes bypass-feature-reader cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer consolidate create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries "
            "genomicsdb-segment-size genomicsdb-shared-posixfs-optimizations genomicsdb-update-workspace-path "
            "genomicsdb-use-gcs-hdfs-connector genomicsdb-vcf-buffer-size genomicsdb-workspace-path header input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-num-intervals-to-import-in-parallel max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "merge-contigs-into-num-partitions merge-input-intervals min-fragment-length min-read-length "
            "minimum-mapping-quality output-cram-version output-interval-list-to-file "
            "overwrite-existing-genomicsdb-workspace platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reader-threads reference sample sample-name-map "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater validate-sample-name-map variant verbosity"
        ).split()
    ),
    "GenotypeGVCFs": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line allele-fraction-error "
            "allow-old-rms-mapping-quality-annotation-data ambig-filter-bases ambig-filter-frac "
            "annotate-with-num-discovered-alleles annotation annotation-group annotations-to-exclude "
            "assembly-complexity-reference-mode black-listed-lanes call-genotypes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 dbsnp denovo-depth-threshold denovo-parent-gq-threshold "
            "dirichlet-keep-prior-in-count dirichlet-prior-pseudo-count disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation disable-tool-default-annotations "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends dont-use-dragstr-priors "
            "enable-all-annotations exclude-intervals filter-too-short flow-order-for-annotations "
            "force-output-intervals founder-id gcs-max-retries genomicsdb-max-alternate-alleles "
            "genomicsdb-shared-posixfs-optimizations genomicsdb-use-bcf-codec genomicsdb-use-gcs-hdfs-connector "
            "genotype-assignment-method heterozygosity heterozygosity-stdev include-non-variant-sites "
            "indel-heterozygosity input input-is-somatic interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-combined-raw-annotations "
            "keep-intervals keep-read-group keep-reverse-strand-only keep-specific-combined-raw-annotation "
            "lenient library mate-too-distant-length max-alternate-alleles max-clipped-bases max-fragment-length "
            "max-genotype-count max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality merge-input-intervals min-fragment-length "
            "min-read-length minimum-mapping-quality num-reference-samples-if-no-call "
            "only-output-calls-starting-in-intervals output output-cram-version pedigree platform-filter-name "
            "population-callset pseudo-count-weight-decay-rate read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample sample-ploidy "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output somatic-quality-threshold "
            "standard-min-confidence-threshold-for-calling tmp-dir tumor-lod-to-emit use-jdk-deflater "
            "use-jdk-inflater use-new-qual-calculator use-posteriors-to-calculate-qual variant "
            "variant-output-filtering verbosity"
        ).split()
    ),
    "GnarlyGenotyper": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "annotate-with-num-discovered-alleles black-listed-lanes call-genotypes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 dbsnp disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends dont-use-dragstr-priors exclude-intervals filter-too-short "
            "gcs-max-retries genomicsdb-max-alternate-alleles genomicsdb-shared-posixfs-optimizations "
            "genomicsdb-use-bcf-codec genomicsdb-use-gcs-hdfs-connector genotype-assignment-method heterozygosity "
            "heterozygosity-stdev indel-heterozygosity input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-all-sites keep-intervals "
            "keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-alternate-alleles max-clipped-bases max-fragment-length max-genotype-count max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality merge-input-intervals min-fragment-length min-read-length "
            "minimum-mapping-quality num-reference-samples-if-no-call only-output-calls-starting-in-intervals "
            "output output-cram-version output-database-name platform-filter-name population-callset read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample sample-ploidy "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output "
            "standard-min-confidence-threshold-for-calling tmp-dir use-jdk-deflater use-jdk-inflater "
            "use-new-qual-calculator use-posteriors-to-calculate-qual variant variant-output-filtering verbosity"
        ).split()
    ),
    "GroupedSVCluster": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line algorithm "
            "alt-allele-summary-strategy ambig-filter-bases ambig-filter-frac black-listed-lanes "
            "breakpoint-summary-strategy cloud-index-prefetch-buffer cloud-prefetch-buffer clustering-config "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "default-no-call disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends enable-cnv exclude-intervals "
            "fast-mode filter-too-short flag-field-logic gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-records-in-ram "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality omit-members "
            "output output-cram-version platform-filter-name ploidy-table read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output stratify-config stratify-num-breakpoint-overlaps "
            "stratify-num-breakpoint-overlaps-interchromosomal stratify-overlap-fraction tmp-dir track-intervals "
            "track-name use-jdk-deflater use-jdk-inflater variant variant-output-filtering variant-prefix "
            "verbosity"
        ).split()
    ),
    "GtfToBed": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries gtf-path input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output sort-by-transcript "
            "tmp-dir use-basic-transcript use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "HaplotypeBasedVariantRecaller": frozenset(
        (
            "QUIET activeregion-alt-multiplier adaptive-pruning adaptive-pruning-initial-error-rate "
            "add-output-sam-program-record add-output-vcf-command-line all-site-pls "
            "allele-informative-reads-overlap-margin alleles alleles-file-vcf allow-non-unique-kmers-in-ref "
            "ambig-filter-bases ambig-filter-frac annotate-with-num-discovered-alleles apply-bqd apply-frd "
            "bam-output bam-writer-type base-quality-score-threshold black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer comparison contamination-fraction-per-sample-file "
            "contamination-fraction-to-filter create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 dbsnp debug-assembly disable-bam-index-caching "
            "disable-cap-base-qualities-to-map-quality disable-optimizations disable-read-filter "
            "disable-sequence-dictionary-validation disable-spanning-event-genotyping "
            "disable-symmetric-hmm-normalizing disable-tool-default-read-filters "
            "do-not-correct-overlapping-quality do-not-run-physical-phasing dont-increase-kmer-sizes-for-cycles "
            "dont-require-soft-clips-both-ends dont-use-dragstr-pair-hmm-scores dont-use-dragstr-priors "
            "dont-use-soft-clipped-bases dragen-378-concordance-mode dragen-mode "
            "dragen-pd-pair-hmm-implementation dragstr-het-hom-ratio dragstr-params-path emit-ref-confidence "
            "enable-dynamic-read-disqualification-for-genotyping exclude-intervals "
            "expected-mismatch-rate-for-read-disqualification filter-too-short floor-blocks "
            "flow-assembly-collapse-partial-mode flow-disallow-probs-larger-than-call flow-fill-empty-bins-value "
            "flow-filter-alleles flow-filter-alleles-qual-threshold flow-filter-alleles-sor-threshold "
            "flow-filter-lone-alleles flow-lump-probs flow-matrix-mods flow-mode flow-probability-scaling-factor "
            "flow-quantization-bins flow-remove-non-single-base-pair-indels flow-remove-one-zero-probs "
            "flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format flow-symmetric-indel-probs "
            "flow-use-t0-tag force-call-filtered-alleles gcs-max-retries genotype-assignment-method graph-output "
            "gvcf-gq-bands haplotypes-file-bam heterozygosity heterozygosity-stdev indel-heterozygosity "
            "indel-size-to-eliminate-in-ref-model input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-boundary-flows keep-intervals "
            "keep-read-group keep-reverse-strand-only kmer-size lenient library likelihood-calculation-engine "
            "linked-de-bruijn-graph mapping-quality-threshold-for-genotyping mate-too-distant-length "
            "matrix-file-csv max-alternate-alleles max-clipped-bases max-effective-depth-adjustment-for-frd "
            "max-fragment-length max-genotype-count max-mnp-distance max-num-haplotypes-in-population "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-unpruned-variants "
            "max-variants-per-shard maximum-mapping-quality min-base-quality-score min-dangling-branch-length "
            "min-fragment-length min-pruning min-read-length minimum-mapping-quality "
            "native-dragen-pd-pair-avx-level native-dragen-pd-pair-hmm-threads native-dragen-pd-pair-memory-limit "
            "native-dragen-pd-pair-open-mp-settings native-pair-hmm-threads native-pair-hmm-use-double-precision "
            "num-pruning-samples num-reference-samples-if-no-call output-cram-version output-mode "
            "pair-hmm-gap-continuation-penalty pair-hmm-implementation pair-hmm-results-file pcr-indel-model "
            "phred-scaled-global-read-mismapping-rate pileup-detection platform-filter-name ploidy-regions "
            "population-callset pruning-lod-threshold pruning-seeding-lod-threshold read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency recover-all-dangling-branches "
            "recover-dangling-heads reference reference-model-deletion-quality sample sample-name sample-ploidy "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output smith-waterman "
            "smith-waterman-dangling-end-gap-extend-penalty smith-waterman-dangling-end-gap-open-penalty "
            "smith-waterman-dangling-end-match-value smith-waterman-dangling-end-mismatch-penalty "
            "smith-waterman-haplotype-to-reference-gap-extend-penalty "
            "smith-waterman-haplotype-to-reference-gap-open-penalty "
            "smith-waterman-haplotype-to-reference-match-value "
            "smith-waterman-haplotype-to-reference-mismatch-penalty "
            "smith-waterman-read-to-haplotype-gap-extend-penalty "
            "smith-waterman-read-to-haplotype-gap-open-penalty smith-waterman-read-to-haplotype-match-value "
            "smith-waterman-read-to-haplotype-mismatch-penalty soft-clip-low-quality-ends "
            "standard-min-confidence-threshold-for-calling tmp-dir transform-dragen-mapping-quality "
            "use-filtered-reads-for-annotations use-jdk-deflater use-jdk-inflater use-new-qual-calculator "
            "use-pdhmm use-pdhmm-overlap-optimization use-posteriors-to-calculate-qual verbosity"
        ).split()
    ),
    "IndexFeatureFile": frozenset(
        (
            "QUIET gcs-max-retries input output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "JointGermlineCNVSegmentation": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line alt-allele-summary-strategy "
            "ambig-filter-bases ambig-filter-frac autosomal-ref-copy-number black-listed-lanes "
            "breakpoint-summary-strategy cloud-index-prefetch-buffer cloud-prefetch-buffer "
            "clustering-breakend-window clustering-interval-overlap clustering-size-similarity "
            "combine-variants-distance create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 defragmentation-padding-fraction disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries "
            "ignore-variants-starting-outside-interval input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases max-distance "
            "max-fragment-length max-read-length max-records-in-ram max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length min-sample-set-fraction-overlap minimum-mapping-quality minimum-qs-score "
            "model-call-intervals output output-cram-version pedigree platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency ref-padding reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater variant variant-output-filtering verbosity"
        ).split()
    ),
    "LearnReadOrientationModel": frozenset(
        (
            "QUIET convergence-threshold gcs-max-retries input max-depth num-em-iterations output tmp-dir "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "LeftAlignIndels": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "MarkDuplicatesWithMateCigar": frozenset(
        (
            "ADD_PG_TAG_TO_READS ASSUME_SORTED ASSUME_SORT_ORDER BLOCK_SIZE COMMENT COMPRESSION_LEVEL "
            "CREATE_INDEX CREATE_MD5_FILE DUPLICATE_SCORING_STRATEGY INPUT MAX_OPTICAL_DUPLICATE_SET_SIZE "
            "MAX_RECORDS_IN_RAM METRICS_FILE MINIMUM_DISTANCE OPTICAL_DUPLICATE_PIXEL_DISTANCE OUTPUT "
            "PROGRAM_GROUP_COMMAND_LINE PROGRAM_GROUP_NAME PROGRAM_GROUP_VERSION PROGRAM_RECORD_ID QUIET "
            "READ_NAME_REGEX REFERENCE_SEQUENCE REMOVE_DUPLICATES SKIP_PAIRS_WITH_NO_MATE_CIGAR TMP_DIR "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "MergeBamAlignment": frozenset(
        (
            "ADD_MATE_CIGAR ADD_PG_TAG_TO_READS ALIGNED_BAM ALIGNED_READS_ONLY ALIGNER_PROPER_PAIR_FLAGS "
            "ATTRIBUTES_TO_REMOVE ATTRIBUTES_TO_RETAIN ATTRIBUTES_TO_REVERSE ATTRIBUTES_TO_REVERSE_COMPLEMENT "
            "CLIP_ADAPTERS CLIP_OVERLAPPING_READS COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "EXPECTED_ORIENTATIONS HARD_CLIP_OVERLAPPING_READS INCLUDE_SECONDARY_ALIGNMENTS IS_BISULFITE_SEQUENCE "
            "JUMP_SIZE MATCHING_DICTIONARY_TAGS MAX_INSERTIONS_OR_DELETIONS MAX_RECORDS_IN_RAM "
            "MIN_UNCLIPPED_BASES OUTPUT PAIRED_RUN PRIMARY_ALIGNMENT_STRATEGY PROGRAM_GROUP_COMMAND_LINE "
            "PROGRAM_GROUP_NAME PROGRAM_GROUP_VERSION PROGRAM_RECORD_ID QUIET READ1_ALIGNED_BAM READ1_TRIM "
            "READ2_ALIGNED_BAM READ2_TRIM REFERENCE_SEQUENCE SORT_ORDER TMP_DIR UNMAPPED_BAM "
            "UNMAPPED_READ_STRATEGY UNMAP_CONTAMINANT_READS USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "MergeMutect2CallsWithMC3": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer confidence "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends evaluation exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir truth "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "MergeSamFiles": frozenset(
        (
            "ASSUME_SORTED COMMENT COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT INTERVALS "
            "MAX_RECORDS_IN_RAM MERGE_SEQUENCE_DICTIONARIES OUTPUT QUIET REFERENCE_SEQUENCE SORT_ORDER TMP_DIR "
            "USE_JDK_DEFLATER USE_JDK_INFLATER USE_THREADING VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "Mutect2": frozenset(
        (
            "QUIET active-probability-threshold adaptive-pruning-initial-error-rate add-output-sam-program-record "
            "add-output-vcf-command-line af-of-alleles-not-in-resource allele-informative-reads-overlap-margin "
            "alleles allow-non-unique-kmers-in-ref allow-old-rms-mapping-quality-annotation-data "
            "ambig-filter-bases ambig-filter-frac annotation annotation-group annotations-to-exclude "
            "assembly-complexity-reference-mode assembly-region-out assembly-region-padding bam-output "
            "bam-writer-type base-qual-correction-factor base-quality-score-threshold black-listed-lanes "
            "callable-depth cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 debug-assembly "
            "denovo-depth-threshold denovo-parent-gq-threshold dirichlet-keep-prior-in-count "
            "dirichlet-prior-pseudo-count disable-adaptive-pruning disable-bam-index-caching "
            "disable-cap-base-qualities-to-map-quality disable-read-filter disable-sequence-dictionary-validation "
            "disable-symmetric-hmm-normalizing disable-tool-default-annotations disable-tool-default-read-filters "
            "dont-increase-kmer-sizes-for-cycles dont-require-soft-clips-both-ends "
            "dont-use-dragstr-pair-hmm-scores dont-use-soft-clipped-bases downsampling-stride "
            "dragen-pd-pair-hmm-implementation dragstr-het-hom-ratio dragstr-params-path emit-ref-confidence "
            "enable-all-annotations enable-dynamic-read-disqualification-for-genotyping exclude-intervals "
            "expected-mismatch-rate-for-read-disqualification f1r2-max-depth f1r2-median-mq f1r2-min-bq "
            "f1r2-tar-gz filter-too-short flow-assembly-collapse-partial-mode "
            "flow-disallow-probs-larger-than-call flow-fill-empty-bins-value flow-filter-alleles "
            "flow-filter-alleles-qual-threshold flow-filter-alleles-sor-threshold flow-filter-lone-alleles "
            "flow-lump-probs flow-matrix-mods flow-mode flow-order-for-annotations "
            "flow-probability-scaling-factor flow-quantization-bins flow-remove-non-single-base-pair-indels "
            "flow-remove-one-zero-probs flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-symmetric-indel-probs flow-use-t0-tag force-active force-call-filtered-alleles founder-id "
            "gcs-max-retries genotype-germline-sites genotype-germline-sites-fraction genotype-pon-sites "
            "germline-resource graph-output gvcf-lod-band ignore-itr-artifacts independent-mates "
            "initial-tumor-lod input interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-boundary-flows keep-intervals keep-read-group "
            "keep-reverse-strand-only kmer-size lenient library likelihood-calculation-engine "
            "linked-de-bruijn-graph mate-too-distant-length max-assembly-region-size max-clipped-bases "
            "max-fragment-length max-mnp-distance max-num-haplotypes-in-population max-population-af "
            "max-prob-propagation-distance max-read-length max-reads-per-alignment-start "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-suspicious-reads-per-alignment-start max-unpruned-variants max-variants-per-shard "
            "maximum-mapping-quality min-assembly-region-size min-base-quality-score min-dangling-branch-length "
            "min-fragment-length min-pruning min-read-length minimum-allele-fraction minimum-mapping-quality "
            "mitochondria-mode native-dragen-pd-pair-avx-level native-dragen-pd-pair-hmm-threads "
            "native-dragen-pd-pair-memory-limit native-dragen-pd-pair-open-mp-settings native-pair-hmm-threads "
            "native-pair-hmm-use-double-precision normal-lod normal-sample num-pruning-samples output "
            "output-cram-version pair-hmm-gap-continuation-penalty pair-hmm-implementation pair-hmm-results-file "
            "panel-of-normals pcr-indel-model pcr-indel-qual pcr-snv-qual pedigree permutect-alt-downsample "
            "permutect-dataset-mode permutect-non-artifact-ratio permutect-ref-downsample permutect-test-dataset "
            "permutect-test-truth permutect-training-dataset permutect-training-truth "
            "phred-scaled-global-read-mismapping-rate pileup-detection platform-filter-name pruning-lod-threshold "
            "pruning-seeding-lod-threshold pseudo-count-weight-decay-rate read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency recover-all-dangling-branches reference "
            "reference-model-deletion-quality sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output smith-waterman smith-waterman-dangling-end-gap-extend-penalty "
            "smith-waterman-dangling-end-gap-open-penalty smith-waterman-dangling-end-match-value "
            "smith-waterman-dangling-end-mismatch-penalty "
            "smith-waterman-haplotype-to-reference-gap-extend-penalty "
            "smith-waterman-haplotype-to-reference-gap-open-penalty "
            "smith-waterman-haplotype-to-reference-match-value "
            "smith-waterman-haplotype-to-reference-mismatch-penalty "
            "smith-waterman-read-to-haplotype-gap-extend-penalty "
            "smith-waterman-read-to-haplotype-gap-open-penalty smith-waterman-read-to-haplotype-match-value "
            "smith-waterman-read-to-haplotype-mismatch-penalty soft-clip-low-quality-ends tmp-dir "
            "tumor-lod-to-emit tumor-sample use-jdk-deflater use-jdk-inflater use-pdhmm "
            "use-pdhmm-overlap-optimization verbosity"
        ).split()
    ),
    "NonNFastaSize": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT INTERVALS MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "NormalizeFasta": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT LINE_LENGTH MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "PathSeqBwaSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes bwa-score-threshold conf create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only library "
            "mate-too-distant-length max-alternate-hits max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio maximum-mapping-quality "
            "microbe-bwa-image microbe-dict microbe-min-seed-length min-fragment-length min-read-length "
            "minimum-mapping-quality num-reducers output-cram-version output-shard-tmp-dir paired-input "
            "paired-output platform-filter-name program-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample sharded-output spark-master "
            "spark-verbosity splitting-index-granularity tmp-dir unpaired-input unpaired-output use-jdk-deflater "
            "use-jdk-inflater use-nio verbosity"
        ).split()
    ),
    "PathSeqFilterSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes conf create-output-bam-index create-output-bam-splitting-index "
            "create-output-variant-index disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends dust-mask-quality dust-t "
            "dust-window exclude-intervals filter-bwa-image filter-bwa-seed-length filter-duplicates "
            "filter-metrics filter-reads-per-partition filter-too-short gcs-max-retries host-kmer-thresh "
            "host-min-identity ignore-alignment-contigs input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter is-host-aligned keep-intervals "
            "keep-read-group keep-reverse-strand-only kmer-file library mate-too-distant-length "
            "max-adapter-mismatches max-clipped-bases max-fragment-length max-masked-bases max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio maximum-mapping-quality "
            "min-adapter-length min-base-quality min-clipped-read-length min-fragment-length min-read-length "
            "minimum-mapping-quality num-reducers output-cram-version output-shard-tmp-dir paired-output "
            "platform-filter-name program-name quality-threshold read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample sharded-output "
            "skip-pre-bwa-repartition skip-quality-filters spark-master spark-verbosity "
            "splitting-index-granularity tmp-dir unpaired-output use-jdk-deflater use-jdk-inflater use-nio "
            "verbosity"
        ).split()
    ),
    "PathSeqPipelineSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes bwa-score-threshold conf create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters divide-by-genome-length "
            "dont-require-soft-clips-both-ends dust-mask-quality dust-t dust-window exclude-intervals "
            "filter-bwa-image filter-bwa-seed-length filter-duplicates filter-metrics filter-reads-per-partition "
            "filter-too-short gcs-max-retries host-kmer-thresh host-min-identity identity-margin "
            "ignore-alignment-contigs input interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter is-host-aligned keep-intervals keep-read-group "
            "keep-reverse-strand-only kmer-file library mate-too-distant-length max-adapter-mismatches "
            "max-alternate-hits max-clipped-bases max-fragment-length max-masked-bases max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio maximum-mapping-quality "
            "microbe-bwa-image microbe-dict microbe-min-seed-length min-adapter-length min-base-quality "
            "min-clipped-read-length min-fragment-length min-read-length min-score-identity "
            "minimum-mapping-quality not-normalized-by-kingdom num-reducers output output-cram-version "
            "output-shard-tmp-dir pipeline-reads-per-partition platform-filter-name program-name "
            "quality-threshold read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency readsPerPartitionOutput reference sample score-metrics "
            "score-reads-per-partition-estimate score-warnings scores-output sharded-output "
            "skip-pre-bwa-repartition skip-quality-filters spark-master spark-verbosity "
            "splitting-index-granularity taxonomy-file tmp-dir use-jdk-deflater use-jdk-inflater use-nio "
            "verbosity"
        ).split()
    ),
    "PathSeqScoreSpark": frozenset(
        (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac bam-partition-size "
            "black-listed-lanes conf create-output-bam-index create-output-bam-splitting-index "
            "create-output-variant-index disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters divide-by-genome-length dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries identity-margin input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio maximum-mapping-quality min-fragment-length min-read-length "
            "min-score-identity minimum-mapping-quality not-normalized-by-kingdom num-reducers output "
            "output-cram-version output-shard-tmp-dir paired-input platform-filter-name program-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample score-metrics "
            "score-reads-per-partition-estimate score-warnings scores-output sharded-output spark-master "
            "spark-verbosity splitting-index-granularity taxonomy-file tmp-dir unpaired-input use-jdk-deflater "
            "use-jdk-inflater use-nio verbosity"
        ).split()
    ),
    "PositionBasedDownsampleSam": frozenset(
        (
            "ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FRACTION "
            "INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET READ_NAME_REGEX REFERENCE_SEQUENCE "
            "REMOVE_DUPLICATE_INFORMATION STOP_AFTER TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "PrintBGZFBlockInformation": frozenset(
        (
            "QUIET bgzf-file gcs-max-retries output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "PrintDistantMates": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "PrintFileDiagnostics": frozenset(
        (
            "QUIET count-limit gcs-max-retries input output tmp-dir use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "PrintReadCounts": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input input-counts interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output-cram-version output-file-list output-prefix platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-index read-name read-validation-stringency reference "
            "sample seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "PrintSVEvidence": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer compression-level "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends evidence-file exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample sample-names "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "ReadAnonymizer": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency ref-base-quality reference "
            "sample seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir "
            "use-jdk-deflater use-jdk-inflater use-simple-cigar verbosity"
        ).split()
    ),
    "ReadsPipelineSpark": frozenset(
        (
            "QUIET active-probability-threshold activeregion-alt-multiplier adaptive-pruning "
            "adaptive-pruning-initial-error-rate add-output-vcf-command-line align all-site-pls "
            "allele-informative-reads-overlap-margin alleles allow-missing-read-group "
            "allow-non-unique-kmers-in-ref allow-old-rms-mapping-quality-annotation-data ambig-filter-bases "
            "ambig-filter-frac annotate-with-num-discovered-alleles annotation annotation-group "
            "annotations-to-exclude apply-bqd apply-frd assembly-complexity-reference-mode "
            "assembly-region-padding bam-output bam-partition-size bam-writer-type base-quality-score-threshold "
            "binary-tag-name black-listed-lanes bqsr-baq-gap-open-penalty bwa-mem-index-image comparison conf "
            "contamination-fraction-per-sample-file contamination-fraction-to-filter create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index dbsnp debug-assembly "
            "default-base-qualities deletions-default-quality denovo-depth-threshold denovo-parent-gq-threshold "
            "dirichlet-keep-prior-in-count dirichlet-prior-pseudo-count disable-cap-base-qualities-to-map-quality "
            "disable-optimizations disable-read-filter disable-sequence-dictionary-validation "
            "disable-spanning-event-genotyping disable-symmetric-hmm-normalizing disable-tool-default-annotations "
            "disable-tool-default-read-filters do-not-correct-overlapping-quality do-not-mark-unmapped-mates "
            "do-not-run-physical-phasing dont-increase-kmer-sizes-for-cycles dont-require-soft-clips-both-ends "
            "dont-use-dragstr-pair-hmm-scores dont-use-dragstr-priors dont-use-soft-clipped-bases "
            "dragen-378-concordance-mode dragen-mode dragen-pd-pair-hmm-implementation dragstr-het-hom-ratio "
            "dragstr-params-path duplicate-scoring-strategy duplicate-tagging-policy emit-original-quals "
            "emit-ref-confidence enable-all-annotations enable-dynamic-read-disqualification-for-genotyping "
            "exclude-intervals expected-mismatch-rate-for-read-disqualification filter-too-short floor-blocks "
            "flow-assembly-collapse-partial-mode flow-disallow-probs-larger-than-call flow-end-pos-uncertainty "
            "flow-fill-empty-bins-value flow-filter-alleles flow-filter-alleles-qual-threshold "
            "flow-filter-alleles-sor-threshold flow-filter-lone-alleles flow-lump-probs flow-matrix-mods "
            "flow-mode flow-order-for-annotations flow-probability-scaling-factor flow-q-is-known-end "
            "flow-quality-sum-strategy flow-quantization-bins flow-remove-non-single-base-pair-indels "
            "flow-remove-one-zero-probs flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-skip-start-homopolymers flow-symmetric-indel-probs flow-use-t0-tag flowbased force-active "
            "force-call-filtered-alleles founder-id gcs-max-retries genotype-assignment-method "
            "global-qscore-prior graph-output gvcf-gq-bands heterozygosity heterozygosity-stdev "
            "indel-heterozygosity indel-size-to-eliminate-in-ref-model indels-context-size input "
            "insertions-default-quality interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-boundary-flows keep-intervals keep-read-group "
            "keep-reverse-strand-only kmer-size known-sites library likelihood-calculation-engine "
            "linked-de-bruijn-graph low-quality-tail mapping-quality-threshold-for-genotyping "
            "mate-too-distant-length max-alternate-alleles max-assembly-region-size max-clipped-bases "
            "max-effective-depth-adjustment-for-frd max-fragment-length max-genotype-count max-mnp-distance "
            "max-num-haplotypes-in-population max-prob-propagation-distance max-read-length "
            "max-reads-per-alignment-start max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-unpruned-variants maximum-cycle-value maximum-mapping-quality min-assembly-region-size "
            "min-base-quality-score min-dangling-branch-length min-fragment-length min-pruning min-read-length "
            "minimum-mapping-quality mismatches-context-size mismatches-default-quality "
            "native-dragen-pd-pair-avx-level native-dragen-pd-pair-hmm-threads native-dragen-pd-pair-memory-limit "
            "native-dragen-pd-pair-open-mp-settings native-pair-hmm-threads native-pair-hmm-use-double-precision "
            "num-pruning-samples num-reducers num-reference-samples-if-no-call output output-bam "
            "output-cram-version output-mode output-shard-tmp-dir pair-hmm-gap-continuation-penalty "
            "pair-hmm-implementation pair-hmm-results-file pcr-indel-model pedigree "
            "phred-scaled-global-read-mismapping-rate pileup-detection platform-filter-name ploidy-regions "
            "population-callset preserve-qscores-less-than program-name pruning-lod-threshold "
            "pruning-seeding-lod-threshold pseudo-count-weight-decay-rate quantize-quals quantizing-levels "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-index read-name read-shard-padding read-shard-size "
            "read-validation-stringency recover-all-dangling-branches recover-dangling-heads reference "
            "reference-model-deletion-quality remove-all-duplicates remove-sequencing-duplicates "
            "round-down-quantized sample sample-name sample-ploidy sharded-output single-end-alignment "
            "single-end-reads-clipping-is-end single-end-reads-end-position-significant smith-waterman "
            "smith-waterman-dangling-end-gap-extend-penalty smith-waterman-dangling-end-gap-open-penalty "
            "smith-waterman-dangling-end-match-value smith-waterman-dangling-end-mismatch-penalty "
            "smith-waterman-haplotype-to-reference-gap-extend-penalty "
            "smith-waterman-haplotype-to-reference-gap-open-penalty "
            "smith-waterman-haplotype-to-reference-match-value "
            "smith-waterman-haplotype-to-reference-mismatch-penalty "
            "smith-waterman-read-to-haplotype-gap-extend-penalty "
            "smith-waterman-read-to-haplotype-gap-open-penalty smith-waterman-read-to-haplotype-match-value "
            "smith-waterman-read-to-haplotype-mismatch-penalty soft-clip-low-quality-ends spark-master "
            "spark-verbosity splitting-index-granularity standard-min-confidence-threshold-for-calling "
            "static-quantized-quals strict tmp-dir transform-dragen-mapping-quality "
            "use-filtered-reads-for-annotations use-jdk-deflater use-jdk-inflater use-new-qual-calculator use-nio "
            "use-original-qualities use-pdhmm use-pdhmm-overlap-optimization use-posteriors-to-calculate-qual "
            "verbosity"
        ).split()
    ),
    "ReblockGVCF": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line add-site-filters-to-genotype "
            "allow-missing-hom-ref-data allow-old-rms-mapping-quality-annotation-data ambig-filter-bases "
            "ambig-filter-frac annotate-with-num-discovered-alleles annotation annotation-group "
            "annotations-to-exclude annotations-to-keep assembly-complexity-reference-mode black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 dbsnp denovo-depth-threshold "
            "denovo-parent-gq-threshold dirichlet-keep-prior-in-count dirichlet-prior-pseudo-count "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-annotations disable-tool-default-read-filters do-qual-score-approximation "
            "dont-require-soft-clips-both-ends dont-use-dragstr-priors drop-low-quals enable-all-annotations "
            "exclude-intervals filter-too-short floor-blocks flow-order-for-annotations "
            "format-annotations-to-remove founder-id gcs-max-retries genotype-assignment-method gvcf-gq-bands "
            "heterozygosity heterozygosity-stdev indel-heterozygosity input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-all-alts keep-intervals keep-read-group keep-reverse-strand-only keep-site-filters lenient "
            "library mate-too-distant-length max-alternate-alleles max-clipped-bases max-fragment-length "
            "max-genotype-count max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality num-reference-samples-if-no-call output output-cram-version pedigree "
            "platform-filter-name population-callset pseudo-count-weight-decay-rate read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference "
            "rgq-threshold-to-no-call sample sample-ploidy seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output standard-min-confidence-threshold-for-calling tmp-dir "
            "tree-score-threshold-to-no-call use-jdk-deflater use-jdk-inflater use-new-qual-calculator "
            "use-posteriors-to-calculate-qual variant variant-output-filtering verbosity"
        ).split()
    ),
    "ReferenceBlockConcordance": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer confidence "
            "confidence-concordance-histogram create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends eval-block-histogram evaluation exclude-intervals filter-too-short "
            "gcs-max-retries input interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output-cram-version platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir truth truth-block-histogram use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "ReorderSam": frozenset(
        (
            "ALLOW_CONTIG_LENGTH_DISCORDANCE ALLOW_INCOMPLETE_DICT_CONCORDANCE COMPRESSION_LEVEL CREATE_INDEX "
            "CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET REFERENCE_SEQUENCE SEQUENCE_DICTIONARY TMP_DIR "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "ReplaceSamHeader": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE HEADER INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "RevertBaseQualityScores": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "RevertOriginalBaseQualitiesAndAddMateCigar": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM MAX_RECORDS_TO_EXAMINE "
            "OUTPUT QUIET REFERENCE_SEQUENCE RESTORE_ORIGINAL_QUALITIES SORT_ORDER TMP_DIR USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "RevertSam": frozenset(
        (
            "ATTRIBUTE_TO_CLEAR ATTRIBUTE_TO_REVERSE ATTRIBUTE_TO_REVERSE_COMPLEMENT COMPRESSION_LEVEL "
            "CREATE_INDEX CREATE_MD5_FILE INPUT KEEP_FIRST_DUPLICATE LIBRARY_NAME MAX_DISCARD_FRACTION "
            "MAX_RECORDS_IN_RAM OUTPUT OUTPUT_BY_READGROUP OUTPUT_BY_READGROUP_FILE_FORMAT OUTPUT_MAP QUIET "
            "REFERENCE_SEQUENCE REMOVE_ALIGNMENT_INFORMATION REMOVE_DUPLICATE_INFORMATION RESTORE_HARDCLIPS "
            "RESTORE_ORIGINAL_QUALITIES SAMPLE_ALIAS SANITIZE SORT_ORDER TMP_DIR USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "SVAnnotate": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-breakend-as-cnv-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length minimum-mapping-quality non-coding-bed output output-cram-version "
            "platform-filter-name promoter-window-length protein-coding-gtf read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater variant "
            "variant-output-filtering verbosity"
        ).split()
    ),
    "SVCluster": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line algorithm "
            "alt-allele-summary-strategy ambig-filter-bases ambig-filter-frac black-listed-lanes "
            "breakpoint-summary-strategy cloud-index-prefetch-buffer cloud-prefetch-buffer "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "default-no-call defrag-padding-fraction defrag-sample-overlap depth-breakend-window "
            "depth-interval-overlap depth-sample-overlap depth-size-similarity disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends enable-cnv exclude-intervals fast-mode filter-too-short "
            "flag-field-logic gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-records-in-ram max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length minimum-mapping-quality mixed-breakend-window mixed-interval-overlap "
            "mixed-sample-overlap mixed-size-similarity omit-members output output-cram-version "
            "pesr-breakend-window pesr-interval-overlap pesr-sample-overlap pesr-size-similarity "
            "platform-filter-name ploidy-table read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater variant variant-output-filtering "
            "variant-prefix verbosity"
        ).split()
    ),
    "SVConcordance": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer clustering-config confidence "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "depth-breakend-window depth-interval-overlap depth-sample-overlap depth-size-similarity "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends evaluation exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-records-in-ram max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-fragment-length "
            "min-read-length minimum-mapping-quality mixed-breakend-window mixed-interval-overlap "
            "mixed-sample-overlap mixed-size-similarity output output-cram-version pesr-breakend-window "
            "pesr-interval-overlap pesr-sample-overlap pesr-size-similarity platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output stratify-config "
            "stratify-num-breakpoint-overlaps stratify-num-breakpoint-overlaps-interchromosomal "
            "stratify-overlap-fraction tmp-dir track-intervals track-name truth use-jdk-deflater use-jdk-inflater "
            "verbosity"
        ).split()
    ),
    "SVStratify": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line allow-multiple-matches "
            "ambig-filter-bases ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality output "
            "output-cram-version output-prefix platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample seconds-between-progress-updates "
            "sequence-dictionary sites-only-vcf-output split-output stratify-config "
            "stratify-num-breakpoint-overlaps stratify-num-breakpoint-overlaps-interchromosomal "
            "stratify-overlap-fraction tmp-dir track-intervals track-name use-jdk-deflater use-jdk-inflater "
            "variant variant-output-filtering verbosity"
        ).split()
    ),
    "SamFormatConverter": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "SamToFastqWithTags": frozenset(
        (
            "CLIPPING_ACTION CLIPPING_ATTRIBUTE CLIPPING_MIN_LENGTH COMPRESSION_LEVEL COMPRESS_OUTPUTS_PER_RG "
            "COMPRESS_OUTPUTS_PER_TAG_GROUP CREATE_INDEX CREATE_MD5_FILE FASTQ INCLUDE_NON_PF_READS "
            "INCLUDE_NON_PRIMARY_ALIGNMENTS INPUT INTERLEAVE MAX_RECORDS_IN_RAM OUTPUT_DIR OUTPUT_PER_RG QUALITY "
            "QUALITY_TAG_GROUP QUIET READ1_MAX_BASES_TO_WRITE READ1_TRIM READ2_MAX_BASES_TO_WRITE READ2_TRIM "
            "REFERENCE_SEQUENCE RE_REVERSE RG_TAG SECOND_END_FASTQ SEQUENCE_TAG_GROUP TAG_GROUP_SEPERATOR TMP_DIR "
            "UNPAIRED_FASTQ USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "ScatterIntervalsByNs": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM MAX_TO_MERGE OUTPUT OUTPUT_TYPE "
            "QUIET REFERENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "ShiftFasta": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-file-name "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library line-width "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality output "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "shift-back-output shift-offset-list sites-only-vcf-output tmp-dir use-jdk-deflater use-jdk-inflater "
            "verbosity"
        ).split()
    ),
    "SimpleMarkDuplicatesWithMateCigar": frozenset(
        (
            "ADD_PG_TAG_TO_READS ASSUME_SORTED ASSUME_SORT_ORDER BARCODE_TAG CLEAR_DT COMMENT COMPRESSION_LEVEL "
            "CREATE_INDEX CREATE_MD5_FILE DUPLEX_UMI DUPLICATE_SCORING_STRATEGY FLOW_DUP_STRATEGY "
            "FLOW_EFFECTIVE_QUALITY_THRESHOLD FLOW_MODE FLOW_Q_IS_KNOWN_END FLOW_SKIP_FIRST_N_FLOWS "
            "FLOW_UNPAIRED_END_UNCERTAINTY FLOW_UNPAIRED_START_UNCERTAINTY FLOW_USE_END_IN_UNPAIRED_READS "
            "FLOW_USE_UNPAIRED_CLIPPED_END INPUT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP "
            "MAX_OPTICAL_DUPLICATE_SET_SIZE MAX_RECORDS_IN_RAM MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP METRICS_FILE "
            "MOLECULAR_IDENTIFIER_TAG OPTICAL_DUPLICATE_PIXEL_DISTANCE OUTPUT PROGRAM_GROUP_COMMAND_LINE "
            "PROGRAM_GROUP_NAME PROGRAM_GROUP_VERSION PROGRAM_RECORD_ID QUIET READ_NAME_REGEX "
            "READ_ONE_BARCODE_TAG READ_TWO_BARCODE_TAG REFERENCE_SEQUENCE REMOVE_DUPLICATES "
            "REMOVE_SEQUENCING_DUPLICATES SORTING_COLLECTION_SIZE_RATIO TAGGING_POLICY TAG_DUPLICATE_SET_MEMBERS "
            "TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "SiteDepthtoBAF": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "baf-evidence-output baf-sites-vcf black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer compression-level create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries input "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio max-std max-variants-per-shard "
            "maximum-mapping-quality min-fragment-length min-het-probability min-read-length min-total-depth "
            "minimum-mapping-quality output-cram-version platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-index read-name read-validation-stringency reference sample sample-names "
            "seconds-between-progress-updates sequence-dictionary site-depth sites-only-vcf-output tmp-dir "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "SortGff": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE SEQUENCE_DICTIONARY TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY nRecordsInMemory"
        ).split()
    ),
    "SortSam": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE SORT_ORDER TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY "
            "VERBOSITY"
        ).split()
    ),
    "SplitNCigarReads": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters do-not-fix-overhangs dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-bases-in-overhang max-clipped-bases max-fragment-length max-mismatches-in-overhang "
            "max-read-length max-reads-in-memory max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name process-secondary-alignments "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-index read-name read-validation-stringency "
            "refactor-cigar-string reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output skip-mapping-quality-transform tmp-dir use-jdk-deflater use-jdk-inflater "
            "verbosity"
        ).split()
    ),
    "SplitSamByLibrary": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT QUIET "
            "REFERENCE_SEQUENCE TMP_DIR USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "SplitSamByNumberOfReads": frozenset(
        (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INPUT MAX_RECORDS_IN_RAM OUTPUT OUT_PREFIX QUIET "
            "REFERENCE_SEQUENCE SPLIT_TO_N_FILES SPLIT_TO_N_READS TMP_DIR TOTAL_READS_IN_INPUT USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "StructuralVariationDiscoveryPipelineSpark": frozenset(
        (
            "QUIET adapter-sequence add-output-vcf-command-line aligner-index-image "
            "allowed-short-fragment-overhang ambig-filter-bases ambig-filter-frac assembled-contigs-output-order "
            "assembly-imprecise-evidence-overlap-uncertainty assembly-to-mapped-size-ratio-guess "
            "bam-partition-size black-listed-lanes breakpoint-evidence-dir breakpoint-intervals "
            "cleaner-max-copy-number cleaner-max-intervals cleaner-min-kmer-count cnv-calls conf contig-sam-file "
            "create-output-bam-index create-output-bam-splitting-index create-output-variant-index "
            "cross-contigs-to-ignore debug-mode disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "exclusion-interval-padding exclusion-intervals exp-interpret expand-assembly-graph external-evidence "
            "external-evidence-uncertainty external-evidence-weight fastq-dir filter-too-short gcs-max-retries "
            "high-coverage-intervals high-depth-coverage-factor high-depth-coverage-peak-factor "
            "imprecise-variant-evidence-threshold include-mapping-location input interval-exclusion-padding "
            "interval-merging-rule interval-only-assembly interval-padding interval-set-rule intervals "
            "inverted-read-filter k-size keep-intervals keep-read-group keep-reverse-strand-only kmer-intervals "
            "kmer-max-dust-score kmers-to-ignore library mate-too-distant-length "
            "max-callable-imprecise-deletion-size max-clipped-bases max-fastq-size max-fragment-length "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-tracked-fragment-length maximum-mapping-quality min-align-length "
            "min-coherent-evidence-coverage-ratio min-evidence-coverage-ratio min-evidence-mapq "
            "min-evidence-match-length min-fragment-length min-kmers-per-interval min-mq min-read-length "
            "minimum-mapping-quality num-reducers output output-cram-version output-shard-tmp-dir "
            "platform-filter-name pop-variant-bubbles program-name qname-intervals-for-assembly "
            "qname-intervals-mapped read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-metadata read-name "
            "read-validation-stringency reference remove-shadowed-contigs run-without-gaps-annotation "
            "run-without-umap-s100-annotation sample sharded-output spark-master spark-verbosity "
            "splitting-index-granularity sv-evidence-filter-model-file sv-evidence-filter-threshold-probability "
            "sv-genome-gaps-file sv-genome-umap-s100-file target-link-file tmp-dir truth-interval-padding "
            "unfiltered-breakpoint-evidence-dir use-jdk-deflater use-jdk-inflater use-nio verbosity write-gfas "
            "z-dropoff"
        ).split()
    ),
    "TransferReadTags": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-tags read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir unmapped-sam "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ).split()
    ),
    "UmiAwareMarkDuplicatesWithMateCigar": frozenset(
        (
            "ADD_PG_TAG_TO_READS ALLOW_MISSING_UMIS ASSUME_SORTED ASSUME_SORT_ORDER BARCODE_TAG CLEAR_DT COMMENT "
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE DUPLEX_UMI DUPLICATE_SCORING_STRATEGY "
            "FLOW_DUP_STRATEGY FLOW_EFFECTIVE_QUALITY_THRESHOLD FLOW_MODE FLOW_Q_IS_KNOWN_END "
            "FLOW_SKIP_FIRST_N_FLOWS FLOW_UNPAIRED_END_UNCERTAINTY FLOW_UNPAIRED_START_UNCERTAINTY "
            "FLOW_USE_END_IN_UNPAIRED_READS FLOW_USE_UNPAIRED_CLIPPED_END INPUT MAX_EDIT_DISTANCE_TO_JOIN "
            "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP MAX_OPTICAL_DUPLICATE_SET_SIZE MAX_RECORDS_IN_RAM "
            "MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP METRICS_FILE MOLECULAR_IDENTIFIER_TAG "
            "OPTICAL_DUPLICATE_PIXEL_DISTANCE OUTPUT PROGRAM_GROUP_COMMAND_LINE PROGRAM_GROUP_NAME "
            "PROGRAM_GROUP_VERSION PROGRAM_RECORD_ID QUIET READ_NAME_REGEX READ_ONE_BARCODE_TAG "
            "READ_TWO_BARCODE_TAG REFERENCE_SEQUENCE REMOVE_DUPLICATES REMOVE_SEQUENCING_DUPLICATES "
            "SORTING_COLLECTION_SIZE_RATIO TAGGING_POLICY TAG_DUPLICATE_SET_MEMBERS TMP_DIR UMI_METRICS_FILE "
            "UMI_TAG_NAME USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ).split()
    ),
    "UnmarkDuplicates": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends exclude-intervals "
            "filter-too-short gcs-max-retries input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output output-cram-version platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-index read-name read-validation-stringency reference sample "
            "seconds-between-progress-updates sequence-dictionary sites-only-vcf-output tmp-dir use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ).split()
    ),
    "ValidateBasicSomaticShortMutations": frozenset(
        (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "annotated-vcf black-listed-lanes cloud-index-prefetch-buffer cloud-prefetch-buffer "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters discovery-sample-name dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries input interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-intervals keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-validation-normal-count max-variants-per-shard maximum-mapping-quality "
            "min-base-quality-cutoff min-fragment-length min-power min-read-length minimum-mapping-quality output "
            "output-cram-version platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list read-index read-name "
            "read-validation-stringency reference sample seconds-between-progress-updates sequence-dictionary "
            "sites-only-vcf-output summary tmp-dir use-jdk-deflater use-jdk-inflater val-case-sample-name "
            "val-control-sample-name variant variant-output-filtering verbosity"
        ).split()
    ),
}
