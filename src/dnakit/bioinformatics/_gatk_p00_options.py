"""Options that ``gatk <Tool> --help`` lists for the GATK tools of pack P00 (GATK 4.6.2.0).

Generated from the help output of each tool. Per tool four groups of long option names (space
separated): ``value`` is passed as given, ``path`` is resolved to an absolute path, ``out`` must
be a new file below an existing folder, ``dir`` is a folder path. ``--help``, ``--version``,
``--showHidden``, ``--arguments_file`` and ``--gatk-config-file`` are not offered.
"""

from __future__ import annotations

_RAW: dict[str, dict[str, str]] = {
    "ASEReadCounter": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer count-overlap-reads-handling create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-base-quality "
            "min-depth-of-non-filtered-base min-fragment-length min-mapping-quality "
            "min-read-length minimum-mapping-quality output-format platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary variant"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "AddFlowBaseQuality": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short flow-disallow-probs-larger-than-call "
            "flow-fill-empty-bins-value flow-lump-probs flow-matrix-mods "
            "flow-probability-scaling-factor flow-quantization-bins "
            "flow-remove-non-single-base-pair-indels flow-remove-one-zero-probs "
            "flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-symmetric-indel-probs flow-use-t0-tag gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-boundary-flows keep-intervals keep-read-group keep-reverse-strand-only lenient "
            "library mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximal-quality-score maximum-mapping-quality "
            "min-fragment-length min-read-length minimal-error-rate minimum-mapping-quality "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency replace-quality-mode sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "AddFlowSNVQuality": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short flow-disallow-probs-larger-than-call "
            "flow-fill-empty-bins-value flow-lump-probs flow-matrix-mods "
            "flow-probability-scaling-factor flow-quantization-bins "
            "flow-remove-non-single-base-pair-indels flow-remove-one-zero-probs "
            "flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-symmetric-indel-probs flow-use-t0-tag gcs-max-retries "
            "gcs-project-for-requester-pays include-qc-failed-read interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-boundary-flows keep-intervals keep-read-group "
            "keep-reverse-strand-only keep-supplementary-alignments lenient library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-phred-score "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output-quality-attribute platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output snvq-mode use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "AnalyzeCovariates": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays ignore-last-modification-times "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("after-report-file before-report-file bqsr-recal-file"),
        "out": ("intermediate-csv-file plots-report-file"),
        "dir": ("tmp-dir"),
    },
    "AnalyzeSaturationMutagenesis": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer codon-translation create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-ignore-disjoint-pairs "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short "
            "find-large-deletions gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality min-alt-length "
            "min-flanking-length min-fragment-length min-length min-mapq min-q min-read-length "
            "min-variant-obs minimum-mapping-quality orf paired-mode platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity "
            "write-rejected-reads"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output-file-prefix"),
        "dir": ("tmp-dir"),
    },
    "AnnotateIntervals": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals feature-query-lookahead filter-too-short gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals "
            "keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": (
            "input mappability-track reference segmental-duplication-track sequence-dictionary"
        ),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "BedToIntervalList": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE KEEP_LENGTH_ZERO_INTERVALS "
            "MAX_RECORDS_IN_RAM QUIET SORT UNIQUE USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "BpmToNormalizationManifestCsv": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("CLUSTER_FILE INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CalculateAverageCombinedAnnotations": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output summed-annotation-to-divide use-jdk-deflater use-jdk-inflater "
            "variant-output-filtering verbosity"
        ),
        "path": ("input reference sequence-dictionary variant"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CalculateFingerprintMetrics": {
        "value": (
            "CALCULATE_BY COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE GENOTYPE_LOD_THRESHOLD "
            "INPUT MAX_RECORDS_IN_RAM NUMBER_OF_SAMPLING QUIET USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("HAPLOTYPE_MAP REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CalculateReadGroupChecksum": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CallCopyRatioSegments": {
        "value": (
            "QUIET calling-copy-ratio-z-score-threshold gcs-max-retries "
            "gcs-project-for-requester-pays neutral-segment-copy-ratio-lower-bound "
            "neutral-segment-copy-ratio-upper-bound "
            "outlier-neutral-segment-copy-ratio-z-score-threshold use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CallableLoci": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short format gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals "
            "keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-depth max-depth-per-sample max-fraction-of-reads-with-low-mapq "
            "max-fragment-length max-low-mapq max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-base-quality min-depth "
            "min-depth-for-low-mapq min-fragment-length min-mapping-quality min-read-length "
            "minimum-mapping-quality platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output summary"),
        "dir": ("tmp-dir"),
    },
    "CheckDuplicateMarking": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM MODE QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CheckFingerprint": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND EXIT_CODE_WHEN_NO_VALID_CHECKS "
            "EXPECTED_SAMPLE_ALIAS GENOTYPES GENOTYPE_LOD_THRESHOLD IGNORE_READ_GROUPS INPUT "
            "MAX_RECORDS_IN_RAM OBSERVED_SAMPLE_ALIAS QUIET USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("HAPLOTYPE_MAP REFERENCE_SEQUENCE"),
        "out": ("DETAIL_OUTPUT OUTPUT SUMMARY_OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CheckPileup": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer continue-after-error create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "ignore-overlaps interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input pileup reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CheckTerminatorBlock": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "dir": ("TMP_DIR"),
    },
    "CollectAlignmentSummaryMetrics": {
        "value": (
            "ADAPTER_SEQUENCE ASSUME_SORTED COLLECT_ALIGNMENT_INFORMATION COMPRESSION_LEVEL "
            "CREATE_INDEX CREATE_MD5_FILE EXPECTED_PAIR_ORIENTATIONS IS_BISULFITE_SEQUENCED "
            "MAX_INSERT_SIZE MAX_RECORDS_IN_RAM METRIC_ACCUMULATION_LEVEL QUIET STOP_AFTER "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("HISTOGRAM_FILE OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectAllelicCounts": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-base-quality minimum-mapping-quality platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CollectArraysVariantCallingMetrics": {
        "value": (
            "CALL_RATE_PF_THRESHOLD COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "MAX_RECORDS_IN_RAM NUM_PROCESSORS QUIET USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("DBSNP INPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectBaseDistributionByCycle": {
        "value": (
            "ALIGNED_READS_ONLY ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "MAX_RECORDS_IN_RAM PF_READS_ONLY QUIET STOP_AFTER USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("CHART_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectF1R2Counts": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals f1r2-max-depth f1r2-median-mq f1r2-min-bq filter-too-short "
            "gcs-max-retries gcs-project-for-requester-pays interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient "
            "library mate-too-distant-length max-clipped-bases max-depth-per-sample "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CollectHiSeqXPfFailMetrics": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE LANE MAX_RECORDS_IN_RAM "
            "NUM_PROCESSORS N_CYCLES PROB_EXPLICIT_READS QUIET USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("BASECALLS_DIR REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectHsMetrics": {
        "value": (
            "ALLELE_FRACTION BAIT_SET_NAME CLIP_OVERLAPPING_READS COMPRESSION_LEVEL COVERAGE_CAP "
            "CREATE_INDEX CREATE_MD5_FILE INCLUDE_INDELS MAX_RECORDS_IN_RAM "
            "METRIC_ACCUMULATION_LEVEL MINIMUM_BASE_QUALITY MINIMUM_MAPPING_QUALITY NEAR_DISTANCE "
            "QUIET SAMPLE_SIZE USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("BAIT_INTERVALS INPUT REFERENCE_SEQUENCE TARGET_INTERVALS"),
        "out": ("OUTPUT PER_BASE_COVERAGE PER_TARGET_COVERAGE THEORETICAL_SENSITIVITY_OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectIlluminaBasecallingMetrics": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE LANE MAX_RECORDS_IN_RAM QUIET "
            "READ_STRUCTURE USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("BARCODES_DIR BASECALLS_DIR INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectIlluminaLaneMetrics": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FILE_EXTENSION MAX_RECORDS_IN_RAM "
            "OUTPUT_PREFIX QUIET READ_STRUCTURE USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("REFERENCE_SEQUENCE RUN_DIRECTORY"),
        "dir": ("OUTPUT_DIRECTORY TMP_DIR"),
    },
    "CollectIndependentReplicateMetrics": {
        "value": (
            "BARCODE_BQ BARCODE_TAG COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "FILTER_UNPAIRED_READS MAX_RECORDS_IN_RAM MINIMUM_BARCODE_BQ MINIMUM_BQ MINIMUM_GQ "
            "MINIMUM_MQ PROGRESS_STEP_INTERVAL QUIET SAMPLE STOP_AFTER USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE VCF"),
        "out": ("MATRIX_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectJumpingLibraryMetrics": {
        "value": (
            "CHIMERA_KB_MIN COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM "
            "MINIMUM_MAPPING_QUALITY QUIET TAIL_LIMIT USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectMultipleMetrics": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE EXTRA_ARGUMENT "
            "FILE_EXTENSION IGNORE_SEQUENCE INCLUDE_UNPAIRED MAX_RECORDS_IN_RAM "
            "METRIC_ACCUMULATION_LEVEL PROGRAM QUIET STOP_AFTER USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("DB_SNP INPUT INTERVALS REFERENCE_SEQUENCE REF_FLAT"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectOxoGMetrics": {
        "value": (
            "COMPRESSION_LEVEL CONTEXTS CONTEXT_SIZE CREATE_INDEX CREATE_MD5_FILE "
            "INCLUDE_NON_PF_READS MAXIMUM_INSERT_SIZE MAX_RECORDS_IN_RAM MINIMUM_INSERT_SIZE "
            "MINIMUM_MAPPING_QUALITY MINIMUM_QUALITY_SCORE QUIET STOP_AFTER USE_JDK_DEFLATER "
            "USE_JDK_INFLATER USE_OQ VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("DB_SNP INPUT INTERVALS REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectQualityYieldMetrics": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FLOW_MODE "
            "INCLUDE_SECONDARY_ALIGNMENTS INCLUDE_SUPPLEMENTAL_ALIGNMENTS MAX_RECORDS_IN_RAM "
            "QUIET STOP_AFTER USE_JDK_DEFLATER USE_JDK_INFLATER USE_ORIGINAL_QUALITIES "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectQualityYieldMetricsFlow": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE INCLUDE_BQ_HISTOGRAM "
            "INCLUDE_SECONDARY_ALIGNMENTS INCLUDE_SUPPLEMENTAL_ALIGNMENTS MAX_RECORDS_IN_RAM "
            "QUIET STOP_AFTER USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY "
            "flow-fill-empty-bins-value flow-ignore-t0-tag"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectQualityYieldMetricsSNVQ": {
        "value": (
            "ALTERNATE_QUALITY_ATTRIBUTE ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX "
            "CREATE_MD5_FILE INCLUDE_BQ_HISTOGRAM INCLUDE_SECONDARY_ALIGNMENTS "
            "INCLUDE_SUPPLEMENTAL_ALIGNMENTS MAX_RECORDS_IN_RAM QUIET STOP_AFTER USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectReadCounts": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short format gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter keep-intervals "
            "keep-read-group keep-reverse-strand-only lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CollectRnaSeqMetrics": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE END_BIAS_BASES "
            "IGNORE_SEQUENCE MAX_RECORDS_IN_RAM METRIC_ACCUMULATION_LEVEL MINIMUM_LENGTH QUIET "
            "RRNA_FRAGMENT_PERCENTAGE STOP_AFTER STRAND_SPECIFICITY USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE REF_FLAT RIBOSOMAL_INTERVALS"),
        "out": ("CHART_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectRrbsMetrics": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE C_QUALITY_THRESHOLD "
            "MAX_MISMATCH_RATE MAX_RECORDS_IN_RAM METRIC_ACCUMULATION_LEVEL MINIMUM_READ_LENGTH "
            "NEXT_BASE_QUALITY_THRESHOLD QUIET SEQUENCE_NAMES USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE"),
        "out": ("METRICS_FILE_PREFIX"),
        "dir": ("TMP_DIR"),
    },
    "CollectSequencingArtifactMetrics": {
        "value": (
            "ASSUME_SORTED COMPRESSION_LEVEL CONTEXTS_TO_PRINT CONTEXT_SIZE CREATE_INDEX "
            "CREATE_MD5_FILE FILE_EXTENSION INCLUDE_DUPLICATES INCLUDE_NON_PF_READS "
            "INCLUDE_UNPAIRED MAXIMUM_INSERT_SIZE MAX_RECORDS_IN_RAM MINIMUM_INSERT_SIZE "
            "MINIMUM_MAPPING_QUALITY MINIMUM_QUALITY_SCORE QUIET STOP_AFTER TANDEM_READS "
            "USE_JDK_DEFLATER USE_JDK_INFLATER USE_OQ VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("DB_SNP INPUT INTERVALS REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectTargetedPcrMetrics": {
        "value": (
            "ALLELE_FRACTION CLIP_OVERLAPPING_READS COMPRESSION_LEVEL COVERAGE_CAP CREATE_INDEX "
            "CREATE_MD5_FILE CUSTOM_AMPLICON_SET_NAME INCLUDE_INDELS MAX_RECORDS_IN_RAM "
            "METRIC_ACCUMULATION_LEVEL MINIMUM_BASE_QUALITY MINIMUM_MAPPING_QUALITY NEAR_DISTANCE "
            "QUIET SAMPLE_SIZE USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("AMPLICON_INTERVALS INPUT REFERENCE_SEQUENCE TARGET_INTERVALS"),
        "out": ("OUTPUT PER_BASE_COVERAGE PER_TARGET_COVERAGE THEORETICAL_SENSITIVITY_OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectUmiPrevalenceMetrics": {
        "value": (
            "BARCODE_BQ BARCODE_TAG COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "FILTER_UNPAIRED_READS MAX_RECORDS_IN_RAM MINIMUM_BARCODE_BQ MINIMUM_MQ "
            "PROGRESS_STEP_INTERVAL QUIET USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY "
            "VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectVariantCallingMetrics": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE GVCF_INPUT MAX_RECORDS_IN_RAM QUIET "
            "THREAD_COUNT USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("DBSNP INPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY TARGET_INTERVALS"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CollectWgsMetrics": {
        "value": (
            "ALLELE_FRACTION COMPRESSION_LEVEL COUNT_UNPAIRED COVERAGE_CAP CREATE_INDEX "
            "CREATE_MD5_FILE INCLUDE_BQ_HISTOGRAM LOCUS_ACCUMULATION_CAP MAX_RECORDS_IN_RAM "
            "MINIMUM_BASE_QUALITY MINIMUM_MAPPING_QUALITY QUIET READ_LENGTH SAMPLE_SIZE "
            "STOP_AFTER USE_FAST_ALGORITHM USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT INTERVALS REFERENCE_SEQUENCE"),
        "out": ("OUTPUT THEORETICAL_SENSITIVITY_OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CombineGenotypingArrayVcfs": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CombineSegmentBreakpoints": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer columns-of-interest create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only labels lenient library mate-too-distant-length "
            "max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality platform-filter-name read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference segments sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CompareBaseQualities": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "VALIDATION_STRINGENCY gcs-max-retries gcs-project-for-requester-pays "
            "round-down-quantized static-quantized-quals throw-on-diff use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("reference"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CompareDuplicatesSpark": {
        "value": (
            "QUIET add-output-vcf-command-line ambig-filter-bases ambig-filter-frac "
            "bam-partition-size black-listed-lanes conf create-output-bam-index "
            "create-output-bam-splitting-index create-output-variant-index disable-read-filter "
            "disable-sequence-dictionary-validation disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends exclude-intervals filter-too-short gcs-max-retries "
            "gcs-project-for-requester-pays input2 interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only library "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "maximum-mapping-quality min-fragment-length min-read-length minimum-mapping-quality "
            "num-reducers platform-filter-name print-summary program-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-name read-validation-stringency sample "
            "sharded-output spark-master spark-verbosity splitting-index-granularity "
            "throw-on-diff use-jdk-deflater use-jdk-inflater use-nio verbosity"
        ),
        "path": ("input reference"),
        "out": ("output output2"),
        "dir": ("output-shard-tmp-dir tmp-dir"),
    },
    "CompareGtcFiles": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("ILLUMINA_BEAD_POOL_MANIFEST_FILE INPUT REFERENCE_SEQUENCE"),
        "dir": ("TMP_DIR"),
    },
    "CompareIntervalLists": {
        "value": (
            "L L2 QUIET gcs-max-retries gcs-project-for-requester-pays use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("reference"),
        "dir": ("tmp-dir"),
    },
    "CompareSAMs": {
        "value": (
            "COMPARE_MQ COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE LENIENT_DUP LENIENT_HEADER "
            "LENIENT_LOW_MQ_ALIGNMENT LENIENT_UNKNOWN_MQ_ALIGNMENT LOW_MQ_THRESHOLD "
            "MAX_RECORDS_IN_RAM QUIET USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY "
            "VERBOSITY"
        ),
        "path": ("REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "ConvertHaplotypeDatabaseToVcf": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "CountBases": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CountReads": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CreateExtendedIlluminaManifest": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE FLAG_DUPLICATES MAX_RECORDS_IN_RAM "
            "QUIET SUPPORTED_BUILD TARGET_BUILD USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": (
            "CLUSTER_FILE DBSNP_FILE INPUT REFERENCE_SEQUENCE SUPPORTED_CHAIN_FILE "
            "SUPPORTED_REFERENCE_FILE"
        ),
        "out": ("BAD_ASSAYS_FILE OUTPUT REPORT_FILE"),
        "dir": ("TMP_DIR"),
    },
    "CreateReadCountPanelOfNormals": {
        "value": (
            "QUIET conf do-impute-zeros extreme-outlier-truncation-percentile "
            "extreme-sample-median-percentile gcs-max-retries gcs-project-for-requester-pays "
            "maximum-chunk-size maximum-zeros-in-interval-percentage "
            "maximum-zeros-in-sample-percentage minimum-interval-median-percentile "
            "number-of-eigensamples program-name spark-master spark-verbosity use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("annotated-intervals input"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "CrosscheckFingerprints": {
        "value": (
            "ALLOW_DUPLICATE_READS CALCULATE_TUMOR_AWARE_RESULTS COMPRESSION_LEVEL CREATE_INDEX "
            "CREATE_MD5_FILE CROSSCHECK_BY CROSSCHECK_MODE EXIT_CODE_WHEN_MISMATCH "
            "EXIT_CODE_WHEN_NO_VALID_CHECKS EXPECT_ALL_GROUPS_TO_MATCH GENOTYPING_ERROR_RATE "
            "INPUT LOD_THRESHOLD LOSS_OF_HET_RATE MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK "
            "MAX_RECORDS_IN_RAM NUM_THREADS OUTPUT_ERRORS_ONLY QUIET REQUIRE_INDEX_FILES "
            "SECOND_INPUT USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": (
            "HAPLOTYPE_MAP INPUT_INDEX_MAP INPUT_SAMPLE_FILE_MAP INPUT_SAMPLE_MAP "
            "REFERENCE_SEQUENCE SAMPLE_INDIVIDUAL_MAP SECOND_INPUT_INDEX_MAP "
            "SECOND_INPUT_SAMPLE_MAP"
        ),
        "out": ("MATRIX_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "DenoiseReadCounts": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays number-of-eigensamples "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("annotated-intervals count-panel-of-normals input"),
        "out": ("denoised-copy-ratios standardized-copy-ratios"),
        "dir": ("tmp-dir"),
    },
    "DepthOfCoverage": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes calculate-coverage-over-genes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer count-type create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "ignore-deletion-sites include-deletions include-ref-n-sites "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-base-quality "
            "max-clipped-bases max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-base-quality min-fragment-length "
            "min-read-length minimum-mapping-quality nBins omit-depth-output-at-each-base "
            "omit-genes-not-entirely-covered-by-traversal omit-interval-statistics "
            "omit-locus-table omit-per-sample-statistics output-format partition-type "
            "platform-filter-name print-base-counts read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output start stop "
            "summary-coverage-threshold use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "DetermineGermlineContigPloidy": {
        "value": (
            "QUIET adamax-beta-1 adamax-beta-2 caller-external-admixing-rate "
            "caller-internal-admixing-rate caller-update-convergence-threshold "
            "convergence-snr-averaging-window convergence-snr-countdown-window "
            "convergence-snr-trigger-threshold disable-annealing disable-caller disable-sampler "
            "exclude-intervals gcs-max-retries gcs-project-for-requester-pays global-psi-scale "
            "initial-temperature input interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals learning-rate "
            "log-emission-samples-per-round log-emission-sampling-median-rel-error "
            "log-emission-sampling-rounds mapping-error-rate max-advi-iter-first-epoch "
            "max-advi-iter-subsequent-epochs max-calling-iters max-training-epochs "
            "mean-bias-standard-deviation min-training-epochs num-thermal-advi-iters "
            "output-prefix sample-psi-scale use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("contig-ploidy-priors model"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "DumpTabixIndex": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("tabix-index"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "ExtractFingerprint": {
        "value": (
            "COMPRESSION_LEVEL CONTAMINATION CREATE_INDEX CREATE_MD5_FILE EXTRACT_CONTAMINATION "
            "EXTRACT_NON_REPRESENTATIVES_TOO INPUT LOCUS_MAX_READS MAX_RECORDS_IN_RAM QUIET "
            "SAMPLE_ALIAS USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("HAPLOTYPE_MAP REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "ExtractIlluminaBarcodes": {
        "value": (
            "BARCODE COMPRESSION_LEVEL COMPRESS_OUTPUTS CREATE_INDEX CREATE_MD5_FILE "
            "DISTANCE_MODE LANE MAX_MISMATCHES MAX_NO_CALLS MAX_RECORDS_IN_RAM "
            "MINIMUM_BASE_QUALITY MINIMUM_QUALITY MIN_MISMATCH_DELTA NUM_PROCESSORS QUIET "
            "READ_STRUCTURE USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("BARCODE_FILE BASECALLS_DIR INPUT_PARAMS_FILE REFERENCE_SEQUENCE"),
        "out": ("METRICS_FILE"),
        "dir": ("OUTPUT_DIR TMP_DIR"),
    },
    "FilterIntervals": {
        "value": (
            "QUIET exclude-intervals extreme-count-filter-maximum-percentile "
            "extreme-count-filter-minimum-percentile extreme-count-filter-percentage-of-samples "
            "gcs-max-retries gcs-project-for-requester-pays interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "low-count-filter-count-threshold low-count-filter-percentage-of-samples "
            "maximum-gc-content maximum-mappability maximum-segmental-duplication-content "
            "minimum-gc-content minimum-mappability minimum-segmental-duplication-content "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("annotated-intervals input"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "FlagStat": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "FlowFeatureMapper": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer copy-attr copy-attr-prefix create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "emit-ref-confidence exclude-intervals exclude-nan-scores filter-too-short "
            "floor-blocks flow-disallow-probs-larger-than-call flow-fill-empty-bins-value "
            "flow-lump-probs flow-matrix-mods flow-probability-scaling-factor "
            "flow-quantization-bins flow-remove-non-single-base-pair-indels "
            "flow-remove-one-zero-probs flow-report-insertion-or-deletion "
            "flow-retain-max-n-probs-base-format flow-symmetric-indel-probs flow-use-t0-tag "
            "gcs-max-retries gcs-project-for-requester-pays gvcf-gq-bands include-dup-reads "
            "include-qc-failed-reads interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-boundary-flows keep-intervals keep-negatives keep-read-group "
            "keep-reverse-strand-only keep-supplementary-alignments lenient library limit-score "
            "mapping-feature mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-read-length max-score max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length min-score minimum-mapping-quality "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output snv-identical-bases snv-identical-bases-after use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "FlowPairHMMAlignReadsToHaplotypes": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line aligner "
            "ambig-filter-bases ambig-filter-frac base-quality-score-threshold black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer concise-output-format "
            "create-output-bam-index create-output-bam-md5 create-output-variant-index "
            "create-output-variant-md5 disable-bam-index-caching "
            "disable-cap-base-qualities-to-map-quality disable-read-filter "
            "disable-sequence-dictionary-validation disable-symmetric-hmm-normalizing "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "dont-use-dragstr-pair-hmm-scores dragstr-het-hom-ratio "
            "enable-dynamic-read-disqualification-for-genotyping exclude-intervals "
            "expected-mismatch-rate-for-read-disqualification filter-too-short "
            "flow-disallow-probs-larger-than-call flow-fill-empty-bins-value flow-lump-probs "
            "flow-matrix-mods flow-probability-scaling-factor flow-quantization-bins "
            "flow-remove-non-single-base-pair-indels flow-remove-one-zero-probs "
            "flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-symmetric-indel-probs flow-use-t0-tag gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-merging-rule "
            "interval-padding interval-set-rule intervals inverted-read-filter "
            "keep-boundary-flows keep-intervals keep-read-group keep-reverse-strand-only lenient "
            "library likelihood-calculation-engine mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality native-pair-hmm-threads "
            "native-pair-hmm-use-double-precision pair-hmm-gap-continuation-penalty "
            "pair-hmm-implementation pcr-indel-model phred-scaled-global-read-mismapping-rate "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency ref-haplotype sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("dragstr-params-path haplotypes input reference sequence-dictionary"),
        "out": ("output pair-hmm-results-file"),
        "dir": ("tmp-dir"),
    },
    "GeneExpressionEvaluation": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals feature-label-key filter-too-short gcs-max-retries "
            "gcs-project-for-requester-pays grouping-type interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-intervals keep-read-group keep-reverse-strand-only lenient "
            "library mate-too-distant-length max-clipped-bases max-fragment-length "
            "max-read-length max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality multi-map-method multi-overlap-method overlap-type "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-strands read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output unspliced use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("gff-file input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "GermlineCNVCaller": {
        "value": (
            "QUIET active-class-padding-hybrid-mode adamax-beta-1 adamax-beta-2 "
            "caller-external-admixing-rate caller-internal-admixing-rate "
            "caller-update-convergence-threshold class-coherence-length cnv-coherence-length "
            "convergence-snr-averaging-window convergence-snr-countdown-window "
            "convergence-snr-trigger-threshold copy-number-posterior-expectation-mode "
            "depth-correction-tau disable-annealing disable-caller disable-sampler "
            "enable-bias-factors exclude-intervals gc-curve-standard-deviation gcs-max-retries "
            "gcs-project-for-requester-pays init-ard-rel-unexplained-variance initial-temperature "
            "input interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-psi-scale interval-set-rule intervals learning-rate "
            "log-emission-samples-per-round log-emission-sampling-median-rel-error "
            "log-emission-sampling-rounds log-mean-bias-standard-deviation mapping-error-rate "
            "max-advi-iter-first-epoch max-advi-iter-subsequent-epochs max-bias-factors "
            "max-calling-iters max-copy-number max-training-epochs min-training-epochs "
            "num-gc-bins num-samples-copy-ratio-approx num-thermal-advi-iters output-prefix "
            "p-active p-alt run-mode sample-psi-scale use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("annotated-intervals contig-ploidy-calls model"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "GetNormalArtifactData": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends error-prob "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality normal-sample platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "GetPileupSummaries": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-depth-per-sample max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality maximum-population-allele-frequency "
            "min-fragment-length min-read-length minimum-mapping-quality "
            "minimum-population-allele-frequency platform-filter-name read-filter "
            "read-filter-expression read-filter-max-hmer read-filter-tag read-filter-tag-comp "
            "read-filter-tag-op read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary variant"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "GetSampleName": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater use-url-encoding verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "GroundTruthReadsBuilder": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac append-sequence base-quality-score-threshold black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-cap-base-qualities-to-map-quality "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-symmetric-hmm-normalizing disable-tool-default-read-filters "
            "discard-non-polyt-softclipped-reads dont-require-soft-clips-both-ends "
            "dont-use-dragstr-pair-hmm-scores dragstr-het-hom-ratio "
            "enable-dynamic-read-disqualification-for-genotyping exclude-intervals "
            "expected-mismatch-rate-for-read-disqualification false-snp-compensation "
            "fill-softclipped-reads fill-trimmed-reads fill-trimmed-reads-Q fill-trimmed-reads-Z "
            "filter-too-short flow-disallow-probs-larger-than-call flow-fill-empty-bins-value "
            "flow-lump-probs flow-matrix-mods flow-probability-scaling-factor "
            "flow-quantization-bins flow-remove-non-single-base-pair-indels "
            "flow-remove-one-zero-probs flow-report-insertion-or-deletion "
            "flow-retain-max-n-probs-base-format flow-symmetric-indel-probs flow-use-t0-tag "
            "gcs-max-retries gcs-project-for-requester-pays gt-no-output "
            "haplotype-output-padding-size include-supp-align interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-boundary-flows keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library likelihood-calculation-engine "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-output-reads "
            "max-read-length max-rq max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-haplotype-score min-haplotype-score-delta min-mq "
            "min-read-length minimum-mapping-quality native-pair-hmm-threads "
            "native-pair-hmm-use-double-precision output-flow-length "
            "pair-hmm-gap-continuation-penalty pair-hmm-implementation pcr-indel-model "
            "phred-scaled-global-read-mismapping-rate platform-filter-name prepend-sequence "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output subsampling-ratio use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": (
            "ancestral-translators-base-path dragstr-params-path input maternal-ref paternal-ref "
            "reference sequence-dictionary"
        ),
        "out": ("output-csv pair-hmm-results-file"),
        "dir": ("tmp-dir"),
    },
    "GroundTruthScorer": {
        "value": (
            "QUIET add-mean-call add-output-sam-program-record add-output-vcf-command-line "
            "ambig-filter-bases ambig-filter-frac base-quality-score-threshold black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-cap-base-qualities-to-map-quality "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-symmetric-hmm-normalizing disable-tool-default-read-filters "
            "dont-require-soft-clips-both-ends dont-use-dragstr-pair-hmm-scores "
            "dragstr-het-hom-ratio enable-dynamic-read-disqualification-for-genotyping "
            "exclude-intervals exclude-zero-flows "
            "expected-mismatch-rate-for-read-disqualification filter-too-short "
            "flow-disallow-probs-larger-than-call flow-fill-empty-bins-value flow-lump-probs "
            "flow-matrix-mods flow-probability-scaling-factor flow-quantization-bins "
            "flow-remove-non-single-base-pair-indels flow-remove-one-zero-probs "
            "flow-report-insertion-or-deletion flow-retain-max-n-probs-base-format "
            "flow-symmetric-indel-probs flow-use-t0-tag gcs-max-retries "
            "gcs-project-for-requester-pays gt-no-output interval-exclusion-padding "
            "interval-merging-rule interval-padding interval-set-rule intervals "
            "inverted-read-filter keep-boundary-flows keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library likelihood-calculation-engine "
            "mate-too-distant-length max-clipped-bases max-fragment-length max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality native-pair-hmm-threads native-pair-hmm-use-double-precision "
            "normalized-score-threshold omit-zeros-from-report pair-hmm-gap-continuation-penalty "
            "pair-hmm-implementation pcr-indel-model phred-scaled-global-read-mismapping-rate "
            "platform-filter-name quality-percentiles read-filter read-filter-expression "
            "read-filter-max-hmer read-filter-tag read-filter-tag-comp read-filter-tag-op "
            "read-group-black-list read-name read-validation-stringency sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater use-softclipped-bases verbosity"
        ),
        "path": (
            "dragstr-params-path features-file genome-prior input reference sequence-dictionary"
        ),
        "out": ("output-csv pair-hmm-results-file report-file"),
        "dir": ("tmp-dir"),
    },
    "GtcToVcf": {
        "value": (
            "ANALYSIS_VERSION_NUMBER COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS EXPECTED_GENDER MAX_RECORDS_IN_RAM "
            "PIPELINE_VERSION QUIET SAMPLE_ALIAS USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": (
            "CLUSTER_FILE EXTENDED_ILLUMINA_MANIFEST FINGERPRINT_GENOTYPES_VCF_FILE GENDER_GTC "
            "ILLUMINA_BEAD_POOL_MANIFEST_FILE INPUT REFERENCE_SEQUENCE"
        ),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "HtsgetReader": {
        "value": (
            "QUIET check-md5 class field format gcs-max-retries gcs-project-for-requester-pays id "
            "intervals notag reader-threads tag url use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "out": ("O"),
        "dir": ("tmp-dir"),
    },
    "IdentifyContaminant": {
        "value": (
            "COMPRESSION_LEVEL CONTAMINATION CREATE_INDEX CREATE_MD5_FILE EXTRACT_CONTAMINATED "
            "INPUT LOCUS_MAX_READS MAX_RECORDS_IN_RAM QUIET SAMPLE_ALIAS USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("HAPLOTYPE_MAP REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "IlluminaBasecallsToFastq": {
        "value": (
            "ADAPTERS_TO_CHECK APPLY_EAMSS_FILTER COMPRESSION_LEVEL COMPRESS_OUTPUTS CREATE_INDEX "
            "CREATE_MD5_FILE DISTANCE_MODE FIRST_TILE FIVE_PRIME_ADAPTER FLOWCELL_BARCODE "
            "FORCE_GC IGNORE_UNEXPECTED_BARCODES INCLUDE_NON_PF_READS LANE MACHINE_NAME "
            "MATCH_BARCODES_INLINE MAX_MISMATCHES MAX_NO_CALLS MAX_READS_IN_RAM_PER_TILE "
            "MAX_RECORDS_IN_RAM MINIMUM_BASE_QUALITY MINIMUM_QUALITY MIN_MISMATCH_DELTA "
            "MIN_TRIMMED_LENGTH NUM_PROCESSORS QUIET READ_NAME_FORMAT READ_STRUCTURE RUN_BARCODE "
            "SORT THREE_PRIME_ADAPTER TILE_LIMIT TRIMMING_QUALITY USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": (
            "BARCODES_DIR BASECALLS_DIR INPUT_PARAMS_FILE MULTIPLEX_PARAMS OUTPUT_PREFIX "
            "REFERENCE_SEQUENCE"
        ),
        "out": ("METRICS_FILE"),
        "dir": ("TMP_DIR"),
    },
    "IlluminaBasecallsToSam": {
        "value": (
            "ADAPTERS_TO_CHECK APPLY_EAMSS_FILTER BARCODE_POPULATION_STRATEGY COMPRESSION_LEVEL "
            "COMPRESS_OUTPUTS CREATE_INDEX CREATE_MD5_FILE DISTANCE_MODE FIRST_TILE "
            "FIVE_PRIME_ADAPTER IGNORE_UNEXPECTED_BARCODES INCLUDE_BARCODE_QUALITY "
            "INCLUDE_BC_IN_RG_TAG INCLUDE_NON_PF_READS LANE LIBRARY_NAME MATCH_BARCODES_INLINE "
            "MAX_MISMATCHES MAX_NO_CALLS MAX_READS_IN_RAM_PER_TILE MAX_RECORDS_IN_RAM "
            "MINIMUM_BASE_QUALITY MINIMUM_QUALITY MIN_MISMATCH_DELTA "
            "MOLECULAR_INDEX_BASE_QUALITY_TAG MOLECULAR_INDEX_TAG NUM_PROCESSORS PLATFORM "
            "PROCESS_SINGLE_TILE QUIET READ_GROUP_ID READ_STRUCTURE RUN_BARCODE RUN_START_DATE "
            "SAMPLE_ALIAS SEQUENCING_CENTER SORT TAG_PER_MOLECULAR_INDEX THREE_PRIME_ADAPTER "
            "TILE_LIMIT USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": (
            "BARCODES_DIR BARCODE_PARAMS BASECALLS_DIR INPUT_PARAMS_FILE LIBRARY_PARAMS "
            "REFERENCE_SEQUENCE"
        ),
        "out": ("METRICS_FILE OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "IntervalListToBed": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET SCORE SORT "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "IntervalListTools": {
        "value": (
            "ACTION BREAK_BANDS_AT_MULTIPLES_OF COMMENT COMPRESSION_LEVEL CREATE_INDEX "
            "CREATE_MD5_FILE DONT_MERGE_ABUTTING INCLUDE_FILTERED INVERT MAX_RECORDS_IN_RAM "
            "OUTPUT_VALUE PADDING QUIET SCATTER_CONTENT SCATTER_COUNT SORT SUBDIVISION_MODE "
            "UNIQUE USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE SECOND_INPUT"),
        "out": ("COUNT_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "LiftOverHaplotypeMap": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("CHAIN INPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "LiftOverIntervalList": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM MIN_LIFTOVER_PCT "
            "QUIET USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("CHAIN INPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY"),
        "out": ("OUTPUT REJECT"),
        "dir": ("TMP_DIR"),
    },
    "LocalAssembler": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac assembly-name black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-gapfill-count min-read-length min-sv-size "
            "min-thin-observations minimum-mapping-quality no-scaffolding pair-padding "
            "platform-filter-name q-min read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output too-many-scaffolds too-many-traversals use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("fasta-file gfa-file"),
        "dir": ("tmp-dir"),
    },
    "MarkIlluminaAdapters": {
        "value": (
            "ADAPTERS ADAPTER_TRUNCATION_LENGTH COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "FIVE_PRIME_ADAPTER MAX_ERROR_RATE_PE MAX_ERROR_RATE_SE MAX_RECORDS_IN_RAM "
            "MIN_MATCH_BASES_PE MIN_MATCH_BASES_SE NUM_ADAPTERS_TO_KEEP PAIRED_RUN "
            "PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN QUIET THREE_PRIME_ADAPTER "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("METRICS OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "MeanQualityByCycle": {
        "value": (
            "ALIGNED_READS_ONLY ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "MAX_RECORDS_IN_RAM PF_READS_ONLY QUIET STOP_AFTER USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("CHART_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "MergeAnnotatedRegions": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference segments sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "MergeAnnotatedRegionsByAnnotation": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac annotations-to-match black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-merge-distance max-read-length "
            "max-soft-clipped-leading-trailing-ratio max-soft-clipped-ratio "
            "max-variants-per-shard maximum-mapping-quality min-fragment-length min-read-length "
            "minimum-mapping-quality output-contig-column output-end-column output-start-column "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference segments sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "MergePedIntoVcf": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY ZCALL_VERSION"
        ),
        "path": ("MAP_FILE ORIGINAL_VCF PED_FILE REFERENCE_SEQUENCE ZCALL_THRESHOLDS_FILE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "ModelSegments": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays genotyping-base-error-rate "
            "genotyping-homozygous-log-ratio-threshold kernel-approximation-dimension "
            "kernel-scaling-allele-fraction kernel-variance-allele-fraction "
            "kernel-variance-copy-ratio maximum-number-of-segments-per-chromosome "
            "maximum-number-of-smoothing-iterations minimum-total-allele-count-case "
            "minimum-total-allele-count-normal minor-allele-fraction-prior-alpha "
            "number-of-burn-in-samples-allele-fraction number-of-burn-in-samples-copy-ratio "
            "number-of-changepoints-penalty-factor number-of-samples-allele-fraction "
            "number-of-samples-copy-ratio number-of-smoothing-iterations-per-fit output-prefix "
            "smoothing-credible-interval-threshold-allele-fraction "
            "smoothing-credible-interval-threshold-copy-ratio use-jdk-deflater use-jdk-inflater "
            "verbosity window-size"
        ),
        "path": ("allelic-counts denoised-copy-ratios normal-allelic-counts segments"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "PathSeqBuildKmers": {
        "value": (
            "QUIET bloom-false-positive-probability gcs-max-retries "
            "gcs-project-for-requester-pays kmer-mask kmer-size kmer-spacing use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("reference"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "PathSeqBuildReferenceTaxonomy": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays genbank-catalog "
            "min-non-virus-contig-length refseq-catalog tax-dump use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("reference"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "PlotDenoisedCopyRatios": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays maximum-copy-ratio "
            "minimum-contig-length output-prefix point-size-copy-ratio use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("denoised-copy-ratios sequence-dictionary standardized-copy-ratios"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "PlotModeledSegments": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays maximum-copy-ratio "
            "minimum-contig-length output-prefix point-size-allele-fraction point-size-copy-ratio "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("allelic-counts denoised-copy-ratios segments sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "PostprocessGermlineCNVCalls": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line allosomal-contig "
            "ambig-filter-bases ambig-filter-frac autosomal-ref-copy-number black-listed-lanes "
            "cloud-index-prefetch-buffer cloud-prefetch-buffer create-output-bam-index "
            "create-output-bam-md5 create-output-variant-index create-output-variant-md5 "
            "disable-bam-index-caching disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "duplication-qs-threshold exclude-intervals filter-too-short gcs-max-retries "
            "gcs-project-for-requester-pays het-deletion-qs-threshold hom-deletion-qs-threshold "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency sample sample-index seconds-between-progress-updates "
            "site-frequency-threshold sites-only-vcf-output use-jdk-deflater use-jdk-inflater "
            "verbosity"
        ),
        "path": (
            "calls-shard-path clustered-breakpoints contig-ploidy-calls input input-intervals-vcf "
            "model-shard-path reference sequence-dictionary"
        ),
        "out": ("output-denoised-copy-ratios output-genotyped-intervals output-genotyped-segments"),
        "dir": ("tmp-dir"),
    },
    "PreprocessIntervals": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac bin-length black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "interval-exclusion-padding interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality padding "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency sample seconds-between-progress-updates "
            "sites-only-vcf-output use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "QualityScoreDistribution": {
        "value": (
            "ALIGNED_READS_ONLY ASSUME_SORTED COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE "
            "INCLUDE_NO_CALLS MAX_RECORDS_IN_RAM PF_READS_ONLY QUIET STOP_AFTER USE_JDK_DEFLATER "
            "USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("CHART_OUTPUT OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "SplitCRAM": {
        "value": (
            "QUIET gcs-max-retries gcs-project-for-requester-pays shard-max-output-count "
            "shard-records use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "SplitIntervals": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-mix-contigs dont-require-soft-clips-both-ends "
            "exclude-intervals extension filter-too-short gcs-max-retries "
            "gcs-project-for-requester-pays interval-exclusion-padding interval-file-num-digits "
            "interval-file-prefix interval-merging-rule interval-padding interval-set-rule "
            "intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-contig-size min-fragment-length min-read-length minimum-mapping-quality "
            "platform-filter-name read-filter read-filter-expression read-filter-max-hmer "
            "read-filter-tag read-filter-tag-comp read-filter-tag-op read-group-black-list "
            "read-name read-validation-stringency sample scatter-count "
            "seconds-between-progress-updates sites-only-vcf-output subdivision-mode "
            "use-jdk-deflater use-jdk-inflater verbosity"
        ),
        "path": ("input reference sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "TagGermlineEvents": {
        "value": (
            "QUIET add-output-sam-program-record add-output-vcf-command-line ambig-filter-bases "
            "ambig-filter-frac black-listed-lanes cloud-index-prefetch-buffer "
            "cloud-prefetch-buffer create-output-bam-index create-output-bam-md5 "
            "create-output-variant-index create-output-variant-md5 disable-bam-index-caching "
            "disable-read-filter disable-sequence-dictionary-validation "
            "disable-tool-default-read-filters dont-require-soft-clips-both-ends endpoint-padding "
            "exclude-intervals filter-too-short gcs-max-retries gcs-project-for-requester-pays "
            "input-call-header interval-exclusion-padding interval-merging-rule interval-padding "
            "interval-set-rule intervals inverted-read-filter keep-intervals keep-read-group "
            "keep-reverse-strand-only lenient library mate-too-distant-length max-clipped-bases "
            "max-fragment-length max-read-length max-soft-clipped-leading-trailing-ratio "
            "max-soft-clipped-ratio max-variants-per-shard maximum-mapping-quality "
            "min-fragment-length min-read-length minimum-mapping-quality platform-filter-name "
            "read-filter read-filter-expression read-filter-max-hmer read-filter-tag "
            "read-filter-tag-comp read-filter-tag-op read-group-black-list read-name "
            "read-validation-stringency reciprocal-threshold sample "
            "seconds-between-progress-updates sites-only-vcf-output use-jdk-deflater "
            "use-jdk-inflater verbosity"
        ),
        "path": ("called-matched-normal-seg-file input reference segments sequence-dictionary"),
        "out": ("output"),
        "dir": ("tmp-dir"),
    },
    "ValidateSamFile": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE IGNORE IGNORE_WARNINGS "
            "INDEX_VALIDATION_STRINGENCY IS_BISULFITE_SEQUENCED MAX_OPEN_TEMP_FILES MAX_OUTPUT "
            "MAX_RECORDS_IN_RAM MODE QUIET SKIP_MATE_VALIDATION USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATE_INDEX VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INPUT REFERENCE_SEQUENCE"),
        "out": ("OUTPUT"),
        "dir": ("TMP_DIR"),
    },
    "VcfToAdpc": {
        "value": (
            "COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE MAX_RECORDS_IN_RAM QUIET "
            "USE_JDK_DEFLATER USE_JDK_INFLATER VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("REFERENCE_SEQUENCE VCF"),
        "out": ("NUM_MARKERS_FILE OUTPUT SAMPLES_FILE"),
        "dir": ("TMP_DIR"),
    },
    "ViewSam": {
        "value": (
            "ALIGNMENT_STATUS COMPRESSION_LEVEL CREATE_INDEX CREATE_MD5_FILE HEADER_ONLY INPUT "
            "MAX_RECORDS_IN_RAM PF_STATUS QUIET RECORDS_ONLY USE_JDK_DEFLATER USE_JDK_INFLATER "
            "VALIDATION_STRINGENCY VERBOSITY"
        ),
        "path": ("INTERVAL_LIST REFERENCE_SEQUENCE"),
        "dir": ("TMP_DIR"),
    },
}

TOOL_OPTIONS: dict[str, dict[str, str]] = {
    tool: {name: kind for kind, names in groups.items() for name in names.split()}
    for tool, groups in _RAW.items()
}
