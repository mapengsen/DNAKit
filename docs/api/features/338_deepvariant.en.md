# DeepVariant: variant-calling runners, pipeline stages, training and helper tools

This group of functions wraps the external command-line programs [DeepVariant](https://github.com/google/deepvariant) (15 functions, module `deepvariant_p12_cli_functions` in `dnakit.variants`). Every function runs one fixed program or sub-command with an explicit, validated argument list (no shell, timeout and output-size limits), writes results to explicit new paths and never overwrites an existing output. Tested versions: DeepVariant 1.10.0. Import them with `from <area> import ...` once the area package exports them, or directly from the module.

The programs are not bundled with DNAKit: install them (DeepVariant: a DeepVariant container image (Docker or Apptainer) whose programs are made executable on the host), then put them on PATH or pass the path option of the function. If a program cannot be found, `BackendUnavailableError` is raised with an install hint; if it exits with an error, `BackendExecutionError` is raised with the tail of its log. Results are JSON-serialisable dictionaries with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and `operation`; each section lists only the extra keys. The programs are not installed in the environment used to produce this page, so every checked example shows behaviour that is real on any machine (a missing program); each is followed by the call to run with the program installed.

## 1) `dnakit.variants.deepvariant_call_variants()`

- **Purpose:** Run the CNN genotype inference on pileup examples with ``call_variants``.
- **Method:** ``examples_path`` names the TFRecord example files written by ``deepvariant_make_examples`` (a sharded path such as ``examples.tfrecord@4.gz`` or a glob; at least one file must exist), ``checkpoint`` the model (a saved-model folder or checkpoint prefix, as understood by DeepVariant) and ``output_path`` the new TFRecord (gzip) file of ``CallVariantsOutput`` protos. ``batch_size``, ``num_readers``, ``writer_threads``, ``execution_hardware`` (for example ``cpu`` or ``accelerator``), ``include_debug_info`` and ``allow_empty_examples`` are passed on; ``extra_flags`` holds further ``call_variants`` flags of DeepVariant 1.10.0.
- **API:** `deepvariant_call_variants(examples_path: 'str', output_path: 'str', checkpoint: 'str', *, batch_size: 'int | None' = None, num_readers: 'int | None' = None, writer_threads: 'int | None' = None, execution_hardware: 'str | None' = None, include_debug_info: 'bool | None' = None, allow_empty_examples: 'bool | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `examples_path`, `output_path`, `checkpoint`; optional keyword-only: `batch_size`, `num_readers`, `writer_threads`, `execution_hardware`, `include_debug_info`, `allow_empty_examples`, `extra_flags`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``input_example_files``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_call_variants

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ex.tfrecord.gz').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_call_variants(
            examples_path=str(root / 'ex.tfrecord.gz'),
            output_path=str(root / 'result.txt'),
            checkpoint='model/ckpt',
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_call_variants

result = deepvariant_call_variants(
    examples_path='ex.tfrecord.gz',
    output_path='result.txt',
    checkpoint='model/ckpt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 2) `dnakit.variants.deepvariant_compare_labelers()`

- **Purpose:** Compare the variant labels of two labeled example sets (two labelers) with ``compare_labelers``.
- **Method:** The example sets (TFRecord, one per labeler) must come from the same candidates. The comparison report is written to ``output_path``. ``compare_labelers`` is a source-only tool of the DeepVariant repository (``//deepvariant/labeler:compare_labelers`` in Bazel); build or wrap it first.
- **API:** `deepvariant_compare_labelers(labeled_examples_a: 'str', labeled_examples_b: 'str', output_path: 'str', *, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `labeled_examples_a`, `labeled_examples_b`, `output_path`; optional keyword-only: `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``lines`` (the report lines, at most 200).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_compare_labelers

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'a.tfrecord.gz').write_text('x\n')
    (root / 'b.tfrecord.gz').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_compare_labelers(
            labeled_examples_a=str(root / 'a.tfrecord.gz'),
            labeled_examples_b=str(root / 'b.tfrecord.gz'),
            output_path=str(root / 'result.txt'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_compare_labelers

result = deepvariant_compare_labelers(
    labeled_examples_a='a.tfrecord.gz',
    labeled_examples_b='b.tfrecord.gz',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 3) `dnakit.variants.deepvariant_labeled_examples_to_vcf()`

- **Purpose:** Convert labeled training examples back to a VCF of their true genotypes with ``labeled_examples_to_vcf``.
- **Method:** ``examples_path`` names the labeled example TFRecord file(s) written by ``deepvariant_make_examples`` in training mode. ``sample_name`` is the sample name written in the VCF, ``max_records`` limits the examples read and ``allow_unlabeled_examples`` also accepts examples without a label.
- **API:** `deepvariant_labeled_examples_to_vcf(reference_path: 'str', examples_path: 'str', output_vcf: 'str', *, sample_name: 'str | None' = None, max_records: 'int | None' = None, allow_unlabeled_examples: 'bool | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `reference_path`, `examples_path`, `output_vcf`; optional keyword-only: `sample_name`, `max_records`, `allow_unlabeled_examples`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result key: ``output_vcf``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_labeled_examples_to_vcf

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ex.tfrecord.gz').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_labeled_examples_to_vcf(
            reference_path=str(root / 'ref.fa'),
            examples_path=str(root / 'ex.tfrecord.gz'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_labeled_examples_to_vcf

result = deepvariant_labeled_examples_to_vcf(
    reference_path='ref.fa',
    examples_path='ex.tfrecord.gz',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 4) `dnakit.variants.deepvariant_make_examples()`

- **Purpose:** Generate pileup examples (and candidate positions) with the DeepVariant ``make_examples`` stage programs.
- **Method:** One function covers six operations of the work list through ``program`` and ``mode``. ``program`` selects the executable and how ``reads`` and ``sample_names`` are read: ``make_examples`` (one BAM/CRAM; germline), ``make_examples_somatic`` (tumor, then optionally normal reads), ``multisample_make_examples`` (several BAMs of one run; joined with commas), ``deeptrio_make_examples`` (child, then parent 1 and parent 2; the program lives in the ``deeptrio`` sub-folder of ``deepvariant_dir``) and ``make_examples_pangenome_aware_dv`` (one BAM plus ``pangenome_path``, a GBZ graph). ``mode`` is ``calling`` (examples for inference), ``training`` (examples with labels; needs ``truth_variants`` and ``confident_regions``) or ``candidate_sweep`` (writes ``candidate_positions`` used to cut the genome into intervals with equally many candidates; run ``calling`` afterwards). ``examples_path`` is the output TFRecord (gzip) path; it may be sharded as ``name.tfrecord@4.gz``, in which case ``task`` (0 to 3) selects the shard written by this call. It is required unless ``mode='candidate_sweep'``. ``gvcf_path`` also writes non-variant sites for gVCF output and ``checkpoint`` points to the model whose example format is used. ``regions`` limits the genome. ``extra_flags`` holds further flags of the chosen program (``--flag=value``); only the documented flags of DeepVariant 1.10.0 for that program are accepted.
- **API:** `deepvariant_make_examples(program: 'str', reference_path: 'str', reads: 'Sequence[str]', *, examples_path: 'str | None' = None, mode: 'str' = 'calling', sample_names: 'Sequence[str] | None' = None, regions: 'Sequence[str] | None' = None, confident_regions: 'str | None' = None, truth_variants: 'str | None' = None, candidate_positions: 'str | None' = None, pangenome_path: 'str | None' = None, gvcf_path: 'str | None' = None, task: 'int | None' = None, checkpoint: 'str | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `program`, `reference_path`, `reads`; optional keyword-only: `examples_path`, `mode`, `sample_names`, `regions`, `confident_regions`, `truth_variants`, `candidate_positions`, `pangenome_path`, `gvcf_path`, `task`, `checkpoint`, `extra_flags`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``examples_files`` (the files written) and ``candidate_positions_files``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_make_examples

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ex.tfrecord.gz').write_text('x\n')
    (root / 'reads.bam').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_make_examples(
            program='make_examples',
            reference_path=str(root / 'ref.fa'),
            reads=[str(root / 'reads.bam')],
            examples_path=str(root / 'ex.tfrecord.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_make_examples

result = deepvariant_make_examples(
    program='make_examples',
    reference_path='ref.fa',
    reads=['reads.bam'],
    examples_path='ex.tfrecord.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 5) `dnakit.variants.deepvariant_postprocess_variants()`

- **Purpose:** Convert the model predictions of ``call_variants`` to a VCF (and gVCF) with ``postprocess_variants``.
- **Method:** ``call_variants_output`` is the TFRecord written by ``deepvariant_call_variants`` (sharded paths and globs are passed on unchanged), ``output_vcf`` the new VCF path (``.vcf.gz`` for compressed output) and ``output_gvcf`` a new gVCF path, which also needs ``nonvariant_site_tfrecord_path`` (the gVCF records written by make_examples with ``gvcf_path``). ``small_model_cvo_records`` are the calls of the small model, ``sample_name`` overrides the sample name, ``cpus`` the worker processes, ``regions`` restricts the output, and ``haploid_contigs``/``par_regions_bed`` handle sex chromosomes. ``extra_flags`` holds further ``postprocess_variants`` flags of DeepVariant 1.10.0.
- **API:** `deepvariant_postprocess_variants(reference_path: 'str', call_variants_output: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, nonvariant_site_tfrecord_path: 'str | None' = None, small_model_cvo_records: 'str | None' = None, sample_name: 'str | None' = None, cpus: 'int | None' = None, regions: 'Sequence[str] | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, process_somatic: 'bool | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `reference_path`, `call_variants_output`, `output_vcf`; optional keyword-only: `output_gvcf`, `nonvariant_site_tfrecord_path`, `small_model_cvo_records`, `sample_name`, `cpus`, `regions`, `haploid_contigs`, `par_regions_bed`, `process_somatic`, `extra_flags`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_vcf`` and ``output_gvcf``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_postprocess_variants

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'cvo.tfrecord.gz').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_postprocess_variants(
            reference_path=str(root / 'ref.fa'),
            call_variants_output=str(root / 'cvo.tfrecord.gz'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_postprocess_variants

result = deepvariant_postprocess_variants(
    reference_path='ref.fa',
    call_variants_output='cvo.tfrecord.gz',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 6) `dnakit.variants.deepvariant_preprocess_truth()`

- **Purpose:** Normalise the alleles and genotypes of a phased truth VCF (star alleles, overlapping variants) with ``preprocess_truth``.
- **Method:** Use it only for truth VCFs that are fully phased, such as the T2T Q100 truth sets. ``truth_vcf`` is the input and ``output_vcf`` the new normalised VCF. ``preprocess_truth`` is a source-only tool of the DeepVariant repository (``tools/preprocess_truth.py``); wrap or build it first.
- **API:** `deepvariant_preprocess_truth(truth_vcf: 'str', output_vcf: 'str', *, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `truth_vcf`, `output_vcf`; optional keyword-only: `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result key: ``output_vcf``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_preprocess_truth

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'truth.vcf').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_preprocess_truth(
            truth_vcf=str(root / 'truth.vcf'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_preprocess_truth

result = deepvariant_preprocess_truth(
    truth_vcf='truth.vcf',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 7) `dnakit.variants.deepvariant_print_f1()`

- **Purpose:** Extract the F1 score of every checkpoint from the evaluation metrics files of a training run with ``print_f1``.
- **Method:** ``metrics_dir`` holds the ``*.metrics`` JSON files that ``train`` writes (names ending ``ckpt-<number>.metrics``). ``print_f1`` prints one tab-separated line per file; it is saved in ``output_path`` and returned as ``records`` (``checkpoint``, ``tps_plus_fns`` and ``f1_all``). ``print_f1`` is a source-only tool of the DeepVariant repository (``tools/print_f1.py``); wrap or build it first.
- **API:** `deepvariant_print_f1(metrics_dir: 'str', output_path: 'str', *, timeout_seconds: 'float' = 3600.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `metrics_dir`, `output_path`; optional keyword-only: `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. No extra keys beyond `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_print_f1

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'metrics').mkdir()
    (root / 'bin').mkdir()
    try:
        deepvariant_print_f1(
            metrics_dir=str(root / 'metrics'),
            output_path=str(root / 'result.txt'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_print_f1

result = deepvariant_print_f1(
    metrics_dir='metrics',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 8) `dnakit.variants.deepvariant_run()`

- **Purpose:** Call germline variants from aligned reads with DeepVariant's ``run_deepvariant`` (make_examples, call_variants, postprocess_variants).
- **Method:** ``model_type`` is WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA, MASSEQ or RNASEQ. ``reference_path`` is the FASTA with its ``.fai`` index and ``reads_path`` an aligned, sorted, indexed BAM or CRAM. ``output_vcf`` (and the optional ``output_gvcf``) must not exist. ``phase_vcf=True`` writes read-backed haplotype phasing into the VCF (``--phase_vcf``, long-read PACBIO and ONT_R104 models only; None keeps the program default, which is on for those models); ``vcf_stats_report=True`` also writes the HTML statistics report ``<output_vcf>.visual_report.html`` (``--vcf_stats_report``, titled ``report_title``). Other options: ``num_shards`` (parallel make_examples shards), ``regions`` (region literals or BED paths), ``sample_name``, ``haploid_contigs`` and ``par_regions_bed`` for sex chromosomes, ``intermediate_results_dir``, ``logging_dir`` with ``runtime_report``, ``customized_model`` (a checkpoint) with ``customized_model_json``, ``disable_small_model``, ``customized_small_model``, ``postprocess_cpus`` and mappings ``make_examples_extra_args``, ``call_variants_extra_args`` and ``postprocess_variants_extra_args`` holding the flags of the corresponding stage (see ``deepvariant_make_examples`` for the names). ``dry_run`` only prints the commands the program would run.
- **API:** `deepvariant_run(model_type: 'str', reference_path: 'str', reads_path: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name: 'str | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, runtime_report: 'bool' = False, customized_model: 'str | None' = None, customized_model_json: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, phase_vcf: 'bool | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, postprocess_cpus: 'int | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `model_type`, `reference_path`, `reads_path`, `output_vcf`; optional keyword-only: `output_gvcf`, `num_shards`, `regions`, `sample_name`, `haploid_contigs`, `par_regions_bed`, `intermediate_results_dir`, `logging_dir`, `runtime_report`, `customized_model`, `customized_model_json`, `disable_small_model`, `customized_small_model`, `phase_vcf`, `vcf_stats_report`, `report_title`, `postprocess_cpus`, `make_examples_extra_args`, `call_variants_extra_args`, `postprocess_variants_extra_args`, `dry_run`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report`` (None when not produced).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'reads.bam').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_run(
            model_type='WGS',
            reference_path=str(root / 'ref.fa'),
            reads_path=str(root / 'reads.bam'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run

result = deepvariant_run(
    model_type='WGS',
    reference_path='ref.fa',
    reads_path='reads.bam',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 9) `dnakit.variants.deepvariant_run_deepsomatic()`

- **Purpose:** Call somatic variants from tumor and normal reads with DeepSomatic (``run_deepsomatic``).
- **Method:** ``model_type`` is WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES or, for tumor-only calling, one of the ``*_TUMOR_ONLY`` types (FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY, WGS_TUMOR_ONLY, WES_TUMOR_ONLY, PACBIO_TUMOR_ONLY, ONT_TUMOR_ONLY). Paired types need ``reads_normal`` and ``sample_name_normal``; tumor-only types must not have them. ``process_somatic`` marks the input as somatic (None keeps the program default), ``pon_filtering`` is a panel-of-normals VCF used to filter calls and ``use_default_pon_filtering`` turns on the default panel for tumor-only models. Other options as in ``deepvariant_run``.
- **API:** `deepvariant_run_deepsomatic(model_type: 'str', reference_path: 'str', reads_tumor: 'str', output_vcf: 'str', sample_name_tumor: 'str', *, reads_normal: 'str | None' = None, sample_name_normal: 'str | None' = None, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model: 'str | None' = None, customized_model_json: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, process_somatic: 'bool | None' = None, pon_filtering: 'str | None' = None, use_default_pon_filtering: 'bool | None' = None, use_candidate_partition: 'bool | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `model_type`, `reference_path`, `reads_tumor`, `output_vcf`, `sample_name_tumor`; optional keyword-only: `reads_normal`, `sample_name_normal`, `output_gvcf`, `num_shards`, `regions`, `intermediate_results_dir`, `logging_dir`, `customized_model`, `customized_model_json`, `disable_small_model`, `customized_small_model`, `process_somatic`, `pon_filtering`, `use_default_pon_filtering`, `use_candidate_partition`, `vcf_stats_report`, `report_title`, `make_examples_extra_args`, `call_variants_extra_args`, `postprocess_variants_extra_args`, `dry_run`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_deepsomatic

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ref.fa').write_text('x\n')
    (root / 'tumor.bam').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_run_deepsomatic(
            model_type='WGS_TUMOR_ONLY',
            reference_path=str(root / 'ref.fa'),
            reads_tumor=str(root / 'tumor.bam'),
            output_vcf=str(root / 'calls.vcf.gz'),
            sample_name_tumor='tumor',
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_deepsomatic

result = deepvariant_run_deepsomatic(
    model_type='WGS_TUMOR_ONLY',
    reference_path='ref.fa',
    reads_tumor='tumor.bam',
    output_vcf='calls.vcf.gz',
    sample_name_tumor='tumor',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 10) `dnakit.variants.deepvariant_run_deeptrio()`

- **Purpose:** Call variants jointly in a trio or a duo from aligned reads with ``run_deeptrio`` (DeepTrio).
- **Method:** ``model_type`` is WGS, WES, PACBIO or ONT. The child's reads, VCF and sample name are required; give the reads, VCF and sample name of one parent for a duo, of both parents for a trio (each parent needs all three). ``*_gvcf`` paths are optional; every output must not exist. ``customized_model_child`` and ``customized_model_parent`` replace the default checkpoints, ``use_candidate_partition`` toggles the candidate-sweep partitioning (None keeps the program default), and the mappings ``make_examples_extra_args``, ``call_variants_extra_args`` and ``postprocess_variants_extra_args`` pass further stage flags (names as in ``deepvariant_make_examples`` for the DeepTrio variant). Other options as in ``deepvariant_run``.
- **API:** `deepvariant_run_deeptrio(model_type: 'str', reference_path: 'str', reads_child: 'str', output_vcf_child: 'str', sample_name_child: 'str', *, reads_parent1: 'str | None' = None, output_vcf_parent1: 'str | None' = None, sample_name_parent1: 'str | None' = None, reads_parent2: 'str | None' = None, output_vcf_parent2: 'str | None' = None, sample_name_parent2: 'str | None' = None, output_gvcf_child: 'str | None' = None, output_gvcf_parent1: 'str | None' = None, output_gvcf_parent2: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model_child: 'str | None' = None, customized_model_parent: 'str | None' = None, disable_small_model: 'bool' = False, use_candidate_partition: 'bool | None' = None, vcf_stats_report: 'bool' = False, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `model_type`, `reference_path`, `reads_child`, `output_vcf_child`, `sample_name_child`; optional keyword-only: `reads_parent1`, `output_vcf_parent1`, `sample_name_parent1`, `reads_parent2`, `output_vcf_parent2`, `sample_name_parent2`, `output_gvcf_child`, `output_gvcf_parent1`, `output_gvcf_parent2`, `num_shards`, `regions`, `intermediate_results_dir`, `logging_dir`, `customized_model_child`, `customized_model_parent`, `disable_small_model`, `use_candidate_partition`, `vcf_stats_report`, `make_examples_extra_args`, `call_variants_extra_args`, `postprocess_variants_extra_args`, `dry_run`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_vcfs`` (a mapping from ``child``, ``parent1``, ``parent2`` to the written VCF) and ``output_gvcfs``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_deeptrio

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'child.bam').write_text('x\n')
    (root / 'p1.bam').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_run_deeptrio(
            model_type='WGS',
            reference_path=str(root / 'ref.fa'),
            reads_child=str(root / 'child.bam'),
            output_vcf_child=str(root / 'child.vcf'),
            sample_name_child='child',
            reads_parent1=str(root / 'p1.bam'),
            output_vcf_parent1=str(root / 'p1.vcf'),
            sample_name_parent1='p1',
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_deeptrio

result = deepvariant_run_deeptrio(
    model_type='WGS',
    reference_path='ref.fa',
    reads_child='child.bam',
    output_vcf_child='child.vcf',
    sample_name_child='child',
    reads_parent1='p1.bam',
    output_vcf_parent1='p1.vcf',
    sample_name_parent1='p1',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 11) `dnakit.variants.deepvariant_run_oracle_inference()`

- **Purpose:** Write the VCF that a perfect classifier would produce from labeled examples (``run_oracle_inference``).
- **Method:** The program runs make_examples in training mode against a truth VCF and then converts the labeled examples to a VCF with ``labeled_examples_to_vcf``; it measures the upper bound that candidate generation allows. ``model_type`` is WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA or MASSEQ; ``truth_variants`` is the truth VCF and ``confident_regions`` the BED file of confident regions; ``labeler_algorithm`` selects the labeler (for example ``haplotype_labeler``). Other options as in ``deepvariant_run``.
- **API:** `deepvariant_run_oracle_inference(model_type: 'str', reference_path: 'str', reads_path: 'str', truth_variants: 'str', confident_regions: 'str', output_vcf: 'str', *, labeler_algorithm: 'str | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `model_type`, `reference_path`, `reads_path`, `truth_variants`, `confident_regions`, `output_vcf`; optional keyword-only: `labeler_algorithm`, `haploid_contigs`, `par_regions_bed`, `num_shards`, `regions`, `sample_name`, `intermediate_results_dir`, `logging_dir`, `make_examples_extra_args`, `dry_run`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result key: ``output_vcf``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_oracle_inference

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'conf.bed').write_text('x\n')
    (root / 'reads.bam').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'truth.vcf').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_run_oracle_inference(
            model_type='WGS',
            reference_path=str(root / 'ref.fa'),
            reads_path=str(root / 'reads.bam'),
            truth_variants=str(root / 'truth.vcf'),
            confident_regions=str(root / 'conf.bed'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_oracle_inference

result = deepvariant_run_oracle_inference(
    model_type='WGS',
    reference_path='ref.fa',
    reads_path='reads.bam',
    truth_variants='truth.vcf',
    confident_regions='conf.bed',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 12) `dnakit.variants.deepvariant_run_pangenome_aware()`

- **Purpose:** Call germline variants with pangenome-aware DeepVariant (``run_pangenome_aware_deepvariant``).
- **Method:** The pileup images also contain haplotypes of a pangenome graph. ``model_type`` is WGS or WES, ``pangenome_path`` the pangenome graph in GBZ format and ``reads_path`` aligned, sorted, indexed short reads. ``sample_name_reads``, ``sample_name_pangenome`` and ``ref_name_pangenome`` name the read sample, the pangenome sample and the reference path inside the graph; ``gbz_shared_memory_name`` and ``gbz_shared_memory_size_gb`` use a graph already loaded into shared memory. Other options as in ``deepvariant_run``.
- **API:** `deepvariant_run_pangenome_aware(model_type: 'str', reference_path: 'str', reads_path: 'str', pangenome_path: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name_reads: 'str | None' = None, sample_name_pangenome: 'str | None' = None, ref_name_pangenome: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, postprocess_cpus: 'int | None' = None, gbz_shared_memory_name: 'str | None' = None, gbz_shared_memory_size_gb: 'int | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `model_type`, `reference_path`, `reads_path`, `pangenome_path`, `output_vcf`; optional keyword-only: `output_gvcf`, `num_shards`, `regions`, `sample_name_reads`, `sample_name_pangenome`, `ref_name_pangenome`, `intermediate_results_dir`, `logging_dir`, `customized_model`, `disable_small_model`, `customized_small_model`, `vcf_stats_report`, `report_title`, `postprocess_cpus`, `gbz_shared_memory_name`, `gbz_shared_memory_size_gb`, `make_examples_extra_args`, `call_variants_extra_args`, `postprocess_variants_extra_args`, `dry_run`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_pangenome_aware

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'graph.gbz').write_text('x\n')
    (root / 'reads.bam').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_run_pangenome_aware(
            model_type='WGS',
            reference_path=str(root / 'ref.fa'),
            reads_path=str(root / 'reads.bam'),
            pangenome_path=str(root / 'graph.gbz'),
            output_vcf=str(root / 'calls.vcf.gz'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_run_pangenome_aware

result = deepvariant_run_pangenome_aware(
    model_type='WGS',
    reference_path='ref.fa',
    reads_path='reads.bam',
    pangenome_path='graph.gbz',
    output_vcf='calls.vcf.gz',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 13) `dnakit.variants.deepvariant_show_examples()`

- **Purpose:** Render pileup examples as images with ``show_examples``.
- **Method:** ``examples_path`` names the example TFRecord file(s) (sharded path or glob), ``example_info_json`` the ``example_info.json`` that describes their channels (needed when it is not next to the examples) and ``vcf`` restricts the output to the variants of a VCF. ``output_prefix`` is the path prefix of the image files; ``image_type`` is ``channels`` (a row of channel images), ``RGB`` (stacked channels), ``both`` or ``none`` (no images). ``regions`` filters by 1-based region or BED path, ``num_records`` limits how many examples are drawn, ``annotation`` and ``truth_labels`` add labels, ``scale`` enlarges the images and ``max_examples_to_scan`` bounds the search.
- **API:** `deepvariant_show_examples(examples_path: 'str', output_prefix: 'str', *, example_info_json: 'str | None' = None, vcf: 'str | None' = None, image_type: 'str' = 'channels', regions: 'Sequence[str] | None' = None, num_records: 'int | None' = None, annotation: 'bool | None' = None, truth_labels: 'bool | None' = None, scale: 'int | None' = None, max_examples_to_scan: 'int | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `examples_path`, `output_prefix`; optional keyword-only: `example_info_json`, `vcf`, `image_type`, `regions`, `num_records`, `annotation`, `truth_labels`, `scale`, `max_examples_to_scan`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``image_files`` (new files that start with the output prefix, at most 1000) and ``image_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_show_examples

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ex.tfrecord.gz').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_show_examples(
            examples_path=str(root / 'ex.tfrecord.gz'),
            output_prefix=str(root / 'img'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_show_examples

result = deepvariant_show_examples(
    examples_path='ex.tfrecord.gz',
    output_prefix='img',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 14) `dnakit.variants.deepvariant_train()`

- **Purpose:** Train a DeepVariant CNN model from labeled examples with ``train``.
- **Method:** ``train_dataset_pbtxt`` and ``tune_dataset_pbtxt`` are the dataset configuration files (``.pbtxt``) that list the labeled example TFRecords of the training and tuning sets; ``experiment_dir`` is the folder for checkpoints and logs (it may exist or be new). ``config_file`` and ``config_name`` select the configuration (``dv_config.py:base`` is the DeepVariant default; other names of that file are also accepted) and ``strategy`` is ``mirrored`` (GPU/CPU) or ``tpu``. ``config_overrides`` sets configuration values (``--config.<name>=<value>``) for the names batch_size, num_epochs, num_validation_examples, warmup_steps, limit, steps_per_iter, shuffle_buffer_elements, init_checkpoint, best_checkpoint_metric, optimizer, momentum, use_ema, ema_momentum, optimizer_weight_decay, weight_decay, early_stopping_patience, learning_rate, learning_rate_num_epochs_per_decay, learning_rate_decay_rate, label_smoothing, backbone_dropout_rate, use_mixed_precision, beta_1, beta_2, adaptive_epsilon and include_snp_indel_metrics.
- **API:** `deepvariant_train(train_dataset_pbtxt: 'str', tune_dataset_pbtxt: 'str', experiment_dir: 'str', *, config_name: 'str' = 'base', config_file: 'str' = 'dv_config.py', strategy: 'str' = 'mirrored', config_overrides: 'Mapping[str, str | int | float | bool] | None' = None, debug: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `train_dataset_pbtxt`, `tune_dataset_pbtxt`, `experiment_dir`; optional keyword-only: `config_name`, `config_file`, `strategy`, `config_overrides`, `debug`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``experiment_dir`` and ``files`` (the files of the experiment folder, relative, at most 200).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_train

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'train.pbtxt').write_text('x\n')
    (root / 'tune.pbtxt').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_train(
            train_dataset_pbtxt=str(root / 'train.pbtxt'),
            tune_dataset_pbtxt=str(root / 'tune.pbtxt'),
            experiment_dir=str(root / 'exp'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_train

result = deepvariant_train(
    train_dataset_pbtxt='train.pbtxt',
    tune_dataset_pbtxt='tune.pbtxt',
    experiment_dir='exp',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 15) `dnakit.variants.deepvariant_vcf_stats_report()`

- **Purpose:** Create the visual HTML statistics report of a single-sample VCF with ``vcf_stats_report`` (as ``run_deepvariant --vcf_stats_report``).
- **Method:** ``input_vcf`` needs exactly one sample and a GT field. The report is written to ``<output_base>.visual_report.html`` (which must not exist). ``title`` is the heading (default: the sample name) and ``max_records`` limits the number of VCF records read (default: all).
- **API:** `deepvariant_vcf_stats_report(input_vcf: 'str', output_base: 'str', *, title: 'str | None' = None, max_records: 'int | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of DeepVariant (1.10.0 is the tested version) on PATH or through `deepvariant_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_vcf`, `output_base`; optional keyword-only: `title`, `max_records`, `timeout_seconds`, `deepvariant_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``report_file``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_vcf_stats_report

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'calls.vcf').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        deepvariant_vcf_stats_report(
            input_vcf=str(root / 'calls.vcf'),
            output_base=str(root / 'report'),
            deepvariant_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_vcf_stats_report

result = deepvariant_vcf_stats_report(
    input_vcf='calls.vcf',
    output_base='report',
    verbose=False,
)
print(result['status'], result['outputs'])
```
