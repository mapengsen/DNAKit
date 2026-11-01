# DeepVariant：变异调用、流程各阶段、训练与辅助工具

这一组函数封装外部命令行程序 [DeepVariant](https://github.com/google/deepvariant)（15 个函数，位于 `dnakit.variants` 的 `deepvariant_p12_cli_functions`）。每个函数运行一个固定的程序或子命令，参数列表经过校验，不经过 shell，带超时和输出大小限制，结果写到明确指定的新路径，绝不覆盖已有输出。测试版本：DeepVariant 1.10.0。导入方式：区域包导出后用 `from <区域包> import ...`，也可直接从模块导入。

这些程序不随 DNAKit 分发：请先安装（DeepVariant：a DeepVariant container image (Docker or Apptainer) whose programs are made executable on the host），再放入 PATH 或传入函数的路径参数。找不到程序时抛出 `BackendUnavailableError` 并给出安装提示；程序出错退出时抛出 `BackendExecutionError` 并附日志末尾。结果是可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 和 `operation`；各小节只列额外的键。生成本页的环境没有安装这些程序，因此每个被检查的示例展示的是在任何机器上都成立的行为（缺少程序）；每个示例后附有安装程序后运行的调用。

## 1) `dnakit.variants.deepvariant_call_variants()`

- **作用：** 在 pileup 样例上运行 CNN 基因型推断。
- **处理方法：** 运行 `call_variants`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_call_variants(examples_path: 'str', output_path: 'str', checkpoint: 'str', *, batch_size: 'int | None' = None, num_readers: 'int | None' = None, writer_threads: 'int | None' = None, execution_hardware: 'str | None' = None, include_debug_info: 'bool | None' = None, allow_empty_examples: 'bool | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`examples_path`、`output_path`、`checkpoint`；可选关键字参数：`batch_size`、`num_readers`、`writer_threads`、`execution_hardware`、`include_debug_info`、`allow_empty_examples`、`extra_flags`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, input_example_files。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 比较两种标注器给出的变异标签。
- **处理方法：** 运行 `compare_labelers`（源码工具，需自行构建）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_compare_labelers(labeled_examples_a: 'str', labeled_examples_b: 'str', output_path: 'str', *, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`labeled_examples_a`、`labeled_examples_b`、`output_path`；可选关键字参数：`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, lines。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 把带标签的训练样例转回真实基因型 VCF。
- **处理方法：** 运行 `labeled_examples_to_vcf`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_labeled_examples_to_vcf(reference_path: 'str', examples_path: 'str', output_vcf: 'str', *, sample_name: 'str | None' = None, max_records: 'int | None' = None, allow_unlabeled_examples: 'bool | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`reference_path`、`examples_path`、`output_vcf`；可选关键字参数：`sample_name`、`max_records`、`allow_unlabeled_examples`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 生成 DeepVariant 的 pileup 样例（推断、训练、候选位点扫描；germline/trio/多样本/体细胞/泛基因组）。
- **处理方法：** 按 `program` 运行对应的 make_examples 可执行文件，`mode` 选择 calling/training/candidate_sweep。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_make_examples(program: 'str', reference_path: 'str', reads: 'Sequence[str]', *, examples_path: 'str | None' = None, mode: 'str' = 'calling', sample_names: 'Sequence[str] | None' = None, regions: 'Sequence[str] | None' = None, confident_regions: 'str | None' = None, truth_variants: 'str | None' = None, candidate_positions: 'str | None' = None, pangenome_path: 'str | None' = None, gvcf_path: 'str | None' = None, task: 'int | None' = None, checkpoint: 'str | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`program`、`reference_path`、`reads`；可选关键字参数：`examples_path`、`mode`、`sample_names`、`regions`、`confident_regions`、`truth_variants`、`candidate_positions`、`pangenome_path`、`gvcf_path`、`task`、`checkpoint`、`extra_flags`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：examples_files, candidate_positions_files。
- **示例代码：**

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

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
```

**安装程序后的用法：**

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

- **作用：** 把模型预测转换为 VCF/gVCF。
- **处理方法：** 运行 `postprocess_variants`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_postprocess_variants(reference_path: 'str', call_variants_output: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, nonvariant_site_tfrecord_path: 'str | None' = None, small_model_cvo_records: 'str | None' = None, sample_name: 'str | None' = None, cpus: 'int | None' = None, regions: 'Sequence[str] | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, process_somatic: 'bool | None' = None, extra_flags: 'Mapping[str, str | int | float | bool] | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`reference_path`、`call_variants_output`、`output_vcf`；可选关键字参数：`output_gvcf`、`nonvariant_site_tfrecord_path`、`small_model_cvo_records`、`sample_name`、`cpus`、`regions`、`haploid_contigs`、`par_regions_bed`、`process_somatic`、`extra_flags`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf, output_gvcf。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 规范化真值 VCF 的等位基因和基因型。
- **处理方法：** 运行 `preprocess_truth`（源码工具，需自行构建）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_preprocess_truth(truth_vcf: 'str', output_vcf: 'str', *, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`truth_vcf`、`output_vcf`；可选关键字参数：`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 从评估指标文件中提取各检查点的 F1。
- **处理方法：** 运行 `print_f1`（源码工具，需自行构建）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_print_f1(metrics_dir: 'str', output_path: 'str', *, timeout_seconds: 'float' = 3600.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`metrics_dir`、`output_path`；可选关键字参数：`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。除 `outputs` 和 `log_tail` 外没有额外的键。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用 `run_deepvariant` 从比对读段调用胚系变异，可输出读段支持的单倍型定相和 VCF 统计报告。
- **处理方法：** 运行 `run_deepvariant`，`phase_vcf`、`vcf_stats_report` 对应同名标志。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_run(model_type: 'str', reference_path: 'str', reads_path: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name: 'str | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, runtime_report: 'bool' = False, customized_model: 'str | None' = None, customized_model_json: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, phase_vcf: 'bool | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, postprocess_cpus: 'int | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`model_type`、`reference_path`、`reads_path`、`output_vcf`；可选关键字参数：`output_gvcf`、`num_shards`、`regions`、`sample_name`、`haploid_contigs`、`par_regions_bed`、`intermediate_results_dir`、`logging_dir`、`runtime_report`、`customized_model`、`customized_model_json`、`disable_small_model`、`customized_small_model`、`phase_vcf`、`vcf_stats_report`、`report_title`、`postprocess_cpus`、`make_examples_extra_args`、`call_variants_extra_args`、`postprocess_variants_extra_args`、`dry_run`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf, output_gvcf, stats_report。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用 DeepSomatic 从肿瘤/正常读段调用体细胞变异。
- **处理方法：** 运行 `run_deepsomatic`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_run_deepsomatic(model_type: 'str', reference_path: 'str', reads_tumor: 'str', output_vcf: 'str', sample_name_tumor: 'str', *, reads_normal: 'str | None' = None, sample_name_normal: 'str | None' = None, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model: 'str | None' = None, customized_model_json: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, process_somatic: 'bool | None' = None, pon_filtering: 'str | None' = None, use_default_pon_filtering: 'bool | None' = None, use_candidate_partition: 'bool | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`model_type`、`reference_path`、`reads_tumor`、`output_vcf`、`sample_name_tumor`；可选关键字参数：`reads_normal`、`sample_name_normal`、`output_gvcf`、`num_shards`、`regions`、`intermediate_results_dir`、`logging_dir`、`customized_model`、`customized_model_json`、`disable_small_model`、`customized_small_model`、`process_somatic`、`pon_filtering`、`use_default_pon_filtering`、`use_candidate_partition`、`vcf_stats_report`、`report_title`、`make_examples_extra_args`、`call_variants_extra_args`、`postprocess_variants_extra_args`、`dry_run`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf, output_gvcf, stats_report。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用 `run_deeptrio` 对三人家系或二人组联合调用变异。
- **处理方法：** 运行 `run_deeptrio`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_run_deeptrio(model_type: 'str', reference_path: 'str', reads_child: 'str', output_vcf_child: 'str', sample_name_child: 'str', *, reads_parent1: 'str | None' = None, output_vcf_parent1: 'str | None' = None, sample_name_parent1: 'str | None' = None, reads_parent2: 'str | None' = None, output_vcf_parent2: 'str | None' = None, sample_name_parent2: 'str | None' = None, output_gvcf_child: 'str | None' = None, output_gvcf_parent1: 'str | None' = None, output_gvcf_parent2: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model_child: 'str | None' = None, customized_model_parent: 'str | None' = None, disable_small_model: 'bool' = False, use_candidate_partition: 'bool | None' = None, vcf_stats_report: 'bool' = False, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`model_type`、`reference_path`、`reads_child`、`output_vcf_child`、`sample_name_child`；可选关键字参数：`reads_parent1`、`output_vcf_parent1`、`sample_name_parent1`、`reads_parent2`、`output_vcf_parent2`、`sample_name_parent2`、`output_gvcf_child`、`output_gvcf_parent1`、`output_gvcf_parent2`、`num_shards`、`regions`、`intermediate_results_dir`、`logging_dir`、`customized_model_child`、`customized_model_parent`、`disable_small_model`、`use_candidate_partition`、`vcf_stats_report`、`make_examples_extra_args`、`call_variants_extra_args`、`postprocess_variants_extra_args`、`dry_run`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcfs, child, output_gvcfs。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用真值 VCF 生成“完美分类器”的 VCF，衡量候选位点上限。
- **处理方法：** 运行 `run_oracle_inference`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_run_oracle_inference(model_type: 'str', reference_path: 'str', reads_path: 'str', truth_variants: 'str', confident_regions: 'str', output_vcf: 'str', *, labeler_algorithm: 'str | None' = None, haploid_contigs: 'str | None' = None, par_regions_bed: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`model_type`、`reference_path`、`reads_path`、`truth_variants`、`confident_regions`、`output_vcf`；可选关键字参数：`labeler_algorithm`、`haploid_contigs`、`par_regions_bed`、`num_shards`、`regions`、`sample_name`、`intermediate_results_dir`、`logging_dir`、`make_examples_extra_args`、`dry_run`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用泛基因组感知的 DeepVariant 调用变异。
- **处理方法：** 运行 `run_pangenome_aware_deepvariant`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_run_pangenome_aware(model_type: 'str', reference_path: 'str', reads_path: 'str', pangenome_path: 'str', output_vcf: 'str', *, output_gvcf: 'str | None' = None, num_shards: 'int' = 1, regions: 'Sequence[str] | None' = None, sample_name_reads: 'str | None' = None, sample_name_pangenome: 'str | None' = None, ref_name_pangenome: 'str | None' = None, intermediate_results_dir: 'str | None' = None, logging_dir: 'str | None' = None, customized_model: 'str | None' = None, disable_small_model: 'bool' = False, customized_small_model: 'str | None' = None, vcf_stats_report: 'bool' = False, report_title: 'str | None' = None, postprocess_cpus: 'int | None' = None, gbz_shared_memory_name: 'str | None' = None, gbz_shared_memory_size_gb: 'int | None' = None, make_examples_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, call_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, postprocess_variants_extra_args: 'Mapping[str, str | int | float | bool] | None' = None, dry_run: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`model_type`、`reference_path`、`reads_path`、`pangenome_path`、`output_vcf`；可选关键字参数：`output_gvcf`、`num_shards`、`regions`、`sample_name_reads`、`sample_name_pangenome`、`ref_name_pangenome`、`intermediate_results_dir`、`logging_dir`、`customized_model`、`disable_small_model`、`customized_small_model`、`vcf_stats_report`、`report_title`、`postprocess_cpus`、`gbz_shared_memory_name`、`gbz_shared_memory_size_gb`、`make_examples_extra_args`、`call_variants_extra_args`、`postprocess_variants_extra_args`、`dry_run`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_vcf, output_gvcf, stats_report。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 把 pileup 样例渲染为图像。
- **处理方法：** 运行 `show_examples`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_show_examples(examples_path: 'str', output_prefix: 'str', *, example_info_json: 'str | None' = None, vcf: 'str | None' = None, image_type: 'str' = 'channels', regions: 'Sequence[str] | None' = None, num_records: 'int | None' = None, annotation: 'bool | None' = None, truth_labels: 'bool | None' = None, scale: 'int | None' = None, max_examples_to_scan: 'int | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`examples_path`、`output_prefix`；可选关键字参数：`example_info_json`、`vcf`、`image_type`、`regions`、`num_records`、`annotation`、`truth_labels`、`scale`、`max_examples_to_scan`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：image_files, image_count。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 用带标签的样例训练 DeepVariant CNN 模型。
- **处理方法：** 运行 `train`，配置项用 `--config.<名称>` 传入。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_train(train_dataset_pbtxt: 'str', tune_dataset_pbtxt: 'str', experiment_dir: 'str', *, config_name: 'str' = 'base', config_file: 'str' = 'dv_config.py', strategy: 'str' = 'mirrored', config_overrides: 'Mapping[str, str | int | float | bool] | None' = None, debug: 'bool' = False, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`train_dataset_pbtxt`、`tune_dataset_pbtxt`、`experiment_dir`；可选关键字参数：`config_name`、`config_file`、`strategy`、`config_overrides`、`debug`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：experiment_dir, files。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

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

- **作用：** 为单样本 VCF 生成 HTML 统计报告。
- **处理方法：** 运行 `vcf_stats_report`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `deepvariant_vcf_stats_report(input_vcf: 'str', output_base: 'str', *, title: 'str | None' = None, max_records: 'int | None' = None, timeout_seconds: 'float' = 86400.0, deepvariant_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DeepVariant（测试版本 1.10.0）位于 PATH 中或通过 `deepvariant_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_vcf`、`output_base`；可选关键字参数：`title`、`max_records`、`timeout_seconds`、`deepvariant_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：report_file。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**安装程序后的用法：**

```python
from dnakit.variants.deepvariant_p12_cli_functions import deepvariant_vcf_stats_report

result = deepvariant_vcf_stats_report(
    input_vcf='calls.vcf',
    output_base='report',
    verbose=False,
)
print(result['status'], result['outputs'])
```
