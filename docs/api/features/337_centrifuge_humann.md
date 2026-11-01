# Centrifuge 与 HUMAnN：宏基因组参考库与功能谱命令行操作

这一组函数封装外部命令行程序 [Centrifuge](https://github.com/DaehwanKimLab/centrifuge)（9 个函数，位于 `dnakit.microbiome` 的 `centrifuge_p12_cli_functions`）、[HUMAnN](https://github.com/biobakery/humann)（13 个函数，位于 `dnakit.microbiome` 的 `humann_p12_cli_functions`）。每个函数运行一个固定的程序或子命令，参数列表经过校验，不经过 shell，带超时和输出大小限制，结果写到明确指定的新路径，绝不覆盖已有输出。测试版本：Centrifuge 1.0.4、HUMAnN 3.9。导入方式：区域包导出后用 `from <区域包> import ...`，也可直接从模块导入。

这些程序不随 DNAKit 分发：请先安装（Centrifuge：`conda install -c bioconda centrifuge`；HUMAnN：`conda install -c bioconda humann` or `pip install humann`），再放入 PATH 或传入函数的路径参数。找不到程序时抛出 `BackendUnavailableError` 并给出安装提示；程序出错退出时抛出 `BackendExecutionError` 并附日志末尾。结果是可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 和 `operation`；各小节只列额外的键。生成本页的环境没有安装这些程序，因此每个被检查的示例展示的是在任何机器上都成立的行为（缺少程序）；每个示例后附有安装程序后运行的调用。

## 1) `dnakit.microbiome.centrifuge_build()`

- **作用：** 用参考 FASTA 和分类树构建 Centrifuge 分类索引，可同时统计不同 k-mer 数。
- **处理方法：** 运行 `centrifuge-build`，`kmer_count` 对应 `--kmer-count`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_build(reference_files: 'Sequence[str]', index_prefix: 'str', conversion_table: 'str', taxonomy_tree: 'str', *, name_table: 'str | None' = None, size_table: 'str | None' = None, threads: 'int' = 1, offrate: 'int | None' = None, ftab_chars: 'int | None' = None, bmax: 'int | None' = None, bmax_divisor: 'int | None' = None, diff_cover_period: 'int | None' = None, seed: 'int | None' = None, kmer_count: 'int | None' = None, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`reference_files`、`index_prefix`、`conversion_table`、`taxonomy_tree`；可选关键字参数：`name_table`、`size_table`、`threads`、`offrate`、`ftab_chars`、`bmax`、`bmax_divisor`、`diff_cover_period`、`seed`、`kmer_count`、`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：index_prefix, index_files。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_build

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'conv.tsv').write_text('x\n')
    (root / 'nodes.dmp').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_build(
            reference_files=[str(root / 'ref.fa')],
            index_prefix=str(root / 'idx'),
            conversion_table=str(root / 'conv.tsv'),
            taxonomy_tree=str(root / 'nodes.dmp'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_build

result = centrifuge_build(
    reference_files=['ref.fa'],
    index_prefix='idx',
    conversion_table='conv.tsv',
    taxonomy_tree='nodes.dmp',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 2) `dnakit.microbiome.centrifuge_compress_genomes()`

- **作用：** 按分类等级合并相近基因组以压缩参考库。
- **处理方法：** 运行 `centrifuge-compress.pl`（需要 Jellyfish 和 MUMmer）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_compress_genomes(download_dir: 'str', taxonomy_dir: 'str', output_dir: 'str', *, output_name: 'str' = 'compressed', level: 'str' = 'species', mapping_file: 'str | None' = None, threads: 'int' = 1, max_genome_size: 'int' = 50000000, no_compress: 'bool' = False, no_dustmasker: 'bool' = False, centrifuge_dir: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`download_dir`、`taxonomy_dir`、`output_dir`；可选关键字参数：`output_name`、`level`、`mapping_file`、`threads`、`max_genome_size`、`no_compress`、`no_dustmasker`、`centrifuge_dir`、`timeout_seconds`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, fasta, taxid_map, size_table, sequence_count, base_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_compress_genomes

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'download').mkdir()
    (root / 'taxonomy').mkdir()
    (root / 'bin').mkdir()
    try:
        centrifuge_compress_genomes(
            download_dir=str(root / 'download'),
            taxonomy_dir=str(root / 'taxonomy'),
            output_dir=str(root / 'out'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_compress_genomes

result = centrifuge_compress_genomes(
    download_dir='download',
    taxonomy_dir='taxonomy',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 3) `dnakit.microbiome.centrifuge_download()`

- **作用：** 下载 RefSeq/GenBank 基因组、污染序列（UniVec/EmVec）或 NCBI 分类（需联网）。
- **处理方法：** 运行 `centrifuge-download`，可选 RNA、dustmasker 低复杂度屏蔽等。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_download(output_dir: 'str', database: 'str', *, domains: 'Sequence[str] | None' = None, assembly_level: 'str' = 'Complete Genome', refseq_category: 'str | None' = None, taxids: 'Sequence[int] | None' = None, include_rna: 'bool' = False, filter_unplaced: 'bool' = False, mask_low_complexity: 'bool' = False, taxid_in_headers: 'bool' = False, threads: 'int' = 1, download_program: 'str | None' = None, mapping_file: 'str | None' = None, timeout_seconds: 'float' = 86400.0, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`output_dir`、`database`；可选关键字参数：`domains`、`assembly_level`、`refseq_category`、`taxids`、`include_rna`、`filter_unplaced`、`mask_low_complexity`、`taxid_in_headers`、`threads`、`download_program`、`mapping_file`、`timeout_seconds`、`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, mapping_file, mapping_rows, files, output_dir。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_download

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'bin').mkdir()
    try:
        centrifuge_download(
            output_dir=str(root / 'out'),
            database='taxonomy',
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_download

result = centrifuge_download(
    output_dir='out',
    database='taxonomy',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 4) `dnakit.microbiome.centrifuge_inspect()`

- **作用：** 从 Centrifuge 索引提取序列、名称、摘要、转换表、大小表或 k-mer 估计。
- **处理方法：** 运行 `centrifuge-inspect`，由 `output` 选择模式。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_inspect(index_prefix: 'str', output_path: 'str', *, output: 'str' = 'fasta', line_width: 'int' = 60, max_records: 'int' = 10000, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`index_prefix`、`output_path`；可选关键字参数：`output`、`line_width`、`max_records`、`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。除 `outputs` 和 `log_tail` 外没有额外的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_inspect

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_inspect(
            index_prefix=str(root / 'idx'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_inspect

result = centrifuge_inspect(
    index_prefix='idx',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 5) `dnakit.microbiome.centrifuge_kreport()`

- **作用：** 由 Centrifuge 分类输出生成 Kraken 风格报告。
- **处理方法：** 运行 `centrifuge-kreport -x <索引>` 并解析报告。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_kreport(index_prefix: 'str', classification_files: 'Sequence[str]', output_path: 'str', *, min_score: 'int | None' = None, min_length: 'int | None' = None, no_lca: 'bool' = False, show_zeros: 'bool' = False, is_count_table: 'bool' = False, max_records: 'int' = 10000, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`index_prefix`、`classification_files`、`output_path`；可选关键字参数：`min_score`、`min_length`、`no_lca`、`show_zeros`、`is_count_table`、`max_records`、`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, percent, clade_reads, taxon_reads, rank, taxid, name, depth, max_records, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_kreport

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'cls.tsv').write_text('x\n')
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_kreport(
            index_prefix=str(root / 'idx'),
            classification_files=[str(root / 'cls.tsv')],
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_kreport

result = centrifuge_kreport(
    index_prefix='idx',
    classification_files=['cls.tsv'],
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 6) `dnakit.microbiome.centrifuge_remove_empty_sequences()`

- **作用：** 删除序列为空的 FASTA 记录。
- **处理方法：** 运行 `centrifuge-RemoveEmptySequence.pl`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_remove_empty_sequences(input_path: 'str', output_path: 'str', *, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count, base_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_empty_sequences

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_remove_empty_sequences(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_empty_sequences

result = centrifuge_remove_empty_sequences(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 7) `dnakit.microbiome.centrifuge_remove_n()`

- **作用：** 删除 FASTA 记录中的 N 并按 80 列重排。
- **处理方法：** 运行 `centrifuge-RemoveN.pl`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_remove_n(input_path: 'str', output_path: 'str', *, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count, base_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_n

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_remove_n(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_n

result = centrifuge_remove_n(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 8) `dnakit.microbiome.centrifuge_shared_sequence()`

- **作用：** 检测多个基因组的共享序列片段并合并。
- **处理方法：** 运行 `centrifuge-BuildSharedSequence.pl`（需要 Jellyfish 和 nucmer）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_shared_sequence(fasta_files: 'Sequence[str]', output_dir: 'str', *, prefix: 'str' = 'shared', kmer_size: 'int' = 53, kmer_portion: 'float' = 0.01, nucmer_identity: 'int' = 99, overlap: 'int' = 250, fragment: 'bool' = False, centrifuge_dir: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`fasta_files`、`output_dir`；可选关键字参数：`prefix`、`kmer_size`、`kmer_portion`、`nucmer_identity`、`overlap`、`fragment`、`centrifuge_dir`、`timeout_seconds`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, files, records, bases, program_output。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_shared_sequence

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'g1.fa').write_text('x\n')
    (root / 'g2.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_shared_sequence(
            fasta_files=[str(root / 'g1.fa'), str(root / 'g2.fa')],
            output_dir=str(root / 'out'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_shared_sequence

result = centrifuge_shared_sequence(
    fasta_files=['g1.fa', 'g2.fa'],
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 9) `dnakit.microbiome.centrifuge_sort_nt()`

- **作用：** 按分类号对核酸 FASTA 记录排序。
- **处理方法：** 运行 `centrifuge-sort-nt.pl`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `centrifuge_sort_nt(sequence_file: 'str', mapping_files: 'Sequence[str]', output_path: 'str', *, taxid_map_output: 'str | None' = None, unmapped_output: 'str | None' = None, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 Centrifuge（测试版本 1.0.4）位于 PATH 中或通过 `centrifuge_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`sequence_file`、`mapping_files`、`output_path`；可选关键字参数：`taxid_map_output`、`unmapped_output`、`centrifuge_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count, base_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_sort_nt

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'a2t.txt').write_text('x\n')
    (root / 'nt.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_sort_nt(
            sequence_file=str(root / 'nt.fa'),
            mapping_files=[str(root / 'a2t.txt')],
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_sort_nt

result = centrifuge_sort_nt(
    sequence_file='nt.fa',
    mapping_files=['a2t.txt'],
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 10) `dnakit.microbiome.humann_barplot()`

- **作用：** 绘制某个基因家族/通路按物种分层的堆叠条形图。
- **处理方法：** 运行 `humann_barplot`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_barplot(input_path: 'str', focal_feature: 'str', output_image: 'str', *, top_taxa: 'int' = 18, as_genera: 'bool' = False, exclude_unclassified: 'bool' = False, remove_zeros: 'bool' = False, sort: 'Sequence[str] | None' = None, focal_metadata: 'str | None' = None, last_metadata: 'str | None' = None, scaling: 'str' = 'original', units: 'str' = 'generic', dimensions: 'Sequence[float]' = (11.0, 6.0), sample_order_output: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`focal_feature`、`output_image`；可选关键字参数：`top_taxa`、`as_genera`、`exclude_unclassified`、`remove_zeros`、`sort`、`focal_metadata`、`last_metadata`、`scaling`、`units`、`dimensions`、`sample_order_output`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：image, image_bytes, sample_order_file。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_barplot

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_barplot(
            input_path=str(root / 'in.fa'),
            focal_feature='UniRef90_X',
            output_image=str(root / 'plot.png'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_barplot

result = humann_barplot(
    input_path='in.fa',
    focal_feature='UniRef90_X',
    output_image='plot.png',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 11) `dnakit.microbiome.humann_build_custom_database()`

- **作用：** 由基因序列（可按属过滤）构建自定义搜索数据库。
- **处理方法：** 运行 `humann_build_custom_database`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_build_custom_database(fasta_path: 'str', output_dir: 'str', *, id_mapping: 'str | None' = None, taxonomic_profile: 'str | None' = None, database_format: 'str' = 'fasta', genus_abundance_threshold: 'float | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`fasta_path`、`output_dir`；可选关键字参数：`id_mapping`、`taxonomic_profile`、`database_format`、`genus_abundance_threshold`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, database_files。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_build_custom_database

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'genes.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_build_custom_database(
            fasta_path=str(root / 'genes.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_build_custom_database

result = humann_build_custom_database(
    fasta_path='genes.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 12) `dnakit.microbiome.humann_expand_cluster()`

- **作用：** 把一个 UniRef90 基因家族扩展到同一 UniRef50 簇的全部成员。
- **处理方法：** 运行 `humann_expand_cluster`（需要 HUMAnN 工具映射数据）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_expand_cluster(input_path: 'str', uniref90: 'str', output_path: 'str', *, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`uniref90`、`output_path`；可选关键字参数：`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_expand_cluster

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_expand_cluster(
            input_path=str(root / 'in.fa'),
            uniref90='UniRef90_X',
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_expand_cluster

result = humann_expand_cluster(
    input_path='in.fa',
    uniref90='UniRef90_X',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 13) `dnakit.microbiome.humann_infer_taxonomy()`

- **作用：** 由 UniRef 簇的最近公共祖先推断未分类基因家族的物种归属。
- **处理方法：** 运行 `humann_infer_taxonomy`，数据库名或 `--devdb` 文件二选一。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_infer_taxonomy(input_path: 'str', output_path: 'str', *, level: 'str' = 'Family', database: 'str | None' = None, database_file: 'str | None' = None, mode: 'str' = 'totals', lca_choice: 'str' = 'humann_lca', threshold: 'float' = 0.001, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`level`、`database`、`database_file`、`mode`、`lca_choice`、`threshold`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows, summary。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_infer_taxonomy

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'tol.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_infer_taxonomy(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            database_file=str(root / 'tol.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_infer_taxonomy

result = humann_infer_taxonomy(
    input_path='in.fa',
    output_path='result.txt',
    database_file='tol.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 14) `dnakit.microbiome.humann_join_tables()`

- **作用：** 把各样本的基因/通路表合并成多样本表。
- **处理方法：** 运行 `humann_join_tables`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_join_tables(input_dir: 'str', output_path: 'str', *, file_name: 'str | None' = None, search_subdirectories: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_dir`、`output_path`；可选关键字参数：`file_name`、`search_subdirectories`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_join_tables

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'input').mkdir()
    (root / 'bin').mkdir()
    try:
        humann_join_tables(
            input_dir=str(root / 'input'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_join_tables

result = humann_join_tables(
    input_dir='input',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 15) `dnakit.microbiome.humann_quantify_alignments()`

- **作用：** 由 SAM/BAM、BLAST m8 比对或基因家族表计算基因家族和通路的丰度与覆盖度。
- **处理方法：** 运行 `humann --input-format ...`，跳过搜索步骤。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_quantify_alignments(input_path: 'str', output_dir: 'str', *, input_format: 'str' = 'blastm8', minpath: 'bool' = True, gap_fill: 'bool' = True, pathways: 'str' = 'metacyc', pathways_database: 'Sequence[str] | None' = None, id_mapping: 'str | None' = None, annotation_gene_index: 'int | None' = None, output_basename: 'str | None' = None, output_format: 'str' = 'tsv', output_max_decimals: 'int | None' = None, remove_stratified_output: 'bool' = False, remove_column_description_output: 'bool' = False, remove_temp_output: 'bool' = False, nucleotide_database: 'str | None' = None, protein_database: 'str | None' = None, threads: 'int' = 1, timeout_seconds: 'float' = 86400.0, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_dir`；可选关键字参数：`input_format`、`minpath`、`gap_fill`、`pathways`、`pathways_database`、`id_mapping`、`annotation_gene_index`、`output_basename`、`output_format`、`output_max_decimals`、`remove_stratified_output`、`remove_column_description_output`、`remove_temp_output`、`nucleotide_database`、`protein_database`、`threads`、`timeout_seconds`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, gene_families, genetable, biom, pathway_abundance, pathway_coverage, file, columns, row_count, rows, biom, file。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_quantify_alignments

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_quantify_alignments(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_quantify_alignments

result = humann_quantify_alignments(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 16) `dnakit.microbiome.humann_reduce_table()`

- **作用：** 用汇总函数把表约简为每个特征一个值。
- **处理方法：** 运行 `humann_reduce_table`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_reduce_table(input_path: 'str', output_path: 'str', *, function: 'str' = 'max', sort_by: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`function`、`sort_by`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_reduce_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_reduce_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_reduce_table

result = humann_reduce_table(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 17) `dnakit.microbiome.humann_regroup_table()`

- **作用：** 把基因家族重新分组为功能组（GO、KO、EC、Pfam、MetaCyc 反应等）。
- **处理方法：** 运行 `humann_regroup_table`，内置或自定义分组二选一。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_regroup_table(input_path: 'str', output_path: 'str', *, groups: 'str | None' = None, custom_groups: 'str | None' = None, reversed_mapping: 'bool' = False, function: 'str' = 'sum', precision: 'int | None' = None, include_ungrouped: 'bool' = True, keep_protected: 'bool' = True, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`groups`、`custom_groups`、`reversed_mapping`、`function`、`precision`、`include_ungrouped`、`keep_protected`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_regroup_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'groups.tsv').write_text('x\n')
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_regroup_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            custom_groups=str(root / 'groups.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_regroup_table

result = humann_regroup_table(
    input_path='in.fa',
    output_path='result.txt',
    custom_groups='groups.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 18) `dnakit.microbiome.humann_rename_table()`

- **作用：** 给基因/通路 ID 加上名称（gloss）。
- **处理方法：** 运行 `humann_rename_table`，内置映射或自定义映射二选一。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_rename_table(input_path: 'str', output_path: 'str', *, names: 'str | None' = None, custom_mapping: 'str | None' = None, simplify: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`names`、`custom_mapping`、`simplify`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_rename_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'names.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_rename_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            custom_mapping=str(root / 'names.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_rename_table

result = humann_rename_table(
    input_path='in.fa',
    output_path='result.txt',
    custom_mapping='names.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 19) `dnakit.microbiome.humann_renorm_table()`

- **作用：** 把基因/通路丰度重新归一化为 CPM 或相对丰度。
- **处理方法：** 运行 `humann_renorm_table`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_renorm_table(input_path: 'str', output_path: 'str', *, units: 'str' = 'cpm', mode: 'str' = 'community', include_special: 'bool' = True, update_sample_names: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`units`、`mode`、`include_special`、`update_sample_names`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_renorm_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_renorm_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_renorm_table

result = humann_renorm_table(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 20) `dnakit.microbiome.humann_split_stratified_table()`

- **作用：** 把分层表拆成分层表和未分层（总量）表。
- **处理方法：** 运行 `humann_split_stratified_table`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_split_stratified_table(input_path: 'str', output_dir: 'str', *, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_dir`；可选关键字参数：`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：stratified, unstratified, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_split_stratified_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_split_stratified_table(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_split_stratified_table

result = humann_split_stratified_table(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 21) `dnakit.microbiome.humann_split_table()`

- **作用：** 把合并表拆成每个样本一张表。
- **处理方法：** 运行 `humann_split_table`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_split_table(input_path: 'str', output_dir: 'str', *, taxonomy_index: 'int | None' = None, taxonomy_level: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_dir`；可选关键字参数：`taxonomy_index`、`taxonomy_level`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：tables, table_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_split_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_split_table(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_split_table

result = humann_split_table(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 22) `dnakit.microbiome.humann_unpack_pathways()`

- **作用：** 把通路丰度拆解为组成基因的丰度。
- **处理方法：** 运行 `humann_unpack_pathways`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `humann_unpack_pathways(gene_table: 'str', pathway_table: 'str', output_path: 'str', *, gene_mapping: 'str | None' = None, pathway_mapping: 'str | None' = None, remove_taxonomy: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 HUMAnN（测试版本 3.9）位于 PATH 中或通过 `humann_dir` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`gene_table`、`pathway_table`、`output_path`；可选关键字参数：`gene_mapping`、`pathway_mapping`、`remove_taxonomy`、`humann_dir`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：table, file, columns, row_count, rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_unpack_pathways

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'genes.tsv').write_text('x\n')
    (root / 'paths.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_unpack_pathways(
            gene_table=str(root / 'genes.tsv'),
            pathway_table=str(root / 'paths.tsv'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_unpack_pathways

result = humann_unpack_pathways(
    gene_table='genes.tsv',
    pathway_table='paths.tsv',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```
