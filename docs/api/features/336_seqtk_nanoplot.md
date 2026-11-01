# seqtk 与 NanoPlot：序列文件和长读长命令行操作

这一组函数封装外部命令行程序 [seqtk](https://github.com/lh3/seqtk)（16 个函数，位于 `dnakit.fasta_fastq` 的 `seqtk_p12_cli_functions`）、[NanoPlot](https://github.com/wdecoster/NanoPlot)（1 个函数，位于 `dnakit.sequencing_qc` 的 `nanoplot_p12_cli_functions`）。每个函数运行一个固定的程序或子命令，参数列表经过校验，不经过 shell，带超时和输出大小限制，结果写到明确指定的新路径，绝不覆盖已有输出。测试版本：seqtk 1.5-r133、NanoPlot 1.48.0。导入方式：区域包导出后用 `from <区域包> import ...`，也可直接从模块导入。

这些程序不随 DNAKit 分发：请先安装（seqtk：`conda install -c bioconda seqtk`；NanoPlot：`conda install -c bioconda nanoplot` or `pip install NanoPlot`），再放入 PATH 或传入函数的路径参数。找不到程序时抛出 `BackendUnavailableError` 并给出安装提示；程序出错退出时抛出 `BackendExecutionError` 并附日志末尾。结果是可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 和 `operation`；各小节只列额外的键。生成本页的环境没有安装这些程序，因此每个被检查的示例展示的是在任何机器上都成立的行为（缺少程序）；每个示例后附有安装程序后运行的调用。

## 1) `dnakit.fasta_fastq.seqtk_cutn()`

- **作用：** 在长 N 段处切断序列，或只列出 N 段。
- **处理方法：** 运行 `seqtk cutN -n -p [-g]`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_cutn(input_path: 'str', output_path: 'str', *, min_n_tract: 'int' = 1000, non_n_penalty: 'int' = 10, gaps_only: 'bool' = False, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`min_n_tract`、`non_n_penalty`、`gaps_only`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。除 `outputs` 和 `log_tail` 外没有额外的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_cutn

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_cutn(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_cutn

result = seqtk_cutn(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 2) `dnakit.fasta_fastq.seqtk_dropse()`

- **作用：** 从交错的双端文件中丢弃没有配对的读段。
- **处理方法：** 运行 `seqtk dropse`；相邻两条记录名（去掉 /1 /2）相同才保留。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_dropse(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_dropse

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_dropse(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_dropse

result = seqtk_dropse(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 3) `dnakit.fasta_fastq.seqtk_famask()`

- **作用：** 把 X 编码的掩蔽 FASTA 应用到源 FASTA。
- **处理方法：** 运行 `seqtk famask <src> <mask>`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_famask(source_path: 'str', mask_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`source_path`、`mask_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_famask

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'mask.fa').write_text('x\n')
    (root / 'src.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_famask(
            source_path=str(root / 'src.fa'),
            mask_path=str(root / 'mask.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_famask

result = seqtk_famask(
    source_path='src.fa',
    mask_path='mask.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 4) `dnakit.fasta_fastq.seqtk_fqchk()`

- **作用：** 按读段位置汇总碱基组成和质量。
- **处理方法：** 运行 `seqtk fqchk -q <阈值>`，解析汇总行和逐位置表格。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_fqchk(input_path: 'str', output_path: 'str', *, quality_threshold: 'int' = 20, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`quality_threshold`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：summary, min_len, max_len, avg_len, distinct_quality_values, columns, rows, max_records, row_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_fqchk

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_fqchk(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_fqchk

result = seqtk_fqchk(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 5) `dnakit.fasta_fastq.seqtk_gap()`

- **作用：** 列出长度不小于给定值的 N 缺口位置。
- **处理方法：** 运行 `seqtk gap -l` 并解析三列表。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_gap(input_path: 'str', output_path: 'str', *, min_length: 'int' = 50, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`min_length`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, start, end, max_records, record_count, output_path。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gap

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_gap(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gap

result = seqtk_gap(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 6) `dnakit.fasta_fastq.seqtk_gc()`

- **作用：** 用 X-dropoff 算法找出高 GC（或高 AT）区域。
- **处理方法：** 运行 `seqtk gc [-w] -f -l -x` 并解析区域表。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_gc(input_path: 'str', output_path: 'str', *, at_rich: 'bool' = False, min_fraction: 'float' = 0.6, min_length: 'int' = 20, x_dropoff: 'float' = 10.0, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`at_rich`、`min_fraction`、`min_length`、`x_dropoff`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, start, end, hit_bases, max_records, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gc

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_gc(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gc

result = seqtk_gc(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 7) `dnakit.fasta_fastq.seqtk_hety()`

- **作用：** 统计含 IUPAC 简并码的 FASTA 的区域杂合度。
- **处理方法：** 运行 `seqtk hety -w -t [-m]` 并解析窗口表。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_hety(input_path: 'str', output_path: 'str', *, window_size: 'int' = 50000, start_positions: 'int' = 5, mask_lowercase: 'bool' = False, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`window_size`、`start_positions`、`mask_lowercase`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, start, end, hets_per_window, called_sites, het_sites, max_records, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hety

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_hety(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hety

result = seqtk_hety(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 8) `dnakit.fasta_fastq.seqtk_hpc()`

- **作用：** 对序列做均聚物压缩。
- **处理方法：** 运行 `seqtk hpc`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_hpc(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hpc

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_hpc(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hpc

result = seqtk_hpc(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 9) `dnakit.fasta_fastq.seqtk_kfreq()`

- **作用：** 统计每条序列中某 k-mer 及其单错配邻居的出现次数。
- **处理方法：** 运行 `seqtk kfreq <kmer> <fa>` 并解析表格。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_kfreq(kmer: 'str', input_path: 'str', output_path: 'str', *, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`kmer`、`input_path`、`output_path`；可选关键字参数：`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, length, strand, neighbour_count, exact_count, max_records, record_count, output_path。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_kfreq

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_kfreq(
            kmer='ACGT',
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_kfreq

result = seqtk_kfreq(
    kmer='ACGT',
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 10) `dnakit.fasta_fastq.seqtk_listhet()`

- **作用：** 列出每个杂合（IUPAC 双碱基）位点的位置。
- **处理方法：** 运行 `seqtk listhet` 并解析表格。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_listhet(input_path: 'str', output_path: 'str', *, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, position, base, max_records, record_count, output_path。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_listhet

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_listhet(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_listhet

result = seqtk_listhet(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 11) `dnakit.fasta_fastq.seqtk_mergefa()`

- **作用：** 合并两个等位基因集合 FASTA（可识别杂合）。
- **处理方法：** 运行 `seqtk mergefa [-q -i -m -r -h]`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_mergefa(first_path: 'str', second_path: 'str', output_path: 'str', *, quality_threshold: 'int' = 0, intersection: 'bool' = False, lowercase_on_n: 'bool' = False, random_allele: 'bool' = False, suppress_hets: 'bool' = False, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`first_path`、`second_path`、`output_path`；可选关键字参数：`quality_threshold`、`intersection`、`lowercase_on_n`、`random_allele`、`suppress_hets`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mergefa

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'a.fa').write_text('x\n')
    (root / 'b.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_mergefa(
            first_path=str(root / 'a.fa'),
            second_path=str(root / 'b.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mergefa

result = seqtk_mergefa(
    first_path='a.fa',
    second_path='b.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 12) `dnakit.fasta_fastq.seqtk_mutfa()`

- **作用：** 按位点列表向 FASTA 引入点突变。
- **处理方法：** 运行 `seqtk mutfa <fa> <snp>`。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_mutfa(input_path: 'str', snp_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`snp_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mutfa

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'snps.txt').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_mutfa(
            input_path=str(root / 'in.fa'),
            snp_path=str(root / 'snps.txt'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mutfa

result = seqtk_mutfa(
    input_path='in.fa',
    snp_path='snps.txt',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 13) `dnakit.fasta_fastq.seqtk_randbase()`

- **作用：** 把 FASTA 中的杂合 IUPAC 碱基随机还原为其中一个等位碱基。
- **处理方法：** 运行 `seqtk randbase`（drand48 未设种子，结果可重复）。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_randbase(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_randbase

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_randbase(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_randbase

result = seqtk_randbase(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 14) `dnakit.fasta_fastq.seqtk_seq()`

- **作用：** 对 FASTA/FASTQ 记录做变换：去除标题注释、按质量屏蔽碱基、转大写、选择配对读段的第 1 或第 2 条。
- **处理方法：** 运行 `seqtk seq`，参数 -C / -q / -n / -Q / -U / -1 / -2 / -l 由对应形参生成。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_seq(input_path: 'str', output_path: 'str', *, drop_comments: 'bool' = False, mask_quality_below: 'int | None' = None, mask_character: 'str | None' = None, quality_shift: 'int' = 33, uppercase: 'bool' = False, mate: 'int | None' = None, line_length: 'int | None' = None, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`drop_comments`、`mask_quality_below`、`mask_character`、`quality_shift`、`uppercase`、`mate`、`line_length`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_file, record_count。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_seq

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_seq(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_seq

result = seqtk_seq(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 15) `dnakit.fasta_fastq.seqtk_size()`

- **作用：** 统计 FASTA/FASTQ 文件的序列数和碱基数。
- **处理方法：** 运行 `seqtk size` 并解析输出的两列。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_size(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。除 `outputs` 和 `log_tail` 外没有额外的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_size

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_size(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_size

result = seqtk_size(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 16) `dnakit.fasta_fastq.seqtk_telo()`

- **作用：** 识别组装序列或长读段两端的端粒重复。
- **处理方法：** 运行 `seqtk telo -m -p -d -s` 并解析区间表和 stderr 中的总数。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `seqtk_telo(input_path: 'str', output_path: 'str', *, motif: 'str' = 'CCCTAA', penalty: 'int' = 1, max_drop: 'int' = 2000, min_score: 'int' = 300, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 seqtk（测试版本 1.5-r133）位于 PATH 中或通过 `seqtk_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_path`、`output_path`；可选关键字参数：`motif`、`penalty`、`max_drop`、`min_score`、`max_records`、`seqtk_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：records, seqid, start, end, length, max_records, record_count, telomere_bases, total_bases。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_telo

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_telo(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
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
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_telo

result = seqtk_telo(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 17) `dnakit.sequencing_qc.nanoplot_plots()`

- **作用：** 绘制 NanoPlot 的全部长读长图（长度直方图、产量、长度-质量散点、比对图、时间/通道图）并返回统计和原始表。
- **处理方法：** 运行一次 `NanoPlot`，按固定名称归类生成的图文件，解析 NanoStats.txt，可选导出 `--raw` 逐读段表。详细参数含义见 API 签名与英文页 Method 说明。
- **API：** `nanoplot_plots(input_paths: 'Sequence[str]', output_dir: 'str', *, input_type: 'str' = 'fastq', plot_names: 'Sequence[str] | None' = None, bivariate_plots: 'Sequence[str] | None' = None, static_images: 'bool' = True, image_formats: 'Sequence[str] | None' = None, log_length: 'bool' = False, percent_quality: 'bool' = False, aligned_length: 'bool' = False, min_length: 'int | None' = None, max_length: 'int | None' = None, min_quality: 'float | None' = None, downsample: 'int | None' = None, drop_outliers: 'bool' = False, no_supplementary: 'bool' = False, show_n50: 'bool | None' = None, color: 'str | None' = None, colormap: 'str | None' = None, title: 'str | None' = None, dpi: 'int | None' = None, raw_table: 'bool' = False, threads: 'int' = 1, timeout_seconds: 'float' = 86400.0, nanoplot_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 NanoPlot（测试版本 1.48.0）位于 PATH 中或通过 `nanoplot_path` 指定；生成本页的环境没有安装它。下面的示例展示不依赖它的行为：缺少程序时抛出的错误。
- **输入：** 必填：`input_paths`、`output_dir`；可选关键字参数：`input_type`、`plot_names`、`bivariate_plots`、`static_images`、`image_formats`、`log_length`、`percent_quality`、`aligned_length`、`min_length`、`max_length`、`min_quality`、`downsample`、`drop_outliers`、`no_supplementary`、`show_n50`、`color`、`colormap`、`title`、`dpi`、`raw_table`、`threads`、`timeout_seconds`、`nanoplot_path`。输入文件必须已存在，输出路径必须尚不存在。
- **输出：** 通用结果字典。额外的键：output_dir, plots, missing_plots, statistics, number_of_reads, number_of_bases, median_read_length, mean_read_length, mean_qual, median_qual, stats_file, report_file, raw_table, raw_table_file, raw_table_columns, raw_table_rows。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.sequencing_qc.nanoplot_p12_cli_functions import nanoplot_plots

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fastq').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        nanoplot_plots(
            input_paths=[str(root / 'in.fastq')],
            output_dir=str(root / 'out'),
            nanoplot_path='/nonexistent/NanoPlot',
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
from dnakit.sequencing_qc.nanoplot_p12_cli_functions import nanoplot_plots

result = nanoplot_plots(
    input_paths=['in.fastq'],
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```
