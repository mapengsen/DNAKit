# P41 序列操作

## 作用

通过专用上游接口补齐 P41 操作，返回 JSON 数据或实际生成的文件路径。

## 处理方法

保留上游算法和默认值，通过明确参数选择算法或输出。拒绝覆盖已有输出，仅在调用时要求可选依赖。

## API

| 模块 | 函数 |
| --- | --- |
| `dnakit.alignment.edlib_p41_functions` | `edlib_prefix_alignment` |
| `dnakit.annotation.ncls_p41_functions` | `ncls_limited_overlaps` |
| `dnakit.phylogeny.pyvolve_p41_functions` | `simulate_pyvolve_sequences`, `pyvolve_state_frequencies` |
| `dnakit.regulation.tangermeme_p41_functions` | `tangermeme_hypothetical_attributions`, `tangermeme_unchunk`, `tangermeme_extract_signal`, `tangermeme_position_entropy` |
| `dnakit.sequence.py2bit_p41_functions` | `twobit_region_composition`, `twobit_masked_blocks` |
| `dnakit.sequence.screed_p41_functions` | `screed_record_at_index`, `screed_record_names`, `screed_is_dna`, `screed_write_record` |
| `dnakit.sequencing_qc.pyslow5_p41_functions` | `slow5_read_records` |

## 输入

输入为明确的本地路径、可 JSON 序列化数组；区间默认是 0-based、左闭右开，特殊上游坐标见函数说明。CLI 程序须另行安装；参数和大小限制以函数文档为准。

## 输出

字典返回状态、解析结果或实际写出的文件路径；绘图接口同时返回序列和特征数据。验证范围见 `paper/write_paper/figure/fig2/task_coverage/coverage_gap/manifests/integration_p41.json`。

### 前缀比对

**示例代码：**

```python
from dnakit.alignment.edlib_p41_functions import edlib_prefix_alignment
r = edlib_prefix_alignment('ACGT', 'ACGTAAA', task='path', verbose=False)
print(r['edit_distance'], r['cigar'])
```

**示例输出：**

```text
0 4=
```

### 限制区间命中数量

**示例代码：**

```python
from dnakit.annotation.ncls_p41_functions import ncls_limited_overlaps
r = ncls_limited_overlaps([0, 5], [4, 9], [1], [3], selection='first', verbose=False)
print(r['query_indices'], r['stored_indices'])
```

**示例输出：**

```text
[0] [0]
```

### FASTX 写出与字符检查

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequence.screed_p41_functions import screed_write_record, screed_is_dna
with TemporaryDirectory() as folder:
    p = Path(folder) / 'read.fa'
    screed_write_record({'name': 'read', 'sequence': 'ACGTN'}, str(p), verbose=False)
    print(p.read_text().strip())
print(screed_is_dna('ACGTN', verbose=False)['is_dna'])
```

**示例输出：**

```text
>read
ACGTN
True
```

### 核苷酸频率

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.phylogeny.pyvolve_p41_functions import pyvolve_state_frequencies
with TemporaryDirectory() as folder:
    p = Path(folder) / 'dna.fa'
    p.write_text('>dna\nACGT\n')
    print(pyvolve_state_frequencies(str(p), verbose=False)['frequencies'])
```

**示例输出：**

```text
{'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
```

### 熵与序列分块还原

**示例代码：**

```python
from dnakit.regulation.tangermeme_p41_functions import tangermeme_position_entropy, tangermeme_unchunk
r = tangermeme_position_entropy([[[0.25], [0.25], [0.25], [0.25]]], verbose=False)
print(round(r['entropy'][0][0], 6))
print(tangermeme_unchunk([[[1, 2]], [[3, 4]]], [4], verbose=False)['sequences'])
```

**示例输出：**

```text
2.0
[[[1.0, 2.0, 3.0, 4.0]]]
```

### SLOW5 读取标识

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequencing_qc.backends import write_slow5_signals
from dnakit.sequencing_qc.pyslow5_p41_functions import slow5_read_records
with TemporaryDirectory() as folder:
    p = Path(folder) / 'read.slow5'
    write_slow5_signals([{'read_id':'r1','signal':[1,2,3],'digitisation':8192,'offset':0,'range':1467,'sampling_rate':4000}], p, verbose=False)
    print(slow5_read_records(str(p), include_signal=False, verbose=False)['records'])
```

**示例输出：**

```text
['r1']
```

### 2bit 掩蔽区间

**示例代码：**

```python
from pathlib import Path
import py2bit
from dnakit.sequence.py2bit_p41_functions import twobit_masked_blocks
# This small fixture is distributed with py2bit; supply your own .2bit in production.
p = Path(py2bit.__file__).parent / 'py2bitTest' / 'foo.2bit'
print(twobit_masked_blocks(str(p), 'chr1', kind='hard', verbose=False)['blocks'])
```

**示例输出：**

```text
[[0, 50], [100, 150]]
```

