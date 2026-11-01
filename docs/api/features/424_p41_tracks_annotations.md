# P41 轨道与注释

## 作用

通过专用上游接口补齐 P41 操作，返回 JSON 数据或实际生成的文件路径。

## 处理方法

保留上游算法和默认值，通过明确参数选择算法或输出。拒绝覆盖已有输出，仅在调用时要求可选依赖。

## API

| 模块 | 函数 |
| --- | --- |
| `dnakit.annotation.prokka_p41_cli_functions` | `annotate_prokka_rnammer` |
| `dnakit.annotation.pybigwig_p41_functions` | `inspect_bigwig_file` |
| `dnakit.databases.gget_p41_functions` | `gget_download_reference`, `gget_isoform_structure_map` |
| `dnakit.epigenome.epigenome_p41_cli_functions` | `pairtools_tile_duplication_stats`, `macs3_bedgraph_cutoff_report`, `macs3_duplicate_statistics`, `bam_region_count_matrix`, `bigwig_region_score_matrix` |
| `dnakit.molecular_lab.dna_features_viewer_p41_functions` | `read_gff_feature_record`, `export_feature_map_record`, `plot_feature_map` |
| `dnakit.variants.cyvcf2_p41_functions` | `cyvcf2_ibd_bins`, `cyvcf2_informative_sites`, `cyvcf2_export_text`, `vcf_call_is_filtered` |

## 输入

输入为明确的本地路径、可 JSON 序列化数组；区间默认是 0-based、左闭右开，特殊上游坐标见函数说明。CLI 程序须另行安装；参数和大小限制以函数文档为准。

## 输出

字典返回状态、解析结果或实际写出的文件路径；绘图接口同时返回序列和特征数据。验证范围见 `paper/write_paper/figure/fig2/task_coverage/coverage_gap/manifests/integration_p41.json`。

GFF 解析需要 bcbio-gff。gget 为联网操作，普通测试不请求真实服务。Prokka/RNAmmer 和 bigBed autoSQL 尚缺真实端到端验证，不能据此声称科学结果已验收。

### 注释特征导出

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from Bio import SeqIO
from dnakit.molecular_lab.dna_features_viewer_p41_functions import export_feature_map_record
with TemporaryDirectory() as folder:
    p = Path(folder) / 'dna.gb'
    export_feature_map_record('ATGAAATAG', [{'start':0,'end':9,'label':'gene'}], str(p), verbose=False)
    r = SeqIO.read(p, 'genbank')
    print(str(r.seq), len(r.features))
```

**示例输出：**

```text
ATGAAATAG 1
```

### 基因型过滤标记

**示例代码：**

```python
from dnakit.variants.cyvcf2_p41_functions import vcf_call_is_filtered
print(vcf_call_is_filtered('0/1', filters=['PASS'], verbose=False)['is_filtered'])
print(vcf_call_is_filtered('0/1', filters=['q10'], verbose=False)['is_filtered'])
```

**示例输出：**

```text
False
True
```

### bigWig 类型检查

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
import pyBigWig
from dnakit.annotation.pybigwig_p41_functions import inspect_bigwig_file
with TemporaryDirectory() as folder:
    p = Path(folder) / 'scores.bw'
    h = pyBigWig.open(str(p), 'w')
    h.addHeader([('chr1', 100)])
    h.addEntries(['chr1'], [0], ends=[5], values=[1.0])
    h.close()
    r = inspect_bigwig_file(str(p), verbose=False)
    print(r['is_bigwig'], r['is_bigbed'], r['autosql'])
```

**示例输出：**

```text
True False None
```

### MACS3 重复统计

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.epigenome.epigenome_p41_cli_functions import macs3_duplicate_statistics
with TemporaryDirectory() as folder:
    p = Path(folder) / 'reads.bed'
    p.write_text('chr1\t0\t10\ta\t0\t+\n' * 2)
    r = macs3_duplicate_statistics(p, file_format='BED', keep_duplicates=1, verbose=False)
    print(r['statistics']['reads_before'], r['statistics']['reads_after'])
```

**示例输出：**

```text
2 1
```

