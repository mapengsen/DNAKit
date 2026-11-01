# BEDOPS Starch 归档与 convert2bed

`dnakit.annotation`（模块 `bedops_p05_cli_functions`）中的这些函数封装 [BEDOPS](https://github.com/bedops/bedops) 2.4.42 的 Starch 归档程序（`starch`、`unstarch`、`starchcat`、`starchstrip`）和 `convert2bed`，后者把 BAM、SAM、GFF3、GTF、GVF、PSL、RepeatMasker 以及 VCF/WIG 文件转成排序后的 BED 或 Starch。Starch 是 BEDOPS 的压缩 BED 格式，含按染色体划分的数据流和元数据。`starch-diff` 没有封装：它在 2.4.42 发布版中是 Python 2 脚本（`unicode`、`iteritems`），在 Python 3 下无法运行。程序把结果写到标准输出；这些函数把它写入明确指定的新文件，绝不覆盖。

BEDOPS 是外部程序：先安装（例如 `conda install -c bioconda bedops`，或从 https://github.com/bedops/bedops/releases 下载发布包），再把其 `bin` 目录放入 PATH 或传入 `bedops_bin_dir`。找不到程序时抛出带安装提示的 `BackendUnavailableError`；程序出错退出时抛出 `BackendExecutionError`（代码 `BEDOPS_EXECUTION_FAILED`）并附日志末尾，且不会留下不完整的输出。BEDOPS 采用 GPL-2.0 许可，DNAKit 不内置它。

所有结果共有这些字段：`status`（`completed`）、`backend`（`bedops`）、`operation`、`program`、`version`（从程序读取）、`command`、`outputs`（写出文件的路径）、`log_tail` 和 `provenance`（工具、版本、测试版本 2.4.42、可执行文件、引用文献）；下面各函数小节只列出额外字段。

## 示例数据

示例从当前目录读取下面的文件，并把输出写到临时目录。

**示例文件：** `ref.bed`

```text
chr1	100	200	r1
chr1	300	400	r2
chr1	1000	1100	r3
chr2	10	60	r4
```

**示例文件：** `map.bed`

```text
chr1	150	250	m1	2	+
chr1	180	350	m2	4	-
chr1	320	330	m3	6	+
chr2	0	20	m4	10	+
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	101	rs1	A	G	50	PASS	.
chr1	51	rs2	C	T	60	PASS	.
```

**示例文件：** `genes.gtf`

```text
chr1	demo	exon	101	200	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	demo	exon	301	400	.	-	.	gene_id "g2"; transcript_id "t2";
```

**示例文件：** `repeats.out`

```text
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  463  10.2  0.0  1.2  chr1        101     200  (900)  +  L1PA2         LINE/L1                 1   100  (10)    1
```

## 1) `dnakit.annotation.starch_compress()`

- **作用：** 把已排序的 BED 文件压缩成 BEDOPS Starch 归档。
- **处理方法：** 运行 `starch [--bzip2|--gzip] [--note=...] [--omit-signature] [--header] [<unique-tag>] <bed>`，把归档写入 `output_path`；随后用 `unstarch --elements` 读取元素数。
- **API：** `starch_compress(input_path: 'str', output_path: 'str', *, compression: 'str' = 'bzip2', note: 'str | None' = None, omit_signature: 'bool' = False, header: 'bool' = False, unique_tag: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `starch`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的、按 `sort-bed` 顺序排序的 BED 文件。`output_path` 必须尚不存在。`compression` 取 `bzip2`（默认，更小）或 `gzip`（解压更快）。`note` 是写入元数据的单行字符串。`omit_signature` 跳过每条染色体的签名。`header` 接受 UCSC track、SAM 或 VCF 头部。`unique_tag` 是可选标识符。
- **输出：** 通用结果字典，另有 `element_count` 和 `size_bytes`。`outputs` 含 `archive`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    result = starch_compress("ref.bed", archive, compression="gzip", note="documentation example", verbose=False)
    print(result["status"], result["element_count"], list(result["outputs"]))
    unstarch_to_bed(archive, str(Path(tmp) / "back.bed"), verbose=False)
    print("round trip equals input:", (Path(tmp) / "back.bed").read_text() == Path("ref.bed").read_text())
```

- **示例输出：**

```text
completed 4 ['archive']
round trip equals input: True
```

## 2) `dnakit.annotation.unstarch_to_bed()`

- **作用：** 把 Starch 归档解压成 BED，可只取一条染色体。
- **处理方法：** 运行 `unstarch [<chromosome>] <archive>`，把 BED 行写入 `output_path`。
- **API：** `unstarch_to_bed(input_path: 'str', output_path: 'str', *, chromosome: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `unstarch`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的 Starch 归档。`output_path` 必须尚不存在。`chromosome` 可选。
- **输出：** 通用结果字典，另有 `line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    whole = unstarch_to_bed(archive, str(Path(tmp) / "all.bed"), verbose=False)
    one = unstarch_to_bed(archive, str(Path(tmp) / "chr1.bed"), chromosome="chr1", verbose=False)
    print(whole["line_count"], one["line_count"], one["preview"][0].replace("\t", " "))
```

- **示例输出：**

```text
4 3 chr1 100 200 r1
```

## 3) `dnakit.annotation.unstarch_summary()`

- **作用：** 不解压 Starch 归档，直接读取一个元数据值。
- **处理方法：** 运行 `unstarch [<chromosome>] --<statistic> <archive>` 并解析输出。`statistic` 取 `elements`（元素数）、`bases`（总碱基数）、`bases-uniq`（不重复碱基数）、`has-duplicate` 和 `has-nested`（返回布尔值）、`list-chromosomes`（返回列表）、`elements-max-string-length` 或 `list-json`（归档元数据字典）。该函数覆盖工作清单中的：Starch 归档的元素数、总碱基数、不重复碱基数、重复元素检测、嵌套元素检测和染色体列表。
- **API：** `unstarch_summary(input_path: 'str', *, statistic: 'str' = 'elements', chromosome: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `unstarch`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的 Starch 归档。`statistic` 取上述值。`chromosome`（可选）把数值限定在该染色体。
- **输出：** 通用结果字典，另有 `statistic`、`chromosome` 和 `value`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_summary

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    for statistic in ("elements", "bases", "bases-uniq", "has-duplicate", "has-nested", "list-chromosomes"):
        print(statistic, unstarch_summary(archive, statistic=statistic, verbose=False)["value"])
    print("chr1 elements", unstarch_summary(archive, statistic="elements", chromosome="chr1", verbose=False)["value"])
```

- **示例输出：**

```text
elements 4
bases 350
bases-uniq 350
has-duplicate False
has-nested False
list-chromosomes ['chr1', 'chr2']
chr1 elements 3
```

## 4) `dnakit.annotation.starchcat_archives()`

- **作用：** 把多个 Starch 归档合并为一个多重集并集归档，或重新压缩单个归档。
- **处理方法：** 运行 `starchcat [--bzip2|--gzip] [--note=...] [--omit-signature] <archives>`，把新归档写入 `output_path`；元素数用 `unstarch --elements` 读取。只有一个输入时，只是重新压缩或更新注释。
- **API：** `starchcat_archives(input_paths: 'Sequence[str]', output_path: 'str', *, compression: 'str' = 'bzip2', note: 'str | None' = None, omit_signature: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `starchcat`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_paths` 是非空的、已存在的已排序且无头部的 Starch 归档列表。`output_path` 必须尚不存在。`compression`、`note`、`omit_signature` 与 `starch_compress` 相同。
- **输出：** 通用结果字典，另有 `element_count` 和 `size_bytes`。`outputs` 含 `archive`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, starchcat_archives, unstarch_summary

with tempfile.TemporaryDirectory() as tmp:
    first, second, merged = (str(Path(tmp) / name) for name in ("ref.starch", "map.starch", "merged.starch"))
    starch_compress("ref.bed", first, verbose=False)
    starch_compress("map.bed", second, verbose=False)
    result = starchcat_archives([first, second], merged, verbose=False)
    print(result["element_count"], unstarch_summary(merged, statistic="list-chromosomes", verbose=False)["value"])
```

- **示例输出：**

```text
8 ['chr1', 'chr2']
```

## 5) `dnakit.annotation.starchstrip_chromosomes()`

- **作用：** 保留或去掉 Starch 归档中的某些染色体，不重新压缩。
- **处理方法：** 运行 `starchstrip --include <list> <archive>` 并把结果写入 `output_path`；归档中不存在的染色体被忽略。BEDOPS 2.4.42 的 `starchstrip --exclude` 在去掉第一条染色体时会写出 `unstarch` 无法读取的归档，因此 `mode="exclude"` 时函数先列出归档的染色体，再传入互补的 `--include` 列表，得到预期的归档。去掉全部染色体会被拒绝。
- **API：** `starchstrip_chromosomes(input_path: 'str', output_path: 'str', chromosomes: 'Sequence[str]', *, mode: 'str' = 'include', bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `starchstrip`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的 Starch 归档。`output_path` 必须尚不存在。`chromosomes` 是不含逗号的非空名称列表。`mode` 取 `include` 或 `exclude`。
- **输出：** 通用结果字典，另有 `mode`、`strategy`、`chromosomes_kept` 和 `element_count`。`outputs` 含 `archive`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, starchstrip_chromosomes, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    result = starchstrip_chromosomes(archive, str(Path(tmp) / "no_chr1.starch"), ["chr1"], mode="exclude", verbose=False)
    print(result["chromosomes_kept"], result["element_count"])
    rows = unstarch_to_bed(result["outputs"]["archive"], str(Path(tmp) / "no_chr1.bed"), verbose=False)
    print(rows["preview"][0].replace("\t", " "))
```

- **示例输出：**

```text
['chr2'] 1
chr2 10 60 r4
```

## 6) `dnakit.annotation.convert_to_bed()`

- **作用：** 把 BAM、SAM、GFF3、GTF、GVF、PSL、RepeatMasker、VCF 或 WIG 文件转换成排序后的 BED 或 Starch 归档。
- **处理方法：** 运行 `convert2bed --input=<format> --output=<bed|starch> [格式选项] < input`，输入文件接到标准输入，结果写入 `output_path`。`input_format` 取 `bam`、`sam`、`gff`（GFF3）、`gtf`、`gvf`、`psl`、`rmsk`（RepeatMasker .out）、`vcf` 或 `wig`；`output_format` 取 `bed` 或 `starch`。一个函数覆盖 `bam2bed`、`sam2bed`、`gff2bed`、`gtf2bed`、`gvf2bed`、`psl2bed`、`rmsk2bed`、`vcf2bed`、`wig2bed` 和对应的 `*2starch` 程序，也包括只是分发同一转换的集群封装 `bam2starch_sge`、`bam2starch_slurm`、`bam2starch_gnuParallel`。`format_options` 是对该输入格式有效的布尔开关（`keep-header`；sam/bam 的 `split`、`all-reads`、`reduced`、`split-with-deletions`；psl 的 `split`；vcf 的 `snvs`、`insertions`、`deletions`、`do-not-split`；wig 的 `zero-indexed`）；`attribute_key` 为 `gff` 和 `gtf` 选择复制到 BED 名称列的 GENCODE 属性。`sam` 和 `bam` 输入需要 PATH 中有 `samtools`，否则抛出 `BackendUnavailableError`。该函数覆盖工作清单中把 BAM、SAM、GFF3、GTF、GVF、PSL、RepeatMasker OUT、VCF 和 WIG 转成 Starch，以及把 SAM、GVF、PSL、RepeatMasker OUT、VCF 和 WIG 转成 BED 的各行。
- **API：** `convert_to_bed(input_path: 'str', output_path: 'str', *, input_format: 'str', output_format: 'str' = 'bed', format_options: 'Sequence[str] | None' = None, attribute_key: 'str | None' = None, sort_output: 'bool' = True, max_memory: 'str | None' = None, temp_dir: 'str | None' = None, starch_compression: 'str' = 'bzip2', starch_note: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `convert2bed`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的、格式为 `input_format` 的文件。`output_path` 必须尚不存在。`format_options` 是该格式允许的选项名列表；`attribute_key` 是 `convert2bed --help-gff` / `--help-gtf` 列出的 GENCODE 键之一，只用于 gff 和 gtf。`sort_output=False` 加 `--do-not-sort`（仅 BED）。`max_memory`（`8G`）和 `temp_dir` 对应 `--max-mem` 和 `--sort-tmpdir`。`starch_compression`（`bzip2` 或 `gzip`）和 `starch_note` 只用于 `output_format="starch"`。
- **输出：** 通用结果字典，另有 `input_format`、`output_format`；BED 输出还有 `line_count` 和 `preview`（`outputs` 含 `output`），Starch 输出还有 `element_count` 和 `size_bytes`（`outputs` 含 `archive`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import convert_to_bed

with tempfile.TemporaryDirectory() as tmp:
    for name, fmt in (("variants.vcf", "vcf"), ("genes.gtf", "gtf"), ("repeats.out", "rmsk")):
        result = convert_to_bed(name, str(Path(tmp) / f"{fmt}.bed"), input_format=fmt, verbose=False)
        print(fmt, "->", result["preview"][0].replace("\t", " ")[:60])
    snvs = convert_to_bed("variants.vcf", str(Path(tmp) / "v.starch"), input_format="vcf", output_format="starch", verbose=False)
    print("vcf -> starch elements:", snvs["element_count"])
    named = convert_to_bed(
        "genes.gtf", str(Path(tmp) / "named.bed"), input_format="gtf", attribute_key="transcript_id", sort_output=False, verbose=False
    )
    print([line.split("\t")[3] for line in named["preview"]])
```

- **示例输出：**

```text
vcf -> chr1 50 51 rs2 60 C T PASS .
gtf -> chr1 100 200 g1 . + demo exon . gene_id "g1"; transcript_id 
rmsk -> chr1 100 200 L1PA2 463 + 10.2 0.0 1.2 (900) LINE/L1 1 100 (1
vcf -> starch elements: 2
['t1', 't2']
```
