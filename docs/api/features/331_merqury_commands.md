# Merqury：基于 k-mer 的组装评估、trio hap-mer 与 k-mer 数据库工具

`dnakit.assembly` 中的这一组函数把 [Merqury](https://github.com/marbl/merqury)（基于 Meryl 的基因组组装 k-mer 评估脚本）封装成带类型的接口：逐序列 QV、k-mer 多重度轨道、假重复计数与轨道、Mash QV、k-mer 长度选择、10x reads 拆分、solid k-mer 过滤、k-mer 数据库集合运算、直方图导出、trio hap-mer 推导、hap-mer 谱与 blob 计数、Hamming 错误、reads 归属与排除、带 switch 错误的相位块、N*/NG* 统计、bed 转 bigWig 以及三个 R 绘图。整个组装的 QV 与带拷贝数谱的完整度（`merqury.sh`）已由 `dnakit.assembly.evaluate_merqury` 提供。每个函数在先把输入链接进新输出目录之后，以 `bash <MERQURY>/<脚本> <参数>` 运行一个脚本，因此不会写入输入所在目录，已存在的输出目录会被拒绝。结果字典含 `status`、`backend`、`version`（Meryl）、`command`、`return_code`、`outputs`、`log_tail`、`operation`、`merqury_home`、`provenance` 以及解析出的值。只有退出码为 0 且所有预期文件都存在时 `status` 才是 `completed`。

测试版本为 Merqury 1.4.1（Meryl 1.4.2）；二者均来自 NIH/NHGRI（Merqury 的许可条款见其仓库），DNAKit 把它们作为独立程序启动，不复制代码。用 `conda install -c bioconda merqury` 安装（会一并安装 Meryl、带 ggplot2 的 R、Java、bedtools 和 samtools），DNAKit 通过环境变量 `MERQURY`、PATH 上的 `merqury.sh` 或 `merqury_home` 参数找到脚本。脚本或程序缺失时，会在创建任何输出之前抛出带安装提示的 `BackendUnavailableError`。上游脚本的已知限制：`count_10x.sh` 通过管道向 `meryl ... -` 输入，Meryl 1.4.2 不接受，`split_ont_fa.sh` 写死了集群的暂存目录，因此二者未封装；个别结尾的 R 绘图在很小的输入上会失败。

## 示例数据

示例在准备步骤中自己生成很小的数据：随机 4000 bp 基因组（分成两条 contig）、带 3% 差异的第二个单倍型、100 bp 的 reads 集、两个亲本 reads 集和一个子代 reads 集，全部计数为 15-mer 的 Meryl 数据库。其他输出都写入新建的临时文件夹。

**示例代码：**

```python
import random
import shutil
import subprocess
from pathlib import Path

rng = random.Random(5)
genome = "".join(rng.choice("ACGT") for _ in range(4000))
other = "".join(base if rng.random() > 0.03 else rng.choice("ACGT") for base in genome)


def write_reads(path, sequence, count, error=0.0):
    records = []
    for index in range(count):
        start = rng.randint(0, len(sequence) - 100)
        read = list(sequence[start:start + 100])
        if rng.random() < error:
            read[rng.randint(0, 99)] = "A"
        records.append(f">r{index}\n{''.join(read)}\n")
    Path(path).write_text("".join(records))


Path("asm.fa").write_text(f">c1\n{genome[:2000]}\n>c2\n{genome[2000:]}\n")
Path("hap2.fa").write_text(f">d1\n{other[:2000]}\n>d2\n{other[2000:]}\n")
write_reads("reads.fa", genome, 600, error=0.1)
write_reads("mat.fa", genome, 1500)
write_reads("pat.fa", other, 1500)
write_reads("child.fa", genome, 700)
with open("child.fa", "a") as handle:
    handle.write(Path("child.fa").read_text().replace(">r", ">p")[:0])
write_reads("child2.fa", other, 700)
Path("child.fa").write_text(Path("child.fa").read_text() + Path("child2.fa").read_text().replace(">r", ">q"))
for name in ("asm", "reads", "mat", "pat", "child"):
    shutil.rmtree(f"{name}.meryl", ignore_errors=True)
    subprocess.run(["meryl", "k=15", "count", "output", f"{name}.meryl", f"{name}.fa"], check=True, capture_output=True)
    print(name, "counted")
```

**示例输出：**

```text
asm counted
reads counted
mat counted
pat counted
child counted
```

准备步骤：生成 FASTA 文件和 15-mer Meryl 数据库（`asm.meryl`、`reads.meryl`、`mat.meryl`、`pat.meryl`、`child.meryl`）。

## 1) `dnakit.assembly.merqury_assign_reads()`

- **作用：** 按 hap-mer 含量把 reads 分配给单倍型 1、单倍型 2 或未知（`trio/assign_reads.sh`）。
- **处理方法：** 输入一个 reads 文件和两个 hap-mer 数据库，输出每条 read 的 k-mer 数、两个单倍型的 k-mer 数及命中数和归属。结果含 `assignments` 与 `counts`。
- **API：** `merqury_assign_reads(read_path: 'str', hap1_name: 'str', hap2_name: 'str', hap1_hapmer_meryl: 'str', hap2_hapmer_meryl: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_assign_reads

work = Path(tempfile.mkdtemp())
r = merqury_assign_reads("child.fa", "mat", "pat", "mat.meryl", "pat.meryl", str(work / "assign"), verbose=False)
print(r["status"], sum(r["counts"].values()))
```

- **示例输出：**

```text
completed 1400
```

## 2) `dnakit.assembly.merqury_bed_to_bigwig()`

- **作用：** 把区间 BED 转换成 bigWig 覆盖度轨道（`util/bed_to_bigwig.sh`）。
- **处理方法：** 先用 `bedtools genomecov` 得到 bedGraph，按序列排序后用 `bedGraphToBigWig` 转换；`convert_pipes` 把序列名中的 `|` 换成 `_`。
- **API：** `merqury_bed_to_bigwig(bed_path: 'str', fasta_index: 'str', output_dir: 'str', *, convert_pipes: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_bed_to_bigwig

work = Path(tempfile.mkdtemp())
__import__('subprocess').run(["samtools", "faidx", "asm.fa"], check=True)
Path("regions.bed").write_text("c1\t10\t200\nc2\t5\t90\n")
r = merqury_bed_to_bigwig("regions.bed", "asm.fa.fai", str(work / "bw"), verbose=False)
print(r["status"], Path(r["bigwig_file"]).name)
```

- **示例输出：**

```text
completed regions.bigwig
```

## 3) `dnakit.assembly.merqury_best_k()`

- **作用：** 根据基因组大小和可容忍的碰撞率选择 k-mer 长度（`best_k.sh`）。
- **处理方法：** k = log4(基因组大小 × (1 − 率) / 率)，默认率 0.001 适用于 Illumina WGS。结果含 `k` 与 `k_rounded_up`。
- **API：** `merqury_best_k(genome_size: 'int', *, tolerable_collision_rate: 'float' = 0.001, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_best_k

work = Path(tempfile.mkdtemp())
r = merqury_best_k(3_000_000_000, verbose=False)
print(r["status"], r["k"], r["k_rounded_up"])
```

- **示例输出：**

```text
completed 20.7233 21
```

## 4) `dnakit.assembly.merqury_block_n_stats()`

- **作用：** 相位块与 contig 的 N*/NG* 统计及曲线（`trio/block_n_stats.sh`）。
- **处理方法：** 输入一或两个组装 FASTA 及各自的 `*.phased_block.bed`，可选基因组大小以得到 NG*；输出尺寸表和 R 绘图。
- **API：** `merqury_block_n_stats(assemblies: 'Sequence[str]', block_beds: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, genome_size: 'int | None' = None, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
该函数没有在示例数据上用真实脚本运行过，因此这里不提供可运行示例。

## 5) `dnakit.assembly.merqury_combine_kmer_dbs()`

- **作用：** 合并 Meryl 数据库：差集、交集或求和并集（`build/diff.sh`、`intersect.sh`、`union_sum.sh`）。
- **处理方法：** `operation` 选择脚本；差集可对结果做 solid k-mer 过滤；并集需要 `k`。结果含 `database`。
- **API：** `merqury_combine_kmer_dbs(operation: 'str', meryl_dbs: 'Sequence[str]', output_dir: 'str', output_name: 'str', *, k: 'int | None' = None, filter_solid: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_combine_kmer_dbs

work = Path(tempfile.mkdtemp())
r = merqury_combine_kmer_dbs("intersect", ["asm.meryl", "reads.meryl"], str(work / "both"), "shared", verbose=False)
print(r["status"], Path(r["database"]).name)
u = merqury_combine_kmer_dbs("union_sum", ["asm.meryl", "reads.meryl"], str(work / "union"), "all", k=15, verbose=False)
print(u["status"], Path(u["database"]).name)
```

- **示例输出：**

```text
completed asm_and_reads.meryl
completed all.k15.meryl
```

## 6) `dnakit.assembly.merqury_exclude_reads()`

- **作用：** 删除含有某个 Meryl 数据库 k-mer 的 reads 及其配对 read（`trio/exclude_reads.sh`）。
- **处理方法：** `read_pairs` 为 `[R1]` 或 `[R1, R2]` 列表；`tenx` 忽略 R1 前 23 个碱基。结果含 `filtered_reads`。
- **API：** `merqury_exclude_reads(meryl_db: 'str', read_pairs: 'Sequence[Sequence[str]]', output_dir: 'str', *, tenx: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_exclude_reads

work = Path(tempfile.mkdtemp())
r = merqury_exclude_reads("mat.meryl", [["child.fa"]], str(work / "excl"), verbose=False)
print(r["status"], len(r["filtered_reads"]))
```

- **示例输出：**

```text
completed 1
```

## 7) `dnakit.assembly.merqury_false_duplications()`

- **作用：** 统计组装中位于单拷贝和双拷贝 reads 峰内的额外 k-mer（`eval/false_duplications.sh`）。
- **处理方法：** 输入 `*.spectra-cn.hist`，截断值为单拷贝峰 read 多重度的 1.5 倍；不写任何文件。
- **API：** `merqury_false_duplications(spectra_cn_hist: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_false_duplications

work = Path(tempfile.mkdtemp())
Path("spectra.hist").write_text("Copies\tkmer_multiplicity\tCount\n1\t1\t50\n1\t2\t10\n1\t5\t400\n2\t5\t30\n")
r = merqury_false_duplications("spectra.hist", verbose=False)
print(r["status"], r["cutoff"], r["copies"]["1"], r["duplicated_percent"])
```

- **示例输出：**

```text
completed 8 460 6.12245
```

## 8) `dnakit.assembly.merqury_false_duplications_track()`

- **作用：** 把假重复 k-mer 在组装上的位置输出为 BED 与 WIG 轨道（`eval/false_duplications_track.sh`）。
- **处理方法：** 需要组装的 Meryl 数据库（与 FASTA 同名的 `.meryl`）和 reads 数据库；结果含 `bed_file`、`wig_file`。
- **API：** `merqury_false_duplications_track(spectra_cn_hist: 'str', assembly_fasta: 'str', assembly_meryl: 'str', reads_meryl: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_false_duplications_track

work = Path(tempfile.mkdtemp())
r = merqury_false_duplications_track("spectra.hist", "asm.fa", "asm.meryl", "reads.meryl", str(work / "fd"), verbose=False)
print(r["status"], r["interval_count"])
```

- **示例输出：**

```text
completed 0
```

## 9) `dnakit.assembly.merqury_filter_solid_kmers()`

- **作用：** 用由倍性推得的计数阈值丢弃 reads 数据库中的错误 k-mer（`build/filt.sh`）。
- **处理方法：** 写出直方图、倍性/阈值文件，并用 `meryl greater-than` 保留高于阈值的 k-mer。结果含 `threshold` 与 `solid_database`。
- **API：** `merqury_filter_solid_kmers(meryl_db: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_filter_solid_kmers

work = Path(tempfile.mkdtemp())
r = merqury_filter_solid_kmers("reads.meryl", str(work / "filt"), verbose=False)
print(r["status"], r["threshold"], Path(r["solid_database"]).name)
```

- **示例输出：**

```text
completed 2 reads.gt2.meryl
```

## 10) `dnakit.assembly.merqury_hamming_error()`

- **作用：** 由每个 contig 的 hap-mer 计数计算组装的 Hamming（相位）错误（`trio/hamming_error.sh`）。
- **处理方法：** 输入 `merqury_hapmer_counts` 的计数表和一或两个组装 FASTA。
- **API：** `merqury_hamming_error(hapmer_count_table: 'str', assemblies: 'Sequence[str]', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hamming_error

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
from dnakit.assembly import merqury_hapmer_counts

hc = merqury_hapmer_counts(
    hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], ["asm.fa", "hap2.fa"], str(work / "counts"), "demo", verbose=False
)
r = merqury_hamming_error(hc["count_file"], ["asm.fa", "hap2.fa"], str(work / "ham"), verbose=False)
print(r["status"], r["columns"])
for row in r["rows"]:
    print(row)
```

- **示例输出：**

```text
completed ['Assembly', 'Hap1', 'Hap2', 'Hap1Error(%)', 'Hap2Error(%)', 'MajorHapError(%)']
['asm', '1025', '0', '0', '100', '0']
['hap2', '0', '1009', '100', '0', '0']
```

## 11) `dnakit.assembly.merqury_hapmer_counts()`

- **作用：** 统计每个 contig 中各单倍型 hap-mer 的数量并绘制 blob 图（`trio/hap_blob.sh`）。
- **处理方法：** 输入两个 hap-mer 数据库和一或两个组装 FASTA，输出 `*.hapmers.count` 与 blob 图（R）。
- **API：** `merqury_hapmer_counts(hap1_meryl: 'str', hap2_meryl: 'str', assemblies: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hapmer_counts

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
hc = merqury_hapmer_counts(
    hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], ["asm.fa", "hap2.fa"], str(work / "counts"), "demo", verbose=False
)
print(hc["status"], hc["columns"])
for row in hc["rows"]:
    print(row)
```

- **示例输出：**

```text
completed ['Assembly', 'Contig', 'mat.inherited.gt1', 'pat.inherited.gt6', 'Size']
['asm', 'c1', '482', '0', '1986']
['asm', 'c2', '543', '0', '1986']
['hap2', 'd1', '0', '472', '1986']
['hap2', 'd2', '0', '537', '1986']
```

## 12) `dnakit.assembly.merqury_hapmer_spectra()`

- **作用：** 一或两个组装的 hap-mer 拷贝数谱与完整度（`trio/spectra-hap.sh`）。
- **处理方法：** 输入 reads 数据库、两个单倍型特异 k-mer 数据库（含 reads 计数）和组装 FASTA，输出完整度表和谱图（R）。
- **API：** `merqury_hapmer_spectra(reads_meryl: 'str', hap1_meryl: 'str', hap2_meryl: 'str', assemblies: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, homopolymer_compressed: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
该函数没有在示例数据上用真实脚本运行过，因此这里不提供可运行示例。

## 13) `dnakit.assembly.merqury_hapmers()`

- **作用：** 从亲本 k-mer 数据库推导单倍型特异 k-mer（hap-mer），并给出子代继承的 hap-mer（`trio/hapmers.sh`）。
- **处理方法：** 输出 `<hap>.only.meryl`、`shrd.meryl` 以及（给出子代时）`<hap>.hapmer.meryl`、`cutoffs.txt` 与谱图。结果含 `hapmer_databases` 与 `cutoffs`。
- **API：** `merqury_hapmers(hap1_meryl: 'str', hap2_meryl: 'str', output_dir: 'str', *, child_meryl: 'str | None' = None, no_filt: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hapmers

work = Path(tempfile.mkdtemp())
hapmers_result = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
print(hapmers_result["status"], hapmers_result["cutoffs"])
hapmers = {Path(p).name: p for p in hapmers_result["hapmer_databases"]}
print(sorted(name for name in hapmers if name.endswith(".hapmer.meryl")))
```

- **示例输出：**

```text
completed {'mat': 1, 'pat': 6, 'shared': 2}
['mat.hapmer.meryl', 'pat.hapmer.meryl']
```

## 14) `dnakit.assembly.merqury_kmer_histogram_table()`

- **作用：** 导出 Meryl 数据库的 k-mer 多重度直方图，格式即 Merqury 绘图脚本读取的表（`plot/to_hist_for_plotting.sh`）。
- **处理方法：** `meryl_dbs` 把标签映射到数据库；输出三列：`kmer`、`kmer_multiplicity`、`Count`。
- **API：** `merqury_kmer_histogram_table(meryl_dbs: 'Mapping[str, str]', output_path: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_kmer_histogram_table

work = Path(tempfile.mkdtemp())
r = merqury_kmer_histogram_table({"reads": "reads.meryl"}, str(work / "hist.tsv"), verbose=False)
print(r["status"], r["columns"], r["rows"][:2])
```

- **示例输出：**

```text
completed ['kmer', 'kmer_multiplicity', 'Count'] [['reads', '1', '487'], ['reads', '2', '4']]
```

## 15) `dnakit.assembly.merqury_kmer_multiplicity_track()`

- **作用：** 沿组装序列输出 k-mer 计数轨道（WIG 与 bigWig）（`eval/asm_multiplicity.sh`、`eval/read_multiplicity.sh`）。
- **处理方法：** `source='assembly'` 给出组装 k-mer 在组装中的拷贝数；`source='reads'` 给出其在 reads 中的多重度。需要 samtools 与 wigToBigWig。
- **API：** `merqury_kmer_multiplicity_track(assembly_fasta: 'str', meryl_db: 'str', output_dir: 'str', *, source: 'str' = 'assembly', prefix: 'str' = 'track', merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_kmer_multiplicity_track

work = Path(tempfile.mkdtemp())
for source, db in (("assembly", "asm.meryl"), ("reads", "reads.meryl")):
    r = merqury_kmer_multiplicity_track("asm.fa", db, str(work / source), source=source, verbose=False)
    print(source, r["status"], Path(r["wig_file"]).name, Path(r["bigwig_file"]).name)
```

- **示例输出：**

```text
assembly completed track.copies.wig track.asm_multiplicity.bigWig
reads completed track.read.wig track.read_multiplicity.bigwig
```

## 16) `dnakit.assembly.merqury_mash_qv()`

- **作用：** 用 Mash 草图快速估计组装 QV（`merqury-mash.sh`）。
- **处理方法：** 对组装做 Mash sketch，再用 `mash screen` 对 reads 打分并转换为 QV；结果与完整 Merqury 不完全一致。
- **API：** `merqury_mash_qv(assembly_fasta: 'str', read_paths: 'Sequence[str]', output_dir: 'str', *, k: 'int' = 21, cpus: 'int' = 1, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_mash_qv

work = Path(tempfile.mkdtemp())
r = merqury_mash_qv("asm.fa", ["reads.fa"], str(work / "mash"), k=15, verbose=False)
print(r["status"], r["columns"], r["rows"][0][:2])
```

- **示例输出：**

```text
completed ['assembly', 'reads', 'qv', 'error_rate'] ['asm', '3965']
```

## 17) `dnakit.assembly.merqury_per_sequence_qv()`

- **作用：** 由 reads 的 k-mer 数据库估计组装中每条序列的 QV（`eval/per_seq_qv.sh`）。
- **处理方法：** 统计只出现在组装中的 k-mer，QV = −10·log10(错误率)；输出 `<prefix>.qv`。
- **API：** `merqury_per_sequence_qv(assembly_fasta: 'str', read_meryl: 'str', output_dir: 'str', *, prefix: 'str' = 'per_seq', merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_per_sequence_qv

work = Path(tempfile.mkdtemp())
r = merqury_per_sequence_qv("asm.fa", "reads.meryl", str(work / "qv"), verbose=False)
print(r["status"], r["columns"])
for row in r["rows"]:
    print(row)
```

- **示例输出：**

```text
completed ['sequence', 'error_kmers', 'total_kmers', 'qv', 'error_rate']
['c1', '0', '1986', '+inf', '0']
['c2', '7', '1986', '36.2826', '0.000235366']
```

## 18) `dnakit.assembly.merqury_phase_blocks()`

- **作用：** 沿组装用 hap-mer 标记构建相位块并统计 switch 错误（`trio/phase_block.sh`、`trio/switch_error.sh`）。
- **处理方法：** 先运行 `phase_block.sh`（固定 100 次切换/20000 bp），参数不同时再对 `<prefix>.sort.bed` 运行 `switch_error.sh`。需要 Java。
- **API：** `merqury_phase_blocks(assembly_fasta: 'str', hap1_meryl: 'str', hap2_meryl: 'str', output_dir: 'str', output_prefix: 'str', *, num_switch: 'int' = 100, short_range: 'int' = 20000, include_gaps: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_phase_blocks

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
r = merqury_phase_blocks(
    "asm.fa", hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], str(work / "phase"), "demo",
    num_switch=50, short_range=10000, verbose=False,
)
print(r["status"], [Path(p).name for p in r["block_files"]])
```

- **示例输出：**

```text
completed ['demo.100_20000.phased_block.bed', 'demo.50_10000.phased_block.bed']
```

## 19) `dnakit.assembly.merqury_plot_blob()`

- **作用：** 把每个 contig 的 hap-mer 计数画成 blob 图（`plot/plot_blob.R`）。
- **处理方法：** 输入含表头的 5 列制表符表（类别、序列 ID、hap1 计数、hap2 计数、序列长度）。需要 Rscript。
- **API：** `merqury_plot_blob(count_table: 'str', output_dir: 'str', output_prefix: 'str', *, width: 'float' = 6.5, height: 'float' = 6.0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_blob

work = Path(tempfile.mkdtemp())
Path("c.count").write_text("category\tseq\thap1\thap2\tsize\nhap1\tc1\t10\t3\t2000\nhap1\tc2\t4\t5\t1000\n")
r = merqury_plot_blob("c.count", str(work / "blob"), "blob", verbose=False)
print(r["status"], Path(r["plot"]).name)
```

- **示例输出：**

```text
completed blob.png
```

## 20) `dnakit.assembly.merqury_plot_block_n()`

- **作用：** 绘制相位块、scaffold、contig 的 N*/NG* 连续性曲线（`plot/plot_block_N.R`）。
- **处理方法：** 输入三列（类型、分组、大小）的 `.sizes` 表，可选基因组大小以得到 NG*。需要 Rscript。
- **API：** `merqury_plot_block_n(output_dir: 'str', output_prefix: 'str', *, block_sizes: 'str | None' = None, scaffold_sizes: 'str | None' = None, contig_sizes: 'str | None' = None, genome_size: 'int' = 0, width: 'float' = 6.0, height: 'float' = 5.0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_block_n

work = Path(tempfile.mkdtemp())
Path("c.sizes").write_text("contig\tasm\t2000\ncontig\tasm\t1000\n")
r = merqury_plot_block_n(str(work / "bn"), "bn", contig_sizes="c.sizes", verbose=False)
print(r["status"], [Path(p).name for p in r["plots"]])
```

- **示例输出：**

```text
completed ['bn.block.N.png', 'bn.continuity.N.png']
```

## 21) `dnakit.assembly.merqury_plot_spectra_cn()`

- **作用：** 由 `*.spectra-cn.hist` 表绘制 k-mer 拷贝数谱（`plot/plot_spectra_cn.R`）。
- **处理方法：** `plot_type` 为 line、fill、stack 或 all，对应 `ln`/`fl`/`st` 三个图文件。需要 Rscript。
- **API：** `merqury_plot_spectra_cn(histogram: 'str', output_dir: 'str', output_prefix: 'str', *, only_histogram: 'str | None' = None, cutoff_file: 'str | None' = None, plot_type: 'str' = 'all', width: 'float' = 6.0, height: 'float' = 5.0, xmax: 'int' = 0, ymax: 'int' = 0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_spectra_cn

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_kmer_histogram_table

merqury_kmer_histogram_table({"reads": "reads.meryl"}, str(work / "h.tsv"), verbose=False)
r = merqury_plot_spectra_cn(str(work / "h.tsv"), str(work / "plots"), "spectra", plot_type="line", verbose=False)
print(r["status"], [Path(p).name for p in r["plots"]])
```

- **示例输出：**

```text
completed ['spectra.ln.png']
```

## 22) `dnakit.assembly.merqury_split_10x_fastq()`

- **作用：** 去掉 10x R1 reads 前 23 个碱基并按 3 亿行拆分、压缩（`build/split_10x.sh`）。
- **处理方法：** 输出到 `split/<name>.NNNN.fq.gz`，并写出每个输入的清单；需要 pigz。结果含 `pieces`。
- **API：** `merqury_split_10x_fastq(read_paths: 'Sequence[str]', output_dir: 'str', *, threads: 'int' = 8, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 Merqury 及其脚本所用的程序（`conda install -c bioconda merqury`；Meryl，视脚本而定还有 Java、R、samtools、bedtools、Mash、pigz）。DNAKit 通过 `merqury_home`、`MERQURY` 或 PATH 上的 `merqury.sh` 找到脚本；缺少任何一项都会在创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_split_10x_fastq

work = Path(tempfile.mkdtemp())
r = merqury_split_10x_fastq(["reads.fa"], str(work / "split"), verbose=False)
print(r["status"], [v for v in r["pieces"].values()])
```

- **示例输出：**

```text
completed [['split/reads.fa.0000.fq.gz']]
```
