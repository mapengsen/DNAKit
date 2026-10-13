# 可选 tskit 与 dadi 群体遗传函数（第 80 组）

这一组函数包装 tskit 与 dadi 中与 DNA 或基因相关、可在本地离线运行的公开函数。tskit 部分在 `dnakit.population`，负责序列与基因组区间导出、文本表读写、PLINK `.fam` 读取、拓扑计数、多叉树拆分与 Newick 输出；dadi 部分在 `dnakit.population`，负责 SNP 与 VCF 数据读取、频率谱的读写与合并、自助抽样、人口模型的期望频率谱、模型拟合与 Fisher 信息不确定性。

这些函数从 `dnakit.population` 导出，示例从该包导入。依赖都是可选的，只有调用对应函数时才会导入，缺少包时报错并给出安装命令。输出文件从不覆盖已有文件；模型只能从固定名单中按名称选择，不执行任意代码。

## 示例数据

第 1 至 12 项为 tskit 示例，第 13 至 30 项为 dadi 示例。下面的文件在运行示例前会写入工作目录。`tables.json` 有 4 个样本、1 棵树和 3 个位点，只用于演示接口，统计值不代表真实群体。

**示例文件：** `tables.json`

```json
{
 "sequence_length": 100,
 "nodes": [
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 1
  },
  {
   "time": 2
  },
  {
   "time": 3
  }
 ],
 "edges": [
  {
   "left": 0,
   "right": 100,
   "parent": 4,
   "child": 0
  },
  {
   "left": 0,
   "right": 100,
   "parent": 4,
   "child": 1
  },
  {
   "left": 0,
   "right": 100,
   "parent": 5,
   "child": 2
  },
  {
   "left": 0,
   "right": 100,
   "parent": 5,
   "child": 3
  },
  {
   "left": 0,
   "right": 100,
   "parent": 6,
   "child": 4
  },
  {
   "left": 0,
   "right": 100,
   "parent": 6,
   "child": 5
  }
 ],
 "sites": [
  {
   "position": 10,
   "ancestral_state": "A"
  },
  {
   "position": 20,
   "ancestral_state": "A"
  },
  {
   "position": 30,
   "ancestral_state": "A"
  }
 ],
 "mutations": [
  {
   "site": 0,
   "node": 0,
   "derived_state": "G"
  },
  {
   "site": 1,
   "node": 4,
   "derived_state": "G"
  },
  {
   "site": 2,
   "node": 5,
   "derived_state": "G"
  }
 ]
}
```

**示例文件：** `nodes.txt`

```text
id	is_sample	time	population	individual	metadata
0	1	0.000000	-1	-1	
1	1	0.000000	-1	-1	
2	1	0.000000	-1	-1	
3	1	0.000000	-1	-1	
4	0	1.000000	-1	-1	
5	0	2.000000	-1	-1	
6	0	3.000000	-1	-1	
```

**示例文件：** `edges.txt`

```text
left	right	parent	child	metadata
0.000000	100.000000	4	0	
0.000000	100.000000	4	1	
0.000000	100.000000	5	2	
0.000000	100.000000	5	3	
0.000000	100.000000	6	4	
0.000000	100.000000	6	5	
```

**示例文件：** `sites.txt`

```text
position	ancestral_state	metadata
10.000000	A	
20.000000	A	
30.000000	A	
```

**示例文件：** `mutations.txt`

```text
site	node	time	derived_state	parent	metadata
0	0	unknown	G	-1	
1	4	unknown	G	-1	
2	5	unknown	G	-1	
```

**示例文件：** `samples.fam`

```text
f1 i1 0 0 1
f1 i2 0 0 2
f2 i3 0 0 0
```

**示例文件：** `snps.txt`

```text
Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid
ACA TCA A 2 1 G 0 1 snp1
GCT GCT C 1 2 T 1 0 snp2
TTG TAG T 2 2 G 0 0 snp3
ACG ACG A 1 1 C 1 1 snp4
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
chr1	100	.	A	G	.	PASS	AA=A	GT	0|0	0|1	1|1	0|0
chr1	200	.	C	T	.	PASS	AA=T	GT	0|1	1|1	0|0	0|0
chr1	300	.	G	A	.	PASS	AA=G	GT	0|0	0|0	0|1	1|1
```

**示例文件：** `popinfo.txt`

```text
s1 popA
s2 popA
s3 popB
s4 popB
```

**示例文件：** `chunks.txt`

```text
Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid
ACG ACG A 2 1 G 0 1 chr1_100
ACG ACG A 1 0 G 1 2 chr1_200
ACG ACG A 0 2 G 2 0 chr1_1100
ACG ACG A 1 1 G 1 1 chr1_1200
ACG ACG A 2 0 G 0 2 chr1_2100
ACG ACG A 0 1 G 2 1 chr1_2200
```

**示例文件：** `sample.exonic_variant_function`

```text
1	nonsynonymous SNV	GENE1	chr1 100
2	synonymous SNV	GENE2	chr1 200
```

**示例文件：** `spectrum.sfs`

```text
# demo
2 2 
1.5 2 3 4
```

**示例文件：** `model.yaml`

```text
description: split with migration
time_units: generations
demes:
  - name: ancestral
    epochs:
      - end_time: 2000
        start_size: 1000
  - name: A
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
  - name: B
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
migrations:
  - source: A
    dest: B
    start_time: 2000
    end_time: 0
    rate: 0.0001
```

## 1) `dnakit.population.export_sample_sequences()`

- **作用：** 把每个样本的单倍型序列写成 FASTA 或 NEXUS 文件。
- **处理方法：** 调用 tskit 的 `write_fasta()` 或 `write_nexus()`。序列由位点和突变推导，没有位点的位置用 `missing_data_character` 填充。输出文件已存在时报错。
- **API：** `dnakit.population.export_sample_sequences(tables, destination, *, file_format='fasta', wrap_width=60, missing_data_character='N', max_rows=100000)`
- **输入：** `tables` 与 `build_tree_sequence()` 相同；`destination` 是新文件路径；`file_format` 为 `fasta`（默认）或 `nexus`；`wrap_width` 是 FASTA 的行宽，0 表示不换行；`missing_data_character` 是单个字符。
- **输出：** 含 `file_format`、`destination`、`bytes`、`num_samples`、`sequence_length`、`num_sites` 和 `provenance` 的字典。`file_format="nexus"` 时同样写出 NEXUS 文件。
- **示例代码：**

```python
import json
from pathlib import Path
from dnakit.population import export_sample_sequences

Path("samples.fa").unlink(missing_ok=True)  # the function never overwrites an existing file
result = export_sample_sequences(json.load(open("tables.json")), "samples.fa", wrap_width=60)
lines = Path("samples.fa").read_text().splitlines()
print(result["num_samples"], lines[0], lines[1][10], lines[1][20], lines[1][30])
```

- **示例输出：**

```text
4 >n0 G G A
```

## 2) `dnakit.population.export_ms_format()`

- **作用：** 把树序列的基因型写成 Hudson `ms` 格式。
- **处理方法：** 调用 tskit 的 `write_ms()`。每行是一个样本，每列是一个分离位点，1 表示衍生等位基因。`print_trees=True` 时同时写出 Newick 树。
- **API：** `dnakit.population.export_ms_format(tables, destination, *, print_trees=False, precision=4, max_rows=100000)`
- **输入：** `tables` 与 `build_tree_sequence()` 相同；`destination` 是新文件路径；`precision` 是位点位置的小数位数。
- **输出：** 含 `destination`、`bytes`、`num_samples`、`num_segregating_sites` 和 `provenance` 的字典。
- **示例代码：**

```python
import json
from pathlib import Path
from dnakit.population import export_ms_format

Path("samples.ms").unlink(missing_ok=True)
result = export_ms_format(json.load(open("tables.json")), "samples.ms", precision=4)
print(result["num_segregating_sites"], Path("samples.ms").read_text().splitlines()[6:10])
```

- **示例输出：**

```text
3 ['110', '010', '001', '001']
```

## 3) `dnakit.population.complement_intervals()`

- **作用：** 返回 `[start, end)` 中没有被给定区间覆盖的半开区间。
- **处理方法：** 调用 tskit 的 `negate_intervals()`。输入区间必须按位置排序、互不重叠，并且位于 `[start, end)` 之内。空列表返回整个区间。
- **API：** `dnakit.population.complement_intervals(intervals, *, start=0.0, end)`
- **输入：** `intervals` 是 `[left, right]` 的列表；`start` 默认 0；`end` 必须给出。
- **输出：** 含 `intervals`（补集区间）、`num_intervals`、`complement_length`、`covered_length` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import complement_intervals

result = complement_intervals([(10, 20), (30, 40)], end=100)
print(result["intervals"], result["complement_length"])
```

- **示例输出：**

```text
[[0.0, 10.0], [20.0, 30.0], [40.0, 100.0]] 80.0
```

## 4) `dnakit.population.read_text_table()`

- **作用：** 读取一个 tskit 文本表文件，返回其中的行。
- **处理方法：** 按 `table` 选择 tskit 的解析器（`parse_nodes`、`parse_edges`、`parse_sites`、`parse_mutations`、`parse_individuals`、`parse_populations`、`parse_migrations`）。元数据列不返回。`strict=True` 要求制表符分隔。
- **API：** `dnakit.population.read_text_table(source, *, table='nodes', strict=True, max_rows=100000)`
- **输入：** `source` 是文本表文件路径；`table` 是表的类型；`max_rows` 限制行数。
- **输出：** 含 `table`、`num_rows`、`columns` 和 `rows`（每行一个字典）的字典。
- **示例代码：**

```python
from dnakit.population import read_text_table

result = read_text_table("sites.txt", table="sites")
print(result["num_rows"], [row["position"] for row in result["rows"]])
```

- **示例输出：**

```text
3 [10.0, 20.0, 30.0]
```

## 5) `dnakit.population.read_text_tree_sequence()`

- **作用：** 由 tskit 文本表文件构建树序列并汇总。
- **处理方法：** 调用 `tskit.load_text()`。`nodes` 和 `edges` 必填，`sites` 与 `mutations` 可选。`sequence_length=0` 时由边推断序列长度。
- **API：** `dnakit.population.read_text_tree_sequence(nodes, edges, *, sites=None, mutations=None, sequence_length=0.0, strict=True)`
- **输入：** `nodes`、`edges` 是必需的文本文件；`sites`、`mutations` 可选；`strict` 与 `read_text_table()` 相同。
- **输出：** 含 `sources`、`num_samples`、`num_trees`、`num_mutations` 等汇总字段和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import read_text_tree_sequence

result = read_text_tree_sequence("nodes.txt", "edges.txt", sites="sites.txt", mutations="mutations.txt")
print(result["num_samples"], result["num_trees"], result["num_mutations"])
```

- **示例输出：**

```text
4 1 3
```

## 6) `dnakit.population.write_text_tree_sequence()`

- **作用：** 把树序列写成 tskit 文本表（节点、边、位点、突变）。
- **处理方法：** 调用 `tskit.dump_text()`，把四张表写入新建的目录。目录已存在时报错。输出可由 `read_text_tree_sequence()` 读回。
- **API：** `dnakit.population.write_text_tree_sequence(tables, destination_dir, *, precision=6, max_rows=100000)`
- **输入：** `tables` 与 `build_tree_sequence()` 相同；`destination_dir` 是不存在的目录；`precision` 是时间的小数位数。
- **输出：** 含 `destination_dir`、`files`、`num_nodes`、`num_edges`、`num_sites`、`num_mutations` 和 `provenance` 的字典。
- **示例代码：**

```python
import json
import shutil
from dnakit.population import write_text_tree_sequence

shutil.rmtree("text_out", ignore_errors=True)  # the directory must not exist yet
result = write_text_tree_sequence(json.load(open("tables.json")), "text_out")
print(result["files"], result["num_mutations"])
```

- **示例输出：**

```text
['nodes.txt', 'edges.txt', 'sites.txt', 'mutations.txt'] 3
```

## 7) `dnakit.population.count_sample_topologies()`

- **作用：** 统计每棵树中命名样本集合的嵌入拓扑计数。
- **处理方法：** 调用 tskit 的 `treeseq_count_topologies()`，对每个样本集合组合返回计数，并按拓扑秩（`shape`、`label`）排序。标签是集合名用 `+` 连接。每个集合至少两个样本。
- **API：** `dnakit.population.count_sample_topologies(tables, sample_sets, *, max_rows=100000)`
- **输入：** `tables` 与 `build_tree_sequence()` 相同；`sample_sets` 是 `{名称: [样本节点]}`。
- **输出：** 含 `sample_sets`、`num_trees` 和 `trees`（每棵树的计数）的字典。
- **示例代码：**

```python
import json
from dnakit.population import count_sample_topologies

result = count_sample_topologies(json.load(open("tables.json")), {"A": [0, 1], "B": [2, 3]})
print(result["num_trees"], sorted(result["trees"][0]["counts"]))
```

- **示例输出：**

```text
1 ['A', 'A+B', 'B']
```

## 8) `dnakit.population.resolve_polytomies()`

- **作用：** 把一棵树中的多叉节点拆成二叉节点，并返回 Newick 字符串。
- **处理方法：** 调用 tskit 的 `split_polytomies()`。`random_seed` 使随机拆分可以复现。已经是二叉的树不变，`added_nodes` 为 0。
- **API：** `dnakit.population.resolve_polytomies(tables, *, tree_index=0, random_seed=None, max_rows=100000)`
- **输入：** `tables`；`tree_index` 选择树；`random_seed` 可选。
- **输出：** 含 `newick`、`added_nodes`、`added_edges` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import resolve_polytomies

star = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}],
    "edges": [{"left": 0, "right": 100, "parent": 4, "child": child} for child in range(4)],
    "sites": [],
    "mutations": [],
}
result = resolve_polytomies(star, random_seed=7)
print(result["added_nodes"], result["newick"].count("("))
```

- **示例输出：**

```text
2 3
```

## 9) `dnakit.population.tree_to_newick()`

- **作用：** 把树序列中的一棵树输出为 Newick 字符串。
- **处理方法：** 调用 tskit 的 `build_newick()`，样本节点用 `taxon_labels` 命名（默认 `n<节点编号>`）。标签必须唯一，且不含 Newick 保留字符。
- **API：** `dnakit.population.tree_to_newick(tables, *, tree_index=0, precision=6, include_branch_lengths=True, taxon_labels=None, max_rows=100000)`
- **输入：** `tables`；`tree_index`；`precision` 是支长的小数位数；`include_branch_lengths` 控制是否写支长。
- **输出：** 含 `tree_index`、`num_samples`、`newick` 和 `provenance` 的字典。
- **示例代码：**

```python
import json
from dnakit.population import tree_to_newick

result = tree_to_newick(json.load(open("tables.json")), taxon_labels=["a", "b", "c", "d"], precision=3)
print(result["newick"])
```

- **示例输出：**

```text
((a:1.000,b:1.000):2.000,(c:2.000,d:2.000):1.000);
```

## 10) `dnakit.population.random_tree_newick()`

- **作用：** 生成一棵给定叶子数的随机或规则树，并返回 Newick 字符串。
- **处理方法：** 调用 tskit 的组合生成器。`random_binary` 在叶子标号的二叉树中均匀抽样，`seed` 使其可复现；`balanced`、`star`、`comb` 为确定的形状。叶子命名为 `n0` 至 `n<叶子数-1>`。
- **API：** `dnakit.population.random_tree_newick(num_leaves, *, shape='random_binary', seed=None, arity=2, branch_length=1.0)`
- **输入：** `num_leaves` 为 2 至 1000；`shape` 为四种形状之一；`seed` 仅用于 `random_binary`；`arity` 仅用于 `balanced`。
- **输出：** 含 `shape`、`num_leaves`、`seed`、`newick` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import random_tree_newick

result = random_tree_newick(5, seed=3)
print(result["newick"])
```

- **示例输出：**

```text
(n0:4.000000,(n3:3.000000,(n4:2.000000,(n1:1.000000,n2:1.000000):1.000000):1.000000):1.000000);
```

## 11) `dnakit.population.map_allele_indices()`

- **作用：** 把一组等位基因映射到另一组等位基因的索引。
- **处理方法：** 调用 tskit 的 `allele_remap()`，每个等位基因取它在 `alleles_to` 中第一次出现的位置。不在 `alleles_to` 中的等位基因获得 `alleles_to` 之后的新索引，并列在 `unmapped_alleles` 中。
- **API：** `dnakit.population.map_allele_indices(alleles_from, alleles_to)`
- **输入：** `alleles_from` 和 `alleles_to` 是等位基因字符串列表；`alleles_to` 不得重复。
- **输出：** 含 `indices`、`num_alleles_to`、`unmapped_alleles` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import map_allele_indices

result = map_allele_indices(["A", "G", "T", "X"], ["A", "G"])
print(result["indices"], result["unmapped_alleles"])
```

- **示例输出：**

```text
[0, 1, 2, 3] ['T', 'X']
```

## 12) `dnakit.population.read_plink_fam()`

- **作用：** 读取 PLINK `.fam` 文件，返回个体列表。
- **处理方法：** 调用 tskit 的 `parse_fam()`。文件有五列：家系号、个体号、父本、母本、性别（1 为男，2 为女，其他为未知）。个体按文件顺序返回。
- **API：** `dnakit.population.read_plink_fam(path)`
- **输入：** `path` 是 `.fam` 文件路径。
- **输出：** 含 `num_individuals`、`num_males`、`num_females`、`num_unknown_sex`、`individuals` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import read_plink_fam

result = read_plink_fam("samples.fam")
print(result["num_individuals"], result["num_males"], result["num_females"], result["num_unknown_sex"])
```

- **示例输出：**

```text
3 1 1 1
```

## 13) `dnakit.population.read_snp_data()`

- **作用：** 读取 dadi 的 SNP 文件并汇总。
- **处理方法：** 调用 dadi 的 `Misc.make_data_dict()`。文件头为 `Context OutgroupContext Allele1 <群体...> Allele2 <编号>`，每行是一个 SNP。外群等位基因属于两个分离等位基因之一时，该 SNP 被极化。
- **API：** `dnakit.population.read_snp_data(path)`
- **输入：** `path` 是 dadi SNP 文件（可为 `.gz` 或 `.zip`）。
- **输出：** 含 `populations`、`num_snps`、`num_biallelic`、`num_polarized` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import read_snp_data

result = read_snp_data("snps.txt")
print(result["num_snps"], result["num_polarized"], result["populations"])
```

- **示例输出：**

```text
4 2 ['pop1', 'pop2']
```

## 14) `dnakit.population.snp_configuration_counts()`

- **作用：** 按 SNP 构型统计指定群体的计数。
- **处理方法：** 调用 dadi 的 `Misc.count_data_dict()`。每行给出各群体的成功检出数与衍生检出数、是否极化，以及该构型的 SNP 数。按构型排序。非双等位基因位点被跳过。
- **API：** `dnakit.population.snp_configuration_counts(path, *, pop_ids)`
- **输入：** `path` 是 SNP 文件；`pop_ids` 是群体编号列表，必须都在文件中。
- **输出：** 含 `pop_ids`、`num_configurations`、`num_snps_counted`、`configurations` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import snp_configuration_counts

result = snp_configuration_counts("snps.txt", pop_ids=["pop1", "pop2"])
print(result["num_configurations"], result["num_snps_counted"], result["configurations"][0])
```

- **示例输出：**

```text
4 4 {'successful_calls': [2, 2], 'derived_calls': [0, 0], 'polarized': False, 'count': 1}
```

## 15) `dnakit.population.read_vcf_data()`

- **作用：** 读取 VCF 文件和群体对应表，得到 dadi 的数据字典并汇总。
- **处理方法：** `parser="builtin"` 用 dadi 自带的 VCF 读取器；`parser="cyvcf2"` 用 cyvcf2（需安装）。两者对同一文件给出相同的字典。`filter=True` 跳过未通过过滤的位点。
- **API：** `dnakit.population.read_vcf_data(vcf_path, popinfo_path, *, parser='builtin', filter=True)`
- **输入：** `vcf_path` 是 VCF 文件；`popinfo_path` 每行为 `样本 群体`。
- **输出：** 含 `parser`、`populations`、`num_snps` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import read_vcf_data

result = read_vcf_data("variants.vcf", "popinfo.txt")
print(result["parser"], result["num_snps"], result["populations"])
```

- **示例输出：**

```text
builtin 3 ['popA', 'popB']
```

## 16) `dnakit.population.bootstrap_snp_spectra()`

- **作用：** 从 SNP 文件的染色体片段中自助抽样，得到频率谱的均值与标准差。
- **处理方法：** 调用 dadi 的 `Misc.fragment_data_dict()` 与 `Misc.bootstraps_from_dd_chunks()`。SNP 编号须以碱基位置结尾（如 `chr1_1500`）。按 `chunk_size` 分段，每次有放回地抽取片段。`seed` 使抽样可复现，并在调用后恢复全局随机状态。
- **API：** `dnakit.population.bootstrap_snp_spectra(path, *, pop_ids, projections, chunk_size, num_boot, seed=0, polarized=True)`
- **输入：** `path` 是 SNP 文件；`pop_ids` 与 `projections`（每个群体的投影样本量）；`chunk_size` 为碱基数；`num_boot` 为重抽次数（最多 1000）。
- **输出：** 含 `num_chunks`、`num_boot`、`mean_spectrum`、`std_spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import bootstrap_snp_spectra

result = bootstrap_snp_spectra(
    "chunks.txt", pop_ids=["pop1", "pop2"], projections=[2, 2], chunk_size=1000, num_boot=20, seed=5
)
print(result["num_chunks"], result["num_boot"])
print([[None if v is None else round(v, 4) for v in row] for row in result["mean_spectrum"]])
```

- **示例输出：**

```text
3 20
[[None, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, None]]
```

## 17) `dnakit.population.pool_sfs_populations()`

- **作用：** 把二或三群体频率谱中的两个群体合并为一个轴。
- **处理方法：** 调用 dadi 的 `Misc.combine_pops()`。被选中的两个轴按 `i + j` 相加，合并后的轴成为结果的第一轴。三群体时只支持 (0, 1)、(0, 2)、(1, 2)。
- **API：** `dnakit.population.pool_sfs_populations(spectrum, *, populations=(0, 1))`
- **输入：** `spectrum` 是二维或三维的列表；`populations` 是两个从 0 开始的轴编号。
- **输出：** 含 `populations`、`sample_sizes`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import pool_sfs_populations

result = pool_sfs_populations([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
print(result["sample_sizes"], result["spectrum"])
```

- **示例输出：**

```text
[4] [None, 6.0, 15.0, 14.0, None]
```

## 18) `dnakit.population.read_sfs_file()`

- **作用：** 读取 dadi 的数组或频率谱文本文件。
- **处理方法：** 调用 dadi 的 `Numerics.array_from_file()`。第一行非注释行给出形状，之后的行为数值。以 `#` 开头的行作为注释返回。
- **API：** `dnakit.population.read_sfs_file(path)`
- **输入：** `path` 是文本文件路径。
- **输出：** 含 `shape`、`spectrum`、`comments` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import read_sfs_file

result = read_sfs_file("spectrum.sfs")
print(result["shape"], result["spectrum"], result["comments"])
```

- **示例输出：**

```text
[2, 2] [[1.5, 2.0], [3.0, 4.0]] ['demo']
```

## 19) `dnakit.population.write_sfs_file()`

- **作用：** 把数组或频率谱写成 dadi 文本文件。
- **处理方法：** 调用 dadi 的 `Numerics.array_to_file()`。文件先写形状，再写数值；`comment_lines` 写在文件开头的 `# ` 之后。不覆盖已有文件。
- **API：** `dnakit.population.write_sfs_file(spectrum, destination, *, precision=16, comment_lines=())`
- **输入：** `spectrum` 为数值列表；`destination` 为新文件路径；`precision` 为数值位数（1 至 17）。
- **输出：** 含 `destination`、`bytes`、`shape`、`num_values` 和 `provenance` 的字典。
- **示例代码：**

```python
from pathlib import Path
from dnakit.population import write_sfs_file

Path("spectrum_out.sfs").unlink(missing_ok=True)  # the function never overwrites an existing file
result = write_sfs_file([[1.5, 2.0], [3.0, 4.0]], "spectrum_out.sfs", comment_lines=["demo"])
print(result["shape"], result["num_values"], Path("spectrum_out.sfs").read_text().splitlines()[0])
```

- **示例输出：**

```text
[2, 2] 4 # demo
```

## 20) `dnakit.population.apply_ancestral_misidentification()`

- **作用：** 对频率谱施加祖先状态误判。
- **处理方法：** 调用 dadi 的 `Numerics.apply_anc_state_misid()`，结果为 `(1 - p) * fs + p * 反转(fs)`，其中 `p` 是祖先状态被误判的位点比例。
- **API：** `dnakit.population.apply_ancestral_misidentification(spectrum, p_misid)`
- **输入：** `spectrum` 为频率谱列表；`p_misid` 在 0 到 1 之间。
- **输出：** 含 `p_misid`、`sample_sizes`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import apply_ancestral_misidentification

result = apply_ancestral_misidentification([0, 1, 2, 3, 0], 0.1)
print([None if v is None else round(v, 4) for v in result["spectrum"]])
```

- **示例输出：**

```text
[None, 1.2, 2.0, 2.8, None]
```

## 21) `dnakit.population.expected_sfs_1d()`

- **作用：** 单群体人口模型的期望位点频率谱。
- **处理方法：** 用 dadi 的一维模型族计算期望频率谱。`model` 可以是 `two_epoch`（代表）、`growth`、`bottlegrowth`、`bottlegrowth_1d`、`three_epoch` 或 `three_epoch_inbreeding`。`grid_points` 是积分网格，增大后结果趋于精确值。
- **API：** `dnakit.population.expected_sfs_1d(model, params, sample_size, *, grid_points=1000)`
- **输入：** `model`；`params` 是该模型的参数（dadi 顺序）；`sample_size` 为 2 至 2000；`grid_points` 为 50 至 10000。
- **输出：** 含 `model`、`params`、`sample_sizes`、`spectrum`（端点为 None）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_1d

result = expected_sfs_1d("two_epoch", [2.0, 0.5], 10, grid_points=50)
print(result["model"], [round(v, 4) for v in result["spectrum"][1:4]])
```

- **示例输出：**

```text
two_epoch [1.6053, 0.6722, 0.3969]
```

## 22) `dnakit.population.expected_sfs_2d()`

- **作用：** 两群体人口模型的期望联合频率谱。
- **处理方法：** 用 dadi 的二维模型族计算。`model` 默认为 `split_mig`（代表），其他名称如 `asym_mig`、`anc_sym_mig`、`IM` 也可用，接口相同 `(params, ns, pts)`。
- **API：** `dnakit.population.expected_sfs_2d(model, params, sample_sizes, *, grid_points=60)`
- **输入：** `model`；`params`；`sample_sizes` 为两个群体的样本量（2 至 200）；`grid_points` 为 20 至 400。
- **输出：** 含 `model`、`params`、`sample_sizes`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_2d

result = expected_sfs_2d("split_mig", [1.0, 1.0, 0.5, 1.0], [4, 4], grid_points=30)
print(result["sample_sizes"], round(result["spectrum"][1][1], 4))
```

- **示例输出：**

```text
[4, 4] 0.1778
```

## 23) `dnakit.population.expected_sfs_3d()`

- **作用：** 三群体人口模型的期望联合频率谱。
- **处理方法：** 用 dadi 的三维模型族计算。`model` 默认为 `split_nomig`（代表），`out_of_africa` 为 Gutenkunst 等（2009）的模型。网格较小，因为三维积分的代价随 `grid_points` 的三次方增长。
- **API：** `dnakit.population.expected_sfs_3d(model, params, sample_sizes, *, grid_points=30)`
- **输入：** `model`；`params`；`sample_sizes` 为三个群体的样本量（2 至 40）；`grid_points` 为 20 至 60。
- **输出：** 含 `model`、`params`、`sample_sizes`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_3d

result = expected_sfs_3d("split_nomig", [1.0, 1.0, 0.5, 1.0, 0.5, 0.5], [2, 2, 2], grid_points=20)
print(result["model"], round(result["spectrum"][1][1][1], 4))
```

- **示例输出：**

```text
split_nomig 0.0069
```

## 24) `dnakit.population.expected_sfs_selection()`

- **作用：** 带选择（DFE 模型）的人口模型期望频率谱。
- **处理方法：** 调用 dadi 的 `DFE.DemogSelModels` 模型。`model` 默认为 `split_mig_sel`（代表）；`*_single_gamma` 变体让各群体共享同一选择系数。`sample_sizes` 为一或两个样本量，视模型而定。
- **API：** `dnakit.population.expected_sfs_selection(model, params, sample_sizes, *, grid_points=40)`
- **输入：** `model`；`params`（如 `nu1, nu2, T, m, gamma1, gamma2`）；`sample_sizes`；`grid_points` 为 20 至 400。
- **输出：** 含 `model`、`params`、`sample_sizes`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_selection

result = expected_sfs_selection("split_mig_sel", [1.0, 1.0, 0.5, 1.0, 1.0, 1.0], [4, 4], grid_points=20)
print(result["model"], round(result["spectrum"][1][1], 4))
```

- **示例输出：**

```text
split_mig_sel 0.2278
```

## 25) `dnakit.population.fit_demographic_model()`

- **作用：** 用最大似然把人口模型拟合到观测频率谱。
- **处理方法：** 调用 dadi 的 `Inference.optimize_log()` 系列。`method` 选择优化器：`bfgs`（对数参数上的 BFGS）、`lbfgsb`、`nelder-mead` 或 `powell`。`multinom=True` 时 dadi 自动缩放 theta，`p0` 不含 theta。结果中的对数似然是拟合参数处的值。
- **API：** `dnakit.population.fit_demographic_model(data, *, model, p0, lower_bound=None, upper_bound=None, method='bfgs', multinom=True, grid_points=40, maxiter=50)`
- **输入：** `data` 为观测频率谱（一或二群体）；`model` 按 `data` 的维数从对应模型族中选择；`p0` 为正的初值；`lower_bound`、`upper_bound` 可选。
- **输出：** 含 `model`、`method`、`params`、`log_likelihood`、`grid_points`、`maxiter` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_1d, fit_demographic_model

expected = expected_sfs_1d("two_epoch", [2.0, 0.5], 10, grid_points=50)["spectrum"]
data = [0.0 if v is None else v * 100 for v in expected]  # counts; the corners are 0
result = fit_demographic_model(
    data, model="two_epoch", p0=[1.0, 0.3], lower_bound=[0.01, 0.001], upper_bound=[10, 5], grid_points=50
)
print([round(v, 3) for v in result["params"]], round(result["log_likelihood"], 3))
```

- **示例输出：**

```text
[1.999, 0.5] -23.198
```

## 26) `dnakit.population.sfs_residuals()`

- **作用：** 模型频率谱与观测计数之间逐项的残差。
- **处理方法：** 调用 dadi 的 `Inference.Anscombe_Poisson_residual()`（默认）或 `linear_Poisson_residual()`，即 `(模型 - 数据) / sqrt(模型)`。两个数组形状必须相同。端点为 None，`sum_of_squares` 只计未屏蔽的项。
- **API：** `dnakit.population.sfs_residuals(model, data, *, residual='anscombe')`
- **输入：** `model` 为模型频率谱；`data` 为观测计数（非负）；`residual` 为 `anscombe` 或 `poisson`。
- **输出：** 含 `residual`、`residuals`、`num_entries`、`sum_of_squares` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import sfs_residuals

result = sfs_residuals([1.0, 4.0, 9.0, 4.0, 1.0], [1.0, 2.0, 2.0, 2.0, 1.0], residual="poisson")
print([None if v is None else round(v, 4) for v in result["residuals"]], round(result["sum_of_squares"], 4))
```

- **示例输出：**

```text
[None, 1.0, 2.3333, 1.0, None] 7.4444
```

## 27) `dnakit.population.parameter_uncertainty_fim()`

- **作用：** 由 Fisher 信息矩阵估计参数的标准误。
- **处理方法：** 调用 dadi 的 `Godambe.FIM_uncert()`。模型先经 `Numerics.make_extrap_func()` 在 `grid_points` 上外推，再由对数似然的数值 Hessian 得到 Fisher 信息。`multinom=True` 时 theta 也被估计，其标准误在最后一项。
- **API：** `dnakit.population.parameter_uncertainty_fim(data, *, model, p0, grid_points=(20, 30, 40), log=False, multinom=True, eps=0.01)`
- **输入：** `data` 为观测计数；`model` 为模型名；`p0` 为正的参数值；`grid_points` 为 1 至 5 个网格大小；`eps` 为有限差分的相对步长。
- **输出：** 含 `params`、`grid_points`、`uncertainties` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_1d, parameter_uncertainty_fim

expected = expected_sfs_1d("two_epoch", [2.0, 0.5], 6, grid_points=50)["spectrum"]
data = [0.0 if v is None else v * 100 for v in expected]
result = parameter_uncertainty_fim(data, model="two_epoch", p0=[2.0, 0.5])
print([round(v, 3) for v in result["uncertainties"]])
```

- **示例输出：**

```text
[0.661, 0.525, 16.189]
```

## 28) `dnakit.population.annotate_snps_from_annovar()`

- **作用：** 用 ANNOVAR 的外显子注释文件保留某一类变异的 SNP。
- **处理方法：** 调用 dadi 的 `Misc.annotate_from_annovar()`。若注释行的变异类型文字包含 `variant_type`，则保留该位点。返回的 SNP 编号已排序。
- **API：** `dnakit.population.annotate_snps_from_annovar(path, annovar_path, *, variant_type)`
- **输入：** `path` 为 SNP 文件；`annovar_path` 为 `.exonic_variant_function` 文件；`variant_type` 如 `nonsynonymous`。
- **输出：** 含 `num_snps_total`、`num_snps_kept`、`kept_snp_ids` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import annotate_snps_from_annovar

result = annotate_snps_from_annovar(
    "chunks.txt", "sample.exonic_variant_function", variant_type="nonsynonymous"
)
print(result["num_snps_total"], result["num_snps_kept"], result["kept_snp_ids"])
```

- **示例输出：**

```text
6 1 ['chr1_100']
```

## 29) `dnakit.population.trinucleotide_frequency_array()`

- **作用：** 把 64 个三核苷酸频率转换为 dadi 的 64 项数组。
- **处理方法：** 调用 dadi 的 `Misc.tri_freq_dict_to_array()`。dadi 的顺序是每个位置按 `C, G, T, A` 排列，第一位变化最慢；`trinucleotides` 给出该顺序。
- **API：** `dnakit.population.trinucleotide_frequency_array(frequencies)`
- **输入：** `frequencies` 是映射，键为 A、C、G、T 组成的 64 个三核苷酸，值为非负频率。
- **输出：** 含 `trinucleotides`、`values` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import trinucleotide_frequency_array

bases = "ACGT"
freqs = {a + b + c: 0.0 for a in bases for b in bases for c in bases}
freqs["GTA"] = 1.0
result = trinucleotide_frequency_array(freqs)
position = result["values"].index(1.0)
print(len(result["values"]), result["trinucleotides"][position])
```

- **示例输出：**

```text
64 GTA
```

## 30) `dnakit.population.expected_sfs_from_demes()`

- **作用：** 由 demes 人口图（YAML 文件）计算指定群体的期望频率谱。
- **处理方法：** 调用 dadi 的 `Demes.SFS()`，图由 `demes` 包读取（需单独安装；它不在 `bio-dadi` 额外项中）。`sampled_demes` 给出一至三个群体名，`sample_sizes` 给出各群体的样本量；`theta` 为缩放系数。
- **API：** `dnakit.population.expected_sfs_from_demes(graph, sampled_demes, sample_sizes, *, grid_points=40, theta=1.0)`
- **输入：** `graph` 是 demes YAML 文件；`sampled_demes` 是群体名列表；`sample_sizes` 与之一一对应；`grid_points` 为 20 至 400（三群体时为 20 至 60）。
- **输出：** 含 `demes`、`sample_sizes`、`theta`、`spectrum` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import expected_sfs_from_demes

result = expected_sfs_from_demes("model.yaml", ["A", "B"], [4, 4], grid_points=40)
print(result["sample_sizes"], round(result["spectrum"][1][1], 4))
```

- **示例输出：**

```text
[4, 4] 0.047
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-tskit]"    # tskit 相关函数
python -m pip install "dnakit[bio-dadi]"     # dadi 相关函数；cyvcf2 读取器另需 cyvcf2
```

项目地址：[tskit](https://tskit.dev/tskit)、[dadi](https://github.com/RyanGutenkunst/dadi)、[cyvcf2](https://github.com/brentp/cyvcf2)。

## 验证方式

测试文件是 `others/tests/unit/test_full_G13.py`，共 30 项，每个包装函数至少一项。
- tskit 的 FASTA、NEXUS、ms 与文本表输出与直接调用 tskit 的字节级结果一致；`complement_intervals()`、`map_allele_indices()` 与手工计算一致；`tree_to_newick()` 的支长与由节点时间算出的值一致；`resolve_polytomies()` 把四叉星树拆为二叉树并新增 2 个节点。
- tskit 的拓扑计数、`random_tree_newick()` 和 PLINK `.fam` 解析与直接调用 tskit 的结果一致。
- dadi 的 SNP 极化数、构型计数、VCF 读取（内置与 cyvcf2 两种读取器）与直接调用 dadi 的结果一致。
- 自助抽样、频率谱合并、模型期望频率谱、模型拟合、残差与 Fisher 信息均与直接调用 dadi 的结果一致；Poisson 残差与手工公式一致；祖先误判与手工公式一致；三核苷酸数组的顺序与 dadi 一致。
- 参数不合法、输出文件已存在、模型名不在名单中等情况均报出 `ConfigurationError`。
测试不访问网络，并在 `tmp_path` 中创建输入；依赖缺失时用 `pytest.importorskip` 跳过。

## 当前范围

这一批包装 tskit 的 12 个函数（其中 6 个为家族代表，另有 `write_nexus`、`tree_count_topologies`、`generate_*` 与 `parse_*` 的其他成员作为家族变体记录在台账中）和 dadi 的 18 个函数。族内的其他模型或解析器通过参数选择，例如 `expected_sfs_2d(model=...)`、`read_text_table(table=...)`、`fit_demographic_model(method=...)`。
未包装的主要部分包括：tskit 的绘图与 HTML 渲染、命令行入口、元数据与 provenance 辅助函数；dadi 的绘图、Triallele 与 TwoLocus 模型（在当前安装的 dadi 2.4.4 与 NumPy 2 下无法运行，前者缺少 `dadi.tridiag`，后者使用了已删除的 `np.float`）、LowPass 低覆盖校正、Godambe 的 GIM 与 bootstrap 变体，以及基于 VCF 子抽样的自助抽样；demes 图的输出、切片与拟合函数（`Demes.output`、`DemesUtil`、`Demes.Inference`）也未包装。
`dnakit.population.backends` 中已有的树序列构建、读写和群体统计函数，以及中性频率谱函数，在台账中标记为 `already-wrapped` 或 `duplicate-of`，本批没有重复包装。
