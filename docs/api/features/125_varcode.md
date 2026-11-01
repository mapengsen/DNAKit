# varcode：变异效应、融合与剪接结果

这一组函数是 `dnakit.molecular_lab` 中基于可选依赖的适配器，封装了 [varcode](https://github.com/openvax/varcode)（Apache-2.0，预测变异对转录本和蛋白的影响）和 [pyensembl](https://github.com/openvax/pyensembl)（Apache-2.0，varcode 依赖的转录本模型数据库），并用 [pysam](https://github.com/pysam-developers/pysam) 读取基因组 FASTA。函数可以预测小变异和结构变异对转录本的影响，把种系与体细胞变异应用到转录本上得到突变 cDNA 和蛋白，枚举并实现剪接方案与剪接结果，导入观察到的 RNA 融合模型（包括 ExACTO 表），列出与变异重叠的基因，从 VCF、MAF 和 JSON 文件加载变异，检查 VCF 样本是否来自同一供体，并对结构变异检出结果分组。它们不做变异检出、不比对读段，也不测量表达量。

大多数函数接收本地注释（`gtf`）和基因组 FASTA（`genome_fasta`），用 `build_local_reference()` 在本地建立 pyensembl 参考，因此不会下载任何数据。参考按内容摘要只建一次索引，并缓存在 `cache_dir` 下（默认 `~/.cache/dnakit/varcode_reference`）；下面的示例使用工作目录中的 `cache_dir="cache"`。变异是带有 `contig`、1 起始的 `start`、`ref` 和 `alt` 的映射（碱基为 A、C、G、T、N；插入缺失变异让 `ref` 与 `alt` 之一为空）。结构变异是带有 `contig`、1 起始且含两端的 `start` 与 `end`，以及 `sv_type`（`DEL`、`DUP`、`INV`、`INS` 或 `CNV`）的映射。转录本用 `transcript_id` 指定。每个结果都是可序列化为 JSON 的字典，含 `operation`、`status`、`backend`、`package` 和 `coordinates`；效应对象用类名和描述概括，并附上变异、转录本 ID 和基因名。参数无效会抛出 `ConfigurationError`；缺少依赖包会抛出带安装命令的 `BackendUnavailableError`。

函数通过 `from dnakit.molecular_lab import ...` 导入。依赖都是可选的，只有调用对应函数时才会导入：varcode 函数需要 `dnakit[bio-varcode]` 额外依赖，使用本地参考的函数还需要 `dnakit[bio-pyensembl]` 和 `dnakit[bio-reads]`（pysam）。根据包元数据，varcode 与 pyensembl 均为 Apache-2.0 许可，因此这两个包没有传染性许可证问题；DNAKit 不复制它们的任何代码。

## 示例数据

示例共用一个很小的合成基因组：一条 450 bp 的序列，名称为 `1`（用 `1` 而不是 `chr1`，因为 varcode 的 VCF 加载器会去掉 `chr` 前缀，基于 VCF 的示例只能匹配没有前缀的参考染色体名）。其上有两个正链基因。`GENEA`（转录本 `T1`）有外显子 1-60、101-150 和 191-240，编码序列从位置 10 到 228-230 的终止密码子，编码一个以 `MPTVYG` 开头的 46 个氨基酸的蛋白。`GENEB`（转录本 `T2`）有外显子 281-330 和 371-420，编码序列从 291 到 409-411 的终止密码子。`T1` 第二个内含子的供体（位置 151-152）是 `AT` 而不是标准的 `GT`，因此剪接规划会看到这里的供体丢失。下面用到的变异是 `1:13 C>G`（`p.P2A`）、`1:20 TC>T`（移码）和 `1:61 G>A`（第一个内含子的供体）。这些文件是为本页合成的，不是下载的；请在包含全部文件的目录中运行示例。

**示例文件：** `genome.fa`

```text
>1
CCGACTGCAATGCCAACAGTCTATGGAGCTGTACATTACTATCGCATACGTTTAGAAGCC
GTTTGAGGGTGAGAGCCCATTCTAGACTTTGAGAGTCGAGGGTCCTGAGAAATGGAAAAA
AGCTGTGACTTTCTCTTATGACTACACAGAATTTAGCGCCCACCCCAGAAGGGGTTGCTC
TAAGCAACAGGCTCCGCAAACCCGATCACACGACCAAGACGAGACACTAAAAGTTTATAG
ATCCTAACATTCTAGAGCCGATAGAGCTAGCTTGACTGTCGACGCCATTAATGAGCCGAC
TGCACCCTACCGTAACAGGCAGACCCGCTCGTGTCTGAGAACTCCTCCTGCGACGGCAAG
CGAGCCGAAGGCCGGACTAAAGAGGTTCGTGCTTTAGAGATACAGGTCTGACGTTCGCTA
AGAACCCTAGTTCAAAGCATCGCAGTAGCT
```

**示例文件：** `genes.gtf`

```text
1	demo	gene	1	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding";
1	demo	transcript	1	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	exon	1	60	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	exon	101	150	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	exon	191	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "3";
1	demo	CDS	10	60	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	CDS	101	150	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	CDS	191	227	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "3";
1	demo	start_codon	10	12	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	stop_codon	228	230	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	gene	281	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding";
1	demo	transcript	281	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
1	demo	exon	281	330	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	exon	371	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	CDS	291	330	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	CDS	371	408	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	start_codon	291	293	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
1	demo	stop_codon	409	411	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
```

这里只保留 varcode 与 pyensembl 所需的最少属性：`gene_id`、`gene_name`、`gene_biotype`、`transcript_id`、`transcript_biotype` 和 `exon_number`。`start_codon` 与 `stop_codon` 行把转录本标记为完整；CDS 行不含终止密码子。

**示例文件：** `phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=450>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	13	.	C	G	.	PASS	.	GT:PS	0|1:13
1	17	.	C	T	.	PASS	.	GT:PS	0|1:13
```

带相位信息的 VCF，样本为 `S1`：两个变异都在相位集 13 的第二条单倍型上（`0|1`）。

**示例文件：** `structures.tsv`

```text
transcript_model_id	reference_transcript_ids	index	read_start	read_end	sequence	type	kind	context	chromosome_1	position_1	strand_1	chromosome_2	position_2	strand_2	transcript_id_1	transcript_id_2
M1	T1,T2	0	0	59	CCGACTGCAATGCCAACAGTCTATGGAGCTGTACATTACTATCGCATACGTTTAGAAGCC	base	match	exon	1	1	+	1	60	+	T1	T1
M1	T1,T2	1	60	109	GGTCCTGAGAAATGGAAAAAAGCTGTGACTTTCTCTTATGACTACACAGA	base	match	exon	1	101	+	1	150	+	T1	T1
M1	T1,T2	2	110	159	GACGCCATTAATGAGCCGACTGCACCCTACCGTAACAGGCAGACCCGCTC	base	match	exon	1	281	+	1	330	+	T2	T2
```

ExACTO 转录本结构表，含一个模型 `M1`，由 `T1` 的外显子 1 和 2（位置 1-60 与 101-150）连接 `T2` 的外显子 1（281-330）构成。各列以制表符分隔。

**示例文件：** `integrated.tsv`

```text
transcript_model_id	reference_transcript_ids	rna_variant_call_id	dna_variant_call_id
M1	T1,T2	R1	D1
```

ExACTO 整合表，把模型 `M1` 的 RNA 检出 `R1` 与 DNA 检出 `D1` 关联起来。

**示例文件：** `donor1.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D1
chr1	1000	.	A	G	.	PASS	.	GT	0/1
chr1	2000	.	C	T	.	PASS	.	GT	1/1
chr1	3000	.	G	A	.	PASS	.	GT	0/0
chr1	4000	.	T	C	.	PASS	.	GT	0/1
chr1	5000	.	A	C	.	PASS	.	GT	0/0
```

**示例文件：** `donor2.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D2
chr1	1000	.	A	G	.	PASS	.	GT	0/1
chr1	2000	.	C	T	.	PASS	.	GT	1/1
chr1	3000	.	G	A	.	PASS	.	GT	0/0
chr1	4000	.	T	C	.	PASS	.	GT	0/1
chr1	5000	.	A	C	.	PASS	.	GT	0/0
```

**示例文件：** `donor3.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D3
chr1	1000	.	A	G	.	PASS	.	GT	1/1
chr1	2000	.	C	T	.	PASS	.	GT	0/0
chr1	3000	.	G	A	.	PASS	.	GT	0/1
chr1	4000	.	T	C	.	PASS	.	GT	0/0
chr1	5000	.	A	C	.	PASS	.	GT	1/1
```

三个种系 VCF 各有五个 SNP。`donor1.vcf` 与 `donor2.vcf` 的基因型相同，`donor3.vcf` 的基因型不同。

**示例文件：** `calls.csv`

```text
call_id,caller,sample,build,chrom,pos,ref,alt,info
c1,callerA,S1,GRCh38,1,10000,N,<DEL>,END=10500;SVTYPE=DEL
c2,callerB,S1,GRCh38,1,10020,N,<DEL>,END=10490;SVTYPE=DEL
c3,callerA,S1,GRCh38,2,5000,N,<DEL>,END=5300;SVTYPE=DEL
```

三条结构变异检出：`c1` 与 `c2` 是染色体 1 上的缺失，断点最多相差 20 个碱基；`c3` 是染色体 2 上的缺失。

## 1) `dnakit.molecular_lab.build_local_reference()`

- **作用：** 把本地 GTF 与基因组 FASTA 建成 varcode 可用的 pyensembl 参考，无需下载。
- **处理方法：** 读取 GTF（每条非注释行必须有九个制表符分隔的字段，exon 与 CDS 行必须有 `transcript_id`；缺失的 CDS `protein_id` 被设为转录本 ID，缺失的 `exon_id` 会自动生成）。用 pysam 从基因组 FASTA 中截取外显子，按转录本链拼接成 cDNA FASTA，并把 CDS 行翻译成蛋白 FASTA；然后用这些文件构建 `pyensembl.Genome` 并调用 `index()`。索引保存在 `cache_dir/<摘要>` 下，摘要是 GTF 内容加上 FASTA 路径、大小和修改时间的 SHA-256。本进程中已建好的索引会被复用，并报告其原始目录。
- **API：** `build_local_reference(gtf: 'str', genome_fasta: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的文件；`cache_dir` 是索引目录。每个转录本必须位于同一条染色体和同一条链上，其染色体必须在 FASTA 中，且至少有一个转录本带有 exon 行（否则抛出 `ConfigurationError`）。
- **输出：** `reference_name`（`dnakit_local`）、`annotation_name`、`cache_directory`、`reused`（本进程中已加载索引时为真），以及计数 `contigs`、`genes`、`transcripts` 和 `exons`。
- **示例代码：**

```python
from dnakit.molecular_lab import build_local_reference

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

first = build_local_reference(**ref)
again = build_local_reference(**ref)
print(first["reference_name"], first["contigs"], first["genes"], first["transcripts"], first["exons"])
print(first["reused"], again["reused"])
```

- **示例输出：**

```text
dnakit_local 1 2 2 5
False True
```

## 2) `dnakit.molecular_lab.varcode_apply_germline_to_transcript()`

- **作用：** 把种系变异应用到转录本上，返回得到的患者转录本。
- **处理方法：** 用种系变异构建 varcode 的 `GermlineContext`，再调用 `varcode.germline.apply_germline_to_transcript`。不给体细胞变异时，应用所有位于该转录本外显子内的种系变异；给出体细胞变异时，只应用落在其窗口（即 `varcode_default_germline_window()` 报告的窗口）内的种系变异。如果没有合格的种系变异，或编辑相互冲突、落在 CDS 之外，varcode 返回空，包装函数不会报错，而是返回一个空的突变记录（空字符串、序列为 `None`、没有编辑）。
- **API：** `varcode_apply_germline_to_transcript(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', germline_variants: 'Sequence[Mapping[str, Any]]', *, somatic_variant: 'Mapping[str, Any] | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`germline_variants` 是最多 1000 个变异映射的列表；`somatic_variant` 是可选的变异映射。
- **输出：** `transcript_id` 与 `mutant_transcript`，后者是含 `transcript_id`、`annotator`、`cdna_sequence`、`mutant_protein_sequence` 和 `edits` 的字典（每个编辑含 `cdna_start`、`cdna_end`、`alt_bases` 和 `source_variant`；cDNA 偏移从 0 起，左闭右开）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_apply_germline_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

somatic = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
near = {"contig": "1", "start": 14, "ref": "C", "alt": "T"}
far = {"contig": "1", "start": 17, "ref": "C", "alt": "T"}

result = varcode_apply_germline_to_transcript(
    transcript_id="T1", germline_variants=[near], somatic_variant=somatic, **ref
)
mutant = result["mutant_transcript"]
print(mutant["annotator"], mutant["mutant_protein_sequence"][:6], len(mutant["edits"]))

outside = varcode_apply_germline_to_transcript(
    transcript_id="T1", germline_variants=[far], somatic_variant=somatic, **ref
)
print(outside["mutant_transcript"]["cdna_sequence"])
```

- **示例输出：**

```text
protein_diff MLTVYG 1
None
```

## 3) `dnakit.molecular_lab.varcode_apply_variant_to_transcript()`

- **作用：** 把一个变异应用到转录本上，返回突变 cDNA 与蛋白。
- **处理方法：** 用变异和转录本模型调用 `varcode.mutant_transcript.apply_variant_to_transcript`，并概括得到的 `MutantTranscript`：剪接后的突变 cDNA、从 CDS 起点翻译出的蛋白，以及已应用的编辑。
- **API：** `varcode_apply_variant_to_transcript(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `variant`、`transcript_id` 和 `mutant_transcript`（含 `transcript_id`、`annotator`、`cdna_sequence`、`mutant_protein_sequence` 和 `edits`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_apply_variant_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variant = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
result = varcode_apply_variant_to_transcript(variant=variant, transcript_id="T1", **ref)
mutant = result["mutant_transcript"]
print(len(mutant["cdna_sequence"]), mutant["mutant_protein_sequence"][:8], len(mutant["mutant_protein_sequence"]))
print(mutant["edits"][0]["cdna_start"], mutant["edits"][0]["cdna_end"], mutant["edits"][0]["alt_bases"])
```

- **示例输出：**

```text
160 MATVYGAV 46
12 13 G
```

## 4) `dnakit.molecular_lab.varcode_apply_variants_to_transcript()`

- **作用：** 把多个变异一起应用到同一个转录本上，返回得到的突变体。
- **处理方法：** 用全部变异一次调用 `varcode.mutant_transcript.apply_variants_to_transcript`，所有编辑合并到同一个突变 cDNA 中，蛋白由该 cDNA 翻译得到（示例中移码使蛋白在替换之后被截短）。
- **API：** `varcode_apply_variants_to_transcript(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variants` 是最多 1000 个变异映射的列表。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `transcript_id`、`variants`（varcode 规范化后的变异）和 `mutant_transcript`（字段与 `varcode_apply_variant_to_transcript()` 相同）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_apply_variants_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = [{"contig": "1", "start": 13, "ref": "C", "alt": "G"}, {"contig": "1", "start": 20, "ref": "TC", "alt": "T"}]
result = varcode_apply_variants_to_transcript(variants=variants, transcript_id="T1", **ref)
mutant = result["mutant_transcript"]
print([(v["start"], v["ref"], v["alt"]) for v in result["variants"]])
print(mutant["mutant_protein_sequence"], len(mutant["edits"]))
```

- **示例输出：**

```text
[(13, 'C', 'G'), (21, 'C', '')]
MATVMELYITIAYV 2
```

## 5) `dnakit.molecular_lab.varcode_candidates_from_effects()`

- **作用：** 预测变异的效应，并把每个效应包装成带来源标签的效应候选。
- **处理方法：** 效应对象无法以 JSON 形式传入，因此包装函数先预测该变异对其重叠的每个转录本的效应（`varcode.effects.predict_variant_effects`），再用 `varcode.candidates_from_effects(effects, source=source)` 包装。
- **API：** `varcode_candidates_from_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', *, source: 'str' = 'varcode', cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`source` 是非空标签（默认 `varcode`）。
- **输出：** `variant`、`source`、`candidates` 和 `count`。每个候选含 `source`、`evidence`（字典）和 `effect`（`class`、`description`、`variant`、`transcript_id`、`gene_name`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_candidates_from_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_candidates_from_effects(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, source="demo", **ref)
print(result["count"], result["source"])
candidate = result["candidates"][0]
print(candidate["source"], candidate["evidence"], candidate["effect"]["class"], candidate["effect"]["description"])
```

- **示例输出：**

```text
1 demo
demo {} Substitution p.P2A
```

## 6) `dnakit.molecular_lab.varcode_classify_from_protein_diff()`

- **作用：** 通过比较参考蛋白与突变蛋白，对变异造成的蛋白改变分类。
- **处理方法：** 用变异、转录本、两条蛋白字符串和 `length_delta` 调用 `varcode.effects.classify.classify_from_protein_diff`。两条蛋白相同则为 `Silent`；突变蛋白不再以 `M` 开头则为 `StartLoss`；否则裁掉共有的两端，对差异分类（替换、插入、缺失、提前终止、终止丢失）。`length_delta` 不是 3 的倍数且蛋白尾部改变时判为移码。
- **API：** `varcode_classify_from_protein_diff(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', ref_protein: 'str', mut_protein: 'str', *, length_delta: 'int' = 0, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`ref_protein` 与 `mut_protein` 是非空字符串；`length_delta` 是编辑造成的核苷酸净长度变化（默认 0）。
- **输出：** `variant`、`transcript_id` 和 `effect`（`class`、`description`、`variant`、`transcript_id`、`gene_name`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_classify_from_protein_diff

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variant = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
same = varcode_classify_from_protein_diff(
    variant=variant, transcript_id="T1", ref_protein="MPTV", mut_protein="MPTV", **ref
)
print(same["effect"]["class"], same["effect"]["description"])
stop = varcode_classify_from_protein_diff(
    variant=variant, transcript_id="T1", ref_protein="MPTV", mut_protein="MP", length_delta=-6, **ref
)
print(stop["effect"]["class"], stop["effect"]["description"])
```

- **示例输出：**

```text
Silent p.M1=
PrematureStop p.T3*
```

## 7) `dnakit.molecular_lab.varcode_codon_table_for_transcript()`

- **作用：** 选择 varcode 为某个转录本使用的密码子表。
- **处理方法：** 调用 `varcode.effects.codon_tables.codon_table_for_transcript`（也以 `varcode.annotators.protein_diff.codon_table_for_transcript` 暴露）。位于染色体名为 `mt`、`chrmt`、`m`、`chrm`、`mtdna` 或 `mito`（不区分大小写）的转录本使用脊椎动物线粒体密码子表（NCBI 表 2），其余转录本使用标准密码子表。
- **API：** `varcode_codon_table_for_transcript(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `transcript_id`、`contig` 和 `codon_table`（含 `name`、`ncbi_id`，以及排序后的 `start_codons` 与 `stop_codons`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_codon_table_for_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_codon_table_for_transcript(transcript_id="T1", **ref)
table = result["codon_table"]
print(result["contig"], table["name"], table["ncbi_id"], table["stop_codons"])
```

- **示例输出：**

```text
1 Standard 1 ['TAA', 'TAG', 'TGA']
```

## 8) `dnakit.molecular_lab.varcode_default_germline_window()`

- **作用：** 返回为寻找可与体细胞变异定相的种系变异而搜索的基因组窗口。
- **处理方法：** 用体细胞变异和转录本调用 `varcode.germline.default_germline_window`。窗口从变异前 3 个碱基到变异后 3 个碱基；如果转录本的某个外显子边界距变异不超过 6 个碱基，窗口会扩展到该边界两侧各 6 个碱基，以覆盖供体和受体信号。
- **API：** `varcode_default_germline_window(gtf: 'str', genome_fasta: 'str', somatic_variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`somatic_variant` 是变异映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** 窗口的 `contig`、`start` 和 `end`（1 起始，含两端）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_default_germline_window

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

inside = varcode_default_germline_window(somatic_variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(inside["contig"], inside["start"], inside["end"])
near_edge = varcode_default_germline_window(
    somatic_variant={"contig": "1", "start": 58, "ref": "C", "alt": "G"}, transcript_id="T1", **ref
)
print(near_edge["contig"], near_edge["start"], near_edge["end"])
```

- **示例输出：**

```text
1 10 16
1 54 66
```

## 9) `dnakit.molecular_lab.varcode_infer_genome()`

- **作用：** 根据发布号或参考名称解析出 Ensembl 发布版本。
- **处理方法：** 调用 `varcode.reference.infer_genome`。整数表示人类的 Ensembl 发布号；字符串是 `GRCh38` 或 `GRCh38:93` 这样的参考名称。`hg19` 等 UCSC 名称会转换为对应的 Ensembl 名称（GRCh37）。发布对象是惰性创建的，因此调用本身不需要数据文件。
- **API：** `varcode_infer_genome(reference: 'str | int', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）。可离线运行。解析名称不会下载任何东西；下载发布数据是另一步（见 `varcode_install_reference_data()`）。
- **输入：** `reference` 是整数或非空字符串。布尔值、浮点数或未知名称会抛出 `ConfigurationError`。
- **输出：** `requested`、`reference_name`、`release`、`species`（拉丁名）和 `converted_from_ucsc`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_infer_genome

ucsc = varcode_infer_genome("hg19", verbose=False)
print(ucsc["reference_name"], ucsc["release"], ucsc["species"], ucsc["converted_from_ucsc"])
number = varcode_infer_genome(93, verbose=False)
print(number["reference_name"], number["release"], number["converted_from_ucsc"])
```

- **示例输出：**

```text
GRCh37 75 homo_sapiens True
GRCh38 93 False
```

## 10) `dnakit.molecular_lab.varcode_infer_genome_for_reference_name()`

- **作用：** 把参考名称（可带发布号）解析为 Ensembl 发布版本。
- **处理方法：** 调用 `varcode.reference.infer_genome_for_reference_name`。名称可以是普通的参考名称（`GRCh37`，对应其已知的最新发布）、带发布号的名称（`GRCh38:93`），或 UCSC 名称（`hg38`，会被转换）。不会下载任何数据。
- **API：** `varcode_infer_genome_for_reference_name(reference_name: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）。可离线运行。
- **输入：** `reference_name` 是非空字符串。
- **输出：** `requested`、`reference_name`、`release`、`species` 和 `converted_from_ucsc`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_infer_genome_for_reference_name

with_release = varcode_infer_genome_for_reference_name("GRCh38:93", verbose=False)
print(with_release["requested"], with_release["reference_name"], with_release["release"])
older = varcode_infer_genome_for_reference_name("GRCh37", verbose=False)
print(older["reference_name"], older["release"], older["converted_from_ucsc"])
```

- **示例输出：**

```text
GRCh38:93 GRCh38 93
GRCh37 75 False
```

## 11) `dnakit.molecular_lab.varcode_left_align_indels()`

- **作用：** 利用本地参考把插入缺失移到最左的等价位置。
- **处理方法：** 在本地参考上构建 `VariantCollection`，调用 `varcode.transforms.left_align_indels`，让每个插入或缺失沿重复序列向左滑动。替换变异不会移动。
- **API：** `varcode_left_align_indels(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variants` 是最多 1000 个变异映射的列表。
- **输出：** `variants`（对齐后的变异，含 `contig`、`start`、`ref`、`alt`）和 `count`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_left_align_indels

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

# The reference has CC at positions 13-14, so deleting the second C is the same as deleting the first.
variants = [{"contig": "1", "start": 14, "ref": "C", "alt": ""}]
result = varcode_left_align_indels(variants=variants, **ref)
print(result["count"], result["variants"])
```

- **示例输出：**

```text
1 [{'contig': '1', 'start': 13, 'ref': 'C', 'alt': ''}]
```

## 12) `dnakit.molecular_lab.varcode_make_rna_outcome()`

- **作用：** 把变异对某个转录本的效应包装成来自 RNA 的效应候选。
- **处理方法：** 先预测变异对该转录本的效应（`predict_variant_effect_on_transcript`），再把它连同来源标签、可选的转录本模型 ID 和可选的读段数传给 `varcode.make_rna_outcome`。不会读取任何 RNA 读段；读段数只作为证据保存，不会被当作似然。
- **API：** `varcode_make_rna_outcome(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, source: 'str' = 'rna', transcript_model_id: 'str | None' = None, read_count: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`source` 是非空标签（默认 `rna`），`transcript_model_id` 是可选标识符，`read_count` 是可选的非负整数。
- **输出：** `variant`、`transcript_id` 和 `candidate`（含 `source`、`evidence`（给出时包含模型 ID 与读段数）和 `effect`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_make_rna_outcome

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_make_rna_outcome(
    variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", source="rna", transcript_model_id="m1", read_count=12, **ref
)
candidate = result["candidate"]
print(candidate["source"], candidate["evidence"])
print(candidate["effect"]["class"], candidate["effect"]["description"])
```

- **示例输出：**

```text
rna {'transcript_model_id': 'm1', 'read_count': 12}
Substitution p.P2A
```

## 13) `dnakit.molecular_lab.varcode_predict_germline_aware_effect()`

- **作用：** 针对应用了种系变异的转录本，预测体细胞变异的效应。
- **处理方法：** 用种系变异构建 `GermlineContext`，再用体细胞变异、转录本、该上下文和注释器（`annotator` 为 None 时使用默认注释器）调用 `varcode.germline.predict_germline_aware_effect`。当某个种系变异与体细胞变异位于同一窗口且二者的相位未知时，varcode 返回 `PhaseCandidateSet`，即跨相位假设的可能效应集合（其简短描述以 `?` 开头）。包装函数只报告该效应的类名和简短描述，不列出各个假设。
- **API：** `varcode_predict_germline_aware_effect(gtf: 'str', genome_fasta: 'str', somatic_variant: 'Mapping[str, Any]', transcript_id: 'str', germline_variants: 'Sequence[Mapping[str, Any]]', *, annotator: 'str | None' = None, max_hypotheses: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`somatic_variant` 是变异映射；`germline_variants` 是变异映射列表；`annotator` 为 None 或 `fast`、`protein_diff`、`transcript_model`、`realized` 之一；`max_hypotheses` 为 None 或 1 到 10000 的整数。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `transcript_id`、`somatic_variant`、`germline_count` 和 `effect`（`class`、`description`、`variant`、`transcript_id`、`gene_name`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_germline_aware_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_predict_germline_aware_effect(
    somatic_variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", germline_variants=[{"contig": "1", "start": 14, "ref": "C", "alt": "T"}], **ref
)
print(result["germline_count"], result["effect"]["class"], result["effect"]["description"])
```

- **示例输出：**

```text
1 PhaseCandidateSet ?p.P2A
```

## 14) `dnakit.molecular_lab.varcode_predict_realized_effect()`

- **作用：** 预测体细胞变异对某个转录本的实际最可能效应，并给出备选结果。
- **处理方法：** 用体细胞变异（一起应用）、转录本、种系变异和 `max_hypotheses` 调用 `varcode.realized_effects.predict_realized_effect`。varcode 能解析相位时，种系变异与体细胞变异定相；它保留的备选结果（最多 `max_hypotheses` 个）作为 `candidates` 返回，其中若干项可能具有相同的描述。包装函数报告最优效应和最多 1000 个备选结果。
- **API：** `varcode_predict_realized_effect(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', transcript_id: 'str', *, germline_variants: 'Sequence[Mapping[str, Any]]' = (), max_hypotheses: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variants` 是最多 1000 个变异映射的列表；`germline_variants` 是变异映射列表（默认为空）；`max_hypotheses` 是 1 到 10000 的整数（默认 64）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `transcript_id`、`top_effect`（`class`、`description` 等）和 `candidates`（备选结果列表，含 `source`、`evidence` 和 `effect`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_realized_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

alone = varcode_predict_realized_effect(variants=[{"contig": "1", "start": 13, "ref": "C", "alt": "G"}], transcript_id="T1", **ref)
print(alone["top_effect"]["description"], len(alone["candidates"]))

phased = varcode_predict_realized_effect(
    variants=[{"contig": "1", "start": 13, "ref": "C", "alt": "G"}], transcript_id="T1", germline_variants=[{"contig": "1", "start": 14, "ref": "C", "alt": "T"}], **ref
)
print(phased["top_effect"]["description"], len(phased["candidates"]))
print(sorted({c["effect"]["description"] for c in phased["candidates"]}))
```

- **示例输出：**

```text
p.P2A 3
p.P2A 6
['p.L2V', 'p.P2A']
```

## 15) `dnakit.molecular_lab.varcode_predict_variant_effect_on_transcript()`

- **作用：** 返回一个变异对一个转录本的效应（例如移码）。
- **处理方法：** 用变异、转录本模型和所选注释器调用 `varcode.effects.predict_variant_effect_on_transcript`。参考碱基不匹配或其他上游 `ValueError` 会变成 `ConfigurationError`。
- **API：** `varcode_predict_variant_effect_on_transcript(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, annotator: 'str | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`annotator` 为 None 或 `fast`、`protein_diff`、`transcript_model`、`realized` 之一。
- **输出：** `variant`（varcode 规范化后的形式）、`transcript_id` 和 `effect`（`class`、`description`、`variant`、`transcript_id`、`gene_name`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_variant_effect_on_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_predict_variant_effect_on_transcript(variant={"contig": "1", "start": 20, "ref": "TC", "alt": "T"}, transcript_id="T1", **ref)
print(result["variant"])
print(result["effect"]["class"], result["effect"]["description"])
```

- **示例输出：**

```text
{'contig': '1', 'start': 21, 'ref': 'C', 'alt': ''}
FrameShift p.Y5fs
```

## 16) `dnakit.molecular_lab.varcode_predict_variant_effect_on_transcript_or_failure()`

- **作用：** 预测一个变异对一个转录本的效应，失败时报告失败而不是抛出异常。
- **处理方法：** 调用 `varcode.effects.predict_variant_effect_on_transcript_or_failure`。当 varcode 返回异常，或返回 `Failure`、`Unresolved` 效应（例如 `ref` 等位基因与参考不符）时，结果中 `status` 为 `failed`、`effect` 为 None，失败文本放在 `failure` 中；否则 `status` 为 `completed`。无效参数仍然抛出 `ConfigurationError`。
- **API：** `varcode_predict_variant_effect_on_transcript_or_failure(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `variant`、`transcript_id`、`status`（`completed` 或 `failed`）、`effect`（失败时为 None）和 `failure`（上游失败文本或 None）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_variant_effect_on_transcript_or_failure

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

good = varcode_predict_variant_effect_on_transcript_or_failure(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(good["status"], good["effect"]["description"], good["failure"])

wrong_ref = {"contig": "1", "start": 13, "ref": "A", "alt": "G"}
bad = varcode_predict_variant_effect_on_transcript_or_failure(variant=wrong_ref, transcript_id="T1", **ref)
print(bad["status"], bad["effect"], "Reference allele mismatch" in bad["failure"])
```

- **示例输出：**

```text
completed p.P2A None
failed None True
```

## 17) `dnakit.molecular_lab.varcode_predict_variant_effects()`

- **作用：** 预测一个变异对其重叠的每个转录本的效应。
- **处理方法：** 在本地参考上构建 varcode 的 `Variant`，调用 `varcode.effects.predict_variant_effects`。不与任何基因重叠的变异得到 `Intergenic` 效应。默认 `raise_on_error=False` 时，注释错误（例如 `ref` 等位基因错误）以 `Failure` 效应返回；`raise_on_error=True` 时则抛出 `ConfigurationError`。
- **API：** `varcode_predict_variant_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', *, raise_on_error: 'bool' = False, annotator: 'str | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`raise_on_error` 是布尔值（默认 False）；`annotator` 为 None 或 `fast`、`protein_diff`、`transcript_model`、`realized` 之一。
- **输出：** `variant`、`effects`（每个效应一条记录，含 `class`、`description`，有转录本时还含 `transcript_id` 和 `gene_name`）和 `count`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_variant_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = {
    "missense": {"contig": "1", "start": 13, "ref": "C", "alt": "G"},
    "splice donor": {"contig": "1", "start": 61, "ref": "G", "alt": "A"},
    "between genes": {"contig": "1", "start": 255, "ref": "G", "alt": "A"},
}
for label, variant in variants.items():
    result = varcode_predict_variant_effects(variant=variant, **ref)
    print(label, result["count"], [(e["class"], e["description"]) for e in result["effects"]])
```

- **示例输出：**

```text
missense 1 [('Substitution', 'p.P2A')]
splice donor 1 [('SpliceOutcomeSet', 'splice-set:intron-retention (donor side, predicted)')]
between genes 1 [('Intergenic', 'intergenic')]
```

## 18) `dnakit.molecular_lab.varcode_apply_phase_resolver_to_effects()`

- **作用：** 预测变异的效应，并用 VCF 样本中记录的相位对效应进行解析。
- **处理方法：** 在本地参考上用 `varcode.load_vcf` 加载 `phased_vcf`，为该样本构建 `VCFPhaseResolver`，预测变异的效应，再交给 `varcode.apply_phase_resolver_to_effects`。在 varcode 10.11 中 `VCFPhaseResolver` 只回答顺式与反式问题，没有突变转录本方法，因此 `apply_phase_resolver_to_effects` 不会改变这些效应：结果与普通的效应预测相同，调用主要用于检查带相位的 VCF 与样本能否加载。
- **API：** `varcode_apply_phase_resolver_to_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', phased_vcf: 'str', sample: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`phased_vcf` 是本地 VCF，其 `sample` 列带有定相基因型（`0|1`）和 `PS` 值；`sample` 是列名。VCF 中的染色体名在 varcode 去掉 `chr` 前缀后必须与参考一致。
- **输出：** `variant`、`sample`、`effects`（含 `class`、`description`、`transcript_id`、`gene_name` 的记录）和 `count`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_apply_phase_resolver_to_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_apply_phase_resolver_to_effects(
    variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, phased_vcf="phased.vcf", sample="S1", **ref
)
print(result["sample"], result["count"], [(e["class"], e["description"]) for e in result["effects"]])
```

- **示例输出：**

```text
S1 1 [('Substitution', 'p.P2A')]
```

## 19) `dnakit.molecular_lab.varcode_build_exon_runs()`

- **作用：** 把转录本的有义链外显子投影到它的基因组布局上。
- **处理方法：** 为转录本构建两侧各带 `flank` 个碱基的 `GenomicLayout`（`GenomicLayout.from_transcript`，序列取自本地基因组 FASTA），并调用 `varcode.transcript_layout.build_exon_runs`。每个 run 是参考外显子在布局中的一次出现。
- **API：** `varcode_build_exon_runs(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, flank: 'int' = 50, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`flank` 是 0 到 100000 的整数（默认 50）；侧翼在染色体起点处被截断。
- **输出：** `transcript_id`、`flank` 和 `runs`。每个 run 含 `key`（外显子编号与出现次数）、`exon_number`、`occurrence`、`layout_start` 与 `layout_end`（布局中的 0 起始偏移）、`reference_start`、`reference_end`、`first_genomic_position` 和 `last_genomic_position`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_build_exon_runs

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_build_exon_runs(transcript_id="T2", flank=10, **ref)
for run in result["runs"]:
    print(run["key"], run["layout_start"], run["layout_end"], run["first_genomic_position"], run["last_genomic_position"])
```

- **示例输出：**

```text
[1, 1] 10 60 281 330
[2, 1] 100 150 371 420
```

## 20) `dnakit.molecular_lab.varcode_enumerate_from_structural_variant()`

- **作用：** 枚举结构变异断点两侧的隐蔽外显子候选。
- **处理方法：** 在本地参考上构建 varcode 的 `StructuralVariant`，用侧翼大小调用 `varcode.cryptic_exons.enumerate_from_structural_variant`。varcode 用默认的剪接位点打分器扫描断点周围序列，把类受体与类供体基序配成候选外显子。包装函数对每个候选调用 `str()`，因此候选只以类名和变异的文本形式报告；varcode 对象中的区间以及供体、受体分数不在结果里。当某个断点没有可用的参考序列时，varcode 会对每个基因组打印一次警告。
- **API：** `varcode_enumerate_from_structural_variant(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', *, flank: 'int' = 500, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structural_variant` 含 `contig`、`start`、`end`（`end` 不小于 `start`）和 `sv_type`（`DEL`、`DUP`、`INV`、`INS` 或 `CNV`，不区分大小写）；`flank` 是 0 到 100000 的整数（默认 500）。
- **输出：** `structural_variant`、`flank`、`candidates`（字符串）和 `count`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_enumerate_from_structural_variant

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

sv = {"contig": "1", "start": 70, "end": 90, "sv_type": "DEL"}
narrow = varcode_enumerate_from_structural_variant(structural_variant=sv, flank=30, **ref)
wide = varcode_enumerate_from_structural_variant(structural_variant=sv, flank=100, **ref)
print(narrow["count"], wide["count"])
print(wide["candidates"][0])
```

- **示例输出：**

```text
0 1
CrypticExonCandidate(variant=StructuralVariant(contig='1', start=70, end=90, sv_type='DEL', alt='<DEL>', reference_name='dnakit_local'))
```

## 21) `dnakit.molecular_lab.varcode_enumerate_layout_splice_plans()`

- **作用：** 枚举转录本布局中满足剪接图约束的剪接方案。
- **处理方法：** 构建转录本的基因组布局，调用 `varcode.splice_graph.enumerate_layout_splice_plans`。参考布局中每个不完整的供体或受体（此处根据布局序列判断）给出一个剪接选择，机制包括外显子跳跃、内含子保留、隐蔽位点和正常剪接；每个方案对每个受损位点固定一个选择，并列出保留的外显子 run。布局完整时只有一个方案。
- **API：** `varcode_enumerate_layout_splice_plans(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, flank: 'int' = 50, max_plans: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`flank` 是 0 到 100000 的整数（默认 50）；`max_plans` 是 1 到 10000 的整数（默认 64）。
- **输出：** `transcript_id`、`flank`、`plans` 和 `count`。每个方案含 `choices`（位点与机制的配对，机制含 `ordinal_rank`、`skipped_runs` 和 `required_runs`）以及 `kept_runs`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_enumerate_layout_splice_plans

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

broken = varcode_enumerate_layout_splice_plans(transcript_id="T1", flank=10, **ref)
print(broken["count"])
for plan in broken["plans"]:
    site, choice = plan["choices"][0]
    print(site["side"], site["exon_number"], choice["mechanism"], plan["kept_runs"])

intact = varcode_enumerate_layout_splice_plans(transcript_id="T2", flank=10, **ref)
print(intact["count"], intact["plans"][0]["choices"])
```

- **示例输出：**

```text
4
donor 2 exon_skip [[1, 1], [3, 1]]
donor 2 intron_retention [[1, 1], [2, 1], [3, 1]]
donor 2 cryptic [[1, 1], [2, 1], [3, 1]]
donor 2 normal [[1, 1], [2, 1], [3, 1]]
1 []
```

## 22) `dnakit.molecular_lab.varcode_predict_structural_variant_effect()`

- **作用：** 对某个转录本上的结构变异（缺失、重复、倒位、插入或拷贝数变异）分类。
- **处理方法：** 在本地参考上构建 varcode 的 `StructuralVariant`，用转录本调用 `varcode.effects.structural.predict_structural_variant_effect`。
- **API：** `varcode_predict_structural_variant_effect(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structural_variant` 含 `contig`、`start`、`end`（`end` 不小于 `start`）和 `sv_type`（`DEL`、`DUP`、`INV`、`INS` 或 `CNV`，不区分大小写）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `structural_variant`（解析后的形式）、`transcript_id` 和 `effect`（`class`、`description`、`variant`、`transcript_id`、`gene_name`）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_predict_structural_variant_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

cases = [
    ("T1", {"contig": "1", "start": 101, "end": 150, "sv_type": "DEL"}),
    ("T1", {"contig": "1", "start": 101, "end": 150, "sv_type": "DUP"}),
    ("T1", {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}),
]
for transcript_id, sv in cases:
    result = varcode_predict_structural_variant_effect(structural_variant=sv, transcript_id=transcript_id, **ref)
    print(sv["sv_type"], sv["start"], sv["end"], result["effect"]["class"], result["effect"]["description"])
```

- **示例输出：**

```text
DEL 101 150 FrameShift p.G18fs
DUP 101 150 FrameShift p.L35fs
DEL 151 280 GeneFusion sv-gene-fusion
```

## 23) `dnakit.molecular_lab.varcode_random_variants()`

- **作用：** 在本地参考的完整编码转录本上随机抽取变异。
- **处理方法：** 以本地基因组作为 `ensembl` 调用 `varcode.util.random_variants`。它在完整编码转录本（有起始密码子、终止密码子和蛋白的转录本）中抽取位置，并根据参考序列构造替换变异，启用时还构造插入和缺失。位置和参考碱基由 `random_seed` 决定；但替代碱基由 varcode 从 Python `set` 中抽取，而 `set` 的顺序在不同 Python 进程之间会变化（字符串哈希是随机化的），所以相同的种子在不同进程中可能给出不同的 `alt`；同一进程内的抽取是可重复的。
- **API：** `varcode_random_variants(gtf: 'str', genome_fasta: 'str', count: 'int', *, deletions: 'bool' = True, insertions: 'bool' = True, random_seed: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`count` 是 1 到 10000 的整数；`deletions` 与 `insertions` 是布尔值（默认 True）；`random_seed` 为 None 或 0 到 2**32 - 1 的整数。
- **输出：** `count`、`random_seed` 和 `variants`（含 `contig`、`start`、`ref`、`alt` 的记录）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_random_variants

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

first = varcode_random_variants(count=4, deletions=False, insertions=False, random_seed=5, **ref)
second = varcode_random_variants(count=4, deletions=False, insertions=False, random_seed=5, **ref)
print(first["count"], first["variants"] == second["variants"])
for variant in first["variants"]:
    print(variant["contig"], variant["start"], variant["ref"], len(variant["alt"]), variant["alt"] != variant["ref"])
```

- **示例输出：**

```text
4 True
1 8 C 1 True
1 150 A 1 True
1 305 C 1 True
1 415 T 1 True
```

## 24) `dnakit.molecular_lab.varcode_realize_exon_path()`

- **作用：** 拼接转录本中选定的外显子 run，并从映射到的 CDS 起点开始翻译。
- **处理方法：** 构建基因组布局并调用 `varcode.transcript_layout.realize_exon_path`，它把选中的外显子 run 连接成 cDNA，并从 CDS 起点翻译。`kept_run_keys` 为 None 时保留所有 run。注意：包装函数只接受非空字符串形式的 `kept_run_keys`，而 `varcode_build_exon_runs()` 报告的每个键是一对整数（例如 `[1, 1]`），varcode 把键当作元组比较；字符串永远不匹配，所以任何字符串列表都选不到 run，得到空的产物。只有默认的 `None` 才能得到可用的产物。
- **API：** `varcode_realize_exon_path(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, kept_run_keys: 'Sequence[str] | None' = None, flank: 'int' = 50, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`kept_run_keys` 为 None 或非空字符串列表（见上面的说明）；`flank` 是 0 到 100000 的整数（默认 50）。
- **输出：** `transcript_id` 与 `product`（含 `transcript`（文本）、`cdna_sequence`、`protein_sequence`、`junction_signature`、`start_codon_present` 和 `evidence`（`start_codon`、`exon_runs`、`genomic_positions`））。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_realize_exon_path

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

everything = varcode_realize_exon_path(transcript_id="T1", flank=10, **ref)
product = everything["product"]
print(len(product["cdna_sequence"]), len(product["protein_sequence"]), product["start_codon_present"])
print(product["evidence"]["exon_runs"], len(product["evidence"]["genomic_positions"]))

by_string = varcode_realize_exon_path(transcript_id="T1", kept_run_keys=["1"], flank=10, **ref)
print(len(by_string["product"]["cdna_sequence"]), repr(by_string["product"]["protein_sequence"]))
```

- **示例输出：**

```text
160 46 True
[[1, 1], [2, 1], [3, 1]] 160
0 ''
```

## 25) `dnakit.molecular_lab.varcode_realize_splice_plan()`

- **作用：** 在转录本布局上实现一个已枚举的剪接方案（外显子跳跃、内含子保留或隐蔽位点选择）。
- **处理方法：** 构建布局，像 `varcode_enumerate_layout_splice_plans()` 一样枚举方案，用 `varcode.splice_graph.disrupted_splice_sites` 计算受损的剪接位点，再对 `plan_index` 指定的方案调用 `varcode.transcript_layout.realize_splice_plan`。产物包含剪接后的 cDNA 和从 CDS 起点翻译的蛋白；对于隐蔽位点选择，会使用 `scan_flank` 个碱基范围内最佳的非标准边界。
- **API：** `varcode_realize_splice_plan(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, plan_index: 'int' = 0, flank: 'int' = 50, scan_flank: 'int' = 50, max_plans: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`plan_index` 是 0 到 `max_plans` 的整数，且必须小于方案数（否则抛出 `ConfigurationError`）；`flank` 与 `scan_flank` 是整数（`scan_flank` 为 0 到 10000，默认 50）；`max_plans` 是 1 到 10000 的整数。
- **输出：** `transcript_id`、`plan_index`、`plan_count`、`product`（与 `varcode_realize_exon_path()` 相同，其 `evidence` 中另有 `splice_choices` 和 `resolved_splice_choices`）以及 `unresolved_mechanisms`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_realize_splice_plan

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

from dnakit.exceptions import ConfigurationError

for index in range(4):
    result = varcode_realize_splice_plan(transcript_id="T1", plan_index=index, flank=10, **ref)
    product = result["product"]
    mechanism = product["evidence"]["splice_choices"][0][1]
    print(index, mechanism, len(product["cdna_sequence"]), len(product["protein_sequence"]), result["plan_count"])

try:
    varcode_realize_splice_plan(transcript_id="T1", plan_index=9, flank=10, **ref)
except ConfigurationError as error:
    print(type(error).__name__, "out of range" in str(error))
```

- **示例输出：**

```text
0 exon_skip 110 33 4
1 intron_retention 200 63 4
2 cryptic 134 41 4
3 normal 160 46 4
ConfigurationError True
```

## 26) `dnakit.molecular_lab.varcode_reference_base()`

- **作用：** 返回 1 起始位置上的正链参考碱基。
- **处理方法：** 在本地参考上调用 `varcode.reference_base(genome, contig, position)`，它会读取基因组 FASTA。
- **API：** `varcode_reference_base(gtf: 'str', genome_fasta: 'str', contig: 'str', position: 'int', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`contig` 是参考中存在的普通标识符；`position` 是不小于 1 的整数。位置超出染色体末端，或染色体不在参考中时，`base` 为空字符串，不会报错。
- **输出：** `contig`、`position` 和 `base`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_reference_base

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_base(contig="1", position=13, **ref)
print(result["contig"], result["position"], result["base"])
beyond = varcode_reference_base(contig="1", position=9999, **ref)
print(repr(beyond["base"]))
```

- **示例输出：**

```text
1 13 C
''
```

## 27) `dnakit.molecular_lab.varcode_reference_range()`

- **作用：** 返回 1 起始、含两端区间上的正链参考序列。
- **处理方法：** 在本地参考上调用 `varcode.reference_range(genome, contig, start, end)`，它会读取基因组 FASTA。即使基因位于负链，返回的序列也总是正链序列。
- **API：** `varcode_reference_range(gtf: 'str', genome_fasta: 'str', contig: 'str', start: 'int', end: 'int', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`contig` 是参考中存在的普通标识符；`start` 是不小于 1 的整数，`end` 是不小于 `start` 的整数。
- **输出：** `contig`、`start`、`end` 和 `sequence`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_reference_range

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_range(contig="1", start=10, end=21, **ref)
print(result["sequence"], len(result["sequence"]))
```

- **示例输出：**

```text
ATGCCAACAGTC 12
```

## 28) `dnakit.molecular_lab.varcode_reference_selenocysteine()`

- **作用：** 返回转录本中已注释硒代半胱氨酸密码子的 cDNA 偏移。
- **处理方法：** 调用 `varcode.effects.selenocysteine.reference_selenocysteine`，它在完整转录本的蛋白序列中查找 `U`，并返回每个这类密码子第一个碱基的 cDNA 偏移。在本地参考中，蛋白 FASTA 由 DNAKit 自己的翻译从 CDS 行得到，该翻译把 `TGA` 读作终止，从不产生 `U`，所以对本地参考，这个函数对每个转录本都返回空列表。
- **API：** `varcode_reference_selenocysteine(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。
- **输出：** `transcript_id` 和 `cdna_offsets`（整数列表）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_reference_selenocysteine

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

for transcript_id in ("T1", "T2"):
    result = varcode_reference_selenocysteine(transcript_id=transcript_id, **ref)
    print(transcript_id, result["cdna_offsets"])
```

- **示例输出：**

```text
T1 []
T2 []
```

## 29) `dnakit.molecular_lab.varcode_make_fusion_outcome()`

- **作用：** 导入结构变异上一条观察到的 RNA 连接或融合模型。
- **处理方法：** 在本地参考上构建 varcode 的 `StructuralVariant`，用 5' 锚定转录本、观察到的序列、模型 ID、可选的伙伴转录本、可选的 `cds_start` 和读段数调用 `varcode.rna_evidence.make_fusion_outcome`。给出 `partner_transcript_id` 时，结果按有义方向连接两个转录本；否则是易位到基因间区。给出 `cds_start` 时，varcode 翻译序列中从起始到终止的完整 ORF（这是从 RNA 得到的预测，不是蛋白观察）；不给时蛋白未确定。
- **API：** `varcode_make_fusion_outcome(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', transcript_id: 'str', sequence: 'str', transcript_model_id: 'str', *, partner_transcript_id: 'str | None' = None, cds_start: 'int | None' = None, source: 'str' = 'rna', read_count: 'int | None' = None, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structural_variant` 含 `contig`、`start`、`end` 和 `sv_type`。`transcript_id` 是已注释的 5' 锚定转录本；`sequence` 是观察到的 5' 到 3' 序列（A、C、G、T 或 N，最长 100000 个碱基）；`transcript_model_id` 是标识符；`partner_transcript_id` 为 None 或转录本 ID；`cds_start` 为 None 或完整 ORF 起点的 0 起始偏移（否则抛出 `ConfigurationError`）；`source` 是标识符（默认 `rna`）；`read_count` 为 None 或非负整数。
- **输出：** `structural_variant`、`transcript_id`、`partner_transcript_id`、`transcript_model_id`、`candidate`（`source`、`evidence` 和 `effect`）以及 `mutant_transcript`（`cdna_sequence`、`mutant_protein_sequence`、`edits` 等）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_make_fusion_outcome

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}
import csv

with open("structures.tsv", newline="") as handle:
    sequence = "".join(row["sequence"] for row in csv.DictReader(handle, delimiter="\t"))

sv = {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}
plain = varcode_make_fusion_outcome(
    structural_variant=sv, transcript_id="T1", sequence=sequence, transcript_model_id="M1",
    partner_transcript_id="T2", read_count=7, **ref
)
print(plain["candidate"]["effect"]["class"], plain["candidate"]["evidence"]["protein_status"], plain["mutant_transcript"]["mutant_protein_sequence"])

orf = varcode_make_fusion_outcome(
    structural_variant=sv, transcript_id="T1", sequence=sequence, transcript_model_id="M1",
    partner_transcript_id="T2", cds_start=22, read_count=7, **ref
)
print(orf["candidate"]["evidence"]["protein_status"], orf["mutant_transcript"]["mutant_protein_sequence"])
print(len(sequence), orf["candidate"]["source"], orf["candidate"]["evidence"]["read_count"])
```

- **示例输出：**

```text
GeneFusion not_determined None
predicted_from_observed_rna MELYITIAYV
160 rna 7
```

## 30) `dnakit.molecular_lab.varcode_apply_rna_evidence_to_effects()`

- **作用：** 把观察到的融合结果挂到结构变异的预测效应上。
- **处理方法：** 对每个条目像 `varcode_make_fusion_outcome()` 一样构建融合结果，把它们放进 varcode 的 `RNAEvidence` 解析器，在本地参考上预测结构变异的效应，再调用 `varcode.apply_rna_evidence_to_effects`。DNA 预测的效应被保留，每个效应列出它新增的候选（varcode 自己的备选加上观察到的 RNA 候选）。
- **API：** `varcode_apply_rna_evidence_to_effects(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', fusion_outcomes: 'list[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structural_variant` 含 `contig`、`start`、`end` 和 `sv_type`。`fusion_outcomes` 是 1 到 200 个映射的列表，含 `transcript_id`、`transcript_model_id` 和 `sequence`，可选 `partner_transcript_id`、`cds_start`、`read_count` 和 `source`。
- **输出：** `structural_variant`、`fusion_count`、`effects`（每个含 `class`、`description`、`transcript_id`、`gene_name` 和 `candidates`）、`count` 以及 `observed_candidate_count`（所有效应的全部候选数，包括 varcode 自己的）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_apply_rna_evidence_to_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}
import csv

with open("structures.tsv", newline="") as handle:
    sequence = "".join(row["sequence"] for row in csv.DictReader(handle, delimiter="\t"))

fusion = {
    "transcript_id": "T1", "transcript_model_id": "M1", "sequence": sequence,
    "partner_transcript_id": "T2", "read_count": 7,
}
result = varcode_apply_rna_evidence_to_effects(structural_variant={"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}, fusion_outcomes=[fusion], **ref)
print(result["fusion_count"], result["count"], result["observed_candidate_count"])
effect = result["effects"][0]
print(effect["class"], effect["transcript_id"])
print([(c["source"], c["effect"]["class"]) for c in effect["candidates"]])
```

- **示例输出：**

```text
1 1 6
GeneFusion T1
[('varcode', 'GeneFusion'), ('varcode', 'StructuralVariantEffect'), ('varcode_splice', 'ExonSkipping'), ('varcode_splice', 'IntronRetention'), ('varcode_splice', 'CrypticDonor'), ('rna', 'GeneFusion')]
```

## 31) `dnakit.molecular_lab.varcode_load_exacto_fusions()`

- **作用：** 导入与结构变异关联的 ExACTO 转录本结构表和整合表。
- **处理方法：** 为每个 DNA 检出 ID 构建 varcode 的 `StructuralVariant`，再用两个 TSV 文件调用 `varcode.exacto.load_exacto_fusions`。只导入 `dna_variant_call_id` 是 `structural_variants` 键的整合行。序列按转录本方向读取。线性的双位点模型成为融合结果，与参考一致的碱基映射到已注释的转录本；环状、多基因和不完整的路径会报错。原生的 primary-structure 肽段和显式 ORF 起点没有暴露，因此蛋白未确定。
- **API：** `varcode_load_exacto_fusions(gtf: 'str', genome_fasta: 'str', structures_path: 'str', integrated_path: 'str', structural_variants: 'Mapping[str, Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structures_path` 与 `integrated_path` 是已存在的制表符分隔文件（接受 gzip）；结构表需要列 `transcript_model_id`、`reference_transcript_ids`、`index`、`read_start`、`read_end`、`sequence`、`type`、`kind`、`context`、`chromosome_1`、`position_1`、`strand_1`、`chromosome_2`、`position_2`、`strand_2`、`transcript_id_1` 和 `transcript_id_2`，整合表需要列 `transcript_model_id`、`reference_transcript_ids`、`rna_variant_call_id` 和 `dna_variant_call_id`。`structural_variants` 把 1 到 200 个 DNA 检出 ID 映射到结构变异映射。
- **输出：** `structural_variants`、`candidates`（记录，其 `source` 为 `exacto`，`evidence` 含 DNA 与 RNA 检出 ID、参考转录本、`partner_status` 和原始表格行，另有 `effect`）、`count`，以及 `structures_path` 与 `integrated_path`（解析后的路径）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_load_exacto_fusions

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_load_exacto_fusions(
    structures_path="structures.tsv", integrated_path="integrated.tsv",
    structural_variants={"D1": {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}}, **ref
)
print(result["count"], list(result["structural_variants"]))
candidate = result["candidates"][0]
evidence = candidate["evidence"]
print(candidate["source"], candidate["effect"]["class"], candidate["effect"]["transcript_id"])
print(evidence["dna_variant_call_id"], evidence["rna_variant_call_ids"], evidence["reference_transcript_ids"], evidence["partner_status"])
print(len(evidence["exacto_structure"]), evidence["protein_status"])
```

- **示例输出：**

```text
1 ['D1']
exacto GeneFusion T1
D1 ['R1'] ['T1', 'T2'] sense
3 not_determined
```

## 32) `dnakit.molecular_lab.varcode_reference_completion_hypotheses()`

- **作用：** 把 ExACTO 融合中已映射的末端补全为明确的参考假设。
- **处理方法：** 先像 `varcode_load_exacto_fusions()` 一样加载观察到的融合，再对每个融合调用 `varcode.reference_completion.reference_completion_hypotheses`。观察模型中未映射的 5' 和/或 3' 末端用已注释转录本的参考碱基延伸，每种转录本选择对应一个假设。这些假设被标为参考补全假设，绝不当作 RNA 观察，读段数也不被用作支持。
- **API：** `varcode_reference_completion_hypotheses(gtf: 'str', genome_fasta: 'str', structures_path: 'str', integrated_path: 'str', structural_variants: 'Mapping[str, Mapping[str, Any]]', *, ends: 'list[str] | None' = None, five_prime_transcript_ids: 'list[str] | None' = None, three_prime_transcript_ids: 'list[str] | None' = None, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`structures_path`、`integrated_path` 和 `structural_variants` 与 `varcode_load_exacto_fusions()` 相同。`ends` 为 None（两端都补全）或由 `five_prime` 和/或 `three_prime` 组成的列表；`five_prime_transcript_ids` 与 `three_prime_transcript_ids` 为 None（使用该末端已注释的转录本）或最多 500 个参考中存在的转录本 ID 的列表。
- **输出：** `structural_variants`、`ends`、`observed_count`、`hypothesis_count` 和 `hypotheses`。每个假设含 `source`（`varcode_reference_completion`）、`evidence`（包括 `observed_evidence`、含观察区间与假设区间的 `reference_completion`、`protein_status`、`protein_completeness`、`cds_start` 和 `cds_end`）以及 `effect`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_reference_completion_hypotheses

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_completion_hypotheses(
    structures_path="structures.tsv", integrated_path="integrated.tsv",
    structural_variants={"D1": {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}}, ends=["three_prime"], **ref
)
print(result["ends"], result["observed_count"], result["hypothesis_count"])
hypothesis = result["hypotheses"][0]
evidence = hypothesis["evidence"]
print(hypothesis["source"], evidence["sequence_status"], evidence["protein_status"], evidence["protein_completeness"])
completion = evidence["reference_completion"]
print(completion["observed_span"], [(s["transcript_id"], s["start"], s["end"], s["label"]) for s in completion["assumed_spans"]])
```

- **示例输出：**

```text
['three_prime'] 1 1
varcode_reference_completion reference_completion_hypothesis predicted_from_reference_completion start_to_stop
[0, 160] [('T2', 160, 210, 'assumed_reference_three_prime')]
```

## 33) `dnakit.molecular_lab.varcode_enumerate_splice_outcomes()`

- **作用：** 枚举某个转录本上剪接信号被破坏背后的剪接机制。
- **处理方法：** 先预测变异对转录本的效应，再调用 `varcode.splice_outcomes.enumerate_splice_outcomes`。剪接供体、受体或剪接位点效应会被包装成 `SpliceOutcomeSet`，其候选就是各种机制（例如内含子保留、隐蔽位点和正常剪接）。`use_genomic_sequence` 为真时，包装函数提供本地参考的基因组序列，因此会扫描隐蔽位点并实际构造内含子保留（描述中会带上蛋白后果）。其他效应原样返回，没有候选。单外显子转录本没有剪接连接。
- **API：** `varcode_enumerate_splice_outcomes(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, use_genomic_sequence: 'bool' = True, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variant` 是含 `contig`、`start`（1 起始整数）、`ref` 和 `alt`（A、C、G、T 或 N，不能同时为空）的映射。`transcript_id` 是普通标识符（字母、数字和 `._:+|-`），必须存在于本地参考中。`use_genomic_sequence` 是布尔值（默认 True）。
- **输出：** `variant`、`transcript_id`、`input_effect`、`outcome_class`、`is_outcome_set`、`genomic_sequence_used`、`candidates`（含 `source`、`evidence` 和 `effect` 的记录）和 `count`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_enumerate_splice_outcomes

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

with_sequence = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 61, "ref": "G", "alt": "A"}, transcript_id="T1", **ref)
print(with_sequence["outcome_class"], with_sequence["is_outcome_set"], with_sequence["genomic_sequence_used"], with_sequence["count"])
for candidate in with_sequence["candidates"]:
    print(candidate["effect"]["class"], "|", candidate["effect"]["description"])

without = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 61, "ref": "G", "alt": "A"}, transcript_id="T1", use_genomic_sequence=False, **ref)
print(without["genomic_sequence_used"], [c["effect"]["description"] for c in without["candidates"]])

missense = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(missense["outcome_class"], missense["is_outcome_set"], missense["count"])
```

- **示例输出：**

```text
SpliceOutcomeSet True True 3
IntronRetention | intron-retention (donor side, p.G18fs*)
CrypticDonor | cryptic-donor ext8bp
NormalSplicing | normal-splicing
False ['intron-retention (donor side, predicted)', 'cryptic-donor (predicted)', 'normal-splicing']
Substitution False 0
```

## 34) `dnakit.molecular_lab.varcode_install_reference_data()`

- **作用：** 为每个参考名称下载并索引 varcode 所需的 Ensembl 参考数据。
- **处理方法：** 对每个不同的名称，pyensembl 通过 `pyensembl.genome_for_reference_name` 选择发布版本：`allow_older_downloaded_release` 为真时选择已下载的最新发布，否则选择最新发布。如果该发布的 GTF 和 FASTA 文件缺失，就从 Ensembl 下载（`genome.download()`）并建立索引（`genome.index()`），与 varcode 命令行的做法相同；已存在的文件不会重复下载。名称在调用 pyensembl 之前先做检查，未知名称会抛出 `ConfigurationError`。
- **API：** `varcode_install_reference_data(reference_names: 'Sequence[str]', *, allow_older_downloaded_release: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（额外依赖 `dnakit[bio-pyensembl]`）；当某个发布的 Ensembl 文件尚不存在时，还需要能访问 Ensembl FTP 服务器的网络，以及人类注释所需的磁盘空间（数百兆字节）。下载失败会抛出 `DownloadError`，建立索引失败会抛出 `BackendExecutionError`。下面的示例不会下载任何内容：它展示未知参考名称的报错，该错误发生在任何下载之前。
- **输入：** `reference_names` 是 1 到 16 个名称的列表，每个名称最长 64 个字符，由字母、数字、`_`、`.` 和 `-` 组成（例如 `GRCh38`）；重复项会被忽略。`allow_older_downloaded_release` 是布尔值（默认 True）。
- **输出：** `references`（每个名称一条记录，含 `reference_name`、`release`、`species`、`cache_directory`、`files_present_before`、`downloaded`、`indexed`）和 `allow_older_downloaded_release`。
- **安装依赖后的用法：**

```python
from dnakit.molecular_lab import varcode_install_reference_data

result = varcode_install_reference_data(["GRCh38"])
for reference in result["references"]:
    print(reference["reference_name"], reference["release"], reference["downloaded"], reference["indexed"])
```

- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import varcode_install_reference_data

for names in (["NOT_A_GENOME"], [], ["bad name"]):
    try:
        varcode_install_reference_data(names, verbose=False)
    except ConfigurationError as error:
        message = str(error)
        print(type(error).__name__, "|", message.split("] ", 1)[-1].split(":")[0])
```

- **示例输出：**

```text
ConfigurationError | Unknown reference name 'NOT_A_GENOME'
ConfigurationError | reference_names must hold 1 to 16 names.
ConfigurationError | Invalid reference name 'bad name'.
```

## 35) `dnakit.molecular_lab.varcode_variant_gene_overlaps()`

- **作用：** 利用本地 GTF 和基因组 FASTA，列出与每个变异重叠的基因。
- **处理方法：** 构建（或复用）本地参考，对每个变异创建 varcode 的 `Variant` 并读取其 `genes` 属性。varcode 无法注释的变异（例如染色体不在 GTF 中）会附上错误信息，`raise_on_error` 为真时则抛出 `BackendExecutionError`。
- **API：** `varcode_variant_gene_overlaps(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, raise_on_error: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `gtf` 与 `genome_fasta` 是已存在的本地文件（注释和基因组 FASTA），`cache_dir` 是索引目录（默认 `~/.cache/dnakit/varcode_reference`）。`variants` 是 1 到 10000 个映射的列表，含 `contig`（GTF 中的名称）、`position`（1 起始整数）、`ref` 和 `alt`（由 1 到 256 个 A、C、G、T 或 N 组成的 DNA 等位基因；与其他 varcode 函数不同，这里不允许空等位基因）。`raise_on_error` 是布尔值（默认 False）。
- **输出：** `variants`（每个输入一条记录，含 `variant_index`、`variant`、`contig`、`position`、`ref`、`alt`、`genes` 和 `error`）、`variant_count`、`annotation_failures` 和 `reference`（摘要、缓存目录以及索引是否被复用）。每个基因含 `gene_id`、`gene_name`、`biotype`、`contig`、`start`、`end` 和 `strand`。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_variant_gene_overlaps

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = [
    {"contig": "1", "position": 13, "ref": "C", "alt": "G"},
    {"contig": "1", "position": 300, "ref": "G", "alt": "A"},
    {"contig": "1", "position": 255, "ref": "G", "alt": "A"},
    {"contig": "2", "position": 5, "ref": "C", "alt": "G"},
]
result = varcode_variant_gene_overlaps(variants=variants, **ref)
print(result["variant_count"], result["annotation_failures"])
for row in result["variants"]:
    print(row["variant"], [g["gene_name"] for g in row["genes"]], row["error"] is not None)
gene = result["variants"][0]["genes"][0]
print(gene["gene_id"], gene["biotype"], gene["start"], gene["end"], gene["strand"])
```

- **示例输出：**

```text
4 1
1:13:C>G ['GENEA'] False
1:300:G>A ['GENEB'] False
1:255:G>A [] False
2:5:C>G [] True
G1 protein_coding 1 240 +
```

## 36) `dnakit.molecular_lab.varcode_check_sample_identity_files()`

- **作用：** 检查多个本地 VCF 文件中的样本是否来自同一供体。
- **处理方法：** 用 `varcode.sample_identity.load_vcf_samples` 加载每个选中的 VCF 样本（作为种系、体细胞或混合样本，以及可选清单里的样本），由 `config` 构建 `SampleCheckConfig`，再调用 `varcode.sample_identity.check_sample_identity`，按基因型一致性和共有体细胞变异比较每一对样本。完整报告以 JSON 写入 `json_output`，样本对表写入 `tsv_output`，报告同时返回。结果是筛查性的启发式判断，不是身份证明。已存在的输出文件不会被覆盖。
- **API：** `varcode_check_sample_identity_files(*, germline_vcfs: 'Sequence[str]' = (), somatic_vcfs: 'Sequence[str]' = (), mixed_vcfs: 'Sequence[str]' = (), manifest_path: 'str | None' = None, assembly: 'str | None' = None, config: 'Mapping[str, Any] | None' = None, json_output: 'str | None' = None, tsv_output: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）。可离线运行。
- **输入：** `germline_vcfs`、`somatic_vcfs` 和 `mixed_vcfs` 是已存在的 VCF 路径列表；`manifest_path` 是可选的 CSV 或 TSV，含 `path` 和 `sample` 列（可选 `kind`、`assembly`、`role`、`label`、`donor_id`、`tumor_id`、`normal_sample`）；VCF 头部未说明时，`assembly` 为 `GRCh37` 或 `GRCh38`；`config` 是 `SampleCheckConfig` 字段的映射（未知字段会抛出 `ConfigurationError`）；`json_output` 与 `tsv_output` 是已存在目录中的新路径，彼此不同，也不同于输入。至少需要两个样本。
- **输出：** `sample_count`、`pair_count`、`identity_conflicts`、`json_output`、`tsv_output` 和 `report`（含 `samples`、`pairs` 与配置；每个样本对含 `germline` 和 `somatic` 结果，其 `status` 可为 `compatible`、`discordant` 或 `inconclusive`，另有 `flags` 与 `identity_conflict`）。报告中的样本标签包含文件路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_check_sample_identity_files

folder = Path(tempfile.mkdtemp())
result = varcode_check_sample_identity_files(
    germline_vcfs=["donor1.vcf", "donor2.vcf", "donor3.vcf"],
    assembly="GRCh38",
    config={"min_shared_snps": 5, "min_nonreference_snps": 1, "min_genomic_bins": 1},
    json_output=str(folder / "report.json"),
    tsv_output=str(folder / "pairs.tsv"),
    verbose=False,
)
print(result["sample_count"], result["pair_count"], result["identity_conflicts"])
for pair in result["report"]["pairs"]:
    germline = pair["germline"]
    print(germline["status"], germline["shared_snps"], germline["genotype_concordance"], germline["ibs0_rate"])
print(sorted(path.name for path in folder.iterdir()))
print((folder / "pairs.tsv").read_text().splitlines()[0].split("\t")[:4])
```

- **示例输出：**

```text
3 3 0
compatible 5 1.0 0.0
discordant 5 0.0 0.4
discordant 5 0.0 0.4
['pairs.tsv', 'report.json']
['sample_a', 'sample_b', 'germline_status', 'shared_snps']
```

## 37) `dnakit.molecular_lab.varcode_compare_sv_calls_directory()`

- **作用：** 对检出结果 CSV 中的结构变异断点分组，并把分组写入一个新目录。
- **处理方法：** 用 `csv.DictReader` 读取 CSV，以 `max_distance` 调用 `varcode.sv_comparison.compare_sv_calls`，并把 `groups.csv`、`members.csv` 和 `comparison.json`（含输入来源信息：varcode 版本、输入的 SHA-256 与路径、检出数和分组数）写入新目录。全链接分组把邻接拓扑相同、且每个对应断点都相距不超过 `max_distance` 个碱基的检出放在一起。组内检出的规范化等位基因不同时状态为 `nearby_candidate`；单个检出为 `exact_reported_allele`。
- **API：** `varcode_compare_sv_calls_directory(calls_csv: 'str', output_dir: 'str', *, max_distance: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）。可离线运行。
- **输入：** `calls_csv` 是已存在的 CSV，含列 `call_id`（唯一）、`caller`、`sample`、`build`（`GRCh37` 或 `GRCh38`）、`chrom`、`pos`、`ref` 和 `alt`，可选 `info`（VCF INFO，例如 `END=10500;SVTYPE=DEL`）与 `end`。`output_dir` 必须尚不存在且其父目录必须存在；`max_distance` 是非负整数（默认 100）。
- **输出：** `call_count`、`group_count`、`output_dir` 和 `files`（`groups.csv`、`members.csv`、`comparison.json`）。
- **示例代码：**

```python
import csv
import json
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_compare_sv_calls_directory

target = Path(tempfile.mkdtemp()) / "comparison"
result = varcode_compare_sv_calls_directory("calls.csv", str(target), max_distance=50, verbose=False)
print(result["call_count"], result["group_count"], result["files"])
print(sorted(path.name for path in target.iterdir()))
with open(target / "groups.csv", newline="") as handle:
    for row in csv.DictReader(handle):
        print(row["status"], row["call_ids"], row["call_count"])
print(json.loads((target / "comparison.json").read_text())["provenance"]["group_count"])
```

- **示例输出：**

```text
3 2 ['groups.csv', 'members.csv', 'comparison.json']
['comparison.json', 'groups.csv', 'members.csv']
exact_reported_allele ["c3"] 1
nearby_candidate ["c1", "c2"] 2
2
```

## 38) `dnakit.molecular_lab.varcode_effects_table()`

- **作用：** 对照本地参考预测已加载变异的编码效应，并写入 CSV 文件。
- **处理方法：** 像 `varcode_load_variant_sources()` 一样加载变异，用相同的 `gtf` 和 `genome_fasta` 把每个变异交给 `varcode_predict_variant_effects()`，每个效应写一行 CSV。设置 `skip_errors` 时，注释失败以 `Failure` 类的行保留；否则第一个失败就会抛出异常。没有任何效应的变异不产生行。varcode 命令行的 `--only-coding` 和 `--one-per-variant` 过滤器不包括在内。
- **API：** `varcode_effects_table(output_csv: 'str', *, gtf: 'str', genome_fasta: 'str', vcf_paths: 'Sequence[str]' = (), maf_paths: 'Sequence[str]' = (), json_paths: 'Sequence[str]' = (), variants: 'Sequence[Mapping[str, Any]]' = (), cache_dir: 'str | None' = None, only_passing: 'bool' = True, skip_errors: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）、pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）。可离线运行，不会下载任何数据。
- **输入：** `output_csv` 是已存在目录中的新文件路径。`gtf` 与 `genome_fasta` 是必需的关键字参数；`vcf_paths`、`maf_paths` 和 `json_paths` 是已存在文件的列表；`variants` 是变异映射（`contig`、`start`、`ref`、`alt`）的列表；至少需要一个来源。`cache_dir` 是索引目录；`only_passing`（默认 True）跳过 FILTER 不是 PASS 的 VCF 记录；`skip_errors` 是布尔值（默认 False）。
- **输出：** `reference`（本地参考的摘要）、`variant_count`、`effect_count`、`class_counts`（各效应类别及其计数）、`columns`（`contig`、`start`、`ref`、`alt`、`class`、`transcript_id`、`gene_name`、`description`）和 `output_csv`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_effects_table

target = Path(tempfile.mkdtemp()) / "effects.csv"
result = varcode_effects_table(
    str(target),
    gtf="genes.gtf",
    genome_fasta="genome.fa",
    vcf_paths=["phased.vcf"],
    variants=[{"contig": "1", "start": 61, "ref": "G", "alt": "A"}],
    cache_dir="cache",
    verbose=False,
)
print(result["variant_count"], result["effect_count"], result["class_counts"])
print(result["columns"])
print(target.read_text().strip())
```

- **示例输出：**

```text
3 3 {'SpliceOutcomeSet': 1, 'Substitution': 2}
['contig', 'start', 'ref', 'alt', 'class', 'transcript_id', 'gene_name', 'description']
contig,start,ref,alt,class,transcript_id,gene_name,description
1,13,C,G,Substitution,T1,GENEA,p.P2A
1,17,C,T,Substitution,T1,GENEA,p.T3I
1,61,G,A,SpliceOutcomeSet,T1,GENEA,"splice-set:intron-retention (donor side, predicted)"
```

## 39) `dnakit.molecular_lab.varcode_load_variant_sources()`

- **作用：** 从 VCF、MAF、JSON 文件和显式记录加载变异并取并集。
- **处理方法：** 用 `varcode.load_vcf` 加载每个 VCF（解析结构变异并应用 `only_passing`），用 `varcode.maf.load_maf` 加载每个 MAF，把显式记录作为参考上的 `varcode.Variant` 对象，用 `VariantCollection.from_json` 加载每个 JSON 文件，然后用 `VariantCollection.union` 合并。参考由 `gtf` 与 `genome_fasta` 构建（两者必须同时给出），或用 `genome` 指定一个已安装的 pyensembl 参考。VCF 的染色体名由 varcode 规范化（去掉 `chr` 前缀），而显式记录保留原样的染色体名，因此同一个变异如果在两个来源中染色体名写法不同，可能出现两次。
- **API：** `varcode_load_variant_sources(vcf_paths: 'Sequence[str]' = (), maf_paths: 'Sequence[str]' = (), json_paths: 'Sequence[str]' = (), variants: 'Sequence[Mapping[str, Any]]' = (), *, gtf: 'str | None' = None, genome_fasta: 'str | None' = None, cache_dir: 'str | None' = None, genome: 'str | None' = None, only_passing: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 varcode（额外依赖 `dnakit[bio-varcode]`）。使用 `gtf` 与 `genome_fasta` 时还需要 pyensembl（`dnakit[bio-pyensembl]`）和 pysam（`dnakit[bio-reads]`）；使用 `genome` 时需要该 Ensembl 发布已经安装。不会下载任何数据。
- **输入：** `vcf_paths`、`maf_paths` 和 `json_paths` 是已存在文件的列表；`variants` 是含 `contig`、`start`（不小于 1 的整数）、`ref` 和 `alt` 的映射列表；至少需要一个来源，显式记录还需要参考。`gtf` 与 `genome_fasta` 要么同时给出，要么都不给；`cache_dir` 是索引目录；`genome` 是 `GRCh38` 这样的参考名称；`only_passing` 是布尔值（默认 True）。
- **输出：** `reference`（`local:<摘要>`、`genome` 名称或 None）、`count`、`sources`（取自 `vcf`、`maf`、`json` 和 `variant` 的文件或记录数）和 `variants`（含 `contig`、`start`、`ref`、`alt` 的记录）。
- **示例代码：**

```python
from dnakit.molecular_lab import varcode_load_variant_sources

result = varcode_load_variant_sources(
    vcf_paths=["phased.vcf"],
    variants=[
        {"contig": "1", "start": 61, "ref": "G", "alt": "A"},
        {"contig": "1", "start": 13, "ref": "C", "alt": "G"},
    ],
    gtf="genes.gtf",
    genome_fasta="genome.fa",
    cache_dir="cache",
    verbose=False,
)
print(result["reference"].startswith("local:"), result["count"], result["sources"])
for variant in result["variants"]:
    print(variant["contig"], variant["start"], variant["ref"], variant["alt"])
```

- **示例输出：**

```text
True 3 {'vcf': 1, 'maf': 0, 'json': 0, 'variant': 2}
1 13 C G
1 17 C T
1 61 G A
```
