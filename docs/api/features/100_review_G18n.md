# G18n 审阅封装：pydna、DnaCauldron 与 DnaChisel

这组审阅封装覆盖账本中仍标为 `other` 的 pydna、DnaCauldron 与 DnaChisel 函数。每个封装位于 `src/dnakit/molecular_lab/` 下的三个模块之一：`pydna_review_functions.py`、`dnacauldron_review_functions.py` 和 `dnachisel_review_functions.py`。封装接入包的工作由负责人之后完成，因此示例直接从模块文件导入；接入后同名函数由 `dnakit.molecular_lab` 导出。

每个函数只在运行时导入对应的包；缺少包时抛出 `BackendUnavailableError` 并给出安装提示。每个函数返回可 JSON 序列化的字典，包含 `status`、`operation`、结果字段和 `provenance`。坐标从 0 开始；区间为半开区间。所有示例离线运行。未封装的函数及原因列在页面末尾。

## 1) `dnakit.molecular_lab.dna_molecules_equal()`

- **作用:** 判断 DNA 字符串是否代表同一个双链分子。
- **处理方法:** 调用 pydna 的 `eq`。序列与其反向互补序列视为相同；`circular=True` 时，环状分子的所有旋转也视为相同。上游函数：`pydna.utils.eq`。
- **API:** `dnakit.molecular_lab.dna_molecules_equal(sequences, *, circular=False, verbose=True)`
- **输入:** 至少两条 DNA 字符串（A、C、G、T；去除空白，忽略大小写）以及布尔值 `circular`。
- **输出:** `equal`、`sequence_count`、`lengths`、`circular` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import dna_molecules_equal
result = dna_molecules_equal(["ATGC", "GCAT"])
print(result["equal"], dna_molecules_equal(["ATGC", "TGCA"])["equal"])
print(dna_molecules_equal(["ATGCAAGTCC", "AAGTCCATGC"], circular=True)["equal"])
```

- **示例输出:**

```text
True False
True
```

## 2) `dnakit.molecular_lab.rare_codon_positions()`

- **作用:** 找出在某种生物密码子使用表中属于稀有的密码子位置。
- **处理方法:** 调用 pydna 的 `rarecodons`，使用内置的 `sce`（酵母）或 `eco`（大肠杆菌）稀有密码子表。从位置 0 开始按框读取密码子；末尾不足一个密码子的碱基被忽略。
- **API:** `dnakit.molecular_lab.rare_codon_positions(sequence, *, organism="sce", verbose=True)`
- **输入:** 一条 DNA 字符串，以及 `organism`（`sce` 或 `eco`）。
- **输出:** `codon_count`、`rare_codon_count`、`positions`（每项含 `start`、`end` 与 `codon`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import rare_codon_positions
result = rare_codon_positions("ATGCGATAGAGAAAAAAAAGCTCGA")
print(result["codon_count"], result["rare_codon_count"])
print(result["positions"])
```

- **示例输出:**

```text
8 1
[{'start': 3, 'end': 6, 'codon': 'CGA'}]
```

## 3) `dnakit.molecular_lab.smallest_circular_rotation()`

- **作用:** 返回环状 DNA 序列的规范形式：字典序最小的旋转。
- **处理方法:** 调用 pydna 的 `smallest_rotation`，基于 Lyndon 词算法（已安装 pydivsufsort 时使用它）。
- **API:** `dnakit.molecular_lab.smallest_circular_rotation(sequence, *, verbose=True)`
- **输入:** 一条 DNA 字符串。
- **输出:** `sequence`（旋转后的序列）、`length` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import smallest_circular_rotation
print(smallest_circular_rotation("CAGT")["sequence"])
print(smallest_circular_rotation("GATTACA")["sequence"])
```

- **示例输出:**

```text
AGTC
ACAGATT
```

## 4) `dnakit.molecular_lab.generate_random_orf()`

- **作用:** 生成随机开放阅读框：ATG、随机有义密码子和一个终止密码子。
- **处理方法:** 调用 pydna 的 `randomORF`。种子只在本次调用内设置；调用前保存全局随机状态，调用后恢复。给定 `max_codons` 时，长度在两者之间按三角分布抽取。
- **API:** `dnakit.molecular_lab.generate_random_orf(codons, *, max_codons=None, seed=0, verbose=True)`
- **输入:** `codons` 为 2 到 10000 的整数，包含起始与终止密码子；可选 `max_codons`，不小于 `codons`。
- **输出:** `sequence`、`length`（为 `codons` 的 3 倍）、`codon_count`、`seed`、`max_codons` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import generate_random_orf
result = generate_random_orf(10, seed=1)
print(result["sequence"], result["codon_count"])
```

- **示例输出:**

```text
ATGACGGACGCAGTGTCTCTGTCGATATAG 10
```

## 5) `dnakit.molecular_lab.annotate_primer_binding_sites()`

- **作用:** 标出两条引物在产物两端的结合位置。
- **处理方法:** 检查正向引物是否为产物的起始序列，反向引物（5' 到 3' 书写）是否为产物末端的反向互补序列。随后调用 pydna 的 `annotate_primer_binding_sites`，在两端各添加一个 `primer_bind` 特征。
- **API:** `dnakit.molecular_lab.annotate_primer_binding_sites(product, forward_primer, reverse_primer, *, verbose=True)`
- **输入:** 产物与两条引物的 DNA 字符串。不匹配时抛出 `ConfigurationError`。
- **输出:** `product_length`、`feature_count`、`features`（每项含 `type`、`start`、`end`、`strand` 和 `label`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import annotate_primer_binding_sites
result = annotate_primer_binding_sites(
    "ATGCGTACGTTAGCATCGGATCCAAAGGGCCC", "ATGCGTACGTTAGCAT", "GGGCCCTTT"
)
for feature in result["features"]:
    print(feature["label"], feature["start"], feature["end"], feature["strand"])
```

- **示例输出:**

```text
forward_primer 0 16 1
reverse_primer 23 32 -1
```

## 6) `dnakit.molecular_lab.suggest_pcr_program()`

- **作用:** 用一对引物模拟 PCR，并返回 pydna 建议的循环程序。
- **处理方法:** 运行 pydna 的 `amplify.pcr`，每条引物须至少有 `min_annealing` 个碱基结合。然后按所选聚合酶应用其 Tm 与退火温度规则：`taq`（`program`）、`dbd`（`dbd_program`）或 `pfu_sso7d`（`pfu_sso7d_program`）。
- **API:** `dnakit.molecular_lab.suggest_pcr_program(template, forward_primer, reverse_primer, *, polymerase="taq", min_annealing=14, verbose=True)`
- **输入:** 模板 DNA、两条引物、`polymerase`，以及 8 到 60 之间的 `min_annealing`。
- **输出:** `product_length`、`program`（文本块）、`program_lines`、`polymerase` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import suggest_pcr_program
A = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGG"
result = suggest_pcr_program(A, A[:22], A[-22:].translate(str.maketrans("ACGT", "TGCA"))[::-1])
print(result["product_length"], result["polymerase"])
for line in result["program_lines"]:
    print(line)
```

- **示例输出:**

```text
104 taq
|95°C|95°C               |    |tmf:58.1
|____|_____          72°C|72°C|tmr:65.4
|3min|30s  \ 51.8°C _____|____|45s/kb
|    |      \______/ 0:30|5min|GC 41%
|    |       30s         |    |104bp
```

## 7) `dnakit.molecular_lab.ligate_blunt_fragments()`

- **作用:** 枚举平末端双链 DNA 片段连接后的产物。
- **处理方法:** 以 `allow_blunt=True` 调用 pydna 的 `ligation_assembly`。不表示酶切产生的粘性末端。产物先按长度、再按序列排序；`product_count` 给出总数。
- **API:** `dnakit.molecular_lab.ligate_blunt_fragments(fragments, *, allow_partial_overlap=False, circular_only=False, max_products=50, verbose=True)`
- **输入:** 2 到 20 条平末端 DNA 字符串；`max_products` 为 1 到 200。
- **输出:** `fragment_count`、`product_count`、`products`（每项含 `sequence`、`length` 和 `circular`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import ligate_blunt_fragments
first = "ATGCGTACGTTAGCATCGGATCCAAAGGGCCCTTTAAAGGAGATCGTTTAGCGCCC"
second = "TTTAGCGCCCAAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGC"
result = ligate_blunt_fragments([first, second])
print(result["product_count"], [row["length"] for row in result["products"]])
print(all(row["circular"] for row in result["products"]))
```

- **示例输出:**

```text
2 [116, 116]
True
```

## 8) `dnakit.molecular_lab.gateway_recombination_products()`

- **作用:** 枚举 Gateway BP 或 LR 重组反应得到的 DNA 产物。
- **处理方法:** 对所选反应调用 pydna 的 `gateway_assembly`。`greedy=True` 使用范围更宽的位点定义，否则使用保守位点。片段的 DNA 中必须含有 att 位点。
- **API:** `dnakit.molecular_lab.gateway_recombination_products(fragments, *, reaction_type="BP", greedy=False, circular_only=False, multi_site_only=False, max_products=50, verbose=True)`
- **输入:** 2 到 20 条含 att 位点的 DNA 字符串；`reaction_type` 为 `BP` 或 `LR`。
- **输出:** `reaction_type`、`site_set`、`product_count`、`products` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import gateway_recombination_products
from Bio.Data.IUPACData import ambiguous_dna_values
from pydna.gateway import recombinase_dict
site = recombinase_dict["BP"]["conservative"].recombinases[0]
def expand(text):
    return "".join(ambiguous_dna_values[ch][0] if ch in ambiguous_dna_values else ch for ch in text.upper())
att_b, att_p = expand(site.site1), expand(site.site2)
linear = att_b + "TTTTACGTACGTAAAGGGCCCAAATTTGGGGCCAAGGCCTTAAGGCCAA" + att_b
plasmid = att_p + "AAAACCCGGGTTTAAACCCGGG" + att_p
result = gateway_recombination_products([linear, plasmid], reaction_type="BP")
print(result["reaction_type"], result["site_set"], result["product_count"])
print([row["length"] for row in result["products"]])
```

- **示例输出:**

```text
BP conservative 7
[96, 157, 166, 227, 324, 350, 411]
```

## 9) `dnakit.molecular_lab.simulate_genome_edit()`

- **作用:** 模拟 DNA 序列上的一次基因组编辑事件。
- **处理方法:** 分派到 pydna 的 `homologous_recombination_*`、`cre_lox_*`、`recombinase_*` 或 `crispr_integration`。整合需要片段；切除与倒位不需要。CRISPR 需要能切割基因组的向导间隔序列。重组酶需要两个位点（大写侧翼、小写同源核心），例如 `("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")`。CRE-lox 整合需要环状片段，因此设 `inserts_circular=True`。`limit`（最小同源长度）只用于同源重组与 CRISPR。
- **API:** `dnakit.molecular_lab.simulate_genome_edit(genome, *, mechanism, operation, inserts=(), guides=(), recombinase_sites=None, genome_circular=False, inserts_circular=False, limit=40, max_products=50, verbose=True)`
- **输入:** 长度不超过 20000 碱基的基因组 DNA；`mechanism`（`homologous_recombination`、`cre_lox`、`recombinase` 或 `crispr`）；`operation`（`integration`、`excision` 或 `excision_or_inversion`；CRISPR 只允许 integration）。
- **输出:** `mechanism`、`edit`、`genome_length`、`insert_count`、`guide_count`、`limit`、`product_count`、`products` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import simulate_genome_edit
result = simulate_genome_edit(
    "CCCCCCATAACTTCGTATAGCATACATTATACGAAGTTATGGGGGATAACTTCGTATAGCATACATTATACGAAGTTATAAAAA",
    mechanism="cre_lox",
    operation="excision",
)
print(result["edit"], result["product_count"], [row["length"] for row in result["products"]])
```

- **示例输出:**

```text
excision 2 [39, 45]
```

- **示例代码:**

```python
from dnakit.molecular_lab import simulate_genome_edit
result = simulate_genome_edit(
    "aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac",
    mechanism="crispr",
    operation="integration",
    inserts=["aaccggttAAAAAAAAAttcaaagcac"],
    guides=["ttcaatgcaaacagtaatga"],
    limit=8,
)
print(result["mechanism"], result["limit"], [row["sequence"] for row in result["products"]])
```

- **示例输出:**

```text
crispr 8 ['AACCGGTTAAAAAAAAATTCAAAGCAC']
```

## 10) `dnakit.molecular_lab.assemble_recombinase_fragments()`

- **作用:** 枚举 DNA 片段经重组酶组装得到的产物。
- **处理方法:** 用由 `recombinase_sites` 构成的双位点重组酶调用 pydna 的 `recombinase_assembly`，格式与 `simulate_genome_edit` 相同。每条片段的 DNA 中含有一个位点，相匹配的位点把片段连接成产物。
- **API:** `dnakit.molecular_lab.assemble_recombinase_fragments(fragments, *, recombinase_sites, circular_only=False, max_products=50, verbose=True)`
- **输入:** 2 到 20 条含位点的 DNA 字符串，以及位点对。
- **输出:** `fragment_count`、`product_count`、`products` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import assemble_recombinase_fragments
sites = ("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
result = assemble_recombinase_fragments(
    ["CCCCATGCCCTAAaaTTGGGGGGGG", "AAaaTTTTTTTCCCTTTTTTTTT"], recombinase_sites=sites
)
print(result["product_count"], [row["length"] for row in result["products"]])
```

- **示例输出:**

```text
2 [14, 34]
```

## 11) `dnakit.molecular_lab.extract_protospacers()`

- **作用:** 提取紧邻向导支架序列之前的 Cas9 原间隔序列。
- **处理方法:** 以 `cas9` 类调用 pydna 的 `protospacer`：20 个碱基的间隔序列后接支架 `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG`。同时搜索两条链；`circular=True` 时还跨越起点搜索。
- **API:** `dnakit.molecular_lab.extract_protospacers(guide_construct, *, cas="cas9", circular=False, verbose=True)`
- **输入:** 向导构建体 DNA 字符串；`cas` 只能为 `cas9`，这是 pydna 唯一定义的核酸酶。
- **输出:** `protospacer_count`、`protospacers`（已排序）、`protospacer_length`、`scaffold`、`circular` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import extract_protospacers
from pydna.crispr import cas9
construct = "ATTA" + "GACGTACGATCGTAGCTAGG" + cas9.scaffold + "TTTT"
result = extract_protospacers(construct)
print(result["protospacer_count"], result["protospacers"])
```

- **示例输出:**

```text
1 ['GACGTACGATCGTAGCTAGG']
```

## 12) `dnakit.molecular_lab.parse_sequence_text()`

- **作用:** 把 EMBL、GenBank 或 FASTA 文本解析为 DNA 记录。
- **处理方法:** 在内存中对文本调用 pydna 的 `parsers.embl_gb_fasta`，无需文件。记录保持输入顺序。
- **API:** `dnakit.molecular_lab.parse_sequence_text(text, *, verbose=True)`
- **输入:** 不超过 1 MB、至少含一条记录的文本。
- **输出:** `record_count`、`records`（每项含 `name`、`length`、`topology` 和 `sequence`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import parse_sequence_text
result = parse_sequence_text(">one\nATGCCC\n>two\nGGGAAA\n")
for record in result["records"]:
    print(record["name"], record["length"], record["topology"], record["sequence"])
```

- **示例输出:**

```text
one 6 linear ATGCCC
two 6 linear GGGAAA
```

## 13) `dnakit.molecular_lab.design_user_cloning_primers()`

- **作用:** 设计 USER 克隆引物的尾巴，使相邻片段带有重叠末端。
- **处理方法:** 用各自的引物通过 pydna 的 `amplify.pcr` 扩增模板，再以 `max_overlap`、`min_overlap` 和 `max_tail` 调用 pydna 的 `user_assembly_design`。尾巴中的尿嘧啶按 pydna 的写法保留为 `U`。
- **API:** `dnakit.molecular_lab.design_user_cloning_primers(fragments, *, max_overlap=15, min_overlap=4, max_tail=50, verbose=True)`
- **输入:** 2 到 20 个映射，每个含 `template`、`forward_primer` 和 `reverse_primer` 的 DNA 字符串；引物必须能扩增模板。
- **输出:** `fragment_count`、`fragments`（每项含 `forward_primer`、`reverse_primer` 和 `product_length`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import design_user_cloning_primers
A = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGG"
B = "GGCGAAGGGATGCTTACCGTTGGTACACACCCCAATCTTTTACTGGAGCGTAAAGGTTACGCATGGCAATTAAATCGCTGGATTGTCGATCCACGCCATCTGG"
def reverse_complement(text):
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]
fragments = [
    {"template": A, "forward_primer": A[:22], "reverse_primer": reverse_complement(A[-22:])},
    {"template": B, "forward_primer": B[:22], "reverse_primer": reverse_complement(B[-22:])},
]
result = design_user_cloning_primers(fragments)
for row in result["fragments"]:
    print(row["forward_primer"], row["reverse_primer"], row["product_length"])
```

- **示例输出:**

```text
ATGGCTAAAGGTGAAGAAC ATCCCTTCGCCCCUTCACCCTCTCCA 104
AGGGGCGAAGGGAUGCT CCAGATGGCGTGGAT 103
```

## 14) `dnakit.molecular_lab.simulate_type2s_assembly()`

- **作用:** 模拟线性部件通过 Golden Gate（Type-2S）组装成环状产物。
- **处理方法:** 构建 DnaCauldron 的 `generate_type2s_restriction_mix`，并请它给出环状组装结果。`enzyme="auto"` 时由 DnaCauldron 的 `autoselect_enzyme` 在 BbsI、AarI、BsaI、BsmBI 和 SapI 中选择。
- **API:** `dnakit.molecular_lab.simulate_type2s_assembly(parts, *, enzyme="auto", max_products=50, verbose=True)`
- **输入:** 从部件名到线性 DNA 字符串的映射，字符串两端带有粘性末端序列。
- **输出:** `enzyme`、`part_names`、`circular_product_count`、`products`（每项含 `sequence`、`length` 和 `circular`）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import simulate_type2s_assembly
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

parts = {"p1": part("AATG", "GCTT", INS1), "p2": part("GCTT", "AATG", INS2)}
result = simulate_type2s_assembly(parts, enzyme="BsmBI")
print(result["enzyme"], result["circular_product_count"], [row["length"] for row in result["products"]])
```

- **示例输出:**

```text
BsmBI 1 [152]
```

## 15) `dnakit.molecular_lab.clone_into_donor_vector()`

- **作用:** 通过替换 Type-2S 填充片段，把插入片段克隆到供体载体中。
- **处理方法:** 调用 DnaCauldron 的 `swap_donor_vector_part`。供体为环状，其填充片段两侧各有一个位点对。插入片段也为环状，位点对位于插入片段两侧，并且必须携带与填充片段相同的两个粘性末端；否则上游失败，本函数抛出 `ConfigurationError`。
- **API:** `dnakit.molecular_lab.clone_into_donor_vector(donor_sequence, insert_sequence, *, enzyme="BsmBI", verbose=True)`
- **输入:** 供体与插入片段的 DNA 字符串，以及 `enzyme`（BbsI、AarI、BsaI、BsmBI 或 SapI）。
- **输出:** `enzyme`、`product_sequence`、`product_length` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import clone_into_donor_vector
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
insert = part("AATG", "GCTT", INS2)
result = clone_into_donor_vector(donor, insert, enzyme="BsmBI")
print(len(donor), result["product_length"])
```

- **示例输出:**

```text
794 794
```

## 16) `dnakit.molecular_lab.assign_parts_to_backbones()`

- **作用:** 按粘性末端匹配，把每个环状 Golden Gate 部件分配到对应的骨架。
- **处理方法:** 以 BsmBI 和 `process_parts_with_backbone=False` 调用 DnaCauldron 的 `insert_parts_on_backbones`。上游的骨架检测只支持 BsmBI，因此本封装只接受 BsmBI。每个部件报告为 `assigned`（附产物）、`no_backbone` 或 `error`；每个部件的状态在 `status_per_part` 中。
- **API:** `dnakit.molecular_lab.assign_parts_to_backbones(parts, backbones, *, min_backbone_length=500, verbose=True)`
- **输入:** 两个从名称到环状 DNA 字符串的映射：部件与骨架。每个骨架含一个 BsmBI 填充片段。
- **输出:** `part_count`、`assigned_count`、`assignments`（每项含 `part`、`status_per_part`、`already_on_backbone`、`backbone`、`product_sequence`、`product_length` 和 `error`）、`enzyme` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import assign_parts_to_backbones
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

backbone = part("AATG", "GCTT", "ACGT" * 10) + "ACGT" * 175
parts = {"p1": part("AATG", "GCTT", INS1), "p2": part("TTAA", "CCGG", INS2)}
result = assign_parts_to_backbones(parts, {"bb": backbone})
print(result["part_count"], result["assigned_count"])
for row in result["assignments"]:
    print(row["part"], row["status_per_part"], row["backbone"])
```

- **示例输出:**

```text
2 1
p1 assigned bb
p2 no_backbone None
```

## 17) `dnakit.molecular_lab.substitute_type2s_overhangs()`

- **作用:** 替换环状载体上的 Type-2S 粘性末端，从而在组装标准中移动部件。
- **处理方法:** 调用 DnaCauldron 的 `utils.substitute_overhangs`。只改变载体酶切后作为粘性末端出现的序列，改变的序列列在 `applied_overhangs` 中；`requested` 列出全部输入的键。
- **API:** `dnakit.molecular_lab.substitute_type2s_overhangs(circular_sequence, substitutions, *, enzyme="BsmBI", verbose=True)`
- **输入:** 环状 DNA 字符串；`substitutions` 把最多 8 个已有的 4 碱基粘性末端映射到新的 4 碱基粘性末端；`enzyme` 为 BbsI、AarI、BsaI、BsmBI 或 SapI。
- **输出:** `applied_overhangs`、`requested`、`product_sequence`、`product_length`、`enzyme` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import substitute_type2s_overhangs
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
result = substitute_type2s_overhangs(donor, {"AATG": "CCCC"}, enzyme="BsmBI")
print(result["applied_overhangs"], result["product_length"])
```

- **示例输出:**

```text
['AATG'] 794
```

## 18) `dnakit.molecular_lab.record_has_backbone()`

- **作用:** 判断环状 BsmBI 构建体除插入片段外是否仍含骨架。
- **处理方法:** 用 BsmBI 切割构建体，把唯一的无位点片段作为插入片段。调用 DnaCauldron 的 `record_contains_backbone`；骨架为环中其余部分，当其长度大于 `min_backbone_length` 时 `has_backbone` 为真。上游忽略 enzyme 参数，因此本封装固定为 BsmBI。
- **API:** `dnakit.molecular_lab.record_has_backbone(circular_sequence, *, min_backbone_length=500, verbose=True)`
- **输入:** 环状 DNA 字符串以及 `min_backbone_length`。
- **输出:** `construct_length`、`insert_length`、`backbone_length`、`has_backbone`、`enzyme` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import record_has_backbone
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
result = record_has_backbone(donor)
print(result["construct_length"], result["insert_length"], result["backbone_length"], result["has_backbone"])
```

- **示例输出:**

```text
794 72 722 True
```

## 19) `dnakit.molecular_lab.dna_window_overlap()`

- **作用:** 返回两个以 (start, end) 表示的基因组窗口的共同区间。
- **处理方法:** 调用 DnaChisel 的 `windows_overlap`。仅在一个坐标处相接的窗口返回零长度重叠（与 DnaChisel 的报告一致）；此时 `overlaps` 为假。
- **API:** `dnakit.molecular_lab.dna_window_overlap(window1, window2, *, verbose=True)`
- **输入:** 两个 (start, end) 对，为非负整数且 start 不大于 end。
- **输出:** `overlap`（[start, end] 对，或 None）、`overlap_length`、`overlaps` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import dna_window_overlap
result = dna_window_overlap([0, 10], [5, 20])
print(result["overlap"], result["overlap_length"], result["overlaps"])
print(dna_window_overlap([0, 3], [5, 9])["overlap"])
```

- **示例输出:**

```text
[5, 10] 5 True
None
```

## 20) `dnakit.molecular_lab.subdivide_dna_window()`

- **作用:** 把 (start, end) 窗口切分为连续的、长度不超过 max_span 碱基的子窗口。
- **处理方法:** 调用 DnaChisel 的 `subdivide_window`。各段首尾相接并覆盖整个窗口；空窗口不产生任何子窗口。
- **API:** `dnakit.molecular_lab.subdivide_dna_window(window, *, max_span, verbose=True)`
- **输入:** 一个 (start, end) 对与正整数 `max_span`。
- **输出:** `window_count`、`windows`（每项为 [start, end] 对）和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import subdivide_dna_window
result = subdivide_dna_window([0, 10], max_span=4)
print(result["window_count"], result["windows"])
```

- **示例输出:**

```text
3 [[0, 4], [4, 8], [8, 10]]
```

## 21) `dnakit.molecular_lab.group_nearby_dna_segments()`

- **作用:** 把起点彼此靠近的 (start, end) 区段分为一组。
- **处理方法:** 先按起点、再按终点排序，然后调用 DnaChisel 的 `group_nearby_segments`。区段与上一区段起点之差小于 `max_start_gap`、且与本组首区段起点之差小于 `max_start_spread` 时加入当前组。`None` 表示不做该检验。
- **API:** `dnakit.molecular_lab.group_nearby_dna_segments(segments, *, max_start_gap=None, max_start_spread=None, verbose=True)`
- **输入:** (start, end) 对的列表，以及可选的正整数限制。
- **输出:** `segment_count`、`group_count`、`groups`（[start, end] 对的列表）、所用限制和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import group_nearby_dna_segments
segments = [[0, 5], [2, 4], [10, 12], [11, 20]]
result = group_nearby_dna_segments(segments, max_start_gap=3)
print(result["group_count"], result["groups"])
```

- **示例输出:**

```text
2 [[[0, 5], [2, 4]], [[10, 12], [11, 20]]]
```

## 22) `dnakit.molecular_lab.codon_backtranslation_table()`

- **作用:** 返回 Biopython 密码子表的回译表，每个键对应排序后的密码子列表。
- **处理方法:** 调用 DnaChisel 的 `get_backtranslation_table`。键为氨基酸字母、`*`（终止）、`START`（起始）以及歧义字母 X、B、J、Z；每个值是排序后的密码子列表。
- **API:** `dnakit.molecular_lab.codon_backtranslation_table(table_name="Standard", *, verbose=True)`
- **输入:** `table_name`，取 Biopython 无歧义 DNA 密码子表名，例如 `Standard` 或 `Vertebrate Mitochondrial`。
- **输出:** `table_name`、`entry_count`、`table` 和 `provenance`。
- **示例代码:**

```python
from dnakit.molecular_lab import codon_backtranslation_table
result = codon_backtranslation_table()
table = result["table"]
print(result["entry_count"], table["M"], table["*"], len(table["X"]))
```

- **示例输出:**

```text
26 ['ATG'] ['TAA', 'TAG', 'TGA'] 61
```

## 未封装的函数

这些行在账本中保留了精确状态和一行原因。

- Biopython 记录辅助函数、写入器与拓扑设置函数（`annotate_record`、`sequence_to_biopython_record`、`change_biopython_record_sequence`、`write_record`、`set_record_topology`、`record_is_linear`）：它们构建或编辑记录，而不是对序列进行计算。
- 文本与显示辅助函数（`assembly2str`、`dict_to_pretty_string`、`round_all_numbers_in_dict`、`score_to_formatted_string`、表格格式化函数、`gbtext_clean`）：输出不是序列结果。
- `reverse_complement_assembly` 与 `partially_digested_seqs_in_restriction_ligation_assembly` 作用于 pydna 的内部组装对象；后者对所试的酶连接产物没有返回任何部分酶切序列，因此没有可记录的示例。
- `pydna.utils.express` 是上游占位函数；`pydna.utils.cai` 需要 `cai2` 包，而 `dnakit-full` 中未安装（`needs-data:cai2-package`）；`parse_snapgene_history` 需要 SnapGene `.dna` 历史文件（`needs-data:snapgene-history-dna-file`）。
- DnaChisel 的 `dna_pattern_to_regexpr` 把 N、R、Y 映射到错误的字符类，因此不封装；它被记为 `dnakit.patterns.scan_motif` 的重复项，后者能正确匹配 IUPAC 字符。
- pandas 报告构建函数（`constraints_before_after_dataframe`、`objectives_before_after_dataframe`、`records_from_breaches_dataframe`）整理已完成的 DnaChisel 问题的结果，而非序列。
