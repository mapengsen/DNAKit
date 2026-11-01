# Biopython 限制性酶切分析、CAPS、BioSQL 与基因组图

本文档涵盖对 [Biopython](https://biopython.org) 1.88 的三类封装：`dnakit.sequence` 中的 `Bio.Restriction` 与 `Bio.CAPS`（使用 Biopython 自带的 REBASE 酶数据；位置从 1 起，为上链切口后的第一个碱基）、`dnakit.databases` 中的 `BioSQL`（在由 BioSQL schema 建立的 SQLite、MySQL 或 PostgreSQL 数据库中存取记录）和 `dnakit.visualization` 中的 `Bio.Graphics`（用 ReportLab 绘制基因组交叉连接图和染色体表意图）。用 `python -m pip install "dnakit[bio-formats,reports]"` 安装；各包只在函数运行时导入。BioSQL 示例需要 BioSQL 的 SQLite schema 文件，Biopython 不自带、DNAKit 也不再分发：请从 https://github.com/biosql/biosql 下载 `biosqldb-sqlite.sql`，并把环境变量 `DNAKIT_BIOSQL_SQLITE_SCHEMA` 设为其路径。示例从模块 `dnakit.sequence.biopython_p19_restriction_functions`、`dnakit.databases.biopython_p19_biosql_functions` 和 `dnakit.visualization.biopython_p19_graphics_functions` 导入。

## 示例数据

`records.gb` 是供 BioSQL 示例使用的双记录 GenBank 文件（`REC1` 30 个碱基、`REC2` 20 个碱基，各含一个 `CDS` 特征）。酶切示例使用 38 个碱基的序列 `GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC`，其中有两个 EcoRI 位点，BamHI、SmaI、KpnI 各一个，没有 HindIII 位点。

**示例文件：** `records.gb`

```text
LOCUS       REC1                      30 bp    DNA              UNK 01-JAN-1980
DEFINITION  test REC1.
ACCESSION   ACC0
VERSION     ACC0.1
KEYWORDS    .
SOURCE      .
  ORGANISM  Testus x
            .
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="g0"
ORIGIN
        1 atggctaaag gtttcccaga actggtttaa
//
LOCUS       REC2                      20 bp    DNA              UNK 01-JAN-1980
DEFINITION  test REC2.
ACCESSION   ACC1
VERSION     ACC1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  Testus x
            .
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="g1"
ORIGIN
        1 ggggccccaa aattttacgt
//
```

## 1) `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzyme_properties()`

- **作用：** 描述限制性内切酶：识别位点、切割模式、末端类型、突出端、出现频率和供应商。
- **处理方法：** 读取 Biopython 自带的 REBASE 数据（`Bio.Restriction.<酶>`）：`elucidate()` 给出切割模式（`^` 为上链切口，`_` 为下链切口），`is_blunt`、`is_5overhang`、`is_3overhang` 给出末端类型，`ovhg` 与 `ovhgseq` 给出突出端，`frequency()` 给出每个预期位点对应的随机碱基数，`suppl`、`supplier_list` 和 `is_comm` 给出供应商。涵盖 Biopython 酶类的末端类型、切割模式、位点频率和供应商查询。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzyme_properties(enzymes)`。
- **输入：** 酶名列表（REBASE 名称，如 `EcoRI`）。
- **输出：** 每个酶一条记录：`site`、`size`、`elucidate`、`end_type`、`overhang`、`overhang_sequence`、`palindromic`、`frequency`、`commercial`、`supplier_codes` 和 `suppliers`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_enzyme_properties

for item in restriction_enzyme_properties(["EcoRI", "SmaI", "KpnI", "BsaI"])["enzymes"]:
    print(item["name"], item["site"], item["elucidate"], item["end_type"], repr(item["overhang_sequence"]), item["palindromic"], item["frequency"], item["commercial"])
print(restriction_enzyme_properties(["EcoRI"])["enzymes"][0]["suppliers"][:3])
```

- **示例输出：**

```text
EcoRI GAATTC G^AATT_C 5_overhang 'AATT' True 4096.0 True
SmaI CCCGGG CCC^_GGG blunt '' True 4096.0 True
KpnI GGTACC G_GTAC^C 3_overhang 'GTAC' True 4096.0 True
BsaI GGTCTC GGTCTCN^NNNN_N 5_overhang 'NNNN' False 4096.0 True
['Thermo Fisher Scientific', 'Minotech Biotechnology', 'SibEnzyme Ltd.']
```

## 2) `dnakit.sequence.biopython_p19_restriction_functions.restriction_related_enzymes()`

- **作用：** 查找与某个酶相关的酶：相容末端、新同切点异源酶、同切点同源酶或同工酶。
- **处理方法：** 调用 `Bio.Restriction.<酶>.<relation>`。`compatible_end` 列出其末端可与该酶末端连接的酶，`neoschizomers` 位点相同但切割位置不同，`equischizomers` 位点与切割位置都相同，`isoschizomers` 为两者之和。默认在全部 REBASE 酶中搜索，`among` 可限定范围。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.restriction_related_enzymes(enzyme, relation='isoschizomers', *, among=None)`。
- **输入：** 酶名、`relation` 和可选的酶名列表 `among`。
- **输出：** `site`、`count` 和排序后的 `enzymes`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_related_enzymes

for enzyme, relation in (("KpnI", "isoschizomers"), ("KpnI", "equischizomers"), ("SmaI", "neoschizomers"), ("EcoRI", "compatible_end")):
    result = restriction_related_enzymes(enzyme, relation)
    print(enzyme, relation, result["count"], result["enzymes"][:4])
print(restriction_related_enzymes("SmaI", "neoschizomers", among=["XmaI", "SmaI", "KpnI"])["enzymes"])
```

- **示例输出：**

```text
KpnI isoschizomers 2 ['Acc65I', 'Asp718I']
KpnI equischizomers 0 []
SmaI neoschizomers 3 ['Cfr9I', 'TspMI', 'XmaI']
EcoRI compatible_end 51 ['AarI', 'Acc36I', 'AceIII', 'AcsI']
['XmaI']
```

## 3) `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzymes_by_supplier()`

- **作用：** 列出某个供应商出售的酶，或列出所有供应商。
- **处理方法：** 使用 `Bio.Restriction.RestrictionBatch.add_supplier` 和 `suppl_codes`。供应商用 REBASE 单字母代码（`N` = New England Biolabs）或全名指定。不给 `supplier` 时返回代码、名称及各供应商出售酶数量的表。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzymes_by_supplier(supplier=None)`。
- **输入：** 可选的 `supplier` 代码或名称。
- **输出：** 给定供应商时：`supplier_code`、`supplier`、`enzyme_count` 和 `enzymes`；否则为 `suppliers` 表。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_enzymes_by_supplier

table = restriction_enzymes_by_supplier()["suppliers"]
print(len(table), table[0])
result = restriction_enzymes_by_supplier("N")
print(result["supplier"], result["enzyme_count"], result["enzymes"][:4])
print(restriction_enzymes_by_supplier("New England Biolabs")["enzyme_count"])
```

- **示例输出：**

```text
15 {'code': 'B', 'name': 'Thermo Fisher Scientific', 'enzyme_count': 211}
New England Biolabs 238 ['AatII', 'AbaSI', 'Acc65I', 'AccI']
238
```

## 4) `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_enzymes()`

- **作用：** 按性质筛选酶：末端类型、是否商品化、供应商、回文位点或位点长度。
- **处理方法：** 对 `enzymes`（默认全部酶）应用 Biopython 酶类的 `is_blunt`、`is_5overhang`、`is_3overhang`、`is_comm`、`suppl`、`is_palindromic`、`size` 和 `is_defined` 判据，所有给定条件须同时满足。Biopython 1.88 自带的 `RestrictionBatch.lambdasplit` 与 `split` 返回空集合，因此直接应用这些判据。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_enzymes(enzymes=None, *, end_type=None, commercial=None, supplier=None, palindromic=None, min_site_size=None, max_site_size=None, defined_site=None)`。
- **输入：** 可选的 `enzymes`、`end_type`（`blunt`、`5_overhang`、`3_overhang`）、`commercial`、`supplier` 代码、`palindromic`、`min_site_size`、`max_site_size` 和 `defined_site`。
- **输出：** `candidate_count`、`count` 和排序后的 `enzymes`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import filter_restriction_enzymes

blunt = filter_restriction_enzymes(end_type="blunt")
print(blunt["count"], "of", blunt["candidate_count"])
print(filter_restriction_enzymes(["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], end_type="5_overhang", palindromic=True, min_site_size=6)["enzymes"])
sold = filter_restriction_enzymes(commercial=True, supplier="N", max_site_size=4, defined_site=True)
print(sold["count"], sold["enzymes"][:5])
```

- **示例输出：**

```text
137 of 1088
['BamHI', 'EcoRI', 'HindIII']
24 ['AciI', 'AluI', 'BfaI', 'BstUI', 'CviAII']
```

## 5) `dnakit.sequence.biopython_p19_restriction_functions.restriction_analysis_report()`

- **作用：** 在序列中搜索限制性酶切位点，并以列表、计数表或图谱形式输出。
- **处理方法：** 对所选酶运行 `Bio.Restriction.Analysis`，用 `print_as`（`list`、`number` 或 `map`）和 `print_that` 输出。位置从 1 起：上链切口后的第一个碱基。环状分子用 `linear=False`。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.restriction_analysis_report(sequence, enzymes, *, linear=True, view='list')`。
- **输入：** DNA `sequence`（IUPAC 字母）、酶名列表、`linear` 和 `view`。
- **输出：** `cuts`（酶到位置）、`no_cut` 和文本 `report`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_analysis_report

result = restriction_analysis_report("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], view="map")
print(result["cuts"], result["no_cut"])
tidy = lambda text: "\n".join(line.rstrip() for line in text.splitlines())
print(tidy(result["report"]))
print(tidy(restriction_analysis_report("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], view="number")["report"]))
```

- **示例输出：**

```text
{'BamHI': [11], 'EcoRI': [2, 20], 'KpnI': [38], 'SmaI': [27]} ['HindIII']
 2 EcoRI
 |
 |        11 BamHI
 |        |
 |        |        20 EcoRI
 |        |        |
 |        |        |      27 SmaI
 |        |        |      |
 |        |        |      |          38 KpnI
 |        |        |      |          |
GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC
||||||||||||||||||||||||||||||||||||||
CTTAAGTTTCCTAGGTTTCTTAAGGGCCCAAACCATGG
1                                   38
   Enzymes which do not cut the sequence.
HindIII
enzymes which cut 1 times :
BamHI      :  11.
KpnI       :  38.
SmaI       :  27.
enzymes which cut 2 times :
EcoRI      :  2, 20.
   Enzymes which do not cut the sequence.
HindIII
```

## 6) `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_results()`

- **作用：** 执行酶切搜索，只保留符合末端类型、切割次数、名称、位点长度或位置条件的结果。
- **处理方法：** 使用 `Bio.Restriction.Analysis` 的过滤方法并依次应用给定条件：`blunt`、`overhang5`、`overhang3`；`with_N_sites` 或 `with_number_list`；`with_sites` 或 `without_site`；酶名和位点长度（直接应用，因为 Biopython 1.88 的 `with_name` 与 `with_site_size` 不能与前一步结果一起使用）；以及带 1 起始 `start`、`end` 的位置过滤 `only_between`、`between`、`show_only_between`、`only_outside`、`outside` 和 `do_not_cut`。与 Biopython 一致，所选末端类型中不切割的酶会以空列表保留，`cuts="cutting"` 可去掉它们。若某一步没有结果则链条终止（Biopython 会重新回到完整结果）。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_results(sequence, enzymes, *, linear=True, end_type=None, cut_count=None, cuts=None, names=None, site_sizes=None, position_filter=None, start=None, end=None)`。
- **输入：** `sequence`、酶名、`linear`，以及至少一个条件：`end_type`、`cut_count`（数或列表）、`cuts`（`cutting` 或 `non_cutting`）、`names`、`site_sizes`，或带 `start` 与 `end` 的 `position_filter`。
- **输出：** 已应用的 `filters`、`match_count` 和 `results`（酶到位置）。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import filter_restriction_results

def show(**criteria):
    print(criteria, filter_restriction_results("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], **criteria)["results"])

show(end_type="blunt")
show(end_type="5_overhang")
show(end_type="5_overhang", cuts="cutting")
show(cut_count=2)
show(cuts="non_cutting")
show(names=["EcoRI", "KpnI"])
show(site_sizes=[6], position_filter="between", start=1, end=20)
show(position_filter="show_only_between", start=15, end=38)
```

- **示例输出：**

```text
{'end_type': 'blunt'} {'SmaI': [27]}
{'end_type': '5_overhang'} {'BamHI': [11], 'EcoRI': [2, 20], 'HindIII': []}
{'end_type': '5_overhang', 'cuts': 'cutting'} {'BamHI': [11], 'EcoRI': [2, 20]}
{'cut_count': 2} {'EcoRI': [2, 20]}
{'cuts': 'non_cutting'} {'HindIII': []}
{'names': ['EcoRI', 'KpnI']} {'EcoRI': [2, 20], 'KpnI': [38]}
{'site_sizes': [6], 'position_filter': 'between', 'start': 1, 'end': 20} {'BamHI': [11], 'EcoRI': [2, 20]}
{'position_filter': 'show_only_between', 'start': 15, 'end': 38} {'EcoRI': [20], 'SmaI': [27]}
```

## 7) `dnakit.sequence.biopython_p19_restriction_functions.caps_differential_cut_sites()`

- **作用：** 查找 CAPS 标记候选：只切割部分比对序列的酶切位点。
- **处理方法：** 调用 `Bio.CAPS.CAPSMap`（Konieczny 和 Ausubel 1993）。对每个酶，当某个比对位置至少有一条序列被切、至少一条不被切时报告该位点。不接受缺口字符，因为 `Bio.Restriction` 无法搜索缺口。
- **API：** `dnakit.sequence.biopython_p19_restriction_functions.caps_differential_cut_sites(alignment, enzymes)`。
- **输入：** 序列名到等长比对 DNA（IUPAC 字母）的映射，以及酶名列表。
- **输出：** `site_count` 和 `sites`：酶、切割处的 0 起始比对位置 `start`，以及被切（`cuts_in`）与不被切（`blocked_in`）的序列。
- **示例代码：**

```python
from dnakit.sequence.biopython_p19_restriction_functions import caps_differential_cut_sites

alignment = {"a": "AAGAATTCGGGAAGGATCCT", "b": "AAGAATTAGGGAAGGATCCT", "c": "AAGAATTCGGGAAGGATCGT"}
result = caps_differential_cut_sites(alignment, ["EcoRI", "BamHI"])
print(result["site_count"])
for site in result["sites"]:
    print(site)
```

- **示例输出：**

```text
2
{'enzyme': 'EcoRI', 'start': 3, 'cuts_in': ['a', 'c'], 'blocked_in': ['b']}
{'enzyme': 'BamHI', 'start': 14, 'cuts_in': ['a', 'b'], 'blocked_in': ['c']}
```

## 8) `dnakit.databases.biopython_p19_biosql_functions.biosql_load_records()`

- **作用：** 把文件中的序列记录载入 BioSQL 子数据库。
- **处理方法：** 用 `Bio.SeqIO` 读文件，子数据库不存在时创建（`BioSeqDatabase.new_database`），用 `BioSeqDatabase.load`（`BioSQL.Loader.DatabaseLoader.load_seqrecord`）载入记录并提交；失败时回滚。不访问 NCBI 分类学。对 SQLite，`schema_sql`（BioSQL 的 `biosqldb-sqlite.sql` 文件，Biopython 不自带）用于创建新的数据库文件；不会重建已存在的数据库。MySQL（`mysql.connector`）和 PostgreSQL（`psycopg2`）数据库通过 `driver`、`host`、`user`、`password` 和 `port` 访问。
- **API：** `dnakit.databases.biopython_p19_biosql_functions.biosql_load_records(database, namespace, records_path, *, record_format='genbank', driver='sqlite3', schema_sql=None, host=None, user=None, password=None, port=None)`。
- **输入：** `database`（SQLite 文件，服务器数据库则为库名）、`namespace`、`records_path` 与 `record_format`（默认 `genbank`）、`driver`，以及可选的 `schema_sql`。示例从环境变量 `DNAKIT_BIOSQL_SQLITE_SCHEMA` 读取 schema 路径（下载 https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-sqlite.sql）。
- **输出：** `records_loaded`、`records_in_namespace`，以及子数据库是否为新建。
- **示例代码：**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
result = biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
print(result["records_loaded"], result["records_in_namespace"], result["created_namespace"])
```

- **示例输出：**

```text
2 2 True
```

## 9) `dnakit.databases.biopython_p19_biosql_functions.biosql_get_record()`

- **作用：** 按标识符从 BioSQL 子数据库取出一条记录。
- **处理方法：** `id_type` 选择标识符类型：`display_id`（`BioSeqDatabase.get_Seq_by_id`）、`accession`（`get_Seq_by_acc`）、`version`（`get_Seq_by_ver`），或 `primary_id`、`name`、`gi`（`BioSeqDatabase.lookup`）。给定 `output_path`（新文件）时还会用 `Bio.SeqIO` 以 `output_format` 写出该记录。
- **API：** `dnakit.databases.biopython_p19_biosql_functions.biosql_get_record(database, namespace, identifier, *, id_type='display_id', include_sequence=True, max_sequence_length=1000000, output_path=None, output_format='genbank', driver='sqlite3', host=None, user=None, password=None, port=None)`。
- **输入：** `database`、`namespace`、`identifier`、`id_type`、`include_sequence`、`max_sequence_length`、可选的 `output_path` 与 `output_format`，以及连接参数。
- **输出：** `record`：id、名称、描述、长度、序列、特征摘要、注释和交叉引用。
- **示例代码：**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_get_record

record = biosql_get_record(database, "demo", "REC1")["record"]
print(record["id"], record["name"], record["length"], record["sequence"], record["feature_types"])
print(biosql_get_record(database, "demo", "ACC1", id_type="accession")["record"]["name"])
print(biosql_get_record(database, "demo", "ACC0.1", id_type="version", include_sequence=False)["record"]["length"])
```

- **示例输出：**

```text
ACC0.1 REC1 30 ATGGCTAAAGGTTTCCCAGAACTGGTTTAA ['CDS']
REC2
30
```

## 10) `dnakit.databases.biopython_p19_biosql_functions.biosql_get_subsequence()`

- **作用：** 不载入整条记录，直接返回已存序列的一个片段。
- **处理方法：** 在数据库内用 SQL `SUBSTR` 截取片段 `[start, end)`（0 起始，半开区间）（`BioSQL.BioSeqDatabase.Adaptor.get_subseq_as_string`）。记录按 `identifier` 与 `id_type` 选择，同 `biosql_get_record`。
- **API：** `dnakit.databases.biopython_p19_biosql_functions.biosql_get_subsequence(database, namespace, identifier, start, end, *, id_type='display_id', driver='sqlite3', host=None, user=None, password=None, port=None)`。
- **输入：** `database`、`namespace`、`identifier`、`start`、`end`、`id_type` 和连接参数。
- **输出：** 截取的 `sequence` 片段和完整的 `sequence_length`。
- **示例代码：**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_get_subsequence

result = biosql_get_subsequence(database, "demo", "REC1", 3, 12)
print(result["sequence"], result["sequence_length"])
print(biosql_get_subsequence(database, "demo", "ACC1.1", 16, 20, id_type="version")["sequence"])
```

- **示例输出：**

```text
GCTAAAGGT 30
ACGT
```

## 11) `dnakit.databases.biopython_p19_biosql_functions.biosql_list_entries()`

- **作用：** 列出 BioSQL 数据库的子数据库，或某个子数据库的条目。
- **处理方法：** 不给 `namespace` 时列出子数据库名及其条目数（`Adaptor.list_biodatabase_names`）。给定 `namespace` 时列出内部条目 id 与显示 id（`Adaptor.list_bioentry_ids`、`list_bioentry_display_ids`），最多 `max_entries` 个。
- **API：** `dnakit.databases.biopython_p19_biosql_functions.biosql_list_entries(database, *, namespace=None, max_entries=10000, driver='sqlite3', host=None, user=None, password=None, port=None)`。
- **输入：** `database`、可选的 `namespace`、`max_entries` 和连接参数。
- **输出：** 带条目数的 `namespaces`，或 `entry_count`、`entry_ids` 和 `display_ids`。
- **示例代码：**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_list_entries

print(biosql_list_entries(database)["namespaces"])
entries = biosql_list_entries(database, namespace="demo")
print(entries["entry_count"], entries["entry_ids"], entries["display_ids"])
```

- **示例输出：**

```text
[{'name': 'demo', 'entry_count': 2}]
2 [1, 2] ['REC1', 'REC2']
```

## 12) `dnakit.visualization.biopython_p19_graphics_functions.genome_diagram_crosslinks()`

- **作用：** 把基因组轨道画成上下排列的条带，并用交叉连接带相连。
- **处理方法：** 构建 `Bio.Graphics.GenomeDiagram.Diagram`，每个基因组一条轨道（奇数轨道层级为连接带留出空间），按链方向把特征画成方框或箭头，用 `CrossLink` 连接不同轨道的区域（`flip=True` 为反向匹配画出翻转的连接带）并写出图形。所有轨道使用最长轨道的比例尺。输出为 `pdf`、`svg` 或 `eps`；`png` 需要 ReportLab 的 Cairo 后端 `rlPyCairo`。输出文件必须不存在。
- **API：** `dnakit.visualization.biopython_p19_graphics_functions.genome_diagram_crosslinks(tracks, links, output_path, *, image_format='pdf', title=None, width_cm=24.0, track_height_cm=3.0)`。
- **输入：** `tracks`（`name`、`length`，`features` 含 `start`、`end`、`strand`、`label`、`color`）、`links`（`track_a`、`start_a`、`end_a`、`track_b`、`start_b`、`end_b`、`color`、`flip`；轨道用名称或序号）、新的 `output_path`、`image_format`、`title`、`width_cm` 和 `track_height_cm`。坐标 0 起始、半开区间；颜色为 `red` 之类名称或 `#RRGGBB`。
- **输出：** `output` 路径、`image_format`、`bytes` 以及所绘的 `tracks` 与 `links`。
- **示例代码：**

```python
import os, tempfile
from dnakit.visualization.biopython_p19_graphics_functions import genome_diagram_crosslinks

tracks = [
    {"name": "genomeA", "length": 1000, "features": [
        {"start": 100, "end": 300, "strand": 1, "label": "geneA", "color": "steelblue"},
        {"start": 500, "end": 700, "strand": -1, "label": "geneB", "color": "#cc3333"}]},
    {"name": "genomeB", "length": 800, "features": [{"start": 50, "end": 250, "strand": -1, "label": "geneA2"}]},
    {"name": "genomeC", "length": 900},
]
links = [
    {"track_a": "genomeA", "start_a": 100, "end_a": 300, "track_b": "genomeB", "start_b": 50, "end_b": 250, "flip": True, "color": "red"},
    {"track_a": 0, "start_a": 500, "end_a": 700, "track_b": 2, "start_b": 100, "end_b": 300},
]
target = os.path.join(tempfile.mkdtemp(), "links.svg")
result = genome_diagram_crosslinks(tracks, links, target, image_format="svg", title="demo")
print(result["image_format"], result["bytes"] > 1000, [t["name"] for t in result["tracks"]], len(result["links"]))
text = open(target).read()
print("geneA" in text, "genomeB" in text)
```

- **示例输出：**

```text
svg True ['genomeA', 'genomeB', 'genomeC'] 2
True True
```

## 13) `dnakit.visualization.biopython_p19_graphics_functions.chromosome_ideogram()`

- **作用：** 把染色体并排画成带端粒帽和标注特征的表意图。
- **处理方法：** 构建 `Bio.Graphics.BasicChromosome.Organism`，每个条目一条 `Chromosome`：两端为 `TelomereSegment` 帽，中间是 `ChromosomeSegment`，有特征时用 `AnnotatedChromosomeSegment`（链方向 +1 画在右侧，-1 画在左侧，0 横跨）。所有染色体共用一个比例尺（以最长者为准）。`telomere_length` 默认为最长染色体的 2%。输出为 `pdf`、`svg` 或 `eps`（`png` 需要 `rlPyCairo`），且文件必须不存在。
- **API：** `dnakit.visualization.biopython_p19_graphics_functions.chromosome_ideogram(chromosomes, output_path, *, title='Chromosomes', image_format='pdf', telomere_length=None, page_width_cm=24.0, page_height_cm=18.0)`。
- **输入：** `chromosomes`（`name`、`length`，`features` 含 `start`、`end`、`strand`、`label`、`color`）、新的 `output_path`、`title`、`image_format`、`telomere_length`、`page_width_cm` 和 `page_height_cm`。
- **输出：** `output` 路径、`image_format`、`bytes`、`telomere_length` 以及所绘的 `chromosomes`。
- **示例代码：**

```python
import os, tempfile
from dnakit.visualization.biopython_p19_graphics_functions import chromosome_ideogram

chromosomes = [
    {"name": "chr1", "length": 1000, "features": [
        {"start": 100, "end": 200, "strand": 1, "label": "a", "color": "red"},
        {"start": 400, "end": 420, "strand": -1, "label": "b"}]},
    {"name": "chr2", "length": 700},
]
target = os.path.join(tempfile.mkdtemp(), "ideogram.svg")
result = chromosome_ideogram(chromosomes, target, title="Demo", image_format="svg")
print(result["image_format"], result["bytes"] > 1000, result["telomere_length"], [c["name"] for c in result["chromosomes"]])
text = open(target).read()
print("chr1" in text, "Demo" in text)
```

- **示例输出：**

```text
svg True 20 ['chr1', 'chr2']
True True
```
