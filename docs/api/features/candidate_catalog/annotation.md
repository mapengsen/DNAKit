# 基因模型、功能与通路注释

本类收录 `annotation` 模块的 **18** 个函数；以下说明涵盖其中 14 个（`convert_annotation_coordinates`、`make_intergenic_regions`、`make_promoter_regions`、`merge_gene_models` 的小节另行说明）。这 14 个中有 12 个无需任何数据库即可由 GFF3/GTF 文本、文件或特征映射计算出真实结果：`parse`、`validate`、`check_annotation_consistency`、`build_gene_model`、`index`、`summarize`、`compare`、`export`、`extract`、`annotate`，以及两个查表函数 `assign_gene_ontology` / `assign_pathways`。其中 `annotate` 需要 `options["annotation"]`，`assign_gene_ontology` 和 `assign_pathways` 需要 `options["mapping"]` 或 `options["gene_sets"]`；缺少时返回 `status: "not_implemented"`，因为 DNAKit 不自带注释、Gene Ontology 或通路数据库。`assign_protein_domains` 和 `find_homologous_genes` 始终返回 `status: "not_implemented"`（它们需要 Pfam/InterPro 加 HMM 搜索，或直系同源数据库 / BLAST）。所有坐标均为 1 起点、闭区间。示例输出均为真实运行结果；`input_summary` 简写为 `{...}`。

## 1) annotate

- **作用：** 用所提供注释中与区域重叠的特征为基因组区域注释；若无重叠，则给出最近的特征。

- **处理方法：** 按序列对所选类型的特征按起点建立索引；当 `feature.start <= region.end` 且 `feature.end >= region.start` 时区域与特征重叠，`overlap_bp` 为重叠区间的长度。若无重叠，则报告左侧或右侧最近的特征及其 `distance`（bp；区域起点减特征终点，或特征起点减区域终点）；`relation` 相对于特征链为 `upstream`/`downstream`，无链特征则为 `before`/`after`；距离超过 `max_distance` 的特征不报告。所有坐标均为 1 起点、闭区间。没有 `options["annotation"]` 时函数返回 `status: "not_implemented"`，因为 DNAKit 不自带注释数据库。

- **API：** `dnakit.annotation.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：形如 `"chr1:100-200"` / `"chr1:150"` 的区域字符串（数字中可含逗号）、`{"seqid"|"chrom", "start", "end", "id"}` 映射、单个区域，或带 `regions` 的映射（最多 100,000 个）。`options`：`annotation`（必填：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`））、`feature_types`（默认 `["gene"]`，不区分大小写）、`max_distance`（默认不限）、`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`feature_types` 和 `annotations`，每个区域一行，含 `region`、`overlapping`（特征的 `id`、`name`、`type`、`seqid`、`start`、`end`、`strand`、`overlap_bp`）和 `nearest`（同样字段外加 `distance` 与 `relation`，或 `null`）。没有 `annotation` 时返回 `status: "not_implemented"`，含 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.annotation import annotate

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = annotate(["chr1:1100-1150", "chr1:6000-6100"], options={"annotation": gff3, "feature_types": ["gene"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_types": [
      "gene"
    ],
    "annotations": [
      {
        "region": {
          "id": "chr1:1100-1150",
          "seqid": "chr1",
          "start": 1100,
          "end": 1150
        },
        "overlapping": [
          {
            "id": "g1",
            "name": "GENE1",
            "type": "gene",
            "seqid": "chr1",
            "start": 1000,
            "end": 5000,
            "strand": "+",
            "overlap_bp": 51
          }
        ],
        "nearest": null
      },
      {
        "region": {
          "id": "chr1:6000-6100",
          "seqid": "chr1",
          "start": 6000,
          "end": 6100
        },
        "overlapping": [],
        "nearest": {
          "id": "g1",
          "name": "GENE1",
          "type": "gene",
          "seqid": "chr1",
          "start": 1000,
          "end": 5000,
          "strand": "+",
          "distance": 1000,
          "relation": "downstream"
        }
      }
    ]
  }
}
```

## 2) assign_gene_ontology

- **作用：** 在调用者提供的基因到 GO 的映射中查找基因的 Gene Ontology 术语；不自带 GO 数据库。

- **处理方法：** 这是查表而不是富集或预测：术语完全来自 `mapping` 和/或 `gene_sets`，二者合并为一张基因 -> 术语表。每个术语必须形如 `GO:` 加七位数字（否则抛 `InputFormatError`，`INVALID_GO_TERM`）；不会查询 GO 层级结构。没有条目的基因列在 `unannotated` 中。没有 `mapping`/`gene_sets` 时函数返回 `status: "not_implemented"`。

- **API：** `dnakit.annotation.assign_gene_ontology(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：基因 id、以逗号/空白分隔的字符串、基因 id 列表，或带 `genes` 的映射。`options`：`mapping`（`{gene: [terms]}`）和/或 `gene_sets`（`{term: [genes]}`），要得到结果至少需要其一；`term_names`（`{term: name}`，可选标签）、`case_insensitive`（默认 false：基因名区分大小写匹配）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `assignments`（基因 -> 排序后的术语列表）、`terms`（术语 -> `name` 及携带它的 `genes`）、`unannotated`（没有任何术语的基因）和 `annotated_fraction`。没有 `mapping` 和 `gene_sets` 时，封装改为 `status: "not_implemented"`，含 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.annotation import assign_gene_ontology

result = assign_gene_ontology(
    "TP53, BRCA1, UNKNOWN",
    options={
        "mapping": {"TP53": ["GO:0006915", "GO:0003700"], "BRCA1": ["GO:0006281"]},
        "term_names": {"GO:0006281": "DNA repair"},
    },
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "assign_gene_ontology",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "assignments": {
      "TP53": [
        "GO:0003700",
        "GO:0006915"
      ],
      "BRCA1": [
        "GO:0006281"
      ],
      "UNKNOWN": []
    },
    "terms": {
      "GO:0003700": {
        "name": null,
        "genes": [
          "TP53"
        ]
      },
      "GO:0006281": {
        "name": "DNA repair",
        "genes": [
          "BRCA1"
        ]
      },
      "GO:0006915": {
        "name": null,
        "genes": [
          "TP53"
        ]
      }
    },
    "unannotated": [
      "UNKNOWN"
    ],
    "annotated_fraction": 0.6666666666666666
  }
}
```

## 3) assign_pathways

- **作用：** 在调用者提供的基因到通路的映射中查找基因的通路；不自带通路数据库。

- **处理方法：** 与 `assign_gene_ontology` 相同的查表，只是通路标识符为自由格式字符串（如 KEGG 或 Reactome 的 id 或名称），不做格式检查。这不是富集检验。没有 `mapping`/`gene_sets` 时函数返回 `status: "not_implemented"`。

- **API：** `dnakit.annotation.assign_pathways(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：基因 id、以逗号/空白分隔的字符串、基因 id 列表，或带 `genes` 的映射。`options`：`mapping`（`{gene: [terms]}`）和/或 `gene_sets`（`{term: [genes]}`），要得到结果至少需要其一；`term_names`（`{term: name}`，可选标签）、`case_insensitive`（默认 false：基因名区分大小写匹配）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `assignments`（基因 -> 排序后的术语列表）、`terms`（术语 -> `name` 及携带它的 `genes`）、`unannotated`（没有任何术语的基因）和 `annotated_fraction`。没有 `mapping` 和 `gene_sets` 时，封装改为 `status: "not_implemented"`，含 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.annotation import assign_pathways

result = assign_pathways(
    ["TP53", "BRCA1", "UNKNOWN"],
    options={"gene_sets": {"p53 signaling": ["TP53", "MDM2"], "DNA repair": ["BRCA1", "TP53"]}},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "assign_pathways",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "assignments": {
      "TP53": [
        "DNA repair",
        "p53 signaling"
      ],
      "BRCA1": [
        "DNA repair"
      ],
      "UNKNOWN": []
    },
    "terms": {
      "DNA repair": {
        "name": null,
        "genes": [
          "TP53",
          "BRCA1"
        ]
      },
      "p53 signaling": {
        "name": null,
        "genes": [
          "TP53"
        ]
      }
    },
    "unannotated": [
      "UNKNOWN"
    ],
    "annotated_fraction": 0.6666666666666666
  }
}
```

## 4) assign_protein_domains

- **作用：** 未实现：蛋白结构域需要 profile-HMM 或结构域数据库（Pfam、InterPro）以及搜索程序。

- **处理方法：** 不做任何计算；既不自带结构域数据库，也不自带 HMM 搜索。

- **API：** `dnakit.annotation.assign_protein_domains(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“Pfam/InterPro 谱数据库以及 hmmscan 或 InterProScan”）和 `message`（结构域注释需要结构域数据库和 HMM 搜索，二者均未自带），没有 `result`。

- **示例代码：**

```python
import json
from dnakit.annotation import assign_protein_domains

result = assign_protein_domains("MKTAYIAKQRQISFVKSHFSRQ")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "assign_protein_domains",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Pfam/InterPro profile databases and hmmscan or InterProScan",
  "message": "Domain assignment needs a domain database and an HMM search; neither is bundled."
}
```

## 5) build_gene_model

- **作用：** 由 GFF3/GTF 文本或特征映射组装基因模型（基因、转录本、外显子、内含子、CDS、UTR）。

- **处理方法：** 转录本为 mRNA/transcript/`*RNA` 特征，通过 GFF3 的 `Parent` 或 GTF 的 `transcript_id` 与其外显子和 CDS 关联。内含子是相邻外显子之间的间隙。若转录本有 CDS（GTF 还需有起始/终止密码子）特征但没有注释的 UTR，则把 CDS 跨度之外的外显子部分推导为 5' 与 3' UTR，并令 `utr_derived` 为 true（GTF 中位于 CDS 之外的终止密码子也会考虑）。当 CDS 长度能被 3 整除时 `protein_length_aa` 为 `cds_length / 3`，否则为 `null`。区间为 `[start, end]`，1 起点闭区间，按基因组升序排列。不读取序列，因此蛋白长度仅是长度。

- **API：** `dnakit.annotation.build_gene_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`gene_count`、`transcript_count` 和 `genes`，每个基因含 `gene_id`、`name`、`type`、`seqid`、`start`、`end`、`strand`、`biotype` 和 `transcripts`（`transcript_id`、`type`、`seqid`、`start`、`end`、`strand`、`name`、`exons`、`introns`、`cds`、`five_prime_utr`、`three_prime_utr`、`utr_derived`、`exon_count`、`spliced_length`、`cds_length`、`protein_length_aa`）。

- **示例代码：**

```python
import json
from dnakit.annotation import build_gene_model

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = build_gene_model(gff3)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "build_gene_model",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "gene_count": 1,
    "transcript_count": 1,
    "genes": [
      {
        "gene_id": "g1",
        "name": "GENE1",
        "type": "gene",
        "seqid": "chr1",
        "start": 1000,
        "end": 5000,
        "strand": "+",
        "biotype": null,
        "transcripts": [
          {
            "transcript_id": "t1",
            "type": "mRNA",
            "seqid": "chr1",
            "start": 1000,
            "end": 5000,
            "strand": "+",
            "name": null,
            "exons": [
              [
                1000,
                1200
              ],
              [
                2000,
                2300
              ]
            ],
            "introns": [
              [
                1201,
                1999
              ]
            ],
            "cds": [
              [
                1101,
                1200
              ],
              [
                2000,
                2199
              ]
            ],
            "five_prime_utr": [
              [
                1000,
                1100
              ]
            ],
            "three_prime_utr": [
              [
                2200,
                2300
              ]
            ],
            "utr_derived": true,
            "exon_count": 2,
            "spliced_length": 502,
            "cds_length": 300,
            "protein_length_aa": 100
          }
        ]
      }
    ]
  }
}
```

## 6) check_annotation_consistency

- **作用：** 检查基因注释的逻辑一致性：子特征在父特征之内、外显子与 CDS 的几何关系以及 CDS 相位。

- **处理方法：** 规则（代码，严重程度）：`DUPLICATE_ID`（GFF3，错误）、`MISSING_PARENT`（错误）、`CHILD_OUTSIDE_PARENT`（转录本在基因之外、外显子/CDS/UTR 在转录本之外；错误）、子父特征间的 `SEQID_MISMATCH` 与 `STRAND_MISMATCH`（错误）、`EXON_OVERLAP`（同一转录本的外显子重叠；错误）、`CDS_OUTSIDE_EXON`（某 CDS 区间未包含于单个外显子；错误）、`CDS_PHASE_MISMATCH`（相位与按转录本顺序由之前 CDS 片段长度推出的相位不符，相对首个相位为 `(3 - 之前长度 mod 3) mod 3`；错误）、`CDS_LENGTH_NOT_MULTIPLE_OF_3`（CDS 加终止密码子的长度；警告）、`GENE_WITHOUT_TRANSCRIPT` 与 `TRANSCRIPT_WITHOUT_EXON`（警告）以及 `BEYOND_SEQUENCE_END`（需要 `contig_lengths`；错误）。它只检查内部一致性，不判断注释在生物学上是否正确。语法由 `validate` 检查。

- **API：** `dnakit.annotation.check_annotation_consistency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）、`contig_lengths`（`{seqid: length}`，启用 `BEYOND_SEQUENCE_END`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `consistent`（无错误）、`issues`（每项 `{severity, code, line, message}`）、`error_count`、`warning_count` 和 `checked`（`format`、`features`、`genes`、`transcripts`）。示例中第二个 CDS 的相位是 0，而按推算应为 2。

- **示例代码：**

```python
import json
from dnakit.annotation import check_annotation_consistency

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t0\tID=c2;Parent=t1\n"
)
result = check_annotation_consistency(gff3)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "check_annotation_consistency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "consistent": false,
    "issues": [
      {
        "severity": "error",
        "code": "CDS_PHASE_MISMATCH",
        "line": 7,
        "message": "transcript 't1': CDS 2000-2199 has phase 0, expected 2"
      }
    ],
    "error_count": 1,
    "warning_count": 0,
    "checked": {
      "format": "gff3",
      "features": 6,
      "genes": 1,
      "transcripts": 1
    }
  }
}
```

## 7) compare

- **作用：** 按特征坐标 `(seqid, type, start, end, strand)` 比较两份注释。

- **处理方法：** 每份注释都转为 `(seqid, type, start, end, strand)` 键的多重集；`shared`、`only_left`、`only_right` 统计交集和各差集中不同键的个数，`jaccard` 为 `shared / union`（两个空集时取 1.0）。属性、id 和分数都被忽略，因此改名但位置相同的特征视为共有。`by_type` 按特征类型拆分计数。

- **API：** `dnakit.annotation.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：第一份注释（GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`））；第二份在 `options["other"]`。也接受 `[left, right]` 对或映射 `{"left", "right"}`；缺少第二份注释会抛 `ConfigurationError`（`OTHER_REQUIRED`）。`options`：`other`、`feature_types`（限定这些类型，不区分大小写）、`max_examples`（默认 20）、`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `left_count`、`right_count`、`shared`、`only_left`、`only_right`、`jaccard`、`by_type` 以及样例列表 `only_left_examples` / `only_right_examples`（`seqid`、`type`、`start`、`end`、`strand`）。

- **示例代码：**

```python
import json
from dnakit.annotation import compare

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
other = gff3.replace("1000\t1200", "1000\t1250", 1)
result = compare(gff3, options={"other": other, "feature_types": ["exon", "CDS"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "left_count": 4,
    "right_count": 4,
    "shared": 3,
    "only_left": 1,
    "only_right": 1,
    "jaccard": 0.6,
    "by_type": {
      "CDS": {
        "shared": 2,
        "only_left": 0,
        "only_right": 0
      },
      "exon": {
        "shared": 1,
        "only_left": 1,
        "only_right": 1
      }
    },
    "only_left_examples": [
      {
        "seqid": "chr1",
        "type": "exon",
        "start": 1000,
        "end": 1200,
        "strand": "+"
      }
    ],
    "only_right_examples": [
      {
        "seqid": "chr1",
        "type": "exon",
        "start": 1000,
        "end": 1250,
        "strand": "+"
      }
    ]
  }
}
```

## 8) convert_annotation_coordinates

- **作用：** 为“基因模型、功能与通路注释”类分析提供可复用的 `convert_annotation_coordinates` 功能。候选编号 CAP-0042，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `convert_annotation_coordinates` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/annotation/local.py`，并通过正式的 DNAKit 命名空间导出。坐标约定、contig 长度和合并规则都通过显式选项给出；不查询注释数据库或参考基因组。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.annotation.convert_annotation_coordinates(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `convert_annotation_coordinates` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的坐标来源信息。它只做确定性的区间运算，不查询注释数据库或参考序列。

- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import convert_annotation_coordinates

result = convert_annotation_coordinates(
    input={"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
    options={"coordinate_system": "1-based-closed", "target_coordinate_system": "0-based-half-open"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0042",
  "operation": "convert_annotation_coordinates",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'convert_annotation_coordinates', 'result': {'records': [{'chrom': 'chr1', 'start': 1, 'end': 5, 'id': 'g1', 'coordinate_system': '0-based-half-open'}], 'record_count': 1, 'source_coordinate_system': '1-based-closed', 'target_coordinate_system': '0-based-half-...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) export

- **作用：** 把注释转换为 GFF3、GTF、BED6、BED12 或 JSON 文本，可选写入文件。

- **处理方法：** GFF3 转 GTF 时通过 `Parent` 链解析 `gene_id`/`transcript_id`；GTF 转 GFF3 时由子特征补出缺失的基因和转录本行（`synthesize_parents`）。BED 使用 0 起点半开起点（`start - 1`）；BED6 导出某一类型的特征；BED12 每个转录本一行（区块为外显子，粗线区域为 CDS 加起始/终止密码子，非编码时为空）。文件为原子写入，且仅在给出 `output` 时才写。

- **API：** `dnakit.annotation.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）。`options`：`format`（目标格式：默认 `"gff3"`，或 `"gtf"`、`"bed"`、`"bed12"`、`"json"`）、`source_format`（输入格式：`"auto"`、`"gff3"`、`"gtf"`）、`feature_type`（仅 BED6，默认 `"gene"`）、`synthesize_parents`（GTF 转 GFF3，默认 true）、`output`（文件路径；不给则不写）、`overwrite`（默认 false；文件已存在会抛 `ConfigurationError`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`text`、`line_count`（数据行数，不含注释）和 `path`（已写入的文件或 `null`）。

- **示例代码：**

```python
import json
from dnakit.annotation import export

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = export(gff3, options={"format": "bed12"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "bed12",
    "text": "chr1\t999\t2300\tt1\t0\t+\t1100\t2199\t0\t2\t201,301,\t0,1000,\n",
    "line_count": 1,
    "path": null
  }
}
```

## 10) extract

- **作用：** 从 GFF3/GTF 注释中提取某一生物学实体的特征，或在 DNA 序列中扫描开放阅读框。

- **处理方法：** `entity`（默认 `gene`）可取 `gene`、`transcript`、`mrna`、`exon`、`cds`、`start_codon`、`stop_codon`、`mirna`（miRNA/pre_miRNA 类型，或具有 miRNA 生物型的基因/转录本）、`intron`（推导为相邻外显子之间的间隙，并按转录本顺序编号）、`utr`/`five_prime_utr`/`three_prime_utr`（已注释的，或像 `build_gene_model` 那样由外显子和 CDS 跨度推导）、其他任意值（不区分大小写地匹配特征类型）或 `orf`。给出 `genome` 时，会为每个特征添加其基因组序列（负链则取反向互补）。`orf` 在所有阅读框中扫描（除非 `both_strands` 为 false，否则含双链）：ORF 从终止密码子之后的第一个起始密码子起，到下一个同框终止密码子止，使用标准遗传密码；坐标为正链上的 1 起点坐标。ORF 仅按此定义报告，除 `min_length` 外不做基因可能性过滤。

- **API：** `dnakit.annotation.extract(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：注释类实体为GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）；`orf` 为 DNA 序列、FASTA 文本、映射 `{id: sequence}` 或列表。`options`：`entity`、`seqid`、`strand`、`start`/`end`（保留与该范围重叠的特征）、`limit`（默认 100000）、`genome`（FASTA 文本/路径或 `{seqid: sequence}`）、`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）；`orf` 另有 `min_length`（含终止密码子的核苷酸数，默认 90）、`start_codons`（默认 `["ATG"]`）、`stop_codons`（默认 `["TAA", "TAG", "TGA"]`）、`require_stop`（默认 true）、`both_strands`（默认 true）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `entity`、`format`、`count`、`truncated` 和 `features`（特征记录，含 `seqid`、`type`、`start`、`end`、`strand`、`phase`、`attributes`、`derived`，给出 `genome` 时还有 `sequence`）；`orf` 的 `result` 含 `entity`、`count` 和 `features`（各 ORF：`seqid`、`strand`、`frame`（+1..+3 / -1..-3）、`start`、`end`、含终止密码子的核苷酸 `length`、`has_stop`、`protein` 和 `protein_length`）。

- **示例代码：**

```python
import json
from dnakit.annotation import extract

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = extract(gff3, options={"entity": "intron"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "extract",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "intron",
    "format": "gff3",
    "count": 1,
    "truncated": false,
    "features": [
      {
        "line": null,
        "seqid": "chr1",
        "source": "dnakit",
        "type": "intron",
        "start": 1201,
        "end": 1999,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "Parent": "t1",
          "intron_number": 1
        },
        "derived": true
      }
    ]
  }
}
```

## 11) find_homologous_genes

- **作用：** 未实现：同源搜索需要直系同源数据库或序列相似性搜索程序。

- **处理方法：** 不做任何计算；既不自带参考数据（Ensembl Compara、OrthoDB），也不自带 BLAST/DIAMOND 及参考蛋白质组。

- **API：** `dnakit.annotation.find_homologous_genes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“直系同源数据库（Ensembl Compara、OrthoDB）或带参考蛋白质组的 BLAST/DIAMOND”）和 `message`（查找同源基因需要参考数据或相似性搜索，二者均未自带），没有 `result`。

- **示例代码：**

```python
import json
from dnakit.annotation import find_homologous_genes

result = find_homologous_genes("BRCA1", options={"species": "mouse"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "find_homologous_genes",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an orthology database (Ensembl Compara, OrthoDB) or BLAST/DIAMOND with reference proteomes",
  "message": "Finding homologous genes needs reference data or a similarity search; neither is bundled."
}
```

## 12) index

- **作用：** 为 GFF3/GTF 注释建立查找表：按 id、按序列（排序的区间）以及父到子。

- **处理方法：** 对特征做一次遍历。`by_id` 保留每个 id 的第一个特征；`by_seqid` 列出按起点、终点排序的 `[start, end, id_or_null, type]`；`children` 把每个父 id 映射到以它为父的特征 id（按文件顺序；没有 id 的特征如 exon/CDS 行以 `"<type>@<line>"` 表示）；`gene_transcripts` 把基因 id 映射到转录本 id。这些是普通字典，不是区间树。

- **API：** `dnakit.annotation.index(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`feature_count`、`by_id`、`by_seqid`、`children`、`gene_transcripts` 和 `type_counts`。

- **示例代码：**

```python
import json
from dnakit.annotation import index

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
)
result = index(gff3)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "index",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 2,
    "by_id": {
      "g1": {
        "seqid": "chr1",
        "type": "gene",
        "start": 1000,
        "end": 5000,
        "strand": "+"
      },
      "t1": {
        "seqid": "chr1",
        "type": "mRNA",
        "start": 1000,
        "end": 5000,
        "strand": "+"
      }
    },
    "by_seqid": {
      "chr1": [
        [
          1000,
          5000,
          "g1",
          "gene"
        ],
        [
          1000,
          5000,
          "t1",
          "mRNA"
        ]
      ]
    },
    "children": {
      "g1": [
        "t1"
      ]
    },
    "gene_transcripts": {
      "g1": [
        "t1"
      ]
    },
    "type_counts": {
      "gene": 1,
      "mRNA": 1
    }
  }
}
```

## 13) make_intergenic_regions

- **作用：** 为“基因模型、功能与通路注释”类分析提供可复用的 `make_intergenic_regions` 功能。候选编号 CAP-0047，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `make_intergenic_regions` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/annotation/local.py`，并通过正式的 DNAKit 命名空间导出。坐标约定、contig 长度和合并规则都通过显式选项给出；不查询注释数据库或参考基因组。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.annotation.make_intergenic_regions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `make_intergenic_regions` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的坐标来源信息。它只做确定性的区间运算，不查询注释数据库或参考序列。

- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import make_intergenic_regions

result = make_intergenic_regions(
    input={"genes": [{"chrom": "chr1", "start": 2, "end": 4}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0047",
  "operation": "make_intergenic_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_intergenic_regions', 'result': {'regions': [{'id': 'intergenic_chr1_1', 'chrom': 'chr1', 'start': 1, 'end': 1, 'feature_type': 'intergenic', 'coordinate_system': '1-based-closed'}, {'id': 'intergenic_chr1_2', 'chrom': 'chr1', 'start': 5, 'end': 10, 'feat...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) make_promoter_regions

- **作用：** 为“基因模型、功能与通路注释”类分析提供可复用的 `make_promoter_regions` 功能。候选编号 CAP-0048，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `make_promoter_regions` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/annotation/local.py`，并通过正式的 DNAKit 命名空间导出。坐标约定、contig 长度和合并规则都通过显式选项给出；不查询注释数据库或参考基因组。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.annotation.make_promoter_regions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `make_promoter_regions` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的坐标来源信息。它只做确定性的区间运算，不查询注释数据库或参考序列。

- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import make_promoter_regions

result = make_promoter_regions(
    input={"genes": [{"id": "g1", "chrom": "chr1", "start": 3, "end": 6, "strand": "+"}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0048",
  "operation": "make_promoter_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_promoter_regions', 'result': {'promoters': [{'id': 'promoter_g1', 'gene_id': 'g1', 'chrom': 'chr1', 'start': 1, 'end': 3, 'strand': '+', 'feature_type': 'promoter', 'coordinate_system': '1-based-closed', 'upstream': 2, 'downstream': 1}], 'promoter_count'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) merge_gene_models

- **作用：** 为“基因模型、功能与通路注释”类分析提供可复用的 `merge_gene_models` 功能。候选编号 CAP-0049，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge_gene_models` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/annotation/local.py`，并通过正式的 DNAKit 命名空间导出。坐标约定、contig 长度和合并规则都通过显式选项给出；不查询注释数据库或参考基因组。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.annotation.merge_gene_models(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `merge_gene_models` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的坐标来源信息。它只做确定性的区间运算，不查询注释数据库或参考序列。

- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import merge_gene_models

result = merge_gene_models(
    input={"segments": [{"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3}, {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5}]} ,
    options={"coordinate_system": "1-based-closed", "merge_adjacent": True}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0049",
  "operation": "merge_gene_models",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'merge_gene_models', 'result': {'models': [{'gene_id': 'g1', 'gene_name': 'g1', 'chrom': 'chr1', 'strand': '+', 'start': 1, 'end': 5, 'segments': [{'start': 1, 'end': 5}], 'segment_count': 1, 'coordinate_system': '1-based-closed'}], 'model_count': 1, 'input_se...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) parse

- **作用：** 把 GFF3 或 GTF 文本（或 `.gff`/`.gtf`/`.gz` 文件）解析为特征记录。

- **处理方法：** 格式由 `##gff-version` 行或属性风格自动判定，除非指定 `format`。GFF3 属性做百分号解码，逗号列表变为列表；GTF 属性为 `key "value";` 对，重复键变为列表。坐标保持 1 起点闭区间。启用 `strict`（默认）时，格式错误的行会抛 `InputFormatError`（`INVALID_ANNOTATION`）；`strict: false` 时跳过该行并列入 `skipped`。特征映射列表只做规范化而不解析。

- **API：** `dnakit.annotation.parse(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：注释文本、文件路径（允许 `.gz`），或特征映射列表。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）、`strict`（默认 true）、`limit`（返回特征的最大数，默认 100000）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`feature_count`、`truncated`、`features`（每个含 `line`、`seqid`、`source`、`type`、`start`、`end`、`score`、`strand`、`phase`、`attributes`）、`directives`（`##` 行）和 `skipped`。

- **示例代码：**

```python
import json
from dnakit.annotation import parse

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
)
result = parse(gff3)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "parse",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 2,
    "truncated": false,
    "features": [
      {
        "line": 2,
        "seqid": "chr1",
        "source": "src",
        "type": "gene",
        "start": 1000,
        "end": 5000,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "ID": "g1",
          "Name": "GENE1"
        }
      },
      {
        "line": 3,
        "seqid": "chr1",
        "source": "src",
        "type": "mRNA",
        "start": 1000,
        "end": 5000,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "ID": "t1",
          "Parent": "g1"
        }
      }
    ],
    "directives": [
      "##gff-version 3"
    ],
    "skipped": []
  }
}
```

## 17) summarize

- **作用：** 汇总 GFF3/GTF 注释：按类型、序列和链统计特征数量与长度，以及基因结构统计量。

- **处理方法：** 长度为 `end - start + 1`。基因结构沿用 `build_gene_model` 的模型：`intron_count` 由外显子之间的间隙推导，`coding_transcripts` 统计含 CDS 的转录本，`multi_exon_transcripts` 统计外显子多于一个的转录本，`transcripts_per_gene` / `exons_per_transcript` 给出均值和最大值（无数据时为 0）。

- **API：** `dnakit.annotation.summarize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：GFF3 或 GTF 文本、GFF3/GTF 文件路径（允许 `.gz`），或特征映射列表（形如 `parse` 的返回；每个需要 `seqid`、整数 `start` <= `end` 且 >= 1，可选 `type`、`strand`、`phase`、`attributes`）。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`feature_count`、`types`（类型 -> `count`、`total_bp`、`mean_length`）、`seqids`、`strands`、`gene_count`、`transcript_count`、`exon_count`、`intron_count`、`coding_transcripts`、`multi_exon_transcripts`、`transcripts_per_gene` 和 `exons_per_transcript`（`mean`、`max`）。

- **示例代码：**

```python
import json
from dnakit.annotation import summarize

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = summarize(gff3)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "summarize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 6,
    "types": {
      "CDS": {
        "count": 2,
        "total_bp": 300,
        "mean_length": 150.0
      },
      "exon": {
        "count": 2,
        "total_bp": 502,
        "mean_length": 251.0
      },
      "gene": {
        "count": 1,
        "total_bp": 4001,
        "mean_length": 4001.0
      },
      "mRNA": {
        "count": 1,
        "total_bp": 4001,
        "mean_length": 4001.0
      }
    },
    "seqids": {
      "chr1": 6
    },
    "strands": {
      "+": 6
    },
    "gene_count": 1,
    "transcript_count": 1,
    "exon_count": 2,
    "intron_count": 1,
    "coding_transcripts": 1,
    "multi_exon_transcripts": 1,
    "transcripts_per_gene": {
      "mean": 1.0,
      "max": 1
    },
    "exons_per_transcript": {
      "mean": 2.0,
      "max": 2
    }
  }
}
```

## 18) validate

- **作用：** 校验 GFF3 或 GTF 文件的语法：列数、坐标、链、相位、属性、id 与父特征。

- **处理方法：** 错误：列数不对、坐标非整数或非正、end < start、链不在 `+ - . ?` 中、相位不在 0/1/2/`.` 中、CDS 缺少相位、属性无法解析；GFF3 另有非 CDS 特征的 `ID` 重复和 `Parent` 指向未定义的 `ID`；GTF 另有缺少 `gene_id`（所有行）或 `transcript_id`（除基因行外的所有行）。警告：GFF3 缺少 `##gff-version 3` 指令。格式错误的行只作报告而不抛出。基因结构的逻辑检查见 `check_annotation_consistency`。

- **API：** `dnakit.annotation.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：注释文本或文件路径（允许 `.gz`）。`options`：`format`（默认 `"auto"`：由 `##gff-version` 行或属性风格判定；`"gff3"`、`"gtf"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `valid`（无错误）、`format`、`feature_count`（成功解析的行数）、`errors` 和 `warnings`（每项 `{line, severity, code, message}`，按行号排序）。

- **示例代码：**

```python
import json
from dnakit.annotation import validate

bad = "chr1\tsrc\tgene\t50\t10\t.\t+\t.\tID=g1\nchr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=nope\n"
result = validate(bad)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "annotation",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "format": "gff3",
    "feature_count": 1,
    "errors": [
      {
        "line": 1,
        "severity": "error",
        "code": "COORDINATE_ORDER",
        "message": "end (10) is smaller than start (50)"
      },
      {
        "line": 2,
        "severity": "error",
        "code": "UNDEFINED_PARENT",
        "message": "Parent 'nope' is not defined"
      }
    ],
    "warnings": [
      {
        "line": 1,
        "severity": "warning",
        "code": "MISSING_VERSION",
        "message": "no '##gff-version 3' directive"
      }
    ]
  }
}
```
