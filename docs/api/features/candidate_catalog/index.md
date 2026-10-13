# DNAKit 837 候选功能设计文档

本目录收录 **837** 项候选：410 个 `new_candidate` 和 427 个 `review_needed`。
这些文档记录正式 DNAKit 命名空间、Agent/MCP JSON 调用计划和本地适配层状态；适配层返回依赖状态时不等于已经完成底层科学后端验证。
不同输入类型、后端、格式、平台、提供方和生物学对象统一作为参数维度记录。

## Categories

| Category | Candidates | Document |
|---|---:|---|
| Alignment | 34 | [alignment.md](alignment.md) |
| Annotation | 18 | [annotation.md](annotation.md) |
| Assembly | 37 | [assembly.md](assembly.md) |
| Breeding | 7 | [breeding.md](breeding.md) |
| Clinical genomics | 18 | [clinical.md](clinical.md) |
| Databases | 8 | [databases.md](databases.md) |
| CRISPR design | 19 | [design-crispr.md](design-crispr.md) |
| Epigenomics | 13 | [epigenome.md](epigenome.md) |
| Experimental design | 9 | [experimental-design.md](experimental-design.md) |
| External adapters | 9 | [external-adapters.md](external-adapters.md) |
| FASTA and FASTQ | 66 | [fasta-fastq.md](fasta-fastq.md) |
| Forensic and reproductive genomics | 9 | [forensic-reproductive.md](forensic-reproductive.md) |
| Three-dimensional genome | 9 | [genome3d.md](genome3d.md) |
| GWAS | 29 | [gwas.md](gwas.md) |
| Haplotypes | 17 | [haplotypes.md](haplotypes.md) |
| Input and output formats | 131 | [io-formats.md](io-formats.md) |
| Knowledge graph | 18 | [knowledge-graph.md](knowledge-graph.md) |
| Microbiome | 21 | [microbiome.md](microbiome.md) |
| Machine learning models | 48 | [ml-models.md](ml-models.md) |
| Molecular laboratory | 10 | [molecular-lab.md](molecular-lab.md) |
| Phylogeny | 8 | [phylogeny.md](phylogeny.md) |
| Population genetics | 8 | [population.md](population.md) |
| Regulation | 11 | [regulation.md](regulation.md) |
| Replication and DNA repair | 10 | [replication-repair.md](replication-repair.md) |
| Sequence analysis | 52 | [sequence.md](sequence.md) |
| Sequencing quality control | 50 | [sequencing-qc.md](sequencing-qc.md) |
| Single-cell and spatial genomics | 54 | [singlecell-spatial.md](singlecell-spatial.md) |
| Transcriptomics | 18 | [transcriptomics.md](transcriptomics.md) |
| Variant analysis | 63 | [variants.md](variants.md) |
| Visualization | 6 | [visualization.md](visualization.md) |
| Workflow and data security | 27 | [workflow-security.md](workflow-security.md) |

## Status counts

- `new_candidate`: 410
- `review_needed`: 427
- `implemented_candidate`: 189
- `partial_candidate`: 5
- `native_facade`: 4
- `adapter_candidate`: 639
- `pending`: 0

## Required implementation follow-up

For production scientific execution, define the concrete backend dependencies, error contract, side effects, tests, resource limits, and MCP exposure policy. Existing public functions and Skills remain the source of workflow-specific backend guidance.
