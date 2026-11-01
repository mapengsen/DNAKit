# 公开 DNA 工具覆盖

`new_vision/2026-10-10/public_dna_tools.md` 是一份标注日期的来源调研。DNAKit
现在把逐项目的对照记录保存在
`new_vision/2026-10-10/public_dna_tools_coverage.{json,csv,md}`，并通过
`dnakit.integrations` 提供同样的记录。

> `new_vision/` 是作者本地保留的审计目录，不在仓库中，所以下面提到的路径在代码检出中找不到。结果已汇总在本页，也可以通过 `dnakit.integrations` 获取。

来源文本声明有 282 条表格记录。仓库中可用的 Markdown 共有 425 处明确提及，
归并后为 320 条规范化来源记录。调研本身声明 282 条；重复的分类、引用和
特定提供方列表作为独立证据保留，拼写别名则保守合并。报告还单独记录了来源中
声明的 169 条有条件可集成记录。

每个纳入范围的候选能力都使用正式的 `dnakit.<module>.<function>()` 公开入口。
`dnakit_candidates` 只是实现登记使用的内部暂存命名空间，不是公开 API。正式入口
可能是本地实现、有限门面或带依赖检查的适配器，具体状态在下文单独记录。

## P00–P45 操作批次接入

46 批清单共 3,512 项操作：1,473 项实现、866 项通过参数合并、1,043 项已有对应、130 项保留跳过原因。清单声明的 2,041 个函数已公开导出并可由 Agent/MCP 发现；另公开 10 个被历史记录引用的已有函数。操作行数与函数数口径不同。

这次接入新增 2,051 个公开入口，当前共有 5,795 个公开 Python 函数，其中 5,792 个兼容 Agent。既有同名入口保留；新增 `dnakit.sequence.biopython_trim_sequence` 与 `dnakit.databases.ensembl_vep_region_query` 用于避免覆盖原 API。可选后端按需安装，外部程序、数据库与模型仍需各自准备。Biotite 2/pyranges1 1.5 要求 Python 3.12，旧 TensorFlow Addons 后端需独立的 Python 3.10/3.11 环境。

原始 `dedup_clusters_all.csv` 和 `dedup_summary.json` 继续作为历史快照，不代表当前公开接口数量。逐项接入与跳过记录位于仓库 `paper/write_paper/figure/fig2/task_coverage/integration/`，全量验证命令见该目录的 `PUBLIC_INTEGRATION_REPORT.md`。

## 查询覆盖矩阵

```python
from dnakit.integrations import (
    describe_third_party_tool,
    list_third_party_tools,
    list_unified_capabilities,
    integration_inventory,
    run_unified_capability,
    third_party_coverage,
    third_party_route,
    unified_capability_summary,
)

print(third_party_coverage()["summary"])
print(integration_inventory())
print([item["name"] for item in list_third_party_tools("biopython")])
print(describe_third_party_tool("Biopython")["citation_dois"])
print(third_party_route("Cas-OFFinder")["api"])
print(unified_capability_summary()["count"])
print(list_unified_capabilities(module="sequence")[0]["api"])
print(run_unified_capability("CAP-0187", {"input": "@r1\\nACGT\\n+\\nIIII\\n"}))
```

路由是刻意显式给出的：

- `analysis_adapter` 对应 `dnakit.integrations.run_analysis`，仍需要单独安装
  可执行程序和参考数据。
- `mcp_registry` 对应 MCP 桥接，仍需要配置好的本地服务器或作者文档中给出的
  服务端点。
- `skill_registry` 对应随包提供的 Skill／工作流层；Skill 自身的许可证、模型和
  数据条件依然适用。
- `pip_adapter` 对应显式的可选 Python 后端，例如 Biopython、pydna、seqme、
  DNA Chisel、Biotite、Edlib、PyRanges、pyfastx、dnaio、mappy、pyBigWig、pybedtools、
  gffutils、pySBOL3、scadnano、polars-bio、DnaCauldron、DNA Features Viewer、gtars、
  GeneBlocks 或 sourmash；对应的包和模型 extra 需要另行安装。
- `native_facade` 把调研中的工具对应到聚焦的 DNAKit 实现，例如
  BWA/Bowtie2/SAMtools/BCFtools/SPAdes/Picard NGS 门面、NCBI
  BLAST+/minimap2/BEDTools 比对门面、MACS3 peak 门面、SeqKit FASTA/FASTQ 函数、
  FastTree 建树门面、SPOA 共识、vcfanno/SURVIVOR 变异门面，以及
  NanoFilt/Filtlong/NanoPlot/Porechop/Trimmomatic/AfterQC/SeqMagick/GATK4/Shasta/Merqury/BUSCO/DeepVariant/cuteSV/Sniffles2/DELLY/HTSlib/seqtk/WhatsHap/BamTools/VSEARCH/BEDOPS/Flye/hifiasm/MEGAHIT/DIAMOND/MMseqs2/Kraken2/Prokka/Bakta/CheckM2/Bracken/Centrifuge/Kaiju/Unicycler/GTDB-Tk/HUMAnN/MethylDackel/GEMINI/Ensembl VEP
  读段处理门面；这并不表示已经重新实现了上游命令的完整命令行。
- `project_evidence_only` 和 `no_local_match` 是发现记录。它们只作为显式的
  外部适配边界路由到 `dnakit.external_adapters.run`，并不表示已复制上游源码，
  也不表示其全部功能都已实现。

`local_match` 表示在 DNAKit 中找到了名称、仓库、依赖、Skill 声明、MCP 条目或
分析入口，并不是科学验证结果。

调用 `integration_inventory()` 可获取当前各接口面的数量。它把 Python 函数、
兼容 Agent 的 Python 函数、MCP 服务器配置、已注册 Skills、调研项目和候选
Schema 作为不同单位分别统计；这些类别不能当作可互换的工具直接相加。
当前原始接口面合计为 **6,419**（5,795 个公开 Python 函数 + 49 个 MCP 服务器
配置 + 575 个已注册 Skills）；这是接口面计数，而不是去重后的生物学工具数，
因为同一能力可能出现在多个接口面中。
清单还包含 `unified_candidate_statuses`，区分已实现、部分实现和适配层候选。
它还保留了调研的 `records_by_kind` 和集成状态计数，避免把调研中的 103 条
MCP 记录误认为 103 个已在本地连接的服务器。
在此快照中，320 条调研记录包括 62 条 `native_facade`、2 条
`analysis_adapter`、26 条 `pip_adapter`、13 条 `mcp_registry`、2 条 `skill_registry`、
44 条 `project_evidence_only` 和 171 条 `no_local_match`。

随包提供的 `unified_capabilities.json` 还记录了 837 个正式候选 API 名称及其
输入／输出 JSON Schema。`list_unified_capabilities()` 和
`describe_unified_capability()` 通过同一套 Python 与 MCP 发现层提供这份目录。
目前有 32 条项目记录与一个或多个能力存在保守的名称／Schema 匹配；这些只是
待审核的关联，并不证明上游项目的完整 API 已经实现。

`run_unified_capability()` 是目录能力的执行桥接。它通过现有的 Agent/MCP
注册表解析正式的 `dnakit.<module>.<operation>` 入口，应用其 JSON 参数 Schema
和序列化限制；涉及写文件、下载、外部程序或模型操作时，必须显式传入
`allow_side_effects=True`。它不会安装上游包，也不会在后台执行未经验证的项目。

调研中采用宽松许可证的 pip 工具都有显式的受限适配器：

- 现有的 `dnakit.bioinformatics` 格式／比对函数在安装 `dnakit[bio-formats]`
  后使用 Biopython，读段函数在安装 `dnakit[bio-reads]` 后使用 pysam，克隆
  函数在安装 `dnakit[bio-assembly]` 后使用 pydna。这些是聚焦的 DNAKit 门面，
  不是对上游完整命名空间的重导出。

- `dnakit.integrations.beta_diversity_scikit_bio()` 和
  `dnakit.integrations.permanova_scikit_bio()` 提供受限的 scikit-bio 统计。
  安装 `dnakit[bio-skibio]`；不会重导出完整的 scikit-bio 命名空间。
- `dnakit.integrations.evaluate_seqme()` 构造白名单内的 `seqme.metrics` 类，
  并序列化评价结果表。安装 `dnakit[bio-seqme]`；模型 extra 和权重需另行处理。
- `dnakit.integrations.optimize_with_dnachisel()` 构造白名单内的一组 DNA Chisel
  约束规范并返回优化后的序列。安装 `dnakit[bio-dnachisel]`；不复制 DNA Chisel 源码。
- `dnakit.integrations.fetch_pyfaidx_region()` 提供有界 FASTA 区间读取，使用
  `dnakit[bio-pyfaidx]`；`summarize_vcf_vcfpy()` 和 `summarize_vcf_cyvcf2()`
  提供有界 VCF 汇总，分别使用对应的可选 extra。
- `dnakit.integrations.align_biotite()`、`align_edlib()` 和 PyRanges 的两个
  BED 适配器提供有界成对比对、编辑距离和区间操作，分别使用
  `bio-biotite`、`bio-edlib` 与 `bio-pyranges`。
- `dnakit.integrations.summarize_pyfastx()` 汇总有界 FASTA/FASTQ 记录，使用
  `bio-pyfastx`；`summarize_pybigwig_region()` 返回有界本地 bigWig 区域统计，使用
  `bio-pybigwig`。二者均不建立可变索引、不写入输入文件，也不访问远程 URL。
- `dnakit.integrations.summarize_dnaio_fastx()` 使用 `bio-dnaio` 只读汇总有界
  FASTA/FASTQ 流；`dnakit.integrations.map_mappy_sequences()` 使用 `bio-mappy`
  和 minimap2/mappy 对本地参考和 reads 做有界映射。二者均不写索引、不接受任意
  上游命令行参数，缺少可选依赖时返回结构化依赖状态。
- `dnakit.integrations.intersect_pybedtools_bed()` 提供有界 BED 交集，使用
  `bio-pybedtools`；它仍要求用户单独安装 `bedtools` 可执行程序，且不接受任意命令行参数。
- `dnakit.integrations.query_gffutils_database()` 只读查询已经存在的 gffutils
  `FeatureDB`，使用 `bio-gffutils`；不会创建数据库或执行任意 SQL。
- `dnakit.integrations.summarize_sbol3_document()` 读取本地 SBOL3 RDF 文档并返回
  有界对象/类型清单，使用 `bio-pysbol3`；不会解析远程 RDF、修改对象或写文件。
- `dnakit.integrations.summarize_scadnano_design()` 读取本地 scadnano `.sc` DNA
  纳米结构设计并返回有界螺旋/链清单，使用 `bio-scadnano`；不会编辑、渲染或下载设计资产。
- `dnakit.integrations.overlap_polars_bio_bed()` 使用 `bio-polarsbio` 对两个本地
  BED/CSV 区间文件做有界交集，输出左侧区间；不会执行任意 SQL 或访问远程 URL。
- `dnakit.integrations.simulate_dnacauldron_assembly()` 使用 `bio-dnacauldron`
  模拟 Gibson 或 Type IIS 组装，返回预测构建体及后端报告的警告/缺陷；不会推断湿实验结果。
- `dnakit.integrations.summarize_dna_features_viewer()` 使用
  `bio-dna-features-viewer` 只读汇总本地 GenBank/GFF 特征；不会写图像文件，也不重导出
  DNA Features Viewer 或 Biopython 的完整 API。
- `dnakit.integrations.summarize_gtars_regions()` 使用 `bio-gtars` 只读汇总本地 BED
  `RegionSet`，返回区间数、碱基覆盖和染色体计数；不接受远程 URL、缓存构建或任意区间表达式。
- `dnakit.integrations.diff_geneblocks_sequences()` 使用 `bio-geneblocks` 对两个有界 DNA
  字符串做确定性的差异块比较；不会调用 BLAST、绘图或特征转移。
- `dnakit.integrations.compare_sourmash_sequences()` 使用 `bio-sourmash` 在内存中比较两条
  有界 A/C/G/T 序列的 MinHash sketch；不会访问数据库、网络服务或分类索引。
- `dnakit.sequencing_qc.run_cutadapt()` 使用用户提供的 Cutadapt 可执行程序，以固定参数白名单
  清理单端或双端 FASTQ 接头；需要显式输出路径，不安装程序，也不覆盖已有文件。
- `dnakit.sequencing_qc.run_fastp()` 使用用户提供的 fastp 可执行程序执行有界 FASTQ 过滤和报告
  生成；支持单端或双端输入，外部进程需要 Agent/MCP 调用者显式授权。
- `dnakit.sequencing_qc.run_fastqc()` 和 `run_multiqc()` 分别运行用户提供的 FastQC 和 MultiQC
  可执行程序，监控 HTML/ZIP 或汇总 HTML 产物；不会捆绑 Java、插件或 GPL 后端。

截图中列出的 UGENE、iLearnPlus、iFeatureOmega 和 `ga` 均已保留为显式
外部边界；它们可以被引用或由用户单独安装后调用，但当前不会被误报为
DNAKit 核心源码已经复制。`third_party_route()` 会返回对应的外部入口和
许可证/引用证据。

## 引用与许可证证据

调研为以下核心工具提供了论文标识；请引用论文中的 DOI，并在任何分发中保留上游软件声明：

| 工具 | 引用 | DNAKit 边界 |
| --- | --- | --- |
| Biopython | Cock et al. (2009), [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) | `dnakit.bioinformatics` 适配器；适用 Biopython 自身的许可证和第三方声明。 |
| SeqKit | Shen, Sipos & Zhao (2024), SeqKit2, [10.1002/imt2.191](https://doi.org/10.1002/imt2.191); 早期 SeqKit 论文为 [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) | 通过外部分析适配器调用单独安装的命令行工具。 |
| scikit-bio | Aton et al. (2026), [10.1038/s41592-025-02981-z](https://doi.org/10.1038/s41592-025-02981-z) | 可选 Python 依赖；BSD-3-Clause 声明随依赖保留。 |
| UGENE | Okonechnikov et al. (2012), [10.1093/bioinformatics/bts091](https://doi.org/10.1093/bioinformatics/bts091) | 外部 GPLv2 程序；不复制进核心包。 |
| pydna | Pereira et al. (2015), [10.1186/s12859-015-0544-x](https://doi.org/10.1186/s12859-015-0544-x) | 可选 DNA 组装依赖；保留 BSD 声明。 |
| DNA Chisel | Zulkower & Rosser (2020), [10.1093/bioinformatics/btaa558](https://doi.org/10.1093/bioinformatics/btaa558) | 可选 MIT 依赖；模型和数据条款另行适用。 |
| seqme | Møller-Larsen et al. (2025), [arXiv:2511.04239](https://arxiv.org/abs/2511.04239) | 可选 BSD-3-Clause 适配器；模型 extra 和权重另行处理。 |
| Biotite | Kunzmann & Hamacher (2018), [10.1186/s12859-018-2367-z](https://doi.org/10.1186/s12859-018-2367-z) | 可选 BSD-3-Clause 适配器；仅限受限的核苷酸比对。 |
| Edlib | Šošić & Šikić (2017), [10.1093/bioinformatics/btw753](https://doi.org/10.1093/bioinformatics/btw753) | 可选 MIT 适配器；仅限受限的编辑距离比对。 |
| PyRanges | Stovner & Sætrom (2020), [10.1093/bioinformatics/btz615](https://doi.org/10.1093/bioinformatics/btz615) | 可选 MIT 适配器；仅限受限的 BED 区间操作。 |
| pyfastx | Du et al. (2021), [10.1093/bib/bbaa368](https://doi.org/10.1093/bib/bbaa368) | 可选 MIT 适配器；仅限受限的 FASTA/FASTQ 汇总。 |
| pyfaidx | Shirley et al. (2015), [10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1)；源码 BSD-3-Clause | 可选 BSD-3-Clause 适配器；受限的只读 FASTA 区间提取。 |
| dnaio | dnaio 软件存档, [10.5281/zenodo.10548864](https://doi.org/10.5281/zenodo.10548864); MIT | 可选 MIT 适配器；受限的只读 FASTA/FASTQ 流汇总。 |
| mappy/minimap2 | Li (2018), [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191); [官方 Python 绑定 README](https://github.com/lh3/minimap2/blob/master/python/README.rst); MIT | 可选 MIT 适配器；使用白名单预设的受限本地映射。 |
| pyBigWig | pyBigWig 软件发布版本, [10.5281/zenodo.45238](https://doi.org/10.5281/zenodo.45238) | 可选 MIT 适配器；仅限受限的本地 bigWig 统计。 |
| pybedtools | Dale et al. (2011), [10.1093/bioinformatics/btr539](https://doi.org/10.1093/bioinformatics/btr539); BEDTools： [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) | 可选 MIT 封装；仍需单独安装 BEDTools 可执行程序。 |
| gffutils | gffutils 文档与 daler 贡献者 | 可选 MIT 适配器；只读查询已有的 FeatureDB；调研未提供 DOI。 |
| pySBOL3 | Mitchell, Beal & Bartley (2022), [10.1021/acssynbio.2c00164](https://doi.org/10.1021/acssynbio.2c00164) | 可选 MIT 适配器；受限的只读 SBOL3 清单。 |
| scadnano | Doty, Lee & Stérin (2020), [arXiv:2005.11841](https://arxiv.org/abs/2005.11841) | 可选 MIT 适配器；受限的只读 `.sc` 设计清单。 |
| polars-bio | polars-bio, [10.1093/bioinformatics/btaf640](https://doi.org/10.1093/bioinformatics/btaf640) | 可选 Apache-2.0 适配器；受限的本地区间交集。 |
| DnaCauldron | DnaCauldron 文档与 Edinburgh Genome Foundry 贡献者；调研未提供 DOI | 可选 MIT 适配器；仅限受限的 Gibson/Type IIS 模拟。 |
| DNA Features Viewer | [官方文档](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/) 与 [Edinburgh Genome Foundry 仓库](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)；调研未提供 DOI | 可选 MIT 适配器；受限的只读 GenBank/GFF 特征清单。 |
| gtars | [官方文档](https://docs.bedbase.org/gtars/python-overview/) 与 [databio 仓库](https://github.com/databio/gtars)；调研未提供 DOI | 可选 BSD-2-Clause 适配器；受限的只读本地 BED RegionSet 清单。 |
| GeneBlocks | [官方文档](https://edinburgh-genome-foundry.github.io/Geneblocks/) 与 [Edinburgh Genome Foundry 仓库](https://github.com/Edinburgh-Genome-Foundry/Geneblocks)；调研未提供 DOI | 可选 MIT 适配器；确定性的受限序列差异比较。 |
| sourmash | [JOSS 论文](https://doi.org/10.21105/joss.06830) 与 [官方文档](https://sourmash.readthedocs.io/); BSD-3-Clause | 可选 BSD-3-Clause 适配器；受限的本地 MinHash 比较。 |
| Cutadapt | Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200); [官方指南](https://cutadapt.readthedocs.io/en/stable/guide.html); MIT | 受限的外部可执行程序适配器，用于 FASTQ 接头修剪。 |
| fastp | Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); [官方仓库](https://github.com/OpenGene/fastp); MIT | 受限的外部可执行程序适配器，用于 FASTQ 预处理和质控。 |
| FastQC | [官方项目页](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); GPL-3.0；调研未提供 DOI | 受限的外部可执行程序适配器，生成 HTML/ZIP 读段质控报告。 |
| MultiQC | Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354); [官方仓库](https://github.com/MultiQC/MultiQC); GPL-3.0-or-later | 受限的外部可执行程序适配器，生成汇总质控报告。 |

截图中的十个项目当前逐项路由如下：

| 项目 | DNAKit 路由 | 当前边界 |
| --- | --- | --- |
| Biopython | `dnakit.bioinformatics` | 已有格式、读取、比对和 motif 相关入口；不是完整命名空间重导出 |
| SeqKit | `dnakit.fasta_fastq` | 已有 FASTA/FASTQ 文件操作；其余 CLI 仍需单独安装 |
| scikit-bio | `dnakit.integrations.beta_diversity_scikit_bio`、`permanova_scikit_bio` | 可选 BSD-3-Clause 后端；只暴露受限统计接口 |
| UGENE | `dnakit.external_adapters.run` | GPLv2 外部程序，不复制进核心包 |
| iLearnPlus | `dnakit.external_adapters.run` | 许可和再分发范围仍需作者确认 |
| iFeatureOmega | `dnakit.external_adapters.run` | 许可和再分发范围仍需作者确认 |
| `ga` | `dnakit.external_adapters.run` | GPLv3 通用遗传算法包，不作为 DNAKit 原生 DNA API |
| seqme | `dnakit.integrations.evaluate_seqme` | 可选 BSD-3-Clause 后端；模型权重另行审查 |
| pydna | `dnakit.bioinformatics.cloning` | 已有构建相关入口；不是完整 pydna API 重导出 |
| DNA Chisel | `dnakit.integrations.optimize_with_dnachisel` | 可选 MIT 后端；仅允许清单内的 specification 类 |
| Biotite | `dnakit.integrations.align_biotite` | 可选 BSD-3-Clause 后端；有界最优核苷酸成对比对 |
| Edlib | `dnakit.integrations.align_edlib` | 可选 MIT 后端；有界编辑距离、位置或 CIGAR 路径 |
| PyRanges | `dnakit.integrations.summarize_pyranges_bed`、`merge_pyranges_bed` | 可选 MIT 后端；有界 BED 区间汇总和合并，超限拒绝 |
| pyfastx | `dnakit.integrations.summarize_pyfastx` | 可选 MIT 后端；有界 FASTA/FASTQ 汇总，超限拒绝 |
| dnaio | `dnakit.integrations.summarize_dnaio_fastx` | 可选 MIT 后端；只读 FASTA/FASTQ 流汇总，超限拒绝 |
| mappy/minimap2 | `dnakit.integrations.map_mappy_sequences` | 可选 MIT 后端；有界本地序列映射，预设值白名单 |
| pyBigWig | `dnakit.integrations.summarize_pybigwig_region` | 可选 MIT 后端；有界本地 bigWig 区域统计 |
| pybedtools | `dnakit.integrations.intersect_pybedtools_bed` | 可选 MIT Python 包；有界 BED 交集，仍需外部 bedtools |
| gffutils | `dnakit.integrations.query_gffutils_database` | 可选 MIT 后端；只读查询已有 FeatureDB，不创建或修改数据库 |
| pySBOL3 | `dnakit.integrations.summarize_sbol3_document` | 可选 MIT 后端；只读汇总本地 SBOL3 文档 |
| scadnano | `dnakit.integrations.summarize_scadnano_design` | 可选 MIT 后端；只读汇总本地 `.sc` 设计 |
| polars-bio | `dnakit.integrations.overlap_polars_bio_bed` | 可选 Apache-2.0 后端；有界本地区间交集 |
| DnaCauldron | `dnakit.integrations.simulate_dnacauldron_assembly` | 可选 MIT 后端；有界 Gibson/Type IIS 组装模拟 |
| DnaFeaturesViewer | `dnakit.integrations.summarize_dna_features_viewer` | 可选 MIT 后端；有界只读 GenBank/GFF 特征清单 |
| gtars | `dnakit.integrations.summarize_gtars_regions` | 可选 BSD-2-Clause 后端；有界只读 BED 区间清单 |
| GeneBlocks | `dnakit.integrations.diff_geneblocks_sequences` | 可选 MIT 后端；有界确定性 DNA 序列差异比较 |
| sourmash | `dnakit.integrations.compare_sourmash_sequences` | 可选 BSD-3-Clause 后端；有界本地 MinHash 比较 |
| Cutadapt | `dnakit.sequencing_qc.run_cutadapt` | MIT 外部可执行程序适配器；固定参数白名单的 FASTQ 接头清理 |
| fastp | `dnakit.sequencing_qc.run_fastp` | MIT 外部可执行程序适配器；固定参数白名单的 FASTQ 预处理和 QC |
| FastQC | `dnakit.sequencing_qc.run_fastqc` | GPL-3.0 外部可执行程序适配器；监控 HTML/ZIP QC 报告 |
| MultiQC | `dnakit.sequencing_qc.run_multiqc` | GPL-3.0-or-later 外部可执行程序适配器；监控汇总 HTML 报告 |
| pyfaidx | `dnakit.integrations.fetch_pyfaidx_region` | 可选 BSD-3-Clause 后端；只读 FASTA 区间 |
| vcfpy | `dnakit.integrations.summarize_vcf_vcfpy` | 可选 MIT 后端；有界 VCF 汇总 |
| cyvcf2 | `dnakit.integrations.summarize_vcf_cyvcf2` | 可选 MIT 后端；有界 VCF 汇总 |

记录的 `citation_status` 为 `missing_in_survey` 时，表示来源文件没有提供 DOI
或推荐引用格式。这是发表前需要补齐的引用缺口：在论文声称使用该工具之前，
应先从上游的 `CITATION.cff`、论文或文档中确认。宽松的代码许可证并不授予再分发
托管 API、模型权重、参考数据库或用户数据的权利。受限、非商业和未经验证的项目
只保留在审计矩阵中，不会复制进 DNAKit。

编辑带日期的来源文件后，用
`new_vision/2026-10-10/build_public_tool_coverage.py` 重新生成证据。该命令离线运行，
不会核实上游当前的许可证，也不会执行科学计算后端。实现阶段的文件清单记录在
`new_vision/2026-10-10/implementation_files_manifest.txt`。
