# public DNA tool coverage

`new_vision/2026-10-10/public_dna_tools.md` is a dated source survey.  DNAKit
now keeps its project-by-project cross-reference in
`new_vision/2026-10-10/public_dna_tools_coverage.{json,csv,md}` and ships the
same records through `dnakit.integrations`.

The source text declares 282 spreadsheet records.  The Markdown that is
available in this repository contains 425 explicit mentions, which collapse
to 320 normalized source records. The survey itself declares 282 records;
repeated classification, citation and provider-specific lists are retained as
separate evidence, while spelling aliases are merged conservatively.  The report also
records the source's declared 169 conditionally integrable records separately.

每个纳入范围的候选能力都使用正式的 `dnakit.<module>.<function>()` 公开入口。
`dnakit_candidates` 只是实现登记使用的内部暂存命名空间，不是公开 API。正式入口
可能是本地实现、有限门面或带依赖检查的适配器，具体状态在下文单独记录。

## Querying the coverage matrix

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

The route is intentionally explicit:

- `analysis_adapter` maps to `dnakit.integrations.run_analysis` and still
  requires the separately installed executable and reference data.
- `mcp_registry` maps to the MCP bridge and still requires a configured local
  server or author-documented endpoint.
- `skill_registry` maps to the bundled Skill/workflow layer; the Skill's own
  license, model and data conditions remain in force.
- `pip_adapter` maps to an explicit optional Python backend such as Biopython,
  pydna, seqme, DNA Chisel, Biotite, Edlib, PyRanges, pyfastx, dnaio, mappy, pyBigWig, pybedtools,
  gffutils, pySBOL3, scadnano, polars-bio, DnaCauldron, DNA Features Viewer, gtars, GeneBlocks or sourmash; the package and any model extras must be
  installed separately.
- `native_facade` maps a surveyed tool to a focused DNAKit implementation,
  such as the BWA/Bowtie2/SAMtools/BCFtools/SPAdes/Picard NGS facades, the NCBI
  BLAST+/minimap2/BEDTools alignment facades, the MACS3 peak facade, the
  SeqKit FASTA/FASTQ functions, the FastTree tree facade, the SPOA consensus,
  vcfanno/SURVIVOR
  variant facades and the NanoFilt/Filtlong/NanoPlot/Porechop/Trimmomatic/AfterQC/SeqMagick/GATK4/Shasta/Merqury/BUSCO/DeepVariant/cuteSV/Sniffles2/DELLY/HTSlib/seqtk/WhatsHap/BamTools/VSEARCH/BEDOPS/Flye/hifiasm/MEGAHIT/DIAMOND/MMseqs2/Kraken2/Prokka/Bakta/CheckM2/Bracken/Centrifuge/Kaiju/Unicycler/GTDB-Tk/HUMAnN/MethylDackel/GEMINI/Ensembl VEP
  read-processing facades; it does not mean that the upstream command's complete CLI has
  been reimplemented.
- `project_evidence_only` and `no_local_match` are discovery records.  They
  route to `dnakit.external_adapters.run` only as an explicit external-adapter
  boundary; they are not claims that the upstream source has been copied or
  that all of its functions are implemented.

`local_match` means a name, repository, dependency, Skill notice, MCP entry or
analysis entry was found in DNAKit.  It is not a scientific validation result.

Call `integration_inventory()` for current surface counts. It reports Python
functions, Agent-compatible Python functions, MCP server recipes, registered
Skills, surveyed projects and candidate schemas as separate units; these
categories must not be added together as if they were interchangeable tools.
The current raw surface sum is **2,432** (1,826 public Python functions + 32
MCP server recipes + 574 registered Skills); it is a surface count rather than
a deduplicated biological-tool count because one capability may appear in
several surfaces.
The inventory also includes `unified_candidate_statuses`, separating
implemented, partial and adapter candidates.
It also keeps the survey's `records_by_kind` and integration-status counts so
the 103 surveyed MCP records are not mistaken for 103 locally connected servers.
At this snapshot the 320 survey records contain 62 `native_facade`, 2
`analysis_adapter`, 26 `pip_adapter`, 13 `mcp_registry`, 2 `skill_registry`,
44 `project_evidence_only` and 171 `no_local_match` records.

The shipped `unified_capabilities.json` also records 837 formal candidate API
names and their input/output JSON Schemas. `list_unified_capabilities()` and
`describe_unified_capability()` make this catalog available through the same
Python and MCP discovery layer. Thirty-two project records currently have a
conservative name/schema match to one or more capabilities; these are review
links, not evidence that the upstream project's complete API is implemented.

`run_unified_capability()` is the execution bridge for a catalog capability. It
resolves the formal `dnakit.<module>.<operation>` entry through the existing
Agent/MCP registry, applies its JSON argument schema and serialization limits,
and requires explicit `allow_side_effects=True` for file writes, downloads,
external programs or model operations. It does not install an upstream package
or silently execute an unverified project.

The survey's permissively licensed pip tools have explicit bounded adapters:

- Existing `dnakit.bioinformatics` format/alignment functions use Biopython
  when `dnakit[bio-formats]` is installed, reads functions use pysam when
  `dnakit[bio-reads]` is installed, and cloning functions use pydna when
  `dnakit[bio-assembly]` is installed. These are focused DNAKit facades, not
  re-exports of the complete upstream namespaces.

- `dnakit.integrations.beta_diversity_scikit_bio()` and
  `dnakit.integrations.permanova_scikit_bio()` expose bounded scikit-bio
  statistics. Install `dnakit[bio-skibio]`; the full scikit-bio namespace is
  not re-exported.
- `dnakit.integrations.evaluate_seqme()` constructs allow-listed
  `seqme.metrics` classes and serializes the evaluation table. Install
  `dnakit[bio-seqme]`; model extras and weights remain separate.
- `dnakit.integrations.optimize_with_dnachisel()` constructs an allow-listed
  set of DNA Chisel specifications and returns the optimized sequence. Install
  `dnakit[bio-dnachisel]`; it does not copy DNA Chisel source code.
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

## Citation and license evidence

The survey provides paper identifiers for the following central tools; use the
DOI in the paper and retain the upstream software notice in any distribution:

| Tool | Citation | DNAKit boundary |
| --- | --- | --- |
| Biopython | Cock et al. (2009), [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) | `dnakit.bioinformatics` adapters; Biopython's own license and third-party notices apply. |
| SeqKit | Shen, Sipos & Zhao (2024), SeqKit2, [10.1002/imt2.191](https://doi.org/10.1002/imt2.191); the earlier SeqKit paper is [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) | Use a separately installed CLI through an external-analysis adapter. |
| scikit-bio | Aton et al. (2026), [10.1038/s41592-025-02981-z](https://doi.org/10.1038/s41592-025-02981-z) | Optional Python dependency; the BSD-3-Clause notice remains with the dependency. |
| UGENE | Okonechnikov et al. (2012), [10.1093/bioinformatics/bts091](https://doi.org/10.1093/bioinformatics/bts091) | External GPLv2 program; it is not copied into the core package. |
| pydna | Pereira et al. (2015), [10.1186/s12859-015-0544-x](https://doi.org/10.1186/s12859-015-0544-x) | Optional DNA-assembly dependency; retain the BSD notice. |
| DNA Chisel | Zulkower & Rosser (2020), [10.1093/bioinformatics/btaa558](https://doi.org/10.1093/bioinformatics/btaa558) | Optional MIT dependency; model and data terms are separate. |
| seqme | Møller-Larsen et al. (2025), [arXiv:2511.04239](https://arxiv.org/abs/2511.04239) | Optional BSD-3-Clause adapter; model extras and weights are separate. |
| Biotite | Kunzmann & Hamacher (2018), [10.1186/s12859-018-2367-z](https://doi.org/10.1186/s12859-018-2367-z) | Optional BSD-3-Clause adapter; bounded nucleotide alignment only. |
| Edlib | Šošić & Šikić (2017), [10.1093/bioinformatics/btw753](https://doi.org/10.1093/bioinformatics/btw753) | Optional MIT adapter; bounded edit-distance alignment only. |
| PyRanges | Stovner & Sætrom (2020), [10.1093/bioinformatics/btz615](https://doi.org/10.1093/bioinformatics/btz615) | Optional MIT adapter; bounded BED interval operations only. |
| pyfastx | Du et al. (2021), [10.1093/bib/bbaa368](https://doi.org/10.1093/bib/bbaa368) | Optional MIT adapter; bounded FASTA/FASTQ summary only. |
| pyfaidx | Shirley et al. (2015), [10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1); source code BSD-3-Clause | Optional BSD-3-Clause adapter; bounded read-only FASTA interval fetching. |
| dnaio | dnaio software archive, [10.5281/zenodo.10548864](https://doi.org/10.5281/zenodo.10548864); MIT | Optional MIT adapter; bounded read-only FASTA/FASTQ stream summary. |
| mappy/minimap2 | Li (2018), [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191); [official Python binding README](https://github.com/lh3/minimap2/blob/master/python/README.rst); MIT | Optional MIT adapter; bounded local mapping with allow-listed presets. |
| pyBigWig | pyBigWig software release, [10.5281/zenodo.45238](https://doi.org/10.5281/zenodo.45238) | Optional MIT adapter; bounded local bigWig statistics only. |
| pybedtools | Dale et al. (2011), [10.1093/bioinformatics/btr539](https://doi.org/10.1093/bioinformatics/btr539); BEDTools: [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) | Optional MIT wrapper; a separately installed BEDTools executable remains required. |
| gffutils | gffutils documentation and daler contributors | Optional MIT adapter; read-only queries against an existing FeatureDB; no DOI was supplied by the survey. |
| pySBOL3 | Mitchell, Beal & Bartley (2022), [10.1021/acssynbio.2c00164](https://doi.org/10.1021/acssynbio.2c00164) | Optional MIT adapter; bounded read-only SBOL3 inventory. |
| scadnano | Doty, Lee & Stérin (2020), [arXiv:2005.11841](https://arxiv.org/abs/2005.11841) | Optional MIT adapter; bounded read-only `.sc` design inventory. |
| polars-bio | polars-bio, [10.1093/bioinformatics/btaf640](https://doi.org/10.1093/bioinformatics/btaf640) | Optional Apache-2.0 adapter; bounded local interval overlap. |
| DnaCauldron | DnaCauldron documentation and Edinburgh Genome Foundry contributors; the survey supplied no DOI | Optional MIT adapter; bounded Gibson/Type IIS simulation only. |
| DNA Features Viewer | [official documentation](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/) and [Edinburgh Genome Foundry repository](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer); the survey supplied no DOI | Optional MIT adapter; bounded read-only GenBank/GFF feature inventory. |
| gtars | [official documentation](https://docs.bedbase.org/gtars/python-overview/) and [databio repository](https://github.com/databio/gtars); the survey supplied no DOI | Optional BSD-2-Clause adapter; bounded read-only local BED RegionSet inventory. |
| GeneBlocks | [official documentation](https://edinburgh-genome-foundry.github.io/Geneblocks/) and [Edinburgh Genome Foundry repository](https://github.com/Edinburgh-Genome-Foundry/Geneblocks); the survey supplied no DOI | Optional MIT adapter; deterministic bounded sequence diff. |
| sourmash | [JOSS paper](https://doi.org/10.21105/joss.06830) and [official documentation](https://sourmash.readthedocs.io/); BSD-3-Clause | Optional BSD-3-Clause adapter; local bounded MinHash comparison. |
| Cutadapt | Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200); [official guide](https://cutadapt.readthedocs.io/en/stable/guide.html); MIT | Bounded external executable adapter for FASTQ adapter trimming. |
| fastp | Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); [official repository](https://github.com/OpenGene/fastp); MIT | Bounded external executable adapter for FASTQ preprocessing and QC. |
| FastQC | [official project page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); GPL-3.0; the survey supplied no DOI | Bounded external executable adapter for HTML/ZIP read-QC reports. |
| MultiQC | Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354); [official repository](https://github.com/MultiQC/MultiQC); GPL-3.0-or-later | Bounded external executable adapter for aggregated QC reports. |

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

The record's `citation_status` is `missing_in_survey` when the source file did
not provide a DOI or preferred citation.  That status is a publication gap to
resolve from the upstream `CITATION.cff`, paper or documentation before a
paper claims to have used that tool.  A permissive code license does not grant
rights to redistribute a hosted API, model weights, reference database or
user data.  Restricted, non-commercial and unverified projects stay in the
audit matrix and are not copied into DNAKit.

Regenerate the evidence after editing the dated source file with
`new_vision/2026-10-10/build_public_tool_coverage.py`.  The command is offline;
it does not verify current upstream licenses or execute scientific backends.
The implementation-pass file list is recorded in
`new_vision/2026-10-10/implementation_files_manifest.txt`.
