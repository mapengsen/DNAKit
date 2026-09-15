"""Build the complete public API reference for the 484-function skill migration.

Run with the source checkout installed: python others/docs/build_native_api_docs.py
Use --check to verify committed pages, signatures, schemas and coverage offline.
No provider is initialized and no scientific operation is executed by this script.
"""

# Bilingual Markdown prose intentionally preserves punctuation and paragraph lines.
# ruff: noqa: E501, RUF001

from __future__ import annotations

import argparse
import html
import inspect
import json
import re
from collections import Counter
from pathlib import Path

from dnakit import databases, encode, workflows
from dnakit.databases._runtime import definitions
from dnakit.download import ena, ncbi
from dnakit.search.ncbi import accessions

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "docs/api/native"

# Explicit source-module membership avoids silently classifying new providers.
GROUPS = {
    "sequences": (
        "基因、序列、测序数据与下载",
        "Genes, sequences, sequencing data and downloads",
        "ena_portal_tool ensembl_compara_tool ensembl_sequence_tool ensembl_tool "
        "hgnc_tool ncbi_datasets_tool ncbi_nucleotide_tool ncbi_sra_tool sra_tool "
        "tark_tool file_download_tool",
    ),
    "variants": (
        "遗传变异、GWAS 与变异评分",
        "Genetic variants, GWAS and variant scores",
        "alphamissense_tool cadd_tool civic_tool clingen_dosage_tool clingen_tool "
        "clinvar_tool cosmic_tool dbsnp_tool ensembl_variation_ext_tool ensembl_vep_tool "
        "eve_tool favor_tool finngen_tool genomics_gene_search_tool gnomad_tool "
        "gwas_tool igsr_tool ldlink_tool mavedb_tool oncokb_tool pheweb_phewas_tool "
        "progenetix_tool spliceai_tool",
    ),
    "phenotypes": (
        "疾病、表型与本体",
        "Diseases, phenotypes and ontologies",
        "disgenet_tool efo_tool gencc_tool hpo_tool marrvel_tool monarch_v3_tool "
        "ols_tool omim_tool orphanet_tool restful_tool",
    ),
    "regulation": (
        "表达、表观调控与肿瘤数据",
        "Expression, regulation and cancer data",
        "cbioportal_tool cellxgene_census_tool chipatlas_tool depmap_tool encode_tool "
        "epigenomics_tool expression_atlas_tool fourdn_tool gdc_tool gtex_tool "
        "gtex_v2_tool harmonizome_tool hpa_tool jaspar_tool regulomedb_tool remap_tool "
        "screen_tool unibind_tool",
    ),
    "organisms": (
        "模式生物、物种与微生物组",
        "Model organisms, taxonomy and microbiomes",
        "alliance_genome_tool gbif_tool gtdb_tool mgnify_expanded_tool mgnify_tool "
        "rgd_strain_tool rgd_tool sgd_tool usda_plants_tool wormbase_tool xenbase_tool",
    ),
    "pathways": (
        "通路、互作、蛋白与结构",
        "Pathways, interactions, proteins and structures",
        "alphafold_tool biogrid_tool blast_tool bvbrc_tool enrichr_ext_tool enrichr_tool "
        "gene_ontology_tool humanbase_tool iedb_tool imgt_tool immport_tool intact_tool "
        "interpro_entry_tool interpro_tool kegg_ext_tool kegg_tool panther_tool "
        "pdbe_sifts_tool pfam_tool plant_reactome_tool proteins_api_tool quickgo_tool "
        "rcsb_data_tool reactome_analysis_tool reactome_content_tool reactome_tool "
        "string_ext_tool string_tool structure_annotation_tool uniprot_tool wikipathways_tool",
    ),
    "targets": (
        "靶点、药物与关联证据",
        "Targets, drugs and association evidence",
        "chem_tool ctd_tool dataset_tool dgidb_tool graphql_tool metabolite_tool "
        "openfda_tool pubchem_tool xml_tool",
    ),
    "literature": (
        "文献、预印本与临床研究检索",
        "Literature, preprints and study discovery",
        "arxiv_tool biorxiv_tool crossref_tool ctg_tool europe_pmc_tool openalex_tool "
        "pubmed_tool semantic_scholar_tool web_search_tool",
    ),
    "models": (
        "已训练模型的预测接口",
        "Prediction with trained models",
        "alphagenome_tool dynamut2_tool esm2_variant_effect_tool esm_tool "
        "evo2_variant_effect_tool nvidia_nim_tool",
    ),
    "analysis": (
        "序列计算、VCF、克隆与遗传模拟",
        "Sequence analysis, VCF, cloning and inheritance",
        "dn_ds_tool dna_tools scientific_calculator_tools sequence_analyze_tool "
        "variant_fraction_tool vcf_stats_tool",
    ),
    "encode": ("ENCODE 查询、下载与数据管理", "ENCODE queries, downloads and data management", ""),
    "execution": (
        "接口发现、工作流与报告",
        "API discovery, workflows and reports",
        "executed_notebook_tool tool_discovery_tools",
    ),
}

OVERRIDES = {
    "bio_rxiv_list_recent_preprints": "literature",
    "fly_mine_search": "organisms",
    "gbif_get_species": "organisms",
    "powo_search_plants": "organisms",
    "ncbigene_get_summary": "sequences",
    "ncbigene_search": "sequences",
    "ncbiprotein_get_summary": "pathways",
    "pdbe_get_uniprot_mappings": "pathways",
    "string_get_interaction_partners": "pathways",
    "string_get_network": "pathways",
    "variant_validator_gene2transcripts": "variants",
    "variant_validator_validate_variant": "variants",
    "e_qtl_get_associations": "regulation",
    "e_qtl_list_datasets": "regulation",
    "my_gene_batch_query": "sequences",
    "my_gene_get_gene_annotation": "sequences",
    "my_gene_query_genes": "sequences",
    "my_variant_get_pathogenicity_scores": "variants",
    "my_variant_get_variant_annotation": "variants",
    "my_variant_query_variants": "variants",
    "pdb_get_structure": "pathways",
    "orphanet_get_disease_genes": "phenotypes",
    "cadd_get_scores": "variants",
    "string_get_functional_enrichment": "pathways",
}

WORKFLOWS = (
    "scripts",
    "resource_path",
    "run",
    "run_support",
    "simulate_inheritance",
    "write_wes_report",
)
HELPERS = ("call", "describe", "execute", "list_operations", "search")
PILOT = {
    "dnakit.search.ncbi.accessions": accessions,
    "dnakit.download.ncbi.sequence": ncbi.sequence,
    "dnakit.download.ena.sequence": ena.sequence,
}

# Explanations that cannot be inferred from type annotations alone.
PARAMETERS = {
    "verbose": (
        "是否显示 DNAKit 运行提示；不会改变返回值。",
        "Show DNAKit progress messages without changing results.",
    ),
    "operation": (
        "本地操作名或登记的上游别名。",
        "Native operation name or registered upstream alias.",
    ),
    "arguments": (
        "数据库调用传参数字典；工作流传 CLI 字符串列表。",
        "Parameter mapping for database calls; CLI string list for workflows.",
    ),
    "request": (
        "操作名，或包含 name、arguments 的字典。",
        "Operation name or a mapping with name and arguments.",
    ),
    "query": (
        "本地目录搜索词，匹配名称与描述；使用英文功能词或函数名。",
        "Local name/description search; use English task terms or a function name.",
    ),
    "limit": (
        "最多返回的匹配条数；数据库目录搜索允许 1–100。",
        "Maximum matches; database catalog search allows 1–100.",
    ),
    "max_result_bytes": (
        "序列化结果的字节上限，超限抛出异常。",
        "Serialized result byte limit; exceeding it raises an exception.",
    ),
    "skill_id": (
        "已收录技能名称；兼容登记过的旧名称。",
        "Bundled skill ID or registered legacy alias.",
    ),
    "resource": (
        "技能包内相对路径，读取前校验哈希。",
        "Relative resource path within the skill; verified before use.",
    ),
    "script": (
        "技能包内 .py、.R 或 .sh 脚本相对路径。",
        "Relative .py, .R or .sh script path inside the bundle.",
    ),
    "workdir": (
        "子进程工作目录；None 使用当前目录。",
        "Subprocess working directory; None uses the current directory.",
    ),
    "timeout_seconds": ("子进程执行超时秒数。", "Subprocess timeout in seconds."),
    "max_output_bytes": (
        "最多捕获的子进程输出字节数。",
        "Maximum captured subprocess output bytes.",
    ),
    "name": (
        "支持流程名称，允许值见本条说明。",
        "Support workflow name; allowed values are described below.",
    ),
    "father": (
        "父方记录，sex 必须为 Male；基因型字段遵循 Genomebook 数据结构。",
        "Father record with sex=Male and Genomebook genotype fields.",
    ),
    "mother": (
        "母方记录，sex 必须为 Female；基因型字段遵循 Genomebook 数据结构。",
        "Mother record with sex=Female and Genomebook genotype fields.",
    ),
    "trait_registry": ("显式提供的性状与位点定义。", "Explicit trait and locus registry."),
    "disease_registry": ("显式提供的疾病位点定义。", "Explicit disease locus registry."),
    "offspring_count": ("模拟后代数量，1–10000。", "Number of simulated offspring, 1–10000."),
    "generation": ("代次标签，正整数。", "Positive integer generation label."),
    "mutation_rate": (
        "每个位点的模拟突变率，范围 0–1。",
        "Per-locus simulated mutation rate, between 0 and 1.",
    ),
    "seed": (
        "局部随机数种子；不修改全局随机状态。",
        "Local random seed; does not modify global RNG state.",
    ),
    "input_report": (
        "已有 UTF-8 Markdown 报告路径，最多 5 MB。",
        "Existing UTF-8 Markdown report path, at most 5 MB.",
    ),
    "output_dir": (
        "报告输出目录，生成 PDF 和来源清单。",
        "Report output directory for the PDF and provenance manifest.",
    ),
    "language": (
        "报告标签语言 en 或 es；不翻译输入正文。",
        "Report label language en or es; input text is not translated.",
    ),
    "title": (
        "自定义报告标题；None 使用语言对应默认标题。",
        "Custom report title; None uses the language default.",
    ),
    "logo_left": ("左侧 Logo 图片路径。", "Left logo image path."),
    "logo_right": ("右侧 Logo 图片路径。", "Right logo image path."),
    "font_path": ("可选 TTF 字体路径。", "Optional TTF font path."),
    "overwrite": ("是否覆盖已存在的报告文件。", "Allow overwriting existing report files."),
    "uids": (
        "NCBI 序列 UID 列表，按每批 200 个解析。",
        "NCBI sequence UIDs, resolved in batches of 200.",
    ),
    "database": (
        "NCBI 序列数据库 nuccore 或 protein。",
        "NCBI sequence database: nuccore or protein.",
    ),
    "config": (
        "查询使用 SearchConfig；下载使用 DownloadConfig。",
        "SearchConfig for queries; DownloadConfig for downloads.",
    ),
    "api_base_url": (
        "服务 API 根地址；默认值见签名。",
        "Service API base URL; see the signature for its default.",
    ),
    "accession": (
        "版本化序列 accession；ENA 使用 INSDC 编号，NCBI 也支持 RefSeq。",
        "Versioned sequence accession; ENA uses INSDC, NCBI also accepts RefSeq.",
    ),
    "output_path": ("目标序列文件路径。", "Destination sequence file path."),
    "format": (
        "fasta、genbank 或 embl；提供方及转换依赖见序列下载教程。",
        "fasta, genbank or embl; see the download guide for provider and conversion support.",
    ),
    "seq_start": (
        "NCBI 核酸子区间起点，1-based 闭区间；须同时指定 seq_stop。",
        "NCBI nucleotide start, 1-based inclusive; requires seq_stop.",
    ),
    "seq_stop": (
        "NCBI 核酸子区间终点，1-based 闭区间；须同时指定 seq_start。",
        "NCBI nucleotide end, 1-based inclusive; requires seq_start.",
    ),
    "strand": (
        "NCBI 核酸序列链方向，1 正向、2 反向互补。",
        "NCBI nucleotide strand: 1 forward, 2 reverse complement.",
    ),
    "search_config": (
        "控制查询超时、返回大小和 NCBI 凭据的 SearchConfig。",
        "SearchConfig for query timeout, response size and NCBI credentials.",
    ),
    "query_progress": ("查询步骤的显式进度回调。", "Explicit query progress callback."),
}

ZH_PURPOSE = {
    "dnakit.databases.call": "根据操作名和参数字典调用数据库接口。",
    "dnakit.databases.describe": "查看单个接口的实际输入 Schema、上游别名和实现信息。",
    "dnakit.databases.execute": "使用操作名或结构化请求分发一次调用。",
    "dnakit.databases.list_operations": "离线列出全部数据库操作名与上游别名。",
    "dnakit.databases.search": "按名称、任务关键词或描述离线搜索数据库接口。",
    "dnakit.databases.execute_tool": "兼容上游工具名的本地分发入口，实际调用 DNAKit 数据库函数。",
    "dnakit.workflows.scripts": "列出技能包中可运行的 Python、R 和 Shell 脚本。",
    "dnakit.workflows.resource_path": "取得已校验的随包资源路径。",
    "dnakit.workflows.run": "传递原脚本参数，执行随包工作流并返回真实退出码和输出。",
    "dnakit.workflows.run_support": "执行组合流程依赖的随包支持脚本。",
    "dnakit.workflows.simulate_inheritance": "根据显式提供的亲本与位点定义模拟等位基因遗传和突变。",
    "dnakit.workflows.write_wes_report": "把已有 Markdown 报告渲染为 PDF，并记录来源哈希。",
    "dnakit.search.ncbi.accessions": "把 NCBI 序列 UID 转为版本化 accession，保留输入顺序和重复项。",
    "dnakit.download.ncbi.sequence": "按 accession 下载 NCBI 序列，可选择格式、核酸子区间和链方向。",
    "dnakit.download.ena.sequence": "按 INSDC accession 下载 ENA 序列文件。",
}

ENCODE_PURPOSES = (
    "按实验类型、物种、组织等条件检索实验。",
    "获取指定实验的完整信息。",
    "列出实验关联的文件。",
    "按文件类型、组装和实验条件检索文件。",
    "下载指定文件并记录元数据。",
    "获取实验的结构化元数据。",
    "批量下载匹配实验的文件。",
    "配置、查看状态或清除 ENCODE 凭据。",
    "查询可用于筛选的字段取值及数量。",
    "获取指定文件的详细元数据。",
    "把实验加入本地追踪目录。",
    "列出本地追踪的实验。",
    "获取追踪实验的论文引用，可导出 BibTeX 或 RIS。",
    "比较两个已追踪实验的兼容性。",
    "登记衍生文件及其来源、工具和参数。",
    "查询衍生文件的来源链。",
    "把本地追踪的实验导出为 CSV、TSV 或 JSON 文本。",
    "按实验、组织等字段汇总本地实验集合。",
    "关联实验与 DOI、PMID、GEO 等外部标识。",
    "查询实验关联的外部引用。",
)
ZH_PURPOSE.update(
    {f"dnakit.encode.{n}": text for n, text in zip(encode.__all__, ENCODE_PURPOSES, strict=True)}
)

# These routes deliberately differ from the retained provider descriptions.
PURPOSE_OVERRIDES = {
    "dnakit.databases.execute_tool": (
        "根据工具名与参数调用本地 DNAKit 数据库操作；禁止递归调用自身。",
        "Dispatch a tool name and arguments to a local DNAKit database operation; recursive self-calls are rejected.",
    ),
    "dnakit.databases.ch_ipatlas_enrichment_analysis": (
        "根据 BED、motif 或基因列表生成 ChIP-Atlas 富集任务的提交说明；本函数不计算富集统计量。",
        "Prepare ChIP-Atlas enrichment submission instructions from BED, motif or gene-list input; this function does not calculate enrichment statistics.",
    ),
    "dnakit.databases.ch_ipatlas_get_peak_data": (
        "查询 ChIP-Atlas 实验峰数据入口；gene 模式先定位 GRCh38 基因区域，提供 peak_file 后计算 BED 重叠，否则返回数据获取计划。",
        "Query ChIP-Atlas experiment peak resources. The gene mode resolves a GRCh38 locus and intersects a supplied peak_file BED, or returns a data retrieval plan when no BED is supplied.",
    ),
}

DATABASE_PARAMETERS = {
    "category": (
        "分类或结果类别，按当前接口的枚举/来源定义填写。",
        "Category or result class defined by this operation's provider/schema.",
    ),
    "query": (
        "检索词或查询表达式；语法由当前数据源定义。",
        "Search term or expression in the current provider's query syntax.",
    ),
    "limit": (
        "当前请求的返回数量上限；不是全库命中总数。",
        "Maximum results for the current request, not a database-wide total.",
    ),
    "gene": (
        "基因标识或基因列表；命名要求和类型见当前接口。",
        "Gene identifier or list; follow this operation's identifier conventions and type.",
    ),
    "tissue": (
        "目标组织名称或标识，按当前数据源填写。",
        "Tissue name or identifier understood by the provider.",
    ),
    "max_cells": (
        "本次返回或处理的细胞数上限。",
        "Maximum cells returned or processed in this request.",
    ),
    "offset": ("分页起始偏移；对应当前查询结果。", "Pagination offset within the current query."),
    "operation": (
        "提供方操作模式；默认值与允许值见本行约束。",
        "Provider operation mode; see this row's default and allowed values.",
    ),
    "experiment_type": (
        "实验类别，按数据源接受的类型筛选。",
        "Experiment type accepted by the source.",
    ),
    "peak_file": (
        "本地 BED 峰文件；必须核对其基因组版本和轨道来源。",
        "Local BED peak file with verified genome assembly and track provenance.",
    ),
    "flank_bp": (
        "基因区域两侧扩展的碱基数。",
        "Number of flanking bases added on each side of the gene locus.",
    ),
    "max_file_bytes": ("允许读取的峰文件字节上限。", "Maximum peak-file size in bytes."),
    "assay_title": ("实验方法名称，例如 ATAC-seq。", "Assay title, for example ATAC-seq."),
    "target": (
        "实验靶标或目标标识，按当前数据源填写。",
        "Assay target or target identifier recognized by the provider.",
    ),
    "organism": (
        "目标生物名称，按提供方命名规则填写。",
        "Organism name using the provider's naming convention.",
    ),
    "status": ("记录的发布或处理状态筛选条件。", "Record release or processing status filter."),
    "biosample": ("生物样本名称或检索条件。", "Biosample name or search filter."),
    "source": (
        "数据来源筛选值；按当前接口允许的来源填写。",
        "Source filter using values accepted by this operation.",
    ),
    "dataset_id": ("提供方的数据集标识。", "Provider dataset identifier."),
    "page": (
        "结果页码；起始页约定按提供方参数说明。",
        "Result page number; the starting page convention follows the provider.",
    ),
    "items_per_page": ("每页返回的记录数。", "Number of records requested per page."),
    "gene_symbol": ("基因符号，例如 BRCA1。", "Gene symbol, for example BRCA1."),
    "size": (
        "当前返回页的数量上限；不自动遍历后续页。",
        "Limit for the returned page; does not automatically retrieve later pages.",
    ),
    "uniprot_id": (
        "UniProt accession；反向结构映射时不要同时传 pdb_id。",
        "UniProt accession; do not also supply pdb_id for the reverse structure lookup.",
    ),
    "protein_ids": (
        "蛋白标识；输入列表或字符串类型见当前接口 Schema。",
        "Protein identifiers; use the list or string type specified by this operation.",
    ),
    "chrom": (
        "染色体标识；须与 start、end 一起传入。",
        "Chromosome identifier; provide with start and end.",
    ),
    "start": (
        "区域起点；本参数分支使用 Ensembl 的 1-based 闭区间。",
        "Region start; this parameter branch uses Ensembl 1-based inclusive coordinates.",
    ),
    "end": (
        "区域终点；本参数分支使用 Ensembl 的 1-based 闭区间。",
        "Region end; this parameter branch uses Ensembl 1-based inclusive coordinates.",
    ),
    "allele": (
        "VEP 区域变异的等位基因参数，遵循 Ensembl VEP 格式。",
        "Allele argument for a region variant in Ensembl VEP format.",
    ),
    "p_value_threshold": (
        "对当前返回页应用的 P 值过滤阈值；不代表已遍历全库。",
        "P-value filter applied to the current page, not an exhaustive database scan.",
    ),
    "format": (
        "返回格式；允许值与服务默认值见当前接口说明。",
        "Output format; see this operation's allowed values and provider defaults.",
    ),
    "species": (
        "物种名或 Taxonomy ID；类型与支持值取决于当前接口。",
        "Species name or taxonomy ID; type and supported values depend on the operation.",
    ),
    "pdb_id": (
        "PDB 结构编号，支持 4 字符编号或 pdb_ 加 8 字符编号。",
        "PDB structure identifier: four characters or pdb_ plus eight characters.",
    ),
    "max_bytes": (
        "坐标文件大小上限，1–100000000 字节。",
        "Coordinate file size limit, 1–100000000 bytes.",
    ),
    "orpha_code": ("Orphanet 的 ORPHA 疾病编号。", "Orphanet ORPHA disorder identifier."),
    "variant": (
        "用于 MyVariant CADD 注释检索的变异查询条件。",
        "Variant query used to retrieve indexed CADD annotations through MyVariant.",
    ),
    "assembly": ("参考基因组版本，例如 hg38。", "Reference genome assembly, for example hg38."),
    "identifiers": (
        "以换行或回车分隔的蛋白标识；将转换成 STRING 标识列表。",
        "Protein identifiers separated by newline or carriage return; converted to a STRING identifier list.",
    ),
}

RETURN_NOTES = {
    "dnakit.databases.pdb_get_structure": (
        "字典含 pdb_id、format、data（坐标文件文本）、url、sha256；不写入磁盘。",
        "Mapping with pdb_id, format, data (coordinate file text), url and sha256; does not write to disk.",
    ),
    "dnakit.databases.orphanet_get_disease_genes": (
        "返回 Orphanet_get_genes 的原始关联记录结构，见 [orphanet_get_genes](phenotypes.md#dnakit.databases.orphanet_get_genes)。",
        "Returns the Orphanet_get_genes association result; see [orphanet_get_genes](phenotypes.md#dnakit.databases.orphanet_get_genes).",
    ),
    "dnakit.databases.cadd_get_scores": (
        "返回 MyVariant 检索结果；字典结果另含 dnakit_score_source，表明分数来自已索引的 CADD 注释。缺失分数保持缺失，不执行新的模型预测。",
        "Returns MyVariant query results; mapping results include dnakit_score_source identifying indexed CADD annotations. Missing scores remain missing; no new model prediction is performed.",
    ),
    "dnakit.databases.string_get_functional_enrichment": (
        "返回 STRING_functional_enrichment 的富集结果，见 [string_functional_enrichment](pathways.md#dnakit.databases.string_functional_enrichment)。",
        "Returns STRING_functional_enrichment results; see [string_functional_enrichment](pathways.md#dnakit.databases.string_functional_enrichment).",
    ),
    "dnakit.databases.search": (
        "列表，每项包含 name、description、upstream_name。",
        "List of name, description and upstream_name mappings.",
    ),
    "dnakit.databases.list_operations": (
        "列表，每项包含 name、upstream_name。",
        "List of name and upstream_name mappings.",
    ),
    "dnakit.databases.describe": (
        "字典；input_schema 为实际输入约束，config 保存提供方配置，upstream_name 为原接口名。",
        "Mapping with actual input_schema, provider config and upstream_name.",
    ),
    "dnakit.workflows.scripts": (
        "技能包中的脚本相对路径列表；可能为空。",
        "List of relative script paths; may be empty.",
    ),
    "dnakit.workflows.resource_path": (
        "pathlib.Path；指向通过哈希校验的本地资源。",
        "pathlib.Path to a checksum-verified local resource.",
    ),
    "dnakit.workflows.run": (
        "字典包含 skill、script、returncode、output、upstream_url、upstream_revision、license；returncode 非零表示脚本失败。",
        "Mapping with skill, script, returncode, output, upstream_url, upstream_revision and license; nonzero returncode indicates failure.",
    ),
    "dnakit.workflows.run_support": (
        "字典包含 name、returncode、output；须检查退出码。",
        "Mapping with name, returncode and output; check the exit code.",
    ),
    "dnakit.workflows.simulate_inheritance": (
        "字典包含 offspring、seed、mutation_rate、scope、upstream；只模拟独立位点遗传，不模拟染色体交叉或真实临床结局。",
        "Mapping with offspring, seed, mutation_rate, scope and upstream; independent-locus simulation, not chromosome crossover or clinical prediction.",
    ),
    "dnakit.search.ncbi.accessions": (
        "tuple[str, ...]；空输入返回空元组，缺失 UID 或版本号抛出 QueryError。",
        "tuple[str, ...]; empty input gives an empty tuple, missing UIDs or versions raise QueryError.",
    ),
}


def tr(zh: str, en: str, language: str) -> str:
    return zh if language == "zh" else en


def cell(value: object) -> str:
    return (
        html.escape(str(value))
        .replace("|", "&#124;")
        .replace("[", "&#91;")
        .replace("]", "&#93;")
        .replace("\n", "<br>")
    )


def annotation(value: object) -> str:
    if value is inspect.Signature.empty:
        return "object"
    text = value if isinstance(value, str) else inspect.formatannotation(value)
    return text.replace("_Object", "object")


def signature(fn: object) -> str:
    """Keep Python defaults intact; avoid quoted future annotations in display."""
    sig = inspect.signature(fn)
    params = [p.replace(annotation=inspect.Signature.empty) for p in sig.parameters.values()]
    return str(sig.replace(parameters=params, return_annotation=inspect.Signature.empty))


def doc_sections(fn: object) -> tuple[str, dict[str, str], str]:
    doc = inspect.getdoc(fn) or ""
    before = re.split(r"\n(?:Args|Arguments|Parameters|Returns):", doc, maxsplit=1)[0]
    args: dict[str, str] = {}
    matched = re.search(r"\nArgs:\n(.*?)(?=\n\w[^\n]*:|\Z)", doc, re.S)
    if matched:
        key = None
        for line in matched[1].splitlines():
            item = re.match(r"    (\w+)(?:\s*\([^)]*\))?:\s*(.*)", line)
            if item:
                key = item[1]
                args[key] = item[2]
            elif key and line.strip():
                args[key] += " " + line.strip()
    returned = re.search(r"\nReturns:\n(.*?)(?=\n\w[^\n]*:|\Z)", doc, re.S)
    return before, args, inspect.cleandoc(returned[1]) if returned else ""


def inventory() -> list[dict]:
    entries = definitions()
    module_group = {
        module: group for group, values in GROUPS.items() for module in values[2].split()
    }
    records = []
    assert len(databases.__all__) == 455 and len(encode.__all__) == 20
    assert set(databases.__all__) == set(entries) | set(HELPERS)
    for name in databases.__all__:
        entry = entries.get(name)
        group = (
            "execution" if entry is None else OVERRIDES.get(name) or module_group[entry["module"]]
        )
        records.append(
            {
                "name": f"dnakit.databases.{name}",
                "fn": getattr(databases, name),
                "group": group,
                "entry": entry,
            }
        )
    for name in encode.__all__:
        records.append(
            {
                "name": f"dnakit.encode.{name}",
                "fn": getattr(encode, name),
                "group": "encode",
                "entry": None,
            }
        )
    for name in WORKFLOWS:
        records.append(
            {
                "name": f"dnakit.workflows.{name}",
                "fn": getattr(workflows, name),
                "group": "analysis" if name == "simulate_inheritance" else "execution",
                "entry": None,
            }
        )
    for name, fn in PILOT.items():
        records.append({"name": name, "fn": fn, "group": "sequences", "entry": None})
    assert len(records) == len({r["name"] for r in records}) == 484
    return sorted(records, key=lambda r: r["name"])


def purpose(record: dict, language: str) -> str:
    if record["name"] in PURPOSE_OVERRIDES:
        return PURPOSE_OVERRIDES[record["name"]][0 if language == "zh" else 1]
    if language == "zh" and record["name"] in ZH_PURPOSE:
        return ZH_PURPOSE[record["name"]]
    entry = record["entry"]
    return entry["config"]["description"] if entry else doc_sections(record["fn"])[0]


def parameter_table(record: dict, language: str) -> str:
    entry = record["entry"]
    schema = entry.get("input_schema", {}) if entry else {}
    properties = schema.get("properties", {})
    required = schema.get("required", [])
    _, arg_docs, _ = doc_sections(record["fn"])
    rows = [
        tr(
            "| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |",
            "| Parameter | Type | Required / default | Description and constraints |",
            language,
        ),
        "| --- | --- | --- | --- |",
    ]
    sig = inspect.signature(record["fn"])
    assert set(properties) <= set(sig.parameters), record["name"]
    for name, p in sig.parameters.items():
        prop = properties.get(name, {})
        kind = prop.get("type", annotation(p.annotation))
        if isinstance(kind, list):
            kind = " | ".join(kind)
        if name in required:
            default = tr("必填（输入 Schema）", "Required by input schema", language)
        elif p.default is inspect.Signature.empty:
            default = tr("必填", "Required", language)
        else:
            default = repr(p.default)
        desc = prop.get("description") or arg_docs.get(name)
        if not desc and entry and name in DATABASE_PARAMETERS:
            desc = DATABASE_PARAMETERS[name][0 if language == "zh" else 1]
        if not desc and not entry:
            pair = PARAMETERS.get(name)
            desc = pair[0 if language == "zh" else 1] if pair else None
        if not desc:
            desc = tr(
                "运行提示开关。" if name == "verbose" else "参数含义见用途与完整输入 Schema。",
                "Progress flag."
                if name == "verbose"
                else "See the purpose and complete input schema.",
                language,
            )
        limits = {k: v for k, v in prop.items() if k not in {"type", "description", "title"}}
        if limits:
            desc += " " + json.dumps(limits, ensure_ascii=False)
        rows.append(f"| `{name}` | {cell(kind)} | {cell(default)} | {cell(desc)} |")
    return "\n".join(rows)


def details(title: str, content: str, kind: str = "json") -> str:
    text = f"```{kind}\n{content}\n```" if kind else content
    return (
        f'??? info "{title}"\n\n'
        + "\n".join("    " + line if line else "" for line in text.splitlines())
        + "\n"
    )


def return_text(record: dict, language: str) -> str:
    name, entry = record["name"], record["entry"]
    if name == "dnakit.databases.ch_ipatlas_enrichment_analysis":
        return tr(
            "字典包含 status=requires_submission、analysis_completed=False 和提交说明；不会返回已完成的富集检验。下面保留提供方原始结果 Schema 供追溯，不能据此认定分析已经完成。",
            "Mapping with status=requires_submission, analysis_completed=False and submission instructions, not a completed enrichment test. The original provider schema is retained below for provenance and does not establish completed analysis.",
            language,
        )
    if name in RETURN_NOTES:
        return RETURN_NOTES[name][0 if language == "zh" else 1]
    if name.startswith("dnakit.download."):
        return tr(
            "DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](../features/15_download.md#accession-sequences)。",
            "DatasetDownloadResult with files, metadata (source, resolved accession and length), and manifest_path. See the [sequence download guide](../features/15_download.md#accession-sequences) for formats and coordinates.",
            language,
        )
    if name == "dnakit.workflows.write_wes_report":
        return tr(
            "返回报告路径和来源记录的字典；输出 wes_report.pdf 与 wes_report.manifest.json。需要 dnakit[reports]。",
            "Returns a mapping of report paths and provenance; writes wes_report.pdf and wes_report.manifest.json. Requires dnakit[reports].",
            language,
        )
    _, _, doc_return = doc_sections(record["fn"])
    if name.startswith("dnakit.encode."):
        return (
            tr("Python 返回类型为 str。", "Python return type: str. ", language) + " " + doc_return
        )
    if entry and entry["config"].get("return_schema"):
        return tr(
            "返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。",
            "Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.",
            language,
        )
    return tr(
        "返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。",
        "Returns object with the called operation's result structure; a common data envelope is not guaranteed. Provider errors raise exceptions.",
        language,
    )


def render_function(record: dict, language: str) -> str:
    name, entry = record["name"], record["entry"]
    output = [
        f"## `{name}` {{#{name}}}",
        "",
        purpose(record, language),
        "",
        "```text",
        name + signature(record["fn"]),
        "```",
        "",
        parameter_table(record, language),
        "",
        "**" + tr("返回值", "Returns", language) + "**",
        "",
        return_text(record, language),
        "",
    ]
    if entry:
        output += [
            details(
                tr(
                    "完整输入 Schema（含别名与条件约束）",
                    "Complete input schema, aliases and conditional constraints",
                    language,
                ),
                json.dumps(entry["input_schema"], ensure_ascii=False, indent=2),
            )
        ]
        schema = entry["config"].get("return_schema")
        if schema:
            output += [
                details(
                    tr("提供方登记的结果 Schema", "Registered provider result schema", language),
                    json.dumps(schema, ensure_ascii=False, indent=2),
                )
            ]
        extension = entry.get("extension")
        if extension:
            output += [
                details(
                    tr("DNAKit 参数扩展规则", "DNAKit argument extension rules", language),
                    json.dumps(extension, ensure_ascii=False, indent=2),
                )
            ]
        output += [
            tr("**执行与来源：**", "**Execution and provenance:**", language)
            + f" `{entry.get('effect', 'network')}`；"
            + tr("上游接口名", "upstream name", language)
            + f" `{entry['upstream_name']}`。",
            "",
        ]
        keys = entry["config"].get("required_api_keys", [])
        if keys:
            output += [
                tr("**提供方凭据：**", "**Provider credentials:**", language)
                + " "
                + ", ".join(f"`{key}`" for key in keys),
                "",
            ]
        output += [
            tr(
                "**查看参数示例（离线，不调用提供方）：**",
                "**Inspect parameters offline without calling the provider:**",
                language,
            ),
            "",
            "```python",
            "import dnakit",
            f'info = dnakit.databases.describe("{name.rsplit(".", 1)[1]}", verbose=False)',
            'print(info["input_schema"])',
            "```",
            "",
        ]
    else:
        before, _, _ = doc_sections(record["fn"])
        if language == "zh" and before != purpose(record, language):
            output += [details("补充说明（源码原文）", before, kind="")]
    return "\n".join(output)


def page_intro(group: str, language: str) -> str:
    common = tr(
        "每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。",
        "Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.",
        language,
    )
    if group == "encode":
        common += "\n\n" + tr(
            "ENCODE 接口返回文本：多数是 JSON，`get_citations` 可返回 BibTeX/RIS，`export_data` 可返回 CSV/TSV。仅 JSON 文本使用 `json.loads()`。追踪与凭据保存在 `~/.dnakit-encode/`；持久化凭据功能需要 `dnakit[encode-credentials]`。源码说明中的 `encode_...` 是保留的上游命名，在 DNAKit 中请使用本页的 `dnakit.encode.*` 函数。",
            "ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials]. Names starting with encode_ in source descriptions are upstream names; use this page's dnakit.encode.* functions in DNAKit.",
            language,
        )
    if group == "models":
        common += "\n\n" + tr(
            "这些入口调用已有模型或预测服务。模型权重、SDK、服务凭据和使用条件继续由原提供方决定；DNAKit 没有重新训练这些模型。",
            "These functions call existing models or prediction services. Weights, SDKs, credentials and service terms remain provider requirements; DNAKit did not retrain these models.",
            language,
        )
    if group == "regulation":
        common += "\n\n" + tr(
            "ChIP-Atlas 按 gene 查峰时，未提供 peak_file 会返回 requires_data 和检索计划；提供已核对组装的 BED 后才计算区间重叠。富集入口返回 requires_submission 提交说明，不能把它当成已计算的富集统计。参数扩展分支的结果可能不同于提供方原始 Schema。",
            "ChIP-Atlas gene-based peak queries without peak_file return requires_data and a retrieval plan; interval overlaps require a BED with verified assembly. The enrichment entry returns requires_submission instructions, not computed enrichment statistics. Extended query branches may return structures different from the original provider schema.",
            language,
        )
    if group == "execution":
        common += "\n\n" + tr(
            "scripts、resource_path、describe、search、list_operations 用于本地发现或资源读取。run、run_support 会执行脚本，可能联网或写文件；解释器、科学软件和参考数据需要按具体任务准备。原有严格 YAML 工作流 run_workflow 与这里的脚本执行入口是不同函数。",
            "scripts, resource_path, describe, search and list_operations support local discovery or resource reads. run and run_support execute scripts which may access networks or write files; prepare their interpreters, scientific software and reference data. The existing strict YAML run_workflow is separate from these script runners.",
            language,
        )
    return common


def render_page(group: str, records: list[dict], language: str) -> str:
    title = GROUPS[group][0 if language == "zh" else 1]
    lines = [
        f"# {title}",
        "",
        tr(
            f"本页收录 **{len(records)} 个新增公开函数**。",
            f"This page documents **{len(records)} newly added public functions**.",
            language,
        )
        + " "
        + tr(
            "[完整索引](index.md) · [调用约定与完整示例](usage.md)",
            "[Complete index](index.md) · [Conventions and examples](usage.md)",
            language,
        ),
        "",
        page_intro(group, language),
        "",
        tr("## 函数目录", "## Function index", language),
        "",
    ]
    lines += [f"- [`{r['name']}`](#{r['name']})" for r in records]
    lines += ["", *[render_function(r, language) for r in records]]
    return "\n".join(lines) + "\n"


def render_index(records: list[dict], language: str) -> str:
    lines = [
        "# " + tr("新增函数完整索引", "Complete index of added functions", language),
        "",
        tr(
            "本轮迁移共新增 **484 个公开 Python 函数**：数据库 455、ENCODE 20、工作流 6、序列检索与下载 3。每个函数在下表中只出现一次，并链接到完整参数与返回值说明。旧函数扩充参数、兼容别名、CLI/MCP 对同一函数的开放不重复计数；484 是接口数量，不是独立科研算法数。",
            "This migration adds **484 public Python functions**: 455 database, 20 ENCODE, 6 workflow, and 3 sequence retrieval/download functions. Each appears once below and links to its parameters and returns. Extended existing functions, aliases and CLI/MCP exposures are not counted again; 484 counts interfaces, not independent scientific algorithms.",
            language,
        ),
        "",
        tr(
            "[调用约定与完整示例](usage.md) · [功能树](../features/function_tree.md) · [API 总览](../index.md)",
            "[Conventions and examples](usage.md) · [Feature tree](../features/function_tree.md) · [API overview](../index.md)",
            language,
        ),
        "",
        tr(
            "可用网页搜索框搜索函数名、数据源或任务关键词，也可以使用浏览器查找。所有接口都能直接通过 Python 使用，不要求先安装 Skill。",
            "Search the site by function, provider or task, or use your browser's Find command. All these APIs can be called from Python without installing a Skill.",
            language,
        ),
        "",
        tr("## 按用途浏览", "## Browse by purpose", language),
        "",
        tr("| 用途分类 | 函数数量 |", "| Purpose | Function count |", language),
        "| --- | ---: |",
    ]
    for group, (zh, en, _) in GROUPS.items():
        members = [r for r in records if r["group"] == group]
        lines.append(f"| [{tr(zh, en, language)}]({group}.md) | {len(members)} |")
    lines += ["| **" + tr("合计", "Total", language) + "** | **484** |", ""]
    for group, (zh, en, _) in GROUPS.items():
        lines += [
            "## " + tr(zh, en, language),
            "",
            tr("| 函数 | 用途摘要 |", "| Function | Purpose summary |", language),
            "| --- | --- |",
        ]
        for r in records:
            if r["group"] == group:
                summary = purpose(r, language).split("\n")[0]
                summary = summary if len(summary) <= 220 else summary[:217] + "…"
                lines.append(f"| [`{r['name']}`]({group}.md#{r['name']}) | {cell(summary)} |")
        lines += [""]
    lines += [
        "## " + tr("来源与维护", "Provenance and maintenance", language),
        "",
        tr(
            "提供方说明、Schema 与上游别名来自当前源码内登记的接口定义，不表示已重新验证每个远程服务。DNAKit 的输入扩展、模型调用边界及仍需提交/下载的步骤见各页说明。上游来源和许可见[第三方声明](../../acknowledgements.md#third-party-notices)。",
            "Provider descriptions, schemas and aliases come from the current bundled definitions; this does not mean every remote service has been revalidated. See each page for argument extensions, model dependencies and steps requiring submission or data. See [third-party notices](../../acknowledgements.md#third-party-notices) for provenance and licensing.",
            language,
        ),
        "",
        tr(
            "文档从源码生成并逐项校验；维护命令在仓库根目录执行：",
            "Generate and check this reference from the source checkout:",
            language,
        ),
        "",
        "```bash",
        "python others/docs/build_native_api_docs.py",
        "python others/docs/build_native_api_docs.py --check",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    records = inventory()
    output: dict[Path, str] = {}
    for language in ("zh", "en"):
        suffix = ".md" if language == "zh" else ".en.md"
        output[DEST / ("index" + suffix)] = render_index(records, language)
        for group in GROUPS:
            output[DEST / (group + suffix)] = render_page(
                group, [r for r in records if r["group"] == group], language
            )
        # Each public name has exactly one detail anchor in each language.
        anchors = Counter(
            re.findall(
                r"\{#(dnakit\.[\w.]+)\}",
                "\n".join(
                    v
                    for p, v in output.items()
                    if p.name.endswith(suffix)
                    and (language == "en" or not p.name.endswith(".en.md"))
                ),
            )
        )
        assert anchors == Counter({r["name"]: 1 for r in records}), anchors
    manifest = {
        "scope": "Public functions added by the skill migration; aliases and existing-function extensions excluded",
        "total": len(records),
        "modules": {"databases": 455, "encode": 20, "workflows": 6, "sequence_retrieval": 3},
        "functions": [
            {
                "name": r["name"],
                "page": f"api/native/{r['group']}.md",
                "anchor": r["name"],
                "signature": signature(r["fn"]),
            }
            for r in records
        ],
    }
    output[ROOT / "registry/native-api-documentation.json"] = (
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    stale = []
    for path, content in output.items():
        if args.check:
            if not path.exists() or path.read_text() != content:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
    if stale:
        raise SystemExit("Regenerate stale documentation: " + ", ".join(stale))
    print(
        json.dumps(
            {
                "functions": 484,
                "languages": 2,
                "generated_pages": 26,
                "coverage": "484/484 unique detail entries per language",
                "mode": "check" if args.check else "write",
                "groups": Counter(r["group"] for r in records),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
