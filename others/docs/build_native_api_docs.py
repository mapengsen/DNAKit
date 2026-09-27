"""Maintain public API details inside the existing functional tutorials.

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

from rich.console import Console
from rich.progress import track
from tutorial_api_examples import OVERRIDES as EXAMPLE_OVERRIDES
from tutorial_api_examples import arguments_for, render_example
from tutorial_api_layout import (
    ENCODE_RESULTS_ZH,
    RESULT_NOTES,
    method_summary,
    result_example,
    success_schema,
)
from tutorial_api_routes import FUNCTION_ROUTES, SECTIONS, page_for, section_for

from dnakit import databases, encode, workflows
from dnakit.databases._runtime import definitions
from dnakit.download import ena, ncbi
from dnakit.search.ncbi import accessions

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "docs/api/features"
REFERENCE = ROOT / "docs/faq/api"
TITLES = json.loads(
    Path(__file__).with_name("tutorial_api_titles.json").read_text(encoding="utf-8")
)

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
        "蛋白标识；输入列表或字符串类型见本行类型说明。",
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
        "返回 Orphanet_get_genes 的原始关联记录结构，见 [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes)。",
        "Returns the Orphanet_get_genes association result; see [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes).",
    ),
    "dnakit.databases.cadd_get_scores": (
        "返回 MyVariant 检索结果；字典结果另含 dnakit_score_source，表明分数来自已索引的 CADD 注释。缺失分数保持缺失，不执行新的模型预测。",
        "Returns MyVariant query results; mapping results include dnakit_score_source identifying indexed CADD annotations. Missing scores remain missing; no new model prediction is performed.",
    ),
    "dnakit.databases.string_get_functional_enrichment": (
        "返回 STRING_functional_enrichment 的富集结果，见 [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment)。",
        "Returns STRING_functional_enrichment results; see [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment).",
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
        if name in required or p.default is inspect.Signature.empty:
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
                "运行提示开关。" if name == "verbose" else "参数用法见本功能的作用与示例代码。",
                "Progress flag."
                if name == "verbose"
                else "See this function's purpose and example code for usage.",
                language,
            )
        limits = {k: v for k, v in prop.items() if k not in {"type", "description", "title"}}
        if limits:
            desc += " " + json.dumps(limits, ensure_ascii=False)
        rows.append(f"| `{name}` | {cell(kind)} | {cell(default)} | {cell(desc)} |")
    for keyword in ("oneOf", "anyOf"):
        alternatives = schema.get(keyword, [])
        if alternatives:
            groups = [
                " + ".join(f"`{name}`" for name in branch["required"]) for branch in alternatives
            ]
            rule = tr(
                "以下参数组必须且只能提供一组："
                if keyword == "oneOf"
                else "以下参数组至少提供一组：",
                "Provide exactly one of these parameter groups: "
                if keyword == "oneOf"
                else "Provide at least one of these parameter groups: ",
                language,
            )
            rows += ["", rule + "; ".join(groups) + tr("。", ".", language)]
    return "\n".join(rows)


def parameter_details(record: dict, language: str) -> str:
    table = parameter_table(record, language)
    keys = (record["entry"] or {}).get("config", {}).get("required_api_keys", [])
    credentials = [key for key in keys if key != "CELLXGENE_CENSUS_PACKAGE_INSTALLED"]
    if credentials:
        table += "\n\n" + tr("运行前配置：", "Configure before running: ", language)
        table += ", ".join(f"`{key}`" for key in credentials) + tr("。", ".", language)
    if "CELLXGENE_CENSUS_PACKAGE_INSTALLED" in keys:
        table += "\n\n" + tr(
            "运行前安装 `cellxgene-census`。",
            "Install `cellxgene-census` before running.",
            language,
        )
    return table + "\n"


def reference_link(record: dict, language: str) -> str:
    label = tr("完整参数与使用说明", "Full parameters and usage notes", language)
    return f"[{label}](../../faq/api/{page_for(record)}.md#{record['name']}-parameters)"


def return_text(record: dict, language: str) -> str:
    name = record["name"]
    if name == "dnakit.databases.download_file":
        return tr(
            "保存模式返回文件路径；`return_content=True` 时返回文件的文本内容。",
            "Returns the file path when saving; with `return_content=True`, returns the file contents as text.",
            language,
        )
    if name == "dnakit.databases.ch_ipatlas_enrichment_analysis":
        return tr(
            "字典包含 `status=requires_submission`、`analysis_completed=False` 和任务提交说明；本次调用尚未完成富集检验。",
            "Mapping with `status=requires_submission`, `analysis_completed=False` and task submission instructions; the enrichment test has not been run.",
            language,
        )
    if name in RETURN_NOTES:
        return RETURN_NOTES[name][0 if language == "zh" else 1]
    short_name = name.rsplit(".", 1)[-1]
    if short_name in RESULT_NOTES:
        return RESULT_NOTES[short_name][0 if language == "zh" else 1]
    if name.startswith("dnakit.download."):
        return tr(
            "DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](15_download.md#accession-sequences)。",
            "DatasetDownloadResult with files, metadata (source, resolved accession and length), and manifest_path. See the [sequence download guide](15_download.md#accession-sequences) for formats and coordinates.",
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
        if language == "zh":
            return ENCODE_RESULTS_ZH[name.rsplit(".", 1)[-1]]
        return (
            tr("Python 返回类型为 str。", "Python return type: str. ", language) + " " + doc_return
        )
    return tr(
        "返回该操作的处理结果，具体内容随输入与数据版本变化。",
        "Returns the operation's result; contents depend on the input and data version.",
        language,
    )


def input_summary(record: dict, language: str) -> str:
    entry = record["entry"] or {}
    required = set(entry.get("input_schema", {}).get("required", []))
    required.update(
        name
        for name, parameter in inspect.signature(record["fn"]).parameters.items()
        if parameter.default is inspect.Parameter.empty
    )
    names = ", ".join(
        f"`{name}`" for name in inspect.signature(record["fn"]).parameters if name in required
    )
    if names:
        return tr(
            f"必须提供 {names}。",
            f"Provide {names}.",
            language,
        )
    return tr(
        "按示例选择查询条件；可选参数使用接口默认值。",
        "Choose filters as shown in the example; omitted options use the operation defaults.",
        language,
    )


def input_requirements(record: dict, language: str) -> str:
    """Keep requirements needed to run the example beside its input."""
    entry = record["entry"] or {}
    schema = entry.get("input_schema", {})
    notes = []
    if record["name"] in {
        "dnakit.databases.sequence_count_residues",
        "dnakit.databases.sequence_stats",
    }:
        notes.append(
            tr(
                "另需提供 `sequence` 或 `uniprot_id`；使用 UniProt 编号时需要联网。",
                "Also provide sequence or uniprot_id; UniProt lookup requires network access.",
                language,
            )
        )
    for keyword in ("oneOf", "anyOf"):
        groups = [
            " + ".join(f"`{name}`" for name in branch.get("required", []))
            for branch in schema.get(keyword, [])
            if branch.get("required")
        ]
        if groups:
            notes.append(
                tr(
                    "参数组恰好选一组：" if keyword == "oneOf" else "参数组至少选一组：",
                    "Select exactly one argument group: "
                    if keyword == "oneOf"
                    else "Select at least one argument group: ",
                    language,
                )
                + " / ".join(groups)
                + tr("。", ".", language)
            )
    keys = entry.get("config", {}).get("required_api_keys", [])
    credentials = [key for key in keys if key != "CELLXGENE_CENSUS_PACKAGE_INSTALLED"]
    if credentials:
        notes.append(
            tr("需配置：", "Configure: ", language)
            + ", ".join(f"`{key}`" for key in credentials)
            + "."
        )
    if "CELLXGENE_CENSUS_PACKAGE_INSTALLED" in keys:
        notes.append(tr("需安装 `cellxgene-census`。", "Requires `cellxgene-census`.", language))
    if entry.get("module") == "vcf_stats_tool":
        notes.append(
            tr(
                "需本地 VCF/BCF 文件，且 PATH 中有 `bcftools`。",
                "Requires a local VCF/BCF file and `bcftools` on PATH.",
                language,
            )
        )
    return " " + " ".join(notes) if notes else ""


def example_notes(record: dict, language: str) -> str:
    name = record["name"]
    if name in {
        "dnakit.databases.hpa_get_protein_interactions_by_gene",
        "dnakit.databases.ctd_get_chemical_gene_interactions",
    }:
        return tr(
            "该数据源入口当前不可用，调用会返回错误；以下代码仅说明调用形式。",
            "This data-source operation is currently unavailable and returns an error; the code only illustrates the call syntax.",
            language,
        )
    notes = []
    if "REPLACE_ME" in repr(arguments_for(record)):
        notes.append(
            tr(
                "先用文件查询功能取得实际编号，替换 `ENCFF_REPLACE_ME`。",
                "Retrieve a real file accession with the file query API and replace `ENCFF_REPLACE_ME`.",
                language,
            )
        )
    if name.startswith("dnakit.encode.") and section_for(record) == "encode-collection":
        notes.append(
            tr(
                "集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。",
                "Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.",
                language,
            )
        )
    if name == "dnakit.databases.read_executed_notebook":
        notes.append(
            tr(
                "准备已经执行并保存输出的 `analysis_executed.ipynb`，或替换为自己的文件路径。",
                "Prepare an executed notebook with saved outputs at `analysis_executed.ipynb`, or supply your own path.",
                language,
            )
        )
    if name == "dnakit.workflows.write_wes_report":
        notes.append(
            tr(
                "先准备 `report.md` 并安装 `dnakit[reports]`；结果将写入 `reports/pdf`。",
                "Prepare `report.md` and install `dnakit[reports]`; output is written to `reports/pdf`.",
                language,
            )
        )
    if name == "dnakit.workflows.simulate_inheritance":
        notes.append(
            tr(
                "使用显式合成位点演示独立遗传；关闭突变，不访问远程数据。",
                "This example uses an explicit synthetic locus with mutation disabled and no remote data.",
                language,
            )
        )
    effect = (record["entry"] or {}).get("effect")
    if effect in {"filesystem", "subprocess", "external_process"}:
        notes.append(
            tr(
                "运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。",
                "Prepare the local files and external programs used by the example, and replace paths with your input or output locations.",
                language,
            )
        )
    return " ".join(notes)


def result_summary(record: dict, language: str) -> str:
    name = record["name"]
    if name in {
        "dnakit.databases.hpa_get_protein_interactions_by_gene",
        "dnakit.databases.ctd_get_chemical_gene_interactions",
    }:
        return tr(
            "此入口当前报告服务不可用，不返回有效互作结果。",
            "This operation currently reports an unavailable service rather than valid interaction results.",
            language,
        )
    schema = success_schema((record["entry"] or {}).get("config", {}).get("return_schema", {}))
    text = return_text(record, language)
    if text.startswith(("返回该操作的处理结果", "Returns the operation's result")):
        if result_example(record):
            return tr(
                "返回结构如下，实际记录随输入和数据源变化。",
                "The response has the following shape; records depend on the input and source.",
                language,
            )
        description = schema.get("description")
        if description and language == "en":
            return description
    return text


def output_summary(record: dict, language: str) -> str:
    """Describe the Python return value before the worked example."""
    name = record["name"]
    if name.startswith("dnakit.download."):
        return tr(
            "返回 `DatasetDownloadResult`，包含文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。",
            "Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.",
            language,
        )
    if name in {
        "dnakit.databases.hpa_get_protein_interactions_by_gene",
        "dnakit.databases.ctd_get_chemical_gene_interactions",
    }:
        return result_summary(record, language)

    description = return_text(record, language)
    if not description.startswith(("返回该操作的处理结果", "Returns the operation's result")):
        return description

    example = result_example(record)
    if example:
        shape = json.loads(example)
        if isinstance(shape, list):
            item = shape[0] if shape else None
            if isinstance(item, dict):
                fields = ", ".join(f"`{key}`" for key in list(item)[:4])
                return tr(
                    f"返回记录列表；每条记录包含 {fields} 等字段。",
                    f"Returns a list of records; each record includes {fields} and related fields.",
                    language,
                )
            title = TITLES[name][language]
            return tr(
                f"返回“{title}”的列表；元素结构见下方示例。",
                f"Returns a list for “{title}”; the item shape is shown below.",
                language,
            )
        if isinstance(shape, dict):
            if not shape:
                return tr(
                    "返回字典；具体字段随请求模式变化。",
                    "Returns a mapping whose fields depend on the request mode.",
                    language,
                )
            if "<record key>" in shape:
                item = shape["<record key>"]
                if isinstance(item, dict):
                    fields = ", ".join(f"`{key}`" for key in list(item)[:4])
                    return tr(
                        f"返回以记录键为索引的字典；每条记录包含 {fields} 等字段。",
                        f"Returns a mapping keyed by record; each record includes {fields} and related fields.",
                        language,
                    )
                return tr(
                    "返回以记录键为索引的字典。",
                    "Returns a mapping keyed by record.",
                    language,
                )
            fields = ", ".join(f"`{key}`" for key in list(shape)[:4])
            return tr(
                f"返回字典，包含 {fields} 等字段。",
                f"Returns a mapping with {fields} and related fields.",
                language,
            )
        if isinstance(shape, str):
            kind = shape.strip("<>")
            kinds = {
                "string": ("字符串", "string"),
                "integer": ("整数", "integer"),
                "number": ("数值", "number"),
                "boolean": ("布尔值", "boolean"),
            }
            if kind in kinds:
                zh, en = kinds[kind]
                return tr(f"返回{zh}。", f"Returns a {en}.", language)

    schema = success_schema((record["entry"] or {}).get("config", {}).get("return_schema", {}))
    if schema.get("description"):
        return schema["description"]
    return tr(
        "返回查询结果；字段随数据源响应变化。",
        "Returns query results; fields depend on the source response.",
        language,
    )


def render_function(record: dict, language: str, number: int) -> str:
    name = record["name"]
    title = TITLES[name][language]
    summary = title + ("。" if language == "zh" else ".")
    if language == "zh" and name in ZH_PURPOSE:
        summary = ZH_PURPOSE[name]
    if name in PURPOSE_OVERRIDES:
        summary = PURPOSE_OVERRIDES[name][0 if language == "zh" else 1]
    if name == "dnakit.databases.alpha_missense_get_protein_scores":
        summary = tr(
            "按 UniProt 编号汇总蛋白质各位点的替换分类、位点平均分数与阈值；可用 max_residues 限制查询位点数。",
            "Aggregate substitution classes, residue mean scores and thresholds across a UniProt protein; max_residues limits the queried positions.",
            language,
        )
    elif name == "dnakit.databases.alpha_missense_get_residue_scores":
        summary = tr(
            "查看指定位点的良性、不确定和致病替换分类，以及位点平均分数；返回值不是每种替换各自的数值分数。",
            "Retrieve benign, ambiguous and pathogenic substitution classes plus residue mean scores; these are not numeric scores for each individual substitution.",
            language,
        )
    elif name == "dnakit.databases.alpha_missense_get_variant_score":
        summary = tr(
            "查询单个错义变异的 AlphaMissense 分类。上游不提供该替换的独立数值分数，pathogenicity_score 为 null，score_available 为 false。",
            "Retrieve an AlphaMissense classification for one missense variant. The upstream operation supplies no numeric substitution score: pathogenicity_score is null and score_available is false.",
            language,
        )
    output = [
        f"## {number}) {title} {{#{name}}}",
        "",
        "- **" + tr("作用：", "Purpose:", language) + "** " + summary,
        "- **" + tr("处理方法：", "Method:", language) + "** " + method_summary(record, language),
        "- **API：** `" + name + signature(record["fn"]) + "`",
        "- **"
        + tr("输入：", "Input:", language)
        + "** "
        + input_summary(record, language)
        + input_requirements(record, language),
    ]
    note = example_notes(record, language)
    if note:
        output[-1] += " " + note
    output += [
        "- **" + tr("输出：", "Output:", language) + "** " + output_summary(record, language),
        "- **" + tr("示例代码：", "Example code:", language) + "**",
        "",
        "```python",
        render_example(record),
        "```",
        "",
        "- **"
        + tr("示例结果：", "Example output:", language)
        + "** "
        + result_summary(record, language),
    ]
    example = result_example(record)
    if example and name not in {
        "dnakit.databases.hpa_get_protein_interactions_by_gene",
        "dnakit.databases.ctd_get_chemical_gene_interactions",
    }:
        output += [
            "",
            tr(
                "返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。",
                "Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.",
                language,
            ),
            "",
            "```json",
            example,
            "```",
        ]
    output += ["", reference_link(record, language), ""]
    return re.sub(r"(?m)^[ \t]+$", "", "\n".join(output))


def page_intro(group: str, language: str) -> str:
    common = tr(
        "各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。",
        "Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.",
        language,
    )
    if group == "encode":
        common += "\n\n" + tr(
            "ENCODE 接口返回文本：多数是 JSON，`get_citations` 可返回 BibTeX/RIS，`export_data` 可返回 CSV/TSV。仅 JSON 文本使用 `json.loads()`。追踪与凭据保存在 `~/.dnakit-encode/`；持久化凭据功能需要 `dnakit[encode-credentials]`。",
            "ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials].",
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
            "ChIP-Atlas 按 gene 查峰时，未提供 peak_file 会返回 requires_data 和检索计划；提供已核对组装的 BED 后才计算区间重叠。富集入口返回 requires_submission 提交说明，不能把它当成已计算的富集统计。",
            "ChIP-Atlas gene-based peak queries without peak_file return requires_data and a retrieval plan; interval overlaps require a BED with verified assembly. The enrichment entry returns requires_submission instructions, not computed enrichment statistics.",
            language,
        )
    if group == "execution":
        common += "\n\n" + tr(
            "scripts、resource_path、describe、search、list_operations 用于本地发现或资源读取。run、run_support 会执行脚本，可能联网或写文件；解释器、科学软件和参考数据需要按具体任务准备。原有严格 YAML 工作流 run_workflow 与这里的脚本执行入口是不同函数。",
            "scripts, resource_path, describe, search and list_operations support local discovery or resource reads. run and run_support execute scripts which may access networks or write files; prepare their interpreters, scientific software and reference data. The existing strict YAML run_workflow is separate from these script runners.",
            language,
        )
    return common


def render_section(
    section: str, records: list[dict], language: str, numbers: dict[str, int]
) -> str:
    lines = [
        f'<span id="{section}"></span>',
        "",
    ]
    lines += [render_function(r, language, numbers[r["name"]]) for r in records]
    return "\n".join(lines).rstrip()


def render_reference(page: str, records: list[dict], language: str) -> str:
    suffix = ".en.md" if language == "en" else ".md"
    title = (DEST / (page + suffix)).read_text(encoding="utf-8").splitlines()[0].lstrip("# ")
    lines = [
        "# " + title + tr("：参数与补充说明", ": parameters and notes", language),
        "",
        tr(
            f"用法和示例见[功能教程](../../api/features/{page}.md)。",
            f"See the [tutorial](../../api/features/{page}.md) for usage and examples.",
            language,
        ),
        "",
    ]
    for group in sorted({record["group"] for record in records}):
        lines += [page_intro(group, language), ""]
    for record in records:
        lines += [
            f"## {TITLES[record['name']][language]} {{#{record['name']}-parameters}}",
            "",
            f"`{record['name']}{signature(record['fn'])}`",
            "",
            parameter_details(record, language),
            "",
        ]
    return "\n".join(lines).rstrip() + "\n"


def update_section(content: str, section: str, rendered: str) -> str:
    """Replace only the managed block, preserving handwritten tutorial prose."""
    start = f"<!-- BEGIN GENERATED API: {section} -->"
    end = f"<!-- END GENERATED API: {section} -->"
    block = f"{start}\n\n{rendered}\n\n{end}"
    if start not in content and end not in content:
        return content.rstrip() + "\n\n" + block + "\n"
    if content.count(start) != 1 or content.count(end) != 1:
        raise ValueError(f"Invalid generated section markers: {section}")
    before, remainder = content.split(start)
    _, after = remainder.split(end)
    return before + block + after


def validate_coverage(output: dict[Path, str], records: list[dict]) -> None:
    """Check all tutorials, including unmanaged prose, for duplicate details."""
    expected = Counter({r["name"]: 1 for r in records})
    for language in ("zh", "en"):
        texts = [
            output.get(path, path.read_text(encoding="utf-8"))
            for path in (ROOT / "docs").rglob("*.md")
            if path.name.endswith(".en.md") == (language == "en")
        ]
        anchors = Counter(
            name
            for name in re.findall(r"\{#(dnakit\.[\w.]+)\}", "\n".join(texts))
            if name in expected
        )
        if anchors != expected:
            raise ValueError(
                f"{language} missing or duplicate API details: {anchors - expected}; missing: {expected - anchors}"
            )


def function_numbers(records: list[dict]) -> dict[str, int]:
    """Number generated APIs where their blocks occur among handwritten tasks."""
    numbers = {}
    entries = re.compile(
        r"^ {0,3}(?P<backticks>`{3,})[^\n]*\n.*?^ {0,3}(?P=backticks)`*[ \t]*$"
        r"|^ {0,3}(?P<tildes>~{3,})[^\n]*\n.*?^ {0,3}(?P=tildes)~*[ \t]*$"
        r"|<!-- BEGIN GENERATED API: (?P<section>[\w-]+) -->.*?"
        r"<!-- END GENERATED API: (?P=section) -->"
        r"|^#{2,3} (?P<number>\d+)[).、]\s+",
        flags=re.M | re.S,
    )
    for page in sorted({page_for(record) for record in records}):
        sections = {
            section: [record for record in records if section_for(record) == section]
            for section, (section_page, _, _) in SECTIONS.items()
            if section_page == page
        }
        bilingual_numbers = []
        for suffix in (".md", ".en.md"):
            text = (DEST / (page + suffix)).read_text(encoding="utf-8")
            page_numbers = {}
            seen = set()
            number = 1
            for entry in entries.finditer(text):
                if entry.group("backticks") is not None or entry.group("tildes") is not None:
                    continue
                section, handwritten_number = entry.group("section", "number")
                if section is None:
                    number = max(number, int(handwritten_number) + 1)
                    continue
                if section not in sections or section in seen:
                    raise ValueError(f"Invalid generated section in {page + suffix}: {section}")
                seen.add(section)
                for record in sections[section]:
                    page_numbers[record["name"]] = number
                    number += 1
            # Newly routed sections are appended in the same order as update_section.
            for section, members in sections.items():
                if section in seen:
                    continue
                for record in members:
                    page_numbers[record["name"]] = number
                    number += 1
            bilingual_numbers.append(page_numbers)
        if bilingual_numbers[0] != bilingual_numbers[1]:
            raise ValueError(f"Tutorial numbering differs between languages: {page}")
        numbers.update(bilingual_numbers[0])
    return numbers


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    records = inventory()
    names = {r["name"] for r in records}
    if set(TITLES) != names:
        raise ValueError(
            f"Title coverage mismatch: missing {names - TITLES.keys()}, extra {TITLES.keys() - names}"
        )
    if set(EXAMPLE_OVERRIDES) - names:
        raise ValueError("Example overrides refer to unknown functions")
    for name, titles in TITLES.items():
        if set(titles) != {"zh", "en"} or any(
            not value.strip() or "dnakit." in value for value in titles.values()
        ):
            raise ValueError(f"Invalid functional title: {name}")
    numbers = function_numbers(records)
    unknown_routes = set(FUNCTION_ROUTES) - {r["name"] for r in records}
    if unknown_routes:
        raise ValueError(f"Routes refer to missing functions: {sorted(unknown_routes)}")
    retired = ROOT / "docs/api/native"
    if retired.exists() and any(retired.glob("*.md")):
        raise ValueError("Move legacy api/native pages into the tutorials before regenerating.")
    output: dict[Path, str] = {}
    for language in ("zh", "en"):
        suffix = ".md" if language == "zh" else ".en.md"
        for section, (page, _, _) in SECTIONS.items():
            members = [r for r in records if section_for(r) == section]
            if not members:
                raise ValueError(f"Empty tutorial section: {section}")
            path = DEST / (page + suffix)
            content = output.get(path)
            if content is None:
                content = path.read_text(encoding="utf-8")
            output[path] = update_section(
                content, section, render_section(section, members, language, numbers)
            )
        reference_links = [
            "# " + tr("接口参数与补充说明", "API parameters and usage notes", language),
            "",
            tr(
                "按功能查阅完整签名、可选参数、凭据与运行约定。直接使用请返回对应教程中的示例。",
                "Browse full signatures, optional parameters, credentials and execution conventions by task. Each page links back to its tutorial examples.",
                language,
            ),
            "",
        ]
        for page in sorted({page_for(record) for record in records}):
            members = [record for record in records if page_for(record) == page]
            output[REFERENCE / (page + suffix)] = render_reference(page, members, language)
        # Preserve handwritten parameter pages alongside the generated references.
        reference_pages = {
            path
            for path in REFERENCE.glob("*.md")
            if path.name not in {"index.md", "index.en.md"}
            and path.name.endswith(".en.md") == (language == "en")
        }
        reference_pages.update(
            path
            for path in output
            if path.parent == REFERENCE
            and path.name not in {"index.md", "index.en.md"}
            and path.name.endswith(".en.md") == (language == "en")
        )
        for path in sorted(reference_pages):
            content = output[path] if path in output else path.read_text(encoding="utf-8")
            title = content.splitlines()[0].lstrip("# ")
            page = path.name.removesuffix(suffix)
            reference_links.append(f"- [{title}]({page}.md)")
        output[REFERENCE / ("index" + suffix)] = "\n".join(reference_links) + "\n"
    validate_coverage(output, records)
    manifest = {
        "scope": "Public functions added by the skill migration; aliases and existing-function extensions excluded",
        "total": len(records),
        "modules": {"databases": 455, "encode": 20, "workflows": 6, "sequence_retrieval": 3},
        "functions": [
            {
                "name": r["name"],
                "page": f"api/features/{page_for(r)}.md",
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
    for path, content in track(
        output.items(),
        description="Check documentation" if args.check else "Write documentation",
        console=Console(stderr=True),
        transient=True,
    ):
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    if stale:
        raise SystemExit("Regenerate stale documentation: " + ", ".join(stale))
    print(
        json.dumps(
            {
                "functions": len(records),
                "languages": 2,
                "tutorial_pages_per_language": len({page_for(r) for r in records}),
                "coverage": "484/484 unique detail entries per language",
                "mode": "check" if args.check else "write",
                "numbered_functional_titles": len(TITLES),
                "validated_call_examples": len(records),
                "tutorials": Counter(page_for(r) for r in records),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
