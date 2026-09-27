"""Concise tutorial methods and explicitly illustrative API result shapes.

Shapes come from registered return schemas or the bundled implementation. They
are documentation, never captured service responses or predicted scientific data.
"""

# Bilingual prose and API labels intentionally retain their punctuation.
# ruff: noqa: E501, RUF001

from __future__ import annotations

import json
import re

METHODS = {
    "call": (
        "校验操作名和参数后调用对应接口，保留提供方返回的字段。",
        "Validate the operation and arguments, then dispatch to its adapter and preserve the returned fields.",
    ),
    "execute": (
        "解析操作名或 name/arguments 请求，再调用对应数据库接口。",
        "Parse an operation name or name/arguments request and dispatch to the database adapter.",
    ),
    "describe": (
        "从本地接口目录读取参数约束与来源信息，不联网。",
        "Read parameter constraints and provenance from the local operation catalog, without network access.",
    ),
    "search": (
        "在本地接口名称和描述中匹配关键词，按匹配数排序。",
        "Match keywords against local operation names and descriptions and rank by match count.",
    ),
    "list_operations": (
        "读取本地登记目录，列出函数名及上游别名。",
        "List function names and upstream aliases from the local registry.",
    ),
    "sequence_gc_content": (
        "统计 G、C 数量，除以序列长度减 N 数量得到 GC 比例。",
        "Divide the G+C count by sequence length minus the N count.",
    ),
    "sequence_reverse_complement": (
        "按 A↔T、C↔G、N↔N 互补后反转顺序。",
        "Complement A↔T, C↔G and N↔N, then reverse the sequence.",
    ),
    "sequence_count_residues": (
        "逐位统计目标残基，输出数量、比例及从 1 开始的位置。",
        "Count the selected residue and report its fraction and 1-based positions.",
    ),
    "sequence_stats": (
        "统计序列长度和字符组成，再按核酸或蛋白类型计算对应性质。",
        "Count sequence length and composition, then calculate nucleotide- or protein-specific properties.",
    ),
    "sequence_dn_ds": (
        "对同框密码子比对使用 Nei–Gojobori 估计及 Jukes–Cantor 校正。",
        "Apply the Nei–Gojobori estimator with Jukes–Cantor correction to in-frame codon-aligned sequences.",
    ),
    "coding_variant_fraction": (
        "以编码区变异为分母，统计满足注释及 VAF 阈值的比例。",
        "Use coding variants as the denominator and count those meeting the annotation and VAF threshold.",
    ),
    "dna_translate_reading_frames": (
        "按标准遗传密码表翻译三个正向阅读框，比较无终止密码子的连续片段。",
        "Translate the three forward frames with the standard genetic code and compare uninterrupted protein segments.",
    ),
    "dna_virtual_digest": (
        "匹配限制酶识别位点，按切割偏移拆分线性或环状序列。",
        "Match restriction sites and split linear or circular sequences at their cleavage offsets.",
    ),
    "dna_gibson_design": (
        "从相邻片段末端提取指定长度的同源区，添加到待组装片段。",
        "Take overlaps of the requested length from neighboring fragment ends and append them to assembly parts.",
    ),
    "dna_golden_gate_design": (
        "分配不回文的唯一四碱基突出端，并添加 Type IIS 识别位点。",
        "Assign unique nonpalindromic four-base overhangs and add Type IIS recognition sites.",
    ),
    "dna_golden_gate_assemble": (
        "模拟 Type IIS 酶切，移除仍含识别位点的片段，再按突出端连接。",
        "Simulate Type IIS digestion, discard fragments retaining recognition sites and ligate matching overhangs.",
    ),
    "vcf_count_variants": (
        "先用 bcftools view 过滤，再用 bcftools stats 统计变异记录。",
        "Filter with bcftools view, then count variant records with bcftools stats.",
    ),
    "vcf_summary_stats": (
        "对筛选后的 VCF 调用 bcftools stats，解析变异类型和样本统计。",
        "Run bcftools stats on the filtered VCF and parse variant types and sample statistics.",
    ),
    "vcf_normalize": (
        "用 bcftools norm 拆分或合并多等位位点；提供参考时左对齐，并比较前后统计。",
        "Use bcftools norm to split or join multiallelic sites, left-align when a reference is supplied, and compare counts.",
    ),
    "simulate_inheritance": (
        "逐独立位点从双亲各抽取一个等位基因，再按给定概率施加突变。",
        "Sample one allele from each parent at each independent locus, then apply the specified mutation probability.",
    ),
    "scripts": (
        "读取随包资源目录，列出指定 Skill 的脚本路径。",
        "Read bundled resources and list script paths for the selected skill.",
    ),
    "resource_path": (
        "定位随包资源并校验校验和，返回本地路径。",
        "Locate a bundled resource, verify its checksum and return its local path.",
    ),
    "run": (
        "使用指定解释器执行随包脚本，收集退出码、日志和来源记录。",
        "Execute the bundled script with the selected interpreter and collect its exit code, output and provenance.",
    ),
    "run_support": (
        "执行指定辅助脚本，收集退出码与日志。",
        "Execute the selected support script and collect its exit code and output.",
    ),
    "write_wes_report": (
        "将准备好的 Markdown 排版为 PDF，并写入来源清单。",
        "Render the prepared Markdown as a PDF and write a provenance manifest.",
    ),
    "read_executed_notebook": (
        "读取本地已执行 notebook 中保存的单元格输出，并按长度上限截取。",
        "Read saved cell outputs from a local executed notebook and apply the output-length limit.",
    ),
    "download_file": (
        "通过 HTTP(S) 分块写入指定路径，或在 return_content=True 时返回文本。",
        "Stream an HTTP(S) response to the destination, or return text when return_content=True.",
    ),
    "cellx_gene_download_h5ad": (
        "通过 CELLxGENE Census 定位数据集，下载对应的 H5AD 文件。",
        "Locate the dataset through CELLxGENE Census and download its H5AD file.",
    ),
}

# Field excerpts checked against the bundled adapters where no useful return
# schema is registered. Bracketed strings are type placeholders, not data.
SHAPES = {
    "describe": {
        "name": "<str>",
        "upstream_name": "<str>",
        "input_schema": "<dict>",
        "config": "<dict>",
    },
    "list_operations": [{"name": "<str>", "upstream_name": "<str>"}],
    "search": [{"name": "<str>", "description": "<str>", "upstream_name": "<str>"}],
    "pdb_get_structure": {
        "pdb_id": "<str>",
        "format": "<str>",
        "data": "<coordinate-file text>",
        "url": "<str>",
        "sha256": "<str>",
    },
    "scripts": ["<relative script path>"],
    "run": {
        "skill": "<str>",
        "script": "<str>",
        "returncode": "<int>",
        "output": "<str>",
        "upstream_url": "<str>",
    },
    "run_support": {"name": "<str>", "returncode": "<int>", "output": "<str>"},
    "simulate_inheritance": {
        "offspring": "<list>",
        "seed": "<int | None>",
        "mutation_rate": "<float>",
        "scope": "<str>",
        "upstream": "<str>",
    },
    "wiki_pathways_search": {
        "status": "success",
        "data": [{"wpid": "<str>", "uri": "<str>", "title": "<str>", "organism": "<str>"}],
        "metadata": {"query": "<str>", "count": "<int>", "backend": "WikiPathways SPARQL"},
    },
    "gdc_get_mutation_frequency": {
        "status": "success",
        "source": "GDC",
        "gene": "<str>",
        "data": {
            "gene_info": "<dict>",
            "total_ssm_occurrences": "<int>",
            "is_cancer_gene_census": "<bool | None>",
        },
    },
    "vcf_count_variants": {
        "status": "success",
        "data": "<variant counts and filters_applied>",
        "metadata": "<dict>",
    },
    "vcf_summary_stats": {
        "status": "success",
        "data": "<variant and sample statistics; filters_applied>",
        "metadata": "<dict>",
    },
    "vcf_normalize": {
        "status": "success",
        "data": {
            "before": "<variant counts>",
            "after": "<variant counts>",
            "left_aligned": "<bool>",
            "multiallelics": "<str>",
        },
        "metadata": "<dict>",
    },
}

ENCODE_RESULTS_ZH = {
    "batch_download": "返回 JSON 文本，包含下载预览（dry_run=True）或逐文件下载结果。",
    "compare_experiments": "返回 JSON 文本，包含兼容性结论、问题、警告和建议。",
    "download_files": "返回 JSON 文本，包含逐文件路径、大小、成功状态与 MD5 校验结果。",
    "export_data": "返回指定 CSV、TSV 或 JSON 格式的表格文本。",
    "get_citations": "返回指定格式的文献引用文本，例如 BibTeX 或 RIS。",
    "get_experiment": "返回 JSON 文本，包含实验详情与文件列表。",
    "get_facets": "返回 JSON 文本，包含筛选字段及各取值的记录数。",
    "get_file_info": "返回 JSON 文本，包含文件格式、大小、下载地址、MD5 和组装版本。",
    "get_metadata": "返回 JSON 文本，列出指定元数据类别的可用取值。",
    "get_provenance": "返回 JSON 文本，包含来源链或衍生文件列表。",
    "get_references": "返回 JSON 文本，列出关联的外部引用。",
    "link_reference": "返回 JSON 文本，报告引用关联结果。",
    "list_files": "返回 JSON 文本，列出文件 accession、格式、大小、下载地址及元数据。",
    "list_tracked": "返回 JSON 文本，包含已追踪实验的元数据表和追踪统计。",
    "log_derived_file": "返回 JSON 文本，包含新建来源记录的 ID。",
    "manage_credentials": "返回 JSON 文本，报告所选凭据操作的结果。",
    "search_experiments": "返回 JSON 文本，包含实验结果、总数和分页信息。",
    "search_files": "返回 JSON 文本，包含文件结果、总数和分页信息。",
    "summarize_collection": "返回 JSON 文本，按多个维度汇总实验数量。",
    "track_experiment": "返回 JSON 文本，包含追踪结果、找到的文献与分析流程信息。",
}

RESULT_NOTES = {
    "call": (
        "返回所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。",
        "Returns the selected operation's response; this example retrieves the Ensembl BRCA1 record.",
    ),
    "execute": (
        "返回请求中所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。",
        "Returns the requested operation's response; this example retrieves the Ensembl BRCA1 record.",
    ),
    "execute_tool": (
        "返回 tool_name 所选接口的结果，字段由该接口定义。",
        "Returns the result of the operation selected by tool_name, using that operation's fields.",
    ),
    "encode_search_experiments": (
        "返回 ENCODE 实验搜索结果；普通查询在 data 中保留响应，biosample 查询返回解析后的 ENCODE 搜索结果。",
        "Returns ENCODE experiment search results; ordinary queries preserve the response in data, while biosample queries return the parsed ENCODE search result.",
    ),
    "enrichr_enrich": (
        "返回富集条目、P 值、z 分数、combined score 和命中基因。",
        "Returns enriched terms, p-values, z-scores, combined scores and overlapping genes.",
    ),
    "gbif_search_species": (
        "字典的 data 包含 GBIF 物种搜索响应及分页信息。",
        "The data field contains the GBIF species search response and pagination information.",
    ),
    "gdc_get_clinical_data": (
        "data.cases 包含病例、诊断和治疗记录，data.pagination 包含分页信息。",
        "data.cases contains case, diagnosis and treatment records; data.pagination contains paging information.",
    ),
    "gdc_get_ssm_by_gene": (
        "data 包含 GDC 体细胞变异记录，gene 标明查询基因。",
        "data contains GDC somatic mutation records; gene identifies the queried gene.",
    ),
    "gdc_get_survival": (
        "data 包含供体随访时间、删失状态、生存估计和队列汇总；donors 最多展示 50 条。",
        "data contains follow-up times, censoring, survival estimates and cohort summaries; donors includes at most 50 records.",
    ),
    "gdc_list_files": (
        "data 中保留 GDC 文件列表及分页响应，记录含文件标识与元数据。",
        "data preserves the GDC file listing and pagination response, including file identifiers and metadata.",
    ),
    "gdc_list_projects": (
        "data 中保留 GDC 项目列表及项目统计。",
        "data preserves the GDC project listing and project statistics.",
    ),
    "gdc_search_cases": (
        "data 中保留 GDC 病例查询响应及分页信息。",
        "data preserves the GDC case search response and pagination information.",
    ),
    "gtex_get_expression_summary": (
        "返回解析到的基因标识及组织表达数据；无表达记录时会附 note。",
        "Returns the resolved gene identifier and tissue expression data; note explains missing expression records.",
    ),
    "gtex_query_eqtl": (
        "data.singleTissueEqtl 包含单组织 eQTL 记录，metadata 记录查询来源与条件。",
        "data.singleTissueEqtl contains single-tissue eQTL records; metadata records the source and query.",
    ),
    "marrvel_get_gene": (
        "返回基因标识交叉引用、名称、位置、别名与摘要；不包含 OMIM 编号。",
        "Returns gene identifier cross-references, names, location, aliases and summary; no OMIM identifier is included.",
    ),
    "mgnify_list_analyses": (
        "data 中包含指定研究的 MGnify 分析记录及服务返回的分页信息。",
        "data contains MGnify analysis records for the study and the service's pagination information.",
    ),
    "mgnify_search_studies": (
        "data 中包含 MGnify 研究检索结果及服务返回的分页信息。",
        "data contains MGnify study search results and the service's pagination information.",
    ),
    "monarch_v3_get_entity": (
        "data 包含实体名称、类别、描述、同义词和交叉引用。",
        "data contains the entity name, category, description, synonyms and cross-references.",
    ),
    "my_variant_get_pathogenicity_scores": (
        "返回 dbNSFP 已有变异评分，包括数据源提供的 REVEL、CADD、SIFT 等；缺失评分不会补算。",
        "Returns existing dbNSFP variant scores such as REVEL, CADD and SIFT where available; missing scores are not calculated.",
    ),
    "nvidia_nim_alphafold2": (
        "返回预测结构的 PDB 坐标及服务提供的置信度和模型指标。",
        "Returns predicted PDB coordinates and service-provided confidence and model metrics.",
    ),
    "openalex_get_author": (
        "返回 OpenAlex 作者记录；接口允许返回 null。",
        "Returns an OpenAlex author record; the declared response also allows null.",
    ),
    "pub_med_get_article": (
        "返回 PMID 对应文章的 XML 文本或已解析元数据，包含标题、摘要、作者和期刊信息。",
        "Returns article XML or parsed metadata for the PMID, including title, abstract, authors and journal information.",
    ),
    "rgd_get_gene": (
        "返回大鼠基因名称、功能描述、基因类型及 Ensembl 交叉引用。",
        "Returns rat gene names, functional description, gene type and Ensembl cross-references.",
    ),
    "usda_plants_get_profile": (
        "返回植物学名、俗名、分类、生长习性和区域原生状态；精细分布并非所有记录均有。",
        "Returns scientific/common names, taxonomy, growth habits and regional native status; fine-grained distributions are not always available.",
    ),
}


def translated(pair: tuple[str, str], language: str) -> str:
    return pair[0 if language == "zh" else 1]


def method_summary(record: dict, language: str) -> str:
    """Describe dispatch or calculation without pretending all APIs calculate."""
    name = record["name"].rsplit(".", 1)[-1]
    if name in METHODS:
        return translated(METHODS[name], language)
    if record["name"].startswith("dnakit.encode."):
        if name in {
            "get_provenance",
            "get_references",
            "get_citations",
            "list_tracked",
            "log_derived_file",
            "link_reference",
            "summarize_collection",
            "compare_experiments",
            "export_data",
            "manage_credentials",
        }:
            return translated(
                (
                    "读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。",
                    "Read or update local ENCODE tracking records and format the selected operation's result.",
                ),
                language,
            )
        return translated(
            (
                "按指定条件访问 ENCODE，解析实验、文件或下载结果。",
                "Query ENCODE with the supplied filters and parse experiment, file or download results.",
            ),
            language,
        )
    if record["name"].startswith("dnakit.download."):
        return translated(
            (
                "解析 accession，向对应序列服务请求数据，写出文件与来源清单。",
                "Resolve the accession, retrieve sequence data, and write the files and provenance manifest.",
            ),
            language,
        )
    if record["name"] == "dnakit.search.ncbi.accessions":
        return translated(
            (
                "用 NCBI ESummary 将 UID 转换成带版本的 accession。",
                "Resolve UIDs to versioned accessions using NCBI ESummary.",
            ),
            language,
        )
    entry = record["entry"] or {}
    config = entry.get("config", {})
    source = config.get("type") or entry.get("module", "").removesuffix("_tool")
    source = re.sub(r"(?:REST)?(?:Tool|API)$", "", source) or "API"
    if entry.get("effect") == "model":
        return translated(
            (
                f"通过 `{source}` 调用已训练模型，使用所选模型和输入返回预测。",
                f"Use `{source}` to run the selected trained model on the supplied input.",
            ),
            language,
        )
    return translated(
        (
            f"通过 `{source}` 按输入条件检索或提交请求，并解析响应。",
            f"Use `{source}` to query or submit the supplied input and parse the response.",
        ),
        language,
    )


def success_schema(schema: dict) -> dict:
    """Select a successful branch for the illustrative output excerpt."""
    for keyword in ("oneOf", "anyOf"):
        if keyword in schema:
            branches = [part for part in schema[keyword] if isinstance(part, dict)]
            successful = [part for part in branches if "error" not in part.get("required", [])]
            return success_schema(next(iter(successful or branches), {}))
    return schema


def schema_shape(schema: dict, depth: int = 0) -> object:
    """Render small, explicitly typed excerpts; never invent response values."""
    schema = success_schema(schema)
    kind = schema.get("type")
    if isinstance(kind, list):
        kind = (
            "array"
            if "array" in kind and "items" in schema
            else next((v for v in kind if v != "null"), "null")
        )
    if depth >= 3:
        return f"<{kind or 'value'}>"
    if kind == "array" or "items" in schema:
        items = schema.get("items", {})
        return [schema_shape(items if isinstance(items, dict) else {}, depth + 1)]
    if "properties" in schema:
        return {
            name: schema_shape(part, depth + 1)
            for name, part in list(schema["properties"].items())[:8]
        }
    if isinstance(schema.get("additionalProperties"), dict):
        return {"<record key>": schema_shape(schema["additionalProperties"], depth + 1)}
    return f"<{kind or 'value'}>"


def result_example(record: dict) -> str | None:
    """Return an illustrative shape only when its structure has a source."""
    name = record["name"].rsplit(".", 1)[-1]
    if name in SHAPES:
        shape = SHAPES[name]
    else:
        schema = success_schema((record["entry"] or {}).get("config", {}).get("return_schema", {}))
        if not schema:
            return None
        shape = schema_shape(schema)
        if shape in ("<object>", "<value>", {}):
            return None
    return json.dumps(shape, ensure_ascii=False, indent=2)
