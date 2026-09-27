# ruff: noqa: E501, RUF001
"""Merge the audited user catalog into reproducible project discovery and FAQs.

This version keeps input claims, observed source evidence and actual local
entrypoints separate; existing functions are credited without counting them as
new implementations.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from rich.progress import track

CATEGORY_DOMAINS = {
    "01": ["breeding", "gwas"],
    "02": ["regulation"],
    "03": ["annotation"],
    "04": ["variants"],
    "05": ["clinical"],
    "06": ["crispr", "functional-assays"],
    "07": ["population", "comparative"],
    "08": ["metagenomics"],
    "09": ["synthetic"],
    "10": ["forensics", "traceability"],
    "11": ["epigenomics"],
    "12": ["sequence"],
    "13": ["models"],
    "14": ["workflows"],
    "15": ["workflows"],
}
EXTRA_DOMAINS = {
    "enformer": ["models"],
    "orca": ["3d-genomics", "models"],
    "basenji": ["models"],
    "borzoi": ["models"],
    "grelu": ["models"],
    "scenic": ["single-cell-dna"],
    "chromvar": ["single-cell-dna"],
    "moods": ["regulation"],
    "edta": ["repeats"],
    "repeatmasker": ["repeats"],
    "repeatmodeler": ["repeats"],
    "jcvi": ["comparative"],
    "orthofinder": ["comparative"],
    "quast": ["assembly"],
    "busco": ["assembly"],
    "ichorcna": ["variants"],
    "coloc": ["gwas"],
    "susier": ["gwas"],
    "saige": ["gwas"],
    "regenie": ["gwas"],
    "twosamplemr": ["gwas"],
    "msprime": ["simulation"],
    "stdpopsim": ["simulation"],
    "slim": ["simulation"],
    "tskit": ["simulation"],
    "dnabert2": ["models"],
    "dna-fountain": ["dna-storage"],
    "scadnano": ["nanotechnology"],
    "oxdna": ["nanotechnology"],
    "mtoolbox": ["organellar"],
    "expansionhunter": ["repeats"],
    "straglr": ["repeats"],
    "strling": ["repeats"],
    "trtools": ["repeats"],
    "genomad": ["viral-genomics"],
    "megahit": ["assembly"],
    "unicycler": ["assembly"],
    "canu": ["assembly"],
    "pggb": ["pangenome"],
    "vg": ["pangenome"],
    "odgi": ["pangenome"],
    "minigraph-cactus": ["pangenome"],
    "ppanggolin": ["pangenome"],
    "panaroo": ["pangenome"],
    "hicexplorer": ["3d-genomics"],
    "pairtools": ["3d-genomics"],
    "cooler": ["3d-genomics"],
    "cooltools": ["3d-genomics"],
    "archr": ["single-cell-dna", "epigenomics"],
    "signac": ["single-cell-dna", "epigenomics"],
    "fastqc": ["sequencing"],
    "multiqc": ["sequencing"],
    "cutadapt": ["sequencing"],
    "trim-galore": ["sequencing"],
    "bwa-mem2": ["alignment"],
    "bowtie2": ["alignment"],
    "minimap2": ["alignment"],
    "winnowmap": ["alignment"],
    "samtools": ["alignment"],
    "htslib": ["alignment"],
    "ncbi-datasets": ["databases"],
    "gget": ["databases"],
}

# These functions predate this import. They expose the specified upstream
# operation; they do not imply that every operation of that project is wrapped.
EXISTING_FUNCTIONS = {
    "biopython": [
        "dnakit.bioinformatics.reads.read_fastx",
        "dnakit.bioinformatics.alignment.read_alignment",
        "dnakit.bioinformatics.motifs.build_motif",
    ],
    "pysam": [
        "dnakit.bioinformatics.reads.alignment_summary",
        "dnakit.bioinformatics.reads.coverage_consensus",
    ],
    "pydna": [
        "dnakit.bioinformatics.cloning.assemble_duplexes",
        "dnakit.bioinformatics.cloning.read_constructs",
    ],
    "mafft": ["dnakit.bioinformatics.alignment.align_sequences"],
    "muscle": ["dnakit.bioinformatics.alignment.align_sequences"],
    "ncbi-blast": [
        "dnakit.bioinformatics.alignment.build_blast_database",
        "dnakit.bioinformatics.alignment.local_blast",
    ],
    "bowtie2": ["dnakit.bioinformatics.ngs.map_reads"],
    "samtools": ["dnakit.bioinformatics.ngs.prepare_bam"],
    "bcftools": ["dnakit.bioinformatics.ngs.call_variants"],
    "spades": ["dnakit.bioinformatics.ngs.assemble_reads"],
    "fastqc": ["dnakit.bioinformatics.ngs.quality_control"],
    "minimap2": ["dnakit.annotation.genes.annotate_genes"],
    "ensembl-vep": [
        "dnakit.annotation.variants.annotate_variant_vep",
        "dnakit.annotation.variants.annotate_rsid_vep",
    ],
    "enformer": [
        "dnakit.predictions.api.predict_sequence_properties",
        "dnakit.predictions.api.predict_enformer_benchmark",
    ],
    "evo2": ["dnakit.predictions.api.predict_properties"],
    "alphagenome-research": ["dnakit.predictions.api.predict_properties"],
    "nucleotide-transformer-agront-segmentnt": [
        "dnakit.predictions.api.predict_sequence_properties"
    ],
}
FUNCTION_NOTES = {
    "mafft": "align_sequences(method='mafft'); external MAFFT required.",
    "muscle": "align_sequences(method='muscle'); MUSCLE 5 required.",
    "bowtie2": "map_reads(method='bowtie2'); bowtie2 and bowtie2-build required.",
    "spades": "Existing SPAdes short-read/hybrid wrapper; not every metaSPAdes mode.",
    "minimap2": "Existing reference-based gene annotation/alignment wrapper; not every minimap2 CLI mode.",
    "ensembl-vep": "Existing Ensembl REST annotation through the configured ToolUniverse provider; not a local VEP CLI installation.",
    "enformer": "Existing configured Enformer prediction backend; local checkpoints required.",
    "evo2": "Existing configured Evo2 prediction backend; local checkpoints required.",
    "alphagenome-research": "Existing AlphaGenome research-weight backend; access and non-commercial model terms apply independently.",
    "nucleotide-transformer-agront-segmentnt": "Existing SegmentNT prediction backend only; not all NT/AgroNT models.",
}


def read(root: Path, name: str) -> dict[str, Any]:
    return json.loads((root / name).read_text(encoding="utf-8"))


def load_user_catalog(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    data = read(root, "registry/user-catalog-335.json")
    evidence = read(root, "registry/user-catalog-evidence.json")
    expected = hashlib.sha256(
        json.dumps(data, ensure_ascii=False, sort_keys=True).encode()
    ).hexdigest()
    if evidence.get("input_sha256") != expected:
        raise ValueError("Source evidence does not match the imported catalog; refresh evidence")
    for records, key in (
        (data["entries"], "input_stable_id"),
        (data["entries"], "project_id"),
        (evidence["pages"], "url"),
        (evidence["relationships"], "index"),
    ):
        if len({record[key] for record in records}) != len(records):
            raise ValueError(f"Duplicate catalog/evidence identity: {key}")
    if {r["index"] for r in evidence["relationships"]} != set(
        range(1, len(data["relationships"]) + 1)
    ):
        raise ValueError("Missing or extra citation relationship evidence")
    indexed = {r["index"]: r for r in evidence["relationships"]}
    for index, original in enumerate(data["relationships"], 1):
        for checked, supplied in (
            ("source", "来源"),
            ("relation_type", "关系类型"),
            ("target", "对象"),
            ("evidence_url", "证据入口"),
        ):
            if indexed[index][checked] != original[supplied]:
                raise ValueError("Citation evidence disagrees with the supplied relationship")
    context = root / "registry/user-catalog-context-review.json"
    if context.exists():
        pages = {page["url"]: page for page in evidence["pages"]}
        relations = {row["index"]: row for row in evidence["relationships"]}
        for review in json.loads(context.read_text())["relationships"]:
            current = relations[review["index"]]
            fields = ("source", "relation_type", "target", "evidence_url")
            if any(current[key] != review[key] for key in fields):
                raise ValueError("Citation relationship identity changed; review must be renewed")
            if pages[review["evidence_url"]].get("sha256") != review["evidence_sha256"]:
                raise ValueError("Citation source changed; context review must be renewed")
            current.update(copy.deepcopy(review))
            current["verification"] = (
                "context_reviewed"
                if review["relationship_verified"]
                else "planned_only_not_implemented"
            )
    return data, evidence


def external_guide(root: Path, *, english: bool) -> str:
    """Refresh only explicitly marked registry tables in the hand-written guide."""
    suffix = ".en" if english else ""
    content = (root / f"docs/api/features/39_external_genomics{suffix}.md").read_text()
    analyses = read(root, "src/dnakit/integrations/analyses.json")["analyses"]
    servers = read(root, "src/dnakit/integrations/mcp_servers.json")["servers"]
    tables = {
        "analyses": [
            "| Analysis ID / 分析 ID | Capability / 功能 | Backend / 后端 | Validation / 验证 |",
            "| --- | --- | --- | --- |",
            *[
                f"| `{a['id']}` | {cell(a['title'] if english else a.get('title_zh', a['title']))} | {cell(a['backend'])} | {cell(a.get('validation', {}).get('backend_execution', 'not_run'))} |"
                for a in analyses
            ],
        ],
        "mcp": [
            "| Server ID / 服务 ID | Upstream / 上游 | Configuration / 配置 |",
            "| --- | --- | --- |",
            *[
                f"| `{s['id']}` | [{cell(s['name'])}]({s['repository']}) | {cell(s.get('installation', 'See list_mcp_servers()'))} |"
                for s in servers
            ],
        ],
    }
    for tag, rows in tables.items():
        start = f"<!-- generated:{tag}:start -->"
        end = f"<!-- generated:{tag}:end -->"
        if content.count(start) != 1 or content.count(end) != 1:
            raise ValueError(f"Guide marker missing or duplicated: {tag}")
        left, rest = content.split(start)
        _, right = rest.split(end)
        content = left + start + "\n" + "\n".join(rows) + "\n" + end + right
    return content


def add_candidates(root: Path, projects: dict[str, dict[str, Any]]) -> None:
    data, evidence = load_user_catalog(root)
    pages = {page["url"]: page for page in evidence["pages"]}
    for row in data["entries"]:
        if row["existing_project_id"]:
            continue
        pid = row["project_id"]
        if pid in projects:
            raise ValueError(f"User catalog project identity collision: {pid}")
        repository = (
            row["project_url"] if row["project_url"].startswith("https://github.com/") else None
        )
        source = pages.get(row["evidence_url"], {})
        projects[pid] = {
            "id": pid,
            "name": row["name"],
            "domains": sorted(
                set(CATEGORY_DOMAINS[row["category"][:2]] + EXTRA_DOMAINS.get(pid, []))
            ),
            "kind": row["input_fields"]["条目类型"] + "; " + row["interface"],
            "homepage": row["project_url"],
            "repository": repository,
            "attribution": (repository or row["name"])
            + " authors and contributors; individual paper authors are not inferred",
            "paper_urls": [],
            "documentation_urls": [row["evidence_url"]],
            "capabilities": [row["purpose"]],
            "capability_evidence": "User-supplied project description; executable operations are separately enumerated below.",
            "license_evidence": {
                "status": "not_cleared_for_redistribution",
                "note": row["input_fields"]["许可限制 / 待审"],
            },
            "verified_urls": [source["retrieved_url"]]
            if source.get("status") == "retrieved"
            else [],
            "availability": "documentation_retrieved"
            if source.get("status") == "retrieved"
            else "unresolved_candidate",
            "validation_status": "source_entrypoint_checked; no_new_scientific_execution_asserted",
            "integration_recommendation": row["input_fields"]["安装路线"],
            "upstream_conditions": row["notes"],
            "research_part": "user-catalog-335.json",
        }


def merge_evidence(root: Path, projects: dict[str, dict[str, Any]]) -> None:
    data, evidence = load_user_catalog(root)
    pages = {page["url"]: page for page in evidence["pages"]}
    relations = {row["index"]: row for row in evidence["relationships"]}
    for row in data["entries"]:
        pid = row["project_id"]
        if pid not in projects:
            raise ValueError(f"Previously catalogued project was lost: {pid}")
        item = projects[pid]
        item["user_catalog"] = {
            "row": row["index"],
            "stable_id": row["input_stable_id"],
            "decision": "merged_existing" if row["existing_project_id"] else "added_project_record",
            "input_project_url": row["project_url"],
            "input_notes": row["notes"],
        }
        item["user_catalog_evidence"] = [
            copy.deepcopy(pages[url])
            for url in dict.fromkeys([row["project_url"], row["evidence_url"]])
            if url in pages
        ]
        proposed = [row["paper_url"]] if row["paper_url"] else []
        proposed += [p["文献入口"] for p in data["publications"] if p["关联项目"] == row["name"]]
        item["user_supplied_publication_urls"] = list(dict.fromkeys(proposed))
        # Preserve the distinction between an independently reviewed citation
        # and a link supplied by the user, even where they have the same DOI.
        item["citation_relationships"] = [
            {**copy.deepcopy(relation), "check": copy.deepcopy(relations[index])}
            for index, relation in enumerate(data["relationships"], 1)
            if any(
                name.strip() == row["name"] or name.strip().startswith(row["name"] + " ")
                for name in relation["来源"].split(";")
            )
        ]
        item["existing_local_functions"] = EXISTING_FUNCTIONS.get(pid, [])
        if pid in FUNCTION_NOTES:
            item["existing_function_scope"] = FUNCTION_NOTES[pid]
    # The supplied BOLE URL adds benchmark evidence, not a runnable framework.
    bole = projects["bole"]
    bole["repository"] = "https://github.com/ChongxiaoLiang/BOLE"
    bole["availability"] = (
        "public_benchmark_evidence; complete_framework_or_callable_api_not_verified"
    )
    bole["integration_recommendation"] = (
        "保留论文、服务与基准证据仓库；尚未核实完整可运行框架或公开 API，未伪造 BOLE 适配器。"
    )
    bole["validation_status"] = "paper_and_benchmark_readme_reviewed; framework_not_executed"
    bole["license_evidence"] = [
        {
            "url": bole["repository"],
            "note": "基准证据仓库与论文许可不等于完整软件框架可再分发或可调用。",
        }
    ]
    bole["verified_urls"] = list(dict.fromkeys([*bole["verified_urls"], bole["repository"]]))
    bole["user_catalog_evidence"].append(copy.deepcopy(pages[bole["repository"]]))


def validate_existing_functions(root: Path) -> None:
    for path in sorted(
        {function for functions in EXISTING_FUNCTIONS.values() for function in functions}
    ):
        module, function = path.rsplit(".", 1)
        source = root / "src" / (module.replace(".", "/") + ".py")
        names = {
            node.name
            for node in ast.parse(source.read_text()).body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        if function not in names:
            raise ValueError(f"Existing local function no longer exists: {path}")


def related_skills(
    root: Path, skills: list[dict[str, Any]], projects: dict[str, dict[str, Any]], repo_key: Any
) -> dict[str, list[str]]:
    """Match explicit upstream repository links, not vague keyword overlap."""
    result: dict[str, list[str]] = {pid: [] for pid in projects}
    keys = {
        pid: repo_key(p.get("repository")) for pid, p in projects.items() if p.get("repository")
    }
    for skill in track(skills, description="Match existing workflow references"):
        path = root / "skills" / skill["id"] / "SKILL.md"
        if not path.exists():
            continue
        body = path.read_text(encoding="utf-8")
        links = {
            repo_key(link.rstrip(".,;"))
            for link in re.findall(r"https://github\.com/[^\s<>)\]`\"']+", body)
        }
        for pid, key in keys.items():
            if any(link == key or link.startswith(key + "/") for link in links):
                result[pid].append(skill["id"])
    return result


def cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def faq(root: Path, survey: dict[str, Any], *, english: bool) -> str:
    data, evidence = load_user_catalog(root)
    projects = {p["id"]: p for p in survey["projects"]}
    local_skills = read(root, "registry/catalog.json")["skills"]
    analyses = read(root, "src/dnakit/integrations/analyses.json")["analyses"]
    servers = read(root, "src/dnakit/integrations/mcp_servers.json")["servers"]
    added_skills = read(root, "registry/user-catalog-skills.json")["skills"]
    lines = [
        "# " + ("Review of the 335 supplied projects" if english else "335 个补充项目的审核与接入"),
        "",
        (
            "All 335 rows were compared with DNAKit. The input contains software, models, workflows and MCP services, not 335 installable Skills."
            if english
            else "已逐条对照 DNAKit 检查全部 335 条。原始清单混合了软件、模型、工作流和 MCP，不能按 335 个可安装 Skill 计算。"
        ),
        "",
        f"- Existing project records / 已有项目记录：**{data['counts']['existing']}**；新增项目记录：**{data['counts']['new']}**。",
        f"- Total project/source records / 合并后项目与来源：**{len(projects)}**。",
        f"- Bundled Skills / 分发 Skills：**{len(local_skills)}**，本次新增 **{len(added_skills)}**。",
        f"- Executable analyses / 外部分析：**{len(analyses)}**（原 19，本次新增 **{len(analyses) - 19}**）；MCP connections / MCP 连接：**{len(servers)}**（原 7，本次新增 **{len(servers) - 7}**）。",
        f"- Source pages retrieved / 取得原始文档响应：**{sum(p['status'] == 'retrieved' for p in evidence['pages'])}/{len(evidence['pages'])}**（不等于运行验证）。",
        "",
        "## 1) " + ("How to use the additions" if english else "如何使用补充功能"),
        "",
        "```python",
        "from dnakit.integrations import list_analyses, describe_analysis, search_projects, list_mcp_servers",
        "print(list_analyses('DIAMOND'))",
        "print(describe_analysis('diamond-blastp'))",
        "print(search_projects('BOLE'))",
        "print(list_mcp_servers())",
        "```",
        "",
        (
            "Use plan_analysis/run_analysis for registered analyses and list_mcp_tools/call_mcp_tool for configured MCP services. Install upstream software and supply databases/checkpoints separately. See the external integrations guide."
            if english
            else "已注册分析通过 `plan_analysis`／`run_analysis` 调用；MCP 通过 `list_mcp_tools`／`call_mcp_tool` 调用。外部软件、参考数据库和模型权重按需配置。详见[外部基因组工具与 MCP](../api/features/39_external_genomics.md)。"
        ),
        "",
        "## 2) " + ("Decisions for every supplied row" if english else "335 条逐项处理结果"),
        "",
        "| # | Project / 项目 | Deduplication / 去重 | Local entrypoints / 本地能力 | Conditions / 条件与限制 |",
        "| ---: | --- | --- | --- | --- |",
    ]
    for row in data["entries"]:
        p = projects[row["project_id"]]
        integration = p["dnakit_integration"]
        interfaces = [
            *integration["native_functions"],
            *integration["analysis_adapters"],
            *integration["mcp_servers"],
            *p.get("existing_local_functions", []),
        ]
        if integration["bundled_skills"]:
            interfaces.append(f"{len(integration['bundled_skills'])} bundled Skills")
        if integration["related_local_skills"]:
            interfaces.append(f"{len(integration['related_local_skills'])} related workflow guides")
        state = "; ".join(dict.fromkeys(interfaces)) or (
            "Source only; no adapter" if english else "仅来源记录；未接入专用执行接口"
        )
        conditions = row["notes"]
        if p.get("existing_function_scope"):
            conditions += "; " + p["existing_function_scope"]
        lines.append(
            f"| {row['index']} | [{cell(row['name'])}]({row['project_url']}) | {'合并已有' if row['existing_project_id'] else '新增记录'} | {cell(state)} | {cell(conditions)} |"
        )
    lines += [
        "",
        "## 3) "
        + ("Citations, wrappers and dependencies" if english else "谁引用了谁：论文、封装与依赖"),
        "",
        (
            "The 86 supplied relationships retain their distinct types. A retrieved page or a matching name alone does not verify a relationship. Context-reviewed relationships are marked separately; other claims remain supplied evidence awaiting review."
            if english
            else "原表 86 条关系保留原始类型：论文引用、文档链接、软件依赖、工具封装、模型封装和数据复用分别记录。文档能访问或出现对象名称，不能单独证明关系成立；未完成上下文核对的关系明确标为待核实。"
        ),
        "",
        "| Source / 来源 | Relationship / 关系 | Target / 对象 | Evidence / 证据 | Verification / 核验 |",
        "| --- | --- | --- | --- | --- |",
    ]
    reviews = {r["index"]: r for r in evidence["relationships"]}
    for index, relation in enumerate(data["relationships"], 1):
        review = reviews[index]
        target = relation["对象"]
        target = f"[{cell(target)}]({target})" if target.startswith("http") else cell(target)
        status = (
            ("Context verified" if english else "上下文已核对")
            if review["relationship_verified"]
            else (
                "Planned description; implementation unverified"
                if english
                else "仅确认计划性描述；未确认实现"
            )
        )
        status += "; " + review.get("review_note", review["verification"])
        lines.append(
            f"| {cell(relation['来源'])} | `{relation['关系类型']}` | {target} | [source]({relation['证据入口']}) | {status} |"
        )
    lines += [
        "",
        "## 4) "
        + ("Added Skills and deferred components" if english else "新增 Skills 与暂不接入的组件"),
        "",
        "| Skill | Upstream / 来源 | Scope / 使用范围 |",
        "| --- | --- | --- |",
        *[
            f"| `{s['id']}` | [{cell(s.get('source_id', s.get('source', 'upstream')))}]({s['upstream_url']}) | {cell(s.get('description', s.get('capability', 'See SKILL.md')))} |"
            for s in added_skills
        ],
        "",
        "| Component / 组件 | Decision / 处理原因 |",
        "| --- | --- |",
        "| [Clair Skills](https://github.com/HKU-BAL/Clair-skills) | 技能内容的再分发许可未明确；保留来源，不复制。 |",
        "| [AWS Kiro microbiology](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/life-sciences-microbiology) | QIIME2 路径只返回 queued，未执行分析；未注册该 MCP 组件。 |",
        "| [AWS Kiro molecular biology](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/life-sciences-molbio) | 酶切位置计算、BLAST/Primer3 成功判断存在问题；未注册该 MCP 组件。 |",
        "| [Precision Medicine mcp-fgbio](https://github.com/lynnlangit/precision-medicine-mcp) | 源码返回固定模拟计数、由 hash 生成假 ENSG；未作为真实分析功能接入。 |",
        "| [Synapse marketplace connector](https://github.com/anthropics/life-sciences) | 需要 OAuth；当前 HTTP 桥接未实现该认证流程，保留来源。 |",
        "| [ClawBio MCP](https://github.com/ClawBio/ClawBio) | 上游已弃用并移除通用 MCP；保留已审查 Skills，不补失效入口。 |",
        "",
        "上述只列影响本轮接入决策的发现。每个外部分析的 validation、MCP 的 limitations，以及 registry/user-catalog-*-review.json 分别保存未执行条件与源码证据。",
        "",
        "## 5) " + ("Evidence boundaries and rebuild" if english else "证据边界与重建"),
        "",
        (
            "The public [BOLE repository](https://github.com/ChongxiaoLiang/BOLE) describes benchmark evidence and conversation logs. A complete runnable framework or public API has not been verified. DNAKit's ridge regression is not presented as a BOLE reproduction."
            if english
            else "BOLE 的[公开仓库](https://github.com/ChongxiaoLiang/BOLE)描述可复现证据、基准数据和对话日志；尚未核实完整可运行框架或公开 API。DNAKit 不将本地岭回归称为 BOLE 复现。"
        ),
        "",
        (
            "Skills require explicit licensing, complete resources and pinned revisions before import. Format adaptations retain the original text and attribution. Model code, weights and database terms apply separately. A launch recipe does not imply installation or scientific validation."
            if english
            else "上游 Skills 只在明确许可、完整资源和版本锁定后导入；格式适配会保留原文与来源。模型代码、权重和数据库许可分别适用。仅登记启动命令不表示在本机安装或完成科学验证。"
        ),
        "",
        (
            "The four supplied files remain unchanged. registry/user-catalog-335.json preserves their SHA-256 digests, all 335 project rows, 62 publication/archive links, 86 relationships, 38 coverage mappings and 52 skill-source claims. registry/user-catalog-evidence.json stores retrieval outcomes and content hashes."
            if english
            else "原始 temp 文件保持不变。`registry/user-catalog-335.json` 保存四文件 SHA-256、335 条原始字段、62 条文献入口、86 条关系、38 条功能映射和 52 条 Skills 来源记录；`registry/user-catalog-evidence.json` 保存读取结果与证据 SHA-256。"
        ),
        "",
        "```bash",
        "python others/skills/import_user_catalog.py --check",
        "python others/skills/build_genomics_survey.py --check",
        "```",
        "",
    ]
    return "\n".join(lines)
