#!/usr/bin/env python3
"""Build the dated project survey, packaged search index and source FAQ.

This version merges primary-source research with pinned existing skill sources,
and derives integration status from actual local registries. It never executes
scientific backends. Run --check to detect stale generated files.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import user_catalog_support as user_catalog
from rich.progress import track

ROOT = Path(__file__).resolve().parents[2]
DATE = "2026-09-20"
DOMAINS = {
    "sequence": "序列处理与格式",
    "sequencing": "测序、质控与覆盖度",
    "alignment": "比对与参考序列",
    "assembly": "基因组组装",
    "pangenome": "泛基因组与图基因组",
    "annotation": "基因识别与功能注释",
    "regulation": "启动子、增强子与调控",
    "variants": "小变异与结构变异",
    "clinical": "疾病、肿瘤与药物基因组",
    "breeding": "育种、性状与种质资源",
    "gwas": "关联分析、遗传力与因果遗传",
    "population": "群体遗传与祖源",
    "comparative": "比较基因组与系统发育",
    "simulation": "遗传与进化模拟",
    "polyploid": "多倍体与单倍型",
    "forensics": "法医、亲缘与 DNA 鉴定",
    "traceability": "物种、食品与环境溯源",
    "ancient-dna": "古 DNA 与损伤",
    "edna": "环境 DNA 与保护遗传",
    "metagenomics": "微生物、宏基因组与耐药",
    "epigenomics": "甲基化与表观遗传",
    "3d-genomics": "染色质空间结构与接触图",
    "single-cell-dna": "单细胞 DNA 与多组学",
    "spatial-genomics": "空间基因组",
    "crispr": "基因编辑与脱靶",
    "synthetic": "合成生物学与序列优化",
    "replication-repair": "复制、损伤修复与端粒",
    "functional-assays": "功能筛选与变异效应实验",
    "repeats": "重复与移动遗传元件",
    "viral-genomics": "病毒与噬菌体基因组",
    "organellar": "线粒体与叶绿体基因组",
    "nanotechnology": "DNA 纳米结构与分子力学",
    "dna-storage": "DNA 数据存储与编码",
    "strand-displacement": "DNA 链置换与分子计算",
    "models": "DNA 基础模型与表示",
    "databases": "数据库与知识检索",
    "workflows": "智能体、Skills、MCP 与复现",
}
ALIASES = {
    "育种与农业": "breeding",
    "种质资源": "breeding",
    "抗病性状": "breeding",
    "DNA鉴定": "forensics",
    "亲子鉴定": "forensics",
    "家系": "forensics",
    "法医学": "forensics",
    "古DNA": "ancient-dna",
    "环境DNA": "edna",
    "保护遗传": "edna",
    "多倍体": "polyploid",
    "局部祖源": "population",
    "祖源": "population",
    "群体遗传": "population",
    "进化仿真": "simulation",
    "遗传仿真": "simulation",
    "食品溯源": "traceability",
    "物种鉴定": "traceability",
    "GWAS": "gwas",
    "变异分析": "variants",
    "序列分析": "sequence",
    "全基因组比对": "alignment",
    "泛基因组": "pangenome",
    "数据库": "databases",
    "疾病与医学": "clinical",
    "medical": "clinical",
    "gene-function": "annotation",
    "系统发育": "comparative",
    "重复序列": "repeats",
    "intervals": "sequence",
}
NATIVE = {
    "vanraden-grm": ["genomic_relationship_matrix"],
    "rrblup": ["genomic_prediction"],
    "pyrodigal": ["predict_prokaryotic_genes"],
    "modkit": ["read_bedmethyl", "summarize_methylation"],
    "chamaeleo": ["encode_dna_storage", "decode_dna_storage"],
    "scadnano": ["design_dna_duplex", "export_oxdna_design"],
}


def read(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def repo_key(value: str | None) -> str:
    key = (
        (value or "")
        .removeprefix("https://github.com/")
        .removesuffix(".git")
        .rstrip("/")
        .casefold()
    )
    return {
        "k-dense-ai/claude-scientific-skills": "k-dense-ai/scientific-agent-skills",
        "weberlab-hhu/helixer": "usadellab/helixer",
        "fenderglass/flye": "mikolmogorov/flye",
    }.get(key, key)


def domains_for_skill(skill: dict[str, Any]) -> set[str]:
    values = skill.get("categories", [skill.get("category", "workflows")])
    result = {ALIASES.get(x, x) for x in values if x}
    text = json.dumps(skill, ensure_ascii=False).casefold()
    for needle, domain in (
        ("mitochond", "organellar"),
        ("chloroplast", "organellar"),
        ("phage", "viral-genomics"),
        ("viral genome", "viral-genomics"),
        ("transpos", "repeats"),
        ("mobile genetic", "repeats"),
        ("pangenome", "pangenome"),
        ("ancient dna", "ancient-dna"),
        ("edna", "edna"),
        ("polyploid", "polyploid"),
    ):
        if needle in text:
            result.add(domain)
    return result & DOMAINS.keys()


def build() -> dict[str, Any]:
    projects: dict[str, dict[str, Any]] = {}
    queries = []
    limitations = []
    parts = sorted((ROOT / "registry/genomics-survey-parts").glob("*.json"))
    for path in track(parts, description="Merge research evidence"):
        part = json.loads(path.read_text(encoding="utf-8"))
        queries.extend({**q, "research_part": path.name} for q in part["queries"])
        limitations.extend(part.get("limitations", []))
        for original in part["candidates"]:
            item = copy.deepcopy(original)
            item["original_domains"] = item["domains"]
            item["domains"] = sorted({ALIASES.get(x, x) for x in item["domains"]})
            item["research_part"] = path.name
            if item["id"] in projects:
                raise ValueError(f"Duplicate project id: {item['id']}")
            projects[item["id"]] = item
    has_user_catalog = (ROOT / "registry/user-catalog-335.json").exists()
    if has_user_catalog:
        user_catalog.add_candidates(ROOT, projects)
    catalog = read("registry/catalog.json")
    by_repository: dict[str, list[dict[str, Any]]] = {}
    for p in projects.values():
        if p.get("repository"):
            by_repository.setdefault(repo_key(p["repository"]), []).append(p)
    for source in catalog["sources"]:
        sid = source["id"]
        if sid == "dnakit":
            continue
        bundled = [
            s for s in catalog["skills"] if s.get("source_id") == sid or s.get("source") == sid
        ]
        if not bundled:
            continue
        key = repo_key(source.get("repository"))
        matches = by_repository.get(key, [])
        if not matches:
            item = {
                "id": "catalog-source-" + sid,
                "name": source["repository"],
                "repository": "https://github.com/" + key,
                "homepage": "https://github.com/" + key,
                "kind": "existing-pinned-skill-source",
                "domains": sorted(set().union(*(domains_for_skill(s) for s in bundled))),
                "attribution": source.get("attribution", source["repository"] + " contributors"),
                "paper_urls": [],
                "capabilities": [],
                "license_evidence": source.get("content_license", "see upstream"),
                "verified_urls": [
                    "https://github.com/" + key + "/tree/" + source.get("revision", "HEAD")
                ],
                "validation": "previously pinned; not re-executed in this survey",
                "research_part": "existing catalog",
            }
            projects[item["id"]] = item
            matches = [item]
            by_repository[key] = matches
        for item in matches:
            item["bundled_skills"] = sorted(
                set(item.get("bundled_skills", [])) | {s["id"] for s in bundled}
            )
            item.setdefault("pinned_source", source)
            item.setdefault("pinned_sources", []).append(source)
            item["domains"] = sorted(
                set(item["domains"]) | set().union(*(domains_for_skill(s) for s in bundled))
            )
    explicit_skill_matches: dict[str, list[str]] = {}
    if has_user_catalog:
        user_catalog.merge_evidence(ROOT, projects)
        user_catalog.validate_existing_functions(ROOT)
        explicit_skill_matches = user_catalog.related_skills(
            ROOT, catalog["skills"], projects, repo_key
        )
    analyses = read("src/dnakit/integrations/analyses.json")["analyses"]
    servers = read("src/dnakit/integrations/mcp_servers.json")["servers"]
    skill_ids = {s["id"] for s in catalog["skills"]}
    for analysis in analyses:
        if analysis.get("project_id") not in projects:
            raise ValueError(f"Analysis {analysis['id']} has no known source project.")
    for server in servers:
        if server.get("project_id"):
            known = server["project_id"] in projects
        else:
            known = any(
                repo_key(p.get("repository")) == repo_key(server["repository"])
                for p in projects.values()
            )
        if not known:
            raise ValueError(f"MCP server {server['id']} has no known source project.")
    for item in projects.values():
        key = repo_key(item.get("repository"))
        text = item["id"]
        native = ["dnakit.genomics." + f for f in NATIVE.get(text, [])]
        recipes = [a["id"] for a in analyses if a["project_id"] == text]
        matched_servers = [
            s
            for s in servers
            if (
                s["project_id"] == text
                if s.get("project_id")
                else key and repo_key(s["repository"]) == key
            )
        ]
        item["integration_sources"] = {
            "analyses": [
                {
                    field: copy.deepcopy(a[field])
                    for field in ("id", "sources", "citations", "license", "validation")
                    if field in a
                }
                for a in analyses
                if a["project_id"] == text
            ],
            "mcp_servers": [
                {
                    field: copy.deepcopy(s[field])
                    for field in (
                        "id",
                        "sources",
                        "license",
                        "attribution",
                        "source_commit",
                        "limitations",
                    )
                    if field in s
                }
                for s in matched_servers
            ],
        }
        item["paper_urls"] = list(
            dict.fromkeys(
                [
                    *item.get("paper_urls", []),
                    *(
                        url
                        for a in analyses
                        if a["project_id"] == text
                        for url in a.get("citations", [])
                    ),
                ]
            )
        )
        related = list(explicit_skill_matches.get(text, []))
        for match in item.get("local_matches", []):
            skill_id = match.get("skill_id") if isinstance(match, dict) else match
            if isinstance(skill_id, str) and skill_id in skill_ids:
                related.append(skill_id)
        if text == "rrblup":
            item["implementation_note"] = (
                "Independent fixed-lambda ridge implementation; does not call R rrBLUP or estimate REML."
            )
        item["dnakit_integration"] = {
            "native_functions": native,
            "existing_local_functions": item.get("existing_local_functions", []),
            "analysis_adapters": recipes,
            "mcp_servers": [s["id"] for s in matched_servers],
            "bundled_skills": item.get("bundled_skills", []),
            "related_local_skills": sorted(set(related)),
            "status": "available"
            if native
            or item.get("existing_local_functions")
            or recipes
            or item.get("bundled_skills")
            or matched_servers
            else "reference_only",
        }
        item["dnakit_integration"]["note"] = (
            "Availability refers to an interface or bundled workflow, not scientific execution or installed external dependencies."
        )
        if "gaps" in item:
            item["initial_gaps"] = item.pop("gaps")
        if "local_matches" in item:
            item["initial_local_matches"] = item.pop("local_matches")
        item["gap_audit_note"] = "Initial fields describe the repository before this expansion."
        if item["dnakit_integration"]["status"] == "reference_only":
            item["remaining_gaps"] = item.get("initial_gaps") or [
                "No project-specific DNAKit executable adapter or bundled workflow is registered."
            ]
        elif (
            native
            or item.get("existing_local_functions")
            or recipes
            or item["dnakit_integration"]["mcp_servers"]
        ):
            item["remaining_gaps"] = [
                "Only the listed entrypoints are integrated; this does not expose every upstream method.",
                "External services, references and optional packages require their own configured environment.",
            ]
        else:
            item["remaining_gaps"] = [
                "Bundled skill instructions/resources are available; project-specific execution is not established by this record."
            ]
    coverage = [
        {
            "id": key,
            "name_zh": name,
            "project_ids": sorted(p["id"] for p in projects.values() if key in p["domains"]),
        }
        for key, name in DOMAINS.items()
    ]
    return {
        "schema_version": 1,
        "snapshot_date": DATE,
        "scope": "All biological DNA directions; the initial 11 fields are not a boundary.",
        "completeness": "Dated, extensible survey, not proof that every DNA paper/project is captured. Existing pinned sources are distinguished from newly reviewed sources.",
        "queries": queries,
        "domains": coverage,
        "limitations": [
            *limitations,
            "No new biological dataset was submitted to external services during the survey.",
            "No claim of scientific equivalence from interface tests.",
            "BOLE benchmark evidence is public; a complete runnable framework or documented callable API was not verified."
            if has_user_catalog
            else "No public BOLE source or documented callable API was verified.",
        ],
        "projects": sorted(projects.values(), key=lambda p: p["id"]),
    }


def cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def faq(data: dict[str, Any], *, english: bool) -> str:
    title = "DNA project sources and coverage" if english else "DNA 项目来源与覆盖范围"
    intro = (
        "This dated survey separates original methods, upstream software, skills, MCP services and DNAKit adapters. Listing a project does not mean its software is installed or its results are validated."
        if english
        else "本页区分原始方法、上游软件、Skills、MCP 服务与 DNAKit 新增适配。收录项目不代表已安装依赖或验证其科学结果。"
    )
    lines = [
        f"# {title}",
        "",
        intro,
        "",
        f"Snapshot: **{DATE}**. Projects/sources: **{len(data['projects'])}**. Logged searches: **{len(data['queries'])}**.",
        "",
        (
            "The search extends beyond the initial 11 areas. It remains a versioned survey, not a claim to have enumerated every DNA paper or project."
            if english
            else "检索不限于最初 11 个领域；这是可持续补充的版本化清单，不能据此宣称已经穷尽全球全部 DNA 论文和项目。"
        ),
        "",
        "## 1) BOLE",
        "",
        (
            "Wang Zishuai, Liang Chongxiao, Zhang Yanlin, Zhou Rong, Zhang Xiaoai, Wei Wenkang and Li Kui describe a breeding agent framework. A public benchmark/evidence repository was found; a complete runnable framework or documented callable API was not verified. DNAKit records these sources without claiming a BOLE reimplementation. The local relationship and ridge-prediction functions independently implement the cited mathematical methods."
            if english
            else "Wang Zishuai、Liang Chongxiao、Zhang Yanlin、Zhou Rong、Zhang Xiaoai、Wei Wenkang、Li Kui 的 BOLE 描述了育种智能体框架。补充清单提供了公开基准与证据仓库，但尚未核实完整可运行框架或公开 API；因此保留来源，不声称复现 BOLE。本地关系矩阵与岭回归预测依据对应数学方法独立实现。"
        ),
        "",
        "[BOLE (2026)](https://doi.org/10.1016/j.jia.2026.07.011) · [VanRaden (2008)](https://doi.org/10.3168/jds.2007-0980) · [Endelman (2011)](https://doi.org/10.3835/plantgenome2011.08.0024)",
        "",
        "[BOLE benchmark evidence](https://github.com/ChongxiaoLiang/BOLE) · [335-project review / 补充清单审核](user_catalog.md)",
        "",
        "## 2) " + ("Coverage" if english else "方向覆盖"),
        "",
        "| Domain / 方向 | Project records / 项目记录数 |",
        "| --- | ---: |",
    ]
    lines += [
        f"| {cell(d['name_zh'])} (`{d['id']}`) | {len(d['project_ids'])} |" for d in data["domains"]
    ]
    lines += [
        "",
        "## 3) " + ("Who is credited by each integration?" if english else "每项功能引用了谁？"),
        "",
        (
            "Authors below are project/paper credits, not DNAKit authorship claims. Repository ownership is shown when individual authors were not verified. Original workflow citations remain in each skill's UPSTREAM.json and SKILL.md."
            if english
            else "以下列出项目或论文作者，不将上游工作归为 DNAKit 原创；未核实个人作者时标明仓库作者与贡献者。原有工作流引用仍保存在各 Skill 的 UPSTREAM.json 与 SKILL.md。"
        ),
        "",
        "| Project / 项目 | Authors and sources / 作者与来源 | DNAKit integration / 本地接口或状态 |",
        "| --- | --- | --- |",
    ]
    for p in data["projects"]:
        authors = (
            p.get("authors")
            or p.get("attribution")
            or (p.get("repository") or p["name"]) + " authors and contributors"
        )
        if isinstance(authors, list):
            authors = ", ".join(authors)
        source = p.get("homepage") or p.get("repository") or p["verified_urls"][0]
        papers = " ".join(
            f"[paper {i + 1}]({url})" for i, url in enumerate(p.get("paper_urls", []))
        )
        supplied = [
            url
            for url in p.get("user_supplied_publication_urls", [])
            if url not in p.get("paper_urls", [])
        ]
        if supplied:
            papers += " " + " ".join(
                f"[supplied source {i + 1}]({url})" for i, url in enumerate(supplied)
            )
        integration = p["dnakit_integration"]
        available = [
            *integration["native_functions"],
            *integration.get("existing_local_functions", []),
            *integration["analysis_adapters"],
            *["MCP: " + s for s in integration["mcp_servers"]],
        ]
        if integration["bundled_skills"]:
            available.append(f"{len(integration['bundled_skills'])} bundled Skills")
        if integration["related_local_skills"]:
            available.append(
                "Related local workflows: " + ", ".join(integration["related_local_skills"])
            )
        note = "; ".join(available) or (
            "Reference only; not integrated" if english else "来源记录；尚未集成"
        )
        if p.get("implementation_note"):
            note += "; " + (
                p["implementation_note"]
                if english
                else "独立实现固定惩罚的岭回归；不调用 R rrBLUP，也不估计 REML"
            )
        lines.append(f"| [{cell(p['name'])}]({source}) | {cell(authors)} {papers} | {cell(note)} |")
    lines += [
        "",
        "## 4) "
        + (
            "How to inspect evidence and dependencies"
            if english
            else "如何查看证据、依赖与具体差距"
        ),
        "",
        "```python",
        "from dnakit.integrations import search_projects, describe_analysis",
        "print(search_projects('BOLE'))",
        "print(search_projects(domain='dna-storage'))",
        "print(describe_analysis('plink2-pca'))",
        "```",
        "",
        (
            "The packaged search returns primary links, license evidence, local skill matches, gaps and actual adapter IDs. The source repository stores the search log in registry/genomics-survey.json and its research parts. Use dnakit skills cite SKILL_ID for pinned workflow attribution. External software/model/data terms apply separately."
            if english
            else "包内检索结果包含原始链接、许可证证据、已有技能匹配、差距与真实适配器 ID。仓库 registry/genomics-survey.json 及其分领域原始记录保存检索日志。使用 dnakit skills cite SKILL_ID 查看固定版本的工作流来源。外部软件、模型和数据适用其各自条款。"
        ),
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = build()
    payload = json.dumps(data, ensure_ascii=False, indent=2) + "\n"
    files = {
        ROOT / "registry/genomics-survey.json": payload,
        ROOT / "src/dnakit/integrations/projects.json": payload,
        ROOT / "docs/faq/genomics_sources.md": faq(data, english=False),
        ROOT / "docs/faq/genomics_sources.en.md": faq(data, english=True),
    }
    if (ROOT / "registry/user-catalog-335.json").exists():
        files[ROOT / "docs/faq/user_catalog.md"] = user_catalog.faq(ROOT, data, english=False)
        files[ROOT / "docs/faq/user_catalog.en.md"] = user_catalog.faq(ROOT, data, english=True)
        files[ROOT / "docs/api/features/39_external_genomics.md"] = user_catalog.external_guide(
            ROOT, english=False
        )
        files[ROOT / "docs/api/features/39_external_genomics.en.md"] = user_catalog.external_guide(
            ROOT, english=True
        )
    stale = []
    for path, content in files.items():
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    if stale:
        raise SystemExit("Stale survey files: " + ", ".join(stale))
    print(
        f"{'Checked' if args.check else 'Wrote'} {len(data['projects'])} project records and {len(data['queries'])} searches."
    )


if __name__ == "__main__":
    main()
