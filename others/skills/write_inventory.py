"""Render the public inventory and review notes from the generated catalog."""
# Chinese documentation intentionally uses CJK punctuation.
# ruff: noqa: RUF001

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    decisions = json.loads((ROOT / "registry/decisions.json").read_text())["decisions"]
    entries = catalog["skills"]
    counts = Counter(entry["source_id"] for entry in entries)
    roles = Counter(entry.get("role", "discovery") for entry in entries)
    checked = sum(entry["validation"]["execution"] == "smoke_passed" for entry in entries)
    rows = [
        "# DNAKit Skills",
        "",
        "DNA／基因组学技能集合：统一搜索、按需安装、保留原作者及许可证。",
        "",
        f"当前快照：{catalog['snapshot_date']}。共 **{len(entries)}** 个技能："
        f"{roles['dna_workflow']} 个 DNA 相关工作流、"
        f"{roles['workflow_dependency']} 个组合流程必需的辅助技能、1 个 DNAKit 入口。",
        "",
        "每个子目录包含一个 `SKILL.md` 及所需辅助文件；原始来源和修改记录保留在 `UPSTREAM.md`。"
        "运行依赖以各技能说明为准。",
        "",
        "```bash",
        "# 在 DNAKit 源码根目录安装当前版本",
        'python -m pip install -e "."',
        'dnakit skills search "变异注释" --limit 5',
        "dnakit skills install dnakit --agent codex",
        "dnakit skills install variant-annotation --agent codex",
        "dnakit skills cite variant-annotation",
        "```",
        "",
        "安装入口后，在 Codex 中使用 `$dnakit` 描述任务。也支持 `--agent claude-code`、"
        "`--agent cursor`；默认安装到当前项目，可添加 `--global`。",
        "",
        "- [使用教程](../docs/llm_tools/skills.md)",
        "- [来源、筛选及复现](../registry/README.md)",
        "- [网络检索范围与补充记录](../registry/discovery.md)",
        "- [论文清单 CSV](../registry/skills.csv)",
        "- [第三方许可与署名](../THIRD_PARTY_NOTICES.md)",
        "",
        "## 来源",
        "",
        "| 来源 | 保留数量 |",
        "| --- | ---: |",
    ]
    rows.extend(f"| {source} | {count} |" for source, count in sorted(counts.items()))
    rows.extend(
        [
            "",
            "同一任务阶段的重复实现只保留一个，替代来源记录在目录中，旧名称可通过别名检索。"
            "输入条件不同的工作流、组合流程和组成步骤可分别保留。",
            "",
            f"元数据、文件完整性与执行验证分别记录。{checked} 个技能通过了小样本离线检查；"
            "其余技能未执行验证。完整流程成功率、Agent 自动选择质量和科学结果仍需专门评估。",
            "",
            "## 技能目录",
            "",
        ]
    )
    for category, labels in catalog["categories"].items():
        selected = [entry for entry in entries if entry["category"] == category]
        rows.extend(
            [
                f"### {labels['zh']} / {labels['en']} ({len(selected)})",
                "",
                "| Skill | 来源 | 用途 |",
                "| --- | --- | --- |",
            ]
        )
        for entry in selected:
            purpose = entry["description"].replace("\n", " ").replace("|", "\\|")
            if len(purpose) > 180:
                purpose = purpose[:177] + "..."
            rows.append(
                f"| [{entry['id']}]({entry['id']}/SKILL.md) | {entry['source_id']} | {purpose} |"
            )
        rows.append("")
    (ROOT / "skills/README.md").write_text("\n".join(rows), encoding="utf-8")

    notes = [
        "# DNAKit Skills 来源与复现记录",
        "",
        f"快照日期：{catalog['snapshot_date']}。Git commit 固定在 [sources.json](sources.json)。",
        "",
        f"本次检查 {len(decisions)} 个 `SKILL.md` 路径（包含聚合副本、示例及范围外内容）。"
        "这些路径数不是不同 DNA 功能的数量。",
        "",
        "| 筛选结果 | 路径数 |",
        "| --- | ---: |",
    ]
    notes.extend(f"| {status} | {count} |" for status, count in catalog["selection"].items())
    notes.extend(
        [
            "",
            f"最终分发 {len(entries)} 个技能，其中 1 个是 DNAKit 自有入口。",
            "",
            "## 选择与去重方法",
            "",
            "按任务、输入输出、分析阶段和必要条件判断是否重复。保留优先级及明确的等价组在 "
            "[curation_rules.py](../others/skills/curation_rules.py)。这一规则集可审查，但不是自动证明语义等价。",
            "第二轮逐来源、逐路径的补充判断见 [expansion.json](expansion.json)，"
            "搜索范围与排除理由见 [discovery.md](discovery.md)。",
            "",
            "例如：FASTQ 初始质控的重复指南保留一份；BWA 的 WGS/WES 与 Bowtie2 的 ChIP/ATAC "
            "输入场景不同，分别保留。WGS 完整流程与其变异注释步骤也分别保留。"
            "多组学与新抗原组合流程的必要辅助技能标记为 `workflow_dependency`。",
            "",
            "[decisions.json](decisions.json) 保留每个候选的来源、路径、筛选状态、替代项与理由。"
            "[catalog.json](catalog.json) 包含已收录内容、完整文件清单与 SHA-256。"
            "保留上游原文，必要的 YAML、技能标识、文档链接和本地环境路径兼容修复逐文件记录。"
            "原始与分发文件分别计算哈希；改动在 `UPSTREAM.md` 中说明。",
            "",
            "## 来源与许可证处理",
            "",
            "| 来源 | 固定版本 | 本版处理 |",
            "| --- | --- | --- |",
        ]
    )
    for source in catalog["sources"]:
        if source["id"] == "dnakit":
            continue
        url = f"https://github.com/{source['repository']}/tree/{source['revision']}"
        notes.append(
            f"| [{source['repository']}]({url}) | `{source['revision'][:12]}` "
            f"| {source['policy']} |"
        )
    notes.extend(
        [
            "",
            "[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) 说明许可范围。SciAgent 的 "
            "README 明确区分 CC-BY-4.0 原创说明内容和 frontmatter 中的底层工具许可；"
            "K-Dense 和 ClawBio 则明确要求遵守各技能自身的许可覆盖规则。",
            "",
            "MIT、Apache-2.0、CC-BY-4.0 和 AGPL-3.0-only 内容按原许可随包分发，不统一改成 MIT。"
            "ENCODE 独立技能保留完整可编辑源码；DeepMind 软件与文档分别适用 Apache 和 CC BY。"
            "暂缓项是对具体来源的判断，不表示 GPL 或 AGPL 一律禁止再分发。",
            "",
            "### 暂不复制的候选",
            "",
            "| 来源 | Skill | 理由 |",
            "| --- | --- | --- |",
        ]
    )
    notes.extend(
        f"| {item['source']} | {item['name']} | {item['reason']} |"
        for item in decisions
        if item["status"] in {"license_review", "quality_review"}
    )
    notes.extend(
        [
            "",
            "## 维护与复现命令",
            "",
            "先安装项目开发依赖。网络下载和中间结果都写入仓库忽略的 `temp/`。 "
            "`--fetch` 获取缺失的固定 Git 版本，不运行上游安装程序。",
            "",
            "```bash",
            'python -m pip install -e ".[dev]"',
            "python others/skills/build_catalog.py --fetch --plan",
            "# 审阅 temp/dnakit-skills-work/ 中的选择与决策后重新生成",
            "python others/skills/build_catalog.py --fetch",
            "python others/skills/build_catalog.py --check",
            "# 更新全部已审阅的 DNAKit 适配",
            "python others/skills/adaptations.py",
            "# --python 指向 Python 3.11+；运行已检查的离线脚本",
            "python others/skills/smoke_workflows.py --python python3",
            "python others/skills/write_inventory.py",
            "python -m build --outdir temp/dnakit-skills-dist",
            "```",
            "",
            "`skills/` 是公开分发目录；原生适配的维护模板在 `others/skills/adaptations/`，"
            "由 `registry/adaptations.json` 绑定上游版本。完整导入也会重放这些适配，"
            "保留旧名称别名、原始元数据、许可和哈希。全部技能都有维护模板；科学后端保留原始实现。"
            "构建时将已登记文件复制进 wheel。用户安装包后不需要 "
            "`temp` 中的原仓库，查询和安装也不需要 Git 或 Node.js。",
            "",
            "## 迁移清单",
            "",
            "[migration.json](migration.json) 逐项列出技能的接入方式、函数和脚本。",
            "[native-providers.json](native-providers.json) 固定数据库封装，",
            "[native-extensions.json](native-extensions.json) 记录 DNAKit 可选参数扩展。",
            "模板重放不运行科研工作流；执行验证须另行记录。",
            "",
            "## 验证范围与论文使用",
            "",
            "[validation.json](validation.json) 记录已安装技能的离线子进程检查：BED/GFF "
            "坐标转换、合成样本的 QC 标记和报告，以及 FASTA 去重、映射报告与异常输入处理。"
            "证据绑定文件清单哈希；文件变化后不继承旧结果。",
            "",
            "这些检查不包含大模型调用、外部数据库、全规模测序任务或完整流程成功率。"
            "元数据有效、文件哈希一致、依赖探测成功和科学结果正确是不同验证层级。",
            "",
            "论文可使用 [skills.csv](skills.csv) 和固定版本清单报告覆盖、来源、许可及去重结果。"
            "DNAKit 的贡献应表述为筛选、整合、检索、封装与工具接入；上游技能、软件、模型和"
            "数据库需分别署名引用。`dnakit skills cite ID --format bibtex` 提供固定版本的软件引用，"
            "不是虚构的论文 DOI。",
            "",
        ]
    )
    (ROOT / "registry/README.md").write_text("\n".join(notes), encoding="utf-8")


if __name__ == "__main__":
    main()
