"""Verify or refresh the manuscript's complete DNAKit feature supplement.

The reviewed capability inventory supplies all counted rows. Supporting objects,
interfaces and development tools are retained separately without extra counts.
Only the named subsection under '# 补充材料' is replaced with --write.
"""

# Chinese punctuation matches the manuscript and reviewed capability inventory.
# ruff: noqa: RUF001

from __future__ import annotations

import argparse
import html
import re
import sys
from collections import Counter
from pathlib import Path

# Support both direct execution and package imports from validation tooling.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from others.validation.audit_table1 import (
    CATEGORIES,
    ROOT,
    audit_inventory,
    resolve_symbol,
)

INVENTORY = ROOT / "docs/planning/table1_dnakit_capabilities.csv"
DEFAULT_DOCUMENT = ROOT / "write_paper/手稿+提示词.md"
TITLE = "## DNAKit 网页功能清单"

# (display name, public API references, repository paths, scope description)
SUPPORT: tuple[tuple[str, tuple[str, ...], tuple[str, ...], str], ...] = (
    (
        "统一 DNA 输入与高级数据对象",
        ("dnakit.DNA", "dnakit.core.DNASequence", "dnakit.core.DNARecord", "dnakit.core.DNASet"),
        (),
        "DNA 接受单条、多条或 record mapping；高级对象表达序列、记录和集合。"
        "ID、描述、metadata、feature、拓扑和链类型是数据结构，不按字段追加功能数。",
    ),
    (
        "下标与切片",
        ("dnakit.DNA.__getitem__",),
        (),
        "下标/切片选择记录后仍返回 DNA；与内存记录选择同属数据访问，不追加计数。",
    ),
    (
        "特征、Gap、字母表与拓扑",
        (
            "dnakit.core.DNAFeature",
            "dnakit.Gap",
            "dnakit.core.DNAAlphabet",
            "dnakit.core.Topology",
            "dnakit.core.Strandedness",
        ),
        (),
        "表示功能区、已知/未知长度 Gap、字母表、线性/环状和单/双链声明；"
        "这些对象的使用已体现在编辑、校验、分析与可视化条目中。",
    ),
    (
        "坐标与复合区间",
        (
            "dnakit.core.Interval",
            "dnakit.core.CompoundLocation",
            "dnakit.core.UnresolvedLocation",
            "dnakit.core.import_location",
            "dnakit.core.export_location",
        ),
        (),
        "内部使用明确的坐标语义，支持复合与未解析区间；转换操作已在 A 列计数。",
    ),
    (
        "普通/流式读写与压缩",
        ("dnakit.read", "dnakit.write", "dnakit.io.iter_chunks"),
        (),
        "read 的 dna/stream 模式、流式写出、批次迭代和支持格式的 gzip 压缩"
        "归入已有读写/批处理能力，不按模式或扩展名重复计数。",
    ),
    (
        "热力学条件与光学修饰参数",
        (
            "dnakit.thermodynamics.ThermodynamicConditions",
            "dnakit.thermodynamics.OpticalModification",
            "dnakit.thermodynamics.LabelAbsorbanceCorrection",
        ),
        (),
        "盐、温度、浓度、修饰与吸光校正参数作用于相应计算；参数对象及不同参数组合不作为新算法。",
    ),
    (
        "图形高亮、样式、Gap 与自定义字符",
        (
            "dnakit.visualization.Highlight",
            "dnakit.visualization.SequencePlotConfig",
            "dnakit.visualization.SVGTheme",
            "dnakit.visualization.plot_sequence",
        ),
        (),
        "配置位置高亮、颜色、字体、布局、Gap 显示和自定义碱基显示字符；"
        "坐标与互补关系仍按真实序列计算，归入 H 列的绘图条目。",
    ),
    (
        "批量、并行、种子与进度",
        (
            "dnakit.batch.run_batch",
            "dnakit.batch.iter_batch",
            "dnakit.ops.evolution_generate",
            "dnakit.datasets.split",
        ),
        (),
        "在各入口支持范围内控制线程、seed、进度、取消和错误处理；"
        "同步/异步或串行/并行不与生物学功能数相乘。",
    ),
    (
        "实现来源、版本、缓存与产物追踪",
        (
            "dnakit.core.ImplementationInfo",
            "dnakit.core.Provenance",
            "dnakit.core.ArtifactRef",
            "dnakit.workflows.RunManifestBuilder",
            "dnakit.cache.JSONCache",
        ),
        (),
        "保存实现来源、参数、依赖、seed、产物与哈希等上下文；"
        "缓存的读写/失效等独立操作已在 I 列计数，构造器和记录字段不再追加。",
    ),
    (
        "结构化错误与警告",
        ("dnakit.exceptions.DNAKitError", "dnakit.core.Issue"),
        (),
        "用错误码、上下文、严重程度和提示表达失败或数据问题；不是独立分析算法。",
    ),
    (
        "Python 与命令行入口",
        ("dnakit.DNA", "dnakit.cli.run"),
        ("pyproject.toml", "docs/api/index.md"),
        "Python API 与 dnakit 命令行用于调用已有能力；CLI 只按实际命令覆盖范围使用，"
        "不宣称每个 Python API 都有单独的 CLI 命令。",
    ),
    (
        "Agent、MCP 与通用 LLM",
        (
            "dnakit.tools.create_server",
            "dnakit.tools.server.main",
            "dnakit.llm.call_llm",
            "dnakit.llm.acall_llm",
        ),
        ("docs/llm_tools/mcp.md", "docs/llm_tools/python_api.md"),
        "已有 Agent 目录/搜索/schema/执行、dnakit-mcp 服务和 LLM 调用能力均在 I 列；"
        "服务商、模型参数及被调用函数不重复累加。",
    ),
    (
        "Galaxy 包装入口",
        (),
        (
            "galaxy/dnakit/dnakit_convert.xml",
            "galaxy/dnakit/dnakit_deduplicate.xml",
            "galaxy/dnakit/dnakit_split.xml",
            "galaxy/dnakit/dnakit_report.xml",
        ),
        "仓库已有 convert、deduplicate、split、report 四个 Galaxy 包装器，"
        "分别复用转换、去重、划分与报告；不能沿用旧文中“不提供 Galaxy”的说法。",
    ),
    (
        "单元测试与一致性验证",
        (),
        ("others/tests", "others/validation/run_validation.py"),
        "项目测试入口为 pytest（测试目录 others/tests）；一致性报告由"
        " others/validation/run_validation.py 生成。这些开发检查不计为 DNA 运行时能力。",
    ),
    (
        "元数据分组与筛选",
        (
            "dnakit.datasets.split",
            "dnakit.datasets.joint_split",
            "dnakit.datasets.exclude_species",
            "dnakit.datasets.exclude_chromosomes",
            "dnakit.io.filter_by_metadata",
        ),
        (),
        "物种、染色体、个体等 metadata key 可用于已有分组/分层划分、联合约束和筛选入口；"
        "这些字段选项已归入 A/D 列，不因更换分组字段新增算法。",
    ),
    (
        "模型比较实验与文档构建",
        (),
        (
            "compared_models/使用方法.md",
            "compared_models/scripts/validate_compared_model_smoke.py",
            "mkdocs.yml",
        ),
        "模型比较实验按 compared_models 下的说明与脚本执行，文档由 MkDocs 构建。"
        "旧 benchmarks.benchmark_core.run_benchmarks 入口已不存在，不能继续列作可调用功能。",
    ),
)


def subsection_span(text: str) -> tuple[int, int]:
    """Locate the exact subsection without consuming adjacent manuscript content."""
    offset = 0
    in_supplement = False
    start: int | None = None
    fence: str | None = None
    for line in text.splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith(("```", "~~~")):
            marker = stripped[:3]
            fence = None if fence == marker else marker if fence is None else fence
        elif fence is None:
            heading = re.fullmatch(r"(#{1,6})\s+(.+?)\s*", stripped)
            if heading:
                level, title = len(heading[1]), heading[2]
                if start is not None and level <= 2:
                    return start, offset
                if level == 1:
                    in_supplement = title == "补充材料"
                if in_supplement and level == 2 and stripped == TITLE:
                    start = offset
        offset += len(line)
    if start is None:
        raise ValueError("Missing '# 补充材料' / '## DNAKit 网页功能清单'.")
    return start, len(text)


def cell(value: str) -> str:
    return html.escape(value, quote=False).replace("|", "&#124;").replace("\n", " ")


def api_cell(references: list[str] | tuple[str, ...]) -> str:
    return "<br>".join(f"`{reference}()`" for reference in references)


def module_label(row: dict[str, str]) -> str:
    labels = {
        "descriptor_field": "固定描述符字段",
        "dna_encoding": "传统 DNA 编码",
        "embedding_model": "DNA embedding 模型",
        "prediction_task": "官方预测任务",
        "ml_algorithm": "机器学习算法",
        "training_method": "自定义深度训练",
    }
    for key, label in labels.items():
        if row.get(key):
            return label
    if row.get("method_bindings"):
        return row["method_bindings"].split(":", 1)[0]
    return CATEGORIES[row["category"]]


def render(rows: list[dict[str, str]], counts: Counter[str], covered: int) -> str:
    """Render every counted entry once, plus verified, uncounted support notes."""
    total = sum(counts.values())
    catalog_counts = {
        key: sum(bool(row.get(key)) for row in rows)
        for key in (
            "descriptor_field",
            "dna_encoding",
            "embedding_model",
            "prediction_task",
            "ml_algorithm",
            "training_method",
        )
    }
    input_count = sum(len(row["task_bindings"].split(";")) for row in rows if row["task_bindings"])
    lines = [
        TITLE,
        "",
        "本清单按当前工作区（含未提交实现）核实，依据"
        "[逐项能力清单](../docs/planning/table1_dnakit_capabilities.csv)和"
        "[命名方法目录](../docs/planning/table1_method_catalog.json)，并对照"
        "[网页功能树](../docs/api/features/function_tree.md)及公开源码。"
        "“网页”指文档网站的功能展示，不代表每项都有独立的交互式 Web 表单。",
        "",
        f"补充表 S1–S9 共 **{total} 项计数能力**，与实验安排 Table 1 的九列对应。"
        "每行保留唯一清单 ID、功能名、全部登记公开入口及范围说明；ID 不连续不表示遗漏。"
        "每个描述符字段、登记模型/任务、算法或明确操作按已审核口径计数，"
        "别名、包装、参数组合及重复依赖不另加计。S10 保留数据对象、调用方式和工程支撑，"
        "这些说明不重复计入总数。",
        "",
        "| 补充表 | 对应列 | 能力数量 |",
        "| --- | --- | ---: |",
    ]
    for index, (category, label) in enumerate(CATEGORIES.items(), 1):
        lines.append(f"| S{index} | {category} {label} | {counts[category]} |")
    lines += [f"| 合计 | A–I | **{total}** |", ""]
    for index, (category, label) in enumerate(CATEGORIES.items(), 1):
        lines += [
            f"### 补充表 S{index}. {label}功能清单（{counts[category]} 项）",
            "",
            "| 清单 ID | 子模块 | 网页展示功能 | 调用函数名 | 函数功能与范围 |",
            "| --- | --- | --- | --- | --- |",
        ]
        for row in rows:
            if row["category"] == category and row["count"] == "1":
                fields = [
                    row["feature_id"],
                    cell(module_label(row)),
                    cell(row["capability"]),
                    api_cell(row["public_api"].split(";")),
                    cell(row["scope_note"]),
                ]
                lines.append("| " + " | ".join(fields) + " |")
        lines.append("")
    lines += [
        "### 补充表 S10. 数据对象、调用方式与工程支撑（不重复计数）",
        "",
        "| 编号 | 支撑内容 | 公开入口或项目文件 | 说明 |",
        "| --- | --- | --- | --- |",
    ]
    for index, (name, references, paths, description) in enumerate(SUPPORT, 1):
        for reference in references:
            resolve_symbol(reference)
        for path in paths:
            if not (ROOT / path).exists():
                raise ValueError(f"Supporting file does not exist: {path}")
        links = [f"[{cell(path)}](../{path})" for path in paths]
        entries = ([api_cell(references)] if references else []) + links
        lines.append(
            f"| P{index:02} | {cell(name)} | {'<br>'.join(entries)} | {cell(description)} |"
        )
    lines += [
        "",
        "### 完整性核对与使用边界",
        "",
        f"- S1–S9 的 **{total} 个唯一 ID** 与计数清单逐项一致，未省略描述符或模型明细。"
        f"公开导出和已审核业务类方法合并别名后覆盖 **{covered} 个实现**。",
        f"- 已核对 {catalog_counts['descriptor_field']} 个固定描述符字段、"
        f"{catalog_counts['dna_encoding']} 种 DNA 编码、"
        f"{catalog_counts['embedding_model']} 个 embedding 模型、"
        f"{catalog_counts['prediction_task']} 个预测任务及其 {input_count} 个输入绑定、"
        f"{catalog_counts['ml_algorithm']} 类机器学习算法、"
        f"{catalog_counts['training_method']} 种深度训练方法，"
        "以及命名方法目录中的逐项归属。科学等价性、全部私有实现与外部后端行为不能仅靠此检查证明。",
        "- 修正旧文的范围说明：变异下载支持 ClinVar VCF/XML/TSV 和 dbSNP VCF（A069）；"
        "提供方 MD5 文件的下载不等于自动校验该 MD5。已存在的 Galaxy 包装器列于 S10，"
        "已删除的旧 benchmark 入口不作为现有功能。",
        "- 可选 Python 依赖、外部程序、网络服务和模型权重按实现需要配置；"
        "固定描述符中的二核苷酸理化统计需要用户提供数值表，缺表会记录不可用原因。"
        "此清单核实接口、范围和计数，不表示所有外部模型均完成端到端实测。",
        "- [清单同步与核对脚本](../others/validation/sync_paper_features.py)会检查源码目录、"
        "逐项内容和计数；源代码变更提示必须人工复核，不能仅修改总数。",
        "",
        "在项目根目录复核本节：",
        "",
        "```bash",
        ".venv/bin/python others/validation/sync_paper_features.py",
        "```",
        "",
        "人工审核能力清单后，仅刷新本节：",
        "",
        "```bash",
        ".venv/bin/python others/validation/sync_paper_features.py --write",
        "```",
        "",
    ]
    return "\n".join(lines) + "\n"


def normalized(text: str) -> list[str]:
    """Ignore Markdown padding and separator widths, while checking all content."""
    result = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("|"):
            cells = [value.strip() for value in stripped.split("|")[1:-1]]
            if all(re.fullmatch(r"[-: ]+", value) for value in cells):
                continue
            stripped = "|".join(cells)
        result.append(stripped)
    return result


def check_section(actual: str, expected: str, rows: list[dict[str, str]]) -> None:
    identifiers = re.findall(r"^\|\s*([A-I]\d{3})\s*\|", actual, re.MULTILINE)
    desired = [row["feature_id"] for row in rows if row["count"] == "1"]
    if Counter(identifiers) != Counter(desired):
        missing = sorted(set(desired) - set(identifiers))
        extra = sorted(set(identifiers) - set(desired))
        duplicates = sorted(key for key, count in Counter(identifiers).items() if count > 1)
        raise ValueError(
            f"Feature IDs differ: missing={missing}, extra={extra}, duplicates={duplicates}"
        )
    if normalized(actual) != normalized(expected):
        raise ValueError(
            "Feature text/counts differ from reviewed content; review, then use --write."
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--document", type=Path, default=DEFAULT_DOCUMENT)
    parser.add_argument(
        "--write", action="store_true", help="Replace only the named supplement subsection."
    )
    args = parser.parse_args()
    try:
        rows, counts, covered = audit_inventory(INVENTORY)
        expected = render(rows, counts, covered)
        original = args.document.read_text(encoding="utf-8")
        start, end = subsection_span(original)
        if args.write and normalized(original[start:end]) != normalized(expected):
            updated = original[:start] + expected + original[end:]
            args.document.write_text(updated, encoding="utf-8")
        else:
            updated = original
        start, end = subsection_span(updated)
        check_section(updated[start:end], expected, rows)
        print(f"Verified {sum(counts.values())} unique capability rows across 9 tables.")
        print(
            f"Verified {len(SUPPORT)} uncounted support entries "
            f"and {covered} public implementations."
        )
        print("Source catalogs, IDs, API references, descriptions and table counts match.")
        return 0
    except (ValueError, KeyError, OSError, ImportError, AttributeError) as exc:
        print(f"Paper feature audit failed: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
