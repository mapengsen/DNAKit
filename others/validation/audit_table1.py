"""Audit the reviewed Table 1 task inventory and optionally update its DNAKit row.

Registered descriptors, encodings, embedding models, prediction tasks and ML
algorithms each count once. Remaining operations use reviewed scope boundaries.
Shared discovery, extraction and training wrappers do not add counts.
This checks source coverage and arithmetic, not scientific or backend execution.
"""

# Chinese punctuation intentionally matches the document's table and headings.
# ruff: noqa: RUF001

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib
import inspect
import json
import re
import sys
from collections import Counter
from functools import cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
METHOD_CATALOG = ROOT / "docs/planning/table1_method_catalog.json"

CATEGORIES = {
    "A": "序列获取与预处理",
    "B": "序列分析与分子设计",
    "C": "描述符、指纹与特征工程",
    "D": "数据集构建",
    "E": "DNA大模型",
    "F": "性质预测与建模（DNA注释）",
    "G": "序列集合评价",
    "H": "可视化",
    "I": "工作流、工程与Agent",
}


def resolve_symbol(name: str) -> Any:
    """Resolve an exported symbol, including a class method, without calling it."""
    parts = name.split(".")
    if parts[0] != "dnakit":
        raise ValueError(f"Only dnakit source references are accepted: {name}")
    # The full name is a symbol: read/normalize/deduplicate can also name modules.
    for length in range(len(parts) - 1, 0, -1):
        module_path = ROOT / "src" / Path(*parts[:length])
        if module_path.with_suffix(".py").is_file() or (module_path / "__init__.py").is_file():
            value: Any = importlib.import_module(".".join(parts[:length]))
            for part in parts[length:]:
                value = getattr(value, part)
            if isinstance(value, property):
                value = value.fget
            if not callable(value):
                raise ValueError(f"Reference is not an operation: {name}")
            return value
    raise ValueError(f"Source module was not found: {name}")


def symbol_key(value: Any) -> str:
    value = inspect.unwrap(getattr(value, "__func__", value))
    return f"{value.__module__}.{value.__qualname__}"


@cache
def reference_info(reference: str) -> tuple[str, str]:
    """Resolve each repeated API once during this read-only source audit."""
    value = resolve_symbol(reference)
    source_value = inspect.unwrap(getattr(value, "__func__", value))
    filename = inspect.getsourcefile(source_value)
    if filename is None:
        raise ValueError("No Python source found")
    relative = Path(filename).resolve().relative_to(ROOT)
    return symbol_key(value), f"{relative.as_posix()}:{inspect.getsourcelines(source_value)[1]}"


def public_functions() -> dict[str, Any]:
    """Discover public exports and directly defined bioinformatics functions."""
    modules = ["dnakit"]
    for path in sorted((ROOT / "src/dnakit").glob("*/__init__.py")):
        if path.parent.name != "cli":
            modules.append(f"dnakit.{path.parent.name}")
    modules.extend(("dnakit.batch", "dnakit.chunking"))
    functions: dict[str, Any] = {}
    for name in modules:
        module = importlib.import_module(name)
        for export in getattr(module, "__all__", ()):
            value = getattr(module, export)
            if inspect.isfunction(inspect.unwrap(value)):
                functions[f"{name}.{export}"] = value
    for path in sorted((ROOT / "src/dnakit/bioinformatics").glob("*.py")):
        if path.stem.startswith("_"):
            continue
        name = f"dnakit.bioinformatics.{path.stem}"
        module = importlib.import_module(name)
        for node in ast.parse(path.read_text(encoding="utf-8")).body:
            if isinstance(
                node, (ast.FunctionDef, ast.AsyncFunctionDef)
            ) and not node.name.startswith("_"):
                functions[f"{name}.{node.name}"] = getattr(module, node.name)
    catalog = json.loads(METHOD_CATALOG.read_text(encoding="utf-8"))
    for name in catalog["reviewed_classes"]:
        cls = resolve_symbol(name)
        for export, value in inspect.getmembers(cls):
            if export.startswith(("_", "to_", "from_")):
                continue
            if inspect.isfunction(value) or inspect.ismethod(value):
                functions[f"{name}.{export}"] = value
    return functions


def method_binding_errors(rows: list[dict[str, str]]) -> list[str]:
    """Check reviewed methods, their owners, and changes to reviewed source code.

    Source digests are review guards, not automatic proof of semantic coverage.
    --write deliberately does not refresh them: changed methods need re-review.
    """
    catalog = json.loads(METHOD_CATALOG.read_text(encoding="utf-8"))
    expected = {
        f"{group['name']}:{item['key']}": (group["category"], item["capability"])
        for group in catalog["groups"]
        for item in group["items"]
    }
    errors: list[str] = []
    seen: set[str] = set()
    for row in rows:
        bindings = list(filter(None, row.get("method_bindings", "").split(";")))
        if len(bindings) > 1:
            errors.append(f"Multiple methods counted as one entry: {row['feature_id']}")
        for binding in bindings:
            if binding in seen:
                errors.append(f"Duplicate method binding: {binding}")
            seen.add(binding)
            if binding not in expected:
                errors.append(f"Unknown method binding: {binding}")
            elif expected[binding][0] != row["category"] or row["count"] != "1":
                errors.append(f"Invalid method category/count: {row['feature_id']}: {binding}")
            elif expected[binding][1] != row["capability"]:
                errors.append(f"Method label differs from reviewed catalog: {binding}")
    if missing := expected.keys() - seen:
        errors.append(f"Unmapped method bindings: {sorted(missing)}")
    for reference, expected_digest in catalog["reviewed_sources"].items():
        value = inspect.unwrap(resolve_symbol(reference))
        digest = hashlib.sha256(inspect.getsource(value).encode("utf-8")).hexdigest()
        if digest != expected_digest:
            errors.append(f"Method source changed; review catalog and inventory: {reference}")
    for filename, expected_digest in catalog["reviewed_files"].items():
        digest = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        if digest != expected_digest:
            errors.append(f"Reviewed source file changed; re-review methods: {filename}")
    return errors


def registry_binding_errors(rows: list[dict[str, str]]) -> list[str]:
    """Check one counted entry per catalog item and task-specific input bindings."""
    from dnakit.bioinformatics.encoding import DNA_ENCODINGS
    from dnakit.descriptors import DESCRIPTOR_NAMES_V1
    from dnakit.predictions import available_machine_learning_models
    from dnakit.predictions.models import TASK_REGISTRY
    from dnakit.representations import available_embedding_models

    catalogs = {
        "descriptor_field": ("C", set(DESCRIPTOR_NAMES_V1)),
        "dna_encoding": ("C", set(DNA_ENCODINGS)),
        "embedding_model": ("E", set(available_embedding_models())),
        "prediction_task": ("F", {f"{model}:{task}" for model, task in TASK_REGISTRY}),
        "ml_algorithm": ("F", {model.name for model in available_machine_learning_models()}),
        "training_method": ("F", {"full", "linear_probe"}),
    }
    errors: list[str] = []
    for field, (category, registered) in catalogs.items():
        seen: set[str] = set()
        for row in rows:
            binding = row.get(field, "")
            if not binding:
                continue
            if binding not in registered:
                errors.append(f"Unknown {field}: {row['feature_id']}: {binding}")
            if binding in seen:
                errors.append(f"Duplicate {field}: {binding}")
            if row["category"] != category or row["count"] != "1":
                errors.append(f"Invalid category/count for {field}: {row['feature_id']}")
            seen.add(binding)
        if missing := registered - seen:
            errors.append(f"Unmapped {field}: {sorted(missing)}")
    for row in rows:
        if sum(bool(row.get(field)) for field in catalogs) > 1:
            errors.append(f"Multiple catalog items in one entry: {row['feature_id']}")
        if row["category"] == "E" and not row.get("embedding_model"):
            errors.append(f"Missing embedding model: {row['feature_id']}")
        task = row.get("prediction_task", "")
        actual = set(filter(None, row.get("task_bindings", "").split(";")))
        model_name, _, task_name = task.partition(":")
        spec = TASK_REGISTRY.get((model_name, task_name))
        expected = {f"{task}:{kind}" for kind in spec.input_kinds} if spec else set()
        if actual != expected:
            errors.append(f"Prediction input bindings differ from task: {row['feature_id']}")
    return errors


def audit_inventory(path: Path) -> tuple[list[dict[str, str]], Counter[str], int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Inventory is empty.")
    counts: Counter[str] = Counter()
    seen_ids: set[str] = set()
    seen_names: set[str] = set()
    covered: set[str] = set()
    bindings: set[str] = set()
    errors: list[str] = []
    owners: dict[str, set[str]] = {}
    for row in rows:
        identifier, category = row["feature_id"], row["category"]
        if identifier in seen_ids or not re.fullmatch(r"[A-IX]\d{3}", identifier):
            errors.append(f"Invalid or duplicate ID: {identifier}")
        if not identifier.startswith(category):
            errors.append(f"ID/category mismatch: {identifier}")
        if row["capability"] in seen_names or not row["scope_note"].strip():
            errors.append(f"Duplicate capability or missing scope: {identifier}")
        seen_ids.add(identifier)
        seen_names.add(row["capability"])
        if category not in {*CATEGORIES, "X"} or row["count"] not in {"0", "1"}:
            errors.append(f"Invalid category/count: {identifier}")
        elif (category == "X") != (row["count"] == "0"):
            errors.append(f"Only explicit exclusions may have count=0: {identifier}")
        else:
            counts[category] += int(row["count"])
        sources = set()
        for reference in row["public_api"].split(";"):
            try:
                key, source = reference_info(reference)
                covered.add(key)
                if row["count"] == "1":
                    owners.setdefault(key, set()).add(category)
                sources.add(source)
            except (AttributeError, ImportError, TypeError, ValueError) as exc:
                errors.append(f"{identifier}: {reference}: {exc}")
        row["source_refs"] = ";".join(sorted(sources))
        for binding in filter(None, row.get("task_bindings", "").split(";")):
            if binding in bindings:
                errors.append(f"Duplicate prediction binding: {binding}")
            bindings.add(binding)
    # Shared callables can serve distinct tasks or registered embedding models.
    # Their explicit task/model bindings determine the inventory entries.
    for key, categories in owners.items():
        if len(categories) > 1:
            errors.append(f"Callable assigned to multiple columns: {key}: {sorted(categories)}")
    discovered = public_functions()
    for name, value in discovered.items():
        if symbol_key(value) not in covered:
            errors.append(f"Public function has no task or exclusion: {name}")
    from dnakit.predictions.models import TASK_REGISTRY

    expected = {
        f"{model}:{task}:{kind}"
        for (model, task), spec in TASK_REGISTRY.items()
        for kind in spec.input_kinds
    }
    if bindings != expected:
        errors.append(f"Unmapped prediction bindings: {sorted(expected - bindings)}")
        errors.append(f"Unknown prediction bindings: {sorted(bindings - expected)}")
    errors.extend(registry_binding_errors(rows))
    errors.extend(method_binding_errors(rows))
    if errors:
        raise ValueError("\n".join(errors))
    return rows, counts, len({symbol_key(value) for value in discovered.values()})


def table_section(text: str) -> tuple[str, str, str]:
    start = text.index("# table1：")
    stop = text.index("## fig2：", start)
    return text[:start], text[start:stop], text[stop:]


def expected_row(counts: Counter[str]) -> str:
    cells = ["**DNAKit**", *(f"✓（{counts[key]}）" for key in CATEGORIES)]
    cells.append(f"（{sum(counts.values())}）")
    return "| " + " | ".join(cells) + " |"


def table_cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.split("|")[1:-1]]


def check_table(section: str, counts: Counter[str]) -> None:
    lines = [line for line in section.splitlines() if line.startswith("|")]
    if len(lines) < 3:
        raise ValueError("Table 1 must contain a header, separator, and DNAKit row.")
    header = table_cells(lines[0])
    expected_header = ["工具", *CATEGORIES.values(), "总共的功能数量"]
    if header != expected_header:
        raise ValueError("Table 1 columns differ from the documented task categories.")
    matches = [line for line in lines if line.split("|")[1].strip() == "**DNAKit**"]
    if len(matches) != 1 or table_cells(matches[0]) != table_cells(expected_row(counts)):
        raise ValueError("DNAKit row does not match the inventory; review, then use --write.")
    arithmetic = " + ".join(str(counts[key]) for key in CATEGORIES)
    if f"当前为 **{arithmetic} = {sum(counts.values())}**" not in section:
        raise ValueError("Review and update the total-count explanation to match the inventory.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inventory", type=Path, default=ROOT / "docs/planning/table1_dnakit_capabilities.csv"
    )
    parser.add_argument("--document", type=Path, default=ROOT / "实验安排.md")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Update source references and only the DNAKit row after a successful audit.",
    )
    args = parser.parse_args()
    try:
        rows, counts, function_count = audit_inventory(args.inventory)
        original = args.document.read_text(encoding="utf-8")
        before, section, after = table_section(original)
        if args.write:
            section, replacements = re.subn(
                r"^\|\s*\*\*DNAKit\*\*\s*\|[^\n]*$",
                lambda match: (
                    match[0]
                    if table_cells(match[0]) == table_cells(expected_row(counts))
                    else expected_row(counts)
                ),
                section,
                flags=re.MULTILINE,
            )
            if replacements != 1:
                raise ValueError("Expected exactly one DNAKit row in Table 1.")
        check_table(section, counts)
        if args.write:
            with args.inventory.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            updated = before + section + after
            if updated != original:
                args.document.write_text(updated, encoding="utf-8")
        for key, label in CATEGORIES.items():
            print(f"{key} {label}: {counts[key]}")
        print(
            f"Total: {' + '.join(str(counts[key]) for key in CATEGORIES)} = {sum(counts.values())}"
        )
        print(f"Public function/method implementations covered (aliases merged): {function_count}")
        print(
            "Source coverage and table arithmetic verified; "
            "external backends/models were not executed."
        )
        return 0
    except (ValueError, KeyError, OSError) as exc:
        print(f"Table 1 audit failed:\n{exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
