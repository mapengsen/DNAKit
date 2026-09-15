"""Commands for finding and installing bundled DNA workflow skills."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Any

import typer
from rich.console import Console
from rich.table import Table

from dnakit.exceptions import DNAKitError
from dnakit.skills import (
    check_skill,
    cite_skill,
    get_skill,
    install_skill,
    list_skill_categories,
    list_skills,
    read_skill_file,
    search_skills,
    skill_catalog_stats,
)

app = typer.Typer(
    help="Search, inspect, cite, and install curated DNA skills.", no_args_is_help=True
)
console = Console()
error_console = Console(stderr=True)
JsonOption = Annotated[bool, typer.Option("--json", help="Return machine-readable JSON.")]


def _json(value: object) -> None:
    typer.echo(json.dumps(value, ensure_ascii=False, indent=2))


def _fail(exc: Exception) -> None:
    error_console.print(str(exc), markup=False)
    raise typer.Exit(2)


def _results(payload: dict[str, Any], as_json: bool) -> None:
    if as_json:
        _json(payload)
        return
    typer.echo(f"DNAKit skills: {payload['total']} matches; showing {len(payload['skills'])}")
    for item in payload["skills"]:
        typer.echo(f"\n{item['id']}")
        typer.echo(f"  {item['source_id']} | {item['license']} | {item['category']}")
        console.print(item["description"][:300], markup=False, highlight=False)
    if not payload["skills"]:
        typer.echo(
            "No matches. Try a shorter term, FASTQ, VCF, GWAS, or a category "
            "from 'dnakit skills categories'."
        )


@app.command("search")
def search_command(
    query: Annotated[
        str, typer.Argument(help="Task keywords, e.g. '变异注释', 'FASTQ QC', or 'GWAS'.")
    ],
    category: Annotated[str | None, typer.Option("--category")] = None,
    source: Annotated[str | None, typer.Option("--source")] = None,
    limit: Annotated[int, typer.Option("--limit", min=1, max=100)] = 10,
    as_json: JsonOption = False,
) -> None:
    """Search the local catalog; no network or model API is used."""
    try:
        _results(search_skills(query, category=category, source=source, limit=limit), as_json)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("list")
def list_command(
    category: Annotated[str | None, typer.Option("--category")] = None,
    source: Annotated[str | None, typer.Option("--source")] = None,
    offset: Annotated[int, typer.Option("--offset", min=0)] = 0,
    limit: Annotated[int, typer.Option("--limit", min=1, max=100)] = 20,
    as_json: JsonOption = False,
) -> None:
    """Browse a page of skills, optionally filtered by category or source."""
    try:
        _results(list_skills(category=category, source=source, offset=offset, limit=limit), as_json)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("categories")
def categories_command(as_json: JsonOption = False) -> None:
    """List searchable DNA task categories."""
    try:
        payload = list_skill_categories()
        if as_json:
            _json(payload)
            return
        table = Table("Category", "中文", "English", "Skills")
        for key, value in payload.items():
            table.add_row(key, value["zh"], value["en"], str(value["count"]))
        console.print(table)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("stats")
def stats_command() -> None:
    """Show catalog counts, licenses, and execution-validation status as JSON."""
    try:
        _json(skill_catalog_stats())
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("info")
def info_command(skill_id: str, as_json: JsonOption = False) -> None:
    """Inspect a skill's purpose, provenance, dependencies, and available files."""
    try:
        entry = get_skill(skill_id)
        if as_json:
            _json(entry)
            return
        typer.echo(entry["id"])
        typer.echo(entry.get("summary_zh") or entry["description"])
        typer.echo(f"Source: {entry['upstream_url']}")
        typer.echo(f"License: {entry['license']}")
        typer.echo(f"Execution validation: {entry['validation']['execution']}")
        typer.echo("Runtime requirements:")
        _json(entry.get("runtime", {}))
        typer.echo("Files:")
        for file in entry["files"]:
            typer.echo(f"  {file['path']}")
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("read")
def read_command(
    skill_id: str,
    resource: Annotated[
        str, typer.Option("--file", help="Relative bundled resource path.")
    ] = "SKILL.md",
    offset: Annotated[int, typer.Option("--offset", min=0)] = 0,
    limit: Annotated[int, typer.Option("--limit", min=1, max=100000)] = 24000,
    as_json: JsonOption = False,
) -> None:
    """Read a workflow or supporting text file, with bounded output."""
    try:
        result = read_skill_file(skill_id, resource=resource, offset=offset, limit=limit)
        if as_json:
            _json(result)
        else:
            typer.echo(result["content"], nl=False)
            if result["next_offset"] is not None:
                error_console.print(
                    f"\nMore content: use --offset {result['next_offset']}", markup=False
                )
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("doctor")
def doctor_command(skill_id: str) -> None:
    """Check declared local dependencies without executing a workflow."""
    try:
        result = check_skill(skill_id)
        _json(result)
        if result["status"] == "missing_dependencies":
            raise typer.Exit(1)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("cite")
def cite_command(
    skill_id: str,
    format: Annotated[str, typer.Option("--format", help="text, bibtex, or json.")] = "text",
) -> None:
    """Print a pinned upstream software citation for a paper or report."""
    try:
        result = cite_skill(skill_id, format=format)
        if isinstance(result, dict):
            _json(result)
        else:
            typer.echo(result)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("scripts")
def scripts_command(skill_id: str) -> None:
    """List bundled workflow scripts without requiring a source checkout."""
    from dnakit.workflows import scripts

    try:
        _json(scripts(skill_id))
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("run", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def run_command(
    ctx: typer.Context,
    skill_id: str,
    script: Annotated[
        str, typer.Option("--script", help="Bundled script shown by skills scripts.")
    ],
    workdir: Annotated[Path | None, typer.Option("--workdir")] = None,
) -> None:
    """Run a verified bundled script; pass its original arguments after --.

    Standard input/output and exit status pass through, including binary pipelines.
    Scientific software and model dependencies must be available in this environment.
    """
    import os
    import subprocess

    from dnakit.workflows.execution import _command

    try:
        _, command = _command(skill_id, script, list(ctx.args))
        result = subprocess.run(
            command, cwd=workdir, check=False, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
        )
        raise typer.Exit(result.returncode)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command(
    "support-run", context_settings={"allow_extra_args": True, "ignore_unknown_options": True}
)
def support_run_command(ctx: typer.Context, name: str) -> None:
    """Run a bundled support stage; pass its original options after --."""
    import os
    import subprocess

    from dnakit.workflows.execution import _support_command

    try:
        command = _support_command(name, list(ctx.args))
        result = subprocess.run(
            command, check=False, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
        )
        raise typer.Exit(result.returncode)
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("report")
def report_command(
    input_report: Annotated[Path, typer.Option("--input", help="WES report Markdown to render.")],
    output_dir: Annotated[Path, typer.Option("--output", "--output-dir")],
    language: Annotated[str, typer.Option("--language")] = "en",
    title: Annotated[str | None, typer.Option("--title")] = None,
    logo_left: Annotated[Path | None, typer.Option("--logo-left")] = None,
    logo_right: Annotated[Path | None, typer.Option("--logo-right")] = None,
    font_path: Annotated[Path | None, typer.Option("--font")] = None,
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Render the supplied report and tables as PDF, without generating new interpretations."""
    from dnakit.workflows import write_wes_report

    try:
        _json(
            write_wes_report(
                input_report,
                output_dir,
                language=language,
                title=title,
                logo_left=logo_left,
                logo_right=logo_right,
                font_path=font_path,
                overwrite=overwrite,
            )
        )
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)


@app.command("install")
def install_command(
    skill_id: str,
    agent: Annotated[
        str | None, typer.Option("--agent", help="codex, claude-code, or cursor.")
    ] = None,
    destination: Annotated[
        Path | None, typer.Option("--destination", help="Custom skill container directory.")
    ] = None,
    project: Annotated[
        Path | None,
        typer.Option("--project", help="Project root; defaults to the current directory."),
    ] = None,
    global_scope: Annotated[
        bool, typer.Option("--global", help="Install to the selected Agent's user directory.")
    ] = False,
    force: Annotated[
        bool,
        typer.Option(
            "--force", help="Replace a DNAKit-managed installation, including local edits."
        ),
    ] = False,
) -> None:
    """Install one complete skill; computational dependencies are configured separately."""
    try:
        _json(
            install_skill(
                skill_id,
                agent=agent,
                destination=destination,
                project_dir=project,
                global_scope=global_scope,
                force=force,
            )
        )
    except (DNAKitError, OSError, ValueError) as exc:
        _fail(exc)
