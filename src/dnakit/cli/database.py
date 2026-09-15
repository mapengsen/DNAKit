"""Command-line access to bundled scientific database functions."""

from __future__ import annotations

import json

import typer

from dnakit import databases
from dnakit.exceptions import DNAKitError

app = typer.Typer(
    help="Search and call DNAKit scientific database functions.", no_args_is_help=True
)


@app.command("search")
def search(query: str, limit: int = 10) -> None:
    typer.echo(json.dumps(databases.search(query, limit=limit), ensure_ascii=False, indent=2))


@app.command("info")
def info(operation: str) -> None:
    typer.echo(json.dumps(databases.describe(operation), ensure_ascii=False, indent=2))


@app.command("call")
def call(operation: str, arguments: str = "{}") -> None:
    """Call an operation with a JSON object; retain original provider result fields."""
    try:
        result = databases.call(operation, json.loads(arguments))
        typer.echo(json.dumps(result, ensure_ascii=False, indent=2))
    except (DNAKitError, ValueError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(2) from exc
