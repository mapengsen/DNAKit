"""CLI for genomic gene annotation with an optional reference assembly."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.exceptions import DNAKitError
from dnakit.io import ReadConfig, read
from dnakit.io._advanced_common import write_text_path


def annotate_genes_command(
    input_path: Annotated[
        Path,
        typer.Option(
            "--input",
            exists=True,
            dir_okay=False,
            help="Single-record genomic FASTA (optionally gzip).",
        ),
    ],
    species: Annotated[
        str, typer.Option("--species", help="Scientific name, alias, or NCBI TaxID.")
    ],
    assembly: Annotated[
        str | None,
        typer.Option(
            "--assembly",
            help="Optional accession/version; omitted selects and pins the species reference.",
        ),
    ] = None,
    output: Annotated[
        Path | None, typer.Option("--output", help="JSON output path; omit for stdout.")
    ] = None,
    reference_fasta: Annotated[
        Path | None, typer.Option("--reference-fasta", exists=True, dir_okay=False)
    ] = None,
    annotation_gff: Annotated[
        Path | None, typer.Option("--annotation-gff", exists=True, dir_okay=False)
    ] = None,
    cache_dir: Annotated[
        Path, typer.Option("--cache-dir", help="Project reference and index cache.")
    ] = Path(".dnakit/gene_annotation"),
    minimap2: Annotated[
        str, typer.Option("--minimap2", help="minimap2 executable name or path.")
    ] = "minimap2",
    preset: Annotated[
        str,
        typer.Option(
            "--preset", help="Genomic minimap2 preset; default assumes long accurate DNA."
        ),
    ] = "asm5",
    threads: Annotated[int, typer.Option("--threads", min=1, max=256)] = 1,
    min_mapq: Annotated[int, typer.Option("--min-mapq", min=1, max=254)] = 20,
    min_identity: Annotated[float, typer.Option("--min-identity", min=0.0, max=1.0)] = 0.90,
    min_alignment_length: Annotated[int, typer.Option("--min-alignment-length", min=1)] = 100,
    allow_download: Annotated[
        bool,
        typer.Option(
            "--allow-download/--no-download",
            help="Allow missing reference resources to be downloaded.",
        ),
    ] = True,
    api_key: Annotated[
        str | None,
        typer.Option(
            "--api-key",
            envvar="NCBI_API_KEY",
            help="Optional NCBI API key (or NCBI_API_KEY environment variable).",
        ),
    ] = None,
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Return named gene intervals, overlaps and intergenic/unknown segments."""

    try:
        if output is not None:
            inputs = [
                path.resolve()
                for path in (input_path, reference_fasta, annotation_gff)
                if path is not None
            ]
            if output.resolve() in inputs:
                raise ValueError("Output must not overwrite an input/reference file.")
            if output.exists() and not overwrite:
                raise FileExistsError(f"Refusing to overwrite {output}; use --overwrite.")
        with read(input_path, format="fasta", config=ReadConfig(max_records=2)) as source:
            records = iter(source)
            record = next(records, None)
            if record is None or next(records, None) is not None:
                raise ValueError("--input must contain exactly one genomic DNA record.")
        result = annotate_genes(
            record.sequence.symbols,
            sequence_id=record.id,
            species=species,
            assembly=assembly,
            reference_fasta=reference_fasta,
            annotation_gff=annotation_gff,
            config=GeneAnnotationConfig(
                cache_dir=str(cache_dir),
                minimap2_executable=minimap2,
                preset=preset,
                threads=threads,
                min_mapq=min_mapq,
                min_identity=min_identity,
                min_alignment_length=min_alignment_length,
                allow_download=allow_download,
                api_key=api_key,
            ),
            progress=lambda message: typer.echo(message, err=True),
        )
        payload = json.dumps(result.to_dict(), ensure_ascii=False, indent=2) + "\n"
        if output is None:
            typer.echo(payload, nl=False)
        else:
            write_text_path(
                output,
                lambda stream: stream.write(payload),
                overwrite=overwrite,
                create_parents=True,
            )
            typer.echo(f"Wrote {len(result.gene_records)} genes to {output}")
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        typer.echo(str(exc), err=True)
        if isinstance(exc, DNAKitError) and exc.context.get("candidates"):
            typer.echo(f"Candidate assemblies: {exc.context['candidates']}", err=True)
        raise typer.Exit(code=2) from exc
