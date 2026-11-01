"""Local genomepy annotation operations and the fixed Boyle Lab blacklist source."""
from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

from datetime import datetime, timezone
from pathlib import Path
import importlib.metadata
import shutil
from typing import Any

from dnakit.bioinformatics._common import require
from dnakit.bioinformatics import _p40_helpers as _h


def _result(**values: Any) -> dict[str, Any]:
    return {"status": "completed", "backend": "genomepy", "version": importlib.metadata.version("genomepy"), **values}


def _annotation(path: str) -> Any:
    module = require("genomepy", "bio-genomepy")
    return module.Annotation(str(_h.existing(path, "annotation")), quiet=True)


def genomepy_list_installed(genomes_dir: str, *, verbose: bool = True) -> dict[str, Any]:
    """List genomes installed under an explicit local directory, with no provider lookup."""
    module = require("genomepy", "bio-genomepy")
    return _result(genomes=module.list_installed_genomes(genomes_dir=str(Path(genomes_dir).resolve())))


def genomepy_import_local_annotation(annotation: str, output_path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Copy a user-supplied local GTF/BED annotation to a new path and parse it with genomepy."""
    source = _h.existing(annotation, "annotation")
    target = _h.new_file(output_path, "output_path")
    obj = _annotation(str(source))
    kind = "bed" if ".bed" in source.suffixes else "gtf"
    frame = getattr(obj, kind)
    shutil.copyfile(source, target)
    return _result(outputs=[str(target)], record_count=len(frame), format=kind)


def genomepy_annotation_attributes(annotation: str, *, verbose: bool = True) -> dict[str, Any]:
    """List GTF attribute names using Annotation.attributes (upstream feature sampling)."""
    return _result(attributes=_annotation(annotation).attributes())


def genomepy_annotation_genes(annotation: str, *, file_format: str = "gtf", verbose: bool = True) -> dict[str, Any]:
    """List unique GTF gene_name values or BED names with Annotation.genes."""
    kind = _h.choice(file_format, ("gtf", "bed"), "file_format")
    return _result(genes=sorted(_annotation(annotation).genes(annot=kind)))


def genomepy_annotation_lengths(annotation: str, *, attribute: str = "gene_name", verbose: bool = True) -> dict[str, Any]:
    """Gene/transcript exon lengths with Annotation.lengths, preserving upstream coordinate arithmetic."""
    _h.choice(attribute, ("gene_name", "gene_id", "transcript_name", "transcript_id"), "attribute")
    lengths = _annotation(annotation).lengths(attribute)
    return _result(lengths={str(k): int(v) for k, v in lengths.items()})


def genomepy_annotation_dictionary(annotation: str, key: str, value: str, *, string_values: bool = True,
                                  verbose: bool = True) -> dict[str, Any]:
    """Make a key-value dictionary from GTF columns/attributes with Annotation.gtf_dict."""
    values = _annotation(annotation).gtf_dict(_h.text(key, "key"), _h.text(value, "value"), string_values=string_values)
    return _result(mapping={str(k): sorted((_h.plain(x) for x in v), key=str) if isinstance(v, (list, tuple, set)) else _h.plain(v)
                            for k, v in values.items()})


def genomepy_sanitize_annotation(annotation: str, genome: str, output_gtf: str, *, verbose: bool = True) -> dict[str, Any]:
    """Match GTF contigs against full FASTA headers, then drop unmatched rows using genomepy sanitize."""
    module = require("genomepy", "bio-genomepy")
    pd = require("pandas", "bio-analysis")
    source = _h.existing(genome, "genome")
    target = _h.new_file(output_gtf, "output_gtf")
    obj = _annotation(annotation)
    records = _h.read_fasta(source)
    obj.genome_file = str(source)
    obj.genome_contigs = [r[0] for r in records]
    # sanitize uses BED solely for the annotation contig names. Supplying this
    # minimal in-memory projection avoids a UCSC-tool conversion dependency.
    obj.bed = pd.DataFrame({"chrom": obj.gtf.seqname})
    before = len(obj.gtf)
    obj.sanitize(match=True, filter=True, overwrite=False)
    utils = require("genomepy.annotation.utils", "bio-genomepy")
    utils.write_annot(obj.gtf, str(target))
    return _result(outputs=[str(target)], input_rows=before, output_rows=len(obj.gtf))


def genomepy_download_blacklist(assembly: str, output_bed: str, *, verbose: bool = True) -> dict[str, Any]:
    """Download and decompress the Boyle Lab blacklist for an assembly supported by genomepy (network)."""
    import gzip
    import urllib.request
    plugin = require("genomepy.plugins.blacklist", "bio-genomepy")
    assembly = _h.choice(assembly, tuple(k for k in plugin.BlacklistPlugin.http_dict if k != "this was a triumph"), "assembly")
    output = _h.new_file(output_bed, "output_bed")
    url = plugin.BlacklistPlugin.http_dict[assembly]
    with urllib.request.urlopen(url, timeout=60) as response:
        content = gzip.decompress(response.read())
    output.write_bytes(content)
    return _result(outputs=[str(output)], source=url, retrieved_at=datetime.now(timezone.utc).isoformat())


__all__ = ["genomepy_list_installed", "genomepy_import_local_annotation", "genomepy_annotation_attributes",
           "genomepy_annotation_genes", "genomepy_annotation_lengths", "genomepy_annotation_dictionary",
           "genomepy_sanitize_annotation", "genomepy_download_blacklist"]

_register_public_verbose(globals(), tuple(__all__))
