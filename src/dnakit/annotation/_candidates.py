"""Annotation candidates: GFF3/GTF handling, gene structure and annotation lookups.

Implemented natively (no database): ``annotate`` (overlap of regions with a supplied annotation), ``assign_gene_ontology`` and
``assign_pathways`` (only from a mapping supplied in ``options``), ``build_gene_model``, ``check_annotation_consistency``,
``compare``, ``export``, ``extract``, ``index``, ``parse``, ``summarize`` and ``validate``.  ``assign_protein_domains`` and
``find_homologous_genes`` need external databases and return ``status: "not_implemented"``.  All coordinates are 1-based, closed.
"""

from __future__ import annotations

import itertools
import json
import re
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import (
    candidate_not_implemented,
    candidate_ok,
    run_candidate_operation,
)
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _impl_gff as _gff
from dnakit.exceptions import ConfigurationError, InputFormatError

_MODULE = "annotation"
_MAX_REGIONS = 100_000


def _opts(options: Mapping[str, object] | None) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping.", code="INVALID_OPTION")
    return dict(options or {})


def _option_bool(settings: Mapping[str, object], key: str, default: bool) -> bool:
    value = settings.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(f"options.{key} must be a boolean.", code="INVALID_OPTION")
    return value


def _overlaps(a: Mapping[str, object], start: int, end: int) -> bool:
    return int(a["start"]) <= end and int(a["end"]) >= start  # type: ignore[call-overload]


def _brief(feature: Mapping[str, object], fmt: str) -> dict[str, object]:
    return {
        "id": _gff.feature_id(feature, fmt),
        "name": _gff.attribute(feature, "Name", "gene_name", "transcript_name"),
        "type": feature["type"],
        "seqid": feature["seqid"],
        "start": feature["start"],
        "end": feature["end"],
        "strand": feature["strand"],
    }


def annotate(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Annotate genomic regions with the overlapping (and, if none overlap, the nearest) features of a supplied annotation.

    ``input``: regions as strings ``"chr1:100-200"`` / ``"chr1:150"`` (1-based, closed), mappings ``{"seqid"|"chrom", "start",
    "end", "id"}``, or a mapping with ``regions``.  Options: ``annotation`` (required: GFF3/GTF text, file path, or list of feature
    mappings), ``feature_types`` (default ``["gene"]``, case-insensitive), ``max_distance`` (nearest features farther than this
    many bp are not reported, default unlimited).  Per region: ``overlapping`` features (with ``overlap_bp``) and ``nearest``
    (``distance`` in bp and ``relation``: upstream/downstream relative to the feature strand, before/after when unstranded).
    Without ``annotation`` the result is ``not_implemented`` (annotating from a gene database is not supported).
    """

    settings = _opts(options)
    if settings.get("annotation") is None:
        return candidate_not_implemented(
            _MODULE,
            "annotate",
            input,
            options=options,
            requires="an annotation file (GFF3/GTF) passed as options['annotation']",
            message="No annotation source was supplied and DNAKit bundles no annotation database, so there is nothing to annotate against.",
        )
    fmt, features = _gff.load_features(settings["annotation"], settings)
    types = {str(item).lower() for item in settings.get("feature_types", ["gene"])}  # type: ignore[union-attr]
    max_distance = settings.get("max_distance")
    regions = _regions(input)
    by_seqid: dict[str, list[dict[str, object]]] = defaultdict(list)
    for feature in features:
        if str(feature["type"]).lower() in types:
            by_seqid[str(feature["seqid"])].append(feature)
    index: dict[str, tuple[list[dict[str, object]], list[int], list[int]]] = {}
    for seqid, items in by_seqid.items():
        items.sort(key=lambda f: (int(f["start"]), int(f["end"])))  # type: ignore[call-overload]
        starts = [int(f["start"]) for f in items]  # type: ignore[call-overload]
        prefix: list[int] = []  # index of the feature with the largest end among items[:i + 1]
        for position, item in enumerate(items):
            best = prefix[-1] if prefix else position
            prefix.append(position if int(item["end"]) > int(items[best]["end"]) else best)  # type: ignore[call-overload]
        index[seqid] = (items, starts, prefix)
    rows = []
    for region in regions:
        items, starts, prefix = index.get(region["seqid"], ([], [], []))  # type: ignore[call-overload]
        cut = bisect_right(starts, region["end"])  # type: ignore[arg-type]
        hits = []
        position = cut - 1
        while position >= 0 and int(items[prefix[position]]["end"]) >= region["start"]:  # type: ignore[operator, call-overload]
            if _overlaps(items[position], region["start"], region["end"]):  # type: ignore[arg-type]
                hits.append(position)
            position -= 1
        overlapping = []
        for position in sorted(hits):
            feature = items[position]
            size = (
                min(int(feature["end"]), region["end"])
                - max(int(feature["start"]), region["start"])
                + 1
            )  # type: ignore[call-overload, operator]
            overlapping.append({**_brief(feature, fmt), "overlap_bp": size})
        nearest = None
        if not overlapping:
            candidates = []
            if cut > 0:
                left = items[prefix[cut - 1]]
                candidates.append((region["start"] - int(left["end"]), left, "left"))  # type: ignore[operator, call-overload]
            if cut < len(items):
                right = items[cut]
                candidates.append((int(right["start"]) - region["end"], right, "right"))  # type: ignore[operator, call-overload]
            if candidates:
                distance, feature, side = min(candidates, key=lambda c: (c[0], c[2]))
                if max_distance is None or distance <= float(max_distance):  # type: ignore[arg-type]
                    strand = str(feature["strand"])
                    if strand in {"+", "-"}:
                        relation = (
                            ("upstream" if side == "right" else "downstream")
                            if strand == "+"
                            else ("downstream" if side == "right" else "upstream")
                        )
                    else:
                        relation = "after" if side == "left" else "before"
                    nearest = {**_brief(feature, fmt), "distance": distance, "relation": relation}
        rows.append({"region": region, "overlapping": overlapping, "nearest": nearest})
    return candidate_ok(
        _MODULE,
        "annotate",
        input,
        {"format": fmt, "feature_types": sorted(types), "annotations": rows},
        options=options,
    )


_REGION = re.compile(r"^\s*([^:\s]+):([\d,]+)(?:-([\d,]+))?\s*$")


def _regions(input: object) -> list[dict[str, object]]:
    items = input.get("regions") if isinstance(input, Mapping) and "regions" in input else input
    if isinstance(items, (str, Mapping)):
        items = [items]
    if not isinstance(items, Sequence):
        raise InputFormatError(
            "input must be a region string, a mapping or a list of them.", code="INVALID_REGION"
        )
    if len(items) > _MAX_REGIONS:
        raise InputFormatError(
            f"At most {_MAX_REGIONS} regions are supported.", code="TOO_MANY_REGIONS"
        )
    regions = []
    for position, item in enumerate(items):
        if isinstance(item, str):
            match = _REGION.match(item)
            if match is None:
                raise InputFormatError(
                    f"Region {item!r} is not 'seqid:start-end'.", code="INVALID_REGION"
                )
            seqid, start, end = (
                match.group(1),
                int(match.group(2).replace(",", "")),
                int((match.group(3) or match.group(2)).replace(",", "")),
            )
            identifier = item.strip()
        elif isinstance(item, Mapping):
            seqid = item.get("seqid", item.get("chrom", item.get("chromosome")))
            start = item.get("start")
            end = item.get("end", start)
            identifier = item.get("id", f"region{position + 1}")
        else:
            raise InputFormatError(
                f"Region {position} must be a string or mapping.", code="INVALID_REGION"
            )
        if (
            not isinstance(seqid, str)
            or any(isinstance(v, bool) or not isinstance(v, int) for v in (start, end))
            or start < 1
            or end < start
        ):
            raise InputFormatError(
                f"Region {item!r} needs a seqid and 1 <= start <= end.", code="INVALID_REGION"
            )
        regions.append({"id": identifier, "seqid": seqid, "start": start, "end": end})
    return regions


def _assign(
    operation: str, input: object, options: Mapping[str, object] | None, *, go: bool
) -> object:
    settings = _opts(options)
    mapping, sets, names = (
        settings.get("mapping"),
        settings.get("gene_sets"),
        settings.get("term_names") or {},
    )
    label = "Gene Ontology" if go else "pathway"
    if mapping is None and sets is None:
        return candidate_not_implemented(
            _MODULE,
            operation,
            input,
            options=options,
            requires=f"a {label} annotation supplied as options['mapping'] ({{gene: [terms]}}) or options['gene_sets'] ({{term: [genes]}})",
            message=f"DNAKit bundles no {label} database; {label} terms can only be assigned from a mapping you provide.",
        )
    genes_input = input.get("genes") if isinstance(input, Mapping) else input
    if isinstance(genes_input, str):
        genes = [item for item in re.split(r"[,\s]+", genes_input) if item]
    elif isinstance(genes_input, Sequence) and all(isinstance(item, str) for item in genes_input):
        genes = [item.strip() for item in genes_input if item.strip()]
    else:
        raise InputFormatError(
            "input must be a gene id, a comma-separated string or a list of gene ids.",
            code="INVALID_INPUT",
        )
    if not genes:
        raise InputFormatError("No gene identifiers were given.", code="INVALID_INPUT")
    insensitive = _option_bool(settings, "case_insensitive", False)
    key = (lambda text: text.casefold()) if insensitive else (lambda text: text)
    table: dict[str, list[str]] = defaultdict(list)
    if mapping is not None:
        if not isinstance(mapping, Mapping):
            raise ConfigurationError(
                "options.mapping must be a mapping gene -> list of terms.", code="INVALID_OPTION"
            )
        for gene, terms in mapping.items():
            for term in [terms] if isinstance(terms, str) else terms:  # type: ignore[union-attr]
                table[key(str(gene))].append(str(term))
    if sets is not None:
        if not isinstance(sets, Mapping):
            raise ConfigurationError(
                "options.gene_sets must be a mapping term -> list of genes.", code="INVALID_OPTION"
            )
        for term, members in sets.items():
            for gene in [members] if isinstance(members, str) else members:  # type: ignore[union-attr]
                table[key(str(gene))].append(str(term))
    if go:
        bad = sorted(
            {
                term
                for terms in table.values()
                for term in terms
                if not re.fullmatch(r"GO:\d{7}", term)
            }
        )
        if bad:
            raise InputFormatError(
                f"Not Gene Ontology identifiers (GO:#######): {bad[:5]}", code="INVALID_GO_TERM"
            )
    assignments = {
        gene: sorted(dict.fromkeys(table.get(key(gene), []))) for gene in dict.fromkeys(genes)
    }
    terms: dict[str, dict[str, object]] = {}
    for gene, found in assignments.items():
        for term in found:
            entry = terms.setdefault(
                term, {"name": names.get(term) if isinstance(names, Mapping) else None, "genes": []}
            )
            entry["genes"].append(gene)  # type: ignore[union-attr]
    unannotated = [gene for gene, found in assignments.items() if not found]
    result = {
        "assignments": assignments,
        "terms": dict(sorted(terms.items())),
        "unannotated": unannotated,
        "annotated_fraction": (len(assignments) - len(unannotated)) / len(assignments),
    }
    return candidate_ok(_MODULE, operation, input, result, options=options)


def assign_gene_ontology(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Look up Gene Ontology terms of genes in a gene-to-GO mapping that the caller supplies (no GO database is bundled).

    ``input``: a gene id, comma-separated string, list of gene ids, or a mapping with ``genes``.  Options: ``mapping``
    (``{gene: [GO:#######, ...]}``) and/or ``gene_sets`` (``{GO term: [genes]}``), ``term_names`` (``{term: name}``,
    optional), ``case_insensitive`` (default false).  Terms must look like ``GO:0008150`` (otherwise ``InputFormatError``).
    Result: ``assignments`` (gene -> sorted terms), ``terms`` (term -> name, genes), ``unannotated`` and
    ``annotated_fraction``.  Without ``mapping``/``gene_sets`` the result is ``not_implemented``.
    """

    return _assign("assign_gene_ontology", input, options, go=True)


def assign_pathways(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Look up pathways of genes in a gene-to-pathway mapping that the caller supplies (no pathway database is bundled).

    Same input and options as ``assign_gene_ontology`` but pathway identifiers are free-form strings (for example KEGG
    or Reactome ids).  Without ``mapping``/``gene_sets`` the result is ``not_implemented``.
    """

    return _assign("assign_pathways", input, options, go=False)


def assign_protein_domains(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: protein domains need profile-HMM or domain databases (Pfam, InterPro) and a search program."""

    return candidate_not_implemented(
        _MODULE,
        "assign_protein_domains",
        input,
        options=options,
        requires="Pfam/InterPro profile databases and hmmscan or InterProScan",
        message="Domain assignment needs a domain database and an HMM search; neither is bundled.",
    )


def build_gene_model(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Assemble gene models (gene -> transcripts -> exons, introns, CDS, UTRs) from GFF3/GTF text or feature mappings.

    ``input``: GFF3/GTF text, a file path (``.gz`` allowed), or a list of feature mappings.  Option ``format`` (auto, gff3,
    gtf).  Transcripts are mRNA/transcript/*RNA features (GFF3 ``Parent`` links, GTF ``transcript_id``).  Introns are the gaps
    between consecutive exons.  If a transcript has CDS (and, for GTF, start/stop codon) features but no annotated UTR, the
    5'/3' UTRs are derived as the exon parts outside the CDS span (``utr_derived`` is true; in GTF a stop codon outside CDS is
    taken into account).  Per transcript: ``exons``, ``introns``, ``cds``, ``five_prime_utr``, ``three_prime_utr`` (1-based
    closed intervals ``[start, end]`` in ascending genomic order), ``exon_count``, ``spliced_length``, ``cds_length`` and
    ``protein_length_aa`` (``cds_length / 3`` when divisible by 3).
    """

    settings = _opts(options)
    fmt, features = _gff.load_features(input, settings)
    models = _gff.gene_models(features, fmt)
    result = {
        "format": fmt,
        "gene_count": len([m for m in models if m["gene_id"]]),
        "transcript_count": sum(len(m["transcripts"]) for m in models),  # type: ignore[arg-type]
        "genes": models,
    }
    return candidate_ok(_MODULE, "build_gene_model", input, result, options=options)


def check_annotation_consistency(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Check the logical consistency of a gene annotation (children inside parents, exon/CDS geometry, phases).

    ``input``: GFF3/GTF text, file path or feature mappings.  Rules (code, severity): ``DUPLICATE_ID`` (GFF3, error), ``MISSING_PARENT``
    (Parent id not defined, error), ``CHILD_OUTSIDE_PARENT`` (transcript outside gene, exon/CDS/UTR outside transcript, error),
    ``SEQID_MISMATCH`` / ``STRAND_MISMATCH`` between child and parent (error), ``EXON_OVERLAP`` (exons of one transcript overlap,
    error), ``CDS_OUTSIDE_EXON`` (a CDS interval not contained in one exon, error), ``CDS_PHASE_MISMATCH`` (phase differs from the
    one implied by the preceding CDS lengths, error), ``CDS_LENGTH_NOT_MULTIPLE_OF_3`` (CDS plus stop codon length, warning),
    ``GENE_WITHOUT_TRANSCRIPT`` and ``TRANSCRIPT_WITHOUT_EXON`` (warning), ``BEYOND_SEQUENCE_END`` (``contig_lengths`` option
    ``{seqid: length}``, error).  Result: ``consistent`` (no error), ``issues`` and ``checked`` counts.
    """

    settings = _opts(options)
    fmt, features = _gff.load_features(input, settings)
    issues: list[dict[str, object]] = []

    def issue(
        severity: str, code: str, message: str, feature: Mapping[str, object] | None = None
    ) -> None:
        issues.append(
            {
                "severity": severity,
                "code": code,
                "line": None if feature is None else feature.get("line"),
                "message": message,
            }
        )

    by_id: dict[str, Mapping[str, object]] = {}
    for feature in features:
        identifier = _gff.feature_id(feature, fmt) if fmt == "gff3" else None
        if identifier is not None:
            if (
                identifier in by_id
                and str(by_id[identifier]["type"]).lower() != "cds"
                and str(feature["type"]).lower() != "cds"
            ):
                issue(
                    "error",
                    "DUPLICATE_ID",
                    f"ID {identifier!r} is used by more than one feature",
                    feature,
                )
            by_id.setdefault(identifier, feature)
    lengths = settings.get("contig_lengths")
    if isinstance(lengths, Mapping):
        for feature in features:
            limit = lengths.get(feature["seqid"])
            if limit is not None and int(feature["end"]) > int(limit):  # type: ignore[call-overload]
                issue(
                    "error",
                    "BEYOND_SEQUENCE_END",
                    f"{feature['type']} ends at {feature['end']} beyond {feature['seqid']} length {limit}",
                    feature,
                )

    def relate(child: Mapping[str, object], parent: Mapping[str, object]) -> None:
        if child["seqid"] != parent["seqid"]:
            issue(
                "error",
                "SEQID_MISMATCH",
                f"{child['type']} on {child['seqid']} but parent on {parent['seqid']}",
                child,
            )
        if int(child["start"]) < int(parent["start"]) or int(child["end"]) > int(parent["end"]):  # type: ignore[call-overload]
            issue(
                "error",
                "CHILD_OUTSIDE_PARENT",
                f"{child['type']} {child['start']}-{child['end']} lies outside its parent {parent['start']}-{parent['end']}",
                child,
            )
        if (
            child["strand"] in {"+", "-"}
            and parent["strand"] in {"+", "-"}
            and child["strand"] != parent["strand"]
        ):
            issue(
                "error",
                "STRAND_MISMATCH",
                f"{child['type']} strand {child['strand']} differs from parent strand {parent['strand']}",
                child,
            )

    if fmt == "gff3":
        for feature in features:
            for parent_id in _gff.parents(feature, fmt):
                if parent_id not in by_id:
                    issue(
                        "error",
                        "MISSING_PARENT",
                        f"Parent {parent_id!r} is not defined in the annotation",
                        feature,
                    )
                else:
                    relate(feature, by_id[parent_id])
    genes = _gff.build_models(features, fmt)
    transcripts = 0
    for gene_id, gene in genes.items():
        gene_feature = gene["feature"]
        if gene_id and gene_feature is not None and not gene["transcripts"]:
            issue(
                "warning",
                "GENE_WITHOUT_TRANSCRIPT",
                f"gene {gene_id!r} has no transcript",
                gene_feature,
            )  # type: ignore[arg-type]
        for entry in gene["transcripts"].values():  # type: ignore[union-attr]
            transcripts += 1
            transcript_feature = entry["feature"]
            if fmt == "gtf":
                if gene_feature is not None and transcript_feature is not None:
                    relate(transcript_feature, gene_feature)  # type: ignore[arg-type]
                if transcript_feature is not None:
                    for key in (
                        "exons",
                        "cds",
                        "five_prime_utr",
                        "three_prime_utr",
                        "utr",
                        "start_codon",
                        "stop_codon",
                    ):
                        for child in entry[key]:  # type: ignore[union-attr]
                            relate(child, transcript_feature)  # type: ignore[arg-type]
            exons = _gff.intervals(entry["exons"])  # type: ignore[arg-type]
            if (
                transcript_feature is not None
                and not exons
                and str(transcript_feature["type"]).lower() == "mrna"
            ):
                issue(
                    "warning",
                    "TRANSCRIPT_WITHOUT_EXON",
                    f"transcript {entry['id']!r} has no exon",
                    transcript_feature,
                )  # type: ignore[arg-type]
            for (a_start, a_end), (b_start, b_end) in itertools.pairwise(exons):
                if b_start <= a_end:
                    issue(
                        "error",
                        "EXON_OVERLAP",
                        f"transcript {entry['id']!r}: exons {a_start}-{a_end} and {b_start}-{b_end} overlap",
                        entry["exons"][0],
                    )  # type: ignore[index]
                    break
            cds_features = sorted(entry["cds"], key=lambda f: int(f["start"]))  # type: ignore[arg-type, call-overload]
            if exons:
                for cds in cds_features:
                    if not any(s <= int(cds["start"]) and int(cds["end"]) <= e for s, e in exons):  # type: ignore[call-overload]
                        issue(
                            "error",
                            "CDS_OUTSIDE_EXON",
                            f"transcript {entry['id']!r}: CDS {cds['start']}-{cds['end']} is not inside one exon",
                            cds,
                        )
            if cds_features:
                strand = _gff.transcript_strand(entry)
                ordered = cds_features[::-1] if strand == "-" else cds_features
                first_phase = ordered[0]["phase"]
                running = 0
                for cds in ordered:
                    if first_phase is not None and cds["phase"] is not None:
                        expected = (3 - ((running - int(first_phase)) % 3)) % 3  # type: ignore[call-overload]
                        if int(cds["phase"]) != expected:  # type: ignore[call-overload]
                            issue(
                                "error",
                                "CDS_PHASE_MISMATCH",
                                f"transcript {entry['id']!r}: CDS {cds['start']}-{cds['end']} has phase {cds['phase']}, expected {expected}",
                                cds,
                            )
                    running += _gff.length(cds)
                total = running + sum(_gff.length(c) for c in entry["stop_codon"])  # type: ignore[union-attr]
                if total % 3:
                    issue(
                        "warning",
                        "CDS_LENGTH_NOT_MULTIPLE_OF_3",
                        f"transcript {entry['id']!r}: CDS length {total} is not a multiple of 3 (partial CDS?)",
                        cds_features[0],
                    )
    errors = [item for item in issues if item["severity"] == "error"]
    result = {
        "consistent": not errors,
        "issues": issues,
        "error_count": len(errors),
        "warning_count": len(issues) - len(errors),
        "checked": {
            "format": fmt,
            "features": len(features),
            "genes": len([g for g in genes if g]),
            "transcripts": transcripts,
        },
    }
    return candidate_ok(_MODULE, "check_annotation_consistency", input, result, options=options)


def compare(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Compare two annotations by feature coordinates ``(seqid, type, start, end, strand)``.

    ``input``: the first annotation (GFF3/GTF text, path or feature mappings); the second is ``options["other"]``.  A pair
    ``[left, right]`` or a mapping ``{"left", "right"}`` is also accepted.  Options: ``feature_types`` (restrict to these
    types), ``max_examples`` (default 20, sample size of the differing features).  Result: ``left_count``, ``right_count``,
    ``shared``, ``only_left``, ``only_right``, ``jaccard`` (shared / union, 1.0 for two empty sets), ``by_type`` (per feature type
    counts) and the example lists ``only_left_examples`` / ``only_right_examples``.
    """

    settings = _opts(options)
    if isinstance(input, Mapping) and "left" in input and "right" in input:
        left_source, right_source = input["left"], input["right"]
    elif isinstance(input, (list, tuple)) and len(input) == 2 and not isinstance(input[0], Mapping):
        left_source, right_source = input
    else:
        left_source, right_source = input, settings.get("other")
        if right_source is None:
            raise ConfigurationError(
                "The second annotation is required.",
                code="OTHER_REQUIRED",
                hint="Pass options={'other': <annotation>} or a [left, right] pair.",
            )
    types = (
        {str(item).lower() for item in settings["feature_types"]}
        if "feature_types" in settings
        else None
    )  # type: ignore[union-attr]

    def keys(source: object) -> Counter[tuple[object, ...]]:
        _, features = _gff.load_features(source, settings)
        return Counter(
            (f["seqid"], str(f["type"]), f["start"], f["end"], f["strand"])
            for f in features
            if types is None or str(f["type"]).lower() in types
        )

    left, right = keys(left_source), keys(right_source)
    shared = left.keys() & right.keys()
    only_left, only_right = left.keys() - right.keys(), right.keys() - left.keys()
    union = len(left.keys() | right.keys())
    limit = int(settings.get("max_examples", 20))

    def example(key: tuple[object, ...]) -> dict[str, object]:
        return dict(zip(("seqid", "type", "start", "end", "strand"), key, strict=False))

    by_type: dict[str, dict[str, int]] = defaultdict(
        lambda: {"shared": 0, "only_left": 0, "only_right": 0}
    )
    for bucket, name in ((shared, "shared"), (only_left, "only_left"), (only_right, "only_right")):
        for key in bucket:
            by_type[str(key[1])][name] += 1
    result = {
        "left_count": len(left),
        "right_count": len(right),
        "shared": len(shared),
        "only_left": len(only_left),
        "only_right": len(only_right),
        "jaccard": len(shared) / union if union else 1.0,
        "by_type": dict(sorted(by_type.items())),
        "only_left_examples": [example(k) for k in sorted(only_left, key=str)[:limit]],
        "only_right_examples": [example(k) for k in sorted(only_right, key=str)[:limit]],
    }
    return candidate_ok(_MODULE, "compare", input, result, options=options)


def export(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Convert an annotation to GFF3, GTF, BED6, BED12 or JSON text (optionally writing a file).

    ``input``: GFF3/GTF text, path or feature mappings.  Options: ``format`` (target: "gff3" default, "gtf", "bed", "bed12",
    "json"), ``source_format`` (format of the input: auto, gff3, gtf), ``feature_type`` (BED6 only: the type to export, default "gene"), ``synthesize_parents`` (GTF to GFF3: create the
    missing gene/transcript lines from their children, default true), ``output`` (file path; nothing is written without it),
    ``overwrite`` (default false).  GFF3 to GTF resolves ``gene_id``/``transcript_id`` through ``Parent`` links.  BED uses 0-based
    half-open starts; BED12 has one line per transcript (blocks = exons, thick region = CDS plus start/stop codons, or empty when
    non-coding).  Result: ``format``, ``text``, ``line_count``, ``path``.
    """

    settings = _opts(options)
    target = str(settings.get("format", "gff3")).lower()
    if target not in {"gff3", "gtf", "bed", "bed12", "json"}:
        raise ConfigurationError(
            "options.format must be gff3, gtf, bed, bed12 or json.", code="INVALID_OPTION"
        )
    fmt, features = _gff.load_features(
        input, {**settings, "format": settings.get("source_format", "auto")}
    )
    if target == "gff3":
        text = _gff.to_gff3(
            features, fmt, synthesize_parents=_option_bool(settings, "synthesize_parents", True)
        )
    elif target == "gtf":
        text = _gff.to_gtf(features, fmt)
    elif target == "bed":
        text = _gff.bed6(features, fmt, str(settings.get("feature_type", "gene")))
    elif target == "bed12":
        text = _gff.bed12(features, fmt)
    else:
        text = json.dumps(features, indent=2, ensure_ascii=False) + "\n"
    path = settings.get("output")
    if path is not None:
        _gff.write_atomic(
            Path(str(path)), text, overwrite=_option_bool(settings, "overwrite", False)
        )
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    return candidate_ok(
        _MODULE,
        "export",
        input,
        {
            "format": target,
            "text": text,
            "line_count": len(lines),
            "path": None if path is None else str(path),
        },
        options=options,
    )


def _extract_features(
    entity: str, features: Sequence[Mapping[str, object]], fmt: str
) -> list[dict[str, object]]:
    if entity in {"intron", "utr", "five_prime_utr", "three_prime_utr"}:
        derived: list[dict[str, object]] = []
        for gene in _gff.build_models(features, fmt).values():
            for entry in gene["transcripts"].values():  # type: ignore[union-attr]
                model = _gff.transcript_model(entry)
                strand, seqid, tid = model["strand"], model["seqid"], model["transcript_id"]
                if entity == "intron":
                    pieces = [("intron", s, e) for s, e in model["introns"]]  # type: ignore[misc]
                    if strand == "-":
                        pieces.reverse()
                    for number, (kind, s, e) in enumerate(pieces, start=1):
                        derived.append(
                            _derived_feature(
                                kind, seqid, s, e, strand, tid, {"intron_number": number}
                            )
                        )
                    continue
                annotated = [*entry["five_prime_utr"], *entry["three_prime_utr"], *entry["utr"]]
                wanted = {
                    "utr": ("five_prime_utr", "three_prime_utr"),
                    "five_prime_utr": ("five_prime_utr",),
                    "three_prime_utr": ("three_prime_utr",),
                }[entity]
                if annotated and not model["utr_derived"]:
                    keep = {
                        "utr": annotated,
                        "five_prime_utr": entry["five_prime_utr"],
                        "three_prime_utr": entry["three_prime_utr"],
                    }[entity]
                    derived.extend({**dict(f), "derived": False} for f in keep)  # type: ignore[union-attr]
                else:
                    for kind in wanted:
                        for s, e in model[kind]:  # type: ignore[union-attr]
                            derived.append(_derived_feature(kind, seqid, s, e, strand, tid, {}))
        return derived
    matchers = {
        "gene": lambda f: str(f["type"]).lower() in _gff.GENE_TYPES,
        "transcript": lambda f: _gff.is_transcript_type(str(f["type"])),
        "mrna": lambda f: str(f["type"]).lower() == "mrna",
        "exon": lambda f: str(f["type"]).lower() == "exon",
        "cds": lambda f: str(f["type"]).lower() == "cds",
        "start_codon": lambda f: str(f["type"]).lower() == "start_codon",
        "stop_codon": lambda f: str(f["type"]).lower() == "stop_codon",
        "mirna": lambda f: (
            str(f["type"]).lower() in _gff.MIRNA_TYPES
            or (
                str(f["type"]).lower() in _gff.GENE_TYPES | {"transcript", "ncrna"}
                and "mirna" in _gff.biotypes(f)
            )
        ),
    }
    match = matchers.get(entity, lambda f: str(f["type"]).lower() == entity)
    return [{**dict(f), "derived": False} for f in features if match(f)]


def _derived_feature(
    kind: str,
    seqid: object,
    start: int,
    end: int,
    strand: object,
    transcript_id: object,
    extra: Mapping[str, object],
) -> dict[str, object]:
    return {
        "line": None,
        "seqid": seqid,
        "source": "dnakit",
        "type": kind,
        "start": start,
        "end": end,
        "score": None,
        "strand": strand,
        "phase": None,
        "attributes": {"Parent": transcript_id, **extra},
        "derived": True,
    }


def extract(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Extract features of one biological entity from GFF3/GTF text, or scan a DNA sequence for ORFs.

    Option ``entity`` (default "gene"): gene, transcript, mrna, exon, cds, start_codon, stop_codon, mirna (miRNA/pre_miRNA types or
    genes/transcripts with a miRNA biotype), intron (derived as the gaps between consecutive exons, numbered in transcript order),
    utr / five_prime_utr / three_prime_utr (annotated, or derived from exons and the CDS span, see ``build_gene_model``), any other
    value is matched against the feature type exactly (case-insensitive), and "orf".  For the annotation entities ``input`` is
    GFF3/GTF text, a path or feature mappings and the filters ``seqid``, ``strand``, ``start``/``end`` (overlap) and ``limit`` (default
    100000) apply; ``genome`` (FASTA text/path or ``{seqid: sequence}``) adds the genomic ``sequence`` of every feature (reverse
    complemented on the minus strand).  For ``orf``, ``input`` is a DNA sequence or FASTA (mapping ``{id: sequence}`` or list also ok)
    and the options are ``min_length`` (nucleotides including the stop codon, default 90), ``start_codons`` (default ["ATG"]),
    ``stop_codons`` (default TAA, TAG, TGA), ``require_stop`` (default true), ``both_strands`` (default true).  An ORF runs from the first
    start codon after a stop to the next in-frame stop codon (standard genetic code); coordinates are 1-based on the forward strand.
    """

    settings = _opts(options)
    entity = str(settings.get("entity", "gene")).lower()
    if entity == "orf":
        starts = [str(c).upper().replace("U", "T") for c in settings.get("start_codons", ["ATG"])]  # type: ignore[union-attr]
        stops = [
            str(c).upper().replace("U", "T")
            for c in settings.get("stop_codons", ["TAA", "TAG", "TGA"])
        ]  # type: ignore[union-attr]
        minimum = int(settings.get("min_length", 90))
        if minimum < 3 or any(len(c) != 3 for c in [*starts, *stops]):
            raise ConfigurationError(
                "min_length must be >= 3 and codons must have 3 letters.", code="INVALID_OPTION"
            )
        orfs = []
        for seqid, sequence in _gff.read_sequences(input):
            orfs += _gff.scan_orfs(
                seqid,
                sequence,
                min_length=minimum,
                start_codons=starts,
                stop_codons=stops,
                require_stop=_option_bool(settings, "require_stop", True),
                both_strands=_option_bool(settings, "both_strands", True),
            )
        return candidate_ok(
            _MODULE,
            "extract",
            input,
            {"entity": "orf", "count": len(orfs), "features": orfs},
            options=options,
        )
    fmt, features = _gff.load_features(input, settings)
    selected = _extract_features(entity, features, fmt)
    if "seqid" in settings:
        selected = [f for f in selected if f["seqid"] == settings["seqid"]]
    if "strand" in settings:
        selected = [f for f in selected if f["strand"] == settings["strand"]]
    if "start" in settings or "end" in settings:
        low, high = int(settings.get("start", 1)), int(settings.get("end", 2**62))  # type: ignore[call-overload]
        selected = [f for f in selected if _overlaps(f, low, high)]
    limit = int(settings.get("limit", 100_000))
    truncated = len(selected) > limit
    selected = selected[:limit]
    if settings.get("genome") is not None:
        genome = dict(_gff.read_sequences(settings["genome"]))
        from dnakit.sequence import reverse_complement_dna

        for feature in selected:
            sequence = genome.get(str(feature["seqid"]))
            if sequence is None or int(feature["end"]) > len(sequence):  # type: ignore[call-overload]
                raise InputFormatError(
                    f"The genome has no {feature['seqid']} covering {feature['start']}-{feature['end']}.",
                    code="GENOME_MISMATCH",
                )
            piece = sequence[int(feature["start"]) - 1 : int(feature["end"])]  # type: ignore[call-overload]
            feature["sequence"] = (
                reverse_complement_dna(piece) if feature["strand"] == "-" else piece
            )
    result = {
        "entity": entity,
        "format": fmt,
        "count": len(selected),
        "truncated": truncated,
        "features": selected,
    }
    return candidate_ok(_MODULE, "extract", input, result, options=options)


def find_homologous_genes(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: homology search needs orthology databases or a sequence-search program (BLAST/DIAMOND)."""

    return candidate_not_implemented(
        _MODULE,
        "find_homologous_genes",
        input,
        options=options,
        requires="an orthology database (Ensembl Compara, OrthoDB) or BLAST/DIAMOND with reference proteomes",
        message="Finding homologous genes needs reference data or a similarity search; neither is bundled.",
    )


def index(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Build lookup tables for a GFF3/GTF annotation: by id, by sequence (sorted intervals) and parent-to-children.

    ``input``: GFF3/GTF text, path or feature mappings.  Result: ``by_id`` (feature id -> seqid, type, start, end, strand),
    ``by_seqid`` (seqid -> list of ``[start, end, id_or_null, type]`` sorted by start then end), ``children`` (parent id -> ids of
    the features that list it, in file order; exon/CDS lines without an id are listed as ``"<type>@<line>"``),
    ``gene_transcripts`` (gene id -> transcript ids) and ``type_counts``.  Option ``format`` (auto, gff3, gtf).
    """

    settings = _opts(options)
    fmt, features = _gff.load_features(input, settings)
    by_id: dict[str, dict[str, object]] = {}
    by_seqid: dict[str, list[list[object]]] = defaultdict(list)
    children: dict[str, list[str]] = defaultdict(list)
    for feature in features:
        identifier = _gff.feature_id(feature, fmt)
        name = identifier or f"{feature['type']}@{feature['line']}"
        if identifier is not None and identifier not in by_id:
            by_id[identifier] = {
                "seqid": feature["seqid"],
                "type": feature["type"],
                "start": feature["start"],
                "end": feature["end"],
                "strand": feature["strand"],
            }
        by_seqid[str(feature["seqid"])].append(
            [feature["start"], feature["end"], identifier, feature["type"]]
        )
        for parent in _gff.parents(feature, fmt):
            children[parent].append(name)
    for rows in by_seqid.values():
        rows.sort(key=lambda row: (row[0], row[1]))  # type: ignore[arg-type, return-value]
    gene_transcripts = {
        gene_id: list(gene["transcripts"])
        for gene_id, gene in _gff.build_models(features, fmt).items()
        if gene_id
    }  # type: ignore[call-overload]
    result = {
        "format": fmt,
        "feature_count": len(features),
        "by_id": by_id,
        "by_seqid": dict(by_seqid),
        "children": dict(children),
        "gene_transcripts": gene_transcripts,
        "type_counts": dict(Counter(str(f["type"]) for f in features)),
    }
    return candidate_ok(_MODULE, "index", input, result, options=options)


def parse(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Parse GFF3 or GTF text (or a ``.gff``/``.gtf``/``.gz`` file) into feature records.

    ``input``: the annotation text or a file path.  Options: ``format`` (auto default: ``##gff-version`` or the attribute style
    decides; gff3, gtf), ``strict`` (default true: a malformed line raises ``InputFormatError``; false skips it and lists it under
    ``skipped``), ``limit`` (maximum features returned, default 100000).  GFF3 attributes are percent-decoded and comma lists become
    lists; GTF attributes are ``key "value";`` pairs (repeated keys become lists).  Each feature has ``line``, ``seqid``,
    ``source``, ``type``, ``start``, ``end`` (1-based closed), ``score``, ``strand``, ``phase`` and ``attributes``.  Result:
    ``format``, ``feature_count``, ``features``, ``directives`` (``##`` lines) and ``skipped``.
    """

    settings = _opts(options)
    strict = _option_bool(settings, "strict", True)
    if isinstance(input, Sequence) and not isinstance(input, (str, bytes)):
        fmt, features = _gff.load_features(input, settings)
        skipped: list[dict[str, object]] = []
        directives: list[str] = []
    else:
        text = _gff.read_text(input)
        fmt = str(settings.get("format", "auto"))
        fmt = _gff.detect_format(text) if fmt == "auto" else ("gff3" if fmt == "gff" else fmt)
        if fmt not in {"gff3", "gtf"}:
            raise ConfigurationError(
                "options.format must be auto, gff3 or gtf.", code="INVALID_OPTION"
            )
        features, skipped = _gff.parse_lines(text, fmt)
        errors = [item for item in skipped if item["severity"] == "error"]
        if strict and errors:
            raise InputFormatError(
                f"Malformed {fmt.upper()} at line {errors[0]['line']}: {errors[0]['message']}",
                code="INVALID_ANNOTATION",
                context={"line": errors[0]["line"], "error_count": len(errors)},
                hint="Pass strict=false to skip malformed lines.",
            )
        directives = [
            line
            for line in text.splitlines()
            if line.startswith("##") and not line.startswith("##FASTA")
        ]
    limit = int(settings.get("limit", 100_000))
    result = {
        "format": fmt,
        "feature_count": len(features),
        "truncated": len(features) > limit,
        "features": features[:limit],
        "directives": directives,
        "skipped": skipped,
    }
    return candidate_ok(_MODULE, "parse", input, result, options=options)


def summarize(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Summarize a GFF3/GTF annotation: feature counts and lengths per type, per sequence and strand, gene structure statistics.

    ``input``: GFF3/GTF text, path or feature mappings (option ``format``).  Result: ``feature_count``; ``types`` (type -> count,
    total_bp, mean_length); ``seqids``; ``strands``; ``gene_count``; ``transcript_count``; ``exon_count``; ``intron_count`` (derived
    from exons); ``coding_transcripts`` (with CDS); ``multi_exon_transcripts``; ``transcripts_per_gene`` and ``exons_per_transcript``
    (mean, max).  Lengths are ``end - start + 1``.
    """

    settings = _opts(options)
    fmt, features = _gff.load_features(input, settings)
    counts: Counter[str] = Counter()
    bases: Counter[str] = Counter()
    for feature in features:
        counts[str(feature["type"])] += 1
        bases[str(feature["type"])] += _gff.length(feature)
    models = _gff.gene_models(features, fmt)
    transcripts = [t for m in models for t in m["transcripts"]]  # type: ignore[union-attr]
    per_gene = [len(m["transcripts"]) for m in models if m["gene_id"]]  # type: ignore[arg-type]
    exon_counts = [int(t["exon_count"]) for t in transcripts]

    def stats(values: Sequence[int]) -> dict[str, object]:
        return {
            "mean": sum(values) / len(values) if values else 0.0,
            "max": max(values) if values else 0,
        }

    result = {
        "format": fmt,
        "feature_count": len(features),
        "types": {
            name: {
                "count": counts[name],
                "total_bp": bases[name],
                "mean_length": bases[name] / counts[name],
            }
            for name in sorted(counts)
        },
        "seqids": dict(sorted(Counter(str(f["seqid"]) for f in features).items())),
        "strands": dict(sorted(Counter(str(f["strand"]) for f in features).items())),
        "gene_count": len([m for m in models if m["gene_id"]]),
        "transcript_count": len(transcripts),
        "exon_count": sum(exon_counts),
        "intron_count": sum(len(t["introns"]) for t in transcripts),
        "coding_transcripts": sum(1 for t in transcripts if t["cds_length"]),
        "multi_exon_transcripts": sum(1 for n in exon_counts if n > 1),
        "transcripts_per_gene": stats(per_gene),
        "exons_per_transcript": stats(exon_counts),
    }
    return candidate_ok(_MODULE, "summarize", input, result, options=options)


def validate(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Validate the syntax of a GFF3 or GTF file (column count, coordinates, strand, phase, attributes, ids, parents).

    ``input``: annotation text or file path; option ``format`` (auto, gff3, gtf).  Errors: wrong column count, non-integer or
    non-positive coordinates, end < start, strand not in ``+ - . ?``, phase not in 0/1/2/'.', CDS without phase, unparsable attributes,
    (GFF3) duplicate ``ID`` of non-CDS features and ``Parent`` pointing to an undefined ``ID``, (GTF) missing ``gene_id`` (all lines) or
    ``transcript_id`` (all but gene lines).  Warning: GFF3 without ``##gff-version 3``.  Logical checks of gene structures are in
    ``check_annotation_consistency``.  Result: ``valid`` (no error), ``format``, ``feature_count``, ``errors``, ``warnings`` (each
    ``{line, severity, code, message}``).
    """

    settings = _opts(options)
    text = _gff.read_text(input)
    fmt = str(settings.get("format", "auto"))
    fmt = _gff.detect_format(text) if fmt == "auto" else ("gff3" if fmt == "gff" else fmt)
    if fmt not in {"gff3", "gtf"}:
        raise ConfigurationError("options.format must be auto, gff3 or gtf.", code="INVALID_OPTION")
    features, issues = _gff.parse_lines(text, fmt)
    if fmt == "gff3":
        if not any(line.startswith("##gff-version") for line in text.splitlines()):
            issues.append(
                {
                    "line": 1,
                    "severity": "warning",
                    "code": "MISSING_VERSION",
                    "message": "no '##gff-version 3' directive",
                }
            )
        seen: dict[str, int] = {}
        for feature in features:
            identifier = _gff.attribute(feature, "ID")
            if identifier is not None:
                if identifier in seen and str(feature["type"]).lower() != "cds":
                    issues.append(
                        {
                            "line": feature["line"],
                            "severity": "error",
                            "code": "DUPLICATE_ID",
                            "message": f"ID {identifier!r} already used at line {seen[identifier]}",
                        }
                    )
                seen.setdefault(identifier, int(feature["line"]))  # type: ignore[call-overload]
        for feature in features:
            for parent in _gff.as_list(feature["attributes"].get("Parent")):  # type: ignore[union-attr]
                if parent not in seen:
                    issues.append(
                        {
                            "line": feature["line"],
                            "severity": "error",
                            "code": "UNDEFINED_PARENT",
                            "message": f"Parent {parent!r} is not defined",
                        }
                    )
    else:
        for feature in features:
            if _gff.attribute(feature, "gene_id") is None:
                issues.append(
                    {
                        "line": feature["line"],
                        "severity": "error",
                        "code": "MISSING_GENE_ID",
                        "message": "GTF lines require a gene_id attribute",
                    }
                )
            if (
                str(feature["type"]).lower() not in _gff.GENE_TYPES
                and _gff.attribute(feature, "transcript_id") is None
            ):
                issues.append(
                    {
                        "line": feature["line"],
                        "severity": "error",
                        "code": "MISSING_TRANSCRIPT_ID",
                        "message": "GTF lines other than gene require a transcript_id attribute",
                    }
                )
    issues.sort(key=lambda item: (int(item["line"]), str(item["code"])))  # type: ignore[call-overload]
    errors = [item for item in issues if item["severity"] == "error"]
    result = {
        "valid": not errors,
        "format": fmt,
        "feature_count": len(features),
        "errors": errors,
        "warnings": [i for i in issues if i["severity"] != "error"],
    }
    return candidate_ok(_MODULE, "validate", input, result, options=options)


def convert_annotation_coordinates(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return run_candidate_operation(
        "annotation", "convert_annotation_coordinates", input, options=options
    )


def make_intergenic_regions(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return run_candidate_operation("annotation", "make_intergenic_regions", input, options=options)


def make_promoter_regions(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return run_candidate_operation("annotation", "make_promoter_regions", input, options=options)


def merge_gene_models(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return run_candidate_operation("annotation", "merge_gene_models", input, options=options)


__all__ = [
    "annotate",
    "assign_gene_ontology",
    "assign_pathways",
    "assign_protein_domains",
    "build_gene_model",
    "check_annotation_consistency",
    "compare",
    "convert_annotation_coordinates",
    "export",
    "extract",
    "find_homologous_genes",
    "index",
    "make_intergenic_regions",
    "make_promoter_regions",
    "merge_gene_models",
    "parse",
    "summarize",
    "validate",
]

_register_public_verbose(globals())
