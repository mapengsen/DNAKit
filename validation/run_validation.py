"""Generate an auditable local correctness-validation report for DNAKit."""

from __future__ import annotations

import argparse
import contextlib
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import tempfile
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final, Literal

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

import dnakit
from dnakit import DNARecord, DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise
from dnakit.datasets import (
    ClusterConfig,
    HierarchicalClusteringConfig,
    cluster_sequences,
    hierarchical_cluster,
)
from dnakit.descriptors import base_composition, gc_at_content, kmer_statistics, length_features
from dnakit.ops import reverse_complement
from dnakit.patterns import scan_restriction_sites
from dnakit.similarity import subsequence_search
from dnakit.thermodynamics import molecular_weight

SCHEMA_VERSION: Final = "dnakit.validation.v1"
ValidationStatus = Literal["pass", "fail", "not_comparable", "not_run"]


@dataclass(frozen=True, slots=True)
class CheckResult:
    id: str
    category: str
    status: ValidationStatus
    method: str
    observed: object
    expected: object
    tolerance: object
    difference: object
    note: str


def _version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _license(distribution: str) -> str | None:
    try:
        metadata = importlib.metadata.metadata(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None
    license_expression: str | None = metadata["License-Expression"]
    legacy_license: str | None = metadata["License"]
    return license_expression or legacy_license


def _manual_checks() -> list[CheckResult]:
    sequence = DNASequence("AACCGGTT")
    composition = base_composition(sequence)
    content = gc_at_content(sequence)
    kmers = kmer_statistics(DNASequence("AAAA"), 2)
    search = subsequence_search(DNASequence("AAA"), DNASequence("AAAA"))
    reverse = reverse_complement(DNASequence("AAGCT"))
    observed_composition = dict(composition.counts)
    observed_kmers = dict(kmers.counts)
    checks = [
        CheckResult(
            "MANUAL-001",
            "manual",
            "pass" if observed_composition == {"A": 2, "C": 2, "G": 2, "T": 2} else "fail",
            "hand-counted canonical composition",
            observed_composition,
            {"A": 2, "C": 2, "G": 2, "T": 2},
            "exact",
            None,
            "Eight bases contain two copies of each canonical symbol.",
        ),
        CheckResult(
            "MANUAL-002",
            "manual",
            "pass" if content.gc_fraction == 0.5 else "fail",
            "hand-counted GC fraction",
            content.gc_fraction,
            0.5,
            "exact rational result",
            None if content.gc_fraction is None else content.gc_fraction - 0.5,
            "Four of eight bases are G or C.",
        ),
        CheckResult(
            "MANUAL-003",
            "manual",
            "pass" if observed_kmers == {"AA": 3} else "fail",
            "overlapping exact k-mer count",
            observed_kmers,
            {"AA": 3},
            "exact",
            None,
            "AAAA contains AA at zero-based starts 0, 1, and 2.",
        ),
        CheckResult(
            "MANUAL-004",
            "manual",
            "pass" if [match.start for match in search.matches] == [0, 1] else "fail",
            "overlapping literal search",
            [match.start for match in search.matches],
            [0, 1],
            "exact",
            None,
            "AAA occurs twice in AAAA when overlaps are enabled.",
        ),
        CheckResult(
            "MANUAL-005",
            "manual",
            "pass" if reverse.symbols == "AGCTT" else "fail",
            "Watson-Crick reverse complement",
            reverse.symbols,
            "AGCTT",
            "exact",
            None,
            "Reverse AAGCT, then complement each canonical symbol.",
        ),
    ]
    return checks


def _boundary_checks() -> list[CheckResult]:
    empty = length_features(DNASequence(""))
    iupac = DNASequence("ACGTRYSWKMBDHVN", alphabet="iupac")
    long_sequence = DNASequence("ACGT" * 50_000)
    long_content = gc_at_content(long_sequence)
    circular_hits = scan_restriction_sites(DNASequence("AATTCG", topology="circular"), ("EcoRI",))
    circular_hit = circular_hits.hits[0] if len(circular_hits.hits) == 1 else None
    return [
        CheckResult(
            "BOUNDARY-001",
            "boundary",
            "pass" if empty.symbol_length == 0 and empty.coordinate_span == 0 else "fail",
            "empty linear sequence length",
            {"symbol_length": empty.symbol_length, "coordinate_span": empty.coordinate_span},
            {"symbol_length": 0, "coordinate_span": 0},
            "exact",
            None,
            "Empty linear DNA is a supported value object; circular empty DNA is rejected.",
        ),
        CheckResult(
            "BOUNDARY-002",
            "boundary",
            "pass" if iupac.symbol_length == 15 and iupac.ambiguity_count == 11 else "fail",
            "complete uppercase IUPAC alphabet construction",
            {"symbol_length": iupac.symbol_length, "ambiguity_count": iupac.ambiguity_count},
            {"symbol_length": 15, "ambiguity_count": 11},
            "exact",
            None,
            "Descriptor-specific ambiguity policies remain explicit at calculation time.",
        ),
        CheckResult(
            "BOUNDARY-003",
            "boundary",
            "pass"
            if long_sequence.symbol_length == 200_000 and long_content.gc_fraction == 0.5
            else "fail",
            "bounded 200000-nt descriptor input",
            {
                "symbol_length": long_sequence.symbol_length,
                "gc_fraction": long_content.gc_fraction,
            },
            {"symbol_length": 200_000, "gc_fraction": 0.5},
            "exact",
            None,
            "This is a boundary regression, not an unlimited-input or performance claim.",
        ),
        CheckResult(
            "BOUNDARY-004",
            "boundary",
            "pass"
            if circular_hit is not None
            and circular_hit.wraps_origin
            and circular_hit.top_cut == 0
            and circular_hit.bottom_cut == 4
            else "fail",
            "circular EcoRI site crossing sequence origin",
            None
            if circular_hit is None
            else {
                "wraps_origin": circular_hit.wraps_origin,
                "top_cut": circular_hit.top_cut,
                "bottom_cut": circular_hit.bottom_cut,
            },
            {"wraps_origin": True, "top_cut": 0, "bottom_cut": 4},
            "exact 0-based circular coordinates",
            None,
            "GAATTC is formed across the origin of circular AATTCG.",
        ),
    ]


def _biopython_checks() -> list[CheckResult]:
    try:
        bio_seq = importlib.import_module("Bio.Seq")
        bio_restriction = importlib.import_module("Bio.Restriction")
        bio_sequtils = importlib.import_module("Bio.SeqUtils")
        bio_align = importlib.import_module("Bio.Align")
    except ImportError:
        return [
            CheckResult(
                "BIOPYTHON-000",
                "external",
                "not_run",
                "optional dependency discovery",
                None,
                "Biopython installed",
                None,
                None,
                "Biopython is optional and was not installed in this environment.",
            )
        ]

    checks: list[CheckResult] = []
    restriction_sequence = "TTTGAATTCAAGGATCCCAAGCTTGGCCGCGGCCGCCCCGGG"
    enzyme_names = ("EcoRI", "BamHI", "HindIII", "HaeIII", "NotI", "SmaI")
    dnakit_hits = scan_restriction_sites(DNASequence(restriction_sequence), enzyme_names)
    dnakit_sites: dict[str, list[int]] = {name: [] for name in enzyme_names}
    for hit in dnakit_hits.hits:
        if hit.top_cut is not None:
            dnakit_sites[hit.enzyme].append(hit.top_cut)
    sequence_object = bio_seq.Seq(restriction_sequence)
    biopython_reported_sites = {
        name: list(getattr(bio_restriction, name).search(sequence_object)) for name in enzyme_names
    }
    # Bio.Restriction reports 1-based cleavage boundaries; DNAKit uses 0-based boundaries.
    biopython_sites = {
        name: [position - 1 for position in positions]
        for name, positions in biopython_reported_sites.items()
    }
    matches = all(dnakit_sites[name] == biopython_sites[name] for name in enzyme_names)
    checks.append(
        CheckResult(
            "BIOPYTHON-001",
            "restriction",
            "pass" if matches else "fail",
            "DNAKit 0-based top cleavage boundary versus Biopython search minus one",
            dnakit_sites,
            biopython_sites,
            "exact integer positions",
            None,
            (
                "Bio.Restriction.search reports 1-based cleavage boundaries, so one is "
                "subtracted before comparison. Recognition-site starts use a different "
                "convention and were not compared."
            ),
        )
    )

    molecular_sequences = tuple(
        ("ACGT" * ((length + 3) // 4))[:length] for length in (1, 2, 4, 10, 14, 30, 60)
    )
    molecular_rows: list[dict[str, object]] = []
    molecular_ok = True
    for text in molecular_sequences:
        dnakit_value = molecular_weight(DNASequence(text)).value_dalton
        biopython_value = bio_sequtils.molecular_weight(
            text, "DNA", double_stranded=False, circular=False
        )
        dnakit_phosphorylated = molecular_weight(
            DNASequence(text), five_prime_phosphorylated=True
        ).value_dalton
        raw_difference = dnakit_value - biopython_value
        phosphorylated_difference = dnakit_phosphorylated - biopython_value
        # DNAKit documents two-decimal residue masses and a 79.0-Da phosphate
        # approximation; Biopython uses higher-precision nucleotide masses.
        molecular_ok = molecular_ok and abs(raw_difference + 79.0) <= 1.0
        molecular_rows.append(
            {
                "sequence": text,
                "dnakit_unphosphorylated_dalton": dnakit_value,
                "dnakit_five_prime_phosphorylated_dalton": dnakit_phosphorylated,
                "biopython_five_prime_phosphate_dalton": biopython_value,
                "raw_difference_dalton": raw_difference,
                "phosphorylated_difference_dalton": phosphorylated_difference,
            }
        )
    checks.append(
        CheckResult(
            "BIOPYTHON-002",
            "molecular_weight",
            "pass" if molecular_ok else "fail",
            "average DNA mass with explicit terminal-convention reconciliation",
            molecular_rows,
            "Biopython 5-prime-phosphate convention and higher-precision residue table",
            "abs((DNAKit unphosphorylated - Biopython) + 79.0) <= 1.0 Da",
            None,
            (
                "Biopython computes sum of average unambiguous nucleotide weights minus "
                "(length-1)*18.0153 Da and states that nucleotide sequences carry a 5-prime "
                "phosphate. DNAKit uses rounded anhydrous residue masses, subtracts 61.96 Da "
                "for hydroxyl termini, and optionally adds an approximate 79.0 Da phosphate. "
                "Both unphosphorylated and phosphorylated errors are recorded across lengths; "
                "the one-Dalton tolerance validates the documented approximate model rather "
                "than claiming formula identity. A higher-precision revision should adopt one "
                "versioned mass table and explicitly test terminal formulas."
            ),
        )
    )

    query, target = "GATTACA", "GCATGCU".replace("U", "T")
    configuration = AlignmentConfig(
        mode="global", match_score=1.0, mismatch_score=-1.0, gap_score=-1.0
    )
    dnakit_alignment = align_pairwise(DNASequence(query), DNASequence(target), config=configuration)
    aligner = bio_align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -1.0
    biopython_score = float(aligner.score(query, target))
    score_difference = dnakit_alignment.score - biopython_score
    checks.append(
        CheckResult(
            "BIOPYTHON-003",
            "alignment",
            "pass" if abs(score_difference) <= 1e-12 else "fail",
            "global linear-gap score; match=1 mismatch=-1 gap=-1",
            dnakit_alignment.score,
            biopython_score,
            "absolute score difference <= 1e-12",
            score_difference,
            "Only optimal score is compared because equally optimal tracebacks may differ.",
        )
    )

    search_query = "AAA"
    search_target = "CAAAAT"
    dnakit_search_starts = [
        match.start
        for match in subsequence_search(
            DNASequence(search_query), DNASequence(search_target), overlapping=True
        ).matches
    ]
    biopython_search_starts = [
        position
        for position, match in bio_seq.Seq(search_target).search([search_query])
        if str(match) == search_query
    ]
    checks.append(
        CheckResult(
            "BIOPYTHON-004",
            "search",
            "pass" if dnakit_search_starts == biopython_search_starts else "fail",
            "overlapping literal subsequence search with zero-based starts",
            dnakit_search_starts,
            biopython_search_starts,
            "exact integer positions",
            None,
            "Both APIs perform literal overlapping search for this canonical query.",
        )
    )
    return checks


def _clustering_check() -> CheckResult:
    records = (
        DNARecord(DNASequence("AAAA"), "a"),
        DNARecord(DNASequence("AAAA"), "b"),
        DNARecord(DNASequence("TTTT"), "c"),
    )
    result = cluster_sequences(
        records,
        config=ClusterConfig(method="identity", threshold=1.0, max_records=10),
    )
    observed = [list(cluster.member_ids) for cluster in result.clusters]
    expected = [["a", "b"], ["c"]]
    return CheckResult(
        "ALGORITHM-001",
        "clustering",
        "pass" if observed == expected else "fail",
        "exhaustive exact-identity threshold graph connected components",
        observed,
        expected,
        "exact deterministic cluster membership and order",
        None,
        "This hand-check validates deterministic graph semantics; it is not a CD-HIT comparison.",
    )


def _biopython_clustering_checks() -> list[CheckResult]:
    try:
        bio_cluster = importlib.import_module("Bio.Cluster")
    except ImportError:
        return [
            CheckResult(
                "BIOCLUSTER-000",
                "external",
                "not_run",
                "optional Biopython clustering discovery",
                None,
                "Bio.Cluster installed",
                None,
                None,
                "Biopython is optional and Bio.Cluster was not available.",
            )
        ]

    records = tuple(
        DNARecord(DNASequence(symbols), record_id)
        for record_id, symbols in (
            ("a", "AAAA"),
            ("b", "AAAT"),
            ("c", "TTTT"),
            ("d", "TTTA"),
        )
    )
    # These are 1 - DNAKit global-alignment identity values. In particular,
    # AAAT versus TTTA has an optimal gapped identity of 0.2, not Hamming 0.5.
    lower_triangle = [[], [0.25], [1.0, 0.75], [0.75, 0.8, 0.25]]
    method_codes = {"single": "s", "complete": "m", "average": "a"}
    rows: list[dict[str, object]] = []
    passed = True
    for linkage, code in method_codes.items():
        dnakit_result = hierarchical_cluster(
            records,
            config=HierarchicalClusteringConfig(
                method="identity",
                linkage=linkage,  # type: ignore[arg-type]
            ),
        )
        reference = bio_cluster.treecluster(
            data=None,
            distancematrix=lower_triangle,
            method=code,
        )
        dnakit_distances = [step.distance for step in dnakit_result.linkage]
        reference_distances = [float(reference[index].distance) for index in range(len(reference))]
        match = all(
            abs(left - right) <= 1e-12
            for left, right in zip(dnakit_distances, reference_distances, strict=True)
        )
        passed = passed and match
        rows.append(
            {
                "linkage": linkage,
                "dnakit_merge_distances": dnakit_distances,
                "biopython_merge_distances": reference_distances,
                "matches": match,
            }
        )
    return [
        CheckResult(
            "BIOCLUSTER-001",
            "clustering",
            "pass" if passed else "fail",
            "hierarchical clustering on the same lower-triangular identity-distance matrix",
            rows,
            "Bio.Cluster single/complete/average linkage merge distances",
            "absolute merge-distance difference <= 1e-12",
            None,
            (
                "Node numbering and left/right order are implementation-specific, so the "
                "deterministic merge-distance sequence is compared for this no-tie fixture."
            ),
        )
    ]


def build_report(*, include_optional: bool = True, show_progress: bool = True) -> dict[str, Any]:
    """Run local checks. This function never imports, probes, installs, or calls NUPACK."""

    groups: list[tuple[str, Callable[[], list[CheckResult]]]] = [
        ("manual", _manual_checks),
        ("boundaries", _boundary_checks),
        ("clustering", lambda: [_clustering_check()]),
    ]
    if include_optional:
        groups.extend(
            (
                ("biopython", _biopython_checks),
                ("biopython-clustering", _biopython_clustering_checks),
            )
        )
    checks: list[CheckResult] = []
    progress = Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=Console(stderr=True),
        disable=not show_progress,
    )
    with progress:
        task = progress.add_task("DNAKit validation", total=len(groups))
        for name, runner in groups:
            progress.update(task, description=f"DNAKit validation: {name}")
            checks.extend(runner())
            progress.advance(task)
    counts = {
        status: sum(check.status == status for check in checks)
        for status in ("pass", "fail", "not_comparable", "not_run")
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "local_correctness_validation_not_paper_reproduction",
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "dnakit": dnakit.__version__,
            "biopython": _version("biopython"),
            "biopython_license": _license("biopython"),
            "primer3_cli": "not automatically discovered or executed",
        },
        "prohibited_backend_audit": {
            "primer3_automatic_discovery_attempted": False,
            "primer3_installation_attempted": False,
            "primer3_import_attempted": False,
            "primer3_call_attempted": False,
            "nupack_installation_attempted": False,
            "nupack_probe_attempted": False,
            "nupack_import_attempted": False,
            "nupack_call_attempted": False,
            "note": (
                "Primer3 and NUPACK are deliberately outside this runner. These booleans "
                "describe runner behavior, not a scan of the environment."
            ),
        },
        "summary": {"total": len(checks), **counts},
        "checks": [asdict(check) for check in checks],
        "interpretation": (
            "A passing check supports only its stated inputs, conventions, model, and tolerance. "
            "not_comparable means the reference was recorded but DNAKit has no equivalent "
            "implemented calculation. No paper-reproduction claim is made."
        ),
    }


def write_report(report: Mapping[str, Any], output: Path, *, force: bool = False) -> None:
    resolved = output.expanduser().resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{resolved.name}.", suffix=".tmp", dir=resolved.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(payload)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if force:
            os.replace(temporary_name, resolved)
        else:
            try:
                os.link(temporary_name, resolved)
            except FileExistsError as exc:
                raise FileExistsError(
                    f"Output already exists: {resolved}. Pass --force to replace it."
                ) from exc
            os.unlink(temporary_name)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(temporary_name)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation/results/local_validation_report.json"),
    )
    parser.add_argument("--skip-optional", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        report = build_report(include_optional=not arguments.skip_optional)
        write_report(report, arguments.output, force=arguments.force)
    except FileExistsError as exc:
        parser.error(str(exc))
    output = arguments.output.expanduser().resolve()
    Console(stderr=True).print(f"[green]Validation report written:[/green] {output}")
    return 1 if report["summary"]["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
