#!/usr/bin/env python3
"""Query the statgen-analysis OKG for datasets related to a trait.

Five query modes:
  --trait <label>      fuzzy match (ILIKE) on trait.attrs->>'label'
  --efo  <EFO:0000384> exact match on the EFO term node
  --seed-dataset <id>  follow dataset->trait edges from a seed dataset,
                       then collect sibling datasets sharing the same
                       trait OR the same EFO term
  --category <name>    filter on the curated trait.attrs->>'category'
                       (e.g. 'immune')

Filters:
  --ancestry EUR,EAS   normalised codes ({EUR, AFR, EAS, AMR, SAS,
                       admixed, unknown}); parses dataset.ancestry_scope
                       string into the controlled vocabulary
  --min-n N            (n_cases + n_controls) >= N  (rows without n_cases
                       are kept iff --min-n is not set; set --min-n 0 to
                       require some N attribute to exist)
  --sumstats-only      only datasets with a non-empty summary_statistics_url

Output:
  TSV to stdout with one row per matching dataset.
  Optional --json <path> writes the same payload as a JSON manifest.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional


# ---------------------------- Ancestry normalisation ----------------------------

# Maps a substring observed in dataset.attrs->>'ancestry_scope' to a
# controlled-vocab code. Match is greedy / order-independent: any free-text
# value can resolve to a set of codes (e.g. "European, African" -> {EUR, AFR}).
# When the eventual `ancestry_codes` derived attribute lands on the graph
# (OpenSpec `normalize-dataset-ancestry-codes`), the skill switches to a
# direct attrs->'ancestry_codes' ? <CODE> match and skips this parser.
ANCESTRY_PATTERNS = [
    (re.compile(r"european|caucasian|^eur\b", re.I),               "EUR"),
    (re.compile(r"\b(african|black)(?! american)|^afr\b", re.I),    "AFR"),
    (re.compile(r"african[ _-]?american", re.I),                    "AFR"),
    (re.compile(r"east[ _-]?asian|chinese|japanese|korean|^eas\b",
                 re.I),                                              "EAS"),
    (re.compile(r"south[ _-]?asian|indian|pakistani|bangladesh"
                r"|sri lankan|^sas\b", re.I),                        "SAS"),
    (re.compile(r"hispanic|latin[oa]|admixed[ _-]?american|^amr\b",
                 re.I),                                              "AMR"),
    (re.compile(r"\bmixed\b|admixed(?![ _-]?american)|multi[ _-]?"
                r"ancestry|trans[ _-]?ancestry", re.I),              "admixed"),
    (re.compile(r"not[ _-]?reported|unknown|^$", re.I),              "unknown"),
]


def ancestry_codes(scope: Optional[str]) -> set[str]:
    """Parse a free-text ancestry_scope into a set of controlled-vocab codes."""
    if not scope:
        return {"unknown"}
    codes: set[str] = set()
    for pat, code in ANCESTRY_PATTERNS:
        if pat.search(scope):
            codes.add(code)
    if not codes:
        codes.add("unknown")
    # If we matched anything other than 'unknown', drop 'unknown' from the set.
    if codes != {"unknown"} and "unknown" in codes:
        codes.discard("unknown")
    return codes


# ---------------------------- Postgres client ----------------------------

DEFAULT_DSN = "postgres://postgres:okg@localhost:5449/statgen_analysis"


def ensure_psycopg() -> None:
    """Install psycopg[binary] into the active env if missing. Idempotent."""
    try:
        import psycopg  # noqa: F401
        return
    except ImportError:
        pass
    # Install into the active env (writable for anaconda + venv setups).
    # Skip --user: that path is invisible to the anaconda interpreter even
    # though pip happily installs there.
    for cmd in [
        [sys.executable, "-m", "pip", "install", "psycopg[binary]"],
        ["uv", "pip", "install", "--python", sys.executable, "psycopg[binary]"],
    ]:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                print(f"[okg-search] installed psycopg via "
                      f"{' '.join(cmd[:2])}", file=sys.stderr)
                # Verify the install is importable from this Python.
                try:
                    import importlib
                    importlib.invalidate_caches()
                    import psycopg  # noqa: F401
                    return
                except ImportError:
                    print(f"[okg-search] install succeeded but psycopg is "
                          f"still not importable; trying next strategy",
                          file=sys.stderr)
                    continue
            else:
                print(f"[okg-search] {cmd[0]} pip install failed "
                      f"(rc={r.returncode}): {r.stderr[:300]}",
                      file=sys.stderr)
        except FileNotFoundError:
            continue
    sys.exit("ERROR: could not install psycopg; "
             f"`{sys.executable} -m pip install 'psycopg[binary]'` manually")


def connect(dsn: Optional[str] = None):
    import psycopg
    return psycopg.connect(dsn or os.environ.get("OKG_DSN", DEFAULT_DSN))


# ---------------------------- Trait / EFO resolution ----------------------------

def resolve_trait_ids_by_label(conn, label: str) -> list[str]:
    """Fuzzy match on trait.attrs->>'label' (ILIKE). Returns every matching
    trait node_id. Tries an exact match first; falls back to substring."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT node_id FROM okg.nodes_live
            WHERE subtype = 'trait' AND valid_to_gen IS NULL
              AND LOWER(attrs->>'label') = LOWER(%s)
        """, (label,))
        rows = [r[0] for r in cur.fetchall()]
        if rows:
            return rows
        cur.execute("""
            SELECT node_id FROM okg.nodes_live
            WHERE subtype = 'trait' AND valid_to_gen IS NULL
              AND LOWER(attrs->>'label') LIKE LOWER(%s)
        """, (f"%{label}%",))
        return [r[0] for r in cur.fetchall()]


def resolve_trait_ids_by_efo(conn, efo: str) -> list[str]:
    """Find traits mapped_to an EFO term. Accepts `EFO:0000384`,
    `term:efo:0000384`, or the bare `0000384`."""
    efo_clean = efo.strip()
    if efo_clean.lower().startswith("term:efo:"):
        term_id = efo_clean.lower()
    elif efo_clean.upper().startswith("EFO:"):
        term_id = "term:efo:" + efo_clean.split(":", 1)[1].strip().lstrip("0").rjust(7, "0")
    else:
        term_id = "term:efo:" + efo_clean.strip().lstrip("0").rjust(7, "0")
    with conn.cursor() as cur:
        cur.execute("""
            SELECT src FROM okg.edges_live
            WHERE valid_to_gen IS NULL
              AND edge_type = 'mapped_to' AND dst = %s
              AND src LIKE 'trait:%%'
        """, (term_id,))
        return [r[0] for r in cur.fetchall()]


def resolve_trait_ids_by_seed_dataset(conn, dataset_id: str) -> list[str]:
    """Follow dataset --references--> trait edges; return all trait IDs."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT dst FROM okg.edges_live
            WHERE valid_to_gen IS NULL
              AND edge_type = 'references'
              AND src = %s AND dst LIKE 'trait:%%'
        """, (dataset_id,))
        return sorted({r[0] for r in cur.fetchall()})


def expand_trait_ids_via_efo(conn, trait_ids: list[str]) -> list[str]:
    """Given a seed set of trait IDs, walk to their EFO terms and pull in
    every sibling trait that maps to the same terms. Catches both the
    curated `trait:*` namespace and the scale-ingest
    `trait:gwas_catalog_scale:*` namespace whose members share an EFO
    target. No-op when the input set has no EFO mappings."""
    if not trait_ids:
        return []
    expanded = set(trait_ids)
    with conn.cursor() as cur:
        cur.execute("""
            SELECT DISTINCT dst FROM okg.edges_live
            WHERE valid_to_gen IS NULL
              AND edge_type = 'mapped_to'
              AND src = ANY(%s) AND dst LIKE 'term:efo:%%'
        """, (list(expanded),))
        efo_terms = [r[0] for r in cur.fetchall()]
        if not efo_terms:
            return sorted(expanded)
        cur.execute("""
            SELECT DISTINCT src FROM okg.edges_live
            WHERE valid_to_gen IS NULL
              AND edge_type = 'mapped_to'
              AND dst = ANY(%s)
              AND src LIKE 'trait:%%'
        """, (efo_terms,))
        for r in cur.fetchall():
            expanded.add(r[0])
    return sorted(expanded)


def resolve_trait_ids_by_category(conn, category: str) -> list[str]:
    """Filter curated traits by attrs->>'category'. Only the canonical
    `trait:*` namespace has category populated (not `trait:gwas_catalog_scale:*`)."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT node_id FROM okg.nodes_live
            WHERE subtype = 'trait' AND valid_to_gen IS NULL
              AND LOWER(attrs->>'category') = LOWER(%s)
        """, (category,))
        return [r[0] for r in cur.fetchall()]


# ---------------------------- Dataset fetch + filter ----------------------------

def fetch_datasets_referencing_traits(conn, trait_ids: list[str]) -> list[dict]:
    """Return one row per dataset_metadata node that references any trait in
    `trait_ids`. Joins the trait label so the caller can show what it matched.
    """
    if not trait_ids:
        return []
    with conn.cursor() as cur:
        cur.execute("""
            SELECT
              d.node_id,
              d.attrs->>'accession' AS gcst,
              d.attrs->>'disease_trait' AS trait,
              COALESCE(d.attrs->>'ancestry_scope', '') AS ancestry_scope,
              COALESCE((d.attrs->>'n_cases')::bigint, 0)
                + COALESCE((d.attrs->>'n_controls')::bigint, 0) AS n_total,
              COALESCE(d.attrs->>'n_cases', '') AS n_cases,
              COALESCE(d.attrs->>'n_controls', '') AS n_controls,
              COALESCE(d.attrs->>'summary_statistics_url', '') AS sumstats_url,
              COALESCE(d.attrs->>'release_date', '') AS release_date,
              COALESCE(d.attrs->>'genome_build', '') AS build,
              COALESCE(d.attrs->>'first_author', '') AS author,
              COALESCE(d.attrs->>'pubmed_id', '') AS pmid,
              array_agg(DISTINCT t.node_id) AS matched_traits
            FROM okg.edges_live e
            JOIN okg.nodes_live d
              ON d.node_id = e.src AND d.subtype = 'dataset_metadata'
                 AND d.valid_to_gen IS NULL
            LEFT JOIN okg.nodes_live t
              ON t.node_id = e.dst AND t.subtype = 'trait'
                 AND t.valid_to_gen IS NULL
            WHERE e.valid_to_gen IS NULL
              AND e.edge_type = 'references'
              AND e.dst = ANY(%s)
            GROUP BY d.node_id, d.attrs
            ORDER BY n_total DESC NULLS LAST, d.node_id
        """, (trait_ids,))
        cols = [c.name for c in cur.description]
        rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    for r in rows:
        r["ancestry_codes"] = sorted(ancestry_codes(r["ancestry_scope"]))
    return rows


# ---------------------------- Filtering ----------------------------

def filter_rows(rows: list[dict], ancestry: Optional[list[str]],
                 min_n: Optional[int], sumstats_only: bool,
                 has_gcst: bool) -> list[dict]:
    out = rows
    if ancestry:
        wanted = set(c.upper() if c.upper() in {"EUR", "AFR", "EAS",
                                                  "AMR", "SAS"} else c
                      for c in ancestry)
        out = [r for r in out if set(r["ancestry_codes"]) & wanted]
    if min_n is not None:
        out = [r for r in out if int(r["n_total"] or 0) >= min_n]
    if sumstats_only:
        out = [r for r in out if r["sumstats_url"]]
    if has_gcst:
        out = [r for r in out if r["gcst"]]
    return out


# ---------------------------- Output ----------------------------

TSV_COLS = ["node_id", "gcst", "trait", "n_total", "ancestry_codes",
            "ancestry_scope", "build", "release_date", "author", "pmid",
            "sumstats_url"]


def emit_tsv(rows: list[dict], out=sys.stdout) -> None:
    out.write("\t".join(TSV_COLS) + "\n")
    for r in rows:
        cells = []
        for c in TSV_COLS:
            v = r.get(c, "")
            if isinstance(v, list):
                v = ",".join(map(str, v))
            cells.append(str(v) if v is not None else "")
        out.write("\t".join(cells) + "\n")


def emit_json_manifest(rows: list[dict], args, path: Path) -> None:
    manifest = {
        "skill": "okg-search", "version": "0.1",
        "query": {
            "trait": args.trait, "efo": args.efo,
            "seed_dataset": args.seed_dataset, "category": args.category,
            "ancestry": args.ancestry, "min_n": args.min_n,
            "sumstats_only": bool(args.sumstats_only),
            "limit": args.limit,
        },
        "n_results": len(rows),
        "results": rows,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    path.write_text(json.dumps(manifest, indent=2, default=str))
    print(f"manifest -> {path}", file=sys.stderr)


# ---------------------------- CLI ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--trait", help="Trait label (fuzzy ILIKE match)")
    grp.add_argument("--efo", help="EFO term id, e.g. EFO:0003767")
    grp.add_argument("--seed-dataset", dest="seed_dataset",
                     help="OKG dataset_metadata node id; find sibling datasets "
                          "sharing the same trait / EFO term")
    grp.add_argument("--category",
                     help="Curated trait category, e.g. 'immune'")

    ap.add_argument("--ancestry", type=lambda s: [x.strip() for x in s.split(",")],
                    help="Comma-separated controlled-vocab codes to keep: "
                         "EUR, AFR, EAS, AMR, SAS, admixed, unknown")
    ap.add_argument("--min-n", type=int, dest="min_n",
                    help="Minimum (n_cases + n_controls) to keep")
    ap.add_argument("--sumstats-only", action="store_true",
                    dest="sumstats_only",
                    help="STRICT: keep only datasets where the OKG has an "
                         "explicit summary_statistics_url attribute. Only "
                         "the curated ~5 hand-imported datasets have this "
                         "today; the bulk 165K scale-ingest namespace "
                         "doesn't. For 'fetchable via gwas-fetch' (every "
                         "GWAS Catalog accession qualifies), use --has-gcst "
                         "instead.")
    ap.add_argument("--has-gcst", action="store_true", dest="has_gcst",
                    help="Keep only datasets with a non-empty GCST "
                         "accession (effectively all dataset_metadata "
                         "nodes; this is the right filter when you plan "
                         "to feed results into `gwas-fetch GCST<X>`).")
    ap.add_argument("--expand-efo", dest="expand_efo",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="After resolving seed traits, walk to their EFO "
                         "terms and pull in every sibling trait mapping "
                         "to the same term. Catches both the curated and "
                         "scale-ingest trait namespaces. Default ON; pass "
                         "--no-expand-efo to disable.")
    ap.add_argument("--limit", type=int, default=50,
                    help="Max rows in stdout output (default 50; manifest "
                         "JSON keeps all)")
    ap.add_argument("--json", type=Path, default=None,
                    help="Write the full result set + query manifest to this path")
    ap.add_argument("--dsn", type=str, default=None,
                    help="Postgres DSN override (default $OKG_DSN or %s)" % DEFAULT_DSN)

    args = ap.parse_args()
    ensure_psycopg()

    with connect(args.dsn) as conn:
        if args.trait:
            trait_ids = resolve_trait_ids_by_label(conn, args.trait)
            mode_desc = f"trait label ~ {args.trait!r}"
        elif args.efo:
            trait_ids = resolve_trait_ids_by_efo(conn, args.efo)
            mode_desc = f"EFO term {args.efo}"
        elif args.seed_dataset:
            trait_ids = resolve_trait_ids_by_seed_dataset(conn, args.seed_dataset)
            mode_desc = f"sibling of {args.seed_dataset}"
        elif args.category:
            trait_ids = resolve_trait_ids_by_category(conn, args.category)
            mode_desc = f"category = {args.category!r}"
        else:
            ap.error("one of --trait / --efo / --seed-dataset / --category required")

        if not trait_ids:
            print(f"[okg-search] no trait nodes matched ({mode_desc})",
                  file=sys.stderr)
            return 1
        n_seed = len(trait_ids)
        if args.expand_efo:
            trait_ids = expand_trait_ids_via_efo(conn, trait_ids)
        print(f"[okg-search] matched {n_seed} seed trait nodes for "
              f"{mode_desc}"
              f"{'; expanded to ' + str(len(trait_ids)) + ' via EFO siblings' if args.expand_efo and len(trait_ids) > n_seed else ''}"
              f"; fetching datasets...", file=sys.stderr)

        rows = fetch_datasets_referencing_traits(conn, trait_ids)
    n_pre = len(rows)
    rows = filter_rows(rows, args.ancestry, args.min_n, args.sumstats_only,
                        args.has_gcst)
    print(f"[okg-search] {len(rows)} / {n_pre} datasets pass filters "
          f"(ancestry={args.ancestry}, min_n={args.min_n}, "
          f"sumstats_only={args.sumstats_only}, "
          f"has_gcst={args.has_gcst})", file=sys.stderr)

    # stdout = top --limit rows
    emit_tsv(rows[: args.limit], sys.stdout)
    if args.json:
        emit_json_manifest(rows, args, args.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
