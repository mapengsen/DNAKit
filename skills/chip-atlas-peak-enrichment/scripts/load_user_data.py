"""
Load user-supplied gene lists.

Accepted layouts
----------------
* Plain text  — one symbol per line, `#` starts a comment.
* CSV / TSV   — auto-detects a gene-name column from common headers.
* DE results  — DESeq2 / edgeR / limma exports work as CSV with a `gene`,
                `gene_name`, `symbol`, or similar column.

Inputs : a path
Outputs: list of gene symbol strings

Example
-------
    python load_user_data.py --file my_genes.txt
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


CANDIDATE_COLUMNS = (
    "gene", "Gene", "GENE", "gene_name", "geneName", "gene_symbol",
    "geneSymbol", "symbol", "Symbol", "id", "ID", "gene_id",
)


def load_user_data(path: str, column_name: Optional[str] = None) -> List[str]:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    extension = os.path.splitext(path)[1].lower()
    if extension == ".txt":
        genes = _read_plain_text(path)
    elif extension in {".csv", ".tsv", ".tab"}:
        sep = "," if extension == ".csv" else "\t"
        genes = _read_delimited(path, sep=sep, column_name=column_name)
    else:
        genes = _infer_format(path, column_name=column_name)
    if not genes:
        raise ValueError(f"No genes loaded from {path}")
    print(f"✓ Loaded {len(genes)} genes from {path}")
    return genes


def _read_plain_text(path: str) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            genes.append(stripped)
    return genes


def _read_delimited(path: str, sep: str, column_name: Optional[str]) -> List[str]:
    frame = pd.read_csv(path, sep=sep)
    if column_name is None:
        for candidate in CANDIDATE_COLUMNS:
            if candidate in frame.columns:
                column_name = candidate
                break
        if column_name is None:
            column_name = frame.columns[0]
            print(f"   inferred gene column: '{column_name}'")
    if column_name not in frame.columns:
        raise ValueError(f"Column '{column_name}' not in {list(frame.columns)}")
    return [
        symbol for symbol in frame[column_name].astype(str).tolist()
        if symbol and symbol.lower() != "nan"
    ]


def _infer_format(path: str, column_name: Optional[str]) -> List[str]:
    for sep in (",", "\t"):
        try:
            return _read_delimited(path, sep=sep, column_name=column_name)
        except Exception:
            continue
    return _read_plain_text(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--file", required=True)
    parser.add_argument("--column", default=None)
    args = parser.parse_args(argv)
    genes = load_user_data(args.file, column_name=args.column)
    print(f"✓ First 10: {genes[:10]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
