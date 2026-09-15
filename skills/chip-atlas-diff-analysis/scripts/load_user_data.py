"""
Load experiment accession lists from user-supplied files.

Supported formats
-----------------
* Plain text: one accession per line; lines starting with `#` are skipped.
* CSV / TSV : auto-detects a column from common headers
              (`experiment_id`, `srx`, `accession`, `sample_id`, `id`).
* Single-file two-group layout via `load_experiment_groups` — pass one
  delimited file with a group column and an accession column.

Inputs : file paths
Outputs: lists of strings; the validator only warns about non-standard IDs.

Example
-------
    python load_user_data.py --file group_a.txt
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import List, Optional, Tuple

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


ACCESSION_RE = re.compile(r"^(SRX|ERX|DRX|GSM)\d+$")
ID_COLUMN_CANDIDATES = (
    "experiment_id", "srx", "srx_id", "SRX", "accession",
    "sample_id", "ID", "id",
)


def load_user_data(path: str, column_name: Optional[str] = None) -> List[str]:
    """Return a list of accession strings from `path`."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    extension = os.path.splitext(path)[1].lower()
    if extension == ".txt":
        ids = _read_plain_text(path)
    elif extension in {".csv", ".tsv", ".tab"}:
        ids = _read_delimited(path, sep="," if extension == ".csv" else "\t", column_name=column_name)
    else:
        ids = _infer_format_and_read(path, column_name=column_name)

    valid, suspect = [], []
    for accession in ids:
        (valid if ACCESSION_RE.match(accession) else suspect).append(accession)

    if suspect:
        print(f"   ! {len(suspect)} ID(s) outside SRX/ERX/DRX/GSM, keeping anyway: {suspect[:5]}")
    if not valid and not suspect:
        raise ValueError(f"No accessions parsed from {path}")

    all_ids = valid + suspect
    print(f"✓ Loaded {len(all_ids)} accession(s) from {path}")
    return all_ids


def load_experiment_groups(
    file_a: Optional[str] = None,
    file_b: Optional[str] = None,
    group_file: Optional[str] = None,
    group_column: str = "group",
    id_column: Optional[str] = None,
) -> Tuple[List[str], List[str]]:
    """Load two groups, either from a pair of files or a single file with a group column."""
    if file_a and file_b:
        return load_user_data(file_a, column_name=id_column), load_user_data(file_b, column_name=id_column)

    if not group_file:
        raise ValueError("Provide either file_a + file_b, or group_file")

    if not os.path.exists(group_file):
        raise FileNotFoundError(group_file)

    extension = os.path.splitext(group_file)[1].lower()
    sep = "," if extension == ".csv" else "\t"
    frame = pd.read_csv(group_file, sep=sep)

    if group_column not in frame.columns:
        raise ValueError(
            f"Group column '{group_column}' not found. Columns: {list(frame.columns)}"
        )

    if id_column is None:
        for candidate in ID_COLUMN_CANDIDATES:
            if candidate in frame.columns:
                id_column = candidate
                break
        if id_column is None:
            for column in frame.columns:
                if column != group_column:
                    id_column = column
                    print(f"   Inferred ID column: '{id_column}'")
                    break

    groups = list(frame[group_column].dropna().unique())
    if len(groups) != 2:
        raise ValueError(f"Expected exactly two groups, found {len(groups)}: {groups}")

    name_a, name_b = sorted(groups)
    list_a = frame.loc[frame[group_column] == name_a, id_column].astype(str).tolist()
    list_b = frame.loc[frame[group_column] == name_b, id_column].astype(str).tolist()
    print(f"✓ Loaded group '{name_a}' ({len(list_a)} IDs) and '{name_b}' ({len(list_b)} IDs)")
    return list_a, list_b


def _read_plain_text(path: str) -> List[str]:
    ids: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            ids.append(stripped)
    if not ids:
        raise ValueError(f"No accessions in {path}")
    return ids


def _read_delimited(path: str, sep: str, column_name: Optional[str]) -> List[str]:
    frame = pd.read_csv(path, sep=sep)
    if column_name is None:
        for candidate in ID_COLUMN_CANDIDATES:
            if candidate in frame.columns:
                column_name = candidate
                break
        if column_name is None:
            column_name = frame.columns[0]
            print(f"   Inferred ID column: '{column_name}'")
    if column_name not in frame.columns:
        raise ValueError(f"Column '{column_name}' not in {list(frame.columns)}")
    values = frame[column_name].astype(str).tolist()
    return [v for v in values if v not in {"", "nan", "None"}]


def _infer_format_and_read(path: str, column_name: Optional[str]) -> List[str]:
    for sep in (",", "\t"):
        try:
            return _read_delimited(path, sep=sep, column_name=column_name)
        except Exception:
            continue
    return _read_plain_text(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--file", help="Single accession list (txt/csv/tsv)")
    parser.add_argument("--column", help="Column name for the accession column")
    parser.add_argument("--group-file", help="Delimited file with a group column")
    parser.add_argument("--group-column", default="group")
    parser.add_argument("--id-column", default=None)
    args = parser.parse_args(argv)

    if args.file:
        ids = load_user_data(args.file, column_name=args.column)
        print(f"✓ First 5: {ids[:5]}")
        return 0

    if args.group_file:
        a, b = load_experiment_groups(group_file=args.group_file, group_column=args.group_column, id_column=args.id_column)
        print(f"✓ Group A first 5: {a[:5]}; Group B first 5: {b[:5]}")
        return 0

    parser.error("Provide --file or --group-file")
    return 1


if __name__ == "__main__":
    sys.exit(main())
