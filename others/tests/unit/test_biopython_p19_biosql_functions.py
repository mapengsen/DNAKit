"""Tests for the P19 BioSQL functions. A BioSQL SQLite schema file is needed (it is not shipped with Biopython):
set DNAKIT_BIOSQL_SQLITE_SCHEMA to the path of biosqldb-sqlite.sql; those tests are skipped without it."""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import pytest

pytest.importorskip("Bio")
pytest.importorskip("BioSQL")

from dnakit.databases import biopython_p19_biosql_functions as m
from dnakit.exceptions import BackendExecutionError, ConfigurationError

SCHEMA = os.environ.get("DNAKIT_BIOSQL_SQLITE_SCHEMA", "")
needs_schema = pytest.mark.skipif(
    not Path(SCHEMA).is_file(), reason="DNAKIT_BIOSQL_SQLITE_SCHEMA is not set"
)


@pytest.fixture()
def genbank(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    records = []
    for index, (name, text) in enumerate(
        [("REC1", "ATGGCTAAAGGTTTCCCAGAACTGGTTTAA"), ("REC2", "GGGGCCCCAAAATTTTACGT")]
    ):
        record = SeqRecord(
            Seq(text),
            id=f"ACC{index}.1",
            name=name,
            description=f"test {name}",
            annotations={"molecule_type": "DNA", "organism": "Testus x"},
        )
        record.features.append(
            SeqFeature(FeatureLocation(0, 9), type="CDS", qualifiers={"gene": [f"g{index}"]})
        )
        records.append(record)
    path = tmp_path / "a.gb"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SeqIO.write(records, str(path), "genbank")
    return path


@pytest.fixture()
def database(tmp_path, genbank):
    path = tmp_path / "bio.db"
    m.biosql_load_records(path, "ns1", genbank, schema_sql=SCHEMA, verbose=False)
    return path


@needs_schema
def test_load_creates_database_and_namespace(tmp_path, genbank) -> None:
    path = tmp_path / "new.db"
    result = m.biosql_load_records(path, "ns1", genbank, schema_sql=SCHEMA, verbose=False)
    assert result["records_loaded"] == 2 and result["created_namespace"] is True
    with pytest.raises(
        BackendExecutionError
    ):  # the same identifiers twice: rolled back, nothing changes
        m.biosql_load_records(path, "ns1", genbank, verbose=False)
    assert m.biosql_list_entries(path, verbose=False)["namespaces"] == [
        {"name": "ns1", "entry_count": 2}
    ]
    other = m.biosql_load_records(path, "ns2", genbank, verbose=False)
    assert other["records_in_namespace"] == 2
    with pytest.raises(ConfigurationError):
        m.biosql_load_records(
            path, "ns3", genbank, schema_sql=SCHEMA, verbose=False
        )  # would recreate


@needs_schema
def test_list_entries(database) -> None:
    namespaces = m.biosql_list_entries(database, verbose=False)["namespaces"]
    assert namespaces == [{"name": "ns1", "entry_count": 2}]
    entries = m.biosql_list_entries(database, namespace="ns1", verbose=False)
    assert entries["display_ids"] == ["REC1", "REC2"] and entries["entry_ids"] == [1, 2]
    assert (
        m.biosql_list_entries(database, namespace="ns1", max_entries=1, verbose=False)["truncated"]
        is True
    )
    with pytest.raises(ConfigurationError):
        m.biosql_list_entries(database, namespace="missing", verbose=False)


@needs_schema
def test_get_record_by_each_identifier(database, tmp_path) -> None:
    by_name = m.biosql_get_record(database, "ns1", "REC1", verbose=False)["record"]
    assert by_name["sequence"] == "ATGGCTAAAGGTTTCCCAGAACTGGTTTAA" and by_name["feature_types"] == [
        "CDS"
    ]
    assert by_name["id"] == "ACC0.1" and by_name["annotations"]["organism"] == "Testus x"
    assert (
        m.biosql_get_record(database, "ns1", "ACC1", id_type="accession", verbose=False)["record"][
            "name"
        ]
        == "REC2"
    )
    assert (
        m.biosql_get_record(database, "ns1", "ACC0.1", id_type="version", verbose=False)["record"][
            "name"
        ]
        == "REC1"
    )
    assert (
        m.biosql_get_record(database, "ns1", "REC2", id_type="name", verbose=False)["record"][
            "length"
        ]
        == 20
    )
    assert (
        m.biosql_get_record(database, "ns1", "ACC0.1", id_type="primary_id", verbose=False)[
            "record"
        ]["name"]
        == "REC1"
    )
    brief = m.biosql_get_record(database, "ns1", "REC1", include_sequence=False, verbose=False)[
        "record"
    ]
    assert brief["sequence"] is None and brief["length"] == 30
    saved = m.biosql_get_record(
        database, "ns1", "REC2", output_path=tmp_path / "o.gb", verbose=False
    )
    assert Path(saved["output"]).read_text().startswith("LOCUS       REC2")
    with pytest.raises(ConfigurationError):
        m.biosql_get_record(database, "ns1", "absent", verbose=False)
    with pytest.raises(ConfigurationError):
        m.biosql_get_record(database, "ns1", "REC1", max_sequence_length=10, verbose=False)
    with pytest.raises(ConfigurationError):
        m.biosql_get_record(database, "ns1", "REC1", id_type="isbn", verbose=False)


@needs_schema
def test_get_subsequence(database) -> None:
    full = "ATGGCTAAAGGTTTCCCAGAACTGGTTTAA"
    result = m.biosql_get_subsequence(database, "ns1", "REC1", 3, 12, verbose=False)
    assert result["sequence"] == full[3:12] and result["sequence_length"] == 30
    assert (
        m.biosql_get_subsequence(database, "ns1", "REC1", 0, 30, verbose=False)["sequence"] == full
    )
    assert (
        m.biosql_get_subsequence(
            database, "ns1", "ACC1.1", 19, 20, id_type="version", verbose=False
        )["sequence"]
        == "T"
    )
    with pytest.raises(ConfigurationError):
        m.biosql_get_subsequence(database, "ns1", "REC1", 5, 31, verbose=False)
    with pytest.raises(ConfigurationError):
        m.biosql_get_subsequence(database, "ns1", "REC1", 5, 5, verbose=False)


def test_validation_without_schema(tmp_path, genbank) -> None:
    with pytest.raises(ConfigurationError):
        m.biosql_list_entries(tmp_path / "missing.db", verbose=False)
    with pytest.raises(ConfigurationError):
        m.biosql_load_records(
            tmp_path / "x.db", "bad name!", genbank, schema_sql=tmp_path / "s.sql", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.biosql_load_records(tmp_path / "x.db", "ns", genbank, driver="oracle", verbose=False)
    with pytest.raises(ConfigurationError):
        m.biosql_load_records(
            tmp_path / "x.db", "ns", genbank, schema_sql=tmp_path / "missing.sql", verbose=False
        )
    assert not (tmp_path / "x.db").exists()  # a failed creation leaves no database file behind
