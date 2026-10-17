"""Real behaviour of the io_formats format operations (map, ped, paf, newick, nexus, phylip, twobit, zarr)."""

from __future__ import annotations

import base64
import hashlib
import http.server
import json
import struct
import threading
from pathlib import Path

import pytest

from dnakit import io_formats
from dnakit.exceptions import ConfigurationError, DownloadError, InputFormatError
from dnakit.io_formats import _impl_zarr

FORMATS = ("map", "newick", "nexus", "paf", "ped", "phylip", "twobit", "zarr")
SPECIFIC_OPERATIONS = (
    "convert",
    "filter",
    "index",
    "merge",
    "read",
    "sort",
    "split",
    "stream",
    "summarize",
    "validate",
    "write",
)
GENERIC_OPERATIONS = ("convert", "filter", "index", "merge", "sort", "stream", "summarize")
ALL_FUNCTIONS = (
    [f"{op}_{fmt}" for op in SPECIFIC_OPERATIONS for fmt in FORMATS]
    + list(GENERIC_OPERATIONS)
    + ["resume_download"]
)


def call(name: str, data: object, **options: object) -> dict:
    """Call a public function and return the ``result`` of its (successful) envelope."""

    envelope = getattr(io_formats, name)(data, options=options or None, verbose=False)
    assert envelope["status"] == "ok", envelope
    assert envelope["module"] == "io_formats"
    assert envelope["operation"] == name
    return envelope["result"]  # type: ignore[no-any-return]


# ---------------------------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------------------------

MAP = "1\trs1\t0\t100\n1\trs2\t0.5\t50\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
PED = "F1 I1 0 0 1 2 A A C G\nF1 I2 0 0 2 1 A G 0 0\nF1 I3 I1 I2 1 -9 G G C C\n"
PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)
NEWICK = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"
PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"
NEXUS = """#NEXUS
[a comment]
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""


def twobit_bytes(records: list[tuple[str, int, list, list, bytes]], endian: str = "<") -> bytes:
    """Build a 2bit file by hand from (name, dna_size, n_blocks, mask_blocks, packed_bytes)."""

    def u32(*values: int) -> bytes:
        return struct.pack(f"{endian}{len(values)}I", *values)

    index_size = sum(1 + len(name.encode()) + 4 for name, *_ in records)
    offset = 16 + index_size
    index, body = b"", b""
    for name, size, n_blocks, mask_blocks, packed in records:
        index += bytes([len(name)]) + name.encode() + u32(offset)
        record = (
            u32(size, len(n_blocks))
            + u32(*[s for s, _ in n_blocks])
            + u32(*[n for _, n in n_blocks])
            + u32(len(mask_blocks))
            + u32(*[s for s, _ in mask_blocks])
            + u32(*[n for _, n in mask_blocks])
            + u32(0)
            + packed
        )
        body += record
        offset += len(record)
    return u32(0x1A412743, 0, len(records), 0) + index + body


# TCAGNNac: T=00 C=01 A=10 G=11 | N,N stored as T=00,00 | a=10 c=01  -> 0x1B 0x09
S1 = ("s1", 8, [(4, 2)], [(6, 2)], bytes([0x1B, 0x09]))
# ACGTACGTAC: A=10 C=01 G=11 T=00 -> 0x9C 0x9C ; AC + 2 padding bases -> 0x90
S2 = ("s2", 10, [], [], bytes([0x9C, 0x9C, 0x90]))
TWOBIT = twobit_bytes([S1, S2])


# ---------------------------------------------------------------------------------------------
# all 96 functions run (smoke matrix), generic dispatch
# ---------------------------------------------------------------------------------------------

SAMPLES = {
    "map": MAP,
    "ped": PED,
    "paf": PAF,
    "newick": NEWICK,
    "phylip": PHYLIP,
    "nexus": NEXUS,
    "twobit": TWOBIT,
}
MERGE_OPTIONS = {
    "map": {"duplicates": "first"},
    "ped": {"duplicates": "first"},
    "twobit": {"duplicates": "first"},
}
SPLIT_OPTIONS = {"phylip": {"site_size": 4}, "nexus": {"site_size": 4}}
CONVERT_TARGET = {
    "map": "bed",
    "ped": "tsv",
    "paf": "bed",
    "newick": "nexus",
    "phylip": "fasta",
    "nexus": "fasta",
    "twobit": "fasta",
}
WRITE_SAMPLE = {
    "map": [["1", "a", 5]],
    "ped": [{"fid": "F", "iid": "1", "genotypes": [["A", "C"]]}],
    "paf": [
        dict(
            query_name="q",
            query_length=10,
            query_start=0,
            query_end=10,
            strand="+",
            target_name="t",
            target_length=20,
            target_start=0,
            target_end=10,
            matches=9,
            block_length=10,
            mapq=5,
        )
    ],
    "newick": ["(A,B);"],
    "phylip": {"A": "ACGT", "B": "ACGA"},
    "nexus": {"A": "ACGT", "B": "ACGA"},
    "twobit": {"s": "ACGT"},
}


@pytest.mark.parametrize("fmt", sorted(SAMPLES))
@pytest.mark.parametrize("op", SPECIFIC_OPERATIONS)
def test_every_specific_function_runs_on_a_sample(op: str, fmt: str, tmp_path: Path) -> None:
    sample = SAMPLES[fmt]
    options: dict[str, object] = {}
    data: object = sample
    if op == "merge":
        data, options = [sample, sample], MERGE_OPTIONS.get(fmt, {})
    elif op == "split":
        options = SPLIT_OPTIONS.get(fmt, {})
    elif op == "convert":
        options = {"to": CONVERT_TARGET[fmt]}
    elif op == "write":
        data = WRITE_SAMPLE[fmt]
    result = call(f"{op}_{fmt}", data, **options)
    assert isinstance(result, dict)


@pytest.mark.parametrize("op", GENERIC_OPERATIONS)
@pytest.mark.parametrize("fmt", sorted(SAMPLES))
def test_generic_functions_detect_the_format(op: str, fmt: str) -> None:
    sample = SAMPLES[fmt]
    options: dict[str, object] = {}
    data: object = sample
    if op == "merge":
        data, options = [sample, sample], MERGE_OPTIONS.get(fmt, {})
    elif op == "convert":
        options = {"to": CONVERT_TARGET[fmt]}
    generic = call(op, data, **options)
    specific = call(f"{op}_{fmt}", data, **options)
    assert generic == specific


def test_generic_detection_covers_files_by_suffix_and_content(tmp_path: Path) -> None:
    path = tmp_path / "hits.paf"
    path.write_text(PAF)
    assert call("summarize", path)["alignment_count"] == 2
    gz = tmp_path / "v.map.gz"
    import gzip

    gz.write_bytes(gzip.compress(MAP.encode()))
    assert call("summarize", str(gz))["variant_count"] == 4
    assert call("summarize", PED, format="ped")["individual_count"] == 3
    with pytest.raises(InputFormatError):
        io_formats.summarize("completely unknown text", verbose=False)


def test_generic_fasta_is_a_convert_source_only() -> None:
    result = call("convert", ">a\nACGT\n>b\nACGA\n", format="fasta", to="phylip")
    assert result["text"] == " 2 4\na          ACGT\nb          ACGA\n"
    with pytest.raises(ConfigurationError):
        io_formats.filter(">a\nACGT\n", options={"format": "fasta"}, verbose=False)


# ---------------------------------------------------------------------------------------------
# map
# ---------------------------------------------------------------------------------------------


def test_map_read_validate_summarize_index() -> None:
    records = call("read_map", MAP)["records"]
    assert records[1] == {"chrom": "1", "id": "rs2", "cm": 0.5, "pos": 50}
    assert call("validate_map", MAP)["valid"] is True
    bad = call("validate_map", "1 rs1 0 100\n1 rs1 0 abc\n2 rs2 0\n1 rs1 0 100\n")
    assert bad["valid"] is False
    assert any("line 2" in e for e in bad["errors"]) and any("line 3" in e for e in bad["errors"])
    duplicate = call("validate_map", "1 a 0 5\n2 a 0 6\n1 b 0 5\n")
    assert duplicate["errors"] == ["line 2: duplicate variant ID 'a' (first on line 1)"]
    summary = call("summarize_map", MAP)
    assert summary["variant_count"] == 4 and summary["chromosome_count"] == 3
    assert summary["chromosomes"]["1"] == {"variant_count": 2, "min_pos": 50, "max_pos": 100}
    assert summary["sorted"] is False
    index = call("index_map", MAP)
    # "1\trs1\t0\t100\n" is 12 bytes, so rs2 starts at byte 12 and its line "1\trs2\t0.5\t50\n" has 13 bytes
    assert index["entries"]["rs2"] == {
        "line": 2,
        "offset": 12,
        "length": 13,
        "chrom": "1",
        "pos": 50,
    }
    assert index["chromosomes"]["1"]["first_line"] == 1


def test_map_filter_sort_merge_split_stream_convert_write() -> None:
    assert call("filter_map", MAP, chromosomes=["1"])["text"] == "1\trs1\t0\t100\n1\trs2\t0.5\t50\n"
    assert call("filter_map", MAP, region="1:60-100")["text"] == "1\trs1\t0\t100\n"
    assert (
        call("filter_map", MAP, exclude_ids="rs1,rs2", max_cm=0.0)["text"]
        == "2\trs3\t0\t10\nX\trs4\t0\t5\n"
    )
    assert (
        call("sort_map", MAP)["text"]
        == "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
    )
    assert call("sort_map", MAP, key="id", reverse=True)["text"].splitlines()[0] == "X\trs4\t0\t5"
    merged = call("merge_map", [MAP, "3 rs9 0 1\n"], sort=True)
    # natural chromosome order: 1, 2, 3, then X
    assert merged["record_count"] == 5 and merged["text"].splitlines()[-2:] == [
        "3\trs9\t0\t1",
        "X\trs4\t0\t5",
    ]
    with pytest.raises(InputFormatError):
        io_formats.merge_map([MAP, MAP], verbose=False)
    parts = call("split_map", MAP)["parts"]
    assert [(p["name"], p["record_count"]) for p in parts] == [("1", 2), ("2", 1), ("X", 1)]
    assert call("split_map", MAP, by="size", size=3)["parts"][1]["text"] == "X\trs4\t0\t5\n"
    page = call("stream_map", MAP, start=1, limit=2)
    assert (page["returned"], page["next_start"], page["total"], page["done"]) == (2, 3, 4, False)
    assert (
        call("convert_map", MAP, to="bed")["text"]
        == "1\t99\t100\trs1\n1\t49\t50\trs2\n2\t9\t10\trs3\nX\t4\t5\trs4\n"
    )
    assert call("convert_map", MAP, to="tsv")["text"].splitlines()[:2] == [
        "chrom\tid\tcm\tpos",
        "1\trs1\t0\t100",
    ]
    assert call("write_map", [["1", "a", 5], {"chrom": "2", "id": "b", "cm": 1.5, "pos": 7}])[
        "text"
    ] == ("1\ta\t0\t5\n2\tb\t1.5\t7\n")


# ---------------------------------------------------------------------------------------------
# ped
# ---------------------------------------------------------------------------------------------


def test_ped_read_validate_summarize_index() -> None:
    result = call("read_ped", PED)
    assert result["variant_count"] == 2
    assert result["records"][1]["genotypes"] == [["A", "G"], ["0", "0"]]
    assert call("validate_ped", PED)["valid"] is True
    ragged = call("validate_ped", "F I1 0 0 1 1 A A\nF I2 0 0 1 1 A A C C\n")
    assert ragged["valid"] is False and "line 2" in ragged["errors"][0]
    assert call("validate_ped", PED, variant_count=3)["errors"] == [
        "file has 2 variants but variant_count=3"
    ]
    odd = call("validate_ped", "F I1 0 0 7 1 A A\nF I1 0 0 1 1 A A\nF I2 X 0 1 1 A A\n")
    assert any("duplicate individual F:I1" in e for e in odd["errors"])
    assert any("sex code 7" in w for w in odd["warnings"]) and any(
        "'X'" in w for w in odd["warnings"]
    )
    summary = call("summarize_ped", PED)
    assert (
        summary["individual_count"] == 3
        and summary["family_count"] == 1
        and summary["founder_count"] == 2
    )
    assert summary["sex"] == {"male": 2, "female": 1, "unknown": 0}
    assert summary["phenotype"] == {"case": 1, "control": 1, "missing": 1, "quantitative": 0}
    assert summary["genotype_missing_rate"] == pytest.approx(
        1 / 6
    )  # one missing call in 3 x 2 genotypes
    assert summary["individual_missing_rate"]["max"] == 0.5
    index = call("index_ped", PED)
    assert index["entries"]["F1:I2"] == {"line": 2, "offset": 22, "length": 22}
    assert index["families"] == {"F1": ["I1", "I2", "I3"]}


def test_ped_filter_options_have_documented_meaning() -> None:
    def kept(**options: object) -> list[str]:
        text = call("filter_ped", PED, **options)["text"]
        return [line.split("\t")[1] for line in text.splitlines()]

    assert kept(sex="male") == ["I1", "I3"]
    assert kept(sex=2) == ["I2"]
    assert kept(phenotype="case") == ["I1"]
    assert kept(phenotype="missing") == ["I3"]
    assert kept(phenotype="observed") == ["I1", "I2"]
    assert kept(founders_only=True) == ["I1", "I2"]
    assert kept(max_missing_rate=0.4) == ["I1", "I3"]  # I2 misses 1 of 2 genotypes
    assert kept(individuals="I2,I3", exclude_individuals=["I3"]) == ["I2"]
    assert kept(families=["nope"]) == []
    columns = call("filter_ped", PED, variant_indices=[1])
    assert columns["variant_count"] == 1 and columns["kept_variant_indices"] == [1]
    assert columns["text"].splitlines()[0] == "F1\tI1\t0\t0\t1\t2\tC\tG"
    clean = call("filter_ped", PED, max_variant_missing_rate=0.0)
    assert (
        clean["variant_count"] == 1 and clean["text"].splitlines()[1] == "F1\tI2\t0\t0\t2\t1\tA\tG"
    )


def test_ped_sort_merge_split_stream_convert_write() -> None:
    order = call("sort_ped", PED, key="sex", reverse=True)["text"].splitlines()
    assert [line.split("\t")[1] for line in order] == ["I2", "I3", "I1"]
    other = "F2 J1 0 0 1 1 T T C C\n"
    merged = call("merge_ped", [PED, other])
    assert (
        merged["record_count"] == 4
        and merged["text"].splitlines()[-1] == "F2\tJ1\t0\t0\t1\t1\tT\tT\tC\tC"
    )
    with pytest.raises(InputFormatError):
        io_formats.merge_ped([PED, "F1 I1 0 0 1 1 A A C C\n"], verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.merge_ped([PED, "F9 I1 0 0 1 1 A A\n"], verbose=False)
    families = call("split_ped", PED + other)["parts"]
    assert [(p["name"], p["record_count"]) for p in families] == [("F1", 3), ("F2", 1)]
    assert [p["record_count"] for p in call("split_ped", PED, by="size", size=2)["parts"]] == [2, 1]
    assert call("stream_ped", PED, limit=2)["next_start"] == 2
    tsv = call("convert_ped", PED, to="tsv")["text"].splitlines()
    assert (
        tsv[0] == "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tV1\tV2"
        and tsv[2] == "F1\tI2\t0\t0\t2\t1\tA/G\t0/0"
    )
    written = call(
        "write_ped", [{"fid": "F", "iid": "1", "sex": 2, "pheno": 1, "genotypes": ["AC", "G/T"]}]
    )
    assert written["text"] == "F\t1\t0\t0\t2\t1\tA\tC\tG\tT\n"


# ---------------------------------------------------------------------------------------------
# paf
# ---------------------------------------------------------------------------------------------


def test_paf_read_validate_summarize_index() -> None:
    first = call("read_paf", PAF)["records"][0]
    assert (first["query_name"], first["target_start"], first["matches"], first["mapq"]) == (
        "q1",
        100,
        80,
        60,
    )
    assert first["tags"] == {"tp": {"type": "A", "value": "P"}, "cg": {"type": "Z", "value": "90M"}}
    assert call("validate_paf", PAF)["valid"] is True
    bad = call(
        "validate_paf",
        "q\t10\t0\t5\t+\tt\t10\t0\t5\t6\t5\t300\nq\t10\t0\t20\t+\tt\t10\t0\t5\t1\t5\t1\nshort\t1\n",
    )
    assert len(bad["errors"]) == 4  # line 1 has two problems (matches > block length, mapq)
    assert "matches" in bad["errors"][0] or "mapping quality" in bad["errors"][0]
    assert call("validate_paf", "q\t10\t0\t5\t+\tt\t10\t0\t5\t5\t5\t1\tbadtag\n")["valid"] is False
    summary = call("summarize_paf", PAF)
    assert (
        summary["alignment_count"] == 2
        and summary["aligned_bases"] == 140
        and summary["matching_bases"] == 120
    )
    assert summary["overall_identity"] == pytest.approx(120 / 140)
    assert summary["mean_identity"] == pytest.approx((80 / 90 + 40 / 50) / 2)
    assert summary["mapq"] == {"10": 1, "60": 1}
    assert summary["strand"] == {"+": 1, "-": 1}
    first_len = len(PAF.splitlines()[0]) + 1
    index = call("index_paf", PAF, by="target")
    assert index["entries"]["t2"] == [
        {"line": 2, "offset": first_len, "length": len(PAF) - first_len, "start": 0, "end": 50}
    ]


def test_paf_filter_sort_merge_split_stream_convert_write() -> None:
    def names(**options: object) -> list[str]:
        text = call("filter_paf", PAF, **options)["text"]
        return [line.split("\t")[0] for line in text.splitlines()]

    assert names(min_mapq=30) == ["q1"]
    assert names(min_identity=0.85) == ["q1"]  # 80/90 = 0.889 passes, 40/50 = 0.8 does not
    assert names(min_query_coverage=0.95) == ["q2"]  # 90/100 vs 50/50
    assert names(strand="-") == ["q2"]
    assert names(min_block_length=60) == ["q1"] and names(min_matches=41) == ["q1"]
    assert names(tags={"tp": "P"}) == ["q1"]
    assert names(target_region="t1:150-300") == ["q1"]
    assert (
        names(target_region="t1:190-300") == []
    )  # [100, 190) does not overlap the half-open [190, 300)
    assert names(query_region="q2:10-20", exclude_target_names=["t1"]) == ["q2"]
    assert names(query_names="q1,q2", exclude_query_names=["q1"]) == ["q2"]
    assert [
        line.split("\t")[0] for line in call("sort_paf", PAF, key="mapq")["text"].splitlines()
    ] == ["q2", "q1"]
    assert [
        line.split("\t")[0]
        for line in call("sort_paf", PAF, key="identity", reverse=True)["text"].splitlines()
    ] == ["q1", "q2"]
    merged = call("merge_paf", [PAF, PAF], deduplicate=True, sort_by="mapq")
    assert merged["record_count"] == 2 and merged["text"].splitlines()[0].startswith("q2")
    assert call("merge_paf", [PAF, PAF])["record_count"] == 4
    assert [p["name"] for p in call("split_paf", PAF, by="target")["parts"]] == ["t1", "t2"]
    assert [p["record_count"] for p in call("split_paf", PAF, by="size", size=1)["parts"]] == [1, 1]
    assert call("stream_paf", PAF, start=1)["records"][0]["query_name"] == "q2"
    assert (
        call("convert_paf", PAF, to="bed")["text"]
        == "t1\t100\t190\tq1\t60\t+\nt2\t0\t50\tq2\t10\t-\n"
    )
    tsv = call("convert_paf", PAF, to="tsv")["text"].splitlines()
    assert (
        tsv[0].endswith("\tcg\ttp") and tsv[1].endswith("\t90M\tP") and tsv[2].endswith("\t10\t\t")
    )
    record = dict(
        query_name="q",
        query_length=10,
        query_start=0,
        query_end=10,
        strand="-",
        target_name="t",
        target_length=20,
        target_start=2,
        target_end=12,
        matches=9,
        block_length=10,
        mapq=5,
        tags={"tp": {"type": "A", "value": "P"}, "NM": 1, "de": 0.5},
    )
    assert (
        call("write_paf", [record])["text"]
        == "q\t10\t0\t10\t-\tt\t20\t2\t12\t9\t10\t5\ttp:A:P\tNM:i:1\tde:f:0.5\n"
    )
    assert call("write_paf", PAF)["text"] == PAF


# ---------------------------------------------------------------------------------------------
# newick
# ---------------------------------------------------------------------------------------------


def test_newick_read_parses_quoting_comments_and_numbers() -> None:
    result = call("read_newick", "('a b':1e-3,[x]c:2,(d,e)f:0.5)root;")
    record = result["records"][0]
    assert record["tips"] == ["a b", "c", "d", "e"]
    assert record["newick"] == "('a b':0.001,c:2,(d,e)f:0.5)root;"
    nodes = record["nodes"]
    assert nodes[0] == {"id": 0, "parent": None, "name": "root", "length": None}
    assert nodes[1] == {"id": 1, "parent": 0, "name": "a b", "length": 0.001}


def test_newick_validate_summarize_index() -> None:
    assert call("validate_newick", NEWICK)["valid"] is True
    for text in ("((A,B);", "(A,B)", "(A,B));", "(A:x,B);", "(A,B)[;"):
        assert call("validate_newick", text)["valid"] is False, text
    dup = call("validate_newick", "(A,A);")
    assert dup["valid"] is False and "duplicated" in dup["errors"][0]
    negative = call("validate_newick", "(A:-1,B:1);")
    assert negative["valid"] is True and negative["warnings"] == [
        "tree 1: 1 negative branch length(s)"
    ]
    trees = call("summarize_newick", NEWICK)["trees"]
    assert (
        trees[0]["tip_count"] == 4
        and trees[0]["internal_node_count"] == 3
        and trees[0]["rooted"] is True
    )
    assert (
        trees[0]["height"] == 2.5
        and trees[0]["total_branch_length"] == 6.0
        and trees[0]["max_depth_edges"] == 2
    )
    assert (
        trees[1]["root_degree"] == 3
        and trees[1]["rooted"] is False
        and trees[1]["total_branch_length"] is None
    )
    assert trees[1]["fully_resolved"] is True
    index = call("index_newick", NEWICK)
    assert index["trees"][0] == {"index": 0, "offset": 0, "length": 30, "tip_count": 4}
    assert index["trees"][1]["offset"] == 31 and index["trees"][1]["length"] == 12
    assert index["tips"]["C"] == [0, 1]


def test_newick_filter_sort_merge_split_stream_convert_write() -> None:
    kept = call("filter_newick", NEWICK, keep_tips=["A", "B", "C"])["text"]
    assert (
        kept == "((A:1,B:2):0.5,C:1.5);\n(A,B,C);\n"
    )  # C inherits the 0.5 branch of the pruned (C,D) clade
    assert (
        call("filter_newick", NEWICK, drop_tips=["A"])["text"]
        == "(B:2.5,(C:1,D:1):0.5);\n(B,(C,D));\n"
    )
    collapsed = call("filter_newick", NEWICK, collapse_below=0.6)["text"]
    assert (
        collapsed == "(A:1,B:2,C:1,D:1);\n(A,B,(C,D));\n"
    )  # both 0.5 branches go; unlabeled-length (C,D) stays
    assert call("filter_newick", NEWICK, min_tips=5)["record_count"] == 0
    assert call("filter_newick", NEWICK, keep_tips=["A"], max_tips=0)["record_count"] == 0
    assert call("sort_newick", "(C,(B,A));", order="alphabetical")["text"] == "((A,B),C);\n"
    assert call("sort_newick", "((A,B,C),D);")["text"] == "(D,(A,B,C));\n"
    assert call("sort_newick", "(D,(A,B,C));", order="descending")["text"] == "((A,B,C),D);\n"
    merged = call(
        "merge_newick", [NEWICK, "(A,B);\n((A:1,B:2):0.5,(C:1,D:1):0.5);"], deduplicate=True
    )
    assert merged["text"] == "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n(A,B);\n"
    assert [p["text"] for p in call("split_newick", NEWICK)["parts"]] == [
        "((A:1,B:2):0.5,(C:1,D:1):0.5);\n",
        "(A,B,(C,D));\n",
    ]
    assert call("split_newick", NEWICK, trees_per_part=2)["part_count"] == 1
    assert call("stream_newick", NEWICK, limit=1)["next_start"] == 1
    nexus = call("convert_newick", "(A:1,B:2);", to="nexus")["text"]
    assert nexus == (
        "#NEXUS\n\nBEGIN TAXA;\n    DIMENSIONS NTAX=2;\n    TAXLABELS A B;\nEND;\n\n"
        "BEGIN TREES;\n    TREE tree1 = (A:1,B:2);\nEND;\n"
    )
    assert call("convert_newick", "(A:1,B:2);", to="tsv")["text"] == (
        "tree\tnode_id\tparent_id\tname\tlength\tis_tip\n0\t0\t\t\t\t0\n0\t1\t0\tA\t1\t1\n0\t2\t0\tB\t2\t1\n"
    )
    nested = {"children": [{"name": "x", "length": 1}, {"name": "a b", "length": 2.5}]}
    assert call("write_newick", nested)["text"] == "(x:1,'a b':2.5);\n"
    table = {
        "nodes": [
            {"id": 0, "parent": None},
            {"id": 1, "parent": 0, "name": "p", "length": 1},
            {"id": 2, "parent": 0, "name": "q"},
        ]
    }
    assert call("write_newick", [table, "(z);"])["text"] == "(p:1,q);\n(z);\n"


def test_newick_handles_very_deep_trees() -> None:
    text = "A"
    for i in range(5000):
        text = f"({text},T{i})"
    summary = call("summarize_newick", text + ";")["trees"][0]
    assert summary["tip_count"] == 5001 and summary["max_depth_edges"] == 5000
    assert call("filter_newick", text + ";", keep_tips=["A", "T4999"])["text"] == "(A,T4999);\n"


# ---------------------------------------------------------------------------------------------
# phylip
# ---------------------------------------------------------------------------------------------

INFORMATIVE = " 4 4\nA ACGT\nB ACGT\nC AGGA\nD AGGA\n"


def test_phylip_read_layouts() -> None:
    sequential = call("read_phylip", PHYLIP)
    assert [(r["name"], r["sequence"]) for r in sequential["records"]] == [
        ("A", "ACGTACGT"),
        ("B", "ACGTACGA"),
        ("C", "ACG-ACGT"),
    ]
    interleaved = call("read_phylip", " 2 8\nA ACGT\nB ACGA\n\nACGT\nTTTT\n")
    assert [r["sequence"] for r in interleaved["records"]] == ["ACGTACGT", "ACGATTTT"]
    spaced = call("read_phylip", " 1 8\nA ACGT ACGT\n")
    assert spaced["records"][0]["sequence"] == "ACGTACGT"
    strict = call("read_phylip", " 1 4\nTaxon one ACGT\n", input_strict_names=True)
    assert strict["records"][0] == {"name": "Taxon one", "sequence": "ACGT"}
    with pytest.raises(InputFormatError):
        io_formats.read_phylip(" 2 8\nA ACGT\nB ACGT\n", verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.read_phylip("not a header\nA ACGT\n", verbose=False)


def test_phylip_validate_summarize_index() -> None:
    assert call("validate_phylip", PHYLIP)["valid"] is True
    assert call("validate_phylip", " 2 8\nA ACGT\nB ACGT\n")["valid"] is False
    assert call("validate_phylip", " 2 4\nA ACGT\nA ACGT\n")["errors"] == [
        "duplicate taxon names: ['A']"
    ]
    protein = " 1 4\nA MKVL\n"
    assert (
        call("validate_phylip", protein)["valid"] is False
    )  # M, K, V, L are not nucleotide codes (default datatype=dna)
    assert call("validate_phylip", protein, datatype="protein")["valid"] is True
    summary = call("summarize_phylip", INFORMATIVE)
    assert (
        summary["constant_sites"],
        summary["variable_sites"],
        summary["parsimony_informative_sites"],
    ) == (2, 2, 2)
    assert summary["gc_content"] == 0.5  # 16 bases: A 6, C 2, G 6, T 2 -> GC = 8/16
    assert summary["composition"] == {"A": 6, "C": 2, "G": 6, "T": 2}
    plain = call("summarize_phylip", PHYLIP)
    assert (
        plain["constant_sites"],
        plain["variable_sites"],
        plain["parsimony_informative_sites"],
    ) == (7, 1, 0)
    assert plain["taxon_missing_fraction"]["max"] == 1 / 8
    assert plain["duplicate_sequence_count"] == 0
    assert call("index_phylip", PHYLIP)["entries"]["C"] == {
        "ordinal": 2,
        "length": 8,
        "residues": 7,
    }


def test_phylip_filter_sort_merge_split_stream_convert_write() -> None:
    assert call("filter_phylip", PHYLIP, exclude=["C"])["text"].splitlines()[0] == " 2 8"
    assert (
        call("filter_phylip", PHYLIP, max_missing_fraction=0.1)["record_count"] == 2
    )  # C is 1/8 = 0.125 missing
    trimmed = call("filter_phylip", PHYLIP, max_site_missing_fraction=0.3)["text"].splitlines()
    assert trimmed[0] == " 3 7" and trimmed[1].endswith(
        "ACGACGT"
    )  # the column with the gap (1/3 missing) is removed
    assert call("filter_phylip", PHYLIP, name_pattern="^[AB]$")["record_count"] == 2
    rows = call("sort_phylip", PHYLIP, key="missing", reverse=True)["text"].splitlines()
    assert [r.split()[0] for r in rows[1:]] == ["C", "B", "A"]
    first, second = " 2 4\nA ACGT\nB ACGT\n", " 2 3\nB TTT\nC GGG\n"
    merged = call("merge_phylip", [first, second])["text"].splitlines()
    assert [r.split()[-1] for r in merged[1:]] == ["ACGT???", "ACGTTTT", "????GGG"] and merged[
        0
    ] == " 3 7"
    assert (
        call("merge_phylip", [" 1 2\nA AC\n", " 1 2\nB GT\n"], mode="append")["record_count"] == 2
    )
    with pytest.raises(InputFormatError):
        io_formats.merge_phylip([first, first], options={"mode": "append"}, verbose=False)
    by_sites = call("split_phylip", PHYLIP, site_size=5)["parts"]
    assert [p["name"] for p in by_sites] == ["sites_1-5", "sites_6-8"]
    assert by_sites[1]["text"].splitlines()[1].endswith("CGT")
    partitions = call("split_phylip", PHYLIP, partitions={"g1": "1-4", "g2": "5-6,8-8"})["parts"]
    assert (
        partitions[1]["text"].splitlines()[1].split()[-1] == "ACT"
    )  # sites 5, 6 and 8 of ACGTACGT
    assert [p["record_count"] for p in call("split_phylip", PHYLIP, taxa_size=2)["parts"]] == [2, 1]
    with pytest.raises(ConfigurationError):
        io_formats.split_phylip(PHYLIP, verbose=False)
    assert call("stream_phylip", PHYLIP, start=2)["records"] == [
        {"name": "C", "sequence": "ACG-ACGT"}
    ]
    assert (
        call("convert_phylip", PHYLIP, to="fasta")["text"]
        == ">A\nACGTACGT\n>B\nACGTACGA\n>C\nACG-ACGT\n"
    )
    nexus = call("convert_phylip", PHYLIP, to="nexus")["text"]
    assert "DIMENSIONS NTAX=3 NCHAR=8;" in nexus and "    C  ACG-ACGT\n" in nexus
    written = call("write_phylip", ">a\nACGTAC\n>b\nACGTAA\n", layout="interleaved", block_size=4)[
        "text"
    ]
    assert written == " 2 6\na          ACGT\nb          ACGT\n\nAC\nAA\n"
    assert call("read_phylip", written)["records"][1]["sequence"] == "ACGTAA"
    assert (
        call("write_phylip", {"A": "AC", "B": "GT"}, strict_names=True)["text"]
        == " 2 2\nA         AC\nB         GT\n"
    )
    with pytest.raises(InputFormatError):
        io_formats.write_phylip({"a name": "AC"}, verbose=False)


# ---------------------------------------------------------------------------------------------
# nexus
# ---------------------------------------------------------------------------------------------


def test_nexus_read_validate_summarize_index() -> None:
    result = call("read_nexus", NEXUS)
    assert result["taxa"] == ["A", "B", "C d"] and result["datatype"] == "dna"
    sequences = [r for r in result["records"] if r["kind"] == "sequence"]
    assert [r["sequence"] for r in sequences] == ["ACGTACGT", "ACGTACGA", "ACG-ACGT"]
    tree = [r for r in result["records"] if r["kind"] == "tree"][0]
    assert (
        tree["name"] == "t1"
        and tree["rooted"] is True
        and tree["newick"] == "((A:1,B:1):1,'C d':2);"
    )
    assert call("validate_nexus", NEXUS)["valid"] is True
    assert call("validate_nexus", "BEGIN DATA;")["valid"] is False
    wrong = NEXUS.replace("NTAX=3", "NTAX=4").replace("NCHAR=8", "NCHAR=9")
    assert call("validate_nexus", wrong)["valid"] is False
    extra = call("validate_nexus", NEXUS.replace("(1:1,2:1)", "(1:1,9:1)"))
    assert extra["valid"] is True and "9" in extra["warnings"][0]
    summary = call("summarize_nexus", NEXUS)
    assert (summary["constant_sites"], summary["variable_sites"]) == (7, 1) and summary[
        "tree_count"
    ] == 1
    assert summary["trees"][0]["tip_count"] == 3
    index = call("index_nexus", NEXUS)
    assert index["entries"]["C d"]["residues"] == 7 and index["trees"] == {"t1": 0}


def test_nexus_matchchar_and_interleave() -> None:
    text = """#NEXUS
BEGIN DATA;
 DIMENSIONS NTAX=2 NCHAR=8;
 FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. INTERLEAVE=YES;
 MATRIX
  A ACGT
  B .GG.
  A ACGT
  B ....
 ;
END;
"""
    assert [r["sequence"] for r in call("read_nexus", text)["records"]] == ["ACGTACGT", "AGGTACGT"]


def test_nexus_filter_sort_merge_split_stream_convert_write() -> None:
    filtered = call("filter_nexus", NEXUS, exclude=["B"])["text"]
    assert "NTAX=2" in filtered and "TREE t1 = [&R] (A:2,'C d':2);" in filtered
    assert call("read_nexus", filtered)["taxa"] == ["A", "C d"]
    assert (
        call("filter_nexus", NEXUS, max_missing_fraction=0.1)["record_count"] == 3
    )  # 2 taxa + 1 tree
    ordered = call("sort_nexus", NEXUS, key="name", reverse=True)["text"]
    assert ordered.index("'C d'") < ordered.index("    B ") < ordered.index("    A ")
    merged = call("merge_nexus", [NEXUS, NEXUS])["text"]
    assert "NCHAR=16" in merged and "A      ACGTACGTACGTACGT" in merged
    assert merged.count("TREE ") == 2 and "TREE t1_2 = " in merged
    parts = call("split_nexus", NEXUS, site_size=4)["parts"]
    assert [p["name"] for p in parts] == ["sites_1-4", "sites_5-8"] and "NCHAR=4" in parts[0][
        "text"
    ]
    only_trees = "#NEXUS\nBEGIN TREES;\nTREE a = (A,B);\nTREE b = (C,D);\nEND;\n"
    tree_parts = call("split_nexus", only_trees, trees_per_part=1)["parts"]
    assert [p["name"] for p in tree_parts] == [
        "trees_0001",
        "trees_0002",
    ] and "TREE b = (C,D);" in tree_parts[1]["text"]
    assert call("stream_nexus", NEXUS, start=3)["records"][0]["kind"] == "tree"
    assert call("convert_nexus", NEXUS, to="newick")["text"] == "((A:1,B:1):1,'C d':2);\n"
    assert call("convert_nexus", NEXUS, to="fasta")["text"].startswith(">A\nACGTACGT\n>B\n")
    assert call("convert_nexus", NEXUS, to="phylip", strict_names=True)["record_count"] == 3
    with pytest.raises(InputFormatError):
        io_formats.convert_nexus(
            NEXUS, options={"to": "phylip"}, verbose=False
        )  # 'C d' has a space
    expected = (
        "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=2 NCHAR=4;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
        "    MATRIX\n    A  ACGT\n    B  ACGA\n    ;\nEND;\n"
    )
    assert call("write_nexus", {"A": "ACGT", "B": "ACGA"})["text"] == expected
    assert (
        call("write_nexus", PHYLIP.replace("A  ", "A ").replace("B  ", "B ").replace("C  ", "C "))[
            "record_count"
        ]
        == 3
    )
    protein = call("write_nexus", [{"name": "p", "sequence": "MKV"}], datatype="protein")["text"]
    assert "DATATYPE=PROTEIN" in protein


# ---------------------------------------------------------------------------------------------
# twobit
# ---------------------------------------------------------------------------------------------


def decode_payload(result: dict) -> bytes:
    return base64.b64decode(result["data_base64"])


def test_twobit_writer_matches_the_hand_built_file() -> None:
    assert decode_payload(call("write_twobit", {"s1": "TCAGNNac", "s2": "ACGTACGTAC"})) == TWOBIT
    assert decode_payload(call("write_twobit", ">s1\nTCAGNNac\n>s2\nACGTAC\nGTAC\n")) == TWOBIT
    assert (
        decode_payload(
            call("write_twobit", [{"name": "s1", "sequence": "TCAGNNac"}, ("s2", "ACGTACGTAC")])
        )
        == TWOBIT
    )
    big = twobit_bytes([S1, S2], endian=">")
    assert (
        decode_payload(
            call("write_twobit", {"s1": "TCAGNNac", "s2": "ACGTACGTAC"}, byte_order="big")
        )
        == big
    )
    assert big[:4] == b"\x1a\x41\x27\x43" and TWOBIT[:4] == b"\x43\x27\x41\x1a"
    # IUPAC ambiguity codes become N blocks (like faToTwoBit)
    iupac = decode_payload(call("write_twobit", {"s": "ACRYGT"}))
    # A C [R Y -> N block] G T: A=10 C=01 T=00 T=00 | G=11 T=00 + 2 padding bases -> 0x90 0xC0
    assert iupac == twobit_bytes([("s", 6, [(2, 2)], [], bytes([0x90, 0xC0]))])
    assert call("read_twobit", iupac)["records"][0]["sequence"] == "ACNNGT"


def test_twobit_read_and_regions() -> None:
    records = call("read_twobit", TWOBIT)["records"]
    assert [(r["name"], r["sequence"]) for r in records] == [
        ("s1", "TCAGNNac"),
        ("s2", "ACGTACGTAC"),
    ]
    assert (
        call("read_twobit", TWOBIT, names=["s1"], soft_mask=False)["records"][0]["sequence"]
        == "TCAGNNAC"
    )
    assert (
        call("read_twobit", TWOBIT, names=["s1"], start=2, end=6)["records"][0]["sequence"]
        == "AGNN"
    )
    assert (
        call("read_twobit", TWOBIT, names=["s2"], start=5, end=9)["records"][0]["sequence"]
        == "CGTA"
    )
    assert (
        call("read_twobit", twobit_bytes([S1, S2], endian=">"))["records"][1]["sequence"]
        == "ACGTACGTAC"
    )
    with pytest.raises(ConfigurationError):
        io_formats.read_twobit(TWOBIT, options={"names": ["zz"]}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.read_twobit(TWOBIT, options={"start": 1}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.read_twobit(TWOBIT, options={"max_bases": 5}, verbose=False)


def test_twobit_validate_summarize_index_stream() -> None:
    ok = call("validate_twobit", TWOBIT)
    assert ok["valid"] is True and ok["sequence_count"] == 2 and ok["byte_order"] == "little"
    assert call("validate_twobit", b"NOTA2BITFILE" * 4)["valid"] is False
    assert call("validate_twobit", TWOBIT[:20])["valid"] is False
    overlapping = twobit_bytes([("x", 10, [(0, 5), (3, 5)], [], bytes(3))])
    assert any("overlap" in e for e in call("validate_twobit", overlapping)["errors"])
    outside = twobit_bytes([("x", 4, [(2, 5)], [], bytes(1))])
    assert any("extends past" in e for e in call("validate_twobit", outside)["errors"])
    assert call("validate_twobit", TWOBIT + b"")["valid"] is True
    bumped = bytearray(TWOBIT)
    bumped[8:12] = struct.pack(
        "<I", 3
    )  # claims 3 sequences, only 2 index entries fit before the data
    assert call("validate_twobit", bytes(bumped))["valid"] is False
    summary = call("summarize_twobit", TWOBIT)
    s1, s2 = summary["sequences"]
    assert (s1["A"], s1["C"], s1["G"], s1["T"], s1["N"], s1["soft_masked"]) == (2, 2, 1, 1, 2, 2)
    assert s1["gc_fraction"] == 0.5  # (C 2 + G 1) / (A 2 + C 2 + G 1 + T 1)
    assert (s2["A"], s2["C"], s2["G"], s2["T"], s2["N"]) == (3, 3, 2, 2, 0) and s2[
        "gc_fraction"
    ] == 0.5
    assert (
        summary["total_length"] == 18
        and summary["n_bases"] == 2
        and summary["gc_fraction"] == pytest.approx(8 / 16)
    )
    # header 16 bytes + two index entries of 1+2+4 bytes -> the first record starts at byte 30
    entries = call("index_twobit", TWOBIT)["entries"]
    assert entries["s1"] == {
        "offset": 30,
        "length": 8,
        "n_block_count": 1,
        "mask_block_count": 1,
        "packed_offset": 62,
    }
    assert (
        entries["s2"]["offset"] == 30 + 4 + 4 + 8 + 4 + 8 + 4 + 2
        and entries["s2"]["n_block_count"] == 0
    )
    page = call("stream_twobit", TWOBIT, window=4, start=1, limit=3)
    assert [(r["name"], r["start"], r["end"], r["sequence"]) for r in page["records"]] == [
        ("s1", 4, 8, "NNac"),
        ("s2", 0, 4, "ACGT"),
        ("s2", 4, 8, "ACGT"),
    ]
    assert page["total"] == 5 and page["next_start"] == 4


def test_twobit_filter_sort_merge_split_convert(tmp_path: Path) -> None:
    def names(result: dict) -> list[str]:
        return [r["name"] for r in call("read_twobit", decode_payload(result))["records"]]

    assert names(call("filter_twobit", TWOBIT, min_length=9)) == ["s2"]
    assert names(call("filter_twobit", TWOBIT, max_n_fraction=0.1)) == ["s2"]
    assert names(call("filter_twobit", TWOBIT, names=["s1"])) == ["s1"]
    assert names(call("filter_twobit", TWOBIT, name_pattern="2$", exclude=["s9"])) == ["s2"]
    assert names(call("sort_twobit", TWOBIT, key="length", reverse=True)) == ["s2", "s1"]
    chrom = decode_payload(call("write_twobit", {"chr10": "AC", "chr2": "GT", "chrX": "A"}))
    assert names(call("sort_twobit", chrom, key="natural")) == ["chr2", "chr10", "chrX"]
    assert names(call("sort_twobit", chrom)) == ["chr10", "chr2", "chrX"]
    # filter/sort copy records bit-exactly
    assert decode_payload(call("filter_twobit", TWOBIT, names=["s1", "s2"])) == TWOBIT
    other = decode_payload(call("write_twobit", {"s3": "GGCC"}))
    merged = call("merge_twobit", [TWOBIT, other])
    assert names(merged) == ["s1", "s2", "s3"]
    with pytest.raises(InputFormatError):
        io_formats.merge_twobit([TWOBIT, TWOBIT], verbose=False)
    assert call("merge_twobit", [TWOBIT, TWOBIT], duplicates="first")["record_count"] == 2
    parts = call("split_twobit", TWOBIT)["parts"]
    assert [p["name"] for p in parts] == ["s1", "s2"]
    assert (
        call("read_twobit", base64.b64decode(parts[1]["data_base64"]))["records"][0]["sequence"]
        == "ACGTACGTAC"
    )
    written = call("split_twobit", TWOBIT, output_dir=str(tmp_path))["parts"]
    assert Path(written[0]["output_path"]).read_bytes() == twobit_bytes([S1])
    assert (
        call("convert_twobit", TWOBIT, to="fasta", line_width=4)["text"]
        == ">s1\nTCAG\nNNac\n>s2\nACGT\nACGT\nAC\n"
    )
    assert (
        call("convert_twobit", TWOBIT, to="fasta", names=["s2"], line_width=0)["text"]
        == ">s2\nACGTACGTAC\n"
    )
    again = call("convert_twobit", twobit_bytes([S1], endian=">"), to="twobit")
    assert decode_payload(again) == twobit_bytes([S1])
    assert call("convert_twobit", TWOBIT, to="json")["text"].count('"n_blocks"') == 2
    assert (
        decode_payload(
            call("convert", ">s1\nTCAGNNac\n>s2\nACGTACGTAC\n", format="fasta", to="twobit")
        )
        == TWOBIT
    )


def test_twobit_works_on_files(tmp_path: Path) -> None:
    path = tmp_path / "genome.2bit"
    out = call("write_twobit", {"s1": "TCAGNNac"}, output=str(path))
    assert out["output_path"] == str(path) and path.read_bytes() == twobit_bytes([S1])
    with pytest.raises(ConfigurationError):
        io_formats.write_twobit({"s1": "A"}, options={"output": str(path)}, verbose=False)
    assert (
        call("write_twobit", {"s1": "A"}, output=str(path), overwrite=True)["bytes"]
        == path.stat().st_size
    )
    path.write_bytes(TWOBIT)
    assert call("read_twobit", str(path), names=["s2"])["records"][0]["sequence"] == "ACGTACGTAC"
    assert call("summarize", path)["sequence_count"] == 2
    assert call("sort_twobit", path, key="length", output=str(tmp_path / "sorted.2bit"))[
        "output_path"
    ]
    assert call("validate_twobit", path)["valid"] is True


def test_twobit_rejects_bad_input() -> None:
    with pytest.raises(InputFormatError):
        io_formats.read_twobit(b"definitely not 2bit data", verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.read_twobit("a text string", verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.write_twobit({"": "ACGT"}, verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.write_twobit({"a": "AC1T"}, verbose=False)
    version1 = bytearray(TWOBIT)
    version1[4:8] = struct.pack("<I", 1)
    with pytest.raises(InputFormatError):
        io_formats.read_twobit(bytes(version1), verbose=False)


# ---------------------------------------------------------------------------------------------
# zarr
# ---------------------------------------------------------------------------------------------

zarr_available = pytest.mark.skipif(
    __import__("importlib").util.find_spec("zarr") is None, reason="zarr is an optional dependency"
)


@pytest.fixture
def store(tmp_path: Path) -> str:
    path = str(tmp_path / "a.zarr")
    call(
        "write_zarr",
        {"x": [[1, 2], [3, 4], [5, 6]], "y": [30, 10, 20], "g/z": [1.5, 2.5, 3.5]},
        output=path,
        chunks=2,
        attrs={"k": "v"},
    )
    return path


@zarr_available
def test_zarr_write_read_summarize_index(store: str) -> None:
    result = call("read_zarr", store)
    assert result["attrs"] == {"k": "v"} and result["groups"] == ["g"]
    assert result["arrays"]["x"]["data"] == [[1, 2], [3, 4], [5, 6]]
    assert result["arrays"]["x"]["shape"] == [3, 2] and result["arrays"]["x"]["chunks"] == [2, 2]
    assert (
        result["arrays"]["g/z"]["data"] == [1.5, 2.5, 3.5]
        and result["arrays"]["g/z"]["dtype"] == "float64"
    )
    assert call("read_zarr", store, array="y", start=1, stop=3)["arrays"]["y"]["data"] == [10, 20]
    assert "data" not in call("read_zarr", store, include_data=False)["arrays"]["x"]
    with pytest.raises(ConfigurationError):
        io_formats.read_zarr(store, options={"max_elements": 3}, verbose=False)
    summary = call("summarize_zarr", store)
    assert (
        summary["array_count"] == 3
        and summary["group_count"] == 1
        and summary["total_elements"] == 6 + 3 + 3
    )
    assert summary["total_nbytes"] == 6 * 8 + 3 * 8 + 3 * 8 and summary["total_chunks"] == 2 + 2 + 2
    index = call("index_zarr", store)["entries"]
    assert (
        index["x"]["chunk_grid"] == [2, 1]
        and index["x"]["nchunks"] == 2
        and index["g/z"]["chunks"] == [2]
    )


@zarr_available
def test_zarr_validate_detects_corruption(store: str, tmp_path: Path) -> None:
    assert call("validate_zarr", store)["valid"] is True
    for chunk in (Path(store) / "x" / "c").rglob("*"):
        if chunk.is_file():
            chunk.write_bytes(b"this is not a zstd frame")
    broken = call("validate_zarr", store)
    assert broken["valid"] is False and broken["errors"][0].startswith(
        "x: first chunk cannot be decoded"
    )
    assert call("validate_zarr", str(tmp_path))["valid"] is False  # a directory that is not a store


@zarr_available
def test_zarr_stream_filter_sort_merge_split_convert(store: str, tmp_path: Path) -> None:
    page = call("stream_zarr", store, array="x", start=1, limit=5)
    assert page["records"] == [[3, 4], [5, 6]] and page["done"] is True and page["total"] == 3
    with pytest.raises(ConfigurationError):
        io_formats.stream_zarr(store, verbose=False)  # several arrays and no 'array'

    def read(path: str) -> dict:
        return {k: v["data"] for k, v in call("read_zarr", path)["arrays"].items()}

    sorted_out = str(tmp_path / "sorted.zarr")
    done = call("sort_zarr", store, by="y", output=sorted_out)
    assert done["permuted_arrays"] == ["g/z", "x", "y"]
    assert read(sorted_out) == {
        "g/z": [2.5, 3.5, 1.5],
        "x": [[3, 4], [5, 6], [1, 2]],
        "y": [10, 20, 30],
    }
    assert read(str(tmp_path / "sorted.zarr"))["y"] == [10, 20, 30]
    descending = str(tmp_path / "desc.zarr")
    call("sort_zarr", store, by="y", reverse=True, output=descending)
    assert read(descending)["y"] == [30, 20, 10]
    filtered = str(tmp_path / "f.zarr")
    assert (
        call("filter_zarr", store, arrays=["x", "y"], rows="1:3", output=filtered)["kept_rows"] == 2
    )
    assert read(filtered) == {"x": [[3, 4], [5, 6]], "y": [10, 20]}
    by_prefix = str(tmp_path / "prefix.zarr")
    call("filter_zarr", store, arrays=["g"], rows=[0, 2], output=by_prefix)
    assert read(by_prefix) == {"g/z": [1.5, 3.5]}
    dropped = str(tmp_path / "dropped.zarr")
    call("filter_zarr", store, exclude=["x", "g"], output=dropped)
    assert read(dropped) == {"y": [30, 10, 20]}
    merged = str(tmp_path / "m.zarr")
    assert call("merge_zarr", [filtered, filtered], output=merged)["input_count"] == 2
    assert read(merged) == {"x": [[3, 4], [5, 6], [3, 4], [5, 6]], "y": [10, 20, 10, 20]}
    with pytest.raises(InputFormatError):
        io_formats.merge_zarr(
            [filtered, dropped], options={"output": str(tmp_path / "bad.zarr")}, verbose=False
        )
    parts = call("split_zarr", store, rows_per_part=2, output_dir=str(tmp_path / "parts"))["parts"]
    assert [p["row_count"] for p in parts] == [2, 1]
    assert read(parts[0]["output_path"])["y"] == [30, 10] and read(parts[1]["output_path"])[
        "x"
    ] == [[5, 6]]
    assert call("convert_zarr", store, to="tsv", array="x")["text"] == "1\t2\n3\t4\n5\t6\n"
    assert call("convert_zarr", store, to="tsv", array="y")["text"] == "30\n10\n20\n"
    as_json = json.loads(call("convert_zarr", store, to="json")["text"])
    assert as_json["arrays"]["y"]["data"] == [30, 10, 20] and as_json["attrs"] == {"k": "v"}
    import io

    import numpy as np

    archive = np.load(io.BytesIO(decode_payload(call("convert_zarr", store, to="npz"))))
    assert archive["x"].tolist() == [[1, 2], [3, 4], [5, 6]] and archive["g__z"].tolist() == [
        1.5,
        2.5,
        3.5,
    ]
    copy = str(tmp_path / "copy.zarr")
    call("convert_zarr", store, to="zarr", output=copy)
    assert read(copy) == read(store)
    assert call("summarize", store)["array_count"] == 3  # generic dispatch recognises the directory


@zarr_available
def test_zarr_safety_and_errors(store: str, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        io_formats.sort_zarr(store, options={"by": "y"}, verbose=False)  # output is required
    with pytest.raises(ConfigurationError):
        io_formats.sort_zarr(store, options={"by": "y", "output": store}, verbose=False)  # exists
    plain = tmp_path / "keep"
    plain.mkdir()
    (plain / "precious.txt").write_text("x")
    with pytest.raises(ConfigurationError):
        io_formats.write_zarr(
            {"a": [1]}, options={"output": str(plain), "overwrite": True}, verbose=False
        )
    assert (plain / "precious.txt").exists()
    again = call("write_zarr", {"a": [1, 2]}, output=store, overwrite=True)
    assert again["output_path"] == store
    assert list(call("read_zarr", store)["arrays"]) == ["a"]
    with pytest.raises(InputFormatError):
        io_formats.write_zarr(
            {"a": ["x", "y"]}, options={"output": str(tmp_path / "s.zarr")}, verbose=False
        )
    with pytest.raises(InputFormatError):
        io_formats.read_zarr(str(tmp_path / "missing.zarr"), verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.sort_zarr(
            store, options={"by": "nope", "output": str(tmp_path / "n.zarr")}, verbose=False
        )


@pytest.mark.parametrize("operation", SPECIFIC_OPERATIONS)
def test_zarr_functions_report_not_implemented_without_the_optional_packages(
    operation: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    def missing() -> None:
        raise _impl_zarr.OptionalBackendMissing(
            "pip install 'zarr>=3' numpy", "zarr is not installed"
        )

    monkeypatch.setattr(_impl_zarr, "_libs", missing)
    data: object = "some.zarr"
    options: dict[str, object] = {"to": "json"} if operation == "convert" else {}
    if operation == "merge":
        data = ["a.zarr", "b.zarr"]
    envelope = getattr(io_formats, f"{operation}_zarr")(
        data, options=options or None, verbose=False
    )
    assert envelope["status"] == "not_implemented"
    assert envelope["required_dependency"] == "pip install 'zarr>=3' numpy"


# ---------------------------------------------------------------------------------------------
# resume_download against a local HTTP server
# ---------------------------------------------------------------------------------------------

PAYLOAD = bytes(range(256)) * 40  # 10240 bytes
ETAG = '"v1"'


class RangeHandler(http.server.BaseHTTPRequestHandler):
    mode = "normal"  # normal | cut | ignore_range
    cut_at = 0
    etag = ETAG
    body = PAYLOAD
    requests: list[dict[str, str | None]] = []

    def log_message(self, *args: object) -> None:  # silence
        pass

    def do_GET(self) -> None:  # noqa: N802
        cls = type(self)
        if self.path == "/redirect":
            self.send_response(302)
            self.send_header("Location", "/file.bin")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        cls.requests.append(
            {"range": self.headers.get("Range"), "if_range": self.headers.get("If-Range")}
        )
        body, start = cls.body, 0
        header = self.headers.get("Range")
        honour = (
            header
            and cls.mode != "ignore_range"
            and self.headers.get("If-Range") in (None, cls.etag)
        )
        if honour:
            start = int(header.split("=")[1].rstrip("-"))
            if start >= len(body):
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{len(body)}")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
        else:
            self.send_response(200)
        data = body[start:]
        self.send_header("Content-Length", str(len(data)))
        self.send_header("ETag", cls.etag)
        self.end_headers()
        if cls.mode == "cut" and not honour:
            self.wfile.write(data[: cls.cut_at])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(data)


@pytest.fixture
def server():
    RangeHandler.mode, RangeHandler.cut_at, RangeHandler.etag, RangeHandler.body = (
        "normal",
        0,
        ETAG,
        PAYLOAD,
    )
    RangeHandler.requests = []
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}/file.bin"
    httpd.shutdown()
    httpd.server_close()


def test_resume_download_full_transfer(server: str, tmp_path: Path) -> None:
    target = tmp_path / "out.bin"
    result = call("resume_download", server, output=str(target))
    assert target.read_bytes() == PAYLOAD and result["bytes"] == len(PAYLOAD)
    assert result["sha256"] == hashlib.sha256(PAYLOAD).hexdigest()
    assert (
        result["resumed_from"] == 0 and result["restarted"] is False and result["complete"] is True
    )
    assert (
        not (tmp_path / "out.bin.part").exists() and not (tmp_path / "out.bin.part.meta").exists()
    )
    assert call("resume_download", {"url": server, "output": str(tmp_path / "b.bin")})[
        "bytes"
    ] == len(PAYLOAD)


def test_resume_download_follows_http_redirects(server: str, tmp_path: Path) -> None:
    target = tmp_path / "out.bin"
    result = call("resume_download", server.replace("/file.bin", "/redirect"), output=str(target))
    assert target.read_bytes() == PAYLOAD and result["bytes"] == len(PAYLOAD)


def test_resume_download_resumes_an_interrupted_transfer(server: str, tmp_path: Path) -> None:
    target = tmp_path / "out.bin"
    RangeHandler.mode, RangeHandler.cut_at = "cut", 3000
    with pytest.raises(DownloadError) as caught:
        io_formats.resume_download(server, options={"output": str(target)}, verbose=False)
    assert caught.value.code in {"DOWNLOAD_LENGTH_MISMATCH", "DOWNLOAD_NETWORK_ERROR"}
    part = tmp_path / "out.bin.part"
    assert part.read_bytes() == PAYLOAD[:3000] and not target.exists()
    RangeHandler.mode = "normal"
    result = call(
        "resume_download",
        server,
        output=str(target),
        expected_sha256=hashlib.sha256(PAYLOAD).hexdigest(),
    )
    assert (
        target.read_bytes() == PAYLOAD
        and result["resumed_from"] == 3000
        and result["restarted"] is False
    )
    assert RangeHandler.requests[-1] == {"range": "bytes=3000-", "if_range": ETAG}


def test_resume_download_restarts_when_range_is_ignored_or_the_resource_changed(
    server: str, tmp_path: Path
) -> None:
    target = tmp_path / "out.bin"
    (tmp_path / "out.bin.part").write_bytes(b"stale-bytes")
    RangeHandler.mode = "ignore_range"
    result = call("resume_download", server, output=str(target))
    assert (
        target.read_bytes() == PAYLOAD
        and result["restarted"] is True
        and result["resumed_from"] == 0
    )
    # the validator saved with the part file no longer matches: the server answers 200 with the full new body
    other = tmp_path / "other.bin"
    RangeHandler.mode, RangeHandler.etag, RangeHandler.body = "normal", '"v2"', PAYLOAD[::-1]
    (tmp_path / "other.bin.part").write_bytes(PAYLOAD[:100])
    (tmp_path / "other.bin.part.meta").write_text(json.dumps({"etag": ETAG}))
    changed = call("resume_download", server, output=str(other))
    assert other.read_bytes() == PAYLOAD[::-1] and changed["restarted"] is True


def test_resume_download_416_means_complete_or_stale(server: str, tmp_path: Path) -> None:
    done = tmp_path / "done.bin"
    (tmp_path / "done.bin.part").write_bytes(PAYLOAD)
    result = call("resume_download", server, output=str(done))
    assert done.read_bytes() == PAYLOAD and result["resumed_from"] == len(PAYLOAD)
    stale = tmp_path / "stale.bin"
    (tmp_path / "stale.bin.part").write_bytes(PAYLOAD + b"extra")
    result = call("resume_download", server, output=str(stale))
    assert stale.read_bytes() == PAYLOAD and result["restarted"] is True


def test_resume_download_verifies_and_bounds(server: str, tmp_path: Path) -> None:
    target = tmp_path / "out.bin"
    with pytest.raises(DownloadError) as caught:
        io_formats.resume_download(
            server, options={"output": str(target), "expected_sha256": "0" * 64}, verbose=False
        )
    assert caught.value.code == "DOWNLOAD_CHECKSUM_MISMATCH"
    assert not target.exists() and not (tmp_path / "out.bin.part").exists()
    with pytest.raises(DownloadError) as caught:
        io_formats.resume_download(
            server, options={"output": str(target), "max_bytes": 1000}, verbose=False
        )
    assert caught.value.code == "DOWNLOAD_TOO_LARGE"
    target.write_bytes(b"existing")
    with pytest.raises(ConfigurationError):
        io_formats.resume_download(server, options={"output": str(target)}, verbose=False)
    assert target.read_bytes() == b"existing"


@pytest.mark.parametrize(
    "url", ["ftp://example.org/a", "file:///etc/passwd", "not a url", "http://"]
)
def test_resume_download_rejects_other_schemes(url: str, tmp_path: Path) -> None:
    with pytest.raises(InputFormatError):
        io_formats.resume_download(url, options={"output": str(tmp_path / "x")}, verbose=False)
    assert not (tmp_path / "x").exists()


def test_resume_download_needs_an_output_and_reports_connection_errors(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        io_formats.resume_download("http://127.0.0.1:9/x", verbose=False)
    with pytest.raises(DownloadError) as caught:
        io_formats.resume_download(
            "http://127.0.0.1:9/x",
            options={"output": str(tmp_path / "x"), "timeout": 2},
            verbose=False,
        )
    assert caught.value.code == "DOWNLOAD_NETWORK_ERROR"


# ---------------------------------------------------------------------------------------------
# invalid input and options
# ---------------------------------------------------------------------------------------------


def test_record_inputs_work_for_validation_and_indexing() -> None:
    assert call("validate_map", [["1", "a", 5], ["1", "a", 6]])["errors"] == [
        "line 2: duplicate variant ID 'a' (first on line 1)"
    ]
    assert call("index_map", [["1", "a", 5]])["entries"]["a"]["offset"] is None
    person = {"fid": "F", "iid": "1", "genotypes": [["A", "C"]]}
    assert "duplicate individual F:1" in call("validate_ped", [person, person])["errors"][0]
    record = dict(
        query_name="q",
        query_length=10,
        query_start=0,
        query_end=10,
        strand="+",
        target_name="t",
        target_length=20,
        target_start=0,
        target_end=10,
        matches=9,
        block_length=10,
        mapq=5,
    )
    assert call("validate_paf", [record])["valid"] is True
    assert call("validate_paf", [{"query_name": "q"}])["valid"] is False


@pytest.mark.parametrize("fmt", ["map", "ped", "paf", "newick", "nexus", "phylip"])
def test_read_rejects_unsupported_input_types(fmt: str) -> None:
    with pytest.raises(InputFormatError):
        getattr(io_formats, f"read_{fmt}")(12345, verbose=False)


@pytest.mark.parametrize(
    ("function", "text"),
    [
        ("read_map", "1 rs1 0\n"),
        ("read_map", "1 rs1 0 -5\n"),
        ("read_ped", "F I 0 0 1\n"),
        ("read_ped", "F I 0 0 x 1 A A\n"),
        ("read_paf", "q\t10\t0\t5\t+\tt\n"),
        ("read_newick", "((A,B);"),
        ("read_newick", "   "),
        ("read_nexus", "(A,B);"),
        ("read_phylip", ""),
    ],
)
def test_malformed_content_raises_an_input_format_error(function: str, text: str) -> None:
    with pytest.raises(InputFormatError) as caught:
        getattr(io_formats, function)(text, verbose=False)
    assert caught.value.code


def test_unknown_options_and_missing_files_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="unknown option"):
        io_formats.read_map(MAP, options={"nonsense": 1}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.filter_map(MAP, options={"min_cm": "high"}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.filter_map(MAP, options={"region": "chr1-100"}, verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.read_paf(str(tmp_path / "absent.paf"), verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.read_paf(tmp_path / "absent.paf", verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.convert_map(MAP, verbose=False)  # 'to' is required
    with pytest.raises(ConfigurationError):
        io_formats.convert_map(MAP, options={"to": "nexus"}, verbose=False)
    with pytest.raises(InputFormatError):
        io_formats.merge_map(MAP, verbose=False)  # merge needs a list
    with pytest.raises(ConfigurationError):
        io_formats.read_map(MAP, options={"max_records": 2}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_formats.read_map(MAP, options={"max_bytes": 10}, verbose=False)


def test_output_option_writes_files_and_refuses_to_overwrite(tmp_path: Path) -> None:
    target = tmp_path / "sub" / "out.map"
    result = call("sort_map", MAP, output=str(target))
    assert "text" not in result and result["output_path"] == str(target)
    assert target.read_text() == "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
    with pytest.raises(ConfigurationError):
        io_formats.sort_map(MAP, options={"output": str(target)}, verbose=False)
    call("sort_map", MAP, output=str(target), overwrite=True)
    src = tmp_path / "in.map"
    src.write_text(MAP)
    written = call("split_map", src, output_dir=str(tmp_path / "parts"))["parts"]
    assert Path(written[0]["output_path"]).read_text() == "1\trs1\t0\t100\n1\trs2\t0.5\t50\n"
    assert Path(written[2]["output_path"]).name == "part_0003_X.map"


# ---------------------------------------------------------------------------------------------
# independent cross-checks with Biopython parsers/writers (optional)
# ---------------------------------------------------------------------------------------------


def test_outputs_are_readable_by_biopython(tmp_path: Path) -> None:
    import io

    pytest.importorskip("Bio")
    from Bio import AlignIO, Phylo, SeqIO

    records = {"alpha": "ACGTACGTAC", "beta": "ACGTACGAAC", "gamma": "AC-TACGTNC"}
    for fmt, function in (("phylip-relaxed", "write_phylip"), ("nexus", "write_nexus")):
        for options in ({}, {"layout": "interleaved", "block_size": 4}):
            text = call(function, records, **options)["text"]
            alignment = AlignIO.read(io.StringIO(text), fmt)
            assert {r.id: str(r.seq) for r in alignment} == records
    tree = call("convert_newick", "((A:1,B:2)0.9:0.5,(C:1,D:1):0.5);", to="nexus")["text"]
    parsed = Phylo.read(io.StringIO(tree), "nexus")
    assert [t.name for t in parsed.get_terminals()] == [
        "A",
        "B",
        "C",
        "D",
    ] and parsed.total_branch_length() == 6.0
    small = tmp_path / "small.2bit"
    small.write_bytes(
        decode_payload(
            call("write_twobit", {"s1": "TCAGNNac", "s2": "ACGTACGTAC", "r": "ACGTRYKMacgtnn"})
        )
    )
    assert {r.id: str(r.seq) for r in SeqIO.parse(small, "twobit")} == {
        "s1": "TCAGNNac",
        "s2": "ACGTACGTAC",
        "r": "ACGTNNNNacgtnn",
    }
    big = tmp_path / "big.2bit"
    big.write_bytes(decode_payload(call("write_twobit", {"s1": "TCAGNNac"}, byte_order="big")))
    assert [str(r.seq) for r in SeqIO.parse(big, "twobit")] == ["TCAGNNac"]


def test_biopython_written_alignments_are_read_back() -> None:
    import io

    pytest.importorskip("Bio")
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = {"alpha": "ACGTACGTAC", "beta": "ACGTACGAAC", "gamma": "AC-TACGTNC"}
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq(s), id=n, annotations={"molecule_type": "DNA"}) for n, s in records.items()]
    )
    for fmt, function, options in (
        ("phylip-relaxed", "read_phylip", {}),
        ("phylip", "read_phylip", {"input_strict_names": True}),
        ("nexus", "read_nexus", {}),
    ):
        buffer = io.StringIO()
        AlignIO.write(alignment, buffer, fmt)
        rows = call(function, buffer.getvalue(), **options)["records"]
        assert {r["name"]: r["sequence"] for r in rows} == records
