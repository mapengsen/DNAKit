"""Regression tests for the io_formats fixes: gzip codecs, validation messages and a few option checks."""

from __future__ import annotations

import base64
import gzip
import hashlib
import inspect
import io
import tarfile
from pathlib import Path

import pytest

from dnakit import io_formats
from dnakit.exceptions import ConfigurationError, InputFormatError

FORMATS = ["map", "newick", "nexus", "paf", "ped", "phylip", "twobit", "zarr"]
ALL_BYTES = bytes(range(256)) * 4


def result_of(function: str, data: object, **options: object) -> dict:
    envelope = getattr(io_formats, function)(data, options=options or None)
    assert envelope["status"] == "ok"
    return envelope["result"]


def restored(result: dict) -> bytes:
    if "data_base64" in result:
        return base64.b64decode(result["data_base64"])
    return result["text"].encode("utf-8")


def roundtrip(fmt: str, data: object) -> bytes:
    packed = result_of(f"compress_{fmt}", data)
    return restored(result_of(f"decompress_{fmt}", packed["data"]))


# ---------------------------------------------------------------------------------------------
# compress / decompress
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
def test_text_and_binary_round_trip_is_bit_exact(fmt: str) -> None:
    text = "1\trs1\t0\t1000\nnaïve – Ω\n"
    assert roundtrip(fmt, text) == text.encode("utf-8")
    assert roundtrip(fmt, text.encode("utf-8")) == text.encode("utf-8")
    assert roundtrip(fmt, ALL_BYTES) == ALL_BYTES
    assert roundtrip(fmt, bytearray(ALL_BYTES)) == ALL_BYTES
    assert roundtrip(fmt, b"") == b""


@pytest.mark.parametrize("fmt", FORMATS)
def test_compress_result_fields_and_deterministic_gzip(fmt: str) -> None:
    payload = ("ACGT" * 500).encode()
    first = result_of(f"compress_{fmt}", payload)
    second = result_of(f"compress_{fmt}", payload)
    assert first == second
    assert first["encoding"] == "gzip+base64"
    packed = base64.b64decode(first["data"])
    assert packed[4:8] == b"\x00\x00\x00\x00"  # no timestamp in the gzip header
    assert gzip.decompress(packed) == payload
    assert first["original_bytes"] == len(payload)
    assert first["compressed_bytes"] == len(packed) < len(payload)
    assert first["sha256_original"] == hashlib.sha256(payload).hexdigest()
    assert "output_path" not in first
    plain = result_of(f"decompress_{fmt}", first["data"])
    assert plain["original_bytes"] == len(payload) and plain["sha256"] == first["sha256_original"]
    assert plain["text"] == payload.decode()


def test_level_option_controls_compression() -> None:
    payload = bytes(range(256)) * 200 + b"ACGT" * 5000
    fast = result_of("compress_map", payload, level=1)
    best = result_of("compress_map", payload, level=9)
    assert fast["compressed_bytes"] >= best["compressed_bytes"]
    assert fast["data"] != best["data"]
    assert gzip.decompress(base64.b64decode(best["data"])) == payload
    for bad in (0, 10, "x", 1.5):
        with pytest.raises(ConfigurationError):
            io_formats.compress_map("x", options={"level": bad})


def test_real_twobit_file_survives_the_round_trip(tmp_path: Path) -> None:
    written = result_of("write_twobit", {"chr1": "ACGTNNacgtACGT" * 9, "chr2": "GGGCCCAAATTT"})
    blob = base64.b64decode(written["data_base64"])
    path = tmp_path / "genome.2bit"
    path.write_bytes(blob)
    assert roundtrip("twobit", blob) == blob
    assert roundtrip("twobit", path) == blob
    assert roundtrip("twobit", str(path)) == blob
    packed = result_of("compress_twobit", path)
    assert packed["original_bytes"] == len(blob)
    # the decompressed 2bit file is still readable by the format engine
    again = tmp_path / "again.2bit"
    result_of("decompress_twobit", packed["data"], output=str(again))
    assert again.read_bytes() == blob
    assert result_of("read_twobit", str(again))["record_count"] == 2


def test_file_inputs_and_missing_or_wrong_paths(tmp_path: Path) -> None:
    source = tmp_path / "x.map"
    source.write_bytes(b"1\trs1\t0\t5\n\xff\xfe")
    assert roundtrip("map", source) == source.read_bytes()
    assert roundtrip("map", str(source)) == source.read_bytes()
    with pytest.raises(InputFormatError) as missing:
        io_formats.compress_map(tmp_path / "absent.map")
    assert missing.value.code == "IO_FORMATS_FILE_NOT_FOUND"
    with pytest.raises(InputFormatError):
        io_formats.compress_map(str(tmp_path / "absent.map"))
    with pytest.raises(InputFormatError) as directory:
        io_formats.compress_map(tmp_path)
    assert directory.value.hint
    with pytest.raises(InputFormatError):
        io_formats.compress_map({"a": 1})
    with pytest.raises(InputFormatError):
        io_formats.compress_map(12)


def test_options_are_validated_and_size_is_capped(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        io_formats.compress_map("x", options={"nope": 1})
    with pytest.raises(ConfigurationError) as too_large:
        io_formats.compress_map("x" * 100, options={"max_bytes": 10})
    assert too_large.value.code == "IO_FORMATS_TOO_LARGE"
    big = tmp_path / "big.map"
    big.write_bytes(b"x" * 100)
    with pytest.raises(ConfigurationError):
        io_formats.compress_map(big, options={"max_bytes": 10})
    bomb = result_of("compress_map", b"\0" * 100_000)["data"]
    with pytest.raises(ConfigurationError):
        io_formats.decompress_map(bomb, options={"max_bytes": 1000})
    assert result_of("decompress_map", bomb, max_bytes=100_000)["original_bytes"] == 100_000


def test_files_are_written_only_with_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    result_of("compress_map", "abc")
    assert list(tmp_path.iterdir()) == []
    target = tmp_path / "out.gz"
    packed = result_of("compress_map", "abc", output=str(target))
    assert packed["output_path"] == str(target)
    assert target.read_bytes() == base64.b64decode(packed["data"])
    assert gzip.decompress(target.read_bytes()) == b"abc"
    with pytest.raises(ConfigurationError) as exists:
        io_formats.compress_map("abc", options={"output": str(target)})
    assert exists.value.code == "IO_FORMATS_EXISTS"
    result_of("compress_map", "different", output=str(target), overwrite=True)
    assert gzip.decompress(target.read_bytes()) == b"different"
    # a .gz path is accepted by decompress, raw gzip bytes too
    assert result_of("decompress_map", target)["text"] == "different"
    assert result_of("decompress_map", str(target))["text"] == "different"
    assert result_of("decompress_map", target.read_bytes())["text"] == "different"
    restored_file = tmp_path / "restored.bin"
    info = result_of("decompress_map", target, output=str(restored_file))
    assert info["output_path"] == str(restored_file) and "text" not in info
    assert restored_file.read_bytes() == b"different"
    with pytest.raises(ConfigurationError):
        io_formats.decompress_map(target, options={"output": str(restored_file)})


def test_binary_decompress_returns_base64_and_text_returns_text() -> None:
    binary = result_of("decompress_map", result_of("compress_map", b"\xff\x00\x80")["data"])
    assert "text" not in binary and base64.b64decode(binary["data_base64"]) == b"\xff\x00\x80"
    textual = result_of("decompress_map", result_of("compress_map", "ok\n")["data"])
    assert textual["text"] == "ok\n" and "data_base64" not in textual


@pytest.mark.parametrize("fmt", FORMATS)
def test_invalid_payloads_raise_instead_of_returning_ok(fmt: str, tmp_path: Path) -> None:
    good = base64.b64decode(result_of(f"compress_{fmt}", "hello world" * 20)["data"])
    bad_inputs = [
        "!!! not base64 !!!",
        base64.b64encode(b"plain text, not gzip").decode(),
        b"plain text, not gzip",
        base64.b64encode(good[: len(good) // 2]).decode(),  # truncated
        good[:-6],  # truncated raw gzip
        base64.b64encode(good[:12] + b"\xff" * 20).decode(),  # corrupt body
        "",
    ]
    for bad in bad_inputs:
        with pytest.raises(InputFormatError) as caught:
            getattr(io_formats, f"decompress_{fmt}")(bad)
        assert caught.value.code and caught.value.hint
    junk = tmp_path / "junk.gz"
    junk.write_bytes(b"not gzip at all")
    with pytest.raises(InputFormatError):
        getattr(io_formats, f"decompress_{fmt}")(junk)
    with pytest.raises(InputFormatError):
        getattr(io_formats, f"decompress_{fmt}")(tmp_path / "missing.gz")
    with pytest.raises(InputFormatError):
        getattr(io_formats, f"decompress_{fmt}")(["x"])


def test_the_functions_no_longer_use_the_shared_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    import dnakit._candidate_runtime as runtime

    def boom(*args: object, **kwargs: object) -> None:
        raise AssertionError("run_candidate_operation must not be used")

    monkeypatch.setattr(runtime, "run_candidate_operation", boom)
    for fmt in FORMATS:
        assert roundtrip(fmt, "abc") == b"abc"
        assert "run_candidate_operation" not in inspect.getsource(getattr(io_formats, f"compress_{fmt}"))
        assert "run_candidate_operation" not in inspect.getsource(getattr(io_formats, f"decompress_{fmt}"))


# ---------------------------------------------------------------------------------------------
# zarr directories
# ---------------------------------------------------------------------------------------------


def make_tree(root: Path) -> None:
    (root / "group" / "deep").mkdir(parents=True)
    (root / "empty_dir").mkdir()
    (root / "zarr.json").write_text('{"zarr_format": 3}\n')
    (root / "group" / "chunk").write_bytes(ALL_BYTES)
    (root / "group" / "deep" / "naïve.bin").write_bytes(b"\x00\x01\x02")


def same_tree(left: Path, right: Path) -> None:
    names = sorted(p.relative_to(left) for p in left.rglob("*"))
    assert names == sorted(p.relative_to(right) for p in right.rglob("*"))
    for name in names:
        if (left / name).is_file():
            assert (left / name).read_bytes() == (right / name).read_bytes()
        else:
            assert (right / name).is_dir()


def test_directory_round_trip_is_exact_and_reproducible(tmp_path: Path) -> None:
    source = tmp_path / "store.zarr"
    make_tree(source)
    first = result_of("compress_zarr", source)
    again = result_of("compress_zarr", str(source))
    assert first == again and first["archive"] == "tar"
    assert first["original_bytes"] > 0
    target = tmp_path / "restored.zarr"
    out = result_of("decompress_zarr", first["data"], output=str(target))
    assert out["output_path"] == str(target) and out["file_count"] == 3
    same_tree(source, target)
    # without 'output' the tar archive is returned and nothing is written
    inline = result_of("decompress_zarr", first["data"])
    assert not (tmp_path / "other").exists() and inline["archive"] == "tar"
    with tarfile.open(fileobj=io.BytesIO(base64.b64decode(inline["data_base64"]))) as archive:
        assert {m.name for m in archive.getmembers()} == {
            "empty_dir",
            "group",
            "group/chunk",
            "group/deep",
            "group/deep/naïve.bin",
            "zarr.json",
        }
    # a non-empty output directory is refused unless overwrite=True
    with pytest.raises(ConfigurationError) as exists:
        io_formats.decompress_zarr(first["data"], options={"output": str(target)})
    assert exists.value.code == "IO_FORMATS_EXISTS"
    result_of("decompress_zarr", first["data"], output=str(target), overwrite=True)
    # the same directory is refused by the other formats, a missing one is an error
    with pytest.raises(InputFormatError):
        io_formats.compress_map(source)
    with pytest.raises(InputFormatError):
        io_formats.compress_zarr(tmp_path / "missing.zarr")


def test_real_zarr_store_survives_the_round_trip(tmp_path: Path) -> None:
    pytest.importorskip("zarr")
    pytest.importorskip("numpy")
    store = tmp_path / "real.zarr"
    result_of("write_zarr", {"a": [[1, 2, 3], [4, 5, 6]], "b": [0.5, 1.5]}, output=str(store))
    packed = result_of("compress_zarr", store)
    target = tmp_path / "copy.zarr"
    result_of("decompress_zarr", packed["data"], output=str(target))
    same_tree(store, target)
    assert result_of("read_zarr", str(target))["arrays"].keys() == {"a", "b"}


def test_directory_size_cap_and_links(tmp_path: Path) -> None:
    source = tmp_path / "s"
    source.mkdir()
    (source / "f").write_bytes(b"x" * 5000)
    with pytest.raises(ConfigurationError):
        io_formats.compress_zarr(source, options={"max_bytes": 1000})
    (source / "link").symlink_to("f")
    with pytest.raises(InputFormatError) as link:
        io_formats.compress_zarr(source)
    assert link.value.code == "IO_FORMATS_UNSAFE_ARCHIVE"


def tar_gz(entries: list[tarfile.TarInfo], payloads: dict[str, bytes] | None = None) -> str:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for info in entries:
            body = (payloads or {}).get(info.name)
            if body is not None:
                info.size = len(body)
                archive.addfile(info, io.BytesIO(body))
            else:
                archive.addfile(info)
    return base64.b64encode(gzip.compress(buffer.getvalue())).decode()


def member(name: str, kind: bytes = tarfile.REGTYPE, link: str = "") -> tarfile.TarInfo:
    info = tarfile.TarInfo(name)
    info.type, info.linkname = kind, link
    return info


@pytest.mark.parametrize(
    "unsafe",
    [
        member("/etc/evil.txt"),
        member("../escape.txt"),
        member("a/../../escape.txt"),
        member("link", tarfile.SYMTYPE, "/etc/passwd"),
        member("link", tarfile.SYMTYPE, "../outside"),
        member("link", tarfile.SYMTYPE, "inside"),
        member("hard", tarfile.LNKTYPE, "/etc/passwd"),
        member("dev", tarfile.CHRTYPE),
        member("fifo", tarfile.FIFOTYPE),
    ],
    ids=lambda info: f"{info.name}-{info.type!r}",
)
def test_unsafe_archive_members_are_refused_before_anything_is_written(
    unsafe: tarfile.TarInfo, tmp_path: Path
) -> None:
    ok = member("fine.txt")
    bodies = {"fine.txt": b"data"}
    if unsafe.type == tarfile.REGTYPE:
        bodies[unsafe.name] = b"x"
    payload = tar_gz([ok, unsafe], bodies)
    target = tmp_path / "out"
    with pytest.raises(InputFormatError) as caught:
        io_formats.decompress_zarr(payload, options={"output": str(target)})
    assert caught.value.code == "IO_FORMATS_UNSAFE_ARCHIVE"
    assert not target.exists()
    assert not (tmp_path / "escape.txt").exists()


def test_archive_structure_errors(tmp_path: Path) -> None:
    clash = tar_gz([member("a"), member("a/b")], {"a": b"1", "a/b": b"2"})
    with pytest.raises(InputFormatError):
        io_formats.decompress_zarr(clash, options={"output": str(tmp_path / "o1")})
    twice = tar_gz([member("a"), member("a")], {"a": b"1"})
    with pytest.raises(InputFormatError):
        io_formats.decompress_zarr(twice, options={"output": str(tmp_path / "o2")})
    # data that is not a tar archive is handled like the other formats (text or bytes, or a file)
    plain = base64.b64encode(gzip.compress(b"just some text, not a tar archive")).decode()
    assert result_of("decompress_zarr", plain)["text"] == "just some text, not a tar archive"
    written = tmp_path / "plain.txt"
    info = result_of("decompress_zarr", plain, output=str(written))
    assert written.read_bytes() == b"just some text, not a tar archive" and "archive" not in info
    file_target = tmp_path / "file"
    file_target.write_text("x")
    good = tar_gz([member("a")], {"a": b"1"})
    with pytest.raises(ConfigurationError):
        io_formats.decompress_zarr(good, options={"output": str(file_target)})


# ---------------------------------------------------------------------------------------------
# validate_map / validate_ped / validate_paf
# ---------------------------------------------------------------------------------------------


def errors_of(function: str, text: object, **options: object) -> list[str]:
    result = result_of(function, text, **options)
    assert result["valid"] is False
    return result["errors"]


def test_validate_map_reports_non_numeric_fields_with_line_field_and_value() -> None:
    errors = errors_of("validate_map", "1\trs1\tx\t100\n1\trs2\t0\tabc\n1\trs3\t0\t-5\n1\trs4\tnan\t7\n")
    assert len(errors) == 4
    assert errors[0].startswith("line 1:") and "'cM'" in errors[0] and "'x'" in errors[0]
    assert errors[1].startswith("line 2:") and "'bp'" in errors[1] and "'abc'" in errors[1]
    assert errors[2].startswith("line 3:") and "'-5'" in errors[2]
    assert errors[3].startswith("line 4:") and "'nan'" in errors[3]
    assert not any("could not convert" in e or "Error" in e for e in errors)


def test_validate_map_reports_every_problem_of_a_line() -> None:
    errors = errors_of("validate_map", "1\trs1\tx\tabc\n")
    assert len(errors) == 2
    assert "'cM'" in errors[0] and "'bp'" in errors[1]
    assert errors_of("validate_map", "1\trs1\n")[0].startswith("line 1: expected 4 columns")
    assert result_of("validate_map", "1\trs1\t0.5\t100\n")["valid"] is True


def test_validate_map_with_record_input_does_not_leak_python_messages() -> None:
    errors = errors_of("validate_map", [{"chrom": "1", "id": "rs1", "cm": "x", "pos": 5}])
    assert "'cm'" in errors[0] and "'x'" in errors[0]
    assert "could not convert" not in errors[0] and "ValueError" not in errors[0]
    assert result_of("validate_map", [{"chrom": "1", "id": "rs1", "cm": 0.5, "pos": 5.0}])["valid"] is True
    with pytest.raises(InputFormatError) as caught:
        io_formats.read_map("1\trs1\tx\t5\n")
    assert "could not convert" not in str(caught.value) and "'cM'" in str(caught.value)


def test_validate_ped_reports_non_numeric_sex_and_phenotype() -> None:
    errors = errors_of("validate_ped", "f1 i1 0 0 x 1 A A\nf1 i2 0 0 1 case A A\nf1 i3 0 0 1.5x 2.5 A A\n")
    assert len(errors) == 3
    assert errors[0].startswith("line 1:") and "'sex'" in errors[0] and "'x'" in errors[0]
    assert errors[1].startswith("line 2:") and "'phenotype'" in errors[1] and "'case'" in errors[1]
    assert errors[2].startswith("line 3:") and "'1.5x'" in errors[2]
    assert not any("could not convert" in e or "is not an integer" == e for e in errors)


def test_validate_ped_reports_every_problem_of_a_row() -> None:
    errors = errors_of("validate_ped", "f1 i1 0 0 x y A A\nf1 i2 0 0 q 1 A A C C\nf1 i3 0 0 1 1 A\n")
    line1 = [e for e in errors if e.startswith("line 1:")]
    assert len(line1) == 2 and "'sex'" in line1[0] and "'phenotype'" in line1[1]
    line2 = [e for e in errors if e.startswith("line 2:")]
    assert len(line2) == 2 and any("variants" in e for e in line2) and any("'sex'" in e for e in line2)
    assert any(e.startswith("line 3:") and "6 + 2k" in e for e in errors)
    errors = errors_of("validate_ped", [{"fid": "f", "iid": "i", "sex": "x", "pheno": "y", "genotypes": []}])
    assert "could not convert" not in errors[0] and "'sex'" in errors[0]


def test_validate_paf_reports_every_problem_of_a_line() -> None:
    row = "q\tx\ty\t10\t*\tt\t20\t0\t10\t9\t10\t300\tbad\n"
    errors = errors_of("validate_paf", row)
    assert all(e.startswith("line 1:") for e in errors)
    joined = " | ".join(errors)
    for expected in ("'query_length'", "'query_start'", "'strand'", "'mapq'", "'bad'"):
        assert expected in joined
    assert len(errors) == 5
    # consistent numbers: all the semantic problems are still listed together
    semantic = errors_of("validate_paf", "q\t10\t0\t20\t*\tt\t10\t0\t5\t9\t5\t300\n")
    assert len(semantic) == 4  # query end, strand, matches, mapq
    assert result_of("validate_paf", "q\t10\t0\t5\t+\tt\t10\t0\t5\t5\t5\t1\ttp:A:P\n")["valid"] is True
    two_lines = errors_of("validate_paf", "q\t10\t0\t5\t*\tt\t10\t0\t5\t5\t5\t999\nq\t10\t0\t5\t+\tt\t10\t0\t5\t5\t5\t1\ta\n")
    assert [e.split(":")[0] for e in two_lines] == ["line 1", "line 1", "line 2"]


# ---------------------------------------------------------------------------------------------
# split_phylip, convert_nexus, index_map, write_phylip, docstrings
# ---------------------------------------------------------------------------------------------

PHYLIP = " 2 4\nA ACGT\nB ACGA\n"
NEXUS = """#NEXUS
BEGIN DATA;
 DIMENSIONS NTAX=2 NCHAR=4;
 FORMAT DATATYPE=DNA MISSING=? GAP=-;
 MATRIX
 A ACGT
 B ACGA
 ;
END;
BEGIN TREES;
 TREE t1 = (A:1,B:2);
END;
"""


def test_split_phylip_documents_no_tree_option_and_rejects_it_clearly() -> None:
    assert "trees_per_part" not in io_formats.split_phylip.__doc__
    assert "trees_per_part" in io_formats.split_nexus.__doc__
    with pytest.raises(ConfigurationError) as caught:
        io_formats.split_phylip(PHYLIP, options={"trees_per_part": 1})
    assert "no trees" in str(caught.value)
    assert result_of("split_phylip", PHYLIP, site_size=2)["part_count"] == 2
    with pytest.raises(ConfigurationError):
        io_formats.split_phylip(PHYLIP, options={"nonsense": 1})


@pytest.mark.parametrize("target", ["twobit", "fasta", "phylip"])
def test_convert_nexus_with_trees_converts_the_matrix_and_warns(target: str) -> None:
    envelope = io_formats.convert_nexus(NEXUS, options={"to": target})
    assert envelope["status"] == "ok"
    result = envelope["result"]
    assert result["record_count"] == 2
    assert len(result["warnings"]) == 1 and "tree" in result["warnings"][0]
    assert envelope["warnings"] == result["warnings"]
    if target == "twobit":
        assert io_formats.read_twobit(restored(result))["result"]["record_count"] == 2
    clean = io_formats.convert_nexus(NEXUS.split("BEGIN TREES")[0], options={"to": target})
    assert "warnings" not in clean["result"] and clean["warnings"] == []


def test_index_map_reports_duplicate_ids() -> None:
    envelope = io_formats.index_map("1\trs1\t0\t10\n1\trs2\t0\t20\n1\trs1\t0\t30\n2\trs2\t0\t5\n")
    result = envelope["result"]
    assert result["record_count"] == 4
    assert result["duplicate_ids"] == ["rs1", "rs2"]
    assert result["entries"]["rs1"]["line"] == 3 and len(result["entries"]) == 2
    assert len(result["warnings"]) == 1 and "rs1" in result["warnings"][0]
    assert envelope["warnings"] == result["warnings"]
    clean = io_formats.index_map("1\trs1\t0\t10\n1\trs2\t0\t20\n")
    assert clean["result"]["duplicate_ids"] == [] and "warnings" not in clean["result"]
    assert clean["warnings"] == []


def test_write_phylip_checks_the_alphabet_of_the_datatype() -> None:
    with pytest.raises(InputFormatError) as caught:
        io_formats.write_phylip({"a": "ACGT", "b": "AC1X"})
    assert "'1'" in str(caught.value) and "'X'" in str(caught.value) and "'b'" in str(caught.value)
    assert caught.value.code == "IO_FORMATS_BAD_ALPHABET" and caught.value.hint
    with pytest.raises(InputFormatError):
        io_formats.write_phylip({"a": "AC>T"})
    # nucleotide IUPAC codes, gaps and missing characters (any case) are accepted
    ok = result_of("write_phylip", {"a": "ACGTNRYacgt-?.", "b": "ACGTNRYacgt-?."})
    assert ok["record_count"] == 2
    assert result_of("write_phylip", {"a": "ACGU"}, datatype="rna")["record_count"] == 1
    assert result_of("write_phylip", {"a": "MKVLX*", "b": "MKV--?"}, datatype="protein")["record_count"] == 2
    assert result_of("write_phylip", {"a": "12 34"}, datatype="any")["record_count"] == 1
    with pytest.raises(InputFormatError):
        io_formats.write_phylip({"a": "MKVLZ1"}, options={"datatype": "protein"})
    # protein text read from PHYLIP/FASTA input is checked against the datatype too
    with pytest.raises(InputFormatError):
        io_formats.write_phylip(">p\nMKVLEE\n")
    assert result_of("write_phylip", ">p\nMKVLEE\n", datatype="protein")["record_count"] == 1
    assert result_of("write_phylip", PHYLIP)["record_count"] == 2


def test_zarr_docstrings_match_the_dependency_hint() -> None:
    expected = "pip install 'zarr>=3' numpy"
    for name in dir(io_formats):
        if name.endswith("_zarr") and not name.startswith("_"):
            doc = getattr(io_formats, name).__doc__ or ""
            if "requires=" in doc:
                assert f'requires="{expected}"' in doc, name
    from dnakit.io_formats import _impl_zarr

    assert expected in inspect.getsource(_impl_zarr)
