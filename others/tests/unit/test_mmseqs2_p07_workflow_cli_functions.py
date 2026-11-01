"""Tests of ``dnakit.alignment.mmseqs2_p07_workflow_cli_functions``: command lines against a stand-in executable, argument
validation, and (when MMseqs2 is installed) results against direct MMseqs2 runs on toy data."""

# ruff: noqa: E501
from __future__ import annotations

import inspect
import json
import shutil
from pathlib import Path

import mmseqs2_p07_toydata as toy
import pytest

from dnakit.alignment import mmseqs2_p07_workflow_cli_functions as module
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

REAL = shutil.which("mmseqs")
requires_mmseqs = pytest.mark.skipif(REAL is None, reason="mmseqs is not installed")

CASES = [
    {
        "name": "mmseqs_createdb",
        "command": "createdb",
        "pos": [["L", "fasta_files"], ["O", "database"]],
        "flags": {
            "dbtype": {"value": 2, "tokens": ["--dbtype", "2"]},
            "shuffle": {"value": True, "tokens": ["--shuffle", "1"]},
            "createdb_mode": {"value": 1, "tokens": ["--createdb-mode", "1"]},
            "write_lookup": {"value": True, "tokens": ["--write-lookup", "1"]},
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
            "id_offset": {"value": 2, "tokens": ["--id-offset", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_map",
        "command": "map",
        "pos": [["D", "query_db"], ["D", "target_db"], ["O", "alignment_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "add_backtrace": {"value": True, "tokens": ["-a", "1"]},
            "sensitivity": {"value": 0.25, "tokens": ["-s", "0.25"]},
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "max_seqs": {"value": 2, "tokens": ["--max-seqs", "2"]},
            "search_type": {"value": 4, "tokens": ["--search-type", "4"]},
        },
        "threads": True,
        "tmp": True,
    },
    {
        "name": "mmseqs_clusterupdate",
        "command": "clusterupdate",
        "pos": [
            ["D", "old_sequence_db"],
            ["D", "new_sequence_db"],
            ["D", "old_cluster_db"],
            ["O", "mapped_sequence_db"],
            ["O", "new_cluster_db"],
        ],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "cluster_mode": {"value": 3, "tokens": ["--cluster-mode", "3"]},
            "sensitivity": {"value": 0.25, "tokens": ["-s", "0.25"]},
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "alignment_mode": {"value": 4, "tokens": ["--alignment-mode", "4"]},
        },
        "threads": True,
        "tmp": True,
    },
    {
        "name": "mmseqs_convertalis",
        "command": "convertalis",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "alignment_db"], ["T", "output"]],
        "flags": {
            "format_mode": {"value": 4, "tokens": ["--format-mode", "4"]},
            "format_output": {
                "value": "query,target,bits",
                "tokens": ["--format-output", "query,target,bits"],
            },
            "search_type": {"value": 4, "tokens": ["--search-type", "4"]},
            "translation_table": {"value": 2, "tokens": ["--translation-table", "2"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_createtsv",
        "command": "createtsv",
        "pos": [["D", "query_db"], ["Dopt", "target_db"], ["D", "result_db"], ["T", "output"]],
        "flags": {
            "first_seq_as_repr": {"value": True, "tokens": ["--first-seq-as-repr", "1"]},
            "target_column": {"value": 2, "tokens": ["--target-column", "2"]},
            "full_header": {"value": True, "tokens": ["--full-header", "1"]},
            "idx_seq_src": {"value": 2, "tokens": ["--idx-seq-src", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_convert2fasta",
        "command": "convert2fasta",
        "pos": [["D", "sequence_db"], ["T", "output"]],
        "flags": {"use_header_file": {"value": True, "tokens": ["--use-header-file", "1"]}},
        "threads": False,
        "tmp": False,
    },
    {
        "name": "mmseqs_clust",
        "command": "clust",
        "pos": [["D", "sequence_db"], ["D", "result_db"], ["O", "cluster_db"]],
        "flags": {
            "cluster_mode": {"value": 3, "tokens": ["--cluster-mode", "3"]},
            "max_iterations": {"value": 2, "tokens": ["--max-iterations", "2"]},
            "similarity_type": {"value": 2, "tokens": ["--similarity-type", "2"]},
            "set_mode": {"value": True, "tokens": ["--set-mode", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_clusthash",
        "command": "clusthash",
        "pos": [["D", "sequence_db"], ["O", "alignment_db"]],
        "flags": {
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_pickconsensusrep",
        "command": "pickconsensusrep",
        "pos": [["D", "sequence_db"], ["D", "cluster_db"], ["O", "new_cluster_db"]],
        "flags": {},
        "threads": False,
        "tmp": True,
    },
    {
        "name": "mmseqs_result2repseq",
        "command": "result2repseq",
        "pos": [["D", "sequence_db"], ["D", "result_db"], ["O", "representative_db"]],
        "flags": {},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_summarizealis",
        "command": "summarizealis",
        "pos": [["D", "alignment_db"], ["O", "summary_db"]],
        "flags": {},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_filterresult",
        "command": "filterresult",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "result_db"], ["O", "filtered_db"]],
        "flags": {
            "max_seq_id": {"value": 0.25, "tokens": ["--max-seq-id", "0.25"]},
            "qid": {"value": "0.5,0.9", "tokens": ["--qid", "0.5,0.9"]},
            "qsc": {"value": 0.25, "tokens": ["--qsc", "0.25"]},
            "cov": {"value": 0.25, "tokens": ["--cov", "0.25"]},
            "diff": {"value": 2, "tokens": ["--diff", "2"]},
            "filter_min_enable": {"value": 2, "tokens": ["--filter-min-enable", "2"]},
            "allow_deletion": {"value": True, "tokens": ["--allow-deletion", "1"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_extractalignedregion",
        "command": "extractalignedregion",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "alignment_db"], ["O", "region_db"]],
        "flags": {"extract_mode": {"value": 2, "tokens": ["--extract-mode", "2"]}},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_gff2db",
        "command": "gff2db",
        "pos": [["L", "gff_files"], ["D", "sequence_db"], ["O", "feature_db"]],
        "flags": {
            "gff_type": {"value": "CDS,gene", "tokens": ["--gff-type", "CDS,gene"]},
            "id_offset": {"value": 2, "tokens": ["--id-offset", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
]


@pytest.fixture
def fake(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    executable = toy.make_fake_mmseqs(bin_dir)
    log = tmp_path / "calls.jsonl"
    monkeypatch.setenv("FAKE_MMSEQS_LOG", str(log))
    work = tmp_path / "work"
    work.mkdir()
    return {"path": str(executable), "log": log, "work": work}


def calls(fake):
    return [json.loads(line) for line in fake["log"].read_text(encoding="utf-8").splitlines()]


def prepare(case, work, *, options):
    """Create the inputs of a case and return (positional arguments, expected positional paths)."""
    args, expected = [], []
    for code, name in case["pos"]:
        path = work / name
        if code in ("D", "Dopt"):
            toy.make_database(path)
        elif code in ("F", "L"):
            path.write_text(">a\nACGT\n", encoding="utf-8")
        if code == "L":
            args.append([str(path)])
        elif code != "Dopt":
            args.append(str(path))
        expected.append(str(path))
    return args, expected


@pytest.mark.parametrize("case", CASES, ids=lambda item: item["name"])
@pytest.mark.parametrize("with_options", [False, True], ids=["defaults", "all_options"])
def test_command_line(case, with_options, fake):
    function = getattr(module, case["name"])
    work = fake["work"]
    args, expected = prepare(case, work, options=with_options)
    keywords = {}
    optional = [name for code, name in case["pos"] if code == "Dopt"]
    for name in optional:
        keywords[name] = str(work / name)
    tokens = []
    if with_options:
        for name, item in case["flags"].items():
            keywords[name] = item["value"]
            tokens += item["tokens"]
    result = function(*args, mmseqs_path=fake["path"], verbose=False, **keywords)
    assert result["status"] == "completed"
    assert result["backend"] == "mmseqs2"
    assert result["version"] == "fake-mmseqs-1"
    [argv] = calls(fake)
    assert argv[0] == case["command"]
    positional = argv[1 : 1 + len(case["pos"]) + (1 if case["tmp"] else 0)]
    if case["tmp"]:
        tmp_folder = positional.pop()
        assert Path(tmp_folder).name.startswith(".mmseqs_tmp_")
        assert not Path(tmp_folder).exists()
    assert positional == expected
    rest = argv[1 + len(positional) + (1 if case["tmp"] else 0) :]
    assert rest == [*tokens, *(["--threads", "1"] if case["threads"] else [])]
    assert result["command"][0] == "mmseqs"
    for code, name in case["pos"]:
        if code in ("O", "T"):
            assert Path(expected[[n for _, n in case["pos"]].index(name)]).exists()


@pytest.mark.parametrize("case", CASES, ids=lambda item: item["name"])
def test_outputs_are_never_overwritten(case, fake):
    if not any(code in ("O", "T") for code, _ in case["pos"]):
        pytest.skip("no new output")
    function = getattr(module, case["name"])
    args, _ = prepare(case, fake["work"], options=False)
    function(*args, mmseqs_path=fake["path"], verbose=False)
    with pytest.raises(ConfigurationError):
        function(*args, mmseqs_path=fake["path"], verbose=False)


@pytest.mark.parametrize("case", CASES, ids=lambda item: item["name"])
def test_missing_input_database_is_refused(case, fake):
    inputs = [name for code, name in case["pos"] if code == "D"]
    if not inputs:
        pytest.skip("no database input")
    function = getattr(module, case["name"])
    args, _ = prepare(case, fake["work"], options=False)
    index = [name for _, name in case["pos"] if _ != "Dopt"].index(inputs[0])
    args[index] = str(fake["work"] / "absent")
    with pytest.raises(ConfigurationError):
        function(*args, mmseqs_path=fake["path"], verbose=False)
    assert not fake["log"].exists()


@pytest.mark.parametrize("case", CASES, ids=lambda item: item["name"])
def test_invalid_option_values_are_refused(case, fake):
    function = getattr(module, case["name"])
    args, _ = prepare(case, fake["work"], options=False)
    for name, item in case["flags"].items():
        wrong = "not a value" if not isinstance(item["value"], str) else 12345
        with pytest.raises(ConfigurationError):
            function(*args, mmseqs_path=fake["path"], verbose=False, **{name: wrong})
    assert not fake["log"].exists()


def test_failed_run_raises_with_log_tail(fake, monkeypatch):
    case = CASES[0]
    function = getattr(module, case["name"])
    args, _ = prepare(case, fake["work"], options=False)
    monkeypatch.setenv("FAKE_MMSEQS_FAIL", "1")
    with pytest.raises(BackendExecutionError) as caught:
        function(*args, mmseqs_path=fake["path"], verbose=False)
    assert "simulated failure" in str(caught.value.context)


def test_missing_program_has_install_hint(fake, tmp_path):
    case = CASES[0]
    function = getattr(module, case["name"])
    args, _ = prepare(case, fake["work"], options=False)
    with pytest.raises(BackendUnavailableError):
        function(*args, mmseqs_path=str(tmp_path / "no-such-mmseqs"), verbose=False)


def test_every_public_function_has_a_keyword_only_verbose():
    assert module.__all__ == sorted(module.__all__)
    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True
        assert inspect.getdoc(getattr(module, name))


OUT = object()


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if REAL is None:
        pytest.skip("mmseqs is not installed")
    return toy.build_world(tmp_path_factory.mktemp("mmseqs_world"))


def sub(tmp_path, name):
    folder = tmp_path / name
    folder.mkdir()
    return folder


def check_database(tmp_path, facade, *direct, threads=True, **kwargs):
    """Run the facade and the plain mmseqs command; the result databases must be byte-identical."""
    out, ref = tmp_path / "facade", tmp_path / "direct"
    result = facade(str(out), **kwargs)
    toy.run_direct(*[ref if item is OUT else item for item in direct], threads=threads)
    assert toy.database_files(out) == toy.database_files(ref)
    assert toy.database_files(out)["entries"] is not None
    return result


def check_file(tmp_path, facade, *direct, threads=True, **kwargs):
    """Run the facade and the plain mmseqs command; the result files must be byte-identical."""
    out, ref = tmp_path / "facade.out", tmp_path / "direct.out"
    result = facade(str(out), **kwargs)
    toy.run_direct(*[ref if item is OUT else item for item in direct], threads=threads)
    assert out.read_bytes() == ref.read_bytes()
    return result


def tmp_for(tmp_path):
    return tmp_path / "scratch"


@requires_mmseqs
def test_createdb_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_createdb([str(world["prot.fa"])], out, verbose=False),
                   "createdb", world["prot.fa"], OUT)  # fmt: skip


@requires_mmseqs
def test_createdb_options_match_direct(world, tmp_path):
    check_database(
        tmp_path,
        lambda out: module.mmseqs_createdb(str(world["contigs.fa"]), out, dbtype=2, shuffle=False, verbose=False),
        "createdb", world["contigs.fa"], OUT, "--dbtype", "2", "--shuffle", "0",
    )  # fmt: skip


@requires_mmseqs
def test_createdb_several_files_and_empty_input(world, tmp_path):
    result = check_database(
        tmp_path,
        lambda out: module.mmseqs_createdb([str(world["set_a.fa"]), str(world["set_b.fa"])], out, verbose=False),
        "createdb", world["set_a.fa"], world["set_b.fa"], OUT,
    )  # fmt: skip
    assert result["entry_counts"]["database"] == 4
    empty = tmp_path / "empty.fa"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(BackendExecutionError):
        module.mmseqs_createdb(str(empty), str(tmp_path / "none"), verbose=False)


@requires_mmseqs
def test_map_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_map(str(world["protDB"]), str(world["protDB"]), out, verbose=False),
                   "map", world["protDB"], world["protDB"], OUT, tmp_path / "tmp")  # fmt: skip


@requires_mmseqs
def test_clusterupdate_matches_direct_and_consumes_the_old_clustering(world, tmp_path):
    for name in ("facade", "direct"):
        toy.run_direct("clust", world["protDB"], world["aln"], tmp_path / f"old_{name}")
    result = module.mmseqs_clusterupdate(
        str(world["protDB"]),
        str(world["protDB"]),
        str(tmp_path / "old_facade"),
        str(tmp_path / "f_map"),
        str(tmp_path / "f_clu"),
        verbose=False,
    )
    toy.run_direct(
        "clusterupdate",
        world["protDB"],
        world["protDB"],
        tmp_path / "old_direct",
        tmp_path / "d_map",
        tmp_path / "d_clu",
        tmp_path / "t",
    )
    assert toy.database_files(tmp_path / "f_clu") == toy.database_files(tmp_path / "d_clu")
    assert toy.database_files(tmp_path / "f_map") == toy.database_files(tmp_path / "d_map")
    assert result["old_cluster_db_removed"] is True
    assert not (tmp_path / "old_direct.dbtype").exists()


@requires_mmseqs
@pytest.mark.parametrize(
    "options", [{}, {"format_output": "query,target,bits"}, {"format_mode": 1}, {"format_mode": 4}]
)
def test_convertalis_matches_direct(world, tmp_path, options):
    direct = []
    if "format_output" in options:
        direct = ["--format-output", options["format_output"]]
    if "format_mode" in options:
        direct = ["--format-mode", str(options["format_mode"])]
    out, ref = tmp_path / "facade.tsv", tmp_path / "direct.tsv"
    result = module.mmseqs_convertalis(
        str(world["protDB"]),
        str(world["protDB"]),
        str(world["aln"]),
        str(out),
        verbose=False,
        **options,
    )
    toy.run_direct("convertalis", world["protDB"], world["protDB"], world["aln"], ref, *direct)
    assert out.read_bytes() == ref.read_bytes()
    if options.get("format_mode") in (None, 0, 4):
        assert result["records"]
        first = result["records"][0]
        assert set(first) == set(
            (
                options.get("format_output")
                or "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
            ).split(",")
        )
    else:
        assert "records" not in result


@requires_mmseqs
@pytest.mark.parametrize("with_target", [True, False], ids=["with_target", "clusters_only"])
def test_createtsv_matches_direct(world, tmp_path, with_target):
    out, ref = tmp_path / "facade.tsv", tmp_path / "direct.tsv"
    if with_target:
        result = module.mmseqs_createtsv(
            str(world["protDB"]),
            str(world["clu"]),
            str(out),
            target_db=str(world["protDB"]),
            verbose=False,
        )
        toy.run_direct("createtsv", world["protDB"], world["protDB"], world["clu"], ref)
    else:
        result = module.mmseqs_createtsv(
            str(world["protDB"]), str(world["clu"]), str(out), target_column=0, verbose=False
        )
        toy.run_direct("createtsv", world["protDB"], world["clu"], ref, "--target-column", "0")
    assert out.read_bytes() == ref.read_bytes()
    assert result["rows"]


@requires_mmseqs
def test_convert2fasta_matches_direct(world, tmp_path):
    result = check_file(tmp_path, lambda out: module.mmseqs_convert2fasta(str(world["protDB"]), out, verbose=False),
                        "convert2fasta", world["protDB"], OUT, threads=False)  # fmt: skip
    assert result["sequence_count"] == 24


@requires_mmseqs
def test_clust_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_clust(str(world["protDB"]), str(world["aln"]), out, verbose=False),
                   "clust", world["protDB"], world["aln"], OUT)  # fmt: skip
    check_database(sub(tmp_path, "connected"),
                   lambda out: module.mmseqs_clust(str(world["protDB"]), str(world["aln"]), out, cluster_mode=1, verbose=False),
                   "clust", world["protDB"], world["aln"], OUT, "--cluster-mode", "1")  # fmt: skip


@requires_mmseqs
def test_clusthash_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_clusthash(str(world["protDB"]), out, verbose=False),
                   "clusthash", world["protDB"], OUT)  # fmt: skip


@requires_mmseqs
def test_pickconsensusrep_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_pickconsensusrep(str(world["protDB"]), str(world["clu"]), out, verbose=False),
                   "pickconsensusrep", world["protDB"], world["clu"], OUT, tmp_path / "t", threads=False)  # fmt: skip


@requires_mmseqs
def test_result2repseq_matches_direct(world, tmp_path):
    result = check_database(tmp_path, lambda out: module.mmseqs_result2repseq(str(world["protDB"]), str(world["clu"]), out, verbose=False),
                            "result2repseq", world["protDB"], world["clu"], OUT)  # fmt: skip
    assert result["entry_counts"]["representative_db"] == 6


@requires_mmseqs
def test_summarizealis_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_summarizealis(str(world["aln"]), out, verbose=False),
                   "summarizealis", world["aln"], OUT)  # fmt: skip


@requires_mmseqs
def test_filterresult_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_filterresult(str(world["protDB"]), str(world["protDB"]), str(world["aln"]), out, qid="0.9", verbose=False),
                   "filterresult", world["protDB"], world["protDB"], world["aln"], OUT, "--qid", "0.9")  # fmt: skip


@requires_mmseqs
@pytest.mark.parametrize("mode", [1, 2])
def test_extractalignedregion_matches_direct(world, tmp_path, mode):
    check_database(tmp_path, lambda out: module.mmseqs_extractalignedregion(str(world["protDB"]), str(world["protDB"]), str(world["aln"]), out, extract_mode=mode, verbose=False),
                   "extractalignedregion", world["protDB"], world["protDB"], world["aln"], OUT, "--extract-mode", str(mode))  # fmt: skip


@requires_mmseqs
def test_gff2db_matches_direct(world, tmp_path):
    result = check_database(tmp_path, lambda out: module.mmseqs_gff2db([str(world["features.gff3"])], str(world["ctgDB"]), out, gff_type="CDS", verbose=False),
                            "gff2db", world["features.gff3"], world["ctgDB"], OUT, "--gff-type", "CDS")  # fmt: skip
    assert result["entry_counts"]["feature_db"] == 8


@requires_mmseqs
def test_maskbygff_translates_sequence_names_to_keys(world, tmp_path):
    lookup = {}
    for line in Path(f"{world['ctgDB']}.lookup").read_text(encoding="utf-8").splitlines():
        key, name, _ = line.split("\t")
        lookup[name] = key
    keyed = tmp_path / "keyed.gff3"
    keyed.write_text(
        "".join(
            (lookup[row.split("\t")[0]] + row[row.index("\t") :])
            if not row.startswith("#")
            else row
            for row in world["features.gff3"].read_text(encoding="utf-8").splitlines(keepends=True)
        ),
        encoding="utf-8",
    )
    check_database(tmp_path, lambda out: module.mmseqs_maskbygff(str(world["features.gff3"]), str(world["ctgDB"]), out, verbose=False),
                   "maskbygff", keyed, world["ctgDB"], OUT)  # fmt: skip
    plain = tmp_path / "plain"
    with pytest.raises(BackendExecutionError):
        module.mmseqs_maskbygff(
            str(world["features.gff3"]),
            str(world["ctgDB"]),
            str(plain),
            names_as_ids=False,
            verbose=False,
        )
    assert not list(tmp_path.glob(".mmseqs_gff_*"))


@requires_mmseqs
def test_linsearch_fasta_and_database_match_direct(world, tmp_path):
    out, ref = tmp_path / "lin.tsv", tmp_path / "ref.tsv"
    result = module.mmseqs_linsearch(
        str(world["contigs.fa"]), str(world["contigs.fa"]), str(out), search_type=3, verbose=False
    )
    toy.run_direct(
        "easy-linsearch",
        world["contigs.fa"],
        world["contigs.fa"],
        ref,
        tmp_path / "t",
        "--search-type",
        "3",
    )
    assert out.read_bytes() == ref.read_bytes()
    assert {row["query"] for row in result["records"]} == {
        "contig0",
        "contig1",
        "contig2",
        "contig3",
    }
    db_out = tmp_path / "lindb"
    with pytest.raises(ConfigurationError):
        module.mmseqs_linsearch(
            str(world["ntDB"]),
            str(world["ntDB"]),
            str(db_out),
            input_kind="database",
            search_type=3,
            verbose=False,
        )
    module.mmseqs_createlinindex(str(world["ntDB"]), search_type=3, verbose=False)
    module.mmseqs_linsearch(
        str(world["ntDB"]),
        str(world["ntDB"]),
        str(db_out),
        input_kind="database",
        search_type=3,
        verbose=False,
    )
    toy.run_direct(
        "linsearch",
        world["ntDB"],
        world["ntDB"],
        tmp_path / "lindb_ref",
        tmp_path / "t2",
        "--search-type",
        "3",
    )
    assert toy.database_files(db_out) == toy.database_files(tmp_path / "lindb_ref")


@requires_mmseqs
def test_linclust_fasta_and_database_match_direct(world, tmp_path):
    result = module.mmseqs_linclust(str(world["prot.fa"]), str(tmp_path / "lc"), verbose=False)
    toy.run_direct("easy-linclust", world["prot.fa"], tmp_path / "ref", tmp_path / "t")
    for suffix in ("_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta"):
        assert (
            Path(f"{tmp_path / 'lc'}{suffix}").read_bytes()
            == Path(f"{tmp_path / 'ref'}{suffix}").read_bytes()
        )
    assert result["sequence_count"] == 24
    assert result["cluster_count"] == len(result["clusters"]) == 12
    check_database(sub(tmp_path, "lcdb"),
                   lambda out: module.mmseqs_linclust(str(world["protDB"]), out, input_kind="database", min_seq_id=0.5, verbose=False),
                   "linclust", world["protDB"], OUT, tmp_path / "t3", "--min-seq-id", "0.5")  # fmt: skip


@requires_mmseqs
def test_rbh_fasta_and_database_match_direct(world, tmp_path):
    out, ref = tmp_path / "rbh.tsv", tmp_path / "ref.tsv"
    result = module.mmseqs_rbh(
        str(world["prot.fa"]), str(world["prot.fa"]), str(out), verbose=False
    )
    toy.run_direct("easy-rbh", world["prot.fa"], world["prot.fa"], ref, tmp_path / "t")
    assert out.read_bytes() == ref.read_bytes()
    assert result["pair_count"] == 24
    check_database(sub(tmp_path, "rbhdb"),
                   lambda dbout: module.mmseqs_rbh(str(world["protDB"]), str(world["protDB"]), dbout, input_kind="database", verbose=False),
                   "rbh", world["protDB"], world["protDB"], OUT, tmp_path / "t4")  # fmt: skip


@requires_mmseqs
def test_createlinindex_matches_direct(world, tmp_path):
    copy_a, copy_b = tmp_path / "a", tmp_path / "b"
    for copy in (copy_a, copy_b):
        toy.run_direct("createdb", world["contigs.fa"], copy, "--dbtype", "2")
    result = module.mmseqs_createlinindex(str(copy_a), search_type=3, verbose=False)
    toy.run_direct("createlinindex", copy_b, tmp_path / "t", "--search-type", "3")
    assert result["entry_counts"]["kmer_index_db"] > 0
    first, second = Path(f"{copy_a}.linidx"), Path(f"{copy_b}.linidx")
    assert first.stat().st_size == second.stat().st_size
    # the tail of the index file holds uninitialised padding that differs between runs
    half = first.stat().st_size // 2
    assert first.read_bytes()[:half] == second.read_bytes()[:half]
    assert (
        first.with_suffix(".linidx.index").read_bytes()
        == second.with_suffix(".linidx.index").read_bytes()
    )
    with pytest.raises(ConfigurationError):
        module.mmseqs_createlinindex(str(copy_a), verbose=False)


@requires_mmseqs
def test_result2msa_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_result2msa(str(world["protDB"]), str(world["protDB"]), str(world["aln"]), out, verbose=False),
                   "result2msa", world["protDB"], world["protDB"], world["aln"], OUT)  # fmt: skip
    sub = tmp_path / "dna"
    sub.mkdir()
    check_database(sub, lambda out: module.mmseqs_result2msa(str(world["ntDB"]), str(world["ntDB"]), str(world["ntaln"]), out, dna=True, verbose=False),
                   "result2dnamsa", world["ntDB"], world["ntDB"], world["ntaln"], OUT)  # fmt: skip


def test_linsearch_modes_use_the_documented_commands(fake):
    work = fake["work"]
    query, target = work / "q.fa", work / "t.fa"
    query.write_text(">a\nACGT\n", encoding="utf-8")
    target.write_text(">b\nACGT\n", encoding="utf-8")
    out = work / "out.tsv"
    result = module.mmseqs_linsearch(str(query), str(target), str(out), search_type=3, e_value=0.5,
                                     format_mode=4, mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    argv = calls(fake)[-1]
    assert argv[0] == "easy-linsearch"
    assert argv[1:3] == [str(query), str(target)]
    assert argv[3] == str(out)
    assert argv[argv.index("--search-type") + 1] == "3"
    assert argv[argv.index("--format-mode") + 1] == "4"
    assert result["outputs"]["output"] == str(out)
    toy.make_database(work / "qdb")
    toy.make_database(work / "tdb")
    with pytest.raises(ConfigurationError):
        module.mmseqs_linsearch(
            str(work / "qdb"),
            str(work / "tdb"),
            str(work / "res0"),
            input_kind="database",
            mmseqs_path=fake["path"],
            verbose=False,
        )
    toy.make_database(work / "tdb.linidx")
    module.mmseqs_linsearch(str(work / "qdb"), str(work / "tdb"), str(work / "res"), input_kind="database",
                            mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    assert calls(fake)[-1][0] == "linsearch"
    with pytest.raises(ConfigurationError):
        module.mmseqs_linsearch(str(work / "qdb"), str(work / "tdb"), str(work / "res2"), input_kind="database",
                                format_mode=0, mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    with pytest.raises(ConfigurationError):
        module.mmseqs_linsearch(
            str(query),
            str(target),
            str(work / "x"),
            input_kind="rows",
            mmseqs_path=fake["path"],
            verbose=False,
        )


def test_linclust_fasta_collects_clusters(fake, monkeypatch):
    work = fake["work"]
    fasta = work / "s.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    prefix = work / "lc"
    monkeypatch.setenv(
        "FAKE_MMSEQS_TOUCH",
        ",".join(f"{prefix}{s}" for s in ("_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta")),
    )
    result = module.mmseqs_linclust(
        str(fasta),
        str(prefix),
        min_seq_id=0.9,
        coverage=0.8,
        mmseqs_path=fake["path"],
        verbose=False,
    )
    argv = calls(fake)[-1]
    assert argv[0] == "easy-linclust"
    assert argv[1] == str(fasta) and argv[2] == str(prefix)
    assert result["clusters"] == {"x": ["y"]}
    assert result["cluster_count"] == 1
    with pytest.raises(ConfigurationError):
        module.mmseqs_linclust(str(fasta), str(prefix), mmseqs_path=fake["path"], verbose=False)


def test_rbh_modes(fake):
    work = fake["work"]
    fasta = work / "s.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    module.mmseqs_rbh(
        str(fasta),
        str(fasta),
        str(work / "o.tsv"),
        sensitivity=7.5,
        mmseqs_path=fake["path"],
        verbose=False,
    )
    argv = calls(fake)[-1]
    assert argv[0] == "easy-rbh" and argv[argv.index("-s") + 1] == "7.5"
    toy.make_database(work / "d")
    module.mmseqs_rbh(
        str(work / "d"),
        str(work / "d"),
        str(work / "res"),
        input_kind="database",
        mmseqs_path=fake["path"],
        verbose=False,
    )
    assert calls(fake)[-1][0] == "rbh"


def test_createlinindex_command(fake, monkeypatch):
    work = fake["work"]
    toy.make_database(work / "db")
    monkeypatch.setenv("FAKE_MMSEQS_TOUCH", f"{work / 'db'}.linidx.dbtype")
    result = module.mmseqs_createlinindex(
        str(work / "db"), search_type=3, kmer_per_seq=40, mmseqs_path=fake["path"], verbose=False
    )
    argv = calls(fake)[-1]
    assert argv[0] == "createlinindex" and argv[1] == str(work / "db")
    assert argv[argv.index("--kmer-per-seq") + 1] == "40"
    assert result["outputs"]["kmer_index_db"] == f"{work / 'db'}.linidx"


def test_result2msa_dna_accepts_only_skip_query(fake):
    work = fake["work"]
    for name in ("q", "t", "r"):
        toy.make_database(work / name)
    module.mmseqs_result2msa(str(work / "q"), str(work / "t"), str(work / "r"), str(work / "m"), dna=True, skip_query=True,
                             mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    argv = calls(fake)[-1]
    assert argv[0] == "result2dnamsa" and argv[argv.index("--skip-query") + 1] == "1"
    with pytest.raises(ConfigurationError):
        module.mmseqs_result2msa(str(work / "q"), str(work / "t"), str(work / "r"), str(work / "m2"), dna=True, cov=0.5,
                                 mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    module.mmseqs_result2msa(str(work / "q"), str(work / "t"), str(work / "r"), str(work / "m3"), msa_format_mode=5,
                             mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    assert calls(fake)[-1][0] == "result2msa"


def test_maskbygff_rewrites_a_temporary_copy_only(fake):
    work = fake["work"]
    toy.make_database(work / "genome")
    Path(f"{work / 'genome'}.lookup").write_text("0\tchrA\t0\n1\tchrB\t0\n", encoding="utf-8")
    gff = work / "f.gff3"
    original = (
        "##gff-version 3\nchrA\ts\tCDS\t1\t3\t.\t+\t0\tID=1\nchrZ\ts\tCDS\t1\t3\t.\t+\t0\tID=2\n"
    )
    gff.write_text(original, encoding="utf-8")
    module.mmseqs_maskbygff(
        str(gff),
        str(work / "genome"),
        str(work / "masked"),
        gff_type="CDS",
        mmseqs_path=fake["path"],
        verbose=False,
    )
    argv = calls(fake)[-1]
    assert argv[0] == "maskbygff"
    assert argv[1] != str(gff)
    assert Path(argv[1]).name == "features.gff3"
    assert gff.read_text(encoding="utf-8") == original
    assert not list(work.glob(".mmseqs_gff_*"))
    module.mmseqs_maskbygff(
        str(gff),
        str(work / "genome"),
        str(work / "masked2"),
        names_as_ids=False,
        mmseqs_path=fake["path"],
        verbose=False,
    )
    assert calls(fake)[-1][1] == str(gff)
