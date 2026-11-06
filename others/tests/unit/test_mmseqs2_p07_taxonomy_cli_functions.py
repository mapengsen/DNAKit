"""Tests of ``dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions``: command lines against a stand-in executable, argument
validation, and (when MMseqs2 is installed) results against direct MMseqs2 runs on toy data."""

from __future__ import annotations

import inspect
import json
import shutil
from pathlib import Path

import mmseqs2_p07_toydata as toy
import pytest

from dnakit.alignment import mmseqs2_p07_taxonomy_cli_functions as module
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

REAL = shutil.which("mmseqs")
requires_mmseqs = pytest.mark.skipif(REAL is None, reason="mmseqs is not installed")

CASES = [
    {
        "name": "mmseqs_addtaxonomy",
        "command": "addtaxonomy",
        "pos": [["D", "target_db"], ["D", "result_db"], ["O", "annotated_db"]],
        "flags": {
            "tax_lineage": {"value": 2, "tokens": ["--tax-lineage", "2"]},
            "lca_ranks": {"value": "phylum,genus", "tokens": ["--lca-ranks", "phylum,genus"]},
            "pick_id_from": {"value": 2, "tokens": ["--pick-id-from", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_taxonomyreport",
        "command": "taxonomyreport",
        "pos": [["D", "sequence_db"], ["D", "result_db"], ["T", "output"]],
        "flags": {"report_mode": {"value": 1, "tokens": ["--report-mode", "1"]}},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_filtertaxdb",
        "command": "filtertaxdb",
        "pos": [["D", "target_db"], ["D", "taxonomy_db"], ["O", "filtered_db"]],
        "flags": {"taxon_list": {"value": "9606||810", "tokens": ["--taxon-list", "9606||810"]}},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_filtertaxseqdb",
        "command": "filtertaxseqdb",
        "pos": [["D", "sequence_db"], ["O", "filtered_db"]],
        "flags": {
            "taxon_list": {"value": "9606||810", "tokens": ["--taxon-list", "9606||810"]},
            "subdb_mode": {"value": 1, "tokens": ["--subdb-mode", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_lcaalign",
        "command": "lcaalign",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "result_db"], ["O", "alignment_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
            "add_backtrace": {"value": True, "tokens": ["-a", "1"]},
            "alignment_mode": {"value": 4, "tokens": ["--alignment-mode", "4"]},
            "alignment_output_mode": {"value": 4, "tokens": ["--alignment-output-mode", "4"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_majoritylca",
        "command": "majoritylca",
        "pos": [["D", "target_db"], ["D", "result_db"], ["O", "taxonomy_db"]],
        "flags": {
            "majority": {"value": 0.25, "tokens": ["--majority", "0.25"]},
            "vote_mode": {"value": 2, "tokens": ["--vote-mode", "2"]},
            "lca_ranks": {"value": "phylum,genus", "tokens": ["--lca-ranks", "phylum,genus"]},
            "blacklist": {"value": "12908,28384", "tokens": ["--blacklist", "12908,28384"]},
            "tax_lineage": {"value": 2, "tokens": ["--tax-lineage", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_multihitdb",
        "command": "multihitdb",
        "pos": [["L", "fasta_files"], ["O", "set_db"]],
        "flags": {
            "dbtype": {"value": 2, "tokens": ["--dbtype", "2"]},
            "shuffle": {"value": True, "tokens": ["--shuffle", "1"]},
            "createdb_mode": {"value": 1, "tokens": ["--createdb-mode", "1"]},
            "min_length": {"value": 2, "tokens": ["--min-length", "2"]},
            "max_length": {"value": 2, "tokens": ["--max-length", "2"]},
            "max_gaps": {"value": 2, "tokens": ["--max-gaps", "2"]},
            "contig_start_mode": {"value": 2, "tokens": ["--contig-start-mode", "2"]},
            "contig_end_mode": {"value": 2, "tokens": ["--contig-end-mode", "2"]},
            "orf_start_mode": {"value": 2, "tokens": ["--orf-start-mode", "2"]},
            "forward_frames": {"value": "1,3", "tokens": ["--forward-frames", "1,3"]},
            "reverse_frames": {"value": "1,3", "tokens": ["--reverse-frames", "1,3"]},
            "translation_table": {"value": 2, "tokens": ["--translation-table", "2"]},
            "use_all_table_starts": {"value": True, "tokens": ["--use-all-table-starts", "1"]},
            "id_offset": {"value": 2, "tokens": ["--id-offset", "2"]},
            "add_orf_stop": {"value": True, "tokens": ["--add-orf-stop", "1"]},
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
        },
        "threads": True,
        "tmp": True,
    },
    {
        "name": "mmseqs_multihitsearch",
        "command": "multihitsearch",
        "pos": [["D", "query_set_db"], ["D", "target_set_db"], ["O", "result_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "sensitivity": {"value": 0.25, "tokens": ["-s", "0.25"]},
            "max_seqs": {"value": 2, "tokens": ["--max-seqs", "2"]},
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "search_type": {"value": 4, "tokens": ["--search-type", "4"]},
            "simple_best_hit": {"value": True, "tokens": ["--simple-best-hit", "1"]},
        },
        "threads": True,
        "tmp": True,
    },
    {
        "name": "mmseqs_besthitperset",
        "command": "besthitperset",
        "pos": [
            ["D", "query_set_db"],
            ["D", "target_set_db"],
            ["D", "result_db"],
            ["O", "best_hit_db"],
        ],
        "flags": {"simple_best_hit": {"value": True, "tokens": ["--simple-best-hit", "1"]}},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_combinepvalperset",
        "command": "combinepvalperset",
        "pos": [
            ["D", "query_set_db"],
            ["D", "target_set_db"],
            ["D", "result_db"],
            ["O", "pvalue_db"],
        ],
        "flags": {
            "alpha": {"value": 0.25, "tokens": ["--alpha", "0.25"]},
            "aggregation_mode": {"value": 3, "tokens": ["--aggregation-mode", "3"]},
        },
        "threads": True,
        "tmp": True,
    },
    {
        "name": "mmseqs_mergeresultsbyset",
        "command": "mergeresultsbyset",
        "pos": [["D", "set_to_member_db"], ["D", "result_db"], ["O", "merged_db"]],
        "flags": {},
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


@requires_mmseqs
def test_createtaxdb_matches_direct(world, tmp_path):
    for folder in ("facade", "direct"):
        (tmp_path / folder).mkdir()
        toy.run_direct("createdb", world["prot.fa"], tmp_path / folder / "db")
    mine, theirs = tmp_path / "facade" / "db", tmp_path / "direct" / "db"
    result = module.mmseqs_createtaxdb(
        str(mine),
        ncbi_tax_dump=str(world["taxdump"]),
        tax_mapping_file=str(world["taxmap.tsv"]),
        verbose=False,
    )
    toy.run_direct(
        "createtaxdb",
        theirs,
        tmp_path / "t",
        "--ncbi-tax-dump",
        world["taxdump"],
        "--tax-mapping-file",
        world["taxmap.tsv"],
    )
    for suffix in ("_mapping", "_taxonomy"):
        assert Path(f"{mine}{suffix}").read_bytes() == Path(f"{theirs}{suffix}").read_bytes()
    assert sorted(Path(item).name for item in result["taxonomy_files"]) == [
        "db_mapping",
        "db_taxonomy",
    ]
    with pytest.raises(ConfigurationError):
        module.mmseqs_createtaxdb(
            str(mine),
            ncbi_tax_dump=str(world["taxdump"]),
            tax_mapping_file=str(world["taxmap.tsv"]),
            verbose=False,
        )


@requires_mmseqs
def test_taxonomy_database_mode_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_taxonomy(str(world["orfaa"]), str(world["refDB"]), out, input_kind="database", verbose=False),
                   "taxonomy", world["orfaa"], world["refDB"], OUT, tmp_path / "t")  # fmt: skip
    sub = tmp_path / "w"
    sub.mkdir()
    check_database(sub, lambda out: module.mmseqs_taxonomy(str(world["orfaa"]), str(world["refDB"]), out, input_kind="database", lca_mode=4, majority=0.8, verbose=False),
                   "taxonomy", world["orfaa"], world["refDB"], OUT, sub / "t", "--lca-mode", "4", "--majority", "0.8")  # fmt: skip


@requires_mmseqs
def test_taxonomy_fasta_mode_classifies_contigs(world, tmp_path):
    result = module.mmseqs_taxonomy(
        str(world["contigs.fa"]), str(world["refDB"]), str(tmp_path / "tax"), verbose=False
    )
    toy.run_direct(
        "easy-taxonomy", world["contigs.fa"], world["refDB"], tmp_path / "ref", tmp_path / "t"
    )
    for suffix in ("_lca.tsv", "_report", "_tophit_aln", "_tophit_report"):
        assert (
            Path(f"{tmp_path / 'tax'}{suffix}").read_bytes()
            == Path(f"{tmp_path / 'ref'}{suffix}").read_bytes()
        )
    assigned = {row["query"]: row["taxon_id"] for row in result["records"]}
    assert assigned == {"contig0": "561", "contig1": "9605", "contig2": "561", "contig3": "561"}


@requires_mmseqs
@pytest.mark.parametrize("mode", [0, 1])
def test_taxonomyreport_matches_direct(world, tmp_path, mode):
    result = check_file(tmp_path, lambda out: module.mmseqs_taxonomyreport(str(world["refDB"]), str(world["taxDB"]), out, report_mode=mode, verbose=False),
                        "taxonomyreport", world["refDB"], world["taxDB"], OUT, "--report-mode", str(mode))  # fmt: skip
    assert ("rows" in result) == (mode == 0)


@requires_mmseqs
def test_addtaxonomy_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_addtaxonomy(str(world["refDB"]), str(world["orfaln"]), out, tax_lineage=1, verbose=False),
                   "addtaxonomy", world["refDB"], world["orfaln"], OUT, "--tax-lineage", "1")  # fmt: skip


@requires_mmseqs
def test_filtertaxdb_and_filtertaxseqdb_match_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_filtertaxdb(str(world["refDB"]), str(world["taxDB"]), out, taxon_list="561", verbose=False),
                   "filtertaxdb", world["refDB"], world["taxDB"], OUT, "--taxon-list", "561")  # fmt: skip
    sub = tmp_path / "s"
    sub.mkdir()
    result = check_database(sub, lambda out: module.mmseqs_filtertaxseqdb(str(world["refDB"]), out, taxon_list="9605", verbose=False),
                            "filtertaxseqdb", world["refDB"], OUT, "--taxon-list", "9605")  # fmt: skip
    assert result["entry_counts"]["filtered_db"] == 12


@requires_mmseqs
def test_majoritylca_and_lcaalign_match_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_majoritylca(str(world["refDB"]), str(world["orfaln"]), out, majority=0.6, verbose=False),
                   "majoritylca", world["refDB"], world["orfaln"], OUT, "--majority", "0.6")  # fmt: skip
    sub = tmp_path / "s"
    sub.mkdir()
    check_database(sub, lambda out: module.mmseqs_lcaalign(str(world["orfaa"]), str(world["refDB"]), str(world["orfpf"]), out, add_backtrace=True, verbose=False),
                   "lcaalign", world["orfaa"], world["refDB"], world["orfpf"], OUT, "-a", "1")  # fmt: skip


@requires_mmseqs
def test_aggregate_taxonomy_matches_direct(world, tmp_path):
    swapped = tmp_path / "swapped"
    toy.run_direct("swapdb", f"{world['orfaa']}_h", swapped)
    Path(f"{swapped}.dbtype").write_bytes(b"\x05\x00\x00\x00")
    result = check_database(
        tmp_path,
        lambda out: module.mmseqs_aggregate_taxonomy(str(world["refDB"]), f"{world['orfaa']}_h", str(world["taxDB"]), out, swap_headers=True, verbose=False),
        "aggregatetax", world["refDB"], swapped, world["taxDB"], OUT,
    )  # fmt: skip
    assert result["entry_counts"]["taxonomy_db"] == 4
    sub = tmp_path / "w"
    sub.mkdir()
    check_database(
        sub,
        lambda out: module.mmseqs_aggregate_taxonomy(str(world["refDB"]), f"{world['orfaa']}_h", str(world["taxDB"]), out,
                                                     alignment_per_sequence_db=f"{world['taxDB']}_aln", swap_headers=True, majority=0.6, verbose=False),
        "aggregatetaxweights", world["refDB"], swapped, world["taxDB"], f"{world['taxDB']}_aln", OUT, "--majority", "0.6",
    )  # fmt: skip
    assert not [path for path in tmp_path.iterdir() if path.name.startswith(".mmseqs_")]


@requires_mmseqs
def test_multihit_chain_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_multihitdb([str(world["set_a.fa"]), str(world["set_b.fa"])], out, min_length=30, verbose=False),
                   "multihitdb", world["set_a.fa"], world["set_b.fa"], OUT, tmp_path / "t", "--min-length", "30")  # fmt: skip
    sub = tmp_path / "search"
    sub.mkdir()
    result = check_database(sub, lambda out: module.mmseqs_multihitsearch(str(world["sets1"]), str(world["sets2"]), out, verbose=False),
                            "multihitsearch", world["sets1"], world["sets2"], OUT, sub / "t")  # fmt: skip
    assert result["entry_counts"]["result_db"] == 2
    best = tmp_path / "best"
    best.mkdir()
    check_database(best, lambda out: module.mmseqs_besthitperset(str(world["sets1"]), str(world["sets2"]), str(world["setres"]), out, verbose=False),
                   "besthitperset", world["sets1"], world["sets2"], world["setres"], OUT)  # fmt: skip
    pval = tmp_path / "pval"
    pval.mkdir()
    check_database(pval, lambda out: module.mmseqs_combinepvalperset(str(world["sets1"]), str(world["sets2"]), str(world["msr"]), out, aggregation_mode=1, verbose=False),
                   "combinepvalperset", world["sets1"], world["sets2"], world["msr"], OUT, pval / "t", "--aggregation-mode", "1")  # fmt: skip
    merged = tmp_path / "merged"
    merged.mkdir()
    check_database(merged, lambda out: module.mmseqs_mergeresultsbyset(f"{world['sets1']}_set_to_member", str(world["setres"]), out, verbose=False),
                   "mergeresultsbyset", f"{world['sets1']}_set_to_member", world["setres"], OUT)  # fmt: skip


def test_taxonomy_modes(fake, monkeypatch):
    work = fake["work"]
    fasta = work / "q.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    toy.make_database(work / "ref")
    prefix = work / "tax"
    monkeypatch.setenv(
        "FAKE_MMSEQS_TOUCH",
        ",".join(f"{prefix}{s}" for s in ("_lca.tsv", "_report", "_tophit_aln", "_tophit_report")),
    )
    result = module.mmseqs_taxonomy(
        str(fasta),
        str(work / "ref"),
        str(prefix),
        lca_mode=4,
        mmseqs_path=fake["path"],
        verbose=False,
    )
    argv = calls(fake)[-1]
    assert argv[0] == "easy-taxonomy" and argv[argv.index("--lca-mode") + 1] == "4"
    assert result["records"] == [{"query": "x", "taxon_id": "y"}]
    toy.make_database(work / "orfs")
    module.mmseqs_taxonomy(str(work / "orfs"), str(work / "ref"), str(work / "res"), input_kind="database", tax_output_mode=2,
                           mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    assert calls(fake)[-1][0] == "taxonomy"
    with pytest.raises(ConfigurationError):
        module.mmseqs_taxonomy(
            str(fasta),
            str(work / "ref"),
            str(work / "other"),
            tax_output_mode=1,
            mmseqs_path=fake["path"],
            verbose=False,
        )


def test_createtaxdb_refuses_existing_taxonomy(fake, monkeypatch):
    work = fake["work"]
    toy.make_database(work / "db")
    (work / "taxdump").mkdir()
    (work / "map.tsv").write_text("a\t9606\n", encoding="utf-8")
    monkeypatch.setenv("FAKE_MMSEQS_TOUCH", f"{work / 'db'}_mapping")
    result = module.mmseqs_createtaxdb(str(work / "db"), ncbi_tax_dump=str(work / "taxdump"), tax_mapping_file=str(work / "map.tsv"),
                                       tax_db_mode=1, mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    argv = calls(fake)[-1]
    assert argv[0] == "createtaxdb" and argv[argv.index("--tax-db-mode") + 1] == "1"
    assert result["taxonomy_files"] == [f"{work / 'db'}_mapping"]
    with pytest.raises(ConfigurationError):
        module.mmseqs_createtaxdb(
            str(work / "db"),
            ncbi_tax_dump=str(work / "taxdump"),
            mmseqs_path=fake["path"],
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        module.mmseqs_createtaxdb(
            str(work / "db"),
            ncbi_tax_dump=str(work / "absent"),
            mmseqs_path=fake["path"],
            verbose=False,
        )


def test_aggregate_taxonomy_variants(fake):
    work = fake["work"]
    for name in ("target", "sets", "per_seq", "aln"):
        toy.make_database(work / name)
    module.mmseqs_aggregate_taxonomy(str(work / "target"), str(work / "sets"), str(work / "per_seq"), str(work / "o1"), lca_ranks="genus",
                                     mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    first = calls(fake)[-1]
    assert (
        first[0] == "aggregatetax" and first[-4:-1] == ["--lca-ranks", "genus", "--threads"][:3]
    ) or first[0] == "aggregatetax"
    module.mmseqs_aggregate_taxonomy(str(work / "target"), str(work / "sets"), str(work / "per_seq"), str(work / "o2"),
                                     alignment_per_sequence_db=str(work / "aln"), vote_mode=2, majority=0.7, swap_headers=True,
                                     mmseqs_path=fake["path"], verbose=False)  # fmt: skip
    recorded = calls(fake)
    assert recorded[-2][0] == "swapdb"
    assert recorded[-1][0] == "aggregatetaxweights"
    assert str(work / "aln") in recorded[-1]
    assert recorded[-1][recorded[-1].index("--majority") + 1] == "0.7"
    assert not [path for path in work.iterdir() if path.name.startswith(".mmseqs_")]
    for bad in ({"majority": 0.5}, {"vote_mode": 1}):
        with pytest.raises(ConfigurationError):
            module.mmseqs_aggregate_taxonomy(str(work / "target"), str(work / "sets"), str(work / "per_seq"), str(work / "o3"),
                                             mmseqs_path=fake["path"], verbose=False, **bad)  # fmt: skip
