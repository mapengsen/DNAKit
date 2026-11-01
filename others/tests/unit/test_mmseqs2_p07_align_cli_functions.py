"""Tests of ``dnakit.alignment.mmseqs2_p07_align_cli_functions``: command lines against a stand-in executable, argument
validation, and (when MMseqs2 is installed) results against direct MMseqs2 runs on toy data."""

# ruff: noqa: E501
from __future__ import annotations

import inspect
import json
import shutil
from pathlib import Path

import mmseqs2_p07_toydata as toy
import pytest

from dnakit.alignment import mmseqs2_p07_align_cli_functions as module
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

REAL = shutil.which("mmseqs")
requires_mmseqs = pytest.mark.skipif(REAL is None, reason="mmseqs is not installed")

CASES = [
    {
        "name": "mmseqs_prefilter",
        "command": "prefilter",
        "pos": [["D", "query_db"], ["D", "target_db"], ["O", "prefilter_db"]],
        "flags": {
            "sensitivity": {"value": 0.25, "tokens": ["-s", "0.25"]},
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "max_seqs": {"value": 2, "tokens": ["--max-seqs", "2"]},
            "min_ungapped_score": {"value": 2, "tokens": ["--min-ungapped-score", "2"]},
            "comp_bias_corr": {"value": True, "tokens": ["--comp-bias-corr", "1"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "spaced_kmer_mode": {"value": True, "tokens": ["--spaced-kmer-mode", "1"]},
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "taxon_list": {"value": "9606||810", "tokens": ["--taxon-list", "9606||810"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_ungappedprefilter",
        "command": "ungappedprefilter",
        "pos": [["D", "query_db"], ["D", "target_db"], ["O", "prefilter_db"]],
        "flags": {
            "max_seqs": {"value": 2, "tokens": ["--max-seqs", "2"]},
            "min_ungapped_score": {"value": 2, "tokens": ["--min-ungapped-score", "2"]},
            "comp_bias_corr": {"value": True, "tokens": ["--comp-bias-corr", "1"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "taxon_list": {"value": "9606||810", "tokens": ["--taxon-list", "9606||810"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_gappedprefilter",
        "command": "gappedprefilter",
        "pos": [["D", "query_db"], ["D", "target_db"], ["O", "prefilter_db"]],
        "flags": {
            "max_seqs": {"value": 2, "tokens": ["--max-seqs", "2"]},
            "min_ungapped_score": {"value": 2, "tokens": ["--min-ungapped-score", "2"]},
            "comp_bias_corr": {"value": True, "tokens": ["--comp-bias-corr", "1"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
            "taxon_list": {"value": "9606||810", "tokens": ["--taxon-list", "9606||810"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_kmermatcher",
        "command": "kmermatcher",
        "pos": [["D", "sequence_db"], ["O", "prefilter_db"]],
        "flags": {
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "kmer_per_seq": {"value": 2, "tokens": ["--kmer-per-seq", "2"]},
            "spaced_kmer_mode": {"value": True, "tokens": ["--spaced-kmer-mode", "1"]},
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
            "adjust_kmer_len": {"value": True, "tokens": ["--adjust-kmer-len", "1"]},
            "hash_shift": {"value": 2, "tokens": ["--hash-shift", "2"]},
            "ignore_multi_kmer": {"value": True, "tokens": ["--ignore-multi-kmer", "1"]},
            "linclust_version": {"value": 2, "tokens": ["--linclust-version", "2"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_kmersearch",
        "command": "kmersearch",
        "pos": [["D", "query_db"], ["D", "kmer_index_db"], ["O", "prefilter_db"]],
        "flags": {
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "kmer_per_seq": {"value": 2, "tokens": ["--kmer-per-seq", "2"]},
            "pick_n_sim_kmer": {"value": 2, "tokens": ["--pick-n-sim-kmer", "2"]},
            "result_direction": {"value": 1, "tokens": ["--result-direction", "1"]},
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_alignall",
        "command": "alignall",
        "pos": [["D", "sequence_db"], ["D", "result_db"], ["O", "alignment_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "add_backtrace": {"value": True, "tokens": ["-a", "1"]},
            "alignment_mode": {"value": 4, "tokens": ["--alignment-mode", "4"]},
            "comp_bias_corr": {"value": True, "tokens": ["--comp-bias-corr", "1"]},
            "score_bias": {"value": 0.25, "tokens": ["--score-bias", "0.25"]},
            "zdrop": {"value": 2, "tokens": ["--zdrop", "2"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_transitivealign",
        "command": "transitivealign",
        "pos": [["D", "sequence_db"], ["D", "alignment_db"], ["O", "transitive_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "add_backtrace": {"value": True, "tokens": ["-a", "1"]},
            "alignment_mode": {"value": 4, "tokens": ["--alignment-mode", "4"]},
            "alignment_output_mode": {"value": 4, "tokens": ["--alignment-output-mode", "4"]},
            "alt_ali": {"value": 2, "tokens": ["--alt-ali", "2"]},
            "max_rejected": {"value": 2, "tokens": ["--max-rejected", "2"]},
            "max_accept": {"value": 2, "tokens": ["--max-accept", "2"]},
            "score_bias": {"value": 0.25, "tokens": ["--score-bias", "0.25"]},
            "realign": {"value": True, "tokens": ["--realign", "1"]},
            "zdrop": {"value": 2, "tokens": ["--zdrop", "2"]},
            "wrapped_scoring": {"value": True, "tokens": ["--wrapped-scoring", "1"]},
            "comp_bias_corr": {"value": True, "tokens": ["--comp-bias-corr", "1"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_rescorediagonal",
        "command": "rescorediagonal",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "prefilter_db"], ["O", "result_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "rescore_mode": {"value": 4, "tokens": ["--rescore-mode", "4"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "add_backtrace": {"value": True, "tokens": ["-a", "1"]},
            "filter_hits": {"value": True, "tokens": ["--filter-hits", "1"]},
            "sort_results": {"value": 1, "tokens": ["--sort-results", "1"]},
            "wrapped_scoring": {"value": True, "tokens": ["--wrapped-scoring", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_fwbw",
        "command": "fwbw",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "alignment_db"], ["O", "result_db"]],
        "flags": {
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "seq_id_mode": {"value": 2, "tokens": ["--seq-id-mode", "2"]},
            "mact": {"value": 0.25, "tokens": ["--mact", "0.25"]},
            "fwbw_gapopen": {"value": 0.25, "tokens": ["--fwbw-gapopen", "0.25"]},
            "fwbw_gapextend": {"value": 0.25, "tokens": ["--fwbw-gapextend", "0.25"]},
            "temperature": {"value": 0.25, "tokens": ["--temperature", "0.25"]},
            "blocklen": {"value": 2, "tokens": ["--blocklen", "2"]},
            "fwbw_backtrace_mode": {"value": 1, "tokens": ["--fwbw-backtrace-mode", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_alignbykmer",
        "command": "alignbykmer",
        "pos": [["D", "query_db"], ["D", "target_db"], ["D", "result_db"], ["O", "alignment_db"]],
        "flags": {
            "kmer_size": {"value": 2, "tokens": ["-k", "2"]},
            "spaced_kmer_mode": {"value": True, "tokens": ["--spaced-kmer-mode", "1"]},
            "e_value": {"value": 0.25, "tokens": ["-e", "0.25"]},
            "min_seq_id": {"value": 0.25, "tokens": ["--min-seq-id", "0.25"]},
            "coverage": {"value": 0.25, "tokens": ["-c", "0.25"]},
            "cov_mode": {"value": 5, "tokens": ["--cov-mode", "5"]},
            "min_aln_len": {"value": 2, "tokens": ["--min-aln-len", "2"]},
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
            "add_self_matches": {"value": True, "tokens": ["--add-self-matches", "1"]},
            "filter_hits": {"value": True, "tokens": ["--filter-hits", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_extractorfs",
        "command": "extractorfs",
        "pos": [["D", "sequence_db"], ["O", "orf_db"]],
        "flags": {
            "min_length": {"value": 2, "tokens": ["--min-length", "2"]},
            "max_length": {"value": 2, "tokens": ["--max-length", "2"]},
            "max_gaps": {"value": 2, "tokens": ["--max-gaps", "2"]},
            "contig_start_mode": {"value": 2, "tokens": ["--contig-start-mode", "2"]},
            "contig_end_mode": {"value": 2, "tokens": ["--contig-end-mode", "2"]},
            "orf_start_mode": {"value": 2, "tokens": ["--orf-start-mode", "2"]},
            "forward_frames": {"value": "1,3", "tokens": ["--forward-frames", "1,3"]},
            "reverse_frames": {"value": "1,3", "tokens": ["--reverse-frames", "1,3"]},
            "translation_table": {"value": 2, "tokens": ["--translation-table", "2"]},
            "translate": {"value": True, "tokens": ["--translate", "1"]},
            "use_all_table_starts": {"value": True, "tokens": ["--use-all-table-starts", "1"]},
            "id_offset": {"value": 2, "tokens": ["--id-offset", "2"]},
            "create_lookup": {"value": True, "tokens": ["--create-lookup", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_extractframes",
        "command": "extractframes",
        "pos": [["D", "sequence_db"], ["O", "frame_db"]],
        "flags": {
            "forward_frames": {"value": "1,3", "tokens": ["--forward-frames", "1,3"]},
            "reverse_frames": {"value": "1,3", "tokens": ["--reverse-frames", "1,3"]},
            "translation_table": {"value": 2, "tokens": ["--translation-table", "2"]},
            "translate": {"value": True, "tokens": ["--translate", "1"]},
            "create_lookup": {"value": True, "tokens": ["--create-lookup", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_orftocontig",
        "command": "orftocontig",
        "pos": [["D", "contig_db"], ["D", "orf_db"], ["O", "alignment_db"]],
        "flags": {},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_recoverlongestorf",
        "command": "recoverlongestorf",
        "pos": [["D", "orf_db"], ["D", "result_db"], ["T", "output"]],
        "flags": {},
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_splitsequence",
        "command": "splitsequence",
        "pos": [["D", "sequence_db"], ["O", "split_db"]],
        "flags": {
            "sequence_overlap": {"value": 2, "tokens": ["--sequence-overlap", "2"]},
            "max_seq_len": {"value": 2, "tokens": ["--max-seq-len", "2"]},
            "sequence_split_mode": {"value": 1, "tokens": ["--sequence-split-mode", "1"]},
            "headers_split_mode": {"value": 1, "tokens": ["--headers-split-mode", "1"]},
            "create_lookup": {"value": True, "tokens": ["--create-lookup", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_masksequence",
        "command": "masksequence",
        "pos": [["D", "sequence_db"], ["O", "masked_db"]],
        "flags": {
            "mask": {"value": True, "tokens": ["--mask", "1"]},
            "mask_prob": {"value": 0.25, "tokens": ["--mask-prob", "0.25"]},
            "mask_lower_case": {"value": True, "tokens": ["--mask-lower-case", "1"]},
            "mask_n_repeat": {"value": 2, "tokens": ["--mask-n-repeat", "2"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_offsetalignment",
        "command": "offsetalignment",
        "pos": [
            ["D", "query_db"],
            ["D", "query_orf_db"],
            ["D", "target_db"],
            ["D", "target_orf_db"],
            ["D", "alignment_db"],
            ["O", "offset_db"],
        ],
        "flags": {
            "search_type": {"value": 4, "tokens": ["--search-type", "4"]},
            "chain_alignments": {"value": True, "tokens": ["--chain-alignments", "1"]},
            "merge_query": {"value": True, "tokens": ["--merge-query", "1"]},
        },
        "threads": True,
        "tmp": False,
    },
    {
        "name": "mmseqs_proteinaln2nucl",
        "command": "proteinaln2nucl",
        "pos": [
            ["D", "nucleotide_query_db"],
            ["D", "nucleotide_target_db"],
            ["D", "protein_query_db"],
            ["D", "protein_target_db"],
            ["D", "alignment_db"],
            ["O", "nucleotide_alignment_db"],
        ],
        "flags": {
            "gap_open": {"value": 2, "tokens": ["--gap-open", "2"]},
            "gap_extend": {"value": 2, "tokens": ["--gap-extend", "2"]},
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


@requires_mmseqs
@pytest.mark.parametrize("function, command, extra", [
    ("mmseqs_prefilter", "prefilter", ["-s", "5.0", "--max-seqs", "50"]),
    ("mmseqs_ungappedprefilter", "ungappedprefilter", []),
    ("mmseqs_gappedprefilter", "gappedprefilter", []),
])  # fmt: skip
def test_prefilters_match_direct(world, tmp_path, function, command, extra):
    options = {"mmseqs_prefilter": {"sensitivity": 5.0, "max_seqs": 50}}.get(function, {})
    check_database(tmp_path, lambda out: getattr(module, function)(str(world["protDB"]), str(world["protDB"]), out, verbose=False, **options),
                   command, world["protDB"], world["protDB"], OUT, *extra)  # fmt: skip


@requires_mmseqs
def test_prefilter_edge_case_single_sequence(world, tmp_path):
    single = tmp_path / "single.fa"
    single.write_text(">only\n" + "M" + "ACDEFGHIK" * 10 + "\n", encoding="utf-8")
    toy.run_direct("createdb", single, tmp_path / "singleDB")
    result = check_database(sub(tmp_path, "edge"),
                            lambda out: module.mmseqs_prefilter(str(tmp_path / "singleDB"), str(tmp_path / "singleDB"), out, verbose=False),
                            "prefilter", tmp_path / "singleDB", tmp_path / "singleDB", OUT)  # fmt: skip
    assert result["entry_counts"]["prefilter_db"] == 1


@requires_mmseqs
def test_kmermatcher_and_kmersearch_match_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_kmermatcher(str(world["protDB"]), out, verbose=False),
                   "kmermatcher", world["protDB"], OUT)  # fmt: skip
    if not Path(f"{world['ntDB']}.linidx.dbtype").exists():
        toy.run_direct("createlinindex", world["ntDB"], tmp_path / "t", "--search-type", "3")
    check_database(sub(tmp_path, "ks"), lambda out: module.mmseqs_kmersearch(str(world["ntDB"]), f"{world['ntDB']}.linidx", out, verbose=False),
                   "kmersearch", world["ntDB"], f"{world['ntDB']}.linidx", OUT)  # fmt: skip


@requires_mmseqs
def test_alignall_transitivealign_fwbw_match_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_alignall(str(world["protDB"]), str(world["clu"]), out, verbose=False),
                   "alignall", world["protDB"], world["clu"], OUT)  # fmt: skip
    a = tmp_path / "a"
    a.mkdir()
    check_database(a, lambda out: module.mmseqs_transitivealign(str(world["protDB"]), str(world["aln"]), out, verbose=False),
                   "transitivealign", world["protDB"], world["aln"], OUT)  # fmt: skip
    b = tmp_path / "b"
    b.mkdir()
    check_database(b, lambda out: module.mmseqs_fwbw(str(world["protDB"]), str(world["protDB"]), str(world["aln"]), out, mact=0.05, verbose=False),
                   "fwbw", world["protDB"], world["protDB"], world["aln"], OUT, "--mact", "0.05")  # fmt: skip


@requires_mmseqs
@pytest.mark.parametrize("rescore_mode", [0, 2, 3])
def test_rescorediagonal_matches_direct(world, tmp_path, rescore_mode):
    check_database(tmp_path, lambda out: module.mmseqs_rescorediagonal(str(world["protDB"]), str(world["protDB"]), str(world["pf"]), out, rescore_mode=rescore_mode, verbose=False),
                   "rescorediagonal", world["protDB"], world["protDB"], world["pf"], OUT, "--rescore-mode", str(rescore_mode))  # fmt: skip


@requires_mmseqs
def test_alignbykmer_matches_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_alignbykmer(str(world["protDB"]), str(world["protDB"]), str(world["pf"]), out, verbose=False),
                   "alignbykmer", world["protDB"], world["protDB"], world["pf"], OUT)  # fmt: skip


@requires_mmseqs
def test_orf_extraction_matches_direct(world, tmp_path):
    result = check_database(tmp_path, lambda out: module.mmseqs_extractorfs(str(world["ctgDB"]), out, translate=True, min_length=30, verbose=False),
                            "extractorfs", world["ctgDB"], OUT, "--translate", "1", "--min-length", "30")  # fmt: skip
    assert result["entry_counts"]["orf_db"] == 56
    a = tmp_path / "frames"
    a.mkdir()
    check_database(a, lambda out: module.mmseqs_extractframes(str(world["ctgDB"]), out, forward_frames="1,2", reverse_frames="1", translate=True, verbose=False),
                   "extractframes", world["ctgDB"], OUT, "--forward-frames", "1,2", "--reverse-frames", "1", "--translate", "1")  # fmt: skip
    b = tmp_path / "contig"
    b.mkdir()
    check_database(b, lambda out: module.mmseqs_orftocontig(str(world["ctgDB"]), str(world["orfaa"]), out, verbose=False),
                   "orftocontig", world["ctgDB"], world["orfaa"], OUT)  # fmt: skip


@requires_mmseqs
def test_recoverlongestorf_matches_direct(world, tmp_path):
    out, ref = tmp_path / "facade.tsv", tmp_path / "direct.tsv"
    module.mmseqs_recoverlongestorf(
        str(world["orfaa"]), str(world["orfaln"]), str(out), verbose=False
    )
    toy.run_direct("recoverlongestorf", world["orfaa"], world["orfaln"], ref)
    assert out.read_bytes() == ref.read_bytes()


@requires_mmseqs
def test_splitsequence_and_masksequence_match_direct(world, tmp_path):
    result = check_database(tmp_path, lambda out: module.mmseqs_splitsequence(str(world["ctgDB"]), out, max_seq_len=300, sequence_overlap=30, verbose=False),
                            "splitsequence", world["ctgDB"], OUT, "--max-seq-len", "300", "--sequence-overlap", "30")  # fmt: skip
    assert result["entry_counts"]["split_db"] == 16
    a = tmp_path / "mask"
    a.mkdir()
    check_database(a, lambda out: module.mmseqs_masksequence(str(world["protDB"]), out, mask_prob=0.5, verbose=False),
                   "masksequence", world["protDB"], OUT, "--mask-prob", "0.5")  # fmt: skip


@requires_mmseqs
def test_offsetalignment_and_proteinaln2nucl_match_direct(world, tmp_path):
    check_database(tmp_path, lambda out: module.mmseqs_offsetalignment(str(world["ctgDB"]), str(world["orfaa"]), str(world["refDB"]), str(world["refDB"]), str(world["orfaln"]), out, verbose=False),
                   "offsetalignment", world["ctgDB"], world["orfaa"], world["refDB"], world["refDB"], world["orfaln"], OUT)  # fmt: skip
    a = tmp_path / "p2n"
    a.mkdir()
    result = check_database(a, lambda out: module.mmseqs_proteinaln2nucl(str(world["cdsDB"]), str(world["cdsDB"]), str(world["prot0DB"]), str(world["prot0DB"]), str(world["aln0"]), out, verbose=False),
                            "proteinaln2nucl", world["cdsDB"], world["cdsDB"], world["prot0DB"], world["prot0DB"], world["aln0"], OUT)  # fmt: skip
    assert result["entry_counts"]["nucleotide_alignment_db"] == 24
