"""Upstream-checked wrappers for selene-sdk, genomic-benchmarks and kipoi (group G17b)."""

from __future__ import annotations

import inspect
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import genomic_benchmarks_functions as gb_api
from dnakit.predictions import kipoi_functions as kipoi_api
from dnakit.predictions import selene_sdk_functions as sel_api
from dnakit.predictions.dna_model_backends import describe_kipoi_model, list_kipoi_models

SELENE_WRAPPERS = [
    sel_api.apply_one_hot_mutations,
    sel_api.decode_one_hot_sequence,
    sel_api.interval_length_probabilities,
    sel_api.list_in_silico_mutants,
    sel_api.load_variant_diff_scores,
    sel_api.order_variant_labels,
    sel_api.read_selene_vcf,
    sel_api.rescale_attribution_scores,
    sel_api.reverse_complement_allele,
    sel_api.reverse_complement_one_hot,
    sel_api.sort_standard_chromosomes,
    sel_api.unpack_bit_packed_sequence,
]
ALL_WRAPPERS = [*SELENE_WRAPPERS, gb_api.locate_sequences_in_reference]
ONE_HOT_ACGTN = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0.25] * 4]


def _rows(sequence: str) -> list[list[float]]:
    """One-hot rows in A, C, G, T order; N is uniform 0.25."""

    table = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]}
    return [list(map(float, table.get(base, [0.25] * 4))) for base in sequence]


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _hide_package(monkeypatch: pytest.MonkeyPatch, package: str) -> None:
    """Make ``package`` and every loaded submodule unimportable for the duration of a test."""

    for name in [key for key in sys.modules if key == package or key.startswith(package + ".")]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, package, None)


def _genome(tmp_path: Path) -> Path:
    return _write(tmp_path / "ref.fa", ">chr1\nACGTACGTAC\n>chr2\nGGGGTTTTAA\n")


def test_every_wrapper_has_native_keyword_only_verbose() -> None:
    for function in ALL_WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_kipoi_module_has_no_offline_dna_wrappers() -> None:
    assert list(kipoi_api.__all__) == []


def test_decode_one_hot_sequence_matches_upstream_and_hand_values(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    import numpy as np
    from selene_sdk.sequences import sequence as upstream

    result = sel_api.decode_one_hot_sequence(ONE_HOT_ACGTN, verbose=False)
    assert result["sequence"] == "ACGTN"
    assert result["operation"] == "decode_one_hot_sequence"
    assert result["provenance"]["package"] == "selene-sdk"
    direct = upstream.encoding_to_sequence(
        np.asarray(ONE_HOT_ACGTN, float), ["A", "C", "G", "T"], "N"
    )
    assert result["sequence"] == direct
    assert (
        sel_api.decode_one_hot_sequence(_rows("TTGA"), unk_base="A", verbose=False)["sequence"]
        == "TTGA"
    )
    with pytest.raises(ConfigurationError, match="one-hot"):
        sel_api.decode_one_hot_sequence([[0.5, 0.5, 0.0, 0.0]], verbose=False)
    with pytest.raises(ConfigurationError, match="finite"):
        sel_api.decode_one_hot_sequence([[True, 0, 0, 0]], verbose=False)


def test_reverse_complement_one_hot_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    import numpy as np
    from selene_sdk.predict import model_predict
    from selene_sdk.sequences.genome import Genome

    rows = _rows("AACG")
    result = sel_api.reverse_complement_one_hot(rows, verbose=False)
    assert sel_api.decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"] == "CGTT"
    direct = model_predict.get_reverse_complement_encoding(
        np.asarray(rows, float), list(Genome.BASES_ARR), dict(Genome.COMPLEMENTARY_BASE_DICT)
    )
    assert result["encoding"] == direct.tolist()
    assert result["length"] == 4


def test_reverse_complement_allele_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk.predict import model_predict
    from selene_sdk.sequences.genome import Genome

    assert sel_api.reverse_complement_allele("AC", verbose=False)["reverse_complement"] == "GT"
    assert sel_api.reverse_complement_allele("ttg", verbose=False)["reverse_complement"] == "CAA"
    assert sel_api.reverse_complement_allele("-", verbose=False)["reverse_complement"] == "*"
    assert sel_api.reverse_complement_allele("*", verbose=False)["reverse_complement"] == "*"
    direct = model_predict.get_reverse_complement("GATTACA", dict(Genome.COMPLEMENTARY_BASE_DICT))
    assert (
        sel_api.reverse_complement_allele("GATTACA", verbose=False)["reverse_complement"] == direct
    )
    with pytest.raises(ConfigurationError):
        sel_api.reverse_complement_allele("A-C", verbose=False)


def test_apply_one_hot_mutations_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    import numpy as np
    from selene_sdk.predict import model_predict
    from selene_sdk.sequences.genome import Genome

    rows = _rows("ACGT")
    result = sel_api.apply_one_hot_mutations(rows, [[0, "T"], [2, "a"]], verbose=False)
    assert sel_api.decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"] == "TCAT"
    assert result["mutation_count"] == 2
    assert result["mutations"] == [{"position": 0, "base": "T"}, {"position": 2, "base": "A"}]
    direct = model_predict.mutate_sequence(np.asarray(rows, float), [(0, "T"), (2, "A")], Genome)
    assert result["encoding"] == direct.tolist()
    with pytest.raises(ConfigurationError, match="repeats"):
        sel_api.apply_one_hot_mutations(rows, [[1, "A"], [1, "G"]], verbose=False)
    with pytest.raises(ConfigurationError, match="0, 3"):
        sel_api.apply_one_hot_mutations(rows, [[4, "A"]], verbose=False)
    with pytest.raises(ConfigurationError, match="A, C, G, T"):
        sel_api.apply_one_hot_mutations(rows, [[1, "N"]], verbose=False)


def test_list_in_silico_mutants_matches_upstream_and_counts() -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk.predict import model_predict

    result = sel_api.list_in_silico_mutants("ACG", verbose=False)
    assert result["count"] == 9
    assert result["mutants"][0] == [[0, "C"]]
    assert result["mutants"][-1] == [[2, "T"]]
    upstream = model_predict.in_silico_mutagenesis_sequences("ACG", 1)
    assert result["mutants"] == [[[int(p), b] for p, b in mutant] for mutant in upstream]
    pairs = sel_api.list_in_silico_mutants("AC", mutate_n_bases=2, verbose=False)
    assert pairs["count"] == 9
    assert pairs["mutants"][0] == [[0, "C"], [1, "A"]]
    shifted = sel_api.list_in_silico_mutants("ACG", start_position=1, verbose=False)
    assert shifted["count"] == 6
    assert {position for mutant in shifted["mutants"] for position, _ in mutant} == {1, 2}
    with pytest.raises(ConfigurationError, match="max_mutants"):
        sel_api.list_in_silico_mutants(
            "ACGTACGTAC", mutate_n_bases=4, max_mutants=10, verbose=False
        )
    with pytest.raises(ConfigurationError, match="shorter"):
        sel_api.list_in_silico_mutants("ACG", mutate_n_bases=2, start_position=2, verbose=False)


def test_read_selene_vcf_matches_upstream_and_flags_reference_mismatch(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk.predict import model_predict
    from selene_sdk.sequences.genome import Genome

    fasta = _genome(tmp_path)
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    vcf = _write(
        tmp_path / "ok.vcf", header + "1\t3\trs1\tG\tT\t.\t.\t.\nchr2\t5\t.\tT\tC\t.\t.\t.\n"
    )
    result = sel_api.read_selene_vcf(str(vcf), str(fasta), verbose=False)
    assert result["variant_count"] == 2
    first = result["variants"][0]
    assert (first["chromosome"], first["position"], first["ref"], first["alt"]) == (
        "chr1",
        3,
        "G",
        "T",
    )
    assert first["ref_matches_fasta"] is True and first["in_bounds"] is True
    upstream = model_predict.read_vcf_file(str(vcf), reference_sequence=Genome(str(fasta)))
    assert [(v[0], v[1], v[2], v[3], v[4]) for v in upstream] == [
        (v["chromosome"], v["position"], v["id"], v["ref"], v["alt"]) for v in result["variants"]
    ]
    bad = _write(
        tmp_path / "bad.vcf", header + "chr1\t3\trs1\tT\tG\t.\t.\t.\nchr1\t99\t.\tA\tC\t.\t.\t.\n"
    )
    flagged = sel_api.read_selene_vcf(str(bad), str(fasta), verbose=False)["variants"]
    assert flagged[0]["ref_matches_fasta"] is False
    assert flagged[1]["in_bounds"] is False and flagged[1]["ref_matches_fasta"] is None
    unknown = _write(tmp_path / "unknown.vcf", header + "chr9\t1\t.\tA\tC\t.\t.\t.\n")
    with pytest.raises(ConfigurationError, match="chr9"):
        sel_api.read_selene_vcf(str(unknown), str(fasta), verbose=False)


def test_unpack_bit_packed_sequence_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    import numpy as np
    from selene_sdk.samplers import dataloader

    # Rows A, C, G, T, N with N as all four bits set: packbits gives one byte per channel.
    packed_rows = [[136, 72, 40, 24]]
    result = sel_api.unpack_bit_packed_sequence(packed_rows, 5, verbose=False)
    assert sel_api.decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"] == "ACGTN"
    direct = dataloader.unpackbits_sequence(np.asarray(packed_rows, dtype=np.uint8), 5)
    assert result["encoding"] == direct.tolist()
    bits = np.asarray([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 1, 1]])
    assert np.packbits(bits.astype(np.uint8), axis=-2).tolist() == packed_rows
    with pytest.raises(ConfigurationError, match="fewer"):
        sel_api.unpack_bit_packed_sequence([[0, 0, 0, 0]], 9, verbose=False)
    with pytest.raises(ConfigurationError, match="0 to 255"):
        sel_api.unpack_bit_packed_sequence([[256, 0, 0, 0]], 1, verbose=False)


def test_interval_length_probabilities_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk import utils

    result = sel_api.interval_length_probabilities([100, 300, 600], verbose=False)
    assert result["indices"] == [0, 1, 2]
    assert result["probabilities"] == pytest.approx([0.1, 0.3, 0.6])
    dropped = sel_api.interval_length_probabilities([0, 300, 600], verbose=False)
    assert dropped["dropped_indices"] == [0]
    assert dropped["probabilities"] == pytest.approx([1 / 3, 2 / 3])
    subset = sel_api.interval_length_probabilities([100, 300, 600], indices=[2, 0], verbose=False)
    assert subset["indices"] == [2, 0]
    assert subset["probabilities"] == pytest.approx([6 / 7, 1 / 7])
    kept, weights = utils.get_indices_and_probabilities([100, 300, 600], [2, 0])
    assert list(kept) == subset["indices"]
    assert list(weights) == pytest.approx(subset["probabilities"])
    with pytest.raises(ConfigurationError, match="positive length"):
        sel_api.interval_length_probabilities([0, 0], verbose=False)
    with pytest.raises(ConfigurationError, match="repeat"):
        sel_api.interval_length_probabilities([1, 2], indices=[0, 0], verbose=False)


def test_rescale_attribution_scores_matches_hand_calculation() -> None:
    pytest.importorskip("selene_sdk")
    scores = [[1, 2, 3, 4], [0, 1, 0, 3]]
    identity = sel_api.rescale_attribution_scores(scores, verbose=False)
    assert identity["scores"] == scores
    both = sel_api.rescale_attribution_scores(
        scores, base_scaling="max_effect", position_scaling="max_effect", verbose=False
    )
    # Position 0: minus min 1 gives 0,1,2,3; sum 6; times max-min 3. Position 1 is analogous.
    expected = [0.0, 0.5, 1.0, 1.5, 0.0, 0.75, 0.0, 2.25]
    assert [value for row in both["scores"] for value in row] == pytest.approx(expected)
    with pytest.raises(ConfigurationError, match="divides by zero"):
        sel_api.rescale_attribution_scores(
            [[0, 0, 0, 0], [1, 2, 3, 4]], position_scaling="max_effect", verbose=False
        )
    with pytest.raises(ConfigurationError, match="base_scaling"):
        sel_api.rescale_attribution_scores(scores, base_scaling="sum", verbose=False)


def test_load_variant_diff_scores_matches_upstream(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk import interpret

    table = _write(
        tmp_path / "scores.tsv",
        "chrom\tpos\tname\tref\talt\tf1\tf2\nchr1\t10\trs1\tA\tG\t0.5\t0.25\nchr2\t20\trs2\tC\tT\t0.1\t0.0\n",
    )
    result = sel_api.load_variant_diff_scores(str(table), verbose=False)
    assert result["features"] == ["f1", "f2"]
    assert result["variants"][1] == {
        "chromosome": "chr2",
        "position": 20,
        "name": "rs2",
        "ref": "C",
        "alt": "T",
    }
    assert result["scores"] == [[0.5, 0.25], [0.1, 0.0]]
    scores, labels, features = interpret.load_variant_abs_diff_scores(str(table))
    assert scores.tolist() == result["scores"] and features == result["features"]
    assert [label[1] for label in labels] == ["10", "20"]
    no_features = _write(tmp_path / "none.tsv", "a\tb\tc\td\te\n")
    with pytest.raises(ConfigurationError, match="no feature"):
        sel_api.load_variant_diff_scores(str(no_features), verbose=False)
    bad_pos = _write(tmp_path / "pos.tsv", "a\tb\tc\td\te\tf1\nchr1\tx\trs\tA\tG\t0.1\n")
    with pytest.raises(ConfigurationError, match="not an integer"):
        sel_api.load_variant_diff_scores(str(bad_pos), verbose=False)


def test_order_variant_labels_is_input_order_independent() -> None:
    pytest.importorskip("selene_sdk")
    labels = [
        ["chr2", 30, "r1", "A", "C"],
        ["chr1", 20, "r2", "G", "T"],
        ["chr1", 10, "r3", "C", "A"],
        ["chrUn_x", 5, "r4", "A", "G"],
        ["chrX", 1, "r5", "A", "G"],
    ]
    result = sel_api.order_variant_labels(labels, verbose=False)
    assert result["order"] == [2, 1, 0, 4, 3]
    assert [v["chromosome"] for v in result["variants"]] == [
        "chr1",
        "chr1",
        "chr2",
        "chrX",
        "chrUn_x",
    ]
    shuffled = [labels[4], labels[0], labels[3], labels[2], labels[1]]
    again = sel_api.order_variant_labels(shuffled, verbose=False)
    assert [shuffled[i][2] for i in again["order"]] == [labels[i][2] for i in result["order"]]
    with pytest.raises(ConfigurationError, match="start with 'chr'"):
        sel_api.order_variant_labels([["1", 5, "r", "A", "G"]], verbose=False)


def test_sort_standard_chromosomes_matches_upstream_key() -> None:
    pytest.importorskip("selene_sdk")
    from selene_sdk.interpret.vis import sort_standard_chrs

    names = ["chr10", "chr2", "chrX", "chr1", "chrM", "chrUn_a", "chrY"]
    result = sel_api.sort_standard_chromosomes(names, verbose=False)
    assert [item["chromosome"] for item in result["chromosomes"]] == [
        "chr1",
        "chr2",
        "chr10",
        "chrX",
        "chrY",
        "chrM",
        "chrUn_a",
    ]
    assert [item["rank"] for item in result["chromosomes"]] == [1, 2, 10, 23, 24, 25, 26]
    assert all(
        item["rank"] == sort_standard_chrs(item["chromosome"]) for item in result["chromosomes"]
    )
    with pytest.raises(ConfigurationError, match="repeat"):
        sel_api.sort_standard_chromosomes(["chr1", "chr1"], verbose=False)


def test_selene_wrappers_raise_backend_unavailable_without_selene(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _hide_package(monkeypatch, "selene_sdk")
    with pytest.raises(BackendUnavailableError, match="bio-selene"):
        sel_api.decode_one_hot_sequence(ONE_HOT_ACGTN, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-selene"):
        sel_api.reverse_complement_allele("AC", verbose=False)


def test_locate_sequences_matches_upstream_and_flags_ambiguity(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks")
    from genomic_benchmarks.seq2loc import fasta2loc

    reference = {"chr1": "ACGTACGTAC", "chr2": "GGGGTTTTAA"}
    query = _write(
        tmp_path / "queries.fa",
        ">unique\nGGGG\n>repeat\nGTAC\n>absent\nAAAAA\n",
    )
    result = gb_api.locate_sequences_in_reference(str(query), reference, verbose=False)
    by_name = {match["sequence"]: match for match in result["matches"]}
    assert by_name["unique"] | {} == {
        "sequence": "unique",
        "chromosome": "chr2",
        "start": 0,
        "end": 4,
        "strand": "+",
        "occurrences": 1,
        "ambiguous": False,
    }
    # GTAC is its own reverse complement and occurs at 2 and 6 on chr1; upstream keeps the
    # last match scanned (6..10) and reports strand '-' for palindromes.
    assert by_name["repeat"]["occurrences"] == 2 and by_name["repeat"]["ambiguous"] is True
    assert (
        by_name["repeat"]["chromosome"],
        by_name["repeat"]["start"],
        by_name["repeat"]["end"],
    ) == (
        "chr1",
        6,
        10,
    )
    assert result["unmatched"] == ["absent"]
    upstream = fasta2loc(str(query), dict(sorted(reference.items())), use_seq_ids=True)
    for name, (chrom, start, end, strand) in upstream.items():
        assert (by_name[name]["chromosome"], by_name[name]["start"]) == (chrom, start)
        assert (by_name[name]["end"], by_name[name]["strand"]) == (end, strand)


def test_locate_sequences_rejects_reverse_complement_collisions(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks")
    query = _write(tmp_path / "pair.fa", ">a\nTTTT\n>b\nAAAA\n")
    with pytest.raises(ConfigurationError, match="reverse complements"):
        gb_api.locate_sequences_in_reference(str(query), {"chr1": "ACGT"}, verbose=False)


def test_locate_sequences_input_validation(tmp_path: Path) -> None:
    pytest.importorskip("genomic_benchmarks")
    lower = _write(tmp_path / "lower.fa", ">a\nacgt\n")
    with pytest.raises(ConfigurationError, match="uppercase"):
        gb_api.locate_sequences_in_reference(str(lower), {"chr1": "ACGT"}, verbose=False)
    dup = _write(tmp_path / "dup.fa", ">a\nAC\n>a\nGT\n")
    with pytest.raises(ConfigurationError, match="more than once"):
        gb_api.locate_sequences_in_reference(str(dup), {"chr1": "ACGT"}, verbose=False)
    many = _write(tmp_path / "many.fa", ">a\nAC\n>b\nGG\n>c\nTT\n")
    with pytest.raises(ConfigurationError, match="max_sequences"):
        gb_api.locate_sequences_in_reference(
            str(many), {"chr1": "ACGT"}, max_sequences=2, verbose=False
        )
    with pytest.raises(ConfigurationError, match="reference"):
        gb_api.locate_sequences_in_reference(str(lower), {}, verbose=False)
    with pytest.raises(ConfigurationError, match="A, C, G, T and N"):
        gb_api.locate_sequences_in_reference(str(lower), {"chr1": "ACGX"}, verbose=False)
    with pytest.raises(ConfigurationError, match="max_reference_length"):
        gb_api.locate_sequences_in_reference(
            str(lower), {"chr1": "ACGT" * 3}, max_reference_length=5, verbose=False
        )


def test_locate_sequences_raises_backend_unavailable_without_genomic_benchmarks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = _write(tmp_path / "q.fa", ">a\nAC\n")
    _hide_package(monkeypatch, "genomic_benchmarks")
    with pytest.raises(BackendUnavailableError, match="bio-genomic-benchmarks"):
        gb_api.locate_sequences_in_reference(str(query), {"chr1": "ACGT"}, verbose=False)


class _FakeKipoi:
    """Stands in for the kipoi package so no model zoo is cloned or listed."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, object]] = []

    def list_models(self, sources=None):
        import pandas as pd

        self.calls.append(("list_models", sources))
        return pd.DataFrame(
            [
                {
                    "model": "fake_model",
                    "version": "1.0",
                    "type": "keras",
                    "license": "MIT",
                    "inputs": "seq",
                    "targets": "y",
                }
            ]
        )

    def get_model_descr(self, name, source="kipoi"):
        self.calls.append(("get_model_descr", (name, source)))
        return SimpleNamespace(
            type="keras",
            info=SimpleNamespace(version="1.0", doc="fake doc", license="MIT", cite_as=None),
            schema=SimpleNamespace(inputs={"seq": 1}, targets={"y": 1}),
        )


def test_existing_kipoi_wrappers_use_a_mocked_kipoi(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("pandas")
    fake = _FakeKipoi()
    monkeypatch.setitem(sys.modules, "kipoi", fake)
    listed = list_kipoi_models(verbose=False)
    assert listed["models"] == [
        {
            "model": "fake_model",
            "version": "1.0",
            "type": "keras",
            "license": "MIT",
            "inputs": "seq",
            "targets": "y",
        }
    ]
    described = describe_kipoi_model("fake_model", verbose=False)
    assert described["description"]["inputs"] == ["seq"]
    assert described["description"]["targets"] == ["y"]
    assert fake.calls == [("list_models", None), ("get_model_descr", ("fake_model", "kipoi"))]


def test_kipoi_wrappers_raise_backend_unavailable_without_kipoi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(sys.modules, "kipoi", None)
    with pytest.raises(BackendUnavailableError, match="bio-kipoi"):
        describe_kipoi_model("fake_model", verbose=False)
