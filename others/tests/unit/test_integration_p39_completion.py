"""P39 wrappers compared to their upstream APIs on small local inputs."""

import json
import shutil
import subprocess
import sys

import pytest

from dnakit.fasta_fastq import cutadapt_p39_functions as cut
from dnakit.phylogeny import toytree_p39_functions as toy


@pytest.mark.parametrize("sequence", ["ACGT", "NNNN", "A"])
def test_cutadapt_filter_matches_cli(tmp_path, sequence):
    pytest.importorskip("cutadapt")
    source, actual, expected = (tmp_path / x for x in ("r.fastq", "a.fastq", "b.fastq"))
    source.write_text(f"@r\n{sequence}\n+\n{'I' * len(sequence)}\n")
    result = cut.cutadapt_filter_reads(str(source), str(actual), max_n=0, verbose=False)
    assert result["status"] == "completed", result
    subprocess.run(
        [sys.executable, "-m", "cutadapt", "--max-n", "0", "-o", str(expected), str(source)],
        capture_output=True,
        check=True,
    )
    assert actual.read_bytes() == expected.read_bytes()


@pytest.mark.parametrize("template", ["{id}", "{id}_changed", "prefix_{id}"])
def test_cutadapt_rename(tmp_path, template):
    pytest.importorskip("cutadapt")
    source, output = tmp_path / "r.fastq", tmp_path / "o.fastq"
    source.write_text("@read\nACGT\n+\nIIII\n")
    result = cut.cutadapt_rewrite_read_names(
        str(source), str(output), rename=template, verbose=False
    )
    assert result["status"] == "completed", result
    assert output.read_text().splitlines()[0] == "@" + template.format(id="read")


def test_cutadapt_demultiplex(tmp_path):
    pytest.importorskip("cutadapt")
    source = tmp_path / "r.fastq"
    source.write_text("@read\nAAAACGT\n+\nIIIIIII\n")
    result = cut.cutadapt_demultiplex_reads(
        str(source), str(tmp_path / "demux"), {"indexA": "^AAA"}, verbose=False
    )
    assert result["status"] == "completed", result
    assert "ACGT" in (tmp_path / "demux/indexA.fastq").read_text()


@pytest.mark.parametrize(
    "first,second",
    [("((a,b),c);", "((a,b),c);"), ("((a,b),(c,d));", "((a,c),(b,d));"), ("(a,b);", "(a,b);")],
)
def test_toytree_rf_matches_upstream(first, second):
    module = pytest.importorskip("toytree")
    got = toy.toytree_rf_distance(first, second, verbose=False)
    assert got["distance"] == module.distance.get_treedist_rf(
        module.tree(first), module.tree(second)
    )


def test_toytree_tree_operations():
    pytest.importorskip("toytree")
    tree = "((a,b),c);"
    assert toy.toytree_mrca(tree, ["a", "b"], verbose=False)["tips"] == ["a", "b"]
    assert toy.toytree_monophyly(tree, ["a", "b"], verbose=False)["monophyletic"]
    assert len(toy.toytree_unique_topologies([tree, tree], verbose=False)["topologies"]) == 1
    assert toy.toytree_consensus([tree, tree], verbose=False)["status"] == "completed"


@pytest.mark.parametrize("distance", [0.0, 0.1, 1.0])
def test_toytree_k80_matches_upstream(distance):
    np = pytest.importorskip("numpy")
    module = pytest.importorskip("toytree.infer.src.maximum_likelihood")
    got = toy.toytree_k80_matrix(distance, verbose=False)["matrix"]
    np.testing.assert_allclose(got, module.K80().get_p_matrix_function()(distance))


@pytest.mark.parametrize("distance", [0.0, 0.1, 1.0])
def test_toytree_tn93_matches_upstream(distance):
    np = pytest.importorskip("numpy")
    module = pytest.importorskip("toytree.infer.src.maximum_likelihood")
    got = toy.toytree_tn93_matrix(distance, verbose=False)["matrix"]
    np.testing.assert_allclose(
        got,
        module.TN93().get_p_matrix_function(kappa_ls=(2.0, 2.0), pi_ls=(0.25, 0.25, 0.25, 0.25))(
            distance
        ),
    )


@pytest.mark.parametrize("size", [3, 4, 5])
@pytest.mark.parametrize("algorithm", ["neighbor_joining", "upgma"])
def test_toytree_distance_tree_matches_upstream(size, algorithm):
    module = pytest.importorskip("toytree")
    pd = pytest.importorskip("pandas")
    names = list("abcde"[:size])
    matrix = [[float(abs(i - j)) for j in range(size)] for i in range(size)]
    method = "neighbor_joining_tree" if algorithm == "neighbor_joining" else "upgma_tree"
    expected = getattr(module.infer, method)(
        pd.DataFrame(matrix, index=names, columns=names)
    ).write()
    result = getattr(toy, "toytree_" + algorithm)(matrix, names, verbose=False)
    assert result["newick"] == expected


def test_toytree_neighbor_joining_rejects_two_tips():
    pytest.importorskip("toytree")
    from dnakit.exceptions import ConfigurationError

    with pytest.raises(ConfigurationError, match="at least three"):
        toy.toytree_neighbor_joining([[0.0, 1.0], [1.0, 0.0]], ["a", "b"], verbose=False)


@pytest.mark.parametrize("selection", [{"index": 0}, {"position": 0}, {"position": 9}])
def test_toytree_local_genealogy_matches_upstream(tmp_path, selection):
    tskit = pytest.importorskip("tskit")
    module = pytest.importorskip("toytree.utils.src.toytree_sequence")
    tables = tskit.TableCollection(sequence_length=10)
    for _ in range(2):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(time=1)
    tables.edges.add_row(left=0, right=10, parent=2, child=0)
    tables.edges.add_row(left=0, right=10, parent=2, child=1)
    tables.sort()
    sequence = tables.tree_sequence()
    source = tmp_path / "trees.trees"
    sequence.dump(source)
    native = module.ToyTreeSequence(sequence)
    expected = (
        native.at_index(selection["index"])
        if "index" in selection
        else native.at_site(selection["position"])
    )
    assert (
        toy.toytree_local_genealogy(str(source), verbose=False, **selection)["newick"]
        == expected.write()
    )


@pytest.mark.parametrize("first,second", [("ACGT", "AGGT"), ("A", "T"), ("ACGT", "ACGT")])
def test_geneblocks_diff_upstream(first, second):
    upstream = pytest.importorskip("geneblocks")
    from dnakit.molecular_lab import geneblocks_p39_functions as module

    got = module.geneblocks_merge_differences(first, second, verbose=False)
    expected = upstream.DiffBlocks.from_sequences(first, second, blast_over=None).merged()
    assert [b["operation"] for b in got["blocks"]] == [b.operation for b in expected.blocks]
    assert (
        module.geneblocks_difference_features(first, second, verbose=False)["status"] == "completed"
    )


@pytest.mark.parametrize("case", ["near_match", "unrelated", "single_base"])
def test_geneblocks_common_pipeline_matches_upstream(tmp_path, case):
    upstream = pytest.importorskip("geneblocks")
    if shutil.which("blastn") is None:
        pytest.skip("BLAST+ blastn is not available")
    import random

    from geneblocks.biotools import annotate_record, sequence_to_record

    from dnakit.molecular_lab import geneblocks_p39_functions as module

    rng = random.Random(108)
    first = "".join(rng.choices("ACGT", k=400))
    second = first[:200] + ("A" if first[200] != "A" else "C") + first[201:]
    if case == "unrelated":
        second = "".join(rng.choices("ACGT", k=400))
    elif case == "single_base":
        first, second = "A", "T"
    sequences = {"first": first, "second": second}
    records = {
        name: sequence_to_record(seq, record_id=name, name=name) for name, seq in sequences.items()
    }
    native = upstream.CommonBlocks.from_sequences(records, min_block_size=80)

    def record_dict(record):
        return {
            "id": record.id,
            "sequence": str(record.seq),
            "features": [
                {
                    "start": int(f.location.start),
                    "end": int(f.location.end),
                    "strand": f.location.strand,
                    "type": f.type,
                    "qualifiers": f.qualifiers,
                }
                for f in record.features
            ],
        }

    got = module.geneblocks_find_common(sequences, verbose=False)
    assert got["blocks"] == native.common_blocks
    json.dumps(got, allow_nan=False)
    assert module.geneblocks_unique_segments(sequences, verbose=False)["segments"] == {
        name: [list(interval) for interval in intervals]
        for name, intervals in native.compute_unique_blocks().items()
    }
    assert module.geneblocks_common_records(sequences, verbose=False)["records"] == [
        record_dict(r) for r in native.common_blocks_records()
    ]
    assert module.geneblocks_unique_records(sequences, verbose=False)["records"] == [
        record_dict(r) for r in native.unique_blocks_records()
    ]
    assert module.geneblocks_annotate_common(sequences, verbose=False)["records"] == {
        name: record_dict(r) for name, r in native.sequences_with_annotated_blocks().items()
    }
    expected, actual = tmp_path / "expected.csv", tmp_path / "actual.csv"
    native.common_blocks_to_csv(str(expected))
    module.geneblocks_export_csv(sequences, str(actual), verbose=False)
    assert actual.read_bytes() == expected.read_bytes()
    # The feature-transfer path additionally verifies strand and coordinates.
    annotate_record(
        records["first"], location=(0, min(20, len(first)), -1), feature_type="misc_feature"
    )
    annotated_native = upstream.CommonBlocks.from_sequences(records, min_block_size=80)
    expected_records = annotated_native.copy_features_between_common_blocks(inplace=False)
    transferred = module.geneblocks_transfer_annotations(
        sequences,
        {"first": [{"start": 0, "end": min(20, len(first)), "strand": -1}]},
        verbose=False,
    )
    assert transferred["records"] == {name: record_dict(r) for name, r in expected_records.items()}
    image = tmp_path / "common.png"
    plotted = module.geneblocks_plot_common(sequences, str(image), verbose=False)
    assert image.is_file() and image.stat().st_size > 0
    assert plotted["blocks"] == native.common_blocks


def test_nupack_optional_dependency():
    from dnakit.exceptions import BackendUnavailableError
    from dnakit.secondary_structure.nupack_p39_functions import nupack_partition_function

    try:
        import nupack
    except ImportError:
        with pytest.raises(BackendUnavailableError, match="NUPACK"):
            nupack_partition_function(["ACGT"], verbose=False)
    else:
        got = nupack_partition_function(["ACGT"], verbose=False)
        assert got["free_energy"] == pytest.approx(
            float(nupack.pfunc(strands=["ACGT"], model=nupack.Model(material="dna"))[1])
        )


@pytest.mark.parametrize(
    "loop,structure", [("AATT", "...."), ("AACCCTT", "(.....)"), ("AA+TT", "((+))")]
)
def test_nupack_public_loop_energies(loop, structure):
    from dnakit.exceptions import BackendUnavailableError
    from dnakit.secondary_structure import nupack_p39_functions as module

    try:
        import nupack
    except ImportError:
        for function in (module.nupack_loop_energy, module.nupack_loop_stacking_energies):
            with pytest.raises(BackendUnavailableError, match="NUPACK"):
                function(loop, structure, verbose=False)
    else:
        model = nupack.Model(material="dna")
        assert module.nupack_loop_energy(loop, structure, verbose=False)[
            "free_energy"
        ] == pytest.approx(float(model.loop_energy(loop=loop, structure=structure)))
        assert module.nupack_loop_stacking_energies(loop, structure, verbose=False)[
            "stacking_energies"
        ] == pytest.approx(
            {
                str(key): float(value)
                for key, value in model.stack_energies(loop=loop, structure=structure).items()
            }
        )


@pytest.mark.parametrize(
    "records", [{"one": "ACGT"}, {"one": "N"}, {"first": "AC--", "second": "ACGT"}]
)
def test_pyhmmer_sequence_alignment_io(tmp_path, records):
    pytest.importorskip("pyhmmer")
    from dnakit.alignment import pyhmmer_p39_functions as module

    fasta, msa = tmp_path / "sequences.fa", tmp_path / "msa.afa"
    rows = [{"name": n, "sequence": s.replace("-", "")} for n, s in records.items()]
    module.pyhmmer_write_sequences(rows, str(fasta), verbose=False)
    assert [
        (r["name"], r["sequence"])
        for r in module.pyhmmer_read_sequences(str(fasta), verbose=False)["records"]
    ] == [(r["name"], r["sequence"]) for r in rows]
    module.pyhmmer_write_alignment(records, str(msa), verbose=False)
    assert [
        (r["name"], r["sequence"])
        for r in module.pyhmmer_read_alignment(str(msa), verbose=False)["records"]
    ] == list(records.items())
