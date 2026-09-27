"""Scientific invariants and hand-computable cases for native bioinformatics APIs."""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from dnakit.bioinformatics import encoding, features, metrics, motifs, optimization, phylogeny
from dnakit.bioinformatics import statistics as statistics_api
from dnakit.core import DNA
from dnakit.exceptions import ConfigurationError

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")


def test_composition_encodings_and_reverse_complement_invariance() -> None:
    sequences = {"first": "AAACGT", "reverse": "ACGTTT"}
    nac = encoding.encode_dna(sequences, "NAC")
    assert nac.ids == ("first", "reverse")
    assert nac.names == ("A", "C", "G", "T")
    np.testing.assert_allclose(nac.values[0], [3 / 6, 1 / 6, 1 / 6, 1 / 6])
    canonical = encoding.encode_dna(sequences, "RCKmer", parameters={"k": 2})
    np.testing.assert_allclose(canonical.values[0], canonical.values[1])
    np.testing.assert_allclose(canonical.values.sum(axis=1), 1)
    with pytest.raises(ValueError):
        nac.values[0, 0] = 0
    anf = encoding.encode_dna({"a": "AACA"}, "ANF")
    np.testing.assert_allclose(anf.values, [[1, 1, 1 / 3, 3 / 4]])
    pairs = encoding.encode_dna({"a": "ACGT"}, "CKSNAP", parameters={"gap": 0})
    assert pairs.values[0, pairs.names.index("gap_0:AC")] == pytest.approx(1 / 3)


def test_native_encoding_accepts_dna_batch_with_generated_ids() -> None:
    result = encoding.encode_dna(DNA(["AACGTG", "CCGTAT"]), "NAC")

    assert result.ids == ("sequence_1", "sequence_2")
    np.testing.assert_allclose(
        result.values,
        [[2 / 6, 1 / 6, 2 / 6, 1 / 6], [1 / 6, 2 / 6, 1 / 6, 2 / 6]],
    )


@pytest.mark.parametrize(
    "descriptor",
    [name for name in encoding.DNA_ENCODINGS if name not in {"PSTNPss", "PSTNPds", "KNN"}],
)
def test_all_native_encodings_have_a_finite_stable_schema(descriptor: str) -> None:
    k = 3 if descriptor in {"TPCP", "TAC", "TCC", "TACC", "PCPseTNC", "SCPseTNC"} else 2
    tables = {
        name: {
            "".join(word): float(index**exponent)
            for index, word in enumerate(itertools.product("ACGT", repeat=k))
        }
        for name, exponent in (("linear", 1), ("square", 2))
    }
    result = encoding.encode_dna(
        {"one": "AACGTGCTACGATGCA", "two": "CCGTACCTAGTACGTA"},
        descriptor,
        parameters={"k": 2}
        if descriptor in {"Kmer", "RCKmer", "Mismatch", "Subsequence", "PseKNC"}
        else {},
        property_tables=tables,
    )
    assert result.values.shape == (2, len(result.names))
    assert len(set(result.names)) == len(result.names)
    assert np.isfinite(result.values).all()


@pytest.mark.parametrize("descriptor", ["PSTNPss", "PSTNPds", "KNN"])
def test_supervised_encodings_exclude_self_and_ignore_test_labels(descriptor: str) -> None:
    sequences = {"p1": "AAAAAA", "p2": "AAAAAT", "n1": "CCCCCC", "n2": "CCCCCG", "test": "AAAAAC"}
    training = ("p1", "p2", "n1", "n2")
    labels = {"p1": 1, "p2": 1, "n1": 0, "n2": 0, "test": 0}
    first = encoding.encode_dna(
        sequences,
        descriptor,
        labels=labels,
        training_ids=training,
        parameters={"neighbors": 1} if descriptor == "KNN" else {},
    )
    labels["test"] = 1
    second = encoding.encode_dna(
        sequences,
        descriptor,
        labels=labels,
        training_ids=training,
        parameters={"neighbors": 1} if descriptor == "KNN" else {},
    )
    np.testing.assert_array_equal(first.values, second.values)
    if descriptor == "KNN":
        assert first.values[-1, 0] == 1
    with pytest.raises(ConfigurationError):
        encoding.encode_dna(sequences, descriptor)


def test_encoding_undefined_inputs_are_not_silent_zero_features() -> None:
    with pytest.raises(ConfigurationError):
        encoding.encode_dna({"a": "AN"}, "NAC")
    with pytest.raises(ConfigurationError):
        encoding.encode_dna({"a": "AC", "b": "ACG"}, "binary")
    with pytest.raises(ConfigurationError):
        encoding.encode_dna({"a": "A"}, "CKSNAP")
    with pytest.raises(ConfigurationError):
        encoding.encode_dna({"a": "ACGTACGT"}, "DAC")


def test_kernel_distances_against_explicit_small_matrices() -> None:
    x, y = np.array([[0.0], [1.0]]), np.array([[2.0], [3.0]])
    kxx = np.exp(-((x - x.T) ** 2) / 2)
    kyy = np.exp(-((y - y.T) ** 2) / 2)
    kxy = np.exp(-((x - y.T) ** 2) / 2)
    expected = kxx[0, 1] + kyy[0, 1] - 2 * kxy.mean()
    assert metrics.kernel_distance(x, y) == pytest.approx(expected)
    assert metrics.kernel_distance(x, x, unbiased=False) == pytest.approx(0)
    assert metrics.kernel_distance(x, x) < 0  # Unbiased estimates must not be clamped.
    expected_kid = (
        ((x @ x.T + 1) ** 3)[0, 1] + ((y @ y.T + 1) ** 3)[0, 1] - 2 * ((x @ y.T + 1) ** 3).mean()
    )
    assert metrics.kernel_distance(x, y, metric="kid") == pytest.approx(expected_kid)
    with pytest.raises(ConfigurationError):
        metrics.kernel_distance(x, y, max_pairs=1)


def test_manifolds_authenticity_and_calibration() -> None:
    reference = [[0.0], [1.0], [2.0], [3.0]]
    assert metrics.manifold_metrics(reference, reference, neighbors=1)["precision"] == 1
    assert metrics.manifold_metrics([[100.0], [101.0]], reference, neighbors=1)["precision"] == 0
    assert metrics.authenticity(reference, reference) == 0
    assert metrics.authenticity([[20.0]], reference) == 1
    clipped = metrics.clipped_manifold_metrics(reference, reference, neighbors=1)
    assert clipped == {"clipped_density": 1.0, "clipped_coverage": 1.0}
    assert metrics.clipped_manifold_metrics(
        [[100.0], [101.0], [102.0], [103.0]], reference, neighbors=1
    ) == {"clipped_density": 0.0, "clipped_coverage": 0.0}


def test_pareto_hypervolume_and_property_statistics() -> None:
    points = [[1, 1], [2, 2], [1, 2], [1, 1]]
    np.testing.assert_array_equal(
        metrics.pareto_ranks(points, maximize=[False, False]), [0, 2, 1, 0]
    )
    assert metrics.hypervolume([[1, 2], [2, 1]], [3, 3], maximize=[False, False]) == pytest.approx(
        3
    )
    assert metrics.hypervolume([[1, 1, 1]], [3, 3, 3], maximize=[False] * 3) == pytest.approx(8)
    assert metrics.property_hit_rate([[1, 2], [3, 4]], lower=[0, 0], upper=[2, 2])["joint"] == 0.5
    np.testing.assert_allclose(metrics.property_kl(points, points), 0)
    conformity = metrics.conformity_score([[1, 1]], points, folds=2)
    assert 0 <= conformity["value"] <= 1


def test_fitted_feature_pipeline_does_not_refit_held_out_data() -> None:
    pytest.importorskip("sklearn")
    fitted = features.analyze_features(
        [[0, 1], [2, 3], [4, 5]], method="none", feature_names=["a", "b"]
    )
    centers = fitted.estimator.named_steps["scale"].mean_.copy()
    transformed = fitted.transform([[100, 101]], feature_names=["a", "b"])
    assert transformed[0, 0] > 10
    np.testing.assert_array_equal(centers, fitted.estimator.named_steps["scale"].mean_)
    with pytest.raises(ConfigurationError):
        fitted.transform([[1, 2]], feature_names=["b", "a"])


def test_native_distance_statistics_match_reference_implementation() -> None:
    skbio = pytest.importorskip("skbio")
    spatial = pytest.importorskip("scipy.spatial.distance")
    pdist, squareform = spatial.pdist, spatial.squareform

    data = squareform(pdist([[0, 0], [1, 0.2], [5, 5], [6, 5.1], [5.5, 6]]))
    groups = ["a", "a", "b", "b", "b"]
    for method in ("permanova", "anosim"):
        own = statistics_api.distance_test(data, groups, method=method, permutations=0)
        other = getattr(skbio.stats.distance, method)(
            skbio.DistanceMatrix(data), groups, permutations=0
        )
        assert own.statistic == pytest.approx(other["test statistic"])
    first = statistics_api.distance_test(data, groups, permutations=19, seed=9)
    assert first == statistics_api.distance_test(data, groups, permutations=19, seed=9)
    assert first.pvalue is not None and 0 < first.pvalue <= 1


def test_native_ordination_and_composition_roundtrips() -> None:
    spatial = pytest.importorskip("scipy.spatial.distance")
    pdist, squareform = spatial.pdist, spatial.squareform

    points = np.array([[0, 0], [1, 1], [2, 0], [3, 2.0]])
    distances = squareform(pdist(points))
    result = statistics_api.pcoa(distances)
    np.testing.assert_allclose(squareform(pdist(result.samples)), distances, atol=1e-8)
    assert not result.negative_eigenvalues
    counts = np.array([[1, 2, 3], [3, 2, 1], [2, 4, 1.0]])
    for method in ("clr", "ilr", "alr"):
        transformed = statistics_api.composition_transform(counts, method=method)
        recovered = statistics_api.composition_transform(transformed, method=method, inverse=True)
        np.testing.assert_allclose(recovered, counts / counts.sum(axis=1, keepdims=True))
    correspondence = statistics_api.ordination(counts)
    assert correspondence.eigenvalues.sum() > 0
    assert statistics_api.alpha_diversity([[1, 1]], metric="shannon")[0] == pytest.approx(np.log(2))


def test_native_tree_and_unifrac(tmp_path: Path) -> None:
    distances = [[0, 2, 4, 4], [2, 0, 4, 4], [4, 4, 0, 2], [4, 4, 2, 0]]
    for method in ("nj", "upgma"):
        tree = phylogeny.build_tree(distances, ["a", "b", "c", "d"], method=method)
        assert set(tree.tips) == {"a", "b", "c", "d"}
        assert all(node.length >= 0 for node in tree.walk())
        pytest.importorskip("Bio")
        path = tmp_path / f"{method}.nwk"
        phylogeny.write_tree(tree, path)
        assert phylogeny.read_tree(path).tips == tree.tips
    simple = phylogeny.TreeNode(children=(phylogeny.TreeNode("a", 1), phylogeny.TreeNode("b", 1)))
    distances = phylogeny.unifrac([[1, 0], [0, 1], [1, 0]], ["a", "b"], simple)
    assert distances[0, 1] == 1
    assert distances[0, 2] == 0


def test_motif_background_threshold_and_circular_coordinates() -> None:
    model = motifs.build_motif(["ACG", "ACG", "ATG"], pseudocount=0.1)
    assert model.consensus == "ACG"
    hits = motifs.scan_motif_model("GTAC", model, threshold=0, circular=True, both_strands=False)
    assert any(hit["start"] == 2 and hit["end"] == 5 for hit in hits)
    threshold = motifs.fpr_threshold(model, fpr=0.1, precision=100)
    false_positive_probability = 0.0
    for word in itertools.product("ACGT", repeat=3):
        if motifs.scan_motif_model("".join(word), model, threshold=threshold, both_strands=False):
            false_positive_probability += 0.25**3
    assert false_positive_probability <= 0.1


def test_native_optimizer_locks_preserves_translation_and_reports_failure(tmp_path: Path) -> None:
    fixed = optimization.specification("avoid_changes", sequence="AAAAAA", positions=[0, 1])
    gc = optimization.specification("gc", minimum=0.5, maximum=0.5)
    result = optimization.optimize_sequence("AAAAAA", constraints=[fixed, gc])
    assert result.constraints_satisfied
    assert result.sequence[:2] == "AA"
    assert sum(base in "GC" for base in result.sequence) == 3
    coding = optimization.optimize_sequence(
        "GCTGCT",
        constraints=[optimization.specification("avoid_pattern", pattern="GCT")],
        preserve_translation=True,
    )
    assert coding.constraints_satisfied
    from dnakit.ops import translate

    assert translate(coding.original) == translate(coding.sequence)
    impossible = optimization.optimize_sequence(
        "AAAA", constraints=[optimization.specification("length", minimum=5)], max_evaluations=2
    )
    assert impossible.status == "search_exhausted"
    assert not impossible.constraints_satisfied
    report = impossible.write_report(tmp_path / "report.json")
    assert '"constraints_satisfied": false' in report.read_text()


def test_mmi_sums_to_positional_information() -> None:
    text = "AAACGCGTTAGCA"
    output = encoding.encode_dna({"sample": text}, "MMI")
    joint = np.zeros((4, 4))
    for a, b in itertools.pairwise(text):
        joint["ACGT".index(a), "ACGT".index(b)] += 1
    joint /= joint.sum()
    expected = 0.0
    for i, j in itertools.product(range(4), repeat=2):
        if joint[i, j]:
            expected += joint[i, j] * np.log2(joint[i, j] / (joint.sum(1)[i] * joint.sum(0)[j]))
    assert output.values[0, :10].sum() == pytest.approx(expected)
    with pytest.raises(ConfigurationError, match="Unknown"):
        encoding.encode_dna({"x": text}, "NAC", parameters={"typo": 1})
    with pytest.raises(ConfigurationError, match="max_output_cells"):
        encoding.encode_dna({"x": text * 1000}, "PS4", max_output_cells=10)


def test_native_mcl_and_tree_topology_operations() -> None:
    groups = features.cluster_features(
        [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]], method="mcl"
    )
    assert groups == ((0, 1), (2, 3))
    node = phylogeny.TreeNode
    abcd = node(
        children=(
            node(children=(node("a", 1), node("b", 1)), length=1),
            node(children=(node("c", 1), node("d", 1)), length=1),
        )
    )
    acbd = node(
        children=(node(children=(node("a"), node("c"))), node(children=(node("b"), node("d"))))
    )
    assert phylogeny.robinson_foulds(abcd, acbd) == 2
    assert phylogeny.robinson_foulds(abcd, acbd, normalized=True) == 1
    majority = phylogeny.consensus_tree([abcd, abcd, acbd])
    assert phylogeny.robinson_foulds(majority, abcd, rooted=True) == 0
    strict = phylogeny.consensus_tree([abcd, acbd], threshold=1)
    assert len(strict.children) == 4
    pruned = phylogeny.prune_tree(abcd, ["a", "c"])
    assert pruned.tips == ("a", "c")
    assert [child.length for child in pruned.children] == [2, 2]


def test_extended_specifications_and_genbank_annotations() -> None:
    from dnakit.ops import translate

    usage = {"".join(word): 1.0 for word in itertools.product("ACGT", repeat=3)}
    usage["GCT"] = 0.01
    rare = optimization.specification("avoid_rare_codons", usage=usage, threshold=0.1)
    assert rare.evaluate("GCTGCC").locations == ((0, 3),)
    optimized = optimization.optimize_sequence(
        "GCTGCC", constraints=[rare], preserve_translation=True
    )
    assert optimized.constraints_satisfied
    assert translate(optimized.sequence) == "AA"
    compatibility = optimization.specification(
        "regions_compatibility", regions=[(0, 4), (4, 8)], match_length=4
    )
    assert compatibility.evaluate("AAAATTTT").violation == 1
    both = optimization.specification("avoid_pattern", pattern="AGG", both_strands=True)
    assert both.evaluate("CCT").violation == 1
    for parameters in ({"kind": "gc", "minimun": 0.4}, {"kind": "sequence"}):
        with pytest.raises(ConfigurationError):
            optimization.specification(**parameters)
    pytest.importorskip("Bio")
    Seq = pytest.importorskip("Bio.Seq").Seq
    bio_feature = pytest.importorskip("Bio.SeqFeature")
    BeforePosition, SeqFeature, SimpleLocation = (
        bio_feature.BeforePosition,
        bio_feature.SeqFeature,
        bio_feature.SimpleLocation,
    )
    SeqRecord = pytest.importorskip("Bio.SeqRecord").SeqRecord

    record = SeqRecord(Seq("ACGTACGT"))
    record.features = [
        SeqFeature(
            SimpleLocation(0, 4),
            qualifiers={"dnakit_constraint": ['{"kind":"gc", "minimum":0.4, "maximum":0.6}']},
        )
    ]
    constraints, objectives = optimization.specifications_from_record(record)
    assert len(constraints) == 1 and not objectives
    assert constraints[0].evaluate(str(record.seq)).violation == 0
    record.features[0].location = SimpleLocation(BeforePosition(0), 4)
    with pytest.raises(ConfigurationError):
        optimization.specifications_from_record(record)


def test_ancom_bc_matches_independent_reference() -> None:
    pytest.importorskip("skbio")
    reference = pytest.importorskip("skbio.stats.composition._ancombc")
    _estimate_bias_em, _estimate_params = reference._estimate_bias_em, reference._estimate_params

    from dnakit.bioinformatics.counts import ancom_bc

    rng = np.random.default_rng(13)
    design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
    data = np.exp(
        rng.normal(4, 0.7, (30, 60))
        + design[:, 1, None] * np.r_[np.full(6, 1.2), np.full(6, -1.0), np.zeros(48)]
    )
    actual = ancom_bc(data, design)
    variance, beta, _, _ = _estimate_params(np.log(data), design)
    expected = np.array(
        [_estimate_bias_em(beta[j], variance[:, j], tol=1e-5, max_iter=1000)[0] for j in range(2)]
    )
    np.testing.assert_allclose(actual.sampling_bias, expected, atol=2e-5)
    np.testing.assert_allclose(actual.standard_errors, np.sqrt(variance), atol=1e-12)
    np.testing.assert_allclose(actual.coefficients, beta.T - expected, atol=2e-5)
    with pytest.raises(ConfigurationError, match="positive"):
        ancom_bc(np.zeros_like(data), design)
    with pytest.raises(ConfigurationError, match="converge"):
        ancom_bc(data, design, max_iterations=1)


def test_nb2_matches_glm_and_reports_untestable_genes() -> None:
    sm = pytest.importorskip("statsmodels.api")
    from dnakit.bioinformatics.counts import adjust_pvalues, differential_expression

    rng = np.random.default_rng(18)
    design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
    mean = np.exp(design @ np.array([[3, 4, 5], [0.8, -1.0, 0.0]]))
    data = rng.negative_binomial(10, 10 / (mean + 10))
    data = np.column_stack([data, np.zeros(30)])
    actual = differential_expression(
        data, design, [0, 1], size_factors=np.ones(30), dispersions=0.1
    )
    assert actual.converged.tolist() == [True, True, True, False]
    assert np.isnan(actual.pvalues[-1]) and np.isnan(actual.qvalues[-1])
    for j in range(3):
        expected = sm.GLM(data[:, j], design, family=sm.families.NegativeBinomial(alpha=0.1)).fit()
        assert actual.log2_fold_change[j] == pytest.approx(expected.params[1] / np.log(2), abs=1e-5)
        assert actual.standard_errors[j] == pytest.approx(expected.bse[1] / np.log(2), abs=1e-5)
    np.testing.assert_allclose(adjust_pvalues([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
    with pytest.raises(ConfigurationError, match="column rank"):
        differential_expression(data, np.ones((30, 2)), [0, 1])


def test_information_gain_and_predictor_bridge() -> None:
    pytest.importorskip("sklearn")
    fitted = features.analyze_features(
        [[0, 1], [0, 1], [1, 1], [1, 1]],
        method="information_gain",
        labels=[0, 0, 1, 1],
        scale="none",
        selected_features=1,
        parameters={"bins": 2},
    )
    assert fitted.output_names == ("feature_0",)
    np.testing.assert_allclose(fitted.estimator.named_steps["analysis"].scores_, [1, 0])
    result = metrics.sequence_properties({"a": "ACGT", "b": "CC"}, {"length": len})
    assert result["ids"] == ("a", "b") and result["names"] == ("length",)
    np.testing.assert_array_equal(result["values"], [[4], [2]])


def test_sequence_properties_accepts_dna_batch_and_preserves_record_ids() -> None:
    sequences = DNA(["ACGT", "GGGG"])
    result = metrics.sequence_properties(
        sequences, {"gc": lambda text: sum(symbol in "GC" for symbol in text) / len(text)}
    )

    assert result["ids"] == ("sequence_1", "sequence_2")
    assert result["names"] == ("gc",)
    np.testing.assert_allclose(result["values"], [[0.5], [1.0]])

    duplicated = DNA([{"id": "same", "sequence": "ACGT"}, {"id": "same", "sequence": "GGGG"}])
    with pytest.raises(ConfigurationError, match="IDs must be unique"):
        metrics.sequence_properties(duplicated, {"length": len})
    raw_sequences: object = ["ACGT", "GGGG"]
    with pytest.raises(ConfigurationError, match="DNA"):
        metrics.sequence_properties(raw_sequences, {"length": len})  # type: ignore[arg-type]


def test_constrained_ordination_explains_fraction_of_total_inertia() -> None:
    x = np.arange(5.0)[:, None]
    y = np.column_stack([x[:, 0], [1.0, 0.0, 3.0, -1.0, 2.0]])
    result = statistics_api.ordination(y, method="rda", constraints=x)
    centered = y - y.mean(axis=0)
    projected = (x - x.mean()) @ np.linalg.lstsq(x - x.mean(), centered, rcond=None)[0]
    expected = (projected**2).sum() / (centered**2).sum()
    assert 0 < result.proportion_explained.sum() < 1
    assert result.proportion_explained.sum() == pytest.approx(expected)


def test_gaussian_mmd_is_translation_invariant_at_large_offsets() -> None:
    left, right = np.array([[0.0], [1.0], [2.0]]), np.array([[1.0], [3.0], [4.0]])
    assert metrics.kernel_distance(left, right) == pytest.approx(
        metrics.kernel_distance(left + 1e12, right + 1e12)
    )
