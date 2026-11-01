"""P40 completion checks against installed upstream implementations, without network."""
import json
from pathlib import Path

import pytest

from dnakit.databases import genomepy_p40_functions as genome
from dnakit.sequence import dnaio_p40_functions as dnaio
from dnakit.sequence import pyfaidx_p40_functions as faidx
from dnakit.sequence import pyfastx_p40_functions as fastx
from dnakit.annotation import genomicranges_p40_functions as gr


@pytest.fixture
def fasta(tmp_path):
    source = tmp_path / "tiny.fa"
    source.write_text(">chr1\nACGTACGT\n>chr2\nNNNN\n")
    return source


@pytest.mark.parametrize("sequence", ["ACGT", "N", "acgtN"])
def test_dnaio_reverse_complement(tmp_path, sequence):
    upstream = pytest.importorskip("dnaio")
    source, output = tmp_path / "r.fa", tmp_path / "rc.fa"
    source.write_text(f">read\n{sequence}\n")
    dnaio.dnaio_reverse_complement_reads(source, output, verbose=False)
    with upstream.open(output) as handle:
        got = next(iter(handle)).sequence
    expected = upstream.SequenceRecord("read", sequence).reverse_complement().sequence
    assert got == expected


@pytest.mark.parametrize("extension", [".gz", ".bz2", ".xz"])
def test_dnaio_compressed_chunks_are_reopenable(tmp_path, extension):
    upstream = pytest.importorskip("dnaio")
    from xopen import xopen

    source = tmp_path / ("reads.fastq" + extension)
    content = b"@read\nACGT\n+\nIIII\n"
    with xopen(source, "wb") as handle:
        handle.write(content)
    result = dnaio.dnaio_split_into_chunks(source, output_dir=tmp_path / "chunks", verbose=False)
    assert result["chunk_count"] == 1
    output = Path(result["outputs"][0])
    with xopen(output, "rb") as handle:
        assert handle.read() == content
    with upstream.open(output) as reader:
        assert next(iter(reader)).sequence == "ACGT"


def test_dnaio_chunk_limit_rejected_before_writing(tmp_path):
    pytest.importorskip("dnaio")
    from dnakit.exceptions import ConfigurationError
    source = tmp_path / "reads.fa"
    source.write_text(">read\nA\n")
    output = tmp_path / "chunks"
    with pytest.raises(ConfigurationError, match="max_chunks"):
        dnaio.dnaio_split_into_chunks(source, output_dir=output, max_chunks=0, verbose=False)
    assert not output.exists()


@pytest.mark.parametrize("sequence", ["A", "N" * 300, "ACGCGTAGTAGCAT" * 50])
def test_pyrodigal_gc_frames_preserve_unknown_and_truncation(tmp_path, sequence):
    upstream = pytest.importorskip("pyrodigal")
    from dnakit.annotation.pyrodigal_p40_functions import pyrodigal_gc_frame_plot
    source = tmp_path / "sequence.fa"
    source.write_text(f">sequence\n{sequence}\n")
    expected = list(upstream.Sequence(sequence).max_gc_frame_plot(10))
    result = pyrodigal_gc_frame_plot(source, window_size=10, max_segments=1, verbose=False)
    assert result["frame_counts"] == [expected.count(frame) for frame in (0, 1, 2)]
    assert result["unknown_frame_count"] == expected.count(-1)
    segment_count = sum(i == 0 or expected[i] != expected[i-1] for i in range(len(expected)))
    assert result["segments_truncated"] == (segment_count > 1)
    assert result["segments"][0]["frame"] == expected[0]


def test_faidx_copy_and_indexes(fasta, tmp_path):
    pytest.importorskip("pyfaidx")
    index = tmp_path / "custom.fai"
    assert faidx.faidx_write_index(fasta, index, verbose=False)["record_count"] == 2
    assert faidx.faidx_read_index(fasta, index_path=index, verbose=False)["record_count"] == 2
    output = tmp_path / "changed.fa"
    faidx.faidx_overwrite_region(fasta, output, "chr1", 2, 3, "TT", verbose=False)
    assert "ATTTACGT" in output.read_text()
    assert "ACGTACGT" in fasta.read_text()


def test_fastx_counts(fasta):
    upstream = pytest.importorskip("pyfastx")
    for length in (1, 5, 100):
        expected = upstream.Fasta(str(fasta)).count(length)
        assert fastx.pyfastx_count_sequences(fasta, min_length=length, verbose=False)["sequence_count"] == expected


@pytest.fixture(params=[("g", "G", 4), ("gene2", "SECOND", 1), ("third", "THIRD", 20)])
def annotation(tmp_path, request):
    gene, name, length = request.param
    source = tmp_path / "annotation.gtf"
    source.write_text(f'chr1\ttest\texon\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; gene_name "{name}"; transcript_id "t";\n')
    return source


def test_genomepy_annotation(annotation, tmp_path):
    upstream = pytest.importorskip("genomepy")
    obj = upstream.Annotation(str(annotation), quiet=True)
    assert genome.genomepy_annotation_attributes(str(annotation), verbose=False)["attributes"] == obj.attributes()
    assert genome.genomepy_annotation_genes(str(annotation), verbose=False)["genes"] == sorted(obj.genes())
    assert genome.genomepy_annotation_lengths(str(annotation), verbose=False)["lengths"] == obj.lengths().to_dict()
    assert genome.genomepy_annotation_dictionary(str(annotation), "gene_id", "gene_name", verbose=False)["mapping"] == obj.gtf_dict("gene_id", "gene_name")
    assert genome.genomepy_list_installed(str(tmp_path), verbose=False)["genomes"] == []
    target = tmp_path / "copy.gtf"
    genome.genomepy_import_local_annotation(str(annotation), str(target), verbose=False)
    assert target.read_bytes() == annotation.read_bytes()


def test_genomepy_sanitize(annotation, fasta, tmp_path):
    pytest.importorskip("genomepy")
    target = tmp_path / "sanitized.gtf"
    result = genome.genomepy_sanitize_annotation(str(annotation), str(fasta), str(target), verbose=False)
    assert result["input_rows"] == result["output_rows"] == 1
    assert target.exists()


@pytest.mark.parametrize("prediction,target", [([1., 2., 3.], [2., 3., 4.]), ([0., 0., 0.], [0., 1., 0.]), ([2.], [1.])])
def test_enformer_loss(prediction, target):
    torch = pytest.importorskip("torch")
    try:
        from enformer_pytorch.modeling_enformer import poisson_loss
    except ImportError as exc:
        pytest.skip(str(exc))
    from dnakit.predictions.enformer_p40_functions import enformer_poisson_loss
    got = enformer_poisson_loss(prediction, target, verbose=False)
    assert got["values"] == pytest.approx(float(poisson_loss(torch.tensor(prediction), torch.tensor(target))))
    json.dumps(got)


def test_genomicranges_restrict(tmp_path):
    upstream = pytest.importorskip("genomicranges")
    import pandas as pd
    source = tmp_path / "r.bed"
    source.write_text("chr1\t0\t8\tfirst\t2\t+\nchr1\t10\t15\tsecond\t3\t+\n")
    for first, last in [(0, 20), (2, 12), (3, 4)]:
        result = gr.genomicranges_restrict_intervals(source, start=first, end=last, verbose=False)
        expected = upstream.GenomicRanges.from_pandas(pd.DataFrame({"seqnames": ["chr1", "chr1"], "starts": [1, 11], "ends": [8, 15], "strand": ["+", "+"]})).restrict(start=first + 1, end=last)
        assert [(r["start"], r["end"]) for r in result["intervals"]] == list(zip([int(v)-1 for v in expected.get_start()], [int(v) for v in expected.get_end()]))


@pytest.mark.parametrize("sequence", ["ACGT" * 8, "N" * 32, "a" * 32])
@pytest.mark.parametrize("kind", ["base", "head", "context", "attention"])
def test_enformer_local_prediction_matches_upstream(tmp_path, sequence, kind):
    torch = pytest.importorskip("torch")
    enformer = pytest.importorskip("enformer_pytorch")
    finetune = pytest.importorskip("enformer_pytorch.finetune")
    from dnakit.predictions import enformer_p40_functions as module
    import numpy as np
    previous_threads = torch.get_num_threads()
    torch.set_num_threads(1)
    try:
        torch.manual_seed(123)
        config = dict(dim=24, depth=1, heads=2, output_heads={"human": 2, "mouse": 1}, target_length=4,
                      num_downsamples=3, dim_divisible_by=4, attn_dim_key=8, dropout_rate=0., attn_dropout=0., pos_dropout=0.)
        model = enformer.Enformer(enformer.EnformerConfig(**config))
        kwargs = {}
        if kind == "base":
            fn = module.enformer_predict_tracks
            kwargs["head"] = "human"
        elif kind == "head":
            model = finetune.HeadAdapterWrapper(enformer=model, num_tracks=2)
            fn = module.enformer_predict_custom_tracks
            kwargs["num_tracks"] = 2
        else:
            cls = finetune.ContextAdapterWrapper if kind == "context" else finetune.ContextAttentionAdapterWrapper
            model = cls(enformer=model, context_dim=3)
            fn = module.enformer_predict_context_tracks if kind == "context" else module.enformer_predict_context_attention
            kwargs.update(context=[[1., 0., 0.], [0., 1., 1.]], context_dim=3)
        model.eval()
        checkpoint, output = tmp_path / "model.pt", tmp_path / "prediction.npy"
        torch.save(model.state_dict(), checkpoint)
        native_args = {"head": "human"} if kind == "base" else {}
        if "context" in kwargs:
            native_args["context"] = torch.tensor(kwargs["context"])
        with torch.no_grad():
            expected = model(enformer.data.str_to_seq_indices([sequence]), **native_args).numpy()
        result = fn(str(checkpoint), [sequence], str(output), config=config, verbose=False, **kwargs)
        np.testing.assert_allclose(np.load(output), expected, rtol=1e-6, atol=1e-6)
        assert result["shape"] == list(expected.shape)
    finally:
        torch.set_num_threads(previous_threads)


@pytest.mark.parametrize("values", [[[1., 2., 3.]], [[3., 1., 0.]], [[0., 0., 0.]]])
def test_enformer_correlations(values):
    torch = pytest.importorskip("torch")
    modeling = pytest.importorskip("enformer_pytorch.modeling_enformer")
    metrics = pytest.importorskip("enformer_pytorch.metrics")
    from dnakit.predictions import enformer_p40_functions as module
    pred = torch.tensor(values).unsqueeze(-1)
    target = torch.tensor([[[0.], [1.], [2.]]])
    expected = modeling.pearson_corr_coef(pred, target).tolist()
    assert module.enformer_pearson_correlation(pred.tolist(), target.tolist(), verbose=False)["values"] == pytest.approx(expected)
    metric = metrics.MeanPearsonCorrCoefPerChannel(1)
    metric.update(pred, target)
    actual = module.enformer_channel_correlation([pred.tolist()], [target.tolist()], verbose=False)["values"]
    expected_channels = metric.compute()
    if torch.isnan(expected_channels).all():
        assert actual is None or actual == [None]
    else:
        assert actual == pytest.approx(expected_channels.tolist())
