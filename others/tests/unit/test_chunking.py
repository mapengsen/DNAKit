from io import StringIO

import pytest

from dnakit import (
    ChunkingConfig,
    ChunkProgress,
    DNARecord,
    DNASequence,
    LengthCurriculum,
    SequenceChunk,
    iter_fasta_chunks,
    iter_sequence_chunks,
    make_length_curriculum,
)
from dnakit.exceptions import InputFormatError


def _chunks(sequence: str, config: ChunkingConfig) -> list[SequenceChunk]:
    return list(iter_sequence_chunks(DNASequence(sequence), config=config))


def test_fixed_chunks_default_to_non_overlapping_train_windows() -> None:
    chunks = _chunks("ACGTACGTAC", ChunkingConfig(length=4))

    assert [chunk.sequence.to_string() for chunk in chunks] == ["ACGT", "ACGT"]
    assert [(chunk.source_start, chunk.source_end) for chunk in chunks] == [(0, 4), (4, 8)]
    assert all(chunk.split == "train" for chunk in chunks)


def test_fixed_chunks_can_keep_a_partial_tail() -> None:
    chunks = _chunks("ACGTACGTAC", ChunkingConfig(length=4, include_partial=True))

    assert [chunk.sequence.to_string() for chunk in chunks] == ["ACGT", "ACGT", "AC"]


def test_sliding_chunks_use_the_requested_step() -> None:
    chunks = _chunks(
        "ACGTACGTAC",
        ChunkingConfig(strategy="sliding", length=4, step=2),
    )

    assert [(chunk.source_start, chunk.source_end) for chunk in chunks] == [
        (0, 4),
        (2, 6),
        (4, 8),
        (6, 10),
    ]


def test_random_chunks_are_reproducible_and_within_range() -> None:
    config = ChunkingConfig(
        strategy="random",
        min_length=2,
        max_length=5,
        num_samples=8,
        seed=17,
    )

    first = _chunks("ACGTACGTACGT", config)
    second = _chunks("ACGTACGTACGT", config)

    assert [(item.source_start, item.source_end) for item in first] == [
        (item.source_start, item.source_end) for item in second
    ]
    assert all(2 <= item.length <= 5 for item in first)


def test_multiscale_chunks_mark_the_length_level() -> None:
    chunks = _chunks(
        "ACGTACGTACGT",
        ChunkingConfig(strategy="multiscale", lengths=(4, 8)),
    )

    assert [(chunk.level_index, chunk.source_start, chunk.source_end) for chunk in chunks] == [
        (0, 0, 4),
        (0, 4, 8),
        (0, 8, 12),
        (1, 0, 8),
    ]


def test_curriculum_is_a_short_to_long_schedule() -> None:
    curriculum = make_length_curriculum((4, 8, 16), stage_steps=(100, 200, 300))
    config = curriculum.to_config()

    assert isinstance(curriculum, LengthCurriculum)
    assert [(stage.length, stage.training_steps) for stage in curriculum.stages] == [
        (4, 100),
        (8, 200),
        (16, 300),
    ]
    assert config.strategy == "curriculum"
    assert config.lengths == (4, 8, 16)
    assert config.stage_steps == (100, 200, 300)


def test_fasta_without_bed_assigns_every_record_to_train() -> None:
    fasta = StringIO(">chr1\nACGTACGT\n>chr2\nTTTT\n")

    chunks = list(
        iter_fasta_chunks(
            fasta,
            config=ChunkingConfig(length=4),
        )
    )

    assert [(chunk.source_id, chunk.split, chunk.sequence.to_string()) for chunk in chunks] == [
        ("chr1", "train", "ACGT"),
        ("chr1", "train", "ACGT"),
        ("chr2", "train", "TTTT"),
    ]


def test_fasta_bed_intervals_select_regions_and_use_bed_names_as_splits() -> None:
    fasta = StringIO(">chr1\nAACCGGTTAACCGGTT\n")
    bed = StringIO("chr1\t0\t8\tvalid\nchr1\t8\t16\ttest\n")

    chunks = list(
        iter_fasta_chunks(
            fasta,
            config=ChunkingConfig(length=4),
            bed=bed,
        )
    )

    assert [(chunk.split, chunk.source_start, chunk.source_end) for chunk in chunks] == [
        ("valid", 0, 4),
        ("valid", 4, 8),
        ("test", 8, 12),
        ("test", 12, 16),
    ]
    assert chunks[0].to_dict()["source_start"] == 0


def test_progress_callback_receives_lifecycle_events() -> None:
    events: list[ChunkProgress] = []
    list(
        iter_sequence_chunks(
            DNARecord(DNASequence("ACGTACGT"), "chr1"),
            config=ChunkingConfig(length=4),
            progress=events.append,
        )
    )

    assert [event.status for event in events] == ["started", "yielded", "yielded", "completed"]
    assert events[0].total == 2
    assert events[-1].processed == 2


def test_fasta_bed_missing_sequence_is_reported() -> None:
    fasta = StringIO(">chr1\nACGT\n")
    bed = StringIO("chr2\t0\t4\ttest\n")

    with pytest.raises(InputFormatError, match="BED contains sequence IDs"):
        list(iter_fasta_chunks(fasta, config=ChunkingConfig(length=2), bed=bed))
