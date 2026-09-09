"""Seeded CPU batch loading shared by custom-property training and evaluation."""

from __future__ import annotations

import random
from collections.abc import Iterator, Sequence
from typing import Any, cast


class _PropertyDataset:
    """Only data, never an encoder, API client, or training configuration, is pickled."""

    def __init__(self, sequences: Sequence[str], targets: Any, representations: Any | None) -> None:
        self.sequences = sequences
        self.targets = targets
        self.representations = representations

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> dict[str, Any]:
        item = {
            "row_index": index,
            "sequence": self.sequences[index],
            "targets": self.targets[index].copy(),
        }
        if self.representations is not None:
            item["representations"] = self.representations[index].copy()
        return item


class _BatchSampler:
    """Replace the sample order between fully consumed train/validation/test passes."""

    def __init__(self) -> None:
        self.order: tuple[tuple[int, ...], ...] = ()

    def __iter__(self) -> Iterator[tuple[int, ...]]:
        return iter(self.order)

    def __len__(self) -> int:
        return len(self.order)


def _seed_worker(worker_id: int) -> None:
    import numpy as np
    import torch  # type: ignore[import-not-found]

    del worker_id
    seed = torch.initial_seed() % (2**32)
    random.seed(seed)
    np.random.seed(seed)


class _PropertyBatchLoader:
    """Use one reusable worker pool for all splits; model execution stays in the caller."""

    def __init__(
        self,
        torch: Any,
        sequences: Sequence[str],
        targets: Any,
        *,
        representations: Any | None,
        batch_size: int,
        workers: int,
        seed: int,
    ) -> None:
        self.batch_size = batch_size
        self.sampler = _BatchSampler()
        self.generator = torch.Generator().manual_seed(seed)
        self.loader = torch.utils.data.DataLoader(
            _PropertyDataset(sequences, targets, representations),
            batch_sampler=self.sampler,
            num_workers=workers,
            persistent_workers=workers > 0,
            # Avoid inheriting a CUDA context initialized by the sequence encoder.
            multiprocessing_context="spawn" if workers > 0 else None,
            generator=self.generator,
            worker_init_fn=_seed_worker,
        )

    def batches(
        self, indices: Sequence[int], *, shuffle: bool, seed: int
    ) -> Iterator[dict[str, Any]]:
        import numpy as np

        values = np.asarray(indices, dtype=np.int64)
        if shuffle:
            values = np.random.default_rng(seed).permutation(values)
        self.sampler.order = tuple(
            tuple(int(value) for value in values[start : start + self.batch_size])
            for start in range(0, len(values), self.batch_size)
        )
        self.generator.manual_seed(seed)
        return cast(Iterator[dict[str, Any]], iter(self.loader))

    def close(self) -> None:
        # PyTorch has no public close() for persistent DataLoader workers. Release
        # them explicitly on normal completion, early stopping, or exceptions.
        iterator = getattr(self.loader, "_iterator", None)
        if iterator is not None:
            iterator._shutdown_workers()
            self.loader._iterator = None
