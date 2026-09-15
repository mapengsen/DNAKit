"""Public, context-managed entry point for text, files and mixed record inputs."""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

from .config import ReadConfig
from .source import RecordSource


def read_inputs(
    value: object,
    *,
    read_config: ReadConfig | None = None,
    sequence_column: str | None = None,
    verbose: bool = True,
) -> RecordSource:
    """Read DNA text, file paths, DNA objects or an iterable through one lazy source.

    File formats are inferred. Uncompressed file limits follow file size; gzip
    has a 20 GB default decoded limit. CSV extras are retained in metadata.
    Use as a context manager when iteration may stop before EOF.
    """
    from dnakit._processing.inputs import iter_inputs

    return RecordSource(
        iter_inputs(value, read_config=read_config, sequence_column=sequence_column)
    )


_register_public_verbose(globals(), ("read_inputs",))
