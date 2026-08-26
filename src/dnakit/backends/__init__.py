"""Lazy optional-backend registration without import-time probing."""

from dnakit.backends.external import (
    EXTERNAL_CLI_ADAPTERS,
    BoundedCommandResult,
    ExternalCLIAdapter,
    execute_bounded_command,
)
from dnakit.backends.registry import BackendRegistration, BackendRegistry, backend_registry
from dnakit.core import BackendInfo


def _probe_primer3_cli() -> BackendInfo:
    from dnakit.thermodynamics import probe_primer3

    return probe_primer3()


def _load_primer3_cli() -> object:
    from dnakit.thermodynamics import Primer3CLIAdapter

    return Primer3CLIAdapter()


def _probe_nupack() -> BackendInfo:
    from dnakit.secondary_structure import probe_nupack

    return probe_nupack()


def _load_nupack() -> object:
    from dnakit.secondary_structure import NupackAdapter

    return NupackAdapter()


def _external_registration(adapter: ExternalCLIAdapter) -> BackendRegistration:
    def load() -> object:
        return adapter

    return BackendRegistration(
        adapter.backend_id,
        adapter.probe,
        load,
        source="builtin:external-cli",
    )


def _register_builtin_backends() -> None:
    existing = set(backend_registry)
    for registration in (
        BackendRegistration(
            "primer3-cli",
            _probe_primer3_cli,
            _load_primer3_cli,
            source="builtin:thermodynamics",
        ),
        BackendRegistration(
            "nupack",
            _probe_nupack,
            _load_nupack,
            source="builtin:secondary-structure",
        ),
        *(_external_registration(adapter) for adapter in EXTERNAL_CLI_ADAPTERS),
    ):
        if registration.backend_id not in existing:
            backend_registry.register(registration)


_register_builtin_backends()

__all__ = [
    "EXTERNAL_CLI_ADAPTERS",
    "BackendRegistration",
    "BackendRegistry",
    "BoundedCommandResult",
    "ExternalCLIAdapter",
    "backend_registry",
    "execute_bounded_command",
]
