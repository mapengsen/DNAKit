"""Tests for lazy optional-backend discovery and resolution."""

import pytest

from dnakit.backends import BackendRegistration, BackendRegistry, backend_registry
from dnakit.backends.registry import _entry_point_registration
from dnakit.core import BackendInfo
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


def test_registry_is_lazy_sorted_and_loads_only_on_request() -> None:
    calls: list[str] = []
    registry = BackendRegistry()

    def probe_zeta() -> BackendInfo:
        calls.append("probe-zeta")
        return BackendInfo("zeta", capabilities={"tm"})

    def load_zeta() -> object:
        calls.append("load-zeta")
        return object()

    registry.register(
        BackendRegistration(
            "zeta",
            probe_zeta,
            load_zeta,
        )
    )
    registry.register(BackendRegistration("alpha", lambda: BackendInfo("alpha", available=False)))

    assert tuple(registry) == ("alpha", "zeta")
    assert calls == []
    registry.probe("zeta")
    assert calls == ["probe-zeta"]
    registry.load("zeta", capability="tm")
    assert calls == ["probe-zeta", "probe-zeta", "load-zeta"]


def test_registry_rejects_duplicates_missing_and_capability_mismatch() -> None:
    registry = BackendRegistry()
    registration = BackendRegistration(
        "example.backend", lambda: BackendInfo("example", capabilities={"search"})
    )
    registry.register(registration)

    with pytest.raises(ConfigurationError) as duplicate:
        registry.register(registration)
    assert duplicate.value.code == "DUPLICATE_BACKEND"
    with pytest.raises(BackendUnavailableError) as missing:
        registry.probe("missing")
    assert missing.value.code == "BACKEND_NOT_REGISTERED"
    with pytest.raises(BackendUnavailableError) as capability:
        registry.load("example.backend", capability="tm")
    assert capability.value.code == "BACKEND_CAPABILITY_UNAVAILABLE"


def test_registry_wraps_invalid_probe_results_and_loader_errors() -> None:
    registry = BackendRegistry()
    registry.register(
        BackendRegistration(
            "bad.probe",
            lambda: object(),  # type: ignore[arg-type,return-value]
        )
    )
    registry.register(
        BackendRegistration(
            "bad.loader",
            lambda: BackendInfo("bad-loader"),
            lambda: (_ for _ in ()).throw(RuntimeError("boom")),
        )
    )

    with pytest.raises(BackendExecutionError) as invalid:
        registry.probe("bad.probe")
    assert invalid.value.code == "INVALID_BACKEND_PROBE_RESULT"
    with pytest.raises(BackendExecutionError) as load_error:
        registry.load("bad.loader")
    assert load_error.value.code == "BACKEND_LOAD_FAILED"


@pytest.mark.parametrize("backend_id", ["Bad", "two words", "../escape", ""])
def test_registration_rejects_unsafe_ids(backend_id: str) -> None:
    with pytest.raises(ConfigurationError):
        BackendRegistration(backend_id, lambda: BackendInfo("safe"))


def test_registry_validates_lookup_arguments_and_entry_point_classes_are_instantiated() -> None:
    registry = BackendRegistry()
    with pytest.raises(ConfigurationError):
        registry.unregister([])  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        registry.probe("Bad")

    class Adapter:
        def probe(self) -> BackendInfo:
            return BackendInfo("adapter")

    class Point:
        name = "adapter"
        value = "tests:Adapter"

        @staticmethod
        def load() -> type[Adapter]:
            return Adapter

    registration = _entry_point_registration(Point())  # type: ignore[arg-type]
    registry.register(registration)
    loaded = registry.load("adapter")

    assert isinstance(loaded, Adapter)
    with pytest.raises(ConfigurationError):
        registry.load("adapter", capability="")


def test_builtin_primer3_cli_adapter_is_registered_without_eager_probe() -> None:
    registrations = {item.backend_id: item for item in backend_registry.registrations()}

    assert {"nupack", "primer3-cli"} <= registrations.keys()
    assert "primer3-py" not in registrations
    assert registrations["primer3-cli"].source == "builtin:thermodynamics"
    info = backend_registry.probe("primer3-cli")
    assert info.name == "primer3-cli"
    assert info.metadata["automatic_path_search"] is False
