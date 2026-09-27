"""Offline transaction recovery, conflict detection and source-boundary checks."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from copy import deepcopy
from pathlib import Path
from typing import Any, cast
from unittest.mock import Mock

import pytest

import dnakit.generated as generated
import dnakit.tools.development_support as support
from dnakit.exceptions import ConfigurationError
from dnakit.tools.development_support import DevelopmentJob, rollback_tool_change


@pytest.fixture
def checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "checkout"
    tools = root / "src/dnakit/tools"
    tools.mkdir(parents=True)
    (root / "pyproject.toml").write_text('[project]\nname = "dnakit"\n', encoding="utf-8")
    (tools / "agent.py").write_text("# source checkout\n", encoding="utf-8")
    source = tools / "development_support.py"
    source.write_text("# transaction support\n", encoding="utf-8")
    monkeypatch.setattr(support, "__file__", str(source))
    monkeypatch.setattr(generated, "_refresh", Mock())
    return root


def _action(revision: int = 1) -> dict[str, Any]:
    return {
        "reason": "Fill a verified offline capability gap",
        "name": "transaction_probe",
        "target": "dnakit.generated.transaction_probe",
        "contract": {"input": "DNA text", "output": "length"},
        "tests": ["assert length('ACGT') == 4"],
        "scenarios": [{"sequence": "ACGT"}],
        "source": f"implementation revision {revision}",
    }


def _job(root: Path, revision: int = 1, *, action: dict[str, Any] | None = None) -> DevelopmentJob:
    return DevelopmentJob(
        root,
        "generate_tool",
        "dnakit.generated.transaction_probe",
        _action(revision) if action is None else action,
        False,
    )


def _files(root: Path) -> tuple[Path, Path]:
    original = root / "src/dnakit/original.py"
    added = root / "src/dnakit/added.py"
    original.write_bytes(b"original implementation\n")
    return original, added


def _activate(root: Path) -> tuple[DevelopmentJob, Path, Path]:
    original, added = _files(root)
    job = _job(root)
    job.watch([original, added])
    job.activate({original: b"new implementation\n", added: b"added implementation\n"})
    return job, original, added


def _state(job: DevelopmentJob) -> dict[str, Any]:
    return cast(
        dict[str, Any], json.loads((job.directory / "state.json").read_text(encoding="utf-8"))
    )


def test_source_checkout_is_accepted(checkout: Path) -> None:
    assert support._source_root() == checkout


def test_installed_module_is_rejected_even_inside_a_directory_with_checkout_markers(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    installed = checkout / "site-packages/dnakit/tools/development_support.py"
    installed.parent.mkdir(parents=True)
    installed.write_text("# installed wheel\n", encoding="utf-8")
    monkeypatch.setattr(support, "__file__", str(installed))

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_SOURCE_REQUIRED"):
        support._source_root()


@pytest.mark.parametrize("field", ["tests", "scenarios", "contract", "target", "name"])
def test_validation_cases_and_contract_cannot_change_between_revisions(
    checkout: Path,
    field: str,
) -> None:
    first = _job(checkout)
    first.reject("Revise the implementation")
    revised = _action(2)
    revised[field] = "weakened validation"

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_FIXED_SCENARIOS"):
        _job(checkout, action=revised)
    assert _state(first)["attempt"] == 1


def test_three_revisions_are_enforced_and_all_candidates_are_preserved(checkout: Path) -> None:
    for revision in range(1, 4):
        job = _job(checkout, revision)
        assert job.attempt == revision
        job.reject("Review rejected this implementation")

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_ROUND_LIMIT"):
        _job(checkout, 4)
    assert sorted(path.name for path in job.directory.glob("candidate-*.json")) == [
        "candidate-1.json",
        "candidate-2.json",
        "candidate-3.json",
    ]


def test_duplicate_activated_receipt_replays_without_writes_or_new_attempt(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    job, original, added = _activate(checkout)
    before = (original.read_bytes(), added.read_bytes(), _state(job))
    writes = Mock(wraps=support._atomic_write)
    monkeypatch.setattr(support, "_atomic_write", writes)

    replay = _job(checkout)

    assert replay.replay == job.state["result"] and replay.attempt == 1
    assert replay.replay and replay.replay["status"] == "activated"
    writes.assert_not_called()
    assert before == (original.read_bytes(), added.read_bytes(), _state(job))


def test_interrupted_validation_does_not_rerun_external_effects(checkout: Path) -> None:
    job = _job(checkout)
    job.phase("executing")

    restored = _job(checkout)

    assert restored.attempt == 1 and restored.replay
    assert restored.replay["code"] == "AGENT_DEVELOPMENT_VALIDATION_UNCERTAIN"
    assert len(list(job.directory.glob("candidate-*.json"))) == 1


def test_activation_exception_restores_all_existing_and_new_targets(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    atomic = support._atomic_write
    failed = False

    def fail_second(path: Path, data: bytes) -> None:
        nonlocal failed
        if path == added and not failed:
            failed = True
            raise OSError("Simulated disk failure")
        atomic(path, data)

    monkeypatch.setattr(support, "_atomic_write", fail_second)
    with pytest.raises(OSError, match="Simulated disk failure"):
        job.activate({original: b"new implementation\n", added: b"added implementation\n"})

    assert original.read_bytes() == b"original implementation\n" and not added.exists()
    assert _state(job)["phase"] == "rejected"


def test_interrupted_partial_activation_rolls_back_before_replay(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    atomic = support._atomic_write

    def interrupt_second(path: Path, data: bytes) -> None:
        if path == added:
            raise KeyboardInterrupt
        atomic(path, data)

    monkeypatch.setattr(support, "_atomic_write", interrupt_second)
    with pytest.raises(KeyboardInterrupt):
        job.activate({original: b"new implementation\n", added: b"added implementation\n"})
    assert original.read_bytes() == b"new implementation\n" and not added.exists()
    assert _state(job)["phase"] == "activating"
    monkeypatch.setattr(support, "_atomic_write", atomic)

    recovered = _job(checkout)

    assert original.read_bytes() == b"original implementation\n" and not added.exists()
    assert recovered.replay and recovered.replay["status"] == "rejected"
    assert recovered.attempt == 1


def test_interrupted_receipt_after_all_writes_is_reconciled_as_activated(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    monkeypatch.setattr(job, "complete", Mock(side_effect=KeyboardInterrupt()))
    with pytest.raises(KeyboardInterrupt):
        job.activate({original: b"new implementation\n", added: b"added implementation\n"})
    assert _state(job)["phase"] == "activating"

    recovered = _job(checkout)

    assert recovered.replay and recovered.replay["status"] == "activated"
    assert recovered.attempt == 1
    assert original.read_bytes() == b"new implementation\n"
    assert added.read_bytes() == b"added implementation\n"


def test_optimistic_conflict_never_overwrites_user_changes(checkout: Path) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    original.write_bytes(b"user edit\n")

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFLICT"):
        job.activate({original: b"new implementation\n", added: b"added implementation\n"})

    assert original.read_bytes() == b"user edit\n" and not added.exists()


def test_retargeted_symlink_cannot_overwrite_an_unwatched_source_file(checkout: Path) -> None:
    original, added = _files(checkout)
    added.write_bytes(b"unrelated user file\n")
    alias = checkout / "src/dnakit/alias.py"
    alias.symlink_to(original)
    job = _job(checkout)
    job.watch([alias])
    alias.unlink()
    alias.symlink_to(added)

    with pytest.raises(ConfigurationError, match=r"AGENT_DEVELOPMENT_(CONFLICT|PATH)"):
        job.activate({alias: b"new implementation\n"})

    assert added.read_bytes() == b"unrelated user file\n"


def test_changed_watched_dependency_blocks_activation_even_if_it_is_not_rewritten(
    checkout: Path,
) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    original.write_bytes(b"user changed validation dependency\n")

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFLICT"):
        job.activate({added: b"candidate validated against obsolete implementation\n"})

    assert not added.exists()
    assert original.read_bytes() == b"user changed validation dependency\n"


def test_activation_outside_checkout_is_rejected_without_writing(checkout: Path) -> None:
    outside = checkout.parent / "outside.py"
    job = _job(checkout)

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_PATH"):
        job.activate({outside: b"unexpected mutation\n"})
    assert not outside.exists()


def test_activation_uses_journaled_targets_when_a_symlink_changes_after_preflight(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original, unrelated = _files(checkout)
    unrelated.write_bytes(b"unrelated user file\n")
    alias = checkout / "src/dnakit/alias.py"
    alias.symlink_to(original)
    job = _job(checkout)
    job.watch([alias])
    phase = job.phase

    def retarget_after_journal(name: str, **data: Any) -> None:
        phase(name, **data)
        if name == "activating":
            alias.unlink()
            alias.symlink_to(unrelated)

    monkeypatch.setattr(job, "phase", retarget_after_journal)
    # Either reject the changed path or write only the journaled target.
    with suppress(ConfigurationError):
        job.activate({alias: b"new implementation\n"})

    assert unrelated.read_bytes() == b"unrelated user file\n"


def test_source_write_lock_excludes_a_second_writer_and_releases_after_interrupt(
    checkout: Path,
) -> None:
    with pytest.raises(KeyboardInterrupt), support._write_lock(checkout):
        with ThreadPoolExecutor(max_workers=1) as workers:
            future = workers.submit(_job, checkout)
            with pytest.raises(OSError):
                future.result(timeout=5)
        raise KeyboardInterrupt
    assert _job(checkout).attempt == 1


def test_conflict_after_journaling_does_not_rollback_over_a_new_user_edit(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original, added = _files(checkout)
    job = _job(checkout)
    job.watch([original, added])
    phase = job.phase

    def user_edits_after_journal(name: str, **data: Any) -> None:
        phase(name, **data)
        if name == "activating":
            original.write_bytes(b"user edit after preflight\n")

    monkeypatch.setattr(job, "phase", user_edits_after_journal)
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFLICT"):
        job.activate({original: b"new implementation\n", added: b"added implementation\n"})

    assert original.read_bytes() == b"user edit after preflight\n" and not added.exists()


def test_stale_job_cannot_replace_a_newer_candidate_receipt(checkout: Path) -> None:
    first = _job(checkout)
    first.reject("First revision failed review")
    second = _job(checkout, 2)
    second.reject("Second revision failed review")
    current = deepcopy(_state(second))

    with pytest.raises(ConfigurationError, match=r"AGENT_DEVELOPMENT_(CONFLICT|STALE)"):
        first.reject("Late response from superseded worker")

    assert _state(second) == current


@pytest.mark.parametrize("phase", ["tested", "validated", "baseline_completed"])
def test_third_revision_survives_repeated_recovery_before_next_phase(
    checkout: Path, phase: str
) -> None:
    for revision in (1, 2):
        _job(checkout, revision).reject("Revise the implementation")
    third = _job(checkout, 3)
    third.phase(phase, validation={"passed": True})

    # Constructing a resumed job may itself be followed by another interruption.
    first_recovery = _job(checkout, 3)
    second_recovery = _job(checkout, 3)

    assert first_recovery.resuming and second_recovery.resuming
    assert first_recovery.attempt == second_recovery.attempt == 3
    assert _state(second_recovery)["phase"] == phase
    assert len(list(second_recovery.directory.glob("candidate-*.json"))) == 3


def test_rollback_restores_previous_bytes_removes_new_files_and_is_idempotent(
    checkout: Path,
) -> None:
    job, original, added = _activate(checkout)

    result = rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)
    original.write_bytes(b"user edits after rollback\n")
    again = rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert result["status"] == again["status"] == "rolled_back"
    assert original.read_bytes() == b"user edits after rollback\n" and not added.exists()
    assert _state(job)["phase"] == "rolled_back"


@pytest.mark.parametrize("failed_cache", ["registry", "catalog"])
def test_rollback_cache_refresh_failure_keeps_recoverable_journal(
    checkout: Path, monkeypatch: pytest.MonkeyPatch, failed_cache: str
) -> None:
    from dnakit.skills.catalog import _catalog
    from dnakit.tools.registry import default_tool_registry

    job, original, added = _activate(checkout)
    cache = default_tool_registry if failed_cache == "registry" else _catalog
    clear = Mock(side_effect=OSError("cache refresh interrupted"))
    monkeypatch.setattr(cache, "cache_clear", clear)

    with pytest.raises(OSError, match="cache refresh interrupted"):
        rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert _state(job)["phase"] == "rolling_back"
    assert original.read_bytes() == b"original implementation\n" and not added.exists()
    clear.side_effect = None
    result = rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert result["status"] == "rolled_back"
    assert _state(job)["phase"] == "rolled_back" and clear.call_count == 2


def test_already_rolled_back_refreshes_caches_without_rewriting_later_edits(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from dnakit.skills.catalog import _catalog
    from dnakit.tools.registry import default_tool_registry

    job, original, added = _activate(checkout)
    rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)
    original.write_bytes(b"user edits after rollback\n")
    registry_clear = Mock()
    catalog_clear = Mock()
    monkeypatch.setattr(default_tool_registry, "cache_clear", registry_clear)
    monkeypatch.setattr(_catalog, "cache_clear", catalog_clear)
    monkeypatch.setattr(
        support, "_refresh_rolled_back_code", Mock(side_effect=AssertionError("already restored"))
    )

    result = rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert result["status"] == "rolled_back"
    assert original.read_bytes() == b"user edits after rollback\n" and not added.exists()
    registry_clear.assert_called_once_with()
    catalog_clear.assert_called_once_with()


def test_post_activation_error_does_not_destroy_the_rollback_receipt(checkout: Path) -> None:
    job, original, added = _activate(checkout)
    try:
        result = job.reject("Post-activation reporting failed")
    except ConfigurationError:
        result = None

    assert _state(job)["phase"] == "activated"
    assert result is None or result["status"] == "activated"
    assert original.read_bytes() == b"new implementation\n" and added.exists()
    rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)
    assert original.read_bytes() == b"original implementation\n" and not added.exists()


def test_rollback_conflict_preserves_all_user_and_activated_files(checkout: Path) -> None:
    job, original, added = _activate(checkout)
    original.write_bytes(b"user edit after activation\n")

    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFLICT"):
        rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert original.read_bytes() == b"user edit after activation\n"
    assert added.read_bytes() == b"added implementation\n"
    assert _state(job)["phase"] == "activated"


def test_interrupted_rollback_resumes_from_inverse_journal(
    checkout: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    job, original, added = _activate(checkout)
    atomic = support._atomic_write

    def interrupt_restore(path: Path, data: bytes) -> None:
        if path == original:
            raise KeyboardInterrupt
        atomic(path, data)

    monkeypatch.setattr(support, "_atomic_write", interrupt_restore)
    with pytest.raises(KeyboardInterrupt):
        rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)
    assert _state(job)["phase"] == "rolling_back" and not added.exists()
    monkeypatch.setattr(support, "_atomic_write", atomic)

    result = rollback_tool_change(job.change_id, allow_side_effects=True, verbose=False)

    assert result["status"] == "rolled_back"
    assert original.read_bytes() == b"original implementation\n" and not added.exists()


def test_rollback_requires_explicit_permission(checkout: Path) -> None:
    job, original, added = _activate(checkout)

    with pytest.raises(ConfigurationError, match="AGENT_TOOL_CONFIRMATION_REQUIRED"):
        rollback_tool_change(job.change_id, verbose=False)

    assert original.read_bytes() == b"new implementation\n" and added.exists()
