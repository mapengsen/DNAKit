"""Check a packaged DNAKit desktop app from the inside.

Run it with the packaged app, which executes it in the bundled Python:

    desktop-app/dist/DNAKit/DNAKit desktop-app/smoke_test.py report.json

It exercises what differs from a pip installation: bundled data, the function registry, and the
helper processes the interface starts through ``sys.executable``. The JSON report is written to
the given path because a windowed app has no terminal; the exit code is 0 only if all checks pass.
"""

from __future__ import annotations

import json
import multiprocessing
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Any

MIN_FUNCTIONS = 5000
MIN_SKILLS = 500


def check_registry() -> Any:
    from dnakit.tools.registry import default_tool_registry

    stats = default_tool_registry(verbose=False).stats()
    assert int(str(stats["total"])) >= MIN_FUNCTIONS, stats
    return stats


def check_skills() -> Any:
    from dnakit.skills import skill_catalog_stats

    count = int(skill_catalog_stats()["count"])
    assert count >= MIN_SKILLS, count
    return count


def check_interface_builds() -> Any:
    from dnakit.gui import GuiSettings, build_app

    with tempfile.TemporaryDirectory() as folder:
        build_app(GuiSettings.local(folder))
    return "built"


def check_cli_subprocess() -> Any:
    from dnakit.gui._backend import CLI_BOOTSTRAP

    done = subprocess.run(
        [sys.executable, "-c", CLI_BOOTSTRAP, "describe", "ACGTACGTAA"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        timeout=300,
        check=False,
    )
    assert done.returncode == 0, done.stderr[-2000:]
    return json.loads(done.stdout)["values"]["symbol_length"]


def check_tool_worker() -> Any:
    from dnakit.gui._worker import ToolWorker

    worker = ToolWorker()
    try:
        with tempfile.TemporaryDirectory() as folder:
            ok, result = worker.call(
                "dnakit.ops.reverse_complement",
                {"value": "AACGT"},
                Path(folder),
                timeout=300,
            )
    finally:
        worker.stop()
    assert ok and "ACGTT" in json.dumps(result), result
    return "ACGTT"


def check_process_pool() -> Any:
    with multiprocessing.Pool(2) as pool:
        values = pool.map(abs, [-1, -2, -3])
    assert values == [1, 2, 3], values
    return values


CHECKS = [
    check_registry,
    check_skills,
    check_interface_builds,
    check_cli_subprocess,
    check_tool_worker,
    check_process_pool,
]


def main() -> int:
    report: dict[str, Any] = {"executable": sys.executable, "frozen": hasattr(sys, "frozen")}
    failed = False
    for check in CHECKS:
        try:
            report[check.__name__] = {"ok": True, "value": check()}
        except Exception:
            failed = True
            report[check.__name__] = {"ok": False, "error": traceback.format_exc()}
    report["ok"] = not failed
    text = json.dumps(report, indent=2, default=str)
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(text, encoding="utf-8")
    print(text)
    return 1 if failed else 0


if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())
