"""Entry point of the downloadable DNAKit desktop app (built with PyInstaller).

Started without arguments (a double click), the app opens the DNAKit window and keeps uploaded
inputs and results in ``~/DNAKit``. With arguments it behaves like the ``dnakit`` command line,
for example ``DNAKit describe ACGT``.

DNAKit starts helper processes with ``sys.executable -m module``, ``sys.executable -c code`` and
``sys.executable script.py``. In the packaged app ``sys.executable`` is the app itself, so this
entry point also accepts those three interpreter forms and runs them in the bundled Python.
"""

from __future__ import annotations

import multiprocessing
import os
import runpy
import sys
from pathlib import Path

INTERPRETER_FLAGS = {"-u", "-I", "-E", "-s", "-S", "-B", "-O", "-OO", "-q"}
LOCAL_HOSTS = ("127.0.0.1", "localhost", "::1")


def _repair_standard_streams() -> bool:
    """Give a windowed app valid stdin, stdout and stderr; return whether any was missing."""

    repaired = False
    for fd, mode in ((0, os.O_RDONLY), (1, os.O_WRONLY), (2, os.O_WRONLY)):
        try:
            os.fstat(fd)
        except OSError:
            null = os.open(os.devnull, mode)
            if null != fd:
                os.dup2(null, fd)
                os.close(null)
            repaired = True
    # These streams stay open for the life of the process.
    text = {"encoding": "utf-8", "errors": "replace", "closefd": False}
    if sys.stdin is None:
        sys.stdin = open(0, **text)  # noqa: SIM115
    if sys.stdout is None:
        sys.stdout = open(1, "w", buffering=1, **text)  # noqa: SIM115
    if sys.stderr is None:
        sys.stderr = open(2, "w", buffering=1, **text)  # noqa: SIM115
    return repaired


def _bypass_proxy_for_local_pages() -> None:
    """The window loads a page served on this computer; never route it through a proxy."""

    current = os.environ.get("NO_PROXY") or os.environ.get("no_proxy") or ""
    hosts = [h.strip() for h in current.split(",") if h.strip()]
    hosts += [h for h in LOCAL_HOSTS if h not in hosts]
    os.environ["NO_PROXY"] = os.environ["no_proxy"] = ",".join(hosts)


def _open_window(windowed: bool) -> None:
    workdir = Path.home() / "DNAKit"
    workdir.mkdir(parents=True, exist_ok=True)
    if windowed:  # no terminal to print to: keep messages in a log file
        log = os.open(workdir / "desktop.log", os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(log, 1)
        os.dup2(log, 2)
        os.close(log)
    os.chdir(workdir)
    _run_cli(["gui", "--workdir", str(workdir)])


def _run_cli(args: list[str]) -> None:
    from dnakit.cli.app import run

    sys.argv = ["dnakit", *args]
    run()


def main() -> None:
    multiprocessing.freeze_support()
    windowed = _repair_standard_streams()
    _bypass_proxy_for_local_pages()
    args = sys.argv[1:]
    while args and args[0] in INTERPRETER_FLAGS:
        args.pop(0)
    if not args or args[0].startswith("-psn_"):  # macOS Finder may pass a process serial number
        _open_window(windowed)
    elif args[0] == "-c" and len(args) > 1:
        sys.argv = ["-c", *args[2:]]
        exec(args[1], {"__name__": "__main__", "__builtins__": __builtins__})
    elif args[0] == "-m" and len(args) > 1:
        sys.argv = args[1:]
        runpy.run_module(args[1], run_name="__main__", alter_sys=True)
    elif args[0].endswith((".py", ".pyw")) and Path(args[0]).is_file():
        sys.argv = args
        sys.path.insert(0, str(Path(args[0]).resolve().parent))
        runpy.run_path(args[0], run_name="__main__")
    else:
        _run_cli(args)


if __name__ == "__main__":
    main()
