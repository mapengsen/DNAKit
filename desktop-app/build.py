"""Build, test and package the downloadable DNAKit desktop app for the current system.

    python desktop-app/build.py --install

``--install`` first installs DNAKit from this checkout with the GUI dependencies and PyInstaller
into the current Python environment (use a fresh virtual environment). The steps are:

1. PyInstaller bundles Python, DNAKit and the interface (``dnakit_desktop.spec``).
2. The packaged app runs ``smoke_test.py`` and opens its window once (``--skip-tests`` omits this).
3. The app is packed into the file users download, written to ``desktop-app/dist/``:
   ``DNAKit-windows-x64.zip``, ``DNAKit-macos-arm64.dmg`` / ``DNAKit-macos-x64.dmg`` or
   ``DNAKit-linux-x64.tar.gz``.

Each system builds only its own app, so the three downloads come from three machines; the
``Desktop app`` GitHub Actions workflow runs this script on all of them.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DIST = HERE / "dist"
WORK = HERE / "build"
SYSTEM = {"win32": "windows", "darwin": "macos"}.get(sys.platform, "linux")
ARCH = {"amd64": "x64", "x86_64": "x64", "aarch64": "arm64"}.get(
    platform.machine().lower(), platform.machine().lower()
)
BUILD_REQUIREMENTS = ["httpx[socks]", "pyinstaller>=6.10,<7", "Pillow"]
QT_WHEELS = "PyQt6,PyQt6-Qt6,PyQt6-sip,PyQt6-WebEngine,PyQt6-WebEngine-Qt6"


def run(command: list[str], **kwargs: object) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, **kwargs)  # type: ignore[call-overload]


def install() -> None:
    """Install DNAKit (not editable, so the Skills are bundled) and the build tools."""

    extra = "gui-qt" if SYSTEM == "linux" else "gui"
    pip = [sys.executable, "-m", "pip", "install", "--upgrade"]
    if SYSTEM == "linux":  # never compile Qt from source on an older Linux
        pip += ["--only-binary", QT_WHEELS]
    run([*pip, f"{ROOT}[{extra}]", *BUILD_REQUIREMENTS])


def bundle() -> Path:
    """Run PyInstaller and return the packaged program."""

    for folder in (DIST, WORK):
        shutil.rmtree(folder, ignore_errors=True)
    run(
        [
            sys.executable,
            "-m",
            "PyInstaller",
            "--noconfirm",
            "--distpath",
            str(DIST),
            "--workpath",
            str(WORK),
            str(HERE / "dnakit_desktop.spec"),
        ]
    )
    if SYSTEM == "macos":
        return DIST / "DNAKit.app" / "Contents" / "MacOS" / "DNAKit"
    return DIST / "DNAKit" / ("DNAKit.exe" if SYSTEM == "windows" else "DNAKit")


def fail(message: str, detail: str = "") -> None:
    """Stop the build; on GitHub Actions also show the reason as an error annotation."""

    if os.environ.get("GITHUB_ACTIONS"):
        text = f"{message}\n{detail[-6000:]}".strip()
        escaped = text.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
        print(f"::error title=DNAKit desktop app ({SYSTEM}-{ARCH})::{escaped}", flush=True)
    sys.exit(f"{message}\n{detail}".strip())


def run_app(
    program: Path, args: list[str], env: dict[str, str], capture: Path | None = None
) -> int:
    """Start the packaged app and return its exit code.

    Without ``capture`` the app starts the way a double click starts it. With ``capture`` its
    messages are written to that file, which is how a failed start is diagnosed.
    """

    print("+", program.name, " ".join(args), flush=True)
    output = None if capture is None else capture.open("w", encoding="utf-8")
    if SYSTEM == "windows" and capture is None:
        # Start-Process gives the app no console and no standard streams.
        quoted = ",".join("'\"" + a.replace("'", "''") + "\"'" for a in args)
        script = (
            f"$p = Start-Process -FilePath '{program}' -ArgumentList @({quoted}) "
            "-Wait -PassThru; exit $p.ExitCode"
        )
        command = ["powershell", "-NoProfile", "-NonInteractive", "-Command", script]
    else:
        command = [str(program), *args]
        if SYSTEM == "linux" and not (env.get("DISPLAY") or env.get("WAYLAND_DISPLAY")):
            if shutil.which("xvfb-run") is None:
                sys.exit("No display: install xvfb or pass --skip-tests.")
            command = ["xvfb-run", "-a", "-s", "-screen 0 1400x1000x24", *command]
    try:
        done = subprocess.run(
            command,
            env=env,
            cwd=WORK,
            stdin=subprocess.DEVNULL if output else None,
            stdout=output,
            stderr=subprocess.STDOUT if output else None,
            timeout=900,
            check=False,
        )
        return done.returncode
    except subprocess.TimeoutExpired:
        return 124
    finally:
        if output is not None:
            output.close()


def test(program: Path) -> None:
    """Run the in-app checks, then open the window and check that the page renders."""

    env = dict(os.environ)
    if SYSTEM == "linux" and os.geteuid() == 0:
        env.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--no-sandbox")
    report = DIST / "smoke-report.json"
    log = DIST / "test-output.log"
    smoke = [str(HERE / "smoke_test.py"), str(report)]
    code = run_app(program, smoke, env)
    checks = json.loads(report.read_text("utf-8")) if report.exists() else {}
    if code != 0 or not checks.get("ok"):
        errors = [
            f"{name}: {value['error']}"
            for name, value in checks.items()
            if isinstance(value, dict) and not value.get("ok")
        ]
        if not errors:  # the app stopped before it wrote the report
            run_app(program, smoke, env, capture=log)
            errors = [log.read_text("utf-8", errors="replace")]
        fail(f"The packaged app failed its checks (exit code {code}).", "\n".join(errors))
    print(json.dumps(checks, indent=2))
    screenshot = DIST / f"window-{SYSTEM}-{ARCH}.png"
    window = ["gui", "--workdir", str(WORK / "gui-work"), "--self-test", "180"]
    code = run_app(program, [*window, "--screenshot", str(screenshot)], env)
    if code != 0:
        run_app(program, window, env, capture=log)
        detail = log.read_text("utf-8", errors="replace")
        fail(f"The packaged app did not open its window (exit code {code}).", detail)


def package() -> Path:
    """Pack the app into the single file that users download."""

    stem = DIST / f"DNAKit-{SYSTEM}-{ARCH}"
    if SYSTEM == "windows":
        return Path(shutil.make_archive(str(stem), "zip", DIST, "DNAKit"))
    if SYSTEM == "linux":
        target = stem.with_name(stem.name + ".tar.gz")
        with tarfile.open(target, "w:gz") as archive:
            archive.add(DIST / "DNAKit", arcname="DNAKit")
        return target
    target = stem.with_name(stem.name + ".dmg")
    staging = WORK / "dmg"
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    run(["ditto", str(DIST / "DNAKit.app"), str(staging / "DNAKit.app")])
    os.symlink("/Applications", staging / "Applications")
    create = ["hdiutil", "create", "-volname", "DNAKit", "-srcfolder", str(staging)]
    for attempt in range(5):  # hdiutil sometimes reports "Resource busy" on build machines
        try:
            run([*create, "-ov", "-format", "UDZO", str(target)])
            break
        except subprocess.CalledProcessError:
            if attempt == 4:
                raise
            time.sleep(10)
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n", maxsplit=1)[0])
    parser.add_argument("--install", action="store_true", help="pip install DNAKit and build tools")
    parser.add_argument("--skip-tests", action="store_true", help="do not run the packaged app")
    options = parser.parse_args()
    if options.install:
        install()
    program = bundle()
    if not options.skip_tests:
        test(program)
    target = package()
    print(f"Built {target} ({target.stat().st_size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
