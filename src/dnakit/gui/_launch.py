"""Start the DNAKit graphical interface in a desktop window, a browser, or as a web server."""

from __future__ import annotations

import os
import sys
import threading
import time
import urllib.request
import webbrowser
from pathlib import Path
from typing import TYPE_CHECKING, Any

from dnakit.exceptions import BackendUnavailableError
from dnakit.gui._settings import GuiSettings

if TYPE_CHECKING:
    import gradio as gr

WINDOW_SIZE = (1280, 900)
MIN_WINDOW_SIZE = (900, 600)


def _require_gradio() -> None:
    os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
    try:
        import gradio  # noqa: F401
    except ImportError as exc:
        raise BackendUnavailableError(
            "The DNAKit graphical interface requires the optional GUI dependencies.",
            code="GUI_BACKEND_UNAVAILABLE",
            hint='Install them with: pip install "dnakit[gui]"',
        ) from exc


def build_app(settings: GuiSettings | None = None) -> gr.Blocks:
    """Build the DNAKit interface as a Gradio ``Blocks`` app without starting a server.

    ``settings`` defaults to ``GuiSettings.local()``. Use ``GuiSettings.public_demo()`` for an
    internet-facing deployment.
    """

    _require_gradio()
    from dnakit.gui._app import build_app as _build

    return _build(GuiSettings.local() if settings is None else settings)


def _has_display() -> bool:
    if not sys.platform.startswith("linux"):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def _page_ready(url: str, timeout: float) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=10) as response:
                if response.status == 200 and b"gradio_config" in response.read():
                    return True
        except OSError:
            pass
        time.sleep(0.5)
    return False


def _window_self_test(
    window: Any, timeout: float, screenshot: Path | None, result: list[bool]
) -> None:
    """Wait until the page shows its header, optionally capture the screen, then close."""

    deadline = time.monotonic() + timeout
    ok = False
    if window.events.loaded.wait(timeout):
        while time.monotonic() < deadline:
            try:
                if window.evaluate_js("document.getElementById('hero') !== null"):
                    ok = True
                    break
            except Exception:  # the page may still be navigating
                pass
            time.sleep(0.5)
    if ok and screenshot is not None:
        time.sleep(2)
        try:
            from PIL import ImageGrab

            ImageGrab.grab().save(screenshot)
        except Exception as exc:  # screenshots are best-effort diagnostics
            print(f"Screenshot failed: {exc}", file=sys.stderr)
    result.append(ok)
    window.destroy()


def launch(
    settings: GuiSettings | None = None,
    *,
    window: bool = True,
    open_browser: bool = True,
    host: str = "127.0.0.1",
    port: int | None = None,
    self_test: float | None = None,
    screenshot: str | Path | None = None,
) -> bool:
    """Start the DNAKit interface and block until it is closed.

    With ``window=True`` the interface opens in a native desktop window (pywebview). Without
    a desktop window backend or display it falls back to the default web browser. With
    ``window=False`` it serves the page at ``host:port`` and, if ``open_browser``, opens a
    browser tab. ``self_test`` (seconds) starts the interface, checks that the page renders,
    closes it and returns whether the check passed; ``screenshot`` saves the screen during that
    check (desktop window only).
    """

    _require_gradio()
    resolved = GuiSettings.local() if settings is None else settings
    from dnakit.gui._app import launch_options

    demo = build_app(resolved)
    demo.queue(default_concurrency_limit=2, max_size=32 if resolved.is_public else None)
    options = launch_options(resolved)
    webview: Any = None
    if window:
        if not _has_display():
            print("No display found; opening the interface in a web browser.", file=sys.stderr)
        else:
            try:
                import webview
            except ImportError:
                print(
                    'pywebview is not installed (pip install "dnakit[gui]"); '
                    "opening the interface in a web browser.",
                    file=sys.stderr,
                )
    if webview is None:
        if self_test is None:
            demo.launch(server_name=host, server_port=port, inbrowser=open_browser, **options)
            return True
        demo.launch(server_name=host, server_port=port, prevent_thread_lock=True, **options)
        try:
            return _page_ready(str(demo.local_url), self_test)
        finally:
            demo.close()

    demo.launch(
        server_name="127.0.0.1", server_port=port, prevent_thread_lock=True, quiet=True, **options
    )
    url = str(demo.local_url)
    try:
        webview.settings["ALLOW_DOWNLOADS"] = True
        win = webview.create_window(
            "DNAKit",
            url,
            width=WINDOW_SIZE[0],
            height=WINDOW_SIZE[1],
            min_size=MIN_WINDOW_SIZE,
            text_select=True,
        )
        result: list[bool] = []
        if self_test is None:
            webview.start()
        else:
            shot = None if screenshot is None else Path(screenshot)
            webview.start(_window_self_test, (win, self_test, shot, result))
        return bool(result[0]) if self_test is not None and result else self_test is None
    except Exception as exc:  # no usable window backend: keep the server and use a browser
        if self_test is not None:
            print(f"Desktop window failed: {exc}", file=sys.stderr)
            return False
        print(f"Desktop window failed ({exc}); opening {url} in a web browser.", file=sys.stderr)
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
        demo.block_thread()
        return True
    finally:
        demo.close()
