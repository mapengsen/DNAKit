# Graphical interface

DNAKit's graphical interface runs the same functions as the Python API, the `dnakit` command line
and the MCP server, without writing code. It comes in two forms built from one code base:

- **Desktop app**: `dnakit gui` opens a native window on your own computer. Files stay local,
  there is no time limit, and functions that write files, download model weights or run external
  programs run after you confirm them.
- **Web app**: the same page hosted on a server for trying DNAKit in a browser. It runs on a
  shared CPU with time limits and does not execute functions with side effects.

## Install and start

```bash
pip install "dnakit[gui]"
dnakit gui
```

The window uses the operating system's web view: WebView2 on Windows 10/11 and WebKit on macOS.
On Linux, pywebview uses GTK WebKit if the system provides it; otherwise install the Qt backend
with `pip install "dnakit[gui-qt]"`. Its wheels need a recent Linux (glibc 2.28 or later) and the
`libxcb-cursor0` system library. Without a display or window backend, `dnakit gui` opens the same
interface in your default web browser.

| Option | Effect |
|---|---|
| `--workdir DIR` | Folder for uploaded inputs and results (default `./dnakit-gui`) |
| `--browser` | Open in the default web browser instead of a window |
| `--no-open` | Only serve the page, e.g. on a remote server |
| `--host`, `--port` | Address and port to serve on (default `127.0.0.1`, first free port from 7860) |
| `--public-demo` | Apply the web-app restrictions described below |

On a remote server without a desktop, run `dnakit gui --no-open --port 7860` and open the page
through SSH port forwarding (`ssh -L 7860:127.0.0.1:7860 server`) or your editor's port
forwarding.

## What the interface contains

| Tab | Contents |
|---|---|
| Sequence | 240 descriptors, six-frame ORFs, IUPAC/regex motif search, sequence comparison, validation |
| Dataset | Deduplication (exact or reverse complement), train/valid/test split, format conversion, a descriptor feature table for every record |
| Function explorer | Search every Agent-callable function, read its schema, run it with JSON arguments |
| Skills | Search and read the bundled DNA workflow Skills |
| Use with agents | Connect the same functions to Claude Code, Codex or other MCP clients |

Every result shows the equivalent command line, Python call and, in the function explorer, the
MCP request an agent sends, so a task tried in the interface can be scripted directly.

## Desktop app and web app

| | Desktop app (`GuiSettings.local`) | Web app (`GuiSettings.public_demo`) |
|---|---|---|
| Files | Saved in the work folder | Per-visitor temporary folder, deleted after 6 hours |
| Time limit per call | None (a **Stop** button ends a call) | 60 s |
| Upload size | Unlimited | 20 MB |
| Functions that write files, download models or run programs | Run after confirmation | Not executed |
| Path arguments | Any path | Only inside the visitor's folder |

Single-sequence tools accept up to 30,000 nt because the sequence is passed on the command line
(Windows limits a command line to 32,767 characters); use the Dataset tab for longer sequences.

## Python

```python
from dnakit.gui import GuiSettings, build_app, launch

launch(GuiSettings.local("results"))  # desktop window, results in ./results
app = build_app(GuiSettings.public_demo())  # Gradio Blocks app for your own server
```

`packaging/README.md` describes the systemd and nginx configuration used for the hosted web app.
