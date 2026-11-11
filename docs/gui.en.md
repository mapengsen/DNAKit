# Graphical interface

DNAKit's graphical interface runs the same functions as the Python API, the `dnakit` command line
and the MCP server, without writing code. It comes in two forms built from one code base:

- **Desktop app**: `dnakit gui` opens a native window on your own computer. Files stay local,
  there is no time limit, and functions that write files, download model weights or run external
  programs run after you confirm them.
- **Web app**: the same page hosted on a server for trying DNAKit in a browser. It runs on a
  shared CPU with time limits and does not execute functions with side effects.

## Download the desktop app

The desktop app can be downloaded as a ready-to-run program that needs no Python:

| System | Download | Start |
|---|---|---|
| Windows 10 / 11 (64-bit) | [DNAKit-windows-x64.zip](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-windows-x64.zip) | Unpack the zip, open the `DNAKit` folder and double-click `DNAKit.exe` |
| macOS, Apple silicon | [DNAKit-macos-arm64.dmg](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-arm64.dmg) | Open the disk image, drag **DNAKit** to Applications and open it |
| macOS, Intel | [DNAKit-macos-x64.dmg](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-x64.dmg) | As above |
| Linux x86-64 (glibc 2.35 or later, e.g. Ubuntu 22.04) | [DNAKit-linux-x64.tar.gz](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-linux-x64.tar.gz) | `tar -xzf DNAKit-linux-x64.tar.gz && ./DNAKit/DNAKit` |

The apps are not code-signed, so the system asks for one confirmation at the first start:

- **Windows:** in the *Windows protected your PC* dialog click *More info*, then *Run anyway*.
- **macOS:** after the message that DNAKit cannot be verified, open *System Settings* →
  *Privacy & Security* and click *Open Anyway* (macOS 14 and earlier: right-click the app and
  choose *Open*). Alternatively run `xattr -dr com.apple.quarantine /Applications/DNAKit.app`.

The downloaded app keeps uploaded inputs and results in the `DNAKit` folder of your home
directory and, when started by double click, writes its messages to `DNAKit/desktop.log` there.
It contains the interface and the functions of a plain `pip install "dnakit[gui]"`. Functions
that need optional packages (deep-learning models, the `bio` extras) or external programs report
what is missing; use the pip installation below for those. Started from a terminal with
arguments, the program is the `dnakit` command line, for example `DNAKit describe ACGT`.

The downloads are rebuilt from the `main` branch by the
[Desktop app workflow](https://github.com/mapengsen/DNAKit/actions/workflows/desktop.yml); the
build recipe is in [`desktop-app/`](https://github.com/mapengsen/DNAKit/tree/main/desktop-app).

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
| Sequence (Example) | Reverse complement, 240 descriptors, six-frame ORFs, IUPAC/regex motif search, sequence comparison, validation |
| Dataset (Example) | Deduplication (exact or reverse complement), train/valid/test split, format conversion, a descriptor feature table for every record |
| Function explorer | Search every Agent-callable function, open its section of this documentation, read its schema, run it with JSON arguments |
| MCP tools | Search the DNA and gene tools of third-party MCP servers; a clicked tool shows its arguments, source repository and the `call_mcp_tool` call, and *docs* opens its entry in the MCP tool catalogue |
| Skills | Search and read the bundled DNA workflow Skills; *docs* opens a Skill's entry in the Skills catalogue |

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
