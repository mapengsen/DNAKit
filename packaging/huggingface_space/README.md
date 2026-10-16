---
title: DNAKit
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: 6.29.1
python_version: "3.12"
app_file: app.py
pinned: false
license: mit
short_description: Unified DNA toolkit for Python, the CLI and AI agents
---

# DNAKit web demo

Browser access to [DNAKit](https://github.com/mapengsen/DNAKit), a DNA sequence toolkit that exposes
the same functions through Python, the `dnakit` command line and an MCP server for AI agents.

| Tab | What it does |
|---|---|
| Sequence | 240 descriptors, six-frame ORFs, IUPAC/regex motif search, sequence comparison, validation |
| Dataset | Deduplication (exact or reverse complement), train/valid/test split, format conversion, a 240-column descriptor feature table |
| Function explorer | Search every Agent-callable function, read its schema, run it with JSON arguments, and copy the equivalent Python and MCP call |
| Skills | Search and read the bundled DNA workflow Skills |
| Use with agents | Install DNAKit and connect it to Claude Code, Codex or other MCP clients |

The page is the `dnakit.gui` interface (the same code as the `dnakit gui` desktop app) in public-demo mode. Functions that write files outside the visitor's session folder,
download model weights or run external programs are listed but not executed; each call stops after
60 s. Install DNAKit locally (`pip install dnakit`) for large data, model training and those tools.
Uploaded files stay in a temporary per-visitor folder and are deleted when the Space restarts.
