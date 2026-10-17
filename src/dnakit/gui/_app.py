"""Page layout of the DNAKit graphical interface, shared by the desktop app and the web demo."""

from __future__ import annotations

import base64
import re
from importlib import resources
from pathlib import Path
from typing import Any

import gradio as gr

import dnakit
from dnakit.gui._backend import EXAMPLE_SEQUENCE, GuiBackend, pick_first_column
from dnakit.gui._settings import GuiSettings

# System fonts only: available on every platform, and no web-font download (Google Fonts is
# unreachable in some regions, which leaves the page in a fallback font).
SANS = [
    "-apple-system",
    "BlinkMacSystemFont",
    "Segoe UI",
    "Roboto",
    "Helvetica Neue",
    "Arial",
    "PingFang SC",
    "Microsoft YaHei",
    "Noto Sans CJK SC",
    "sans-serif",
]
MONO = ["SFMono-Regular", "Menlo", "Consolas", "Liberation Mono", "DejaVu Sans Mono", "monospace"]

CSS = """
.gradio-container {max-width: 1180px !important; margin: auto;}
#hero h1 {margin-bottom: 0.2rem; display: flex; align-items: center; gap: 0.45em;}
#hero h1 img.dk-logo {height: 1.5em; width: auto; margin: 0;}
#hero p {margin-top: 0; color: var(--body-text-color-subdued);}
.dk-plot {display: flex; justify-content: center;}
.dk-plot svg {width: 100%; max-width: 480px; height: auto;}
"""


def _gui_released() -> bool:
    """Whether this DNAKit release ships ``dnakit gui`` (first in the release after 0.2.0)."""

    parts = tuple(int(p) for p in re.findall(r"\d+", dnakit.__version__)[:3])
    return parts > (0, 2, 0)


def _agent_guide(settings: GuiSettings) -> str:
    desktop = ""
    if settings.is_public and _gui_released():
        desktop = (
            "\n**Desktop app.** Run this interface on your own computer, without the demo limits:\n"
            '```bash\npip install "dnakit[gui]"\ndnakit gui\n```\n'
        )
    return f"""
### Use the same functions from an AI agent

Everything in this interface also runs from Python, the command line, or an MCP-capable agent
(Claude Code, Codex, Cursor and others). The agent finds functions with `search_dnakit_tools`,
reads the schema with `describe_dnakit_tool` and runs them with `call_dnakit_tool`, exactly like
the *Function explorer* tab.
{desktop}
**1. Install**
```bash
python3 -m venv ~/.venvs/dnakit-mcp
~/.venvs/dnakit-mcp/bin/python -m pip install "dnakit[agent]=={dnakit.__version__}"
```

**2. Connect the MCP server**
```bash
claude mcp add dnakit -- ~/.venvs/dnakit-mcp/bin/dnakit-mcp     # Claude Code
codex mcp add dnakit -- ~/.venvs/dnakit-mcp/bin/dnakit-mcp      # Codex
```
Other clients use the `mcpServers` JSON format:
```json
{{"mcpServers": {{"dnakit": {{"command": "/absolute/path/to/dnakit-mcp"}}}}}}
```

**3. Optional: install the planning and workflow Skills**
```bash
dnakit skills install dnakit --agent claude-code
dnakit skills install dnakit-plan --agent codex
```

Documentation: <https://mapengsen.github.io/DNAKit/en/> ·
Source: <https://github.com/mapengsen/DNAKit>
"""


def logo_path() -> Path:
    """Return the bundled DNAKit logo (also the window and browser-tab icon)."""

    return Path(str(resources.files("dnakit.gui").joinpath("assets/dnakit-icon.png")))


def _logo_html() -> str:
    data = base64.b64encode(logo_path().read_bytes()).decode("ascii")
    return f'<img class="dk-logo" src="data:image/png;base64,{data}" alt="DNAKit logo">'


def _hero(backend: GuiBackend) -> str:
    settings = backend.settings
    head = (
        f"# {_logo_html()} DNAKit\n"
        "A unified DNA toolkit for Python, the command line and AI agents: "
        f"**{backend.stats['agent_compatible']:,} functions** in {backend.stats['categories']} "
        f"modules and **{backend.skill_total} workflow Skills**. Version {dnakit.__version__}. "
    )
    if settings.is_public:
        return head + (
            "This public demo runs on a shared CPU; install DNAKit locally for large data, model "
            "training and file-writing tools."
        )
    return head + f"Running on this computer; inputs and results are saved in `{settings.workdir}`."


def build_app(settings: GuiSettings) -> gr.Blocks:
    """Build the interface for ``settings`` without starting a server."""

    backend = GuiBackend(settings)
    public = settings.is_public
    limit = f" (≤ {settings.max_upload_mb} MB)" if settings.max_upload_mb else ""
    records = settings.max_table_records
    table_note = f" (≤ {records:,} records)" if records else ""
    files_label = (
        "Files written in your session folder" if public else f"Files written in {settings.workdir}"
    )

    demo = gr.Blocks(title="DNAKit")
    with demo:
        state = gr.State(None)
        gr.Markdown(_hero(backend), elem_id="hero")

        with gr.Tab("Sequence"):
            seq_in = gr.Textbox(
                label="DNA sequence (plain text or one FASTA record)",
                value=EXAMPLE_SEQUENCE,
                lines=4,
            )
            with gr.Tabs():
                with gr.Tab("Descriptors (240)"):
                    b = gr.Button("Compute descriptors", variant="primary")
                    d_head = gr.Markdown()
                    d_table = gr.Dataframe(headers=["field", "value"], max_height=420, wrap=True)
                    d_file = gr.File(label="Download CSV")
                    d_code = gr.Markdown()
                    b.click(
                        backend.seq_descriptors,
                        [seq_in, state],
                        [state, d_head, d_table, d_file, d_code],
                    )
                with gr.Tab("ORFs"):
                    orf_len = gr.Number(label="Minimum ORF length (nt)", value=60, precision=0)
                    b = gr.Button("Scan six frames", variant="primary")
                    o_head = gr.Markdown()
                    o_table = gr.Dataframe(
                        headers=["strand", "frame", "start", "end", "length (nt)", "translation"],
                        wrap=True,
                    )
                    o_code = gr.Markdown()
                    b.click(
                        backend.seq_orfs, [seq_in, orf_len, state], [state, o_head, o_table, o_code]
                    )
                with gr.Tab("Motif search"):
                    with gr.Row():
                        motif = gr.Textbox(label="Motif", value="CAYCA")
                        motif_mode = gr.Dropdown(
                            ["iupac", "exact", "regex"], value="iupac", label="Mode"
                        )
                    b = gr.Button("Search both strands", variant="primary")
                    m_head = gr.Markdown()
                    m_table = gr.Dataframe(headers=["strand", "start", "end", "match"])
                    m_code = gr.Markdown()
                    b.click(
                        backend.seq_motif,
                        [seq_in, motif, motif_mode, state],
                        [state, m_head, m_table, m_code],
                    )
                with gr.Tab("Compare"):
                    other = gr.Textbox(
                        label="Second sequence", value=EXAMPLE_SEQUENCE[5:110], lines=2
                    )
                    with gr.Row():
                        cmp_method = gr.Dropdown(
                            [
                                ("Edit distance (Levenshtein)", "edit"),
                                ("Hamming", "hamming"),
                                ("k-mer Jaccard", "kmer_jaccard"),
                                ("Exact / reverse complement", "exact"),
                            ],
                            value="edit",
                            label="Method",
                        )
                        cmp_k = gr.Number(label="k (k-mer Jaccard)", value=3, precision=0)
                    b = gr.Button("Compare", variant="primary")
                    c_head = gr.Markdown()
                    c_json = gr.JSON(label="Result")
                    c_code = gr.Markdown()
                    b.click(
                        backend.seq_compare,
                        [seq_in, other, cmp_method, cmp_k, state],
                        [state, c_head, c_json, c_code],
                    )
                with gr.Tab("Validate"):
                    b = gr.Button("Normalise and validate", variant="primary")
                    v_head = gr.Markdown()
                    v_json = gr.JSON(label="Audit")
                    v_code = gr.Markdown()
                    b.click(backend.seq_validate, [seq_in, state], [state, v_head, v_json, v_code])

        with gr.Tab("Dataset"):
            gr.Markdown(
                f"Upload FASTA, FASTQ, GenBank, CSV, TSV or JSON(L){limit}, or use the example "
                "(26 sequences including one exact and one reverse-complement duplicate)."
            )
            with gr.Row():
                ds_file = gr.File(label="Sequence file", file_count="single", type="filepath")
                ds_example = gr.Checkbox(label="Use the example file", value=True)
            with gr.Tabs():
                with gr.Tab("Deduplicate"):
                    eq = gr.Dropdown(
                        [("Exact", "exact"), ("Exact or reverse complement", "reverse_complement")],
                        value="reverse_complement",
                        label="Duplicates are",
                    )
                    b = gr.Button("Deduplicate", variant="primary")
                    dd_head = gr.Markdown()
                    dd_table = gr.Dataframe(
                        headers=["kept", "group size", "members", "orientations"], wrap=True
                    )
                    dd_files = gr.File(label="Download", file_count="multiple")
                    dd_code = gr.Markdown()
                    b.click(
                        backend.ds_deduplicate,
                        [ds_file, ds_example, eq, state],
                        [state, dd_head, dd_table, dd_files, dd_code],
                    )
                with gr.Tab("Train/valid/test split"):
                    with gr.Row():
                        sp_method = gr.Dropdown(
                            [("Seeded random", "random"), ("Stable hash of record ID", "hash")],
                            value="random",
                            label="Method",
                        )
                        sp_ratios = gr.Textbox(label="Ratios", value="train=0.8,valid=0.1,test=0.1")
                        sp_seed = gr.Number(label="Seed", value=42, precision=0)
                    b = gr.Button("Split", variant="primary")
                    sp_head = gr.Markdown()
                    sp_table = gr.Dataframe(headers=["record", "subset"], max_height=320)
                    sp_files = gr.File(label="Download", file_count="multiple")
                    sp_code = gr.Markdown()
                    b.click(
                        backend.ds_split,
                        [ds_file, ds_example, sp_method, sp_ratios, sp_seed, state],
                        [state, sp_head, sp_table, sp_files, sp_code],
                    )
                with gr.Tab("Convert format"):
                    fmt = gr.Dropdown(
                        ["fasta", "genbank", "csv", "tsv", "json", "jsonl"],
                        value="csv",
                        label="Output format",
                    )
                    b = gr.Button("Convert", variant="primary")
                    cv_head = gr.Markdown()
                    cv_files = gr.File(label="Download", file_count="multiple")
                    cv_code = gr.Markdown()
                    b.click(
                        backend.ds_convert,
                        [ds_file, ds_example, fmt, state],
                        [state, cv_head, cv_files, cv_code],
                    )
                with gr.Tab("Descriptor feature table"):
                    gr.Markdown(
                        "All 240 descriptor fields for every record, ready for machine learning"
                        f"{table_note}."
                    )
                    b = gr.Button("Build feature table", variant="primary")
                    ft_head = gr.Markdown()
                    ft_table = gr.Dataframe(max_height=360)
                    ft_files = gr.File(label="Download", file_count="multiple")
                    ft_code = gr.Markdown()
                    b.click(
                        backend.ds_descriptor_table,
                        [ds_file, ds_example, state],
                        [state, ft_head, ft_table, ft_files, ft_code],
                    )

        with gr.Tab(f"Function explorer ({backend.stats['agent_compatible']:,})"):
            gr.Markdown(
                "Search every Agent-callable DNAKit function (the same catalog an AI agent sees "
                "through the MCP server), then inspect its schema and run it with JSON arguments."
            )
            with gr.Row():
                fx_query = gr.Textbox(label="What do you want to do?", value="cpg island", scale=3)
                fx_cat = gr.Dropdown(backend.category_choices, value="", label="Module", scale=2)
                fx_btn = gr.Button("Search", variant="primary", scale=1)
            fx_head = gr.Markdown()
            fx_table = gr.Dataframe(
                headers=["function", "module", "effect", "description"],
                max_height=300,
                wrap=True,
                interactive=False,
            )
            fx_name = gr.Dropdown(label="Function", choices=[], allow_custom_value=True)
            fx_doc = gr.Markdown()
            with gr.Accordion("Input JSON schema", open=False):
                fx_schema = gr.JSON()
            with gr.Row():
                with gr.Column(scale=3):
                    fx_args = gr.Code(label="Arguments (JSON)", language="json", lines=8)
                with gr.Column(scale=2):
                    fx_file = gr.File(label="Optional input file", type="filepath")
                    fx_file_note = gr.Markdown()
            fx_consent = gr.Checkbox(
                label="Allow this function to write files, download model weights or run "
                "external programs on this computer",
                value=False,
                visible=False,
            )
            with gr.Row():
                fx_run_btn = gr.Button("Run", variant="primary", interactive=False, scale=4)
                fx_stop_btn = gr.Button("Stop", variant="stop", scale=1)
            fx_status = gr.Markdown()
            fx_plot = gr.HTML()
            fx_result = gr.JSON(label="Result")
            fx_files = gr.File(label=files_label, file_count="multiple", visible=False)
            fx_code = gr.Markdown()
            fx_btn.click(backend.fx_search, [fx_query, fx_cat], [fx_head, fx_table, fx_name])
            fx_query.submit(backend.fx_search, [fx_query, fx_cat], [fx_head, fx_table, fx_name])
            fx_table.select(pick_first_column, [fx_table], [fx_name])
            fx_name.change(
                backend.fx_describe,
                [fx_name],
                [fx_doc, fx_args, fx_schema, fx_run_btn, fx_consent],
            )
            fx_file.upload(backend.fx_upload, [fx_file, state], [state, fx_file_note])
            run_event = fx_run_btn.click(
                backend.fx_run,
                [fx_name, fx_args, fx_consent, state],
                [state, fx_status, fx_result, fx_plot, fx_files, fx_code],
            )
            fx_stop_btn.click(backend.fx_stop, None, [fx_status], cancels=[run_event])

        with gr.Tab(f"Skills ({backend.skill_total})"):
            gr.Markdown(
                "Curated DNA workflow Skills that guide an AI agent through multi-step analyses."
            )
            with gr.Row():
                sk_query = gr.Textbox(label="Search Skills", value="variant annotation", scale=4)
                sk_btn = gr.Button("Search", variant="primary", scale=1)
            sk_head = gr.Markdown()
            sk_table = gr.Dataframe(
                headers=["skill", "category", "source", "description"],
                max_height=300,
                wrap=True,
                interactive=False,
            )
            sk_name = gr.Dropdown(label="Skill", choices=[], allow_custom_value=True)
            sk_doc = gr.Markdown(max_height=640, container=True, padding=True)
            sk_btn.click(backend.sk_search, [sk_query], [sk_head, sk_table, sk_name])
            sk_query.submit(backend.sk_search, [sk_query], [sk_head, sk_table, sk_name])
            sk_table.select(pick_first_column, [sk_table], [sk_name])
            sk_name.change(backend.sk_read, [sk_name], [sk_doc])

        with gr.Tab("Use with agents"):
            gr.Markdown(_agent_guide(settings))

        demo.load(backend.fx_search, [fx_query, fx_cat], [fx_head, fx_table, fx_name])
        demo.load(backend.sk_search, [sk_query], [sk_head, sk_table, sk_name])

    return demo


def launch_options(settings: GuiSettings) -> dict[str, Any]:
    """Return the ``Blocks.launch`` keyword arguments that belong to the page itself."""

    options: dict[str, Any] = {
        "theme": gr.themes.Soft(primary_hue="blue", font=SANS, font_mono=MONO),
        "css": CSS,
        "allowed_paths": [str(settings.workdir)],
        "footer_links": ["api", "gradio"],
        "favicon_path": str(logo_path()),
    }
    if settings.max_upload_mb:
        options["max_file_size"] = f"{settings.max_upload_mb}mb"
    return options
