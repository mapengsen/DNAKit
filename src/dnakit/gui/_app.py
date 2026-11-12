"""Page layout of the DNAKit graphical interface, shared by the desktop app and the web demo."""

from __future__ import annotations

import base64
import json
import re
from importlib import resources
from pathlib import Path
from typing import Any

import gradio as gr

import dnakit
from dnakit.gui._backend import (
    DOCS_URL,
    EXAMPLE_SEQUENCE,
    FX_HEADERS,
    MCP_HEADERS,
    SK_HEADERS,
    GuiBackend,
    mark_selected,
    pick_first_column,
)
from dnakit.gui._chat import build_chat_tab
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


def logo_path() -> Path:
    """Return the bundled DNAKit logo (also the window and browser-tab icon)."""

    return Path(str(resources.files("dnakit.gui").joinpath("assets/dnakit-icon.png")))


def _logo_html() -> str:
    data = base64.b64encode(logo_path().read_bytes()).decode("ascii")
    return f'<img class="dk-logo" src="data:image/png;base64,{data}" alt="DNAKit logo">'


def catalog_counts() -> dict[str, Any] | None:
    """Deduplicated catalogue counts (operations, MCP tools, Skills) shipped with the interface."""

    path = resources.files("dnakit.gui").joinpath("assets/catalog_counts.json")
    if not path.is_file():
        return None
    counts: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return counts


def _hero(backend: GuiBackend) -> str:
    settings = backend.settings
    counts = catalog_counts()
    if counts:
        summary = (
            f"**{counts['operations_total']:,} DNA operations** ({counts['own_functions']:,} own "
            f"functions, {counts['mcp_bridge']:,} via the MCP bridge and "
            f"{counts['agent_skills']:,} Agent Skills), with **{counts['mcp_tools']:,} tools** "
            "callable through MCP"
        )
    else:
        summary = (
            f"**{backend.stats['agent_compatible']:,} functions** in {backend.stats['categories']} "
            f"modules and **{backend.skill_total} workflow Skills**"
        )
    head = (
        f"# {_logo_html()} DNAKit\n"
        "A unified DNA toolkit for Python, the command line and AI agents: "
        f"{summary}. Version {dnakit.__version__}. "
    )
    links = f" [Documentation]({DOCS_URL}) · [Source](https://github.com/mapengsen/DNAKit)"
    if settings.is_public:
        desktop = (
            ' Run it on your own computer with `pip install "dnakit[gui]"` and `dnakit gui`.'
            if _gui_released()
            else ""
        )
        return head + (
            "This public demo runs on a shared CPU; install DNAKit locally for large data, model "
            f"training and file-writing tools.{desktop}{links}"
        )
    return (
        head + f"Running on this computer; inputs and results are saved in `{settings.workdir}`."
        f"{links}"
    )


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

        with gr.Tab("Sequence (Example)"):
            seq_in = gr.Textbox(
                label="DNA sequence (plain text or one FASTA record)",
                value=EXAMPLE_SEQUENCE,
                lines=4,
            )
            with gr.Tabs():
                with gr.Tab("Reverse complement"):
                    b = gr.Button("Reverse-complement", variant="primary")
                    rc_head = gr.Markdown()
                    rc_out = gr.Textbox(
                        label="Reverse complement (5'→3')", lines=4, buttons=["copy"]
                    )
                    rc_file = gr.File(label="Download FASTA")
                    rc_code = gr.Markdown()
                    b.click(
                        backend.seq_reverse_complement,
                        [seq_in, state],
                        [state, rc_head, rc_out, rc_file, rc_code],
                    )
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

        with gr.Tab("Dataset (Example)"):
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

        with gr.Tab("Chat (AI agent)"):
            build_chat_tab(backend, state, files_label)

        counts = catalog_counts()
        listed = backend.stats["agent_compatible"]
        with gr.Tab(f"Function explorer ({counts['own_functions'] if counts else listed:,})"):
            note = (
                f" It lists all {listed:,} Agent-callable function names; names that implement the "
                f"same operation (aliases, batch versions) count once in the "
                f"{counts['own_functions']:,} functions."
                if counts
                else ""
            )
            gr.Markdown(
                "Search every Agent-callable DNAKit function (the same catalog an AI agent sees "
                "through the MCP server), then inspect its schema and run it with JSON arguments."
                + note
                + " Click a row to show the function below; *docs* opens its documentation page."
            )
            with gr.Row():
                fx_query = gr.Textbox(label="What do you want to do?", value="cpg island", scale=3)
                fx_cat = gr.Dropdown(backend.category_choices, value="", label="Module", scale=2)
                fx_btn = gr.Button("Search", variant="primary", scale=1)
            fx_head = gr.Markdown()
            fx_table = gr.Dataframe(
                headers=FX_HEADERS,
                datatype=["str", "markdown", "str", "str", "str"],
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
            fx_name.change(mark_selected, [fx_table, fx_name], [fx_table])
            fx_file.upload(backend.fx_upload, [fx_file, state], [state, fx_file_note])
            run_event = fx_run_btn.click(
                backend.fx_run,
                [fx_name, fx_args, fx_consent, state],
                [state, fx_status, fx_result, fx_plot, fx_files, fx_code],
            )
            fx_stop_btn.click(backend.fx_stop, None, [fx_status], cancels=[run_event])

        mcp_total = counts["mcp_tools"] if counts else None
        with gr.Tab(f"MCP tools ({mcp_total:,})" if mcp_total else "MCP tools"):
            if counts:
                own = counts["mcp_tools"] - counts["mcp_bridge"]
                gr.Markdown(
                    f"An AI agent connected to the DNAKit MCP server can call **{mcp_total:,} "
                    f"tools**: {own:,} DNAKit functions (searchable in the *Function explorer*) "
                    f"and {counts['mcp_bridge']:,} operations of third-party MCP servers that "
                    "DNAKit reaches through its MCP bridge. Search the indexed third-party tools "
                    "below and click a row for its arguments, source and how to call it; *docs* "
                    "opens its entry in the [MCP tool catalogue ↗]"
                    f"({DOCS_URL}llm_tools/mcp_tools/). They run once their server is installed "
                    f"locally (see the [MCP guide ↗]({DOCS_URL}llm_tools/mcp/)), not in this demo."
                )
            with gr.Row():
                mcp_query = gr.Textbox(label="Search third-party MCP tools", value="gwas", scale=4)
                mcp_btn = gr.Button("Search", variant="primary", scale=1)
            mcp_head = gr.Markdown()
            mcp_table = gr.Dataframe(
                headers=MCP_HEADERS,
                datatype=["str", "markdown", "str", "str", "str"],
                max_height=300,
                wrap=True,
                interactive=False,
            )
            mcp_detail = gr.Markdown(max_height=640, container=True, padding=True)
            with gr.Accordion("Indexed MCP servers", open=False):
                gr.Dataframe(
                    value=backend.mcp_servers(),
                    headers=["server", "DNA/gene tools"],
                    max_height=300,
                    interactive=False,
                )
            mcp_btn.click(backend.mcp_search, [mcp_query], [mcp_head, mcp_table, mcp_detail])
            mcp_query.submit(backend.mcp_search, [mcp_query], [mcp_head, mcp_table, mcp_detail])
            mcp_table.select(backend.mcp_select, [mcp_table], [mcp_detail, mcp_table])

        with gr.Tab(f"Skills ({backend.skill_total})"):
            gr.Markdown(
                "Curated DNA workflow Skills that guide an AI agent through multi-step analyses. "
                "Click a row to read the Skill below; *docs* opens its entry in the "
                f"[Skills catalogue ↗]({DOCS_URL}llm_tools/skills_catalog/), and the "
                f"[Skills guide ↗]({DOCS_URL}llm_tools/skills/) explains installing them."
            )
            with gr.Row():
                sk_query = gr.Textbox(label="Search Skills", value="variant annotation", scale=4)
                sk_btn = gr.Button("Search", variant="primary", scale=1)
            sk_head = gr.Markdown()
            sk_table = gr.Dataframe(
                headers=SK_HEADERS,
                datatype=["str", "markdown", "str", "str", "str"],
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
            sk_name.change(mark_selected, [sk_table, sk_name], [sk_table])

        demo.load(backend.fx_search, [fx_query, fx_cat], [fx_head, fx_table, fx_name])
        demo.load(backend.sk_search, [sk_query], [sk_head, sk_table, sk_name])
        demo.load(backend.mcp_search, [mcp_query], [mcp_head, mcp_table, mcp_detail])

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
