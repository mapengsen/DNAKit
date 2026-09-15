---
name: bioinformatics-methods-results-writer
description: >
  Write Methods and Results sections for bioinformatics manuscripts from code,
  notebooks, logs, figures, and tables. Use when drafting or revising
  computational biology prose, especially to preserve reproducibility details,
  figure/table references, and concise results reporting. Triggers on: "methods
  section", "results section", "manuscript", "write methods", "write results",
  "draft methods", "draft results", "publication", "paper", "scientific prose",
  "reproducibility", "methods and results".
---

## DNAKit execution

Use `dnakit skills doctor bioinformatics-methods-results-writer` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Bioinformatics Methods & Results Writer

Draft manuscript **Methods** and **Results** sections from code, notebooks, console output, and generated figures/tables. Optimize for reproducibility, clarity, and section-appropriate tone.

## Core objective

Turn computational artifacts into publication-ready prose that a field expert could use to understand, reproduce, and verify the analysis.

## Input assumptions

Expect some combination of:

* code cells, scripts, or notebooks
* package lists, environment files, or session info
* parameter settings and random seeds
* tables, plots, and summary outputs
* notes on data provenance, preprocessing, or filtering
* study constraints such as IRB, consent, privacy, or licensing

## Output expectations

Produce one or more of the following, as requested:

* a Methods paragraph or subsection
* a Results paragraph or subsection
* parallel Methods/Results drafts with aligned structure
* a short checklist of missing details
* a revision of existing prose to better match the evidence

## Operating rules

### 1) Separate Methods from Results

* **Methods**: explain how the work was done.
* **Results**: report what was found.
* Do not mix interpretation into Results.
* Do not turn Methods into a recipe; keep it readable and structured.

### 2) Methods should be reproducible

Always capture, when available:

* data source and version or accession
* inclusion/exclusion or filtering criteria
* preprocessing steps and transformations
* software, package, and model names
* version numbers
* key parameters and non-default settings
* random seeds for stochastic steps
* hardware or runtime environment when relevant
* code/data availability, repository location, and license when known
* ethics approval, consent, privacy handling, and legal constraints when relevant

### 3) Methods style

* Write in **past tense**.
* Use active or passive voice as appropriate; prefer the clearest option.
* Order content **logically and chronologically**.
* Use subheadings when the analysis has distinct stages.
* Mirror the structure of the Results section when possible.
* Write for a technically trained reader, not a lay audience.

### 4) Results style

* Write in **past tense**.
* Be concise, objective, and information-dense.
* Report the key findings, not every intermediate number.
* Include statistics, effect sizes, confidence intervals, and p-values when relevant.
* Mention **negative or null results** when they matter.
* Reference every figure or table that appears in the Results, and do not leave visuals unexplained.
* Do not repeat in text what is already obvious from a figure or table.

### 5) Notebook-to-narrative mapping

Translate artifacts into prose using these defaults:

* `read_csv` / `load` / `import` → dataset acquisition or import sentence
* `set.seed(...)` / RNG initialization → reproducibility sentence
* model constructor / fit call → model specification and training sentence
* transformation code → preprocessing sentence
* printed metric or summary output → Results sentence with the metric and context
* plot generation → Results sentence that points the reader to the figure
* dataframe display → Results sentence that summarizes the table

### 6) Preserve provenance

When multiple steps are involved, describe:

* what was taken in as input
* what was changed
* what was computed
* what was retained for downstream analysis

If a detail is missing but likely needed for reproducibility, flag it explicitly rather than inventing it.

### 7) Prefer concrete wording

Use specific nouns and values over generic prose.

* Say "pandas v2.1" rather than "a data library."
* Say "100 trees, max depth 5" rather than "default settings."
* Say "Figure 2" rather than "the plot above."

## Drafting workflow

1. Extract the analysis stages from the input.
2. Classify each stage as Methods or Results content.
3. Build a paragraph outline that follows the analysis order.
4. Convert each step into plain scientific prose.
5. Insert figure/table references where outputs are presented.
6. Check for missing reproducibility details.
7. Remove interpretation from Results and remove outcome claims from Methods.
8. Return the draft plus a brief list of any unresolved gaps.

## Suggested section templates

### Methods template

Start with the data and study design, then move through preprocessing, analysis, and statistical or computational methods.

> Data were obtained from [source]. We [filtered/transformed/normalized] the input by [criterion]. Analyses were performed in [language] using [software/packages] version [X]. We fit [model/test] with [key parameters]. Randomness was controlled by setting the seed to [N]. [Ethics/privacy/availability statement, if applicable].

### Results template

Start with the primary outcome, then add supporting metrics and comparisons.

> [Primary finding]. Compared with [reference], [group/model] showed [direction and magnitude] (metric = [value], p = [value]). Figure [N] summarizes [what the figure shows]. Table [N] provides [what the table contains].

## Quality checks before finalizing

* Methods and Results should use different verbs and different purposes.
* Every reported computational result should have a clear source in the notebook or code.
* Every figure/table mentioned in Results should be described in text.
* All non-default computational choices should be stated.
* No speculative interpretation should leak into Results.
* No missing provenance should be silently fabricated.

## Default refusal boundary

Do not invent methods details that are not supported by the supplied code, notebook, or notes. When details are absent, mark them as unknown and request them or provide a bracketed placeholder.
