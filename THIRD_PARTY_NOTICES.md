# DNAKit Third-Party Notices

Last reviewed: 2026-08-23.

This file is an engineering compliance inventory, not legal advice. DNAKit itself is licensed under MIT. A dependency, external executable, database, paper, or user-supplied parameter table keeps its own terms; citing a paper does not by itself grant permission to redistribute its software or data.

## Distribution boundary

The DNAKit wheel and sdist do **not** bundle DiProDB numerical tables, Primer3, NUPACK, DSSR/3DNA, Dashing, BLAST, MMseqs2, sourmash, REBASE, JASPAR, DNA foundation-model checkpoints, NCBI, Ensembl, ENA, ENCODE, UCSC, or their databases. External programs are never downloaded automatically. Public-data adapters access providers only when explicitly called. Where an executable adapter exists, the user must install the tool separately and provide an explicit executable or installation path. A model checkpoint is downloaded to the local `ckpt/` directory only when the user explicitly invokes the representation API; it is not redistributed in DNAKit artifacts.

Python dependencies declared in `pyproject.toml` are resolved as separate distributions by the installer; their source code is not vendored into the DNAKit wheel. This inventory covers direct dependencies. Users and redistributors must also inspect resolved transitive dependencies for their actual environment and versions.

## Direct Python dependencies

| Scope | Component | License identifier or upstream label | DNAKit use |
| --- | --- | --- | --- |
| Runtime | PyYAML | MIT | YAML configuration |
| Runtime | Rich | MIT | terminal rendering and progress |
| Runtime | Typer | MIT | command-line interface |
| Runtime, Python <3.11 | tomli | MIT | TOML parsing |
| `viz` extra | CairoSVG | LGPL-3.0-or-later | optional SVG conversion |
| `viz` extra | Pillow | MIT-CMU | optional raster export |
| `io` extra | PyArrow | Apache-2.0 | optional Parquet I/O |
| `validation` extra | Biopython | Biopython License Agreement | validation only |
| `neural` extra | Hugging Face Hub, NumPy, scikit-learn, PyTorch, Transformers | Apache-2.0, BSD-family, and Apache-2.0 labels; verify resolved versions | checkpoint retrieval, inference, and k-means |
| `neural-caduceus` extra | mamba-ssm | Apache-2.0 label; verify resolved version | optional Caduceus runtime |
| `neural-enformer` extra | enformer-pytorch | MIT label; verify resolved version | optional Enformer runtime |
| `neural-evo2` extra | evo2 and transitive dependencies | verify the installed version and official repository | optional Evo 2 runtime |
| `docs` extra | MkDocs | BSD-2-Clause | documentation build |
| `docs` extra | Material for MkDocs | MIT | documentation theme |
| `docs` extra | mkdocstrings | ISC | API documentation |
| `docs` extra | nbmake | Apache-2.0 | notebook test gate |
| `dev` extra | build | MIT | local package build |
| `dev` extra | Hypothesis | MPL-2.0 | tests |
| `dev` extra | mypy | MIT | type checking |
| `dev` extra | pytest | MIT | tests |
| `dev` extra | pytest-cov | MIT | coverage |
| `dev` extra | Ruff | MIT | lint and formatting |
| `dev` extra | twine | Apache-2.0 | local distribution checks |
| `dev` extra | types-PyYAML | Apache-2.0 | type stubs |

Exact versions and license texts must be checked in the distributions actually installed. Useful upstream locations include [PyPI](https://pypi.org/), [Biopython licensing](https://biopython.org/DIST/docs/LICENSE), and each project repository.

## User-installed tools and supplied data

| Resource | DNAKit integration and distribution | Terms/cost warning |
| --- | --- | --- |
| Primer3 CLI (`primer3_core`, `oligotm`, `ntthal`) | Pure external CLI adapter; explicit paths; no Python binding and no bundled binary | Primer3 source files state GPL-2.0-or-later; verify the installed release and comply with its license. See the [official repository](https://github.com/primer3-org/primer3) and [manual](https://github.com/primer3-org/primer3/blob/main/src/primer3_manual.htm). |
| NUPACK | User installs separately; DNAKit does not download, bundle, or expose it as a hosted service | Separate NUPACK terms/subscription apply and access may require payment. Academic and commercial terms differ. Check the [official download and licensing page](https://www.nupack.org/download/overview) before use. |
| DSSR/3DNA | User installs separately; DNAKit only parses explicitly supplied output | Academic and commercial permissions differ; free academic availability is still subject to the issued license. Check the [official DSSR licensing notice](https://home.x3dna.org/highlights/x3dna-dssr-is-funded-and-dssr-basic-academic-is-free). |
| Dashing | User supplies an explicit executable; no binary/source is included in distributions | GPL-3.0; verify the installed version in the [official repository](https://github.com/dnbaker/dashing). |
| DiProDB-derived or other dinucleotide tables | DNAKit supplies only a JSON schema and loader. It contains no DiProDB numerical values. | No redistribution permission for the public table has been assumed. The user must obtain or create a table they are entitled to use and cite its actual source. See [DiProDB](https://diprodb.fli-leibniz.de/). |
| BLAST, MMseqs2, sourmash and external databases | Not installed or bundled; registered capabilities are conditional and version-specific | Verify each program and database license separately. Database access does not imply redistribution rights. |
| NCBI, Ensembl, ENA, ENCODE and UCSC public services/data | Queried or downloaded only on explicit user calls; source URLs and checksums are recorded; no remote dataset is embedded in distributions | Check provider usage, attribution, privacy, controlled-access and redistribution terms. dbGaP and other controlled data still require formal authorization. |
| Hugging Face DNA-model checkpoints | Downloaded into `ckpt/` after explicit model selection with a source manifest; weights are not bundled in wheel/sdist | DNABERT-2, NTv2, HyenaDNA, Caduceus, GROVER, LucaOne, GENERator, Enformer, and Evo 2 may have different terms. Review the official links on the [clustering feature page](docs/api/features/10_clustering.md#data-027-neural-clustering). |
| AlphaGenome | Loaded only after the user accepts access terms, obtains the checkpoint, and installs the official research code | The checkpoint is gated and carries non-commercial model terms. Review the [official repository](https://github.com/google-deepmind/alphagenome_research) and current model page. |
| JanusDNA | Official Harvard Dataverse file is downloaded with MD5 verification; the user must also supply the official source environment | Review the [official repository](https://github.com/Qihao-Duan/JanusDNA), [Dataverse DOI](https://doi.org/10.7910/DVN/HDT0RN), and applicable code/data terms. |

The DNABERT-2, NTv2, HyenaDNA, Caduceus, LucaOne, and GENERator Transformers
checkpoints include Python code required by their loaders. DNAKit rejects that execution by
default. It is enabled only after the user reviews the source and explicitly sets
`allow_remote_code=True`. This opt-in does not replace a code audit, dependency lock, or
license review.

Keeping a tool in a separate process or asking the user to install it reduces bundling and automatic-distribution risk, but it is not a legal conclusion that every intended use is permitted. Commercial deployment, hosted services, redistribution, or journal supplementary archives require a new review of the exact artifacts being distributed.

## Scientific attribution

Methods, papers, databases, and websites used as scientific foundations are listed in `docs/acknowledgements.md`. Those citations document provenance; they do not replace software/data licenses. Runtime provenance should record the actual tool/table name, version, source, parameters, and checksum.

## No assurance of zero cost or zero risk

DNAKit does not promise that any optional tool, database, or future version is free of charge. Terms, availability, and fees can change. Before each release, archive, hosted deployment, or commercial use, repeat the dependency and artifact audit and retain copies of the applicable license texts.
