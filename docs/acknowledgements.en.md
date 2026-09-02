# Acknowledgments

This page summarizes DNAKit's main sources, licenses, third-party use boundaries and disclaimers for centralized review and traceability.

## Acknowledgments and Primary Sources {#acknowledgements}

DNAKit would like to thank the paper authors, database maintainers, standards organizations and open source project contributors who provided the theoretical foundation, public parameters, format specifications, software interfaces and structural data for this project.

This section provides a unified summary of the methods, papers, databases and websites used in the project. The function page will still retain the PMID or usage boundary directly related to the specific parameters to facilitate item-by-item tracing. Listing of a source does not imply that DNAKit has obtained permission to redistribute its data or software.

### Method, paper and format basis {#methods-and-references}

| Function | DNAKit Type | Primary Basis Should Be Cited |
| --- | --- | --- |
| FASTA/FASTQ, GenBank, GFF3, BED, AGP | `reimplementation` / Format integration | [NCBI GenBank sample records ](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html), [Sequence Ontology GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), [UCSC BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), [NCBI AGP 2.1](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/); GenBank only supports documented subsets |
| Coordinate and sequence object design | `native`, refer to the existing interface | Biopython `Seq`/`SeqRecord` and scikit-bio `DNA` are only for design comparison |
| Shannon entropy | `reimplementation` | C. E. Shannon, *A Mathematical Theory of Communication* (1948), DOI [`10.1002/j.1538-7305.1948.tb01338.x`](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x) |
| linguistic complexity | `reimplementation` | Published definition of observed/possible k-word product; result saves exact formula |
| LZ76 complexity | `reimplementation` | Lempel & Ziv, 1976, DOI [`10.1109/TIT.1976.1055501`](https://doi.org/10.1109/TIT.1976.1055501) |
| 15 sets of dinucleotide attributes | User-provided parameter table | DNAKit only defines 15×16 JSON schema and does not have built-in DiProDB values. If users obtain tables from DiProDB or other sources by themselves, they should record the [DiProDB paper](https://doi.org/10.1093/nar/gkn597), actual table pages, original paper, version and SHA-256; for field mapping, see [240-item descriptor table](api/features/05_all_descriptors.md) |
| Standard genetic code table | `reimplementation` / Standard data | [NCBI Standard Genetic Code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) |
| Needleman–Wunsch global alignment | `reimplementation` | Needleman & Wunsch, 1970, DOI [`10.1016/0022-2836(70)90057-4`](https://doi.org/10.1016/0022-2836(70)90057-4) |
| Smith–Waterman local alignment | `reimplementation` | Smith & Waterman, 1981, DOI [`10.1016/0022-2836(81)90087-5`](https://doi.org/10.1016/0022-2836(81)90087-5) |
| Levenshtein distance | `reimplementation` | V. I. Levenshtein, 1966 |
| MinHash | `reimplementation` | A. Z. Broder, 1997, resemblance/containment sketching |
| EvoAug sequence generation | `reimplementation` | Lee et al., 2023, [EvoAug paper ](https://doi.org/10.1186/s13059-023-02941-w); [Official PyTorch code ](https://github.com/p-koo/evoaug); DNAKit adopts an independent `DNASequence` level implementation and does not introduce PyTorch or Gaussian noise |
| k-mer remains shuffled | `reimplementation` | [uShuffle paper ](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/); DNAKit is implemented using random Euler paths on de Bruijn multigraphs to accurately maintain specified overlapping k-mer counts |
| Dashing Jaccard adapter | `adapter` | Baker & Langmead, 2019, DOI [`10.1186/s13059-019-1875-0`](https://doi.org/10.1186/s13059-019-1875-0); Actual Dashing version, k, exact/HLL mode and sketch size |
| DNA basic model rep and k-means | `adapter` | Reference the selected model and checkpoint during actual runtime: [DNABERT-2](https://github.com/MAGICS-LAB/DNABERT_2), [Nucleotide Transformer](https://github.com/instadeepai/nucleotide-transformer), [HyenaDNA](https://github.com/HazyResearch/hyena-dna), [Caduceus](https://github.com/kuleshov-group/caduceus), [GROVER](https://huggingface.co/PoetschLab/GROVER), [LucaOne](https://github.com/LucaOne/LucaOne), [GENERator](https://github.com/GenerTeam/GENERator), [Enformer PyTorch](https://github.com/lucidrains/enformer-pytorch), [AlphaGenome research](https://github.com/google-deepmind/alphagenome_research), [JanusDNA](https://github.com/Qihao-Duan/JanusDNA) or [Evo 2](https://github.com/ArcInstitute/evo2); clustering record actual scikit-learn Version, pooling, PCA, seed and checkpoint sources |
| NT Revised/Genomic Benchmarks task classification | `adapter` | Cite the [Enformer](https://doi.org/10.1038/s41592-021-01252-x) backbone, the [NT Revised](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised) or [Genomic Benchmarks](https://doi.org/10.1186/s12863-023-01123-8) task, and the [GENERator Appendix C.4](https://arxiv.org/abs/2502.07272) full-fine-tuning protocol; also record the actual checkpoint file, label mapping, and source link |
| Fréchet DNA distance | `reimplementation` + `adapter` | Fréchet Gaussian distance follows the mathematical form of [FCD original paper ](https://doi.org/10.1021/acs.jcim.8b00234); the representation separately refers to the actual model and checkpoint, and the default model refers to [LucaOne](https://doi.org/10.1038/s42256-025-01044-4). This indicator is not ChemNet FCD |
| DNA Frag / SNN | `reimplementation` + `adapter` | Index formula reference [MOSES](https://doi.org/10.3389/fphar.2020.565644); DNA adaptation uses fixed-length k-mer to replace BRICS fragments, and hashed k-mer bit fingerprint to replace Morgan fingerprint. Values cannot be directly compared with molecular indicators |
| DNA 260 nm extinction coefficient | `reimplementation` | Warshaw & Tinoco, 1966, DOI [`10.1016/0022-2836(66)90115-X`](https://doi.org/10.1016/0022-2836(66)90115-X); Cantor, Warshaw & Shapiro, 1970, DOI [`10.1002/bip.1970.360090909`](https://doi.org/10.1002/bip.1970.360090909); actual parameter set version |
| DNA nearest-neighbor/Tm | `reimplementation` | SantaLucia, 1998, DOI [`10.1073/pnas.95.4.1460`](https://doi.org/10.1073/pnas.95.4.1460) |
| Primer3 thermodynamic/design adapter | `adapter` | Untergasser et al., 2012, DOI [`10.1093/nar/gks596`](https://doi.org/10.1093/nar/gks596); Actual Primer3 CLI version, executable path and parameter directory |
| NUPACK secondary structure adapter | `adapter` | Zadeh et al., 2011, DOI [`10.1002/jcc.21596`](https://doi.org/10.1002/jcc.21596); Actual NUPACK version, model, temperature and salt conditions |
| 3DNA/DSSR parameter analysis | `adapter` | 3DNA DOI [`10.1093/nar/gkg680`](https://doi.org/10.1093/nar/gkg680); DSSR DOI [`10.1093/nar/gkv716`](https://doi.org/10.1093/nar/gkv716); actual output version |
| PDB three-dimensional structure sample | Data record | Actual RCSB PDB ID, experimental method, model number and coordinate file SHA-256; current samples are 1BNA, 1AC7, 139D |
| Restriction enzymes | Versioned built-in mini-directory/optional database | Actual enzyme definition; reference REBASE version for external use |
| TF PWM | Caller provided matrix | Actual motif database (e.g. JASPAR) versions and entries; DNAKit does not have a built-in activity model |
| novelty/memorization | `native` Transparency rules | Actual reference library version, digest, filters and selected similarity method |
| synthesis-risk | `native` Transparency rules | DNAKit configuration and hit rules; must not be described as experimental success rate |

Each run only needs to reference the entries that are actually used; backends that are not called should not be counted against.

#### Original papers that user dinucleotide tables may cite {#dinucleotide-primary-references}

The table below only helps users trace historical sources and is not a list of DNAKit's built-in data, nor does it grant permission to copy or redistribute the data in these papers.

| Parameters | Original paper |
| --- | --- |
| Twist, Tilt, Roll, Shift, Slide, Rise | Perez et al. 2004, [PMID 15562006](https://pubmed.ncbi.nlm.nih.gov/15562006/) |
| Bend, Inclination, Major groove width, Minor groove width | Karas et al. 1996, [PMID 8996793](https://pubmed.ncbi.nlm.nih.gov/8996793/) |
| Direction | Shpigelman et al. 1993, [PMID 8402210](https://pubmed.ncbi.nlm.nih.gov/8402210/) |
| Propeller twist | Gorin et al. 1995, [PMID 7897660](https://pubmed.ncbi.nlm.nih.gov/7897660/) |
| Persistence length | Hogan & Austin 1987, [PMID 3627268](https://pubmed.ncbi.nlm.nih.gov/3627268/) |
| Stacking energy | Sponer et al. 1997, [PMID 9199773](https://pubmed.ncbi.nlm.nih.gov/9199773/) |
| Free energy | Sugimoto et al. 1996, [PMID 8948641](https://pubmed.ncbi.nlm.nih.gov/8948641/) |

### Websites, databases and public data {#websites-and-data}

| Resources | Purpose of this project | Website |
| --- | --- | --- |
| DiProDB | User-checkable DNA dinucleotide parameters and original paper index; DNAKit does not contain their values | [Database paper](https://doi.org/10.1093/nar/gkn597), [Public parameter table](https://diprodb.fli-leibniz.de/ShowTable.php) |
| NCBI Genetic Codes | Standard genetic code table 1 | [NCBI Standard Genetic Code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) |
| NCBI Datasets | Classification, assembly, gene, virus query and data package | [REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/), [Genome download](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/), [Virus package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) |
| NCBI Entrez, BLAST, ClinVar, dbSNP, GEO | accession/project/sample/variant/expression/literature query, BLAST jobs and public file downloads | [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/), [BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html), [GEO access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html) |
| Ensembl REST | Coordinate sequences, transcripts, regions/variants/regulations, homologous and comparative genomes | [Ensembl REST endpoints](https://rest.ensembl.org/) |
| European Nucleotide Archive | Study/Sample/Experiment/Run/Analysis Metadata and public files | [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api) |
| ENCODE Portal | Epigenomic experiments, peaks, signals and public documents | [ENCODE REST API](https://www.encodeproject.org/help/rest-api/) |
| UCSC Genome Browser | Chromosomes, coordinate sequences, annotations/duplicates/conserved tracks and download directories | [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html) |
| IDT and Sigma-Aldrich | OD260/oligonucleotide quantitative instructions and public ACGT calculation examples | [IDT quantitative instructions](https://sg.idtdna.com/page/support-and-education/decoded-plus/oligo-quantification-getting-it-right), [Sigma-Aldrich parameters and calculation examples](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/pcr/quantitation-of-oligos) |
| NUPACK | Analysis interface, model conditions and permission boundaries of conditional secondary structure backend | [Analysis documentation](https://docs.nupack.org/analysis/), [Utility functions](https://docs.nupack.org/utilities/), [Model documentation](https://docs.nupack.org/model/), [Download and licensing](https://www.nupack.org/download/overview) |
| 3DNA/DSSR | Three-dimensional structure parameter definition and DSSR JSON parsing basis | [3DNA parameter representation](https://x3dna.org/highlights/schematic-diagrams-of-base-pair-parameters), [DSSR JSON document](https://x3dna.org/highlights/dssr-output-in-json-format) |
| RCSB PDB | DNA three-dimensional structure records used for document verification | [1BNA](https://www.rcsb.org/structure/1BNA), [1AC7](https://www.rcsb.org/structure/1AC7), [139D](https://www.rcsb.org/structure/139D) |

## Third Party Statement {#third-party-notices}

Last review: 2026-09-01. This page is a construction compliance checklist, not legal advice.

DNAKit itself is licensed under the MIT license. Third-party packages, external programs, databases, papers, and user-provided parameter sheets retain their respective terms; citing a paper does not constitute a license to redistribute its software or data.

### Distribution Boundary {#distribution-boundary}

DNAKit's wheel and sdist do not include DiProDB numerical tables, Primer3, NUPACK, DSSR/3DNA, Dashing, BLAST, MMseqs2, sourmash, REBASE, JASPAR, DNA base-model checkpoints, or their databases, nor do they automatically download external programs. When a CLI adapter is present, the user must still install it separately and provide an explicit path. The representation API downloads permitted models to local `ckpt/` only after an explicit call; users download and place the new 27 task checkpoints from the [shared Google Drive folder](https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing). Neither path redistributes weights inside DNAKit.

Python dependencies in `pyproject.toml` are resolved by the installer as standalone distribution packages and are not copied into the DNAKit wheel. NCBI, Ensembl, ENA, ENCODE, and UCSC data are also not distributed with wheel/sdist; the adapter only accesses the public interface when explicitly called by the user. The main permissions for direct dependencies are as follows; locked versions and transitive dependencies must still be checked before actual release.

| Scope | Direct Dependencies | License Flag |
| --- | --- | --- |
| Core | PyYAML, Rich, Typer, tomli (Python <3.11) | MIT |
| `viz` | CairoSVG; Pillow | LGPL-3.0-or-later; MIT-CMU |
| `io` | PyArrow | Apache-2.0 |
| `validation` | Biopython | Biopython License Agreement |
| `external-tools` | ToolUniverse | Apache-2.0; subject to the installed release and transitive dependencies |
| `agent` | FastMCP | Apache-2.0; subject to the installed release and transitive dependencies |
| `neural` | Hugging Face Hub, NumPy, scikit-learn, PyTorch, Transformers | Apache-2.0, BSD series, Apache-2.0; subject to the actual installed version |
| `neural-caduceus` | mamba-ssm | Apache-2.0; subject to the actual installed version |
| `neural-enformer` | encoder-pytorch | MIT; subject to the actual installed version |
| `neural-evo2` | evo2 and its transitive dependencies | The actual installed version and official repository shall prevail |
| `docs` | MkDocs; Material for MkDocs; mkdocstrings; nbmake | BSD-2-Clause; MIT; ISC; Apache-2.0 |
| `dev` | build, Hypothesis, mypy, pytest, pytest-cov, Ruff, twine, types-PyYAML | MIT, MPL-2.0, MIT, MIT, MIT, MIT, Apache-2.0, Apache-2.0 |

### External tools and user data {#external-tools-and-user-data}

| Resources | DNAKit processing methods | Must confirm before use |
| --- | --- | --- |
| ToolUniverse | Installed as an optional Python dependency and loaded only for explicitly called allowlisted functions; source is not vendored | Code is Apache-2.0; remote Ensembl, NCBI, gnomAD, and EMBL-EBI services retain their own terms and availability limits. See the [official repository](https://github.com/mims-harvard/ToolUniverse). |
| Primer3 CLI | Only called via explicit `primer3_core`, `oligotm`, `ntthal` paths; no `primer3-py` or binary | Official source tag GPL-2.0-or-later; check [ repository ](https://github.com/primer3-org/primer3) with actual version |
| NUPACK | Obtained and installed separately by user; not downloaded, not packaged, not provided as an online service | Separate terms/subscriptions, fees may apply, academic and commercial conditions differ; see [official license page](https://www.nupack.org/download/overview) |
| DSSR/3DNA | Separate user installation; DNAKit only parses explicit output | Academic and commercial permissions are different; free academic Basic is still subject to license; see [Official Notes](https://home.x3dna.org/highlights/x3dna-dssr-is-funded-and-dssr-basic-academic-is-free) |
| Dashing | User provides explicit executable file; distribution package does not contain source code/binary | GPL-3.0; see [official repository](https://github.com/dnbaker/dashing) |
| DiProDB or other dinucleotide tables | DNAKit only provides JSON schema/loader, without any DiProDB values | Users must have access to and record the true source, version and SHA-256; see [DiProDB](https://diprodb.fli-leibniz.de/) |
| Other CLI/database | No automatic installation or packaging | Check program and database terms separately; accessible does not mean redistributable |
| NCBI/Ensembl/ENA/ENCODE/UCSC public interface and data | Only download according to user query, save source and checksum; do not embed remote data into the distribution package | Check the use policy, attribution, privacy and redistribution terms of each provider; controlled data such as dbGaP still requires formal authorization |
| Hugging Face DNA model checkpoint | After the user selects the model, download it to `ckpt/` and record the source manifest; wheel/sdist does not contain weights | DNABERT-2, NTv2, HyenaDNA, Caduceus, GROVER, LucaOne, GENERator, Enformer and Evo 2 may have different terms; check item by item [ Official link in the sequence characterization page ](api/features/08_fingerprints.md#neural-representations) |
| AlphaGenome | Only loaded after the user accepts the access terms, obtains checkpoint, and installs the official research code | The checkpoint page indicates restricted access and non-commercial model terms; check [official repository](https://github.com/google-deepmind/alphagenome_research) and current model page |
| JanusDNA | Download and verify MD5 from the official files of Harvard Dataverse; users are also required to provide official source code environments | Check [official repositories](https://github.com/Qihao-Duan/JanusDNA), [Dataverse DOI](https://doi.org/10.7910/DVN/HDT0RN) and their code/data terms |

Transformers for DNABERT-2, NTv2, HyenaDNA, Caduceus, LucaOne and GENERator
The checkpoint contains the Python code required to load. DNAKit denies execution by default; only the user reviews the source and
Enabled after explicitly setting `allow_remote_code=True`. This confirmation cannot replace code auditing and dependency locking
or license check.

Process isolation and "installed by user" reduce bundle distribution risks but do not automatically justify all uses. Commercial deployments, hosting services, source/binary redistributions, or archiving of journal supplementary materials should be re-examined against the actual submission.

The complete English list is saved with the source code and installation package in the "Third Party Statement" section of the root directory `DISCLAIMER.md`; for scientific methods, papers and websites, see [Acknowledgments and Primary Sources](#acknowledgements), and for general risks of use, see [Disclaimer](#disclaimer).

## Disclaimer {#disclaimer}

DNAKit is software for research and teaching purposes, not a medical device, and is not intended for clinical diagnosis, treatment, patient management, experimental safety decision-making, or other high-risk uses. Output may be incomplete, inaccurate, or inappropriate for a particular experiment; users must independently verify calculations, parameters, input data, and experimental conclusions.

DNAKit is provided "as is" under the MIT License, without warranties express or implied, and liability is limited to the maximum extent permitted by applicable law. The `LICENSE` in the root directory of the repository is the official license text; this page only describes the boundaries of use and does not replace or modify the MIT license, nor does it constitute legal, medical, or compliance advice.

"This project is mainly developed for academic research" is not an exemption clause. MIT allows both research and commercial use, so DNAKit cannot be described as "non-commercial use only" at the same time, unless it is changed to another reviewed license in the future.

Users are responsible for:

- Confirm that the data, parameter tables, databases and external programs it provides have applicable access, use and redistribution rights;
- Do not use ClinVar, VEP, BLAST or other remote provider records directly as clinical diagnosis, source determination or experimental safety conclusion;
- Comply with third-party licenses, subscriptions, attribution/citations, privacy, ethics, biosecurity, export controls and local laws;
- Disclose actual methods, versions, parameters, limitations and sources in the paper and software;
- Independent verification and necessary professional review before experimentation, publication, deployment or commercial use.

DNAKit does not guarantee that third-party resources are free, continuously available, suitable for commercial use, or permit redistribution. Terms and fees may change and you should check the official page when using. For specific distribution boundaries, see [Third Party Statement](#third-party-notices).

## License {#licenses}

### Project License {#project-license}

DNAKit itself adopts the MIT license, and the copyright is Pengsen Ma. The complete text can be found in `LICENSE` in the root directory of the repository. The SPDX identifier is `MIT`. The installation package contains `DISCLAIMER.md` which incorporates third-party notices; for the web version, see [Third-Party Notices](#third-party-notices) and [Disclaimer](#disclaimer).

The current version of the project is `0.1.3`, which is still in the early development stage; dependency compatibility and external backends still need to be reviewed according to the actual environment.

A license for a third-party package does not automatically become a license for DNAKit.

### Implement tag {#implementation-labels}

| Label | Meaning | Current Example |
| --- | --- | --- |
| `native` | DNAKit own objects, combinational logic, or explicit simple deterministic logic | Core objects, auditing, cache, workflow, scorecard |
| `reimplementation` | Independently implemented based on public algorithms, papers or format specifications | alignment, SantaLucia NN, MinHash, format codec |
| `adapter` | Calling independently installed packages, CLIs or databases | Primer3CLIAdapter, PyArrow Parquet, DashingAdapter, external CLI metadata handles |
| `novel` | New methods confirmed by definition, retrieval, baseline and ablation | There are currently no confirmed projects |

Mixed/multi-scale fingerprints are currently `native` combinations and are not marked as `novel`.


### Current local dependent license snapshot {#dependency-license-snapshot}

The following table is from the installed distribution metadata of the `dnakit-dev` environment or a third party in the repository `LICENSE` and is for local audit purposes and does not substitute for legal advice.

| Components | Local Version/Status | License Tags | DNAKit Boundaries |
| --- | --- | --- | --- |
| PyYAML | 6.0.3 | MIT | Core configuration parsing |
| Rich | 15.0.0 | MIT | CLI/Script Progress Display |
| Typer | 0.27.1 | MIT | CLI |
| tomli | 2.4.1 | MIT | Python <3.11 Conditional TOML parsing |
| CairoSVG | 2.9.0 | LGPL-3.0-or-later | `viz` extra; SVG to PNG/PDF |
| Pillow | 12.3.0 | MIT-CMU | `viz` extra; TIFF/resolution metadata |
| PyArrow | 25.0.1 | Apache-2.0 | `io` extra; Parquet table reading and writing |
| NumPy | 2.2.6 | BSD-3-Clause | Local rep matrix and clustering tests |
| scikit-learn | 1.7.2 | BSD-3-Clause | Local k-means/PCA/silhouette testing |
| Biopython | 1.88 | `LicenseRef-Biopython-License-Agreement` | Differential control for validation extra only |
| MkDocs | 1.6.1 | BSD-2-Clause | Documentation Development |
| MkDocs Material | 9.7.7 | MIT | Document Development |
| mkdocstrings | 0.30.1 | ISC | API documentation generation |
| nbmake | 1.5.5 | Apache-2.0 | Notebook Access Control |
| Dashing local copy | Out-of-package third-party directory under project root | GPL-3.0 (this directory LICENSE); runtime adapter record `GPL-3.0-only` | Excluded from wheel/sdist; explicit adapter does not automatically select or package this copy, real scientific difference is not completed |

Primer3 is no longer a Python dependency and is not in any extras. Currently `dnakit-dev` is not installed yet
Complete PyTorch/Transformers model inference stack; GROVER real checkpoint smoke in standalone compatibility
Model environment execution, the above table cannot be interpreted to mean that all 11 models have been verified. The above table only represents the currently directly dependent
Local metadata snapshot; the license list must be regenerated after version upgrade, repackaging or changing distribution method.
For a complete list of direct dependencies, see [Third Party Statement](#third-party-notices).

### License prompt for registered external CLI handle {#registered-cli-license-notices}

These values are metadata in `BackendInfo.license_expression` to prompt the user to review, and are not a legal conclusion about DNAKit's suitability for third-party licenses. The path discovery of BLAST/MMseqs2/sourmash does not execute the program, and only the user explicitly calls `.version()` to run the restricted version command; Dashing also has a scientific computing adapter that must explicitly provide the executable path and call it explicitly.

| Backend ID | License hint for DNAKit records | Current capabilities |
| --- | --- | --- |
| `blastn` | `LicenseRef-NCBI` | Passive path positioning, explicit version query |
| `mmseqs2` | `GPL-3.0-or-later` | Passive path positioning, explicit version query |
| `sourmash` | `BSD-3-Clause` | Passive path positioning, explicit version query |
| `dashing` | `GPL-3.0-only` | registry passive path/version handle; exact/HLL Jaccard/Top-k adapter with explicit path; project root third-party copy is not included in the release archive |

Before actual use, the respective versions, license files and usage conditions of the binary/source code and database provided by the user must prevail.

### Restricted backends and databases {#restricted-backends-and-databases}

| Project | Currently Processing |
| --- | --- |
| NUPACK | Provides passive detection and explicit adapter, but does not automatically install/download, does not package, and does not serve as a web backend; independent subscription terms may be charged, the current environment is not available and the real value difference is subject to licensing and installation conditions |
| Primer3 | Users install the CLI separately and explicitly provide `oligotm`, `ntthal` or `primer3_core` paths; DNAKit no longer relies on `primer3-py`, does not search according to `PATH`, and does not package; the license prompt is `GPL-2.0-or-later`, and must still be reviewed according to the actual version |
| DSSR/3DNA | Users obtain a separate license and install it; DNAKit only parses the output provided by the user, and the academic free terms cannot be extrapolated to commercial free or arbitrary redistribution |
| Dashing | Local GPLv3 source code/binary does not enter the release archive; explicit path adapter can do bounded exact/HLL Jaccard and Top-k, but there is no real scientific difference and the local copy will not be automatically selected |
| BLAST/MMseqs2/sourmash | No automatic installation or download; passive metadata/explicit version handles registered, but no search, clustering or sketch executors |
| RepeatMasker/TRF/FIMO | No automatic installation or download; currently no unified adapter |
| Databases such as REBASE/JASPAR | Not redistributed with packages; provided by user under license and recorded version, date, filter and checksum |
| DiProDB/dinucleotide parameter table | No built-in DiProDB values; users provide the 15×16 table they have permission to use according to the fixed JSON schema, and the result records the table name, version, source statement and SHA-256 |
| DNA basic model checkpoint | is only downloaded to `ckpt/` when the user explicitly calls the rep API, and is not distributed with the package; checkpoint/source code terms are checked model by model, and remote code is rejected by default; currently only GROVER real smoke |
