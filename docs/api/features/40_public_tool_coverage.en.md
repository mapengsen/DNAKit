# Public DNA tool coverage

`new_vision/2026-10-10/public_dna_tools.md` is a dated source survey. DNAKit
keeps the project-by-project cross-reference in
`new_vision/2026-10-10/public_dna_tools_coverage.{json,csv,md}` and exposes the
same records through `dnakit.integrations`.

> `new_vision/` is a local audit directory kept by the author. It is not part of the repository, so the paths below point to files you will not find in a checkout. The results are summarised on this page and available through `dnakit.integrations`.

The source text declares 282 spreadsheet records. The Markdown available in
this repository contains 425 explicit mentions, which collapse to 320
conservatively normalized source records. Classification, citation and
provider-specific entries remain evidence; only spelling aliases are merged.
The report also preserves the source's declared 169 conditionally integrable
records as a separate figure.

Every in-scope candidate capability has a formal public callable name in the
form `dnakit.<module>.<function>()`. `dnakit_candidates` is an internal staging
namespace used by the implementation registry and is never the public API.
The formal entry may be a native implementation, a bounded facade, or an
explicit dependency-checked adapter; the implementation status is recorded
separately below.

## Query the coverage matrix

```python
from dnakit.integrations import (
    describe_third_party_tool,
    list_third_party_tools,
    list_unified_capabilities,
    integration_inventory,
    run_unified_capability,
    third_party_coverage,
    third_party_route,
    unified_capability_summary,
)

print(third_party_coverage()["summary"])
print(integration_inventory())
print([item["name"] for item in list_third_party_tools("biopython")])
print(describe_third_party_tool("Biopython")["citation_dois"])
print(third_party_route("Cas-OFFinder")["api"])
print(unified_capability_summary()["count"])
print(list_unified_capabilities(module="sequence")[0]["api"])
print(run_unified_capability("CAP-0187", {"input": "@r1\\nACGT\\n+\\nIIII\\n"}))
```

Routes are explicit. `analysis_adapter` uses
`dnakit.integrations.run_analysis`; `mcp_registry` uses the MCP bridge;
`skill_registry` uses the bundled Skill/workflow layer; and
`pip_adapter` uses an explicit optional Python backend such as Biopython,
 pydna, seqme, DNA Chisel, Biotite, Edlib, PyRanges, pyfastx, dnaio, mappy, pyBigWig, pybedtools,
 gffutils, pySBOL3, scadnano, polars-bio, DnaCauldron, DNA Features Viewer, gtars, GeneBlocks or sourmash; the package and any model extras remain separate.
`native_facade` maps a surveyed tool to a focused DNAKit implementation, such
as the BWA/Bowtie2/SAMtools/BCFtools/SPAdes/Picard NGS facades, the NCBI
BLAST+/minimap2/BEDTools alignment facades, the MACS3 peak facade, SeqKit
FASTA/FASTQ functions, the FastTree tree facade, the SPOA consensus and the
vcfanno/SURVIVOR variant
facades and the NanoFilt/Filtlong/NanoPlot/Porechop/Trimmomatic/AfterQC/SeqMagick/GATK4/Shasta/Merqury/BUSCO/DeepVariant/cuteSV/Sniffles2/DELLY/HTSlib/seqtk/WhatsHap/BamTools/VSEARCH/BEDOPS/Flye/hifiasm/MEGAHIT/GEMINI/Ensembl VEP external
facades; it does not mean that the upstream command's complete CLI has been
reimplemented.
`project_evidence_only` or `no_local_match` is an external-adapter boundary,
not a claim that upstream source was copied or that every upstream function is
implemented. All external programs, databases, models and credentials remain
separate requirements.

`local_match` means that a name, repository, dependency, Skill notice, MCP
entry or analysis entry was found in DNAKit. It is not scientific validation.

Call `integration_inventory()` for current surface counts. It reports Python
functions, Agent-compatible Python functions, MCP server recipes, registered
Skills, surveyed projects and candidate schemas as separate units; these
categories should not be added together as if they were interchangeable tools.
The current raw surface sum is **2,432** (1,826 public Python functions + 32
MCP server recipes + 574 registered Skills). This is a surface count rather
than a deduplicated biological-tool count because one capability can appear in
several surfaces.
The inventory also includes `unified_candidate_statuses`, separating
implemented, partial and adapter candidates.
It also keeps the survey's kind and integration-status counts, so the 103
surveyed MCP records are not mistaken for 103 locally connected servers.
At this snapshot the 320 survey records contain 62 `native_facade`, 2
`analysis_adapter`, 26 `pip_adapter`, 13 `mcp_registry`, 2 `skill_registry`,
44 `project_evidence_only` and 171 `no_local_match` records.

The shipped `unified_capabilities.json` records 837 formal candidate API names
and their input/output JSON Schemas. `list_unified_capabilities()` and
`describe_unified_capability()` expose that catalog through the same Python and
MCP discovery layer. Thirty-two project records currently have conservative
name/schema matches to one or more capabilities; these are review links, not
proof that an upstream project's full API is implemented.

`run_unified_capability()` executes a catalog capability through the existing
Agent/MCP registry. It applies the registered JSON argument coercion, output
serialization limits and side-effect confirmation rules. It never installs an
upstream package or silently executes an unverified project.

The permissively licensed pip entries have bounded adapters:

- Existing `dnakit.bioinformatics` format/alignment functions use Biopython
  with `dnakit[bio-formats]`, reads functions use pysam with
  `dnakit[bio-reads]`, and cloning functions use pydna with
  `dnakit[bio-assembly]`. These are focused DNAKit facades, not complete
  upstream namespace re-exports.

- `dnakit.integrations.beta_diversity_scikit_bio()` and
  `dnakit.integrations.permanova_scikit_bio()` expose bounded scikit-bio
  statistics. Install `dnakit[bio-skibio]`; the full scikit-bio namespace is
  not re-exported.
- `dnakit.integrations.evaluate_seqme()` constructs allow-listed
  `seqme.metrics` classes and serializes the evaluation table. Install
  `dnakit[bio-seqme]`; model extras and weights remain separate.
- `dnakit.integrations.optimize_with_dnachisel()` constructs an allow-listed
  set of DNA Chisel specifications and returns the optimized sequence. Install
  `dnakit[bio-dnachisel]`; it does not copy DNA Chisel source code.
- `dnakit.integrations.fetch_pyfaidx_region()` provides bounded FASTA interval
  reads with `dnakit[bio-pyfaidx]`; `summarize_vcf_vcfpy()` and
  `summarize_vcf_cyvcf2()` provide bounded VCF summaries with their respective
  optional extras.
- `dnakit.integrations.align_biotite()`, `align_edlib()` and the two PyRanges
  BED adapters provide bounded pairwise alignment, edit-distance alignment and
  interval operations through `bio-biotite`, `bio-edlib` and `bio-pyranges`.
- `dnakit.integrations.summarize_pyfastx()` provides bounded FASTA/FASTQ
  aggregation through `bio-pyfastx`; `summarize_pybigwig_region()` provides
  bounded local bigWig statistics through `bio-pybigwig`. Neither writes input
  files, builds a mutable index, or accesses remote URLs.
- `dnakit.integrations.summarize_dnaio_fastx()` uses `bio-dnaio` for bounded,
  read-only FASTA/FASTQ stream summaries; `dnakit.integrations.map_mappy_sequences()`
  uses `bio-mappy` and minimap2/mappy for bounded local mapping. Neither writes an
  index or accepts arbitrary upstream command-line flags; missing extras produce a
  structured dependency status.
- `dnakit.integrations.intersect_pybedtools_bed()` provides bounded BED
  intersection through `bio-pybedtools`; it still requires a separately
  installed `bedtools` executable and accepts no arbitrary command-line flags.
- `dnakit.integrations.query_gffutils_database()` performs read-only queries
  against an existing gffutils `FeatureDB` through `bio-gffutils`; it does not
  create databases or execute arbitrary SQL.
- `dnakit.integrations.summarize_sbol3_document()` inventories a local SBOL3 RDF
  document through `bio-pysbol3`; it does not resolve remote RDF, mutate objects,
  or write files.
- `dnakit.integrations.summarize_scadnano_design()` inventories a local scadnano
  `.sc` DNA nanostructure design through `bio-scadnano`; it does not edit, render,
  or download design assets.
- `dnakit.integrations.overlap_polars_bio_bed()` performs bounded overlap of two
  local BED/CSV interval files through `bio-polarsbio`; it does not execute
  arbitrary SQL or access remote URLs.
- `dnakit.integrations.simulate_dnacauldron_assembly()` runs bounded Gibson or
  Type IIS assembly simulation through `bio-dnacauldron`, returning predicted
  constructs and backend warnings/flaws without inferring a wet-lab result.
- `dnakit.integrations.summarize_dna_features_viewer()` uses
  `bio-dna-features-viewer` to inventory local GenBank/GFF features without
  writing a plot file or re-exporting the full DNA Features Viewer or Biopython API.
- `dnakit.integrations.summarize_gtars_regions()` uses `bio-gtars` to inventory
  a local BED `RegionSet`, returning region count, covered bases and chromosome
  counts; it rejects remote URLs, cache construction and arbitrary interval expressions.
- `dnakit.integrations.diff_geneblocks_sequences()` uses `bio-geneblocks` for a
  deterministic diff of two bounded DNA strings; it does not invoke BLAST,
  plotting or feature transfer.
- `dnakit.integrations.compare_sourmash_sequences()` uses `bio-sourmash` for an
  in-memory MinHash comparison of two bounded A/C/G/T sequences; it does not
  access databases, network services or taxonomy indexes.
- `dnakit.sequencing_qc.run_cutadapt()` invokes a caller-provided Cutadapt
  executable through a fixed option allowlist for single-end or paired FASTQ
  trimming; it requires explicit output paths and never overwrites them.
- `dnakit.sequencing_qc.run_fastp()` invokes a caller-provided fastp executable
  through a fixed option allowlist for bounded FASTQ filtering and reports; an
  external process requires explicit Agent/MCP authorization.
- `dnakit.sequencing_qc.run_fastqc()` and `run_multiqc()` invoke caller-provided
  FastQC and MultiQC executables and monitor HTML/ZIP or aggregated HTML
  artifacts; Java, plugins and GPL backends are not bundled.

UGENE, iLearnPlus, iFeatureOmega and `ga` from the screenshot are kept as
explicit external boundaries. They may be cited or called after the user
installs the corresponding backend, but they are not reported as copied into
DNAKit core. `third_party_route()` returns the external entry point together
with the available license and citation evidence.

## Citation and license evidence

The survey supplies paper identifiers for these central tools:

| Tool | Citation | DNAKit boundary |
| --- | --- | --- |
| Biopython | Cock et al. (2009), [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) | `dnakit.bioinformatics` adapters; upstream license and notices apply. |
| SeqKit | Shen, Sipos & Zhao (2024), SeqKit2, [10.1002/imt2.191](https://doi.org/10.1002/imt2.191); the earlier SeqKit paper is [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) | Separately installed CLI through an external-analysis adapter. |
| scikit-bio | Aton et al. (2026), [10.1038/s41592-025-02981-z](https://doi.org/10.1038/s41592-025-02981-z) | Optional Python dependency; retain its BSD-3-Clause notice. |
| UGENE | Okonechnikov et al. (2012), [10.1093/bioinformatics/bts091](https://doi.org/10.1093/bioinformatics/bts091) | External GPLv2 program; not copied into the core package. |
| pydna | Pereira et al. (2015), [10.1186/s12859-015-0544-x](https://doi.org/10.1186/s12859-015-0544-x) | Optional DNA-assembly dependency; retain its BSD notice. |
| DNA Chisel | Zulkower & Rosser (2020), [10.1093/bioinformatics/btaa558](https://doi.org/10.1093/bioinformatics/btaa558) | Optional MIT dependency; model and data terms are separate. |
| seqme | Møller-Larsen et al. (2025), [arXiv:2511.04239](https://arxiv.org/abs/2511.04239) | Optional BSD-3-Clause adapter; model extras and weights are separate. |
| Biotite | Kunzmann & Hamacher (2018), [10.1186/s12859-018-2367-z](https://doi.org/10.1186/s12859-018-2367-z) | Optional BSD-3-Clause adapter; bounded nucleotide alignment only. |
| Edlib | Šošić & Šikić (2017), [10.1093/bioinformatics/btw753](https://doi.org/10.1093/bioinformatics/btw753) | Optional MIT adapter; bounded edit-distance alignment only. |
| PyRanges | Stovner & Sætrom (2020), [10.1093/bioinformatics/btz615](https://doi.org/10.1093/bioinformatics/btz615) | Optional MIT adapter; bounded BED interval operations only. |
| pyfastx | Du et al. (2021), [10.1093/bib/bbaa368](https://doi.org/10.1093/bib/bbaa368) | Optional MIT adapter; bounded FASTA/FASTQ summary only. |
| pyfaidx | Shirley et al. (2015), [10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1); source code BSD-3-Clause | Optional BSD-3-Clause adapter; bounded read-only FASTA interval fetching. |
| dnaio | dnaio software archive, [10.5281/zenodo.10548864](https://doi.org/10.5281/zenodo.10548864); MIT | Optional MIT adapter; bounded read-only FASTA/FASTQ stream summary. |
| mappy/minimap2 | Li (2018), [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191); [official Python binding README](https://github.com/lh3/minimap2/blob/master/python/README.rst); MIT | Optional MIT adapter; bounded local mapping with allow-listed presets. |
| pyBigWig | pyBigWig software release, [10.5281/zenodo.45238](https://doi.org/10.5281/zenodo.45238) | Optional MIT adapter; bounded local bigWig statistics only. |
| pybedtools | Dale et al. (2011), [10.1093/bioinformatics/btr539](https://doi.org/10.1093/bioinformatics/btr539); BEDTools: [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) | Optional MIT wrapper; a separately installed BEDTools executable remains required. |
| gffutils | gffutils documentation and daler contributors | Optional MIT adapter; read-only queries against an existing FeatureDB; the survey supplied no DOI. |
| pySBOL3 | Mitchell, Beal & Bartley (2022), [10.1021/acssynbio.2c00164](https://doi.org/10.1021/acssynbio.2c00164) | Optional MIT adapter; bounded read-only SBOL3 inventory. |
| scadnano | Doty, Lee & Stérin (2020), [arXiv:2005.11841](https://arxiv.org/abs/2005.11841) | Optional MIT adapter; bounded read-only `.sc` design inventory. |
| polars-bio | polars-bio, [10.1093/bioinformatics/btaf640](https://doi.org/10.1093/bioinformatics/btaf640) | Optional Apache-2.0 adapter; bounded local interval overlap. |
| DnaCauldron | DnaCauldron documentation and Edinburgh Genome Foundry contributors; the survey supplied no DOI | Optional MIT adapter; bounded Gibson/Type IIS simulation only. |
| DNA Features Viewer | [official documentation](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/) and [Edinburgh Genome Foundry repository](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer); the survey supplied no DOI | Optional MIT adapter; bounded read-only GenBank/GFF feature inventory. |
| gtars | [official documentation](https://docs.bedbase.org/gtars/python-overview/) and [databio repository](https://github.com/databio/gtars); the survey supplied no DOI | Optional BSD-2-Clause adapter; bounded read-only local BED RegionSet inventory. |
| GeneBlocks | [official documentation](https://edinburgh-genome-foundry.github.io/Geneblocks/) and [Edinburgh Genome Foundry repository](https://github.com/Edinburgh-Genome-Foundry/Geneblocks); the survey supplied no DOI | Optional MIT adapter; deterministic bounded sequence diff. |
| sourmash | [JOSS paper](https://doi.org/10.21105/joss.06830) and [official documentation](https://sourmash.readthedocs.io/); BSD-3-Clause | Optional BSD-3-Clause adapter; local bounded MinHash comparison. |
| Cutadapt | Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200); [official guide](https://cutadapt.readthedocs.io/en/stable/guide.html); MIT | Bounded external executable adapter for FASTQ adapter trimming. |
| fastp | Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); [official repository](https://github.com/OpenGene/fastp); MIT | Bounded external executable adapter for FASTQ preprocessing and QC. |
| FastQC | [official project page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); GPL-3.0; the survey supplied no DOI | Bounded external executable adapter for HTML/ZIP read-QC reports. |
| MultiQC | Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354); [official repository](https://github.com/MultiQC/MultiQC); GPL-3.0-or-later | Bounded external executable adapter for aggregated QC reports. |

The ten projects shown in the screenshot have the following explicit routes:

| Project | DNAKit route | Current boundary |
| --- | --- | --- |
| Biopython | `dnakit.bioinformatics` | Existing format, reads, alignment and motif entry points; not a full namespace re-export |
| SeqKit | `dnakit.fasta_fastq` | Existing FASTA/FASTQ operations; other CLI commands still require a separate install |
| scikit-bio | `dnakit.integrations.beta_diversity_scikit_bio`, `permanova_scikit_bio` | Optional BSD-3-Clause backend; bounded statistics only |
| UGENE | `dnakit.external_adapters.run` | GPLv2 external program; not copied into the core package |
| iLearnPlus | `dnakit.external_adapters.run` | Redistribution scope still requires author/license confirmation |
| iFeatureOmega | `dnakit.external_adapters.run` | Redistribution scope still requires author/license confirmation |
| `ga` | `dnakit.external_adapters.run` | GPLv3 general genetic-algorithm package, not a native DNAKit DNA API |
| seqme | `dnakit.integrations.evaluate_seqme` | Optional BSD-3-Clause backend; model weights require separate review |
| pydna | `dnakit.bioinformatics.cloning` | Existing construct entry points; not a full pydna API re-export |
| DNA Chisel | `dnakit.integrations.optimize_with_dnachisel` | Optional MIT backend; only allow-listed specification classes are constructible |
| Biotite | `dnakit.integrations.align_biotite` | Optional BSD-3-Clause backend; bounded optimal nucleotide pairwise alignment |
| Edlib | `dnakit.integrations.align_edlib` | Optional MIT backend; bounded edit distance, locations or CIGAR path |
| PyRanges | `dnakit.integrations.summarize_pyranges_bed`, `merge_pyranges_bed` | Optional MIT backend; bounded BED interval summary and merge, rejecting over-limit input |
| pyfastx | `dnakit.integrations.summarize_pyfastx` | Optional MIT backend; bounded FASTA/FASTQ summary, rejecting over-limit input |
| dnaio | `dnakit.integrations.summarize_dnaio_fastx` | Optional MIT backend; read-only FASTA/FASTQ stream summary, rejecting over-limit input |
| mappy/minimap2 | `dnakit.integrations.map_mappy_sequences` | Optional MIT backend; bounded local sequence mapping with allow-listed presets |
| pyBigWig | `dnakit.integrations.summarize_pybigwig_region` | Optional MIT backend; bounded local bigWig region statistics |
| pybedtools | `dnakit.integrations.intersect_pybedtools_bed` | Optional MIT Python package; bounded BED intersection, external bedtools still required |
| gffutils | `dnakit.integrations.query_gffutils_database` | Optional MIT backend; read-only queries against an existing FeatureDB, without database creation or mutation |
| pySBOL3 | `dnakit.integrations.summarize_sbol3_document` | Optional MIT backend; read-only inventory of a local SBOL3 document |
| scadnano | `dnakit.integrations.summarize_scadnano_design` | Optional MIT backend; read-only inventory of a local `.sc` design |
| polars-bio | `dnakit.integrations.overlap_polars_bio_bed` | Optional Apache-2.0 backend; bounded local interval overlap |
| DnaCauldron | `dnakit.integrations.simulate_dnacauldron_assembly` | Optional MIT backend; bounded Gibson/Type IIS assembly simulation |
| DnaFeaturesViewer | `dnakit.integrations.summarize_dna_features_viewer` | Optional MIT backend; bounded read-only GenBank/GFF feature inventory |
| gtars | `dnakit.integrations.summarize_gtars_regions` | Optional BSD-2-Clause backend; bounded read-only BED RegionSet inventory |
| GeneBlocks | `dnakit.integrations.diff_geneblocks_sequences` | Optional MIT backend; bounded deterministic DNA sequence diff |
| sourmash | `dnakit.integrations.compare_sourmash_sequences` | Optional BSD-3-Clause backend; bounded local MinHash comparison |
| Cutadapt | `dnakit.sequencing_qc.run_cutadapt` | MIT external executable adapter; fixed-option FASTQ adapter trimming |
| fastp | `dnakit.sequencing_qc.run_fastp` | MIT external executable adapter; fixed-option FASTQ preprocessing and QC |
| FastQC | `dnakit.sequencing_qc.run_fastqc` | GPL-3.0 external executable adapter; monitors HTML/ZIP QC reports |
| MultiQC | `dnakit.sequencing_qc.run_multiqc` | GPL-3.0-or-later external executable adapter; monitors aggregated HTML reports |
| pyfaidx | `dnakit.integrations.fetch_pyfaidx_region` | Optional BSD-3-Clause backend; read-only FASTA intervals |
| vcfpy | `dnakit.integrations.summarize_vcf_vcfpy` | Optional MIT backend; bounded VCF summaries |
| cyvcf2 | `dnakit.integrations.summarize_vcf_cyvcf2` | Optional MIT backend; bounded VCF summaries |

`citation_status` is `missing_in_survey` when the source did not provide a DOI
or preferred citation. Resolve that gap from the upstream `CITATION.cff`, paper
or documentation before claiming use in a paper. A permissive code license
does not grant redistribution rights for a hosted API, model weights,
reference database or user data. Restricted, non-commercial and unverified
projects remain in the audit matrix and are not copied into DNAKit.

Regenerate the evidence after editing the dated source with
`new_vision/2026-10-10/build_public_tool_coverage.py`. The command is offline;
it does not verify current upstream licenses or execute scientific backends.
The implementation-pass file list is recorded in
`new_vision/2026-10-10/implementation_files_manifest.txt`.
