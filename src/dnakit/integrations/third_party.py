"""Discover the public DNA-tool coverage matrix and its unified routes.

The matrix is intentionally separate from the Skill catalog.  It records
third-party projects named by the dated survey, the license evidence stated
there, the corresponding DNAKit registry match and the safe next route.  A
``reference_only`` route never executes or downloads an upstream project.
"""

from __future__ import annotations

import copy
import json
import re
from collections import Counter
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError

_CATALOG = Path(__file__).with_name("third_party_tools.json")
_CAPABILITIES = Path(__file__).with_name("unified_capabilities.json")
_MCP_CATALOG = Path(__file__).with_name("mcp_servers.json")
_SKILL_CATALOG = Path(__file__).parents[3] / "registry/catalog.json"


def _name_key(value: str) -> str:
    """Normalize package spellings while keeping provider names distinct."""

    normalized = re.sub(r"[^a-z0-9]+", "-", value.casefold()).strip("-")
    return {
        "dna-chisel": "dnachisel",
        "dna-weaver": "dnaweaver",
        "dna-cauldron": "dnacauldron",
        "dnafeaturesviewer": "dna-features-viewer",
        "pysbol3": "sbol3",
    }.get(normalized, normalized)

_NATIVE_FACADES: dict[str, dict[str, Any]] = {
    "mafft": {
        "api": "dnakit.alignment.align_sequences",
        "functions": ["dnakit.alignment.align_sequences"],
        "sources": ["https://mafft.cbrc.jp/alignment/software/"],
        "citation_dois": ["10.1093/molbev/mst010"],
        "citation_references": ["Katoh & Standley (2013), DOI 10.1093/molbev/mst010"],
        "license_evidence": ["BSD"],
        "boundary": (
            "Runs a caller-installed MAFFT executable through the bounded multiple-alignment "
            "facade; the full MAFFT option surface is not re-exported."
        ),
    },
    "muscle": {
        "api": "dnakit.alignment.align_sequences",
        "functions": ["dnakit.alignment.align_sequences"],
        "sources": ["https://github.com/rcedgar/muscle"],
        "citation_dois": ["10.1038/s41467-022-34630-w"],
        "citation_references": [
            "Edgar (2022), DOI 10.1038/s41467-022-34630-w"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed MUSCLE 5 executable through the bounded multiple-"
            "alignment facade; the full MUSCLE command surface is not re-exported."
        ),
    },
    "macs3": {
        "api": "dnakit.ngs.call_chip_peaks",
        "functions": ["dnakit.ngs.call_chip_peaks"],
        "sources": ["https://github.com/macs3-project/MACS"],
        "citation_references": [
            "MACS3 project authors and contributors; the survey did not supply a DOI"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "Runs caller-installed MACS3 callpeak with fixed treatment/control BAM inputs and "
            "an explicit effective genome size; arbitrary MACS3 subcommands are not exposed."
        ),
    },
    "ncbi-blast": {
        "api": "dnakit.alignment.build_blast_database / dnakit.alignment.local_blast",
        "functions": [
            "dnakit.alignment.build_blast_database",
            "dnakit.alignment.local_blast",
        ],
        "sources": ["https://blast.ncbi.nlm.nih.gov/doc/blast-help/"],
        "citation_dois": ["10.1186/1471-2105-10-421"],
        "citation_references": [
            "Camacho et al. (2009), DOI 10.1186/1471-2105-10-421"
        ],
        "license_evidence": ["Public domain (NCBI BLAST+ software)"],
        "boundary": (
            "Builds a local nucleotide BLAST database and runs fixed-format blastn XML; "
            "the caller supplies the BLAST+ executables and database, and arbitrary BLAST+ "
            "subcommands are not re-exported."
        ),
    },
    "minimap2": {
        "api": "dnakit.alignment.map_long_reads",
        "functions": ["dnakit.alignment.map_long_reads"],
        "sources": ["https://github.com/lh3/minimap2"],
        "citation_dois": ["10.1093/bioinformatics/bty191"],
        "citation_references": ["Li (2018), DOI 10.1093/bioinformatics/bty191"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed minimap2 with a bounded preset allow-list and SAM/PAF "
            "output; the full minimap2 command surface and index-management modes are not "
            "re-exported."
        ),
    },
    "bedtools": {
        "api": "dnakit.alignment.intersect_bed",
        "functions": ["dnakit.alignment.intersect_bed"],
        "sources": ["https://github.com/arq5x/bedtools2"],
        "citation_dois": ["10.1093/bioinformatics/btq033"],
        "citation_references": ["Quinlan & Hall (2010), DOI 10.1093/bioinformatics/btq033"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed BEDTools for the fixed local BED intersection operation; "
            "the complete BEDTools subcommand and option surface is not re-exported."
        ),
    },
    "bedops": {
        "api": "dnakit.alignment.sort_bed_bedops",
        "functions": ["dnakit.alignment.sort_bed_bedops"],
        "sources": ["https://github.com/bedops/bedops"],
        "citation_dois": ["10.1093/bioinformatics/bts277"],
        "citation_references": [
            "Neph et al. (2012), DOI 10.1093/bioinformatics/bts277"
        ],
        "license_evidence": ["GPL-2.0"],
        "boundary": (
            "Runs caller-installed BEDOPS sort-bed for a bounded BED sort operation; "
            "the GPL-2.0 code, other BEDOPS commands and format converters remain separate."
        ),
    },
    "bwa": {
        "api": "dnakit.ngs.map_reads",
        "functions": ["dnakit.ngs.map_reads"],
        "sources": ["https://github.com/lh3/bwa"],
        "citation_dois": ["10.1093/bioinformatics/btp324"],
        "citation_references": ["Li & Durbin (2009), DOI 10.1093/bioinformatics/btp324"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed BWA executable through the bounded map_reads facade; "
            "the facade does not bundle BWA or expose every BWA command."
        ),
    },
    "bwa-mem2": {
        "api": "dnakit.ngs.map_reads",
        "functions": ["dnakit.ngs.map_reads"],
        "sources": ["https://github.com/bwa-mem2/bwa-mem2"],
        "citation_dois": ["10.1109/IPDPS.2019.00041"],
        "citation_references": [
            "Vasimuddin et al. (2019), DOI 10.1109/IPDPS.2019.00041"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed BWA-MEM2 executable through the bounded map_reads facade; "
            "the facade does not bundle BWA-MEM2 or expose every BWA-MEM2 command."
        ),
    },
    "iq-tree2": {
        "api": "dnakit.phylogeny.infer_ml_tree",
        "functions": ["dnakit.phylogeny.infer_ml_tree"],
        "sources": ["https://github.com/iqtree/iqtree2"],
        "citation_dois": ["10.1093/molbev/msaa015"],
        "citation_references": [
            "Minh et al. (2020), DOI 10.1093/molbev/msaa015"
        ],
        "license_evidence": ["GPL-2.0"],
        "boundary": (
            "Runs caller-installed IQ-TREE 2 for bounded maximum-likelihood tree inference; "
            "the full IQ-TREE command surface and GPL code are not bundled into DNAKit."
        ),
    },
    "fasttree": {
        "api": "dnakit.phylogeny.infer_fast_tree",
        "functions": ["dnakit.phylogeny.infer_fast_tree"],
        "sources": ["https://github.com/morgannprice/fasttree"],
        "citation_dois": ["10.1371/journal.pone.0009490"],
        "citation_references": [
            "Price, Dehal & Arkin (2010), DOI 10.1371/journal.pone.0009490"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed FastTree executable for bounded approximate maximum-"
            "likelihood inference; the complete FastTree CLI and GPL code are not bundled."
        ),
    },
    "shasta": {
        "api": "dnakit.assembly.assemble_shasta",
        "functions": ["dnakit.assembly.assemble_shasta"],
        "sources": ["https://github.com/paoloshasta/shasta"],
        "citation_dois": ["10.1038/s41587-020-0503-6"],
        "citation_references": [
            "Shafin et al. (2020), DOI 10.1038/s41587-020-0503-6"
        ],
        "license_evidence": ["MIT (current Shasta source LICENSE)"],
        "boundary": (
            "Runs a caller-installed Shasta long-read assembler with bounded input, "
            "configuration, thread and output-directory controls; the complete Shasta "
            "CLI, bundled dependencies and source are not copied into DNAKit."
        ),
    },
    "merqury": {
        "api": "dnakit.assembly.evaluate_merqury",
        "functions": ["dnakit.assembly.evaluate_merqury"],
        "sources": ["https://github.com/marbl/merqury"],
        "citation_dois": ["10.1186/s13059-020-02134-9"],
        "citation_references": [
            "Rhie et al. (2020), DOI 10.1186/s13059-020-02134-9"
        ],
        "license_evidence": ["Public domain / United States Government Work"],
        "boundary": (
            "Runs a caller-installed Merqury evaluation with an existing Meryl read database "
            "and assembly FASTA; only the documented QV/completeness report boundary is exposed, "
            "and Meryl databases and all upstream scripts remain caller responsibilities."
        ),
    },
    "busco": {
        "api": "dnakit.assembly.run_busco",
        "functions": ["dnakit.assembly.run_busco"],
        "sources": ["https://gitlab.com/ezlab/busco"],
        "citation_dois": ["10.1093/molbev/msab199"],
        "citation_references": [
            "Manni et al. (2021), DOI 10.1093/molbev/msab199"
        ],
        "license_evidence": ["MIT (BUSCO software; lineage datasets have separate terms)"],
        "boundary": (
            "Runs a caller-installed BUSCO completeness assessment with bounded input, lineage, "
            "mode, output-name and thread controls; BUSCO lineage datasets, databases and the "
            "broader workflow remain separate dependencies."
        ),
    },
    "flye": {
        "api": "dnakit.assembly.assemble_flye",
        "functions": ["dnakit.assembly.assemble_flye"],
        "sources": ["https://github.com/mikolmogorov/Flye"],
        "citation_dois": ["10.1038/s41587-019-0072-8"],
        "citation_references": [
            "Kolmogorov et al. (2019), DOI 10.1038/s41587-019-0072-8"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "Runs caller-installed Flye for bounded long-read assembly with fixed platform, "
            "genome-size, thread and output controls; metagenome and polishing subcommands "
            "remain separate upstream responsibilities."
        ),
    },
    "hifiasm": {
        "api": "dnakit.assembly.assemble_hifiasm",
        "functions": ["dnakit.assembly.assemble_hifiasm"],
        "sources": ["https://github.com/chhylp123/hifiasm"],
        "citation_dois": ["10.1038/s41592-020-01056-5"],
        "citation_references": [
            "Cheng et al. (2021), DOI 10.1038/s41592-020-01056-5"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed hifiasm for bounded HiFi assembly graph generation; "
            "the full haplotype, trio and polyploid option surface is not re-exported."
        ),
    },
    "megahit": {
        "api": "dnakit.assembly.assemble_megahit",
        "functions": ["dnakit.assembly.assemble_megahit"],
        "sources": ["https://github.com/voutcn/megahit"],
        "citation_dois": ["10.1093/bioinformatics/btv033"],
        "citation_references": [
            "Li et al. (2015), DOI 10.1093/bioinformatics/btv033"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed MEGAHIT for bounded single-end or paired-end assembly; "
            "GPL-3.0 code and the complete memory/k-mer option surface remain upstream."
        ),
    },
    "freebayes": {
        "api": "dnakit.variants.call_freebayes",
        "functions": ["dnakit.variants.call_freebayes"],
        "sources": ["https://github.com/freebayes/freebayes"],
        "citation_references": [
            "Garrison & Marth (2012), arXiv:1207.3907"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed FreeBayes against caller-prepared BAM/CRAM inputs with "
            "bounded ploidy and alternate-support filters; arbitrary FreeBayes options are not "
            "exposed."
        ),
    },
    "gemini": {
        "api": "dnakit.variants.load_gemini_database / dnakit.variants.query_gemini_database",
        "functions": [
            "dnakit.variants.load_gemini_database",
            "dnakit.variants.query_gemini_database",
        ],
        "sources": [
            "https://github.com/arq5x/gemini",
            "https://gemini.readthedocs.io/en/latest/content/quick_start.html",
        ],
        "citation_dois": ["10.1371/journal.pcbi.1003153"],
        "citation_references": [
            "Paila et al. (2013), DOI 10.1371/journal.pcbi.1003153"
        ],
        "license_evidence": ["MIT (upstream repository)"],
        "boundary": (
            "Loads caller-prepared VCF files into a local GEMINI database and executes a "
            "bounded SQL query with tabular output; GEMINI annotation assets, database "
            "construction inputs and the complete command surface remain external."
        ),
    },
    "ensembl vep": {
        "api": "dnakit.variants.annotate_vcf_ensembl_vep",
        "functions": ["dnakit.variants.annotate_vcf_ensembl_vep"],
        "sources": ["https://github.com/Ensembl/ensembl-vep"],
        "citation_dois": ["10.1186/s13059-016-0974-4"],
        "citation_references": [
            "McLaren et al. (2016), DOI 10.1186/s13059-016-0974-4"
        ],
        "license_evidence": ["Apache-2.0 (upstream repository)"],
        "boundary": (
            "Runs caller-installed Ensembl VEP on a caller-prepared VCF with bounded species, "
            "assembly, cache, FASTA, output-format and fork options; VEP caches, transcript "
            "resources, plugins and the complete CLI remain external."
        ),
    },
    "gatk4": {
        "api": "dnakit.variants.call_gatk_haplotypecaller",
        "functions": ["dnakit.variants.call_gatk_haplotypecaller"],
        "sources": ["https://github.com/broadinstitute/gatk"],
        "citation_dois": ["10.1101/gr.107524.110"],
        "citation_references": [
            "McKenna et al. (2010), DOI 10.1101/gr.107524.110"
        ],
        "license_evidence": ["Apache-2.0 (current GATK source LICENSE.TXT)"],
        "boundary": (
            "Runs a caller-installed GATK4 HaplotypeCaller with bounded reference, BAM/CRAM, "
            "VCF/GVCF and confidence controls; the broader GATK tool collection, Java runtime "
            "and reference resources are not bundled."
        ),
    },
    "deepvariant": {
        "api": "dnakit.variants.call_deepvariant",
        "functions": ["dnakit.variants.call_deepvariant"],
        "sources": ["https://github.com/google/deepvariant"],
        "citation_dois": ["10.1038/nbt.4235"],
        "citation_references": [
            "Poplin et al. (2018), DOI 10.1038/nbt.4235"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "Runs a caller-installed DeepVariant variant-calling entry point with bounded "
            "reference, BAM/CRAM, model and VCF/GVCF controls; models, containers and the "
            "broader pipeline are not bundled into DNAKit."
        ),
    },
    "cutesv": {
        "api": "dnakit.variants.call_cutesv",
        "functions": ["dnakit.variants.call_cutesv"],
        "sources": ["https://github.com/tjiangHIT/cuteSV"],
        "citation_dois": ["10.1186/s13059-020-02107-y"],
        "citation_references": [
            "Jiang et al. (2020), DOI 10.1186/s13059-020-02107-y"
        ],
        "license_evidence": ["MIT (upstream repository)"],
        "boundary": (
            "Runs caller-installed cuteSV with fixed long-read alignment, reference, VCF, "
            "work-directory, thread and sample controls; the complete cuteSV option surface "
            "and Python implementation are not copied into DNAKit."
        ),
    },
    "sniffles2": {
        "api": "dnakit.variants.call_sniffles2",
        "functions": ["dnakit.variants.call_sniffles2"],
        "sources": ["https://github.com/fritzsedlazeck/Sniffles"],
        "citation_dois": ["10.1038/s41587-023-02024-y"],
        "citation_references": [
            "Smolka et al. (2024), DOI 10.1038/s41587-023-02024-y"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed Sniffles2 with bounded BAM/CRAM, VCF, reference, thread, "
            "sample and minimum-SV-length controls; population, mosaic and auxiliary modes "
            "remain upstream responsibilities."
        ),
    },
    "delly": {
        "api": "dnakit.variants.call_delly",
        "functions": ["dnakit.variants.call_delly"],
        "sources": ["https://github.com/dellytools/delly"],
        "citation_dois": ["10.1093/bioinformatics/bts378"],
        "citation_references": [
            "Rausch et al. (2012), DOI 10.1093/bioinformatics/bts378"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "Runs caller-installed DELLY short-read or long-read SV discovery with a bounded "
            "reference, alignment and mode boundary; genotyping, CNV, filtering and the full "
            "DELLY command surface remain upstream."
        ),
    },
    "vcftools": {
        "api": "dnakit.variants.filter_vcf_vcftools",
        "functions": ["dnakit.variants.filter_vcf_vcftools"],
        "sources": ["https://github.com/vcftools/vcftools"],
        "citation_dois": ["10.1093/bioinformatics/btr330"],
        "citation_references": [
            "Danecek et al. (2011), DOI 10.1093/bioinformatics/btr330"
        ],
        "license_evidence": ["LGPL-3.0"],
        "boundary": (
            "Runs caller-installed VCFtools for a fixed recoding and filtering subset; "
            "the full VCFtools option surface is not re-exported."
        ),
    },
    "vcfanno": {
        "api": "dnakit.variants.annotate_vcf_vcfanno",
        "functions": ["dnakit.variants.annotate_vcf_vcfanno"],
        "sources": ["https://github.com/brentp/vcfanno"],
        "citation_dois": ["10.1186/s13059-016-0973-5"],
        "citation_references": [
            "Pedersen et al. (2016), DOI 10.1186/s13059-016-0973-5"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed vcfanno with a caller-provided TOML configuration and "
            "fixed VCF output; the full vcfanno Lua and annotation option surface is not "
            "re-exported as Python parameters."
        ),
    },
    "spoa": {
        "api": "dnakit.alignment.consensus_spoa",
        "functions": ["dnakit.alignment.consensus_spoa"],
        "sources": ["https://github.com/rvaser/spoa"],
        "citation_dois": ["10.1093/bioinformatics/btg109"],
        "citation_references": [
            "Lee (2003), DOI 10.1093/bioinformatics/btg109"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed SPOA for a bounded consensus-FASTA operation; its full "
            "alignment, graph and scoring option surface is not re-exported."
        ),
    },
    "survivor": {
        "api": "dnakit.variants.merge_structural_variants_survivor",
        "functions": ["dnakit.variants.merge_structural_variants_survivor"],
        "sources": ["https://github.com/fritzsedlazeck/SURVIVOR"],
        "citation_dois": ["10.1038/ncomms14061"],
        "citation_references": [
            "Jeffares et al. (2017), DOI 10.1038/NCOMMS14061"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed SURVIVOR's bounded VCF merge mode; simulation, conversion, "
            "filtering and comparison subcommands remain outside this facade."
        ),
    },
    "picard": {
        "api": "dnakit.ngs.mark_duplicates_picard",
        "functions": ["dnakit.ngs.mark_duplicates_picard"],
        "sources": ["https://github.com/broadinstitute/picard"],
        "citation_references": ["Picard Toolkit (2019), Broad Institute"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed Picard MarkDuplicates with fixed input, output, metrics and "
            "duplicate-policy parameters; the complete Picard command surface is not bundled."
        ),
    },
    "bowtie2": {
        "api": "dnakit.ngs.map_reads",
        "functions": ["dnakit.ngs.map_reads"],
        "sources": ["https://github.com/BenLangmead/bowtie2"],
        "citation_dois": ["10.1038/nmeth.1923"],
        "citation_references": ["Langmead & Salzberg (2012), DOI 10.1038/nmeth.1923"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed Bowtie2 and Bowtie2-build executables through the bounded "
            "map_reads facade; the complete Bowtie2 CLI is not re-exported."
        ),
    },
    "samtools": {
        "api": "dnakit.ngs.prepare_bam",
        "functions": ["dnakit.ngs.prepare_bam"],
        "sources": ["https://github.com/samtools/samtools"],
        "citation_dois": ["10.1093/gigascience/giab008"],
        "citation_references": ["Danecek et al. (2021), DOI 10.1093/gigascience/giab008"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed SAMtools executable through bounded sort, duplicate "
            "marking, indexing and summary stages; arbitrary SAMtools subcommands are not exposed."
        ),
    },
    "bcftools": {
        "api": "dnakit.ngs.call_variants",
        "functions": ["dnakit.ngs.call_variants"],
        "sources": ["https://github.com/samtools/bcftools"],
        "citation_dois": ["10.1093/gigascience/giab008"],
        "citation_references": ["Danecek et al. (2021), DOI 10.1093/gigascience/giab008"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed BCFtools executable through fixed mpileup, call, norm "
            "and index stages; arbitrary BCFtools plugins and commands are not exposed."
        ),
    },
    "htslib": {
        "api": "dnakit.ngs.compress_index_vcf_htslib",
        "functions": ["dnakit.ngs.compress_index_vcf_htslib"],
        "sources": ["https://github.com/samtools/htslib"],
        "citation_dois": ["10.1093/gigascience/giab007"],
        "citation_references": [
            "Bonfield et al. (2021), DOI 10.1093/gigascience/giab007"
        ],
        "license_evidence": ["MIT/BSD (HTSlib distribution)"],
        "boundary": (
            "Runs caller-installed HTSlib bgzip and tabix utilities for bounded VCF compression "
            "and indexing; the C library API and broader samtools/bcftools surfaces remain separate."
        ),
    },
    "seqtk": {
        "api": "dnakit.ngs.sample_sequences_seqtk",
        "functions": ["dnakit.ngs.sample_sequences_seqtk"],
        "sources": ["https://github.com/lh3/seqtk"],
        "citation_references": [
            "seqtk project authors and contributors; the project does not provide a DOI"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed seqtk for bounded FASTA/FASTQ random sampling; conversion, "
            "trimming, subsequence and telomere subcommands remain upstream responsibilities."
        ),
    },
    "whatshap": {
        "api": "dnakit.ngs.phase_variants_whatshap",
        "functions": ["dnakit.ngs.phase_variants_whatshap"],
        "sources": ["https://github.com/whatshap/whatshap"],
        "citation_dois": ["10.1101/085050", "10.1089/cmb.2014.0157"],
        "citation_references": [
            "Martin et al. (2016), DOI 10.1101/085050; Patterson et al. (2015), DOI 10.1089/cmb.2014.0157"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed WhatsHap phase with bounded VCF, BAM/CRAM, reference and "
            "output controls; pedigree, polyploid, genotype and auxiliary modes remain upstream."
        ),
    },
    "bamtools": {
        "api": "dnakit.ngs.sort_bam_bamtools",
        "functions": ["dnakit.ngs.sort_bam_bamtools"],
        "sources": ["https://github.com/pezmaster31/bamtools"],
        "citation_dois": ["10.1093/bioinformatics/btr174"],
        "citation_references": [
            "Barnett et al. (2011), DOI 10.1093/bioinformatics/btr174"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed BamTools for bounded BAM sorting; filtering, merging, "
            "random access and the C++ API remain outside this facade."
        ),
    },
    "vsearch": {
        "api": "dnakit.fasta_fastq.cluster_vsearch",
        "functions": ["dnakit.fasta_fastq.cluster_vsearch"],
        "sources": ["https://github.com/torognes/vsearch"],
        "citation_dois": ["10.7717/peerj.2584"],
        "citation_references": [
            "Rognes et al. (2016), DOI 10.7717/peerj.2584"
        ],
        "license_evidence": ["BSD-2-Clause or GPL-3.0 (dual license)"],
        "boundary": (
            "Runs caller-installed VSEARCH greedy FASTA clustering with bounded identity, "
            "thread and centroid/UC output controls; search, dereplication, chimera and the "
            "complete option surface remain upstream."
        ),
    },
    "spades": {
        "api": "dnakit.ngs.assemble_reads",
        "functions": ["dnakit.ngs.assemble_reads"],
        "sources": ["https://github.com/ablab/spades"],
        "citation_dois": ["10.1089/cmb.2012.0021"],
        "citation_references": ["Bankevich et al. (2012), DOI 10.1089/cmb.2012.0021"],
        "license_evidence": ["GPL-2.0"],
        "boundary": (
            "Runs a caller-installed SPAdes executable through a bounded short-read or hybrid "
            "assembly facade; the complete SPAdes command suite is not bundled."
        ),
    },
    "diamond": {
        "api": "dnakit.alignment.search_diamond_blastp / dnakit.alignment.search_diamond_blastx",
        "functions": [
            "dnakit.alignment.search_diamond_blastp",
            "dnakit.alignment.search_diamond_blastx",
        ],
        "sources": [
            "https://github.com/bbuchfink/diamond",
            "https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options",
        ],
        "citation_dois": ["10.1038/s41592-021-01101-x"],
        "citation_references": ["Buchfink et al. (2021), DOI 10.1038/s41592-021-01101-x"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Builds a caller-provided protein database and exposes fixed DIAMOND blastp or "
            "blastx tabular searches; the full CLI, database distribution and functional "
            "interpretation remain upstream responsibilities."
        ),
    },
    "mmseqs2": {
        "api": "dnakit.alignment.search_mmseqs_protein / dnakit.alignment.cluster_mmseqs_protein",
        "functions": [
            "dnakit.alignment.search_mmseqs_protein",
            "dnakit.alignment.cluster_mmseqs_protein",
        ],
        "sources": ["https://github.com/soedinglab/MMseqs2"],
        "citation_dois": ["10.1038/nbt.3988"],
        "citation_references": ["Steinegger & Söding (2017), DOI 10.1038/nbt.3988"],
        "license_evidence": ["MIT"],
        "boundary": (
            "Exposes fixed MMseqs2 easy-search and easy-cluster protein FASTA workflows with "
            "explicit memory and identity/coverage controls; temporary databases, the full "
            "command surface and biological interpretation remain upstream."
        ),
    },
    "kraken2": {
        "api": "dnakit.microbiome.classify_kraken2",
        "functions": ["dnakit.microbiome.classify_kraken2"],
        "sources": ["https://github.com/DerrickWood/kraken2/wiki/Manual"],
        "citation_dois": [
            "10.1186/s13059-019-1891-0",
            "10.1186/s13059-018-1568-0",
        ],
        "citation_references": [
            "Wood et al. (2019), DOI 10.1186/s13059-019-1891-0",
            "Wood et al. (2018), DOI 10.1186/s13059-018-1568-0",
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed Kraken2 against a caller-prepared database and returns the "
            "classification and report files; database construction, abundance re-estimation "
            "and taxonomic interpretation are not bundled."
        ),
    },
    "prokka": {
        "api": "dnakit.annotation.annotate_prokka",
        "functions": ["dnakit.annotation.annotate_prokka"],
        "sources": ["https://github.com/tseemann/prokka"],
        "citation_dois": ["10.1093/bioinformatics/btu153"],
        "citation_references": ["Seemann (2014), DOI 10.1093/bioinformatics/btu153"],
        "license_evidence": ["GPL-2.0"],
        "boundary": (
            "Runs caller-installed Prokka with bounded output prefix and CPU controls and "
            "requires its separately prepared databases; the complete annotation option "
            "surface is not re-exported."
        ),
    },
    "bakta": {
        "api": "dnakit.annotation.annotate_bakta",
        "functions": ["dnakit.annotation.annotate_bakta"],
        "sources": ["https://github.com/oschwengers/bakta"],
        "citation_dois": ["10.1099/mgen.0.000685"],
        "citation_references": ["Schwengers et al. (2021), DOI 10.1099/mgen.0.000685"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed Bakta with an explicitly supplied database, translation "
            "table and minimum contig length; database preparation, AMRFinder resources and "
            "the complete annotation surface remain caller responsibilities."
        ),
    },
    "checkm2": {
        "api": "dnakit.annotation.predict_checkm2",
        "functions": ["dnakit.annotation.predict_checkm2"],
        "sources": ["https://github.com/chklovski/CheckM2"],
        "citation_dois": ["10.1038/s41592-023-01940-w"],
        "citation_references": ["Chklovski et al. (2023), DOI 10.1038/s41592-023-01940-w"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed CheckM2 against a bins directory and compatible DIAMOND "
            "database; model files, database preparation and scientific interpretation of "
            "completeness/contamination estimates remain external."
        ),
    },
    "bracken": {
        "api": "dnakit.microbiome.estimate_bracken_abundance",
        "functions": ["dnakit.microbiome.estimate_bracken_abundance"],
        "sources": ["https://ccb.jhu.edu/software/bracken/"],
        "citation_dois": ["10.7717/peerj-cs.104", "10.1038/s41596-022-00738-y"],
        "citation_references": [
            "Lu et al. (2017), DOI 10.7717/peerj-cs.104",
            "Lu et al. (2022), DOI 10.1038/s41596-022-00738-y",
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed Bracken against a Kraken2 report and prepared database "
            "with bounded read-length, rank and threshold controls; database construction and "
            "taxonomic interpretation remain external."
        ),
    },
    "centrifuge": {
        "api": "dnakit.microbiome.classify_centrifuge",
        "functions": ["dnakit.microbiome.classify_centrifuge"],
        "sources": ["https://github.com/infphilo/centrifuge"],
        "citation_dois": ["10.1101/gr.210641.116"],
        "citation_references": ["Kim et al. (2016), DOI 10.1101/gr.210641.116"],
        "license_evidence": ["GPL-3.0 (Centrifuge source)"],
        "boundary": (
            "Runs caller-installed Centrifuge against a prebuilt index and writes fixed read "
            "assignments and a report; index construction and taxonomic interpretation are "
            "not bundled."
        ),
    },
    "kaiju": {
        "api": "dnakit.microbiome.classify_kaiju",
        "functions": ["dnakit.microbiome.classify_kaiju"],
        "sources": ["https://github.com/bioinformatics-centre/kaiju"],
        "citation_dois": ["10.1038/ncomms11257"],
        "citation_references": ["Menzel et al. (2016), DOI 10.1038/ncomms11257"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed Kaiju with explicit protein-database and NCBI taxonomy "
            "files; database construction, preprocessing and downstream abundance analysis "
            "remain caller responsibilities."
        ),
    },
    "unicycler": {
        "api": "dnakit.assembly.assemble_unicycler",
        "functions": ["dnakit.assembly.assemble_unicycler"],
        "sources": ["https://github.com/rrwick/Unicycler"],
        "citation_dois": ["10.1371/journal.pcbi.1005595"],
        "citation_references": ["Wick et al. (2017), DOI 10.1371/journal.pcbi.1005595"],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed Unicycler for bounded short-read, long-read or hybrid "
            "bacterial assembly; the underlying SPAdes, Pilon, minimap and complete pipeline "
            "options remain external."
        ),
    },
    "gtdb-tk": {
        "api": "dnakit.microbiome.classify_gtdbtk",
        "functions": ["dnakit.microbiome.classify_gtdbtk"],
        "sources": ["https://github.com/Ecogenomics/GTDBTk"],
        "citation_dois": ["10.1093/bioinformatics/btz848", "10.1093/bioinformatics/btac672"],
        "citation_references": [
            "Chaumeil et al. (2020), DOI 10.1093/bioinformatics/btz848",
            "Chaumeil et al. (2022), DOI 10.1093/bioinformatics/btac672",
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs caller-installed GTDB-Tk classify_wf with GTDBTK_DATA_PATH set to the "
            "caller-supplied reference-data directory, plus genome directory, extension and "
            "CPU limit; the GTDB reference package and full workflow remain external."
        ),
    },
    "humann": {
        "api": "dnakit.microbiome.profile_humann",
        "functions": ["dnakit.microbiome.profile_humann"],
        "sources": ["https://github.com/biobakery/humann"],
        "citation_dois": ["10.7554/eLife.65088"],
        "citation_references": ["Beghini et al. (2021), DOI 10.7554/eLife.65088"],
        "license_evidence": ["GPL-3.0 (upstream repository)"],
        "boundary": (
            "Runs caller-installed HUMAnN with fixed input/output and optional nucleotide or "
            "protein database paths; database preparation, normalization choices and the full "
            "biobakery workflow remain external."
        ),
    },
    "methyldackel": {
        "api": "dnakit.epigenome.extract_methylation_methyldackel",
        "functions": ["dnakit.epigenome.extract_methylation_methyldackel"],
        "sources": ["https://github.com/dpryan79/MethylDackel"],
        "citation_references": [
            "MethylDackel project authors; the official repository requests URL citation and supplies no DOI"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs caller-installed MethylDackel extract against an indexed reference and "
            "coordinate-sorted alignment, exposing CpG/CHG/CHH bedGraph outputs; indexing, "
            "alignment quality and downstream methylation interpretation remain external."
        ),
    },
    "cutadapt": {
        "api": "dnakit.sequencing_qc.run_cutadapt",
        "functions": ["dnakit.sequencing_qc.run_cutadapt"],
        "boundary": (
            "Runs a caller-installed Cutadapt executable through a bounded argv adapter; "
            "the adapter does not install Cutadapt or provide its complete CLI."
        ),
        "sources": [
            "https://github.com/marcelm/cutadapt",
            "https://cutadapt.readthedocs.io/en/stable/guide.html",
        ],
        "citation_dois": ["10.14806/ej.17.1.200"],
        "citation_references": ["Martin (2011), DOI 10.14806/ej.17.1.200"],
        "license_evidence": ["MIT"],
    },
    "fastp": {
        "api": "dnakit.sequencing_qc.run_fastp",
        "functions": ["dnakit.sequencing_qc.run_fastp"],
        "boundary": (
            "Runs a caller-installed fastp executable through a bounded argv adapter; "
            "the adapter does not install fastp or provide its complete CLI."
        ),
        "sources": ["https://github.com/OpenGene/fastp"],
        "citation_dois": ["10.1002/imt2.70078"],
        "citation_references": ["Chen (2025), DOI 10.1002/imt2.70078"],
        "license_evidence": ["MIT"],
    },
    "fastqc": {
        "api": "dnakit.sequencing_qc.run_fastqc",
        "functions": ["dnakit.sequencing_qc.run_fastqc"],
        "boundary": (
            "Runs a caller-installed FastQC executable through a bounded argv adapter; "
            "the adapter does not bundle Java or reimplement the complete report engine."
        ),
        "sources": ["https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"],
        "citation_references": ["FastQC documentation; no DOI supplied by the survey"],
        "license_evidence": ["GPL-3.0"],
    },
    "multiqc": {
        "api": "dnakit.sequencing_qc.run_multiqc",
        "functions": ["dnakit.sequencing_qc.run_multiqc"],
        "boundary": (
            "Runs a caller-installed MultiQC executable through a bounded argv adapter; "
            "the adapter does not bundle the GPL backend or its plugin collection."
        ),
        "sources": ["https://github.com/MultiQC/MultiQC"],
        "citation_dois": ["10.1093/bioinformatics/btw354"],
        "citation_references": ["Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354"],
        "license_evidence": ["GPL-3.0-or-later"],
    },
    "nanofilt": {
        "api": "dnakit.sequencing_qc.filter_nanopore_reads",
        "functions": ["dnakit.sequencing_qc.filter_nanopore_reads"],
        "sources": ["https://github.com/wdecoster/nanofilt"],
        "citation_dois": ["10.1093/bioinformatics/bty149"],
        "citation_references": [
            "De Coster et al. (2018), DOI 10.1093/bioinformatics/bty149"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed NanoFilt executable for bounded quality, length and "
            "end-cropping filters; NanoFilt is no longer actively updated upstream and the "
            "complete CLI is not re-exported."
        ),
    },
    "filtlong": {
        "api": "dnakit.sequencing_qc.filter_long_reads",
        "functions": ["dnakit.sequencing_qc.filter_long_reads"],
        "sources": ["https://github.com/rrwick/Filtlong"],
        "citation_references": [
            "Filtlong project authors; the official project supplies no DOI"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed Filtlong executable with bounded numeric length and "
            "retention thresholds; reference-driven, unit-suffixed and trimming modes are "
            "not re-exported by this facade."
        ),
    },
    "nanoplot": {
        "api": "dnakit.sequencing_qc.summarize_long_reads",
        "functions": ["dnakit.sequencing_qc.summarize_long_reads"],
        "sources": ["https://github.com/wdecoster/NanoPlot"],
        "citation_dois": ["10.1093/bioinformatics/btad311"],
        "citation_references": [
            "De Coster & Rademakers (2023), DOI 10.1093/bioinformatics/btad311"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed NanoPlot executable in bounded FASTQ summary mode with "
            "static plots disabled; the complete plotting and input-format surface is not "
            "re-exported."
        ),
    },
    "porechop": {
        "api": "dnakit.sequencing_qc.trim_porechop",
        "functions": ["dnakit.sequencing_qc.trim_porechop"],
        "sources": ["https://github.com/rrwick/Porechop"],
        "citation_dois": ["10.1099/mgen.0.000132"],
        "citation_references": [
            "Wick et al. (2017), DOI 10.1099/mgen.0.000132"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed Porechop executable for bounded adapter trimming and "
            "optional no-split processing; the upstream project is officially unsupported and "
            "the complete CLI is not re-exported."
        ),
    },
    "trimmomatic": {
        "api": "dnakit.sequencing_qc.trim_trimmomatic",
        "functions": ["dnakit.sequencing_qc.trim_trimmomatic"],
        "sources": ["https://github.com/usadellab/Trimmomatic"],
        "citation_dois": ["10.1093/bioinformatics/btu170"],
        "citation_references": [
            "Bolger, Lohse & Usadel (2014), DOI 10.1093/bioinformatics/btu170"
        ],
        "license_evidence": ["GPL-3.0 (adapter sequences have separate Illumina terms)"],
        "boundary": (
            "Runs a caller-installed Trimmomatic executable through fixed single-end or "
            "paired-end quality and adapter steps; adapter sequence terms and the complete "
            "Trimmomatic CLI remain upstream responsibilities."
        ),
    },
    "afterqc": {
        "api": "dnakit.sequencing_qc.run_afterqc",
        "functions": ["dnakit.sequencing_qc.run_afterqc"],
        "sources": ["https://github.com/OpenGene/AfterQC"],
        "citation_dois": ["10.1186/s12859-017-1469-3"],
        "citation_references": [
            "Chen et al. (2017), DOI 10.1186/s12859-017-1469-3"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "Runs a caller-installed AfterQC script for bounded single-end or paired-end "
            "filtering and QC output; the complete option surface and Python source are not "
            "bundled."
        ),
    },
    "biopython": {
        "api": "dnakit.bioinformatics",
        "functions": [
            "dnakit.bioinformatics.reads.read_fastx",
            "dnakit.bioinformatics.alignment.read_alignment",
            "dnakit.bioinformatics.motifs.build_motif",
        ],
        "boundary": (
            "Related DNAKit entry points; this is not a claim that all Biopython APIs are "
            "re-exported."
        ),
    },
    "biotite": {
        "api": "dnakit.integrations.align_biotite",
        "functions": ["dnakit.integrations.align_biotite"],
        "citation_dois": ["10.1186/s12859-018-2367-z"],
        "citation_references": ["Kunzmann & Hamacher (2018), DOI 10.1186/s12859-018-2367-z"],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "The optional BSD-3-Clause backend exposes bounded optimal nucleotide pairwise "
            "alignment; Biotite's broader structure and database APIs are not re-exported."
        ),
    },
    "edlib": {
        "api": "dnakit.integrations.align_edlib",
        "functions": ["dnakit.integrations.align_edlib"],
        "citation_dois": ["10.1093/bioinformatics/btw753"],
        "citation_references": ["Šošić & Šikić (2017), DOI 10.1093/bioinformatics/btw753"],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded edit-distance alignment; arbitrary "
            "native extension calls are not available through the JSON boundary."
        ),
    },
    "seqkit": {
        "api": "dnakit.fasta_fastq",
        "functions": [
            "dnakit.fasta_fastq.calculate_gc_profile",
            "dnakit.fasta_fastq.calculate_length_distribution",
            "dnakit.fasta_fastq.convert_to_fasta",
            "dnakit.fasta_fastq.count_records",
            "dnakit.fasta_fastq.count_reads",
        ],
        "sources": [
            "https://doi.org/10.1002/imt2.191",
            "https://github.com/shenwei356/seqkit",
        ],
        "citation_dois": [
            "10.1002/imt2.191",
            "10.1371/journal.pone.0163962",
        ],
        "citation_references": [
            "Shen, Sipos & Zhao (2024), SeqKit2, DOI 10.1002/imt2.191",
            "Shen et al. (2016), SeqKit, DOI 10.1371/journal.pone.0163962",
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "SeqKit remains a separately installed CLI for operations not covered by these "
            "file-oriented functions."
        ),
    },
    "seqmagick": {
        "api": "dnakit.fasta_fastq.convert_seqmagick",
        "functions": ["dnakit.fasta_fastq.convert_seqmagick"],
        "sources": ["https://github.com/fhcrc/seqmagick"],
        "citation_references": [
            "Matsen Group seqmagick project; the official project supplies no DOI"
        ],
        "license_evidence": ["GPL-3.0"],
        "boundary": (
            "Runs a caller-installed SeqMagick executable for bounded format conversion and "
            "common sequence transformations; GPL-3.0 code is not bundled and the complete "
            "SeqMagick CLI is not re-exported."
        ),
    },
    "scikit-bio": {
        "api": (
            "dnakit.integrations.beta_diversity_scikit_bio / "
            "dnakit.integrations.permanova_scikit_bio"
        ),
        "functions": [
            "dnakit.integrations.beta_diversity_scikit_bio",
            "dnakit.integrations.permanova_scikit_bio",
        ],
        "boundary": (
            "The adapter exposes two bounded statistics operations; scikit-bio's full "
            "namespace is not copied."
        ),
    },
    "pydna": {
        "api": "dnakit.bioinformatics.cloning",
        "functions": [
            "dnakit.bioinformatics.cloning.assemble_duplexes",
            "dnakit.bioinformatics.cloning.read_constructs",
        ],
        "boundary": (
            "DNAKit adapters cover the documented construct operations; pydna's complete "
            "API is not implied."
        ),
    },
    "pysam": {
        "api": "dnakit.bioinformatics.reads",
        "functions": [
            "dnakit.bioinformatics.reads.index_reference",
            "dnakit.bioinformatics.reads.alignment_summary",
            "dnakit.bioinformatics.reads.coverage_consensus",
        ],
        "boundary": (
            "DNAKit exposes focused FASTA/BAM indexing and summary operations through pysam; "
            "pysam's complete API is not re-exported."
        ),
    },
    "pyranges": {
        "api": (
            "dnakit.integrations.summarize_pyranges_bed / "
            "dnakit.integrations.merge_pyranges_bed"
        ),
        "functions": [
            "dnakit.integrations.summarize_pyranges_bed",
            "dnakit.integrations.merge_pyranges_bed",
        ],
        "citation_dois": ["10.1093/bioinformatics/btz615"],
        "citation_references": ["Stovner & Sætrom (2020), DOI 10.1093/bioinformatics/btz615"],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded BED interval summaries and merges; "
            "the full PyRanges/pandas API is not re-exported."
        ),
    },
    "pybigwig": {
        "api": "dnakit.integrations.summarize_pybigwig_region",
        "functions": ["dnakit.integrations.summarize_pybigwig_region"],
        "citation_dois": ["10.5281/zenodo.45238"],
        "citation_references": ["pyBigWig software release, DOI 10.5281/zenodo.45238"],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded, read-only local bigWig region "
            "statistics; file creation, bigBed entries and remote URLs are not re-exported."
        ),
    },
    "pybedtools": {
        "api": "dnakit.integrations.intersect_pybedtools_bed",
        "functions": ["dnakit.integrations.intersect_pybedtools_bed"],
        "citation_dois": [
            "10.1093/bioinformatics/btr539",
            "10.1093/bioinformatics/btq033",
        ],
        "citation_references": [
            "Dale et al. (2011), DOI 10.1093/bioinformatics/btr539",
            "Quinlan & Hall (2010), DOI 10.1093/bioinformatics/btq033",
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT Python wrapper exposes bounded local BED intersection; it still "
            "requires a separately installed BEDTools executable and does not expose arbitrary "
            "flags."
        ),
    },
    "pyfastx": {
        "api": "dnakit.integrations.summarize_pyfastx",
        "functions": ["dnakit.integrations.summarize_pyfastx"],
        "citation_dois": ["10.1093/bib/bbaa368"],
        "citation_references": [
            "Du et al. (2021), DOI 10.1093/bib/bbaa368"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded FASTA/FASTQ aggregate statistics and a "
            "small preview; pyfastx random-access indexes and write operations are not re-exported."
        ),
    },
    "dnaio": {
        "api": "dnakit.integrations.summarize_dnaio_fastx",
        "functions": ["dnakit.integrations.summarize_dnaio_fastx"],
        "sources": [
            "https://doi.org/10.5281/zenodo.10548864",
            "https://github.com/marcelm/dnaio",
        ],
        "citation_dois": ["10.5281/zenodo.10548864"],
        "citation_references": ["dnaio software archive, Zenodo DOI 10.5281/zenodo.10548864"],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded read-only FASTA/FASTQ stream statistics; "
            "paired-end writers, format conversion and the complete dnaio API are not re-exported."
        ),
    },
    "mappy": {
        "api": "dnakit.integrations.map_mappy_sequences",
        "functions": ["dnakit.integrations.map_mappy_sequences"],
        "sources": ["https://github.com/lh3/minimap2"],
        "citation_dois": ["10.1093/bioinformatics/bty191"],
        "citation_references": [
            "Li (2018), minimap2, DOI 10.1093/bioinformatics/bty191"
        ],
        "license_evidence": ["MIT (minimap2/mappy upstream)"],
        "boundary": (
            "The optional MIT mappy binding exposes bounded local mapping of a FASTA/FASTQ stream "
            "against a caller-prepared reference; index writing and arbitrary minimap2 flags remain external."
        ),
    },
    "gffutils": {
        "api": "dnakit.integrations.query_gffutils_database",
        "functions": ["dnakit.integrations.query_gffutils_database"],
        "citation_references": [
            "gffutils documentation and daler contributors (MIT)"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes read-only queries against an existing gffutils "
            "FeatureDB; database creation, SQL and feature mutation are not re-exported."
        ),
    },
    "sbol3": {
        "api": "dnakit.integrations.summarize_sbol3_document",
        "functions": ["dnakit.integrations.summarize_sbol3_document"],
        "citation_dois": ["10.1021/acssynbio.2c00164"],
        "citation_references": [
            "Mitchell, Beal & Bartley (2022), DOI 10.1021/acssynbio.2c00164"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes a bounded, read-only inventory of a "
            "local SBOL3 document; it does not resolve remote RDF, mutate objects, or write files."
        ),
    },
    "scadnano": {
        "api": "dnakit.integrations.summarize_scadnano_design",
        "functions": ["dnakit.integrations.summarize_scadnano_design"],
        "citation_references": [
            "Doty, Lee & Stérin (2020), scadnano, arXiv:2005.11841"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes a bounded, read-only inventory of a local "
            "scadnano .sc design; it does not edit, render, or download design assets."
        ),
    },
    "polars-bio": {
        "api": "dnakit.integrations.overlap_polars_bio_bed",
        "functions": ["dnakit.integrations.overlap_polars_bio_bed"],
        "citation_dois": ["10.1093/bioinformatics/btaf640"],
        "citation_references": [
            "polars-bio, Bioinformatics, DOI 10.1093/bioinformatics/btaf640"
        ],
        "license_evidence": ["Apache-2.0"],
        "boundary": (
            "The optional Apache-2.0 backend exposes bounded overlap of local BED/CSV interval "
            "files; arbitrary SQL, Python expressions, remote URLs and the full DataFrame API "
            "are not re-exported."
        ),
    },
    "dnacauldron": {
        "api": "dnakit.integrations.simulate_dnacauldron_assembly",
        "functions": ["dnakit.integrations.simulate_dnacauldron_assembly"],
        "citation_references": [
            "DnaCauldron documentation and Edinburgh Genome Foundry contributors; "
            "no DOI supplied by the survey"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes bounded Gibson and Type IIS assembly simulation; "
            "it reports predicted constructs and flaws but does not infer a wet-lab result."
        ),
    },
    "dna-features-viewer": {
        "api": "dnakit.integrations.summarize_dna_features_viewer",
        "functions": ["dnakit.integrations.summarize_dna_features_viewer"],
        "citation_references": [
            "DNA Features Viewer documentation and Edinburgh Genome Foundry contributors; "
            "no DOI supplied by the survey"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes a bounded, read-only feature inventory for a "
            "local GenBank or GFF record; plotting, writing files and the full Biopython "
            "translator API are not re-exported."
        ),
    },
    "gtars": {
        "api": "dnakit.integrations.summarize_gtars_regions",
        "functions": ["dnakit.integrations.summarize_gtars_regions"],
        "citation_references": [
            "gtars documentation and databio contributors; no DOI supplied by the survey"
        ],
        "license_evidence": ["BSD-2-Clause"],
        "boundary": (
            "The optional BSD-2-Clause backend exposes a bounded, read-only inventory of a "
            "local BED RegionSet; URLs, binary caches, tokenizers and arbitrary interval "
            "expressions are not re-exported."
        ),
    },
    "geneblocks": {
        "api": "dnakit.integrations.diff_geneblocks_sequences",
        "functions": ["dnakit.integrations.diff_geneblocks_sequences"],
        "citation_references": [
            "GeneBlocks documentation and Edinburgh Genome Foundry contributors; "
            "no DOI supplied by the survey"
        ],
        "license_evidence": ["MIT"],
        "boundary": (
            "The optional MIT backend exposes a bounded, deterministic DNA sequence diff; "
            "BLAST acceleration, plotting, feature transfer and the full GeneBlocks API "
            "are not re-exported."
        ),
    },
    "sourmash": {
        "api": "dnakit.integrations.compare_sourmash_sequences",
        "functions": ["dnakit.integrations.compare_sourmash_sequences"],
        "citation_dois": ["10.21105/joss.06830"],
        "citation_references": [
            "Pierce-Ward et al. (2024), sourmash v4, DOI 10.21105/joss.06830"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "The optional BSD-3-Clause backend exposes a bounded, local in-memory DNA "
            "MinHash comparison; databases, taxonomy indexes, plotting and network queries "
            "are not re-exported."
        ),
    },
    "pyfaidx": {
        "api": "dnakit.integrations.fetch_pyfaidx_region",
        "functions": ["dnakit.integrations.fetch_pyfaidx_region"],
        "sources": [
            "https://doi.org/10.7287/peerj.preprints.970v1",
            "https://github.com/mdshw5/pyfaidx",
        ],
        "citation_dois": ["10.7287/peerj.preprints.970v1"],
        "citation_references": [
            "Shirley et al. (2015), DOI 10.7287/peerj.preprints.970v1"
        ],
        "license_evidence": ["BSD-3-Clause"],
        "boundary": (
            "The optional BSD-3-Clause backend exposes bounded, read-only FASTA interval "
            "fetching; pyfaidx's complete API is not re-exported."
        ),
    },
    "vcfpy": {
        "api": "dnakit.integrations.summarize_vcf_vcfpy",
        "functions": ["dnakit.integrations.summarize_vcf_vcfpy"],
        "boundary": (
            "The optional MIT backend exposes bounded VCF summary parsing; arbitrary record "
            "mutation is not available through the JSON boundary."
        ),
    },
    "cyvcf2": {
        "api": "dnakit.integrations.summarize_vcf_cyvcf2",
        "functions": ["dnakit.integrations.summarize_vcf_cyvcf2"],
        "boundary": (
            "The optional MIT backend exposes bounded VCF summary parsing; cyvcf2's complete "
            "API and genotype arrays are not re-exported."
        ),
    },
    "dna chisel": {
        "api": "dnakit.integrations.optimize_with_dnachisel",
        "functions": ["dnakit.integrations.optimize_with_dnachisel"],
        "boundary": (
            "DNAKit also has native codon/rule optimization, but this route invokes DNA Chisel "
            "only when its optional backend is installed."
        ),
    },
    "dnachisel": {
        "api": "dnakit.integrations.optimize_with_dnachisel",
        "functions": ["dnakit.integrations.optimize_with_dnachisel"],
        "boundary": (
            "The MIT DNA Chisel package is an optional backend; this adapter allow-lists "
            "built-in specifications and does not copy its source."
        ),
    },
    "primer3-py": {
        "api": "dnakit.molbio.Primer3CLIDesignAdapter",
        "functions": ["dnakit.molbio.Primer3CLIDesignAdapter"],
        "boundary": "Primer3 is an external GPL backend; DNAKit does not bundle primer3-py.",
    },
    "seqme": {
        "api": "dnakit.integrations.evaluate_seqme",
        "functions": ["dnakit.integrations.evaluate_seqme"],
        "boundary": (
            "The BSD-3-Clause seqme package is an optional backend; model extras, weights and "
            "their licenses remain separate."
        ),
    },
    "ilearnplus": {
        "api": "dnakit.external_adapters.run",
        "functions": [],
        "boundary": (
            "Academic-use language and redistribution terms are not sufficiently clear for a "
            "bundled adapter; use an explicitly installed external backend."
        ),
    },
    "ifeatureomega-cli": {
        "api": "dnakit.external_adapters.run",
        "functions": [],
        "boundary": (
            "The source survey does not establish a redistribution license; use an explicit "
            "external backend after author/license review."
        ),
    },
    "ifeatureomega": {
        "api": "dnakit.external_adapters.run",
        "functions": [],
        "boundary": (
            "The source survey does not establish a redistribution license; use an explicit "
            "external backend after author/license review."
        ),
    },
    "ga": {
        "api": "dnakit.external_adapters.run",
        "functions": [],
        "boundary": (
            "The GPLv3 ga package is a general-purpose genetic-algorithm dependency, not a "
            "DNAKit-native DNA analysis API; use it only as a separately installed backend."
        ),
    },
    "ugene": {
        "api": "dnakit.external_adapters.run",
        "functions": [],
        "boundary": "UGENE remains a separately installed GPLv2 external program.",
    },
}


@lru_cache(maxsize=1)
def _data() -> dict[str, Any]:
    if not _CATALOG.is_file():
        return {"summary": {}, "records": []}
    return json.loads(_CATALOG.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _capability_data() -> dict[str, Any]:
    if not _CAPABILITIES.is_file():
        return {"count": 0, "capabilities": []}
    value = json.loads(_CAPABILITIES.read_text(encoding="utf-8"))
    capabilities = value.get("capabilities", [])
    if isinstance(capabilities, dict):
        capabilities = list(capabilities.values())
    return {**value, "capabilities": capabilities}


def _matches(record: dict[str, Any], terms: list[str]) -> bool:
    text = json.dumps(record, ensure_ascii=False).casefold()
    return all(term in text for term in terms)


def list_third_party_tools(
    query: str = "",
    *,
    kind: str | None = None,
    decision: str | None = None,
    integration_status: str | None = None,
    include_unverified: bool = True,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """List survey projects and their current DNAKit coverage evidence.

    ``local_match`` means that a DNAKit registry, dependency, Skill notice or
    analysis entry was found.  It does not mean that every upstream method is
    implemented locally, nor that a backend, model or database is installed.
    """

    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    if kind is not None and not isinstance(kind, str):
        raise ConfigurationError("kind must be a string or None.")
    if decision is not None and not isinstance(decision, str):
        raise ConfigurationError("decision must be a string or None.")
    if integration_status is not None and not isinstance(integration_status, str):
        raise ConfigurationError("integration_status must be a string or None.")
    terms = query.casefold().split()
    result = []
    for record in _data().get("records", []):
        if kind is not None and kind not in record.get("kind", []):
            continue
        if decision is not None and record.get("primary_decision") != decision:
            continue
        if (
            integration_status is not None
            and record.get("integration_status") != integration_status
        ):
            continue
        if not include_unverified and record.get("primary_decision") in {
            "blocked_or_unverified",
            "review_needed",
        }:
            continue
        if _matches(record, terms):
            result.append(copy.deepcopy(record))
    return result


def describe_third_party_tool(tool: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return one survey record by name, key, URL or normalized name."""

    if not isinstance(tool, str) or not tool.strip():
        raise ConfigurationError("tool must be a non-empty project name, key or URL.")
    query = tool.casefold().strip()
    records = _data().get("records", [])
    normalized = _name_key(query)
    for record in records:
        candidates = {
            str(record.get("key", "")).casefold(),
            str(record.get("name", "")).casefold(),
            _name_key(str(record.get("name", ""))),
            *(str(url).casefold() for url in record.get("urls", [])),
        }
        if query in candidates or normalized in candidates:
            return copy.deepcopy(record)
    raise ConfigurationError(f"Unknown surveyed DNA tool: {tool!r}.")


def third_party_coverage(*, verbose: bool = True) -> dict[str, Any]:
    """Return counts for the dated public-tool survey and local cross-reference."""

    data = _data()
    return {
        "schema_version": data.get("schema_version"),
        "source": data.get("source"),
        "summary": copy.deepcopy(data.get("summary", {})),
    }


def integration_inventory(*, verbose: bool = True) -> dict[str, Any]:
    """Count DNAKit's Python, MCP, Skill and survey surfaces separately.

    A Python function, an MCP server recipe, a Skill and a survey project are
    different integration units. Candidate schemas are reported separately
    from executable functions, so this inventory does not imply full parity
    with every upstream project's namespace.
    """

    from dnakit.tools.registry import default_tool_registry

    python_stats = default_tool_registry(verbose=False).stats()
    mcp_servers = []
    if _MCP_CATALOG.is_file():
        mcp_servers = json.loads(_MCP_CATALOG.read_text(encoding="utf-8")).get("servers", [])
    skills = []
    if _SKILL_CATALOG.is_file():
        skills = json.loads(_SKILL_CATALOG.read_text(encoding="utf-8")).get("skills", [])
    coverage = third_party_coverage(verbose=False)
    summary = coverage.get("summary", {})
    candidate_data = _capability_data()
    candidate_statuses = dict(
        Counter(str(item.get("status", "")) for item in candidate_data.get("capabilities", []))
    )
    python_functions = python_stats["total"]
    mcp_server_recipes = len(mcp_servers)
    registered_skills = len(skills)
    return {
        "python_functions": python_functions,
        "agent_compatible_python_functions": python_stats["agent_compatible"],
        "python_only_functions": python_stats["python_only"],
        "python_categories": python_stats["categories"],
        "mcp_server_recipes": mcp_server_recipes,
        "registered_skills": registered_skills,
        "tool_surface_total": python_functions + mcp_server_recipes + registered_skills,
        "tool_surface_definition": (
            "Raw sum of public Python functions, MCP server recipes and registered Skills; "
            "the same upstream capability can appear in more than one surface."
        ),
        "survey_unique_projects": summary.get("unique_records", 0),
        "survey_local_matches": summary.get("local_matches", 0),
        "survey_records_by_kind": summary.get("kinds", {}),
        "survey_records_by_integration_status": summary.get("integration_status", {}),
        "unified_candidate_capabilities": candidate_data.get(
            "count", len(candidate_data.get("capabilities", []))
        ),
        "unified_candidate_statuses": candidate_statuses,
        "source": {
            "python": "dnakit.tools.catalog / AgentToolRegistry",
            "mcp": str(_MCP_CATALOG),
            "skills": str(_SKILL_CATALOG),
            "survey": _data().get("source"),
            "candidates": candidate_data.get("source"),
        },
    }


def list_unified_capabilities(
    query: str = "",
    *,
    module: str | None = None,
    status: str | None = None,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """List the 837 candidate capabilities with their formal DNAKit API names."""

    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    if module is not None and not isinstance(module, str):
        raise ConfigurationError("module must be a string or None.")
    if status is not None and not isinstance(status, str):
        raise ConfigurationError("status must be a string or None.")
    terms = query.casefold().split()
    result = []
    for capability in _capability_data().get("capabilities", []):
        if module is not None and capability.get("module") != module:
            continue
        if status is not None and capability.get("status") != status:
            continue
        if _matches(capability, terms):
            result.append(copy.deepcopy(capability))
    return result


def describe_unified_capability(capability_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return one candidate capability, including input and output schemas."""

    if not isinstance(capability_id, str) or not capability_id.strip():
        raise ConfigurationError("capability_id must be a non-empty string.")
    for capability in _capability_data().get("capabilities", []):
        if capability.get("capability_id") == capability_id:
            return copy.deepcopy(capability)
    raise ConfigurationError(f"Unknown unified capability: {capability_id!r}.")


def unified_capability_summary(*, verbose: bool = True) -> dict[str, Any]:
    """Return the count and source of the shipped unified capability schemas."""

    data = _capability_data()
    capabilities = data.get("capabilities", [])
    modules: dict[str, int] = {}
    statuses: dict[str, int] = {}
    for capability in capabilities:
        modules[capability.get("module", "")] = modules.get(capability.get("module", ""), 0) + 1
        statuses[capability.get("status", "")] = statuses.get(capability.get("status", ""), 0) + 1
    return {
        "source": data.get("source"),
        "count": data.get("count", len(capabilities)),
        "modules": modules,
        "statuses": statuses,
    }


def run_unified_capability(
    capability_id: str,
    arguments: Mapping[str, Any] | None = None,
    *,
    allow_side_effects: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Execute one registered unified capability through the DNAKit Agent registry.

    The capability catalog is a discovery contract, while the registry is the
    execution authority.  Reusing the registry means the same JSON argument
    coercion, output bounds, denied-tool rules and side-effect confirmation are
    applied to Python, MCP and Skills callers.  ``arguments`` normally follows
    the actual ``dnakit.<module>.<operation>`` signature.  Candidate adapters
    use ``{"input": ..., "options": ...}``; when a native compatibility
    function exposes one different required parameter (for example ``source``),
    ``input`` is mapped to that parameter without changing the public route.
    """

    if not isinstance(arguments, Mapping) and arguments is not None:
        raise ConfigurationError("arguments must be a mapping or None.")
    capability = describe_unified_capability(capability_id, verbose=False)
    supplied = {} if arguments is None else dict(arguments)
    if any(not isinstance(key, str) for key in supplied):
        raise ConfigurationError("arguments must use string keys.")

    from dnakit.tools.registry import default_tool_registry

    registry = default_tool_registry(verbose=False)
    tool = registry.resolve(capability["api"])
    parameter_names = {parameter.name for parameter in tool.call_plan.parameters}
    required_names = {
        parameter.name
        for parameter in tool.call_plan.parameters
        if parameter.required and parameter.name != "verbose"
    }
    call_arguments = dict(supplied)
    if "input" in call_arguments and "input" not in parameter_names:
        if len(required_names) != 1:
            raise ConfigurationError(
                "This native compatibility route does not have one unambiguous input parameter.",
                context={"api": capability["api"], "parameters": sorted(parameter_names)},
            )
        input_name = next(iter(required_names))
        call_arguments[input_name] = call_arguments.pop("input")
    options = call_arguments.get("options")
    if "options" not in parameter_names and options not in (None, {}):
        raise ConfigurationError(
            "This native compatibility route does not accept an options mapping.",
            context={"api": capability["api"]},
        )
    if "options" not in parameter_names:
        call_arguments.pop("options", None)
    if "verbose" in parameter_names and "verbose" not in call_arguments:
        call_arguments["verbose"] = verbose

    result = registry.execute(
        capability["api"],
        call_arguments,
        allow_side_effects=allow_side_effects,
    )
    return {
        "status": "completed",
        "capability_id": capability["capability_id"],
        "api": capability["api"],
        "provenance": {
            "catalog": _capability_data().get("source"),
            "implementation_status": capability.get("status"),
            "execution": "dnakit.tools.registry.AgentToolRegistry",
        },
        "result": result,
    }


def third_party_route(tool: str, *, verbose: bool = True) -> dict[str, Any]:
    """Describe the unified DNAKit route for a surveyed project.

    The result is a route description, not an implicit execution request.  A
    caller must explicitly invoke ``run_analysis`` or an MCP/Skill API after
    checking dependencies, credentials, service terms and provenance.
    """

    record = describe_third_party_tool(tool, verbose=False)
    native = _NATIVE_FACADES.get(record["name"].casefold())
    if native is None:
        native = _NATIVE_FACADES.get(_name_key(record["name"]))
    matches = record.get("local_matches", [])
    route: dict[str, Any] = {
        "tool": record["name"],
        "status": record.get("integration_status", "no_local_match"),
        "decision": record.get("primary_decision", "review_needed"),
        "sources": record.get("urls", []),
        "citation_dois": record.get("citation_dois", []),
        "citation_references": record.get("citation_references", record.get("citation_dois", [])),
        "license_evidence": record.get("licenses", []),
        "provenance": {
            "survey_lines": record.get("source_lines", []),
            "survey_sections": record.get("sections", []),
        },
        "related_capabilities": list(record.get("related_capabilities", [])),
    }
    analysis_ids = [item["value"] for item in matches if item.get("kind") == "analysis_adapter"]
    mcp_ids = [item["value"] for item in matches if item.get("kind") == "mcp_registry"]
    skill_ids = [item["value"] for item in matches if item.get("kind") == "skill_registry"]
    if analysis_ids:
        route.update(api="dnakit.integrations.run_analysis", analysis_ids=analysis_ids)
    elif mcp_ids:
        route.update(api="dnakit.integrations.call_mcp_tool", mcp_server_ids=mcp_ids)
    elif skill_ids:
        route.update(api="dnakit.skills / dnakit.workflows", skill_ids=skill_ids)
    elif record.get("integration_status") == "project_evidence_only":
        route.update(api="dnakit.external_adapters.run", status="reference_only")
    else:
        route.update(api="dnakit.external_adapters.run", status="unregistered_backend")
    if native is not None:
        capability_by_api = {
            item.get("api"): item.get("capability_id")
            for item in _capability_data().get("capabilities", [])
        }
        facade_status = "native_facade"
        if native["api"] == "dnakit.external_adapters.run":
            facade_status = "external_facade"
        elif native["api"].startswith("dnakit.integrations"):
            facade_status = "optional_backend"
        route.update(
            api=native["api"],
            sources=list(native.get("sources", route["sources"])),
            native_functions=list(native["functions"]),
            boundary=native["boundary"],
            status=facade_status,
            citation_dois=list(native.get("citation_dois", route["citation_dois"])),
            citation_references=list(
                native.get("citation_references", route["citation_references"])
            ),
            license_evidence=list(native.get("license_evidence", route["license_evidence"])),
            related_capabilities=[
                capability_by_api[function]
                for function in native["functions"]
                if function in capability_by_api
            ],
        )
    return route


__all__ = [
    "describe_third_party_tool",
    "describe_unified_capability",
    "integration_inventory",
    "list_third_party_tools",
    "list_unified_capabilities",
    "run_unified_capability",
    "third_party_coverage",
    "third_party_route",
    "unified_capability_summary",
]
