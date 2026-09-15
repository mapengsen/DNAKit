---
name: plant-genome-guided-primer-design
version: 2.0.0
description: >
  A genome-guided workflow for plant gene identity validation,
  homolog analysis, and experimental primer design.
  The workflow prevents incorrect primer design caused by
  gene ID misannotation, duplicated genes, paralogs, and polyploid genomes.

author: Li-Chia-ching
language: en

domain:
  - plant molecular biology
  - bioinformatics
  - genome annotation
  - primer design
  - gene expression analysis

supported_tasks:
  - gene_identity_validation
  - gene_sequence_confirmation
  - gene_family_context_analysis
  - full_length_cds_cloning_primer_design
  - infusion_cloning_primer_design
  - sqRT-PCR_primer_design
  - qRT-PCR_primer_design
  - homolog_specific_primer_design
  - polyploid_specificity_evaluation
  - PCR_simulation
  - scientific_report_generation

supported_species:
  - diploid_plants
  - autotetraploid_plants
  - allopolyploid_plants

input_requirements:
  required_files:
    - genome_fasta
    - annotation_gff3
    - cds_fasta
    - transcript_fasta

  optional_files:
    - protein_fasta
    - haplotype_sequences
    - gene_family_sequences
    - reference_species_sequences

input_types:
  - gene_id
  - gene_list
  - gene_family_name
  - nucleotide_sequence
  - protein_sequence

mandatory_workflow:
  - target_identity_validation
  - transcript_selection
  - homolog_analysis
  - experimental_objective_confirmation
  - primer_design
  - specificity_validation
  - pcr_simulation
  - report_generation

design_modes:
  cloning:
    enabled: true
    purpose: full_length_orf_amplification

  sqRT-PCR:
    enabled: true
    purpose: semi_quantitative_expression_detection

  qRT-PCR:
    enabled: true
    purpose: quantitative_expression_analysis

  homolog_specific:
    enabled: true
    purpose: distinguish_related_genes

critical_rules:
  - never_design_primers_directly_from_gene_id
  - always_validate_gene_identity_before_primer_design
  - always_check_gene_family_specificity_for_duplicate_genes
  - always_consider_polyploid_homologs_when_available
  - never_apply_qPCR_parameters_to_sqRT-PCR
  - stop_workflow_if_target_identity_fails

recommended_tools:
  sequence_analysis:
    - Biopython
    - BLAST+
    - HMMER
    - MAFFT

  primer_design:
    - Primer3
    - primer3-py

  report_generation:
    - Python
    - matplotlib
    - reportlab

outputs:
  - Gene_identity_validation.tsv
  - Gene_relationship_analysis.tsv
  - Validated_target_sequence.fasta
  - Cloning_primers.tsv
  - sqRT_PCR_primers.tsv
  - qRT_PCR_primers.tsv
  - Primer_specificity.tsv
  - PCR_simulation.tsv
  - Primer_Design_Report.pdf

failure_conditions:
  - incorrect_gene_annotation
  - missing_sequence_information
  - failed_identity_validation
  - unresolved_target_conflict

interaction_policy:
  ask_before_execution:
    - experimental_objective
    - target_gene_source
    - confidence_of_gene_identity

priority:
  biological_validation: highest
  primer_design: secondary
  automation: after_validation

---

## DNAKit execution

Use `dnakit skills doctor plant-genome-guided-primer-design` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


<!-- DNAKit metadata portability edit, 2026-09-12; see UPSTREAM.md. -->

# Plant Genome-Guided Gene Validation and Primer Design Skill

## 1. Skill Overview

This Skill provides a genome-guided workflow for reliable plant gene validation and molecular biology primer design.

The primary purpose is to prevent incorrect primer design caused by:

- incorrect gene identifiers
- outdated genome annotations
- duplicated genes
- paralog confusion
- polyploid genome complexity
- incorrect transcript selection

The workflow follows a mandatory validation-first strategy:

```

Candidate gene information
↓
Gene identity validation
↓
Sequence confirmation
↓
Gene family / homolog analysis
↓
Experimental objective definition
↓
Primer design
↓
Specificity validation
↓
PCR simulation
↓
Experimental report

```

The Skill supports:

- full-length CDS cloning primer design
- In-Fusion/Gibson assembly primer design
- semi-quantitative RT-PCR (sqRT-PCR) primer design
- qRT-PCR primer design
- homolog/paralog-specific primer design
- polyploid genome primer evaluation


---

# 2. Supported Species and Data Requirements

The workflow is optimized for plant genomes with:

Required:

```

Genome/

├── genome.fa
├── annotation.gff3
├── cds.fa
└── transcript.fa

```

Recommended:

```

├── protein.fa
├── haplotype_sequences/
├── gene_family_sequences/
└── reference_species_sequences/

```

Compatible with:

- diploid plants
- autotetraploid plants
- allopolyploid plants
- chromosome-level assemblies
- annotated draft genomes


---

# 3. User Input and Interactive Confirmation

Before starting analysis, the Agent should collect necessary information.

## 3.1 Required user information

Ask the user:

### Question 1

What is the experimental purpose?

Options:

```

A. Full-length gene cloning
B. Expression analysis
C. qRT-PCR quantification
D. Functional validation
E. Other

```


### Question 2

How was the target selected?

Options:

```

A. Gene ID from genome annotation
B. Gene family analysis
C. Transcriptome result
D. Published gene
E. User-provided sequence

```


### Question 3

Is the target gene identity already experimentally confirmed?

Options:

```

A. Yes
B. No
C. Unknown

```


If the target identity is uncertain:

The workflow MUST perform additional validation.


---

# 4. Input Files and Target Definition

The user may provide:

## Option A

Gene ID list

Example:

```

Gene001
Gene002

```


## Option B

Gene family name

Example:

```

ABC transporter
WRKY transcription factor

```

## Option C

Sequence file

Example:

```

target.fasta

```


Important:

A gene ID is treated as a candidate identifier only.

The Agent MUST NOT assume that the gene ID is biologically correct.


---

# 5. Mandatory Gene Identity Validation Module

## 5.1 Purpose

Confirm:

```

Is this genomic locus actually the biological target?

```

before primer design.


---

## 5.2 Annotation validation

Extract:

- Gene ID
- Transcript ID
- chromosome/scaffold
- strand
- exon number
- CDS length
- protein length
- annotation description


Generate:

```

Gene_identity_validation.tsv

```


Columns:

```

Gene_ID
Transcript_ID
Location
Strand
CDS_length
Protein_length
Annotation
Validation_status

```


---

# 6. Sequence-Based Identity Validation

The Agent should evaluate target identity using available evidence.


## 6.1 Protein/domain validation

If protein sequences exist:

Analyze:

- conserved domains
- functional motifs
- protein architecture
- similarity to known proteins


Methods:

Recommended:

- HMMER
- BLASTP
- local alignment
- motif scanning


The workflow should adapt to the target gene.

Examples:

- enzyme:
  catalytic domains

- transcription factor:
  DNA-binding domains

- transporter:
  membrane domains

- structural proteins:
  conserved regions


---

## 6.2 Homology validation

Compare against:

- curated reference proteins
- published homologs
- related species databases


Report:

```

Identity evidence:
Strong
Moderate
Weak
Failed

```


---

# 7. Validation Decision Gate

Primer design cannot start until validation is completed.


## Accepted

```

VALIDATED

```

Evidence:

- correct annotation
- sequence agreement
- expected functional characteristics


## Warning

```

UNCERTAIN

```

Possible causes:

- incomplete annotation
- multiple isoforms
- weak functional evidence


The Agent should ask user whether to continue.


## Failed

```

FAILED

```

Possible causes:

- incorrect gene ID
- missing sequence
- inconsistent annotation


The Agent must stop primer design.


---

# 8. Transcript Selection Module

For genes with multiple transcripts:

Evaluate:

- longest complete CDS
- canonical transcript
- experimentally supported isoform


Report:

```

Selected transcript:

Reason:

```


Do not automatically choose the first transcript in annotation.


---

# 9. Gene Family and Homolog Analysis

Mandatory for:

- gene families
- duplicated genes
- polyploid species
- paralog-rich genomes


Identify:

```

Target gene

├── Allelic variants
├── Homeologs
├── Paralogs
└── Related family members

```


Generate:

```

Gene_relationship_analysis.tsv

```


Columns:

```

Gene_ID
Similarity
Relationship
Protein_identity
Risk_level

```


---

# 10. Experimental Objective Selection

Primer design strategy depends on the purpose.


## A. Full-length cloning

Goal:

Obtain complete ORF.


## B. sqRT-PCR

Goal:

Detect transcript expression by agarose gel.


## C. qRT-PCR

Goal:

Quantitative expression measurement.


## D. Specific homolog detection

Goal:

Distinguish highly similar genes.


The Agent MUST select different primer rules.


---

# 11. Full-Length CDS Cloning Primer Design

Requirements:

Template:

validated CDS


Primer structure:

Forward:

```

Vector overlap
+
5' CDS-specific sequence

```


Reverse:

```

Vector overlap
+
3' CDS-specific sequence

```


Gene-specific region:

```

20-30 bp

```


Parameters:

```

Length:
18-30 nt

Tm:
58-65°C

GC:
40-60%

```


Mandatory checks:

- complete ORF coverage
- correct orientation
- no internal stop codon
- expected PCR product size


Output:

```

Cloning_primers.tsv

```


---

# 12. sqRT-PCR Primer Design

Purpose:

Semi-quantitative RT-PCR.


Amplicon:

```

200-500 bp

```

Preferred:

```

250-350 bp

```


Primer parameters:

```

Length:
18-25 nt

Tm:
58-62°C

GC:
40-60%

```


Prefer:

- unique transcript regions
- exon boundaries
- non-conserved regions


Avoid:

- conserved family domains
- repetitive sequences


Output:

```

sqRT_PCR_primers.tsv

```


---

# 13. qRT-PCR Primer Design

Amplicon:

```

80-200 bp

```


Preferred:

```

100-150 bp

```


Additional requirements:

- primer efficiency consideration
- minimal secondary structures
- consistent Tm


Output:

```

qRT_PCR_primers.tsv

```


---

# 14. Primer Specificity Validation

Mandatory for all primer types.


Databases:

Required:

```

cds.fa
transcript.fa

```


Recommended:

```

genome.fa

```


Check:

- exact matches
- mismatches
- off-target genes
- family members


Classification:

## Class I

Target-specific


## Class II

Possible homolog amplification


## Class III

Non-specific


Output:

```

Primer_specificity.tsv

```


---

# 15. Polyploid Genome Evaluation

For polyploid plants:

Evaluate:

- allele variation
- homeolog similarity
- SNP/InDel in primer binding sites


Report:

```

Universal
Allele-biased
Homeolog-specific
High-risk

```


---

# 16. PCR Simulation

Perform in silico PCR.


Check:

- primer binding position
- orientation
- expected fragment size
- multiple amplification products


Output:

```

PCR_simulation.tsv

```


---

# 17. Primer Ranking

Each gene should provide:

2-3 candidate primer pairs.


Suggested scoring:


|Category|Weight|
|-|-|
|Specificity|35%|
|Sequence validation confidence|25%|
|Primer quality|20%|
|Experimental suitability|20%|


---

# 18. Output Files

Generate:


```

Gene_identity_validation.tsv

Gene_relationship_analysis.tsv

Validated_target_sequence.fasta

Cloning_primers.tsv

sqRT_PCR_primers.tsv

qRT_PCR_primers.tsv

Primer_specificity.tsv

PCR_simulation.tsv

Primer_Design_Report.pdf

```


---

# 19. Final Report Requirements

Generate an A4 scientific report.


Include:


## 1. Project information

- species
- genome version
- objective


## 2. Target validation

- annotation evidence
- sequence evidence
- confidence level


## 3. Gene relationship analysis

- homologs
- paralogs
- polyploid risks


## 4. Primer design

- sequences
- parameters
- expected products


## 5. Specificity evaluation


## 6. PCR simulation


## 7. Experimental recommendations


---

# 20. Fail-Safe Rules

The Agent MUST NOT:


❌ Use gene ID directly for primer design

```

Gene_ID
↓
Primer

```


❌ Assume annotation is correct


❌ Ignore homologs in polyploid genomes


❌ Use qPCR parameters for sqRT-PCR


Correct workflow:


```

Gene candidate

↓

Identity validation

↓

Sequence confirmation

↓

Experimental objective

↓

Primer design

↓

Specificity validation

↓

PCR simulation

↓

Report

```
