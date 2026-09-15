
---

# Polyploid Gene Family Analysis & Seamless Cloning Pipeline

### AI Agent Skill · v3.0.0

A genome-aware, polyploid-aware workflow for **gene family identification, candidate prioritization, validated full-length CDS cloning, homoeolog-aware primer design, and end-to-end construct verification**.

---

## Overview

This repository provides a specialized AI Agent Skill for plant genomics and molecular cloning workflows, particularly those involving **polyploid genomes such as *Medicago sativa* L.**

Version 3.0.0 inherits the sequence-validation framework established in v2.0.0 and substantially extends it from **primer-level QC to construct-level validation**.

The pipeline is designed to maintain sequence fidelity across the entire workflow:

```text
Gene Identity
      ↓
Genome Locus
      ↓
Transcript
      ↓
CDS
      ↓
Protein / ORF
      ↓
Construction Intent
      ↓
Vector Architecture
      ↓
Cloning Primers
      ↓
PCR Product
      ↓
Recombinant Construct
      ↓
Final ORF
      ↓
Final Protein
```

The central objective is not simply to generate a primer pair that can amplify a target sequence.

> **The objective is to computationally demonstrate that the proposed primers will generate the intended biological molecule in the intended vector architecture before primer synthesis.**

---

# What's New in v3.0.0?

Version 2.0.0 established a strict **Gene → CDS → Primer** validation chain.

Version 3.0.0 extends this into a complete:

> **Gene → CDS → Construction → Primer → PCR Product → Construct → Protein**

validation framework.

### Major upgrades

* **Construction Intent Lock**
* **Terminal ORF Integrity Gate**
* **Reverse Primer Boundary Audit**
* **Partial Codon Detection**
* **STOP Codon Integrity Control**
* **Actual Fusion Architecture Verification**
* **PCR Product Reconstruction**
* **Recombinant Construct Reconstruction**
* **Final Protein Reconstruction**
* **Unexpected Vector-Derived Tail Detection**
* **Homoeolog Specificity Feasibility Analysis**
* **Explicit “Physically Unresolvable” Classification**
* **Current vs. Obsolete Design Detection**
* **Mandatory GO / NO-GO / HOLD decision**

These additions are specifically designed to prevent a class of errors in which a primer pair appears structurally valid but produces a **truncated, frame-altered, or unintended protein product**.

---

# Core Principles

## 1. Never design directly from a gene name

The Agent must never design cloning primers directly from:

* gene symbols
* annotation descriptions
* phylogenetic rankings
* protein names
* old primer names
* historical experimental notes
* plasmid filenames

Instead, the identity must be resolved through:

```text
Gene Symbol
    ↓
Genome Gene ID
    ↓
Transcript ID
    ↓
CDS Sequence
    ↓
Protein Sequence
    ↓
Validated ORF
```

If the identity chain cannot be established, primer design stops.

---

## 2. Sequence is the ultimate reference

When conflicting information exists, the Skill prioritizes actual sequence evidence over historical design descriptions.

The default evidence hierarchy is:

```text
Raw nucleotide sequence
        ↓
GFF3-derived sequence
        ↓
Validated feature annotation
        ↓
Current primer record
        ↓
Current simulated construct
        ↓
Historical primer design
        ↓
Historical reports / notes
        ↓
Filename or project label
```

This prevents obsolete `.dna` files, primer names, or project descriptions from silently overriding the actual sequence.

---

## 3. Full-length cloning means full-length CDS

For a task explicitly defined as **full-length CDS cloning**, the default invariant is:

```text
PCR gene-derived sequence
        =
Complete validated CDS
```

This includes:

* complete coding sequence
* complete final amino-acid codon
* native STOP codon

The Agent must never truncate the CDS merely to improve primer convenience or apparent homoeolog specificity.

---

# Six Mandatory QC Gates

Version 3.0.0 introduces six sequential validation gates.

| Gate                                                  | Purpose                                                   |
| ----------------------------------------------------- | --------------------------------------------------------- |
| **GATE 1 — Gene Identity Lock**                       | Confirm exact gene, locus, and transcript                 |
| **GATE 2 — CDS / ORF Integrity Lock**                 | Confirm complete and valid coding sequence                |
| **GATE 3 — Construction Intent Lock**                 | Determine the intended protein architecture               |
| **GATE 4 — Primer Boundary & Strand Lock**            | Verify primer orientation and exact CDS boundaries        |
| **GATE 5 — PCR Product Reconstruction**               | Reconstruct the actual expected amplicon                  |
| **GATE 6 — Final Construct / Protein Reconstruction** | Verify recombinant construct and final translated protein |

A primer cannot be released for synthesis until all mandatory gates pass.

---

# Key Features

## Gene Identity Lock

The Agent generates a structured identity record containing:

* Gene ID
* Transcript ID
* chromosome / scaffold
* genomic coordinates
* strand
* CDS length
* protein length
* source sequence

This is particularly important in polyploid genomes where multiple homoeologs may have highly similar names and sequences.

---

## CDS Validation Gate

The Skill validates:

* CDS length
* coding frame
* start codon
* terminal codon
* STOP codon
* internal premature STOP codons
* ORF integrity

When both FASTA and GFF3 are available, they are independently checked for consistency.

---

## Construction Intent Lock

Version 3.0.0 explicitly separates:

### Standalone protein expression

```text
Promoter
   ↓
Target CDS
   ↓
STOP
```

### N-terminal fusion

```text
Promoter
   ↓
Fusion partner
   ↓
Target CDS
   ↓
STOP
```

### C-terminal fusion

```text
Promoter
   ↓
Target CDS
   ↓
Linker
   ↓
Fusion partner
   ↓
STOP
```

The Agent must verify the **actual vector architecture** rather than infer fusion from a plasmid filename or historical note.

In particular:

> The presence of GFP elsewhere in a plasmid does not prove that the target CDS is intended to be fused with GFP.

---

## Terminal ORF Integrity Gate

This is a major v3.0.0 addition.

Every reverse cloning primer is explicitly mapped to the reference CDS.

The Agent must determine:

```text
Reverse primer
      ↓
Reverse complement
      ↓
Exact CDS coordinates
      ↓
PCR product endpoint
      ↓
Final codon
      ↓
STOP
```

The system specifically detects **partial final codons**, a failure mode that may otherwise remain hidden during routine primer checking.

---

## Seamless Cloning Boundary Verification

For Gibson / In-Fusion-style cloning, every primer is decomposed into:

```text
[Vector Homology Arm]
+
[Gene-Specific Annealing Region]
```

The Agent verifies:

* vector linearization site
* left and right vector boundaries
* homology-arm orientation
* effective homology length
* gene-specific binding coordinates
* restriction-site conflicts
* final junction sequences

---

## Homoeolog-Aware Primer Design

For polyploid genomes, the Skill evaluates primer specificity against relevant homoeologs.

It distinguishes between:

* target-specific
* partially specific
* non-specific
* physically unresolvable

### Important design principle

Homoeolog specificity and CDS completeness are treated as **separate constraints**.

If no primer can simultaneously provide:

```text
Complete CDS
+
Target homoeolog specificity
```

the Agent must report the limitation rather than damage the CDS to manufacture specificity.

A valid alternative is:

```text
Full-length amplification
        ↓
Cloning
        ↓
Sanger sequencing
        ↓
Homoeolog identification
        ↓
Select correct clone
```

---

## PCR Product Reconstruction

Version 3.0.0 does not stop at primer analysis.

The Agent reconstructs the expected PCR product from:

```text
Forward primer
+
Reference CDS
+
Reverse primer
```

It then reports:

* PCR product length
* CDS coverage
* 5′ truncation
* 3′ truncation
* internal mismatches
* indels
* terminal codon status

For full-length cloning:

> The gene-derived region must equal the validated CDS.

---

## Recombinant Construct Reconstruction

The Agent reconstructs the predicted recombinant plasmid:

```text
Vector Left Fragment
        +
PCR Product
        +
Vector Right Fragment
        ↓
Predicted Recombinant Construct
```

If a SnapGene or equivalent simulated construct is available, the independent reconstruction is compared against it.

This creates an additional safeguard against errors hidden inside simulated construct files.

---

## Final Protein Reconstruction

The final recombinant ORF is translated computationally.

The Agent compares:

```text
Reference Protein
        vs.
Predicted Construct Protein
```

It detects:

* altered amino acids
* missing amino acids
* premature STOP
* missing STOP
* frame shifts
* vector-derived residues
* unintended linker residues
* unintended C-terminal tails
* accidental fusion

For standalone full-length expression:

```text
Predicted Protein == Reference Protein
```

must hold unless an intentional modification has been explicitly specified.

---

# STOP Codon Policy

For standalone full-length cloning:

> **Native STOP codon is required.**

STOP may be intentionally omitted only when the Agent independently verifies a genuine downstream fusion architecture:

```text
Target CDS
   ↓
Linker
   ↓
Fusion Partner
   ↓
STOP
```

The Agent must verify that the fusion partner:

* actually exists downstream
* is in the correct orientation
* is in-frame
* has the intended linker
* contains an appropriate terminal STOP

An old annotation such as:

```text
"STOP omitted for GFP fusion"
```

is not sufficient evidence.

---

# Error-Resistant Primer Release

Before recommending synthesis, the Agent performs a final release audit.

### Mandatory checks

```text
[✓] Gene identity verified
[✓] Transcript verified
[✓] CDS verified
[✓] ORF verified
[✓] Construction intent verified
[✓] Vector architecture verified
[✓] Fwd primer verified
[✓] Rev primer verified
[✓] Full CDS coverage verified
[✓] Final amino-acid codon verified
[✓] STOP / fusion architecture verified
[✓] Homology arms verified
[✓] Homoeolog specificity assessed
[✓] PCR product reconstructed
[✓] Recombinant construct reconstructed
[✓] Final protein reconstructed
[✓] No unexpected non-native sequence
```

If a mandatory check fails:

> **The Agent must not release the primer for synthesis.**

---

# Decision System

Every final cloning task receives exactly one decision:

### GO

```text
GO — Primer is sequence-validated and ready for synthesis
```

All mandatory QC gates pass.

### NO-GO

```text
NO-GO — Primer design contains a correctable design error
```

The design must be revised before synthesis.

### HOLD

```text
HOLD — Source sequence, construct intent, or vector architecture remains unresolved
```

The Agent must not make assumptions to force a GO decision.

---

# Standard Output

The Skill generates structured outputs for reproducibility and downstream laboratory use.

### Core records

1. `Gene_ID_Record.tsv`
2. `CDS_validation_report.tsv`
3. `Seamless_Primers.tsv`
4. `Primer_Audit_Report.tsv`

### v3.0.0 extended QC

5. `<Gene>_Cloning_Primer_Final_QC_Report.md`
6. `PCR_Product_Reconstruction`
7. `Construct_Reconstruction`
8. `Final_Protein_QC`

Where computational scripts are used, they should be preserved in a dedicated `qc/` directory together with machine-readable results where practical.

---

# Usage

Import `prompt.md` into a compatible local Agent framework as a **System Prompt / AI Skill**.

Typical environments include:

* Dify
* FastGPT
* Coze
* LangChain-based Agent systems
* local AI Agent frameworks

The Skill can be supplied with:

* genome FASTA
* CDS FASTA
* protein FASTA
* GFF3/GTF annotation
* vector `.dna` / sequence files
* simulated constructs
* primer records
* experimental design notes

For best results, provide both the **reference sequence resources** and the **actual vector/construct files** used for the planned experiment.

---

# Recommended Workflow

A typical cloning task should proceed as:

```text
1. Identify candidate gene
        ↓
2. Lock genome locus and transcript
        ↓
3. Validate CDS
        ↓
4. Validate protein / ORF
        ↓
5. Lock construction intent
        ↓
6. Inspect actual vector architecture
        ↓
7. Analyze homoeolog specificity
        ↓
8. Design cloning primers
        ↓
9. Audit primer boundaries
        ↓
10. Reconstruct PCR product
        ↓
11. Reconstruct recombinant construct
        ↓
12. Translate final ORF
        ↓
13. Compare final protein
        ↓
14. GO / NO-GO / HOLD
        ↓
15. Synthesis
```

---

# What v3.0.0 Prevents

The upgraded architecture is specifically designed to catch errors such as:

* selecting the wrong homoeolog
* selecting the wrong transcript
* amplifying an incomplete CDS
* truncating the final amino-acid codon
* unintentionally deleting the STOP codon
* assuming GFP is a fusion partner when it is independently expressed
* introducing vector-derived C-terminal residues
* creating an unintended frame shift
* using an obsolete primer design
* using the wrong homology arm
* accepting a simulated construct without independently reconstructing it
* sacrificing CDS integrity to obtain homoeolog specificity

---

# Design Philosophy

The Skill follows a simple principle:

> **Sequence correctness must be demonstrated at every level of the construct, not assumed from the primer design.**

For a standard full-length cloning experiment:

```text
Reference CDS
      =
PCR Gene-Derived Sequence
      =
Recombinant Target ORF
      =
Intended Protein
```

If these four layers do not agree, the design is not ready for synthesis.

---

# Version History

## v3.0.0 — Construct-Level Validation

### Inherited from v2.0.0

* Gene Identity Lock
* CDS validation
* genome-aware sequence analysis
* polyploid / homoeolog-aware analysis
* In-Fusion boundary verification
* vector compatibility assessment
* primer audit
* separation of cloning and qPCR workflows

### Introduced in v3.0.0

* Construction Intent Lock
* Terminal ORF Integrity Gate
* reverse-primer boundary audit
* partial-codon detection
* STOP-codon integrity control
* actual fusion architecture verification
* PCR product reconstruction
* recombinant construct reconstruction
* final protein reconstruction
* unexpected vector-derived tail detection
* homoeolog specificity feasibility analysis
* physically unresolvable specificity classification
* current/obsolete design detection
* mandatory GO / NO-GO / HOLD release system
* end-to-end construct-level QC report

---

# Scope and Limitations

This Skill provides **in silico sequence and construct validation**.

It does not guarantee:

* PCR amplification success
* cloning efficiency
* transformation efficiency
* plasmid stability
* protein expression
* protein function

Wet-lab validation remains necessary.

For polyploid organisms, computational specificity may also be fundamentally limited by sequence identity between homoeologs. When sequence-level discrimination is impossible, the Skill should explicitly recommend an experimentally appropriate downstream identification strategy rather than artificially compromising the target CDS.

---

## Disclaimer

AI-generated primers, homology arms, construct reconstructions, and sequence interpretations must be independently verified against the **actual laboratory vector, reference sequence, transcript annotation, and experimental objective** before synthesis or experimental use.

**Never synthesize a primer solely because an AI Agent reports that it is valid.**

