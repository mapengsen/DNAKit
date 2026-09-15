# Plant Genome-Guided Primer Design Agent Skill


## Automated Gene Validation, Primer Design, and Expression Analysis Workflow


# Overview

Reliable molecular experiments in plants depend on accurate target gene identification and primer design.

However, primer design failures frequently originate from upstream problems:

- incorrect gene identifiers
- outdated genome annotations
- incomplete transcript models
- duplicated genes
- paralog interference
- polyploid genome complexity

This Skill provides a genome-guided AI-assisted workflow that prioritizes:

```
Target validation
        ↓
Sequence confirmation
        ↓
Homolog analysis
        ↓
Primer design
        ↓
Specificity validation
        ↓
PCR simulation
        ↓
Experimental report
```

The workflow is designed for plant functional genomics projects using genome assemblies, gene annotations, and molecular validation experiments.


---

# Features


## 1. Genome-guided gene identity validation

Unlike traditional primer design pipelines, this Skill does not directly convert a gene ID into primers.

Before primer design, the workflow verifies:

- gene annotation
- transcript structure
- CDS completeness
- protein characteristics
- sequence consistency

The system treats user-provided gene IDs as candidate identifiers rather than confirmed biological targets.


Example workflow:

```
Gene ID

↓

Genome annotation validation

↓

Sequence confirmation

↓

Functional evidence evaluation

↓

Primer design
```


This prevents errors caused by:

- incorrect gene annotation
- misidentified family members
- obsolete gene models


---

# 2. Multiple primer design modes


The Skill supports different experimental purposes.


## Full-length CDS cloning

For:

- In-Fusion cloning
- Gibson assembly
- overexpression vector construction

Features:

- complete ORF amplification
- start/stop codon coverage
- cloning-compatible primer architecture


---

## Semi-quantitative RT-PCR (sqRT-PCR)

For:

- transcript detection
- preliminary expression validation
- gel-based expression comparison


Design preference:

```
Amplicon:

200–500 bp

Preferred:

250–350 bp
```


---

## qRT-PCR

For:

- quantitative expression analysis
- RNA-seq validation
- relative expression studies


Design preference:

```
Amplicon:

80–200 bp

Preferred:

100–150 bp
```


---

## Homolog-specific primer design

For:

- duplicated genes
- gene families
- polyploid genomes


The workflow evaluates whether primers can distinguish:

- target gene
- paralogs
- homeologs
- allelic variants


---

# 3. Gene family and polyploid genome support


Many plant genomes contain highly similar duplicated genes.

This Skill performs family-level specificity evaluation.

Supported scenarios:

- autotetraploid crops
- allopolyploid crops
- duplicated gene families


Examples:

- Medicago sativa
- wheat
- cotton
- soybean
- other polyploid plants


The workflow evaluates:

```
Target gene

├── Allelic variants

├── Homeologs

├── Paralogs

└── Related family members
```


---

# 4. Primer specificity validation


Primer candidates are evaluated against:

Required:

```
CDS database

Transcript database
```


Recommended:

```
Genome sequence
```


The workflow classifies primers:


## Class I

Target-specific

Suitable for individual gene analysis


## Class II

Potential homolog amplification

Requires experimental consideration


## Class III

Non-specific

Not recommended


---

# 5. Experimental objective-aware workflow


Primer parameters are selected according to experimental purpose.

The Skill does not apply qPCR rules to cloning PCR or sqRT-PCR.


Supported objectives:

- full-length cloning
- expression detection
- quantitative expression analysis
- homolog discrimination
- candidate gene validation


---

# Recommended Workflow


```
Genome resources

        |

        ↓

Target gene identity validation

        |

        ↓

Transcript and sequence confirmation

        |

        ↓

Gene family / homolog analysis

        |

        ↓

Experimental objective selection

        |

        ↓

Primer design

        |

        ↓

Specificity validation

        |

        ↓

PCR simulation

        |

        ↓

Primer report

        |

        ↓

Experimental validation

        |

        ↓

Expression analysis
```


---

# Required Input


Users should provide:


```
1. Genome FASTA

2. GFF3 annotation

3. CDS FASTA

4. Transcript FASTA

5. Target information:

   - Gene ID
   - Gene list
   - Gene family name
   - Sequence

Optional:

   - Protein FASTA
   - Haplotype sequences
   - Reference homolog sequences
   - Previous primers
```


---

# Interactive Confirmation


Before execution, the Agent may request:


## Experimental purpose

Examples:

```
- cloning
- sqRT-PCR
- qRT-PCR
- homolog-specific detection
```


## Target source

Examples:

```
- genome annotation
- gene family analysis
- transcriptome analysis
- published gene
```


## Confidence of target identity

Examples:

```
confirmed
uncertain
unknown
```


This prevents automated design based on incorrect assumptions.


---

# Output Files


The Skill generates:


## Target validation


```
Gene_identity_validation.tsv

Gene_relationship_analysis.tsv
```


Including:

- annotation evidence
- transcript selection
- homolog information
- confidence level


---

## Primer results


```
Cloning_primers.tsv

sqRT_PCR_primers.tsv

qRT_PCR_primers.tsv
```


Including:

- primer sequence
- Tm
- GC content
- product size
- specificity evaluation


---

## Validation results


```
Primer_specificity.tsv

PCR_simulation.tsv
```


Including:

- off-target analysis
- expected amplification products
- primer reliability


---

## Final report


```
Primer_Design_Report.pdf
```


Including:

- target validation
- sequence evidence
- primer design strategy
- specificity assessment
- experimental recommendations


---

# Recommended Use Cases


Suitable for:


## Functional genomics

- candidate gene validation
- overexpression experiments
- gene characterization


## Gene family studies

- duplicated gene analysis
- homolog-specific expression


## Transcriptome studies

- RNA-seq validation
- stress response experiments
- developmental expression analysis


## Polyploid crop research

Examples:

- Medicago sativa
- wheat
- cotton
- soybean


---

# Limitations


This Skill does not:


- replace biological interpretation
- guarantee experimental PCR success
- determine gene function automatically
- replace laboratory validation


Experimental confirmation is still required:


For cloning:

- colony PCR
- Sanger sequencing
- ORF verification


For expression analysis:

- PCR validation
- qRT-PCR efficiency testing
- biological replication


---

# Version


Version:

```
2.0
```


Designed for:

Plant molecular biology laboratories


Workflow type:

```
AI-assisted genome-guided molecular biology pipeline
```


Core principle:

```
Validate the gene first.
Design primers second.
```
