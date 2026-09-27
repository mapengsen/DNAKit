---
inclusion: manual
---

# Bioinformatics File Formats

A practical guide to the most common file formats in life sciences data analysis.

## Sequence Formats

### FASTA (.fasta, .fa, .fna, .faa)
- **Description**: Plain-text format for nucleotide or amino acid sequences with a header line starting with `>`.
- **When to use**: Storing reference genomes, protein sequences, query sequences for BLAST.
- **Common tools**: BLAST, Clustal Omega, MUSCLE, samtools faidx.
- **Example**:
```
>sp|P04637|P53_HUMAN Cellular tumor antigen p53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
```

### FASTQ (.fastq, .fq, .fastq.gz)
- **Description**: Sequence format that includes per-base quality scores (Phred+33 encoding).
- **When to use**: Raw sequencing reads from Illumina, PacBio, or Nanopore instruments.
- **Common tools**: FastQC, Trimmomatic, cutadapt, BWA, STAR.
- **Example**:
```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

## Alignment Formats

### BAM / SAM (.bam, .sam, .cram)
- **Description**: Binary (BAM) or text (SAM) format for aligned sequencing reads. CRAM is a compressed alternative.
- **When to use**: Storing read alignments after mapping to a reference genome.
- **Common tools**: samtools, Picard, GATK, IGV, mosdepth.
- **Key fields**: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQ, QUAL.

## Variant Formats

### VCF (.vcf, .vcf.gz)
- **Description**: Variant Call Format for storing SNPs, indels, and structural variants.
- **When to use**: Variant calling output, population genetics, clinical annotation.
- **Common tools**: bcftools, GATK, VEP, SnpEff, ClinVar.
- **Example header**: `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE`

## Annotation Formats

### GFF / GTF (.gff, .gff3, .gtf)
- **Description**: General Feature Format for genomic annotations (genes, exons, CDS).
- **When to use**: Gene annotation, RNA-seq quantification, genome browsers.
- **Common tools**: bedtools, featureCounts, StringTie, gffread.

### BED (.bed)
- **Description**: Tab-delimited format for genomic intervals (chrom, start, end).
- **When to use**: Defining regions of interest, peak calling, coverage analysis.
- **Common tools**: bedtools, UCSC Genome Browser, deepTools.

## Imaging Formats

### DICOM (.dcm)
- **Description**: Digital Imaging and Communications in Medicine standard for medical images.
- **When to use**: Radiology, pathology, and clinical imaging data exchange.
- **Common tools**: pydicom, 3D Slicer, OsiriX, OHIF Viewer.

## Systems Biology Formats

### SBML (.sbml, .xml)
- **Description**: Systems Biology Markup Language for computational models of biological processes.
- **When to use**: Metabolic network modeling, kinetic simulations, pathway analysis.
- **Common tools**: COPASI, libSBML, CellDesigner, tellurium.

### BioPAX (.owl, .rdf)
- **Description**: Biological Pathway Exchange format for pathway data in RDF/OWL.
- **When to use**: Pathway data integration, knowledge graphs, semantic web applications.
- **Common tools**: Paxtools, PathVisio, Cytoscape.

## Experimental Metadata Formats

### MAGE-TAB (.idf.txt, .sdrf.txt)
- **Description**: MicroArray Gene Expression Tabular format for experiment metadata.
- **When to use**: Submitting microarray or RNA-seq experiments to ArrayExpress.
- **Common tools**: ArrayExpress submission tool, Annotare.

### ISA-Tab (.txt)
- **Description**: Investigation/Study/Assay framework for describing multi-omics experiments.
- **When to use**: Multi-assay experiment descriptions, data sharing, FAIR compliance.
- **Common tools**: ISAcreator, isatools Python library.

## Quick Reference

| Format | Extension | Binary? | Indexed? | Typical Size |
|--------|-----------|---------|----------|-------------|
| FASTA  | .fa       | No      | .fai     | MB–GB       |
| FASTQ  | .fq.gz    | Compressed | No   | GB–TB       |
| BAM    | .bam      | Yes     | .bai     | GB–TB       |
| VCF    | .vcf.gz   | Compressed | .tbi | MB–GB       |
| GFF3   | .gff3     | No      | .tbi     | MB          |
| BED    | .bed      | No      | No       | KB–MB       |
