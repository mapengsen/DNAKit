# varcode: variant effects, fusions and splice outcomes

This group of functions in `dnakit.molecular_lab` is a set of adapters over [varcode](https://github.com/openvax/varcode) (Apache-2.0), which predicts the effect of variants on transcripts and proteins, and [pyensembl](https://github.com/openvax/pyensembl) (Apache-2.0), the transcript-model database that varcode builds on. [pysam](https://github.com/pysam-developers/pysam) reads the genome FASTA. The functions predict the effect of small variants and structural variants on a transcript, apply germline and somatic variants to a transcript to get the mutant cDNA and protein, enumerate and realise splice plans and splice outcomes, import observed RNA fusion models (including ExACTO tables), list the genes that overlap variants, load variants from VCF, MAF and JSON files, check whether VCF samples belong to the same donors, and group structural-variant calls. They do not call variants, align reads or measure expression.

Most functions take a local annotation (`gtf`) and a genome FASTA (`genome_fasta`) and build a pyensembl reference from them with `build_local_reference()`, so nothing is downloaded. The reference is indexed once per content digest and cached under `cache_dir` (default `~/.cache/dnakit/varcode_reference`); the examples below use `cache_dir="cache"` in the working folder. Variants are mappings with `contig`, 1-based `start`, `ref` and `alt` (letters A, C, G, T, N; an indel leaves one of `ref` and `alt` empty). Structural variants are mappings with `contig`, 1-based inclusive `start` and `end`, and `sv_type` (`DEL`, `DUP`, `INV`, `INS` or `CNV`). A transcript is named by its `transcript_id`. Every result is a JSON-serialisable dictionary with `operation`, `status`, `backend`, `package` and `coordinates`; effect objects are summarised by their class name and description, plus the variant, the transcript ID and the gene name. Invalid arguments raise `ConfigurationError`; a missing package raises `BackendUnavailableError` with the install command.

The functions are imported with `from dnakit.molecular_lab import ...`. The dependencies are optional and are imported only when a function is called: the varcode functions need the `dnakit[bio-varcode]` extra, the local-reference functions also `dnakit[bio-pyensembl]` and `dnakit[bio-reads]` (pysam). varcode and pyensembl are Apache-2.0 licensed according to their package metadata, so there is no copyleft caveat for these two; DNAKit copies none of their code.

## Example data

The examples share one small synthetic genome: a single 450 bp contig named `1` (the name is `1`, not `chr1`, because the varcode VCF loader strips a `chr` prefix, so VCF-based examples only match a reference whose contigs have no prefix). It carries two plus-strand genes. `GENEA` (transcript `T1`) has exons 1-60, 101-150 and 191-240 and a coding sequence from position 10 to the stop codon at 228-230, which encodes a 46-residue protein starting `MPTVYG`. `GENEB` (transcript `T2`) has exons 281-330 and 371-420 and a coding sequence from 291 to the stop codon at 409-411. The donor of the second intron of `T1` (positions 151-152) is `AT` instead of the canonical `GT`, so splice planning sees a lost donor there. The variants used below are `1:13 C>G` (`p.P2A`), `1:20 TC>T` (a frameshift) and `1:61 G>A` (the donor of the first intron). The files are synthetic and were generated for this page, not downloaded; run the examples in a folder that contains all of them.

**Example file:** `genome.fa`

```text
>1
CCGACTGCAATGCCAACAGTCTATGGAGCTGTACATTACTATCGCATACGTTTAGAAGCC
GTTTGAGGGTGAGAGCCCATTCTAGACTTTGAGAGTCGAGGGTCCTGAGAAATGGAAAAA
AGCTGTGACTTTCTCTTATGACTACACAGAATTTAGCGCCCACCCCAGAAGGGGTTGCTC
TAAGCAACAGGCTCCGCAAACCCGATCACACGACCAAGACGAGACACTAAAAGTTTATAG
ATCCTAACATTCTAGAGCCGATAGAGCTAGCTTGACTGTCGACGCCATTAATGAGCCGAC
TGCACCCTACCGTAACAGGCAGACCCGCTCGTGTCTGAGAACTCCTCCTGCGACGGCAAG
CGAGCCGAAGGCCGGACTAAAGAGGTTCGTGCTTTAGAGATACAGGTCTGACGTTCGCTA
AGAACCCTAGTTCAAAGCATCGCAGTAGCT
```

**Example file:** `genes.gtf`

```text
1	demo	gene	1	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding";
1	demo	transcript	1	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	exon	1	60	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	exon	101	150	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	exon	191	240	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "3";
1	demo	CDS	10	60	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	CDS	101	150	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	CDS	191	227	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding"; exon_number "3";
1	demo	start_codon	10	12	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	stop_codon	228	230	.	+	.	gene_id "G1"; gene_name "GENEA"; gene_biotype "protein_coding"; transcript_id "T1"; transcript_biotype "protein_coding";
1	demo	gene	281	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding";
1	demo	transcript	281	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
1	demo	exon	281	330	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	exon	371	420	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	CDS	291	330	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "1";
1	demo	CDS	371	408	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding"; exon_number "2";
1	demo	start_codon	291	293	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
1	demo	stop_codon	409	411	.	+	.	gene_id "G2"; gene_name "GENEB"; gene_biotype "protein_coding"; transcript_id "T2"; transcript_biotype "protein_coding";
```

The attributes are the minimum that varcode and pyensembl need here: `gene_id`, `gene_name`, `gene_biotype`, `transcript_id`, `transcript_biotype` and `exon_number`. `start_codon` and `stop_codon` lines mark the transcripts as complete; the CDS lines do not include the stop codon.

**Example file:** `phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=450>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	13	.	C	G	.	PASS	.	GT:PS	0|1:13
1	17	.	C	T	.	PASS	.	GT:PS	0|1:13
```

A phased VCF with sample `S1`: both variants are on the second haplotype (`0|1`) of phase set 13.

**Example file:** `structures.tsv`

```text
transcript_model_id	reference_transcript_ids	index	read_start	read_end	sequence	type	kind	context	chromosome_1	position_1	strand_1	chromosome_2	position_2	strand_2	transcript_id_1	transcript_id_2
M1	T1,T2	0	0	59	CCGACTGCAATGCCAACAGTCTATGGAGCTGTACATTACTATCGCATACGTTTAGAAGCC	base	match	exon	1	1	+	1	60	+	T1	T1
M1	T1,T2	1	60	109	GGTCCTGAGAAATGGAAAAAAGCTGTGACTTTCTCTTATGACTACACAGA	base	match	exon	1	101	+	1	150	+	T1	T1
M1	T1,T2	2	110	159	GACGCCATTAATGAGCCGACTGCACCCTACCGTAACAGGCAGACCCGCTC	base	match	exon	1	281	+	1	330	+	T2	T2
```

An ExACTO transcript-structure table with one model, `M1`, built from `T1` exons 1 and 2 (positions 1-60 and 101-150) joined to `T2` exon 1 (281-330). The columns are tab separated.

**Example file:** `integrated.tsv`

```text
transcript_model_id	reference_transcript_ids	rna_variant_call_id	dna_variant_call_id
M1	T1,T2	R1	D1
```

The ExACTO integration table that links the RNA call `R1` of model `M1` to the DNA call `D1`.

**Example file:** `donor1.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D1
chr1	1000	.	A	G	.	PASS	.	GT	0/1
chr1	2000	.	C	T	.	PASS	.	GT	1/1
chr1	3000	.	G	A	.	PASS	.	GT	0/0
chr1	4000	.	T	C	.	PASS	.	GT	0/1
chr1	5000	.	A	C	.	PASS	.	GT	0/0
```

**Example file:** `donor2.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D2
chr1	1000	.	A	G	.	PASS	.	GT	0/1
chr1	2000	.	C	T	.	PASS	.	GT	1/1
chr1	3000	.	G	A	.	PASS	.	GT	0/0
chr1	4000	.	T	C	.	PASS	.	GT	0/1
chr1	5000	.	A	C	.	PASS	.	GT	0/0
```

**Example file:** `donor3.vcf`

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	D3
chr1	1000	.	A	G	.	PASS	.	GT	1/1
chr1	2000	.	C	T	.	PASS	.	GT	0/0
chr1	3000	.	G	A	.	PASS	.	GT	0/1
chr1	4000	.	T	C	.	PASS	.	GT	0/0
chr1	5000	.	A	C	.	PASS	.	GT	1/1
```

The three germline VCFs have five SNPs each. `donor1.vcf` and `donor2.vcf` carry the same genotypes; `donor3.vcf` carries different ones.

**Example file:** `calls.csv`

```text
call_id,caller,sample,build,chrom,pos,ref,alt,info
c1,callerA,S1,GRCh38,1,10000,N,<DEL>,END=10500;SVTYPE=DEL
c2,callerB,S1,GRCh38,1,10020,N,<DEL>,END=10490;SVTYPE=DEL
c3,callerA,S1,GRCh38,2,5000,N,<DEL>,END=5300;SVTYPE=DEL
```

Three structural-variant calls: `c1` and `c2` are deletions on contig 1 whose breakpoints differ by at most 20 bases, and `c3` is a deletion on contig 2.

## 1) `dnakit.molecular_lab.build_local_reference()`

- **Purpose:** Index a local GTF and genome FASTA as a pyensembl reference for varcode, without downloading.
- **Method:** Reads the GTF (every non-comment line needs nine tab-separated fields, and exon and CDS lines need a `transcript_id`; a missing CDS `protein_id` is set to the transcript ID and a missing `exon_id` is generated). With pysam it cuts the exons out of the genome FASTA, joins them on the transcript strand into a cDNA FASTA, and translates the CDS lines into a protein FASTA. It then builds a `pyensembl.Genome` from these files and calls `index()`. The index is stored under `cache_dir/<digest>`, where the digest is a SHA-256 of the GTF content plus the FASTA path, size and modification time. An index already built in this process is reused and its original folder is reported.
- **API:** `build_local_reference(gtf: 'str', genome_fasta: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing files; `cache_dir` is the index folder. Every transcript must lie on one contig and strand, its contig must be in the FASTA, and at least one transcript needs exon lines (otherwise `ConfigurationError`).
- **Output:** `reference_name` (`dnakit_local`), `annotation_name`, `cache_directory`, `reused` (true when the index was already loaded in this process) and the counts `contigs`, `genes`, `transcripts` and `exons`.
- **Example code:**

```python
from dnakit.molecular_lab import build_local_reference

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

first = build_local_reference(**ref)
again = build_local_reference(**ref)
print(first["reference_name"], first["contigs"], first["genes"], first["transcripts"], first["exons"])
print(first["reused"], again["reused"])
```

- **Example output:**

```text
dnakit_local 1 2 2 5
False True
```

## 2) `dnakit.molecular_lab.varcode_apply_germline_to_transcript()`

- **Purpose:** Apply germline variants to a transcript and return the resulting patient transcript.
- **Method:** Builds a varcode `GermlineContext` from the germline variants and calls `varcode.germline.apply_germline_to_transcript`. Without a somatic variant, every germline variant that lies in an exon of the transcript is applied. With a somatic variant, only the germline variants inside its window (the one reported by `varcode_default_germline_window()`) are applied. When no germline variant qualifies, or the edits conflict or fall outside the CDS, varcode returns nothing and the wrapper returns an empty mutant record (empty strings, `None` sequences, no edits) instead of raising.
- **API:** `varcode_apply_germline_to_transcript(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', germline_variants: 'Sequence[Mapping[str, Any]]', *, somatic_variant: 'Mapping[str, Any] | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `germline_variants` is a list of up to 1000 variant mappings; `somatic_variant` is an optional variant mapping.
- **Output:** `transcript_id` and `mutant_transcript`, a dictionary with `transcript_id`, `annotator`, `cdna_sequence`, `mutant_protein_sequence` and `edits` (each edit has `cdna_start`, `cdna_end`, `alt_bases` and `source_variant`; the cDNA offsets are 0-based and half-open).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_apply_germline_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

somatic = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
near = {"contig": "1", "start": 14, "ref": "C", "alt": "T"}
far = {"contig": "1", "start": 17, "ref": "C", "alt": "T"}

result = varcode_apply_germline_to_transcript(
    transcript_id="T1", germline_variants=[near], somatic_variant=somatic, **ref
)
mutant = result["mutant_transcript"]
print(mutant["annotator"], mutant["mutant_protein_sequence"][:6], len(mutant["edits"]))

outside = varcode_apply_germline_to_transcript(
    transcript_id="T1", germline_variants=[far], somatic_variant=somatic, **ref
)
print(outside["mutant_transcript"]["cdna_sequence"])
```

- **Example output:**

```text
protein_diff MLTVYG 1
None
```

## 3) `dnakit.molecular_lab.varcode_apply_variant_to_transcript()`

- **Purpose:** Apply one variant to a transcript and return the mutant cDNA and protein.
- **Method:** Calls `varcode.mutant_transcript.apply_variant_to_transcript` with the variant and the transcript model, and summarises the resulting `MutantTranscript`: the spliced mutant cDNA, the protein translated from the CDS start, and the edits that were applied.
- **API:** `varcode_apply_variant_to_transcript(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `variant`, `transcript_id` and `mutant_transcript` with `transcript_id`, `annotator`, `cdna_sequence`, `mutant_protein_sequence` and `edits`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_apply_variant_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variant = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
result = varcode_apply_variant_to_transcript(variant=variant, transcript_id="T1", **ref)
mutant = result["mutant_transcript"]
print(len(mutant["cdna_sequence"]), mutant["mutant_protein_sequence"][:8], len(mutant["mutant_protein_sequence"]))
print(mutant["edits"][0]["cdna_start"], mutant["edits"][0]["cdna_end"], mutant["edits"][0]["alt_bases"])
```

- **Example output:**

```text
160 MATVYGAV 46
12 13 G
```

## 4) `dnakit.molecular_lab.varcode_apply_variants_to_transcript()`

- **Purpose:** Apply several variants together to one transcript and return the resulting mutant.
- **Method:** Calls `varcode.mutant_transcript.apply_variants_to_transcript` once with all variants,, so the edits are combined in a single mutant cDNA and the protein is translated from that cDNA (in the example the frameshift truncates the protein after the substitution).
- **API:** `varcode_apply_variants_to_transcript(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variants` is a list of up to 1000 variant mappings. `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `transcript_id`, `variants` (the variants as varcode normalised them) and `mutant_transcript` (same fields as in `varcode_apply_variant_to_transcript()`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_apply_variants_to_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = [{"contig": "1", "start": 13, "ref": "C", "alt": "G"}, {"contig": "1", "start": 20, "ref": "TC", "alt": "T"}]
result = varcode_apply_variants_to_transcript(variants=variants, transcript_id="T1", **ref)
mutant = result["mutant_transcript"]
print([(v["start"], v["ref"], v["alt"]) for v in result["variants"]])
print(mutant["mutant_protein_sequence"], len(mutant["edits"]))
```

- **Example output:**

```text
[(13, 'C', 'G'), (21, 'C', '')]
MATVMELYITIAYV 2
```

## 5) `dnakit.molecular_lab.varcode_candidates_from_effects()`

- **Purpose:** Predict a variant's effects and wrap each one as an effect candidate with a source label.
- **Method:** Effect objects cannot be passed in as JSON, so the wrapper first predicts the variant's effects on every transcript it overlaps (`varcode.effects.predict_variant_effects`) and then wraps them with `varcode.candidates_from_effects(effects, source=source)`.
- **API:** `varcode_candidates_from_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', *, source: 'str' = 'varcode', cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `source` is a non-empty label (default `varcode`).
- **Output:** `variant`, `source`, `candidates` and `count`. Each candidate has `source`, `evidence` (a dictionary) and `effect` (`class`, `description`, `variant`, `transcript_id`, `gene_name`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_candidates_from_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_candidates_from_effects(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, source="demo", **ref)
print(result["count"], result["source"])
candidate = result["candidates"][0]
print(candidate["source"], candidate["evidence"], candidate["effect"]["class"], candidate["effect"]["description"])
```

- **Example output:**

```text
1 demo
demo {} Substitution p.P2A
```

## 6) `dnakit.molecular_lab.varcode_classify_from_protein_diff()`

- **Purpose:** Classify a variant's protein change by comparing a reference protein with a mutant protein.
- **Method:** Calls `varcode.effects.classify.classify_from_protein_diff` with the variant, the transcript, the two protein strings and `length_delta`. Identical proteins give `Silent`; a mutant that no longer starts with `M` gives `StartLoss`; otherwise the shared flanks are trimmed and the difference is classified (substitution, insertion, deletion, premature stop, stop loss). A `length_delta` that is not a multiple of 3 marks a frameshift when the tail of the protein changes.
- **API:** `varcode_classify_from_protein_diff(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', ref_protein: 'str', mut_protein: 'str', *, length_delta: 'int' = 0, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `ref_protein` and `mut_protein` are non-empty strings; `length_delta` is the net nucleotide length change of the edit (default 0).
- **Output:** `variant`, `transcript_id` and `effect` (`class`, `description`, `variant`, `transcript_id`, `gene_name`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_classify_from_protein_diff

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variant = {"contig": "1", "start": 13, "ref": "C", "alt": "G"}
same = varcode_classify_from_protein_diff(
    variant=variant, transcript_id="T1", ref_protein="MPTV", mut_protein="MPTV", **ref
)
print(same["effect"]["class"], same["effect"]["description"])
stop = varcode_classify_from_protein_diff(
    variant=variant, transcript_id="T1", ref_protein="MPTV", mut_protein="MP", length_delta=-6, **ref
)
print(stop["effect"]["class"], stop["effect"]["description"])
```

- **Example output:**

```text
Silent p.M1=
PrematureStop p.T3*
```

## 7) `dnakit.molecular_lab.varcode_codon_table_for_transcript()`

- **Purpose:** Select the codon table that varcode uses for a transcript.
- **Method:** Calls `varcode.effects.codon_tables.codon_table_for_transcript` (exposed as `varcode.annotators.protein_diff.codon_table_for_transcript`). Transcripts on a contig whose name is one of `mt`, `chrmt`, `m`, `chrm`, `mtdna` or `mito` (case-insensitive) get the vertebrate mitochondrial table (NCBI table 2); every other transcript gets the standard table.
- **API:** `varcode_codon_table_for_transcript(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `transcript_id`, `contig` and `codon_table` with `name`, `ncbi_id`, and the sorted `start_codons` and `stop_codons`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_codon_table_for_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_codon_table_for_transcript(transcript_id="T1", **ref)
table = result["codon_table"]
print(result["contig"], table["name"], table["ncbi_id"], table["stop_codons"])
```

- **Example output:**

```text
1 Standard 1 ['TAA', 'TAG', 'TGA']
```

## 8) `dnakit.molecular_lab.varcode_default_germline_window()`

- **Purpose:** Return the genomic window that is searched for germline variants that can phase with a somatic variant.
- **Method:** Calls `varcode.germline.default_germline_window` with the somatic variant and the transcript. The window runs from three bases before to three bases after the somatic variant; when an exon boundary of the transcript lies within 6 bases of the variant, the window is widened to 6 bases on both sides of that boundary so that donor and acceptor signals are covered.
- **API:** `varcode_default_germline_window(gtf: 'str', genome_fasta: 'str', somatic_variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `somatic_variant` is a variant mapping. `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `contig`, `start` and `end` of the window (1-based inclusive).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_default_germline_window

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

inside = varcode_default_germline_window(somatic_variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(inside["contig"], inside["start"], inside["end"])
near_edge = varcode_default_germline_window(
    somatic_variant={"contig": "1", "start": 58, "ref": "C", "alt": "G"}, transcript_id="T1", **ref
)
print(near_edge["contig"], near_edge["start"], near_edge["end"])
```

- **Example output:**

```text
1 10 16
1 54 66
```

## 9) `dnakit.molecular_lab.varcode_infer_genome()`

- **Purpose:** Resolve an Ensembl release from a release number or a reference name.
- **Method:** Calls `varcode.reference.infer_genome`. An integer is an Ensembl release number for human; a string is a reference name such as `GRCh38` or `GRCh38:93`. UCSC names such as `hg19` are converted to the Ensembl equivalent (GRCh37). The release object is created lazily, so the call itself needs no data files.
- **API:** `varcode_infer_genome(reference: 'str | int', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package varcode (extra `dnakit[bio-varcode]`). It works offline. Resolving the name downloads nothing; downloading the release data is a separate step (see `varcode_install_reference_data()`).
- **Input:** `reference` is an integer or a non-empty string. A boolean, a float or an unknown name raises `ConfigurationError`.
- **Output:** `requested`, `reference_name`, `release`, `species` (Latin name) and `converted_from_ucsc`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_infer_genome

ucsc = varcode_infer_genome("hg19", verbose=False)
print(ucsc["reference_name"], ucsc["release"], ucsc["species"], ucsc["converted_from_ucsc"])
number = varcode_infer_genome(93, verbose=False)
print(number["reference_name"], number["release"], number["converted_from_ucsc"])
```

- **Example output:**

```text
GRCh37 75 homo_sapiens True
GRCh38 93 False
```

## 10) `dnakit.molecular_lab.varcode_infer_genome_for_reference_name()`

- **Purpose:** Resolve a reference name, optionally with a release, to an Ensembl release.
- **Method:** Calls `varcode.reference.infer_genome_for_reference_name`. The name may be a plain reference name (`GRCh37`, which maps to the newest release known for it), a name with a release (`GRCh38:93`) or a UCSC name (`hg38`), which is converted. Nothing is downloaded.
- **API:** `varcode_infer_genome_for_reference_name(reference_name: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package varcode (extra `dnakit[bio-varcode]`). It works offline.
- **Input:** `reference_name` is a non-empty string.
- **Output:** `requested`, `reference_name`, `release`, `species` and `converted_from_ucsc`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_infer_genome_for_reference_name

with_release = varcode_infer_genome_for_reference_name("GRCh38:93", verbose=False)
print(with_release["requested"], with_release["reference_name"], with_release["release"])
older = varcode_infer_genome_for_reference_name("GRCh37", verbose=False)
print(older["reference_name"], older["release"], older["converted_from_ucsc"])
```

- **Example output:**

```text
GRCh38:93 GRCh38 93
GRCh37 75 False
```

## 11) `dnakit.molecular_lab.varcode_left_align_indels()`

- **Purpose:** Shift indels to their leftmost equivalent position using the local reference.
- **Method:** Builds a `VariantCollection` on the local reference and calls `varcode.transforms.left_align_indels`, which slides each insertion or deletion to the left through repeated sequence. Substitutions are not moved.
- **API:** `varcode_left_align_indels(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variants` is a list of up to 1000 variant mappings.
- **Output:** `variants` (the aligned variants with `contig`, `start`, `ref`, `alt`) and `count`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_left_align_indels

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

# The reference has CC at positions 13-14, so deleting the second C is the same as deleting the first.
variants = [{"contig": "1", "start": 14, "ref": "C", "alt": ""}]
result = varcode_left_align_indels(variants=variants, **ref)
print(result["count"], result["variants"])
```

- **Example output:**

```text
1 [{'contig': '1', 'start': 13, 'ref': 'C', 'alt': ''}]
```

## 12) `dnakit.molecular_lab.varcode_make_rna_outcome()`

- **Purpose:** Wrap a variant's effect on one transcript as an RNA-sourced effect candidate.
- **Method:** Predicts the effect of the variant on the transcript (`predict_variant_effect_on_transcript`) and passes it to `varcode.make_rna_outcome` with the source label, the optional transcript model ID and the optional read count. No RNA reads are read; the read count is only stored as evidence and is not used as a likelihood.
- **API:** `varcode_make_rna_outcome(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, source: 'str' = 'rna', transcript_model_id: 'str | None' = None, read_count: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `source` is a non-empty label (default `rna`), `transcript_model_id` an optional identifier and `read_count` an optional non-negative integer.
- **Output:** `variant`, `transcript_id` and `candidate` with `source`, `evidence` (the model ID and read count when given) and `effect`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_make_rna_outcome

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_make_rna_outcome(
    variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", source="rna", transcript_model_id="m1", read_count=12, **ref
)
candidate = result["candidate"]
print(candidate["source"], candidate["evidence"])
print(candidate["effect"]["class"], candidate["effect"]["description"])
```

- **Example output:**

```text
rna {'transcript_model_id': 'm1', 'read_count': 12}
Substitution p.P2A
```

## 13) `dnakit.molecular_lab.varcode_predict_germline_aware_effect()`

- **Purpose:** Predict a somatic variant's effect against the germline-applied transcript.
- **Method:** Builds a `GermlineContext` from the germline variants and calls `varcode.germline.predict_germline_aware_effect` with the somatic variant, the transcript, the context and the annotator (the default annotator when `annotator` is None). When a germline variant shares the window with the somatic variant and the phase between them is unknown, varcode returns a `PhaseCandidateSet`, the set of plausible effects across the phase hypotheses (its short description starts with `?`). The wrapper reports only the class and the short description of that effect, so the individual hypotheses are not listed.
- **API:** `varcode_predict_germline_aware_effect(gtf: 'str', genome_fasta: 'str', somatic_variant: 'Mapping[str, Any]', transcript_id: 'str', germline_variants: 'Sequence[Mapping[str, Any]]', *, annotator: 'str | None' = None, max_hypotheses: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `somatic_variant` is a variant mapping; `germline_variants` a list of variant mappings; `annotator` is None or one of `fast`, `protein_diff`, `transcript_model`, `realized`; `max_hypotheses` is None or an integer from 1 to 10000. `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `transcript_id`, `somatic_variant`, `germline_count` and `effect` (`class`, `description`, `variant`, `transcript_id`, `gene_name`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_germline_aware_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_predict_germline_aware_effect(
    somatic_variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", germline_variants=[{"contig": "1", "start": 14, "ref": "C", "alt": "T"}], **ref
)
print(result["germline_count"], result["effect"]["class"], result["effect"]["description"])
```

- **Example output:**

```text
1 PhaseCandidateSet ?p.P2A
```

## 14) `dnakit.molecular_lab.varcode_predict_realized_effect()`

- **Purpose:** Predict the realised top effect of somatic variants on one transcript, with the alternative outcomes.
- **Method:** Calls `varcode.realized_effects.predict_realized_effect` with the somatic variants (applied together), the transcript, the germline variants and `max_hypotheses`. Germline variants are phased against the somatic ones when varcode can resolve phase; the alternative outcomes it keeps (up to `max_hypotheses`) are returned as `candidates`, and several of them can carry the same description. The wrapper reports the top effect and at most 1000 alternatives.
- **API:** `varcode_predict_realized_effect(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', transcript_id: 'str', *, germline_variants: 'Sequence[Mapping[str, Any]]' = (), max_hypotheses: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variants` is a list of up to 1000 variant mappings; `germline_variants` a list of variant mappings (default none); `max_hypotheses` an integer from 1 to 10000 (default 64). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `transcript_id`, `top_effect` (`class`, `description`, ...) and `candidates`, a list of alternative outcomes with `source`, `evidence` and `effect`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_realized_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

alone = varcode_predict_realized_effect(variants=[{"contig": "1", "start": 13, "ref": "C", "alt": "G"}], transcript_id="T1", **ref)
print(alone["top_effect"]["description"], len(alone["candidates"]))

phased = varcode_predict_realized_effect(
    variants=[{"contig": "1", "start": 13, "ref": "C", "alt": "G"}], transcript_id="T1", germline_variants=[{"contig": "1", "start": 14, "ref": "C", "alt": "T"}], **ref
)
print(phased["top_effect"]["description"], len(phased["candidates"]))
print(sorted({c["effect"]["description"] for c in phased["candidates"]}))
```

- **Example output:**

```text
p.P2A 3
p.P2A 6
['p.L2V', 'p.P2A']
```

## 15) `dnakit.molecular_lab.varcode_predict_variant_effect_on_transcript()`

- **Purpose:** Return the effect (for example a frameshift) of one variant on one transcript.
- **Method:** Calls `varcode.effects.predict_variant_effect_on_transcript` with the variant, the transcript model and the selected annotator. A reference mismatch or another upstream `ValueError` becomes a `ConfigurationError`.
- **API:** `varcode_predict_variant_effect_on_transcript(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, annotator: 'str | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `annotator` is None or one of `fast`, `protein_diff`, `transcript_model`, `realized`.
- **Output:** `variant` (as varcode normalised it), `transcript_id` and `effect` (`class`, `description`, `variant`, `transcript_id`, `gene_name`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_variant_effect_on_transcript

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_predict_variant_effect_on_transcript(variant={"contig": "1", "start": 20, "ref": "TC", "alt": "T"}, transcript_id="T1", **ref)
print(result["variant"])
print(result["effect"]["class"], result["effect"]["description"])
```

- **Example output:**

```text
{'contig': '1', 'start': 21, 'ref': 'C', 'alt': ''}
FrameShift p.Y5fs
```

## 16) `dnakit.molecular_lab.varcode_predict_variant_effect_on_transcript_or_failure()`

- **Purpose:** Predict a variant's effect on one transcript and report a failure instead of raising.
- **Method:** Calls `varcode.effects.predict_variant_effect_on_transcript_or_failure`. When varcode returns an exception or a `Failure` or `Unresolved` effect (for example because the `ref` allele does not match the reference), the result has `status` set to `failed`, `effect` set to None and the failure text in `failure`; otherwise `status` is `completed`. Invalid arguments still raise `ConfigurationError`.
- **API:** `varcode_predict_variant_effect_on_transcript_or_failure(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `variant`, `transcript_id`, `status` (`completed` or `failed`), `effect` (None on failure) and `failure` (the upstream failure text, or None).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_variant_effect_on_transcript_or_failure

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

good = varcode_predict_variant_effect_on_transcript_or_failure(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(good["status"], good["effect"]["description"], good["failure"])

wrong_ref = {"contig": "1", "start": 13, "ref": "A", "alt": "G"}
bad = varcode_predict_variant_effect_on_transcript_or_failure(variant=wrong_ref, transcript_id="T1", **ref)
print(bad["status"], bad["effect"], "Reference allele mismatch" in bad["failure"])
```

- **Example output:**

```text
completed p.P2A None
failed None True
```

## 17) `dnakit.molecular_lab.varcode_predict_variant_effects()`

- **Purpose:** Predict the effect of one variant on every transcript it overlaps.
- **Method:** Builds a varcode `Variant` on the local reference and calls `varcode.effects.predict_variant_effects`. A variant that overlaps no gene gets an `Intergenic` effect. With the default `raise_on_error=False`, an annotation error (for example a wrong `ref` allele) is returned as a `Failure` effect; with `raise_on_error=True` it raises `ConfigurationError`.
- **API:** `varcode_predict_variant_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', *, raise_on_error: 'bool' = False, annotator: 'str | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `raise_on_error` is a boolean (default False); `annotator` is None or one of `fast`, `protein_diff`, `transcript_model`, `realized`.
- **Output:** `variant`, `effects` (one record per effect with `class`, `description`, and `transcript_id` and `gene_name` when a transcript applies) and `count`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_variant_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = {
    "missense": {"contig": "1", "start": 13, "ref": "C", "alt": "G"},
    "splice donor": {"contig": "1", "start": 61, "ref": "G", "alt": "A"},
    "between genes": {"contig": "1", "start": 255, "ref": "G", "alt": "A"},
}
for label, variant in variants.items():
    result = varcode_predict_variant_effects(variant=variant, **ref)
    print(label, result["count"], [(e["class"], e["description"]) for e in result["effects"]])
```

- **Example output:**

```text
missense 1 [('Substitution', 'p.P2A')]
splice donor 1 [('SpliceOutcomeSet', 'splice-set:intron-retention (donor side, predicted)')]
between genes 1 [('Intergenic', 'intergenic')]
```

## 18) `dnakit.molecular_lab.varcode_apply_phase_resolver_to_effects()`

- **Purpose:** Predict a variant's effects and resolve them with the phase recorded in a VCF sample.
- **Method:** Loads `phased_vcf` with `varcode.load_vcf` on the local reference, builds a `VCFPhaseResolver` for the sample, predicts the variant's effects and passes them to `varcode.apply_phase_resolver_to_effects`. In varcode 10.11 the `VCFPhaseResolver` only answers cis and trans questions and has no mutant-transcript method, so `apply_phase_resolver_to_effects` leaves the effects unchanged: the result equals a plain effect prediction, and the call mainly checks that the phased VCF and the sample can be loaded.
- **API:** `varcode_apply_phase_resolver_to_effects(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', phased_vcf: 'str', sample: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `phased_vcf` is a local VCF whose `sample` column has phased genotypes (`0|1`) and `PS` values; `sample` is the column name. The contig names in the VCF must match the reference after varcode drops a `chr` prefix.
- **Output:** `variant`, `sample`, `effects` (records with `class`, `description`, `transcript_id`, `gene_name`) and `count`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_apply_phase_resolver_to_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_apply_phase_resolver_to_effects(
    variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, phased_vcf="phased.vcf", sample="S1", **ref
)
print(result["sample"], result["count"], [(e["class"], e["description"]) for e in result["effects"]])
```

- **Example output:**

```text
S1 1 [('Substitution', 'p.P2A')]
```

## 19) `dnakit.molecular_lab.varcode_build_exon_runs()`

- **Purpose:** Project a transcript's sense-strand exons onto its genomic layout.
- **Method:** Builds a `GenomicLayout` for the transcript with `flank` bases on each side (`GenomicLayout.from_transcript`, with the sequence read from the local genome FASTA) and calls `varcode.transcript_layout.build_exon_runs`. Each run is one occurrence of a reference exon in the layout.
- **API:** `varcode_build_exon_runs(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, flank: 'int' = 50, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `flank` is an integer from 0 to 100000 (default 50); the flank is clipped at the start of the contig.
- **Output:** `transcript_id`, `flank` and `runs`. Each run has `key` (exon number and occurrence), `exon_number`, `occurrence`, `layout_start` and `layout_end` (0-based offsets in the layout), `reference_start`, `reference_end`, `first_genomic_position` and `last_genomic_position`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_build_exon_runs

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_build_exon_runs(transcript_id="T2", flank=10, **ref)
for run in result["runs"]:
    print(run["key"], run["layout_start"], run["layout_end"], run["first_genomic_position"], run["last_genomic_position"])
```

- **Example output:**

```text
[1, 1] 10 60 281 330
[2, 1] 100 150 371 420
```

## 20) `dnakit.molecular_lab.varcode_enumerate_from_structural_variant()`

- **Purpose:** Enumerate cryptic-exon candidates in the flanks of a structural variant breakpoint.
- **Method:** Builds a varcode `StructuralVariant` on the local reference and calls `varcode.cryptic_exons.enumerate_from_structural_variant` with the flank size. varcode scans the sequence around the breakpoints with its default splice-site scorer and pairs acceptor-like and donor-like motifs into candidate exons. The wrapper converts each candidate with `str()`, so a candidate is reported as the text of its class and variant; the interval and the donor and acceptor scores that the varcode object holds are not in the result. varcode prints a warning once per genome when a breakpoint has no reference sequence available.
- **API:** `varcode_enumerate_from_structural_variant(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', *, flank: 'int' = 500, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structural_variant` has `contig`, `start`, `end` (`end` at least `start`) and `sv_type` (`DEL`, `DUP`, `INV`, `INS` or `CNV`, case-insensitive); `flank` is an integer from 0 to 100000 (default 500).
- **Output:** `structural_variant`, `flank`, `candidates` (strings) and `count`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_enumerate_from_structural_variant

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

sv = {"contig": "1", "start": 70, "end": 90, "sv_type": "DEL"}
narrow = varcode_enumerate_from_structural_variant(structural_variant=sv, flank=30, **ref)
wide = varcode_enumerate_from_structural_variant(structural_variant=sv, flank=100, **ref)
print(narrow["count"], wide["count"])
print(wide["candidates"][0])
```

- **Example output:**

```text
0 1
CrypticExonCandidate(variant=StructuralVariant(contig='1', start=70, end=90, sv_type='DEL', alt='<DEL>', reference_name='dnakit_local'))
```

## 21) `dnakit.molecular_lab.varcode_enumerate_layout_splice_plans()`

- **Purpose:** Enumerate the graph-valid splice plans of a transcript layout.
- **Method:** Builds the genomic layout of the transcript and calls `varcode.splice_graph.enumerate_layout_splice_plans`. Each donor or acceptor of the reference layout that is not intact (here judged from the sequence of the layout) gives a splice choice with the mechanisms exon skip, intron retention, cryptic site and normal; each plan fixes one choice per disrupted site and lists the exon runs it keeps. An intact layout has one plan.
- **API:** `varcode_enumerate_layout_splice_plans(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, flank: 'int' = 50, max_plans: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `flank` is an integer from 0 to 100000 (default 50); `max_plans` is an integer from 1 to 10000 (default 64).
- **Output:** `transcript_id`, `flank`, `plans` and `count`. Each plan has `choices` (pairs of a site and a mechanism with its `ordinal_rank`, `skipped_runs` and `required_runs`) and `kept_runs`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_enumerate_layout_splice_plans

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

broken = varcode_enumerate_layout_splice_plans(transcript_id="T1", flank=10, **ref)
print(broken["count"])
for plan in broken["plans"]:
    site, choice = plan["choices"][0]
    print(site["side"], site["exon_number"], choice["mechanism"], plan["kept_runs"])

intact = varcode_enumerate_layout_splice_plans(transcript_id="T2", flank=10, **ref)
print(intact["count"], intact["plans"][0]["choices"])
```

- **Example output:**

```text
4
donor 2 exon_skip [[1, 1], [3, 1]]
donor 2 intron_retention [[1, 1], [2, 1], [3, 1]]
donor 2 cryptic [[1, 1], [2, 1], [3, 1]]
donor 2 normal [[1, 1], [2, 1], [3, 1]]
1 []
```

## 22) `dnakit.molecular_lab.varcode_predict_structural_variant_effect()`

- **Purpose:** Classify a structural variant (deletion, duplication, inversion, insertion or CNV) on one transcript.
- **Method:** Builds a varcode `StructuralVariant` on the local reference and calls `varcode.effects.structural.predict_structural_variant_effect` with the transcript.
- **API:** `varcode_predict_structural_variant_effect(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structural_variant` has `contig`, `start`, `end` (`end` at least `start`) and `sv_type` (`DEL`, `DUP`, `INV`, `INS` or `CNV`, case-insensitive). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `structural_variant` (as parsed), `transcript_id` and `effect` (`class`, `description`, `variant`, `transcript_id`, `gene_name`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_predict_structural_variant_effect

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

cases = [
    ("T1", {"contig": "1", "start": 101, "end": 150, "sv_type": "DEL"}),
    ("T1", {"contig": "1", "start": 101, "end": 150, "sv_type": "DUP"}),
    ("T1", {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}),
]
for transcript_id, sv in cases:
    result = varcode_predict_structural_variant_effect(structural_variant=sv, transcript_id=transcript_id, **ref)
    print(sv["sv_type"], sv["start"], sv["end"], result["effect"]["class"], result["effect"]["description"])
```

- **Example output:**

```text
DEL 101 150 FrameShift p.G18fs
DUP 101 150 FrameShift p.L35fs
DEL 151 280 GeneFusion sv-gene-fusion
```

## 23) `dnakit.molecular_lab.varcode_random_variants()`

- **Purpose:** Draw random variants that overlap complete coding transcripts of the local reference.
- **Method:** Calls `varcode.util.random_variants` with the local genome as `ensembl`. It samples positions in complete coding transcripts (those with a start codon, a stop codon and a protein) and builds substitutions and, if enabled, insertions and deletions from the reference. The positions and reference bases follow `random_seed`. The alternative base, however, is drawn by varcode from a Python `set`, whose order changes between Python processes (string hashing is randomised), so the same seed can give different `alt` values in different processes; within one process the draws repeat.
- **API:** `varcode_random_variants(gtf: 'str', genome_fasta: 'str', count: 'int', *, deletions: 'bool' = True, insertions: 'bool' = True, random_seed: 'int | None' = None, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `count` is an integer from 1 to 10000; `deletions` and `insertions` are booleans (default True); `random_seed` is None or an integer from 0 to 2**32 - 1.
- **Output:** `count`, `random_seed` and `variants` (records with `contig`, `start`, `ref`, `alt`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_random_variants

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

first = varcode_random_variants(count=4, deletions=False, insertions=False, random_seed=5, **ref)
second = varcode_random_variants(count=4, deletions=False, insertions=False, random_seed=5, **ref)
print(first["count"], first["variants"] == second["variants"])
for variant in first["variants"]:
    print(variant["contig"], variant["start"], variant["ref"], len(variant["alt"]), variant["alt"] != variant["ref"])
```

- **Example output:**

```text
4 True
1 8 C 1 True
1 150 A 1 True
1 305 C 1 True
1 415 T 1 True
```

## 24) `dnakit.molecular_lab.varcode_realize_exon_path()`

- **Purpose:** Splice selected exon runs of a transcript and translate from its mapped CDS start.
- **Method:** Builds the genomic layout and calls `varcode.transcript_layout.realize_exon_path`, which joins the selected exon runs into a cDNA and translates it from the CDS start. With `kept_run_keys` None all runs are kept. Be aware that the wrapper accepts `kept_run_keys` only as non-empty strings, whereas `varcode_build_exon_runs()` reports each key as a pair of integers (for example `[1, 1]`) and varcode compares the keys as tuples; no string matches, so any list of strings selects no run and gives an empty product. Only the default `None` returns a usable product.
- **API:** `varcode_realize_exon_path(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, kept_run_keys: 'Sequence[str] | None' = None, flank: 'int' = 50, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `kept_run_keys` is None or a list of non-empty strings (see the note above); `flank` is an integer from 0 to 100000 (default 50).
- **Output:** `transcript_id` and `product` with `transcript` (text), `cdna_sequence`, `protein_sequence`, `junction_signature`, `start_codon_present` and `evidence` (`start_codon`, `exon_runs`, `genomic_positions`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_realize_exon_path

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

everything = varcode_realize_exon_path(transcript_id="T1", flank=10, **ref)
product = everything["product"]
print(len(product["cdna_sequence"]), len(product["protein_sequence"]), product["start_codon_present"])
print(product["evidence"]["exon_runs"], len(product["evidence"]["genomic_positions"]))

by_string = varcode_realize_exon_path(transcript_id="T1", kept_run_keys=["1"], flank=10, **ref)
print(len(by_string["product"]["cdna_sequence"]), repr(by_string["product"]["protein_sequence"]))
```

- **Example output:**

```text
160 46 True
[[1, 1], [2, 1], [3, 1]] 160
0 ''
```

## 25) `dnakit.molecular_lab.varcode_realize_splice_plan()`

- **Purpose:** Realise one enumerated splice plan (exon skip, intron retention or a cryptic choice) on the transcript layout.
- **Method:** Builds the layout, enumerates the plans as `varcode_enumerate_layout_splice_plans()` does, computes the disrupted splice sites with `varcode.splice_graph.disrupted_splice_sites`, and calls `varcode.transcript_layout.realize_splice_plan` for the plan at `plan_index`. The product holds the spliced cDNA and the protein translated from the CDS start; for a cryptic choice the best non-canonical boundary within `scan_flank` bases is used.
- **API:** `varcode_realize_splice_plan(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, plan_index: 'int' = 0, flank: 'int' = 50, scan_flank: 'int' = 50, max_plans: 'int' = 64, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `plan_index` is an integer from 0 up to `max_plans` and must be smaller than the number of plans (otherwise `ConfigurationError`); `flank` and `scan_flank` are integers (`scan_flank` from 0 to 10000, default 50); `max_plans` is an integer from 1 to 10000.
- **Output:** `transcript_id`, `plan_index`, `plan_count`, `product` (as in `varcode_realize_exon_path()`, with `splice_choices` and `resolved_splice_choices` in its `evidence`) and `unresolved_mechanisms`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_realize_splice_plan

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

from dnakit.exceptions import ConfigurationError

for index in range(4):
    result = varcode_realize_splice_plan(transcript_id="T1", plan_index=index, flank=10, **ref)
    product = result["product"]
    mechanism = product["evidence"]["splice_choices"][0][1]
    print(index, mechanism, len(product["cdna_sequence"]), len(product["protein_sequence"]), result["plan_count"])

try:
    varcode_realize_splice_plan(transcript_id="T1", plan_index=9, flank=10, **ref)
except ConfigurationError as error:
    print(type(error).__name__, "out of range" in str(error))
```

- **Example output:**

```text
0 exon_skip 110 33 4
1 intron_retention 200 63 4
2 cryptic 134 41 4
3 normal 160 46 4
ConfigurationError True
```

## 26) `dnakit.molecular_lab.varcode_reference_base()`

- **Purpose:** Return the forward-strand reference base at a 1-based position.
- **Method:** Calls `varcode.reference_base(genome, contig, position)` on the local reference, which reads the genome FASTA.
- **API:** `varcode_reference_base(gtf: 'str', genome_fasta: 'str', contig: 'str', position: 'int', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `contig` is a plain identifier that exists in the reference; `position` is an integer of at least 1. A position beyond the end of the contig, or a contig that is not in the reference, gives an empty `base` and no error.
- **Output:** `contig`, `position` and `base`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_reference_base

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_base(contig="1", position=13, **ref)
print(result["contig"], result["position"], result["base"])
beyond = varcode_reference_base(contig="1", position=9999, **ref)
print(repr(beyond["base"]))
```

- **Example output:**

```text
1 13 C
''
```

## 27) `dnakit.molecular_lab.varcode_reference_range()`

- **Purpose:** Return the forward-strand reference sequence over an inclusive 1-based interval.
- **Method:** Calls `varcode.reference_range(genome, contig, start, end)` on the local reference, which reads the genome FASTA. The sequence is always on the forward strand, even for a gene on the minus strand.
- **API:** `varcode_reference_range(gtf: 'str', genome_fasta: 'str', contig: 'str', start: 'int', end: 'int', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `contig` is a plain identifier that exists in the reference; `start` is an integer of at least 1 and `end` an integer of at least `start`.
- **Output:** `contig`, `start`, `end` and `sequence`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_reference_range

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_range(contig="1", start=10, end=21, **ref)
print(result["sequence"], len(result["sequence"]))
```

- **Example output:**

```text
ATGCCAACAGTC 12
```

## 28) `dnakit.molecular_lab.varcode_reference_selenocysteine()`

- **Purpose:** Return the cDNA offsets of annotated selenocysteine codons in a transcript.
- **Method:** Calls `varcode.effects.selenocysteine.reference_selenocysteine`, which looks for `U` in the protein sequence of a complete transcript and returns the cDNA offset of the first base of each such codon. In a local reference the protein FASTA is derived from the CDS lines by DNAKit's own translation, which reads `TGA` as a stop and never produces `U`, so on a local reference this function returns an empty list for every transcript.
- **API:** `varcode_reference_selenocysteine(gtf: 'str', genome_fasta: 'str', transcript_id: 'str', *, cache_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference.
- **Output:** `transcript_id` and `cdna_offsets` (a list of integers).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_reference_selenocysteine

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

for transcript_id in ("T1", "T2"):
    result = varcode_reference_selenocysteine(transcript_id=transcript_id, **ref)
    print(transcript_id, result["cdna_offsets"])
```

- **Example output:**

```text
T1 []
T2 []
```

## 29) `dnakit.molecular_lab.varcode_make_fusion_outcome()`

- **Purpose:** Import one observed RNA junction or fusion model on a structural variant.
- **Method:** Builds a varcode `StructuralVariant` on the local reference and calls `varcode.rna_evidence.make_fusion_outcome` with the 5' anchor transcript, the observed sequence, the model ID, the optional partner transcript, the optional `cds_start` and the read count. When `partner_transcript_id` is given the outcome joins the two transcripts in sense orientation; otherwise it is a translocation to an intergenic region. With `cds_start` varcode translates a complete start-to-stop ORF of the sequence (a prediction from RNA, not a protein observation); without it the protein is not determined.
- **API:** `varcode_make_fusion_outcome(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', transcript_id: 'str', sequence: 'str', transcript_model_id: 'str', *, partner_transcript_id: 'str | None' = None, cds_start: 'int | None' = None, source: 'str' = 'rna', read_count: 'int | None' = None, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structural_variant` has `contig`, `start`, `end` and `sv_type`. `transcript_id` is the annotated 5' anchor; `sequence` is the observed 5'-to-3' sequence (A, C, G, T or N, up to 100000 bases); `transcript_model_id` is an identifier; `partner_transcript_id` is None or a transcript ID; `cds_start` is None or a 0-based offset where a complete ORF starts (otherwise `ConfigurationError`); `source` is an identifier (default `rna`); `read_count` is None or a non-negative integer.
- **Output:** `structural_variant`, `transcript_id`, `partner_transcript_id`, `transcript_model_id`, `candidate` (`source`, `evidence` and `effect`) and `mutant_transcript` (`cdna_sequence`, `mutant_protein_sequence`, `edits`, ...).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_make_fusion_outcome

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}
import csv

with open("structures.tsv", newline="") as handle:
    sequence = "".join(row["sequence"] for row in csv.DictReader(handle, delimiter="\t"))

sv = {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}
plain = varcode_make_fusion_outcome(
    structural_variant=sv, transcript_id="T1", sequence=sequence, transcript_model_id="M1",
    partner_transcript_id="T2", read_count=7, **ref
)
print(plain["candidate"]["effect"]["class"], plain["candidate"]["evidence"]["protein_status"], plain["mutant_transcript"]["mutant_protein_sequence"])

orf = varcode_make_fusion_outcome(
    structural_variant=sv, transcript_id="T1", sequence=sequence, transcript_model_id="M1",
    partner_transcript_id="T2", cds_start=22, read_count=7, **ref
)
print(orf["candidate"]["evidence"]["protein_status"], orf["mutant_transcript"]["mutant_protein_sequence"])
print(len(sequence), orf["candidate"]["source"], orf["candidate"]["evidence"]["read_count"])
```

- **Example output:**

```text
GeneFusion not_determined None
predicted_from_observed_rna MELYITIAYV
160 rna 7
```

## 30) `dnakit.molecular_lab.varcode_apply_rna_evidence_to_effects()`

- **Purpose:** Attach observed fusion outcomes to the predicted effects of a structural variant.
- **Method:** Builds one fusion outcome per entry as `varcode_make_fusion_outcome()` does, wraps them in a varcode `RNAEvidence` resolver, predicts the effects of the structural variant on the local reference and calls `varcode.apply_rna_evidence_to_effects`. The DNA-predicted effect is kept, and each effect lists the candidates it gained (varcode's own alternatives plus the observed RNA candidates).
- **API:** `varcode_apply_rna_evidence_to_effects(gtf: 'str', genome_fasta: 'str', structural_variant: 'Mapping[str, Any]', fusion_outcomes: 'list[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structural_variant` has `contig`, `start`, `end` and `sv_type`. `fusion_outcomes` is a list of 1 to 200 mappings with `transcript_id`, `transcript_model_id` and `sequence`, and optionally `partner_transcript_id`, `cds_start`, `read_count` and `source`.
- **Output:** `structural_variant`, `fusion_count`, `effects` (each with `class`, `description`, `transcript_id`, `gene_name` and `candidates`), `count` and `observed_candidate_count` (all candidates over all effects, including varcode's own).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_apply_rna_evidence_to_effects

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}
import csv

with open("structures.tsv", newline="") as handle:
    sequence = "".join(row["sequence"] for row in csv.DictReader(handle, delimiter="\t"))

fusion = {
    "transcript_id": "T1", "transcript_model_id": "M1", "sequence": sequence,
    "partner_transcript_id": "T2", "read_count": 7,
}
result = varcode_apply_rna_evidence_to_effects(structural_variant={"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}, fusion_outcomes=[fusion], **ref)
print(result["fusion_count"], result["count"], result["observed_candidate_count"])
effect = result["effects"][0]
print(effect["class"], effect["transcript_id"])
print([(c["source"], c["effect"]["class"]) for c in effect["candidates"]])
```

- **Example output:**

```text
1 1 6
GeneFusion T1
[('varcode', 'GeneFusion'), ('varcode', 'StructuralVariantEffect'), ('varcode_splice', 'ExonSkipping'), ('varcode_splice', 'IntronRetention'), ('varcode_splice', 'CrypticDonor'), ('rna', 'GeneFusion')]
```

## 31) `dnakit.molecular_lab.varcode_load_exacto_fusions()`

- **Purpose:** Import ExACTO transcript-structure and integration tables linked to structural variants.
- **Method:** Builds a varcode `StructuralVariant` for each DNA call ID and calls `varcode.exacto.load_exacto_fusions` with the two TSV files. Only the integration rows whose `dna_variant_call_id` is a key of `structural_variants` are imported. Sequences are taken in transcript orientation. A linear two-locus model becomes a fusion outcome whose reference-matching bases map to the annotated transcripts; circular, multi-gene and incomplete paths raise an error. Native primary-structure peptides and explicit ORF starts are not exposed, so the protein is not determined.
- **API:** `varcode_load_exacto_fusions(gtf: 'str', genome_fasta: 'str', structures_path: 'str', integrated_path: 'str', structural_variants: 'Mapping[str, Mapping[str, Any]]', *, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structures_path` and `integrated_path` are existing tab-separated files (gzip is accepted); the structure table needs the columns `transcript_model_id`, `reference_transcript_ids`, `index`, `read_start`, `read_end`, `sequence`, `type`, `kind`, `context`, `chromosome_1`, `position_1`, `strand_1`, `chromosome_2`, `position_2`, `strand_2`, `transcript_id_1` and `transcript_id_2`, and the integration table `transcript_model_id`, `reference_transcript_ids`, `rna_variant_call_id` and `dna_variant_call_id`. `structural_variants` maps 1 to 200 DNA call IDs to structural-variant mappings.
- **Output:** `structural_variants`, `candidates` (records with `source` set to `exacto`, `evidence` with the DNA and RNA call IDs, the reference transcripts, `partner_status` and the original table rows, and `effect`), `count`, `structures_path` and `integrated_path` (the resolved paths).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_load_exacto_fusions

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_load_exacto_fusions(
    structures_path="structures.tsv", integrated_path="integrated.tsv",
    structural_variants={"D1": {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}}, **ref
)
print(result["count"], list(result["structural_variants"]))
candidate = result["candidates"][0]
evidence = candidate["evidence"]
print(candidate["source"], candidate["effect"]["class"], candidate["effect"]["transcript_id"])
print(evidence["dna_variant_call_id"], evidence["rna_variant_call_ids"], evidence["reference_transcript_ids"], evidence["partner_status"])
print(len(evidence["exacto_structure"]), evidence["protein_status"])
```

- **Example output:**

```text
1 ['D1']
exacto GeneFusion T1
D1 ['R1'] ['T1', 'T2'] sense
3 not_determined
```

## 32) `dnakit.molecular_lab.varcode_reference_completion_hypotheses()`

- **Purpose:** Complete the mapped ends of ExACTO fusions as explicit reference hypotheses.
- **Method:** Loads the observed fusions as `varcode_load_exacto_fusions()` does, then calls `varcode.reference_completion.reference_completion_hypotheses` for each. The unmapped 5' and/or 3' end of the observed model is extended with reference bases of an annotated transcript, one hypothesis per choice of transcript. The hypotheses are labelled as reference-completion hypotheses, never as RNA observations, and the read counts are not used as support.
- **API:** `varcode_reference_completion_hypotheses(gtf: 'str', genome_fasta: 'str', structures_path: 'str', integrated_path: 'str', structural_variants: 'Mapping[str, Mapping[str, Any]]', *, ends: 'list[str] | None' = None, five_prime_transcript_ids: 'list[str] | None' = None, three_prime_transcript_ids: 'list[str] | None' = None, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `structures_path`, `integrated_path` and `structural_variants` as in `varcode_load_exacto_fusions()`. `ends` is None (both ends) or a list of `five_prime` and/or `three_prime`; `five_prime_transcript_ids` and `three_prime_transcript_ids` are None (the annotated transcript of that end) or lists of up to 500 transcript IDs that exist in the reference.
- **Output:** `structural_variants`, `ends`, `observed_count`, `hypothesis_count` and `hypotheses`. Each hypothesis has `source` (`varcode_reference_completion`), `evidence` (including `observed_evidence`, `reference_completion` with the observed and assumed spans, `protein_status`, `protein_completeness`, `cds_start` and `cds_end`) and `effect`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_reference_completion_hypotheses

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

result = varcode_reference_completion_hypotheses(
    structures_path="structures.tsv", integrated_path="integrated.tsv",
    structural_variants={"D1": {"contig": "1", "start": 151, "end": 280, "sv_type": "DEL"}}, ends=["three_prime"], **ref
)
print(result["ends"], result["observed_count"], result["hypothesis_count"])
hypothesis = result["hypotheses"][0]
evidence = hypothesis["evidence"]
print(hypothesis["source"], evidence["sequence_status"], evidence["protein_status"], evidence["protein_completeness"])
completion = evidence["reference_completion"]
print(completion["observed_span"], [(s["transcript_id"], s["start"], s["end"], s["label"]) for s in completion["assumed_spans"]])
```

- **Example output:**

```text
['three_prime'] 1 1
varcode_reference_completion reference_completion_hypothesis predicted_from_reference_completion start_to_stop
[0, 160] [('T2', 160, 210, 'assumed_reference_three_prime')]
```

## 33) `dnakit.molecular_lab.varcode_enumerate_splice_outcomes()`

- **Purpose:** Enumerate the splice mechanisms behind a splice-signal disruption on one transcript.
- **Method:** Predicts the variant's effect on the transcript, then calls `varcode.splice_outcomes.enumerate_splice_outcomes`. A splice donor, acceptor or splice-site effect is wrapped in a `SpliceOutcomeSet` whose candidates are the mechanisms (for example intron retention, a cryptic site and normal splicing). With `use_genomic_sequence` true the wrapper supplies the genome sequence of the local reference, so cryptic sites are scanned and intron retention is materialised (the descriptions then carry protein consequences). Any other effect is returned unchanged with no candidates. A single-exon transcript has no splice junction.
- **API:** `varcode_enumerate_splice_outcomes(gtf: 'str', genome_fasta: 'str', variant: 'Mapping[str, Any]', transcript_id: 'str', *, use_genomic_sequence: 'bool' = True, cache_dir: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variant` is a mapping with `contig`, `start` (1-based integer), `ref` and `alt` (A, C, G, T or N; not both empty). `transcript_id` is a plain identifier (letters, digits and `._:+|-`) that must exist in the local reference. `use_genomic_sequence` is a boolean (default True).
- **Output:** `variant`, `transcript_id`, `input_effect`, `outcome_class`, `is_outcome_set`, `genomic_sequence_used`, `candidates` (records with `source`, `evidence` and `effect`) and `count`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_enumerate_splice_outcomes

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

with_sequence = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 61, "ref": "G", "alt": "A"}, transcript_id="T1", **ref)
print(with_sequence["outcome_class"], with_sequence["is_outcome_set"], with_sequence["genomic_sequence_used"], with_sequence["count"])
for candidate in with_sequence["candidates"]:
    print(candidate["effect"]["class"], "|", candidate["effect"]["description"])

without = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 61, "ref": "G", "alt": "A"}, transcript_id="T1", use_genomic_sequence=False, **ref)
print(without["genomic_sequence_used"], [c["effect"]["description"] for c in without["candidates"]])

missense = varcode_enumerate_splice_outcomes(variant={"contig": "1", "start": 13, "ref": "C", "alt": "G"}, transcript_id="T1", **ref)
print(missense["outcome_class"], missense["is_outcome_set"], missense["count"])
```

- **Example output:**

```text
SpliceOutcomeSet True True 3
IntronRetention | intron-retention (donor side, p.G18fs*)
CrypticDonor | cryptic-donor ext8bp
NormalSplicing | normal-splicing
False ['intron-retention (donor side, predicted)', 'cryptic-donor (predicted)', 'normal-splicing']
Substitution False 0
```

## 34) `dnakit.molecular_lab.varcode_install_reference_data()`

- **Purpose:** Download and index the Ensembl reference data that varcode needs, for each reference name.
- **Method:** For every distinct name pyensembl picks a release through `pyensembl.genome_for_reference_name`: with `allow_older_downloaded_release` true the newest release that is already downloaded, otherwise the newest release. When the GTF and FASTA files of that release are missing, they are downloaded from Ensembl (`genome.download()`) and the index is built (`genome.index()`), as the varcode command line does; present files are never downloaded again. Names are checked before pyensembl is called, and an unknown name raises `ConfigurationError`.
- **API:** `varcode_install_reference_data(reference_names: 'Sequence[str]', *, allow_older_downloaded_release: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and, when the Ensembl files of a release are not yet present, internet access to the Ensembl FTP server and disk space for the human annotation (hundreds of megabytes). A download failure raises `DownloadError`; an indexing failure raises `BackendExecutionError`. The example below does not download anything: it shows the error for an unknown reference name, which is raised before any download.
- **Input:** `reference_names` is a list of 1 to 16 names of up to 64 characters made of letters, digits, `_`, `.` and `-` (for example `GRCh38`); duplicates are ignored. `allow_older_downloaded_release` is a boolean (default True).
- **Output:** `references`, a list with one record per name (`reference_name`, `release`, `species`, `cache_directory`, `files_present_before`, `downloaded`, `indexed`), and `allow_older_downloaded_release`.
- **Usage with the dependency installed:**

```python
from dnakit.molecular_lab import varcode_install_reference_data

result = varcode_install_reference_data(["GRCh38"])
for reference in result["references"]:
    print(reference["reference_name"], reference["release"], reference["downloaded"], reference["indexed"])
```

- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import varcode_install_reference_data

for names in (["NOT_A_GENOME"], [], ["bad name"]):
    try:
        varcode_install_reference_data(names, verbose=False)
    except ConfigurationError as error:
        message = str(error)
        print(type(error).__name__, "|", message.split("] ", 1)[-1].split(":")[0])
```

- **Example output:**

```text
ConfigurationError | Unknown reference name 'NOT_A_GENOME'
ConfigurationError | reference_names must hold 1 to 16 names.
ConfigurationError | Invalid reference name 'bad name'.
```

## 35) `dnakit.molecular_lab.varcode_variant_gene_overlaps()`

- **Purpose:** List the genes that overlap each variant, using a local GTF and genome FASTA.
- **Method:** Builds (or reuses) the local reference and, for each variant, creates a varcode `Variant` and reads its `genes` property. Variants that varcode cannot annotate (for example a contig that is not in the GTF) are reported with an error message, or raise `BackendExecutionError` when `raise_on_error` is true.
- **API:** `varcode_variant_gene_overlaps(gtf: 'str', genome_fasta: 'str', variants: 'Sequence[Mapping[str, Any]]', *, cache_dir: 'str | None' = None, raise_on_error: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `gtf` and `genome_fasta` are existing local files (the annotation and the genome FASTA) and `cache_dir` is the folder for the index (default `~/.cache/dnakit/varcode_reference`). `variants` is a list of 1 to 10000 mappings with `contig` (a name as in the GTF), `position` (1-based integer), `ref` and `alt` (DNA alleles of 1 to 256 letters A, C, G, T or N; unlike the other varcode functions, empty alleles are not allowed). `raise_on_error` is a boolean (default False).
- **Output:** `variants` (one record per input with `variant_index`, `variant`, `contig`, `position`, `ref`, `alt`, `genes` and `error`), `variant_count`, `annotation_failures` and `reference` (the digest, cache folder and whether the index was reused). Each gene has `gene_id`, `gene_name`, `biotype`, `contig`, `start`, `end` and `strand`.
- **Example code:**

```python
from dnakit.molecular_lab import varcode_variant_gene_overlaps

ref = {"gtf": "genes.gtf", "genome_fasta": "genome.fa", "cache_dir": "cache", "verbose": False}

variants = [
    {"contig": "1", "position": 13, "ref": "C", "alt": "G"},
    {"contig": "1", "position": 300, "ref": "G", "alt": "A"},
    {"contig": "1", "position": 255, "ref": "G", "alt": "A"},
    {"contig": "2", "position": 5, "ref": "C", "alt": "G"},
]
result = varcode_variant_gene_overlaps(variants=variants, **ref)
print(result["variant_count"], result["annotation_failures"])
for row in result["variants"]:
    print(row["variant"], [g["gene_name"] for g in row["genes"]], row["error"] is not None)
gene = result["variants"][0]["genes"][0]
print(gene["gene_id"], gene["biotype"], gene["start"], gene["end"], gene["strand"])
```

- **Example output:**

```text
4 1
1:13:C>G ['GENEA'] False
1:300:G>A ['GENEB'] False
1:255:G>A [] False
2:5:C>G [] True
G1 protein_coding 1 240 +
```

## 36) `dnakit.molecular_lab.varcode_check_sample_identity_files()`

- **Purpose:** Check whether the samples in several local VCF files are the same donors.
- **Method:** Loads every selected VCF sample with `varcode.sample_identity.load_vcf_samples` (as germline, somatic or mixed samples, plus the samples of an optional manifest), builds a `SampleCheckConfig` from `config`, and calls `varcode.sample_identity.check_sample_identity`, which compares every pair of samples by genotype concordance and shared somatic variants. The full report is written as JSON to `json_output` and a pair table to `tsv_output`; the report is also returned. The results are screening heuristics, not proof of identity. Existing output files are never overwritten.
- **API:** `varcode_check_sample_identity_files(*, germline_vcfs: 'Sequence[str]' = (), somatic_vcfs: 'Sequence[str]' = (), mixed_vcfs: 'Sequence[str]' = (), manifest_path: 'str | None' = None, assembly: 'str | None' = None, config: 'Mapping[str, Any] | None' = None, json_output: 'str | None' = None, tsv_output: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package varcode (extra `dnakit[bio-varcode]`). It works offline.
- **Input:** `germline_vcfs`, `somatic_vcfs` and `mixed_vcfs` are lists of existing VCF paths; `manifest_path` is an optional CSV or TSV with `path` and `sample` columns (and optionally `kind`, `assembly`, `role`, `label`, `donor_id`, `tumor_id`, `normal_sample`); `assembly` is `GRCh37` or `GRCh38` when the VCF header does not say; `config` is a mapping of `SampleCheckConfig` fields (unknown fields raise `ConfigurationError`); `json_output` and `tsv_output` are new paths in an existing folder, different from each other and from the inputs. At least two samples are needed.
- **Output:** `sample_count`, `pair_count`, `identity_conflicts`, `json_output`, `tsv_output` and `report` (with `samples`, `pairs` and the configuration; each pair has `germline` and `somatic` results with a `status` such as `compatible`, `discordant` or `inconclusive`, `flags` and `identity_conflict`). Sample labels in the report contain the file path.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_check_sample_identity_files

folder = Path(tempfile.mkdtemp())
result = varcode_check_sample_identity_files(
    germline_vcfs=["donor1.vcf", "donor2.vcf", "donor3.vcf"],
    assembly="GRCh38",
    config={"min_shared_snps": 5, "min_nonreference_snps": 1, "min_genomic_bins": 1},
    json_output=str(folder / "report.json"),
    tsv_output=str(folder / "pairs.tsv"),
    verbose=False,
)
print(result["sample_count"], result["pair_count"], result["identity_conflicts"])
for pair in result["report"]["pairs"]:
    germline = pair["germline"]
    print(germline["status"], germline["shared_snps"], germline["genotype_concordance"], germline["ibs0_rate"])
print(sorted(path.name for path in folder.iterdir()))
print((folder / "pairs.tsv").read_text().splitlines()[0].split("\t")[:4])
```

- **Example output:**

```text
3 3 0
compatible 5 1.0 0.0
discordant 5 0.0 0.4
discordant 5 0.0 0.4
['pairs.tsv', 'report.json']
['sample_a', 'sample_b', 'germline_status', 'shared_snps']
```

## 37) `dnakit.molecular_lab.varcode_compare_sv_calls_directory()`

- **Purpose:** Group structural-variant breakpoints from a CSV of calls and write the groups to a new directory.
- **Method:** Reads the CSV with `csv.DictReader`, calls `varcode.sv_comparison.compare_sv_calls` with `max_distance`, and writes `groups.csv`, `members.csv` and `comparison.json` (with input provenance: varcode version, SHA-256 and path of the input, call and group counts) to the new directory. Complete-link grouping puts calls with the same adjacency topology together when every corresponding breakpoint lies within `max_distance` bases. A group whose calls have different normalised alleles is `nearby_candidate`; a single call is `exact_reported_allele`.
- **API:** `varcode_compare_sv_calls_directory(calls_csv: 'str', output_dir: 'str', *, max_distance: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package varcode (extra `dnakit[bio-varcode]`). It works offline.
- **Input:** `calls_csv` is an existing CSV with the columns `call_id` (unique), `caller`, `sample`, `build` (`GRCh37` or `GRCh38`), `chrom`, `pos`, `ref` and `alt`, and optionally `info` (VCF INFO, for example `END=10500;SVTYPE=DEL`) and `end`. `output_dir` must not exist yet and its parent folder must exist; `max_distance` is a non-negative integer (default 100).
- **Output:** `call_count`, `group_count`, `output_dir` and `files` (`groups.csv`, `members.csv`, `comparison.json`).
- **Example code:**

```python
import csv
import json
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_compare_sv_calls_directory

target = Path(tempfile.mkdtemp()) / "comparison"
result = varcode_compare_sv_calls_directory("calls.csv", str(target), max_distance=50, verbose=False)
print(result["call_count"], result["group_count"], result["files"])
print(sorted(path.name for path in target.iterdir()))
with open(target / "groups.csv", newline="") as handle:
    for row in csv.DictReader(handle):
        print(row["status"], row["call_ids"], row["call_count"])
print(json.loads((target / "comparison.json").read_text())["provenance"]["group_count"])
```

- **Example output:**

```text
3 2 ['groups.csv', 'members.csv', 'comparison.json']
['comparison.json', 'groups.csv', 'members.csv']
exact_reported_allele ["c3"] 1
nearby_candidate ["c1", "c2"] 2
2
```

## 38) `dnakit.molecular_lab.varcode_effects_table()`

- **Purpose:** Predict the coding effects of loaded variants against a local reference and write them to a CSV file.
- **Method:** Loads the variants as `varcode_load_variant_sources()` does, passes each one to `varcode_predict_variant_effects()` with the same `gtf` and `genome_fasta`, and writes one CSV row per effect. With `skip_errors` the annotation failures are kept as rows of class `Failure`; otherwise the first failure raises. A variant without any effect gets no row. The `--only-coding` and `--one-per-variant` filters of the varcode command are not included.
- **API:** `varcode_effects_table(output_csv: 'str', *, gtf: 'str', genome_fasta: 'str', vcf_paths: 'Sequence[str]' = (), maf_paths: 'Sequence[str]' = (), json_paths: 'Sequence[str]' = (), variants: 'Sequence[Mapping[str, Any]]' = (), cache_dir: 'str | None' = None, only_passing: 'bool' = True, skip_errors: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages varcode (extra `dnakit[bio-varcode]`), pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`). It works offline; nothing is downloaded.
- **Input:** `output_csv` is a new file path in an existing folder. `gtf` and `genome_fasta` are required keyword arguments; `vcf_paths`, `maf_paths` and `json_paths` are lists of existing files; `variants` is a list of variant mappings (`contig`, `start`, `ref`, `alt`); at least one source is needed. `cache_dir` is the index folder; `only_passing` (default True) skips VCF records whose FILTER is not PASS; `skip_errors` is a boolean (default False).
- **Output:** `reference` (the digest of the local reference), `variant_count`, `effect_count`, `class_counts` (effect classes with their counts), `columns` (`contig`, `start`, `ref`, `alt`, `class`, `transcript_id`, `gene_name`, `description`) and `output_csv`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import varcode_effects_table

target = Path(tempfile.mkdtemp()) / "effects.csv"
result = varcode_effects_table(
    str(target),
    gtf="genes.gtf",
    genome_fasta="genome.fa",
    vcf_paths=["phased.vcf"],
    variants=[{"contig": "1", "start": 61, "ref": "G", "alt": "A"}],
    cache_dir="cache",
    verbose=False,
)
print(result["variant_count"], result["effect_count"], result["class_counts"])
print(result["columns"])
print(target.read_text().strip())
```

- **Example output:**

```text
3 3 {'SpliceOutcomeSet': 1, 'Substitution': 2}
['contig', 'start', 'ref', 'alt', 'class', 'transcript_id', 'gene_name', 'description']
contig,start,ref,alt,class,transcript_id,gene_name,description
1,13,C,G,Substitution,T1,GENEA,p.P2A
1,17,C,T,Substitution,T1,GENEA,p.T3I
1,61,G,A,SpliceOutcomeSet,T1,GENEA,"splice-set:intron-retention (donor side, predicted)"
```

## 39) `dnakit.molecular_lab.varcode_load_variant_sources()`

- **Purpose:** Load variants from VCF, MAF and JSON files and explicit records, and union them.
- **Method:** Loads each VCF with `varcode.load_vcf` (structural variants parsed, `only_passing` applied), each MAF with `varcode.maf.load_maf`, the explicit records as `varcode.Variant` objects on the reference, and each JSON file with `VariantCollection.from_json`, then combines the collections with `VariantCollection.union`. The reference is built from `gtf` and `genome_fasta` (the two must be given together), or `genome` names an installed pyensembl reference. VCF contig names are normalised by varcode (a `chr` prefix is dropped), while explicit records keep their contig names as given, so the same variant can appear twice if its contig is spelled differently in the two sources.
- **API:** `varcode_load_variant_sources(vcf_paths: 'Sequence[str]' = (), maf_paths: 'Sequence[str]' = (), json_paths: 'Sequence[str]' = (), variants: 'Sequence[Mapping[str, Any]]' = (), *, gtf: 'str | None' = None, genome_fasta: 'str | None' = None, cache_dir: 'str | None' = None, genome: 'str | None' = None, only_passing: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package varcode (extra `dnakit[bio-varcode]`). With `gtf` and `genome_fasta` it also needs pyensembl (extra `dnakit[bio-pyensembl]`) and pysam (extra `dnakit[bio-reads]`); with `genome` it needs that Ensembl release to be installed already. Nothing is downloaded.
- **Input:** `vcf_paths`, `maf_paths` and `json_paths` are lists of existing files; `variants` is a list of mappings with `contig`, `start` (integer of at least 1), `ref` and `alt`; at least one source is needed, and explicit records need a reference. `gtf` and `genome_fasta` are given together or not at all; `cache_dir` is the index folder; `genome` is a reference name such as `GRCh38`; `only_passing` is a boolean (default True).
- **Output:** `reference` (`local:<digest>`, the `genome` name or None), `count`, `sources` (the number of files or records taken from `vcf`, `maf`, `json` and `variant`) and `variants` (records with `contig`, `start`, `ref`, `alt`).
- **Example code:**

```python
from dnakit.molecular_lab import varcode_load_variant_sources

result = varcode_load_variant_sources(
    vcf_paths=["phased.vcf"],
    variants=[
        {"contig": "1", "start": 61, "ref": "G", "alt": "A"},
        {"contig": "1", "start": 13, "ref": "C", "alt": "G"},
    ],
    gtf="genes.gtf",
    genome_fasta="genome.fa",
    cache_dir="cache",
    verbose=False,
)
print(result["reference"].startswith("local:"), result["count"], result["sources"])
for variant in result["variants"]:
    print(variant["contig"], variant["start"], variant["ref"], variant["alt"])
```

- **Example output:**

```text
True 3 {'vcf': 1, 'maf': 0, 'json': 0, 'variant': 2}
1 13 C G
1 17 C T
1 61 G A
```
