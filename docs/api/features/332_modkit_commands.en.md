# modkit: modified-base pileups, modBAM tools, motifs and differential methylation

This group of functions in `dnakit.epigenome` is a set of command-line facades over [modkit](https://github.com/nanoporetech/modkit), the Oxford Nanopore tool for base-modification calls stored in a modBAM (the `MM` and `ML` tags). They cover the bedMethyl pileup (`modkit pileup`, `pileup-hemi`), read-level tables (`extract`), summaries and quality checks of a modBAM (`modbam summary`, `sample-probs`, `check-tags`), editing a modBAM (`adjust-mods`, `call-mods`, `update-tags`, `repair`), bedMethyl utilities (`bedmethyl merge`, `tobigwig`, `map-to-genome`, `stats`, `localize`), sequence motifs (`motif bed`, `search`, `evaluate`, `refine`), methylation entropy, validation against ground truth, differential methylation (`dmr`) and 6mA open-chromatin prediction. Each function runs one fixed subcommand with an explicit argument list (no shell), refuses to overwrite an existing output, and returns a dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `operation` and `provenance`; tables also get `row_count` and a short `rows` preview. The `status` is `completed` only when modkit exits with code 0 and every expected output exists; a modkit error gives `failed` and the reason is in `log_tail`.

modkit 0.6.4 is the tested version; it is released under the Oxford Nanopore Technologies Public License v1.0 and DNAKit starts it as a separate program without copying any of its code. Install a release binary from the project's GitHub releases page (or `conda install -c bioconda modkit`) and put it on PATH, or pass `modkit_path`; if the program is missing the function raises `BackendUnavailableError` with an install hint before any output is created. An option left at `None` or `False` is not passed, so modkit's own defaults apply. Positional arguments are always placed before options because modkit's multi-value options (such as `--modified-bases`) would otherwise swallow them. The examples call every function with `verbose=False` and build their small modBAM files from the SAM text below with pysam.

## Example data

The examples use the small files below: a two-contig reference, three sets of 100 bp reads with CpG 5mC calls (`MM:Z:C+m?` and `ML` tags; `other.sam` has few modified calls, `duplex.sam` also has the opposite-strand `G-m?` calls that hemi-methylation needs), a two-transcript reference with two transcriptome-aligned read sets, and small BED, chromosome-size and GTF files. The preparation step below sorts and indexes the BAM files and writes bgzip-compressed, tabix-indexed bedMethyl files into a `data` folder. Every other output goes to a fresh temporary folder because the functions never overwrite a file.

**Example file:** `ref.fa`

```text
>chr1
GCGCTCCTCCGGATTCGAAAATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAACGTGCATCGATGCCGGAAGCG
>chr2
ATTCGGCACGGACGTCGATCGATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGATGACGCGCCAGCTTCGAAAA
```

Reference with two contigs (240 bp and 160 bp).

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:240
@SQ	SN:chr2	LN:160
a001	16	chr1	1	60	100M	*	0	0	GCGCTCCTCCGGATTCGAAAATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCG	*	NM:i:0	MM:Z:C+m?,0,0,1,0,0,1,0,0,0,1,2,0,0,0,0,0,0,1,0;	ML:B:C,254,4,231,230,250,231,253,238,200,217,214,206,1,241,224,227,233,228,35
a002	0	chr1	21	60	100M	*	0	0	ATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCC	*	NM:i:0	MM:Z:C+m?,1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,0;	ML:B:C,214,248,218,26,35,6,246,7,246,32,27,242,218,231,32,254,215,225
a003	16	chr1	41	60	100M	*	0	0	GTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCG	*	NM:i:0	MM:Z:C+m?,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,2,0;	ML:B:C,211,244,247,228,206,233,23,201,219,239,225,210,200,234,35,232,254,229
a004	0	chr1	61	60	100M	*	0	0	TCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTC	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,0,0,0,1,0,4,1,1,1,0,1,1;	ML:B:C,35,246,250,32,33,213,203,223,212,26,222,200,239,221,201,40,237
a005	16	chr1	81	60	100M	*	0	0	GCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTAT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,3,0,0,1,1,0,0,0,0,1,0;	ML:B:C,205,251,16,243,255,200,217,207,211,204,216,210,241,229,231
a006	0	chr1	101	60	100M	*	0	0	AACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACG	*	NM:i:0	MM:Z:C+m?,1,0,4,1,1,1,0,1,1,5,2,0,0,1,0,2;	ML:B:C,201,221,212,216,32,38,201,225,246,28,243,253,240,233,233,225
a007	16	chr1	121	60	100M	*	0	0	GGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAAC	*	NM:i:0	MM:Z:C+m?,0,1,1,0,0,2,0,0,0,1,0,3,0,0,1;	ML:B:C,251,240,247,213,19,204,219,226,208,254,252,236,252,39,224
a008	0	chr1	141	60	100M	*	0	0	CGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAACGTGCATCGATGCCGGAAGCG	*	NM:i:0	MM:Z:C+m?,0,1,1,5,2,0,0,1,0,2,3,0,1,1,1,0;	ML:B:C,206,243,37,206,24,231,239,18,212,36,221,217,253,235,253,231
a009	16	chr2	1	60	100M	*	0	0	ATTCGGCACGGACGTCGATCGATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGA	*	NM:i:0	MM:Z:C+m?,1,0,0,1,0,1,0,0,3,0,0,0,2,0,0,1,0,0,1,1;	ML:B:C,234,246,208,234,248,232,23,207,255,31,235,220,204,250,208,239,224,235,205,23
a010	0	chr2	21	60	100M	*	0	0	GATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCG	*	NM:i:0	MM:Z:C+m?,0,1,0,0,1,0,0,0,0,1,0,1,2,0,2,0,0,3;	ML:B:C,36,207,217,202,0,200,207,2,250,26,228,215,254,224,34,18,245,206
a011	16	chr2	41	60	100M	*	0	0	TGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGA	*	NM:i:0	MM:Z:C+m?,1,1,3,0,1,2,0,0,1,0,1,0,0,3,0,0,0;	ML:B:C,220,200,218,220,220,204,38,7,250,234,30,216,213,215,252,248,241
a012	0	chr2	61	60	100M	*	0	0	AGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGATGACGCGCCAGCTTCGAAAA	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,2,0,2,0,0,3,0,0,0,0,3;	ML:B:C,221,24,2,220,237,19,206,237,5,201,25,235,4,200,250,230
```

Eight plus-strand and reverse-strand reads on chr1 and four on chr2; 80 percent of the CpG calls are modified.

**Example file:** `other.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:240
@SQ	SN:chr2	LN:160
b001	16	chr1	1	60	100M	*	0	0	GCGCTCCTCCGGATTCGAAAATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCG	*	NM:i:0	MM:Z:C+m?,0,0,1,0,0,1,0,0,0,1,2,0,0,0,0,0,0,1,0;	ML:B:C,3,223,19,13,37,27,32,34,32,2,23,20,27,33,211,1,211,232,32
b002	0	chr1	21	60	100M	*	0	0	ATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCC	*	NM:i:0	MM:Z:C+m?,1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,0;	ML:B:C,11,28,33,23,22,28,248,29,15,31,22,29,29,35,31,20,10,39
b003	16	chr1	41	60	100M	*	0	0	GTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCG	*	NM:i:0	MM:Z:C+m?,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,2,0;	ML:B:C,30,32,32,37,13,23,39,21,12,6,241,237,6,8,15,3,2,223
b004	0	chr1	61	60	100M	*	0	0	TCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTC	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,0,0,0,1,0,4,1,1,1,0,1,1;	ML:B:C,243,207,1,201,8,10,33,24,15,202,239,7,31,228,38,16,39
b005	16	chr1	81	60	100M	*	0	0	GCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTAT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,3,0,0,1,1,0,0,0,0,1,0;	ML:B:C,30,5,20,1,8,25,20,221,38,1,8,16,210,229,32
b006	0	chr1	101	60	100M	*	0	0	AACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACG	*	NM:i:0	MM:Z:C+m?,1,0,4,1,1,1,0,1,1,5,2,0,0,1,0,2;	ML:B:C,2,14,4,37,39,16,17,0,224,7,5,6,248,13,242,19
b007	16	chr1	121	60	100M	*	0	0	GGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAAC	*	NM:i:0	MM:Z:C+m?,0,1,1,0,0,2,0,0,0,1,0,3,0,0,1;	ML:B:C,24,13,27,1,3,33,6,30,33,7,18,23,26,219,249
b008	0	chr1	141	60	100M	*	0	0	CGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAACGTGCATCGATGCCGGAAGCG	*	NM:i:0	MM:Z:C+m?,0,1,1,5,2,0,0,1,0,2,3,0,1,1,1,0;	ML:B:C,1,3,31,37,0,23,4,31,236,29,222,7,7,254,25,6
b009	16	chr2	1	60	100M	*	0	0	ATTCGGCACGGACGTCGATCGATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGA	*	NM:i:0	MM:Z:C+m?,1,0,0,1,0,1,0,0,3,0,0,0,2,0,0,1,0,0,1,1;	ML:B:C,242,2,31,29,223,33,26,18,10,16,5,36,222,234,226,5,2,218,21,33
b010	0	chr2	21	60	100M	*	0	0	GATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCG	*	NM:i:0	MM:Z:C+m?,0,1,0,0,1,0,0,0,0,1,0,1,2,0,2,0,0,3;	ML:B:C,9,27,233,32,10,15,22,36,29,1,24,11,3,25,26,30,35,5
b011	16	chr2	41	60	100M	*	0	0	TGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGA	*	NM:i:0	MM:Z:C+m?,1,1,3,0,1,2,0,0,1,0,1,0,0,3,0,0,0;	ML:B:C,14,12,24,40,20,29,6,213,38,13,6,35,12,37,31,8,243
b012	0	chr2	61	60	100M	*	0	0	AGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGATGACGCGCCAGCTTCGAAAA	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,2,0,2,0,0,3,0,0,0,0,3;	ML:B:C,16,11,13,4,31,4,31,21,32,5,22,248,25,40,8,231
```

The same layout with only about 20 percent modified CpG calls (a second condition).

**Example file:** `duplex.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:240
@SQ	SN:chr2	LN:160
d001	16	chr1	1	60	100M	*	0	0	GCGCTCCTCCGGATTCGAAAATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCG	*	NM:i:0	MM:Z:C+m?,0,0,1,0,0,1,0,0,0,1,2,0,0,0,0,0,0,1,0;G-m?,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,1,1,0,4;	ML:B:C,234,238,237,200,30,214,245,253,225,9,209,33,200,4,237,249,255,238,24,27,251,208,6,231,243,240,26,224,234,237,221,1,38,10,20,36,245,40
d002	0	chr1	21	60	100M	*	0	0	ATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCC	*	NM:i:0	MM:Z:C+m?,1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,0;G-m?,0,0,0,0,0,0,2,1,0,0,0,1,0,0,1,0,0,0;	ML:B:C,36,207,254,30,251,209,227,7,239,24,21,217,202,204,234,212,18,9,21,208,24,255,241,35,251,227,15,27,33,221,226,20,239,208,240,222
d003	16	chr1	41	60	100M	*	0	0	GTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCG	*	NM:i:0	MM:Z:C+m?,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,2,0;G-m?,0,1,1,1,4,0,1,0,0,0,1,0,1,0,0,0,0,1;	ML:B:C,22,17,1,243,223,229,238,223,248,238,250,249,36,8,214,17,211,6,220,14,254,5,241,236,214,2,212,253,11,21,5,22,226,250,29,226
d004	0	chr1	61	60	100M	*	0	0	TCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTC	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,0,0,0,1,0,4,1,1,1,0,1,1;G-m?,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,3,0;	ML:B:C,236,226,200,253,32,245,247,248,33,34,214,37,7,2,251,12,27,200,207,219,201,226,239,208,255,251,22,207,252,215,8,0,236,3
d005	16	chr1	81	60	100M	*	0	0	GCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTAT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,3,0,0,1,1,0,0,0,0,1,0;G-m?,0,2,5,1,1,0,1,1,1,4,0,1,0,0,0;	ML:B:C,15,233,203,249,3,2,204,202,5,231,220,224,37,216,221,208,245,5,202,238,34,251,227,240,244,244,201,234,237,227
d006	0	chr1	101	60	100M	*	0	0	AACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACG	*	NM:i:0	MM:Z:C+m?,1,0,4,1,1,1,0,1,1,5,2,0,0,1,0,2;G-m?,0,0,0,1,1,0,0,3,0,1,0,0,0,2,0,0;	ML:B:C,208,20,235,16,244,252,33,6,36,251,200,215,2,236,224,201,254,254,230,18,251,249,233,6,206,3,255,252,211,29,25,238
d007	16	chr1	121	60	100M	*	0	0	GGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAAC	*	NM:i:0	MM:Z:C+m?,0,1,1,0,0,2,0,0,0,1,0,3,0,0,1;G-m?,2,0,3,2,0,1,0,0,2,5,1,1,0,1,1;	ML:B:C,222,205,215,26,10,236,33,9,9,231,33,37,11,248,237,214,34,18,26,237,217,19,230,12,223,230,9,244,238,237
d008	0	chr1	141	60	100M	*	0	0	CGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAACGTGCATCGATGCCGGAAGCG	*	NM:i:0	MM:Z:C+m?,0,1,1,5,2,0,0,1,0,2,3,0,1,1,1,0;G-m?,0,3,0,1,0,0,0,2,0,0,1,1,0,1,1,2;	ML:B:C,35,246,225,2,214,4,16,212,8,245,2,16,202,227,251,205,218,228,243,201,227,204,237,31,234,217,227,233,6,244,249,228
d009	16	chr2	1	60	100M	*	0	0	ATTCGGCACGGACGTCGATCGATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGA	*	NM:i:0	MM:Z:C+m?,1,0,0,1,0,1,0,0,3,0,0,0,2,0,0,1,0,0,1,1;G-m?,0,2,0,2,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1;	ML:B:C,243,16,221,34,242,222,218,36,40,223,29,235,221,38,219,29,250,204,211,34,225,217,203,2,217,249,221,28,249,231,7,237,207,246,244,226,25,8,209,12
d010	0	chr2	21	60	100M	*	0	0	GATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCG	*	NM:i:0	MM:Z:C+m?,0,1,0,0,1,0,0,0,0,1,0,1,2,0,2,0,0,3;G-m?,2,0,0,2,0,0,0,3,0,0,1,0,1,0,0,2,1,0;	ML:B:C,210,212,23,1,28,24,237,231,19,30,212,0,248,31,240,212,13,232,28,218,9,250,203,239,204,201,20,222,8,253,2,21,13,12,248,220
d011	16	chr2	41	60	100M	*	0	0	TGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGA	*	NM:i:0	MM:Z:C+m?,1,1,3,0,1,2,0,0,1,0,1,0,0,3,0,0,0;G-m?,0,0,0,3,0,0,2,0,2,1,0,1,0,0,0,0,1;	ML:B:C,202,254,210,230,34,13,19,229,25,229,16,223,223,214,13,223,229,210,210,232,201,238,231,253,228,202,3,202,201,14,224,210,207,207
d012	0	chr2	61	60	100M	*	0	0	AGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGATGACGCGCCAGCTTCGAAAA	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,2,0,2,0,0,3,0,0,0,0,3;G-m?,3,0,0,1,0,1,0,0,2,1,0,3,1,2,0,1;	ML:B:C,251,17,29,19,235,252,222,203,5,200,243,210,230,241,221,248,251,222,2,248,25,218,226,224,234,25,247,24,35,211,226,228
```

Reads that carry the top-strand `C+m?` and the bottom-strand `G-m?` calls on the same molecule.

**Example file:** `tx.fa`

```text
>tx1
ACGTTGTCTTCGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGCGCTAGATTCG
>tx2
ACGTTATTCGAACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACTACGTGTAATG
```

Transcriptome reference with transcripts `tx1` (160 nt) and `tx2` (120 nt).

**Example file:** `tx_a.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:tx1	LN:160
@SQ	SN:tx2	LN:120
ta001	16	tx1	1	60	60M	*	0	0	ACGTTGTCTTCGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTT	*	NM:i:0	MM:Z:C+m?,0,0,0,2,2,1,0,0,1;	ML:B:C,206,230,204,235,3,234,249,206,201
ta002	0	tx1	11	60	60M	*	0	0	CGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTC	*	NM:i:0	MM:Z:C+m?,0,2,4,1,1,0,0,0,1;	ML:B:C,251,217,219,255,254,223,238,224,211
ta003	16	tx1	21	60	60M	*	0	0	GAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTA	*	NM:i:0	MM:Z:C+m?,0,0,0,0,0,0,2,2,1;	ML:B:C,217,252,35,200,36,219,232,227,227
ta004	0	tx1	31	60	60M	*	0	0	TCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGG	*	NM:i:0	MM:Z:C+m?,1,1,0,0,0,1,0,1,0;	ML:B:C,214,252,205,240,33,230,209,12,212
ta005	16	tx1	41	60	60M	*	0	0	CGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAA	*	NM:i:0	MM:Z:C+m?,1,2,0,0,0,0,0,0,2;	ML:B:C,228,222,237,235,220,203,217,239,7
ta006	0	tx1	51	60	60M	*	0	0	CGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACAT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,0,1,0;	ML:B:C,211,201,244,218,220,218,209,226
ta007	16	tx1	61	60	60M	*	0	0	TCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGG	*	NM:i:0	MM:Z:C+m?,4,0,1,2,0,0,0;	ML:B:C,39,204,212,18,224,210,200
ta008	0	tx1	71	60	60M	*	0	0	GAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGA	*	NM:i:0	MM:Z:C+m?,1,0,1,0,4,0,1;	ML:B:C,229,250,218,228,27,7,203
ta009	16	tx1	81	60	60M	*	0	0	ATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAG	*	NM:i:0	MM:Z:C+m?,1,0,0,1,0,4,0,1,2;	ML:B:C,238,209,234,215,202,233,226,212,215
ta010	0	tx1	91	60	60M	*	0	0	TTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGC	*	NM:i:0	MM:Z:C+m?,0,0,4,0,1,0,0,0,1;	ML:B:C,231,226,241,13,202,216,213,226,9
ta011	16	tx1	101	60	60M	*	0	0	ACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGCGCTAGATTCG	*	NM:i:0	MM:Z:C+m?,0,1,1,2,0,0,1,0,4,0;	ML:B:C,220,236,241,250,202,205,255,210,242,241
ta012	0	tx2	1	60	60M	*	0	0	ACGTTATTCGAACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCC	*	NM:i:0	MM:Z:C+m?,0,0,2,1,0,0,0,0,1,1,2;	ML:B:C,226,241,251,225,4,17,12,239,217,203,40
ta013	16	tx2	11	60	60M	*	0	0	AACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCA	*	NM:i:0	MM:Z:C+m?,0,0,2,3,0,0,0,0,1,0,0;	ML:B:C,230,224,200,210,216,251,224,35,212,238,235
ta014	0	tx2	21	60	60M	*	0	0	ACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTA	*	NM:i:0	MM:Z:C+m?,0,0,0,0,0,1,1,2,3,2,1;	ML:B:C,230,28,202,207,255,249,223,233,19,252,234
ta015	16	tx2	31	60	60M	*	0	0	CGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGA	*	NM:i:0	MM:Z:C+m?,0,0,1,1,0,2,3,0,0,0;	ML:B:C,213,219,248,225,207,23,227,245,204,212
ta016	0	tx2	41	60	60M	*	0	0	ATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCA	*	NM:i:0	MM:Z:C+m?,1,2,3,2,1,0,2,1;	ML:B:C,253,207,210,210,209,255,234,0
ta017	16	tx2	51	60	60M	*	0	0	CTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,1,0,2;	ML:B:C,200,207,231,238,204,226,222,211
ta018	0	tx2	61	60	60M	*	0	0	AGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACTACGTGTAATG	*	NM:i:0	MM:Z:C+m?,0,2,1,0,2,1,3,1;	ML:B:C,40,203,220,229,236,255,228,228
```

Reads aligned to the transcripts, mostly modified (condition A).

**Example file:** `tx_b.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:tx1	LN:160
@SQ	SN:tx2	LN:120
tb001	16	tx1	1	60	60M	*	0	0	ACGTTGTCTTCGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTT	*	NM:i:0	MM:Z:C+m?,0,0,0,2,2,1,0,0,1;	ML:B:C,22,33,229,15,10,230,24,36,13
tb002	0	tx1	11	60	60M	*	0	0	CGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTC	*	NM:i:0	MM:Z:C+m?,0,2,4,1,1,0,0,0,1;	ML:B:C,11,24,4,239,8,255,249,10,18
tb003	16	tx1	21	60	60M	*	0	0	GAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTA	*	NM:i:0	MM:Z:C+m?,0,0,0,0,0,0,2,2,1;	ML:B:C,12,40,12,24,23,9,21,38,38
tb004	0	tx1	31	60	60M	*	0	0	TCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGG	*	NM:i:0	MM:Z:C+m?,1,1,0,0,0,1,0,1,0;	ML:B:C,21,222,30,11,11,201,22,1,26
tb005	16	tx1	41	60	60M	*	0	0	CGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAA	*	NM:i:0	MM:Z:C+m?,1,2,0,0,0,0,0,0,2;	ML:B:C,37,28,211,12,215,29,22,16,6
tb006	0	tx1	51	60	60M	*	0	0	CGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACAT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,0,1,0;	ML:B:C,23,2,5,32,9,34,243,11
tb007	16	tx1	61	60	60M	*	0	0	TCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGG	*	NM:i:0	MM:Z:C+m?,4,0,1,2,0,0,0;	ML:B:C,40,244,30,3,238,25,247
tb008	0	tx1	71	60	60M	*	0	0	GAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGA	*	NM:i:0	MM:Z:C+m?,1,0,1,0,4,0,1;	ML:B:C,16,26,240,231,13,246,8
tb009	16	tx1	81	60	60M	*	0	0	ATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAG	*	NM:i:0	MM:Z:C+m?,1,0,0,1,0,4,0,1,2;	ML:B:C,26,227,3,18,239,33,31,20,18
tb010	0	tx1	91	60	60M	*	0	0	TTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGC	*	NM:i:0	MM:Z:C+m?,0,0,4,0,1,0,0,0,1;	ML:B:C,9,252,11,31,11,217,35,23,250
tb011	16	tx1	101	60	60M	*	0	0	ACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGCGCTAGATTCG	*	NM:i:0	MM:Z:C+m?,0,1,1,2,0,0,1,0,4,0;	ML:B:C,37,19,17,16,18,21,37,30,16,17
tb012	0	tx2	1	60	60M	*	0	0	ACGTTATTCGAACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCC	*	NM:i:0	MM:Z:C+m?,0,0,2,1,0,0,0,0,1,1,2;	ML:B:C,32,22,22,26,11,28,21,33,12,30,5
tb013	16	tx2	11	60	60M	*	0	0	AACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCA	*	NM:i:0	MM:Z:C+m?,0,0,2,3,0,0,0,0,1,0,0;	ML:B:C,26,37,26,35,40,1,37,11,11,2,14
tb014	0	tx2	21	60	60M	*	0	0	ACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTA	*	NM:i:0	MM:Z:C+m?,0,0,0,0,0,1,1,2,3,2,1;	ML:B:C,8,211,35,229,39,213,23,225,26,35,2
tb015	16	tx2	31	60	60M	*	0	0	CGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGA	*	NM:i:0	MM:Z:C+m?,0,0,1,1,0,2,3,0,0,0;	ML:B:C,19,234,37,22,26,36,23,10,36,12
tb016	0	tx2	41	60	60M	*	0	0	ATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCA	*	NM:i:0	MM:Z:C+m?,1,2,3,2,1,0,2,1;	ML:B:C,38,252,0,20,39,20,36,27
tb017	16	tx2	51	60	60M	*	0	0	CTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACT	*	NM:i:0	MM:Z:C+m?,0,0,1,0,1,1,0,2;	ML:B:C,31,24,10,16,19,26,219,18
tb018	0	tx2	61	60	60M	*	0	0	AGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACTACGTGTAATG	*	NM:i:0	MM:Z:C+m?,0,2,1,0,2,1,3,1;	ML:B:C,201,239,15,8,0,35,30,15
```

Reads aligned to the transcripts, mostly unmodified (condition B).

**Example file:** `sizes.txt`

```text
chr1	240
chr2	160
```

Chromosome sizes of `ref.fa`.

**Example file:** `regions.bed`

```text
chr1	0	120	reg1
chr1	120	240	reg2
chr2	0	160	reg3
```

Three regions to summarise or compare.

**Example file:** `features.bed`

```text
chr1	60	61	f1	0	+
chr1	180	181	f2	0	+
```

Two point features for `localize`.

**Example file:** `truth.bed`

```text
chr1	1	2	m	0	+
chr1	9	10	-	0	+
chr1	15	16	m	0	+
chr1	26	27	-	0	+
chr1	30	31	m	0	+
chr1	35	36	-	0	+
chr1	39	40	m	0	+
chr1	42	43	-	0	+
chr1	46	47	m	0	+
chr1	53	54	-	0	+
chr1	58	59	m	0	+
chr1	61	62	-	0	+
chr1	68	69	m	0	+
chr1	74	75	-	0	+
chr1	79	80	m	0	+
chr1	84	85	-	0	+
chr1	87	88	m	0	+
chr1	92	93	-	0	+
chr1	98	99	m	0	+
chr1	104	105	-	0	+
chr1	108	109	m	0	+
chr1	119	120	-	0	+
chr1	126	127	m	0	+
chr1	132	133	-	0	+
chr1	138	139	m	0	+
chr1	140	141	-	0	+
chr1	147	148	m	0	+
chr1	151	152	-	0	+
chr1	167	168	m	0	+
chr1	173	174	-	0	+
chr1	181	182	m	0	+
chr1	184	185	-	0	+
chr1	190	191	m	0	+
chr1	192	193	-	0	+
chr1	198	199	m	0	+
chr1	206	207	-	0	+
chr1	211	212	m	0	+
chr1	219	220	-	0	+
chr1	226	227	m	0	+
chr1	232	233	-	0	+
chr1	238	239	m	0	+
```

Ground-truth modified (`m`) and canonical (`-`) positions for `validate`.

**Example file:** `genes.gtf`

```text
chr1	test	transcript	1001	1280	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	exon	1001	1080	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	exon	1201	1280	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	transcript	1001	1360	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
chr1	test	exon	1001	1060	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
chr1	test	exon	1301	1360	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
```

Gene `g1` with the two transcripts `tx1` and `tx2` (exon lengths 160 and 120).

**Example code:**

```python
import shutil
from pathlib import Path

import pysam

from dnakit.epigenome import modkit_pileup

for name in ("reads", "other", "duplex", "tx_a", "tx_b"):
    pysam.sort("-o", f"{name}.bam", f"{name}.sam")
    pysam.index(f"{name}.bam")
pysam.faidx("ref.fa")
pysam.faidx("tx.fa")
shutil.rmtree("data", ignore_errors=True)
Path("data").mkdir()
for name, bam, reference in (("a", "reads", "ref"), ("b", "other", "ref"), ("tx_a", "tx_a", "tx"), ("tx_b", "tx_b", "tx")):
    result = modkit_pileup(f"{bam}.bam", f"data/{name}.bed.gz", reference=f"{reference}.fa", cpg=True,
                           modified_bases=["5mC"], bgzf=True, verbose=False)
    pysam.tabix_index(f"data/{name}.bed.gz", preset="bed", force=True, keep_original=True)
    print(name, result["status"], result["record_count"])
```

**Example output:**

```text
a completed 116
b completed 119
tx_a completed 77
tx_b completed 79
```

Preparation: sorted and indexed BAM files, and CpG bedMethyl pileups of the genome reads (`a`, `b`) and the transcript reads (`tx_a`, `tx_b`) as bgzip files with tabix indexes in `data`.

## 1) `dnakit.epigenome.modkit_bedmethyl_map_to_genome()`

- **Purpose:** Lift the records of a transcriptome-aligned bedMethyl file for one transcript onto genome coordinates.
- **Method:** Run `modkit bedmethyl map-to-genome <bedmethyl> <output> --transcript-id ... --gtf ...` (plus `--transcript-version`, `--ignore-version`, `--header` when requested). modkit finds the transcript's exons in the GTF and converts each transcript position (on the minus strand also the strand) to a genome position. The status is `completed` when the output exists.
- **API:** `modkit_bedmethyl_map_to_genome(bedmethyl_path: 'str', output_path: 'str', gtf_path: 'str', transcript_id: 'str', *, transcript_version: 'int | None' = None, ignore_version: 'bool' = False, header: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index whose contigs are transcript ids; `output_path` a new file; `gtf_path` the annotation; `transcript_id` the transcript to map; `transcript_version` and `ignore_version` handle version suffixes; `header` writes a header line.
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_bedmethyl_map_to_genome

work = Path(tempfile.mkdtemp())
result = modkit_bedmethyl_map_to_genome(
    "data/tx_a.bed.gz", str(work / "genome.bed"), "genes.gtf", "tx1", verbose=False
)
print(result["status"], result["record_count"])
first = (work / "genome.bed").read_text().splitlines()[0].split("\t")
print(first[0], first[1], first[2], first[5])
```

- **Example output:**

```text
completed 43
chr1 1010 1011 +
```

## 2) `dnakit.epigenome.modkit_bedmethyl_merge()`

- **Purpose:** Combine two or more bedMethyl files into one by an outer join that sums the counts of overlapping records.
- **Method:** Run `modkit bedmethyl merge <file> <file>... --genome-sizes ... --out-bed ...` (plus `--header`, `--min-samples`, `--min-sample-coverage`, `--threads`). Records at the same position, strand and modification code are added together; records seen in fewer than `min_samples` samples or with too little coverage are dropped.
- **API:** `modkit_bedmethyl_merge(bedmethyl_paths: 'Sequence[str]', genome_sizes: 'str', output_path: 'str', *, header: 'bool' = False, min_samples: 'int | None' = None, min_sample_coverage: 'int | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_paths` lists at least two bgzip bedMethyl files with tabix indexes; `genome_sizes` is a two-column chromosome sizes file; `output_path` a new file; `min_samples` and `min_sample_coverage` are optional integers.
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_bedmethyl_merge

work = Path(tempfile.mkdtemp())
result = modkit_bedmethyl_merge(
    ["data/a.bed.gz", "data/b.bed.gz"], "sizes.txt", str(work / "merged.bed"), verbose=False
)
print(result["status"], result["record_count"])
print((work / "merged.bed").read_text().splitlines()[0].split("\t")[:6])
```

- **Example output:**

```text
completed 122
['chr1', '2', '3', 'm', '2', '-']
```

## 3) `dnakit.epigenome.modkit_bedmethyl_to_bigwig()`

- **Purpose:** Turn a bedMethyl table into a BigWig track of the fraction modified.
- **Method:** Run `modkit bedmethyl tobigwig <bedmethyl> <output> --sizes ...` (or `--header <bam>`) with one `--mod-codes` per code, plus `--negative-strand-values`, `--nzooms` and `--uncompressed` when requested. With several codes modkit writes one track per code and appends the code to the file name.
- **API:** `modkit_bedmethyl_to_bigwig(bedmethyl_path: 'str', output_path: 'str', mod_codes: 'Sequence[str]', *, chrom_sizes: 'str | None' = None, header_bam: 'str | None' = None, negative_strand_values: 'bool' = False, nzooms: 'int | None' = None, uncompressed: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is an uncompressed bedMethyl file (a bgzip file is not accepted by modkit here); `output_path` a new file; `mod_codes` the modification codes (for example `['m']`); exactly one of `chrom_sizes` (two-column file) or `header_bam` supplies the chromosome sizes.
- **Output:** The common result dictionary.
- **Example code:**

```python
import gzip
import shutil
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_bedmethyl_to_bigwig

work = Path(tempfile.mkdtemp())
with gzip.open("data/a.bed.gz", "rb") as source, open(work / "a.bed", "wb") as target:
    shutil.copyfileobj(source, target)
result = modkit_bedmethyl_to_bigwig(
    str(work / "a.bed"), str(work / "frac.bw"), ["m"], chrom_sizes="sizes.txt", verbose=False
)
print(result["status"], Path(result["outputs"][0]).name, Path(result["outputs"][0]).stat().st_size > 0)
```

- **Example output:**

```text
completed frac.bw True
```

## 4) `dnakit.epigenome.modkit_dmr()`

- **Purpose:** Test regions for differential methylation between two samples (pair) or between all pairs of several samples (multi).
- **Method:** `mode='pair'` runs `modkit dmr pair -a <control> -b <experiment> --regions-bed ... --ref ... --out-path ...`; `mode='multi'` runs `modkit dmr multi --sample <file> <name> ...` for every sample and writes one BED per pair into a new directory. Both take `--base` (default C), optional `--assign-code`, `--single-code`, `--mask`, `--min-valid-coverage`, `--header`, and, in pair mode, `--segment`. modkit scores each region with a Bayesian model of the modified fraction; the score column is the magnitude of the difference.
- **API:** `modkit_dmr(samples: 'dict[str, str]', regions_bed: 'str', reference: 'str', output_path: 'str', *, mode: 'str' = 'pair', base: 'Sequence[str]' = ('C',), assign_code: 'Sequence[str] | None' = None, single_code: 'str | None' = None, mask: 'bool' = False, min_valid_coverage: 'int | None' = None, segment: 'str | None' = None, reference_sample: 'str | None' = None, prefix: 'str | None' = None, header: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `samples` maps a sample name to a bgzip bedMethyl with a tabix index (name: letters, digits, `_`, `.`, `-`); pair mode needs exactly two (the first is the control `-a`), multi mode at least two. `regions_bed` and `reference` (FASTA with .fai) are required. `output_path` is a new BED file (pair) or a new directory (multi, with an optional `prefix` and `reference_sample`).
- **Output:** The common result dictionary plus `mode` and `record_count` (pair) or `files` (multi).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_dmr

work = Path(tempfile.mkdtemp())
samples = {"a": "data/a.bed.gz", "b": "data/b.bed.gz"}
pair = modkit_dmr(samples, "regions.bed", "ref.fa", str(work / "dmr.bed"), verbose=False)
print(pair["status"], pair["mode"], pair["record_count"])
for line in (work / "dmr.bed").read_text().splitlines():
    fields = line.split("\t")
    print(fields[0], fields[1], fields[2], fields[4], fields[5])
multi = modkit_dmr(samples, "regions.bed", "ref.fa", str(work / "multi"), mode="multi", verbose=False)
print(multi["status"], [Path(p).name for p in multi["files"]])
```

- **Example output:**

```text
completed pair 3
chr1 0 120 25.70579668561382 .
chr1 120 240 25.18549683798696 .
chr2 0 160 22.296727753583554 .
completed ['a_b.bed']
```

## 5) `dnakit.epigenome.modkit_dmr_compare_tx_sites()`

- **Purpose:** Compare the modification level of every transcript site between two conditions.
- **Method:** Run `modkit dmr compare-tx-sites -a <A> -b <B> --out ... --gtf ...` (alias `compare-transcriptome-sites`) with the optional `--full`, `--ignore-version`, `--min-valid-coverage`, `--single-code`, `--plot`, `--top-k`, `--min-effect-size`, `--plot-title`. Each site is tested with a differential model and the table lists score and p-value with the gene of the transcript.
- **API:** `modkit_dmr_compare_tx_sites(condition_a: 'str', condition_b: 'str', gtf_path: 'str', output_path: 'str', *, min_valid_coverage: 'int | None' = None, single_code: 'str | None' = None, full: 'bool' = False, ignore_version: 'bool' = False, plot_path: 'str | None' = None, top_k: 'int | None' = None, min_effect_size: 'float | None' = None, plot_title: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `condition_a` and `condition_b` are transcriptome-aligned bgzip bedMethyl files with tabix indexes; `gtf_path` the annotation; `output_path` a new table; `plot_path` an optional new plot file (with `top_k`, `min_effect_size`, `plot_title`).
- **Output:** The common result dictionary plus `columns`, `row_count`, `rows` and `plot`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_dmr_compare_tx_sites

work = Path(tempfile.mkdtemp())
result = modkit_dmr_compare_tx_sites(
    "data/tx_a.bed.gz", "data/tx_b.bed.gz", "genes.gtf", str(work / "sites.tsv"), verbose=False
)
print(result["status"], result["row_count"])
print(result["columns"])
print(result["rows"][0][:4], result["rows"][0][-2:])
```

- **Example output:**

```text
completed 72
['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'p_value', 'gene_id', 'gene_name']
['chr1', '1001', '1002', 'm'] ['g1.0', 'GENE1']
```

## 6) `dnakit.epigenome.modkit_dmr_isoform()`

- **Purpose:** Compare methylation between the transcript isoforms of a gene from one transcriptome-aligned bedMethyl file.
- **Method:** Run `modkit dmr isoform <bedmethyl> <output> --gtf ...` with the optional `--gene-id`, `--gene-name`, `--gene-strand`, `--gene-chrom`, `--min-valid-coverage`, `--single-code`, `--full`, `--ignore-version`, `--plot`, `--min-score`, `--max-pval`. The sites that map to the same genome position in different isoforms are compared with a differential methylation model.
- **API:** `modkit_dmr_isoform(bedmethyl_path: 'str', gtf_path: 'str', output_path: 'str', *, gene_id: 'str | None' = None, gene_name: 'str | None' = None, gene_strand: 'str | None' = None, gene_chrom: 'str | None' = None, min_valid_coverage: 'int | None' = None, single_code: 'str | None' = None, full: 'bool' = False, ignore_version: 'bool' = False, plot_path: 'str | None' = None, min_score: 'float | None' = None, max_pval: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a transcriptome-aligned bgzip bedMethyl with a tabix index, `gtf_path` the annotation, `output_path` a new table; `gene_id`/`gene_name` select one gene (otherwise all genes are tested); `plot_path` is an optional new SVG plot (one gene).
- **Output:** The common result dictionary plus `columns`, `row_count`, `rows` and `plot`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_dmr_isoform

work = Path(tempfile.mkdtemp())
result = modkit_dmr_isoform("data/tx_a.bed.gz", "genes.gtf", str(work / "isoform.tsv"), verbose=False)
print(result["status"], result["row_count"])
print(result["columns"][:7], result["columns"][-3:])
```

- **Example output:**

```text
completed 4
['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'p_value'] ['n_transcripts', 'gene_id', 'gene_name']
```

## 7) `dnakit.epigenome.modkit_entropy()`

- **Purpose:** Compute methylation entropy over windows of consecutive modified sites.
- **Method:** Run `modkit entropy --in-bam <bam> --ref ... --out-bed ...` with the site selection (`--cpg`, `--motif` or `--base`), `--num-positions` (default 4 sites), `--window-size` (default 50 bp), `--min-coverage`, `--combine-strands`, `--regions` and the filtering options. For each window of `num_positions` modified sites the reads are grouped by their methylation pattern and the Shannon entropy of the pattern frequencies is reported; zero entropy means every read has the same pattern.
- **API:** `modkit_entropy(bam_paths: 'Sequence[str]', reference: 'str', output_path: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, base: 'str | None' = None, regions_bed: 'str | None' = None, num_positions: 'int | None' = None, window_size: 'int | None' = None, min_coverage: 'int | None' = None, max_filtered_positions: 'int | None' = None, combine_strands: 'bool' = False, mask: 'bool' = False, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, num_reads: 'int | None' = None, header: 'bool' = False, drop_zeros: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_paths` lists sorted, indexed modBAM files (pooled when several); `reference` is a FASTA with .fai; one of `cpg`, `motif` (+`motif_offset`) or `base` selects the sites; `output_path` is a new BED file, or with `regions_bed` a new directory that receives one BED per region (files start with `prefix`).
- **Output:** The common result dictionary plus `record_count` (BED mode) or `files` (directory mode).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_entropy

work = Path(tempfile.mkdtemp())
result = modkit_entropy(
    ["reads.bam"], "ref.fa", str(work / "entropy.bed"), cpg=True, num_positions=3, window_size=40,
    header=True, verbose=False,
)
print(result["status"], result["record_count"])
lines = (work / "entropy.bed").read_text().splitlines()
print(lines[0].split("\t")[:6])
print(lines[1].split("\t")[:6])
```

- **Example output:**

```text
completed 24
['#chrom', 'start', 'end', 'entropy', 'strand', 'num_reads']
['chr1', '42', '55', '0.3060986', '+', '3']
```

## 8) `dnakit.epigenome.modkit_extract()`

- **Purpose:** Extract read-level modification information from a modBAM as a table: probabilities, hard calls or per-read counts.
- **Method:** `table='full'` runs `modkit extract full <bam> <output>` (probability of every modification call per read position); `table='calls'` runs `modkit extract calls` (calls made with the same pass threshold as the pileup; `pass_only` keeps only passing calls); `table='read-stats'` runs `modkit extract read-stats` (one CSV row per read with modified and unmodified counts per code, `--mod-codes` required). Options that belong to another table raise `ConfigurationError` instead of being ignored.
- **API:** `modkit_extract(bam_path: 'str', output_path: 'str', *, table: 'str' = 'full', reference: 'str | None' = None, region: 'str | None' = None, include_bed: 'str | None' = None, exclude_bed: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, annotate_motifs: 'bool' = False, mask: 'bool' = False, ignore: 'str | None' = None, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, mapped_only: 'bool' = False, allow_non_primary: 'bool' = False, num_reads: 'int | None' = None, kmer_size: 'int | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, pass_only: 'bool' = False, sampling_frac: 'float | None' = None, mod_codes: 'Sequence[str] | None' = None, ignore_implicit: 'bool' = False, ignore_index: 'bool' = False, bgzf: 'bool' = False, no_headers: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM; `output_path` a new file; `table` is `full`, `calls` or `read-stats`; `reference` (needed for `cpg`, `motif`, `annotate_motifs`), `region`, `include_bed`, `exclude_bed`, `mapped_only`, `edge_filter`, `kmer_size`, `bgzf`, `no_headers` select and format rows; `filter_percentile`, `filter_threshold`, `mod_thresholds`, `no_filtering`, `pass_only` apply to `calls`; `mod_codes` such as `['C:m']` is required for `read-stats`.
- **Output:** The common result dictionary plus `table`, `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_extract

work = Path(tempfile.mkdtemp())
for table, extra in (("full", {}), ("calls", {"pass_only": True}), ("read-stats", {"mod_codes": ["C:m"]})):
    result = modkit_extract("reads.bam", str(work / f"{table}.tsv"), table=table, verbose=False, **extra)
    print(table, result["status"], result["row_count"], result["columns"][:3])
```

- **Example output:**

```text
full completed 205 ['read_id', 'forward_read_position', 'ref_position']
calls completed 185 ['read_id', 'forward_read_position', 'ref_position']
read-stats completed 12 ['read_id', 'chrom', 'aln_start']
```

## 9) `dnakit.epigenome.modkit_localize()`

- **Purpose:** Aggregate modification levels in windows centred on genomic features, such as promoters or binding sites.
- **Method:** Run `modkit localize <bedmethyl> --regions ... --genome-sizes ... --out-file ...` (alias `localise`) with `--window` (default 2000 bp on each side), `--stranded`, `--stranded-features`, `--min-coverage`, `--chart` and `--name`. Window positions are measured from the centre of each feature and the counts of all features are added up per position.
- **API:** `modkit_localize(bedmethyl_path: 'str', regions_bed: 'str', genome_sizes: 'str', output_path: 'str', *, window: 'int | None' = None, stranded: 'str | None' = None, stranded_features: 'bool' = False, min_coverage: 'int | None' = None, chart_path: 'str | None' = None, chart_name: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index; `regions_bed` the features; `genome_sizes` a two-column chromosome sizes file; `output_path` a new table; `window` (bp on each side), `stranded` (`true`, `false` or `both`), `min_coverage`; `chart_path` an optional new HTML chart.
- **Output:** The common result dictionary plus `columns`, `row_count`, `rows` and `chart`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_localize

work = Path(tempfile.mkdtemp())
result = modkit_localize(
    "data/a.bed.gz", "features.bed", "sizes.txt", str(work / "localize.tsv"), window=30,
    chart_path=str(work / "localize.html"), verbose=False,
)
print(result["status"], result["row_count"], Path(result["chart"]).name)
print(result["columns"][:6])
```

- **Example output:**

```text
completed 20 localize.html
['mod_code', 'offset', 'n_valid', 'n_mod', 'percent_modified']
```

## 10) `dnakit.epigenome.modkit_modbam_adjust_mods()`

- **Purpose:** Rewrite the base-modification calls of a modBAM: drop a code, convert one code to another, or keep or discard sequence contexts.
- **Method:** Run `modkit modbam adjust-mods <in> <out>` with `--ignore`, `--convert <from> <to>`, `--cpg` or `--motif` (with `--discard-motifs`), `--edge-filter`, `--invert-edge-filter`, `--filter-probs`, `--mapped-only`, `--num-reads`, the pass-threshold options and `--output-sam`. Converting `h` to `m` collapses 5hmC calls into 5mC; ignoring a code removes its calls and renormalises the probabilities.
- **API:** `modkit_modbam_adjust_mods(bam_path: 'str', output_bam: 'str', *, ignore: 'str | None' = None, convert_from: 'str | None' = None, convert_to: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, discard_motifs: 'bool' = False, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, filter_probs: 'bool' = False, mapped_only: 'bool' = False, num_reads: 'int | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, output_sam: 'bool' = False, fast_forward: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM and `output_bam` a new file; `ignore` is one modification code; `convert_from` and `convert_to` are given together; `cpg` / `motif` (+`motif_offset`) and `discard_motifs` select contexts (a reference is not needed); `edge_filter` is `N` or `start,end`; `output_sam` writes SAM instead of BAM.
- **Output:** The common result dictionary.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_adjust_mods

work = Path(tempfile.mkdtemp())
import pysam

result = modkit_modbam_adjust_mods(
    "reads.bam", str(work / "converted.bam"), convert_from="m", convert_to="h", verbose=False
)
print(result["status"])
with pysam.AlignmentFile("reads.bam") as before, pysam.AlignmentFile(str(work / "converted.bam")) as after:
    print(next(iter(before)).get_tag("MM")[:8], "->", next(iter(after)).get_tag("MM")[:8])
```

- **Example output:**

```text
completed
C+m?,0,0 -> C+h?,0,0
```

## 11) `dnakit.epigenome.modkit_modbam_call_mods()`

- **Purpose:** Threshold the modification probabilities of a modBAM and write a new modBAM of hard calls.
- **Method:** Run `modkit modbam call-mods <in> <out>` with the pass-threshold options (`--filter-percentile` default 10, `--filter-threshold`, `--mod-threshold`, `--no-filtering`), sampling (`--num-reads`, `--sampling-frac`, `--seed`), `--edge-filter`, `--cpg` / `--motif` (+`--discard-motifs`) and `--output-sam`. Every base is set to probability 100 percent when a modification is called and to 0 percent when it is called canonical.
- **API:** `modkit_modbam_call_mods(bam_path: 'str', output_bam: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, discard_motifs: 'bool' = False, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, output_sam: 'bool' = False, fast_forward: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM and `output_bam` a new file; `filter_threshold` and `filter_percentile` set the pass threshold; `mod_thresholds` entries look like `m:0.8`; `output_sam` writes SAM.
- **Output:** The common result dictionary.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_call_mods

work = Path(tempfile.mkdtemp())
import pysam

result = modkit_modbam_call_mods("reads.bam", str(work / "called.bam"), filter_threshold=0.8, verbose=False)
print(result["status"])
with pysam.AlignmentFile(str(work / "called.bam")) as bam:
    record = next(iter(bam))
    print(sorted(set(record.get_tag("ML"))))
```

- **Example output:**

```text
completed
[0, 255]
```

## 12) `dnakit.epigenome.modkit_modbam_check_tags()`

- **Purpose:** Validate the MM and ML modification tags of a BAM.
- **Method:** Run `modkit modbam check-tags <bam> --out-dir <new directory>` with `--permissive`, `--region`, `--head`, `--mapped-only`, `--allow-non-primary`, `--ignore-index`, `--prefix` and `--threads`. modkit checks that the number of ML values matches the MM entries of each read and writes tables of the modifications and MM headers it found plus any failure tables. Without `permissive` it exits with code 1 when invalid records are found, which gives the status `failed`.
- **API:** `modkit_modbam_check_tags(bam_path: 'str', output_dir: 'str', *, permissive: 'bool' = False, region: 'str | None' = None, num_reads: 'int | None' = None, head: 'int | None' = None, mapped_only: 'bool' = False, allow_non_primary: 'bool' = False, ignore_index: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM; `output_dir` is a new directory for the tables; `permissive`, `region`, `head` (number of reads), `mapped_only`, `allow_non_primary` and `prefix` are optional.
- **Output:** The common result dictionary plus `files`, the tables that modkit wrote.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_check_tags

work = Path(tempfile.mkdtemp())
result = modkit_modbam_check_tags("reads.bam", str(work / "checks"), verbose=False)
print(result["status"], [Path(p).name for p in result["files"]])
print((work / "checks" / "modified_bases.tsv").read_text().splitlines())
```

- **Example output:**

```text
completed ['modified_bases.tsv', 'valid_mm_headers.tsv']
['strand\tprimary_base\tmod_code\tmode', '+\tC\tm\t?']
```

## 13) `dnakit.epigenome.modkit_modbam_repair()`

- **Purpose:** Copy correct MM/ML tags from a donor modBAM onto the reads of an acceptor BAM.
- **Method:** Run `modkit modbam repair --donor-bam <donor> --acceptor-bam <acceptor> --output-bam <out>`. Both BAMs must be sorted by read name and the donor must contain every acceptor read (extra donor reads are allowed). Typical use: restore tags that a trimming or alignment step removed.
- **API:** `modkit_modbam_repair(donor_bam: 'str', acceptor_bam: 'str', output_bam: 'str', *, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `donor_bam` carries the tags, `acceptor_bam` the reads that lack them, `output_bam` is a new file; `threads` is optional.
- **Output:** The common result dictionary.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_repair

work = Path(tempfile.mkdtemp())
import pysam

pysam.sort("-n", "-o", str(work / "donor.bam"), "reads.sam")
with pysam.AlignmentFile("reads.sam") as source, pysam.AlignmentFile(str(work / "stripped.unsorted.bam"), "wb", template=source) as target:
    for record in source:
        record.set_tag("MM", None)
        record.set_tag("ML", None)
        target.write(record)
pysam.sort("-n", "-o", str(work / "stripped.bam"), str(work / "stripped.unsorted.bam"))
result = modkit_modbam_repair(str(work / "donor.bam"), str(work / "stripped.bam"), str(work / "repaired.bam"), verbose=False)
print(result["status"])
with pysam.AlignmentFile(str(work / "repaired.bam"), check_sq=False) as bam:
    print(sum(1 for record in bam if record.has_tag("MM")), "reads have MM again")
```

- **Example output:**

```text
completed
12 reads have MM again
```

## 14) `dnakit.epigenome.modkit_modbam_sample_probs()`

- **Purpose:** Estimate the distribution of base-modification probabilities in a modBAM and the pass thresholds derived from it.
- **Method:** Run `modkit modbam sample-probs <bam> --out-dir <new directory>` with `--percentiles` (default 10,50,90), `--quantiles`, `--hist`, `--reference`, `--motif`, `--region`, `--include-bed`, `--edge-filter`, `--use-non-primary`, sampling options and `--prefix`. modkit samples reads (10,042 by default), writes `probabilities.tsv` and `thresholds.tsv` and, with `--hist`, HTML histograms.
- **API:** `modkit_modbam_sample_probs(bam_path: 'str', output_dir: 'str', *, percentiles: 'str | None' = None, quantiles: 'str | None' = None, histogram: 'bool' = False, reference: 'str | None' = None, motif: 'str | None' = None, motif_offset: 'int' = 0, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, use_non_primary: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, no_sampling: 'bool' = False, seed: 'int | None' = None, ignore_index: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM and `output_dir` a new directory; `percentiles` / `quantiles` are comma-separated numbers; `histogram` adds the HTML histograms; `no_sampling` reads every read.
- **Output:** The common result dictionary plus `files`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_sample_probs

work = Path(tempfile.mkdtemp())
result = modkit_modbam_sample_probs("reads.bam", str(work / "probs"), percentiles="10,50,90", verbose=False)
print(result["status"], [Path(p).name for p in result["files"]])
for line in (work / "probs" / "thresholds.tsv").read_text().splitlines():
    print(line.rstrip())
```

- **Example output:**

```text
completed ['thresholds.tsv']
 base  percentile  threshold  mod_qual
 C     10          0.8066406  206
 C     50          0.8925781  228
 C     90          0.9824219  251
```

## 15) `dnakit.epigenome.modkit_modbam_summary()`

- **Purpose:** Summarise the modification tags present in a modBAM with basic call statistics.
- **Method:** Run `modkit modbam summary <bam>` with `--tsv`, `--combine-mods`, `--cpg` / `--motif` (with `--reference`), `--region`, `--include-bed`, `--edge-filter`, `--mapped-only`, `--matched-only`, sampling and pass-threshold options. The default output is a totals block (lines starting with `#`) followed by a table of calls per base and code; `--tsv` writes one tab-separated table. Standard output is streamed into the new file.
- **API:** `modkit_modbam_summary(bam_path: 'str', output_path: 'str', *, tsv: 'bool' = False, combine_mods: 'bool' = False, cpg: 'bool' = False, reference: 'str | None' = None, motif: 'str | None' = None, motif_offset: 'int' = 0, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, mapped_only: 'bool' = False, matched_only: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, filter_quantile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, ignore_index: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM and `output_path` a new file; `tsv`, `combine_mods`, `cpg`/`motif` (needs `reference`), `region`, `edge_filter`, `num_reads` (default sampling 10,042 reads), `filter_threshold`, `mod_thresholds` are optional.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_summary

work = Path(tempfile.mkdtemp())
result = modkit_modbam_summary("reads.bam", str(work / "summary.txt"), verbose=False)
print(result["status"])
for line in (work / "summary.txt").read_text().splitlines()[:8]:
    print(line.rstrip())
```

- **Example output:**

```text
completed
# bases                     C
# total_reads_used          12
# count_reads_C             12
# pass_threshold_C          0.8066406
# modification_codes_for_C  m
 base  code  pass_count  pass_frac   all_count  all_frac
 C     -     45          0.24324325  45         0.2195122
 C     m     140         0.7567568   160        0.7804878
```

## 16) `dnakit.epigenome.modkit_modbam_update_tags()`

- **Purpose:** Rename lower-case `Mm`/`Ml` tags to `MM`/`ML` and set the mode flag of the modification tags.
- **Method:** Run `modkit modbam update-tags <in> <out>` with `--mode explicit|implicit`, `--no-implicit-probs`, `--output-sam` and `--threads`. `explicit` writes the `?` flag, `implicit` the `.` flag (unlisted bases are canonical); modkit by default refuses implicit-mode tags in a pileup, so this step makes older BAMs usable.
- **API:** `modkit_modbam_update_tags(bam_path: 'str', output_bam: 'str', *, mode: 'str | None' = None, no_implicit_probs: 'bool' = False, output_sam: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a BAM with modification tags and `output_bam` a new file; `mode` is `explicit` or `implicit`; `no_implicit_probs` skips adding implicit canonical calls; `output_sam` writes SAM.
- **Output:** The common result dictionary.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_update_tags

work = Path(tempfile.mkdtemp())
import pysam

result = modkit_modbam_update_tags("reads.bam", str(work / "implicit.bam"), mode="implicit", verbose=False)
print(result["status"])
with pysam.AlignmentFile(str(work / "implicit.bam")) as bam:
    print(next(iter(bam)).get_tag("MM")[:8])
```

- **Example output:**

```text
completed
C+m.,0,0
```

## 17) `dnakit.epigenome.modkit_motif_bed()`

- **Purpose:** Write a BED file with every location of a sequence motif in a FASTA.
- **Method:** Run `modkit motif bed <fasta> <motif> <offset>` (and `--mask`) and stream its standard output into the new file. The motif may use IUPAC codes; `offset` is the 0-based position inside the motif that the interval marks, so `CG` with offset 0 marks the C of each CpG; the same site is written a second time with strand `-` at the position of the G (the cytosine of the opposite strand).
- **API:** `modkit_motif_bed(fasta_path: 'str', motif: 'str', offset: 'int', output_path: 'str', *, mask: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `fasta_path` is a FASTA, `motif` letters only (IUPAC allowed), `offset` an integer smaller than the motif length, `output_path` a new file; `mask=True` respects soft-masking.
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_motif_bed

work = Path(tempfile.mkdtemp())
result = modkit_motif_bed("ref.fa", "CG", 0, str(work / "cpg.bed"), verbose=False)
print(result["status"], result["record_count"])
lines = (work / "cpg.bed").read_text().splitlines()
print(lines[0])
print(lines[1])
```

- **Example output:**

```text
completed 138
chr1	1	2	.	.	+
chr1	2	3	.	.	-
```

## 18) `dnakit.epigenome.modkit_motif_evaluate()`

- **Purpose:** Calculate enrichment statistics of given motifs from a bedMethyl table.
- **Method:** Run `modkit motif evaluate --in-bedmethyl ... --ref ... --out ...` with `--known-motif <sequence> <offset> <mod_code>` (repeatable), `--known-motifs-table`, `--contig`, `--min-coverage`, `--context-size` and the `--low-thresh` / `--high-thresh` fraction-modified cut-offs. For each motif modkit counts the genome sites that match it with a high, middle or low fraction modified and reports the log-odds of enrichment.
- **API:** `modkit_motif_evaluate(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, contig: 'str | None' = None, context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index, `reference` a FASTA with .fai, `output_table` a new TSV; at least one of `known_motifs` (`[sequence, offset, mod_code]` triples) or `known_motifs_table` is required.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_motif_evaluate

work = Path(tempfile.mkdtemp())
result = modkit_motif_evaluate(
    "data/a.bed.gz", "ref.fa", str(work / "evaluate.tsv"), known_motifs=[["CG", "0", "m"]], min_coverage=1,
    verbose=False,
)
print(result["status"], result["columns"])
print(result["rows"][0])
```

- **Example output:**

```text
completed ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count', 'log_odds']
['m', 'CG', '0', '0.85057473', '74', '13', '23', '0']
```

## 19) `dnakit.epigenome.modkit_motif_refine()`

- **Purpose:** Refine a list of known motifs with a bedMethyl table.
- **Method:** Run `modkit motif refine --in-bedmethyl ... --ref ... --out ...` with `--known-motif`, `--known-motifs-table`, `--min-refine-frac-mod`, `--min-refine-sites`, `--context-size`, `--min-coverage`, `--min-sites`, `--min-log-odds`, `--min-frac-mod` and the fraction-modified thresholds. Each known motif is extended or narrowed with the sequence context that raises its enrichment.
- **API:** `modkit_motif_refine(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, contig: 'str | None' = None, context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, min_sites: 'int | None' = None, min_log_odds: 'float | None' = None, min_frac_mod: 'float | None' = None, min_refine_frac_mod: 'float | None' = None, min_refine_sites: 'int | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index, `reference` a FASTA with .fai, `output_table` a new TSV; at least one of `known_motifs` or `known_motifs_table` is required.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_motif_refine

work = Path(tempfile.mkdtemp())
result = modkit_motif_refine(
    "data/a.bed.gz", "ref.fa", str(work / "refine.tsv"), known_motifs=[["C", "0", "m"]], min_coverage=1,
    min_sites=3, min_refine_sites=3, verbose=False,
)
print(result["status"], result["columns"])
print(result["rows"][0])
```

- **Example output:**

```text
completed ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count', 'log_odds']
['m', 'C', '0', '0.85057473', '74', '13', '23', '0']
```

## 20) `dnakit.epigenome.modkit_motif_search()`

- **Purpose:** Search a genome for modification-enriched sequence motifs from a bedMethyl table.
- **Method:** Run `modkit motif search --in-bedmethyl ... --ref ... --out-table ...` with `--contig`, `--mod-code`, `--known-motif`, `--known-motifs-table`, `--context-size`, `--init-context-size`, `--min-coverage`, `--min-sites`, `--min-log-odds`, `--min-frac-mod` and the thresholds. modkit starts from short seeds whose sites are mostly modified, grows them into longer motifs and keeps those with enough sites and log-odds.
- **API:** `modkit_motif_search(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, contig: 'str | None' = None, mod_codes: 'Sequence[str] | None' = None, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, context_size: 'Sequence[int] | None' = None, init_context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, min_sites: 'int | None' = None, min_log_odds: 'float | None' = None, min_frac_mod: 'float | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index, `reference` a FASTA with .fai, `output_table` a new TSV; `mod_codes` restricts the codes; `known_motifs` are `[sequence, offset, mod_code]` triples to report next to the discovered ones.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_motif_search

work = Path(tempfile.mkdtemp())
result = modkit_motif_search(
    "data/a.bed.gz", "ref.fa", str(work / "search.tsv"), min_coverage=1, min_sites=3, verbose=False
)
print(result["status"], result["row_count"], result["columns"])
print(result["rows"][0][:4])
```

- **Example output:**

```text
completed 1 ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count']
['m', 'C', '0', '0.85057473']
```

## 21) `dnakit.epigenome.modkit_open_chromatin_predict()`

- **Purpose:** Predict open-chromatin regions from exogenous 6mA signal in a modBAM (Fiber-seq style data).
- **Method:** Run `modkit open-chromatin predict <bam> --model <model directory> --output <bed>` with `--region`, `--include-bed`, `--min-coverage`, `--threshold`, `--batch-size`, `--step-size`, `--no-header` and `--dryrun`. The model is a directory that contains `model_config.json` and the weights (the modkit release archive ships one under `models/`). With `dry_run=True` modkit loads the model, checks the inputs and exits without prediction, leaving an empty output file.
- **API:** `modkit_open_chromatin_predict(bam_path: 'str', model_path: 'str', output_path: 'str', *, region: 'str | None' = None, include_bed: 'str | None' = None, min_coverage: 'int | None' = None, threshold: 'float | None' = None, batch_size: 'int | None' = None, step_size: 'int | None' = None, no_header: 'bool' = False, dry_run: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a modBAM with 6mA calls, `model_path` an existing model directory, `output_path` a new BED file; `threshold` is the score cut-off, `min_coverage` the minimum read coverage; `dry_run` only checks inputs and the model.
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_open_chromatin_predict

work = Path(tempfile.mkdtemp())
import shutil

models = Path(shutil.which("modkit")).resolve().parent / "models"
model = sorted(models.iterdir())[0]
result = modkit_open_chromatin_predict("reads.bam", str(model), str(work / "open.bed"), dry_run=True, verbose=False)
print(result["status"], result["record_count"])
```

- **Example output:**

```text
completed 0
```

## 22) `dnakit.epigenome.modkit_pileup()`

- **Purpose:** Tabulate base-modification calls per genomic position into a bedMethyl file.
- **Method:** Run `modkit pileup <bam> <output>` with `--cpg` or `--motif` (with `--reference`), `--modified-bases`, `--combine-mods`, `--combine-strands`, `--region`, `--include-bed`, `--edge-filter`, `--max-depth`, sampling options, the pass-threshold options (`--filter-percentile` default 10, `--filter-threshold`, `--mod-threshold`, `--no-filtering`), `--duplex`, `--header` and `--bgzf`. A call counts as modified or canonical only above the pass threshold, otherwise it is counted as filtered; the 18-column bedMethyl has the coverage, the percent modified and these counts per position, strand and code.
- **API:** `modkit_pileup(bam_path: 'str', output_path: 'str', *, modified_bases: 'Sequence[str] | None' = None, reference: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, combine_mods: 'bool' = False, combine_strands: 'bool' = False, mask: 'bool' = False, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, max_depth: 'int | None' = None, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, duplex: 'bool' = False, header: 'bool' = False, bgzf: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a sorted, indexed modBAM, `output_path` a new file; `cpg` or `motif` (+`motif_offset`) need `reference` (FASTA with .fai); `modified_bases` such as `['5mC']` or `['C:m']`; `bgzf=True` writes bgzip output that you index with tabix before using it as input of `modkit_stats`, `modkit_localize`, `modkit_bedmethyl_merge` or `modkit_dmr`; unset options keep modkit's defaults.
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_pileup

work = Path(tempfile.mkdtemp())
result = modkit_pileup(
    "reads.bam", str(work / "cpg.bed"), reference="ref.fa", cpg=True, modified_bases=["5mC"], verbose=False
)
print(result["status"], result["record_count"])
first = (work / "cpg.bed").read_text().splitlines()[0].split("\t")
print(first[:6], first[9:12])
```

- **Example output:**

```text
completed 116
['chr1', '2', '3', 'm', '1', '-'] ['1', '0.00', '0']
```

## 23) `dnakit.epigenome.modkit_pileup_hemi()`

- **Purpose:** Tabulate double-stranded modification patterns (such as hemi-methylation) across palindromic motif positions.
- **Method:** Run `modkit pileup-hemi <bam> --ref ... --out-bed ...` with `--cpg` (shorthand for motif `CG` offset 0) or `--motif`, plus `--combine-mods`, `--ignore`, `--region`, `--include-bed`, `--edge-filter`, sampling and pass-threshold options. Each molecule must report calls on both strands (duplex reads with `C+m?` and `G-m?`); the output lists for every position the pattern counts such as `m,m,C`, `m,-,C` and `-,m,C` (first base, second base, primary base).
- **API:** `modkit_pileup_hemi(bam_path: 'str', reference: 'str', output_path: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, combine_mods: 'bool' = False, mask: 'bool' = False, ignore: 'str | None' = None, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, max_depth: 'int | None' = None, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_path` is a sorted, indexed duplex modBAM, `reference` a FASTA with .fai, `output_path` a new file; exactly one of `cpg=True` or `motif` (+`motif_offset`, reverse-complement palindromic).
- **Output:** The common result dictionary plus `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_pileup_hemi

work = Path(tempfile.mkdtemp())
result = modkit_pileup_hemi("duplex.bam", "ref.fa", str(work / "hemi.bed"), cpg=True, verbose=False)
print(result["status"], result["record_count"])
patterns = sorted({line.split("\t")[3] for line in (work / "hemi.bed").read_text().splitlines()})
print(patterns)
```

- **Example output:**

```text
completed 133
['-,-,C', '-,m,C', 'm,-,C', 'm,m,C']
```

## 24) `dnakit.epigenome.modkit_stats()`

- **Purpose:** Report the modification level (fraction modified) over each region of a BED file.
- **Method:** Run `modkit stats <bedmethyl> --regions ... --out-table ...` with `--mod-codes`, `--min-coverage` (default 1), `--no-header` and `--threads`. The valid coverage and the modified counts of all bedMethyl records inside a region are added up per modification code and divided.
- **API:** `modkit_stats(bedmethyl_path: 'str', regions_bed: 'str', output_path: 'str', *, mod_codes: 'str | None' = None, min_coverage: 'int | None' = None, no_header: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bedmethyl_path` is a bgzip bedMethyl with a tabix index, `regions_bed` the regions, `output_path` a new table; `mod_codes` is a comma-separated string such as `'m,h'`.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_stats

work = Path(tempfile.mkdtemp())
result = modkit_stats("data/a.bed.gz", "regions.bed", str(work / "stats.tsv"), verbose=False)
print(result["status"], result["row_count"])
print(result["columns"])
for row in result["rows"]:
    print(row)
```

- **Example output:**

```text
completed 3
['#chrom', 'start', 'end', 'name', 'strand', 'count_m', 'count_valid_m', 'percent_m']
['chr1', '0', '120', 'reg1', '.', '49', '65', '75.38461']
['chr1', '120', '240', 'reg2', '.', '46', '57', '80.70176']
['chr2', '0', '160', 'reg3', '.', '45', '63', '71.42857']
```

## 25) `dnakit.epigenome.modkit_validate()`

- **Purpose:** Validate modified-base calls against ground-truth positions and report the confusion between calls and truth.
- **Method:** Run `modkit validate --bam-and-bed <bam> <bed> ... --out-filepath ...` with `--canonical-base`, `--ignore`, `--edge-filter`, `--invert-edge-filter`, `--min-identity`, `--min-length`, `--filter-quantile` and `--filter-threshold`. The name column of each truth BED holds the modification code, or `-` for a canonical position. modkit needs the `NM` tag in the BAM, builds contingency tables and writes a machine-readable summary.
- **API:** `modkit_validate(bam_and_bed: 'Sequence[Sequence[str]]', output_path: 'str', *, canonical_base: 'str | None' = None, ignore: 'str | None' = None, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, min_alignment_identity: 'float | None' = None, min_alignment_length: 'int | None' = None, filter_quantile: 'float | None' = None, filter_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the modkit program (a release binary, or `conda install -c bioconda modkit`). DNAKit looks for it at `modkit_path`, then on PATH; if it is missing, `BackendUnavailableError` is raised before the program is started and before anything is created.
- **Input:** `bam_and_bed` is a list of `[modBAM, truth BED]` pairs (one per sample), `output_path` a new file; `canonical_base` (A, C, G or T) is needed when the truth has only canonical sites or ChEBI codes.
- **Output:** The common result dictionary plus `columns`, `row_count` and `rows` (the first lines of the summary).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_validate

work = Path(tempfile.mkdtemp())
result = modkit_validate([["reads.bam", "truth.bed"]], str(work / "validate.tsv"), verbose=False)
print(result["status"], result["row_count"])
for line in (work / "validate.tsv").read_text().splitlines()[:3]:
    print(line[:90])
```

- **Example output:**

```text
completed 6
full_contingency_table: [["ground_truth_label","C","m","No Call"],["C",8,26,1],["m",11,22,
raw_accuracy: 45.454544
raw_contingency_table: [["ground_truth_label","C","m"],["C",8,25],["m",11,22]]
```
