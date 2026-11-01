# Picard SAM/BAM, Illumina-run, VCF and interval tools

This group of functions in `dnakit.alignment and dnakit.variants` is a set of command-line facades over [Picard](https://github.com/broadinstitute/picard) (Broad Institute, MIT licence) tools that read or rewrite SAM/BAM/CRAM files, Illumina run folders, VCF/BCF files, interval lists, FASTA and GFF3 files. Each function runs one fixed Picard tool with an explicit `KEY=value` argument list (no shell), refuses to overwrite an existing output, and returns a dictionary with `status`, `backend`, `tool`, `version`, `command`, `return_code`, `outputs` and `log_tail` plus a few parsed facts. `status` is `completed` only when the exit code is accepted and every expected output exists. Arguments left as `None` are not passed, so Picard's own defaults apply. Tested with Picard 3.5.0. Picard's JavaScript filters (`FilterSamReads includeJavascript`, `FilterVcf JAVASCRIPT_FILE`) are not exposed because they would run caller-supplied scripts.

Picard is an external program. Install it (`conda install -c bioconda picard`, or `picard.jar` plus Java 17) and put `picard` on PATH, pass `picard_path` (a launcher or a `picard.jar`) or set `PICARD_JAR`; if it is missing the function raises `BackendUnavailableError` with an install hint. The examples import from the modules directly and call every function with `verbose=False`.

## Example data

The examples read the small files below. Every output is written to a fresh temporary directory because the functions never overwrite a file.

**Example file:** `example.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1200
@SQ	SN:chr2	LN:800
@RG	ID:rg1	SM:sample1	LB:lib1	PL:ILLUMINA	PU:FC1.1
M1:7:FC1:1:1101:1481:2143	99	chr1	125	60	36M	=	267	178	ATTTTTATTACACTCAGAAACAGAACTCGGGTAATT	BH?AFHBIEAFEAFGI?IFGG?DH?E@FFFI?HIF?	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1407:2121	99	chr1	201	60	36M	=	382	217	GACACTCGCTATGAATCTCTGATTTACCCACTCTGC	H@?IBIGCD@HCCAC@FCBA55?GHID5AEDEI@BG	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1481:2143	147	chr1	267	60	36M	=	125	-178	GTAGTCAATAGAGCGAACGCATTATTCGGTTACTTA	EHDF@AGE?@CG5DA5C55IGAEB?I@D?GAG?IC@	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1407:2121	147	chr1	382	60	36M	=	201	-217	TGACTTGTCAGGCGGCAAATAGAGTTTCTCCGACGC	GH?@IA@@EHD?5EEAAIC55GFD@B?H5FIDC?E5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1333:2099	99	chr1	445	60	36M	=	615	206	GCCAAGACTATAGGCACTGTCGCATCACAAACGATT	CIFIIHHEF5HIADHF@EA5IFB@F@HAFBA5@CCD	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1000:2000	99	chr1	516	60	36M	=	704	224	GACTGGTTTACGATAGTATGTCCAACGGCGAGCTTT	FB???D@FGAA@HGEID@5HDIDGGF5D5CCDACID	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1001:2001	99	chr1	516	60	36M	=	704	224	GACTGGTTTACGATAGTATGTCCAACGGCGAGCTTT	FB???D@FGAA@HGEID@5HDIDGGF5D5CCDACID	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1185:2055	99	chr1	537	60	36M	=	669	168	CCAACGGCGAGCTTTACATTTGCTGTGAGAGGTACA	GADG?@@FF?5??@FEEGD5H5HGC@IACB@5BH?G	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1148:2044	99	chr1	598	60	36M	=	731	169	CAATTCGTACCTTGGGGGTCGTTACCACTCTGTTCC	ADEDEGGCIFI?@CCC?BF@?HBICFDH@FBFAFAD	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1222:2066	99	chr1	609	60	36M	=	751	178	TTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCA	DBDIE5A?@@CD@5BDFC?CFDCDH??DCFADAEBC	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1074:2022	99	chr1	610	60	36M	=	774	200	TGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCAT	5C@ACFCEBBAB55@F?GCEH5FDECI?FAHI@DCH	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1333:2099	147	chr1	615	60	36M	=	445	-206	CAGAAATGCCGCTCGTGGGAACAGAGTGGTAACGAC	?AHB@@HIEHEAIA5FIE@HCIB@I@GGACH?FD@H	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1037:2011	99	chr1	646	60	36M	=	794	184	TTCTGGATGGCCAGCTTTTGACATTTAATTTCACCC	HGBHF5G?B?F5DA5B?BCH@?5GFB?EGF@E?F@B	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1185:2055	147	chr1	669	60	36M	=	537	-168	CTTGCAGCTTTACGCTGGTTTATGGGTGAAATTAAA	?CAEGD55ADG5E5GAAA5@G@C5EBDGBE?AHDHI	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1370:2110	99	chr1	693	60	36M	=	826	169	TAAAGCTGCAAGTGGCTCCATGAACTTAGCTGCTAG	GEDA?IB5CEA55BBA?IBE?@CEEGCB@F?55EE?	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1444:2132	99	chr1	693	60	36M	=	829	172	TAAAGCTGCAAGTGGCTCCATGAACTTAGCTGCTAG	I@DB5EGHCGAE?FCFEDFH@DGG?5IHCGHBGGDC	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1259:2077	99	chr1	702	60	36M	=	832	166	AAGTGGCTCCATGAACTTAGCTGCTAGTGTCAGACT	IBGEGCFADGFA@?HF?FBID5HIG@B5DI?I@ACA	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1000:2000	147	chr1	704	60	36M	=	516	-224	CGAGTCTGACACTAGCAGCTAAGTTCATGGAGCCAC	GCDF5CF@HCADHH5C?F@?CDAFH5A@DDEH555H	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1001:2001	147	chr1	704	60	36M	=	516	-224	CGAGTCTGACACTAGCAGCTAAGTTCATGGAGCCAC	GCDF5CF@HCADHH5C?F@?CDAFH5A@DDEH555H	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1111:2033	99	chr1	720	60	36M	=	892	208	AGCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTA	BICBCDFFGDHC5IEDEB@FB@DGDGA?CCGACAD5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1004:2004	99	chr1	720	60	36M	=	892	208	AGCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTA	BICBCDFFGDHC5IEDEB@FB@DGDGA?CCGACAD5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1296:2088	99	chr1	727	60	36M	=	900	209	AGTGTCAGACTCGCCTCGGATCCTTACTACACTAAC	ADCH?@B?BGIAIH5D5G@DAB@DI5FBHH@GAGEI	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1148:2044	147	chr1	731	60	36M	=	598	-169	TCAAGTTAGTGTAGTAAGGATCCGAGGCGAGTCTGA	@5HGG?CGHHI5ID55BIIF5BD?G5H5A@EFGBHF	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1222:2066	147	chr1	751	60	36M	=	609	-178	ACTCTTTGACCACTAGGCGTTCAAGTTAGTGTAGTA	AD5BH5C@AI@?ABF@FEEA@CCAIDDHGABEFAAE	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1074:2022	147	chr1	774	60	36M	=	610	-200	TATAGATACCGACGATTACCAGTACTCTTTGACCAC	C@HAGGBF?IIEBHIHI@D?5DFG?EDG@DBGG?DE	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1037:2011	147	chr1	794	60	36M	=	646	-184	TGTTTCCCCTCCCCTGCTTATATAGATACCGACGAT	DGBBAI?IFBEGIGAHDAFICEFBGEEB5ACAAFFD	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1370:2110	147	chr1	826	60	36M	=	693	-169	CTTAGCATTAGGAGTCACCGGCTGAGAACAAATGTT	IICIGB?HEDEAFC5C?HBHGIHIBHA?@I55D@BC	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1444:2132	147	chr1	829	60	36M	=	693	-172	TGTCTTAGCATTAGGAGTCACCGGCTGAGAACAAAT	EHH@BCFH5AAHIEI?@HGCFGDCFAGEDB?A@AFI	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1259:2077	147	chr1	832	60	36M	=	702	-166	AAATGTCTTAGCATTAGGAGTCACCGGCTGAGAACA	H??FCFBA?IB?AB@IDBCDEHH@B@5CHHICD5HI	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1111:2033	147	chr1	892	60	36M	=	720	-208	GTGCCTGCTTCAGCTGGTTGCTCAGATTTATGGCAT	55BGEBFBFGDFFIHDDEC5GHCE5H?FA?DCFDHF	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1004:2004	147	chr1	892	60	36M	=	720	-208	GTGCCTGCTTCAGCTGGTTGCTCAGATTTATGGCAT	55BGEBFBFGDFFIHDDEC5GHCE5H?FA?DCFDHF	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1296:2088	147	chr1	900	60	36M	=	727	-209	GCACTGTCGTGCCTGCTTCAGCTGGTTGCTCAGATT	FBDHHGC5?HB5GGI5AH?5CACI?DIIDIGABF?C	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1666:2198	99	chr2	194	60	36M	=	354	196	CCCATTGGATCTCGTTTATCGATTAAGCCCGATCTA	DCFHGF@HGGCAIGBIE5HBHFIEFBCFFB@B5FE?	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1592:2176	99	chr2	196	60	36M	=	375	215	CATTGGATCTCGTTTATCGATTAAGCCCGATCTAGG	HACFIDC5ICHCEFCAAC@@A5HEDEDGB@G?@BIB	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1518:2154	99	chr2	198	60	36M	=	355	193	TTGGATCTCGTTTATCGATTAAGCCCGATCTAGGTT	IH5B?ACFIFCIE5GC?CFCG?5HABCAIE5GE?5E	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1629:2187	99	chr2	215	60	36M	=	386	207	ATTAAGCCCGATCTAGGTTCCTAGAGGTTAAATTGG	AGBFH?AAI5@G5??GCI@5ABFH5HC5ACCI5HED	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1703:2209	99	chr2	246	60	36M	=	387	177	ATTGGACGTCTTCCCACTCCGTTGCTGCGTGTCTAG	C5?HI5?HDHC5AGDDDHHA5B5BIDAACACDHBBE	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1666:2198	147	chr2	354	60	36M	=	194	-196	CCGCTAGCCCAACCCTCACACGACCTCGAGTGAATC	HC@HAD?5G@?5FFAF@BGCI@@I@F5CIAEEAHCD	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1518:2154	147	chr2	355	60	36M	=	198	-193	GCCGCTAGCCCAACCCTCACACGACCTCGAGTGAAT	??B@@FBHHD@GBFIBE55C@EFE55?@GHHGDE@I	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1592:2176	147	chr2	375	60	36M	=	196	-215	GATGTGATAGTTTCATAATTGCCGCTAGCCCAACCC	BIGFHC?AG?G@BGCECIDI?EC@BBF5@HBAI5A5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1629:2187	147	chr2	386	60	36M	=	215	-207	CCCGCTTATGTGATGTGATAGTTTCATAATTGCCGC	GHC@5D5?HGCEGDBE55CGHC5DGIIC@?5@A@F?	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1703:2209	147	chr2	387	60	36M	=	246	-177	GCCCGCTTATGTGATGTGATAGTTTCATAATTGCCG	AG@EB@BB?C5EA@CHGGEAG5AIC5E@D@BH5?@5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1555:2165	99	chr2	443	60	36M	=	614	207	TCCATAAAACACTAGCTCAGCAGTTGAAAAAATGGC	DAGF?CCFAB@GG5A@CIECGEDCC5CGECA5AEG5	RG:Z:rg1	MC:Z:36M	NM:i:0
M1:7:FC1:1:1101:1555:2165	147	chr2	614	60	36M	=	443	-207	GCTGGTGATGGAGAACTACCCATACGTCGTTACTTT	H@IH@BDB?FBCGGFG@I5F?ADHGH?CBA@H?BCI	RG:Z:rg1	MC:Z:36M	NM:i:0
```

**Example file:** `ref.fa`

```text
>chr1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC
TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC
GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA
CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT
ACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAA
GACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATTTGCCGCCTGACAAGTCAATG
CGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATT
AACTGATAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAA
CGGCGAGCTTTACATTTGCTGTGAGAGGTACAGGGATTAGTGAGAAGCCGTGCGTATCAA
TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGC
TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTA
GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCA
AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAG
CCGGTGACTCCTAATGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAA
ATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTA
GCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGATACCTCTCCATCTGACCCAAG
ATTGTGCTTGTTCAATTCTTCTTAACGTGATAACAGAATCAAACCTGCCAGGCGGTCGTC
GCGGACCTCGGTCGAAGTAGTGGTGCGGATCCAGGGGAACCGTTGACTCAAAAGGAGCTG
CCGTCCACCTAACGTGAAGTTCCAAAATCCCAAACCTCTCGAGATATTTATCCAGCAAGG
>chr2
AGTGGCAACGCCCGCTGCTTTAATCGCTACCAAAACGCAAACAAAAGCATACCCAAAAGT
ACACGGGTGAGGGAGGTGATATAGTACAGCTACGAAGTATCTGGCGCCTCAATAGGATTA
TAGCGGTCTCTCAGGCTGCTTGCCGTCCGGCCCGGCCGCGACACTCCGGTGCAAGCTTAA
TTCGTACGTACTTCCCATTGGATCTCGTTTATCGATTAAGCCCGATCTAGGTTCCTAGAG
GTTAAATTGGACGTCTTCCCACTCCGTTGCTGCGTGTCTAGGCGGTTTAGCGTAAGCGAA
CAGGACCCTGCCTCAGCTCATAAGTCCTTATTCTCTCACGTTGTGTTACGAAAGATTCAC
TCGAGGTCGTGTGAGGGTTGGGCTAGCGGCAATTATGAAACTATCACATCACATAAGCGG
GCTAGATATAATTTAATCTTAATCCATAAAACACTAGCTCAGCAGTTGAAAAAATGGCTA
GGTTCCAGCTTTTGGGGAGACGTCTTTCTGAGGGTCAGCCGTGATTCCGATTCGATTAGA
CTGGTCCCCACGGGTCCATGAGTACGAGGAAACTCGGTATCGAGCCTAAAAGTTATAAGG
CATCTCGCCCAGGAAAGTAACGACGTATGGGTAGTTCTCCATCACCAGCTATAATGGCTA
GCGCACTCTCGTTCCAGGGCGTAGTTACACTGAGCGTGCCATGTCAGCATGCTAGCGTAT
CGCCCCCCAATGCCCCGCAATAGGGTAATTCGCCGACGAGTAAGCGTAGATTACACACCC
AGGAAACGATCTAGACAGAT
```

**Example file:** `ref.fa.fai`

```text
chr1	1200	6	60	61
chr2	800	1232	60	61
```

**Example file:** `ref.dict`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1200
@SQ	SN:chr2	LN:800
```

**Example file:** `hap.map`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:1200
@SQ	SN:chr2	LN:800
#CHROMOSOME	POSITION	NAME	MAJOR_ALLELE	MINOR_ALLELE	MAF	ANCHOR_SNP	PANELS
chr1	518	snp518	C	G	0.3	snp518	
chr1	525	snp525	A	C	0.3	snp525	
chr1	532	snp532	T	A	0.3	snp532	
chr1	539	snp539	A	C	0.3	snp539	
chr1	546	snp546	A	C	0.3	snp546	
chr1	609	snp609	T	A	0.3	snp609	
chr1	616	snp616	T	A	0.3	snp616	
chr1	623	snp623	C	G	0.3	snp623	
```

**Example file:** `genotypes.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1200>
##contig=<ID=chr2,length=800>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	518	.	C	G	50	.	DP=10;QD=5;FS=1.0	GT:DP:GQ	0/0:8:30
chr1	525	.	A	C	50	.	DP=10;QD=6;FS=1.0	GT:DP:GQ	0/0:9:31
chr1	532	.	T	A	50	.	DP=10;QD=7;FS=1.0	GT:DP:GQ	0/0:10:32
chr1	539	.	A	C	50	.	DP=10;QD=8;FS=1.0	GT:DP:GQ	0/0:11:33
chr1	546	.	A	C	50	.	DP=10;QD=9;FS=1.0	GT:DP:GQ	0/0:12:34
chr1	609	.	T	A	50	.	DP=10;QD=10;FS=1.0	GT:DP:GQ	0/0:13:35
chr1	616	.	T	A	50	.	DP=10;QD=11;FS=1.0	GT:DP:GQ	0/0:14:36
chr1	623	.	C	G	50	.	DP=10;QD=12;FS=1.0	GT:DP:GQ	0/0:15:37
```

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1200>
##contig=<ID=chr2,length=800>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	518	.	C	G	50	.	DP=10;QD=5;FS=1.0	GT:DP:GQ	0/1:8:30
chr1	525	.	A	C	50	.	DP=10;QD=6;FS=1.0	GT:DP:GQ	0/1:9:31
chr1	532	.	T	A	50	.	DP=10;QD=7;FS=1.0	GT:DP:GQ	0/1:10:32
chr1	539	.	A	C	50	.	DP=10;QD=8;FS=1.0	GT:DP:GQ	0/1:11:33
chr1	546	.	A	C	50	.	DP=10;QD=9;FS=1.0	GT:DP:GQ	0/1:12:34
chr1	609	.	T	A	50	.	DP=10;QD=10;FS=1.0	GT:DP:GQ	0/1:13:35
chr1	616	.	T	A	50	.	DP=10;QD=11;FS=1.0	GT:DP:GQ	0/1:14:36
chr1	623	.	C	G	50	.	DP=10;QD=12;FS=1.0	GT:DP:GQ	0/1:15:37
```

The next step derives the remaining inputs (indexed VCFs, a BAM, a truncated BAM, tagged reads, annotation and estimator files, an array VCF and a tiny Illumina run folder) from the files above.

- **Example code:**

```python
import struct
from pathlib import Path

import pysam

for name in ("sample", "genotypes"):
    pysam.tabix_index(f"{name}.vcf", preset="vcf", force=True, keep_original=True)
Path("a.interval_list").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1200\n@SQ\tSN:chr2\tLN:800\n"
    "chr1\t100\t200\t+\ta\nchr1\t300\t400\t+\tb\nchr2\t10\t50\t+\tc\n"
)
Path("b.interval_list").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1200\n@SQ\tSN:chr2\tLN:800\nchr1\t150\t180\t+\tx\n"
)
Path("messy.fa").write_text(">s1 desc\n" + "ACGT" * 10 + "\nAC\n>s2\nACGTAC\nGT\n")
sequence = "ACGT" * 25 + "N" * 20 + "TTGCA" * 20
Path("gap.fa").write_text(f">chrN\n{sequence}\n")
Path("gap.fa.fai").write_text(f"chrN\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n")
Path("gap.dict").write_text(f"@HD\tVN:1.6\n@SQ\tSN:chrN\tLN:{len(sequence)}\n")
Path("x.gff3").write_text(
    "##gff-version 3\nchr2\tsrc\tgene\t10\t50\t.\t+\t.\tID=g2\n"
    "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g1b\nchr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1a\n"
)
Path("a_to_b.chain").write_text(
    "chain 1000 chr1 1200 + 0 1200 chrA 1200 + 0 1200 1\n1200\n\n"
    "chain 1000 chr2 800 + 0 800 chrB 800 + 0 800 2\n800\n\n"
)
Path("target.fa").write_text(Path("ref.fa").read_text().replace(">chr1", ">chrA").replace(">chr2", ">chrB"))
Path("target.fa.fai").write_text(Path("ref.fa.fai").read_text().replace("chr1", "chrA").replace("chr2", "chrB"))
Path("target.dict").write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:1200\n@SQ\tSN:chrB\tLN:800\n")
Path("swapped.dict").write_text("@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:800\n@SQ\tSN:chr1\tLN:1200\n")
Path("mixed.vcf").write_text(
    Path("sample.vcf").read_text()
    + "chr1\t1000\t.\tA\tAT\t50\t.\tDP=10;QD=5;FS=1.0\tGT:DP:GQ\t0/1:8:30\n"
)
clusters, base = 4, Path("run/Data/Intensities")
(base / "BaseCalls" / "L001" / "C1.1").mkdir(parents=True, exist_ok=True)
(base / "BaseCalls" / "L001" / "C2.1").mkdir(parents=True, exist_ok=True)
(base / "L001").mkdir(exist_ok=True)
Path("run/InterOp").mkdir(exist_ok=True)
Path("run/InterOp/TileMetricsOut.bin").write_bytes(bytes([2, 10]) + struct.pack("<HHHf", 1, 1101, 100, 1000.0))
(base / "BaseCalls" / "L001" / "s_1_1101.filter").write_bytes(struct.pack("<III", 0, 3, clusters) + bytes([1, 1, 0, 1]))
for cycle, calls in ((1, (0, 1, 2, 3)), (2, (3, 2, 1, 0))):
    path = base / "BaseCalls" / "L001" / f"C{cycle}.1" / "s_1_1101.bcl"
    path.write_bytes(struct.pack("<I", clusters) + bytes((30 << 2) | call for call in calls))
(base / "L001" / "s_1_1101.locs").write_bytes(
    struct.pack("<IfI", 1, 1.0, clusters) + b"".join(struct.pack("<ff", 100.0 + i * 10, 200.0 + i * 10) for i in range(clusters))
)
print(sorted(path.name for path in Path(".").iterdir() if path.suffix in {".gz", ".chain", ".interval_list", ".gff3"} or path.name == "run"))
```

- **Example output:**

```text
['a.interval_list', 'a_to_b.chain', 'b.interval_list', 'genotypes.vcf.gz', 'run', 'sample.vcf.gz', 'x.gff3']
```

## 1) `dnakit.alignment.picard_downsample_sam()`

- **Purpose:** Keep a deterministic random subset of the reads of a SAM/BAM/CRAM file, by query name (`DownsampleSam`).
- **Method:** Reads of one template (pairs, secondary and supplementary alignments) are kept or dropped together. `probability` is the chance of keeping a template (Picard's default is 1.0, which keeps everything); `strategy` is `ConstantMemory` (default), `HighAccuracy` or `Chained`; `random_seed` makes the choice reproducible. `metrics_file` receives a `QualityYieldMetrics` file of the downsampled output.
- **API:** `dnakit.alignment.picard_downsample_sam(input_path, output_path, *, probability=None, strategy=None, random_seed=None, accuracy=None, metrics_file=None, reference_sequence=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and `probability`; optional `strategy`, `random_seed`, `metrics_file`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` of the downsampled file when `metrics_file` is given.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_downsample_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


kept = picard_downsample_sam("example.sam", str(out / "half.sam"), probability=0.5, random_seed=11, metrics_file=str(out / "half.txt"), verbose=False)
names = [read.query_name for read in reads(out / "half.sam")]
print(kept["status"], len(names), all(names.count(n) == 2 for n in names), kept["metrics"][0]["TOTAL_READS"] == len(names))
again = picard_downsample_sam("example.sam", str(out / "half2.sam"), probability=0.5, random_seed=11, verbose=False)
print([read.query_name for read in reads(out / "half2.sam")] == names)
```

- **Example output:**

```text
completed 22 True True
True
```

## 2) `dnakit.alignment.picard_position_based_downsample_sam()`

- **Purpose:** Downsample a SAM/BAM/CRAM file by flowcell position, keeping approximately `fraction` of the templates (`PositionBasedDownsampleSam`).
- **Method:** The choice is deterministic from each read's tile and x/y cluster coordinates (parsed from the read name with `read_name_regex`), not random, so optical-duplicate structure is kept. Reads of one template stay together. `fraction` must be between 0 and 1.
- **API:** `dnakit.alignment.picard_position_based_downsample_sam(input_path, output_path, fraction, *, read_name_regex=None, remove_duplicate_information=None, allow_multiple_downsampling_despite_warnings=None, stop_after=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and `fraction` between 0 and 1.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the downsampled file in `outputs`.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_position_based_downsample_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


for name in ("a", "b"):
    picard_position_based_downsample_sam("example.sam", str(out / f"{name}.sam"), 0.5, verbose=False)
first = [read.query_name for read in reads(out / "a.sam")]
print(len(first) < 44, first == [read.query_name for read in reads(out / "b.sam")])
```

- **Example output:**

```text
True True
```

## 3) `dnakit.alignment.picard_filter_sam_reads()`

- **Purpose:** Subset the reads of a SAM/BAM/CRAM file by one filter (`FilterSamReads`).
- **Method:** `filter` is one of `includeAligned` / `excludeAligned` (the input must be queryname sorted and both mates must be aligned to be included), `includeReadList` / `excludeReadList` (`read_list_file` holds one read name per line), `includePairedIntervals` (`interval_list`) and `includeTagValues` / `excludeTagValues` (`tag` and string `tag_values`). The JavaScript filter is not offered. `sort_order` defaults to the input's.
- **API:** `dnakit.alignment.picard_filter_sam_reads(input_path, output_path, filter, *, read_list_file=None, interval_list=None, tag=None, tag_values=None, sort_order=None, write_reads_files=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and `filter` (`includeAligned`, `excludeAligned`, `includeReadList`, `excludeReadList`, `includePairedIntervals`, `includeTagValues`, `excludeTagValues`) with its companion argument.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the filtered file in `outputs`.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_filter_sam_reads


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


names = [read.query_name for read in reads(Path("example.sam"))]
(out / "names.txt").write_text(names[0] + "\n")
one = picard_filter_sam_reads("example.sam", str(out / "one.sam"), "includeReadList", read_list_file=str(out / "names.txt"), verbose=False)
print(one["status"], len(reads(out / "one.sam")))
tag = picard_filter_sam_reads("example.sam", str(out / "tag.sam"), "includeTagValues", tag="RG", tag_values=["rg2"], verbose=False)
print(tag["status"], len(reads(out / "tag.sam")))
```

- **Example output:**

```text
completed 2
completed 0
```

## 4) `dnakit.alignment.picard_split_sam_by_number_of_reads()`

- **Purpose:** Split a query-grouped SAM/BAM/CRAM file into several files with the same sort order (`SplitSamByNumberOfReads`).
- **Method:** Give exactly one of `split_to_n_files` or `split_to_n_reads` (about that many reads per file). `output_dir` is an existing directory; the files `<out_prefix>_N.<ext>` (`shard` by default, extension of the input) are written into it and must not exist yet. The new files are listed in `files`. `total_reads_in_input` saves a first pass over the input.
- **API:** `dnakit.alignment.picard_split_sam_by_number_of_reads(input_path, output_dir, *, split_to_n_files=None, split_to_n_reads=None, out_prefix=None, total_reads_in_input=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, an existing output directory and exactly one of `split_to_n_files` and `split_to_n_reads`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `files`: the new shard files.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_revert_sam, picard_split_sam_by_number_of_reads


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


picard_revert_sam("example.sam", str(out / "unmapped.sam"), verbose=False)
shards = out / "shards"
shards.mkdir()
result = picard_split_sam_by_number_of_reads(str(out / "unmapped.sam"), str(shards), split_to_n_files=3, verbose=False)
print([Path(p).name for p in result["files"]], sum(len(reads(Path(p))) for p in result["files"]))
```

- **Example output:**

```text
['shard_0001.sam', 'shard_0002.sam', 'shard_0003.sam'] 44
```

## 5) `dnakit.alignment.picard_reorder_sam()`

- **Purpose:** Reorder the reads of a SAM/BAM/CRAM file to match the contig order of a new reference (`ReorderSam`).
- **Method:** Contigs are matched by exact name; reads on contigs absent from the new reference are unmapped. `sequence_dictionary` may be a SAM/BAM/CRAM, VCF/BCF, interval list, FASTA or `.dict` file. Both `allow_*` flags are dangerous relaxations and default to Picard's false.
- **API:** `dnakit.alignment.picard_reorder_sam(input_path, output_path, sequence_dictionary, *, allow_contig_length_discordance=None, allow_incomplete_dict_concordance=None, reference_sequence=None, validation_stringency=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and a sequence dictionary (SAM/BAM/CRAM, VCF, interval list, FASTA or `.dict`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the reordered file in `outputs`.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_reorder_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


result = picard_reorder_sam("example.sam", str(out / "re.sam"), "swapped.dict", verbose=False)
contigs = [read.reference_name for read in reads(out / "re.sam")]
print(result["status"], contigs[0], contigs[-1])
```

- **Example output:**

```text
completed chr2 chr1
```

## 6) `dnakit.alignment.picard_revert_sam()`

- **Purpose:** Revert an aligned SAM/BAM/CRAM file to an unaligned file, restoring original qualities (`RevertSam`).
- **Method:** With Picard's defaults the alignment and duplicate flags are removed, original qualities (`OQ`) are restored, hard-clipped bases are restored and the output is queryname sorted. `attribute_to_clear` replaces nothing: its tags are added to Picard's default list (`NM UQ PG MD MQ SA MC AS`). With `output_by_readgroup=True`, `output_path` is a new directory that receives one file per read group (`output_by_readgroup_file_format` is `sam`, `bam`, `cram` or `dynamic`, by default the format of the input); its files are listed in `files`. `sanitize` is destructive and discards inconsistent reads.
- **API:** `dnakit.alignment.picard_revert_sam(input_path, output_path, *, output_by_readgroup=None, output_by_readgroup_file_format=None, remove_alignment_information=None, remove_duplicate_information=None, restore_original_qualities=None, restore_hardclips=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file and a new output path (or an existing directory with `output_by_readgroup=True`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `files` when writing one file per read group.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_revert_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


result = picard_revert_sam("example.sam", str(out / "u.sam"), verbose=False)
reverted = reads(out / "u.sam")
print(result["status"], len(reverted), all(read.is_unmapped for read in reverted), any(read.has_tag("MC") for read in reverted))
folder = out / "per_rg"
folder.mkdir()
print([Path(p).name for p in picard_revert_sam("example.sam", str(folder), output_by_readgroup=True, verbose=False)["files"]])
```

- **Example output:**

```text
completed 44 True False
['rg1.sam']
```

## 7) `dnakit.alignment.picard_mark_illumina_adapters()`

- **Purpose:** Find Illumina adapter sequences in queryname-sorted reads and write the `XT` trimming tag (`MarkIlluminaAdapters`).
- **Method:** `metrics_path` receives a histogram of bases clipped per read (parsed in `histogram`). Without `output_path` only the metrics are produced. `adapters` are Picard adapter pairs (`PAIRED_END`, `INDEXED`, `SINGLE_END`, `NEXTERA_V1`, `NEXTERA_V2`, `DUAL_INDEXED`, `FLUIDIGM`, `TRUSEQ_SMALLRNA`, `ALTERNATIVE_SINGLE_END`); by default `INDEXED`, `DUAL_INDEXED` and `PAIRED_END`. Existing `XT` tags are cleared.
- **API:** `dnakit.alignment.picard_mark_illumina_adapters(input_path, metrics_path, *, output_path=None, adapters=None, five_prime_adapter=None, three_prime_adapter=None, max_error_rate_pe=None, max_error_rate_se=None, ..., verbose=True)`.
- **Input:** A queryname-sorted SAM/BAM/CRAM file and a new metrics path; optional `output_path` and `adapters`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `outputs` and, when adapters are found, `histogram`.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_mark_illumina_adapters, picard_revert_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


picard_revert_sam("example.sam", str(out / "u.sam"), verbose=False)
result = picard_mark_illumina_adapters(str(out / "u.sam"), str(out / "adapters.txt"), output_path=str(out / "marked.sam"), adapters=["PAIRED_END"], verbose=False)
print(result["status"], sorted(result["outputs"]), len(reads(out / "marked.sam")))
```

- **Example output:**

```text
completed ['metrics', 'output'] 44
```

## 8) `dnakit.alignment.picard_merge_bam_alignment()`

- **Purpose:** Merge aligner output with the original unmapped BAM, restoring the tags the aligner dropped (`MergeBamAlignment`).
- **Method:** `unmapped_bam` must be queryname sorted and contain only unmapped reads. Give the alignments either as `aligned_bams` or as separate `read1_aligned_bams` and `read2_aligned_bams` (not both). The merged file is coordinate sorted unless `sort_order` says otherwise. `attributes_to_reverse` and `attributes_to_reverse_complement` replace Picard's defaults (`OQ U2` and `E2 SQ`) only when given. `primary_alignment_strategy` is `BestMapq` (Picard's default), `EarliestFragment`, `BestEndMapq` or `MostDistant`.
- **API:** `dnakit.alignment.picard_merge_bam_alignment(unmapped_bam, reference_sequence, output_path, *, aligned_bams=None, read1_aligned_bams=None, read2_aligned_bams=None, sort_order=None, add_mate_cigar=None, clip_adapters=None, ..., verbose=True)`.
- **Input:** The unmapped BAM, the reference FASTA, a new output path and `aligned_bams` (or `read1_aligned_bams` and `read2_aligned_bams`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the merged file in `outputs`.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_merge_bam_alignment, picard_revert_sam


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


picard_revert_sam("example.sam", str(out / "u.sam"), verbose=False)
result = picard_merge_bam_alignment(str(out / "u.sam"), "ref.fa", str(out / "merged.sam"), aligned_bams=["example.sam"], verbose=False)
merged = reads(out / "merged.sam")
print(result["status"], len(merged), sum(not read.is_unmapped for read in merged), {read.get_tag("RG") for read in merged})
```

- **Example output:**

```text
completed 44 44 {'rg1'}
```

## 9) `dnakit.alignment.picard_mark_duplicates_with_mate_cigar()`

- **Purpose:** Flag duplicate molecules in coordinate-sorted, mate-cigar-tagged reads (`SimpleMarkDuplicatesWithMateCigar` or `UmiAwareMarkDuplicatesWithMateCigar`).
- **Method:** Reads need `MC` mate-cigar tags. `umi_aware=False` runs `SimpleMarkDuplicatesWithMateCigar` (Picard marks it EXPERIMENTAL); `umi_aware=True` runs `UmiAwareMarkDuplicatesWithMateCigar`, which also groups reads by their UMI (tag `umi_tag_name`, `RX` by default; reads within `max_edit_distance_to_join` edits are joined) and needs `umi_metrics_file`. The duplication metrics are parsed in `metrics` (and `umi_metrics`). `tagging_policy` is `DontTag`, `OpticalOnly` or `All`; `duplicate_scoring_strategy` is `SUM_OF_BASE_QUALITIES` (default), `TOTAL_MAPPED_REFERENCE_LENGTH` or `RANDOM`.
- **API:** `dnakit.alignment.picard_mark_duplicates_with_mate_cigar(input_paths, output_path, metrics_file, *, umi_aware=False, umi_metrics_file=None, umi_tag_name=None, max_edit_distance_to_join=None, allow_missing_umis=None, remove_duplicates=None, ..., verbose=True)`.
- **Input:** Coordinate-sorted files with `MC` tags, a new output path and a new metrics path; with `umi_aware=True` also `umi_metrics_file`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` (and `umi_metrics`).
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_mark_duplicates_with_mate_cigar


def reads(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read for read in handle.fetch(until_eof=True)]


result = picard_mark_duplicates_with_mate_cigar(["example.sam"], str(out / "marked.sam"), str(out / "dup.txt"), verbose=False)
row = result["metrics"][0]
print(row["READ_PAIRS_EXAMINED"], row["READ_PAIR_DUPLICATES"], sum(read.is_duplicate for read in reads(out / "marked.sam")))
lines = [line if line.startswith("@") else line + "\tRX:Z:ACGT" for line in Path("example.sam").read_text().splitlines()]
(out / "umi.sam").write_text("\n".join(lines) + "\n")
umi = picard_mark_duplicates_with_mate_cigar([str(out / "umi.sam")], str(out / "umi_marked.sam"), str(out / "umi_dup.txt"), umi_aware=True, umi_metrics_file=str(out / "umi_metrics.txt"), verbose=False)
print(umi["status"], umi["umi_metrics"][0]["MEAN_UMI_LENGTH"])
```

- **Example output:**

```text
22 2 4
completed 4
```

## 10) `dnakit.alignment.picard_illumina_basecalls_to_fastq()`

- **Purpose:** Generate FASTQ files from the raw base calls of an Illumina run folder (`IlluminaBasecallsToFastq`).
- **Method:** `basecalls_dir` is `Data/Intensities/BaseCalls` of the run, `lanes` the lane numbers and `read_structure` the description of the clusters (for example `76T8B76T`). Give `output_prefix` for a run without demultiplexing (files `<prefix>.1.fastq`, `<prefix>.2.fastq`, `<prefix>.barcode_1.fastq` ...) or `multiplex_params`, a tab-separated file (`OUTPUT_PREFIX`, `BARCODE_1` ...) that demultiplexes by barcode. Read names are in Casava 1.8 style unless `read_name_format='ILLUMINA'`; that style needs `machine_name` and `flowcell_barcode`. The FASTQ files that exist afterwards are listed in `files`.
- **API:** `dnakit.alignment.picard_illumina_basecalls_to_fastq(basecalls_dir, lanes, read_structure, run_barcode, *, output_prefix=None, multiplex_params=None, machine_name=None, flowcell_barcode=None, read_name_format=None, compress_outputs=None, ..., verbose=True)`.
- **Input:** The BaseCalls directory, lane numbers, the read structure, the run barcode and `output_prefix` (or `multiplex_params`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `files`: the FASTQ files written.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_illumina_basecalls_to_fastq

(out / "fastq").mkdir()
result = picard_illumina_basecalls_to_fastq("run/Data/Intensities/BaseCalls", [1], "2T", "run1", output_prefix=str(out / "fastq" / "s"), machine_name="M1", flowcell_barcode="FC1", verbose=False)
records = Path(result["files"][0]).read_text().splitlines()
print([Path(p).name for p in result["files"]], records[1::4], records[0].rsplit(":", 2)[0])
```

- **Example output:**

```text
['s.1.fastq'] ['AT', 'CG', 'GC', 'TA'] @M1:run1:FC1:1:1101:2000:3000 :N
```

## 11) `dnakit.alignment.picard_illumina_basecalls_to_sam()`

- **Purpose:** Transform raw Illumina base calls into unmapped SAM/BAM/CRAM files, demultiplexed by barcode (`IlluminaBasecallsToSam`).
- **Method:** Single-sample run: give `output_path` (a new file whose extension picks the format), `sample_alias` and `library_name`. Multi-sample run: give `library_params`, a tab-separated file with the columns `OUTPUT`, `SAMPLE_ALIAS`, `LIBRARY_NAME` and `BARCODE_1` ...; the output files it names are listed in `files`. `read_structure` describes the clusters (for example `76T8B76T`); `sequencing_center` fills the `CN` tag of the read group; `molecular_index_tag` is `RX` by default; `run_start_date` is a date such as `2024-05-01`.
- **API:** `dnakit.alignment.picard_illumina_basecalls_to_sam(basecalls_dir, lanes, read_structure, run_barcode, sequencing_center, *, output_path=None, sample_alias=None, library_name=None, library_params=None, platform=None, read_group_id=None, ..., verbose=True)`.
- **Input:** The BaseCalls directory, lane numbers, the read structure, the run barcode, the sequencing centre and `output_path` with `sample_alias` (or `library_params`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `files`: the unmapped files written.
- **Example code:**

```python
import pysam
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.picard_sam_p08_cli_functions import picard_illumina_basecalls_to_sam

(out / "bam").mkdir()
result = picard_illumina_basecalls_to_sam("run/Data/Intensities/BaseCalls", [1], "2T", "run1", "XX", output_path=str(out / "bam" / "u.bam"), sample_alias="s1", library_name="l1", verbose=False)
with pysam.AlignmentFile(result["files"][0], check_sq=False) as handle:
    group = handle.header.to_dict()["RG"][0]
    sequences = sorted(read.query_sequence for read in handle.fetch(until_eof=True))
print(group["SM"], group["LB"], group["CN"], sequences)
```

- **Example output:**

```text
s1 l1 XX ['AT', 'CG', 'GC', 'TA']
```

## 12) `dnakit.variants.picard_filter_vcf()`

- **Purpose:** Apply hard filters to the genotypes and sites of a VCF/BCF (`FilterVcf`).
- **Method:** `min_dp` and `min_gq` filter individual genotypes (Picard default 0 = off), `min_ab` is the minimum allele balance of heterozygotes, `min_qd` the minimum QD and `max_fs` the maximum Fisher strand value that a site may have before it is filtered. The JavaScript filter is not offered.
- **API:** `dnakit.variants.picard_filter_vcf(input_path, output_path, *, min_dp=None, min_gq=None, min_ab=None, min_qd=None, max_fs=None, picard_path=None, ..., verbose=True)`.
- **Input:** A VCF/BCF, a new output path and the thresholds.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `records` and the `index` path.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_filter_vcf

result = picard_filter_vcf("sample.vcf", str(out / "f.vcf"), min_gq=34, verbose=False)
rows = [line.split("\t") for line in (out / "f.vcf").read_text().splitlines() if not line.startswith("#")]
print(result["records"], [("LowGQ" in row[9]) for row in rows])
site = picard_filter_vcf("sample.vcf", str(out / "g.vcf"), min_qd=8.5, verbose=False)
print([line.split("\t")[6] for line in (out / "g.vcf").read_text().splitlines() if not line.startswith("#")])
```

- **Example output:**

```text
8 [True, True, True, True, False, False, False, False]
['LowQD', 'LowQD', 'LowQD', 'LowQD', 'PASS', 'PASS', 'PASS', 'PASS']
```

## 13) `dnakit.variants.picard_liftover_vcf()`

- **Purpose:** Lift a VCF/BCF over to another reference build with a UCSC chain file (`LiftoverVcf`).
- **Method:** `reference_sequence` is the target FASTA (with a `.dict` and `.fai` next to it). The output is sorted and indexed in the target coordinates unless `disable_sort` is set; variants that cannot be lifted go to `reject_path`. `liftover_min_match` defaults to 1.0. With `recover_swapped_ref_alt` a bi-allelic SNP whose alternative allele is the new reference allele is swapped; `tags_to_reverse` (default `AF`) are then transformed to 1-x and `tags_to_drop` (default `MAX_AF`) removed.
- **API:** `dnakit.variants.picard_liftover_vcf(input_path, output_path, reject_path, chain, reference_sequence, *, liftover_min_match=None, recover_swapped_ref_alt=None, allow_missing_fields_in_header=None, tags_to_drop=None, tags_to_reverse=None, write_original_position=None, ..., verbose=True)`.
- **Input:** A VCF/BCF, new output and reject paths, a chain file and the target reference FASTA (with `.fai` and `.dict`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `lifted_records`, `rejected_records` and `index`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_liftover_vcf

result = picard_liftover_vcf("sample.vcf", str(out / "lifted.vcf"), str(out / "rej.vcf"), "a_to_b.chain", "target.fa", write_original_position=True, verbose=False)
rows = [line.split("\t") for line in (out / "lifted.vcf").read_text().splitlines() if not line.startswith("#")]
print(result["lifted_records"], result["rejected_records"], {row[0] for row in rows}, "OriginalContig=chr1" in rows[0][7])
```

- **Example output:**

```text
8 0 {'chrA'} True
```

## 14) `dnakit.variants.picard_make_sites_only_vcf()`

- **Purpose:** Remove all genotype columns from a VCF/VCF.gz/BCF, keeping the site-level information and annotations such as AN and AF (`MakeSitesOnlyVcf`).
- **Method:** `samples` keeps only the named samples' contribution to the site-level annotations (Picard's `SAMPLE`).
- **API:** `dnakit.variants.picard_make_sites_only_vcf(input_path, output_path, *, samples=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A VCF/VCF.gz/BCF and a new output path; optional `samples`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `records` and the `index` path.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_make_sites_only_vcf

result = picard_make_sites_only_vcf("sample.vcf.gz", str(out / "sites.vcf"), verbose=False)
rows = [line.split("\t") for line in (out / "sites.vcf").read_text().splitlines() if not line.startswith("#")]
print(result["records"], {len(row) for row in rows})
```

- **Example output:**

```text
8 {8}
```

## 15) `dnakit.variants.picard_make_vcf_sample_name_map()`

- **Purpose:** Write a tab-separated map from sample name to VCF/GVCF path, one line per input (`MakeVcfSampleNameMap`).
- **Method:** Each input must be a single-sample VCF/GVCF. The map is also returned as `entries` (a list of `sample` and `path` pairs). Picard 3.5.0 writes the path in the first column and the sample name in the second, although its documentation says the reverse; `entries` tells the two apart by the input paths.
- **API:** `dnakit.variants.picard_make_vcf_sample_name_map(input_paths, output_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** Single-sample VCF/GVCF files and a new output path.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `entries`: sample/path pairs.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_make_vcf_sample_name_map

result = picard_make_vcf_sample_name_map(["sample.vcf.gz", "genotypes.vcf.gz"], str(out / "map.tsv"), verbose=False)
print([(entry["sample"], Path(entry["path"]).name) for entry in result["entries"]])
```

- **Example output:**

```text
[('sample1', 'sample.vcf.gz'), ('sample1', 'genotypes.vcf.gz')]
```

## 16) `dnakit.variants.picard_rename_sample_in_vcf()`

- **Purpose:** Rename the sample of a single-sample VCF/BCF (`RenameSampleInVcf`).
- **Method:** `old_sample_name`, when given, asserts that this is the existing name.
- **API:** `dnakit.variants.picard_rename_sample_in_vcf(input_path, output_path, new_sample_name, *, old_sample_name=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A single-sample VCF/BCF, a new output path and `new_sample_name`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the renamed file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_rename_sample_in_vcf

result = picard_rename_sample_in_vcf("sample.vcf", str(out / "renamed.vcf"), "NEW", old_sample_name="sample1", verbose=False)
header = [line for line in (out / "renamed.vcf").read_text().splitlines() if line.startswith("#CHROM")][0]
print(result["status"], header.split("\t")[-1])
wrong = picard_rename_sample_in_vcf("sample.vcf", str(out / "bad.vcf"), "NEW", old_sample_name="other", verbose=False)
print(wrong["status"])
```

- **Example output:**

```text
completed NEW
failed
```

## 17) `dnakit.variants.picard_split_vcfs()`

- **Purpose:** Split a VCF/BCF into a SNP file and an indel file with identical headers (`SplitVcfs`).
- **Method:** Both outputs are indexed. With Picard's default `strict=True` a record that is neither a SNP nor an indel is an error; `strict=False` drops such records.
- **API:** `dnakit.variants.picard_split_vcfs(input_path, snp_output, indel_output, *, strict=None, sequence_dictionary=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A VCF/BCF and new SNP and indel output paths.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `snp_records` and `indel_records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_split_vcfs

result = picard_split_vcfs("mixed.vcf", str(out / "snp.vcf"), str(out / "indel.vcf"), verbose=False)
print(result["snp_records"], result["indel_records"])
```

- **Example output:**

```text
8 1
```

## 18) `dnakit.variants.picard_liftover_haplotype_map()`

- **Purpose:** Lift a fingerprinting haplotype database over to another reference with a UCSC chain file (`LiftOverHaplotypeMap`).
- **Method:** `sequence_dictionary` (any file a dictionary can be read from) is written into the output header.
- **API:** `dnakit.variants.picard_liftover_haplotype_map(input_path, output_path, chain, sequence_dictionary, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A haplotype map, a new output path, a chain file and a sequence dictionary.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `snps`: the number of lifted SNPs.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_liftover_haplotype_map

result = picard_liftover_haplotype_map("hap.map", str(out / "hap2.map"), "a_to_b.chain", "target.dict", verbose=False)
rows = [line.split("\t") for line in (out / "hap2.map").read_text().splitlines() if not line.startswith(("@", "#"))]
print(result["snps"], {row[0] for row in rows})
```

- **Example output:**

```text
8 {'chrA'}
```

## 19) `dnakit.variants.picard_interval_list_to_bed()`

- **Purpose:** Convert a Picard interval list into a BED file (`IntervalListToBed`).
- **Method:** `score` (0 to 1000, Picard's default 500) is written in the BED score column; `sort` (default true) sorts the intervals first. Interval lists are 1-based closed, BED is 0-based half-open; Picard converts.
- **API:** `dnakit.variants.picard_interval_list_to_bed(input_path, output_path, *, score=None, sort=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** An interval list and a new BED path; optional `score` and `sort`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `intervals`: the number of intervals.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_interval_list_to_bed

result = picard_interval_list_to_bed("a.interval_list", str(out / "a.bed"), score=7, verbose=False)
print(result["intervals"])
print((out / "a.bed").read_text().strip().replace("\t", " "))
```

- **Example output:**

```text
3
chr1 99 200 a 7 +
chr1 299 400 b 7 +
chr2 9 50 c 7 +
```

## 20) `dnakit.variants.picard_interval_list_tools()`

- **Purpose:** Combine and manipulate interval lists (`IntervalListTools`).
- **Method:** `action` picks the set operation: `CONCAT` (all intervals of all inputs, no sorting or merging; Picard's default), `UNION` (merged, sorted), `INTERSECT` (loci present in every input), `SUBTRACT` (loci of the inputs not in `second_input_paths`), `SYMDIFF` (loci in exactly one side) and `OVERLAPS` (the whole intervals of the inputs that overlap any interval of `second_input_paths`). Inputs may be interval lists or VCFs. `padding` widens (or, when negative, shrinks) each interval before the action; `invert` returns the complement; `count_value` (`BASES` or `INTERVALS`) adds the number of bases or intervals of the output as `count`.
- **API:** `dnakit.variants.picard_interval_list_tools(input_paths, output_path, *, action='CONCAT', second_input_paths=None, padding=None, invert=None, sort=None, unique=None, ..., verbose=True)`.
- **Input:** Interval lists (or VCFs), a new output path, `action` and, for `SUBTRACT`/`SYMDIFF`/`OVERLAPS`, `second_input_paths`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `intervals` and, with `count_value`, `count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_interval_list_tools

both = picard_interval_list_tools(["a.interval_list", "b.interval_list"], str(out / "c.il"), verbose=False)
print("CONCAT", both["intervals"])
hit = picard_interval_list_tools(["a.interval_list"], str(out / "o.il"), action="OVERLAPS", second_input_paths=["b.interval_list"], count_value="INTERVALS", verbose=False)
print("OVERLAPS", hit["intervals"], hit["count"])
rest = picard_interval_list_tools(["a.interval_list"], str(out / "s.il"), action="SUBTRACT", second_input_paths=["b.interval_list"], count_value="BASES", verbose=False)
print("SUBTRACT bases", rest["count"])
```

- **Example output:**

```text
CONCAT 4
OVERLAPS 1 1
SUBTRACT bases 212
```

## 21) `dnakit.variants.picard_scatter_intervals_by_ns()`

- **Purpose:** Split a reference into alternating `N` and `ACGT` intervals, for scattering a variant-calling pipeline (`ScatterIntervalsByNs`).
- **Method:** The FASTA needs an index (`.fai`) and a dictionary (`.dict`). `output_type` is `BOTH` (Picard's default), `N` or `ACGT`; `max_to_merge` (default 1) is the number of contiguous N bases tolerated inside an ACGT interval.
- **API:** `dnakit.variants.picard_scatter_intervals_by_ns(reference_sequence, output_path, *, output_type=None, max_to_merge=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A reference FASTA (indexed, with dictionary) and a new interval-list path; optional `output_type` and `max_to_merge`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `intervals`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_scatter_intervals_by_ns

result = picard_scatter_intervals_by_ns("gap.fa", str(out / "n.il"), verbose=False)
print([tuple(line.split("\t")[1:3]) + (line.split("\t")[4],) for line in (out / "n.il").read_text().splitlines() if not line.startswith("@")])
```

- **Example output:**

```text
[('1', '100', 'ACGTmer'), ('101', '120', 'Nmer'), ('121', '220', 'ACGTmer')]
```

## 22) `dnakit.variants.picard_normalize_fasta()`

- **Purpose:** Rewrite a FASTA so that every sequence line has the same length, except the last of each record (`NormalizeFasta`).
- **Method:** `line_length` defaults to Picard's 100 bases; `truncate_sequence_names_at_whitespace` shortens record names at the first whitespace. Picard writes the sequence as it reads it (it does not change case).
- **API:** `dnakit.variants.picard_normalize_fasta(input_path, output_path, *, line_length=None, truncate_sequence_names_at_whitespace=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A FASTA, a new output path and `line_length`; optional `truncate_sequence_names_at_whitespace`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_normalize_fasta

result = picard_normalize_fasta("messy.fa", str(out / "n.fa"), line_length=16, truncate_sequence_names_at_whitespace=True, verbose=False)
print(result["records"])
print((out / "n.fa").read_text().strip())
```

- **Example output:**

```text
2
>s1
ACGTACGTACGTACGT
ACGTACGTACGTACGT
ACGTACGTAC
>s2
ACGTACGT
```

## 23) `dnakit.variants.picard_sort_gff()`

- **Purpose:** Sort a GFF3 file by coordinates so that it can be indexed, adding flush directives where possible (`SortGff`).
- **Method:** Contigs are ordered by `sequence_dictionary` when given, otherwise lexicographically. `n_records_in_memory` (Picard's `nRecordsInMemory`, default 50000) is the number of records held in RAM before spilling to disk.
- **API:** `dnakit.variants.picard_sort_gff(input_path, output_path, *, sequence_dictionary=None, n_records_in_memory=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A GFF3 and a new output path; optional `sequence_dictionary` and `n_records_in_memory`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `features`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_sort_gff

result = picard_sort_gff("x.gff3", str(out / "sorted.gff3"), verbose=False)
print([line.split("\t")[8] for line in (out / "sorted.gff3").read_text().splitlines() if line and not line.startswith("#")])
ordered = picard_sort_gff("x.gff3", str(out / "ordered.gff3"), sequence_dictionary="swapped.dict", verbose=False)
print([line.split("\t")[8] for line in (out / "ordered.gff3").read_text().splitlines() if line and not line.startswith("#")])
```

- **Example output:**

```text
['ID=g1a', 'ID=g1b', 'ID=g2']
['ID=g2', 'ID=g1a', 'ID=g1b']
```
