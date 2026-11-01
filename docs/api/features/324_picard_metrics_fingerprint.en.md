# Picard metrics, quality-control and fingerprint tools

This group of functions in `dnakit.sequencing_qc` is a set of command-line facades over [Picard](https://github.com/broadinstitute/picard) (Broad Institute, MIT licence) metrics, quality-control and fingerprinting tools. Each function runs one fixed Picard tool with an explicit `KEY=value` argument list (no shell), refuses to overwrite an existing output, and returns a dictionary with `status`, `backend`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail` and, where Picard writes a metrics file, the parsed `metrics`. `status` is `completed` only when the exit code is accepted and every expected output exists. Arguments left as `None` are not passed, so Picard's own defaults apply. Tested with Picard 3.5.0.

Picard is an external program. Install it (`conda install -c bioconda picard`, or `picard.jar` plus Java 17) and put `picard` on PATH, pass `picard_path` (a launcher or a `picard.jar`) or set `PICARD_JAR`; if it is missing the function raises `BackendUnavailableError` with an install hint. Charts (`chart_output`) need `Rscript`. The examples import from the module directly and call every function with `verbose=False`.

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

**Example file:** `dbsnp.vcf`

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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	dbsnp
chr1	518	rs518	C	G	50	.	DP=10;QD=5;FS=1.0	GT:DP:GQ	0/0:8:30
chr1	525	rs525	A	C	50	.	DP=10;QD=6;FS=1.0	GT:DP:GQ	0/0:9:31
chr1	532	rs532	T	A	50	.	DP=10;QD=7;FS=1.0	GT:DP:GQ	0/0:10:32
chr1	539	rs539	A	C	50	.	DP=10;QD=8;FS=1.0	GT:DP:GQ	0/0:11:33
chr1	546	rs546	A	C	50	.	DP=10;QD=9;FS=1.0	GT:DP:GQ	0/0:12:34
chr1	609	rs609	T	A	50	.	DP=10;QD=10;FS=1.0	GT:DP:GQ	0/0:13:35
chr1	616	rs616	T	A	50	.	DP=10;QD=11;FS=1.0	GT:DP:GQ	0/0:14:36
chr1	623	rs623	C	G	50	.	DP=10;QD=12;FS=1.0	GT:DP:GQ	0/0:15:37
```

The next step derives the remaining inputs (indexed VCFs, a BAM, a truncated BAM, tagged reads, annotation and estimator files, an array VCF and a tiny Illumina run folder) from the files above.

- **Example code:**

```python
import struct
from pathlib import Path

import pysam

for name in ("genotypes", "sample", "dbsnp"):
    pysam.tabix_index(f"{name}.vcf", preset="vcf", force=True, keep_original=True)
sam = pysam.AlignmentFile("example.sam")
bam = pysam.AlignmentFile("example.bam", "wb", template=sam)
for read in sam:
    bam.write(read)
bam.close()
Path("example.bam").read_bytes()
Path("cut.bam").write_bytes(Path("example.bam").read_bytes()[:1000])
tagged = []
for line in Path("example.sam").read_text().splitlines():
    if not line.startswith("@"):
        quality = line.split("\t")[10]
        line += "".join(f"\tq{base}:Z:{quality}" for base in "acgt")
    tagged.append(line)
Path("snvq.sam").write_text("\n".join(tagged) + "\n")
Path("genes.refflat").write_text(
    "geneA\ttxA\tchr1\t+\t100\t900\t100\t900\t2\t100,500,\t300,900,\n"
    "geneB\ttxB\tchr2\t-\t50\t700\t50\t700\t1\t50,\t700,\n"
)
Path("bafregress.txt").write_text(
    "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\nS1\t0.012\t0.004\t3.0\t0.01\t0.99\t900\n"
)
Path("verifyidintensity.txt").write_text("ID\t%Mix\tLLK\tLLK0\n----\n1\t0.05\t-120.5\t-125.5\n")
controls = [
    ("DNP(High)", "Staining"), ("DNP(Bgnd)", "Staining"), ("Biotin(High)", "Staining"),
    ("Biotin(Bgnd)", "Staining"), ("Extension(A)", "Extension"), ("Extension(T)", "Extension"),
    ("Extension(C)", "Extension"), ("Extension(G)", "Extension"), ("TargetRemoval", "TargetRemoval"),
    ("Hyb(High)", "Hybridization"), ("Hyb(Medium)", "Hybridization"), ("Hyb(Low)", "Hybridization"),
    ("String(PM)", "Stringency"), ("String(MM)", "Stringency"), ("NSB(Bgnd)Red", "Non-SpecificBinding"),
    ("NSB(Bgnd)Purple", "Non-SpecificBinding"), ("NSB(Bgnd)Blue", "Non-SpecificBinding"),
    ("NSB(Bgnd)Green", "Non-SpecificBinding"), ("NP(A)", "Non-Polymorphic"), ("NP(T)", "Non-Polymorphic"),
    ("NP(C)", "Non-Polymorphic"), ("NP(G)", "Non-Polymorphic"), ("Restore", "Restoration"),
]
header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1200>", "##contig=<ID=chr2,length=800>"]
header += [f"##{name}={name}|{category}|1|2" for name, category in controls]
header += [
    "##sampleAlias=sample1", "##arrayType=TestChip", "##autocallGender=F", "##autocallVersion=1.0",
    "##autocallDate=09/21/2016 20:40", "##imagingDate=8/15/2015 7:28:52 AM",
    "##extendedIlluminaManifestVersion=1.4", "##clusterFile=test.egt", "##p95Green=1000", "##p95Red=900",
    "##scannerName=iScan1",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
]
sites = [line.split("\t") for line in Path("sample.vcf").read_text().splitlines() if not line.startswith("#")]
rows = [f"{s[0]}\t{s[1]}\trs{s[1]}\t{s[3]}\t{s[4]}\t.\t.\t.\tGT\t{'0/0' if i % 2 else '0/1'}" for i, s in enumerate(sites)]
Path("arrays.vcf").write_text("\n".join(header + rows) + "\n")
pysam.tabix_index("arrays.vcf", preset="vcf", force=True, keep_original=True)
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
print(sorted(path.name for path in Path(".").iterdir() if path.suffix in {".bam", ".gz", ".tbi", ".refflat"} or path.name in {"run", "snvq.sam"}))
```

- **Example output:**

```text
['arrays.vcf.gz', 'arrays.vcf.gz.tbi', 'cut.bam', 'dbsnp.vcf.gz', 'dbsnp.vcf.gz.tbi', 'example.bam', 'genes.refflat', 'genotypes.vcf.gz', 'genotypes.vcf.gz.tbi', 'run', 'sample.vcf.gz', 'sample.vcf.gz.tbi', 'snvq.sam']
```

## 1) `dnakit.sequencing_qc.picard_accumulate_metrics()`

- **Purpose:** Merge sharded Picard metric files into one (`AccumulateQualityYieldMetrics` or `AccumulateVariantCallingMetrics`).
- **Method:** `metric_type="quality_yield"` merges `QualityYieldMetrics` files (`input_paths` are the files, `output_path` is the merged file). `metric_type="variant_calling"` merges Variant Calling metrics shards: `input_paths` and `output_path` are file prefixes without the `.variant_calling_detail_metrics` / `.variant_calling_summary_metrics` extensions, and both merged files are written.
- **API:** `dnakit.sequencing_qc.picard_accumulate_metrics(input_paths, output_path, *, metric_type='quality_yield', picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A list of metric files (`metric_type='quality_yield'`) or file prefixes without extension (`'variant_calling'`) and a new output path or prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `outputs` and the merged `metrics` (or `summary_metrics` and `detail_metrics`).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_accumulate_metrics, picard_variant_calling_metrics

for name in ("a", "b"):
    picard_variant_calling_metrics("sample.vcf.gz", "dbsnp.vcf.gz", str(out / name), verbose=False)
merged = picard_accumulate_metrics([str(out / "a"), str(out / "b")], str(out / "all"), metric_type="variant_calling", verbose=False)
print(merged["status"], merged["summary_metrics"][0]["TOTAL_SNPS"])
```

- **Example output:**

```text
completed 16
```

## 2) `dnakit.sequencing_qc.picard_format_contamination_metrics()`

- **Purpose:** Format a contamination estimator's output as a Picard metrics file.
- **Method:** `source="bafregress"` runs `CreateBafRegressMetricsFile` (writes `<prefix>.bafregress_metrics`); `source="verifyidintensity"` runs `CreateVerifyIDIntensityContaminationMetricsFile` (writes `<prefix>.verifyidintensity_metrics`). `input_path` is the estimator's captured standard output.
- **API:** `dnakit.sequencing_qc.picard_format_contamination_metrics(input_path, output_prefix, *, source='bafregress', picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** The captured standard output of the estimator and a new output prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the metrics file `<prefix>.bafregress_metrics` or `.verifyidintensity_metrics`, parsed in `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_format_contamination_metrics

baf = picard_format_contamination_metrics("bafregress.txt", str(out / "baf"), verbose=False)
print(baf["status"], baf["metrics"][0]["SAMPLE"], baf["metrics"][0]["ESTIMATE"], baf["metrics"][0]["NHOM"])
vid = picard_format_contamination_metrics("verifyidintensity.txt", str(out / "vid"), source="verifyidintensity", verbose=False)
print(vid["status"], vid["metrics"][0]["PCT_MIX"], vid["metrics"][0]["LLK"])
```

- **Example output:**

```text
completed S1 0.012 900
completed 0.05 -120.5
```

## 3) `dnakit.sequencing_qc.picard_convert_artifact_to_oxog()`

- **Purpose:** Convert `CollectSequencingArtifactMetrics` output into `CollectOxoGMetrics` format (`ConvertSequencingArtifactToOxoG`).
- **Method:** `input_base` is the prefix used for the artifact metrics (`<base>.pre_adapter_detail_metrics` and `<base>.bait_bias_detail_metrics` must exist); `output_base` is the prefix of the new `<base>.oxog_metrics` file.
- **API:** `dnakit.sequencing_qc.picard_convert_artifact_to_oxog(input_base, output_base, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** The prefix of an artifact-metrics run and a new output prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `<prefix>.oxog_metrics`, parsed in `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_convert_artifact_to_oxog, picard_sequencing_artifact_metrics

picard_sequencing_artifact_metrics("example.sam", str(out / "art"), "ref.fa", verbose=False)
converted = picard_convert_artifact_to_oxog(str(out / "art"), str(out / "oxo"), verbose=False)
print(converted["status"], len(converted["metrics"]) > 0, sorted(converted["metrics"][0])[:3])
```

- **Example output:**

```text
completed True ['ALT_NONOXO_BASES', 'ALT_OXO_BASES', 'CONTEXT']
```

## 4) `dnakit.sequencing_qc.picard_base_distribution_by_cycle()`

- **Purpose:** Nucleotide distribution per read cycle of a SAM/BAM/CRAM file (`CollectBaseDistributionByCycle`).
- **Method:** Picard's defaults apply (all reads, `assume_sorted` true). `chart_output` is an optional new `.pdf` and needs R.
- **API:** `dnakit.sequencing_qc.picard_base_distribution_by_cycle(input_path, output_path, *, aligned_reads_only=None, pf_reads_only=None, chart_output=None, assume_sorted=None, stop_after=None, validation_stringency=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file and a new output path; optional `aligned_reads_only`, `pf_reads_only`, `chart_output` (needs R).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `cycles` (one row per read and cycle) and `cycle_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_base_distribution_by_cycle

result = picard_base_distribution_by_cycle("example.sam", str(out / "bd.txt"), verbose=False)
print(result["status"], result["cycle_count"], round(sum(result["cycles"][0][k] for k in ("PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N"))))
```

- **Example output:**

```text
completed 72 100
```

## 5) `dnakit.sequencing_qc.picard_duplicate_metrics()`

- **Purpose:** Duplication metrics read from an already duplicate-marked SAM/BAM/CRAM file (`CollectDuplicateMetrics`).
- **Method:** Picard notes that the optical-duplicate count is always 0 here, so `ESTIMATED_LIBRARY_SIZE` is slightly off.
- **API:** `dnakit.sequencing_qc.picard_duplicate_metrics(input_path, metrics_file, *, assume_sorted=None, stop_after=None, validation_stringency=None, reference_sequence=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A duplicate-marked SAM/BAM/CRAM file and a new metrics path.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` with READ_PAIRS_EXAMINED, READ_PAIR_DUPLICATES, PERCENT_DUPLICATION.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_duplicate_metrics

row = picard_duplicate_metrics("example.sam", str(out / "dup.txt"), verbose=False)["metrics"][0]
print(row["LIBRARY"], row["READ_PAIRS_EXAMINED"], row["READ_PAIR_DUPLICATES"])
```

- **Example output:**

```text
lib1 22 0
```

## 6) `dnakit.sequencing_qc.picard_oxog_metrics()`

- **Purpose:** 8-oxoguanine (CpCG to CpCA) oxidation error metrics per sequence context (`CollectOxoGMetrics`).
- **Method:** Picard's defaults apply (context size 1, mapping quality 30, base quality 20, insert size 60 to 600).
- **API:** `dnakit.sequencing_qc.picard_oxog_metrics(input_path, output_path, reference_sequence, *, context_size=None, contexts=None, dbsnp=None, intervals=None, minimum_mapping_quality=None, minimum_quality_score=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and the reference FASTA; optional `contexts`, `dbsnp`, `intervals` and quality filters.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics`: one row per sequence context.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_oxog_metrics

result = picard_oxog_metrics("example.sam", str(out / "oxog.txt"), "ref.fa", verbose=False)
first = result["metrics"][0]
print(result["status"], first["SAMPLE_ALIAS"], len(first["CONTEXT"]), first["TOTAL_SITES"] > 0)
```

- **Example output:**

```text
completed sample1 3 True
```

## 7) `dnakit.sequencing_qc.picard_snvq_yield_metrics()`

- **Purpose:** Yield of bases passing quality thresholds from SNV-quality reads (`CollectQualityYieldMetricsSNVQ`).
- **Method:** The input reads must carry the per-base SNV quality tags `qa`, `qc`, `qg` and `qt` (as written by Ultima Genomics pipelines); without them Picard stops with an error and the result has `status` `failed`.
- **API:** `dnakit.sequencing_qc.picard_snvq_yield_metrics(input_path, output_path, *, alternate_quality_attribute=None, include_bq_histogram=None, include_secondary_alignments=None, include_supplemental_alignments=None, assume_sorted=None, stop_after=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file with the `qa`/`qc`/`qg`/`qt` tags and a new output path.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics`: total and passing-filter reads and bases, Q20/Q30/Q40 bases.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_snvq_yield_metrics

print(picard_snvq_yield_metrics("example.sam", str(out / "plain.txt"), verbose=False)["status"])
row = picard_snvq_yield_metrics("snvq.sam", str(out / "snvq.txt"), verbose=False)["metrics"][0]
print(row["TOTAL_READS"], row["TOTAL_BASES"], row["PF_READS"])
```

- **Example output:**

```text
failed
44 1584 44
```

## 8) `dnakit.sequencing_qc.picard_rnaseq_metrics()`

- **Purpose:** RNA-seq coverage and transcript-region metrics from aligned reads (`CollectRnaSeqMetrics`).
- **Method:** `ref_flat` is a UCSC refFlat gene annotation. `strand_specificity` is `NONE`, `FIRST_READ_TRANSCRIPTION_STRAND` or `SECOND_READ_TRANSCRIPTION_STRAND`. `metric_accumulation_level` takes `ALL_READS` (Picard's default), `SAMPLE`, `LIBRARY` and `READ_GROUP`. `chart_output` needs R.
- **API:** `dnakit.sequencing_qc.picard_rnaseq_metrics(input_path, output_path, ref_flat, strand_specificity, *, ribosomal_intervals=None, chart_output=None, end_bias_bases=None, minimum_length=None, rrna_fragment_percentage=None, ignore_sequence=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path, a refFlat annotation file and `strand_specificity` (`NONE`, `FIRST_READ_TRANSCRIPTION_STRAND` or `SECOND_READ_TRANSCRIPTION_STRAND`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` with the base counts per region and strand statistics.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_rnaseq_metrics

row = picard_rnaseq_metrics("example.sam", str(out / "rna.txt"), "genes.refflat", "NONE", verbose=False)["metrics"][0]
print(row["PF_BASES"], row["CODING_BASES"] + row["UTR_BASES"] + row["INTRONIC_BASES"] + row["INTERGENIC_BASES"] == row["PF_ALIGNED_BASES"])
```

- **Example output:**

```text
1584 True
```

## 9) `dnakit.sequencing_qc.picard_sam_error_metrics()`

- **Purpose:** Base-error metrics stratified by context, cycle, quality and other covariates (`CollectSamErrorMetrics`).
- **Method:** One file `<output_prefix>.<error>_by_<stratifier>` is written per `ERROR[:STRATIFIER...]` request; they are listed in `files`. `error_metrics` entries are added to Picard's default list (`ERROR`, `ERROR:BASE_QUALITY`, `ERROR:INSERT_LENGTH`, ... `INDEL_ERROR`), so repeating a default entry is rejected by Picard.
- **API:** `dnakit.sequencing_qc.picard_sam_error_metrics(input_path, output_prefix, reference_sequence, *, error_metrics=None, vcf=None, intervals=None, interval_iterator=None, min_base_q=None, min_mapping_q=None, ..., verbose=True)`.
- **Input:** A SAM/BAM file, a new output prefix and the reference FASTA; optional `error_metrics`, `vcf` of known variants, `intervals`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `files`: every metrics file written.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from pathlib import Path

from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_sam_error_metrics

result = picard_sam_error_metrics("example.sam", str(out / "err"), "ref.fa", verbose=False)
names = sorted(Path(path).name for path in result["files"])
print(result["status"], len(names) > 20, names[:3])
```

- **Example output:**

```text
completed True ['err.error_by_all', 'err.error_by_base_quality', 'err.error_by_binned_length_homopolymer_and_following_ref_base']
```

## 10) `dnakit.sequencing_qc.picard_sequencing_artifact_metrics()`

- **Purpose:** Library-preparation (pre-adapter) and bait-bias substitution artifacts (`CollectSequencingArtifactMetrics`).
- **Method:** Writes `<prefix>.pre_adapter_detail_metrics`, `.pre_adapter_summary_metrics`, `.bait_bias_detail_metrics`, `.bait_bias_summary_metrics` and `.error_summary_metrics`. Picard's defaults apply (context size 1, mapping quality 30, base quality 20, insert size 60 to 600, duplicates and non-PF reads excluded).
- **API:** `dnakit.sequencing_qc.picard_sequencing_artifact_metrics(input_path, output_prefix, reference_sequence, *, context_size=None, contexts_to_print=None, dbsnp=None, intervals=None, minimum_mapping_quality=None, minimum_quality_score=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output prefix and the reference FASTA.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus five files `<prefix>.*_metrics`; `pre_adapter_summary` parsed.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_sequencing_artifact_metrics

result = picard_sequencing_artifact_metrics("example.sam", str(out / "art"), "ref.fa", verbose=False)
print(sorted(result["outputs"]))
print(len(result["pre_adapter_summary"]) > 0)
```

- **Example output:**

```text
['bait_bias_detail_metrics', 'bait_bias_summary_metrics', 'error_summary_metrics', 'pre_adapter_detail_metrics', 'pre_adapter_summary_metrics']
True
```

## 11) `dnakit.sequencing_qc.picard_wgs_metrics_nonzero_coverage()`

- **Purpose:** Whole-genome coverage metrics, also restricted to bases with non-zero coverage (`CollectWgsMetricsWithNonZeroCoverage`).
- **Method:** Picard marks this tool EXPERIMENTAL. Defaults: coverage cap 250, mapping quality 20, base quality 20, read length 150, sample size 10000. `chart_output` needs R. The `metrics` rows are `WHOLE_GENOME` and `NON_ZERO_REGIONS`.
- **API:** `dnakit.sequencing_qc.picard_wgs_metrics_nonzero_coverage(input_path, output_path, reference_sequence, *, chart_output=None, theoretical_sensitivity_output=None, allele_fraction=None, coverage_cap=None, count_unpaired=None, include_bq_histogram=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a new output path and the reference FASTA; optional `coverage_cap`, `intervals`, `chart_output` (needs R).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics`: the two categories.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_wgs_metrics_nonzero_coverage

rows = {r["CATEGORY"]: r for r in picard_wgs_metrics_nonzero_coverage("example.sam", str(out / "wgs.txt"), "ref.fa", verbose=False)["metrics"]}
print(sorted(rows), rows["WHOLE_GENOME"]["GENOME_TERRITORY"], rows["NON_ZERO_REGIONS"]["MEAN_COVERAGE"] > rows["WHOLE_GENOME"]["MEAN_COVERAGE"])
```

- **Example output:**

```text
['NON_ZERO_REGIONS', 'WHOLE_GENOME'] 2000 True
```

## 12) `dnakit.sequencing_qc.picard_variant_calling_metrics()`

- **Purpose:** Per-sample and aggregate variant-calling metrics of a VCF (`CollectVariantCallingMetrics`).
- **Method:** `dbsnp` is a reference dbSNP file (dbSNP or VCF format, indexed). Writes `<prefix>.variant_calling_detail_metrics` and `<prefix>.variant_calling_summary_metrics`.
- **API:** `dnakit.sequencing_qc.picard_variant_calling_metrics(input_path, dbsnp, output_prefix, *, gvcf_input=None, sequence_dictionary=None, target_intervals=None, thread_count=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** An indexed VCF, an indexed dbSNP VCF and a new output prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `summary_metrics` and `detail_metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_variant_calling_metrics

result = picard_variant_calling_metrics("sample.vcf.gz", "dbsnp.vcf.gz", str(out / "vc"), verbose=False)
row = result["summary_metrics"][0]
print(row["TOTAL_SNPS"], row["NUM_IN_DB_SNP"], row["PCT_DBSNP"], result["detail_metrics"][0]["SAMPLE_ALIAS"])
```

- **Example output:**

```text
8 8 1 sample1
```

## 13) `dnakit.sequencing_qc.picard_arrays_variant_calling_metrics()`

- **Purpose:** Summary and per-sample metrics of a genotyping-array VCF written by GtcToVcf (`CollectArraysVariantCallingMetrics`).
- **Method:** Writes `<prefix>.arrays_variant_calling_detail_metrics`, `.arrays_variant_calling_summary_metrics` and `.arrays_control_code_summary_metrics`. The call-rate threshold for an autocall pass defaults to 0.98.
- **API:** `dnakit.sequencing_qc.picard_arrays_variant_calling_metrics(input_path, dbsnp, output_prefix, *, call_rate_pf_threshold=None, sequence_dictionary=None, num_processors=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** An indexed array VCF, an indexed dbSNP VCF and a new output prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `summary_metrics`, `detail_metrics` and the control-code summary file.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_arrays_variant_calling_metrics

result = picard_arrays_variant_calling_metrics("arrays.vcf.gz", "dbsnp.vcf.gz", str(out / "arr"), verbose=False)
row = result["summary_metrics"][0]
print(row["NUM_ASSAYS"], row["NUM_CALLS"], row["CALL_RATE"], sorted(result["outputs"]))
```

- **Example output:**

```text
8 8 1 ['arrays_control_code_summary_metrics', 'arrays_variant_calling_detail_metrics', 'arrays_variant_calling_summary_metrics']
```

## 14) `dnakit.sequencing_qc.picard_read_group_checksum()`

- **Purpose:** Hash code of the read groups in a SAM/BAM header (`CalculateReadGroupChecksum`).
- **Method:** Adding or removing read groups changes the code. The hash is written to `output_path` and returned as `checksum`.
- **API:** `dnakit.sequencing_qc.picard_read_group_checksum(input_path, output_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A SAM/BAM file and a new output path.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `checksum` (32 hexadecimal characters).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from pathlib import Path

from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_read_group_checksum

first = picard_read_group_checksum("example.sam", str(out / "one.txt"), verbose=False)["checksum"]
Path("renamed.sam").write_text(Path("example.sam").read_text().replace("SM:sample1", "SM:sample2"))
second = picard_read_group_checksum("renamed.sam", str(out / "two.txt"), verbose=False)["checksum"]
print(len(first), first != second)
```

- **Example output:**

```text
32 True
```

## 15) `dnakit.sequencing_qc.picard_check_terminator_block()`

- **Purpose:** Check that the last BGZF block of a block-compressed file such as a BAM is the well-formed terminator (`CheckTerminatorBlock`).
- **Method:** Picard exits with code 100 when the terminator is missing or malformed; that is a verdict, not a failure, so `status` is `completed` and `valid` is False.
- **API:** `dnakit.sequencing_qc.picard_check_terminator_block(input_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A block-compressed file such as a BAM.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `valid` (true when the exit code is 0).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_check_terminator_block

good = picard_check_terminator_block("example.bam", verbose=False)
bad = picard_check_terminator_block("cut.bam", verbose=False)
print(good["valid"], good["return_code"], bad["valid"], bad["return_code"], bad["status"])
```

- **Example output:**

```text
True 0 False 100 completed
```

## 16) `dnakit.sequencing_qc.picard_check_illumina_directory()`

- **Purpose:** Check an Illumina BaseCalls directory for the files every tile and cycle needs (`CheckIlluminaDirectory`).
- **Method:** `read_structure` is the integer/letter description of a run (`B` sample barcode, `M` molecular barcode, `T` template, `S` skip), for example `76T8B76T`. `valid` is true when Picard finds no missing or empty file; otherwise `status` is `failed` and `log_tail` lists what is missing. Files are never created (Picard's `FAKE_FILES` and `LINK_LOCS` are not exposed).
- **API:** `dnakit.sequencing_qc.picard_check_illumina_directory(basecalls_dir, lanes, read_structure, *, data_types=None, tile_numbers=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** The BaseCalls directory, lane numbers and the read structure (for example `76T8B76T`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `valid` (true when no file is missing or empty).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_check_illumina_directory

result = picard_check_illumina_directory("run/Data/Intensities/BaseCalls", [1], "2T", data_types=["BaseCalls"], verbose=False)
print(result["status"], result["valid"])
missing = picard_check_illumina_directory("run/Data/Intensities/BaseCalls", [2], "2T", data_types=["BaseCalls"], verbose=False)
print(missing["status"], missing["valid"])
```

- **Example output:**

```text
completed True
failed False
```

## 17) `dnakit.sequencing_qc.picard_fingerprint_metrics()`

- **Purpose:** Viability statistics of single fingerprints (`CalculateFingerprintMetrics`).
- **Method:** Reports p-values under the null that each fingerprint is a non-contaminated diploid human, plus the discriminatory power. `calculate_by` is `READGROUP` (Picard's default), `LIBRARY`, `SAMPLE` or `FILE`.
- **API:** `dnakit.sequencing_qc.picard_fingerprint_metrics(input_paths, haplotype_map, output_path, *, calculate_by=None, genotype_lod_threshold=None, number_of_sampling=None, reference_sequence=None, validation_stringency=None, picard_path=None, ..., verbose=True)`.
- **Input:** SAM/BAM/CRAM or VCF files, a haplotype map and a new output path; `calculate_by` is `READGROUP`, `LIBRARY`, `SAMPLE` or `FILE`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics`: one row per fingerprint unit.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_fingerprint_metrics

row = picard_fingerprint_metrics(["example.sam"], "hap.map", str(out / "fp.txt"), calculate_by="SAMPLE", verbose=False)["metrics"][0]
print(row["SAMPLE_ALIAS"], row["HAPLOTYPES"], row["HAPLOTYPES_WITH_EVIDENCE"] <= row["HAPLOTYPES"])
```

- **Example output:**

```text
sample1 8 True
```

## 18) `dnakit.sequencing_qc.picard_extract_fingerprint()`

- **Purpose:** Extract the fingerprint genotype likelihoods of a file into a VCF (`ExtractFingerprint`).
- **Method:** The likelihoods are the PL values at the fingerprint sites. `reference_sequence` is required by Picard for this tool.
- **API:** `dnakit.sequencing_qc.picard_extract_fingerprint(input_path, haplotype_map, output_path, reference_sequence, *, sample_alias=None, contamination=None, extract_contamination=None, extract_non_representatives_too=None, locus_max_reads=None, validation_stringency=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM file, a haplotype map, a new output path and the reference FASTA.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus the fingerprint VCF and `sites` (number of sites).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_extract_fingerprint

result = picard_extract_fingerprint("example.sam", "hap.map", str(out / "fp.vcf"), "ref.fa", sample_alias="S9", verbose=False)
print(result["status"], result["sites"], "PL" in (out / "fp.vcf").read_text())
```

- **Example output:**

```text
completed 8 True
```

## 19) `dnakit.sequencing_qc.picard_check_fingerprint()`

- **Purpose:** Check that a SAM/BAM/CRAM or VCF matches the expected sample's known genotypes (`CheckFingerprint`).
- **Method:** `genotypes` is a VCF of known genotypes; `haplotype_map` lists the fingerprint SNPs. Writes `<prefix>.fingerprinting_summary_metrics` and `<prefix>.fingerprinting_detail_metrics`. Exit code 1 (expected sample not found in the genotypes) and 2 (no valid checks) are reported as `failed`.
- **API:** `dnakit.sequencing_qc.picard_check_fingerprint(input_path, genotypes, haplotype_map, output_prefix, *, expected_sample_alias=None, observed_sample_alias=None, genotype_lod_threshold=None, ignore_read_groups=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`.
- **Input:** A SAM/BAM/CRAM or VCF file, a genotypes VCF, a haplotype map and a new output prefix.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `summary_metrics` (LOD of the expected sample) and `detail_metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_check_fingerprint

result = picard_check_fingerprint("example.sam", "genotypes.vcf.gz", "hap.map", str(out / "chk"), verbose=False)
row = result["summary_metrics"][0]
print(row["SAMPLE"], row["HAPLOTYPES_WITH_GENOTYPES"], isinstance(row["LOD_EXPECTED_SAMPLE"], float))
missing = picard_check_fingerprint("example.sam", "genotypes.vcf.gz", "hap.map", str(out / "chk2"), expected_sample_alias="nobody", verbose=False)
print(missing["status"], missing["return_code"])
```

- **Example output:**

```text
sample1 8 True
failed 1
```

## 20) `dnakit.sequencing_qc.picard_crosscheck_fingerprints()`

- **Purpose:** Check that groups of reads or variants come from the same individual (`CrosscheckFingerprints`).
- **Method:** `crosscheck_by` rolls fingerprints up to `READGROUP` (Picard's default), `LIBRARY`, `SAMPLE` or `FILE`. `second_input_paths` compares each sample of the first set with the same sample of the second set (`crosscheck_mode` is `CHECK_SAME_SAMPLE` by default, or `CHECK_ALL_OTHERS`). Picard exits with code 1 when a group relates to the others unexpectedly; the metrics file is still written, so `status` stays `completed` and `mismatch_detected` is true. Exit code 2 (no valid checks) is a failure.
- **API:** `dnakit.sequencing_qc.picard_crosscheck_fingerprints(input_paths, haplotype_map, output_path, *, crosscheck_by=None, second_input_paths=None, crosscheck_mode=None, lod_threshold=None, expect_all_groups_to_match=None, output_errors_only=None, ..., verbose=True)`.
- **Input:** SAM/BAM/CRAM or VCF files, a haplotype map and a new output path; `crosscheck_by` is `READGROUP`, `LIBRARY`, `SAMPLE` or `FILE`.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` (LOD score and RESULT per pair of groups) and `mismatch_detected`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_crosscheck_fingerprints

result = picard_crosscheck_fingerprints(["example.sam"], "hap.map", str(out / "x.txt"), crosscheck_by="SAMPLE", verbose=False)
print(result["status"], result["mismatch_detected"], {row["RESULT"] for row in result["metrics"]})
```

- **Example output:**

```text
completed False {'EXPECTED_MATCH'}
```

## 21) `dnakit.sequencing_qc.picard_cluster_crosscheck_metrics()`

- **Purpose:** Cluster the groups of a `CrosscheckFingerprints` metrics file by LOD score (`ClusterCrosscheckMetrics`).
- **Method:** Two groups are connected when their LOD score exceeds `lod_threshold` (Picard's default is 0.0); the resulting clusters are written to `output_path` and parsed in `metrics`.
- **API:** `dnakit.sequencing_qc.picard_cluster_crosscheck_metrics(input_path, output_path, *, lod_threshold=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`.
- **Input:** A crosscheck metrics file and a new output path.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`, plus `metrics` with a `CLUSTER` and `CLUSTER_SIZE` for each group.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_cluster_crosscheck_metrics, picard_crosscheck_fingerprints

picard_crosscheck_fingerprints(["example.sam"], "hap.map", str(out / "x.txt"), crosscheck_by="SAMPLE", verbose=False)
result = picard_cluster_crosscheck_metrics(str(out / "x.txt"), str(out / "clusters.txt"), verbose=False)
print(result["status"], [(row["CLUSTER"], row["CLUSTER_SIZE"]) for row in result["metrics"]])
```

- **Example output:**

```text
completed [(0, 1)]
```
