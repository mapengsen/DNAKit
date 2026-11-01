# Picard 的 SAM/BAM、Illumina 运行目录、VCF 与区间工具

这一组函数是 `dnakit.alignment and dnakit.variants` 中针对 [Picard](https://github.com/broadinstitute/picard)（Broad Institute，MIT 许可）的命令行封装，涉及读写 SAM/BAM/CRAM、Illumina 运行目录、VCF/BCF、区间列表、FASTA 和 GFF3 的工具。每个函数用明确的 `KEY=value` 参数列表（不经过 shell）运行一个固定的 Picard 工具，拒绝覆盖已存在的输出，并返回包含 `status`、`backend`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail` 以及少量解析结果的字典。只有退出码被接受且所有预期输出都存在时 `status` 才是 `completed`。值为 `None` 的参数不会传给 Picard，所以使用 Picard 自己的默认值。测试版本为 Picard 3.5.0。Picard 的 JavaScript 过滤（`FilterSamReads includeJavascript`、`FilterVcf JAVASCRIPT_FILE`）会运行调用者提供的脚本，因此不开放。

Picard 是外部程序。请先安装（`conda install -c bioconda picard`，或 `picard.jar` 加 Java 17），并把 `picard` 放到 PATH、传入 `picard_path`（启动脚本或 `picard.jar`）或设置 `PICARD_JAR`；找不到时函数抛出带安装提示的 `BackendUnavailableError`。示例直接从模块导入，所有函数都以 `verbose=False` 调用。

## 示例数据

示例读取下面这些很小的文件。函数从不覆盖已有文件，所以每个输出都写入新建的临时目录。

**示例文件：** `example.sam`

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

**示例文件：** `ref.fa`

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

**示例文件：** `ref.fa.fai`

```text
chr1	1200	6	60	61
chr2	800	1232	60	61
```

**示例文件：** `ref.dict`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1200
@SQ	SN:chr2	LN:800
```

**示例文件：** `hap.map`

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

**示例文件：** `genotypes.vcf`

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

**示例文件：** `sample.vcf`

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

下一步由上面的文件派生出其余输入（已建索引的 VCF、BAM、被截断的 BAM、带标签的读段、注释与估计程序输出文件、芯片 VCF 和一个很小的 Illumina 运行目录）。

- **示例代码：**

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

- **示例输出：**

```text
['a.interval_list', 'a_to_b.chain', 'b.interval_list', 'genotypes.vcf.gz', 'run', 'sample.vcf.gz', 'x.gff3']
```

## 1) `dnakit.alignment.picard_downsample_sam()`

- **作用：** 按模板名保留确定性的随机读段子集。
- **处理方法：** 运行 `DownsampleSam`；同一模板（配对、次级和补充比对）一起保留或丢弃。
- **API：** `dnakit.alignment.picard_downsample_sam(input_path, output_path, *, probability=None, strategy=None, random_seed=None, accuracy=None, metrics_file=None, reference_sequence=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和 `probability`；可选 `strategy`、`random_seed`、`metrics_file`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有给出 `metrics_file` 时的下采样结果指标 `metrics`。
- **示例代码：**

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

- **示例输出：**

```text
completed 22 True True
True
```

## 2) `dnakit.alignment.picard_position_based_downsample_sam()`

- **作用：** 按流动槽位置对读段做确定性的下采样。
- **处理方法：** 运行 `PositionBasedDownsampleSam`，由读名中的 tile 与 x/y 坐标决定保留哪些模板（不是随机的）。
- **API：** `dnakit.alignment.picard_position_based_downsample_sam(input_path, output_path, fraction, *, read_name_regex=None, remove_duplicate_information=None, allow_multiple_downsampling_despite_warnings=None, stop_after=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和 0 到 1 之间的 `fraction`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 中的下采样文件。
- **示例代码：**

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

- **示例输出：**

```text
True True
```

## 3) `dnakit.alignment.picard_filter_sam_reads()`

- **作用：** 按一种过滤条件筛选读段。
- **处理方法：** 运行 `FilterSamReads`；支持已比对/未比对、读名列表、区间和标签值过滤，不提供 JavaScript 过滤。
- **API：** `dnakit.alignment.picard_filter_sam_reads(input_path, output_path, filter, *, read_list_file=None, interval_list=None, tag=None, tag_values=None, sort_order=None, write_reads_files=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和 `filter`，以及该过滤所需的 `read_list_file`、`interval_list` 或 `tag`/`tag_values`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 中的筛选结果文件。
- **示例代码：**

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

- **示例输出：**

```text
completed 2
completed 0
```

## 4) `dnakit.alignment.picard_split_sam_by_number_of_reads()`

- **作用：** 把按查询名分组的文件拆成多个读段数相近的文件。
- **处理方法：** 运行 `SplitSamByNumberOfReads`；`split_to_n_files` 与 `split_to_n_reads` 二选一。
- **API：** `dnakit.alignment.picard_split_sam_by_number_of_reads(input_path, output_dir, *, split_to_n_files=None, split_to_n_reads=None, out_prefix=None, total_reads_in_input=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、已存在的输出目录，以及 `split_to_n_files` 与 `split_to_n_reads` 之一。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`files`：新生成的分片文件。
- **示例代码：**

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

- **示例输出：**

```text
['shard_0001.sam', 'shard_0002.sam', 'shard_0003.sam'] 44
```

## 5) `dnakit.alignment.picard_reorder_sam()`

- **作用：** 按新参考序列的 contig 顺序重排读段。
- **处理方法：** 运行 `ReorderSam`，按名称精确匹配 contig；新参考中不存在的 contig 上的读段变为未比对。
- **API：** `dnakit.alignment.picard_reorder_sam(input_path, output_path, sequence_dictionary, *, allow_contig_length_discordance=None, allow_incomplete_dict_concordance=None, reference_sequence=None, validation_stringency=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和序列字典文件（SAM/BAM/CRAM、VCF、区间列表、FASTA 或 `.dict`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 中重排后的文件。
- **示例代码：**

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

- **示例输出：**

```text
completed chr2 chr1
```

## 6) `dnakit.alignment.picard_revert_sam()`

- **作用：** 把已比对的文件还原为未比对状态。
- **处理方法：** 运行 `RevertSam`：去掉比对信息和重复标记、恢复原始质量（OQ）、默认按查询名排序。
- **API：** `dnakit.alignment.picard_revert_sam(input_path, output_path, *, output_by_readgroup=None, output_by_readgroup_file_format=None, remove_alignment_information=None, remove_duplicate_information=None, restore_original_qualities=None, restore_hardclips=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件和新的输出路径（`output_by_readgroup=True` 时为已存在的目录）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有按读组输出时的 `files`。
- **示例代码：**

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

- **示例输出：**

```text
completed 44 True False
['rg1.sam']
```

## 7) `dnakit.alignment.picard_mark_illumina_adapters()`

- **作用：** 在按查询名排序的读段中查找 Illumina 接头并写入 `XT` 标签。
- **处理方法：** 运行 `MarkIlluminaAdapters`；写出被剪切碱基数的直方图。
- **API：** `dnakit.alignment.picard_mark_illumina_adapters(input_path, metrics_path, *, output_path=None, adapters=None, five_prime_adapter=None, three_prime_adapter=None, max_error_rate_pe=None, max_error_rate_se=None, ..., verbose=True)`。
- **输入：** 按查询名排序的 SAM/BAM/CRAM 文件和新的指标路径；可选 `output_path`、`adapters`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs`，找到接头时还有 `histogram`。
- **示例代码：**

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

- **示例输出：**

```text
completed ['metrics', 'output'] 44
```

## 8) `dnakit.alignment.picard_merge_bam_alignment()`

- **作用：** 把比对器输出与未比对 BAM 合并，恢复比对器丢失的标签。
- **处理方法：** 运行 `MergeBamAlignment`；未比对文件须按查询名排序且只含未比对读段。
- **API：** `dnakit.alignment.picard_merge_bam_alignment(unmapped_bam, reference_sequence, output_path, *, aligned_bams=None, read1_aligned_bams=None, read2_aligned_bams=None, sort_order=None, add_mate_cigar=None, clip_adapters=None, ..., verbose=True)`。
- **输入：** 未比对 BAM、参考 FASTA、新的输出路径，以及 `aligned_bams`（或 `read1_aligned_bams` 与 `read2_aligned_bams`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 中合并后的文件。
- **示例代码：**

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

- **示例输出：**

```text
completed 44 44 {'rg1'}
```

## 9) `dnakit.alignment.picard_mark_duplicates_with_mate_cigar()`

- **作用：** 在带 `MC` 标签的坐标排序读段中标记重复分子。
- **处理方法：** `umi_aware=False` 运行 `SimpleMarkDuplicatesWithMateCigar`（Picard 标注为实验性），`umi_aware=True` 运行 `UmiAwareMarkDuplicatesWithMateCigar` 并按 UMI 分组。
- **API：** `dnakit.alignment.picard_mark_duplicates_with_mate_cigar(input_paths, output_path, metrics_file, *, umi_aware=False, umi_metrics_file=None, umi_tag_name=None, max_edit_distance_to_join=None, allow_missing_umis=None, remove_duplicates=None, ..., verbose=True)`。
- **输入：** 带 `MC` 标签的坐标排序文件、新的输出路径和指标路径；`umi_aware=True` 时还需 `umi_metrics_file`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`（以及 `umi_metrics`）。
- **示例代码：**

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

- **示例输出：**

```text
22 2 4
completed 4
```

## 10) `dnakit.alignment.picard_illumina_basecalls_to_fastq()`

- **作用：** 从 Illumina 运行目录的原始碱基调用生成 FASTQ。
- **处理方法：** 运行 `IlluminaBasecallsToFastq`；`output_prefix` 不拆分样本，`multiplex_params` 按条形码拆分。
- **API：** `dnakit.alignment.picard_illumina_basecalls_to_fastq(basecalls_dir, lanes, read_structure, run_barcode, *, output_prefix=None, multiplex_params=None, machine_name=None, flowcell_barcode=None, read_name_format=None, compress_outputs=None, ..., verbose=True)`。
- **输入：** BaseCalls 目录、泳道编号、读段结构、运行条形码，以及 `output_prefix`（或 `multiplex_params`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`files`：写出的 FASTQ 文件。
- **示例代码：**

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

- **示例输出：**

```text
['s.1.fastq'] ['AT', 'CG', 'GC', 'TA'] @M1:run1:FC1:1:1101:2000:3000 :N
```

## 11) `dnakit.alignment.picard_illumina_basecalls_to_sam()`

- **作用：** 把 Illumina 原始碱基调用转换成按条形码拆分的未比对 SAM/BAM/CRAM。
- **处理方法：** 运行 `IlluminaBasecallsToSam`；单样本用 `output_path`、`sample_alias`、`library_name`，多样本用 `library_params`。
- **API：** `dnakit.alignment.picard_illumina_basecalls_to_sam(basecalls_dir, lanes, read_structure, run_barcode, sequencing_center, *, output_path=None, sample_alias=None, library_name=None, library_params=None, platform=None, read_group_id=None, ..., verbose=True)`。
- **输入：** BaseCalls 目录、泳道编号、读段结构、运行条形码、测序中心，以及 `output_path` 与 `sample_alias`（或 `library_params`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`files`：写出的未比对文件。
- **示例代码：**

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

- **示例输出：**

```text
s1 l1 XX ['AT', 'CG', 'GC', 'TA']
```

## 12) `dnakit.variants.picard_filter_vcf()`

- **作用：** 对 VCF 的基因型和位点施加硬过滤。
- **处理方法：** 运行 `FilterVcf`；`min_dp`、`min_gq` 过滤基因型，`min_ab`、`min_qd`、`max_fs` 过滤位点。不提供 JavaScript 过滤。
- **API：** `dnakit.variants.picard_filter_vcf(input_path, output_path, *, min_dp=None, min_gq=None, min_ab=None, min_qd=None, max_fs=None, picard_path=None, ..., verbose=True)`。
- **输入：** VCF/BCF 文件、新的输出路径和阈值。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`records` 与索引路径 `index`。
- **示例代码：**

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

- **示例输出：**

```text
8 [True, True, True, True, False, False, False, False]
['LowQD', 'LowQD', 'LowQD', 'LowQD', 'PASS', 'PASS', 'PASS', 'PASS']
```

## 13) `dnakit.variants.picard_liftover_vcf()`

- **作用：** 用 UCSC chain 文件把 VCF 转换到另一个参考版本。
- **处理方法：** 运行 `LiftoverVcf`；输出按目标坐标排序并建索引，无法转换的记录写入 `reject_path`。
- **API：** `dnakit.variants.picard_liftover_vcf(input_path, output_path, reject_path, chain, reference_sequence, *, liftover_min_match=None, recover_swapped_ref_alt=None, allow_missing_fields_in_header=None, tags_to_drop=None, tags_to_reverse=None, write_original_position=None, ..., verbose=True)`。
- **输入：** VCF/BCF、新的输出和拒绝文件路径、chain 文件、目标参考 FASTA（需 `.fai` 与 `.dict`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`lifted_records`、`rejected_records` 与 `index`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_liftover_vcf

result = picard_liftover_vcf("sample.vcf", str(out / "lifted.vcf"), str(out / "rej.vcf"), "a_to_b.chain", "target.fa", write_original_position=True, verbose=False)
rows = [line.split("\t") for line in (out / "lifted.vcf").read_text().splitlines() if not line.startswith("#")]
print(result["lifted_records"], result["rejected_records"], {row[0] for row in rows}, "OriginalContig=chr1" in rows[0][7])
```

- **示例输出：**

```text
8 0 {'chrA'} True
```

## 14) `dnakit.variants.picard_make_sites_only_vcf()`

- **作用：** 去掉 VCF 的全部基因型列，保留位点级信息。
- **处理方法：** 运行 `MakeSitesOnlyVcf`，AN、AF 等注释保留。
- **API：** `dnakit.variants.picard_make_sites_only_vcf(input_path, output_path, *, samples=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** VCF/VCF.gz/BCF 和新的输出路径；可选 `samples`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`records` 与 `index`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_make_sites_only_vcf

result = picard_make_sites_only_vcf("sample.vcf.gz", str(out / "sites.vcf"), verbose=False)
rows = [line.split("\t") for line in (out / "sites.vcf").read_text().splitlines() if not line.startswith("#")]
print(result["records"], {len(row) for row in rows})
```

- **示例输出：**

```text
8 {8}
```

## 15) `dnakit.variants.picard_make_vcf_sample_name_map()`

- **作用：** 生成样本名到 VCF 路径的映射表。
- **处理方法：** 运行 `MakeVcfSampleNameMap`；Picard 3.5.0 在第一列写路径、第二列写样本名，`entries` 按输入路径区分两列。
- **API：** `dnakit.variants.picard_make_vcf_sample_name_map(input_paths, output_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 单样本 VCF/GVCF 文件列表和新的输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`entries`：样本名/路径对。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_make_vcf_sample_name_map

result = picard_make_vcf_sample_name_map(["sample.vcf.gz", "genotypes.vcf.gz"], str(out / "map.tsv"), verbose=False)
print([(entry["sample"], Path(entry["path"]).name) for entry in result["entries"]])
```

- **示例输出：**

```text
[('sample1', 'sample.vcf.gz'), ('sample1', 'genotypes.vcf.gz')]
```

## 16) `dnakit.variants.picard_rename_sample_in_vcf()`

- **作用：** 重命名单样本 VCF/BCF 中的样本。
- **处理方法：** 运行 `RenameSampleInVcf`；给出 `old_sample_name` 时会断言现有名称与之相同。
- **API：** `dnakit.variants.picard_rename_sample_in_vcf(input_path, output_path, new_sample_name, *, old_sample_name=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 单样本 VCF/BCF、新的输出路径和 `new_sample_name`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 中改名后的文件。
- **示例代码：**

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

- **示例输出：**

```text
completed NEW
failed
```

## 17) `dnakit.variants.picard_split_vcfs()`

- **作用：** 把 VCF/BCF 拆成表头相同的 SNP 文件和插入缺失文件。
- **处理方法：** 运行 `SplitVcfs`；两个输出都建索引，`strict=True`（默认）时遇到其它类型记录报错。
- **API：** `dnakit.variants.picard_split_vcfs(input_path, snp_output, indel_output, *, strict=None, sequence_dictionary=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** VCF/BCF 以及新的 SNP、indel 输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`snp_records` 与 `indel_records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_split_vcfs

result = picard_split_vcfs("mixed.vcf", str(out / "snp.vcf"), str(out / "indel.vcf"), verbose=False)
print(result["snp_records"], result["indel_records"])
```

- **示例输出：**

```text
8 1
```

## 18) `dnakit.variants.picard_liftover_haplotype_map()`

- **作用：** 用 chain 文件把指纹单倍型数据库转换到另一个参考。
- **处理方法：** 运行 `LiftOverHaplotypeMap`，把序列字典写入输出表头。
- **API：** `dnakit.variants.picard_liftover_haplotype_map(input_path, output_path, chain, sequence_dictionary, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 单倍型图、新的输出路径、chain 文件和序列字典。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`snps`：转换后的 SNP 数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_liftover_haplotype_map

result = picard_liftover_haplotype_map("hap.map", str(out / "hap2.map"), "a_to_b.chain", "target.dict", verbose=False)
rows = [line.split("\t") for line in (out / "hap2.map").read_text().splitlines() if not line.startswith(("@", "#"))]
print(result["snps"], {row[0] for row in rows})
```

- **示例输出：**

```text
8 {'chrA'}
```

## 19) `dnakit.variants.picard_interval_list_to_bed()`

- **作用：** 把 Picard 区间列表转换成 BED。
- **处理方法：** 运行 `IntervalListToBed`；区间列表为 1 起始闭区间，BED 为 0 起始半开区间。
- **API：** `dnakit.variants.picard_interval_list_to_bed(input_path, output_path, *, score=None, sort=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 区间列表和新的 BED 路径；可选 `score`、`sort`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`intervals`：区间数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_interval_list_to_bed

result = picard_interval_list_to_bed("a.interval_list", str(out / "a.bed"), score=7, verbose=False)
print(result["intervals"])
print((out / "a.bed").read_text().strip().replace("\t", " "))
```

- **示例输出：**

```text
3
chr1 99 200 a 7 +
chr1 299 400 b 7 +
chr2 9 50 c 7 +
```

## 20) `dnakit.variants.picard_interval_list_tools()`

- **作用：** 对区间列表做并、交、差、对称差、重叠选择与拼接。
- **处理方法：** 运行 `IntervalListTools`；`action` 为 `CONCAT`、`UNION`、`INTERSECT`、`SUBTRACT`、`SYMDIFF` 或 `OVERLAPS`。
- **API：** `dnakit.variants.picard_interval_list_tools(input_paths, output_path, *, action='CONCAT', second_input_paths=None, padding=None, invert=None, sort=None, unique=None, ..., verbose=True)`。
- **输入：** 区间列表（或 VCF）、新的输出路径、`action`，`SUBTRACT`/`SYMDIFF`/`OVERLAPS` 还需 `second_input_paths`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`intervals`，给出 `count_value` 时还有 `count`。
- **示例代码：**

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

- **示例输出：**

```text
CONCAT 4
OVERLAPS 1 1
SUBTRACT bases 212
```

## 21) `dnakit.variants.picard_scatter_intervals_by_ns()`

- **作用：** 把参考序列按 N 与 ACGT 交替拆成区间。
- **处理方法：** 运行 `ScatterIntervalsByNs`；FASTA 需要 `.fai` 和 `.dict`。
- **API：** `dnakit.variants.picard_scatter_intervals_by_ns(reference_sequence, output_path, *, output_type=None, max_to_merge=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 参考 FASTA（有索引和字典）和新的区间列表路径；可选 `output_type`、`max_to_merge`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`intervals`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_scatter_intervals_by_ns

result = picard_scatter_intervals_by_ns("gap.fa", str(out / "n.il"), verbose=False)
print([tuple(line.split("\t")[1:3]) + (line.split("\t")[4],) for line in (out / "n.il").read_text().splitlines() if not line.startswith("@")])
```

- **示例输出：**

```text
[('1', '100', 'ACGTmer'), ('101', '120', 'Nmer'), ('121', '220', 'ACGTmer')]
```

## 22) `dnakit.variants.picard_normalize_fasta()`

- **作用：** 把 FASTA 每条序列的行长统一（每条记录最后一行除外）。
- **处理方法：** 运行 `NormalizeFasta`；默认每行 100 个碱基。
- **API：** `dnakit.variants.picard_normalize_fasta(input_path, output_path, *, line_length=None, truncate_sequence_names_at_whitespace=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** FASTA、新的输出路径和 `line_length`；可选 `truncate_sequence_names_at_whitespace`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.picard_vcf_intervals_p08_cli_functions import picard_normalize_fasta

result = picard_normalize_fasta("messy.fa", str(out / "n.fa"), line_length=16, truncate_sequence_names_at_whitespace=True, verbose=False)
print(result["records"])
print((out / "n.fa").read_text().strip())
```

- **示例输出：**

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

- **作用：** 按坐标排序 GFF3，使其可以建索引。
- **处理方法：** 运行 `SortGff`；给出序列字典时按字典排 contig，否则按字典序。
- **API：** `dnakit.variants.picard_sort_gff(input_path, output_path, *, sequence_dictionary=None, n_records_in_memory=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** GFF3 和新的输出路径；可选 `sequence_dictionary`、`n_records_in_memory`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`features`。
- **示例代码：**

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

- **示例输出：**

```text
['ID=g1a', 'ID=g1b', 'ID=g2']
['ID=g2', 'ID=g1a', 'ID=g1b']
```
