# Picard 指标、质控与指纹工具

这一组函数是 `dnakit.sequencing_qc` 中针对 [Picard](https://github.com/broadinstitute/picard)（Broad Institute，MIT 许可）指标、质控和指纹工具的命令行封装。每个函数用明确的 `KEY=value` 参数列表（不经过 shell）运行一个固定的 Picard 工具，拒绝覆盖已存在的输出，并返回包含 `status`、`backend`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典；Picard 写出指标文件时还给出解析后的 `metrics`。只有退出码被接受且所有预期输出都存在时 `status` 才是 `completed`。值为 `None` 的参数不会传给 Picard，所以使用 Picard 自己的默认值。测试版本为 Picard 3.5.0。

Picard 是外部程序。请先安装（`conda install -c bioconda picard`，或 `picard.jar` 加 Java 17），并把 `picard` 放到 PATH、传入 `picard_path`（启动脚本或 `picard.jar`）或设置 `PICARD_JAR`；找不到时函数抛出带安装提示的 `BackendUnavailableError`。画图（`chart_output`）需要 `Rscript`。示例直接从模块导入，所有函数都以 `verbose=False` 调用。

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

**示例文件：** `dbsnp.vcf`

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

下一步由上面的文件派生出其余输入（已建索引的 VCF、BAM、被截断的 BAM、带标签的读段、注释与估计程序输出文件、芯片 VCF 和一个很小的 Illumina 运行目录）。

- **示例代码：**

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

- **示例输出：**

```text
['arrays.vcf.gz', 'arrays.vcf.gz.tbi', 'cut.bam', 'dbsnp.vcf.gz', 'dbsnp.vcf.gz.tbi', 'example.bam', 'genes.refflat', 'genotypes.vcf.gz', 'genotypes.vcf.gz.tbi', 'run', 'sample.vcf.gz', 'sample.vcf.gz.tbi', 'snvq.sam']
```

## 1) `dnakit.sequencing_qc.picard_accumulate_metrics()`

- **作用：** 合并分片的 Picard 指标文件（质量产出或变异识别指标）。
- **处理方法：** 对 QualityYield 文件调用 `AccumulateQualityYieldMetrics`，对变异识别指标前缀调用 `AccumulateVariantCallingMetrics`，把多个分片合成一个文件。
- **API：** `dnakit.sequencing_qc.picard_accumulate_metrics(input_paths, output_path, *, metric_type='quality_yield', picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 指标文件列表（`quality_yield`）或不含扩展名的前缀列表（`variant_calling`），以及新的输出路径或前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`outputs` 和合并后的 `metrics`（或 `summary_metrics`、`detail_metrics`）。
- **示例代码：**

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

- **示例输出：**

```text
completed 16
```

## 2) `dnakit.sequencing_qc.picard_format_contamination_metrics()`

- **作用：** 把污染估计程序（bafRegress 或 VerifyIDIntensity）的输出整理成 Picard 指标文件。
- **处理方法：** `source='bafregress'` 运行 `CreateBafRegressMetricsFile`，`source='verifyidintensity'` 运行 `CreateVerifyIDIntensityContaminationMetricsFile`。
- **API：** `dnakit.sequencing_qc.picard_format_contamination_metrics(input_path, output_prefix, *, source='bafregress', picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 估计程序的标准输出文件和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有指标文件 `<前缀>.bafregress_metrics` 或 `.verifyidintensity_metrics`，解析结果在 `metrics`。
- **示例代码：**

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

- **示例输出：**

```text
completed S1 0.012 900
completed 0.05 -120.5
```

## 3) `dnakit.sequencing_qc.picard_convert_artifact_to_oxog()`

- **作用：** 把测序假象指标转换成 8-oxoG 指标格式。
- **处理方法：** 运行 `ConvertSequencingArtifactToOxoG`，读取 `CollectSequencingArtifactMetrics` 的 pre-adapter 和 bait-bias 详细指标。
- **API：** `dnakit.sequencing_qc.picard_convert_artifact_to_oxog(input_base, output_base, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 假象指标的前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.oxog_metrics`，解析结果在 `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_convert_artifact_to_oxog, picard_sequencing_artifact_metrics

picard_sequencing_artifact_metrics("example.sam", str(out / "art"), "ref.fa", verbose=False)
converted = picard_convert_artifact_to_oxog(str(out / "art"), str(out / "oxo"), verbose=False)
print(converted["status"], len(converted["metrics"]) > 0, sorted(converted["metrics"][0])[:3])
```

- **示例输出：**

```text
completed True ['ALT_NONOXO_BASES', 'ALT_OXO_BASES', 'CONTEXT']
```

## 4) `dnakit.sequencing_qc.picard_base_distribution_by_cycle()`

- **作用：** 按测序循环统计碱基组成。
- **处理方法：** 运行 `CollectBaseDistributionByCycle`，输出每个循环 A/C/G/T/N 的百分比。
- **API：** `dnakit.sequencing_qc.picard_base_distribution_by_cycle(input_path, output_path, *, aligned_reads_only=None, pf_reads_only=None, chart_output=None, assume_sorted=None, stop_after=None, validation_stringency=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件和新的输出路径；可选 `aligned_reads_only`、`pf_reads_only`、`chart_output`（需要 R）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`cycles`（每个读段和循环一行）与 `cycle_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_base_distribution_by_cycle

result = picard_base_distribution_by_cycle("example.sam", str(out / "bd.txt"), verbose=False)
print(result["status"], result["cycle_count"], round(sum(result["cycles"][0][k] for k in ("PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N"))))
```

- **示例输出：**

```text
completed 72 100
```

## 5) `dnakit.sequencing_qc.picard_duplicate_metrics()`

- **作用：** 从已标记重复的文件收集重复指标。
- **处理方法：** 运行 `CollectDuplicateMetrics`；Picard 注明光学重复数恒为 0。
- **API：** `dnakit.sequencing_qc.picard_duplicate_metrics(input_path, metrics_file, *, assume_sorted=None, stop_after=None, validation_stringency=None, reference_sequence=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 已标记重复的 SAM/BAM/CRAM 文件和新的指标文件路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`，含 READ_PAIRS_EXAMINED、READ_PAIR_DUPLICATES、PERCENT_DUPLICATION。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_duplicate_metrics

row = picard_duplicate_metrics("example.sam", str(out / "dup.txt"), verbose=False)["metrics"][0]
print(row["LIBRARY"], row["READ_PAIRS_EXAMINED"], row["READ_PAIR_DUPLICATES"])
```

- **示例输出：**

```text
lib1 22 0
```

## 6) `dnakit.sequencing_qc.picard_oxog_metrics()`

- **作用：** 统计 8-oxoG（CpCG 到 CpCA）氧化错误率。
- **处理方法：** 运行 `CollectOxoGMetrics`，按序列上下文汇报氧化错误率；需要参考序列。
- **API：** `dnakit.sequencing_qc.picard_oxog_metrics(input_path, output_path, reference_sequence, *, context_size=None, contexts=None, dbsnp=None, intervals=None, minimum_mapping_quality=None, minimum_quality_score=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和参考 FASTA；可选 `contexts`、`dbsnp`、`intervals` 与质量过滤参数。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`：每个序列上下文一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_oxog_metrics

result = picard_oxog_metrics("example.sam", str(out / "oxog.txt"), "ref.fa", verbose=False)
first = result["metrics"][0]
print(result["status"], first["SAMPLE_ALIAS"], len(first["CONTEXT"]), first["TOTAL_SITES"] > 0)
```

- **示例输出：**

```text
completed sample1 3 True
```

## 7) `dnakit.sequencing_qc.picard_snvq_yield_metrics()`

- **作用：** 统计带 SNV 质量标签的读段的质量产出。
- **处理方法：** 运行 `CollectQualityYieldMetricsSNVQ`；读段必须带 `qa`、`qc`、`qg`、`qt` 标签，否则 Picard 报错，`status` 为 `failed`。
- **API：** `dnakit.sequencing_qc.picard_snvq_yield_metrics(input_path, output_path, *, alternate_quality_attribute=None, include_bq_histogram=None, include_secondary_alignments=None, include_supplemental_alignments=None, assume_sorted=None, stop_after=None, ..., verbose=True)`。
- **输入：** 带 `qa`/`qc`/`qg`/`qt` 标签的 SAM/BAM/CRAM 文件和新的输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`：总读段/碱基、通过过滤的读段/碱基以及 Q20/Q30/Q40 碱基。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_snvq_yield_metrics

print(picard_snvq_yield_metrics("example.sam", str(out / "plain.txt"), verbose=False)["status"])
row = picard_snvq_yield_metrics("snvq.sam", str(out / "snvq.txt"), verbose=False)["metrics"][0]
print(row["TOTAL_READS"], row["TOTAL_BASES"], row["PF_READS"])
```

- **示例输出：**

```text
failed
44 1584 44
```

## 8) `dnakit.sequencing_qc.picard_rnaseq_metrics()`

- **作用：** 统计 RNA-seq 比对在转录本区域的分布。
- **处理方法：** 运行 `CollectRnaSeqMetrics`，使用 refFlat 注释，把碱基归入编码、UTR、内含子、基因间或核糖体区域。
- **API：** `dnakit.sequencing_qc.picard_rnaseq_metrics(input_path, output_path, ref_flat, strand_specificity, *, ribosomal_intervals=None, chart_output=None, end_bias_bases=None, minimum_length=None, rrna_fragment_percentage=None, ignore_sequence=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径、refFlat 注释文件和 `strand_specificity`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`，含各区域碱基数和链特异性统计。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_rnaseq_metrics

row = picard_rnaseq_metrics("example.sam", str(out / "rna.txt"), "genes.refflat", "NONE", verbose=False)["metrics"][0]
print(row["PF_BASES"], row["CODING_BASES"] + row["UTR_BASES"] + row["INTRONIC_BASES"] + row["INTERGENIC_BASES"] == row["PF_ALIGNED_BASES"])
```

- **示例输出：**

```text
1584 True
```

## 9) `dnakit.sequencing_qc.picard_sam_error_metrics()`

- **作用：** 按上下文、循环、质量等协变量统计碱基错误率。
- **处理方法：** 运行 `CollectSamErrorMetrics`，每个 `ERROR[:分层]` 请求写一个 `<前缀>.<error>_by_<stratifier>` 文件。
- **API：** `dnakit.sequencing_qc.picard_sam_error_metrics(input_path, output_prefix, reference_sequence, *, error_metrics=None, vcf=None, intervals=None, interval_iterator=None, min_base_q=None, min_mapping_q=None, ..., verbose=True)`。
- **输入：** SAM/BAM 文件、新的输出前缀和参考 FASTA；可选 `error_metrics`、已知变异 `vcf`、`intervals`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`files`：写出的所有指标文件。
- **示例代码：**

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

- **示例输出：**

```text
completed True ['err.error_by_all', 'err.error_by_base_quality', 'err.error_by_binned_length_homopolymer_and_following_ref_base']
```

## 10) `dnakit.sequencing_qc.picard_sequencing_artifact_metrics()`

- **作用：** 统计文库制备（pre-adapter）和 bait-bias 的单碱基替换假象。
- **处理方法：** 运行 `CollectSequencingArtifactMetrics`，写出 pre-adapter、bait-bias 的详细和汇总指标以及总错误汇总。
- **API：** `dnakit.sequencing_qc.picard_sequencing_artifact_metrics(input_path, output_prefix, reference_sequence, *, context_size=None, contexts_to_print=None, dbsnp=None, intervals=None, minimum_mapping_quality=None, minimum_quality_score=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出前缀和参考 FASTA。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有五个 `<前缀>.*_metrics` 文件；`pre_adapter_summary` 为解析结果。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_sequencing_artifact_metrics

result = picard_sequencing_artifact_metrics("example.sam", str(out / "art"), "ref.fa", verbose=False)
print(sorted(result["outputs"]))
print(len(result["pre_adapter_summary"]) > 0)
```

- **示例输出：**

```text
['bait_bias_detail_metrics', 'bait_bias_summary_metrics', 'error_summary_metrics', 'pre_adapter_detail_metrics', 'pre_adapter_summary_metrics']
True
```

## 11) `dnakit.sequencing_qc.picard_wgs_metrics_nonzero_coverage()`

- **作用：** 统计全基因组覆盖度，并同时给出覆盖度非零区域的指标。
- **处理方法：** 运行 `CollectWgsMetricsWithNonZeroCoverage`（Picard 标注为实验性），指标行为 `WHOLE_GENOME` 和 `NON_ZERO_REGIONS`。
- **API：** `dnakit.sequencing_qc.picard_wgs_metrics_nonzero_coverage(input_path, output_path, reference_sequence, *, chart_output=None, theoretical_sensitivity_output=None, allele_fraction=None, coverage_cap=None, count_unpaired=None, include_bq_histogram=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、新的输出路径和参考 FASTA；可选 `coverage_cap`、`intervals`、`chart_output`（需要 R）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`：两个类别各一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_wgs_metrics_nonzero_coverage

rows = {r["CATEGORY"]: r for r in picard_wgs_metrics_nonzero_coverage("example.sam", str(out / "wgs.txt"), "ref.fa", verbose=False)["metrics"]}
print(sorted(rows), rows["WHOLE_GENOME"]["GENOME_TERRITORY"], rows["NON_ZERO_REGIONS"]["MEAN_COVERAGE"] > rows["WHOLE_GENOME"]["MEAN_COVERAGE"])
```

- **示例输出：**

```text
['NON_ZERO_REGIONS', 'WHOLE_GENOME'] 2000 True
```

## 12) `dnakit.sequencing_qc.picard_variant_calling_metrics()`

- **作用：** 统计 VCF 的每样本和总体变异识别指标。
- **处理方法：** 运行 `CollectVariantCallingMetrics`，需要已建索引的 dbSNP 参考文件。
- **API：** `dnakit.sequencing_qc.picard_variant_calling_metrics(input_path, dbsnp, output_prefix, *, gvcf_input=None, sequence_dictionary=None, target_intervals=None, thread_count=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 已建索引的 VCF、已建索引的 dbSNP VCF 和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary_metrics` 与 `detail_metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_variant_calling_metrics

result = picard_variant_calling_metrics("sample.vcf.gz", "dbsnp.vcf.gz", str(out / "vc"), verbose=False)
row = result["summary_metrics"][0]
print(row["TOTAL_SNPS"], row["NUM_IN_DB_SNP"], row["PCT_DBSNP"], result["detail_metrics"][0]["SAMPLE_ALIAS"])
```

- **示例输出：**

```text
8 8 1 sample1
```

## 13) `dnakit.sequencing_qc.picard_arrays_variant_calling_metrics()`

- **作用：** 统计基因分型芯片（GtcToVcf 产生）VCF 的汇总和每样本指标。
- **处理方法：** 运行 `CollectArraysVariantCallingMetrics`，VCF 头部必须有 Infinium 控制信息行。
- **API：** `dnakit.sequencing_qc.picard_arrays_variant_calling_metrics(input_path, dbsnp, output_prefix, *, call_rate_pf_threshold=None, sequence_dictionary=None, num_processors=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 已建索引的芯片 VCF、已建索引的 dbSNP VCF 和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary_metrics`、`detail_metrics` 以及控制码汇总文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_arrays_variant_calling_metrics

result = picard_arrays_variant_calling_metrics("arrays.vcf.gz", "dbsnp.vcf.gz", str(out / "arr"), verbose=False)
row = result["summary_metrics"][0]
print(row["NUM_ASSAYS"], row["NUM_CALLS"], row["CALL_RATE"], sorted(result["outputs"]))
```

- **示例输出：**

```text
8 8 1 ['arrays_control_code_summary_metrics', 'arrays_variant_calling_detail_metrics', 'arrays_variant_calling_summary_metrics']
```

## 14) `dnakit.sequencing_qc.picard_read_group_checksum()`

- **作用：** 计算 SAM/BAM 头部读组的哈希值。
- **处理方法：** 运行 `CalculateReadGroupChecksum`；增删读组会改变哈希值。
- **API：** `dnakit.sequencing_qc.picard_read_group_checksum(input_path, output_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** SAM/BAM 文件和新的输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`checksum`（32 位十六进制）。
- **示例代码：**

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

- **示例输出：**

```text
32 True
```

## 15) `dnakit.sequencing_qc.picard_check_terminator_block()`

- **作用：** 检查块压缩文件（如 BAM）最后一个 BGZF 块是否是完整的终止块。
- **处理方法：** 运行 `CheckTerminatorBlock`；退出码 100 表示终止块缺失或损坏，这是判定结果而不是失败。
- **API：** `dnakit.sequencing_qc.picard_check_terminator_block(input_path, *, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 块压缩文件，例如 BAM。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`valid`（退出码为 0 时为真）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_check_terminator_block

good = picard_check_terminator_block("example.bam", verbose=False)
bad = picard_check_terminator_block("cut.bam", verbose=False)
print(good["valid"], good["return_code"], bad["valid"], bad["return_code"], bad["status"])
```

- **示例输出：**

```text
True 0 False 100 completed
```

## 16) `dnakit.sequencing_qc.picard_check_illumina_directory()`

- **作用：** 检查 Illumina BaseCalls 目录中每个 tile 和循环所需的文件是否齐全。
- **处理方法：** 运行 `CheckIlluminaDirectory`；不暴露会创建文件的 `FAKE_FILES` 与 `LINK_LOCS`。
- **API：** `dnakit.sequencing_qc.picard_check_illumina_directory(basecalls_dir, lanes, read_structure, *, data_types=None, tile_numbers=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** BaseCalls 目录、泳道编号和读段结构（例如 `76T8B76T`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`valid`（没有缺失或空文件时为真）。
- **示例代码：**

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

- **示例输出：**

```text
completed True
failed False
```

## 17) `dnakit.sequencing_qc.picard_fingerprint_metrics()`

- **作用：** 统计单个指纹的可用性。
- **处理方法：** 运行 `CalculateFingerprintMetrics`，报告指纹相对于“无污染二倍体人类”零假设的 p 值与区分能力。
- **API：** `dnakit.sequencing_qc.picard_fingerprint_metrics(input_paths, haplotype_map, output_path, *, calculate_by=None, genotype_lod_threshold=None, number_of_sampling=None, reference_sequence=None, validation_stringency=None, picard_path=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 或 VCF 文件列表、单倍型图和新的输出路径；`calculate_by` 取 `READGROUP`、`LIBRARY`、`SAMPLE`、`FILE`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`：每个指纹单位一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_fingerprint_metrics

row = picard_fingerprint_metrics(["example.sam"], "hap.map", str(out / "fp.txt"), calculate_by="SAMPLE", verbose=False)["metrics"][0]
print(row["SAMPLE_ALIAS"], row["HAPLOTYPES"], row["HAPLOTYPES_WITH_EVIDENCE"] <= row["HAPLOTYPES"])
```

- **示例输出：**

```text
sample1 8 True
```

## 18) `dnakit.sequencing_qc.picard_extract_fingerprint()`

- **作用：** 把文件在指纹位点的基因型似然（PL）提取成 VCF。
- **处理方法：** 运行 `ExtractFingerprint`；Picard 要求提供参考序列。
- **API：** `dnakit.sequencing_qc.picard_extract_fingerprint(input_path, haplotype_map, output_path, reference_sequence, *, sample_alias=None, contamination=None, extract_contamination=None, extract_non_representatives_too=None, locus_max_reads=None, validation_stringency=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 文件、单倍型图、新的输出路径和参考 FASTA。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有指纹 VCF 与位点数 `sites`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_extract_fingerprint

result = picard_extract_fingerprint("example.sam", "hap.map", str(out / "fp.vcf"), "ref.fa", sample_alias="S9", verbose=False)
print(result["status"], result["sites"], "PL" in (out / "fp.vcf").read_text())
```

- **示例输出：**

```text
completed 8 True
```

## 19) `dnakit.sequencing_qc.picard_check_fingerprint()`

- **作用：** 检查文件的身份是否与预期样本的已知基因型一致。
- **处理方法：** 运行 `CheckFingerprint`；退出码 1（预期样本不在基因型文件中）和 2（无有效检查）按失败处理。
- **API：** `dnakit.sequencing_qc.picard_check_fingerprint(input_path, genotypes, haplotype_map, output_prefix, *, expected_sample_alias=None, observed_sample_alias=None, genotype_lod_threshold=None, ignore_read_groups=None, reference_sequence=None, validation_stringency=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 或 VCF 文件、基因型 VCF、单倍型图和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary_metrics`（预期样本的 LOD）与 `detail_metrics`。
- **示例代码：**

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

- **示例输出：**

```text
sample1 8 True
failed 1
```

## 20) `dnakit.sequencing_qc.picard_crosscheck_fingerprints()`

- **作用：** 检查多组读段或变异是否来自同一个体。
- **处理方法：** 运行 `CrosscheckFingerprints`；退出码 1 表示有组之间的关系不符合预期，指标文件仍会写出，所以 `status` 仍为 `completed`，同时 `mismatch_detected` 为真。
- **API：** `dnakit.sequencing_qc.picard_crosscheck_fingerprints(input_paths, haplotype_map, output_path, *, crosscheck_by=None, second_input_paths=None, crosscheck_mode=None, lod_threshold=None, expect_all_groups_to_match=None, output_errors_only=None, ..., verbose=True)`。
- **输入：** SAM/BAM/CRAM 或 VCF 文件列表、单倍型图和新的输出路径；`crosscheck_by` 取 `READGROUP`、`LIBRARY`、`SAMPLE`、`FILE`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`（每对组的 LOD 与 RESULT）和 `mismatch_detected`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_crosscheck_fingerprints

result = picard_crosscheck_fingerprints(["example.sam"], "hap.map", str(out / "x.txt"), crosscheck_by="SAMPLE", verbose=False)
print(result["status"], result["mismatch_detected"], {row["RESULT"] for row in result["metrics"]})
```

- **示例输出：**

```text
completed False {'EXPECTED_MATCH'}
```

## 21) `dnakit.sequencing_qc.picard_cluster_crosscheck_metrics()`

- **作用：** 按 LOD 对 `CrosscheckFingerprints` 的结果聚类。
- **处理方法：** 运行 `ClusterCrosscheckMetrics`；LOD 超过 `lod_threshold`（默认 0）的两组相连。
- **API：** `dnakit.sequencing_qc.picard_cluster_crosscheck_metrics(input_path, output_path, *, lod_threshold=None, picard_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 交叉检查指标文件和新的输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`metrics`，每个组带有 `CLUSTER` 与 `CLUSTER_SIZE`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.sequencing_qc.picard_metrics_p08_cli_functions import picard_cluster_crosscheck_metrics, picard_crosscheck_fingerprints

picard_crosscheck_fingerprints(["example.sam"], "hap.map", str(out / "x.txt"), crosscheck_by="SAMPLE", verbose=False)
result = picard_cluster_crosscheck_metrics(str(out / "x.txt"), str(out / "clusters.txt"), verbose=False)
print(result["status"], [(row["CLUSTER"], row["CLUSTER_SIZE"]) for row in result["metrics"]])
```

- **示例输出：**

```text
completed [(0, 1)]
```
