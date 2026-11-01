# modkit：修饰碱基 pileup、modBAM 工具、motif 与差异甲基化

`dnakit.epigenome` 中的这一组函数是对 [modkit](https://github.com/nanoporetech/modkit) 的命令行封装。modkit 是 Oxford Nanopore 处理 modBAM（`MM`、`ML` 标签）中碱基修饰信息的工具。函数涵盖：bedMethyl pileup（`modkit pileup`、`pileup-hemi`）、read 级别表格（`extract`）、modBAM 的汇总与质检（`modbam summary`、`sample-probs`、`check-tags`）、modBAM 编辑（`adjust-mods`、`call-mods`、`update-tags`、`repair`）、bedMethyl 工具（`bedmethyl merge`、`tobigwig`、`map-to-genome`、`stats`、`localize`）、序列 motif（`motif bed`、`search`、`evaluate`、`refine`）、甲基化熵、对真值的验证、差异甲基化（`dmr`）以及 6mA 开放染色质预测。每个函数用明确的参数列表（不经过 shell）运行一个固定子命令，拒绝覆盖已有输出，返回含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`operation` 和 `provenance` 的字典；表格类结果还有 `row_count` 和很短的 `rows` 预览。只有 modkit 退出码为 0 且所有预期输出都存在时 `status` 才是 `completed`；modkit 报错时为 `failed`，原因在 `log_tail` 里。

测试版本为 modkit 0.6.4。它采用 Oxford Nanopore Technologies Public License v1.0 许可，DNAKit 把它作为独立程序启动，不复制其代码。请从项目的 GitHub releases 页面安装发布版二进制（或 `conda install -c bioconda modkit`）并放到 PATH 上，也可以传入 `modkit_path`；程序缺失时函数会在创建任何输出之前抛出带安装提示的 `BackendUnavailableError`。值为 `None` 或 `False` 的选项不会传给 modkit，因此使用 modkit 自己的默认值。位置参数总是放在选项之前，因为 modkit 的多值选项（如 `--modified-bases`）会把后面的位置参数吞掉。示例对所有函数都使用 `verbose=False`，并用 pysam 把下面的 SAM 文本转换成小的 modBAM 文件。

## 示例数据

示例使用下面这些很小的文件：含两条序列的参考、三组带 CpG 5mC 判读的 100 bp reads（`MM:Z:C+m?` 与 `ML` 标签；`other.sam` 修饰判读很少，`duplex.sam` 还带有半甲基化分析所需的对侧链 `G-m?` 判读）、含两条转录本的参考以及两组比对到转录组的 reads，还有很小的 BED、染色体长度和 GTF 文件。下面的准备步骤会排序并索引 BAM 文件，并把 bgzip 压缩、带 tabix 索引的 bedMethyl 文件写入 `data` 文件夹。其他输出都写入新建的临时文件夹，因为这些函数从不覆盖文件。

**示例文件：** `ref.fa`

```text
>chr1
GCGCTCCTCCGGATTCGAAAATTCAACGCACGTTTCGACCGTCGTCCGGTAGTCGAGACGTCGTAATTCGAACACGTTGCGCTACGTCGAGTCGATTTCGAACACGTACGTCCTATCTCCGGATTCCGGACACGTAACCGCGCAGGGCGCACGCAATGTCATCTCCTCGACCACGTTTATTCGACGCTGGCGCGCACACGTAGCCCCGGATCGATCAAACGTGCATCGATGCCGGAAGCG
>chr2
ATTCGGCACGGACGTCGATCGATATGCGCACGATCGAGTTTGCGCATCGATCGATTTCGAAGAGGACGTCGCCGGTCGTCCGGCTCACGTCGCCACGTGAAGTTCGAAATGCGCCTTCCGGGTGATTCGAAGTTCGAAGATGACGCGCCAGCTTCGAAAA
```

含两条序列（240 bp 和 160 bp）的参考。

**示例文件：** `reads.sam`

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

chr1 上 8 条、chr2 上 4 条正反链 reads；80% 的 CpG 判读为已修饰。

**示例文件：** `other.sam`

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

同样的布局，但只有约 20% 的 CpG 判读为已修饰（第二个条件）。

**示例文件：** `duplex.sam`

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

同一分子上同时带有正链 `C+m?` 与负链 `G-m?` 判读的 reads。

**示例文件：** `tx.fa`

```text
>tx1
ACGTTGTCTTCGAACCAATCGAACCCCCGGTCCGGTAGTCCGGGCGAATTCGAATACGTTTCACGTTTTCGAAACACGTAATTCGCATGGTTTCGATGAAACGTATACATACCCCGGGGGACGTCCGTGACGTCGCGCAGGCGCCCCAGCGCTAGATTCG
>tx2
ACGTTATTCGAACCTCGCTAACGTCGTGCGCGCGCTTCGAATCCGGAGATCTGCCGGCCCAGCGCCCGCAAGACGTGTTACGTCCTTCGAACAAAGCGCATCCTCGAACTACGTGTAATG
```

转录组参考，含转录本 `tx1`（160 nt）和 `tx2`（120 nt）。

**示例文件：** `tx_a.sam`

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

比对到转录本的 reads，大多已修饰（条件 A）。

**示例文件：** `tx_b.sam`

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

比对到转录本的 reads，大多未修饰（条件 B）。

**示例文件：** `sizes.txt`

```text
chr1	240
chr2	160
```

`ref.fa` 的染色体长度。

**示例文件：** `regions.bed`

```text
chr1	0	120	reg1
chr1	120	240	reg2
chr2	0	160	reg3
```

三个用于汇总或比较的区域。

**示例文件：** `features.bed`

```text
chr1	60	61	f1	0	+
chr1	180	181	f2	0	+
```

用于 `localize` 的两个点状特征。

**示例文件：** `truth.bed`

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

用于 `validate` 的真值位置：已修饰（`m`）与未修饰（`-`）。

**示例文件：** `genes.gtf`

```text
chr1	test	transcript	1001	1280	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	exon	1001	1080	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	exon	1201	1280	.	+	.	gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";
chr1	test	transcript	1001	1360	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
chr1	test	exon	1001	1060	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
chr1	test	exon	1301	1360	.	+	.	gene_id "g1"; transcript_id "tx2"; gene_name "GENE1";
```

基因 `g1` 及其两个转录本 `tx1` 和 `tx2`（外显子总长 160 和 120）。

**示例代码：**

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

**示例输出：**

```text
a completed 116
b completed 119
tx_a completed 77
tx_b completed 79
```

准备步骤：得到排序并索引的 BAM 文件，以及基因组 reads（`a`、`b`）和转录本 reads（`tx_a`、`tx_b`）的 CpG bedMethyl pileup（bgzip 文件，带 tabix 索引，位于 `data`）。

## 1) `dnakit.epigenome.modkit_bedmethyl_map_to_genome()`

- **作用：** 把比对到转录组的 bedMethyl 中某个转录本的记录换算到基因组坐标。
- **处理方法：** 运行 `modkit bedmethyl map-to-genome <bedmethyl> <输出> --transcript-id ... --gtf ...`（按需加 `--transcript-version`、`--ignore-version`、`--header`）。modkit 在 GTF 中找到该转录本的外显子，把每个转录本位置（负链时还有链方向）换算为基因组位置。输出存在时 `status` 为 `completed`。
- **API：** `modkit_bedmethyl_map_to_genome(bedmethyl_path: 'str', output_path: 'str', gtf_path: 'str', transcript_id: 'str', *, transcript_version: 'int | None' = None, ignore_version: 'bool' = False, header: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是 bgzip 压缩、带 tabix 索引、序列名为转录本 ID 的 bedMethyl；`output_path` 是新文件；`gtf_path` 是注释；`transcript_id` 是要换算的转录本；`transcript_version` 与 `ignore_version` 处理版本后缀；`header` 写表头。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 43
chr1 1010 1011 +
```

## 2) `dnakit.epigenome.modkit_bedmethyl_merge()`

- **作用：** 对两个或更多 bedMethyl 文件做外连接合并，重叠记录的计数相加。
- **处理方法：** 运行 `modkit bedmethyl merge <文件> <文件>... --genome-sizes ... --out-bed ...`（以及 `--header`、`--min-samples`、`--min-sample-coverage`、`--threads`）。同一位置、链和修饰代码的记录计数相加；出现在少于 `min_samples` 个样本中或覆盖度太低的记录被丢弃。
- **API：** `modkit_bedmethyl_merge(bedmethyl_paths: 'Sequence[str]', genome_sizes: 'str', output_path: 'str', *, header: 'bool' = False, min_samples: 'int | None' = None, min_sample_coverage: 'int | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_paths` 列出至少两个带 tabix 索引的 bgzip bedMethyl 文件；`genome_sizes` 是两列的染色体长度文件；`output_path` 是新文件；`min_samples` 与 `min_sample_coverage` 是可选整数。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 122
['chr1', '2', '3', 'm', '2', '-']
```

## 3) `dnakit.epigenome.modkit_bedmethyl_to_bigwig()`

- **作用：** 把 bedMethyl 表转换成修饰比例的 BigWig 轨道。
- **处理方法：** 运行 `modkit bedmethyl tobigwig <bedmethyl> <输出> --sizes ...`（或 `--header <bam>`），每个代码一个 `--mod-codes`，按需加 `--negative-strand-values`、`--nzooms`、`--uncompressed`。给出多个代码时 modkit 每个代码写一个轨道，并把代码追加到文件名上。
- **API：** `modkit_bedmethyl_to_bigwig(bedmethyl_path: 'str', output_path: 'str', mod_codes: 'Sequence[str]', *, chrom_sizes: 'str | None' = None, header_bam: 'str | None' = None, negative_strand_values: 'bool' = False, nzooms: 'int | None' = None, uncompressed: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是未压缩的 bedMethyl 文件（modkit 在这里不接受 bgzip 文件）；`output_path` 是新文件；`mod_codes` 是修饰代码（如 `['m']`）；染色体长度由 `chrom_sizes`（两列文件）或 `header_bam` 提供，二者给且只给一个。
- **输出：** 通用结果字典。
- **示例代码：**

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

- **示例输出：**

```text
completed frac.bw True
```

## 4) `dnakit.epigenome.modkit_dmr()`

- **作用：** 检验区域在两个样本之间（pair）或多个样本两两之间（multi）的差异甲基化。
- **处理方法：** `mode='pair'` 运行 `modkit dmr pair -a <对照> -b <实验> --regions-bed ... --ref ... --out-path ...`；`mode='multi'` 对每个样本给出 `--sample <文件> <名称>`，运行 `modkit dmr multi`，并把每一对的 BED 写入新目录。两者都使用 `--base`（默认 C），可选 `--assign-code`、`--single-code`、`--mask`、`--min-valid-coverage`、`--header`，pair 模式还有 `--segment`。modkit 用修饰比例的贝叶斯模型给每个区域打分，score 列表示差异大小。
- **API：** `modkit_dmr(samples: 'dict[str, str]', regions_bed: 'str', reference: 'str', output_path: 'str', *, mode: 'str' = 'pair', base: 'Sequence[str]' = ('C',), assign_code: 'Sequence[str] | None' = None, single_code: 'str | None' = None, mask: 'bool' = False, min_valid_coverage: 'int | None' = None, segment: 'str | None' = None, reference_sample: 'str | None' = None, prefix: 'str | None' = None, header: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `samples` 把样本名映射到带 tabix 索引的 bgzip bedMethyl（名称只能含字母、数字、`_`、`.`、`-`）；pair 模式必须恰好两个（第一个是对照 `-a`），multi 模式至少两个。`regions_bed` 与 `reference`（带 .fai 的 FASTA）必填。`output_path` 在 pair 模式是新 BED 文件，在 multi 模式是新目录（可选 `prefix` 和 `reference_sample`）。
- **输出：** 通用结果字典，另有 `mode`，以及 `record_count`（pair）或 `files`（multi）。
- **示例代码：**

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

- **示例输出：**

```text
completed pair 3
chr1 0 120 25.70579668561382 .
chr1 120 240 25.18549683798696 .
chr2 0 160 22.296727753583554 .
completed ['a_b.bed']
```

## 5) `dnakit.epigenome.modkit_dmr_compare_tx_sites()`

- **作用：** 比较两个条件之间每个转录本位点的修饰水平。
- **处理方法：** 运行 `modkit dmr compare-tx-sites -a <A> -b <B> --out ... --gtf ...`（别名 `compare-transcriptome-sites`），可选 `--full`、`--ignore-version`、`--min-valid-coverage`、`--single-code`、`--plot`、`--top-k`、`--min-effect-size`、`--plot-title`。每个位点用差异模型检验，表中给出得分、p 值以及转录本所属基因。
- **API：** `modkit_dmr_compare_tx_sites(condition_a: 'str', condition_b: 'str', gtf_path: 'str', output_path: 'str', *, min_valid_coverage: 'int | None' = None, single_code: 'str | None' = None, full: 'bool' = False, ignore_version: 'bool' = False, plot_path: 'str | None' = None, top_k: 'int | None' = None, min_effect_size: 'float | None' = None, plot_title: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `condition_a` 与 `condition_b` 是比对到转录组、带 tabix 索引的 bgzip bedMethyl；`gtf_path` 是注释；`output_path` 是新表格；`plot_path` 是可选的新绘图文件（配合 `top_k`、`min_effect_size`、`plot_title`）。
- **输出：** 通用结果字典，另有 `columns`、`row_count`、`rows` 与 `plot`。
- **示例代码：**

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

- **示例输出：**

```text
completed 72
['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'p_value', 'gene_id', 'gene_name']
['chr1', '1001', '1002', 'm'] ['g1.0', 'GENE1']
```

## 6) `dnakit.epigenome.modkit_dmr_isoform()`

- **作用：** 用一个比对到转录组的 bedMethyl 文件比较基因各转录本异构体之间的甲基化。
- **处理方法：** 运行 `modkit dmr isoform <bedmethyl> <输出> --gtf ...`，可选 `--gene-id`、`--gene-name`、`--gene-strand`、`--gene-chrom`、`--min-valid-coverage`、`--single-code`、`--full`、`--ignore-version`、`--plot`、`--min-score`、`--max-pval`。映射到同一基因组位置的不同异构体的位点用差异甲基化模型比较。
- **API：** `modkit_dmr_isoform(bedmethyl_path: 'str', gtf_path: 'str', output_path: 'str', *, gene_id: 'str | None' = None, gene_name: 'str | None' = None, gene_strand: 'str | None' = None, gene_chrom: 'str | None' = None, min_valid_coverage: 'int | None' = None, single_code: 'str | None' = None, full: 'bool' = False, ignore_version: 'bool' = False, plot_path: 'str | None' = None, min_score: 'float | None' = None, max_pval: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是比对到转录组、带 tabix 索引的 bgzip bedMethyl，`gtf_path` 是注释，`output_path` 是新表格；`gene_id`/`gene_name` 选择一个基因（否则检验所有基因）；`plot_path` 是可选的新 SVG 图（针对一个基因）。
- **输出：** 通用结果字典，另有 `columns`、`row_count`、`rows` 与 `plot`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_dmr_isoform

work = Path(tempfile.mkdtemp())
result = modkit_dmr_isoform("data/tx_a.bed.gz", "genes.gtf", str(work / "isoform.tsv"), verbose=False)
print(result["status"], result["row_count"])
print(result["columns"][:7], result["columns"][-3:])
```

- **示例输出：**

```text
completed 4
['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'p_value'] ['n_transcripts', 'gene_id', 'gene_name']
```

## 7) `dnakit.epigenome.modkit_entropy()`

- **作用：** 计算连续修饰位点窗口上的甲基化熵。
- **处理方法：** 运行 `modkit entropy --in-bam <bam> --ref ... --out-bed ...`，位点选择用 `--cpg`、`--motif` 或 `--base`，另有 `--num-positions`（默认 4 个位点）、`--window-size`（默认 50 bp）、`--min-coverage`、`--combine-strands`、`--regions` 及过滤选项。对每个含 `num_positions` 个修饰位点的窗口，按甲基化模式给 reads 分组，报告模式频率的香农熵；熵为 0 表示所有 reads 模式相同。
- **API：** `modkit_entropy(bam_paths: 'Sequence[str]', reference: 'str', output_path: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, base: 'str | None' = None, regions_bed: 'str | None' = None, num_positions: 'int | None' = None, window_size: 'int | None' = None, min_coverage: 'int | None' = None, max_filtered_positions: 'int | None' = None, combine_strands: 'bool' = False, mask: 'bool' = False, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, num_reads: 'int | None' = None, header: 'bool' = False, drop_zeros: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_paths` 列出已排序并索引的 modBAM（多个时合并）；`reference` 是带 .fai 的 FASTA；用 `cpg`、`motif`（+`motif_offset`）或 `base` 之一选择位点；`output_path` 是新的 BED 文件，给出 `regions_bed` 时是新目录，每个区域一个 BED（文件名以 `prefix` 开头）。
- **输出：** 通用结果字典，另有 `record_count`（BED 模式）或 `files`（目录模式）。
- **示例代码：**

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

- **示例输出：**

```text
completed 24
['#chrom', 'start', 'end', 'entropy', 'strand', 'num_reads']
['chr1', '42', '55', '0.3060986', '+', '3']
```

## 8) `dnakit.epigenome.modkit_extract()`

- **作用：** 把 modBAM 中 read 级别的修饰信息提取为表格：概率、硬判读或每条 read 的计数。
- **处理方法：** `table='full'` 运行 `modkit extract full <bam> <输出>`（每条 read 每个位置的修饰概率）；`table='calls'` 运行 `modkit extract calls`（用与 pileup 相同的通过阈值做出判读，`pass_only` 只保留通过的判读）；`table='read-stats'` 运行 `modkit extract read-stats`（每条 read 一行 CSV，给出每种代码的已修饰与未修饰计数，必须给 `--mod-codes`）。属于其他 table 的选项会抛出 `ConfigurationError`，不会被忽略。
- **API：** `modkit_extract(bam_path: 'str', output_path: 'str', *, table: 'str' = 'full', reference: 'str | None' = None, region: 'str | None' = None, include_bed: 'str | None' = None, exclude_bed: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, annotate_motifs: 'bool' = False, mask: 'bool' = False, ignore: 'str | None' = None, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, mapped_only: 'bool' = False, allow_non_primary: 'bool' = False, num_reads: 'int | None' = None, kmer_size: 'int | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, pass_only: 'bool' = False, sampling_frac: 'float | None' = None, mod_codes: 'Sequence[str] | None' = None, ignore_implicit: 'bool' = False, ignore_index: 'bool' = False, bgzf: 'bool' = False, no_headers: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM；`output_path` 是新文件；`table` 为 `full`、`calls` 或 `read-stats`；`reference`（`cpg`、`motif`、`annotate_motifs` 需要）、`region`、`include_bed`、`exclude_bed`、`mapped_only`、`edge_filter`、`kmer_size`、`bgzf`、`no_headers` 用于选择和格式化行；`filter_percentile`、`filter_threshold`、`mod_thresholds`、`no_filtering`、`pass_only` 适用于 `calls`；`read-stats` 必须给 `mod_codes`，如 `['C:m']`。
- **输出：** 通用结果字典，另有 `table`、`columns`、`row_count` 与 `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_extract

work = Path(tempfile.mkdtemp())
for table, extra in (("full", {}), ("calls", {"pass_only": True}), ("read-stats", {"mod_codes": ["C:m"]})):
    result = modkit_extract("reads.bam", str(work / f"{table}.tsv"), table=table, verbose=False, **extra)
    print(table, result["status"], result["row_count"], result["columns"][:3])
```

- **示例输出：**

```text
full completed 205 ['read_id', 'forward_read_position', 'ref_position']
calls completed 185 ['read_id', 'forward_read_position', 'ref_position']
read-stats completed 12 ['read_id', 'chrom', 'aln_start']
```

## 9) `dnakit.epigenome.modkit_localize()`

- **作用：** 对以基因组特征（如启动子、结合位点）为中心的窗口汇总修饰水平。
- **处理方法：** 运行 `modkit localize <bedmethyl> --regions ... --genome-sizes ... --out-file ...`（别名 `localise`），可选 `--window`（每侧默认 2000 bp）、`--stranded`、`--stranded-features`、`--min-coverage`、`--chart` 与 `--name`。窗口位置从每个特征的中心起算，各特征的计数按位置累加。
- **API：** `modkit_localize(bedmethyl_path: 'str', regions_bed: 'str', genome_sizes: 'str', output_path: 'str', *, window: 'int | None' = None, stranded: 'str | None' = None, stranded_features: 'bool' = False, min_coverage: 'int | None' = None, chart_path: 'str | None' = None, chart_name: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是带 tabix 索引的 bgzip bedMethyl；`regions_bed` 是特征；`genome_sizes` 是两列的染色体长度文件；`output_path` 是新表格；`window`（每侧 bp 数）、`stranded`（`true`、`false` 或 `both`）、`min_coverage`；`chart_path` 是可选的新 HTML 图。
- **输出：** 通用结果字典，另有 `columns`、`row_count`、`rows` 与 `chart`。
- **示例代码：**

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

- **示例输出：**

```text
completed 20 localize.html
['mod_code', 'offset', 'n_valid', 'n_mod', 'percent_modified']
```

## 10) `dnakit.epigenome.modkit_modbam_adjust_mods()`

- **作用：** 改写 modBAM 中的碱基修饰判读：丢弃某个代码、把一个代码转换成另一个，或保留/丢弃某些序列上下文。
- **处理方法：** 运行 `modkit modbam adjust-mods <输入> <输出>`，可选 `--ignore`、`--convert <原> <新>`、`--cpg` 或 `--motif`（配合 `--discard-motifs`）、`--edge-filter`、`--invert-edge-filter`、`--filter-probs`、`--mapped-only`、`--num-reads`、通过阈值选项和 `--output-sam`。把 `h` 转换为 `m` 会把 5hmC 判读并入 5mC；忽略某个代码会移除其判读并重新归一化概率。
- **API：** `modkit_modbam_adjust_mods(bam_path: 'str', output_bam: 'str', *, ignore: 'str | None' = None, convert_from: 'str | None' = None, convert_to: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, discard_motifs: 'bool' = False, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, filter_probs: 'bool' = False, mapped_only: 'bool' = False, num_reads: 'int | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, output_sam: 'bool' = False, fast_forward: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM，`output_bam` 是新文件；`ignore` 是一个修饰代码；`convert_from` 与 `convert_to` 必须同时给出；`cpg` / `motif`（+`motif_offset`）与 `discard_motifs` 选择上下文（不需要参考）；`edge_filter` 为 `N` 或 `start,end`；`output_sam` 写 SAM 而不是 BAM。
- **输出：** 通用结果字典。
- **示例代码：**

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

- **示例输出：**

```text
completed
C+m?,0,0 -> C+h?,0,0
```

## 11) `dnakit.epigenome.modkit_modbam_call_mods()`

- **作用：** 对 modBAM 的修饰概率设阈值，写出只含硬判读的新 modBAM。
- **处理方法：** 运行 `modkit modbam call-mods <输入> <输出>`，可选通过阈值选项（`--filter-percentile` 默认 10、`--filter-threshold`、`--mod-threshold`、`--no-filtering`）、抽样（`--num-reads`、`--sampling-frac`、`--seed`）、`--edge-filter`、`--cpg` / `--motif`（+`--discard-motifs`）和 `--output-sam`。判读为修饰的碱基概率设为 100%，判读为未修饰的设为 0%。
- **API：** `modkit_modbam_call_mods(bam_path: 'str', output_bam: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, discard_motifs: 'bool' = False, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, output_sam: 'bool' = False, fast_forward: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM，`output_bam` 是新文件；`filter_threshold` 与 `filter_percentile` 设置通过阈值；`mod_thresholds` 的元素形如 `m:0.8`；`output_sam` 写 SAM。
- **输出：** 通用结果字典。
- **示例代码：**

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

- **示例输出：**

```text
completed
[0, 255]
```

## 12) `dnakit.epigenome.modkit_modbam_check_tags()`

- **作用：** 校验 BAM 中的 MM 与 ML 修饰标签。
- **处理方法：** 运行 `modkit modbam check-tags <bam> --out-dir <新目录>`，可选 `--permissive`、`--region`、`--head`、`--mapped-only`、`--allow-non-primary`、`--ignore-index`、`--prefix`、`--threads`。modkit 检查每条 read 的 ML 值个数与 MM 条目是否一致，写出所发现的修饰和 MM 头部的表格以及失败记录表。未设置 `permissive` 时发现无效记录会以退出码 1 结束，此时状态为 `failed`。
- **API：** `modkit_modbam_check_tags(bam_path: 'str', output_dir: 'str', *, permissive: 'bool' = False, region: 'str | None' = None, num_reads: 'int | None' = None, head: 'int | None' = None, mapped_only: 'bool' = False, allow_non_primary: 'bool' = False, ignore_index: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM；`output_dir` 是存放表格的新目录；`permissive`、`region`、`head`（read 数）、`mapped_only`、`allow_non_primary` 和 `prefix` 为可选。
- **输出：** 通用结果字典，另有 `files`，即 modkit 写出的表格。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import modkit_modbam_check_tags

work = Path(tempfile.mkdtemp())
result = modkit_modbam_check_tags("reads.bam", str(work / "checks"), verbose=False)
print(result["status"], [Path(p).name for p in result["files"]])
print((work / "checks" / "modified_bases.tsv").read_text().splitlines())
```

- **示例输出：**

```text
completed ['modified_bases.tsv', 'valid_mm_headers.tsv']
['strand\tprimary_base\tmod_code\tmode', '+\tC\tm\t?']
```

## 13) `dnakit.epigenome.modkit_modbam_repair()`

- **作用：** 把供体 modBAM 中正确的 MM/ML 标签复制到受体 BAM 的 reads 上。
- **处理方法：** 运行 `modkit modbam repair --donor-bam <供体> --acceptor-bam <受体> --output-bam <输出>`。两个 BAM 都必须按 read 名排序，供体必须包含受体的所有 reads（供体可以有多余的 reads）。典型用途：恢复被修剪或比对步骤丢掉的标签。
- **API：** `modkit_modbam_repair(donor_bam: 'str', acceptor_bam: 'str', output_bam: 'str', *, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `donor_bam` 带有标签，`acceptor_bam` 是缺少标签的 reads，`output_bam` 是新文件；`threads` 可选。
- **输出：** 通用结果字典。
- **示例代码：**

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

- **示例输出：**

```text
completed
12 reads have MM again
```

## 14) `dnakit.epigenome.modkit_modbam_sample_probs()`

- **作用：** 估计 modBAM 中碱基修饰概率的分布以及由此得到的通过阈值。
- **处理方法：** 运行 `modkit modbam sample-probs <bam> --out-dir <新目录>`，可选 `--percentiles`（默认 10,50,90）、`--quantiles`、`--hist`、`--reference`、`--motif`、`--region`、`--include-bed`、`--edge-filter`、`--use-non-primary`、抽样选项和 `--prefix`。modkit 抽样 reads（默认 10042 条），写出 `probabilities.tsv` 与 `thresholds.tsv`，使用 `--hist` 时还有 HTML 直方图。
- **API：** `modkit_modbam_sample_probs(bam_path: 'str', output_dir: 'str', *, percentiles: 'str | None' = None, quantiles: 'str | None' = None, histogram: 'bool' = False, reference: 'str | None' = None, motif: 'str | None' = None, motif_offset: 'int' = 0, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, use_non_primary: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, no_sampling: 'bool' = False, seed: 'int | None' = None, ignore_index: 'bool' = False, prefix: 'str | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM，`output_dir` 是新目录；`percentiles` / `quantiles` 是逗号分隔的数字；`histogram` 额外输出 HTML 直方图；`no_sampling` 读取所有 reads。
- **输出：** 通用结果字典，另有 `files`。
- **示例代码：**

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

- **示例输出：**

```text
completed ['thresholds.tsv']
 base  percentile  threshold  mod_qual
 C     10          0.8066406  206
 C     50          0.8925781  228
 C     90          0.9824219  251
```

## 15) `dnakit.epigenome.modkit_modbam_summary()`

- **作用：** 汇总 modBAM 中存在的修饰标签及基本判读统计。
- **处理方法：** 运行 `modkit modbam summary <bam>`，可选 `--tsv`、`--combine-mods`、`--cpg` / `--motif`（需 `--reference`）、`--region`、`--include-bed`、`--edge-filter`、`--mapped-only`、`--matched-only`、抽样和通过阈值选项。默认输出是汇总块（以 `#` 开头的行）加上按碱基和代码统计判读的表格；`--tsv` 写成一张制表符分隔表。标准输出被流式写入新文件。
- **API：** `modkit_modbam_summary(bam_path: 'str', output_path: 'str', *, tsv: 'bool' = False, combine_mods: 'bool' = False, cpg: 'bool' = False, reference: 'str | None' = None, motif: 'str | None' = None, motif_offset: 'int' = 0, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, mapped_only: 'bool' = False, matched_only: 'bool' = False, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, filter_quantile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, ignore_index: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是 modBAM，`output_path` 是新文件；`tsv`、`combine_mods`、`cpg`/`motif`（需要 `reference`）、`region`、`edge_filter`、`num_reads`（默认抽样 10042 条）、`filter_threshold`、`mod_thresholds` 为可选。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`。
- **示例代码：**

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

- **示例输出：**

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

- **作用：** 把小写的 `Mm`/`Ml` 标签改名为 `MM`/`ML`，并设置修饰标签的模式标志。
- **处理方法：** 运行 `modkit modbam update-tags <输入> <输出>`，可选 `--mode explicit|implicit`、`--no-implicit-probs`、`--output-sam`、`--threads`。`explicit` 写 `?` 标志，`implicit` 写 `.` 标志（未列出的碱基视为未修饰）；modkit 默认拒绝对隐式模式标签做 pileup，因此这一步能让旧 BAM 可用。
- **API：** `modkit_modbam_update_tags(bam_path: 'str', output_bam: 'str', *, mode: 'str | None' = None, no_implicit_probs: 'bool' = False, output_sam: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是带修饰标签的 BAM，`output_bam` 是新文件；`mode` 为 `explicit` 或 `implicit`；`no_implicit_probs` 不添加隐式的未修饰判读；`output_sam` 写 SAM。
- **输出：** 通用结果字典。
- **示例代码：**

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

- **示例输出：**

```text
completed
C+m.,0,0
```

## 17) `dnakit.epigenome.modkit_motif_bed()`

- **作用：** 写出 FASTA 中某个序列 motif 的所有位置的 BED 文件。
- **处理方法：** 运行 `modkit motif bed <fasta> <motif> <offset>`（可加 `--mask`），把标准输出流式写入新文件。motif 可使用 IUPAC 简并码；`offset` 是区间所标记的 motif 内 0 起始位置，因此 `CG` 加偏移 0 标记每个 CpG 的 C；同一位点会以链 `-` 在 G 的位置（对侧链的胞嘧啶）再写一次。
- **API：** `modkit_motif_bed(fasta_path: 'str', motif: 'str', offset: 'int', output_path: 'str', *, mask: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `fasta_path` 是 FASTA，`motif` 仅含字母（允许 IUPAC），`offset` 是小于 motif 长度的整数，`output_path` 是新文件；`mask=True` 尊重软屏蔽。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 138
chr1	1	2	.	.	+
chr1	2	3	.	.	-
```

## 18) `dnakit.epigenome.modkit_motif_evaluate()`

- **作用：** 根据 bedMethyl 表计算给定 motif 的富集统计。
- **处理方法：** 运行 `modkit motif evaluate --in-bedmethyl ... --ref ... --out ...`，可选 `--known-motif <序列> <偏移> <修饰代码>`（可重复）、`--known-motifs-table`、`--contig`、`--min-coverage`、`--context-size` 及 `--low-thresh` / `--high-thresh` 修饰比例阈值。对每个 motif，modkit 统计匹配该 motif 的基因组位点中修饰比例高、中、低的数量，并报告富集的对数几率。
- **API：** `modkit_motif_evaluate(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, contig: 'str | None' = None, context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是带 tabix 索引的 bgzip bedMethyl，`reference` 是带 .fai 的 FASTA，`output_table` 是新 TSV；`known_motifs`（`[序列, 偏移, 修饰代码]` 三元组）与 `known_motifs_table` 至少给一个。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`。
- **示例代码：**

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

- **示例输出：**

```text
completed ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count', 'log_odds']
['m', 'CG', '0', '0.85057473', '74', '13', '23', '0']
```

## 19) `dnakit.epigenome.modkit_motif_refine()`

- **作用：** 用 bedMethyl 表细化一组已知 motif。
- **处理方法：** 运行 `modkit motif refine --in-bedmethyl ... --ref ... --out ...`，可选 `--known-motif`、`--known-motifs-table`、`--min-refine-frac-mod`、`--min-refine-sites`、`--context-size`、`--min-coverage`、`--min-sites`、`--min-log-odds`、`--min-frac-mod` 以及修饰比例阈值。每个已知 motif 用能提高其富集度的序列上下文进行扩展或收窄。
- **API：** `modkit_motif_refine(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, contig: 'str | None' = None, context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, min_sites: 'int | None' = None, min_log_odds: 'float | None' = None, min_frac_mod: 'float | None' = None, min_refine_frac_mod: 'float | None' = None, min_refine_sites: 'int | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是带 tabix 索引的 bgzip bedMethyl，`reference` 是带 .fai 的 FASTA，`output_table` 是新 TSV；`known_motifs` 与 `known_motifs_table` 至少给一个。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`。
- **示例代码：**

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

- **示例输出：**

```text
completed ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count', 'log_odds']
['m', 'C', '0', '0.85057473', '74', '13', '23', '0']
```

## 20) `dnakit.epigenome.modkit_motif_search()`

- **作用：** 根据 bedMethyl 表在基因组中搜索修饰富集的序列 motif。
- **处理方法：** 运行 `modkit motif search --in-bedmethyl ... --ref ... --out-table ...`，可选 `--contig`、`--mod-code`、`--known-motif`、`--known-motifs-table`、`--context-size`、`--init-context-size`、`--min-coverage`、`--min-sites`、`--min-log-odds`、`--min-frac-mod` 及阈值。modkit 从位点多半被修饰的短种子出发，逐步延长成更长的 motif，并保留位点数和对数几率足够的结果。
- **API：** `modkit_motif_search(bedmethyl_path: 'str', reference: 'str', output_table: 'str', *, contig: 'str | None' = None, mod_codes: 'Sequence[str] | None' = None, known_motifs: 'Sequence[Sequence[str]] | None' = None, known_motifs_table: 'str | None' = None, context_size: 'Sequence[int] | None' = None, init_context_size: 'Sequence[int] | None' = None, min_coverage: 'int | None' = None, min_sites: 'int | None' = None, min_log_odds: 'float | None' = None, min_frac_mod: 'float | None' = None, low_threshold: 'float | None' = None, high_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是带 tabix 索引的 bgzip bedMethyl，`reference` 是带 .fai 的 FASTA，`output_table` 是新 TSV；`mod_codes` 限定修饰代码；`known_motifs` 是 `[序列, 偏移, 修饰代码]` 三元组，会与发现的 motif 一并报告。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`。
- **示例代码：**

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

- **示例输出：**

```text
completed 1 ['mod_code', 'motif', 'offset', 'frac_mod', 'high_count', 'low_count', 'mid_count']
['m', 'C', '0', '0.85057473']
```

## 21) `dnakit.epigenome.modkit_open_chromatin_predict()`

- **作用：** 根据 modBAM 中外源 6mA 信号预测开放染色质区域（Fiber-seq 类数据）。
- **处理方法：** 运行 `modkit open-chromatin predict <bam> --model <模型目录> --output <bed>`，可选 `--region`、`--include-bed`、`--min-coverage`、`--threshold`、`--batch-size`、`--step-size`、`--no-header`、`--dryrun`。模型是包含 `model_config.json` 与权重的目录（modkit 发布包的 `models/` 下带有一个）。`dry_run=True` 时 modkit 加载模型、检查输入后直接退出，不做预测，留下空输出文件。
- **API：** `modkit_open_chromatin_predict(bam_path: 'str', model_path: 'str', output_path: 'str', *, region: 'str | None' = None, include_bed: 'str | None' = None, min_coverage: 'int | None' = None, threshold: 'float | None' = None, batch_size: 'int | None' = None, step_size: 'int | None' = None, no_header: 'bool' = False, dry_run: 'bool' = False, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是带 6mA 判读的 modBAM，`model_path` 是已有的模型目录，`output_path` 是新 BED 文件；`threshold` 为得分阈值，`min_coverage` 为最小 read 覆盖度；`dry_run` 只检查输入和模型。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 0
```

## 22) `dnakit.epigenome.modkit_pileup()`

- **作用：** 把各基因组位置的碱基修饰判读汇总成 bedMethyl 文件。
- **处理方法：** 运行 `modkit pileup <bam> <输出>`，可选 `--cpg` 或 `--motif`（配合 `--reference`）、`--modified-bases`、`--combine-mods`、`--combine-strands`、`--region`、`--include-bed`、`--edge-filter`、`--max-depth`、抽样选项、通过阈值选项（`--filter-percentile` 默认 10、`--filter-threshold`、`--mod-threshold`、`--no-filtering`）、`--duplex`、`--header` 与 `--bgzf`。只有高于通过阈值的判读才计为已修饰或未修饰，否则计为被过滤；18 列的 bedMethyl 给出每个位置、链和代码的覆盖度、修饰百分比以及这些计数。
- **API：** `modkit_pileup(bam_path: 'str', output_path: 'str', *, modified_bases: 'Sequence[str] | None' = None, reference: 'str | None' = None, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, combine_mods: 'bool' = False, combine_strands: 'bool' = False, mask: 'bool' = False, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, max_depth: 'int | None' = None, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, duplex: 'bool' = False, header: 'bool' = False, bgzf: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是已排序并索引的 modBAM，`output_path` 是新文件；`cpg` 或 `motif`（+`motif_offset`）需要 `reference`（带 .fai 的 FASTA）；`modified_bases` 如 `['5mC']` 或 `['C:m']`；`bgzf=True` 写 bgzip 输出，需先用 tabix 建索引才能作为 `modkit_stats`、`modkit_localize`、`modkit_bedmethyl_merge` 或 `modkit_dmr` 的输入；未设置的选项保持 modkit 默认值。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 116
['chr1', '2', '3', 'm', '1', '-'] ['1', '0.00', '0']
```

## 23) `dnakit.epigenome.modkit_pileup_hemi()`

- **作用：** 汇总回文 motif 位置上的双链修饰模式（如半甲基化）。
- **处理方法：** 运行 `modkit pileup-hemi <bam> --ref ... --out-bed ...`，使用 `--cpg`（即 motif `CG`、偏移 0 的简写）或 `--motif`，另有 `--combine-mods`、`--ignore`、`--region`、`--include-bed`、`--edge-filter`、抽样和通过阈值选项。每个分子必须在两条链上都有判读（带 `C+m?` 和 `G-m?` 的 duplex reads）；输出列出每个位置上的模式计数，如 `m,m,C`、`m,-,C` 和 `-,m,C`（第一个碱基、第二个碱基、主碱基）。
- **API：** `modkit_pileup_hemi(bam_path: 'str', reference: 'str', output_path: 'str', *, cpg: 'bool' = False, motif: 'str | None' = None, motif_offset: 'int' = 0, combine_mods: 'bool' = False, mask: 'bool' = False, ignore: 'str | None' = None, region: 'str | None' = None, include_bed: 'str | None' = None, edge_filter: 'str | None' = None, filter_percentile: 'float | None' = None, filter_threshold: 'float | None' = None, mod_thresholds: 'Sequence[str] | None' = None, no_filtering: 'bool' = False, max_depth: 'int | None' = None, num_reads: 'int | None' = None, sampling_frac: 'float | None' = None, seed: 'int | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_path` 是已排序并索引的 duplex modBAM，`reference` 是带 .fai 的 FASTA，`output_path` 是新文件；`cpg=True` 与 `motif`（+`motif_offset`，必须是反向互补回文）二选一。
- **输出：** 通用结果字典，另有 `record_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed 133
['-,-,C', '-,m,C', 'm,-,C', 'm,m,C']
```

## 24) `dnakit.epigenome.modkit_stats()`

- **作用：** 报告 BED 文件中每个区域的修饰水平（修饰比例）。
- **处理方法：** 运行 `modkit stats <bedmethyl> --regions ... --out-table ...`，可选 `--mod-codes`、`--min-coverage`（默认 1）、`--no-header`、`--threads`。把区域内所有 bedMethyl 记录的有效覆盖度和已修饰计数按修饰代码累加后相除。
- **API：** `modkit_stats(bedmethyl_path: 'str', regions_bed: 'str', output_path: 'str', *, mod_codes: 'str | None' = None, min_coverage: 'int | None' = None, no_header: 'bool' = False, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bedmethyl_path` 是带 tabix 索引的 bgzip bedMethyl，`regions_bed` 是区域，`output_path` 是新表格；`mod_codes` 是逗号分隔字符串，如 `'m,h'`。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`。
- **示例代码：**

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

- **示例输出：**

```text
completed 3
['#chrom', 'start', 'end', 'name', 'strand', 'count_m', 'count_valid_m', 'percent_m']
['chr1', '0', '120', 'reg1', '.', '49', '65', '75.38461']
['chr1', '120', '240', 'reg2', '.', '46', '57', '80.70176']
['chr2', '0', '160', 'reg3', '.', '45', '63', '71.42857']
```

## 25) `dnakit.epigenome.modkit_validate()`

- **作用：** 用真值位置验证修饰碱基判读，并报告判读与真值之间的混淆情况。
- **处理方法：** 运行 `modkit validate --bam-and-bed <bam> <bed> ... --out-filepath ...`，可选 `--canonical-base`、`--ignore`、`--edge-filter`、`--invert-edge-filter`、`--min-identity`、`--min-length`、`--filter-quantile`、`--filter-threshold`。每个真值 BED 的 name 列是修饰代码，未修饰位置用 `-`。modkit 要求 BAM 带 `NM` 标签，会构建列联表并写出机器可读的汇总。
- **API：** `modkit_validate(bam_and_bed: 'Sequence[Sequence[str]]', output_path: 'str', *, canonical_base: 'str | None' = None, ignore: 'str | None' = None, edge_filter: 'str | None' = None, invert_edge_filter: 'bool' = False, min_alignment_identity: 'float | None' = None, min_alignment_length: 'int | None' = None, filter_quantile: 'float | None' = None, filter_threshold: 'float | None' = None, threads: 'int | None' = None, modkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 modkit 程序（发布版二进制，或 `conda install -c bioconda modkit`）。DNAKit 先看 `modkit_path`，再看 PATH；找不到时会在启动程序、创建任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `bam_and_bed` 是 `[modBAM, 真值 BED]` 对的列表（每个样本一对），`output_path` 是新文件；真值只含未修饰位点或 ChEBI 代码时需要 `canonical_base`（A、C、G 或 T）。
- **输出：** 通用结果字典，另有 `columns`、`row_count` 与 `rows`（汇总的前几行）。
- **示例代码：**

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

- **示例输出：**

```text
completed 6
full_contingency_table: [["ground_truth_label","C","m","No Call"],["C",8,26,1],["m",11,22,
raw_accuracy: 45.454544
raw_contingency_table: [["ground_truth_label","C","m"],["C",8,25],["m",11,22]]
```
