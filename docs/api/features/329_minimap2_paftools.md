# minimap2 与 paftools：比对、比对格式转换、组装与注释评估

`dnakit.alignment` 中的这组函数是对 [minimap2](https://github.com/lh3/minimap2)（用 2.31 测试）及其配套脚本的带类型适配器。`minimap2_map_reads()` 运行比对器，包括剪接 RNA/cDNA 预设和读段全对全重叠预设。其余函数各运行 `paftools.js`（由 [k8](https://github.com/attractivechaos/k8) JavaScript shell 执行）或 `pafcluster.js` 的一个子命令：在 PAF、SAM、BED12、GFF3、MUMmer delta 和模拟读段格式之间转换，报告比对、组装、基因完整性和错误连接统计，转换坐标，检出和比较变异，评估比对、剪接连接点、外显子和读段重叠。每个函数以固定的、带类型的选项列表不经 shell 运行其程序。它们没有暴露全部选项，`paftools.js bedcov` 也没有封装。

minimap2 和 paftools.js 使用 MIT 许可。DNAKit 把它们作为独立进程运行，不复制其代码。

这些函数位于模块 `dnakit.alignment.minimap2_p09_cli_functions`。这些程序不是 Python 包：用 `conda install -c bioconda minimap2 k8` 安装（会把 `minimap2`、`paftools.js` 和 `k8` 放进 PATH），或用发布的压缩包（`misc/paftools.js`、`misc/pafcluster.js`），不在 PATH 中时传入 `minimap2_path`、`paftools_path`、`k8_path` 或 `pafcluster_path`。程序缺失时在创建任何输出目录之前抛出带安装提示的 `BackendUnavailableError`。每个函数都写入新的 `output_directory`（该目录不能已存在）：输出到标准输出的程序被直接写入命名的结果文件，标准错误写入日志文件，程序失败时抛出带日志末尾的 `BackendExecutionError`。输入会被校验（`ConfigurationError`）。可选参数默认为 `None`，此时不传给程序，因此沿用工具的默认值。`paftools.js call` 和 `ov-eval` 需要按目标名称和起点排序的 PAF；这两个函数会自己排序并保留排序后的副本。

## 示例数据

所有示例都从工作目录读取这几个小文件：`ref.fa` 含参考 contig `chrA`（3000 碱基，300-450、800-950、1300-1450 处有三个外显子，内含子为 GT...AG）和 `chrB`（2000 碱基）；`cdna.fa` 含剪接转录本 `tx1`（450 碱基）和带一个替换的变体；`asm.fa` 是该参考序列的一个组装，`chrA` 中有四处替换和一个 12 碱基缺失，`chrB` 中有一个 10 碱基插入。更小的输入（BED、VCF、GTF、delta、模拟器文件）由示例自己写出，每个示例都把输出写入新的临时目录。示例假设 `minimap2`、`k8` 和 `paftools.js`（最后一个示例还需要 `pafcluster.js`）在 PATH 中。

**示例文件：** `ref.fa`

```text
>chrA
CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGA
CCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCCTAGCCTAGCTACGCGTTTGCGCATCAGGCTGT
CCCATACATCAAGCGGTTCCCCTCAAATTATCCGGACTCGGTAAGGGCAGCGAGTAAATATTTTACAATA
CGTTTCTTGTCAATCTGCTGCTTTGTACGCGTCACAGTTACTCGGCGAAGGCCCGTCTTTTTGCTGACCA
GGAAATTTCACAGCTGAGCCTAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCC
TAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGA
TCAGTGCCCTCTAAGTCTCTAAGCTGTGTTGTAGGTACGAGCCCGAGCCCTTCAGGACCGAGTAAACTTG
TAGCGTTTCTCATCAGTCCAGGGGCATCCCACCCACATAACCAACCACCTATGGGTATATTCAAGTGCGG
GTGTGAAGATGCCGGTAGTCAGTATCGCATGGTCATCCACCCGACTCGTCGCGTCGGCGAACGGTCTAGG
CCAACTCTCCTTGCTACAACTATAAGACGTGTTAGGATGTGGGCGGCCAGCAGACGCAAACGCCGCCACG
TGGCTTGACGGCGTCATTCCTATTATCAAAGCAATATGTTTGCGCGACCTGGGTAGAACCTGTGCTGCGG
TTCGCCCACGTTGCGAAGACCACTTTGCAGAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACC
CTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGA
CCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGTACATTGTAGAGGGCACCAGCTCAAGTAG
TTTGATGCGCACTTCAATGACCGACCTCCGTGACCGGTGGGTGCCGAGCGTCGGCAAATCTTTCTTCCAG
TATGGTAGAGCGTCAGGTGTCTGTCCCTATCAGCTGTTGACTGCAGGGGACATGTACAACCTATCCATAA
TTACTTCCAAGCCGAGCTGTAGTTGATCACGTTAGGTAAGCTGTAGACGTACGGATTACGACCCCTTAGT
CCAGATAGGGACGCTGCTGAGCTCTTAGGCTTTCCATTGTAGAGCCAACAGATAGACGTAGTCTACTCTC
GCCCCCGGCTTAACTATAGCGTCTCAGAATACACGGGGAGGGTAGATTGAGAGGATTGTCAACAACAATT
GTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTA
CCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTGAAGGGACTCAGTCGACTATA
CAGAATCATCCTCGTGAAAAATGTCTCATCTCGTGACTTAAACGAGCAGTAACCGCCGAACGACAAGGAT
TGACGACACCGTCTCGCTTGTGCTCGAGGTTTGCCCAGCATTCCCTAAACTCCCTTCGTACAGTGGCGCT
AAGCTTAACGCTTGAGCTACGGATTACAATATCGTGAGCCCATTTGAACAGCGCGTGGCTCGATCCGCGA
TGGAGGTTTGCGGTACCCGCAATGCATCGCCAAGCGTTCCCATCCTTCGAGACCTGAGTGAAGCTTCCTG
GCGGCCTATATATCTAAGAGAGGGGAAAGGAATGCCAACCTGGACTGCCGGGAGGGCAGTTATCAGGTTA
TTGGTTGCAACATTAACCGCAGTCGGAGTCTTAACGTAATCATATTACGAAGTATCGCATAATGCGACAA
ATAATGTCGGAGTCTTTACCAGGCTTCACCTATACACGCTACAGTTATCTGAGAATCTGCCAGCTCTAGA
GTTTCCCATCGCTCCTGCCACTTAGACCGAATCTGCTTAATGTGAGGGTGCTAGTGCAGTTGGGAGACCT
GGTGTCTACCCTGTGGTAAAAGCAGTAGTCGATAAGGCACAGTAGGCCCGAGGACAACCCCCCCATACTC
CGGAGGTGGAGGCTAACGATATCATTTCGCATATAGATCATTCCTGGTAGGAGCGCGAAATCGTCGCTTC
TTGTAGGCAATTCATGCACAACTAAAGAAATTATACAGGCCCAAGCTGTAGAGTTCACGTCGCGTGGCAC
CACTCTCTAATTAACAGTATAATTTTTTTCTCTATTGAGGTGAAGCGATGTGCACTCCACTACTCCACAA
GTTCAGGCGACAATATGCATTCGTAATCCTGTTGCTCGTGGGAACCGGCGCCTTGCATTTTAGTCAGAAG
GGTTGGAAAACAAGGGGTCCTGCACATGTGACGTTAATTCTGTCTCCCACATCTCGACTCATAATGGTGC
CTGTTTTGTGGACTGTGCGAAAGGCTCTTGCTGGCATATGATTCAAACCTATCACCAGATCTCAGCTGAT
GCGACTGTGCTAGTCCATCTGCCCAAGGGCATTTCCCAATACTGATACGCTTACAGCTTGTTGAAAGTAA
AGTGCTAACATACTAGCGGCACGACACAAGGGCGGAAAAACGTTGTTGTCGTTCACTCCCGTACCTACTG
GGGCTTCTAGCCCCAACTAGGGTTGTAGCCGCATGGCCGTTCTCCCGAGTTATTAGCTGGACCATCGGAG
CATCAGACCCAATCATCAAACCAGTACGCGAATAAGACGGGATCGAAGGTTGGTGTACGATTCTGTTCCG
TACGATATACACGACGAGAAAGCTCAGCGGCAACGGTTAACCTGTGGTAGCGTCGAGGGTTTTAGGGCCC
TGTGGAATCATGTAATGTTATCAAGAACAAATGATCCACAAGTGAGAAAACCCGTGAAGTCTACTCCTCT
GTCACAGAAGGCCCCAGGATAAAAAGTCTAGTACCACGCACCAGGTTGGTGGGGGGTGAG
>chrB
CTGGACTTTTATCTGAGGATGCGGTAGCCCGCTCACAGCAATCTATTCTCTAGCGGTCGAAGTGGCTGGA
GCCGCCTTGCTTGTTAGTCGTGATGCCCTCTTTCAGAGTGTCCAAGTAAAAACCAAGAGCACCTTATTCA
TCTTACCCGTTGCGGGCCGACTTCAATTCGTGGTTCTCAGCAATATTTCAGTGGGTTATTTTTTGAGCGG
TTATGGGACCGTGCCGAAAGAACGAATTTAGCGAGTTCGCGCGACTCACAAGGACCAAGATCTGCAAATC
CGTTAAGATGAGCGGAAGGAAGGTCTTGGATGGAGCTTAAGGGTGGGAGAGGGTTGTTACTAGAAGGAGA
TGTCTCCGTGGTTCGCTTGCTAAAGGGCCAGACATGTTATCACGCGACCAGGGCTGCGGCCCAACACTGT
AAACCCTGCAAATACCATGATGGCTCCGTGAGGCTAATGACATTACAAATCAGTATATGTCGTGCTATTT
CGTGTAGTCCATTGACCACGTAGAGAAATGGATGATTAGTACTAGTTTTCTCTGTTCCGACCAGCTTGGC
CCACTTTCGATCGCCCAACAGAGAGGTGATGATATCGAGGTAGGTGGCGGATCTGATCTGTACAATTGGA
CTGTTGAATAACTAGTCAGTAAGTCTACAAAGGAACCTCGAGTCGAGTACGAGCCAACCGTGTATCATGC
AGGGCTAACAGCTACTAGCTTTCCCCGTCGACACTTTTACAATGCCTCAGACCCGGCCGGTCGAGCAGCG
GATATGAGGAGGGGATCGTTTACTGAAACTCACGGATTAAATGTTTTGATATTAGGCGTGCATCTGCCTG
GTCAGCAAAACAACTTAAACGAGAAGCTGGTGCGCAGATAACCCCTCGTTATGAATAGCGAGCACCTCCG
ACAAGTTGAATAGCCGCCTGAAACTGTAGAGAACACTAGTTACGTGGTAGGGGATTTAAGGGCACCGGTC
ATCGAACATCGTAAATGCAGGAATGTTAAATGAAGTACTCTGAGTTTCCAGCCGACATGAGAGCTCGTCC
GGCTTAATGTTTAACCATCGCGGTGGCTTTGTCGAACAGACGTCTCATCGGACCGGTGAGCACCAAAATC
TGTTTGTATACCTAATCTAGACTGGTCTGCAGCGCCACTCTGAGAATAAAGAGAAATCGTTGAGGACAAG
TCTAGAAGTGATGGAGTACTATAGTTCTGTAGATGGAGAGGTTGAAGGACTATAGGATACATTGGCGAAA
TATAGTAAGTGTATGGGAGACAAGGAAAACGGGACGACTTCACCACCAAGGCGTACCTCCCACATTTCTC
GCAACTAGTCTCGGTTCTAAAGCTTGTCCCATAGATCTACGATGGTACGGACACGCGGGGGTAACTTATT
TTTTGATCATGCCTATGGCCACCACCCGTTCTATTGGTCCATAAGATTGGGAATTCTCCGTTTGTGGGAC
ATAGGCGCTCCGCCATTTAAGGCTCCAATGCTCACTAAGATATAAGACGTCAGGACTGCGGTAATGGTGA
CGCTATAGGGGCCCTTCTAACTCGTAAAACGCTGTACCGCGCCGGTACATGCTTCTGGAAAAATCTCTTA
CGACATAAGATCAACACTCGCATCAATATGTAGTCGGACTGATGCGTCACCCATGGGCACATCGCGATAG
AATTTGATCGAAGTACGAGTCTATCTCATGTCTGTGTGGATTGAATAGATACCCCCAGTATGAAGTGGAT
CTCGGATTAGCTGAGTCGTCTGAGCAGCTAACCCAGGGAGTGGCCCCAGGTTATGGTTGTCCCCCCCGGT
ATTCCACCATACCTCTGTCATGAGTTATGGAGCTAAGTATGGTAGCCTTACATACTCATTCGCCACGACG
TGGGAAGTCGTCCGAAGCCTTGAGAACGTAACGGTCTATTCCGACACTTAGCAATTTCTCCAGTGGGCGT
CTAGTGCGTTAGTAACTAGGCCGTACGCGTGGCCCAGGGG
```

**示例文件：** `cdna.fa`

```text
>tx1
TAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCCTAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGATCAGTGCCCTCTAAGTCTCTAAGCTGTGTTAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACCCTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGACCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGGTAGATTGAGAGGATTGTCAACAACAATTGTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTACCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTG
>tx1_variant
TAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCCTAGCCAGCAGATCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGATCAGTGCCCTCTAAGTCTCTAAGCTGTGTTAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACCCTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGACCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGGTAGATTGAGAGGATTGTCAACAACAATTGTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTACCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTG
```

**示例文件：** `asm.fa`

```text
>ctgA
CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGA
CCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCCTAGCCTAGCTACGCGTTTGCGCATCAGGCTGT
CCCATACATCAAGCGGTTCCCCTCAAATTATCCGGACTCGGTAAGGGCAGCGAGTAAATATTTTACAATA
CGTTTCTTGTCAATCTGCTGCTTTGTACGCGTCACAGTTACTCGGCGAAGGCCCGTCTTTTTGCTGACCA
GGAAATTTCACAGCTGAGCCTAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCC
TAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGA
TCAGTGCCCTCTAAGTCTCTAAGCTGTGTTGTAGGTACGAGCCCGAGCCCTTCAGGACCGAGTAAACTTG
TAGCGTTTCTAATCAGTCCAGGGGCATCCCACCCACATAACCAACCACCTATGGGTATATTCAAGTGCGG
GTGTGAAGATGCCGGTAGTCAGTATCGCATGGTCATCCACCCGACTCGTCGCGTCGGCGAACGGTCTAGG
CCAACTCTCCTTGCTACAACTATAAGACGTGTTAGGATGTGGGCGGCCAGCAGACGCAAACGCCGCCACG
TGGCTTGACGGCGTCATTCCTATTATCAAAGCAATATGTTTGCGCGACCTGGGTAGAACCTGTGCTGCGG
TTCGCCCACGTTGCGAAGACCACTTTGCAGAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACC
CTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGA
CCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGTACATTGTAGAGGGCACCAGCTCAAGTAG
TTTGATGCGCACTTCAATGACCGACCTCCGTGACCGGTGGGTGCCGAGCGTCGGCAAATCTTTCTTCCAG
TATGGTAGAGCGTCAGGTGTCTGTCCCTATCAGCTGTTGACTGCAGGGGAAATGTACAACCTATCCATAA
TTACTTCCAAGCCGAGCTGTAGTTGATCACGTTAGGTAAGCTGTAGACGTACGGATTACGACCCCTTAGT
CCAGATAGGGACGCTGCTGAGCTCTTAGGCTTTCCATTGTAGAGCCAACAGATAGACGTAGTCTACTCTC
GCCCCCGGCTTAACTATAGCGTCTCAGAATACACGGGGAGGGTAGATTGAGAGGATTGTCAACAACAATT
GTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTA
CCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTGAAGGGACTCAGTCGACTATA
CAGAATCATCCTCGTGAAAAATGTCTCATCTCGTGACTTAAACGAGCAGTAACCGCCGAACGACAAGGAT
TGACGACACCGTCTCGCTTGTGCTCGAGGTTTGCCCAGCATTCCCTAAACTCCCTTCGTACAGTGGCGCT
AAGCTTAACGCTTGAGCTACGGATTACAATATCGTGAGCCCATTTGAACAGCGCGTGGCTCGATCCGCGA
TGGAGGTTTGCGGTACCCGCCATGCATCGCCAAGCGTTCCCATCCTTCGAGACCTGAGTGAAGCTTCCTG
GCGGCCTATATATCTAAGAGAGGGGAAAGGAATGCCAACCTGGACTGCCGGGAGGGCAGTTATCAGGTTA
TTGGTTGCAACATTAACCGCAGTCGGAGTCTTAACGTAATCATATTACGAAGTATCGCATAATGCGACAA
ATAATGTCGGAGTCTTTACCAGGCTTCACCTATACACGCTACAGTTATCTGAGAATCTGCCAGCTCTAGA
GTTTCCCATCGCTCCTGCCACTTAGACCGAATCTGCTTAAAGTGCAGTTGGGAGACCTGGTGTCTACCCT
GTGGTAAAAGCAGTAGTCGATAAGGCACAGTAGGCCCGAGGACAACCCCCCCATACTCCGGAGGTGGAGG
CTAACGATATCATTTCGCATATAGATCATTCCTGGTAGGAGCGCGAAATCGTCGCTTCTTGTAGGCAATT
CATGCACAACTAAAGAAATTATACAGGCCCAAGCTGTAGAGTTCACGTCGCGTGGCACCACTCTCTAATT
AACAGTATAATTTTTTTCTCTATTGAGGTGAAGCGATGTGCACTCCACTACTCCACAAGTTCAGGCGACA
ATATGCATTCGTAATCCTGTTGCTCGTGGGAACCGGCGCCTTGCATTTTAGTCAGAAGGGTTGGAAAACA
AGGGGTCCAGCACATGTGACGTTAATTCTGTCTCCCACATCTCGACTCATAATGGTGCCTGTTTTGTGGA
CTGTGCGAAAGGCTCTTGCTGGCATATGATTCAAACCTATCACCAGATCTCAGCTGATGCGACTGTGCTA
GTCCATCTGCCCAAGGGCATTTCCCAATACTGATACGCTTACAGCTTGTTGAAAGTAAAGTGCTAACATA
CTAGCGGCACGACACAAGGGCGGAAAAACGTTGTTGTCGTTCACTCCCGTACCTACTGGGGCTTCTAGCC
CCAACTAGGGTTGTAGCCGCATGGCCGTTCTCCCGAGTTATTAGCTGGACCATCGGAGCATCAGACCCAA
TCATCAAACCAGTACGCGAATAAGACGGGATCGAAGGTTGGTGTACGATTCTGTTCCGTACGATATACAC
GACGAGAAAGCTCAGCGGCAACGGTTAACCTGTGGTAGCGTCGAGGGTTTTAGGGCCCTGTGGAATCATG
TAATGTTATCAAGAACAAATGATCCACAAGTGAGAAAACCCGTGAAGTCTACTCCTCTGTCACAGAAGGC
CCCAGGATAAAAAGTCTAGTACCACGCACCAGGTTGGTGGGGGGTGAG
>ctgB
CTGGACTTTTATCTGAGGATGCGGTAGCCCGCTCACAGCAATCTATTCTCTAGCGGTCGAAGTGGCTGGA
GCCGCCTTGCTTGTTAGTCGTGATGCCCTCTTTCAGAGTGTCCAAGTAAAAACCAAGAGCACCTTATTCA
TCTTACCCGTTGCGGGCCGACTTCAATTCGTGGTTCTCAGCAATATTTCAGTGGGTTATTTTTTGAGCGG
TTATGGGACCGTGCCGAAAGAACGAATTTAGCGAGTTCGCGCGACTCACAAGGACCAAGATCTGCAAATC
CGTTAAGATGAGCGGAAGGAAGGTCTTGGATGGAGCTTAAGGGTGGGAGAGGGTTGTTACTAGAAGGAGA
TGTCTCCGTGGTTCGCTTGCTAAAGGGCCAGACATGTTATCACGCGACCAGGGCTGCGGCCCAACACTGT
AAACCCTGCAAATACCATGATGGCTCCGTGAGGCTAATGACATTACAAATCAGTATATGTCGTGCTATTT
CGTGTAGTCCATTGACCACGTAGAGAAATGGATGATTAGTACTAGTTTTCTCTGTTCCGACCAGCTTGGC
CCACTTTCGATCGCCCAACAGAGAGGTGATGATATCGAGGTAGGTGGCGGATCTGATCTGTACAATTGGA
CTGTTGAATAACTAGTCAGTAAGTCTACAAAGGAACCTCGAGTCGAGTACGAGCCAACCGTGTATCATGC
AGGGCTAACAGCTACTAGCTTTCCCCGTCGACACTTTTACAATGCCTCAGACCCGGCCGGTCGAGCAGCG
GATATGAGGAGGGGATCGTTTACTGAAACTCACGGATTAAATGTTTTGATATTAGGCGTGCATCTGCCTG
GTCAGCAAAACAACTTAAACGAGAAGCTGGTGCGCAGATAACCCCTCGTTATGAATAGCGAGCACCTCCG
ACAAGTTGAATAGCCGCCTGAAACTGTAGAGAACACTAGTTACGTGGTAGGGGATTTAAGGGCACCGGTC
ATCGAACATCGTAAATGCAGACGTACGTACGAATGTTAAATGAAGTACTCTGAGTTTCCAGCCGACATGA
GAGCTCGTCCGGCTTAATGTTTAACCATCGCGGTGGCTTTGTCGAACAGACGTCTCATCGGACCGGTGAG
CACCAAAATCTGTTTGTATACCTAATCTAGACTGGTCTGCAGCGCCACTCTGAGAATAAAGAGAAATCGT
TGAGGACAAGTCTAGAAGTGATGGAGTACTATAGTTCTGTAGATGGAGAGGTTGAAGGACTATAGGATAC
ATTGGCGAAATATAGTAAGTGTATGGGAGACAAGGAAAACGGGACGACTTCACCACCAAGGCGTACCTCC
CACATTTCTCGCAACTAGTCTCGGTTCTAAAGCTTGTCCCATAGATCTACGATGGTACGGACACGCGGGG
GTAACTTATTTTTTGATCATGCCTATGGCCACCACCCGTTCTATTGGTCCATAAGATTGGGAATTCTCCG
TTTGTGGGACATAGGCGCTCCGCCATTTAAGGCTCCAATGCTCACTAAGATATAAGACGTCAGGACTGCG
GTAATGGTGACGCTATAGGGGCCCTTCTAACTCGTAAAACGCTGTACCGCGCCGGTACATGCTTCTGGAA
AAATCTCTTACGACATAAGATCAACACTCGCATCAATATGTAGTCGGACTGATGCGTCACCCATGGGCAC
ATCGCGATAGAATTTGATCGAAGTACGAGTCTATCTCATGTCTGTGTGGATTGAATAGATACCCCCAGTA
TGAAGTGGATCTCGGATTAGCTGAGTCGTCTGAGCAGCTAACCCAGGGAGTGGCCCCAGGTTATGGTTGT
CCCCCCCGGTATTCCACCATACCTCTGTCATGAGTTATGGAGCTAAGTATGGTAGCCTTACATACTCATT
CGCCACGACGTGGGAAGTCGTCCGAAGCCTTGAGAACGTAACGGTCTATTCCGACACTTAGCAATTTCTC
CAGTGGGCGTCTAGTGCGTTAGTAACTAGGCCGTACGCGTGGCCCAGGGG
```

## 1) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_map_reads()`

- **作用：** 用 minimap2 把读段或序列比对到参考序列，包括剪接 RNA/cDNA 比对和读段间的全对全重叠检测。
- **处理方法：** 运行 `minimap2 -x <preset>`。`preset` 可为 `map-ont`、`map-hifi`、`map-pb`、`lr:hq`、`sr`、`asm5`、`asm10`、`asm20`，剪接预设 `splice`、`splice:hq`、`splice:sr`（配合 `splice_strand` = `-u`、`max_intron_length` = `-G`、`junction_bed` = `--junc-bed`），以及重叠预设 `ava-pb`、`ava-ont`（`reference` 和 `query` 给同一个文件）。`mate` 增加双端短读段的第二个文件，`output_format='sam'` 加 `-a`，`base_level_alignment` 加 `-c`，`cs_tag`、`md_tag`、`eqx` 加 `--cs`、`--MD`、`--eqx`；`kmer_size`、`window_size`、`secondary`、`max_secondary`、`soft_clip_supplementary`、`threads` 对应 `-k`、`-w`、`--secondary`、`-N`、`-Y`、`-t`。标准输出直接写入结果文件，标准错误写入 `minimap2.log`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_map_reads(reference: 'str | Path', query: 'str | Path', output_directory: 'str | Path', *, preset: "Literal['map-ont', 'map-hifi', 'map-pb', 'lr:hq', 'splice', 'splice:hq', 'splice:sr', 'sr', 'asm5', 'asm10', 'asm20', 'ava-pb', 'ava-ont']" = 'map-ont', mate: 'str | Path | None' = None, output_format: "Literal['paf', 'sam']" = 'paf', base_level_alignment: 'bool' = False, cs_tag: "Literal['short', 'long'] | None" = None, md_tag: 'bool' = False, eqx: 'bool' = False, kmer_size: 'int | None' = None, window_size: 'int | None' = None, max_intron_length: 'int | str | None' = None, splice_strand: "Literal['f', 'r', 'b', 'n'] | None" = None, junction_bed: 'str | Path | None' = None, secondary: 'bool | None' = None, max_secondary: 'int | None' = None, soft_clip_supplementary: 'bool' = False, threads: 'int | None' = None, minimap2_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 参考和查询文件（FASTA/FASTQ）以及新的输出目录。程序从 `PATH` 中查找；其他位置请传 `minimap2_path`、`paftools_path`、`k8_path`。
- **输出：** `outputs`（`alignments` 为 `alignments.paf` 或 `.sam`，`log`）、`alignment_count`、`preset`、`version`，PAF 时还有 `records`（12 个标准列加 `tags` 字典）。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads

result = minimap2_map_reads("ref.fa", "cdna.fa", tempfile.mkdtemp() + "/run", preset="splice:hq", base_level_alignment=True, verbose=False)
first = result["records"][0]
print(result["alignment_count"], first["query"], first["target"], first["target_start"], first["target_end"])
print(first["tags"]["cg"])
```

- **示例输出：**

```text
2 tx1 chrA 300 1450
150M350N150M350N150M
```

## 2) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_view_paf()`

- **作用：** 把 PAF 比对渲染为类 BLAST 文本、MAF 或 lastz-cigar。
- **处理方法：** 通过 `k8` 运行 `paftools.js view -f <output_format>`。PAF 需要长格式 `cs` 标签。`line_length` 是类 BLAST 排版的 `-l`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_view_paf(paf: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['aln', 'maf', 'lastz-cigar']" = 'aln', line_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 带 `--cs=long` 的 PAF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `view.txt`、`log`）和 `output_format`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_view_paf

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="long", verbose=False)["outputs"]["alignments"]
result = minimap2_view_paf(asm, work + "/view", line_length=60, verbose=False)
lines = open(result["outputs"]["output"]).read().splitlines()
print(lines[1][:20], lines[3][:20])
```

- **示例输出：**

```text
Ref+:          1 CTT Qry+:          1 CTT
```

## 3) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_sam_to_paf()`

- **作用：** 把 SAM 比对转换为 PAF。
- **处理方法：** 运行 `paftools.js sam2paf`。`alignments` 为 `all`、`primary_and_supplementary`（`-p`）或 `primary`（`-P`）；`include_unmapped` 为 `-U`，`long_cs` 为 `-L`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_sam_to_paf(sam: 'str | Path', output_directory: 'str | Path', *, alignments: "Literal['all', 'primary_and_supplementary', 'primary']" = 'all', include_unmapped: 'bool' = False, long_cs: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** SAM 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `converted.paf`）、`rows`（前几行）和 `row_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_sam_to_paf

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
result = minimap2_sam_to_paf(splice_sam, work + "/paf", alignments="primary", verbose=False)
print(result["row_count"], result["rows"][0][:9])
```

- **示例输出：**

```text
2 ['tx1', '450', '0', '450', '+', 'chrA', '3000', '300', '1450']
```

## 4) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_delta_to_paf()`

- **作用：** 把 MUMmer delta 比对转换为 PAF。
- **处理方法：** 运行 `paftools.js delta2paf`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_delta_to_paf(delta: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** `.delta` 文件（nucmer 输出）和新的输出目录。
- **输出：** `outputs`（`output` 为 `converted.paf`）、`rows` 和 `row_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_delta_to_paf

work = tempfile.mkdtemp()
delta = work + "/toy.delta"
open(delta, "w").write("ref.fa qry.fa\nNUCMER\n>chrX qryX 1000 900\n1 900 1 900 0 0 0\n0\n")
result = minimap2_delta_to_paf(delta, work + "/paf", verbose=False)
print(result["rows"][0][:12], result["rows"][0][12:])
```

- **示例输出：**

```text
['qryX', '900', '0', '900', '+', 'chrX', '1000', '0', '900', '900', '900', '0'] ['NM:i:0', 'cg:Z:900M']
```

## 5) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_spliced_alignments_to_bed12()`

- **作用：** 把 PAF 或 SAM 中的剪接比对转换为 BED12 基因结构。
- **处理方法：** 运行 `paftools.js splice2bed`；`keep_multiple` 为 `-m`（保留次要比对）。PAF 输入需要 `cg` 标签（`minimap2 -c`）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_spliced_alignments_to_bed12(alignments: 'str | Path', output_directory: 'str | Path', *, keep_multiple: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 剪接 PAF 或 SAM 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `structures.bed`）、`records`（BED12 的 12 列）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_spliced_alignments_to_bed12

work = tempfile.mkdtemp()
splice = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_spliced_alignments_to_bed12(splice, work + "/bed", verbose=False)
record = result["records"][0]
print(record["name"], record["start"], record["end"], record["block_count"], record["block_sizes"], record["block_starts"])
```

- **示例输出：**

```text
tx1 300 1450 3 150,150,150, 0,500,1000,
```

## 6) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_annotation_to_bed12()`

- **作用：** 把 GTF 或 GENCODE 风格的 GFF3 基因注释转换为 BED12。
- **处理方法：** 运行 `paftools.js gff2bed`。`junctions_only` 为 `-j`，`short_names` 为 `-s`，`fai` 为 `-u`，`ensembl_canonical_only` 为 `-e`，`output_gff` 为 `-g`（需 `fai`）。GFF3 输入需要 GENCODE/Ensembl 风格的 `transcript_id=` 属性；其他 GFF3 会被该工具拒绝并抛出 `BackendExecutionError`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_annotation_to_bed12(annotation: 'str | Path', output_directory: 'str | Path', *, junctions_only: 'bool' = False, short_names: 'bool' = False, fai: 'str | Path | None' = None, ensembl_canonical_only: 'bool' = False, output_gff: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** GTF 或 GFF3 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `annotation.bed` 或 `annotation.gff`）、`rows` 和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_annotation_to_bed12

work = tempfile.mkdtemp()
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_annotation_to_bed12(gtf, work + "/bed", verbose=False)
print(result["rows"][0])
junctions = minimap2_annotation_to_bed12(gtf, work + "/junctions", junctions_only=True, verbose=False)
print(junctions["rows"])
```

- **示例输出：**

```text
['chrA', '300', '1450', 't1||g1', '1000', '+', '300', '1450', '196,196,196', '3', '150,150,150,', '0,500,1000,']
[['chrA', '450', '800', 't1||g1', '1000', '+'], ['chrA', '950', '1300', 't1||g1', '1000', '+']]
```

## 7) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_gff3_to_junctions()`

- **作用：** 从 GFF3 文件中提取内含子连接点并输出为 BED。
- **处理方法：** 运行 `paftools.js gff2junc`。内含子是共享同一个 `Parent=` 的相邻 `feature` 记录（`-f`，工具默认 `CDS`）之间的间隔。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_gff3_to_junctions(gff3: 'str | Path', output_directory: 'str | Path', *, feature: 'str | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** GFF3 文件、新的输出目录和可选的 `feature`。
- **输出：** `outputs`（`output` 为 `junctions.bed`）、`records`（chrom、start、end、name、score、strand）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_gff3_to_junctions

work = tempfile.mkdtemp()
gff = work + "/genes.gff3"
rows = ["##gff-version 3", "chrA\ttest\tmRNA\t301\t1450\t.\t+\t.\tID=t1"]
rows += [f"chrA\ttest\tCDS\t{s + 1}\t{e}\t.\t+\t0\tParent=t1" for s, e in [(300, 450), (800, 950), (1300, 1450)]]
open(gff, "w").write("\n".join(rows) + "\n")
result = minimap2_gff3_to_junctions(gff, work + "/junctions", verbose=False)
print([(row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
[(450, 800), (950, 1300)]
```

## 8) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_longcs_to_sequences()`

- **作用：** 从长格式 cs 的 PAF 重建比对上的目标序列或查询序列。
- **处理方法：** 运行 `paftools.js longcs2seq`；`query=True` 为 `-q`（输出查询序列而不是目标序列）。PAF 需要 `--cs=long`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_longcs_to_sequences(paf: 'str | Path', output_directory: 'str | Path', *, query: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 长格式 cs 的 PAF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `sequences.fasta`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_longcs_to_sequences

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="long", verbose=False)["outputs"]["alignments"]
result = minimap2_longcs_to_sequences(asm, work + "/seq", query=True, verbose=False)
lines = open(result["outputs"]["output"]).read().splitlines()
print(result["record_count"], lines[0], len(lines[1]))
```

- **示例输出：**

```text
2 >ctgA_0_2988 2988
```

## 9) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_paf_to_gff3()`

- **作用：** 把 miniprot 蛋白到基因组的 PAF 转换为 GFF3 基因模型。
- **处理方法：** 运行 `paftools.js paf2gff`；`amino_acid` 为 `-a`。PAF 需要 `cg` 和 `AS` 标签。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_paf_to_gff3(paf: 'str | Path', output_directory: 'str | Path', *, amino_acid: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** miniprot 的 PAF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `models.gff3`）、`rows` 和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_paf_to_gff3

work = tempfile.mkdtemp()
paf = work + "/miniprot.paf"
open(paf, "w").write("prot1\t100\t0\t100\t+\tchrA\t3000\t300\t600\t290\t300\t60\tAS:i:100\tcg:Z:300M\n")
result = minimap2_paf_to_gff3(paf, work + "/gff", verbose=False)
for row in result["rows"]:
    print(row[2], row[3], row[4])
```

- **示例输出：**

```text
transcript 301 600
CDS 301 600
```

## 10) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_convert_simulated_reads()`

- **作用：** 把 mason2、PBSIM 或 Badread 模拟的读段转换为名称中带真实坐标的 FASTQ/FASTA。
- **处理方法：** `simulator='mason2'` 运行 `paftools.js mason2fq`（一个 SAM 文件，输出 FASTQ），`'pbsim'` 运行 `pbsim2fq`（一个或多个 MAF 文件加 `reference_fai`），`'badread'` 运行 `badread2fa`（一个 FASTQ 加 `reference_fai`）。读段名称变为 `name!chrom!start!end!strand`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_convert_simulated_reads(simulator: "Literal['mason2', 'pbsim', 'badread']", input_paths: 'Sequence[str | Path]', output_directory: 'str | Path', *, reference_fai: 'str | Path | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** `simulator`、模拟器输出文件、新的输出目录，`pbsim` 和 `badread` 还需要参考序列的 FASTA 索引。
- **输出：** `outputs`（`output` 为 `reads.fastq` 或 `reads.fasta`）、`record_count` 和 `simulator`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_convert_simulated_reads

work = tempfile.mkdtemp()
fai = work + "/ref.fa.fai"
open(fai, "w").write("chrA\t3000\t6\t70\t71\nchrB\t2000\t3100\t70\t71\n")
mason = work + "/mason.sam"
open(mason, "w").write("simulated.r1\t0\tchrB\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXE:i:0\tXS:i:0\tXI:i:0\n")
result = minimap2_convert_simulated_reads("mason2", [mason], work + "/mason", verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
maf = work + "/pbsim.maf"
open(maf, "w").write("a\ns ref 100 8 + 5000 ACGTACGT\ns S1_1 0 8 - 8 ACGT-ACGT\n")
result = minimap2_convert_simulated_reads("pbsim", [maf], work + "/pbsim", reference_fai=fai, verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
bad = work + "/badread.fq"
open(bad, "w").write("@id1 chrA,+strand,100-108 length=8 read_identity=91.5%\nACGTACGT\n+\nIIIIIIII\n")
result = minimap2_convert_simulated_reads("badread", [bad], work + "/badread", reference_fai=fai, verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
```

- **示例输出：**

```text
@r1!chrB!100!110!+ 0:0:0
>S1_1!chrA!100!108!-
>S1!chrA!100!108!+	ri:f:91.5
```

## 11) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_simulated_names_to_bed()`

- **作用：** 把模拟读段名称转换为 BED 格式的真实坐标。
- **处理方法：** 对每行一个读段名称 `name!chrom!start!end!strand` 的文件运行 `paftools.js sim2bed`（也接受 mason2 双端读段的格式）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_simulated_names_to_bed(names_file: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 读段名称文本文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `truth.bed`）、`records`（chrom、start、end、name、score、strand）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_simulated_names_to_bed

work = tempfile.mkdtemp()
names = work + "/names.txt"
open(names, "w").write("r1!chrB!100!250!+\nr2!chrA!900!1050!-\n")
result = minimap2_simulated_names_to_bed(names, work + "/bed", verbose=False)
for row in result["records"]:
    print(row["chrom"], row["start"], row["end"], row["strand"])
```

- **示例输出：**

```text
chrB 100 250 +
chrA 900 1050 -
```

## 12) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_alignment_stats()`

- **作用：** PAF 或 SAM 比对的汇总统计。
- **处理方法：** 运行 `paftools.js stat`：比对上的序列数、主要和次要比对数、碱基数、替换数和按大小分组的插入缺失数。`per_alignment_cs` 为 `-c`，`gap_output_length` 为 `-l`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_alignment_stats(alignments: 'str | Path', output_directory: 'str | Path', *, per_alignment_cs: 'bool' = False, gap_output_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** PAF（带 `cs` 标签）或 SAM 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `stats.txt`）和 `metrics`，字典，如 `substitutions`、`insertions_in_0_50`、`deletions_in_0_50`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_alignment_stats

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_alignment_stats(asm, work + "/stat", verbose=False)
metrics = result["metrics"]
print(metrics["mapped_bases"], metrics["substitutions"], metrics["insertions_in_0_50"], metrics["deletions_in_0_50"])
```

- **示例输出：**

```text
4998 4 1 1
```

## 13) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_stats()`

- **作用：** 由组装到参考序列的比对得到组装的连续性统计。
- **处理方法：** 用参考序列的 FASTA 索引和每个组装一个 PAF 运行 `paftools.js asmstat`。`min_query_length` 为 `-q`，`min_block_length` 为 `-l`（默认 10000），`max_divergence` 为 `-d`（0.01）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_stats(fai: 'str | Path', alignments: 'Sequence[str | Path]', output_directory: 'str | Path', *, min_query_length: 'int | None' = None, min_block_length: 'int | None' = None, max_divergence: 'float | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** `.fai` 文件、PAF 文件列表和新的输出目录。
- **输出：** `outputs`（`output` 为 `assembly_stats.tsv`）和 `metrics`，把每个指标（`Length`、`NG50`、`NGA50` 等）映射到每个组装一个值。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_assembly_stats

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
fai = work + "/ref.fa.fai"
open(fai, "w").write("chrA\t3000\nchrB\t2000\n")
result = minimap2_assembly_stats(fai, [asm], work + "/asmstat", min_block_length=1000, verbose=False)
print(result["metrics"]["Length"], result["metrics"]["NG50"], result["metrics"]["Rcov"])
```

- **示例输出：**

```text
[4998] [2988] ['100.00%']
```

## 14) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_gene_completeness()`

- **作用：** 由同一组转录本的剪接比对评估组装的基因完整性。
- **处理方法：** 用比对到参考序列和各个组装的转录本（`minimap2 -cx splice`）运行 `paftools.js asmgene`。`min_identity` 为 `-i`（0.99），`min_coverage` 为 `-c`（0.99），`autosomes_only` 为 `-a`，`show_incomplete` 为 `-e`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_gene_completeness(reference_paf: 'str | Path', assembly_pafs: 'Sequence[str | Path]', output_directory: 'str | Path', *, min_identity: 'float | None' = None, min_coverage: 'float | None' = None, autosomes_only: 'bool' = False, show_incomplete: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 参考序列的 PAF、组装 PAF 列表和新的输出目录。
- **输出：** `outputs`（`output` 为 `completeness.tsv`）和 `metrics`，把 `full_sgl`、`full_dup`、`frag`、`part50+` 等映射到每个输入一个计数。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_assembly_gene_completeness

work = tempfile.mkdtemp()
splice = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_assembly_gene_completeness(splice, [splice], work + "/asmgene", verbose=False)
print(result["metrics"]["full_sgl"], result["metrics"]["frag"])
```

- **示例输出：**

```text
[2, 2] [0, 0]
```

## 15) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_detect_misjoins()`

- **作用：** 在组装到参考序列的比对中检测大尺度错误连接。
- **处理方法：** 运行 `paftools.js misjoin`。`centromere_bed` 为 `-c`，`centromere_overlap_ratio` 为 `-r`，`min_block_length` 为 `-l`，`max_gap` 为 `-g`（数字或 `1m` 形式），`show_non_centromeric` 为 `-e`，`show_long_blocks` 为 `-p`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_detect_misjoins(paf: 'str | Path', output_directory: 'str | Path', *, centromere_bed: 'str | Path | None' = None, centromere_overlap_ratio: 'float | None' = None, min_block_length: 'int | str | None' = None, max_gap: 'int | str | None' = None, show_non_centromeric: 'bool' = False, show_long_blocks: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 组装到参考序列的 PAF 和新的输出目录。
- **输出：** `outputs`（`output` 为 `misjoins.txt`）和 `metrics`（染色体间错误连接、间隔和候选倒位的计数）。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_detect_misjoins

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_detect_misjoins(asm, work + "/misjoin", verbose=False)
print(result["metrics"])
```

- **示例输出：**

```text
{'inter_chromosomal_misjoins': [0, 0], 'intra_chromosomal_gaps': [0, 0], 'candidate_inversions_in_the_middle': [0, 0], 'candidate_inversions_at_contig_ends': [0, 0]}
```

## 16) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_liftover_coordinates()`

- **作用：** 通过 PAF 比对在序列之间转换 BED 区间坐标。
- **处理方法：** 运行 `paftools.js liftover`（简易 liftOver）。PAF 需要 `cg` 标签，`query_bed` 给出查询序列上的区间。`min_mapping_quality` 为 `-q`，`min_alignment_length` 为 `-l`（默认 50000），`max_divergence` 为 `-d`，`merge_intervals` 为 `-m`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_liftover_coordinates(paf: 'str | Path', query_bed: 'str | Path', output_directory: 'str | Path', *, min_mapping_quality: 'int | None' = None, min_alignment_length: 'int | None' = None, max_divergence: 'float | None' = None, merge_intervals: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** PAF、BED 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `lifted.bed`）、`rows` 和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_liftover_coordinates

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
bed = work + "/query.bed"
open(bed, "w").write("ctgA\t100\t200\nctgA\t2500\t2600\n")
result = minimap2_liftover_coordinates(asm, bed, work + "/lift", min_alignment_length=1000, verbose=False)
for row in result["rows"]:
    print(row[:3])
```

- **示例输出：**

```text
['chrA', '100', '200']
['chrA', '2512', '2612']
```

## 17) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_call_assembly_variants()`

- **作用：** 由带 cs 标签的组装到参考序列的比对检出变异。
- **处理方法：** 按工具要求对 PAF 按目标名称和起点排序（`sort -k6,6 -k8,8n`，保存为 `sorted.paf`），然后运行 `paftools.js call`。给出 `reference_fasta`（`-f`）时调用结果写成 VCF（`sample_name` 为 `-s`），否则为工具的文本格式。`min_coverage_length` 为 `-l`（10000），`min_call_length` 为 `-L`（50000），`min_mapping_quality` 为 `-q`，`gap_threshold` 为 `-g`。汇总行解析到 `summary`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_call_assembly_variants(paf: 'str | Path', output_directory: 'str | Path', *, reference_fasta: 'str | Path | None' = None, sample_name: 'str | None' = None, min_coverage_length: 'int | None' = None, min_call_length: 'int | None' = None, min_mapping_quality: 'int | None' = None, gap_threshold: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 用 `--cs` 生成的 PAF、可选的参考 FASTA 和新的输出目录。
- **输出：** `outputs`（`output` 为 `variants.vcf` 或 `variants.txt`、`sorted_paf`）和 `summary`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_call_assembly_variants

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_call_assembly_variants(asm, work + "/call", reference_fasta="ref.fa", min_coverage_length=1000, min_call_length=1000, verbose=False)
records = [line.split("\t") for line in open(result["outputs"]["output"]) if not line.startswith("#")]
print(result["summary"]["substitutions"], len(records))
print([row[:5] for row in records[:2]])
```

- **示例输出：**

```text
4 6
[['chrA', '501', '.', 'C', 'A'], ['chrA', '1101', '.', 'C', 'A']]
```

## 18) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_pair_haplotype_vcf()`

- **作用：** 把两个单倍型的 VCF 记录合并为一个定相的二倍体 VCF。
- **处理方法：** 运行 `paftools.js vcfpair`。`male` 为 `-m`，并需要 `genome_version='37'`（`-g`）来处理 chrX 的拟常染色体区；`sample_name` 为 `-s`（默认 `syndip`）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_pair_haplotype_vcf(pair_vcf: 'str | Path', output_directory: 'str | Path', *, male: 'bool' = False, genome_version: "Literal['37'] | None" = None, sample_name: 'str | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 每个单倍型一个样本列的 VCF 和新的输出目录。
- **输出：** `outputs`（`output` 为 `diploid.vcf`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_pair_haplotype_vcf

work = tempfile.mkdtemp()
pair = work + "/pair.vcf"
open(pair, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thap1\thap2\nchr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,5\t0/0:5,0\n")
result = minimap2_pair_haplotype_vcf(pair, work + "/diploid", sample_name="demo", verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[-1])
```

- **示例输出：**

```text
chr1	100	.	A	G	30	.	.	GT:AD	1|0:5,5
```

## 19) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_vcf_stats()`

- **作用：** VCF 的汇总统计：替换、转换/颠换以及按大小分组的插入缺失。
- **处理方法：** 运行 `paftools.js vcfstat`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_vcf_stats(vcf: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** VCF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `vcf_stats.txt`）和 `metrics`（`substitutions`、`ts_tv`、`insertions`、`deletions` 等）。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_vcf_stats

work = tempfile.mkdtemp()
vcf = work + "/toy.vcf"
open(vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\t.\t.\nchr1\t200\t.\tC\tA\t.\t.\t.\nchr1\t300\t.\tT\tTAAA\t.\t.\t.\nchr1\t400\t.\tGCCC\tG\t.\t.\t.\n")
result = minimap2_vcf_stats(vcf, work + "/stats", verbose=False)
metrics = result["metrics"]
print(metrics["substitutions"], metrics["ts_tv"], metrics["insertions"], metrics["deletions"])
```

- **示例输出：**

```text
2 1.0 1 1
```

## 20) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_select_vcf_by_length()`

- **作用：** 保留变异长度落在某个范围内的 VCF 记录。
- **处理方法：** 运行 `paftools.js vcfsel`；`min_length` 为 `-l`，`max_length` 为 `-L`。保留头部。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_select_vcf_by_length(vcf: 'str | Path', output_directory: 'str | Path', *, min_length: 'int | None' = None, max_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** VCF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `selected.vcf`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_select_vcf_by_length

work = tempfile.mkdtemp()
vcf = work + "/toy.vcf"
open(vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\t.\t.\nchr1\t300\t.\tT\tTAAA\t.\t.\t.\nchr1\t400\t.\tGCCC\tG\t.\t.\t.\n")
result = minimap2_select_vcf_by_length(vcf, work + "/selected", min_length=2, verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
2
```

## 21) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_sv_callsets()`

- **作用：** 比较两个 VCF 格式的结构变异调用集。
- **处理方法：** 运行 `paftools.js sveval base.vcf call.vcf`。`confident_regions_bed` 为 `-r`，`min_discard_length` 为 `-f`（30），`min_sv_length` 为 `-i`（50），`max_sv_length` 为 `-x`（100000），`window_size` 为 `-w`（500），`max_allele_difference` 为 `-d`（0.5），`show_errors` 为 `-e`，`show_matches` 为 `-p`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_sv_callsets(base_vcf: 'str | Path', call_vcf: 'str | Path', output_directory: 'str | Path', *, confident_regions_bed: 'str | Path | None' = None, min_discard_length: 'int | None' = None, min_sv_length: 'int | None' = None, max_sv_length: 'int | None' = None, window_size: 'int | None' = None, max_allele_difference: 'float | None' = None, show_errors: 'bool' = False, show_matches: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 基准 VCF、待评估 VCF 和新的输出目录。
- **输出：** `outputs`（`output` 为 `sveval.txt`）、`rows`（工具的报告行）和 `row_count`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_compare_sv_callsets

work = tempfile.mkdtemp()
header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
base, call = work + "/base.vcf", work + "/call.vcf"
open(base, "w").write(header + "chr1\t1000\t.\tA\tA" + "G" * 99 + "\t.\t.\t.\nchr1\t5000\t.\tA" + "C" * 99 + "\tA\t.\t.\t.\n")
open(call, "w").write(header + "chr1\t1010\t.\tA\tA" + "G" * 99 + "\t.\t.\t.\n")
result = minimap2_compare_sv_callsets(base, call, work + "/sveval", verbose=False)
print(result["row_count"])
print(open(result["outputs"]["output"]).read().splitlines()[0][:40])
```

- **示例输出：**

```text
3
SN	2	1	0.500000
```

## 22) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_mapping()`

- **作用：** 对照模拟读段名称中编码的真实位置评估比对准确度。
- **处理方法：** 运行 `paftools.js mapeval`。读段名称必须是 `name!chrom!start!end!strand`（见 `minimap2_convert_simulated_reads`）或 mason2 双端读段格式。`min_overlap_ratio` 为 `-r`（0.1），`print_wrong_mapq` 为 `-Q`，`mode` 为 `-m`（0 最长比对，1 第一个，2 所有主要比对）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_mapping(alignments: 'str | Path', output_directory: 'str | Path', *, min_overlap_ratio: 'float | None' = None, print_wrong_mapq: 'int | None' = None, mode: 'Literal[0, 1, 2] | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 模拟读段的 PAF 或 SAM 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `mapeval.tsv`）、`records`（每个比对质量阈值一项：`mapq`、`reads`、`wrong`、`cumulative_error_rate`、`cumulative_reads`）和 `unmapped`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_mapping

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/sim.fq"
open(reads, "w").write("".join(f"@r{i}!chrB!{s}!{s + 100}!+\n{chrB[s:s + 100]}\n+\n{'I' * 100}\n" for i, s in enumerate(range(0, 1800, 150))))
aln = minimap2_map_reads("ref.fa", reads, work + "/aln", preset="sr", base_level_alignment=True, verbose=False)["outputs"]["alignments"]
result = minimap2_evaluate_mapping(aln, work + "/mapeval", verbose=False)
print(result["records"][-1], result["unmapped"])
```

- **示例输出：**

```text
{'mapq': 60, 'reads': 12, 'wrong': 0, 'cumulative_error_rate': 0.0, 'cumulative_reads': 12} None
```

## 23) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_paf()`

- **作用：** 逐碱基比较两个 PAF 比对文件。
- **处理方法：** 运行 `paftools.js pafcmp base.paf test.paf`；`min_mapping_quality` 为 `-q`（10）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_paf(base_paf: 'str | Path', test_paf: 'str | Path', output_directory: 'str | Path', *, min_mapping_quality: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 两个 PAF 文件和新的输出目录。
- **输出：** `outputs`（`output` 为 `pafcmp.txt`）和 `metrics`，统计测试文件中正确、错误、缺失或多出的碱基比对。只比较长度至少 5000 碱基且只有一个主要比对的查询序列。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_compare_paf

work = tempfile.mkdtemp()
base, test = work + "/base.paf", work + "/test.paf"


def line(name, start):
    return f"{name}\t6000\t0\t6000\t+\tchr1\t1000000\t{start}\t{start + 6000}\t5900\t6000\t60\n"


open(base, "w").write(line("q1", 100) + line("q2", 20000) + line("q3", 40000))
open(test, "w").write(line("q1", 150) + line("q3", 90000))
result = minimap2_compare_paf(base, test, work + "/cmp", verbose=False)
print(result["metrics"])
```

- **示例输出：**

```text
{'base_alignments_with_mapq_10': 3, 'base_alignments_correctly_mapped_by_test': 1, 'wrong_test_alignment': 1, 'base_alignments_missing': 1, 'additional_test_alignments_with_mapq_10': 0}
```

## 24) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_junctions()`

- **作用：** 对照基因注释评估剪接比对的剪接连接点。
- **处理方法：** 运行 `paftools.js junceval annotation.gtf alignments`。`input_format` 为 `sam`（默认）、`miniprot_paf`（`-a`）或 `bed`（`-b`）；`tolerance` 为 `-l`，`standard_chromosomes_only` 为 `-c`，`first_alignment_only` 为 `-1`，`print_overlapping` 为 `-p`，`print_errors` 为 `-e`。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_junctions(annotation_gtf: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['sam', 'miniprot_paf', 'bed']" = 'sam', tolerance: 'int | None' = None, standard_chromosomes_only: 'bool' = False, first_alignment_only: 'bool' = False, print_overlapping: 'bool' = False, print_errors: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** GTF 文件、剪接比对和新的输出目录。
- **输出：** `outputs`（`output` 为 `junceval.txt`）、`metrics`（`predicted_introns`、`correct_introns`、`correct_introns_percent` 等）、`report_text` 和 `rows`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_junctions

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_evaluate_junctions(gtf, splice_sam, work + "/junc", verbose=False)
print(result["metrics"])
```

- **示例输出：**

```text
{'unmapped_reads': 0, 'mapped_reads': 2, 'primary_alignments': 2, 'singletons': 0, 'predicted_introns': 4, 'non_overlapping_introns': 0, 'correct_introns': 4, 'correct_introns_percent': 100.0}
```

## 25) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_exons()`

- **作用：** 对照基因注释评估剪接比对的外显子。
- **处理方法：** 运行 `paftools.js exoneval annotation.gtf alignments`。选项与剪接连接点评估相同，另有 `coding_only`（`-d`）、`skip_first_exon`（`-f`）、`skip_terminal_exons`（`-t`）和 `base_level`（`-s`，碱基水平的灵敏度和特异度）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_exons(annotation_gtf: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['sam', 'miniprot_paf', 'bed']" = 'sam', tolerance: 'int | None' = None, coding_only: 'bool' = False, standard_chromosomes_only: 'bool' = False, first_alignment_only: 'bool' = False, skip_first_exon: 'bool' = False, skip_terminal_exons: 'bool' = False, base_level: 'bool' = False, print_overlapping: 'bool' = False, print_errors: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** GTF 文件、剪接比对和新的输出目录。
- **输出：** `outputs`（`output` 为 `exoneval.txt`）、`metrics`（`predicted_exons`、`correct_exons` 等）、`report_text` 和 `rows`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_exons

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_evaluate_exons(gtf, splice_sam, work + "/exon", verbose=False)
print(result["metrics"])
```

- **示例输出：**

```text
{'unmapped_reads': 0, 'mapped_reads': 2, 'primary_alignments': 2, 'predicted_exons': 6, 'non_overlapping_exons': 0, 'correct_exons': 6, 'correct_exons_percent': 100.0}
```

## 26) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_overlaps()`

- **作用：** 对照读段到参考序列的比对评估读段重叠检测的灵敏度。
- **处理方法：** 按目标名称和起点对读段到参考序列的 PAF 排序（`sort -k6,6 -k8,8n`，保存为 `sorted_to_ref.paf`），然后运行 `paftools.js ov-eval`。`min_overlap_length` 为 `-l`（2000），`min_mapping_quality` 为 `-q`（10），`min_mapped_fraction` 为 `-f`（0.95）。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_overlaps(read_to_reference_paf: 'str | Path', overlap_paf: 'str | Path', output_directory: 'str | Path', *, min_overlap_length: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_mapped_fraction: 'float | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 读段比对到参考序列的 PAF、全对全重叠 PAF（`ava-ont`/`ava-pb`）和新的输出目录。
- **输出：** `outputs`（`output` 为 `ov_eval.txt`、`sorted_to_ref_paf`）、`metrics`（`inferred_overlaps`、`missed_overlaps`、`sensitivity_percent`）和 `report_text`。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_overlaps

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/long.fa"
open(reads, "w").write("".join(f">l{i}\n{chrB[s:s + 800]}\n" for i, s in enumerate(range(0, 1000, 150))))
to_ref = minimap2_map_reads("ref.fa", reads, work + "/ref", preset="map-ont", verbose=False)["outputs"]["alignments"]
overlaps = minimap2_map_reads(reads, reads, work + "/ava", preset="ava-ont", verbose=False)["outputs"]["alignments"]
result = minimap2_evaluate_overlaps(to_ref, overlaps, work + "/ov", min_overlap_length=200, verbose=False)
print(result["report_text"])
```

- **示例输出：**

```text
15 overlaps inferred from the reference mapping
0 missed by the read overlapper
100.00% sensitivity
```

## 27) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_cluster_overlap_reads()`

- **作用：** 把全对全重叠 PAF 中的读段聚成读段组。
- **处理方法：** 通过 `k8` 运行 `pafcluster.js`（在 minimap2 源码树的 `misc/` 目录中）。`min_coverage` 为 `-c`（0.9），`max_divergence` 为 `-d`（0.015），`max_difference` 为 `-e`（20000）。脚本从 `PATH` 或 `pafcluster_path` 中查找。
- **API：** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_cluster_overlap_reads(overlap_paf: 'str | Path', output_directory: 'str | Path', *, min_coverage: 'float | None' = None, max_divergence: 'float | None' = None, max_difference: 'float | None' = None, pafcluster_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 全对全重叠 PAF 和新的输出目录。
- **输出：** `outputs`（`output` 为 `clusters.txt`）和工具记录格式的 `rows`（`SD` 种子、`CL` 簇长度、以 `//` 结束的块）。
- **示例代码：**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_cluster_overlap_reads

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/long.fa"
open(reads, "w").write("".join(f">l{i}\n{chrB[s:s + 1000]}\n" for i, s in enumerate(range(0, 800, 200))))
overlaps = minimap2_map_reads(reads, reads, work + "/ava", preset="ava-ont", verbose=False)["outputs"]["alignments"]
result = minimap2_cluster_overlap_reads(overlaps, work + "/clusters", verbose=False)
print(result["row_count"], result["rows"][:3])
```

- **示例输出：**

```text
12 [['SD', 'l0', '1'], ['CL', 'l0', '1000'], ['//']]
```
