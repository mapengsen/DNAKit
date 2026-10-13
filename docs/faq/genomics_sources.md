# DNA 项目来源与覆盖范围

本页区分原始方法、上游软件、Skills、MCP 服务与 DNAKit 新增适配。收录项目不代表已安装依赖或验证其科学结果。

Snapshot: **2026-09-20**. Projects/sources: **419**. Logged searches: **88**.

针对 `new_vision/2026-10-10/public_dna_tools.md` 的逐项清单，另有[公开 DNA 工具覆盖矩阵](../api/features/40_public_tool_coverage.md)：原文声明 282 条表格记录，当前 Markdown 可解析出 425 次明确提及，保守去重后得到 320 条来源记录；分类、引用和不同提供方的记录会保留为证据，拼写别名才会合并。`dnakit.integrations` 提供离线查询、引用和统一路由，并为 scikit-bio、seqme、DNA Chisel、Biotite、Edlib、GeneBlocks、sourmash、PyRanges、pyfastx、dnaio、mappy、pyBigWig、pybedtools、gffutils、pySBOL3、scadnano、polars-bio、DnaCauldron、DNA Features Viewer、gtars、pyfaidx、vcfpy 与 cyvcf2 提供可选 pip 适配器；Cutadapt、fastp、FastQC、MultiQC、BWA、Bowtie2、SAMtools、BCFtools、SPAdes、SPOA、NanoFilt、Filtlong、NanoPlot、Porechop、Trimmomatic、AfterQC、SeqMagick、GATK4、GEMINI、Ensembl VEP、Shasta、Merqury、BUSCO、DeepVariant、cuteSV、Sniffles2、DELLY、BEDOPS、Flye、hifiasm、MEGAHIT、vcfanno、SURVIVOR 和 Picard 通过固定参数外部适配器接入；不把“发现了项目”误标为“已复制全部功能”。

检索不限于最初 11 个领域；这是可持续补充的版本化清单，不能据此宣称已经穷尽全球全部 DNA 论文和项目。

## 1) BOLE

Wang Zishuai、Liang Chongxiao、Zhang Yanlin、Zhou Rong、Zhang Xiaoai、Wei Wenkang、Li Kui 的 BOLE 描述了育种智能体框架。补充清单提供了公开基准与证据仓库，但尚未核实完整可运行框架或公开 API；因此保留来源，不声称复现 BOLE。本地关系矩阵与岭回归预测依据对应数学方法独立实现。

[BOLE (2026)](https://doi.org/10.1016/j.jia.2026.07.011) · [VanRaden (2008)](https://doi.org/10.3168/jds.2007-0980) · [Endelman (2011)](https://doi.org/10.3835/plantgenome2011.08.0024)

[BOLE benchmark evidence](https://github.com/ChongxiaoLiang/BOLE) · [335-project review / 补充清单审核](user_catalog.md)

## 2) 方向覆盖

| Domain / 方向 | Project records / 项目记录数 |
| --- | ---: |
| 序列处理与格式 (`sequence`) | 57 |
| 测序、质控与覆盖度 (`sequencing`) | 21 |
| 比对与参考序列 (`alignment`) | 10 |
| 基因组组装 (`assembly`) | 15 |
| 泛基因组与图基因组 (`pangenome`) | 10 |
| 基因识别与功能注释 (`annotation`) | 36 |
| 启动子、增强子与调控 (`regulation`) | 26 |
| 小变异与结构变异 (`variants`) | 41 |
| 疾病、肿瘤与药物基因组 (`clinical`) | 42 |
| 育种、性状与种质资源 (`breeding`) | 21 |
| 关联分析、遗传力与因果遗传 (`gwas`) | 19 |
| 群体遗传与祖源 (`population`) | 39 |
| 比较基因组与系统发育 (`comparative`) | 32 |
| 遗传与进化模拟 (`simulation`) | 4 |
| 多倍体与单倍型 (`polyploid`) | 9 |
| 法医、亲缘与 DNA 鉴定 (`forensics`) | 13 |
| 物种、食品与环境溯源 (`traceability`) | 8 |
| 古 DNA 与损伤 (`ancient-dna`) | 6 |
| 环境 DNA 与保护遗传 (`edna`) | 6 |
| 微生物、宏基因组与耐药 (`metagenomics`) | 42 |
| 甲基化与表观遗传 (`epigenomics`) | 39 |
| 染色质空间结构与接触图 (`3d-genomics`) | 7 |
| 单细胞 DNA 与多组学 (`single-cell-dna`) | 6 |
| 空间基因组 (`spatial-genomics`) | 1 |
| 基因编辑与脱靶 (`crispr`) | 27 |
| 合成生物学与序列优化 (`synthetic`) | 18 |
| 复制、损伤修复与端粒 (`replication-repair`) | 5 |
| 功能筛选与变异效应实验 (`functional-assays`) | 18 |
| 重复与移动遗传元件 (`repeats`) | 11 |
| 病毒与噬菌体基因组 (`viral-genomics`) | 4 |
| 线粒体与叶绿体基因组 (`organellar`) | 2 |
| DNA 纳米结构与分子力学 (`nanotechnology`) | 4 |
| DNA 数据存储与编码 (`dna-storage`) | 8 |
| DNA 链置换与分子计算 (`strand-displacement`) | 5 |
| DNA 基础模型与表示 (`models`) | 21 |
| 数据库与知识检索 (`databases`) | 15 |
| 智能体、Skills、MCP 与复现 (`workflows`) | 48 |

## 3) 每项功能引用了谁？

以下列出项目或论文作者，不将上游工作归为 DNAKit 原创；未核实个人作者时标明仓库作者与贡献者。原有工作流引用仍保存在各 Skill 的 UPSTREAM.json 与 SKILL.md。

| Project / 项目 | Authors and sources / 作者与来源 | DNAKit integration / 本地接口或状态 |
| --- | --- | --- |
| [ABC enhancer-gene model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) | https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction authors and contributors; individual paper authors are not inferred  | Related local workflows: causal-genomics-effector-gene-prioritization |
| [ABRicate](https://github.com/tseemann/abricate) | https://github.com/tseemann/abricate authors and contributors; individual paper authors are not inferred  | abricate-screen |
| [ADMIXTOOLS 2](https://github.com/uqrmaie1/admixtools) | uqrmaie1 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: population-genetics-population-structure |
| [ADMIXTURE](https://dalexander.github.io/admixture/) | ADMIXTURE authors and contributors; individual paper authors are not inferred | 来源记录；尚未集成 |
| [Aldy](https://github.com/0xTCG/aldy) | https://github.com/0xTCG/aldy authors and contributors; individual paper authors are not inferred | 来源记录；尚未集成 |
| [AlphaGenome API](https://github.com/google-deepmind/alphagenome) | https://github.com/google-deepmind/alphagenome authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41586-025-10014-0) | 来源记录；尚未集成 |
| [AlphaGenome Research](https://github.com/google-deepmind/alphagenome_research) | https://github.com/google-deepmind/alphagenome_research authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41586-025-10014-0) | dnakit.predictions.api.predict_properties |
| [AlphaMate](https://github.com/AlphaGenes/AlphaMate) | AlphaGenes 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: ecological-genomics-conservation-genetics |
| [AlphaSimR](https://github.com/gaynorr/AlphaSimR) | R. Chris Gaynor, Gregor Gorjanc, John M. Hickey [paper 1](https://doi.org/10.1093/g3journal/jkaa017) [supplied source 1](https://doi.org/10.1101/gr.083634.108) | Related local workflows: plant-genomics |
| [AMRFinderPlus](https://github.com/ncbi/amr) | ncbi项目作者与贡献者；terms of the United States Copyright Act.  It was written as part of；thus cannot be copyrighted.  This software/database is freely available [paper 1](https://doi.org/10.1038/s41598-021-91456-0) [paper 2](https://pubmed.ncbi.nlm.nih.gov/34135355/) | amrfinder-nucleotide; Related local workflows: metagenomics-amr-detection |
| [ANGSD](https://github.com/ANGSD/angsd) | https://github.com/ANGSD/angsd authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [AnnotSV](https://github.com/lgmgeo/AnnotSV) | https://github.com/lgmgeo/AnnotSV authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Anthropic Life Sciences](https://github.com/anthropics/life-sciences) | https://github.com/anthropics/life-sciences authors and contributors; individual paper authors are not inferred  | MCP: pubmed-anthropic; MCP: opentargets-platform; MCP: 10x-genomics; 3 bundled Skills |
| [arcasHLA](https://github.com/RabadanLab/arcasHLA) | https://github.com/RabadanLab/arcasHLA authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ArchR](https://github.com/GreenleafLab/ArchR) | https://github.com/GreenleafLab/ArchR authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ARIBA](https://github.com/sanger-pathogens/ariba) | https://github.com/sanger-pathogens/ariba authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ASCAT](https://github.com/VanLoo-lab/ascat) | https://github.com/VanLoo-lab/ascat authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ASTRAL](https://github.com/smirarab/ASTRAL) | https://github.com/smirarab/ASTRAL authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) | https://github.com/Gaius-Augustus/Augustus authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/bioinformatics/btn013) | augustus-genes |
| [AutoBA](https://github.com/Panomica/AutoBA) | https://github.com/Panomica/AutoBA authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1002/advs.202407094) | 来源记录；尚未集成 |
| [AutoPVS1](https://github.com/JiguangPeng/autopvs1) | https://github.com/JiguangPeng/autopvs1 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [AWS Kiro Life Sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences) | https://github.com/aws-samples/sample-kiro-power-life-sciences authors and contributors; individual paper authors are not inferred  | MCP: kiro-genomics; MCP: kiro-model-organisms; MCP: kiro-epigenomics; MCP: kiro-agriculture; MCP: kiro-ecology; MCP: kiro-clinical; MCP: kiro-pipelines; 1 bundled Skills |
| [Azimuth](https://github.com/MicrosoftResearch/Azimuth) | https://github.com/MicrosoftResearch/Azimuth authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BAGEL2](https://github.com/hart-lab/bagel) | https://github.com/hart-lab/bagel authors and contributors; individual paper authors are not inferred  | Related local workflows: crispr-screens-bagel-essentiality, crispr-screens-screen-qc |
| [Bakta](https://github.com/oschwengers/bakta) | https://github.com/oschwengers/bakta authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1099/mgen.0.000685) | `dnakit.annotation.annotate_bakta`（需用户提供 Bakta、数据库） |
| [Barrnap](https://github.com/tseemann/barrnap) | tseemann/barrnap authors and contributors  | barrnap-rrna |
| [Basenji](https://github.com/calico/basenji) | https://github.com/calico/basenji authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BCFtools](https://github.com/samtools/bcftools) | Danecek et al. (2021), [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008); MIT | `dnakit.ngs.call_variants`（需用户提供 bcftools） |
| [BE-Hive bystander model](https://github.com/maxwshen/be_predict_bystander) | https://github.com/maxwshen/be_predict_bystander authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BEDTools](https://github.com/arq5x/bedtools2) | Quinlan & Hall (2010), [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033); MIT | `dnakit.alignment.intersect_bed`（需用户提供 bedtools） |
| [BGLR](https://github.com/gdlc/BGLR-R) | Paulino Pérez, Gustavo de los Campos [paper 1](https://doi.org/10.1534/genetics.114.164442) [paper 2](https://doi.org/10.1093/genetics/iyac112) | Related local workflows: plant-genomics |
| [BioCoder](https://github.com/gersteinlab/BioCoder) | https://github.com/gersteinlab/BioCoder authors and contributors; individual paper authors are not inferred  [supplied source 1](https://arxiv.org/abs/2308.16458) | 来源记录；尚未集成 |
| [BioContextAI Knowledgebase MCP](https://github.com/biocontext-ai/knowledgebase-mcp) | https://github.com/biocontext-ai/knowledgebase-mcp authors and contributors; individual paper authors are not inferred  | MCP: biocontext-kb |
| [BioContextAI Registry](https://github.com/biocontext-ai/registry) | https://github.com/biocontext-ai/registry authors and contributors; individual paper authors are not inferred  [supplied source 1](https://www.nature.com/articles/s41587-025-02900-9) | 来源记录；尚未集成 |
| [bioframe](https://github.com/open2c/bioframe) | https://github.com/open2c/bioframe authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BioinfoMCP](https://github.com/florensiawidjaja/BioinfoMCP) | florensiawidjaja/BioinfoMCP authors and contributors [paper 1](https://arxiv.org/abs/2510.02139) | 来源记录；尚未集成 |
| [Bioinformatics MCP](https://github.com/raonyguimaraes/bioinformatics-mcp) | raonyguimaraes/bioinformatics-mcp authors and contributors  | 来源记录；尚未集成 |
| [BioMCP](https://github.com/genomoncology/biomcp) | genomoncology项目作者与贡献者；Copyright (c) 2025 Ian Maurer；The above copyright notice and this permission notice shall be included in all  | MCP: biomcp; 1 bundled Skills; Related local workflows: clinical-databases-clinvar-lookup, variant-interpretation |
| [BioMCP TypeScript](https://github.com/yeyuan98/biomcp-ts) | yeyuan98项目作者与贡献者；"Licensor" shall mean the copyright owner or entity authorized by；the copyright owner that is granting the License. [paper 1](https://doi.org/10.5281/zenodo.22178957) | Related local workflows: clinical-databases-clinvar-lookup |
| [Biopython](https://github.com/biopython/biopython) | https://github.com/biopython/biopython authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/bioinformatics/btp163) | dnakit.bioinformatics.reads.read_fastx; dnakit.bioinformatics.alignment.read_alignment; dnakit.bioinformatics.motifs.build_motif |
| [bioSkills](https://github.com/GPTomics/bioSkills) | GPTomics项目作者与贡献者；Copyright (c) 2026 GPTomics；The above copyright notice and this permission notice shall be included in all  | 331 bundled Skills; Related local workflows: metagenomics-kraken, methylation-calling, variant-calling-deepvariant |
| [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | Biostrings authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BioThings MCP](https://github.com/longevity-genie/biothings-mcp) | https://github.com/longevity-genie/biothings-mcp authors and contributors; individual paper authors are not inferred  | MCP: biothings |
| [BioTools MCP](https://github.com/neuralx-dev/biotools-mcp-server) | neuralx-dev 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: blast-searches, phylo-tree-io |
| [bioutils](https://github.com/biocommons/bioutils) | https://github.com/biocommons/bioutils authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Bismark](https://github.com/FelixKrueger/Bismark) | FelixKrueger项目作者与贡献者  | Related local workflows: methylation-bismark-alignment, methylation-calling |
| [BOLD Data Portal API](https://www.boldsystems.org/data/api/) | BOLD Data Portal API authors and contributors [paper 1](https://arxiv.org/abs/2404.05696) | Related local workflows: ecological-genomics-edna-metabarcoding |
| [BOLE](https://bole.zishuailab.com/) | Zishuai Wang, Chongxiao Liang, Yanlin Zhang, Rong Zhou, Xiaoai Zhang, Wenkang Wei, Kui Li [paper 1](https://doi.org/10.1016/j.jia.2026.07.011) [paper 2](https://www.chinaagrisci.com/Jwk_zgnykxen/EN/10.1016/j.jia.2026.07.011) [paper 3](https://www.biorxiv.org/content/10.64898/2026.02.04.703673v1) [supplied source 1](https://www.sciencedirect.com/science/article/pii/S2095311926003011) | Related local workflows: causal-genomics-heritability-partitioning, plant-genomics, workflows-gwas-pipeline |
| [Borzoi](https://github.com/calico/borzoi) | https://github.com/calico/borzoi authors and contributors; individual paper authors are not inferred  [supplied source 1](https://www.biorxiv.org/content/10.1101/2023.08.30.555582v1) | 来源记录；尚未集成 |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) | Langmead & Salzberg (2012), [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923); GPL-3.0 | `dnakit.ngs.map_reads`（需用户提供 bowtie2 和 bowtie2-build） |
| [BPNet](https://github.com/kundajelab/bpnet) | https://github.com/kundajelab/bpnet authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Bracken](https://github.com/jenniferlu717/Bracken) | jenniferlu717项目作者与贡献者；GPL-3.0；[paper 1](https://doi.org/10.7717/peerj-cs.104) [paper 2](https://doi.org/10.1038/s41596-022-00738-y) | `dnakit.microbiome.estimate_bracken_abundance`（需用户提供 Bracken、Kraken2 报告和数据库）; Related local workflows: metagenomics-abundance |
| [BRAKER](https://github.com/Gaius-Augustus/BRAKER) | https://github.com/Gaius-Augustus/BRAKER authors and contributors [paper 1](https://doi.org/10.1093/nar/gki937) [paper 2](https://doi.org/10.1101/550145) [paper 3](https://doi.org/10.1093/nargab/lqaa026) [paper 4](https://doi.org/10.48550/arXiv.2403.19416) | Related local workflows: genome-annotation-eukaryotic-gene-prediction |
| [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) | https://github.com/Gaius-Augustus/BRAKER4 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BrAPI MCP](https://github.com/cyanheads/brapi-mcp-server) | https://github.com/cyanheads/brapi-mcp-server authors and contributors; individual paper authors are not inferred  | MCP: brapi |
| [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) | BSgenome authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [bsseq](https://github.com/hansenlab/bsseq) | https://github.com/hansenlab/bsseq authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [BUSCO](https://busco.ezlab.org/) | BUSCO authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1093/molbev/msab199) | busco-proteins; busco-genome |
| [BV-BRC MCP](https://github.com/BV-BRC/BV-BRC-MCP) | BV-BRC项目作者与贡献者  | Related local workflows: microbial-genome-characterization |
| [BWA](https://github.com/lh3/bwa) | Li & Durbin (2009), [10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324); GPL-3.0 | `dnakit.ngs.map_reads`（需用户提供 bwa） |
| [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) | Vasimuddin et al. (2019), [10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041); MIT | `dnakit.ngs.map_reads`（`method="bwa-mem2"`，需用户提供 bwa-mem2） |
| [IQ-TREE 2](https://github.com/iqtree/iqtree2) | Minh et al. (2020), [10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015); GPL-2.0 | `dnakit.phylogeny.infer_ml_tree`（需用户提供 iqtree2） |
| [FastTree](https://github.com/morgannprice/fasttree) | Price, Dehal & Arkin (2010), [10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490); GPL-3.0 | `dnakit.phylogeny.infer_fast_tree`（需用户提供 FastTree） |
| [FreeBayes](https://github.com/freebayes/freebayes) | Garrison & Marth (2012), [arXiv:1207.3907](https://arxiv.org/abs/1207.3907); MIT | `dnakit.variants.call_freebayes`（需用户提供 freebayes） |
| [GEMINI](https://github.com/arq5x/gemini) | Paila et al. (2013), [10.1371/journal.pcbi.1003153](https://doi.org/10.1371/journal.pcbi.1003153)；MIT | `dnakit.variants.load_gemini_database`、`dnakit.variants.query_gemini_database`（需用户提供 GEMINI 与注释资源） |
| [BUSCO](https://gitlab.com/ezlab/busco) | Manni et al. (2021), [10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199)；BUSCO 软件 MIT，lineage 数据集另有条款 | `dnakit.assembly.run_busco`（需用户提供 BUSCO 与 lineage 数据） |
| [DeepVariant](https://github.com/google/deepvariant) | Poplin et al. (2018), [10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235)；BSD-3-Clause | `dnakit.variants.call_deepvariant`（需用户提供 DeepVariant 与模型） |
| [BEDOPS](https://github.com/bedops/bedops) | Neph et al. (2012), [10.1093/bioinformatics/bts277](https://doi.org/10.1093/bioinformatics/bts277)；GPL-2.0 | `dnakit.alignment.sort_bed_bedops`（需用户提供 sort-bed） |
| [Flye](https://github.com/mikolmogorov/Flye) | Kolmogorov et al. (2019), [10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)；BSD-3-Clause | `dnakit.assembly.assemble_flye`（需用户提供 Flye） |
| [hifiasm](https://github.com/chhylp123/hifiasm) | Cheng et al. (2021), [10.1038/s41592-020-01056-5](https://doi.org/10.1038/s41592-020-01056-5)；MIT | `dnakit.assembly.assemble_hifiasm`（需用户提供 hifiasm） |
| [MEGAHIT](https://github.com/voutcn/megahit) | Li et al. (2015), [10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)；GPL-3.0 | `dnakit.assembly.assemble_megahit`（需用户提供 MEGAHIT） |
| [Merqury](https://github.com/marbl/merqury) | Rhie et al. (2020), [10.1186/s13059-020-02134-9](https://doi.org/10.1186/s13059-020-02134-9)；Public domain / United States Government Work | `dnakit.assembly.evaluate_merqury`（需用户提供 Merqury 与 Meryl 数据库） |
| [Shasta](https://github.com/paoloshasta/shasta) | Shafin et al. (2020), [10.1038/s41587-020-0503-6](https://doi.org/10.1038/s41587-020-0503-6)；当前源码 MIT | `dnakit.assembly.assemble_shasta`（需用户提供 Shasta） |
| [GATK4](https://github.com/broadinstitute/gatk) | McKenna et al. (2010), [10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110)；当前源码 `LICENSE.TXT` 为 Apache-2.0 | `dnakit.variants.call_gatk_haplotypecaller`（需用户提供 GATK4） |
| [VCFtools](https://github.com/vcftools/vcftools) | Danecek et al. (2011), [10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330); LGPL-3.0 | `dnakit.variants.filter_vcf_vcftools`（需用户提供 vcftools） |
| [vcfanno](https://github.com/brentp/vcfanno) | Pedersen et al. (2016), [10.1186/s13059-016-0973-5](https://doi.org/10.1186/s13059-016-0973-5); MIT | `dnakit.variants.annotate_vcf_vcfanno`（需用户提供 vcfanno） |
| [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) | Jeffares et al. (2017), [10.1038/NCOMMS14061](https://doi.org/10.1038/NCOMMS14061); MIT | `dnakit.variants.merge_structural_variants_survivor`（需用户提供 SURVIVOR） |
| [Picard](https://github.com/broadinstitute/picard) | Broad Institute, *Picard Toolkit* (2019); MIT | `dnakit.ngs.mark_duplicates_picard`（需用户提供 Picard） |
| [SPOA](https://github.com/rvaser/spoa) | Christopher Lee (2003), [10.1093/bioinformatics/btg109](https://doi.org/10.1093/bioinformatics/btg109); MIT | `dnakit.alignment.consensus_spoa`（需用户提供 spoa） |
| [NanoFilt](https://github.com/wdecoster/nanofilt) | De Coster et al. (2018), [10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149); GPL-3.0 | `dnakit.sequencing_qc.filter_nanopore_reads`（需用户提供 NanoFilt；上游已停止更新） |
| [Filtlong](https://github.com/rrwick/Filtlong) | Filtlong 官方项目未提供 DOI；GPL-3.0 | `dnakit.sequencing_qc.filter_long_reads`（需用户提供 filtlong） |
| [NanoPlot](https://github.com/wdecoster/NanoPlot) | De Coster & Rademakers (2023), [10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311); MIT | `dnakit.sequencing_qc.summarize_long_reads`（需用户提供 NanoPlot） |
| [Porechop](https://github.com/rrwick/Porechop) | Wick et al. (2017), [10.1099/mgen.0.000132](https://doi.org/10.1099/mgen.0.000132); GPL-3.0 | `dnakit.sequencing_qc.trim_porechop`（需用户提供 Porechop；官方已停止维护） |
| [Trimmomatic](https://github.com/usadellab/Trimmomatic) | Bolger, Lohse & Usadel (2014), [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170); GPL-3.0；adapter 序列另有 Illumina 条款 | `dnakit.sequencing_qc.trim_trimmomatic`（需用户提供 Trimmomatic） |
| [AfterQC](https://github.com/OpenGene/AfterQC) | Chen et al. (2017), [10.1186/s12859-017-1469-3](https://doi.org/10.1186/s12859-017-1469-3)；MIT | `dnakit.sequencing_qc.run_afterqc`（需用户提供 after.py；官方 README 建议维护中的工作流使用 fastp） |
| [SeqMagick](https://github.com/fhcrc/seqmagick) | Matsen Group seqmagick 项目；官方项目未提供 DOI；GPL-3.0 | `dnakit.fasta_fastq.convert_seqmagick`（需用户提供 seqmagick） |
| [bwa-meth](https://github.com/brentp/bwa-meth) | https://github.com/brentp/bwa-meth authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CADD](https://github.com/kircherlab/CADD-scripts) | https://github.com/kircherlab/CADD-scripts authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/nar/gkad989) [supplied source 2](https://doi.org/10.1186/s13073-021-00835-9) [supplied source 3](https://doi.org/10.1038/ng.2892) | 来源记录；尚未集成 |
| [Caduceus](https://github.com/kuleshov-group/caduceus) | https://github.com/kuleshov-group/caduceus authors and contributors; individual paper authors are not inferred  | Related local workflows: caduceus |
| [CalicoST](https://github.com/raphael-group/CalicoST) | raphael-group项目作者与贡献者；Copyright (c) 2023, Princeton University；1. Redistributions of source code must retain the above copyright notice, this  | Related local workflows: copy-number-allele-specific-copy-number |
| [Canu](https://github.com/marbl/canu) | https://github.com/marbl/canu authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Cas-OFFinder](https://github.com/snugel/cas-offinder) | https://github.com/snugel/cas-offinder authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btu048) | cas-offinder; Related local workflows: crispr-screens-library-design |
| [Cas-OFFinder bulge branch](https://github.com/pnucolab/cas-offinder-2) | https://github.com/pnucolab/cas-offinder-2 authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btu048) | Related local workflows: crispr-screens-library-design |
| [hubayirp/agentic-science](https://github.com/hubayirp/agentic-science) | hubayirp/agentic-science contributors  | 3 bundled Skills |
| [NVIDIA-BioNeMo/bionemo-agent-toolkit](https://github.com/nvidia-bionemo/bionemo-agent-toolkit) | NVIDIA Corporation and contributors  | 3 bundled Skills |
| [ClawBio/ClawBio](https://github.com/clawbio/clawbio) | ClawBio contributors  | 23 bundled Skills; Related local workflows: amplicon-qc, analyze-fasta, ancestry-risk-profiler, article-data-fetcher, bioconductor-bridge, crispr-screen-triage, eqtl-catalogue-region-fetch, equity-scorer, galaxy-bridge, gwas-catalog-region-fetch, illumina-bridge, ld-1000g-region-compute, locuscompare-region-render, multiqc-reporter, nfcore-sarek-wrapper, rare-high-impact-variants, recombinator, sample-qc-triage, seq-wrangler, ukb-ppp-region-fetch, wgs-prs |
| [DamLabResources/damlab-skills](https://github.com/damlabresources/damlab-skills) | DamLabResources/damlab-skills contributors  | 2 bundled Skills |
| [google-deepmind/science-skills](https://github.com/google-deepmind/science-skills) | Copyright 2026 Google LLC  | 9 bundled Skills |
| [companion-inc/feynman](https://github.com/companion-inc/feynman) | companion-inc/feynman contributors  | 1 bundled Skills |
| [muwenbo/genetic-counseling-skills](https://github.com/muwenbo/genetic-counseling-skills) | muwenbo/genetic-counseling-skills contributors  | 6 bundled Skills |
| [juanesarango/isabl-mcp](https://github.com/juanesarango/isabl-mcp) | juanesarango/isabl-mcp contributors  | 8 bundled Skills |
| [tanaylab/misha](https://github.com/tanaylab/misha) | tanaylab/misha contributors  | 1 bundled Skills |
| [KalarisLabs/Myriad](https://github.com/kalarislabs/myriad) | MYRIAD by Srihari Muralikrishnan and Kalaris Labs. https://github.com/KalarisLabs/myriad  | 1 bundled Skills |
| [htlin222/ngs-tertiary-analysis-skills](https://github.com/htlin222/ngs-tertiary-analysis-skills) | htlin222/ngs-tertiary-analysis-skills contributors  | 1 bundled Skills |
| [Li-Chia-ching/Plant-qPCR-Agent-Skill](https://github.com/li-chia-ching/plant-qpcr-agent-skill) | Li-Chia-ching/Plant-qPCR-Agent-Skill contributors  | 1 bundled Skills |
| [Xuzhen-Li/bioinfo-agent-skills](https://github.com/xuzhen-li/bioinfo-agent-skills) | Xuzhen-Li/bioinfo-agent-skills contributors  | 6 bundled Skills |
| [Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline](https://github.com/li-chia-ching/polyploid-genefamily-ai-pipeline) | Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline contributors  | 1 bundled Skills |
| [Ynzuses/research-skills-suite](https://github.com/ynzuses/research-skills-suite) | Ynzuses/research-skills-suite contributors  | 2 bundled Skills |
| [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/sciagent-skills) | jaechang-hits; credited K-Dense-AI and Biomni contributors  | 21 bundled Skills |
| [mims-harvard/ToolUniverse](https://github.com/mims-harvard/tooluniverse) | ToolUniverse contributors, mims-harvard  | MCP: tooluniverse; 35 bundled Skills |
| [ccsmeth](https://github.com/PengNi/ccsmeth) | https://github.com/PengNi/ccsmeth authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Cello](https://github.com/CIDARLAB/cello) | https://github.com/CIDARLAB/cello authors and contributors  | 来源记录；尚未集成 |
| [Cello v2](https://github.com/CIDARLAB/Cello-v2) | https://github.com/CIDARLAB/Cello-v2 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) | https://github.com/DaehwanKimLab/centrifuge authors and contributors; GPL-3.0；[paper 1](https://doi.org/10.1101/gr.210641.116) | `dnakit.microbiome.classify_centrifuge`（需用户提供 Centrifuge 索引） |
| [Chamaeleo](https://github.com/ntpz870817/Chamaeleo) | ntpz870817/Chamaeleo authors and contributors [paper 1](https://doi.org/10.1038/nature11875) | dnakit.genomics.encode_dna_storage; dnakit.genomics.decode_dna_storage; Related local workflows: dnakit-dna-storage |
| [ChAMP](https://bioconductor.org/packages/release/bioc/html/ChAMP.html) | ChAMP authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CheckM2](https://github.com/chklovski/CheckM2) | https://github.com/chklovski/CheckM2 authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1038/s41592-023-01940-w) | `dnakit.annotation.predict_checkm2`（需用户提供 CheckM2、数据库）; Related local workflows: checkm2 |
| [chewBBACA](https://github.com/B-UMMI/chewBBACA) | https://github.com/B-UMMI/chewBBACA authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CHISEL](https://github.com/raphael-group/chisel) | raphael-group项目作者与贡献者；Copyright (c) 2019, Princeton University；1. Redistributions of source code must retain the above copyright notice, this [paper 1](https://doi.org/10.1038/s41587-020-0661-6) | Related local workflows: copy-number-allele-specific-copy-number |
| [Chorus](https://github.com/pinellolab/chorus) | https://github.com/pinellolab/chorus authors and contributors  | MCP: chorus |
| [ChromBPNet](https://github.com/kundajelab/chrombpnet) | https://github.com/kundajelab/chrombpnet authors and contributors [paper 1](https://www.biorxiv.org/content/10.1101/2024.12.25.630221v1) [paper 2](https://www.biorxiv.org/content/early/2024/12/25/2024.12.25.630221) [paper 3](https://www.biorxiv.org/content/early/2024/12/25/2024.12.25.630221.full.pdf) [supplied source 1](https://doi.org/10.1101/2024.12.25.630221) | Related local workflows: atac-seq-footprinting |
| [chromVAR](https://bioconductor.org/packages/release/bioc/html/chromVAR.html) | chromVAR authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Clair Skills](https://github.com/HKU-BAL/Clair-skills) | HKU-BAL项目作者与贡献者  | Related local workflows: long-read-sequencing-clair3-variants |
| [Clair3](https://github.com/HKU-BAL/Clair3) | HKU-BAL项目作者与贡献者；Copyright 2021 The University of Hong Kong, Department of Computer Science；1. Redistributions of source code must retain the above copyright notice, this [paper 1](https://www.biorxiv.org/content/10.64898/2026.02.13.705285v1) [paper 2](https://www.biorxiv.org/content/10.1101/2021.12.29.474431v2) [paper 3](https://doi.org/10.1093/bioinformatics/btag181) [paper 4](https://doi.org/10.1093/bioinformatics/btag376) | Related local workflows: long-read-sequencing-clair3-variants |
| [ClassifyCNV](https://github.com/Genotek/ClassifyCNV) | https://github.com/Genotek/ClassifyCNV authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CNVkit](https://github.com/etal/cnvkit) | etal项目作者与贡献者；Copyright (c) 2013-2018 Eric Talevich；Copyright (c) 2013-2018 University of California  | Related local workflows: copy-number-cnvkit-analysis |
| [coloc](https://github.com/chr1swallace/coloc) | https://github.com/chr1swallace/coloc authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1371/journal.pgen.1011697) [supplied source 2](https://doi.org/10.1371/journal.pgen.1009440) | 来源记录；尚未集成 |
| [Control-FREEC](https://github.com/BoevaLab/FREEC) | https://github.com/BoevaLab/FREEC authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Cooler](https://github.com/open2c/cooler) | open2c项目作者与贡献者；Copyright (c) 2015-2023, Cooler developers；1. Redistributions of source code must retain the above copyright notice, this [paper 1](https://doi.org/10.1093/bioinformatics/btz540) [paper 2](https://doi.org/10.1093/bioinformatics/btz540}) | Related local workflows: hi-c-analysis-hic-data-io, hi-c-analysis-matrix-operations |
| [Cooltools](https://github.com/open2c/cooltools) | open2c项目作者与贡献者；Copyright (c) 2017 Cooltools developers；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.1101/2022.10.31.514564) | Related local workflows: hi-c-analysis-compartment-analysis, hi-c-analysis-tad-detection |
| [CRABS](https://github.com/gjeunen/reference_database_creator) | G.-J. Jeunen, E. Dowle, J. Edgecombe, U. von Ammon, N. J. Gemmell, H. Cross [paper 1](https://doi.org/10.1111/1755-0998.13741) | Related local workflows: ecological-genomics-edna-metabarcoding |
| [CREsted](https://github.com/aertslab/CREsted) | https://github.com/aertslab/CREsted authors and contributors [paper 1](https://www.nature.com/articles/s41592-026-03057-2) [paper 2](https://doi.org/10.5281/zenodo.15045960) | 来源记录；尚未集成 |
| [CRISPOR](https://github.com/maximilianh/crisporWebsite) | https://github.com/maximilianh/crisporWebsite authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CRISPR-GPT public light](https://github.com/cong-lab/crispr-gpt-pub) | https://github.com/cong-lab/crispr-gpt-pub authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41551-025-01463-z) | 来源记录；尚未集成 |
| [crisprDesign](https://github.com/crisprVerse/crisprDesign) | https://github.com/crisprVerse/crisprDesign authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [CRISPResso2](https://github.com/pinellolab/CRISPResso2) | https://github.com/pinellolab/CRISPResso2 authors and contributors [paper 1](https://doi.org/10.1038/s41587-019-0032-3) | crispresso2-amplicon; Related local workflows: crispr-screens-crispresso-editing |
| [CRISPRitz](https://github.com/pinellolab/CRISPRitz) | https://github.com/pinellolab/CRISPRitz authors and contributors  | Related local workflows: crispr-screens-library-design |
| [CRISPRme](https://github.com/pinellolab/CRISPRme) | https://github.com/pinellolab/CRISPRme authors and contributors [paper 1](https://doi.org/10.1038/s41588-022-01257-y) | 来源记录；尚未集成 |
| [crisprScore](https://github.com/crisprVerse/crisprScore) | https://github.com/crisprVerse/crisprScore authors and contributors [paper 1](https://www.biorxiv.org/content/10.1101/2022.04.21.488824v3) [paper 2](https://doi.org/10.1101/2022.06.27.497780) | Related local workflows: crispr-screens-library-design |
| [Cutadapt](https://github.com/marcelm/cutadapt) | Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200); MIT | `dnakit.sequencing_qc.run_cutadapt`（需用户提供 Cutadapt 可执行程序） |
| [cuteSV](https://github.com/tjiangHIT/cuteSV) | Jiang et al. (2020), [10.1186/s13059-020-02107-y](https://doi.org/10.1186/s13059-020-02107-y)；MIT | `dnakit.variants.call_cutesv`（需用户提供 cuteSV） |
| [cyvcf2](https://github.com/brentp/cyvcf2) | https://github.com/brentp/cyvcf2 authors and contributors; individual paper authors are not inferred  | `dnakit.integrations.summarize_vcf_cyvcf2`（需 `dnakit[bio-cyvcf2]`） |
| [DADA2](https://github.com/benjjneb/dada2) | https://github.com/benjjneb/dada2 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DAS Tool](https://github.com/cmks/DAS_Tool) | https://github.com/cmks/DAS_Tool authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41564-018-0171-1) | 来源记录；尚未集成 |
| [deepDNAshape](https://github.com/JinsenLi/deepDNAshape) | JinsenLi/deepDNAshape authors and contributors  | 来源记录；尚未集成 |
| [DeepSignal](https://github.com/bioinfomaticsCSU/deepsignal) | https://github.com/bioinfomaticsCSU/deepsignal authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DeepSTARR](https://github.com/bernardo-de-almeida/DeepSTARR) | https://github.com/bernardo-de-almeida/DeepSTARR authors and contributors [paper 1](https://www.nature.com/articles/s41588-022-01048-5) [paper 2](https://doi.org/10.1002/cpmb.105) [paper 3](https://doi.org/10.1038/s41586-019-1210-7) [paper 4](http://dx.doi.org/10.1038/nature13994) [paper 5](https://doi.org/10.5281/zenodo.5502060) | 来源记录；尚未集成 |
| [deepTools](https://github.com/deeptools/deepTools) | https://github.com/deeptools/deepTools authors and contributors; individual paper authors are not inferred  | Related local workflows: cite-encode |
| [DeepVariant](https://github.com/google/deepvariant) | google项目作者与贡献者；Copyright 2020 Google LLC.；1. Redistributions of source code must retain the above copyright notice, this [paper 1](https://doi.org/10.1038/nbt.4235) [paper 2](https://doi.org/10.1093/bioinformatics/btaa1081) | Related local workflows: variant-calling-deepvariant |
| [DELLY](https://github.com/dellytools/delly) | Rausch et al. (2012), [10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)；BSD-3-Clause | `dnakit.variants.call_delly`（需用户提供 DELLY） |
| [DIAMOND](https://github.com/bbuchfink/diamond) | https://github.com/bbuchfink/diamond authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1038/s41592-021-01101-x) | `dnakit.alignment.search_diamond_blastp`; `dnakit.alignment.search_diamond_blastx`（需用户提供 DIAMOND） |
| [DMRcate](https://bioconductor.org/packages/release/bioc/html/DMRcate.html) | DMRcate authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [dmrseq](https://bioconductor.org/packages/release/bioc/html/dmrseq.html) | dmrseq authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DNA-Aeon](https://github.com/MW55/DNA-Aeon) | MW55/DNA-Aeon authors and contributors [paper 1](https://doi.org/10.1038/s41467-023-36297-3) | 来源记录；尚未集成 |
| [DNA Cauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron) | https://github.com/Edinburgh-Genome-Foundry/DnaCauldron authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) | https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer authors and contributors; individual paper authors are not inferred  | `dnakit.integrations.summarize_dna_features_viewer` |
| [DNA Fountain](https://github.com/TeamErlich/dna-fountain) | https://github.com/TeamErlich/dna-fountain authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DNABERT2](https://github.com/MAGICS-LAB/DNABERT_2) | https://github.com/MAGICS-LAB/DNABERT_2 authors and contributors; individual paper authors are not inferred  | Related local workflows: dnabert2 |
| [DNA Chisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) | https://github.com/Edinburgh-Genome-Foundry/DnaChisel authors and contributors  [supplied source 1](https://academic.oup.com/bioinformatics/article/36/16/4508/5869515) | `dnakit.integrations.optimize_with_dnachisel`（需 `dnakit[bio-dnachisel]`） |
| [DNAplotlib](https://github.com/VoigtLab/dnaplotlib) | https://github.com/VoigtLab/dnaplotlib authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DNAscent](https://github.com/MBoemo/DNAscent) | https://github.com/MBoemo/DNAscent authors and contributors [paper 1](https://doi.org/10.1038/s41467-025-63168-w) [paper 2](https://doi.org/10.1093/nar/gkad093) [paper 3](https://doi.org/10.1186/s12864-021-07736-6) [paper 4](https://www.nature.com/articles/s41592-019-0394-y) | 来源记录；尚未集成 |
| [DNAshapeR](https://github.com/TsuPeiChiu/DNAshapeR) | TsuPeiChiu/DNAshapeR authors and contributors  | 来源记录；尚未集成 |
| [DNASpiderWeb](https://github.com/HaolingZHANG/DNASpiderWeb) | HaolingZHANG/DNASpiderWeb authors and contributors  | 来源记录；尚未集成 |
| [dnastore](https://github.com/ihh/dnastore) | ihh/dnastore authors and contributors [paper 1](https://arxiv.org/abs/1606.01799) | 来源记录；尚未集成 |
| [DNA Weaver](https://github.com/Edinburgh-Genome-Foundry/DnaWeaver) | https://github.com/Edinburgh-Genome-Foundry/DnaWeaver authors and contributors  | 来源记录；尚未集成 |
| [Dorado](https://github.com/nanoporetech/dorado) | https://github.com/nanoporetech/dorado authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [DR4DNA](https://github.com/umr-ds/DR4DNA) | umr-ds/DR4DNA authors and contributors  | 来源记录；尚未集成 |
| [DrugZ](https://github.com/hart-lab/drugz) | https://github.com/hart-lab/drugz authors and contributors; individual paper authors are not inferred  | drugz-unpaired; Related local workflows: crispr-screens-drugz-chemogenomic |
| [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html) | DSS authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [EBI Tools MCP (QuentinCody)](https://github.com/QuentinCody/ebi-tools-mcp-server) | https://github.com/QuentinCody/ebi-tools-mcp-server authors and contributors; individual paper authors are not inferred  | MCP: ebi-tools |
| [EDTA](https://github.com/oushujun/EDTA) | https://github.com/oushujun/EDTA authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | https://github.com/eggnogdb/eggnog-mapper authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1093/molbev/msab293) [paper 2](https://doi.org/10.1093/nar/gkaf1249) | eggnog-proteins |
| [EIGENSOFT](https://github.com/DReichLab/EIG) | https://github.com/DReichLab/EIG authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ENCODE Toolkit](https://github.com/ammawla/encode-toolkit) | ammawla项目作者与贡献者；Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>；with two steps: (1) assert copyright on the software, and (2) offer [paper 1](https://doi.org/10.5281/zenodo.18917511) [paper 2](https://doi.org/10.5281/zenodo.18917511}) | MCP: encode-toolkit; 16 bundled Skills; Related local workflows: accessibility-aggregation, hic-aggregation, methylation-aggregation |
| [Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) | https://github.com/google-deepmind/deepmind-research/tree/master/enformer authors and contributors; individual paper authors are not inferred  | dnakit.predictions.api.predict_sequence_properties; dnakit.predictions.api.predict_enformer_benchmark |
| [Enrich2](https://github.com/FowlerLab/Enrich2) | https://github.com/FowlerLab/Enrich2 authors and contributors [paper 1](https://doi.org/10.5281/zenodo.14681278) [paper 2](https://doi.org/10.1186/s13059-017-1272-5) | 来源记录；尚未集成 |
| [Ensembl MCP server](https://github.com/cyanheads/ensembl-mcp-server) | https://github.com/cyanheads/ensembl-mcp-server authors and contributors  | MCP: ensembl |
| [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) | McLaren et al. (2016), [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4); Apache-2.0 | `dnakit.variants.annotate_vcf_ensembl_vep`（需用户提供 VEP、缓存、参考序列和插件） |
| [EpiDISH](https://bioconductor.org/packages/release/bioc/html/EpiDISH.html) | EpiDISH authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [EuroForMix](https://github.com/oyvble/euroformix) | oyvble 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: variant-calling |
| [Evo2](https://github.com/ArcInstitute/evo2) | https://github.com/ArcInstitute/evo2 authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41586-026-10176-5) | dnakit.predictions.api.predict_properties; Related local workflows: evo2 |
| [Exomiser](https://github.com/exomiser/Exomiser) | exomiser项目作者与贡献者  | Related local workflows: clinical-databases-variant-prioritization, rare-disease-diagnosis |
| [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | https://github.com/Illumina/ExpansionHunter authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [f5c](https://github.com/hasindu2008/f5c) | https://github.com/hasindu2008/f5c authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [FACETS](https://github.com/mskcc/facets) | https://github.com/mskcc/facets authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [FastANI](https://github.com/ParBLiSS/FastANI) | ParBLiSS/FastANI authors and contributors [paper 1](https://doi.org/10.1038/s41467-018-07641-9) | fastani |
| [fastp](https://github.com/OpenGene/fastp) | Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); MIT | `dnakit.sequencing_qc.run_fastp`（需用户提供 fastp 可执行程序） |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | FastQC authors and contributors; individual paper authors are not inferred; GPL-3.0 | `dnakit.sequencing_qc.run_fastqc`（需用户提供 FastQC 可执行程序） |
| [fastSTRUCTURE](https://github.com/rajanil/fastStructure) | https://github.com/rajanil/fastStructure authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [FDSTools](https://github.com/Jerrythafast/FDSTools) | Jerrythafast 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: variant-calling |
| [FlashFry](https://github.com/mckennalab/FlashFry) | https://github.com/mckennalab/FlashFry authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Flye](https://github.com/fenderglass/Flye) | fenderglass项目作者与贡献者；Copyright (c) 2016, The Regents of the University of California；* Redistributions of source code must retain the above copyright [paper 1](https://doi.org/10.5281/zenodo.5228989) [paper 2](https://doi.org/10.1038/s41592-020-00971-x) [paper 3](https://doi.org/10.1038/s41587-019-0072-8) [paper 4](https://www.doi.org/10.1073/pnas.1604560113) | flye-assembly; Related local workflows: genome-assembly-long-read-assembly, genome-assembly-metagenome-assembly |
| [FORECasT / SelfTarget](https://github.com/felicityallen/SelfTarget) | https://github.com/felicityallen/SelfTarget authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [forrel / pedsuite](https://github.com/magnusdv/forrel) | Magnus Dehli Vigeland / pedsuite contributors; cite Pedigree Analysis in R (2021), ISBN 9780128244302, as requested by pedsuite.  | Related local workflows: genomics-study-design |
| [FreeBayes](https://github.com/freebayes/freebayes) | https://github.com/freebayes/freebayes authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Galaxy](https://github.com/galaxyproject/galaxy) | https://github.com/galaxyproject/galaxy authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Galaxy Claude Skills (nekrut)](https://github.com/nekrut/galaxy-claude-skills) | https://github.com/nekrut/galaxy-claude-skills authors and contributors; individual paper authors are not inferred  | 1 bundled Skills |
| [Galaxy MCP](https://github.com/galaxyproject/galaxy-mcp) | galaxyproject/galaxy-mcp authors and contributors  | 来源记录；尚未集成 |
| [GAPIT](https://github.com/jiabowang/GAPIT) | Jiabo Wang, Zhiwu Zhang [paper 1](https://doi.org/10.1093/molbev/msag107) [supplied source 1](https://doi.org/10.1093/bioinformatics/bts444) [supplied source 2](https://doi.org/10.1016/j.gpb.2021.08.005) [supplied source 3](https://doi.org/10.1093/gigascience/giy154) [supplied source 4](https://doi.org/10.1371/journal.pgen.1005767) | Related local workflows: workflows-gwas-pipeline |
| [GATK](https://github.com/broadinstitute/gatk) | broadinstitute项目作者与贡献者  [supplied source 1](https://genome.cshlp.org/content/20/9/1297.short) | Related local workflows: copy-number-gatk-cnv, gatk-variant-calling |
| [GATK-SV](https://github.com/broadinstitute/gatk-sv) | https://github.com/broadinstitute/gatk-sv authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [GCTA](https://github.com/jianyangqt/gcta) | jianyangqt/gcta authors and contributors  | gcta-grm; gcta-reml |
| [GEMMA](https://github.com/genetics-statistics/GEMMA) | https://github.com/genetics-statistics/GEMMA authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1038/ng.2310) | gemma-lmm |
| [GeneChat MCP](https://github.com/natecostello/genechat-mcp) | natecostello项目作者与贡献者；Copyright (c) 2026 GeneChat Contributors；The above copyright notice and this permission notice shall be included in all  | MCP: genechat; Related local workflows: clinical-databases-pharmacogenomics, vcf-basics |
| [GeneGPT](https://github.com/ncbi/GeneGPT) | https://github.com/ncbi/GeneGPT authors and contributors; individual paper authors are not inferred  [supplied source 1](https://academic.oup.com/bioinformatics/article/40/2/btae075/7606338) | 来源记录；尚未集成 |
| [GenePattern MCP](https://github.com/genepattern/genepattern-mcp) | genepattern/genepattern-mcp authors and contributors  | MCP: genepattern |
| [geNomad](https://github.com/apcamargo/genomad) | https://github.com/apcamargo/genomad authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1038/s41587-023-01953-y) | genomad-mobile-elements; Related local workflows: genomad |
| [Genomic Intelligence skill](https://github.com/K-Dense-AI/scientific-agent-skills) | https://github.com/K-Dense-AI/scientific-agent-skills authors and contributors  | 9 bundled Skills; Related local workflows: genomic-intelligence |
| [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) | GenomicRanges authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [gffutils](https://github.com/daler/gffutils) | gffutils documentation and daler contributors; MIT；来源调查未提供 DOI | `dnakit.integrations.query_gffutils_database`（需 `dnakit[bio-gffutils]`） |
| [GeneBlocks](https://github.com/Edinburgh-Genome-Foundry/Geneblocks) | [official documentation](https://edinburgh-genome-foundry.github.io/Geneblocks/); Edinburgh Genome Foundry contributors; MIT；来源调查未提供 DOI | `dnakit.integrations.diff_geneblocks_sequences`（需 `dnakit[bio-geneblocks]`） |
| [sourmash](https://github.com/sourmash-bio/sourmash) | [JOSS paper](https://doi.org/10.21105/joss.06830); BSD-3-Clause | `dnakit.integrations.compare_sourmash_sequences`（需 `dnakit[bio-sourmash]`） |
| [pySBOL3](https://github.com/SynBioDex/pySBOL3) | Mitchell, Beal & Bartley (2022), [paper](https://doi.org/10.1021/acssynbio.2c00164); MIT | `dnakit.integrations.summarize_sbol3_document`（需 `dnakit[bio-pysbol3]`） |
| [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano) | Doty, Lee & Stérin (2020), [arXiv:2005.11841](https://arxiv.org/abs/2005.11841); MIT | `dnakit.integrations.summarize_scadnano_design`（需 `dnakit[bio-scadnano]`） |
| [polars-bio](https://github.com/biodatageeks/polars-bio) | polars-bio, [paper](https://doi.org/10.1093/bioinformatics/btaf640); Apache-2.0 | `dnakit.integrations.overlap_polars_bio_bed`（需 `dnakit[bio-polarsbio]`） |
| [DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron) | DnaCauldron documentation and Edinburgh Genome Foundry contributors; survey supplied no DOI; MIT | `dnakit.integrations.simulate_dnacauldron_assembly`（需 `dnakit[bio-dnacauldron]`） |
| [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) | DNA Features Viewer documentation and Edinburgh Genome Foundry contributors; survey supplied no DOI; MIT | `dnakit.integrations.summarize_dna_features_viewer`（需 `dnakit[bio-dna-features-viewer]`） |
| [dnaio](https://github.com/marcelm/dnaio) | dnaio 软件归档，[Zenodo DOI](https://doi.org/10.5281/zenodo.10548864)；MIT；这是软件 DOI，不是同行评审论文 | `dnakit.integrations.summarize_dnaio_fastx`（需 `dnakit[bio-dnaio]`） |
| [gtars](https://github.com/databio/gtars) | [official documentation](https://docs.bedbase.org/gtars/python-overview/) and databio contributors; survey supplied no DOI; BSD-2-Clause | `dnakit.integrations.summarize_gtars_regions`（需 `dnakit[bio-gtars]`） |
| [gget](https://github.com/scverse/gget) | https://github.com/scverse/gget authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [gget MCP](https://github.com/longevity-genie/gget-mcp) | https://github.com/longevity-genie/gget-mcp authors and contributors; individual paper authors are not inferred  | MCP: gget |
| [Gigwa-MCP](https://github.com/gkanogiannis/Gigwa-MCP) | gkanogiannis / Gigwa-MCP contributors; separate attribution to Gigwa and underlying methods.  | MCP: gigwa-mcp; 5 bundled Skills; Related local workflows: gigwa-diversity-report, gigwa-explore-instance, gigwa-import-and-qc, gigwa-qc-triage, gigwa-region-scan, plant-genomics, population-genetics-scikit-allel-analysis |
| [GrEBI](https://github.com/EBISPOT/GrEBI) | https://github.com/EBISPOT/GrEBI authors and contributors; individual paper authors are not inferred  | MCP: grebi |
| [gReLU](https://github.com/Genentech/gReLU) | https://github.com/Genentech/gReLU authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | Ecogenomics项目作者与贡献者；Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>；(1) assert copyright on the software, and (2) offer you this License [paper 1](https://www.nature.com/articles/s41564-021-00918-8) [paper 2](https://doi.org/10.1038/s41587-020-0501-8) [paper 3](https://www.nature.com/articles/nbt.4229) [paper 4](http://dx.doi.org/10.1038/nbt.4229) [paper 5](https://www.nature.com/articles/s41592-023-02018-3) | `dnakit.microbiome.classify_gtdbtk`（需用户提供 GTDB-Tk 与参考数据） |
| [GuideScan CLI](https://github.com/pritykinlab/guidescan-cli) | https://github.com/pritykinlab/guidescan-cli authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [GWAS Catalog MCP](https://github.com/gwas-catalog/mcp) | gwas-catalog 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: gwas-database |
| [GWAS MCP (muslus)](https://github.com/muslus/gwas-mcp) | https://github.com/muslus/gwas-mcp authors and contributors; individual paper authors are not inferred  | MCP: gwas-mcp |
| [GWASpoly](https://github.com/jendelman/GWASpoly) | jendelman 项目及贡献者；论文作者见逐项原始引用 [paper 1](https://doi.org/10.3835/plantgenome2015.08.0073) | Related local workflows: workflows-gwas-pipeline |
| [Hail](https://hail.is/) | Hail authors and contributors; individual paper authors are not inferred  | 来源记录；已接入可选适配器（`bio-hail`） |
| [hap.py](https://github.com/Illumina/hap.py) | https://github.com/Illumina/hap.py authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [haploqtl](https://github.com/tayloranthonyanderson/haploqtl) | Taylor A. Anderson, Martha A. Sudermann, Darlene M. DeJong, David M. Francis, Christine D. Smart, Martha A. Mutschler [paper 1](https://doi.org/10.1111/tpj.16495) | 1 bundled Skills; Related local workflows: plant-genomics, population-genetics-population-structure |
| [Hartwig HMFtools](https://github.com/hartwigmedical/hmftools) | https://github.com/hartwigmedical/hmftools authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Helixer](https://github.com/usadellab/Helixer) | https://github.com/usadellab/Helixer authors and contributors [paper 1](https://doi.org/10.1038/s41592-025-02939-1) [paper 2](https://doi.org/10.1093/bioinformatics/btaa1044) | Related local workflows: genome-annotation-eukaryotic-gene-prediction |
| [hgvs](https://github.com/biocommons/hgvs) | https://github.com/biocommons/hgvs authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [HiCExplorer](https://github.com/deeptools/HiCExplorer) | https://github.com/deeptools/HiCExplorer authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/gigascience/giac061) | 来源记录；尚未集成 |
| [hifiasm](https://github.com/chhylp123/hifiasm) | chhylp123项目作者与贡献者；Copyright (c) 2019；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.1038/s41592-020-01056-5) [paper 2](https://doi.org/10.1038/s41587-022-01261-x) [paper 3](https://doi.org/10.1038/s41592-024-02269-8) | hifiasm-assembly; Related local workflows: genome-assembly-hifi-assembly, optimize-hifi-diploid-assembly |
| [HipSTR](https://github.com/HipSTR-Tool/HipSTR) | Thomas Willems (upstream README author); cite Genome-wide profiling of heritable and de novo STR variations, Nature Methods, doi:10.1038/nmeth.4267. [paper 1](https://www.nature.com/articles/nmeth.4267) | Related local workflows: variant-calling |
| [HTSlib](https://github.com/samtools/htslib) | Bonfield et al. (2021), [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007)；MIT/BSD | `dnakit.ngs.compress_index_vcf_htslib`（需用户提供 bgzip/tabix） |
| [HUMAnN](https://github.com/biobakery/humann) | https://github.com/biobakery/humann authors and contributors; individual paper authors are not inferred | `dnakit.microbiome.profile_humann`（需用户提供 HUMAnN 数据库） |
| [HyenaDNA](https://github.com/HazyResearch/hyena-dna) | https://github.com/HazyResearch/hyena-dna authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ichorCNA](https://github.com/broadinstitute/ichorCNA) | https://github.com/broadinstitute/ichorCNA authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [IGV](https://github.com/igvteam/igv) | https://github.com/igvteam/igv authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [inDelphi](https://github.com/maxwshen/inDelphi-model) | https://github.com/maxwshen/inDelphi-model authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [inStrain](https://github.com/MrOlm/inStrain) | https://github.com/MrOlm/inStrain authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [InterProScan](https://github.com/ebi-pf-team/interproscan) | https://github.com/ebi-pf-team/interproscan authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [InterVar](https://github.com/WGLab/InterVar) | https://github.com/WGLab/InterVar authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [IQ-TREE 2](https://github.com/iqtree/iqtree2) | https://github.com/iqtree/iqtree2 authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/molbev/msaa015) [supplied source 2](https://doi.org/10.1093/molbev/msx281) [supplied source 3](https://doi.org/10.1038/nmeth.4285) [supplied source 4](https://doi.org/10.1093/bioinformatics/btad540) | 来源记录；尚未集成 |
| [IQ-TREE 3](https://github.com/iqtree/iqtree3) | iqtree 项目及贡献者；论文作者见逐项原始引用 [paper 1](https://doi.org/10.1093/molbev/msag117) | Related local workflows: phylo-modern-tree-inference |
| [JCVI](https://github.com/tanghaibao/jcvi) | https://github.com/tanghaibao/jcvi authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Juicebox MCP](https://github.com/aidenlab/juicebox-mcp) | aidenlab项目作者与贡献者；Copyright (c) 2016-2019 The Regents of the University of California；The above copyright notice and this permission notice shall be included in all  | MCP: juicebox; Related local workflows: hi-c-analysis-hic-visualization |
| [Just DNA Agents](https://github.com/dna-seq/just-dna-agents) | https://github.com/dna-seq/just-dna-agents authors and contributors  | 来源记录；尚未集成 |
| [Kaiju](https://github.com/bioinformatics-centre/kaiju) | https://github.com/bioinformatics-centre/kaiju authors and contributors; GPL-3.0；[paper 1](https://doi.org/10.1038/ncomms11257) | `dnakit.microbiome.classify_kaiju`（需用户提供 Kaiju 蛋白库和 taxonomy） |
| [KinDA](https://github.com/DNA-and-Natural-Algorithms-Group/KinDA) | DNA-and-Natural-Algorithms-Group/KinDA authors and contributors [paper 1](https://doi.org/10.1098/rsif.2018.0107) | 来源记录；尚未集成 |
| [Kleborate](https://github.com/klebgenomics/Kleborate) | https://github.com/klebgenomics/Kleborate authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [KMCP](https://github.com/shenwei356/kmcp) | https://github.com/shenwei356/kmcp authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/bioinformatics/btac845) | 来源记录；尚未集成 |
| [Kraken2](https://github.com/DerrickWood/kraken2) | DerrickWood项目作者与贡献者；Copyright (c) 2017-2023 Derrick Wood；The above copyright notice and this permission notice shall be included in all [paper 1](https://dx.doi.org/10.1186/gb-2014-15-3-r46) [paper 2](https://dx.doi.org/10.1186/s13059-018-1568-0) [paper 3](https://dx.doi.org/10.1186/s13059-019-1891-0) [paper 4](https://dx.doi.org/10.1038/s41596-022-00738-y) [supplied source 1](https://doi.org/10.1186/s13059-019-1891-0) | `dnakit.microbiome.classify_kraken2`（需用户提供 Kraken2、数据库）; Related local workflows: metagenomics-kraken |
| [Liftoff](https://github.com/agshumate/Liftoff) | https://github.com/agshumate/Liftoff authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btaa1016) | liftoff-annotation; Related local workflows: genome-annotation-annotation-transfer |
| [Lindel](https://github.com/shendurelab/Lindel) | https://github.com/shendurelab/Lindel authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [LIRICAL](https://github.com/TheJacksonLaboratory/LIRICAL) | https://github.com/TheJacksonLaboratory/LIRICAL authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Locus](https://github.com/mattsmcnulty/locus) | Matthew McNulty (MIT copyright) and Locus contributors; underlying models/data have separate citations.  | Related local workflows: ancestry-risk-profiler |
| [LUMPY](https://github.com/arq5x/lumpy-sv) | https://github.com/arq5x/lumpy-sv authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MACS3](https://github.com/macs3-project/MACS) | MACS3 authors and contributors；调查表记录 BSD-3-Clause，未提供 DOI | `dnakit.ngs.call_chip_peaks` |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Katoh & Standley (2013), [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010); BSD | `dnakit.alignment.align_sequences`；`dnakit.bioinformatics.alignment.align_sequences` |
| [maftools](https://github.com/PoisonAlien/maftools) | https://github.com/PoisonAlien/maftools authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) | MAGeCK authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1186/s13059-014-0554-4) | mageck-test; mageck-count |
| [Manta](https://github.com/Illumina/manta) | https://github.com/Illumina/manta authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [mapDamage](https://github.com/ginolhac/mapDamage) | ginolhac 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: grapevine-adna |
| [MAPpoly](https://github.com/mmollina/mappoly) | mmollina 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: polyploid-gene-family-cloning |
| [MaveDB API](https://github.com/VariantEffect/mavedb-api) | https://github.com/VariantEffect/mavedb-api authors and contributors [paper 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1845-6) | 来源记录；尚未集成 |
| [MaveDB Link](https://github.com/berntpopp/mavedb-link) | https://github.com/berntpopp/mavedb-link authors and contributors [paper 1](https://doi.org/10.5281/zenodo.11201736) [paper 2](https://doi.org/10.1186/s13059-019-1845-6) | MCP: mavedb-link |
| [MAVE-NN](https://github.com/jbkinney/mavenn) | https://github.com/jbkinney/mavenn authors and contributors [paper 1](https://doi.org/10.1186/s13059-022-02661-7) | 来源记录；尚未集成 |
| [Medea](https://github.com/mims-harvard/Medea) | https://github.com/mims-harvard/Medea authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MEGAHIT](https://github.com/voutcn/megahit) | https://github.com/voutcn/megahit authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Megalodon](https://github.com/nanoporetech/megalodon) | https://github.com/nanoporetech/megalodon authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MEME Suite / FIMO](https://meme-suite.org/meme/) | Charles E. Grant, Timothy L. Bailey, William S. Noble; MEME Suite contributors [paper 1](https://doi.org/10.1093/bioinformatics/btr064) | fimo-motifs |
| [memesuite-lite](https://github.com/jmschrei/memesuite-lite) | https://github.com/jmschrei/memesuite-lite authors and contributors [paper 1](https://www.biorxiv.org/content/10.1101/2025.05.27.656386v1) | Related local workflows: chipseq-motif-analysis, motif-search |
| [Metagenomics Skills](https://github.com/PKU-EMBL/Metagenomics-Skills) | PKU-EMBL项目作者与贡献者；Copyright (c) 2026 PKU-EMBL；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.1038/s41467-024-46539-7) [paper 2](https://doi.org/10.1186/s40168-018-0541-1) [paper 3](https://doi.org/10.1038/s41564-018-0171-1) [paper 4](https://doi.org/10.1038/s41587-023-01688-w) [paper 5](https://doi.org/10.7554/eLife.65088) [paper 6](https://doi.org/10.1038/s41587-020-0603-3) [paper 7](https://doi.org/10.1038/s41587-020-0501-8) [paper 8](https://doi.org/10.1093/nar/gkac1080) | 53 bundled Skills; Related local workflows: basalt, checkm2, metabat2, mgnify |
| [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) | biobakery项目作者与贡献者 [paper 1](https://doi.org/10.1038/s41587-023-01688-w) [paper 2](https://doi.org/10.1101/2024.02.19.580813) | Related local workflows: biobakery-databases, metagenomics-metaphlan, metagenomics-strain-tracking |
| [metaWRAP](https://github.com/bxlab/metaWRAP) | https://github.com/bxlab/metaWRAP authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [methplotlib](https://github.com/wdecoster/methplotlib) | https://github.com/wdecoster/methplotlib authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MethylDackel](https://github.com/dpryan79/MethylDackel) | https://github.com/dpryan79/MethylDackel authors and contributors; individual paper authors are not inferred | `dnakit.epigenome.extract_methylation_methyldackel`（需用户提供 MethylDackel） |
| [MethylGPT](https://github.com/albert-ying/MethylGPT) | albert-ying项目作者与贡献者；"Licensor" shall mean the copyright owner or entity authorized by；the copyright owner that is granting the License. [paper 1](https://www.biorxiv.org/content/10.1101/2024.10.30.621013v2) | Related local workflows: methylation-epigenetic-clocks |
| [methylKit](https://github.com/al2na/methylKit) | al2na项目作者与贡献者  | Related local workflows: methylation-dmr-detection, methylation-methylkit |
| [MGnify Genomes Catalogue Pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) | EBI-Metagenomics项目作者与贡献者；Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>；(1) assert copyright on the software, and (2) offer you this License [paper 1](https://doi.org/10.1016/j.jmb.2023.168016) | Related local workflows: mgnify, microbial-genome-characterization |
| [MGnify MCP](https://github.com/MGS-sails/mgnify-mcp) | MGS-sails项目作者与贡献者  | Related local workflows: mgnify |
| [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html) | minfi authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [minigraph](https://github.com/lh3/minigraph) | lh3/minigraph authors and contributors [paper 1](https://doi.org/10.1186/s13059-020-02168-z) | minigraph-pangenome |
| [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) | ComparativeGenomicsToolkit 项目及贡献者；论文作者见逐项原始引用 [paper 1](https://doi.org/10.1038/s41587-023-01793-w) | Related local workflows: comparative-genomics-pangenome-analysis, vitis-pangenome |
| [Minimap2](https://github.com/lh3/minimap2) | Li (2018), [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191); MIT | `dnakit.alignment.map_long_reads`；`dnakit.annotation.genes.annotate_genes` |
| [miniprot](https://github.com/lh3/miniprot) | https://github.com/lh3/miniprot authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btad014) | miniprot-align; Related local workflows: genome-annotation-annotation-transfer |
| [mlst](https://github.com/tseemann/mlst) | https://github.com/tseemann/mlst authors and contributors; individual paper authors are not inferred  | mlst-typing |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | https://github.com/soedinglab/MMseqs2 authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1038/nbt.3988) | `dnakit.alignment.search_mmseqs_protein`; `dnakit.alignment.cluster_mmseqs_protein`（需用户提供 MMseqs2） |
| [MMSplice / MTSplice](https://github.com/gagneurlab/MMSplice_MTSplice) | https://github.com/gagneurlab/MMSplice_MTSplice authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MoChA](https://github.com/freeseek/mocha) | https://github.com/freeseek/mocha authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Modkit](https://github.com/nanoporetech/modkit) | nanoporetech项目作者与贡献者；non-exclusive license under Contributor copyrights Licensable by such；Section 2.1 above, no copyright license is granted:  | dnakit.genomics.read_bedmethyl; dnakit.genomics.summarize_methylation; modkit-pileup; Related local workflows: long-read-sequencing-nanopore-methylation, methylation-calling |
| [moments](https://github.com/MomentsLD/moments) | https://github.com/MomentsLD/moments authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MOODS](https://github.com/jhkorhonen/MOODS) | https://github.com/jhkorhonen/MOODS authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [mosdepth](https://github.com/brentp/mosdepth) | brentp/mosdepth authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btx699) | mosdepth-coverage |
| [motifmatchr](https://bioconductor.org/packages/release/bioc/html/motifmatchr.html) | motifmatchr authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [mOTUs](https://github.com/motu-tool/mOTUs) | https://github.com/motu-tool/mOTUs authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MPRAflow](https://github.com/shendurelab/MPRAflow) | https://github.com/shendurelab/MPRAflow authors and contributors  | Related local workflows: functional-screen-analysis |
| [MrBayes](https://github.com/NBISweden/MrBayes) | https://github.com/NBISweden/MrBayes authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MSMC2](https://github.com/stschiff/msmc2) | https://github.com/stschiff/msmc2 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [msprime](https://github.com/tskit-dev/msprime) | tskit-dev 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: population-genetics-scikit-allel-analysis |
| [MToolBox](https://github.com/mitoNGS/MToolBox) | https://github.com/mitoNGS/MToolBox authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [MultiQC](https://github.com/MultiQC/MultiQC) | Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354); GPL-3.0-or-later | `dnakit.sequencing_qc.run_multiqc`（需用户提供 MultiQC 可执行程序） |
| [Multistrand](https://github.com/DNA-and-Natural-Algorithms-Group/multistrand) | DNA-and-Natural-Algorithms-Group/multistrand authors and contributors [paper 1](https://doi.org/10.1007/978-3-319-21999-8_13) | 来源记录；尚未集成 |
| [MUSCLE](https://github.com/rcedgar/muscle) | Edgar (2022), [10.1038/s41467-022-34630-w](https://doi.org/10.1038/s41467-022-34630-w); GPL-3.0 | `dnakit.alignment.align_sequences`；`dnakit.bioinformatics.alignment.align_sequences` |
| [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) | https://github.com/UMCUGenetics/MutationalPatterns authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [NanoMethViz](https://github.com/Shians/NanoMethViz) | https://github.com/Shians/NanoMethViz authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Nanopolish](https://github.com/jts/nanopolish) | https://github.com/jts/nanopolish authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) | Camacho et al. (2009), [10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421); software public-domain statement | `dnakit.alignment.build_blast_database`；`dnakit.alignment.local_blast` |
| [NCBI Database MCP (hpend2373)](https://github.com/hpend2373/NCBI-Database-MCP) | https://github.com/hpend2373/NCBI-Database-MCP authors and contributors; individual paper authors are not inferred  | MCP: ncbi-database |
| [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) | NCBI Datasets authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [NCBI EGAPx](https://github.com/ncbi/egapx) | https://github.com/ncbi/egapx authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Nextflow](https://github.com/nextflow-io/nextflow) | https://github.com/nextflow-io/nextflow authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/nbt.3820) | 来源记录；尚未集成 |
| [nf-core/atacseq](https://nf-co.re/atacseq) | nf-core/atacseq authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [nf-core/bacass](https://nf-co.re/bacass) | nf-core/bacass authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [nf-core/chipseq](https://nf-co.re/chipseq) | nf-core/chipseq authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [nf-core/mag](https://nf-co.re/mag) | nf-core/mag authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [nf-core/methylong](https://github.com/nf-core/methylong) | nf-core项目作者与贡献者；Copyright (c) The nf-core/methylong team；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.5281/zenodo.15366449) [paper 2](https://doi.org/10.5281/zenodo.15366448) [paper 3](https://dx.doi.org/10.1038/s41587-020-0439-x) | Related local workflows: long-read-sequencing-nanopore-methylation |
| [nf-core/methylseq](https://github.com/nf-core/methylseq) | nf-core项目作者与贡献者；Copyright (c) The nf-core/methylseq team；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.5281/zenodo.1343417) [paper 2](https://dx.doi.org/10.1038/s41587-020-0439-x) | Related local workflows: methylation-bismark-alignment |
| [nf-core/sarek](https://github.com/nf-core/sarek) | nf-core项目作者与贡献者；Copyright (c) The nf-core/sarek team；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.5281/zenodo.3476425) [paper 2](https://doi.org/10.1093/nargab/lqae031) [paper 3](http://dx.doi.org/10.12688/f1000research.16665.2) [paper 4](https://dx.doi.org/10.1038/s41587-020-0439-x) | Related local workflows: nfcore-sarek-wrapper, workflows-fastq-to-variants |
| [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler) | nf-core项目作者与贡献者；Copyright (c) The nf-core/taxprofiler team；The above copyright notice and this permission notice shall be included in all [paper 1](https://doi.org/10.5281/zenodo.7728364) [paper 2](https://doi.org/10.1101/2023.10.20.563221) [paper 3](https://dx.doi.org/10.1038/s41587-020-0439-x) | Related local workflows: metagenomics-kraken, metagenomics-metaphlan |
| [nf-core/eager](https://github.com/nf-core/eager) | nf-core 项目及贡献者；论文作者见逐项原始引用 [paper 1](https://doi.org/10.7717/peerj.10947) | Related local workflows: grapevine-adna |
| [NOREC4DNA](https://github.com/umr-ds/NOREC4DNA) | umr-ds/NOREC4DNA authors and contributors [paper 1](https://doi.org/10.1186/s12859-021-04318-x) | 来源记录；尚未集成 |
| [Nucleotide Transformer / AgroNT / SegmentNT](https://github.com/instadeepai/nucleotide-transformer) | https://github.com/instadeepai/nucleotide-transformer authors and contributors; individual paper authors are not inferred  [supplied source 1](https://www.nature.com/articles/s41592-024-02523-z) | dnakit.predictions.api.predict_sequence_properties; Related local workflows: nucleotide-transformer |
| [Nuskell](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell) | DNA-and-Natural-Algorithms-Group/nuskell authors and contributors  | 来源记录；尚未集成 |
| [OBITools4](https://github.com/metabarcoding/obitools4) | metabarcoding 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: ecological-genomics-edna-metabarcoding, workflows-edna-pipeline |
| [ODGI](https://github.com/pangenome/odgi) | https://github.com/pangenome/odgi authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [OKseqHMM](https://github.com/CL-CHEN-Lab/OK-Seq) | https://github.com/CL-CHEN-Lab/OK-Seq authors and contributors [paper 1](https://doi.org/10.5281/zenodo.7428883) [paper 2](https://doi.org/10.1038/ncomms10208) [paper 3](https://doi.org/10.1038/s41467-020-17858-2) [paper 4](https://doi.org/10.1093/nar/gkac1239) [paper 5](https://doi.org/10.1038/s41596-022-00793-5) | 来源记录；尚未集成 |
| [Omics Skills](https://github.com/fmschulz/omics-skills) | fmschulz项目作者与贡献者；Copyright (c) 2026 Frederik Schulz；The above copyright notice and this permission notice shall be included in all  | 10 bundled Skills; Related local workflows: microbial-genome-characterization, viromics |
| [OmicsClaw](https://github.com/TianGzlab/OmicsClaw) | TianGzlab/OmicsClaw authors and contributors  | 来源记录；尚未集成 |
| [Open Targets MCP](https://github.com/nickzren/opentargets-mcp) | https://github.com/nickzren/opentargets-mcp authors and contributors; individual paper authors are not inferred  | MCP: opentargets |
| [OpenCRAVAT](https://github.com/KarchinLab/open-cravat) | https://github.com/KarchinLab/open-cravat authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Operon](https://github.com/swaruplab/operon) | https://github.com/swaruplab/operon authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [OptiType](https://github.com/FRED-2/OptiType) | https://github.com/FRED-2/OptiType authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Orca](https://github.com/jzhoulab/orca) | https://github.com/jzhoulab/orca authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1038/s41588-022-01065-4) | 来源记录；尚未集成 |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | https://github.com/davidemms/OrthoFinder authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [oxDNA](https://github.com/lorenzo-rovigatti/oxDNA) | lorenzo-rovigatti/oxDNA authors and contributors  | 来源记录；尚未集成 |
| [pairtools](https://github.com/open2c/pairtools) | https://github.com/open2c/pairtools authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Panaroo](https://github.com/gtonkinhill/panaroo) | https://github.com/gtonkinhill/panaroo authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Pangolin splice predictor](https://github.com/tkzeng/Pangolin) | https://github.com/tkzeng/Pangolin authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) | https://github.com/PacificBiosciences/pb-CpG-tools authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [PCAngsd](https://github.com/Rosemeis/pcangsd) | https://github.com/Rosemeis/pcangsd authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Peppercorn](https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator) | DNA-and-Natural-Algorithms-Group/peppercornenumerator authors and contributors  | 来源记录；尚未集成 |
| [PGGB](https://github.com/pangenome/pggb) | pangenome 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: comparative-genomics-pangenome-analysis, vitis-pangenome |
| [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) | https://github.com/PGScatalog/pgsc_calc authors and contributors; individual paper authors are not inferred  | Related local workflows: clinical-databases-polygenic-risk |
| [PharmCAT](https://github.com/PharmGKB/PharmCAT) | PharmGKB项目作者与贡献者；applicable copyright doctrines of fair use, fair dealing, or other；(including copyright notices, patent notices, disclaimers of warranty, [paper 1](https://dx.doi.org/10.1002/cpt.928) | Related local workflows: clinical-databases-pharmacogenomics, clinpgx |
| [Piperine](https://github.com/DNA-and-Natural-Algorithms-Group/piperine) | DNA-and-Natural-Algorithms-Group/piperine authors and contributors  | 来源记录；尚未集成 |
| [Plant Genomics MCP](https://github.com/musharna/plant-genomics-mcp) | https://github.com/musharna/plant-genomics-mcp authors and contributors; individual paper authors are not inferred  | MCP: plant-genomics |
| [PLINK 2](https://github.com/chrchang/plink-ng) | chrchang/plink-ng authors and contributors  | plink2-pca; plink2-king; plink2-ld-prune; plink2-gwas |
| [polyRAD](https://github.com/lvclark/polyRAD) | Lindsay V. Clark, Alexander E. Lipka, Erik J. Sacks [paper 1](https://doi.org/10.1534/g3.118.200913) | Related local workflows: polyploid-gene-family-cloning |
| [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) | https://github.com/labgem/PPanGGOLiN authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Precision Medicine MCP](https://github.com/lynnlangit/precision-medicine-mcp) | https://github.com/lynnlangit/precision-medicine-mcp authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Predmoter](https://github.com/weberlab-hhu/Predmoter) | https://github.com/weberlab-hhu/Predmoter authors and contributors [paper 1](https://doi.org/10.1038/nmeth.2688) [paper 2](https://doi.org/10.1038/nature01080) [paper 3](https://doi.org/10.1186/GB-2008-9-9-R137) [paper 4](https://doi.org/10.1093/bioadv/vbae074) | 来源记录；尚未集成 |
| [PrimeDesign](https://github.com/pinellolab/PrimeDesign) | https://github.com/pinellolab/PrimeDesign authors and contributors  | Related local workflows: crispr-screens-prime-editing-screens |
| [Primer3](https://github.com/primer3-org/primer3) | https://github.com/primer3-org/primer3 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [primer3-py](https://github.com/libnano/primer3-py) | https://github.com/libnano/primer3-py authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ProCapNet](https://github.com/kundajelab/ProCapNet) | https://github.com/kundajelab/ProCapNet authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Prodigal](https://github.com/hyattpd/Prodigal) | https://github.com/hyattpd/Prodigal authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1186/1471-2105-11-119) | 来源记录；尚未集成 |
| [Prokka](https://github.com/tseemann/prokka) | tseemann/prokka authors and contributors [paper 1](https://doi.org/10.1093/bioinformatics/btu153) | `dnakit.annotation.annotate_prokka`（需用户提供 Prokka、数据库） |
| [PSMC](https://github.com/lh3/psmc) | https://github.com/lh3/psmc authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [PubMed MCP (cyanheads)](https://github.com/cyanheads/pubmed-mcp-server) | https://github.com/cyanheads/pubmed-mcp-server authors and contributors; individual paper authors are not inferred  | MCP: pubmed |
| [PureCN](https://github.com/lima1/PureCN) | https://github.com/lima1/PureCN authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [pybedtools](https://github.com/daler/pybedtools) | Dale et al. (2011), [10.1093/bioinformatics/btr539](https://doi.org/10.1093/bioinformatics/btr539); BEDTools [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033); MIT | `dnakit.integrations.intersect_pybedtools_bed`（需 `dnakit[bio-pybedtools]` 和外部 bedtools） |
| [pyBigWig](https://github.com/deeptools/pyBigWig) | pyBigWig software release, [DOI](https://doi.org/10.5281/zenodo.45238); MIT | `dnakit.integrations.summarize_pybigwig_region`（需 `dnakit[bio-pybigwig]`） |
| [PyBrOpS](https://github.com/rzshrote/pybrops) | Robert Shrote and PyBrOpS contributors (MIT copyright); full paper author list at linked PMC article. [paper 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11457082/) | Related local workflows: plant-genomics |
| [pydna](https://github.com/pydna-group/pydna) | https://github.com/pydna-group/pydna authors and contributors; individual paper authors are not inferred  | dnakit.bioinformatics.cloning.assemble_duplexes; dnakit.bioinformatics.cloning.read_constructs |
| [pyfaidx](https://github.com/mdshw5/pyfaidx) | Shirley et al. (2015), [10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1)；源代码 BSD-3-Clause；预印本另按 CC-BY-4.0 说明 | `dnakit.integrations.fetch_pyfaidx_region`（需 `dnakit[bio-pyfaidx]`） |
| [pyfastx](https://github.com/lmdu/pyfastx) | Du et al. (2021), [paper](https://doi.org/10.1093/bib/bbaa368); MIT | `dnakit.integrations.summarize_pyfastx`（需 `dnakit[bio-pyfastx]`） |
| [PyPGx](https://github.com/sbslee/pypgx) | https://github.com/sbslee/pypgx authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [PyRanges0](https://github.com/pyranges/pyranges0) | https://github.com/pyranges/pyranges0 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Pyrodigal](https://github.com/althonos/pyrodigal) | https://github.com/althonos/pyrodigal authors and contributors [paper 1](https://doi.org/10.21105/joss.04296) [paper 2](https://doi.org/10.1186/1471-2105-11-119) | dnakit.genomics.predict_prokaryotic_genes; Related local workflows: genome-annotation-prokaryotic-annotation |
| [pysam](https://github.com/pysam-developers/pysam) | https://github.com/pysam-developers/pysam authors and contributors; individual paper authors are not inferred  | dnakit.bioinformatics.reads.alignment_summary; dnakit.bioinformatics.reads.coverage_consensus |
| [pySBOL3](https://github.com/SynBioDex/pySBOL3) | https://github.com/SynBioDex/pySBOL3 authors and contributors [paper 1](https://doi.org/10.1021/acssynbio.2c00249) | 来源记录；尚未集成 |
| [QIIME2 / Rachis framework](https://github.com/rachis-org/rachis) | https://github.com/rachis-org/rachis authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [QUAST](https://github.com/ablab/quast) | https://github.com/ablab/quast authors and contributors; individual paper authors are not inferred [paper 1](https://doi.org/10.1093/bioinformatics/btt086) | quast-assembly |
| [R/qtl](https://rqtl.org/) | R/qtl authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [R/qtl2](https://rqtl.org/qtl2/) | R/qtl2 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [RAxML-NG](https://github.com/amkozlov/raxml-ng) | amkozlov 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: phylo-modern-tree-inference |
| [REGENIE](https://github.com/rgcgithub/regenie) | https://github.com/rgcgithub/regenie authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Remora](https://github.com/nanoporetech/remora) | https://github.com/nanoporetech/remora authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) | https://github.com/Dfam-consortium/RepeatMasker authors and contributors; individual paper authors are not inferred  | repeatmasker-library |
| [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) | https://github.com/Dfam-consortium/RepeatModeler authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [ResFinder](https://github.com/genomicepidemiology/resfinder) | https://github.com/genomicepidemiology/resfinder authors and contributors; individual paper authors are not inferred  | Related local workflows: disinfinder |
| [RFMix](https://github.com/slowkoni/rfmix) | slowkoni 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: ancestry-risk-profiler, population-genetics-population-structure |
| [RGI](https://github.com/arpcard/rgi) | https://github.com/arpcard/rgi authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [rMVP](https://github.com/xiaolei-lab/rMVP) | https://github.com/xiaolei-lab/rMVP authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1016/j.gpb.2020.10.007) | 来源记录；尚未集成 |
| [rrBLUP](https://cran.r-project.org/package=rrBLUP) | Jeffrey B. Endelman [paper 1](https://doi.org/10.3835/plantgenome2011.08.0024) [paper 2](https://doi.org/10.1534/g3.112.004259) [paper 3](https://doi.org/10.3168/jds.2007-0980) | dnakit.genomics.genomic_prediction; Related local workflows: workflows-gwas-pipeline; 独立实现固定惩罚的岭回归；不调用 R rrBLUP，也不估计 REML |
| [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) | Rsamtools authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Rule Set 3](https://github.com/gpp-rnd/rs3) | https://github.com/gpp-rnd/rs3 authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SAIGE](https://github.com/saigegit/SAIGE) | https://github.com/saigegit/SAIGE authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SAMtools](https://github.com/samtools/samtools) | Danecek et al. (2021), [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008); MIT | `dnakit.ngs.prepare_bam`（需用户提供 samtools）；Related local workflows: cite-encode, samtools-bam-processing |
| [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package) | UC-Davis-molecular-computing/scadnano-python-package authors and contributors [paper 1](https://doi.org/10.4230/LIPIcs.DNA.2020.9) | dnakit.genomics.design_dna_duplex; dnakit.genomics.export_oxdna_design; Related local workflows: dnakit-dna-nanotechnology |
| [scarHRD](https://github.com/sztup/scarHRD) | https://github.com/sztup/scarHRD authors and contributors [paper 1](https://www.nature.com/articles/s41523-018-0066-6) | 来源记录；尚未集成 |
| [SCENIC+](https://github.com/aertslab/scenicplus) | https://github.com/aertslab/scenicplus authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [schmutzi](https://github.com/grenaud/schmutzi) | grenaud 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: grapevine-adna |
| [scikit-allel](https://github.com/cggh/scikit-allel) | https://github.com/cggh/scikit-allel authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Sei](https://github.com/FunctionLab/sei-framework) | https://github.com/FunctionLab/sei-framework authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SemiBin](https://github.com/BigDataBiology/SemiBin) | https://github.com/BigDataBiology/SemiBin authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SeqKit](https://github.com/shenwei356/seqkit) | shenwei356项目作者与贡献者；Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.；The above copyright notice and this permission notice shall be included in [paper 1](https://doi.org/10.1002/imt2.191) | Related local workflows: read-sequences |
| [seqtk](https://github.com/lh3/seqtk) | seqtk 项目作者；MIT；官方未提供 DOI | `dnakit.ngs.sample_sequences_seqtk`（需用户提供 seqtk） |
| [SeSAMe](https://bioconductor.org/packages/release/bioc/html/sesame.html) | SeSAMe authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Signac](https://github.com/stuart-lab/signac) | https://github.com/stuart-lab/signac authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment) | https://github.com/SigProfilerSuite/SigProfilerAssignment authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) | AlexandrovLab项目作者与贡献者；Copyright (c) 2019, S M Ashiqul Islam (Mishu) [Alexandrov Lab]；1. Redistributions of source code must retain the above copyright notice, this [paper 1](https://doi.org/10.1016/j.xgen.2022.100179) | Related local workflows: clinical-databases-somatic-signatures |
| [SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) | https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SimBiology MCP](https://github.com/Sepanta-Yalameha/SIMBIOLOGY-MCP) | https://github.com/Sepanta-Yalameha/SIMBIOLOGY-MCP authors and contributors  | 来源记录；尚未集成 |
| [SLiM](https://github.com/MesserLab/SLiM) | Benjamin C. Haller (copyright) and MesserLab/SLiM contributors; use version-specific citation from SLiM manual.  | Related local workflows: ecological-genomics-conservation-genetics |
| [SMC++](https://github.com/popgenmethods/smcpp) | https://github.com/popgenmethods/smcpp authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Snakemake](https://github.com/snakemake/snakemake) | https://github.com/snakemake/snakemake authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SnapATAC2](https://github.com/scverse/snapatac2) | scverse项目作者与贡献者；Copyright (c) 2022-2026 Kai Zhang；The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. [paper 1](https://doi.org/10.1038/s41592-023-02139-9) | Related local workflows: atac-seq-single-cell-atac |
| [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | Smolka et al. (2024), [10.1038/s41587-023-02024-y](https://doi.org/10.1038/s41587-023-02024-y)；MIT | `dnakit.variants.call_sniffles2`（需用户提供 Sniffles2） |
| [Snippy](https://github.com/tseemann/snippy) | https://github.com/tseemann/snippy authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SnpEff / SnpSift](https://github.com/pcingola/SnpEff) | https://github.com/pcingola/SnpEff authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [sommer](https://github.com/covaruber/sommer) | Giovanny Covarrubias-Pazaran [paper 1](https://doi.org/10.1371/journal.pone.0156744) | Related local workflows: causal-genomics-heritability-partitioning |
| [SPAdes](https://github.com/ablab/spades) | Bankevich et al. (2012), [10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021); GPL-2.0 | `dnakit.ngs.assemble_reads`（需用户提供 spades.py）；Related local workflows: genome-assembly-short-read-assembly |
| [SpliceAI](https://github.com/Illumina/SpliceAI) | https://github.com/Illumina/SpliceAI authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Statistical Genetics Skills](https://github.com/svdorn/statgen-skills) | svdorn/statgen-skills authors and contributors  | 10 bundled Skills |
| [stdpopsim](https://github.com/popsim-consortium/stdpopsim) | https://github.com/popsim-consortium/stdpopsim authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Straglr](https://github.com/BirolLab/straglr) | https://github.com/BirolLab/straglr authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [STRait Razor](https://github.com/Ahhgust/STRaitRazor) | Ahhgust 项目及贡献者；论文作者见逐项原始引用  | Related local workflows: variant-calling |
| [Strelka2](https://github.com/Illumina/strelka) | https://github.com/Illumina/strelka authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [STRling](https://github.com/quinlan-lab/STRling) | https://github.com/quinlan-lab/STRling authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1186/s13059-022-02826-4) | 来源记录；尚未集成 |
| [susieR](https://github.com/stephenslab/susieR) | https://github.com/stephenslab/susieR authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1111/rssb.12388) [supplied source 2](https://doi.org/10.1371/journal.pgen.1010299) | 来源记录；尚未集成 |
| [SVIM](https://github.com/eldariont/svim) | https://github.com/eldariont/svim authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [SynBioHub](https://github.com/SynBioHub/synbiohub) | https://github.com/SynBioHub/synbiohub authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [tangermeme](https://github.com/jmschrei/tangermeme) | https://github.com/jmschrei/tangermeme authors and contributors [paper 1](https://www.biorxiv.org/content/10.1101/2025.08.08.669296v2) | 1 bundled Skills |
| [TASSEL](https://tassel.bitbucket.io/) | TASSEL authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Telomerecat](https://github.com/cancerit/telomerecat) | https://github.com/cancerit/telomerecat authors and contributors  | 来源记录；尚未集成 |
| [TelomereHunter2](https://github.com/ferdinand-popp/telomerehunter2) | https://github.com/ferdinand-popp/telomerehunter2 authors and contributors [paper 1](https://doi.org/10.1093/bioadv/vbag187) [paper 2](https://www.biorxiv.org/content/10.1101/2024.08.28.609339v2) [paper 3](https://doi.org/10.1186/s12859-019-2851-0) | 来源记录；尚未集成 |
| [TF-MoDISco](https://github.com/kundajelab/tfmodisco) | https://github.com/kundajelab/tfmodisco authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [TOBIAS](https://github.com/loosolab/TOBIAS) | https://github.com/loosolab/TOBIAS authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Tombo](https://github.com/nanoporetech/tombo) | https://github.com/nanoporetech/tombo authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Trim Galore](https://github.com/FelixKrueger/TrimGalore) | https://github.com/FelixKrueger/TrimGalore authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [TRTools](https://github.com/gymrek-lab/TRTools) | https://github.com/gymrek-lab/TRTools authors and contributors; individual paper authors are not inferred  [supplied source 1](https://doi.org/10.1093/bioinformatics/btaa736) | 来源记录；尚未集成 |
| [Truvari](https://github.com/ACEnglish/truvari) | https://github.com/ACEnglish/truvari authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [tskit](https://github.com/tskit-dev/tskit) | https://github.com/tskit-dev/tskit authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | https://github.com/MRCIEU/TwoSampleMR authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [TxAgent](https://github.com/mims-harvard/TxAgent) | https://github.com/mims-harvard/TxAgent authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [UCSC Genome MCP](https://github.com/hlydecker/ucsc-genome-mcp) | https://github.com/hlydecker/ucsc-genome-mcp authors and contributors; individual paper authors are not inferred  | MCP: ucsc-genome |
| [UNACORM](https://github.com/AAnzel/UNACORM) | AAnzel/UNACORM authors and contributors [paper 1](https://arxiv.org/abs/2608.09673) | 来源记录；尚未集成 |
| [Unicycler](https://github.com/rrwick/Unicycler) | https://github.com/rrwick/Unicycler authors and contributors; GPL-3.0；[paper 1](https://doi.org/10.1371/journal.pcbi.1005595) | `dnakit.assembly.assemble_unicycler`（需用户提供 Unicycler） |
| [updog](https://github.com/dcgerard/updog) | David Gerard, Luís Felipe Ventorim Ferrão, Antonio Augusto Franco Garcia, Matthew Stephens [paper 1](https://doi.org/10.1534/genetics.118.301468) | Related local workflows: polyploid-gene-family-cloning |
| [VAMB](https://github.com/RasmussenLab/vamb) | https://github.com/RasmussenLab/vamb authors and contributors; individual paper authors are not inferred  | Related local workflows: taxvamb |
| [VanRaden genomic relationship method 1](https://doi.org/10.3168/jds.2007-0980) | Paul M. VanRaden [paper 1](https://doi.org/10.3168/jds.2007-0980) | dnakit.genomics.genomic_relationship_matrix |
| [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) | VariantAnnotation authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [Bioinformatics Agent Skills](https://github.com/variomeanalytics/bioinformatics-agent-skills) | variomeanalytics/bioinformatics-agent-skills authors and contributors  | 来源记录；尚未集成 |
| [vg](https://github.com/vgteam/vg) | https://github.com/vgteam/vg authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |
| [VSEARCH](https://github.com/torognes/vsearch) | Rognes et al. (2016), [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)；BSD-2-Clause 或 GPL-3.0 | `dnakit.fasta_fastq.cluster_vsearch`（需用户提供 VSEARCH） |
| [BamTools](https://github.com/pezmaster31/bamtools) | Barnett et al. (2011), [10.1093/bioinformatics/btr174](https://doi.org/10.1093/bioinformatics/btr174)；MIT | `dnakit.ngs.sort_bam_bamtools`（需用户提供 BamTools） |
| [WhatsHap](https://github.com/whatshap/whatshap) | Martin et al. (2016), [10.1101/085050](https://doi.org/10.1101/085050)；MIT | `dnakit.ngs.phase_variants_whatshap`（需用户提供 WhatsHap） |
| [Winnowmap](https://github.com/marbl/Winnowmap) | https://github.com/marbl/Winnowmap authors and contributors; individual paper authors are not inferred  | 来源记录；尚未集成 |

## 4) 如何查看证据、依赖与具体差距

```python
from dnakit.integrations import search_projects, describe_analysis
print(search_projects('BOLE'))
print(search_projects(domain='dna-storage'))
print(describe_analysis('plink2-pca'))
```

包内检索结果包含原始链接、许可证证据、已有技能匹配、差距与真实适配器 ID。仓库 registry/genomics-survey.json 及其分领域原始记录保存检索日志。使用 dnakit skills cite SKILL_ID 查看固定版本的工作流来源。外部软件、模型和数据适用其各自条款。
