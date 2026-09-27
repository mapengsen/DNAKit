# Review of the 335 supplied projects

All 335 rows were compared with DNAKit. The input contains software, models, workflows and MCP services, not 335 installable Skills.

- Existing project records / 已有项目记录：**82**；新增项目记录：**253**。
- Total project/source records / 合并后项目与来源：**419**。
- Bundled Skills / 分发 Skills：**573**，本次新增 **13**。
- Executable analyses / 外部分析：**43**（原 19，本次新增 **24**）；MCP connections / MCP 连接：**32**（原 7，本次新增 **25**）。
- Source pages retrieved / 取得原始文档响应：**340/340**（不等于运行验证）。

## 1) How to use the additions

```python
from dnakit.integrations import list_analyses, describe_analysis, search_projects, list_mcp_servers
print(list_analyses('DIAMOND'))
print(describe_analysis('diamond-blastp'))
print(search_projects('BOLE'))
print(list_mcp_servers())
```

Use plan_analysis/run_analysis for registered analyses and list_mcp_tools/call_mcp_tool for configured MCP services. Install upstream software and supply databases/checkpoints separately. See the external integrations guide.

## 2) Decisions for every supplied row

| # | Project / 项目 | Deduplication / 去重 | Local entrypoints / 本地能力 | Conditions / 条件与限制 |
| ---: | --- | --- | --- | --- |
| 1 | [AlphaMate](https://github.com/AlphaGenes/AlphaMate) | 合并已有 | 1 related workflow guides | 不要与实测表型预测混为一谈 |
| 2 | [AlphaSimR](https://gaynorr.github.io/AlphaSimR/) | 合并已有 | 1 related workflow guides |  |
| 3 | [BGLR](https://github.com/gdlc/BGLR-R) | 合并已有 | 1 related workflow guides |  |
| 4 | [GAPIT](https://github.com/jiabowang/GAPIT) | 合并已有 | 1 related workflow guides | README含GAPIT4与历史版本引文；不同版本分别引用 |
| 5 | [GEMMA](https://github.com/genetics-statistics/GEMMA) | 新增记录 | gemma-lmm |  |
| 6 | [PyBrOpS](https://github.com/rzshrote/pybrops) | 合并已有 | 1 related workflow guides |  |
| 7 | [R/qtl](https://rqtl.org/) | 新增记录 | Source only; no adapter |  |
| 8 | [R/qtl2](https://rqtl.org/qtl2/) | 新增记录 | Source only; no adapter |  |
| 9 | [rMVP](https://github.com/xiaolei-lab/rMVP) | 新增记录 | Source only; no adapter |  |
| 10 | [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html) | 合并已有 | dnakit.genomics.genomic_prediction; 1 related workflow guides |  |
| 11 | [sommer](https://cran.r-project.org/web/packages/sommer/index.html) | 合并已有 | 1 related workflow guides |  |
| 12 | [TASSEL](https://tassel.bitbucket.io/) | 新增记录 | Source only; no adapter | 需Java运行时 |
| 13 | [ABC enhancer-gene model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) | 新增记录 | 1 related workflow guides |  |
| 14 | [Basenji](https://github.com/calico/basenji) | 新增记录 | Source only; no adapter |  |
| 15 | [Borzoi](https://github.com/calico/borzoi) | 新增记录 | Source only; no adapter | 依赖baskerville；模型权重和训练数据单独管理 |
| 16 | [BPNet](https://github.com/kundajelab/bpnet) | 新增记录 | Source only; no adapter |  |
| 17 | [ChromBPNet](https://github.com/kundajelab/chrombpnet) | 合并已有 | 1 related workflow guides |  |
| 18 | [chromVAR](https://bioconductor.org/packages/release/bioc/html/chromVAR.html) | 新增记录 | Source only; no adapter |  |
| 19 | [DeepSTARR](https://github.com/bernardo-de-almeida/DeepSTARR) | 合并已有 | Source only; no adapter | 原始模型针对果蝇S2细胞；不能直接泛化到所有物种 |
| 20 | [Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) | 新增记录 | dnakit.predictions.api.predict_sequence_properties; dnakit.predictions.api.predict_enformer_benchmark | 官方研究仓库子目录；旧TensorFlow环境需隔离；使用Basenji数据和输出格式; Existing configured Enformer prediction backend; local checkpoints required. |
| 21 | [gReLU](https://github.com/Genentech/gReLU) | 新增记录 | Source only; no adapter |  |
| 22 | [MEME Suite](https://meme-suite.org/meme/) | 合并已有 | fimo-motifs | 套件按一项计数；子命令不冒充独立项目 |
| 23 | [MOODS](https://github.com/jhkorhonen/MOODS) | 新增记录 | Source only; no adapter |  |
| 24 | [motifmatchr](https://bioconductor.org/packages/release/bioc/html/motifmatchr.html) | 新增记录 | Source only; no adapter |  |
| 25 | [Orca](https://github.com/jzhoulab/orca) | 新增记录 | Source only; no adapter | 需源码、特定Selene分支、权重和参考资源；README未提供普通pip发行包；非商业和学术使用条款 |
| 26 | [ProCapNet](https://github.com/kundajelab/ProCapNet) | 新增记录 | Source only; no adapter |  |
| 27 | [SCENIC+](https://github.com/aertslab/scenicplus) | 新增记录 | Source only; no adapter |  |
| 28 | [Sei](https://github.com/FunctionLab/sei-framework) | 新增记录 | Source only; no adapter |  |
| 29 | [Tangermeme](https://github.com/jmschrei/tangermeme) | 合并已有 | 1 bundled Skills |  |
| 30 | [TF-MoDISco](https://github.com/kundajelab/tfmodisco) | 新增记录 | Source only; no adapter |  |
| 31 | [TOBIAS](https://github.com/loosolab/TOBIAS) | 新增记录 | Source only; no adapter |  |
| 32 | [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) | 新增记录 | augustus-genes |  |
| 33 | [Bakta](https://github.com/oschwengers/bakta) | 新增记录 | bakta-annotation |  |
| 34 | [BRAKER1-3](https://github.com/Gaius-Augustus/BRAKER) | 合并已有 | 1 related workflow guides | 原仓库标为deprecated，推荐BRAKER4；依赖组件另查授权 |
| 35 | [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) | 新增记录 | Source only; no adapter |  |
| 36 | [BUSCO](https://busco.ezlab.org/) | 新增记录 | busco-proteins; busco-genome |  |
| 37 | [DIAMOND](https://github.com/bbuchfink/diamond) | 新增记录 | diamond-blastp; diamond-blastx |  |
| 38 | [EDTA](https://github.com/oushujun/EDTA) | 新增记录 | Source only; no adapter |  |
| 39 | [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | 新增记录 | eggnog-proteins |  |
| 40 | [Helixer](https://github.com/usadellab/Helixer) | 合并已有 | 1 related workflow guides | 旧weberlab-hhu地址已迁移至usadellab |
| 41 | [InterProScan](https://github.com/ebi-pf-team/interproscan) | 新增记录 | Source only; no adapter | 属于DNA下游辅助，不是DNA序列模型 |
| 42 | [JCVI](https://github.com/tanghaibao/jcvi) | 新增记录 | Source only; no adapter |  |
| 43 | [Liftoff](https://github.com/agshumate/Liftoff) | 合并已有 | liftoff-annotation; 1 related workflow guides |  |
| 44 | [miniprot](https://github.com/lh3/miniprot) | 合并已有 | miniprot-align; 1 related workflow guides |  |
| 45 | [NCBI EGAPx](https://github.com/ncbi/egapx) | 新增记录 | Source only; no adapter |  |
| 46 | [OrthoFinder](https://github.com/davidemms/OrthoFinder) | 新增记录 | Source only; no adapter |  |
| 47 | [Prodigal](https://github.com/hyattpd/Prodigal) | 新增记录 | Source only; no adapter |  |
| 48 | [Prokka](https://github.com/tseemann/prokka) | 合并已有 | prokka-annotation |  |
| 49 | [Pyrodigal](https://github.com/althonos/pyrodigal) | 合并已有 | dnakit.genomics.predict_prokaryotic_genes; 1 related workflow guides |  |
| 50 | [QUAST](https://github.com/ablab/quast) | 新增记录 | quast-assembly |  |
| 51 | [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) | 新增记录 | repeatmasker-library |  |
| 52 | [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) | 新增记录 | Source only; no adapter |  |
| 53 | [BCFtools](https://github.com/samtools/bcftools) | 新增记录 | dnakit.bioinformatics.ngs.call_variants |  |
| 54 | [CADD](https://github.com/kircherlab/CADD-scripts) | 新增记录 | Source only; no adapter | 评分/数据库使用许可需单独核对 |
| 55 | [Clair3](https://github.com/HKU-BAL/Clair3) | 合并已有 | 1 related workflow guides |  |
| 56 | [CNVkit](https://github.com/etal/cnvkit) | 合并已有 | 1 related workflow guides |  |
| 57 | [Control-FREEC](https://github.com/BoevaLab/FREEC) | 新增记录 | Source only; no adapter |  |
| 58 | [cuteSV](https://github.com/tjiangHIT/cuteSV) | 新增记录 | Source only; no adapter |  |
| 59 | [DeepVariant](https://github.com/google/deepvariant) | 合并已有 | 1 related workflow guides | 通常用官方容器；不是普通pip包 |
| 60 | [DELLY](https://github.com/dellytools/delly) | 新增记录 | Source only; no adapter |  |
| 61 | [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) | 合并已有 | dnakit.annotation.variants.annotate_variant_vep; dnakit.annotation.variants.annotate_rsid_vep; 1 related workflow guides | ; Existing Ensembl REST annotation through the configured ToolUniverse provider; not a local VEP CLI installation. |
| 62 | [FreeBayes](https://github.com/freebayes/freebayes) | 新增记录 | Source only; no adapter |  |
| 63 | [GATK-SV](https://github.com/broadinstitute/gatk-sv) | 新增记录 | Source only; no adapter | 不是GATK4单个命令的简单同义名 |
| 64 | [GATK4](https://github.com/broadinstitute/gatk) | 合并已有 | 2 related workflow guides | 子工具不重复计为独立项目 |
| 65 | [hap.py](https://github.com/Illumina/hap.py) | 新增记录 | Source only; no adapter |  |
| 66 | [LUMPY](https://github.com/arq5x/lumpy-sv) | 新增记录 | Source only; no adapter |  |
| 67 | [Manta](https://github.com/Illumina/manta) | 新增记录 | Source only; no adapter | 历史工具；部署前核对维护与旧依赖 |
| 68 | [MMSplice / MTSplice](https://github.com/gagneurlab/MMSplice_MTSplice) | 新增记录 | Source only; no adapter |  |
| 69 | [MoChA](https://github.com/freeseek/mocha) | 新增记录 | Source only; no adapter | BCFtools插件；README明确不适用于WES或低覆盖WGS |
| 70 | [Pangolin splice predictor](https://github.com/tkzeng/Pangolin) | 新增记录 | Source only; no adapter | 不是同名病毒谱系分型工具 |
| 71 | [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | 合并已有 | sniffles-sv; 2 related workflow guides |  |
| 72 | [SnpEff / SnpSift](https://github.com/pcingola/SnpEff) | 新增记录 | Source only; no adapter | 同一项目的组件合并计数 |
| 73 | [SpliceAI](https://github.com/Illumina/SpliceAI) | 新增记录 | Source only; no adapter | 代码与模型授权分别检查 |
| 74 | [Strelka2](https://github.com/Illumina/strelka) | 新增记录 | Source only; no adapter | 历史工具；部署前核对维护与旧依赖 |
| 75 | [SVIM](https://github.com/eldariont/svim) | 新增记录 | Source only; no adapter |  |
| 76 | [Truvari](https://github.com/ACEnglish/truvari) | 新增记录 | Source only; no adapter |  |
| 77 | [WhatsHap](https://github.com/whatshap/whatshap) | 新增记录 | Source only; no adapter |  |
| 78 | [Aldy](https://github.com/0xTCG/aldy) | 新增记录 | Source only; no adapter |  |
| 79 | [AnnotSV](https://github.com/lgmgeo/AnnotSV) | 新增记录 | Source only; no adapter |  |
| 80 | [arcasHLA](https://github.com/RabadanLab/arcasHLA) | 新增记录 | Source only; no adapter | RNA输入；作为DNA分型辅助收录 |
| 81 | [ASCAT](https://github.com/VanLoo-lab/ascat) | 新增记录 | Source only; no adapter |  |
| 82 | [AutoPVS1](https://github.com/JiguangPeng/autopvs1) | 新增记录 | Source only; no adapter |  |
| 83 | [ClassifyCNV](https://github.com/Genotek/ClassifyCNV) | 新增记录 | Source only; no adapter |  |
| 84 | [Exomiser](https://github.com/exomiser/Exomiser) | 合并已有 | 2 related workflow guides | 不是自动确诊；参考数据库需固定版本 |
| 85 | [FACETS](https://github.com/mskcc/facets) | 新增记录 | Source only; no adapter |  |
| 86 | [Hartwig HMFtools](https://github.com/hartwigmedical/hmftools) | 新增记录 | Source only; no adapter | 多个组件；按一个项目计数 |
| 87 | [ichorCNA](https://github.com/broadinstitute/ichorCNA) | 新增记录 | Source only; no adapter |  |
| 88 | [InterVar](https://github.com/WGLab/InterVar) | 新增记录 | Source only; no adapter | ANNOVAR等依赖另核授权；不能替代人工判读 |
| 89 | [LIRICAL](https://github.com/TheJacksonLaboratory/LIRICAL) | 新增记录 | Source only; no adapter |  |
| 90 | [maftools](https://github.com/PoisonAlien/maftools) | 新增记录 | Source only; no adapter |  |
| 91 | [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) | 新增记录 | Source only; no adapter |  |
| 92 | [OpenCRAVAT](https://github.com/KarchinLab/open-cravat) | 新增记录 | Source only; no adapter |  |
| 93 | [OptiType](https://github.com/FRED-2/OptiType) | 新增记录 | Source only; no adapter |  |
| 94 | [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) | 新增记录 | 1 related workflow guides | 根据README准备标准化基因型输入；不是直接从原始WGS reads计算临床风险 |
| 95 | [PharmCAT](https://github.com/PharmGKB/PharmCAT) | 合并已有 | 2 related workflow guides | 不是仅凭单SNP即可完整分型的替代品 |
| 96 | [PureCN](https://github.com/lima1/PureCN) | 新增记录 | Source only; no adapter |  |
| 97 | [PyPGx](https://github.com/sbslee/pypgx) | 新增记录 | Source only; no adapter |  |
| 98 | [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment) | 新增记录 | Source only; no adapter |  |
| 99 | [SigProfilerExtractor](https://github.com/SigProfilerSuite/SigProfilerExtractor) | 合并已有 | 1 related workflow guides |  |
| 100 | [SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) | 新增记录 | Source only; no adapter |  |
| 101 | [Azimuth](https://github.com/MicrosoftResearch/Azimuth) | 新增记录 | Source only; no adapter | 旧模型与依赖需隔离 |
| 102 | [BAGEL2](https://github.com/hart-lab/bagel) | 新增记录 | 2 related workflow guides |  |
| 103 | [BE-Hive bystander model](https://github.com/maxwshen/be_predict_bystander) | 新增记录 | Source only; no adapter |  |
| 104 | [Cas-OFFinder](https://github.com/snugel/cas-offinder) | 合并已有 | cas-offinder; 1 related workflow guides |  |
| 105 | [CRISPOR](https://github.com/maximilianh/crisporWebsite) | 新增记录 | Source only; no adapter |  |
| 106 | [crisprDesign](https://github.com/crisprVerse/crisprDesign) | 新增记录 | Source only; no adapter |  |
| 107 | [CRISPResso2](https://github.com/pinellolab/CRISPResso2) | 合并已有 | crispresso2-amplicon; 1 related workflow guides |  |
| 108 | [CRISPRitz](https://github.com/pinellolab/CRISPRitz) | 合并已有 | 1 related workflow guides |  |
| 109 | [crisprScore](https://github.com/crisprVerse/crisprScore) | 合并已有 | 1 related workflow guides |  |
| 110 | [DrugZ](https://github.com/hart-lab/drugz) | 新增记录 | drugz-unpaired; 1 related workflow guides |  |
| 111 | [FlashFry](https://github.com/mckennalab/FlashFry) | 新增记录 | Source only; no adapter |  |
| 112 | [FORECasT / SelfTarget](https://github.com/felicityallen/SelfTarget) | 新增记录 | Source only; no adapter |  |
| 113 | [GuideScan CLI](https://github.com/pritykinlab/guidescan-cli) | 新增记录 | Source only; no adapter |  |
| 114 | [inDelphi](https://github.com/maxwshen/inDelphi-model) | 新增记录 | Source only; no adapter |  |
| 115 | [Lindel](https://github.com/shendurelab/Lindel) | 新增记录 | Source only; no adapter |  |
| 116 | [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) | 新增记录 | mageck-test; mageck-count |  |
| 117 | [PrimeDesign](https://github.com/pinellolab/PrimeDesign) | 合并已有 | 1 related workflow guides |  |
| 118 | [Rule Set 3](https://github.com/gpp-rnd/rs3) | 新增记录 | Source only; no adapter |  |
| 119 | [ADMIXTOOLS 2](https://github.com/uqrmaie1/admixtools) | 合并已有 | 1 related workflow guides |  |
| 120 | [ADMIXTURE](https://dalexander.github.io/admixture/) | 新增记录 | Source only; no adapter | 公开工具；源码及再分发授权待核 |
| 121 | [ANGSD](https://github.com/ANGSD/angsd) | 新增记录 | Source only; no adapter |  |
| 122 | [ASTRAL](https://github.com/smirarab/ASTRAL) | 新增记录 | Source only; no adapter |  |
| 123 | [coloc](https://github.com/chr1swallace/coloc) | 新增记录 | Source only; no adapter |  |
| 124 | [EIGENSOFT](https://github.com/DReichLab/EIG) | 新增记录 | Source only; no adapter |  |
| 125 | [fastSTRUCTURE](https://github.com/rajanil/fastStructure) | 新增记录 | Source only; no adapter | 旧Python依赖需核对 |
| 126 | [IQ-TREE 2](https://github.com/iqtree/iqtree2) | 新增记录 | Source only; no adapter |  |
| 127 | [MAFFT](https://mafft.cbrc.jp/alignment/software/) | 新增记录 | dnakit.bioinformatics.alignment.align_sequences | ; align_sequences(method='mafft'); external MAFFT required. |
| 128 | [moments](https://github.com/MomentsLD/moments) | 新增记录 | Source only; no adapter |  |
| 129 | [MrBayes](https://github.com/NBISweden/MrBayes) | 新增记录 | Source only; no adapter |  |
| 130 | [MSMC2](https://github.com/stschiff/msmc2) | 新增记录 | Source only; no adapter |  |
| 131 | [msprime](https://github.com/tskit-dev/msprime) | 合并已有 | 1 related workflow guides |  |
| 132 | [MUSCLE](https://github.com/rcedgar/muscle) | 新增记录 | dnakit.bioinformatics.alignment.align_sequences | ; align_sequences(method='muscle'); MUSCLE 5 required. |
| 133 | [PCAngsd](https://github.com/Rosemeis/pcangsd) | 新增记录 | Source only; no adapter |  |
| 134 | [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | 合并已有 | plink2-pca; plink2-king; plink2-ld-prune; plink2-gwas |  |
| 135 | [PSMC](https://github.com/lh3/psmc) | 新增记录 | Source only; no adapter |  |
| 136 | [RAxML-NG](https://github.com/amkozlov/raxml-ng) | 合并已有 | 1 related workflow guides |  |
| 137 | [REGENIE](https://github.com/rgcgithub/regenie) | 新增记录 | Source only; no adapter |  |
| 138 | [SAIGE](https://github.com/saigegit/SAIGE) | 新增记录 | Source only; no adapter |  |
| 139 | [scikit-allel](https://github.com/cggh/scikit-allel) | 新增记录 | Source only; no adapter |  |
| 140 | [SLiM](https://github.com/MesserLab/SLiM) | 合并已有 | 1 related workflow guides |  |
| 141 | [SMC++](https://github.com/popgenmethods/smcpp) | 新增记录 | Source only; no adapter |  |
| 142 | [stdpopsim](https://github.com/popsim-consortium/stdpopsim) | 新增记录 | Source only; no adapter |  |
| 143 | [susieR](https://github.com/stephenslab/susieR) | 新增记录 | Source only; no adapter |  |
| 144 | [tskit](https://github.com/tskit-dev/tskit) | 新增记录 | Source only; no adapter |  |
| 145 | [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | 新增记录 | Source only; no adapter | 部分数据接口依赖OpenGWAS；统计证据不等于已证实因果 |
| 146 | [ABRicate](https://github.com/tseemann/abricate) | 新增记录 | abricate-screen | 多个外部数据库各自有条款 |
| 147 | [AMRFinderPlus](https://github.com/ncbi/amr) | 合并已有 | amrfinder-nucleotide; 1 related workflow guides | 用于检测和监测；结果不替代药敏试验 |
| 148 | [ARIBA](https://github.com/sanger-pathogens/ariba) | 新增记录 | Source only; no adapter | README说明缺少资源提供支持；旧流程复现需隔离环境 |
| 149 | [Bracken](https://github.com/jenniferlu717/Bracken) | 合并已有 | 1 related workflow guides |  |
| 150 | [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) | 新增记录 | Source only; no adapter |  |
| 151 | [CheckM2](https://github.com/chklovski/CheckM2) | 新增记录 | checkm2-predict; 1 related workflow guides |  |
| 152 | [chewBBACA](https://github.com/B-UMMI/chewBBACA) | 新增记录 | Source only; no adapter |  |
| 153 | [DADA2](https://github.com/benjjneb/dada2) | 新增记录 | Source only; no adapter |  |
| 154 | [DAS Tool](https://github.com/cmks/DAS_Tool) | 新增记录 | Source only; no adapter |  |
| 155 | [geNomad](https://github.com/apcamargo/genomad) | 新增记录 | genomad-mobile-elements; 1 related workflow guides |  |
| 156 | [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | 合并已有 | 1 related workflow guides | 必须固定GTDB参考库版本 |
| 157 | [HUMAnN](https://github.com/biobakery/humann) | 新增记录 | 1 related workflow guides |  |
| 158 | [inStrain](https://github.com/MrOlm/inStrain) | 新增记录 | Source only; no adapter |  |
| 159 | [Kaiju](https://github.com/bioinformatics-centre/kaiju) | 新增记录 | Source only; no adapter |  |
| 160 | [Kleborate](https://github.com/klebgenomics/Kleborate) | 新增记录 | Source only; no adapter | 旧katholt地址已迁移 |
| 161 | [KMCP](https://github.com/shenwei356/kmcp) | 新增记录 | Source only; no adapter |  |
| 162 | [Kraken2](https://github.com/DerrickWood/kraken2) | 合并已有 | kraken2-classify; 1 related workflow guides |  |
| 163 | [MEGAHIT](https://github.com/voutcn/megahit) | 新增记录 | Source only; no adapter |  |
| 164 | [MetaPhlAn / StrainPhlAn](https://github.com/biobakery/MetaPhlAn) | 合并已有 | 3 related workflow guides | StrainPhlAn为同项目组件，不重复计数 |
| 165 | [metaWRAP](https://github.com/bxlab/metaWRAP) | 新增记录 | Source only; no adapter |  |
| 166 | [mlst](https://github.com/tseemann/mlst) | 新增记录 | mlst-typing |  |
| 167 | [mOTUs](https://github.com/motu-tool/mOTUs) | 新增记录 | Source only; no adapter |  |
| 168 | [Panaroo](https://github.com/gtonkinhill/panaroo) | 新增记录 | Source only; no adapter |  |
| 169 | [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) | 新增记录 | Source only; no adapter |  |
| 170 | [QIIME2 / Rachis framework](https://github.com/rachis-org/rachis) | 新增记录 | Source only; no adapter | 原qiime2/qiime2框架仓库已更名Rachis；QIIME2扩增子分发仍独立存在 |
| 171 | [ResFinder](https://github.com/genomicepidemiology/resfinder) | 新增记录 | 1 related workflow guides | 需另外下载相应参考库 |
| 172 | [RGI](https://github.com/arpcard/rgi) | 新增记录 | Source only; no adapter | CARD数据库许可与软件分别核对 |
| 173 | [SemiBin](https://github.com/BigDataBiology/SemiBin) | 新增记录 | Source only; no adapter |  |
| 174 | [Snippy](https://github.com/tseemann/snippy) | 新增记录 | Source only; no adapter |  |
| 175 | [SPAdes / metaSPAdes](https://github.com/ablab/spades) | 合并已有 | dnakit.bioinformatics.ngs.assemble_reads; 1 related workflow guides | 同一套件组件合并计数; Existing SPAdes short-read/hybrid wrapper; not every metaSPAdes mode. |
| 176 | [VAMB](https://github.com/RasmussenLab/vamb) | 新增记录 | 1 related workflow guides |  |
| 177 | [VSEARCH](https://github.com/torognes/vsearch) | 新增记录 | vsearch-dereplicate; vsearch-cluster |  |
| 178 | [Cello v2](https://github.com/CIDARLAB/Cello-v2) | 新增记录 | Source only; no adapter | 安装与适用宿主按原文档核对 |
| 179 | [DNA Cauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron) | 新增记录 | Source only; no adapter |  |
| 180 | [DNA Chisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) | 合并已有 | 1 related workflow guides |  |
| 181 | [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) | 新增记录 | Source only; no adapter |  |
| 182 | [DNA Fountain](https://github.com/TeamErlich/dna-fountain) | 新增记录 | Source only; no adapter | 原论文代码使用旧Python；不应直接装入主环境 |
| 183 | [DNA Weaver](https://github.com/Edinburgh-Genome-Foundry/DnaWeaver) | 合并已有 | Source only; no adapter |  |
| 184 | [DNAplotlib](https://github.com/VoigtLab/dnaplotlib) | 新增记录 | Source only; no adapter | README采用PYTHONPATH导入；不要假定包名即PyPI包 |
| 185 | [oxDNA](https://github.com/lorenzo-rovigatti/oxDNA) | 合并已有 | Source only; no adapter | 包括Python绑定；编译和运行环境独立 |
| 186 | [Primer3](https://github.com/primer3-org/primer3) | 新增记录 | Source only; no adapter |  |
| 187 | [primer3-py](https://github.com/libnano/primer3-py) | 新增记录 | Source only; no adapter | 不是独立引物算法；需引用底层Primer3 |
| 188 | [pydna](https://github.com/pydna-group/pydna) | 新增记录 | dnakit.bioinformatics.cloning.assemble_duplexes; dnakit.bioinformatics.cloning.read_constructs |  |
| 189 | [pySBOL3](https://github.com/SynBioDex/pySBOL3) | 合并已有 | Source only; no adapter |  |
| 190 | [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano) | 合并已有 | dnakit.genomics.design_dna_duplex; dnakit.genomics.export_oxdna_design; 1 related workflow guides | 属于DNA纳米技术扩展方向 |
| 191 | [SynBioHub](https://github.com/SynBioHub/synbiohub) | 新增记录 | Source only; no adapter |  |
| 192 | [EuroForMix](https://euroformix.com/) | 合并已有 | 1 related workflow guides | 需实验室验证、正确群体频率与专家复核 |
| 193 | [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | 新增记录 | Source only; no adapter | 当前仓库标注PolyForm Strict 1.0.0；不是无条件可商用开源 |
| 194 | [FDSTools](https://github.com/Jerrythafast/FDSTools) | 合并已有 | 1 related workflow guides | 计算结果不能自动生成司法鉴定结论 |
| 195 | [forrel](https://github.com/magnusdv/forrel) | 合并已有 | 1 related workflow guides | pedsuite生态；需要相应标记频率与家系假设，结果不是自动司法鉴定意见 |
| 196 | [mapDamage](https://github.com/ginolhac/mapDamage) | 合并已有 | 1 related workflow guides |  |
| 197 | [MToolBox](https://github.com/mitoNGS/MToolBox) | 新增记录 | Source only; no adapter | 用于谱系及医学辅助；不能单独判定身份 |
| 198 | [schmutzi](https://github.com/grenaud/schmutzi) | 合并已有 | 1 related workflow guides |  |
| 199 | [Straglr](https://github.com/BirolLab/straglr) | 新增记录 | Source only; no adapter | 旧bcgsc地址已迁移 |
| 200 | [STRling](https://github.com/quinlan-lab/STRling) | 新增记录 | Source only; no adapter |  |
| 201 | [TRTools](https://github.com/gymrek-lab/TRTools) | 新增记录 | Source only; no adapter |  |
| 202 | [Bismark](https://github.com/FelixKrueger/Bismark) | 合并已有 | 2 related workflow guides |  |
| 203 | [bsseq](https://github.com/hansenlab/bsseq) | 新增记录 | Source only; no adapter |  |
| 204 | [bwa-meth](https://github.com/brentp/bwa-meth) | 新增记录 | Source only; no adapter | 作者建议提取甲基化使用MethylDackel |
| 205 | [ccsmeth](https://github.com/PengNi/ccsmeth) | 新增记录 | Source only; no adapter |  |
| 206 | [ChAMP](https://bioconductor.org/packages/release/bioc/html/ChAMP.html) | 新增记录 | Source only; no adapter |  |
| 207 | [DeepSignal](https://github.com/bioinfomaticsCSU/deepsignal) | 新增记录 | Source only; no adapter | 此原版依赖旧TensorFlow与Tombo；应与新版区分 |
| 208 | [DMRcate](https://bioconductor.org/packages/release/bioc/html/DMRcate.html) | 新增记录 | Source only; no adapter |  |
| 209 | [dmrseq](https://bioconductor.org/packages/release/bioc/html/dmrseq.html) | 新增记录 | Source only; no adapter |  |
| 210 | [Dorado](https://github.com/nanoporetech/dorado) | 新增记录 | Source only; no adapter | 匹配测序化学版本与模型 |
| 211 | [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html) | 新增记录 | Source only; no adapter |  |
| 212 | [EpiDISH](https://bioconductor.org/packages/release/bioc/html/EpiDISH.html) | 新增记录 | Source only; no adapter |  |
| 213 | [f5c](https://github.com/hasindu2008/f5c) | 新增记录 | Source only; no adapter | 重实现部分Nanopolish模块 |
| 214 | [Megalodon](https://github.com/nanoporetech/megalodon) | 新增记录 | Source only; no adapter | 旧生态兼容性应专门核对 |
| 215 | [methplotlib](https://github.com/wdecoster/methplotlib) | 新增记录 | Source only; no adapter |  |
| 216 | [MethylDackel](https://github.com/dpryan79/MethylDackel) | 新增记录 | Source only; no adapter |  |
| 217 | [methylKit](https://github.com/al2na/methylKit) | 合并已有 | 2 related workflow guides |  |
| 218 | [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html) | 新增记录 | Source only; no adapter |  |
| 219 | [Modkit](https://github.com/nanoporetech/modkit) | 合并已有 | dnakit.genomics.read_bedmethyl; dnakit.genomics.summarize_methylation; modkit-pileup; 2 related workflow guides |  |
| 220 | [NanoMethViz](https://github.com/Shians/NanoMethViz) | 新增记录 | Source only; no adapter |  |
| 221 | [Nanopolish](https://github.com/jts/nanopolish) | 新增记录 | Source only; no adapter |  |
| 222 | [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) | 新增记录 | Source only; no adapter | 官方预编译包；不是普通pip包 |
| 223 | [Remora](https://github.com/nanoporetech/remora) | 新增记录 | Source only; no adapter |  |
| 224 | [SeSAMe](https://bioconductor.org/packages/release/bioc/html/sesame.html) | 新增记录 | Source only; no adapter |  |
| 225 | [Tombo](https://github.com/nanoporetech/tombo) | 新增记录 | Source only; no adapter | 旧模型依赖链；不建议盲目升级安装 |
| 226 | [ArchR](https://github.com/GreenleafLab/ArchR) | 新增记录 | Source only; no adapter |  |
| 227 | [BEDTools](https://github.com/arq5x/bedtools2) | 新增记录 | 1 related workflow guides |  |
| 228 | [bioframe](https://github.com/open2c/bioframe) | 新增记录 | Source only; no adapter |  |
| 229 | [Biopython](https://github.com/biopython/biopython) | 新增记录 | dnakit.bioinformatics.reads.read_fastx; dnakit.bioinformatics.alignment.read_alignment; dnakit.bioinformatics.motifs.build_motif |  |
| 230 | [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | 新增记录 | Source only; no adapter |  |
| 231 | [bioutils](https://github.com/biocommons/bioutils) | 新增记录 | Source only; no adapter |  |
| 232 | [Bowtie2](https://github.com/BenLangmead/bowtie2) | 新增记录 | dnakit.bioinformatics.ngs.map_reads | ; map_reads(method='bowtie2'); bowtie2 and bowtie2-build required. |
| 233 | [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) | 新增记录 | Source only; no adapter |  |
| 234 | [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) | 新增记录 | Source only; no adapter |  |
| 235 | [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) | 合并已有 | 2 related workflow guides |  |
| 236 | [Canu](https://github.com/marbl/canu) | 新增记录 | Source only; no adapter |  |
| 237 | [Cooler](https://github.com/open2c/cooler) | 合并已有 | 2 related workflow guides |  |
| 238 | [cooltools](https://github.com/open2c/cooltools) | 合并已有 | 2 related workflow guides |  |
| 239 | [Cutadapt](https://github.com/marcelm/cutadapt) | 新增记录 | Source only; no adapter |  |
| 240 | [cyvcf2](https://github.com/brentp/cyvcf2) | 新增记录 | Source only; no adapter |  |
| 241 | [deepTools](https://github.com/deeptools/deepTools) | 新增记录 | 1 related workflow guides |  |
| 242 | [fastp](https://github.com/OpenGene/fastp) | 合并已有 | 1 related workflow guides |  |
| 243 | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 新增记录 | dnakit.bioinformatics.ngs.quality_control |  |
| 244 | [Flye](https://github.com/mikolmogorov/Flye) | 合并已有 | flye-assembly; 2 related workflow guides | 旧fenderglass地址已迁移 |
| 245 | [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) | 新增记录 | Source only; no adapter |  |
| 246 | [gffutils](https://github.com/daler/gffutils) | 新增记录 | 1 related workflow guides |  |
| 247 | [gget](https://github.com/scverse/gget) | 新增记录 | Source only; no adapter | 由pachterlab迁移至scverse |
| 248 | [Hail](https://hail.is/) | 新增记录 | Source only; no adapter | 需要独立计算环境；安装包不等于配置完毕 |
| 249 | [hgvs](https://github.com/biocommons/hgvs) | 新增记录 | Source only; no adapter |  |
| 250 | [HiCExplorer](https://github.com/deeptools/HiCExplorer) | 新增记录 | Source only; no adapter |  |
| 251 | [hifiasm](https://github.com/chhylp123/hifiasm) | 合并已有 | hifiasm-assembly; 2 related workflow guides |  |
| 252 | [HTSlib](https://github.com/samtools/htslib) | 新增记录 | Source only; no adapter | 供pysam、BCFtools等使用 |
| 253 | [IGV](https://github.com/igvteam/igv) | 新增记录 | Source only; no adapter |  |
| 254 | [MACS3](https://github.com/macs3-project/MACS) | 新增记录 | 1 related workflow guides |  |
| 255 | [minigraph](https://github.com/lh3/minigraph) | 合并已有 | minigraph-pangenome |  |
| 256 | [Minimap2](https://github.com/lh3/minimap2) | 新增记录 | dnakit.annotation.genes.annotate_genes | ；官方另提供mappy Python绑定，需独立核对版本; Existing reference-based gene annotation/alignment wrapper; not every minimap2 CLI mode. |
| 257 | [MMseqs2](https://github.com/soedinglab/MMseqs2) | 新增记录 | mmseqs-protein-search; mmseqs-protein-cluster |  |
| 258 | [MultiQC](https://github.com/MultiQC/MultiQC) | 新增记录 | Source only; no adapter |  |
| 259 | [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) | 新增记录 | dnakit.bioinformatics.alignment.build_blast_database; dnakit.bioinformatics.alignment.local_blast |  |
| 260 | [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) | 新增记录 | Source only; no adapter | CLI与API是数据获取入口，不是新增预测算法 |
| 261 | [ODGI](https://github.com/pangenome/odgi) | 新增记录 | Source only; no adapter |  |
| 262 | [pairtools](https://github.com/open2c/pairtools) | 新增记录 | Source only; no adapter |  |
| 263 | [PGGB](https://github.com/pangenome/pggb) | 合并已有 | 2 related workflow guides |  |
| 264 | [pybedtools](https://github.com/daler/pybedtools) | 新增记录 | Source only; no adapter | 多数操作需要外部bedtools可执行程序 |
| 265 | [pyBigWig](https://github.com/deeptools/pyBigWig) | 新增记录 | 1 related workflow guides |  |
| 266 | [pyfaidx](https://github.com/mdshw5/pyfaidx) | 新增记录 | Source only; no adapter |  |
| 267 | [pyfastx](https://github.com/lmdu/pyfastx) | 新增记录 | Source only; no adapter |  |
| 268 | [PyRanges0](https://github.com/pyranges/pyranges0) | 新增记录 | Source only; no adapter | 旧pyranges地址已重定向；新旧API需区分 |
| 269 | [pysam](https://github.com/pysam-developers/pysam) | 新增记录 | dnakit.bioinformatics.reads.alignment_summary; dnakit.bioinformatics.reads.coverage_consensus |  |
| 270 | [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) | 新增记录 | Source only; no adapter |  |
| 271 | [SAMtools](https://github.com/samtools/samtools) | 新增记录 | dnakit.bioinformatics.ngs.prepare_bam; 2 related workflow guides |  |
| 272 | [SeqKit](https://github.com/shenwei356/seqkit) | 合并已有 | 1 related workflow guides |  |
| 273 | [seqtk](https://github.com/lh3/seqtk) | 新增记录 | Source only; no adapter |  |
| 274 | [Signac](https://github.com/stuart-lab/signac) | 新增记录 | Source only; no adapter |  |
| 275 | [Trim Galore](https://github.com/FelixKrueger/TrimGalore) | 新增记录 | Source only; no adapter |  |
| 276 | [Unicycler](https://github.com/rrwick/Unicycler) | 新增记录 | Source only; no adapter |  |
| 277 | [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) | 新增记录 | Source only; no adapter |  |
| 278 | [vg](https://github.com/vgteam/vg) | 新增记录 | Source only; no adapter |  |
| 279 | [Winnowmap](https://github.com/marbl/Winnowmap) | 新增记录 | Source only; no adapter |  |
| 280 | [AlphaGenome API](https://github.com/google-deepmind/alphagenome) | 新增记录 | Source only; no adapter | API路径涉及数据外发；与本地研究代码分开 |
| 281 | [AlphaGenome Research](https://github.com/google-deepmind/alphagenome_research) | 新增记录 | dnakit.predictions.api.predict_properties | 代码Apache2；权重需接受非商业条款；要求较高的GPU资源; Existing AlphaGenome research-weight backend; access and non-commercial model terms apply independently. |
| 282 | [Caduceus](https://github.com/kuleshov-group/caduceus) | 新增记录 | 1 related workflow guides |  |
| 283 | [DNABERT2](https://github.com/MAGICS-LAB/DNABERT_2) | 新增记录 | 1 related workflow guides |  |
| 284 | [Evo2](https://github.com/ArcInstitute/evo2) | 新增记录 | dnakit.predictions.api.predict_properties; 1 related workflow guides | 按模型卡限定任务与资源；不将生成序列视为已验证生物功能; Existing configured Evo2 prediction backend; local checkpoints required. |
| 285 | [HyenaDNA](https://github.com/HazyResearch/hyena-dna) | 新增记录 | Source only; no adapter |  |
| 286 | [Nucleotide Transformer / AgroNT / SegmentNT](https://github.com/instadeepai/nucleotide-transformer) | 新增记录 | dnakit.predictions.api.predict_sequence_properties; 1 related workflow guides | 一个项目的模型家族合并计数；不同权重适用物种和任务不同; Existing SegmentNT prediction backend only; not all NT/AgroNT models. |
| 287 | [Anthropic Life Sciences](https://github.com/anthropics/life-sciences) | 新增记录 | pubmed-anthropic; opentargets-platform; 10x-genomics; 3 bundled Skills | 含nextflow-development；合作服务不必然可本地开源部署 |
| 288 | [AutoBA](https://github.com/Panomica/AutoBA) | 新增记录 | Source only; no adapter | GitHub搜索确认仓库；README本次未成功读取，作者关联与完整运行需复核 |
| 289 | [AWS Kiro Life Sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences) | 新增记录 | kiro-genomics; kiro-model-organisms; kiro-epigenomics; kiro-agriculture; kiro-ecology; kiro-clinical; kiro-pipelines; 1 bundled Skills | 外部云服务可能需账户与费用 |
| 290 | [BioCoder](https://github.com/gersteinlab/BioCoder) | 新增记录 | Source only; no adapter | 评测数据与抽取方法，不算新增分析功能 |
| 291 | [BioContextAI Knowledgebase MCP](https://github.com/biocontext-ai/knowledgebase-mcp) | 新增记录 | biocontext-kb | 含Ensembl和文献接口；非所有后端无需认证 |
| 292 | [BioContextAI Registry](https://github.com/biocontext-ai/registry) | 新增记录 | Source only; no adapter | 这是注册目录，不是一个DNA分析算法 |
| 293 | [BioinfoMCP](https://github.com/florensiawidjaja/BioinfoMCP) | 合并已有 | Source only; no adapter | 论文报告的测试不等同在用户环境测试 |
| 294 | [BioMCP (GenomOncology)](https://github.com/genomoncology/biomcp) | 合并已有 | biomcp; 1 bundled Skills; 2 related workflow guides | 正确PyPI发行包为biomcp-cli；biomcp为无关同名包 |
| 295 | [BioThings MCP](https://github.com/longevity-genie/biothings-mcp) | 新增记录 | biothings |  |
| 296 | [BOLE](https://bole.zishuailab.com/) | 合并已有 | 3 related workflow guides | 机构新闻确认发布；站点和完整安装入口本次未核实。GitHub ChongxiaoLiang/BOLE的README仅说明基准证据与日志，不能据此认定完整代码已获得 |
| 297 | [BrAPI MCP](https://github.com/cyanheads/brapi-mcp-server) | 新增记录 | brapi | 社区实现；需连接BrAPI后端，不能替代基因组选择算法 |
| 298 | [Chorus](https://github.com/pinellolab/chorus) | 合并已有 | chorus | 作者明确不在PyPI；pip install chorus会装到无关同名包；应固定源码版本 |
| 299 | [Clair Skills](https://github.com/HKU-BAL/Clair-skills) | 合并已有 | 1 related workflow guides | 是操作指导层，不是新变异检测算法 |
| 300 | [ClawBio](https://github.com/ClawBio/ClawBio) | 合并已有 | 23 bundled Skills; 21 related workflow guides | catalog中存在planned/mvp且未benchmark_validated的项目；不能把所有条目视为生产可用 |
| 301 | [CRISPR-GPT public light](https://github.com/cong-lab/crispr-gpt-pub) | 新增记录 | Source only; no adapter | README明确完整版未公开；不能按完整版功能接入 |
| 302 | [EBI Tools MCP (QuentinCody)](https://github.com/QuentinCody/ebi-tools-mcp-server) | 新增记录 | ebi-tools | 社区候选，不标为EBI官方 |
| 303 | [ENCODE Toolkit](https://github.com/ammawla/encode-toolkit) | 合并已有 | encode-toolkit; 16 bundled Skills; 3 related workflow guides | AGPL3；Zenodo条目是软件存档，不冒充同行评议论文 |
| 304 | [Ensembl MCP (cyanheads)](https://github.com/cyanheads/ensembl-mcp-server) | 合并已有 | ensembl | 社区封装，不标为Ensembl官方发布 |
| 305 | [Galaxy Claude Skills (nekrut)](https://github.com/nekrut/galaxy-claude-skills) | 新增记录 | 1 bundled Skills | 搜索发现；本次未成功读取README |
| 306 | [GeneGPT](https://github.com/ncbi/GeneGPT) | 新增记录 | Source only; no adapter | 旧API/模型接口需适配；研究用途 |
| 307 | [GenePattern MCP](https://github.com/genepattern/genepattern-mcp) | 合并已有 | genepattern | 需GenePattern后端及相应模块 |
| 308 | [gget MCP](https://github.com/longevity-genie/gget-mcp) | 新增记录 | gget | 底层gget单列；接口数不等于独立算法数 |
| 309 | [Google DeepMind Science Skills](https://github.com/google-deepmind/science-skills) | 合并已有 | 9 bundled Skills | 只选择DNA相关子集；不同Skill与第三方数据分别授权 |
| 310 | [GPTomics bioSkills](https://github.com/GPTomics/bioSkills) | 合并已有 | 331 bundled Skills; 3 related workflow guides | 仓库于2026-08-15归档；保留为可复现参考，非持续维护承诺 |
| 311 | [GrEBI](https://github.com/EBISPOT/GrEBI) | 新增记录 | grebi | 项目提供开发服务器；不等同生产服务保证 |
| 312 | [GWAS MCP (muslus)](https://github.com/muslus/gwas-mcp) | 新增记录 | gwas-mcp | 搜索发现；本次未成功读取README |
| 313 | [K-Dense Scientific Agent Skills](https://github.com/K-Dense-AI/scientific-agent-skills) | 合并已有 | 9 bundled Skills; 1 related workflow guides | 原claude-scientific-skills名称已迁移；不把全科学Skill总数当DNA功能数 |
| 314 | [Medea](https://github.com/mims-harvard/Medea) | 新增记录 | Source only; no adapter | DNA相关性限基因/患者数据部分；研究输出非诊疗建议 |
| 315 | [NCBI Database MCP (hpend2373)](https://github.com/hpend2373/NCBI-Database-MCP) | 新增记录 | ncbi-database | 搜索发现；本次未成功读取README |
| 316 | [Open Targets MCP](https://github.com/nickzren/opentargets-mcp) | 新增记录 | opentargets | 社区封装；查询不等于临床决策 |
| 317 | [Operon](https://github.com/swaruplab/operon) | 新增记录 | Source only; no adapter | 工作台、协议与分析工具分别验收 |
| 318 | [Plant Genomics MCP](https://github.com/musharna/plant-genomics-mcp) | 新增记录 | plant-genomics | 社区项目；批量/单条接口不重复当独立算法 |
| 319 | [Precision Medicine MCP](https://github.com/lynnlangit/precision-medicine-mcp) | 新增记录 | Source only; no adapter | 示例/研究假说不代表临床验证；必须专业审核 |
| 320 | [PubMed MCP (cyanheads)](https://github.com/cyanheads/pubmed-mcp-server) | 新增记录 | pubmed | 文献基础设施，不是DNA算法 |
| 321 | [SciAgent Skills](https://github.com/jaechang-hits/SciAgent-Skills) | 合并已有 | 21 bundled Skills | DNA任务是其中子集；非全库均属于DNA |
| 322 | [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) | 合并已有 | tooluniverse; 35 bundled Skills | 跨学科聚合器；不能把其所有工具算DNA工具；有tooluniverse-doctor |
| 323 | [TxAgent](https://github.com/mims-harvard/TxAgent) | 新增记录 | Source only; no adapter | 属于基因/疾病证据辅助，不是DNA序列分析器 |
| 324 | [UCSC Genome MCP](https://github.com/hlydecker/ucsc-genome-mcp) | 新增记录 | ucsc-genome | 社区封装 |
| 325 | [Galaxy](https://github.com/galaxyproject/galaxy) | 新增记录 | Source only; no adapter | 平台和平台内工具分别计数 |
| 326 | [Nextflow](https://github.com/nextflow-io/nextflow) | 新增记录 | Source only; no adapter | 本身不是DNA算法 |
| 327 | [nf-core/atacseq](https://nf-co.re/atacseq) | 新增记录 | Source only; no adapter |  |
| 328 | [nf-core/bacass](https://nf-co.re/bacass) | 新增记录 | Source only; no adapter |  |
| 329 | [nf-core/chipseq](https://nf-co.re/chipseq) | 新增记录 | Source only; no adapter |  |
| 330 | [nf-core/eager](https://nf-co.re/eager) | 合并已有 | 1 related workflow guides |  |
| 331 | [nf-core/mag](https://nf-co.re/mag) | 新增记录 | Source only; no adapter |  |
| 332 | [nf-core/methylseq](https://nf-co.re/methylseq) | 合并已有 | 1 related workflow guides |  |
| 333 | [nf-core/sarek](https://nf-co.re/sarek) | 合并已有 | 2 related workflow guides |  |
| 334 | [nf-core/taxprofiler](https://nf-co.re/taxprofiler) | 合并已有 | 2 related workflow guides |  |
| 335 | [Snakemake](https://github.com/snakemake/snakemake) | 新增记录 | Source only; no adapter | 本身不是DNA算法 |

## 3) Citations, wrappers and dependencies

The 86 supplied relationships retain their distinct types. A retrieved page or a matching name alone does not verify a relationship. Context-reviewed relationships are marked separately; other claims remain supplied evidence awaiting review.

| Source / 来源 | Relationship / 关系 | Target / 对象 | Evidence / 证据 | Verification / 核验 |
| --- | --- | --- | --- | --- |
| AlphaSimR | `documentation_links_publication` | [https://doi.org/10.1093/g3journal/jkaa017](https://doi.org/10.1093/g3journal/jkaa017) | [source](https://gaynorr.github.io/AlphaSimR/) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| GAPIT | `documentation_links_publication` | [https://doi.org/10.1093/molbev/msag107](https://doi.org/10.1093/molbev/msag107) | [source](https://github.com/jiabowang/GAPIT) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| rMVP | `documentation_links_publication` | [https://doi.org/10.1016/j.gpb.2020.10.007](https://doi.org/10.1016/j.gpb.2020.10.007) | [source](https://github.com/xiaolei-lab/rMVP) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Borzoi | `documentation_links_publication` | [https://www.biorxiv.org/content/10.1101/2023.08.30.555582v1](https://www.biorxiv.org/content/10.1101/2023.08.30.555582v1) | [source](https://github.com/calico/borzoi) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| ChromBPNet | `documentation_links_publication` | [https://doi.org/10.1101/2024.12.25.630221](https://doi.org/10.1101/2024.12.25.630221) | [source](https://github.com/kundajelab/chrombpnet) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Orca | `documentation_links_publication` | [https://doi.org/10.1038/s41588-022-01065-4](https://doi.org/10.1038/s41588-022-01065-4) | [source](https://github.com/jzhoulab/orca) | Context verified; README 手稿已发表链接 https://www.nature.com/articles/s41588-022-01065-4 与输入 DOI 对应。 |
| AUGUSTUS | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btn013](https://doi.org/10.1093/bioinformatics/btn013) | [source](https://github.com/Gaius-Augustus/Augustus) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta | `documentation_links_publication` | [https://doi.org/10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta documentation | `documentation_links_publication` | [https://doi.org/10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta documentation | `documentation_links_publication` | [https://doi.org/10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Prokka | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153) | [source](https://github.com/tseemann/prokka) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta documentation | `documentation_links_publication` | [https://doi.org/10.21105/joss.04296](https://doi.org/10.21105/joss.04296) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| CADD | `documentation_links_publication` | [https://doi.org/10.1093/nar/gkad989](https://doi.org/10.1093/nar/gkad989) | [source](https://github.com/kircherlab/CADD-scripts) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| cuteSV | `documentation_links_publication` | [https://doi.org/10.1186/s13059-020-02107-y](https://doi.org/10.1186/s13059-020-02107-y) | [source](https://github.com/tjiangHIT/cuteSV) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| DeepVariant | `documentation_links_publication` | [https://doi.org/10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235) | [source](https://github.com/google/deepvariant) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| DELLY | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378) | [source](https://github.com/dellytools/delly) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| coloc | `documentation_links_publication` | [https://doi.org/10.1371/journal.pgen.1011697](https://doi.org/10.1371/journal.pgen.1011697) | [source](https://github.com/chr1swallace/coloc) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| IQ-TREE 2 | `documentation_links_publication` | [https://doi.org/10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015) | [source](https://github.com/iqtree/iqtree2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| susieR | `documentation_links_publication` | [https://doi.org/10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388) | [source](https://github.com/stephenslab/susieR) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta documentation | `documentation_links_publication` | [https://doi.org/10.1038/s41598-021-91456-0](https://doi.org/10.1038/s41598-021-91456-0) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bracken | `documentation_links_publication` | [https://peerj.com/articles/cs-104/](https://peerj.com/articles/cs-104/) | [source](https://github.com/jenniferlu717/Bracken) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| DAS Tool | `documentation_links_publication` | [https://doi.org/10.1038/s41564-018-0171-1](https://doi.org/10.1038/s41564-018-0171-1) | [source](https://github.com/cmks/DAS_Tool) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| KMCP | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btac845](https://doi.org/10.1093/bioinformatics/btac845) | [source](https://github.com/shenwei356/kmcp) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Kraken2 | `documentation_links_publication` | [https://doi.org/10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0) | [source](https://github.com/DerrickWood/kraken2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| DNA Chisel | `documentation_links_publication` | [https://academic.oup.com/bioinformatics/article/36/16/4508/5869515](https://academic.oup.com/bioinformatics/article/36/16/4508/5869515) | [source](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| STRling | `documentation_links_publication` | [https://doi.org/10.1186/s13059-022-02826-4](https://doi.org/10.1186/s13059-022-02826-4) | [source](https://github.com/quinlan-lab/STRling) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| TRTools | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btaa736](https://doi.org/10.1093/bioinformatics/btaa736) | [source](https://github.com/gymrek-lab/TRTools) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Biopython | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) | [source](https://github.com/biopython/biopython) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| HiCExplorer | `documentation_links_publication` | [https://doi.org/10.1093/gigascience/giac061](https://doi.org/10.1093/gigascience/giac061) | [source](https://github.com/deeptools/HiCExplorer) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Minimap2 | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191) | [source](https://github.com/lh3/minimap2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| MMseqs2 | `documentation_links_publication` | [https://doi.org/10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988) | [source](https://github.com/soedinglab/MMseqs2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Bakta documentation | `documentation_links_publication` | [https://doi.org/10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421) | [source](https://github.com/oschwengers/bakta) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| SeqKit | `documentation_links_publication` | [https://doi.org/10.1002/imt2.191](https://doi.org/10.1002/imt2.191) | [source](https://github.com/shenwei356/seqkit) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| AlphaGenome API; AlphaGenome Research | `documentation_links_publication` | [https://doi.org/10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0) | [source](https://github.com/google-deepmind/alphagenome) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Evo2 | `documentation_links_publication` | [https://doi.org/10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5) | [source](https://github.com/ArcInstitute/evo2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Nucleotide Transformer / AgroNT / SegmentNT | `documentation_links_publication` | [https://www.nature.com/articles/s41592-024-02523-z](https://www.nature.com/articles/s41592-024-02523-z) | [source](https://github.com/instadeepai/nucleotide-transformer) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| BioinfoMCP paper | `paper_cites_paper` | [https://doi.org/10.1002/advs.202407094](https://doi.org/10.1002/advs.202407094) | [source](https://arxiv.org/html/2510.02139v1) | Context verified; 人工核对 BioinfoMCP arXiv HTML 参考文献条目；作者、题名和目标文献对应。 |
| BioCoder | `documentation_links_publication` | [https://arxiv.org/abs/2308.16458](https://arxiv.org/abs/2308.16458) | [source](https://github.com/gersteinlab/BioCoder) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| BioContextAI Registry | `documentation_links_publication` | [https://www.nature.com/articles/s41587-025-02900-9](https://www.nature.com/articles/s41587-025-02900-9) | [source](https://github.com/biocontext-ai/registry) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| BioinfoMCP | `documentation_links_publication` | [https://arxiv.org/abs/2510.02139](https://arxiv.org/abs/2510.02139) | [source](https://github.com/florensiawidjaja/BioinfoMCP) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| CRISPR-GPT public light | `documentation_links_publication` | [https://doi.org/10.1038/s41551-025-01463-z](https://doi.org/10.1038/s41551-025-01463-z) | [source](https://github.com/cong-lab/crispr-gpt-pub) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| ENCODE Toolkit | `software_archive` | [https://doi.org/10.5281/zenodo.18917511](https://doi.org/10.5281/zenodo.18917511) | [source](https://github.com/ammawla/encode-toolkit) | Context verified; README 的 Zenodo 链接指向软件存档；不归为同行评议论文。 |
| GeneGPT | `documentation_links_publication` | [https://academic.oup.com/bioinformatics/article/40/2/btae075/7606338](https://academic.oup.com/bioinformatics/article/40/2/btae075/7606338) | [source](https://github.com/ncbi/GeneGPT) | Context verified; README Citation 的 GeneGPT 题名、Jin/Yang/Chen/Lu、Bioinformatics 40(2):btae075 (2024) 与输入文献一致；README 没有该完整网页 URL。 |
| Nextflow | `documentation_links_publication` | [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) | [source](https://github.com/nextflow-io/nextflow) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| GAPIT | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/bts444](https://doi.org/10.1093/bioinformatics/bts444) | [source](https://github.com/jiabowang/GAPIT) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| GAPIT | `documentation_links_publication` | [https://doi.org/10.1016/j.gpb.2021.08.005](https://doi.org/10.1016/j.gpb.2021.08.005) | [source](https://github.com/jiabowang/GAPIT) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| GAPIT | `documentation_links_publication` | [https://doi.org/10.1093/gigascience/giy154](https://doi.org/10.1093/gigascience/giy154) | [source](https://github.com/jiabowang/GAPIT) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| GAPIT | `documentation_links_publication` | [https://doi.org/10.1371/journal.pgen.1005767](https://doi.org/10.1371/journal.pgen.1005767) | [source](https://github.com/jiabowang/GAPIT) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| AlphaSimR | `documentation_links_publication` | [https://doi.org/10.1101/gr.083634.108](https://doi.org/10.1101/gr.083634.108) | [source](https://gaynorr.github.io/AlphaSimR/) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| CADD | `documentation_links_publication` | [https://doi.org/10.1186/s13073-021-00835-9](https://doi.org/10.1186/s13073-021-00835-9) | [source](https://github.com/kircherlab/CADD-scripts) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| CADD | `documentation_links_publication` | [https://doi.org/10.1038/ng.2892](https://doi.org/10.1038/ng.2892) | [source](https://github.com/kircherlab/CADD-scripts) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| DeepVariant | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btaa1081](https://doi.org/10.1093/bioinformatics/btaa1081) | [source](https://github.com/google/deepvariant) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| cuteSV | `documentation_links_publication` | [https://doi.org/10.1186/s13059-025-03642-2](https://doi.org/10.1186/s13059-025-03642-2) | [source](https://github.com/tjiangHIT/cuteSV) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| susieR | `documentation_links_publication` | [https://doi.org/10.1371/journal.pgen.1010299](https://doi.org/10.1371/journal.pgen.1010299) | [source](https://github.com/stephenslab/susieR) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| coloc | `documentation_links_publication` | [https://doi.org/10.1371/journal.pgen.1009440](https://doi.org/10.1371/journal.pgen.1009440) | [source](https://github.com/chr1swallace/coloc) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| Minimap2 | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btab705](https://doi.org/10.1093/bioinformatics/btab705) | [source](https://github.com/lh3/minimap2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| IQ-TREE 2 | `documentation_links_publication` | [https://doi.org/10.1093/molbev/msx281](https://doi.org/10.1093/molbev/msx281) | [source](https://github.com/iqtree/iqtree2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| IQ-TREE 2 | `documentation_links_publication` | [https://doi.org/10.1038/nmeth.4285](https://doi.org/10.1038/nmeth.4285) | [source](https://github.com/iqtree/iqtree2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| IQ-TREE 2 | `documentation_links_publication` | [https://doi.org/10.1093/bioinformatics/btad540](https://doi.org/10.1093/bioinformatics/btad540) | [source](https://github.com/iqtree/iqtree2) | Context verified; 人工检查原始项目文档中的引用/DOI/论文链接上下文；只证明文档列出该文献，不证明另一论文引用它。 |
| BioinfoMCP paper | `paper_cites_paper` | [https://genome.cshlp.org/content/20/9/1297.short](https://genome.cshlp.org/content/20/9/1297.short) | [source](https://arxiv.org/html/2510.02139v1) | Context verified; 人工核对 BioinfoMCP arXiv HTML 参考文献条目；作者、题名和目标文献对应。 |
| BioinfoMCP paper | `paper_cites_paper` | [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/) | [source](https://arxiv.org/html/2510.02139v1) | Context verified; 人工核对 BioinfoMCP arXiv HTML 参考文献条目；作者、题名和目标文献对应。 |
| gget MCP | `wraps_or_uses_tool` | gget | [source](https://github.com/longevity-genie/gget-mcp) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| BioThings MCP | `wraps_or_uses_tool` | BioThings APIs | [source](https://github.com/longevity-genie/biothings-mcp) | Context verified; README 明确列出 MyGene/MyVariant/MyChem 数据源及 biothings-typed-client GeneTools/VariantTools 封装；没有字面字符串 BioThings APIs。 |
| BrAPI MCP | `wraps_or_uses_tool` | BrAPI | [source](https://github.com/cyanheads/brapi-mcp-server) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| Ensembl MCP (cyanheads) | `wraps_or_uses_tool` | Ensembl REST | [source](https://github.com/cyanheads/ensembl-mcp-server) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| UCSC Genome MCP | `wraps_or_uses_tool` | UCSC Genome Browser API | [source](https://github.com/hlydecker/ucsc-genome-mcp) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| Open Targets MCP | `wraps_or_uses_tool` | Open Targets Platform | [source](https://github.com/nickzren/opentargets-mcp) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| GenePattern MCP | `wraps_or_uses_tool` | GenePattern | [source](https://github.com/genepattern/genepattern-mcp) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| Pyrodigal | `wraps_or_uses_tool` | Prodigal | [source](https://github.com/althonos/pyrodigal) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| crisprDesign | `wraps_or_uses_tool` | crisprScore | [source](https://github.com/crisprVerse/crisprDesign) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| TxAgent | `wraps_or_uses_tool` | ToolUniverse | [source](https://github.com/mims-harvard/TxAgent) | Context verified; 人工核对 README 中对底层工具或 API 的封装说明；没有推断 DNAKit 已运行该科学分析。 |
| Medea | `wraps_or_uses_tool` | ToolUniverse | [source](https://github.com/mims-harvard/Medea) | Context verified; README 说明 Medea 工具可注册到 ToolUniverse；属于可选兼容集成，不能推断所有 Medea 路径都以 ToolUniverse 为必需依赖。 |
| Chorus | `wraps_model` | AlphaGenome | [source](https://github.com/pinellolab/chorus) | Context verified; Chorus README 明确列出提供统一接口的模型；各自权重/运行环境仍需独立配置。 |
| Chorus | `wraps_model` | Enformer | [source](https://github.com/pinellolab/chorus) | Context verified; Chorus README 明确列出提供统一接口的模型；各自权重/运行环境仍需独立配置。 |
| Chorus | `wraps_model` | Borzoi | [source](https://github.com/pinellolab/chorus) | Context verified; Chorus README 明确列出提供统一接口的模型；各自权重/运行环境仍需独立配置。 |
| Chorus | `wraps_model` | ChromBPNet | [source](https://github.com/pinellolab/chorus) | Context verified; Chorus README 明确列出提供统一接口的模型；各自权重/运行环境仍需独立配置。 |
| Chorus | `wraps_model` | Sei | [source](https://github.com/pinellolab/chorus) | Context verified; Chorus README 明确列出提供统一接口的模型；各自权重/运行环境仍需独立配置。 |
| Bakta | `uses_component` | Pyrodigal | [source](https://github.com/oschwengers/bakta) | Context verified; 人工核对 README 的依赖/功能说明或 skills/catalog.json 对应具体组件。 |
| Bakta | `uses_component` | DIAMOND | [source](https://github.com/oschwengers/bakta) | Context verified; 人工核对 README 的依赖/功能说明或 skills/catalog.json 对应具体组件。 |
| Bakta | `uses_component` | AMRFinderPlus | [source](https://github.com/oschwengers/bakta) | Context verified; 人工核对 README 的依赖/功能说明或 skills/catalog.json 对应具体组件。 |
| DELLY | `uses_component` | HTSlib | [source](https://github.com/dellytools/delly) | Context verified; 人工核对 README 的依赖/功能说明或 skills/catalog.json 对应具体组件。 |
| Enformer | `reuses_data_or_schema` | Basenji | [source](https://github.com/google-deepmind/deepmind-research/blob/master/enformer/README.md) | Context verified; 人工核对数据来源/输出格式说明；未误标成必须的运行时软件依赖。 |
| AlphaGenome Research | `reuses_data_or_schema` | Borzoi | [source](https://github.com/google-deepmind/alphagenome_research) | Context verified; 人工核对数据来源/输出格式说明；未误标成必须的运行时软件依赖。 |
| ClawBio analyze-fasta | `uses_component` | Biopython | [source](https://github.com/ClawBio/ClawBio/blob/main/skills/catalog.json) | Context verified; 人工核对 README 的依赖/功能说明或 skills/catalog.json 对应具体组件。 |
| ClawBio bioqc-mcp | `wraps_or_uses_tool` | FastQC / MultiQC | [source](https://github.com/ClawBio/ClawBio/blob/main/skills/catalog.json) | Planned description; implementation unverified; catalog.json 描述 FastQC/MultiQC，但 bioqc-mcp 的 status=planned；只确认计划性描述，未确认实际封装已实现。 |
| Anthropic Life Sciences nextflow-development | `skill_guides_workflow` | nf-core/sarek; nf-core/atacseq | [source](https://github.com/anthropics/life-sciences) | Context verified; nextflow-development 段明确列出 nf-core sarek 和 atacseq；这是指导关系。 |

## 4) Added Skills and deferred components

| Skill | Upstream / 来源 | Scope / 使用范围 |
| --- | --- | --- |
| `bioconductor-bridge` | [user335-clawbio](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/bioconductor-bridge) | Bioconductor package discovery, workflow recommendation, setup inspection, and starter code generation grounded in official Bioconductor containers and BiocManager. |
| `galaxy-bioblend` | [user335-galaxy](https://github.com/nekrut/galaxy-claude-skills/tree/0f833c2f11fa624d9fd87c6d139dac404d6dd4ff) | Use Galaxy REST or the BioBlend Python SDK to manage genomic analysis histories, upload or inspect datasets, discover tools and invoke existing workflows on a user-selected Galaxy server. Includes the original REST and SDK guides; server tools, account access and scientific workflows are external. |
| `galaxy-bridge` | [user335-clawbio](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/galaxy-bridge) | Galaxy tool discovery, intelligent recommendation, and execution — 8,000+ bioinformatics tools from usegalaxy.org with multi-signal scoring and workflow suggestions |
| `instrument-data-to-allotrope` | [user335-anthropic](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/instrument-data-to-allotrope) | Convert laboratory instrument output files (PDF, CSV, Excel, TXT) to Allotrope Simple Model (ASM) JSON format or flattened 2D CSV. Use this skill when scientists need to standardize instrument data for LIMS systems, data lakes, or downstream analysis. Supports auto-detection of instrument types. Outputs include full ASM JSON, flattened CSV for easy import, and exportable Python code for data engineers. Common triggers include converting instrument files, standardizing lab data, preparing data for upload to LIMS/ELN systems, or generating parser code for production pipelines. |
| `interpro-database` | [user335-deepmind](https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0/skills/interpro_database) | Identify domains, families, and sites in proteins; find all proteins in a family or sharing a domain; explore species distribution for a domain; annotate genomes with protein families and GO terms. InterPro combines 14 databases (e.g., Pfam, CDD) into one searchable resource. InterPro-N significantly expands annotation and sequence coverage with deep learning. Includes domain architecture (IDA) search. |
| `kiro-healthomics-genomics` | [user335-aws-kiro](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/kiro-life-sciences) | Prepare DNA sequencing and variant workflows for AWS HealthOmics using selected upstream Kiro guidance. Use for WDL/Nextflow import, genomic file-format checks, container references and workflow/run inspection; requires a separately configured official AWS HealthOmics MCP server or AWS CLI. |
| `nextflow-development` | [user335-anthropic](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/nextflow-development) | Run nf-core bioinformatics pipelines (rnaseq, sarek, atacseq) on sequencing data. Use when analyzing RNA-seq, WGS/WES, or ATAC-seq data—either local FASTQs or public datasets from GEO/SRA. Triggers on nf-core, Nextflow, FASTQ analysis, variant calling, gene expression, differential expression, GEO reanalysis, GSE/GSM/SRR accessions, or samplesheet creation. |
| `onekgpd-individual-queries` | [user335-kdense](https://github.com/K-Dense-AI/scientific-agent-skills/tree/330c8e764435a731eff571e3efdda70b363d0792/skills/onekgpd) | Query the 1000 Genomes Project dataset (3,202 whole-genome-sequenced individuals, GRCh38) at the level of individual participants. Use when a question is about individuals or variants in the 1000 Genomes Project cohort: which individuals carry variants matching specific criteria in a gene or region, which individuals are homozygous-reference at a position, which variants exist in the dataset or carried by specified individuals in a gene or region, the relatedness between two specified individuals. Variants are returned with 1000 Genomes allele frequencies (AF), gnomAD v4.1 exome and genome AF, AlphaMissense score, and HGVSp annotations. |
| `outlier-splicing-detection` | [user335-bioskills](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alternative-splicing/outlier-splicing-detection) | Detects aberrant splicing in single rare-disease patients vs a control panel using FRASER 2.0 (Bioconductor; Beta-binomial autoencoder on Intron Jaccard Index, default delta cutoff 0.1, q hyperparameter), OUTRIDER (gene-level outlier expression via autoencoder denoising), LeafcutterMD (Dirichlet-multinomial outlier mode of LeafCutter for annotation-free junctions), and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression for clinical diagnostics). The statistical model is fundamentally different from differential splicing — single-sample-vs-cohort outlier detection rather than two-group comparison. Standard tool in EU rare-disease (Solve-RD) and NIH UDN programs. Use when applying RNA-seq to undiagnosed Mendelian disease, validating predicted splice variants in clinical samples, or detecting cryptic splicing in disease tissue. |
| `pdb-database` | [user335-deepmind](https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0/skills/pdb_database) | Use when you want to search for or download experimentally-determined 3D structures for biomolecules (proteins, nucleic acids, bound ligands). Supports searching by sequence similarity, structure similarity, chemical and other attributes. Also use to get metadata about biomolecular structure experiments. |
| `scvi-tools` | [user335-anthropic](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/scvi-tools) | Deep learning for single-cell analysis using scvi-tools. This skill should be used when users need (1) data integration and batch correction with scVI/scANVI, (2) ATAC-seq analysis with PeakVI, (3) CITE-seq multi-modal analysis with totalVI, (4) multiome RNA+ATAC analysis with MultiVI, (5) spatial transcriptomics deconvolution with DestVI, (6) label transfer and reference mapping with scANVI/scArches, (7) RNA velocity with veloVI, or (8) any deep learning-based single-cell method. Triggers include mentions of scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI, sysVI, scArches, variational autoencoder, VAE, batch correction, data integration, multi-modal, CITE-seq, multiome, reference mapping, latent space. |
| `single-cell-cnv-inference` | [user335-bioskills](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/cnv-inference) | Infer large-scale copy-number alterations from tumor single-cell or single-nucleus RNA-seq to separate malignant from normal cells and call subclones, using inferCNV, copyKAT, Numbat, and SCEVAN. Use when separating malignant from normal cells in a tumor scRNA-seq dataset, inferring chromosome-arm CNVs or aneuploidy from expression, calling tumor subclones from single cells, choosing a CNV-inference method (reference-based vs reference-free, expression-only vs allele-aware), or deciding which cells are tumor before downstream analysis. |
| `ukb-ppp-region-fetch` | [user335-clawbio](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/ukb-ppp-region-fetch) | Fetch a regional slice of plasma pQTL summary statistics from the UK Biobank Pharma Proteomics Project (UKB-PPP; Sun 2023 Nature) for a specific (protein, ancestry) measurement. Use when an agent needs per-variant beta / SE / p-value around a coloc-lead variant for downstream colocalisation, Mendelian randomisation, or regional plotting against a pQTL exposure. The canonical use case is the cis-window around the protein's coding gene TSS, but UKB-PPP releases full-genome summary stats per protein so any GRCh38 window (including trans loci) is supported when the user supplies an explicit (chromosome, start_bp, end_bp). Input: protein_label (HGNC or UniProt), ancestry, chromosome, start_bp, end_bp. Output: harmonised TSV slice + manifest + human-readable report. |

| Component / 组件 | Decision / 处理原因 |
| --- | --- |
| [Clair Skills](https://github.com/HKU-BAL/Clair-skills) | 技能内容的再分发许可未明确；保留来源，不复制。 |
| [AWS Kiro microbiology](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/life-sciences-microbiology) | QIIME2 路径只返回 queued，未执行分析；未注册该 MCP 组件。 |
| [AWS Kiro molecular biology](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/life-sciences-molbio) | 酶切位置计算、BLAST/Primer3 成功判断存在问题；未注册该 MCP 组件。 |
| [Precision Medicine mcp-fgbio](https://github.com/lynnlangit/precision-medicine-mcp) | 源码返回固定模拟计数、由 hash 生成假 ENSG；未作为真实分析功能接入。 |
| [Synapse marketplace connector](https://github.com/anthropics/life-sciences) | 需要 OAuth；当前 HTTP 桥接未实现该认证流程，保留来源。 |
| [ClawBio MCP](https://github.com/ClawBio/ClawBio) | 上游已弃用并移除通用 MCP；保留已审查 Skills，不补失效入口。 |

上述只列影响本轮接入决策的发现。每个外部分析的 validation、MCP 的 limitations，以及 registry/user-catalog-*-review.json 分别保存未执行条件与源码证据。

## 5) Evidence boundaries and rebuild

The public [BOLE repository](https://github.com/ChongxiaoLiang/BOLE) describes benchmark evidence and conversation logs. A complete runnable framework or public API has not been verified. DNAKit's ridge regression is not presented as a BOLE reproduction.

Skills require explicit licensing, complete resources and pinned revisions before import. Format adaptations retain the original text and attribution. Model code, weights and database terms apply separately. A launch recipe does not imply installation or scientific validation.

The four supplied files remain unchanged. registry/user-catalog-335.json preserves their SHA-256 digests, all 335 project rows, 62 publication/archive links, 86 relationships, 38 coverage mappings and 52 skill-source claims. registry/user-catalog-evidence.json stores retrieval outcomes and content hashes.

```bash
python others/skills/import_user_catalog.py --check
python others/skills/build_genomics_survey.py --check
```
