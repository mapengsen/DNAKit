# 致谢与引用 {#acknowledgements-and-citations}

本页合并 DNAKit 的主要来源、项目引用方式、许可证、第三方使用边界和免责声明，便于集中查阅与追溯。

## 致谢与主要来源 {#acknowledgements}

DNAKit 感谢为本项目提供理论基础、公开参数、格式规范、软件接口和结构数据的论文作者、数据库维护者、标准组织与开源项目贡献者。

本节统一汇总项目使用的方法、论文、数据库和网站。功能页仍会保留与具体参数直接相关的 PMID 或使用边界，便于逐项追溯。列出来源不表示 DNAKit 已获得其数据或软件的再分发许可。

### 方法、论文与格式依据 {#methods-and-references}

| 功能 | DNAKit 类型 | 应引用的主要依据 |
| --- | --- | --- |
| FASTA/FASTQ、GenBank、GFF3、BED、AGP | `reimplementation` / 格式整合 | [NCBI GenBank 样例记录](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)、[Sequence Ontology GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)、[UCSC BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)、[NCBI AGP 2.1](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/)；GenBank 仅支持文档化子集 |
| 坐标、序列对象设计 | `native`，参考既有接口 | Biopython `Seq`/`SeqRecord` 与 scikit-bio `DNA` 仅作设计对照 |
| Shannon entropy | `reimplementation` | C. E. Shannon, *A Mathematical Theory of Communication* (1948), DOI [`10.1002/j.1538-7305.1948.tb01338.x`](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x) |
| linguistic complexity | `reimplementation` | 已公开的 observed/possible k-word product 定义；结果保存精确公式 |
| LZ76 complexity | `reimplementation` | Lempel & Ziv, 1976, DOI [`10.1109/TIT.1976.1055501`](https://doi.org/10.1109/TIT.1976.1055501) |
| 15 组二核苷酸属性 | 用户提供的参数表 | DNAKit 只定义 15×16 JSON schema，不内置 DiProDB 数值。用户若自行取得 DiProDB 或其他来源的表，应记录[DiProDB 论文](https://doi.org/10.1093/nar/gkn597)、实际表页、原始论文、版本和 SHA-256；字段映射见[240 项描述符表](api/features/05_all_descriptors.md) |
| 标准遗传密码表 | `reimplementation` / 标准数据 | [NCBI Standard Genetic Code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) |
| Needleman–Wunsch global alignment | `reimplementation` | Needleman & Wunsch, 1970, DOI [`10.1016/0022-2836(70)90057-4`](https://doi.org/10.1016/0022-2836(70)90057-4) |
| Smith–Waterman local alignment | `reimplementation` | Smith & Waterman, 1981, DOI [`10.1016/0022-2836(81)90087-5`](https://doi.org/10.1016/0022-2836(81)90087-5) |
| Levenshtein distance | `reimplementation` | V. I. Levenshtein, 1966 |
| MinHash | `reimplementation` | A. Z. Broder, 1997, resemblance/containment sketching |
| EvoAug 序列生成 | `reimplementation` | Lee et al., 2023, [EvoAug 论文](https://doi.org/10.1186/s13059-023-02941-w)；[官方 PyTorch 代码](https://github.com/p-koo/evoaug)；DNAKit 采用独立的 `DNASequence` 级实现，不引入 PyTorch 或 Gaussian noise |
| k-mer 保持打乱 | `reimplementation` | [uShuffle 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/)；DNAKit 采用 de Bruijn 多重图上的随机 Euler 路径实现，精确保持指定重叠 k-mer 计数 |
| Dashing Jaccard adapter | `adapter` | Baker & Langmead, 2019, DOI [`10.1186/s13059-019-1875-0`](https://doi.org/10.1186/s13059-019-1875-0)；实际 Dashing 版本、k、exact/HLL 模式和 sketch size |
| DNA 基础模型 rep 与 k-means | `adapter` | 实际运行时引用所选模型及 checkpoint： [DNABERT-2](https://github.com/MAGICS-LAB/DNABERT_2)、[Nucleotide Transformer](https://github.com/instadeepai/nucleotide-transformer)、[HyenaDNA](https://github.com/HazyResearch/hyena-dna)、[Caduceus](https://github.com/kuleshov-group/caduceus)、[GROVER](https://huggingface.co/PoetschLab/GROVER)、[LucaOne](https://github.com/LucaOne/LucaOne)、[GENERator](https://github.com/GenerTeam/GENERator)、[Enformer PyTorch](https://github.com/lucidrains/enformer-pytorch)、[AlphaGenome research](https://github.com/google-deepmind/alphagenome_research)、[JanusDNA](https://github.com/Qihao-Duan/JanusDNA)或[Evo 2](https://github.com/ArcInstitute/evo2)；聚类记录实际 scikit-learn 版本、pooling、PCA、seed 和 checkpoint 来源 |
| Fréchet DNA distance | `reimplementation` + `adapter` | Fréchet 高斯距离沿用 [FCD 原始论文](https://doi.org/10.1021/acs.jcim.8b00234)的数学形式；表征另行引用实际模型和 checkpoint，默认模型引用 [LucaOne](https://doi.org/10.1038/s42256-025-01044-4)。该指标不是 ChemNet FCD |
| DNA Frag / SNN | `reimplementation` + `adapter` | 指标公式参考 [MOSES](https://doi.org/10.3389/fphar.2020.565644)；DNA 适配分别用 fixed-length k-mer 代替 BRICS 片段、用 hashed k-mer 位指纹代替 Morgan 指纹，数值不可与分子指标直接比较 |
| DNA 260 nm 消光系数 | `reimplementation` | Warshaw & Tinoco, 1966, DOI [`10.1016/0022-2836(66)90115-X`](https://doi.org/10.1016/0022-2836(66)90115-X)；Cantor, Warshaw & Shapiro, 1970, DOI [`10.1002/bip.1970.360090909`](https://doi.org/10.1002/bip.1970.360090909)；实际参数集版本 |
| DNA nearest-neighbor/Tm | `reimplementation` | SantaLucia, 1998, DOI [`10.1073/pnas.95.4.1460`](https://doi.org/10.1073/pnas.95.4.1460) |
| Primer3 thermodynamic/design adapter | `adapter` | Untergasser et al., 2012, DOI [`10.1093/nar/gks596`](https://doi.org/10.1093/nar/gks596)；实际 Primer3 CLI 版本、可执行路径和参数目录 |
| NUPACK 二级结构 adapter | `adapter` | Zadeh et al., 2011, DOI [`10.1002/jcc.21596`](https://doi.org/10.1002/jcc.21596)；实际 NUPACK 版本、模型、温度和盐条件 |
| 3DNA/DSSR 参数解析 | `adapter` | 3DNA DOI [`10.1093/nar/gkg680`](https://doi.org/10.1093/nar/gkg680)；DSSR DOI [`10.1093/nar/gkv716`](https://doi.org/10.1093/nar/gkv716)；实际输出版本 |
| PDB 三维结构样本 | 数据记录 | 实际 RCSB PDB ID、实验方法、模型号和坐标文件 SHA-256；当前样本为 1BNA、1AC7、139D |
| 限制酶 | 版本化内置小目录 / 可选数据库 | 实际酶定义；外部使用时引用 REBASE 版本 |
| TF PWM | 调用方提供矩阵 | 实际 motif 数据库（如 JASPAR）版本与条目；DNAKit 不内置活性模型 |
| novelty/memorization | `native` 透明规则 | 实际参考库版本、digest、筛选条件和所选相似度方法 |
| synthesis-risk | `native` 透明规则 | DNAKit 配置与命中规则；不得描述为实验成功率 |

每个运行只需引用真正使用的条目，不应把未调用的后端列为计算依据。

#### 用户二核苷酸表可能引用的原始论文 {#dinucleotide-primary-references}

下表仅帮助用户追溯历史来源，不是 DNAKit 内置数据清单，也不授予复制或再分发这些论文数据的权限。

| 参数 | 原始论文 |
| --- | --- |
| Twist、Tilt、Roll、Shift、Slide、Rise | Perez et al. 2004，[PMID 15562006](https://pubmed.ncbi.nlm.nih.gov/15562006/) |
| Bend、Inclination、Major groove width、Minor groove width | Karas et al. 1996，[PMID 8996793](https://pubmed.ncbi.nlm.nih.gov/8996793/) |
| Direction | Shpigelman et al. 1993，[PMID 8402210](https://pubmed.ncbi.nlm.nih.gov/8402210/) |
| Propeller twist | Gorin et al. 1995，[PMID 7897660](https://pubmed.ncbi.nlm.nih.gov/7897660/) |
| Persistence length | Hogan & Austin 1987，[PMID 3627268](https://pubmed.ncbi.nlm.nih.gov/3627268/) |
| Stacking energy | Sponer et al. 1997，[PMID 9199773](https://pubmed.ncbi.nlm.nih.gov/9199773/) |
| Free energy | Sugimoto et al. 1996，[PMID 8948641](https://pubmed.ncbi.nlm.nih.gov/8948641/) |

### 网站、数据库与公开数据 {#websites-and-data}

| 资源 | 本项目用途 | 网址 |
| --- | --- | --- |
| DiProDB | 用户可自行核对的 DNA 二核苷酸参数及原始论文索引；DNAKit 不含其数值 | [数据库论文](https://doi.org/10.1093/nar/gkn597)、[公共参数表](https://diprodb.fli-leibniz.de/ShowTable.php) |
| NCBI Genetic Codes | 标准遗传密码表 1 | [NCBI Standard Genetic Code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) |
| NCBI Datasets | 分类、组装、基因、病毒查询及数据包 | [REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)、[Genome download](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/)、[Virus package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) |
| NCBI Entrez、BLAST、ClinVar、dbSNP、GEO | accession/项目/样本/变异/表达/文献查询、BLAST 作业和公开文件下载 | [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)、[BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html)、[GEO access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html) |
| Ensembl REST | 坐标序列、转录本、区域/变异/调控、同源和比较基因组 | [Ensembl REST endpoints](https://rest.ensembl.org/) |
| European Nucleotide Archive | Study/Sample/Experiment/Run/Analysis 元数据和公开文件 | [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api) |
| ENCODE Portal | 表观基因组实验、峰、信号和公开文件 | [ENCODE REST API](https://www.encodeproject.org/help/rest-api/) |
| UCSC Genome Browser | 染色体、坐标序列、注释/重复/保守性轨道和下载目录 | [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html) |
| IDT 与 Sigma-Aldrich | OD260/寡核苷酸定量说明和公开 ACGT 算例 | [IDT 定量说明](https://sg.idtdna.com/page/support-and-education/decoded-plus/oligo-quantification-getting-it-right)、[Sigma-Aldrich 参数与算例](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/pcr/quantitation-of-oligos) |
| NUPACK | 条件二级结构后端的分析接口、模型条件及许可边界 | [分析文档](https://docs.nupack.org/analysis/)、[实用函数](https://docs.nupack.org/utilities/)、[模型文档](https://docs.nupack.org/model/)、[下载与许可](https://www.nupack.org/download/overview) |
| 3DNA/DSSR | 三维结构参数定义与 DSSR JSON 解析依据 | [3DNA 参数示意](https://x3dna.org/highlights/schematic-diagrams-of-base-pair-parameters)、[DSSR JSON 文档](https://x3dna.org/highlights/dssr-output-in-json-format) |
| RCSB PDB | 文档验证所用 DNA 三维结构记录 | [1BNA](https://www.rcsb.org/structure/1BNA)、[1AC7](https://www.rcsb.org/structure/1AC7)、[139D](https://www.rcsb.org/structure/139D) |

## 第三方声明 {#third-party-notices}

最后复核：2026-08-23。本页是工程合规清单，不是法律意见。

DNAKit 自身采用 MIT 许可证。第三方包、外部程序、数据库、论文和用户提供的参数表保留各自条款；引用论文不等于获得其软件或数据的再分发许可。

### 分发边界 {#distribution-boundary}

DNAKit 的 wheel 和 sdist 不包含 DiProDB 数值表、Primer3、NUPACK、DSSR/3DNA、Dashing、BLAST、MMseqs2、sourmash、REBASE、JASPAR、DNA 基础模型 checkpoint 或其数据库，也不会自动下载外部程序。存在 CLI adapter 时，用户仍须单独安装并提供显式路径。模型 checkpoint 只在用户显式调用表征 API 时下载到本地 `ckpt/`，不等于随 DNAKit 再分发。

`pyproject.toml` 中的 Python 依赖由安装器作为独立发行包解析，不被复制进 DNAKit wheel。NCBI、Ensembl、ENA、ENCODE 和 UCSC 数据也不随 wheel/sdist 分发；adapter 只在用户显式调用时访问公开接口。直接依赖的主要许可如下；实际发布前仍须检查锁定版本及传递依赖。

| 范围 | 直接依赖 | 许可标记 |
| --- | --- | --- |
| 核心 | PyYAML、Rich、Typer、tomli（Python <3.11） | MIT |
| `viz` | CairoSVG；Pillow | LGPL-3.0-or-later；MIT-CMU |
| `io` | PyArrow | Apache-2.0 |
| `validation` | Biopython | Biopython License Agreement |
| `neural` | Hugging Face Hub、NumPy、scikit-learn、PyTorch、Transformers | Apache-2.0、BSD 系列、Apache-2.0；以实际安装版本为准 |
| `neural-caduceus` | mamba-ssm | Apache-2.0；以实际安装版本为准 |
| `neural-enformer` | enformer-pytorch | MIT；以实际安装版本为准 |
| `neural-evo2` | evo2 及其传递依赖 | 以实际安装版本和官方仓库为准 |
| `docs` | MkDocs；Material for MkDocs；mkdocstrings；nbmake | BSD-2-Clause；MIT；ISC；Apache-2.0 |
| `dev` | build、Hypothesis、mypy、pytest、pytest-cov、Ruff、twine、types-PyYAML | MIT、MPL-2.0、MIT、MIT、MIT、MIT、Apache-2.0、Apache-2.0 |

### 外部工具与用户数据 {#external-tools-and-user-data}

| 资源 | DNAKit 处理方式 | 使用前必须确认 |
| --- | --- | --- |
| Primer3 CLI | 仅通过显式 `primer3_core`、`oligotm`、`ntthal` 路径调用；不含 `primer3-py` 或二进制 | 官方源码标记 GPL-2.0-or-later；核对[仓库](https://github.com/primer3-org/primer3)与实际版本 |
| NUPACK | 用户单独取得并安装；不下载、不打包、不作为在线服务提供 | 独立条款/订阅，可能收费，学术和商业条件不同；见[官方许可页](https://www.nupack.org/download/overview) |
| DSSR/3DNA | 用户单独安装；DNAKit 只解析显式输出 | 学术与商业权限不同；免费 academic Basic 仍受许可约束；见[官方说明](https://home.x3dna.org/highlights/x3dna-dssr-is-funded-and-dssr-basic-academic-is-free) |
| Dashing | 用户提供显式可执行文件；发行包不含源码/二进制 | GPL-3.0；见[官方仓库](https://github.com/dnbaker/dashing) |
| DiProDB 或其他二核苷酸表 | DNAKit 仅提供 JSON schema/loader，不含任何 DiProDB 数值 | 用户须有权使用并记录真实来源、版本和 SHA-256；见[DiProDB](https://diprodb.fli-leibniz.de/) |
| 其他 CLI/数据库 | 不自动安装或打包 | 分别核对程序与数据库条款；可访问不等于可再分发 |
| NCBI/Ensembl/ENA/ENCODE/UCSC 公共接口与数据 | 仅按用户查询下载，保存来源和 checksum；不把远程数据嵌入发行包 | 核对各提供方使用政策、署名、隐私和再分发条款；dbGaP 等受控数据仍需正式授权 |
| Hugging Face DNA 模型 checkpoint | 用户选择模型后下载到 `ckpt/` 并记录来源 manifest；wheel/sdist 不包含权重 | DNABERT-2、NTv2、HyenaDNA、Caduceus、GROVER、LucaOne、GENERator、Enformer 和 Evo 2 各自条款可能不同；逐项核对[序列表征页中的官方链接](api/features/08_fingerprints.md#neural-representations) |
| AlphaGenome | 仅在用户接受访问条款、取得 checkpoint 并安装官方 research 代码后加载 | checkpoint 页面标明受限访问和非商业模型条款；核对[官方仓库](https://github.com/google-deepmind/alphagenome_research)及当前模型页面 |
| JanusDNA | 从 Harvard Dataverse 的官方文件下载并校验 MD5；还需用户提供官方源码环境 | 核对[官方仓库](https://github.com/Qihao-Duan/JanusDNA)、[Dataverse DOI](https://doi.org/10.7910/DVN/HDT0RN)及其代码/数据条款 |

DNABERT-2、NTv2、HyenaDNA、Caduceus、LucaOne 和 GENERator 的 Transformers
checkpoint 包含加载所需 Python 代码。DNAKit 默认拒绝执行；只有用户审查来源并
显式设置 `allow_remote_code=True` 后才启用。该确认不能替代代码审计、依赖锁定
或许可证核对。

进程隔离和“由用户安装”可减少捆绑分发风险，但不能自动证明所有用途合法。商业部署、托管服务、源码/二进制再分发或期刊补充材料归档时，应按实际提交的文件重新审查。

完整英文清单随源码和安装包保存为 `THIRD_PARTY_NOTICES.md`；科学方法、论文和网站见[致谢与主要来源](#acknowledgements)，一般使用风险见[免责声明](#disclaimer)。

## 免责声明 {#disclaimer}

DNAKit 是研究与教学用途的软件，不是医疗器械，也不用于临床诊断、治疗、患者管理、实验安全决策或其他高风险用途。输出可能不完整、不准确或不适用于特定实验；使用者必须独立验证计算、参数、输入数据和实验结论。

DNAKit 按 MIT 许可证“按原样”提供，不作明示或默示担保，并在适用法律允许的最大范围内限制责任。仓库根目录的 `LICENSE` 是正式许可文本；本页只说明使用边界，不替代或修改 MIT 许可证，也不构成法律、医疗或合规意见。

“本项目以学术研究为主要开发目的”不是免责条款。MIT 允许研究和商业使用，因此不能同时把 DNAKit 描述成“仅限非商业使用”，除非未来改用另一份经过审查的许可证。

使用者负责：

- 确认其提供的数据、参数表、数据库和外部程序具有适用的访问、使用与再分发权限；
- 不将 ClinVar、VEP、BLAST 或其他远程提供方记录直接作为临床诊断、来源定论或实验安全结论；
- 遵守第三方许可证、订阅、署名/引用、隐私、伦理、生物安全、出口管制及当地法律；
- 在论文和软件中披露实际方法、版本、参数、限制和来源；
- 在实验、发表、部署或商业使用前进行独立验证及必要的专业审查。

DNAKit 不保证第三方资源免费、持续可用、适合商业用途或允许再分发。条款和费用可能变化，使用时应查看官方页面。具体分发边界见[第三方声明](#third-party-notices)。

## 引用与许可证 {#citations-and-license}

### 引用 DNAKit {#citing-dnakit}

仓库根目录的 `CITATION.cff` 提供机器可读元数据。DNAKit 尚未发布论文或 DOI；当前开发预览版可引用为：

```text
DNAKit contributors. DNAKit 0.1.0.dev0 (development preview), 2026.
```

引用 DNAKit 不能替代对实际算法、参数集、后端和参考数据库的引用。项目使用的方法、论文、数据库和网站统一列在[致谢与主要来源](#acknowledgements)中；结果中的 `Provenance`、`BackendInfo`、`ReferenceLibrary` 和 `RunManifest` 用于保留这些信息。

### 项目许可证 {#project-license}

DNAKit 自身采用 MIT 许可证，版权人为 Pengsen Ma，完整文本见仓库根目录的 `LICENSE`。SPDX 标识符为 `MIT`。安装包同时包含 `THIRD_PARTY_NOTICES.md` 和 `DISCLAIMER.md`；网页版本见[第三方声明](#third-party-notices)和[免责声明](#disclaimer)。

项目仍处于 `0.1.0.dev0` 开发预览阶段，依赖兼容性和外部后端仍需按实际环境复核。

第三方包的许可证不会自动成为 DNAKit 的许可证。

### 实现标签 {#implementation-labels}

| 标签 | 含义 | 当前示例 |
| --- | --- | --- |
| `native` | DNAKit 自身对象、组合逻辑或明确的简单确定性逻辑 | 核心对象、审计、缓存、workflow、scorecard |
| `reimplementation` | 根据公开算法、论文或格式规范独立实现 | alignment、SantaLucia NN、MinHash、格式 codec |
| `adapter` | 调用独立安装的包、CLI 或数据库 | Primer3CLIAdapter、PyArrow Parquet、DashingAdapter、外部 CLI metadata 句柄 |
| `novel` | 经过定义、检索、基线和消融确认的新方法 | 当前没有已确认项目 |

混合/多尺度指纹目前是 `native` 组合，不标记为 `novel`。


### 当前本地依赖许可证快照 {#dependency-license-snapshot}

下表来自 `dnakit-dev` 环境的已安装 distribution metadata 或仓库内第三方 `LICENSE`，用于本地审计，不替代法律意见。

| 组件 | 本地版本/状态 | 许可证标记 | DNAKit 边界 |
| --- | --- | --- | --- |
| PyYAML | 6.0.3 | MIT | 核心配置解析 |
| Rich | 15.0.0 | MIT | CLI/脚本进度显示 |
| Typer | 0.27.1 | MIT | CLI |
| tomli | 2.4.1 | MIT | Python <3.11 条件 TOML 解析 |
| CairoSVG | 2.9.0 | LGPL-3.0-or-later | `viz` extra；SVG 转 PNG/PDF |
| Pillow | 12.3.0 | MIT-CMU | `viz` extra；TIFF/分辨率 metadata |
| PyArrow | 25.0.1 | Apache-2.0 | `io` extra；Parquet 表读写 |
| NumPy | 2.2.6 | BSD-3-Clause | 本地 rep 矩阵和聚类测试 |
| scikit-learn | 1.7.2 | BSD-3-Clause | 本地 k-means/PCA/silhouette 测试 |
| Biopython | 1.88 | `LicenseRef-Biopython-License-Agreement` | 仅 validation extra 的差分对照 |
| MkDocs | 1.6.1 | BSD-2-Clause | 文档开发 |
| MkDocs Material | 9.7.7 | MIT | 文档开发 |
| mkdocstrings | 0.30.1 | ISC | API 文档生成 |
| nbmake | 1.5.5 | Apache-2.0 | Notebook 门禁 |
| Dashing 本地副本 | 项目根下的包外第三方目录 | GPL-3.0（该目录 LICENSE）；运行时 adapter 记录 `GPL-3.0-only` | 从 wheel/sdist 排除；显式 adapter 不自动选择或打包该副本，真实科学差分未完成 |

Primer3 不再是 Python 依赖，也不在任何 extra 中。当前 `dnakit-dev` 尚未安装
完整 PyTorch/Transformers 模型推理栈；GROVER 真实 checkpoint smoke 在独立兼容
模型环境执行，不能把上表解释为 11 种模型均已验证。上表只表示当前直接依赖的
本地 metadata 快照；版本升级、重新打包或改变分发方式后必须重新生成许可清单，
完整直接依赖清单见[第三方声明](#third-party-notices)。

### 已注册外部 CLI 句柄的许可证提示 {#registered-cli-license-notices}

这些值是 `BackendInfo.license_expression` 中用于提醒使用者复核的 metadata，不是 DNAKit 对第三方许可适用性的法律结论。BLAST/MMseqs2/sourmash 的路径发现不执行程序，只有用户显式调用 `.version()` 才运行受限版本命令；Dashing 另有必须显式提供可执行路径并显式调用的科学计算 adapter。

| 后端 ID | DNAKit 记录的许可证提示 | 当前能力 |
| --- | --- | --- |
| `blastn` | `LicenseRef-NCBI` | 被动路径定位、显式版本查询 |
| `mmseqs2` | `GPL-3.0-or-later` | 被动路径定位、显式版本查询 |
| `sourmash` | `BSD-3-Clause` | 被动路径定位、显式版本查询 |
| `dashing` | `GPL-3.0-only` | registry 被动路径/版本句柄；另有显式路径的 exact/HLL Jaccard/Top-k adapter；项目根第三方副本不入发行归档 |

实际使用前必须以用户提供二进制/源码和数据库各自的版本、许可证文件与使用条件为准。

### 受限后端与数据库 {#restricted-backends-and-databases}

| 项目 | 当前处理 |
| --- | --- |
| NUPACK | 提供被动探测和显式 adapter，但不自动安装/下载、不打包、不作为 Web 后端；独立订阅条款可能收费，当前环境不可用且真实数值差分受许可与安装条件约束 |
| Primer3 | 用户单独安装 CLI，并显式提供 `oligotm`、`ntthal` 或 `primer3_core` 路径；DNAKit 不再依赖 `primer3-py`，不按 `PATH` 搜索、不打包；许可证提示为 `GPL-2.0-or-later`，仍须按实际版本复核 |
| DSSR/3DNA | 用户单独取得许可证并安装；DNAKit 只解析用户提供的输出，学术免费条款不能外推为商业免费或任意再分发 |
| Dashing | 本地 GPLv3 源码/二进制不进入发行归档；显式路径 adapter 可做有界 exact/HLL Jaccard 与 Top-k，但尚无真实科学差分且不会自动选择本地副本 |
| BLAST/MMseqs2/sourmash | 不自动安装或下载；已注册被动 metadata/显式版本句柄，但没有搜索、聚类或 sketch 执行器 |
| RepeatMasker/TRF/FIMO | 不自动安装或下载；当前没有统一 adapter |
| REBASE/JASPAR 等数据库 | 不随包再分发；由用户按许可提供并记录版本、日期、筛选和 checksum |
| DiProDB/二核苷酸参数表 | 不内置任何 DiProDB 数值；用户按固定 JSON schema 提供其有权使用的 15×16 表，结果记录表名、版本、来源声明和 SHA-256 |
| DNA 基础模型 checkpoint | 仅在用户显式调用 rep API 时下载到 `ckpt/`，不随包分发；逐模型核对 checkpoint/源码条款，远程代码默认拒绝；当前只完成 GROVER 真实 smoke |

更详细的可用性判定见[条件与不可用功能](planning/04_conditional_and_unavailable_features.md)。
