# DNA 相关 Skills：网络补充检索记录

检索与快照日期：2026-09-12。目标是在原有 DNAKit 技能包上补齐 DNA／基因组学能力，保留具备适用再分发许可的内容，并按任务阶段与输入输出去重。

本轮新增 **143** 个技能，保留原有 393 个，合计 **536** 个。累计锁定 77 个来源仓库，枚举 7114 个 `SKILL.md` 路径；最终从 22 个上游来源收录内容。路径数包含镜像、旧版、示例及非 DNA 内容，不代表独立 DNA 功能数。

这是公开网络在指定日期的可复查快照。私人仓库、未被索引、已删除及未来新增的内容无法穷尽。未把“搜索到仓库”“元数据可解析”和“完整科学流程运行成功”等同。

## 查找与筛选方法

1. 从 `agent-source.md`、BioTender 索引、GoekeLab genomic-skills 清单、Claw4Science、OpenClaw Medical 和 SkillsMP 等线索追踪原始来源。
2. 对 DNA 模型、变异、甲基化、泛基因组、长读长、古 DNA、线粒体、质粒、移动遗传元件、植物及合成生物学等方向做网页检索，并使用 GitHub 仓库检索补查。
3. 下载固定 Git 版本到 `temp/`，检查技能定义、支持文件、仓库与个别技能的许可证、作者声明、引用和镜像来源。未运行第三方安装器或下载模型／参考数据库。
4. 按任务、输入输出、阶段与必要条件判断重复。相同任务只保留一份；不同分析阶段、完整流程及组成步骤可以分别保留。名称相似不构成等价证明。
5. 将选中的完整技能资源、许可证、来源版本和文件哈希加入 `skills/`，记录兼容修复。验证目录完整性、搜索、安装与分发包；仅对少量已检查脚本做离线执行验证。

本文件列出 24 条保留了原始查询字符串的网页查询和 17 个 GitHub 查询／分页请求。bioinformatics 仓库查询取回当时的全部 5 页；DNA、variant 和宽泛 readme 查询为候选线索抽样，未声称翻遍所有搜索结果。完整查询与来源结果元数据见 [discovery.json](discovery.json)。

## 本轮主要补充

| 方向 | 示例 |
| --- | --- |
| DNA 模型与调控预测 | AlphaGenome、DNABERT-2、Nucleotide Transformer、Caduceus、Evo 2、Borzoi |
| 微生物／病毒与基因组挖掘 | 分箱、MAG 质控、病毒识别、代谢功能、生物合成与抗噬菌体防御 |
| 质粒与移动遗传元件 | MOB-suite、IntegronFinder、ISEScan、mobileOG、MobileElementFinder |
| 调控与多组学 | ENCODE 聚合、DNase-seq、ChIP-Atlas、染色质功能筛选 |
| 植物与遗传研究 | 古 DNA、植物遗传分析、基因组指导的引物设计、多倍体基因家族克隆、家系图与研究报告 |
| 实验与可复现性支持 | 研究设计、DNA 图表、FASTA 数据库整理、工作流记录与 Isabl 操作指南 |

## 许可与资源处理

每个收录来源的版本和依据见 [sources.json](sources.json)，逐路径选择见 [expansion.json](expansion.json) 与 [decisions.json](decisions.json)。许可范围另见 [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md)。

- DeepMind：README 区分 Apache-2.0 软件和 CC-BY-4.0 其他材料，两份完整许可和署名随技能分发。
- ENCODE Toolkit：独立技能保持 AGPL-3.0-only，分发完整可编辑文件、许可和引用；外部 MCP 服务另行配置。不能把这些文件改称 MIT。
- PKU-EMBL：所选 MIT 指南与模型权重、数据库条款分开；仓库公共文档复制到技能内部并调整链接。
- POD5：原仓库的 `bin` 符号链接指向作者本机虚拟环境，已移除该链接并把说明改为使用单独安装的 `pod5` 命令；没有附带虚拟环境。
- 少数 YAML 描述、技能标识和文档路径做了可追踪兼容修复，保留原始与分发文件的校验值。
- NGS tertiary、Isabl、ToolUniverse、ENCODE 等技能仍需要原项目运行环境或外部服务。安装指南不会自动提供这些计算工具。

## 重要未收录项

| 来源 | 原因 |
| --- | --- |
| [HKU-BAL/Clair-skills](https://github.com/HKU-BAL/Clair-skills/tree/c95b825f246d7484e0385c05aaffb42cecbaaf6e) | No applicable license for the skill repository found; the separate Clair caller licenses do not establish permission to redistribute this guide. |
| [JiaqiLi1024/s2f-agent](https://github.com/JiaqiLi1024/s2f-agent/tree/6d2d70dea2ae051bd17d88b4be09fdbfb96310ad) | DNA model workflows found, but no applicable content redistribution license found. Licenses of referenced model implementations do not license this skill collection. |
| [Pavel-Kravchenko/Bioinformatics](https://github.com/Pavel-Kravchenko/Bioinformatics/tree/1dc244624720bb43ad7bd17f3c3c49784ce8c65c) | No applicable redistribution license found at the pinned revision; retain metadata and original-source links only. |
| [openai/plugins](https://github.com/openai/plugins/tree/1dc195897af4161d039b80d8471ec0a10c9bbc89) | Reviewed the life-science-research subtree only. No applicable license declaration was found there or at the root; licenses in unrelated plugin subtrees do not cover these 50 skills. |
| [stjudecab/CAB-aiSkills](https://github.com/stjudecab/CAB-aiSkills/tree/04520240060c10e8937dcd945a01fc6ebbdc731e) | Reference only in this public general-use distribution: content has NonCommercial and ShareAlike conditions; it is not MIT and is not prohibited from all redistribution. |
| [Kinprove/dna-research-recipes](https://github.com/Kinprove/dna-research-recipes/tree/277bb68d89142ca51be52d2391713215e6180cbc) | Reference only: current FSL competing-use restriction and unfinished third-party data notice; the future MIT grant is not yet effective for the 2026 snapshot. |
| operon 中的 MethylGPT、gene-panel-design | 引用的仓库自有 Python 程序不存在于固定版本，暂缓为 `quality_review`。 |
| plant-skills 中的 grapevine-chip | 指向外部能力的兼容入口，当前仓库未提供完整实现，暂缓收录。 |

带 NonCommercial／竞争用途限制的内容不进入本次面向通用用户的公开分发包；这不等于宣称其一切再分发均被禁止。没有发现明确许可也不等于可以直接复制。镜像优先追溯原作者，不采用镜像根目录许可证替代原许可。

## 网页查询字符串

- `"DNA" "agent skills" genomics GitHub -writing -story -site:reddit.com`
- `"SKILL.md" "epigenomics" -site:reddit.com`
- `"SKILL.md" "CRISPR" -site:reddit.com`
- `"SKILL.md" "metagenomics" -site:reddit.com`
- `"SKILL.md" "AlphaGenome"`
- `"SKILL.md" "nucleotide transformer"`
- `"SKILL.md" "Evo2"`
- `"SKILL.md" "Borzoi"`
- `"SKILL.md" "mitochondrial"`
- `"SKILL.md" "ancient DNA"`
- `"SKILL.md" "plasmid"`
- `"SKILL.md" "DNA storage"`
- `"SKILL.md" "methylation" github`
- `"SKILL.md" "pangenome" github`
- `"SKILL.md" "long read" github`
- `"SKILL.md" "GWAS" github`
- `site:github.com genomics skills in:name`
- `"SKILL.md" "DNA" "license: MIT" -site:reddit.com`
- `"SKILL.md" "Oxford Nanopore" -site:reddit.com`
- `"SKILL.md" "synthetic biology" -site:reddit.com`
- `"genomics" "SKILL.md" "plant"`
- `"genetics" "SKILL.md" "license"`
- `"DNA" "SKILL.md" "origami"`
- `"SKILL.md" "STR" "genome"`

## 固定来源与处理结果

| 仓库 | Skill 路径 | 收录 | 许可／依据 |
| --- | ---: | ---: | --- |
| [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a) | 562 | 329 | MIT |
| [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb) | 208 | 21 | CC-BY-4.0 |
| [ClawBio/ClawBio](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d) | 99 | 20 | MIT |
| [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34) | 165 | 8 | MIT |
| [TianGzlab/OmicsClaw](https://github.com/TianGzlab/OmicsClaw/tree/6fbd79f262e3fd3552efe9f675574ffb165019ab) | 124 | 0 | Apache-2.0 |
| [NVIDIA-BioNeMo/bionemo-agent-toolkit](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit/tree/0e67a612e4045f007e38fa77adc8f3ebfc5616b6) | 62 | 3 | Apache-2.0 AND CC-BY-4.0 |
| [BioTender-max/awesome-bio-agent-skills](https://github.com/BioTender-max/awesome-bio-agent-skills/tree/8cbdd18837aa6296c4e77616a03b323bd69b57b1) | 1705 | 0 | CC0-1.0 |
| [FreedomIntelligence/OpenClaw-Medical-Skills](https://github.com/FreedomIntelligence/OpenClaw-Medical-Skills/tree/b1f9b6e3306f8d870e23b2023057bf730090036f) | 896 | 0 | 未建立适用再分发许可／仅作线索 |
| [variomeanalytics/bioinformatics-agent-skills](https://github.com/variomeanalytics/bioinformatics-agent-skills/tree/b3d959475e9b3c5a830368ff9a44319ce97e2549) | 0 | 0 | MIT |
| [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7) | 465 | 35 | Apache-2.0 |
| [AlterLab-IEU/AlterLab-Academic-Skills](https://github.com/AlterLab-IEU/AlterLab-Academic-Skills/tree/4a5b75358026b33d3e53101bf551331e12113bee) | 240 | 0 | MIT |
| [ArcInstitute/SRAgent](https://github.com/ArcInstitute/SRAgent/tree/67b60c24a49385d3eed8e5121be58d13c82839cd) | 1 | 0 | MIT |
| [BioTender-max/awesome-ai4bio-skills](https://github.com/BioTender-max/awesome-ai4bio-skills/tree/2d04b66d5831fb6274e1f34c6b8b5a7da80a47f3) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [Bioconductor/ai-agent-skills](https://github.com/Bioconductor/ai-agent-skills/tree/f87f5a70e228d98b62ab6f85628bd7fb9ba52147) | 15 | 0 | MIT |
| [Bloomed-Health/cancer-predisposition-variant-analyst](https://github.com/Bloomed-Health/cancer-predisposition-variant-analyst/tree/0cace7c9dfe9c7a49f6fac28896495f9f1ffe8e0) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [ChevyLiang/Geng-BEVA-agent-skills](https://github.com/ChevyLiang/Geng-BEVA-agent-skills/tree/3a26ad93a33c82c860cd7c65adbe7ff55a1055fa) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [DamLabResources/damlab-skills](https://github.com/DamLabResources/damlab-skills/tree/4ad096e3b9b0116271d10bdb746641b60a1b2211) | 16 | 2 | MIT |
| [Dddddong-higher/acr_skill](https://github.com/Dddddong-higher/acr_skill/tree/53f3710178f8a2c7184dc49313409b979df400f7) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [GoekeLab/awesome-genomic-skills](https://github.com/GoekeLab/awesome-genomic-skills/tree/c11b8d6b6fa0676229e8728d89ce78db3fa9d6aa) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [Golden2002/population-genetics-skill](https://github.com/Golden2002/population-genetics-skill/tree/c8fec68b6522f22befb5eaaf467750394209f722) | 1 | 0 | MIT |
| [HKU-BAL/Clair-skills](https://github.com/HKU-BAL/Clair-skills/tree/c95b825f246d7484e0385c05aaffb42cecbaaf6e) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [Jason-0409-G/vivarium](https://github.com/Jason-0409-G/vivarium/tree/a915249440e52879bac4013630533d78505e87d1) | 6 | 0 | MIT |
| [JiaqiLi1024/s2f-agent](https://github.com/JiaqiLi1024/s2f-agent/tree/6d2d70dea2ae051bd17d88b4be09fdbfb96310ad) | 17 | 0 | 未建立适用再分发许可／仅作线索 |
| [Jiawang1209/GeneFam-Pipeline-skills](https://github.com/Jiawang1209/GeneFam-Pipeline-skills/tree/aa87c4fae5f5e2844f401700186164e90a514f40) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [JimLiu/science-skills](https://github.com/JimLiu/science-skills/tree/fb309c32ee9a54dad169fa845638057d1cfac77f) | 29 | 0 | 未建立适用再分发许可／仅作线索 |
| [KalarisLabs/Myriad](https://github.com/KalarisLabs/Myriad/tree/5efbc5e9940fd37608f1d75cccac2befebad1215) | 3 | 1 | Apache-2.0 AND CC-BY-4.0 |
| [Kinprove/dna-research-recipes](https://github.com/Kinprove/dna-research-recipes/tree/277bb68d89142ca51be52d2391713215e6180cbc) | 1 | 0 | FSL-1.1-MIT |
| [Li-Chia-ching/Plant-qPCR-Agent-Skill](https://github.com/Li-Chia-ching/Plant-qPCR-Agent-Skill/tree/23acf68cc4cefe67f8f3b5f5057b0eabf0ac9a4e) | 1 | 1 | Apache-2.0 |
| [Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline](https://github.com/Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline/tree/ecce861e5215e25cc180f1381d73a8751c5257bd) | 1 | 1 | MIT |
| [Linxira-OS/linxira-bio-sdk](https://github.com/Linxira-OS/linxira-bio-sdk/tree/dcfd512a9a33272c05a5101aa2724abd62658924) | 35 | 0 | AGPL-3.0-only |
| [NathanSkene/claude-genome-skill](https://github.com/NathanSkene/claude-genome-skill/tree/ffac2ebe4d085830422fe73a580d4e8e099e7e92) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [Omicsio/gene2ai-openclaw-skill](https://github.com/Omicsio/gene2ai-openclaw-skill/tree/a972ea8f7f4e35c7e0f2db9dd23d46d6f83f3726) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [PKU-EMBL/Metagenomics-Skills](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258) | 163 | 53 | MIT |
| [Pavel-Kravchenko/Bioinformatics](https://github.com/Pavel-Kravchenko/Bioinformatics/tree/1dc244624720bb43ad7bd17f3c3c49784ce8c65c) | 213 | 0 | 未建立适用再分发许可／仅作线索 |
| [PengPPPP/MetaGLens-skills](https://github.com/PengPPPP/MetaGLens-skills/tree/f4e7f615723fc38d99ba06fe56cf057bb26c9e79) | 9 | 0 | 未建立适用再分发许可／仅作线索 |
| [Pratikrishi97/sciagent-skills](https://github.com/Pratikrishi97/sciagent-skills/tree/c9e6dd5ad2389af5c7933385bd8f6a65adc9be25) | 40 | 0 | MIT |
| [Runchuan-BU/BioClaw](https://github.com/Runchuan-BU/BioClaw/tree/50a7d9b1dfa536cfe7e047d63be1de368b47a814) | 52 | 0 | MIT |
| [SiYangming/bioskills](https://github.com/SiYangming/bioskills/tree/4d8327be29c11b1056bc8f68cfaa20f055a479e9) | 0 | 0 | MIT |
| [Skillwager/agent-skill-profiles](https://github.com/Skillwager/agent-skill-profiles/tree/660a444ebc249a1e51e675d9e8da0a71cd581cf8) | 0 | 0 | MIT |
| [Xuzhen-Li/bioinfo-agent-skills](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83) | 8 | 6 | MIT |
| [Ynzuses/research-skills-suite](https://github.com/Ynzuses/research-skills-suite/tree/a5cb708252a58bcfb771d46d39f5e38d64b8b2f6) | 8 | 2 | Apache-2.0 |
| [YuliaNuzhnenko/bioinformatics-agent-skills](https://github.com/YuliaNuzhnenko/bioinformatics-agent-skills/tree/9b9dbfe1802ef721a76c1e04efd305242042d1ca) | 12 | 0 | MIT |
| [agenticbioinformatics/skills](https://github.com/agenticbioinformatics/skills/tree/45b5e0c21d3a8d90ff870f994240627495f0dabc) | 4 | 0 | MIT |
| [ammawla/encode-toolkit](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b) | 94 | 16 | AGPL-3.0-only |
| [ankurgenomics/genomics-skills](https://github.com/ankurgenomics/genomics-skills/tree/c6f1adfc45fb2791c0051ea74fcf4b3a96c5de01) | 8 | 0 | MIT |
| [anthropics/life-sciences](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c) | 6 | 0 | 未建立适用再分发许可／仅作线索 |
| [aws-samples/sample-kiro-power-life-sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [awslabs/hcls-agent-skills](https://github.com/awslabs/hcls-agent-skills/tree/0cdbb370b2c53ea9af86eab186232751cbb2a2c0) | 41 | 0 | MIT-0 |
| [brycewang-stanford/Auto-Research-Skills](https://github.com/brycewang-stanford/Auto-Research-Skills/tree/1c277adf6833b1b7ac184ac5653f7b6b3c81d5c3) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [cameronehrlich/genex](https://github.com/cameronehrlich/genex/tree/ade282ff21f418987ba166af2066e29e2f1116cf) | 0 | 0 | MIT |
| [companion-inc/feynman](https://github.com/companion-inc/feynman/tree/dfdcb7cf2c73183cff7b10aa8ea8ce370c8b152c) | 48 | 1 | MIT |
| [cultivarium/full-circle-labs-skill](https://github.com/cultivarium/full-circle-labs-skill/tree/72eac4aea8810d9e424c95376cab24408cedae54) | 1 | 0 | MIT |
| [fmschulz/omics-skills](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637) | 34 | 10 | MIT |
| [gilleswainrib-ext-boop/bioinfo-skills](https://github.com/gilleswainrib-ext-boop/bioinfo-skills/tree/b632655f93a38df2c70f4d3fc42fc4d24b8fc7cd) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [google-deepmind/science-skills](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d) | 40 | 7 | Apache-2.0 AND CC-BY-4.0 |
| [hantao23/dna-storage-paper-reading-skill](https://github.com/hantao23/dna-storage-paper-reading-skill/tree/070636bc2cb283d93eafa020ee8e739555656ecc) | 1 | 0 | 未建立适用再分发许可／仅作线索 |
| [htlin222/ngs-tertiary-analysis-skills](https://github.com/htlin222/ngs-tertiary-analysis-skills/tree/d483c21fddd14a034e7f2f038f237ad7bcd05b1a) | 3 | 1 | MIT |
| [hubayirp/agentic-science](https://github.com/hubayirp/agentic-science/tree/213588a92bb74e705d795ba106e0119975877e8c) | 39 | 3 | MIT |
| [hydra-genetics/hydra-genetics-skills](https://github.com/hydra-genetics/hydra-genetics-skills/tree/a2c7acec765f3ed995eed68147e6cdd31abf9813) | 3 | 0 | 未建立适用再分发许可／仅作线索 |
| [jarxunlai/My-skills](https://github.com/jarxunlai/My-skills/tree/23fd1af77f8791f851c76bd6c7394d38d15763fe) | 16 | 0 | 未建立适用再分发许可／仅作线索 |
| [juanesarango/isabl-mcp](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e) | 8 | 8 | MIT |
| [lyfar/genomic-variant-research](https://github.com/lyfar/genomic-variant-research/tree/1a0e4f7f9161e6922a03539be14e42f945e47784) | 1 | 0 | MIT |
| [maratgaliev/claude-skill-genomic-pipelines](https://github.com/maratgaliev/claude-skill-genomic-pipelines/tree/b1fd8c09bdd257efe7962ecd3f771147ac35f74d) | 0 | 0 | 未建立适用再分发许可／仅作线索 |
| [masaomi/GenomicsChain_SkillSets](https://github.com/masaomi/GenomicsChain_SkillSets/tree/983bd50d61f97a92005fc4b9d900b95ab100d0f7) | 5 | 0 | MIT |
| [mcap91/bioinfo-agent-toolkit](https://github.com/mcap91/bioinfo-agent-toolkit/tree/ac16a7f57f67d586cda5158129b90a3ac72e8f6b) | 10 | 0 | MIT |
| [muwenbo/genetic-counseling-skills](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4) | 10 | 6 | MIT |
| [muwenbo/noxia-skills](https://github.com/muwenbo/noxia-skills/tree/8941b532a5c0cf3f81f94d3db2bae56752f85104) | 6 | 0 | 未建立适用再分发许可／仅作线索 |
| [nekrut/galaxy-claude-skills](https://github.com/nekrut/galaxy-claude-skills/tree/0f833c2f11fa624d9fd87c6d139dac404d6dd4ff) | 0 | 0 | MIT |
| [openai/plugins](https://github.com/openai/plugins/tree/1dc195897af4161d039b80d8471ec0a10c9bbc89) | 536 | 0 | 未建立适用再分发许可／仅作线索 |
| [stjudecab/CAB-aiSkills](https://github.com/stjudecab/CAB-aiSkills/tree/04520240060c10e8937dcd945a01fc6ebbdc731e) | 12 | 0 | CC-BY-NC-SA-4.0 |
| [swaruplab/operon](https://github.com/swaruplab/operon/tree/77f5361af5a2433e38ecfc652c454c49baec4ac3) | 690 | 0 | MIT |
| [tanaylab/misha](https://github.com/tanaylab/misha/tree/b13d4e2593d1883acd43efb30f9bf1132891fc3b) | 1 | 1 | MIT |
| [waldronlab/ai-agent-skills](https://github.com/waldronlab/ai-agent-skills/tree/01739035bb9522328ee84b64d7c3de4298ce1d17) | 10 | 0 | MIT |
| [wu-yc/LabClaw](https://github.com/wu-yc/LabClaw/tree/df378022e8b0764e5ee07df3e3a03e2fa7c7dd65) | 240 | 0 | 未建立适用再分发许可／仅作线索 |
| [zamushwani2/biomedical-ai-skills](https://github.com/zamushwani2/biomedical-ai-skills/tree/71c9884f8b1a2931d1b68d3ac3eafb9337da8d6d) | 16 | 0 | MIT |
| [zoebischuribe-cloud/postgwas-skills](https://github.com/zoebischuribe-cloud/postgwas-skills/tree/cc7ad3dbb68bb7c7b13e3ffde8366824f0a56a4f) | 2 | 0 | MIT |
| [zongtingwei/Bioclaw_Skills_Hub](https://github.com/zongtingwei/Bioclaw_Skills_Hub/tree/d056daa153685cc59d24e9c432e4804d9cd122fe) | 66 | 0 | MIT |

各来源的范围和个别许可例外以 [sources.json](sources.json) 的 `policy`、所选技能的 `license_basis` 及原始声明为准；表中仓库许可不能自动覆盖所有子目录。

## 本地复现与论文记录

维护命令见 [README.md](README.md)。已收录文件不依赖 `temp/`，PyPI 分发包包含目录、完整资源及许可。原始下载、Git 缓存和搜索响应位于 `temp/dna-skills-discovery-2026-09-12/` 与 `temp/skill-sources/`，不随包分发。

论文可报告该固定快照的覆盖、来源、去重方法、许可处理和检索／安装验证；不能据此声称已穷尽网络上的 DNA 技能或全部科学流程验证通过。小样本执行证据见 [validation.json](validation.json)。
