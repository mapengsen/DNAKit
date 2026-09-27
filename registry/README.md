# DNAKit Skills 来源与复现记录

快照日期：2026-09-25。Git commit 固定在 [sources.json](sources.json)。

下表保留原始扫描的 7114 个 `SKILL.md` 路径决策（包含聚合副本、示例及范围外内容）。`selected` 是原始扫描入选数量，不含后续独立清单的扩展，也不等于当前分发总量。这些路径数不是不同 DNA 功能的数量。

| 筛选结果 | 路径数 |
| --- | ---: |
| out_of_scope | 2153 |
| selected | 535 |
| duplicate | 519 |
| reference_only | 3849 |
| invalid_metadata | 42 |
| license_review | 13 |
| quality_review | 3 |

后续新增 18 个上游技能见 [genomics-new-skills.json](genomics-new-skills.json)，6 个原生流程见 [genomics-native-skills.json](genomics-native-skills.json)。用户提供的 335 项清单另补 13 个技能，见 [user-catalog-skills.json](user-catalog-skills.json)；[选择与排除审查](user-catalog-skills-review.json)记录逐来源处理及验证范围。最终分发 574 个技能，其中包含 DNAKit 检索入口 `dnakit` 和规划入口 `dnakit-plan`。

## 选择与去重方法

按任务、输入输出、分析阶段和必要条件判断是否重复。保留优先级及明确的等价组在 [curation_rules.py](../others/skills/curation_rules.py)。这一规则集可审查，但不是自动证明语义等价。
第二轮逐来源、逐路径的补充判断见 [expansion.json](expansion.json)，搜索范围与排除理由见 [discovery.md](discovery.md)。

例如：FASTQ 初始质控的重复指南保留一份；BWA 的 WGS/WES 与 Bowtie2 的 ChIP/ATAC 输入场景不同，分别保留。WGS 完整流程与其变异注释步骤也分别保留。多组学与新抗原组合流程的必要辅助技能标记为 `workflow_dependency`。

[decisions.json](decisions.json) 保留每个候选的来源、路径、筛选状态、替代项与理由。[catalog.json](catalog.json) 包含已收录内容、完整文件清单与 SHA-256。保留上游原文，必要的 YAML、技能标识、文档链接和本地环境路径兼容修复逐文件记录。原始与分发文件分别计算哈希；改动在 `UPSTREAM.md` 中说明。

## 来源与许可证处理

| 来源 | 固定版本 | 本版处理 |
| --- | --- | --- |
| [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a) | `d91ed3d56301` | Original-source DNA/genomics workflows; MIT notice retained. |
| [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb) | `fe505cae14d2` | Original content CC-BY-4.0; README identifies frontmatter licenses as underlying tool licenses. Retain referenced third-party notices. |
| [ClawBio/ClawBio](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d) | `1fcecb7e9531` | Individual skill license overrides apply; proprietary skills excluded; GPL skill deferred for separate review. |
| [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34) | `c1ed16d97dd6` | Only skills explicitly licensed MIT enter this snapshot; other individual terms require a separate review. |
| [TianGzlab/OmicsClaw](https://github.com/TianGzlab/OmicsClaw/tree/6fbd79f262e3fd3552efe9f675574ffb165019ab) | `6fbd79f262e3` | Root Apache and individual MIT scope needs clarification; record candidates without copying. |
| [NVIDIA-BioNeMo/bionemo-agent-toolkit](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit/tree/0e67a612e4045f007e38fa77adc8f3ebfc5616b6) | `0e67a612e404` | Retain Apache and CC-BY licenses plus NOTICE, with distinct code/content scopes and individual AND declarations. |
| [BioTender-max/awesome-bio-agent-skills](https://github.com/BioTender-max/awesome-bio-agent-skills/tree/8cbdd18837aa6296c4e77616a03b323bd69b57b1) | `8cbdd18837aa` | Discovery index only: aggregator license does not relicense original third-party skills; prefer original sources. |
| [FreedomIntelligence/OpenClaw-Medical-Skills](https://github.com/FreedomIntelligence/OpenClaw-Medical-Skills/tree/b1f9b6e3306f8d870e23b2023057bf730090036f) | `b1f9b6e3306f` | Discovery index only: root LICENSE absent at snapshot; a README badge is insufficient to relicense third-party content. |
| [variomeanalytics/bioinformatics-agent-skills](https://github.com/variomeanalytics/bioinformatics-agent-skills/tree/b3d959475e9b3c5a830368ff9a44319ce97e2549) | `b3d959475e9b` | Service/reference only: this repository contains no SKILL.md, and redistribution of the remote SkillGraph content is not established. |
| [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7) | `7a0ceb2f3a2a` | Root Apache-2.0; conflicting individual declarations deferred. ToolUniverse runtime remains a separate dependency. |
| [AlterLab-IEU/AlterLab-Academic-Skills](https://github.com/AlterLab-IEU/AlterLab-Academic-Skills/tree/4a5b75358026b33d3e53101bf551331e12113bee) | `4a5b75358026` | Reference/deduplication only. NOTICE and THIRD_PARTY_NOTICES distinguish K-Dense derivatives and four CC-BY-NC skills; do not infer MIT for every folder or ignore individual declarations. |
| [ArcInstitute/SRAgent](https://github.com/ArcInstitute/SRAgent/tree/67b60c24a49385d3eed8e5121be58d13c82839cd) | `67b60c24a493` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [BioTender-max/awesome-ai4bio-skills](https://github.com/BioTender-max/awesome-ai4bio-skills/tree/2d04b66d5831fb6274e1f34c6b8b5a7da80a47f3) | `2d04b66d5831` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [Bioconductor/ai-agent-skills](https://github.com/Bioconductor/ai-agent-skills/tree/f87f5a70e228d98b62ab6f85628bd7fb9ba52147) | `f87f5a70e228` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Bloomed-Health/cancer-predisposition-variant-analyst](https://github.com/Bloomed-Health/cancer-predisposition-variant-analyst/tree/0cace7c9dfe9c7a49f6fac28896495f9f1ffe8e0) | `0cace7c9dfe9` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [ChevyLiang/Geng-BEVA-agent-skills](https://github.com/ChevyLiang/Geng-BEVA-agent-skills/tree/3a26ad93a33c82c860cd7c65adbe7ff55a1055fa) | `3a26ad93a33c` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [DamLabResources/damlab-skills](https://github.com/DamLabResources/damlab-skills/tree/4ad096e3b9b0116271d10bdb746641b60a1b2211) | `4ad096e3b9b0` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Dddddong-higher/acr_skill](https://github.com/Dddddong-higher/acr_skill/tree/53f3710178f8a2c7184dc49313409b979df400f7) | `53f3710178f8` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [GoekeLab/awesome-genomic-skills](https://github.com/GoekeLab/awesome-genomic-skills/tree/c11b8d6b6fa0676229e8728d89ce78db3fa9d6aa) | `c11b8d6b6fa0` | Discovery/reference source; individual content terms require review before copying. |
| [Golden2002/population-genetics-skill](https://github.com/Golden2002/population-genetics-skill/tree/c8fec68b6522f22befb5eaaf467750394209f722) | `c8fec68b6522` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [HKU-BAL/Clair-skills](https://github.com/HKU-BAL/Clair-skills/tree/c95b825f246d7484e0385c05aaffb42cecbaaf6e) | `c95b825f246d` | No applicable license for the skill repository found; the separate Clair caller licenses do not establish permission to redistribute this guide. |
| [Jason-0409-G/vivarium](https://github.com/Jason-0409-G/vivarium/tree/a915249440e52879bac4013630533d78505e87d1) | `a915249440e5` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [JiaqiLi1024/s2f-agent](https://github.com/JiaqiLi1024/s2f-agent/tree/6d2d70dea2ae051bd17d88b4be09fdbfb96310ad) | `6d2d70dea2ae` | DNA model workflows found, but no applicable content redistribution license found. Licenses of referenced model implementations do not license this skill collection. |
| [Jiawang1209/GeneFam-Pipeline-skills](https://github.com/Jiawang1209/GeneFam-Pipeline-skills/tree/aa87c4fae5f5e2844f401700186164e90a514f40) | `aa87c4fae5f5` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [JimLiu/science-skills](https://github.com/JimLiu/science-skills/tree/fb309c32ee9a54dad169fa845638057d1cfac77f) | `fb309c32ee9a` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [KalarisLabs/Myriad](https://github.com/KalarisLabs/Myriad/tree/5efbc5e9940fd37608f1d75cccac2befebad1215) | `5efbc5e9940f` | LICENSE.md assigns procedural skill implementations to Apache-2.0 and explanatory/data materials to CC-BY-4.0. Preserve both full licenses, NOTICE.md, authorship and citation; do not describe the task graph as 10,000 installed skills. |
| [Kinprove/dna-research-recipes](https://github.com/Kinprove/dna-research-recipes/tree/277bb68d89142ca51be52d2391713215e6180cbc) | `277bb68d8914` | Reference only: current FSL competing-use restriction and unfinished third-party data notice; the future MIT grant is not yet effective for the 2026 snapshot. |
| [Li-Chia-ching/Plant-qPCR-Agent-Skill](https://github.com/Li-Chia-ching/Plant-qPCR-Agent-Skill/tree/23acf68cc4cefe67f8f3b5f5057b0eabf0ac9a4e) | `23acf68cc4ce` | Apache-2.0 notice retained. Import only the reviewed DNA-related scope and retain individual notices. |
| [Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline](https://github.com/Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline/tree/ecce861e5215e25cc180f1381d73a8751c5257bd) | `ecce861e5215` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Linxira-OS/linxira-bio-sdk](https://github.com/Linxira-OS/linxira-bio-sdk/tree/dcfd512a9a33272c05a5101aa2724abd62658924) | `dcfd512a9a33` | AGPL source collection; review scientific overlap and complete source requirements before importing. |
| [NathanSkene/claude-genome-skill](https://github.com/NathanSkene/claude-genome-skill/tree/ffac2ebe4d085830422fe73a580d4e8e099e7e92) | `ffac2ebe4d08` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [Omicsio/gene2ai-openclaw-skill](https://github.com/Omicsio/gene2ai-openclaw-skill/tree/a972ea8f7f4e35c7e0f2db9dd23d46d6f83f3726) | `a972ea8f7f4e` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [PKU-EMBL/Metagenomics-Skills](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258) | `e0e19a600d99` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Pavel-Kravchenko/Bioinformatics](https://github.com/Pavel-Kravchenko/Bioinformatics/tree/1dc244624720bb43ad7bd17f3c3c49784ce8c65c) | `1dc244624720` | No applicable redistribution license found at the pinned revision; retain metadata and original-source links only. |
| [PengPPPP/MetaGLens-skills](https://github.com/PengPPPP/MetaGLens-skills/tree/f4e7f615723fc38d99ba06fe56cf057bb26c9e79) | `f4e7f615723f` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [Pratikrishi97/sciagent-skills](https://github.com/Pratikrishi97/sciagent-skills/tree/c9e6dd5ad2389af5c7933385bd8f6a65adc9be25) | `c9e6dd5ad238` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Runchuan-BU/BioClaw](https://github.com/Runchuan-BU/BioClaw/tree/50a7d9b1dfa536cfe7e047d63be1de368b47a814) | `50a7d9b1dfa5` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [SiYangming/bioskills](https://github.com/SiYangming/bioskills/tree/4d8327be29c11b1056bc8f68cfaa20f055a479e9) | `4d8327be29c1` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Skillwager/agent-skill-profiles](https://github.com/Skillwager/agent-skill-profiles/tree/660a444ebc249a1e51e675d9e8da0a71cd581cf8) | `660a444ebc24` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Xuzhen-Li/bioinfo-agent-skills](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83) | `6dea4bee7aad` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [Ynzuses/research-skills-suite](https://github.com/Ynzuses/research-skills-suite/tree/a5cb708252a58bcfb771d46d39f5e38d64b8b2f6) | `a5cb708252a5` | Apache-2.0 notice retained. Import only the reviewed DNA-related scope and retain individual notices. |
| [YuliaNuzhnenko/bioinformatics-agent-skills](https://github.com/YuliaNuzhnenko/bioinformatics-agent-skills/tree/9b9dbfe1802ef721a76c1e04efd305242042d1ca) | `9b9dbfe1802e` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [agenticbioinformatics/skills](https://github.com/agenticbioinformatics/skills/tree/45b5e0c21d3a8d90ff870f994240627495f0dabc) | `45b5e0c21d3a` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [ammawla/encode-toolkit](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b) | `b1d90bc9e7ea` | Independent skill source files remain AGPL-3.0-only, as declared in pyproject.toml; complete editable skill sources, copyright, license and citation supplied. The separate ENCODE MCP runtime is not bundled or linked into DNAKit. |
| [ankurgenomics/genomics-skills](https://github.com/ankurgenomics/genomics-skills/tree/c6f1adfc45fb2791c0051ea74fcf4b3a96c5de01) | `c6f1adfc45fb` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [anthropics/life-sciences](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c) | `e96556b637b5` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [aws-samples/sample-kiro-power-life-sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) | `95b67e763df4` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [awslabs/hcls-agent-skills](https://github.com/awslabs/hcls-agent-skills/tree/0cdbb370b2c53ea9af86eab186232751cbb2a2c0) | `0cdbb370b2c5` | Amazon MIT-0 license retained; evaluate exact workflow overlap before importing. |
| [brycewang-stanford/Auto-Research-Skills](https://github.com/brycewang-stanford/Auto-Research-Skills/tree/1c277adf6833b1b7ac184ac5653f7b6b3c81d5c3) | `1c277adf6833` | Discovery/reference source; individual content terms require review before copying. |
| [cameronehrlich/genex](https://github.com/cameronehrlich/genex/tree/ade282ff21f418987ba166af2066e29e2f1116cf) | `ade282ff21f4` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [companion-inc/feynman](https://github.com/companion-inc/feynman/tree/dfdcb7cf2c73183cff7b10aa8ea8ce370c8b152c) | `dfdcb7cf2c73` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [cultivarium/full-circle-labs-skill](https://github.com/cultivarium/full-circle-labs-skill/tree/72eac4aea8810d9e424c95376cab24408cedae54) | `72eac4aea881` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [fmschulz/omics-skills](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637) | `7a4eb5bace3f` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [gilleswainrib-ext-boop/bioinfo-skills](https://github.com/gilleswainrib-ext-boop/bioinfo-skills/tree/b632655f93a38df2c70f4d3fc42fc4d24b8fc7cd) | `b632655f93a3` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [google-deepmind/science-skills](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d) | `28b8482603a4` | README separately licenses software Apache-2.0 and other materials CC-BY-4.0. Retain both full licenses, attribution and SKILL_LICENSES.md; third-party APIs/data have separate terms. |
| [hantao23/dna-storage-paper-reading-skill](https://github.com/hantao23/dna-storage-paper-reading-skill/tree/070636bc2cb283d93eafa020ee8e739555656ecc) | `070636bc2cb2` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [htlin222/ngs-tertiary-analysis-skills](https://github.com/htlin222/ngs-tertiary-analysis-skills/tree/d483c21fddd14a034e7f2f038f237ad7bcd05b1a) | `d483c21fddd1` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [hubayirp/agentic-science](https://github.com/hubayirp/agentic-science/tree/213588a92bb74e705d795ba106e0119975877e8c) | `213588a92bb7` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [hydra-genetics/hydra-genetics-skills](https://github.com/hydra-genetics/hydra-genetics-skills/tree/a2c7acec765f3ed995eed68147e6cdd31abf9813) | `a2c7acec765f` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [jarxunlai/My-skills](https://github.com/jarxunlai/My-skills/tree/23fd1af77f8791f851c76bd6c7394d38d15763fe) | `23fd1af77f87` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [juanesarango/isabl-mcp](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e) | `85b3e298d9cf` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [lyfar/genomic-variant-research](https://github.com/lyfar/genomic-variant-research/tree/1a0e4f7f9161e6922a03539be14e42f945e47784) | `1a0e4f7f9161` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [maratgaliev/claude-skill-genomic-pipelines](https://github.com/maratgaliev/claude-skill-genomic-pipelines/tree/b1fd8c09bdd257efe7962ecd3f771147ac35f74d) | `b1fd8c09bdd2` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [masaomi/GenomicsChain_SkillSets](https://github.com/masaomi/GenomicsChain_SkillSets/tree/983bd50d61f97a92005fc4b9d900b95ab100d0f7) | `983bd50d61f9` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [mcap91/bioinfo-agent-toolkit](https://github.com/mcap91/bioinfo-agent-toolkit/tree/ac16a7f57f67d586cda5158129b90a3ac72e8f6b) | `ac16a7f57f67` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [muwenbo/genetic-counseling-skills](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4) | `eda3cc5dace5` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [muwenbo/noxia-skills](https://github.com/muwenbo/noxia-skills/tree/8941b532a5c0cf3f81f94d3db2bae56752f85104) | `8941b532a5c0` | No applicable redistribution license established, or no portable SKILL.md at this revision; reference only. |
| [nekrut/galaxy-claude-skills](https://github.com/nekrut/galaxy-claude-skills/tree/0f833c2f11fa624d9fd87c6d139dac404d6dd4ff) | `0f833c2f11fa` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [openai/plugins](https://github.com/openai/plugins/tree/1dc195897af4161d039b80d8471ec0a10c9bbc89) | `1dc195897af4` | Reviewed the life-science-research subtree only. No applicable license declaration was found there or at the root; licenses in unrelated plugin subtrees do not cover these 50 skills. |
| [stjudecab/CAB-aiSkills](https://github.com/stjudecab/CAB-aiSkills/tree/04520240060c10e8937dcd945a01fc6ebbdc731e) | `04520240060c` | Reference only in this public general-use distribution: content has NonCommercial and ShareAlike conditions; it is not MIT and is not prohibited from all redistribution. |
| [swaruplab/operon](https://github.com/swaruplab/operon/tree/77f5361af5a2433e38ecfc652c454c49baec4ac3) | `77f5361af5a2` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [tanaylab/misha](https://github.com/tanaylab/misha/tree/b13d4e2593d1883acd43efb30f9bf1132891fc3b) | `b13d4e2593d1` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [waldronlab/ai-agent-skills](https://github.com/waldronlab/ai-agent-skills/tree/01739035bb9522328ee84b64d7c3de4298ce1d17) | `01739035bb95` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [wu-yc/LabClaw](https://github.com/wu-yc/LabClaw/tree/df378022e8b0764e5ee07df3e3a03e2fa7c7dd65) | `df378022e8b0` | Discovery mirror: no repository-level license at this revision. Resolve K-Dense, ToolUniverse and other original sources instead of relicensing copied skills. |
| [zamushwani2/biomedical-ai-skills](https://github.com/zamushwani2/biomedical-ai-skills/tree/71c9884f8b1a2931d1b68d3ac3eafb9337da8d6d) | `71c9884f8b1a` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [zoebischuribe-cloud/postgwas-skills](https://github.com/zoebischuribe-cloud/postgwas-skills/tree/cc7ad3dbb68bb7c7b13e3ffde8366824f0a56a4f) | `cc7ad3dbb68b` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [zongtingwei/Bioclaw_Skills_Hub](https://github.com/zongtingwei/Bioclaw_Skills_Hub/tree/d056daa153685cc59d24e9c432e4804d9cd122fe) | `d056daa15368` | Original repository MIT notice retained. Import only the reviewed DNA-related scope; preserve any individual license overrides. |
| [gkanogiannis/Gigwa-MCP](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf) | `d6ccdd1a5eb1` | Reviewed DNA workflows; complete skill resources and license notices retained. |
| [tayloranthonyanderson/haploqtl](https://github.com/tayloranthonyanderson/haploqtl/tree/210eace2912c7bf11d119cb503c89f8f6cad479f) | `210eace2912c` | Reviewed DNA workflows; complete skill resources and license notices retained. |
| [jmschrei/tangermeme](https://github.com/jmschrei/tangermeme/tree/d5c51b83d374576b64effb9284fc79a71a604734) | `d5c51b83d374` | Reviewed DNA workflows; complete skill resources and license notices retained. |
| [svdorn/statgen-skills](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0) | `dc5f0c60d299` | Reviewed DNA workflows; complete skill resources and license notices retained. |
| [genomoncology/biomcp](https://github.com/genomoncology/biomcp/tree/3b1d6452ab63741562d2ef97495a97281b304729) | `3b1d6452ab63` | Reviewed DNA workflows; complete skill resources and license notices retained. |
| [anthropics/life-sciences](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c) | `e96556b637b5` | Only the three selected skill subdirectories have verified Apache-2.0 LICENSE.txt; no repository-wide licensing claim. |
| [aws-samples/sample-kiro-power-life-sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) | `95b67e763df4` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |
| [nekrut/galaxy-claude-skills](https://github.com/nekrut/galaxy-claude-skills/tree/0f833c2f11fa624d9fd87c6d139dac404d6dd4ff) | `0f833c2f11fa` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |
| [ClawBio/ClawBio](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17) | `4848cc280c44` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |
| [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a) | `d91ed3d56301` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |
| [google-deepmind/science-skills](https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0) | `68832757cbbf` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |
| [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills/tree/330c8e764435a731eff571e3efdda70b363d0792) | `330c8e764435` | Additive reviewed DNA skill leaves; original baseline source metadata remains unchanged. Underlying tools/data retain separate terms. |

[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) 说明许可范围。SciAgent 的 README 明确区分 CC-BY-4.0 原创说明内容和 frontmatter 中的底层工具许可；K-Dense 和 ClawBio 则明确要求遵守各技能自身的许可覆盖规则。

MIT、Apache-2.0、CC-BY-4.0 和 AGPL-3.0-only 内容按原许可随包分发，不统一改成 MIT。ENCODE 独立技能保留完整可编辑源码；DeepMind 软件与文档分别适用 Apache 和 CC BY。暂缓项是对具体来源的判断，不表示 GPL 或 AGPL 一律禁止再分发。

### 暂不复制的候选

| 来源 | Skill | 理由 |
| --- | --- | --- |
| clawbio | wes-clinical-report-en | Explicit proprietary skill; no redistribution permission established. |
| clawbio | wes-clinical-report-es | Explicit proprietary skill; no redistribution permission established. |
| omicsclaw | genomics-alignment | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-assembly | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-cnv-calling | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-epigenomics | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-phasing | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-qc | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-sv-detection | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-variant-annotation | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-variant-calling | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| omicsclaw | genomics-vcf-operations | Root Apache-2.0 and per-skill MIT declarations need their scope clarified. |
| operon | epigenomics-methylgpt-agent | Guide invokes a repository-specific Python executable that is absent from the skill and pinned repository; a scientific runtime has not been established. |
| operon | gene-panel-design-agent | Guide invokes a repository-specific Python executable that is absent from the skill and pinned repository; a scientific runtime has not been established. |
| plant-skills | grapevine-chip | Compatibility stub refers to grapeancestry/chip outside this repository; complete canonical workflow not supplied. |
| tooluniverse | tooluniverse-gwas-study-explorer | Workflow body states MIT while the repository declares Apache-2.0; review pending. |

## 维护与复现命令

先安装项目开发依赖。网络下载和中间结果都写入仓库忽略的 `temp/`。 `--fetch` 获取缺失的固定 Git 版本，不运行上游安装程序。

```bash
python -m pip install -e ".[dev]"
python others/skills/build_catalog.py --fetch --plan
# 审阅 temp/dnakit-skills-work/ 中的选择与决策后重新生成
python others/skills/build_catalog.py --fetch
python others/skills/build_catalog.py --check
# 验证新增上游技能的固定版本资源；不运行科学工作流
python others/skills/import_genomics_skills.py --check
python others/skills/sync_genomics_expansion.py --check
# 更新全部已审阅的 DNAKit 适配
python others/skills/adaptations.py
# --python 指向 Python 3.11+；运行已检查的离线脚本
python others/skills/smoke_workflows.py --python python3
python others/skills/write_inventory.py
python -m build --outdir temp/dnakit-skills-dist
```

`skills/` 是公开分发目录；原生适配的维护模板在 `others/skills/adaptations/`，由 `registry/adaptations.json` 绑定上游版本。完整导入也会重放这些适配，保留旧名称别名、原始元数据、许可和哈希。新增上游技能由 `genomics-new-skills.json` 固定完整资源；新增原生入口由 `genomics-native-skills.json` 维护 API 与类别。`sync_genomics_expansion.py` 合并这些清单，完整目录构建也会保留它们。科学后端保留原始实现。构建时将已登记文件复制进 wheel。用户安装包后不需要 `temp` 中的原仓库，查询和安装也不需要 Git 或 Node.js。

## 迁移清单

[migration.json](migration.json) 逐项列出技能的接入方式、函数和脚本。
[native-providers.json](native-providers.json) 固定数据库封装，
[native-extensions.json](native-extensions.json) 记录 DNAKit 可选参数扩展。
模板重放不运行科研工作流；执行验证须另行记录。

## 验证范围与论文使用

[validation.json](validation.json) 记录已安装技能的离线子进程检查：BED/GFF 坐标转换、合成样本的 QC 标记和报告，以及 FASTA 去重、映射报告与异常输入处理。证据绑定文件清单哈希；文件变化后不继承旧结果。

这些检查不包含大模型调用、外部数据库、全规模测序任务或完整流程成功率。元数据有效、文件哈希一致、依赖探测成功和科学结果正确是不同验证层级。

论文可使用 [skills.csv](skills.csv) 和固定版本清单报告覆盖、来源、许可及去重结果。DNAKit 的贡献应表述为筛选、整合、检索、封装与工具接入；上游技能、软件、模型和数据库需分别署名引用。`dnakit skills cite ID --format bibtex` 提供固定版本的软件引用，不是虚构的论文 DOI。
