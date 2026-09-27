# DNAKit Skills：搜索、安装与使用

DNAKit 汇总 DNA／基因组学工作流，提供本地检索、按需安装、来源追溯和 MCP 读取。一个 Skill 是一份操作指南及其脚本、参考资料；安装 Skill 后，还需要准备该流程要求的软件、参考基因组或 API。

建议先安装 **`dnakit` 统一入口**。用户提出自然语言任务，Agent 根据入口说明搜索目录，读取合适的工作流，再使用可用工具执行。无需预先记住几百个技能名称。

## 1. 安装 DNAKit

以下功能需要包含 Skills 目录的新版本。当前源码可在 DNAKit 仓库根目录安装：

```bash
python -m pip install -e "."
dnakit skills stats
```

发布包含本功能的版本后，用户可以通过 `python -m pip install -U dnakit` 获取。仅搜索、读取和安装 Skill 不需要 MCP 或模型 API。

## 2. 安装统一入口

在自己的分析项目目录执行，三种客户端选择一种：

```bash
dnakit skills install dnakit --agent codex
# 或
dnakit skills install dnakit --agent claude-code
# 或
dnakit skills install dnakit --agent cursor
```

默认只安装到当前项目。需要在所有项目使用时，添加 `--global`。

| 客户端 | 项目目录 | 用户目录（`--global`） |
| --- | --- | --- |
| Codex | `.agents/skills/dnakit/` | `~/.agents/skills/dnakit/` |
| Claude Code | `.claude/skills/dnakit/` | `~/.claude/skills/dnakit/` |
| Cursor | `.cursor/skills/dnakit/` | `~/.cursor/skills/dnakit/` |

目录依据 [Codex Skills 文档](https://learn.chatgpt.com/docs/build-skills)、[Claude Code Skills 文档](https://code.claude.com/docs/en/skills)及 [Cursor Skills 文档](https://cursor.com/docs/skills)。重新打开客户端会话，让客户端发现新技能；客户端运行环境应能访问安装了 DNAKit 的 Python 环境。

然后在 Codex 中输入，例如：

```text
$dnakit 把 BED 区间 chr1 999 1000 转成 GFF 坐标。
请先找合适的技能，检查依赖，运行转换并解释坐标约定。
```

预期过程：检索 `genomic-coordinates` → 阅读说明和依赖 → 使用附带脚本转换 → 返回 `chr1 1000 1000`，区间长度为 1 bp。Agent 在当前任务授权范围内执行，具体工具取决于所选工作流。

### 用 `dnakit-plan` 先规划复杂任务

仓库在 `skills/dnakit-plan/` 收录 Codex 专用的规划 Skill，用于在执行前调查 DNAKit MCP 的可用工具、询问影响结果的关键问题，并在用户确认后保存可供下次会话执行的 Markdown 计划。仓库自身的 `.agents/skills/dnakit-plan/` 是指向该正式目录的相对软链接；修改正式目录后，本项目中的 Codex 会读取新内容。其他项目可以安装目录中的版本：

```bash
dnakit skills install dnakit-plan --agent codex
# 要在所有项目使用，添加 --global
```

安装后重新打开 Codex 会话。安装 Skill 不会自动配置 DNAKit MCP；请按 [MCP 教程](mcp.md)连接服务。已在本仓库根目录使用 Codex 的用户无需重复安装。

需要确保每次进入规划流程时，请显式使用 `$dnakit-plan`；仅安装 Skill 不会强制 Codex 对所有请求先规划。

先按 [MCP 教程](mcp.md)接入 DNAKit，并在 Codex 中用 `/mcp` 确认工具可用。之后只需输入任务目标，例如：

```text
$dnakit-plan 下载 hg38 参考基因组并准备训练数据。
```

Skill 要求先通过 DNAKit MCP 搜索所需工具并核对参数，在 Markdown 中逐步列出具体工具名、关键参数、输入输出和验收方法；关键问题、确认后保存以及停止执行也由 Skill 规定。按提示回答问题并确认最终方案即可；在允许写文件的模式中，Codex 会保存到 `plans/`，无需重复说明这些规则。下次只需说“按 `plans/<任务名>.md` 执行”。`/plan` 是可选的额外规划模式；若该模式不允许写文件，确认后切回普通模式再保存。规划 Skill 提供流程指引；若必须由程序保证“未经确认绝不执行”，还需在调用 Codex 的入口设置独立的确认门禁。

## 3. 自己搜索并安装某个技能

```bash
dnakit skills search "变异注释" --limit 5
dnakit skills search "FASTQ quality" --limit 5
dnakit skills search "GWAS" --category gwas --limit 5
dnakit skills search "整合子" --limit 5
dnakit skills search "DNABERT" --category models --limit 5
dnakit skills categories

dnakit skills info variant-annotation
dnakit skills doctor variant-annotation
dnakit skills install variant-annotation --agent codex
```

搜索使用本地关键词索引和中文查询别名，不调用模型 API。Agent 搜索时返回的摘要、之后读取的正文仍会占用对话 token。底层工具检索与 Skill 检索是两件事：前者找可调用函数，后者找如何完成任务的流程。

`doctor` 只检查已登记的 Python 模块、命令、环境变量是否存在；完整版本要求、R 包、模型／数据库许可、网络接口及参考文件还需阅读说明。检查通过不代表完整流程已验证。它不会输出 API 密钥值。

安装会复制所选技能的全部已收录资源和许可证。已有同名非 DNAKit 目录会报错；已安装内容相同时不重复写入；修改过的 DNAKit 技能只有显式 `--force` 才会替换。安装不会自动安装计算依赖。

## 4. 一个无需模型 API 的本地例子

`genomic-coordinates` 的上游兼容要求为 Python 3.11+。将完整技能装到自己的工作目录：

```bash
dnakit skills install genomic-coordinates --destination ./analysis-skills
python3 analysis-skills/genomic-coordinates/scripts/convert_coords.py \
  --from bed --to gff chr1 999 1000 --format json
```

结果中的 `start=1000`、`end=1000`、`length=1` 表示从 BED 的 0-based 半开区间转换为 GFF 的 1-based 闭区间。这里改变的是坐标约定；参考基因组版本转换还需要对应参考文件或映射方法。

## 5. 通过 MCP 按需读取

按 [MCP 教程](mcp.md)接入 DNAKit 后，Agent 可以直接调用：

```text
search_dnakit_skills(query="变异注释", limit=5)
read_dnakit_skill(skill_id="variant-annotation")
```

读取辅助文档时，传入响应 `files` 中的相对路径作为 `resource`。响应出现 `next_offset` 时，以该值继续读取。MCP 文本读取不会把脚本复制到客户端。服务环境可以直接调用 `dnakit.workflows.run` 运行包内脚本，也可以使用 `dnakit skills install ... --destination ...` 导出整个技能。

全部收录技能都提供 `DNAKIT.md` 接入说明。数据库和集合封装已迁入 DNAKit；实际科学软件、模型、参考数据和凭据继续按任务准备。详见[统一接口教程](#unified-skill-api)。
已经迁移的 `sequence-retrieval` 直接调用 DNAKit 的 NCBI/ENA 接口，无需 ToolUniverse。
FASTA、NCBI GenBank 和 ENA EMBL 使用基础包；跨格式转换需要 `dnakit[bio-formats]`。
详见[序列检索](../api/features/22_database_query_sequences.md)和[序列下载](../api/features/15_download.md#accession-sequences)。

```bash
dnakit skills search "sequence retrieval" --limit 5
dnakit skills install sequence-retrieval --agent codex
dnakit skills read sequence-retrieval --file examples.md
```

旧名称 `tooluniverse-sequence-retrieval` 和 `public-db-lookup` 仍作为别名解析到新技能。
原有客户端安装目录不会自动删除或移动；安装新名称后，可按需清理旧目录以免重复发现。
其他技能的接入状态以 `skills info` 中的运行依赖和技能说明为准。

## 6. Skill 统一接口 {#unified-skill-api}

直接使用函数时，请按任务进入对应教程：[数据库查询](../api/features/22_database_query_discovery.md)、[下载](../api/features/15_download.md)、[变异注释](../api/features/35_variant_annotation.md)、[变异效应评分](../api/features/36_variant_effect_scores.md)、[自定义模型训练](../api/features/23_prediction_custom_training.md)、[调控轨道预测](../api/features/23_prediction_regulatory_tracks.md)、[蛋白质结构预测](../api/features/23_prediction_protein_features.md)或[工作流与报告](../api/features/14_engineering.md)。参数、返回值和调用示例均在对应教程中。

DNAKit 为整个技能目录提供统一的搜索、安装和执行入口。集合封装中的数据库访问、共享脚本、报告与调用代码迁入 DNAKit；Biopython、samtools、AlphaGenome 等实际科学软件或模型继续使用原实现。

Skill 名称去掉集合前缀，旧名称保留为别名；真实软件名称继续保留。例如 `tooluniverse-sequence-retrieval` 对应 `sequence-retrieval`，`virulencefinder` 和 `alphagenome-atlas-website-links` 保持原名。

### 数据库函数

先搜索所需接口，再查看参数：

```bash
dnakit database search "gene lookup"
dnakit database info ensembl_lookup_gene
dnakit database call ensembl_lookup_gene '{"id":"BRCA1","species":"homo_sapiens"}'
```

Python 直接使用能力函数，无需初始化 ToolUniverse 或加载全部工具：

```python
import dnakit

gene = dnakit.databases.ensembl_lookup_gene(
    id="BRCA1",                # 也接受原参数 gene_id
    species="homo_sapiens",
)
variants = dnakit.databases.gnomad_get_sv_by_region(
    chrom="17", start=43044295, end=43070295,  # end 也可写成 stop
)
print(gene)
print(variants)
```

`dnakit.databases.describe()` 的 `input_schema` 描述 DNAKit 接受的参数，`config.parameter` 保留原服务参数。别名出现冲突、参数不合法、提供方报错或结果超出上限时抛出异常。结果保留提供方的字典、列表和字段结构，不能假定每个接口都有相同的 `data` 层。

可选参数会实际参与处理。基因名称可按接口要求解析为 ID，坐标别名会转换后调用；Open Targets 的 `size` 限制当前返回页，并标记剩余结果，GWAS 的 `p_value_threshold` 对当前页实施过滤。当前页结果不能被视为完整数据库检索结果。

### 随包脚本

用户不需要克隆原技能仓库来定位已打包的脚本：

```bash
dnakit skills scripts genomic-coordinates
dnakit skills run genomic-coordinates --script scripts/convert_coords.py -- \
  --from bed --to gff chr1 999 1000 --format json
```

Python 调用同一份脚本，保留原来的可选参数：

```python
import dnakit

result = dnakit.workflows.run(
    "genomic-coordinates",
    "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60,
)
print(result["returncode"])
print(result["output"])
```

该示例把 BED 的 0-based 半开区间转换成 GFF 的 1-based 闭区间；不是参考基因组版本转换。CLI 转发标准输入、输出和退出码；Python 接口限制执行时间和捕获输出大小。脚本及其辅助资源在执行前核对已登记的 SHA-256。流程依赖的科学软件仍须位于运行环境中。

### ENCODE、组合流程与报告

ENCODE 的查询、下载和本地追踪使用同一个 DNAKit 包：

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)
```

ENCODE 函数保留原接口选项和 JSON 文本返回形式。本地追踪及凭据使用 `~/.dnakit-encode/`；持久化凭据存储需要 `dnakit[encode-credentials]`。源代码许可为 AGPL-3.0-only，详见第三方声明。

nf-core 等组合流程需要的已许可支持脚本可通过 `dnakit skills support-run NAME -- ...` 或 `dnakit.workflows.run_support()` 执行。支持步骤包括 `clinical-variant-reporter`、`clinical-trial-finder` 、`omics-target-evidence-mapper`、`gwas-lookup` 和 `fine-mapping`。这些函数运行原流程，不保证其结论已经通过临床验证。

DNAKit 自己的报告函数接收已有 Markdown 内容，生成 PDF 和来源哈希，不自动生成新的医学解释：

```bash
python -m pip install -e ".[reports]"
dnakit skills report --input wes_report.md --output report-pdf --language en
```

对应函数为 `dnakit.workflows.write_wes_report()`，提供语言标签、标题、左右 Logo、字体和覆盖选项。没有复制上游标为专有许可的英文或西班牙文报告实现。

### MCP 与验证范围

按 [MCP 教程](mcp.md)安装 `dnakit[agent]` 后，Agent 使用 `search_dnakit_tools` → `describe_dnakit_tool` → `call_dnakit_tool` 发现并调用上述函数。脚本执行、模型计算、文件写入等操作带有相应的操作类型；目录增加函数不会把全部定义放入模型的初始上下文。

每个 Skill 的 `DNAKIT.md`、运行依赖和验证状态可以通过 `dnakit skills info ID --json` 查看。`registry/migration.json` 记录逐项接入方式，`registry/adaptations.json` 绑定维护模板和上游版本。函数签名、文件完整性、单元测试、真实服务调用及完整科研流程属于不同验证层级；未经执行的流程继续标为 `not_run`。

移入的代码保留原作者、许可及原始文件哈希。DNAKit 新增的参数转换、执行入口和报告函数另行记录，不把上游模型或代码改称原创。RCSB 坐标文件接口依据 [RCSB 官方下载说明](https://www.rcsb.org/docs/programmatic-access/file-download-services)，CELLxGENE 有界查询依据 [官方 Python API](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html)。

ChIP-Atlas 的基因区域查询先返回 `requires_data` 状态及峰文件下载链接；提供核对过的 `peak_file` 后才计算区域重叠。原富集接口只提供提交说明，标为 `requires_submission`，不能当成已算出的统计结果。区域重叠也不能直接证明 motif 被破坏，详见 [ChIP-Atlas 官方 API](https://chip-atlas.org/agents)。
