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

全部收录技能都提供 `DNAKIT.md` 接入说明。数据库和集合封装已迁入 DNAKit；实际科学软件、模型、参考数据和凭据继续按任务准备。详见[统一接口教程](native_skills.md)。
已经迁移的 `sequence-retrieval` 直接调用 DNAKit 的 NCBI/ENA 接口，无需 ToolUniverse。
FASTA、NCBI GenBank 和 ENA EMBL 使用基础包；跨格式转换需要 `dnakit[bio-formats]`。
详见[序列检索](../api/features/22_database_query.md)和[序列下载](../api/features/15_download.md#accession-sequences)。

```bash
dnakit skills search "sequence retrieval" --limit 5
dnakit skills install sequence-retrieval --agent codex
dnakit skills read sequence-retrieval --file examples.md
```

旧名称 `tooluniverse-sequence-retrieval` 和 `public-db-lookup` 仍作为别名解析到新技能。
原有客户端安装目录不会自动删除或移动；安装新名称后，可按需清理旧目录以免重复发现。
其他技能的接入状态以 `skills info` 中的运行依赖和技能说明为准。

## 6. 来源、去重、许可证与论文引用

当前目录收录 bioSkills、ToolUniverse、ClawBio、SciAgent、K-Dense、NVIDIA BioNeMo，以及 DeepMind science-skills、PKU-EMBL Metagenomics-Skills、omics-skills、ENCODE Toolkit 等来源的工作流，并提供 DNAKit 自有入口。内容覆盖 DNA 模型、变异、调控、宏基因组、病毒、质粒和移动遗传元件等任务。同一任务阶段的重复实现按已审核规则保留一个；输入条件不同的流程、完整工作流及其组成步骤可以分别保留。旧名称无歧义时，已登记的重复项别名会指向保留版本；同名但任务不同的情况，请从搜索结果选择具体技能 ID。

查询当前数量和具体来源：

```bash
dnakit skills stats
dnakit skills list --source tooluniverse --limit 20
dnakit skills info genomic-coordinates --json
dnakit skills cite genomic-coordinates
dnakit skills cite genomic-coordinates --format bibtex
```

每个上游技能都带有 `UPSTREAM.md`、原始许可证、固定 Git commit 和文件校验值。完整清单及去重／暂缓理由见仓库的 [`registry/`](https://github.com/mapengsen/DNAKit/tree/main/registry)，论文整理可使用 [`registry/skills.csv`](https://github.com/mapengsen/DNAKit/blob/main/registry/skills.csv)。

仓库核心代码保持 MIT；随包分发的第三方内容按各自 MIT、Apache-2.0、CC-BY-4.0、AGPL-3.0-only 条款保留。ENCODE 的独立技能保持 AGPL 并附完整可编辑源码；DeepMind 的软件和说明内容分别适用 Apache-2.0 与 CC-BY-4.0。此处收录的是技能指南及辅助文件，外部 MCP 服务、模型权重和数据库仍需按各自条件准备。

检索范围、关键词、固定来源和未收录原因见 [`registry/discovery.md`](https://github.com/mapengsen/DNAKit/blob/main/registry/discovery.md)。没有找到适用许可证的来源、当前带非商业或竞争用途限制的内容，以及缺失关键脚本的候选，没有复制进此通用分发包。发现平台只作为线索，优先回到原作者来源。该清单是公开网络在指定日期的可复查快照，不能保证穷尽所有 DNA 技能。

论文应分别说明 DNAKit 的筛选、去重、检索、封装与集成贡献，以及上游工作流贡献。除 DNAKit 外，还要引用实际使用的上游项目、软件、模型和数据集。目录纳入、元数据校验、文件完整性检查不等于整套流程经过科学验证；具体执行验证范围记录在目录和 `registry/validation.json` 中。
