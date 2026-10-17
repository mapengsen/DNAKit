# DNAKit 接入说明

正式技能名称：`dnakit-tool-generation`。在 Codex 中安装：

```bash
dnakit skills install dnakit-tool-generation --agent codex
```

本 Skill 指导 Agent 先用 DNAKit MCP 确认已有工具能否完成需求，确实缺少能力时，再在隔离的候选代码目录中实现并检查新的候选工具。它不附带脚本，也不会自行登记 manifest、重启 MCP 或提交代码；这些步骤由独立评测程序和运行控制器决定。

使用前需要：连接 DNAKit MCP（或等效的 DNAKit CLI/API）、一个已 clone 并固定 commit 的候选代码目录，以及预先配置好的 `admet_ai` conda 环境。安装 Skill 不会配置这些环境。完整流程见 [SKILL.md](SKILL.md)。
