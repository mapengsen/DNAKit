# 基因组育种与农业

本类共 **1** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`dnakit-genomic-breeding`](#dnakit-genomic-breeding) | dnakit | Compute diploid genomic relationships and fixed-lambda RR-BLUP predictions with DNAKit. Use for genotype-to-phenotype breeding research when aligned dosage mat… |

## `dnakit-genomic-breeding` {#dnakit-genomic-breeding}

Compute diploid genomic relationships and fixed-lambda RR-BLUP predictions with DNAKit. Use for genotype-to-phenotype breeding research when aligned dosage matrices and training phenotypes are available.

- **来源：**[DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-genomic-breeding)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `numpy`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit-genomic-breeding
dnakit skills install dnakit-genomic-breeding --agent claude-code   # 或 --agent codex / cursor
```
