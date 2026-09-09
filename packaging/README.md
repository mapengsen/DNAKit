# 平台打包适配

这里保存 DNAKit 在 PyPI 之外的平台适配源文件。它们不会改变
`src/dnakit/` 中的核心代码。

## Bioconda

[`bioconda/meta.yaml`](bioconda/meta.yaml) 是 `dnakit 0.1.1` 的 Bioconda
配方源稿，使用 PyPI 已发布的源码包和固定 SHA-256。

本地渲染：

```bash
conda-build packaging/bioconda --output
```

Bioconda 已收录 DNAKit，0.1.1 和 0.1.3 均已发布；本目录仍保留
Galaxy 当前依赖的 0.1.1 配方。现在可以使用：

```bash
conda install -c conda-forge -c bioconda dnakit
```

## GNU Guix

[`guix/dnakit.scm`](guix/dnakit.scm) 是可独立加载的 Guix 包定义，包名为
`python-dnakit`，使用 GitHub 版本标签源码和固定哈希。

在安装了 Guix 的系统中验证和安装：

```bash
guix build -f packaging/guix/dnakit.scm
guix package -f packaging/guix/dnakit.scm
```

正式进入 GNU Guix 仍需等待官方构建通过并完成审核合并。
合并后可直接运行 `guix install python-dnakit`。

当前提交：[Guix PR #10903](https://codeberg.org/guix/guix/pulls/10903)。
改动已整理为单个提交 `6edc2d2` 并推送到原 PR。
2026-09-09 修复了官方 CI 的五项测试失败：测试运行时生成的脚本写死了
`#!/bin/sh`，而 Guix 构建沙箱没有 `/bin/sh`。配方在解包后将两个测试文件
中的脚本解释器替换为构建环境提供的 `sh`，并保留 pytest 的 `pythonpath=.`
设置以加载源码中的 `benchmarks` 和 `validation`；没有禁用测试。

本地使用 Python 3.12 和 Guix pytest 插件，在缺少 `/bin/sh` 的隔离环境中
先复现了同样的五项失败，再执行实际的 Guix `substitute*` 修复阶段，
完整测试结果为 **1001 passed、15 skipped**。已验证此独立配方与 PR 中的
包定义一致；这不等同于完整 Guix 构建，后者仍需由 Cuirass 确认。

## Galaxy

Galaxy Tool Shed 包装器位于 [`galaxy/dnakit/`](../galaxy/dnakit/) 中。Galaxy
依赖 Bioconda 提供 `dnakit=0.1.1`，该版本现已发布。
2026-09-07 已发布到
[Main Tool Shed](https://toolshed.g2.bx.psu.edu/view/mapengsen/dnakit/9f6dc4f3006b)，
修订号为 `9f6dc4f3006b`，包含四个 `0.1.1+galaxy0` 工具。
本地 `planemo lint`、`planemo shed_lint` 及 Galaxy 24.2 的四项工具测试
全部通过。运行测试时将独立安装的 Bioconda 0.1.1 环境加入 `PATH`，
使用 `--no_dependency_resolution`；未测试 Galaxy 自动创建 Conda 环境。
验证和上传命令见 [`galaxy/dnakit/README.rst`](../galaxy/dnakit/README.rst)。
