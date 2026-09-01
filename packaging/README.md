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

正式发布时，需要把配方复制到 `bioconda-recipes` 仓库的
`recipes/dnakit/meta.yaml`，并向 Bioconda 提交 Pull Request。合并后才可使用：

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

正式进入 GNU Guix 仍需按 Guix 项目的贡献流程提交补丁并等待审核。
合并后可直接运行 `guix install python-dnakit`。

## Galaxy

Galaxy Tool Shed 包装器位于 [`galaxy/dnakit/`](../galaxy/dnakit/) 中。Galaxy
依赖 Bioconda 提供 `dnakit=0.1.1`，因此应先完成 Bioconda 发布。
