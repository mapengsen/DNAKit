# DNAKit

**DNAKit: A Comprehensive Toolkit for Efficient DNA Research**

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases=true&label=PyPI&logo=pypi&cacheSeconds=300" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda">
    <img src="https://img.shields.io/badge/Bioconda-recipe-43B02A?logo=anaconda&logoColor=white" alt="Bioconda recipe">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit">
    <img src="https://img.shields.io/badge/Galaxy-wrapper-2C3143?logo=galaxy&logoColor=white" alt="Galaxy wrapper">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/guix">
    <img src="https://img.shields.io/badge/GNU%20Guix-package-A42E2B?logo=gnu&logoColor=white" alt="GNU Guix package">
  </a>
  <a href="https://dnakit.mapengsen.asia">
    <img src="https://img.shields.io/badge/Web%20demo-online-2EA44F?logo=googlechrome&logoColor=white" alt="Web demo">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Windows-0078D6?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDBoMTEuNHYxMS40SDB6TTEyLjYgMEgyNHYxMS40SDEyLjZ6TTAgMTIuNmgxMS40VjI0SDB6TTEyLjYgMTIuNkgyNFYyNEgxMi42eiIvPjwvc3ZnPg==" alt="GUI for Windows">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-macOS-000000?logo=apple&logoColor=white" alt="GUI for macOS">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Linux-FCC624?logo=linux&logoColor=black" alt="GUI for Linux">
  </a>
</p>

DNAKit 是面向 DNA/基因方便操作的可复现分析工具包。它提供：

1. **数据准备：** DNA 序列与注释文件读写、数据库下载与查询、合法性检查。
2. **数据处理：** 序列标准化与编辑、转录与翻译、去重与聚类、数据集划分。
3. **数据分析：** DNA 描述符、表征与指纹、序列与模式搜索、理化及热力学性质、深度学习性质预测。
4. **数据评价：** 有效性、唯一性、多样性与新颖性评价，以及相似度、距离和序列比对。
5. **可视化：** 序列文字图、位置高亮、线性与环状 DNA 图、序列比对图。
6. **大语言模型工具：** 提供 [Python API](llm_tools/python_api.md)、[MCP 接入](llm_tools/mcp.md)、[Skills 工作流](llm_tools/skills.md)与内置 [Agent](llm_tools/agent.md)，支持通用大模型调用与自然语言任务执行。
7. **图形界面与网页版：** [网页版](https://dnakit.mapengsen.asia)无需安装，打开浏览器即可使用；[桌面版](gui.md)（`dnakit gui`）支持 Windows、macOS 和 Linux，提供序列与数据集工具、全部函数的检索与运行以及 Skills 目录。

# 支持与社区

**Github仓库：**[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

<div id="wechat-community" class="wechat-community" aria-labelledby="wechat-title">
  <strong id="wechat-title">微信社区：</strong>
  <img src="assets/images/community/wechat-community-qr.png" alt="DNAKit 微信社区联系人二维码" loading="lazy">
</div>

# 引用 DNAKit

DNAKit 尚未发布论文或 DOI；当前版本的中文引用信息保存在仓库 `README_CN.md` 中，可引用为：

```text
DNAKit contributors. DNAKit 0.2.0, 2026.
```

引用 DNAKit 不能替代对实际算法、参数集、后端和参考数据库的引用。项目使用的方法、论文、数据库和网站统一列在[致谢与主要来源](acknowledgements.md#acknowledgements)中；结果中的 `Provenance`、`BackendInfo`、`ReferenceLibrary` 和 `RunManifest` 用于保留这些信息。
