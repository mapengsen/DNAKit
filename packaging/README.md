# 平台打包适配

这里保存 DNAKit 在 PyPI 之外的平台适配源文件。它们不会改变
`src/dnakit/` 中的核心代码。

## PyPI 发布

将 `src/dnakit/_version.py`、中英文安装说明、引用和更新日志同步至新版本，
通过本地测试后推送 `main`，再推送与版本号一致的 `vX.Y.Z` 标签。
标签将触发 TestPyPI 流程：测试、Notebook、文档、后端验证、构建、
源码与 Skills 校验和审查、隔离安装检查，以及 TestPyPI 上传。
成功后自动调度现有正式发布流程；正式流程核对标签及 TestPyPI 上
同一组文件的 SHA-256，再上传到 PyPI 并建立附带发行包的 GitHub Release。

两个上传任务继续使用各自的 `testpypi` / `pypi` 环境和 Trusted Publishing。
手动调度仍保留原有的显式确认输入。

## Bioconda

[`bioconda/meta.yaml`](bioconda/meta.yaml) 是 `dnakit 0.1.1` 的 Bioconda
配方源稿，使用 PyPI 已发布的源码包和固定 SHA-256。

本地渲染：

```bash
conda-build packaging/bioconda --output
```

截至 2026-09-11，[Bioconda 已发布 DNAKit 0.1.4](https://anaconda.org/bioconda/dnakit)，
0.1.1 和 0.1.3 也仍可下载；
[0.1.4 更新 PR #68962](https://github.com/bioconda/bioconda-recipes/pull/68962)
已合并并通过六项检查。本目录仍保留 Galaxy 当前依赖的 0.1.1 配方。
现在可以使用：

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
截至 2026-09-11，PR 仍未合并；其中 `6edc2d2` 已通过
[Cuirass 官方 x86_64-linux 构建 #2196787](https://pulls.ci.guix.gnu.org/build/2196787/details)，
官方日志为 **1001 passed、15 skipped**。
[审核者也已成功构建并报告闭包大小为 495.0 MiB](https://codeberg.org/guix/guix/pulls/10903#issuecomment-22572178)。

根据审核意见，配方已移除重复的 `python-wheel` 直接构建依赖；
Guix 的 `python-setuptools` 已通过 `propagated-inputs` 提供它。
该调整已合入单个提交 `ddf6209` 并推送到原 PR。
删除重复声明后的提交仍需官方复验，不能将旧提交的成功结果当作新提交的结果。
本地已使用 setuptools 80.9.0 重新构建、安装干净的 `v0.1.1` 源码，并在
缺少 `/bin/sh` 的隔离环境中验证 **1001 passed、15 skipped**；依赖一致性和
CLI 冒烟检查通过。该环境没有单独安装 `wheel` 发行包，而是使用 PyPI
setuptools 提供的副本；这不等同于再次执行完整 Guix 构建。

2026-09-09 修复了官方 CI 的五项测试失败：测试运行时生成的脚本写死了
`#!/bin/sh`，而 Guix 构建沙箱没有 `/bin/sh`。配方在解包后将两个测试文件
中的脚本解释器替换为构建环境提供的 `sh`，并保留 pytest 的 `pythonpath=.`
设置以加载源码中的 `benchmarks` 和 `validation`；没有禁用测试。

本地使用 Python 3.12 和 Guix pytest 插件，在缺少 `/bin/sh` 的隔离环境中
先复现了同样的五项失败，再执行实际的 Guix `substitute*` 修复阶段，
完整测试结果为 **1001 passed、15 skipped**。已验证此独立配方与 PR 中的
包定义一致。

官方评估 #9292 中另外两个失败包为 `guile2.0-haunt` 和 `guile2.2-haunt`，
日志均报告配置阶段缺少 `(htmlprag)`；它们既不是本 PR 修改的包，也不是
DNAKit 的构建依赖。评估整体仍有失败，不应描述为所有 CI 均通过。

## Galaxy

Galaxy Tool Shed 包装器位于 [`galaxy/dnakit/`](../galaxy/dnakit/) 中。Galaxy
依赖 Bioconda 提供 `dnakit=0.1.1`，该版本现已发布。
2026-09-07 已发布到
[Main Tool Shed](https://toolshed.g2.bx.psu.edu/view/mapengsen/dnakit/9f6dc4f3006b)，
修订号为 `9f6dc4f3006b`，包含四个 `0.1.1+galaxy0` 工具。
2026-09-11 再次核对 Main Tool Shed：仓库公开、未弃用、修订可下载，
四个工具均有效，未报告缺失测试组件。Galaxy 仍固定依赖 Bioconda 0.1.1，
并未自动升级为 0.1.4。
本地 `planemo lint`、`planemo shed_lint` 及 Galaxy 24.2 的四项工具测试
全部通过。运行测试时将独立安装的 Bioconda 0.1.1 环境加入 `PATH`，
使用 `--no_dependency_resolution`；未测试 Galaxy 自动创建 Conda 环境。
验证和上传命令见 [`galaxy/dnakit/README.rst`](../galaxy/dnakit/README.rst)。

## 网页演示（Hugging Face Spaces 或自建服务器）

网页演示与桌面程序 `dnakit gui` 共用 [`src/dnakit/gui/`](../src/dnakit/gui/) 中的界面代码，
修改界面只改这一处。[`huggingface_space/app.py`](huggingface_space/app.py) 只是以
`GuiSettings.public_demo()` 启动同一界面：不执行写文件、下载模型或调用外部程序的函数，单次调用
限时 60 秒，路径参数只能指向访问者自己的临时目录。

`dnakit.gui` 从 0.2.0 之后的版本开始随 PyPI 发布。在此之前，`requirements.txt` 仍固定
`dnakit==0.2.0`，由 `deploy.sh` 把源码中的 `src/dnakit/gui` 复制到服务器；发布新版本后把
`requirements.txt` 改为 `dnakit[gui,bio-analysis,bio-formats,bio-plots]==<新版本>`。

本地预览（与服务器相同的限制）：

```bash
python -m pip install -e ".[gui]"
dnakit gui --public-demo --browser
```

### Hugging Face Spaces

Hugging Face 现在要求 PRO 订阅才能在 cpu-basic 上运行 Gradio Space。有订阅时：

```bash
pip install -U huggingface_hub
hf auth login
hf repo create <用户名>/DNAKit --repo-type space --space-sdk gradio
hf upload <用户名>/DNAKit packaging/huggingface_space . --repo-type space
```

### 自建服务器（systemd + nginx）

[`dnakit-web.service`](huggingface_space/dnakit-web.service) 以无登录权限的 `dnakit` 用户运行网页，
程序在 `/opt/dnakit-web/app`，虚拟环境在 `/opt/dnakit-web/venv`，限制为 2 核 CPU、1 GB 内存，
系统目录只读、`/tmp` 私有，只监听 `127.0.0.1:7860`。前面由 nginx
（[`nginx-dnakit.conf`](huggingface_space/nginx-dnakit.conf)）监听 80 端口：HTML 到达后立即显示
加载页，Gradio 前端文件直接从磁盘发送并缓存一年，文本内容 gzip 压缩。

首次安装（在服务器上，以 root 执行）：

```bash
apt-get install -y python3-venv nginx
useradd --system --home-dir /var/lib/dnakit-web --shell /usr/sbin/nologin dnakit
mkdir -p /opt/dnakit-web/app
python3 -m venv /opt/dnakit-web/venv
/opt/dnakit-web/venv/bin/pip install "dnakit[bio-analysis,bio-formats,bio-plots]==0.2.0" "gradio==6.29.1"
ln -sfn "$(/opt/dnakit-web/venv/bin/python -c 'import gradio, pathlib; print(pathlib.Path(gradio.__file__).parent / "templates/frontend")')" /opt/dnakit-web/frontend
rm -f /etc/nginx/sites-enabled/default
ln -sf /etc/nginx/sites-available/dnakit /etc/nginx/sites-enabled/dnakit
systemctl enable nginx dnakit-web  # 服务文件由下面的 deploy.sh 写入
```

之后每次更新界面，在本机源码目录执行一条命令即可（复制界面代码和服务文件、重启并检查）：

```bash
packaging/huggingface_space/deploy.sh root@<服务器>
```

升级 gradio 后重新执行上面的 `ln -sfn` 一行，使 nginx 发送新版前端文件。在中国大陆服务器上用域名对外
提供网站需先完成 ICP 备案。
