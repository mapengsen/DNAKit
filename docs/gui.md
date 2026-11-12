# 图形界面

DNAKit 的图形界面调用的函数与 Python API、`dnakit` 命令行和 MCP 服务器相同，用户不用写代码。
同一套代码提供两种形式：

- **桌面程序**：`dnakit gui` 在本机打开一个原生窗口。文件留在本地，运行不限时；会写文件、下载模型权重
  或调用外部程序的函数，经用户确认后执行。
- **网页版**：同一页面部署在服务器上，用浏览器就能试用。共享 CPU、单次调用限时，不执行有副作用的函数。

## 下载桌面程序

桌面程序可以直接下载使用，不需要安装 Python：

| 系统 | 下载 | 打开方式 |
|---|---|---|
| Windows 10 / 11（64 位） | [DNAKit-windows-x64.zip](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-windows-x64.zip) | 解压后进入 `DNAKit` 文件夹，双击 `DNAKit.exe` |
| macOS（Apple 芯片） | [DNAKit-macos-arm64.dmg](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-arm64.dmg) | 打开磁盘映像，把 **DNAKit** 拖进“应用程序”后打开 |
| macOS（Intel） | [DNAKit-macos-x64.dmg](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-x64.dmg) | 同上 |
| Linux x86-64（glibc 2.35 及以上，如 Ubuntu 22.04） | [DNAKit-linux-x64.tar.gz](https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-linux-x64.tar.gz) | `tar -xzf DNAKit-linux-x64.tar.gz && ./DNAKit/DNAKit` |

程序没有代码签名，第一次打开时系统会要求确认一次：

- **Windows：** 在“Windows 已保护你的电脑”对话框中点“更多信息”，再点“仍要运行”。
- **macOS：** 出现“无法验证 DNAKit”的提示后，打开“系统设置”→“隐私与安全性”，点“仍要打开”
  （macOS 14 及更早版本：右键点击程序，选“打开”）。也可以运行
  `xattr -dr com.apple.quarantine /Applications/DNAKit.app`。

下载版把上传的文件和结果保存在用户主目录下的 `DNAKit` 文件夹；双击启动时，运行信息写入其中的
`DNAKit/desktop.log`。下载版包含界面和 `pip install "dnakit[gui]"` 所带的函数。依赖可选软件包
（深度学习模型、`bio` 系列依赖组）或外部程序的函数会提示缺少什么，这些函数请使用下文的 pip 安装。
在终端中带参数启动时，该程序就是 `dnakit` 命令行，例如 `DNAKit describe ACGT`。

下载文件由 [Desktop app 工作流](https://github.com/mapengsen/DNAKit/actions/workflows/desktop.yml)
从 `main` 分支构建，构建方法见
[`desktop-app/`](https://github.com/mapengsen/DNAKit/tree/main/desktop-app)。

## 安装与启动

```bash
pip install "dnakit[gui]"
dnakit gui
```

窗口使用操作系统自带的网页组件：Windows 10/11 用 WebView2，macOS 用 WebKit。Linux 上如果系统装有
GTK WebKit，pywebview 会直接使用；否则安装 Qt 后端 `pip install "dnakit[gui-qt]"`，它的安装包要求较新的
Linux（glibc 2.28 及以上）和系统库 `libxcb-cursor0`。没有显示器或窗口组件时，`dnakit gui` 会在默认浏览器
中打开同一界面。

| 选项 | 作用 |
|---|---|
| `--workdir DIR` | 上传文件和结果的保存目录（默认 `./dnakit-gui`） |
| `--browser` | 用默认浏览器打开，不弹窗口 |
| `--no-open` | 只启动服务，例如在远程服务器上 |
| `--host`、`--port` | 服务地址和端口（默认 `127.0.0.1`，从 7860 起第一个空闲端口） |
| `--public-demo` | 使用下文网页版的限制 |

在没有桌面的远程服务器上，运行 `dnakit gui --no-open --port 7860`，再通过 SSH 端口转发
（`ssh -L 7860:127.0.0.1:7860 server`）或编辑器的端口转发在本机浏览器打开。

## 界面内容

| 标签页 | 内容 |
|---|---|
| Sequence (Example) | 反向互补链、240 个描述符、六框 ORF、IUPAC/正则 motif 搜索、序列比较、序列校验 |
| Dataset (Example) | 去重（完全相同或反向互补）、训练/验证/测试集划分、格式转换、每条记录的描述符特征表 |
| Chat (AI agent) | 用自然语言描述任务，由所选平台的大模型检索并运行 DNAKit 函数和 Skills（见[与大模型对话](#与大模型对话)） |
| Function explorer | 检索全部 Agent 可调用函数，跳转到本文档中对应的章节，查看参数说明，用 JSON 参数直接运行 |
| MCP tools | 检索第三方 MCP 服务器中的 DNA 与基因工具；点击某个工具可查看参数、源码仓库和 `call_mcp_tool` 调用写法，*docs* 链接打开它在 MCP 工具目录中的条目 |
| Skills | 检索和阅读内置 DNA 工作流 Skills；*docs* 链接打开它在 Skills 目录中的条目 |

每个结果都会给出等价的命令行和 Python 调用；函数检索页还给出 agent 发送的 MCP 请求，
在界面里试过的操作可以直接写成脚本。

## 与大模型对话

**Chat (AI agent)** 标签页提供对话式的 DNAKit：输入一句请求，例如“找出这条序列中不短于 60 nt 的全部
ORF，并翻译最长的一个”，大模型会自己规划步骤、调用 DNAKit，并根据运行结果回答。

1. 选择**平台**：DeepSeek、Grok（xAI）、通义千问（Qwen）、Kimi、OpenRouter，
   或任何兼容 OpenAI chat-completions 格式的接口（填写它的 base URL）。GPT、Gemini、Claude 等模型通过
   OpenRouter 使用。
2. 填入自己在该平台的 **API key**。填好后界面会从平台读取模型列表，从中选择一个模型，也可以直接输入
   模型名。
3. 输入请求（或点击一个示例），按 **Send**。

模型使用的就是 [DNAKit MCP 服务器](llm_tools/mcp.md)的工具：`search_dnakit_tools`、
`describe_dnakit_tool`、`call_dnakit_tool`、`search_dnakit_skills`、`read_dnakit_skill`，另有
`list_session_files` 用来查看用户添加的文件。每次工具调用连同参数和返回结果都显示在对话中，可以看到
答案出自哪些函数；函数返回的图直接显示为图片。在对话下方添加的文件保存在工作目录中，请求里写文件名
即可；目录中始终有一个 `example.fasta`（26 条记录）。

| 示例请求 | 模型的做法 |
|---|---|
| 求一条序列的反向互补链，并报告 GC 含量和熔解温度 | 调用序列与热力学函数 |
| 在六个读码框中找出不短于 60 nt 的 ORF，并翻译最长的一个 | 六框 ORF 扫描和翻译函数 |
| 查找 CpG 岛以及 EcoRI / BamHI / HindIII 酶切位点 | `dnakit.patterns.find_cpg_islands`、`dnakit.patterns.scan_restriction_sites` |
| 为目标序列设计 CRISPR-Cas9 向导 RNA | `dnakit.design_crispr.design` 及排序 |
| 统计 `example.fasta` 的记录数，找出互为反向互补的重复序列 | `dnakit.fasta_fastq.count_records`、`dnakit.datasets.deduplicate` |
| 哪些 Skills 与变异注释有关 | 检索 Skills 目录并概括其中一个 Skill |

使用须知：

- **key 和数据。** 请求由运行界面的那台电脑（网页版即服务器）发往所选平台，使用你的 key，费用由你的
  账户承担。key 只用于当前这条消息的请求，不保存，也不写入日志。对话内容（包括序列和工具返回的结果）
  会发送给该平台。
- **限制。** 工具调用遵守所在界面的规则（见下一节）。网页版只运行没有副作用的函数，单次调用限时
  60 秒，自定义接口地址必须是公网 HTTPS 地址。桌面程序中，勾选 *Allow tools with side effects* 后才
  运行会写文件、下载模型权重或调用外部程序的函数；也可以使用本机部署的模型（例如
  `http://localhost:11434/v1`），此时不需要 key。
- **用量。** 一条消息最多向模型请求 16 次，最后一次必须给出回答。输入框下方显示所用模型、工具调用次数、
  耗时，以及平台返回的 token 用量。对话过长时先丢弃较早的工具结果。
- **连通性。** 所选平台必须能从运行界面的那台电脑访问。在线网页版可以访问 DeepSeek、通义千问、Kimi 和
  OpenRouter，不能直接访问 Grok 的接口；在网页版使用 Grok 请通过 OpenRouter 或兼容 OpenAI 格式的中转地址
  （*Other OpenAI-compatible API*），或者使用桌面程序。

## 桌面程序与网页版的区别

| | 桌面程序（`GuiSettings.local`） | 网页版（`GuiSettings.public_demo`） |
|---|---|---|
| 文件 | 保存在工作目录 | 每位访问者一个临时目录，6 小时后删除 |
| 单次调用时间 | 不限（可用 **Stop** 按钮中止） | 60 秒 |
| 上传大小 | 不限 | 20 MB |
| 写文件、下载模型、调用外部程序的函数 | 确认后执行 | 不执行 |
| 对话：接口地址 | 任意 HTTP(S) 地址，包括本机部署的模型 | 预设平台，或公网 HTTPS 地址 |
| 路径参数 | 任意路径 | 只能指向访问者自己的目录 |

单序列工具最多接受 30,000 nt，因为序列通过命令行参数传递（Windows 命令行上限为 32,767 个字符）；
更长的序列请用 Dataset 标签页上传文件。

## Python 调用

```python
from dnakit.gui import GuiSettings, build_app, launch

launch(GuiSettings.local("results"))  # 打开桌面窗口，结果保存在 ./results
app = build_app(GuiSettings.public_demo())  # 部署到自己服务器用的 Gradio Blocks 应用
```

网页版所用的 systemd 和 nginx 配置见 `packaging/README.md`。
