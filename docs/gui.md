# 图形界面

DNAKit 的图形界面调用的函数与 Python API、`dnakit` 命令行和 MCP 服务器相同，用户不用写代码。
同一套代码提供两种形式：

- **桌面程序**：`dnakit gui` 在本机打开一个原生窗口。文件留在本地，运行不限时；会写文件、下载模型权重
  或调用外部程序的函数，经用户确认后执行。
- **网页版**：同一页面部署在服务器上，用浏览器就能试用。共享 CPU、单次调用限时，不执行有副作用的函数。

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
| Sequence | 240 个描述符、六框 ORF、IUPAC/正则 motif 搜索、序列比较、序列校验 |
| Dataset | 去重（完全相同或反向互补）、训练/验证/测试集划分、格式转换、每条记录的描述符特征表 |
| Function explorer | 检索全部 Agent 可调用函数，查看参数说明，用 JSON 参数直接运行 |
| Skills | 检索和阅读内置 DNA 工作流 Skills |
| Use with agents | 把同一套函数接入 Claude Code、Codex 等 MCP 客户端 |

每个结果都会给出等价的命令行和 Python 调用；函数检索页还给出 agent 发送的 MCP 请求，
在界面里试过的操作可以直接写成脚本。

## 桌面程序与网页版的区别

| | 桌面程序（`GuiSettings.local`） | 网页版（`GuiSettings.public_demo`） |
|---|---|---|
| 文件 | 保存在工作目录 | 每位访问者一个临时目录，6 小时后删除 |
| 单次调用时间 | 不限（可用 **Stop** 按钮中止） | 60 秒 |
| 上传大小 | 不限 | 20 MB |
| 写文件、下载模型、调用外部程序的函数 | 确认后执行 | 不执行 |
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
