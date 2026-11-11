# DNAKit 桌面程序（可直接下载的版本）

这里是把 `dnakit gui` 打包成可下载程序的全部文件。用户下载后解压、双击即可使用，不需要安装 Python。

| 文件 | 作用 |
|---|---|
| `build.py` | 一条命令完成打包、自检和压缩，三个系统通用 |
| `dnakit_desktop.spec` | PyInstaller 打包配置（包含哪些模块和数据、图标、macOS 的 `.app`） |
| `dnakit_desktop.py` | 程序入口：双击时打开窗口；带参数时等同 `dnakit` 命令行，并代替 Python 解释器运行界面启动的子进程 |
| `smoke_test.py` | 在打包后的程序内部运行的自检（函数注册表、Skills、界面、子进程、进程池） |
| `RELEASE_NOTES.md` | 下载页面（GitHub `desktop` release）上显示的说明 |

## 更新下载文件

平时不需要手动打包。把改动推送到 GitHub 的 `main` 分支后，
[`Desktop app` 工作流](../.github/workflows/desktop.yml)会在 Windows、macOS（Apple 芯片和 Intel）
和 Linux 上各打包一次，运行自检，四个都通过后替换 `desktop` release 里的文件。README 中的下载按钮
指向这个 release 的固定地址，所以链接不用改。

- 改动涉及 `desktop-app/`、`src/dnakit/gui/`、`src/dnakit/cli/gui.py` 或版本号时，工作流自动运行。
- 其他改动（例如新增了函数）想更新下载文件时，在 GitHub 的 Actions 页面选择 `Desktop app`，
  点 `Run workflow`；或运行 `gh workflow run desktop.yml`。
- 一次运行约 15–20 分钟。任何一个系统失败都不会替换下载文件，旧文件保持可用。
- 每次运行的自检报告和窗口截图保存在该次运行的 `test-<系统>` 附件中。

下载地址：

```text
https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-windows-x64.zip
https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-arm64.dmg
https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-macos-x64.dmg
https://github.com/mapengsen/DNAKit/releases/download/desktop/DNAKit-linux-x64.tar.gz
```

## 在自己的电脑上打包

每个系统只能打出自己的包：Windows 的包要在 Windows 上打，macOS 的包要在 macOS 上打。
在仓库根目录、一个新建的 Python 3.12 虚拟环境里运行：

```bash
python -m venv build-env
source build-env/bin/activate        # Windows: build-env\Scripts\activate
python desktop-app/build.py --install
```

`--install` 会把当前源码和打包工具装进这个环境；之后重复打包可以去掉它。结果在 `desktop-app/dist/`：

| 系统 | 生成的文件 |
|---|---|
| Windows | `DNAKit-windows-x64.zip` |
| macOS | `DNAKit-macos-arm64.dmg` 或 `DNAKit-macos-x64.dmg` |
| Linux | `DNAKit-linux-x64.tar.gz` |

`build.py` 打包后会运行 `smoke_test.py` 并打开一次窗口检查页面能否显示；加 `--skip-tests` 可跳过。
Linux 上打包需要 Qt 的系统库和（没有显示器时）`xvfb`，软件包列表见工作流文件。

## 常见改动

- **让下载版带上更多可选依赖**（例如 `bio` 依赖组）：在 `build.py` 的 `install()` 里把依赖组加到安装
  命令中。安装包会变大；加入 `torch` 会超过 2 GB。
- **换图标**：替换 `docs/assets/images/DNAKit-icon.png`（`dnakit_desktop.spec` 中的 `ICON`）。
- **打包后缺少某个模块或数据文件**：在 `dnakit_desktop.spec` 的 `hiddenimports` 或 `datas` 中补上，
  并在 `smoke_test.py` 里加一项检查。
- **消除首次打开的系统提示**：需要购买代码签名证书（Apple 开发者账号、Windows 代码签名证书），
  目前没有签名。
