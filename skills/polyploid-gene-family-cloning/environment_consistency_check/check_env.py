import os
import sys
import subprocess
import shutil

# 定义学术调色板与符号标识
COLOR_GREEN = "\033[92m"
COLOR_RED = "\033[91m"
COLOR_YELLOW = "\033[93m"
COLOR_CYAN = "\033[96m"
COLOR_RESET = "\033[0m"

TICK_OK = f"{COLOR_GREEN}[✔]{COLOR_RESET}"
CROSS_FAIL = f"{COLOR_RED}[✘]{COLOR_RESET}"
WARN_ATTN = f"{COLOR_YELLOW}[!]{COLOR_RESET}"

def print_header(title):
    print("\n" + "=" * 60)
    print(f"{COLOR_CYAN}{title.center(60)}{COLOR_RESET}")
    print("=" * 60)

def check_windows_pip():
    print_header("1. 正在检测 Windows 环境中的 Python 第三方库")
    
    required_libs = {
        "pandas": "用于结构化表格解析与数据清洗",
        "openpyxl": "用于生成 Excel 格式文件",
        "matplotlib": "用于出版级图表绘制 (Figure 1)",
        "Bio": "Biopython，用于序列操作与系统树进化解析",
        "snapgene_reader": "用于直接解析 .dna 载体文件 (可选)"
    }
    
    missing_any = False
    for lib, desc in required_libs.items():
        try:
            if lib == "Bio":
                import Bio
            elif lib == "snapgene_reader":
                import snapgene_reader
            else:
                __import__(lib)
            print(f"  {TICK_OK} {lib:<16} -> 已就绪 ({desc})")
        except ImportError:
            if lib == "snapgene_reader":
                print(f"  {WARN_ATTN} {lib:<16} -> 未安装 ({desc}，若不存在可用其他库替代)")
            else:
                print(f"  {CROSS_FAIL} {lib:<16} -> 缺失! ({desc})")
                missing_any = True
                
    if missing_any:
        print(f"\n{WARN_ATTN} 修复提示：请在 Windows 终端中运行以下命令安装缺失的库：")
        print(f"  {COLOR_YELLOW}pip install pandas openpyxl matplotlib biopython snapgene_reader{COLOR_RESET}")
        return False
    return True

def run_wsl_cmd(cmd):
    """安全包装 WSL 命令调用"""
    try:
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "--"] + cmd if "Ubuntu-24.04" in get_wsl_distros() else ["wsl", "--"] + cmd,
            capture_output=True,
            text=True,
            shell=True
        )
        return result.returncode == 0, result.stdout.strip()
    except Exception:
        return False, ""

def get_wsl_distros():
    """获取 Windows 中注册的 WSL 列表"""
    try:
        res = subprocess.run(["wsl", "-l", "-q"], capture_output=True, text=True, errors='ignore')
        return res.stdout.split()
    except:
        return []

def check_wsl_environment():
    print_header("2. 正在检测 WSL2 (Ubuntu) 系统状态")
    
    # 检查 Windows 侧是否有 wsl 可执行文件
    if not shutil.which("wsl"):
        print(f"  {CROSS_FAIL} Windows 系统未找到 wsl 命令。请先启用 Windows 的 WSL2 功能。")
        return False
    print(f"  {TICK_OK} Windows WSL2 基础组件底层支持良好。")
    
    # 检查分发版
    distros = get_wsl_distros()
    if not distros:
        print(f"  {CROSS_FAIL} 未检测到任何已安装的 WSL 子系统。")
        return False
        
    print(f"  {TICK_OK} 检测到已安装的子系统: {', '.join(distros)}")
    
    # 测试能否成功唤醒 WSL 执行简单指令
    success, out = run_wsl_cmd(["uname", "-a"])
    if not success:
        print(f"  {CROSS_FAIL} 无法唤醒 WSL 或内部执行出错。请确保 WSL 子系统能够正常开机。")
        return False
        
    print(f"  {TICK_OK} WSL 桥接连通成功！内核信息: {out}")
    return True

def check_bio_software():
    print_header("3. 正在检测 WSL 内的生信工具链")
    
    tools = {
        "hmmsearch": "HMMER 3.x 家族成员鉴定工具",
        "mafft": "MAFFT 多序列比对工具",
        "fasttree": "FastTree 最大似然进化树构建工具"
    }
    
    missing_tools = []
    for tool, desc in tools.items():
        # 在 WSL 内部使用 which 检测
        success, _ = run_wsl_cmd(["which", tool])
        if success:
            print(f"  {TICK_OK} {tool:<12} -> 已部署 ({desc})")
        else:
            print(f"  {CROSS_FAIL} {tool:<12} -> 未找到! ({desc})")
            missing_tools.append(tool)
            
    if missing_tools:
        print(f"\n{WARN_ATTN} 修复提示：请进入 WSL2 (Ubuntu) 终端中运行以下命令安装缺失的生信软件：")
        print(f"  {COLOR_YELLOW}sudo apt-get update && sudo apt-get install -y hmmer mafft fasttree{COLOR_RESET}")
        return False
    return True

def check_workspace():
    print_header("4. 正在验证工作区与参考基因组数据结构")
    
    cwd = os.getcwd()
    print(f"  当前运行脚本的根目录: {cwd}")
    
    # 检查关键文件夹结构
    required_paths = [
        "ZM4_V2",
        "ZM4_HAP1",
        "Reference_Anchors"
    ]
    
    missing_path = False
    for path in required_paths:
        full_path = os.path.join(cwd, path)
        if os.path.exists(full_path) and os.path.isdir(full_path):
            print(f"  {TICK_OK} 目录检测通过: {path}/")
        else:
            print(f"  {CROSS_FAIL} 缺失核心目录: {path}/")
            missing_path = True
            
    if missing_path:
        print(f"\n{WARN_ATTN} 警告：请将此自检脚本放置在标准项目根目录下运行。")
        print("  工作区应包含 `ZM4_V2/`, `ZM4_HAP1/` 和 `Reference_Anchors/` 文件夹。")
        return False
    return True

def main():
    print(f"{COLOR_CYAN}=" * 60)
    print("多倍体基因家族分析及克隆自动化流水线 —— 环境自检引擎".center(48))
    print(f"={COLOR_RESET}" * 60)
    
    win_ok = check_windows_pip()
    wsl_ok = check_wsl_environment()
    
    bio_ok = False
    if wsl_ok:
        bio_ok = check_bio_software()
    else:
        print(f"\n{CROSS_FAIL} 由于 WSL 未就绪，跳过生信工具链的检测。")
        
    space_ok = check_workspace()
    
    print_header("自检报告最终结论")
    if win_ok and wsl_ok and bio_ok and space_ok:
        print(f"{COLOR_GREEN}恭喜！所有干、湿实验自动化环境组件全部就绪。{COLOR_RESET}")
        print("请直接将你的目标参数和引物设计方案输入 AI Agent，开启全自动分析流水线。")
    else:
        print(f"{COLOR_RED}检测到环境不完整。{COLOR_RESET}")
        print("请按照上方红色 {CROSS_FAIL} 或黄色 {WARN_ATTN} 的提示修复对应的模块，然后再运行 AI 自动化 Agent。")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()
