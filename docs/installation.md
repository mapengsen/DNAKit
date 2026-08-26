# 安装

## 系统范围

- Python 3.10 或更高版本；
- 当前认证环境：Linux；

## 创建完整开发环境

在项目根目录执行：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
export PYTHONNOUSERSITE=1
python -m pip install -e ".[dev,docs,io,validation,viz,neural]"
```

`PYTHONNOUSERSITE=1` 用于避免用户级 `site-packages` 污染测试。当前已验证环境使用 Python 3.10.20。

## 可选后端

### Primer3

DNAKit 不依赖或安装 `primer3-py`。用户须从 Primer3 官方渠道单独取得 CLI，并把 `oligotm`、`ntthal` 或 `primer3_core` 的绝对/显式相对路径传给 adapter；普通命令名和 `PATH` 自动发现会被拒绝。构造 adapter 只检查路径，不执行版本探测；只有显式调用计算方法才启动有超时和输出上限的子进程，且不经过 shell。

`Primer3CLIAdapter` 提供 Tm、hairpin、self-dimer 和 heterodimer；`Primer3CLIDesignAdapter` 执行 Boulder-IO 引物设计。结果记录 `GPL-2.0-or-later` 提示、实际路径、请求条件和 provenance，但使用者仍须以安装版本中的许可证为准。用户若提供独立 thermodynamic parameter 目录，也应同时确认其来源和分发权限。

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter

adapter = Primer3CLIAdapter(
    oligotm_path="/opt/primer3/src/oligotm",
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.hairpin(DNASequence("GCGTTTTTCGC"), output_structure=True)
print(result.tm_celsius, result.delta_g_kcal_per_mol)
```

### 可视化格式

`viz` extra 提供 CairoSVG/Pillow。SVG 原生生成不需要这些依赖；PNG、TIFF 和 PDF 转换才需要。像素数受 `ImageExportConfig.max_output_pixels` 限制。

### Parquet

`io` extra 安装 PyArrow。`export_table(..., format="parquet")` 和 `read_table(..., format="parquet")` 会先做被动探测；缺少后端时抛出结构化 `PARQUET_BACKEND_UNAVAILABLE`。当前环境使用 PyArrow 25.0.1 完成了 DNAKit 原子写出、DNAKit 读回和 schema metadata 复核。四种表格式都要求调用方提供显式 `TableSchema`；读取限制行数、列数、单元格、输入文件和解码字节，写出限制行数、列数、单元格和输出字节并在超限时原子回滚。

### DNA 基础模型 rep 与 k-means 聚类 {#dna-rep-k-means}

`neural` extra 安装 checkpoint 下载、Transformers/PyTorch 推理、NumPy 和
scikit-learn 聚类依赖：

```bash
python -m pip install -e ".[neural]"
```

Enformer、Caduceus 和 Evo 2 还需要各自的可选 extra：

```bash
python -m pip install -e ".[neural,neural-enformer]"
python -m pip install -e ".[neural,neural-caduceus]"
python -m pip install -e ".[neural,neural-evo2]"
```

Evo 2、Caduceus、AlphaGenome 和 JanusDNA 对 Python/CUDA/JAX、编译依赖或官方
源码环境有额外要求，不能保证与基础 `dnakit-dev` 环境一次解析成功。应按对应
官方仓库建立兼容环境后再安装 DNAKit。AlphaGenome 和 JanusDNA 的源码不由
DNAKit extra 安装；JanusDNA 需把官方源码路径传给 `model_source_path`。

模型 checkpoint 不随安装包分发。默认模型为 LucaOne，对应 checkpoint
`LucaGroup/LucaOne-gene-step36.8M`；首次显式调用表征 API 时，默认下载到当前
工作目录的 `ckpt/lucaone-gene-step36-8m/`，以后检测到完整文件便直接复用。
LucaOne 需要执行 checkpoint 自带 Python 代码，DNAKit 不会隐式放开此权限；
调用方审查来源后须显式设置 `allow_remote_code=True`。模型列表、checkpoint 来源和完整示例见
[序列表征：神经网络表征](api/features/08_fingerprints.md#neural-representations)。

### NUPACK 与外部搜索工具

DNAKit 提供 `probe_nupack()` 和 `NupackAdapter`，但 NUPACK 采用单独下载/许可流程，不在 extras 中，也不会由项目自动安装或下载。用户在适用许可范围内独立安装后，才可显式调用二级结构和 tube analysis；当前开发环境没有 NUPACK。BLAST、MMseqs2 和 sourmash 已注册被动路径探测及显式版本查询句柄，但没有序列搜索、聚类或相似度执行器。

`DashingAdapter` 是独立的显式科学计算 adapter：调用方必须提供一个现存、可执行的 Dashing 文件路径；DNAKit 不从 `PATH` 自动发现、不安装，也不选择项目根下的第三方副本。`matrix()` 可运行 exact k-mer set 或 HLL-sketch Jaccard，`top_k()` 对已验证矩阵做稳定排序。adapter 固定子命令/flag 白名单且不经 shell，限制输入项、输入/输出/捕获字节、sketch 内存、线程和超时，并记录输入/原始输出 SHA-256、后端版本、GPL 标记和 provenance。受控协议/失败/安全测试及本地 Dashing `v1.0.2-4-g0635` 的两序列 exact 文档示例 smoke 已通过；它不是科学差分，因此状态保持 `conditional`。

任何外部程序和数据库都不会自动下载或随包分发。NUPACK 可能需要付费订阅，DSSR 的学术与商业授权也不同；不能把“单独安装”理解为“免费或当然获准使用”。完整边界见[第三方声明](acknowledgements.md#third-party-notices)和[条件与不可用功能](planning/04_conditional_and_unavailable_features.md)。
