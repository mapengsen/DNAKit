# 正确性验证与 benchmark

## 当前结果

本地机器可读工件：

- `validation/results/local_validation_report.json`：生成于 `2026-08-15`；15 checks，15 pass，0 fail，0 not-run，0 not-comparable；
- `benchmarks/results/local_benchmark_report.json`：生成于 `2026-08-13T11:36:49Z`；2 个输入规模，7 个 DNAKit 任务、3 个有 Biopython 对等实现的任务，共 20 个 case。

环境快照：Python 3.10.20、DNAKit 0.1.0.dev0、Biopython 1.88、Linux/WSL2。Primer3 不是 Python 依赖，验证器明确记录其未被自动发现或执行。

!!! warning "结论范围"
    验证只覆盖表中明确列出的比较。它不证明所有输入、所有后端或实验正确性，也不包含真实 Primer3、NUPACK、CD-HIT、MMseqs2、BLAST、Dashing 或论文复现实验。

## 验证项目与允许误差

| ID/类别 | 对照 | 通过条件 | 当前结果 |
| --- | --- | --- | --- |
| `MANUAL-001..005` | 人工计数 | A/C/G/T、GC、overlap k-mer、overlap search、reverse complement 精确相等 | pass |
| `BOUNDARY-001` | 人工期望 | 空线性序列长度/坐标跨度为 0；空环状拒绝 | pass |
| `BOUNDARY-002` | 人工期望 | 完整 IUPAC 字母表长度与 ambiguity 数精确相等 | pass |
| `BOUNDARY-003` | 人工期望 | 200,000 nt 序列长度与 GC=0.5 精确相等 | pass |
| `BOUNDARY-004` | 人工期望 | 环状 EcoRI 跨原点切割坐标精确相等 | pass |
| `ALGORITHM-001` | 人工小图 | threshold graph cluster membership/order 精确相等 | pass |
| `BIOPYTHON-001` | Biopython Restriction | 统一为 0-based cleavage boundary 后整数精确相等 | pass |
| `BIOPYTHON-002` | Biopython molecular weight | `abs((DNAKit未磷酸化−Biopython)+79.0) <= 1.0 Da` | pass |
| `BIOPYTHON-003` | Biopython PairwiseAligner | 同一 global linear-gap 参数下最优分数绝对差 `<=1e-12` | pass |
| `BIOPYTHON-004` | Biopython Seq.search | overlap literal search 0-based start 精确相等 | pass |
| `BIOCLUSTER-001` | Biopython Bio.Cluster | 同一 DNAKit global-alignment identity distance 矩阵下 single/complete/average linkage 合并距离逐项绝对差 `<=1e-12` | pass |

### 分子量约定

Biopython 使用更高精度平均 nucleotide weights，并按其接口约定包含 5′ phosphate。DNAKit 当前使用两位小数的无水残基表、`-61.96 Da` 羟基端修正和可选 `+79.0 Da` 磷酸近似。报告逐长度保存原始误差和显式磷酸化误差；`1.0 Da` 容差验证的是文档化近似模型，不宣称公式完全相同。

### Primer3 审计 {#primer3}

当前正式验证器没有 `PRIMER3-*` 项，也不搜索、安装、导入或调用 Primer3。`Primer3CLIAdapter` 和 `Primer3CLIDesignAdapter` 由单元测试临时创建受控假可执行文件，验证参数白名单、`oligotm`/`ntthal` 文本解析、Boulder-IO、序列/条件/候选绑定、失败、超时及输出上限。这些测试只证明接口契约，不能写成真实 Primer3 数值验证。

用户若在自己的许可环境中做真实科学对照，应记录 Primer3 版本、可执行文件 checksum、thermodynamic parameter 目录、完整条件与容差，并单独保存不受再分发限制的验证摘要。

### NUPACK 审计

当前报告的禁止行为审计明确记录：验证器没有安装、探测、导入或调用 NUPACK；同一审计也记录没有自动发现、安装、导入或调用 Primer3。这些布尔值只描述验证器行为，不是对整台机器的扫描。`NupackAdapter` 的字段与异常边界由受控替身覆盖；当前环境没有真实 MFE、配对概率或 tube concentration 数值差分结论。

## 运行验证

默认拒绝覆盖已有报告。要生成新报告，应先选择新的输出路径；若脚本版本支持显式覆盖，再按 `--help` 使用对应选项。

```bash
PYTHONNOUSERSITE=1 python -m validation.run_validation \
  --output validation/results/local_validation_report.refresh.json
```

查看参数：

```bash
PYTHONNOUSERSITE=1 python -m validation.run_validation --help
```

## 本地 microbenchmark

固定参数：seed `20260813`，size 100/1000 nt，1 次重复，1 次预热。DNAKit 运行 construct、normalize、GC、k-mer fingerprint、MinHash、subsequence search 和 reverse complement；Biopython 1.88 只运行存在直接公开对等入口的 construct、GC 和 reverse complement。其余任务明确不伪造对等比较。

| 实现/任务 | 100 nt median | 1000 nt median | 1000 nt Python allocator peak |
| --- | ---: | ---: | ---: |
| DNAKit construct | 0.089710 ms | 0.482677 ms | 592 B |
| Biopython construct | 0.043011 ms | 0.014829 ms | 1,181 B |
| DNAKit normalize | 4.373388 ms | 25.123024 ms | 358,508 B |
| DNAKit GC content | 0.158554 ms | 0.265315 ms | 1,516 B |
| Biopython GC content | 0.040655 ms | 0.030125 ms | 396 B |
| DNAKit k-mer fingerprint | 3.408180 ms | 6.061041 ms | 42,842 B |
| DNAKit MinHash | 0.899842 ms | 9.534002 ms | 84,604 B |
| DNAKit subsequence search | 0.061717 ms | 0.137280 ms | 1,692 B |
| DNAKit reverse complement | 0.023589 ms | 0.439863 ms | 2,834 B |
| Biopython reverse complement | 0.007116 ms | 0.004219 ms | 2,154 B |

报告还按统一 `inspect.getsourcelines` 口径统计所选公开 callable 的非空、非注释行：DNAKit 637 行、Biopython 241 行。该值只说明本次选择的 callable body，不含传递依赖，也不是质量、可维护性或工作量指标。

这些数值来自一次本机短跑，重复数很少；1000 nt 个别任务比 100 nt 快是计时噪声示例，不能解释为反向 scaling。`tracemalloc` 只统计 Python allocator，不是进程 RSS 或原生库总内存；本报告不支持跨机器性能排名。

运行一个新的报告：

```bash
PYTHONNOUSERSITE=1 python -m benchmarks.benchmark_core \
  --sizes 100,1000,10000 \
  --repeats 3 \
  --warmups 1 \
  --seed 20260813 \
  --tasks construct,normalize,gc_content,kmer_fingerprint,minhash,subsequence_search \
  --implementations dnakit,biopython \
  --output benchmarks/results/local_benchmark_report.refresh.json
```

## 尚未完成的验证

- NUPACK adapter 契约已测试；真实二级结构/配对概率/tube 数值差分仍受独立许可和安装条件阻断；
- 光学/浓度、双链平衡、dot-bracket、条件 NUPACK adapter、PDB 几何及 3DNA/DSSR 解析已纳入全量回归；当前结果为 `945 passed, 1 skipped`；
- Biopython linkage 已对照；Dashing adapter 已完成固定协议、解析、错误与资源边界契约测试，并通过本地 Dashing `v1.0.2-4-g0635` 两序列 exact 文档示例 smoke，但尚无科学差分；CD-HIT/MMseqs2/BLAST/sourmash 的大库搜索/聚类差分仍因执行后端或数据库未集成而缺失；
- Parquet 已完成 PyArrow 25.0.1 下的 DNAKit 原子写出与 `read_table()` 往返；尚未做其他 Parquet 引擎/版本的交叉矩阵；
- 多机器、多线程 scaling、进程 RSS 和除已配置 Biopython 对等任务外的外部工具性能排名；源码行口径已记录，但不作为质量指标；
- 论文基准数据、消融、统计检验和论文复现实验。

对应状态在追踪矩阵中保持 `partial`/`conditional`，真实后端对照不会由替身单元测试通过替代。
