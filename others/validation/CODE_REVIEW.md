# 本地验证与 benchmark 代码审查

审查范围：`benchmarks/**`、`others/validation/**`、`others/tests/test_benchmarks.py`、`others/tests/validation/**`。

结论：未发现阻止本地使用的代码问题。输出使用同目录临时文件与原子链接/替换，默认拒绝覆盖；JSON 禁止 NaN；输入规模、重复、预热、案例数和总 nucleotide-runs 均有硬上限；随机输入使用显式 seed 和 SHA-256 摘要；脚本记录环境、参数、计时器和内存计量定义。

已核对的关键边界：

- benchmark 的 `tracemalloc` 只代表 Python allocator 峰值，报告中已禁止解释为进程总内存；
- Primer3 已从正式验证器移除；报告和单元测试确认验证器不自动发现、安装、导入或调用 Primer3。CLI adapter 只由临时受控假可执行文件验证命令白名单、解析、结果绑定和资源上限，不写成真实科学差分；
- Biopython restriction search 的 1-based cleavage boundary 在减一后才与 DNAKit 0-based boundary 对照；
- 分子量质量表和端基约定不完全相同，报告保留逐长度原始误差和显式磷酸化误差；
- NUPACK 没有安装、探测、导入或调用尝试；Primer3 同样没有自动发现、安装、导入或调用尝试；
- 未生成或声称任何论文复现实验。

剩余限制：microbenchmark 只适用于本机；可选后端升级后必须重新生成报告；
single/complete/average linkage merge distance 已与 Biopython `Bio.Cluster` 在相同
identity-distance matrix 上对照，最大绝对误差不超过 `1e-12`，但 CD-HIT、MMseqs2、
Dashing 等外部聚类工具没有在本验证器中执行。
