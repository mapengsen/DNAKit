# DNAKit 本地正确性验证

`run_validation.py` 生成机器可读 JSON，包含：

- 人工可核验的组成、GC、k-mer、重叠搜索和反向互补；
- 小型确定性聚类语义；
- 安装可用时，与 Biopython 对照限制酶切位点、分子量终端约定和全局比对得分；
- 安装可用时，将同一 identity 距离矩阵的 single、complete、average 层次聚类合并距离与 Biopython `Bio.Cluster` 对照。

当前容差与约定：

- 限制酶切位点和重叠搜索：坐标换算后整数精确相等；
- 全局比对：相同线性 gap 参数下最优得分绝对差不超过 `1e-12`；
- 层次聚类：无并列小样上按相同距离矩阵比较合并距离序列，绝对误差不超过 `1e-12`；
- 分子量：Biopython 使用更高精度平均核苷酸质量并声明 5' phosphate；DNAKit 使用取至两位小数的无水残基表。未磷酸化误差和显式磷酸化误差均逐长度保存，以 `abs((DNAKit未磷酸化−Biopython)+79.0) <= 1.0 Da` 验证当前近似模型，不宣称公式完全相同。
- 边界：IUPAC、空线性序列、200,000 nt 输入和跨原点环状限制酶位点均有人工期望值。

Primer3 与 NUPACK 均不在本验证器中：不自动发现、不安装、不导入、不调用。Primer3 CLI adapter 的协议与防御边界由临时受控替身单元测试覆盖，不能替代真实科学差分。此目录也不声称完成论文复现实验。

示例：

```bash
python -m others.validation.run_validation \
  --output others/validation/results/local_validation_report.json
```
