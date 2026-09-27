# 遗传模拟：参数与补充说明

用法和示例见[功能教程](../../api/features/34_inheritance.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 模拟独立位点的等位基因遗传 {#dnakit.workflows.simulate_inheritance-parameters}

`dnakit.workflows.simulate_inheritance(father, mother, trait_registry, disease_registry, *, offspring_count=3, generation=1, mutation_rate=0.001, seed=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `father` | dict&#91;str, Any&#93; | 必填 | 父方记录，sex 必须为 Male；基因型字段遵循 Genomebook 数据结构。 |
| `mother` | dict&#91;str, Any&#93; | 必填 | 母方记录，sex 必须为 Female；基因型字段遵循 Genomebook 数据结构。 |
| `trait_registry` | dict&#91;str, Any&#93; | 必填 | 显式提供的性状与位点定义。 |
| `disease_registry` | dict&#91;str, Any&#93; | 必填 | 显式提供的疾病位点定义。 |
| `offspring_count` | int | 3 | 模拟后代数量，1–10000。 |
| `generation` | int | 1 | 代次标签，正整数。 |
| `mutation_rate` | float | 0.001 | 每个位点的模拟突变率，范围 0–1。 |
| `seed` | int &#124; None | None | 局部随机数种子；不修改全局随机状态。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
