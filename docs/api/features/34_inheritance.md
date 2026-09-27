# 遗传模拟

<!-- BEGIN GENERATED API: inheritance -->

<span id="inheritance"></span>

## 1) 模拟独立位点的等位基因遗传 {#dnakit.workflows.simulate_inheritance}

- **作用：** 根据显式提供的亲本与位点定义模拟等位基因遗传和突变。
- **处理方法：** 逐独立位点从双亲各抽取一个等位基因，再按给定概率施加突变。
- **API：** `dnakit.workflows.simulate_inheritance(father, mother, trait_registry, disease_registry, *, offspring_count=3, generation=1, mutation_rate=0.001, seed=None, verbose=True)`
- **输入：** 必须提供 `father`, `mother`, `trait_registry`, `disease_registry`。 使用显式合成位点演示独立遗传；关闭突变，不访问远程数据。
- **输出：** 字典包含 offspring、seed、mutation_rate、scope、upstream；只模拟独立位点遗传，不模拟染色体交叉或真实临床结局。
- **示例代码：**

```python
import dnakit

example_father = {'id': 'father',
 'sex': 'Male',
 'loci': {'L1': {'chromosome': '1',
                 'position': 1,
                 'ref': 'A',
                 'alt': 'C',
                 'dominance': 'additive',
                 'effect_size': 1,
                 'alleles': ['A', 'A']}}}

example_mother = {'id': 'mother',
 'sex': 'Female',
 'loci': {'L1': {'chromosome': '1',
                 'position': 1,
                 'ref': 'A',
                 'alt': 'C',
                 'dominance': 'additive',
                 'effect_size': 1,
                 'alleles': ['C', 'C']}}}

example_trait_registry = {'traits': {'dosage': {'loci': [{'chromosome': '1',
                                 'position': 1,
                                 'ref': 'A',
                                 'alt': 'C',
                                 'dominance': 'additive',
                                 'effect_size': 1,
                                 'id': 'L1',
                                 'effect': 1}]}}}

result = dnakit.workflows.simulate_inheritance(
    father=example_father,
    mother=example_mother,
    trait_registry=example_trait_registry,
    disease_registry={'diseases': {}},
    offspring_count=2,
    mutation_rate=0,
    seed=13,
    verbose=False,
)
print(result)
```

- **示例结果：** 字典包含 offspring、seed、mutation_rate、scope、upstream；只模拟独立位点遗传，不模拟染色体交叉或真实临床结局。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "offspring": "<list>",
  "seed": "<int | None>",
  "mutation_rate": "<float>",
  "scope": "<str>",
  "upstream": "<str>"
}
```

[完整参数与使用说明](../../faq/api/34_inheritance.md#dnakit.workflows.simulate_inheritance-parameters)

<!-- END GENERATED API: inheritance -->
