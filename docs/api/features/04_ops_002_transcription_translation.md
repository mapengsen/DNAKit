# OPS-002 转录与翻译

将 DNA 序列转录为 RNA，或按照指定阅读框和密码表翻译为蛋白质序列。

- **作用：** 按链方向把 DNA 转录为 RNA，或按指定阅读框和遗传密码表翻译为蛋白质，用于检查转录产物、编码区和氨基酸序列。
- **API：** `dnakit.ops.transcribe(sequence[必须], strand[可选])`、`dnakit.ops.translate(sequence[必须], frame[可选], table[可选], strand[可选], stop_policy[可选], unknown_policy[可选], incomplete_policy[可选])`。
- **输入：** 转录必填 `DNASequence`；翻译可传 `DNASequence` 或已规范化 DNA/RNA 字符串；可选 `strand`、`frame`、`table`、终止/模糊/不完整密码子策略。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import transcribe, translate

seq = DNASequence("ATGGGCTAA")
print(transcribe(seq))
print(translate(seq))
```

- **示例结果：**

```text
AUGGGCUAA
MG*
```

- **限制：** 当前只支持 NCBI genetic code table 1 和 frame `0/1/2`；显式 Gap 会被拒绝。该功能只做序列规则转换，不预测表达或蛋白质功能。
