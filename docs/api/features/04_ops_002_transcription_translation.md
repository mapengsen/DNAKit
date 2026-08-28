# OPS-002 转录与翻译

将 DNA 序列转录为 RNA，或按照指定阅读框和密码表翻译为蛋白质序列。

## 1) OPS-002.1 转录

- **作用：** 按指定链方向将 DNA 中的胸腺嘧啶 `T` 转换为尿嘧啶 `U`，得到大写 RNA 序列，用于检查转录产物。
- **API：** `dnakit.ops.transcribe(sequence[必须], strand[可选])`。
- **输入：** 必填一条 `DNASequence`；`strand` 可选 `"forward"` 或 `"reverse"`，默认使用正向链。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import transcribe

seq = DNASequence("ATGGGCTAA")
print(transcribe(seq))
```

- **示例结果：**

```text
AUGGGCUAA
```

## 2) OPS-002.2 翻译

- **作用：** 按指定链方向、阅读框和遗传密码表将 DNA 或 RNA 翻译为蛋白质序列，用于检查编码区和氨基酸序列。
- **API：** `dnakit.ops.translate(sequence[必须], frame[可选], table[可选], strand[可选], stop_policy[可选], unknown_policy[可选], incomplete_policy[可选])`。
- **输入：** 必填 `DNASequence` 或已规范化的大写 DNA/RNA 字符串；可选 `strand`、`frame`、`table`、终止密码子策略、模糊密码子策略和不完整密码子策略。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import translate

seq = DNASequence("ATGGGCTAA")
print(translate(seq))
```

- **示例结果：**

```text
MG*
```
