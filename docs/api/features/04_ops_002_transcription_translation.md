# 转录与翻译

将 DNA 序列转录为 RNA，或按照指定阅读框和密码表翻译为蛋白质序列。

## 1) 转录 {#1-ops-0021}

- **作用：** 按指定链方向将 DNA 中的胸腺嘧啶 `T` 转换为尿嘧啶 `U`，得到大写 RNA 序列，用于检查转录产物。
- **处理方法：** 正向链将 T 替换为 U；反向链先反向互补再转录。
- **API：** `dnakit.ops.transcribe(sequence[必须], strand[可选])`。
- **输入：** 必填一条 `DNASequence`；`strand` 可选 `"forward"` 或 `"reverse"`，默认使用正向链。
- **输出：** 返回转录后的 RNA 字符串。
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

## 2) 翻译 {#2-ops-0022}

- **作用：** 按指定链方向、阅读框和遗传密码表将 DNA 或 RNA 翻译为蛋白质序列，用于检查编码区和氨基酸序列。
- **处理方法：** 按阅读框将序列分为密码子，再查指定遗传密码表并应用终止和模糊密码子策略。
- **API：** `dnakit.ops.translate(sequence[必须], frame[可选], table[可选], strand[可选], stop_policy[可选], unknown_policy[可选], incomplete_policy[可选])`。
- **输入：** 必填 `DNASequence` 或已规范化的大写 DNA/RNA 字符串；可选 `strand`、`frame`、`table`、终止密码子策略、模糊密码子策略和不完整密码子策略。
- **输出：** 返回翻译后的氨基酸序列字符串。
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

<!-- BEGIN GENERATED API: translation -->

<span id="translation"></span>

## 3) 翻译三个正向阅读框并查找最长开放阅读框 {#dnakit.databases.dna_translate_reading_frames}

- **作用：** 翻译三个正向阅读框并查找最长开放阅读框。
- **处理方法：** 按标准遗传密码表翻译三个正向阅读框，比较无终止密码子的连续片段。
- **API：** `dnakit.databases.dna_translate_reading_frames(*, operation=None, sequence, frame=None, genetic_code=None, verbose=True)`
- **输入：** 必须提供 `sequence`。
- **输出：** 返回字典，包含 `sequence_length`, `frames`, `best_frame`, `best_orf` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dna_translate_reading_frames(
    operation='translate_reading_frames',
    sequence='ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "sequence_length": "<integer>",
  "frames": "<object>",
  "best_frame": "<integer>",
  "best_orf": "<string>",
  "best_orf_length_aa": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/04_ops_002_transcription_translation.md#dnakit.databases.dna_translate_reading_frames-parameters)

<!-- END GENERATED API: translation -->
