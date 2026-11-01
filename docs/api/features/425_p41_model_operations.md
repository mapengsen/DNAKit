# P41 模型操作

Kipoi 0.8.6、genomic-benchmarks 1.0.0、PyTorch 2.14.1 的公开检查已通过。Kipoi 和 CNN 使用真实本地小模型与直接上游比较；Evo2 仅做契约测试。TensorFlow 导出和真实 Evo2 权重推理尚未运行，尚未进行独立验收。

## Kipoi 激活和输入梯度

**作用：** 提取中间层激活，或计算模型输入梯度。

**处理方法：** `kipoi.get_model(..., with_dataloader=False)` 加载模型，再调用 `predict_activation_on_batch` 或 `input_grad`。`input_dtype='float32'` 支持 JSON 浮点列表；设为 `None` 保留类型。`pre_nonlinearity` 是否可用取决于上游模型。

**API：**

- `dnakit.predictions.kipoi_p41_functions.kipoi_layer_activations(model, inputs, layer, *, source='dir', pre_nonlinearity=False, input_dtype='float32', verbose=True)`
- `dnakit.predictions.kipoi_p41_functions.kipoi_input_gradients(model, inputs, *, source='dir', filter_idx=None, avg_func=None, layer=None, final_layer=True, pre_nonlinearity=False, input_dtype='float32', verbose=True)`

**输入：** 可信的本地 Kipoi 模型目录、模型要求的数组或嵌套字典和层名。修改 source 可触发上游模型查询。

**输出：** JSON 激活或梯度数组，以及 Kipoi 版本。测试使用固定权重线性 PyTorchModel，对独热、均匀和全零三组输入与上游比较。

## Kipoi 数据加载批次

**作用：** 有界读取 raw、train 或 predict 批次。

**处理方法：** 调用原生批次迭代器；train 关闭无限循环。读取一批额外数据判断 truncated。

**API：** `dnakit.predictions.kipoi_p41_functions.kipoi_dataloader_batches(dataloader, *, source='dir', dataloader_args=None, mode='predict', batch_size=32, max_batches=1, verbose=True)`

**输入：** 可信的本地加载器目录和初始化参数；batch_size 为 1–10000，max_batches 为 1–100。

**输出：** batches、batch_count 和 truncated。测试涵盖 1、4、5 条记录和三种模式。

## BED 预测结果和逐碱基轨道

**作用：** 写出 BED6 加预测列，或 bedGraph/bigWig 逐碱基轨道。

**处理方法：** 使用原生 BedBatchWriter、BedGraphWriter、BigWigWriter。二维 bedGraph 值由上游按通道求和。

**API：**

- `dnakit.predictions.kipoi_p41_functions.kipoi_write_bed_predictions(ranges, predictions, output_path, *, header=True, verbose=True)`
- `dnakit.predictions.kipoi_p41_functions.kipoi_write_region_track(regions, values, output_path, *, output_format='bedgraph', chrom_sizes=None, verbose=True)`

**输入：** BED 并行数组 chr/start/end/id/strand，或 chr/start/end 区间与等长有限数值。bigWig 需要 chrom_sizes 和排序区间；输出文件必须不存在。

**输出：** 实际输出路径。文本与上游逐字比较，bigWig 按读取值比较。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.predictions.kipoi_p41_functions import kipoi_write_region_track
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / 'track.bedgraph'
    result = kipoi_write_region_track(
        [{'chr': 'chr1', 'start': 0, 'end': 2}],
        [[0.1, 0.2]], str(path), verbose=False)
    print(result['status'], result['region_count'])
    print(len(path.read_text().splitlines()))
```

**示例输出：**

```text
completed 1
2
```

## genomic-benchmarks CNN

**作用：** 训练、评估和预测本地 DNA token CNN，或导出 TensorFlow 初始化架构。

**处理方法：** PyTorch 调用原生 CNN.fit/test/forward；CPU 训练保存 state_dict，加载使用 weights_only=True。保留上游二分类 sigmoid 与 binary_cross_entropy_with_logits 的组合。TensorFlow 调用 get_basic_cnn_model_v0，导出 .keras，不训练。

**API：** 模块 `dnakit.predictions.genomic_benchmarks_p41_functions` 提供：

- `train_genomic_benchmark_cnn(token_ids, labels, config, output_checkpoint, *, epochs=1, batch_size=32, seed=0, verbose=True)`
- `evaluate_genomic_benchmark_cnn(token_ids, labels, config, checkpoint, *, batch_size=32, verbose=True)`
- `predict_genomic_benchmark_cnn(token_ids, config, checkpoint, *, verbose=True)`
- `build_genomic_benchmark_tf_cnn(number_of_classes, vocab_size, output_path, *, input_len=100, verbose=True)`

**输入：** config 包含 number_of_classes、vocab_size、embedding_dim、input_len；token 为明确词表索引。CNN 架构要求足够长的输入，示例采用 100 碱基；更短的参数范围不保证上游卷积可执行。

**输出：** 权重文件、上游 accuracy/F1/六位小数损失，或预测分数。固定种子的混合、同聚物和多类三组训练结果与上游比较；不代表生物学预测质量。TensorFlow 依赖缺失，实际导出尚未验证。

**示例代码：**

```python
import tempfile
from pathlib import Path
import torch
from dnakit.predictions.genomic_benchmarks_p41_functions import train_genomic_benchmark_cnn, predict_genomic_benchmark_cnn
torch.set_num_threads(1)
config = {'number_of_classes': 2, 'vocab_size': 5, 'embedding_dim': 2, 'input_len': 100}
tokens = [[0] * 100, [1] * 100, [2] * 100]
with tempfile.TemporaryDirectory() as folder:
    checkpoint = str(Path(folder) / 'cnn.pth')
    trained = train_genomic_benchmark_cnn(tokens, [0, 1, 0], config, checkpoint, batch_size=3, verbose=False)
    predicted = predict_genomic_benchmark_cnn(tokens, config, checkpoint, verbose=False)
    print(trained['status'], predicted['status'])
    print(len(predicted['scores']), len(predicted['scores'][0]))
```

**示例输出：**

```text
completed completed
3 1
```

## 本地 Evo2 推理

**作用：** 序列似然、前向 logits/嵌入、逐位置熵或困惑度。

**处理方法：** 显式使用本地权重，不下载。保留大小写。评分可平均正反链；逐位置指标先平均熵，再取指数得到困惑度。前向输出有 max_output_values 限制。

**API：** 模块 `dnakit.predictions.evo2_p41_functions` 提供：

- `score_evo2_sequences(sequences, checkpoint, *, model_name='evo2_7b', batch_size=1, prepend_bos=False, reduce_method='mean', average_reverse_complement=False, verbose=True)`
- `forward_evo2_sequences(sequences, checkpoint, *, model_name='evo2_7b', layer_names=None, prepend_bos=False, device='cuda:0', max_output_values=1000000, verbose=True)`
- `evo2_positional_uncertainty(sequence, checkpoint, *, model_name='evo2_7b', metric='entropy', reverse_complement=False, prepend_bos=False, device='cuda:0', verbose=True)`

**输入：** 2–100000 碱基 A/C/G/T/N 序列和本地 checkpoint，兼容的 Evo2/Vortex/CUDA。

**输出：** JSON 数组和版本。仅参数、调用协议、反链位置和输出上限有契约测试；未使用真实权重验证科学结果。

**示例代码：**

```python
from dnakit.predictions.evo2_p41_functions import score_evo2_sequences
from dnakit.exceptions import ConfigurationError
try:
    score_evo2_sequences(['AC-X'], 'unused.ckpt', verbose=False)
except ConfigurationError:
    print('invalid DNA rejected before checkpoint loading')
```

**示例输出：**

```text
invalid DNA rejected before checkpoint loading
```
