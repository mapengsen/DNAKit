# DNAKit 本地性能基准

`benchmark_core.py` 测量原生核心操作的耗时和 Python 分配器峰值内存。报告包含完整参数、随机种子、输入摘要、Python/DNAKit/Rich 版本与平台信息。

结果只适用于当前机器与软件环境，不是跨机器性能声明。`tracemalloc` 不统计所有原生库内存。脚本有输入长度、重复次数、预热次数、案例数和总工作量硬上限，默认拒绝覆盖已有报告。

示例：

```bash
python -m benchmarks.benchmark_core \
  --sizes 100,1000,10000 \
  --repeats 3 \
  --warmups 1 \
  --seed 20260813 \
  --tasks construct,normalize,gc_content,kmer_fingerprint,minhash,subsequence_search \
  --output benchmarks/results/local_benchmark_report.json
```

