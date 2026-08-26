# 固定输入只读演示

<div class="demo-banner">
  本页是静态展示：没有 DNA 序列输入、文件上传、数据库查询或在线计算后端。
</div>

## 固定输入

| ID | 序列 |
| --- | --- |
| `seq-a` | `ACGTACGT` |
| `seq-b` | `ACGTACGT` |
| `seq-c` | `ACGTTCGT` |

## 人工可核验结果

<div class="demo-grid">
  <div class="demo-card">
    <h3>seq-a 基础结果</h3>
    <dl>
      <dt>符号长度</dt><dd>8</dd>
      <dt>GC 比例</dt><dd>0.5</dd>
      <dt>反向互补</dt><dd><code>ACGTACGT</code></dd>
    </dl>
  </div>
  <div class="demo-card">
    <h3>seq-a 高级确定性结果</h3>
    <dl>
      <dt>Wallace Tm</dt><dd>24.0 °C</dd>
      <dt>linguistic complexity（k=1…3）</dt><dd>0.38095238095238093</dd>
      <dt>exact tandem-repeat fraction</dt><dd>1.0</dd>
      <dt>透明 synthesis-risk level</dt><dd>low</dd>
    </dl>
  </div>
  <div class="demo-card">
    <h3>seq-a 的 2-mer 计数</h3>
    <dl>
      <dt>AC</dt><dd>2</dd>
      <dt>CG</dt><dd>2</dd>
      <dt>GT</dt><dd>2</dd>
      <dt>TA</dt><dd>1</dd>
    </dl>
  </div>
  <div class="demo-card">
    <h3>数据集检查</h3>
    <dl>
      <dt>记录数</dt><dd>3</dd>
      <dt>唯一原始序列数</dt><dd>2</dd>
      <dt>完全重复组</dt><dd>seq-a, seq-b</dd>
    </dl>
  </div>
</div>

基础组成、反向互补与 2-mer 值先人工计算；高级值由固定参数的当前本地 DNAKit API 生成。`tests/test_static_demo.py` 直接复核页面展示的长度、GC、反向互补、2-mer、去重、Wallace Tm、linguistic complexity、repeat fraction 和 synthesis-risk。synthesis-risk 是透明序列规则，不是实验成功率。页面仍保留为静态夹具，不会在浏览器中计算或接受 DNA 序列输入。机器可读副本见 [`fixed_demo.json`](data/fixed_demo.json)，源码夹具和复核命令见[示例页面](../examples/index.md)。

## 安全边界

- 页面不会读取访问者数据；
- 页面不会执行 NUPACK、Primer3、BLAST 或其他外部后端；
- MkDocs 已禁用站内搜索输入和外部 Web 字体；固定演示无运行时网络依赖；
- 所有展示内容在构建文档前已经固定。
