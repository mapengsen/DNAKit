---
name: bioinfo-pitfalls
description: >-
  Bioinformatics troubleshooting: HWE, multiple testing, aDNA damage,
  VCF merge traps, Slurm, array QC. Use when a pipeline looks significant
  but wrong, a tool errors on macOS/HPC, or the user asks 踩坑 / 排错.
---

## DNAKit execution

Use `dnakit skills doctor bioinfo-pitfalls` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Bioinformatics pitfalls

When debugging, match the symptom first, then apply the fix. Add new entries in the same shape: 现象 / 原因 / 解决 / 预防 / 来源.

# 生信踩坑与排错汇总

> 全库「踩坑速查总表」——按领域组织（统计/生信工具/古DNA/集群/群体遗传/芯片/数据库），每条含现象、原因、解决、预防、来源。踩新坑按格式追加。
> 2026-07-10 从根目录 `0-pitfalls.md` 移入此处（内容属生信技术，非知识库运维）。

---

## 统计与方法论

### HWE 检验未做导致假阳性位点混入

- **现象**：GWAS 或群体遗传分析中出现了大量"显著"位点，但这些位点在校正后消失；或者观察到杂合度偏离预期的位点在高频出现。
- **原因**：未对 SNP 做哈代-温伯格平衡（HWE）检验。HWE 偏离的位点通常是测序错误、比对错误或基因型 calling 错误导致的假阳性。这些假阳性位点在群体中表现出系统性的杂合子过量或缺失。
- **解决**：
  ```bash
  # PLINK 方法
  plink --bfile input --hwe 1e-6 --make-bed --out hwe_filtered
  # 或 vcftools 方法
  vcftools --vcf input.vcf --hwe 0.001 --recode --out hwe_filtered
  ```
  标准阈值：PLINK 用 `--hwe 1e-6`（严格）或 `1e-4`（宽松）。
- **预防**：在变异过滤流水线中，HWE 过滤应放在硬过滤之后、MAF 过滤之前。所有群体遗传分析前的 SNP 集都应通过 HWE 过滤。
- **来源**：GATK Best Practices

---

### 多重检验未校正导致 p 值被严重低估

- **现象**：基因组范围内的关联分析/差异表达分析中，出现了成千上万个"显著"结果（p < 0.05），但其中绝大多数是假阳性。结果无法在独立数据集中复现。
- **原因**：同时对大量假设做检验（如 20000 个基因 × 100 个代谢物 = 200 万次检验），如果单次检验的 alpha = 0.05，则期望有 10 万个假阳性。未做 Bonferroni 或 Benjamini-Hochberg (BH) 多重检验校正。
- **解决**：
  ```r
  # R 中 BH 校正
  data$p_adj <- p.adjust(data$p.value, method = "BH")
  # Dsuite 常用此校正
  ```
  ```python
  # Python 中
  from statsmodels.stats.multitest import multipletests
  reject, p_adj, _, _ = multipletests(p_values, method='fdr_bh')
  ```
  - Bonferroni：非常保守，适合严格发表要求
  - BH (FDR)：较为平衡，推荐用于探索性分析
  - 在多组比较时（如所有群体对的 Fst），也需要校正
- **预防**：在分析计划阶段就确定多重检验校正策略。报告时同时报告原始 p 值和校正后 p 值。Dsuite 用 `adjustP.r` 做 BH 校正，TreeMix p 值也应校正。
- **来源**：`notes/r/r-patterns-ggplot2.md`（adjustP.r），`notes/stats/mixed-effect-model-python.md`

---

### 混合效应模型只有 2 个分组水平导致随机效应估计不可靠

- **现象**：LMM 拟合后，随机效应的方差分量估计值极端（要么趋近 0，要么非常大），且标准误异常大。模型 summary 中的 `Group Var` 不可信。
- **原因**：混合效应模型中的随机效应需要足够的组数来估计组间方差。只有 2 个分组水平时（如"处理组 vs 对照组"），随机截距的方差本质上只有 1 个自由度，估计极不可靠。这不是代码 bug，是方法学限制。
- **解决**：
  - 如果只有 2 个分组：放弃随机效应，改用固定效应模型（或简单的 t 检验）
  - 如果分组 >= 5：LMM 开始有意义
  - 如果分组 3-4：谨慎使用，报告随机效应方差估计的不确定性
- **预防**：在设计实验时就确保每个随机效应至少有 5 个水平。如果数据已经收集完毕、只有 2 个分组——认命，不要强行用 LMM。
- **来源**：`notes/stats/mixed-effect-model-python.md`（奶牛代谢物案例——`Group` 字段只有 2 个水平，`ALL_UFA ~ 1 + (1 | Group)` 为空模型 + 不可靠的随机效应）

---

### Pearson 相关阈值 0.1 无生物学意义

- **现象**：代谢物-转录组关联分析中，筛选出大量 `|r| > 0.1` 的候选基因，但生物学验证发现几乎全是噪音。后续分析无法复现任何候选基因。
- **原因**：Pearson 相关系数 0.1 在生物学上代表极弱的相关（决定系数 R² = 0.01，即仅解释 1% 的方差）。即使统计显著（p < 0.05），生物学意义也极弱。18 个样本的转录组-代谢组配对中，|r| > 0.1 可能大部分是随机波动。
- **解决**：
  - 将阈值提高到 `|r| > 0.6` 或 `|r| > 0.7`（至少中等强度相关）
  - 同时要求 p 值校正后显著（BH 校正）
  - 用 Spearman 替代 Pearson（对异常值更稳健）
  - 用 WGCNA 做共表达模块分析而非单个基因-代谢物配对
- **预防**：任何"全转录组 × 代谢组"的大规模关联分析都应设合理的效应量阈值（r、fold-change），不能仅靠 p 值。Pearson r 的"大效应"阈值通常是 0.5-0.7。
- **来源**：`ideas/melon-metabolome-correlation.md`（甜瓜关联分析脚本中 `abs(correlation) > 0.1`），Cohen 1988 效应量标准

---

## 生信工具

### EDTA Helitron 假阳性 30-50%

- **现象**：EDTA 输出的 Helitron 序列中，30-50% 缺少 Helitron 的保守末端结构特征（5' TC 和 3' CTRR）。这些假阳性被当作真正的 Helitron 放入金标准库，污染了下游的全基因组注释。
- **原因**：EDTA 的 Helitron 鉴定算法基于结构特征而非序列保守性。在没有明显 Helitron 末端信号的区域，EDTA 会过度预测。尤其在 GC 含量不均或重复密集区，假阳性率飙升。
- **解决**：
  ```bash
  # 提取 Helitron 序列末端 10bp 检查保守性
  seqkit subseq -r 1:10 helitron_seqs.fa | grep -v '^>' | sort | uniq -c | sort -rn
  seqkit subseq -r -10:-1 helitron_seqs.fa | grep -v '^>' | sort | uniq -c | sort -rn
  ```
  移除缺少 5' TC 和 3' CTRR 的序列。在金标准库构建前必须做这一步。
- **预防**：EDTA 的 Helitron 输出默认不可信。每次构建金标准库前都手工核查 Helitron 保守末端。质检时关注 Helitron/LTR 比是否 > 30%。
- **来源**：`notes/bioinfo/annotation/te-annotation-pipeline.md`（EDTA 十大缺陷第 2 条）

---

### EDTA LINE/SINE 灵敏度极低

- **现象**：EDTA 注释的 LINE 和 SINE 占基因组比例远低于预期（如植物基因组 LINE 通常 1-5%，但 EDTA 只检出 < 0.5%）。这静默导致全基因组 TE 占比被低估。
- **原因**：EDTA 的核心算法（LTR 结构识别 + TIR 重复识别）对 nonLTR 元件（LINE/SINE，无末端重复结构）不敏感。LINE/SINE 的复制机制（靶向引发逆转录，TPRT）不产生可被 EDTA 识别为 TE 的结构特征。
- **解决**：加跑 RepeatModeler（对 nonLTR 敏感），将 RepeatModeler 输出中分类为 LINE/SINE 的序列补充进金标准库。或使用 NeuralTE（深度学习方法的 TE 发现工具）。
- **预防**：在 EDTA 运行前须知：对于 nonLTR 主导的基因组（如人类，LINE 占 ~21%），仅靠 EDTA 远远不够。质检时检查 LINE 占比是否过低。
- **来源**：`notes/bioinfo/annotation/te-annotation-pipeline.md`（EDTA 十大缺陷第 1 条）

---

### GATK `--min-meanDP 10` 忘记过滤导致低深度位点混入

- **现象**：群体遗传分析中发现大量"稀有变异"（MAF < 0.01），但 Sanger 验证发现约半数位点为假阳性。PCA 中低深度样本出现异常的离群位置。
- **原因**：GATK 的硬过滤不能完全排除低深度导致的假阳性。平均深度不足的位点，基因型 calling 的不确定性高，容易系统性地错误 calling 为杂合子。这些假位点在 PCA 中形成人为的"群体分化"信号。
- **解决**：
  ```bash
  vcftools --vcf input.vcf --min-meanDP 10 --recode --out dp_filtered
  ```
  这个过滤必须在 MAF 过滤之前执行——低深度位点的假阳性会抬高 MAF 分布的低端。
- **预防**：VCF 过滤流水线的标准顺序：硬过滤 → min-meanDP → MAF → Missing。min-meanDP 阈值根据平均测序深度设定（通常为平均深度的 1/3 到 1/2）。
- **来源**：`notes/bioinfo/ngs/gatk-variant-calling.md`（步骤 5：二次过滤），GATK Best Practices

---

### Beagle 18 号染色体 OOM（需增大 Xmx）

- **现象**：Beagle 5 对第 18 号染色体做 phasing 时，进程被 `Killed`（OOM）。其他染色体都正常。
- **原因**：18 号染色体可能特别大（如人类 chr18 ~80Mb）或 LD 区块特别长（如着丝粒附近的长单倍型块），导致 HMM 状态空间爆炸。Beagle 默认的 Java 堆内存不足以处理。
- **解决**：增大 Beagle 的 `-Xmx` 参数——专门为 18 号染色体分配更大的堆内存：
  ```bash
  java -Xmx128g -jar beagle.22Jul22.46e.jar gt=chr18.vcf.gz out=chr18_phased
  ```
  如果 128GB 还不够，将大染色体拆分为小段分别 phasing，再合并。
- **预防**：在投递作业前检查各染色体的 SNP 数量和 LD 衰减速度，预测内存需求。对已知的大染色体/高 LD 区域，提前分配更多内存。
- **来源**：`ideas/finestructure-troubleshooting.md`，用户实际排错记录

---

### Purge_dups 6 参数手调错误导致过度去冗余

- **现象**：Purge_dups 去重后，BUSCO 单拷贝完整度 (S%) 从 95% 暴跌到 75%，而多拷贝完整度 (D%) 从 15% 降到 2%。"去重"变成了"去基因"。
- **原因**：手工确定的 6 个 cutoff 参数中，`transition`（单倍型峰与双倍型峰之间的低谷）设得过低，导致双倍型 contig 被分类为"冗余单倍型"而被错误删除。看图时误把单倍型峰的"肩膀"当成了 transition。
- **解决**：
  1. 重新绘制覆盖度直方图（`hist_plot.py`）
  2. 重新确定 transition——它应该是单倍型峰和双倍型峰之间**最深**的点，不是峰的边缘
  3. 如果无法确定明确的 transition（两个峰重叠严重），调低期望——只去明显的冗余
  4. 验证：去重前后各跑一次 BUSCO，D% 下降的同时 S% 不应下降超过 5%
- **预防**：Purge_dups 的 6 个参数确定需要在看懂覆盖度直方图的基础上进行，不能凭"感觉"。建议先在测试数据集上调参验证，再去批量处理全部样本。如果覆盖度直方图没有清晰的双峰——跳过 Purge_dups，强行去重只会破坏组装。
- **来源**：`notes/bioinfo/assembly/assembly-postprocessing.md`（Purge_dups 6 参数图解 + Vacer 葡萄实例 `"5 5 22 32 44 60"`）

---

### TGS-GapCloser FASTQ 转 FASTA 格式问题导致一个洞都补不上

- **现象**：TGS-GapCloser 运行完成，输出文件存在，但 gap 数量和使用 HiFi reads 补闭的 gap 数均为 0。所有 gap 原样保留。
- **原因**：`seqtk seq -A` 输出的 FASTA 可能不是 TGS-GapCloser 期望的标准格式。常见细节问题：(1) 每行碱基数不固定；(2) 序列名中出现特殊字符；(3) FASTA 文件的换行符不标准 (CRLF vs LF)；(4) 同时可能 `-K 100M` minimap2 批次大小参数不匹配。
- **解决**：
  ```bash
  # 用 seqkit 替代 seqtk，输出更标准
  seqkit fq2fa hifi_reads.fq.gz -o hifi_reads.fa

  # 验证 FASTA 格式
  grep -c '^>' hifi_reads.fa
  awk 'NR%4==2' hifi_reads.fa | awk '{print length}' | sort -n | head -5

  # 调小 -K 参数：从 100M 降至 10M 或 50M
  tgsgapcloser --scaff assembly.fa --reads hifi_reads.fa \
    --ne --thread 36 --tgstype pb \
    --minmap_arg '-x asm20 -K 10M' \
    --out filled.fa
  ```
  如果仍然不行：(1) HiFi reads 长度太短，无法跨越 gap；(2) gap 侧翼序列在 reads 中没有覆盖（高杂合区）。
- **预防**：在补洞前用 `samtools faidx hifi_reads.fa` 验证 FASTA 格式。记录 HiFi reads 的 N50——如果 N50 < 最小 gap 长度，补洞不可能成功。
- **来源**：`notes/bioinfo/assembly/telomere-completion-t2t.md`（排错小节），`00_pan&genome/ex_gap_for_pcr.md`

---

## 古DNA 特异性

### 双端 BAM 喂给 mapDamage 导致长度分布错误

- **现象**：mapDamage 的 `Length_plot.pdf` 中片段长度分布出现双峰——第一个峰在 35bp，第二个峰在 70bp。而 fastp 报告显示实际片段分布是单峰（~50bp）。
- **原因**：mapDamage 对双端（paired-end）BAM 无法正确计算片段长度。它会错误地将双端 read 各自独立计算长度，产生两个虚假的长度分布峰。这个 bug 在 mapDamage2 文档中有明确记载。
- **解决**：在运行 mapDamage 前，将双端 BAM 转换为单端——只保留第一个 mate：
  ```bash
  samtools view -b -f 64 -F 1 paired.bam > single_end.bam
  samtools index single_end.bam
  mapDamage -i single_end.bam -r ref.fa -d output_dir
  ```
  `-f 64` 保留第一个 mate，`-F 1` 过滤掉第二个 mate。
- **预防**：所有 mapDamage 输入都使用单端 BAM——无论是原始单端数据还是双端数据的第一个 mate。在用 AdapterRemoval `--collapse` 选项时输出已经是单端，不需要再转换。
- **来源**：`notes/bioinfo/ancient-dna/adna-damage-assessment.md`（双端 BAM 陷阱），mapDamage2 官方文档

---

### 损伤水平过低（< 0.01）导致 Bayesian 计算被禁用——这不是报错

- **现象**：mapDamage 输出中显示 "Bayesian estimation disabled"，并且 `Stats_out_MCMC_hist.pdf` 未生成。新人看到 "disabled" 字样以为是错误。
- **原因**：当损伤水平低于 0.01 时，mapDamage 内部算法自动禁用贝叶斯估计——因为信号太弱，贝叶斯模型无法可靠收敛。这是 mapDamage 的**设计行为**，不是 bug。损伤水平用 `5pCtoT_freq.txt` 中的末端替换率评估——如果接近背景水平（如 < 0.02），贝叶斯就被禁用。
- **解决**：什么都不用做。检查 `Fragmisincorporation_plot.pdf` 确认损伤模式。如果曲线显示有轻微但可辨别的 5' C→T 升高，损伤就是真实存在的——只是水平低到无法贝叶斯建模。
- **预防**：在解读 mapDamage 输出前阅读其文档，理解哪些输出是可选的。如果损伤水平极低（< 0.01）+ 片段长度 > 100bp——应怀疑样本可能不是真古DNA（见下一条）。
- **来源**：`notes/bioinfo/ancient-dna/adna-damage-assessment.md`（常见报错表），mapDamage2 文档

---

### DNA 损伤水平极低可能意味着不是真古DNA

- **现象**：mapDamage 损伤曲线平坦，`5pCtoT_freq.txt` 中末端替换率与背景水平无异（< 0.02）。同时片段长度中位数 > 120bp。研究者仍按"古DNA"处理数据。
- **原因**：真古DNA 的 defiing 特征就是 C→T 末端损伤升高（胞嘧啶脱氨基）。如果损伤水平极低且片段长度偏长，最可能的解释是（a）样本中有大量现代 DNA 污染，稀释了古DNA 信号；或（b）样本根本不是古DNA。保存条件极好的样本（如永久冻土）可能有低损伤，但其片段长度仍应 < 100bp。
- **解决**：
  1. 用 DamageProfiler 快速交叉验证（比 mapDamage 快，不需要 Bayesian 计算）
  2. 做 BLAST nt / Kraken2 验证物种来源
  3. 如果片段长度 > 100bp + 损伤平坦 → 考虑作为"现代样本"处理
  4. 如果片段长度 < 100bp + 损伤平坦 → 可能保存条件特殊，记录并标记为"待确认古DNA"
- **预防**：古DNA 认证三要素（损伤模式 + 片段长度 + 物种来源）必须全部满足，不能只看其中一项。如果三项中有两项不及格——按现代 DNA 标准流程处理。
- **来源**：`notes/bioinfo/ancient-dna/adna-damage-assessment.md`（古DNA 认证三要素），`notes/bioinfo/ngs/fastp-qc-interpretation.md`（片段长度解读）

---

## 集群与运维

### Slurm 作业 Pending 不跑 — QOS 限额/分区不对

- **现象**：`squeue` 显示作业状态为 `PENDING`，`Reason` 列为 `QOSMaxCpuPerUserLimit` 或 `AssocGrpCPUMinutesLimit`。作业排了几个小时甚至几天都不跑。
- **原因**：已提交的作业总数超过了你在该 QOS (Quality of Service) 下的 CPU 限额。也可能是选错了分区——某些分区有自己的 QOS 限制（如 GPU 分区限制更严）。还有可能是 `Reason=NONE` 但队列拥堵。
- **解决**：
  ```bash
  # 查看队列状态
  squeue -o "%.18i %.9P %.12j %.12u %.12T %.12M %.16l %.6D %R" -u $USER

  # 查看分区限制
  sacctmgr show qos format=name,maxcpuperuser,maxtresperuser

  # 取消一些低优先级 PENDING 作业腾出限额
  scancel -t PENDING -u $USER --partition=GPU

  # 改到空闲分区
  sbatch --partition=C032M0128G --qos=low job.sh
  ```
- **预防**：提交前用 `squeue -u $USER | grep PENDING | wc -l` 检查积压情况。熟悉你的集群分区和 QOS 限制。对大任务使用 `--array` 分批提交而非一次全部提交。
- **来源**：`notes/tools/hpc/slurm-usage-guide.md`（常见问题表）

---

### Conda 环境改名失败 — conda 不支持 rename

- **现象**：尝试 `conda rename old_name new_name` 时报错 `conda: command not found`（但 conda 本身可用）。或直接找不到 rename 子命令。
- **原因**：conda **从来就不支持 rename 操作**。这不是 bug，是设计如此——conda 环境内部有大量硬编码路径，直接重命名会导致环境损坏。网上有些帖子的 `conda rename` 实际来自第三方插件或错误写法。
- **解决**：用 clone + remove 两步操作：
  ```bash
  conda create --name new_name --clone old_name
  conda remove --name old_name --all
  ```
  验证新环境可正常使用后，再删除旧环境。
- **预防**：创建环境时就用最终名称。如果需要改名——clone + remove。不要尝试手动 `mv` 环境目录（会破坏硬链接和路径）。
- **来源**：`notes/tools/conda/conda-environment-management.md`（环境改名 = clone + remove old），Conda 官方文档

---

### VCFtools 临时文件报错 — TMPDIR 未设置

- **现象**：vcftools 运行中突然报错 `Error: Could not open temporary file.` 或 `Write error`，无法完成过滤/统计。
- **原因**：vcftools 需要临时目录存储中间数据。如果 `TMPDIR` 未设置或指向的空间不足，vcftools 会尝试使用系统默认的 `/tmp`——而 `/tmp` 可能空间有限或被挂载为 tmpfs 内存文件系统（重启后清空）。
- **解决**：
  ```bash
  export TMPDIR=/path/to/large/scratch
  mkdir -p $TMPDIR
  vcftools --vcf input.vcf --recode --out output
  ```
  [版本相关] 新版 vcftools 中 `--temp-prefix` 参数已被移除，只能用 `TMPDIR` 环境变量。
- **预防**：在 Slurm 作业脚本开头 export TMPDIR。每个样本的临时文件使用独立子目录避免冲突。
- **来源**：`notes/bioinfo/ngs/vcf-tools-cheatsheet.md`（临时文件问题）

---

### R 包 00LOCK — unlink 或 kill 占用进程

- **现象**：`install.packages("xxx")` 报错：`ERROR: failed to lock directory '/path/to/R/library' for modifying. Try removing '/path/to/R/library/00LOCK-xxx'`
- **原因**：上一次 R 包安装被异常中断（Ctrl+C、崩溃、网络中断），残留的 `00LOCK-xxx` 目录未被清理。R 通过这个锁文件/目录机制防止并发安装冲突。
- **解决**：
  ```bash
  # 方案 1：手动删除 LOCK 目录
  unlink("/path/to/library/00LOCK-xxx", recursive=TRUE)
  # 或在 shell 中
  rm -rf /path/to/R/library/00LOCK-*

  # 方案 2：如果删除失败（有 R 进程占用）
  lsof | grep 00LOCK       # 找占用进程
  kill <pid>               # 终止占用进程
  rm -rf /path/to/R/library/00LOCK-*
  ```
- **预防**：安装大 R 包时不要按 Ctrl+C。使用 `nohup R CMD INSTALL ... &` 在后台安装大型包。安装前确保磁盘空间充足。
- **来源**：`notes/r/r-patterns-ggplot2.md`（R 包手动安装与 00LOCK 排错）

---

## 群体遗传学

### fineSTRUCTURE 样本名以数字开头导致数据静默不一致

- **现象**：fineSTRUCTURE 报错 `Data 0 is inconsistent: 97!=96`，或者更糟糕——**静默**产生错误结果而不报错。数据对不齐但看起来"运行成功"。
- **原因**：fineSTRUCTURE 内部对样本 ID 的处理中，以数字开头的 ID（如 `001_sample`）会被某些子程序解析为纯数字，导致字符串排序和数值排序不一致。这种不一致只在一部分子步骤中触发，所以表面上看"程序跑完了"，但数据在内部已经错位。
- **解决**：
  ```bash
  # 方案 1：所有 ID 前加字母前缀
  awk '{print "A"$1}' ids.txt > fixed_ids.txt
  # 方案 2：用 seq 生成字母前缀 ID
  seq -f "A%.0f" 1 96 > fixed_ids.txt
  ```
  确保 VCF 中的样本名与 fineSTRUCTURE 输入中的样本名都以字母开头。
- **预防**：任何以数字开头的生物样本 ID 都是定时炸弹——不仅 fineSTRUCTURE，PLINK、vcftools 等工具也有类似问题。样本命名规范：必须以字母开头，仅含字母、数字、下划线。
- **来源**：`ideas/finestructure-troubleshooting.md`（fineSTRUCTURE 最关键的踩坑记录）

---

### TreeMix 树中只有一个群体 — 外群名称拼错

- **现象**：TreeMix 输出的树中只有一个节点，没有预期的多群体分化结构。或者所有群体被堆在一起不分枝。
- **原因**：`-root` 参数指定的外群名称与 `pop.cluster` 文件中的群体名不一致——大小写不对、多了/少了空格、或者群体在输入文件中根本不存在。TreeMix 找不到外群时不会报错，而是静默生成一棵无根（或全部堆积的）树。
- **解决**：
  ```bash
  # 检查 pop.cluster 文件中的群体名
  cut -f3 pop.cluster | sort -u

  # 确保 -root 参数与上面的某个群体名完全一致
  treemix -i input.treemix.gz -root Outgroup -o output
  # 注意：上面用 "Outgroup" 而 pop.cluster 中写的是 "outgroup" → 不匹配！
  ```
- **预防**：在运行 TreeMix 前 `grep` 确认 `-root` 参数值与 `pop.cluster` 中的某一行完全一致（含大小写）。用脚本自动校验而非人工检查。
- **来源**：`notes/bioinfo/population/treemix-gene-flow.md`（常见问题与排错）

---

### easySFS 输出浮点数 SFS — 需 `--dtype int`

- **现象**：easySFS 输出的折叠 SFS (Site Frequency Spectrum) 中每个频率 bin 的值是小数（如 15.7 个 SNP），而下游工具（如 FitCoal）期望整数计数值。导致 FitCoal 报错或输出异常。
- **原因**：easySFS 默认对低覆盖度位点使用比例投影（proportional projection），导致 SFS 条目变成浮点数。这对于某些下游工具是合法的（如 dadi），但 FitCoal 和部分推断方法要求整数。
- **解决**：在 easySFS 命令中加 `--dtype int` 标志：
  ```bash
  easySFS.py -i input.vcf -p pop.txt --proj 20 --dtype int -o output.sfs
  ```
- **预防**：在运行 easySFS 前确认下游工具对 SFS 格式的要求。如果下游工具文档没有明确说明，优先用 `--dtype int`（大多数工具接受整数）。
- **来源**：`ideas/finestructure-troubleshooting.md`（FitCoal SFS 浮点数问题）

---

## 芯片设计

### DP10 过滤后补 — 早期数据需重跑

- **现象**：芯片探针设计时使用了含 `--min-meanDP 10` 过滤的 VCF。但后来发现部分早期测序数据（2018-2019 年）平均深度仅 6-8x，被该过滤排除了大量真实位点。需用未做 DP 过滤的原始 VCF 重新提取。
- **原因**：芯片探针设计中的位点来源有多批数据，各批的测序深度不同。统一的 DP10 过滤对低深度批次过于严格，丢弃了在低深度批次中出现的真实变异——这些变异恰好是"稀有"等位基因（在不同群体中频率差异大的位点），对芯片区分品种最有价值。
- **解决**：
  1. 回退到 DP 过滤前的 VCF（或原始 gVCF）
  2. 对低深度批次单独设定更宽松的 DP 阈值（如 4x），或在低深度位点用基因型似然度方法重新 calling
  3. 重新提取 MAF 和位点信息
- **预防**：多批次数据的过滤策略不应一刀切。按批次统计深度分布后按批次设过滤阈值。芯片设计对位点多样性的需求优先于位点 calling 的置信度。
- **来源**：`ideas/chip-design-160k-panel.md`（芯片设计全流程，方法学讨论）

---

### MinimalMarkers 大样本内存溢出 — 需迂回方案

- **现象**：使用 MinimalMarkers 工具从数万个候选位点中选取最小标记集时，样本数 > 3000 直接导致内存溢出（> 256GB）。
- **原因**：MinimalMarkers 的算法需要构建全两两样本间的差异矩阵，复杂度和内存消耗随样本数平方增长（O(n²)）。3000+ 样本的差异矩阵在内存中需要 > 70GB。
- **解决**：
  1. 方案 A：分群体抽样（每个群体随机选 200-300 个代表样本），用子集做 MinimalMarkers 筛选，再在全集中验证
  2. 方案 B：先用 LD 过滤大幅降低候选位点数，减少矩阵维度
  3. 方案 C：用增量贪心算法替代 MinimalMarkers 的全局优化（牺牲最优解换取内存可行性）
- **预防**：在大样本项目中，提前评估工具的内存复杂度。对于 O(n²) 内存的工具，样本数 1000 是实用上限。超过此阈值需设计迂回方案。
- **来源**：`ideas/chip-design-160k-panel.md`（实际工程踩坑）

---

### 重测序过估计杂合 — Hetero → Homo 差异达 49.99%

- **现象**：用芯片数据与 WGS 重测序数据做基因型一致性对比时，发现约 50% 的位点在 WGS 中被 calling 为杂合子，在芯片中被 calling 为纯合子。差异不是少数几个点，而是系统性的 ~49.99%。
- **原因**：重测序数据在某些区域因比对错误、GC 偏倚或低深度，系统性地过估计杂合度。尤其在重复区（如着丝粒周边），reads 的错配被 GATK HaplotypeCaller 解释为杂合信号。芯片探针因设计时避开了重复区，这些位点显示为纯合。
- **解决**：
  1. 在重复区过滤位点（RepeatMasker 掩蔽区内的位点）
  2. 对"Hetero→Homo"差异位点做 Sanger 验证确认哪一方正确
  3. 芯片数据通常更可靠（探针设计在 unique 区），但需确认
  4. WGS 数据加 GATK VQSR（variant quality score recalibration）过滤
- **预防**：基因型一致性对比分析中，按基因组区域分层（unique / repetitive / centromeric），分别统计一致性。不要仅用全基因组平均值。
- **来源**：`ideas/chip-design-160k-panel.md`（基因型一致性分析结果）

---

## 数据库

### `--curatedlib` 放错一条 TE 导致全基因组污染

- **现象**：EDTA `--curatedlib` 运行后，某个基因家族的所有成员（几十到几百个）全部被注释为 TE 并被 RepeatMasker 掩蔽。该基因家族的后续功能分析无法进行。
- **原因**：金标准库中混入了一条基因编码序列（如某个 LTR 逆转座子的 gag 蛋白在 BLAST 中偶然匹配到一个真实基因的保守结构域）。EDTA `--curatedlib` 对金标准库中的序列 100% 信任——不会做任何二次验证。一条误放序列会让 EDTA 在全基因组范围内对该家族的每个拷贝做注释。
- **解决**：
  1. 从金标准库中移除该污染序列
  2. 重新运行 `EDTA --curatedlib`（去掉该序列）
  3. 对被掩蔽的基因区域做 `blastn` 到 nt/nr 库，验证它们是真实基因
  4. 修复被污染的 GFF 注释
- **预防**：金标准库构建时的 CDS 过滤（`blastn -evalue 1e-10` 到植物 CDS 数据库）是**必需的**，不是可选的。同时随机抽取金标准库中 10 条序列 BLAST 到 nt 库做人工验证。切勿在未经验证的情况下将任何序列放入 `--curatedlib`。
- **来源**：`notes/bioinfo/annotation/te-annotation-pipeline.md`（`--curatedlib` 双刃剑警示）

---

### DTH vs DHH 混淆 — 只差一个字母但机制完全不同

- **现象**：在论文或报告中，将 Harbinger 转座子（DTH）讨论为"滚环复制"机制，或者将 Helitron（DHH）描述为"cut-and-paste"机制。审稿人指出分类错误。
- **原因**：Wicker 2007 三字母命名体系中，DTH 和 DHH 只差第三个字母（T vs H），但分属完全不同的转座子超家族：
  - **DHH** = Class II / **H**elitron — 滚环复制（rolling-circle replication），无 TIR
  - **DTH** = Class II / TIR / **H**arbinger — cut-and-paste 转座，有 TIR

  混淆通常发生在：(1) 自动化分类工具将 Helitron 误标为 DTH；(2) 人工阅读时因字母相似而看错；(3) TEsorter 输出中两个分类在相邻行。
- **解决**：
  - DHH: Helitron, rolling-circle, 5' TC + 3' CTRR 保守末端
  - DTH: TIR / Harbinger, cut-and-paste, 有 TIR 末端重复
- **预防**：在涉及 Helitron 和 Harbinger 的图表/文字中，反复核对三字母缩写。建立三字母速查表贴在实验室/显示器旁边（见 TE 注释 SOP 的附录）。
- **来源**：`notes/bioinfo/annotation/te-annotation-pipeline.md`（Wicker 2007 三字母命名速查 + DTH vs DHH 警示）

---

## 贡献指南

每次踩新坑后，按以下格式追加到对应分类下：

```markdown
### 坑名（简短描述）

- **现象：** <你看到了什么——具体的报错信息、异常结果、或奇怪的行为>
- **原因：** <为什么会发生——技术原理或操作失误>
- **解决：** <具体怎么做——可执行命令或操作步骤>
- **预防：** <下次怎么避免——工作流调整或检查点>
- **来源：** <哪个笔记/哪个脚本/哪个对话>
```

如果不确定属于哪个分类，先放在"其他"分类下，后续整理时再归类。

---

## 版本历史

| 版本 | 日期 | 变更内容 |
|------|------|---------|
| 1.0 | 2026-06-02 | 初始版本，含 7 大类 19 条踩坑记录 |

## 来源

本文件中的每条记录均标有来源，汇总自：
- `notes/` 下所有原子笔记中的"常见报错"、"常见问题"小节
- `ideas/` 下的半成品归档中的排错记录
- 用户实际脚本中的注释和排错日志（`debug&shell/`、`XTBG/` 等目录）
