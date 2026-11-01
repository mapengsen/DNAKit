# GATK 注释评分模型与体细胞/线粒体过滤

`dnakit.variants.gatk_p02_cli_functions` 中的这些函数通过固定参数白名单（不经 shell）调用 GATK 命令行工具（测试版本 4.7.0.0）。需要安装 GATK 与 Java 17（`conda install -c bioconda gatk4` 或发布包），并把 `gatk` 放入 PATH 或传入 `gatk_path`。输出写入新路径（不覆盖已有文件），GATK 生成的索引文件列于 `index_files`。结果为字典，含 `status`（仅当退出码为 0 且输出存在时为 `completed`）、`version`、`command`、`outputs`、`log_tail` 及各工具特有的值。未设置的参数沿用 GATK 默认值。找不到 `gatk` 时抛出带安装提示的 `BackendUnavailableError`。

## 示例数据

示例在同一目录中按顺序运行，后面的示例会使用前面生成的文件。需要 `gatk` 在 PATH 中。

**示例文件：** `ref.fa`

```text
>chr1
TTTTTTTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATT
TGGTCTTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCG
AGGTTTCGATTTCCGCTATTTATCCTGATTTCCTGACCTGATTGTCGTAATTTTGTAGTCTAAAACAAGC
AAGATATTAAATCGTCAGCATCATCAGTGCAATAGTACATCATTATAACGTATTTAAAGTTTTAACTGTT
GCAACTCATTTCTTTTACAAAATCGGATAACGTATTGTGTCTAAGTATCTGGTTATACAATCGTTGATGC
TAATTGCTTAAAAAAGGTTTAAATTGGACCACAGCCGGATATCGGCAAGCTGGCTTTCCCGTTCTATCCT
ACTTACTAGTATTGGTCCAGACTCTCTGTAAGCATCCGCACTTCTTTGATAACAATATATTTCTTTTTAA
CGGAAACGTAGGATTTAATCAATTGTTATTAGAGGACCCTCTACTTAAGTTATCATGACAAAAGCCCTGC
CTGTAATCTTTATATATATGTGAATACGCTCGTATATTAT
>chr2
GTGCTGGCTTACTTGCGCATTGAGGGACAACGGTAGTACCCCCGGGGATCGTGGCATTACGGGCCGTTAG
GCTTCCTTGCTGTCATATACCGTGCTGACCTTCCAATCAGGTACAGATCTTCTAGAGGCTTTCATCACGG
CGACAACGGAGCTGAGGGGTAGCGGGCATTAACCGCGCGTCCACGGTAAGTTACAGGTTATGCGAGGATT
ATCGCTGAGTAGATTCAGAGACCCCGAAGTGAACTGAGGTCAGATAGGACTGACGCCGTATCTTGAACAA
ATGCTGGCGTCCGATTGACA
```

**示例文件：** `ref.fa.fai`

```text
chr1	600	6	70	71
chr2	300	621	70	71
```

**示例文件：** `ref.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:300
```

**示例文件：** `tumor.vcf`

```text
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
##FORMAT=<ID=FAD,Number=R,Type=Integer,Description="Count of fragments supporting each allele.">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 --output m2.vcf --intervals chr1:1-30000 --input rg.bam --reference genome.fa --f1r2-median-mq 50 --f1r2-min-bq 20 --f1r2-max-depth 200 --flow-likelihood-parallel-threads 0 --flow-likelihood-optimized-comp false --trim-to-haplotype true --exact-matching false --flow-use-t0-tag false --flow-remove-non-single-base-pair-indels false --flow-remove-one-zero-probs false --flow-quantization-bins 121 --flow-fill-empty-bins-value 0.001 --flow-symmetric-indel-probs false --flow-report-insertion-or-deletion false --flow-disallow-probs-larger-than-call false --flow-lump-probs false --flow-retain-max-n-probs-base-format false --flow-probability-scaling-factor 10 --flow-order-cycle-length 4 --keep-boundary-flows false --genotype-pon-sites false --genotype-germline-sites false --genotype-germline-sites-fraction 1.0 --af-of-alleles-not-in-resource -1.0 --mitochondria-mode false --permutect-ref-downsample 10 --permutect-alt-downsample 20 --permutect-non-artifact-ratio 1 --permutect-dataset-mode ILLUMINA --tumor-lod-to-emit 3.0 --initial-tumor-lod 2.0 --pcr-snv-qual 40 --pcr-indel-qual 40 --base-qual-correction-factor 5 --max-population-af 0.01 --downsampling-stride 1 --callable-depth 10 --max-suspicious-reads-per-alignment-start 0 --normal-lod 2.2 --ignore-itr-artifacts false --gvcf-lod-band -2.5 --gvcf-lod-band -2.0 --gvcf-lod-band -1.5 --gvcf-lod-band -1.0 --gvcf-lod-band -0.5 --gvcf-lod-band 0.0 --gvcf-lod-band 0.5 --gvcf-lod-band 1.0 --minimum-allele-fraction 0.0 --independent-mates false --flow-mode NONE --disable-adaptive-pruning false --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --min-dangling-branch-length 4 --recover-all-dangling-branches false --max-num-haplotypes-in-population 128 --min-pruning 2 --adaptive-pruning-initial-error-rate 0.001 --pruning-lod-threshold 2.302585092994046 --pruning-seeding-lod-threshold 9.210340371976184 --max-unpruned-variants 100 --linked-de-bruijn-graph false --disable-artificial-haplotype-recovery false --enable-legacy-graph-cycle-detection false --debug-assembly false --debug-graph-transformations false --capture-assembly-failure-bam false --num-matching-bases-in-dangling-end-to-recover -1 --error-correction-log-odds -Infinity --error-correct-reads false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --dragstr-het-hom-ratio 2 --dont-use-dragstr-pair-hmm-scores false --pair-hmm-gap-continuation-penalty 10 --expected-mismatch-rate-for-read-disqualification 0.02 --pair-hmm-implementation FASTEST_AVAILABLE --dragen-pd-pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --disable-symmetric-hmm-normalizing false --disable-cap-base-qualities-to-map-quality false --enable-dynamic-read-disqualification-for-genotyping false --dynamic-read-disqualification-threshold 1.0 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --native-dragen-pd-pair-hmm-threads 8 --native-dragen-pd-pair-avx-level FASTEST_AVAILABLE --native-dragen-pd-pair-memory-limit 1024 --native-dragen-pd-pair-open-mp-settings FASTEST_AVAILABLE --flow-hmm-engine-min-indel-adjust 6 --flow-hmm-engine-flat-insertion-penatly 45 --flow-hmm-engine-flat-deletion-penatly 45 --pileup-detection false --use-pdhmm false --use-pdhmm-overlap-optimization false --make-determined-haps-from-pd-code false --print-pileupcalling-status false --fallback-gga-if-pdhmm-fails true --pileup-detection-enable-indel-pileup-calling false --pileup-detection-active-region-phred-threshold 0.0 --num-artificial-haplotypes-to-add-per-allele 5 --artifical-haplotype-filtering-kmer-size 10 --pileup-detection-snp-alt-threshold 0.1 --pileup-detection-indel-alt-threshold 0.1 --pileup-detection-absolute-alt-depth 0.0 --pileup-detection-snp-adjacent-to-assembled-indel-range 5 --pileup-detection-snp-basequality-filter 12 --pileup-detection-bad-read-tolerance 0.0 --pileup-detection-proper-pair-read-badness true --pileup-detection-edit-distance-read-badness-threshold 0.08 --pileup-detection-chimeric-read-badness true --pileup-detection-template-mean-badness-threshold 0.0 --pileup-detection-template-std-badness-threshold 0.0 --pileup-detection-filter-assembly-alt-bad-read-tolerance 0.0 --pileup-detection-edit-distance-read-badness-for-assembly-filtering-threshold 0.12 --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --override-fragment-softclip-check false --min-base-quality-score 10 --smith-waterman FASTEST_AVAILABLE --emit-ref-confidence NONE --max-mnp-distance 1 --force-call-filtered-alleles false --reference-model-deletion-quality 30 --soft-clip-low-quality-ends false --allele-informative-reads-overlap-margin 2 --smith-waterman-dangling-end-match-value 25 --smith-waterman-dangling-end-mismatch-penalty -50 --smith-waterman-dangling-end-gap-open-penalty -110 --smith-waterman-dangling-end-gap-extend-penalty -6 --smith-waterman-haplotype-to-reference-match-value 200 --smith-waterman-haplotype-to-reference-mismatch-penalty -150 --smith-waterman-haplotype-to-reference-gap-open-penalty -260 --smith-waterman-haplotype-to-reference-gap-extend-penalty -11 --smith-waterman-read-to-haplotype-match-value 10 --smith-waterman-read-to-haplotype-mismatch-penalty -15 --smith-waterman-read-to-haplotype-gap-open-penalty -30 --smith-waterman-read-to-haplotype-gap-extend-penalty -5 --flow-assembly-collapse-hmer-size 0 --flow-assembly-collapse-partial-mode false --flow-filter-alleles false --flow-filter-alleles-qual-threshold 30.0 --flow-filter-alleles-sor-threshold 3.0 --flow-filter-lone-alleles false --flow-filter-alleles-debug-graphs false --min-assembly-region-size 50 --max-assembly-region-size 300 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --force-active false --assembly-region-padding 100 --padding-around-indels 75 --padding-around-snps 20 --padding-around-strs 75 --max-extension-into-assembly-region-padding-legacy 25 --max-reads-per-alignment-start 50 --enable-legacy-assembly-region-trimming false --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --output-cram-version 3.1 --create-output-variant-index true --create-output-variant-md5 false --max-variants-per-shard 0 --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false --minimum-mapping-quality 20 --max-read-length 2147483647 --min-read-length 30 --disable-tool-default-annotations false --enable-all-annotations false",Version="4.7.0.0",Date="October 20, 2026 at 10:25:25 AM CST">
##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.">
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
##INFO=<ID=CONTQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to contamination">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Number of potential somatic events in the assembly region">
##INFO=<ID=ECNTH,Number=A,Type=Integer,Description="Number of somatic events in best supporting haplotype for each alt allele">
##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality by allele">
##INFO=<ID=MFRL,Number=R,Type=Integer,Description="median fragment length by allele">
##INFO=<ID=MMQ,Number=R,Type=Integer,Description="median mapping quality by allele">
##INFO=<ID=MPOS,Number=A,Type=Integer,Description="median distance from end of read">
##INFO=<ID=NALOD,Number=A,Type=Float,Description="Log 10 odds of artifact in normal with same allele fraction as tumor">
##INFO=<ID=NCount,Number=1,Type=Integer,Description="Count of N bases in the pileup">
##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal log 10 likelihood ratio of diploid het or hom alt genotypes">
##INFO=<ID=OCM,Number=1,Type=Integer,Description="Number of alt reads whose original alignment doesn't match the current contig.">
##INFO=<ID=PON,Number=0,Type=Flag,Description="site found in panel of normals">
##INFO=<ID=POPAF,Number=A,Type=Float,Description="negative log 10 population allele frequencies of alt alleles">
##INFO=<ID=ROQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to read orientation artifact">
##INFO=<ID=RPA,Number=R,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=SEQQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not sequencing errors">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description="Phred-scaled quality of strand bias artifact">
##INFO=<ID=STRQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors">
##INFO=<ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">
##MutectVersion=2.2
##contig=<ID=chr1,length=600>
##contig=<ID=chr2,length=300>
##filtering_status=Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives.
##source=Mutect2
##tumor_sample=tumor
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tumor
chr1	179	.	T	C	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=4;ECNTH=2;MBQ=37,38;MFRL=0,0;MMQ=60,60;MPOS=17;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	0|1:2,1:0.400:3:0,0:2,1:2,1:0|1:179_T_C:179:2,0,1,0
chr1	180	.	G	T	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=4;ECNTH=2;MBQ=33,33;MFRL=0,0;MMQ=60,60;MPOS=28;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	1|0:2,1:0.400:3:0,0:2,1:2,1:1|0:179_T_C:179:2,0,1,0
chr1	189	.	A	C	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=4;ECNTH=2;MBQ=35,35;MFRL=0,0;MMQ=60,60;MPOS=37;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	1|0:2,1:0.400:3:0,0:2,1:2,1:1|0:179_T_C:179:2,0,1,0
chr1	191	.	T	G	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=4;ECNTH=2;MBQ=33,37;MFRL=0,0;MMQ=60,60;MPOS=5;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	0|1:2,1:0.400:3:0,0:2,1:2,1:0|1:179_T_C:179:2,0,1,0
chr1	434	.	G	A	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=2;ECNTH=2;MBQ=36,36;MFRL=0,0;MMQ=60,60;MPOS=58;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	0|1:2,1:0.400:3:0,0:2,1:2,1:0|1:434_G_A:434:2,0,1,0
chr1	442	.	C	G	.	.	AS_SB_TABLE=2,0|1,0;DP=3;ECNT=2;ECNTH=2;MBQ=31,40;MFRL=0,0;MMQ=60,60;MPOS=66;POPAF=7.30;TLOD=3.42	GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB	0|1:2,1:0.400:3:0,0:2,1:2,1:0|1:434_G_A:434:2,0,1,0
```

**示例文件：** `tumor.vcf.stats`

```text
statistic	value
callable	635.0
```

## 1) `extract_variant_annotations()`

- **作用：** 把 VCF 的位点级注释、标签和元数据提取为 HDF5 文件。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk ExtractVariantAnnotations`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.extract_variant_annotations(vcf: 'str', annotations: 'list[str]', output_prefix: 'str', *, resources: 'list[str | dict[str, Any]] | None' = None, modes: 'list[str] | None' = None, maximum_number_of_unlabeled_variants: 'int | None' = None, reservoir_sampling_random_seed: 'int | None' = None, resource_matching_strategy: 'str | None' = None, omit_alleles_in_hdf5: 'bool' = False, do_not_gzip_vcf_output: 'bool' = False, do_not_trust_all_polymorphic: 'bool' = False, ignore_all_filters: 'bool' = False, ignore_filters: 'list[str] | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF、注释名称、带标签的标注资源与变异类型。
- **输出：** `<prefix>.annot.hdf5`、`<prefix>.vcf.gz` 及未标注文件，列于 `prefix_outputs`。
- **示例代码：**

```python
from pathlib import Path
import random
import subprocess

rng = random.Random(7)
header = (
    "##fileformat=VCFv4.2\n"
    + "".join(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">\n' for t in ("QD", "FS", "MQ", "SOR"))
    + "##contig=<ID=chr1,length=70000>\n"
)
call = [header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
truth = ["##fileformat=VCFv4.2\n##contig=<ID=chr1,length=70000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
for i in range(900):
    pos, good = 100 + 75 * i, i % 4 != 0
    qd = rng.gauss(18, 3) if good else rng.gauss(4, 2)
    fs = abs(rng.gauss(2, 2)) if good else abs(rng.gauss(30, 10))
    mq = rng.gauss(60, 1) if good else rng.gauss(40, 5)
    sor = abs(rng.gauss(1, 0.5)) if good else abs(rng.gauss(4, 1))
    call.append(f"chr1\t{pos}\t.\tA\tC\t50\t.\tQD={qd:.3f};FS={fs:.3f};MQ={mq:.3f};SOR={sor:.3f}\n")
    if good and i % 3 == 0:
        truth.append(f"chr1\t{pos}\t.\tA\tC\t.\t.\t.\n")
Path("scored.vcf").write_text("".join(call))
Path("truth.vcf").write_text("".join(truth))
for name in ("scored.vcf", "truth.vcf"):
    subprocess.run(["gatk", "IndexFeatureFile", "-I", name], check=True, capture_output=True)
from dnakit.variants import gatk_p02_cli_functions as g

labels = [{"path": "truth.vcf", "name": "truth", "tags": {"training": True, "calibration": True}}]
names = ["QD", "FS", "MQ", "SOR"]
r = g.extract_variant_annotations("scored.vcf", names, "ext", resources=labels, modes=["SNP"])
print(r["status"], sorted(Path(p).name for p in r["prefix_outputs"]))
```

- **示例输出：**

```text
completed ['ext.annot.hdf5', 'ext.vcf.gz', 'ext.vcf.gz.tbi']
```

## 2) `train_variant_annotations_model()`

- **作用：** 基于提取的注释训练变异评分模型（默认孤立森林）。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk TrainVariantAnnotationsModel`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.train_variant_annotations_model(annotations_hdf5: 'str', output_prefix: 'str', *, unlabeled_annotations_hdf5: 'str | None' = None, model_backend: 'str' = 'PYTHON_IFOREST', hyperparameters_json: 'str | None' = None, python_script: 'str | None' = None, modes: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** 带标签的 HDF5、输出前缀、模型后端与类型。需要 GATK Python 环境（scikit-learn、h5py、dill）。
- **输出：** 评分器与训练/校准分数文件，列于 `prefix_outputs`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

import os
import sys

os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
labels = [{"path": "truth.vcf", "name": "truth", "tags": {"training": True, "calibration": True}}]
names = ["QD", "FS", "MQ", "SOR"]
r = g.train_variant_annotations_model("ext.annot.hdf5", "trn", modes=["SNP"])
print(r["status"], sorted(Path(p).name for p in r["prefix_outputs"]))
```

- **示例输出：**

```text
completed ['trn.snp.calibrationScores.hdf5', 'trn.snp.scorer.pkl', 'trn.snp.trainingScores.hdf5']
```

## 3) `score_variant_annotations()`

- **作用：** 用已训练模型为变异打分，并可按校准灵敏度硬过滤。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk ScoreVariantAnnotations`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.score_variant_annotations(vcf: 'str', annotations: 'list[str]', model_prefix: 'str', output_prefix: 'str', *, resources: 'list[str | dict[str, Any]] | None' = None, modes: 'list[str] | None' = None, model_backend: 'str' = 'PYTHON_IFOREST', python_script: 'str | None' = None, snp_calibration_sensitivity_threshold: 'float | None' = None, indel_calibration_sensitivity_threshold: 'float | None' = None, score_key: 'str | None' = None, calibration_sensitivity_key: 'str | None' = None, low_score_filter_name: 'str | None' = None, snp_key: 'str | None' = None, double_format: 'str | None' = None, resource_matching_strategy: 'str | None' = None, omit_alleles_in_hdf5: 'bool' = False, do_not_gzip_vcf_output: 'bool' = False, do_not_trust_all_polymorphic: 'bool' = False, ignore_all_filters: 'bool' = False, ignore_filters: 'list[str] | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF、注释名称、模型前缀、校准资源与阈值。
- **输出：** `scored_vcf`（SCORE、CALIBRATION_SENSITIVITY、LOW_SCORE 过滤）、`summary` 与 HDF5 文件。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

import os
import sys

os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
labels = [{"path": "truth.vcf", "name": "truth", "tags": {"training": True, "calibration": True}}]
names = ["QD", "FS", "MQ", "SOR"]
r = g.score_variant_annotations(
    "scored.vcf", names, "trn", "sc", resources=labels, modes=["SNP"],
    snp_calibration_sensitivity_threshold=0.95,
)
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **示例输出：**

```text
completed 900 {'.': 648, 'LOW_SCORE': 252}
```

## 4) `filter_mutect_calls()`

- **作用：** 过滤 Mutect2 调出的体细胞 SNV 与 indel。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk FilterMutectCalls`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.filter_mutect_calls(vcf: 'str', reference: 'str', output_path: 'str', *, stats: 'str | None' = None, filtering_stats_path: 'str | None' = None, contamination_tables: 'list[str] | None' = None, contamination_estimate: 'float | None' = None, tumor_segmentation: 'list[str] | None' = None, orientation_bias_priors: 'list[str] | None' = None, threshold_strategy: 'str | None' = None, false_discovery_rate: 'float | None' = None, f_score_beta: 'float | None' = None, initial_threshold: 'float | None' = None, min_allele_fraction: 'float | None' = None, unique_alt_read_count: 'int | None' = None, min_median_base_quality: 'int | None' = None, min_median_mapping_quality: 'int | None' = None, min_median_read_position: 'int | None' = None, max_events_in_region: 'int | None' = None, max_alt_allele_count: 'int | None' = None, normal_p_value_threshold: 'float | None' = None, microbial_mode: 'bool' = False, mitochondria_mode: 'bool' = False, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** Mutect2 VCF 及其 .stats 文件、参考序列，可选污染表、方向偏倚先验与阈值。
- **输出：** 过滤后的 VCF、`summary`，以及可选的过滤统计文件。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.filter_mutect_calls("tumor.vcf", "ref.fa", "filtered.vcf", filtering_stats_path="filtered.stats")
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **示例输出：**

```text
completed 6 {'PASS': 2, 'clustered_events;haplotype': 4}
```

## 5) `create_somatic_panel_of_normals()`

- **作用：** 由正常样本的 GenomicsDB 工作区生成 Mutect2 的 panel of normals。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk CreateSomaticPanelOfNormals`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.create_somatic_panel_of_normals(genomicsdb_workspace: 'str', output_path: 'str', *, reference: 'str | None' = None, germline_resource: 'str | None' = None, min_sample_count: 'int | None' = None, max_germline_probability: 'float | None' = None, call_genotypes: 'bool' = False, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** GenomicsDB 工作区目录（GenomicsDBImport 产出）、最少样本数，可选种系资源。
- **输出：** panel VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

import subprocess

text = Path("tumor.vcf").read_text()
for sample in ("normal1", "normal2"):
    Path(f"{sample}.vcf").write_text(text.replace("\ttumor\n", f"\t{sample}\n"))
    subprocess.run(["gatk", "IndexFeatureFile", "-I", f"{sample}.vcf"], check=True, capture_output=True)
subprocess.run(
    ["gatk", "GenomicsDBImport", "-R", "ref.fa", "--genomicsdb-workspace-path", "pon_db",
     "-L", "chr1:1-600", "-V", "normal1.vcf", "-V", "normal2.vcf"],
    check=True, capture_output=True,
)
r = g.create_somatic_panel_of_normals("pon_db", "pon.vcf", reference="ref.fa", min_sample_count=2)
print(r["status"], r["summary"]["records"], r["summary"]["samples"])
```

- **示例输出：**

```text
completed 6 []
```

## 6) `filter_alignment_artifacts()`

- **作用：** 通过对每个变异周围的 reads 重比对来过滤比对假象。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk FilterAlignmentArtifacts`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.filter_alignment_artifacts(vcf: 'str', reads: 'str', reference: 'str', bwa_mem_index_image: 'str', output_path: 'str', *, bam_output_path: 'str | None' = None, dont_skip_filtered_variants: 'bool' = False, kmer_size: 'int | None' = None, drop_ratio: 'float | None' = None, min_seed_length: 'int | None' = None, seed_split_factor: 'float | None' = None, indel_start_tolerance: 'int | None' = None, max_reasonable_fragment_length: 'int | None' = None, min_aligner_score_difference_per_base: 'float | None' = None, min_mismatch_difference_per_base: 'float | None' = None, num_regular_contigs: 'int | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF、其来源 BAM、参考序列及其 BWA-MEM 索引镜像。
- **输出：** 过滤后的 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

import random
import subprocess

import pysam

seq = "".join(x.strip() for x in Path("ref.fa").read_text().split(">chr2")[0].splitlines()[1:])
variants = [
    x.split("\t") for x in Path("filtered.vcf").read_text().splitlines()
    if not x.startswith("#") and x.split("\t")[6] == "PASS"
]
rng = random.Random(1)
head = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 600}, {"SN": "chr2", "LN": 300}],
    "RG": [{"ID": "g1", "SM": "tumor", "PL": "ILLUMINA", "LB": "l", "PU": "u"}],
}
with pysam.AlignmentFile("unsorted.bam", "wb", header=head) as out:
    for i, start in enumerate(range(0, 500, 5)):
        bases = list(seq[start : start + 100])
        for v in variants:
            offset = int(v[1]) - 1 - start
            if 0 <= offset < 100 and rng.random() < 0.4:
                bases[offset] = v[4]
        read = pysam.AlignedSegment(out.header)
        read.query_name, read.flag, read.reference_id = f"r{i}", 0, 0
        read.reference_start, read.mapping_quality, read.cigarstring = start, 60, "100M"
        read.query_sequence = "".join(bases)
        read.query_qualities = pysam.qualitystring_to_array("I" * 100)
        read.set_tag("RG", "g1")
        out.write(read)
pysam.sort("-o", "tumor.bam", "unsorted.bam")
pysam.index("tumor.bam")
subprocess.run(["gatk", "BwaMemIndexImageCreator", "-I", "ref.fa", "-O", "ref.img"], check=True, capture_output=True)
r = g.filter_alignment_artifacts("filtered.vcf", "tumor.bam", "ref.fa", "ref.img", "realigned.vcf")
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **示例输出：**

```text
completed 2 {'PASS': 2}
```

## 7) `filter_mt_low_heteroplasmy()`

- **作用：** 当通过其他过滤的低异质性位点过多时，过滤全部低异质性位点。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk MTLowHeteroplasmyFilterTool`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.filter_mt_low_heteroplasmy(vcf: 'str', output_path: 'str', *, low_het_threshold: 'float | None' = None, max_allowed_low_hets: 'int | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** Mutect2 线粒体 VCF、低异质性阈值与允许数量。
- **输出：** 应用 `mt_many_low_hets` 过滤后的 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.filter_mt_low_heteroplasmy("filtered.vcf", "lowhet.vcf", low_het_threshold=0.5, max_allowed_low_hets=1)
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **示例输出：**

```text
completed 6 {'PASS': 2, 'clustered_events;haplotype': 4}
```

## 8) `filter_numt_variants()`

- **作用：** 依据等位基因深度与常染色体覆盖度标记疑似 NuMT 的等位基因。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk NuMTFilterTool`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.filter_numt_variants(vcf: 'str', output_path: 'str', *, autosomal_coverage: 'float | None' = None, max_numt_autosomal_copies: 'float | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** Mutect2 线粒体 VCF 与常染色体中位覆盖度。
- **输出：** 带 `possible_numt` 过滤的 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.filter_numt_variants("filtered.vcf", "numt.vcf", autosomal_coverage=100.0)
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **示例输出：**

```text
completed 6 {'clustered_events;haplotype;possible_numt': 4, 'possible_numt': 2}
```
