# GATK 变异、GVCF 与结构变异工具（GATK，P01 包）

这组 30 个函数位于 `dnakit.variants.gatk_p01_variants_cli_functions`，调用 Broad Institute 的 [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 命令行（含其内置的 Picard 工具）。每个函数运行一个固定工具，或由模式参数选择的一族工具；不经过 shell，带超时和输出大小限制。`options` 只接受 `gatk <Tool> --help` 列出的长选项，已由参数设置的选项不能重复，GATK 会执行脚本的 `JAVASCRIPT_FILE` 不可用。输出写到显式的新路径，不会覆盖已有文件。

GATK 是外部程序：请安装（例如 `conda install -c bioconda gatk4` 或发布包）并准备 Java 17，然后把 `gatk` 启动器放进 PATH、设置 `DNAKIT_GATK` 或传入 `gatk_path`。找不到时抛出带安装提示的 `BackendUnavailableError`（代码 `GATK_UNAVAILABLE`）；工具运行失败抛出 `BackendExecutionError`（代码 `GATK_EXECUTION_FAILED`，附日志尾部）；工具返回 0 却未写出所请求的输出时抛出 `GATK_OUTPUT_MISSING`。DNAKit 不捆绑 GATK。Spark 工具以本地模式运行（可通过 `options` 设置 `spark-master`）。参考 FASTA 旁需要 `.fai` 和 `.dict`，BAM/VCF 需要各自索引，与 GATK 本身要求一致。结果为字典，每个函数的小节列出其额外字段。

## 1) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reblock_gvcf()`

- **作用:** 把单样本 GVCF 中的纯合参考块合并为更少更宽的块。
- **处理方法:** 运行 `gatk ReblockGVCF`。
- **API:** `gatk_reblock_gvcf(gvcf_path: 'str', reference_path: 'str', output_path: 'str', *, gvcf_gq_bands: 'Sequence[str] | None' = None, drop_low_quals: 'bool | None' = None, do_qual_score_approximation: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ReblockGVCF。
- **输入:** `gvcf_path` (必需): 单样本 GVCF（需已索引）。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的重新分块 GVCF。; `gvcf_gq_bands` (可选): GQ 分块边界。; `drop_low_quals` (可选): 丢弃低质量变异和非参考块。; `do_qual_score_approximation` (可选): 由等位基因特异注释估算 QUAL。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为重新分块的 GVCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reblock_gvcf

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_path = folder / "gvcf_path.dat"; gvcf_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_reblock_gvcf(
            gvcf_path=str(gvcf_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 2) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_calibrate_dragstr_model()`

- **作用:** 由 BAM 估计 DRAGstr 模型（感知 STR 的插入缺失错误模型）的参数。
- **处理方法:** 运行 `gatk CalibrateDragstrModel`。
- **API:** `gatk_calibrate_dragstr_model(input_path: 'str', reference_path: 'str', str_table_path: 'str', output_path: 'str', *, threads: 'int | None' = None, min_loci_count: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CalibrateDragstrModel。
- **输入:** `input_path` (必需): 输入 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `str_table_path` (必需): `gatk_compose_str_table` 生成的 STR 表。; `output_path` (必需): 新建的 DRAGstr 参数文件。; `threads` (可选): 线程数。; `min_loci_count` (可选): 每个周期/重复组合的最少位点数。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为参数文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_calibrate_dragstr_model

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    str_table_path = folder / "str_table_path.dat"; str_table_path.write_bytes(b"x")
    try:
        gatk_calibrate_dragstr_model(
            input_path=str(input_path),
            reference_path=str(reference_path),
            str_table_path=str(str_table_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 3) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_combine_gvcfs()`

- **作用:** 把一个或多个 HaplotypeCaller GVCF 文件合并为一个 GVCF。
- **处理方法:** 运行 `gatk CombineGVCFs`。
- **API:** `gatk_combine_gvcfs(gvcf_paths: 'Sequence[str]', reference_path: 'str', output_path: 'str', *, convert_to_base_pair_resolution: 'bool | None' = None, break_bands_at_multiples_of: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CombineGVCFs。
- **输入:** `gvcf_paths` (必需): 要合并的 GVCF 文件（需已索引）。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的合并 GVCF。; `convert_to_base_pair_resolution` (可选): 每个碱基写一条记录。; `break_bands_at_multiples_of` (可选): 在该值的倍数处断开参考带。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为合并后的 GVCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_combine_gvcfs

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_paths = folder / "gvcf_paths.dat"; gvcf_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_combine_gvcfs(
            gvcf_paths=[str(gvcf_paths)],
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 4) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_genomicsdb_import()`

- **作用:** 把单样本 GVCF 导入 GenomicsDB 工作区，供联合分型使用。
- **处理方法:** 运行 `gatk GenomicsDBImport`，创建工作区（`mode=create`）或向已有工作区追加样本（`mode=update`）。 `mode`: `create` = 新建工作区，`update` = 向已有工作区追加样本。
- **API:** `gatk_genomicsdb_import(*, mode: 'str' = 'create', workspace_path: 'str | None' = None, existing_workspace_path: 'str | None' = None, gvcf_paths: 'Sequence[str] | None' = None, sample_name_map: 'str | None' = None, intervals: 'Sequence[str] | None' = None, batch_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：GenomicsDBImport。
- **输入:** `workspace_path` (部分模式必需): 新建的 GenomicsDB 工作区目录（由 GATK 创建）。; `existing_workspace_path` (部分模式必需): 要追加样本的已有工作区。; `gvcf_paths` (可选): 要导入的 GVCF 文件（需已索引）；与 `sample_name_map` 二选一。; `sample_name_map` (可选): 两列文件 `样本<TAB>GVCF 路径`。; `intervals` (部分模式必需): 要导入的区间（`create` 必需）。; `batch_size` (可选): 每批读取的样本数（0 为全部）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；新建工作区时返回 `outputs.workspace_path` 与 `files`。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_genomicsdb_import

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_paths = folder / "gvcf_paths.dat"; gvcf_paths.write_bytes(b"x")
    try:
        gatk_genomicsdb_import(
            workspace_path=str(folder / "workspace_path"),
            gvcf_paths=[str(gvcf_paths)],
            intervals=["chr1"],
            mode='create',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 5) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_genotype_gvcfs()`

- **作用:** 对 GVCF 做联合分型（标准模型或快速的 GnarlyGenotyper 模型）。
- **处理方法:** 运行 `gatk GenotypeGVCFs`（`engine=standard`）或 beta 的 `gatk GnarlyGenotyper`（`engine=gnarly`，需 GenomicsDB 工作区或带等位基因特异注释的合并 GVCF）；区间通过 `intervals` 选项给出。 `engine`: `standard` = GenotypeGVCFs，`gnarly` = GnarlyGenotyper。
- **API:** `gatk_genotype_gvcfs(variant_source: 'str', reference_path: 'str', output_path: 'str', *, engine: 'str' = 'standard', dbsnp: 'str | None' = None, include_non_variant_sites: 'bool | None' = None, keep_all_sites: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：GenotypeGVCFs, GnarlyGenotyper。
- **输入:** `variant_source` (必需): 合并后的 GVCF，或 GenomicsDB 工作区 `gendb://<workspace>`。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的分型 VCF。; `dbsnp` (可选): 用于 ID 列的 dbSNP VCF。; `include_non_variant_sites` (可选): 同时输出非变异位点。; `keep_all_sites` (可选): 保留所有位点，包括没有备选等位基因的位点。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_genotype_gvcfs

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_source = folder / "variant_source.dat"; variant_source.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_genotype_gvcfs(
            variant_source=str(variant_source),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            engine='standard',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 6) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_haplotype_based_variant_recaller()`

- **作用:** 针对读段计算 VCF 中每个等位基因的似然，范围限定在一组单倍型内。
- **处理方法:** 运行 `gatk HaplotypeBasedVariantRecaller`（实验性）。
- **API:** `gatk_haplotype_based_variant_recaller(input_path: 'str', reference_path: 'str', alleles_vcf: 'str', haplotypes_bam: 'str', matrix_csv: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：HaplotypeBasedVariantRecaller。
- **输入:** `input_path` (必需): 读段 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `alleles_vcf` (必需): 要评估的等位基因 VCF。; `haplotypes_bam` (必需): 单倍型 BAM。; `matrix_csv` (必需): 保存似然矩阵的新建 CSV。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.matrix_csv` 为似然矩阵。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_haplotype_based_variant_recaller

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    alleles_vcf = folder / "alleles_vcf.dat"; alleles_vcf.write_bytes(b"x")
    haplotypes_bam = folder / "haplotypes_bam.dat"; haplotypes_bam.write_bytes(b"x")
    try:
        gatk_haplotype_based_variant_recaller(
            input_path=str(input_path),
            reference_path=str(reference_path),
            alleles_vcf=str(alleles_vcf),
            haplotypes_bam=str(haplotypes_bam),
            matrix_csv=str(folder / "matrix_csv.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 7) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_learn_read_orientation_model()`

- **作用:** 估计 FilterMutectCalls 使用的读段方向偏倚混合模型的伪迹先验的最大似然值。
- **处理方法:** 运行 `gatk LearnReadOrientationModel`。
- **API:** `gatk_learn_read_orientation_model(f1r2_tar_gz_files: 'Sequence[str]', output_path: 'str', *, num_em_iterations: 'int | None' = None, convergence_threshold: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：LearnReadOrientationModel。
- **输入:** `f1r2_tar_gz_files` (必需): Mutect2 写出的 F1R2 计数压缩包（`f1r2-tar-gz`）。; `output_path` (必需): 新建的方向偏倚模型压缩包（.tar.gz）。; `num_em_iterations` (可选): EM 最大迭代次数。; `convergence_threshold` (可选): EM 收敛阈值。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为模型压缩包。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_learn_read_orientation_model

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    f1r2_tar_gz_files = folder / "f1r2_tar_gz_files.dat"; f1r2_tar_gz_files.write_bytes(b"x")
    try:
        gatk_learn_read_orientation_model(
            f1r2_tar_gz_files=[str(f1r2_tar_gz_files)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 8) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_mutect2()`

- **作用:** 通过单倍型局部重组装检出体细胞 SNV 和 indel（Mutect2）。
- **处理方法:** 运行 `gatk Mutect2`；区域通过 `intervals` 选项给出。
- **API:** `gatk_mutect2(input_paths: 'Sequence[str]', reference_path: 'str', output_path: 'str', *, tumor_sample: 'str | None' = None, normal_samples: 'Sequence[str] | None' = None, germline_resource: 'str | None' = None, panel_of_normals: 'str | None' = None, f1r2_tar_gz: 'str | None' = None, bam_output: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：Mutect2。
- **输入:** `input_paths` (必需): 肿瘤（及正常）BAM/CRAM 文件。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的体细胞 VCF。; `tumor_sample` (可选): 肿瘤样本名（读组 SM）。; `normal_samples` (可选): 正常样本名。; `germline_resource` (可选): 人群等位基因频率 VCF（如 gnomAD）。; `panel_of_normals` (可选): 正常样本面板 VCF。; `f1r2_tar_gz` (可选): 供 `gatk_learn_read_orientation_model` 使用的 F1R2 计数新建压缩包。; `bam_output` (可选): 局部重组装读段的新建 BAM。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_mutect2

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_paths = folder / "input_paths.dat"; input_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_mutect2(
            input_paths=[str(input_paths)],
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 9) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reads_pipeline_spark()`

- **作用:** 在一个 Spark 流程中对读段运行 BWA（可选）、MarkDuplicates、BQSR 和 HaplotypeCaller，生成 VCF/GVCF。
- **处理方法:** 运行 `gatk ReadsPipelineSpark`（beta，Spark 本地模式）。
- **API:** `gatk_reads_pipeline_spark(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, output_bam: 'str | None' = None, align: 'bool | None' = None, emit_ref_confidence: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ReadsPipelineSpark。
- **输入:** `input_path` (必需): 未比对或已比对的 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `known_sites` (必需): BQSR 使用的已知变异 VCF 文件（需已索引）。; `output_path` (必需): 新建的 VCF/GVCF。; `output_bam` (可选): 处理后读段的新建 BAM。; `align` (可选): 先运行 BWA-MEM。; `emit_ref_confidence` (可选): 变异检出的参考置信模式。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为变异文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reads_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_reads_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 10) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_collect_sv_evidence()`

- **作用:** 从 BAM 中收集 GATK-SV 流程使用的双端、剪切读段和深度证据文件。
- **处理方法:** 运行 `gatk CollectSVEvidence`（beta）；写出所请求的证据文件，至少需请求一个。
- **API:** `gatk_collect_sv_evidence(input_path: 'str', sample_name: 'str', *, pe_file: 'str | None' = None, sr_file: 'str | None' = None, sd_file: 'str | None' = None, site_depth_locs_vcf: 'str | None' = None, depth_evidence_file: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CollectSVEvidence。
- **输入:** `input_path` (必需): 按坐标排序的 BAM/CRAM。; `sample_name` (必需): 写入证据文件的样本名。; `pe_file` (可选): 新建的双端（不一致配对）证据文件。; `sr_file` (可选): 新建的剪切读段证据文件。; `sd_file` (可选): 新建的位点深度证据文件。; `site_depth_locs_vcf` (可选): 位点深度计数所用位点的 VCF。; `depth_evidence_file` (可选): 新建的读深度证据文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs` 列出证据文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_collect_sv_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_collect_sv_evidence(
            input_path=str(input_path),
            sample_name="S1",
            pe_file=str(folder / "pe_file.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 11) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_condense_depth_evidence()`

- **作用:** 合并 GATK-SV 深度文件中相邻的深度证据记录。
- **处理方法:** 运行 `gatk CondenseDepthEvidence`（实验性）。
- **API:** `gatk_condense_depth_evidence(depth_evidence: 'str', output_path: 'str', *, max_interval_size: 'int | None' = None, min_interval_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CondenseDepthEvidence。
- **输入:** `depth_evidence` (必需): 深度证据文件（`.rd.txt`）。; `output_path` (必需): 新建的合并后深度证据文件。; `max_interval_size` (可选): 合并区间的最大长度。; `min_interval_size` (可选): 区间最小长度。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为合并后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_condense_depth_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    depth_evidence = folder / "depth_evidence.dat"; depth_evidence.write_bytes(b"x")
    try:
        gatk_condense_depth_evidence(
            depth_evidence=str(depth_evidence),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 12) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_convert_counts_to_depth_file()`

- **作用:** 把 GATK 读计数文件（`*.counts.tsv`）转换为 GATK-SV 深度文件（`.rd.txt`）。
- **处理方法:** 运行 `gatk ConvertCountsToDepthFile`（实验性）。
- **API:** `gatk_convert_counts_to_depth_file(counts_path: 'str', output_path: 'str', *, sample_name: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ConvertCountsToDepthFile。
- **输入:** `counts_path` (必需): 读计数文件（`*.counts.tsv`）。; `output_path` (必需): 新建的深度证据文件（`.rd.txt`）。; `sample_name` (可选): 写入文件的样本名。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为深度文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_convert_counts_to_depth_file

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    counts_path = folder / "counts_path.dat"; counts_path.write_bytes(b"x")
    try:
        gatk_convert_counts_to_depth_file(
            counts_path=str(counts_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 13) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_cluster()`

- **作用:** 按类型、大小和断点对结构变异聚类（可按分层和轨道重叠分组）。
- **处理方法:** 运行 `gatk SVCluster`（`mode=cluster`）或 `GroupedSVCluster`（`mode=grouped`，同时需要聚类与分层配置）；两者均为 beta 工具。 `mode`: `cluster` = SVCluster，`grouped` = GroupedSVCluster。
- **API:** `gatk_sv_cluster(variant_paths: 'Sequence[str]', reference_path: 'str', ploidy_table: 'str', output_path: 'str', *, mode: 'str' = 'cluster', clustering_config: 'str | None' = None, stratify_config: 'str | None' = None, algorithm: 'str | None' = None, variant_prefix: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：GroupedSVCluster, SVCluster。
- **输入:** `variant_paths` (必需): 结构变异 VCF 文件。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `ploidy_table` (必需): 制表符分隔的样本/contig 倍性表。; `output_path` (必需): 新建的聚类 VCF。; `clustering_config` (部分模式必需): 聚类配置表（按变异分组）。; `stratify_config` (部分模式必需): 分层配置表。; `algorithm` (可选): 聚类算法。; `variant_prefix` (可选): 新变异 ID 的前缀。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为聚类后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_cluster

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_paths = folder / "variant_paths.dat"; variant_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    ploidy_table = folder / "ploidy_table.dat"; ploidy_table.write_bytes(b"x")
    try:
        gatk_sv_cluster(
            variant_paths=[str(variant_paths)],
            reference_path=str(reference_path),
            ploidy_table=str(ploidy_table),
            output_path=str(folder / "output_path.out"),
            mode='cluster',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 14) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_joint_germline_cnv_segmentation()`

- **作用:** 把多个样本分段的种系拷贝数（gCNV）VCF 合并为一个队列 VCF。
- **处理方法:** 运行 `gatk JointGermlineCNVSegmentation`（beta）。
- **API:** `gatk_joint_germline_cnv_segmentation(variant_paths: 'Sequence[str]', pedigree: 'str', reference_path: 'str', output_path: 'str', *, model_call_intervals: 'str | None' = None, min_sample_set_fraction_overlap: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：JointGermlineCNVSegmentation。
- **输入:** `variant_paths` (必需): 各样本分段后的 gCNV VCF。; `pedigree` (必需): 家系（PED）文件。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的合并 VCF。; `model_call_intervals` (可选): gCNV 模型调用的区间列表。; `min_sample_set_fraction_overlap` (可选): 合并片段所需的最小重叠比例。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为合并后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_joint_germline_cnv_segmentation

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_paths = folder / "variant_paths.dat"; variant_paths.write_bytes(b"x")
    pedigree = folder / "pedigree.dat"; pedigree.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_joint_germline_cnv_segmentation(
            variant_paths=[str(variant_paths)],
            pedigree=str(pedigree),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 15) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_print_read_counts()`

- **作用:** 从读计数文件打印 CNV 判定用的各样本计数文件。
- **处理方法:** 运行 `gatk PrintReadCounts`（实验性）。输出文件名由 GATK 根据 `output_prefix` 决定。
- **API:** `gatk_print_read_counts(counts_path: 'str', *, output_prefix: 'str | None' = None, output_file_list: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PrintReadCounts。
- **输入:** `counts_path` (必需): 输入读计数文件。; `output_prefix` (可选): 输出路径前缀 `<前缀><样本>.counts.tsv`；所含目录必须已存在。; `output_file_list` (可选): 新建的两列列表 `样本<TAB>输出路径`。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；请求时含 `outputs.output_file_list`。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_print_read_counts

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    counts_path = folder / "counts_path.dat"; counts_path.write_bytes(b"x")
    try:
        gatk_print_read_counts(
            counts_path=str(counts_path),
            output_prefix="x",
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 16) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_print_sv_evidence()`

- **作用:** 把 GATK-SV 证据记录合并或按样本取子集，写入一个输出文件。
- **处理方法:** 运行 `gatk PrintSVEvidence`（实验性）。
- **API:** `gatk_print_sv_evidence(evidence_file: 'str', output_path: 'str', *, sample_names: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PrintSVEvidence。
- **输入:** `evidence_file` (必需): SV 证据文件（`.pe.txt`、`.sr.txt`、`.rd.txt`、`.sd.txt`、`.baf.txt`）。; `output_path` (必需): 新建的证据文件。; `sample_names` (可选): 要保留的样本。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为证据文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_print_sv_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evidence_file = folder / "evidence_file.dat"; evidence_file.write_bytes(b"x")
    try:
        gatk_print_sv_evidence(
            evidence_file=str(evidence_file),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 17) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_annotate()`

- **作用:** 给 SV VCF 添加基因重叠和变异后果预测注释。
- **处理方法:** 运行 `gatk SVAnnotate`。
- **API:** `gatk_sv_annotate(variant_path: 'str', output_path: 'str', *, protein_coding_gtf: 'str | None' = None, non_coding_bed: 'str | None' = None, promoter_window_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SVAnnotate。
- **输入:** `variant_path` (必需): GATK-SV 流程产生的结构变异 VCF。; `output_path` (必需): 新建的注释 VCF。; `protein_coding_gtf` (可选): 蛋白编码基因 GTF。; `non_coding_bed` (可选): 非编码元件 BED。; `promoter_window_length` (可选): 基因上游启动子窗口（bp）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为注释后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_annotate

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_sv_annotate(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 18) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_concordance()`

- **作用:** 给结构变异添加与真值集的基因型一致性注释。
- **处理方法:** 运行 `gatk SVConcordance`（beta）。
- **API:** `gatk_sv_concordance(evaluation: 'str', truth: 'str', output_path: 'str', *, clustering_config: 'str | None' = None, stratify_config: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SVConcordance。
- **输入:** `evaluation` (必需): 待评估 SV VCF。; `truth` (必需): 真值 SV VCF。; `output_path` (必需): 新建的、含一致性注释的 VCF。; `clustering_config` (可选): 聚类配置表。; `stratify_config` (可选): 分层配置表。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为注释后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_sv_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 19) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_stratify()`

- **作用:** 按 SV 类型、大小和参考轨道重叠对结构变异分层注释。
- **处理方法:** 运行 `gatk SVStratify`（beta）。
- **API:** `gatk_sv_stratify(variant_path: 'str', stratify_config: 'str', output_path: 'str', *, track_name: 'Sequence[str] | None' = None, track_intervals: 'Sequence[str] | None' = None, split_output: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SVStratify。
- **输入:** `variant_path` (必需): SV VCF。; `stratify_config` (必需): 分层配置表。; `output_path` (必需): 新建的分层 VCF。; `track_name` (可选): 参考轨道名。; `track_intervals` (可选): 参考轨道的区间文件。; `split_output` (可选): 每个分层写一个文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为分层后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_stratify

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    stratify_config = folder / "stratify_config.dat"; stratify_config.write_bytes(b"x")
    try:
        gatk_sv_stratify(
            variant_path=str(variant_path),
            stratify_config=str(stratify_config),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 20) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_site_depth_to_baf()`

- **作用:** 把位点深度证据转换为 GATK-SV 使用的 B 等位基因频率（BAF）证据。
- **处理方法:** 运行 `gatk SiteDepthtoBAF`（实验性）。
- **API:** `gatk_site_depth_to_baf(site_depth_paths: 'Sequence[str]', baf_sites_vcf: 'str', baf_evidence_output: 'str', *, sample_names: 'Sequence[str] | None' = None, min_het_probability: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SiteDepthtoBAF。
- **输入:** `site_depth_paths` (必需): 位点深度证据文件（`.sd.txt`）。; `baf_sites_vcf` (必需): 候选杂合位点 VCF。; `baf_evidence_output` (必需): 新建的 B 等位基因频率证据文件。; `sample_names` (可选): 输入对应的样本名。; `min_het_probability` (可选): 杂合位点的最小概率。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.baf_evidence_output` 为 BAF 文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_site_depth_to_baf

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    site_depth_paths = folder / "site_depth_paths.dat"; site_depth_paths.write_bytes(b"x")
    baf_sites_vcf = folder / "baf_sites_vcf.dat"; baf_sites_vcf.write_bytes(b"x")
    try:
        gatk_site_depth_to_baf(
            site_depth_paths=[str(site_depth_paths)],
            baf_sites_vcf=str(baf_sites_vcf),
            baf_evidence_output=str(folder / "baf_evidence_output.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 21) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_structural_variation_discovery_pipeline_spark()`

- **作用:** 对单个样本运行结构变异发现流程（证据收集、局部组装、断点检出）。
- **处理方法:** 运行 `gatk StructuralVariationDiscoveryPipelineSpark`（beta，Spark 本地模式）。
- **API:** `gatk_structural_variation_discovery_pipeline_spark(input_path: 'str', reference_path: 'str', aligner_index_image: 'str', kmers_to_ignore: 'str', contig_sam_file: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：StructuralVariationDiscoveryPipelineSpark。
- **输入:** `input_path` (必需): 单个样本的 BAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `aligner_index_image` (必需): 参考序列的 BWA-MEM 索引镜像。; `kmers_to_ignore` (必需): `gatk_find_bad_genomic_kmers_spark` 生成的 k-mer 列表。; `contig_sam_file` (必需): 新建的、重组装 contig 比对 SAM。; `output_path` (必需): 新建的 SV VCF。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 SV VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_structural_variation_discovery_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    aligner_index_image = folder / "aligner_index_image.dat"; aligner_index_image.write_bytes(b"x")
    kmers_to_ignore = folder / "kmers_to_ignore.dat"; kmers_to_ignore.write_bytes(b"x")
    try:
        gatk_structural_variation_discovery_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            aligner_index_image=str(aligner_index_image),
            kmers_to_ignore=str(kmers_to_ignore),
            contig_sam_file=str(folder / "contig_sam_file.out"),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 22) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_calculate_genotype_posteriors()`

- **作用:** 利用家系和/或已知群体基因型计算基因型后验概率。
- **处理方法:** 运行 `gatk CalculateGenotypePosteriors`。
- **API:** `gatk_calculate_genotype_posteriors(variant_path: 'str', output_path: 'str', *, pedigree: 'str | None' = None, supporting_callsets: 'Sequence[str] | None' = None, de_novo_prior: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CalculateGenotypePosteriors。
- **输入:** `variant_path` (必需): 含基因型似然的输入 VCF。; `output_path` (必需): 新建的、含基因型后验概率的 VCF。; `pedigree` (可选): 用于家系先验的家系（PED）文件。; `supporting_callsets` (可选): 提供等位基因计数先验的群体 VCF。; `de_novo_prior` (可选): 新生突变先验概率。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_calculate_genotype_posteriors

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_calculate_genotype_posteriors(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 23) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_concordance()`

- **作用:** 评估 VCF 与已验证真值 VCF 的一致性（真/假阳性、假阴性）。
- **处理方法:** 运行 `gatk Concordance`；汇总表解析到 `tables.summary_path`。区域通过 `intervals` 选项给出。
- **API:** `gatk_concordance(evaluation: 'str', truth: 'str', summary_path: 'str', *, true_positives_and_false_negatives: 'str | None' = None, true_positives_and_false_positives: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：Concordance。
- **输入:** `evaluation` (必需): 待评估 VCF。; `truth` (必需): 已验证的真值 VCF。; `summary_path` (必需): 新建的汇总表（灵敏度、精确度）。; `true_positives_and_false_negatives` (可选): 新建的 TP 与 FN 记录 VCF。; `true_positives_and_false_positives` (可选): 新建的 TP 与 FP 记录 VCF。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `tables.summary_path`（按变异类型的行）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            summary_path=str(folder / "summary_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 24) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_count_false_positives()`

- **作用:** 统计检出集中 PASS 变异的数量（样本无真变异时即假阳性数）。
- **处理方法:** 运行 `gatk CountFalsePositives`（beta）。
- **API:** `gatk_count_false_positives(variant_path: 'str', intervals: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CountFalsePositives。
- **输入:** `variant_path` (必需): 对真值样本的检出 VCF（所有检出均视为假阳性）。; `intervals` (必需): 限定计数范围的区间或区间文件。; `output_path` (必需): 保存计数的新建文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 保存计数。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_count_false_positives

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_count_false_positives(
            variant_path=str(variant_path),
            intervals=["chr1"],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 25) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_evaluate_info_field_concordance()`

- **作用:** 评估 VCF 的某个 INFO 字段与真值 VCF 对应字段的一致性。
- **处理方法:** 运行 `gatk EvaluateInfoFieldConcordance`（beta）。
- **API:** `gatk_evaluate_info_field_concordance(evaluation: 'str', truth: 'str', eval_info_key: 'str', truth_info_key: 'str', summary_path: 'str', *, epsilon: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：EvaluateInfoFieldConcordance。
- **输入:** `evaluation` (必需): 待评估 VCF。; `truth` (必需): 真值 VCF。; `eval_info_key` (必需): 待评估 VCF 中的 INFO 键。; `truth_info_key` (必需): 真值 VCF 中的 INFO 键。; `summary_path` (必需): 新建的汇总文件。; `epsilon` (可选): 比较数值时的容差。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.summary_path` 为汇总文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_evaluate_info_field_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_evaluate_info_field_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            eval_info_key="AF",
            truth_info_key="AF",
            summary_path=str(folder / "summary_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 26) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reference_block_concordance()`

- **作用:** 评估 GVCF 与真值 GVCF 的参考块一致性。
- **处理方法:** 运行 `gatk ReferenceBlockConcordance`。
- **API:** `gatk_reference_block_concordance(evaluation: 'str', truth: 'str', eval_block_histogram: 'str', truth_block_histogram: 'str', confidence_concordance_histogram: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ReferenceBlockConcordance。
- **输入:** `evaluation` (必需): 待评估 GVCF。; `truth` (必需): 真值 GVCF。; `eval_block_histogram` (必需): 新建的待评估参考块直方图。; `truth_block_histogram` (必需): 新建的真值参考块直方图。; `confidence_concordance_histogram` (必需): 新建的置信度一致性直方图。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs` 列出三个直方图。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reference_block_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_reference_block_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            eval_block_histogram=str(folder / "eval_block_histogram.out"),
            truth_block_histogram=str(folder / "truth_block_histogram.out"),
            confidence_concordance_histogram=str(folder / "confidence_concordance_histogram.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 27) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_filter_funcotations()`

- **作用:** 按临床意义相关的 Funcotator 注释（ClinVar、等位基因频率等）过滤变异。
- **处理方法:** 运行 `gatk FilterFuncotations`（实验性）。
- **API:** `gatk_filter_funcotations(variant_path: 'str', output_path: 'str', ref_version: 'str', allele_frequency_data_source: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FilterFuncotations。
- **输入:** `variant_path` (必需): Funcotator 注释后的 VCF。; `output_path` (必需): 新建的过滤后 VCF。; `ref_version` (必需): 注释所用的参考基因组版本。; `allele_frequency_data_source` (必需): 提供等位基因频率的群体数据库。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为过滤后的 VCF。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_filter_funcotations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_filter_funcotations(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            ref_version="hg38",
            allele_frequency_data_source="gnomad",
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 28) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_find_mendelian_violations()`

- **作用:** 在含三联体的 VCF 中找出所有类型的孟德尔违例。
- **处理方法:** 运行 Picard `gatk FindMendelianViolations`；指标文件解析到 `tables.output_path`。
- **API:** `gatk_find_mendelian_violations(input_path: 'str', pedigree_path: 'str', output_path: 'str', *, min_gq: 'int | None' = None, min_dp: 'int | None' = None, skip_chroms: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FindMendelianViolations。
- **输入:** `input_path` (必需): 含三联体成员的 VCF。; `pedigree_path` (必需): 三联体 PED 文件（无基因型列）。; `output_path` (必需): 新建的指标文件。; `min_gq` (可选): 最低基因型质量。; `min_dp` (可选): 最低深度。; `skip_chroms` (可选): 跳过的 contig。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `tables.output_path`（Picard 指标行）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_find_mendelian_violations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    pedigree_path = folder / "pedigree_path.dat"; pedigree_path.write_bytes(b"x")
    try:
        gatk_find_mendelian_violations(
            input_path=str(input_path),
            pedigree_path=str(pedigree_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 29) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_funcotate()`

- **作用:** 用 Funcotator 对变异或拷贝数片段做功能注释（基因、后果、ClinVar、gnomAD 等）。
- **处理方法:** 运行 `gatk Funcotator`（`input_kind=variants`）或 `gatk FuncotateSegments`（`input_kind=segments`，beta；上游输出格式尚未定型）。 `input_kind`: `variants` = Funcotator，`segments` = FuncotateSegments。
- **API:** `gatk_funcotate(input_path: 'str', reference_path: 'str', data_sources_path: 'str', ref_version: 'str', output_path: 'str', output_file_format: 'str', *, input_kind: 'str' = 'variants', transcript_selection_mode: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FuncotateSegments, Funcotator。
- **输入:** `input_path` (必需): 待注释的 VCF（`variants`）或片段文件（`segments`）。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `data_sources_path` (必需): Funcotator 数据源目录（见 FuncotatorDataSourceDownloader）。; `ref_version` (必需): 输入所用参考基因组版本，如 hg19 或 hg38。; `output_path` (必需): 新建的注释文件。; `output_file_format` (必需): 输出格式。; `transcript_selection_mode` (可选): 变异所属转录本的选择方式。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为注释后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_funcotate

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    data_sources_path = folder / "data_sources_path"; data_sources_path.mkdir()
    try:
        gatk_funcotate(
            input_path=str(input_path),
            reference_path=str(reference_path),
            data_sources_path=str(data_sources_path),
            ref_version="hg38",
            output_path=str(folder / "output_path.out"),
            output_file_format="VCF",
            input_kind='variants',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 30) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_validate_basic_somatic_short_mutations()`

- **作用:** 用同一批样本另行测序的肿瘤-正常 BAM 验证体细胞短变异（验证功效）。
- **处理方法:** 运行 `gatk ValidateBasicSomaticShortMutations`（实验性）。
- **API:** `gatk_validate_basic_somatic_short_mutations(variant_path: 'str', input_path: 'str', reference_path: 'str', discovery_sample_name: 'str', val_case_sample_name: 'str', val_control_sample_name: 'str', output_path: 'str', *, summary_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ValidateBasicSomaticShortMutations。
- **输入:** `variant_path` (必需): 待验证的发现检出 VCF。; `input_path` (必需): 验证用肿瘤与正常样本的 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `discovery_sample_name` (必需): 发现 VCF 中的样本名。; `val_case_sample_name` (必需): BAM 中的验证病例（肿瘤）样本名。; `val_control_sample_name` (必需): BAM 中的验证对照（正常）样本名。; `output_path` (必需): 新建的验证结果表。; `summary_path` (可选): 新建的汇总文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为验证结果表。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_validate_basic_somatic_short_mutations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_validate_basic_somatic_short_mutations(
            variant_path=str(variant_path),
            input_path=str(input_path),
            reference_path=str(reference_path),
            discovery_sample_name="T",
            val_case_sample_name="T2",
            val_control_sample_name="N2",
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```
