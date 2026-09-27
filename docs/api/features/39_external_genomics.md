# 外部基因组分析与 MCP

通过统一 Python API 检索已登记功能、调用外部分析软件和连接第三方 MCP。安装、完整参数与运行约定见[参数 FAQ](../../faq/api/39_external_genomics.md)，作者、论文和许可证见[来源 FAQ](../../faq/genomics_sources.md)。

希望通过 pip 一并安装计算依赖，可先使用[基因组教程中的 `genomics` 安装选项](38_genomics_extensions.md#pip-genomics)。它提供进程内的基因型 PCA、ANI、原核 CDS 预测等接口；本页 `run_analysis()` 继续使用所登记的外部可执行程序，二者的输入和算法范围分别说明。43 项功能的逐项迁移策略见[FAQ](../../faq/external_analysis_migration.md)。

## 1) 检索项目与可执行功能

- **作用：** 查找相关项目与 DNAKit 已登记的分析接口，读取所需输入和后端。
- **计算方法：** 对本地目录记录进行关键词匹配，并按领域筛选；`describe_analysis` 根据分析 ID 返回接口定义。检索过程不运行分析。
- **API：** `dnakit.integrations.search_projects(query[可选], domain[可选], verbose[可选])`；`dnakit.integrations.list_analyses(query[可选], domain[可选], verbose[可选])`；`dnakit.integrations.describe_analysis(analysis[必须], verbose[可选])`。
- **输入：** 项目关键词、领域名称或分析 ID，例如 `plink2`、`population`、`plink2-pca`。
- **输出：** 返回项目列表、可用分析列表和分析定义字典。
- **示例代码：**

```python
from dnakit.integrations import search_projects, list_analyses, describe_analysis

projects = search_projects("plink2", verbose=False)
analyses = list_analyses("plink2", verbose=False)
definition = describe_analysis("plink2-pca", verbose=False)
print([project["id"] for project in projects["projects"]])
print([analysis["id"] for analysis in analyses])
print(definition["backend"], sorted(definition["inputs"]))
```

- **示例结果：** 当前目录的离线查询结果。

```text
['plink2']
['plink2-pca', 'plink2-king', 'plink2-ld-prune', 'plink2-gwas']
plink2 ['genotypes']
```

## 2) 外部软件调用

- **作用：** 根据已登记的分析定义生成运行计划，并在准备好后端和输入后执行分析。
- **计算方法：** 校验输入文件、参数和输出目录，按登记模板构造命令参数。`plan_analysis` 返回计划；`run_analysis` 调用后端并检查输出文件，保存运行记录。
- **API：** `dnakit.integrations.plan_analysis(analysis[必须], inputs[必须], output_directory[必须], parameters[可选], threads[可选], executable_paths[可选], verbose[可选])`；`dnakit.integrations.run_analysis(analysis[必须], inputs[必须], output_directory[必须], parameters[可选], threads[可选], executable_paths[可选], timeout_seconds[可选], max_artifact_bytes[可选], verbose[可选])`。
- **输入：** 分析 ID、输入文件映射和尚不存在的输出目录。以下示例需已有 `data/cohort.bed`、`data/cohort.bim`、`data/cohort.fam`。
- **输出：** `plan_analysis` 返回执行计划；实际调用返回状态、步骤和产物路径。
- **示例代码：**

```python
from dnakit.integrations import plan_analysis, run_analysis

inputs = {"genotypes": "data/cohort"}
plan = plan_analysis(
    "plink2-pca", inputs, "results/cohort-pca",
    parameters={"components": 3}, threads=4, verbose=False,
)
print(plan["analysis"], plan["status"])
print(plan["parameters"])
print(plan["steps"][0]["expected_outputs"])

# 已安装 PLINK 2，并确认输入与计划后执行：
# result = run_analysis(
#     "plink2-pca", inputs, "results/cohort-pca",
#     parameters={"components": 3}, threads=4,
# )
```

- **示例结果：** 输入文件就绪时的计划输出；此步骤不执行 PLINK，也不写结果文件。

```text
plink2-pca planned
{'components': 3}
['pca.eigenvec', 'pca.eigenval']
```

<details markdown="1">
<summary>已登记分析与后端</summary>

<!-- generated:analyses:start -->
| Analysis ID / 分析 ID | Capability / 功能 | Backend / 后端 | Validation / 验证 |
| --- | --- | --- | --- |
| `plink2-pca` | 群体结构主成分 | plink2 | passed_synthetic_fixture |
| `plink2-king` | 样本亲缘系数 | plink2 | passed_synthetic_fixture |
| `plink2-ld-prune` | 连锁不平衡剪枝 | plink2 | passed_synthetic_fixture |
| `plink2-gwas` | 带协变量的全基因组关联分析 | plink2 | passed_synthetic_fixture |
| `gcta-grm` | 基因组关系矩阵文件 | gcta64 | not_run |
| `gcta-reml` | 遗传力估计 | gcta64 | not_run |
| `minigraph-pangenome` | 参考基因组与组装图构建 | minigraph | not_run |
| `hifiasm-assembly` | HiFi 单样本组装 | hifiasm | not_run |
| `flye-assembly` | Nanopore 长读长组装 | flye | not_run |
| `sniffles-sv` | 长读长结构变异检测 | sniffles | not_run |
| `modkit-pileup` | 修饰碱基覆盖统计 | modkit | not_run |
| `mosdepth-coverage` | 测序深度统计 | mosdepth | not_run |
| `fastani` | 微生物平均核苷酸一致性 | fastANI | not_run |
| `kraken2-classify` | 微生物读段分类 | kraken2 | not_run |
| `amrfinder-nucleotide` | 组装序列耐药基因注释 | amrfinder | not_run |
| `prokka-annotation` | 原核基因组注释 | prokka | not_run |
| `barrnap-rrna` | 基因组核糖体 RNA 位点识别 | barrnap | not_run |
| `fimo-motifs` | DNA 模体位点搜索 | fimo | not_run |
| `cas-offinder` | CRISPR 脱靶候选枚举 | cas-offinder | not_run |
| `bakta-annotation` | 细菌、MAG与质粒标准化注释 | bakta | not_run |
| `busco-proteins` | 蛋白集完整性评估 | busco | not_run |
| `busco-genome` | 基因组完整性评估 | busco | not_run |
| `augustus-genes` | 真核基因结构预测 | augustus | not_run |
| `diamond-blastp` | 蛋白相似性搜索 | diamond | passed_synthetic_smoke |
| `diamond-blastx` | 翻译核酸相似性搜索 | diamond | passed_synthetic_smoke |
| `eggnog-proteins` | 直系同源支持的蛋白功能注释 | emapper.py | not_run |
| `miniprot-align` | 含剪接和移码的蛋白到基因组比对 | miniprot | passed_synthetic_smoke |
| `liftoff-annotation` | 组装之间的基因注释迁移 | liftoff | not_run |
| `checkm2-predict` | MAG完整性与污染度预测 | checkm2 | not_run |
| `quast-assembly` | 参考支持的组装质量评估 | quast.py | not_run |
| `vsearch-dereplicate` | 核酸序列去重复 | vsearch | passed_synthetic_smoke |
| `vsearch-cluster` | 核酸相似性聚类 | vsearch | passed_synthetic_smoke |
| `mmseqs-protein-search` | 蛋白序列快速搜索 | mmseqs | passed_synthetic_smoke |
| `mmseqs-protein-cluster` | 蛋白序列聚类 | mmseqs | passed_synthetic_smoke |
| `abricate-screen` | 组装序列耐药与毒力标记筛查 | abricate | not_run |
| `mlst-typing` | PubMLST方案多位点序列分型 | mlst | not_run |
| `genomad-mobile-elements` | 病毒和质粒等移动元件识别 | genomad | not_run |
| `repeatmasker-library` | 指定重复序列库的核酸屏蔽 | RepeatMasker | not_run |
| `mageck-test` | CRISPR筛选基因效应排序 | mageck | not_run |
| `mageck-count` | CRISPR筛选sgRNA读数统计 | mageck | not_run |
| `crispresso2-amplicon` | 单扩增子编辑结果定量 | CRISPResso | not_run |
| `drugz-unpaired` | 非配对CRISPR药物交互筛选 | drugz.py | passed_synthetic_smoke |
| `gemma-lmm` | 亲缘校正的混合模型关联检验 | gemma | not_run |
<!-- generated:analyses:end -->

</details>

## 3) 第三方 MCP 连接

- **作用：** 查看已登记服务，并通过实时工具清单选择和调用第三方能力。
- **计算方法：** `list_mcp_servers` 读取本地目录；`list_mcp_tools` 连接指定服务并获取工具 schema；`call_mcp_tool` 校验参数后调用一个工具。
- **API：** `dnakit.integrations.list_mcp_servers(verbose[可选])`；`dnakit.integrations.list_mcp_tools(server_id[必须], transport[可选], executable_path[可选], workdir[可选])`；`dnakit.integrations.call_mcp_tool(server_id[必须], tool_name[必须], arguments[必须], transport[可选], executable_path[可选], workdir[可选])`。连接超时、环境变量和响应限制见[参数 FAQ](../../faq/api/39_external_genomics.md#mcp-connection)。
- **输入：** 目录中的服务 ID；实际调用还需已配置的服务，以及实时 schema 中的工具名称和参数。
- **输出：** 返回 MCP 服务和工具列表；每项包含服务 ID、命令或工具定义。
- **示例代码：**

```python
from dnakit.integrations import list_mcp_servers, list_mcp_tools, call_mcp_tool

server = next(
    item for item in list_mcp_servers(verbose=False) if item["id"] == "biomcp"
)
print(server["id"], server["command"])
print(server["arguments"])

# 配置好 BioMCP 后，查询实时工具；这一步会启动服务：
# tools = list_mcp_tools("biomcp")
# 按 tools 返回的工具名称和 inputSchema 准备参数，再调用：
# result = call_mcp_tool("biomcp", tool_name, arguments)
```

- **示例结果：** 离线服务目录结果；未启动服务或调用远端工具。

```text
biomcp biomcp
['serve']
```

<details markdown="1">
<summary>已登记 MCP 服务与上游配置</summary>

<!-- generated:mcp:start -->
| Server ID / 服务 ID | Upstream / 上游 | Configuration / 配置 |
| --- | --- | --- |
| `gigwa-mcp` | [Gigwa MCP](https://github.com/gkanogiannis/Gigwa-MCP) | Install gigwa-mcp separately into a Python >=3.10 environment; configure GIGWA_URL and optional paired GIGWA_USER/GIGWA_PASS. |
| `chorus` | [Chorus](https://github.com/pinellolab/chorus) | Install the pinellolab/chorus source project and provision its oracle environments, model weights and references separately. The PyPI package named chorus is unrelated. Set executable_path to the installed chorus-mcp. |
| `mavedb-link` | [MaveDB-Link](https://github.com/berntpopp/mavedb-link) | For stdio, install the source project's dependencies in Python >=3.12 and set workdir to that checkout and executable_path to its Python. The author-documented HTTP service is available with transport='http'. |
| `biomcp` | [BioMCP](https://github.com/genomoncology/biomcp) | Install genomoncology/biomcp separately; its PyPI distribution is biomcp-cli, not the unrelated biomcp distribution. The installed executable is biomcp. |
| `juicebox` | [Juicebox MCP](https://github.com/aidenlab/juicebox-mcp) | Prepare the official source checkout and npm dependencies separately, then set workdir to the checkout and executable_path to Node if necessary. |
| `genechat` | [GeneChat MCP](https://github.com/natecostello/genechat-mcp) | Install and initialize GeneChat and its local genome/annotation databases separately. GENECHAT_CONFIG optionally selects an existing TOML configuration. |
| `genepattern` | [GenePattern MCP](https://github.com/genepattern/genepattern-mcp) | Install the official source checkout requirements separately with Python >=3.11. Set workdir to the checkout, executable_path to that Python, and GENEPATTERN_KEY to an existing API token. |
| `brapi` | [BrAPI MCP](https://github.com/cyanheads/brapi-mcp-server) | Separately prepare the pinned @cyanheads/brapi-mcp-server release/source with its dist build and dependencies; set workdir to the package root and executable_path to Node >=24. The bridge never runs npm, npx, bunx or an installer. |
| `ensembl` | [Ensembl MCP (cyanheads)](https://github.com/cyanheads/ensembl-mcp-server) | Separately prepare the pinned @cyanheads/ensembl-mcp-server build and dependencies; set workdir to the package root and executable_path to Node >=24. |
| `biothings` | [BioThings MCP](https://github.com/longevity-genie/biothings-mcp) | Install the pinned longevity-genie/biothings-mcp project separately in Python >=3.10; use its installed biothings-mcp entry point. |
| `gget` | [gget MCP](https://github.com/longevity-genie/gget-mcp) | Install gget-mcp and its separately specified gget dependencies in Python >=3.10; the verified command is gget-mcp stdio. |
| `ucsc-genome` | [UCSC Genome MCP](https://github.com/hlydecker/ucsc-genome-mcp) | Prepare the pinned source checkout and its requirements in Python >=3.10, set workdir to the checkout and executable_path to that Python. Use the existing hyphenated ucsc-genome-mcp.py script; upstream README/packaging also mention nonexistent underscore filenames. |
| `opentargets` | [Open Targets MCP](https://github.com/nickzren/opentargets-mcp) | Install nickzren/opentargets-mcp separately in Python >=3.10; set executable_path to its interpreter. Uses the verified Python module entry point with explicit stdio. |
| `plant-genomics` | [Plant Genomics MCP](https://github.com/musharna/plant-genomics-mcp) | Install the pinned plant-genomics-mcp Python distribution and dependencies separately. The plant-genomics-mcp entry point selects stdio; no HTTP token is needed for this local transport. |
| `biocontext-kb` | [BioContextAI Knowledgebase MCP](https://github.com/biocontext-ai/knowledgebase-mcp) | Install biocontext_kb separately in Python >=3.12. DNAKit fixes MCP_ENVIRONMENT=DEVELOPMENT for the verified stdio entry point. |
| `grebi` | [GrEBI](https://github.com/EBISPOT/GrEBI) | Select transport='http' for the author-documented development service. Local deployment is a separate Java/API plus graph-database and HPC data-pipeline task; no local startup recipe is claimed here. |
| `ebi-tools` | [EBI Tools MCP (QuentinCody)](https://github.com/QuentinCody/ebi-tools-mcp-server) | Select transport='http' for the repository author's published Cloudflare Worker MCP service. No standalone Node stdio entry point is present; no mcp-remote/npx helper is required by DNAKit. |
| `gwas-mcp` | [GWAS MCP (muslus)](https://github.com/muslus/gwas-mcp) | Prepare the pinned muslus/gwas-mcp source and its dependencies in an isolated Python >=3.10 environment, then point executable_path at that interpreter. Avoid assuming that the separately published gwas-mcp distribution has identical code. |
| `ncbi-database` | [NCBI Database MCP (hpend2373)](https://github.com/hpend2373/NCBI-Database-MCP) | Prepare the pinned NCBI-Database-MCP checkout and its requirements, set workdir to the checkout and executable_path to its Python. Uses the README-supported standard MCP script rather than executing its shell installer/launcher. |
| `pubmed` | [PubMed MCP (cyanheads)](https://github.com/cyanheads/pubmed-mcp-server) | Prepare the pinned @cyanheads/pubmed-mcp-server build and dependencies separately, then set workdir to the package root and executable_path to Node >=24. |
| `encode-toolkit` | [ENCODE Toolkit](https://github.com/ammawla/encode-toolkit) | Install the Python encode-toolkit distribution separately and set executable_path to its interpreter. DNAKit invokes the verified Python module directly; it does not use the npm shim which can bootstrap uv/Python. |
| `tooluniverse` | [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| `kiro-genomics` | [AWS Kiro Life Sciences: genomics](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-genomics component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-genomics console script. No uvx or automatic install is used. |
| `kiro-model-organisms` | [AWS Kiro Life Sciences: model-organisms](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-model-organisms component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-model-organisms console script. No uvx or automatic install is used. |
| `kiro-epigenomics` | [AWS Kiro Life Sciences: epigenomics](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-epigenomics component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-epigenomics console script. No uvx or automatic install is used. |
| `kiro-agriculture` | [AWS Kiro Life Sciences: agriculture](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-agriculture component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-agriculture console script. No uvx or automatic install is used. |
| `kiro-ecology` | [AWS Kiro Life Sciences: ecology](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-ecology component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-ecology console script. No uvx or automatic install is used. |
| `kiro-clinical` | [AWS Kiro Life Sciences: clinical](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-clinical component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-clinical console script. No uvx or automatic install is used. |
| `kiro-pipelines` | [AWS Kiro Life Sciences: pipelines](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-pipelines component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-pipelines console script. No uvx or automatic install is used. |
| `pubmed-anthropic` | [PubMed (Anthropic Life Sciences hosted connector)](https://github.com/anthropics/life-sciences) | Select transport=http. The provider-published endpoint supports anonymous MCP initialization and tools/list as verified on 2026-09-20. No local installer is run. |
| `opentargets-platform` | [Open Targets Platform official MCP](https://github.com/anthropics/life-sciences) | Select transport=http. The provider-published endpoint supports anonymous MCP initialization and tools/list as verified on 2026-09-20. No local installer is run. |
| `10x-genomics` | [10x Genomics Cloud MCP](https://github.com/10XGenomics/txg-mcp) | Separately prepare and unpack the upstream 10x Genomics MCP bundle, including its server/, node_modules/ and platform bin/ files. Set workdir to that unpacked bundle, executable_path to installed Node.js >=18, and the named TXG_CLI_ACCESS_TOKEN environment variable. DNAKit does not download a bundle or TXG binary. |
<!-- generated:mcp:end -->

</details>

## 4) 同一会话连续调用

- **作用：** 在同一 MCP 会话内依次执行多个工具，保留模型加载等服务内状态。
- **计算方法：** 每步重新读取 schema、校验显式参数并按顺序调用；遇到远端错误即停止，返回已执行步骤及错误状态。
- **API：** `dnakit.integrations.call_mcp_workflow(server_id[必须], calls[必须], transport[可选], executable_path[可选], workdir[可选], timeout_seconds[可选])`；其余连接参数见[参数 FAQ](../../faq/api/39_external_genomics.md#mcp-workflow)。
- **输入：** 服务 ID 和调用列表；每项包含实时 schema 中的 `tool_name` 与 `arguments`。下例的加载和预测工具名称、参数需由调用者传入。
- **输出：** 返回同一会话中各工具调用的结果字典。
- **示例代码：**

```python
from dnakit.integrations import call_mcp_workflow


def load_then_predict(load_tool, load_arguments, predict_tool, predict_arguments):
    calls = [
        {"tool_name": load_tool, "arguments": load_arguments},
        {"tool_name": predict_tool, "arguments": predict_arguments},
    ]
    result = call_mcp_workflow("chorus", calls, timeout_seconds=300)
    print(result["status"], result["requested_steps"], result["executed_steps"])
    return result
```

- **示例结果：** 两步均成功时的输出形式；此处未连接 Chorus，实际结果取决于服务返回。

```text
completed 2 2
```

## 5) 通过 DNAKit MCP 暴露

- **作用：** 将基因组分析和外部工具接口提供给支持 MCP 的客户端，并在调用前查询其输入与执行类型。
- **计算方法：** DNAKit 从公开函数生成工具目录、签名和 JSON schema；客户端检索目录后调用指定工具。下面离线查看同一工具目录中的运行类型。
- **API：** 启动入口 `dnakit-mcp`；MCP 工具 `search_dnakit_tools(query[必须], category[可选], limit[可选])`、`describe_dnakit_tool(tool_name[必须])`；本地目录入口 `dnakit.tools.default_tool_registry(verbose[可选])`。
- **输入：** 检索词 `genomics` 或 `integrations`，以及工具完整名称，如 `dnakit.integrations.run_analysis`。
- **输出：** 返回 DNAKit MCP 工具说明或执行结果。
- **示例代码：**

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry(verbose=False)
info = registry.describe("dnakit.integrations.run_analysis")
print(info["name"])
print(info["effect"])
```

- **示例结果：** 离线读取工具信息；未执行分析。

```text
dnakit.integrations.run_analysis
external_process
```
