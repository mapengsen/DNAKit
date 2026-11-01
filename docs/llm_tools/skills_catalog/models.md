# DNA 模型与加速

本类共 **15** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`alphagenome-atlas-website-links`](#alphagenome-atlas-website-links) | deepmind | Constructs deep-links and URLs for the AlphaGenome Atlas website. Supports generating single-variant exploration links (1-based chr:pos:ref&gt;alt), genomic locus… |
| [`alphagenome-single-variant-analysis`](#alphagenome-single-variant-analysis) | deepmind | Analyzes genetic variant effects on gene expression (RNA-seq), chromatin accessibility (DNASE), histone marks (ChIP), and transcription factors using the Alpha… |
| [`alphagenome-variant-impact-score`](#alphagenome-variant-impact-score) | deepmind | Score, annotate, and analyze the functional impact of genetic variants using AlphaGenome Variant Impact (AVI) scores. Query variants in chr:pos:ref&gt;alt format,… |
| [`borzoi`](#borzoi) | feynman | Use Borzoi-style regulatory genomics models for sequence-to-expression or variant-effect analysis. Use when the task asks for noncoding variant impact, regulat… |
| [`caduceus`](#caduceus) | pku-metagenomics | Caduceus — bi-directional equivariant long-range DNA language model (Mamba / RC-aware) for long genomic contexts. Use when: (1) Long DNA sequences where transf… |
| [`dnabert2`](#dnabert2) | pku-metagenomics | DNABERT-2 — efficient multi-species genome language model (BPE + ALiBi) with the Genome Understanding Evaluation (GUE) benchmark. Use when: (1) DNA embeddings… |
| [`evo2`](#evo2) | pku-metagenomics | Evo 2 — long-context DNA language model for genome modeling and design across all domains of life (Arc Institute; StripedHyena 2 / OpenGenome2). Use when: (1)… |
| [`evo2-nim`](#evo2-nim) | bionemo | Generate and analyze DNA sequences using NVIDIA's Evo 2 BioNeMo NIM microservice. Use for Evo2/Evo 2, DNA generation, genomic sequence generation, hosted gener… |
| [`genome-language-model`](#genome-language-model) | pku-metagenomics | Route genome language models (gLMs) for metagenomics and microbial genomics: representation models (Nucleotide Transformer / NTv3, DNABERT-2, Caduceus), genera… |
| [`genomic-intelligence`](#genomic-intelligence) | kdense | Predict regulatory features, gene structure, and expression directly from DNA sequence using Genomic Intelligence's hosted transformer DNA language models — no… |
| [`genomics-workflow-acceleration`](#genomics-workflow-acceleration) | bionemo | Use when accelerating existing genomics workflows with NVIDIA Parabricks, improving runtime or price/performance, converting pipeline steps to GPUs, or compari… |
| [`nucleotide-transformer`](#nucleotide-transformer) | pku-metagenomics | InstaDeep Nucleotide Transformer family — DNA foundation models from NT/NT-v2 through NTv3 (multi-species, long-context seq–function and generation). Use when:… |
| [`omg`](#omg) | pku-metagenomics | Open MetaGenomic (OMG) corpus and gLM2 mixed-modality genomic language model — large pretraining data from IMG + MGnify (Tbp-scale DNA + billions of CDS). Use… |
| [`parabricks`](#parabricks) | bionemo | Route NVIDIA Parabricks pbrun tools, assess GPU/runtime readiness, and provide version-aware command guidance for FASTQ/BAM processing, RNA-seq, variant callin… |
| [`tangermeme`](#tangermeme) | tangermeme | Use for any task involving the tangermeme library — post-training analysis of genomic sequence-to-function (S2F) deep learning models. Triggers on predictions,… |

## `alphagenome-atlas-website-links` {#alphagenome-atlas-website-links}

Constructs deep-links and URLs for the AlphaGenome Atlas website. Supports generating single-variant exploration links (1-based chr:pos:ref&gt;alt), genomic locus views (1-based closed chr:start-end), candidate summary tables, and AlphaGenome reference vs. alternate predictions. Use whenever visualizing, exploring, charting, or linking genetic variants and genomic loci on the AlphaGenome Atlas, or when asked to inspect, view, or link predictions for a genomic variant.

- **来源：**[Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/alphagenome_atlas_website_links)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alphagenome-atlas-website-links
dnakit skills install alphagenome-atlas-website-links --agent claude-code   # 或 --agent codex / cursor
```

## `alphagenome-single-variant-analysis` {#alphagenome-single-variant-analysis}

Analyzes genetic variant effects on gene expression (RNA-seq), chromatin accessibility (DNASE), histone marks (ChIP), and transcription factors using the AlphaGenome API. Use when the user asks about non-coding variant effects, pathogenicity, clinical significance, disease associations, functional effects, gene expression changes, splicing disruption, or regulatory effects in promoters and enhancers. Also use for resolving biological terms to tissue/cell-type ontologies (UBERON/CL) or analyzing variants in chr:pos:ref&gt;alt format.

- **来源：**[Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/alphagenome_single_variant_analysis)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alphagenome-single-variant-analysis
dnakit skills install alphagenome-single-variant-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `alphagenome-variant-impact-score` {#alphagenome-variant-impact-score}

Score, annotate, and analyze the functional impact of genetic variants using AlphaGenome Variant Impact (AVI) scores. Query variants in chr:pos:ref&gt;alt format, annotate VCF/tabular callsets, perform saturation mutagenesis window scans (1-based closed chr:start-end), and extract GENCODE v46 GTF gene/exon/junction coordinates all via the AlphaGenome Atlas API.

- **来源：**[Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/alphagenome_variant_impact_score)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alphagenome-variant-impact-score
dnakit skills install alphagenome-variant-impact-score --agent claude-code   # 或 --agent codex / cursor
```

## `borzoi` {#borzoi}

Use Borzoi-style regulatory genomics models for sequence-to-expression or variant-effect analysis. Use when the task asks for noncoding variant impact, regulatory sequence design, or expression prediction.

- **来源：**[companion-inc/feynman contributors @ `dfdcb7c`](https://github.com/companion-inc/feynman/tree/dfdcb7cf2c73183cff7b10aa8ea8ce370c8b152c/skills/borzoi)
- **许可证：**`MIT`
- **别名：**`alterlab-borzoi`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info borzoi
dnakit skills install borzoi --agent claude-code   # 或 --agent codex / cursor
```

## `caduceus` {#caduceus}

Caduceus — bi-directional equivariant long-range DNA language model (Mamba / RC-aware) for long genomic contexts. Use when: (1) Long DNA sequences where transformer gLMs are costly, (2) Reverse-complement equivariant embeddings (Caduceus-PS) or RC-augmented Ph models, (3) Comparing long-range DNA LMs to DNABERT-2 / NT. Upstream: https://github.com/kuleshov-group/caduceus. HF: kuleshov-group/caduceus-*. Route via genome-language-model.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/caduceus)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info caduceus
dnakit skills install caduceus --agent claude-code   # 或 --agent codex / cursor
```

## `dnabert2` {#dnabert2}

DNABERT-2 — efficient multi-species genome language model (BPE + ALiBi) with the Genome Understanding Evaluation (GUE) benchmark. Use when: (1) DNA embeddings or fine-tunes on multi-species tasks, (2) Comparing representation gLMs without NT/Evo-scale hardware, (3) GUE-style evaluation. Upstream: https://github.com/MAGICS-LAB/DNABERT_2. HF: zhihan1996/DNABERT-2-117M. Related: DNABERT (v1), DNABERT-S. Route via genome-language-model.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/dnabert2)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnabert2
dnakit skills install dnabert2 --agent claude-code   # 或 --agent codex / cursor
```

## `evo2` {#evo2}

Evo 2 — long-context DNA language model for genome modeling and design across all domains of life (Arc Institute; StripedHyena 2 / OpenGenome2). Use when: (1) Autoregressive DNA likelihoods or embeddings at up to ~1 Mb context, (2) Sequence generation / design, (3) Contrasting generative gLMs with Nucleotide Transformer or gLM2. Upstream: https://github.com/ArcInstitute/evo2. Predecessor: Evo (Science 2024). Route via genome-language-model. Heavy GPU requirements.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/evo2)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info evo2
dnakit skills install evo2 --agent claude-code   # 或 --agent codex / cursor
```

## `evo2-nim` {#evo2-nim}

Generate and analyze DNA sequences using NVIDIA's Evo 2 BioNeMo NIM microservice. Use for Evo2/Evo 2, DNA generation, genomic sequence generation, hosted generation, local Docker deployment, local forward passes, layer outputs, logits, sampled probabilities, and BioNeMo NIM workflows.

- **来源：**[NVIDIA Corporation and contributors @ `0e67a61`](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit/tree/0e67a612e4045f007e38fa77adc8f3ebfc5616b6/nim-skills/evo2-nim)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info evo2-nim
dnakit skills install evo2-nim --agent claude-code   # 或 --agent codex / cursor
```

## `genome-language-model` {#genome-language-model}

Route genome language models (gLMs) for metagenomics and microbial genomics: representation models (Nucleotide Transformer / NTv3, DNABERT-2, Caduceus), generative models (Evo / Evo 2), and mixed-modality metagenomic LMs (OMG + gLM2). Use when choosing DNA embeddings, generation, or functional-track prediction — not for classical profilers or Foldseek. Child skills: omg, nucleotide-transformer, dnabert2, caduceus, evo2. Protein LMs via protein-language-model. Trees via phylogenomics / lucaphylo.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/genome-language-model)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info genome-language-model
dnakit skills install genome-language-model --agent claude-code   # 或 --agent codex / cursor
```

## `genomic-intelligence` {#genomic-intelligence}

Predict regulatory features, gene structure, and expression directly from DNA sequence using Genomic Intelligence's hosted transformer DNA language models — no local GPU or model weights. Six tasks over a REST API and a hosted MCP server (keyless public demo): promoter regions, splice donor/acceptor sites, enhancer activity, chromatin state, sequence-to-expression (log TPM), and de-novo gene annotation, plus a composite find-genes-then-predict-expression workflow. Use when the user has a gene symbol, a genomic region, or a DNA/FASTA sequence and wants any of these predictions, mentions Genomic Intelligence, genomicintelligence.ai, api.genomicintelligence.ai, or mcp.genomicintelligence.ai.

- **来源：**[Genomic Intelligence; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/genomic-intelligence)
- **许可证：**`MIT`
- **别名：**`gi-annotation`, `gi-chromatin`, `gi-enhancer`, `gi-expression`, `gi-promoter`, `gi-splice`, `query-genomic-intelligence`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info genomic-intelligence
dnakit skills install genomic-intelligence --agent claude-code   # 或 --agent codex / cursor
```

## `genomics-workflow-acceleration` {#genomics-workflow-acceleration}

Use when accelerating existing genomics workflows with NVIDIA Parabricks, improving runtime or price/performance, converting pipeline steps to GPUs, or comparing CPU and GPU workflow outputs. Adds optional GPU steps in-place with runtime toggles (default off). Do NOT use for individual pbrun command routing — use parabricks.

- **来源：**[NVIDIA Corporation and contributors @ `0e67a61`](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit/tree/0e67a612e4045f007e38fa77adc8f3ebfc5616b6/library-skills/genomics-workflow-acceleration)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info genomics-workflow-acceleration
dnakit skills install genomics-workflow-acceleration --agent claude-code   # 或 --agent codex / cursor
```

## `nucleotide-transformer` {#nucleotide-transformer}

InstaDeep Nucleotide Transformer family — DNA foundation models from NT/NT-v2 through NTv3 (multi-species, long-context seq–function and generation). Use when: (1) DNA embeddings or MLM scores on contigs/MAGs, (2) Functional-track / genome-annotation style prediction with NTv3, (3) Contrasting representation gLMs vs generative Evo 2 or metagenomic gLM2. Upstream: https://github.com/instadeepai/nucleotide-transformer. HF: InstaDeepAI. Route via genome-language-model. Check CC BY-NC-SA license for your use case.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/nucleotide-transformer)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info nucleotide-transformer
dnakit skills install nucleotide-transformer --agent claude-code   # 或 --agent codex / cursor
```

## `omg` {#omg}

Open MetaGenomic (OMG) corpus and gLM2 mixed-modality genomic language model — large pretraining data from IMG + MGnify (Tbp-scale DNA + billions of CDS). Use when: (1) Pretraining or evaluating genomic LMs, (2) Using gLM2 embeddings for function / coevolution / regulatory syntax, (3) Distinguishing ML corpora from study reanalysis. Hugging Face: tattabio/OMG, tattabio/gLM2_650M. Upstream: https://github.com/TattaBio/OMG. Route via genome-language-model.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/omg)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info omg
dnakit skills install omg --agent claude-code   # 或 --agent codex / cursor
```

## `parabricks` {#parabricks}

Route NVIDIA Parabricks pbrun tools, assess GPU/runtime readiness, and provide version-aware command guidance for FASTQ/BAM processing, RNA-seq, variant calling, BAM QC, and GVCF workflows. Do NOT use for inspecting or accelerating whole pipelines — use genomics-workflow-acceleration.

- **来源：**[NVIDIA Corporation and contributors @ `0e67a61`](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit/tree/0e67a612e4045f007e38fa77adc8f3ebfc5616b6/library-skills/parabricks)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info parabricks
dnakit skills install parabricks --agent claude-code   # 或 --agent codex / cursor
```

## `tangermeme` {#tangermeme}

Use for any task involving the tangermeme library — post-training analysis of genomic sequence-to-function (S2F) deep learning models. Triggers on predictions, DeepLIFT/SHAP attributions, marginalization/ablation/spacing of motifs, saturation mutagenesis (ISM), variant effect scoring, sequence design, seqlet calling, loading loci/FASTA/bigWig/MEME/VCF, or wrapping a PyTorch genomic model for these analyses. This is a router skill — read the relevant file under references/ for details and footguns before writing tangermeme code.

- **来源：**[Copyright 2024-2026 Jacob Schreiber; tangermeme contributors @ `d5c51b8`](https://github.com/jmschrei/tangermeme/tree/d5c51b83d374576b64effb9284fc79a71a604734/tangermeme/_skills/data)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `tangermeme`, `torch`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info tangermeme
dnakit skills install tangermeme --agent claude-code   # 或 --agent codex / cursor
```
