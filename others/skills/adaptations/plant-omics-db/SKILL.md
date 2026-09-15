---
name: plant-omics-db
description: >-
  Build a journal-ready plant multi-omics database (genome, germplasm,
  variome, transcriptome, GWAS, tools). Use when the user wants a NAR/MP/PC
  class plant genomics web database, module list, schema, or Flask/FastAPI prototype.
---

## DNAKit execution

Use `dnakit skills doctor plant-omics-db` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# 植物多组学数据库构建 Skill

> 调用方式: 将此文件内容作为系统提示词，或在 Claude Code 中使用 `/skill` 加载。
> Playbook for building a genus-level plant multi-omics database aimed at Molecular Plant / NAR / Plant Communications.
> Distilled from shipping a grapevine (*Vitis*) database. Do not copy unpublished sample counts or private schemas from any one project.

---

## 一、启动流程

当用户说要"构建一个XX物种的多组学数据库"时，按以下顺序执行：

### Step 1: 需求摸底（5 个必问问题）
1. **数据类型与规模**: 基因组多少？重测序多少份？转录组？代谢组？表型？芯片？特殊数据（单细胞/空间/古DNA）？
2. **目标期刊**: MP (~25) / NAR (~16) / PC (~10) / HR (~8)？决定了框架复杂度
3. **团队情况**: 有网页工程师吗？有生信人员吗？还是全靠 AI 实现？
4. **现有数据**: 哪些已有？哪些需要整合公开数据？哪些需要新生成？
5. **差异化卖点**: 和同类数据库比，你的独特优势是什么？

### Step 2: 框架设计
- 输出一份 `frame_work.md`（参考已发表植物数据库的模块拆法）
- 列出所有模块、每个模块的数据内容、交互功能、参考来源
- 设计架构图（ASCII art）
- 明确"王牌模块"（2-3 个独有/碾压级优势）

### Step 3: 任务拆分
- 创建 `tasks/` 目录，每个模块一个 md 文件
- 写给网页工程师看：数据格式、API 端点、可视化类型、参考网站
- 写给生信人员看：SQL 建表、数据预处理 pipeline、文件格式要求

### Step 4: Schema 设计
- 输出 `schema/` 目录：完整建表 SQL + 数据准备指南（CSV 模板）+ 跨表联动文档
- 核心表: species, assembly, gene, gene_function, accession, variant, genotype, phenotype_value, expression_value, metabolite, gwas_result

### Step 5: 原型搭建
- Flask/FastAPI + SQLite (开发) / PostgreSQL (生产)
- 复古学术风格 CSS（复古学术风，主色可换）
- ECharts 本地化（1MB，无 CDN 依赖）
- JBrowse2 嵌入

### Step 6: 文档同步
- 每次代码改动，同步更新 frame_work.md 和 tasks/*.md

---

## 二、数据库标准模块清单

### 必做模块（所有数据库都需要）
1. **Genome** — 多基因组列表 + 基因浏览器 + 泛基因组 PAV + 注释页
2. **Variome** — 变异查询 + 群体遗传学指标 + PCA/ADMIXTURE  
3. **Accessions** — 品种卡片（每份种质一个详情页）
4. **Transcriptome** — 表达热图 + 共表达网络
5. **GWAS / Multi-Omics** — Manhattan plot + 单倍型 + 选择信号
6. **Tools** — BLAST + Primer3 + GO/KEGG + SeqFetch
7. **Download / API**
8. **Search** — 跨模块全局搜索（NCBI 风格）

### 加分模块（有数据就做，没数据预留）
9. **Epigenome** — ATAC-seq / ChIP-seq / WGBS tracks
10. **Metabolome** — 代谢物分类浏览 + mGWAS
11. **Phenotype** — 性状定义 + 表型值 + 品种对比
12. **VitisChip / SNP Array** — 芯片页面（QC指标 + 标记列表 + 填充服务）
13. **Smart Breeding** — 基因组预测 + 亲本优化
14. **Single-Cell & Spatial** — UMAP + 空间热图
15. **Ancient DNA** — 时间轴 + 损伤模式 + PCA 投影
16. **TE (Transposon)** — 超家族分布 + TIP 浏览器

### 工具集成
17. **LocusZoom** — 区域 GWAS 三面板联动
18. **JBrowse2** — 嵌入式基因组浏览器
19. **专业工具包**: Heatmap+Clustering, Volcano Plot, Phylo Tree, Network Viewer, Dot Plot, GO DAG, SeqLogo, UpSet Plot

---

## 三、数据库表结构模板（核心 15 表）

```sql
-- 元数据
species (id, scientific_name, common_name, chinese_name, taxon_id, genus, subgenus)
assembly (id, species_id, name, level, is_t2t, total_size_bp, num_genes, busco_score)
chromosome (id, assembly_id, name, length_bp)

-- 基因
gene (id, assembly_id, chromosome_id, gene_id, gene_name, gene_symbol, gene_type,
      start_pos, end_pos, strand, cds_length_bp, exon_count, protein_length_aa,
      protein_sequence, gc_content, pav_category, is_transcription_factor, tf_family,
      subcellular_local, description, mature_fruit_tpm, leaf_tpm, root_tpm,
      tissue_specificity, ortholog_arabidopsis, ortholog_rice)
gene_function (id, gene_id, annotation_db, accession, term_name, term_type, evidence_code)
gene_alias (id, gene_id, alias, source_db)

-- 种质 & 变异
accession (id, accession_name, accession_code, species_id, country, region,
           category, usage_type, domestication_status, dual_origin_group,
           clonal_group_id, is_core_accession, berry_color, description)
variant (id, assembly_id, chromosome_id, rs_id, pos, ref_allele, alt_allele,
         variant_type, consequence, impact, gene_id, gene_name, global_maf)
genotype (variant_id, accession_id, genotype) -- 分区表, 数亿行

-- 表型 & 表达
phenotype_trait (id, trait_name, trait_name_cn, trait_category, unit, data_type)
phenotype_value (id, trait_id, accession_id, value, environment, year, replicate)
expression_value (id, gene_id, sample_name, tpm, fpkm, tissue, developmental_stage)

-- GWAS
gwas_result (id, trait_id, variant_id, chromosome_id, pos, p_value, neg_log10_p,
             effect_size_beta, maf, model, environment, nearest_gene_id, significant)
```

---

## 四、前端架构模式

### CSS 配色（葡萄学术风 — 可替换物种主色调）
```css
--bg:           #FAF8F4;      /* 暖羊皮纸 */
--bg-header:    #3E1A2E;      /* 深酒红（可换成物种主题色）*/
--accent:       #4A235A;      /* 葡萄紫 */
--link:         #5C2A6E;
--success:      #2D5A27;      /* 葡萄叶绿 */
--warning:      #8B1A1A;      /* 酒红警告 */
--font-mono:    "SF Mono", "Courier New", monospace;
--font-serif:   "Georgia", "Times New Roman", serif;
```

### ECharts 加载策略
```html
<!-- 本地 1MB, 无 CDN 依赖, 国内可用 -->
<script src="/static/js/echarts.min.js"></script>
```

### 图表模板（最常用 5 种）
```javascript
// 1. 热图 (Expression)
echarts.init(dom).setOption({
  visualMap: {min:0, max:200, inRange:{color:['#FAF8F4','#F5DEB3','#D35400','#6B2D5B']}},
  series: [{type:'heatmap', data: heatData, label:{show:true}}]
});
// 2. Manhattan (GWAS)
echarts.init(dom).setOption({
  series: [{type:'scatter', data: points, large:true,
    itemStyle:{color:function(p){return p.data[2]>8?'#8B1A1A':'#4A235A';}}}]
});
// 3. PCA scatter
echarts.init(dom).setOption({
  series: groups.map(g=>({name:g.name, type:'scatter', data:g.points, itemStyle:{color:g.color}}))
});
// 4. Bar (多指标对比)  
echarts.init(dom).setOption({
  series: [{type:'bar', data: vals, barMaxWidth:20}]
});
// 5. Pie (组成)
echarts.init(dom).setOption({
  series: [{type:'pie', radius:['35%','65%'], data: pieData}]
});
```

---

## 五、后端架构

```python
# Flask 路由模式
@app.route('/gene/<gene_id>')    # Gene Card (核心页面)
@app.route('/accession/<id>')    # 品种卡片
@app.route('/variome')           # 变异浏览器
@app.route('/gwas')              # GWAS
@app.route('/search')            # 跨模块搜索 (6 表联动)
@app.route('/tools')             # 工具列表
@app.route('/upload', methods=['GET','POST'])  # 数据上传+校验
```

### 跨模块搜索路由模式
```python
@app.route('/search')
def search():
    q = request.args.get('q','')
    kw = f'%{q}%'
    results = {
        'genes': db.execute("SELECT ... FROM gene WHERE ... LIKE ? LIMIT 40", (kw,)),
        'accessions': db.execute("SELECT ... FROM accession WHERE ... LIKE ? LIMIT 30", (kw,)),
        'variants': db.execute("SELECT ... FROM variant WHERE ... LIKE ? LIMIT 30", (kw,)),
        'gwas': db.execute("SELECT ... FROM gwas_result WHERE ... LIKE ? LIMIT 15", (kw,)),
        'traits': db.execute("SELECT ... FROM phenotype_trait WHERE ... LIKE ? LIMIT 15", (kw,)),
        'metabolites': db.execute("SELECT ... FROM metabolite WHERE ... LIKE ? LIMIT 15", (kw,)),
    }
    return render_template('search.html', results=results)
```

---

## 六、部署清单

### 必需文件
```
project/
  frame_work.md              # 学术框架设计
  tasks/                     # 开发任务文档
    00_README.md
    module_01_genome.md
    module_02_variome.md
    ...
  schema/
    full_schema.sql          # 建表 SQL
    data_preparation_guide.md # CSV 模板 + 填写说明
    cross_table_linkage.md   # 表关联规则 + 导入顺序
  deploy/
    docker-compose.yml       # PostgreSQL + Redis + ES + MinIO + Nginx
    Dockerfile.backend       # FastAPI + samtools + blast + snpEff
    nginx.conf
    TOOLS_DEPLOY.md          # 30+ 工具安装指南
    DEPLOY.md                # 部署手册
  processor/
    auto_process.py          # 原始数据自动处理器
  app/
    app.py                   # Flask/FastAPI 主程序
    static/css/style.css     # 学术风 CSS（主色可换）
    static/js/echarts.min.js # ECharts 本地 (1MB)
    templates/               # Jinja2 模板 (20+ 页面)
    templates_data/          # 可下载的 CSV 模板 (10 种)
```

### Docker 一键启动
```bash
cd deploy && docker compose up -d
# PostgreSQL:5432, Redis:6379, ES:9200, MinIO:9000, Backend:8000, Nginx:80
```

---

## 七、关键设计原则（写进 frame_work.md）

1. **不设数据上限** — 自产不够就整合公开数据，表结构先建好
2. **功能可见，不极简** — 重要功能多处入口，冗余入口 > 找不到入口
3. **每个品种一张卡片** — 不是表格行，是带照片/表型/群体的完整页面
4. **跨模块联动** — 搜基因→自动显示其变异/表达/GWAS/育种信息
5. **图表本地化** — ECharts 1MB 放本地，不依赖任何 CDN
6. **MAF 不过滤** — Reference panel 构建时不过滤稀有变异
7. **工具要全** — BLAST/Primer3/CRISPR/GO/LD/LocusZoom/JBrowse2/Heatmap/Volcano/PhyloTree/Network/DotPlot/GO DAG/SeqLogo/UpSet
8. **批次溯源** — 每个数据文件记录来源/处理pipeline/参考版本

---

## 八、参考数据库速查表

| 数据库 | 物种 | 期刊 | 年份 | 模块数 | 核心借鉴 |
|--------|------|------|------|--------|---------|
| TomOmics | 番茄 | MP | 2026 | 8 | Multi-Omics框架 (Haplotype+GWAS+eQTL) + Gene Card |
| SSBP | 大豆 | MP | 2026 | 4 | SoyRefHap + Meta-GWAS + 基因组预测 + 亲本优化 |
| BnIR | 油菜 | MP | 2023 | 6 | GWAS-eQTL-TWAS-SMR-COLOC 五重关联 |
| SoyMD | 大豆 | NAR | 2024 | 5 | 38基因组 + DL变异效应预测 + 14在线工具 |
| iFish | 88鱼类 | NAR | 2026 | 6 | 共表达(2.89亿对) + TF调控网络 + 单细胞 |
| AMIR | 菊科 | NAR | 2025 | 科级 | 74物种132基因组 |
| Gramene | 139物种 | NAR | 2026 | 综合 | 葡萄泛基因组门户 + rsID标准化 |
| DeepPGDB | 50+植物 | PC | 2025 | AI | LLM+RAG+Summarize自然语言查询 |
| OrchidMD | 兰科 | PBJ | 2025 | 5 | CNN品种识别(86.5%准确率) |
| WheatOmics | 小麦 | MP | 2021 | 综合 | GeneHub基因卡片 + Mutants Track |
| Vitis super-pangenome | 葡萄 | NG | 2025 | — | 72基因组 + graph pangenome + SV-eQTL |
| Rice3KGS | 水稻 | PC | 2025 | — | 8模型 + LGBMY + SNP/SV/PAV多标记 |
| SoyDNGP | 大豆 | Brief Bioinform | 2023 | — | 3D CNN web GP服务器 |

---

## 九、常用命令速查

```bash
# 初始化数据库
python3 app.py init

# 启动开发服务器
python3 app.py                    # Flask
uvicorn app:app --port 8000       # FastAPI

# 原始数据处理
python3 processor/auto_process.py --input /raw/files/ --batch BATCH_001 --lab "LabName"

# Docker
docker compose up -d
docker compose exec backend python -c "import app; app.init_db()"

# VCF 预处理
bgzip variants.vcf && tabix -p vcf variants.vcf.gz

# 基因组索引
samtools faidx genome.fa
bwa index genome.fa

# BLAST 数据库
makeblastdb -in genome.fa -dbtype nucl

# 大文件导入 PostgreSQL
\copy gene FROM '/data/genes.csv' CSV HEADER;
```
