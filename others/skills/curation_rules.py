"""Reviewed scope and functional-equivalence rules for the first DNA skill snapshot.

Changing an upstream revision requires reviewing these rules and the generated decisions.
Names alone never establish scientific equivalence; groups below describe the same task
at the same workflow stage. Composite workflows and their individual stages are distinct.
"""

from __future__ import annotations

CATEGORIES = {
    "sequence": {"zh": "DNA 序列处理", "en": "DNA sequence processing"},
    "sequencing": {"zh": "测序与质控", "en": "Sequencing and quality control"},
    "alignment": {"zh": "比对与比对文件", "en": "Alignment and alignment files"},
    "intervals": {"zh": "基因组区间", "en": "Genomic intervals"},
    "assembly": {"zh": "基因组组装", "en": "Genome assembly"},
    "annotation": {"zh": "基因组注释", "en": "Genome annotation"},
    "variants": {"zh": "变异检测与处理", "en": "Variant calling and processing"},
    "clinical": {"zh": "遗传变异与疾病研究", "en": "Variant interpretation and disease research"},
    "gwas": {"zh": "GWAS 与因果遗传学", "en": "GWAS and causal genetics"},
    "population": {"zh": "群体与进化遗传学", "en": "Population and evolutionary genetics"},
    "epigenomics": {"zh": "调控与表观基因组", "en": "Regulatory and epigenomics"},
    "crispr": {"zh": "CRISPR 与分子设计", "en": "CRISPR and molecular design"},
    "comparative": {"zh": "比较基因组与系统发育", "en": "Comparative genomics and phylogenetics"},
    "metagenomics": {"zh": "微生物与环境基因组", "en": "Microbial and environmental genomics"},
    "databases": {"zh": "基因组数据库", "en": "Genomic databases"},
    "models": {"zh": "DNA 模型与加速", "en": "DNA models and acceleration"},
    "workflows": {"zh": "组合流程与技能入口", "en": "Composite workflows and skill discovery"},
}

QUERY_ALIASES = {
    "变异注释": "variant annotation",
    "变异检测": "variant calling",
    "变异解读": "variant interpretation",
    "结构变异": "structural variant",
    "拷贝数": "copy number",
    "胚系": "germline",
    "体细胞": "somatic",
    "变异": "variant",
    "突变": "variant",
    "基因注释": "gene annotation",
    "基因组注释": "genome annotation",
    "序列检索": "sequence retrieval",
    "序列获取": "sequence retrieval",
    "序列下载": "sequence retrieval",
    "序列处理": "sequence",
    "序列": "sequence",
    "反向互补": "reverse complement",
    "格式转换": "format conversion",
    "参考基因组": "reference genome",
    "基因组组装": "genome assembly",
    "组装": "assembly",
    "短读长比对": "bwa alignment",
    "长读长比对": "long read alignment",
    "比对": "alignment",
    "质量控制": "quality control",
    "质量报告": "quality reports",
    "质控": "quality",
    "接头去除": "adapter trimming",
    "引物设计": "primer design",
    "引物": "primer",
    "甲基化": "methylation",
    "表观遗传": "epigenomics",
    "表观基因组": "epigenomics",
    "染色质": "chromatin",
    "峰检测": "peak calling",
    "峰注释": "peak annotation",
    "群体遗传": "population genetics",
    "群体结构": "population structure",
    "系统发育": "phylogenetic tree",
    "比较基因组": "comparative genomics",
    "泛基因组": "pangenome",
    "宏基因组": "metagenomics",
    "微生物组": "microbiome",
    "环境DNA": "edna",
    "环境dna": "edna",
    "古DNA": "ancient dna",
    "古dna": "ancient dna",
    "古基因组": "ancient dna",
    "DNA模型": "dna model",
    "dna模型": "dna model",
    "基因组大模型": "genome language model",
    "基因组语言模型": "genome language model",
    "质粒": "plasmid",
    "整合子": "integron",
    "插入序列": "insertion sequence",
    "转座子": "transposon",
    "移动遗传元件": "mobile genetic element",
    "病毒基因组": "viral genome",
    "噬菌体": "phage",
    "谱系追踪": "lineage tracing",
    "遗传咨询": "genetic counseling",
    "家系图": "pedigree",
    "多倍体": "polyploid",
    "组装基因组质控": "mag quality",
    "单倍型": "haplotype",
    "基因型填补": "genotype imputation",
    "填补": "imputation",
    "精细定位": "fine mapping",
    "孟德尔随机化": "mendelian randomization",
    "多基因风险": "polygenic risk",
    "遗传风险": "polygenic risk",
    "关联分析": "association",
    "遗传力": "heritability",
    "共定位": "colocalization",
    "基因编辑": "genome editing",
    "脱靶": "off target",
    "限制酶": "restriction enzyme",
    "罕见病": "rare disease",
    "药物基因组": "pharmacogenomics",
    "癌症": "cancer",
    "肿瘤": "tumor",
    "基因组区间": "genomic interval",
    "坐标转换": "genomic coordinates",
    "免疫受体": "immune repertoire",
    "测序": "sequencing",
    "检索": "search",
    "分析": "analysis",
    "流程": "pipeline",
}

BIO_DOMAINS = {
    "alignment-files": "alignment",
    "alignment": "alignment",
    "atac-seq": "epigenomics",
    "causal-genomics": "gwas",
    "chip-seq": "epigenomics",
    "clinical-databases": "clinical",
    "comparative-genomics": "comparative",
    "copy-number": "variants",
    "crispr-screens": "crispr",
    "database-access": "databases",
    "ecological-genomics": "metagenomics",
    "epidemiological-genomics": "metagenomics",
    "genome-annotation": "annotation",
    "genome-assembly": "assembly",
    "genome-engineering": "crispr",
    "genome-intervals": "intervals",
    "hi-c-analysis": "epigenomics",
    "liquid-biopsy": "clinical",
    "long-read-sequencing": "sequencing",
    "metagenomics": "metagenomics",
    "methylation-analysis": "epigenomics",
    "microbiome": "metagenomics",
    "phasing-imputation": "population",
    "phylogenetics": "comparative",
    "population-genetics": "population",
    "primer-design": "crispr",
    "read-alignment": "alignment",
    "read-qc": "sequencing",
    "restriction-analysis": "crispr",
    "sequence-io": "sequence",
    "sequence-manipulation": "sequence",
    "tcr-bcr-analysis": "clinical",
    "workflow-management": "workflows",
    "variant-calling": "variants",
}
BIO_EXCLUDE = {
    "alignment/structural-alignment",
    "read-alignment/hisat2-alignment",
    "read-alignment/star-alignment",
    "read-qc/rnaseq-qc",
    "long-read-sequencing/isoseq-analysis",
    "database-access/remote-homology",
    "database-access/uniprot-access",
    "database-access/interaction-databases",
}
BIO_SUPPORT = {
    "single-cell/data-io": "workflows",
    "single-cell/preprocessing": "workflows",
    "single-cell/clustering": "workflows",
    "single-cell/multimodal-integration": "workflows",
    "single-cell/scatac-analysis": "epigenomics",
    "immunoinformatics/mhc-binding-prediction": "clinical",
    "immunoinformatics/mhc-class-ii-prediction": "clinical",
    "immunoinformatics/immunogenicity-scoring": "clinical",
    "immunoinformatics/epitope-prediction": "clinical",
}
BIO_EXTRA = {
    **BIO_SUPPORT,
    "gene-regulatory-networks/multiomics-grn": "epigenomics",
    "gene-regulatory-networks/scenic-regulons": "epigenomics",
    "immunoinformatics/neoantigen-prediction": "clinical",
}
BIO_WORKFLOWS = {
    "atacseq-pipeline",
    "causal-genomics-pipeline",
    "chipseq-pipeline",
    "cnv-pipeline",
    "crispr-editing-pipeline",
    "crispr-screen-pipeline",
    "edna-pipeline",
    "fastq-to-variants",
    "genome-annotation-pipeline",
    "genome-assembly-pipeline",
    "gwas-pipeline",
    "hic-pipeline",
    "liquid-biopsy-pipeline",
    "longread-sv-pipeline",
    "metagenomics-pipeline",
    "methylation-pipeline",
    "microbiome-pipeline",
    "multiome-pipeline",
    "neoantigen-pipeline",
    "outbreak-pipeline",
    "somatic-variant-pipeline",
    "tcr-pipeline",
}

TOOLUNIVERSE = {
    "sequence-retrieval": "databases",
    "sequence-analysis": "sequence",
    "fastq-qc": "sequencing",
    "primer-design": "crispr",
    "molecular-cloning": "crispr",
    "variant-analysis": "variants",
    "variant-interpretation": "clinical",
    "acmg-variant-classification": "clinical",
    "structural-variant-analysis": "clinical",
    "variant-functional-annotation": "variants",
    "cancer-variant-interpretation": "clinical",
    "cancer-genomics-tcga": "clinical",
    "rare-disease-genomics": "clinical",
    "regulatory-genomics": "epigenomics",
    "regulatory-variant-analysis": "epigenomics",
    "variant-to-mechanism": "epigenomics",
    "epigenomics": "epigenomics",
    "epigenomics-chromatin": "epigenomics",
    "gene-regulatory-networks": "epigenomics",
    "crispr-screen-analysis": "crispr",
    "functional-genomics-screens": "crispr",
    "gwas-drug-discovery": "gwas",
    "gwas-finemapping": "gwas",
    "gwas-snp-interpretation": "gwas",
    "gwas-study-explorer": "gwas",
    "gwas-trait-to-gene": "gwas",
    "mendelian-randomization": "gwas",
    "phewas": "gwas",
    "polygenic-risk-score": "gwas",
    "pathway-disease-genetics": "gwas",
    "population-genetics": "population",
    "population-genetics-1000genomes": "databases",
    "comparative-genomics": "comparative",
    "phylogenetics": "comparative",
    "plant-genomics": "comparative",
    "microbial-genome-characterization": "metagenomics",
    "metagenomics-analysis": "metagenomics",
    "model-organism-genetics": "databases",
    "gene-disease-association": "clinical",
    "hla-immunogenomics": "clinical",
    "immune-repertoire-analysis": "clinical",
    "pharmacogenomics": "clinical",
    "microbiome-research": "metagenomics",
    "variant-predictor-dms-validation": "variants",
    "rare-disease-diagnosis": "clinical",
}

CLAWBIO = {
    "analyze-fasta": "sequence",
    "ancestry-risk-profiler": "population",
    "archaic-introgression": "comparative",
    "article-data-fetcher": "databases",
    "bioqc-mcp": "sequencing",
    "busco-assessor": "assembly",
    "claw-amplicon-qc": "sequencing",
    "claw-ancestry-pca": "population",
    "claw-metagenomics": "metagenomics",
    "clinical-variant-prioritizer": "clinical",
    "clinical-variant-reporter": "clinical",
    "clinpgx": "clinical",
    "cnv-acmg-classifier": "clinical",
    "crispr-screen-triage": "crispr",
    "dnasp": "population",
    "eqtl-catalogue-region-fetch": "databases",
    "equity-scorer": "population",
    "fastreer": "comparative",
    "fine-mapping": "gwas",
    "gwas-catalog-region-fetch": "databases",
    "gwas-lookup": "gwas",
    "gwas-pipeline": "workflows",
    "gwas-prs": "gwas",
    "hla-typing": "clinical",
    "illumina-bridge": "workflows",
    "just-prs-mcp": "gwas",
    "ld-1000g-region-compute": "population",
    "locuscompare-region-render": "gwas",
    "mendelian-randomisation": "gwas",
    "methylation-clock": "epigenomics",
    "multiqc-reporter": "sequencing",
    "ncbi-datasets": "databases",
    "nfcore-sarek-wrapper": "workflows",
    "pharmgx-reporter": "clinical",
    "phylogenetics-builder": "comparative",
    "polars-bio": "intervals",
    "rare-high-impact-variants": "clinical",
    "recombinator": "population",
    "sample-qc-triage": "sequencing",
    "seq-wrangler": "workflows",
    "variant-annotation": "variants",
    "vcf-annotator": "variants",
    "wes-clinical-report-en": "clinical",
    "wes-clinical-report-es": "clinical",
    "wgs-prs": "workflows",
    **{
        f"gi-{task}": "models"
        for task in ("annotation", "chromatin", "enhancer", "expression", "promoter", "splice")
    },
}
KDENSE = {
    "biopython": "sequence",
    "bioservices": "databases",
    "depmap": "crispr",
    "dnanexus-integration": "workflows",
    "etetoolkit": "comparative",
    "folklore-variant-evidence": "clinical",
    "geniml": "intervals",
    "genomic-coordinates": "intervals",
    "genomic-intelligence": "models",
    "gget": "databases",
    "gtars": "intervals",
    "onekgpd": "databases",
    "pathogen-variant-surveillance": "metagenomics",
    "phylogenetics": "comparative",
    "polars-bio": "intervals",
    "pysam": "alignment",
    "scikit-bio": "sequence",
    "tiledbvcf": "variants",
}
SCI_EXCLUDE = {
    "star-rna-seq-aligner",
    "arboreto-grn-inference",
    "archs4-database",
    "kegg-database",
    "quickgo-database",
    "viennarna-structure-prediction",
}

# Candidate names sharing a complete task/stage. The first matching source priority wins,
# unless PREFERRED selects a specific implementation. Distinct data modalities remain separate.
EQUIVALENT = {
    "single-cell-atac-analysis": [
        "bio-atac-seq-single-cell-atac",
        "bio-single-cell-scatac-analysis",
    ],
    "fastq-quality-reports": [
        "bio-read-qc-quality-reports",
        "tooluniverse-fastq-qc",
        "bioqc-mcp",
        "genomics-qc",
    ],
    "fastq-preprocessing-fastp": ["bio-read-qc-fastp-workflow", "fastp-fastq-preprocessing"],
    "multiqc-report-aggregation": ["multiqc-qc-reports", "multiqc-reporter"],
    "short-read-wgs-alignment": ["bio-read-alignment-bwa-alignment", "bwa-mem2-dna-aligner"],
    "prokaryotic-genome-annotation": [
        "bio-genome-annotation-prokaryotic-annotation",
        "bakta-genome-annotation",
        "prokka-genome-annotation",
    ],
    "pangenome-analysis": ["bio-comparative-genomics-pangenome-analysis", "roary-pangenome"],
    "python-genomic-files": ["pysam", "pysam-genomic-files"],
    "biopython-sequence-analysis": [
        "biopython",
        "biopython-molecular-biology",
        "biopython-sequence-analysis",
        "tooluniverse-sequence-analysis",
    ],
    "vcf-annotation": [
        "bio-variant-annotation",
        "variant-annotation",
        "vcf-annotator",
        "snpeff-variant-annotation",
        "tooluniverse-variant-functional-annotation",
    ],
    "vcf-manipulation": [
        "bio-vcf-manipulation",
        "bcftools-variant-manipulation",
        "genomics-vcf-operations",
    ],
    "vcf-filtering": ["bio-variant-calling-filtering-best-practices", "vcf-variant-filtering"],
    "vcf-statistics": ["bio-vcf-statistics", "genomics-variant-calling"],
    "germline-variant-calling-gatk": ["bio-gatk-variant-calling", "gatk-variant-calling"],
    "germline-acmg-classification": [
        "bio-clinical-databases-acmg-classification",
        "tooluniverse-acmg-variant-classification",
        "clinical-variant-reporter",
        "bio-variant-calling-clinical-interpretation",
    ],
    "rare-variant-prioritization": [
        "bio-clinical-databases-variant-prioritization",
        "clinical-variant-prioritizer",
    ],
    "germline-cnv-interpretation": [
        "bio-copy-number-germline-cnv-interpretation",
        "cnv-acmg-classifier",
        "tooluniverse-structural-variant-analysis",
    ],
    "cnvkit-analysis": ["bio-copy-number-cnvkit-analysis", "cnvkit-copy-number"],
    "copy-ratio-segmentation": ["bio-copy-number-copy-ratio-segmentation", "genomics-cnv-calling"],
    "clinvar-lookup": ["bio-clinical-databases-clinvar-lookup", "clinvar-database"],
    "dbsnp-lookup": ["bio-clinical-databases-dbsnp-queries", "dbsnp-database"],
    "gnomad-frequencies": ["bio-clinical-databases-gnomad-frequencies", "gnomad-database"],
    "ensembl-api": ["bio-ensembl-rest", "ensembl-database"],
    "geo-data-retrieval": ["bio-geo-data", "geo-database"],
    "ncbi-datasets": ["bio-ncbi-datasets-cli", "ncbi-datasets"],
    "hla-typing": ["bio-clinical-databases-hla-typing", "hla-typing"],
    "pharmacogenomics": [
        "bio-clinical-databases-pharmacogenomics",
        "tooluniverse-pharmacogenomics",
        "pharmgx-reporter",
    ],
    "clinpgx-api": ["clinpgx", "clinpgx-database"],
    "polygenic-risk-score": [
        "bio-clinical-databases-polygenic-risk",
        "tooluniverse-polygenic-risk-score",
        "gwas-prs",
        "just-prs-mcp",
    ],
    "gwas-fine-mapping": [
        "bio-causal-genomics-fine-mapping",
        "tooluniverse-gwas-finemapping",
        "fine-mapping",
    ],
    "mendelian-randomization": [
        "bio-causal-genomics-mendelian-randomization",
        "tooluniverse-mendelian-randomization",
        "mendelian-randomisation",
    ],
    "gwas-association-testing": [
        "bio-population-genetics-association-testing",
        "plink2-gwas-analysis",
    ],
    "gwas-pipeline": ["bio-workflows-gwas-pipeline", "gwas-pipeline"],
    "gwas-variant-evidence": ["tooluniverse-gwas-snp-interpretation", "gwas-lookup"],
    "1000genomes-queries": ["tooluniverse-population-genetics-1000genomes", "onekgpd"],
    "population-structure": ["bio-population-genetics-population-structure", "claw-ancestry-pca"],
    "introgression-detection": [
        "bio-comparative-genomics-introgression-detection",
        "archaic-introgression",
    ],
    "population-genetics-statistics": [
        "bio-population-genetics-scikit-allel-analysis",
        "tooluniverse-population-genetics",
        "dnasp",
    ],
    "genome-assembly-qc": ["bio-genome-assembly-assembly-qc", "genomics-assembly"],
    "genomic-interval-algebra": [
        "bio-genome-intervals-interval-arithmetic",
        "bedtools-genomic-intervals",
        "gtars",
        "polars-bio",
    ],
    "genomic-interval-embedding": ["geniml"],
    "genomics-database-gget": ["gget", "gget-genomic-databases"],
    "phylogenetic-tree-manipulation": ["bio-phylo-tree-manipulation", "etetoolkit"],
    "phylogenetic-tree-inference": [
        "bio-phylo-modern-tree-inference",
        "tooluniverse-phylogenetics",
        "phylogenetics-builder",
        "phylogenetics",
    ],
    "epigenetic-clock": ["bio-methylation-epigenetic-clocks", "methylation-clock"],
    "primer-design": ["bio-primer-design-primer-basics", "tooluniverse-primer-design"],
    "grna-design": ["bio-genome-engineering-grna-design", "sgrna-design-guide"],
    "depmap-crispr-essentiality": ["depmap", "depmap-crispr-essentiality"],
    "crispr-screen-pipeline": [
        "bio-workflows-crispr-screen-pipeline",
        "tooluniverse-crispr-screen-analysis",
    ],
    "chipseq-motif-analysis": ["bio-chipseq-motif-analysis", "homer-motif-analysis"],
    "chipseq-peak-calling": ["bio-chipseq-peak-calling", "macs3-peak-calling"],
    "metagenomics-pipeline": [
        "bio-workflows-metagenomics-pipeline",
        "tooluniverse-metagenomics-analysis",
        "claw-metagenomics",
    ],
    "genomic-intelligence": [
        "genomic-intelligence",
        "gi-annotation",
        "gi-chromatin",
        "gi-enhancer",
        "gi-expression",
        "gi-promoter",
        "gi-splice",
    ],
    "evo2-nim": ["evo2-nim"],
    "parabricks": ["parabricks"],
    "genomics-workflow-acceleration": ["genomics-workflow-acceleration"],
}
PRIORITY = {
    "bioskills": 0,
    "tooluniverse": 1,
    "kdense": 2,
    "clawbio": 3,
    "sciagent": 4,
    "omicsclaw": 5,
    "bionemo": 6,
}
PREFERRED = {
    "biopython-sequence-analysis": "tooluniverse-sequence-analysis",
    "genomic-intelligence": "genomic-intelligence",
    "genomic-interval-algebra": "bio-genome-intervals-interval-arithmetic",
    "evo2-nim": "nim-skills/evo2-nim/SKILL.md",
    "parabricks": "library-skills/parabricks/SKILL.md",
    "genomics-workflow-acceleration": "library-skills/genomics-workflow-acceleration/SKILL.md",
}


def classify(source: str, path: str, name: str) -> str | None:
    """Return the reviewed task category or exclude a non-DNA/legacy/example definition."""
    stem = path.removesuffix("/SKILL.md")
    if source == "bioskills":
        if stem in BIO_EXCLUDE:
            return None
        if stem in BIO_EXTRA:
            return BIO_EXTRA[stem]
        if stem.startswith("workflows/"):
            return "workflows" if stem.split("/")[1] in BIO_WORKFLOWS else None
        return BIO_DOMAINS.get(stem.split("/")[0])
    if source == "tooluniverse":
        return (
            TOOLUNIVERSE.get(name.removeprefix("tooluniverse-")) if path.count("/") == 2 else None
        )
    if source == "clawbio":
        return CLAWBIO.get(name) if path.count("/") == 2 else None
    if source == "kdense":
        return KDENSE.get(name)
    if source == "sciagent":
        if name in SCI_EXCLUDE or any(part in path for part in ("/rnaseq/", "/single-cell/")):
            return None
        if path.startswith("skills/molecular-biology/"):
            return "crispr"
        if not path.startswith("skills/genomics-bioinformatics/"):
            return None
        if "/databases/" in path:
            return "databases"
        if "/variant/" in path:
            return "variants"
        if "/alignment/" in path:
            return "alignment"
        if "/interval-ops/" in path:
            return "intervals"
        if "/annotation/" in path:
            return "annotation"
        if "/qc/" in path:
            return "sequencing"
        return "sequence"
    if source == "omicsclaw" and path.startswith("skills/genomics/"):
        return (
            "variants"
            if "variant" in name or "cnv" in name or "vcf" in name or "sv-" in name
            else "sequencing"
        )
    if source == "bionemo" and any(
        part in path for part in ("/evo2-nim/", "/parabricks/", "/genomics-workflow-acceleration/")
    ):
        return "models"
    return None
