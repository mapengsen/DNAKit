---
name: isabl-submit-data
description: Submit sequencing data to Isabl by creating individuals, samples, and experiments. Use when registering new data files or importing from external sources.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-submit-data` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Submitting Data to Isabl

You are helping the user submit new sequencing data to the Isabl platform.

## Checklist

Work through these steps systematically:

1. [ ] **Understand data hierarchy** (Individual → Sample → Experiment)
2. [ ] **Check required reference data** (Center, Technique, Platform, Disease)
3. [ ] **Prepare individual records**
4. [ ] **Prepare sample records**
5. [ ] **Prepare experiment records** with raw_data
6. [ ] **Validate before submission**
7. [ ] **Submit to API**
8. [ ] **Verify submission**

## Data Hierarchy

```
Individual (patient/subject)
├── identifier: unique ID from source system
├── species: HUMAN, MOUSE, etc.
├── gender: MALE, FEMALE, UNKNOWN
├── center: where samples originated
│
└── Sample (tissue specimen)
    ├── identifier: sample ID
    ├── category: TUMOR, NORMAL, GERMLINE, etc.
    ├── disease: disease classification
    │
    └── Experiment (sequencing run)
        ├── identifier: library/run ID
        ├── technique: WGS, WES, RNA, etc.
        ├── platform: ILLUMINA, ONT, etc.
        ├── center: where sequenced
        ├── projects: [project memberships]
        └── raw_data: [file paths with metadata]
```

## Step 1: Check Reference Data

```python
import isabl_cli as ii

# Check available centers
centers = list(ii.get_instances("centers"))
print("Centers:", [c.slug for c in centers])

# Check available techniques
techniques = list(ii.get_instances("techniques"))
print("Techniques:", [(t.slug, t.method) for t in techniques])

# Check available platforms
platforms = list(ii.get_instances("platforms"))
print("Platforms:", [p.slug for p in platforms])

# Check available diseases
diseases = list(ii.get_instances("diseases"))
print("Diseases:", [d.acronym for d in diseases[:10]])

# Check target project
project = ii.get_instance("projects", PROJECT_PK)
print(f"Target project: {project.title}")
```

## Step 2: Create Reference Data (if needed)

```python
# Create a new center
center = ii.create_instance("centers", {
    "name": "Memorial Sloan Kettering",
    "acronym": "MSK",
})
print(f"Created center: {center.slug}")

# Create a new technique
technique = ii.create_instance("techniques", {
    "name": "WGS 30X",
    "method": "WGS",
    "category": "DNA",
})
print(f"Created technique: {technique.slug}")

# Create a new disease
disease = ii.create_instance("diseases", {
    "name": "Acute Myeloid Leukemia",
    "acronym": "AML",
})
print(f"Created disease: {disease.slug}")
```

## Step 3: Prepare Submission Data

```python
# Example: submitting tumor-normal pair

submission_data = {
    "individual": {
        "identifier": "PATIENT_001",
        "species": "HUMAN",
        "gender": "FEMALE",
        "center": {"slug": "msk"},
    },
    "samples": [
        {
            "identifier": "PATIENT_001_TUMOR",
            "category": "TUMOR",
            "disease": {"acronym": "AML"},
        },
        {
            "identifier": "PATIENT_001_NORMAL",
            "category": "NORMAL",
            "disease": {"acronym": "AML"},
        },
    ],
    "experiments": [
        {
            "sample_identifier": "PATIENT_001_TUMOR",
            "identifier": "LIB_001_T_WGS",
            "technique": {"slug": "wgs-30x"},
            "platform": {"slug": "illumina-novaseq"},
            "center": {"slug": "msk"},
            "raw_data": [
                {
                    "file_url": "/data/fastq/LIB_001_T_R1.fastq.gz",
                    "file_type": "FASTQ",
                    "file_data": {"read": "R1", "lane": "L001"},
                },
                {
                    "file_url": "/data/fastq/LIB_001_T_R2.fastq.gz",
                    "file_type": "FASTQ",
                    "file_data": {"read": "R2", "lane": "L001"},
                },
            ],
        },
        {
            "sample_identifier": "PATIENT_001_NORMAL",
            "identifier": "LIB_001_N_WGS",
            "technique": {"slug": "wgs-30x"},
            "platform": {"slug": "illumina-novaseq"},
            "center": {"slug": "msk"},
            "raw_data": [
                {
                    "file_url": "/data/fastq/LIB_001_N_R1.fastq.gz",
                    "file_type": "FASTQ",
                    "file_data": {"read": "R1", "lane": "L001"},
                },
                {
                    "file_url": "/data/fastq/LIB_001_N_R2.fastq.gz",
                    "file_type": "FASTQ",
                    "file_data": {"read": "R2", "lane": "L001"},
                },
            ],
        },
    ],
}
```

## Step 4: Validate Before Submission

```python
import os

# Check all raw_data files exist
errors = []

for exp in submission_data["experiments"]:
    for rd in exp.get("raw_data", []):
        path = rd["file_url"]
        if not os.path.exists(path):
            errors.append(f"Missing: {path}")

if errors:
    print("Validation errors:")
    for e in errors:
        print(f"  {e}")
else:
    print("All files validated!")

# Check identifiers don't already exist
existing = list(ii.get_experiments(
    identifier__in=[e["identifier"] for e in submission_data["experiments"]]
))
if existing:
    print(f"Warning: {len(existing)} experiments already exist")
```

## Step 5: Submit Using API

### Option A: Create individually (more control)

```python
# 1. Create individual
individual = ii.create_instance("individuals", {
    "identifier": submission_data["individual"]["identifier"],
    "species": submission_data["individual"]["species"],
    "gender": submission_data["individual"]["gender"],
    "center": submission_data["individual"]["center"],
})
print(f"Created individual: {individual.system_id}")

# 2. Create samples
samples = {}
for sample_data in submission_data["samples"]:
    sample = ii.create_instance("samples", {
        "identifier": sample_data["identifier"],
        "category": sample_data["category"],
        "disease": sample_data["disease"],
        "individual": {"pk": individual.pk},
    })
    samples[sample_data["identifier"]] = sample
    print(f"Created sample: {sample.system_id}")

# 3. Create experiments
for exp_data in submission_data["experiments"]:
    sample = samples[exp_data["sample_identifier"]]
    experiment = ii.create_instance("experiments", {
        "identifier": exp_data["identifier"],
        "sample": {"pk": sample.pk},
        "technique": exp_data["technique"],
        "platform": exp_data["platform"],
        "center": exp_data["center"],
        "raw_data": exp_data["raw_data"],
        "projects": [{"pk": PROJECT_PK}],
    })
    print(f"Created experiment: {experiment.system_id}")
```

### Option B: Use nested creation (simpler)

```python
# Create experiment with nested individual/sample
# API will get-or-create the hierarchy

experiment = ii.create_instance("experiments", {
    "identifier": "LIB_001_T_WGS",
    "sample": {
        "identifier": "PATIENT_001_TUMOR",
        "category": "TUMOR",
        "disease": {"acronym": "AML"},
        "individual": {
            "identifier": "PATIENT_001",
            "species": "HUMAN",
            "gender": "FEMALE",
            "center": {"slug": "msk"},
        },
    },
    "technique": {"slug": "wgs-30x"},
    "platform": {"slug": "illumina-novaseq"},
    "center": {"slug": "msk"},
    "projects": [{"pk": PROJECT_PK}],
    "raw_data": [
        {"file_url": "/data/fastq/R1.fastq.gz", "file_type": "FASTQ"},
        {"file_url": "/data/fastq/R2.fastq.gz", "file_type": "FASTQ"},
    ],
})
```

### Option C: Use submission form (bulk)

```python
# Create submission form for bulk import
submission = ii.create_instance("submissions", {
    "title": "Batch import 2024-01",
    "data": {
        "experiments": [...],  # List of experiment dicts
    },
})

# Process the submission
ii.process_submission(submission.pk, commit=True)
```

## Step 6: Verify Submission

```python
# Check experiments were created
experiments = list(ii.get_experiments(
    projects=PROJECT_PK,
    identifier__startswith="LIB_001"
))

print(f"Created {len(experiments)} experiments:")
for exp in experiments:
    print(f"  {exp.system_id}")
    print(f"    Sample: {exp.sample.identifier} ({exp.sample.category})")
    print(f"    Individual: {exp.sample.individual.identifier}")
    print(f"    Raw data files: {len(exp.raw_data or [])}")
```

## Raw Data Format

The `raw_data` field is a list of file records:

```python
raw_data = [
    {
        "file_url": "/path/to/file.fastq.gz",  # Required
        "file_type": "FASTQ",                   # FASTQ, BAM, CRAM, etc.
        "file_data": {                          # Optional metadata
            "read": "R1",                       # R1 or R2 for paired
            "lane": "L001",                     # Flowcell lane
            "PU": "FLOWCELL:LANE:SAMPLE",       # Platform unit
            "LB": "library_id",                 # Library
            "PL": "ILLUMINA",                   # Platform
        },
        "hash_value": "abc123...",              # Optional: file checksum
        "hash_method": "MD5",                   # MD5, SHA256, etc.
    }
]
```

## Sample Categories

| Category | Description |
|----------|-------------|
| `TUMOR` | Tumor tissue |
| `NORMAL` | Normal tissue (for somatic calling) |
| `GERMLINE` | Germline sample |
| `PRIMARY` | Primary tumor |
| `METASTASIS` | Metastatic site |
| `RELAPSE` | Relapse sample |
| `XENOGRAFT` | PDX model |
| `UNKNOWN` | Unknown category |

## Common Patterns

### Import from CSV

```python
import pandas as pd

df = pd.read_csv("samples.csv")

for _, row in df.iterrows():
    experiment = ii.create_instance("experiments", {
        "identifier": row["library_id"],
        "sample": {
            "identifier": row["sample_id"],
            "category": row["category"],
            "disease": {"acronym": row["disease"]},
            "individual": {
                "identifier": row["patient_id"],
                "species": "HUMAN",
                "center": {"slug": "msk"},
            },
        },
        "technique": {"slug": row["technique"]},
        "platform": {"slug": row["platform"]},
        "center": {"slug": "msk"},
        "projects": [{"pk": PROJECT_PK}],
        "raw_data": [
            {"file_url": row["r1_path"], "file_type": "FASTQ"},
            {"file_url": row["r2_path"], "file_type": "FASTQ"},
        ],
    })
    print(f"Created: {experiment.system_id}")
```

### Update existing experiment

```python
# Add raw_data to existing experiment
exp = ii.get_instance("experiments", system_id="ISB_H000001_T01_WGS01")

ii.patch_instance("experiments", exp.pk, {
    "raw_data": [
        {"file_url": "/new/path/R1.fastq.gz", "file_type": "FASTQ"},
        {"file_url": "/new/path/R2.fastq.gz", "file_type": "FASTQ"},
    ],
})
```

### Add experiment to project

```python
# Add existing experiment to a new project
exp = ii.get_instance("experiments", system_id="ISB_H000001_T01_WGS01")
current_projects = [{"pk": p.pk} for p in exp.projects]
current_projects.append({"pk": NEW_PROJECT_PK})

ii.patch_instance("experiments", exp.pk, {
    "projects": current_projects,
})
```
