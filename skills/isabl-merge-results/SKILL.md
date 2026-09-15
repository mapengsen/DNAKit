---
name: isabl-merge-results
description: Aggregate results from multiple Isabl analyses into combined datasets. Use when merging VCFs, TSVs, or other outputs across samples or cohorts.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-merge-results` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Merging Results from Multiple Analyses

You are helping the user aggregate results from multiple Isabl analyses.

## Checklist

Work through these steps systematically:

1. [ ] **Identify target analyses** (by project, application, status)
2. [ ] **Understand result schema** (what files/keys are available)
3. [ ] **Collect result paths** from storage_url
4. [ ] **Verify files exist** before merging
5. [ ] **Choose merge strategy** (concat, join, custom)
6. [ ] **Execute merge** and save output
7. [ ] **Validate merged output**

## Step 1: Identify Target Analyses

```python
import isabl_cli as ii

# Get succeeded analyses for an application in a project
analyses = ii.get_analyses(
    projects=PROJECT_PK,
    application__name="MUTECT",
    status="SUCCEEDED"
)

print(f"Found {len(analyses)} analyses to merge")
for a in analyses[:5]:
    target = a.targets[0].system_id if a.targets else "N/A"
    print(f"  [{a.pk}] {target}")
```

## Step 2: Understand Result Schema

```python
# Check what results are available
sample_analysis = analyses[0]
print(f"Available results: {list(sample_analysis.results.keys())}")

# Check application's result definitions
app = ii.get_instance("applications", sample_analysis.application.pk)
print(f"Result schema: {app.results}")
```

Common result keys:
- `vcf` - Variant calls
- `bam` - Aligned reads
- `tsv` / `csv` - Tabular data
- `summary` - Summary statistics
- `qc_metrics` - Quality metrics

## Step 3: Collect Result Paths

```python
from pathlib import Path

result_key = "tsv"  # The result key to merge

paths = []
for a in analyses:
    # Option 1: From results dict (if path stored there)
    if result_key in a.results:
        path = a.results[result_key]
        if isinstance(path, dict):
            path = path.get("path") or path.get("url")
        paths.append({"analysis_pk": a.pk, "path": path})

    # Option 2: Construct from storage_url
    else:
        expected_path = Path(a.storage_url) / f"{result_key}.tsv"
        paths.append({"analysis_pk": a.pk, "path": str(expected_path)})

print(f"Collected {len(paths)} paths")
```

## Step 4: Verify Files Exist

```python
import os

valid_paths = []
missing = []

for item in paths:
    if os.path.exists(item["path"]):
        valid_paths.append(item)
    else:
        missing.append(item)

print(f"Valid: {len(valid_paths)}, Missing: {len(missing)}")

if missing:
    print("Missing files:")
    for m in missing[:5]:
        print(f"  Analysis {m['analysis_pk']}: {m['path']}")
```

## Step 5: Choose Merge Strategy

### Strategy A: Simple Concatenation (same schema)

```python
import pandas as pd

dfs = []
for item in valid_paths:
    df = pd.read_csv(item["path"], sep="\t")
    df["analysis_pk"] = item["analysis_pk"]
    dfs.append(df)

merged = pd.concat(dfs, ignore_index=True)
print(f"Merged shape: {merged.shape}")
```

### Strategy B: Wide Format (pivot by sample)

```python
dfs = []
for a in analyses:
    target_id = a.targets[0].system_id if a.targets else f"analysis_{a.pk}"
    df = pd.read_csv(a.results["tsv"], sep="\t")
    df = df.rename(columns={"value": target_id})
    dfs.append(df)

# Merge on common key
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="key_column", how="outer")
```

### Strategy C: Extract Specific Metrics

```python
# Build summary table from individual results
data = []
for a in analyses:
    target = a.targets[0] if a.targets else None
    row = {
        "analysis_pk": a.pk,
        "sample_id": target.system_id if target else None,
        "category": target.sample.category if target else None,
    }

    # Extract metrics from results
    if "summary" in a.results:
        summary = a.results["summary"]
        if isinstance(summary, dict):
            row.update(summary.get("data", {}))

    data.append(row)

merged = pd.DataFrame(data)
```

## Step 6: Execute Merge

```python
# Full example: merge mutation counts from MUTECT analyses

import pandas as pd
from pathlib import Path

analyses = ii.get_analyses(
    projects=PROJECT_PK,
    application__name="MUTECT",
    status="SUCCEEDED"
)

data = []
for a in analyses:
    target = a.targets[0]
    reference = a.references[0] if a.references else None

    # Read the mutations TSV
    muts_path = Path(a.storage_url) / "mutations.tsv"
    if muts_path.exists():
        muts = pd.read_csv(muts_path, sep="\t")

        data.append({
            "sample_id": target.system_id,
            "tumor_category": target.sample.category,
            "normal_id": reference.system_id if reference else None,
            "total_mutations": len(muts),
            "snvs": len(muts[muts["type"] == "SNV"]),
            "indels": len(muts[muts["type"] == "INDEL"]),
            "analysis_pk": a.pk,
        })

merged = pd.DataFrame(data)
merged.to_csv("mutation_summary.csv", index=False)
print(f"Saved {len(merged)} samples to mutation_summary.csv")
```

## Step 7: Validate Merged Output

```python
# Check for issues
print(f"Shape: {merged.shape}")
print(f"Columns: {list(merged.columns)}")
print(f"Missing values:\n{merged.isnull().sum()}")
print(f"Duplicates: {merged.duplicated().sum()}")

# Sample preview
print(merged.head())
```

## Common Merge Patterns

### Merge VCFs with bcftools

```python
# First, generate VCF list file from Python
vcf_paths = []
for a in analyses:
    vcf_path = Path(a.storage_url) / "output.vcf.gz"
    if vcf_path.exists():
        vcf_paths.append(str(vcf_path))

# Write to file for bcftools
with open("vcf_files.txt", "w") as f:
    f.write("\n".join(vcf_paths))

print(f"Listed {len(vcf_paths)} VCF files")
```

```bash
# Then merge with bcftools
bcftools merge -l vcf_files.txt -o merged.vcf.gz -O z
```

### Merge BAM coverage

```python
# Collect coverage metrics from multiple analyses
coverage_data = []
for a in analyses:
    metrics_file = Path(a.storage_url) / "coverage_metrics.txt"
    if metrics_file.exists():
        # Parse coverage metrics (format varies by tool)
        with open(metrics_file) as f:
            for line in f:
                if line.startswith("MEAN_COVERAGE"):
                    coverage_data.append({
                        "sample": a.targets[0].system_id,
                        "mean_coverage": float(line.split("\t")[1])
                    })
```

### Project-level aggregation

```python
# Use project-level analysis if available
project_analyses = ii.get_analyses(
    project_level_analysis=PROJECT_PK,
    application__name="COHORT_SUMMARY",
    status="SUCCEEDED"
)

if project_analyses:
    # Results already aggregated
    merged = project_analyses[0].results
else:
    # Need to aggregate manually
    pass
```

## Handling Large Datasets

```python
# Process in chunks to avoid memory issues
import pandas as pd

chunk_size = 100
output_file = "merged_large.csv"
first_chunk = True

for i in range(0, len(valid_paths), chunk_size):
    chunk_paths = valid_paths[i:i+chunk_size]

    dfs = [pd.read_csv(p["path"], sep="\t") for p in chunk_paths]
    chunk_df = pd.concat(dfs, ignore_index=True)

    chunk_df.to_csv(
        output_file,
        mode="w" if first_chunk else "a",
        header=first_chunk,
        index=False
    )
    first_chunk = False
    print(f"Processed {min(i+chunk_size, len(valid_paths))}/{len(valid_paths)}")
```

## Result File Conventions

Standard files in `analysis.storage_url`:

| File | Content |
|------|---------|
| `head_job.log` | Job stdout |
| `head_job.err` | Job stderr |
| `head_job.sh` | Job script |
| `*.vcf.gz` | Variant calls |
| `*.bam` | Aligned reads |
| `*.tsv` / `*.csv` | Tabular results |
| `summary.json` | Summary metrics |
