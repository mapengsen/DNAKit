---
name: isabl-query-data
description: Help construct queries to retrieve data from Isabl using the Python SDK. Use when searching for experiments, analyses, samples, or other data.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-query-data` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Querying Isabl Data

You are helping the user query data from the Isabl platform using the Python SDK.

## Checklist

Work through these steps systematically:

1. [ ] **Understand the query goal** (what data is needed)
2. [ ] **Identify the entity type** (experiments, analyses, samples, individuals)
3. [ ] **Determine filter criteria** (project, status, dates, etc.)
4. [ ] **Choose appropriate operators** (exact match, contains, ranges)
5. [ ] **Build the query** using ii.get_* functions
6. [ ] **Execute and verify results**
7. [ ] **Format output** as needed (DataFrame, JSON, etc.)
8. [ ] **When showing IDs for experiments/samples/individuals, prefer `system_id`** (or `identifier` if `system_id` is not available); avoid `pk` unless explicitly requested
9. [ ] **For analysis summaries, choose scope-aware identifiers**:
       - `individual_level_analysis` -> show individual `system_id` (e.g., `IID_H210514`)
       - `project_level_analysis` -> show project `pk`
       - avoid listing long `targets`/`references` experiment arrays unless explicitly requested
10. [ ] **When building frontend links from `ISABL_API_URL`, remove `/api/v1/`** before appending query params (e.g., `?analysis=<pk>`)
11. [ ] **Before any table, add a one-line plain-language summary** of what the table contains (scope, filters, and count if known)

## Step 1: Understand the Query Goal

Ask clarifying questions if needed:
- What entity type? (experiments, analyses, samples, individuals)
- What filters? (project, application, status, dates)
- What output format? (list, DataFrame, count only)

## Step 2: Identify Entity Type

```python
import isabl_cli as ii

# Main entity types and their query functions
experiments = ii.get_experiments(...)      # Sequencing data
analyses = ii.get_analyses(...)            # Pipeline results
samples = ii.get_instances("samples", ...) # Tissue specimens
individuals = ii.get_instances("individuals", ...)  # Patients/subjects
applications = ii.get_instances("applications", ...)  # Pipelines
projects = ii.get_instances("projects", ...)  # Project groupings
```

## Step 3: Determine Filter Criteria

```python
# Filter by project
experiments = ii.get_experiments(projects=102)

# Filter by status
analyses = ii.get_analyses(status="SUCCEEDED")

# Filter by application
analyses = ii.get_analyses(application__name="MUTECT")

# Filter by date range
analyses = ii.get_analyses(created__gte="2024-01-01")

# Combine multiple filters
experiments = ii.get_experiments(
    projects=102,
    sample__category="TUMOR",
    technique__method="WGS"
)
```

## Step 4: Choose Appropriate Operators

| Operator | Example | Purpose |
|----------|---------|---------|
| (none) | `status="SUCCEEDED"` | Exact match |
| `__contains` | `name__contains="tumor"` | Substring |
| `__icontains` | `name__icontains="TUMOR"` | Case-insensitive substring |
| `__startswith` | `name__startswith="PT"` | Starts with |
| `__in` | `status__in=["SUCCEEDED","FAILED"]` | One of values |
| `__gt`, `__gte` | `created__gt="2024-01-01"` | Greater than |
| `__lt`, `__lte` | `total__lte=100` | Less than |
| `__isnull` | `results__isnull=False` | Is/isn't null |
| `!` prefix | `status!="FAILED"` | Negation |

### Traverse Relationships

Use double underscore to traverse related objects:

```python
import isabl_cli as ii

# Experiments for a specific individual
experiments = ii.get_experiments(sample__individual__pk=500)

# Analyses by application name
analyses = ii.get_analyses(application__name="variant_caller")

# Experiments with specific technique
experiments = ii.get_experiments(technique__method="WGS")
```

## Step 5: Build the Query

### Find experiments by project

```python
import isabl_cli as ii

experiments = ii.get_experiments(
    projects=102,
    sample__category="TUMOR"
)
for exp in experiments:
    print(f"{exp.system_id}: {exp.sample.identifier}")
```

### ID display convention (important)

When users ask for an experiment/sample/individual "ID", prefer the human-readable identifier (`system_id` or `identifier`) and omit numeric `pk` unless asked.

```python
import isabl_cli as ii

experiments = ii.get_experiments(projects=102, fields=["system_id", "sample"])
for exp in experiments:
    print(f"experiment_id={exp.system_id}")
```

For MCP `isabl_query`, prefer:

```python
isabl_query(
    endpoint="experiments",
    filters={"projects": 102},
    output_fields=["system_id"]
)
```

### Analysis scope convention (important)

When listing analyses, use identifiers that match the analysis scope:

- For individual-level analyses (for example `WGTS_REPORT`), display the individual id (`IID_H...`).
- For project-level analyses, display the project `pk`.
- Do not show large `targets`/`references` experiment lists unless the user asks for them.

### Table response convention (important)

Before rendering a table, include one concise summary line that states what is being shown.

Examples:
- "Showing 3 failed WGTS_REPORT analyses in project 267, with individual IDs and stderr links."
- "Showing the latest 10 succeeded MUTECT analyses with target experiment IDs."

### Find successful analyses for an application

```python
import isabl_cli as ii

analyses = ii.get_analyses(
    application__name="variant_caller",
    application__version="2.0.0",
    status="SUCCEEDED"
)
for a in analyses:
    print(f"Analysis {a.pk}: {a.results}")
```

### Find failed analyses

```python
import isabl_cli as ii

failed = ii.get_analyses(
    status="FAILED",
    created__gt="2024-01-01"
)
for a in failed:
    print(f"{a.pk}: {a.application.name}")
```

### Get single instance by ID

```python
import isabl_cli as ii

# By system_id
exp = ii.Experiment("SAMPLE_001")

# By primary key
analysis = ii.Analysis(12345)
```

### Count records without fetching

```python
import isabl_cli as ii

count = ii.get_instances_count("experiments", projects=102)
print(f"Total experiments: {count}")
```

### Get individual's full tree

```python
import isabl_cli as ii

individual = ii.get_tree(individual_pk)
for sample in individual.sample_set:
    print(f"Sample: {sample.identifier}")
    for exp in sample.experiment_set:
        print(f"  Experiment: {exp.system_id}")
```

## Step 6: Execute and Verify Results

### Get results from analyses

```python
import isabl_cli as ii
from isabl_cli import utils

# Get specific result from an experiment
vcf_path, analysis = utils.get_result(
    experiment=experiment,
    application_key=123,
    result_key="vcf"
)
print(f"VCF: {vcf_path}")

# Get all results for an application
results = utils.get_results(
    experiment=experiment,
    application_name="variant_caller"
)
```

### Performance Tips

```python
import isabl_cli as ii

# Only fetch needed fields
experiments = ii.get_experiments(
    projects=102,
    fields=["pk", "system_id", "sample"]
)

# Stop after N results
experiments = ii.get_experiments(
    projects=102,
    count_limit=100
)

# Use cursor pagination for very large datasets
experiments = ii.get_instances(
    "experiments",
    projects=102,
    paginator="cursor"
)
```

## Step 7: Format Output

### Convert to DataFrame

```python
import isabl_cli as ii
import pandas as pd

experiments = ii.get_experiments(projects=102)
df = pd.DataFrame([
    {
        "system_id": e.system_id,
        "sample": e.sample.identifier,
        "category": e.sample.category,
        "technique": e.technique.method,
    }
    for e in experiments
])
df.head()
```

## Data Model Reference

```
Individual
├── identifier (unique)
├── species (HUMAN, MOUSE, etc.)
├── gender
└── sample_set → [Sample]

Sample
├── identifier
├── individual → Individual
├── category (TUMOR, NORMAL, etc.)
└── experiment_set → [Experiment]

Experiment
├── system_id (unique)
├── sample → Sample
├── technique → Technique
├── platform → Platform
├── raw_data (dict)
├── bam_files (dict)
└── analysis_set → [Analysis]

Analysis
├── pk
├── application → Application
├── targets → [Experiment]
├── references → [Experiment]
├── status (CREATED, STAGED, SUBMITTED, STARTED, SUCCEEDED, FAILED)
├── storage_url (path to results)
└── results (dict)

Application
├── name
├── version
├── assembly
├── species
└── settings (dict)
```
