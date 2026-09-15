---
name: isabl-project-report
description: Generate a status report for an Isabl project. Use when summarizing project progress, identifying issues, or preparing updates.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-project-report` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Generating a Project Report

You are helping generate a status report for an Isabl project.

## Checklist

Work through these steps systematically:

1. [ ] **Get project info** (title, description, team)
2. [ ] **Count experiments** by technique and category
3. [ ] **Count analyses** by status and application
4. [ ] **Identify failures** that need attention
5. [ ] **Calculate storage** usage
6. [ ] **Generate summary** with key metrics
7. [ ] **List action items** if any

## Step 1: Get Project Info

```python
import isabl_cli as ii

project = ii.get_instance("projects", PROJECT_PK)

print(f"Project: {project.title}")
print(f"Short title: {project.short_title}")
print(f"PI: {project.principal_investigator}")
print(f"Analyst: {project.analyst}")
print(f"Description: {project.description}")
```

## Step 2: Count Experiments

```python
from collections import Counter

experiments = ii.get_experiments(projects=PROJECT_PK)

# By technique
techniques = Counter(e.technique.method for e in experiments)
print("Experiments by technique:")
for tech, count in techniques.most_common():
    print(f"  {tech}: {count}")

# By sample category
categories = Counter(e.sample.category for e in experiments)
print("\nExperiments by category:")
for cat, count in categories.most_common():
    print(f"  {cat}: {count}")

print(f"\nTotal experiments: {len(experiments)}")
```

## Step 3: Count Analyses

```python
analyses = ii.get_analyses(projects=PROJECT_PK)

# By status
statuses = Counter(a.status for a in analyses)
print("Analyses by status:")
for status in ["SUCCEEDED", "FAILED", "STARTED", "STAGED", "CREATED"]:
    count = statuses.get(status, 0)
    print(f"  {status}: {count}")

# By application
apps = Counter(a.application.name for a in analyses)
print("\nAnalyses by application (top 10):")
for app, count in apps.most_common(10):
    print(f"  {app}: {count}")

print(f"\nTotal analyses: {len(analyses)}")
```

## Step 4: Identify Failures

```python
failed = ii.get_analyses(
    projects=PROJECT_PK,
    status="FAILED"
)

if failed:
    print(f"FAILED analyses ({len(failed)}):")
    for a in failed[:10]:  # Show first 10
        target = a.targets[0].system_id if a.targets else "N/A"
        print(f"  [{a.pk}] {a.application.name}: {target}")

    if len(failed) > 10:
        print(f"  ... and {len(failed) - 10} more")
else:
    print("No failed analyses!")
```

## Step 5: Calculate Storage

```python
# Get storage usage from project
print(f"Project storage: {project.storage_usage / 1e9:.2f} GB")

# Or calculate from analyses
total_storage = sum(a.storage_usage or 0 for a in analyses)
print(f"Total analysis storage: {total_storage / 1e9:.2f} GB")
```

## Step 6: Generate Summary

Create a summary like this:

```markdown
# Project Report: {project.title}

**Date**: {today}
**PI**: {project.principal_investigator}
**Analyst**: {project.analyst}

## Overview

| Metric | Count |
|--------|-------|
| Individuals | X |
| Samples | X |
| Experiments | X |
| Analyses | X |

## Analysis Status

| Status | Count | % |
|--------|-------|---|
| SUCCEEDED | X | X% |
| FAILED | X | X% |
| IN PROGRESS | X | X% |

## Top Applications

| Application | Succeeded | Failed |
|-------------|-----------|--------|
| MUTECT | X | X |
| BATTENBERG | X | X |

## Issues Requiring Attention

- X failed analyses need investigation
- [List specific failures if any]

## Storage

Total: X.XX GB
```

## Step 7: Action Items

Based on the report, identify:

- [ ] Failed analyses to re-run or debug
- [ ] Missing analyses (samples without expected apps)
- [ ] Stale analyses (STARTED for too long)

## Common Patterns

### Find stale analyses (running too long)

```python
from datetime import datetime, timedelta

stale_cutoff = datetime.now() - timedelta(days=7)

stale = ii.get_analyses(
    projects=PROJECT_PK,
    status="STARTED",
    modified__lt=stale_cutoff.isoformat()
)

if stale:
    print(f"Stale analyses (started > 7 days ago): {len(stale)}")
```

### Find samples without a specific app

```python
# All experiments in project
all_exps = set(e.pk for e in ii.get_experiments(projects=PROJECT_PK))

# Experiments with MUTECT analysis
mutect_exps = set()
for a in ii.get_analyses(projects=PROJECT_PK, application__name="MUTECT"):
    mutect_exps.update(t.pk for t in a.targets)

# Missing
missing = all_exps - mutect_exps
print(f"Experiments without MUTECT: {len(missing)}")
```

### Export to CSV

```python
import pandas as pd

data = []
for a in analyses:
    data.append({
        "pk": a.pk,
        "application": a.application.name,
        "status": a.status,
        "target": a.targets[0].system_id if a.targets else None,
        "created": a.created,
    })

df = pd.DataFrame(data)
df.to_csv("project_report.csv", index=False)
```
