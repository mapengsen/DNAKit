---
name: isabl-run-pipeline
description: Submit and run multiple Isabl applications as a pipeline. Use when chaining apps together or running a sequence of analyses on samples.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-run-pipeline` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Running Isabl Pipelines

You are helping the user run multiple Isabl applications as a coordinated pipeline.

## Checklist

Work through these steps systematically:

0. [ ] **Find Available Applications** (apps by assemblies)
1. [ ] **Define pipeline steps** (apps and order)
2. [ ] **Identify target experiments**
3. [ ] **Check prerequisites** (dependencies, raw data)
4. [ ] **Stage analyses** (create without running)
5. [ ] **Submit pipeline** (commit to run)
6. [ ] **Monitor execution**
7. [ ] **Verify completion**


## Setp 0: Find Avaiable Apps

### Find Apps

Find all apps per assembly on the format `isabl apps-<assembly> --help`

```bash
isabl --help

# i.e 
isabl apps-grch37 --help
```

## Step 1: Define Pipeline Steps

### Common Pipelines

**WGS Tumor-Normal Pipeline:**
```python
PIPELINE = [
    {"app": "bwa-0-7-17", "type": "TARGETS"},           # Align tumor
    {"app": "bwa-0-7-17", "type": "REFERENCES"},        # Align normal
    {"app": "mutect-2-0", "type": "PAIRS"},             # Somatic SNVs
    {"app": "strelka-2-9", "type": "PAIRS"},            # Somatic variants
    {"app": "battenberg-2-2", "type": "PAIRS"},         # CNV calling
    {"app": "brass-6-0", "type": "PAIRS"},              # Structural variants
]
```

**RNA-Seq Pipeline:**
```python
PIPELINE = [
    {"app": "star-2-7", "type": "TARGETS"},             # Align
    {"app": "rsem-1-3", "type": "TARGETS"},             # Quantification
    {"app": "arriba-2-0", "type": "TARGETS"},           # Fusion calling
]
```

**Single Sample Pipeline:**
```python
PIPELINE = [
    {"app": "bwa-0-7-17", "type": "TARGETS"},           # Align
    {"app": "gatk-haplotype-4-0", "type": "TARGETS"},   # Germline variants
    {"app": "coverage-stats", "type": "TARGETS"},       # QC metrics
]
```

## Step 2: Identify Target Experiments

```python
import isabl_cli as ii

# Get experiments to process
PROJECT_PK = 102

# Option A: All WGS experiments in project
experiments = list(ii.get_experiments(
    projects=PROJECT_PK,
    technique__method="WGS"
))

# Option B: Specific experiments by ID
experiment_ids = ["ISB_H000001_T01_WGS01", "ISB_H000002_T01_WGS01"]
experiments = [ii.Experiment(eid) for eid in experiment_ids]

# Option C: Tumor-Normal pairs
tumors = list(ii.get_experiments(
    projects=PROJECT_PK,
    sample__category="TUMOR",
    technique__method="WGS"
))

# Find matching normals
pairs = []
for tumor in tumors:
    normals = list(ii.get_experiments(
        sample__individual=tumor.sample.individual.pk,
        sample__category="NORMAL",
        technique__method="WGS"
    ))
    if normals:
        pairs.append((tumor, normals[0]))

print(f"Found {len(pairs)} tumor-normal pairs")
```

## Step 3: Check Prerequisites

```python
from pathlib import Path

def check_prerequisites(experiments, pipeline):
    """Check if experiments are ready for pipeline."""
    issues = []

    for exp in experiments:
        # Check raw data exists
        if not exp.raw_data:
            issues.append(f"{exp.system_id}: No raw data")
            continue

        # Check files exist
        for rd in exp.raw_data:
            if not Path(rd["file_url"]).exists():
                issues.append(f"{exp.system_id}: Missing {rd['file_url']}")

    # Check required apps are registered
    for step in pipeline:
        apps = list(ii.get_instances("applications", name=step["app"].rsplit("-", 1)[0]))
        if not apps:
            issues.append(f"Application not found: {step['app']}")

    return issues

issues = check_prerequisites(experiments, PIPELINE)
if issues:
    print("Prerequisites not met:")
    for issue in issues:
        print(f"  - {issue}")
else:
    print("All prerequisites met!")
```

## Step 4: Stage Analyses (Dry Run)

```python
def stage_pipeline(experiments, pipeline, pairs=None):
    """Stage pipeline without committing."""
    staged = []

    for step in pipeline:
        app_name = step["app"]
        app_type = step["type"]

        print(f"\n=== Staging {app_name} ({app_type}) ===")

        if app_type == "TARGETS":
            # Run on each experiment
            for exp in experiments:
                cmd = f"isabl {app_name} --targets {exp.system_id}"
                print(f"  {cmd}")
                staged.append({"app": app_name, "targets": [exp.system_id]})

        elif app_type == "REFERENCES":
            # Run on reference experiments (normals)
            if pairs:
                for tumor, normal in pairs:
                    cmd = f"isabl {app_name} --targets {normal.system_id}"
                    print(f"  {cmd}")
                    staged.append({"app": app_name, "targets": [normal.system_id]})

        elif app_type == "PAIRS":
            # Run on tumor-normal pairs
            if pairs:
                for tumor, normal in pairs:
                    cmd = f"isabl {app_name} --targets {tumor.system_id} --references {normal.system_id}"
                    print(f"  {cmd}")
                    staged.append({
                        "app": app_name,
                        "targets": [tumor.system_id],
                        "references": [normal.system_id]
                    })

    return staged

staged = stage_pipeline(experiments, PIPELINE, pairs=pairs)
print(f"\nTotal analyses to stage: {len(staged)}")
```

## Step 5: Submit Pipeline

### Option A: Sequential Submission (Simple)

```python
def run_pipeline_sequential(experiments, pipeline, pairs=None, commit=False):
    """Run pipeline steps sequentially."""
    commit_flag = "--commit" if commit else ""

    for step in pipeline:
        app_name = step["app"]
        app_type = step["type"]

        print(f"\n=== Running {app_name} ===")

        if app_type == "TARGETS":
            for exp in experiments:
                cmd = f"isabl {app_name} --targets {exp.system_id} {commit_flag}"
                print(f"Executing: {cmd}")
                if commit:
                    os.system(cmd)

        elif app_type == "PAIRS" and pairs:
            for tumor, normal in pairs:
                cmd = f"isabl {app_name} --targets {tumor.system_id} --references {normal.system_id} {commit_flag}"
                print(f"Executing: {cmd}")
                if commit:
                    os.system(cmd)

# Dry run first
run_pipeline_sequential(experiments, PIPELINE, pairs=pairs, commit=False)

# Then commit
# run_pipeline_sequential(experiments, PIPELINE, pairs=pairs, commit=True)
```

### Option B: Batch Submission (Recommended)

```python
import subprocess

def run_pipeline_batch(experiments, pipeline, pairs=None, commit=False):
    """Submit all pipeline steps at once."""
    commands = []
    commit_flag = "--commit" if commit else ""

    for step in pipeline:
        app_name = step["app"]
        app_type = step["type"]

        if app_type == "TARGETS":
            # Submit all targets at once
            targets = ",".join(e.system_id for e in experiments)
            commands.append(f"isabl {app_name} --targets {targets} {commit_flag}")

        elif app_type == "PAIRS" and pairs:
            # Submit all pairs at once
            for tumor, normal in pairs:
                commands.append(
                    f"isabl {app_name} --targets {tumor.system_id} "
                    f"--references {normal.system_id} {commit_flag}"
                )

    # Execute commands
    for cmd in commands:
        print(f"Executing: {cmd}")
        if commit:
            subprocess.run(cmd, shell=True)

    return commands

commands = run_pipeline_batch(experiments, PIPELINE, pairs=pairs, commit=False)
```

### Option C: Using Python API

```python
from isabl_cli import api

def run_pipeline_api(experiments, pipeline, pairs=None, commit=False):
    """Run pipeline using Python API."""

    for step in pipeline:
        app_name = step["app"]
        app_type = step["type"]

        # Get application class
        app_class = api.get_application(app_name)
        app = app_class()

        print(f"\n=== {app_name} ===")

        if app_type == "TARGETS":
            tuples = [([exp], []) for exp in experiments]
        elif app_type == "PAIRS" and pairs:
            tuples = [([tumor], [normal]) for tumor, normal in pairs]
        else:
            continue

        # Run application
        app.run(tuples=tuples, commit=commit)

# run_pipeline_api(experiments, PIPELINE, pairs=pairs, commit=True)
```

## Step 6: Monitor Execution

```python
from datetime import datetime
import time

def monitor_pipeline(project_pk, pipeline_apps, interval=300):
    """Monitor pipeline progress."""

    while True:
        print(f"\n=== Pipeline Status ({datetime.now()}) ===")

        all_done = True
        for app_name in pipeline_apps:
            # Get analyses for this app
            analyses = list(ii.get_analyses(
                projects=project_pk,
                application__name=app_name.rsplit("-", 1)[0]
            ))

            if not analyses:
                continue

            # Count statuses
            from collections import Counter
            statuses = Counter(a.status for a in analyses)

            succeeded = statuses.get("SUCCEEDED", 0)
            failed = statuses.get("FAILED", 0)
            running = statuses.get("STARTED", 0)
            pending = statuses.get("SUBMITTED", 0) + statuses.get("STAGED", 0)

            print(f"  {app_name}: ✓{succeeded} ✗{failed} ▶{running} ⏳{pending}")

            if running > 0 or pending > 0:
                all_done = False

        if all_done:
            print("\nPipeline complete!")
            break

        time.sleep(interval)

pipeline_apps = [step["app"] for step in PIPELINE]
# monitor_pipeline(PROJECT_PK, pipeline_apps)
```

## Step 7: Verify Completion

```python
def verify_pipeline(project_pk, pipeline, experiments, pairs=None):
    """Verify all pipeline steps completed successfully."""
    results = {}

    for step in pipeline:
        app_name = step["app"].rsplit("-", 1)[0]  # Remove version

        analyses = list(ii.get_analyses(
            projects=project_pk,
            application__name=app_name,
            status="SUCCEEDED"
        ))

        expected = len(experiments) if step["type"] == "TARGETS" else len(pairs or [])
        actual = len(analyses)

        results[step["app"]] = {
            "expected": expected,
            "succeeded": actual,
            "complete": actual >= expected
        }

    print("\n=== Pipeline Verification ===")
    all_complete = True
    for app, stats in results.items():
        status = "✓" if stats["complete"] else "✗"
        print(f"  {status} {app}: {stats['succeeded']}/{stats['expected']}")
        if not stats["complete"]:
            all_complete = False

    return all_complete

verify_pipeline(PROJECT_PK, PIPELINE, experiments, pairs=pairs)
```

## Common Pipeline Patterns

### Run app only on samples missing results

```python
def run_missing(app_name, experiments, commit=False):
    """Run app only on experiments that don't have results yet."""
    # Get experiments with succeeded analyses
    succeeded = set()
    for a in ii.get_analyses(application__name=app_name, status="SUCCEEDED"):
        for t in a.targets:
            succeeded.add(t.pk)

    # Filter to missing
    missing = [e for e in experiments if e.pk not in succeeded]
    print(f"Running {app_name} on {len(missing)} missing experiments")

    if missing and commit:
        targets = ",".join(e.system_id for e in missing)
        os.system(f"isabl {app_name} --targets {targets} --commit")

run_missing("bwa-0-7-17", experiments, commit=False)
```

### Retry failed analyses

```python
def retry_failed(app_name, project_pk, commit=False):
    """Retry all failed analyses for an app."""
    failed = list(ii.get_analyses(
        projects=project_pk,
        application__name=app_name.rsplit("-", 1)[0],
        status="FAILED"
    ))

    print(f"Retrying {len(failed)} failed {app_name} analyses")

    for a in failed:
        targets = ",".join(t.system_id for t in a.targets)
        refs = ",".join(r.system_id for r in a.references)
        cmd = f"isabl {app_name} --targets {targets}"
        if refs:
            cmd += f" --references {refs}"
        cmd += " --force"
        if commit:
            cmd += " --commit"
        print(f"  {cmd}")
        if commit:
            os.system(cmd)

retry_failed("mutect-2-0", PROJECT_PK, commit=False)
```

### Chain dependent apps

```python
def run_after_dependency(app_name, dep_app_name, project_pk, commit=False):
    """Run app on experiments that have dependency completed."""
    # Get experiments with dependency succeeded
    ready = set()
    for a in ii.get_analyses(
        projects=project_pk,
        application__name=dep_app_name,
        status="SUCCEEDED"
    ):
        for t in a.targets:
            ready.add(t.system_id)

    # Get experiments already processed by this app
    done = set()
    for a in ii.get_analyses(
        projects=project_pk,
        application__name=app_name
    ):
        for t in a.targets:
            done.add(t.system_id)

    # Run on ready but not done
    to_run = ready - done
    print(f"Running {app_name} on {len(to_run)} experiments")

    if to_run and commit:
        targets = ",".join(to_run)
        os.system(f"isabl {app_name} --targets {targets} --commit")

# Run mutect after bwa completes
run_after_dependency("mutect-2-0", "bwa", PROJECT_PK)
```

### Generate pipeline script

```python
def generate_pipeline_script(experiments, pipeline, pairs=None, output="run_pipeline.sh"):
    """Generate a shell script to run the pipeline."""
    with open(output, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("set -e\n\n")
        f.write(f"# Pipeline generated on {datetime.now()}\n\n")

        for step in pipeline:
            app = step["app"]
            f.write(f"echo '=== Running {app} ==='\n")

            if step["type"] == "TARGETS":
                targets = ",".join(e.system_id for e in experiments)
                f.write(f"isabl {app} --targets {targets} --commit\n")
            elif step["type"] == "PAIRS" and pairs:
                for tumor, normal in pairs:
                    f.write(f"isabl {app} --targets {tumor.system_id} --references {normal.system_id} --commit\n")

            f.write("\n")

        f.write("echo 'Pipeline complete!'\n")

    print(f"Script written to {output}")
    print(f"Run with: bash {output}")

generate_pipeline_script(experiments, PIPELINE, pairs=pairs)
```
