---
name: isabl-monitor-analyses
description: Monitor and track the status of Isabl analyses. Use when checking job progress, finding failed analyses, or tracking pipeline execution.
tools: Read, Bash
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-monitor-analyses` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Monitoring Isabl Analyses

You are helping the user monitor and track the status of analyses.

## Checklist

Work through these steps systematically:

1. [ ] **Define scope** (project, application, time range)
2. [ ] **Check status counts** (overall health)
3. [ ] **Find problematic analyses** (failed, stale)
4. [ ] **Investigate failures** (logs, errors)
5. [ ] **Track running jobs** (progress, runtime)
6. [ ] **Generate status report**
7. [ ] **Suggest actions** (retry, debug, escalate)

## Step 1: Define Scope

```python
import isabl_cli as ii
from datetime import datetime, timedelta

# Define what to monitor
PROJECT_PK = 102
APPLICATION_NAME = "MUTECT"  # or None for all
HOURS_BACK = 24  # Look at last N hours

cutoff = datetime.now() - timedelta(hours=HOURS_BACK)
```

## Step 2: Check Status Counts

```python
from collections import Counter

# Get all analyses in scope
filters = {"projects": PROJECT_PK}
if APPLICATION_NAME:
    filters["application__name"] = APPLICATION_NAME

analyses = list(ii.get_analyses(**filters))

# Count by status
status_counts = Counter(a.status for a in analyses)

print(f"=== Analysis Status (Project {PROJECT_PK}) ===")
print(f"Total: {len(analyses)}")
for status in ["SUCCEEDED", "FAILED", "STARTED", "SUBMITTED", "STAGED", "CREATED"]:
    count = status_counts.get(status, 0)
    pct = (count / len(analyses) * 100) if analyses else 0
    print(f"  {status}: {count} ({pct:.1f}%)")
```

## Step 3: Find Problematic Analyses

### Failed Analyses

```python
failed = ii.get_analyses(
    projects=PROJECT_PK,
    status="FAILED",
    modified__gte=cutoff.isoformat()
)

if failed:
    print(f"\n=== FAILED ({len(failed)}) ===")
    for a in failed[:10]:
        target = a.targets[0].system_id if a.targets else "N/A"
        print(f"  [{a.pk}] {a.application.name}: {target}")
        print(f"       Failed: {a.finished}")
```

### Stale Analyses (Running Too Long)

```python
stale_cutoff = datetime.now() - timedelta(hours=48)

stale = ii.get_analyses(
    projects=PROJECT_PK,
    status="STARTED",
    started__lt=stale_cutoff.isoformat()
)

if stale:
    print(f"\n=== STALE (running > 48h) ({len(stale)}) ===")
    for a in stale:
        hours = (datetime.now() - a.started).total_seconds() / 3600
        print(f"  [{a.pk}] {a.application.name}: running {hours:.1f}h")
```

### Stuck in Queue

```python
queue_cutoff = datetime.now() - timedelta(hours=6)

stuck = ii.get_analyses(
    projects=PROJECT_PK,
    status="SUBMITTED",
    submitted__lt=queue_cutoff.isoformat()
)

if stuck:
    print(f"\n=== STUCK IN QUEUE (> 6h) ({len(stuck)}) ===")
    for a in stuck:
        hours = (datetime.now() - a.submitted).total_seconds() / 3600
        print(f"  [{a.pk}] {a.application.name}: waiting {hours:.1f}h")
```

## Step 4: Investigate Failures

```python
from pathlib import Path

def investigate_failure(analysis):
    """Get failure details for an analysis."""
    info = {
        "pk": analysis.pk,
        "application": f"{analysis.application.name} v{analysis.application.version}",
        "target": analysis.targets[0].system_id if analysis.targets else None,
        "storage_url": analysis.storage_url,
        "run_time": f"{analysis.run_time:.1f} min" if analysis.run_time else "N/A",
    }

    # Try to read last lines of error log
    if analysis.storage_url:
        err_path = Path(analysis.storage_url) / "head_job.err"
        if err_path.exists():
            with open(err_path) as f:
                lines = f.readlines()
                info["last_error"] = "".join(lines[-10:])

    return info

# Investigate recent failures
for a in failed[:5]:
    details = investigate_failure(a)
    print(f"\n--- Analysis {details['pk']} ---")
    print(f"App: {details['application']}")
    print(f"Target: {details['target']}")
    print(f"Runtime: {details['run_time']}")
    if details.get("last_error"):
        print(f"Last error:\n{details['last_error']}")
```

## Step 5: Track Running Jobs

```python
running = ii.get_analyses(
    projects=PROJECT_PK,
    status="STARTED",
    ordering="-started"
)

if running:
    print(f"\n=== RUNNING ({len(running)}) ===")
    for a in running[:20]:
        hours = (datetime.now() - a.started).total_seconds() / 3600
        target = a.targets[0].system_id if a.targets else "N/A"
        print(f"  [{a.pk}] {a.application.name}: {target} ({hours:.1f}h)")
```

### Check Job Progress (if supported by app)

```python
from pathlib import Path
import json

def check_progress(analysis):
    """Check progress file if application writes one."""
    progress_file = Path(analysis.storage_url) / "progress.json"
    if progress_file.exists():
        with open(progress_file) as f:
            return json.load(f)
    return None

for a in running[:5]:
    progress = check_progress(a)
    if progress:
        print(f"  [{a.pk}] Progress: {progress.get('percent', 'unknown')}%")
```

## Step 6: Generate Status Report

```python
from datetime import datetime

def generate_status_report(project_pk, app_name=None):
    """Generate a status report for a project."""
    filters = {"projects": project_pk}
    if app_name:
        filters["application__name"] = app_name

    analyses = list(ii.get_analyses(**filters))
    status_counts = Counter(a.status for a in analyses)

    # Calculate metrics
    total = len(analyses)
    succeeded = status_counts.get("SUCCEEDED", 0)
    failed = status_counts.get("FAILED", 0)
    running = status_counts.get("STARTED", 0)
    queued = status_counts.get("SUBMITTED", 0) + status_counts.get("STAGED", 0)

    # Success rate
    completed = succeeded + failed
    success_rate = (succeeded / completed * 100) if completed else 0

    # Average runtime
    runtimes = [a.run_time for a in analyses if a.run_time and a.status == "SUCCEEDED"]
    avg_runtime = sum(runtimes) / len(runtimes) if runtimes else 0

    report = f"""
# Analysis Status Report
**Generated**: {datetime.now().strftime("%Y-%m-%d %H:%M")}
**Project**: {project_pk}
**Application**: {app_name or "All"}

## Summary
| Metric | Value |
|--------|-------|
| Total | {total} |
| Succeeded | {succeeded} ({success_rate:.1f}%) |
| Failed | {failed} |
| Running | {running} |
| Queued | {queued} |
| Avg Runtime | {avg_runtime:.1f} min |

## Status Breakdown
"""
    for status, count in sorted(status_counts.items()):
        pct = count / total * 100 if total else 0
        report += f"- {status}: {count} ({pct:.1f}%)\n"

    return report

print(generate_status_report(PROJECT_PK, APPLICATION_NAME))
```

## Step 7: Suggest Actions

```python
def suggest_actions(project_pk):
    """Suggest actions based on analysis status."""
    actions = []

    # Check for failures
    failed = list(ii.get_analyses(projects=project_pk, status="FAILED"))
    if failed:
        actions.append(f"- Debug {len(failed)} failed analyses (use /isabl-debug-analysis)")

    # Check for stale
    stale_cutoff = datetime.now() - timedelta(days=2)
    stale = list(ii.get_analyses(
        projects=project_pk,
        status="STARTED",
        started__lt=stale_cutoff.isoformat()
    ))
    if stale:
        actions.append(f"- Investigate {len(stale)} stale analyses (running > 48h)")

    # Check for queue issues
    queue_cutoff = datetime.now() - timedelta(hours=12)
    stuck = list(ii.get_analyses(
        projects=project_pk,
        status="SUBMITTED",
        submitted__lt=queue_cutoff.isoformat()
    ))
    if stuck:
        actions.append(f"- Check batch system for {len(stuck)} stuck jobs")

    if not actions:
        actions.append("- No immediate actions required")

    print("\n## Suggested Actions")
    for action in actions:
        print(action)

suggest_actions(PROJECT_PK)
```

## Common Monitoring Queries

### Analyses by application today

```python
today = datetime.now().date().isoformat()
analyses = ii.get_analyses(
    projects=PROJECT_PK,
    created__date=today
)
apps = Counter(a.application.name for a in analyses)
for app, count in apps.most_common():
    print(f"{app}: {count}")
```

### Failed analyses with specific error

```python
import os

error_pattern = "OutOfMemory"
matches = []

for a in ii.get_analyses(projects=PROJECT_PK, status="FAILED"):
    err_file = os.path.join(a.storage_url, "head_job.err")
    if os.path.exists(err_file):
        with open(err_file) as f:
            if error_pattern in f.read():
                matches.append(a)

print(f"Found {len(matches)} analyses with '{error_pattern}' error")
```

### Track a specific analysis

```python
def watch_analysis(pk, interval=60):
    """Watch an analysis until completion."""
    import time

    while True:
        a = ii.get_instance("analyses", pk)
        print(f"[{datetime.now()}] Analysis {pk}: {a.status}")

        if a.status in ["SUCCEEDED", "FAILED"]:
            print(f"Completed with status: {a.status}")
            break

        time.sleep(interval)

# watch_analysis(12345)
```

### Bulk retry failed analyses

```python
# Get failed analyses to retry
failed = list(ii.get_analyses(
    projects=PROJECT_PK,
    application__name="MUTECT",
    status="FAILED"
))

# Generate retry commands
print("# Retry commands:")
for a in failed:
    target = a.targets[0].system_id if a.targets else ""
    print(f"isabl mutect-2-0 --targets {target} --force --commit")
```
