---
name: isabl-debug-analysis
description: Systematically debug a failed Isabl analysis. Use when an analysis has FAILED status or unexpected behavior.
tools: Read, Bash, Glob, Grep
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-debug-analysis` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Debugging an Isabl Analysis

You are helping debug a failed or problematic Isabl analysis.

## Checklist

Work through these steps systematically:

1. [ ] **Get analysis details** (status, application, targets)
2. [ ] **Check the command script** (head_job.sh)
3. [ ] **Review stdout log** (head_job.log)
4. [ ] **Review stderr log** (head_job.err)
5. [ ] **Check input data exists** (BAM files, dependencies)
6. [ ] **Verify settings** (application_settings values)
7. [ ] **Identify root cause**
8. [ ] **Suggest fix**

## Step 1: Get Analysis Details

```python
import isabl_cli as ii

# Get the analysis
analysis = ii.Analysis(ANALYSIS_PK)

print(f"Status: {analysis.status}")
print(f"Application: {analysis.application.name} v{analysis.application.version}")
print(f"Storage URL: {analysis.storage_url}")
print(f"Targets: {[t.system_id for t in analysis.targets]}")
print(f"References: {[r.system_id for r in analysis.references]}")
```

## Step 2: Check Command Script

The generated command is in `head_job.sh`:

```bash
# View the command that was executed
cat {analysis.storage_url}/head_job.sh
```

Look for:
- Incorrect paths
- Missing environment variables
- Malformed command syntax

## Step 3: Review Logs

```bash
# Standard output
cat {analysis.storage_url}/head_job.log

# Standard error (usually more informative)
cat {analysis.storage_url}/head_job.err
```

Common error patterns:

| Error | Likely Cause |
|-------|--------------|
| "File not found" | Input file path incorrect or missing |
| "Permission denied" | Storage directory permissions |
| "Command not found" | Tool path in settings is wrong |
| "Memory allocation" | Job needs more memory (batch system) |
| "Timeout" | Job exceeded time limit |

## Step 4: Check Input Data

```python
# Verify target experiments have data
for target in analysis.targets:
    print(f"{target.system_id}:")
    print(f"  Raw data: {target.raw_data}")
    print(f"  BAM files: {target.bam_files}")
```

```bash
# Check files exist
ls -la {path_to_input_file}
```

## Step 5: Check Dependencies

If the app depends on other analyses:

```python
# Find dependency analyses
from isabl_cli import utils

# Check if dependency result exists
result, dep_analysis = utils.get_result(
    experiment=analysis.targets[0],
    application_key=DEPENDENCY_APP_PK,
    result_key="expected_result"
)
print(f"Dependency analysis: {dep_analysis}")
print(f"Dependency result: {result}")
```

## Step 6: Verify Settings

```python
# Get application settings
app = ii.get_instance("applications", analysis.application.pk)
print(f"Settings: {app.settings}")

# Check specific setting
from my_apps import MyApplication
instance = MyApplication()
print(f"Tool path: {instance.settings.tool_path}")
```

## Step 7: Common Fixes

### Fix 1: Re-run with --force

```bash
# Wipe failed analysis and restart
isabl myapp-1-0 --analyses pk={ANALYSIS_PK} --force --commit
```

### Fix 2: Re-run with --restart

```bash
# Resume from checkpoint (if supported)
isabl myapp-1-0 --analyses pk={ANALYSIS_PK} --restart --commit
```

### Fix 3: Run locally for debugging

```bash
# Execute locally instead of batch system
isabl myapp-1-0 --targets EXPERIMENT_ID --local --commit
```

### Fix 4: Manually patch status

```python
# If analysis completed but status is wrong
ii.patch_instance("analyses", analysis.pk, status="SUCCEEDED")
```

### Fix 5: Check batch system logs

```bash
# For SLURM
sacct -j JOB_ID --format=JobID,State,ExitCode,MaxRSS,Elapsed

# For LSF
bjobs -l JOB_ID
bhist -l JOB_ID
```

## Signals and Automations

If the analysis SUCCEEDED but something else failed:

```python
# Check if signals are configured
from isabl_cli.settings import system_settings
print(f"ON_STATUS_CHANGE: {system_settings.ON_STATUS_CHANGE}")

# Manually trigger signals
from isabl_cli import signals
signals.run_on_status_change(analysis)
```

## Reporting the Issue

When reporting, include:
1. Analysis PK and status
2. Application name and version
3. Relevant portion of head_job.err
4. Input experiment IDs
5. Steps to reproduce
