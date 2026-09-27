---
inclusion: auto
fileMatchPattern: "**/*.{wdl,nf,cwl}"
---

# Variant Calling Pipeline with HealthOmics

Step-by-step guide to create, configure, and run a variant calling workflow on AWS HealthOmics.

## Step 1: Prepare Workflow Definition

1. Choose a variant calling pipeline (e.g., GATK Best Practices germline short variant).
2. Lint the workflow definition to catch syntax errors:
   - Use `LintWorkflowDefinition` with `workflow_format: "wdl"` (or `"cwl"`).
3. Package the workflow files into a ZIP:
   - Use `PackageWorkflow` with the main workflow file and any imported sub-workflows.

## Step 2: Set Up Container Access

1. Check if required containers are available in ECR:
   - Use `CheckContainerAvailability` for each container image (e.g., `broadinstitute/gatk:4.5.0.0`).
2. If containers are not available, clone them to ECR:
   - Use `CloneContainerToECR` with the source image reference.
3. Grant HealthOmics access to ECR repositories:
   - Use `GrantHealthOmicsRepository` for each repository.
4. Create a container registry map:
   - Use `CreateContainerRegistryMap` to auto-discover pull-through caches.

## Step 3: Create HealthOmics Workflow

1. Create the workflow:
   - Use `CreateWorkflow` with the packaged ZIP, name, and parameter template.
2. Verify the workflow was created successfully:
   - Use `GetWorkflow` with the returned workflow ID.
   - Wait for status to become `ACTIVE`.

## Step 4: Prepare Input Data

1. Locate input FASTQ or BAM files in S3.
2. Locate the reference genome (e.g., GRCh38) in S3 or a reference store.
3. Prepare the parameters JSON matching the workflow's parameter template.

## Step 5: Start the Run

1. Start the workflow run:
   - Use `StartRun` with workflow_id, role_arn, output_uri, and parameters.
   - Use `storage_type: "DYNAMIC"` unless inputs exceed 1 TiB.
2. Note the returned run ID for monitoring.

## Step 6: Monitor the Run

1. Check run status:
   - Use `GetRun` with the run ID. Look for status: `RUNNING`, `COMPLETED`, or `FAILED`.
2. List task-level progress:
   - Use `ListRunTasks` to see individual task statuses.
3. If a run fails:
   - Use `DiagnoseRunFailure` for actionable error messages and recommendations.

## Step 7: Analyze Results

1. Check run outputs in the S3 output URI.
2. Generate a timeline visualization:
   - Use `GenerateRunTimeline` to see task parallelism and identify bottlenecks.
3. Analyze performance for optimization:
   - Use `AnalyzeRunPerformance` to get resource utilization and cost recommendations.

## Expected Outputs
- VCF file with called variants.
- BAM file with aligned and processed reads.
- Run metrics and quality reports.
