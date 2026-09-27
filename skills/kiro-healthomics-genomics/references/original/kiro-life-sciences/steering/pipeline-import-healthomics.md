---
inclusion: auto
fileMatchPattern: "**/*.{wdl,nf,cwl},**/nextflow.config"
---

# Import Pipeline to AWS HealthOmics

Step-by-step guide to import an nf-core or WDL pipeline into AWS HealthOmics.

## Step 1: Select a Pipeline

1. Browse available pipelines:
   - Use `life-sciences-pipelines` → `pipeline_search` tool to find pipelines by keyword.
   - Or browse by language: `pipeline_by_language` with "Nextflow", "WDL", or "CWL".
2. Get pipeline details including repository URL and HealthOmics import instructions.

## Step 2: Prepare the Workflow Definition

1. Clone or download the pipeline source from its repository.
2. For nf-core pipelines:
   - Ensure `nextflow.config` has HealthOmics-compatible settings.
   - Remove any local executor configurations.
3. For WDL pipelines:
   - Verify WDL version compatibility (HealthOmics supports WDL 1.0 and 1.1).
4. Lint the workflow:
   - Use `LintWorkflowDefinition` or `LintWorkflowBundle` for multi-file workflows.

## Step 3: Set Up Containers

1. Identify all container images used by the pipeline.
2. For each container, check ECR availability:
   - Use `CheckContainerAvailability` with the repository name and tag.
3. Clone missing containers to ECR:
   - Use `CloneContainerToECR` for each missing image.
4. Create a container registry map:
   - Use `CreateContainerRegistryMap` to generate mappings.

## Step 4: Package and Create Workflow

1. Package the workflow files:
   - Use `PackageWorkflow` with the main file and additional imports.
2. Create the HealthOmics workflow:
   - Use `CreateWorkflow` with the packaged definition, name, and container registry map.
3. Wait for workflow to become ACTIVE:
   - Use `GetWorkflow` to check status.

## Step 5: Configure and Run

1. Prepare input parameters matching the workflow's parameter template.
2. Start a test run with a small dataset:
   - Use `StartRun` with the workflow ID and test parameters.
3. Monitor the run:
   - Use `GetRun` and `ListRunTasks` to track progress.
4. If the run fails, diagnose:
   - Use `DiagnoseRunFailure` for detailed error analysis.

## Step 6: Optimize

1. After a successful run, analyze performance:
   - Use `AnalyzeRunPerformance` to identify optimization opportunities.
2. Adjust resource allocations based on recommendations.
3. Create a new workflow version with optimized settings:
   - Use `CreateWorkflowVersion` with updated definitions.

## Expected Outputs
- HealthOmics workflow ID ready for production runs.
- Container registry map for all pipeline containers.
- Performance baseline from test run.
