---
inclusion: manual
---

# Genomics Pipeline Best Practices

Practical guidance for designing, optimizing, and running genomics workflows in WDL, Nextflow, and CWL.

## Workflow Design Patterns

### Modular Task Design
- Break workflows into small, reusable tasks with clear inputs and outputs.
- Each task should do one thing well (e.g., alignment, sorting, variant calling).
- Use sub-workflows for common patterns like "align + sort + markdup".

### Input Validation
- Validate file formats and required inputs at workflow start before expensive compute.
- Check reference genome compatibility with input reads (species, assembly version).
- Validate sample sheet structure and required columns early.

```wdl
task ValidateInputs {
  input { File fastq_r1; File fastq_r2; File reference }
  command <<<
    samtools quickcheck ~{fastq_r1} ~{fastq_r2}
    samtools faidx ~{reference}
  >>>
}
```

## Resource Optimization

### CPU and Memory
- Profile tasks to right-size CPU and memory requests. Over-provisioning wastes cost.
- Alignment tasks (BWA, STAR): 8–16 CPUs, 32–64 GB RAM.
- Variant calling (GATK HaplotypeCaller): 2–4 CPUs, 8–16 GB RAM.
- Sorting and indexing: 4–8 CPUs, 16–32 GB RAM.
- Use HealthOmics `AnalyzeRunPerformance` to identify over-provisioned tasks.

### Storage
- Use DYNAMIC storage on HealthOmics unless inputs exceed 1 TiB.
- Compress intermediate files (BAM → CRAM, VCF → VCF.gz + tabix).
- Clean up temporary files within tasks to reduce storage footprint.

### Parallelism
- Scatter by chromosome or genomic interval for variant calling.
- Scatter by sample for multi-sample workflows.
- Use `maxRetries` for transient failures in cloud environments.

## Error Handling

### Retry Strategy
- Set `maxRetries: 2` for tasks that may fail due to transient cloud issues.
- Use preemptible/spot instances for fault-tolerant tasks with retry enabled.
- Log stderr to a file for post-mortem debugging.

### Checkpointing
- Design workflows so completed tasks are not re-run on restart.
- Use HealthOmics run caching (`CACHE_ON_FAILURE`) to avoid re-running successful tasks.
- Store intermediate outputs in persistent storage (S3) for long-running pipelines.

```nextflow
process alignment {
  cache 'lenient'
  errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'finish' }
  maxRetries 2
  // ...
}
```

## Reproducibility

### Version Pinning
- Pin all tool versions in container images (e.g., `bwa:0.7.17`, not `bwa:latest`).
- Use container digests (`sha256:...`) for maximum reproducibility.
- Record the workflow version, input parameters, and reference genome version in run metadata.

### Container Best Practices
- Use ECR pull-through caches for HealthOmics workflows.
- Grant HealthOmics access to ECR repositories with `GrantHealthOmicsRepository`.
- Use multi-stage Docker builds to minimize image size.

### Workflow Testing
- Test with a small downsampled dataset before running on full data.
- Validate outputs against known truth sets (e.g., Genome in a Bottle).
- Use `LintWorkflowDefinition` to catch syntax errors before deployment.

## HealthOmics-Specific Tips

- Use `CreateWorkflow` to register workflows, then `StartRun` to execute.
- Monitor runs with `GetRun` and `ListRunTasks` for task-level status.
- Use `DiagnoseRunFailure` for failed runs to get actionable error messages.
- Use `GenerateRunTimeline` to visualize task parallelism and identify bottlenecks.
- Set up run groups with `CreateRunGroup` to limit concurrent resource usage.

## Common Pipeline Patterns

| Pattern | Use Case | Key Tools |
|---------|----------|-----------|
| Germline variant calling | WGS/WES SNP/indel detection | BWA + GATK HaplotypeCaller |
| Somatic variant calling | Tumor-normal pair analysis | BWA + Mutect2 |
| RNA-seq quantification | Gene expression analysis | STAR + featureCounts |
| ChIP-seq peak calling | Transcription factor binding | BWA + MACS2 |
| Methylation analysis | Bisulfite sequencing | Bismark + methylKit |
