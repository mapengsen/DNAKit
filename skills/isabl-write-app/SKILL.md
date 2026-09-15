---
name: isabl-write-app
description: Guide through creating a new Isabl bioinformatics application. Use when building pipelines that integrate with the Isabl platform.
tools: Read, Write, Edit, Glob, Grep
model: inherit
---

## DNAKit execution

Use `dnakit skills doctor isabl-write-app` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Writing an Isabl Application

You are guiding the user through creating an Isabl bioinformatics application.

## Checklist

Work through these steps systematically:

1. [ ] **Define application metadata** (NAME, VERSION, ASSEMBLY, SPECIES)
2. [ ] **Choose CLI options** (TARGETS, REFERENCES, PAIRS)
3. [ ] **Define application_settings** for configurable paths/parameters
4. [ ] **Implement validate_experiments()** to check input validity
5. [ ] **Implement get_dependencies()** if app depends on other app results
6. [ ] **Implement get_command()** to generate the shell command
7. [ ] **Implement get_analysis_results()** to extract outputs
8. [ ] **Register in INSTALLED_APPLICATIONS**
9. [ ] **Write tests** using pytest fixtures

## Step 1: Define Application Metadata

Create the application class with required metadata:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    """
    Brief description of what this application does.
    """

    # Required metadata
    NAME = "my_application"
    VERSION = "1.0.0"

    # Optional: restrict to specific assembly/species
    ASSEMBLY = "GRCh37"  # or "GRCh38", None for any
    SPECIES = "HUMAN"    # or None for any
```

## Step 2: Choose CLI Options

Configure how the application receives input experiments:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"

    # CLI configuration
    cli_help = "Run my application on experiments"
    cli_options = [options.TARGETS]  # or REFERENCES, PAIRS
```

### CLI Option Types

- `options.TARGETS` - Single or multiple target experiments
- `options.REFERENCES` - Reference experiments (e.g., normals)
- `options.PAIRS` - Tumor-normal pairs

## Step 3: Define Application Settings

Configure settings that can be overridden in the database:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    # Configurable settings (can be overridden in database)
    application_settings = {
        "tool_path": "/usr/bin/mytool",
        "threads": 4,
        "memory_gb": 16,
    }

    # Define expected results
    application_results = {
        "output_file": {
            "frontend_type": "text-file",
            "description": "Main output file",
            "verbose_name": "Output",
        }
    }
```

## Step 4: Implement validate_experiments()

Validate input experiments before creating analyses:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    def validate_experiments(self, targets, references):
        """
        Raise AssertionError if experiments are invalid for this app.
        Called before creating analyses.
        """
        assert len(targets) == 1, "Requires exactly one target experiment"
        assert targets[0].technique.method == "WGS", "Only WGS supported"
```

### Tumor-Normal Pair Validation

```python
from isabl_cli import AbstractApplication, options

class MyPairApplication(AbstractApplication):
    NAME = "my_pair_application"
    VERSION = "1.0.0"
    cli_options = [options.PAIRS]

    def validate_experiments(self, targets, references):
        assert len(targets) == 1, "One tumor per analysis"
        assert len(references) == 1, "One normal per analysis"
        assert targets[0].sample.category == "TUMOR"
        assert references[0].sample.category == "NORMAL"
```

## Step 5: Implement get_dependencies()

Define dependencies on other application results:

```python
from isabl_cli import AbstractApplication, options, utils

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    application_settings = {
        "alignment_app_pk": 10,  # PK of alignment application
    }

    def get_dependencies(self, targets, references, settings):
        """
        Return (dependency_analyses, inputs_dict) if this app needs
        results from other applications.
        """
        bam, analysis_key = utils.get_result(
            experiment=targets[0],
            application_key=settings.alignment_app_pk,
            result_key="bam"
        )
        return [analysis_key], {"input_bam": bam}
```

## Step 6: Implement get_command()

Generate the shell command to execute:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    application_settings = {
        "tool_path": "/usr/bin/mytool",
        "threads": 4,
    }

    def get_command(self, analysis, inputs, settings):
        """
        Return the shell command to execute.
        This is the core of the application.
        """
        target = analysis.targets[0]
        output_dir = analysis.storage_url

        return f"""
        {settings.tool_path} \\
            --input {target.bam_files["GRCh37"]["url"]} \\
            --output {output_dir}/result.txt \\
            --threads {settings.threads}
        """
```

## Step 7: Implement get_analysis_results()

Extract output paths after successful completion:

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    NAME = "my_application"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    application_results = {
        "output_file": {
            "frontend_type": "text-file",
            "description": "Main output file",
            "verbose_name": "Output",
        }
    }

    def get_analysis_results(self, analysis):
        """
        Return dict of result paths after successful completion.
        Keys should match application_results.
        """
        return {
            "output_file": f"{analysis.storage_url}/result.txt"
        }
```

## Step 8: Register in INSTALLED_APPLICATIONS

Add the application to the appropriate assembly apps module (e.g., `apps_grch38`):

```python
INSTALLED_APPLICATIONS = [
    "my_apps.MyApplication",
    "my_apps.MyPairApplication",
]
```

### Running the Application

Applications are invoked via the assembly-specific CLI command:

```bash
# Dry run (no --commit)
isabl apps-<assembly-name> <app-name>-<version> --targets EXPERIMENT_SYSTEM_ID

# Paired tumor-normal
isabl apps-<assembly-name> <app-name>-<version> --pairs TUMOR_SYSTEM_ID NORMAL_SYSTEM_ID

# Actually submit
isabl apps-<assembly-name> <app-name>-<version> --pairs TUMOR_SYSTEM_ID NORMAL_SYSTEM_ID --commit
```

For example, a GRCh38 app called `my_application` version `1.0.0`:

```bash
isabl apps-grch38 my-application-1.0.0 --pairs IID_H209923_T01_01_WG01 IID_H209923_N01_01_WG01 --commit
```

## Step 9: Write Tests

Test the application using pytest fixtures:

```python
def test_my_application(tmpdir, commit):
    import isabl_cli as ii
    from isabl_cli import api, factories
    from my_apps import MyApplication

    # Create test experiment
    experiment = api.create_instance(
        "experiments",
        **factories.ExperimentFactory()
    )

    # Run application
    app = MyApplication()
    app.run(
        tuples=[([experiment], [])],
        commit=commit
    )
```

## Complete Application Template

```python
from isabl_cli import AbstractApplication, options

class MyApplication(AbstractApplication):
    """
    Brief description of what this application does.
    """

    # Required metadata
    NAME = "my_application"
    VERSION = "1.0.0"

    # Optional: restrict to specific assembly/species
    ASSEMBLY = "GRCh37"  # or "GRCh38", None for any
    SPECIES = "HUMAN"    # or None for any

    # CLI configuration
    cli_help = "Run my application on experiments"
    cli_options = [options.TARGETS]  # or REFERENCES, PAIRS

    # Configurable settings (can be overridden in database)
    application_settings = {
        "tool_path": "/usr/bin/mytool",
        "threads": 4,
    }

    # Define expected results
    application_results = {
        "output_file": {
            "frontend_type": "text-file",
            "description": "Main output file",
            "verbose_name": "Output",
        }
    }

    def validate_experiments(self, targets, references):
        """
        Raise AssertionError if experiments are invalid for this app.
        Called before creating analyses.
        """
        assert len(targets) == 1, "Requires exactly one target experiment"
        assert targets[0].technique.method == "WGS", "Only WGS supported"

    def get_dependencies(self, targets, references, settings):
        """
        Return (dependency_analyses, inputs_dict) if this app needs
        results from other applications.
        """
        return [], {}

    def get_command(self, analysis, inputs, settings):
        """
        Return the shell command to execute.
        This is the core of the application.
        """
        target = analysis.targets[0]
        output_dir = analysis.storage_url

        return f"""
        {settings.tool_path} \\
            --input {target.bam_files["GRCh37"]["url"]} \\
            --output {output_dir}/result.txt \\
            --threads {settings.threads}
        """

    def get_analysis_results(self, analysis):
        """
        Return dict of result paths after successful completion.
        Keys should match application_results.
        """
        return {
            "output_file": f"{analysis.storage_url}/result.txt"
        }
```

## Advanced Patterns

### Multiple Targets (Cohort Analysis)

```python
from isabl_cli import AbstractApplication, options, api

class CohortApplication(AbstractApplication):
    NAME = "cohort_analysis"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]
    unique_analysis_per_individual = False  # Allow multiple targets

    def get_experiments_from_cli_options(self, **cli_options):
        # Custom logic to group experiments
        targets = api.get_instances("experiments", **filters)
        return [(targets, [])]  # Single analysis for all targets
```

### Project-Level Merge

```python
from isabl_cli import AbstractApplication, options

class MergeApplication(AbstractApplication):
    NAME = "merge_results"
    VERSION = "1.0.0"
    cli_options = [options.TARGETS]

    application_project_level_results = {
        "merged_output": {
            "frontend_type": "text-file",
            "description": "Merged output across project",
        }
    }

    def merge_project_analyses(self, analysis, analyses):
        # Combine results from all analyses in project
        pass

    def get_project_analysis_results(self, analysis):
        return {"merged_output": f"{analysis.storage_url}/merged.txt"}
```

## Common Issues

- **AssertionError in validate_experiments**: Check your validation logic matches actual data
- **Command fails silently**: Check `head_job.log` and `head_job.err` in storage_url
- **Dependencies not found**: Ensure dependency app has SUCCEEDED status
- **Results not extracted**: Verify paths in get_analysis_results() exist
