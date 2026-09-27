# Contributing tools and compatible fixes

DNAKit functions share one Python implementation with the MCP discovery and
execution endpoints. A new function does not need its own MCP server.

## Native conventions

- Keep public functions typed, documented and explicitly exported through
  `__all__`. Follow the existing domain module's input and result conventions.
- Accept keyword-only `verbose: bool = True`. Use `dnakit._progress.track` or
  `Progress` for long work; public entrypoints use DNAKit's verbosity wrapper.
- Prefer existing DNAKit functions. Use DNAKit exceptions for input,
  configuration and backend failures. Import optional backends inside functions.
- Document dependencies, file/network/process effects and Python/MCP usage.
- Generated exports live under `dnakit.generated`. Versioned implementation
  filenames do not change the public Python function or MCP tool name.

## Compatible repairs

Repairs fill missing behavior or correct a reproducible defect. Preserve the
function name, public import and MCP tool name. Keep all existing parameter
names, order, kinds, types and default values. Keep the return interface.

Only add optional parameters with defaults; new Python parameters must be
keyword-only so positional calls retain their meaning. Preserve old behavior
when the new option is omitted, apart from the documented defect being fixed.
Renaming/deleting parameters or adding required parameters is not a compatible
repair. Keep commands and scientific assumptions outside the identified defect
unchanged.

Before a repair, freeze examples showing both the defect and successful old
calls. Execute baseline and candidate implementations on the same examples.
Keep old successes, demonstrate the fix and add the examples as regression tests.
Interface checks and model review supplement these tests; they do not establish
correctness for every scientific input.

## Source and provenance

Use the existing DNAKit checkout for native code. External MCP implementations
need an identifiable source repository/revision or a caller-supplied source
directory. An MCP description does not provide the server's source code.
Keep upstream attribution and license files with adapted code; record the exact
revision and changes. Editing a local copy does not change a running remote MCP
service. Validate against a configured local service before claiming that a
server repair works.

## Reviewable submissions

Include the implementation, exports/metadata, reproducible tests and usage
documentation. Agent development can prepare these files and a contribution
record locally. Review fixtures for user data, credentials, machine-specific
paths and external-service dependencies before submitting them.

`.dnakit/development` holds local candidate code, logs, backups and execution
state. It is ignored by Git and is not a submission directory. Public manifests
must not contain local audit paths. Effectful generated tests require an explicit
environment opt-in as documented in their contribution record.

Install the development environment explicitly:

```bash
python -m pip install -e '.[dev]'
```

Check formatting, lint, types and relevant tests before submitting:

```bash
python -m ruff format --check src others/tests
python -m ruff check src others/tests
python -m mypy
python -m pytest others/tests/unit
```

Build an installable package when exports or resources change:

```bash
python -m build
```

Commit the selected contribution files on your branch, push to your fork and
open a pull request. Explain the original problem, the resulting behavior and
the validation actually performed. The maintainer decides whether to merge.
Agent development does not automatically commit, push, open a PR or merge.
