# DNAKit execution: instrument-data-to-allotrope

Use scripts/convert_to_asm.py for supported instrument exports and scripts/validate_asm.py for schema checks. DNA scope includes NanoDrop concentration, DNA integrity and QuantStudio/qPCR exports. Install the upstream allotropy parser and format-specific optional dependencies. PDF/Excel support depends on the actual parser; preserving this example skill does not validate every instrument format.

Install this complete folder with `dnakit skills install instrument-data-to-allotrope --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor instrument-data-to-allotrope`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
