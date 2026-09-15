---
name: pod5
description: Inspect, merge, filter, subset, repack, and convert nanopore POD5
  files using the pod5 CLI. Use when working with .pod5 files, converting fast5
  to pod5 or pod5 to fast5, extracting read summaries, filtering reads by
  read_id, subsetting by barcode or other metadata columns, merging multiple pod5
  files, or updating pod5 files to the latest schema. Triggers on tasks involving
  .pod5, .fast5 conversion, nanopore raw signal files, read_id filtering, basecall
  input preparation, or Oxford Nanopore sequencing data.
---

## DNAKit execution

Use `dnakit skills doctor pod5` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


<!-- DNAKit metadata portability edit, 2026-09-12; see UPSTREAM.md. -->

# pod5

## Environment

Install the POD5 CLI separately using the upstream `environment.yaml` or the official POD5 installation instructions. The `pod5` executable must be available on PATH.

Before issuing any commands, resolve the full absolute path for this machine:
```bash
command -v pod5
```
If the command prints no path, install and activate the POD5 environment before continuing. Use the printed output literally as the first token in every command.
In examples below, `$POD5` is a readable placeholder for that resolved path.

## Subcommands

**Inspection**
- `view` — write a tab-separated summary table of all reads (replaces `inspect reads`; 200× faster)
- `inspect reads` — (deprecated) print csv table of all reads
- `inspect read` — print detailed info for a single read by read_id
- `inspect summary` — print a high-level summary of pod5 file contents
- `inspect debug` — print debugging information

**Manipulation**
- `merge` — merge multiple pod5 files into one; checks for duplicate read_ids by default
- `filter` — extract reads by read_id list into a single output file
- `subset` — subset reads into one or more outputs using a CSV mapping or table columns
- `repack` — repack pod5 files one-to-one (useful for defragmentation)
- `recover` — attempt to recover corrupted pod5 files (writes `_recovered.pod5`)
- `update` — update old pod5 files to the latest schema version

**Conversion**
- `convert fast5` — convert fast5 file(s) to pod5
- `convert to_fast5` — convert pod5 file(s) to fast5 (4000 reads/file by default)

## Common patterns

**Summarise reads in a pod5 file:**
```bash
$POD5 view input.pod5
# Only read_ids, no header (pipe-friendly):
$POD5 view input.pod5 --ids --no-header
# Select specific fields:
$POD5 view input.pod5 --include "read_id, channel, num_samples, end_reason"
# Write to file:
$POD5 view *.pod5 --output summary.tsv
```

**Merge pod5 files:**
```bash
$POD5 merge *.pod5 --output merged.pod5
# Allow duplicate read_ids (e.g. after re-basecalling):
$POD5 merge *.pod5 --output merged.pod5 --duplicate-ok
```

**Filter reads by read_id list:**
```bash
# read_ids.txt: one UUID per line, comments (#) allowed
$POD5 filter input.pod5 --ids read_ids.txt --output filtered.pod5
# Allow some requested IDs to be absent:
$POD5 filter input.pod5 --ids read_ids.txt --output filtered.pod5 --missing-ok
```

**Random subsample of 1000 reads:**
```bash
$POD5 view all.pod5 --ids --no-header --output all_ids.txt
sort --random-sort all_ids.txt | head -1000 > 1k_ids.txt
$POD5 filter all.pod5 --ids 1k_ids.txt --output 1k.pod5
```

**Subset by metadata column (e.g. barcode):**
```bash
# Generate summary table first, then subset on barcode column
$POD5 view input.pod5 --output summary.tsv
$POD5 subset input.pod5 --output barcode_subset/ --summary summary.tsv --columns barcode
# Custom output filename template:
$POD5 subset input.pod5 --output barcode_subset/ --summary summary.tsv --columns barcode \
    --template "{barcode}.pod5"
```

**Convert fast5 → pod5:**
```bash
# All fast5 into one pod5:
$POD5 convert fast5 ./input/*.fast5 --output converted.pod5
# One-to-one (preserves directory structure):
$POD5 convert fast5 ./input/*.fast5 --output output_pod5s/ --one-to-one ./input/
```

**Convert pod5 → fast5:**
```bash
$POD5 convert to_fast5 input.pod5 --output pod5_to_fast5/
# Control reads per file (default 4000):
$POD5 convert to_fast5 input.pod5 --output pod5_to_fast5/ --file-read-count 2000
```

**Update old pod5 files to latest schema:**
```bash
$POD5 update old/*.pod5 --output updated/
```

## Allowlist entries

Resolve and add to your terminal command allowlist (Cursor: Settings → Features → Terminal):
```bash
command -v pod5
```

## Full flag reference

To look up all flags for a specific subcommand:
```bash
grep -A 50 "^### \`subcommand\`" "$(dirname <path-to-this-SKILL.md>)/reference.md"
```
Full reference: [reference.md](reference.md)

## Patterns

Reusable real-world patterns accumulated over time. To search:
```bash
grep -A 20 "keyword" "$(dirname <path-to-this-SKILL.md>)/patterns.md"
```
[patterns.md](patterns.md)
