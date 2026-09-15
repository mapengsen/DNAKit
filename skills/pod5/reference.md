# pod5 — Full Reference

Binary: `~/.cursor/skills/pod5/bin/bin/pod5`

Each entry contains the verbatim `--help` output. Grep for a subcommand:
```bash
grep -A 50 "^### \`subcommand\`" ~/.cursor/skills/pod5/reference.md
```
Increase `-A` if output appears truncated.

---

## Top-level

### `pod5`

```
usage: pod5 [-h] [-v]
            {convert,inspect,merge,repack,subset,filter,recover,update,view}
            ...

**********      POD5 Tools      **********

Tools for inspecting, converting, subsetting and formatting POD5 files

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show pod5 version and exit.

Example: pod5 convert fast5 input.fast5 --output output.pod5
```

---

## Inspection

### `view`

```
usage: pod5 view [-h] [-o OUTPUT] [-r] [-f] [-t THREADS] [-H]
                 [--separator SEPARATOR] [-I] [-i INCLUDE] [-x EXCLUDE] [-L]
                 [inputs [inputs ...]]

    Write contents of some pod5 file(s) as a table to stdout or --output if given.
    The default separator is <tab>.
    The column order is always as shown in -L/--list-fields"
    

positional arguments:
  inputs                Input pod5 file(s) to view (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -f, --force-overwrite
                        Overwrite destination files (default: False)
  -t THREADS, --threads THREADS
                        Set the number of reader workers (default: 4)

Formatting:
  -H, --no-header       Omit the header line (default: False)
  --separator SEPARATOR
                        Table separator character (e.g. ',') (default: )

Selection:
  -I, --ids             Only write 'read_id' field (default: False)
  -i INCLUDE, --include INCLUDE
                        Include a double-quoted comma-separated list of fields
                        (default: None)
  -x EXCLUDE, --exclude EXCLUDE
                        Exclude a double-quoted comma-separated list of
                        fields. (default: None)

List Fields:
  -L, --list-fields     List all groups and fields available for selection and
                        exit (default: False)

Example: pod5 view input.pod5
```

### `inspect`

```
usage: pod5 inspect [-h] {summary,reads,read,debug} ...

Inspect the contents of a pod5 file

optional arguments:
  -h, --help            show this help message and exit

command:
  {summary,reads,read,debug}

Example: pod5 inspect reads input.pod5
```

### `inspect reads`

```
usage: pod5 inspect reads [-h] [-r] input_files [input_files ...]

Print read information on all reads as a csv table

positional arguments:
  input_files

optional arguments:
  -h, --help       show this help message and exit
  -r, --recursive  Search for input files recursively matching `*.pod5`

Example: pod5 inspect reads input.pod5
```

### `inspect read`

```
usage: pod5 inspect read [-h] input_files read_id

Print detailed read information for a named read id

positional arguments:
  input_files
  read_id

optional arguments:
  -h, --help   show this help message and exit

Example: pod5 inspect read input.pod5 0000173c-bf67-44e7-9a9c-1ad0bc728e74
```

### `inspect summary`

```
usage: pod5 inspect summary [-h] input_files [input_files ...]

Print a summary of the contents of pod5 files

positional arguments:
  input_files

optional arguments:
  -h, --help   show this help message and exit

Example: pod5 inspect summary input.pod5
```

### `inspect debug`

```
usage: pod5 inspect debug [-h] input_files

Print debugging information

positional arguments:
  input_files

optional arguments:
  -h, --help   show this help message and exit

Example: pod5 inspect debug input.pod5
```

---

## Manipulation

### `merge`

```
usage: pod5 merge [-h] -o OUTPUT [-r] [-f] [-t THREADS] [-R READERS] [-D]
                  inputs [inputs ...]

Merge multiple pod5 files

positional arguments:
  inputs                Pod5 filepaths to use as inputs

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filepath (default: None)
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -f, --force-overwrite
                        Overwrite destination files (default: False)
  -t THREADS, --threads THREADS
                        Number of workers (default: 4)
  -R READERS, --readers READERS
                        number of merge readers TESTING ONLY (default: 20)
  -D, --duplicate-ok    Allow duplicate read_ids (default: False)

Example: pod5 merge inputs/*.pod5 merged.pod5
```

### `filter`

```
usage: pod5 filter [-h] [-r] [-f] -i IDS -o OUTPUT [-t THREADS] [-M] [-D]
                   inputs [inputs ...]

Take a subset of reads using a list of read_ids from one or more inputs

positional arguments:
  inputs                Pod5 filepaths to use as inputs

optional arguments:
  -h, --help            show this help message and exit
  -r, --recursive       Search for input files recursively matching `*.pod5`
  -f, --force-overwrite
                        Overwrite destination files
  -t THREADS, --threads THREADS
                        Number of workers

required arguments:
  -i IDS, --ids IDS     A file containing a list of only valid read ids to
                        filter from inputs
  -o OUTPUT, --output OUTPUT
                        Destination output filename

content settings:
  -M, --missing-ok      Allow missing read_ids
  -D, --duplicate-ok    Allow duplicate read_ids

Example: pod5 filter inputs*.pod5 --ids read_ids.txt --output filtered.pod5
```

### `subset`

```
usage: pod5 subset [-h] [-o OUTPUT] [-r] [-f] [-t THREADS] [--csv CSV]
                   [-s TABLE] [-R READ_ID_COLUMN] [-c COLUMNS [COLUMNS ...]]
                   [--template TEMPLATE] [-T] [-M] [-D]
                   inputs [inputs ...]

Given one or more pod5 input files, take subsets of reads into one or more pod5 output files by a user-supplied mapping.

positional arguments:
  inputs                Pod5 filepaths to use as inputs

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Destination directory to write outputs
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -f, --force-overwrite
                        Overwrite destination files (default: False)
  -t THREADS, --threads THREADS
                        Number of subsetting workers (default: 4)

direct mapping:
  --csv CSV             CSV file mapping output filename to read ids (default:
                        None)

table mapping:
  -s TABLE, --summary TABLE, --table TABLE
                        Table filepath (csv or tsv) (default: None)
  -R READ_ID_COLUMN, --read-id-column READ_ID_COLUMN
                        Name of the read_id column in the summary (default:
                        read_id)
  -c COLUMNS [COLUMNS ...], --columns COLUMNS [COLUMNS ...]
                        Names of --summary / --table columns to subset on
                        (default: None)
  --template TEMPLATE   template string to generate output filenames (e.g.
                        "mux-{mux}_barcode-{barcode}.pod5"). default is to
                        concatenate all columns to values. (default: None)
  -T, --ignore-incomplete-template
                        Suppress the exception raised if the --template string
                        does not contain every --columns key (default: None)

content settings:
  -M, --missing-ok      Allow missing read_ids (default: False)
  -D, --duplicate-ok    Allow duplicate read_ids (default: False)

Example: pod5 subset inputs.pod5 --output subset_mux/ --summary summary.tsv --columns mux
```

### `repack`

```
usage: pod5 repack [-h] [-o OUTPUT] [-r] [-f] [-t THREADS] inputs [inputs ...]

Repack a pod5 files into a single output

positional arguments:
  inputs                Input pod5 file(s) to repack

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory for pod5 files
  -r, --recursive       Search for input files recursively matching `*.pod5`
  -f, --force-overwrite
                        Overwrite destination files
  -t THREADS, --threads THREADS
                        Number of repacking workers

Example: pod5 repack inputs/*.pod5 repacked/
```

### `recover`

```
usage: pod5 recover [-h] [-r] [-f] inputs [inputs ...]

Attempt to recover pod5 files. Recovered files are written to sibling files
with the '_recovered.pod5` suffix

positional arguments:
  inputs                Input pod5 file(s) to update

optional arguments:
  -h, --help            show this help message and exit
  -r, --recursive       Search for input files recursively matching `*.pod5`
  -f, --force-overwrite
                        Overwrite destination files
```

### `update`

```
usage: pod5 update [-h] -o OUTPUT [-r] [-f] inputs [inputs ...]

Update a pod5 files to the latest available version

positional arguments:
  inputs                Input pod5 file(s) to update

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory for updated pod5 files (default:
                        None)
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -f, --force-overwrite
                        Overwrite destination files (default: False)

Example: pod5 update my.pod5 --output updated/
```

---

## Conversion

### `convert fast5`

```
usage: pod5 convert fast5 [-h] -o OUTPUT [-r] [-t THREADS] [--strict]
                          [-O ONE_TO_ONE] [-f]
                          [--signal-chunk-size SIGNAL_CHUNK_SIZE]
                          inputs [inputs ...]

Convert fast5 file(s) into a pod5 file(s)

positional arguments:
  inputs                Input path for fast5 file

optional arguments:
  -h, --help            show this help message and exit
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -t THREADS, --threads THREADS
                        Set the number of threads to use (default: 4)
  --strict              Immediately quit if an exception is encountered during
                        conversion instead of continuing with remaining inputs
                        after issuing a warning (default: False)

required arguments:
  -o OUTPUT, --output OUTPUT
                        Output path for the pod5 file(s). This can be an
                        existing directory (creating 'output.pod5' within it)
                        or a new named file path. A directory must be given
                        when using --one-to-one. (default: None)

output control arguments:
  -O ONE_TO_ONE, --one-to-one ONE_TO_ONE
                        Output files are written 1:1 to inputs. 1:1 output
                        files are written to the output directory in a new
                        directory structure relative to the directory path
                        provided to this argument. This directory path must be
                        a relative parent of all inputs. (default: None)
  -f, --force-overwrite
                        Overwrite destination files (default: False)
  --signal-chunk-size SIGNAL_CHUNK_SIZE
                        Chunk size to use for signal data set (default:
                        102400)
```

### `convert to_fast5`

```
usage: pod5 convert to_fast5 [-h] -o OUTPUT [-r] [-t THREADS] [-f]
                             [--file-read-count FILE_READ_COUNT]
                             inputs [inputs ...]

Convert pod5 file(s) into fast5 file(s)

positional arguments:
  inputs

optional arguments:
  -h, --help            show this help message and exit
  -r, --recursive       Search for input files recursively matching `*.pod5`
                        (default: False)
  -t THREADS, --threads THREADS
                        How many file writers to keep active (default: 4)

required arguments:
  -o OUTPUT, --output OUTPUT
                        Output path for the pod5 file(s). This can be an
                        existing directory (creating 'output.pod5' within it)
                        or a new named file path. (default: None)

output control arguments:
  -f, --force-overwrite
                        Overwrite destination files (default: False)
  --file-read-count FILE_READ_COUNT
                        Number of reads to write per file. (default: 4000)
```
