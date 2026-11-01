# genomepy and gget: genome catalogues, annotations and database queries

This group of functions in `dnakit.databases` is a set of adapters over two optional packages and a few external programs: [genomepy](https://github.com/vanheeringen-lab/genomepy) (genome catalogues of the GENCODE, Ensembl, UCSC and NCBI providers, annotation download and conversion, NCBI assembly reports, aligner index commands) and [gget](https://github.com/pachterlab/gget) (Ensembl reference links, NCBI Virus and GenBank downloads, cBioPortal data and heatmaps, CELLxGENE Census queries, COSMIC reference files, and the MUSCLE 5 and DIAMOND alignment programs). Five further functions are offline transforms of genomepy catalogue records that need no package at all. The gget functions return the service name and the access date next to the data; DNAKit itself caches nothing (genomepy keeps its own on-disk cache of provider catalogues).

genomepy is MIT-licensed and gget is BSD-2-Clause (from the package metadata). MUSCLE, DIAMOND, the aligners and the UCSC conversion tools are separate programs with their own licences; DNAKit starts them as external processes and copies none of their code. The COSMIC download functions need your own COSMIC account and use it only for the request.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises `BackendUnavailableError` that includes the install command; the aligner and MUSCLE/DIAMOND wrappers raise it when the program is not found. Install the extras with `dnakit[bio-genomepy]`, `dnakit[bio-gget]` and `dnakit[bio-cellxgene]`. The examples import from `dnakit.databases`.

## Example data

No example file is declared. The examples that need inputs write tiny files into a temporary directory, and the functions that need internet access, an account or an external program that is not installed here show the error they raise for invalid input or a missing program. Each of those sections adds an unexecuted "Usage with the dependency installed" block with a realistic call.

## 1) `dnakit.databases.genomepy_bowtie2_index_command()`

- **Purpose:** Run a `bowtie2-build` index command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `bowtie2-build`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `bowtie2` index plugin step as an argument vector.
- **API:** `genomepy_bowtie2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bowtie2-build` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: The example passes a command that does not start with `bowtie2-build`; the check fails before any program is looked up or run, so it gives the same output on every machine.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `bowtie2-build`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`bowtie2_index_command`), `backend` (`bowtie2`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
from dnakit.databases import genomepy_bowtie2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_bowtie2_index_command(
        ["samtools", "faidx", "genome.fa"], working_directory=".", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError The bowtie2 index command must start with 'bowtie2-build'.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_bowtie2_index_command

result = genomepy_bowtie2_index_command(
    ["bowtie2-build", "genome.fa", "genome"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 2) `dnakit.databases.genomepy_bwa_index_command()`

- **Purpose:** Run a `bwa index` command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `bwa`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `bwa` index plugin step as an argument vector.
- **API:** `genomepy_bwa_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bwa` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: bwa is not installed here, so the example writes a tiny shell script `fake_bwa.sh` that creates two index files and passes it as `executable_path`. The function itself runs for real (argument check, subprocess, detection of new files); only the aligner is replaced.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `bwa`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`bwa_index_command`), `backend` (`bwa`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_bwa_index_command

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("genome.fa").write_text(">chr1\nACGTACGT\n")
    fake = Path("fake_bwa.sh")
    fake.write_text('#!/bin/sh\ntouch "$2.bwt" "$2.sa"\n')
    fake.chmod(0o755)
    result = genomepy_bwa_index_command(
        ["bwa", "index", "genome.fa"],
        working_directory=".",
        executable_path=str(fake.resolve()),
        verbose=False,
    )
    print(result["status"], result["operation"], result["backend"])
    print([Path(path).name for path in result["outputs"].values()])
```

- **Example output:**

```text
complete bwa_index_command bwa
['genome.fa.bwt', 'genome.fa.sa']
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_bwa_index_command

result = genomepy_bwa_index_command(
    ["bwa", "index", "genome.fa"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 3) `dnakit.databases.genomepy_download_annotation()`

- **Purpose:** Download or copy an annotation file into a genome folder and convert it to GTF and BED files, as genomepy's base and GENCODE providers do.
- **Method:** Validates the arguments, checks that the UCSC tool `genePredToGtf` is on PATH (otherwise raises `BackendUnavailableError`), then calls `genomepy.providers.base.download_annotation`. genomepy writes `<genomes_dir>/<localname>/<localname>.annotation.gtf` and `.annotation.bed` (the conversion uses `gtfToGenePred`, `gff3ToGenePred`, `bedToGenePred`, `genePredToGtf` and `genePredToBed`, depending on the input type). DNAKit then lists every file under that folder.
- **API:** `genomepy_download_annotation(genomes_dir: 'str', annotation_url: 'str', localname: 'str', *, n: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`), the UCSC conversion tools on PATH (for example from bioconda) and, for a remote URL, internet access. The UCSC tools are not installed here, so the example shows the validation error. Example: The example uses an empty `localname`, which is rejected before any tool is checked or any download starts.
- **Input:** `genomes_dir`: base folder; `annotation_url`: URL (or an existing local path, which genomepy copies); `localname`: non-empty folder name; `n`: optional positive integer, download only the first lines. Text arguments are at most 512 characters.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings). `outputs` holds `directory` and `file_0`, `file_1`, ... (sorted paths of all files in the folder).
- **Example code:**

```python
from dnakit.databases import genomepy_download_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_download_annotation("genomes", "https://example.org/a.gtf.gz", "", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError localname must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_download_annotation

result = genomepy_download_annotation(
    "genomes",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz",
    "GRCh38_gencode",
)
print(result["outputs"]["directory"])
```

## 4) `dnakit.databases.genomepy_download_annotation_head()`

- **Purpose:** Download only the first lines of a (gzipped) annotation file into a new file.
- **Method:** Validates the URL and the new output path (its parent folder is created), then calls `genomepy.providers.base.download_head`. That function opens the URL with `urllib`, decompresses it when the URL ends in `.gz`, writes every line, and stops after `n + 2` lines that do not start with `#`; comment lines are written but not counted. The returned `lines` is the requested `n`, not the number of lines written.
- **API:** `genomepy_download_annotation_head(annotation_url: 'str', annotation_file: 'str', *, n: 'int' = 5, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`). A remote URL needs internet access; a `file://` URL works offline, which the example uses. Example: The example reads a local GTF through a `file://` URL, so nothing touches the network. With `n=2` the file gets the comment line plus four feature lines (`n + 2`), five lines in total.
- **Input:** `annotation_url`: URL string of at most 512 characters; `annotation_file`: output path that must not exist yet; `n`: positive integer (default 5).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings). `outputs` holds `head` (the new file); `lines` is the requested `n`.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_download_annotation_head

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    rows = ["#!genome-build GRCh38"]
    for i in range(1, 8):
        rows.append(f'chr1\tsrc\tgene\t{i}\t{i + 9}\t.\t+\t.\tgene_id "g{i}";')
    Path("annotation.gtf").write_text("\n".join(rows) + "\n")
    url = Path("annotation.gtf").resolve().as_uri()
    result = genomepy_download_annotation_head(url, "head/head.gtf", n=2, verbose=False)
    print(result["operation"], result["lines"], list(result["outputs"]))
    print(len(Path("head/head.gtf").read_text().splitlines()))
```

- **Example output:**

```text
annotation_head 2 ['head']
5
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_download_annotation_head

result = genomepy_download_annotation_head(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz",
    "head/gencode_head.gtf",
    n=5,
)
print(result["outputs"]["head"])
```

## 5) `dnakit.databases.genomepy_download_assembly_report()`

- **Purpose:** Retrieve the NCBI assembly report (chromosome names, roles and lengths) of an assembly accession.
- **Method:** Calls `genomepy.providers.ncbi.download_assembly_report`. Without `output_path` the report rows are returned; with `output_path` (a new path) genomepy writes the report as a TSV file. upstream returns `None` without any request when the accession does not start with `GCA` or `GCF` or when NCBI has no such assembly; DNAKit turns that into `BackendExecutionError`. DNAKit itself only checks that the accession is a non-empty string.
- **API:** `genomepy_download_assembly_report(accession: 'str', *, output_path: 'str | None' = None, quiet: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. The example uses an accession that upstream rejects before any request, so it runs offline.
- **Input:** `accession`: string, normally `GCA_...` or `GCF_...`; `output_path`: optional new file path; `quiet`: silence genomepy warnings (default True).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `accession`, and either `rows` (a list of dictionaries with the report columns) and `row_count`, or, with `output_path`, `outputs` holding `report`.
- **Example code:**

```python
from dnakit.databases import genomepy_download_assembly_report
from dnakit.exceptions import BackendExecutionError

try:
    genomepy_download_assembly_report("XYZ_1", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
BackendExecutionError No NCBI assembly report for XYZ_1.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_download_assembly_report

result = genomepy_download_assembly_report("GCF_000001405.40")
print(result["row_count"], result["rows"][0]["Sequence-Name"])
```

## 6) `dnakit.databases.genomepy_ensembl_genomes()`

- **Purpose:** Build the Ensembl genome catalogue from the Ensembl REST API.
- **Method:** Calls `genomepy.providers.ensembl.get_genomes(rest_url)`, which asks the REST service for the divisions and the genomes in each division (EnsemblBacteria is skipped), keeps a few summary fields per genome and adds the fixed GRCh37 record. genomepy caches the result on disk.
- **API:** `genomepy_ensembl_genomes(rest_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. Example: The example passes an empty URL, which is rejected before any request is made.
- **Input:** `rest_url`: non-empty URL string of at most 512 characters, for example `https://rest.ensembl.org/`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `genomes` (assembly name mapped to its record) and `genome_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_ensembl_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ensembl_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError rest_url must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_ensembl_genomes

result = genomepy_ensembl_genomes("https://rest.ensembl.org/")
print(result["genome_count"], "GRCh37" in result["genomes"])
```

## 7) `dnakit.databases.genomepy_find_annotation_urls()`

- **Purpose:** Look for GTF or GFF3 annotation files in the remote folder that holds a genome file.
- **Method:** Calls `genomepy.providers.url.search_url_for_annotations(url)`: it takes the folder part of the genome URL, reads that folder listing and matches file names against the genome's local name with a fuzzy search. The links of the matching files are returned (an empty list when none matches).
- **API:** `genomepy_find_annotation_urls(url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. Example: The example passes an empty URL, which is rejected before any request is made.
- **Input:** `url`: non-empty URL string (at most 512 characters) of a genome file.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `urls` (list of annotation links) and `url_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_find_annotation_urls
from dnakit.exceptions import ConfigurationError

try:
    genomepy_find_annotation_urls("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError url must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_find_annotation_urls

result = genomepy_find_annotation_urls(
    "https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz"
)
print(result["urls"])
```

## 8) `dnakit.databases.genomepy_gencode_genomes()`

- **Purpose:** Build the GENCODE genome catalogue: the latest release of each major assembly.
- **Method:** Calls `genomepy.providers.gencode.get_genomes(ftp_link)`. It connects to the FTP site, lists the human and mouse release folders, takes the releases in order (highest first, using the same release rule as `genomepy_gencode_releases`), finds the primary-assembly file of each release and records the first release found for each assembly, then adds the GRCh37 record for human.
- **API:** `genomepy_gencode_genomes(ftp_link: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and FTP access to the GENCODE server. Example: The example passes an empty link, which is rejected before any connection is made.
- **Input:** `ftp_link`: non-empty FTP link string of at most 512 characters, for example `ftp://ftp.ebi.ac.uk/pub/databases/gencode`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `genomes` (assembly name mapped to a record with `annotations`, `taxonomy_id`, `species`, `text_search`) and `genome_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_gencode_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_gencode_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError ftp_link must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_gencode_genomes

result = genomepy_gencode_genomes("ftp://ftp.ebi.ac.uk/pub/databases/gencode")
print(sorted(result["genomes"]))
```

## 9) `dnakit.databases.genomepy_gencode_releases()`

- **Purpose:** Extract the GENCODE release numbers from a listing of release folders, highest first.
- **Method:** Checks the arguments and calls `genomepy.providers.gencode.get_releases(listing, species)`. For every entry whose last two characters are digits the number is taken; numbers of 21 or below are dropped, the rest are sorted from high to low, and for `mouse` an `M` is put in front of each. Entries such as `README.txt` or `release_M9` are ignored. The function is pure and makes no request.
- **API:** `genomepy_gencode_releases(listing: 'Sequence[str]', species: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) for the helper it calls; no internet access.
- **Input:** `listing`: a list of strings (a single string is refused), usually paths ending in `release_NN` or `release_MNN`; `species`: `human` or `mouse`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `releases` (list of strings, highest first) and `release_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_gencode_releases

human = genomepy_gencode_releases(
    ["release_44", "release_45", "release_19", "README.txt", "release_43"],
    "human",
    verbose=False,
)
print(human["releases"], human["release_count"])
mouse = genomepy_gencode_releases(
    ["release_M33", "release_M34", "release_M9"], "mouse", verbose=False
)
print(mouse["releases"])
```

- **Example output:**

```text
['45', '44', '43'] 3
['M34', 'M33']
```

## 10) `dnakit.databases.genomepy_gencode_to_ucsc_names()`

- **Purpose:** Map GENCODE (Ensembl-style) assembly names to UCSC names for a catalogue.
- **Method:** Calls `genomepy.providers.gencode.get_gencode2ucsc(genomes)`. The mapping starts with `GRCh37` to `hg19` and `GRCm38` to `mm10`; for every other assembly the species field decides the prefix (`hg` if it starts with `H`, otherwise `mm`) and the digits of the assembly name give the number. No request is made. A record without a `species` field makes genomepy raise `KeyError`, which DNAKit does not convert.
- **API:** `genomepy_gencode_to_ucsc_names(genomes: 'Mapping[str, Any]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`); no internet access.
- **Input:** `genomes`: a non-empty mapping of assembly name to a record that has a `species` string, for example a catalogue from `genomepy_gencode_genomes`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `mapping` (GENCODE name mapped to UCSC name) and `entry_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_gencode_to_ucsc_names

catalogue = {
    "GRCh38": {"species": "Homo sapiens"},
    "GRCm39": {"species": "Mus musculus"},
    "GRCh37": {"species": "Homo sapiens"},
}
result = genomepy_gencode_to_ucsc_names(catalogue, verbose=False)
print(result["mapping"], result["entry_count"])
```

- **Example output:**

```text
{'GRCh37': 'hg19', 'GRCm38': 'mm10', 'GRCh38': 'hg38', 'GRCm39': 'mm39'} 4
```

## 11) `dnakit.databases.genomepy_generate_annotation()`

- **Purpose:** Convert a GTF annotation to BED, or a BED annotation to GTF, with the UCSC tools.
- **Method:** Requires `template` to be an existing file and refuses an existing `target` unless `overwrite` is True; the folder of `target` is created. It checks that the UCSC tool `genePredToGtf` is on PATH (otherwise `BackendUnavailableError`) and calls `genomepy.annotation.utils.generate_annot`, which pipes the template through `gtfToGenePred` and `genePredToBed` (GTF to BED) or `bedToGenePred` and `genePredToGtf` (BED to GTF). A template whose name contains neither `gtf` nor `bed` is reported as `ConfigurationError`; a missing output file as `BackendExecutionError`.
- **API:** `genomepy_generate_annotation(template: 'str', target: 'str', *, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and the UCSC conversion tools on PATH (for example `conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-genepredtogtf ucsc-bedtogenepred`). No internet access. The UCSC tools are not installed here, so the example shows the input check. Example: The example names a template that does not exist; the text after the colon in the message is a machine-specific path, so only the part before it is printed.
- **Input:** `template`: path of an existing GTF or BED file; `target`: output path; `overwrite`: allow replacing `target` (default False).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings). `outputs` holds `annotation` (the converted file).
- **Example code:**

```python
from dnakit.databases import genomepy_generate_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_generate_annotation("no_such_file.gtf", "out.bed", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1].split(":")[0])
```

- **Example output:**

```text
ConfigurationError template does not exist
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_generate_annotation

result = genomepy_generate_annotation("annotation.gtf", "annotation.bed")
print(result["outputs"]["annotation"])
```

## 12) `dnakit.databases.genomepy_genome_size()`

- **Purpose:** Return the total sequence length of an NCBI assembly, summed from its assembly report.
- **Method:** Calls `genomepy.providers.ncbi.get_genome_size(accession)`, which downloads the assembly report and sums the `Sequence-Length` column; the same implementation serves the ncbi, gencode and ucsc providers. Upstream returns -1 when there is no report (for example for an accession that does not start with `GCA` or `GCF`, which is decided without any request); DNAKit reports that as `genome_size` of `None`.
- **API:** `genomepy_genome_size(accession: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. The example uses an accession that upstream rejects before any request, so it runs offline.
- **Input:** `accession`: non-empty string of at most 512 characters, normally `GCA_...` or `GCF_...`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `accession` and `genome_size` (integer, or `None` when no report exists).
- **Example code:**

```python
from dnakit.databases import genomepy_genome_size

result = genomepy_genome_size("XYZ_1", verbose=False)
print(result["operation"], result["accession"], result["genome_size"])
```

- **Example output:**

```text
genome_size XYZ_1 None
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_genome_size

result = genomepy_genome_size("GCF_000001405.40")
print(result["genome_size"])
```

## 13) `dnakit.databases.genomepy_gmap_index_command()`

- **Purpose:** Run a `gmap_build` index command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `gmap_build`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `gmap` index plugin step as an argument vector.
- **API:** `genomepy_gmap_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `gmap_build` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: The example points `working_directory` at a folder that does not exist; the check fails before the program is looked up. The text after the colon is a machine-specific path, so it is cut off.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `gmap_build`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`gmap_index_command`), `backend` (`gmap`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
from dnakit.databases import genomepy_gmap_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_gmap_index_command(
        ["gmap_build", "-d", "genome", "genome.fa"],
        working_directory="no_such_directory",
        verbose=False,
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1].split(":")[0])
```

- **Example output:**

```text
ConfigurationError working_directory does not exist
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_gmap_index_command

result = genomepy_gmap_index_command(
    ["gmap_build", "-D", "indexes", "-d", "genome", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 14) `dnakit.databases.genomepy_head_annotations()`

- **Purpose:** Show the first annotation lines of a genome, for each annotation style the providers offer.
- **Method:** Validates `name` and `lines`, checks that the UCSC tool `genePredToGtf` is on PATH (otherwise `BackendUnavailableError`), and runs `genomepy.functions.head_annotations(name, provider, n=lines)` with standard output captured. genomepy visits the online providers that list the genome, downloads the first `lines` genes of each annotation into a temporary folder, converts them to GTF and prints up to `lines` non-comment GTF lines. The non-empty printed lines are returned.
- **API:** `genomepy_head_annotations(name: 'str', *, provider: 'str | None' = None, lines: 'int' = 2, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`), the UCSC conversion tools on PATH and internet access. The UCSC tools are not installed here, so the example shows the input check. Example: The example passes an empty genome name, which is rejected before any tool check or download.
- **Input:** `name`: genome name such as `hg38`; `provider`: optional provider name; `lines`: positive integer (default 2).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `genome` and `lines` (the list of printed lines).
- **Example code:**

```python
from dnakit.databases import genomepy_head_annotations
from dnakit.exceptions import ConfigurationError

try:
    genomepy_head_annotations("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError name must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_head_annotations

result = genomepy_head_annotations("hg38", provider="UCSC", lines=2)
print(result["lines"])
```

## 15) `dnakit.databases.genomepy_hisat2_index_command()`

- **Purpose:** Run a `hisat2-build` index command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `hisat2-build`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `hisat2` index plugin step as an argument vector.
- **API:** `genomepy_hisat2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `hisat2-build` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: The example passes the command as one shell string, which the function refuses; no program is looked up or run.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `hisat2-build`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`hisat2_index_command`), `backend` (`hisat2`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
from dnakit.databases import genomepy_hisat2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_hisat2_index_command(
        "hisat2-build genome.fa genome", working_directory=".", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError command must be a non-empty list of strings (no shell string).
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_hisat2_index_command

result = genomepy_hisat2_index_command(
    ["hisat2-build", "genome.fa", "genome"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 16) `dnakit.databases.genomepy_list_available_genomes()`

- **Purpose:** List the genomes offered by genomepy's online providers.
- **Method:** Calls `genomepy.functions.list_available_genomes(provider, size)`, which asks each online provider (GENCODE, Ensembl, UCSC, NCBI) for its catalogue, and converts the rows to plain lists. Each row is `[name, provider, description...]`. DNAKit adds no validation: an unknown provider name makes genomepy raise `ValueError`, which is not converted.
- **API:** `genomepy_list_available_genomes(*, provider: 'str | None' = None, include_size: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. The example uses an unknown provider name, which genomepy rejects before any request.
- **Input:** `provider`: optional provider name (default: all online providers); `include_size`: add the absolute genome size to each row (default False).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `rows` (list of lists) and `row_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_list_available_genomes

try:
    genomepy_list_available_genomes(provider="nonexistent_provider", verbose=False)
except ValueError as error:
    print(type(error).__name__, str(error).split(".")[0])
```

- **Example output:**

```text
ValueError Unknown provider 'nonexistent_provider'
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_list_available_genomes

result = genomepy_list_available_genomes(provider="UCSC")
print(result["row_count"], result["rows"][0][:2])
```

## 17) `dnakit.databases.genomepy_minimap2_index_command()`

- **Purpose:** Run a `minimap2` index command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `minimap2`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `minimap2` index plugin step as an argument vector.
- **API:** `genomepy_minimap2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `minimap2` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: The example passes an empty command list, which the function refuses; no program is looked up or run.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `minimap2`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`minimap2_index_command`), `backend` (`minimap2`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
from dnakit.databases import genomepy_minimap2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_minimap2_index_command([], working_directory=".", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError command must be a non-empty list of strings (no shell string).
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_minimap2_index_command

result = genomepy_minimap2_index_command(
    ["minimap2", "-d", "genome.mmi", "genome.fa"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 18) `dnakit.databases.genomepy_ncbi_genomes()`

- **Purpose:** Parse NCBI assembly summary files into a genome catalogue.
- **Method:** Calls `genomepy.providers.ncbi.get_genomes(assembly_url)`. genomepy opens four files under the URL with `urllib` (`assembly_summary_genbank_historical.txt`, `assembly_summary_refseq_historical.txt`, `assembly_summary_genbank.txt`, `assembly_summary_refseq.txt`), skips the comment line, reads the header line, and keeps eight columns of each row. Rows whose `ftp_path` is `na` are dropped. Later files overwrite earlier ones with the same assembly name (column 15, made safe for use as a name).
- **API:** `genomepy_ncbi_genomes(assembly_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`). The real NCBI summaries are large and need internet access; a `file://` URL to local files works offline, which the example uses. Example: The example writes four tiny summary files (one assembly in the RefSeq file) and passes their folder as a `file://` URL, so nothing touches the network.
- **Input:** `assembly_url`: non-empty URL string (at most 512 characters) of the folder that holds the four summary files, for example `https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `genomes` (assembly name mapped to a record with `assembly_accession`, `taxid`, `species_taxid`, `organism_name`, `asm_submitter`, `gbrs_paired_asm`, `paired_asm_comp`, `ftp_path`) and `genome_count`.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_ncbi_genomes

COLUMNS = [
    "assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category",
    "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate",
    "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date",
    "asm_name", "asm_submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
]
FILES = [
    "assembly_summary_genbank_historical.txt",
    "assembly_summary_refseq_historical.txt",
    "assembly_summary_genbank.txt",
    "assembly_summary_refseq.txt",
]

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("summaries").mkdir()
    record = dict.fromkeys(COLUMNS, "na")
    record.update(
        assembly_accession="GCF_000001405.40", taxid="9606", species_taxid="9606",
        organism_name="Homo sapiens", asm_name="GRCh38.p14", asm_submitter="Genome Reference Consortium",
        gbrs_paired_asm="GCA_000001405.29", paired_asm_comp="identical",
        ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
    )
    for name in FILES:
        text = "# comment line\n#" + "\t".join(COLUMNS) + "\n"
        if name == "assembly_summary_refseq.txt":
            text += "\t".join(record[column] for column in COLUMNS) + "\n"
        Path("summaries", name).write_text(text)
    url = Path("summaries").resolve().as_uri()
    result = genomepy_ncbi_genomes(url, verbose=False)
    genome = result["genomes"]["GRCh38.p14"]
    print(result["genome_count"], list(result["genomes"]))
    print(genome["assembly_accession"], genome["organism_name"], genome["taxid"])
```

- **Example output:**

```text
1 ['GRCh38.p14']
GCF_000001405.40 Homo sapiens 9606
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_ncbi_genomes

result = genomepy_ncbi_genomes("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS")
print(result["genome_count"])
```

## 19) `dnakit.databases.genomepy_run_index_cmd()`

- **Purpose:** Run the index command of one named aligner (bowtie2, bwa, gmap, hisat2, minimap2 or star) as an argument list.
- **Method:** Checks that `name` is one of the six aligners and then runs the same step as the dedicated `genomepy_<aligner>_index_command` wrapper: the first element of `command` must be named after the aligner's index program (`bowtie2-build`, `bwa`, `gmap_build`, `hisat2-build`, `minimap2`, `STAR`), the program is resolved from `executable_path` or PATH, run without a shell in `working_directory`, and the new files are listed. Unlike genomepy's own `run_index_cmd`, which hands a string to a shell, it accepts only an argument list.
- **API:** `genomepy_run_index_cmd(name: 'str', command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program of the chosen aligner on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: The example asks for an aligner name that is not in the list; the check fails before any program is looked up.
- **Input:** `name`: one of `bowtie2`, `bwa`, `gmap`, `hisat2`, `minimap2`, `star`; `command`, `working_directory` and `executable_path` as in the dedicated wrappers.
- **Output:** A dictionary with `status` (`complete`), `operation` (`<name>_index_command`), `backend` (the aligner name), `version`, `outputs` (`index_0`, ... for the new files), `executable` and `log_tail`.
- **Example code:**

```python
from dnakit.databases import genomepy_run_index_cmd
from dnakit.exceptions import ConfigurationError

try:
    genomepy_run_index_cmd("blast", ["makeblastdb"], working_directory=".", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError name must be one of ['bowtie2', 'bwa', 'gmap', 'hisat2', 'minimap2', 'star'].
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_run_index_cmd

result = genomepy_run_index_cmd(
    "star",
    ["STAR", "--runMode", "genomeGenerate", "--genomeDir", "star_index",
     "--genomeFastaFiles", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 20) `dnakit.databases.genomepy_search_genomes()`

- **Purpose:** Search genome names, descriptions, accession IDs and taxonomy IDs across the online providers.
- **Method:** Calls `genomepy.providers.search(term, provider, exact, size)`. genomepy searches the names and descriptions case-insensitively (a regular expression is allowed; several words in any order), matches a term that starts with `GCA_` or `GCF_` against the assembly accession, and a number against the taxonomy ID. Each row is `[name, provider, description...]`.
- **API:** `genomepy_search_genomes(term: 'str', *, provider: 'str | None' = None, exact: 'bool' = False, include_size: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. Example: The example passes an empty search term, which is rejected before any provider is contacted.
- **Input:** `term`: non-empty string of at most 512 characters; `provider`: optional provider name; `exact`: require an exact match (default False); `include_size`: add the genome size to each row (default False).
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `rows` (list of lists) and `row_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_search_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_search_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError term must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_search_genomes

result = genomepy_search_genomes("zebrafish", provider="Ensembl")
print(result["row_count"], result["rows"][0][:2])
```

## 21) `dnakit.databases.genomepy_star_index_command()`

- **Purpose:** Run a STAR genome-generation command given as an argument list and report the new index files.
- **Method:** Checks that `command` is a list and that its first element is named `STAR`, resolves the program from `executable_path` or PATH, runs the remaining arguments in `working_directory` without a shell (bounded by DNAKit's runner: one thread, 3600 s timeout, 1 MB output cap), and lists the files that appeared. It raises `BackendExecutionError` when the program wrote no new file. It does not import the genomepy package; it reproduces genomepy's `star` index plugin step as an argument vector.
- **API:** `genomepy_star_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `STAR` on PATH, or its path in `executable_path`. The genomepy package is not needed. Example: STAR is not installed here, so the example writes a tiny shell script `fake_star.sh` that creates three files and passes it as `executable_path`. The function itself runs for real; only the aligner is replaced.
- **Input:** `command`: a non-empty list of strings (a shell string is refused) whose first element has the file name `STAR`; `working_directory`: an existing directory in which the program runs; `executable_path`: optional path of the program (default: found on PATH).
- **Output:** A dictionary with `status` (`complete`), `operation` (`star_index_command`), `backend` (`star`), `version`, `outputs` (`index_0`, `index_1`, ... mapped to the paths of the files that are new in `working_directory`, searched recursively), `executable` (resolved program path) and `log_tail` (last 2000 characters of the program output).
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_star_index_command

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("genome.fa").write_text(">chr1\nACGTACGT\n")
    fake = Path("fake_star.sh")
    fake.write_text('#!/bin/sh\ntouch Genome SA SAindex\n')
    fake.chmod(0o755)
    result = genomepy_star_index_command(
        ["STAR", "--runMode", "genomeGenerate", "--genomeFastaFiles", "genome.fa"],
        working_directory=".",
        executable_path=str(fake.resolve()),
        verbose=False,
    )
    print(result["operation"], result["backend"])
    print([Path(path).name for path in result["outputs"].values()])
```

- **Example output:**

```text
star_index_command star
['Genome', 'SA', 'SAindex']
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_star_index_command

result = genomepy_star_index_command(
    ["STAR", "--runMode", "genomeGenerate", "--genomeDir", "star_index",
     "--genomeFastaFiles", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 22) `dnakit.databases.genomepy_ucsc_download_annotation()`

- **Purpose:** Download a UCSC gene annotation (extended genePred) from the UCSC MySQL server and convert it to GTF and BED files.
- **Method:** Validates the text arguments, checks that the UCSC tool `genePredToGtf` is on PATH (otherwise `BackendUnavailableError`), and calls `genomepy.providers.ucsc.download_annotation`. genomepy reads the column list of the annotation table, queries the table rows from the MySQL server, writes an extended genePred file and converts it with `genePredToGtf` and `genePredToBed` into `<genomes_dir>/<localname>/`. DNAKit lists every file in that folder.
- **API:** `genomepy_ucsc_download_annotation(name: 'str', annotation: 'str', genomes_dir: 'str', localname: 'str', *, n: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`), the UCSC conversion tools on PATH and internet access to the UCSC MySQL server. The UCSC tools are not installed here, so the example shows the input check. Example: The example passes an empty genome name, which is rejected before any tool check or download.
- **Input:** `name`: UCSC genome (database) name such as `hg38`; `annotation`: table name such as `ncbiRefSeq`, `refGene`, `ensGene` or `knownGene`; `genomes_dir`: base folder; `localname`: output folder name; `n`: optional positive integer limiting the rows. Text arguments are non-empty and at most 512 characters.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings). `outputs` holds `directory` and `file_0`, `file_1`, ... (sorted paths of the files in the folder).
- **Example code:**

```python
from dnakit.databases import genomepy_ucsc_download_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ucsc_download_annotation("", "refGene", "genomes", "hg38", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError name must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_ucsc_download_annotation

result = genomepy_ucsc_download_annotation("hg38", "refGene", "genomes", "hg38_refgene")
print(result["outputs"]["directory"])
```

## 23) `dnakit.databases.genomepy_ucsc_genomes()`

- **Purpose:** Build the UCSC genome catalogue from the UCSC REST API.
- **Method:** Calls `genomepy.providers.ucsc.get_genomes(rest_url)`, which requests the JSON list of UCSC genomes, keeps five fields per genome (`taxId`, `description`, `scientificName`, `sourceName`, `htmlPath`) and adds the annotation tables and assembly accessions found on the UCSC MySQL server (the record transforms are available offline as `genomepy_ucsc_assembly_accessions`, `genomepy_ucsc_sourcename_accessions` and `genomepy_ucsc_annotation_links`).
- **API:** `genomepy_ucsc_genomes(rest_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomepy (extra `dnakit[bio-genomepy]`) and internet access. Example: The example passes an empty URL, which is rejected before any request is made.
- **Input:** `rest_url`: non-empty URL string of at most 512 characters, for example `https://api.genome.ucsc.edu/list/ucscGenomes`.
- **Output:** A dictionary with the common genomepy result fields `status` (`complete`), `operation`, `backend`, `version` (installed genomepy version, or `None`) and `outputs` (names mapped to path strings), `genomes` (UCSC assembly name mapped to its record) and `genome_count`.
- **Example code:**

```python
from dnakit.databases import genomepy_ucsc_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ucsc_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError rest_url must be a non-empty string of at most 512 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import genomepy_ucsc_genomes

result = genomepy_ucsc_genomes("https://api.genome.ucsc.edu/list/ucscGenomes")
print(result["genome_count"], "hg38" in result["genomes"])
```

## 24) `dnakit.databases.genomepy_ucsc_assembly_accessions()`

- **Purpose:** Set `assembly_accession` on a UCSC catalogue from assembly-equivalence rows.
- **Method:** Offline version of the step `add_accessions1` in genomepy's UCSC provider, which reads the `asmEquivalent` table from the UCSC MySQL server; here the rows are passed in. For each source assembly only the rows with the highest match count are kept; the `GCA_`/`GCF_` accession is cut from the destination name; when an assembly keeps several rows, only `GCA_` accessions are used. Sources that are not in `genomes` are ignored, and a kept row whose destination holds no accession assigns nothing (upstream stores NaN). The input is copied, not changed.
- **API:** `genomepy_ucsc_assembly_accessions(genomes: 'Mapping[str, Mapping[str, Any]]', equivalences: 'Sequence[Sequence[Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`.
- **Requirements:** Needs no optional package and no internet access; works after plain pip install. Example: `hg19` keeps two rows with the top count 900, so only the `GCA_` accession is used; `zzz` is not in the catalogue and is ignored; the input catalogue stays empty.
- **Input:** `genomes`: mapping of UCSC assembly name to a record (a mapping); `equivalences`: rows `(source, destination, matchCount)` with exactly three fields. Wrong types raise `ConfigurationError`.
- **Output:** A new catalogue (dictionary) with the same keys; records that received an accession have an `assembly_accession` string.
- **Example code:**

```python
from dnakit.databases import genomepy_ucsc_assembly_accessions

catalogue = {"hg19": {}, "hg38": {}, "mm10": {}}
rows = [
    ("hg19", "GCA_000001405.1", 900),
    ("hg19", "GCF_000001405.13", 900),
    ("hg19", "GCA_000000000.1", 100),
    ("hg38", "GCF_000001405.39", 5000),
    ("mm10", "GCA_000001635.2", 10),
    ("zzz", "GCA_000000001.1", 3),
]
result = genomepy_ucsc_assembly_accessions(catalogue, rows, verbose=False)
for name, record in result.items():
    print(name, record["assembly_accession"])
print(catalogue)
```

- **Example output:**

```text
hg19 GCA_000001405.1
hg38 GCF_000001405.39
mm10 GCA_000001635.2
{'hg19': {}, 'hg38': {}, 'mm10': {}}
```

## 25) `dnakit.databases.genomepy_ucsc_sourcename_accessions()`

- **Purpose:** Set `assembly_accession` on a UCSC catalogue from the accession found in each record's `sourceName`.
- **Method:** Offline version of the step `add_accessions2` in genomepy's UCSC provider. For every record the first `GCA_`/`GCF_` accession in `sourceName` (nine digits and a version) is stored. A record without `sourceName` raises `ConfigurationError` (upstream: `KeyError`); a `sourceName` of `None` is skipped (upstream: `TypeError`); a non-string value raises `ConfigurationError`. The input is copied, not changed.
- **API:** `genomepy_ucsc_sourcename_accessions(genomes: 'Mapping[str, Mapping[str, Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`.
- **Requirements:** Needs no optional package and no internet access; works after plain pip install.
- **Input:** `genomes`: mapping of assembly name to a record that has a `sourceName` field (string or `None`).
- **Output:** A new catalogue (dictionary); records whose `sourceName` contains an accession gain `assembly_accession`.
- **Example code:**

```python
from dnakit.databases import genomepy_ucsc_sourcename_accessions
from dnakit.exceptions import ConfigurationError

catalogue = {
    "mm39": {"sourceName": "Genome Reference Consortium GCA_000001635.9"},
    "hg38": {"sourceName": "GRCh38 GCF_000001405.40"},
    "unknown": {"sourceName": None},
}
result = genomepy_ucsc_sourcename_accessions(catalogue, verbose=False)
for name, record in result.items():
    print(name, record.get("assembly_accession"))
print("assembly_accession" in catalogue["mm39"])
try:
    genomepy_ucsc_sourcename_accessions({"bad": {}}, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
mm39 GCA_000001635.9
hg38 GCF_000001405.40
unknown None
False
ConfigurationError genome record 'bad' has no sourceName field.
```

## 26) `dnakit.databases.genomepy_ucsc_annotation_links()`

- **Purpose:** Append the annotation table names that UCSC lists for each genome to the `annotations` list of its record.
- **Method:** Offline version of the step `add_annotation_links` in genomepy's UCSC provider. Upstream queries `information_schema` of the UCSC MySQL server for the tables `ensGene`, `ncbiRefSeq`, `knownGene` and `refGene`; here the `(database, table)` rows of that result are passed in, and every table name is appended to the `annotations` list of the genome with that database name, in row order. Rows for genomes that are not in the catalogue are ignored. A matching record without an `annotations` list raises `ConfigurationError` (upstream: `KeyError`). The input is copied, not changed.
- **API:** `genomepy_ucsc_annotation_links(genomes: 'Mapping[str, Mapping[str, Any]]', table_rows: 'Sequence[Sequence[Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`.
- **Requirements:** Needs no optional package and no internet access; works after plain pip install.
- **Input:** `genomes`: mapping of UCSC database name to a record with an `annotations` list; `table_rows`: rows `(database, table)` with exactly two fields.
- **Output:** A new catalogue (dictionary) whose `annotations` lists include the table names.
- **Example code:**

```python
from dnakit.databases import genomepy_ucsc_annotation_links

catalogue = {"hg38": {"annotations": []}, "mm39": {"annotations": ["existing"]}}
rows = [("hg38", "ensGene"), ("hg38", "refGene"), ("mm39", "knownGene"), ("other", "ensGene")]
result = genomepy_ucsc_annotation_links(catalogue, rows, verbose=False)
print(result)
print(catalogue)
```

- **Example output:**

```text
{'hg38': {'annotations': ['ensGene', 'refGene']}, 'mm39': {'annotations': ['existing', 'knownGene']}}
{'hg38': {'annotations': []}, 'mm39': {'annotations': ['existing']}}
```

## 27) `dnakit.databases.genomepy_ensembl_add_grch37()`

- **Purpose:** Add the fixed GRCh37 record that genomepy lists next to the Ensembl release catalogue.
- **Method:** Offline version of `add_grch37` in genomepy's Ensembl provider. The record is a fixed literal (Homo sapiens, assembly accession `GCA_000001405.14`, taxonomy ID 9606, division `vertebrates`, Ensembl GRCh37 release 87 GTF) and is inserted under the key `GRCh37`; an existing `GRCh37` entry is replaced. The `genome` and `annotation` URLs hold a `{}` placeholder for the release number, exactly as upstream. The input is copied, not changed.
- **API:** `genomepy_ensembl_add_grch37(genomes: 'Mapping[str, Mapping[str, Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`.
- **Requirements:** Needs no optional package and no internet access; works after plain pip install.
- **Input:** `genomes`: Ensembl catalogue, a mapping of assembly name to a record (mapping).
- **Output:** A new catalogue (dictionary) with the extra key `GRCh37`; its record has `genome`, `annotation`, `assembly_accession`, `taxonomy_id`, `name`, `scientific_name`, `url_name`, `assembly_name`, `division`, `base_count`, `display_name` and `genebuild`.
- **Example code:**

```python
from dnakit.databases import genomepy_ensembl_add_grch37

catalogue = {"GRCh38": {"name": "Homo_sapiens"}}
result = genomepy_ensembl_add_grch37(catalogue, verbose=False)
record = result["GRCh37"]
print(list(result), list(catalogue))
print(record["assembly_accession"], record["taxonomy_id"], record["scientific_name"])
print(record["annotation"])
```

- **Example output:**

```text
['GRCh38', 'GRCh37'] ['GRCh38']
GCA_000001405.14 9606 Homo sapiens
http://ftp.ensembl.org/pub/grch37/release-{}/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
```

## 28) `dnakit.databases.genomepy_gencode_add_grch37()`

- **Purpose:** Add the GRCh37 annotation record that GENCODE files in its own mapping folder.
- **Method:** Offline version of `add_grch37` in genomepy's GENCODE provider. Among the keys that look like `GRCh<digits>` it takes the highest, reads the release number from the `release_` segment of that assembly's first annotation link (the last two characters of the segment) and builds the GRCh37 lift-over URL `<ftp_link>/Gencode_human/release_<N>/GRCh37_mapping/gencode.v<N>lift37.annotation.gtf.gz`. Unlike upstream, which takes the maximum over every key, only `GRCh` keys are considered; a catalogue without a usable `GRCh` key, annotation link or `release_` segment raises `ConfigurationError`. The input is copied, not changed.
- **API:** `genomepy_gencode_add_grch37(genomes: 'Mapping[str, Mapping[str, Any]]', ftp_link: 'str', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`.
- **Requirements:** Needs no optional package and no internet access; works after plain pip install.
- **Input:** `genomes`: GENCODE catalogue, a mapping of assembly name to a record with an `annotations` list; `ftp_link`: non-empty base link of the GENCODE listing, for example `ftp://ftp.ebi.ac.uk/pub/databases/gencode`.
- **Output:** A new catalogue (dictionary) with the extra key `GRCh37`; its record has `annotations` (one URL), `taxonomy_id` (9606), `species` and `text_search`.
- **Example code:**

```python
from dnakit.databases import genomepy_gencode_add_grch37

base = "ftp://ftp.ebi.ac.uk/pub/databases/gencode"
catalogue = {
    "GRCh38": {"annotations": [base + "/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"]},
    "GRCm39": {"annotations": []},
}
result = genomepy_gencode_add_grch37(catalogue, base, verbose=False)
record = result["GRCh37"]
print(list(result), list(catalogue))
print(record["annotations"][0])
print(record["taxonomy_id"], record["species"], record["text_search"])
```

- **Example output:**

```text
['GRCh38', 'GRCm39', 'GRCh37'] ['GRCh38', 'GRCm39']
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz
9606 Homo sapiens human
```

## 29) `dnakit.databases.gget_alphainfluenza_download()`

- **Purpose:** Download influenza A (alphainfluenza) sequences from NCBI's cached data packages into a folder.
- **Method:** Creates `outdir` if needed, then calls `gget.gget_virus.download_alphainfluenza_optimized` with the optional filters (host, complete-only, annotated, accession and whether to use the accession). gget returns the ZIP path, the filters NCBI could apply and the filters it could not.
- **API:** `gget_alphainfluenza_download(*, outdir: 'str', host: 'str | None' = None, complete_only: 'bool | None' = None, annotated: 'bool | None' = None, accession: 'str | None' = None, use_accession: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes an empty folder path, which is rejected before gget is imported or any request is made.
- **Input:** `outdir`: output folder path (created if missing); `host`, `accession`: optional strings; `complete_only`, `annotated`: optional booleans (`None` leaves the filter unset); `use_accession`: boolean (default False). A non-boolean flag raises `ConfigurationError`.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `outdir`, `zip_path`, `applied_filters` and `missing_filters`.
- **Example code:**

```python
from dnakit.databases import gget_alphainfluenza_download
from dnakit.exceptions import ConfigurationError

try:
    gget_alphainfluenza_download(outdir="", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError outdir must be a folder path.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_alphainfluenza_download

result = gget_alphainfluenza_download(outdir="influenza_a", host="human", complete_only=True)
print(result["zip_path"], result["applied_filters"])
```

## 30) `dnakit.databases.gget_cbio_download()`

- **Purpose:** Download cBioPortal study archives (mutations, copy number, structural variants, clinical data) into a folder.
- **Method:** Validates the study IDs, creates `out_dir` if needed, then calls `gget.gget_cbio.download_cbioportal_data(study_ids, verbose=False, out_dir, confirm_download)`. When gget reports that not all files were downloaded, the function raises `BackendExecutionError`. The `verbose` argument of DNAKit is not passed to gget.
- **API:** `gget_cbio_download(study_ids: 'Sequence[str]', *, out_dir: 'str', confirm_download: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes a study ID with spaces, which fails the identifier check before any folder is created.
- **Input:** `study_ids`: list of 1 to 10000 identifiers (letters, digits and `._:/+-`, at most 64 characters each); `out_dir`: output folder (created if missing); `confirm_download`: passed to gget (default False).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `study_ids` and `out_dir`.
- **Example code:**

```python
from dnakit.databases import gget_cbio_download
from dnakit.exceptions import ConfigurationError

try:
    gget_cbio_download(["bad study id"], out_dir="cbio", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError study_ids[] has an unsupported format: 'bad study id'.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_cbio_download

result = gget_cbio_download(["msk_impact_2017"], out_dir="cbio_data")
print(result["out_dir"])
```

## 31) `dnakit.databases.gget_cbio_plot()`

- **Purpose:** Download cBioPortal study data and draw a heatmap of gene alterations per stratum.
- **Method:** Validates the study IDs, up to 200 gene symbols, `stratification`, `variation_type`, a plain `figure_filename` (no path separators) and `dpi` (30 to 1200); creates `figure_dir` and `data_dir` if needed; then calls `gget.cbio_plot(..., show=False)`. When gget returns a false value the function raises `BackendExecutionError`. The figure path is returned only if the PNG exists afterwards.
- **API:** `gget_cbio_plot(study_ids: 'Sequence[str]', genes: 'Sequence[str]', *, figure_dir: 'str', figure_filename: 'str', data_dir: 'str', stratification: 'str' = 'tissue', variation_type: 'str' = 'mutation_occurrences', confirm_download: 'bool' = False, dpi: 'int' = 100, figure_title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes an empty study list, which is rejected before any folder is created.
- **Input:** `study_ids`: list of study identifiers; `genes`: list of 1 to 200 gene symbols; `figure_dir`, `figure_filename`, `data_dir`: explicit output locations; `stratification` (default `tissue`) and `variation_type` (default `mutation_occurrences`): passed to gget; `confirm_download`, `dpi` (default 100), `figure_title`.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `study_ids`, `genes`, `figure` (path of the PNG, or `None`) and `data_dir`.
- **Example code:**

```python
from dnakit.databases import gget_cbio_plot
from dnakit.exceptions import ConfigurationError

try:
    gget_cbio_plot(
        [], ["TP53"], figure_dir="figs", figure_filename="a.png", data_dir="data", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError study_ids must hold between 1 and 10000 identifiers.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_cbio_plot

result = gget_cbio_plot(
    ["msk_impact_2017"],
    ["TP53", "KRAS"],
    figure_dir="figures",
    figure_filename="alterations.png",
    data_dir="cbio_data",
)
print(result["figure"])
```

## 32) `dnakit.databases.gget_diamond_align()`

- **Purpose:** Align query sequences to a reference with DIAMOND and return the tabular hits.
- **Method:** Checks the query and reference files, the new output path, the sensitivity name and the thread count, then finds the `diamond` program (from `diamond_path` or PATH; a missing program raises `BackendUnavailableError` with an install hint) and reads its version. In a temporary folder next to `output_path` it runs `diamond makedb` on the reference, then `diamond blastp` (or `blastx` when `translated` is True) with `--outfmt 6` and the fixed 14 columns, and removes the folder. The output file is parsed into one dictionary per hit.
- **API:** `gget_diamond_align(query: 'str', reference: 'str', output_path: 'str', *, translated: 'bool' = False, sensitivity: 'str' = 'very-sensitive', threads: 'int' = 1, diamond_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program DIAMOND (2.1 or later) on PATH, or its path in `diamond_path` (gget ships a Linux build under `gget/bins/Linux/diamond`). The gget package itself is not imported. DIAMOND is not on PATH here, so the example shows the missing-program error. Example: The example points `diamond_path` at a file that does not exist, so it fails with `BackendUnavailableError` whether or not DIAMOND is installed.
- **Input:** `query`, `reference`: existing FASTA files; `output_path`: new file path whose folder exists; `translated`: False for protein against protein, True for nucleotide against protein (default False); `sensitivity`: `fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive` or `ultra-sensitive` (default `very-sensitive`); `threads`: 1 to 256; `diamond_path`: optional program path.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `program` (name, version, path), `mode` (`blastp` or `blastx`), `sensitivity`, `output`, `columns` (qseqid, sseqid, pident, qlen, slen, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore), `hits`, `rows` (one dictionary per hit, at most 100000) and `log_tail`.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_diamond_align
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("proteins.fa").write_text(">p1\nMKTAYIAKQRQISFVKSHFSRQ\n")
    try:
        gget_diamond_align(
            "proteins.fa", "proteins.fa", "hits.tsv",
            diamond_path="no_such_diamond", verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split("] ", 1)[-1].split(" Hint:")[0])
```

- **Example output:**

```text
BackendUnavailableError diamond_path does not exist: no_such_diamond
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_diamond_align

result = gget_diamond_align("query.fa", "reference.fa", "hits.tsv", sensitivity="sensitive")
print(result["hits"], result["rows"][0]["sseqid"])
```

## 33) `dnakit.databases.gget_genbank_metadata()`

- **Purpose:** Fetch GenBank records for a list of accessions through NCBI E-utilities and write the XML and a metadata CSV.
- **Method:** Validates the accessions and requires `xml_path`, `csv_path` and the optional `failed_log_path` to be new files in existing folders. It then calls `gget.gget_virus.fetch_genbank_metadata` with the batch size, the pause between requests and the optional NCBI API key; the key is passed to gget only and is never returned.
- **API:** `gget_genbank_metadata(accessions: 'Sequence[str]', *, xml_path: 'str', csv_path: 'str', batch_size: 'int' = 200, delay: 'float' = 0.5, failed_log_path: 'str | None' = None, api_key: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. An NCBI API key is optional and raises the request rate. Example: The example passes an empty accession list, which is rejected before any file path is checked.
- **Input:** `accessions`: list of 1 to 10000 identifiers; `xml_path`, `csv_path`: new output files; `batch_size`: 1 to 1000 (default 200); `delay`: seconds between requests, 0 to 60 (default 0.5); `failed_log_path`: optional new file; `api_key`: optional NCBI key.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `accessions` (count), `xml_path`, `csv_path` and `summary` (the value returned by gget).
- **Example code:**

```python
from dnakit.databases import gget_genbank_metadata
from dnakit.exceptions import ConfigurationError

try:
    gget_genbank_metadata([], xml_path="records.xml", csv_path="records.csv", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError accessions must hold between 1 and 10000 identifiers.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_genbank_metadata

result = gget_genbank_metadata(
    ["NC_045512.2"], xml_path="records.xml", csv_path="records.csv"
)
print(result["summary"])
```

## 34) `dnakit.databases.gget_muscle_align()`

- **Purpose:** Align nucleotide or amino-acid sequences with MUSCLE 5 and write an aligned FASTA file.
- **Method:** Checks that `fasta` exists and `output_path` is new, finds the `muscle` program (from `muscle_path` or PATH; a missing program raises `BackendUnavailableError` with an install hint), runs `muscle -version` and requires version 5 (otherwise `BackendVersionError`). It then runs `muscle -align <fasta> -output <output_path>` (or `-super5` when `super5` is True), the same command as gget's, and summarises the aligned FASTA (number of records, aligned length, identifiers).
- **API:** `gget_muscle_align(fasta: 'str', output_path: 'str', *, super5: 'bool' = False, muscle_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program MUSCLE 5 on PATH, or its path in `muscle_path` (gget ships a Linux build under `gget/bins/Linux/muscle`). The gget package itself is not imported. MUSCLE is not on PATH here, so the example shows the missing-program error. Example: The example points `muscle_path` at a file that does not exist, so it fails with `BackendUnavailableError` whether or not MUSCLE is installed.
- **Input:** `fasta`: existing FASTA file; `output_path`: new file path whose folder exists; `super5`: use the Super5 algorithm for large inputs (default False); `muscle_path`: optional program path.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `program` (name, version, path), `command`, `output`, `records`, `aligned_length` (`None` when the sequences differ in length), `ids` and `log_tail`.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_muscle_align
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("seqs.fa").write_text(">a\nACGTACGT\n>b\nACGTACGA\n")
    try:
        gget_muscle_align(
            "seqs.fa", "aligned.fa", muscle_path="no_such_muscle", verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split("] ", 1)[-1].split(" Hint:")[0])
```

- **Example output:**

```text
BackendUnavailableError muscle_path does not exist: no_such_muscle
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_muscle_align

result = gget_muscle_align("seqs.fa", "aligned.fa")
print(result["records"], result["aligned_length"], result["ids"])
```

## 35) `dnakit.databases.gget_ref()`

- **Purpose:** Return the Ensembl FTP links for a species' reference files, without downloading them.
- **Method:** Validates the species (lowercase letters, digits and underscores) and `which`, then calls `gget.ref(species, which, release, ftp, save=False)`. gget looks up the links on the Ensembl FTP site and returns them; nothing is saved. With `ftp=True` gget returns only the links.
- **API:** `gget_ref(species: 'str', *, which: 'Sequence[str] | str' = 'all', release: 'int | None' = None, ftp: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes a species name with capital letters and a space, which fails the format check before gget is imported.
- **Input:** `species`: `genus_species` such as `homo_sapiens`, or a shortcut such as `human`; `which`: `all` or a subset of `gtf`, `cdna`, `dna`, `cds`, `cdrna`, `pep` (string or list); `release`: optional Ensembl release number, 1 to 1000 (default: latest); `ftp`: return only the links (default False).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `species`, `which` (list), `release`, `ftp_only` and `links` (the value returned by gget).
- **Example code:**

```python
from dnakit.databases import gget_ref
from dnakit.exceptions import ConfigurationError

try:
    gget_ref("Homo Sapiens", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError species has an unsupported format: 'Homo Sapiens'.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_ref

result = gget_ref("homo_sapiens", which=["gtf", "dna"], ftp=True)
print(result["links"])
```

## 36) `dnakit.databases.gget_sars_cov2_download()`

- **Purpose:** Download SARS-CoV-2 sequences from NCBI's cached data packages into a folder.
- **Method:** Creates `outdir` if needed, then calls `gget.gget_virus.download_sars_cov2_optimized` with the optional filters (host, complete-only, annotated, lineage, accession and whether to use the accession). gget returns the ZIP path, the filters NCBI could apply and the filters it could not.
- **API:** `gget_sars_cov2_download(*, outdir: 'str', host: 'str | None' = None, complete_only: 'bool | None' = None, annotated: 'bool | None' = None, lineage: 'str | None' = None, accession: 'str | None' = None, use_accession: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes an empty folder path, which is rejected before gget is imported or any request is made.
- **Input:** `outdir`: output folder path (created if missing); `host`, `lineage`, `accession`: optional strings; `complete_only`, `annotated`: optional booleans; `use_accession`: boolean (default False).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `outdir`, `zip_path`, `applied_filters` and `missing_filters`.
- **Example code:**

```python
from dnakit.databases import gget_sars_cov2_download
from dnakit.exceptions import ConfigurationError

try:
    gget_sars_cov2_download(outdir="", lineage="B.1.1.7", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError outdir must be a folder path.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_sars_cov2_download

result = gget_sars_cov2_download(outdir="sars_cov2", lineage="B.1.1.7", complete_only=True)
print(result["zip_path"], result["missing_filters"])
```

## 37) `dnakit.databases.gget_virus()`

- **Purpose:** Download a virus genome dataset (sequences and metadata) from NCBI Virus into a folder.
- **Method:** Validates the taxon name (unless `is_accession` is True), the sequence-length bounds and the dates, creates `outfolder`, and calls `gget.virus(...)` with the filters unchanged; unset filters stay `None`. The value returned by gget is converted to plain Python data and returned as `result`.
- **API:** `gget_virus(virus: 'str', *, outfolder: 'str', is_accession: 'bool' = False, host: 'str | None' = None, min_seq_length: 'int | None' = None, max_seq_length: 'int | None' = None, nuc_completeness: 'str | None' = None, annotated: 'bool | None' = None, geographic_location: 'str | None' = None, min_release_date: 'str | None' = None, max_release_date: 'str | None' = None, lineage: 'str | None' = None, genbank_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes a taxon name with characters that are not allowed, which fails before any folder is created.
- **Input:** `virus`: taxon name such as `Zika virus` (letters, digits and common punctuation, at most 128 characters), or an accession list with `is_accession=True`; `outfolder`: output folder; filters `host`, `nuc_completeness`, `annotated`, `geographic_location`, `lineage`; `min_seq_length` and `max_seq_length` (non-negative integers, minimum not above maximum); `min_release_date` and `max_release_date` (`YYYY-MM-DD`); `genbank_metadata`: also fetch GenBank metadata.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `virus`, `outfolder` and `result` (the value returned by gget).
- **Example code:**

```python
from dnakit.databases import gget_virus
from dnakit.exceptions import ConfigurationError

try:
    gget_virus("<bad>", outfolder="virus_data", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError virus must be a taxon name (letters, digits and common punctuation).
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_virus

result = gget_virus(
    "Zika virus",
    outfolder="virus_data",
    host="human",
    min_seq_length=10000,
    min_release_date="2020-01-01",
)
print(result["outfolder"])
```

## 38) `dnakit.databases.gget_virus_metadata()`

- **Purpose:** Fetch virus genome metadata from the NCBI Datasets API without downloading sequences.
- **Method:** Validates the name and the optional date, then calls `gget.gget_virus.fetch_virus_metadata` with the host, location, annotation, completeness and RefSeq filters. The returned records are converted to plain data; `count` is their number when gget returns a list.
- **API:** `gget_virus_metadata(virus: 'str', *, accession: 'bool' = False, host: 'str | None' = None, geographic_location: 'str | None' = None, annotated: 'bool | None' = None, complete_only: 'bool | None' = None, min_release_date: 'str | None' = None, refseq_only: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes an empty virus name, which is rejected before gget is imported or any request is made.
- **Input:** `virus`: non-empty name (at most 256 characters); `accession`: treat `virus` as accession(s) (boolean, default False); `host`, `geographic_location`: optional strings; `annotated`, `complete_only`, `refseq_only`: optional booleans; `min_release_date`: optional `YYYY-MM-DD`.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `virus`, `records` and `count`.
- **Example code:**

```python
from dnakit.databases import gget_virus_metadata
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_metadata("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError virus must be a non-empty string of at most 256 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_virus_metadata

result = gget_virus_metadata("Zika virus", host="human", complete_only=True)
print(result["count"])
```

## 39) `dnakit.databases.gget_virus_metadata_chunked()`

- **Purpose:** Fetch metadata for a very large virus dataset in release-date chunks.
- **Method:** Validates the name and the dates, then calls `gget.gget_virus.fetch_virus_metadata_chunked`, which splits the query into release-date ranges to keep each request small. The returned records are converted to plain data; `count` is their number when gget returns a list.
- **API:** `gget_virus_metadata_chunked(virus: 'str', *, accession: 'bool' = False, host: 'str | None' = None, geographic_location: 'str | None' = None, annotated: 'bool | None' = None, complete_only: 'bool' = False, min_release_date: 'str | None' = None, max_release_date: 'str | None' = None, refseq_only: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. Example: The example passes a number instead of a name, which is rejected before gget is imported or any request is made.
- **Input:** `virus`: non-empty name (at most 256 characters); `accession`: boolean (default False); `host`, `geographic_location`: optional strings; `annotated`: optional boolean; `complete_only`, `refseq_only`: booleans (default False; `None` is refused); `min_release_date`, `max_release_date`: optional `YYYY-MM-DD`.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `virus`, `records` and `count`.
- **Example code:**

```python
from dnakit.databases import gget_virus_metadata_chunked
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_metadata_chunked(123, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError virus must be a non-empty string of at most 256 characters.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_virus_metadata_chunked

result = gget_virus_metadata_chunked(
    "Severe acute respiratory syndrome coronavirus 2",
    min_release_date="2023-01-01",
    max_release_date="2023-03-31",
)
print(result["count"])
```

## 40) `dnakit.databases.gget_virus_sequences_by_accession()`

- **Purpose:** Download virus genome sequences for a list of accession numbers from NCBI, in batches.
- **Method:** Validates the accessions, creates `outdir`, and calls `gget.gget_virus.download_sequences_by_accessions` with the batch size and the optional NCBI API key (passed to gget only, never returned). gget returns the path of the FASTA file it wrote.
- **API:** `gget_virus_sequences_by_accession(accessions: 'Sequence[str]', *, outdir: 'str', batch_size: 'int' = 200, api_key: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gget (extra `dnakit[bio-gget]`) and internet access. An NCBI API key is optional and raises the request rate. Example: The example includes an accession with a space, which fails the identifier check before any folder is created.
- **Input:** `accessions`: list of 1 to 10000 identifiers (letters, digits and `._:/+-`, at most 64 characters each); `outdir`: output folder; `batch_size`: 1 to 1000 (default 200); `api_key`: optional NCBI key.
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `accessions` (count), `outdir` and `fasta_path`.
- **Example code:**

```python
from dnakit.databases import gget_virus_sequences_by_accession
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_sequences_by_accession(
        ["NC_045512.2", "bad accession"], outdir="virus_fasta", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError accessions[] has an unsupported format: 'bad accession'.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_virus_sequences_by_accession

result = gget_virus_sequences_by_accession(["NC_045512.2"], outdir="virus_fasta")
print(result["fasta_path"])
```

## 41) `dnakit.databases.gget_cellxgene_query()`

- **Purpose:** Query CELLxGENE Census expression data for genes and cell-metadata filters and write the result to a new file.
- **Method:** Checks the arguments first: `output_path` must be new and end in `.h5ad` (or `.csv` when `meta_only` is True), the species must be one of five, and genes, columns and the five label filters must match the allowed patterns. It then checks that `cellxgene_census` can be imported, and calls `gget.gget_cellxgene.cellxgene(...)` with `out=output_path`. gget writes the expression AnnData (`.h5ad`) or, with `meta_only`, the cell-metadata table (`.csv`); no data returned by gget raises `BackendExecutionError`.
- **API:** `gget_cellxgene_query(output_path: 'str', *, species: 'str' = 'homo_sapiens', genes: 'list[str] | None' = None, ensembl: 'bool' = False, tissue: 'list[str] | None' = None, cell_type: 'list[str] | None' = None, disease: 'list[str] | None' = None, sex: 'list[str] | None' = None, assay: 'list[str] | None' = None, column_names: 'list[str] | None' = None, is_primary_data: 'bool' = True, meta_only: 'bool' = False, census_version: 'str' = 'stable', verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages gget (extra `dnakit[bio-gget]`) and cellxgene-census (extra `dnakit[bio-cellxgene]`), and internet access to the Census. cellxgene-census is not installed here, so the example shows an input check that happens before the package is imported. Example: The example asks for an expression query (`meta_only` False) but names a `.csv` output; the extension check fails before the Census package is imported, so the output does not depend on whether it is installed.
- **Input:** `output_path`: new file path; `species`: `homo_sapiens`, `mus_musculus`, `macaca_mulatta`, `callithrix_jacchus` or `pan_troglodytes`; `genes`: up to 500 symbols (Ensembl IDs when `ensembl` is True); `tissue`, `cell_type`, `disease`, `sex`, `assay`: lists of up to 50 Census labels; `column_names`; `is_primary_data` (default True); `meta_only` (default False); `census_version` (default `stable`).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `species`, `census_version`, `genes`, `ensembl`, `filters` (only the filters that were set), `is_primary_data`, `meta_only`, `output_path`, `counts` (`rows` for metadata, or `cells` and `genes`), `obs_columns` and `var_names` (each at most 200 entries).
- **Example code:**

```python
from dnakit.databases import gget_cellxgene_query
from dnakit.exceptions import ConfigurationError

try:
    gget_cellxgene_query("expression.csv", genes=["TP53"], tissue=["lung"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError output_path must end in .h5ad for this query.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_cellxgene_query

result = gget_cellxgene_query(
    "lung_tp53.h5ad", genes=["TP53"], tissue=["lung"], cell_type=["macrophage"]
)
print(result["counts"])
```

## 42) `dnakit.databases.gget_cosmic_download_reference()`

- **Purpose:** Download one COSMIC reference archive with the user's own account and extract it.
- **Method:** Checks that the link is an https address on `cancer.sanger.ac.uk` whose `path` parameter names a `.tar` file, and checks the e-mail and password. The first request sends the account as HTTP basic authentication and receives a signed download address; the archive is then fetched from that address without the credentials, saved as `<output_folder>/<name>.tar`, and extracted into `<output_folder>/<name>/` after a member check that rejects absolute paths, `..` and links. With `unzip` True each extracted `.gz` file is also decompressed next to itself. The password is never stored or returned. Only the Python standard library is used for the transfer; gget is not imported.
- **API:** `gget_cosmic_download_reference(download_link: 'str', output_folder: 'str', *, email: 'str', password: 'str', unzip: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Works after plain pip install (gget is not needed), but needs internet access and a COSMIC account (e-mail and password) that is allowed to download the file. Example: The example uses a link on another host; it is refused before any request is made, so the example sends no credentials.
- **Input:** `download_link`: https link on `cancer.sanger.ac.uk` with a `path=` parameter that ends in `.tar`; `output_folder`: target folder (created if missing); `email`, `password`: account (keyword-only); `unzip`: decompress `.gz` files (default False).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `download_link`, `tar_path`, `extracted_folder`, `files` (at most 200 extracted paths) and `unzipped` (paths of decompressed files).
- **Example code:**

```python
from dnakit.databases import gget_cosmic_download_reference
from dnakit.exceptions import ConfigurationError

try:
    gget_cosmic_download_reference(
        "https://example.org/data.tar",
        "cosmic",
        email="user@example.org",
        password="not-a-real-password",
        verbose=False,
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **Example output:**

```text
ConfigurationError download_link must be an https address on cancer.sanger.ac.uk.
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_cosmic_download_reference

link = (
    "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"
    "?path=grch38/cosmic/v99/Cosmic_MutantCensus_Tsv_v99_GRCh38.tar&bucket=downloads"
)
result = gget_cosmic_download_reference(
    link, "cosmic", email="you@example.org", password="your-cosmic-password", unzip=True
)
print(result["extracted_folder"])
```

## 43) `dnakit.databases.gget_cosmic_select_reference()`

- **Purpose:** Download one COSMIC project release for a genome build, or reuse the file if it already exists.
- **Method:** Checks the project name (`cancer`, `cell_line`, `census`, `resistance`, `genome_screen` or `targeted_screen`), the GRCh version (37 or 38; `cancer` only exists for 37) and the release number (1 to 200), builds the archive name and the expected file name, and looks in `output_folder`. If `<output_folder>/<archive name>/<file>` exists (the `.tsv`, or `.tsv.gz` when `unzip` is False) and `overwrite` is False, it returns at once with `downloaded` False. Otherwise it calls `gget_cosmic_download_reference` with the account and returns the file path.
- **API:** `gget_cosmic_select_reference(output_folder: 'str', cosmic_project: 'str', grch_version: 'int', cosmic_version: 'int', *, email: 'str', password: 'str', unzip: 'bool' = True, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Works after plain pip install (gget is not needed). Reusing an existing file needs no network; a download needs internet access and a COSMIC account (e-mail and password). Example: The example creates the file that the function expects for the `census` project, GRCh38, release 99, so the function reuses it and sends no request. The account arguments are placeholders and are not returned.
- **Input:** `output_folder`: folder for the files; `cosmic_project`: project name; `grch_version`: 37 or 38; `cosmic_version`: release number; `email`, `password`: account (keyword-only; not used when the existing file is reused); `unzip`: default True; `overwrite`: download again even when the file exists (default False).
- **Output:** A dictionary with the common gget result fields `operation`, `status` (`completed`), `backend` (`gget`), `package` (`gget <version>`), `service` (name of the service that was used) and `accessed` (today's date), `cosmic_project`, `grch_version`, `cosmic_version`, `file_path` and `downloaded`. When reused it also has `overwritten` (False); after a download it has `download_link`, `overwritten`, `tar_path` and `file_count` instead.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_cosmic_select_reference

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    release = Path("cosmic") / "Cosmic_MutantCensus_Tsv_v99_GRCh38"
    release.mkdir(parents=True)
    (release / "Cosmic_MutantCensus_v99_GRCh38.tsv").write_text(
        "GENE_SYMBOL\tMUTATION_AA\nTP53\tp.R175H\n"
    )
    result = gget_cosmic_select_reference(
        "cosmic", "census", 38, 99,
        email="user@example.org", password="not-a-real-password", verbose=False,
    )
    print(result["operation"], result["cosmic_project"], result["grch_version"], result["cosmic_version"])
    print(Path(result["file_path"]).name, result["downloaded"], result["overwritten"])
    print("password" in str(result))
```

- **Example output:**

```text
select_reference census 38 99
Cosmic_MutantCensus_v99_GRCh38.tsv False False
False
```

- **Usage with the dependency installed:**

```python
from dnakit.databases import gget_cosmic_select_reference

result = gget_cosmic_select_reference(
    "cosmic", "census", 38, 99, email="you@example.org", password="your-cosmic-password"
)
print(result["file_path"], result["downloaded"])
```
