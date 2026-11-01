# Variant and annotation I/O: bed-reader, bgen-reader, hgvs, polars-bio and pyensembl

This group of functions in `dnakit.variants` and `dnakit.annotation` is a set of adapters over five optional packages: [bed-reader](https://github.com/fastlmm/bed-reader) (PLINK 1 `.bed` files, also at a URL), [bgen-reader](https://github.com/limix/bgen-reader-py) (BGEN sample files), [hgvs](https://github.com/biocommons/hgvs) (HGVS variant strings and sequence-context text), [polars-bio](https://github.com/biodatageeks/polars-bio) (BGEN, PLINK 2 PGEN, VCF Zarr and CRAM tables, and an interval plot) and [pyensembl](https://github.com/openvax/pyensembl) (Ensembl release caches and annotation indexes). They read genotype files into plain dictionaries, format HGVS text, map one variant to its g., c., n. and p. forms, list Ensembl annotation dates, and build or inspect pyensembl annotation caches. They do not call variants or predict their effect.

The functions are imported with `from dnakit.variants import ...` (sections 1 to 10) and `from dnakit.annotation import ...` (sections 11 to 33). bed-reader, hgvs, polars-bio and pyensembl are Apache-2.0 and bgen-reader is MIT, as stated in the installed package metadata; none of them is copyleft. The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises `BackendUnavailableError` with an install hint.

Several functions need something that is not a plain Python package: the two bed-reader sample functions and the three pyensembl listing functions use the network, `pyensembl_install_ensembl_release()` downloads from the Ensembl FTP server, and the two HGVS mapper functions need an hgvs mapper built on a UTA database. For these the example shows the function's real behaviour on invalid input, and an unexecuted block shows the call with the dependency in place.

## Example data

The BGEN and PGEN files are binary, so they are declared as hex text and written to disk by the examples. `cohort.bgen.hex` is a 150-byte BGEN v1.2 file (layout 2, 8-bit probabilities, no compression) with three samples (S1, S2, S3) and two variants (`rs1` at chromosome 1 position 100 and `rs2` at position 250). `cohort.pgen.hex` is a PLINK 1-mode PGEN with the same three samples and two variants: the three magic bytes `6c1b01` and one packed byte per variant. `cohort.pvar` and `cohort.psam` are its companion files. `make_vcz.py` writes a two-variant, two-sample VCF Zarr store with the `zarr` package. `a.bed` and `b.bed` are two small BED files with overlapping intervals, and `mini.gtf` is a three-exon GTF annotation.

**Example file:** `cohort.bgen.hex`

```text
280000001400000002000000030000006267656e08000080140000000300
000002005331020053320200533304007661723003007273310100316400
000002000100000041010000004713000000030000000200020202020200
08ff0000ff00000400766172310300727332010031fa0000000200010000
0043010000005413000000030000000200020202020200080000ff0000ff
```

**Example file:** `cohort.pgen.hex`

```text
6c1b010b1e
```

**Example file:** `cohort.pvar`

```text
#CHROM	POS	ID	REF	ALT
1	100	rs1	A	G
1	250	rs2	C	T
```

**Example file:** `cohort.psam`

```text
#IID	SEX
S1	NA
S2	NA
S3	NA
```

**Example file:** `make_vcz.py`

```text
import numpy as np
import zarr
from zarr.core.dtype import VariableLengthUTF8


def make_vcz(path):
    """Write a 2-variant, 2-sample VCF Zarr store (Zarr v2 layout)."""
    root = zarr.open_group(path, mode="w", zarr_format=2)
    root.attrs["vcf_zarr_version"] = "0.4"
    root.attrs["vcf_header"] = (
        "##fileformat=VCFv4.3\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    )
    root.attrs["contigs"] = [{"ID": "chr1", "length": 1000}]

    def put(name, data, dtype=None):
        data = np.asarray(data, dtype=dtype)
        zdtype = VariableLengthUTF8() if data.dtype == object else data.dtype
        array = root.create_array(name, shape=data.shape, dtype=zdtype, chunks=data.shape)
        array[...] = data

    put("contig_id", ["chr1"], object)
    put("sample_id", ["S1", "S2"], object)
    put("variant_contig", [0, 0], np.int32)
    put("variant_position", [100, 250], np.int32)
    put("variant_id", ["rs1", "rs2"], object)
    put("variant_allele", [["A", "G"], ["C", "T"]], object)
    put("variant_quality", [50.0, 60.0], np.float32)
    put("variant_filter", np.zeros((2, 1), dtype=bool))
    put("filter_id", ["PASS"], object)
    put("variant_DP", [14, 20], np.int32)
    genotype = np.array([[[0, 0], [0, 1]], [[1, 1], [0, 0]]], dtype=np.int8)
    put("call_genotype", genotype)
    put("call_genotype_mask", genotype < 0)
    put("call_genotype_phased", np.zeros((2, 2), dtype=bool))
```

**Example file:** `a.bed`

```text
chr1	10	20	a
chr1	40	60	b
```

**Example file:** `b.bed`

```text
chr1	15	30	c
chr1	50	55	d
```

**Example file:** `mini.gtf`

```text
chr1	toy	exon	100	300	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number "1"; gene_name "GENEA";
chr1	toy	exon	600	900	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number "2"; gene_name "GENEA";
chr1	toy	exon	2000	2500	.	-	.	gene_id "G2"; transcript_id "T2"; exon_number "1"; gene_name "GENEB";
```

## 1) `dnakit.variants.check_cloud_location()`

- **Purpose:** Report whether bed-reader can read the object at a URL.
- **Method:** Validate the URL, then call bed-reader's native `check_file_cloud(url, options)`. A `ValueError` from the native check (a missing or unreadable object) is reported as `False`; any other error is raised.
- **API:** `check_cloud_location(location: 'str', *, cloud_options: 'dict[str, str] | None' = None, verbose: bool = True) -> 'bool'`.
- **Requirements:** Needs the optional package bed-reader (extra `dnakit[bio-bed-reader]`). A `file://` URL reads a local file; `https://` and `s3://` URLs need network access and, for private buckets, credentials in `cloud_options`.
- **Input:** `location` must be a string of at most 4096 characters with a scheme (`https://`, `s3://`, `file://`); a plain local path is refused with `ConfigurationError`. `cloud_options` is an optional dictionary of strings passed to the storage client.
- **Output:** A bool, `True` when the object can be opened for reading.
- **Example code:**

```python
from pathlib import Path

import numpy as np
from bed_reader import to_bed

from dnakit.variants import check_cloud_location

to_bed("tiny.bed", np.array([[1, 0, np.nan], [2, 0, 2]], dtype=np.float32))
url = Path("tiny.bed").resolve().as_uri()
print(check_cloud_location(url, verbose=False))
print(check_cloud_location(url.replace("tiny", "absent"), verbose=False))
```

- **Example output:**

```text
True
False
```

## 2) `dnakit.variants.read_cloud_genotypes()`

- **Purpose:** Read the genotype matrix of a PLINK 1 `.bed` file given by URL.
- **Method:** Open the URL with bed-reader's `open_bed(url, cloud_options=...)` and call `read(index=(iid, sid), dtype=dtype)`, which uses the native `read_cloud_i8`, `read_cloud_f32` or `read_cloud_f64` reader. Missing calls (-127 in `int8`, NaN in floats) become `None`. The selected cell count is checked against `max_cells` before the read.
- **API:** `read_cloud_genotypes(location: 'str', *, dtype: 'str' = 'float32', iid_index: 'list[int] | None' = None, sid_index: 'list[int] | None' = None, cloud_options: 'dict[str, str] | None' = None, max_cells: 'int' = 5000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages bed-reader and NumPy (extra `dnakit[bio-bed-reader]`). The `.fam` and `.bim` files must sit beside the `.bed` object. Remote URLs need network access.
- **Input:** `location` as above. `dtype` is `"int8"`, `"float32"` or `"float64"`. `iid_index` (samples) and `sid_index` (variants) are non-empty lists of 0-based positions in range, or `None` for all. `max_cells` is bounded above by 100,000,000; a larger selection raises `ConfigurationError`.
- **Output:** A dictionary with `status`, `operation`, `location`, `dtype`, `shape` (samples by variants), `missing_calls`, `values` (nested lists, `None` for missing) and `provenance`.
- **Example code:**

```python
from pathlib import Path

import numpy as np
from bed_reader import to_bed

from dnakit.variants import read_cloud_genotypes

to_bed("tiny.bed", np.array([[1, 0, np.nan], [2, 0, 2]], dtype=np.float32))
url = Path("tiny.bed").resolve().as_uri()
result = read_cloud_genotypes(url, verbose=False)
print(result["shape"], result["missing_calls"], result["values"])
subset = read_cloud_genotypes(url, dtype="int8", iid_index=[1], sid_index=[2, 0], verbose=False)
print(subset["dtype"], subset["shape"], subset["values"])
```

- **Example output:**

```text
[2, 3] 1 [[1.0, 0.0, None], [2.0, 0.0, 2.0]]
int8 [1, 2] [[2, 2]]
```

## 3) `dnakit.variants.sample_bed_file()`

- **Purpose:** Download one of bed-reader's sample `.bed` files (if not cached) and return its local path.
- **Method:** Validate the name, then call `bed_reader.sample_file(name)`, which fetches the file with pooch into bed-reader's cache. For a name ending in `.bed` the paths of the matching `.fam` and `.bim` files are added when they exist. A missing pooch package is reported as `BackendUnavailableError` with the `bed-reader[samples]` install hint.
- **API:** `sample_bed_file(filepath: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages bed-reader and pooch (extra `dnakit[bio-bed-reader]`) and internet access the first time a sample is fetched. The example below shows the name check, which runs before any download; with the packages and network in place, use the call in the usage block.
- **Input:** `filepath` is the name of a bed-reader sample file such as `small.bed`: a non-empty string of at most 255 characters without control characters.
- **Output:** `service`, `source` (the name) and `result` with `bed`, `fam` and `bim` (absolute paths, `fam` and `bim` may be `None`).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import sample_bed_file

try:
    sample_bed_file("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
filepath must be the name of a bed-reader sample file, such as small.bed.
```

Usage with the dependency installed:

```python
from dnakit.variants import sample_bed_file

result = sample_bed_file("small.bed")
print(result["result"]["bed"], result["result"]["fam"], result["result"]["bim"])
```

## 4) `dnakit.variants.sample_bed_url()`

- **Purpose:** Return the URL of one of bed-reader's sample `.bed` files.
- **Method:** Validate the name, then call `bed_reader.sample_url(name)`, which fetches the file with pooch first when it is not cached and returns a `file://` URL. A missing pooch package is reported as `BackendUnavailableError`.
- **API:** `sample_bed_url(filepath: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages bed-reader and pooch (extra `dnakit[bio-bed-reader]`) and internet access the first time a sample is fetched. The example shows the name check; the usage block shows the call with network access.
- **Input:** `filepath`, validated as in `sample_bed_file()`.
- **Output:** `service`, `source` and `result` with the single key `url`. The URL can be passed to `read_cloud_genotypes()`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import sample_bed_url

try:
    sample_bed_url("bad\nname.bed", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
filepath must not contain control characters.
```

Usage with the dependency installed:

```python
from dnakit.variants import read_cloud_genotypes, sample_bed_url

url = sample_bed_url("small.bed")["result"]["url"]
print(read_cloud_genotypes(url)["shape"])
```

## 5) `dnakit.variants.fetch_bgen_example()`

- **Purpose:** Return the local path of one of bgen-reader's sample BGEN files, downloading it on first use.
- **Method:** Check the name against the list bgen-reader keeps in `bgen_reader._example`, then call `bgen_reader.example_filepath(name)`, which downloads the file from the `fastlmm/bgen-sample-files` repository on GitHub into the bgen-reader cache and checks its hash. The result adds the file size and the source URL.
- **API:** `fetch_bgen_example(filename: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package bgen-reader (extra `dnakit[bio-bgen-reader]`, which also installs `requests`) and internet access the first time a file is fetched. The example shows the refusal of an unknown name, which happens before any download.
- **Input:** `filename` must be one of the names bgen-reader knows, for example `example.bgen` or `complex.bgen`; anything else raises `ConfigurationError` listing the known names.
- **Output:** `status`, `backend`, `name`, `path` (absolute local path), `size_bytes` and `source` (the download URL).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import fetch_bgen_example

try:
    fetch_bgen_example("nope.bgen", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__)
    print(error.message.split(" Known names")[0])
```

- **Example output:**

```text
ConfigurationError
Unknown bgen_reader example file 'nope.bgen'.
```

Usage with the dependency installed:

```python
from dnakit.variants import fetch_bgen_example

result = fetch_bgen_example("example.bgen")
print(result["path"], result["size_bytes"])
```

## 6) `dnakit.variants.format_sequence_blocks()`

- **Purpose:** Lay out a DNA sequence in groups of bases with a coordinate ruler.
- **Method:** Check the arguments and call hgvs' `hgvs.utils.context.format_sequence(sequence, start, end, group_size=...)`, joining the lines it yields into one string. Each row has 100 columns: a coordinate column, a ruler line that marks the last base of every group, and the bases. Displayed coordinates are 1-based.
- **API:** `format_sequence_blocks(sequence: 'str', *, start: 'int | None' = None, end: 'int | None' = None, group_size: 'int' = 3, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package hgvs (extra `dnakit[bio-hgvs]`). It uses only the hgvs formatting code, with no database.
- **Input:** `sequence` is a non-empty string. `start` and `end` are 0-based slice offsets with `0 <= start < end <= len(sequence)` (defaults: the whole sequence); the displayed numbers are 1-based. `group_size` is an integer from 1 to 88.
- **Output:** One multi-line string. Rows are separated by blank lines.
- **Example code:**

```python
from dnakit.variants import format_sequence_blocks

text = format_sequence_blocks("ACGTACGTACGTACGTAC", verbose=False)
part = format_sequence_blocks("ACGTACGTACGTACGTAC", start=3, end=12, group_size=4, verbose=False)
for block in (text, part):
    for line in block.splitlines():
        if line.strip():
            print("|" + line.rstrip())
```

- **Example output:**

```text
|          :   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60
|        1 : ACG TAC GTA CGT ACG TAC
|          :    7   11   15   19   23   27   31   35   39   43   47   51   55   59   63   67   71
|        4 : TACG TACG T
```

## 7) `dnakit.variants.hgvs_full_house()`

- **Purpose:** Express one variant as its g., c., n. and p. HGVS strings.
- **Method:** Parse the variant with the hgvs parser and call `hgvs.utils.context.full_house(mapper, variant, tx_ac=transcript)`. For a genomic variant without `transcript` it asks the mapper for `relevant_transcripts` and requires exactly one; it then calls `g_to_n`, `n_to_c` and `c_to_p` (for n. and c. input it calls `n_to_g`, `c_to_g` and `c_to_n` as needed). A `RuntimeError` from hgvs becomes `ConfigurationError`.
- **API:** `hgvs_full_house(mapper: 'Any', variant: 'str', *, transcript: 'str | None' = None, verbose: bool = True) -> 'dict[str, str]'`.
- **Requirements:** Needs the optional package hgvs (extra `dnakit[bio-hgvs]`) and an hgvs mapper built on a UTA database connection and a sequence store; hgvs ships neither, and the UTA database is reached over the network or from a local copy. The docstring calls the mapper an `AlignmentMapper`, but the methods used (`relevant_transcripts`, `g_to_n`, `n_to_c`, `c_to_p`) belong to hgvs' `AssemblyMapper`. The example uses a small stand-in object that lists two transcripts, to show the refusal; the usage block shows a real mapper.
- **Input:** `mapper` (required, not `None`), `variant` as an HGVS string such as `NC_000001.11:g.100A>G`, and `transcript`, which is needed when the mapper reports more than one overlapping transcript for a genomic variant.
- **Output:** A dictionary of four strings keyed `g`, `c`, `n` and `p`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import hgvs_full_house


class TwoTranscripts:
    """Stand-in mapper: two transcripts overlap every variant."""

    def relevant_transcripts(self, variant):
        return ["NM_000001.1", "NM_000002.1"]


try:
    hgvs_full_house(TwoTranscripts(), "NC_000001.11:g.100A>G", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError
2 relevant transcripts for NC_000001.11; you need to pick one
```

Usage with the dependency installed:

```python
import hgvs.assemblymapper
import hgvs.dataproviders.uta

from dnakit.variants import hgvs_full_house

hdp = hgvs.dataproviders.uta.connect()  # needs the UTA database and a sequence store
mapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
print(hgvs_full_house(mapper, "NC_000007.13:g.117199646_117199648del", transcript="NM_000492.3"))
```

## 8) `dnakit.variants.hgvs_pointer_line()`

- **Purpose:** Draw the marker line that points at a variant inside a displayed sequence span.
- **Method:** Parse the variant and call `hgvs.utils.context.pointer_line(variant, (start, end))`. Insertions are marked with `><`; other edits get one `*` per base. The line ends with the variant text.
- **API:** `hgvs_pointer_line(variant: 'str', span: 'Sequence[int]', *, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package hgvs (extra `dnakit[bio-hgvs]`). It uses only the hgvs parser, with no database.
- **Input:** `variant` is a genomic, coding or non-coding HGVS variant with a position (g., c. or n.; other types raise `ConfigurationError`). `span` is a pair of integers with `0 <= start <= end`.
- **Output:** One text line.
- **Example code:**

```python
from dnakit.variants import hgvs_pointer_line

for variant in ("NM_000399.3:c.5A>G", "NM_000399.3:c.5_6insA"):
    line = hgvs_pointer_line(variant, (0, 10), verbose=False)
    print("|" + line.rstrip())
```

- **Example output:**

```text
|                                    *                                         NM_000399.3:c.5A>G
|                                    ><                                        NM_000399.3:c.5_6insA
```

## 9) `dnakit.variants.hgvs_sequence_line()`

- **Purpose:** Format one sequence line of a variant-context display.
- **Method:** Parse the variant and call `hgvs.utils.context.seq_line_fmt(var=..., span=..., content=..., dir=...)`. The line shows the accession and variant type at the left, the start coordinate, the strand marker, the sequence text, the strand marker again, the end coordinate, and the variant string as a trailing comment.
- **API:** `hgvs_sequence_line(variant: 'str', span: 'Sequence[int]', content: 'str', *, direction: 'str' = '', verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package hgvs (extra `dnakit[bio-hgvs]`). It uses only the hgvs parser, with no database.
- **Input:** `variant` is any HGVS variant string that parses. `span` is a pair of integers `(start, end)` with `0 <= start <= end`. `content` is the sequence text to print. `direction` is `""`, `">"` or `"<"`.
- **Output:** One text line.
- **Example code:**

```python
from dnakit.variants import hgvs_sequence_line

line = hgvs_sequence_line("NM_000399.3:c.5A>G", (1, 10), "ACGTACGTA", direction=">", verbose=False)
print("|" + line.rstrip())
```

- **Example output:**

```text
|   NM_000399.3  c          1 > ACGTACGTA                                     >       10 NM_000399.3:c.5A>G
```

## 10) `dnakit.variants.hgvs_variant_context_with_alignment()`

- **Purpose:** Show a coding variant against the genome-to-transcript alignment, with genomic and cDNA context.
- **Method:** Check that the `uta_align` package imports, parse the variant, and call hgvs' experimental `variant_context_w_alignment(mapper, variant, margin=..., tx_ac=...)`. It fetches the genomic and transcript sequences from the mapper's data provider, aligns them globally with `uta_align`, and returns a multi-line text figure. It is limited to coding variants.
- **API:** `hgvs_variant_context_with_alignment(mapper: 'Any', variant: 'str', *, margin: 'int' = 20, transcript: 'str | None' = None, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional packages hgvs and uta-align (extra `dnakit[bio-hgvs-alignment]`) and an hgvs mapper on a UTA database connection and a sequence store. Without `uta_align` it raises `BackendUnavailableError` with an install hint, which is what the example shows on a machine without the package.
- **Input:** `mapper`, `variant` and `transcript` as in `hgvs_full_house()`. `margin` is the number of bases added on each side and must not be negative.
- **Output:** A multi-line string.
- **Example code:**

```python
from dnakit.exceptions import BackendUnavailableError
from dnakit.variants import hgvs_variant_context_with_alignment

try:
    hgvs_variant_context_with_alignment(None, "NM_000399.3:c.5A>G", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.code)
    print("uta-align" in error.hint)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
True
```

Usage with the dependency installed:

```python
import hgvs.assemblymapper
import hgvs.dataproviders.uta

from dnakit.variants import hgvs_variant_context_with_alignment

hdp = hgvs.dataproviders.uta.connect()  # needs the UTA database and a sequence store
mapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
print(hgvs_variant_context_with_alignment(mapper, "NM_000492.3:c.1521_1523del", margin=10))
```

## 11) `dnakit.annotation.polars_bio_describe_bgen()`

- **Purpose:** List the columns a local BGEN file produces, with their types and layout.
- **Method:** Call polars-bio's `describe_bgen(path, sample_path=..., bgi_path=...)` and convert each row of the result to a plain dictionary.
- **API:** `polars_bio_describe_bgen(path: 'str', *, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local file ending in `.bgen`; URLs (`://`) are refused. `sample_path` and `bgi_path` are optional local companion files.
- **Output:** `status`, `operation`, `field_count`, `fields` (one dictionary per output column with `name`, `type`, `layout`, `index` and coordinate-system fields) and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_describe_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_describe_bgen("cohort.bgen", verbose=False)
print(result["field_count"])
print([field["name"] for field in result["fields"]])
print(result["fields"][0]["layout"], result["fields"][0]["coordinate_system_zero_based"])
```

- **Example output:**

```text
9
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
2 false
```

## 12) `dnakit.annotation.polars_bio_describe_pgen()`

- **Purpose:** List the columns a local PLINK 2 PGEN fileset produces.
- **Method:** Call polars-bio's `describe_pgen(path, pvar_path=..., psam_path=..., pgi_path=...)` and convert each row to a plain dictionary.
- **API:** `polars_bio_describe_pgen(path: 'str', *, pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local file ending in `.pgen`. The `.pvar` and `.psam` files are looked up beside it unless `pvar_path`, `psam_path` or `pgi_path` give them.
- **Output:** `status`, `operation`, `field_count`, `fields` (one dictionary per output column with `name`, `type`, `storage_mode`, `index` and coordinate-system fields) and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_describe_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_describe_pgen("cohort.pgen", verbose=False)
print(result["field_count"], result["fields"][0]["storage_mode"])
print([field["name"] for field in result["fields"]])
```

- **Example output:**

```text
7 0x01
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
```

## 13) `dnakit.annotation.polars_bio_describe_vcf_zarr()`

- **Purpose:** List the INFO and FORMAT fields of a local VCF Zarr store.
- **Method:** Call polars-bio's `describe_vcf_zarr(path)` and convert each row to a plain dictionary.
- **API:** `polars_bio_describe_vcf_zarr(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`). Building the example store also needs `zarr`.
- **Input:** `path` is an existing local directory (the Zarr store); URLs are refused.
- **Output:** `status`, `operation`, `field_count`, `fields` (each with `name`, `field_type` of `INFO` or `FORMAT`, `data_type` and `description`) and `provenance`.
- **Example code:**

```python
from dnakit.annotation import polars_bio_describe_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_describe_vcf_zarr("cohort.vcz", verbose=False)
print(result["field_count"])
for field in result["fields"]:
    print(field["name"], field["field_type"], field["data_type"])
```

- **Example output:**

```text
2
DP INFO Integer
genotypes FORMAT Struct
```

## 14) `dnakit.annotation.polars_bio_read_bgen()`

- **Purpose:** Read a local BGEN file into one record per variant, with probabilities or dosages.
- **Method:** Call polars-bio's `read_bgen(path, genotype_output=..., samples=..., sample_path=..., bgi_path=...)`, convert the polars frame to dictionaries, and refuse to return more than `max_rows` rows.
- **API:** `polars_bio_read_bgen(path: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.bgen` file. `genotype_output` is `"probability"` (field `GP`) or `"dosage"` (field `DS`). `samples` is an optional non-empty list of sample names to keep. `max_rows` is an integer from 1 to 1,000,000; a file with more variants raises `ConfigurationError`.
- **Output:** `status`, `operation`, `genotype_output`, `variant_count`, `variants` (each with `chrom`, `start`, `end`, `id`, `rsid`, `alleles`, `phased`, `bits` and a `genotypes` mapping) and `provenance`. Coordinates are 1-based.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_read_bgen("cohort.bgen", verbose=False)
print(result["variant_count"], result["genotype_output"])
first = result["variants"][0]
print(first["rsid"], first["start"], first["alleles"], first["genotypes"]["GP"])
dosage = polars_bio_read_bgen("cohort.bgen", genotype_output="dosage", verbose=False)
print(dosage["variants"][1]["genotypes"]["DS"])
```

- **Example output:**

```text
2 probability
rs1 100 ['A', 'G'] [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
[2.0, 0.0, 1.0]
```

## 15) `dnakit.annotation.polars_bio_read_bgen_matrix()`

- **Purpose:** Read the ALT dosages of a local BGEN file as a dense variants-by-samples matrix.
- **Method:** Call polars-bio's `read_bgen_matrix(path, samples=..., missing=..., ...)` and convert the matrix, its positions and its sample names to lists. The matrix size is limited to 5,000,000 cells.
- **API:** `polars_bio_read_bgen_matrix(path: 'str', *, samples: 'list[str] | None' = None, missing: 'float | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.bgen` file. `samples` selects and orders samples by name. `missing` is a number used for missing dosages, or `None` for the polars-bio default.
- **Output:** `status`, `operation`, `variant_count`, `sample_count`, `positions`, `sample_names`, `values` (a list of rows, one per variant) and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_bgen_matrix

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_read_bgen_matrix("cohort.bgen", verbose=False)
print(result["variant_count"], result["sample_count"], result["positions"])
print(result["sample_names"], result["values"])
subset = polars_bio_read_bgen_matrix("cohort.bgen", samples=["S3", "S1"], verbose=False)
print(subset["sample_names"], subset["values"])
```

- **Example output:**

```text
2 3 [100, 250]
['S1', 'S2', 'S3'] [[0.0, 1.0, 2.0], [2.0, 0.0, 1.0]]
['S3', 'S1'] [[2.0, 0.0], [1.0, 2.0]]
```

## 16) `dnakit.annotation.polars_bio_read_pgen()`

- **Purpose:** Read a local PLINK 2 PGEN fileset into one record per variant.
- **Method:** Call polars-bio's `read_pgen(path, genotype_fields=..., samples=..., psam_id_mode=..., ...)`, convert the frame to dictionaries, and refuse more than `max_rows` rows. The field list defaults to `["GT"]`.
- **API:** `polars_bio_read_pgen(path: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.pgen` file with `.pvar` and `.psam` beside it. `genotype_fields` is a non-empty list of field names (default `GT`). `psam_id_mode` is `"iid"`, `"fid_iid"` or `"fid_iid_sid"`. `samples` keeps named samples. `max_rows` is an integer from 1 to 1,000,000.
- **Output:** `status`, `operation`, `variant_count`, `variants` (each with `chrom`, `start`, `end`, `id`, `ref`, `alt` and a `genotypes` mapping, where a missing call is `None`) and `provenance`. Coordinates are 1-based.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_read_pgen("cohort.pgen", verbose=False)
print(result["variant_count"])
for variant in result["variants"]:
    print(variant["id"], variant["start"], variant["ref"], variant["alt"], variant["genotypes"]["GT"])
```

- **Example output:**

```text
2
rs1 100 A ['G'] [[0, 0], [0, 1], [1, 1]]
rs2 250 C ['T'] [[0, 1], [0, 0], None]
```

## 17) `dnakit.annotation.polars_bio_read_pgen_matrix()`

- **Purpose:** Read one genotype field of a PGEN fileset as a dense variants-by-samples matrix.
- **Method:** Call polars-bio's `read_pgen_matrix(path, field=..., samples=..., missing=..., psam_id_mode=..., ...)` and convert the matrix, positions and sample names to lists. The matrix size is limited to 5,000,000 cells.
- **API:** `polars_bio_read_pgen_matrix(path: 'str', *, field: 'str' = 'ALT_COUNT', samples: 'list[str] | None' = None, missing: 'float | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.pgen` file. `field` is `"ALT_COUNT"` or `"DS"`. `missing` is a number used for missing calls, or `None` for the polars-bio default (-9 for `ALT_COUNT`). The other arguments are as in `polars_bio_read_pgen()`.
- **Output:** `status`, `operation`, `variant_count`, `sample_count`, `positions`, `sample_names`, `values` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_pgen_matrix

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
default = polars_bio_read_pgen_matrix("cohort.pgen", verbose=False)
print(default["sample_names"], default["values"])
marked = polars_bio_read_pgen_matrix("cohort.pgen", missing=-1, verbose=False)
print(marked["values"])
```

- **Example output:**

```text
['S1', 'S2', 'S3'] [[0, 1, 2], [1, 0, -9]]
[[0, 1, 2], [1, 0, -1]]
```

## 18) `dnakit.annotation.polars_bio_read_vcf_zarr()`

- **Purpose:** Read a local VCF Zarr store into one record per variant.
- **Method:** Call polars-bio's `read_vcf_zarr(path, info_fields=..., format_fields=..., samples=...)`, convert the frame to dictionaries, and refuse more than `max_rows` rows.
- **API:** `polars_bio_read_vcf_zarr(path: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`). Building the example store also needs `zarr`.
- **Input:** `path` is an existing local directory. `info_fields`, `format_fields` and `samples` are optional non-empty lists of names. `max_rows` is an integer from 1 to 1,000,000.
- **Output:** `status`, `operation`, `variant_count`, `variants` (each with `chrom`, `start`, `end`, `id`, `ref`, `alt`, `qual`, `filter`, the requested INFO fields and a `genotypes` mapping) and `provenance`. Coordinates are 1-based.
- **Example code:**

```python
from dnakit.annotation import polars_bio_read_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_read_vcf_zarr("cohort.vcz", info_fields=["DP"], verbose=False)
print(result["variant_count"])
for variant in result["variants"]:
    print(variant["id"], variant["start"], variant["ref"], variant["alt"], variant["DP"], variant["genotypes"]["GT"])
```

- **Example output:**

```text
2
rs1 100 A G 14 [[0, 0], [0, 1]]
rs2 250 C T 20 [[1, 1], [0, 0]]
```

## 19) `dnakit.annotation.polars_bio_register_bgen()`

- **Purpose:** Load a local BGEN file as a temporary SQL table and return all of its rows.
- **Method:** Register the file in the polars-bio session with `register_bgen(path, table_name, ...)`, run `SELECT * FROM <table> LIMIT max_rows+1`, convert the rows, and always deregister the table afterwards. A table of the same name already in the session raises `ConfigurationError`.
- **API:** `polars_bio_register_bgen(path: 'str', table_name: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.bgen` file. `table_name` is lowercase letters, digits and underscores, starting with a letter or underscore, at most 63 characters. The other arguments are as in `polars_bio_read_bgen()`.
- **Output:** `status`, `operation`, `table_name`, `columns`, `row_count`, `rows` and `provenance`. The table does not exist after the call.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_register_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_register_bgen("cohort.bgen", "cohort_bgen", verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([row["rsid"] for row in result["rows"]])
```

- **Example output:**

```text
cohort_bgen 2
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
['rs1', 'rs2']
```

## 20) `dnakit.annotation.polars_bio_register_cram()`

- **Purpose:** Load a local CRAM file as a temporary SQL table and return all of its rows.
- **Method:** Register the file with polars-bio's `register_cram(path, table_name, tag_fields=...)`, run `SELECT * FROM <table> LIMIT max_rows+1`, convert the rows, and deregister the table. polars-bio takes no reference FASTA here, so only CRAM files with embedded sequences decode.
- **API:** `polars_bio_register_cram(path: 'str', table_name: 'str', *, tag_fields: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`). The example writes its CRAM file with `pysam`.
- **Input:** `path` is a local CRAM file. `table_name` follows the rule given for `polars_bio_register_bgen()`. `tag_fields` is an optional non-empty list of alignment tags to add as columns. `max_rows` is an integer from 1 to 1,000,000.
- **Output:** `status`, `operation`, `table_name`, `columns`, `row_count`, `rows` and `provenance`.
- **Example code:**

```python
import pysam

from dnakit.annotation import polars_bio_register_cram

header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
with pysam.AlignmentFile("unmapped.cram", "wc", header=header) as out:
    for index, sequence in enumerate(["ACGTACGTAC", "GGGGCCCCAA"]):
        read = pysam.AlignedSegment()
        read.query_name = f"r{index}"
        read.query_sequence = sequence
        read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        read.flag = 4
        read.reference_id = -1
        read.reference_start = -1
        read.set_tag("NM", 0)
        out.write(read)

result = polars_bio_register_cram("unmapped.cram", "reads_cram", tag_fields=["NM"], verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([(row["name"], row["sequence"], row["NM"]) for row in result["rows"]])
```

- **Example output:**

```text
reads_cram 2
['name', 'chrom', 'start', 'end', 'flags', 'cigar', 'mapping_quality', 'mate_chrom', 'mate_start', 'sequence', 'quality_scores', 'template_length', 'NM']
[('r0', 'ACGTACGTAC', 0), ('r1', 'GGGGCCCCAA', 0)]
```

## 21) `dnakit.annotation.polars_bio_register_pgen()`

- **Purpose:** Load a local PGEN fileset as a temporary SQL table and return all of its rows.
- **Method:** Register the fileset with polars-bio's `register_pgen(path, table_name, ...)`, run `SELECT * FROM <table> LIMIT max_rows+1`, convert the rows, and deregister the table.
- **API:** `polars_bio_register_pgen(path: 'str', table_name: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** `path` is a local `.pgen` file; `table_name` follows the rule given for `polars_bio_register_bgen()`; the other arguments are as in `polars_bio_read_pgen()`.
- **Output:** `status`, `operation`, `table_name`, `columns`, `row_count`, `rows` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_register_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_register_pgen("cohort.pgen", "cohort_pgen", verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([row["genotypes"]["GT"] for row in result["rows"]])
```

- **Example output:**

```text
cohort_pgen 2
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
[[[0, 0], [0, 1], [1, 1]], [[0, 1], [0, 0], None]]
```

## 22) `dnakit.annotation.polars_bio_register_vcf_zarr()`

- **Purpose:** Load a local VCF Zarr store as a temporary SQL table and return all of its rows.
- **Method:** Register the store with polars-bio's `register_vcf_zarr(path, table_name, ...)`, run `SELECT * FROM <table> LIMIT max_rows+1`, convert the rows, and deregister the table.
- **API:** `polars_bio_register_vcf_zarr(path: 'str', table_name: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`). Building the example store also needs `zarr`.
- **Input:** `path` is an existing local directory; `table_name` follows the rule given for `polars_bio_register_bgen()`; the other arguments are as in `polars_bio_read_vcf_zarr()`.
- **Output:** `status`, `operation`, `table_name`, `columns`, `row_count`, `rows` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import polars_bio_register_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_register_vcf_zarr("cohort.vcz", "cohort_vcz", info_fields=["DP"], verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([(row["id"], row["DP"]) for row in result["rows"]])
```

- **Example output:**

```text
cohort_vcz 2
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'qual', 'filter', 'DP', 'genotypes']
[('rs1', 14), ('rs2', 20)]
```

## 23) `dnakit.annotation.polars_bio_scan_bgen()`

- **Purpose:** Scan a local BGEN file lazily, count all variants, and return at most `max_rows` of them.
- **Method:** Open the file with polars-bio's `scan_bgen(path, ...)`, count rows with a lazy `pl.len()`, then collect only the first `max_rows` rows. The result reports the column schema and whether the output was truncated.
- **API:** `polars_bio_scan_bgen(path: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** As in `polars_bio_read_bgen()`. Unlike the read function, a file with more than `max_rows` variants is cut, not refused.
- **Output:** `status`, `operation`, `columns` (name and dtype pairs), `variant_count` (all variants), `returned_count`, `truncated`, `variants` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_scan_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_scan_bgen("cohort.bgen", max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["rsid"])
```

- **Example output:**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
rs1
```

## 24) `dnakit.annotation.polars_bio_scan_pgen()`

- **Purpose:** Scan a local PGEN fileset lazily, count all variants, and return at most `max_rows` of them.
- **Method:** Open the fileset with polars-bio's `scan_pgen(path, ...)`, count rows lazily, then collect only the first `max_rows` rows.
- **API:** `polars_bio_scan_pgen(path: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`).
- **Input:** As in `polars_bio_read_pgen()`.
- **Output:** `status`, `operation`, `columns`, `variant_count`, `returned_count`, `truncated`, `variants` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_scan_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_scan_pgen("cohort.pgen", max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["id"])
```

- **Example output:**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
rs1
```

## 25) `dnakit.annotation.polars_bio_scan_vcf_zarr()`

- **Purpose:** Scan a local VCF Zarr store lazily, count all variants, and return at most `max_rows` of them.
- **Method:** Open the store with polars-bio's `scan_vcf_zarr(path, ...)`, count rows lazily, then collect only the first `max_rows` rows.
- **API:** `polars_bio_scan_vcf_zarr(path: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package polars-bio (extra `dnakit[bio-polarsbio]`). Building the example store also needs `zarr`.
- **Input:** As in `polars_bio_read_vcf_zarr()`.
- **Output:** `status`, `operation`, `columns`, `variant_count`, `returned_count`, `truncated`, `variants` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import polars_bio_scan_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_scan_vcf_zarr("cohort.vcz", info_fields=["DP"], max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["id"], result["variants"][0]["DP"])
```

- **Example output:**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'qual', 'filter', 'DP', 'genotypes']
rs1 14
```

## 26) `dnakit.annotation.polars_bio_visualize_intervals()`

- **Purpose:** Draw the overlapping interval pairs of two BED files into one PNG, one panel per pair.
- **Method:** Read both files with polars-bio's `read_bed(..., use_zero_based=True)`, compute `overlap(...)`, and for each of the first `max_pairs` pairs call polars-bio's `visualize_intervals`, which opens its own matplotlib figure. The figures are captured and stacked into one image written to `output_path`. If the files have no overlap, `ConfigurationError` is raised.
- **API:** `polars_bio_visualize_intervals(left_bed: 'str', right_bed: 'str', output_path: 'str', *, label: 'str' = 'overlapping pair', max_pairs: 'int' = 20, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages polars-bio, matplotlib (extra `dnakit[bio-plots]`) and bioframe (installed by `polars-bio[viz]`; a missing bioframe raises `BackendUnavailableError`).
- **Input:** `left_bed` and `right_bed` are local BED files. `output_path` is a new `.png` path (an existing file is not overwritten; missing parent folders are created). `label` is a non-empty string of at most 128 characters. `max_pairs` is an integer from 1 to 50.
- **Output:** `path` (the PNG) and `data` with `label`, `pair_count` (all overlaps), `plotted_count`, `pairs` (the drawn pairs, with `_1` and `_2` suffixed BED columns) and `provenance`. Intervals are 0-based half-open.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import polars_bio_visualize_intervals

target = Path(tempfile.mkdtemp()) / "pairs.png"
result = polars_bio_visualize_intervals("a.bed", "b.bed", str(target), verbose=False)
print(result["data"]["pair_count"], result["data"]["plotted_count"])
for pair in result["data"]["pairs"]:
    print(pair["name_1"], pair["start_1"], pair["end_1"], pair["name_2"], pair["start_2"], pair["end_2"])
print(target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **Example output:**

```text
2 2
a 10 20 c 15 30
b 40 60 d 50 55
True
```

## 27) `dnakit.annotation.cached_ensembl_release()`

- **Purpose:** Describe an Ensembl release that is already installed in the local pyensembl cache.
- **Method:** Call `pyensembl.ensembl_release.cached_release(release, species=species)`, then `inspect_data()`, which reads only the cache. If the release is installed, the result lists gene and contig counts from the annotation index and the cached files; if not, `BackendUnavailableError` is raised with the `pyensembl install` command as its hint. Nothing is downloaded.
- **API:** `cached_ensembl_release(release: 'int | str', *, species: 'str' = 'human', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and an installed release (for example from `pyensembl install --release 110 --species human`, or from `pyensembl_install_ensembl_release()`). The cache is `$PYENSEMBL_CACHE_DIR` or `~/.cache/pyensembl`. The example points the cache at an empty folder and shows the not-installed error; the usage block shows the call on an installed release.
- **Input:** `release` is a positive integer or an annotation date string such as `2026_04`. `species` is a non-empty species name known to pyensembl; an unknown name raises `ConfigurationError`.
- **Output:** `status`, `operation`, `service`, `backend`, `version`, `request`, `retrieved_utc` and `data` with `species`, `reference_name`, `release`, `annotation_name`, `annotation_version`, `cache_directory`, `installed`, `annotation_status`, `reference_dna_status`, `gene_count`, `contig_count` and `files`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["PYENSEMBL_CACHE_DIR"] = str(Path("empty_cache").resolve())

from dnakit.annotation import cached_ensembl_release
from dnakit.exceptions import BackendUnavailableError

try:
    cached_ensembl_release(110, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.code)
    print(error.hint.split("; ")[0])
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
Install the annotation with the pyensembl command line first, for example: pyensembl install --release 110 --species homo_sapiens
```

Usage with the dependency installed:

```python
from dnakit.annotation import cached_ensembl_release

result = cached_ensembl_release(110)  # after: pyensembl install --release 110 --species human
print(result["data"]["gene_count"], result["data"]["contig_count"])
```

## 28) `dnakit.annotation.ensembl_release_for_reference_name()`

- **Purpose:** Describe the installed Ensembl release that pyensembl picks for a reference genome name such as GRCh38.
- **Method:** Call `pyensembl.reference_name.genome_for_reference_name(name, allow_older_downloaded_release=...)`. With the flag true the newest installed release is chosen, otherwise the newest release with downloaded files. The chosen release is then described exactly as in `cached_ensembl_release()`, and `BackendUnavailableError` is raised if it is not installed. Nothing is downloaded.
- **API:** `ensembl_release_for_reference_name(reference_name: 'str', *, allow_older_downloaded_release: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and an installed release for the reference. The example shows the refusal of an unknown reference name; the usage block shows an installed one.
- **Input:** `reference_name` is a non-empty string such as `GRCh38`; an unknown name raises `ConfigurationError`. `allow_older_downloaded_release` must be a bool.
- **Output:** The same dictionary as `cached_ensembl_release()`, with `request` holding `reference_name` and `allow_older_downloaded_release`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["PYENSEMBL_CACHE_DIR"] = str(Path("empty_cache").resolve())

from dnakit.annotation import ensembl_release_for_reference_name
from dnakit.exceptions import ConfigurationError

try:
    ensembl_release_for_reference_name("GRCh99", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
pyensembl genome_for_reference_name rejected the input: Reference genome 'GRCh99' not found
```

Usage with the dependency installed:

```python
from dnakit.annotation import ensembl_release_for_reference_name

result = ensembl_release_for_reference_name("GRCh38")  # needs an installed GRCh38 release
print(result["data"]["release"], result["data"]["gene_count"])
```

## 29) `dnakit.annotation.fetch_all_ensembl_dated_releases()`

- **Purpose:** Collect the annotation dates of every species that has dated Ensembl releases.
- **Method:** Call `pyensembl.dated_releases.fetch_all_dated_releases(expire_after=..., timeout=...)`. Cached listings younger than `expire_after_days` are used as they are; the human listing is fetched first, and if it cannot be read every species falls back to its cached dates. Species with unknown dates map to `None`. If nothing can be fetched and nothing is cached, `BackendExecutionError` is raised.
- **API:** `fetch_all_ensembl_dated_releases(*, expire_after_days: 'float' = 1.0, timeout_seconds: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and internet access for listings that are missing or expired. The example shows the argument check, which runs before any network access.
- **Input:** `expire_after_days` is a finite number of zero or more (a negative value raises `ConfigurationError`). `timeout_seconds` is a positive integer.
- **Output:** The standard listing envelope (`status`, `operation`, `service`, `backend`, `version`, `request`, `retrieved_utc`) and `data` with `dates_by_species`, `species_count`, `species_with_dates` and `species_without_dates`.
- **Example code:**

```python
from dnakit.annotation import fetch_all_ensembl_dated_releases
from dnakit.exceptions import ConfigurationError

try:
    fetch_all_ensembl_dated_releases(expire_after_days=-1, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
expire_after_days must be a finite number of zero or more.
```

Usage with the dependency installed:

```python
from dnakit.annotation import fetch_all_ensembl_dated_releases

result = fetch_all_ensembl_dated_releases(expire_after_days=7)
print(result["data"]["species_count"], result["data"]["dates_by_species"]["mus_musculus"])
```

## 30) `dnakit.annotation.list_ensembl_annotation_dates()`

- **Purpose:** List the annotation dates Ensembl publishes for one assembly accession.
- **Method:** Validate the accession and provider, then call `pyensembl.dated_releases.available_annotation_dates(accession, provider, refresh=...)`. The listing is read from Ensembl the first time and reused from the cache afterwards unless `refresh` is true. Dates are returned oldest first in `YYYY_MM` form.
- **API:** `list_ensembl_annotation_dates(assembly_accession: 'str', *, provider: 'str' = 'ensembl', refresh: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and internet access for a listing that is not cached. The example shows the accession check, which runs before any network access.
- **Input:** `assembly_accession` is a versioned GenBank or RefSeq accession matching `GCA_` or `GCF_` plus nine digits, a dot and a version, for example `GCA_000001405.29`. `provider` is a lower-case directory name (letters and underscores, at most 32). `refresh` is a bool.
- **Output:** The standard listing envelope and `data` with `dates` (strings) and `count`. A listing that cannot be read raises `BackendExecutionError`; an unknown accession raises `ConfigurationError`.
- **Example code:**

```python
from dnakit.annotation import list_ensembl_annotation_dates
from dnakit.exceptions import ConfigurationError

try:
    list_ensembl_annotation_dates("GRCh38", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
assembly_accession must be a versioned GCA or GCF accession, for example GCA_000001405.29.
```

Usage with the dependency installed:

```python
from dnakit.annotation import list_ensembl_annotation_dates

result = list_ensembl_annotation_dates("GCA_000001405.29")
print(result["data"]["count"], result["data"]["dates"][:3])
```

## 31) `dnakit.annotation.list_ensembl_dated_releases()`

- **Purpose:** List the annotation dates published for a species' current assembly.
- **Method:** Resolve the species with pyensembl's registry (`check_species_object`); a species without dated releases is refused before any network access. Otherwise call `pyensembl.dated_releases.available_dated_releases(species, refresh=...)`, which reads the Ensembl listing for the species' assembly accession and provider, or reuses the cached listing.
- **API:** `list_ensembl_dated_releases(species: 'str' = 'human', *, refresh: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and internet access unless a cached listing exists and `refresh` is false. The example shows the refusal of a species without dated releases, which happens offline.
- **Input:** `species` is a species name known to pyensembl (an unknown name raises `ConfigurationError`). `refresh` is a bool.
- **Output:** The standard listing envelope and `data` with `latin_name`, `assembly_accession`, `provider`, `dates` and `count`.
- **Example code:**

```python
from dnakit.annotation import list_ensembl_dated_releases
from dnakit.exceptions import ConfigurationError

try:
    list_ensembl_dated_releases("toxoplasma_gondii", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
toxoplasma_gondii has no dated Ensembl releases.
```

Usage with the dependency installed:

```python
from dnakit.annotation import list_ensembl_dated_releases

result = list_ensembl_dated_releases("mouse")
print(result["data"]["assembly_accession"], result["data"]["dates"])
```

## 32) `dnakit.annotation.pyensembl_install_local_annotation()`

- **Purpose:** Build a pyensembl index from a local GTF file (and optional local FASTA files) so it can be queried like an Ensembl release.
- **Method:** Create a `pyensembl.Genome` from the local paths, then call its `download()` (which leaves local files where they are) and `index()`, as the `pyensembl install` action does with `--gtf`. It never uses the network. The result counts the indexed genes and transcripts.
- **API:** `pyensembl_install_local_annotation(reference_name: 'str', annotation_name: 'str', gtf_path: 'str | Path', *, annotation_version: 'int | None' = None, transcript_fasta_paths: 'Sequence[str | Path]' = (), protein_fasta_paths: 'Sequence[str | Path]' = (), genome_fasta_path: 'str | Path | None' = None, cache_directory_path: 'str | Path | None' = None, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`). It writes the index into the cache folder.
- **Input:** `reference_name` and `annotation_name` are names (at most 128 characters) under which pyensembl keeps the genome. `gtf_path` is an existing local file. `annotation_version` is a positive integer or `None`. The FASTA arguments are lists (at most 10 files) of existing local files, or a single file for `genome_fasta_path`. `cache_directory_path` must be an existing folder, otherwise pyensembl's default cache is used. `overwrite=True` rebuilds an existing index.
- **Output:** `status` (`installed` or `not_installed`), `operation`, `reference_name`, `annotation_name`, `annotation_version`, `gtf` (file name), `cache_directory`, `genes`, `transcripts` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import pyensembl_install_local_annotation

Path("toy_cache").mkdir(exist_ok=True)
result = pyensembl_install_local_annotation(
    "toyref", "toyann", "mini.gtf", cache_directory_path="toy_cache", verbose=False
)
print(result["status"], result["gtf"], result["genes"], result["transcripts"])
print(sorted(path.name for path in Path("toy_cache").iterdir()))
```

- **Example output:**

```text
installed mini.gtf 2 2
['mini.db']
```

## 33) `dnakit.annotation.pyensembl_install_ensembl_release()`

- **Purpose:** Download and index one Ensembl release for a species into the pyensembl cache.
- **Method:** Create `pyensembl.EnsemblRelease(release=..., species=...)` and call `download()` and `index()`, as `pyensembl install --release` does. pyensembl fetches the GTF annotation and the transcript and protein FASTA files from the Ensembl FTP server and builds the gene index. The result records the release, the species, the GTF URL, and the gene and transcript counts.
- **API:** `pyensembl_install_ensembl_release(release: 'int | str', species: 'str' = 'human', *, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyensembl (extra `dnakit[bio-pyensembl]`) and internet access to the Ensembl FTP server; the download is large for release annotations and genome FASTA sets. The example shows the argument check that runs before any download.
- **Input:** `release` is an integer from 1 to 1000 or a string such as `2026_04`. `species` is a pyensembl species name such as `human` or `caenorhabditis_elegans`. `overwrite` is a bool. A species or release that pyensembl rejects raises `ConfigurationError`.
- **Output:** `status` (`installed` or `not_installed`), `operation`, `ensembl_release`, `species`, `gtf_url`, `cache_directory`, `genes`, `transcripts` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import pyensembl_install_ensembl_release
from dnakit.exceptions import ConfigurationError

try:
    pyensembl_install_ensembl_release(0, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
release must be between 1 and 1000.
```

Usage with the dependency installed:

```python
from dnakit.annotation import pyensembl_install_ensembl_release

result = pyensembl_install_ensembl_release(110, species="caenorhabditis_elegans")
print(result["status"], result["genes"], result["transcripts"])
```
