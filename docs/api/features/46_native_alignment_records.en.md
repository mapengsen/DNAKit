# Native SAM records and alignment statistics

These functions process SAM text or materialized SAM-like records. They do not read binary BAM/CRAM, build a reference index, or replace BWA, Bowtie2, minimap2, or another external aligner. A methods reference for the file and record ecosystem is SAMtools/BCFtools/HTSlib: <https://doi.org/10.1093/gigascience/giab008>.

## 1) `dnakit.alignment.parse_sam_records()`

- **Purpose:** Parse SAM text into serializable records.
- **Method:** Read SAM records with at least 11 columns and parse flags, coordinates, CIGAR, sequence, quality, and tags.
- **API:** `dnakit.alignment.parse_sam_records(source)`.
- **Input:** SAM text or a local SAM path.
- **Output:** A list of record dictionaries with 0-based start positions.
- **Example code:** `dnakit.alignment.parse_sam_records("r1\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tATGC\\tIIII")`.
- **Example output:** `[{"query_name": "r1", "reference": "chr1", "position": 0, "cigar": "4M"}]`.

## 2) `dnakit.alignment.calculate_coverage()`

- **Purpose:** Calculate covered bases, total depth, and mean depth for SAM-like records.
- **Method:** Accumulate depth for CIGAR `M`, `=`, `X`, and reference-consuming operations.
- **API:** `dnakit.alignment.calculate_coverage(input)`.
- **Input:** `{"sam": "..."}` or `{"records": [...]}`.
- **Output:** `covered_bases`, `total_depth`, `mean_depth`, and depth by reference.
- **Example code:** `dnakit.alignment.calculate_coverage({"sam": sam_text})`.
- **Example output:** `{"covered_bases": 4, "total_depth": 4, "mean_depth": 1.0}`.

## 3) `dnakit.alignment.calculate_mapping_rate()`

- **Purpose:** Calculate the proportion of mapped records.
- **Method:** Use the SAM unmapped flag or a `*` reference.
- **API:** `dnakit.alignment.calculate_mapping_rate(input)`.
- **Input:** SAM text, a path, or materialized records.
- **Output:** A floating-point fraction between 0 and 1.
- **Example code:** `dnakit.alignment.calculate_mapping_rate({"records": records})`.
- **Example output:** `0.75`.

## 4) `dnakit.alignment.generate_pileup()`

- **Purpose:** Generate depth records ordered by reference and position.
- **Method:** Reuse local CIGAR depth calculation without calling samtools.
- **API:** `dnakit.alignment.generate_pileup(input)`.
- **Input:** SAM-like records.
- **Output:** `[{"reference": ..., "position": ..., "depth": ...}]`.
- **Example code:** `dnakit.alignment.generate_pileup({"sam": sam_text})`.
- **Example output:** `[{"reference": "chr1", "position": 0, "depth": 2}]`.

## 5) `dnakit.alignment.mark_duplicates()`

- **Purpose:** Mark duplicate records by reference, position, orientation, and sequence.
- **Method:** Keep the first occurrence and add a `DUP=1` tag to later records with the same key.
- **API:** `dnakit.alignment.mark_duplicates(input)`.
- **Input:** SAM-like records.
- **Output:** A list of records with duplicate tags.
- **Example code:** `dnakit.alignment.mark_duplicates({"records": records})`.
- **Example output:** `[{"tags": {"DUP": "1"}}]`.

