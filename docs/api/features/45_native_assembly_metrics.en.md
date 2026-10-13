# Native assembly contiguity metrics

These functions calculate contiguity metrics from explicit contig lengths. They do not run an assembler or query a reference genome or BUSCO database. For a methods reference, see Assemblathon 2: <https://doi.org/10.1038/nmeth.1935>.

## 1) `dnakit.assembly.calculate_n50()`

- **Purpose:** Calculate N50 and its corresponding L50.
- **Method:** Sort contigs in descending order and select the contig crossing 50% of total length.
- **API:** `dnakit.assembly.calculate_n50(input, options=None)`.
- **Input:** `[100, 60, 40]` or `{"lengths": [100, 60, 40]}`.
- **Output:** A JSON-serializable dictionary with `status`, `result`, and `provenance`.
- **Example code:** `dnakit.assembly.calculate_n50([100, 60, 40])`.
- **Example output:** `{"result": {"n50": 100, "l50": 1, "total_length": 200}}`.

## 2) `dnakit.assembly.calculate_nx()`

- **Purpose:** Calculate Nx at an arbitrary percentile.
- **Method:** Use `options.x` or `options.percentile` as the cumulative threshold.
- **API:** `dnakit.assembly.calculate_nx(input, options={"x": 90})`.
- **Input:** A list of positive contig lengths.
- **Output:** A result containing `nx`, the threshold rank, and total length.
- **Example code:** `dnakit.assembly.calculate_nx([100, 60, 40], options={"x": 90})`.
- **Example output:** `{"result": {"percentile": 90.0, "nx": 40, "contigs_to_threshold": 3}}`.

## 3) `dnakit.assembly.calculate_l50()`

- **Purpose:** Calculate the number of contigs needed to reach a cumulative percentile.
- **Method:** Use the same descending sort and cumulative rule as Nx, defaulting to 50%.
- **API:** `dnakit.assembly.calculate_l50(input, options=None)`.
- **Input:** A list of positive contig lengths.
- **Output:** `l50` or the rank at the requested percentile.
- **Example code:** `dnakit.assembly.calculate_l50([100, 60, 40])`.
- **Example output:** `{"result": {"l50": 1, "n_at_percentile": 100}}`.

## 4) `dnakit.assembly.calculate_ng50()`

- **Purpose:** Calculate NG50 with an expected genome-size denominator.
- **Method:** Use `options.genome_size` instead of assembled length as the cumulative denominator.
- **API:** `dnakit.assembly.calculate_ng50(input, options={"genome_size": 250})`.
- **Input:** Contig lengths and an explicit expected genome size.
- **Output:** A result containing `ng50`, the expected genome size, and assembled length.
- **Example code:** `dnakit.assembly.calculate_ng50([100, 60, 40], options={"genome_size": 250})`.
- **Example output:** `{"result": {"ng50": 60, "expected_genome_size": 250}}`.

## 5) `dnakit.assembly.calculate_stats()`

- **Purpose:** Return common assembly length, N50/L50, and N90/L90 metrics together.
- **Method:** Compute descriptive statistics over the same explicit contig lengths.
- **API:** `dnakit.assembly.calculate_stats(input, options=None)`.
- **Input:** Positive contig lengths; add `genome_size` to include NG50.
- **Output:** Counts, total/min/max/mean lengths, N50/L50, N90/L90, and sorted lengths.
- **Example code:** `dnakit.assembly.calculate_stats([100, 60, 40], options={"genome_size": 250})`.
- **Example output:** `{"result": {"contig_count": 3, "n50": 100, "ng50": 60}}`.

