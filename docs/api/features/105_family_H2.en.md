# H2 features: family-variant check of pydna, DnaCauldron, DnaChisel and primer3

This page records the check of 101 `family-variant` ledger rows of pydna, DnaCauldron, DnaChisel and primer3, and the 10 wrapper functions that the check added. They live in three modules under `dnakit.molecular_lab`: `pydna_family_functions`, `dnacauldron_family_functions` and `dnachisel_family_functions`. primer3 gets no new function and so no module.

- Offline: the examples use only call arguments or example files and never use the network.
- Method: `duplicate-of` compares the two upstream objects with `is`, or compares the result of the same code; a kept `family-variant` row is called through the parameter of the existing wrapper and compared with the direct upstream call; the other rows get a more precise status or a new wrapper. All checks are in `others/tests/unit/test_family_H2.py`.
- Examples import from the module that holds the function, for example `from dnakit.molecular_lab import simulate_pcr_assembly`.

## Ledger status (101 rows)

| Package | family-variant before | duplicate-of | family-variant (verified) | wrapped | Other statuses |
| --- | ---: | ---: | ---: | --- | --- |
| pydna | 59 | 34 | 13 | 7 rows (6 functions) | deprecated 3, broken-upstream 1, needs-network 1 |
| dnacauldron | 21 | 18 | 0 | 3 rows (2 functions) | none |
| dnachisel | 13 | 11 | 0 | 2 rows (2 functions) | none |
| primer3 | 8 | 8 | 0 | 0 | none |
| Total | 101 | 71 | 13 | 12 rows (10 functions) | 5 |

## Main findings

- `pydna.readers.read` and two other rows were called variants of `parse_primer_sequences`. That is wrong: `read` returns a DNA record, not a primer; it is `parse` plus a check for exactly one record, and `parse_sequence_text` already covers it.
- `pydna.assembly2.common_sub_strings` is not an alias of the string version: it takes Dseqrecords and can cross the origin of a circular molecule (the example gives `(4, 0, 11)` instead of `(4, 0, 7)`), so it got a wrapper.
- `pydna.assembly2.pcr_assembly` differs from `amplify.pcr`: it allows primer mismatches, supports circular templates and can give several products.
- `tm_dbd`, `tm_product`, `ta_default`, `ta_dbd` and `tmbresluc` in `pydna.tm` were filed under `primer_melting_temperature`, but that wrapper has no parameter that reaches them, so each got a wrapper. `tm_neb` needs the NEB web service (`needs-network`), and the body of `Q5` is only `raise NotImplementedError` (`broken-upstream`).
- `homologous_recombination_excision`, `cre_lox_excision` and `recombinase_excision` are deprecated in pydna 5.5.16 (they emit a DeprecationWarning and point to `*_excision_or_inversion`), so they are now `deprecated`; `simulate_genome_edit(operation="excision")` still reaches them.
- `pfu_sso7d_program` is only an alias of `dbd_program`, so `polymerase="dbd"` and `"pfu_sso7d"` of `suggest_pcr_program` give the same text.
- Most other rows were re-exports of the same object (the `snapgene_history_parser`, `all` and `design` modules), so they are now `duplicate-of`.
- DnaCauldron's `load_records_from_files`, `load_records_from_file` and `string_to_records` were said to be covered by `read_dnacauldron_parts`; that wrapper reads one record per file by extension and cannot read multi-record files, archives, folders or bare sequence text.
- `annotate_loxP_sites` and `annotate_gateway_sites` only write the sites found by `find_*_sites` as features; their positions equal the result of `find_recombination_sites`.
- The six primer3 `calc_*_tm` functions are `calc_*(...).tm`; the two `_b` functions are bytes forms, and `sanitize_sequence` calls `sanitize_sequence_b` internally.

## 1) `find_common_substrings_with_topology`

- **Purpose:** Find common substrings of two DNA molecules, across the origin of circular ones.
- **Method:** Calls pydna.assembly2.common_sub_strings. A circular molecule is read as if its ends were joined, so a match can run from the end to the start and is never longer than the shorter molecule; for two linear molecules the result equals `find_common_substrings`. Each match is `(x_start, y_start, length)` with 0-based starts, sorted and capped at 1000 rows.
- **API:** `find_common_substrings_with_topology(sequence_x[required], sequence_y[required], min_length[optional], circular_x[optional], circular_y[optional], verbose[optional])`
- **Input:** Two A/C/G/T sequences, a minimum match length (default 25) and two booleans, `circular_x` and `circular_y`.
- **Output:** A dictionary with `min_length`, `circular_x`, `circular_y`, `match_count`, `matches` (`x_start`, `y_start`, `length`) and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import find_common_substrings_with_topology

x = "ACGTTTGCAAC"
y = "TTGCAACACGTT"
circular = find_common_substrings_with_topology(x, y, min_length=6, circular_x=True, verbose=False)
linear = find_common_substrings_with_topology(x, y, min_length=6, verbose=False)
print(circular["matches"])
print(linear["matches"])
```

**Example output:**

```text
[{'x_start': 4, 'y_start': 0, 'length': 11}]
[{'x_start': 4, 'y_start': 0, 'length': 7}]
```

## 2) `simulate_pcr_assembly`

- **Purpose:** Enumerate every PCR product of a primer pair, allowing primer mismatches.
- **Method:** Calls pydna.assembly2.pcr_assembly. Each primer must anneal over at least `min_annealing` bases, and `mismatches` allows that many mismatched bases in the annealing region, which `simulate_pcr_amplicon` does not. A circular template can give a product across the origin, and several binding sites give several products. Products are sorted by length and sequence; an empty list means no product.
- **API:** `simulate_pcr_assembly(template[required], forward_primer[required], reverse_primer[required], min_annealing[optional], mismatches[optional], template_circular[optional], max_products[optional], verbose[optional])`
- **Input:** A template, a forward and a reverse primer (A/C/G/T), a minimum annealing length (default 14), the allowed mismatches (0 to 5), whether the template is circular and the most products to list.
- **Output:** A dictionary with `template_length`, `template_circular`, `min_annealing`, `mismatches`, `product_count`, `products` (`sequence`, `length`, `circular`) and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import simulate_pcr_assembly

template = "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGGATGATCAGTGGGTAAAGG"
forward = "AGAGTTTTTCGAACTCGTGTTG"
reverse = "CCCACGACTAAAAGCCCTGCCA"
result = simulate_pcr_assembly(template, forward, reverse, verbose=False)
print(result["product_count"], result["products"][0]["length"], result["products"][0]["sequence"] == template[20:102])
changed = forward[:10] + ("A" if forward[10] != "A" else "C") + forward[11:]
strict = simulate_pcr_assembly(template, changed, reverse, verbose=False)
tolerant = simulate_pcr_assembly(template, changed, reverse, mismatches=1, verbose=False)
print(strict["product_count"], tolerant["product_count"], [row["length"] for row in tolerant["products"]])
```

**Example output:**

```text
1 82 True
0 1 [82]
```

## 3) `primer_melting_temperature_dbd`

- **Purpose:** Compute a primer melting temperature with the rule for DNA-binding-domain polymerases such as Phusion and Pfu-Sso7d.
- **Method:** Calls pydna.tm.tm_dbd: nearest-neighbour method, SantaLucia 1998 table, 50 mM Na+, 1.5 mM Mg2+ and 250 nM of each strand. The value differs from `primer_melting_temperature`, which uses the Taq buffer. It is rounded to three decimals.
- **API:** `primer_melting_temperature_dbd(sequence[required], verbose[optional])`
- **Input:** A primer of 6 to 200 A/C/G/T bases.
- **Output:** A dictionary with `sequence`, `length`, `tm_celsius` and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import primer_melting_temperature
from dnakit.molecular_lab import primer_melting_temperature_dbd

primer = "AGAGTTTTTCGAACTCGTGTTG"
print(primer_melting_temperature(primer, verbose=False)["tm_celsius"])
print(primer_melting_temperature_dbd(primer, verbose=False)["tm_celsius"])
```

**Example output:**

```text
59.445
55.982
```

## 4) `primer_melting_temperature_breslauer`

- **Purpose:** Compute a primer melting temperature with the Breslauer parameters.
- **Method:** Calls pydna.tm.tmbresluc: Breslauer nearest-neighbour enthalpy and entropy values, corrected for the primer concentration `primer_nm` (nM, default 500) and the monovalent cation concentration `salt_mm` (mM, default 50). Doubling the salt adds 16.6 times log10(2), about 5.0 degrees. The value is rounded to three decimals.
- **API:** `primer_melting_temperature_breslauer(sequence[required], primer_nm[optional], salt_mm[optional], verbose[optional])`
- **Input:** A primer of 6 to 200 A/C/G/T bases, a primer concentration and a salt concentration (both positive).
- **Output:** A dictionary with `sequence`, `length`, `primer_nm`, `salt_mm`, `tm_celsius` and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import primer_melting_temperature_breslauer

primer = "AGAGTTTTTCGAACTCGTGTTG"
base = primer_melting_temperature_breslauer(primer, verbose=False)
salty = primer_melting_temperature_breslauer(primer, salt_mm=100.0, verbose=False)
print(base["tm_celsius"], salty["tm_celsius"])
print(round(salty["tm_celsius"] - base["tm_celsius"], 2))
```

**Example output:**

```text
59.082 64.079
5.0
```

## 5) `amplicon_melting_temperature`

- **Purpose:** Estimate the melting temperature of a PCR product from its length and GC content.
- **Method:** Calls pydna.tm.tm_product, the Rychlik et al. (1990) formula: 81.5 + 0.41 times GC percent + 16.6 times log10(K) - 675 divided by the length, where K is the monovalent cation concentration `potassium_m` (mol/L, default 0.05). The value is rounded to three decimals.
- **API:** `amplicon_melting_temperature(sequence[required], potassium_m[optional], verbose[optional])`
- **Input:** A product of at least 20 A/C/G/T bases and a cation concentration (positive).
- **Output:** A dictionary with `length`, `gc_fraction`, `potassium_m`, `tm_celsius` and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import amplicon_melting_temperature

product = "AGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGG"
result = amplicon_melting_temperature(product, verbose=False)
print(result["length"], result["gc_fraction"], result["tm_celsius"])
print(amplicon_melting_temperature(product, potassium_m=0.5, verbose=False)["tm_celsius"])
```

**Example output:**

```text
82 0.463415 70.671
87.271
```

## 6) `pcr_annealing_temperature`

- **Purpose:** Suggest a PCR annealing temperature for a primer pair and its product.
- **Method:** `polymerase="taq"` (default) calls pydna.tm.ta_default: 0.3 times the lower primer Tm + 0.7 times the product Tm - 14.9. `polymerase="dbd"` calls pydna.tm.ta_dbd: the lower primer Tm (DBD rule) plus 3 degrees. The Tm of both primers is reported too; all values are rounded to three decimals.
- **API:** `pcr_annealing_temperature(forward_primer[required], reverse_primer[required], product[required], polymerase[optional], verbose[optional])`
- **Input:** The annealing parts of two primers (6 to 200 bases), a product of at least 20 bases and the polymerase type.
- **Output:** A dictionary with `polymerase`, `forward_tm_celsius`, `reverse_tm_celsius`, `product_length`, `annealing_celsius` and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import pcr_annealing_temperature

forward = "AGAGTTTTTCGAACTCGTGTTG"
reverse = "CCCACGACTAAAAGCCCTGCCA"
product = "AGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGG"
taq = pcr_annealing_temperature(forward, reverse, product, verbose=False)
dbd = pcr_annealing_temperature(forward, reverse, product, polymerase="dbd", verbose=False)
print(taq["forward_tm_celsius"], taq["reverse_tm_celsius"], taq["annealing_celsius"])
print(dbd["forward_tm_celsius"], dbd["reverse_tm_celsius"], dbd["annealing_celsius"])
```

**Example output:**

```text
59.445 68.485 52.403
55.982 64.16 58.982
```

## 7) `read_dnacauldron_records`

- **Purpose:** Read every record of local sequence files with DnaCauldron.
- **Method:** Calls `load_records_from_files` and `load_records_from_file` of dnacauldron.biotools.sequence_io. Unlike `read_dnacauldron_parts`, a file may hold several records, the format is detected from the content, `.zip` archives are opened, and `folder` reads a whole directory (sorted by path). A record without an ID takes its file name, with a number when a file holds several records; with `use_file_names_as_ids=True` a file with one record gives that record the file name as ID.
- **API:** `read_dnacauldron_records(sources[optional], folder[optional], use_file_names_as_ids[optional], max_records[optional], verbose[optional])`
- **Input:** A list of file paths `sources` or a directory `folder` (give at least one), whether file names become IDs and the most records (default 500).
- **Output:** A dictionary with `file_count`, `record_count`, `files` (`file`, `format`, `record_count`; the format of an archive is `zip`), `records` (`id`, `file`, `length`, `topology`, `sequence`) and `provenance`.

**Example file:** `parts.fa`

```text
>partA
ATGCGTACGTTAGC
>partB
GGCATCGATCGGA
```

**Example file:** `backbone.fasta`

```text
>
ACGTACGTAGCTAGCTAGC
```

**Example code:**

```python
from dnakit.molecular_lab import read_dnacauldron_records

result = read_dnacauldron_records(["parts.fa", "backbone.fasta"], verbose=False)
print(result["file_count"], result["record_count"])
print([(item["file"], item["format"], item["record_count"]) for item in result["files"]])
print([(row["id"], row["length"]) for row in result["records"]])
```

**Example output:**

```text
2 3
[('parts.fa', 'fasta', 2), ('backbone.fasta', 'fasta', 1)]
[('partA', 14), ('partB', 13), ('backbone', 19)]
```

## 8) `parse_dnacauldron_sequence_text`

- **Purpose:** Detect the format of sequence text and read its records.
- **Method:** Calls dnacauldron.biotools.sequence_io.string_to_records. The text is a bare ACGT string (format `ATGC`, one record without a name), FASTA or GenBank. Unlike `parse_sequence_text`, a bare sequence is accepted and the detected format is reported; any other format raises `ConfigurationError`.
- **API:** `parse_dnacauldron_sequence_text(text[required], verbose[optional])`
- **Input:** Non-empty text of at most 1 MB.
- **Output:** A dictionary with `format`, `record_count`, `records` (`id`, `length`, `topology`, `sequence`) and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import parse_dnacauldron_sequence_text

bare = parse_dnacauldron_sequence_text("ACGTACGT", verbose=False)
print(bare["format"], bare["records"][0]["sequence"])
fasta = parse_dnacauldron_sequence_text(">x\nACGT\n>y\nGGCC\n", verbose=False)
print(fasta["format"], [row["id"] for row in fasta["records"]])
```

**Example output:**

```text
ATGC ACGTACGT
fasta ['x', 'y']
```

## 9) `find_sequence_difference_positions`

- **Purpose:** Return a per-position difference vector for two equal-length DNA sequences.
- **Method:** Calls dnachisel.biotools.sequences_differences.sequences_differences_array. `difference_array` has one value per position: 1 where the sequence differs from the reference and 0 where it agrees; `positions` lists the differing 0-based positions. Unlike `find_sequence_difference_segments`, it does not merge positions into segments. No alignment is performed.
- **API:** `find_sequence_difference_positions(sequence[required], reference[required], verbose[optional])`
- **Input:** Two A/C/G/T sequences of the same length.
- **Output:** A dictionary with `sequence_length`, `difference_count`, `positions`, `difference_array` and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import find_sequence_difference_positions

result = find_sequence_difference_positions("ACGTACGTAC", "ACCTACGAAC", verbose=False)
print(result["difference_count"], result["positions"])
print(result["difference_array"])
```

**Example output:**

```text
2 [2, 7]
[0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
```

## 10) `annotate_sequence_differences`

- **Purpose:** Annotate where a DNA sequence differs from a reference as labelled features.
- **Method:** Calls dnachisel.biotools.annotate_differences. Each run of adjacent differing positions becomes one feature (`start` inclusive, `end` exclusive) labelled with `prefix` followed by the reference bases at that location, for example `#A` or `#TT`, which `find_sequence_difference_segments` does not give. The sequences must have the same length; identical sequences give zero features (upstream raises IndexError in that case).
- **API:** `annotate_sequence_differences(sequence[required], reference[required], feature_type[optional], prefix[optional], verbose[optional])`
- **Input:** Two A/C/G/T sequences of the same length, a feature type (default `misc_feature`) and a label prefix (default `#`, at most 20 characters).
- **Output:** A dictionary with `sequence_length`, `feature_count`, `features` (`start`, `end`, `type`, `label`) and `provenance`.

**Example code:**

```python
from dnakit.molecular_lab import annotate_sequence_differences

result = annotate_sequence_differences("AATTGGCCAC", "AACCGGCCGC", verbose=False)
for row in result["features"]:
    print(row["start"], row["end"], row["type"], row["label"])
same = annotate_sequence_differences("ACGT", "ACGT", verbose=False)
print(same["feature_count"])
```

**Example output:**

```text
2 4 misc_feature #CC
8 9 misc_feature #G
0
```
