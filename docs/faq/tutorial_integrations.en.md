# Tutorial details and integrations

Installation, parameters and operational details for the feature tutorials.

## Reverse-complement files and performance {#reverse-complement}

Strings containing a dot or path separator are treated as paths; pass `Path` for extensionless files. Objects preserve their return types and map features, strands and per-base annotations.


Large uncompressed FASTA/FNA inputs (at least 16 MiB) with default reader settings
use file-region processing. Workers read, validate, transform, and write disjoint
regions directly; only paths and offsets cross process boundaries. Record order is
preserved and the temporary output is installed only after successful validation.
`workers=1` runs the same block algorithm serially; `workers=-1` uses available CPUs
subject to memory limits. Multiprocessing scripts need an `if __name__ == "__main__":`
guard. Irregular wrapping, compressed files, other formats, explicit `read_config`,
and output byte limits retain the general processing path and its validation rules.


DNA strings and sequence files are also supported directly:

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AATCGC"))  # GCGATT; accepts lowercase and whitespace
result = reverse_complement("hg38.fna")  # writes hg38_reverse_complement.fna
result = reverse_complement("sequences.txt", "rc.txt")
result = reverse_complement("sequences.csv", "rc.csv", sequence_column="sequence")
```

Strings containing a dot or path separator are treated as paths; use `Path` for extensionless files.
Files infer FASTA/FNA, TXT, CSV/TSV and other DNAKit I/O formats, including gzip, and stream one record
at a time. The result is a `WriteResult` with `record_count` and `target_artifact.relative_path`
(the output path relative to the current working directory).
File mode displays start/completion messages and processed record counts; set `show_progress=False`
to silence it. DNA strings return uppercase IUPAC strings; existing DNA object return types remain unchanged.

Without `read_config`, uncompressed input limits follow the actual file size, so hg38 requires no
manual size settings. Gzip input has a 20 GB decompressed reading ceiling. Pass `read_config=ReadConfig(...)` to set explicit limits. Input and output must differ; existing outputs require
`write_config=WriteConfig(overwrite=True)` to replace them.

TXT defaults to one sequence per non-empty line and writes one sequence per line without IDs.
Use `read_config=ReadConfig(text_mode="sequence")` for one wrapped sequence spanning the whole file.
CSV reads the `sequence` column by default, with `sequence_column` allowing another name, and preserves
the input column names, order, labels and other extra columns. Record annotations follow the reverse-complement coordinate mapping.

## Descriptor fields, conditions and sources {#descriptors}

The result object can be converted into a serializable dictionary via `to_dict()`. Where the interface provides `cross_gaps`, `True` is only allowed to span spanable Gap; the position of `Gap(crossable=False)` explicitly set is always a hard boundary. Descriptors that accept circular sequences still use the current origin as the scan boundary, and do not automatically complement adjacent pairs or windows from the end to the beginning; when cross-origin semantics are required, the sequence should be rotated or explicitly linearized first.

A total of 13 basic functions that can be called individually are listed below, including ambiguous base statistics, length and composition; when you need to output a more complete fixed feature vector at one time, please continue to view [ "All descriptor calculations" ](../api/features/05_all_descriptors.md#all-descriptors) later in this page.


### 14) All descriptor calculations (240 items)

**Full 240 descriptors**

Compute the compositional, statistical, physicochemical, and dinucleotide descriptors of DNA sequences at once to obtain 240 feature results in a fixed order.

`all_descriptors()` Use the fixed version `descriptor_schema_v1` and return the non-computable reason, calculation conditions and source information at the same time. The order of the fields will not change with the input; the first 180 items are calculated according to the input and their respective applicable fields. DNAKit does not have a built-in dinucleotide value table, so the last 60 items default to `None`.

**Shortest usage**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

The CLI outputs the same complete set of fields by default; if you need the old four sets of condensed results, you can use `--compact` explicitly:

```bash
dnakit describe ACGT
dnakit describe ACGT --compact
```

**API**

- `dnakit.descriptors.all_descriptors(value[required], ambiguity_policy[optional], conditions[optional], dinucleotide_property_table[optional])`: Returns all fields; the last 60 are calculated only if the table is explicitly provided.
- `dnakit.descriptors.load_dinucleotide_property_table(path[required])`: Bounded read and strict validation of user JSON table, record file SHA-256.
- `dnakit.descriptors.descriptor_schema_v1()`: Returns 240 immutable `DescriptorField`, each containing `index`, `name`, `category`, `unit`, `formula`, and `source`.
- `dnakit.descriptors.DESCRIPTOR_SCHEMA_V1`: constant form of the same fixed schema.
- `dnakit.descriptors.DESCRIPTOR_NAMES_V1`: Contains only ordered field names.

Result object structure:

| Properties | Meaning |
| ----------------------- | --------------------------------------------------------------- |
| `schema_version` | Fixed to `descriptor_schema_v1` |
| `sequence_id` | Keep record ID if input is `DNARecord`, otherwise `None` |
| `values` | 240 values strictly arranged by schema |
| `unavailable_reasons` | Contains only fields with value `None`, and each field has exactly one reason |
| `conditions` | IUPAC, Gap, k-mer, ORF, replicate, thermodynamic and parameter table conditions |
| `provenance` | DNAKit version, implementation type, and user table declaration name, version, origin, and SHA-256 |

**Unified calculation caliber**

- Only A/C/G/T counts in the canonical denominator; defaults to `ambiguity_policy="ignore"`, can also be set to `"error"`.
- k-mers, CpG/GpC, user dinucleotide tables, and other neighbor word statistics all allow overlap, but do not span IUPAC ambiguous symbols or explicit gaps.
- All scales returning `None` when the denominator is 0 are not replaced with fake `0`; the reason is written to `unavailable_reasons`.
- User table statistics `sd` use the population standard deviation; the table must provide a specified set of 15 attributes and a finite number of 16 DNA dinucleotides per set.
- frame 0 uses NCBI standard genetic code 1; the complete ORF uses `ATG` to start and `TAA/TAG/TGA` to stop, and scans the forward and reverse strands for a total of six reading frames.
- Thermodynamic default conditions are 37 °C, Na⁺ 0.05 M, chain concentration 250 nM; Wallace only 2–13 nt, SantaLucia NN only 2–60 nt.
- Molecular weight and ε260 are only for linear, no-gap, no-modification, A/C/G/T DNA; ε260 is a theoretical coefficient, not experimental A260.
- Currently LZ76 calculates sequences that are linear, have no gaps, no ambiguous symbols and do not exceed 10,000 nt; longer sequences return `None` to avoid runaway quadratic time complexity.

**Category and field range**

| Serial number | Category | Number of fields |
| -------- | ---------------------------------------- | -----: |
| 1–12 | Base Length and Data Quality | 12 |
| 13–28 | Chemical grouping of bases | 16 |
| 29–112 | 1/2/3-mer frequency | 84 |
| 113–128 | skew, CpG/GpC and Chargaff metrics | 16 |
| 129–148 | Entropy, complexity, homopolymers and repetition | 20 |
| 149–164 | Codons and six reading frame ORF | 16 |
| 165–180 | Molecular weight, ε260, Tm and NN Thermodynamics | 16 |
| 181–240 | 15 sets of user dinucleotide parameters × mean/sd/min/max | 60 |

**240 complete field table**

| # | Field | Category | Unit | Formula | Source |
| --: | -------------------------------------------------- | ---------------------------------- | --------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| 1 | `symbol_length` | `basic` | `nt` | `number of nucleotide symbols; explicit gaps excluded` | DNAKit descriptor_schema_v1 |
| 2 | `coordinate_span` | `basic` | `nt` | `symbol_length + sum(known gap lengths); undefined with an unknown gap` | DNAKit descriptor_schema_v1 |
| 3 | `canonical_base_count` | `basic` | `count` | `count(A,C,G,T)` | DNAKit descriptor_schema_v1 |
| 4 | `ambiguity_symbol_count` | `basic` | `count` | `symbol_length - canonical_base_count` | DNAKit descriptor_schema_v1 |
|   5 | `gap_object_count`                         | `basic`                 | `count`       | `number of explicit Gap objects`                                                           | DNAKit descriptor_schema_v1                                        |
|   6 | `known_gap_nt`                             | `basic`                 | `nt`          | `sum(length of known explicit gaps)`                                                       | DNAKit descriptor_schema_v1                                        |
|   7 | `unknown_gap_count`                        | `basic`                 | `count`       | `number of explicit gaps with unknown length`                                              | DNAKit descriptor_schema_v1                                        |
|   8 | `canonical_symbol_fraction`                | `basic`                 | `fraction`    | `canonical_base_count / symbol_length`                                                     | DNAKit descriptor_schema_v1                                        |
|   9 | `ambiguity_symbol_fraction`                | `basic`                 | `fraction`    | `ambiguity_symbol_count / symbol_length`                                                   | DNAKit descriptor_schema_v1                                        |
|  10 | `known_gap_fraction`                       | `basic`                 | `fraction`    | `known_gap_nt / coordinate_span`                                                           | DNAKit descriptor_schema_v1                                        |
|  11 | `canonical_run_count`                      | `basic`                 | `count`       | `number of uninterrupted A/C/G/T runs split by ambiguity or gaps`                          | DNAKit descriptor_schema_v1                                        |
|  12 | `longest_canonical_run_nt`                 | `basic`                 | `nt`          | `maximum uninterrupted A/C/G/T run length`                                                 | DNAKit descriptor_schema_v1                                        |
|  13 | `purine_count`                             | `composition`           | `count`       | `count(A)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  14 | `purine_fraction`                          | `composition`           | `fraction`    | `purine_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  15 | `pyrimidine_count`                         | `composition`           | `count`       | `count(C)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  16 | `pyrimidine_fraction`                      | `composition`           | `fraction`    | `pyrimidine_count / canonical_base_count`                                                  | DNAKit descriptor_schema_v1                                        |
|  17 | `amino_count`                              | `composition`           | `count`       | `count(A)+count(C)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  18 | `amino_fraction`                           | `composition`           | `fraction`    | `amino_count / canonical_base_count`                                                       | DNAKit descriptor_schema_v1                                        |
|  19 | `keto_count`                               | `composition`           | `count`       | `count(G)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  20 | `keto_fraction`                            | `composition`           | `fraction`    | `keto_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  21 | `weak_count`                               | `composition`           | `count`       | `count(A)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  22 | `weak_fraction`                            | `composition`           | `fraction`    | `weak_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  23 | `strong_count`                             | `composition`           | `count`       | `count(C)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  24 | `strong_fraction`                          | `composition`           | `fraction`    | `strong_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  25 | `purine_pyrimidine_skew`                   | `composition`           | `ratio`       | `(purine_count-pyrimidine_count)/(purine_count+pyrimidine_count)`                          | DNAKit descriptor_schema_v1                                        |
|  26 | `amino_keto_skew`                          | `composition`           | `ratio`       | `(amino_count-keto_count)/(amino_count+keto_count)`                                        | DNAKit descriptor_schema_v1                                        |
|  27 | `weak_strong_skew`                         | `composition`           | `ratio`       | `(weak_count-strong_count)/(weak_count+strong_count)`                                      | DNAKit descriptor_schema_v1                                        |
|  28 | `gc_at_ratio`                              | `composition`           | `ratio`       | `strong_count / weak_count`                                                                | DNAKit descriptor_schema_v1                                        |
|  29 | `k1_A_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(A) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  30 | `k1_C_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(C) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  31 | `k1_G_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(G) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  32 | `k1_T_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(T) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  33 | `k2_AA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  34 | `k2_AC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  35 | `k2_AG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  36 | `k2_AT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  37 | `k2_CA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  38 | `k2_CC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  39 | `k2_CG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  40 | `k2_CT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  41 | `k2_GA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  42 | `k2_GC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  43 | `k2_GG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  44 | `k2_GT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  45 | `k2_TA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  46 | `k2_TC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  47 | `k2_TG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  48 | `k2_TT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  49 | `k3_AAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  50 | `k3_AAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  51 | `k3_AAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  52 | `k3_AAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  53 | `k3_ACA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  54 | `k3_ACC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  55 | `k3_ACG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  56 | `k3_ACT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  57 | `k3_AGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  58 | `k3_AGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  59 | `k3_AGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  60 | `k3_AGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  61 | `k3_ATA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  62 | `k3_ATC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  63 | `k3_ATG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  64 | `k3_ATT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  65 | `k3_CAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  66 | `k3_CAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  67 | `k3_CAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  68 | `k3_CAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  69 | `k3_CCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  70 | `k3_CCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  71 | `k3_CCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  72 | `k3_CCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  73 | `k3_CGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  74 | `k3_CGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  75 | `k3_CGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  76 | `k3_CGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  77 | `k3_CTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  78 | `k3_CTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  79 | `k3_CTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  80 | `k3_CTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  81 | `k3_GAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  82 | `k3_GAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  83 | `k3_GAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  84 | `k3_GAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  85 | `k3_GCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  86 | `k3_GCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  87 | `k3_GCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  88 | `k3_GCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  89 | `k3_GGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  90 | `k3_GGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  91 | `k3_GGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  92 | `k3_GGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  93 | `k3_GTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  94 | `k3_GTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  95 | `k3_GTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  96 | `k3_GTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  97 | `k3_TAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  98 | `k3_TAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  99 | `k3_TAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 100 | `k3_TAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 101 | `k3_TCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 102 | `k3_TCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 103 | `k3_TCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 104 | `k3_TCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 105 | `k3_TGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 106 | `k3_TGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 107 | `k3_TGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 108 | `k3_TGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 109 | `k3_TTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 110 | `k3_TTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 111 | `k3_TTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 112 | `k3_TTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 113 | `gc_skew`                                  | `skew_cpg`              | `ratio`       | `(count(G)-count(C))/(count(G)+count(C))`                                                  | DNAKit descriptor_schema_v1                                        |
| 114 | `at_skew`                                  | `skew_cpg`              | `ratio`       | `(count(A)-count(T))/(count(A)+count(T))`                                                  | DNAKit descriptor_schema_v1                                        |
| 115 | `cpg_count`                                | `skew_cpg`              | `count`       | `overlapping count(CG)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 116 | `cpg_density`                              | `skew_cpg`              | `fraction`    | `count(CG) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 117 | `cpg_observed_expected`                    | `skew_cpg`              | `ratio`       | `count(CG)*canonical_base_count/(count(C)*count(G))`                                       | DNAKit descriptor_schema_v1                                        |
| 118 | `gpc_count`                                | `skew_cpg`              | `count`       | `overlapping count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 119 | `gpc_density`                              | `skew_cpg`              | `fraction`    | `count(GC) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 120 | `cpg_gpc_ratio`                            | `skew_cpg`              | `ratio`       | `count(CG) / count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 121 | `cumulative_gc_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 122 | `cumulative_gc_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 123 | `cumulative_gc_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_gc_skew_max - cumulative_gc_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 124 | `cumulative_at_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 125 | `cumulative_at_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 126 | `cumulative_at_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_at_skew_max - cumulative_at_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 127 | `dinucleotide_rc_total_variation`          | `skew_cpg`              | `fraction`    | `0.5*sum_xy(abs(f_xy-f_reverse_complement(xy)))`                                           | DNAKit descriptor_schema_v1                                        |
| 128 | `mono_chargaff_l1_distance`                | `skew_cpg`              | `fraction`    | `abs(f_A-f_T)+abs(f_C-f_G)`                                                                | DNAKit descriptor_schema_v1                                        |
| 129 | `shannon_entropy_k1_bits`                  | `complexity`            | `bits`        | `-sum(p(1-mer)*log2(p(1-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 130 | `shannon_entropy_k2_bits`                  | `complexity`            | `bits`        | `-sum(p(2-mer)*log2(p(2-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 131 | `shannon_entropy_k3_bits`                  | `complexity`            | `bits`        | `-sum(p(3-mer)*log2(p(3-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 132 | `normalized_entropy_k1`                    | `complexity`            | `fraction`    | `shannon_entropy_k1_bits / log2(4**1)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 133 | `normalized_entropy_k2`                    | `complexity`            | `fraction`    | `shannon_entropy_k2_bits / log2(4**2)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 134 | `normalized_entropy_k3`                    | `complexity`            | `fraction`    | `shannon_entropy_k3_bits / log2(4**3)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 135 | `linguistic_complexity_k2`                 | `complexity`            | `fraction`    | `unique 2-mers / min(4**2, valid 2-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 136 | `linguistic_complexity_k3`                 | `complexity`            | `fraction`    | `unique 3-mers / min(4**3, valid 3-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 137 | `linguistic_complexity_k4`                 | `complexity`            | `fraction`    | `unique 4-mers / min(4**4, valid 4-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 138 | `linguistic_complexity_k5`                 | `complexity`            | `fraction`    | `unique 5-mers / min(4**5, valid 5-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 139 | `linguistic_complexity_k6`                 | `complexity`            | `fraction`    | `unique 6-mers / min(4**6, valid 6-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 140 | `linguistic_complexity_product_k1_k6`      | `complexity`            | `fraction`    | `product of defined linguistic_complexity_k values for k=1..6`                             | Observed/possible k-word linguistic complexity                     |
| 141 | `lz76_complexity`                          | `complexity`            | `count`       | `number of phrases in exhaustive Lempel-Ziv 1976 parsing`                                  | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 142 | `normalized_lz76_complexity`               | `complexity`            | `ratio`       | `lz76_complexity*log_base4(canonical_base_count)/canonical_base_count`                     | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 143 | `longest_homopolymer_nt`                   | `complexity`            | `nt`          | `max canonical homopolymer run`                                                            | DNAKit descriptor_schema_v1                                        |
| 144 | `longest_homopolymer_a_nt`                 | `complexity`            | `nt`          | `max homopolymer run of A`                                                                 | DNAKit descriptor_schema_v1                                        |
| 145 | `longest_homopolymer_c_nt`                 | `complexity`            | `nt`          | `max homopolymer run of C`                                                                 | DNAKit descriptor_schema_v1                                        |
| 146 | `longest_homopolymer_g_nt`                 | `complexity`            | `nt`          | `max homopolymer run of G`                                                                 | DNAKit descriptor_schema_v1                                        |
| 147 | `longest_homopolymer_t_nt`                 | `complexity`            | `nt`          | `max homopolymer run of T`                                                                 | DNAKit descriptor_schema_v1                                        |
| 148 | `exact_tandem_repeat_coverage_fraction`    | `complexity`            | `fraction`    | `union bases covered by exact tandem repeats / canonical_base_count`                       | DNAKit exact tandem repeat scanner; units 1..20; minimum repeats 2 |
| 149 | `frame0_codon_count`                       | `coding`                | `count`       | `number of valid forward frame-0 codons`                                                   | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 150 | `frame0_unique_codon_count`                | `coding`                | `count`       | `number of distinct forward frame-0 codons`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 151 | `frame0_start_codon_count`                 | `coding`                | `count`       | `count(ATG) in forward frame 0`                                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 152 | `frame0_stop_codon_count`                  | `coding`                | `count`       | `count(TAA,TAG,TGA) in forward frame 0`                                                    | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 153 | `frame0_start_codon_fraction`              | `coding`                | `fraction`    | `frame0_start_codon_count / frame0_codon_count`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 154 | `frame0_stop_codon_fraction`               | `coding`                | `fraction`    | `frame0_stop_codon_count / frame0_codon_count`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 155 | `frame0_codon_entropy_bits`                | `coding`                | `bits`        | `-sum(frame0 codon frequency*log2(frequency))`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 156 | `frame0_effective_number_of_codons`        | `coding`                | `count`       | `2**frame0_codon_entropy_bits`                                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 157 | `frame0_gc1_fraction`                      | `coding`                | `fraction`    | `GC bases at position 1 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 158 | `frame0_gc2_fraction`                      | `coding`                | `fraction`    | `GC bases at position 2 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 159 | `frame0_gc3_fraction`                      | `coding`                | `fraction`    | `GC bases at position 3 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 160 | `six_frame_complete_orf_count`             | `coding`                | `count`       | `complete start-to-next-stop ORFs across three frames on both strands`                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 161 | `six_frame_forward_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three forward frames`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 162 | `six_frame_reverse_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three reverse-complement frames`                                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 163 | `six_frame_longest_complete_orf_nt`        | `coding`                | `nt`          | `maximum complete six-frame ORF length including terminal stop`                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 164 | `six_frame_complete_orf_coverage_fraction` | `coding`                | `fraction`    | `union of symbol positions covered by complete six-frame ORFs / canonical_base_count`      | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 165 | `mw_ss_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime OH`                                       | DNAKit anhydrous DNA residue mass table v1                         |
| 166 | `mw_ss_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime phosphate`                                | DNAKit anhydrous DNA residue mass table v1                         |
| 167 | `mw_ds_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement with 5-prime OH`              | DNAKit anhydrous DNA residue mass table v1                         |
| 168 | `mw_ds_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement; both 5-prime phosphorylated` | DNAKit anhydrous DNA residue mass table v1                         |
| 169 | `epsilon260_ss_m_inverse_cm_inverse`       | `physicochemical`       | `M^-1 cm^-1`  | `nearest-neighbor epsilon260 pair sum minus internal-base sum`                             | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 170 | `nmol_per_a260_1ml_1cm`                    | `physicochemical`       | `nmol`        | `1e6 / epsilon260 for A260=1, volume=1 mL, path=1 cm`                                      | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 171 | `ug_per_a260_1ml_1cm`                      | `physicochemical`       | `ug`          | `1000*mw_ss_oh_da/epsilon260 for A260=1, volume=1 mL, path=1 cm`                           | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 172 | `tm_wallace_c`                             | `physicochemical`       | `degree C`    | `2*(A+T)+4*(G+C)`                                                                          | Wallace short-oligo 2AT+4GC rule                                   |
| 173 | `stacking_delta_h_kcal_per_mol`            | `physicochemical`       | `kcal/mol`    | `sum SantaLucia nearest-neighbor stacking delta H`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 174 | `stacking_delta_s_cal_per_k_mol_k`         | `physicochemical`       | `cal/(K mol)` | `sum SantaLucia nearest-neighbor stacking delta S`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 175 | `stacking_delta_g37_kcal_per_mol`          | `physicochemical`       | `kcal/mol`    | `stacking delta H - 310.15*stacking delta S/1000`                                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 176 | `nn_delta_h_kcal_per_mol`                  | `physicochemical`       | `kcal/mol`    | `SantaLucia complete-duplex delta H with initiation and symmetry`                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 177 | `nn_delta_s_cal_per_mol_k`                 | `physicochemical`       | `cal/(K mol)` | `SantaLucia complete-duplex delta S with initiation, symmetry, and salt`                   | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 178 | `nn_delta_g37_kcal_per_mol`                | `physicochemical`       | `kcal/mol`    | `complete-duplex delta H - 310.15*delta S/1000`                                            | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 179 | `nn_tm_c`                                  | `physicochemical`       | `degree C`    | `SantaLucia concentration- and sodium-adjusted Tm`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 180 | `self_complementary`                       | `physicochemical`       | `boolean`     | `sequence == reverse_complement(sequence)`                                                 | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 181 | `diprodb_twist_mean`                       | `dinucleotide_property` | `degree`      | `population mean of Twist values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 182 | `diprodb_twist_sd`                         | `dinucleotide_property` | `degree`      | `population sd of Twist values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 183 | `diprodb_twist_min`                        | `dinucleotide_property` | `degree`      | `population min of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 184 | `diprodb_twist_max`                        | `dinucleotide_property` | `degree`      | `population max of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 185 | `diprodb_tilt_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Tilt values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 186 | `diprodb_tilt_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Tilt values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 187 | `diprodb_tilt_min`                         | `dinucleotide_property` | `degree`      | `population min of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 188 | `diprodb_tilt_max`                         | `dinucleotide_property` | `degree`      | `population max of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 189 | `diprodb_roll_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Roll values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 190 | `diprodb_roll_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Roll values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 191 | `diprodb_roll_min`                         | `dinucleotide_property` | `degree`      | `population min of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 192 | `diprodb_roll_max`                         | `dinucleotide_property` | `degree`      | `population max of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 193 | `diprodb_shift_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Shift values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 194 | `diprodb_shift_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Shift values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 195 | `diprodb_shift_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 196 | `diprodb_shift_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 197 | `diprodb_slide_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Slide values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 198 | `diprodb_slide_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Slide values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 199 | `diprodb_slide_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 200 | `diprodb_slide_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 201 | `diprodb_rise_mean`                        | `dinucleotide_property` | `angstrom`    | `population mean of Rise values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 202 | `diprodb_rise_sd`                          | `dinucleotide_property` | `angstrom`    | `population sd of Rise values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 203 | `diprodb_rise_min`                         | `dinucleotide_property` | `angstrom`    | `population min of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 204 | `diprodb_rise_max`                         | `dinucleotide_property` | `angstrom`    | `population max of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 205 | `diprodb_bend_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Bend values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 206 | `diprodb_bend_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Bend values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 207 | `diprodb_bend_min`                         | `dinucleotide_property` | `degree`      | `population min of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 208 | `diprodb_bend_max`                         | `dinucleotide_property` | `degree`      | `population max of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 209 | `diprodb_inclination_mean`                 | `dinucleotide_property` | `degree`      | `population mean of Inclination values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 210 | `diprodb_inclination_sd`                   | `dinucleotide_property` | `degree`      | `population sd of Inclination values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 211 | `diprodb_inclination_min`                  | `dinucleotide_property` | `degree`      | `population min of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 212 | `diprodb_inclination_max`                  | `dinucleotide_property` | `degree`      | `population max of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 213 | `diprodb_direction_mean`                   | `dinucleotide_property` | `degree`      | `population mean of Direction values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 214 | `diprodb_direction_sd`                     | `dinucleotide_property` | `degree`      | `population sd of Direction values over valid overlapping dinucleotides`                   | Caller-supplied table; DNAKit bundles no numerical values          |
| 215 | `diprodb_direction_min`                    | `dinucleotide_property` | `degree`      | `population min of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 216 | `diprodb_direction_max`                    | `dinucleotide_property` | `degree`      | `population max of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 217 | `diprodb_propeller_twist_mean`             | `dinucleotide_property` | `degree`      | `population mean of Propeller twist values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 218 | `diprodb_propeller_twist_sd`               | `dinucleotide_property` | `degree`      | `population sd of Propeller twist values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 219 | `diprodb_propeller_twist_min`              | `dinucleotide_property` | `degree`      | `population min of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 220 | `diprodb_propeller_twist_max`              | `dinucleotide_property` | `degree`      | `population max of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 221 | `diprodb_major_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Major groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 222 | `diprodb_major_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Major groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 223 | `diprodb_major_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 224 | `diprodb_major_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 225 | `diprodb_minor_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Minor groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 226 | `diprodb_minor_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Minor groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 227 | `diprodb_minor_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 228 | `diprodb_minor_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 229 | `diprodb_persistence_length_mean`          | `dinucleotide_property` | `nanometer`   | `population mean of Persistence length values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 230 | `diprodb_persistence_length_sd`            | `dinucleotide_property` | `nanometer`   | `population sd of Persistence length values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 231 | `diprodb_persistence_length_min`           | `dinucleotide_property` | `nanometer`   | `population min of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 232 | `diprodb_persistence_length_max`           | `dinucleotide_property` | `nanometer`   | `population max of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 233 | `diprodb_stacking_energy_mean`             | `dinucleotide_property` | `kcal/mol`    | `population mean of Stacking energy values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 234 | `diprodb_stacking_energy_sd`               | `dinucleotide_property` | `kcal/mol`    | `population sd of Stacking energy values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 235 | `diprodb_stacking_energy_min`              | `dinucleotide_property` | `kcal/mol`    | `population min of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 236 | `diprodb_stacking_energy_max` | `dinucleotide_property` | `kcal/mol` | `population max of Stacking energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 237 | `diprodb_free_energy_mean` | `dinucleotide_property` | `kcal/mol` | `population mean of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 238 | `diprodb_free_energy_sd` | `dinucleotide_property` | `kcal/mol` | `population sd of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 239 | `diprodb_free_energy_min` | `dinucleotide_property` | `kcal/mol` | `population min of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 240 | `diprodb_free_energy_max` | `dinucleotide_property` | `kcal/mol` | `population max of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |

Methods, papers, databases and URLs have been moved to [Acknowledgments and Primary Sources](../acknowledgements.md#methods-and-references); permissions and user responsibilities can be found in [Third Party Statements](../acknowledgements.md#third-party-notices).

## Stable hash and K-fold splits {#dataset-splits}

### 1) Random and stable hash partitioning {#1-data-012-random-and-stable-hash-partitioning}

**Function:** Partition records into named subsets such as `train / valid / test` by ratios, or generate a fixed K-fold partition and select one fold as `valid` and the rest as `train`. Return the subsets and assignment manifest.


**API:** `dnakit.datasets.split(records[required], config[required])`; use `dnakit.datasets.SplitConfig` with `method="random"`, `method="hash"`, or `method="kfold"`.

**Input:** Records are required. Ratio splitting uses fractions summing to 1, defaulting to `train=0.8, valid=0.1, test=0.1`. K-fold uses `n_splits` for the fold count and `fold_id` for the validation fold. Defaults are `seed=42` and `shuffle=True`. Hash splitting requires stable, unique record IDs.

The default seed for `SplitConfig` and `dnakit split` changed from `0` to `42`. Pass `seed=0` explicitly to reproduce earlier ratio splits that omitted the seed.

**Sample code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

**Example results:**

```text
{'train': 8, 'test': 2}
```

**Order-independent `hash` partitioning**

The `hash` pattern uses versioned SHA-256 to calculate the stable sort key for `seed + record.id` when `shuffle=True` and the stable endianness for `record.id` when `shuffle=False` and allocates the records proportionally. Therefore, even if the input order of the same batch of records changes, the classification will still be exactly the same when viewed by `record.id`; it does not use Python's built-in process randomization `hash()`.

**Sample code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
config = SplitConfig(
    method="hash",
    ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
    seed=17,
    preserve_order=False,
)
first = split(records, config=config)
second = split(list(reversed(records)), config=config)
first_by_id = {item.record_id: item.split for item in first.assignments}
second_by_id = {item.record_id: item.split for item in second.assignments}
print(first_by_id == second_by_id)
print(first.get("train").ids == second.get("train").ids)
```

**Example results:**

```text
True
True
```

#### Fixed K-fold partitioning: train / valid

Set `method="kfold"` to use [scikit-learn KFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html) semantics:

```python
KFold(n_splits=10, shuffle=True, random_state=42)
```

**One random seed → one fixed K-fold partition → K validation folds selected in turn.** Keep the records and their input order, `seed`, `shuffle`, and `n_splits` unchanged; change only `fold_id`. For 10 folds, use `fold_id=0..9`: each sample serves as validation exactly once. Do not use `seed + fold_id`. Each call reconstructs the same partition and returns the selected pair.

| `SplitConfig` argument | Default | Meaning |
| --- | --- | --- |
| `method` | `"random"` | Select `"kfold"` for K-fold; `"random"` / `"hash"` retain ratio splitting |
| `n_splits` | `10` in K-fold mode | Fold count, such as 5 or 10; requires `2 <= n_splits <= sample count` |
| `fold_id` | `0` in K-fold mode | Validation fold index, from `0` through `n_splits - 1` |
| `seed` | `42` | Fixed seed for the whole partition; K-fold requires `0 <= seed < 2**32` |
| `shuffle` | `True` | Shuffle before partitioning; `False` gives consecutive validation blocks and seed has no effect |
| `ratios` | `0.8 / 0.1 / 0.1` in ratio mode | Omit for K-fold; target train/valid ratios are derived from the fold count |

K-fold always returns subsets named `train` and `valid`. Records within each subset retain input order, matching scikit-learn, even with `preserve_order=False`. If the sample count is not divisible by K, the first `sample count % K` validation folds contain one extra sample; fold sizes differ by at most 1. Results and exported manifests record `seed`, `shuffle`, `n_splits`, `fold_id`, and each record's assignment.

This mode requires scikit-learn: install `python -m pip install -e ".[bio-analysis]"`. Environments with `.[ml]` or `.[training]` already include it. Ratio splitting does not require this dependency.

**Example: hold out a fixed test set, then apply 10-fold splitting to train.**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("ACGT"), f"r{i}") for i in range(100)]
seed = 42
n_splits = 10  # Change to 5 for five folds.

## Split off the fixed test set once.
outer = split(
    records,
    config=SplitConfig(method="random", ratios={"train": 0.8, "test": 0.2}, seed=seed),
)
training_pool = outer.get("train")
test = outer.get("test")

## Only fold_id changes; never pass the test set into K-fold.
for fold_id in range(n_splits):
    fold = split(
        training_pool,
        config=SplitConfig(
            method="kfold", n_splits=n_splits, fold_id=fold_id,
            shuffle=True, seed=seed,
        ),
    )
    train = fold.get("train")
    valid = fold.get("valid")
    print(f"fold={fold_id}: train={len(train)}, valid={len(valid)}, test={len(test)}")
    # Create and train a fresh model here, then evaluate on valid.
```

Folds 0 through 9 each contain `train=72, valid=8, test=20`; validation samples rotate while the 20 test samples stay fixed. If you already have training data, use it directly as `training_pool`. `n_splits` sets the total fold count, and `fold_id` selects one fold. Iterate all folds for a complete K-fold procedure; `split()` partitions data and does not train models or aggregate metrics.

To use `train_custom_property_model()`, write the current train/valid subsets and fixed test set using the training API's supported split labels, and map that column with `DataMapping.split_field`. Values such as `fold_0` are not training split labels. See [custom model training](../api/features/23_prediction_custom_training.md).

**CLI:** Pass the existing training pool to export the selected train/valid pair and manifest:

```bash
dnakit split train.fasta runs/fold_0 --method kfold --n-splits 10 --fold-id 0 --seed 42
```

Keep the input and seed fixed for subsequent runs, changing `--fold-id` and the output directory. Use `--no-shuffle` to disable shuffling. Agent/MCP calls to `dnakit.datasets.split` accept the same `config`. JSON/YAML workflow `split.params` also accepts `method`, `n_splits`, and `fold_id`; its seed comes from the workflow's global `seed`, which should be set to `42` for this example.

## Downloads: formats, progress, retries and species names {#downloads}

- **Agent arguments:** `query` is an assembly identifier or taxon name, not tool-search keywords. For human hg38, pass `query="hg38"`, not `"human hg38 reference genome FASTA"`. The tool schema includes this constraint. The downloader's `api_key` is an NCBI key, not a DeepSeek key.


Download files, reference genomes, sequences, and annotations from URLs, species names, database identifiers, or query results.

This page includes parameters and examples for NCBI/ENA sequences, URL downloads and ENCODE batch downloads. See [data queries](../api/features/22_database_query_encode.md#encode) for experiments and [dataset preparation](../api/features/10_datasets.md#encode-collection) for tracking, citations and provenance.

### Download sequences by accession

`dnakit.download.ncbi.sequence(accession, output_path, *, database="nuccore", format="fasta",
seq_start=None, seq_stop=None, strand=1, config=None, search_config=None)` uses NCBI
EFetch directly. `dnakit.download.ena.sequence(accession, output_path, *, format="fasta",
config=None, search_config=None)` uses ENA Browser. Neither requires ToolUniverse.

```python
import dnakit

result = dnakit.download.ncbi.sequence("NC_000913.3", "ecoli_k12.fasta", format="fasta")
print(result.files[0].path)
print(result.metadata["resolved_accession"], result.metadata["sequence_length"])
print(result.manifest_path)
```

| Provider | Native formats with base DNAKit | Optional conversion with `dnakit[bio-formats]` |
| --- | --- | --- |
| NCBI nuccore (DNA/RNA) | fasta, genbank | embl, converted from GenBank with Biopython |
| NCBI protein | fasta, genbank (GenPept) | EMBL protein output is unsupported |
| ENA INSDC nucleotide records | fasta, embl | genbank, converted from EMBL with Biopython |

NCBI accepts sequence accessions, including RefSeq NZ_. ENA requires INSDC
accessions such as U00096.3. Conversion does not imply cross-database identity.
Use existing `genome_package()` for GCA_/GCF_ assembly accessions.

The returned `DatasetDownloadResult` identifies the sequence file and sibling
`.manifest.json`: requested/resolved accession, sequence length, redacted URL,
source/output format, conversion flag and checksums. Responses are checked for
accession/version, format and sequence length before installation. Hashes record
local integrity; they are not independent provider signatures. Existing outputs
are protected unless `DownloadConfig(overwrite=True)` is supplied.

Optional seq_start/seq_stop must be supplied together and are **1-based inclusive**;
strand=2 requests the reverse complement. These options apply to NCBI nucleotide
records; proteins use full records. Ensembl coordinate downloads retain their
0-based half-open convention.

DownloadConfig controls files, overwrite and byte limits. SearchConfig controls
timeout, response size (20 MB default, 100 MB maximum), and NCBI api_key/email.
An explicit search_config takes precedence for credentials; otherwise
DownloadConfig.api_key is used. Credentials are redacted in persisted URLs.
Use genome packages for larger assemblies. Install optional conversion support
with `python -m pip install 'dnakit[bio-formats]'`.

MCP names and flat Python aliases are `dnakit.download.ncbi_sequence` and
`dnakit.download.ena_sequence`. The existing `dnakit.download.sequence` continues
to retrieve Ensembl coordinate sequences.

#### Download progress

Network downloads print a start message, progress bar, and completion path to stderr by default.
Reference genomes also show extraction progress. If the server omits the file size, the display
shows downloaded bytes without estimating a percentage. Completion is reported only after the
download, requested checksums, and output installation succeed. Disable the display with
`download_genome(..., show_progress=False)`, or `DownloadConfig(show_progress=False)` for network
download APIs in `dnakit.download` that accept `config`. An explicit `progress` callback (`query_progress` for sequence
downloads) takes over the display to avoid duplicate progress bars. NCBI packages, ENA, ENCODE,
UCSC, and generic file downloads share this behavior.

#### Download failures and retries

Genome packages, generic HTTPS file downloads and shared remote GET queries make at most three
attempts for timeouts, interrupted connections and HTTP 408/429/500/502/503/504, with delays of
0.25 and 0.5 seconds. HTTP 404, invalid parameters, checksum failures and invalid ZIP archives
are not retried. POST requests are not automatically repeated.

Retries restart the staged transfer; persistent HTTP Range resume after process exit is not
implemented. When `Content-Length` is supplied, the received byte count must match. Genome
packages also require a matching assembly accession and FASTA checksum. Files are installed
only after complete download and validation. Overwrite conflicts retain `FileExistsError`, and directory permission failures retain their filesystem exceptions;
and failed transfers do not replace complete existing outputs. Exhausted transient failures
set `context["attempts"]` to 3 and guidance to retry after recovery.

**Common species names**

This section lists common NCBI species names that can be copied directly into
`download_genome(query, ...)` or a `dnakit download-genome` command.

- **API**: `dnakit.references.download_genome(query[required], output_dir[required], overwrite[optional], keep_package[optional], timeout[optional], api_key[optional], api_base_url[optional], chunk_size[optional], max_download_bytes[optional], progress[optional])`, `dnakit download-genome query[required] output_dir[required] --overwrite[optional] --keep-package[optional] --api-key[optional] --progress/--no-progress[optional]`.

**Humans and common model animals**

| Common name | `query` that can be copied directly | Current reference assembly accession |
| -------- | -------------------------- | -------------------------- |
| Human | `Homo sapiens` | `GCF_000001405.40` |
| Mouse | `Mus musculus` | `GCF_000001635.27` |
| Rat | `Rattus norvegicus` | `GCF_036323735.1` |
| Pig | `Sus scrofa` | `GCF_054392235.1` |
| Cow | `Bos taurus` | `GCF_002263795.3` |
| Sheep | `Ovis aries` | `GCF_016772045.2` |
| Goat | `Capra hircus` | `GCF_001704415.2` |
| Horse | `Equus caballus` | `GCF_041296265.1` |
| Dog | `Canis lupus familiaris` | `GCF_011100685.1` |
| Macaque | `Macaca mulatta` | `GCF_049350105.2` |
| Chicken | `Gallus gallus` | `GCF_016699485.2` |
| Zebrafish | `Danio rerio` | `GCF_049306965.2` |
| Xenopus | `Xenopus tropicalis` | `GCF_000004195.4` |

**Common crops and plants**

| Common name | `query` that can be copied directly | Current reference assembly accession |
| -------------------------- | ------------------------------- | ----------------------- |
| Arabidopsis | `Arabidopsis thaliana` | `GCF_000001735.4` |
| Rice | `Oryza sativa Japonica Group` | `GCF_034140825.1` |
| Corn | `Zea mays` | `GCF_902167145.1` |
| Wheat | `Triticum aestivum` | `GCF_018294505.1` |
| Soybeans | `Glycine max` | `GCF_000004515.6` |
| Tomato | `Solanum lycopersicum` | `GCF_036512215.1` |
| Potatoes | `Solanum tuberosum` | `GCF_000226075.1` |
| Upland cotton | `Gossypium hirsutum` | `GCF_007990345.1` |
| Tobacco | `Nicotiana tabacum` | `GCF_000715075.1` |
| Populus trichocarpa | `Populus trichocarpa` | `GCF_000002775.5` |
| Arabica coffee | `Coffea arabica` | `GCF_036785885.1` |

**Other commonly used model creatures**

| Common name | `query` that can be copied directly | Current reference assembly accession |
| --------------------- | ---------------------------------- | ----------------------- |
| Drosophila | `Drosophila melanogaster` | `GCF_000001215.4` |
| Caenorhabditis elegans | `Caenorhabditis elegans` | `GCF_000002985.6` |
| Bombyx mori | `Bombyx mori` | `GCF_030269925.1` |
| Bee | `Apis mellifera` | `GCF_003254395.2` |
| Saccharomyces cerevisiae | `Saccharomyces cerevisiae S288C` | `GCF_000146045.2` |


### 22) Download from a URL {#22-dbd-022-download-from-a-url}

`download_url(url, output_dir)` needs only a direct HTTPS file URL and a destination directory.
Missing directories are created automatically. The default filename comes from the last URL
path component, with percent escapes decoded and query parameters excluded.

```python
from dnakit.download import download_url

result = download_url(
    "https://example.org/files/data.fa.gz",  # Replace with your direct file URL
    "downloads/my_data",
)
print(result.path)       # Absolute path to the downloaded file
print(result.byte_size)  # File size in bytes
```

Pass `filename="reference.fa.gz"` to choose a different name; this parameter is required if
the URL has no filename. Existing files raise `FileExistsError` by default. To replace one,
pass `config=DownloadConfig(overwrite=True)` and import `DownloadConfig` from `dnakit.download`.
`DownloadConfig` also controls the timeout, file size limit, and progress display.
Use `expected_md5` / `expected_sha256` to verify the contents, or `headers` for request headers.

Any file type is supported and response bytes are saved unchanged, so `.gz` files stay compressed.
This function accepts direct file URLs; it does not resolve web or sharing pages or extract archives.
Agents and MCP can call `dnakit.download.download_url` with the same two required arguments,
`url` and `output_dir`. For multiple files with a provenance manifest, use
`dataset([RemoteFile(url1), RemoteFile(url2)], output_dir)`.


### 23) Search and download a sequence {#downloads-sequence-download-example}

This uses the extended existing `entrez` plus the new `accessions` and NCBI `sequence` functions. It queries a remote service and writes a file in the current directory.

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(accessions[0], "ecoli_k12.fasta", format="fasta")
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("No sequence matches; review the query.")
```

For ENA use an INSDC accession, for example `dnakit.download.ena.sequence("U00096.3", "ecoli_ena.fasta")`. Download functions save files and manifests, whereas database sequence wrappers return provider records. See [formats, dependencies and coordinates](../api/features/15_download.md#accession-sequences).

## Gene annotation: installation, coordinates and reference selection {#gene-annotation}

`dnakit.annotation.annotate_genes()` maps genomic DNA to a reference and transfers named gene intervals from its GFF3. Supply a species name or NCBI TaxID; `assembly` is optional. On the first request without an assembly, DNAKit selects the species' current RefSeq reference assembly. The exact downloaded assembly and annotation are then pinned in the project cache.

This transfers existing reference annotations. It does not perform de novo gene prediction, name novel genes, or annotate spliced cDNA/transcripts. Missing gene names remain null; stable IDs or locus tags remain available.

### Installation and command line

DNAKit and the external minimap2 executable are required. Importing the module does not download references or run external programs.

Install the working copy:

```bash
python -m pip install -e .
```

Install the external mapper:

```bash
conda install -c conda-forge -c bioconda minimap2
```

Supply only a species to select and cache the default assembly automatically. `query.fasta` should contain one genomic DNA sequence:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

Specify an exact assembly:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --assembly GCF_000001405.40 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

Run offline with a matched local reference/annotation pair:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --reference-fasta reference.fna \
  --annotation-gff annotation.gff3 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json \
  --no-download
```

Use `--minimap2 /path/to/minimap2` for an explicit executable, `--threads 4` for concurrency, and `--overwrite` to replace an existing output. FASTA/GFF3 may be gzip-compressed. Each invocation processes one query; Python callers can loop over records and reuse the cache.

### Python

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",  # NCBI TaxIDs such as 9606 are also accepted
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(cache_dir=".dnakit/gene_annotation", threads=4),
)
print(result.reference["assembly"])
for gene in result.gene_records:
    print(gene.to_dict())
for segment in result.segments:
    print(segment.to_dict())
```

`assembly` accepts versioned `GCA_/GCF_` accessions, existing fixed aliases, or unambiguous species-scoped assembly names. Explicit old accessions are not replaced with a newer version. Multiple matches (including paired GenBank/RefSeq assemblies sharing a name) require an exact accession. Missing references or gene annotations produce actionable errors.

### Results

- `reference`: actual assembly accession and name, annotation metadata, file SHA256 checksums, source, and selection method. Automatic downloads validate the species/assembly relationship; strains may match their parent species.
- `gene_records`: overlapping gene instances, names/IDs, strand, reference intervals, and projected query parts. A missing `gene_name` is `null`; `display_name` falls back to locus tag/ID. `mapping_id` preserves mapping identity at repeated loci.
- `segments`: a disjoint partition of the whole query into `gene`, `intergenic`, and `unknown`. A segment can refer to several overlapping genes.
- `alignments`: mapping coordinates, strand, MAPQ, identity and acceptance status.

Python `Interval` objects are **0-based half-open**. JSON from `to_dict()` and the CLI is **1-based closed** throughout. For example, input positions 1–12000 are `Interval(0, 12000)` internally. `gene_records.start/end` is an envelope in the input; `parts` describes the confidently mapped bases and may have holes. Reverse strands and indels are projected using CIGAR blocks. Circular origin-spanning genes are supported when the reference is explicitly marked circular in GFF3.

For example, genes A:1–12000, B:18001–45000 and C:43000–52000 produce three `gene_records` and the following `segments`:

```text
1–12000      gene         A
12001–18000  intergenic
18001–42999  gene         B
43000–45000  gene         B + C
45001–52000  gene         C
```

### Default assembly and reliability

When `assembly` is omitted, DNAKit queries current RefSeq reference assemblies for the species and its strains. Missing references, multiple candidates, species mismatches, or absent GFF3 gene annotations produce explicit errors; an arbitrary alternative assembly is never selected. Supply local files for species or assemblies not yet catalogued.

Identical species/assembly requests reuse a pinned reference without querying NCBI. SHA256 protects the cached FASTA, GFF3 and minimap2 index; indexes are keyed by reference content, preset and mapper version. Use a new `cache_dir` to refresh reference/annotation selection. `--no-download` permits only pinned references or local file pairs.

Downloaded assemblies are checked against the requested taxonomy (including strain-to-species ancestry). Offline local species are caller-declared (`species_validated=false`). Local assembly declarations, sequence IDs and bounds are checked when available; callers remain responsible for providing a matching pair. GFF3 needs gene-level rows such as gene, pseudogene, or immune gene segments, not only transcripts/CDS. Contig identifiers must match exactly; `chr1` and accession aliases are not guessed.

The default minimap2 preset is `asm5` for long, accurate genomic DNA; adjust `preset` for your data. Defaults accept primary mappings with MAPQ ≥20, identity ≥0.90 and ≥100 aligned bases; MAPQ=255 is unknown. These are configurable filters, not calibrated accuracy guarantees.

`intergenic` means reliably mapped bases outside genes in the selected annotation. Unmapped/low-confidence regions, insertions, overlapping trusted primary mappings and contigs missing annotation coverage information are `unknown`. Absence of annotation is not evidence of biological absence of a gene.

Downloading a large reference and building a single-part index can require substantial storage/RAM. Configuration bounds reference/query length, downloads, GFF records, index/alignment size, and external process duration. A single-batch index avoids unreliable MAPQ from multi-batch indexing.

### Validation and sources

Run coordinate, default assembly resolution, cache integrity and CLI tests; installing minimap2 also enables the real-process test:

```bash
PYTHONPATH=src python -m pytest -q others/tests/unit/test_gene_annotation.py
```

Sources: [NCBI reference selection](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genome-processing/reference-selection/), [genome package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/), [GFF3 conventions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/), [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html).

This version adds multi-species gene annotation with optional assembly selection, overlapping intervals and reproducible reference caches.

## Optional backend dependencies and return values {#optional-bioinformatics}

The `annotation`, `comparative` and `molbio` APIs return `ProviderResult`. The `dnakit.databases` functions for [variant annotation](../api/features/35_variant_annotation.md#variant-annotation), [selection pressure](../api/features/30_phylogeny.md#selection), and [enrichment](../api/features/37_functional_enrichment.md#enrichment) retain their own result structures. Follow each function's dependencies and coordinates; see the [query guide](../api/features/26_gene_annotation.md#dnakit.databases.ensembl_lookup_gene) for database calling conventions.

These notes describe functions backed by optional scientific tools and exposed
under capability-oriented namespaces: `annotation`, `comparative`, and `molbio`. Users do not
call an `integrations.tooluniverse` namespace.

Install the optional backend with:

```bash
python -m pip install "dnakit[external-tools]"
```

No function trains or fine-tunes a model. Ensembl, ClinVar, dbSNP, and gnomAD calls require network
access; dN/dS and Golden Gate run locally. One backend tool
is loaded per call. Every function returns `ProviderResult`; use `result.data` for the provider
payload; `result.parameters`, `result.metadata`, and `result.provenance` contain the call parameters, backend version, and source, respectively.


### Boundaries

- Missing optional dependencies raise `BackendUnavailableError` with a hint to install `dnakit[external-tools]`.
- Remote providers may rate-limit requests, enter maintenance, or update database versions; retain `result.provenance`.
- ClinVar, VEP, dbSNP, and gnomAD results are research annotations, not clinical diagnoses.
- These functions use bundled DNAKit adapters without importing ToolUniverse. Additional APIs are described in the [native skill interface guide](../llm_tools/skills.md#unified-skill-api).

## Prediction dependencies, checkpoints and training options {#property-prediction}

This page covers custom training and prediction, plus the parameters, returns and credentials of trained-model APIs, including AlphaGenome, ESM, Evo 2 and NVIDIA NIM. For examples, see [custom model training](../api/features/23_prediction_custom_training.md), [regulatory tracks](../api/features/23_prediction_regulatory_tracks.md), [variant effects](../api/features/23_prediction_variant_effects.md), and [protein features and structure](../api/features/23_prediction_protein_features.md). See [variant effect scores and experimental data](../api/features/36_variant_effect_scores.md) for recorded CADD and AlphaMissense scores, and [variant splicing effects](../api/features/23_prediction_variant_effects.md#variant-splicing) for SpliceAI and Pangolin.

If a prediction dependency is missing or cannot be imported, DNAKit prints installation guidance to stderr before raising the existing exception. The exception also retains the installation commands. Dependency error diagnostics remain visible with `verbose=False`; DNAKit does not install packages automatically.

From the **DNAKit project root**, install the extras for your workflow:

| Workflow | Installation command |
| --- | --- |
| Enformer, including its 27 classification tasks | `python -m pip install -e ".[neural,neural-enformer]"` |
| General PyTorch / Transformers backbones | `python -m pip install -e ".[neural]"` |
| Custom deep-learning training | `python -m pip install -e ".[training]"` |
| Basic machine learning | `python -m pip install -e ".[ml]"` |
| XGBoost / LightGBM / CatBoost | Use `.[xgboost]` / `.[lightgbm]` / `.[catboost]`, respectively |

For a published package, replace `-e ".[extras]"` with `"dnakit[extras]"`; no source directory is needed. Evo 2, Caduceus and source-based backends receive their own guidance. The Evo 2 extra requires Python 3.11–3.12; AlphaGenome, LucaOneTasks and JanusDNA also require their official environments.

Models can share an installed PyTorch or TensorFlow version when it satisfies all dependency constraints; the two frameworks are not interchangeable. These commands do not force upgrades or reinstallations of compatible packages. Missing or incompatible dependencies are resolved by pip. The current `neural` / `training` extras specify a PyTorch range of `>=2.5,<3`, rather than different exact versions for every model; third-party backends may impose tighter constraints. See [pip installation behavior](https://pip.pypa.io/en/stable/cli/pip_install/).

This page covers three paths:

1. **User-defined properties:** train classification or regression models from your own files, with full fine-tuning or linear probing.
2. **Official pretrained properties:** run inference with released task heads without retraining locally.
3. **Machine-learning properties:** train your own classification or regression model with scikit-learn or optional XGBoost, LightGBM, or CatBoost.

The existing [neural representations](../api/features/08_fingerprints.md) API remains embedding-only. Custom training is introduced first; numbered sections 1–54 describe the 54 official direct-inference functions.

Each function below describes its purpose, prediction method, API, input, and output. Placeholders in example results indicate output structure only. Actual values, label counts, and track counts depend on the input, species, tissue filters, and checkpoint.

Model papers, official repositories, and weight sources are collected in [FAQ: references for deep-learning property prediction](../faq.md#deep-learning-property-prediction-references).

### Input and output error handling

Direct prediction rejects empty sequences and invalid symbols with `record_id` in the sequence
error. Enformer benchmark inputs may be 1, 199, 200 or 201 bp; the default upper bound is
196608 bp. `max_length` may tighten but cannot relax that architecture limit. Inputs are never
silently truncated, and IUPAC symbols follow `ambiguity_policy`.

`PredictionOutput` rejects empty, non-finite, complex or unconvertible values and duplicate
output names. Direct classification/segmentation outputs must be probabilities in `[0,1]`.
Enformer classification additionally requires the class dimension to match and probabilities
to sum to one within `1e-5`. Regression, scores and tracks retain their own numeric semantics.
Custom backends must follow the selected task's output contract.

Direct prediction and embedding APIs report memory exhaustion, including wrapped backend OOM
errors, as `MODEL_OUT_OF_MEMORY`, with batch/device guidance instead of substitute predictions.
An unavailable explicit CUDA device or device index is rejected before model loading.

### Checkpoint reuse and automatic downloads {#property-prediction-checkpoint-cache}

Pretrained weights are checked locally before prediction: missing weights are downloaded only for the selected model/task, with start/completion messages and absolute local paths. Existing weights are reused. Set `verbose=False` to silence feedback. Enformer task weights default to `./ckpt/enformer-benchmarks/`; other backbones use `./ckpt/<model>/`. ML embedding backbones follow the same policy. User-trained estimators (`estimator.joblib`) and neural heads (`head.npz`) must be retained in `model_dir`; they have no registered public download source. LucaOneTasks uses its official downloader in the supplied source checkout.

For the 27 Enformer classification tasks, no checkpoint path is required. Use `checkpoint_dir` to select a cache directory: existing files are also checked in its `enformer-benchmarks/` and task-family (`nt-revised/` or `genomic-benchmarks/`) subdirectories; a missing file is downloaded directly into the selected directory. Alternatively, `checkpoint_path` selects an exact file, downloading the selected Enformer task there if absent. The two parameters are mutually exclusive.

Install with `python -m pip install "dnakit[neural,neural-enformer]"`. Google Drive downloads use `gdown`, included in these extras. Downloads are written to temporary files and checked before installation; failed or interrupted files are not reused as checkpoints. Prediction validates the task and model parameters. Read class meanings from `output.output_names`.

### Machine-learning classification and regression

Use the existing `CustomPropertyTrainingConfig`, `train_custom_property_model()` and `predict_custom_properties()` APIs with `finetune_method="machine_learning"`. Here `model` selects the estimator family; `embedding_model` separately selects an optional frozen DNA backbone.

Install ML and tabular readers without PyTorch:

```bash
python -m pip install -e ".[ml]"
```

Install the three optional boosting libraries:

```bash
python -m pip install -e ".[xgboost,lightgbm,catboost]"
```

Embedding inputs additionally require `.[neural]` and the selected backbone's extras. Missing optional libraries affect only their own algorithms and produce an installation hint.

#### Available algorithms

`available_machine_learning_models()` lists all 37 families, supported tasks and dependencies without importing optional libraries.

| Tasks | `model` names |
| --- | --- |
| Classification and regression | `ridge`, `sgd`, `svm`, `linear_svm`, `nu_svm`, `knn`, `decision_tree`, `extra_tree`, `random_forest`, `extra_trees`, `gradient_boosting`, `hist_gradient_boosting`, `adaboost`, `bagging`, `gaussian_process`, `mlp`, `dummy`, `xgboost`, `lightgbm`, `catboost` |
| Classification only | `logistic_regression`, `perceptron`, `gaussian_nb`, `multinomial_nb`, `bernoulli_nb`, `complement_nb`, `lda`, `qda` |
| Regression only | `linear_regression`, `lasso`, `elastic_net`, `bayesian_ridge`, `ard`, `huber`, `ransac`, `theil_sen`, `kernel_ridge` |

`mlp` is scikit-learn's shallow neural network; `dummy` is a baseline. This interface covers supervised prediction, not clustering, anomaly detection or automatic model search. Individual estimators retain their sample-size, parameter and computational constraints; Gaussian processes, for example, are costly on large datasets.

#### Feature sources

| `feature_source` | Behavior |
| --- | --- |
| `sequence` (default) | Sequence length, GC fraction, ambiguous-base fraction and fixed-order k-mer frequencies. `kmer_k=3`, allowed range 1–6; dimension is `3 + 4**kmer_k`. Ambiguous windows are skipped; GC uses only A/C/G/T in its denominator. Undefined frequencies become zero. No deep-learning checkpoint is required. |
| `embedding` | Frozen vectors from `embedding_model` (default `dnabert2`). Checkpoint, pooling, device and batching settings follow the representation API. |
| `numeric` | Explicit, ordered `DataMapping.feature_fields`. Set `sequence_field=None` for files without sequences. Missing, nonnumeric or infinite features fail validation; there is no automatic imputation. |

Numeric columns must be explicitly selected and must not overlap labels, ID, split or sequence fields. Prediction requires the same feature names and order. Avoid target-derived features or other leakage.

Train a random forest from sequence data:

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

config = CustomPropertyTrainingConfig(
    task_type="regression",  # or classification
    finetune_method="machine_learning",
    model="random_forest",
    model_params={"n_estimators": 300, "max_depth": 12},
    feature_source="sequence",
    kmer_k=3,
    mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
    output_dir="runs/random-forest",
    workers=6,
    seed=42,
)
result = train_custom_property_model("my_labels.csv", config=config)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir,
    trust_model=True,  # only for artifacts from a trusted source
)
print(predictions.predictions)
```

For numeric input use `feature_source="numeric"` and, for example, `DataMapping(sequence_field=None, feature_fields=("gc", "length"), label_fields=("target",))`. For embeddings use `feature_source="embedding", embedding_model="dnabert2"` and the checkpoint's required loading options, including `allow_remote_code` where applicable.

#### Training and persistence

Defaults remain seed 42 and an 80%/10%/10% split; a mapped split column takes precedence. Target semantics and metrics match the deep-learning path. Multilabel, multioutput classification and multi-target regression fit an independent estimator per target. Constant binary targets use constant classifiers. Class predictions use the highest probability; multilabel predictions use `threshold`.

`model_params` supplies serializable estimator constructor parameters, not arbitrary estimator objects. `workers` controls default `n_jobs` and numerical-library threads (0 means 1); CatBoost uses `thread_count`. Explicit estimator parameters can override their parallelism and seeds. `standardize_features=True` fits scaling only on training rows. Nonnegative-only `multinomial_nb` and `complement_nb` use MaxAbsScaler; negative inputs are rejected instead of silently shifted.

Classifiers without `predict_proba` use sigmoid calibration with up to three stratified folds inside the training split. Each class needs at least two training samples. Scaling is also fitted inside each calibration fold. Validation/test rows are excluded, and calibration is recorded in result warnings. See the [scikit-learn calibration API](https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html).

The ML path does not use neural `epochs`, AdamW settings or gradient accumulation. Set iterations/tree counts through `model_params`. Validation is evaluation-only; there is no automatic tuning or best-model selection.

Regression targets use training-set mean/std normalization, as in the existing API; reports and predictions return original units. Classification `class_weight` / `class_weights` dictionaries use zero-based encoded class indices, ordered as saved in `model.json` under `target_schema.class_names`. Options such as `"balanced"` remain available where the estimator supports them.

Artifacts are `model.json`, `estimator.joblib`, `metrics.json`, `data_inspection.json`, `splits.csv` and `test_predictions.jsonl`. Fitted transforms are stored with estimators; the manifest records feature order, target schema, settings and Python/library versions. The returned `history` is empty because there is no neural training history.

Loading joblib can execute Python code, so ML inference requires explicit `trust_model=True` for artifacts you trained or otherwise trust. Reloading checks the recorded Python and dependency versions and requires the original environment if they differ. See [scikit-learn persistence guidance](https://scikit-learn.org/stable/model_persistence.html).

Edit `others/examples/train_machine_learning_properties.py` to use your own data, or use `--demo` for synthetic training and reload verification:

```bash
python others/examples/train_machine_learning_properties.py --demo
```

### User-defined deep-learning classification and regression

Install the training and common tabular-format dependencies first: `python -m pip install "dnakit[training]"`.

`train_custom_property_model()` accepts one local data file and maps an ID, a DNA/RNA sequence, one or more targets, and an optional existing split. The default train/validation/test ratios are 0.8/0.1/0.1, with seed 42. A mapped existing split column takes precedence. The caller must explicitly set `task_type` to `classification` or `regression`; DNAKit does not infer the scientific task merely from how values look.

Built-in readers cover CSV, TSV, other single-character delimiters, JSON, JSONL, YAML, XML, Excel, Parquet, SQLite, FASTA, FASTQ, and GenBank. Dot paths address nested records, and `key=value` FASTA/FASTQ description tokens become `header.key`. A caller can implement the `DataReader` protocol for any other format. An LLM may propose only a validated declarative `DataMapping`; DNAKit never executes model-generated Python, SQL, regular expressions, or shell commands.

#### Label semantics

- One categorical label column: binary or multiclass classification.
- Multiple binary label columns: automatically detected as one-hot multiclass or multilabel classification; `classification_mode` can select the interpretation explicitly.
- Multiple nonbinary categorical columns: multioutput classification, with an independent class table for each column.
- One or more numeric label columns: single-target or multi-target regression.

Classification output includes accuracy, balanced accuracy, precision, recall, F1, Jaccard, MCC, Cohen's kappa, log loss, AUROC, average precision, confusion matrices, and per-class/per-label results. Multilabel output also includes subset accuracy, Hamming loss, LRAP, and coverage error. Regression output includes MAE, MSE, RMSE, R², explained variance, median absolute error, max error, bias, MAPE, SMAPE, Pearson, and Spearman metrics, both aggregate and per target. Undefined metrics are `null`; classification metrics such as F1 are not reported for regression.

#### Models and fine-tuning methods

The default is `finetune_method="full"`. All 11 models provide full backbone fine-tuning and frozen `linear_probe`: DNABERT-2, NT-v2, HyenaDNA, Caduceus, GROVER, LucaOne, GENERator, Enformer, Evo 2, JanusDNA, and AlphaGenome. Full training updates the backbone and new task head; probing updates only the linear task head. Query `available_custom_training_models()` for the support matrix. A registered capability does not establish that optional dependencies are installed or that every official large checkpoint has been trained locally.

Checkpoint-provided code still requires `allow_remote_code=True`. Enformer requires `enformer-pytorch`; Evo 2 requires the official Evo2/Vortex CUDA environment; JanusDNA requires its official source environment; AlphaGenome uses local research code, JAX/Haiku/Orbax, and PyTorch with a differentiable JAX bridge, not a remote inference API. Set `model_source_path` (`--model-source-path` in the CLI) for a local source checkout. Prediction accepts the same option to relocate the saved source path. Original `.bin` checkpoints require PyTorch 2.6 or newer under Transformers' safe-loading rules.

Older AlphaGenome source exposes only inference normalization statistics. Full fine-tuning keeps those running statistics fixed while updating learnable scale/offset and other backbone weights; newer source with training-mode support also updates running statistics. The backbone manifest records the actual policy.

Enformer full training defaults to 4096 bp chunks with uncropped embeddings, variable-length positional encoding, valid-position pooling, and averaging across chunks. Set `max_length` to change the chunk size. Frozen Enformer probing explicitly configured with a shorter context uses the same short-context processing. Evo 2 full training uses the final normalized representation so all backbone blocks receive gradients; frozen probing retains the registered intermediate layer. Record these choices in comparisons. Evo 2's official constructor uses all visible GPUs: set `CUDA_VISIBLE_DEVICES` to expose only the intended training GPU and select `device="cuda:0"`. The trainer remains single-device.

Full models save the backbone and restore the best validation epoch before evaluation. Reloading requires the model runtime but not the original pretrained weights. Probing saves the task head and original backbone location, so inference needs the same pretrained checkpoint. Training and prediction show progress; `show_progress=False` disables it. Validation includes reduced official architectures, gradient updates, frozen weights, reloads, and end-to-end tests; it does not claim training runs for all 11 official large checkpoints.

AlphaGenome additionally passed a single-GPU check using the local official pretrained weights: a 2048 bp input completed one forward pass, backward pass, and AdamW parameter update. This does not establish complete training, multi-epoch convergence, or full-context validation.

#### Training parameters (Python API)

Configure training directly with `CustomPropertyTrainingConfig`; no command-line arguments are needed. Full fine-tuning and linear probing share these defaults:

Use `dtype="auto"` or `"float32"` for full training; compatible CUDA backends also accept `"bfloat16"`. Pure `"float16"` optimizer parameters are rejected because AdamW can produce non-finite weights. AlphaGenome accepts `"auto"` or `"float32"`.

| Python parameter | Default | Meaning |
| --- | --- | --- |
| `workers` | `6` | CPU data-loading processes; `0` loads in the main process |
| `batch_size` | `32` | Samples per training device per batch, before gradient accumulation |
| `learning_rate` | `1e-4` | AdamW learning rate (lr) |
| `beta1` | `0.9` | AdamW first-moment decay |
| `beta2` | `0.95` | AdamW second-moment decay |
| `weight_decay` | `0.1` | AdamW weight decay |

`workers` fetch records and assemble CPU batches for training, validation, and testing; encoder/GPU operations remain in the main process. This setting does not parallelize initial file parsing or LLM requests. Changing the worker count does not change the seeded split or epoch sample order. The trainer currently runs on one selected `device` (for example `"cuda:0"`), not multi-GPU DDP; available GPU count does not multiply the batch size. `gradient_accumulation_steps` defaults to `1` and is also configurable through Python.

Scripts using `workers>0` need the `if __name__ == "__main__":` guard below. Use `workers=0` in notebooks or environments without subprocess support. Worker processes are released on completion, early stopping, or errors.

```python
import os

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    LLMSchemaConfig,
    train_custom_property_model,
)

def main():
    result = train_custom_property_model(
        "my_labels.xlsx",
        config=CustomPropertyTrainingConfig(
            task_type="regression",  # or "classification"
            model="dnabert2",
            finetune_method="full",  # or "linear_probe"
            output_dir="runs/my-property",
            workers=6,
            batch_size=32,
            learning_rate=1e-4,
            beta1=0.9,
            beta2=0.95,
            weight_decay=0.1,
            seed=42,
            train_ratio=0.8,
            validation_ratio=0.1,
            test_ratio=0.1,
            allow_remote_code=True,
            llm=LLMSchemaConfig(
                provider="deepseek",
                model="deepseek-v4-flash",
                api_key=os.environ["DEEPSEEK_API_KEY"],
            ),  # Omit llm when schema assistance is unnecessary.
        ),
    )
    print(result.metrics["test"])


if __name__ == "__main__":
    main()
```

The output directory contains `model.json`, `head.npz`, `metrics.json`, `history.json`, `data_inspection.json`, `splits.csv`, and `test_predictions.jsonl`; full fine-tuning also writes `backbone/`. The `training` section of `model.json` records the worker count, batch size, learning rate, both betas, and weight decay. `predict_custom_properties()` reloads those artifacts for later inference. The editable script `others/examples/train_custom_properties.py` can be run from an IDE after configuring the source file; LLM assistance is off by default in that script.

#### LLM file analysis and privacy boundary

`LLMSchemaConfig(provider=..., model=...)` supports the following connection presets through OpenAI-compatible `POST /chat/completions`, without vendor SDK dependencies. Query `available_schema_llm_providers()` offline to list them. Except for the legacy DeepSeek default, explicitly select a model ID available in your account; DNAKit does not automatically substitute a model or provider.

| `provider` | Default `base_url` | API-key environment variable |
| --- | --- | --- |
| `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| `custom` | Required explicitly | `LLM_API_KEY` |

For example, use Google Gemini to analyze only field names and type structure:

```python
from dnakit.predictions import LLMSchemaConfig, inspect_property_data

llm = LLMSchemaConfig(
    provider="gemini",
    model="gemini-3.8-flash",  # Replace with an available text/JSON model if needed.
    include_sample_values=False,
)
inspection = inspect_property_data("my_labels.xlsx", llm=llm)
print(inspection.mapping.to_dict())
print(inspection.recommendations)
## Pass the same llm into CustomPropertyTrainingConfig when training.
```

Override `api_key_env` for your own environment-variable name or pass `api_key` in Python without committing it to source control. Explicit `base_url` overrides support regional or private gateways. HTTP redirects are refused so that keys and sampled profiles are not forwarded to another endpoint. Gemini uses its [official compatibility endpoint](https://ai.google.dev/gemini-api/docs/openai), not native `generateContent` or Vertex AI. Qwen defaults to the still-supported Beijing endpoint; a workspace-specific endpoint such as `https://<workspaceID>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` can also be used. Endpoints and API keys must match the region, as described in [Alibaba Cloud's documentation](https://help.aliyun.com/zh/model-studio/compatibility-of-openai-with-dashscope). For Kimi's international endpoint, set `base_url="https://api.moonshot.ai/v1"`, per its [quickstart](https://platform.kimi.ai/docs/overview). OpenRouter requires its complete model ID; check the chosen model's supported parameters and downstream provider/data policies in the [quickstart](https://openrouter.ai/docs/quickstart) and [structured-output documentation](https://openrouter.ai/docs/guides/features/structured-outputs).

Optional Python settings:

- `temperature=None` (default) omits temperature for models that reject custom sampling parameters.
- `token_limit_parameter="auto"` sends `max_completion_tokens` for OpenAI and `max_tokens` for other compatible services; either field can also be selected explicitly. This follows [OpenAI Docs](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create), which notes that the older `max_tokens` field is incompatible with some reasoning models.
- `json_mode=True` requests JSON Object output. Set it to `False` for a model without this parameter; the prompt and mapping validation still require a JSON object.
- `extra_body` passes model-specific fields, such as `{"enable_thinking": False}` for supported Qwen models or `{"thinking": {"type": "disabled"}}` for supported Kimi models. Check the [Qwen](https://help.aliyun.com/zh/model-studio/qwen-structured-output) or [Kimi parameter reference](https://platform.kimi.ai/docs/api/models-overview) before applying these options; DNAKit does not force thinking off for all models. Requests are non-streaming, so models requiring streamed thinking must disable thinking or be replaced with a compatible model.

`extra_body` cannot replace authentication, messages, model selection, token limits, or other core request fields. Truncated responses fail explicitly: adjust `max_output_tokens` or thinking settings yourself; DNAKit never increases the budget or retries another provider automatically. Presets describe interface support, not universal model compatibility. Provider contracts are covered by mocked tests, not paid live tests against every service.

`LLMSchemaConfig` uses OpenAI-compatible `POST /chat/completions` with JSON Output, including DeepSeek and compatible providers. By default, it shares at most 8 sampled records, 12,000 characters, and a structural profile—not the whole source. Even a file with fewer than 8 records retains at least one complete unshared record. Set `include_sample_values=False` to share only field/type structure. `DataInspection` records how many sample rows and characters were sent. API keys are excluded from request bodies, results, exception displays, and training artifacts; the Python example reads the key from an environment variable.

## Enformer task checkpoint paths {#task-checkpoints}

- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H2AFZ.ckpt`; the default path is `./ckpt/enformer-benchmarks/H2AFZ.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K27ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27ac.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K27me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27me3.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K36me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K36me3.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K4me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me1.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K4me2.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me2.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K4me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me3.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K9ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9ac.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H3K9me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9me3.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `H4K20me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H4K20me1.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `enhancers.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `enhancers_types.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers_types.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `promoter_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_all.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `promoter_no_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `promoter_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_tata.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `splice_sites_acceptors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `splice_sites_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `splice_sites_donors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `demo_coding_vs_intergenomic_seqs.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `demo_human_or_worm.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `drosophila_enhancers_stark.ckpt`; the default path is `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `dummy_mouse_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `human_enhancers_cohn.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `human_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `human_ensembl_regulatory.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `human_nontata_promoters.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`; a custom directory or exact path can also be used.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#property-prediction-checkpoint-cache)) `human_ocr_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`; a custom directory or exact path can also be used.
## Script execution, report rendering and validation {#engineering}

### Run a bundled coordinate conversion {#engineering-bundled-script-example}

This runs a packaged Python script and converts BED 0-based half-open coordinates to GFF 1-based inclusive coordinates, not between genome assemblies.

```python
import dnakit

print(dnakit.workflows.scripts("genomic-coordinates", verbose=False))
print(dnakit.workflows.resource_path(
    "genomic-coordinates", "scripts/convert_coords.py", verbose=False,
))
result = dnakit.workflows.run(
    "genomic-coordinates", "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60, verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

`run_support()` accepts `clinical-variant-reporter`, `clinical-trial-finder`, `omics-target-evidence-mapper`, `gwas-lookup` and `fine-mapping`. Use `arguments=["--help"]` to inspect their CLI arguments and prepare the scientific dependencies imported by each script. See [workflow parameters and resource verification](../api/features/14_engineering.md#execution).

### Render an existing report {#engineering-report-example}

```bash
python -m pip install -e ".[reports]"
```

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import dnakit

with TemporaryDirectory() as directory:
    root = Path(directory)
    report = root / "report.md"
    report.write_text("# Results\nUser-supplied report content.\n", encoding="utf-8")
    result = dnakit.workflows.write_wes_report(report, root / "pdf", verbose=False)
    print(result["output"])
    print(result["manifest"])
    print(Path(result["output"]).exists())  # True
```

Use your own output directory to keep the files. This renderer preserves supplied content without generating medical interpretation. `language` selects English or Spanish labels without translating the text. See [title, logos, font and overwrite options](../api/features/14_engineering.md#dnakit.workflows.write_wes_report).

### MCP and validation scope {#engineering-mcp-validation}

The same functions are available through tool discovery, description and invocation in the [MCP guide](../llm_tools/mcp.md). Python, CLI and MCP exposures do not count as separate capabilities.

Documentation checks verify 484 unique function entries, current signatures and parameters, and build both languages. Documentation coverage and offline example checks do not establish end-to-end validation of all remote services, models or scientific workflows.

## ENCODE query and download previews {#encode-download}

### ENCODE queries and download previews {#encode-download-encode-download-example}

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)
preview = json.loads(dnakit.encode.batch_download(
    download_dir="downloads/encode", assay_title="ATAC-seq",
    assembly="GRCh38", file_format="bed", limit=3, dry_run=True,
))
print(preview)
```

`dry_run=True` queries the file list without downloading the files. Set it to `False` to download. ENCODE functions return text: mostly JSON, but citation and table exports also support BibTeX, RIS, CSV and TSV. Do not parse those formats with `json.loads()`. Tracking, credentials and downloads can write local state. See [experiment queries](../api/features/22_database_query_encode.md#encode), [collection management](../api/features/10_datasets.md#encode-collection) and [credential management](../api/features/14_engineering.md#encode-credentials) for the related APIs.
