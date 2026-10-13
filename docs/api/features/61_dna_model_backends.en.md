# Data utilities for optional DNA model packages

This group of functions is a set of optional-dependency adapters in `dnakit.predictions`. They cover DNA-level data operations of [alphagenome](https://github.com/google-deepmind/alphagenome), [grelu](https://github.com/Genentech/gReLU), [selene-sdk](https://github.com/FunctionLab/selene), [enformer-pytorch](https://github.com/lucidrains/enformer-pytorch), [genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks), [evo2](https://github.com/ArcInstitute/evo2), and [kipoi](https://github.com/kipoi/kipoi). They include interval and variant parsing, GENCODE gene and transcript extraction, sequence retrieval and transformation, one-hot encoding, BED windows, benchmark dataset summaries, local Evo 2 generation, and kipoi model descriptions. These functions do not run model inference or download weights. Model inference and hosted prediction APIs are already provided by `dnakit.predictions` and `dnakit.databases`, so they are not repeated here.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## 1) `dnakit.predictions.parse_alphagenome_interval()`

- **Purpose:** Parse an interval string of the form `chr:start-end[:strand]`.
- **Method:** Call alphagenome's `Interval.from_str`. The start is used as given and treated as 0-based. The end is exclusive. No 1-based conversion is applied.
- **API:** `dnakit.predictions.parse_alphagenome_interval(text)`.
- **Input:** An interval string, for example `"chr1:101-200:+"`.
- **Output:** An `interval` field with `chromosome`, `start`, `end`, `strand`, and `width`. The coordinate convention is recorded in `provenance`.
- **Example code:**

```python
from dnakit.predictions import parse_alphagenome_interval

result = parse_alphagenome_interval("chr1:101-200:+")
print(result["interval"])
```

- **Example output:**

```text
{'chromosome': 'chr1', 'start': 101, 'end': 200, 'strand': '+', 'width': 99}
```

## 2) `dnakit.predictions.parse_alphagenome_variant()`

- **Purpose:** Parse a variant string of the form `chr:position:REF>ALT` and report its variant type.
- **Method:** Call alphagenome's `Variant.from_str`. `position` is read as a 1-based VCF position, so `start` is `position - 1`. `end` is the exclusive end computed by alphagenome for the variant type. For example, `chr1:10:AT>A` has end 11. `variant_format` can be `default`, `gtex`, `open_targets`, `open_targets_bigquery`, or `gnomad`.
- **API:** `dnakit.predictions.parse_alphagenome_variant(text, *, variant_format="default")`.
- **Input:** A variant string and an optional `variant_format`.
- **Output:** A `variant` field with `chromosome`, `position`, `reference`, `alternate`, `start`, `end`, and the flags `is_snv`, `is_indel`, `is_insertion`, `is_deletion`, and `is_structural`.
- **Example code:**

```python
from dnakit.predictions import parse_alphagenome_variant

variant = parse_alphagenome_variant("chr1:10:A>G")["variant"]
print(variant["start"], variant["end"], variant["is_snv"])
```

- **Example output:**

```text
9 10 True
```

## 3) `dnakit.predictions.merge_alphagenome_intervals()`

- **Purpose:** Merge overlapping intervals on the same chromosome.
- **Method:** Parse each string with `Interval.from_str`, then call alphagenome's `merge_overlapping_intervals`.
- **API:** `dnakit.predictions.merge_alphagenome_intervals(intervals)`.
- **Input:** A list of interval strings, up to 10,000 entries.
- **Output:** `input_count`, `merged_count`, and the merged `intervals` list.
- **Example code:**

```python
from dnakit.predictions import merge_alphagenome_intervals

result = merge_alphagenome_intervals(["chr1:0-10", "chr1:5-20", "chr1:30-40"])
print([(item["start"], item["end"]) for item in result["intervals"]])
```

- **Example output:**

```text
[(0, 20), (30, 40)]
```

## 4) `dnakit.predictions.read_gencode_gene_intervals()`

- **Purpose:** Return the intervals of selected genes from a GENCODE GTF.
- **Method:** This module parses the `gene` rows of the GTF in pure Python. It converts the 1-based closed GTF start to 0-based by subtracting 1 and keeps the end as the exclusive end. It then passes the table to alphagenome's `get_gene_intervals`. Exactly one of `gene_symbols` and `gene_ids` must be given.
- **API:** `dnakit.predictions.read_gencode_gene_intervals(gtf_path, *, gene_symbols=None, gene_ids=None, max_records=200000)`.
- **Input:** A GTF path, which can be `.gz`; a list of gene symbols or Ensembl gene IDs; and `max_records`, the limit on selected rows.
- **Output:** A `genes` list. Each item has the interval coordinates, `strand`, and `name`.
- **Example file:** `genes.gtf`

```text
chr1	test	gene	101	500	.	+	.	gene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";
chr1	test	gene	700	800	.	-	.	gene_id "G2"; gene_name "GENE2"; gene_type "lincRNA";
```

- **Example code:**

```python
from dnakit.predictions import read_gencode_gene_intervals

result = read_gencode_gene_intervals("genes.gtf", gene_symbols=["GENE1"])
print(result["genes"])
```

- **Example output:**

```text
[{'chromosome': 'chr1', 'start': 100, 'end': 500, 'strand': '+', 'width': 400, 'name': 'GENE1'}]
```

## 5) `dnakit.predictions.extract_gencode_transcripts()`

- **Purpose:** Return transcripts that overlap an interval, with their exons and CDS.
- **Method:** This module parses the `transcript`, `exon`, `CDS`, `start_codon`, and `stop_codon` rows of the GTF in pure Python, using the same coordinate conversion as the previous function. It then calls alphagenome's `TranscriptExtractor`. The CDS intervals come from alphagenome's `Transcript` object and can differ from the raw GTF rows.
- **API:** `dnakit.predictions.extract_gencode_transcripts(gtf_path, region, *, max_records=200000)`.
- **Input:** A GTF path and a 0-based half-open interval string, for example `"chr1:0-1000"`.
- **Output:** A `transcripts` list. Each item has `transcript_id`, `gene_id`, `strand`, `exons`, `cds`, and `is_coding`. All coordinates are 0-based half-open.
- **Example file:** `tx.gtf`

```text
chr1	test	transcript	101	500	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	exon	101	201	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	exon	301	500	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	CDS	151	198	.	+	0	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	stop_codon	199	201	.	+	0	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
```

- **Example code:**

```python
from dnakit.predictions import extract_gencode_transcripts

result = extract_gencode_transcripts("tx.gtf", "chr1:0-1000")
for item in result["transcripts"]:
    print(item["transcript_id"], item["exons"], item["cds"])
```

- **Example output:**

```text
T1 [[100, 201], [300, 500]] [[150, 198]]
```

## 6) `dnakit.predictions.reverse_complement_with_grelu()`

- **Purpose:** Reverse-complement one or more DNA sequences.
- **Method:** Call grelu's `reverse_complement`. Only `A`, `C`, `G`, `T`, and `N` are accepted. Input is converted to upper case.
- **API:** `dnakit.predictions.reverse_complement_with_grelu(sequences)`.
- **Input:** A DNA string or a list of DNA strings, up to 10,000 entries.
- **Output:** A `sequences` list. A single string input also returns a list.
- **Example code:**

```python
from dnakit.predictions import reverse_complement_with_grelu

print(reverse_complement_with_grelu(["AAC", "GGTN"])["sequences"])
```

- **Example output:**

```text
['GTT', 'NACC']
```

## 7) `dnakit.predictions.shuffle_dinucleotides()`

- **Purpose:** Generate random shuffles of a sequence that keep its dinucleotide counts.
- **Method:** Call grelu's `dinuc_shuffle`, which runs on tangermeme. If several shuffles are requested and all of them come out identical, upstream raises an error. This function converts that case to `ConfigurationError`.
- **API:** `dnakit.predictions.shuffle_dinucleotides(sequence, *, n_shuffles=1, seed=None)`.
- **Input:** A DNA sequence. `n_shuffles` is an integer from 1 to 1000. `seed` is optional and must be from 0 to 2^32 - 1.
- **Output:** A `shuffles` list. Each shuffle has the same length as the input.
- **Example code:**

```python
from dnakit.predictions import shuffle_dinucleotides

result = shuffle_dinucleotides("AACCGGTTAGCTAGGATCCAGTAC", n_shuffles=2, seed=7)
print(len(result["shuffles"]), result["shuffles"][0][:8])
```

- **Example output:**

```text
2 AGGTAGGA
```

## 8) `dnakit.predictions.resize_dna_sequences()`

- **Purpose:** Trim or pad DNA sequences with `N` to a target length.
- **Method:** Call grelu's `resize`. `end` can be `left`, `right`, or `both`. Padding uses `N`.
- **API:** `dnakit.predictions.resize_dna_sequences(sequences, seq_len, *, end="both")`.
- **Input:** A DNA string or a list of DNA strings, the target length `seq_len`, and `end`, which selects the side to change.
- **Output:** A `sequences` list. `seq_len` and `end` are also returned.
- **Example code:**

```python
from dnakit.predictions import resize_dna_sequences

print(resize_dna_sequences(["ACGT", "ACGTAC"], 6)["sequences"])
```

- **Example output:**

```text
['NACGTN', 'ACGTAC']
```

## 9) `dnakit.predictions.mutate_dna_sequence()`

- **Purpose:** Replace the base at one position with another base.
- **Method:** Call grelu's `mutate`. `position` is a 0-based index. The original base is returned so the change can be checked.
- **API:** `dnakit.predictions.mutate_dna_sequence(sequence, allele, *, position)`.
- **Input:** A DNA sequence, a single-base `allele`, and a 0-based `position` inside the sequence.
- **Output:** `sequence` (the mutated sequence), `reference_base`, `position`, and `allele`.
- **Example code:**

```python
from dnakit.predictions import mutate_dna_sequence

result = mutate_dna_sequence("ACGT", "T", position=1)
print(result["sequence"], result["reference_base"])
```

- **Example output:**

```text
ATGT C
```

## 10) `dnakit.predictions.fetch_genome_sequences()`

- **Purpose:** Retrieve sequences for given intervals from a FASTA file.
- **Method:** Call selene-sdk's `Genome` and `get_sequence_from_coords`. Intervals are 0-based half-open. With `pad=False`, an out-of-bounds interval raises an error, because selene returns an empty string for it. With `pad=True`, the sequence is padded with `N`. pyfaidx writes a `.fai` index file next to the FASTA.
- **API:** `dnakit.predictions.fetch_genome_sequences(fasta_path, regions, *, pad=False, max_regions=1000, max_region_length=100000)`.
- **Input:** A FASTA path and a list of regions. Each region is a dictionary with `chrom`, `start`, `end`, and an optional `strand`.
- **Output:** A `regions` list. Each item has `chrom`, `start`, `end`, `strand`, and `sequence`.
- **Example file:** `genome.fa`

```text
>chr1
ACGTACGTNN
>chr2
GGGCCCAATT
```

- **Example code:**

```python
from dnakit.predictions import fetch_genome_sequences

result = fetch_genome_sequences(
    "genome.fa",
    [{"chrom": "chr1", "start": 0, "end": 4}, {"chrom": "chr1", "start": 8, "end": 12}],
    pad=True,
)
print([item["sequence"] for item in result["regions"]])
```

- **Example output:**

```text
['ACGT', 'NNNN']
```

## 11) `dnakit.predictions.read_contig_lengths()`

- **Purpose:** List the name and length of each sequence in a FASTA file.
- **Method:** Read the index with selene-sdk's `Genome` and call `get_chr_lens`.
- **API:** `dnakit.predictions.read_contig_lengths(fasta_path, *, max_contigs=100000)`.
- **Input:** A FASTA path.
- **Output:** A `contigs` list. Each item has `name` and `length`. `contig_count` is also returned.
- **Example code:**

```python
from dnakit.predictions import read_contig_lengths

print(read_contig_lengths("genome.fa")["contigs"])
```

- **Example output:**

```text
[{'name': 'chr1', 'length': 10}, {'name': 'chr2', 'length': 10}]
```

## 12) `dnakit.predictions.encode_one_hot_selene()`

- **Purpose:** One-hot encode DNA with selene-sdk's encoding.
- **Method:** Call `Genome.sequence_to_encoding`. The channel order is A, C, G, T. `N` is 0.25 in every channel, which differs from enformer-pytorch.
- **API:** `dnakit.predictions.encode_one_hot_selene(sequence, *, max_length=10000)`.
- **Input:** A DNA sequence no longer than `max_length`. The upper limit is 100,000.
- **Output:** A `one_hot` matrix with one row per base and four columns, plus `length` and `channels`.
- **Example code:**

```python
from dnakit.predictions import encode_one_hot_selene

result = encode_one_hot_selene("ACGTN")
print(result["one_hot"][4])
```

- **Example output:**

```text
[0.25, 0.25, 0.25, 0.25]
```

## 13) `dnakit.predictions.encode_one_hot_enformer()`

- **Purpose:** One-hot encode DNA with enformer-pytorch's encoding.
- **Method:** Call `enformer_pytorch.data.str_to_one_hot`. The channel order is A, C, G, T. `N` is all zeros in every channel, which differs from selene-sdk.
- **API:** `dnakit.predictions.encode_one_hot_enformer(sequence, *, max_length=10000)`.
- **Input:** A DNA sequence no longer than `max_length`. The upper limit is 100,000.
- **Output:** A `one_hot` matrix with one row per base and four columns, plus `length` and `channels`.
- **Example code:**

```python
from dnakit.predictions import encode_one_hot_enformer

result = encode_one_hot_enformer("ACGTN")
print(result["one_hot"][4])
```

- **Example output:**

```text
[0.0, 0.0, 0.0, 0.0]
```

## 14) `dnakit.predictions.read_enformer_bed_windows()`

- **Purpose:** Build enformer-pytorch input windows for BED intervals and count the padding and unknown bases in each window.
- **Method:** Call `enformer_pytorch.data.GenomeIntervalDataset`. Intervals shorter than `context_length` are extended symmetrically. The extension positions are padding `.` characters, which encode as 0.25 in every channel. Intervals longer than `context_length` keep their BED length. The function returns counts only and does not return the full one-hot matrices.
- **API:** `dnakit.predictions.read_enformer_bed_windows(bed_path, fasta_path, *, context_length=196608, max_windows=100)`.
- **Input:** A BED file with at least three tab-separated columns, a FASTA path, and limits for `context_length` and `max_windows`.
- **Output:** A `windows` list. Each item has `chromosome`, `start`, `end`, `length`, `padding_positions`, `unknown_positions`, and `acgt_positions`.
- **Example file:** `windows.bed`

```text
chr1	0	4
chr1	6	10
```

- **Example code:**

```python
from dnakit.predictions import read_enformer_bed_windows

result = read_enformer_bed_windows("windows.bed", "genome.fa", context_length=8)
for window in result["windows"]:
    print(window["chromosome"], window["length"], window["padding_positions"], window["unknown_positions"])
```

- **Example output:**

```text
chr1 8 2 0
chr1 8 2 2
```

## 15) `dnakit.predictions.list_genomic_benchmark_datasets()`

- **Purpose:** List genomic-benchmarks dataset names.
- **Method:** Without `dataset_path`, the function reads the dataset index that ships with the package and makes no network access. With a directory, it lists the subdirectories of that directory.
- **API:** `dnakit.predictions.list_genomic_benchmark_datasets(*, dataset_path=None)`.
- **Input:** An optional local directory.
- **Output:** `source` (`package_index` or `local_directory`), `count`, and a sorted `datasets` list.
- **Example code:**

```python
from dnakit.predictions import list_genomic_benchmark_datasets

result = list_genomic_benchmark_datasets()
print(result["count"], result["datasets"][:3])
```

- **Example output:**

```text
9 ['demo_coding_vs_intergenomic_seqs', 'demo_human_or_worm', 'drosophila_enhancers_stark']
```

## 16) `dnakit.predictions.describe_genomic_benchmark()`

- **Purpose:** Count the sequences in each class and split of a local genomic-benchmarks dataset.
- **Method:** Call `genomic_benchmarks.data_check.info` with `local_repo=True`, so there is no network access. The summary lines that upstream prints are captured in `messages` and are not written to standard output.
- **API:** `dnakit.predictions.describe_genomic_benchmark(dataset_path, *, version=None)`.
- **Input:** A dataset directory that contains `metadata.yaml`, `train/`, and `test/`. `version` is optional, and a mismatch with the metadata raises an error.
- **Output:** `classes` (a mapping from each class to its `train` and `test` counts), `total_sequences`, and `messages`.
- **Example code:**

```python
import gzip
from pathlib import Path

from dnakit.predictions import describe_genomic_benchmark

root = Path("demo_ds")
for split in ("train", "test"):
    (root / split).mkdir(parents=True, exist_ok=True)
    for label in ("neg", "pos"):
        with gzip.open(root / split / f"{label}.csv.gz", "wt") as handle:
            handle.write("id,region,start,end\n1,chr1,10,20\n2,chr1,30,45\n")
(root / "metadata.yaml").write_text("version: 0\nclasses:\n  - neg\n  - pos\n")

result = describe_genomic_benchmark(root, version=0)
print(result["classes"], result["total_sequences"])
```

- **Example output:**

```text
{'neg': {'train': 2, 'test': 2}, 'pos': {'train': 2, 'test': 2}} 8
```

## 17) `dnakit.predictions.check_genomic_benchmark_cache()`

- **Purpose:** Check whether a genomic-benchmarks dataset exists in the local cache folder.
- **Method:** Call `genomic_benchmarks.data_check.is_downloaded`, which only checks whether the folder exists. It does not download anything. Without `cache_path`, the package default `~/.genomic_benchmarks` is used.
- **API:** `dnakit.predictions.check_genomic_benchmark_cache(dataset, *, cache_path=None)`.
- **Input:** A dataset name that may contain only letters, digits, underscores, and hyphens, plus an optional cache folder.
- **Output:** `downloaded` (a boolean), `dataset`, and `cache_path`.
- **Example code:**

```python
from dnakit.predictions import check_genomic_benchmark_cache

result = check_genomic_benchmark_cache("demo_ds", cache_path="empty_cache")
print(result["downloaded"])
```

- **Example output:**

```text
False
```

## 18) `dnakit.predictions.generate_evo2_sequences()`

- **Purpose:** Continue DNA prompts with a local Evo 2 checkpoint.
- **Method:** Construct `evo2.Evo2(model_name=..., local_path=...)` and call `generate`. `weights_path` must be an existing local file, and the function never triggers the Hugging Face download. `model_name` must be one of the names that evo2 supports.
- **API:** `dnakit.predictions.generate_evo2_sequences(prompts, *, weights_path, model_name="evo2_7b", n_tokens=100, temperature=1.0, top_k=4, top_p=1.0, max_prompts=8, max_prompt_length=10000)`.
- **Input:** A list of DNA prompts, a local weights file, and sampling parameters.
- **Output:** `sequences` and `scores`. `scores` holds the score that upstream returns for each sequence.
- **Example code (not run locally):**

```python
from dnakit.predictions import generate_evo2_sequences

result = generate_evo2_sequences(
    ["ACGTACGTAC"],
    weights_path="/path/to/evo2_7b.pt",
    n_tokens=32,
)
print(result["sequences"])
```

- **Example output (not run locally):**

```text
(not run locally; no output)
```

## 19) `dnakit.predictions.list_kipoi_models()`

- **Purpose:** List models from kipoi sources with their name, version, type, and input and target names.
- **Method:** Call `kipoi.list_models`. kipoi reads git-backed model sources, and its first use can clone a repository. Without `sources`, every configured source is read.
- **API:** `dnakit.predictions.list_kipoi_models(*, sources=None, limit=1000)`.
- **Input:** An optional list of source names, and `limit`, the maximum number of rows returned.
- **Output:** `total_models`, `returned`, and a `models` list.
- **Example code (not run locally):**

```python
from dnakit.predictions import list_kipoi_models

result = list_kipoi_models(sources=["kipoi"], limit=5)
print(result["total_models"])
```

- **Example output (not run locally):**

```text
(not run locally; no output)
```

## 20) `dnakit.predictions.describe_kipoi_model()`

- **Purpose:** Return the type, version, documentation, license, citation, and input and target names of one kipoi model.
- **Method:** Call `kipoi.get_model_descr`. The input and target names are read from a dictionary, a list, or a single array object, which are the three forms kipoi accepts.
- **API:** `dnakit.predictions.describe_kipoi_model(model, *, source="kipoi")`.
- **Input:** A model name such as `"DeepSEA/variantEffects"`, and a source name that defaults to `"kipoi"`.
- **Output:** A `description` dictionary with `type`, `version`, `doc`, `license`, `cite_as`, `inputs`, and `targets`.
- **Example code (not run locally):**

```python
from dnakit.predictions import describe_kipoi_model

result = describe_kipoi_model("DeepSEA/variantEffects")
print(result["description"]["inputs"])
```

- **Example output (not run locally):**

```text
(not run locally; no output)
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-alphagenome]"          # alphagenome: intervals, variants, GENCODE helpers
python -m pip install "dnakit[bio-grelu]"                # grelu: DNA string utilities
python -m pip install "dnakit[bio-selene]"               # selene-sdk: FASTA sequences and one-hot encoding
python -m pip install "dnakit[bio-enformer]"             # enformer-pytorch: one-hot encoding and BED windows
python -m pip install "dnakit[bio-genomic-benchmarks]"   # genomic-benchmarks: dataset summaries
python -m pip install "dnakit[bio-evo2]"                 # evo2: local generation, needs a GPU and weights
python -m pip install "dnakit[bio-kipoi]"                # kipoi: model descriptions, needs network access
```

Project pages: [alphagenome](https://github.com/google-deepmind/alphagenome), [grelu](https://github.com/Genentech/gReLU), [selene-sdk](https://github.com/FunctionLab/selene), [enformer-pytorch](https://github.com/lucidrains/enformer-pytorch), [genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks), [evo2](https://github.com/ArcInstitute/evo2), [kipoi](https://github.com/kipoi/kipoi).

## How these functions are verified

The test file is `others/tests/unit/test_dna_model_backends.py`, with 43 tests. 42 pass. One is a kipoi live test that is skipped unless `DNAKIT_LIVE_TESTS=1` is set. The tests check the following:

- 17 functions call the real upstream packages: 5 alphagenome, 4 grelu, 3 selene-sdk, 2 enformer-pytorch, and 3 genomic-benchmarks. Hand-calculated results are compared with direct upstream calls. Examples include GTF coordinates converted from 1-based to 0-based, the encoding of `N`, and window padding counts.
- `generate_evo2_sequences`, `list_kipoi_models`, and `describe_kipoi_model` are tested only with fake modules injected through `sys.modules`. These tests check arguments, return structure, and error mapping. They do not show what the upstream models output. evo2 cannot be imported in this environment because `flash_attn_2_cuda` is missing. A separate test checks the install hint for this real import failure.
- For each of the 7 optional packages, a test checks that a missing package raises `BackendUnavailableError` with the matching extra name.
- Input checks reject invalid DNA characters, booleans, out-of-range positions, lengths above the limits, empty strings, and NUL characters.

## Current scope

This batch covers DNA-level data operations from 7 packages, 20 functions in total. The following are not wrapped:

- alphagenome hosted prediction and local JAX inference: already provided by `dnakit.databases` and `dnakit.predictions`.
- evo2 sequence scoring (`score_sequences`): already used internally by the Evo 2 variant-effect backend in `dnakit.predictions`. This batch wraps only local generation.
- enformer-pytorch `from_pretrained` and track prediction: already covered by `dnakit.predictions`.
- selene-sdk model training, in-silico mutagenesis, and model prediction: these need a trained model, so they are not wrapped in this batch.
- grelu model training, interpretation, and variant-effect prediction: these need a model, so they are not wrapped in this batch.
- genomic-benchmarks `loc2seq`, `seq2loc`, and dataset downloads: these need network access to download genomes, so they are not wrapped in this batch.
- kipoi model download, dependency installation, and model execution: these need network access and model weights, so they are not wrapped in this batch.
