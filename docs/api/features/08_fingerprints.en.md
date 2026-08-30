# DNACharacterization

Characterize DNA sequence embeddings.

## 1) FP-001 integer encoding

- **Function:** Convert each DNA character into an integer according to the explicit code table, while retaining the original position order, which is convenient for word segmentation, model input and reversible checking of encoding rules.
- **API:** `dnakit.fingerprints.integer_encode(value[required], ambiguity_policy[optional], gap_policy[optional], max_output_length[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional IUPAC policy, Gap policy and maximum output length.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import integer_encode

result = integer_encode(DNASequence("ACGT"))
print(result.values)
```

- **Example results:**

```text
(0, 1, 2, 3)
```

## 2) FP-002 One-hot encoding

- **Function:** Convert each base into a four-dimensional vector corresponding to A/C/G/T, forming a "sequence length × 4" position retention matrix, which is convenient for use in statistical models and neural networks.
- **API:** `dnakit.fingerprints.one_hot_encode(value[required], ambiguity_policy[optional], gap_policy[optional], base_order[optional], max_output_length[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional IUPAC, Gap, Column Order and Output Length Policy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import one_hot_encode

result = one_hot_encode(DNASequence("ACGT"))
print(result.values[0], result.values[-1])
```

- **Example results:**

```text
(1.0, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 1.0)
```

## 3) FP-003 k-mer features

- **Function:** Combine the number, frequency or presence of all k-mers into a feature vector in a fixed order for sequence comparison and traditional machine learning; reverse complementary k-mers can be merged into the same feature through `canonical=True`.
- **API:** `dnakit.fingerprints.kmer(value[required], k[required], canonical[optional], mode[optional], representation[optional], ambiguity_policy[optional], overlapping[optional], cross_gaps[optional], max_dimension[optional])` is recommended; the old name `kmer_fingerprint(...)` is retained as a compatible alias only.
- **Input:** Required sequence and `k`; optional `canonical`, `mode="count"|"frequency"|"binary"`, storage representation, overlapping policy, `ambiguity_policy="error"|"ignore"` and `cross_gaps`.
- **Term:** `count`, `frequency`, and `binary` are classified as k-mer features; `binary` is an existential feature. When a fixed length DNA fingerprint is required, use [Hashed k-mer bit fingerprint](08_feature_engineering.md#1-hashed-k-mer).
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import kmer

sequence = DNASequence("ACGT")
standard = kmer(sequence, k=2, representation="sparse")
canonical = kmer(sequence, k=2, canonical=True, representation="sparse")
print(dict(standard.values), standard.dimension)
print(dict(canonical.values), canonical.dimension)
```

- **Example results:**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 16
{'AC': 2, 'CG': 1} 10
```

## 4) FP-005 k-mer Sketch (MinHash/FracMinHash)

- **Function:** Perform Bottom-k or threshold sampling on the k-mer hash of the sequence to generate a smaller and reproducible Sketch for approximating Jaccard similarity and retrieving large batches of sequences.
- **API:** `dnakit.fingerprints.minhash(value[required], k[required], num_hashes[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, `dnakit.fingerprints.fracminhash(value[required], k[required], scaled[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`
- **Input:** Required sequence and `k`; optional `num_hashes` for MinHash, optional `scaled` for FracMinHash, optional seed/canonical/upper limit for both.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash

sequence = DNASequence("ACGTACGT")
first = minhash(sequence, k=3, num_hashes=3, seed=7)
second = minhash(sequence, k=3, num_hashes=3, seed=7)
print(first.hashes == second.hashes, len(first.hashes))
```

- **Example results:**

```text
True 2
```

## 5) Neural network representation

- **Function:** Use the pre-trained DNA basic model to compress each sequence into a fixed-length floating point vector for clustering, retrieval or downstream machine learning; the meaning of the vector depends on the model, checkpoint and pooling configuration.

The neural network representation is the rep extraction part of `DATA-027`: converting each DNA sequence into a fixed-length float32 vector using a pretrained DNA base model. See [DATA-027 Neural Network Clustering ](10_clustering.md#data-027-neural-clustering) for the process of using rep with k-means.

### 5.1 Representation API

- **API:** `dnakit.representations.extract_representations(records[required], config[optional])`.
- **Configuration:** `dnakit.representations.RepresentationConfig`, optional model, checkpoint, pooling, device, dtype, batch size, maximum length and progress bar.
- **Input:** One or more `DNARecord`; does not accept empty sequence or explicit `Gap`.
- **Output:** `RepresentationResult`, where `representations` is a read-only two-dimensional matrix, with each row corresponding to one input record; the actual dimensions are stored in `embedding_dimension`.
- **Aggregation rules:** Support `mean`, `cls`, `max`, `last` pooling. Sequences that exceed the model context window are divided into blocks, and the pooled vectors of each block are finally equally weighted and averaged.

### 5.2 Rep of 11 types of neural networks supported

The following `model` names come directly from the current `MODEL_REGISTRY`:

| Serial number | `model` | rep source | Default official checkpoint | Loading conditions |
| ---: | ------------------- | -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| 1 | `alphagenome` | AlphaGenome all folds | [google/alphagenome-all-folds](https://huggingface.co/google/alphagenome-all-folds) | Requires acceptance of model terms and installation of [official research code](https://github.com/google-deepmind/alphagenome_research) |
| 2 | `caduceus` | Caduceus-Ph 131k | [kuleshov-group/caduceus-ph 131k](https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16) | `neural`, `neural-caduceus` extras; must be explicit `allow_remote_code=True` |
| 3 | `dnabert2` | DNABERT-2-117M | [zhihan1996/DNABERT-2-117M](https://huggingface.co/zhihan1996/DNABERT-2-117M) | `neural` extra; must be explicit `allow_remote_code=True` |
| 4 | `enformer` | Enformer PyTorch | [EleutherAI/enformer-official-rough](https://huggingface.co/EleutherAI/enformer-official-rough) | `neural`, `neural-enformer` extras |
| 5 | `evo2` | Evo 2 7B | [arcinstitute/evo2_7b](https://huggingface.co/arcinstitute/evo2_7b) | `neural`, `neural-evo2` extras; need to be compatible with GPU environment |
| 6 | `generator` | GENERator v2 eukaryote 1.2B | [GenerTeam/GENERator-v2-eukaryote-1.2b-base](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base) | `neural` extra; must be explicit `allow_remote_code=True` |
| 7 | `grover` | GROVER | [PoetschLab/GROVER](https://huggingface.co/PoetschLab/GROVER) | `neural` extra; completed real checkpoint rep + k-means smoke |
| 8 | `hyenadna` | HyenaDNA medium 450k | [LongSafari/hyenadna-medium-450k-seqlen-hf](https://huggingface.co/LongSafari/hyenadna-medium-450k-seqlen-hf) | `neural` extra; must be explicit `allow_remote_code=True` |
| 9 | `janusdna` | JanusDNA 131k | [Harvard Dataverse DOI 10.7910/DVN/HDT0RN](https://doi.org/10.7910/DVN/HDT0RN) | checkpoint to verify the official MD5; [official source code ](https://github.com/Qihao-Duan/JanusDNA) and `model_source_path` must also be provided |
| 10 | `lucaone` (default) | LucaOne gene step 36.8M | [LucaGroup/LucaOne-gene-step36.8M](https://huggingface.co/LucaGroup/LucaOne-gene-step36.8M) | `neural` extra; must be explicit `allow_remote_code=True` |
| 11 | `ntv2` | Nucleotide Transformer v2 500M multi-species | [InstaDeepAI/nucleotide-transformer-v2-500m-multi-species](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) | `neural` extra; must be explicit `allow_remote_code=True` |

Query all currently registered models, checkpoint will not be downloaded:

```python
from dnakit.representations import available_embedding_models

print(available_embedding_models())
```

```text
('alphagenome', 'caduceus', 'dnabert2', 'enformer', 'evo2', 'generator', 'grover', 'hyenadna', 'janusdna', 'lucaone', 'ntv2')
```

### 5.3 Extract rep example

```python
from dnakit import DNARecord, DNASequence
from dnakit.representations import RepresentationConfig, extract_representations

records = [
    DNARecord(DNASequence("ACGTACGT"), "seq-1"),
    DNARecord(DNASequence("AACCGGTT"), "seq-2"),
]
result = extract_representations(
    records,
    config=RepresentationConfig(
        model="lucaone",
        pooling="mean",
        allow_remote_code=True,
    ),
)
print(result.model_name)
print(result.representations.shape)
```

By default, checkpoint is downloaded to `ckpt/lucaone-gene-step36-8m/` in the current working directory, and is reused directly if the complete file already exists. Downloads and sequence-by-sequence extractions display a progress bar by default; this can be adjusted explicitly with `checkpoint_dir`, `checkpoint_path`, or `show_progress`.
