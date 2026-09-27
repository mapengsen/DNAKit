# Sequence representations and feature processing

Encode DNA sequences as numerical vectors, then scale, reduce or select features.

## 1) integer encoding {#1-fp-001-integer-encoding}

- **Purpose:** Convert each DNA character into an integer according to the explicit code table, while retaining the original position order, which is convenient for word segmentation, model input and reversible checking of encoding rules.
- **Method:** Encode characters with A=0, C=1, G=2 and T=3 by default; handle other symbols through explicit policies.
- **API:** `dnakit.fingerprints.integer_encode(value[required], ambiguity_policy[optional], gap_policy[optional], max_output_length[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional IUPAC policy, Gap policy and maximum output length.
- **Output:** Returns an `IntegerEncodingResult`; the example shows `values`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import integer_encode

result = integer_encode(DNASequence("ACGT"))
print(result.values)
```

- **Example output:**

```text
(0, 1, 2, 3)
```

## 2) One-hot encoding {#2-fp-002-one-hot-encoding}

- **Purpose:** Convert each base into a four-dimensional vector corresponding to A/C/G/T, forming a "sequence length × 4" position retention matrix, which is convenient for use in statistical models and neural networks.
- **Method:** Create four columns per position in base_order; the matching canonical base is 1 and the other columns are 0.
- **API:** `dnakit.fingerprints.one_hot_encode(value[required], ambiguity_policy[optional], gap_policy[optional], base_order[optional], max_output_length[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional IUPAC, Gap, Column Order and Output Length Policy.
- **Output:** Returns a `OneHotEncodingResult`; the example shows `values`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import one_hot_encode

result = one_hot_encode(DNASequence("ACGT"))
print(result.values[0], result.values[-1])
```

- **Example output:**

```text
(1.0, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 1.0)
```

## 3) k-mer features {#3-fp-003-k-mer-features}

- **Purpose:** Combine the number, frequency or presence of all k-mers into a feature vector in a fixed order for sequence comparison and traditional machine learning; reverse complementary k-mers can be merged into the same feature through `canonical=True`.
- **Method:** Count k-mers in windows; count retains counts, frequency normalizes them and binary records presence.
- **API:** `dnakit.fingerprints.kmer(value[required], k[required], canonical[optional], mode[optional], representation[optional], ambiguity_policy[optional], overlapping[optional], cross_gaps[optional], max_dimension[optional])` is recommended; the old name `kmer_fingerprint(...)` is retained as a compatible alias only.
- **Input:** Required sequence and `k`; optional `canonical`, `mode="count"|"frequency"|"binary"`, storage representation, overlapping policy, `ambiguity_policy="error"|"ignore"` and `cross_gaps`.
- **Output:** Returns k-mer features with `values` and `dimension`.
- **Term:** `count`, `frequency`, and `binary` are classified as k-mer features; `binary` is an existential feature. When a fixed length DNA fingerprint is required, use [Hashed k-mer bit fingerprint](08_feature_engineering.md#1-hashed-k-mer).
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import kmer

sequence = DNASequence("ACGT")
standard = kmer(sequence, k=2, representation="sparse")
canonical = kmer(sequence, k=2, canonical=True, representation="sparse")
print(dict(standard.values), standard.dimension)
print(dict(canonical.values), canonical.dimension)
```

- **Example output:**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 16
{'AC': 2, 'CG': 1} 10
```

## 4) k-mer Sketch (MinHash/FracMinHash) {#4-fp-005-k-mer-sketch-minhashfracminhash}

- **Purpose:** Perform Bottom-k or threshold sampling on the k-mer hash of the sequence to generate a smaller and reproducible Sketch for approximating Jaccard similarity and retrieving large batches of sequences.
- **Method:** MinHash keeps the smallest num_hashes unique hashes; FracMinHash retains hashes below the scaled threshold.
- **API:** `dnakit.fingerprints.minhash(value[required], k[required], num_hashes[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, `dnakit.fingerprints.fracminhash(value[required], k[required], scaled[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`
- **Input:** Required sequence and `k`; optional `num_hashes` for MinHash, optional `scaled` for FracMinHash, optional seed/canonical/upper limit for both.
- **Output:** Returns `SketchResult` with hash values and sketch parameters.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash

sequence = DNASequence("ACGTACGT")
first = minhash(sequence, k=3, num_hashes=3, seed=7)
second = minhash(sequence, k=3, num_hashes=3, seed=7)
print(first.hashes == second.hashes, len(first.hashes))
```

- **Example output:**

```text
True 2
```

<a id="neural-representations"></a>

## 5) Neural network representation

- **Purpose:** Convert each DNA sequence to a fixed-length floating-point vector.
- **Method:** Pool token representations with mean/cls/max/last; split long sequences at the model context limit and average the chunk vectors.
- **API:** `dnakit.representations.extract_representations(records, config=RepresentationConfig(...))`.
- **Input:** Nonempty DNA records without explicit gaps; select a model, checkpoint, pooling and device. The default LucaOne backend requires `allow_remote_code=True` and model dependencies.
- **Output:** Returns a `RepresentationResult`; the example shows `model_name`, `representations`.
- **Example code:**

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

- **Example output:** `model_name` is `lucaone`; the matrix has two rows and a checkpoint-dependent number of columns. This example has not been run with model weights here.

[FAQ: details](../../faq/tutorial_basics.md#embedding-models) · [FAQ: details](../../faq/tutorial_basics.md#embedding-loading)

<span id="51-representation-api"></span>
<span id="52-rep-of-11-types-of-neural-networks-supported"></span>
<span id="53-extract-rep-example"></span>

## 6) Encode DNA sequences as a feature matrix {#bio-23}

- **Purpose:** Create numeric features for classification, clustering or statistics.
- **Method:** Compute the requested descriptor; NAC gives A/C/G/T frequencies and RCKmer merges reverse-complement k-mer frequencies.
- **API:** `dnakit.bioinformatics.encoding.encode_dna(sequences, descriptor, parameters=None, property_tables=None, labels=None, training_ids=None)`
- **Input:** `DNA([...])` or an ID→A/C/G/T sequence mapping; requires `bio-analysis`. Positional encodings require equal lengths; physicochemical encodings require explicit property tables; see [all encodings and parameters](../../faq/bioinformatics.md#encoding).
- **Output:** Returns `EncodingResult` with sequence `ids` and numeric feature matrix `values`.
- **Example code:**

```python
from dnakit import DNA
from dnakit.bioinformatics.encoding import encode_dna

sequences = DNA(["AACGTG", "CCGTAT"])
nac = encode_dna(sequences, "NAC")
rckmer = encode_dna(sequences, "RCKmer", parameters={"k": 2})
print(nac.ids)
print(nac.values.round(3).tolist())
print(rckmer.values.shape)
```

- **Example output:**

```text
('sequence_1', 'sequence_2')
[[0.333, 0.167, 0.333, 0.167], [0.167, 0.333, 0.167, 0.333]]
(2, 10)
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#encoding).

## 7) Scale, reduce and select features {#bio-24}

- **Purpose:** Fit a feature transformation on training data and reuse it for new samples.
- **Method:** Scale features, then apply a method such as PCA; supervised selection uses `labels`, and new samples use `transform`.
- **API:** `dnakit.bioinformatics.features.analyze_features(values, method='pca', labels=None, feature_names=None, scale='standard', components=2, selected_features=10)`
- **Input:** A samples×features matrix; requires `bio-analysis`. Training and new data must share column order; t-SNE has no `transform`; see [available methods](../../faq/bioinformatics.md#encoding).
- **Output:** Returns `FeatureAnalysis` with transformed `values`; `transform` handles new samples.
- **Example code:**

```python
from dnakit.bioinformatics.features import analyze_features

fitted = analyze_features([[1, 2], [2, 1], [3, 4]],
                          feature_names=["x", "y"], method="pca", components=1)
held_out = fitted.transform([[2, 3]], feature_names=["x", "y"])
print(fitted.values.shape, held_out.shape)
```

- **Example output:**

```text
(3, 1) (1, 1)
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#encoding).
