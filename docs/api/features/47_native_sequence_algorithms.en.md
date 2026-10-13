# Native k-mer, motif, and sequence-complexity algorithms

These functions operate only on explicit DNA strings. They do not query a reference database or load a trained model. Methods references include SeqKit2 <https://doi.org/10.1002/imt2.191>, MEME Suite <https://doi.org/10.1093/nar/gkp335>, and VSEARCH sequence clustering <https://doi.org/10.7717/peerj.2584>.

## 1) `dnakit.sequence.count_kmers()`

- **Purpose:** Count each k-mer in a sequence.
- **Method:** Scan overlapping windows over normalized DNA and return an order-independent count mapping.
- **API:** `dnakit.sequence.count_kmers(sequence, k=5)`.
- **Input:** A DNA string and positive integer `k`.
- **Output:** A mapping from k-mer to count.
- **Example code:** `dnakit.sequence.count_kmers("ATAT", k=2)`.
- **Example output:** `{"AT": 2, "TA": 1}`.

## 2) `dnakit.sequence.count_canonical_kmers()`

- **Purpose:** Combine forward and reverse-complement k-mer counts.
- **Method:** Use the lexicographically smaller forward or reverse-complement representation for each window.
- **API:** `dnakit.sequence.count_canonical_kmers(sequence, k=5)`.
- **Input:** A DNA string and positive integer `k`.
- **Output:** A mapping from canonical k-mer to count.
- **Example code:** `dnakit.sequence.count_canonical_kmers("ATGC", k=2)`.
- **Example output:** `{"AT": 1, "CA": 1, "GC": 1}`.

## 3) `dnakit.sequence.build_minimizer_index()`

- **Purpose:** Build a bounded minimizer index for a sequence.
- **Method:** Select the lexicographically smallest k-mer in each window and record its start positions.
- **API:** `dnakit.sequence.build_minimizer_index(sequence, k=15, window=10)`.
- **Input:** A DNA string, k-mer length, and minimizer window size.
- **Output:** A mapping from minimizer to position list.
- **Example code:** `dnakit.sequence.build_minimizer_index("ACGTACGT", k=3, window=2)`.
- **Example output:** `{"ACG": [0, 4]}`.

## 4) `dnakit.sequence.generate_kmer_table()`

- **Purpose:** Generate a k-mer table with sequence indices for multiple sequences.
- **Method:** Apply deterministic k-mer counting to each sequence and expand counts into records.
- **API:** `dnakit.sequence.generate_kmer_table(sequences, k=5)`.
- **Input:** A list of DNA strings and `k`.
- **Output:** Records containing `sequence_index`, `kmer`, and `count`.
- **Example code:** `dnakit.sequence.generate_kmer_table(["ATAT"], k=2)`.
- **Example output:** `[{"sequence_index": 0, "kmer": "AT", "count": 2}]`.

## 5) `dnakit.sequence.scan_pwm()`

- **Purpose:** Scan a sequence and return PWM hits above a threshold.
- **Method:** Sum per-position base scores for each window and compare with `threshold`.
- **API:** `dnakit.sequence.scan_pwm(sequence, pwm, threshold=0.0)`.
- **Input:** A DNA string, per-position base-score mappings, and a threshold.
- **Output:** Hits containing position, window sequence, and score.
- **Example code:** `dnakit.sequence.scan_pwm("ATGC", [{"A": 1}, {"T": 1}], threshold=1.5)`.
- **Example output:** `[{"start": 0, "end": 2, "sequence": "AT", "score": 2.0}]`.

## 6) `dnakit.sequence.mask_low_complexity()`

- **Purpose:** Replace low-entropy windows with a mask character.
- **Method:** Calculate sliding-window Shannon entropy and replace windows below the threshold with `N`.
- **API:** `dnakit.sequence.mask_low_complexity(sequence, window=10, entropy_threshold=1.2)`.
- **Input:** A DNA string, window size, and entropy threshold.
- **Output:** The masked DNA string.
- **Example code:** `dnakit.sequence.mask_low_complexity("AAAAAAAA", window=4, entropy_threshold=1.2)`.
- **Example output:** `"NNNNNNNN"`.

## 7) `dnakit.sequence.cluster_sequences_by_identity()`

- **Purpose:** Create deterministic sequence-identity clusters.
- **Method:** Use equal-length positional identity and assign each new sequence to the first representative cluster meeting the threshold.
- **API:** `dnakit.sequence.cluster_sequences_by_identity(sequences, threshold=0.9)`.
- **Input:** A list of DNA strings and an identity threshold from 0 to 1.
- **Output:** A list of input-index lists, one per cluster.
- **Example code:** `dnakit.sequence.cluster_sequences_by_identity(["ACGT", "ACGA"], threshold=0.75)`.
- **Example output:** `[[0, 1]]`.

