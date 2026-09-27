# Phylogenetics and Molecular Evolution: parameters and notes

See the [tutorial](../../api/features/30_phylogeny.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Calculate dN/dS between two coding sequences {#dnakit.databases.sequence_dn_ds-parameters}

`dnakit.databases.sequence_dn_ds(*, seq1=None, seq2=None, fasta1_path=None, fasta2_path=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `seq1` | string &#124; null | None | First coding sequence (in-frame, codon-aligned to seq2). DNA/RNA letters. |
| `seq2` | string &#124; null | None | Second coding sequence (same length/frame as seq1). |
| `fasta1_path` | string &#124; null | None | Alternative to seq1: path to a single-record FASTA |
| `fasta2_path` | string &#124; null | None | Alternative to seq2: path to a single-record FASTA |
| `verbose` | bool | True | Progress flag. |
