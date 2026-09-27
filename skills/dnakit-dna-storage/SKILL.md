---
name: dnakit-dna-storage
description: Encode bytes into DNA payloads and decode intact packets through the optional Chamaeleo Goldman codec. Use for computational DNA storage demonstrations and codec integration.
license: MIT
---

# DNA storage payloads

Install `pip install "dnakit[dna-storage]"`.

Use `dnakit.genomics.encode_dna_storage(data_hex, segment_bytes=32)`, then preserve the complete returned JSON packet. Decode with `dnakit.genomics.decode_dna_storage(packet)`. The return value contains `data_hex` and an integrity result. Test a round trip on representative bytes before using a new backend version.

The implementation selects Chamaeleo's fixed-Huffman Goldman codec. Segment indices, lengths and digest live in JSON, not in the DNA sequences. Losing that metadata makes DNA-only payloads insufficient for this API. SHA-256 detects accidental corruption; it provides neither error correction nor authenticity. Missing/modified segments are rejected rather than repaired. This does not implement the complete redundant overlapping archive protocol of the original paper or a wet-lab storage workflow.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite [Goldman et al. (2013)](https://doi.org/10.1038/nature11875) and [Chamaeleo: Zhi Ping, Hao-Ling Zhang, Sha Zhu and contributors](https://github.com/ntpz870817/Chamaeleo). DNAKit supplies the packet validation and optional-package adapter.
