# minimap2 and paftools: mapping, alignment conversion, assembly and annotation evaluation

This group of functions in `dnakit.alignment` is a set of typed adapters over [minimap2](https://github.com/lh3/minimap2) (tested with 2.31) and its companion scripts. `minimap2_map_reads()` runs the aligner, including the spliced RNA/cDNA presets and the all-vs-all read overlap presets. The other functions run one sub-command each of `paftools.js` (executed by the [k8](https://github.com/attractivechaos/k8) JavaScript shell) or `pafcluster.js`: they convert between PAF, SAM, BED12, GFF3, MUMmer delta and simulated-read formats, report alignment, assembly, gene-completeness and misjoin statistics, lift coordinates, call and compare variants, and evaluate mappings, splice junctions, exons and read overlaps. Each function runs its program with a fixed, typed list of options and without a shell. They do not expose every option, and `paftools.js bedcov` is not wrapped.

minimap2 and paftools.js are distributed under the MIT licence. DNAKit runs them as separate processes and copies none of their code.

The functions are in the module `dnakit.alignment.minimap2_p09_cli_functions`. The programs are not Python packages: install them with `conda install -c bioconda minimap2 k8` (this puts `minimap2`, `paftools.js` and `k8` on PATH) or from the release tarball (`misc/paftools.js`, `misc/pafcluster.js`), and pass `minimap2_path`, `paftools_path`, `k8_path` or `pafcluster_path` when they are not on PATH. A missing program raises `BackendUnavailableError` with an install hint before any output directory is created. Every function writes into a new `output_directory` (it must not exist yet): programs that print to standard output are streamed into a named result file, standard error goes to a log file, and a failing program raises `BackendExecutionError` with the log tail. Inputs are validated (`ConfigurationError`). Optional parameters default to `None` and are then not passed, so the tool defaults apply. `paftools.js call` and `ov-eval` need a PAF sorted by target name and start; the two facades sort it themselves and keep the sorted copy.

## Example data

All examples read these small files from the working directory: `ref.fa` has the reference contigs `chrA` (3000 bases, with three exons at 300-450, 800-950 and 1300-1450 and GT...AG introns) and `chrB` (2000 bases), `cdna.fa` has the spliced transcript `tx1` (450 bases) and a variant of it with one substitution, and `asm.fa` is an assembly of the reference with four substitutions and a 12-base deletion in `chrA` and a 10-base insertion in `chrB`. Smaller inputs (BED, VCF, GTF, delta, simulator files) are written by the examples themselves, and every example puts its output into a new temporary directory. The examples assume that `minimap2`, `k8` and `paftools.js` (and, for the last example, `pafcluster.js`) are on PATH.

**Example file:** `ref.fa`

```text
>chrA
CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGA
CCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCCTAGCCTAGCTACGCGTTTGCGCATCAGGCTGT
CCCATACATCAAGCGGTTCCCCTCAAATTATCCGGACTCGGTAAGGGCAGCGAGTAAATATTTTACAATA
CGTTTCTTGTCAATCTGCTGCTTTGTACGCGTCACAGTTACTCGGCGAAGGCCCGTCTTTTTGCTGACCA
GGAAATTTCACAGCTGAGCCTAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCC
TAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGA
TCAGTGCCCTCTAAGTCTCTAAGCTGTGTTGTAGGTACGAGCCCGAGCCCTTCAGGACCGAGTAAACTTG
TAGCGTTTCTCATCAGTCCAGGGGCATCCCACCCACATAACCAACCACCTATGGGTATATTCAAGTGCGG
GTGTGAAGATGCCGGTAGTCAGTATCGCATGGTCATCCACCCGACTCGTCGCGTCGGCGAACGGTCTAGG
CCAACTCTCCTTGCTACAACTATAAGACGTGTTAGGATGTGGGCGGCCAGCAGACGCAAACGCCGCCACG
TGGCTTGACGGCGTCATTCCTATTATCAAAGCAATATGTTTGCGCGACCTGGGTAGAACCTGTGCTGCGG
TTCGCCCACGTTGCGAAGACCACTTTGCAGAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACC
CTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGA
CCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGTACATTGTAGAGGGCACCAGCTCAAGTAG
TTTGATGCGCACTTCAATGACCGACCTCCGTGACCGGTGGGTGCCGAGCGTCGGCAAATCTTTCTTCCAG
TATGGTAGAGCGTCAGGTGTCTGTCCCTATCAGCTGTTGACTGCAGGGGACATGTACAACCTATCCATAA
TTACTTCCAAGCCGAGCTGTAGTTGATCACGTTAGGTAAGCTGTAGACGTACGGATTACGACCCCTTAGT
CCAGATAGGGACGCTGCTGAGCTCTTAGGCTTTCCATTGTAGAGCCAACAGATAGACGTAGTCTACTCTC
GCCCCCGGCTTAACTATAGCGTCTCAGAATACACGGGGAGGGTAGATTGAGAGGATTGTCAACAACAATT
GTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTA
CCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTGAAGGGACTCAGTCGACTATA
CAGAATCATCCTCGTGAAAAATGTCTCATCTCGTGACTTAAACGAGCAGTAACCGCCGAACGACAAGGAT
TGACGACACCGTCTCGCTTGTGCTCGAGGTTTGCCCAGCATTCCCTAAACTCCCTTCGTACAGTGGCGCT
AAGCTTAACGCTTGAGCTACGGATTACAATATCGTGAGCCCATTTGAACAGCGCGTGGCTCGATCCGCGA
TGGAGGTTTGCGGTACCCGCAATGCATCGCCAAGCGTTCCCATCCTTCGAGACCTGAGTGAAGCTTCCTG
GCGGCCTATATATCTAAGAGAGGGGAAAGGAATGCCAACCTGGACTGCCGGGAGGGCAGTTATCAGGTTA
TTGGTTGCAACATTAACCGCAGTCGGAGTCTTAACGTAATCATATTACGAAGTATCGCATAATGCGACAA
ATAATGTCGGAGTCTTTACCAGGCTTCACCTATACACGCTACAGTTATCTGAGAATCTGCCAGCTCTAGA
GTTTCCCATCGCTCCTGCCACTTAGACCGAATCTGCTTAATGTGAGGGTGCTAGTGCAGTTGGGAGACCT
GGTGTCTACCCTGTGGTAAAAGCAGTAGTCGATAAGGCACAGTAGGCCCGAGGACAACCCCCCCATACTC
CGGAGGTGGAGGCTAACGATATCATTTCGCATATAGATCATTCCTGGTAGGAGCGCGAAATCGTCGCTTC
TTGTAGGCAATTCATGCACAACTAAAGAAATTATACAGGCCCAAGCTGTAGAGTTCACGTCGCGTGGCAC
CACTCTCTAATTAACAGTATAATTTTTTTCTCTATTGAGGTGAAGCGATGTGCACTCCACTACTCCACAA
GTTCAGGCGACAATATGCATTCGTAATCCTGTTGCTCGTGGGAACCGGCGCCTTGCATTTTAGTCAGAAG
GGTTGGAAAACAAGGGGTCCTGCACATGTGACGTTAATTCTGTCTCCCACATCTCGACTCATAATGGTGC
CTGTTTTGTGGACTGTGCGAAAGGCTCTTGCTGGCATATGATTCAAACCTATCACCAGATCTCAGCTGAT
GCGACTGTGCTAGTCCATCTGCCCAAGGGCATTTCCCAATACTGATACGCTTACAGCTTGTTGAAAGTAA
AGTGCTAACATACTAGCGGCACGACACAAGGGCGGAAAAACGTTGTTGTCGTTCACTCCCGTACCTACTG
GGGCTTCTAGCCCCAACTAGGGTTGTAGCCGCATGGCCGTTCTCCCGAGTTATTAGCTGGACCATCGGAG
CATCAGACCCAATCATCAAACCAGTACGCGAATAAGACGGGATCGAAGGTTGGTGTACGATTCTGTTCCG
TACGATATACACGACGAGAAAGCTCAGCGGCAACGGTTAACCTGTGGTAGCGTCGAGGGTTTTAGGGCCC
TGTGGAATCATGTAATGTTATCAAGAACAAATGATCCACAAGTGAGAAAACCCGTGAAGTCTACTCCTCT
GTCACAGAAGGCCCCAGGATAAAAAGTCTAGTACCACGCACCAGGTTGGTGGGGGGTGAG
>chrB
CTGGACTTTTATCTGAGGATGCGGTAGCCCGCTCACAGCAATCTATTCTCTAGCGGTCGAAGTGGCTGGA
GCCGCCTTGCTTGTTAGTCGTGATGCCCTCTTTCAGAGTGTCCAAGTAAAAACCAAGAGCACCTTATTCA
TCTTACCCGTTGCGGGCCGACTTCAATTCGTGGTTCTCAGCAATATTTCAGTGGGTTATTTTTTGAGCGG
TTATGGGACCGTGCCGAAAGAACGAATTTAGCGAGTTCGCGCGACTCACAAGGACCAAGATCTGCAAATC
CGTTAAGATGAGCGGAAGGAAGGTCTTGGATGGAGCTTAAGGGTGGGAGAGGGTTGTTACTAGAAGGAGA
TGTCTCCGTGGTTCGCTTGCTAAAGGGCCAGACATGTTATCACGCGACCAGGGCTGCGGCCCAACACTGT
AAACCCTGCAAATACCATGATGGCTCCGTGAGGCTAATGACATTACAAATCAGTATATGTCGTGCTATTT
CGTGTAGTCCATTGACCACGTAGAGAAATGGATGATTAGTACTAGTTTTCTCTGTTCCGACCAGCTTGGC
CCACTTTCGATCGCCCAACAGAGAGGTGATGATATCGAGGTAGGTGGCGGATCTGATCTGTACAATTGGA
CTGTTGAATAACTAGTCAGTAAGTCTACAAAGGAACCTCGAGTCGAGTACGAGCCAACCGTGTATCATGC
AGGGCTAACAGCTACTAGCTTTCCCCGTCGACACTTTTACAATGCCTCAGACCCGGCCGGTCGAGCAGCG
GATATGAGGAGGGGATCGTTTACTGAAACTCACGGATTAAATGTTTTGATATTAGGCGTGCATCTGCCTG
GTCAGCAAAACAACTTAAACGAGAAGCTGGTGCGCAGATAACCCCTCGTTATGAATAGCGAGCACCTCCG
ACAAGTTGAATAGCCGCCTGAAACTGTAGAGAACACTAGTTACGTGGTAGGGGATTTAAGGGCACCGGTC
ATCGAACATCGTAAATGCAGGAATGTTAAATGAAGTACTCTGAGTTTCCAGCCGACATGAGAGCTCGTCC
GGCTTAATGTTTAACCATCGCGGTGGCTTTGTCGAACAGACGTCTCATCGGACCGGTGAGCACCAAAATC
TGTTTGTATACCTAATCTAGACTGGTCTGCAGCGCCACTCTGAGAATAAAGAGAAATCGTTGAGGACAAG
TCTAGAAGTGATGGAGTACTATAGTTCTGTAGATGGAGAGGTTGAAGGACTATAGGATACATTGGCGAAA
TATAGTAAGTGTATGGGAGACAAGGAAAACGGGACGACTTCACCACCAAGGCGTACCTCCCACATTTCTC
GCAACTAGTCTCGGTTCTAAAGCTTGTCCCATAGATCTACGATGGTACGGACACGCGGGGGTAACTTATT
TTTTGATCATGCCTATGGCCACCACCCGTTCTATTGGTCCATAAGATTGGGAATTCTCCGTTTGTGGGAC
ATAGGCGCTCCGCCATTTAAGGCTCCAATGCTCACTAAGATATAAGACGTCAGGACTGCGGTAATGGTGA
CGCTATAGGGGCCCTTCTAACTCGTAAAACGCTGTACCGCGCCGGTACATGCTTCTGGAAAAATCTCTTA
CGACATAAGATCAACACTCGCATCAATATGTAGTCGGACTGATGCGTCACCCATGGGCACATCGCGATAG
AATTTGATCGAAGTACGAGTCTATCTCATGTCTGTGTGGATTGAATAGATACCCCCAGTATGAAGTGGAT
CTCGGATTAGCTGAGTCGTCTGAGCAGCTAACCCAGGGAGTGGCCCCAGGTTATGGTTGTCCCCCCCGGT
ATTCCACCATACCTCTGTCATGAGTTATGGAGCTAAGTATGGTAGCCTTACATACTCATTCGCCACGACG
TGGGAAGTCGTCCGAAGCCTTGAGAACGTAACGGTCTATTCCGACACTTAGCAATTTCTCCAGTGGGCGT
CTAGTGCGTTAGTAACTAGGCCGTACGCGTGGCCCAGGGG
```

**Example file:** `cdna.fa`

```text
>tx1
TAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCCTAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGATCAGTGCCCTCTAAGTCTCTAAGCTGTGTTAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACCCTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGACCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGGTAGATTGAGAGGATTGTCAACAACAATTGTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTACCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTG
>tx1_variant
TAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCCTAGCCAGCAGATCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGATCAGTGCCCTCTAAGTCTCTAAGCTGTGTTAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACCCTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGACCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGGTAGATTGAGAGGATTGTCAACAACAATTGTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTACCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTG
```

**Example file:** `asm.fa`

```text
>ctgA
CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGA
CCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCCTAGCCTAGCTACGCGTTTGCGCATCAGGCTGT
CCCATACATCAAGCGGTTCCCCTCAAATTATCCGGACTCGGTAAGGGCAGCGAGTAAATATTTTACAATA
CGTTTCTTGTCAATCTGCTGCTTTGTACGCGTCACAGTTACTCGGCGAAGGCCCGTCTTTTTGCTGACCA
GGAAATTTCACAGCTGAGCCTAGCTTCCTAAATCCATTTGCGCGGGAAACACGGGACATGTCAACGGTCC
TAGCCAGCAGTTCTAGACAGTCTGAGCGATCCTCCGTGACTCGGCATACACGGACCTTTCCGCTTCTTGA
TCAGTGCCCTCTAAGTCTCTAAGCTGTGTTGTAGGTACGAGCCCGAGCCCTTCAGGACCGAGTAAACTTG
TAGCGTTTCTAATCAGTCCAGGGGCATCCCACCCACATAACCAACCACCTATGGGTATATTCAAGTGCGG
GTGTGAAGATGCCGGTAGTCAGTATCGCATGGTCATCCACCCGACTCGTCGCGTCGGCGAACGGTCTAGG
CCAACTCTCCTTGCTACAACTATAAGACGTGTTAGGATGTGGGCGGCCAGCAGACGCAAACGCCGCCACG
TGGCTTGACGGCGTCATTCCTATTATCAAAGCAATATGTTTGCGCGACCTGGGTAGAACCTGTGCTGCGG
TTCGCCCACGTTGCGAAGACCACTTTGCAGAGTTCGTTGCAGGGGTAGCCCAGCCCGAAATCCTTGTACC
CTTCTTGACCAAGAATGTTGATAGCGTAGATTTTTGATCAGCTGGTCTTGCGGGGTAAGAATTGTGTCGA
CCACGTTTACATTCTCCGCACCAGTCATGAGACTTTGAGTGTACATTGTAGAGGGCACCAGCTCAAGTAG
TTTGATGCGCACTTCAATGACCGACCTCCGTGACCGGTGGGTGCCGAGCGTCGGCAAATCTTTCTTCCAG
TATGGTAGAGCGTCAGGTGTCTGTCCCTATCAGCTGTTGACTGCAGGGGAAATGTACAACCTATCCATAA
TTACTTCCAAGCCGAGCTGTAGTTGATCACGTTAGGTAAGCTGTAGACGTACGGATTACGACCCCTTAGT
CCAGATAGGGACGCTGCTGAGCTCTTAGGCTTTCCATTGTAGAGCCAACAGATAGACGTAGTCTACTCTC
GCCCCCGGCTTAACTATAGCGTCTCAGAATACACGGGGAGGGTAGATTGAGAGGATTGTCAACAACAATT
GTTTTCCACACTGCCAAATCTTACGTCCATTCACAATTAATATAGCACCAGGTCCAATGTTTCGGGTGTA
CCGAGGCTGAATCGTTGACACTATCATGTTCTAAGGATCGGACACAACTGAAGGGACTCAGTCGACTATA
CAGAATCATCCTCGTGAAAAATGTCTCATCTCGTGACTTAAACGAGCAGTAACCGCCGAACGACAAGGAT
TGACGACACCGTCTCGCTTGTGCTCGAGGTTTGCCCAGCATTCCCTAAACTCCCTTCGTACAGTGGCGCT
AAGCTTAACGCTTGAGCTACGGATTACAATATCGTGAGCCCATTTGAACAGCGCGTGGCTCGATCCGCGA
TGGAGGTTTGCGGTACCCGCCATGCATCGCCAAGCGTTCCCATCCTTCGAGACCTGAGTGAAGCTTCCTG
GCGGCCTATATATCTAAGAGAGGGGAAAGGAATGCCAACCTGGACTGCCGGGAGGGCAGTTATCAGGTTA
TTGGTTGCAACATTAACCGCAGTCGGAGTCTTAACGTAATCATATTACGAAGTATCGCATAATGCGACAA
ATAATGTCGGAGTCTTTACCAGGCTTCACCTATACACGCTACAGTTATCTGAGAATCTGCCAGCTCTAGA
GTTTCCCATCGCTCCTGCCACTTAGACCGAATCTGCTTAAAGTGCAGTTGGGAGACCTGGTGTCTACCCT
GTGGTAAAAGCAGTAGTCGATAAGGCACAGTAGGCCCGAGGACAACCCCCCCATACTCCGGAGGTGGAGG
CTAACGATATCATTTCGCATATAGATCATTCCTGGTAGGAGCGCGAAATCGTCGCTTCTTGTAGGCAATT
CATGCACAACTAAAGAAATTATACAGGCCCAAGCTGTAGAGTTCACGTCGCGTGGCACCACTCTCTAATT
AACAGTATAATTTTTTTCTCTATTGAGGTGAAGCGATGTGCACTCCACTACTCCACAAGTTCAGGCGACA
ATATGCATTCGTAATCCTGTTGCTCGTGGGAACCGGCGCCTTGCATTTTAGTCAGAAGGGTTGGAAAACA
AGGGGTCCAGCACATGTGACGTTAATTCTGTCTCCCACATCTCGACTCATAATGGTGCCTGTTTTGTGGA
CTGTGCGAAAGGCTCTTGCTGGCATATGATTCAAACCTATCACCAGATCTCAGCTGATGCGACTGTGCTA
GTCCATCTGCCCAAGGGCATTTCCCAATACTGATACGCTTACAGCTTGTTGAAAGTAAAGTGCTAACATA
CTAGCGGCACGACACAAGGGCGGAAAAACGTTGTTGTCGTTCACTCCCGTACCTACTGGGGCTTCTAGCC
CCAACTAGGGTTGTAGCCGCATGGCCGTTCTCCCGAGTTATTAGCTGGACCATCGGAGCATCAGACCCAA
TCATCAAACCAGTACGCGAATAAGACGGGATCGAAGGTTGGTGTACGATTCTGTTCCGTACGATATACAC
GACGAGAAAGCTCAGCGGCAACGGTTAACCTGTGGTAGCGTCGAGGGTTTTAGGGCCCTGTGGAATCATG
TAATGTTATCAAGAACAAATGATCCACAAGTGAGAAAACCCGTGAAGTCTACTCCTCTGTCACAGAAGGC
CCCAGGATAAAAAGTCTAGTACCACGCACCAGGTTGGTGGGGGGTGAG
>ctgB
CTGGACTTTTATCTGAGGATGCGGTAGCCCGCTCACAGCAATCTATTCTCTAGCGGTCGAAGTGGCTGGA
GCCGCCTTGCTTGTTAGTCGTGATGCCCTCTTTCAGAGTGTCCAAGTAAAAACCAAGAGCACCTTATTCA
TCTTACCCGTTGCGGGCCGACTTCAATTCGTGGTTCTCAGCAATATTTCAGTGGGTTATTTTTTGAGCGG
TTATGGGACCGTGCCGAAAGAACGAATTTAGCGAGTTCGCGCGACTCACAAGGACCAAGATCTGCAAATC
CGTTAAGATGAGCGGAAGGAAGGTCTTGGATGGAGCTTAAGGGTGGGAGAGGGTTGTTACTAGAAGGAGA
TGTCTCCGTGGTTCGCTTGCTAAAGGGCCAGACATGTTATCACGCGACCAGGGCTGCGGCCCAACACTGT
AAACCCTGCAAATACCATGATGGCTCCGTGAGGCTAATGACATTACAAATCAGTATATGTCGTGCTATTT
CGTGTAGTCCATTGACCACGTAGAGAAATGGATGATTAGTACTAGTTTTCTCTGTTCCGACCAGCTTGGC
CCACTTTCGATCGCCCAACAGAGAGGTGATGATATCGAGGTAGGTGGCGGATCTGATCTGTACAATTGGA
CTGTTGAATAACTAGTCAGTAAGTCTACAAAGGAACCTCGAGTCGAGTACGAGCCAACCGTGTATCATGC
AGGGCTAACAGCTACTAGCTTTCCCCGTCGACACTTTTACAATGCCTCAGACCCGGCCGGTCGAGCAGCG
GATATGAGGAGGGGATCGTTTACTGAAACTCACGGATTAAATGTTTTGATATTAGGCGTGCATCTGCCTG
GTCAGCAAAACAACTTAAACGAGAAGCTGGTGCGCAGATAACCCCTCGTTATGAATAGCGAGCACCTCCG
ACAAGTTGAATAGCCGCCTGAAACTGTAGAGAACACTAGTTACGTGGTAGGGGATTTAAGGGCACCGGTC
ATCGAACATCGTAAATGCAGACGTACGTACGAATGTTAAATGAAGTACTCTGAGTTTCCAGCCGACATGA
GAGCTCGTCCGGCTTAATGTTTAACCATCGCGGTGGCTTTGTCGAACAGACGTCTCATCGGACCGGTGAG
CACCAAAATCTGTTTGTATACCTAATCTAGACTGGTCTGCAGCGCCACTCTGAGAATAAAGAGAAATCGT
TGAGGACAAGTCTAGAAGTGATGGAGTACTATAGTTCTGTAGATGGAGAGGTTGAAGGACTATAGGATAC
ATTGGCGAAATATAGTAAGTGTATGGGAGACAAGGAAAACGGGACGACTTCACCACCAAGGCGTACCTCC
CACATTTCTCGCAACTAGTCTCGGTTCTAAAGCTTGTCCCATAGATCTACGATGGTACGGACACGCGGGG
GTAACTTATTTTTTGATCATGCCTATGGCCACCACCCGTTCTATTGGTCCATAAGATTGGGAATTCTCCG
TTTGTGGGACATAGGCGCTCCGCCATTTAAGGCTCCAATGCTCACTAAGATATAAGACGTCAGGACTGCG
GTAATGGTGACGCTATAGGGGCCCTTCTAACTCGTAAAACGCTGTACCGCGCCGGTACATGCTTCTGGAA
AAATCTCTTACGACATAAGATCAACACTCGCATCAATATGTAGTCGGACTGATGCGTCACCCATGGGCAC
ATCGCGATAGAATTTGATCGAAGTACGAGTCTATCTCATGTCTGTGTGGATTGAATAGATACCCCCAGTA
TGAAGTGGATCTCGGATTAGCTGAGTCGTCTGAGCAGCTAACCCAGGGAGTGGCCCCAGGTTATGGTTGT
CCCCCCCGGTATTCCACCATACCTCTGTCATGAGTTATGGAGCTAAGTATGGTAGCCTTACATACTCATT
CGCCACGACGTGGGAAGTCGTCCGAAGCCTTGAGAACGTAACGGTCTATTCCGACACTTAGCAATTTCTC
CAGTGGGCGTCTAGTGCGTTAGTAACTAGGCCGTACGCGTGGCCCAGGGG
```

## 1) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_map_reads()`

- **Purpose:** Align reads or sequences to a reference with minimap2, including spliced RNA/cDNA alignment and all-vs-all read overlaps.
- **Method:** Run `minimap2 -x <preset>`. `preset` is one of `map-ont`, `map-hifi`, `map-pb`, `lr:hq`, `sr`, `asm5`, `asm10`, `asm20`, the spliced presets `splice`, `splice:hq`, `splice:sr` (with `splice_strand` = `-u`, `max_intron_length` = `-G`, `junction_bed` = `--junc-bed`) and the overlap presets `ava-pb`, `ava-ont` (give the same file as `reference` and `query`). `mate` adds the second file of paired short reads, `output_format='sam'` adds `-a`, `base_level_alignment` adds `-c`, `cs_tag`, `md_tag`, `eqx` add `--cs`, `--MD`, `--eqx`; `kmer_size`, `window_size`, `secondary`, `max_secondary`, `soft_clip_supplementary` and `threads` map to `-k`, `-w`, `--secondary`, `-N`, `-Y`, `-t`. Standard output is streamed into the result file and standard error into `minimap2.log`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_map_reads(reference: 'str | Path', query: 'str | Path', output_directory: 'str | Path', *, preset: "Literal['map-ont', 'map-hifi', 'map-pb', 'lr:hq', 'splice', 'splice:hq', 'splice:sr', 'sr', 'asm5', 'asm10', 'asm20', 'ava-pb', 'ava-ont']" = 'map-ont', mate: 'str | Path | None' = None, output_format: "Literal['paf', 'sam']" = 'paf', base_level_alignment: 'bool' = False, cs_tag: "Literal['short', 'long'] | None" = None, md_tag: 'bool' = False, eqx: 'bool' = False, kmer_size: 'int | None' = None, window_size: 'int | None' = None, max_intron_length: 'int | str | None' = None, splice_strand: "Literal['f', 'r', 'b', 'n'] | None" = None, junction_bed: 'str | Path | None' = None, secondary: 'bool | None' = None, max_secondary: 'int | None' = None, soft_clip_supplementary: 'bool' = False, threads: 'int | None' = None, minimap2_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Reference and query files (FASTA/FASTQ) and a new output directory. The tools are called from `PATH`; pass `minimap2_path`, `paftools_path` and `k8_path` for other locations.
- **Output:** `outputs` (`alignments` = `alignments.paf` or `.sam`, `log`), `alignment_count`, `preset`, `version` and, for PAF, `records` (12 standard columns plus a `tags` dictionary).
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads

result = minimap2_map_reads("ref.fa", "cdna.fa", tempfile.mkdtemp() + "/run", preset="splice:hq", base_level_alignment=True, verbose=False)
first = result["records"][0]
print(result["alignment_count"], first["query"], first["target"], first["target_start"], first["target_end"])
print(first["tags"]["cg"])
```

- **Example output:**

```text
2 tx1 chrA 300 1450
150M350N150M350N150M
```

## 2) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_view_paf()`

- **Purpose:** Render PAF alignments as BLAST-like text, MAF or lastz-cigar.
- **Method:** Run `paftools.js view -f <output_format>` through `k8`. The PAF needs the long `cs` tag. `line_length` is `-l` for the BLAST-like layout.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_view_paf(paf: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['aln', 'maf', 'lastz-cigar']" = 'aln', line_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A PAF file with `--cs=long` and a new output directory.
- **Output:** `outputs` (`output` = `view.txt`, `log`) and `output_format`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_view_paf

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="long", verbose=False)["outputs"]["alignments"]
result = minimap2_view_paf(asm, work + "/view", line_length=60, verbose=False)
lines = open(result["outputs"]["output"]).read().splitlines()
print(lines[1][:20], lines[3][:20])
```

- **Example output:**

```text
Ref+:          1 CTT Qry+:          1 CTT
```

## 3) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_sam_to_paf()`

- **Purpose:** Convert SAM alignments to PAF.
- **Method:** Run `paftools.js sam2paf`. `alignments` is `all`, `primary_and_supplementary` (`-p`) or `primary` (`-P`); `include_unmapped` is `-U` and `long_cs` is `-L`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_sam_to_paf(sam: 'str | Path', output_directory: 'str | Path', *, alignments: "Literal['all', 'primary_and_supplementary', 'primary']" = 'all', include_unmapped: 'bool' = False, long_cs: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A SAM file and a new output directory.
- **Output:** `outputs` (`output` = `converted.paf`), `rows` (first lines) and `row_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_sam_to_paf

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
result = minimap2_sam_to_paf(splice_sam, work + "/paf", alignments="primary", verbose=False)
print(result["row_count"], result["rows"][0][:9])
```

- **Example output:**

```text
2 ['tx1', '450', '0', '450', '+', 'chrA', '3000', '300', '1450']
```

## 4) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_delta_to_paf()`

- **Purpose:** Convert MUMmer delta alignments to PAF.
- **Method:** Run `paftools.js delta2paf`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_delta_to_paf(delta: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A `.delta` file (nucmer output) and a new output directory.
- **Output:** `outputs` (`output` = `converted.paf`), `rows` and `row_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_delta_to_paf

work = tempfile.mkdtemp()
delta = work + "/toy.delta"
open(delta, "w").write("ref.fa qry.fa\nNUCMER\n>chrX qryX 1000 900\n1 900 1 900 0 0 0\n0\n")
result = minimap2_delta_to_paf(delta, work + "/paf", verbose=False)
print(result["rows"][0][:12], result["rows"][0][12:])
```

- **Example output:**

```text
['qryX', '900', '0', '900', '+', 'chrX', '1000', '0', '900', '900', '900', '0'] ['NM:i:0', 'cg:Z:900M']
```

## 5) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_spliced_alignments_to_bed12()`

- **Purpose:** Convert spliced alignments in PAF or SAM to BED12 gene structures.
- **Method:** Run `paftools.js splice2bed`; `keep_multiple` is `-m` (keep secondary mappings). PAF input needs the `cg` tag (`minimap2 -c`).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_spliced_alignments_to_bed12(alignments: 'str | Path', output_directory: 'str | Path', *, keep_multiple: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A spliced PAF or SAM file and a new output directory.
- **Output:** `outputs` (`output` = `structures.bed`), `records` (the 12 BED12 columns) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_spliced_alignments_to_bed12

work = tempfile.mkdtemp()
splice = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_spliced_alignments_to_bed12(splice, work + "/bed", verbose=False)
record = result["records"][0]
print(record["name"], record["start"], record["end"], record["block_count"], record["block_sizes"], record["block_starts"])
```

- **Example output:**

```text
tx1 300 1450 3 150,150,150, 0,500,1000,
```

## 6) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_annotation_to_bed12()`

- **Purpose:** Convert GTF or GENCODE-style GFF3 gene annotation to BED12.
- **Method:** Run `paftools.js gff2bed`. `junctions_only` is `-j`, `short_names` is `-s`, `fai` is `-u`, `ensembl_canonical_only` is `-e` and `output_gff` is `-g` (with `fai`). GFF3 input needs `transcript_id=` attributes as GENCODE/Ensembl write them; other GFF3 is rejected by the tool and raises `BackendExecutionError`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_annotation_to_bed12(annotation: 'str | Path', output_directory: 'str | Path', *, junctions_only: 'bool' = False, short_names: 'bool' = False, fai: 'str | Path | None' = None, ensembl_canonical_only: 'bool' = False, output_gff: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A GTF or GFF3 file and a new output directory.
- **Output:** `outputs` (`output` = `annotation.bed`, or `annotation.gff`), `rows` and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_annotation_to_bed12

work = tempfile.mkdtemp()
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_annotation_to_bed12(gtf, work + "/bed", verbose=False)
print(result["rows"][0])
junctions = minimap2_annotation_to_bed12(gtf, work + "/junctions", junctions_only=True, verbose=False)
print(junctions["rows"])
```

- **Example output:**

```text
['chrA', '300', '1450', 't1||g1', '1000', '+', '300', '1450', '196,196,196', '3', '150,150,150,', '0,500,1000,']
[['chrA', '450', '800', 't1||g1', '1000', '+'], ['chrA', '950', '1300', 't1||g1', '1000', '+']]
```

## 7) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_gff3_to_junctions()`

- **Purpose:** Extract intron junctions from a GFF3 file as BED.
- **Method:** Run `paftools.js gff2junc`. Introns are the gaps between consecutive `feature` records (`-f`, tool default `CDS`) that share a `Parent=`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_gff3_to_junctions(gff3: 'str | Path', output_directory: 'str | Path', *, feature: 'str | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A GFF3 file, a new output directory and optionally `feature`.
- **Output:** `outputs` (`output` = `junctions.bed`), `records` (chrom, start, end, name, score, strand) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_gff3_to_junctions

work = tempfile.mkdtemp()
gff = work + "/genes.gff3"
rows = ["##gff-version 3", "chrA\ttest\tmRNA\t301\t1450\t.\t+\t.\tID=t1"]
rows += [f"chrA\ttest\tCDS\t{s + 1}\t{e}\t.\t+\t0\tParent=t1" for s, e in [(300, 450), (800, 950), (1300, 1450)]]
open(gff, "w").write("\n".join(rows) + "\n")
result = minimap2_gff3_to_junctions(gff, work + "/junctions", verbose=False)
print([(row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
[(450, 800), (950, 1300)]
```

## 8) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_longcs_to_sequences()`

- **Purpose:** Reconstruct the aligned target or query sequences from a long-cs PAF.
- **Method:** Run `paftools.js longcs2seq`; `query=True` is `-q` (query instead of target). The PAF needs `--cs=long`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_longcs_to_sequences(paf: 'str | Path', output_directory: 'str | Path', *, query: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A long-cs PAF file and a new output directory.
- **Output:** `outputs` (`output` = `sequences.fasta`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_longcs_to_sequences

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="long", verbose=False)["outputs"]["alignments"]
result = minimap2_longcs_to_sequences(asm, work + "/seq", query=True, verbose=False)
lines = open(result["outputs"]["output"]).read().splitlines()
print(result["record_count"], lines[0], len(lines[1]))
```

- **Example output:**

```text
2 >ctgA_0_2988 2988
```

## 9) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_paf_to_gff3()`

- **Purpose:** Convert miniprot protein-to-genome PAF to GFF3 gene models.
- **Method:** Run `paftools.js paf2gff`; `amino_acid` is `-a`. The PAF needs the `cg` and `AS` tags.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_paf_to_gff3(paf: 'str | Path', output_directory: 'str | Path', *, amino_acid: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A miniprot PAF file and a new output directory.
- **Output:** `outputs` (`output` = `models.gff3`), `rows` and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_paf_to_gff3

work = tempfile.mkdtemp()
paf = work + "/miniprot.paf"
open(paf, "w").write("prot1\t100\t0\t100\t+\tchrA\t3000\t300\t600\t290\t300\t60\tAS:i:100\tcg:Z:300M\n")
result = minimap2_paf_to_gff3(paf, work + "/gff", verbose=False)
for row in result["rows"]:
    print(row[2], row[3], row[4])
```

- **Example output:**

```text
transcript 301 600
CDS 301 600
```

## 10) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_convert_simulated_reads()`

- **Purpose:** Convert reads simulated with mason2, PBSIM or Badread to FASTQ/FASTA whose names carry the true coordinates.
- **Method:** `simulator='mason2'` runs `paftools.js mason2fq` (one SAM file, FASTQ output), `'pbsim'` runs `pbsim2fq` (one or more MAF files plus `reference_fai`) and `'badread'` runs `badread2fa` (one FASTQ plus `reference_fai`). Read names become `name!chrom!start!end!strand`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_convert_simulated_reads(simulator: "Literal['mason2', 'pbsim', 'badread']", input_paths: 'Sequence[str | Path]', output_directory: 'str | Path', *, reference_fai: 'str | Path | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `simulator`, the simulator output file(s), a new output directory and, for `pbsim` and `badread`, the FASTA index of the reference.
- **Output:** `outputs` (`output` = `reads.fastq` or `reads.fasta`), `record_count` and `simulator`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_convert_simulated_reads

work = tempfile.mkdtemp()
fai = work + "/ref.fa.fai"
open(fai, "w").write("chrA\t3000\t6\t70\t71\nchrB\t2000\t3100\t70\t71\n")
mason = work + "/mason.sam"
open(mason, "w").write("simulated.r1\t0\tchrB\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXE:i:0\tXS:i:0\tXI:i:0\n")
result = minimap2_convert_simulated_reads("mason2", [mason], work + "/mason", verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
maf = work + "/pbsim.maf"
open(maf, "w").write("a\ns ref 100 8 + 5000 ACGTACGT\ns S1_1 0 8 - 8 ACGT-ACGT\n")
result = minimap2_convert_simulated_reads("pbsim", [maf], work + "/pbsim", reference_fai=fai, verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
bad = work + "/badread.fq"
open(bad, "w").write("@id1 chrA,+strand,100-108 length=8 read_identity=91.5%\nACGTACGT\n+\nIIIIIIII\n")
result = minimap2_convert_simulated_reads("badread", [bad], work + "/badread", reference_fai=fai, verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[0])
```

- **Example output:**

```text
@r1!chrB!100!110!+ 0:0:0
>S1_1!chrA!100!108!-
>S1!chrA!100!108!+	ri:f:91.5
```

## 11) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_simulated_names_to_bed()`

- **Purpose:** Convert simulated read names to BED truth coordinates.
- **Method:** Run `paftools.js sim2bed` on a file with one read name `name!chrom!start!end!strand` per line (the mason2 paired-end form is also accepted).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_simulated_names_to_bed(names_file: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A text file of read names and a new output directory.
- **Output:** `outputs` (`output` = `truth.bed`), `records` (chrom, start, end, name, score, strand) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_simulated_names_to_bed

work = tempfile.mkdtemp()
names = work + "/names.txt"
open(names, "w").write("r1!chrB!100!250!+\nr2!chrA!900!1050!-\n")
result = minimap2_simulated_names_to_bed(names, work + "/bed", verbose=False)
for row in result["records"]:
    print(row["chrom"], row["start"], row["end"], row["strand"])
```

- **Example output:**

```text
chrB 100 250 +
chrA 900 1050 -
```

## 12) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_alignment_stats()`

- **Purpose:** Summary statistics of PAF or SAM alignments.
- **Method:** Run `paftools.js stat`: mapped sequences, primary and secondary alignments, bases, substitutions and indels by size. `per_alignment_cs` is `-c` and `gap_output_length` is `-l`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_alignment_stats(alignments: 'str | Path', output_directory: 'str | Path', *, per_alignment_cs: 'bool' = False, gap_output_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A PAF (with the `cs` tag) or SAM file and a new output directory.
- **Output:** `outputs` (`output` = `stats.txt`) and `metrics`, a dictionary such as `substitutions`, `insertions_in_0_50`, `deletions_in_0_50`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_alignment_stats

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_alignment_stats(asm, work + "/stat", verbose=False)
metrics = result["metrics"]
print(metrics["mapped_bases"], metrics["substitutions"], metrics["insertions_in_0_50"], metrics["deletions_in_0_50"])
```

- **Example output:**

```text
4998 4 1 1
```

## 13) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_stats()`

- **Purpose:** Assembly contiguity statistics from assembly-to-reference alignments.
- **Method:** Run `paftools.js asmstat` with the FASTA index of the reference and one PAF per assembly. `min_query_length` is `-q`, `min_block_length` `-l` (default 10000) and `max_divergence` `-d` (0.01).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_stats(fai: 'str | Path', alignments: 'Sequence[str | Path]', output_directory: 'str | Path', *, min_query_length: 'int | None' = None, min_block_length: 'int | None' = None, max_divergence: 'float | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A `.fai` file, a list of PAF files and a new output directory.
- **Output:** `outputs` (`output` = `assembly_stats.tsv`) and `metrics`, mapping each metric (`Length`, `NG50`, `NGA50`, ...) to one value per assembly.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_assembly_stats

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
fai = work + "/ref.fa.fai"
open(fai, "w").write("chrA\t3000\nchrB\t2000\n")
result = minimap2_assembly_stats(fai, [asm], work + "/asmstat", min_block_length=1000, verbose=False)
print(result["metrics"]["Length"], result["metrics"]["NG50"], result["metrics"]["Rcov"])
```

- **Example output:**

```text
[4998] [2988] ['100.00%']
```

## 14) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_gene_completeness()`

- **Purpose:** Gene completeness of assemblies from spliced alignments of the same transcripts.
- **Method:** Run `paftools.js asmgene` with the transcripts aligned to the reference and to each assembly (`minimap2 -cx splice`). `min_identity` is `-i` (0.99), `min_coverage` `-c` (0.99), `autosomes_only` `-a` and `show_incomplete` `-e`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_assembly_gene_completeness(reference_paf: 'str | Path', assembly_pafs: 'Sequence[str | Path]', output_directory: 'str | Path', *, min_identity: 'float | None' = None, min_coverage: 'float | None' = None, autosomes_only: 'bool' = False, show_incomplete: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The reference PAF, a list of assembly PAFs and a new output directory.
- **Output:** `outputs` (`output` = `completeness.tsv`) and `metrics`, mapping `full_sgl`, `full_dup`, `frag`, `part50+`, ... to one count per input.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_assembly_gene_completeness

work = tempfile.mkdtemp()
splice = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_assembly_gene_completeness(splice, [splice], work + "/asmgene", verbose=False)
print(result["metrics"]["full_sgl"], result["metrics"]["frag"])
```

- **Example output:**

```text
[2, 2] [0, 0]
```

## 15) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_detect_misjoins()`

- **Purpose:** Detect large-scale misjoins in assembly-to-reference alignments.
- **Method:** Run `paftools.js misjoin`. `centromere_bed` is `-c`, `centromere_overlap_ratio` `-r`, `min_block_length` `-l` and `max_gap` `-g` (numbers or `1m`-style sizes), `show_non_centromeric` `-e`, `show_long_blocks` `-p`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_detect_misjoins(paf: 'str | Path', output_directory: 'str | Path', *, centromere_bed: 'str | Path | None' = None, centromere_overlap_ratio: 'float | None' = None, min_block_length: 'int | str | None' = None, max_gap: 'int | str | None' = None, show_non_centromeric: 'bool' = False, show_long_blocks: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** An assembly-to-reference PAF and a new output directory.
- **Output:** `outputs` (`output` = `misjoins.txt`) and `metrics` (the counts of inter-chromosomal misjoins, gaps and candidate inversions).
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_detect_misjoins

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_detect_misjoins(asm, work + "/misjoin", verbose=False)
print(result["metrics"])
```

- **Example output:**

```text
{'inter_chromosomal_misjoins': [0, 0], 'intra_chromosomal_gaps': [0, 0], 'candidate_inversions_in_the_middle': [0, 0], 'candidate_inversions_at_contig_ends': [0, 0]}
```

## 16) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_liftover_coordinates()`

- **Purpose:** Lift BED intervals between sequences through PAF alignments.
- **Method:** Run `paftools.js liftover` (a simplistic liftOver). The PAF needs the `cg` tag and `query_bed` names intervals on the query sequences. `min_mapping_quality` is `-q`, `min_alignment_length` `-l` (default 50000), `max_divergence` `-d` and `merge_intervals` `-m`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_liftover_coordinates(paf: 'str | Path', query_bed: 'str | Path', output_directory: 'str | Path', *, min_mapping_quality: 'int | None' = None, min_alignment_length: 'int | None' = None, max_divergence: 'float | None' = None, merge_intervals: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A PAF, a BED file and a new output directory.
- **Output:** `outputs` (`output` = `lifted.bed`), `rows` and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_liftover_coordinates

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
bed = work + "/query.bed"
open(bed, "w").write("ctgA\t100\t200\nctgA\t2500\t2600\n")
result = minimap2_liftover_coordinates(asm, bed, work + "/lift", min_alignment_length=1000, verbose=False)
for row in result["rows"]:
    print(row[:3])
```

- **Example output:**

```text
['chrA', '100', '200']
['chrA', '2512', '2612']
```

## 17) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_call_assembly_variants()`

- **Purpose:** Call variants from assembly-to-reference alignments that carry the cs tag.
- **Method:** Sort the PAF by target name and start as the tool requires (`sort -k6,6 -k8,8n`, kept as `sorted.paf`) and run `paftools.js call`. With `reference_fasta` (`-f`) the calls are written as VCF (`sample_name` is `-s`); otherwise as the tool's text. `min_coverage_length` is `-l` (10000), `min_call_length` `-L` (50000), `min_mapping_quality` `-q` and `gap_threshold` `-g`. The summary lines are parsed into `summary`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_call_assembly_variants(paf: 'str | Path', output_directory: 'str | Path', *, reference_fasta: 'str | Path | None' = None, sample_name: 'str | None' = None, min_coverage_length: 'int | None' = None, min_call_length: 'int | None' = None, min_mapping_quality: 'int | None' = None, gap_threshold: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A PAF made with `--cs`, optionally the reference FASTA, and a new output directory.
- **Output:** `outputs` (`output` = `variants.vcf` or `variants.txt`, `sorted_paf`) and `summary`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_call_assembly_variants

work = tempfile.mkdtemp()
asm = minimap2_map_reads("ref.fa", "asm.fa", work + "/asm", preset="asm5", base_level_alignment=True, cs_tag="short", verbose=False)["outputs"]["alignments"]
result = minimap2_call_assembly_variants(asm, work + "/call", reference_fasta="ref.fa", min_coverage_length=1000, min_call_length=1000, verbose=False)
records = [line.split("\t") for line in open(result["outputs"]["output"]) if not line.startswith("#")]
print(result["summary"]["substitutions"], len(records))
print([row[:5] for row in records[:2]])
```

- **Example output:**

```text
4 6
[['chrA', '501', '.', 'C', 'A'], ['chrA', '1101', '.', 'C', 'A']]
```

## 18) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_pair_haplotype_vcf()`

- **Purpose:** Combine the VCF records of two haplotypes into one phased diploid VCF.
- **Method:** Run `paftools.js vcfpair`. `male` is `-m` and requires `genome_version='37'` (`-g`) for the pseudo-autosomal regions of chrX; `sample_name` is `-s` (default `syndip`).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_pair_haplotype_vcf(pair_vcf: 'str | Path', output_directory: 'str | Path', *, male: 'bool' = False, genome_version: "Literal['37'] | None" = None, sample_name: 'str | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A VCF with one sample column per haplotype and a new output directory.
- **Output:** `outputs` (`output` = `diploid.vcf`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_pair_haplotype_vcf

work = tempfile.mkdtemp()
pair = work + "/pair.vcf"
open(pair, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thap1\thap2\nchr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,5\t0/0:5,0\n")
result = minimap2_pair_haplotype_vcf(pair, work + "/diploid", sample_name="demo", verbose=False)
print(open(result["outputs"]["output"]).read().splitlines()[-1])
```

- **Example output:**

```text
chr1	100	.	A	G	30	.	.	GT:AD	1|0:5,5
```

## 19) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_vcf_stats()`

- **Purpose:** Summary statistics of a VCF: substitutions, transitions/transversions and indels by size.
- **Method:** Run `paftools.js vcfstat`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_vcf_stats(vcf: 'str | Path', output_directory: 'str | Path', *, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A VCF file and a new output directory.
- **Output:** `outputs` (`output` = `vcf_stats.txt`) and `metrics` (`substitutions`, `ts_tv`, `insertions`, `deletions`, ...).
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_vcf_stats

work = tempfile.mkdtemp()
vcf = work + "/toy.vcf"
open(vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\t.\t.\nchr1\t200\t.\tC\tA\t.\t.\t.\nchr1\t300\t.\tT\tTAAA\t.\t.\t.\nchr1\t400\t.\tGCCC\tG\t.\t.\t.\n")
result = minimap2_vcf_stats(vcf, work + "/stats", verbose=False)
metrics = result["metrics"]
print(metrics["substitutions"], metrics["ts_tv"], metrics["insertions"], metrics["deletions"])
```

- **Example output:**

```text
2 1.0 1 1
```

## 20) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_select_vcf_by_length()`

- **Purpose:** Keep the VCF records whose variant length lies in a range.
- **Method:** Run `paftools.js vcfsel`; `min_length` is `-l` and `max_length` is `-L`. The header is kept.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_select_vcf_by_length(vcf: 'str | Path', output_directory: 'str | Path', *, min_length: 'int | None' = None, max_length: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A VCF file and a new output directory.
- **Output:** `outputs` (`output` = `selected.vcf`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_select_vcf_by_length

work = tempfile.mkdtemp()
vcf = work + "/toy.vcf"
open(vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\t.\t.\nchr1\t300\t.\tT\tTAAA\t.\t.\t.\nchr1\t400\t.\tGCCC\tG\t.\t.\t.\n")
result = minimap2_select_vcf_by_length(vcf, work + "/selected", min_length=2, verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
2
```

## 21) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_sv_callsets()`

- **Purpose:** Compare two structural-variant callsets in VCF.
- **Method:** Run `paftools.js sveval base.vcf call.vcf`. `confident_regions_bed` is `-r`, `min_discard_length` `-f` (30), `min_sv_length` `-i` (50), `max_sv_length` `-x` (100000), `window_size` `-w` (500), `max_allele_difference` `-d` (0.5), `show_errors` `-e` and `show_matches` `-p`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_sv_callsets(base_vcf: 'str | Path', call_vcf: 'str | Path', output_directory: 'str | Path', *, confident_regions_bed: 'str | Path | None' = None, min_discard_length: 'int | None' = None, min_sv_length: 'int | None' = None, max_sv_length: 'int | None' = None, window_size: 'int | None' = None, max_allele_difference: 'float | None' = None, show_errors: 'bool' = False, show_matches: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The base and the called VCF and a new output directory.
- **Output:** `outputs` (`output` = `sveval.txt`), `rows` (the tool's report lines) and `row_count`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_compare_sv_callsets

work = tempfile.mkdtemp()
header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
base, call = work + "/base.vcf", work + "/call.vcf"
open(base, "w").write(header + "chr1\t1000\t.\tA\tA" + "G" * 99 + "\t.\t.\t.\nchr1\t5000\t.\tA" + "C" * 99 + "\tA\t.\t.\t.\n")
open(call, "w").write(header + "chr1\t1010\t.\tA\tA" + "G" * 99 + "\t.\t.\t.\n")
result = minimap2_compare_sv_callsets(base, call, work + "/sveval", verbose=False)
print(result["row_count"])
print(open(result["outputs"]["output"]).read().splitlines()[0][:40])
```

- **Example output:**

```text
3
SN	2	1	0.500000
```

## 22) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_mapping()`

- **Purpose:** Evaluate mapping accuracy against the truth encoded in simulated read names.
- **Method:** Run `paftools.js mapeval`. Read names must be `name!chrom!start!end!strand` (see `minimap2_convert_simulated_reads`) or the mason2 paired-end form. `min_overlap_ratio` is `-r` (0.1), `print_wrong_mapq` `-Q` and `mode` `-m` (0 longest alignment, 1 first, 2 all primary).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_mapping(alignments: 'str | Path', output_directory: 'str | Path', *, min_overlap_ratio: 'float | None' = None, print_wrong_mapq: 'int | None' = None, mode: 'Literal[0, 1, 2] | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A PAF or SAM file of the simulated reads and a new output directory.
- **Output:** `outputs` (`output` = `mapeval.tsv`), `records` (per mapping-quality threshold: `mapq`, `reads`, `wrong`, `cumulative_error_rate`, `cumulative_reads`) and `unmapped`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_mapping

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/sim.fq"
open(reads, "w").write("".join(f"@r{i}!chrB!{s}!{s + 100}!+\n{chrB[s:s + 100]}\n+\n{'I' * 100}\n" for i, s in enumerate(range(0, 1800, 150))))
aln = minimap2_map_reads("ref.fa", reads, work + "/aln", preset="sr", base_level_alignment=True, verbose=False)["outputs"]["alignments"]
result = minimap2_evaluate_mapping(aln, work + "/mapeval", verbose=False)
print(result["records"][-1], result["unmapped"])
```

- **Example output:**

```text
{'mapq': 60, 'reads': 12, 'wrong': 0, 'cumulative_error_rate': 0.0, 'cumulative_reads': 12} None
```

## 23) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_paf()`

- **Purpose:** Compare two PAF alignment files base by base.
- **Method:** Run `paftools.js pafcmp base.paf test.paf`; `min_mapping_quality` is `-q` (10).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_compare_paf(base_paf: 'str | Path', test_paf: 'str | Path', output_directory: 'str | Path', *, min_mapping_quality: 'int | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Two PAF files and a new output directory.
- **Output:** `outputs` (`output` = `pafcmp.txt`) and `metrics` counting the base alignments that are correct, wrong, missing or additional in the test file. Only queries of at least 5000 bases with a single primary alignment are compared.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_compare_paf

work = tempfile.mkdtemp()
base, test = work + "/base.paf", work + "/test.paf"


def line(name, start):
    return f"{name}\t6000\t0\t6000\t+\tchr1\t1000000\t{start}\t{start + 6000}\t5900\t6000\t60\n"


open(base, "w").write(line("q1", 100) + line("q2", 20000) + line("q3", 40000))
open(test, "w").write(line("q1", 150) + line("q3", 90000))
result = minimap2_compare_paf(base, test, work + "/cmp", verbose=False)
print(result["metrics"])
```

- **Example output:**

```text
{'base_alignments_with_mapq_10': 3, 'base_alignments_correctly_mapped_by_test': 1, 'wrong_test_alignment': 1, 'base_alignments_missing': 1, 'additional_test_alignments_with_mapq_10': 0}
```

## 24) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_junctions()`

- **Purpose:** Evaluate the splice junctions of spliced alignments against a gene annotation.
- **Method:** Run `paftools.js junceval annotation.gtf alignments`. `input_format` is `sam` (default), `miniprot_paf` (`-a`) or `bed` (`-b`); `tolerance` is `-l`, `standard_chromosomes_only` `-c`, `first_alignment_only` `-1`, `print_overlapping` `-p` and `print_errors` `-e`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_junctions(annotation_gtf: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['sam', 'miniprot_paf', 'bed']" = 'sam', tolerance: 'int | None' = None, standard_chromosomes_only: 'bool' = False, first_alignment_only: 'bool' = False, print_overlapping: 'bool' = False, print_errors: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A GTF file, spliced alignments and a new output directory.
- **Output:** `outputs` (`output` = `junceval.txt`), `metrics` (`predicted_introns`, `correct_introns`, `correct_introns_percent`, ...), `report_text` and `rows`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_junctions

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_evaluate_junctions(gtf, splice_sam, work + "/junc", verbose=False)
print(result["metrics"])
```

- **Example output:**

```text
{'unmapped_reads': 0, 'mapped_reads': 2, 'primary_alignments': 2, 'singletons': 0, 'predicted_introns': 4, 'non_overlapping_introns': 0, 'correct_introns': 4, 'correct_introns_percent': 100.0}
```

## 25) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_exons()`

- **Purpose:** Evaluate the exons of spliced alignments against a gene annotation.
- **Method:** Run `paftools.js exoneval annotation.gtf alignments`. Options as for the junction evaluation plus `coding_only` (`-d`), `skip_first_exon` (`-f`), `skip_terminal_exons` (`-t`) and `base_level` (`-s`, base-level sensitivity and specificity).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_exons(annotation_gtf: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['sam', 'miniprot_paf', 'bed']" = 'sam', tolerance: 'int | None' = None, coding_only: 'bool' = False, standard_chromosomes_only: 'bool' = False, first_alignment_only: 'bool' = False, skip_first_exon: 'bool' = False, skip_terminal_exons: 'bool' = False, base_level: 'bool' = False, print_overlapping: 'bool' = False, print_errors: 'bool' = False, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A GTF file, spliced alignments and a new output directory.
- **Output:** `outputs` (`output` = `exoneval.txt`), `metrics` (`predicted_exons`, `correct_exons`, ...), `report_text` and `rows`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_exons

work = tempfile.mkdtemp()
splice_sam = minimap2_map_reads("ref.fa", "cdna.fa", work + "/splice", preset="splice", output_format="sam", verbose=False)["outputs"]["alignments"]
gtf = work + "/genes.gtf"
open(gtf, "w").write("".join(f'chrA\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for s, e in [(300, 450), (800, 950), (1300, 1450)]))
result = minimap2_evaluate_exons(gtf, splice_sam, work + "/exon", verbose=False)
print(result["metrics"])
```

- **Example output:**

```text
{'unmapped_reads': 0, 'mapped_reads': 2, 'primary_alignments': 2, 'predicted_exons': 6, 'non_overlapping_exons': 0, 'correct_exons': 6, 'correct_exons_percent': 100.0}
```

## 26) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_overlaps()`

- **Purpose:** Evaluate read-overlap sensitivity against the mapping of the reads to a reference.
- **Method:** Sort the read-to-reference PAF by target name and start (`sort -k6,6 -k8,8n`, kept as `sorted_to_ref.paf`) and run `paftools.js ov-eval`. `min_overlap_length` is `-l` (2000), `min_mapping_quality` `-q` (10) and `min_mapped_fraction` `-f` (0.95).
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_evaluate_overlaps(read_to_reference_paf: 'str | Path', overlap_paf: 'str | Path', output_directory: 'str | Path', *, min_overlap_length: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_mapped_fraction: 'float | None' = None, paftools_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The PAF of the reads mapped to the reference, the all-vs-all overlap PAF (`ava-ont`/`ava-pb`) and a new output directory.
- **Output:** `outputs` (`output` = `ov_eval.txt`, `sorted_to_ref_paf`), `metrics` (`inferred_overlaps`, `missed_overlaps`, `sensitivity_percent`) and `report_text`.
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_evaluate_overlaps

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/long.fa"
open(reads, "w").write("".join(f">l{i}\n{chrB[s:s + 800]}\n" for i, s in enumerate(range(0, 1000, 150))))
to_ref = minimap2_map_reads("ref.fa", reads, work + "/ref", preset="map-ont", verbose=False)["outputs"]["alignments"]
overlaps = minimap2_map_reads(reads, reads, work + "/ava", preset="ava-ont", verbose=False)["outputs"]["alignments"]
result = minimap2_evaluate_overlaps(to_ref, overlaps, work + "/ov", min_overlap_length=200, verbose=False)
print(result["report_text"])
```

- **Example output:**

```text
15 overlaps inferred from the reference mapping
0 missed by the read overlapper
100.00% sensitivity
```

## 27) `dnakit.alignment.minimap2_p09_cli_functions.minimap2_cluster_overlap_reads()`

- **Purpose:** Cluster reads from an all-vs-all overlap PAF into read groups.
- **Method:** Run `pafcluster.js` (in the `misc/` folder of the minimap2 source tree) through `k8`. `min_coverage` is `-c` (0.9), `max_divergence` `-d` (0.015) and `max_difference` `-e` (20000). The script is found on `PATH` or by `pafcluster_path`.
- **API:** `dnakit.alignment.minimap2_p09_cli_functions.minimap2_cluster_overlap_reads(overlap_paf: 'str | Path', output_directory: 'str | Path', *, min_coverage: 'float | None' = None, max_divergence: 'float | None' = None, max_difference: 'float | None' = None, pafcluster_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** An all-vs-all overlap PAF and a new output directory.
- **Output:** `outputs` (`output` = `clusters.txt`) and `rows` in the tool's record format (`SD` seed, `CL` cluster length, `//` terminated blocks).
- **Example code:**

```python
import tempfile
from dnakit.alignment.minimap2_p09_cli_functions import minimap2_map_reads, minimap2_cluster_overlap_reads

work = tempfile.mkdtemp()
chrB = "".join(open("ref.fa").read().split(">chrB")[1].split())
reads = work + "/long.fa"
open(reads, "w").write("".join(f">l{i}\n{chrB[s:s + 1000]}\n" for i, s in enumerate(range(0, 800, 200))))
overlaps = minimap2_map_reads(reads, reads, work + "/ava", preset="ava-ont", verbose=False)["outputs"]["alignments"]
result = minimap2_cluster_overlap_reads(overlaps, work + "/clusters", verbose=False)
print(result["row_count"], result["rows"][:3])
```

- **Example output:**

```text
12 [['SD', 'l0', '1'], ['CL', 'l0', '1000'], ['//']]
```
