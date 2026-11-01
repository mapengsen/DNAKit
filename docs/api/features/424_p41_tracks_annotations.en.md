# P41 tracks and annotations

## Purpose

Complete selected P41 operations through dedicated upstream APIs, returning JSON data or actual file paths.

## Method

Upstream algorithms and defaults are preserved. Choose algorithm/output variants with explicit parameters. Existing outputs are refused. Optional dependencies are required only when called.

## API

| Module | Functions |
| --- | --- |
| `dnakit.annotation.prokka_p41_cli_functions` | `annotate_prokka_rnammer` |
| `dnakit.annotation.pybigwig_p41_functions` | `inspect_bigwig_file` |
| `dnakit.databases.gget_p41_functions` | `gget_download_reference`, `gget_isoform_structure_map` |
| `dnakit.epigenome.epigenome_p41_cli_functions` | `pairtools_tile_duplication_stats`, `macs3_bedgraph_cutoff_report`, `macs3_duplicate_statistics`, `bam_region_count_matrix`, `bigwig_region_score_matrix` |
| `dnakit.molecular_lab.dna_features_viewer_p41_functions` | `read_gff_feature_record`, `export_feature_map_record`, `plot_feature_map` |
| `dnakit.variants.cyvcf2_p41_functions` | `cyvcf2_ibd_bins`, `cyvcf2_informative_sites`, `cyvcf2_export_text`, `vcf_call_is_filtered` |

## Input

Use explicit local file paths, JSON-compatible arrays and 0-based half-open intervals unless an API documents a native alternative. External CLI programs must be installed. Refer to each function docstring for limits and parameter choices.

## Output

Dictionaries contain status and parsed results, or paths of files actually written. Plot interfaces also return the input sequence and features. Tests and limitations: `paper/write_paper/figure/fig2/task_coverage/coverage_gap/manifests/integration_p41.json`.

GFF parsing requires bcbio-gff. gget functions perform network operations and are excluded from ordinary live tests. Prokka/RNAmmer and bigBed autoSQL lack real end-to-end fixtures in this environment; no successful scientific validation is claimed.

### Feature record export

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from Bio import SeqIO
from dnakit.molecular_lab.dna_features_viewer_p41_functions import export_feature_map_record
with TemporaryDirectory() as folder:
    p = Path(folder) / 'dna.gb'
    export_feature_map_record('ATGAAATAG', [{'start':0,'end':9,'label':'gene'}], str(p), verbose=False)
    r = SeqIO.read(p, 'genbank')
    print(str(r.seq), len(r.features))
```

**Example output:**

```text
ATGAAATAG 1
```

### Genotype filtering

**Example code:**

```python
from dnakit.variants.cyvcf2_p41_functions import vcf_call_is_filtered
print(vcf_call_is_filtered('0/1', filters=['PASS'], verbose=False)['is_filtered'])
print(vcf_call_is_filtered('0/1', filters=['q10'], verbose=False)['is_filtered'])
```

**Example output:**

```text
False
True
```

### bigWig type inspection

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
import pyBigWig
from dnakit.annotation.pybigwig_p41_functions import inspect_bigwig_file
with TemporaryDirectory() as folder:
    p = Path(folder) / 'scores.bw'
    h = pyBigWig.open(str(p), 'w')
    h.addHeader([('chr1', 100)])
    h.addEntries(['chr1'], [0], ends=[5], values=[1.0])
    h.close()
    r = inspect_bigwig_file(str(p), verbose=False)
    print(r['is_bigwig'], r['is_bigbed'], r['autosql'])
```

**Example output:**

```text
True False None
```

### MACS3 duplicate statistics

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.epigenome.epigenome_p41_cli_functions import macs3_duplicate_statistics
with TemporaryDirectory() as folder:
    p = Path(folder) / 'reads.bed'
    p.write_text('chr1\t0\t10\ta\t0\t+\n' * 2)
    r = macs3_duplicate_statistics(p, file_format='BED', keep_duplicates=1, verbose=False)
    print(r['statistics']['reads_before'], r['statistics']['reads_after'])
```

**Example output:**

```text
2 1
```

