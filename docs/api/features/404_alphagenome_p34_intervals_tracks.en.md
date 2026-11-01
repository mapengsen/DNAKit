# AlphaGenome intervals, junctions and tracks

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: Interval/track operations are local. Hosted output/scorer metadata need the caller credentials and network access; ordinary tests do not call them.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_adjust_interval`

**Purpose:** Pad or shift the boundaries of a genomic interval (alphagenome ``Interval`` methods).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_adjust_interval(interval: str | dict[str, Any], *, operation: str='pad', start: int=0, end: int=0, use_strand: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Pad or shift the boundaries of a genomic interval (alphagenome ``Interval`` methods).

``operation`` is ``pad`` (``Interval.pad`` / ``pad_inplace``: ``start`` and ``end`` bases are
ADDED on each side, negative values trim) or ``boundary_shift`` (``Interval.boundary_shift``:
the start boundary moves by ``start`` bases and the end boundary by ``end`` bases, so a
positive ``start`` moves the start to the right). With ``use_strand`` the two sides swap for a
negative-strand interval, so ``start`` always means the 5' side.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_interval_center`

**Purpose:** Centre coordinate of an interval (alphagenome ``Interval.center``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_interval_center(interval: str | dict[str, Any], *, use_strand: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Centre coordinate of an interval (alphagenome ``Interval.center``).

With ``use_strand`` the rounding of an even width follows the strand, as alphagenome does.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_interval_within_reference`

**Purpose:** Test whether an interval lies inside ``[0, reference_length]`` (``Interval.within_reference``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_interval_within_reference(interval: str | dict[str, Any], *, reference_length: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Test whether an interval lies inside ``[0, reference_length]`` (``Interval.within_reference``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_interval_strand`

**Purpose:** Flip or remove the strand of an interval.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_interval_strand(interval: str | dict[str, Any], *, mode: str='swap', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Flip or remove the strand of an interval.

``mode`` is ``swap`` (``Interval.swap_strand``: ``+`` becomes ``-`` and the reverse) or
``drop`` (``Interval.as_unstranded``: strand becomes ``.``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_interval_relation`

**Purpose:** Test whether an interval overlaps or contains another one (same chromosome required).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_interval_relation(interval: str | dict[str, Any], other: str | dict[str, Any], *, relation: str='overlaps', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Test whether an interval overlaps or contains another one (same chromosome required).

``relation`` is ``overlaps`` (``Interval.overlaps``: they share at least one base) or
``contains`` (``Interval.contains``: ``interval`` completely covers ``other``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_interval_coverage`

**Purpose:** Coverage of an interval by other intervals, as a count track or a boolean mask.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_interval_coverage(interval: str | dict[str, Any], intervals: list[Any], *, bin_size: int=1, stranded: bool=False, mask: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Coverage of an interval by other intervals, as a count track or a boolean mask.

``stranded=False, mask=False`` is ``Interval.coverage`` (number of intervals per bin);
``mask=True`` is ``Interval.binary_mask`` (coverage above zero); ``stranded=True`` is
``Interval.coverage_stranded`` / ``binary_mask_stranded``: two columns, same-strand and
opposite-strand coverage relative to ``interval``. ``bin_size`` must divide the interval width;
counts are summed within a bin.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_overlap_ranges`

**Purpose:** Relative sub-ranges of an interval that other intervals overlap (``Interval.overlap_ranges``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_overlap_ranges(interval: str | dict[str, Any], intervals: list[Any], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Relative sub-ranges of an interval that other intervals overlap (``Interval.overlap_ranges``).

Each range is ``[start, end]`` relative to the start of ``interval``, clipped to it.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_variant_overlap`

**Purpose:** Test whether the reference or alternate allele of a variant overlaps an interval.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_variant_overlap(variant: str | dict[str, Any], interval: str | dict[str, Any], *, allele: str='reference', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Test whether the reference or alternate allele of a variant overlaps an interval.

``allele`` is ``reference`` (``Variant.reference_overlaps``) or ``alternate``
(``Variant.alternate_overlaps``: the span the ALT bases would occupy from the variant start).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_split_variant`

**Purpose:** Split a variant into an upstream and a downstream part at a 0-based anchor (``Variant.split``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_split_variant(variant: str | dict[str, Any], *, anchor: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Split a variant into an upstream and a downstream part at a 0-based anchor (``Variant.split``).

If the anchor lies outside the reference span the whole variant is on one side and the other
part is None.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_junction_regions`

**Purpose:** Regions around the splice sites of a junction.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_junction_regions(junction: str | dict[str, Any], *, region: str='dinucleotide', overhang: list[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Regions around the splice sites of a junction.

``region`` is ``dinucleotide`` (``Junction.dinucleotide_region``: the first and last two bases
of the junction), ``acceptor`` (``Junction.acceptor_region``) or ``donor``
(``Junction.donor_region``); the last two are the site padded by ``overhang`` = [upstream,
downstream] bases (default [250, 250]). The junction must be stranded (``+`` or ``-``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_transcript_cds_offset`

**Purpose:** Offset of a 0-based genome position within the coding sequence (``Transcript.offset_in_cds``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_transcript_cds_offset(transcript: dict[str, Any], *, position: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Offset of a 0-based genome position within the coding sequence (``Transcript.offset_in_cds``).

The offset counts coding bases (CDS plus stop codon) from the start of the CDS along the
transcript direction, so it accounts for strand and introns. None when the position is not in
the CDS.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_fix_transcript_truncation`

**Purpose:** Repair the CDS frame of a transcript whose coding sequence is cut (``Transcript.fix_truncation``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_fix_transcript_truncation(transcript: dict[str, Any], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Repair the CDS frame of a transcript whose coding sequence is cut (``Transcript.fix_truncation``).

Each CDS and stop-codon entry needs a frame: ``[start, end, frame]``. The result has the CDS
trimmed to whole codons and new start and stop codon intervals; a transcript with 6 or fewer
coding bases is returned unchanged.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_transcripts_from_gtf`

**Purpose:** Build transcript models from the rows of a GENCODE GTF (``Transcript.from_gtf_df``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_transcripts_from_gtf(gtf_path: str, *, transcript_ids: list[str], fix_truncation: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build transcript models from the rows of a GENCODE GTF (``Transcript.from_gtf_df``).

``gtf_path`` is a GTF or GTF.gz file; ``transcript_ids`` lists the transcripts to build (their
exon, CDS, start_codon, stop_codon and Selenocysteine rows are used). With ``fix_truncation`` the
CDS frame is repaired as ``alphagenome_fix_transcript_truncation`` does. Coordinates are
0-based half-open (GTF start minus one).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_slice_track_data`

**Purpose:** Cut track data down to a genomic region.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_slice_track_data(track: dict[str, Any], *, interval: str | dict[str, Any] | None=None, start: int | None=None, end: int | None=None, match_resolution: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Cut track data down to a genomic region.

Give either ``interval`` (``TrackData.slice_by_interval``: genome coordinates, which the track's
own interval must fully contain; ``match_resolution`` widens it to whole bins) or ``start`` and
``end`` (``TrackData.slice_by_positions``: 0-based base positions relative to the track start).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_track_bin_index`

**Purpose:** Bin index of a position relative to the track start (``TrackData.bin_index``): position // resolution.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_track_bin_index(track: dict[str, Any], *, relative_position: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Bin index of a position relative to the track start (``TrackData.bin_index``): position // resolution.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_resize_track_data`

**Purpose:** Crop or zero-pad track data to ``width`` bases around its centre.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_resize_track_data(track: dict[str, Any], *, width: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Crop or zero-pad track data to ``width`` bases around its centre.

This is ``TrackData.resize`` (the model-output ``Output.resize`` does the same on its data).
``width`` must be divisible by the resolution when growing.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_change_track_resolution`

**Purpose:** Change the resolution of track data, finer (upsample) or coarser (downsample).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_change_track_resolution(track: dict[str, Any], *, resolution: int, aggregation: str='sum', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Change the resolution of track data, finer (upsample) or coarser (downsample).

A target below the current resolution is ``TrackData.upsample`` (values repeated; with ``sum``
they are divided by the repeat so the total is kept), a target above it is
``TrackData.downsample`` (values pooled with ``sum`` or ``max``); the target must divide or be
divisible by the current resolution. The result reports which direction ran.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_reverse_complement_track_data`

**Purpose:** Reverse-complement track data (``TrackData.reverse_complement``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_reverse_complement_track_data(track: dict[str, Any], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Reverse-complement track data (``TrackData.reverse_complement``).

Positions are reversed and the strands of stranded tracks swap; the interval, if present, must
be stranded (it is flipped too).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_normalize_junction_counts`

**Purpose:** Scale splice-junction read counts so that all values sum to ``total_k`` (``JunctionData.normalize_values``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_normalize_junction_counts(junction_data: dict[str, Any], *, total_k: float=10.0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Scale splice-junction read counts so that all values sum to ``total_k`` (``JunctionData.normalize_values``).

``junction_data`` is ``{"junctions": [intervals], "values": junctions x tracks, "names": [...]}``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_junctions_in_interval`

**Purpose:** Keep the splice junctions that overlap an interval (``JunctionData.intersect_with_interval``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_junctions_in_interval(junction_data: dict[str, Any], *, interval: str | dict[str, Any], verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Keep the splice junctions that overlap an interval (``JunctionData.intersect_with_interval``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_extract_fasta_interval`

**Purpose:** Read the reference sequence of an interval from an indexed FASTA (alphagenome ``FastaExtractor``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_extract_fasta_interval(fasta_path: str, *, interval: str | dict[str, Any], verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read the reference sequence of an interval from an indexed FASTA (alphagenome ``FastaExtractor``).

The FASTA needs a ``.fai`` index. Bases outside the chromosome are padded with ``N``, the
sequence is upper-case, and a negative-strand interval returns the reverse complement.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_output_metadata`

**Purpose:** List the output tracks of the hosted AlphaGenome DNA model (``DnaClient.output_metadata``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_output_metadata(api_key: str, *, organism: str='homo_sapiens', max_rows: int=50, timeout: float | None=None, address: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List the output tracks of the hosted AlphaGenome DNA model (``DnaClient.output_metadata``).

Network call with an API key (sent only to the service). For each output type (atac, cage,
dnase, rna_seq, chip_histone, chip_tf, splice_sites, splice_site_usage, splice_junctions,
contact_maps, procap) the result gives the number of tracks, the columns and the first
``max_rows`` rows. ``organism`` is ``homo_sapiens`` or ``mus_musculus``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_scorer_metadata`

**Purpose:** List the variant scorers of the AlphaGenome Atlas service (``AtlasClient.scorer_metadata``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_scorer_metadata(api_key: str, *, max_rows: int=50, timeout: float | None=None, address: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List the variant scorers of the AlphaGenome Atlas service (``AtlasClient.scorer_metadata``).

Network call with an API key. For each scorer: its name, whether scores are signed, and the
track metadata (row count, columns, first ``max_rows`` rows).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_plot_contact_map_difference`

**Purpose:** Draw the difference between two contact maps (ALT minus REF) with a red-blue scale centred on 0.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_plot_contact_map_difference(contact_maps: dict[str, Any], output_path: str, *, vmin: float=-1.0, vmax: float=1.0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Draw the difference between two contact maps (ALT minus REF) with a red-blue scale centred on 0.

This is alphagenome's ``ContactMapsDiff`` component. ``contact_maps`` is track data whose values
are bins x bins x tracks (one square map per track) and that has an ``interval``. The scale is
clipped to ``[vmin, vmax]``. Writes a PNG and returns its path with the data.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_plot_overlaid_tracks`

**Purpose:** Draw several track sets on top of each other, such as REF and ALT predictions.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_plot_overlaid_tracks(tracks: dict[str, Any], output_path: str, *, interval: str | dict[str, Any] | None=None, shared_y_scale: bool=False, alpha: float=0.8, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Draw several track sets on top of each other, such as REF and ALT predictions.

This is alphagenome's ``OverlaidTracks`` component. ``tracks`` maps a label (``"REF"``,
``"ALT"``) to track data; the sets need the same track names and an ``interval``. ``interval``
is the plotted window (default: the interval of the first set). Writes a PNG.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.alphagenome_p34_functions.alphagenome_plot_annotations`

**Purpose:** Mark genomic intervals and/or variants across a plot window.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def alphagenome_plot_annotations(interval: str | dict[str, Any], output_path: str, *, intervals: list[Any] | None=None, variants: list[Any] | None=None, track: dict[str, Any] | None=None, labels: list[str] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Mark genomic intervals and/or variants across a plot window.

``intervals`` are drawn with alphagenome's ``IntervalAnnotation`` (shaded rectangles) and
``variants`` with ``VariantAnnotation`` (thin lines with labels). ``track`` (optional) is track
data drawn underneath; without it an empty panel is used. ``labels`` replaces the default
labels of the intervals (when given) or of the variants. Writes a PNG.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.predictions.alphagenome_p34_functions import alphagenome_interval_center
r = alphagenome_interval_center({'chromosome': 'chr1', 'start': 10, 'end': 20}, verbose=False)
print(r['center'])
```

**Example output:**

```text
15
```
