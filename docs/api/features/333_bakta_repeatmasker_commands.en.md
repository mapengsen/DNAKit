# Bakta and RepeatMasker: genome annotation, repeat masking and repeat-annotation tools

This group of functions in `dnakit.annotation` is a set of command-line facades over [Bakta](https://github.com/oschwengers/bakta) and [RepeatMasker](https://www.repeatmasker.org/). Bakta: `annotate_bakta_features` annotates a bacterial genome and lets the caller choose the detection steps (CDS with their functional annotation by database look-ups, AMR genes and expert proteins, pseudogenes, sORFs, tRNA, tmRNA, rRNA, ncRNA, ncRNA regions, CRISPR arrays, assembly gaps, oriC/oriV/oriT; each skipped step is a `--skip-*` switch), `export_bakta_annotation` re-writes a Bakta JSON in every format (TSV, inference TSV, hypotheticals, GenBank, EMBL, GFF3, FASTA nucleotide/protein, summary) and `plot_bakta_genome` draws the circular maps (including the COG map). RepeatMasker: `mask_repeats`, `clip_bacterial_insertion_sequences`, `process_repeat_hits`, `date_repeats`, `annotate_segmental_duplications`, `mask_repeat_proteins`, `mask_tandem_repeats` run the programs, and the `util/` scripts are covered by `convert_repeats_to_gff3`, `convert_repeats_to_bed`, `export_repeats_ucsc_tables`, `build_repeat_track_hub`, `extract_repeat_sequences`, `mask_sequences_from_annotation`, `summarize_repeat_annotation`, `recalculate_repeat_divergence`, `plot_repeat_landscape`, `plot_duplicon_map` and `combine_repeat_annotations`. Not provided: `wublastToCrossmatch.pl` (reformats output of the discontinued WU-BLAST engine) and `addRepBase.pl` (needs the licensed RepBase RepeatMasker Edition).

Bakta 1.12.2 (GPL-3.0) and RepeatMasker 4.2.6 (Open Software License) are the tested versions and are started as separate programs. Install Bakta with `conda install -c bioconda bakta` (plus its database, `bakta_db download`) and RepeatMasker with `conda install -c bioconda repeatmasker` (then `perl configure`), put the programs on PATH or pass `bakta_path`, `bakta_io_path`, `bakta_plot_path` or `repeatmasker_dir`. Every function writes into a new `output_directory`; RepeatMasker programs write next to their input, so inputs are copied there. The full `bakta` annotation needs the database and several tools and was verified only up to argument parsing; `bakta_io`, `bakta_plot` and the RepeatMasker functions (except the ones that need FamDB: `clip_bacterial_insertion_sequences`, `date_repeats`, `annotate_segmental_duplications`, `mask_repeat_proteins`) were run for real. The examples use `verbose=False`.

## Example data and examples

**Example file:** `example.out`

```text
   SW   perc perc perc  query     position in query    matching repeat       position in repeat
score   div. del. ins.  sequence  begin end   (left)   repeat   class/family begin  end    (left)  ID

  734   27.6  9.9  4.0  seq1       1483  1806  (326) + L1M5     LINE/L1      5368   5710  (414)   1
 2134   10.7  3.3  1.3  seq1       1809  2111   (21) + AluSx    SINE/Alu        3    311    (1)   2
```

**Example code:**

```python
from pathlib import Path
import tempfile
from dnakit.annotation import repeatmasker_p11_cli_functions as rm

work = Path(tempfile.mkdtemp())
gff = rm.convert_repeats_to_gff3("example.out", work / "gff", verbose=False)
print(gff["features"])
print(Path(gff["outputs"]["repeats.gff3"]).read_text().splitlines()[2])
bed = rm.convert_repeats_to_bed("example.out", work / "bed", verbose=False)
print(bed["intervals"])
```

**Example output:**

```text
2
seq1	RepeatMasker	dispersed_repeat	1483	1806	734	+	.	ID=1;Target=L1M5 5368 5710
2
```

**Example code:**

```python
from pathlib import Path
import json, random, tempfile
from dnakit.annotation import bakta_p11_cli_functions as bk

rng = random.Random(2)
seq = "".join(rng.choice("ACGT") for _ in range(3000))
cds = {"sequence": "c1", "strand": "+", "nt": seq[100:166], "aa": "MKTAYIAKQRQISFVKSHFSRQ", "aa_digest": None,
       "start_type": "ATG", "rbs_motif": None, "frame": 2, "functions": [], "genes": []}
data = {
    "genome": {"genus": "E", "species": "coli", "strain": "", "complete": False, "gram": "?", "translation_table": 11},
    "stats": {"no_sequences": 1, "size": 3000, "gc": 0.5, "n_ratio": 0.0, "n50": 3000, "n90": 3000, "coding_ratio": 0.1},
    "features": [
        {**cds, "type": "cds", "start": 101, "stop": 166, "gene": "abcD", "product": "test protein", "locus": "T_00005", "db_xrefs": ["COG:J"]},
        {"type": "tRNA", "sequence": "c1", "start": 501, "stop": 576, "strand": "+", "gene": "tRNA-Ala", "product": "tRNA-Ala(tgc)",
         "locus": "T_00015", "nt": seq[500:576], "amino_acid": "Ala", "anti_codon": "tgc", "score": 80.5, "db_xrefs": ["SO:0001272"]},
    ],
    "sequences": [{"id": "c1", "description": "", "nt": seq, "length": 3000, "complete": False, "type": "contig",
                   "topology": "linear", "name": None, "orig_id": "c1", "orig_description": ""}],
    "version": {"bakta": "1.12.2", "db": {"type": "light", "version": "6.0", "doi": "10.5281/zenodo.0"}},
    "run": {"start": "2026-10-20 10:00:00", "end": "2026-10-20 10:01:00"},
}
work = Path(tempfile.mkdtemp())
(work / "demo.json").write_text(json.dumps(data))
export = bk.export_bakta_annotation(work / "demo.json", work / "export", verbose=False)
print(export["version"], sorted(k for k in export["outputs"] if k in ("gbff", "embl", "gff3", "tsv", "faa", "ffn", "txt")))
print(export["summary"]["CDSs"], export["summary"]["tRNAs"])
plot = bk.plot_bakta_genome(work / "demo.json", work / "plot", plot_type="cog", size=4, dpi=150, verbose=False)
print(sorted(k for k in plot["outputs"] if k in ("png", "svg")))
```

**Example output:**

```text
1.12.2 ['embl', 'faa', 'ffn', 'gbff', 'gff3', 'tsv', 'txt']
1 1
['png', 'svg']
```


## 1) `dnakit.annotation.bakta_p11_cli_functions.annotate_bakta_features()`

- **Purpose:** Annotate a bacterial genome, MAG or plasmid with Bakta, choosing which feature types are detected.
- **Method:** ``features`` selects the detection steps to run (default: all, as Bakta does): ``cds`` (protein-coding genes with their functional annotation: UniRef100/UniParc/PSC/PSCC database look-ups, AMRFinderPlus resistance genes and the curated expert protein sequences matched with DIAMOND), ``pseudogenes``, ``small_orfs`` (sORFs below 30 aa), ``trna`` (tRNAscan-SE), ``tmrna`` (ARAGORN), ``rrna`` and ``ncrna`` and ``ncrna_regions`` (Infernal cmscan against Rfam), ``crispr`` (PILER-CR), ``gaps`` (runs of N in the assembly) and ``origins`` (oriC/oriV/oriT, BLASTn against the origin databases). Steps that are not listed are skipped (``--skip-*``); ``pseudogenes`` needs ``cds``. ``database`` is a Bakta database folder. Further options are Bakta's own: organism (``genus``, ``species``, ``strain``, ``plasmid``), ``complete`` replicons, ``translation_table`` (11, 4 or 25), ``gram`` (``+``, ``-`` or ``?`` for signal peptide prediction), locus naming (``locus``, ``locus_tag``, ``locus_tag_increment``), ``compliant`` INSDC output, ``replicons`` table, pre-annotated ``regions`` (GFF3 or GenBank), ``trusted_proteins`` (FASTA) and ``trusted_hmms`` for CDS annotation, ``metagenome`` mode, ``partial`` genes at sequence ends and ``skip_overlap_filters`` (``--skip-filter``). ``plot`` also draws the circular genome plots (off by default; use ``plot_bakta_genome`` for them).
- **API:** `dnakit.annotation.bakta_p11_cli_functions.annotate_bakta_features(genome: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, features: 'Sequence[str] | None' = None, prefix: 'str' = 'annotation', min_contig_length: 'int | None' = None, genus: 'str | None' = None, species: 'str | None' = None, strain: 'str | None' = None, plasmid: 'str | None' = None, complete: 'bool' = False, prodigal_training_file: 'str | Path | None' = None, translation_table: 'int | None' = None, gram: 'str | None' = None, locus: 'str | None' = None, locus_tag: 'str | None' = None, locus_tag_increment: 'int | None' = None, keep_contig_headers: 'bool' = False, compliant: 'bool' = False, replicons: 'str | Path | None' = None, regions: 'str | Path | None' = None, trusted_proteins: 'str | Path | None' = None, trusted_hmms: 'str | Path | None' = None, metagenome: 'bool' = False, partial: 'bool' = False, skip_overlap_filters: 'bool' = False, plot: 'bool' = False, threads: 'int | None' = None, bakta_path: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `genome`, `database`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``annotation_directory``, ``outputs`` (keyed by extension: ``json``, ``tsv``, ``gff3``, ``gbff``, ``embl``, ``fna``, ``ffn``, ``faa``, ``hypotheticals.tsv``, ``hypotheticals.faa``, ``inference.tsv``, ``txt``), ``features`` (the steps run), ``feature_counts`` (per type, from the JSON) and ``summary`` (the report in ``<prefix>.txt``).

## 2) `dnakit.annotation.bakta_p11_cli_functions.export_bakta_annotation()`

- **Purpose:** Write a Bakta JSON annotation in all output formats (``bakta_io``), without re-running the annotation.
- **Method:** ``annotation_json`` is the ``.json`` file written by ``annotate_bakta_features`` (or by ``bakta``). The export contains the human readable feature table (``tsv``), the feature inference table (``inference.tsv``), the hypothetical CDS table (``hypotheticals.tsv``), GenBank (``gbff``), EMBL (``embl``) and GFF3 (``gff3``) flat files, the genome sequences (``fna``), the CDS/feature nucleotide sequences (``ffn``), the translated CDS sequences (``faa``), the hypothetical CDS proteins (``hypotheticals.faa``), the circular genome plot (``png`` and ``svg``) and the genome and annotation summary report (``txt``). ``prefix`` defaults to the stem of the input file.
- **API:** `dnakit.annotation.bakta_p11_cli_functions.export_bakta_annotation(annotation_json: 'str | Path', output_directory: 'str | Path', *, prefix: 'str | None' = None, bakta_io_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `annotation_json`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``export_directory``, ``outputs`` (keyed by extension) and ``summary`` (the report).

## 3) `dnakit.annotation.bakta_p11_cli_functions.plot_bakta_genome()`

- **Purpose:** Draw circular genome maps from a Bakta JSON annotation (``bakta_plot``).
- **Method:** ``plot_type`` ``features`` (default) colours features by type; ``cog`` colours CDS by COG functional category. ``sequences`` selects the sequences to plot by one-based number or name (default all); ``config`` is a YAML plot configuration; ``label`` is the centre label (``|`` makes a line break); ``size`` is 4, 8 or 16 inches (Bakta default 8) and ``dpi`` 150, 300 or 600 (Bakta default 300). The images are ``outputs['png']`` and ``outputs['svg']``.
- **API:** `dnakit.annotation.bakta_p11_cli_functions.plot_bakta_genome(annotation_json: 'str | Path', output_directory: 'str | Path', *, plot_type: 'str' = 'features', prefix: 'str | None' = None, sequences: 'Sequence[str] | None' = None, config: 'str | Path | None' = None, label: 'str | None' = None, size: 'int | None' = None, dpi: 'int | None' = None, bakta_plot_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `annotation_json`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``plot_directory`` and ``plot_type``.

## 4) `dnakit.annotation.repeatmasker_p11_cli_functions.annotate_segmental_duplications()`

- **Purpose:** Annotate segmental duplications (duplicons) by searching DNA against the duplicon library (DupMasker).
- **Method:** DupMasker repeat-masks the sequence (human libraries), searches the masked sequence against the duplicon library, builds realignment regions for each duplicon seed, realigns them on the unmasked DNA and joins the duplicons; the annotation is ``<name>.duplicons``. An existing RepeatMasker annotation can be given as ``repeatmasker_out`` so that the masking step is not repeated. ``max_divergence`` (``-maxDiv``) filters duplicon seeds, ``max_gap_width`` (``-maxWidth``, default 300) is the largest non-seed gap that is realigned, ``force_search`` (``-forceSearch``) repeats searches whose results exist, ``alignments`` (``-align``) stores alignments, ``gff_output`` (``-gff``) adds a GFF file, ``search_engine`` is ``abblast``, ``wublast`` or ``ncbi``, ``parallel`` is ``-pa`` and ``seeds_only`` (``-dupout``) stops after the seed search.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.annotate_segmental_duplications(fasta: 'str | Path', output_directory: 'str | Path', *, repeatmasker_out: 'str | Path | None' = None, max_divergence: 'float | None' = None, max_gap_width: 'int | None' = None, force_search: 'bool' = False, alignments: 'bool' = False, gff_output: 'bool' = False, search_engine: 'str | None' = None, parallel: 'int | None' = None, seeds_only: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``duplicons`` (path of the annotation or None).

## 5) `dnakit.annotation.repeatmasker_p11_cli_functions.build_repeat_track_hub()`

- **Purpose:** Build a UCSC bigRmsk track hub (bigBed files and hub text files) from a repeat annotation (``util/rmToTrackHub.pl``).
- **Method:** ``genome`` is the UCSC genome (assembly) name, ``hub_name`` the hub folder name (no spaces), created inside ``output_directory``. The script needs the UCSC program ``bedToBigBed`` on PATH. The chromosome sizes are downloaded from UCSC (network access) unless ``chrom_sizes`` (a two-column file) is given. ``alignments`` (``.align``) adds the alignment bigBed. ``only_bigbed`` (``-only_beds``) writes only the bigBed files, not ``hub.txt``, ``genomes.txt``, ``trackDB.txt`` and the HTML stubs.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.build_repeat_track_hub(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, genome: 'str', hub_name: 'str', alignments: 'str | Path | None' = None, chrom_sizes: 'str | Path | None' = None, only_bigbed: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``hub_directory`` and ``hub_files`` (file names in the hub folder).

## 6) `dnakit.annotation.repeatmasker_p11_cli_functions.clip_bacterial_insertion_sequences()`

- **Purpose:** Clip bacterial insertion-sequence (IS) elements out of DNA, for example E. coli contamination (``RepeatMasker -is_only``).
- **Method:** RepeatMasker searches the sequences for E. coli IS elements (the IS library is part of FamDB), removes the elements that are flanked by a short duplicated target site and writes the remaining sequence to ``<name>.withoutIS``; the matches that are only reported (not clipped) go to ``<name>.alert``. Nothing else is masked or annotated.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.clip_bacterial_insertion_sequences(fasta: 'str | Path', output_directory: 'str | Path', *, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``clipped_fasta`` (path or None when no IS element was found) and ``outputs``.

## 7) `dnakit.annotation.repeatmasker_p11_cli_functions.combine_repeat_annotations()`

- **Purpose:** Merge two RepeatMasker annotation sets (``.out`` and ``.align``) into one sorted set (``util/combineRMFiles.pl``).
- **Method:** Each set is a RepeatMasker ``.out`` annotation with its ``.align`` alignment file (a run with ``-a``). The second set must come from searching sequence pre-masked with the first, so the annotations do not overlap; element IDs are renumbered and the entries are placed in sequence and position order. The result is ``<output_prefix>.out`` and ``<output_prefix>.align``.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.combine_repeat_annotations(first_annotation: 'str | Path', first_alignments: 'str | Path', second_annotation: 'str | Path', second_alignments: 'str | Path', output_directory: 'str | Path', *, output_prefix: 'str' = 'combined', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `first_annotation`, `first_alignments`, `second_annotation`, `second_alignments`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``annotations`` and ``class_counts`` of the combined ``.out``.

## 8) `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_bed()`

- **Purpose:** Convert a repeat annotation (``.out`` or ``.align``) to BED with class, length and divergence filters (``util/RM2Bed.py``).
- **Method:** The BED file(s) are named ``<prefix>_rm.bed`` (or ``<prefix>_<group>.bed`` with ``split_by``). Options: ``output_prefix`` (``-p``, default the first field of the input name), ``split_by`` (``-sp``: ``name``, ``family``, ``class`` or ``subclass``), ``sort_by`` (``-s``: ``size``, ``name``, ``family``, ``class``, ``subclass`` or ``diverge``), ``min_length`` (``-m``), ``min_hits`` (``-n``, with ``split_by``), ``max_divergence`` (``-dmax``) and ``min_divergence`` (``-dmin``), and ``overlap_resolution`` (``-o``: ``higher_score``, ``longer_element`` or ``lower_divergence``, which assigns every base to one element). The script needs the pandas package.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_bed(repeat_file: 'str | Path', output_directory: 'str | Path', *, output_prefix: 'str | None' = None, split_by: 'str | None' = None, sort_by: 'str | None' = None, min_length: 'int | None' = None, min_hits: 'int | None' = None, max_divergence: 'float | None' = None, min_divergence: 'float | None' = None, overlap_resolution: 'str | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_file`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``bed_files`` (paths) and ``intervals`` (total rows).

## 9) `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_gff3()`

- **Purpose:** Convert a RepeatMasker annotation (``.out``) to GFF3 (``util/rmOutToGFF3.pl``).
- **Method:** Every annotation becomes a ``dispersed_repeat`` feature (source ``RepeatMasker``, score = Smith-Waterman score, ``ID`` = the element ID, ``Target`` = the repeat family and its coordinates).
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_gff3(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'repeats.gff3', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``features``.

## 10) `dnakit.annotation.repeatmasker_p11_cli_functions.date_repeats()`

- **Purpose:** Label repeats as lineage-specific or ancestral relative to other species, and optionally mask the lineage-specific ones.
- **Method:** ``repeat_annotation`` is the RepeatMasker ``.out`` file of the analysed sequence from ``query_species``; ``comparison_species`` are the other (mammalian) species or clades, resolved with FamDB. DateRepeats writes a copy of the annotation with one added column per comparison species telling whether the repeat is expected to be present there (``<name>.out_<species>...``). With ``mask_species`` (one of the comparison species) it also writes the sequence file with all repeats that are in the query but absent from that species masked, for genome alignments; this needs ``sequence_file``, the FASTA that the annotation was made from. ``mask_uncertain`` (``-aggressive``) also masks the repeats whose age is unclear; ``mask_low_complexity=False`` (``-nolow``) leaves satellites, simple repeats and low-complexity DNA unmasked. FamDB must be configured in the RepeatMasker installation.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.date_repeats(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, query_species: 'str', comparison_species: 'Sequence[str]', sequence_file: 'str | Path | None' = None, mask_species: 'str | None' = None, mask_uncertain: 'bool' = False, mask_low_complexity: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 11) `dnakit.annotation.repeatmasker_p11_cli_functions.export_repeats_ucsc_tables()`

- **Purpose:** Export a repeat annotation as UCSC ``rmsk`` tables (``util/rmToUCSCTables.pl``).
- **Method:** Writes ``<name>.out.tsv`` (the annotation as the UCSC ``rmsk`` table) and ``<name>.join.tsv`` (the table joined with the element IDs) and, when ``alignments`` (the ``.align`` file, optionally gzipped) is given, ``<name>.align.tsv``; the alignment table needs ``intersectBed`` (bedtools) on PATH.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.export_repeats_ucsc_tables(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, alignments: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``tables`` (file name to row count).

## 12) `dnakit.annotation.repeatmasker_p11_cli_functions.extract_repeat_sequences()`

- **Purpose:** Extract the sequence of every annotated repeat element as FASTA (``util/rmOut2Fasta.pl``).
- **Method:** ``fasta`` is the file with all sequences named in the RepeatMasker ``.out`` annotation. Each record is named ``<sequence>:<start>-<end> ( <repeat>:<repeat start>-<repeat end> orient=<+|-> )`` and is the genomic sequence of that annotation, in the orientation of the genome.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.extract_repeat_sequences(fasta: 'str | Path', repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'repeats.fa', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``sequences``.

## 13) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeat_proteins()`

- **Purpose:** Mask repeat-derived protein-coding fragments in DNA by a translated search against repeat proteins (RepeatProteinMask).
- **Method:** The DNA is compared (BLASTX) with the RepeatPeps database of repeat-encoded proteins; matching regions are masked and annotated. Low-complexity and simple repeats are always searched first. ``search_engine`` is ``ncbi``, ``abblast`` or ``wublast``; ``p_value`` (``-pvalue``) and ``min_score`` (``-minscore``) are thresholds (none by default); ``word_size`` is ``-wordsize`` (default 3), ``query_statistics_length`` is ``-queryStatLen``, ``max_amino_acid_distance`` is ``-maxAADist`` (default 333) and ``mask_low_complexity=False`` is ``-noLowSimple``. The annotation is ``<name>.annot`` and the masked sequence ``<name>.masked``; RepeatProteinMask reads ``Libraries/RepeatPeps.lib`` of the installation.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeat_proteins(fasta: 'str | Path', output_directory: 'str | Path', *, search_engine: 'str | None' = None, p_value: 'float | None' = None, min_score: 'int | None' = None, word_size: 'int | None' = None, query_statistics_length: 'int | None' = None, max_amino_acid_distance: 'int | None' = None, mask_low_complexity: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 14) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeats()`

- **Purpose:** Screen DNA for interspersed repeats and low-complexity sequence: annotation and masked FASTA (RepeatMasker).
- **Method:** ``species`` (``-species``, a clade or species name resolved with FamDB; RepeatMasker's default is human) or ``library`` (``-lib``, a FASTA of repeat consensus sequences, ``>name#class/subclass``) selects the repeats searched; ``search_engine`` is ``-engine`` (``crossmatch``, ``wublast``, ``abblast``, ``ncbi``, ``rmblast`` or ``hmmer``). ``speed`` is ``quick`` (``-q``), ``rush`` (``-qq``) or ``slow`` (``-s``); ``parallel`` is ``-pa``; ``cutoff`` is the score cut-off for custom libraries (``-cutoff``, default 225). What is masked: ``mask_low_complexity=False`` is ``-nolow``, ``mask_interspersed=False`` is ``-noint`` (only simple repeats and low complexity), ``screen_small_rnas=False`` is ``-norna``, ``alu_only`` is ``-alu`` and ``max_divergence`` is ``-div``. Search details: ``gc_level`` (``-gc``), ``calculate_gc`` (``-gccalc``), ``max_fragment_length`` (``-frag``), ``skip_excision`` (``-nocut``). ``insertion_sequences`` controls the bacterial insertion-element check: ``report`` (default), ``clip`` (``-is_clip``, clip them out before the analysis) or ``skip`` (``-no_is``). Output: ``mask_style`` ``n`` (repeats replaced by N, default), ``x`` (``-x``), ``lowercase_all`` (``-small``) or ``soft`` (``-xsmall``, repeats lower case); ``lowercase_ambiguous`` (``-lcambig``); extra files: ``alignments`` (``-a``, ``.align``), ``alignments_in_repeat_orientation`` (``-inv``), ``polymorphic_simple_repeats`` (``-poly``), ``ace_output``, ``gff_output`` (``-gff``), ``html_output``, ``untouched_annotation`` (``-u``), ``crossmatch_output`` (``-xm``), ``no_id_column``, ``exclude_n_runs`` (``-excln``), ``uncurated_families``, ``cache_dir`` and ``famdb_data_dir``.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeats(fasta: 'str | Path', output_directory: 'str | Path', *, species: 'str | None' = None, library: 'str | Path | None' = None, search_engine: 'str | None' = None, speed: 'str | None' = None, parallel: 'int | None' = None, cutoff: 'int | None' = None, mask_low_complexity: 'bool' = True, mask_interspersed: 'bool' = True, screen_small_rnas: 'bool' = True, alu_only: 'bool' = False, max_divergence: 'float | None' = None, gc_level: 'int | None' = None, calculate_gc: 'bool' = False, max_fragment_length: 'int | None' = None, skip_excision: 'bool' = False, insertion_sequences: 'str' = 'report', mask_style: 'str' = 'n', lowercase_ambiguous: 'bool' = False, alignments: 'bool' = False, alignments_in_repeat_orientation: 'bool' = False, polymorphic_simple_repeats: 'bool' = False, ace_output: 'bool' = False, gff_output: 'bool' = False, html_output: 'bool' = False, untouched_annotation: 'bool' = False, crossmatch_output: 'bool' = False, no_id_column: 'bool' = False, exclude_n_runs: 'bool' = False, uncurated_families: 'bool' = False, cache_dir: 'str | Path | None' = None, famdb_data_dir: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``summary`` (sequences, total length, GC, masked bases and percent from ``.tbl``), ``annotations`` and ``class_counts`` (from ``.out``) and ``outputs`` (``<name>.out``, ``<name>.masked``, ``<name>.tbl``, ``<name>.cat`` and the optional files).

## 15) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_sequences_from_annotation()`

- **Purpose:** Hard- or soft-mask a FASTA file from a RepeatMasker annotation (``util/maskFile.pl``).
- **Method:** Repeats listed in the ``.out`` file are replaced by N (default) or written in lower case (``soft_mask``, the script's ``-softmask``). ``min_divergence`` and ``max_divergence`` (``-minDivergence`` / ``-maxDivergence``, percent) mask only the repeats within that divergence range. The result is ``<name>.masked`` in the output folder.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_sequences_from_annotation(fasta: 'str | Path', repeat_annotation: 'str | Path', output_directory: 'str | Path', *, soft_mask: 'bool' = False, min_divergence: 'float | None' = None, max_divergence: 'float | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``masked_fasta``, ``masked_bases`` (N or lower-case bases in the output).

## 16) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_tandem_repeats()`

- **Purpose:** Mask simple tandem repeats in DNA with Tandem Repeats Finder (RepeatMasker's ``util/trfMask``).
- **Method:** TRF (configured in the RepeatMasker installation) finds the tandem repeats; the masked sequence is written to ``<name>.masked`` (repeats replaced by N, or lower case with ``soft_mask=True`` which is ``-xsmall``), with the annotation in ``<name>.out`` and the alignments in ``<name>.align``.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_tandem_repeats(fasta: 'str | Path', output_directory: 'str | Path', *, soft_mask: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``summary`` (the line trfMask prints: sequence count, total length, GC level and masked bases).

## 17) `dnakit.annotation.repeatmasker_p11_cli_functions.plot_duplicon_map()`

- **Purpose:** Draw the DupMasker duplicon map of every sequence as SVG (``util/dupliconToSVG.pl``).
- **Method:** ``duplicons_file`` is the ``.duplicons`` annotation written by ``annotate_segmental_duplications`` (the script works with the human duplicon library). One SVG file per sequence is written, named ``<file>.<n>.svg``.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.plot_duplicon_map(duplicons_file: 'str | Path', output_directory: 'str | Path', *, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `duplicons_file`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``svg_files`` (paths).

## 18) `dnakit.annotation.repeatmasker_p11_cli_functions.plot_repeat_landscape()`

- **Purpose:** Plot the repeat landscape: genome fraction occupied by each repeat class at each divergence level (``util/createRepeatLandscape.pl``).
- **Method:** ``divergence_summary`` is the ``.divsum`` file of ``recalculate_repeat_divergence``. The genome size for the percentages is ``genome_size`` (``-g``, bp) or taken from ``two_bit_file`` (``-twoBit``, needs the UCSC program ``twoBitInfo``). The plot is a stacked-area Google Charts page written to ``output_name`` (HTML; the page loads the charting library from google.com when opened) or only its JavaScript with ``javascript_only`` (``-j``); ``title`` is the graph title.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.plot_repeat_landscape(divergence_summary: 'str | Path', output_directory: 'str | Path', *, genome_size: 'int | None' = None, two_bit_file: 'str | Path | None' = None, title: 'str | None' = None, javascript_only: 'bool' = False, output_name: 'str' = 'landscape.html', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `divergence_summary`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``classes`` (repeat classes) and ``landscape`` (one row per divergence level with ``divergence`` and the percentage of the genome per class), which is the data behind the plot.

## 19) `dnakit.annotation.repeatmasker_p11_cli_functions.process_repeat_hits()`

- **Purpose:** Post-process raw RepeatMasker search hits (``.cat``) into the repeat annotation (``.out``) and summary (``.tbl``).
- **Method:** ``cat_file`` is the ``.cat`` (or ``.cat.gz``) file of a RepeatMasker run. ProcessRepeats resolves overlapping hits, joins fragments of one element and writes the annotation again; use it to recreate or reformat results without repeating the search. ``species`` (``-species``, RepeatMasker's default is human) or ``library`` (``-lib``) must match the original run. ``mask_low_complexity=False`` is ``-nolow``, ``mask_interspersed=False`` is ``-noint``; ``lowercase_ambiguous`` is ``-lcambig``; extra files: ``untouched_annotation`` (``-u``), ``crossmatch_output`` (``-xm``), ``ace_output``, ``gff_output`` (``-gff``), ``polymorphic_simple_repeats`` (``-poly``), ``alignments`` (``-a``, needs the alignment-bearing ``.cat``); ``no_id_column`` (``-no_id``), ``exclude_n_runs`` (``-excln``), ``align_l1_orf2`` (``-orf2``). With ``mask_source`` (a FASTA file, ``-maskSource``) the sequence is masked in ``mask_style`` ``n``, ``x`` or ``soft``.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.process_repeat_hits(cat_file: 'str | Path', output_directory: 'str | Path', *, species: 'str | None' = None, library: 'str | Path | None' = None, mask_low_complexity: 'bool' = True, mask_interspersed: 'bool' = True, lowercase_ambiguous: 'bool' = False, untouched_annotation: 'bool' = False, crossmatch_output: 'bool' = False, ace_output: 'bool' = False, gff_output: 'bool' = False, polymorphic_simple_repeats: 'bool' = False, no_id_column: 'bool' = False, exclude_n_runs: 'bool' = False, align_l1_orf2: 'bool' = False, alignments: 'bool' = False, mask_source: 'str | Path | None' = None, mask_style: 'str' = 'n', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `cat_file`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``annotations`` and ``class_counts`` (from ``.out``) and ``summary`` (from ``.tbl``).

## 20) `dnakit.annotation.repeatmasker_p11_cli_functions.recalculate_repeat_divergence()`

- **Purpose:** Recalculate the divergence of every repeat alignment and summarise it per family (``util/calcDivergenceFromAlign.pl``).
- **Method:** ``alignment_file`` is the ``.align`` (optionally gzipped) file of a RepeatMasker run made with ``-a``. The divergence is the Kimura 2-parameter distance, with CpG sites of the consensus counted specially (two transitions as one, one transition as 1/10) unless ``adjust_cpg=False`` (``-noCpGMod``). The summary table (``-s``, the input of ``plot_repeat_landscape``) is written to ``summary_name``; ``write_alignments`` (``-a``) also writes ``<name>.with_divergence.align`` with the divergence on each alignment.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.recalculate_repeat_divergence(alignment_file: 'str | Path', output_directory: 'str | Path', *, summary_name: 'str' = 'repeats.divsum', write_alignments: 'bool' = False, adjust_cpg: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignment_file`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``families``, one row per repeat family with ``class``, ``repeat``, ``absolute_length``, ``well_characterised_length`` and ``kimura_percent``.

## 21) `dnakit.annotation.repeatmasker_p11_cli_functions.summarize_repeat_annotation()`

- **Purpose:** Summarise a repeat annotation (``.out``) per family and class: counts, bases and genome fraction (``util/buildSummary.pl``).
- **Method:** ``min_divergence`` and ``max_divergence`` (``-minDiv`` / ``-maxDiv``, a fraction 0-1) filter the annotations. ``species`` (``-species``) splits families into lineage-specific and ancestral (for runs made with ``-species``; needs FamDB). ``genome`` (``-genome``) is a ``.2bit`` file or a two-column ``sequence<TAB>non-N size`` table that sets the genome size and the sequences counted; ``use_absolute_genome_size`` is ``-useAbsoluteGenomeSize``. ``library_dir`` is ``-libdir``. buildSummary.pl requires the UCSC program ``twoBitInfo`` on PATH.
- **API:** `dnakit.annotation.repeatmasker_p11_cli_functions.summarize_repeat_annotation(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'summary.txt', min_divergence: 'float | None' = None, max_divergence: 'float | None' = None, species: 'str | None' = None, genome: 'str | Path | None' = None, use_absolute_genome_size: 'bool' = False, library_dir: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `repeat_annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``report`` (the first 4000 characters of the summary).

