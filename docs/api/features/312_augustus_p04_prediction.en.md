# AUGUSTUS: gene prediction, training and pipelines

The functions in `dnakit.annotation` listed here are typed facades over AUGUSTUS 3.5.0 programs and helper scripts ([AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)). Each one runs a single fixed program with an explicit argument list (no shell). Install AUGUSTUS (`conda install -c bioconda augustus`) and put its `bin` folder on PATH, or pass `augustus_bin_dir`. If the program is missing the function raises `BackendUnavailableError` with the install hint. The species parameter sets are read from `AUGUSTUS_CONFIG_PATH` (or `config_path`); training functions write there, so the folder must be writable. Output paths must be new; nothing is overwritten. Each result is a dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and `provenance`, plus small parsed counts.

## 1) `dnakit.annotation.augustus_auto_predict()`

- **Purpose:** Chunked genome-wide gene prediction with joined output (``autoAugPred.pl --singleCPU``).
- **Method:** Runs `augustus_options`, `cpus`, `extrinsic_config_path`, `hints_path` through the DNAKit bounded executor. The genome is split into overlapping chunks that are predicted with AUGUSTUS (``cpus`` at a time) and the chunk predictions are joined into one non-redundant set. ``utr`` and ``hints_path`` (with ``extrinsic_config_path``) switch on UTR prediction and hint use. ``augustus_options`` are further ``--name=value`` AUGUSTUS options from the list ``--strand --genemodel --alternatives-from-evidence --alternatives-from-sampling --sample --softmasking --mea --protein --codingseq --gff3``. The results are written under the new ``work_dir``; the joined prediction files are listed in ``prediction_files``.
- **API:** `dnakit.annotation.augustus_auto_predict(genome_path, species, work_dir, *, hints_path=None, extrinsic_config_path=None, utr=False, augustus_options=(), use_existing=False, cpus=1, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `species`, `work_dir`. Optional keyword-only: `hints_path`, `extrinsic_config_path`, `utr`, `augustus_options`, `use_existing`, `cpus`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 2) `dnakit.annotation.augustus_auto_train_and_predict()`

- **Purpose:** Automated training and annotation pipeline (``autoAug.pl --singleCPU``).
- **Method:** Runs `cdna_path`, `cpus`, `est_alignment_path`, `genome_path` through the DNAKit bounded executor. Trains ``species`` on ``training_path`` (genes in GenBank, GFF or protein FASTA format) on ``genome_path``, optionally builds UTR models from ``est_alignment_path`` (PSL) or ``cdna_path`` (cDNA FASTA, aligned with BLAT), and predicts genes (with ``hints_path`` when given). The cluster mode of the script is not exposed; the pipeline always runs on this machine with ``cpus`` processes. Results and temporary files are written under the new ``work_dir``; the prediction files are listed in ``prediction_files``.
- **API:** `dnakit.annotation.augustus_auto_train_and_predict(genome_path, species, work_dir, *, training_path=None, hints_path=None, est_alignment_path=None, cdna_path=None, no_utr=False, optimization_rounds=None, max_intron_length=None, use_existing=False, cpus=1, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `species`, `work_dir`. Optional keyword-only: `training_path`, `hints_path`, `est_alignment_path`, `cdna_path`, `no_utr`, `optimization_rounds`, `max_intron_length`, `use_existing`, `cpus`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 3) `dnakit.annotation.augustus_choose_parameter_set()`

- **Purpose:** Choose the AUGUSTUS parameter set per sequence and predict genes (``augustify.py``).
- **Method:** Runs `augustify.py`, `classification_path`, `prediction_path`, `species` through the DNAKit bounded executor. ``mode="dominant_species"`` returns the majority-vote species (``species`` in the result) for the whole input. ``mode="metagenomic"`` writes a tab-separated table to the new ``classification_path`` that assigns every sequence to its best species with a probability, and with ``prediction_path`` also writes the ab initio gene predictions made with those parameter sets. ``species_list`` are the species parameter sets to try.
- **API:** `dnakit.annotation.augustus_choose_parameter_set(genome_path, species_list, *, mode='dominant_species', classification_path=None, prediction_path=None, threads=1, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `species_list`. Optional keyword-only: `mode`, `classification_path`, `prediction_path`, `threads`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 4) `dnakit.annotation.augustus_fix_in_frame_stop_genes()`

- **Purpose:** Replace genes with in-frame stop codons by genes re-predicted with AUGUSTUS MEA (``fix_in_frame_stop_codon_genes.py``).
- **Method:** Runs `.gff3`, `additional_augustus_arguments`, `annotation_path`, `bad_genes_path` through the DNAKit bounded executor. ``annotation_path`` is a GTF (or GFF3 with ``annotation_format="gff3"``) file, ``bad_genes_path`` lists the transcript identifiers of the genes with in-frame stops. The fixed annotation is written to ``<output_stem>.gtf`` (or ``.gff3``); ``output_stem`` must not exist as such a file yet. ``additional_augustus_arguments`` is a quoted string of further ``--name=value`` AUGUSTUS options (not ``--exonnames`` and not ``--alternatives-from-evidence=true``).
- **API:** `dnakit.annotation.augustus_fix_in_frame_stop_genes(genome_path, annotation_path, bad_genes_path, output_stem, species, *, annotation_format='gtf', hints_path=None, extrinsic_config_path=None, softmasking=False, utr=False, print_utr=False, additional_augustus_arguments=None, cores=1, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `annotation_path`, `bad_genes_path`, `output_stem`, `species`. Optional keyword-only: `annotation_format`, `hints_path`, `extrinsic_config_path`, `softmasking`, `utr`, `print_utr`, `additional_augustus_arguments`, `cores`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 5) `dnakit.annotation.augustus_homologous_gene_mapping()`

- **Purpose:** Map gene models across genomes with a HAL alignment (``homGeneMapping``).
- **Method:** Runs `gene_files_table_path`, `gtffilenames.tbl`, `halLiftover`, `hal_bin_dir` through the DNAKit bounded executor. ``gene_files_table_path`` is the ``gtffilenames.tbl`` text file (genome name, gene GTF, optional hints GTF per line). For every gene the program counts how many exons and introns agree with the genes of the other genomes and which genes are homologous; the annotated GTF files are written to the new ``output_dir``. ``homologs_path`` is a new file for the disjunct homolog sets. The HAL tool ``halLiftover`` must be on PATH or in ``hal_bin_dir``.
- **API:** `dnakit.annotation.augustus_homologous_gene_mapping(hal_path, gene_files_table_path, output_dir, *, no_dupes=False, details=False, unmapped=False, homologs_path=None, hal_bin_dir=None, cpus=1, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `hal_path`, `gene_files_table_path`, `output_dir`. Optional keyword-only: `no_dupes`, `details`, `unmapped`, `homologs_path`, `hal_bin_dir`, `cpus`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 6) `dnakit.annotation.augustus_join_gene_sets()`

- **Purpose:** Merge several GTF gene sets into one, joining and selecting transcripts (``joingenes``).
- **Method:** Runs `alternatives`, `bacterium`, `error_distance`, `eukaryote` through the DNAKit bounded executor. ``gene_set_paths`` are GTF files; ``priorities`` (positive integers, one per file) decide conflicts, bigger wins. ``error_distance`` (default 1000 in the program, negative disables) marks predictions near a sequence-range border; ``gene_model`` is ``eukaryote`` or ``bacterium``; ``alternatives`` joins genes whose transcripts are alternative variants; ``no_join`` and ``no_selection`` switch off the joining and the final selection. The merged GTF is written to the new ``output_path``.
- **API:** `dnakit.annotation.augustus_join_gene_sets(gene_set_paths, output_path, *, priorities=None, error_distance=None, gene_model=None, alternatives=False, suppress_priorities=None, stop_in_coding=False, no_join=False, no_selection=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gene_set_paths`, `output_path`. Optional keyword-only: `priorities`, `error_distance`, `gene_model`, `alternatives`, `suppress_priorities`, `stop_in_coding`, `no_join`, `no_selection`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_join_gene_sets([str(w / "genes.gtf")], str(w / "merged.gtf"), verbose=False)
print(r["genes"], r["transcripts"])
```

- **Example output:**

```text
1 1
```

## 7) `dnakit.annotation.augustus_join_segment_predictions()`

- **Purpose:** Join AUGUSTUS predictions of overlapping sequence segments into one renumbered set (``join_aug_pred.pl``).
- **Method:** Runs `drop_list_path`, `join_aug_pred.pl`, `prediction_paths` through the DNAKit bounded executor. ``prediction_paths`` are the AUGUSTUS outputs of the segments, in order; they are concatenated (the ``# This output was generated with AUGUSTUS`` header lines that the script needs are kept). One copy of every gene seen in two overlapping segments is dropped; ``drop_list_path`` lists gene identifiers to exclude.
- **API:** `dnakit.annotation.augustus_join_segment_predictions(prediction_paths, output_path, *, drop_list_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `prediction_paths`, `output_path`. Optional keyword-only: `drop_list_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
a = m.augustus_predict_genes(str(w / "genome.fa"), "human", str(w / "a.gtf"), verbose=False)
j = m.augustus_join_segment_predictions([a["outputs"]["predictions"]] * 2, str(w / "joined.gtf"), verbose=False)
print(j["input_files"], j["genes"] == a["genes"])
```

- **Example output:**

```text
2 True
```

## 8) `dnakit.annotation.augustus_optimize_parameters()`

- **Purpose:** Train AUGUSTUS and optimize its meta-parameters by k-fold cross-validated accuracy (``optimize_augustus.pl``).
- **Method:** Runs `augustus_train_species`, `only_train_path`, `optimize_augustus.pl`, `species` through the DNAKit bounded executor. ``species`` must already exist (see ``augustus_train_species``); ``training_path`` is a GenBank file of bona fide genes. Defaults are those of the script: 5 rounds, 8-fold cross validation, 5 test points per parameter. ``only_train_path`` adds genes used for etraining only. The optimized parameters are written into the species folder of the config directory. The script's temporary ``tmp_opt_<species>`` folder is created inside the new ``work_dir`` (and removed with ``cleanup=True``).
- **API:** `dnakit.annotation.augustus_optimize_parameters(species, training_path, work_dir, *, rounds=5, cpus=1, kfold=8, pstep=5, metapars_path=None, only_train_path=None, utr=False, sens_spec_bias=None, train_only_utr=False, no_train_pars=False, cleanup=False, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `species`, `training_path`, `work_dir`. Optional keyword-only: `rounds`, `cpus`, `kfold`, `pstep`, `metapars_path`, `only_train_path`, `utr`, `sens_spec_bias`, `train_only_utr`, `no_train_pars`, `cleanup`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 9) `dnakit.annotation.augustus_predict_comparative()`

- **Purpose:** Comparative gene prediction over a multi-species genome alignment (``augustus --treefile --alnfile``, CGP).
- **Method:** Runs `alignment_path`, `dd_iterations`, `dual_decomposition`, `output_path` through the DNAKit bounded executor. ``tree_path`` is a Newick tree, ``alignment_path`` a MAF alignment and ``species_filenames_path`` a two-column text file (species identifier, genome FASTA). ``species`` is the model species whose parameters (and CGP parameters) are used. ``dual_decomposition`` switches on the dual-decomposition decoder (``--/CompPred/dualdecomp=true``) with ``dd_iterations`` iterations. The predictions of all species are written to the new ``output_path``.
- **API:** `dnakit.annotation.augustus_predict_comparative(species, tree_path, alignment_path, species_filenames_path, output_path, *, dual_decomposition=False, dd_iterations=None, only_species=None, softmasking=None, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `species`, `tree_path`, `alignment_path`, `species_filenames_path`, `output_path`. Optional keyword-only: `dual_decomposition`, `dd_iterations`, `only_species`, `softmasking`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 10) `dnakit.annotation.augustus_predict_genes()`

- **Purpose:** Predict gene structures on genomic DNA with AUGUSTUS (HMM, optional hints; ``augustus --species=...``).
- **Method:** Runs `augustus_extract_annotation_fasta`, `coding_sequences`, `extrinsic_config_path`, `hints_path` through the DNAKit bounded executor. One function for the AUGUSTUS prediction modes: ab initio or hint-guided prediction (``hints_path`` with ``extrinsic_config_path``), alternative transcripts with posterior probabilities by sampling (``alternatives_from_sampling=True`` with ``sample=n``), maximum-expected-accuracy decoding (``mea=True``) and gene prediction guided by a protein-family profile (``protein_profile_path``, AUGUSTUS-PPX). ``softmasking=None`` keeps the species default. ``protein_sequences`` and ``coding_sequences`` add the ``# protein sequence`` and ``# coding sequence`` comment blocks (read them with ``augustus_extract_annotation_fasta``). The prediction (GTF-like text, or GFF3 with ``gff3=True``) is written to the new ``output_path``.
- **API:** `dnakit.annotation.augustus_predict_genes(sequence_path, species, output_path, *, strand='both', gene_model='partial', single_strand=False, hints_path=None, extrinsic_config_path=None, utr=False, softmasking=None, alternatives_from_evidence=None, alternatives_from_sampling=False, sample=None, min_exon_intron_prob=None, min_mean_exon_intron_prob=None, max_tracks=None, mea=False, protein_profile_path=None, gff3=False, protein_sequences=False, coding_sequences=False, prediction_start=None, prediction_end=None, no_in_frame_stop=False, unique_gene_id=False, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `sequence_path`, `species`, `output_path`. Optional keyword-only: `strand`, `gene_model`, `single_strand`, `hints_path`, `extrinsic_config_path`, `utr`, `softmasking`, `alternatives_from_evidence`, `alternatives_from_sampling`, `sample`, `min_exon_intron_prob`, `min_mean_exon_intron_prob`, `max_tracks`, `mea`, `protein_profile_path`, `gff3`, `protein_sequences`, `coding_sequences`, `prediction_start`, `prediction_end`, `no_in_frame_stop`, `unique_gene_id`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_predict_genes(str(w / "genome.fa"), "human", str(w / "pred.gtf"), mea=True, verbose=False)
print(r["status"], r["version"], r["genes"] >= 0)
print(r["command"][1:5])
```

- **Example output:**

```text
completed 3.5.0 True
['--species=human', '--strand=both', '--genemodel=partial', '--UTR=off']
```

## 11) `dnakit.annotation.augustus_predict_utrs_from_rnaseq()`

- **Purpose:** Predict UTRs of annotated genes from RNA-Seq coverage (``utrrnaseq``).
- **Method:** Runs `coding_regions_path`, `introns_path`, `output_path`, `repeats_path` through the DNAKit bounded executor. Starting from the start and stop codons in ``coding_regions_path`` (GTF/GFF), the coverage in ``wiggle_path`` downstream/upstream is monitored for a drastic drop; spliced alignments (``introns_path``) are incorporated as UTR introns and repeat-masked drops (``repeats_path``) are ignored. Parameter defaults are those of the program. The UTR models with RNA-Seq evidence are written to the new ``output_path`` (GFF).
- **API:** `dnakit.annotation.augustus_predict_utrs_from_rnaseq(genome_path, coding_regions_path, introns_path, wiggle_path, output_path, *, repeats_path=None, smoothing_window_size=150, read_length=150, limit=5000, drop_window_size=50, minimum_length=2, minimum_average_coverage=10, percent_window=0.6, percent_intron=0.5, percent_multiplicity=0.1, splice_sites='GT_AG', zero_coverage=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `coding_regions_path`, `introns_path`, `wiggle_path`, `output_path`. Optional keyword-only: `repeats_path`, `smoothing_window_size`, `read_length`, `limit`, `drop_window_size`, `minimum_length`, `minimum_average_coverage`, `percent_window`, `percent_intron`, `percent_multiplicity`, `splice_sites`, `zero_coverage`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 12) `dnakit.annotation.augustus_protein_profile_block_search()`

- **Purpose:** Locate the blocks of a protein-family profile in genomic DNA (``fastBlockSearch``, AUGUSTUS-PPX).
- **Method:** Runs `.prfl`, `fastBlockSearch`, `msa2prfl.pl`, `output_path` through the DNAKit bounded executor. ``profile_path`` is a ``.prfl`` block profile (see ``msa2prfl.pl``); ``species`` supplies the ``/ProteinModel/`` thresholds from its parameter file. The block hits that AUGUSTUS prints are written to the new ``output_path``.
- **API:** `dnakit.annotation.augustus_protein_profile_block_search(sequence_path, profile_path, species, output_path, *, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `sequence_path`, `profile_path`, `species`, `output_path`. Optional keyword-only: `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 13) `dnakit.annotation.augustus_train_species()`

- **Purpose:** Train gene-model parameters of a new species from training genes (GenBank, GFF or protein FASTA).
- **Method:** Runs `autoAugTrain.pl`, `cpus`, `crf`, `est_alignment_path` through the DNAKit bounded executor. ``method="etraining"`` creates the species folder with ``new_species.pl`` (``prokaryotic=True`` uses the prokaryotic template) and trains it with ``etraining --species`` on a GenBank training file; with ``set_stop_codon_frequencies=True`` it then runs ``setStopCodonFreqs.pl`` on the etraining log. ``method="autoAugTrain"`` runs ``autoAugTrain.pl`` (training set in GenBank, or GFF / protein FASTA with ``genome_path``; ``utr`` with ``est_alignment_path`` trains UTR models; ``optimization_rounds``, ``crf``, ``flanking_dna``, ``use_existing`` and ``cpus`` map to the script options; ``work_dir`` is a new folder). The parameter files are written into the species folder of the config directory, which must be writable.
- **API:** `dnakit.annotation.augustus_train_species(species, training_path, *, method='etraining', prokaryotic=False, set_stop_codon_frequencies=False, genome_path=None, utr=False, est_alignment_path=None, work_dir=None, optimization_rounds=None, crf=False, flanking_dna=None, use_existing=False, cpus=None, config_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `species`, `training_path`. Optional keyword-only: `method`, `prokaryotic`, `set_stop_codon_frequencies`, `genome_path`, `utr`, `est_alignment_path`, `work_dir`, `optimization_rounds`, `crf`, `flanking_dna`, `use_existing`, `cpus`, `config_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

