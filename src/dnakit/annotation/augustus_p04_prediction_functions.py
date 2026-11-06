"""AUGUSTUS gene prediction, training and pipeline facades (integration pack P04, part 1).

Each public function runs one AUGUSTUS 3.5.0 program or helper script with an explicit argument list: ``augustus``
(single genome, hints, MEA decoding, sampling, protein profiles and comparative mode), ``etraining`` and
``new_species.pl``, ``optimize_augustus.pl``, ``augustify.py``, ``autoAug.pl``, ``autoAugPred.pl``,
``autoAugTrain.pl``, ``joingenes`` and ``join_aug_pred.pl``, ``fix_in_frame_stop_codon_genes.py``,
``homGeneMapping``, ``utrrnaseq`` and ``fastBlockSearch``. The user installs AUGUSTUS (``conda install -c bioconda
augustus``) and either puts its ``bin`` folder on PATH or passes ``augustus_bin_dir``. The species parameter sets live in
the AUGUSTUS config folder (``config_path`` or the ``AUGUSTUS_CONFIG_PATH`` environment variable); training functions
write into it, so it must be writable. Outputs go to new explicit paths and nothing is overwritten. Results are
JSON-serialisable dictionaries with the status, backend, version, command, outputs, small parsed counts and a bounded
log tail. Defaults are the defaults of the AUGUSTUS command-line options.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _augustus_p04_cli as _a
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError

_STRANDS = ("both", "forward", "backward")
_GENE_MODELS = ("partial", "intronless", "complete", "atleastone", "exactlyone")
_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def _species(value: object, name: str = "species") -> str:
    return _a.name_token(value, name)


def _gtf_counts(path: Path) -> dict[str, int]:
    """Genes and transcripts in an AUGUSTUS GTF/GFF3 output (comment lines are ignored).

    Uses the ``gene``/``transcript``/``mRNA`` feature lines; when a file has none (``joingenes`` writes
    exon-level GTF only) the distinct ``gene_id`` and ``transcript_id`` attributes are counted.
    """
    genes = transcripts = 0
    gene_ids: set[str] = set()
    transcript_ids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "gene":
                genes += 1
            elif fields[2] in {"transcript", "mRNA"}:
                transcripts += 1
            gene_match = _GENE_ID.search(fields[8])
            transcript_match = _TRANSCRIPT_ID.search(fields[8])
            if gene_match:
                gene_ids.add(gene_match.group(1))
            if transcript_match:
                transcript_ids.add(transcript_match.group(1))
    return {"genes": genes or len(gene_ids), "transcripts": transcripts or len(transcript_ids)}


def _work_directory(value: object, name: str = "work_dir") -> Path:
    return _cli.new_directory(str(value), name)


def _list_file(values: Sequence[str], path: Path, name: str) -> Path:
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    for value in values:
        _a.name_token(value, name)
    path.write_text("\n".join(values) + "\n", encoding="utf-8")
    return path


def augustus_predict_genes(
    sequence_path: str,
    species: str,
    output_path: str,
    *,
    strand: str = "both",
    gene_model: str = "partial",
    single_strand: bool = False,
    hints_path: str | None = None,
    extrinsic_config_path: str | None = None,
    utr: bool = False,
    softmasking: bool | None = None,
    alternatives_from_evidence: bool | None = None,
    alternatives_from_sampling: bool = False,
    sample: int | None = None,
    min_exon_intron_prob: float | None = None,
    min_mean_exon_intron_prob: float | None = None,
    max_tracks: int | None = None,
    mea: bool = False,
    protein_profile_path: str | None = None,
    gff3: bool = False,
    protein_sequences: bool = False,
    coding_sequences: bool = False,
    prediction_start: int | None = None,
    prediction_end: int | None = None,
    no_in_frame_stop: bool = False,
    unique_gene_id: bool = False,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict gene structures on genomic DNA with AUGUSTUS (HMM, optional hints; ``augustus --species=...``).

    One function for the AUGUSTUS prediction modes: ab initio or hint-guided prediction (``hints_path`` with
    ``extrinsic_config_path``), alternative transcripts with posterior probabilities by sampling
    (``alternatives_from_sampling=True`` with ``sample=n``), maximum-expected-accuracy decoding (``mea=True``) and gene
    prediction guided by a protein-family profile (``protein_profile_path``, AUGUSTUS-PPX). ``softmasking=None``
    keeps the species default. ``protein_sequences`` and ``coding_sequences`` add the ``# protein sequence`` and
    ``# coding sequence`` comment blocks (read them with ``augustus_extract_annotation_fasta``). The prediction
    (GTF-like text, or GFF3 with ``gff3=True``) is written to the new ``output_path``.
    """
    source = _cli.existing_file(sequence_path, "sequence_path")
    target = _cli.new_file(output_path, "output_path")
    name = _species(species)
    arguments = [f"--species={name}", f"--strand={_cli.choice(strand, _STRANDS, 'strand')}"]
    arguments.append(f"--genemodel={_cli.choice(gene_model, _GENE_MODELS, 'gene_model')}")
    arguments += _a.switch("--singlestrand=true", _a.flag(single_strand, "single_strand"))
    if hints_path is not None:
        arguments.append(f"--hintsfile={_cli.existing_file(hints_path, 'hints_path')}")
    if extrinsic_config_path is not None:
        arguments.append(
            f"--extrinsicCfgFile={_cli.existing_file(extrinsic_config_path, 'extrinsic_config_path')}"
        )
    arguments.append(f"--UTR={'on' if _a.flag(utr, 'utr') else 'off'}")
    arguments += _a.boolean_word("--softmasking", softmasking, ("1", "0"))
    arguments += _a.boolean_word("--alternatives-from-evidence", alternatives_from_evidence)
    if _a.flag(alternatives_from_sampling, "alternatives_from_sampling"):
        arguments.append("--alternatives-from-sampling=true")
    elif sample is not None:
        raise ConfigurationError("sample needs alternatives_from_sampling=True.")
    arguments += _a.eq("--sample", _a.optional_int(sample, "sample", minimum=1))
    arguments += _a.eq(
        "--minexonintronprob",
        _a.optional_number(min_exon_intron_prob, "min_exon_intron_prob", minimum=0, maximum=1),
    )
    arguments += _a.eq(
        "--minmeanexonintronprob",
        _a.optional_number(
            min_mean_exon_intron_prob, "min_mean_exon_intron_prob", minimum=0, maximum=1
        ),
    )
    arguments += _a.eq("--maxtracks", _a.optional_int(max_tracks, "max_tracks", minimum=-1))
    arguments += _a.switch("--mea=1", _a.flag(mea, "mea"))
    if protein_profile_path is not None:
        arguments.append(
            f"--proteinprofile={_cli.existing_file(protein_profile_path, 'protein_profile_path')}"
        )
    arguments.append(f"--gff3={'on' if _a.flag(gff3, 'gff3') else 'off'}")
    arguments += _a.switch("--protein=on", _a.flag(protein_sequences, "protein_sequences"))
    arguments += _a.switch("--codingseq=on", _a.flag(coding_sequences, "coding_sequences"))
    arguments += _a.eq(
        "--predictionStart", _a.optional_int(prediction_start, "prediction_start", minimum=1)
    )
    arguments += _a.eq(
        "--predictionEnd", _a.optional_int(prediction_end, "prediction_end", minimum=1)
    )
    arguments += _a.switch("--noInFrameStop=true", _a.flag(no_in_frame_stop, "no_in_frame_stop"))
    arguments += _a.switch("--uniqueGeneId=true", _a.flag(unique_gene_id, "unique_gene_id"))
    arguments.append(str(source))
    run = _a.run(
        "augustus",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        stdout=target,
        config_path=config_path,
    )
    return _a.result("predict_genes", run, outputs={"predictions": target}, **_gtf_counts(target))


def augustus_predict_comparative(
    species: str,
    tree_path: str,
    alignment_path: str,
    species_filenames_path: str,
    output_path: str,
    *,
    dual_decomposition: bool = False,
    dd_iterations: int | None = None,
    only_species: int | None = None,
    softmasking: bool | None = None,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Comparative gene prediction over a multi-species genome alignment (``augustus --treefile --alnfile``, CGP).

    ``tree_path`` is a Newick tree, ``alignment_path`` a MAF alignment and ``species_filenames_path`` a two-column
    text file (species identifier, genome FASTA). ``species`` is the model species whose parameters (and CGP
    parameters) are used. ``dual_decomposition`` switches on the dual-decomposition decoder
    (``--/CompPred/dualdecomp=true``) with ``dd_iterations`` iterations. The predictions of all species are written
    to the new ``output_path``.
    """
    name = _species(species)
    tree = _cli.existing_file(tree_path, "tree_path")
    alignment = _cli.existing_file(alignment_path, "alignment_path")
    genomes = _cli.existing_file(species_filenames_path, "species_filenames_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        f"--species={name}",
        f"--treefile={tree}",
        f"--alnfile={alignment}",
        f"--speciesfilenames={genomes}",
    ]
    arguments += _a.switch(
        "--/CompPred/dualdecomp=true", _a.flag(dual_decomposition, "dual_decomposition")
    )
    arguments += _a.eq(
        "--/CompPred/dd_iterations", _a.optional_int(dd_iterations, "dd_iterations", minimum=1)
    )
    arguments += _a.eq(
        "--/CompPred/only_species", _a.optional_int(only_species, "only_species", minimum=0)
    )
    arguments += _a.boolean_word("--softmasking", softmasking, ("1", "0"))
    run = _a.run(
        "augustus",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        stdout=target,
        config_path=config_path,
    )
    return _a.result(
        "predict_comparative", run, outputs={"predictions": target}, **_gtf_counts(target)
    )


def augustus_protein_profile_block_search(
    sequence_path: str,
    profile_path: str,
    species: str,
    output_path: str,
    *,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate the blocks of a protein-family profile in genomic DNA (``fastBlockSearch``, AUGUSTUS-PPX).

    ``profile_path`` is a ``.prfl`` block profile (see ``msa2prfl.pl``); ``species`` supplies the ``/ProteinModel/``
    thresholds from its parameter file. The block hits that AUGUSTUS prints are written to the new ``output_path``.
    """
    source = _cli.existing_file(sequence_path, "sequence_path")
    profile = _cli.existing_file(profile_path, "profile_path")
    name = _species(species)
    target = _cli.new_file(output_path, "output_path")
    run = _a.run(
        "fastBlockSearch",
        [f"--species={name}", str(source), str(profile)],
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        stdout=target,
        config_path=config_path,
    )
    return _a.result(
        "protein_profile_block_search", run, outputs={"hits": target}, **_a.line_summary(target)
    )


def augustus_train_species(
    species: str,
    training_path: str,
    *,
    method: str = "etraining",
    prokaryotic: bool = False,
    set_stop_codon_frequencies: bool = False,
    genome_path: str | None = None,
    utr: bool = False,
    est_alignment_path: str | None = None,
    work_dir: str | None = None,
    optimization_rounds: int | None = None,
    crf: bool = False,
    flanking_dna: int | None = None,
    use_existing: bool = False,
    cpus: int | None = None,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train gene-model parameters of a new species from training genes (GenBank, GFF or protein FASTA).

    ``method="etraining"`` creates the species folder with ``new_species.pl`` (``prokaryotic=True`` uses the
    prokaryotic template) and trains it with ``etraining --species`` on a GenBank training file; with
    ``set_stop_codon_frequencies=True`` it then runs ``setStopCodonFreqs.pl`` on the etraining log.
    ``method="autoAugTrain"`` runs ``autoAugTrain.pl`` (training set in GenBank, or GFF / protein FASTA with
    ``genome_path``; ``utr`` with ``est_alignment_path`` trains UTR models; ``optimization_rounds``, ``crf``,
    ``flanking_dna``, ``use_existing`` and ``cpus`` map to the script options; ``work_dir`` is a new folder).
    The parameter files are written into the species folder of the config directory, which must be writable.
    """
    name = _species(species)
    training = _cli.existing_file(training_path, "training_path")
    chosen = _cli.choice(method, ("etraining", "autoAugTrain"), "method")
    directory = _a.bin_dir(augustus_bin_dir)
    if chosen == "etraining":
        scratch = Path(tempfile.mkdtemp(prefix="dnakit-augustus-train-"))
        try:
            new_species_args = [f"--species={name}", "--silent"]
            new_species_args += _a.switch("--prokaryotic", _a.flag(prokaryotic, "prokaryotic"))
            _a.run(
                "new_species.pl",
                new_species_args,
                directory=directory,
                cwd=scratch,
                config_path=config_path,
            )
            log_path = scratch / "etraining.log"
            run = _a.run(
                "etraining",
                [f"--species={name}", str(training)],
                directory=directory,
                cwd=scratch,
                stdout=log_path,
                config_path=config_path,
            )
            if _a.flag(set_stop_codon_frequencies, "set_stop_codon_frequencies"):
                config = _a.config_directory(config_path, _a.executable("etraining", directory))
                import os

                folder = config or os.environ.get("AUGUSTUS_CONFIG_PATH")
                if folder is None:
                    raise ConfigurationError("set_stop_codon_frequencies needs a config folder.")
                _a.run(
                    "setStopCodonFreqs.pl",
                    [
                        f"--AUGUSTUS_CONFIG_PATH={folder}",
                        f"--species={name}",
                        f"--etrainOut={log_path}",
                    ],
                    directory=directory,
                    cwd=scratch,
                    config_path=config_path,
                )
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        species_files = _species_files(directory, config_path, name)
        frequencies = {
            stop: float(value)
            for stop, value in re.findall(
                r"^(tag|taa|tga):\s*\d+\s*\(([0-9.]+)\)", log_text, re.MULTILINE
            )
        }
        return _a.result(
            "train_species",
            run,
            outputs={"species_folder": species_files["folder"]},
            method="etraining",
            species=name,
            parameter_files=species_files["files"],
            stop_codon_frequencies=frequencies,
            training_genes=_count_training_genes(training),
        )
    working = (
        _work_directory(work_dir)
        if work_dir is not None
        else Path(tempfile.mkdtemp(prefix="dnakit-augustus-"))
    )
    arguments = [f"--species={name}", f"--trainingset={training}", f"--workingdir={working}"]
    if genome_path is not None:
        arguments.append(f"--genome={_cli.existing_file(genome_path, 'genome_path')}")
    if _a.flag(utr, "utr"):
        arguments.append("--utr")
        if est_alignment_path is None:
            raise ConfigurationError("utr=True needs est_alignment_path (PSL alignments of cDNAs).")
        arguments.append(f"--estali={_cli.existing_file(est_alignment_path, 'est_alignment_path')}")
    arguments += _a.eq(
        "--optrounds", _a.optional_int(optimization_rounds, "optimization_rounds", minimum=0)
    )
    arguments += _a.switch("--CRF", _a.flag(crf, "crf"))
    arguments += _a.eq("--flanking_DNA", _a.optional_int(flanking_dna, "flanking_dna", minimum=0))
    arguments += _a.switch("--useexisting", _a.flag(use_existing, "use_existing"))
    arguments += _a.eq("--cpus", _a.optional_int(cpus, "cpus", minimum=1))
    log = working / "autoAugTrain.stdout"
    run = _a.run(
        "autoAugTrain.pl",
        arguments,
        directory=directory,
        cwd=working,
        stdout=log,
        config_path=config_path,
    )
    species_files = _species_files(directory, config_path, name)
    return _a.result(
        "train_species",
        run,
        outputs={"species_folder": species_files["folder"], "work_dir": working, "log": log},
        method="autoAugTrain",
        species=name,
        parameter_files=species_files["files"],
    )


def _species_files(directory: str | None, config_path: str | None, name: str) -> dict[str, Any]:
    import os

    config = _a.config_directory(
        config_path, _a.executable("augustus", directory)
    ) or os.environ.get("AUGUSTUS_CONFIG_PATH", "")
    folder = Path(config) / "species" / name
    files = sorted(path.name for path in folder.iterdir()) if folder.is_dir() else []
    return {"folder": folder, "files": files}


def _count_training_genes(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(
            1 for line in handle if line.lstrip().startswith("CDS") or line.startswith("     CDS")
        )


def augustus_optimize_parameters(
    species: str,
    training_path: str,
    work_dir: str,
    *,
    rounds: int = 5,
    cpus: int = 1,
    kfold: int = 8,
    pstep: int = 5,
    metapars_path: str | None = None,
    only_train_path: str | None = None,
    utr: bool = False,
    sens_spec_bias: float | None = None,
    train_only_utr: bool = False,
    no_train_pars: bool = False,
    cleanup: bool = False,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train AUGUSTUS and optimize its meta-parameters by k-fold cross-validated accuracy (``optimize_augustus.pl``).

    ``species`` must already exist (see ``augustus_train_species``); ``training_path`` is a GenBank file of bona
    fide genes. Defaults are those of the script: 5 rounds, 8-fold cross validation, 5 test points per parameter.
    ``only_train_path`` adds genes used for etraining only. The optimized parameters are written into the species
    folder of the config directory. The script's temporary ``tmp_opt_<species>`` folder is created inside the new
    ``work_dir`` (and removed with ``cleanup=True``).
    """
    name = _species(species)
    training = _cli.existing_file(training_path, "training_path")
    working = _work_directory(work_dir)
    arguments = [f"--species={name}"]
    arguments.append(f"--rounds={_cli.integer(rounds, 'rounds', minimum=1)}")
    arguments.append(f"--cpus={_cli.integer(cpus, 'cpus', minimum=1)}")
    arguments.append(f"--kfold={_cli.integer(kfold, 'kfold', minimum=2)}")
    arguments.append(f"--pstep={_cli.integer(pstep, 'pstep', minimum=1)}")
    if metapars_path is not None:
        arguments.append(f"--metapars={_cli.existing_file(metapars_path, 'metapars_path')}")
    if only_train_path is not None:
        arguments.append(f"--onlytrain={_cli.existing_file(only_train_path, 'only_train_path')}")
    arguments += _a.switch("--UTR=on", _a.flag(utr, "utr"))
    arguments += _a.eq(
        "--sens_spec_bias", _a.optional_number(sens_spec_bias, "sens_spec_bias", minimum=0)
    )
    arguments += _a.switch("--trainOnlyUtr=1", _a.flag(train_only_utr, "train_only_utr"))
    arguments += _a.switch("--noTrainPars=1", _a.flag(no_train_pars, "no_train_pars"))
    arguments += _a.switch("--cleanup=1", _a.flag(cleanup, "cleanup"))
    arguments.append(str(training))
    log = working / "optimize_augustus.log"
    directory = _a.bin_dir(augustus_bin_dir)
    run = _a.run(
        "optimize_augustus.pl",
        arguments,
        directory=directory,
        cwd=working,
        stdout=log,
        config_path=config_path,
    )
    species_files = _species_files(directory, config_path, name)
    text = log.read_text(encoding="utf-8", errors="replace")
    accuracies = [
        float(value)
        for value in re.findall(r"accuracy[^\d-]*([0-9]*\.?[0-9]+)", text, re.IGNORECASE)
    ]
    return _a.result(
        "optimize_parameters",
        run,
        outputs={"species_folder": species_files["folder"], "work_dir": working, "log": log},
        species=name,
        parameter_files=species_files["files"],
        accuracy_values=accuracies[-10:],
    )


def augustus_choose_parameter_set(
    genome_path: str,
    species_list: Sequence[str],
    *,
    mode: str = "dominant_species",
    classification_path: str | None = None,
    prediction_path: str | None = None,
    threads: int = 1,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Choose the AUGUSTUS parameter set per sequence and predict genes (``augustify.py``).

    ``mode="dominant_species"`` returns the majority-vote species (``species`` in the result) for the whole input.
    ``mode="metagenomic"`` writes a tab-separated table to the new ``classification_path`` that assigns every
    sequence to its best species with a probability, and with ``prediction_path`` also writes the ab initio gene
    predictions made with those parameter sets. ``species_list`` are the species parameter sets to try.
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    chosen = _cli.choice(mode, ("dominant_species", "metagenomic"), "mode")
    directory = _a.bin_dir(augustus_bin_dir)
    threads = _cli.integer(threads, "threads", minimum=1)
    if chosen == "dominant_species" and (
        classification_path is not None or prediction_path is not None
    ):
        raise ConfigurationError("classification_path and prediction_path need mode='metagenomic'.")
    if chosen == "metagenomic" and classification_path is None:
        raise ConfigurationError("mode='metagenomic' needs classification_path.")
    if prediction_path is not None and classification_path is None:
        raise ConfigurationError("prediction_path needs classification_path.")
    with tempfile.TemporaryDirectory(prefix="dnakit-augustify-") as scratch:
        scratch_dir = Path(scratch)
        names = _list_file(species_list, scratch_dir / "species.lst", "species_list")
        arguments = ["-g", str(genome), "-p", str(names), "-t", str(threads)]
        outputs: dict[str, Path] = {}
        if chosen == "dominant_species":
            arguments.append("-s")
            run = _a.run(
                "augustify.py",
                arguments,
                directory=directory,
                cwd=scratch_dir,
                capture=True,
                config_path=config_path,
            )
            text = _a.captured_text(run)
            lines = [line.strip() for line in text.splitlines() if line.strip()]
            return _a.result(
                "choose_parameter_set",
                run,
                outputs=outputs,
                mode=chosen,
                species=lines[-1].split(":", 1)[-1].strip() if lines else None,
            )
        classification = _cli.new_file(classification_path, "classification_path")  # type: ignore[arg-type]
        arguments += ["-m", str(classification)]
        outputs["classification"] = classification
        if prediction_path is not None:
            prediction = _cli.new_file(prediction_path, "prediction_path")
            arguments += ["-P", str(prediction)]
            outputs["predictions"] = prediction
        run = _a.run(
            "augustify.py",
            arguments,
            directory=directory,
            cwd=classification.parent,
            capture=True,
            expected=list(outputs.values()),
            config_path=config_path,
        )
        _a.captured_text(run)
    rows = [
        line.split("\t")
        for line in classification.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return _a.result(
        "choose_parameter_set",
        run,
        outputs=outputs,
        mode=chosen,
        sequences_classified=len(rows),
        species_assigned=sorted({row[1] for row in rows if len(row) > 1}),
    )


def augustus_auto_train_and_predict(
    genome_path: str,
    species: str,
    work_dir: str,
    *,
    training_path: str | None = None,
    hints_path: str | None = None,
    est_alignment_path: str | None = None,
    cdna_path: str | None = None,
    no_utr: bool = False,
    optimization_rounds: int | None = None,
    max_intron_length: int | None = None,
    use_existing: bool = False,
    cpus: int = 1,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Automated training and annotation pipeline (``autoAug.pl --singleCPU``).

    Trains ``species`` on ``training_path`` (genes in GenBank, GFF or protein FASTA format) on ``genome_path``,
    optionally builds UTR models from ``est_alignment_path`` (PSL) or ``cdna_path`` (cDNA FASTA, aligned with BLAT),
    and predicts genes (with ``hints_path`` when given). The cluster mode of the script is not exposed; the pipeline
    always runs on this machine with ``cpus`` processes. Results and temporary files are written under the new
    ``work_dir``; the prediction files are listed in ``prediction_files``.
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    name = _species(species)
    working = _work_directory(work_dir)
    arguments = [
        f"--species={name}",
        f"--genome={genome}",
        f"--workingdir={working}",
        "--singleCPU",
    ]
    arguments.append(f"--cpus={_cli.integer(cpus, 'cpus', minimum=1)}")
    if training_path is not None:
        arguments.append(f"--trainingset={_cli.existing_file(training_path, 'training_path')}")
    if hints_path is not None:
        arguments.append(f"--hints={_cli.existing_file(hints_path, 'hints_path')}")
    if est_alignment_path is not None:
        arguments.append(f"--estali={_cli.existing_file(est_alignment_path, 'est_alignment_path')}")
    if cdna_path is not None:
        arguments.append(f"--cdna={_cli.existing_file(cdna_path, 'cdna_path')}")
    arguments += _a.switch("--noutr", _a.flag(no_utr, "no_utr"))
    arguments += _a.eq(
        "--optrounds", _a.optional_int(optimization_rounds, "optimization_rounds", minimum=0)
    )
    arguments += _a.eq(
        "--maxIntronLen", _a.optional_int(max_intron_length, "max_intron_length", minimum=1)
    )
    arguments += _a.switch("--useexisting", _a.flag(use_existing, "use_existing"))
    log = working / "autoAug.stdout"
    run = _a.run(
        "autoAug.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=working,
        stdout=log,
        config_path=config_path,
    )
    predictions = sorted(str(path) for path in working.rglob("augustus.gff"))
    return _a.result(
        "auto_train_and_predict",
        run,
        outputs={"work_dir": working, "log": log},
        species=name,
        prediction_files=predictions,
    )


def augustus_auto_predict(
    genome_path: str,
    species: str,
    work_dir: str,
    *,
    hints_path: str | None = None,
    extrinsic_config_path: str | None = None,
    utr: bool = False,
    augustus_options: Sequence[str] = (),
    use_existing: bool = False,
    cpus: int = 1,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Chunked genome-wide gene prediction with joined output (``autoAugPred.pl --singleCPU``).

    The genome is split into overlapping chunks that are predicted with AUGUSTUS (``cpus`` at a time) and the
    chunk predictions are joined into one non-redundant set. ``utr`` and ``hints_path`` (with
    ``extrinsic_config_path``) switch on UTR prediction and hint use. ``augustus_options`` are further
    ``--name=value`` AUGUSTUS options from the list ``--strand --genemodel --alternatives-from-evidence
    --alternatives-from-sampling --sample --softmasking --mea --protein --codingseq --gff3``. The results are
    written under the new ``work_dir``; the joined prediction files are listed in ``prediction_files``.
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    name = _species(species)
    working = _work_directory(work_dir)
    allowed = {
        "--strand",
        "--genemodel",
        "--alternatives-from-evidence",
        "--alternatives-from-sampling",
        "--sample",
        "--softmasking",
        "--mea",
        "--protein",
        "--codingseq",
        "--gff3",
    }
    options = []
    for option in augustus_options:
        key = str(option).split("=", 1)[0]
        if key not in allowed or any(ch in str(option) for ch in " \t\n;|&$`"):
            raise ConfigurationError(
                f"augustus_options entry {option!r} is not allowed; choose from {sorted(allowed)}."
            )
        options.append(str(option))
    arguments = [
        f"--species={name}",
        f"--genome={genome}",
        f"--workingdir={working}",
        "--singleCPU",
    ]
    arguments.append(f"--cpus={_cli.integer(cpus, 'cpus', minimum=1)}")
    if hints_path is not None:
        arguments.append(f"--hints={_cli.existing_file(hints_path, 'hints_path')}")
    if extrinsic_config_path is not None:
        arguments.append(
            f"--extrinsiccfg={_cli.existing_file(extrinsic_config_path, 'extrinsic_config_path')}"
        )
    arguments += _a.switch("--utr", _a.flag(utr, "utr"))
    arguments += _a.switch("--useexisting", _a.flag(use_existing, "use_existing"))
    if options:
        arguments.append("--augustusoptions=" + " ".join(options))
    log = working / "autoAugPred.stdout"
    run = _a.run(
        "autoAugPred.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=working,
        stdout=log,
        config_path=config_path,
    )
    predictions = sorted(str(path) for path in working.rglob("augustus.gff"))
    return _a.result(
        "auto_predict",
        run,
        outputs={"work_dir": working, "log": log},
        species=name,
        prediction_files=predictions,
    )


def augustus_join_segment_predictions(
    prediction_paths: Sequence[str],
    output_path: str,
    *,
    drop_list_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join AUGUSTUS predictions of overlapping sequence segments into one renumbered set (``join_aug_pred.pl``).

    ``prediction_paths`` are the AUGUSTUS outputs of the segments, in order; they are concatenated (the
    ``# This output was generated with AUGUSTUS`` header lines that the script needs are kept). One copy of every
    gene seen in two overlapping segments is dropped; ``drop_list_path`` lists gene identifiers to exclude.
    """
    if isinstance(prediction_paths, str) or not prediction_paths:
        raise ConfigurationError(
            "prediction_paths must be a non-empty list of AUGUSTUS output files."
        )
    files = [_cli.existing_file(path, "prediction_paths") for path in prediction_paths]
    target = _cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    if drop_list_path is not None:
        arguments.append(f"--droplist={_cli.existing_file(drop_list_path, 'drop_list_path')}")
    with tempfile.TemporaryDirectory(prefix="dnakit-augustus-join-") as scratch:
        merged = Path(scratch) / "augustus.concat"
        with merged.open("wb") as out:
            for file in files:
                with file.open("rb") as handle:
                    shutil.copyfileobj(handle, out)
        run = _a.run(
            "join_aug_pred.pl",
            arguments,
            directory=_a.bin_dir(augustus_bin_dir),
            cwd=target.parent,
            stdin=merged,
            stdout=target,
        )
    return _a.result(
        "join_segment_predictions",
        run,
        outputs={"joined": target},
        input_files=len(files),
        **_gtf_counts(target),
    )


def augustus_join_gene_sets(
    gene_set_paths: Sequence[str],
    output_path: str,
    *,
    priorities: Sequence[int] | None = None,
    error_distance: int | None = None,
    gene_model: str | None = None,
    alternatives: bool = False,
    suppress_priorities: Sequence[int] | None = None,
    stop_in_coding: bool = False,
    no_join: bool = False,
    no_selection: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several GTF gene sets into one, joining and selecting transcripts (``joingenes``).

    ``gene_set_paths`` are GTF files; ``priorities`` (positive integers, one per file) decide conflicts, bigger
    wins. ``error_distance`` (default 1000 in the program, negative disables) marks predictions near a sequence-range
    border; ``gene_model`` is ``eukaryote`` or ``bacterium``; ``alternatives`` joins genes whose transcripts are
    alternative variants; ``no_join`` and ``no_selection`` switch off the joining and the final selection. The merged
    GTF is written to the new ``output_path``.
    """
    if isinstance(gene_set_paths, str) or not gene_set_paths:
        raise ConfigurationError("gene_set_paths must be a non-empty list of GTF files.")
    files = [str(_cli.existing_file(path, "gene_set_paths")) for path in gene_set_paths]
    target = _cli.new_file(output_path, "output_path")
    arguments = ["--genesets=" + ",".join(files), f"--output={target}"]
    if priorities is not None:
        if len(priorities) != len(files):
            raise ConfigurationError("priorities needs one positive integer per gene set.")
        arguments.append(
            "--priorities="
            + ",".join(str(_cli.integer(p, "priorities", minimum=1)) for p in priorities)
        )
    arguments += _a.eq("--errordistance", _a.optional_int(error_distance, "error_distance"))
    if gene_model is not None:
        arguments.append(
            f"--genemodel={_cli.choice(gene_model, ('eukaryote', 'bacterium'), 'gene_model')}"
        )
    arguments += _a.switch("--alternatives", _a.flag(alternatives, "alternatives"))
    if suppress_priorities:
        arguments.append(
            "--suppress="
            + ",".join(
                str(_cli.integer(p, "suppress_priorities", minimum=1)) for p in suppress_priorities
            )
        )
    arguments += _a.switch("--stopincoding", _a.flag(stop_in_coding, "stop_in_coding"))
    arguments += _a.switch("--nojoin", _a.flag(no_join, "no_join"))
    arguments += _a.switch("--noselection", _a.flag(no_selection, "no_selection"))
    run = _a.run(
        "joingenes",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        expected=[target],
    )
    return _a.result(
        "join_gene_sets",
        run,
        outputs={"merged": target},
        input_files=len(files),
        **_gtf_counts(target),
    )


def augustus_fix_in_frame_stop_genes(
    genome_path: str,
    annotation_path: str,
    bad_genes_path: str,
    output_stem: str,
    species: str,
    *,
    annotation_format: str = "gtf",
    hints_path: str | None = None,
    extrinsic_config_path: str | None = None,
    softmasking: bool = False,
    utr: bool = False,
    print_utr: bool = False,
    additional_augustus_arguments: str | None = None,
    cores: int = 1,
    config_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace genes with in-frame stop codons by genes re-predicted with AUGUSTUS MEA (``fix_in_frame_stop_codon_genes.py``).

    ``annotation_path`` is a GTF (or GFF3 with ``annotation_format="gff3"``) file, ``bad_genes_path`` lists the
    transcript identifiers of the genes with in-frame stops. The fixed annotation is written to ``<output_stem>.gtf``
    (or ``.gff3``); ``output_stem`` must not exist as such a file yet. ``additional_augustus_arguments`` is a
    quoted string of further ``--name=value`` AUGUSTUS options (not ``--exonnames`` and not
    ``--alternatives-from-evidence=true``).
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    annotation = _cli.existing_file(annotation_path, "annotation_path")
    bad = _cli.existing_file(bad_genes_path, "bad_genes_path")
    fmt = _cli.choice(annotation_format, ("gtf", "gff3"), "annotation_format")
    stem = Path(output_stem).expanduser().resolve()
    if not stem.parent.is_dir():
        raise ConfigurationError(f"The folder for output_stem does not exist: {stem.parent}")
    target = stem.with_name(stem.name + (".gtf" if fmt == "gtf" else ".gff3"))
    if target.exists():
        raise ConfigurationError(f"output_stem would overwrite {target}.")
    name = _species(species)
    arguments = [
        "-g",
        str(genome),
        "-t" if fmt == "gtf" else "-3",
        str(annotation),
        "-b",
        str(bad),
        "-o",
        str(stem),
        "-s",
        name,
    ]
    if hints_path is not None:
        arguments += ["-H", str(_cli.existing_file(hints_path, "hints_path"))]
    if extrinsic_config_path is not None:
        arguments += ["-e", str(_cli.existing_file(extrinsic_config_path, "extrinsic_config_path"))]
    arguments += ["-m", "on" if _a.flag(softmasking, "softmasking") else "off"]
    arguments += ["-u", "on" if _a.flag(utr, "utr") else "off"]
    arguments += ["-U", "on" if _a.flag(print_utr, "print_utr") else "off"]
    if additional_augustus_arguments is not None:
        extra = _cli.text(additional_augustus_arguments, "additional_augustus_arguments")
        for token in extra.split():
            if not token.startswith("--") or any(ch in token for ch in ";|&$`<>"):
                raise ConfigurationError(
                    "additional_augustus_arguments must be --name=value tokens."
                )
        arguments += [f"--additional_aug_args={extra} "]
    arguments += ["-c", str(_cli.integer(cores, "cores", minimum=1))]
    directory = _a.bin_dir(augustus_bin_dir)
    config = _a.config_directory(config_path, _a.executable("augustus", directory))
    if config is not None:
        arguments += ["-a", config]
    arguments += ["-A", str(Path(_a.executable("augustus", directory)).parent)]
    run = _a.run(
        "fix_in_frame_stop_codon_genes.py",
        arguments,
        directory=directory,
        cwd=stem.parent,
        capture=True,
        expected=[target],
        config_path=config_path,
    )
    _a.captured_text(run)
    return _a.result(
        "fix_in_frame_stop_genes", run, outputs={"annotation": target}, **_gtf_counts(target)
    )


def augustus_homologous_gene_mapping(
    hal_path: str,
    gene_files_table_path: str,
    output_dir: str,
    *,
    no_dupes: bool = False,
    details: bool = False,
    unmapped: bool = False,
    homologs_path: str | None = None,
    hal_bin_dir: str | None = None,
    cpus: int = 1,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map gene models across genomes with a HAL alignment (``homGeneMapping``).

    ``gene_files_table_path`` is the ``gtffilenames.tbl`` text file (genome name, gene GTF, optional hints GTF per
    line). For every gene the program counts how many exons and introns agree with the genes of the other genomes
    and which genes are homologous; the annotated GTF files are written to the new ``output_dir``.
    ``homologs_path`` is a new file for the disjunct homolog sets. The HAL tool ``halLiftover`` must be on PATH or in
    ``hal_bin_dir``.
    """
    hal = _cli.existing_file(hal_path, "hal_path")
    table = _cli.existing_file(gene_files_table_path, "gene_files_table_path")
    out = _work_directory(output_dir, "output_dir")
    arguments = [
        f"--gtfs={table}",
        f"--halfile={hal}",
        f"--outdir={out}",
        f"--tmpdir={out / 'tmp'}",
    ]
    arguments.append(f"--cpus={_cli.integer(cpus, 'cpus', minimum=1)}")
    arguments += _a.switch("--noDupes", _a.flag(no_dupes, "no_dupes"))
    arguments += _a.switch("--details", _a.flag(details, "details"))
    arguments += _a.switch("--unmapped", _a.flag(unmapped, "unmapped"))
    outputs: dict[str, Path] = {"output_dir": out}
    if homologs_path is not None:
        homologs = _cli.new_file(homologs_path, "homologs_path")
        arguments.append(f"--printHomologs={homologs}")
        outputs["homologs"] = homologs
    if hal_bin_dir is not None:
        folder = Path(hal_bin_dir).expanduser()
        if not folder.is_dir():
            raise ConfigurationError(f"hal_bin_dir is not a folder: {folder}")
        arguments.append(f"--halLiftover_exec_dir={folder.resolve()}")
    (out / "tmp").mkdir()
    log = out / "homGeneMapping.stdout"
    run = _a.run(
        "homGeneMapping", arguments, directory=_a.bin_dir(augustus_bin_dir), cwd=out, stdout=log
    )
    produced = sorted(
        path.name
        for path in out.iterdir()
        if path.suffix in {".gtf", ".gff"} or path.name.endswith(".out")
    )
    return _a.result("homologous_gene_mapping", run, outputs=outputs, annotated_files=produced)


def augustus_predict_utrs_from_rnaseq(
    genome_path: str,
    coding_regions_path: str,
    introns_path: str,
    wiggle_path: str,
    output_path: str,
    *,
    repeats_path: str | None = None,
    smoothing_window_size: int = 150,
    read_length: int = 150,
    limit: int = 5000,
    drop_window_size: int = 50,
    minimum_length: int = 2,
    minimum_average_coverage: int = 10,
    percent_window: float = 0.6,
    percent_intron: float = 0.5,
    percent_multiplicity: float = 0.1,
    splice_sites: str = "GT_AG",
    zero_coverage: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict UTRs of annotated genes from RNA-Seq coverage (``utrrnaseq``).

    Starting from the start and stop codons in ``coding_regions_path`` (GTF/GFF), the coverage in ``wiggle_path``
    downstream/upstream is monitored for a drastic drop; spliced alignments (``introns_path``) are incorporated as UTR
    introns and repeat-masked drops (``repeats_path``) are ignored. Parameter defaults are those of the program. The UTR
    models with RNA-Seq evidence are written to the new ``output_path`` (GFF).
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    coding = _cli.existing_file(coding_regions_path, "coding_regions_path")
    introns = _cli.existing_file(introns_path, "introns_path")
    wiggle = _cli.existing_file(wiggle_path, "wiggle_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        "--in-scaffold-file",
        str(genome),
        "--in-coding-region-file",
        str(coding),
        "--in-intron-file",
        str(introns),
        "--in-wiggle-file",
        str(wiggle),
        "--out-file-name",
        str(target),
    ]
    if repeats_path is not None:
        arguments += ["--in-repeat-file", str(_cli.existing_file(repeats_path, "repeats_path"))]
    arguments += [
        "--smoothing-window-size",
        str(_cli.integer(smoothing_window_size, "smoothing_window_size", minimum=1)),
    ]
    arguments += ["--read-length", str(_cli.integer(read_length, "read_length", minimum=1))]
    arguments += ["--limit", str(_cli.integer(limit, "limit", minimum=1))]
    arguments += [
        "--drop-window-size",
        str(_cli.integer(drop_window_size, "drop_window_size", minimum=1)),
    ]
    arguments += [
        "--minimum-length",
        str(_cli.integer(minimum_length, "minimum_length", minimum=1)),
    ]
    arguments += [
        "--minimum-average-coverage",
        str(_cli.integer(minimum_average_coverage, "minimum_average_coverage", minimum=0)),
    ]
    arguments += [
        "--percent-window",
        str(_cli.number(percent_window, "percent_window", minimum=0, maximum=1)),
    ]
    arguments += [
        "--percent-intron",
        str(_cli.number(percent_intron, "percent_intron", minimum=0, maximum=1)),
    ]
    arguments += [
        "--percent-multiplicity",
        str(_cli.number(percent_multiplicity, "percent_multiplicity", minimum=0, maximum=1)),
    ]
    arguments += [
        "--splices-sites",
        _cli.choice(splice_sites, ("GT_AG", "GC_AG", "all"), "splice_sites"),
    ]
    arguments += _a.switch("--zero-coverage", _a.flag(zero_coverage, "zero_coverage"))
    run = _a.run(
        "utrrnaseq",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        capture=True,
        expected=[target],
    )
    _a.captured_text(run)
    return _a.result(
        "predict_utrs_from_rnaseq", run, outputs={"utrs": target}, **_a.line_summary(target)
    )


__all__ = [
    "augustus_auto_predict",
    "augustus_auto_train_and_predict",
    "augustus_choose_parameter_set",
    "augustus_fix_in_frame_stop_genes",
    "augustus_homologous_gene_mapping",
    "augustus_join_gene_sets",
    "augustus_join_segment_predictions",
    "augustus_optimize_parameters",
    "augustus_predict_comparative",
    "augustus_predict_genes",
    "augustus_predict_utrs_from_rnaseq",
    "augustus_protein_profile_block_search",
    "augustus_train_species",
]
_register_public_verbose(globals())
