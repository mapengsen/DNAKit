"""DendroPy alignment, population-genetic and GenBank operations that DNAKit did not yet expose (pack P38).

Multi-block PHYLIP reading, HKY85 rate matrix and correction factor, stationary root sequences, simulation of
DNA sequences from fragmented populations, the folded site frequency spectrum, two-population summary
statistics, offline parsing of GenBank INSDSeq XML (records, FASTA, feature table, metadata),
building alignments from those records, concatenation of alignments into a partitioned supermatrix and
extraction of alignment columns.

Install with ``python -m pip install "dnakit[bio-dendropy]"``. Nothing here uses the network. Output files are
new files (an existing path is refused). Tested with dendropy 5.1.
"""

from __future__ import annotations

import random
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p38_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dendropy"
_MAX_SEQS = 5000
_MAX_LEN = 5_000_000
_MAX_FILE = 200_000_000
_READ_FORMATS = ("fasta", "nexus", "nexml", "phylip", "phylip-relaxed")
_WRITE_FORMATS = ("fasta", "nexus", "nexml", "phylip", "phylip-relaxed")
_DNA = frozenset("ACGTRYSWKMBDHVN-?.")


def _dendropy() -> Any:
    return require("dendropy", _EXTRA)


def _file_text(path: Path, name: str) -> str:
    if path.stat().st_size > _MAX_FILE:
        raise ConfigurationError(f"{name} is larger than {_MAX_FILE} bytes.")
    return path.read_text(encoding="utf-8", errors="replace")


def _schema(fmt: str) -> tuple[str, dict[str, Any]]:
    if fmt == "phylip":
        return "phylip", {"strict": True}
    if fmt == "phylip-relaxed":
        return "phylip", {"strict": False}
    return fmt, {}


def _read_matrix(dendropy: Any, path: Path, fmt: str, data_type: str, namespace: Any = None) -> Any:
    schema, options = _schema(h.choice(fmt, "input_format", _READ_FORMATS))
    kind = h.choice(data_type, "data_type", ("dna", "rna"))
    cls = dendropy.DnaCharacterMatrix if kind == "dna" else dendropy.RnaCharacterMatrix
    if schema in ("fasta", "phylip"):
        options["data_type"] = kind
    kwargs = {} if namespace is None else {"taxon_namespace": namespace}
    try:
        if schema in ("nexus", "nexml"):
            dataset = dendropy.DataSet.get(data=_file_text(path, "alignment"), schema=schema, **kwargs)
            if not dataset.char_matrices:
                raise ConfigurationError("The file holds no character matrix.")
            return dataset.char_matrices[0]
        return cls.get(data=_file_text(path, "alignment"), schema=schema, **kwargs, **options)
    except ConfigurationError:
        raise
    except Exception as exc:  # dendropy raises DataParseError and its subclasses
        raise ConfigurationError(f"dendropy could not read {path.name}: {exc}") from exc


def _matrix_from_mapping(dendropy: Any, sequences: Mapping[str, str], name: str) -> Any:
    if not isinstance(sequences, Mapping) or len(sequences) < 2 or len(sequences) > _MAX_SEQS:
        raise ConfigurationError(f"{name} must map at least two IDs to aligned DNA strings.")
    lengths = {len(str(v)) for v in sequences.values()}
    if len(lengths) != 1 or 0 in lengths or max(lengths) > _MAX_LEN:
        raise ConfigurationError(f"{name} sequences must have the same non-zero length.")
    cleaned = {}
    for key, value in sequences.items():
        text = str(value).upper()
        if set(text) - _DNA:
            raise ConfigurationError(f"{name} sequence {key!r} has symbols outside the IUPAC DNA alphabet.")
        cleaned[h.text(key, "sequence ID")] = text
    return dendropy.DnaCharacterMatrix.from_dict(cleaned)


def _matrix_rows(matrix: Any) -> dict[str, str]:
    return {str(t.label): matrix[t].symbols_as_string() for t in matrix.taxon_namespace}


def _states_list(matrix: Any) -> list[Any]:
    return [matrix[t] for t in matrix.taxon_namespace]


def _write(matrix_or_dataset: Any, target: Path, fmt: str) -> None:
    schema, options = _schema(h.choice(fmt, "output_format", _WRITE_FORMATS))
    with target.open("x", encoding="utf-8") as handle:
        matrix_or_dataset.write(file=handle, schema=schema, **options)


def read_multi_phylip_alignments(
    path: str,
    *,
    data_type: str = "dna",
    strict: bool = False,
    interleaved: bool = False,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a multi-block PHYLIP file (several alignments in one file) with dendropy's ``multiphylip`` reader.

    Every block starts with its own ``ntax nchar`` header line. ``data_type`` is ``dna`` or ``rna``; ``strict``
    selects strict (10-character names) rather than relaxed PHYLIP; ``interleaved`` reads interleaved blocks.
    Returns each alignment's taxa, length and sequences (at most ``limit`` blocks).
    """

    source = h.existing(path, "path")
    kind = h.choice(data_type, "data_type", ("dna", "rna"))
    cap = h.integer(limit, "limit", 1, 10_000)
    dendropy = _dendropy()
    try:
        dataset = dendropy.DataSet.get(
            data=_file_text(source, "path"), schema="multiphylip", data_type=kind,
            strict=h.flag(strict, "strict"), interleaved=h.flag(interleaved, "interleaved"),
        )
    except Exception as exc:  # dendropy raises DataParseError and TypeError
        raise ConfigurationError(f"dendropy could not read the multi-block PHYLIP file: {exc}") from exc
    blocks = []
    for index, matrix in enumerate(dataset.char_matrices[:cap]):
        rows = _matrix_rows(matrix)
        blocks.append({
            "index": index, "taxa": list(rows), "taxon_count": len(rows),
            "length": int(matrix.max_sequence_size), "sequences": rows,
        })
    return {
        "status": "completed",
        "operation": "read_multi_phylip_alignments",
        "data_type": kind,
        "alignment_count": len(dataset.char_matrices),
        "returned": len(blocks),
        "alignments": blocks,
        "provenance": h.provenance("dendropy", "dendropy.DataSet.read(schema='multiphylip')"),
    }


def _base_freqs(value: Sequence[float]) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 4:
        raise ConfigurationError("base_frequencies must be four numbers in the order A, C, G, T.")
    freqs = [h.number(v, "base frequency", 0.0, 1.0) for v in value]
    if abs(sum(freqs) - 1.0) > 1e-6 or min(freqs) <= 0:
        raise ConfigurationError("base_frequencies must be positive and sum to 1.")
    return freqs


def calculate_hky85_rate_matrix(
    *,
    base_frequencies: Sequence[float] = (0.25, 0.25, 0.25, 0.25),
    kappa: float = 1.0,
    rate: float = 1.0,
    branch_length: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """HKY85 instantaneous nucleotide rate matrix Q and the branch-length correction factor, using
    ``dendropy.model.discrete.Hky85.qmatrix`` and ``corrected_substitution_rate``.

    ``base_frequencies`` are ordered A, C, G, T; ``kappa`` is the transition/transversion rate ratio; ``rate`` is
    the relative rate. ``corrected_substitution_rate`` is the factor that converts a branch length to expected
    substitutions per site (``rate * 0.5 / f``). With ``branch_length`` the matrix of transition probabilities
    over that length (``pij``) is returned too.
    """

    dendropy = _dendropy()
    discrete = require("dendropy.model.discrete", _EXTRA)
    freqs = _base_freqs(base_frequencies)
    k = h.number(kappa, "kappa", 0.0, 1e6)
    r = h.number(rate, "rate", 0.0, 1e6)
    model = discrete.Hky85(kappa=k, base_freqs=freqs)
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "calculate_hky85_rate_matrix",
        "states": ["A", "C", "G", "T"],
        "base_frequencies": freqs,
        "kappa": k,
        "rate": r,
        "qmatrix": [[h.plain(float(v)) for v in row] for row in model.qmatrix(r)],
        "corrected_substitution_rate": h.plain(float(model.corrected_substitution_rate(r))),
        "provenance": h.provenance("dendropy", "dendropy.model.discrete.Hky85"),
    }
    if branch_length is not None:
        length = h.number(branch_length, "branch_length", 0.0, 1e6)
        result["branch_length"] = length
        result["transition_probabilities"] = [
            [h.plain(float(model.pij(i, j, length, r))) for j in range(4)] for i in range(4)
        ]
    del dendropy
    return result


def sample_stationary_nucleotide_sequence(
    length: int,
    *,
    base_frequencies: Sequence[float] = (0.25, 0.25, 0.25, 0.25),
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample a root nucleotide sequence from the stationary base distribution, using
    ``dendropy.model.discrete.NucleotideCharacterEvolutionModel.stationary_sample``.

    ``base_frequencies`` are ordered A, C, G, T; ``seed`` makes the draw reproducible.
    """

    _dendropy()
    discrete = require("dendropy.model.discrete", _EXTRA)
    n = h.integer(length, "length", 1, _MAX_LEN)
    freqs = _base_freqs(base_frequencies)
    rng = random.Random(None if seed is None else h.integer(seed, "seed", 0, 2**31 - 1))
    model = discrete.Hky85(kappa=1.0, base_freqs=freqs, rng=rng)
    states = model.stationary_sample(n, rng=rng)
    text = "".join(str(s.symbol) for s in states)
    return {
        "status": "completed",
        "operation": "sample_stationary_nucleotide_sequence",
        "length": n,
        "base_frequencies": freqs,
        "sequence": text,
        "provenance": h.provenance("dendropy", "NucleotideCharacterEvolutionModel.stationary_sample"),
    }


def simulate_fragmented_population_sequences(
    species_name: str,
    *,
    div_time_gens: int,
    num_desc_pops: int = 2,
    mutrate_per_site_per_generation: float = 1e-7,
    desc_pop_size: int = 10_000,
    samples_per_pop: int = 10,
    seq_len: int = 2000,
    seed: int | None = None,
    output_path: str | None = None,
    output_format: str = "fasta",
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate DNA sequences from fragmented (structured) populations, using
    ``dendropy.simulate.popgensim.FragmentedPopulations.generate_sequences``.

    A population tree with ``num_desc_pops`` descendant populations of ``desc_pop_size`` diverging
    ``div_time_gens`` generations ago is built, a gene tree is drawn with ``samples_per_pop`` genes per
    population, and ``seq_len`` sites are evolved on the mutation tree (dendropy's built-in HKY85, equal base
    frequencies). The sequences are returned and, with ``output_path``, written in ``output_format`` (fasta,
    nexus, nexml or phylip). ``seed`` makes the run reproducible.
    """

    _dendropy()
    popgensim = require("dendropy.simulate.popgensim", _EXTRA)
    name = h.text(species_name, "species_name")
    target = None if output_path is None else h.new_path(output_path, "output_path")
    params = {
        "div_time_gens": h.integer(div_time_gens, "div_time_gens", 1, 10**9),
        "num_desc_pops": h.integer(num_desc_pops, "num_desc_pops", 2, 20),
        "mutrate_per_site_per_generation": h.number(mutrate_per_site_per_generation, "mutrate_per_site_per_generation", 1e-12, 1.0),
        "desc_pop_size": h.integer(desc_pop_size, "desc_pop_size", 1, 10**9),
    }
    per_pop = h.integer(samples_per_pop, "samples_per_pop", 1, 200)
    sites = h.integer(seq_len, "seq_len", 1, 1_000_000)
    rng = random.Random(None if seed is None else h.integer(seed, "seed", 0, 2**31 - 1))
    try:
        simulator = popgensim.FragmentedPopulations(rng=rng, **params)
        dataset = simulator.generate_sequences(name, samples_per_pop=per_pop, seq_len=sites)
    except Exception as exc:  # dendropy raises ValueError, KeyError, AttributeError for bad parameters
        raise ConfigurationError(f"dendropy could not simulate sequences: {exc}") from exc
    matrix = dataset.char_matrices[0]
    rows = _matrix_rows(matrix)
    if target is not None:
        _write(matrix, target, output_format)
    return {
        "status": "completed",
        "operation": "simulate_fragmented_population_sequences",
        "species_name": name,
        "population_count": params["num_desc_pops"],
        "sequence_count": len(rows),
        "length": sites,
        "sequences": rows,
        "output_path": None if target is None else str(target),
        "provenance": h.provenance("dendropy", "dendropy.simulate.popgensim.FragmentedPopulations.generate_sequences"),
    }


def calculate_folded_site_frequency_spectrum(
    sequences: Mapping[str, str],
    *,
    pad_to_unfolded_length: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Folded (minor-allele) site frequency spectrum of a DNA alignment, using
    ``DiscreteCharacterMatrix.folded_site_frequency_spectrum``.

    ``sequences`` maps IDs to aligned DNA strings. Entry i counts the sites whose minor allele occurs in i
    sequences; the vector has ceil(N/2) entries, or N+1 with ``pad_to_unfolded_length``, whose first entry then
    counts monomorphic sites (as dendropy does).
    """

    dendropy = _dendropy()
    matrix = _matrix_from_mapping(dendropy, sequences, "sequences")
    pad = h.flag(pad_to_unfolded_length, "pad_to_unfolded_length")
    spectrum = [int(v) for v in matrix.folded_site_frequency_spectrum(pad)]
    return {
        "status": "completed",
        "operation": "calculate_folded_site_frequency_spectrum",
        "sequence_count": len(sequences),
        "alignment_length": int(matrix.max_sequence_size),
        "spectrum": spectrum,
        "padded": pad,
        "provenance": h.provenance("dendropy", "DiscreteCharacterMatrix.folded_site_frequency_spectrum"),
    }


def calculate_population_pair_statistics(
    population1: Mapping[str, str],
    population2: Mapping[str, str],
    *,
    ignore_uncertain: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Two-population pairwise summary statistics, using
    ``dendropy.calculate.popgenstat.PopulationPairSummaryStatistics`` (Hickerson et al. 2006, Wakeley 1996).

    Each population maps IDs to aligned DNA strings of the same length. Returns average pairwise differences
    (overall, between, within, net), segregating sites, Watterson's theta, Wakeley's psi and Tajima's D. Needs
    at least two sequences per population.
    """

    dendropy = _dendropy()
    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    first = dict(population1)
    second = dict(population2)
    if set(first) & set(second):
        raise ConfigurationError("The two populations must not share sequence IDs.")
    combined = _matrix_from_mapping(dendropy, {**first, **second}, "populations")
    rows = {k: v for k, v in _matrix_rows(combined).items()}
    one = [combined[t] for t in combined.taxon_namespace if str(t.label) in first]
    two = [combined[t] for t in combined.taxon_namespace if str(t.label) in second]
    if len(one) < 2 or len(two) < 2:
        raise ConfigurationError("Each population needs at least two sequences.")
    try:
        stats = popgen.PopulationPairSummaryStatistics(one, two, ignore_uncertain=h.flag(ignore_uncertain, "ignore_uncertain"))
    except (ZeroDivisionError, ValueError) as exc:
        raise ConfigurationError(f"Statistics are undefined for these sequences: {exc}") from exc
    del rows
    return {
        "status": "completed",
        "operation": "calculate_population_pair_statistics",
        "population_sizes": [len(one), len(two)],
        "average_number_of_pairwise_differences": h.plain(float(stats.average_number_of_pairwise_differences)),
        "average_number_of_pairwise_differences_between": h.plain(float(stats.average_number_of_pairwise_differences_between)),
        "average_number_of_pairwise_differences_within": h.plain(float(stats.average_number_of_pairwise_differences_within)),
        "average_number_of_pairwise_differences_net": h.plain(float(stats.average_number_of_pairwise_differences_net)),
        "num_segregating_sites": int(stats.num_segregating_sites),
        "wattersons_theta": h.plain(float(stats.wattersons_theta)),
        "wakeleys_psi": h.plain(float(stats.wakeleys_psi)),
        "tajimas_d": h.plain(float(stats.tajimas_d)),
        "provenance": h.provenance("dendropy", "dendropy.calculate.popgenstat.PopulationPairSummaryStatistics"),
    }


def _gbseq_store(dendropy: Any, xml_path: Path, kind: str) -> Any:
    genbank = require("dendropy.interop.genbank", _EXTRA)
    store = genbank.GenBankDna() if kind == "dna" else genbank.GenBankRna()
    try:
        store.read_xml_string(_file_text(xml_path, "xml_path"))
    except Exception as exc:  # xml.etree ParseError, AttributeError for non-GBSeq input
        raise ConfigurationError(f"dendropy could not parse the GenBank XML: {exc}") from exc
    return store


def parse_genbank_xml_records(
    xml_path: str,
    *,
    molecule: str = "dna",
    include_features: bool = True,
    fasta_labels: str = "defline",
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse GenBank INSDSeq XML records offline, using dendropy's GenBank interop classes.

    ``xml_path`` is an NCBI efetch ``rettype=gbc, retmode=xml`` file (INSDSet/INSDSeq XML, the format dendropy's
    GenBank interop reads; plain ``gb`` GBSeq XML is not read by dendropy). For each record the function returns
    the metadata (accession, version, organism, definition, taxonomy, dates, other sequence IDs), the sequence,
    its FASTA text (``GenBankAccessionRecord.as_fasta``; ``fasta_labels`` is ``defline`` for the NCBI-style
    ``gi|..|gb|..|`` header or ``accession_organism`` for a composed label) and, with ``include_features``, the
    feature table (key, location, qualifiers). ``molecule`` is ``dna`` or ``rna``.
    """

    source = h.existing(xml_path, "xml_path")
    kind = h.choice(molecule, "molecule", ("dna", "rna"))
    cap = h.integer(limit, "limit", 1, 10_000)
    mode = h.choice(fasta_labels, "fasta_labels", ("defline", "accession_organism"))
    features_wanted = h.flag(include_features, "include_features")
    dendropy = _dendropy()
    store = _gbseq_store(dendropy, source, kind)
    out = []
    for record in list(store)[:cap]:
        features = []
        if features_wanted:
            for feature in record.feature_table:
                features.append({
                    "key": feature.key,
                    "location": feature.location,
                    "qualifiers": {str(q.name): q.value for q in feature.qualifiers},
                })
        out.append({
            "gi": record.gi,
            "primary_accession": record.primary_accession,
            "accession_version": record.accession_version,
            "organism": record.organism,
            "definition": record.definition,
            "taxonomy": record.taxonomy,
            "moltype": record.moltype,
            "topology": record.topology,
            "length": record.length,
            "create_date": record.create_date,
            "update_date": record.update_date,
            "other_seq_ids": dict(record.other_seq_ids),
            "sequence": (record.sequence_text or "").upper(),
            "fasta": record.as_fasta(generate_new_label=(mode == "accession_organism")),
            "features": features,
        })
    return {
        "status": "completed",
        "operation": "parse_genbank_xml_records",
        "record_count": len(store),
        "returned": len(out),
        "records": out,
        "provenance": h.provenance("dendropy", "dendropy.interop.genbank GenBankAccessionRecord (offline)"),
    }


def genbank_xml_to_alignment(
    xml_path: str,
    output_path: str,
    *,
    output_format: str = "nexml",
    molecule: str = "dna",
    label_components: Sequence[str] = ("accession", "organism"),
    label_separator: str = "_",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a nucleotide character matrix from GenBank INSDSeq XML records and write it, using
    ``GenBankResourceStore.generate_char_matrix``.

    Taxon labels are composed from ``label_components`` (record attributes such as ``accession``, ``organism``,
    ``gi``) joined by ``label_separator``. The records need not be aligned: write ``fasta`` for unaligned sequences.
    (dendropy 5.1's GenBank annotation helpers ``as_annotation`` and ``as_reference_annotation`` raise ``TypeError``
    and are not used; metadata and feature tables are available from ``parse_genbank_xml_records``.)
    """

    source = h.existing(xml_path, "xml_path")
    target = h.new_path(output_path, "output_path")
    kind = h.choice(molecule, "molecule", ("dna", "rna"))
    fmt = h.choice(output_format, "output_format", _WRITE_FORMATS)
    components = [h.text(c, "label component") for c in label_components]
    separator = h.text(label_separator, "label_separator")
    dendropy = _dendropy()
    store = _gbseq_store(dendropy, source, kind)
    if len(store) == 0:
        raise ConfigurationError("The XML file holds no GenBank records.")
    try:
        matrix = store.generate_char_matrix(
            label_components=components, label_component_separator=separator,
        )
    except (ValueError, AttributeError) as exc:
        raise ConfigurationError(f"dendropy could not build the matrix: {exc}") from exc
    _write(matrix, target, fmt)
    lengths = {str(t.label): len(matrix[t]) for t in matrix.taxon_namespace}
    return {
        "status": "completed",
        "operation": "genbank_xml_to_alignment",
        "sequence_count": len(lengths),
        "sequence_lengths": lengths,
        "output_format": fmt,
        "output_path": str(target),
        "provenance": h.provenance("dendropy", "GenBankResourceStore.generate_char_matrix"),
    }


def concatenate_alignments(
    alignment_paths: Sequence[str],
    output_path: str,
    *,
    input_format: str = "fasta",
    output_format: str = "nexus",
    data_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate alignments into a partitioned supermatrix, using
    ``dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatenate`` (and ``concatenate_from_paths``).

    Every alignment must contain the same taxa (labels are matched by name). The parts are recorded as
    character subsets named after the input files; the ``nexus`` output keeps them as ``charset`` blocks. The
    result reports each part's column range (0-based, half-open).
    """

    if isinstance(alignment_paths, (str, bytes)) or len(alignment_paths) < 2:
        raise ConfigurationError("alignment_paths must list at least two alignment files.")
    paths = [h.existing(p, "alignment path") for p in alignment_paths]
    target = h.new_path(output_path, "output_path")
    fmt = h.choice(output_format, "output_format", _WRITE_FORMATS)
    dendropy = _dendropy()
    namespace = dendropy.TaxonNamespace()
    matrices = [_read_matrix(dendropy, p, input_format, data_type, namespace) for p in paths]
    try:
        merged = matrices[0].__class__.concatenate(matrices)
    except Exception as exc:  # dendropy raises ValueError for missing taxa
        raise ConfigurationError(f"dendropy could not concatenate the alignments: {exc}") from exc
    parts = []
    offset = 0
    for path, matrix in zip(paths, matrices):
        width = int(matrix.max_sequence_size)
        parts.append({"name": path.name, "start": offset, "end": offset + width})
        offset += width
    _write(merged, target, fmt)
    return {
        "status": "completed",
        "operation": "concatenate_alignments",
        "taxon_count": len(merged.taxon_namespace),
        "alignment_length": int(merged.max_sequence_size),
        "parts": parts,
        "output_format": fmt,
        "output_path": str(target),
        "provenance": h.provenance("dendropy", "CharacterMatrix.concatenate"),
    }


def extract_alignment_columns(
    alignment_path: str,
    output_path: str,
    *,
    columns: Sequence[int] | None = None,
    ranges: Sequence[Sequence[int]] | None = None,
    input_format: str = "fasta",
    output_format: str = "fasta",
    data_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract alignment columns by index or range, using ``CharacterMatrix.export_character_indices`` (and
    ``export_character_subset`` for named subsets).

    ``columns`` lists 0-based column indices; ``ranges`` lists ``[start, end)`` intervals (0-based, half-open);
    both can be given and are merged in alignment order. When the alignment is NEXUS with ``charset`` blocks,
    ``input_format="nexus"`` keeps them and ``ranges`` can be replaced by naming nothing: the function
    extracts only the requested columns.
    """

    source = h.existing(alignment_path, "alignment_path")
    target = h.new_path(output_path, "output_path")
    fmt = h.choice(output_format, "output_format", _WRITE_FORMATS)
    dendropy = _dendropy()
    matrix = _read_matrix(dendropy, source, input_format, data_type)
    width = int(matrix.max_sequence_size)
    chosen: set[int] = set()
    for item in columns or ():
        chosen.add(h.integer(item, "column index", 0, width - 1))
    for pair in ranges or ():
        if isinstance(pair, (str, bytes)) or len(pair) != 2:
            raise ConfigurationError("Each range must be [start, end).")
        start = h.integer(pair[0], "range start", 0, width)
        end = h.integer(pair[1], "range end", start, width)
        chosen.update(range(start, end))
    if not chosen:
        raise ConfigurationError("Give columns or ranges that select at least one column.")
    indices = sorted(chosen)
    subset = matrix.export_character_indices(indices)
    _write(subset, target, fmt)
    return {
        "status": "completed",
        "operation": "extract_alignment_columns",
        "original_length": width,
        "selected_columns": len(indices),
        "taxon_count": len(subset.taxon_namespace),
        "output_format": fmt,
        "output_path": str(target),
        "provenance": h.provenance("dendropy", "CharacterMatrix.export_character_indices"),
    }


__all__ = [
    "calculate_folded_site_frequency_spectrum",
    "calculate_hky85_rate_matrix",
    "calculate_population_pair_statistics",
    "concatenate_alignments",
    "extract_alignment_columns",
    "genbank_xml_to_alignment",
    "parse_genbank_xml_records",
    "read_multi_phylip_alignments",
    "sample_stationary_nucleotide_sequence",
    "simulate_fragmented_population_sequences",
]

_register_public_verbose(globals(), tuple(__all__))
