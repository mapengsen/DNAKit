"""Contract tests for optional sequence backends: Biopython, scikit-bio, pyfastx, screed, py2bit."""

from __future__ import annotations

import itertools
import struct
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequence.backends import (
    build_screed_database,
    calculate_kmer_frequencies,
    convert_sequence_file,
    detect_fastq_quality_encoding,
    expand_degenerate_sequences,
    extract_feature_sequences,
    fetch_screed_records,
    fetch_twobit_region,
    read_sequence_records,
    summarize_twobit,
    translate_dna_with_genetic_code,
)

# Forward CDS 1..18, minus-strand CDS 21..32, and a join() gene on 34..36 + 41..44 (1-based).
GENE_SEQUENCE = "ATGGCCAAATTTGGGCCCTAAGGATCCATGAAATTTTAACCGGTA"
FASTA_TEXT = ">r1 first read\nACGTNACGT\n>r2\nggcc\n"
FASTQ_TEXT = "@q1 desc\nACGTAC\n+\nIIII#!\n@q2\nGGCCAA\n+\nABCDEF\n"


def _gene_record():
    """Build the GenBank/EMBL fixture with Biopython; the caller has already imported Bio."""

    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq(GENE_SEQUENCE),
        id="AB000001.1",
        name="AB000001",
        description="test record for dnakit",
        annotations={"molecule_type": "DNA"},
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 18, strand=1),
            type="CDS",
            qualifiers={"gene": ["geneA"]},
        )
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(20, 32, strand=-1),
            type="CDS",
            qualifiers={"gene": ["geneB"]},
        )
    )
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [FeatureLocation(33, 36, strand=1), FeatureLocation(40, 44, strand=1)]
            ),
            type="gene",
            qualifiers={"gene": ["geneC"]},
        )
    )
    return record


def _write_gene_file(tmp_path: Path, file_format: str) -> Path:
    from Bio import SeqIO

    suffix = {"embl": "embl", "genbank": "gb"}[file_format]
    path = tmp_path / f"gene.{suffix}"
    SeqIO.write([_gene_record()], str(path), file_format)
    return path


def _write_text(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


_CODES = {"T": 0, "C": 1, "A": 2, "G": 3}


def _write_twobit(path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a little-endian UCSC .2bit file; lowercase is soft-masked and N is hard-masked."""

    def runs(flags: list[bool]) -> tuple[list[int], list[int]]:
        starts: list[int] = []
        sizes: list[int] = []
        index = 0
        while index < len(flags):
            if flags[index]:
                stop = index
                while stop < len(flags) and flags[stop]:
                    stop += 1
                starts.append(index)
                sizes.append(stop - index)
                index = stop
            else:
                index += 1
        return starts, sizes

    header = struct.pack("<IIII", 0x1A412743, 0, len(records), 0)
    index_bytes = b""
    body = b""
    offset = len(header) + sum(1 + len(name.encode()) + 4 for name, _ in records)
    for name, sequence in records:
        upper = sequence.upper()
        n_starts, n_sizes = runs([symbol == "N" for symbol in upper])
        mask_starts, mask_sizes = runs([symbol.islower() for symbol in sequence])
        packed = bytearray((len(sequence) + 3) // 4)
        for position, symbol in enumerate(upper):
            packed[position // 4] |= _CODES.get(symbol, 0) << (6 - 2 * (position % 4))
        record = struct.pack("<I", len(sequence))
        record += struct.pack("<I", len(n_starts))
        record += struct.pack(f"<{len(n_starts)}I", *n_starts)
        record += struct.pack(f"<{len(n_sizes)}I", *n_sizes)
        record += struct.pack("<I", len(mask_starts))
        record += struct.pack(f"<{len(mask_starts)}I", *mask_starts)
        record += struct.pack(f"<{len(mask_sizes)}I", *mask_sizes)
        record += struct.pack("<I", 0) + bytes(packed)
        raw_name = name.encode()
        index_bytes += bytes([len(raw_name)]) + raw_name + struct.pack("<I", offset)
        body += record
        offset += len(record)
    path.write_bytes(header + index_bytes + body)
    return path


TWOBIT_RECORDS = [("chr1", "NNNNACGTacgtGGCCN"), ("chr2", "TTTTAAAACCCCGGGG")]


# Biopython -----------------------------------------------------------------------------


def test_read_sequence_records_parses_embl_records_and_features(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    path = _write_gene_file(tmp_path, "embl")

    result = read_sequence_records(path, format="embl", verbose=False)

    assert result["record_count"] == 1
    record = result["records"][0]
    # 45 bases; two CDS features and one gene feature.
    assert record["id"] == "AB000001.1"
    assert record["description"] == "test record for dnakit"
    assert record["length"] == len(GENE_SEQUENCE) == 45
    assert record["sequence"] == GENE_SEQUENCE
    assert record["feature_counts"] == {"CDS": 2, "gene": 1}
    assert result["total_bases"] == 45


def test_read_sequence_records_enforces_explicit_limits(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    path = _write_gene_file(tmp_path, "genbank")
    with pytest.raises(ConfigurationError, match="max_bases=10"):
        read_sequence_records(path, format="genbank", max_bases=10, verbose=False)
    with pytest.raises(ConfigurationError, match="format must be one of"):
        read_sequence_records(path, format="fasta", verbose=False)
    with pytest.raises(ConfigurationError, match="max_records must be an integer"):
        read_sequence_records(path, format="genbank", max_records=True, verbose=False)


def test_extract_feature_sequences_reverse_complements_minus_strand(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    path = _write_gene_file(tmp_path, "genbank")

    result = extract_feature_sequences(path, format="genbank", verbose=False)

    first, second = result["features"]
    assert result["feature_count"] == 2
    assert (first["gene"], first["start"], first["end"], first["strand"]) == ("geneA", 1, 18, 1)
    assert first["sequence"] == "ATGGCCAAATTTGGGCCC"
    # Hand calculation: GENE_SEQUENCE[20:32] is AGGATCCATGAA, whose reverse complement is
    # TTCATGGATCCT.
    assert (second["gene"], second["start"], second["end"], second["strand"]) == (
        "geneB",
        21,
        32,
        -1,
    )
    assert second["sequence"] == "TTCATGGATCCT"


def test_extract_feature_sequences_joins_multi_part_locations(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    path = _write_gene_file(tmp_path, "embl")

    result = extract_feature_sequences(path, format="embl", feature_type="gene", verbose=False)

    (gene,) = result["features"]
    # join(34..36, 41..44): TTT from 34-36 followed by CGGT from 41-44.
    assert gene["parts"] == 2
    assert (gene["start"], gene["end"], gene["length"]) == (34, 44, 7)
    assert gene["sequence"] == "TTTCGGT"


def test_extract_feature_sequences_enforces_feature_limit(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    path = _write_gene_file(tmp_path, "genbank")
    with pytest.raises(ConfigurationError, match="max_features=1"):
        extract_feature_sequences(path, format="genbank", max_features=1, verbose=False)


def test_convert_sequence_file_writes_fasta_and_refuses_overwrite(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    source = _write_gene_file(tmp_path, "embl")
    target = tmp_path / "converted.fasta"

    result = convert_sequence_file(
        source, target, input_format="embl", output_format="fasta", verbose=False
    )

    assert result["record_count"] == 1
    assert target.read_text(encoding="utf-8").splitlines() == [
        ">AB000001.1 test record for dnakit",
        GENE_SEQUENCE,
    ]
    with pytest.raises(ConfigurationError, match="already exists"):
        convert_sequence_file(source, target, input_format="embl", verbose=False)

    # Round trip through GenBank keeps the sequence and both feature types.
    genbank = tmp_path / "converted.gb"
    written = convert_sequence_file(
        source, genbank, input_format="embl", output_format="genbank", verbose=False
    )
    reread = read_sequence_records(genbank, format="genbank", verbose=False)
    assert written["record_count"] == reread["record_count"] == 1
    assert reread["records"][0]["sequence"] == GENE_SEQUENCE
    assert reread["records"][0]["feature_counts"] == {"CDS": 2, "gene": 1}


def test_convert_sequence_file_enforces_record_limit_and_removes_partial_output(
    tmp_path: Path,
) -> None:
    pytest.importorskip("Bio")
    from Bio import SeqIO

    source = tmp_path / "two.embl"
    SeqIO.write([_gene_record(), _gene_record()], str(source), "embl")
    target = tmp_path / "partial.gb"
    with pytest.raises(ConfigurationError, match="max_records=1"):
        convert_sequence_file(
            source,
            target,
            input_format="embl",
            output_format="genbank",
            max_records=1,
            verbose=False,
        )
    assert not target.exists()


def test_biopython_missing_package_names_the_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_text(tmp_path, "empty.embl", "")
    monkeypatch.setitem(sys.modules, "Bio.SeqIO", None)
    with pytest.raises(BackendUnavailableError, match="bio-formats"):
        read_sequence_records(path, verbose=False)


# scikit-bio ----------------------------------------------------------------------------


def test_translate_uses_selected_ncbi_genetic_code() -> None:
    pytest.importorskip("skbio")
    # Standard code: ATG=M, TGA=stop, ATA=I, AGA=R.
    standard = translate_dna_with_genetic_code("ATGTGAATAAGA", verbose=False)
    assert standard["protein"] == "M*IR"
    assert standard["genetic_code"] == {"table_id": 1, "name": "Standard"}
    # Vertebrate mitochondrial code: TGA=W, ATA=M, AGA=stop.
    mitochondrial = translate_dna_with_genetic_code("ATGTGAATAAGA", genetic_code=2, verbose=False)
    assert mitochondrial["protein"] == "MWM*"
    assert mitochondrial["genetic_code"]["name"] == "Vertebrate Mitochondrial"


def test_translate_reverse_frame_uses_reverse_complement() -> None:
    pytest.importorskip("skbio")
    # Reverse complement of TTATTCACAGCATCC is GGATGCTGTGAATAA: GGA TGC TGT GAA TAA.
    result = translate_dna_with_genetic_code("TTATTCACAGCATCC", reading_frame=-1, verbose=False)
    assert result["protein"] == "GCCE*"
    assert result["reading_frame"] == -1


def test_translate_start_and_stop_policies() -> None:
    pytest.importorskip("skbio")
    stopped = translate_dna_with_genetic_code("ATGGCCTAA", stop="require", verbose=False)
    assert stopped["protein"] == "MA"
    with pytest.raises(ConfigurationError, match="could not translate"):
        translate_dna_with_genetic_code("CCCGGG", start="require", verbose=False)


def test_translate_rejects_gaps_bools_and_unassigned_tables() -> None:
    pytest.importorskip("skbio")
    with pytest.raises(ConfigurationError, match="unsupported symbols"):
        translate_dna_with_genetic_code("AC-GT", verbose=False)
    with pytest.raises(ConfigurationError, match="genetic_code must be an integer"):
        translate_dna_with_genetic_code("ATG", genetic_code=True, verbose=False)
    with pytest.raises(ConfigurationError, match="reading_frame must be one of"):
        translate_dna_with_genetic_code("ATG", reading_frame=4, verbose=False)
    with pytest.raises(ConfigurationError, match="could not translate"):
        translate_dna_with_genetic_code("ATG", genetic_code=7, verbose=False)


def test_kmer_frequencies_non_overlapping_relative_matches_hand_count() -> None:
    pytest.importorskip("skbio")
    # Non-overlapping 3-mers of ACACATTTATTA: ACA, CAT, TTA, TTA.
    result = calculate_kmer_frequencies(
        "ACACATTTATTA", k=3, overlap=False, relative=True, verbose=False
    )
    assert result["kmers"] == {"ACA": 0.25, "CAT": 0.25, "TTA": 0.5}
    assert result["distinct_kmers"] == 3


def test_kmer_frequencies_overlapping_counts_match_direct_scikit_bio() -> None:
    skbio = pytest.importorskip("skbio")
    # Overlapping windows: ACA, CAC, ACA, CAT, ATT, TTT, TTA, TAT, ATT, TTA.
    result = calculate_kmer_frequencies("ACACATTTATTA", k=3, verbose=False)
    expected = dict(sorted(skbio.DNA("ACACATTTATTA").kmer_frequencies(3).items()))
    assert result["kmers"] == expected
    assert result["kmers"]["ACA"] == 2
    assert result["kmers"]["TTA"] == 2
    assert sum(result["kmers"].values()) == 10


def test_expand_degenerate_sequences_enumerates_all_definite_sequences() -> None:
    pytest.importorskip("skbio")
    # R is A or G, so TRG expands to TAG and TGG.
    assert expand_degenerate_sequences("TRG", verbose=False)["sequences"] == ["TAG", "TGG"]
    # N is any base: TNR has 4 * 2 = 8 definite sequences.
    result = expand_degenerate_sequences("TNR", verbose=False)
    expected = sorted("T" + first + second for first, second in itertools.product("ACGT", "AG"))
    assert result["sequences"] == expected
    assert result["expansion_count"] == 8
    assert result["degenerate_positions"] == 2


def test_expand_degenerate_sequences_bounds_the_expansion() -> None:
    pytest.importorskip("skbio")
    # NNN expands to 64 sequences, above the explicit limit of 63.
    with pytest.raises(ConfigurationError, match="raise max_sequences"):
        expand_degenerate_sequences("NNN", max_sequences=63, verbose=False)


def test_scikit_bio_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "skbio", None)
    with pytest.raises(BackendUnavailableError, match="bio-skibio"):
        translate_dna_with_genetic_code("ATG", verbose=False)


# pyfastx -------------------------------------------------------------------------------


def test_detect_fastq_quality_encoding_reports_phred33_range(tmp_path: Path) -> None:
    pytest.importorskip("pyfastx")
    path = _write_text(tmp_path, "reads.fq", FASTQ_TEXT)

    result = detect_fastq_quality_encoding(path, verbose=False)

    # Quality characters: '!' (33) to 'I' (73), so scores 0 to 40 under offset 33.
    assert result["phred_offset"] == 33
    assert result["possible_phred_offsets"] == [33]
    assert (result["min_quality_ascii"], result["max_quality_ascii"]) == (33, 73)
    assert (result["min_quality_score"], result["max_quality_score"]) == (0, 40)
    assert (result["min_length"], result["max_length"]) == (6, 6)


def test_detect_fastq_quality_encoding_reports_phred64(tmp_path: Path) -> None:
    pytest.importorskip("pyfastx")
    path = _write_text(tmp_path, "illumina.fq", "@r1\nACGT\n+\nhhhh\n")

    result = detect_fastq_quality_encoding(path, verbose=False)

    # 'h' is 104; under offset 64 it is Phred 40.
    assert result["phred_offset"] == 64
    assert (result["min_quality_score"], result["max_quality_score"]) == (40, 40)


def test_detect_fastq_quality_encoding_reports_ambiguous_offset(tmp_path: Path) -> None:
    pytest.importorskip("pyfastx")
    path = _write_text(tmp_path, "ambiguous.fq", "@r1\nACGT\n+\nIIII\n")

    result = detect_fastq_quality_encoding(path, verbose=False)

    # 'I' (73) fits both offsets, so the wrapper reports both and no score range.
    assert result["phred_offset"] is None
    assert result["possible_phred_offsets"] == [33, 64]
    assert result["min_quality_score"] is None


def test_detect_fastq_quality_encoding_rejects_fasta_and_large_inputs(tmp_path: Path) -> None:
    pytest.importorskip("pyfastx")
    fasta = _write_text(tmp_path, "reads.fa", FASTA_TEXT)
    with pytest.raises(ConfigurationError, match="could not read"):
        detect_fastq_quality_encoding(fasta, verbose=False)
    fastq = _write_text(tmp_path, "reads.fq", FASTQ_TEXT)
    with pytest.raises(ConfigurationError, match="max_bytes=10"):
        detect_fastq_quality_encoding(fastq, max_bytes=10, verbose=False)


def test_pyfastx_missing_package_names_the_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_text(tmp_path, "reads.fq", FASTQ_TEXT)
    monkeypatch.setitem(sys.modules, "pyfastx", None)
    with pytest.raises(BackendUnavailableError, match="bio-pyfastx"):
        detect_fastq_quality_encoding(path, verbose=False)


# screed --------------------------------------------------------------------------------


def test_screed_database_round_trips_fasta_records(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write_text(tmp_path, "reads.fa", FASTA_TEXT)
    built = build_screed_database(source, database=tmp_path / "reads", verbose=False)

    assert built["format"] == "fasta"
    assert built["record_count"] == 2
    assert built["database"].endswith("reads_screed")

    fetched = fetch_screed_records(built["database"], names=["r2", "r1"], verbose=False)
    records = fetched["records"]
    assert [record["name"] for record in records] == ["r2", "r1"]
    assert records[0] == {"name": "r2", "description": "", "length": 4, "sequence": "GGCC"}
    assert records[1] == {
        "name": "r1",
        "description": "first read",
        "length": 9,
        "sequence": "ACGTNACGT",
    }


def test_screed_database_keeps_fastq_quality(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write_text(tmp_path, "reads.fq", FASTQ_TEXT)
    built = build_screed_database(source, database=tmp_path / "reads_fq", verbose=False)

    assert built["format"] == "fastq"
    (record,) = fetch_screed_records(built["database"], names=["q1"], verbose=False)["records"]
    assert record["sequence"] == "ACGTAC"
    assert record["quality"] == "IIII#!"
    assert record["description"] == "desc"


def test_screed_build_refuses_existing_database_and_removes_partial_output(
    tmp_path: Path,
) -> None:
    pytest.importorskip("screed")
    source = _write_text(tmp_path, "reads.fa", FASTA_TEXT)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        build_screed_database(source, database=tmp_path / "limited", max_records=1, verbose=False)
    assert not (tmp_path / "limited_screed").exists()

    built = build_screed_database(source, database=tmp_path / "once", verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        build_screed_database(source, database=tmp_path / "once", verbose=False)
    assert built["record_count"] == 2


def test_screed_fetch_enforces_limits_and_missing_names(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write_text(tmp_path, "reads.fa", FASTA_TEXT)
    built = build_screed_database(source, database=tmp_path / "limits", verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        fetch_screed_records(built["database"], max_records=1, verbose=False)
    with pytest.raises(ConfigurationError, match="no records named"):
        fetch_screed_records(built["database"], names=["r1", "missing"], verbose=False)


def test_screed_missing_package_names_the_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _write_text(tmp_path, "reads.fa", FASTA_TEXT)
    monkeypatch.setitem(sys.modules, "screed", None)
    with pytest.raises(BackendUnavailableError, match="bio-screed"):
        build_screed_database(source, database=tmp_path / "unused", verbose=False)


# py2bit --------------------------------------------------------------------------------


def test_summarize_twobit_counts_bases_and_masks(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    path = _write_twobit(tmp_path / "toy.2bit", TWOBIT_RECORDS)

    result = summarize_twobit(path, verbose=False)

    chr1, chr2 = result["chromosomes"]
    # chr1 is NNNNACGTacgtGGCCN: 17 bases, 5 N, and 12 ACGT with G+C = 8 (66.7%).
    assert (chr1["name"], chr1["length"]) == ("chr1", 17)
    assert (chr1["hard_masked_bases"], chr1["acgt_bases"]) == (5, 12)
    assert chr1["gc_percent_of_acgt"] == pytest.approx(100 * 8 / 12)
    # chr2 is TTTTAAAACCCCGGGG: 16 ACGT bases and G+C = 8 (50%).
    assert (chr2["length"], chr2["hard_masked_bases"], chr2["acgt_bases"]) == (16, 0, 16)
    assert chr2["gc_percent_of_acgt"] == pytest.approx(50.0)
    assert result["total_length"] == 33
    assert result["total_hard_masked"] == 5


def test_fetch_twobit_region_is_zero_based_and_half_open(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    path = _write_twobit(tmp_path / "toy.2bit", TWOBIT_RECORDS)

    chr1 = fetch_twobit_region(path, "chr1", 4, 12, verbose=False)
    # Positions 4 to 11 of NNNNACGTacgtGGCCN are ACGTacgt; soft-masking is returned uppercase.
    assert chr1["sequence"] == "ACGTACGT"
    assert (chr1["start"], chr1["end"], chr1["length"]) == (4, 12, 8)
    assert chr1["hard_masked_bases"] == 0

    chr2 = fetch_twobit_region(path, "chr2", 2, 6, verbose=False)
    assert chr2["sequence"] == "TTAA"


def test_fetch_twobit_region_rejects_invalid_intervals(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    path = _write_twobit(tmp_path / "toy.2bit", TWOBIT_RECORDS)
    with pytest.raises(ConfigurationError, match="is not present"):
        fetch_twobit_region(path, "chrX", 0, 4, verbose=False)
    with pytest.raises(ConfigurationError, match="beyond"):
        fetch_twobit_region(path, "chr2", 0, 17, verbose=False)
    with pytest.raises(ConfigurationError, match="max_length=3"):
        fetch_twobit_region(path, "chr1", 0, 8, max_length=3, verbose=False)
    with pytest.raises(ConfigurationError, match="end"):
        fetch_twobit_region(path, "chr1", 5, 5, verbose=False)


def test_py2bit_rejects_a_file_that_is_not_2bit(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    path = _write_text(tmp_path, "not_twobit.bin", "this is not a 2bit file\n")
    with pytest.raises(ConfigurationError, match="could not open"):
        summarize_twobit(path, verbose=False)


def test_py2bit_missing_package_names_the_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_text(tmp_path, "toy.2bit", "")
    monkeypatch.setitem(sys.modules, "py2bit", None)
    with pytest.raises(BackendUnavailableError, match="bio-py2bit"):
        summarize_twobit(path, verbose=False)
