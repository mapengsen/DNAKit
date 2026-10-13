"""Family group H3: Biopython, scikit-bio, screed, cutadapt, umi_tools and khmer.

The first half tests the ten new Biopython wrappers, ``correct_jc69_distance`` and
``query_saved_khmer_table`` against the real upstream call or a hand calculation. The second half
checks, for the ``family-variant`` rows that were kept, that the claimed parameter of an existing
wrapper reaches the variant, and records why the other claims were changed (same object, deprecated,
broken upstream). Every input is built in ``tmp_path``; nothing needs the network.
"""

from __future__ import annotations

import inspect
import io
import math
import random
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import khmer_family_functions as khmer_family
from dnakit.sequence import Bio_family_functions as bio_family
from dnakit.sequence import skbio_family_functions as skbio_family

WRAPPERS = [
    bio_family.bootstrap_replicate_trees,
    bio_family.parse_primer3_records,
    bio_family.read_alignment_blocks,
    bio_family.read_alignment_file,
    bio_family.read_blast_records,
    bio_family.read_genbank_record_fields,
    bio_family.read_goa_gene_product_file,
    bio_family.read_phd_records,
    bio_family.read_phylogeny_collection,
    bio_family.read_substitution_matrix_file,
    skbio_family.correct_jc69_distance,
    khmer_family.query_saved_khmer_table,
]

MAF_TEXT = (
    "##maf version=1 scoring=none\n"
    "a score=23262.0\n"
    "s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n"
    "s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG\n"
    "s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG\n"
    "\n"
    "a score=5062.0\n"
    "s hg18.chr7    27699739 6 + 158545518 TAAAGA\n"
    "s panTro1.chr6 28862317 6 + 161576975 TAAAGA\n"
    "s baboon         241163 6 +   4622798 TAAAGA\n"
    "\n"
)

M10_TEXT = """FASTA searches a protein or DNA sequence data bank
 version 36.3.8 Aug, 2015

1>>>query_seq - 20 nt
 vs  demo_lib library
>>>query_seq, 20 nt vs demo_lib library
; pg_name: FASTA
; pg_ver: 36.3.8
; pg_matrix: +5/-4
>>sbj1 demo subject one (22 nt)
; fa_frame: f
; fa_opt: 100
; fa_expect:  0.001
; sw_score: 100
; sw_overlap: 20
>query_seq ..
; sq_len: 20
; sq_offset: 1
; sq_type: D
; al_start: 1
; al_stop: 20
; al_display_start: 1
ACGTACGGTTAACCGGATAC
>sbj1 ..
; sq_len: 22
; sq_offset: 1
; sq_type: D
; al_start: 3
; al_stop: 22
; al_display_start: 3
ACGTACGGTTAACCGGATAC
>--
; fa_frame: f
; fa_opt: 80
; fa_expect:  0.01
; sw_score: 80
; sw_overlap: 10
>query_seq ..
; sq_len: 20
; sq_offset: 1
; sq_type: D
; al_start: 1
; al_stop: 10
; al_display_start: 1
ACGTACGGTT
>sbj1 ..
; sq_len: 22
; sq_offset: 1
; sq_type: D
; al_start: 11
; al_stop: 20
; al_display_start: 11
ACGTACGGTT
>>><<<
"""

BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.13.0+</BlastOutput_version>
  <BlastOutput_reference>Reference text</BlastOutput_reference>
  <BlastOutput_db>demo_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query_one demo query</BlastOutput_query-def>
  <BlastOutput_query-len>20</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
      <Parameters_gap-open>5</Parameters_gap-open>
      <Parameters_gap-extend>2</Parameters_gap-extend>
      <Parameters_filter>L;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>query_one demo query</Iteration_query-def>
      <Iteration_query-len>20</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>subject_one demo subject</Hit_def>
          <Hit_accession>0</Hit_accession>
          <Hit_len>60</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>37.4</Hsp_bit-score>
              <Hsp_score>20</Hsp_score>
              <Hsp_evalue>1.2e-05</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>11</Hsp_hit-from>
              <Hsp_hit-to>30</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_positive>20</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_qseq>ACGTACGGTTAACCGGATAC</Hsp_qseq>
              <Hsp_hseq>ACGTACGGTTAACCGGATAC</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>2</Statistics_db-num>
          <Statistics_db-len>120</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.46</Statistics_kappa>
          <Statistics_lambda>1.28</Statistics_lambda>
          <Statistics_entropy>0.85</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

BLAST_XML2 = """<?xml version="1.0"?>
<BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd">
<BlastOutput2>
<report>
<Report>
<program>blastn</program>
<version>BLASTN 2.13.0+</version>
<reference>Reference text</reference>
<search-target>
<Target>
<db>demo_db</db>
</Target>
</search-target>
<params>
<Parameters>
<expect>10</expect>
<sc-match>1</sc-match>
<sc-mismatch>-3</sc-mismatch>
<gap-open>5</gap-open>
<gap-extend>2</gap-extend>
<filter>L;m;</filter>
</Parameters>
</params>
<results>
<Results>
<search>
<Search>
<query-id>Query_1</query-id>
<query-title>query_one demo query</query-title>
<query-len>20</query-len>
<hits>
<Hit>
<num>1</num>
<description>
<HitDescr>
<id>gnl|BL_ORD_ID|0</id>
<accession>0</accession>
<title>subject_one demo subject</title>
</HitDescr>
</description>
<len>60</len>
<hsps>
<Hsp>
<num>1</num>
<bit-score>37.4</bit-score>
<score>20</score>
<evalue>1.2e-05</evalue>
<identity>20</identity>
<positive>20</positive>
<query-from>1</query-from>
<query-to>20</query-to>
<hit-from>11</hit-from>
<hit-to>30</hit-to>
<query-strand>Plus</query-strand>
<hit-strand>Plus</hit-strand>
<align-len>20</align-len>
<gaps>0</gaps>
<qseq>ACGTACGGTTAACCGGATAC</qseq>
<hseq>ACGTACGGTTAACCGGATAC</hseq>
<midline>||||||||||||||||||||</midline>
</Hsp>
</hsps>
</Hit>
</hits>
<stat>
<Statistics>
<db-num>2</db-num>
<db-len>120</db-len>
<hsp-len>0</hsp-len>
<eff-space>0</eff-space>
<kappa>0.46</kappa>
<lambda>1.28</lambda>
<entropy>0.85</entropy>
</Statistics>
</stat>
</Search>
</search>
</Results>
</results>
</Report>
</report>
</BlastOutput2>
</BlastXML2>
"""  # noqa: E501

PRIMER3_TEXT = (
    "# EPRIMER3 RELEASE 5.0.0\n"
    "#                      Start  Len   Tm     GC%   Sequence\n"
    "\n"
    "   1 PRODUCT SIZE: 60\n"
    "     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT\n"
    "     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA\n"
    "\n"
    "# EPRIMER3 RELEASE 5.0.0\n"
    "#                      Start  Len   Tm     GC%   Sequence\n"
    "\n"
    "   1 PRODUCT SIZE: 80\n"
    "     FORWARD PRIMER    5  18  59.50  55.00  GGCCTTAAGGCCTTAAGG\n"
    "     REVERSE PRIMER   80  18  59.80  50.00  CCGGAATTCCGGAATTCC\n"
)

GENBANK_TEXT = """LOCUS       SEQ1                      18 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  Demo sequence one.
ACCESSION   SEQ1
VERSION     SEQ1.1
KEYWORDS    demo; test.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria.
FEATURES             Location/Qualifiers
     source          1..18
                     /organism="Escherichia coli"
     gene            1..12
                     /gene="abcA"
ORIGIN
        1 atgaaacgta ttgcgtaa
//
LOCUS       SEQ2                      12 bp    DNA     linear   BCT 02-JAN-2020
DEFINITION  Demo sequence two.
ACCESSION   SEQ2
VERSION     SEQ2.1
KEYWORDS    .
SOURCE      Bacillus subtilis
  ORGANISM  Bacillus subtilis
            Bacteria; Firmicutes; Bacilli.
FEATURES             Location/Qualifiers
     CDS             1..12
                     /gene="xyzB"
ORIGIN
        1 atggcgtacg tt
//
"""

PHD_TEXT = (
    "BEGIN_SEQUENCE read1\nBEGIN_COMMENT\nCHROMAT_FILE: read1.ab1\n"
    "PHRED_VERSION: 0.980528.c\nEND_COMMENT\nBEGIN_DNA\nA 10 1\nC 20 4\nEND_DNA\nEND_SEQUENCE\n\n"
    "BEGIN_SEQUENCE read2\nBEGIN_COMMENT\nCHROMAT_FILE: read2.ab1\n"
    "PHRED_VERSION: 0.980528.c\nEND_COMMENT\nBEGIN_DNA\nG 30 2\nT 40 3\nA 50 5\nEND_DNA\n"
    "END_SEQUENCE\n"
)

NEWICK_TREES = (
    "((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);\n((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06);\n"
)

NEXUS_TREES = (
    "#NEXUS\nBEGIN TAXA;\nDIMENSIONS NTAX=4;\nTAXLABELS A B C D;\nEND;\nBEGIN TREES;\n"
    "TREE t1 = ((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);\n"
    "TREE t2 = ((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06);\nEND;\n"
)

ALIGNED = {
    "a": "ACGTACGTAC",
    "b": "ACGTACGAAC",
    "c": "ACGAACGTTC",
    "d": "TCGTACGTAC",
}


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# --- conventions ------------------------------------------------------------------------------


def test_every_wrapper_has_native_keyword_only_verbose() -> None:
    for function in WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_modules_export_only_new_names_and_no_collisions() -> None:
    import dnakit.microbiome as microbiome
    import dnakit.sequence as sequence

    names = list(bio_family.__all__) + list(skbio_family.__all__) + list(khmer_family.__all__)
    assert len(names) == len(set(names)) == len(WRAPPERS)
    for name in bio_family.__all__ + skbio_family.__all__:
        assert name not in sequence.__all__ or inspect.unwrap(
            getattr(sequence, name)
        ) is inspect.unwrap(getattr(bio_family, name, None) or getattr(skbio_family, name))
    for name in khmer_family.__all__:
        assert name in microbiome.__all__
        assert inspect.unwrap(getattr(microbiome, name)) is inspect.unwrap(
            getattr(khmer_family, name)
        )


def test_public_signatures_have_no_integer_keyed_mappings() -> None:
    for function in WRAPPERS:
        for parameter in inspect.signature(function).parameters.values():
            assert "Mapping[int" not in str(parameter.annotation), function.__name__


# --- read_alignment_file ----------------------------------------------------------------------


@pytest.mark.parametrize("file_format", ["bed", "clustal", "fasta", "maf", "phylip", "psl", "sam"])
def test_read_alignment_file_matches_bio_align_parse(tmp_path: Path, file_format: str) -> None:
    from Bio import Align

    from dnakit.sequence.Bio_review_functions import write_pairwise_alignment_file

    path = tmp_path / f"pair.{file_format}"
    write_pairwise_alignment_file(
        "ACGTACGGTTAACCGGATAC",
        "TTACGTACGTTAACCGATACGG",
        path,
        file_format=file_format,
        verbose=False,
    )
    result = bio_family.read_alignment_file(path, file_format=file_format, verbose=False)
    upstream = list(Align.parse(str(path), file_format))
    assert result["alignment_count"] == len(upstream) == 1
    item = result["alignments"][0]
    assert item["coordinates"] == upstream[0].coordinates.tolist()
    assert item["sequence_ids"] == ["query", "target"]
    assert item["aligned_columns"] == upstream[0].length
    if file_format in {"clustal", "fasta", "maf", "phylip"}:
        assert item["aligned_rows"] == ["--ACGTACGGTTAACCGGATAC--", "TTACGTAC-GTTAACC-GATACGG"]
    else:
        assert item["aligned_rows"] is None
    json_ready = __import__("json").dumps(result)
    assert "NaN" not in json_ready


def test_read_alignment_file_exactly_one_follows_bio_align_read(tmp_path: Path) -> None:
    from Bio import Align

    path = tmp_path / "pair.psl"
    aligner = Align.PairwiseAligner()
    first = aligner.align("ACGTACGT", "ACGTACGT")[0]
    second = aligner.align("ACGTAC", "ACGTAC")[0]
    with path.open("w", encoding="utf-8") as handle:
        Align.write([first, second], handle, "psl")
    assert bio_family.read_alignment_file(path, verbose=False)["alignment_count"] == 2
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_file(path, exactly_one=True, verbose=False)
    with pytest.raises(Exception, match=r"than one"):
        Align.read(str(path), "psl")
    single = tmp_path / "single.psl"
    with single.open("w", encoding="utf-8") as handle:
        Align.write([first], handle, "psl")
    one = bio_family.read_alignment_file(single, exactly_one=True, verbose=False)
    assert one["alignment_count"] == 1
    assert one["provenance"]["parser"] == "Bio.Align.read"


def test_read_alignment_file_validates_inputs(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_file(tmp_path / "missing.psl", verbose=False)
    path = _write(tmp_path / "x.psl", "not a psl line\n")
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_file(path, verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_file(path, file_format="exonerate", verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_file(path, exactly_one="yes", verbose=False)  # type: ignore[arg-type]


# --- read_alignment_blocks --------------------------------------------------------------------


def test_read_alignment_blocks_reads_every_maf_block(tmp_path: Path) -> None:
    from Bio import AlignIO

    from dnakit.bioinformatics.alignment import read_alignment

    path = _write(tmp_path / "blocks.maf", MAF_TEXT)
    result = bio_family.read_alignment_blocks(path, verbose=False)
    upstream = list(AlignIO.parse(str(path), "maf"))
    assert result["block_count"] == len(upstream) == 2
    for block, reference in zip(result["blocks"], upstream, strict=True):
        assert [row["sequence"] for row in block["rows"]] == [str(r.seq) for r in reference]
        assert [row["id"] for row in block["rows"]] == [r.id for r in reference]
    first = result["blocks"][0]["rows"][0]["annotations"]
    assert first == {"size": 38, "srcSize": 158545518, "start": 27578828, "strand": 1}
    # The single-alignment reader of dnakit rejects a file with several blocks.
    with pytest.raises(ValueError, match="More than one"):
        read_alignment(path, format="maf", verbose=False)
    with pytest.raises(ConfigurationError, match="max_blocks=1"):
        bio_family.read_alignment_blocks(path, max_blocks=1, verbose=False)


def test_read_alignment_blocks_reads_fasta_m10_hits(tmp_path: Path) -> None:
    from Bio import AlignIO

    path = _write(tmp_path / "hits.m10", M10_TEXT)
    result = bio_family.read_alignment_blocks(path, file_format="fasta-m10", verbose=False)
    upstream = list(AlignIO.parse(str(path), "fasta-m10"))
    assert result["block_count"] == len(upstream) == 2
    assert [block["alignment_length"] for block in result["blocks"]] == [20, 10]
    assert result["blocks"][0]["annotations"]["fa_expect"] == "0.001"
    assert result["blocks"][1]["rows"][1]["sequence"] == "ACGTACGGTT"
    assert result["blocks"][0]["rows"][1]["annotations"]["original_length"] == 22
    broken = _write(tmp_path / "broken.m10", "garbage\n")
    assert (
        bio_family.read_alignment_blocks(broken, file_format="fasta-m10", verbose=False)[
            "block_count"
        ]
        == 0
    )
    with pytest.raises(ConfigurationError):
        bio_family.read_alignment_blocks(path, file_format="clustal", verbose=False)


# --- read_substitution_matrix_file ------------------------------------------------------------


def test_read_substitution_matrix_file_matches_upstream(tmp_path: Path) -> None:
    from Bio.Align import substitution_matrices

    text = (
        "# demo nucleotide matrix\n   A  C  G  T\n"
        "A  2 -3 -3 -3\nC -3  2 -3 -3\nG -3 -3  2 -3\nT -3 -3 -3  2\n"
    )
    path = _write(tmp_path / "nuc.mat", text)
    result = bio_family.read_substitution_matrix_file(path, verbose=False)
    upstream = substitution_matrices.read(io.StringIO(text))
    assert result["alphabet"] == str(upstream.alphabet) == "ACGT"
    assert result["matrix"][0] == [2.0, -3.0, -3.0, -3.0]
    assert result["matrix"][2][1] == float(upstream[("G", "C")])
    assert result["symmetric"] is True
    assert result["comments"] == ["demo nucleotide matrix"]
    asymmetric = _write(tmp_path / "asym.mat", "   A  C\nA  1  2\nC  3  4\n")
    assert bio_family.read_substitution_matrix_file(asymmetric, verbose=False)["symmetric"] is False


def test_read_substitution_matrix_file_rejects_protein_alphabets(tmp_path: Path) -> None:
    protein = _write(tmp_path / "prot.mat", "   L  E\nL  4 -3\nE -3  5\n")
    with pytest.raises(ConfigurationError, match="non-nucleotide"):
        bio_family.read_substitution_matrix_file(protein, verbose=False)
    bad = _write(tmp_path / "bad.mat", "   A  C\nA  1\n")
    with pytest.raises(ConfigurationError):
        bio_family.read_substitution_matrix_file(bad, verbose=False)


# --- read_blast_records -----------------------------------------------------------------------


def test_read_blast_records_matches_bio_blast_for_xml_and_xml2(tmp_path: Path) -> None:
    from Bio import Blast, SearchIO

    xml = _write(tmp_path / "blast.xml", BLAST_XML)
    xml2 = _write(tmp_path / "blast2.xml", BLAST_XML2)
    for path in (xml, xml2):
        result = bio_family.read_blast_records(path, verbose=False)
        record = Blast.read(str(path))
        row = result["records"][0]
        assert result["record_count"] == 1
        assert row["query_id"] == record.query.id == "Query_1"
        assert row["query_length"] == 20
        assert row["program"] == "blastn"
        assert row["database"] == "demo_db"
        hit = row["hits"][0]
        assert hit["target_description"] == "subject_one demo subject"
        assert hit["target_length"] == 60
        hsp = hit["hsps"][0]
        assert hsp["evalue"] == pytest.approx(1.2e-05)
        assert hsp["bit_score"] == pytest.approx(37.4)
        assert hsp["target_coordinates"] == [10, 30]
        assert hsp["query_coordinates"] == [0, 20]
        assert hsp["identity"] == 20
    # The SearchIO parser behind parse_search_results returns nothing for the XML2 layout.
    assert len(list(SearchIO.parse(str(xml), "blast-xml"))) == 1
    assert list(SearchIO.parse(str(xml2), "blast-xml")) == []


def test_read_blast_records_exactly_one_and_errors(tmp_path: Path) -> None:
    xml = _write(tmp_path / "blast.xml", BLAST_XML)
    one = bio_family.read_blast_records(xml, exactly_one=True, verbose=False)
    assert one["provenance"]["parser"] == "Bio.Blast.read"
    assert one["records"][0]["hit_count"] == 1
    with pytest.raises(ConfigurationError, match="max_hits"):
        bio_family.read_blast_records(xml, max_hits=0, verbose=False)
    not_xml = _write(tmp_path / "not.xml", "this is not xml\n")
    with pytest.raises(ConfigurationError):
        bio_family.read_blast_records(not_xml, verbose=False)


# --- parse_primer3_records --------------------------------------------------------------------


def test_parse_primer3_records_reads_every_record(tmp_path: Path) -> None:
    from Bio.Emboss import Primer3

    from dnakit.sequence.Bio_functions import parse_primer3_output

    path = _write(tmp_path / "multi.txt", PRIMER3_TEXT)
    result = bio_family.parse_primer3_records(path, verbose=False)
    with path.open(encoding="utf-8") as handle:
        upstream = list(Primer3.parse(handle))
    assert result["record_count"] == len(upstream) == 2
    assert result["pair_count"] == 2
    second = result["records"][1]["primer_pairs"][0]
    assert second["product_size"] == 80
    assert second["forward"] == {
        "sequence": "GGCCTTAAGGCCTTAAGG",
        "start": 5,
        "length": 18,
        "tm": 59.5,
        "gc_percent": 55.0,
    }
    assert second["reverse"]["start"] == 80
    assert second["internal"] is None
    # parse_primer3_output reads one record only.
    with pytest.raises(ConfigurationError, match="single Primer3 record"):
        parse_primer3_output(path, verbose=False)
    single = _write(
        tmp_path / "single.txt", PRIMER3_TEXT.split("# EPRIMER3", 2)[1].join(["# EPRIMER3", ""])
    )
    assert parse_primer3_output(single, verbose=False)["pair_count"] == 1


def test_parse_primer3_records_refuses_empty_files_and_limits(tmp_path: Path) -> None:
    empty = _write(tmp_path / "empty.txt", "\n\n")
    with pytest.raises(ConfigurationError, match="empty"):
        bio_family.parse_primer3_records(empty, verbose=False)
    path = _write(tmp_path / "multi.txt", PRIMER3_TEXT)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        bio_family.parse_primer3_records(path, max_records=1, verbose=False)
    headless = _write(tmp_path / "headless.txt", "FORWARD PRIMER 1\n")
    with pytest.raises(ConfigurationError):
        bio_family.parse_primer3_records(headless, verbose=False)


# --- read_genbank_record_fields ---------------------------------------------------------------


def test_read_genbank_record_fields_matches_bio_genbank(tmp_path: Path) -> None:
    from Bio import GenBank

    from dnakit.sequence.backends import read_sequence_records

    path = _write(tmp_path / "two.gb", GENBANK_TEXT)
    result = bio_family.read_genbank_record_fields(path, verbose=False)
    with path.open(encoding="utf-8") as handle:
        upstream = list(GenBank.parse(handle))
    assert result["record_count"] == len(upstream) == 2
    first = result["records"][0]
    assert (first["locus"], first["size"], first["data_file_division"]) == ("SEQ1", 18, "BCT")
    assert first["date"] == "01-JAN-2020"
    assert first["version"] == "SEQ1.1"
    assert first["keywords"] == ["demo", "test"]
    assert first["organism"] == "Escherichia coli"
    assert first["taxonomy"] == ["Bacteria", "Proteobacteria", "Gammaproteobacteria"]
    assert first["sequence"] == str(upstream[0].sequence).upper()
    assert [feature["type"] for feature in first["features"]] == ["source", "gene"]
    assert first["features"][1]["location"] == "1..12"
    assert first["features"][1]["qualifiers"] == [{"key": "gene", "value": '"abcA"'}]
    assert result["records"][1]["keywords"] == []
    assert result["total_bases"] == 30
    # read_sequence_records keeps none of the header fields.
    kept = read_sequence_records(path, format="genbank", verbose=False)["records"][0]
    assert "organism" not in kept and "taxonomy" not in kept and "date" not in kept


def test_read_genbank_record_fields_exactly_one_follows_bio_genbank_read(tmp_path: Path) -> None:
    path = _write(tmp_path / "two.gb", GENBANK_TEXT)
    with pytest.raises(ConfigurationError, match="More than one"):
        bio_family.read_genbank_record_fields(path, exactly_one=True, verbose=False)
    single = _write(tmp_path / "one.gb", GENBANK_TEXT.split("//\n")[0] + "//\n")
    result = bio_family.read_genbank_record_fields(single, exactly_one=True, verbose=False)
    assert result["record_count"] == 1
    assert result["provenance"]["parser"] == "Bio.GenBank.read"
    with pytest.raises(ConfigurationError, match="max_bases"):
        bio_family.read_genbank_record_fields(path, max_bases=20, verbose=False)


# --- read_phylogeny_collection ----------------------------------------------------------------


@pytest.mark.parametrize("file_format", ["newick", "nexus", "phyloxml", "nexml", "cdao"])
def test_read_phylogeny_collection_reads_every_tree(tmp_path: Path, file_format: str) -> None:
    from Bio import Phylo

    from dnakit.bioinformatics.phylogeny import read_tree

    trees = list(Phylo.parse(io.StringIO(NEWICK_TREES), "newick"))
    path = tmp_path / f"trees.{file_format}"
    if file_format == "nexus":
        _write(path, NEXUS_TREES)
    else:
        Phylo.write(trees, str(path), file_format)
    result = bio_family.read_phylogeny_collection(path, file_format=file_format, verbose=False)
    upstream = list(Phylo.parse(str(path), file_format))
    assert result["tree_count"] == len(upstream) == 2
    for item, reference in zip(result["trees"], upstream, strict=True):
        assert item["tips"] == ["A", "B", "C", "D"]
        assert item["newick"] == reference.format("newick").strip()
        assert item["total_branch_length"] == pytest.approx(reference.total_branch_length())
    if file_format in {"newick", "nexus", "phyloxml"}:
        assert result["trees"][0]["total_branch_length"] == pytest.approx(1.11)
    # read_tree takes a single tree and raises on a replicate set.
    with pytest.raises(ValueError, match="multiple trees"):
        read_tree(path, format=file_format, verbose=False)


def test_read_phylogeny_collection_names_and_limits(tmp_path: Path) -> None:
    path = _write(tmp_path / "trees.nex", NEXUS_TREES)
    result = bio_family.read_phylogeny_collection(path, file_format="nexus", verbose=False)
    assert [tree["name"] for tree in result["trees"]] == ["t1", "t2"]
    with pytest.raises(ConfigurationError, match="max_trees=1"):
        bio_family.read_phylogeny_collection(path, file_format="nexus", max_trees=1, verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.read_phylogeny_collection(path, file_format="mrbayes", verbose=False)
    garbage = _write(tmp_path / "garbage.xml", "<not-a-tree/>")
    with pytest.raises(ConfigurationError):
        bio_family.read_phylogeny_collection(garbage, file_format="phyloxml", verbose=False)


# --- bootstrap_replicate_trees ----------------------------------------------------------------


def test_bootstrap_replicate_trees_matches_upstream_and_restores_random_state() -> None:
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.Consensus import bootstrap_trees
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    random.seed(99)
    expected_next = random.Random(99).random()
    result = bio_family.bootstrap_replicate_trees(ALIGNED, replicates=4, seed=7, verbose=False)
    assert random.random() == expected_next  # the caller's random state is untouched
    alignment = MultipleSeqAlignment([SeqRecord(Seq(v), id=k) for k, v in ALIGNED.items()])
    constructor = DistanceTreeConstructor(DistanceCalculator("identity"), "nj")
    state = random.getstate()
    random.seed(7)
    try:
        direct = [t.format("newick").strip() for t in bootstrap_trees(alignment, 4, constructor)]
    finally:
        random.setstate(state)
    assert result["trees"] == direct
    assert result["ids"] == ["a", "b", "c", "d"]
    again = bio_family.bootstrap_replicate_trees(ALIGNED, replicates=4, seed=7, verbose=False)
    assert again["trees"] == result["trees"]
    other = bio_family.bootstrap_replicate_trees(ALIGNED, replicates=4, seed=8, verbose=False)
    assert other["trees"] != result["trees"]


def test_bootstrap_replicate_trees_validates_inputs() -> None:
    with pytest.raises(ConfigurationError):
        bio_family.bootstrap_replicate_trees({"a": "ACGT", "b": "ACGT"}, verbose=False)
    with pytest.raises(ConfigurationError, match="same length"):
        bio_family.bootstrap_replicate_trees({"a": "ACGT", "b": "ACG", "c": "ACGT"}, verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.bootstrap_replicate_trees(ALIGNED, replicates=0, verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.bootstrap_replicate_trees(ALIGNED, method="ml", verbose=False)
    with pytest.raises(ConfigurationError):
        bio_family.bootstrap_replicate_trees({**ALIGNED, "e": "ACGN" * 2 + "AC"}, verbose=False)


def test_bootstrap_replicate_trees_feed_branch_support() -> None:
    from dnakit.sequence.Bio_review_functions import bootstrap_branch_support

    result = bio_family.bootstrap_replicate_trees(ALIGNED, replicates=5, seed=3, verbose=False)
    target = "((a:0.1,b:0.1):0.1,c:0.1,d:0.1);"
    support = bootstrap_branch_support(target, result["trees"], verbose=False)
    assert support["status"] == "completed"


# --- read_phd_records -------------------------------------------------------------------------


def test_read_phd_records_matches_bio_phd_parse(tmp_path: Path) -> None:
    from Bio.Sequencing import Phd

    from dnakit.sequence.Bio_functions import read_phd_sequence

    path = _write(tmp_path / "reads.phd", PHD_TEXT)
    result = bio_family.read_phd_records(path, verbose=False)
    upstream = list(Phd.parse(str(path)))
    assert result["record_count"] == len(upstream) == 2
    second = result["records"][1]
    assert second["file_name"] == "read2"
    assert second["sequence"] == "GTA"
    assert second["qualities"] == [30, 40, 50]
    assert second["peak_locations"] == [2, 3, 5]
    assert second["comments"]["chromat_file"] == "read2.ab1"
    assert [site[0] for site in upstream[1].sites] == list(second["sequence"])
    with pytest.raises(ConfigurationError, match="PHD"):
        read_phd_sequence(path, verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        bio_family.read_phd_records(path, max_records=1, verbose=False)
    empty = _write(tmp_path / "empty.phd", "")
    assert bio_family.read_phd_records(empty, verbose=False)["record_count"] == 0


# --- read_goa_gene_product_file ---------------------------------------------------------------


GPA_TEXT = (
    "!gpa-version: 1.1\n"
    "UniProtKB\tP12345\tenables\tGO:0003700\tPMID:1234567\tECO:0000314\t\t\t20240101\tUniProt\t\t\n"
    "UniProtKB\tP67890\tinvolved_in\tGO:0006355\tPMID:7654321|PMID:1\tECO:0000501\tGO:1\t\t"
    "20240102\tUniProt\t\t\n"
)
GPI_TEXT = (
    "!gpi-version: 1.2\n"
    "UniProtKB\tP12345\tGENEA\tTranscription factor A\tGENEA|TFA\tprotein\ttaxon:9606\t\t\t\n"
)


def test_read_goa_gene_product_file_matches_upstream_iterators(tmp_path: Path) -> None:
    from Bio.UniProt import GOA

    gpa = _write(tmp_path / "prod.gpa", GPA_TEXT)
    result = bio_family.read_goa_gene_product_file(gpa, verbose=False)
    with gpa.open(encoding="utf-8") as handle:
        upstream = list(GOA.gpa_iterator(handle))
    assert result["version"] == "1.1"
    assert result["record_count"] == len(upstream) == 2
    second = result["records"][1]
    assert second["DB_Object_ID"] == "P67890"
    assert second["DB:Reference"] == ["PMID:7654321", "PMID:1"]
    assert second["With"] == ["GO:1"]
    assert second["ECO_Evidence_code"] == "ECO:0000501"
    assert result["records"][0]["With"] == []

    gpi = _write(tmp_path / "prod.gpi", GPI_TEXT)
    info = bio_family.read_goa_gene_product_file(gpi, file_format="gpi", verbose=False)
    assert info["version"] == "1.2"
    record = info["records"][0]
    assert record["DB_Object_Symbol"] == "GENEA"
    assert record["DB_Object_Synonym"] == ["GENEA", "TFA"]
    assert record["DB_Object_Type"] == "protein"
    assert record["Taxon"] == "taxon:9606"


def test_read_goa_gene_product_file_rejects_wrong_header(tmp_path: Path) -> None:
    gaf = _write(tmp_path / "x.gaf", "!gaf-version: 2.2\n")
    with pytest.raises(ConfigurationError):
        bio_family.read_goa_gene_product_file(gaf, verbose=False)
    gpi = _write(tmp_path / "x.gpi", GPI_TEXT)
    with pytest.raises(ConfigurationError):
        bio_family.read_goa_gene_product_file(gpi, file_format="gpa", verbose=False)
    two = _write(tmp_path / "two.gpa", GPA_TEXT)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        bio_family.read_goa_gene_product_file(two, max_records=1, verbose=False)


# --- correct_jc69_distance --------------------------------------------------------------------


def test_correct_jc69_distance_matches_closed_form_and_flags_saturation() -> None:
    result = skbio_family.correct_jc69_distance([0.0, 0.1, 0.5, 0.75, 0.9], verbose=False)
    # D = -3/4 ln(1 - 4/3 p)
    assert result["corrected"][0] == 0.0
    assert result["corrected"][1] == pytest.approx(-0.75 * math.log(1 - 0.1 / 0.75))
    assert result["corrected"][2] == pytest.approx(-0.75 * math.log(1 - 0.5 / 0.75))
    assert result["corrected"][3:] == [None, None]
    assert result["status"] == "partial"
    assert result["saturated_count"] == 2
    assert result["saturation_limit"] == 0.75
    protein = skbio_family.correct_jc69_distance([0.2], alphabet_size=20, verbose=False)
    assert protein["corrected"][0] == pytest.approx(-0.95 * math.log(1 - 0.2 / 0.95))
    assert protein["status"] == "completed"


def test_correct_jc69_distance_agrees_with_nucleotide_distance_model_jc69() -> None:
    from dnakit.sequence.skbio_functions import nucleotide_distance

    left, right = "ACGTACGTAAGTCCGTAGCTAGCATCGA", "ACGTTCGTACGTCCATAGCTAGGATCGT"
    p = nucleotide_distance(left, right, model="pdist", verbose=False)["distance"]
    direct = nucleotide_distance(left, right, model="jc69", verbose=False)["distance"]
    assert skbio_family.correct_jc69_distance([p], verbose=False)["corrected"][0] == pytest.approx(
        direct
    )


def test_correct_jc69_distance_validates_inputs() -> None:
    for bad in ([], "0.1", [-0.1], [1.5], [float("nan")], [True], ["0.1"]):
        with pytest.raises(ConfigurationError):
            skbio_family.correct_jc69_distance(bad, verbose=False)  # type: ignore[arg-type]
    for size in (1, 0, True, 2.5):
        with pytest.raises(ConfigurationError):
            skbio_family.correct_jc69_distance([0.1], alphabet_size=size, verbose=False)  # type: ignore[arg-type]


# --- query_saved_khmer_table ------------------------------------------------------------------

KHMER_SEQ = "ACGTACGGATCGATCGGCTAGCTAGGATCGATTTAGC"


def _save_tables(tmp_path: Path) -> dict[str, Path]:
    khmer = pytest.importorskip("khmer")
    paths = {}
    for name, cls in (
        ("counttable", khmer.Counttable),
        ("nodetable", khmer.Nodetable),
        ("nodegraph", khmer.Nodegraph),
        ("countgraph", khmer.Countgraph),
    ):
        table = cls(5, 1000, 4)
        table.consume(KHMER_SEQ)
        table.consume("ACGTA")
        paths[name] = tmp_path / name
        table.save(str(paths[name]))
    return paths


def test_query_saved_khmer_table_matches_khmer_loaders(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    paths = _save_tables(tmp_path)
    queries = ["ACGTA", "TACGT", "GGGGG", "acgta"]
    for name, loader in (
        ("counttable", khmer.load_counttable),
        ("nodetable", khmer.load_nodetable),
        ("nodegraph", khmer.load_nodegraph),
    ):
        result = bio = khmer_family.query_saved_khmer_table(
            paths[name], queries, table_type=name, verbose=False
        )
        table = loader(str(paths[name]))
        assert result["values"] == {
            kmer: int(table.get(kmer)) for kmer in sorted({q.upper() for q in queries})
        }
        assert bio["ksize"] == 5
    counts = khmer_family.query_saved_khmer_table(paths["counttable"], queries, verbose=False)
    assert counts["value_kind"] == "count"
    assert counts["values"]["ACGTA"] == 2  # added twice
    assert counts["values"]["TACGT"] == 2  # reverse complement shares the entry
    presence = khmer_family.query_saved_khmer_table(
        paths["nodegraph"], queries, table_type="nodegraph", verbose=False
    )
    assert presence["value_kind"] == "presence"
    assert presence["values"]["GGGGG"] == 0


def test_query_saved_countgraph_does_not_read_a_counttable(tmp_path: Path) -> None:
    """Why the load_counttable claim was wrong: the hash layouts differ and no error is raised."""

    from dnakit.microbiome.khmer_functions import query_saved_countgraph

    paths = _save_tables(tmp_path)
    wrong = query_saved_countgraph(paths["counttable"], ["ACGTA"], verbose=False)
    right = khmer_family.query_saved_khmer_table(paths["counttable"], ["ACGTA"], verbose=False)
    assert wrong["counts"]["ACGTA"] == 0
    assert right["values"]["ACGTA"] == 2
    with pytest.raises(ConfigurationError):
        query_saved_countgraph(paths["nodegraph"], ["ACGTA"], verbose=False)


def test_query_saved_khmer_table_validates_inputs(tmp_path: Path) -> None:
    paths = _save_tables(tmp_path)
    with pytest.raises(ConfigurationError, match="table_type"):
        khmer_family.query_saved_khmer_table(paths["counttable"], ["ACGTA"], table_type="x")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="not a khmer nodegraph"):
        khmer_family.query_saved_khmer_table(
            paths["counttable"], ["ACGTA"], table_type="nodegraph", verbose=False
        )
    with pytest.raises(ConfigurationError, match="ksize"):
        khmer_family.query_saved_khmer_table(paths["counttable"], ["ACGTAC"], verbose=False)
    with pytest.raises(ConfigurationError):
        khmer_family.query_saved_khmer_table(paths["counttable"], ["ACGNA"], verbose=False)
    with pytest.raises(ConfigurationError):
        khmer_family.query_saved_khmer_table(paths["counttable"], [], verbose=False)
    with pytest.raises(ConfigurationError):
        khmer_family.query_saved_khmer_table(tmp_path / "missing", ["ACGTA"], verbose=False)


# =============================== family-variant rows: the parameter reaches the variant ========


def test_write_tree_format_reaches_phylo_codec_writers(tmp_path: Path) -> None:
    from Bio.Phylo import CDAOIO, BaseTree, NeXMLIO, NexusIO, PhyloXMLIO

    from dnakit.bioinformatics.phylogeny import read_tree, write_tree

    source = _write(tmp_path / "t.nwk", "((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);\n")
    tree = read_tree(source, verbose=False)

    def clade(node):
        return BaseTree.Clade(
            branch_length=node.length, name=node.name, clades=[clade(c) for c in node.children]
        )

    reference = BaseTree.Tree(root=clade(tree))
    for fmt, module in (
        ("nexus", NexusIO),
        ("nexml", NeXMLIO),
        ("phyloxml", PhyloXMLIO),
        ("cdao", CDAOIO),
    ):
        written = tmp_path / f"w.{fmt}"
        write_tree(tree, written, format=fmt, verbose=False)
        buffer = io.StringIO()
        module.write([reference], buffer)
        assert written.read_text(encoding="utf-8") == buffer.getvalue(), fmt
        back = read_tree(written, format=fmt, verbose=False)
        assert sorted(back.tips) == ["A", "B", "C", "D"]


def test_read_alignment_format_nexus_reaches_nexus_iterator(tmp_path: Path) -> None:
    from Bio.AlignIO import NexusIO

    from dnakit.bioinformatics.alignment import read_alignment

    text = (
        "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=3 NCHAR=8;\nFORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
        "MATRIX\ns1 ACGTACGT\ns2 ACGTTCGT\ns3 ACG-ACGA\n;\nEND;\n"
    )
    path = _write(tmp_path / "aln.nex", text)
    with path.open(encoding="utf-8") as handle:
        blocks = list(NexusIO.NexusIterator(handle))
    assert len(blocks) == 1
    alignment = read_alignment(path, format="nexus", verbose=False)
    assert list(alignment.ids) == [r.id for r in blocks[0]]
    assert list(alignment.sequences) == [str(r.seq) for r in blocks[0]]


def test_parse_search_results_blast_xml_reaches_ncbixml_data(tmp_path: Path) -> None:
    from Bio.Blast import NCBIXML

    from dnakit.bioinformatics.alignment import parse_search_results

    path = _write(tmp_path / "blast.xml", BLAST_XML)
    with path.open(encoding="utf-8") as handle:
        legacy = list(NCBIXML.parse(handle))
    with path.open(encoding="utf-8") as handle:
        single = NCBIXML.read(handle)
    result = list(parse_search_results(path, format="blast-xml", verbose=False))
    assert len(legacy) == len(result) == 1
    hsp_legacy = legacy[0].alignments[0].hsps[0]
    hsp = result[0][0][0]
    assert legacy[0].query_length == result[0].seq_len == single.query_length == 20
    assert hsp_legacy.expect == hsp.evalue
    assert hsp_legacy.align_length == hsp.aln_span
    assert hsp_legacy.identities == hsp.ident_num
    assert hsp_legacy.query == str(hsp.query.seq)
    assert hsp_legacy.sbjct == str(hsp.hit.seq)
    assert legacy[0].alignments[0].length == result[0][0].seq_len


def test_convert_search_results_output_format_reaches_searchio_and_blast_writers(
    tmp_path: Path,
) -> None:
    from Bio import Blast, SearchIO

    from dnakit.sequence.Bio_review_functions import convert_search_results

    source = _write(tmp_path / "blast.xml", BLAST_XML)
    converted = tmp_path / "converted.xml"
    result = convert_search_results(source, "blast-xml", converted, "blast-xml", verbose=False)
    assert result["query_count"] == result["hsp_count"] == 1
    direct = tmp_path / "direct.xml"
    SearchIO.write(list(SearchIO.parse(str(source), "blast-xml")), str(direct), "blast-xml")
    assert converted.read_text(encoding="utf-8") == direct.read_text(encoding="utf-8")
    written = tmp_path / "bio_blast.xml"
    with Blast.parse(str(source)) as records:
        Blast.write(records, str(written), "XML")

    def summary(path: Path) -> list:
        with Blast.parse(str(path)) as records:
            return [
                (
                    r.query.id,
                    [
                        (h.target.id, [(s.score, s.annotations["evalue"], s.coordinates.tolist())])
                        for h in r
                        for s in h
                    ],
                )
                for r in records
            ]

    assert summary(written) == summary(converted)


def test_convert_quality_score_direction_reaches_solexa_quality_from_phred() -> None:
    from Bio.SeqIO.QualityIO import solexa_quality_from_phred

    from dnakit.sequence.Bio_functions import convert_quality_score

    for value in (0.0, 3.0, 10.0, 30.0, 40.5):
        wrapped = convert_quality_score(value, direction="phred_to_solexa", verbose=False)
        assert wrapped["output"] == pytest.approx(solexa_quality_from_phred(value))


def test_sequence_checksum_method_reaches_crc32_crc64_and_gcg() -> None:
    from Bio.SeqUtils import CheckSum

    from dnakit.sequence.Bio_functions import sequence_checksum

    text = "ACGTACGGTTAACC"
    for method in ("crc32", "crc64", "gcg"):
        wrapped = sequence_checksum(text, method=method, verbose=False)
        assert wrapped["value"] == getattr(CheckSum, method)(text)


def test_local_composition_complexity_full_window_reaches_lcc_simp() -> None:
    from Bio.SeqUtils import lcc

    from dnakit.sequence.Bio_functions import local_composition_complexity

    for text in ("ACGTACGTAAGTCCGTAGC", "AAAAAAAAAA", "ACGGTTTACAG"):
        wrapped = local_composition_complexity(text, window=len(text), verbose=False)
        assert wrapped["values"] == [pytest.approx(lcc.lcc_simp(text))]


def test_mcdonald_kreitman_test_species_grouping_reaches_codonalign_mktest() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio.codonalign import mktest as codonalign_mktest
        from Bio.codonalign.codonalignment import CodonAlignment
        from Bio.codonalign.codonalignment import mktest as module_mktest
        from Bio.codonalign.codonseq import CodonSeq
        from Bio.SeqRecord import SeqRecord

    from dnakit.sequence.Bio_functions import mcdonald_kreitman_test

    assert codonalign_mktest is module_mktest
    # Every codon column holds at most two distinct codons, so the upstream count does not depend
    # on set iteration order (with three or more codons per column both mktest functions do).
    ancestor = ["GCT", "GCT", "CTG", "CTG", "AAA", "AAA", "TTT", "GGG"]
    diverged = ["GCC", "GGT", "CTG", "CTG", "AAA", "AAA", "TTT", "GGG"]
    sequences: dict[str, str] = {}
    species: dict[str, str] = {}
    for index in range(4):
        within = list(ancestor)
        if index == 1:
            within[2] = "CTC"  # synonymous polymorphism
        if index == 2:
            within[3] = "GTG"  # nonsynonymous polymorphism
        sequences[f"A{index}"] = "".join(within)
        species[f"A{index}"] = "A"
        sequences[f"B{index}"] = "".join(diverged)
        species[f"B{index}"] = "B"
    alignments = [
        CodonAlignment(
            [SeqRecord(CodonSeq(sequences[k]), id=k) for k in sorted(sequences) if species[k] == s]
        )
        for s in ("A", "B")
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = codonalign_mktest(alignments)
    wrapped = mcdonald_kreitman_test(sequences, species, verbose=False)
    assert wrapped["p_value"] == pytest.approx(direct)


MOTIF_FILES = {
    "alignace": (
        "alignace.out",
        "AlignACE 4.0 Compiled on Sep 15 2004\n./AlignACE -i seqs.fa -numcols 6\n"
        "Input Sequences:\n#1\tseqA\n#2\tseqB\n#3\tseqC\nParameters:\nexpect  = 10\n\n"
        "Motif 1\nACGTAC\t1\t+\tseqA\nACGTAT\t1\t-\tseqB\nACCTAC\t2\t+\tseqC\n"
        "******\nMAP Score: 8.5\n",
    ),
    "clusterbuster": (
        "cb.txt",
        ">C1\n3\t0\t1\t0\n0\t4\t0\t0\n0\t0\t1\t3\n>C2\n1\t1\t1\t1\n2\t0\t2\t0\n",
    ),
    "minimal": (
        "mini.meme",
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies (from uniform background):\n"
        "A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\nMOTIF M1\n"
        "letter-probability matrix: alength= 4 w= 3 nsites= 10 E= 1.5e-3\n"
        "0.700000 0.100000 0.100000 0.100000\n0.100000 0.700000 0.100000 0.100000\n"
        "0.100000 0.100000 0.700000 0.100000\n\n",
    ),
    "xms": (
        "motif.xms",
        '<?xml version="1.0" encoding="UTF-8"?>\n<motifset>\n<motif>\n<name>X1</name>\n'
        '<column pos="0"><weight symbol="adenine">3</weight><weight symbol="cytosine">0</weight>'
        '<weight symbol="guanine">1</weight><weight symbol="thymine">0</weight></column>\n'
        '<column pos="1"><weight symbol="adenine">0</weight><weight symbol="cytosine">4</weight>'
        '<weight symbol="guanine">0</weight><weight symbol="thymine">0</weight></column>\n'
        "</motif>\n</motifset>\n",
    ),
    "pfm-four-columns": ("four_cols.txt", "3\t0\t0\t1\n0\t4\t0\t0\n0\t1\t3\t0\n"),
    "pfm-four-rows": ("four_rows.txt", "3 0 0 1\n0 4 1 0\n0 0 3 0\n1 0 0 4\n"),
}


@pytest.mark.parametrize("file_format", sorted(MOTIF_FILES))
def test_read_motifs_format_reaches_motif_codec_readers(tmp_path: Path, file_format: str) -> None:
    import numpy
    from Bio import motifs

    from dnakit.bioinformatics.motifs import read_motifs

    name, text = MOTIF_FILES[file_format]
    path = _write(tmp_path / name, text)
    with path.open(encoding="utf-8") as handle:
        direct = motifs.parse(handle, file_format)
    wrapped = read_motifs(path, format=file_format, verbose=False)
    assert len(wrapped) == len(direct)
    for model, motif in zip(wrapped, direct, strict=True):
        expected = numpy.array([motif.counts[base] for base in "ACGT"], dtype=float)
        assert numpy.allclose(model.counts, expected)


def test_write_motifs_file_format_reaches_motif_codec_writers() -> None:
    from Bio import motifs
    from Bio.motifs import clusterbuster, jaspar, pfm, transfac

    from dnakit.sequence.Bio_review_functions import write_motifs

    instances = {"M1": ["ACGTA", "ACGTT", "ACCTA", "TCGTA"], "M2": ["GGGC", "GGCC", "GAGC"]}
    built = []
    for name in sorted(instances):
        motif = motifs.create(instances[name])
        motif.name = name
        motif.matrix_id = name
        built.append(motif)
    assert write_motifs(instances, "jaspar", verbose=False)["text"] == jaspar.write(built, "jaspar")
    assert write_motifs(instances, "transfac", verbose=False)["text"] == transfac.write(built)
    cluster = write_motifs(instances, "clusterbuster", verbose=False)["text"]
    assert cluster == clusterbuster.write(built)
    # Bio.motifs.pfm.write is NOT the JASPAR "pfm" layout: its text equals the ClusterBuster layout.
    assert pfm.write(built) == cluster
    single = write_motifs({"M1": instances["M1"]}, "pfm", verbose=False)["text"]
    assert single == jaspar.write(built[:1], "pfm")
    assert single != pfm.write(built[:1])


def test_nucleotide_distance_model_reaches_every_skbio_distance() -> None:
    skbio = pytest.importorskip("skbio")
    from skbio.sequence import distance

    from dnakit.sequence.skbio_functions import nucleotide_distance

    left, right = "ACGTACGTAAGTCCGTAGCTAGCATCGA", "ACGTTCGTACGTCCATAGCTAGGATCGT"
    for model in ("f81", "f84", "jc69", "k2p", "logdet", "paralin", "tn93"):
        wrapped = nucleotide_distance(left, right, model=model, verbose=False)["distance"]
        direct = float(getattr(distance, model)(skbio.DNA(left), skbio.DNA(right)))
        assert wrapped == pytest.approx(direct), model


def test_nucleotide_transition_matrix_model_reaches_every_skbio_matrix() -> None:
    pytest.importorskip("skbio")
    from skbio.sequence import transition

    from dnakit.sequence.skbio_functions import nucleotide_transition_matrix

    freqs = [0.3, 0.2, 0.2, 0.3]
    cases = {
        "f81": ({"base_frequencies": freqs}, lambda: transition.f81(0.1, freqs, dtype="float64")),
        "hky85": (
            {"base_frequencies": freqs, "kappa": 0.5},
            lambda: transition.hky85(0.1, freqs, 0.5, dtype="float64"),
        ),
        "k2p": ({"kappa": 0.5}, lambda: transition.k2p(0.1, 0.5, dtype="float64")),
        "tn93": (
            {"base_frequencies": freqs, "kappa_purine": 0.4, "kappa_pyrimidine": 0.6},
            lambda: transition.tn93(0.1, freqs, 0.4, 0.6, dtype="float64"),
        ),
    }
    for model, (kwargs, direct) in cases.items():
        wrapped = nucleotide_transition_matrix(0.1, model=model, verbose=False, **kwargs)
        assert wrapped["matrix"] == [[float(x) for x in row] for row in direct().data.tolist()]


def test_multiple_nucleotide_alignment_parameters_reach_multi_align() -> None:
    skbio = pytest.importorskip("skbio")
    from skbio.alignment import multi_align

    from dnakit.sequence.skbio_functions import multiple_nucleotide_alignment

    sequences = {"s1": "ACGTACGTAC", "s2": "ACGACGTAC", "s3": "ACGTTCGTAC"}
    wrapped = multiple_nucleotide_alignment(
        sequences,
        match_score=2.0,
        mismatch_score=-3.0,
        gap_cost=4.0,
        free_ends=False,
        verbose=False,
    )
    result = multi_align(
        [skbio.DNA(text) for text in sequences.values()],
        sub_score=(2.0, -3.0),
        gap_cost=4.0,
        free_ends=False,
        ids=list(sequences),
    )
    aligned = [
        str(item) for item in result.path.to_aligned([skbio.DNA(t) for t in sequences.values()])
    ]
    assert list(wrapped["alignment"].values()) == aligned


def test_skbio_registry_names_are_the_same_objects() -> None:
    import skbio
    import skbio.io

    assert skbio.io.registry.read is skbio.io.read is skbio.read
    assert skbio.io.registry.write is skbio.io.write is skbio.write
    assert skbio.io.registry.sniff is skbio.io.sniff


def test_screed_fastq_parameters_reach_tofastq_and_fastq_iter(tmp_path: Path) -> None:
    import screed
    import screed.dump_fastq

    from dnakit.sequence import screed_functions
    from dnakit.sequence.backends import build_screed_database

    assert screed.ToFastq is screed.conversion.ToFastq is screed.dump_fastq.ToFastq
    assert screed.fastq.fastq_iter is screed.openscreed.fastq_iter
    text = "@q1 desc\nACGTAC\n+\nIIII#!\n@q2 \nGGCCAA\n+\nABCDEF\n"
    source = _write(tmp_path / "reads.fq", text)
    built = build_screed_database(source, database=tmp_path / "reads_fq", verbose=False)
    wrapped_out = tmp_path / "wrapped.fq"
    screed_functions.export_screed_database(
        built["database"], wrapped_out, output_format="fastq", verbose=False
    )
    direct_out = tmp_path / "direct.fq"
    screed.ToFastq(built["database"], str(direct_out))
    assert wrapped_out.read_text(encoding="utf-8") == direct_out.read_text(encoding="utf-8")
    records = screed_functions.read_screed_fastx_records(source, file_format="fastq", verbose=False)
    with source.open(encoding="utf-8") as handle:
        direct = list(screed.fastq.fastq_iter(handle))
    assert [r["name"] for r in records["records"]] == [str(r["name"]).split()[0] for r in direct]
    assert [r["sequence"] for r in records["records"]] == [str(r["sequence"]) for r in direct]


def test_cutadapt_aliases_and_parameters(tmp_path: Path) -> None:
    from cutadapt import adapters, align, cli, kmer_heuristic, modifiers, parser, qualtrim

    from dnakit.sequencing_qc.cutadapt_functions import hamming_neighbours, parse_adapter_spec

    assert adapters.edit_environment is align.edit_environment
    assert adapters.hamming_sphere is align.hamming_sphere
    assert adapters.create_positions_and_kmers is kmer_heuristic.create_positions_and_kmers
    assert adapters.kmer_probability_analysis is kmer_heuristic.kmer_probability_analysis
    assert modifiers.nextseq_trim_index is qualtrim.nextseq_trim_index
    assert modifiers.poly_a_trim_index is qualtrim.poly_a_trim_index
    assert cli.make_adapters_from_specifications is parser.make_adapters_from_specifications
    # exact=False reaches hamming_environment.
    wrapped = hamming_neighbours("ACGTT", 2, exact=False, verbose=False)
    direct = sorted((str(t), int(e), int(m)) for t, e, m in align.hamming_environment("ACGTT", 2))
    assert [(r["sequence"], r["distance"], r["matches"]) for r in wrapped["neighbours"]] == sorted(
        direct, key=lambda row: (row[1], row[0])
    )
    # parse_adapter_spec for a non-file specification is parser.make_adapter.
    search = {"max_errors": 0.1, "min_overlap": 3, "indels": True}
    adapter = parser.make_adapter("name=ACGTACGT", "back", search)
    row = parse_adapter_spec("name=ACGTACGT", adapter_type="back", verbose=False)["adapters"][0]
    assert (row["name"], row["sequence"], row["adapter_class"]) == (
        adapter.name,
        adapter.sequence,
        type(adapter).__name__,
    )


def test_cutadapt_reference_edit_environments_equal_edit_neighbours() -> None:
    from cutadapt import align

    from dnakit.sequencing_qc.cutadapt_functions import edit_neighbours

    wrapped = edit_neighbours("ACGTT", 2, verbose=False)
    pairs = {(r["sequence"], r["distance"]) for r in wrapped["neighbours"]}
    assert pairs == {(str(t), int(e)) for t, e, _ in align.edit_environment("ACGTT", 2)}
    assert pairs == {(str(t), int(e)) for t, e, _ in align.slow_edit_environment("ACGTT", 2)}
    assert pairs == {(str(t), int(e)) for t, e, _ in align.py_edit_environment("ACGTT", 2)}
    naive = list(align.naive_edit_environment("ACGTT", 2))
    assert set(naive) == {sequence for sequence, _ in pairs}
    assert len(naive) > len(set(naive))  # the naive form repeats strings and yields no distances


def test_umi_tools_read_name_and_knee_parameters() -> None:
    from umi_tools import extract_methods, network, sam_methods, whitelist_methods

    from dnakit.sequencing_qc.umi_tools_functions import (
        call_cell_whitelist,
        mask_low_quality_umi,
        parse_umi_read_name,
        umi_connected_component,
    )

    read = SimpleNamespace(qname="read7_AAAC_GGTT")
    for cell, sep in ((False, "_"), (True, "_")):
        umi, barcode = sam_methods.get_barcode_read_id(read, cell_barcode=cell, sep=sep)
        wrapped = parse_umi_read_name(read.qname, separator=sep, cell_barcode=cell, verbose=False)
        assert wrapped["umi"] == umi.decode()
        assert wrapped["cell_barcode"] == (barcode.decode() if barcode else None)
    assert sam_methods.get_umi_read_string("read7_AAAC_GGTT") == (None, b"GGTT")
    assert parse_umi_read_name("read7_AAAC_GGTT", verbose=False)["umi"] == "GGTT"
    assert sam_methods.get_cell_umi_read_string("read7_AAAC_GGTT") == (b"AAAC", b"GGTT")

    rng = random.Random(1)
    counts = {}
    for index in range(40):
        barcode = "".join(rng.choice("ACGT") for _ in range(8))
        counts[barcode] = rng.randint(300, 500) if index < 6 else rng.randint(1, 6)
    import collections

    knee_distance = whitelist_methods.getKneeEstimateDistance(collections.Counter(counts))
    knee_density = whitelist_methods.getKneeEstimateDensity(collections.Counter(counts))
    assert sorted(
        call_cell_whitelist(counts, knee_method="distance", verbose=False)["whitelist"]
    ) == sorted(knee_distance)
    assert sorted(
        call_cell_whitelist(counts, knee_method="density", verbose=False)["whitelist"]
    ) == sorted(knee_density)

    # get_below_threshold and umi_below_threshold are contained in mask_low_quality_umi.
    for umi, quals in (("ACGTAC", "II#IJ5"), ("ACGTAC", "IIIIII"), ("AAAA", "!!!!")):
        flags = extract_methods.get_below_threshold(quals, "phred33", 10)
        wrapped = mask_low_quality_umi(umi, quals, threshold=10, verbose=False)
        assert wrapped["masked_positions"] == [i for i, flag in enumerate(flags) if flag]
        assert (wrapped["masked_count"] > 0) is extract_methods.umi_below_threshold(
            quals, "phred33", 10
        )

    # The recursive searches return the component that umi_connected_component returns.
    graph = {"a": ["b"], "b": ["a", "c"], "c": ["b"], "d": []}
    component = umi_connected_component(graph, "a", verbose=False)
    assert component["component"] == sorted(network.breadth_first_search_recursive("a", graph))
    assert set(network.recursive_search("a", graph)) | {"a"} == set(component["component"])


def test_khmer_parameters_reach_the_variant_functions(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    from khmer import khmer_args

    from dnakit.microbiome.khmer_functions import (
        encode_kmer_2bit,
        estimate_kmer_table_size,
        hash_kmer_murmur3,
        inspect_khmer_graph,
    )

    for kmer in ("ACGTA", "TACGT", "TTTTT", "GGGAC"):
        assert encode_kmer_2bit(kmer, canonical=False, verbose=False)["code"] == int(
            khmer.forward_hash_no_rc(kmer, len(kmer))
        )
        assert hash_kmer_murmur3(kmer, canonical=False, verbose=False)["hash"] == int(
            khmer.hash_no_rc_murmur3(kmer)
        )
    by_memory = estimate_kmer_table_size(1_000_000, mem_cap=4_000_000, verbose=False)
    direct_memory = khmer_args.estimate_optimal_with_K_and_M(1_000_000, 4_000_000)
    assert (by_memory["num_htables"], by_memory["htable_size"]) == (
        direct_memory.num_htables,
        direct_memory.htable_size,
    )
    by_rate = estimate_kmer_table_size(1_000_000, fp_rate=0.01, verbose=False)
    direct_rate = khmer_args.estimate_optimal_with_K_and_f(1_000_000, 0.01)
    assert (by_rate["num_htables"], by_rate["htable_size"]) == (
        direct_rate.num_htables,
        direct_rate.htable_size,
    )
    graph = khmer.Nodegraph(5, 1000, 4)
    graph.consume(KHMER_SEQ)
    path = tmp_path / "graph.ng"
    graph.save(str(path))
    info = inspect_khmer_graph(path, kind="nodegraph", verbose=False)
    assert tuple(khmer.extract_nodegraph_info(str(path))) == (
        info["ksize"],
        info["table_size"],
        info["n_tables"],
        info["format_version"],
        info["table_type_flag"],
        info["n_occupied"],
    )


# =============================== claims changed from family-variant ============================


def test_biopython_deprecated_formatters_warn() -> None:
    from Bio.Seq import Seq
    from Bio.SeqIO import QualityIO, TabIO
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq("ACGT"), id="r", description="", letter_annotations={"phred_quality": [30] * 4}
    )
    for function in (QualityIO.as_fastq_illumina, QualityIO.as_fastq_solexa, TabIO.as_tab):
        with pytest.warns(DeprecationWarning, match="deprecated"):
            function(record)
    assert format(record, "tab") == "r\tACGT\n"


def test_gafbyproteiniterator_drops_the_last_protein_group() -> None:
    from Bio.UniProt import GOA

    from dnakit.sequence.Bio_functions import read_gaf_annotations

    rows = []
    for protein, term in (
        ("P1", "GO:0000001"),
        ("P1", "GO:0000002"),
        ("P2", "GO:0000003"),
        ("P3", "GO:0000004"),
    ):
        fields = ["UniProtKB", protein, "S", "", term, "PMID:1", "IEA", "", "F", "N", "", "protein"]
        fields += ["taxon:9606", "20240101", "UniProt"]
        rows.append("\t".join(fields))
    text = "!gaf-version: 2.2\n" + "\n".join(rows) + "\n"
    groups = [
        [r["DB_Object_ID"] for r in group] for group in GOA.gafbyproteiniterator(io.StringIO(text))
    ]
    assert groups == [["P1", "P1"], ["P2"]]  # P3 is lost: upstream never yields the last group
    flat = [r["DB_Object_ID"] for r in GOA.gafiterator(io.StringIO(text))]
    assert flat == ["P1", "P1", "P2", "P3"]
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        path = _write(Path(directory) / "x.gaf", text)
        wrapped = read_gaf_annotations(path, verbose=False)
    assert [r["db_object_id"] for r in wrapped["records"]] == flat


def test_ace_parse_contigs_are_the_contigs_of_read_ace_assembly(tmp_path: Path) -> None:
    from Bio.Sequencing import Ace

    from dnakit.sequence.Bio_functions import read_ace_assembly

    text = (
        "AS 1 2\n\nCO Contig1 12 2 4 U\nACGTACGTACGT\n\nBQ\n"
        + " ".join(["30"] * 12)
        + "\n\nAF read1 U 1\nAF read2 U 2\n\n"
        "RD read1 12 0 0\nACGTACGTACGT\n\nQA 1 12 1 12\n"
        "DS CHROMAT_FILE read1 PHD_FILE read1.phd.1 TIME Thu Jan 1 00:00:00 2024\n\n"
        "RD read2 12 0 0\nACGTACGTACGT\n\nQA 1 12 1 12\n"
        "DS CHROMAT_FILE read2 PHD_FILE read2.phd.1 TIME Thu Jan 1 00:00:00 2024\n\n"
    )
    path = _write(tmp_path / "contig.ace", text)
    with path.open(encoding="utf-8") as handle:
        contigs = list(Ace.parse(handle))
    wrapped = read_ace_assembly(path, verbose=False)
    assert [c.name for c in contigs] == [c["name"] for c in wrapped["contigs"]]
    assert [str(c.sequence) for c in contigs] == [c["sequence"] for c in wrapped["contigs"]]
    assert "list(parse(handle))" in inspect.getsource(Ace.read)


def test_goa_writers_are_called_by_write_gaf_file(tmp_path: Path) -> None:
    from Bio.UniProt import GOA

    from dnakit.sequence.Bio_review_functions import write_gaf_file

    record = {
        "DB": "UniProtKB",
        "DB_Object_ID": "P12345",
        "DB_Object_Symbol": "GENEA",
        "Qualifier": ["enables"],
        "GO_ID": "GO:0003700",
        "DB:Reference": ["PMID:1"],
        "Evidence": "IDA",
        "With": [],
        "Aspect": "F",
        "DB_Object_Name": "TF A",
        "Synonym": ["TFA"],
        "DB_Object_Type": "protein",
        "Taxon_ID": ["taxon:9606"],
        "Date": "20240101",
        "Assigned_By": "UniProt",
        "Annotation_Extension": "",
        "Gene_Product_Form_ID": "",
    }
    path = tmp_path / "out.gaf"
    write_gaf_file([record, record], path, verbose=False)
    direct = io.StringIO()
    GOA.writebyproteinrec([dict(record, With="") for _ in range(2)], direct, fields=list(record))
    assert path.read_text(encoding="utf-8") == direct.getvalue()
    single = io.StringIO()
    GOA.writerec(dict(record, With=""), single, fields=list(record))
    assert single.getvalue() == direct.getvalue().splitlines(keepends=True)[0]


def test_seqio_read_and_to_dict_return_the_records_of_read_sequence_records(tmp_path: Path) -> None:
    from Bio import SeqIO

    from dnakit.sequence.backends import read_sequence_records

    path = _write(tmp_path / "one.gb", GENBANK_TEXT.split("//\n")[0] + "//\n")
    wrapped = read_sequence_records(path, format="genbank", verbose=False)["records"]
    record = SeqIO.read(str(path), "genbank")
    assert [wrapped[0]["id"], wrapped[0]["sequence"]] == [record.id, str(record.seq).upper()]
    both = _write(tmp_path / "two.gb", GENBANK_TEXT)
    table = SeqIO.to_dict(SeqIO.parse(str(both), "genbank"))
    all_records = read_sequence_records(both, format="genbank", verbose=False)["records"]
    assert sorted(table) == sorted(r["id"] for r in all_records)


def test_to_igraph_numbering_differs_from_tree_to_graph_edges() -> None:
    pytest.importorskip("igraph")
    from Bio import Phylo

    from dnakit.sequence.Bio_review_functions import tree_to_graph_edges

    newick = "((A:0.1,B:0.2)X:0.05,(C:0.3,D:0.4):0.06);"
    graph = Phylo.to_igraph(
        Phylo.read(io.StringIO(newick), "newick"),
        vertex_attributes=("name",),
        edge_attributes=("branch_length",),
    )
    wrapped = tree_to_graph_edges(newick, verbose=False)
    assert graph.ecount() == wrapped["edge_count"] and graph.vcount() == wrapped["node_count"]
    labels = graph.vs["name"]
    igraph_edges = sorted(
        (str(labels[a]), str(labels[b]), length)
        for (a, b), length in zip(graph.get_edgelist(), graph.es["branch_length"], strict=True)
    )
    names = {node["id"]: str(node["name"]) for node in wrapped["nodes"]}
    wrapped_edges = sorted(
        (names[e["source"]], names[e["target"]], e["branch_length"]) for e in wrapped["edges"]
    )
    # The named (tip and X) edges agree; the vertex ids of the two layouts differ.
    named = lambda rows: [row for row in rows if row[0] != "None" and row[1] != "None"]  # noqa: E731
    assert named(igraph_edges) == named(wrapped_edges)
    assert graph.get_edgelist() != [(e["source"], e["target"]) for e in wrapped["edges"]]
