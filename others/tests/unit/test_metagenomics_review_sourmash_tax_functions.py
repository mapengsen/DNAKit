"""Tests for dnakit.metagenomics_review_sourmash_tax_functions (sourmash tax metagenome output formats)."""

from __future__ import annotations

import inspect
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit import metagenomics_review_sourmash_tax_functions as tax
from dnakit.exceptions import ConfigurationError

_SOURMASH = Path(sys.executable).parent / "sourmash"
pytestmark = pytest.mark.skipif(
    not _SOURMASH.is_file() and shutil.which("sourmash") is None, reason="sourmash is not installed"
)

_LINEAGES = (
    "ident,superkingdom,phylum,class,order,family,genus,species,strain\n"
    "g1,d__Bacteria,p__Firmicutes,c__Bacilli,o__Bacillales,f__Bacillaceae,g__Bacillus,s__Bacillus one,\n"
    "g2,d__Bacteria,p__Proteobacteria,c__Gamma,o__Enterob,f__Enterobacteriaceae,g__Escherichia,s__Escherichia two,\n"
    "g3,d__Bacteria,p__Proteobacteria,c__Gamma,o__Enterob,f__Enterobacteriaceae,g__Klebsiella,s__Klebsiella three,\n"
)
_TAXPATH = (
    "ident,superkingdom,phylum,class,order,family,genus,species,strain,taxpath\n"
    "g1,Bacteria,Firmicutes,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus one,,2|1239|91061|1385|186817|1386|1386|\n"
    "g2,Bacteria,Proteobacteria,Gamma,Enterob,Enterobacteriaceae,Escherichia,Escherichia two,,2|1224|1236|91347|543|561|562|\n"
    "g3,Bacteria,Proteobacteria,Gamma,Enterob,Enterobacteriaceae,Klebsiella,Klebsiella three,,2|1224|1236|91347|543|570|571|\n"
)


def _sm(*args: str) -> None:
    subprocess.run([str(_SOURMASH), *args], capture_output=True, text=True, timeout=300, check=True)


@pytest.fixture(scope="module")
def gather(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("sourmash_review_tax")
    genomes = {}
    for number in (1, 2, 3):
        rng = random.Random(number)
        seq = "".join(rng.choice("ACGT") for _ in range(30000))
        path = root / f"g{number}.fa"
        path.write_text(f">g{number}\n{seq}\n")
        genomes[number] = seq
    header = "name,genome_filename,protein_filename\n"
    csv = root / "genomes.csv"
    csv.write_text(header + "".join(f"g{i},{root / f'g{i}.fa'},\n" for i in (1, 2, 3)))
    _sm("sketch", "fromfile", str(csv), "--output-signatures", str(root / "db.zip"), "-p", "dna,k=31,scaled=100")
    _sm("index", str(root / "db.sbt.zip"), str(root / "db.zip"), "-k", "31", "--scaled", "100")
    # query: half of genome 1 plus random sequence, so part of the query is unassigned
    rng = random.Random(99)
    query = genomes[1][:15000] + "".join(rng.choice("ACGT") for _ in range(15000))
    (root / "query.fa").write_text(f">q\n{query}\n")
    _sm("sketch", "dna", str(root / "query.fa"), "-p", "k=31,scaled=100", "-o", str(root / "query.sig"))
    _sm(
        "gather", str(root / "query.sig"), str(root / "db.sbt.zip"), "-k", "31", "--scaled", "100",
        "-o", str(root / "gather.csv"), "--threshold-bp", "0",
    )
    (root / "lineages.csv").write_text(_LINEAGES)
    (root / "taxpath.csv").write_text(_TAXPATH)
    blank = "ident,superkingdom,phylum,class,order,family,genus,species,strain,taxpath\n"
    for line in _LINEAGES.splitlines()[1:]:
        blank += line + ",\n"
    (root / "blank_taxpath.csv").write_text(blank)
    lin = "ident,superkingdom,phylum,class,order,family,genus,species,strain,lin\n"
    for number, line in zip((1, 2, 3), _LINEAGES.splitlines()[1:]):
        lin += line + f",{number}\n"
    (root / "lin.csv").write_text(lin)
    (root / "lingroups.csv").write_text("name,lin\nBacillus,1\n")
    return {"root": root, "gather": root / "gather.csv", "lineages": root / "lineages.csv",
            "taxpath": root / "taxpath.csv", "blank": root / "blank_taxpath.csv",
            "lin": root / "lin.csv", "lingroups": root / "lingroups.csv"}


def test_krona_profile_requires_a_rank(gather: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        tax.krona_taxonomy_profile([str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "k"), rank="")


def test_krona_profile_writes_fractions_by_lineage(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.krona_taxonomy_profile(
        [str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "sample"), rank="genus"
    )
    assert result["status"] == "completed", result["log_tail"]
    lines = (tmp_path / "sample.krona.tsv").read_text().splitlines()
    assert lines[0].split("\t")[0] == "fraction"
    assert any("g__Bacillus" in line for line in lines[1:])


def test_csv_summary_has_per_rank_rows(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.csv_taxonomy_summary([str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "csv"))
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "csv.summarized.csv").read_text()
    assert text.splitlines()[0].startswith("query_name,rank,fraction,lineage")
    assert ",genus," in text


def test_human_summary_and_rank_option(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.human_taxonomy_summary(
        [str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "human"), rank="genus"
    )
    assert result["status"] == "completed", result["log_tail"]
    assert "g__Bacillus" in (tmp_path / "human.human.txt").read_text()


def test_lineage_summary_requires_valid_rank(gather: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        tax.lineage_taxonomy_summary(
            [str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "l"), rank="kingdomish"
        )


def test_lineage_summary_at_rank(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.lineage_taxonomy_summary(
        [str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "lin"), rank="genus"
    )
    assert result["status"] == "completed", result["log_tail"]
    assert (tmp_path / "lin.lineage_summary.tsv").read_text().startswith("lineage")


def test_bioboxes_needs_taxpath_column(gather: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        tax.bioboxes_taxonomy_profile([str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "b"))


def test_bioboxes_profile_with_ncbi_taxids(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.bioboxes_taxonomy_profile([str(gather["gather"])], [str(gather["taxpath"])], str(tmp_path / "bb"))
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "bb.bioboxes.profile").read_text()
    assert "@Ranks:superkingdom|phylum|class|order|family|genus|species|strain" in text
    assert "1386\tgenus\t" in text


def test_bioboxes_blank_taxids_fail_upstream(gather: dict[str, Path], tmp_path: Path) -> None:
    # sourmash 4.9.4 raises TypeError when the reported lineage has no taxid; the function reports it as failed
    result = tax.bioboxes_taxonomy_profile([str(gather["gather"])], [str(gather["blank"])], str(tmp_path / "bb2"))
    assert result["status"] == "failed"
    assert "TypeError" in result["log_tail"]


def test_kreport_has_unclassified_row(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.kreport_taxonomy_report([str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "kr"))
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "kr.kreport.txt").read_text()
    assert "unclassified" in text
    assert "d__Bacteria" in text


def test_lingroup_report_needs_lin_column(gather: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        tax.lingroup_taxonomy_report(
            [str(gather["gather"])], [str(gather["lineages"])], str(gather["lingroups"]), str(tmp_path / "lg")
        )


def test_lingroup_report_with_lin_taxonomy(gather: dict[str, Path], tmp_path: Path) -> None:
    result = tax.lingroup_taxonomy_report(
        [str(gather["gather"])], [str(gather["lin"])], str(gather["lingroups"]), str(tmp_path / "lg")
    )
    assert result["status"] == "completed", result["log_tail"]
    assert (tmp_path / "lg.lingroup.tsv").is_file()


def test_existing_output_is_refused(gather: dict[str, Path], tmp_path: Path) -> None:
    (tmp_path / "taken.summarized.csv").write_text("x")
    with pytest.raises(ConfigurationError):
        tax.csv_taxonomy_summary([str(gather["gather"])], [str(gather["lineages"])], str(tmp_path / "taken"))


def test_public_functions_take_keyword_only_verbose() -> None:
    public = [
        tax.krona_taxonomy_profile, tax.csv_taxonomy_summary, tax.human_taxonomy_summary,
        tax.lineage_taxonomy_summary, tax.bioboxes_taxonomy_profile,
        tax.kreport_taxonomy_report, tax.lingroup_taxonomy_report,
    ]
    for function in public:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__
