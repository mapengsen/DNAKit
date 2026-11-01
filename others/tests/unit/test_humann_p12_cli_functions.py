"""Tests of the HUMAnN facades: validation, argument vectors (stand-in programs) and real HUMAnN utility runs."""

from __future__ import annotations

import json
import os
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.microbiome import humann_p12_cli_functions as hm

REAL = (
    shutil.which("humann_join_tables") is not None
    or (Path(sys.executable).parent / "humann_join_tables").is_file()
)
needs_humann = pytest.mark.skipif(not REAL, reason="HUMAnN is not installed")
EXAMPLES = os.environ.get("HUMANN_EXAMPLES_DIR", "")
needs_demo = pytest.mark.skipif(
    not REAL or not (Path(EXAMPLES) / "demo.m8").is_file(),
    reason="HUMAnN or its examples folder (HUMANN_EXAMPLES_DIR) is not available",
)

TABLE_A = "# Gene Family\tA_Abundance-RPKs\nUNMAPPED\t1.0\nUniRef90_X\t4.0\nUniRef90_X|g__Bacteroides.s__Bacteroides_vulgatus\t3.0\nUniRef90_X|unclassified\t1.0\nUniRef90_Y\t2.0\nUniRef90_Y|unclassified\t2.0\n"
TABLE_B = "# Gene Family\tB_Abundance-RPKs\nUNMAPPED\t2.0\nUniRef90_X\t6.0\nUniRef90_X|g__Bacteroides.s__Bacteroides_vulgatus\t1.0\nUniRef90_X|unclassified\t5.0\nUniRef90_Y\t0.0\n"


@pytest.fixture()
def tables(tmp_path: Path) -> dict[str, Path]:
    folder = tmp_path / "tables"
    folder.mkdir()
    (folder / "A_genefamilies.tsv").write_text(TABLE_A)
    (folder / "B_genefamilies.tsv").write_text(TABLE_B)
    return {"dir": folder, "root": tmp_path}


def _stand_in(folder: Path, name: str, create: tuple[str, ...] = (), stdout: str = "") -> Path:
    record = folder / "record.jsonl"
    script = folder / "bin" / name
    script.parent.mkdir(exist_ok=True)
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, pathlib, sys\n"
        f"open({str(record)!r}, 'a').write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:] == ['--version']:\n"
        "    print('humann v9.9'); sys.exit(0)\n"
        f"for path in {list(create)!r}:\n"
        "    pathlib.Path(path).write_text('# Gene Family\\tS_Abundance\\nUniRef90_X\\t1.0\\n')\n"
        f"sys.stderr.write({stdout!r})\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _calls(folder: Path) -> list[list[str]]:
    return [json.loads(line) for line in (folder / "record.jsonl").read_text().splitlines()]


def test_public_names_have_verbose() -> None:
    import inspect

    assert len(hm.__all__) == 13
    for name in hm.__all__:
        assert inspect.signature(getattr(hm, name)).parameters["verbose"].default is True


def test_missing_humann_dir_and_programs(tables: dict[str, Path]) -> None:
    out = tables["root"] / "joined.tsv"
    with pytest.raises(ConfigurationError):
        hm.humann_join_tables(
            str(tables["dir"]), str(out), humann_dir="/nonexistent", verbose=False
        )
    empty = tables["root"] / "emptybin"
    empty.mkdir()
    with pytest.raises(BackendUnavailableError):
        hm.humann_join_tables(str(tables["dir"]), str(out), humann_dir=str(empty), verbose=False)


def test_validation_of_choices_and_pairs(tables: dict[str, Path]) -> None:
    source = tables["dir"] / "A_genefamilies.tsv"
    out = tables["root"] / "o.tsv"
    with pytest.raises(ConfigurationError):
        hm.humann_quantify_alignments(
            str(source), str(tables["root"] / "q"), input_format="fastq", verbose=False
        )
    with pytest.raises(ConfigurationError):
        hm.humann_quantify_alignments(
            str(source), str(tables["root"] / "q2"), output_basename="../x", verbose=False
        )
    with pytest.raises(ConfigurationError, match="exactly one"):
        hm.humann_rename_table(str(source), str(out), verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        hm.humann_rename_table(
            str(source), str(out), names="ec", custom_mapping=str(source), verbose=False
        )
    with pytest.raises(ConfigurationError):
        hm.humann_rename_table(str(source), str(out), names="not-a-mapping", verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        hm.humann_regroup_table(str(source), str(out), verbose=False)
    with pytest.raises(ConfigurationError):
        hm.humann_regroup_table(
            str(source), str(out), groups="uniref90_rxn", function="median", verbose=False
        )
    with pytest.raises(ConfigurationError):
        hm.humann_renorm_table(str(source), str(out), units="tpm", verbose=False)
    with pytest.raises(ConfigurationError):
        hm.humann_reduce_table(str(source), str(out), sort_by="size", verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        hm.humann_infer_taxonomy(str(source), str(out), verbose=False)
    with pytest.raises(ConfigurationError):
        hm.humann_expand_cluster(str(source), "bad id!", str(out), verbose=False)
    with pytest.raises(ConfigurationError, match="extension"):
        hm.humann_barplot(
            str(source), "UniRef90_X", str(tables["root"] / "plot.txt"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        hm.humann_barplot(
            str(source),
            "UniRef90_X",
            str(tables["root"] / "p.png"),
            sort=["alphabetical"],
            verbose=False,
        )


def test_quantify_argument_vector(tables: dict[str, Path]) -> None:
    source = tables["root"] / "sample.m8"
    source.write_text("x\n")
    out = tables["root"] / "q"
    base = out / "sample"
    script = _stand_in(
        tables["root"],
        "humann",
        create=(),
    )
    # the stand-in writes the three outputs after the output folder exists: pre-create through a second helper
    script.write_text(
        script.read_text().replace(
            "sys.stderr.write(",
            "[pathlib.Path(p).write_text('# Gene Family\\tS\\nUniRef90_X\\t1.0\\n') for p in "
            f"[{str(base) + '_genefamilies.tsv'!r}, {str(base) + '_pathabundance.tsv'!r}, {str(base) + '_pathcoverage.tsv'!r}]]\n"
            "sys.stderr.write(",
        )
    )
    result = hm.humann_quantify_alignments(
        str(source),
        str(out),
        input_format="blastm8",
        minpath=False,
        gap_fill=False,
        pathways="unipathway",
        output_max_decimals=4,
        remove_stratified_output=True,
        threads=3,
        humann_dir=str(script.parent),
        verbose=False,
    )
    call = _calls(tables["root"])[0]
    assert call[:6] == ["--input", str(source), "--output", str(out), "--input-format", "blastm8"]
    for pair in (["--minpath", "off"], ["--gap-fill", "off"], ["--pathways", "unipathway"], ["--threads", "3"],
                 ["--output-max-decimals", "4"], ["--output-basename", "sample"]):  # fmt: skip
        assert any(call[i : i + 2] == pair for i in range(len(call) - 1)), pair
    assert "--remove-stratified-output" in call and "--xipe" not in call
    assert result["gene_families"]["row_count"] == 1 and result["version"] == "9.9"
    json.dumps(result)


def test_genetable_input_has_no_gene_family_output(tables: dict[str, Path]) -> None:
    source = tables["root"] / "genes.tsv"
    source.write_text(TABLE_A)
    out = tables["root"] / "q"
    base = out / "genes"
    script = _stand_in(tables["root"], "humann")
    script.write_text(
        script.read_text().replace(
            "sys.stderr.write(",
            "[pathlib.Path(p).write_text('# Pathway\\tS\\nPWY\\t1.0\\n') for p in "
            f"[{str(base) + '_pathabundance.tsv'!r}, {str(base) + '_pathcoverage.tsv'!r}]]\n"
            "sys.stderr.write(",
        )
    )
    result = hm.humann_quantify_alignments(
        str(source),
        str(out),
        input_format="genetable",
        humann_dir=str(script.parent),
        verbose=False,
    )
    assert "gene_families" not in result and result["pathway_abundance"]["rows"] == [["PWY", 1.0]]


def test_failure_is_reported(tables: dict[str, Path]) -> None:
    script = tables["root"] / "humann_join_tables"
    script.write_text(
        f"#!{sys.executable}\nimport sys\nsys.stderr.write('no tables')\nsys.exit(1)\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    with pytest.raises(BackendExecutionError) as caught:
        hm.humann_join_tables(
            str(tables["dir"]),
            str(tables["root"] / "j.tsv"),
            humann_dir=str(tables["root"]),
            verbose=False,
        )
    assert "no tables" in caught.value.context["log_tail"]


@needs_humann
def test_join_split_and_renorm(tables: dict[str, Path]) -> None:
    root = tables["root"]
    joined = hm.humann_join_tables(str(tables["dir"]), str(root / "joined.tsv"), verbose=False)
    assert joined["table"]["columns"] == ["Gene Family", "A_Abundance-RPKs", "B_Abundance-RPKs"]
    assert joined["table"]["row_count"] == 6
    only_a = hm.humann_join_tables(
        str(tables["dir"]), str(root / "only_a.tsv"), file_name="A_", verbose=False
    )
    assert len(only_a["table"]["columns"]) == 2 and only_a["table"]["columns"][1].startswith("A_")
    split = hm.humann_split_table(str(root / "joined.tsv"), str(root / "split"), verbose=False)
    assert split["table_count"] == 2
    renormed = hm.humann_renorm_table(
        str(root / "joined.tsv"), str(root / "relab.tsv"), units="relab", verbose=False
    )
    first = dict((row[0], row[1:]) for row in renormed["table"]["rows"])
    assert first["UniRef90_X"][0] == pytest.approx(
        4.0 / 7.0, rel=1e-3
    )  # community total of sample A = 1 + 4 + 2 (UNMAPPED included)
    cpm = hm.humann_renorm_table(
        str(root / "joined.tsv"), str(root / "cpm.tsv"), include_special=False, verbose=False
    )
    assert "UNMAPPED" not in [row[0] for row in cpm["table"]["rows"]]
    stratified = hm.humann_split_stratified_table(
        str(root / "joined.tsv"), str(root / "strata"), verbose=False
    )
    assert (
        stratified["unstratified"]["row_count"] == 3 and stratified["stratified"]["row_count"] == 3
    )


@needs_humann
def test_rename_regroup_reduce(tables: dict[str, Path]) -> None:
    root = tables["root"]
    hm.humann_join_tables(str(tables["dir"]), str(root / "joined.tsv"), verbose=False)
    (root / "names.tsv").write_text("UniRef90_X\tgene x\n")
    renamed = hm.humann_rename_table(
        str(root / "joined.tsv"),
        str(root / "renamed.tsv"),
        custom_mapping=str(root / "names.tsv"),
        verbose=False,
    )
    assert renamed["table"]["rows"][1][0] == "UniRef90_X: gene x"
    (root / "groups.tsv").write_text("GROUP1\tUniRef90_X\tUniRef90_Y\n")
    regrouped = hm.humann_regroup_table(
        str(root / "joined.tsv"),
        str(root / "regrouped.tsv"),
        custom_groups=str(root / "groups.tsv"),
        verbose=False,
    )
    rows = {row[0]: row[1:] for row in regrouped["table"]["rows"]}
    assert rows["GROUP1"] == [6.0, 6.0] and rows["UNMAPPED"] == [1.0, 2.0]
    no_ungrouped = hm.humann_regroup_table(
        str(root / "joined.tsv"),
        str(root / "regrouped2.tsv"),
        custom_groups=str(root / "groups.tsv"),
        include_ungrouped=False,
        function="mean",
        verbose=False,
    )
    assert {row[0]: row[1:] for row in no_ungrouped["table"]["rows"]}["GROUP1"] == [3.0, 3.0]
    assert "UNGROUPED" not in [row[0] for row in no_ungrouped["table"]["rows"]]
    reduced = hm.humann_reduce_table(
        str(root / "joined.tsv"), str(root / "reduced.tsv"), function="sum", verbose=False
    )
    rows = {row[0]: row[1] for row in reduced["table"]["rows"]}
    assert rows["UniRef90_X"] == 10.0 and rows["UniRef90_Y"] == 2.0
    maximum = hm.humann_reduce_table(str(root / "joined.tsv"), str(root / "max.tsv"), verbose=False)
    assert {row[0]: row[1] for row in maximum["table"]["rows"]}["UniRef90_X"] == 6.0


@needs_humann
def test_infer_taxonomy_with_a_development_database(tables: dict[str, Path]) -> None:
    root = tables["root"]
    hm.humann_join_tables(str(tables["dir"]), str(root / "joined.tsv"), verbose=False)
    database = root / "tol.tsv"
    database.write_text(
        "TAXID\tNAME\tRANK\tPARENT\nRoot\tRoot\tRoot\tRoot\n2\tBacteria\tKingdom\tRoot\n200\tProteobacteria\tPhylum\t2\n"
        "FAMILY\tSOURCE_TAX\tUNIREF_LCA\tHUMANN_LCA\nUniRef90_X\t2\t200\t200\n"
    )
    result = hm.humann_infer_taxonomy(
        str(root / "joined.tsv"),
        str(root / "inferred.tsv"),
        database_file=str(database),
        level="Phylum",
        verbose=False,
    )
    assert "p__Proteobacteria" in "".join(row[0] for row in result["table"]["rows"])
    assert result["summary"].startswith("Summary: Of")


@needs_humann
def test_build_custom_database_and_barplot(tables: dict[str, Path]) -> None:
    root = tables["root"]
    fasta = root / "genes.fa"
    fasta.write_text(">Bacteroides_vulgatus|g1\nACGTACGT\n>Escherichia_coli|g2\nACGTTTGG\n")
    mapping = root / "ids.tsv"
    mapping.write_text("g1\tBacteroides\ng2\tEscherichia\n")
    built = hm.humann_build_custom_database(
        str(fasta), str(root / "db"), id_mapping=str(mapping), verbose=False
    )
    assert built["database_files"] == [str(root / "db" / "genes.filtered.fa")]
    assert (root / "db" / "genes.filtered.fa").read_text().count(">") == 2
    hm.humann_join_tables(str(tables["dir"]), str(root / "joined.tsv"), verbose=False)
    plot = hm.humann_barplot(
        str(root / "joined.tsv"),
        "UniRef90_X",
        str(root / "bar.png"),
        top_taxa=2,
        sort=["sum"],
        sample_order_output=str(root / "order.txt"),
        verbose=False,
    )
    assert plot["image_bytes"] > 1000
    assert sorted((root / "order.txt").read_text().split()) == [
        "A_Abundance-RPKs",
        "B_Abundance-RPKs",
    ]


@needs_demo
def test_quantify_demo_alignments_and_unpack(tmp_path: Path) -> None:
    demo = Path(EXAMPLES) / "demo.m8"
    result = hm.humann_quantify_alignments(
        str(demo), str(tmp_path / "q"), minpath=False, verbose=False
    )
    assert result["status"] == "completed" and result["version"]
    assert (
        result["gene_families"]["row_count"] > 50 and result["pathway_abundance"]["row_count"] > 5
    )
    again = hm.humann_quantify_alignments(
        result["gene_families"]["file"],
        str(tmp_path / "q2"),
        input_format="genetable",
        minpath=False,
        verbose=False,
    )
    assert again["pathway_abundance"]["row_count"] == result["pathway_abundance"]["row_count"]
    unpacked = hm.humann_unpack_pathways(
        result["gene_families"]["file"],
        result["pathway_abundance"]["file"],
        str(tmp_path / "up.tsv"),
        verbose=False,
    )
    assert unpacked["table"]["row_count"] > result["pathway_abundance"]["row_count"]
