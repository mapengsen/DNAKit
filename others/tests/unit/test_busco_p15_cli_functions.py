"""BUSCO facade tests: BUSCO and its lineage datasets are not available here, so a stand-in program that follows BUSCO's command line
and output file layout checks argument assembly and result parsing."""

from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

import pytest

from dnakit.assembly import busco_p15_cli_functions as bu
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

STUB = r'''#!{python}
import json, os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("BUSCO 5.8.2"); sys.exit(0)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "args.json"), "w").write(json.dumps(args))
def opt(n): return args[args.index(n) + 1]
run = os.path.join(opt("--out_path"), opt("-o"))
os.makedirs(os.path.join(run, "run_bacteria_odb10"))
open(os.path.join(run, "short_summary.specific.bacteria_odb10.%s.txt" % opt("-o")), "w").write("""# BUSCO version is: 5.8.2
\tC:98.4%[S:97.6%,D:0.8%],F:1.0%,M:0.6%,n:124\t
\t122\tComplete BUSCOs (C)
\t121\tComplete and single-copy BUSCOs (S)
\t1\tComplete and duplicated BUSCOs (D)
\t1\tFragmented BUSCOs (F)
\t1\tMissing BUSCOs (M)
\t124\tTotal BUSCO groups searched

Assembly Statistics:
\t1\tNumber of scaffolds
\t1\tNumber of contigs
\t4641652\tTotal length
\t0.000%\tPercent gaps
\t4 MB\tScaffold N50
""")
open(os.path.join(run, "run_bacteria_odb10", "full_table.tsv"), "w").write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n# comment\nb1\tComplete\tchr\t1\t99\t+\t50.0\t99\n")
if "--scaffold_composition" in args:
    open(os.path.join(run, "scaffold_composition.txt"), "w").write("Scaffold\tA\tC\tG\tT\tN\tothers\nchr\t10\t10\t10\t10\t0\t0\n")
'''


def test_missing_program(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    fasta.write_text(">c\nACGT\n")
    with pytest.raises(BackendUnavailableError):
        bu.busco_assess_completeness(str(fasta), "bacteria_odb10", str(tmp_path / "o"), busco_path=str(tmp_path / "missing"), verbose=False)
    with pytest.raises(ConfigurationError):
        bu.busco_assess_completeness(str(fasta), "bacteria_odb10", str(tmp_path / "o2"), gene_predictor="sapling", verbose=False)
    with pytest.raises(ConfigurationError):
        bu.busco_assess_completeness(str(fasta), "bacteria_odb10", str(tmp_path / "o3"), mode="transcriptome", gene_predictor="augustus", verbose=False)
    with pytest.raises(ConfigurationError):
        bu.busco_assess_completeness(str(fasta), "bacteria_odb10", str(tmp_path / "o4"), augustus_species="human", verbose=False)


def test_arguments_and_parsing(tmp_path: Path) -> None:
    stub = tmp_path / "bin" / "busco"
    stub.parent.mkdir()
    stub.write_text(STUB.format(python=sys.executable))
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    fasta = tmp_path / "g.fa"
    fasta.write_text(">c\nACGT\n")
    result = bu.busco_assess_completeness(
        str(fasta), "eukaryota_odb10", str(tmp_path / "out"), gene_predictor="miniprot", scaffold_composition=True, offline=True, cpus=2,
        busco_path=str(stub), verbose=False,
    )
    args = json.loads((stub.parent / "args.json").read_text())
    assert "--miniprot" in args and "--scaffold_composition" in args and "--offline" in args and args[args.index("-l") + 1] == "eukaryota_odb10"
    assert result["scores_percent"]["complete"] == 98.4 and result["marker_count"] == 124
    assert result["assembly_statistics"]["total_length"] == 4641652 and result["assembly_statistics"]["percent_gaps"] == 0.0
    assert result["counts"]["complete_buscos_c"] == 122
    assert result["full_table"][0]["Busco id"] == "b1" and result["scaffold_composition"][0]["Scaffold"] == "chr"
