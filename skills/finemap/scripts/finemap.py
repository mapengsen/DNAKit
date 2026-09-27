#!/usr/bin/env python3
"""Fine-mapping wrapper around mancusolab/sushie (+ SuSiE-inf, SuSiE-ash / SuSiE 2.0).

Subcommands:
    finemap.py susie     --vcf <one.vcf> --pheno <one.pheno> [--covar ...] --out <prefix>
    finemap.py sushie    --vcf <k.vcf ...> --pheno <k.pheno ...> [--covar ...] --out <prefix>
    finemap.py sumstats  --z <k.z.tsv ...> --ld <k.ld ...> --n <n1 n2 ...> --out <prefix>
                         [--method {sushie,susie_inf,auto}]
                         [--susie-inf-impl {susie2,finucane}]
    finemap.py recommend --okg-dataset-id <node> [--n-inputs N]
    finemap.py --verify-install   # run the bundled 3-ancestry tutorial

On first use the script clones https://github.com/mancusolab/sushie to
~/.cache/sushie/repo and pip-installs it, then runs the bundled tutorial
under <repo>/data/ as a smoke test. The first time --method susie_inf is
used the script installs the chosen SuSiE-inf implementation — by default
StatFunGen/susieR (SuSiE 2.0; Wang 2025) via R+remotes into
~/.cache/finemap/r_library/. With `--susie-inf-impl finucane` it instead
clones FinucaneLab/fine-mapping-inf (Cui 2024) into the same py310 venv
sushie uses. Records OKG provenance to a .finemap.json sidecar when
$OKG_REPO is set; the manifest carries `susie_inf_impl` so you can tell
the two apart downstream.

The `sumstats` mode dispatches by --method:
  - `sushie` (default) — sushie's infer_sushie_ss; K=1 SuSiE-RSS-equivalent
    or K>=2 multi-ancestry SuShiE.
  - `susie_inf` — K=1 with an infinitesimal polygenic background. Two
    impls available via `--susie-inf-impl`:
      * `susie2` (default) — StatFunGen/susieR (SuSiE 2.0; cleaner
        convergence on real GWAS). Cites paper:susie_ash_2025.
      * `finucane` — FinucaneLab/fine-mapping-inf (Cui 2024 Python; kept
        for regression-testing). Cites paper:susie_inf_2024.
  - `auto` — query the OKG (via --okg-dataset-id) and pick one. See
    `recommend` for the standalone decision-tree view.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


CACHE_ROOT = Path.home() / ".cache" / "sushie"
FINEMAP_CACHE = Path.home() / ".cache" / "finemap"
DEFAULT_REPO = "https://github.com/mancusolab/sushie"
SUSIE_INF_REPO = "https://github.com/FinucaneLab/fine-mapping-inf"
SUSIE_ASH_REPO = "https://github.com/StatFunGen/susieR"

OKG_NODES_SUSHIE = {
    "method": "method:sushie",
    "software": "software:sushie",
    "paper": "paper:sushie_2025",
}
OKG_NODES_SUSIE = {
    "method": "method:susie_finemapping",
    "software_operational": "software:sushie",
    "software_canonical": "software:susie",
    "paper_canonical": "paper:susie_2020",
}
# SuSiE-inf model has two implementations available in this skill:
#   `susie2`    = StatFunGen/susieR (SuSiE 2.0; Wang 2025) — DEFAULT.
#                Cleaner convergence on real GWAS; uses `susie_rss(...,
#                unmappable_effects = "inf")`.
#   `finucane`  = FinucaneLab/fine-mapping-inf (Cui 2024) — Python; kept
#                for continuity / regression-testing against the original
#                impl. Selected via `--susie-inf-impl finucane`.
# Both implement the same model (`method:susie_inf` in OKG); software +
# paper citations differ accordingly.
OKG_NODES_SUSIE_INF_FINUCANE = {
    "method": "method:susie_inf",
    "software": "software:susie_inf",
    "paper": "paper:susie_inf_2024",
}
OKG_NODES_SUSIE_INF_SUSIE2 = {
    "method": "method:susie_inf",
    "software": "software:susie_ash",
    "paper": "paper:susie_ash_2025",
}
# Default impl citation; the recommender + dispatch use SUSIE2 unless
# the user opts in to `--susie-inf-impl finucane`.
OKG_NODES_SUSIE_INF = OKG_NODES_SUSIE_INF_SUSIE2

# Methods the `sumstats` subcommand can dispatch to. `sushie` covers both
# K=1 (SuSiE-RSS-equivalent) and K>=2 (multi-ancestry SuShiE); `susie_inf`
# dispatches to one of two impls via `--susie-inf-impl`.
SUMSTATS_METHODS = ("sushie", "susie_inf", "auto")
SUSIE_INF_IMPLS = ("susie2", "finucane")


# ---------------------------- Install ----------------------------

def ensure_repo(repo_url: str, commit: Optional[str], cache_root: Path,
                refresh: bool = False) -> Path:
    repo_dir = cache_root / "repo"
    if refresh and repo_dir.exists():
        shutil.rmtree(repo_dir)
    if not repo_dir.exists():
        cache_root.mkdir(parents=True, exist_ok=True)
        print(f"cloning {repo_url} -> {repo_dir}", file=sys.stderr)
        subprocess.check_call(
            ["git", "clone", "--depth", "1", repo_url, str(repo_dir)],
            stdout=sys.stderr,
        )
    if commit:
        subprocess.check_call(
            ["git", "fetch", "--unshallow"], cwd=repo_dir,
            stdout=sys.stderr, stderr=sys.stderr,
        )
        subprocess.check_call(
            ["git", "checkout", commit], cwd=repo_dir, stdout=sys.stderr,
        )
    return repo_dir


def ensure_sushie_installed(repo_dir: Path) -> Path:
    """Install sushie into an isolated `uv` venv when the host Python is
    incompatible with sushie's pinned deps. Returns the python binary to
    invoke (which `--vcf`/`--gwas`/etc. will then be passed to via
    `python -m sushie`). Reuses the venv on subsequent runs.

    sushie pins glimix_core==3.1.13 which is Python 3.7-3.10 only. On
    Python 3.11+ host we can't install into the host env directly, so
    we use uv to spin up a Python 3.10 venv keyed to the sushie cache.
    """
    venv_dir = repo_dir.parent / "venv-py310"
    venv_py = venv_dir / "bin" / "python"
    if venv_py.exists():
        # Confirm sushie loads.
        r = subprocess.run([str(venv_py), "-c", "import sushie"],
                            capture_output=True, text=True)
        if r.returncode == 0:
            return venv_py

    if not shutil.which("uv"):
        sys.exit("ERROR: `uv` is required to install sushie into a "
                 "Python 3.10 venv (sushie's deps don't support 3.11+). "
                 "Install with `brew install uv` or "
                 "`pipx install uv`.")
    print(f"[sushie] creating Python 3.10 venv at {venv_dir} ...",
          file=sys.stderr)
    r = subprocess.run(["uv", "venv", "--python", "3.10", str(venv_dir)],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: uv venv creation failed: {r.stderr}")
    print(f"[sushie] installing from {repo_dir} into {venv_dir} ...",
          file=sys.stderr)
    r = subprocess.run(["uv", "pip", "install", "--python", str(venv_py),
                         str(repo_dir)],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: uv pip install of sushie failed:\n{r.stderr}")

    # cyvcf2 ships a prebuilt macOS wheel linked against an older htslib;
    # on hosts where the local htslib differs we hit
    # `symbol not found in flat namespace '_bcf_float_missing'`. Force a
    # source rebuild of cyvcf2 against the local htslib install
    # (homebrew/conda) so the dynamic symbols match.
    r = subprocess.run(
        [str(venv_py), "-c",
         "from cyvcf2 import VCF"],
        capture_output=True, text=True)
    if r.returncode != 0 and "bcf_float_missing" in r.stderr:
        print("[sushie] cyvcf2 dylib mismatch detected; rebuilding from "
              "source against local htslib ...", file=sys.stderr)
        # cyvcf2's source build calls `autoreconf -i`; ensure it's there.
        if not shutil.which("autoreconf"):
            if sys.platform == "darwin" and shutil.which("brew"):
                print("[sushie] installing autoconf/automake/libtool ...",
                      file=sys.stderr)
                subprocess.run(["brew", "install", "autoconf", "automake",
                                  "libtool"], capture_output=True, text=True)
            elif shutil.which("apt-get"):
                subprocess.run(["sudo", "apt-get", "install", "-y",
                                  "autoconf", "automake", "libtool"],
                                capture_output=True, text=True)
        env = os.environ.copy()
        # Homebrew prefix on Apple Silicon vs Intel macs.
        prefix = "/opt/homebrew" if Path("/opt/homebrew").exists() \
            else "/usr/local"
        # cyvcf2 1.x ships its own htslib by default; tell it to use the
        # system install (homebrew or conda) instead so the linked
        # symbols match the dylib actually on the host.
        env["CYVCF2_HTSLIB_MODE"] = "EXTERNAL"
        env["HTSLIB_LIBRARY_DIR"] = f"{prefix}/lib"
        env["HTSLIB_INCLUDE_DIR"] = f"{prefix}/include"
        env["CFLAGS"] = (env.get("CFLAGS", "") +
                          f" -I{prefix}/include").strip()
        env["LDFLAGS"] = (env.get("LDFLAGS", "") +
                           f" -L{prefix}/lib").strip()
        # Sushie pins pandas==1.5.0 which doesn't support numpy>=2.0
        # ('numpy.dtype size changed' ABI break). Pin numpy<2 in the
        # cyvcf2 rebuild so the venv's pandas keeps working.
        r2 = subprocess.run(
            ["uv", "pip", "install", "--python", str(venv_py),
             "--force-reinstall", "--no-binary", "cyvcf2",
             "numpy<2", "cyvcf2"],
            capture_output=True, text=True, env=env)
        if r2.returncode != 0:
            sys.exit(f"ERROR: cyvcf2 source rebuild failed:\n{r2.stderr}")
    return venv_py


def get_repo_commit(repo_dir: Path) -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_dir,
                       capture_output=True, text=True)
    return r.stdout.strip()


def which_sushie_cli(venv_py: Optional[Path] = None) -> str:
    """Return the command (as a space-separated string for shell-equivalent
    invocation, split before subprocess.call) that launches sushie.

    Prefers (1) the `sushie` console-script in the venv's bin/ directory
    when an isolated venv was created, (2) a `sushie` on PATH, (3) the
    host `python -m sushie` as last resort (rarely works since sushie
    doesn't ship a __main__).
    """
    if venv_py is not None:
        venv_sushie = venv_py.parent / "sushie"
        if venv_sushie.exists():
            return str(venv_sushie)
    p = shutil.which("sushie")
    return p if p else f"{sys.executable} -m sushie"


# ---------------------------- Verify install ----------------------------

def verify_install(repo_dir: Path, venv_py: Optional[Path] = None) -> int:
    """Run the bundled 3-ancestry tutorial to verify the install works."""
    sentinel = CACHE_ROOT / "verify_install" / ".verified"
    data_dir = repo_dir / "data"
    out_dir = CACHE_ROOT / "verify_install"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "test_result"
    if not data_dir.exists():
        sys.exit(f"ERROR: bundled tutorial data not found at {data_dir}")
    cli = which_sushie_cli(venv_py)
    cmd = (cli.split() + [
        "finemap",
        "--pheno", "EUR.pheno", "AFR.pheno", "EAS.pheno",
        "--vcf",
        str(data_dir / "vcf" / "EUR.vcf"),
        str(data_dir / "vcf" / "AFR.vcf"),
        str(data_dir / "vcf" / "EAS.vcf"),
        "--covar", "EUR.covar", "AFR.covar", "EAS.covar",
        "--output", str(out_prefix),
    ])
    print(f"[verify-install] running 3-ancestry tutorial in {data_dir}",
          file=sys.stderr)
    print(f"  cmd: {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd, cwd=str(data_dir))
    if rc != 0:
        sys.exit(f"ERROR: sushie finemap tutorial failed (exit {rc}); "
                 f"see logs above. Sentinel NOT written.")
    # Check that the expected output exists.
    found_any = any(out_dir.glob("test_result*"))
    if not found_any:
        sys.exit(f"ERROR: sushie finemap returned 0 but no test_result* "
                 f"files were created under {out_dir}")
    sentinel.write_text(json.dumps({
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "tutorial": "3-ancestry EUR+AFR+EAS HapMap",
        "outputs": [str(p) for p in sorted(out_dir.glob("test_result*"))],
    }, indent=2))
    print(f"[verify-install] OK; sentinel -> {sentinel}", file=sys.stderr)
    return 0


def is_verified() -> bool:
    return (CACHE_ROOT / "verify_install" / ".verified").exists()


# ---------------------------- OKG resolver ----------------------------

def resolve_okg(okg_repo: Optional[Path], node_ids: list[str]) -> dict:
    if okg_repo is None:
        return {}
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
        return {}
    env = os.environ.copy()
    env.setdefault("OKG_DSN",
                   "postgres://postgres:okg@localhost:5449/statgen_analysis")
    proc = subprocess.Popen(
        ["uv", "run", "--extra", "mcp", "python",
         "deployments/statgen-analysis/server.py"],
        cwd=str(okg_repo), env=env,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, bufsize=1,
    )
    def send(m): proc.stdin.write(json.dumps(m) + "\n"); proc.stdin.flush()
    def read(): return json.loads(proc.stdout.readline())
    found = {}
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "finemap-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        for nid in node_ids:
            send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                  "params": {"name": "get_node",
                             "arguments": {"node_id": nid}}})
            resp = read()
            sc = resp.get("result", {}).get("structuredContent") or {}
            if sc.get("node_id") == nid:
                found[nid] = nid
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return found


# ---------------------------- Run sushie ----------------------------

def run_finemap(args, mode: str, repo_dir: Path,
                 venv_py: Optional[Path] = None) -> int:
    """Invoke `sushie finemap` with the user's flags. Mode is 'susie' or
    'sushie' — chosen by number of VCF args."""
    cli = which_sushie_cli(venv_py).split()
    cmd = cli + ["finemap",
                 "--vcf", *[str(p) for p in args.vcf],
                 "--pheno", *[str(p) for p in args.pheno],
                 "--output", str(args.out)]
    if args.covar:
        cmd += ["--covar", *[str(p) for p in args.covar]]
    if args.extra:
        cmd += args.extra
    print(f"[finemap {mode}] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)

    # OKG provenance.
    if mode == "sushie":
        nids = list(OKG_NODES_SUSHIE.values())
        found = resolve_okg(args.okg_repo, nids)
        okg_node_ids = {k: v for k, v in OKG_NODES_SUSHIE.items()
                         if v in found}
    else:
        nids = list(OKG_NODES_SUSIE.values())
        found = resolve_okg(args.okg_repo, nids)
        okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE.items()
                         if v in found}

    # Summary from any test_result.cs.tsv style output.
    summary = _parse_outputs(args.out)
    manifest = {
        "subcommand": mode,
        "output_prefix": str(args.out),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": len(args.vcf),
        "inputs": {
            "vcf":   [str(p) for p in args.vcf],
            "pheno": [str(p) for p in args.pheno],
            "covar": [str(p) for p in (args.covar or [])],
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(str(args.out) + ".finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    if rc == 0 and summary:
        print(f"summary: {summary.get('n_cs', '?')} CSs, "
              f"mean size {summary.get('mean_cs_size', '?')}, "
              f"max PIP {summary.get('max_pip', '?')}")
    return rc


def _find_sushie_output(out_prefix: Path, suffix: str) -> Optional[Path]:
    """Sushie writes `<out_prefix>.sushie.<suffix>` (e.g. .sushie.cs.tsv)
    for the `region` subcommand and `<out_prefix>.<suffix>` for the
    individual-level tutorial path. Return the existing path, or None.
    """
    for candidate in (Path(f"{out_prefix}.sushie.{suffix}"),
                       Path(f"{out_prefix}.{suffix}")):
        if candidate.exists():
            return candidate
    return None


def _parse_outputs(out_prefix: Path) -> dict:
    """Parse sushie's per-CS + weights tables. Handles both naming
    conventions: the `region` subcommand produces `<prefix>.sushie.cs.tsv`
    + `<prefix>.sushie.weights.tsv`; the individual-level tutorial
    produces `<prefix>.cs.tsv` + `<prefix>.weight.tsv`."""
    summary: dict = {}
    cs_path = _find_sushie_output(out_prefix, "cs.tsv")
    if cs_path is not None:
        try:
            with open(cs_path) as f:
                header = f.readline().rstrip("\n").split("\t")
                rows = [line.rstrip("\n").split("\t") for line in f]
            if rows and "CSIndex" in header:
                idx = header.index("CSIndex")
                cs_groups = {}
                for r in rows:
                    cs_groups.setdefault(r[idx], []).append(r)
                summary["n_cs"] = len(cs_groups)
                if cs_groups:
                    sizes = [len(v) for v in cs_groups.values()]
                    summary["mean_cs_size"] = round(sum(sizes) / len(sizes), 2)
                    summary["cs_sizes"] = sorted(sizes)
            if rows and "pip_all" in header:
                pip_idx = header.index("pip_all")
                pips = []
                for r in rows:
                    try:
                        pips.append(float(r[pip_idx]))
                    except (ValueError, IndexError):
                        continue
                if pips:
                    summary["max_pip"] = round(max(pips), 4)
        except Exception as e:
            summary["cs_parse_error"] = str(e)
    weight_path = _find_sushie_output(out_prefix, "weights.tsv") or \
                   _find_sushie_output(out_prefix, "weight.tsv")
    if weight_path is not None:
        try:
            with open(weight_path) as f:
                header = f.readline().rstrip("\n").split("\t")
                # The region subcommand uses `sushie_pip_all`; the
                # tutorial uses `PIP`. Prefer whichever exists.
                for pip_name in ("sushie_pip_all", "PIP", "pip_all"):
                    if pip_name in header:
                        pip_idx = header.index(pip_name)
                        break
                else:
                    pip_idx = None
                if pip_idx is not None:
                    max_pip = 0.0
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        try:
                            v = float(parts[pip_idx])
                            if v > max_pip:
                                max_pip = v
                        except (ValueError, IndexError):
                            continue
                    summary["max_pip"] = round(max_pip, 4)
        except Exception as e:
            summary["weight_parse_error"] = str(e)
    corr_path = _find_sushie_output(out_prefix, "corr.tsv")
    if corr_path is not None:
        summary["cross_ancestry_correlation_file"] = str(corr_path)
    return summary


# ---------------------------- htslib tool install + tabix slicing ----------------------------

_HTSLIB_BINS = ("tabix", "bgzip", "bcftools")
_FINEMAP_REFCACHE = Path.home() / ".cache" / "finemap"

# Per-chromosome 1000G phase3 VCF URL pattern (build hg19).
_1000G_HG19_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/"
                    "20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_"
                    "integrated_v5b.20130502.genotypes.vcf.gz")
# Build GRCh38 NYGC re-call:
_1000G_HG38_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
                    "data_collections/1000_genomes_project/release/"
                    "20190312_biallelic_SNV_and_INDEL/"
                    "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_"
                    "27022019.GRCh38.phased.vcf.gz")
# Sample → population integrated panel file.
_1000G_PANEL_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/"
                     "20130502/integrated_call_samples_v3.20130502."
                     "ALL.panel")

_POP_TO_SUPERPOP = {
    "eur": "EUR", "eas": "EAS", "afr": "AFR", "sas": "SAS", "amr": "AMR",
}


def ensure_htslib_tools() -> dict:
    """Make sure `tabix`, `bgzip`, `bcftools` are on PATH. If missing, try
    a platform-appropriate install (brew on macOS, apt on Debian/Ubuntu,
    conda if no system path works). Returns {bin: path} or sys.exits.
    """
    paths = {b: shutil.which(b) for b in _HTSLIB_BINS}
    missing = [b for b, p in paths.items() if not p]
    if not missing:
        return paths

    print(f"[finemap] installing missing htslib tools {missing} ...",
          file=sys.stderr)
    # Install htslib/bcftools + autoconf/automake/libtool — the latter
    # three are needed when cyvcf2 has to rebuild from source against
    # the local htslib (its CMake calls `autoreconf -i`).
    installers = []
    if sys.platform == "darwin" and shutil.which("brew"):
        installers.append(["brew", "install", "htslib", "bcftools",
                            "autoconf", "automake", "libtool"])
    if shutil.which("apt-get"):
        installers.append(["sudo", "apt-get", "install", "-y",
                            "tabix", "bcftools", "autoconf", "automake",
                            "libtool"])
    if shutil.which("conda"):
        installers.append(["conda", "install", "-y", "-c", "bioconda",
                            "htslib", "bcftools"])
    if shutil.which("mamba"):
        installers.append(["mamba", "install", "-y", "-c", "bioconda",
                            "htslib", "bcftools"])

    for cmd in installers:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                print(f"[finemap] installed via: {' '.join(cmd[:3])}",
                      file=sys.stderr)
                break
        except FileNotFoundError:
            continue

    paths = {b: shutil.which(b) for b in _HTSLIB_BINS}
    still_missing = [b for b, p in paths.items() if not p]
    if still_missing:
        sys.exit(f"ERROR: could not install {still_missing}. Install "
                 f"manually: `brew install htslib bcftools` (macOS) / "
                 f"`apt-get install tabix bcftools` (Debian) / "
                 f"`conda install -c bioconda htslib bcftools`.")
    return paths


def _fetch_1000g_panel(cache_root: Path) -> Path:
    """Cache the 1000G integrated panel file (sample → population)."""
    p = cache_root / "1000g_integrated_panel.tsv"
    if p.exists() and p.stat().st_size > 0:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    print(f"[finemap] downloading 1000G integrated panel -> {p}",
          file=sys.stderr)
    urllib.request.urlretrieve(_1000G_PANEL_URL, p)
    return p


def _samples_for_superpop(panel_file: Path, superpop: str) -> list[str]:
    """Return the list of 1000G sample IDs in a given super-population."""
    samples = []
    with open(panel_file) as f:
        header = f.readline().rstrip("\n").split()
        try:
            i_sample = header.index("sample")
            i_super = header.index("super_pop")
        except ValueError:
            return []
        for line in f:
            r = line.rstrip("\n").split()
            if len(r) > i_super and r[i_super] == superpop:
                samples.append(r[i_sample])
    return samples


def ensure_1000g_slice(chrom: int, start: int, end: int,
                        population: str, build: str = "hg19",
                        cache_root: Path = _FINEMAP_REFCACHE) -> Path:
    """Tabix-slice the 1000G phase3 VCF to (chrom:start-end), subset to a
    super-population, bgzip+index, and cache. Returns the local cached path.

    Args:
        chrom: 1-22 (numeric)
        start, end: bp coordinates on `build`
        population: eur, eas, afr, sas, amr (case-insensitive)
        build: hg19 (default) or hg38

    Reuses the cached file when present. Total fetch on first run is
    typically a few MB for a 1 Mb window.
    """
    pop = population.lower()
    superpop = _POP_TO_SUPERPOP.get(pop)
    if not superpop:
        sys.exit(f"REFUSED: unknown population {population!r}; "
                 f"allowed: {list(_POP_TO_SUPERPOP)}")
    if build not in ("hg19", "hg38"):
        sys.exit(f"REFUSED: unknown build {build!r}; allowed: hg19, hg38")

    ensure_htslib_tools()
    cache_root = cache_root / f"1000g_{build}_{pop}"
    cache_root.mkdir(parents=True, exist_ok=True)
    out_vcf = cache_root / f"chr{chrom}_{start}_{end}.vcf.gz"
    if out_vcf.exists() and out_vcf.stat().st_size > 0:
        idx = Path(str(out_vcf) + ".tbi")
        if idx.exists():
            return out_vcf

    url_template = _1000G_HG19_URL if build == "hg19" else _1000G_HG38_URL
    url = url_template.format(chrom=chrom)
    region = f"{chrom}:{start}-{end}"

    # 1) Cache the sample → population panel.
    panel = _fetch_1000g_panel(_FINEMAP_REFCACHE)
    samples = _samples_for_superpop(panel, superpop)
    if not samples:
        sys.exit(f"REFUSED: no {superpop} samples found in 1000G panel "
                 f"at {panel}")
    keep_file = cache_root / f"{pop}.samples.txt"
    keep_file.write_text("\n".join(samples) + "\n")

    # 2) tabix-slice the remote VCF to the window, pipe into bcftools to
    # subset samples, write bgzipped output, then tabix-index.
    print(f"[finemap] tabix-slicing {url} @ {region} ({len(samples)} "
          f"{superpop} samples) -> {out_vcf}", file=sys.stderr)
    tmp_raw = cache_root / f"chr{chrom}_{start}_{end}.raw.vcf.gz"
    # Stream tabix output through bgzip directly so we never materialise the
    # full-population intermediate as plain VCF.
    with open(tmp_raw, "wb") as out:
        p_tab = subprocess.Popen(["tabix", "-h", url, region],
                                  stdout=subprocess.PIPE)
        p_bz = subprocess.Popen(["bgzip", "-c"],
                                 stdin=p_tab.stdout, stdout=out)
        p_tab.stdout.close()
        p_bz.communicate()
        p_tab.wait()
        if p_tab.returncode != 0 or p_bz.returncode != 0:
            sys.exit(f"ERROR: tabix/bgzip pipe failed for {region} from "
                     f"{url}")
    # bcftools view -S to subset samples + annotate ID = chr:pos. The
    # GWAS converter writes IDs in the same `chr:pos` form, so the join
    # is position-only — sushie's internal logic does allele alignment
    # from the genotype dosages vs the GWAS effect/other alleles.
    # Filter to biallelic SNPs (-m2 -M2 -v snps) to drop multi-allelic
    # sites that would create duplicate IDs.
    tmp_sub = out_vcf.with_suffix(".sub.vcf.gz")
    rc = subprocess.call(["bcftools", "view",
                           "-S", str(keep_file),
                           "-m2", "-M2", "-v", "snps",
                           "-Oz", "-o", str(tmp_sub),
                           str(tmp_raw)])
    if rc != 0:
        sys.exit(f"ERROR: bcftools view failed (rc={rc})")
    tmp_raw.unlink(missing_ok=True)
    rc = subprocess.call([
        "bcftools", "annotate",
        "--set-id", "%CHROM:%POS",
        "-Oz", "-o", str(out_vcf), str(tmp_sub)])
    if rc != 0:
        sys.exit(f"ERROR: bcftools annotate failed (rc={rc})")
    tmp_sub.unlink(missing_ok=True)

    rc = subprocess.call(["tabix", "-p", "vcf", str(out_vcf)])
    if rc != 0:
        sys.exit(f"ERROR: tabix index failed for {out_vcf}")
    return out_vcf


# ---------------------------- GWAS region mode (sumstats + ref VCF/PLINK) ----------------------------

# sushie's `finemap --summary` mode reads a sumstats TSV with this exact
# column order; we auto-convert from GWAS-Catalog harmonised or LDSC
# munged input.
SUSHIE_GWAS_HEADER = ["chrom", "snp", "pos", "a1", "a0", "z"]

_NA_TOKENS = {"", "NA", "nan", "NaN", "N/A", ".", "null", "None"}


def _open_text(path: Path):
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _detect_ref_snp_id_format(ref_plink: Optional[list] = None,
                                ref_vcf: Optional[list] = None,
                                ref_bgen: Optional[list] = None,
                                sample_lines: int = 500) -> str:
    """Peek at the reference panel's variant IDs to choose between rsid and
    chr:pos GWAS-side IDs. Returns 'rsid' if >50% of sampled IDs start with
    'rs<digits>', else 'chrpos'. Defaults to 'chrpos' on any failure — it's
    the safer fallback because chr:pos is constructible from any reference
    whereas rsids aren't always populated."""
    import re
    sample = []
    if ref_plink:
        bim = Path(str(ref_plink[0]) + ".bim")
        if bim.exists():
            with open(bim) as f:
                for i, line in enumerate(f):
                    if i >= sample_lines: break
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        sample.append(parts[1].strip())
    elif ref_vcf:
        try:
            r = subprocess.run(
                ["bcftools", "view", "-H", "--no-header", str(ref_vcf[0])],
                capture_output=True, text=True, timeout=30, check=False,
            )
            for line in r.stdout.splitlines()[:sample_lines]:
                parts = line.split("\t")
                if len(parts) >= 3:
                    sample.append(parts[2].strip())
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    if not sample:
        return "chrpos"
    rsid_re = re.compile(r"^rs\d+$")
    n_rs = sum(1 for s in sample if rsid_re.match(s))
    return "rsid" if n_rs / len(sample) > 0.5 else "chrpos"


def detect_and_convert_to_sushie_gwas(input_path: Path,
                                       chrom: int, start: int, end: int,
                                       out_path: Path,
                                       extras_path: Optional[Path] = None,
                                       id_format: str = "chrpos"
                                       ) -> int:
    """Read a GWAS sumstats file (GWAS Catalog harmonised TSV, LDSC munged
    `.sumstats.gz`, or any TSV with chrom/snp/pos/effect_allele/other_allele
    + signed effect or z columns), filter to the (chrom, start, end) window,
    and write a sushie-compatible TSV with header
    `chrom snp pos a1 a0 z`.

    Returns the number of rows written.
    """
    with _open_text(input_path) as f:
        header = f.readline().rstrip("\n").split("\t")
    cols = {c: i for i, c in enumerate(header)}
    cols_l = {c.lower(): i for i, c in enumerate(header)}

    def col(*candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
            if c.lower() in cols_l:
                return cols_l[c.lower()]
        return None

    # Prefer harmonised columns; alleles get uppercased on the way out.
    chr_i = col("CHR", "chr", "chromosome", "hm_chrom")
    pos_i = col("POS", "pos", "base_pair_location", "hm_pos", "bp")
    snp_i = col("SNP", "snp", "hm_rsid", "rsid", "variant_id", "snpid")
    a1_i  = col("A1", "a1", "hm_effect_allele", "effect_allele")
    a0_i  = col("A2", "a0", "a2", "hm_other_allele", "other_allele")
    beta_i = col("BETA", "beta", "hm_beta", "b")
    se_i   = col("SE", "se", "standard_error", "stderr")
    z_i    = col("Z", "z", "zscore")
    or_i   = col("OR", "hm_odds_ratio", "odds_ratio")
    p_i    = col("P", "p", "p_value", "pval", "Pvalue", "PValue")

    missing = [name for name, idx in
                [("chrom", chr_i), ("pos", pos_i), ("snp", snp_i),
                 ("a1", a1_i), ("a0", a0_i)] if idx is None]
    if missing:
        sys.exit(f"REFUSED: input header missing required columns "
                 f"{missing}; header was {header[:25]}")
    has_signed = z_i is not None or (beta_i is not None and se_i is not None) \
                  or (or_i is not None and se_i is not None)
    if not has_signed:
        sys.exit("REFUSED: cannot derive a signed effect from the input "
                 "(need Z, or BETA+SE, or OR+SE)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    # Sidecar `<out>.gwas_extras.tsv` carries the per-SNP rsid + beta + se
    # + p mapped by chr:pos, so the post-run CS table can join back to
    # rsids and effect-size details (sushie's CS output uses chr:pos as
    # the SNP id and doesn't preserve rsids).
    extras_handle = None
    if extras_path is not None:
        extras_path.parent.mkdir(parents=True, exist_ok=True)
        extras_handle = open(extras_path, "w")
        extras_handle.write("chrpos\trsid\tbeta\tse\tp\n")
    with _open_text(input_path) as f, open(out_path, "w") as g:
        f.readline()  # skip header
        g.write("\t".join(SUSHIE_GWAS_HEADER) + "\n")
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) < max(i for i in (chr_i, pos_i, snp_i, a1_i, a0_i,
                                          beta_i, se_i, z_i, or_i)
                              if i is not None) + 1:
                continue
            # Normalise chrom to int; skip X/Y/MT (sushie only does 1-22).
            chrom_s = str(r[chr_i]).lstrip("chr").strip()
            if not chrom_s.isdigit():
                continue
            try:
                pos_int = int(r[pos_i])
            except ValueError:
                continue
            chrom_i = int(chrom_s)
            if chrom_i != int(chrom) or pos_int < int(start) or pos_int > int(end):
                continue
            snp = r[snp_i]
            a1 = r[a1_i].upper(); a0 = r[a0_i].upper()
            if snp in _NA_TOKENS or a1 in _NA_TOKENS or a0 in _NA_TOKENS:
                continue
            # signed effect -> z
            try:
                if z_i is not None and r[z_i] not in _NA_TOKENS:
                    z = float(r[z_i])
                elif beta_i is not None and se_i is not None:
                    b = float(r[beta_i]); se = float(r[se_i])
                    if se == 0:
                        continue
                    z = b / se
                elif or_i is not None and se_i is not None:
                    import math
                    b = math.log(float(r[or_i])); se = float(r[se_i])
                    if se == 0:
                        continue
                    z = b / se
                else:
                    continue
            except (ValueError, IndexError):
                continue
            # SNP-ID matching against the reference panel — the join is by
            # exact string equality, so the GWAS-side ID must use the same
            # convention as the panel's .bim / VCF ID column.
            #   id_format='rsid':    use the harmonised rsid when present
            #                        (matches HM3 / 1000G-with-rsids panels).
            #   id_format='chrpos':  use 'chr:pos' (matches 1000G Phase3 VCFs
            #                        shipped with ID='.', or panels that
            #                        already use chr:pos IDs).
            # The chr:pos fallback also kicks in when id_format='rsid' but
            # the GWAS row's rsid is missing (rare in harmonised files but
            # possible at the edges).
            if id_format == "rsid" and snp not in _NA_TOKENS:
                snp_id = snp
            else:
                snp_id = f"{chrom_i}:{pos_int}"
            g.write(f"{chrom_i}\t{snp_id}\t{pos_int}\t{a1}\t{a0}\t{z:.6g}\n")
            if extras_handle is not None:
                # Capture beta + se + p if present, else NA — used by the
                # post-run summary table to reconstruct effect-size detail.
                def _safe(idx):
                    if idx is None:
                        return "NA"
                    try:
                        v = r[idx]
                        return v if v not in _NA_TOKENS else "NA"
                    except IndexError:
                        return "NA"
                p_str = _safe(p_i)
                if beta_i is not None and se_i is not None:
                    b_str = _safe(beta_i); se_str = _safe(se_i)
                elif or_i is not None and se_i is not None:
                    try:
                        b_str = f"{math.log(float(r[or_i])):.6g}"
                    except (ValueError, IndexError):
                        b_str = "NA"
                    se_str = _safe(se_i)
                else:
                    b_str = "NA"; se_str = _safe(se_i)
                extras_handle.write(f"{snp_id}\t{snp}\t{b_str}\t{se_str}\t{p_str}\n")
            n_written += 1
    if extras_handle is not None:
        extras_handle.close()
    return n_written


def _compute_aligned_ld_from_plink(sushie_gwas_path: Path,
                                     plink_prefix: Path,
                                     chrom: int, start: int, end: int,
                                     out_z: Path, out_ld: Path,
                                     venv_py: Path) -> int:
    """Use pandas_plink (already in the py310 venv) to:
        1. Read the GWAS-formatted TSV (chrom snp pos a1 a0 z).
        2. Read the PLINK ref panel, slice to the (chrom, start, end) window.
        3. Intersect SNPs by ID; allele-align z so the effect allele matches
           the PLINK A1; compute LD = corr(genotypes).
        4. Write `<out_z>` (SNP + Z) and `<out_ld>` (.npy).

    Returns the number of SNPs retained. This is the same LD-computation
    step sushie does internally for its `region --summary` mode; we lift
    it out so the susie_inf backends (susie2 / finucane) can share the
    same aligned matrix.
    """
    helper = r'''
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin
ap = argparse.ArgumentParser()
ap.add_argument("--gwas", required=True); ap.add_argument("--plink", required=True)
ap.add_argument("--chrom", type=int, required=True)
ap.add_argument("--start", type=int, required=True); ap.add_argument("--end", type=int, required=True)
ap.add_argument("--out-z", required=True); ap.add_argument("--out-ld", required=True)
a = ap.parse_args()

gwas = pd.read_csv(a.gwas, sep="\t")
gwas = gwas[(gwas["chrom"].astype(str) == str(a.chrom))
             & (gwas["pos"] >= a.start) & (gwas["pos"] <= a.end)]
print(f"[ld-helper] GWAS rows in window: {len(gwas)}", file=sys.stderr)

G = read_plink1_bin(a.plink + ".bed", verbose=False)
keep = (G.chrom.values.astype(str) == str(a.chrom)) & \
        (G.pos.values >= a.start) & (G.pos.values <= a.end)
G = G.isel(variant=np.where(keep)[0])
print(f"[ld-helper] PLINK variants in window: {G.shape[1]}", file=sys.stderr)

# Intersect by rsid (PLINK .bim 'snp' column).
common = set(gwas["snp"]) & set(G.snp.values)
if not common:
    sys.exit("ERROR: no overlap between GWAS rsids and PLINK panel")
print(f"[ld-helper] common SNPs: {len(common)}", file=sys.stderr)

gwas = gwas[gwas["snp"].isin(common)].set_index("snp")
G = G.isel(variant=np.where(np.isin(G.snp.values, list(common)))[0])

# Reorder both to the PLINK order (preserves chr/pos sorting).
snp_order = list(G.snp.values)
gwas = gwas.loc[snp_order]

# Allele-align z to PLINK A1; if GWAS A1 differs from PLINK A1, flip sign.
plink_a1 = np.array([str(x).upper() for x in G.a1.values])
gwas_a1  = np.array([str(x).upper() for x in gwas["a1"].values])
flip = plink_a1 != gwas_a1
n_flip = int(flip.sum())
print(f"[ld-helper] sign-flipped {n_flip} of {len(snp_order)} z-scores to "
       f"PLINK A1", file=sys.stderr)
z = gwas["z"].values.astype(float).copy()
z[flip] = -z[flip]

# Compute LD = standardised genotype correlation, with monomorphic-SNP
# drop and diagonal clamp so the matrix is properly conditioned for
# SuSiE-RSS / SuSiE-inf / SuSiE-ash (all need R[i,i] = 1 and R psd).
g = np.asarray(G.values, dtype=np.float64)  # samples x variants
col_mean = np.nanmean(g, axis=0)
nan_mask = np.isnan(g)
if nan_mask.any():
    g[nan_mask] = np.take(col_mean, np.where(nan_mask)[1])
sd = g.std(axis=0)
keep_mask = sd > 1e-6
n_drop = int((~keep_mask).sum())
if n_drop > 0:
    print(f"[ld-helper] dropping {n_drop} monomorphic SNPs", file=sys.stderr)
    g = g[:, keep_mask]
    snp_order = [s for s, k in zip(snp_order, keep_mask) if k]
    z = z[keep_mask]
g -= g.mean(axis=0, keepdims=True)
g /= g.std(axis=0, keepdims=True)
ld = (g.T @ g) / g.shape[0]
np.fill_diagonal(ld, 1.0)  # ensure R[i,i] == 1 exactly
# Tiny ridge for numerical stability (PRS / fine-mapping convention).
ld = 0.999 * ld + 0.001 * np.eye(ld.shape[0])
np.save(a.out_ld, ld.astype(np.float64))
# Also write a .tsv companion (no header) so the R-side wrapper
# (StatFunGen/susieR; no numpy reader) can consume the same matrix.
tsv_path = a.out_ld.rsplit(".", 1)[0] + ".tsv"
np.savetxt(tsv_path, ld, fmt="%.6f", delimiter="\t")

pd.DataFrame({"SNP": snp_order, "Z": z}
              ).to_csv(a.out_z, sep="\t", index=False)
print(f"[ld-helper] wrote z ({len(z)} SNPs) -> {a.out_z}", file=sys.stderr)
print(f"[ld-helper] wrote LD ({ld.shape}) -> {a.out_ld}", file=sys.stderr)
print(len(snp_order))  # stdout: n_snps for the caller
'''
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False,
                                       dir=str(out_z.parent)) as hf:
        hf.write(helper)
        helper_path = hf.name
    cmd = [str(venv_py), helper_path,
           "--gwas", str(sushie_gwas_path),
           "--plink", str(plink_prefix),
           "--chrom", str(chrom),
           "--start", str(start), "--end", str(end),
           "--out-z", str(out_z), "--out-ld", str(out_ld)]
    print(f"[finemap region] computing aligned LD from PLINK panel via "
          f"pandas_plink: {plink_prefix}", file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True)
    try:
        os.unlink(helper_path)
    except OSError:
        pass
    if r.returncode != 0:
        sys.exit(f"ERROR: LD computation failed:\n{r.stderr[-1500:]}")
    # Echo helper's stderr for the user's diagnostics.
    if r.stderr:
        print(r.stderr.rstrip(), file=sys.stderr)
    try:
        return int(r.stdout.strip().split("\n")[-1])
    except ValueError:
        return 0


def run_region(args, repo_dir: Path,
                venv_py: Optional[Path] = None) -> int:
    """Single-locus fine-mapping with sushie's `finemap --summary` (default
    method) or, when `--method susie_inf` is given, via the SuSiE-inf
    model — using either StatFunGen/susieR (default, `--susie-inf-impl
    susie2`) or FinucaneLab/fine-mapping-inf (`--susie-inf-impl finucane`)
    after a shared pandas_plink-based LD-computation step against the
    local PLINK ref panel.
    """
    sushie_gwas = Path(str(args.out) + ".gwas.tsv")
    extras_sidecar = Path(str(args.out) + ".gwas_extras.tsv")
    # Pick the SNP-ID convention BEFORE writing the GWAS so the IDs match
    # the reference panel — getting this wrong is sushie's #1 "no common
    # SNPs" cause (see GH issue thread). Detect from the panel actually
    # passed in; fall back to chr:pos when no panel is available (sumstats
    # mode upstream).
    id_format = args.snp_id_format
    if id_format == "auto":
        id_format = _detect_ref_snp_id_format(
            ref_plink=args.ref_plink, ref_vcf=args.ref_vcf,
            ref_bgen=args.ref_bgen)
        print(f"[finemap region] auto-detected reference SNP-ID format: "
              f"{id_format} (override via --snp-id-format)", file=sys.stderr)
    n = detect_and_convert_to_sushie_gwas(
        args.gwas_sumstats, args.chrom, args.start, args.end, sushie_gwas,
        extras_path=extras_sidecar, id_format=id_format)
    print(f"[finemap region] converted GWAS -> {sushie_gwas} ({n} rows "
          f"in chr{args.chrom}:{args.start}-{args.end}, "
          f"SNP-ID={id_format})", file=sys.stderr)
    if n < 50:
        sys.exit(f"REFUSED: only {n} GWAS rows in the requested window; "
                 f"need at least 50 SNPs to fine-map")

    # ------ Method dispatch for region mode (--method susie_inf) ------
    method = getattr(args, "method", "sushie")
    if method == "auto":
        rec = recommend_method(args.okg_repo,
                                getattr(args, "okg_dataset_id", None),
                                n_inputs=1)
        method = rec["method"]
        print(f"[finemap region auto] recommended method: {method} "
              f"({rec['rationale']})", file=sys.stderr)
        args.method = method
        args.recommendation = rec
    if method == "susie_inf":
        if not args.ref_plink:
            sys.exit("REFUSED: --method susie_inf in region mode "
                     "currently requires --ref-plink <prefix> (the "
                     "pandas_plink-based LD-computation step needs PLINK "
                     "format). Convert your VCF to PLINK or use the "
                     "default --method sushie path which accepts VCF/BGEN.")
        # Compute aligned z + LD once from the local PLINK panel.
        out_z  = Path(str(args.out) + ".aligned_z.tsv")
        out_ld = Path(str(args.out) + ".aligned_ld.npy")
        m = _compute_aligned_ld_from_plink(
            sushie_gwas, args.ref_plink[0],
            args.chrom, args.start, args.end,
            out_z, out_ld, venv_py)
        if m < 50:
            sys.exit(f"REFUSED: only {m} SNPs after PLINK intersection; "
                     f"need at least 50")

        # Build an args-like for the sumstats runner.
        class _SumstatsArgs:
            pass
        sa = _SumstatsArgs()
        sa.z = [out_z]
        sa.ld = [out_ld]
        sa.n = [args.N]
        sa.out = args.out
        # L=5 is a conservative default for a fine-mapping window (~2-3k
        # SNPs, typically ≤2 independent causals). susie_inf's MoM
        # update and susie_rss both degenerate when many of the L slots
        # have no signal — lower L = fewer dead effects = stabler updates.
        sa.L = 5
        sa.max_iter = 500
        sa.threshold = 0.95
        sa.purity = 0.5
        sa.okg_repo = args.okg_repo
        sa.method = method
        sa.okg_dataset_id = getattr(args, "okg_dataset_id", None)

        # SuSiE-inf impl dispatch (default: susie2 / StatFunGen).
        impl = getattr(args, "susie_inf_impl", "susie2")
        if impl == "susie2":
            # The R wrapper reads LD from a whitespace TSV (no numpy);
            # the LD helper writes a `.tsv` companion alongside `.npy`.
            sa.ld = [Path(str(out_ld).rsplit(".", 1)[0] + ".tsv")]
            rscript, r_lib = ensure_susie_inf_susie2_installed()
            rc = run_susie_inf_susie2_sumstats(sa, rscript=rscript,
                                                 r_lib=r_lib)
        elif impl == "finucane":
            si_repo = ensure_susie_inf_installed(venv_py)
            rc = run_susie_inf_sumstats(sa, venv_py=venv_py,
                                          susie_inf_repo=si_repo)
        else:
            sys.exit(f"REFUSED: unknown --susie-inf-impl {impl!r}; "
                     f"expected one of {SUSIE_INF_IMPLS}")

        # Augment the just-written manifest with region-mode context.
        mpath = Path(f"{args.out}.finemap.json")
        if mpath.exists():
            try:
                m_data = json.loads(mpath.read_text())
                m_data["subcommand"] = "region"
                m_data["locus"] = {"chrom": int(args.chrom),
                                    "start": int(args.start),
                                    "end": int(args.end)}
                m_data["aligned_inputs"] = {"z": str(out_z),
                                             "ld": str(out_ld),
                                             "n_snps": m}
                m_data["inputs"]["gwas_sumstats"] = str(args.gwas_sumstats)
                m_data["inputs"]["ref_plink"] = [str(args.ref_plink[0])]
                if args.okg_dataset_id:
                    m_data.setdefault("okg_node_ids", {})["dataset"] = \
                        args.okg_dataset_id
                if args.okg_ld_panel_id:
                    m_data.setdefault("okg_node_ids", {})["ld_panel"] = \
                        args.okg_ld_panel_id
                mpath.write_text(json.dumps(m_data, indent=2))
            except Exception as e:
                print(f"[finemap region] manifest annotation failed: {e}",
                      file=sys.stderr)
        return rc

    cli = which_sushie_cli(venv_py).split()
    cmd = cli + ["finemap",
                  "--summary",
                  "--gwas", str(sushie_gwas),
                  "--gwas-header", *SUSHIE_GWAS_HEADER,
                  "--sample-size", str(args.N),
                  "--chrom", str(args.chrom),
                  "--start", str(args.start),
                  "--end", str(args.end),
                  "--output", str(args.out)]
    # LD source resolution (in priority order):
    #   1. Explicit --ref-vcf / --ref-plink / --ref-bgen / --ld
    #   2. --ref-1000g <pop>: tabix-slice from 1000G FTP + subset samples
    ref_vcf_paths = list(args.ref_vcf or [])
    if args.ref_1000g and not ref_vcf_paths and not args.ref_plink and \
            not args.ref_bgen and not args.ld:
        sliced = ensure_1000g_slice(args.chrom, args.start, args.end,
                                      population=args.ref_1000g,
                                      build=args.ref_1000g_build)
        ref_vcf_paths = [sliced]
        print(f"[finemap region] 1000G {args.ref_1000g.upper()} "
              f"slice -> {sliced}", file=sys.stderr)
    if ref_vcf_paths:
        cmd += ["--vcf", *[str(p) for p in ref_vcf_paths]]
    elif args.ref_plink:
        cmd += ["--plink", *[str(p) for p in args.ref_plink]]
    elif args.ref_bgen:
        cmd += ["--bgen", *[str(p) for p in args.ref_bgen]]
    elif args.ld:
        cmd += ["--ld", *[str(p) for p in args.ld]]
    else:
        sys.exit("REFUSED: provide one of --ref-vcf / --ref-plink / "
                 "--ref-bgen / --ld / --ref-1000g for the LD source")
    if args.extra:
        # Drop the literal `--` separator that argparse REMAINDER captures
        # but the underlying sushie CLI doesn't accept.
        cmd += [a for a in args.extra if a != "--"]

    print(f"[finemap region] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)

    # OKG provenance (single-ancestry SuSiE-RSS path).
    nids = list(OKG_NODES_SUSIE.values())
    found = resolve_okg(args.okg_repo, nids)
    okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE.items() if v in found}
    if args.okg_ld_panel_id:
        okg_node_ids["ld_panel"] = args.okg_ld_panel_id
    if args.okg_dataset_id:
        okg_node_ids["dataset"] = args.okg_dataset_id

    summary = _parse_outputs(args.out)
    summary["n_snps_in_window"] = n
    manifest = {
        "subcommand": "region",
        "output_prefix": str(args.out),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": 1,
        "locus": {"chrom": int(args.chrom),
                   "start": int(args.start),
                   "end": int(args.end)},
        "inputs": {
            "gwas_sumstats": str(args.gwas_sumstats),
            "ref_vcf":   [str(p) for p in (args.ref_vcf or [])],
            "ref_plink": [str(p) for p in (args.ref_plink or [])],
            "ref_bgen":  [str(p) for p in (args.ref_bgen or [])],
            "ld":        [str(p) for p in (args.ld or [])],
        },
        "sample_size": args.N,
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{args.out}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    if rc == 0 and summary:
        print(f"summary: {summary.get('n_cs', '?')} CSs, "
              f"mean size {summary.get('mean_cs_size', '?')}, "
              f"max PIP {summary.get('max_pip', '?')}")
        # Print the rsID-enriched credible-set table, optionally with
        # gnomAD + GTEx annotations appended.
        annot_sources: set = set()
        if args.annotate in ("gnomad", "both"):
            annot_sources.add("gnomad")
        if args.annotate in ("gtex", "both"):
            annot_sources.add("gtex")
        _print_cs_summary_table(Path(args.out), extras_sidecar,
                                  annot_sources=annot_sources)
    return rc


_ANNOT_CACHE = Path.home() / ".cache" / "variant-annotate"


def _load_variant_annotate():
    """Import variant-annotate/scripts/annotate.py as a module. The
    variant-annotate skill lives as a sibling directory under
    statgen-skills/. Returns the module, or None if not findable."""
    # finemap.py is at: statgen-skills/finemap/scripts/finemap.py
    # variant-annotate annotate.py is at:
    #         statgen-skills/variant-annotate/scripts/annotate.py
    here = Path(__file__).resolve().parent  # statgen-skills/finemap/scripts
    sibling = here.parent / "references" / "variant-annotate" / "scripts"
    if not (sibling / "annotate.py").exists():
        return None
    if str(sibling) not in sys.path:
        sys.path.insert(0, str(sibling))
    try:
        import annotate as va_mod  # noqa
        return va_mod
    except Exception as e:
        sys.stderr.write(f"[finemap annotate] failed to import "
                          f"variant-annotate: {type(e).__name__}: {e}\n")
        return None


def _annotate_variants(rsids: list, sources: set) -> dict:
    """Fetch gnomAD and/or GTEx annotations for rsIDs by delegating to
    the sibling `variant-annotate` skill. Returns
    {rsid: {gnomad: dict, gtex: dict}}, empty dict if the skill isn't
    available. Per-rsID JSON is cached by variant-annotate under
    ~/.cache/variant-annotate/."""
    va = _load_variant_annotate()
    if va is None:
        sys.stderr.write("[finemap annotate] variant-annotate skill not "
                          "found alongside finemap; skipping annotations\n")
        return {}
    records = va.annotate_rsids(list(rsids), set(sources), _ANNOT_CACHE)
    out: dict = {}
    for rec in records:
        rsid = rec.get("rsid")
        bundle: dict = {}
        if "gnomad" in sources and rec.get("gnomad"):
            bundle["gnomad"] = rec["gnomad"]
        if "gtex" in sources and rec.get("gtex"):
            bundle["gtex"] = rec["gtex"]
        if rsid and bundle:
            out[rsid] = bundle
    return out


def _print_cs_summary_table(out_prefix: Path,
                             extras_sidecar: Path,
                             annot_sources: Optional[set] = None) -> None:
    """Print a human-readable rsID-enriched credible-set table to stdout.

    Reads sushie's `<prefix>.sushie.cs.tsv` (which carries CSIndex, snp,
    pos, a0, a1, pip_all) and joins by chr:pos with the `<prefix>.
    gwas_extras.tsv` sidecar (chrpos, rsid, beta, se, p) the converter
    wrote. Output: one row per CS member, sorted by CS index then by
    descending pip_all within each CS."""
    cs_path = _find_sushie_output(out_prefix, "cs.tsv")
    if cs_path is None or not extras_sidecar.exists():
        return
    extras = {}
    with open(extras_sidecar) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_chrpos = header.index("chrpos")
        idx_rsid = header.index("rsid")
        idx_beta = header.index("beta")
        idx_se = header.index("se")
        idx_p = header.index("p")
        for line in f:
            r = line.rstrip("\n").split("\t")
            extras[r[idx_chrpos]] = (r[idx_rsid], r[idx_beta], r[idx_se],
                                       r[idx_p])
    rows = []
    with open(cs_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_cs = header.index("CSIndex")
            i_snp = header.index("snp")
            i_pos = header.index("pos")
            i_a0 = header.index("a0")
            i_a1 = header.index("a1")
            i_pip = header.index("pip_all")
        except ValueError:
            return
        for line in f:
            r = line.rstrip("\n").split("\t")
            cs = int(r[i_cs])
            try:
                pip = float(r[i_pip])
            except ValueError:
                pip = 0.0
            chrpos = r[i_snp]
            ex = extras.get(chrpos, ("NA", "NA", "NA", "NA"))
            rows.append((cs, -pip, chrpos, r[i_pos], r[i_a0], r[i_a1],
                          pip, *ex))
    if not rows:
        return
    rows.sort(key=lambda x: (x[0], x[1]))

    # Optional: fetch gnomAD + GTEx annotations for every CS-member rsid.
    annots = {}
    if annot_sources:
        rsids = sorted({r[7] for r in rows
                        if r[7] and r[7] != "NA" and r[7].startswith("rs")})
        if rsids:
            print(f"[finemap annotate] querying "
                  f"{','.join(sorted(annot_sources))} for {len(rsids)} "
                  f"variants (cached to {_ANNOT_CACHE}) ...",
                  file=sys.stderr)
            annots = _annotate_variants(rsids, annot_sources)

    # Render markdown-style table.
    print()
    print("Credible-set members (sorted by CS, descending PIP within CS):")
    print()
    base_cols = ["CS", "PIP", "rsID", "chr:pos", "Eff", "NonEff",
                  "Beta", "SE", "P"]
    extra_cols = []
    if "gnomad" in (annot_sources or set()):
        extra_cols += ["Csq", "Gene", "AF_nfe"]
    if "gtex" in (annot_sources or set()):
        extra_cols += ["top eQTL", "top sQTL"]
    cols = base_cols + extra_cols
    print("| " + " | ".join(cols) + " |")
    print("|" + "|".join("---" for _ in cols) + "|")
    for cs, _, chrpos, pos, a0, a1, pip, rsid, beta, se, p in rows:
        row = [str(cs), f"{pip:.4f}", rsid, chrpos, a1, a0,
                beta, se, p]
        rec = annots.get(rsid, {})
        if "gnomad" in (annot_sources or set()):
            g = rec.get("gnomad") or {}
            csq = (g.get("consequence") or "").replace("_variant", "")
            gene = g.get("gene") or ""
            af_nfe = g.get("af", {}).get("nfe")
            af_str = f"{af_nfe:.4g}" if isinstance(af_nfe, (int, float)) \
                else ""
            row += [csq, gene, af_str]
        if "gtex" in (annot_sources or set()):
            t = rec.get("gtex") or {}
            eq = t.get("top_eqtl") or {}
            sq = t.get("top_sqtl") or {}
            eq_str = (f"{eq.get('gene')}@{eq.get('tissue')} "
                       f"p={eq.get('p'):.1e}") if eq.get("gene") else ""
            sq_str = (f"{sq.get('gene')}@{sq.get('tissue')} "
                       f"p={sq.get('p'):.1e}") if sq.get("gene") else ""
            row += [eq_str, sq_str]
        print("| " + " | ".join(row) + " |")
    print()


# ---------------------------- Sumstats mode ----------------------------

def _load_z(path: Path) -> tuple[list[str], list[float]]:
    """Load a z-score file. Accept TSV/CSV/whitespace with a header. If a
    column named Z/z is present, use it. Else compute z = beta / se from
    BETA/beta + SE/se. Returns (snp_ids, z) — snp_ids may be index strings
    if no SNP/snp/rsid/SNPID column is found."""
    import csv
    with open(path) as f:
        first = f.readline()
    sep = "\t" if "\t" in first else (
        "," if "," in first and ";" not in first else None)
    rows = []
    with open(path) as f:
        if sep:
            reader = csv.reader(f, delimiter=sep)
        else:
            reader = ([x for x in line.split()] for line in f)
        header = next(iter(reader))
        for r in reader:
            if r and any(x.strip() for x in r):
                rows.append(r)
    cols = {c.lower(): i for i, c in enumerate(header)}
    snp_col = next((cols[k] for k in ("snp", "snpid", "rsid", "snp_id",
                                       "variant_id", "id") if k in cols), None)
    z_col = next((cols[k] for k in ("z", "zscore", "z_score") if k in cols),
                  None)
    if z_col is not None:
        snps = [r[snp_col] if snp_col is not None else str(i)
                for i, r in enumerate(rows)]
        zs = [float(r[z_col]) for r in rows]
        return snps, zs
    beta_col = next((cols[k] for k in ("beta", "effect_size", "b")
                      if k in cols), None)
    se_col = next((cols[k] for k in ("se", "standard_error", "stderr")
                    if k in cols), None)
    if beta_col is None or se_col is None:
        sys.exit(f"ERROR: {path} has no Z/z column and no BETA+SE columns "
                 f"to compute it; header was {header}")
    snps = [r[snp_col] if snp_col is not None else str(i)
            for i, r in enumerate(rows)]
    zs = [float(r[beta_col]) / float(r[se_col]) for r in rows]
    return snps, zs


def _load_ld(path: Path):
    """Load an LD matrix. Supports .npy (numpy save) and whitespace TSV."""
    import numpy as np
    p = Path(path)
    if p.suffix == ".npy":
        ld = np.load(p)
    else:
        ld = np.loadtxt(p)
    if ld.ndim != 2 or ld.shape[0] != ld.shape[1]:
        sys.exit(f"ERROR: LD at {p} is not a square matrix "
                 f"(shape={ld.shape})")
    return ld


def run_sumstats(args, repo_dir: Path) -> int:
    """Programmatic SuShiE/SuSiE sumstats fine-mapping via infer_sushie_ss."""
    import numpy as np
    try:
        from sushie.infer_ss import infer_sushie_ss
    except ImportError:
        sys.exit("ERROR: sushie package not importable; "
                 "run --verify-install first")

    K = len(args.z)
    if len(args.ld) != K:
        sys.exit(f"ERROR: --z count ({K}) != --ld count ({len(args.ld)})")
    if len(args.n) != K:
        sys.exit(f"ERROR: --z count ({K}) != --n count ({len(args.n)})")

    # Load + align by SNP if all files carry SNP IDs; else by row order.
    z_arrs, snp_lists, lds = [], [], []
    for i in range(K):
        snps_i, zs_i = _load_z(Path(args.z[i]))
        ld_i = _load_ld(Path(args.ld[i]))
        if ld_i.shape[0] != len(zs_i):
            sys.exit(f"ERROR: ancestry {i+1}: z file has {len(zs_i)} rows "
                     f"but LD is {ld_i.shape[0]}x{ld_i.shape[0]}")
        z_arrs.append(np.asarray(zs_i, dtype=np.float64))
        snp_lists.append(snps_i)
        lds.append(ld_i.astype(np.float64))

    # Sanity check: all ancestries must have same number of SNPs.
    ms = {len(z) for z in z_arrs}
    if len(ms) != 1:
        sys.exit(f"ERROR: ancestries have mismatched SNP counts: "
                 f"{[len(z) for z in z_arrs]}")
    m = next(iter(ms))

    # Run infer_sushie_ss.
    ns = np.asarray([int(x) for x in args.n], dtype=np.float64)
    print(f"[sumstats] K={K} ancestries, m={m} SNPs, L={args.L}, "
          f"max_iter={args.max_iter}", file=sys.stderr)
    result = infer_sushie_ss(
        lds=np.asarray(lds), ns=ns, zs=np.asarray(z_arrs),
        L=args.L, max_iter=args.max_iter, min_tol=args.min_tol,
        threshold=args.threshold, purity=args.purity,
        min_snps=max(10, args.min_snps),
    )

    # Write outputs in a finemap-friendly format.
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cs_path = Path(f"{out_prefix}.cs.tsv")
    weight_path = Path(f"{out_prefix}.weight.tsv")

    cs_df = result.cs  # pd.DataFrame with SNPIndex, CSIndex, alpha, c_alpha
    snps0 = snp_lists[0]
    if hasattr(cs_df, "to_csv"):
        # Resolve SNP IDs from ancestry 1 (all should match).
        cs_out = cs_df.copy()
        try:
            cs_out["SNP"] = [snps0[int(i)] for i in cs_out["SNPIndex"].values]
        except Exception:
            pass
        cs_out.to_csv(cs_path, sep="\t", index=False)
    pip_arr = np.asarray(result.pip)
    with open(weight_path, "w") as f:
        f.write("SNPIndex\tSNP\tPIP\n")
        for i in range(m):
            f.write(f"{i}\t{snps0[i]}\t{pip_arr[i]:.6g}\n")

    # OKG provenance.
    if K == 1:
        node_map = OKG_NODES_SUSIE
    else:
        node_map = OKG_NODES_SUSHIE
    found = resolve_okg(args.okg_repo, list(node_map.values()))
    okg_node_ids = {k: v for k, v in node_map.items() if v in found}

    summary = _parse_outputs(out_prefix)
    if "max_pip" not in summary:
        summary["max_pip"] = round(float(pip_arr.max()), 4)
    manifest = {
        "subcommand": "sumstats",
        "output_prefix": str(out_prefix),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": K,
        "n_snps": m,
        "inputs": {
            "z":  [str(p) for p in args.z],
            "ld": [str(p) for p in args.ld],
            "n":  [int(x) for x in args.n],
        },
        "infer_params": {
            "L": args.L, "max_iter": args.max_iter,
            "min_tol": args.min_tol, "threshold": args.threshold,
            "purity": args.purity, "min_snps": args.min_snps,
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{out_prefix}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    print(f"summary: {summary.get('n_cs', '?')} CSs, max PIP "
          f"{summary.get('max_pip', '?')}")
    return 0


# ---------------------------- SuSiE-inf (FinucaneLab/fine-mapping-inf) ----------------------------

def ensure_susie_inf_installed(venv_py: Path,
                                cache_root: Path = FINEMAP_CACHE,
                                refresh: bool = False) -> Path:
    """Clone FinucaneLab/fine-mapping-inf and install its two packages
    (`susieinf`, `finemapinf`) plus dependencies into the same py310 venv
    sushie already uses. Returns the path to the cloned repo (used to
    invoke `run_fine_mapping.py`). On subsequent calls, just verifies the
    import works.
    """
    repo_dir = cache_root / "fine-mapping-inf"
    if refresh and repo_dir.exists():
        shutil.rmtree(repo_dir)
    if not repo_dir.exists():
        cache_root.mkdir(parents=True, exist_ok=True)
        print(f"[susie_inf] cloning {SUSIE_INF_REPO} -> {repo_dir}",
              file=sys.stderr)
        subprocess.check_call(
            ["git", "clone", "--depth", "1", SUSIE_INF_REPO, str(repo_dir)],
            stdout=sys.stderr)

    # Smoke-import; install if needed.
    r = subprocess.run([str(venv_py), "-c",
                        "import susieinf, finemapinf"],
                        capture_output=True, text=True)
    if r.returncode == 0:
        return repo_dir
    print(f"[susie_inf] installing into {venv_py.parent.parent} ...",
          file=sys.stderr)
    # The repo has two setup.py packages under subdirs `susieinf/` and
    # `finemapinf/`. finemapinf's setup.py does `import numpy` without
    # declaring numpy as a build dep, and susieinf's setup.py needs
    # setuptools available at build time. The py310 venv already has
    # numpy (sushie depends on it). Pre-seed setuptools + pip so the
    # legacy `import setuptools` lines work under --no-build-isolation.
    r = subprocess.run(
        ["uv", "pip", "install", "--python", str(venv_py),
         "setuptools", "wheel", "pip"],
        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: failed to seed setuptools/wheel/pip into "
                 f"{venv_py}:\n{r.stderr[-500:]}")
    for pkg in ("susieinf", "finemapinf"):
        cmd = ["uv", "pip", "install", "--python", str(venv_py),
               "--no-build-isolation", str(repo_dir / pkg)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            # Try pip directly via the venv's python.
            print(f"[susie_inf] uv pip install of {pkg} failed; "
                  f"falling back to pip ...", file=sys.stderr)
            r2 = subprocess.run(
                [str(venv_py), "-m", "pip", "install",
                 "--no-build-isolation", str(repo_dir / pkg)],
                capture_output=True, text=True)
            if r2.returncode != 0:
                sys.exit(f"ERROR: install of {pkg} failed via both uv and "
                         f"pip:\nuv stderr={r.stderr[-500:]}\n"
                         f"pip stderr={r2.stderr[-500:]}")
    r = subprocess.run(
        ["uv", "pip", "install", "--python", str(venv_py), "bgzip"],
        capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[susie_inf] warning: bgzip install failed: {r.stderr}",
              file=sys.stderr)
    # Final verification.
    r = subprocess.run([str(venv_py), "-c",
                        "import susieinf, finemapinf; print('ok')"],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: susie_inf import still failing after install:\n"
                 f"{r.stderr}")
    return repo_dir


def run_susie_inf_sumstats(args, venv_py: Path,
                            susie_inf_repo: Path) -> int:
    """Run FinucaneLab/fine-mapping-inf on K=1 summary statistics.

    Calls the package's Python API directly (rather than the shell-out
    via run_fine_mapping.py) so we can hand it numpy arrays from the
    same z/ld/n inputs the sushie sumstats path accepts.
    """
    if len(args.z) != 1:
        sys.exit("ERROR: --method susie_inf requires K=1 "
                 "(single-ancestry); got {} --z files".format(len(args.z)))
    if len(args.ld) != 1 or len(args.n) != 1:
        sys.exit("ERROR: --method susie_inf requires exactly one "
                 "--ld and one --n value")
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # Run via subprocess so we use the py310 venv that has susieinf
    # installed; pass z, ld, n via temp .npy files for speed.
    import tempfile
    runner = textwrap_dedent_runner_susie_inf()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False,
                                       dir=str(out_prefix.parent)) as rf:
        rf.write(runner)
        runner_path = rf.name
    cmd = [
        str(venv_py), runner_path,
        "--z", str(args.z[0]),
        "--ld", str(args.ld[0]),
        "--n", str(args.n[0]),
        "--out", str(out_prefix),
        "--L", str(args.L),
        "--coverage", str(args.threshold),
        "--purity", str(args.purity),
    ]
    print(f"[finemap susie_inf] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)
    try:
        os.unlink(runner_path)
    except OSError:
        pass
    if rc != 0:
        print(f"[finemap susie_inf] runner exit code {rc}", file=sys.stderr)
        return rc

    # Build manifest.
    found = resolve_okg(args.okg_repo,
                         list(OKG_NODES_SUSIE_INF_FINUCANE.values()))
    okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE_INF_FINUCANE.items()
                     if v in found}
    summary = _parse_outputs(out_prefix)
    manifest = {
        "subcommand": "sumstats",
        "method": "susie_inf",
        "susie_inf_impl": "finucane",
        "output_prefix": str(out_prefix),
        "susie_inf_repo": SUSIE_INF_REPO,
        "susie_inf_commit": get_repo_commit(susie_inf_repo),
        "n_ancestries": 1,
        "inputs": {
            "z":  [str(args.z[0])],
            "ld": [str(args.ld[0])],
            "n":  [int(args.n[0])],
        },
        "infer_params": {
            "L": args.L, "max_iter": args.max_iter,
            "coverage": args.threshold, "purity": args.purity,
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{out_prefix}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    print(f"summary: {summary.get('n_cs', '?')} CSs, max PIP "
          f"{summary.get('max_pip', '?')}")
    return 0


def textwrap_dedent_runner_susie_inf() -> str:
    """Return a self-contained Python script that runs susie_inf on the
    given z/ld/n inputs and writes <prefix>.cs.tsv + <prefix>.weights.tsv
    in the same schema the sushie sumstats path uses."""
    return r'''
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from susieinf import susie as susieinf_susie, cred as susieinf_cred

def _load_z(path):
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}
    snp_c = next((cols[k] for k in ("snp","snpid","rsid","snp_id","variant_id","id") if k in cols), None)
    z_c = next((cols[k] for k in ("z","zscore","z_score") if k in cols), None)
    if z_c is not None:
        z = df[z_c].astype(float).values
    else:
        bc = next((cols[k] for k in ("beta","effect_size","b","effect") if k in cols), None)
        sc = next((cols[k] for k in ("se","standard_error","stderr","sderr") if k in cols), None)
        if bc is None or sc is None:
            sys.exit(f"ERROR: {path} has no Z col and no BETA+SE cols (header={list(df.columns)})")
        z = (df[bc].astype(float).values / df[sc].astype(float).values)
    snps = df[snp_c].astype(str).values if snp_c else np.array([str(i) for i in range(len(z))])
    return snps, z

def _load_ld(path):
    p = Path(path)
    if p.suffix == ".npy":
        ld = np.load(p)
    else:
        ld = np.loadtxt(p)
    if ld.ndim != 2 or ld.shape[0] != ld.shape[1]:
        sys.exit(f"ERROR: LD at {p} is not square (shape={ld.shape})")
    return ld

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", required=True); ap.add_argument("--ld", required=True)
    ap.add_argument("--n", type=int, required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--L", type=int, default=10); ap.add_argument("--coverage", type=float, default=0.95)
    ap.add_argument("--purity", type=float, default=0.5)
    a = ap.parse_args()
    snps, z = _load_z(a.z)
    LD = _load_ld(a.ld)
    if LD.shape[0] != len(z):
        sys.exit(f"ERROR: LD shape {LD.shape} != z length {len(z)}")
    # MLE empirical-Bayes (bounded optimisation) is more numerically
    # robust than MoM on real GWAS data — MoM solves a 2x2 system that
    # can NaN out when some L-slot effects fail to converge. The test
    # fixtures in FinucaneLab/fine-mapping-inf/susieinf/test_susieinf.py
    # demonstrate MLE as the recommended path.
    out = susieinf_susie(z=z, meansq=1.0, n=a.n, L=a.L, LD=LD,
                          est_sigmasq=True, est_tausq=True, method="MLE")
    cs = susieinf_cred(out["PIP"], coverage=a.coverage, purity=a.purity, LD=LD, n=a.n)
    pip = 1 - np.prod(1 - out["PIP"], axis=1)
    # weights.tsv: SNPIndex, SNP, PIP (matches sushie sumstats path)
    pd.DataFrame({"SNPIndex": np.arange(len(snps)), "SNP": snps, "PIP": pip}
                  ).to_csv(f"{a.out}.weights.tsv", sep="\t", index=False)
    # cs.tsv: CSIndex, SNPIndex, SNP, pip_all
    rows = []
    if cs is not None:
        for ci, members in enumerate(cs):
            for idx in members:
                idx = int(idx)
                rows.append({"CSIndex": ci+1, "SNPIndex": idx, "SNP": snps[idx], "pip_all": float(pip[idx])})
    cs_df = pd.DataFrame(rows or [], columns=["CSIndex","SNPIndex","SNP","pip_all"])
    cs_df.to_csv(f"{a.out}.cs.tsv", sep="\t", index=False)
    # log
    Path(f"{a.out}.log").write_text(
        f"susie_inf done\\nN={a.n} L={a.L}\\nm={len(z)} K_cs={cs_df['CSIndex'].nunique() if len(cs_df) else 0}\\n"
        f"sigmasq={out.get('sigmasq', 'na')} tausq={out.get('tausq', 'na')}\\n"
        f"max_pip={float(pip.max()):.6g}\\n"
    )
    print(f"susie_inf: m={len(z)} N={a.n} CSs={cs_df['CSIndex'].nunique() if len(cs_df) else 0} max_PIP={float(pip.max()):.4f}", file=sys.stderr)

if __name__ == "__main__":
    main()
'''


# ---------------------------- SuSiE-ash / SuSiE 2.0 (StatFunGen/susieR via R) ----------------------------

def ensure_r_installed() -> str:
    """Return path to Rscript, installing R if missing (brew on macOS,
    apt on Debian/Ubuntu). sys.exit on failure."""
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
    print("[susie_inf:susie2] R not found; installing ...",
          file=sys.stderr)
    if sys.platform == "darwin" and shutil.which("brew"):
        r = subprocess.run(["brew", "install", "r"], capture_output=True,
                            text=True)
        if r.returncode != 0:
            print(f"[susie_inf:susie2] brew install r warnings: "
                  f"{r.stderr[:500]}", file=sys.stderr)
    elif shutil.which("apt-get"):
        subprocess.run(["sudo", "apt-get", "install", "-y",
                          "r-base", "r-base-dev"],
                        capture_output=True, text=True)
    rscript = shutil.which("Rscript")
    if not rscript:
        sys.exit("ERROR: could not install R. Install manually: "
                 "`brew install r` (macOS) or `apt-get install r-base` "
                 "(Debian).")
    return rscript


def ensure_susie_inf_susie2_installed(cache_root: Path = FINEMAP_CACHE,
                                refresh: bool = False) -> tuple[str, Path]:
    """Install StatFunGen/susieR (the SuSiE 2.0 fork) into a project-local
    R library so we don't touch the user's system library. Returns
    (Rscript path, R library path).
    """
    rscript = ensure_r_installed()
    r_lib = cache_root / "r_library"
    r_lib.mkdir(parents=True, exist_ok=True)
    if refresh:
        for child in r_lib.iterdir():
            shutil.rmtree(child, ignore_errors=True)
    # Smoke-check the install.
    probe = subprocess.run(
        [rscript, "-e",
         f'.libPaths(c("{r_lib}", .libPaths())); '
         f'suppressMessages(library(susieR)); '
         f'cat(packageVersion("susieR"),"\\n")'],
        capture_output=True, text=True)
    if probe.returncode == 0 and probe.stdout.strip():
        return rscript, r_lib

    print(f"[susie_inf:susie2] installing StatFunGen/susieR (SuSiE 2.0) "
          f"into {r_lib} ...", file=sys.stderr)
    install_script = (
        f'.libPaths(c("{r_lib}", .libPaths()));\n'
        f'options(repos = c(CRAN = "https://cloud.r-project.org"));\n'
        f'if (!requireNamespace("remotes", quietly=TRUE)) '
        f'install.packages("remotes", lib="{r_lib}");\n'
        f'remotes::install_github("StatFunGen/susieR", lib="{r_lib}", '
        f'upgrade="never", quiet=TRUE);\n'
        f'cat("susieR=",as.character(packageVersion("susieR")),"\\n", sep="");\n'
    )
    r = subprocess.run([rscript, "-e", install_script],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: R install of StatFunGen/susieR failed:\n"
                 f"stdout={r.stdout[-1000:]}\nstderr={r.stderr[-1000:]}")
    return rscript, r_lib


def run_susie_inf_susie2_sumstats(args, rscript: str, r_lib: Path) -> int:
    """Run the SuSiE-inf model via StatFunGen/susieR (SuSiE 2.0) on K=1
    summary statistics. Implements `susie_rss(z, R, n,
    unmappable_effects = "inf")` — the SuSiE 2.0 codepath for the same
    spike+slab+infinitesimal model as the Cui 2024 Python impl, but
    with cleaner convergence on real GWAS data. The OKG provenance cites
    paper:susie_ash_2025 / software:susie_ash (the canonical "SuSiE 2.0"
    package) and method:susie_inf (the model).
    """
    if len(args.z) != 1:
        sys.exit("ERROR: --method susie_inf --susie-inf-impl susie2 "
                 "requires K=1 (single-ancestry); got {} --z files"
                 .format(len(args.z)))
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    import tempfile
    r_script = _susie_inf_susie2_r_wrapper()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False,
                                       dir=str(out_prefix.parent)) as rf:
        rf.write(r_script)
        r_path = rf.name
    cmd = [
        rscript, r_path,
        "--z", str(args.z[0]),
        "--ld", str(args.ld[0]),
        "--n", str(args.n[0]),
        "--out", str(out_prefix),
        "--L", str(args.L),
        "--coverage", str(args.threshold),
        "--purity", str(args.purity),
        "--lib", str(r_lib),
    ]
    print(f"[finemap susie_inf:susie2] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)
    try:
        os.unlink(r_path)
    except OSError:
        pass
    if rc != 0:
        print(f"[finemap susie_inf (susie2)] Rscript exit code {rc}",
              file=sys.stderr)
        return rc

    # Read susieR package version from a sentinel the R script wrote.
    susie_ver = None
    ver_path = Path(f"{out_prefix}.susier_version.txt")
    if ver_path.exists():
        susie_ver = ver_path.read_text().strip()
    found = resolve_okg(args.okg_repo,
                         list(OKG_NODES_SUSIE_INF_SUSIE2.values()))
    okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE_INF_SUSIE2.items()
                     if v in found}
    summary = _parse_outputs(out_prefix)
    manifest = {
        "subcommand": "sumstats",
        "method": "susie_inf",
        "susie_inf_impl": "susie2",
        "output_prefix": str(out_prefix),
        "susie2_repo": SUSIE_ASH_REPO,
        "susier_version": susie_ver,
        "r_library": str(r_lib),
        "n_ancestries": 1,
        "inputs": {
            "z":  [str(args.z[0])],
            "ld": [str(args.ld[0])],
            "n":  [int(args.n[0])],
        },
        "infer_params": {
            "L": args.L, "max_iter": args.max_iter,
            "coverage": args.threshold, "purity": args.purity,
            "unmappable_effects": "inf",
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{out_prefix}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    print(f"summary: {summary.get('n_cs', '?')} CSs, max PIP "
          f"{summary.get('max_pip', '?')}")
    return 0


def _susie_inf_susie2_r_wrapper() -> str:
    """Return an R script that loads z + LD, runs SuSiE 2.0 (susie_rss
    with unmappable_effects='inf'), and writes the skill's standard
    `<prefix>.cs.tsv` + `<prefix>.weights.tsv` + `<prefix>.log` outputs."""
    return r'''
suppressMessages({
  args <- commandArgs(trailingOnly=TRUE)
  parse_arg <- function(name, default=NULL) {
    i <- match(paste0("--", name), args)
    if (is.na(i)) return(default)
    args[i+1]
  }
  z_path <- parse_arg("z"); ld_path <- parse_arg("ld")
  n <- as.integer(parse_arg("n")); out_pref <- parse_arg("out")
  L <- as.integer(parse_arg("L", "10"))
  coverage <- as.numeric(parse_arg("coverage", "0.95"))
  purity <- as.numeric(parse_arg("purity", "0.5"))
  rlib <- parse_arg("lib", NULL)
  if (!is.null(rlib)) .libPaths(c(rlib, .libPaths()))
  library(susieR)

  # Read z: TSV/CSV with header. Detect Z col or BETA+SE.
  read_z <- function(p) {
    df <- read.table(p, header=TRUE, sep="\t", check.names=FALSE,
                      stringsAsFactors=FALSE)
    if (ncol(df) <= 1) df <- read.table(p, header=TRUE, sep=",",
                                          check.names=FALSE,
                                          stringsAsFactors=FALSE)
    if (ncol(df) <= 1) df <- read.table(p, header=TRUE,
                                          check.names=FALSE,
                                          stringsAsFactors=FALSE)
    cols <- tolower(colnames(df))
    snp_c <- intersect(c("snp","snpid","rsid","snp_id","variant_id","id"), cols)[1]
    if (is.na(snp_c)) snp <- as.character(seq_len(nrow(df))) else snp <- as.character(df[[colnames(df)[match(snp_c, cols)]]])
    z_c <- intersect(c("z","zscore","z_score"), cols)[1]
    if (!is.na(z_c)) {
      z <- as.numeric(df[[colnames(df)[match(z_c, cols)]]])
    } else {
      bc <- intersect(c("beta","effect_size","b","effect"), cols)[1]
      sc <- intersect(c("se","standard_error","stderr","sderr"), cols)[1]
      if (is.na(bc) || is.na(sc)) stop(sprintf("No Z or BETA+SE in %s (cols=%s)", p, paste(cols, collapse=",")))
      z <- as.numeric(df[[colnames(df)[match(bc, cols)]]]) /
            as.numeric(df[[colnames(df)[match(sc, cols)]]])
    }
    list(snp=snp, z=z)
  }
  zd <- read_z(z_path)

  read_ld <- function(p) {
    if (endsWith(tolower(p), ".npy")) stop(".npy LD not supported in R; convert to TSV or .rds")
    if (endsWith(tolower(p), ".rds")) return(readRDS(p))
    as.matrix(read.table(p, header=FALSE))
  }
  R <- read_ld(ld_path)
  if (nrow(R) != length(zd$z)) stop(sprintf("LD %dx%d != z len %d", nrow(R), ncol(R), length(zd$z)))

  fit <- susie_rss(z = zd$z, R = R, n = n, L = L,
                    unmappable_effects = "inf",
                    coverage = coverage, min_abs_corr = purity)

  pip <- susie_get_pip(fit)
  cs_list <- susie_get_cs(fit, coverage = coverage, min_abs_corr = purity)$cs
  # weights.tsv
  weights_df <- data.frame(SNPIndex = seq_along(zd$z) - 1L,
                            SNP = zd$snp,
                            PIP = pip)
  write.table(weights_df, file = paste0(out_pref, ".weights.tsv"),
               sep = "\t", quote = FALSE, row.names = FALSE)
  # cs.tsv
  cs_rows <- list()
  if (!is.null(cs_list)) {
    for (ci in seq_along(cs_list)) {
      for (idx in cs_list[[ci]]) {
        cs_rows[[length(cs_rows)+1]] <- data.frame(
          CSIndex = ci, SNPIndex = idx - 1L,
          SNP = zd$snp[idx], pip_all = pip[idx])
      }
    }
  }
  cs_df <- if (length(cs_rows)) do.call(rbind, cs_rows) else
            data.frame(CSIndex=integer(), SNPIndex=integer(), SNP=character(), pip_all=numeric())
  write.table(cs_df, file = paste0(out_pref, ".cs.tsv"),
               sep = "\t", quote = FALSE, row.names = FALSE)

  # version sentinel
  writeLines(as.character(packageVersion("susieR")), paste0(out_pref, ".susier_version.txt"))
  # log
  log_lines <- c(
    "susie_inf done (susie_rss + unmappable_effects=inf, SuSiE 2.0)",
    sprintf("N=%d L=%d m=%d", n, L, length(zd$z)),
    sprintf("n_cs=%d", if (is.null(cs_list)) 0 else length(cs_list)),
    sprintf("max_pip=%.6g", max(pip)),
    sprintf("susieR_version=%s", packageVersion("susieR"))
  )
  writeLines(log_lines, paste0(out_pref, ".log"))
  message(sprintf("susie_inf:susie2: m=%d N=%d CSs=%d max_PIP=%.4f",
                   length(zd$z), n, if (is.null(cs_list)) 0 else length(cs_list), max(pip)))
})
'''


# ---------------------------- OKG-driven method recommender ----------------------------

def _okg_get_node(okg_repo: Path, node_id: str) -> Optional[dict]:
    """Single-node lookup via MCP. Returns the node's attrs dict or None."""
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
        return None
    env = os.environ.copy()
    env.setdefault("OKG_DSN",
                   "postgres://postgres:okg@localhost:5449/statgen_analysis")
    proc = subprocess.Popen(
        ["uv", "run", "--extra", "mcp", "python",
         "deployments/statgen-analysis/server.py"],
        cwd=str(okg_repo), env=env,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, bufsize=1,
    )
    def send(m): proc.stdin.write(json.dumps(m) + "\n"); proc.stdin.flush()
    def read(): return json.loads(proc.stdout.readline())
    out: Optional[dict] = None
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "finemap-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": "get_node",
                         "arguments": {"node_id": node_id}}})
        resp = read()
        sc = resp.get("result", {}).get("structuredContent") or {}
        if sc.get("node_id") == node_id:
            out = sc
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return out


def _okg_list_neighbors(okg_repo: Path, node_id: str,
                         edge_type: Optional[str] = None,
                         direction: str = "out") -> list[dict]:
    """List neighbor nodes via MCP. Returns the list of neighbor dicts."""
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
        return []
    env = os.environ.copy()
    env.setdefault("OKG_DSN",
                   "postgres://postgres:okg@localhost:5449/statgen_analysis")
    proc = subprocess.Popen(
        ["uv", "run", "--extra", "mcp", "python",
         "deployments/statgen-analysis/server.py"],
        cwd=str(okg_repo), env=env,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, bufsize=1,
    )
    def send(m): proc.stdin.write(json.dumps(m) + "\n"); proc.stdin.flush()
    def read(): return json.loads(proc.stdout.readline())
    out: list[dict] = []
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "finemap-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        args = {"node_id": node_id, "direction": direction, "limit": 50}
        if edge_type is not None:
            args["edge_type"] = edge_type
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": "list_neighbors", "arguments": args}})
        resp = read()
        sc = resp.get("result", {}).get("structuredContent") or {}
        out = sc.get("neighbors", []) or sc.get("results", []) or []
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return out


def recommend_method(okg_repo: Optional[Path],
                      dataset_id: Optional[str],
                      n_inputs: int) -> dict:
    """Choose a sumstats fine-mapping method given OKG context.

    Rules (in priority order):
      1. K=2+ inputs supplied directly --> 'sushie' (multi-ancestry).
      2. dataset has multi-ancestry sibling datasets in the OKG --> 'sushie'
         (recommend; user must pass per-ancestry sumstats themselves).
      3. dataset N (n_cases+n_controls or n_initial) >= 100_000 -->
         'susie_inf' (SuSiE-inf model; default impl is SuSiE 2.0 /
         StatFunGen susieR — see --susie-inf-impl). The infinitesimal
         polygenic background helps at large N where standard SuSiE-RSS
         can be biased by uncaptured polygenicity.
      4. else --> 'sushie' (K=1 SuSiE-RSS-equivalent; sushie's default).

    Returns: {"method": str, "rationale": str,
              "signals": {...graph-derived facts...}}
    """
    signals: dict = {"n_inputs": n_inputs}
    if n_inputs >= 2:
        return {"method": "sushie",
                "rationale": ("multi-ancestry inputs supplied (K=%d) -> use "
                              "sushie" % n_inputs),
                "signals": signals}

    if okg_repo is None or dataset_id is None:
        return {"method": "sushie",
                "rationale": ("no OKG dataset id supplied; defaulting to "
                              "sushie K=1 (SuSiE-RSS-equivalent)"),
                "signals": signals}

    node = _okg_get_node(okg_repo, dataset_id)
    if node is None:
        return {"method": "sushie",
                "rationale": (f"dataset {dataset_id} not found in OKG; "
                              f"defaulting to sushie K=1"),
                "signals": signals}
    attrs = node.get("attrs") or {}
    # Resolve total N.
    n_cases = attrs.get("n_cases")
    n_controls = attrs.get("n_controls")
    n_initial = attrs.get("n_initial")
    if n_cases is not None and n_controls is not None:
        n_total = int(n_cases) + int(n_controls)
        signals["n_cases"] = int(n_cases)
        signals["n_controls"] = int(n_controls)
    elif n_initial is not None:
        n_total = int(n_initial)
        signals["n_initial"] = int(n_initial)
    else:
        n_total = None
    signals["n_total"] = n_total
    signals["ancestry_scope"] = attrs.get("ancestry_scope")

    # Check for multi-ancestry sibling datasets for the same trait.
    trait_neighbors = _okg_list_neighbors(
        okg_repo, dataset_id, edge_type="references", direction="out")
    trait_ids = [n.get("node_id") for n in trait_neighbors
                  if (n.get("node_id") or "").startswith("trait:")]
    sibling_ancestries = []
    for tid in trait_ids:
        # Reverse: which datasets reference this trait?
        ds_neighbors = _okg_list_neighbors(
            okg_repo, tid, edge_type="references", direction="in")
        for ds in ds_neighbors:
            nid = ds.get("node_id") or ""
            if not nid.startswith("dataset:") or nid == dataset_id:
                continue
            ds_node = _okg_get_node(okg_repo, nid)
            if ds_node is None:
                continue
            ds_anc = (ds_node.get("attrs") or {}).get("ancestry_scope")
            if ds_anc:
                sibling_ancestries.append({"dataset_id": nid,
                                             "ancestry_scope": ds_anc})
    signals["sibling_datasets"] = sibling_ancestries[:5]

    own_anc = (attrs.get("ancestry_scope") or "").lower()
    sib_ancs = [s["ancestry_scope"].lower() for s in sibling_ancestries]
    has_other_ancestry = any(
        own_anc not in a and a not in own_anc and "european" not in (own_anc + a).lower()
        for a in sib_ancs
    ) or any(
        "east_asian" in a or "african" in a or "south_asian" in a
        for a in sib_ancs if a not in own_anc
    )

    # Decision tree.
    if has_other_ancestry:
        return {"method": "sushie",
                "rationale": (f"OKG knows {len(sibling_ancestries)} "
                              f"sibling dataset(s) for the same trait in "
                              f"other ancestries; recommend multi-ancestry "
                              f"SuShiE if you can supply their sumstats"),
                "signals": signals}
    if n_total is not None and n_total >= 500_000:
        return {"method": "susie_inf",
                "rationale": (f"very large N ({n_total:,}) - susie_inf "
                              f"(SuSiE 2.0) models a polygenic "
                              f"infinitesimal tail, recommended for high-N "
                              f"GWAS where standard SuSiE-RSS can be biased "
                              f"by uncaptured polygenicity"),
                "signals": signals}
    if n_total is not None and n_total >= 100_000:
        return {"method": "susie_inf",
                "rationale": (f"moderate-large N ({n_total:,}) - susie_inf "
                              f"(default impl SuSiE 2.0 / StatFunGen "
                              f"susieR) offers a richer prior than vanilla "
                              f"SuSiE-RSS"),
                "signals": signals}
    n_str = f"{n_total:,}" if n_total else "unknown"
    return {"method": "sushie",
            "rationale": (f"small-moderate N ({n_str}); sushie K=1 "
                          f"(SuSiE-RSS-equivalent) is the conservative "
                          f"default"),
            "signals": signals}


# ---------------------------- CLI ----------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--vcf", type=Path, nargs="+", required=True,
                   help="One or more VCF files (one per ancestry)")
    p.add_argument("--pheno", type=Path, nargs="+", required=True,
                   help="One or more phenotype files (sushie format)")
    p.add_argument("--covar", type=Path, nargs="+", default=None)
    p.add_argument("--out", type=Path, required=True,
                   help="Output prefix")
    p.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    p.add_argument("--sushie-commit", dest="commit", type=str, default=None,
                   help="Pin sushie to a specific commit")
    p.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    p.add_argument("--refresh", action="store_true",
                   help="Re-clone sushie and re-run verify-install")
    p.add_argument("--re-verify", action="store_true",
                   help="Force re-running verify-install even if sentinel exists")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None,
                   help="OKG repo path (honors $OKG_REPO)")
    p.add_argument("extra", nargs=argparse.REMAINDER,
                   help="Extra flags forwarded to `sushie finemap`")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--verify-install", action="store_true",
                   help="Run the bundled 3-ancestry tutorial and exit")
    sub = p.add_subparsers(dest="subcmd")

    psusie = sub.add_parser("susie", help="Single-ancestry fine-mapping")
    _add_common(psusie)
    psushie = sub.add_parser("sushie", help="Multi-ancestry fine-mapping")
    _add_common(psushie)

    preg = sub.add_parser("region",
                           help="Single-locus fine-mapping from GWAS sumstats "
                                "+ a reference VCF/PLINK (sushie computes "
                                "LD internally; no precomputed matrix needed)")
    preg.add_argument("--gwas-sumstats", dest="gwas_sumstats", type=Path,
                      required=True,
                      help="GWAS sumstats file (GWAS Catalog harmonised "
                           "TSV, LDSC munged .sumstats.gz, or any TSV "
                           "with chrom/snp/pos/effect_allele/other_allele "
                           "+ z or beta/se)")
    preg.add_argument("--chrom", type=int, required=True, choices=range(1, 23),
                      help="Chromosome (1-22)")
    preg.add_argument("--start", type=int, required=True,
                      help="Window start, bp")
    preg.add_argument("--end", type=int, required=True,
                      help="Window end, bp")
    preg.add_argument("--N", type=int, required=True,
                      help="GWAS total sample size for this trait")
    preg.add_argument("--ref-vcf", dest="ref_vcf", type=Path, nargs="+",
                      help="Reference genotype VCF (sushie computes LD "
                           "from this; usually 1000G or UKB on a build "
                           "matching the GWAS)")
    preg.add_argument("--ref-plink", dest="ref_plink", type=Path, nargs="+",
                      help="Reference genotype PLINK1.9 prefix(es) "
                           "(alternative to --ref-vcf)")
    preg.add_argument("--ref-bgen", dest="ref_bgen", type=Path, nargs="+",
                      help="Reference genotype BGEN 1.3 file(s) "
                           "(alternative to --ref-vcf)")
    preg.add_argument("--ref-1000g", dest="ref_1000g", type=str,
                      choices=list(_POP_TO_SUPERPOP),
                      help="Convenience: auto-fetch the 1000G phase3 "
                           "VCF slice for the given super-population "
                           "(eur/eas/afr/sas/amr) at --chrom:--start-end, "
                           "subset to that population, and use it as the "
                           "LD reference. Installs htslib + bcftools on "
                           "first use; caches the slice under "
                           "~/.cache/finemap/1000g_<build>_<pop>/.")
    preg.add_argument("--ref-1000g-build", dest="ref_1000g_build",
                      type=str, choices=["hg19", "hg38"], default="hg19",
                      help="Build for --ref-1000g (default: hg19)")
    preg.add_argument("--ld", type=Path, nargs="+",
                      help="Pre-computed LD matrix (tsv/tsv.gz) — bypasses "
                           "the reference-genotype LD computation")
    preg.add_argument("--out", type=Path, required=True,
                      help="Output prefix")
    preg.add_argument("--okg-dataset-id", dest="okg_dataset_id", type=str,
                      help="OKG dataset_metadata node ID for the GWAS "
                           "(recorded in the manifest's okg_node_ids)")
    preg.add_argument("--okg-ld-panel-id", dest="okg_ld_panel_id", type=str,
                      help="OKG ld_panel node ID for the reference "
                           "genotypes (recorded in the manifest)")
    preg.add_argument("--method", type=str, default="sushie",
                      choices=list(SUMSTATS_METHODS),
                      help="Fine-mapping backend for the region. `sushie` "
                           "(default; SuShiE region path with internal LD), "
                           "`susie_inf` (SuSiE-inf model, K=1; requires "
                           "--ref-plink; impl picked via --susie-inf-impl), "
                           "or `auto` (OKG-recommended given "
                           "--okg-dataset-id).")
    preg.add_argument("--susie-inf-impl", dest="susie_inf_impl", type=str,
                      default="susie2", choices=list(SUSIE_INF_IMPLS),
                      help="Backend impl for --method susie_inf. `susie2` "
                           "(default) = StatFunGen/susieR (SuSiE 2.0; Wang "
                           "2025); `finucane` = FinucaneLab/fine-mapping-inf "
                           "(Cui 2024 Python). Same statistical model; "
                           "different convergence behavior. `susie2` is "
                           "preferred on real GWAS (cleaner convergence + "
                           "newer paper citation).")
    preg.add_argument("--snp-id-format", dest="snp_id_format", type=str,
                      default="auto", choices=["auto", "rsid", "chrpos"],
                      help="SNP-ID convention emitted in the GWAS file "
                           "passed to sushie. Must match the reference "
                           "panel's variant IDs (sushie joins by exact "
                           "string equality). 'auto' (default) peeks at "
                           "the panel's .bim / VCF and picks 'rsid' if "
                           ">50%% of IDs are rs<digits>, else 'chrpos'. "
                           "Override to force one format.")
    preg.add_argument("--annotate", type=str, default="both",
                      choices=["none", "gnomad", "gtex", "both"],
                      help="Post-fine-mapping variant annotation, delegated "
                           "to the sibling `variant-annotate` skill. "
                           "`gnomad` adds most-severe consequence + gene + "
                           "NFE allele frequency. `gtex` adds top eQTL + "
                           "sQTL across tissues. `both` (default) adds both "
                           "column groups to the CS summary table. Network "
                           "calls; responses cached at "
                           "~/.cache/variant-annotate/.")
    preg.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    preg.add_argument("--sushie-commit", dest="commit", type=str, default=None)
    preg.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    preg.add_argument("--refresh", action="store_true")
    preg.add_argument("--re-verify", action="store_true")
    preg.add_argument("--okg-repo", type=Path,
                      default=Path(os.environ["OKG_REPO"])
                              if os.environ.get("OKG_REPO") else None)
    preg.add_argument("extra", nargs=argparse.REMAINDER,
                      help="Extra flags forwarded to `sushie finemap`")

    pss = sub.add_parser("sumstats",
                          help="Sumstats fine-mapping via infer_sushie_ss "
                               "(K=1 SuSiE-RSS or K>=2 SuShiE)")
    pss.add_argument("--z", type=Path, nargs="+", required=True,
                     help="One or more z-score files (TSV with Z, or "
                          "BETA+SE columns) — one per ancestry")
    pss.add_argument("--ld", type=Path, nargs="+", required=True,
                     help="One or more LD matrices (.npy or whitespace TSV) "
                          "— one per ancestry, same SNP order as --z")
    pss.add_argument("--n", type=int, nargs="+", required=True,
                     help="Sample size per ancestry (one int per --z file)")
    pss.add_argument("--out", type=Path, required=True,
                     help="Output prefix")
    pss.add_argument("--L", type=int, default=10,
                     help="Max number of single effects (default 10)")
    pss.add_argument("--max-iter", type=int, default=500)
    pss.add_argument("--min-tol", type=float, default=1e-4)
    pss.add_argument("--threshold", type=float, default=0.95,
                     help="Credible set coverage threshold (default 0.95)")
    pss.add_argument("--purity", type=float, default=0.5,
                     help="Min CS purity (default 0.5)")
    pss.add_argument("--min-snps", type=int, default=100,
                     help="Min SNPs to fine-map (default 100)")
    pss.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    pss.add_argument("--sushie-commit", dest="commit", type=str, default=None)
    pss.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    pss.add_argument("--refresh", action="store_true")
    pss.add_argument("--re-verify", action="store_true")
    pss.add_argument("--okg-repo", type=Path,
                     default=Path(os.environ["OKG_REPO"])
                             if os.environ.get("OKG_REPO") else None)
    pss.add_argument("--method", type=str, default="sushie",
                     choices=list(SUMSTATS_METHODS),
                     help="Fine-mapping method: `sushie` (default; "
                          "K=1 SuSiE-RSS-equivalent or K>=2 SuShiE), "
                          "`susie_inf` (SuSiE-inf model, K=1; impl picked "
                          "via --susie-inf-impl), or `auto` (recommend "
                          "via OKG given --okg-dataset-id).")
    pss.add_argument("--susie-inf-impl", dest="susie_inf_impl", type=str,
                     default="susie2", choices=list(SUSIE_INF_IMPLS),
                     help="Backend impl for --method susie_inf. `susie2` "
                          "(default) = StatFunGen/susieR (SuSiE 2.0; Wang "
                          "2025); `finucane` = FinucaneLab/fine-mapping-inf "
                          "(Cui 2024 Python). Same model; SuSiE 2.0 "
                          "converges more cleanly on real GWAS.")
    pss.add_argument("--okg-dataset-id", dest="okg_dataset_id", type=str,
                     help="OKG dataset_metadata node ID for the input "
                          "(enables --method auto, recorded in manifest)")

    # `recommend` subcommand: prints a method choice + rationale, no run.
    prec = sub.add_parser("recommend",
                           help="Recommend a sumstats fine-mapping method "
                                "from OKG context (no run)")
    prec.add_argument("--okg-dataset-id", dest="okg_dataset_id", type=str,
                      required=True,
                      help="OKG dataset_metadata node ID to inspect")
    prec.add_argument("--n-inputs", type=int, default=1,
                      help="Number of ancestry-specific inputs you plan "
                           "to supply (default 1)")
    prec.add_argument("--okg-repo", type=Path,
                      default=Path(os.environ["OKG_REPO"])
                              if os.environ.get("OKG_REPO") else None)

    args, _ = p.parse_known_args()

    # Pure verify-install mode.
    if args.verify_install:
        repo_dir = ensure_repo(DEFAULT_REPO, None,
                                CACHE_ROOT, refresh=False)
        venv_py = ensure_sushie_installed(repo_dir)
        return verify_install(repo_dir, venv_py=venv_py)

    if args.subcmd not in ("susie", "sushie", "sumstats", "region",
                            "recommend"):
        p.print_help(sys.stderr)
        return 2

    # Re-parse for the subcommand's flags.
    args = p.parse_args()

    # `recommend` is read-only and doesn't need a sushie install.
    if args.subcmd == "recommend":
        rec = recommend_method(args.okg_repo, args.okg_dataset_id,
                                args.n_inputs)
        print(json.dumps(rec, indent=2))
        return 0

    # Common setup.
    repo_dir = ensure_repo(args.repo_url, args.commit,
                            Path(args.repo_cache), refresh=args.refresh)
    venv_py = ensure_sushie_installed(repo_dir)

    # First-run verify-install gate.
    if args.refresh or args.re_verify or not is_verified():
        rc = verify_install(repo_dir, venv_py=venv_py)
        if rc != 0:
            sys.exit("verify-install failed; refusing to proceed with user run")

    if args.subcmd == "sumstats":
        # Resolve --method auto via OKG context.
        method = args.method
        if method == "auto":
            rec = recommend_method(args.okg_repo,
                                    getattr(args, "okg_dataset_id", None),
                                    n_inputs=len(args.z))
            method = rec["method"]
            print(f"[finemap auto] recommended method: {method} "
                  f"({rec['rationale']})", file=sys.stderr)
            args.method = method  # so the runner records the resolved choice
            args.recommendation = rec
        if method == "susie_inf":
            impl = getattr(args, "susie_inf_impl", "susie2")
            if impl == "susie2":
                rscript, r_lib = ensure_susie_inf_susie2_installed()
                return run_susie_inf_susie2_sumstats(args, rscript=rscript,
                                                       r_lib=r_lib)
            elif impl == "finucane":
                si_repo = ensure_susie_inf_installed(venv_py)
                return run_susie_inf_sumstats(args, venv_py=venv_py,
                                               susie_inf_repo=si_repo)
            else:
                sys.exit(f"REFUSED: unknown --susie-inf-impl {impl!r}; "
                         f"expected one of {SUSIE_INF_IMPLS}")
        # default: sushie
        return run_sumstats(args, repo_dir)
    if args.subcmd == "region":
        return run_region(args, repo_dir, venv_py=venv_py)

    # Validate ancestry count vs subcommand.
    if args.subcmd == "susie" and len(args.vcf) != 1:
        p.error("--vcf must be a single file in `susie` mode "
                "(got {} files)".format(len(args.vcf)))
    if args.subcmd == "sushie" and len(args.vcf) < 2:
        p.error("`sushie` mode requires K>=2 ancestries; "
                "got {} VCF".format(len(args.vcf)))
    if len(args.vcf) != len(args.pheno):
        p.error("number of --vcf files must match number of --pheno files "
                "({} vs {})".format(len(args.vcf), len(args.pheno)))
    if args.covar and len(args.covar) != len(args.vcf):
        p.error("number of --covar files must match number of --vcf files")

    return run_finemap(args, args.subcmd, repo_dir, venv_py=venv_py)


if __name__ == "__main__":
    sys.exit(main())
