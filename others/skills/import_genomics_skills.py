"""Import the reviewed, commit-pinned genomics skills; --check is offline.

This additive importer does not edit the main catalog, install scientific
dependencies, execute imported workflows, or overwrite existing skill folders.
It needs PyYAML in the documentation/build environment.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import posixpath
import re
import sys
import tarfile
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
DATE = "2026-09-20"
MANIFEST = Path("registry/genomics-new-skills.json")
SOURCES = [
    {
        "id": "gigwa-mcp",
        "repository": "gkanogiannis/Gigwa-MCP",
        "revision": "d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf",
        "archive_sha256": "c74191ff79b8adf0f085158514c25e268f8eb8ce01d0e557236b914f5c30b1c8",
        "content_license": "Apache-2.0",
        "attribution": "Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors",
        "prefix": "skills",
    },
    {
        "id": "haploqtl",
        "repository": "tayloranthonyanderson/haploqtl",
        "revision": "210eace2912c7bf11d119cb503c89f8f6cad479f",
        "archive_sha256": "e4259f0c4ae53e2a7e607e8cf9cd70dcfdbbf23a39d6c74b3a17972b472f3a66",
        "content_license": "MIT",
        "attribution": "Copyright 2026 Taylor A. Anderson; haploqtl contributors",
        "prefix": "skills",
    },
    {
        "id": "tangermeme",
        "repository": "jmschrei/tangermeme",
        "revision": "d5c51b83d374576b64effb9284fc79a71a604734",
        "archive_sha256": "110e8bd48e5827a27b0a501e39b1e8e0c52c5d0aa8405a2d62e9fe6aa2302715",
        "content_license": "MIT",
        "attribution": "Copyright 2024-2026 Jacob Schreiber; tangermeme contributors",
        "prefix": "tangermeme/_skills/data",
    },
    {
        "id": "statgen-skills",
        "repository": "svdorn/statgen-skills",
        "revision": "dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0",
        "archive_sha256": "755eb7296e5aa7c1d22fd45fb09c46758309da69b2152782c6ff55a3b9ff8bd4",
        "content_license": "MIT",
        "attribution": "Copyright 2026 Stephen Dorn; statgen-skills contributors",
        "prefix": "",
    },
    {
        "id": "biomcp",
        "repository": "genomoncology/biomcp",
        "revision": "3b1d6452ab63741562d2ef97495a97281b304729",
        "archive_sha256": "b8c849562219b87dfdba7774e041ce574e74512d4b3bedccf160296209befa60",
        "content_license": "MIT",
        "attribution": "GenomOncology and BioMCP contributors; see retained LICENSE/CITATION.cff",
        "prefix": "skills",
    },
]

# category, concise Chinese purpose, required Python modules, commands, environment.
RUNTIMES: dict[str, tuple[str, str, list[str], list[str], list[str]]] = {
    "gigwa-diversity-report": (
        "population",
        "Gigwa 群体多样性与关系报告",
        ["gigwa_mcp"],
        ["gigwa-mcp"],
        [],
    ),
    "gigwa-explore-instance": (
        "databases",
        "Gigwa 种质基因型数据发现",
        ["gigwa_mcp"],
        ["gigwa-mcp"],
        [],
    ),
    "gigwa-import-and-qc": ("variants", "Gigwa 基因型导入与质控", ["gigwa_mcp"], ["gigwa-mcp"], []),
    "gigwa-qc-triage": ("variants", "Gigwa 基因型质控评估", ["gigwa_mcp"], ["gigwa-mcp"], []),
    "gigwa-region-scan": ("population", "Gigwa 区域变异与多样性", ["gigwa_mcp"], ["gigwa-mcp"], []),
    "qtl-candidate-gene": (
        "annotation",
        "番茄 QTL 候选基因与辅助选择标记",
        ["numpy", "allel"],
        ["python3"],
        [],
    ),
    "tangermeme": ("models", "DNA 模型归因、扰动与序列设计", ["tangermeme", "torch"], [], []),
    "finemap": ("gwas", "多后端 GWAS 精细定位与变异注释", ["numpy", "pandas"], ["python3"], []),
    "gwas-fetch": ("databases", "GWAS Catalog 汇总统计下载与溯源", [], ["python3"], []),
    "ldsc": (
        "gwas",
        "LD Score 遗传力与遗传相关估计",
        ["numpy", "pandas", "scipy"],
        ["python3", "git"],
        [],
    ),
    "liftover": ("intervals", "GWAS 汇总统计坐标转换", ["pyliftover"], ["python3"], []),
    "locuszoom": ("gwas", "GWAS 区域关联与精细定位绘图", [], ["python3", "Rscript"], []),
    "manhattan-qq": (
        "gwas",
        "GWAS Manhattan、QQ 与质控报告",
        ["numpy", "pandas", "scipy", "matplotlib"],
        ["python3"],
        [],
    ),
    "okg-search": (
        "databases",
        "统计遗传知识库 GWAS 数据集检索",
        ["psycopg"],
        ["python3"],
        ["OKG_DSN"],
    ),
    "prs": (
        "gwas",
        "SBayesRC 多基因评分权重估计",
        [],
        ["python3", "Rscript", "uv"],
        ["OKG_REPO", "OKG_DSN"],
    ),
    "twas": (
        "gwas",
        "FUSION 转录组关联与 GTEx 权重",
        [],
        ["python3", "Rscript", "git"],
        ["OKG_REPO", "OKG_DSN"],
    ),
    "variant-annotate": ("variants", "gnomAD 与 GTEx 变异证据注释", [], ["python3"], []),
    "biomcp": ("clinical", "BioMCP 基因、变异、药物基因组与文献查询", [], ["biomcp"], []),
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def progress(done: int, total: int, text: str) -> None:
    width = 20
    filled = width * done // total
    print(f"[{'#' * filled}{'-' * (width - filled)}] {done}/{total} {text}", file=sys.stderr)


def source_files(source: dict[str, str], cache: Path) -> dict[str, tuple[bytes, int]]:
    archive = cache / (source["repository"].replace("/", "__") + ".tar.gz")
    if not archive.exists():
        url = f"https://codeload.github.com/{source['repository']}/tar.gz/{source['revision']}"
        with urllib.request.urlopen(url, timeout=90) as response:
            blob = response.read(100_000_001)
        if len(blob) > 100_000_000 or digest(blob) != source["archive_sha256"]:
            raise ValueError(f"Unexpected archive for {source['repository']}")
        cache.mkdir(parents=True, exist_ok=True)
        with archive.open("xb") as handle:
            handle.write(blob)
    if digest(archive.read_bytes()) != source["archive_sha256"]:
        raise ValueError(f"Archive SHA-256 mismatch: {archive}")
    result = {}
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            parts = PurePosixPath(member.name).parts
            if (
                not parts
                or not parts[0].endswith(source["revision"])
                or ".." in parts
                or member.name.startswith("/")
            ):
                raise ValueError(f"Unsafe or unpinned archive path: {member.name}")
            relative = str(PurePosixPath(*parts[1:]))
            handle = tar.extractfile(member)
            assert handle is not None
            if relative in result:
                raise ValueError(f"Duplicate archive path: {relative}")
            result[relative] = (handle.read(), 0o755 if member.mode & 0o111 else 0o644)
    if "LICENSE" not in result:
        raise ValueError(f"Missing reviewed license: {source['repository']}")
    return result


def metadata(blob: bytes) -> tuple[dict[str, Any], int]:
    match = re.match(rb"\A---\r?\n(.*?)\r?\n---(?:\r?\n|$)", blob, re.S)
    if match is None:
        raise ValueError("SKILL.md has no YAML frontmatter")
    data = yaml.safe_load(match[1])
    if not isinstance(data, dict) or not data.get("name") or not data.get("description"):
        raise ValueError("Skill name and description are required")
    if not re.fullmatch(r"[a-z0-9-]{1,63}", data["name"]):
        raise ValueError(f"Invalid skill name: {data['name']}")
    return data, match.end()


def guide(skill_id: str, source: dict[str, str]) -> str:
    repo, revision = source["repository"], source["revision"]
    text = f"""# DNAKit execution: {skill_id}

This skill is curated from `{repo}` at `{revision}`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor {skill_id}` after catalog integration.
Install the complete folder with `dnakit skills install {skill_id} --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

"""
    if source["id"] == "gigwa-mcp":
        text += """Install `gigwa-mcp` separately, connect its MCP server, and configure a local or
remote Gigwa instance. `GIGWA_URL` selects the instance; credentials are optional for
public reads and required for private data/imports. Inputs must be readable by the server.
Import only into the dataset/server selected by the user. `DNAKit` does not install Gigwa.
Read [NOTICE](NOTICE) for the distinction between Gigwa and this MCP wrapper.

Scientific interpretation: VanRaden G is the additive relationship matrix, whereas
coancestry is conventionally G/2. No universal G threshold proves sibling/clone status.
Complete call rate alone does not establish a bad import. UPGMA from IBS is a sample
distance dendrogram; PCA plus K-means is a structure heuristic, not ADMIXTURE inference.
Apply these qualifications when reading the upstream example wording.
"""
    elif skill_id == "qtl-candidate-gene":
        text += """The full upstream script/reference/data folder is retained. Interval lookup uses
the bundled ITAG4.1 slice offline. Supply your own matching full GFF outside EB-5/EB-9.
MAS analysis additionally requires NumPy/scikit-allel and an input VCF with stated groups;
the disease-group defaults are only the published EB-9 example. UniProt/PubMed steps
use network; optional citation-support checking needs anthropic and ANTHROPIC_API_KEY.
The core interval scripts do not require installing the whole haploqtl project.
Cite Anderson et al. (2024), DOI 10.1111/tpj.16495, and the SGN/UniProt/NCBI data used.
See the retained [CITATION.cff](CITATION.cff).
"""
    elif skill_id == "tangermeme":
        text += f"""Install `tangermeme`/PyTorch in a separate model environment.
The source snapshot
can be installed with `python -m pip install git+https://github.com/{repo}.git@{revision}`.
Use the real `tangermeme.predict`, `tangermeme.deep_lift_shap`,
`tangermeme.variant_effect`, and `tangermeme.design` modules as documented in references.
Supply a compatible trained model and its preprocessing; no pretrained model or weights
are included here. GPU, memelite and model-specific dependencies are optional by task.
"""
    elif skill_id == "biomcp":
        text += """Install the upstream BioMCP CLI separately according to its README, then inspect
`biomcp list` and `biomcp skill list`. The pinned project may use a different runtime
from older Python packages with a similar name; do not substitute an unverified package.
Gene/variant/PGx/literature functions query upstream services and may use downloaded caches
or credentials for particular providers. Use the genomic/biomedical research playbooks
actually relevant to the request; source data, evidence and model limitations remain visible.
SDLC issue-filing instructions are excluded. Cite the retained CITATION.cff and data sources.
"""
    else:
        text += """Run the bundled script from this skill's directory in a dedicated analysis
environment. Upstream helpers may install Python/R packages, clone scientific tools or
download large reference panels on first ANALYSIS use; importing this skill performs none
of these operations. Inspect arguments and reference paths before an analysis run.

OKG integration is a separate service: when using it, supply your own `OKG_REPO`
(containing deployments/statgen-analysis/server.py) and `OKG_DSN`. The repository's
example localhost DSN is not a supplied database. The OKG schema/data/server are not
bundled here. `okg-search` and the current PRS/TWAS metadata paths depend on this
configuration; other workflows have standalone paths documented upstream.
R/PLINK/LDSC/SBayesRC/FUSION/SuSiE and LD/GTEx/reference files are separate dependencies,
not implementations supplied by the skill package. Preserve build, allele and ancestry
matching; report actual solver/reference versions. Runtime downloads can be unpinned
upstream even though this skill is pinned, so record resolved versions for reproducibility.
"""
        if skill_id == "finemap":
            text += """
The optional annotation dependency is bundled under `references/variant-annotate/`.
One resource-loading path is adapted so installing this skill alone retains annotation.
Read that dependency's SKILL.md for its gnomAD/GTEx assumptions and API limitations.
"""
        if skill_id == "locuszoom":
            text += "\nLDlink access additionally requires the user's own API token.\n"
    return text


def make_skill(
    source: dict[str, str], folder: str, upstream: dict[str, tuple[bytes, int]]
) -> tuple[dict[str, Any], dict[str, tuple[bytes, int]]]:
    prefix = folder + "/"
    files = {
        name[len(prefix) :]: value for name, value in upstream.items() if name.startswith(prefix)
    }
    meta, end = metadata(files["SKILL.md"][0])
    skill_id = meta["name"]
    if skill_id not in RUNTIMES:
        raise ValueError(f"Unreviewed skill: {skill_id}")
    upstream_paths = {name: prefix + name for name in files}
    for name in ("LICENSE", "NOTICE", "CITATION.cff"):
        if name in upstream:
            files[name] = upstream[name]
            upstream_paths[name] = name
    files["references/upstream-readme.md"] = upstream["README.md"]
    upstream_paths["references/upstream-readme.md"] = "README.md"
    if skill_id == "finemap":
        for name, value in upstream.items():
            if name.startswith("variant-annotate/"):
                target = "references/" + name
                files[target] = value
                upstream_paths[target] = name
    originals = [
        {"path": name, "upstream_path": original, "sha256": digest(upstream[original][0])}
        for name, original in sorted(upstream_paths.items())
    ]
    modifications = ["DNAKIT.md", "SKILL.md", "UPSTREAM.md"]
    # Keep relative resources local when bundled; point repository-external
    # documentation links at the pinned source so a single skill is portable.
    for name, original_path in upstream_paths.items():
        if not name.endswith(".md"):
            continue
        original, mode = files[name]

        def link(match: re.Match[str], name: str = name, original_path: str = original_path) -> str:
            target = match[1]
            if target.startswith(("https://", "http://", "mailto:", "#")):
                return match[0]
            path, separator, fragment = target.partition("#")
            local = posixpath.normpath(posixpath.join(posixpath.dirname(name), path))
            if local in files or any(key.startswith(local.rstrip("/") + "/") for key in files):
                return match[0]
            remote = posixpath.normpath(posixpath.join(posixpath.dirname(original_path), path))
            url = f"https://github.com/{source['repository']}/blob/{source['revision']}/{remote}"
            return "](" + url + (separator + fragment if separator else "") + ")"

        revised = re.sub(r"\]\(([^\s)]+)\)", link, original.decode()).encode()
        if revised != original:
            files[name] = revised, mode
            modifications.append(name)
    blob, mode = files["SKILL.md"]
    local_meta = dict(meta)
    local_meta["metadata"] = dict(meta.get("metadata", {}))
    if "compatibility" in local_meta:
        local_meta["metadata"]["compatibility"] = local_meta.pop("compatibility")
    descriptions = {
        "finemap": (
            "Fine-map GWAS loci with SuShiE/SuSiE or SuSiE-inf using individual genotypes "
            "or summary statistics and matched LD panels. Supports single/multiple ancestry, "
            "credible sets, optional gnomAD/GTEx annotation and provenance manifests. "
            "Use for posterior inclusion probabilities, credible sets or fine-mapping "
            "method selection; solver/reference installation is separate."
        ),
        "locuszoom": (
            "Draw LocusZoom-style regional GWAS association plots with locuszoomr, Ensembl "
            "gene tracks and LDlink LD. Select a gene, region or lead SNP and optionally "
            "add fine-mapping posterior probabilities and credible sets. Use for regional "
            "association/credible-set figures; requires R and the user's LDlink API token."
        ),
    }
    if skill_id in descriptions:
        local_meta["description"] = descriptions[skill_id]
    local_meta["description"] = local_meta["description"].replace("<", "(").replace(">", ")")
    # Recompute the frontmatter boundary after portable Markdown-link edits.
    _, end = metadata(blob)
    blob = (
        b"---\n"
        + yaml.safe_dump(local_meta, sort_keys=False, allow_unicode=True).encode()
        + b"---\n"
        + blob[end:]
    )
    _, end = metadata(blob)
    intro = (
        b"\n## DNAKit execution\n\nRead [DNAKIT.md](DNAKIT.md) first for actual runtime, "
        b"resource paths and validation status.\n\n"
    )
    files["SKILL.md"] = (blob[:end] + intro + blob[end:], mode)
    if skill_id == "finemap":
        name = "scripts/finemap.py"
        old = b'here.parent.parent / "variant-annotate" / "scripts"'
        new = b'here.parent / "references" / "variant-annotate" / "scripts"'
        blob, mode = files[name]
        if blob.count(old) != 1:
            raise ValueError("Reviewed finemap portability edit no longer matches")
        files[name] = (blob.replace(old, new), mode)
        modifications.append(name)
    files["DNAKIT.md"] = (guide(skill_id, source).encode(), 0o644)
    portability = (
        "The optional variant-annotate dependency is retained inside this folder; "
        "its loader path is adapted."
        if skill_id == "finemap"
        else "Upstream scientific scripts are retained without algorithm changes."
    )
    notice = f"""# Upstream provenance

- Repository: https://github.com/{source["repository"]}
- Commit: `{source["revision"]}`
- Skill folder: `{folder}`
- Attribution: {source["attribution"]}
- Content license: {source["content_license"]}; see [LICENSE](LICENSE).
- Import date: {DATE}
- Archive SHA-256: `{source["archive_sha256"]}`

All skill-folder files are retained, including scripts, data and relative references.
DNAKit adds its execution/provenance notes and links them from the entrypoint.
Original and packaged hashes are recorded in registry/genomics-new-skills.json.
{portability}
Underlying software, model weights, databases, and services retain their own terms.
Cite the upstream project at this commit and the methods, data and papers actually used;
retaining a workflow is not validation of its scientific or clinical conclusions.
"""
    files["UPSTREAM.md"] = (notice.encode(), 0o644)
    category, summary, modules, commands, environment = RUNTIMES[skill_id]
    for name, (blob, _) in files.items():
        if name.endswith(".py"):
            ast.parse(blob.decode(), filename=f"{skill_id}/{name}")
    record = {
        "id": skill_id,
        "role": "dna_workflow",
        "description": local_meta["description"],
        "summary_zh": summary,
        "category": category,
        "capability": skill_id,
        "source_id": source["id"],
        "license": source["content_license"],
        "license_basis": (
            "Pinned repository LICENSE retained; underlying tools/data have separate terms."
        ),
        "upstream_license_field": meta.get("license"),
        "upstream_metadata": meta,
        "upstream_name": skill_id,
        "upstream_url": f"https://github.com/{source['repository']}/tree/{source['revision']}/{folder}",
        "upstream_path": folder,
        "upstream_revision": source["revision"],
        "attribution": source["attribution"],
        "citation_guidance": (
            "Cite the pinned upstream project and underlying methods/data actually used; "
            "see UPSTREAM.md and retained citation files."
        ),
        "aliases": [],
        "alternatives": [],
        "keywords": [category, summary],
        "runtime": {
            "python_modules": modules,
            "commands": commands,
            "environment": environment,
            "manual_checks": [
                "Read DNAKIT.md and SKILL.md for optional dependencies, versions, models, "
                "reference panels and service access."
            ],
        },
        "required_skills": [],
        "validation": {
            "metadata": "passed",
            "resources": "checksummed",
            "python_syntax": "passed",
            "execution": "not_run",
        },
        "modifications": modifications,
        "excluded_upstream_resources": {},
        "upstream_files": originals,
        "files": [
            {"path": name, "sha256": digest(blob), "bytes": len(blob), "mode": mode}
            for name, (blob, mode) in sorted(files.items())
        ],
    }
    return record, files


def check(root: Path) -> None:
    manifest = json.loads((root / MANIFEST).read_text())
    records = manifest["skills"]
    expected_sources = {s["id"]: s for s in SOURCES}
    actual_sources = {s["id"]: s for s in manifest["sources"]}
    if len(actual_sources) != len(manifest["sources"]) or set(actual_sources) != set(
        expected_sources
    ):
        raise ValueError("Unexpected or duplicate pinned source identities")
    for source_id, expected in expected_sources.items():
        for key in ("repository", "revision", "archive_sha256", "content_license", "attribution"):
            if actual_sources[source_id].get(key) != expected[key]:
                raise ValueError(f"Pinned source metadata mismatch: {source_id}/{key}")
    seen = set()
    for i, record in enumerate(records, 1):
        skill_id = record["id"]
        if skill_id in seen or skill_id not in RUNTIMES:
            raise ValueError(f"Unexpected or duplicate skill: {skill_id}")
        seen.add(skill_id)
        source = expected_sources.get(record["source_id"])
        if source is None or record["upstream_revision"] != source["revision"]:
            raise ValueError(f"Unpinned source: {skill_id}")
        expected_url = (
            f"https://github.com/{source['repository']}/tree/"
            f"{source['revision']}/{record['upstream_path']}"
        )
        if (
            record["license"] != source["content_license"]
            or record["attribution"] != source["attribution"]
            or record["upstream_url"] != expected_url
        ):
            raise ValueError(f"Pinned skill provenance mismatch: {skill_id}")
        folder = root / "skills" / skill_id
        if folder.is_symlink():
            raise ValueError(f"Unsafe skill directory: {folder}")
        expected = {item["path"] for item in record["files"]}
        actual = {
            str(p.relative_to(folder))
            for p in folder.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts
        }
        if expected != actual:
            raise ValueError(f"Unexpected/missing resources in {skill_id}: {expected ^ actual}")
        for item in record["files"]:
            relative = PurePosixPath(item["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe manifest resource: {item['path']}")
            path = folder / item["path"]
            if path.is_symlink() or not path.resolve().is_relative_to(folder.resolve()):
                raise ValueError(f"Unsafe resource: {path}")
            blob = path.read_bytes()
            if digest(blob) != item["sha256"] or len(blob) != item["bytes"]:
                raise ValueError(f"Resource checksum mismatch: {path}")
            if path.stat().st_mode & 0o777 != item["mode"]:
                raise ValueError(f"Resource mode mismatch: {path}")
            if path.suffix == ".py":
                ast.parse(blob.decode(), filename=str(path))
        meta, _ = metadata((folder / "SKILL.md").read_bytes())
        if meta["name"] != skill_id:
            raise ValueError(f"Skill name mismatch: {skill_id}")
        progress(i, len(records), f"checked {skill_id}")
    if seen != set(RUNTIMES):
        raise ValueError(f"Missing reviewed skills: {set(RUNTIMES) - seen}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Check installed files offline; no writes."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="Destination DNAKit checkout or empty review directory.",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=ROOT / "temp/genomics-skills-pinned",
        help="Commit-pinned tarball cache.",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    if args.check:
        check(root)
        return 0
    manifest_path = root / MANIFEST
    if manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite {manifest_path}; use --check.")
    plans = []
    for i, source in enumerate(SOURCES, 1):
        upstream = source_files(source, args.cache)
        prefix = source["prefix"]
        entries = [
            name
            for name in upstream
            if name.endswith("/SKILL.md") and (not prefix or name.startswith(prefix + "/"))
        ]
        for entry in sorted(entries):
            folder = entry.rsplit("/", 1)[0]
            record, files = make_skill(source, folder, upstream)
            target = root / "skills" / record["id"]
            if target.exists():
                raise FileExistsError(f"Refusing to overwrite existing skill: {target}")
            plans.append((record, files, target))
        progress(i, len(SOURCES), f"verified pinned {source['repository']}")
    ids = [record["id"] for record, _, _ in plans]
    if len(ids) != len(set(ids)) or set(ids) != set(RUNTIMES):
        raise ValueError("Imported set differs from reviewed skill selection")
    for i, (record, files, folder) in enumerate(plans, 1):
        folder.mkdir(parents=True, exist_ok=False)
        for name, (blob, mode) in sorted(files.items()):
            target = folder / name
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as handle:
                handle.write(blob)
            target.chmod(mode)
        progress(i, len(plans), f"imported {record['id']}")
    sources = [
        {
            **source,
            "cache_path": "temp/genomics-skills-pinned",
            "handling": "curated-import",
            "policy": (
                "Reviewed DNA workflows; complete skill resources and license notices retained."
            ),
            "citation_guidance": (
                "Cite pinned project and the scientific methods and data actually used."
            ),
        }
        for source in SOURCES
    ]
    manifest = {
        "schema_version": 1,
        "snapshot_date": DATE,
        "sources": sources,
        "skills": [record for record, _, _ in sorted(plans, key=lambda p: p[0]["id"])],
        "excluded_candidates": [
            {
                "source_id": "biomcp",
                "path": "sdlc/skills/biomcp-research",
                "reason": (
                    "Repository-specific developer issue-filing workflow; genomic research "
                    "is covered by the portable biomcp skill."
                ),
            }
        ],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    check(root)
    print(f"Imported {len(plans)} skills; merge {manifest_path} through the main catalog workflow.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
