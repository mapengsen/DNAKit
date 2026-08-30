# Three-dimensional structure and mechanical properties

Read DNA 3D structure files and analyze coordinate geometry, NMR multi-model flexibility, and 3DNA/DSSR structural parameters.

Three-dimensional properties must come from explicit coordinates or trajectories, without forging unique three-dimensional structures from ordinary sequences. The following examples are run from the project root directory, and the results are actually calculated from the current local source code and downloaded samples.

For the paper sources and internal formulas of various calculations, see [FAQ: Calculation Basis and References for Three-Dimensional Structure and Mechanical Properties ](../../faq.md#structure3d-references).

## 1) Read a single PDB model

- **Function:** Read a DNA model specified in the PDB, extract atomic coordinates, residues, chains and deduced sequences, and generate structural objects used in subsequent three-dimensional geometry calculations.
- **Calculation method:** Parse `HEADER`, `TITLE`, `MODEL` and `ATOM` according to legacy PDB fixed column specifications, and only retain canonical DNA residues; alternative conformations only accept blanks or `A`, repeat atoms are prioritized according to the main conformation, and then selected by occupancy, and finally the A/C/G/T sequence is mapped in residue order. This is file parsing, not structure prediction.
- **API:** `dnakit.structure3d.load_pdb(path[required], model_index[optional])`
- **Input:** Required human-readable legacy PDB path; optional positive integer model number.
- **Sample code:**

```python
from pathlib import Path

from dnakit.structure3d import load_pdb

structure = load_pdb(Path("temp/dna_structures/1BNA.pdb"))
print(structure.pdb_id, structure.model_index)
print(len(structure.atoms), len(structure.residues), structure.chain_ids)
print(structure.sequence_by_chain)
```

- **Example results:**

```text
1BNA 1
486 24 ('A', 'B')
('CGCGAATTCGCG', 'CGCGAATTCGCG')
```

## 2) Read PDB multi-model collection

- **Function:** Read multiple DNA conformations from the same PDB or model collection, unify atomic correspondences, and serve as input for RMSF and conformational change comparisons.
- **Calculation method:** Use the same fixed column parsing rules as single model, group by `MODEL` number and return all models. This step only establishes chain/residue/atom identities and does not perform rotation, translation alignment or RMSF calculations.
- **API:** `dnakit.structure3d.load_pdb_ensemble(path[required])`
- **Input:** Required The path to the legacy PDB containing one or more DNA models.
- **Sample code:**

```python
from pathlib import Path

from dnakit.structure3d import load_pdb_ensemble

models = load_pdb_ensemble(Path("temp/dna_structures/1AC7.pdb"))
print(len(models), models[0].pdb_id)
print(len(models[0].residues), models[0].sequence_by_chain)
```

- **Example results:**

```text
10 1AC7
16 ('ATCCTAGTTATAGGAT',)
```

## 3) Explicit coordinate geometry analysis

- **Function:** Calculate indicators such as radius of gyration, SASA, volume, shape, backbone dihedral angle, and geometric hydrogen bonds from explicit DNA atomic coordinates to quantitatively describe a three-dimensional conformation.
- **Calculation method:** Center of mass, radius of gyration, and tensor eigenvalues are calculated from atomic-mass-weighted coordinates; SASA is summed from deterministic sampling-point exposure proportions on the `vdW radius + probe` sphere, and volume is estimated from the union of voxels occupied by vdW spheres. Backbone angles use four-atom dihedrals, and hydrogen bonds are filtered with explicit distance/angle thresholds. The double-stranded helix metric is a DNAKit approximation based on paired `C1′` centers and global axes.
- **API:** `dnakit.structure3d.analyze_structure(structure[required], sasa_probe_radius_angstrom[optional], sasa_points_per_atom[optional], volume_grid_spacing_angstrom[optional], progress[optional])`
- **Input:** Required `DNA3DStructure`; optional 0–5 Å SASA probe, 24–4096 sample points per atom, 0.25–5 Å voxel spacing, and progress callbacks.
- **Sample code:**

```python
from pathlib import Path

from dnakit.structure3d import analyze_structure, load_pdb

structure = load_pdb(Path("temp/dna_structures/1BNA.pdb"))
result = analyze_structure(
    structure,
    sasa_points_per_atom=24,
    volume_grid_spacing_angstrom=1.5,
)
assert result.helix is not None
print(result.atom_count, result.residue_count)
print(
    round(result.shape.radius_of_gyration_angstrom, 3),
    round(result.solvent_accessible_surface_area_angstrom2, 3),
    round(result.molecular_volume_angstrom3, 3),
)
print(
    round(result.helix.mean_rise_angstrom, 3),
    round(result.helix.mean_twist_degree, 3),
)
```

- **Example results:**

```text
486 24
13.228 4650.051 4890.375
3.351 35.797
```

## 4) NMR/Multi-model RMSF Flexible

- **Function:** Calculate atom-by-atom and residue-by-residue RMSF across multiple 3D models with consistent atomic correspondence, and quantify the degree of fluctuation of different positions in the conformation ensemble.
- **Calculation method:** First take the atomic bonds common to all models, and subtract the geometric center of the shared atoms for each model to remove translation; calculate `RMSF_i = sqrt[Σm |r_im − <r_i>|²/M]` for atoms `i`, and then perform root mean square aggregation of the residue RMSF on its atomic RMSF. Implementation without Kabsch rotation fitting.
- **API:** `dnakit.structure3d.analyze_ensemble_flexibility(structures[required])`
- **Input:** Required 2–10000 `DNA3DStructure` models with at least 3 common atoms.
- **Sample code:**

```python
from pathlib import Path

from dnakit.structure3d import analyze_ensemble_flexibility, load_pdb_ensemble

models = load_pdb_ensemble(Path("temp/dna_structures/1AC7.pdb"))
result = analyze_ensemble_flexibility(models)
print(result.model_count, result.common_atom_count)
print(
    round(result.mean_atomic_rmsf_angstrom, 3),
    round(result.max_atomic_rmsf_angstrom, 3),
)
```

- **Example results:**

```text
10 510
0.986 2.281
```

## 5) 3DNA `bp_step.par` standard parameter analysis

- **Function:** Parse external 3DNA parameter files and summarize base pair steps and helix parameters such as shift, slide, rise, twist, etc. to facilitate structural statistics and comparison between models.
- **Calculation method:** Analyze each row using the 12 rigid-body parameter columns of 3DNA `bp_step.par`, checking the row count against the declared base-pair count; average rise and twist from the second row onward, then calculate `bp/turn = 360/mean twist` and `pitch = mean rise × bp/turn`. DNAKit does not recalculate the 3DNA local reference frame.
- **API:** `dnakit.structure3d.read_3dna_bp_step(path[required])`
- **Input:** Required An existing, UTF-8 decodable file that conforms to the 3DNA `bp_step.par` column structure.
- **Sample code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.structure3d import read_3dna_bp_step

with TemporaryDirectory() as directory:
    path = Path(directory) / "bp_step.par"
    path.write_text(
        "2 # base-pairs\n"
        "# Shear Stretch Stagger Buckle Prop-T Opening Shift Slide Rise Tilt Roll Twist\n"
        "A-T 0.1 0.2 0.3 1 2 3 0.4 0.5 3.4 4 5 36\n"
        "T-A 0.2 0.3 0.4 2 3 4 0.5 0.6 3.5 5 6 35\n",
        encoding="utf-8",
    )
    result = read_3dna_bp_step(path)
    print(result.base_pair_count)
    print(result.mean_rise_angstrom, result.mean_twist_degree)
    print(round(result.base_pairs_per_turn or 0.0, 4))
```

- **Example results:**

```text
2
3.5 35.0
10.2857
```

## 6) DSSR JSON summary parsing

- **Function:** Parse external DSSR JSON, extract base, pairing, helix, stem-loop and hydrogen bond annotations, and convert them into structured results that can be queried, summarized and reported.
- **Calculation method:** Read the program version and `nts`, `pairs`, `helices`, `stems`, `hairpins`, `hbonds` and `ntParams/ntPars` fields in the DSSR JSON, and return the summary after first verifying the declaration count and array length. Structure identification is done by an external DSSR that generates this JSON, DNAKit does not implement the DSSR algorithm.
- **API:** `dnakit.structure3d.read_dssr_json(path[required])`
- **Input:** Required Existing DSSR JSON output file.
- **Sample code:**

```python
from pathlib import Path

from dnakit.structure3d import read_dssr_json

result = read_dssr_json(Path("temp/dna_structures/1ehz-dssr-example.json"))
print(
    result.nucleotide_count,
    result.pair_count,
    result.helix_count,
    result.stem_count,
    result.hairpin_count,
)
print(result.hydrogen_bond_count, result.backbone_torsion_record_count)
```

- **Example results:**

```text
76 34 2 4 3
116 76
```

<span id="7"></span>
