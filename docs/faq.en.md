# FAQ

## What matching strategy is used for sequence function search? {#pattern-matching-strategy}

[The 10 types of functions in sequence function search](api/features/06_patterns.md) are all **deterministic rule algorithms**: when the algorithm version, input sequence and parameters are the same, the results must be the same. They do not use machine learning models, random sampling, or inference from experimental data. The PWM score is also a matching score calculated according to a fixed formula, not a model predicted probability.

These functions return candidate sequences that match the rules. Hits on promoter motifs, TF motifs, PAMs, or inverted repeats do not necessarily mean that the position must have true promoter activity, transcription factor binding, CRISPR editing efficiency, or hairpin structure.

### Universal matching convention

- **Strand direction:** Scan the forward strand directly; generate the reverse complementary sequence first and then scan the reverse strand; `strand="both"` scan both strands.
- **IUPAC Match:** Each symbol represents a set of possible bases, such as `N=ACGT`, `R=AG`, `Y=CT`, `W=AT`. When the base set represented by the target symbol and the rule symbol intersect, the position is compatible; the entire window is hit only when all positions are compatible.
- **Gap:** Explicit Gap is a scan boundary and candidates will not cross the Gap.
- **Circular sequence:** exact/IUPAC motif, PWM, promoter motif, restriction enzyme site and PAM can span the circular origin under supported conditions; canonical motif, codon and repeat structure scans cannot span the circular origin.
- **Resource upper limit:** `max_matches` Limit the number of returns. When the upper limit is reached, the result will be marked as truncated; when the number of comparisons or scanning units exceeds the corresponding upper limit, an error will be reported and calculations will not be silently omitted.

### 1. Common motif and PWM

- **exact:** Scan bit by bit with a window of the same length as the motif, and the window string must be exactly the same as the motif to hit.
- **IUPAC:** Also scan window by window, but each position is judged according to the IUPAC set compatibility rules above.
- **regex:** Use restricted DNA regular expression lookups; whether overlapping hits are allowed is controlled by `overlapping`.
- **PWM:** Compute a fixed log2-odds score for each window of the same length as the matrix: `sum over positions of log2(PWM probability / background probability)`. The score is greater than or equal to `threshold` to hit; the background default A, C, G, T are 0.25 each, and the target window containing ambiguous bases will be skipped.

The motif, canonical, PWM and threshold are provided by the caller and are not built-in biological function prediction models.

### 2. Start and stop codons

The program starts from the 0th, 1st, and 2nd bases of each selected strand respectively, and advances 3 bp each time to form three reading frames; a total of six reading frames are obtained when scanning both strands. Each triplet is a hit as long as it belongs to the start or stop codon set.

- NCBI Genetic Code Table 1: The start codon is `ATG`, and the stop codon is `TAA`, `TAG`, `TGA`.
- NCBI Genetic Code Table 11: The start codons are `ATG`, `GTG`, `TTG`, and the stop codons are `TAA`, `TAG`, `TGA`.
- When passing in `start_codons` or `stop_codons`, use the collection provided by the caller instead of the corresponding default collection.

This feature only reports individual codon positions, does not require an in-frame stop codon following the start codon, and is not equivalent to full ORF detection.

### 3. Promoter motif

The program treats each promoter consensus sequence as an IUPAC motif and matches it window by window on the selected strand. The only built-in rules are:

- Eukaryotic TATA box: `TATAWAWR`;
- Bacteria -Zone 10: `TATAAT`;
- Bacteria -Zone 35: `TTGACA`.

When `motifs` is passed in, the caller-supplied motif mapping is used instead of the built-in mapping. The program does not examine the distance between the -10/-35 regions, the transcription start site, or other regulatory elements and therefore only returns promoter motif candidates and does not predict promoter activity.

### 4.TF motif

The caller must provide the PWM and threshold of the transcription factor. The program scores windows by windows using the same log2-odds formula as the general-purpose PWM, returning only windows with scores greater than or equal to the threshold. DNAKit does not have built-in TF motif databases such as JASPAR, nor does it interpret scores as true binding strengths or binding probabilities.

### 5. Restriction enzyme sites

The program matches the restriction enzyme recognition sequence on both front and back strands according to IUPAC rules. After a hit, the cut point is calculated based on the upper and lower strand cutting offsets in the definition of the enzyme. The forward and reverse strand repeat results of palindromic recognition sites are combined.

The built-in mini directories include `BamHI`, `EcoRI`, `HaeIII`, `HindIII`, `NotI`, and `SmaI`; a custom `RestrictionEnzyme` can also be passed in. This rule does not include full REBASE and does not take into account the impact of methylation on digestion.

### 6. CRISPR PAM

The program first searches for the PAM bit by bit according to IUPAC rules on the selected chain, and then extracts a fixed-length guide from the adjacent position based on whether the PAM is located on the 3' side or the 5' side:

- `SpCas9`: `NGG`, PAM is located on the 3' side of the guide, the default guide length is 20 bp;
- `SaCas9`: `NNGRRT`, PAM is located on the 3' side of the guide, the default guide length is 21 bp;
- `AsCas12a`: `TTTV`, PAM is located on the 5' side of the guide, the default guide length is 20 bp.

The guide does not hit when it exceeds the sequence boundary; guides containing ambiguous bases are excluded by default. Then do deterministic filtering by `min_gc`, `max_gc` and `exclude_motifs`. Custom `PAMRule` can also be provided. This feature does not predict editing efficiency and off-target risk.

### 7. Palindrome sequence

The program enumerates each starting point and each length in the range `min_length` to `max_length`, comparing the candidate fragment bit by bit to its own reverse complement; only hits are made if all positions meet IUPAC compatibility rules. `maximal_per_start=True`, only the longest hit is retained for each starting point.

### 8. Repeat upside down

The program enumerates the starting point, arm length and loop length of the left arm, and then calculates the reverse complementary sequence of the left arm; it hits when the right arm and this sequence satisfy the IUPAC compatibility rules bit by bit. Currently no mismatch or indel is allowed, so it is only a sequence candidate for potential hairpin structures, no folding or free energy predictions are made.

### 9. Tandem repetition

The program enumerates repeating units from each starting point in ascending order of unit length, checking whether subsequent adjacent units are exactly the same as the first unit string, and extending to the first different location. After reaching the corresponding minimum number of repetitions, the smallest qualified unit from the starting point and its longest continuous repetition interval are reported.

The exact character equality rule applies here, and the IUPAC ambiguous symbols must also be the same character; no mismatch, indel, or break is allowed. When set to `overlapping=False`, scanning will continue from the end of the repeat interval after a hit; when set to `True`, scanning will continue from the next base.

### 10. Microsatellites

Microsatellite invokes the same tandem repeat algorithm, but fixes the repeat unit length at 1 to 6 bp. The default minimum number of repetitions is:

- 1 bp unit repeated at least 6 times;
- The 2-6 bp unit is repeated at least 3 times.

Matches must be consecutive and identical, no breaks, mismatches, or indels allowed. `min_repeats_by_unit` The threshold can be modified, but a threshold must be provided for each unit length from 1 to 6 bp.

## What calculation methods and references are available for Diversity and Novelty? {#diversity-novelty-references}

DNAKit retains its normalized-similarity methods and adds the paper's raw Levenshtein-distance methods:

| Metric | Default method | Second method |
| --- | --- | --- |
| Diversity | `calculation="similarity"`: distance is `1 - pair_similarity`; `score` is mean nearest-neighbor distance, with mean pair distance and threshold-cluster summaries also reported. | `calculation="levenshtein"`: `Σ(i≠j) Levenshtein(xᵢ,xⱼ) / [n(n−1)]`, equivalent to the mean raw edit distance over all unordered sequence pairs. |
| Novelty | `novelty_calculation="similarity"`: each query receives `1 - nearest_reference_similarity`. | `novelty_calculation="levenshtein"`: `meanᵢ minₛ Levenshtein(queryᵢ, referenceₛ)`. |

The second methods follow Cherednichenko & Poptsova, *Data augmentation with generative models improves detection of Non-B DNA structures*, **Computers in Biology and Medicine** 184 (2025) 109440, [DOI 10.1016/j.compbiomed.2024.109440](https://doi.org/10.1016/j.compbiomed.2024.109440). Section 2.8, equations (20) and (21), specifies Levenshtein distance. The terminology is inherited from Jain et al., *Biological Sequence Design with GFlowNets*, ICML 2022, [PMLR paper page](https://proceedings.mlr.press/v162/jain22a.html).

The article has an [official GitHub repository](https://github.com/powidla/nonB-DNA-structures-generation), with the relevant code in [`seq_analysis.ipynb`](https://github.com/powidla/nonB-DNA-structures-generation/blob/ea61a37f95c5a1effe64324af366c781755fe4c8/notebooks/seq_analysis.ipynb). As of 2026-09-02, the repository declares no open-source license, and its notebook uses fixed 100-bp preprocessing, flattened one-hot/KDTree operations, and chunk-local calculations that do not fully match the published equations. DNAKit therefore does not copy that code: it independently implements the published formulas with its bounded Levenshtein routine, performs no implicit padding or truncation, and does not claim exact reproduction of the notebook values in Table 2.

Levenshtein results are unnormalized counts of edit operations. Larger values mean greater diversity or greater distance from the reference library; compare datasets with similar sequence-length distributions.

## What are the calculation basis and references for physical and chemical properties? {#physicochemical-references}

The functions in [Physical and Chemical Properties](api/features/07_physicochemical.md) do not use machine learning models. They are divided into theoretical formulations, public empirical parametric models, DNAKit internal transparency rules, and external Primer3 thermodynamic structure predictions. The table below itemizes the actual basis; internal rules for unsourced papers are clearly marked, without supplementary packaging with non-existent citations.

| Function | Basis for calculation | References or sources |
| --- | --- | --- |
| THERMO-001 Molecular Weight | Anhydrous DNA residue mass summation, end and 5' phosphate correction | The current implementation records the standard oligonucleotide mass formula and includes a local Biopython numerical comparison, but this comparison cannot replace a primary scientific citation. |
| THERMO-014 260 nm extinction coefficient | Sum of adjacent dinucleotide hypochromicity parameters minus internal single base parameter | Warshaw & Tinoco, 1966, [DOI 10.1016/0022-2836(66)90115-X](https://doi.org/10.1016/0022-2836(66)90115-X); Cantor, Warshaw & Shapiro, 1970, [DOI 10.1002/bip.1970.360090909](https://doi.org/10.1002/bip.1970.360090909). |
| THERMO-002 Tm (Wallace) | `2 × (A+T) + 4 × (G+C)` Rule of thumb for short oligonucleotides | Wallace et al., 1979, [DOI 10.1093/nar/6.11.3543](https://doi.org/10.1093/nar/6.11.3543). The current implementation does not write this DOI into the result provenance, and the FAQ is here to supplement the literature relationship. |
| THERMO-002 Tm (nearest-neighbor) | Neighbor packing, terminal, symmetry, salt concentration and chain concentration models | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). |
| THERMO-003 salt concentration correction | SantaLucia unit price salt entropy correction | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). |
| THERMO-012 local Tm | Sliding window repeatedly calls Wallace or nearest-neighbor Tm | Does not introduce a new scientific model; references inherit the chosen Wallace 1979 or SantaLucia 1998 method. |
| EVAL-013 Synthesis Risk | Five-component equal weighted average of GC, isobase contiguous, tandem repeats, and inverted repeats | DNAKit internal transparent heuristic rules; no external paper or vendor rule sets, scores are not experimental synthesis success probabilities. |
| THERMO-004 Thermodynamic Parameters | ΔH, ΔS, ΔG and Tm for a perfectly complementary duplex | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). |
| THERMO-005 Nearest-neighbor | Same as THERMO-004, but returns stack-by-stack step details | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). |
| THERMO-006 Duplex stability | Native full complementation path using SantaLucia; optional mismatch/dangling path calling Primer3 | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460); Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596). The stable Boolean criterion `Tm > configured temperature` is a DNAKit result-interpretation rule. |
| THERMO-007 Base stacking | Query ΔH/ΔS parameters and calculate ΔG by adjacent dinucleotide | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). |
| THERMO-008 Hairpin | User-installed Primer3 `ntthal` Thermodynamic Hairpin Structure Prediction | Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) and [Primer3 Official Manual ](https://primer3.org/manual.html); specific values also depend on the actual Primer3 version and parameter directory. |
| THERMO-009 Self-dimer | User-installed Primer3 `ntthal` thermodynamic self-dimer prediction | Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) and [Primer3 official manual ](https://primer3.org/manual.html); specific values also depend on the actual Primer3 version and parameter directory. |
| THERMO-010 Heterodimer | User-installed Primer3 `ntthal` thermodynamic heterodimer prediction | Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) and [Primer3 Official Manual ](https://primer3.org/manual.html); exact values also depend on the actual Primer3 version and parameter catalog. |

These references illustrate the sources of algorithms and parameters and do not imply that any calculated results are equivalent to experimental measurements. Tm, ΔG, and structural results are still affected by salt concentration, chain concentration, temperature, sequence length, chemical modifications, and the actual backend version.

## What are the calculation basis and references for double-chain thermodynamic expansion? {#duplex-thermodynamics-references}

[Double-stranded thermodynamics extension](api/features/19_duplex_thermodynamics.md) does not use machine learning. Native double-strand energies are derived from SantaLucia public parameters, other results are standard thermodynamic relations, DNAKit transparent combination rules, or external Primer3 calculations.

| Function | Actual calculation basis | References or sources |
| --- | --- | --- |
| 1. Perfectly complementary duplex `ΔH/ΔS/ΔG/Tm` | Adjacent stacking, termini, symmetry, monovalent salt and chain concentration formulas | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). This DOI has written native results provenance. |
| 2. Unified double-chain stability | `native` follows item 1; `primer3-cli` calls the mismatch/dangling-end structural model of `ntthal` | The native path refers to SantaLucia 1998; the Primer3 path refers to Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) and [Primer3 official manual](https://primer3.org/manual.html). `Tm > configured temperature` is a DNAKit result-interpretation rule and does not have a separate paper. |
| 3. Contribution of adjacent base pairs to steps | Check the adjacent stacking `ΔH/ΔS` table and calculate `ΔG=ΔH−TΔS/1000` | SantaLucia, 1998, [DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460). No additional references are made to the new model. |
| 4. Conditions and Na⁺/K⁺ | The condition object is recorded as a numerical checksum; the native model is corrected with `0.368(N−1)ln([Na⁺]+[K⁺])` `ΔS` | The salt correction refers to SantaLucia 1998. Object checksum `Na⁺+K⁺` field merging is a DNAKit implementation rule and does not have a separate paper. |
| 5. `Ka`, `Kd` and double-strand ratio | `Ka=exp(−ΔG/RT)`, `Kd=1/Ka`, and then solve the ideal two-state mass conservation equation | `ΔG` cited SantaLucia 1998; the equilibrium constant relationship and quadratic equation are standard thermodynamics and DNAKit transparent algebra implementation, not bound to separate papers or fitted data. |
| 6. Theoretical melting curve | Repeat items 1 and 5 at each temperature, and do linear interpolation for the intersection point of double-strand ratio 0.5 | Reference inheritance SantaLucia 1998; temperature scan and linear interpolation are DNAKit internal deterministic rules, no new paper. |
| 7. 5′/3′ end stability | Calculate the nearest neighbor `ΔG` for two equal-length end windows separately and then compare | Energy reference SantaLucia 1998; window truncation and "higher `ΔG` is the less stable end" are DNAKit internal rules. |
| 8. DMSO/formamide correction | Linear empirical summation formula of Primer3 manual | [Primer3 official manual ](https://primer3.org/manual.html). The default DMSO factor of 0.6 is quoted in the manual from Musielski et al., 1981; the formamide formula is quoted from Blake & Delcourt, 1996, [DOI 10.1093/nar/24.11.2095](https://doi.org/10.1093/nar/24.11.2095). This function is not a mechanistic free energy model. |
| 9. Primer3 CLI extension | `oligotm` calculates Tm, `ntthal` calculates hairpins and dimers; DNAKit only does adapter | Untergasser et al., 2012, [DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) and [Primer3 official manual ](https://primer3.org/manual.html). The actual Primer3 version and thermodynamic parameter catalog should also be documented. |

## What are the calculation basis and references for secondary structure properties? {#secondary-structure-references}

[Secondary structure properties](api/features/20_secondary_structure.md) Contains DNAKit internal parsing/formulas and conditional NUPACK adapter. Structural energy and equilibrium results are NUPACK calculations only if items 6 and 7 are run explicitly.

| Function | Actual calculation basis | References or sources |
| --- | --- | --- |
| 1. Dot-bracket parsing | Stack pairing, continuous nested pair merging and fixed structure classification | DNAKit internal deterministic parser, no paper parameters or training models. For the symbol convention of Dot-parens-plus, please refer to Zadeh et al., 2011, [DOI 10.1002/jcc.21596](https://doi.org/10.1002/jcc.21596); the extended bracket parsing of `() [] {} <>` is implemented by DNAKit. |
| 2. Pairing probability and window accessibility | Marginal probability derivation of NUPACK style dense probability matrix | For matrix semantics, please refer to [NUPACK analysis document ](https://docs.nupack.org/analysis/) and Zadeh et al. 2011. Symmetry/rowsum checks and "arithmetic mean of window edge probabilities" are DNAKit internal rules, not NUPACK joint open probabilities. |
| 3. Normalized ensemble defect | The expected proportion of each base that does not adopt the target pairing state | Zadeh, Wolfe & Pierce, 2011, [DOI 10.1002/jcc.21633](https://doi.org/10.1002/jcc.21633). DNAKit only evaluates this definition for the provided probability matrix. |
| 4. Target structure probability | `exp[−(Gtarget−Gensemble)/(RT)]` Boltzmann relation | Standard statistical thermodynamic formula, DNAKit is not bound to a separate paper or fitting parameters. If the free energy comes from NUPACK, the actual NUPACK model and version should be quoted. |
| 5. NUPACK passive detection | Python package location and metadata check | Engineering environment check, does not perform NUPACK, so no scientific paper citations; the results must not be written as NUPACK predictions. |
| 6. NUPACK single complex | Partition functions, MFE, pairing probabilities, suboptimal structures, and Boltzmann sampling for external NUPACK 4 | DNAKit provenance Current record Zadeh et al. 2011, [DOI 10.1002/jcc.21596](https://doi.org/10.1002/jcc.21596). The NUPACK 4 official citation page also specifies Fornace, Porubsky & Pierce, 2020, [DOI 10.1021/acssynbio.9b00523](https://doi.org/10.1021/acssynbio.9b00523); citations should be selected according to the actual calling function and version, see [NUPACK official citation page ](https://docs.nupack.org/). |
| 7. NUPACK tube equilibration | Multi-chain complex and tube equilibration analysis of external NUPACK 4 | Dirks et al., 2007, [DOI 10.1137/060651100](https://doi.org/10.1137/060651100); Fornace, Porubsky & Pierce, 2020, [DOI 10.1021/acssynbio.9b00523](https://doi.org/10.1021/acssynbio.9b00523); and [NUPACK official citation page ](https://docs.nupack.org/). DNAKit's target/non-target ratio is an internal summation of the concentrations returned by the adapter. |

## What are the calculation basis and references for three-dimensional structure and mechanical properties? {#structure3d-references}

[Three-dimensional structure and mechanical properties](api/features/21_structure3d.md)Do not predict three-dimensional structures from DNA sequences. Items 1–4 are deterministic analytical/geometric calculations of explicit coordinates, items 5 and 6 are adapters to external 3DNA/DSSR results.

| Function | Actual calculation basis | References or sources |
| --- | --- | --- |
| 1. Read a single PDB model | legacy PDB `MODEL/ATOM` fixed columns, DNA residue table and alternative conformation selection | [wwPDB PDB v3.3 coordinate record specification ](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html). This is a file format source, not a scientific prediction paper. |
| 2. Read PDB multi-model | Same as item 1, then group by `MODEL` sequence number | [wwPDB PDB v3.3 coordinate record specification ](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html). Multi-model organizations are format rules and there are no individual predictive models. |
| 3. Explicit Coordinate Geometry | Center of mass, radius of gyration, tensor eigenvalues, spherical point sampling SASA, voxel volume, dihedral angles, distance/angle hydrogen bond filtering and global helix axis approximation | Current results provenance marks the entire native geometry implementation as a DNAKit internal method, unbound to the unity paper. Among them, SASA is the Shrake–Rupley style spherical point exposure method, please refer to Shrake & Rupley, 1973, [DOI 10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9); voxel volume, hydrogen bond threshold and `C1′` global axis are DNAKit specific rules and should not be written as 3DNA/DSSR standard parameters. |
| 4. NMR/Multi-model RMSF | Square root of mean square displacement after centering of common atomic translations | DNAKit internal implementation of the standard RMSF mathematical definition, not bound to a separate paper. The implementation does not do a Kabsch rotation fit and does not infer persistence lengths or mechanical moduli. |
| 5. 3DNA `bp_step.par` Analysis | Read the 12 rigid body parameters calculated by 3DNA. DNAKit only finds the average, `360/twist` and pitch | Lu & Olson, 2003, [DOI 10.1093/nar/gkg680](https://doi.org/10.1093/nar/gkg680) and [ Official description of 3DNA parameters ](https://x3dna.org/highlights/schematic-diagrams-of-base-pair-parameters). |
| 6. DSSR JSON summary | Read DSSR identified nucleotide, pair, helix, stem, hairpin, and hydrogen bond counts | Lu, Bussemaker & Olson, 2015, [DOI 10.1093/nar/gkv716](https://doi.org/10.1093/nar/gkv716) and [DSSR JSON official description ](https://x3dna.org/highlights/dssr-output-in-json-format). DSSR completes structure recognition, and DNAKit only parses existing JSON. |

## What are the calculation basis and references for optical and concentration conversion? {#optics-concentration-references}

[Optics and concentration conversion](api/features/18_optics_concentration.md) does not use machine learning. It uses the public single-stranded `ε260` empirical parameters, Beer–Lambert's law, molecular weight and unit conversions, and double-stranded/modified parameters explicitly provided by the caller.

| Function | Actual calculation basis | References or sources |
| --- | --- | --- |
| 1. Single-chain molar extinction coefficient at 260 nm | `Σ adjacent-dinucleotide coefficients − Σ internal single-base coefficients` | Warshaw & Tinoco, 1966, [DOI 10.1016/0022-2836(66)90115-X](https://doi.org/10.1016/0022-2836(66)90115-X); Cantor, Warshaw & Shapiro, 1970, [DOI 10.1002/bip.1970.360090909](https://doi.org/10.1002/bip.1970.360090909). Two articles have been written about the provenance of this result. |
| 2. Theoretical optical properties of single chain/double chain | Use item 1 for single chain; use `13200 × bp` average formula or `(ε1+ε2)(1−h)` for double chain; use the sum of anhydrous residues for molecular weight | Single chain `ε` cite the above two papers. `13200 M⁻¹·cm⁻¹/bp` is the traditional `1 OD260 ≈ 50 µg/mL` double-chain average conversion used in the implementation. The current source code is not bound to the original paper; `h`, `Δε` and `ΔMW` are provided by the caller, and their actual experiments, literature or manufacturer data sheets should be cited. The anhydrous residue molecular weight formula is not yet bound to the original paper. |
| 3. 1 OD260 corresponding to nmol/mass | `10⁶/ε` and `1000×MW/ε` | Directly derived from OD260 definition and unit conversion, no new paper; citations should inherit from selected `ε260` and `MW` sources. |
| 4. A260 to concentration and total amount | `c=A/(εl)`, plus dilution factor, `m=nMW`, `n=cV` and explicit dye subtraction | Beer–Lambert standard laws of physics with DNAKit transparent unit conversion, currently implemented unbound to a specific modern paper. `ε`, `MW` and dye factors should each cite their actual source. |
| 5. Concentration/quantity/mass interchange of substances | `m=nMW`, `n=m/MW`, `c=n/V`, `n=cV` | DNAKit internal dimension conversion, not an empirical model, no separate paper. |
| 6. Dye and modification correction | `Σcount×Δε`, `Σcount×ΔMW`, and `Σlabel-peak×correction-factor` | DNAKit does not have a built-in dye parameter table and only does deterministic addition and subtraction, so there is no unified paper. The caller must document authentic experimental, literature, or manufacturer sources for each parameter. |

The above citations only indicate the sources of algorithms, parameters or external tools. For free energy, hypochromicity, dye correction, or 3DNA/DSSR/NUPACK outputs provided by the caller, the actual parameter table, software version, input conditions, and data source must also be documented.

## What are the references for deep-learning property prediction? {#deep-learning-property-prediction-references}

The following table covers the primary model papers, task datasets, and training protocols used by all 54 functions integrated in [deep-learning property prediction](api/features/23_deep_learning_property_prediction.md). Multiple output heads share one paper, so the same citation is not repeated for every function.

| Functions | Reference |
| --- | --- |
| RNA-seq, CAGE, PRO-cap, ATAC-seq, DNase-seq, ChIP-seq, splicing, and contact maps | Avsec et al., *Advancing regulatory variant effect prediction with AlphaGenome*, **Nature** 649, 1206–1218 (2026), [DOI 10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0). |
| Human and mouse regulatory tracks | Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x). |
| 18 NT Revised classification tasks | Backbone: Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x); task definitions: Dalla-Torre et al., *Nucleotide Transformer: building and evaluating robust foundation models for human genomics*, **Nature Methods** 22, 287–297 (2025), [DOI 10.1038/s41592-024-02523-z](https://doi.org/10.1038/s41592-024-02523-z), and the [revised dataset card](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised); full-fine-tuning protocol: Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272), Appendix C.4. |
| 9 Genomic Benchmarks classification tasks | Backbone: Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x); task datasets: Grešová et al., *Genomic benchmarks: a collection of datasets for genomic sequence classification*, **BMC Genomic Data** 24, 25 (2023), [DOI 10.1186/s12863-023-01123-8](https://doi.org/10.1186/s12863-023-01123-8); full-fine-tuning protocol: Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272), Appendix C.4. |
| Foundation encoder for 14-class single-nucleotide genome segmentation | Dalla-Torre et al., *Nucleotide Transformer: building and evaluating robust foundation models for human genomics*, **Nature Methods** 22, 287–297 (2025), [DOI 10.1038/s41592-024-02523-z](https://doi.org/10.1038/s41592-024-02523-z). |
| Trained 14-class single-nucleotide genome segmentation head | de Almeida et al., *Annotating the genome at single-nucleotide resolution with DNA foundation models*, **Nature Methods** (2025), [DOI 10.1038/s41592-025-02881-2](https://doi.org/10.1038/s41592-025-02881-2). |
| Long-context zero-shot variant effects and exon probability | Brixi et al., *Genome modelling and design across all domains of life with Evo 2*, **Nature** (2026), [DOI 10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5). |
| Allele-conditional-probability variant effects | Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272); Li et al., *GENERator-v2: Reconciling Coarse Tokenization with Single-Nucleotide Resolution in Genomic Language Modeling*, **bioRxiv** (2026), [DOI 10.64898/2026.01.27.702015](https://doi.org/10.64898/2026.01.27.702015). |
| Central dogma, taxonomy, species, protein localization/stability, ncRNA family, and pair-interaction tasks | He et al., *Generalized biological foundation model with unified nucleic acid and protein language*, **Nature Machine Intelligence** 7, 942–953 (2025), [DOI 10.1038/s42256-025-01044-4](https://doi.org/10.1038/s42256-025-01044-4). |

Direct inference also depends on the following released code or trained weights:

- The 11 genomic-track tasks: [AlphaGenome research](https://github.com/google-deepmind/alphagenome_research) and [`google/alphagenome-all-folds`](https://huggingface.co/google/alphagenome-all-folds).
- Human/mouse regulatory tracks: the [official Enformer implementation](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) and [`EleutherAI/enformer-official-rough`](https://huggingface.co/EleutherAI/enformer-official-rough).
- The 27 task-classification checkpoints are downloaded from the [shared Google Drive checkpoint folder](https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing); task definitions come from the [NT Revised dataset](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised) and [Genomic Benchmarks repository](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks), and DNAKit does not redistribute these weights inside the package.
- Single-nucleotide genome segmentation: the [Nucleotide Transformer repository](https://github.com/instadeepai/nucleotide-transformer) and [`InstaDeepAI/segment_nt`](https://huggingface.co/InstaDeepAI/segment_nt).
- Zero-shot variant effects/exon probability: the [Evo 2 repository](https://github.com/ArcInstitute/evo2), [`arcinstitute/evo2_7b`](https://huggingface.co/arcinstitute/evo2_7b), [`arcinstitute/evo2_7b_base`](https://huggingface.co/arcinstitute/evo2_7b_base), and [`schmojo/evo2-exon-classifier`](https://huggingface.co/schmojo/evo2-exon-classifier).
- Allele-conditional-probability variant effects: the [GENERator repository](https://github.com/GenerTeam/GENERator) and [`GenerTeam/GENERator-v2-eukaryote-1.2b-base`](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base).
- Nucleic-acid/protein downstream tasks: [LucaOne](https://github.com/LucaOne/LucaOne), [LucaOneTasks](https://github.com/LucaOne/LucaOneTasks), and [Zenodo 10.5281/zenodo.15171943](https://doi.org/10.5281/zenodo.15171943).

These citations identify the model, trained task-head, and checkpoint sources. They do not make a prediction equivalent to an experimental measurement or a clinical conclusion.
