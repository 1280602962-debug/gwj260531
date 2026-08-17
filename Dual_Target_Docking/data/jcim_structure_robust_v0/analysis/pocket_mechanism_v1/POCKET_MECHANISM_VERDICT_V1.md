# Pocket mechanism verdict v1 — why PIK3CA-end receptor swap collapses summary_min

> Script: `pocket_superposition_v1.py`; raw output: `pocket_superposition_v1_output.txt`.
> Inputs: already-committed crystal coordinates only (`*_protein.pdb`, `*_crystal.pdb`/cognate ligand files
> under `jcim_structure_robust_v0/` and `pik3ca_mtor_panel48_v0/tables/`). **No new docking.**
> Question: `STRUCTURE_ROBUSTNESS_VERDICT_V1.md` showed that swapping PIK3CA (4L23→4JPS/5DXT) collapses
> PM48 `summary_min` from 0.692 to 0.486/0.505, while swapping mTOR (4JT6→4JSX) only mildly reduces it
> (0.639, Δ≈−0.05). Do these deposited structures differ in Cα geometry in a way that is *consistent with*
> that asymmetry, and at what resolution (whole-domain vs local pocket)? Consistency is not causation.

## Method

1. Extract the longest protein chain Cα coordinates from each receptor's `*_protein.pdb`.
2. Match Cα atoms between the frozen main-panel structure (4L23 or 4JT6) and each alternate (4JPS, 5DXT,
   4JSX) by residue number **and** residue name (exact matches only; any mismatch is reported, none found).
3. Kabsch-superpose the alternate onto the main-panel structure using **all** matched Cα atoms
   (Bio.PDB `Superimposer`) → **global Cα RMSD** (one rigid-body fit for the whole matched chain).
4. Define the pocket-residue set from the main-panel structure's **own cognate ligand** (heavy-atom
   distance ≤ 5 Å to any residue atom, not just Cα) — X6K/PI-103 for 4L23, X6K for 4JT6.
5. Apply the **same** global transform (no separate local fit) to the pocket-residue Cα atoms only →
   **local pocket Cα RMSD**. This tests whether the pocket moves more or less than the rest of the domain
   under one rigid-body superposition, not whether it can be locally re-fit to look similar.
6. Apply the same transform to the alternate structure's own cognate ligand (1LT for 4JPS, 5H5 for 5DXT,
   17G/Torin2 for 4JSX) and report the **centroid distance** to the main-panel structure's own cognate
   ligand centroid — tests whether the two ligands occupy the same general site after aligning the
   proteins, independent of docking.

## Results

| ref (main panel) | alt (swap) | matched Cα | mismatched Cα | **global Cα RMSD (Å)** | pocket n | **local pocket Cα RMSD (Å)** | cognate-ligand centroid distance (Å) |
|---|---|---:|---:|---:|---:|---:|---:|
| 4L23 (PIK3CA) | 4JPS | 982 | 0 | **1.486** | 20 | **0.867** | 2.566 |
| 4L23 (PIK3CA) | 5DXT | 862 | 0 | **1.441** | 20 | **0.343** | 2.072 |
| 4JT6 (mTOR)   | 4JSX | 1054 | 0 | **0.454** | 18 | **0.467** | 2.196 |

Zero residue-identity mismatches at matched positions: compared residues are the same amino acids.
Unmatched positions were not compared (5DXT has 862 matched Cα vs 982 for 4JPS), so this is not a
full-construct mutant screen.

Pocket residues used (heavy-atom ≤5 Å of the reference cognate ligand):
- PIK3CA (4L23, cognate X6K): Met772, Trp780, Ile800, Lys802, Leu807, Asp810, Leu814, Tyr836, Cys838,
  Ile848, Glu849, Val850, Val851, Ser854, Thr856, Gln859, Met922, Phe930, Ile932, Asp933.
- mTOR (4JT6, cognate X6K/PI-103): Ile2163, Pro2169, Leu2185, Lys2187, Glu2190, Leu2192, Asp2195,
  Tyr2225, Val2227, Ile2237, Gly2238, Trp2239, Val2240, Met2345, Leu2354, Ile2356, Asp2357, Phe2358.

## Interpretation

1. **Whole-domain Cα variability is larger for these PIK3CA crystal forms (1.44–1.49 Å) than for these
   mTOR crystal forms (0.45 Å).** That difference is **consistent in direction** with the score-level
   asymmetry (PIK3CA-end swap Δ ≈ −0.19 to −0.21; mTOR-end swap Δ ≈ −0.05). It does **not** quantitatively
   explain or prove causation. Limits: n = 2 PIK3CA alternates and n = 1 mTOR alternate; 5DXT matched only
   862 Cα versus 982 for 4JPS, so the two global RMSDs are not equal-coverage comparisons. Receptor
   dependence remains an experimental docking result; this superposition is a descriptive structural
   correlate in this structure set, not a general crystal-form rule.

2. **Local pocket Cα geometry does not track the AUROC collapse.** For 5DXT, the local pocket Cα RMSD
   (0.343 Å) is *smaller* than the global RMSD (1.441 Å) — the ATP-site backbone itself is well conserved
   — yet `summary_min` on 5DXT still fell to 0.505, essentially matching 4JPS (local RMSD 0.867 Å,
   summary_min 0.486). **Cα-level pocket conservation alone does not guarantee that pocket-matched
   discrimination transfers.** Side-chain rotamers, protonation, or docking search-space sensitivity are
   untested; there is no PLIF- or rotamer-level evidence in this round.

3. **Cognate ligands occupy the same general site across crystal forms** (centroid distances 2.1–2.6 Å).
   This rules out docking into a grossly unrelated pocket. It does **not** prove that the cognate binding
   modes are identical.

Zero residue-identity mismatches at matched positions mean the compared residues are the same amino acids.
Unmatched residues (especially the extra 120 Cα on 4L23 vs 5DXT) were not compared and are not a mutant
screen of the full construct.

## Claim implication

The manuscript may state: *in this structure set, PIK3CA inter-crystal-form global Cα RMSD (1.44–1.49 Å)
is larger than mTOR (0.45 Å), consistent in direction with greater PIK3CA-end sensitivity of
pocket-matched discrimination; local pocket Cα conservation (5DXT 0.343 Å) is not sufficient to preserve
that discrimination.* Do **not** claim a quantitative mechanistic explanation, a solved mechanism, or a
validated general rule from n=2 / n=1 alternates.
