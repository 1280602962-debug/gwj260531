# Pocket mechanism verdict v1 — why PIK3CA-end receptor swap collapses summary_min

> Script: `pocket_superposition_v1.py`; raw output: `pocket_superposition_v1_output.txt`.
> Inputs: already-committed crystal coordinates only (`*_protein.pdb`, `*_crystal.pdb`/cognate ligand files
> under `jcim_structure_robust_v0/` and `pik3ca_mtor_panel48_v0/tables/`). **No new docking.**
> Question: `STRUCTURE_ROBUSTNESS_VERDICT_V1.md` showed that swapping PIK3CA (4L23→4JPS/5DXT) collapses
> PM48 `summary_min` from 0.692 to 0.486/0.505, while swapping mTOR (4JT6→4JSX) only mildly reduces it
> (0.639, Δ≈−0.05). Is this asymmetry explained by a real structural difference at the ATP pocket, and if
> so, at what resolution (whole-domain vs local pocket)?

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

Zero residue-identity mismatches at matched positions: these are genuinely the same protein (no
crystallization mutants/tags at pocket-adjacent positions detected), so the differences below are
conformational, not sequence-level.

Pocket residues used (heavy-atom ≤5 Å of the reference cognate ligand):
- PIK3CA (4L23, cognate X6K): Met772, Trp780, Ile800, Lys802, Leu807, Asp810, Leu814, Tyr836, Cys838,
  Ile848, Glu849, Val850, Val851, Ser854, Thr856, Gln859, Met922, Phe930, Ile932, Asp933.
- mTOR (4JT6, cognate X6K/PI-103): Ile2163, Pro2169, Leu2185, Lys2187, Glu2190, Leu2192, Asp2195,
  Tyr2225, Val2227, Ile2237, Gly2238, Trp2239, Val2240, Met2345, Leu2354, Ile2356, Asp2357, Phe2358.

## Interpretation

1. **Whole-domain conformational variability across crystal forms is markedly higher for PIK3CA
   (1.44–1.49 Å) than for mTOR (0.45 Å) in this structure set.** This alone is consistent with, and
   quantitatively explains, the asymmetry already observed at the score level: swapping the PIK3CA end
   collapses `summary_min` (Δ ≈ −0.19 to −0.21), while swapping the mTOR end barely moves it (Δ ≈ −0.05).
   The "receptor dependence" reported in `STRUCTURE_ROBUSTNESS_VERDICT_V1.md` is not an isolated docking
   artifact; it tracks a real, measurable difference in how much these deposited PIK3CA structures differ
   from one another relative to how much these mTOR structures differ.
2. **Local pocket Cα geometry is not always the limiting factor.** For 5DXT, the local pocket Cα RMSD
   (0.343 Å) is *smaller* than the global RMSD (1.441 Å) — the ATP-site backbone itself is well conserved
   even though the rest of the domain has moved substantially — yet `summary_min` on 5DXT still collapsed
   to 0.505, essentially matching 4JPS (where the pocket itself is comparably or more perturbed, local
   RMSD 0.867 Å). This means **Cα-level pocket conservation alone does not guarantee that pocket-matched
   discrimination transfers**; side-chain rotamers, protonation, or docking search-space sensitivity not
   captured by a Cα-only metric likely also contribute. We do not have PLIF-level or rotamer-level evidence
   for this round and do not claim it; it is flagged as the natural next step, not resolved here.
3. **Cognate ligands occupy the same general site across crystal forms** (centroid distances 2.1–2.6 Å for
   all three alternates, all inside a single ATP-competitive pocket) — the AUROC collapse is not explained
   by docking into a grossly different, unrelated site.

## Claim implication

The manuscript may state, with this evidence: *the PIK3CA end of the PM pair shows substantially larger
inter-crystal-form Cα variability than the mTOR end (1.44–1.49 Å vs 0.45 Å in this structure set), which is
consistent with its greater sensitivity of pocket-matched discrimination to receptor choice; however, local
pocket Cα conservation (as in 5DXT) is not sufficient to preserve the discrimination, indicating that
finer-grained (side-chain/rotamer or search-space) factors beyond Cα geometry are also at play and were not
resolved in this round.* Do **not** claim a fully solved mechanism or a validated general rule from n=2
alternates on one target end.
