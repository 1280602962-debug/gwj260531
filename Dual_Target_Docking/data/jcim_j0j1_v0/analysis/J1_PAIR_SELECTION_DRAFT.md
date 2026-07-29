# J1 — Docking pair selection draft (select, do not dock)

**Status:** DRAFT — **awaiting user approval before any J2+ docking**  
**Inputs:** J0 supply audit + existing frozen receptors + `pdb_holo_counts.csv`  
**Machine-readable:** `../../protocols/PAIR_ROLES_DRAFT_J1.yaml`

## K=4 recommended roster

| Seat | Pair | Role | Strict supply | Prep / docking status | Why |
|------|------|------|---------------|----------------------|-----|
| 1 | **PIK3CA / mTOR** | development | Y (80/81) | panel48 LigPrep — **needs J2 RDKit** | Best supply + existing signal vs volume |
| 2 | **EGFR / HER2** | **case / supply-limited** | N (39/7) | EH110 unified RDKit **done** | Direction inversion + baseline fail; **not** a thick strict panel |
| 3 | **AChE / BChE** | development (new) | Y (189/78) | **not docked** — J3/J4 | Passes strict; conventional pockets; holo rich |
| 4 | **PIK3CA / PIK3CB** | development / isoform control | Y (56/67) | **not docked** — J3/J4 | Only remaining non-metal Y pair; “too close” narrative = feature |

### Tier T (not in K=4 sampling budget unless approved)

| Pair | min HN | Role |
|------|-------:|------|
| Mcl-1 / Bcl-xL | 12 | PPI pose-gold case; thin |
| PIK3CB / mTOR | 25 | thin pathway |
| HDAC1 / HDAC6 | 93 | **reject as primary** despite Y — Zn metal + isoform |

### Explicitly not selected

- NLRP3/JNK1 — private holdout  
- BRD4/HDAC*, JAK2/HDAC1 — metal + failed supply  
- Expanding EGFR for significance — forbidden by Stage M / route doc  

## Suggested panel N (when approved)

| Pair | Suggested N | Strict quota sketch | Structure candidates |
|------|------------:|---------------------|----------------------|
| PIK3CA/mTOR | 110 (expand from 48) or keep 48 + report | dual≥25, A≥20, B≥20, neither≥10 under **strict** | A: **4L23**; B: **4JT6** (frozen) |
| EGFR/HER2 | 110 (keep) | θ=6 panel OK; **do not** claim strict-thick | A: **3POZ**; B: **3RCD** |
| AChE/BChE | 100–120 | strict dual/A/B/neither quotas | A: 4EY4 / 6O5V; B: 6ZWI / 1P0I (freeze in J3) |
| PIK3CA/PIK3CB | 100–120 | strict quotas; isoform-matched boxes | A: 4L23; B: pick β holo in J3 cognate QC |

## Docking budget estimate (not executed)

| Task | Vina jobs (approx) | Notes |
|------|-------------------:|-------|
| J2 PM48 → RDKit | ~96 | mandatory for cross-pair prep unity |
| Optional PM expand → ~110 | ~220 | if strict quotas need more ligands |
| J3 receptor freeze AChE/BChE + PIK3CB | 0 | judgment + cognate QC |
| J4 AChE/BChE panel ~110 ×2 | ~220 | |
| J4 PIK3CA/PIK3CB ~110 ×2 | ~220 | PIK3CA end may reuse 4L23 |
| J5 GNINA CNN rescore all poses | 0 new sampling | ~800–1000 rescores |
| **Total new Vina if full plan** | **~540–760** | matches `JCIM_ROUTE_ASSESSMENT_V1` §4.3 |

**EGFR unified prep: 0 new jobs** (already assembled).

## Risks

| Risk | Mitigation |
|------|------------|
| Only 3 non-metal strict-Y pairs in 49 | K=4 must include EGFR as *case*, not as 4th thick pair |
| PIK3CA/PIK3CB too close | Frame as isoform control; report directional metrics separately |
| AChE/BChE homologous pockets | Expect weaker “end asymmetry” story; still valid supply test |
| ChEMBL fetch queue incomplete | Re-run J0 after API recovery before locking forever |

## Gate

**J1 = Draft complete.** No docking performed. Proceed to J2/J3 **only after explicit user approval.**
