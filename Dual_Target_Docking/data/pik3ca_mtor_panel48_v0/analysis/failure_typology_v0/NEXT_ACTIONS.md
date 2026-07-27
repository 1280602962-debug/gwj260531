# Next actions after failure typology v0

## Done
- [x] Export RTM-best poses for 7 case ligands
- [x] Pose-level hinge/clash/cognate-occupancy QC
- [x] Case reports + cross-pair table
- [x] Shortfall penalty pre-test (negative for T2)

## Do next (priority order)

### 1. Protocol write-back (docs only, no new docking)
Update SOP / master plan:
- Pair-specific `exhaustiveness_v0_1`: EGFR=8, PIK3CA/mTOR=16
- Primary **reporting** arms: `vina_mean` **and** `rtm_min_z` (not RTM alone)
- Limitations: PM48_34 8-mode; T2/T5 boundaries; clash gate negative result

### 2. Chemotype warning layer (no rank brushing)
Script to emit `warning_flags.csv` for full panel48 + panel40:
- amino-triazine-like
- morpholine-ATP / high MCS-to-cognate
- (EGFR) anilinoquinazoline / warhead flags already exist for EH40_23
Flags are **diagnostic columns**, not score modifiers in v0.

### 3. Frozen decision-rule ablation (on existing poses only)
Arms to compare:
1. vina_mean
2. rtm_min_z
3. rtm_min_z + shortfall (expect ~null)
4. consensus candidate (e.g. mean of ranks, or require both vina_min and rtm_min above threshold)

Freeze thresholds before looking at case-level wins. Success = hardneg Top10↓ without destroying Torin1/Omipalisib **or** document impossibility.

### 4. Optional 1-pager
Full 9-mode hinge/occupancy for Torin1 & Omipalisib on 4L23 → sampling vs RTM chemotype gap.

## Later (P2)
- Expand panels to 120–200
- Third pair only after T2/T5 language is in the paper outline

## Do not
- Rerun full panel / raise E for cosmetics
- Retune clash to kill PM48_26
- Claim successful C4 extrapolation
- Restart moiety cover story
