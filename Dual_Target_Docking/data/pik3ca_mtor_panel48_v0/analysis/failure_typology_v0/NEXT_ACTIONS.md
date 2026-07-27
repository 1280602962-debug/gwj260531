# Next actions after failure typology v0 / P1

## Done
- [x] Export RTM-best poses for 7 case ligands
- [x] Pose-level hinge/clash/cognate-occupancy QC
- [x] Case reports + cross-pair table
- [x] Shortfall penalty pre-test (negative for T2)
- [x] Protocol write-back (E=8 vs E=16; dual reporting arms; Limitations)
- [x] Chemotype `warning_flags.csv` (panel48 + panel40; flags do not enter score)
- [x] Frozen decision ablation v0 → **cannot jointly satisfy** hardneg↓ + T5 spare

## Do next (optional / P2)
1. Optional 1-pager: full 9-mode hinge/occupancy for Torin1 & Omipalisib on 4L23
2. Expand panels to 120–200 (design only until user green-lights docking)
3. Third pair only after T2/T5 language is in the paper outline

## Do not
- Rerun full panel / raise E for cosmetics
- Retune clash to kill PM48_26
- Claim successful C4 extrapolation
- Restart moiety cover story
- Feed warning flags into gated scores in v0
