# Dual-dock merge / Pareto (P2 only)

This folder no longer stores docking scores.

Historical **Glide XP** dual-target merges and 6-molecule shortlists were deleted. Do not treat any leftover notes as production ranks.

Generate P2 tables locally:

```bash
JOBS=8 bash scripts/run_funnel_p2.sh
```

Outputs: `results/repurposing/pareto_merged_scores.csv` and `pareto_shortlist.csv`.
