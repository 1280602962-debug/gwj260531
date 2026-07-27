# 4JT6 PM48_34 mode count note

- Requested: num_modes=9, seed=20260727, exhaustiveness=16
- Vina stdout lists 9 rows, but mode 9 affinity = **+75.84** (physically nonsense)
- Output PDBQT contains only **8 MODEL** blocks (modes 1–8, affinities −9.85 … −8.58)
- Rerun with identical config is **bit-identical** → not a random crash
- Treatment: accept 8 valid poses; do not invent a 9th pose; flag in job_status/scores
