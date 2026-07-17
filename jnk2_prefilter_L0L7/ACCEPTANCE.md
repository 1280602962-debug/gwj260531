# Prefilter acceptance (docking-front)

- [x] L5 calibration + thresholds.json exists
- [x] core vs full Sim Jaccard = 0.2676 (<1)
- [ ] 56d LOO not systematically discarded — see L5_calibration.md / L5/loo.json
- [ ] Novel random 50 manual review — pending human

## Delivered dock_ready: 10000
{
  "sim_yl": 4200,
  "sim_56d": 1800,
  "novel": 3500,
  "pan": 500
}

Track ID intersections (must be 0): {'sim_yl∩sim_56d': 0, 'sim_yl∩novel': 0, 'sim_yl∩pan': 0, 'sim_56d∩novel': 0, 'sim_56d∩pan': 0, 'novel∩pan': 0}
