# L5 calibration

- Feature rows: **239,990**
- Spearman(max_tc_core, max_tc_full): **ρ=0.662** (p=0.00e+00)
- Jaccard(Sim by core-window vs full-window): **0.268** (<1 ⇒ warhead-strip matters)
- ErG percentiles: p50=0.691, p75=0.734, p90=0.767

## Working thresholds (written to config/thresholds.json)
```json
{
  "sim_tc_core_lo": 0.22,
  "sim_tc_core_hi": 0.55,
  "near_duplicate_tc_full": 0.7,
  "novel_tc_core_max": 0.22,
  "novel_tc_core_pref_lo": 0.08,
  "novel_erg_min": 0.7344,
  "novel_erg_rule": "erg_max >= p75 OR hinge_hits>=1",
  "erg_p50": 0.6912,
  "erg_p75": 0.7344,
  "erg_p90": 0.767,
  "watch_quota_frac": 0.2,
  "quotas": {
    "sim_yl": 4200,
    "sim_56d": 1800,
    "novel": 3500,
    "pan": 500
  },
  "spearman_core_vs_full": {
    "rho": 0.6625,
    "p": 0.0
  },
  "jaccard_sim_core_vs_full_window": 0.2676,
  "n_features": 239990,
  "chosen_reason": "ARS initial window 0.22\u20130.55; ErG Novel gate at p75; confirmed by grid stability of Sim size"
}
```

## Leave-one-anchor-out (nearest_core preference recovery into Sim via other anchors)

- **YL5084**: {'n_pref': 31285, 'recover_sim_frac_at_0.22_0.55': 0.0101}
- **YL2056**: {'n_pref': 70045, 'recover_sim_frac_at_0.22_0.55': 0.0217}
- **JNK-IN-8**: {'n_pref': 81292, 'recover_sim_frac_at_0.22_0.55': 0.0233}
- **56d**: {'n_pref': 57368, 'recover_sim_frac_at_0.22_0.55': 0.0211}

## Tc grid (Sim counts)

- lo=0.15, hi=0.5: sim=112,843, near_dup=0, novel_lt_lo=127,092
- lo=0.15, hi=0.55: sim=112,868, near_dup=0, novel_lt_lo=127,092
- lo=0.15, hi=0.6: sim=112,882, near_dup=0, novel_lt_lo=127,092
- lo=0.15, hi=0.7: sim=112,898, near_dup=0, novel_lt_lo=127,092
- lo=0.18, hi=0.5: sim=53,917, near_dup=0, novel_lt_lo=186,018
- lo=0.18, hi=0.55: sim=53,942, near_dup=0, novel_lt_lo=186,018
- lo=0.18, hi=0.6: sim=53,956, near_dup=0, novel_lt_lo=186,018
- lo=0.18, hi=0.7: sim=53,972, near_dup=0, novel_lt_lo=186,018
- lo=0.2, hi=0.5: sim=32,132, near_dup=0, novel_lt_lo=207,803
- lo=0.2, hi=0.55: sim=32,157, near_dup=0, novel_lt_lo=207,803
- lo=0.2, hi=0.6: sim=32,171, near_dup=0, novel_lt_lo=207,803
- lo=0.2, hi=0.7: sim=32,187, near_dup=0, novel_lt_lo=207,803
- lo=0.22, hi=0.5: sim=17,242, near_dup=0, novel_lt_lo=222,693
- lo=0.22, hi=0.55: sim=17,267, near_dup=0, novel_lt_lo=222,693
- lo=0.22, hi=0.6: sim=17,281, near_dup=0, novel_lt_lo=222,693
- lo=0.22, hi=0.7: sim=17,297, near_dup=0, novel_lt_lo=222,693
- lo=0.25, hi=0.5: sim=7,497, near_dup=0, novel_lt_lo=232,438
- lo=0.25, hi=0.55: sim=7,522, near_dup=0, novel_lt_lo=232,438
- lo=0.25, hi=0.6: sim=7,536, near_dup=0, novel_lt_lo=232,438
- lo=0.25, hi=0.7: sim=7,552, near_dup=0, novel_lt_lo=232,438

## Anchor self-check (full molecules as if products)

- YL5084 full-Tc to anchors: YL5084=1.000, YL2056=0.805, JNK-IN-8=0.377, 56d=0.151
- YL2056 full-Tc to anchors: YL5084=0.805, YL2056=1.000, JNK-IN-8=0.381, 56d=0.153
- JNK-IN-8 full-Tc to anchors: YL5084=0.377, YL2056=0.381, JNK-IN-8=1.000, 56d=0.289
- 56d full-Tc to anchors: YL5084=0.151, YL2056=0.153, JNK-IN-8=0.289, 56d=1.000