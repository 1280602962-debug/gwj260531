# PM48 RDKit reprep — prep sensitivity

- Protocol: 4L23/4JT6, E=16, seed=20260727, n_modes=9, RTM best-of-9
- Prep: RDKit ETKDG + meeko

## Directional by prep

```
       prep       arm  auroc_D_vs_A  auroc_D_vs_B  summary_min
ligprep_old vina_mean      0.698413      0.597222     0.597222
ligprep_old rtm_min_z      0.611111      0.791667     0.611111
  rdkit_new vina_mean      0.722222      0.671296     0.671296
  rdkit_new rtm_min_z      0.519841      0.671296     0.519841
```

- |Δvina_mean| median=0.272
- |Δrtm_min_z| median=0.342
