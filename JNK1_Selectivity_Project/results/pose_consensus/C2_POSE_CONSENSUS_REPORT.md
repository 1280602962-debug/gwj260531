# C2 Pose Consensus — AutoDock Vina Multi-Seed (Open Fallback)

**Scope:** geometry consensus across seeds, not Glide score validation.

## Box definitions (cognate ligand)

```json
{
  "JNK1": {
    "pdb": "3ELJ",
    "cognate": "GS7:365",
    "center": [
      22.779866666666667,
      8.341866666666666,
      31.16666666666668
    ],
    "size": [
      26.697,
      27.503999999999998,
      22.849
    ]
  },
  "JNK2": {
    "pdb": "3E7O",
    "cognate": "35F:1",
    "center": [
      -30.941892857142854,
      -32.67992857142857,
      23.934321428571426
    ],
    "size": [
      30.57,
      20.822,
      22.094
    ]
  },
  "JNK3": {
    "pdb": "3TTI",
    "cognate": "KBI:465",
    "center": [
      -21.579218750000003,
      9.343593749999998,
      -32.186875
    ],
    "size": [
      27.223999999999997,
      25.836,
      22.913
    ]
  }
}
```

## Consensus summary

|   compound_id | isoform   |   mean_vina_score |   std_vina_score |   pairwise_rmsd_mean |   fraction_pairs_rmsd_le_2A | pose_consensus_pass   |
|--------------:|:----------|------------------:|-----------------:|---------------------:|----------------------------:|:----------------------|
|           690 | JNK1      |          -8.27167 |        0.0264491 |             0.989939 |                    1        | True                  |
|           690 | JNK2      |          -8.34567 |        0.094454  |             0.622061 |                    1        | True                  |
|           690 | JNK3      |          -8.93467 |        0.221095  |             1.22072  |                    1        | True                  |
|          2231 | JNK1      |          -8.09533 |        0.0747277 |             1.02507  |                    1        | True                  |
|          2231 | JNK2      |          -7.77133 |        0.155341  |             2.31319  |                    0.333333 | False                 |
|          2231 | JNK3      |          -8.272   |        0.034535  |             0.432173 |                    1        | True                  |

## Pairwise seed RMSD

|   compound_id | isoform   |   seed_a |   seed_b |   heavy_atom_rmsd_A | consensus_le_2A   |
|--------------:|:----------|---------:|---------:|--------------------:|:------------------|
|           690 | JNK1      |        1 |        2 |           1.0063    | True              |
|           690 | JNK1      |        1 |        3 |           0.723224  | True              |
|           690 | JNK1      |        2 |        3 |           1.24029   | True              |
|           690 | JNK2      |        1 |        2 |           0.347459  | True              |
|           690 | JNK2      |        1 |        3 |           0.723342  | True              |
|           690 | JNK2      |        2 |        3 |           0.795384  | True              |
|           690 | JNK3      |        1 |        2 |           1.3876    | True              |
|           690 | JNK3      |        1 |        3 |           0.723916  | True              |
|           690 | JNK3      |        2 |        3 |           1.55065   | True              |
|          2231 | JNK1      |        1 |        2 |           0.284851  | True              |
|          2231 | JNK1      |        1 |        3 |           1.36396   | True              |
|          2231 | JNK1      |        2 |        3 |           1.42639   | True              |
|          2231 | JNK2      |        1 |        2 |           2.65265   | False             |
|          2231 | JNK2      |        1 |        3 |           2.82083   | False             |
|          2231 | JNK2      |        2 |        3 |           1.46607   | True              |
|          2231 | JNK3      |        1 |        2 |           0.0926249 | True              |
|          2231 | JNK3      |        1 |        3 |           0.579393  | True              |
|          2231 | JNK3      |        2 |        3 |           0.624501  | True              |

## Scores

|   compound_id | isoform   | pdb   |   seed |   vina_score | status   |   n_modes |
|--------------:|:----------|:------|-------:|-------------:|:---------|----------:|
|           690 | JNK1      | 3ELJ  |      1 |       -8.255 | ok       |         5 |
|           690 | JNK1      | 3ELJ  |      2 |       -8.309 | ok       |         5 |
|           690 | JNK1      | 3ELJ  |      3 |       -8.251 | ok       |         5 |
|           690 | JNK2      | 3E7O  |      1 |       -8.286 | ok       |         5 |
|           690 | JNK2      | 3E7O  |      2 |       -8.479 | ok       |         5 |
|           690 | JNK2      | 3E7O  |      3 |       -8.272 | ok       |         5 |
|           690 | JNK3      | 3TTI  |      1 |       -9.093 | ok       |         5 |
|           690 | JNK3      | 3TTI  |      2 |       -8.622 | ok       |         5 |
|           690 | JNK3      | 3TTI  |      3 |       -9.089 | ok       |         5 |
|          2231 | JNK1      | 3ELJ  |      1 |       -8.093 | ok       |         5 |
|          2231 | JNK1      | 3ELJ  |      2 |       -8.005 | ok       |         5 |
|          2231 | JNK1      | 3ELJ  |      3 |       -8.188 | ok       |         5 |
|          2231 | JNK2      | 3E7O  |      1 |       -7.659 | ok       |         5 |
|          2231 | JNK2      | 3E7O  |      2 |       -7.664 | ok       |         5 |
|          2231 | JNK2      | 3E7O  |      3 |       -7.991 | ok       |         5 |
|          2231 | JNK3      | 3TTI  |      1 |       -8.319 | ok       |         5 |
|          2231 | JNK3      | 3TTI  |      2 |       -8.26  | ok       |         5 |
|          2231 | JNK3      | 3TTI  |      3 |       -8.237 | ok       |         5 |
