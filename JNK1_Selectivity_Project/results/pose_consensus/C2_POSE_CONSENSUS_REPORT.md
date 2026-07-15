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
|          2157 | JNK1      |          -6.57267 |        0.0283941 |             0.27197  |                    1        | True                  |
|          2157 | JNK2      |          -6.21733 |        0.111236  |             1.72742  |                    0.333333 | False                 |
|          2157 | JNK3      |          -6.49033 |        0.102021  |             1.40784  |                    1        | True                  |

## Pairwise seed RMSD

|   compound_id | isoform   |   seed_a |   seed_b |   heavy_atom_rmsd_A | consensus_le_2A   |
|--------------:|:----------|---------:|---------:|--------------------:|:------------------|
|           690 | JNK1      |        1 |        2 |            1.0063   | True              |
|           690 | JNK1      |        1 |        3 |            0.723224 | True              |
|           690 | JNK1      |        2 |        3 |            1.24029  | True              |
|           690 | JNK2      |        1 |        2 |            0.347459 | True              |
|           690 | JNK2      |        1 |        3 |            0.723342 | True              |
|           690 | JNK2      |        2 |        3 |            0.795384 | True              |
|           690 | JNK3      |        1 |        2 |            1.3876   | True              |
|           690 | JNK3      |        1 |        3 |            0.723916 | True              |
|           690 | JNK3      |        2 |        3 |            1.55065  | True              |
|          2157 | JNK1      |        1 |        2 |            0.337481 | True              |
|          2157 | JNK1      |        1 |        3 |            0.340471 | True              |
|          2157 | JNK1      |        2 |        3 |            0.137958 | True              |
|          2157 | JNK2      |        1 |        2 |            2.32088  | False             |
|          2157 | JNK2      |        1 |        3 |            2.23454  | False             |
|          2157 | JNK2      |        2 |        3 |            0.62683  | True              |
|          2157 | JNK3      |        1 |        2 |            0.875166 | True              |
|          2157 | JNK3      |        1 |        3 |            1.5635   | True              |
|          2157 | JNK3      |        2 |        3 |            1.78486  | True              |

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
|          2157 | JNK1      | 3ELJ  |      1 |       -6.56  | ok       |         5 |
|          2157 | JNK1      | 3ELJ  |      2 |       -6.612 | ok       |         5 |
|          2157 | JNK1      | 3ELJ  |      3 |       -6.546 | ok       |         5 |
|          2157 | JNK2      | 3E7O  |      1 |       -6.361 | ok       |         5 |
|          2157 | JNK2      | 3E7O  |      2 |       -6.201 | ok       |         5 |
|          2157 | JNK2      | 3E7O  |      3 |       -6.09  | ok       |         5 |
|          2157 | JNK3      | 3TTI  |      1 |       -6.524 | ok       |         5 |
|          2157 | JNK3      | 3TTI  |      2 |       -6.595 | ok       |         5 |
|          2157 | JNK3      | 3TTI  |      3 |       -6.352 | ok       |         5 |
