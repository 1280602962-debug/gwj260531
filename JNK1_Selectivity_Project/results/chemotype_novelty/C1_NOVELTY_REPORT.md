# C1 Chemotype Novelty Audit

- Fingerprint: Morgan/ECFP4 r=2, 2048-bit
- ChEMBL/JNK pool size: **1835** unique SMILES

## Summary

|   compound_id | smiles                                     | murcko                            |   maxTc_vs_literature_refs | nearest_literature_ref   |   maxTc_vs_chembl_jnk_pool | nearest_chembl_id   | interpretation                                  |
|--------------:|:-------------------------------------------|:----------------------------------|---------------------------:|:-------------------------|---------------------------:|:--------------------|:------------------------------------------------|
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O |                     0.2329 | Q63                      |                     0.2716 | CHEMBL101035        | chemically distant from curated JNK set (ECFP4) |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       |                     0.2254 | Q63                      |                     0.2676 | CHEMBL1761572       | chemically distant from curated JNK set (ECFP4) |

## vs literature references (full matrix)

|   compound_id | hit_smiles                                 | query_murcko                      | ref_name   | ref_murcko                                      |   ecfp4_tanimoto | same_murcko   |
|--------------:|:-------------------------------------------|:----------------------------------|:-----------|:------------------------------------------------|-----------------:|:--------------|
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | SP600125   | O=C1c2ccccc2-c2n[nH]c3cccc1c23                  |           0.1077 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | CC-90001   | c1cc(NC2CCCCC2)ncn1                             |           0.1974 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | CC-930     | c1ccc(Nc2nc3cnc(NC4CCCCC4)nc3n2[C@H]2CCOC2)cc1  |           0.1236 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | E1         | c1ccc(Nc2ccnc(Nc3ccccc3)n2)cc1                  |           0.2    | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | Q63        | c1ccc(Nc2ccncn2)cc1                             |           0.2329 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | TCS_JNK_6O | O=C(Cc1ccccc1)Nc1ccccn1                         |           0.1529 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | AS602801   | c1ccc2sc(Cc3ccnc(OCc4ccc(CN5CCOCC5)cc4)n3)nc2c1 |           0.0729 | False         |
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | JNK-IN-8   | O=C(Nc1ccc(Nc2nccc(-c3cccnc3)n2)cc1)c1ccccc1    |           0.1809 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | SP600125   | O=C1c2ccccc2-c2n[nH]c3cccc1c23                  |           0.0952 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | CC-90001   | c1cc(NC2CCCCC2)ncn1                             |           0.1429 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | CC-930     | c1ccc(Nc2nc3cnc(NC4CCCCC4)nc3n2[C@H]2CCOC2)cc1  |           0.1412 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | E1         | c1ccc(Nc2ccnc(Nc3ccccc3)n2)cc1                  |           0.1512 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | Q63        | c1ccc(Nc2ccncn2)cc1                             |           0.2254 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | TCS_JNK_6O | O=C(Cc1ccccc1)Nc1ccccn1                         |           0.131  | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | AS602801   | c1ccc2sc(Cc3ccnc(OCc4ccc(CN5CCOCC5)cc4)n3)nc2c1 |           0.1236 | False         |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | JNK-IN-8   | O=C(Nc1ccc(Nc2nccc(-c3cccnc3)n2)cc1)c1ccccc1    |           0.1613 | False         |

## Nearest ChEMBL JNK neighbor

|   compound_id | hit_smiles                                 | query_murcko                      | nearest_chembl_id   | nearest_smiles                                            | nearest_murcko                               | nearest_source_isoform   |   max_ecfp4_tanimoto_chembl_jnk | same_murcko_as_nearest   |
|--------------:|:-------------------------------------------|:----------------------------------|:--------------------|:----------------------------------------------------------|:---------------------------------------------|:-------------------------|--------------------------------:|:-------------------------|
|           690 | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 | O=C(Nc1ccccc1)Nc1ncc2c(n1)CCCC2=O | CHEMBL101035        | Cc1ccc(-n2nc(C(C)(C)C)cc2NC(=O)Nc2ccc(OCCN3CCOCC3)cc2)cc1 | O=C(Nc1ccc(OCCN2CCOCC2)cc1)Nc1ccnn1-c1ccccc1 | JNK2                     |                          0.2716 | False                    |
|          2157 | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              | c1cnc(NCc2cccc3c2OCCCO3)nc1       | CHEMBL1761572       | Cc1ccc(-c2cccc3cnc(Nc4ccc5c(c4)OCO5)nc23)cn1              | c1cncc(-c2cccc3cnc(Nc4ccc5c(c4)OCO5)nc23)c1  | JNK3                     |                          0.2676 | False                    |

## Interpretation guide

- Tc ≥ 0.55: treat as known-like / near-neighbor risk for novelty claim.
- 0.35–0.55: moderate; discuss scaffold relationship carefully.
- < 0.35: ECFP4-distant from curated set; still may share hinge-binder pharmacophore.
- Same Murcko as E1/CC-90001/SP600125 would be a red flag even at modest Tc.
