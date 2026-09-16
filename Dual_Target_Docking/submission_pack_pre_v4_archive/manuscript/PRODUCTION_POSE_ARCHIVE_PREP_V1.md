# Production pose archive

Track B production poses are the original saved-pose tree used by the score table, not a new docking run.

Path: `data/jcim_chembl_universe_v0/local_track_b_v0/poses`

| Item | Status |
|---|---|
| Files | 9742 `*.pdbqt` (1094 `mode_01.pdbqt`) |
| Independent replay | 1094 match, 6 recorded skips, exit 0 |
| Git | committed on `cursor/jcim-language-polish-0b1a` |
| Zenodo | **not uploaded**; no public DOI is claimed |
| Checksum list | `data/jcim_chembl_universe_v0/local_track_b_v0/tables/production_pose_sha256_manifest_v1.csv` |

The committed replay CSV remains a local check record. Public statistics still rebuild from the score long tables; the pose tree now also lets a clone reread `REMARK VINA RESULT` for the 1094 scored jobs.
