# Production pose archive preparation (not published)

Local Track B production `mode_01.pdbqt` files were found outside git:

`repo_Dual_Target_Docking/Dual_Target_Docking/data/jcim_chembl_universe_v0/local_track_b_v0/poses`

This is the original saved-pose tree used by the production score table, not a new docking run.

| Item | Status |
|---|---|
| Files found | 1094 `mode_01.pdbqt` |
| Independent temp replay | 1094 match, 6 recorded skips, exit 0 |
| Git | poses are not committed |
| Zenodo | **not uploaded**; no public DOI is claimed |
| Checksum list | `data/jcim_chembl_universe_v0/local_track_b_v0/tables/production_pose_sha256_manifest_v1.csv` |

The committed replay CSV is a local check record. After this pass, the same 1094/6 counts were also obtained from the pose files themselves. Public reproduction still rebuilds statistics from score long tables and does not require these pose files.
