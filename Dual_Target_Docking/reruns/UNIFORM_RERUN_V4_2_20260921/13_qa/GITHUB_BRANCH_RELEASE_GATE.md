# GITHUB_BRANCH_RELEASE_GATE

GITHUB_BRANCH_RELEASE_GATE = PASS

reason: remote tree verified after ordinary push

local_orphan_commit: af328e3e52a6bdd4e1c49575c1f6f515e227ce10
parent_count: 0
ncommits_at_first_push: 1
worktree: /tmp/v4_2_final_audited_release_20260930
branch: cursor/v4-2-final-audited-results-20260930
remote: origin https://github.com/1280602962-debug/gwj260531.git

remote_push: SUCCESS
push_force: NO
ls-remote_sha: af328e3e52a6bdd4e1c49575c1f6f515e227ce10
sha_match_local_head: YES
github_branch_sha: af328e3e52a6bdd4e1c49575c1f6f515e227ce10
github_commit_parents: 0
pack_written: 161367 objects, 129.81 MiB

remote_path_ok:
- README.md
- .gitignore
- Dual_Target_Docking/results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv
- Dual_Target_Docking/results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv
- Dual_Target_Docking/scripts/rerun_v4_2/formal_metrics_lib.py
- Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRE_COMMIT_LOCAL_RELEASE_GATE.md

remote_ligands: sdf=785 pdbqt=785
remote_absent_404:
- JNK1_Selectivity_Project
- Dual_Target_Docking/current_score_master.csv
- Dual_Target_Docking/scripts/analysis

source_worktree: /tmp/pr39_fiveseed
source_branch: cursor/scientific-freeze-c7cc
source_02_ligands_still_symlink: YES

This Gate file is a second ordinary commit. The first orphan commit was not amended.

LOCAL_RELEASE_COMPLETE
REMOTE_RELEASE_COMPLETE
