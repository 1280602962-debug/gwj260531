# Config

| File | Role |
|------|------|
| `campaign_c1.yaml` | **C1 科学锁**（过关线、轨道、强制回收名单）。不要传给 `run_gnina_batch.py` |
| `docking_c1.yaml` | **C1 gnina 引擎**（`num_modes: 9`，seed 42）。batch `--config` 用这个 |
| `docking_production_p2.yaml` | **冻结生产** gnina P2, 9DKB + 7ALV（`num_modes: 1`） |
| `docking_open_source.yaml` | Protocol-comparison / redock boxes (Vina + gnina) |
| `targets.yaml` | ChEMBL IDs and curation rules |
| `docking_ensemble.yaml` | UNUSED (three-state / Glide) |
| `dual_path.yaml` | UNUSED (generative Path B) |
| `model_hierarchy.yaml` | UNUSED (MiniMol / fusion) |
