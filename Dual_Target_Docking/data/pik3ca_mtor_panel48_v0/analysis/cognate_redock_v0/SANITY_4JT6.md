# SANITY_4JT6 — 轻量检查（升 E=16 前）

日期：2026-07-27T17:30:50+08:00  
主机：LAPTOP-3GOC1J6E  
工作根：`/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0`

**总评：S1–S4 全 PASS → 进入 Step 1（书面冻结 exhaustiveness=16）。**  
未做 Protein Prep Wizard 重跑；未改配体/面板。

---

| ID | 检查 | 结果 | 备注 |
|----|------|------|------|
| S1 | 受体是 mTOR 激酶 ATP 位点，不是 FRB/FKBP | **PASS** | `HEADER`/`TITLE`：`STRUCTURE OF MTORDELTAN-MLST8-PI-103 COMPLEX`（4JT6）。链 A，残基约 1385–2549；含激酶域区段。HETATM 仅 **X6K**（PI-103），无 rapamycin/FKBP 配体故事。 |
| S2 | 盒子来自同坐标系 X6K | **PASS** | `tables/4JT6_cocrystal_X6K.pdb` 重原子质心 ≈ (51.87, −0.00, −46.93)；`boxes/4JT6_box.json` center (51.949, 0.065, −47.707)；质心在盒子半边长内。X6K 质心到蛋白最近重原子 ≈ 2.69 Å（在口袋内）。 |
| S3 | PM48_01 = PI-103 | **PASS** | 面板：CHEMBL573339 / PI-103 / InChIKey `TUVCWJQQGGETHL-UHFFFAOYSA-N` / heavy_atoms=26。`ligands_pdbqt/PM48_01.pdbqt` 重原子 26，与 X6K 重原子数一致。 |
| S4 | 盒子尺寸合理 | **PASS** | size (20.332, 20.0, 20.0)：AABB(X6K)+5 Å，min edge≥20。 |
| S5 | E=16 诊断可复现线索 | **PASS** | 已有 `analysis/cognate_redock_v0/poses_E16/`（18 mode files）与 `tables/pm48_01_rmsd_E16_diag.csv`（4JT6 best9≈0.445）。本轮 Phase B′ 将按规范重跑并写入 `poses/cognate_E16/`（不覆盖 E=8）。 |

## 分支决策

S1–S4 全 PASS → **不回退 E=8 全面板**；进入书面 E=16 + Phase B′。
