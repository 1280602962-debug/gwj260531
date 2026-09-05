# C5 短名单已冻结（零新对接）

冻结日期：2026-09-05。  
脚本：`scripts/freeze_c5_shortlist.py`。  
机器读出：`shortlist_freeze_summary.json`。

这些名字是 **putative dual-node 候选假说**，不是 dual inhibitors。对接分不是亲和力。没有湿实验。

**筛选过程逐步审计 + 单槽 MD 分子预登记：`docs/C5_SCREENING_AUDIT_AND_MD_PICK.md`。**  
读表前必看两条：A1 门 LR+ 只有 2.20，且**挡掉了同池中的 lesinurad / verinurad / puliginurad**；A2 门 LR+ = 1.00（无判别力），所以备份 21 个在 URAT1 侧基本等于未筛，**不得称"次优候选"**。

---

## 冻结规则（不再改名单）

1. **成员资格 = 已算好的 W3 机械交**。不因为 W2 IFP 把谁升上去或踢下来。
2. **主名单** = tier1 去掉 GSK-3008348（结构对照，不是候选）。
3. **备份** = tier2 去掉 3 个头孢 β-内酰胺（Cefetrizole / Cefazedone / Cefoxazole）。
4. W2 IFP 只是注释。A1 仍是 URAT1 判别门；IFP 是更严的确认子集。
5. lesinurad / verinurad / puliginurad 在备份里，是**已知羧酸药语境**，不是这次筛出来的新 hit。

缺 SDF：0。临床 9DKB 三种子姿态都在。

---

## 主名单（n=12，可写进正文候选表）

| 名字 | ligand_id | 最高临床阶段 | A2 结构门种子数 | W2 IFP 过门种子 |
|---|---|---:|---:|---:|
| Lanifibranor | REP_00940 | 3 | 3 | 3 |
| Tonapofylline | REP_00950 | 3 | 2 | 1 |
| Caficrestat | REP_00992 | 3 | 3 | 3 |
| Admilparant | REP_01147 | 3 | 2 | 2 |
| Lintitript | REP_05141 | 2 | 3 | 0 |
| Spiroglumide | REP_05302 | 2 | 3 | 1 |
| Cavosonstat | REP_06119 | 2 | 3 | 1 |
| Runcaciguat | REP_06643 | 2 | 3 | 2 |
| Fulimetibant | REP_06875 | 2 | 3 | 3 |
| PF-03882845 | REP_07580 | 1 | 3 | 3 |
| PSI-697 | REP_07704 | 1 | 3 | 3 |
| CR-3465 Free Acid | REP_08210 | 1 | 3 | 3 |

其中 8 个 W2 IFP ≥2/3（Lanifibranor, Caficrestat, Admilparant, Runcaciguat, Fulimetibant, PF-03882845, PSI-697, CR-3465）。另外 4 个仍在主名单：它们过了 A1，只是 IFP 更严。**不要因此降级。**

---

## 不进候选表

| 名字 | 原因 |
|---|---|
| GSK-3008348 | 结构对照，不是候选 |
| Cefetrizole / Cefazedone / Cefoxazole | β-内酰胺；不得进可报告备份或 MD |

---

## 备份（n=21，可报告；不是第二套主名单）

Terutroban, Glycocholic acid, Aleglitazar, Tenivastatin, Abitesartan, Aseripide, Florantyrone, GW590735, Deferitrin, Posenacaftor, Pocenbrodib, TAS-119, GSK-2018682, ORE-1001, AGG-523, AZD-8075, Cholylsarcosine, ASP-7657，以及已知 URAT1 酸 **lesinurad / verinurad / puliginurad**。

备份里只有 Posenacaftor、Pocenbrodib 的 W2 IFP ≥2/3。其余多数是 A2 过、A1 不过（CNNscore Top-1 酸根远离 Arg477）。这正是 A2 无判别力的体现，所以它们不是主名单。

---

## 本地还要不要算

**选出候选之前：不用。** 必须的 gnina 已经齐了。

| 项 | 还要本机 gnina？ |
|---|---|
| 点这 12 个名字 | 否 |
| 重试 15 个空 SDF 诱饵 / W1 其余 29 格 | 可选 SI，不当解锁 |
| 临床 156 重对接 / Rank 轨 | 禁止 |
| W5 MD | 冻结已完成，但仍要先改正文并显式授权；`md_authorized` 仍为 false。只许 2–3 个主名单分子 + 晶体对照 |

下一步是写正文，不是再对接。
