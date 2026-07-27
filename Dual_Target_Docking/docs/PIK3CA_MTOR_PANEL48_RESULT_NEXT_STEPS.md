# PIK3CA/mTOR panel48 v0 — 结果解读与下一步

> 数据：`data/pik3ca_mtor_panel48_v0/`（E=16，seed=20260727，RTM best-of-9）  
> 对照：EGFR/HER2 panel40 v0  
> **姿态分型已完成：** `analysis/failure_typology_v0/`

---

## 1. 这一枪完成了什么

| 项 | 状态 |
|----|------|
| Cognate QC (PI-103) | **Go @ E=16** |
| 全面板 + RTM 消融 | 完成；主臂 `rtm_min_z` AUROC 0.685 |
| RTM-best pose 导出 | `analysis/rtm_best_pose_export_v1/` |
| **失败分型 v0** | **完成**（T1/T2/T5 + 跨对靶表） |

---

## 2. 分型结论（决定下一步）

| 类型 | 代表 | 含义 |
|------|------|------|
| T2 化学型同源假双靶 | PM48_26/20/21 | 弱端 pose 也干净（hinge+占用），RTM 两端都高；**clash/shortfall 压不住** |
| T1 被救 B_only | WYE-132 | Vina Top10 → RTM #40 |
| T5 误伤真 dual | Torin1 / Omipalisib | Vina #1/#3 → RTM #31/#30；4L23 RTM-best 离 hinge |
| 金标准 | PI-103 | 两端 OK；4JT6 用 mode3 |

**一句话：** 重打分必要，但第二对 **规则未闭环**；主文卖诊断+分型，不卖“尺子已外推成功”。

---

## 3. 接下来该干什么（更新后优先级）

### P0 — 已完成
- [x] 失败分型个案  
- [x] 跨对靶对照表  
- [x] shortfall 预实验（阴性：λ 怎么调，PM48_26 仍 #1）

### P1 — 立即执行（无需重对接）
1. **协议定稿回写**：E=16 仅本对；并列报告 `vina_mean` + `rtm_min_z`；写入 T2/T5 / PM48_34 8-mode Limitations  
2. **化学型警告层**（只出旗标，不改分）：与 EH40_23 警告层对齐  
3. **冻结阈值的决策消融**：vina / rtm_min_z / shortfall / consensus；证明增益是否来自“决策”而非换分  

### P2 — 分型与消融写入大纲之后
- 扩面板 120–200；第三对靶  

### 明确不做
- 全面板重跑 / E=32；拧 clash 打 PM48_26；乘客主线；宣称 C4 成功  

---

## 4. 给本地 / 下一任 agent

> 分型已落地。下一步做 **协议回写 + 警告层 + 冻结消融**；仍不要开第三对对接。细节见 `analysis/failure_typology_v0/NEXT_ACTIONS.md`。
