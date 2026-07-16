# 非强制 OA 期刊硬筛选（按你当前实际证据）

**结论先说：**  
内容达不到 *JCIM* 近三年同类门槛，**不等于发不出非 OA**。  
以你现在的材料（选择性过滤器负校准 + 漏斗案例 + 采购解耦；酶活待测），现实可投档位是 **计算设计/分子建模中档 hybrid 刊**，不是 *JCIM* / *JMC* / *EJMC*。

**“非 OA”操作定义：** 选 **Hybrid**，投稿时勾 **Subscription / 不付 APC**。不是 MDPI 全 OA。

---

## 0. 用你现有证据做硬过滤

| 你有什么 | 期刊能不能吃 |
|----------|--------------|
| Δsel / Gly87 / ML selective 失败表 | ✅ 方法/负结果刊吃这个 |
| ML→Glide→ADMET→MD→采购漏斗 | ✅ 作案例 OK；❌ 不当“发现药物”主贡献 |
| 开源确认层（Vina/Gnina/MD） | ✅ 复现性加分；未做完前投会挨打 |
| 仅购 2 分子、尚无 IC50 | ❌ 药化应用刊不够；✅ 纯方法刊勉强够 |
| 有 JNK1/2/3 IC50（哪怕无选择性） | ✅ 可升到 *ChemMedChem* / *BMC* |
| 无 kinome / 无合成 SAR / 无细胞 | ❌ *JMC*、*EJMC*、多数高分药化主刊 |

**一句话匹配：** 你是 **“便宜 isoform 选择性过滤器不可靠”的负结果方法文**，外加一条采购校准故事；不是选择性抑制剂发现文。

---

## 1. 通过筛选：建议只认这 5 本（按投稿顺序）

### 路线 A — 不等酶活、近期就投（纯计算/负结果）

| 排序 | 期刊 | IF 量级（近年） | Hybrid 订阅？ | 与你匹配度 | 投什么、别投什么 |
|------|------|-----------------|---------------|------------|------------------|
| **1** | **Journal of Computer-Aided Molecular Design (JCAMD)** | ~3 | ✅ | **最高** | 主打 RQ-C 负结果 + 采购解耦；漏斗作案例。别写 hit discovery |
| **2** | **Journal of Molecular Graphics and Modelling (JMGM)** | ~2.5–3 | ✅ 订阅免费 | **高** | docking/MD/VS 应用方法最常见；补 redock + 无约束 MD 副本后再投更稳 |
| **3** | **Molecular Informatics** | ~3–3.5 | ✅ | **高** | 吃 ML F1=0、家族门高召回低特异；少写“成药” |

**路线 A 通过条件（缺一审稿会被卡）：**
1. 失败指标表清晰（Δsel 方向准确率、Gly87 不判别、ML F1=0）  
2. Selection（历史 Glide）与 Confirmation（开源）叙事分开  
3. 不把 690/2231 写成已验证选择性 hit  
4. 至少有：多 seed 对接共识 + 无约束 MD 协议/部分结果  

### 路线 B — 等货测完 JNK1/2/3 IC50 再投（更划算）

| 排序 | 期刊 | Hybrid？ | 匹配度 | 说明 |
|------|------|----------|--------|------|
| **1** | **ChemMedChem** | ✅ | **有活性后首选** | 计算+酶活应用文常见；µM 级、无 SI 也可诚实写“家族富集/校准失败” |
| **2** | **Bioorganic & Medicinal Chemistry (BMC)** 或 **BMCL** | ✅ | 中高 | 偏短名单+IC50；方法负结果作 Discussion，不宜当唯一卖点 |
| **3** | **JCAMD** | ✅ | 仍高 | 若活性弱/阴性，退回方法刊反而更贴切 |

**路线 B 最低湿实验门槛：**  
同批 **E1、CC-90001 + 690 + 2231** 的 JNK1/2/3 IC50；用锁定的 SI 规则报告；**禁止**用 n=2 吹 hit-rate。

---

## 2. 保底（能发，但档次明显下一档）

仅在路线 A/B 拒稿或你要“先发再说”时用：

| 期刊 | Hybrid？ | 备注 |
|------|----------|------|
| **Computational Biology and Chemistry** | ✅ | VS+MD 流程文接受度高 |
| **Journal of Biomolecular Structure and Dynamics (JBSD)** | ✅ | docking/MD 文多；声誉参差 |
| **Current Computer-Aided Drug Design** | ✅（Bentham） | 门槛低、影响力一般 |

---

## 3. 明确筛掉（别浪费时间）

| 期刊 | 原因 |
|------|------|
| *JCIM* | 近三年同类要大基准/强方法增量；你当前规模偏紧 |
| *J. Med. Chem.* / *ACS Med. Chem. Lett.* | 无/弱活性 + 无 SAR/kinome |
| *Eur. J. Med. Chem.* | 偏合成优化与强活性 |
| *Molecules* / *Pharmaceuticals* / *ACS Omega* 等全 OA 或你不想对标档 | 违背“非 OA”或档位策略 |
| *ChemistrySelect* | 偏弱应用刊 |

---

## 4. 决策树（照此选，不要犹豫）

```text
现在有没有 JNK1/2/3 IC50？
│
├─ 没有 → 只投路线 A
│         1) JCAMD
│         2) JMGM
│         3) Molecular Informatics
│         （正文禁止：selective hit / enrichment proven）
│
└─ 有了 → 看活性
          ├─ 至少一个新分子 µM 级活性（有无 SI 都行）
          │     → 1) ChemMedChem  2) BMC/BMCL  3) JCAMD
          └─ 全阴或极弱
                → 仍投 JCAMD/JMGM（主贡献=方法负结果；活性作阴性校准）
```

---

## 5. 对你课题的最终推荐（可直接执行）

1. **默认策略：等 IC50** → 首投 ***ChemMedChem***（订阅通道）。  
2. **若必须先发、不等货** → 首投 ***JCAMD***；备选 ***JMGM***。  
3. **不要**把 *JCIM* 当当前目标；证据规模上来后再谈。  
4. 投稿前最低补齐：C3 无约束 MD 副本 + 开源对接确认层（已在 wait-window 计划里）。

**一句话：** 你的内容适合 **JCAMD / JMGM / Molecular Informatics**；有三亚型酶活后升到 ***ChemMedChem***。这些都是 hybrid，可以不付 APC。
