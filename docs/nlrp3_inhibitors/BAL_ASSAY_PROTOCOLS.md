# BAL 位点抑制剂测活方法汇总

> 配套课题汇报 Slide 15；完整叙事见 [presentation/BAL_PROJECT_PRESENTATION_GUIDE.md](./presentation/BAL_PROJECT_PRESENTATION_GUIDE.md)。

---

## 核心区别：BAL vs MCC950 测活

| 特征 | BAL 系列 | MCC950 |
|------|----------|--------|
| ATP 酶 assay | **不抑制（阴性）** | 抑制 |
| 机制分型实验 | nanoDSF / SPR binning | 同左 |
| 细胞金标准 | THP-1 LPS + nigericin → IL-1β | 相同 |
| 种属 | 人/灵长类有效；鼠 WT 弱 | 人鼠均有效 |

---

## 推荐测活金字塔

见 `presentation/images/assay_pyramid.png`

| Tier | 方法 | 参考 IC50/KD |
|------|------|--------------|
| 1 结合 | SPR (NACHT) | BAL-1516 KD=14.2 nM |
| 1 结合 | nanoDSF ΔTm | BAL-0028 ~+4°C |
| 2 机制 | ATP 酶 (ADP-Glo) | BAL：**无抑制** |
| 2 机制 | SPR binning vs MCC950 | 可加合 |
| 3 细胞 | THP-1 IL-1β ELISA | BAL-0028 IC50=57.5 nM |
| 3 细胞 | ASC speck 流式 | BAL-1516 IC50=14.5 nM |
| 3 细胞 | iPSC 小胶质细胞 IL-1β | BAL-1516 IC50=11 nM |
| 4 体内 | 人源化 NLRP3 小鼠腹膜炎 | BAL-0598 口服有效 |

---

## THP-1 金标准流程（初筛推荐）

```
PMA 分化 2-3h → LPS 100 ng/mL 致敏 3h → Opti-MEM 换液
→ 化合物预孵育 30 min → nigericin 10 µM 1h → 上清 IL-1β ELISA
```

平行：LDH（毒性）、CellTiter-Glo（活力）、MCC950 阳性对照

文献：[Wilhelmsen 2025 JEM](https://doi.org/10.1084/jem.20242403) Methods

---

## 参考文献

| 主题 | 链接 |
|------|------|
| BAL-0028 完整测活 | https://doi.org/10.1084/jem.20242403 |
| BAL-1516 ASC speck / 小胶质细胞 | https://doi.org/10.1101/2025.07.01.662566 |
| DEL 发现 | https://doi.org/10.1016/j.bmcl.2024.129675 |
