# 汇报图片资源清单与引用规范

> 本目录 [`images/`](./images/) 包含可直接插入 PPT 的图片，以及需从原始文献获取的论文原图索引。  
> **使用前请遵守各来源的版权与引用要求。**

---

## 一、本地图片（可直接用于 PPT）

| 文件名 | 内容 | 建议用于 Slide | 来源/许可 |
|--------|------|----------------|-----------|
| `binding_sites_schematic.png` | NLRP3 双位点示意图（MCC950 vs BAL） | Slide 3 | **自制**，可自由使用 |
| `7PZC_assembly.png` | NLRP3 decamer 冷冻电镜结构（含 MCC950） | Slide 6, 7 | [RCSB PDB 7PZC](https://www.rcsb.org/structure/7PZC)，引用见下 |
| `7ALV_assembly.png` | NLRP3 NACHT 晶体结构 | Slide 7 | [RCSB PDB 7ALV](https://www.rcsb.org/structure/7ALV) |
| `7VTP_assembly.png` | 人源 NLRP3 hexamer | 附录 | [RCSB PDB 7VTP](https://www.rcsb.org/structure/7VTP) |
| `8ETR_assembly.png` | NLRP3 NACHT + G2394 | 附录 | [RCSB PDB 8ETR](https://www.rcsb.org/structure/8ETR) |
| `murcko_top1.png` | 主 Murcko 骨架结构式 | Slide 9 | **RDKit 绘制**，基于专利数据 |
| `murcko_top2.png` | 第二 Murcko 骨架结构式 | Slide 9 | **RDKit 绘制** |
| `patent_data_summary_chart.png` | 专利数据统计柱状图 | Slide 8 | **自制** |
| `assay_pyramid.png` | 测活金字塔四层体系 | Slide 15 | **自制** |
| `four_phase_workflow.png` | 四阶段实施路线图 | Slide 12, 16 | **自制** |

### RCSB PDB 图片引用格式（贴于 PPT 图注）

```
Image from RCSB PDB (https://www.rcsb.org/structure/XXXX)
Hochheiser IV et al. (2022) Nature 610:374-379 [7PZC]
```

| PDB | 引用 DOI |
|-----|----------|
| 7PZC | https://doi.org/10.1038/s41586-022-04467-w |
| 7ALV | https://doi.org/10.1016/j.jmb.2021.167309 |
| 7VTP | https://doi.org/10.1073/pnas.2121353119 |
| 8ETR | https://doi.org/10.1021/acs.jmedchem.2c01250 |

---

## 二、需从论文获取的原图（高价值，建议截图）

> 以下图片**未包含在仓库中**（版权归属出版社/作者）。请从链接下载 preprint/PDF，按 Figure 编号截取，并在 PPT 注明出处。学术汇报通常属于合理使用，正式发表需获得许可。

### BAL 系列核心（优先）

| Figure | 内容 | 用于 Slide | 获取链接 |
|--------|------|------------|----------|
| Hartman 2024 **Fig. 1** | BAL-0028 结构式 + DEL 筛选流程 | Slide 6 | https://doi.org/10.1016/j.bmcl.2024.129675 |
| Wilhelmsen 2025 **Fig. 1B** | THP-1 IL-1β 剂量响应（BAL vs MCC950） | Slide 6, 15 | https://doi.org/10.1084/jem.20242403 |
| Wilhelmsen 2025 **Fig. 2B** | ASC speck 成像抑制 | Slide 15 | 同上 |
| Wilhelmsen 2025 **Fig. 4A** | ATP 酶 assay（BAL 不抑制） | Slide 11, 15 | 同上 |
| Wilhelmsen 2025 **Fig. 4C-D** | nanoDSF 热稳定化 + MCC950 竞争 | Slide 11 | 同上 |
| Torp 2025 **Fig. 1c-e** | SPR KD + binning（可加合） | Slide 6, 11 | https://doi.org/10.1101/2025.07.01.662566 |
| Torp 2025 **Fig. 2a-d** | BAL-1516 cryo-EM 结合模式（**最重要**） | Slide 6, 14 | 同上 |
| Torp 2025 **Fig. 4a** | ASC speck 流式 IC50 | Slide 15 | 同上 |
| Torp 2025 **Fig. 4d** | 小胶质细胞 IL-1β 剂量响应 | Slide 15 | 同上 |

### 结构生物学背景

| Figure | 内容 | 用于 Slide | 获取链接 |
|--------|------|------------|----------|
| Hochheiser 2022 **Fig. 1** | NLRP3 decamer 整体架构 | Slide 7 | https://doi.org/10.1038/s41586-022-04467-w |
| Ohto 2022 **Fig. 2** | NLRP3 寡聚化调控 | 附录 | https://doi.org/10.1073/pnas.2121353119 |

### AI 方法学

| Figure | 内容 | 用于 Slide | 获取链接 |
|--------|------|------------|----------|
| FoldBench 2025 **Fig. 2-3** | AF3/Boltz/Chai 蛋白-配体成功率对比 | Slide 13 | https://doi.org/10.1038/s41467-025-67127-3 |
| Allosteric Paradox 2026 **Fig. 1** | 变构位点预测短板 | Slide 13 | https://doi.org/10.64898/2026.02.24.707829 |

### NLRP3 背景

| Figure | 内容 | 用于 Slide | 获取链接 |
|--------|------|------------|----------|
| Broz & Dixit 2016 **Fig. 1** | 炎性小体组装机制 | Slide 2 | https://doi.org/10.1038/nri.2016.58 |
| Coll 2019 **Fig. 1** | MCC950 结合机制 | Slide 3 | https://doi.org/10.1038/s41589-019-0277-7 |

---

## 三、9IHN / 9Q8V（BAL-1516 共晶）特别说明

- PDB 状态：**HPUB**（尚未公开下载 coordinates）
- RCSB CDN **暂无**预览图（已尝试，404）
- **替代方案**：
  1. 使用 Torp 2025 bioRxiv **Fig. 2** 的 cryo-EM 密度图（推荐）
  2. 使用本地 `7PZC_assembly.png` 作为 NLRP3 背景
  3. 9IHN 发布后从 RCSB 下载并渲染：
     - https://www.rcsb.org/structure/9IHN
     - https://www.rcsb.org/structure/9Q8V

---

## 四、PPT 图注模板

### 自制图
```
Source: Generated for this project (2026)
```

### RCSB PDB
```
Image: RCSB PDB 7PZC — Hochheiser et al. Nature 2022
https://www.rcsb.org/structure/7PZC
```

### 论文原图
```
Adapted from Torp et al. bioRxiv 2025, Fig. 2
https://doi.org/10.1101/2025.07.01.662566
```

---

## 五、在线渲染工具（可自行生成更多结构图）

| 工具 | 链接 | 用途 |
|------|------|------|
| RCSB Mol* | https://molstar.org/viewer/ | 上传 PDB 渲染结合位点 |
| RCSB Images | https://images.rcsb.org | 在线生成结构图片 |
| PyMOL / ChimeraX | 本地软件 | 高质量出版级渲染 |

### 推荐渲染场景（9IHN 发布后）
1. NLRP3 NACHT + BAL-1516 结合沟槽表面图
2. Y258/H260 氢键近距离图
3. MCC950 口袋空腔 vs BAL 结合对比图

---

## 六、图片文件大小参考

```
7PZC_assembly.png          ~85 KB
7ALV_assembly.png          ~67 KB
binding_sites_schematic.png ~52 KB
four_phase_workflow.png    ~46 KB
assay_pyramid.png          ~45 KB
patent_data_summary_chart.png ~39 KB
murcko_top1/2.png          ~16 KB each
```

所有图片均为 PNG 格式，适合直接插入 PowerPoint / Keynote / Google Slides。

---

*最后更新：2026-07-10*
