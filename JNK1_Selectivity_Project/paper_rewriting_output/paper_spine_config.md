# PaperSpine 配置摘要

| 字段 | 值 |
|------|-----|
| workflow | `build_from_materials` |
| scene | `journal` |
| tier | `flash` |
| output_language | `en`（另出中文翻译包） |
| target_name | JCIM / Molecules / ChemMedChem（贡献确认后再定刊） |
| materials_dir | `JNK1_Selectivity_Project/` |
| draft_path | `docs/JNK1_PROJECT_REPORT.md`（素材主报告，非定稿） |
| reference_mode | `local_first` |
| citation_target_count | 20 |
| humanize_tier | `medium` |
| ui_language | `zh` |

## 用户动机（待 Research 后确认）

诚实的端到端 ML–对接–MD–采购管线，用于从商业库富集 **JNK 家族结合剂**；并对常见计算 isoform 选择性预测方法做关键 benchmark（Δsel / Gly87 / ML 选择性），证明在近同源 ATP 口袋上失效。湿实验（690+2157 + E1/CC-90001）验证**家族活性富集**，而非宣称选择性发现。

## 特殊约束

- 不得在无 IC50 证据时宣称「发现 JNK1 选择性抑制剂」
- 必须把「计算选择性预测失败」推到前台作为方法学贡献之一
- n=2 新分子 + 已有阳性；kinome 选择性未测
- isoform 选择性 ≠ kinome 选择性，写作必须分开
