# JNK1 项目汇报 PPT

## 文件

| 文件 | 说明 |
|------|------|
| `JNK1_Project_Presentation_v3.0.pptx` | 可编辑汇报幻灯片（25 主 + 1 分隔 + 4 备用） |
| `SPEAKER_NOTES.md` | 主幻灯片逐页演讲备注（与 PPTX notes 同步） |
| `SPEAKER_NOTES_BACKUP.md` | 备用幻灯片演讲备注 |
| `OUTLINE.md` | 幻灯片大纲与版式建议 |

## 重新生成

```bash
pip install python-pptx
python scripts/build_project_pptx.py
```

输出路径：`docs/presentation/JNK1_Project_Presentation_v3.0.pptx`

## 内容来源

基于 `docs/JNK1_PROJECT_REPORT.md` v3.0，叙事主线：ML → Glide → MD → 湿实验；§1.2 设计分析 + 选择性探索失败记录。
