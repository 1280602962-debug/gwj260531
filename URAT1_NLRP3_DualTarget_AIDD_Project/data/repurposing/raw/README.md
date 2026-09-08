# ChEMBL 原始导出

将以下三个文件放在本目录：

- `Phase1_2_3.xls`（或 `.xlsx`）
- `Level1 ATC.xls`
- `Level 2 ATC.xls`

然后运行：

```bash
python3 scripts/build_repurposing_library.py --input-dir data/repurposing/raw
```

输出见 `data/repurposing/repurposing_manifest.csv`。
