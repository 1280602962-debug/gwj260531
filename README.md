# gwj260531

**完整项目均在 `JNK1_Selectivity_Project/` 目录内。**

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py --library data/libraries/your_library.csv --output results/screening_v2
```

详见 [JNK1_Selectivity_Project/README.md](JNK1_Selectivity_Project/README.md)。
