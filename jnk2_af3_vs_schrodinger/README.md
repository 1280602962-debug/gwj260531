# JNK2 Phase 0: AF3 covalent vs Schrödinger covalent docking

See **[AF3_vs_薛定谔共价对接对比总结.md](./AF3_vs_薛定谔共价对接对比总结.md)** for the Chinese summary.

**Pipeline (post Phase 0):**
- [SCREENING_SOP.md](./SCREENING_SOP.md) — frozen ranking rules (AF3 mPAE final; Glide coarse filter)
- [LIBRARY_FUNNEL_COMPLETE.md](./LIBRARY_FUNNEL_COMPLETE.md) — full library funnel (ChEMBL amine→acrylamide, dual-track, Glide→AF3)
- [NEXT_STEPS.md](./NEXT_STEPS.md) — action list

**Data:**
- `af3_analysis/` — AF3 confidence metrics & discrimination tables  
- `schrodinger_analysis/` — Glide covalent docking scores & AF3 correlation tables  

Raw AF3/Glide pose archives are kept locally under `D:\CADD paper exercise\JNK2\` (not vendored here).
