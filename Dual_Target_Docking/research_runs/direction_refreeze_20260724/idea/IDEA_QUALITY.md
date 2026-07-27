# Idea quality: three direction cards

**Skill:** ResearchStudio `idea_quality` (absolute track)  
**Date:** 2026-07-24

## Idea A — Dual-VSDS-Decision (adopted)

### Title
Architecture-agnostic dual-target docking decision ruler

### Motivation
Naive dual-end docking fusion is near-random on Dual vs A-only/B-only (AUROC≈0.55). Dual ligands span merged/linked/other architectures, so a mechanism-specific passenger story cannot cover the task.

### Method
1. Fixed dual-end docking + top-K poses  
2. Optional ML rescoring (RTMScore) + pose/interaction gates  
3. Per-target calibration  
4. Weak-end min/shortfall decision  
5. Evaluate on leakage-controlled four-class panels; stratify by architecture  

### Scores
- **A Problem = 4** — four-class hard-negative gap is real and under-served  
- **B Method = 3** — components known; novelty is disciplined protocol + benchmark  
- **C Fit = 5** — directly solves “cover all dual architectures”  
- **Overall ≈ 67 (strong-borderline); A/C gate passes**

---

## Idea B — Passenger/moiety (rejected)

### Motivation
Second pharmacophore contaminates single-pocket scores.

### Method
Annotate moieties; mask passenger; re-fuse.

### Scores
- **A = 2** — only covers bipartite/linked-like ligands  
- **B = 3**  
- **C = 1** — fails user constraint to cover all duals  
- **Verdict: weak/borderline (C gate)**

---

## Idea C — Benchmark-only (fallback)

### Motivation
Field lacks Dual/A_only/B_only/neither docking panels with architecture labels.

### Method
Curate multi-pair public resource; report failure of naive fusion without claiming a new decision rule.

### Scores
- **A = 4, B = 2, C = 3 → ~50 borderline**  
Use if Phase-1 protocol ablations fail Kill criteria.
