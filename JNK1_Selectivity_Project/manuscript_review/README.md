# Manuscript Review Package（PaperSpine + ARS）

**Date:** 2026-07-17  
**Skills invoked from GitHub:**
- [WUBING2023/PaperSpine](https://github.com/WUBING2023/PaperSpine) → `humanize` + `logic-transfer-audit`（draft diagnosis）
- [Imbad0202/academic-research-skills-codex](https://github.com/Imbad0202/academic-research-skills-codex) → `academic-paper` Writing Quality Check + `academic-paper-reviewer` full mode

**Draft reviewed:** polished English Results & Discussion extract from  
`1_Results_and_Discussion_translated_checked.md`（PLK1–NLRP3 dual-target CADD pipeline）

> Scope caveat: only R&D was available (no Abstract / Introduction / Methods / full references). Verdicts apply to this excerpt plus submission-completeness blockers.

## Deliverables

| File | Content |
|------|---------|
| [`01_language_optimization_report.md`](01_language_optimization_report.md) | PaperSpine humanize matrix + ARS writing-quality checklist + Top-10 before/after rewrites |
| [`02_logic_and_peer_review_report.md`](02_logic_and_peer_review_report.md) | 5-persona peer review, rubrics, logic diagnosis, Must-Fix roadmap, editorial letter |

## Bottom line

| Track | Verdict |
|-------|---------|
| **Language** | Near-ready for internal draft circulation; needs **1 focused language pass** before journal submission |
| **Logic / peer review** | **Major Revision — Incomplete Manuscript**（R&D-only weighted ≈57.7） |

### Shared Must-Fixes (language + science)

1. Typo: `and and PC1`
2. Cys133 vs Cys97 hinge numbering inconsistency
3. Soften overclaims (`inhibitor`, micromolar potency, “outperforming MCC950”)
4. Fix enrichment ≈147-fold and intersection null-model framing
5. Resolve Williams-plot \(p=270\) vs 50-PC model inconsistency
6. Clarify whether ESM-2 embeddings enter the scoring function or are exploratory only
7. Reduce template cadence (16× sentence-initial `To …`; `champion` / `outclasses` / “self-confidence awareness”)

## Skill install paths (this Cloud Agent)

```text
~/.cursor/skills/paper-spine/
~/.cursor/skills/academic-research-suite/
~/.cursor/skills-cursor/{paper-spine,academic-research-suite}/
~/.codex/skills/{paper-spine,academic-research-suite}/
/workspace/.cursor/skills/{paper-spine,academic-research-suite}/   # gitignored
```
