# FUTURE_APPLICATION_SELECTION_RULES

Status: principles only. No pair is selected in Stage 1.

## Role

If a new unlabeled library is ever screened, that exercise is a **benchmark-informed application**.

It is not:

- blind prospective validation;
- independent external validation;
- experimental discovery of a dual inhibitor.

## What this phase forbids

- naming a shortlist of pairs;
- scoring or ranking the eight formal pairs;
- choosing a winner;
- docking a new library;
- writing a candidate list;
- running MD.

Knowing the existing benchmark results does not authorize a pair choice in this document.

## Selection principles for a later freeze

A later, separate freeze may choose one pair, and at most one failure-control pair, using a written checklist. Allowed criteria, none of which may be silently reduced to “highest AUROC”:

1. Dual-target biological rationale.
2. Formal sample supply (complete both-pocket scores; neither supply is not required for the application itself).
3. Whether both directional arms are interpretable, including known anti-discrimination or method-worsening cases if the application is explicitly a negative control.
4. Consistency or interpretability of M0–M3 on that pair, used transparently rather than as a hidden AUROC sort.
5. Receptor-robustness information, if the pair is in the frozen alternative-receptor subpanel.
6. Holdout readiness, if a reserved-sample check is desired.
7. Suitability for computational candidate ranking without wet experiments.

If several pairs pass the later checklist, the tie-break must be written in that later freeze and must not be AUROC.

## Evidence boundary for any future library ranking

Without wet experiments, outputs may be called computationally prioritized dual-target candidates.

Forbidden:

- “we discovered new dual-target inhibitors”;
- treating MD stability as activity validation;
- inventing a new consensus score for the library after seeing the list;
- changing the joint-ranking formula after seeing candidates.

The ranking rule, if used, must be the already frozen minimax joint rank on M0, with seed/receptor/method-consistency checks used only as robustness annotations.
