# Figure visual QA

Date: 2026-09-20

Scientific values were not changed. Only layout: size, spacing, annotation offset, legend room, margins.

## Figure 2 = PASS

File: `figures/jcim_article/Fig2_negative_class_formulation.pdf` and `.png`

Fixes:

- Height 8.50 → 9.00 in (ACS depth limit 9.167 in).
- `hspace` 0.86 → 0.98 so panel legends do not sit on the next axis label.
- Panel C `n=` annotation moved from `point + 0.025` to `CI_hi + 0.028`, so the PIK3CA neither n=4 label is outside the interval.
- x-limit on panel C extended to 1.12 to keep that annotation inside the axes.
- Right/top/bottom margins adjusted; k= labels on panel D remain inside the left margin.

Checked at publication size: axis labels, panel letters, legends, CIs, n=4, long pair names, and tick labels. PDF export uses the existing `savefig.bbox=standard` path (no clipping flag).

## Figure S4 = PASS

File: `figures/jcim_article/FigS4_label_source_robustness.pdf` and `.png`

Fixes:

- Size 7 × 3.80 → 7 × 4.35 in.
- Larger left margin and wspace so JAK1/TYK2 and document labels are not cut.
- Colorbar pad increased; dagger moved from y=8.55 (over the heatmap/colorbar) to an annotation under the colorbar.
- Panel letters offset slightly outward.

Cluster/document labels no longer collide. Caption still matches panels A (θ grid) and B (cluster Δ).

## Other figures

Figure 3 (ECFP4) and Figure S1 (descriptive descriptor forest) were regenerated with the same numeric sources. Nested probability-scale results are tabulated in SI S5, not plotted as a replacement forest.
