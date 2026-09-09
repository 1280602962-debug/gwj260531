# Review of the pasted Chinese JCIM manuscript (2026-09-08)

This note reviews the Chinese Introduction–Conclusion text provided in
the query. It is **not** a rewrite of `docs/MANUSCRIPT_JCIM_ZH.md`.
Nature-skills used: `nature-writing` (paper-review, research paper type,
stance, nature-introduction, nature-results-discussion, zh-to-en /
Chinese-author workflow) and `nature-polishing` (failure-modes,
consistency-sweep, main-text-discipline, discussion-argument-language).
Journal axis: JCIM (not Nature Portfolio); treat Nature-skill checks as
discipline, not Nature house style.

## 1. Overall quality

The draft is scientifically clearer than the current assembled English
Results. The flagship claim is now correctly placed: **the same docking
scores change interpretation when the experimental-state comparison
changes**, and that change is pair- and arm-dependent. The four-state
design, pocket-matched scoring, ligand-only baseline, and the explicit
refusal to treat Dual-versus-neither as “the same task, harder negatives”
are all the right scientific spine.

It is a **research/evaluation paper**, not a methods paper that claims a
new scoring function. That is the right paper type. The remaining
problems are: leftover SI numbering from the pre-merge manuscript;
one Methods selection-rule inaccuracy; AI-template cadence in Chinese;
Results/Discussion overlap; and a few over-compressed or slightly
overclaimed sentences.

Quality verdict: **usable scientific draft, not yet submission-ready
Chinese/English JCIM text.**

## 2. AI-writing traces (moderate, not extreme)

This is not generic chatbot sludge. The argument is real. The AI
fingerprint is **template cadence + parallel scaffolding**, which is
exactly what `nature-polishing/failure-modes.md` flags.

Frequent patterns:

- Section openers restating the same claim (“考察…是否改变…判断”).
- Roadmap lists: 首先 / 随后 / 进一步 / 最后, repeated in Intro, 2.1,
  and the start of Results.
- Hedging stacks that are scientifically honest but rhythmically
  identical (“并不能…不能据此…不宜将…”).
- English leftovers in a Chinese main text: `docking`, `unused-pool
  holdout`, `pair`, unformatted `summarymin`.
- Topic–comment sentences that name the analysis before the finding.

Chinese-academic fix is not “more flowery Chinese.” It is: one claim per
paragraph, finding first, then the test. Cut the protocol-list endings.
Define `docking` / AUROC / \(\mathrm{summary}_{\min}\) once, then use
Chinese terms (`对接评分`、`判别`、`较弱方向汇总`).

Highest-priority prose cuts:

1. Intro last paragraph: keep the question, delete the four-step protocol.
2. Delete or shrink §2.1 to four sentences; it currently reprints the
   Intro.
3. §3.2: lead with EGFR/HER2 pocket A (0.430 vs 0.808), then JAK1/TYK2
   as replication, then the pairs where the gap is compatible with zero.
   Do not open with the eight-pair `summary_min` range.
4. §4.1: stop repeating AUROC numbers. State the task distinction, then
   one sentence that the magnitude is pair-dependent (point to Table 3 /
   Table S4).
5. Conclusion: independent GNINA tests **formulation-gap persistence**,
   not “matched-pocket / structure correspondence.” Do not fold it into
   the same sentence as pocket-swap evidence.

## 3. Framework

The 3.1 supply → 3.2 formulation → 3.3 chemistry → 3.4 structure →
3.5 labels/membership → 3.6 external-impossibility order is the correct
diagnostic chain. Keep it.

What to change structurally (Nature main-text discipline):

| Keep in main text | Move to SI pointer only |
|---|---|
| Four-state definition and pocket-matched scoring | Full bootstrap algorithm |
| Negative-class swap on a fixed score | Exhaustiveness grid, five-seed table |
| ECFP4 increment (max \|Δ\| = 0.023) | Per-descriptor AUROC grid |
| Matched vs mismatched pocket; unused-pool | Document-cluster CIs, power simulation |
| One PIK3CA/mTOR crystal swap | Full swap numeric table |
| Independent GNINA on two pairs | RTMScore/CNN per-pair grid |
| BindingDB gate = 0; time-split impossible | Full gate tables |

§2.2 mentions a high-confidence human single-protein relabel, but §3.5
only reports max-versus-median. Either report the high-confidence result
in one sentence or drop it from Methods.

§2.3 currently says bidirectional supply used **A-only and B-only ≥ 50**
as the operational gate, then lists all eight pairs. That was the original
thick-panel freeze. The five census pairs were **not** admitted by that
≥50 strict gate; EGFR/HER2 also failed it and was retained as
supply-limited. Rewrite the selection paragraph: thick-panel pairs versus
census pairs, and say EGFR/HER2 is the supply-limited exception.

## 4. Data and conclusion audit

Headline numbers in the pasted Results **match the frozen CSVs**
(`unified_threshold_sensitivity_v2.csv`,
`five_pair_stack_v1/table2_comparable_theta6_v1.csv`,
`formulation_equal_score_negative_v1.csv`, SI Tables S2–S13). Census
funnel 2,164,618 / 63,790 / 5,253 / 86 matches
`universe_census_summary_v1.csv`.

| Claim | Verdict |
|---|---|
| `summary_min` range 0.345–0.692 | Correct (F2/F10, PIK3CA/mTOR) |
| Only PPARG/PPARA CI entirely >0.5: 0.649 [0.504, 0.751] | Correct |
| Other CIs contain 0.5 or lie entirely below | Correct (F2/F10 [0.211, 0.477] is entirely <0.5) |
| EGFR pocket A: 0.430 vs 0.808, Δ 0.378 [0.205, 0.547] | Correct |
| JAK1/TYK2 pocket A Δ 0.444 [0.263, 0.620] | Correct |
| EGFR / JAK1/TYK2 Dual-vs-neither 0.756 / 0.770 | Correct (`vina_mean`) |
| EGFR Top-10 1 Dual / 9 selectives / 0 neither | Correct |
| AND filter: 14 Dual + 33 selectives, precision 0.298 | Correct |
| AChE TPSA 0.733 / 0.801; PIK3CA/mTOR heavy-atom 0.463 | Correct |
| ECFP4+docking, 16 arms, max \|Δ\| = 0.023 | Correct |
| EGFR 3POZ top-1 9.505 Å, best-of-9 0.760 Å | Correct |
| Matched−mismatched 0.170 and 0.161; holdout −0.079 to +0.150 | Correct |
| 4JPS `summary_min` 0.486 [0.259, 0.692] | Correct; also state 5DXT 0.505 and 4JSX 0.639 or point to Table S8 |
| Independent GNINA EGFR 0.783 vs 0.220; JAK1/TYK2 0.705 vs 0.317 | Correct |
| PPARG RTMScore 0.369; CNN 0.500; holdout 0.535 | Correct |
| EGFR strict B-only n = 7, `summary_min` 0.324 | Correct |
| Max vs median agreement 93.6%–100%, max \|Δ\| = 0.023 | Correct |
| Holdout points 0.618 / 0.765 / 0.392 / 0.475 / 0.619 / 0.535 / 0.445 | Correct |
| Document-cluster PIK3CA/mTOR Dual-vs-B-only [0.000, 0.818] | Correct |
| BindingDB gate 0 pairs; 2018 time-split no bidirectional test set | Correct |

Mechanical errors to fix before any submission:

1. **SI numbers are the old series.** After the SI merge, Tables S25, S28,
   S29, S30, S32, S34, S39–S41, S46, S48–S49, S54 do not exist. Map to:
   S4 (negative-class), S13 (EGFR ranking/AND), S2 (cognate RMSD; boxes
   and RMSD are now the same table), S3 (θ and aggregation), S5 (ECFP4),
   S6 (matched pocket), S7 (holdout), S8 (crystal swap), S9 (GNINA /
   rescore / seeds), S10 (document-cluster / power), S11 (BindingDB),
   S12 (time-split). Methods §2.7 still splits boxes = S2 and cognate = S3;
   both are Table S2.
2. “JAK1/TYK2 holdout 0.475，略有上升” understates 0.365 → 0.475.
3. §3.6 leftover English “pair”.
4. Conclusion over-packs independent GNINA into “结构对应优势”.

Scientific conclusions in §§4–5 are **directionally correct** and match
the locked interpretation: Dual-versus-neither is a different estimand;
ligand chemistry can carry the signal; matched-pocket advantage is not
stable across pairs or unused-pool ligands; retrospective AUROC is not
prospective hit rate. Do not strengthen those claims.

## 5. JCIM: is it worth trying?

**Yes, as a tightly scoped evaluation Article. No, as a new docking
method paper.**

JCIM will not treat this as a scoring-function advance: there is no new
algorithm, no superiority table, and directional AUROCs are mostly
compatible with 0.5. The publishable object is the **evaluation finding**:
dual-target docking evidence is not invariant to experimental-state
definition, and apparent AUROC is not automatically pocket information.

What improved since the 2026-08-26 four-pair internal note: eight pairs,
ligand-chemistry baseline, unused-pool, BindingDB gate, independent
GNINA on two pairs. What has not: no wet-lab, no independent external
pair that passed the gate, Vina-centric main analysis, Zenodo DOI still
intentionally not minted.

Desk-reject risk remains real if the cover letter sells “a dual-target
docking benchmark that works.” It is lower if the cover letter sells
“how dual-target docking evidence should be defined, shown on eight
pairs.” Do not submit until:

- English JCIM manuscript exists (the journal is English);
- SI citations match Tables S1–S13;
- Zenodo/tag deposit is ready;
- claim language stays inside the locked bounds (no PIK3CA/PIK3CB
  performance; 9V8H peptide retained; BindingDB = not docked).

Honest one-line: **可以按 evaluation article 投稿试一次，但不要按新对接方法投；中文稿只适合组内定稿，外投必须是英文稿加归档。**
