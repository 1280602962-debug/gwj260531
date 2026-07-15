# Software license note (for Methods / compliance)

## Policy for this project

- Any Methods statement that names **Schrödinger / Glide / Prime / Desmond / QikProp** must correspond to runs performed under a **valid license** (institutional academic or commercial).
- Unlicensed copies must **not** be used to generate publishable results.
- If a legal Schrödinger environment is unavailable:
  1. Keep archived figures only if they were generated under license historically **and** the PI can attest; otherwise re-run with open tools; or
  2. Re-dock / re-simulate key molecules with **AutoDock Vina / Gnina / OpenMM / GROMACS** and rewrite Methods to those engines.
- **Never** relabel Glide poses or scores as another program.

## What reviewers usually ask

Reviewers typically ask for **version, precision mode, and parameters**, not a license PDF. Compliance risk is institutional/ethical, not primarily peer-review detection.

## Recommended Methods wording (when licensed)

> Molecular docking was performed with Glide (Schrödinger Suite release X.Y) using XP precision. Protein preparation and grid generation followed the Schrödinger Protein Preparation Wizard defaults unless noted. Licenses were provided by [institution].

## Recommended wording (open-tool re-run)

> Binding poses for purchased candidates were re-evaluated with AutoDock Vina / Gnina (version …) under the protocol in `docs/protocols/C2_C3_pose_md_replica_protocol.md`. Archived Glide ranks guided shortlisting but were not re-used as the sole pose evidence in the final manuscript.
