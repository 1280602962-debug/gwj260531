# LOCAL_ENVIRONMENT_PROVENANCE_ONLY

`FORMAL_AUTHORITY_PATHS_RESOLVED.json` records the absolute paths actually resolved on the original execution machine (`/tmp/pr39_fiveseed/...`).

It is **not** a current normative execution config.

Future analysis must resolve paths from:

- `00_protocol/FORMAL_AUTHORITY_PATHS.yaml` (PROJECT_ROOT-relative)
- `--project-root` or script-derived Dual_Target_Docking PROJECT_ROOT

Do not rewrite the historical absolute values in the JSON.
