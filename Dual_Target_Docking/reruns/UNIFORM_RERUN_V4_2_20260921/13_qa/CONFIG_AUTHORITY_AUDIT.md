# CONFIG_AUTHORITY_AUDIT

{
  "primary_science_authority": "SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml + ablation_config.py",
  "alt_science_authority": "ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml",
  "analysis_freeze_method_fields": "STALE_PROVENANCE_ONLY",
  "current_runners_reading_analysis_freeze_yaml": [],
  "CRITICAL_CONFIG_DRIFT": false,
  "findings": [
    {
      "item": "analysis_freeze gnina.independent_docking_in_v4_1 false",
      "status": "STALE_PROVENANCE_ONLY",
      "drives_current_M3": false
    },
    {
      "item": "analysis_freeze allowed_pdb 4JPS/5DXT/4JSX",
      "status": "STALE_PROVENANCE_ONLY",
      "drives_formal_alt": false
    },
    {
      "item": "ablation freeze forbids the 8 PDBs that alt experiment executed",
      "status": "CHANNEL_SEPARATION_REQUIRED",
      "note": "alt channel uses ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL; primary M0-M3 use UNIQUE_14 only"
    },
    {
      "item": "ablation_config UNIQUE_14",
      "status": "ACTIVE_PRIMARY",
      "matches_scoring_freeze": true
    },
    {
      "item": "RTM checkpoint rtmscore_model1.pth",
      "status": "ACTIVE",
      "hardcoded_in": "ablation_config.py"
    },
    {
      "item": "GNINA v1.3.2 f23dd2b",
      "status": "ACTIVE",
      "hardcoded_in": "ablation_config.py + scoring freeze"
    }
  ]
}

Primary M0–M3 runners do not read analysis_freeze.yaml method fields.
Alt experiment is a separate authoritative channel.
No CRITICAL_CONFIG_DRIFT detected for receptor/ligand/method selection of current runners.
