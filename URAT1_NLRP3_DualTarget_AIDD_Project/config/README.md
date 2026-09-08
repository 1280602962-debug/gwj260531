# Current manuscript configuration

- [campaign_final.yaml](campaign_final.yaml): current scientific route, evidence and interpretation.
- [docking_final.yaml](docking_final.yaml): consolidated GNINA parameters; original logs retain per-job authority.
- [chemistry_final.yaml](chemistry_final.yaml): actual historical exclusions, SMARTS and annotation roles.
- [urat1_gate_definition.yaml](urat1_gate_definition.yaml): A1/A2/W2, atom scope and residue mapping.
- [nlrp3_gate_definition.yaml](nlrp3_gate_definition.yaml): structural gate and validation denominators.
- [md_protocol.yaml](md_protocol.yaml): unexecuted plan and unresolved run parameters.
- [archive](archive/README.md): original C1–C5/development configurations and retirement registry.

Scientific policy YAML is not a batch-runner configuration. Only docking_final follows the engine schema; replicate seeds require explicit per-run overrides. Current evidence reconstruction does not execute docking or MD.
