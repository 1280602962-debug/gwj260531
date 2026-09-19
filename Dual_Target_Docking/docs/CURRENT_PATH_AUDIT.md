# Current-path audit

Date: 2026-09-19
Scanner: `scripts/qa/audit_current_paths.py`
Rule: every current-facing source/protocol/script path must exist, or be marked historical / archived / unavailable.
Checksum manifests are integrity records, not scientific PASS/FAIL authority.

- files scanned: 33
- path citations: 349
- exists: 325
- historical_or_unavailable: 24
- unresolved broken current paths: **0**

No unresolved broken current-facing paths.

## Verdict

path audit unresolved = 0

Machine table: `results/qa/current_path_audit.csv`.
