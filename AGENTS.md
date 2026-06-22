# AGENTS.md

Guidance for AI agents working in this repository.

## Repository overview

This repository hosts a **single academic manuscript artifact**:

- `1_Results_and_Discussion_translated_checked.md` — scientific review, bilingual alignment, and polished English for a CADD paper (dual-target PLK1 + NLRP3 inhibitors).

There is **no application source code**, package manifest, Docker Compose, database, or CI configuration in this repo. The full computational pipeline (Python/RDKit, docking, AMBER MD, etc.) is described in the document but lives in separate local/HPC projects, not here.

## Cursor Cloud specific instructions

### What to run

| Goal | Command / approach |
|------|-------------------|
| Read or edit the manuscript | Open `1_Results_and_Discussion_translated_checked.md` in the editor |
| Validate structure (9 sections 3.1–3.9) | `grep '^### Section' 1_Results_and_Discussion_translated_checked.md` |
| Optional Markdown style check | `npx --yes markdownlint-cli2 "1_Results_and_Discussion_translated_checked.md"` (expect MD013/MD030 on long academic lines; not a repo gate) |
| HTML preview (local only) | See **Manuscript preview** below |

No dev servers, databases, or `npm install` / `pip install -r` steps are required from files in this repo.

### Manuscript preview

To render and preview without changing tracked files:

```bash
pip install --user markdown
export PATH="$HOME/.local/bin:$PATH"
python3 -c "
import markdown
from pathlib import Path
p = Path('1_Results_and_Discussion_translated_checked.md')
html = markdown.markdown(p.read_text(encoding='utf-8'), extensions=['tables', 'fenced_code'])
Path('/tmp/manuscript-preview').mkdir(exist_ok=True)
Path('/tmp/manuscript-preview/index.html').write_text(
    '<!DOCTYPE html><html><head><meta charset=utf-8></head><body>' + html + '</body></html>',
    encoding='utf-8')
print('Wrote /tmp/manuscript-preview/index.html')
"
```

Serve with tmux (Cloud Agent VMs):

```bash
SESSION_NAME=manuscript-preview
tmux -f /exec-daemon/tmux.portal.conf new-session -d -s "$SESSION_NAME" -c /tmp/manuscript-preview -- bash -l
tmux -f /exec-daemon/tmux.portal.conf send-keys -t "$SESSION_NAME:0.0" \
  'python3 -m http.server 8765 --bind 127.0.0.1' C-m
curl -s -o /dev/null -w 'HTTP %{http_code}\n' http://127.0.0.1:8765/
```

Preview output is under `/tmp/manuscript-preview/` and is not part of the git tree.

### Lint / test / build

- **Lint**: No project linter is configured. Optional `markdownlint-cli2` via `npx` (see table above).
- **Test**: None in-repo.
- **Build**: None in-repo.

### If you expected the CADD codebase

The manuscript references paths such as `Document_PLK1 and NLRP3` on the author’s machine. Clone or attach that project separately if you need to run QSAR, docking, or MD workflows.
