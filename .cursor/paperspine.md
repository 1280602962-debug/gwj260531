# PaperSpine (Cursor)

Installed from https://github.com/WUBING2023/PaperSpine (v4.0.0).

- Skill path: `.cursor/skills/paper-spine/` and `~/.cursor/skills/paper-spine/`
- Entry: ask the agent to use the `paper-spine` skill / PaperSpine workflow
- Official hosts in upstream installer: Claude Code, Codex, OpenClaw, Hermes
- This Cursor install mirrors the Claude skill package into Cursor skill directories

To update:
```bash
cd /tmp && rm -rf PaperSpine && git clone --depth 1 https://github.com/WUBING2023/PaperSpine.git
bash PaperSpine/install.sh
cp -a ~/.claude/skills/paper-spine ~/.cursor/skills/
cp -a ~/.claude/skills/paper-spine /workspace/.cursor/skills/
```
