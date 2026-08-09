# Host acceptance (Codex / Claude Code / MCP)

## Goal

Acceptance is **not** “the CLI binary can be invoked”.

An Agent host must **discriminate durable lifecycle states**:

| Agent state | Meaning | Must not |
|-------------|---------|----------|
| `running` | Unit still executing | Report final metrics |
| `completed` | Unit finished | Upgrade to public claim |
| `failed` | Started then errored/timeout | Hide error / invent success |
| `unavailable` | Dependency/adapter missing | Impute predictions |
| `claim_ready` | Release gate green | Apply outside frozen plan |
| `incomplete_not_claim_ready` | Gate blocked | Publish as claim-ready |
| `cancelled` | Cooperative stop | Report as completed |
| `resumed` | Completed via checkpoint resume | Drop provenance |

Commands under test: **`doctor`**, **`inspect`**, **`benchmark-run`** (plus discovery of observe / host-acceptance / MCP tools).

## What is checked

1. **Discovery**
   - Claude: `.claude-plugin/plugin.json`, `skills/celltypepilot/SKILL.md`
   - Codex: `.codex-plugin/plugin.json`, `skills/.../agents/openai.yaml`, `AGENTS.md`
   - MCP: `.mcp.json`, `tool_*` inventory + optional FastMCP runtime
2. **Parameter surface** — `benchmark-run --help` exposes fold worker flags
3. **JSON returns** — `doctor --json`, `inspect --json` serializable
4. **Lifecycle discrimination** — synthetic fixtures for all core states
5. **Timeout / cancel / resume** semantics (classifier + status files)
6. **Isolation** — fixtures under `scratch/host_acceptance/`; optional independent git worktree

## Run

```bash
# In-repo isolated fixtures
celltypepilot host-acceptance --skip-worktree --json

# Independent worktree (sibling clone)
python scripts/run_host_acceptance_worktree.py --json

# Explicit worktree path
celltypepilot host-acceptance --worktree ../celltypepilot-host-acceptance --json
```

Report path: `scratch/host_acceptance/host_acceptance_report.json` (or under the worktree).

## Agent host protocol (expected behavior)

```text
1. doctor --json
   → route on capabilities (available vs unavailable extras)
2. inspect --json
   → confirm species support / cluster keys (routing only)
3. benchmark-run ...
   → poll checkpoints OR observe --json
4. For each fold×method status:
   - running     → wait / surface ETA; no final table
   - completed   → may merge predictions for that unit
   - failed      → retain negative result; do not invent labels
   - unavailable → negative result; do not claim comparator ran
5. release readiness:
   - claim_ready only if release_manifest says so
   - incomplete_not_claim_ready blocks public robustness claims
```

MCP tool: `tool_agent_lifecycle_status(output_dir)` returns the same discrimination
vocabulary without mutating predictions.

## Pass criteria

- All discovery checks green
- Core five states map distinctly: running, completed, failed, unavailable, claim_ready
- Timeout classifies as **failed** (not unavailable)
- Cancel classifies as **cancelled**
- Resume classifies as **resumed** (or completed with resume flag)
- Report `overall_status == passed`
