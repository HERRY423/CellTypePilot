---
description: Check CellTypePilot environment — dependencies, capabilities, MCP status.
---

# Doctor

Verify that CellTypePilot is properly installed and report the capability tier.

## Instructions

1. **Run doctor**:
   ```
   celltypepilot doctor
   ```

2. **Interpret results**:
   - **Core dependencies**: Python version, scanpy, anndata, matplotlib — MUST all pass
   - **Optional dependencies**: Flask (web inspector), rpy2 (Seurat .rds), etc.
   - **Capability tier**: `full` (all deps), `degraded` (core only), or `broken` (core missing)
   - **MCP / Literature**: PubMed direct access, MCP server status
   - **License**: Current license tier (free/academic/commercial/trial)

3. **If core deps are missing**:
   - Tell the user to run `pip install celltypepilot` or `pip install -e .`
   - Do NOT attempt annotation without passing the doctor check

4. **If optional deps are missing**:
   - Inform the user which features are unavailable
   - Suggest `pip install celltypepilot[web]` for Web Inspector
   - Suggest `pip install celltypepilot[seurat]` for Seurat .rds support

## Notes

- Run this FIRST before any annotation workflow
- The doctor check is non-destructive and fast
