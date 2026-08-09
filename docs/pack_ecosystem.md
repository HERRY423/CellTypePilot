# Evidence / reference pack ecosystem

Packs are the sustainable contribution surface for CellTypePilot: **signed,
versioned, data-only** tissue/disease knowledge that still passes ordinary
runtime gates.

## Hard rules

| Rule | Detail |
|------|--------|
| Data only | Only `pack.json`, `marker_atlas.json`, `state_atlas.json`, `LICENSE`, `README.md`, `ontology_map.json`, `reference_manifest.json`, `pack.sig.json` |
| No code | `.py`, `.so`, `.sh`, notebooks, nested trees → install/sign fail closed |
| Versioned | `name` + `version` + content SHA-256 fingerprint |
| License / provenance / ontology | Declared in `pack.json`; edges carry sources + verification_status |
| Runtime gates unchanged | Every pack marker still hits DE evidence, critic, abstention, conflict detection |
| Trust tiers | `atlas` (full provenance) vs `hypothesis` (draft; cannot silently accept identity) |

## Pack kinds

- `evidence` — marker/state knowledge for scoring  
- `reference` — curated reference-oriented panels (still data-only)  
- `mixed` — both  

Tissue + disease labels (`tissues`, `diseases`) enable specialized packs without
forking the core plugin.

## Contributor workflow

```bash
# 1) Scaffold
celltypepilot pack scaffold my-ibd-pack -o ./my-ibd-pack \
  --tissue gut --disease IBD --kind evidence --license CC-BY-4.0

# 2) Edit marker_atlas.json (add cl_id + marker_evidence per edge)

# 3) Sign (HMAC dev default, or RSA private key)
celltypepilot pack sign ./my-ibd-pack --signer "lab-curator@example.org"

# 4) Verify
celltypepilot pack verify ./my-ibd-pack --require-signature

# 5) Install (fail closed unless --trust hypothesis)
celltypepilot pack install ./my-ibd-pack

# 6) Use — still ordinary gates
celltypepilot annotate -i data.h5ad -k leiden --pack my-ibd-pack -o out/
```

## Runtime guarantee

Installing a pack never injects Python. Annotation merge tags incomplete packs as
context/hypothesis so the critic and abstention path cannot “trust-bypass” them
into high-confidence identity without evidence.

## Human review after pack-assisted annotation

See Web Inspector **Evidence** panel: Identity × State × Novelty, supporting vs
opposing markers, neighbors, donor/batch strata, literature hooks.

```bash
# After apply overrides (append-only audit)
celltypepilot review-resign -o out/ --signer "Dr. Reviewer"
```

Stale derived artifacts clear only after regenerate + re-sign.
