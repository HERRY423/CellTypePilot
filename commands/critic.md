---
description: Deep-review a specific cluster with CellTypePilot's Annotation Critic.
---

# Critic

Investigate a flagged or uncertain cluster using the Annotation Critic.

## Instructions

1. **Identify the target** — The user should specify which cluster to review. If not specified, check the evidence table for clusters with `needs_review` or `low` confidence.

2. **Run critic**:
   ```
   celltypepilot critic \
     --input <path> \
     --cluster-key <key> \
     --focus <cluster_id>
   ```

3. **Analyze the output**:
   - Top-5 candidate cell types with scores
   - Which markers support each candidate
   - Which negative markers conflict
   - Doublet likelihood assessment
   - Ontology consistency check

4. **Present findings**:
   - Show the evidence for and against the current annotation
   - Suggest alternative cell types if the current one is uncertain
   - Recommend actions: accept, re-annotate, sub-cluster, or mark as doublet

## When to use

- Critic flagged a cluster with `LOW_EVIDENCE`, `NEG_MARKER_CONFLICT`, or `POSSIBLE_DOUBLET`
- User questions a specific annotation
- Confidence is `low` or `needs_review`
- Two lineage signatures are co-expressed
