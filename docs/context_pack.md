# Governed Context Pack

CellTypePilot uses a governed context interface for study-specific knowledge. Free text is
provenance, not evidence. Only explicit marker hypotheses can affect candidate scoring, and those
markers pass through the same statistical and critic gates as bundled atlas markers.

## JSON schema

Use `--context-file context.json` with schema `celltypepilot.context.v1`:

```json
{
  "schema_version": "celltypepilot.context.v1",
  "species": "human",
  "tissue": "kidney",
  "condition": "post-ischemic reperfusion injury",
  "timepoint": "24 h",
  "anatomical_region": "cortex",
  "free_text": "Pay attention to damaged proximal tubule states.",
  "review_status": "draft",
  "identity_hypotheses": [
    {
      "cell_type": "proximal tubule epithelial cell",
      "cl_id": "CL:0002306",
      "positive_markers": ["LRP2", "CUBN"],
      "negative_markers": ["PECAM1"],
      "source": "lab panel v2",
      "review_status": "reviewed"
    }
  ],
  "state_hypotheses": [
    {
      "state": "injury_responsive",
      "parent_cell_types": ["proximal tubule epithelial cell"],
      "positive_markers": ["VCAM1", "HAVCR1"],
      "negative_markers": [],
      "source": "lab panel v2",
      "review_status": "draft"
    }
  ]
}
```

`species` and `tissue` must match the run. `review_status` is `draft` or `reviewed`. Identity
ontology identifiers, when supplied, must match `CL:` followed by seven digits. A conflicting ID
for an existing atlas label is rejected.

## Custom marker CSV

`--custom-markers markers.csv` uses one marker per row:

```csv
axis,label,gene,polarity,cl_id,parent_cell_types,source,review_status
identity,proximal tubule epithelial cell,LRP2,positive,CL:0002306,,lab panel v2,reviewed
state,injury_responsive,VCAM1,positive,,proximal tubule epithelial cell,lab panel v2,draft
```

Required columns are `axis`, `label`, `gene`, and `polarity`. Allowed axes are `identity` and
`state`; allowed polarities are `positive` and `negative`. Separate multiple parent cell types
with semicolons.

## Trust boundary

- `--context` text alone cannot add candidates or unlock an unsupported tissue.
- A structured marker does not count merely because it is present in the matrix. Supporting
  evidence requires positive direction, log2FC and BH-FDR thresholds, and the configured cluster
  expression fraction.
- Expected-marker coverage includes every declared marker. Genes absent from the matrix and genes
  present but silent are reported separately.
- Draft context-only identity evidence must abstain. Reviewed context-only evidence is capped at
  medium confidence and remains visibly tagged.
- `manifest.json` records the schema, canonical context hash, source-file hashes, review status,
  scope counts, and whether free text was present. The optional normalized pack is written as
  `context_pack.normalized.json` for local audit.
