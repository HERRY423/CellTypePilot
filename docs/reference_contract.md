# Reference safety contract

Plugin-facing annotation fails closed when a reference cannot be matched to the query.

Custom `.h5ad` references must contain:

```python
reference.uns["celltypepilot_reference"] = {
    "species": "human",
    "tissues": ["blood"],
    "source": "GEO:GSE...",
    "version": "2026-08-locked",
    "label_ontology": "Cell Ontology 2026-07",
    "training_studies": ["GSE...", "GSE..."],
}
```

Explicit CellTypist model files require an adjacent `<model>.pkl.json` sidecar containing
`species`, `tissues`, `source`, `version`, `labels`, and the model `sha256`.

Automatic model selection is intentionally narrow. At present only the registered human
blood `Immune_All_Low.pkl` scope is eligible. No immune model is silently applied to lung,
brain, tumor, or another unmatched tissue. `--allow-unverified-reference` is an explicit
escape hatch and is recorded in `manifest.json`; it does not turn the reference into verified
evidence.
