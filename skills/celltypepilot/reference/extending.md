# CellTypePilot — Extending the Atlas

> How to add custom markers, new tissues, or expand the knowledge graph.

## 1. Adding Custom Markers

### 1.1 Via JSON

Create a custom marker file following the atlas schema:

```json
{
  "custom_tissue": {
    "name": "My Custom Tissue",
    "organ_system": "custom",
    "cell_types": {
      "My Cell Type": {
        "cl_id": "CL:9999999",
        "synonyms": ["alias1", "alias2"],
        "positive_markers": ["GENE1", "GENE2", "GENE3"],
        "negative_markers": ["GENE4", "GENE5"],
        "subtypes": {}
      }
    }
  }
}
```

### 1.2 Merging with Built-in Atlas

To use custom markers alongside the built-in atlas, merge them programmatically:

```python
import json
from celltypepilot.data_adapter import load_marker_atlas

# Load built-in atlas
atlas = load_marker_atlas("human")

# Load custom markers
with open("my_markers.json") as f:
    custom = json.load(f)

# Merge
atlas["tissues"].update(custom)

# Now use the merged atlas for scoring
```

## 2. Adding New Tissues

Add a new tissue entry to the atlas JSON. Each tissue needs:
- `name`: Display name
- `organ_system`: Organ system category
- `cell_types`: Dict of cell type definitions

Each cell type needs:
- `cl_id`: Cell Ontology ID (format: CL:XXXXXXX)
- `synonyms`: List of alternative names
- `positive_markers`: Genes whose expression supports this cell type
- `negative_markers`: Genes whose expression argues against this cell type
- `subtypes`: Optional nested dict of subtype definitions (same schema)

## 3. Marker Selection Guidelines

### Good positive markers:
- Expressed in >50% of the target cell type
- Relatively specific (not broadly expressed across many types)
- Well-established in literature
- Detectable by scRNA-seq (membrane/cytosolic proteins preferred)

### Good negative markers:
- Known to be absent from the target type
- Expressed in closely related types that might be confused
- Well-established in literature

### Marker quality tiers:
1. **Gold**: Canonical markers universally accepted (e.g., CD3D for T cells)
2. **Silver**: Well-supported but context-dependent (e.g., tissue-specific subtypes)
3. **Bronze**: Literature-supported but may vary by condition/species

## 4. Species Extension

To add support for a new species:

1. Create gene symbol mappings in the `mouse_gene_map` format
2. Add species-specific marker variants if needed
3. Update the `detect_species()` function in `data_adapter.py`

## 5. Contributing to the Built-in Atlas

The built-in atlas prioritizes:
1. **Breadth**: Cover common tissues across human and mouse
2. **Quality**: Only well-established markers from curated databases
3. **Usability**: Markers that work well with scRNA-seq data

Contributions should include:
- Marker gene symbols (HGNC for human)
- Source reference (database/literature)
- Evidence level (canonical / well-supported / emerging)
- Tissue and cell type context
