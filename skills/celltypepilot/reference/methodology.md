# CellTypePilot — Methodology Reference

> Detailed algorithm description for the marker scoring and critic pipeline.

## 1. Marker Scoring Algorithm

### 1.1 Differential Expression

For each cluster, CellTypePilot computes differential expression using the Wilcoxon rank-sum
test (scanpy's `rank_genes_groups` with `method="wilcoxon"`). This produces per-gene:
- Log2 fold change (logfoldchange)
- P-value and adjusted p-value

### 1.2 Marker Score Components

For each cluster × cell_type combination, five metrics are computed:

**a) Overlap (pct_overlap)**
```
pct_overlap = n_DE_markers / n_detected_markers
```
Fraction of the cell type's positive markers (that exist in the data) that are also
differentially expressed in this cluster. Range: [0, 1].

**b) Expression magnitude (mean_log2fc)**
```
mean_log2fc = mean(log2fc of DE positive markers)
```
Average fold change of the detected positive markers that are DE. Normalized to [0, 1]
by dividing by 2.0 and clipping.

**c) Expression breadth (pct_expressed)**
```
pct_expressed = n_expressed_markers / n_detected_markers
```
Fraction of positive markers expressed in ≥25% of cells within the cluster.

**d) Specificity**
```
specificity = mean(n_cells_expressing_in_cluster / n_cells_expressing_total) per marker
```
For each marker gene, what fraction of all expressing cells belong to this cluster?
Averaged across all positive markers. High specificity means the marker is selective
for this cluster.

**e) Negative conflict (neg_conflict)**
```
neg_conflict = n_neg_expressed / n_neg_detected
```
Fraction of negative markers (that should NOT be expressed) that ARE expressed in >15%
of cluster cells. This is a penalty term.

### 1.3 Combined Score

```
combined = 0.35 × pct_overlap
         + 0.25 × min(mean_log2fc / 2.0, 1.0)
         + 0.20 × specificity
         + 0.20 × pct_expressed
         - 0.30 × neg_conflict
```

Clipped to [0, 1]. The weights prioritize marker overlap and expression magnitude,
while penalizing negative marker violations.

### 1.4 Confidence Assignment

| Confidence | Conditions |
|---|---|
| high | score ≥ 0.7, overlap ≥ 0.5, neg_conflict < 0.1, specificity ≥ 0.5 |
| medium | score ≥ 0.5, overlap ≥ 0.3, neg_conflict < 0.2 |
| low | score ≥ 0.3 |
| needs_review | score < 0.3 |

## 2. Annotation Critic

The Critic performs four independent checks on each annotation:

### 2.1 Evidence Sufficiency

Counts how many of the cell type's positive markers are:
1. Present in the data (gene exists in var_names)
2. Expressed in ≥25% of cluster cells

Coverage < 20% → `LOW_EVIDENCE` flag
Coverage 20-50% → `PARTIAL_EVIDENCE` flag

### 2.2 Negative Marker Conflict

For each negative marker defined for the cell type:
- Check if it's expressed in >15% of cluster cells
- If any are, flag as `NEG_MARKER_CONFLICT`

This catches cases where the annotation is contradicted by markers that should be absent.

### 2.3 Doublet / Mixed Signal Heuristic

For the cluster being annotated, check ALL cell types in the atlas:
- Compute marker coverage for each cell type
- If two distinct lineage types both have ≥30% marker coverage, flag as `POSSIBLE_DOUBLET`

This detects cases where two mutually exclusive lineage programs are simultaneously active,
suggesting either a doublet or a transitional state.

### 2.4 Ontology Consistency

Basic checks:
- CL ID format: must match `CL:\d{7}`
- CL ID presence: flag if no Cell Ontology ID assigned

### 2.5 Confidence Recalibration

The critic can only **downgrade** confidence, never upgrade it:
- `LOW_EVIDENCE` → cap at `needs_review`
- `PARTIAL_EVIDENCE` → cap at `low`
- `NEG_MARKER_CONFLICT` → cap at `needs_review`
- `POSSIBLE_DOUBLET` → cap at `needs_review`

## 3. Marker Knowledge Graph

### 3.1 Data Model

```
CellType(cell_ontology_id, name, synonyms[])
Marker(gene_symbol, species)
Tissue(name)

Edge: CellType --marked_by(polarity: pos/neg, specificity_score)--> Marker
Edge: CellType --observed_in--> Tissue
Provenance: source, evidence_snippet_ref, confidence_tier, kg_version
```

### 3.2 Sources

The built-in atlas integrates markers from:
- **PanglaoDB**: Curated markers across tissues and species
- **CellMarker 2.0**: Comprehensive cell type marker database
- **Cell Ontology (CL)**: Standardized cell type identifiers
- **Literature curation**: Additional markers from key publications

### 3.3 Version Management

Each atlas release is versioned as `mkg-YYYY.MM`. Every run records which MKG version
was used in the manifest.json, ensuring reproducibility.

## 4. Species Support

### 4.1 Auto-Detection

Gene naming conventions differ between species:
- Human: ALL CAPS (e.g., CD3D, MS4A1)
- Mouse: First letter capitalized (e.g., Cd3d, Ms4a1)

CellTypePilot samples the first 500 gene names and classifies based on the dominant pattern.

### 4.2 Mouse Gene Conversion

Most human-to-mouse conversions follow a simple rule: capitalize first letter, lowercase rest.
Key exceptions (CD molecules, etc.) are handled by a lookup table in the atlas.
