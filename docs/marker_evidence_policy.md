# Marker-edge evidence policy

PMID/DOI completeness and edge verification are different properties.

| Status | Meaning | Required fields |
|---|---|---|
| `aggregate_source_only_not_edge_verified` | The database publication is cited; this exact relation has not been traced | standard provenance fields |
| `database_record_verified` | A curator checked a stable database record for this gene/cell/tissue relation | record ID/URL, curator, verification date |
| `primary_source_verified` | A curator checked a primary source at a figure/table/supplement locator | primary source, evidence locator, curator, verification date |

Annotation policies:

- `database`: exploratory draft labels may use all structured relations;
- `edge_verified`: excludes aggregate-only relations;
- `primary`: uses only primary-source-verified relations.

The default does not upgrade the current atlas: its 733 bundled relationships remain
aggregate-only until curated. The critic emits `AGGREGATE_PROVENANCE_ONLY`, and reports must
not describe those relations as marker-specific experimental validation.

Adding an PMID found by an automated search is insufficient. A curator must verify the
cell type, gene polarity, species, tissue/state context, and the cited evidence locator.
