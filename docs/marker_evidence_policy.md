# Marker-edge evidence policy

PMID/DOI completeness and edge verification are different properties.

| Status | Meaning | Required fields |
|---|---|---|
| `aggregate_source_only_not_edge_verified` | The database publication is cited; this exact relation has not been traced | standard provenance fields |
| `literature_cooccurrence_supported` | PubMed title/abstract co-occurrence (>=2 hits) supports the gene/cell-type association; NOT a claim of marker specificity or primary-experiment review | evidence locator (query, hit count, top PMIDs), curator, verification date |
| `database_record_verified` | A curator checked a stable database record for this gene/cell/tissue relation | record ID/URL, curator, verification date |
| `primary_source_verified` | A curator checked a primary source at a figure/table/supplement locator | primary source, evidence locator, curator, verification date |

Annotation policies:

- `database`: exploratory draft labels may use all structured relations;
- `literature`: requires at least literature co-occurrence support;
- `edge_verified`: excludes aggregate-only and co-occurrence-only relations;
- `primary`: uses only primary-source-verified relations.

Co-occurrence upgrades are produced by the auditable sweep
(`celltypepilot curate`): every upgrade records the exact PubMed query, hit
count, top PMIDs, timestamp, and curator identity, and the merged atlas must
pass provenance validation before it is written (fail closed).

Cell Ontology identifiers are checked against the live Cell Ontology
(`celltypepilot ontology update` then `celltypepilot ontology check`):
unknown or obsolete CL identifiers are errors, lexical label mismatches are
warnings (the atlas key may be an intentional refinement).

Adding an PMID found by an automated search is insufficient for edge or
primary verification. A curator must verify the cell type, gene polarity,
species, tissue/state context, and the cited evidence locator.
