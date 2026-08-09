# Azimuth reference governance

Azimuth is an external-reference method, not a fold-trained comparator. A successful
`RunAzimuth` call is therefore insufficient for the independent primary track.

Before any prediction is accepted, `reference_audit.json` must freeze:

- official reference name, version or publication date, source URL, retrieval UTC, and
  SHA-256/byte size for `ref.Rds`, `idx.annoy`, and accompanying metadata;
- species, tissue, assay/technology, gene identifier space, label levels, and software
  versions required to load the reference;
- every disclosed training study, DOI/accession, donor or sample identifier when
  available, and whether provenance is complete enough to audit;
- exact and possible overlap with each evaluation study, dataset version, donor, and
  sample;
- license/redistribution status and whether the frozen files may be shipped or only
  referenced by hash;
- a primary-track status and reason using the states below.

Primary-track states:

- `eligible_external_reference`: tissue-compatible and training provenance is complete,
  with no evaluation-study/donor/sample overlap found;
- `ineligible_training_study_overlap`: the evaluation study contributed to the reference;
- `ineligible_donor_or_sample_overlap`: a held-out donor/sample contributed to the
  reference even if study names differ;
- `ineligible_provenance_incomplete`: independence cannot be established from disclosed
  provenance;
- `unavailable_no_tissue_reference`: no official tissue-compatible reference exists;
- `unavailable_asset_or_dependency_failure`: required files or software cannot be loaded.

Anything other than `eligible_external_reference` is retained as a negative result and
must not be imputed, replaced by a generic reference, or described as fold-trained. An
eligible Azimuth result belongs to a separately named external-reference track; it does
not establish leave-one-donor-out training equivalence with CellTypePilot, CellTypist,
SingleR, or popV.
