# Lung marker-edge promotion packet

Status: `pending_human_review`; no runtime atlas was changed.

This deliberately small packet contains two high-value edge upgrades:

1. `CA4 -> Capillary endothelial cell`: proposed as
   `primary_source_verified` from adult human lung endothelial single-cell work
   with serial immunofluorescence localization to capillaries.
2. `MYH11 -> Bronchial smooth muscle cell`: proposed as
   `database_record_verified` from the HPA lung single-cell record `c-29`.

Both proposals were drafted by an automated search and therefore remain
candidate-only. `evidence_promotion.add_promotion_review` requires two distinct
human reviewers; only `apply_approved_promotion` may create a new versioned
atlas after those approvals.

`LUM -> Pulmonary interstitial fibroblast` was reviewed but not proposed. The
located evidence supports lung/alveolar/fibroblast association, not the exact
CellTypePilot identity label, so promoting it would overstate specificity.

Benchmark leakage boundary: the Travaglini Smart-seq2 rerun must use the source
atlas, not these unapplied proposals. This packet is not benchmark evidence and
must not be used to tune labels on the held-out evaluation dataset.
