# Minimum technical verification: travaglini_lung_smartseq2_2020

Status: **minimum_verification_complete**

This is a real single-study, single-platform, three-donor technical verification. It is not the public multi-cohort benchmark release.

## Donor-weighted results

| method | metric | estimate | ci_lower | ci_upper | se | n_donors | n_studies | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| celltypepilot | accuracy | 0.0 | 0.0 | 0.0 | 0.0 | 3 | 1 | estimated |
| celltypepilot | macro_f1 | 0.0 | 0.0 | 0.0 | 0.0 | 3 | 1 | estimated |
| celltypepilot | balanced_accuracy | 0.0 | 0.0 | 0.0 | 0.0 | 3 | 1 | estimated |
| celltypepilot | coverage | 0.0 | 0.0 | 0.0 | 0.0 | 3 | 1 | estimated |
| celltypepilot | selective_accuracy |  |  |  |  | 0 | 0 | not_estimable |
| celltypist | accuracy | 0.9067883807968582 | 0.7975921745673439 | 0.978612893369997 | 0.04458487101440844 | 3 | 1 | estimated |
| celltypist | macro_f1 | 0.776129753671052 | 0.7280780643635026 | 0.817626907784181 | 0.021191844998786707 | 3 | 1 | estimated |
| celltypist | balanced_accuracy | 0.8014606105329546 | 0.7263162329547362 | 0.8622632529365927 | 0.03172309829864837 | 3 | 1 | estimated |
| celltypist | coverage | 1.0 | 1.0 | 1.0 | 0.0 | 3 | 1 | estimated |
| celltypist | selective_accuracy | 0.9067883807968582 | 0.7975921745673439 | 0.978612893369997 | 0.04458487101440844 | 3 | 1 | estimated |
| popv | accuracy | 0.927800162440677 | 0.8620516679207424 | 0.9804460739382829 | 0.028106713123496017 | 3 | 1 | estimated |
| popv | macro_f1 | 0.8537159853002564 | 0.819232394109393 | 0.8849319231710173 | 0.015500780028337064 | 3 | 1 | estimated |
| popv | balanced_accuracy | 0.8741682384466909 | 0.8319920077559733 | 0.9070133958937013 | 0.017626297130228695 | 3 | 1 | estimated |
| popv | coverage | 1.0 | 1.0 | 1.0 | 0.0 | 3 | 1 | estimated |
| popv | selective_accuracy | 0.927800162440677 | 0.8620516679207424 | 0.9804460739382829 | 0.028106713123496017 | 3 | 1 | estimated |
| singler | accuracy | 0.9225594277846582 | 0.8753448708301982 | 0.9709746410021388 | 0.022460687250738118 | 3 | 1 | estimated |
| singler | macro_f1 | 0.8658173155303417 | 0.8043132322659298 | 0.906100582911905 | 0.02521861727770684 | 3 | 1 | estimated |
| singler | balanced_accuracy | 0.8751599674778564 | 0.7949981697212801 | 0.9164228318891481 | 0.032182019098467116 | 3 | 1 | estimated |
| singler | coverage | 0.9765830446898658 | 0.9651000465332714 | 0.9874732661167126 | 0.005270904390917585 | 3 | 1 | estimated |
| singler | selective_accuracy | 0.9445865287903392 | 0.8957905544147845 | 0.9832920792079208 | 0.020798686374168766 | 3 | 1 | estimated |

## Paired donor-level comparisons

| method_a | method_b | metric | n_paired_donors | mean_difference | median_difference | difference_ci_lower | difference_ci_upper | win_fraction | p_value | status | p_value_bh | reject_fdr_005 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| celltypepilot | celltypist | macro_f1 | 3 | -0.776129753671052 | -0.7826842888654724 | -0.817626907784181 | -0.7280780643635026 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | popv | macro_f1 | 3 | -0.8537159853002564 | -0.8569836386203594 | -0.8849319231710173 | -0.819232394109393 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | singler | macro_f1 | 3 | -0.8658173155303417 | -0.8870381314131904 | -0.906100582911905 | -0.8043132322659298 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | celltypist | balanced_accuracy | 3 | -0.8014606105329546 | -0.8158023457075351 | -0.8622632529365927 | -0.7263162329547362 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | popv | balanced_accuracy | 3 | -0.8741682384466909 | -0.8834993116903984 | -0.9070133958937013 | -0.8319920077559733 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | singler | balanced_accuracy | 3 | -0.8751599674778564 | -0.9140589008231407 | -0.9164228318891481 | -0.7949981697212801 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | celltypist | coverage | 3 | -1.0 | -1.0 | -1.0 | -1.0 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | popv | coverage | 3 | -1.0 | -1.0 | -1.0 | -1.0 | 0.0 | 0.25 | estimated | 0.25 | False |
| celltypepilot | singler | coverage | 3 | -0.9765830446898658 | -0.9771758214196138 | -0.9874732661167126 | -0.9651000465332714 | 0.0 | 0.25 | estimated | 0.25 | False |

## Batch sensitivity

| axis | method | metric | status | n_supported_levels |
| --- | --- | --- | --- | --- |
| platform | celltypepilot | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| platform | celltypist | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| platform | popv | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| platform | singler | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| plate_batch | celltypepilot | macro_f1 | not_estimable_lt2_supported_levels | 0 |
| plate_batch | celltypist | macro_f1 | not_estimable_lt2_supported_levels | 0 |
| plate_batch | popv | macro_f1 | not_estimable_lt2_supported_levels | 0 |
| plate_batch | singler | macro_f1 | not_estimable_lt2_supported_levels | 0 |
| condition | celltypepilot | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| condition | celltypist | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| condition | popv | macro_f1 | not_estimable_lt2_supported_levels | 1 |
| condition | singler | macro_f1 | not_estimable_lt2_supported_levels | 1 |

## QC-stratified performance

| diagnostic | method | stratum | donor_unit | metric | estimate | n_cells | status | identity_effect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| low_quality | celltypepilot | not_flagged | travaglini_2020::travaglini_2020::1 | macro_f1 | 0.0 | 3987.0 | estimated_donor_stratum | none |
| low_quality | celltypepilot | not_flagged | travaglini_2020::travaglini_2020::2 | macro_f1 | 0.0 | 3273.0 | estimated_donor_stratum | none |
| low_quality | celltypepilot | not_flagged | travaglini_2020::travaglini_2020::3 | macro_f1 | 0.0 | 2149.0 | estimated_donor_stratum | none |
| low_quality | celltypist | not_flagged | travaglini_2020::travaglini_2020::1 | macro_f1 | 0.7826842888654724 | 3987.0 | estimated_donor_stratum | none |
| low_quality | celltypist | not_flagged | travaglini_2020::travaglini_2020::2 | macro_f1 | 0.817626907784181 | 3273.0 | estimated_donor_stratum | none |
| low_quality | celltypist | not_flagged | travaglini_2020::travaglini_2020::3 | macro_f1 | 0.7280780643635026 | 2149.0 | estimated_donor_stratum | none |
| low_quality | popv | not_flagged | travaglini_2020::travaglini_2020::1 | macro_f1 | 0.8569836386203594 | 3987.0 | estimated_donor_stratum | none |
| low_quality | popv | not_flagged | travaglini_2020::travaglini_2020::2 | macro_f1 | 0.8849319231710173 | 3273.0 | estimated_donor_stratum | none |
| low_quality | popv | not_flagged | travaglini_2020::travaglini_2020::3 | macro_f1 | 0.8192323941093929 | 2149.0 | estimated_donor_stratum | none |
| low_quality | singler | not_flagged | travaglini_2020::travaglini_2020::1 | macro_f1 | 0.8870381314131904 | 3987.0 | estimated_donor_stratum | none |
| low_quality | singler | not_flagged | travaglini_2020::travaglini_2020::2 | macro_f1 | 0.906100582911905 | 3273.0 | estimated_donor_stratum | none |
| low_quality | singler | not_flagged | travaglini_2020::travaglini_2020::3 | macro_f1 | 0.8043132322659298 | 2149.0 | estimated_donor_stratum | none |
| doublet |  | not_assessed |  | macro_f1 |  |  | not_assessed_missing_predeclared_column | none |
| ambient_rna |  | not_assessed |  | macro_f1 |  |  | not_assessed_missing_predeclared_column | none |
| low_rna |  | not_assessed |  | macro_f1 |  |  | not_assessed_missing_predeclared_column | none |
| high_mito |  | not_assessed |  | macro_f1 |  |  | not_assessed_missing_predeclared_column | none |

The locked low-quality rule was `nGene <= 500`. The dataset minimum was 501, so no retained cell entered the flagged stratum; this is evidence of upstream filtering in this asset, not evidence that low-quality robustness was tested.

## Sample-enrichment flags

| method | cluster | n_cells | n_samples | dominant_sample | dominant_sample_fraction | normalized_sample_entropy | flag | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| celltypist | b_cell | 127 | 6 | blood 1 | 0.8740157480314961 | 0.3027284478036828 | SAMPLE_ENRICHED | estimated_descriptive |
| popv | b_cell | 139 | 7 | blood 1 | 0.841726618705036 | 0.34765062125400814 | SAMPLE_ENRICHED | estimated_descriptive |
| singler | b_cell | 126 | 6 | blood 1 | 0.8888888888888888 | 0.28001924150834445 | SAMPLE_ENRICHED | estimated_descriptive |

## QC and retained negative findings

| category | finding | severity | detail |
| --- | --- | --- | --- |
| scope | single_study_single_platform_technical_verification | claim_boundary | No multi-study or cross-platform robustness claim is supported. |
| azimuth | ineligible_donor_or_sample_overlap | retained_negative_result | Official HLCA sample-level metadata includes Stanford/Krasnow samples from Travaglini et al. with donor 1, donor 2, and donor 3, the same three donors used by the evaluation cohort. Platform separation does not restore donor independence. |
| celltypepilot | complete_abstention_on_all_out_of_fold_cells | retained_negative_result | Coverage was 0/9409; this result is retained and was not replaced by candidate labels. |
| batch_sensitivity:platform | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypepilot |
| batch_sensitivity:platform | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypist |
| batch_sensitivity:platform | not_estimable_lt2_supported_levels | diagnostic_unassessed | popv |
| batch_sensitivity:platform | not_estimable_lt2_supported_levels | diagnostic_unassessed | singler |
| batch_sensitivity:plate_batch | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypepilot |
| batch_sensitivity:plate_batch | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypist |
| batch_sensitivity:plate_batch | not_estimable_lt2_supported_levels | diagnostic_unassessed | popv |
| batch_sensitivity:plate_batch | not_estimable_lt2_supported_levels | diagnostic_unassessed | singler |
| batch_sensitivity:condition | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypepilot |
| batch_sensitivity:condition | not_estimable_lt2_supported_levels | diagnostic_unassessed | celltypist |
| batch_sensitivity:condition | not_estimable_lt2_supported_levels | diagnostic_unassessed | popv |
| batch_sensitivity:condition | not_estimable_lt2_supported_levels | diagnostic_unassessed | singler |
| doublet | not_assessed_missing_predeclared_column | diagnostic_unassessed | Missing diagnostics are not evidence of absence. |
| ambient_rna | not_assessed_missing_predeclared_column | diagnostic_unassessed | Missing diagnostics are not evidence of absence. |
| low_rna | not_assessed_missing_predeclared_column | diagnostic_unassessed | Missing diagnostics are not evidence of absence. |
| high_mito | not_assessed_missing_predeclared_column | diagnostic_unassessed | Missing diagnostics are not evidence of absence. |
| sample_enrichment | predicted_label_concentrated_in_one_sample | diagnostic_flag | celltypist:b_cell -> blood 1 (0.874) |
| sample_enrichment | predicted_label_concentrated_in_one_sample | diagnostic_flag | popv:b_cell -> blood 1 (0.842) |
| sample_enrichment | predicted_label_concentrated_in_one_sample | diagnostic_flag | singler:b_cell -> blood 1 (0.889) |
| method_comparison | minimal_inference_only_three_paired_donors | underpowered | Exact sign-flip tests and donor bootstrap intervals are reported, but n=3 donors cannot support a superiority claim. |

## Claim boundary

Cell-level metrics are descriptive. Donors are the independent units. Three donors provide only minimal inferential support; one study and one platform cannot establish study, platform, tumor, inflamed-tissue, or general batch robustness. Missing doublet or ambient-RNA metadata is unassessed, not negative. Confidence values from each tool retain their method-specific, non-probabilistic semantics.
