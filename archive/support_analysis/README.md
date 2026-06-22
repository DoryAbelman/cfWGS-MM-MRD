# Archived Support Analyses

This folder contains scripts that are preserved for auditability but are not
part of the routine command-line manuscript regeneration pipeline.

## `3_3_Plot_optimal_cutoff_tumor_naive_calls_and_clinical_concordance.R`

This script explores tumor-naive blood cfWGS clinical-concordance behavior.
It is archived because no current final manuscript figure, extended-data
figure, or supplementary table is mapped to it in
`docs/manuscript_artifact_source_map.tsv`.

Routine manuscript regeneration uses:

- `3_1_Optimize_cfWGS_thresholds.R` for model training/metric artifacts.
- `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` for dilution-series
  model application and Supplementary Table 7.
- `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` for the final
  tumor-informed clinical-concordance figures and Supplementary Tables 8 and 10.

Keep this script only if tumor-naive support/sensitivity review is needed. It
should not be required for reviewers to regenerate the manuscript figures and
tables.

## `4_3_cfWGS_vs_EasyM_Proteomic_MRD_Comparison.R`

This script explores direct cfWGS-vs-EasyM proteomic MRD comparisons and writes
support figures/tables under `Output_figures_2025/` and
`Output_tables_2025/cfWGS_vs_EasyM_comparison/`.

It is archived because no current final manuscript figure, extended-data
figure, or supplementary table is mapped to it in
`docs/manuscript_artifact_source_map.tsv`. The active manuscript workflow uses:

- `3_1_A_Process_and_optimize_EasyM.R` to process EasyM data, choose/report the
  optimized EasyM calls, and export `EasyM_all_samples_with_optimized_calls.csv`.
- `4_1_Survival_Analysis.R` to merge those EasyM calls into the survival and
  relapse-detection analyses that create the final Figure 3F/Figure 4E and
  Extended Data Figure 6/8 outputs.

Keep this script only for supplemental EasyM review or future manuscript
extensions. It should not be required for reviewers to regenerate the current
manuscript figures and tables.
