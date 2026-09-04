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

Keep this script only if tumor-naive support/sensitivity analysis is needed. It
is not required to regenerate the manuscript figures and tables. If run
manually, it writes only to `Output_figures_2025/tumor_naive_support/`
and `Output_tables_2025/tumor_naive_support/`, not to
`final_manuscript_objects/`.

The script applies a legacy fixed blood-probability threshold of 0.457 and also
derives a separate high-specificity threshold. These are exploratory support
rules and are not the current mapped manuscript calling rules.

Older duplicate copies of the pre-rerun `3_3` code are treated as local-only
legacy material and are ignored by Git. The copy in this `support_analysis/`
folder is the documented optional support version.

## `4_3_cfWGS_vs_EasyM_Proteomic_MRD_Comparison.R`

This script explores direct cfWGS-vs-EasyM proteomic MRD comparisons and writes
support figures under `Output_figures_2025/`. It creates
`Output_tables_2025/cfWGS_vs_EasyM_comparison/` for a local cache but does not
write analysis tables or source-data files there.

The EasyM binary status in this archived script comes from the supplied legacy
positive/negative workbook export. It does not calculate the manuscript's
isotype-specific reference-threshold call from residual immunoglobulin as a
percentage of baseline bone-marrow immunoglobulin. That calculation is in
`3_1_A_Process_and_optimize_EasyM.R`.

It is archived because no current final manuscript figure, extended-data
figure, or supplementary table is mapped to it in
`docs/manuscript_artifact_source_map.tsv`. The active manuscript workflow uses:

- `3_1_A_Process_and_optimize_EasyM.R` to process EasyM data, calculate the
  isotype-specific reference-threshold call, and export
  `EasyM_all_samples_with_optimized_calls.csv`.
- `4_1_Survival_Analysis.R` to merge those EasyM calls into the survival and
  relapse-detection analyses that create the final Figure 3F/Figure 4E and
  Extended Data Figure 6/8 outputs.

Keep this script only for legacy EasyM exploration. It is not required to
regenerate the current manuscript figures and tables.
