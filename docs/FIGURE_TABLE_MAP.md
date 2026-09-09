# Figure and table map

This map identifies the script that calculates or draws each result used in the
paper. It does not include scripts that only copy files into local manuscript
folders.

## Main figures

| Panel | Script | Notes |
| --- | --- | --- |
| Figure 1A | `2_1_Part2_Cohort_Swim_Plot.R` | Treatment timeline; also writes Supplementary Table 1. |
| Figure 1B | `1_6_Identify_High_Quality_Patient_Pairs.R` | Exports the sample-flow counts; the visual panel is assembled manually. |
| Figure 2A | `2_4_Longitudinal_features_analysis.R` | Longitudinal example patients. |
| Figure 2B-E | `2_4B_Build_all_evaluable_longitudinal_panels.R` | All-evaluable longitudinal summaries. |
| Figure 3A | `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | BM-informed 50-repeat patient-grouped nested-CV ROC panel. |
| Figure 3B | `3_1_Optimize_cfWGS_thresholds.R` | Preserved-model training and test-cohort performance display. |
| Figure 3C | `3_1C_Summarize_dilution_correlations_across_patients.R`, using scored inputs from `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Seven series from four patients; technical replicates are averaged within patient, then the four patients are weighted equally. |
| Figure 3D-E | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | BM-informed clinical concordance. |
| Figure 3F | `4_1_Survival_Analysis.R` | One-year-maintenance landmark survival analysis. |
| Figure 4A | `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | Blood-informed 50-repeat patient-grouped nested-CV ROC panel. |
| Figure 4B | `3_1_Optimize_cfWGS_thresholds.R` | Preserved-model training and test-cohort performance display. |
| Figure 4C-D | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Blood-informed clinical concordance. |
| Figure 4E | `4_1_Survival_Analysis.R` | Blood-informed survival analysis. |

## Extended Data figures

| Panel | Script | Notes |
| --- | --- | --- |
| Extended Data Figure 1 | `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Baseline WGS alteration heatmap. |
| Extended Data Figure 2A-C and 2E-F | `2_3_Feature_Concordance_And_Mutation_Counts.R` | Mutation, CNA, FISH, and BM/cfDNA concordance summaries. |
| Extended Data Figure 2D | Manual figure assembly; `5_1_Export_Locked_Figure_Source_Data.R` exports the numerical source | VA-09 chromosome 1 copy-number profile with the 1q FISH-probe interval. |
| Extended Data Figure 2G | `1_2_Part2_Get_Mutation_Overlap.R` | All-evaluable patient-level mutation-set Jaccard overlap. |
| Extended Data Figure 3A-C | `2_4B_Build_all_evaluable_longitudinal_panels.R` | All-evaluable longitudinal companion panels. |
| Extended Data Figure 3D-E | `1_8C_Analyze_MRDetect_Healthy_Control_Platform_Calibration.R`; `1_8D_Build_ED3DE_MRDetect_Platform_Calibration.R` | MRDetect healthy-control platform comparison and matched-reference calibration. |
| Extended Data Figure 3F-G | `1_7D_Compare_CHARM_Xplus_HC_controls.R`; `1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R` | Paired fragmentomics platform shift and leave-one-control-out calibration validation. |
| Extended Data Figure 4 | `2_4_Longitudinal_features_analysis.R` | Longitudinal feature analysis. |
| Extended Data Figure 5A | `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | BM grouped nested-CV operating-point panel. |
| Extended Data Figure 5B-C | `3_1_Optimize_cfWGS_thresholds.R` | BM model performance support panels. |
| Extended Data Figure 5D | `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | BM dilution-series panel. |
| Extended Data Figure 5E-H | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | BM clinical-comparator panels. |
| Extended Data Figure 6A-K | `4_1_Survival_Analysis.R` | BM-informed survival and prospective relapse-detection panels. |
| Extended Data Figure 7A | `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Blood grouped nested-CV operating-point panel. |
| Extended Data Figure 7B-C and 7E | `3_1_Optimize_cfWGS_thresholds.R` | Blood model performance support panels. |
| Extended Data Figure 7D | `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Blood dilution-series panel. |
| Extended Data Figure 7F-I | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Blood clinical-comparator panels. |
| Extended Data Figure 8A-F | `4_1_Survival_Analysis.R` | Blood-informed survival and prospective relapse-detection panels. |
| Extended Data Figure 9A-B | `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Full model-library grouped nested-CV panels. |
| Extended Data Figure 9C-F | `3_1_Optimize_cfWGS_thresholds.R` | Preserved-model performance panels; panel F is calculated from sample-level predictions for 37 samples from 26 patients. |
| Extended Data Figure 10 | External ichorCNA plotting workflow; `4_2_Compare_subclonal_evolution.R` | The repository script provides the high-risk CNA event summary, not the final genome-wide tracks. |

## Main table

| Table | Script | Notes |
| --- | --- | --- |
| Table 1 | `2_1_Clinical_Demographics_Table.R` | Builds the cohort clinical and demographic summary; final DOCX/PDF formatting includes a document-export step. |

## Supplementary tables

| Table | Script | Current status |
| --- | --- | --- |
| Supplementary Table 1 | `2_1_Part2_Cohort_Swim_Plot.R` | Final CSV matches the retained script output. |
| Supplementary Table 2 | `2_2_Baseline_demographics_by_WGS_heatmap_updated.R`; `2_3_Feature_Concordance_And_Mutation_Counts.R` | Six-sheet assembly requires reconciliation before regeneration. |
| Supplementary Table 3 | `2_3_Feature_Concordance_And_Mutation_Counts.R` | Final CSV matches the retained correlation export. |
| Supplementary Table 4 | `6_12_Patient_Grouped_Repeated_Nested_CV.R`; `6_13_Assemble_All_Model_Grouped_CV_Results.R` | Uses the 50-repeat grouped nested-CV result; final CSV applies rounding and removes one descriptive column. |
| Supplementary Table 5 | `3_1_Optimize_cfWGS_thresholds.R` | Full-training refit metrics; distinct from held-out nested-CV performance. |
| Supplementary Table 6 | `3_1C_Expanded_test_clustered_sensitivity.R` | Final workbook removes one presentation-only `resampling_unit` column. |
| Supplementary Table 7 | `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Pooled correlations across 48 scored libraries plus the selected 73-column scored-data export; this weighting differs from Figure 3C. |
| Supplementary Table 8 | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | All four sheets match the reviewed generated values. |
| Supplementary Table 9 | `4_1_Survival_Analysis.R` | Final workbook uses the 16-row-per-sheet prospective result before later comparator filtering. |
| Supplementary Table 10 | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Values match; the historical generated filename incorrectly calls this Table 9. |

The retained Supplementary Table 2 workbook should not be regenerated until its
six sheets have been reconciled with the current `2_2` and `2_3` outputs.

## Figure source-data workbooks

| File | Script | Notes |
| --- | --- | --- |
| `Source_Data_Main_Figures.xlsx` | `5_1_Export_Locked_Figure_Source_Data.R`; `5_2_Build_Figure_Source_Data_Workbooks.R` | The final manuscript copy has 17 sheets. The generated workbook has 18 and additionally includes `Fig1B`, the manually assembled sample-flow panel. |
| `Source_Data_Extended_Data_Figures.xlsx` | `5_1_Export_Locked_Figure_Source_Data.R`; `5_2_Build_Figure_Source_Data_Workbooks.R` | The final manuscript copy has 54 sheets. The generated workbook has 56 and additionally includes `ED10A` and `ED10B`; the final genome-wide tracks come from the separate ichorCNA plotting workflow. |

`5_1` reconstructs some source tables from retained models, thresholds, and
analysis objects when a complete panel CSV is unavailable. It does not refit
models or redraw figures, but those reconstructed tables should be compared
with the final figure components before replacing the reviewed workbooks. The
retained final workbooks and generated workbooks were not changed during this
review.

## Grouped nested-CV settings used in the paper

`6_12_Patient_Grouped_Repeated_Nested_CV.R` uses five outer folds repeated 50
times, five inner folds repeated five times, and 2,000 patient-clustered
bootstrap replicates. Patients remain grouped in every split. Inner
out-of-fold predictions select the threshold; outer-held-out patients estimate
performance. `6_13` assembles the 32-model result, and `6_17` and the
`6_14_Generate_*` script draw the panels listed above.
