# Script guide

This guide explains how the analysis scripts connect to one another and to the
paper. It uses the current filenames; no files have been renamed for this
release.

Run scripts from the main project directory so their relative input and output
paths resolve correctly.

## Workflow overview

```text
clinical data ─┐
mutations ─────┤
copy number ───┼─> integrated WGS features ─┐
translocations ┘                            │
fragmentomics ──────────────────────────────┼─> integrated analysis table
MRDetect ───────────────────────────────────┤             │
MFC / clonoSEQ / EasyM ─────────────────────┘             ├─> figures and tables
                                                          ├─> model scoring
prespecified model definitions ─> grouped nested CV ──────┘
```

The 50-repeat grouped nested-CV workflow is separate from the historical model
scoring script. `3_1_Optimize_cfWGS_thresholds.R` loads preserved models for
test-cohort scoring by default. `6_12` and `6_13` perform and assemble the
patient-grouped repeated nested CV used in the paper.

### Shared setup and helper files

- `helpers.R` contains study-specific metadata corrections, identity aliases,
  cohort rules, MRDetect parsing/control definitions, revision mutation counts,
  dilution metadata, and FISH recovery used by multiple numbered scripts. Some
  helper calls also write audit CSVs under `Output_tables_2025/`.
- `next_event_endpoint_helpers.R` constructs the sample-relative progression or
  censor endpoint used by `2_4`, `2_4B`, `4_1`, and `4_1B`. It admits a
  progression up to 30 days before a sample as a day-zero event and otherwise
  selects the first later progression or latest available censor date.
- `config.R` supplies the package inventory used by
  `run_pipeline.R --check-packages`. Its three path variables are historical;
  they do not redirect the numbered scripts.
- `setup_packages.R` is an optional common-package loader, not an installer or
  the complete environment definition.
- `publication_export_helpers.R` relabels internal cohort terminology as
  Training or Testing in copied table and workbook data. It does not change
  numeric values, sample inclusion, or model predictions.

## 1. Clinical data and cohort definitions

### `1_0_Process_clinical_metadata.R`

- **Purpose:** harmonize clinical and sample metadata across SPORE, M4,
  IMMAGINE, and MyC and construct relapse, PFS, censoring, and follow-up tables.
- **Main inputs:** cohort clinical spreadsheets, sample identifiers, treatment
  and relapse records, and the reviewed cohort-assignment table.
- **Main outputs:** `combined_clinical_data_updated_April2025.csv`, cleaned M4
  laboratory and relapse tables, and PFS/follow-up objects under
  `Exported_data_tables_clinical/`.
- **Downstream use:** almost every feature-processing, longitudinal, model, and
  survival script.
- **Important detail:** the file also contains older support calculations. The
  current endpoint exports are produced near the end of the script.

### `1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R`

- **Purpose:** process transplant dates, clinical progression, OS follow-up,
  and cohort-specific FISH information.
- **Main inputs:** SPORE, M4, and IMMAGINE clinical workbooks plus outputs from
  `1_0`.
- **Main outputs:** active FISH and clinical helper files under `Clinical data/`,
  including `Clinical data/SPORE/tidy_fish.csv`.
- **Downstream use:** `1_1B`, the integrated analysis table, and Supplementary
  Table 2.
- **Important detail:** the old file `SPORE_fish_flags.csv` is a stale legacy
  export and is not used by the current workflow.

### `1_1B_Process_clinical_labs.R`

- **Purpose:** combine laboratory, staging, demographic, and FISH fields across
  cohorts.
- **Main inputs:** M4, SPORE, and IMMAGINE clinical workbooks and the FISH files
  produced by `1_1A`.
- **Main output:**
  `Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv`.
- **Downstream use:** `2_0_Assemble_Table_With_All_Features.R`.
- **Important detail:** date-based joins use a ±60-day window.

### `1_6_Identify_High_Quality_Patient_Pairs.R`

- **Purpose:** define WGS eligibility and cohort assignments and export the
  counts used for Figure 1B.
- **Main inputs:** integrated WGS features, processing logs, clinical metadata,
  and reviewed cohort assignments.
- **Main outputs:** baseline and longitudinal high-quality patient lists,
  cohort-assignment files, and the Figure 1B source/count tables.
- **Paper use:** Figure 1B source data.
- **Important detail:** despite its filename, the script does not construct
  matched sample pairs.

## 2. WGS alterations

### `1_2_Process_Mutation_Data.R`

- **Purpose:** read BM and cfDNA MAF files, calculate VAFs, annotate sample
  metadata, and create mutation intermediates.
- **Main inputs:** BM and cfDNA MAF files and the clinical table from `1_0`.
- **Main outputs:** BM and blood mutation RDS files plus narrower temporary MAF
  files used by WGS integration.
- **Downstream use:** mutation overlap, Extended Data Figure 1, WGS feature
  integration, and MRDetect-related analyses.
- **Important detail:** the broader RDS outputs and narrower temporary MAFs have
  different eligibility scopes and are not interchangeable.

### `1_2_Part2_Get_Mutation_Overlap.R`

- **Purpose:** calculate patient-level Jaccard overlap between BM and cfDNA
  mutation sets.
- **Main inputs:** the BM and blood mutation RDS files, clinical metadata,
  cohort assignments, and the manuscript display-ID map.
- **Main outputs:** all-evaluable and Frontline overlap tables and lollipop
  plots.
- **Paper use:** the 41-patient all-evaluable plot is Extended Data Figure 2G.

### `1_3_Process_Ig_Translocation_Info.R`

- **Purpose:** parse IgCaller structural-variant calls, annotate cytobands,
  prepare IGV review material, and export manually confirmed Ig translocations.
- **Main inputs:** IgCaller filtered TSV files or the preserved parsed cache,
  cytobands, blacklist regions, and the reviewed IGV workbook.
- **Main output:**
  `Jan2025_exported_data/translocation_data_cytoband_updated.rds` and its text
  companion.
- **Downstream use:** WGS integration, Extended Data Figure 1, and
  Supplementary Table 2.
- **Important detail:** only rows with `Looks_real == 1` are promoted to the
  final positive-call matrix.

### `1_4_Process_CNA_Data.R`

- **Purpose:** convert ichorCNA segments into arm-level del1p, amp1q, del13q,
  del17p, hyperdiploidy, and FISH-probe features.
- **Main inputs:** ichorCNA segment files or the retained combined segment
  cache, cytobands, and probe locations.
- **Main outputs:** `cna_data_ichorCNA.rds` and
  `FISH_probe_calls_bin_cytoband_ichorCNA.rds` under
  `Jan2025_exported_data/`.
- **Downstream use:** WGS integration, baseline heatmaps, FISH concordance, and
  Supplementary Table 2.
- **Important detail:** arm calls use the fraction of evaluated segments, not
  base-pair coverage.

### `1_4A_Process_sequenza_CNA_Data.R`

- **Purpose:** build ploidy-aware Sequenza CNA calls, FISH-probe calls, purity,
  and ploidy estimates.
- **Main inputs:** Sequenza segment and confints files, clinical metadata,
  cytobands, and probe locations.
- **Main outputs:** `cna_data_from_sequenza_400_updated.rds`,
  `FISH_data_from_sequenza_400_updated.rds`, and
  `Sample_ploidy_from_sequenza_400.rds`.
- **Downstream use:** WGS integration, FISH/WGS concordance, and ploidy
  summaries.
- **Important detail:** the direct-coordinate FISH output is the version read
  by `1_5`; the later cytoband-expanded calculation is a comparison.

### `1_5_Integrate_WGS_Feature_Data.R`

- **Purpose:** merge mutation, ichorCNA, Sequenza, Ig translocation,
  FISH-probe, tumour-fraction, and clinical sample features.
- **Main inputs:** outputs from `1_2`, `1_3`, `1_4`, and `1_4A`, tumour-fraction
  data, clinical metadata, and the reviewed IGV workbook.
- **Main outputs:** `All_feature_data_Sep2025_updated2.rds`,
  `CNA_translocation_Sep2025_updated2.rds`, and mutation/CNA helper tables under
  `Jan2025_exported_data/`.
- **Downstream use:** `1_6`, `2_0`, baseline figures, concordance, models, and
  longitudinal analyses.
- **Current limitation:** the presently staged inputs expose an IMG-159 sample
  key mismatch, so a clean top-to-bottom rerun must be checked against the
  retained integrated table before replacing it.

## 3. Fragmentomics and platform calibration

### `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R`

- **Purpose:** calculate per-site nucleosome coverage, midpoint, amplitude,
  z-scores, and MM-DAR summaries.
- **Main inputs:** historical and XPlus nucleosome-distance tables for patients
  and healthy controls.
- **Main outputs:** `MM_DARs_chromatin_activation_data.csv/.rds` and per-site
  metric/statistic tables under `Results_Fragmentomics/`.
- **Downstream use:** `1_7B`, platform calibration, and integrated models.
- **Important detail:** fold change is an ordinary ratio, not log2 fold change.
  Continuous coverage features are the model inputs.

### `1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R`

- **Purpose:** combine fragment score, short-fragment proportion, and MM-DAR
  coverage and harmonize XPlus continuous features to the historical scale.
- **Main inputs:** fragment-score and insert-size summaries, `1_7A` outputs,
  healthy controls, and clinical sample metadata.
- **Main outputs:** `Key_fragmentomics_data_updated2.csv/.rds` and the platform
  reference-parameter table.
- **Downstream use:** `2_0`, patient scoring, and fragmentomics figures.
- **Important detail:** the transformation is estimated from 19 matched control
  identities. It changes continuous XPlus predictors before scoring and is not
  only a plotting adjustment.

### `1_7D_Compare_CHARM_Xplus_HC_controls.R`

- **Purpose:** prepare the 19 matched controls and describe platform shifts
  across scalar and high-dimensional fragmentomics domains.
- **Main inputs:** historical and XPlus healthy-control fragmentomics outputs.
- **Main output used downstream:**
  `Results_Fragmentomics/CHARM_Xplus_HC_comparison/sample_level_feature_values.csv`.
- **Downstream use:** `1_7E` and `1_8E`.
- **Important detail:** its exploratory Wilcoxon tests are unpaired even though
  the controls are identity matched. The transformation itself is defined in
  `1_7B`.

### `1_7E_Harmonize_Fragmentomics_Sequencing_Platforms.R`

- **Purpose:** plot the fixed mean/SD transformation before and after mapping
  the same 19 XPlus controls to the historical scale.
- **Inputs:** the `1_7D` matched-control table, `1_7B` parameter table, and
  `1_7A` regulatory-coverage values.
- **Outputs:** optional QC tables and plots under
  `Results_Fragmentomics/Healthy_control_platform_harmonization/`.
- **Important detail:** this is an in-sample implementation check, not the
  leave-one-control-out validation used in the paper.

### `1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R`

- **Purpose:** validate the fragmentomics mapping while leaving each control
  identity out of the parameter estimation.
- **Inputs:** matched-control values and the historical feature definitions.
- **Outputs:** panel data and plots for the platform-calibration analysis.
- **Paper use:** Extended Data Figure 3F-G.

## 4. MRDetect

### `1_8_Process_Cumulative_VAFs_MRDetect.R`

- **Purpose:** process patient MRDetect output, calculate detection-rate and
  mutation-panel-specific healthy-control z-score features, and attach sample
  metadata.
- **Main inputs:** Winter 2025 and Spring 2026 MRDetect CSV exports, combined
  clinical metadata, the integrated baseline feature table, cohort assignments,
  and the optional audited VA15 rerun replacement table.
- **Control reference:** historical samples use the historical CHARM controls;
  Spring 2026 XPlus samples request the allowlisted XPlus CHARM controls.
- **Main outputs:** complete raw/z-scored tables, BM-informed and blood-informed
  longitudinal CSVs, the historical CHARM control export, mutation-source
  inventories, and audit files under
  `MRDetect_output_winter_2025/Processed_R_outputs/` and
  `Output_tables_2025/clinical_support/`.
- **Downstream use:** `2_0` adds the BM- and blood-informed features to the
  integrated analysis table; `2_4` uses the complete z-scored data; `2_1_Part2`
  uses the all-patient blood-informed table; and `1_10`/`3_1_part2` use the
  historical control export.
- **Important boundary:** the retained `sites_rate_zscore_charm > 4.5` screen is
  an upstream feature definition. Final model-based cfWGS calls are defined in
  `3_1`, not in this script.

### `1_8A_Process_Cumulative_VAFs_for_dilution_series.R`

- **Purpose:** process the experimental dilution-series MRDetect results and
  calculate the z-score features used by the downstream dilution analysis.
- **Control references:** historical NovaSeq 6000 dilution samples use 26
  CHARM libraries. Spring 2026 NovaSeq XPlus dilution samples use the complete
  allowlisted 22-library XPlus reference, representing 21 control identities.
  The 19 controls matched one-to-one across platforms are reserved for the
  platform-comparison and reference-sensitivity analyses.
- **Main outputs:** the dilution-series raw and z-scored RDS/text tables plus
  reference-completeness, mutation-list, and input-availability audits under
  `MRDetect_output_winter_2025/Processed_R_outputs/`.
- **Downstream use:** `3_1_part2` uses the z-scored table for Figure 3C,
  Extended Data Figures 5D and 7D, and Supplementary Table 7.

### `1_8C_Analyze_MRDetect_Healthy_Control_Platform_Calibration.R`

- **Purpose:** compare historical and XPlus healthy-control MRDetect values
  using the 19 shared control identities as the inferential unit.
- **Important detail:** the eight mutation/VCF panels are repeated assay
  contexts, not 152 independent controls.
- **Downstream use:** `1_8D`.

### `1_8D_Build_ED3DE_MRDetect_Platform_Calibration.R`

- **Purpose:** plot the MRDetect platform shift and platform-matched reference
  calibration.
- **Paper use:** Extended Data Figure 3D-E.

`1_8B_Export_MRDetect_patient_feature_table.R` is a derived export not read by a
main manuscript script. Keep it outside the required sequence unless that table
is specifically needed.

## 5. Integrated tables and descriptive figures

### `2_0_Assemble_Table_With_All_Features.R`

- **Purpose:** assemble the sample-level table containing clinical data,
  genomic features, fragmentomics, MRDetect, MFC, clonoSEQ, and EasyM fields.
- **Main inputs:** the processed clinical, WGS, fragmentomics, and MRDetect
  tables from stages 1-4.
- **Main output:**
  `Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds`.
- **Downstream use:** all final descriptive, model, concordance, longitudinal,
  and survival analyses.
- **Important detail:** the script currently retains several historical
  fallback values to reproduce the submitted table; do not simplify those
  joins without row-level comparisons.

### Manuscript output scripts

| Script | Manuscript use |
| --- | --- |
| `2_1_Clinical_Demographics_Table.R` | Main Table 1 |
| `2_1_Part2_Cohort_Swim_Plot.R` | Figure 1A and Supplementary Table 1 |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Extended Data Figure 1 and Supplementary Table 2 sheets A-B source data |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Extended Data Figure 2A-C and 2E-F, Supplementary Table 2 sheets C-F source data, and Supplementary Table 3 |
| `2_4_Longitudinal_features_analysis.R` | Figure 2A and Extended Data Figure 4 |
| `2_4B_Build_all_evaluable_longitudinal_panels.R` | Figure 2B-E and Extended Data Figure 3A-C |

Supplementary Table 2 combines outputs from `2_2` and `2_3`. Its final six-sheet
assembly is not currently implemented in one script. Do not replace the
retained final workbook without reconciling all six sheets against the reviewed
`2_2` and `2_3` outputs.

Extended Data Figure 2D is a manually assembled VA-09 chromosome 1 copy-number
panel. `5_1_Export_Locked_Figure_Source_Data.R` exports its numerical Sequenza
segment and FISH-probe source data; `2_3` does not draw that panel.

## 6. Model scoring, assay comparisons, and survival

### `3_1_Optimize_cfWGS_thresholds.R`

- **Default role:** load the preserved February 2026 model library and
  thresholds and score the current cohort.
- **Outputs:** scored BM-informed and blood-informed call tables, performance
  panels including the sample-level Extended Data Figure 9F confusion matrix,
  and the full-training refit metrics used in Supplementary Table 5.
- **Important detail:** the older model-development branch behind
  `CFWGS_RETRAIN_MODELS=1` is not the final 50-repeat grouped nested CV.

### `3_1_A_Process_and_optimize_EasyM.R`

- **Required role:** calculate the prespecified isotype-specific EasyM calls:
  negative when residual IgG is at most 1% of the baseline BM immunoglobulin or
  residual IgA/light-chain signal is at most 0.05% of the baseline BM
  immunoglobulin.
- **Downstream use:** `3_2` and `4_1`.
- **Important detail:** the same script also contains exploratory optimized
  thresholds. Those are separate from the prespecified calls used in the paper.

### `3_1C_Expanded_test_clustered_sensitivity.R`

- **Purpose:** calculate patient-clustered bootstrap and earliest-sample-per-
  patient sensitivity summaries without retraining the models.
- **Paper use:** Supplementary Table 6.

### `3_2_Plot_optimal_cutoff_and_clinical_concordance.R`

- **Purpose:** build clinical assay concordance plots, contingency tables,
  calibration support, and decision-curve support.
- **Paper use:** Figures 3D-E and 4C-D, Extended Data Figures 5E-H and 7F-I,
  and Supplementary Tables 8 and 10.

### `4_1_Survival_Analysis.R`

- **Purpose:** build landmark survival analyses and prospective relapse-
  detection summaries.
- **Inputs:** scored cfWGS and clinical-assay calls, PFS event/censor tables,
  patient follow-up, and progression dates.
- **Paper use:** Figures 3F and 4E, Extended Data Figures 6 and 8, and
  Supplementary Table 9.

`4_1B_Build_all_evaluable_first_nonbaseline_KM.R` contains alternative survival
anchors and is not the generator of the submitted Figure 3F.

## 7. Dilution series

| Script | Role |
| --- | --- |
| `1_7C_Process_fragmentomics_data_dilution_series_updated.R` | Process dilution-series fragmentomics. |
| `1_8A_Process_Cumulative_VAFs_for_dilution_series.R` | Process dilution-series MRDetect output. |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Apply preserved models; build Extended Data Figures 5D and 7D and Supplementary Table 7; write the point-level inputs for Figure 3C. |
| `3_1C_Summarize_dilution_correlations_across_patients.R` | Build Figure 3C from seven series across four patients, with technical replicates averaged within patient and patients weighted equally. |

Supplementary Table 7 treats the 48 dilution libraries as observations. Figure
3C instead shows the seven series-level correlations and an equal-patient mean
across four patients. Other three-patient and patient-mean displays in
`3_1C_Summarize_dilution_correlations_across_patients.R` are sensitivity
analyses.

`1_9_Create_dilution_series_eligibility_table.R` and
`1_10_Estimate_MRDetect_LOD_proxy.R` are experiment-planning/QC scripts and are
not required to reproduce the submitted dilution results. The manuscript
analysis begins with the completed physical dilution measurements processed by
`1_7C` and `1_8A`; the values calculated by `1_9` are candidate mixing proxies,
not measured tumour fractions or a formal analytical LOD.

## 8. Patient-grouped repeated nested cross-validation

### `6_12_Patient_Grouped_Repeated_Nested_CV.R`

- **Purpose:** estimate performance for pre-specified model definitions while
  keeping all samples from each patient together in every outer and inner fold.
- **Design used in the paper:** five outer folds repeated 50 times, five inner
  folds repeated five times, and 2,000 patient-clustered bootstrap replicates.
- **Models:** one-predictor specifications use centered/scaled logistic
  regression; multi-predictor specifications tune elastic-net alpha and lambda
  within the grouped inner resamples.
- **Threshold selection:** inner out-of-fold predictions only.
- **Evaluation:** outer-held-out patients only.
- **Model library used in the paper:** `all`, which is also the current script
  default. Passing `--model-library all` explicitly records that choice in the
  run command.
- **Outputs:** versioned run directories containing folds, predictions, tuning,
  thresholds, metrics, bootstrap results, warnings, checksums, and QC files.

### `6_13_Assemble_All_Model_Grouped_CV_Results.R`

- **Purpose:** combine the completed BM, blood, and full-cohort fragmentomics
  blocks and verify the expected 32-model inventory.
- **Paper use:** assembled values for Supplementary Table 4 and all final
  grouped-CV panels.

### Grouped-CV figure scripts

| Script | Paper use |
| --- | --- |
| `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | Figures 3A and 4A |
| `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Extended Data Figures 5A, 7A, and 9A-B |

The tracked `6_14_Build_Corrected_Grouped_CV_Panels.R` produces an earlier
four-model display and is not the final manuscript panel generator.

`6_13_Temporal_Validation_Subset_Audit.R` separately compares fixed-model
performance in the original seven-patient hold-out and later-accrued patients.
It is a sensitivity audit, not an input to a manuscript figure or table.

## 9. CNA change analysis

### `4_2_Compare_subclonal_evolution.R`

- **Purpose:** compare longitudinal del1p, amp1q, del13q, and del17p calls and
  summarize emergent events.
- **Paper use:** event summaries supporting Extended Data Figure 10.
- **Boundary:** the final genome-wide ichorCNA tracks were produced by an
  external plotting workflow.

Patient-specific mutation-tracking scripts `4_2B` through `4_2H` are not used in
the paper and are not part of the GitHub release.

## 10. Figure source-data workbooks

### `5_1_Export_Locked_Figure_Source_Data.R`

- **Purpose:** prepare 74 panel-level CSVs covering Figures 1-4 and Extended
  Data Figures 1-10.
- **Inputs:** final figure PDFs, available panel source CSVs, retained
  figure/model inputs, and `id_map.rds`.
- **Outputs:** panel CSVs, a manifest, audit table, and schema contract under
  `Output_tables_2025/Figure_Source_Data/`.
- **Important detail:** despite its historical name, this script reconstructs
  some panel tables when a complete panel CSV is unavailable. It does not
  redraw figures or refit models, but the reconstructed values require
  comparison with the final figure components before release.

### `5_2_Build_Figure_Source_Data_Workbooks.R`

- **Purpose:** assemble the `5_1` panel CSVs into the main-figure and Extended
  Data Excel workbooks.
- **Representation changes:** suitable text columns are converted to logical or
  numeric values, calendar-date columns are removed, and cohort labels are
  changed to Training or Testing.
- **Safety:** the script checks sheet counts and original patient identifiers,
  but overwrites the workbook destinations.

Neither source-data script is called by `run_pipeline.R`.

The workbooks currently stored with the final manuscript contain 17 main-
figure sheets and 54 Extended Data sheets. The generated 18/56-sheet versions
add `Fig1B`, `ED10A`, and `ED10B`. No other sheet names differ. Figure 1B is
assembled from exported counts, and the final Extended Data Figure 10 tracks
come from the separate ichorCNA plotting workflow. The retained final workbooks
were not replaced during this review.

## 11. Outputs and table assembly

Final figures and tables are stored separately from the scripts in the working
project.

At present:

- Supplementary Tables 1, 3, 5, and 8 reconcile with their reviewed script
  outputs;
- Supplementary Tables 4, 6, 7, 9, and 10 need explicit final formatting or
  filename/staging steps; and
- Supplementary Table 2 requires a reviewed scientific reconciliation before
  it is regenerated.

The local manuscript-copy, package-building, Code Ocean, and audit scripts are
not analysis entry points and are omitted from this guide.

Two tracked maintenance helpers also sit outside the analysis workflow:

- `run_manuscript_workflow.R` combines the incomplete configured source plan
  with local manuscript-output refresh and packaging steps; it is not the
  complete paper workflow.
- `shorten_manuscript_figure_filenames.R` renames files in place without a dry
  run and is not required for manuscript reproduction.
