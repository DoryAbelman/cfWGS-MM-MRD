# cfWGS-MM-MRD

Analysis code for *"Cell-free DNA Whole Genome Sequencing for Non-Invasive
Minimal Residual Disease Detection in Multiple Myeloma"*.

The repository contains the R scripts used to process the clinical and cfWGS
data, construct the analysis tables, fit and evaluate the MRD classifiers, and
produce the main and supplementary figures and tables. The scripts retain their
original analysis-based numbering so that saved objects and manuscript records
can be traced back to the code that generated them.

The files are not one simple linear pipeline: several scripts branch into
mutation, copy-number, fragmentomics, MRDetect, dilution-series, longitudinal,
and survival analyses. The tables below describe the purpose, inputs, outputs,
and manuscript contribution of each script. The
[`docs/SCRIPT_GUIDE.md`](docs/SCRIPT_GUIDE.md) explains how the scripts connect,
and [`docs/FIGURE_TABLE_MAP.md`](docs/FIGURE_TABLE_MAP.md) provides a concise
panel-level crosswalk.

The analysis should currently be run as the documented individual stages. The
tracked `run_pipeline.R` is an incomplete internal helper: it requires local
metadata that are not included in this repository, ends at `4_2`,
and does not include the final 50-repeat nested-CV or figure-source-workbook
steps. It is therefore not the public entry point for reproducing the paper.

**Raw data files are not included in this repository** (see [Data availability](#data-availability) below).

---

## Quick start

Run commands from the repository root. Before running a script, stage the
protected inputs listed in its header and confirm that the preceding stage has
produced the required objects. Use [Pipeline scripts](#pipeline-scripts) for
the source-data and figure stages and
[Patient-grouped repeated nested cross-validation](#patient-grouped-repeated-nested-cross-validation)
for the analysis used in Figures 3A and 4A, Extended Data Figures 5A, 7A and
9A-B, and Supplementary Table 4.

Scripts that call `manuscript_output_helpers.R` also place manuscript-labeled
copies of selected outputs under:

```text
final_manuscript_objects/
```

That folder is organized by final manuscript item, for example:

```text
final_manuscript_objects/
├── 01_main_figures/Figure_1/Figure_1A/
├── 02_extended_data_figures/
├── 03_main_tables/
├── 04_supplementary_tables/Supplementary_Table_1/
└── manuscript_direct_output_manifest.tsv
```

Each copied/saved file is recorded in
`final_manuscript_objects/manuscript_direct_output_manifest.tsv` with the
artifact ID, final figure/table label, source path, destination path, checksum,
and generating script.

The local manuscript-number helper writes section-organized review files under:

```text
manuscript_writing/
```

The key files are:

```text
manuscript_numbers_by_section.tsv
manuscript_numbers_by_section.xlsx
manuscript_numbers_by_section.md
manuscript_draft_paragraph_index.tsv
manuscript_numeric_paragraph_review.tsv
manuscript_paragraph_metric_audit.tsv
by_section/
```

`manuscript_numbers_by_section.*` lists
the section, statistic label, formatted value, numerator/denominator when
relevant, source file/script, related figure or table, update trigger, and
caveats.

The paragraph review files are internal manuscript-checking aids: they extract
number-containing paragraphs from local DOCX drafts and link them to analysis
outputs where possible. They are not required to run the scientific analyses.

The following internal maintenance command refreshes and then validates the
organized manuscript-output tree:

```sh
Rscript validate_manuscript_outputs.R
```

It requires the local manuscript artifact map, which is not included here. It
is also not read-only: it first refreshes mapped copies and
indexes under `final_manuscript_objects/`, and can therefore replace files in
that generated output tree. It then checks the source map, direct-output
manifest, organized paths, script headers, helper calls, and README coverage,
and writes:

```text
final_manuscript_objects/manuscript_output_validation_report.tsv
```

After the refresh, these indexes summarize what each mapped script produces:

```text
final_manuscript_objects/manuscript_output_index.tsv
final_manuscript_objects/script_output_index.tsv
```

## Manuscript figure/table ownership

The numbered scripts are still organized by analysis question, not by final
figure number. Each script header lists the same outputs shown below. When a
script runs, the files used in the manuscript are copied or saved into
`final_manuscript_objects/` under the final figure/table label.

| Script | Final manuscript outputs |
| --- | --- |
| `1_2_Part2_Get_Mutation_Overlap.R` | Extended Data Figure 2G |
| `1_8C_Analyze_MRDetect_Healthy_Control_Platform_Calibration.R`; `1_8D_Build_ED3DE_MRDetect_Platform_Calibration.R` | Extended Data Figure 3D-E (MRDetect healthy-control platform shift and platform-matched z-score calibration) |
| `1_7D_Compare_CHARM_Xplus_HC_controls.R`; `1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R` | Extended Data Figure 3F-G (paired fragmentomics platform shift and leave-one-control-out calibration validation) |
| `1_6_Identify_High_Quality_Patient_Pairs.R` | Figure 1B source table |
| `2_1_Clinical_Demographics_Table.R` | Table 1 |
| `2_1_Part2_Cohort_Swim_Plot.R` | Figure 1A; Supplementary Table 1 |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Extended Data Figure 1; Supplementary Table 2 sheets A-B |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Extended Data Figure 2A-C and 2E-F; Supplementary Table 2 sheets C-F; Supplementary Table 3 |
| `2_4_Longitudinal_features_analysis.R` | Figure 2A; Extended Data Figure 4 |
| `2_4B_Build_all_evaluable_longitudinal_panels.R` | Figure 2B-E; Extended Data Figure 3A-C |
| `3_1_Optimize_cfWGS_thresholds.R` | Figure 3B; Figure 4B; Extended Data Figures 5B-C, 7B-C/7E and 9C-F; Supplementary Table 5; saved full-training models and thresholds used downstream |
| `3_1C_Expanded_test_clustered_sensitivity.R` | Patient-clustered bootstrap and deterministic one-sample-per-patient sensitivity outputs for the revision-inclusive Figure 3B/Figure 4B test analyses; final Supplementary Table 6 workbook combining those results with the preserved expanded-test classifier metrics and exact sample manifest |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Scored dilution inputs for Figure 3C; Extended Data Figure 5D; Extended Data Figure 7D; Supplementary Table 7 |
| `3_1C_Summarize_dilution_correlations_across_patients.R` | Figure 3C patient-weighted correlation panel |
| `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Figure 3D-E; Figure 4C-D; Extended Data Figure 5E-H; Extended Data Figure 7F-I; Supplementary Tables 8 and 10 |
| `4_1_Survival_Analysis.R` | Figure 3F; Figure 4E; Extended Data Figures 6A-K and 8A-F; Supplementary Table 9. The time-window results use prospective labels that require progression on or after the sample date, or adequate follow-up for a non-event call. |
| `4_2_Compare_subclonal_evolution.R` | Supporting event summary for Extended Data Figure 10; the final genome-wide tracks come from the separate ichorCNA plotting workflow |
| `6_12_Patient_Grouped_Repeated_Nested_CV.R` | Final 50-repeat patient-grouped nested-CV fitting and held-out predictions used for Supplementary Table 4 |
| `6_13_Assemble_All_Model_Grouped_CV_Results.R` | Combines the 32 model-specific grouped-CV runs and writes the final performance table |
| `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | Figure 3A and Figure 4A |
| `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Extended Data Figures 5A, 7A, 9A and 9B |

`run_manuscript_workflow.R` is an internal orchestration helper. It runs the
configured source plan and then performs several local refresh, validation,
manuscript-number, and table-packaging steps. It does not run the final
50-repeat grouped-CV scripts or the complete figure-source workbook workflow,
and it currently runs
`3_1C_Expanded_test_clustered_sensitivity.R` twice during a full execution.
It is therefore not the recommended entry point for reproducing the paper.

When used in the full internal project, `run_pipeline.R` writes timestamped
logs and `run_manifest.tsv` under `pipeline_logs/`. It skips
`3_1_Optimize_cfWGS_thresholds.R` by default because that script can refit the
historical model library and update saved training-derived objects. The final
50-repeat validation is a separate, explicitly versioned analysis described
below.

Each numbered script header lists the manuscript figure/table outputs that the
script creates or updates. Where a final manuscript object is generated, the
script calls `manuscript_output_helpers.R` to place a labeled copy into
`final_manuscript_objects/`. Use
[`docs/FIGURE_TABLE_MAP.md`](docs/FIGURE_TABLE_MAP.md) for the current,
manually reviewed panel-to-script map.

---

## Pipeline scripts

The stage numbers describe the broad order of the analysis. Scripts within a
stage may form separate branches, and optional/support scripts are identified
below. Run commands from the repository root.

### Stage 0 - Cohort selection

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `0_1_select_additional_samples_for_cohort_expansion.R` | Identify possible test-cohort expansion samples by linking IMMAGINE inventory records to patient IDs, classifying samples as Baseline / MRD / Treatment, and reporting blood-MRD opportunities. This is a planning analysis and is not required to reproduce the current paper. | IMMAGINE patient ID mapping, sample inventory, MRD date, diagnosis date, and treatment-line Excel files (`Clinical data/IMMAGINE/`). | Cohort-expansion planning tables under `Output_tables_2025/cohort_expansion/`; matching, duplicate, and blood-MRD checks under `Output_tables_2025/cohort_expansion/support_qc/`. |

### Stage 1 - Data processing

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `1_0_Process_clinical_metadata.R` | Harmonise clinical metadata across SPORE, M4, IMMAGINE, and MyC cohorts; annotate timepoints; compute relapse/PFS date helpers; and export the main clinical metadata table consumed by downstream mutation, CNA, fragmentomics, swim-plot, model, and survival scripts. Support-only count summaries, BAM review lists, and QC plots are separated from manuscript outputs. | Clinical metadata spreadsheets from `Clinical data/`, `M4_CMRG_Data/`, and reviewed `cohort_assignment_table_updated.rds`. | Downstream inputs retained at historical paths, especially `combined_clinical_data_updated_April2025.csv`, `M4_labs_cleaned.csv`, `Relapse_dates_M4_clean.csv`, `Relapse dates cfWGS updated.csv`, and PFS/censor-date files under `Exported_data_tables_clinical/`; support-only summaries under `Output_tables_2025/clinical_metadata_support/`. |
| `1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R` | Parse transplant dates, relapse/progression events, OS follow-up, and FISH flags for SPORE, M4, and IMMAGINE. This is an upstream clinical-helper step, not a direct manuscript figure/table owner. | SPORE, M4, and IMMAGINE clinical Excel/CSV workbooks with ASCT, progression, FISH, and OS follow-up data; clinical and relapse-date outputs from `1_0`. | Active downstream clinical inputs in `Clinical data/SPORE/tidy_*.csv` and `Clinical data/Exported clinical data April 2025/*.csv`; support-only relapse/OS helper tables and QC plots under `Output_tables_2025/clinical_followup_support/`. |
| `1_1B_Process_clinical_labs.R` | Compile cross-cohort clinical laboratory, staging, demographic, and FISH/cytogenetic helper fields for M4, SPORE, and IMMAGINE. This is an upstream clinical-helper step, not a direct manuscript figure/table owner. | M4 clinical workbooks (`M4_COHORT_*.xlsx`), curated SPORE/IMMAGINE clinical sheets, and FISH helper files generated by `1_1A`. | Active downstream input `Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv`, consumed by `2_0_Assemble_Table_With_All_Features.R`. |
| `1_2_Process_Mutation_Data.R` | Ingest bone-marrow and cfDNA MAF files, compute VAFs, annotate with clinical metadata, and create the mutation-call intermediates used by downstream mutation-overlap, heatmap, feature-integration, and MRDetect analyses. VAF/read-depth plots are support-only QC, not manuscript panels. | BM and cfDNA WGS MAF files (`*.maf`), `combined_clinical_data_updated_April2025.csv`. | Core intermediates `combined_maf_bm_all_muts.rds`, `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts_updated.rds`, and temporary MAF exports used by downstream scripts; support-only QC plots under `Final Tables and Figures/mutation_processing_support/`. |
| `1_2_Part2_Get_Mutation_Overlap.R` | Compute patient-level Jaccard overlap between baseline/diagnosis bone-marrow and cfDNA mutation sets. The 41-patient all-evaluable manuscript-ID lollipop panel is used for Extended Data Figure 2G; the 28-patient Frontline version is a cohort-specific support analysis. Support-only Venn diagrams and QA barplots are written under `Final Tables and Figures/mutation_overlap_support/`. | `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts_updated.rds`, `combined_clinical_data_updated_April2025.csv`, `cohort_assignment_table_updated.rds`, `id_map.rds` | Extended Data Figure 2G all-evaluable panel PNG plus source-data and summary CSVs; cohort-specific and per-patient support plots are retained separately. |
| `1_3_Process_Ig_Translocation_Info.R` | Parse Ig-caller structural variant calls, map breakpoints to cytobands, flag MM-relevant translocations, generate IGV review helpers, and export the manually IGV-confirmed translocation feature matrix. This is an upstream feature-generation step, not a direct manuscript figure/table owner. | Ig-caller `*_filtered.tsv` files when available, otherwise the preserved parsed Ig-caller cache; cytoband and ENCODE blacklist references; `Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm` with `Looks_real` review calls. | Active downstream inputs `Jan2025_exported_data/translocation_data_cytoband_updated.rds/txt`; support-only parsing caches, pre-IGV matrices, and QC plots under `Output_tables_2025/translocation_processing_support/` and `Output_figures_2025/translocation_processing_support/`. |
| `1_4_Process_CNA_Data.R` | Active upstream ichorCNA processing step. Summarises raw ichorCNA segment calls into arm-level del1p, amp1q, del13q, del17p, and hyperdiploidy features; also creates ichorCNA calls at clinical FISH probe loci. It does not directly own a final figure/table, but its outputs feed WGS feature integration, baseline feature heatmaps, FISH/WGS concordance, and Supplementary Table 2 support. | Raw ichorCNA `*.seg` files in `Oct 2024 data/Ichor_CNA/`; if those are unavailable in a lightweight review bundle, the preserved combined segment cache under `Output_tables_2025/ichor_cna_processing_support/` is used. | Active downstream inputs `Jan2025_exported_data/cna_data_ichorCNA.rds/txt` and `Jan2025_exported_data/FISH_probe_calls_bin_cytoband_ichorCNA.rds/txt`; support/QC segment cache and hyperdiploidy audit tables under `Output_tables_2025/ichor_cna_processing_support/`. |
| `1_4A_Process_sequenza_CNA_Data.R` | Active upstream Sequenza CNA processing step. Converts raw Sequenza segments and confints files into ploidy-aware arm-level CNA calls, direct-coordinate FISH-probe helper calls, and sample-level purity/ploidy estimates. It does not directly own a final figure/table, but its outputs feed WGS feature integration, FISH/WGS concordance, ploidy summaries, and Supplementary Table 2 support. | Raw Sequenza `*_segments.txt(.gz)` files in `Oct 2024 data/Sequenza/All_Segments_400/`, matching `*_confints_CP.txt` files in `Oct 2024 data/Sequenza/All_confints_400/`, clinical metadata, cytobands, and FISH probe locations. | Active downstream inputs `Jan2025_exported_data/cna_data_from_sequenza_400_updated.rds/txt`, `Jan2025_exported_data/FISH_data_from_sequenza_400_updated.rds/txt`, and `Jan2025_exported_data/Sample_ploidy_from_sequenza_400.rds`; support/QC segment, probe, cytoband-expanded FISH, ploidy, and summary caches under `Output_tables_2025/sequenza_cna_processing_support/`. |
| `1_5_Integrate_WGS_Feature_Data.R` | Active upstream WGS feature-integration step. Merges mutation, ichorCNA/Sequenza CNA, IG translocation, FISH-probe CNA, tumor-fraction, clinical sample metadata, and reviewed IGV translocation overrides into the final integrated WGS feature tables. It does not directly own a final figure/table, but it is the main bridge into aggregate clinical table assembly, baseline heatmaps, concordance analyses, model training/application, and longitudinal analyses. | Current outputs from `1_2`, `1_3`, `1_4`, and `1_4A`; `Oct 2024 data/tumor_fraction_cfWGS.txt`; reviewed IGV translocation spreadsheet; clinical metadata. | Active downstream inputs `Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds/txt`, `Jan2025_exported_data/CNA_translocation_Sep2025_updated2.rds/txt`, `Jan2025_exported_data/CNA_at_FISH_sites_combined.rds/txt`, `Jan2025_exported_data/mutation_export_updated2.rds/txt`, and `Jan2025_exported_data/mutation_export_updated_more_info2.rds/txt`; support/QC checks under `Output_tables_2025/feature_integration_support/`. |
| `1_6_Identify_High_Quality_Patient_Pairs.R` | Define high-quality baseline sample eligibility, cohort assignment, and the Figure 1B sample-flow source table. The final Figure 1B visual panel remains manually assembled from this exported table, but the underlying counts/QC annotations are regenerated from the command line. | `All_feature_data_Sep2025_updated2.rds`, processing log XLSX, clinical metadata, sample-flow workbook, reviewed `cohort_assignment_table_updated.rds`. | Figure 1B source table in `Final Tables and Figures/` and `final_manuscript_objects/`; downstream intermediates `Output_tables_2025/high_quality_patients_list_for_baseline_mut_calling2.*` and `patient_cohort_assignment.*`; support-only processing log under `Output_tables_2025/sample_qc_support/`. |
| `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R` | Compute nucleosome-accessibility coverage, midpoint, amplitude, and platform-specific healthy-control z-scores. Fold change is an ordinary tumour/healthy ratio: values >1 use the upper reference boundary and values ≤1 use the lower boundary. The resulting `Threshold.*` fields are descriptive binary flags; downstream models use continuous coverage features. | Historical and XPlus patient/control nucleosome-distance TSV files, clinical metadata, and the shared helpers. | Per-site metrics/statistics, `MM_DARs_chromatin_activation_data.csv/rds`, platform reference thresholds, and a corrected-versus-legacy flag audit under `Results_Fragmentomics/`. |
| `1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R` | Compute fragment-size (short-fragment) scores and merge with nucleosome accessibility and clinical data. | `insert_size_summary.tsv`, `fragment_scores.tsv`, `MM_DARs_chromatin_activation_data.csv`. | `Key_fragmentomics_data_updated.csv` |
| `1_7C_Process_fragmentomics_data_dilution_series_updated.R` | Apply the same nucleosome-accessibility and fragment-size processing to the physical dilution-series samples. It uses the same ordinary-ratio direction boundary at 1 and records the historical `>0` flags only for sensitivity comparison. | Dilution-series nucleosome-distance files, fragment-size summaries, dilution metadata, and platform-specific healthy controls. | `key_fragmentomics_info_dilution_series.csv/rds`, per-site metrics/statistics, platform thresholds, and a corrected-versus-legacy flag audit used upstream of the dilution model application. |
| `1_7D_Compare_CHARM_Xplus_HC_controls.R` | Assemble the 19 paired healthy-control profiles used to compare fragmentomics features across NovaSeq 6000 and NovaSeq XPlus. | Matched-control fragment scores, insert-size summaries, nucleosome-accessibility files, and control identity map. | `sample_level_feature_values.csv`, comparison summaries, and QC plots under `Results_Fragmentomics/CHARM_Xplus_HC_comparison/`. |
| `1_7E_Harmonize_Fragmentomics_Sequencing_Platforms.R` | Reconstruct and plot the fixed location/scale mapping used to place XPlus fragmentomics features on the NovaSeq 6000 reference scale. This is an in-sample implementation check and does not retrain a model. | Paired-control table from `1_7D` and the platform mean/SD parameters used by `1_7B`. | Platform summaries and before/after harmonization figures under `Results_Fragmentomics/Healthy_control_platform_harmonization/`. |
| `1_8_Process_Cumulative_VAFs_MRDetect.R` | Process MRDetect cfDNA mutation-detection output; compute z-scores relative to CHARM healthy controls; filter to patient timepoints. | MRDetect CSV output files (`MRDetect_outputs/*.csv`). | `cfWGS_Winter2025All_MRDetect_with_Zscore.rds/txt`, BM and Blood processed CSVs. |
| `1_8A_Process_Cumulative_VAFs_for_dilution_series.R` | Process dilution-series MRDetect output and calculate platform- and mutation-panel-specific healthy-control z-scores. Historical NovaSeq 6000 samples use the preserved 26-library reference; Spring 2026 XPlus samples use the complete allowlisted 22-library reference. The 19 one-to-one paired controls are used only in the separate platform comparison and sensitivity analysis. | Historical and Spring 2026 dilution MRDetect CSV files, dilution metadata, clinical metadata, and the preserved historical healthy-control RDS. | `cfWGS_Winter2025Dilution_series_May2025_with_zscore.rds/txt` plus availability, mutation-list, and healthy-reference audit CSVs used upstream of Figure 3C, Extended Data Figures 5D and 7D, and Supplementary Table 7. |
| `1_8B_Export_MRDetect_patient_feature_table.R` | Create an analysis-ready matched-patient MRDetect table containing both BM-derived and cfDNA-derived mutation-source features. This is a derived export and does not modify the processed MRDetect input. | Processed MRDetect z-score RDS from `1_8` and the cfWGS sample identity map. | `MRDetect_patient_features_BM_and_cfDNA_baselines_Feb2026.csv/rds` under `MRDetect_output_winter_2025/Processed_R_outputs/Derived_exports/`. |
| `1_8C_Analyze_MRDetect_Healthy_Control_Platform_Calibration.R` | Compare paired healthy-control MRDetect results across sequencing platforms and evaluate platform-matched z-score calibration at the level of the 19 matched controls. | Processed CHARM and XPlus healthy-control MRDetect outputs and the control identity map. | Paired-control summaries and source tables used by Extended Data Figure 3D-E. |
| `1_8D_Build_ED3DE_MRDetect_Platform_Calibration.R` | Plot the MRDetect platform comparison and matched-reference calibration analysis prepared by `1_8C`. | The paired-control and calibration tables written by `1_8C`. | Extended Data Figure 3D-E panels and plotted source-data files. |
| `1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R` | Validate fragmentomics platform calibration by leaving out each of the 19 matched controls in turn, estimating the mapping from the other 18, and evaluating the held-out control. | Paired fragmentomics control table from `1_7D`. | Extended Data Figure 3F-G panels, paired-test summaries, and plotted source data. |
| `1_9_Create_dilution_series_eligibility_table.R` | Optional experimental-planning analysis that selects a high-signal (≥0.5% detection rate) and low-signal (≤0.05% and z-score <2) physical sample from the same patient/mutation-panel group and calculates candidate mixing fractions for nominal MRDetect detection-rate targets from 10⁻¹ to 10⁻⁶. These targets are planning proxies, not measured tumour fractions or a formal analytical LOD. The manuscript analysis starts from the completed dilution-series measurements instead. | `cfWGS_MRDetect_BM_data_updated_Feb2026.csv` from `1_8`. | Candidate-pair, mixing-plan, and feasibility CSVs; none are read by the manuscript figure or table scripts. |
| `1_10_Estimate_MRDetect_LOD_proxy.R` | Add read-denominator-based one-read and 95% Poisson LOD proxy fields for MRDetect QC and interpretation. This does not rerun MRDetect or change detection calls. | Baseline/control MRDetect detection-rate table and its informative-read denominator (`reads_checked` by default). | `All_detection_rates_baseline_and_controls_Feb2026_with_mrdetect_lod.csv`. |

The platform-calibration branches have explicit dependencies that are not fully
represented in the current runner: run `1_7D` before `1_7E` or `1_8E`, and run
`1_8C` before `1_8D`. The `1_7E` plot is an optional implementation check;
`1_8D` and `1_8E` generate the platform-calibration panels used in Extended
Data Figure 3.

### Stage 2 - Summary statistics and baseline figures

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `2_0_Assemble_Table_With_All_Features.R` | Merge MRD assay results (MFC, clonoSEQ, EasyM), cfWGS metrics, clinical data, and fragmentomics into a single master table. | Outputs from all Stage 1 scripts. | `Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated*.rds` |
| `2_1_Clinical_Demographics_Table.R` | Build Table 1 (patient demographics, disease characteristics by cohort). | Master feature table from `2_0`. | `table1_categorical_updated_final.docx`, `cohort_assignment_table.rds` |
| `2_1_Part2_Cohort_Swim_Plot.R` | Generate the treatment-timeline swim plot and manuscript event table using the approved display IDs and days from baseline. Support-only cohort/legend QA exports are written under `Final Tables and Figures/swim_plot_support/`. | `tidy_treatments.csv`, M4 and IMMAGINE chemotherapy tables. | Figure 1A component; Supplementary Table 1 event table. Support QA exports are not copied to `final_manuscript_objects/`. |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Create the baseline integrated alteration heatmap (BM overlaid with cfDNA; mutations, CNAs, translocations) and disease-feature catalog. | `Final_aggregate_table*.rds`, cohort assignment, CNA/translocation/mutation RDS files. | Extended Data Figure 1 component; feature catalogue used as sheet A of Supplementary Table 2. |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Compute FISH-WGS concordance; summarise baseline mutation counts by cohort; build concordance and feature-correlation exports. BAM archive/unarchive helper tables are support-only and run only when `CFWGS_RUN_BAM_ARCHIVE_DIAGNOSTICS=true`. | Master feature table, cohort assignments, mutation export RDS. Optional BAM diagnostic mode also reads `All_bam_storage_locations.xlsx`. | Extended Data Figure 2A-C and 2E-F; Supplementary Tables 2 and 3. Panel 2D is assembled separately, with its numerical source exported by `5_1`. Optional BAM diagnostic files are written under `Output_tables_2025_updated/support_only_bam_archive_diagnostics/` and are not copied to `final_manuscript_objects/`. |
| `2_4_Longitudinal_features_analysis.R` | Generate the selected longitudinal examples and the broader all-evaluable feature-change analysis. Earlier frontline-only panels retained in this file have been replaced by `2_4B`. | Master feature table and cohort assignments. | Figure 2A and Extended Data Figure 4; additional historical/support panels. |
| `2_4B_Build_all_evaluable_longitudinal_panels.R` | Build the all-evaluable longitudinal summaries used for the current multi-patient panels. | Master feature table, event dates, and cohort assignments. | Figure 2B-E and Extended Data Figure 3A-C, with their source tables. |

### Stage 3 - Model training, model application, and validation

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `3_1_Optimize_cfWGS_thresholds.R` | Fit or load the historical elastic-net model library, full-training refits, and thresholds used for downstream sample scoring. Its older nested-CV exports are not the final 50-repeat patient-grouped validation; that analysis is performed by `6_12` and assembled by `6_13`. | Master feature table, clonoSEQ/MFC ground-truth labels. | Preserved models and thresholds; Figures 3B and 4B; Extended Data Figures 5B-C, 7B-C/7E, and 9C-F; full-training metrics used in Supplementary Table 5; inputs to downstream scoring and concordance analyses. |
| `3_1_A_Process_and_optimize_EasyM.R` | Process EasyM residual immunoglobulin as a percentage of the baseline BM immunoglobulin and calculate the prespecified isotype-specific call: negative at ≤1% for IgG and ≤0.05% for IgA/light-chain disease. The file also retains exploratory optimized thresholds separately. | EasyM CSV files, clinical/isotype information, and the cfWGS call table from `3_1`. | `EasyM_all_samples_with_optimized_calls.csv`, including the named isotype-specific reference call, plus threshold and QC tables used by `3_2` and `4_1`. |
| `3_1C_Expanded_test_clustered_sensitivity.R` | Evaluate the saved classifier calls in the expanded test cohort while accounting for repeated samples by patient-clustered bootstrap and one-sample-per-patient sensitivity analyses. It does not refit the classifiers or change their thresholds. | Current scored test-cohort table and preserved model thresholds. | Supplementary Table 6, exact scored-sample manifest, clustered-bootstrap summaries, and one-sample-per-patient sensitivity results supporting Figure 3B and Figure 4B. |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Apply saved models and thresholds to the experimental dilution-series libraries and describe detection performance across the tested dilution range. This is not a formal analytical limit-of-detection study. | Saved models/thresholds from `3_1`; fragmentomics and MRDetect dilution-series outputs. | Scored patient/replicate inputs for Figure 3C; Extended Data Figures 5D and 7D; pooled 48-library correlations and scored data for Supplementary Table 7. |
| `3_1C_Summarize_dilution_correlations_across_patients.R` | Calculate dilution-series Spearman correlations by patient and technical replicate, average replicates within patient, and then give each of the four patients equal weight. | Point-level dilution source tables written by `3_1_part2`. | Figure 3C and its patient-, series-, and equal-patient source data; alternate three-patient and patient-mean displays. |
| `3_1D_Build_Dilution_cVAF_vs_controls.R` | Build additional cVAF and MRDetect z-score plots comparing four dilution series with matched healthy-control measurements and unrelated plasma BAM/mutation-panel pairings. | Dilution and healthy-control source tables, all-by-all MRDetect object, high-quality mutation-panel list, and sample-scoring manifest. | Additional descriptive plots and source CSVs; not used in a final manuscript panel. |
| `3_1D_Audit_MRDetect_XPlus_reference_update.R` | Compare sample scores before and after the XPlus MRDetect reference update and report changed values and calls. | Two scored RDS files supplied as command-line arguments. | Sample-level impact CSV and call-flip summary; audit only. |
| `3_1E_Audit_MRDetect_XPlus_dilution_reference_sensitivity.R` | Compare dilution scores calculated with the paired 19-control and primary 22-library XPlus references. | Two dilution CSV files supplied as command-line arguments. | Sample-level impact, call-flip, and correlation-summary CSVs; audit only. |
| `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Generate tumour-informed cfWGS clinical-concordance figures: assay positivity, model-vs-clinical assay comparisons, calibration/decision-curve support, and contingency tables. | `all_patients_with_BM_and_blood_calls_updated*.rds`, threshold table. | Figure 3D-E and Figure 4C-D; Extended Data Figures 5E-H and 7F-I; Supplementary Tables 8 and 10. |
| `archive/support_analysis/3_3_Plot_optimal_cutoff_tumor_naive_calls_and_clinical_concordance.R` | Archived tumour-naive blood cfDNA support/sensitivity analysis. This is not required for routine manuscript regeneration because the final clinical-concordance figures/tables are produced by `3_2`, and no current final manuscript artifact is mapped to `3_3`. | `all_patients_with_BM_and_blood_calls_updated*.rds`, threshold table. | Support-only tumour-naive review outputs; not copied to `final_manuscript_objects/`. |

Figure 3D and Figure 4C use the same cohort assignment but different frontline
eligibility rules. Figure 3D includes all 42 landmark samples with an evaluable
BM-informed cfWGS call; 39 of those also have MFC or clonoSEQ data. Figure 4C
first requires an evaluable blood-informed cfWGS call and at least one of MFC
or clonoSEQ, leaving 41 of 46 blood-call-evaluable landmark samples.

### Stage 4 - Clinical outcome analyses

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `4_1_Survival_Analysis.R` | Kaplan-Meier PFS curves stratified by MRD status at landmark timepoints (Post-ASCT, 1yr Maintenance). Calculate sensitivity of each assay for detecting future relapse. The current manuscript time-window outputs are regenerated under `detection_progression_updated6` from prospective labels that require either a future progression event or adequate follow-up through each prediction window. PFS event/censor dates remain separate from patient-level last-follow-up dates; the prospective labels use `patient_followup_dates_updated.*` generated by `1_0_Process_clinical_metadata.R` when available. | `all_patients_with_BM_and_blood_calls_updated*.rds`, EasyM calls, `Censor_dates_per_patient_for_PFS_updated.rds`, `patient_followup_dates_updated.rds`, curated progression-date CSVs. | Figure 3F and Figure 4E components; Extended Data Figures 6 and 8; Supplementary Table 9; main and prospective QC outputs under `Output_tables_2025/detection_progression_updated6/`. |
| `4_1B_Build_all_evaluable_first_nonbaseline_KM.R` | Optional descriptive all-evaluable Kaplan–Meier companion analysis using one prospective non-baseline assessment per patient. It is not a submitted manuscript panel. | Current cfWGS/clinical call table, PFS and relapse-date objects, follow-up tables, and cohort-specific treatment records. | Versioned KM figures, source data, denominator audits, and Cox-model summaries under `final_manuscript_objects/additional_all_evaluable_first_nonbaseline_km/`. |
| `4_2_Compare_subclonal_evolution.R` | Identify emergent CNA events between baseline and relapse cfDNA samples. The final Extended Data Figure 10 genome-wide CNA tracks come from an external ichorCNA plotting workflow; this script provides repo-side supporting event outputs. | `All_feature_data_Sep2025_updated2.rds`, cohort assignments. | `Emergent_CNA_events.csv` and subclonal-evolution support outputs for Extended Data Figure 10. |
| `archive/support_analysis/4_3_cfWGS_vs_EasyM_Proteomic_MRD_Comparison.R` | Archived cfWGS-vs-EasyM support/sensitivity analysis. This is not required for routine manuscript regeneration because final EasyM call generation is handled by `3_1_A`, and final survival/relapse EasyM panels are handled by `4_1`. | cfWGS call table, EasyM quantitative and binary CSVs, PFS censor dates. | Support-only EasyM comparison figures/tables; not copied to `final_manuscript_objects/`. |

### Stage 5 - Patient-grouped model validation and final ROC panels

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `6_12_Patient_Grouped_Repeated_Nested_CV.R` | Fit and evaluate the complete model library using patient-grouped repeated nested cross-validation. The manuscript run uses 50 outer repeats, five outer folds, five inner repeats, five inner folds, and 2,000 patient-clustered bootstrap replicates. | Preserved training data and model/validation objects, current patient/sample identities, and cohort assignments. | Versioned fold assignments, held-out predictions, tuning results, thresholds, repeat-level metrics, clustered-bootstrap summaries, QC files, and `RUN_COMPLETE`. |
| `6_13_Assemble_All_Model_Grouped_CV_Results.R` | Combine the separately run BM, blood, and full-cohort fragmentomics blocks and verify that all 32 model/cohort specifications are present exactly once. It does not fit models. | Three completed `6_12` run directories. | Combined 32-model result used by the final plots and Supplementary Table 4. |
| `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | Plot mean empirical ROC curves from the 50 outer repeats for the complete BM-informed and blood-informed model sets. It does not refit models or select new thresholds. | Completed combined result from `6_13`. | Figure 3A and Figure 4A PNG/PDF panels and plotted source-data tables. |
| `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Plot the grouped nested-CV operating-point and fragmentomics validation panels from the completed result. | Completed combined result from `6_13`. | Extended Data Figures 5A, 7A and 9A-B, with plotted ROC and operating-point source tables. |

`6_13_Temporal_Validation_Subset_Audit.R` is a separate sensitivity analysis.
It compares fixed-model performance in the original seven-patient hold-out and
the patients accrued later. It does not refit a model or change a threshold,
and its results are not used in a manuscript figure or table.

### Stage 6 - Figure source-data workbooks

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `5_1_Export_Locked_Figure_Source_Data.R` | Prepare the 74 panel-level CSVs for Figures 1-4 and Extended Data Figures 1-10. Despite its historical name, it reconstructs some panel tables from retained models, thresholds, and analysis objects when a complete panel CSV is unavailable; it does not redraw figures or refit models. | Final figure PDFs, generated panel CSVs, retained figure/model inputs, and `id_map.rds`. | Panel CSVs, a 74-row manifest, audit table, and schema contract under `Output_tables_2025/Figure_Source_Data/`. |
| `5_2_Build_Figure_Source_Data_Workbooks.R` | Assemble the panel CSVs from `5_1` into Excel workbooks, convert suitable columns to logical or numeric values, remove calendar-date columns, and check for original patient identifiers. | The `5_1` manifest, audit table, panel CSVs, and `id_map.rds`. | Main-figure and Extended Data source-data workbooks. Existing workbook destinations are overwritten. |

These scripts are not called by `run_pipeline.R`. The panel reconstructions in
`5_1` should be compared with the final figures before regenerated workbooks
replace the reviewed manuscript copies. The currently retained final workbooks
contain 17 main-figure sheets and 54 Extended Data sheets. The generated
18/56-sheet workbooks additionally contain `Fig1B`, `ED10A`, and `ED10B`, which
are absent from the retained final copies. Figure 1B is assembled from exported
counts, and the final Extended Data Figure 10 tracks come from the separate
ichorCNA plotting workflow.

---

## Helper files

| File | Role |
|------|------|
| `setup_packages.R` | Optionally checks and loads a commonly used subset of packages. Individual scripts still load their own dependencies, and `config.R` contains the broader package inventory. |
| `config.R` | Package inventory used by `run_pipeline.R --check-packages`, plus three historical path constants. The numbered scripts use their own project-relative paths rather than these constants. |
| `helpers.R` | Shared loaders and study-specific rules for metadata corrections, sample identities and aliases, cohort assignment, MRDetect parsing and controls, revision mutation counts, dilution metadata, and baseline FISH calls. Some functions write audit CSVs under `Output_tables_2025/`. |
| `next_event_endpoint_helpers.R` | Shared sample-relative progression/censor endpoint construction used by the longitudinal and survival analyses. It applies the 30-day event grace rule and two documented relapse-sample day-zero overrides by default. |
| `publication_export_helpers.R` | Relabels internal cohort terms as Training or Testing in copied table/workbook data. It does not change numeric measurements, sample inclusion, model probabilities, or the caller's input object in place. |
| `prepare_local_inputs.R` | Local author utility that checks a limited input manifest and can recursively copy missing paths from another local mirror. It is not a complete clean-clone input manifest; even `--check` writes a report. |
| `manuscript_output_helpers.R` | Shared utilities used by numbered scripts to copy/save selected figure, table, and source-data components under `final_manuscript_objects/`. Some functions replace mapped copies or remove stale aliases, so they should not be called as read-only inspection tools. |
| `validate_manuscript_outputs.R` | Refreshes mapped manuscript-output copies and indexes, then runs structural checks and writes `manuscript_output_validation_report.tsv`. A passing report confirms mapped-file organization, not scientific equivalence to the manuscript. |
| `run_manuscript_workflow.R` | Internal wrapper around the incomplete configured source plan and local manuscript-output checks. It is not the complete paper workflow and is not the recommended entry point. |
| `5_0_Build_Manuscript_Text_Number_Exports.R` | Local manuscript-checking helper that indexes working DOCX drafts and writes section-organized number tables under `manuscript_writing/`. It is not part of the scientific analysis. |
| `pipeline_metadata.R` | Reads `docs/manuscript_artifact_source_map.tsv` and builds script-to-artifact crosswalks used by the command-line runner and documentation. |
| `shorten_manuscript_figure_filenames.R` | Historical in-place filename migrator for `final_manuscript_objects/`. It has no dry-run and can invalidate links; it is not required to reproduce the paper and should not be run routinely. |

---

## Directory structure

```
Final_Scripts/
├── 0_1_*.R                 ← optional cohort-expansion planning analysis
├── 1_0_*.R … 4_2_*.R      ← clinical, feature, figure, model, and outcome scripts
├── 6_12_*.R … 6_17_*.R    ← 50-repeat grouped nested-CV and plotting scripts
├── archive/support_analysis/ ← analyses not required for the submitted paper
├── config.R                ← package inventory and historical path constants
├── helpers.R               ← shared utility functions
├── manuscript_output_helpers.R ← direct manuscript-output helper for numbered scripts
├── validate_manuscript_outputs.R ← refresh-and-validate output helper
├── pipeline_metadata.R      ← shared script/output metadata helpers
├── final_manuscript_objects/ ← manuscript-labeled outputs created at run time
├── docs/                    ← script guide and figure/table crosswalk
├── run_pipeline.R           ← incomplete internal source-plan runner
├── setup_packages.R        ← optional common-package checker/loader
└── README.md               ← this file
```

---

## Running the analysis

**Prerequisites**

- R ≥ 4.2
- Key libraries include **tidyverse**, **readxl**, **data.table**, **lubridate**, **ComplexHeatmap**, **circlize**, **ChromHeatMap**, **maftools**, **GenomicRanges**, **gtsummary**, **officer**, **flextable**, **gt**, **pROC**, **patchwork**, **rmda**, **exact2x2**, **survival**, **survminer**, **tableone**, **ggbreak**, **ggridges**, **ggpubr**, **GGally**, **pbapply**, **PRROC**, **VennDiagram**, **caret**, **glmnet**, **Matrix**, **DescTools**, **GeneCycle**, **RColorBrewer**, **fuzzyjoin**, **rstatix**, **openxlsx**, **scales**, **viridis**, and **writexl**. See each script header and `config.R` for dependencies.

**Source-script workflow**

1. Obtain the protected inputs described in [Data availability](#data-availability)
   and place them at the project-relative paths listed in each script header.
   `prepare_local_inputs.R` may help authors copy an existing local mirror, but
   it does not yet cover every input required for a clean-clone reproduction.
2. Run from the repository root. Most scripts use their own project-relative paths;
   the historical constants in `config.R` do not redirect all inputs or outputs.
3. Follow the stage table above and the input/output header in each script.
   Several branches are independent, so only run the stages needed for the
   intended figure or table.
4. Run the repeated nested-CV workflow below separately when reproducing its
   figures and Supplementary Table 4.

### Patient-grouped repeated nested cross-validation

The manuscript's internal-validation analysis assigns patients—not individual
sample rows—to both outer and inner folds. Thus, all repeated samples belonging
to a patient remain on the same side of every split. This workflow is separate
from the saved full-training models and thresholds used to score the test
cohort; it does not overwrite those test-scoring objects.

#### What each script does and why

| Step | Script | What it does | Why it is separate |
|---|---|---|---|
| 1 | `6_12_Patient_Grouped_Repeated_Nested_CV.R` | Reconstructs the historical training frame, creates repeated outer and inner folds grouped by patient, tunes `glmnet` inside each outer-training set, derives a Youden threshold from inner out-of-fold predictions, and scores outer-held-out patients. It then reports repeat-pooled metrics with patient-clustered bootstrap intervals. | This is the only step that fits validation models. Keeping it versioned protects the frozen test-scoring models and makes every fold, prediction, threshold, and warning auditable. |
| 2 | `6_13_Assemble_All_Model_Grouped_CV_Results.R` | Combines the separately executed BM, blood, and full-cohort fragmentomics blocks; verifies that all 32 model/cohort combinations are present once and checks prediction and fold integrity. | Splitting the computation into three blocks makes the long run manageable; this assembly step does not fit models. |
| 3 | `6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R` | Calculates and plots mean empirical ROC curves across the 50 outer repeats. | This creates the main Figure 3A and Figure 4A panels without refitting a model or selecting a new threshold. |
| 4 | `6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R` | Summarizes outer-test-fold operating points and generates the remaining grouped-CV panels. | This creates Extended Data Figures 5A, 7A and 9A-B and their plotted source tables from the same completed result. Supplementary Table 4 comes from the assembled `6_13` performance table. |

#### Manuscript settings

The paper uses the complete 32-model library, five outer folds repeated 50
times, five inner folds repeated five times, 2,000 patient-clustered bootstrap
replicates, and base seed `20260731`. These are now the script defaults. The
commands below still specify them explicitly so the saved command records the
complete resampling design.

Run the three model blocks from the repository root with new, unused run IDs:

```sh
Rscript 6_12_Patient_Grouped_Repeated_Nested_CV.R \
  --run-id paper_50rep_bm_v1 --model-library all \
  --models BM_Sites,BM_cVAF,BM_Raw_cVAF,BM_All_Mutation_Features,BM_Combined_Mutation_Zscores,BM_Mutation_Fragmentomics_Full,BM_Mutation_Fragmentomics_Min,Fragmentomics_BMRestricted_Full,Fragmentomics_BMRestricted_Min,Fragmentomics_BMRestricted_FS,Fragmentomics_BMRestricted_Mean_Coverage,Fragmentomics_BMRestricted_Proportion_Short,Fragmentomics_BMRestricted_Tumor_Fraction \
  --outer-repeats 50 --inner-repeats 5 --bootstrap-reps 2000 \
  --seed 20260731

Rscript 6_12_Patient_Grouped_Repeated_Nested_CV.R \
  --run-id paper_50rep_blood_v1 --model-library all \
  --models Blood_Sites,Blood_cVAF,Blood_Raw_cVAF,Blood_All_Mutation_Features,Blood_Combined_Mutation_Zscores,Blood_Mutation_Fragmentomics_Full,Blood_Mutation_Fragmentomics_Min,Fragmentomics_BloodRestricted_Full,Fragmentomics_BloodRestricted_Min,Fragmentomics_BloodRestricted_FS,Fragmentomics_BloodRestricted_Mean_Coverage,Fragmentomics_BloodRestricted_Proportion_Short,Fragmentomics_BloodRestricted_Tumor_Fraction \
  --outer-repeats 50 --inner-repeats 5 --bootstrap-reps 2000 \
  --seed 20260731

Rscript 6_12_Patient_Grouped_Repeated_Nested_CV.R \
  --run-id paper_50rep_fullfrag_v1 --model-library all \
  --models Fragmentomics_FullCohort_Full,Fragmentomics_FullCohort_Min,Fragmentomics_FullCohort_FS,Fragmentomics_FullCohort_Mean_Coverage,Fragmentomics_FullCohort_Proportion_Short,Fragmentomics_FullCohort_Tumor_Fraction \
  --outer-repeats 50 --inner-repeats 5 --bootstrap-reps 2000 \
  --seed 20260731
```

Each run writes all fold assignments, held-out predictions, fitted settings,
thresholds, warnings, summaries, session information, and QC results under
`Output_tables_2025/patient_grouped_repeated_nested_cv/<run-id>/`. It refuses
to overwrite an existing run directory and writes `RUN_COMPLETE` only after
its checks pass. The preserved model/validation input files listed in the
`6_12` header must be present and non-writable.

After all three blocks complete, assemble the 32-model result:

```sh
Rscript 6_13_Assemble_All_Model_Grouped_CV_Results.R \
  --source-runs paper_50rep_bm_v1,paper_50rep_blood_v1,paper_50rep_fullfrag_v1 \
  --output-run-id paper_50rep_combined_v1
```

Generate the main ROC panels and the Extended Data operating-point panels from
that completed combined result:

```sh
Rscript 6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R \
  --input-run-id paper_50rep_combined_v1 \
  --output-run-id paper_50rep_main_roc_v1

Rscript 6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R \
  --input-run-id=paper_50rep_combined_v1 \
  --output-run-id=paper_50rep_extended_data_v1
```

The `6_14_Generate_*` script uses `--name=value` arguments; the preceding
`6_12`, `6_13`, and `6_17` commands use `--name value` arguments.

#### Interpretation boundaries

- The grouped-CV AUC estimates internal discrimination for a previously unseen
  patient, conditional on the already-selected feature specification.
- Inner folds tune `glmnet` and generate out-of-fold predictions for threshold
  selection; outer-held-out patients are used only for evaluation.
- Patient-clustered bootstrap intervals preserve all repeated samples from a
  resampled patient.
- The model specifications were selected before this repeated validation; the
  grouped-CV AUC therefore evaluates those fixed specifications rather than a
  new feature-discovery process.
- None of these scripts overwrites the saved full-training models, thresholds,
  older CV objects, or test-cohort predictions.

Outputs are written to new versioned directories under
`Output_tables_2025/patient_grouped_repeated_nested_cv/` and
`Output_figures_2025/patient_grouped_repeated_nested_cv/`.

---

## Data availability

Because the input files contain protected patient information they are **not distributed with this repository**. This includes:

- Clinical metadata spreadsheets (SPORE, M4, IMMAGINE cohorts)
- WGS MAF files, ichorCNA segmentation outputs, and Ig-caller structural variant calls
- MRDetect cumulative-VAF CSV files
- Fragmentomics nucleosome-distance and fragment-score files
- EasyM proteomic MRD data

To reproduce the analysis, request the de-identified data from the study authors.

---

## Output

The pipeline produces:

- **CSV / XLSX tables** - cleaned metadata, mutation counts, MRD feature matrices, concordance statistics, sensitivity/specificity tables (written to `Output_tables_2025/`)
- **PDF / PNG figures** - the manuscript panels generated in R (primarily under `Final Tables and Figures/` and `Output_figures_2025/`); Figure 1B is assembled from an exported source table, and the final Extended Data Figure 10 genome-wide tracks come from a separate ichorCNA plotting workflow
- **Support-only QA outputs** - non-manuscript review exports are kept in explicitly named support folders such as `Final Tables and Figures/swim_plot_support/` and `Output_tables_2025_updated/support_only_bam_archive_diagnostics/`
- **RDS objects** - intermediate R data objects passed between scripts

---

## Contact

For questions about the analysis, contact **Dory Abelman** or open an issue on this repository.
