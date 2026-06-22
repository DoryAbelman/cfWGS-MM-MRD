# cfWGS-MM-MRD

Analysis code for **Abelman et al. (2025)**: *"Cell-free DNA Whole Genome Sequencing for Non-Invasive MRD Detection in Multiple Myeloma"*.

The pipeline is a series of numbered R scripts. Each script reads the outputs of the previous step and can be run from the command line through the included runner.

This directory is the **primary GitHub/Code Ocean analysis pipeline**. The numbered scripts remain organized by analysis stage and scientific question, matching the original workflow: raw/clinical processing, WGS feature processing, dilution-series analyses, model training/application, concordance analyses, and survival analyses. The scripts directly create both the analysis outputs and organized manuscript-facing outputs.

The normal user-facing workflow is `run_pipeline.R` or `run_manuscript_workflow.R`. Supporting provenance and validation files are included to make the script-to-manuscript mapping explicit, but the numbered scripts are the scientific source of truth.

The workflow is intended to be end-to-end from staged raw/protected inputs through final tables, source-data files, and manuscript-facing figure components. The deliberate exception is routine test-cohort expansion: final training-derived model objects, nested-CV fold objects, metric summaries, and related cached intermediates are treated as preserved reproducibility artifacts unless `--include-cache-sensitive` is supplied. New test-cohort samples should be processed through the upstream feature scripts and then scored with preserved training-derived artifacts, so test-cohort metrics can update without unintentionally recomputing the training-stage model outputs.

**Raw data files are not included in this repository** (see [Data availability](#data-availability) below).

---

## Quick start

Check that required local inputs are staged inside the project before running
the analysis:

```sh
Rscript Scripts_2025/Final_Scripts/prepare_local_inputs.R --check
```

If the required files are available in another local mirror, copy missing inputs
into the project while preserving the relative paths expected by the scripts:

```sh
Rscript Scripts_2025/Final_Scripts/prepare_local_inputs.R --source-dir /path/to/staged_inputs --copy-missing
```

Dry-run the numbered source pipeline from the project root:

```sh
Rscript Scripts_2025/Final_Scripts/run_pipeline.R
```

Execute the pipeline from the project root:

```sh
Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute
```

Numbered scripts also write manuscript-labeled copies of their final
figure/table components to:

```text
Scripts_2025/Final_Scripts/final_manuscript_objects/
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

Validate the organized manuscript-output tree:

```sh
Rscript Scripts_2025/Final_Scripts/validate_manuscript_outputs.R
```

This validator checks the audited source map, direct-output manifest,
human-readable output indexes, final easy-use folders, generated component
mirrors, per-script manuscript helper calls, manuscript-output headers in the
numbered scripts, and README coverage for every script that owns a mapped final
artifact. It writes:

```text
Scripts_2025/Final_Scripts/final_manuscript_objects/manuscript_output_validation_report.tsv
```

The most useful files for seeing what each numbered script produces are:

```text
Scripts_2025/Final_Scripts/final_manuscript_objects/manuscript_output_index.tsv
Scripts_2025/Final_Scripts/final_manuscript_objects/script_output_index.tsv
```

For GitHub or Code Ocean packaging, use
[`CODE_OCEAN_UPLOAD_CHECKLIST.md`](CODE_OCEAN_UPLOAD_CHECKLIST.md).

## Manuscript figure/table ownership

The numbered scripts are still organized by analysis question, not by final
figure number. Each script header lists the same outputs shown below. When a
script runs, the final manuscript-facing files are copied or saved into
`final_manuscript_objects/` under the final figure/table label.

| Script | Final manuscript outputs |
| --- | --- |
| `1_2_Part2_Get_Mutation_Overlap.R` | Extended Data Figure 2G |
| `1_6_Identify_High_Quality_Patient_Pairs.R` | Figure 1B source table |
| `2_1_Clinical_Demographics_Table.R` | Table 1 |
| `2_1_Part2_Cohort_Swim_Plot.R` | Figure 1A; Supplementary Table 1 |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Extended Data Figure 1; Supplementary Table 1A |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Extended Data Figure 2A-F; Supplementary Tables 2-3 |
| `2_4_Longitudinal_features_analysis.R` | Figure 2A-D; Extended Data Figure 3A-C; Extended Data Figure 4 |
| `3_1_Optimize_cfWGS_thresholds.R` | Figure 3A-B; Figure 4A-B; Extended Data Figure 5A-C; Extended Data Figure 7A-C and 7E; Extended Data Figure 9A-F; Supplementary Tables 4-6 |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Figure 3C; Extended Data Figure 5D; Extended Data Figure 7D; Supplementary Table 7 |
| `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Figure 3D-E; Figure 4C-D; Extended Data Figure 5E-G; Extended Data Figure 7F-H; Supplementary Tables 8 and 10 |
| `4_1_Survival_Analysis.R` | Figure 3F; Figure 4E; Extended Data Figure 6A-K; Extended Data Figure 8A-D plus bottom panels; Supplementary Table 9. The submitted time-window outputs for Extended Data Figure 6I, Extended Data Figure 8D, and Supplementary Table 9 are reproduced from preserved `detection_progression_updated5` source CSVs; a separate prospective robust time-window QC analysis is written for future test-cohort updates. |
| `4_2_Compare_subclonal_evolution.R` | Extended Data Figure 10A-B |

The source of truth for exact source filenames, companion CSVs, and historical
name corrections is `docs/manuscript_artifact_source_map.tsv`. The validator
uses that file plus the direct-output manifest to confirm that every mapped
artifact lands in the organized manuscript-output tree.

Run the full manuscript command sequence through one entry point:

```sh
Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R
Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R --execute
Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R --execute --skip-source
```

`run_manuscript_workflow.R` dry-runs by default. With `--execute`, it runs the
numbered `Final_Scripts` pipeline through `run_pipeline.R`, refreshes the
stage-ordered script/artifact map, and validates the direct manuscript outputs
under `Scripts_2025/Final_Scripts/final_manuscript_objects/`.

The source pipeline skips cache-sensitive nested-CV/model training unless
`--include-cache-sensitive` is explicitly supplied.

Run a subset of the numbered source pipeline:

```sh
Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute --from 2_0 --to 2_4
Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute --only 2_1_part2
```

Check packages without running analysis:

```sh
Rscript Scripts_2025/Final_Scripts/prepare_local_inputs.R --check
Rscript Scripts_2025/Final_Scripts/run_pipeline.R --check-packages
```

The runner writes timestamped logs and `run_manifest.tsv` under `Scripts_2025/Final_Scripts/pipeline_logs/`. It skips `3_1_Optimize_cfWGS_thresholds.R` by default because that script contains cache-sensitive nested-CV model training and manuscript values should remain stable during routine test-cohort expansion. Add `--include-cache-sensitive` only when intentional model recomputation is expected; regenerated artifacts may differ unless the exact same inputs, seeds, package versions, and execution environment are preserved.

Each numbered script header lists the manuscript figure/table outputs that the
script creates or updates. Where a final manuscript object is generated, the
script calls `manuscript_output_helpers.R` to place a labeled copy into
`final_manuscript_objects/`.

The runner also refreshes `script_index.tsv`, a table of every numbered script, its stage, run policy, and mapped manuscript figures/tables. See `MANUSCRIPT_OUTPUT_MAP.md` for the concise script-to-figure/table map.

It also refreshes the stage-ordered crosswalk:

```text
../../docs/stage_ordered_script_artifact_map.tsv
../../docs/stage_ordered_script_artifact_map.md
```

Those files keep scripts grouped by analysis stage/question and show the
manuscript figure/table, generated output path, regeneration status, and
validation status for each mapped artifact.

---

## Pipeline scripts

Scripts are numbered to indicate execution order. Run them sequentially from a working directory set to the project root.

### Stage 0 - Cohort selection

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `0_1_select_additional_samples_for_cohort_expansion.R` | Identify candidate test-cohort expansion samples by linking IMMAGINE inventory records to patient IDs, classifying samples as Baseline / MRD / Treatment, and reporting blood-MRD expansion opportunities. These are review/operations outputs, not final manuscript tables. | IMMAGINE patient ID mapping, sample inventory, MRD date, diagnosis date, and treatment-line Excel files (`Clinical data/IMMAGINE/`). | Cohort-expansion review tables under `Output_tables_2025/cohort_expansion/`; matching/duplicate/blood-MRD diagnostics under `Output_tables_2025/cohort_expansion/support_qc/`. |

### Stage 1 - Data processing

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `1_0_Process_clinical_metadata.R` | Harmonise clinical metadata across SPORE, M4, IMMAGINE, and MyC cohorts; annotate timepoints; compute relapse/PFS date helpers; and export the canonical clinical metadata table consumed by downstream mutation, CNA, fragmentomics, swim-plot, model, and survival scripts. Support-only count summaries, BAM review lists, and QC plots are separated from manuscript outputs. | Clinical metadata spreadsheets from `Clinical data/`, `M4_CMRG_Data/`, and reviewed `cohort_assignment_table_updated.rds`. | Canonical downstream inputs retained at historical paths, especially `combined_clinical_data_updated_April2025.csv`, `M4_labs_cleaned.csv`, `Relapse_dates_M4_clean.csv`, `Relapse dates cfWGS updated.csv`, and PFS/censor-date files under `Exported_data_tables_clinical/`; support-only summaries under `Output_tables_2025/clinical_metadata_support/`. |
| `1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R` | Parse transplant dates, relapse/progression events, OS follow-up, and FISH flags for SPORE, M4, and IMMAGINE. This is an upstream clinical-helper step, not a direct manuscript figure/table owner. | SPORE, M4, and IMMAGINE clinical Excel/CSV workbooks with ASCT, progression, FISH, and OS follow-up data; canonical clinical and relapse-date outputs from `1_0`. | Active downstream clinical inputs in `Clinical data/SPORE/tidy_*.csv` and `Clinical data/Exported clinical data April 2025/*.csv`; support-only relapse/OS helper tables and QC plots under `Output_tables_2025/clinical_followup_support/`. |
| `1_1B_Process_clinical_labs.R` | Compile cross-cohort clinical laboratory, staging, demographic, and FISH/cytogenetic helper fields for M4, SPORE, and IMMAGINE. This is an upstream clinical-helper step, not a direct manuscript figure/table owner. | M4 clinical workbooks (`M4_COHORT_*.xlsx`), curated SPORE/IMMAGINE clinical sheets, and FISH helper files generated by `1_1A`. | Active downstream input `Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv`, consumed by `2_0_Assemble_Table_With_All_Features.R`. |
| `1_2_Process_Mutation_Data.R` | Ingest bone-marrow and cfDNA MAF files, compute VAFs, annotate with clinical metadata, and create the mutation-call intermediates used by downstream mutation-overlap, heatmap, feature-integration, and MRDetect analyses. VAF/read-depth plots are support-only QC, not manuscript panels. | BM and cfDNA WGS MAF files (`*.maf`), `combined_clinical_data_updated_April2025.csv`. | Core intermediates `combined_maf_bm_all_muts.rds`, `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts_updated.rds`, and temporary MAF exports used by downstream scripts; support-only QC plots under `Final Tables and Figures/mutation_processing_support/`. |
| `1_2_Part2_Get_Mutation_Overlap.R` | Compute per-patient overlap between baseline bone-marrow and cfDNA mutation calls; generate the protected-ID lollipop panel used for Extended Data Figure 2G. Support-only Venn diagrams and QA barplots are written under `Final Tables and Figures/mutation_overlap_support/`. | `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts_updated.rds`, `combined_clinical_data_updated_April2025.csv`, `cohort_assignment_table_updated.rds`, `id_map.rds` | Extended Data Figure 2G panel PNG plus source-data and summary CSVs; support-only Venn/barplot outputs are not copied to `final_manuscript_objects/`. |
| `1_3_Process_Ig_Translocation_Info.R` | Parse Ig-caller structural variant calls, map breakpoints to cytobands, flag MM-relevant translocations, generate IGV review helpers, and export the manually IGV-confirmed translocation feature matrix. This is an upstream feature-generation step, not a direct manuscript figure/table owner. | Ig-caller `*_filtered.tsv` files when available, otherwise the preserved parsed Ig-caller cache; cytoband and ENCODE blacklist references; `Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm` with `Looks_real` review calls. | Active downstream inputs `Jan2025_exported_data/translocation_data_cytoband_updated.rds/txt`; support-only parsing caches, pre-IGV matrices, and QC plots under `Output_tables_2025/translocation_processing_support/` and `Output_figures_2025/translocation_processing_support/`. |
| `1_4_Process_CNA_Data.R` | Active upstream ichorCNA processing step. Summarises raw ichorCNA segment calls into arm-level del1p, amp1q, del13q, del17p, and hyperdiploidy features; also creates ichorCNA calls at clinical FISH probe loci. It does not directly own a final figure/table, but its outputs feed WGS feature integration, baseline feature heatmaps, FISH/WGS concordance, and Supplementary Table 2 support. | Raw ichorCNA `*.seg` files in `Oct 2024 data/Ichor_CNA/`; if those are unavailable in a lightweight review bundle, the preserved combined segment cache under `Output_tables_2025/ichor_cna_processing_support/` is used. | Active downstream inputs `Jan2025_exported_data/cna_data_ichorCNA.rds/txt` and `Jan2025_exported_data/FISH_probe_calls_bin_cytoband_ichorCNA.rds/txt`; support/QC segment cache and hyperdiploidy audit tables under `Output_tables_2025/ichor_cna_processing_support/`. |
| `1_4A_Process_sequenza_CNA_Data.R` | Active upstream Sequenza CNA processing step. Converts raw Sequenza segments and confints files into ploidy-aware arm-level CNA calls, direct-coordinate FISH-probe helper calls, and sample-level purity/ploidy estimates. It does not directly own a final figure/table, but its outputs feed WGS feature integration, FISH/WGS concordance, ploidy summaries, and Supplementary Table 2 support. | Raw Sequenza `*_segments.txt(.gz)` files in `Oct 2024 data/Sequenza/All_Segments_400/`, matching `*_confints_CP.txt` files in `Oct 2024 data/Sequenza/All_confints_400/`, clinical metadata, cytobands, and FISH probe locations. | Active downstream inputs `Jan2025_exported_data/cna_data_from_sequenza_400_updated.rds/txt`, `Jan2025_exported_data/FISH_data_from_sequenza_400_updated.rds/txt`, and `Jan2025_exported_data/Sample_ploidy_from_sequenza_400.rds`; support/QC segment, probe, cytoband-expanded FISH, ploidy, and summary caches under `Output_tables_2025/sequenza_cna_processing_support/`. |
| `1_5_Integrate_WGS_Feature_Data.R` | Active upstream WGS feature-integration step. Merges mutation, ichorCNA/Sequenza CNA, IG translocation, FISH-probe CNA, tumor-fraction, clinical sample metadata, and reviewed IGV translocation overrides into the final integrated WGS feature tables. It does not directly own a final figure/table, but it is the main bridge into aggregate clinical table assembly, baseline heatmaps, concordance analyses, model training/application, and longitudinal analyses. | Current outputs from `1_2`, `1_3`, `1_4`, and `1_4A`; `Oct 2024 data/tumor_fraction_cfWGS.txt`; reviewed IGV translocation spreadsheet; clinical metadata. | Active downstream inputs `Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds/txt`, `Jan2025_exported_data/CNA_translocation_Sep2025_updated2.rds/txt`, `Jan2025_exported_data/CNA_at_FISH_sites_combined.rds/txt`, `Jan2025_exported_data/mutation_export_updated2.rds/txt`, and `Jan2025_exported_data/mutation_export_updated_more_info2.rds/txt`; support/QC checks under `Output_tables_2025/feature_integration_support/`. |
| `1_6_Identify_High_Quality_Patient_Pairs.R` | Define high-quality baseline sample eligibility, cohort assignment, and the Figure 1B sample-flow source table. The final Figure 1B visual panel remains manually assembled from this exported table, but the underlying counts/QC annotations are regenerated from the command line. | `All_feature_data_August2025.rds`, processing log XLSX, clinical metadata, sample-flow workbook, reviewed `cohort_assignment_table_updated.rds`. | Figure 1B source table in `Final Tables and Figures/` and `final_manuscript_objects/`; downstream intermediates `Output_tables_2025/high_quality_patients_list_for_baseline_mut_calling2.*` and `patient_cohort_assignment.*`; support-only processing log under `Output_tables_2025/sample_qc_support/`. |
| `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R` | Compute nucleosome-accessibility metrics (coverage, amplitude, z-scores vs healthy controls) from cfWGS data. | Nucleosome-distance `.tsv` files from `Fragmentomics_data/` and `Normals/`. | `griffin_per_site_metrics.tsv`, `griffin_per_site_stats.tsv`, `MM_DARs_chromatin_activation_data.csv` |
| `1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R` | Compute fragment-size (short-fragment) scores and merge with nucleosome accessibility and clinical data. | `insert_size_summary.tsv`, `fragment_scores.tsv`, `MM_DARs_chromatin_activation_data.csv`. | `Key_fragmentomics_data_updated.csv` |
| `1_7C_Process_fragmentomics_data_dilution_series_updated.R` | Repeat fragmentomics processing (steps 1.7A–B) specifically for the experimental dilution-series samples (physically diluted in the lab). | Dilution-series nucleosome-distance files and fragmentomics summaries. | `key_fragmentomics_info_dilution_series.csv/rds` |
| `1_8_Process_Cumulative_VAFs_MRDetect.R` | Process MRDetect cfDNA mutation-detection output; compute z-scores relative to CHARM healthy controls; filter to patient timepoints. | MRDetect CSV output files (`MRDetect_outputs/*.csv`). | `cfWGS_Winter2025All_MRDetect_with_Zscore.rds/txt`, BM and Blood processed CSVs. |
| `1_8A_Process_Cumulative_VAFs_for_dilution_series.R` | Identical processing pipeline to `1_8`, applied to the dilution-series MRDetect outputs. | MRDetect CSV files for dilution-series samples. | `cfWGS_Winter2025Dilution_series_with_zscore.rds` |
| `1_9_Create_dilution_series_eligibility_table.R` | Identify eligible patient pairs for the experimental dilution series: one "tumor-high" (≥0.5% detection rate) and one "tumor-low" (≤0.05% AND z-score < 2 vs healthy controls) timepoint per patient. Computes physical mixing fractions for target tumor fractions 10⁻¹ to 10⁻⁶. | `cfWGS_MRDetect_BM_data_updated_Feb2026.csv` (from step 1.8). | `eligible_dilution_pairs_Feb2026.csv`, `dilution_plan_raw_Feb2026.csv`, `dilution_plan_diff_vs_hc_Feb2026.csv` |

### Stage 2 - Summary statistics and baseline figures

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `2_0_Assemble_Table_With_All_Features.R` | Merge MRD assay results (MFC, clonoSEQ, EasyM), cfWGS metrics, clinical data, and fragmentomics into a single master table. | Outputs from all Stage 1 scripts. | `Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated*.rds` |
| `2_1_Clinical_Demographics_Table.R` | Build Table 1 (patient demographics, disease characteristics by cohort). | Master feature table from `2_0`. | `table1_categorical_updated_final.docx`, `cohort_assignment_table.rds` |
| `2_1_Part2_Cohort_Swim_Plot.R` | Generate the treatment-timeline swim plot and privacy-protected event table. Support-only cohort/legend QA exports are written under `Final Tables and Figures/swim_plot_support/`. | `tidy_treatments.csv`, M4 and IMMAGINE chemotherapy tables. | Figure 1A component; Supplementary Table 1 event table. Support QA exports are not copied to `final_manuscript_objects/`. |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Create the baseline integrated alteration heatmap (BM overlaid with cfDNA; mutations, CNAs, translocations) and disease-feature catalog. | `Final_aggregate_table*.rds`, cohort assignment, CNA/translocation/mutation RDS files. | Extended Data Figure 1 component; Supplementary Table 1A. |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Compute FISH-WGS concordance; summarise baseline mutation counts by cohort; build concordance and feature-correlation exports. BAM archive/unarchive helper tables are support-only and run only when `CFWGS_RUN_BAM_ARCHIVE_DIAGNOSTICS=true`. | Master feature table, cohort assignments, mutation export RDS. Optional BAM diagnostic mode also reads `All_bam_storage_locations.xlsx`. | Extended Data Figure 2 components; Supplementary Tables 2 and 3. Optional BAM diagnostic files are written under `Output_tables_2025_updated/support_only_bam_archive_diagnostics/` and are not copied to `final_manuscript_objects/`. |
| `2_4_Longitudinal_features_analysis.R` | Summarise how cfWGS features change over time; generate longitudinal mutation, CNA, and fragmentomics panels. | Master feature table, cohort assignments. | Figure 2 components; Extended Data Figures 3 and 4 source panels. |

### Stage 3 - Model training, model application, and validation

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `3_1_Optimize_cfWGS_thresholds.R` | **Core model analysis.** Train or load elastic-net classifiers using 5x5 nested cross-validation on BM-derived, blood-derived, and fragmentomics features. Evaluate ROC/AUC, calibration, and sensitivity/specificity at fixed operating points. Skipped by default by the command-line runners because these artifacts are cache-sensitive preserved reproducibility outputs during routine test-cohort expansion. | Master feature table, clonoSEQ/MFC ground-truth labels. | Figure 3A-B and Figure 4A-B components; Extended Data Figures 5A-C, 7A-C/E, and 9A-F; Supplementary Tables 4-6; preserved model and metric RDS files. |
| `3_1_A_Process_and_optimize_EasyM.R` | Process EasyM (proteomic MRD) quantitative data and export the optimized EasyM call/threshold tables required by `3_2` and `4_1`. | EasyM CSV files from clinical collaborators; cfWGS call table from `3_1`. | Required intermediate CSVs: `EasyM_all_samples_with_optimized_calls.csv` and `EasyM_threshold_values_by_timepoint.csv`; support-only landmark-analysis tables/figures. |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Apply trained models and thresholds to the dilution-series samples to establish the LOD (limit of detection) for each feature and combined model. | Saved models/thresholds from `3_1`; fragmentomics and MRDetect dilution-series outputs. | Figure 3C component; Extended Data Figures 5D and 7D; Supplementary Table 7. |
| `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Generate tumour-informed cfWGS clinical-concordance figures: assay positivity, model-vs-clinical assay comparisons, calibration/decision-curve support, and contingency tables. | `all_patients_with_BM_and_blood_calls_updated*.rds`, threshold table. | Figure 3D-E and Figure 4C-D components; Extended Data Figures 5E-G and 7F-H; Supplementary Tables 8 and 10. |
| `archive/support_analysis/3_3_Plot_optimal_cutoff_tumor_naive_calls_and_clinical_concordance.R` | Archived tumour-naive blood cfDNA support/sensitivity analysis. This is not required for routine manuscript regeneration because the final clinical-concordance figures/tables are produced by `3_2`, and no current final manuscript artifact is mapped to `3_3`. | `all_patients_with_BM_and_blood_calls_updated*.rds`, threshold table. | Support-only tumour-naive review outputs; not copied to `final_manuscript_objects/`. |

### Stage 4 - Clinical outcome analyses

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `4_1_Survival_Analysis.R` | Kaplan-Meier PFS curves stratified by MRD status at landmark timepoints (Post-ASCT, 1yr Maintenance). Calculate sensitivity of each assay for detecting future relapse. For exact manuscript reproduction, the submitted time-window outputs use preserved `detection_progression_updated5` source CSVs. For future test-cohort additions, the script also exports a prospective robust QC analysis that requires either a future progression event or adequate follow-up through each prediction window. PFS event/censor dates remain separate from patient-level last-follow-up dates; the prospective QC uses `patient_followup_dates_updated.*` generated by `1_0_Process_clinical_metadata.R` when available. | `all_patients_with_BM_and_blood_calls_updated*.rds`, EasyM calls, `Censor_dates_per_patient_for_PFS_updated.rds`, `patient_followup_dates_updated.rds`, curated progression-date CSVs. | Figure 3F and Figure 4E components; Extended Data Figures 6 and 8; Supplementary Table 9; prospective QC outputs under `Output_tables_2025/detection_progression_updated6/prospective_timewindow_qc/`. |
| `4_2_Compare_subclonal_evolution.R` | Identify emergent CNA events between baseline and relapse cfDNA samples. The final Extended Data Figure 10 genome-wide CNA tracks come from an external ichorCNA plotting workflow; this script provides repo-side supporting event outputs. | `All_feature_data_August2025.rds`, cohort assignments. | `Emergent_CNA_events.csv` and subclonal-evolution support outputs for Extended Data Figure 10. |
| `archive/support_analysis/4_3_cfWGS_vs_EasyM_Proteomic_MRD_Comparison.R` | Archived cfWGS-vs-EasyM support/sensitivity analysis. This is not required for routine manuscript regeneration because final EasyM call generation is handled by `3_1_A`, and final survival/relapse EasyM panels are handled by `4_1`. | cfWGS call table, EasyM quantitative and binary CSVs, PFS censor dates. | Support-only EasyM comparison figures/tables; not copied to `final_manuscript_objects/`. |

---

## Helper files

| File | Role |
|------|------|
| `setup_packages.R` | Checks and loads R packages required by the pipeline. It stops with a package-specific error if anything is missing. |
| `config.R` | Complete package list and common directory path constants (`clinical_data_dir`, `wgs_results_dir`, `output_tables_dir`). Update paths here before running. |
| `helpers.R` | Shared utility functions used across scripts (e.g. `clean_sample_id()` for standardising sample identifiers). |
| `prepare_local_inputs.R` | Local/Code Ocean preflight utility. It checks that required protected/de-identified inputs are staged under project-relative paths and can copy missing files/directories from a local source mirror with `--source-dir ... --copy-missing`. |
| `manuscript_output_helpers.R` | Shared utilities used inside numbered scripts to copy/save final manuscript figure, table, and source-data components into `final_manuscript_objects/` with audited manuscript labels. |
| `validate_manuscript_outputs.R` | Command-line validator for the direct manuscript-output tree. It checks coverage against `docs/manuscript_artifact_source_map.tsv`, verifies organized paths, and writes `manuscript_output_validation_report.tsv`. |
| `pipeline_metadata.R` | Reads `docs/manuscript_artifact_source_map.tsv` and builds script-to-artifact crosswalks used by the command-line runner and documentation. |

---

## Directory structure

```
Final_Scripts/
├── 0_1_*.R … 4_2_*.R      ← active numbered analysis scripts (primary manuscript pipeline)
├── archive/support_analysis/ ← non-required support/sensitivity scripts retained for audit
├── config.R                ← package list and path configuration
├── helpers.R               ← shared utility functions
├── manuscript_output_helpers.R ← direct manuscript-output helper for numbered scripts
├── validate_manuscript_outputs.R ← direct manuscript-output validator
├── pipeline_metadata.R      ← shared script/output metadata helpers
├── final_manuscript_objects/ ← manuscript-labeled outputs created at run time
├── run_manuscript_workflow.R ← one-command manuscript workflow orchestrator
├── setup_packages.R        ← one-shot package checker/loader
├── prepare_local_inputs.R  ← local input staging/check utility
└── README.md               ← this file
```

---

## Running the analysis

**Prerequisites**

- R ≥ 4.2
- Run `Rscript Scripts_2025/Final_Scripts/run_pipeline.R --check-packages` to check all required packages. Key libraries: **tidyverse**, **readxl**, **data.table**, **lubridate**, **ComplexHeatmap**, **circlize**, **ChromHeatMap**, **maftools**, **GenomicRanges**, **gtsummary**, **officer**, **flextable**, **gt**, **pROC**, **patchwork**, **rmda**, **exact2x2**, **survival**, **survminer**, **tableone**, **ggbreak**, **ggridges**, **ggpubr**, **GGally**, **pbapply**, **PRROC**, **VennDiagram**, **caret**, **glmnet**, **Matrix**, **DescTools**, **GeneCycle**, **RColorBrewer**, **fuzzyjoin**, **rstatix**, **openxlsx**, **scales**, **viridis**, **writexl**, and others. See `config.R` for the full list.

**Steps**

1. Obtain the staged raw/protected input files (see [Data availability](#data-availability)) and place them in the directory structure expected by each script. Paths and key inputs are documented in the relevant script headers and in the tables above. Use `prepare_local_inputs.R --check` to verify the local staging before running the analysis, or `prepare_local_inputs.R --source-dir /path/to/staged_inputs --copy-missing` to copy missing inputs from a local mirror.
2. Update `config.R` with your local path for `clinical_data_dir`, `wgs_results_dir`, and `output_tables_dir`.
3. Dry-run the command-line sequence:
   ```sh
   Rscript Scripts_2025/Final_Scripts/run_pipeline.R
   ```
4. Execute scripts in numbered order through the runner:
   ```sh
   Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute
   ```
Each script runs in a fresh R process from the project root and writes its outputs to the historical analysis locations such as `Output_tables_2025/` and `Final Tables and Figures/`. Scripts under `Scripts_2025/Final_Scripts/archive/support_analysis/` are retained for audit/support review and are not selected by `run_pipeline.R`.
5. Inspect direct manuscript-facing exports:
   ```sh
   find Scripts_2025/Final_Scripts/final_manuscript_objects -maxdepth 4 -type f
   ```
   This folder contains manuscript-labeled figures, tables, source-data files,
   manifests, and validation reports produced by the numbered scripts.
6. Validate that the manuscript-facing output tree is complete:
   ```sh
   Rscript Scripts_2025/Final_Scripts/validate_manuscript_outputs.R
   ```

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
- **PDF / PNG figures** - all manuscript and supplementary figures (written to `Final Tables and Figures/`)
- **Support-only QA outputs** - non-manuscript review exports are kept in explicitly named support folders such as `Final Tables and Figures/swim_plot_support/` and `Output_tables_2025_updated/support_only_bam_archive_diagnostics/`
- **RDS objects** - intermediate R data objects passed between scripts

---

## Contact

For questions about the analysis, contact **Dory Abelman** or open an issue on this repository.
