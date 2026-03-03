# cfWGS-MM-MRD

Analysis code for **Abelman et al. (2025)**: *"Cell-free DNA Whole Genome Sequencing for Non-Invasive MRD Detection in Multiple Myeloma"*.

The pipeline is a series of numbered R scripts. Each script reads the outputs of the previous step — run them in order.  
**Raw data files are not included in this repository** (see [Data availability](#data-availability) below).

---

## Quick start

```r
# 1. Install dependencies (once)
source("setup_packages.R")

# 2. Set data directory paths for your environment
source("config.R")

# 3. Load shared utility functions
source("helpers.R")

# 4. Run scripts in numbered order (1_0 → 1_1A → … → 4_3)
```

---

## Pipeline scripts

Scripts are numbered to indicate execution order. Run them sequentially from a working directory set to the project root.

### Stage 0 — Cohort selection

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `0_1_select_additional_samples_for_cohort_expansion.R` | Identifies patients with a sample collected within ±30 days of diagnosis; classifies all samples as Baseline / MRD / Treatment. | IMMAGINE patient ID mapping and sample inventory Excel files (`Clinical data/IMMAGINE/`). | Sample eligibility classification table. |

### Stage 1 — Data processing

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `1_0_Process_clinical_metadata.R` | Harmonise clinical metadata across SPORE, M4, IMMAGINE, and MyC cohorts. Annotate timepoints, compute PFS metrics, export master metadata table. | Clinical metadata spreadsheets from `Clinical data/` and `M4_CMRG_Data/`. | `combined_clinical_data_updated_April2025.csv` |
| `1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R` | Parse transplant dates, relapse events, OS follow-up, and FISH flags. | SPORE, M4, and IMMAGINE clinical Excel workbooks with ASCT and progression data. | `tidy_timepoints.csv`, `tidy_treatments.csv`, `tidy_progression.csv`, `ASCT_relapse_summary.csv` |
| `1_1B_Process_clinical_labs.R` | Compile laboratory measurements (CBC, serum protein, staging) for M4 cohort. | M4 clinical lab spreadsheets (`M4_COHORT_*.xlsx`). | Master clinical labs table. |
| `1_2_Process_Mutation_Data.R` | Ingest bone-marrow and cfDNA MAF files, compute VAFs, filter to diagnosis/baseline samples. | WGS MAF files (`*.maf`). | `combined_maf_bm_all_muts.rds`, `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts.rds` |
| `1_2_Part2_Get_Mutation_Overlap.R` | Compute per-patient overlap between bone-marrow and cfDNA mutation calls; generate Venn diagrams and overlap barplots. | `combined_maf_bm_dx.rds`, `combined_maf_blood_all_muts.rds` | VennDiagram PNGs, `percent_overlap_barplot.png` |
| `1_3_Process_Ig_Translocation_Info.R` | Parse Ig-caller structural variant calls, map to cytobands, flag MM-relevant translocations (t(11;14), t(4;14) etc.), and filter using manually IGV-verified calls. | Ig-caller `*_filtered.tsv` files; `Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm` (IGV-verified calls, `Looks_real` confidence column). | `translocation_data_cytoband_updated.rds`, translocation count summary plots. |
| `1_4_Process_CNA_Data.R` | Summarise arm-level copy-number alterations from ichorCNA segments. Produces binary matrix for del1p, amp1q, del13q, del17p, and hyperdiploidy. | ichorCNA `*.seg` files (`Oct 2024 data/Ichor_CNA/`). | `Oct_2024_combined_corrected_calls.rds/csv`, `myeloma_CNA_matrix_with_HRD` |
| `1_4A_Process_sequenza_CNA_Data.R` | Ploidy-aware arm-level CNA calling from Sequenza segments, with realistic hyperdiploidy criteria for myeloma. | Sequenza `*_segments.txt` files. | `cna_data_from_sequenza_400.rds`, per-sample ploidy estimates CSV. |
| `1_5_Integrate_WGS_Feature_Data.R` | Merge CNA, translocation, tumour fraction, and mutation features into a single `All_feature_data` table. | Outputs from steps 1.2–1.4 and ichorCNA tumour fraction. | `All_feature_data_August2025.rds/txt`, `mutation_export.rds/txt` |
| `1_6_Identify_High_Quality_Patient_Pairs.R` | Flag patients with both BM and cfDNA data; produce patient-level sample availability and cohort assignment tables. | `All_feature_data_August2025.rds`, processing log XLSX, clinical metadata. | `patient_cohort_assignment.rds/csv`, `high_quality_patients_list.rds/csv` |
| `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R` | Compute nucleosome-accessibility metrics (coverage, amplitude, z-scores vs healthy controls) from cfWGS data. | Nucleosome-distance `.tsv` files from `Fragmentomics_data/` and `Normals/`. | `griffin_per_site_metrics.tsv`, `griffin_per_site_stats.tsv`, `MM_DARs_chromatin_activation_data.csv` |
| `1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R` | Compute fragment-size (short-fragment) scores and merge with nucleosome accessibility and clinical data. | `insert_size_summary.tsv`, `fragment_scores.tsv`, `MM_DARs_chromatin_activation_data.csv`. | `Key_fragmentomics_data_updated.csv` |
| `1_7C_Process_fragmentomics_data_dilution_series_updated.R` | Repeat fragmentomics processing (steps 1.7A–B) specifically for the experimental dilution-series samples (physically diluted in the lab). | Dilution-series nucleosome-distance files and fragmentomics summaries. | `key_fragmentomics_info_dilution_series.csv/rds` |
| `1_8_Process_Cumulative_VAFs_MRDetect.R` | Process MRDetect cfDNA mutation-detection output; compute z-scores relative to CHARM healthy controls; filter to patient timepoints. | MRDetect CSV output files (`MRDetect_outputs/*.csv`). | `cfWGS_Winter2025All_MRDetect_with_Zscore.rds/txt`, BM and Blood processed CSVs. |
| `1_8A_Process_Cumulative_VAFs_for_dilution_series.R` | Identical processing pipeline to `1_8`, applied to the dilution-series MRDetect outputs. | MRDetect CSV files for dilution-series samples. | `cfWGS_Winter2025Dilution_series_with_zscore.rds` |
| `1_9_Create_dilution_series_eligibility_table.R` | Identify eligible patient pairs for the experimental dilution series: one "tumor-high" (≥0.5% detection rate) and one "tumor-low" (≤0.05% AND z-score < 2 vs healthy controls) timepoint per patient. Computes physical mixing fractions for target tumor fractions 10⁻¹ to 10⁻⁶. | `cfWGS_MRDetect_BM_data_updated_Feb2026.csv` (from step 1.8). | `eligible_dilution_pairs_Feb2026.csv`, `dilution_plan_raw_Feb2026.csv`, `dilution_plan_diff_vs_hc_Feb2026.csv` |

### Stage 2 — Summary statistics and baseline figures

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `2_0_Assemble_Table_With_All_Features.R` | Merge MRD assay results (MFC, clonoSEQ, EasyM), cfWGS metrics, clinical data, and fragmentomics into a single master table. | Outputs from all Stage 1 scripts. | `Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated*.rds` |
| `2_1_Clinical_Demographics_Table.R` | Build Table 1 (patient demographics, disease characteristics by cohort). | Master feature table from `2_0`. | `table1_categorical_updated_final.docx`, `cohort_assignment_table.rds` |
| `2_1_Part2_Cohort_Swim_Plot.R` | Generate treatment-timeline swim plot for the full cohort. | `tidy_treatments.csv`, M4 and IMMAGINE chemotherapy tables. | Swim-plot PDF/PNG. |
| `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Create Figure 1 integrated alteration heatmap (BM overlaid with cfDNA; mutations, CNAs, translocations). | `Final_aggregate_table*.rds`, cohort assignment, CNA/translocation/mutation RDS files. | Figure 1 heatmap PDFs/PNGs. |
| `2_3_Feature_Concordance_And_Mutation_Counts.R` | Compute FISH–WGS concordance; summarise baseline mutation counts by cohort; scatterplots of mutation burden vs. tumour fraction. | Master feature table, cohort assignments, mutation export RDS. | Concordance tables (CSV/XLSX), concordance and mutation-count figures. |
| `2_4_Longitudinal_features_analysis.R` | Summarise how MRD features change over time; generate spaghetti and violin plots of longitudinal cfWGS measurements. | Master feature table, cohort assignments. | Longitudinal comparison CSVs and figures. |

### Stage 3 — Model training, threshold optimisation, and validation

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `3_1_Optimize_cfWGS_thresholds.R` | **Core analysis.** Train elastic-net classifiers using 5×5 nested cross-validation on BM-derived, blood-derived, and fragmentomics features. Evaluate ROC/AUC, calibration, and sensitivity/specificity at fixed thresholds. | Master feature table, clonoSEQ/MFC ground-truth labels. | `selected_combo_models_*.rds`, `selected_combo_thresholds_*.rds`, `all_patients_with_BM_and_blood_calls_updated*.rds` |
| `3_1_A_Process_and_optimize_EasyM.R` | Process EasyM (proteomic MRD) quantitative data; optimise clearance thresholds at each timepoint; generate EasyM-vs-cfWGS comparison tables and Kaplan-Meier plots. | EasyM CSV files from clinical collaborators; cfWGS call table from `3_1`. | EasyM optimised-call CSV, landmark-analysis figures, concordance tables. |
| `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Apply trained models and thresholds to the dilution-series samples to establish the LOD (limit of detection) for each feature and combined model. | Saved models/thresholds from `3_1`; fragmentomics and MRDetect dilution-series outputs. | LOD figures (Extended Data 5D, 7D), source-data CSVs. |
| `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Generate publication figures using the tumour-informed cfWGS model: ROC curves, density plots, waterfall plots, calibration, decision-curve analysis, and contingency tables. | `all_patients_with_BM_and_blood_calls_updated*.rds`, threshold table. | Main manuscript figure panels (PNG, 300 dpi). |
| `3_3_Plot_optimal_cutoff_tumor_naive_calls_and_clinical_concordance.R` | Same figure set as `3_2` but for the **tumour-naïve** (z-score only, no patient-specific mutation sites) blood cfDNA model. | `all_patients_with_BM_and_blood_calls_updated2.rds`, threshold table. | Tumour-naïve figure panels (PNG). |

### Stage 4 — Clinical outcome analyses

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `4_1_Survival_Analysis.R` | Kaplan-Meier PFS curves stratified by MRD status at landmark timepoints (Post-ASCT, 1yr Maintenance). Calculate sensitivity of each assay for detecting future relapse. | `all_patients_with_BM_and_blood_calls_updated*.rds`, EasyM calls, `Censor_dates_per_patient_for_PFS_updated.rds`. | KM-curve PNGs, sensitivity tables (CSV). |
| `4_2_Compare_subclonal_evolution.R` | Identify emergent CNA events between baseline and relapse cfDNA samples (subclonal evolution). | `All_feature_data_August2025.rds`, cohort assignments. | `Subclonal_evolution_plots.pdf`, `Emergent_CNA_events.csv`. |
| `4_3_cfWGS_vs_EasyM_Proteomic_MRD_Comparison.R` | Generate figures directly comparing cfWGS tumour-informed calls vs. EasyM at paired post-ASCT and 1-year maintenance timepoints. | cfWGS call table, EasyM quantitative and binary CSVs, PFS censor dates. | Comparison figures and concordance tables. |

---

## Helper files

| File | Role |
|------|------|
| `setup_packages.R` | Installs/loads all R packages required by the pipeline. Run once before using any numbered script. |
| `config.R` | Complete package list and common directory path constants (`clinical_data_dir`, `wgs_results_dir`, `output_tables_dir`). Update paths here before running. |
| `helpers.R` | Shared utility functions used across scripts (e.g. `clean_sample_id()` for standardising sample identifiers). |

---

## Directory structure

```
Final_Scripts/
├── 0_1_*.R … 4_3_*.R      ← numbered analysis scripts (primary pipeline)
├── config.R                ← package list and path configuration
├── helpers.R               ← shared utility functions
├── setup_packages.R        ← one-shot package loader
└── README.md               ← this file
```

---

## Running the analysis

**Prerequisites**

- R ≥ 4.2
- Run `source("setup_packages.R")` once to install all required packages. Key libraries: **tidyverse**, **readxl**, **data.table**, **lubridate**, **ComplexHeatmap**, **circlize**, **ChromHeatMap**, **maftools**, **GenomicRanges**, **gtsummary**, **officer**, **flextable**, **gt**, **pROC**, **patchwork**, **rmda**, **exact2x2**, **survival**, **survminer**, **tableone**, **ggbreak**, **ggridges**, **ggpubr**, **GGally**, **pbapply**, **PRROC**, **VennDiagram**, **caret**, **glmnet**, **Matrix**, **DescTools**, **GeneCycle**, **RColorBrewer**, **fuzzyjoin**, **rstatix**, **openxlsx**, **scales**, **viridis**, **writexl**, and others. See `config.R` for the full list.

**Steps**

1. Obtain raw data files (see [Data availability](#data-availability)) and place them in the directory structure expected by each script (paths are documented in each file's header).
2. Update `config.R` with your local path for `clinical_data_dir`, `wgs_results_dir`, and `output_tables_dir`.
3. Source the helpers at the start of your R session:
   ```r
   source("setup_packages.R")
   source("config.R")
   source("helpers.R")
   ```
4. Execute scripts in numbered order (`0_1` → `1_0` → `1_1A` → … → `4_3`). Each script writes its outputs to `Output_tables_2025/` or `Final Tables and Figures/` (created automatically if absent).

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

- **CSV / XLSX tables** — cleaned metadata, mutation counts, MRD feature matrices, concordance statistics, sensitivity/specificity tables (written to `Output_tables_2025/`)
- **PDF / PNG figures** — all manuscript and supplementary figures (written to `Final Tables and Figures/`)
- **RDS objects** — intermediate R data objects passed between scripts

---

## Contact

For questions about the analysis, contact **Dory Abelman** or open an issue on this repository.
