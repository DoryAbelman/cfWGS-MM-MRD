# cfWGS-MM-MRD

This repository contains the analysis scripts used in **Abelman et al. (2025)** for the study *"Cell-free DNA Whole Genome Sequencing for Non‑Invasive MRD Detection in Multiple Myeloma"*.

The code is organised as a series of numbered R scripts that process raw clinical data, integrate whole genome sequencing features, and generate the figures found in the manuscript. The raw data files referenced by these scripts are not included in this repository.

## Repository layout

| Step | Script | Purpose | Expected input |
| --- | --- | --- | --- |
| 1.0 | `1_0_Process_clinical_metadata.R` | Harmonise clinical metadata from multiple cohorts and export cleaned tables. | Clinical metadata spreadsheets (`*.xlsx`, `*.csv`). |
| 1.1A | `1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R` | Parse transplant dates, relapse events, and survival outcomes. | SPORE, M4 and IMMAGINE clinical Excel workbooks with ASCT and progression data. |
| 1.1B | `1_1B_Process_clinical_labs.R` | Compile laboratory measurements for the cohort. | M4 clinical lab spreadsheets (`M4_COHORT_*.xlsx`). |
| 1.2 | `1_2_Process_Mutation_Data.R` | Prepare mutation calls for baseline and follow‑up samples. | WGS mutation files (MAF). |
| 1.2B | `1_2_Part2_Get_Mutation_Overlap.R` | Compare bone‑marrow and cfDNA mutations. | Pair of bone‑marrow and cfDNA MAF files per patient. |
| 1.3 | `1_3_Process_Ig_Translocation_Info.R` | Extract immunoglobulin translocation calls from WGS data. | Structural variant calls from WGS. |
| 1.4 | `1_4_Process_CNA_Data.R` | Summarise copy number features from WGS results. | Copy number segmentation files. |
| 1.5 | `1_5_Integrate_WGS_Feature_Data.R` | Combine mutation, CNA, and translocation features. | Outputs from steps 1.2–1.4. |
| 1.6 | `1_6_Identify_High_Quality_Patient_Pairs.R` | Annotate samples that have matched bone‑marrow and cfDNA data. | Integrated feature tables plus clinical metadata. |
| 1.7A | `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R` | Derive nucleosome accessibility metrics. | Nucleosome‑distance tables for samples and healthy controls. |
| 1.7B | `1_7B_Process_fragment_score_and_integrate_fragmentomics_data.R` | Compute fragment‑size scores and merge fragmentomics metrics. | `insert_size_summary.tsv`, `fragment_scores.tsv`, and MM‑DARs activation data. |
| 1.7C | `1_7C_Process_fragmentomics_data_dilution_series_updated.R` | Process fragmentomics for dilution‑series samples. | Dilution‑series nucleosome‑distance files and fragmentomics summaries. |
| 1.8 | `1_8_Process_Cumulative_VAFs_MRDetect.R` | Process MRDetect output and calculate z‑scores. | MRDetect CSV output for patient samples. |
| 1.8A | `1_8A_Process_Cumulative_VAFs_for_dilution_series.R` | Process MRDetect output for dilution series. | MRDetect CSV output for dilution‑series samples. |
| 2.0 | `2_0_Assemble_Table_With_All_Features.R` | Merge all processed features into a single table for downstream analyses. | Outputs from steps 1.*. |
| 2.1 | `2_1_Clinical_Demographics_Table.R` | Produce patient demographic tables used in the manuscript. | Cleaned clinical tables. |
| 2.1B | `2_1_Part2_Cohort_Swim_Plot.R` | Generate cohort treatment timeline swim plot. | `tidy_treatments.csv` and M4 chemotherapy Excel data. |
| 2.2 | `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Create integrated alteration heatmaps. | Feature table from step 2.0. |
| 2.3 | `2_3_Feature_Concordance_And_Mutation_Counts.R` | Calculate concordance between bone‑marrow and cfDNA calls. | Aggregate feature tables and cohort assignments. |
| 2.4 | `2_4_Longitudinal_features_analysis.R` | Summarise longitudinal MRD features and paired changes. | Final aggregated table and cohort assignments. |
| 3.1 | `3_1_Optimize_cfWGS_thresholds.R` | Determine optimal cfWGS tumour‑fraction thresholds. | Feature matrix with clonoSEQ/MFC labels. |
| 3.1B | `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` | Score dilution‑series samples with trained models. | Saved models, thresholds and dilution‑series fragmentomics. |
| 3.2 | `3_2_Plot_optimal_cutoff_and_clinical_concordance.R` | Plot thresholds and clinical concordance panels. | `all_patients_with_BM_and_blood_calls.rds` from step 3.1. |
| 3.3 | `3_3_Plot_optimal_cutoff_tumor_naive_calls_and_clinical_concordance.R` | Tumour‑naïve ROC and concordance figures. | Tumour‑naïve call dataset with clinical labels. |
| 4.1 | `4_1_Survival_Analysis.R` | Perform survival analyses. | `Censor_dates_per_patient_for_PFS_updated.rds` and call table. |
| 4.2 | `4_2_Compare_subclonal_evolution.R` | Assess emergent subclones in longitudinal samples. | `All_feature_data_*.rds` and cohort assignments. |

Each script is intended to be run individually after placing the appropriate input files in the paths referenced within the code. Figures and intermediate tables are written to the working directory.

A small set of helper scripts is included:
`setup_packages.R` loads all required libraries, `config.R` defines common directory locations, and `helpers.R` provides utility functions such as `clean_sample_id()`.

## Running the analysis

1. Ensure you have an R installation (\>=4.2) with all packages used throughout the scripts. Key libraries include **tidyverse** (dplyr, ggplot2, readr, tidyr, tibble, purrr, stringr, forcats), **readxl**, **data.table**, **lubridate**, **ComplexHeatmap**, **circlize**, **ChromHeatMap**, **maftools**, **GenomicRanges**, **gtsummary**, **officer**, **flextable**, **gt**, **pROC**, **patchwork**, **rmda**, **exact2x2**, **survival**, **survminer**, **tableone**, **ggbreak**, **ggridges**, **ggpubr**, **GGally**, **pbapply**, **PRROC**, **VennDiagram**, **caret**, **glmnet**, **Matrix**, **DescTools**, **GeneCycle**, **RColorBrewer**, **fuzzyjoin**, **rstatix**, **openxlsx**, **scales**, **viridis**, **viridisLite**, **writexl**, and others. See `config.R` for the complete package list.
2. Place the raw data files in the directory structure expected by the scripts (see comments at the top of each file for details).
3. Execute each script in numerical order. Many steps read the outputs produced by earlier scripts.

## Configuration

Update `config.R` to point to your local directories for the clinical metadata,
WGS result files and desired output folder. By default it assumes a structure
like:

```
clinical_data_dir/     # spreadsheets and text files
wgs_results_dir/       # variant calls, CNAs, fragmentomics
output_tables_dir/     # where tables and figures will be written
```

Edit these paths to match your environment before running the scripts.

Before running any script, source the helper files:

```r
source("setup_packages.R")  # loads required packages
source("config.R")           # defines directory locations
source("helpers.R")          # utility functions
```
If any packages are missing, install them with `install.packages(packages)` before sourcing.

Because some of the input spreadsheets contain protected patient information, they are not distributed with this repository. You will need to request them from the study authors.

## Output

The pipeline generates several CSV tables summarising the processed clinical metadata, mutation data, and MRD results, along with PDF/PNG figures. These outputs are used directly in the manuscript figures and supplementary tables.

## Contact

For questions about the analysis please contact Dory Abelman or open an issue on this repository.
