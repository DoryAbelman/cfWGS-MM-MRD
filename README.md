# cfWGS-MM-MRD

This repository contains the analysis scripts used in **Abelman et al. (2025)** for the study *"Cell-free DNA Whole Genome Sequencing for Non‑Invasive MRD Detection in Multiple Myeloma"*.

The code is organised as a series of numbered R scripts that process raw clinical data, integrate whole genome sequencing features, and generate the figures found in the manuscript. The raw data files referenced by these scripts are not included in this repository.

## Repository layout

| Step | Script | Purpose |
| --- | --- | --- |
| 1.0 | `1_0_Process_clinical_metadata.R` | Harmonise clinical metadata from multiple cohorts and export cleaned tables. |
| 1.1 | `1_1_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R` | Tidy transplant dates, relapse information, and overall survival metadata. |
| 1.2 | `1_2_Process_Mutation_Data.R` | Prepare mutation calls for baseline and follow‑up samples. |
| 1.3 | `1_3_Process_Ig_Translocation_Info.R` | Extract immunoglobulin translocation calls from WGS data. |
| 1.4 | `1_4_Process_CNA_Data.R` | Summarise copy number features from WGS results. |
| 1.5 | `1_5_Integrate_WGS_Feature_Data.R` | Combine mutation, CNA, and translocation features. |
| 1.6 | `1_6_Identify_High_Quality_Patient_Pairs.R` | Annotate samples that have matched bone‑marrow and cfDNA data. |
| 1.7A | `1_7A_Process_fragmentomics_data_nucleosome_accessibility.R` | Derive nucleosome accessibility metrics from fragmentomics. |
| 1.7B | `1_7B_Process_fragment_score_and_integrate_fragmentomics_data.R` | Compute fragmentation scores and integrate with the main table. |
| 1.8 | `1_8_Process_Cumulative_VAFs_MRDetect.R` | Generate cumulative VAF plots using MRDetect results. |
| 2.0 | `2_0_Assemble_Table_With_All_Features.R` | Merge all processed features into a single table for downstream analyses. |
| 2.1 | `2_1_Clinical_Demographics_Table.R` | Produce patient demographic tables used in the manuscript. |
| 2.2 | `2_2_Baseline_demographics_by_WGS_heatmap_updated.R` | Create integrated alteration heatmaps. |
| 2.3 | `2_3_Feature_Concordance_And_Mutation_Counts.R` | Calculate concordance between bone‑marrow and cfDNA calls. |
| 3.1 | `3_1_Optimize_cfWGS_thresholds.R` | Determine optimal cfWGS tumour‑fraction thresholds. |
| 3.2 | `3_2_Plot_optimal_cutoff.R` | Plot the ROC curves for threshold selection. |
| 3.3 | `3_3_Plot_optimal_cutoff_tumor_naive_calls.R` | ROC analysis using tumour‑naïve mutation calls. |

Each script is intended to be run individually after placing the appropriate input files in the paths referenced within the code. Figures and intermediate tables are written to the working directory.

A small set of helper scripts is included:
`setup_packages.R` loads all required libraries, `config.R` defines common directory locations, and `helpers.R` provides utility functions such as `clean_sample_id()`.

## Running the analysis

1. Ensure you have an R installation (\>=4.2) with the packages used throughout the scripts. The main packages are **tidyverse**, **readxl**, **lubridate**, **ComplexHeatmap**, **maftools**, and **janitor**.
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
