# =============================================================================
# config.R
# Project:  cfWGS MRD Detection in Multiple Myeloma
# Author:   Dory Abelman
#
# Purpose:
#   Central configuration file for the cfWGS-MM-MRD analysis pipeline.
#   Contains two things:
#     1. PACKAGE LIST  - complete list of every R package required by any
#        numbered analysis script. Sourcing config.R (or setup_packages.R,
#        which also references this list) ensures the full environment is
#        loaded before running the pipeline.
#     2. DIRECTORY PATHS - common path variables used throughout the scripts.
#        Update the variables below to match your local environment before
#        running any analysis script.
#
# Usage:
#   source("config.R")   # called at the top of each numbered script
#
# Directory variables (update these for your system):
#   clinical_data_dir  <- path to clinical metadata spreadsheets
#   wgs_results_dir    <- path to WGS variant call, CNA, and fragmentomics files
#   output_tables_dir  <- path where processed tables and figures are written
# =============================================================================

# Packages required across the analysis scripts.
packages <- c(
  "BoutrosLab.plotting.general",
  "ChromHeatMap",
  "ComplexHeatmap",
  "DescTools",
  "GGally",
  "GeneCycle",
  "GenomicRanges",
  "Hmisc",
  "Matrix",
  "PRROC",
  "RColorBrewer",
  "VennDiagram",
  "broom",
  "caret",
  "circlize",
  "conflicted",
  "data.table",
  "doParallel",
  "dplyr",
  "exact2x2",
  "flextable",
  "forcats",
  "fuzzyjoin",
  "ggbreak",
  "ggplot2",
  "ggpubr",
  "ggridges",
  "glmnet",
  "glue",
  "gridExtra",
  "gt",
  "gtsummary",
  "janitor",
  "lubridate",
  "maftools",
  "officer",
  "openxlsx",
  "pROC",
  "patchwork",
  "pbapply",
  "purrr",
  "readr",
  "readxl",
  "reshape2",
  "rlang",
  "rmda",
  "rstatix",
  "scales",
  "stringdist",
  "stringr",
  "survival",
  "survminer",
  "tableone",
  "tibble",
  "tidyr",
  "tidyverse",
  "timeROC",
  "viridis",
  "viridisLite",
  "writexl",
  "yardstick"
)

# Central configuration of directory locations.
# Update these variables to match your local environment.

clinical_data_dir <- "Clinical data"
wgs_results_dir   <- "M4_CMRG_Data"
output_tables_dir <- "Output_tables_2025"
