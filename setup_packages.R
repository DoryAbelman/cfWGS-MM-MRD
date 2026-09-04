# =============================================================================
# setup_packages.R
# Project:  cfWGS MRD Detection in Multiple Myeloma
# Author:   Dory Abelman
#
# Purpose:
#   Optional loader for a commonly used subset of the packages in this project.
#   The numbered scripts load their own dependencies and do not require this
#   file to be sourced first.
#
#   The script checks its package list and stops at the first missing package.
#   It is not a complete dependency installer: additional packages used by
#   specific scripts are listed in config.R and in those scripts' headers.
#
# Usage:
#   source("setup_packages.R")
#
# Installation note:
#   CRAN, Bioconductor, and non-CRAN dependencies require their corresponding
#   installation methods; a single install.packages() call is not sufficient
#   for every package used by the project.
#
# Manuscript outputs created/updated:
#   - None directly. This support file standardizes package loading for
#     reproducible reruns.
# =============================================================================

# Shared package loader
packages <- c(
  # Core tidyverse packages for data manipulation and plotting
  "tidyverse",
  "readxl",
  "lubridate",

  # Genomics / plotting utilities
  "ComplexHeatmap",
  "circlize",
  "maftools",
  "GenomicRanges",

  # Table generation and reporting
  "gtsummary",
  "officer",
  "flextable",
  "gt",

  # Statistical helpers
  "pROC",
  "patchwork",
  "rmda",
  "exact2x2",

  # General utilities
  "janitor",
  "Hmisc",
  "broom",
  "glue",
  "writexl",
  "ggridges",
  "viridis",
  "pbapply",
  "scales"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed", pkg))
  }
  library(pkg, character.only = TRUE)
}
