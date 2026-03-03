# =============================================================================
# setup_packages.R
# Project:  cfWGS MRD Detection in Multiple Myeloma
# Author:   Dory Abelman
#
# Purpose:
#   One-shot package installer and loader for the cfWGS-MM-MRD pipeline.
#   Run this file once per R session (or whenever setting up a new environment)
#   before executing any of the numbered analysis scripts.
#
#   The script checks whether each required package is installed and throws
#   an informative error if any are missing, instead of failing silently mid-
#   analysis. Install missing packages with install.packages() before re-running.
#
# Usage:
#   source("setup_packages.R")
#
# To install all dependencies at once:
#   install.packages(packages)  # where `packages` is the vector defined below
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
