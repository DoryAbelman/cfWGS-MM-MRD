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
