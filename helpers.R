# =============================================================================
# helpers.R
# Project:  cfWGS MRD Detection in Multiple Myeloma
# Author:   Dory Abelman
#
# Purpose:
#   Shared utility functions used across the numbered analysis scripts.
#   Source this file at the start of any script that calls these functions.
#
# Functions defined here:
#   clean_sample_id(x)  - strips whitespace and converts sample IDs to
#                         uppercase for consistent joining across data sources.
#
# Usage:
#   source("helpers.R")  # called at the top of each numbered script
# =============================================================================

# Collection of helper functions used across scripts

clean_sample_id <- function(x) {
  # remove spaces and convert to uppercase
  toupper(gsub("\\s+", "", x))
}
