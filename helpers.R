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
#   spring2026_*()      - revision-cohort helpers. These are deliberately kept
#                         in one file so every downstream script uses the same
#                         metadata authority, exclusion logic, alias rules, and
#                         dilution-series parsing assumptions.
#
# Usage:
#   source("helpers.R")  # called at the top of each numbered script
#
# Manuscript outputs created/updated:
#   - None directly. This support file centralizes reusable helper functions
#     for consistent sample-ID handling.
# =============================================================================

# Collection of helper functions used across scripts

clean_sample_id <- function(x) {
  # remove spaces and convert to uppercase
  toupper(gsub("\\s+", "", x))
}

if (!exists("%>%", mode = "function")) {
  `%>%` <- magrittr::`%>%`
}

project_file <- function(...) {
  file.path(...)
}

spring2026_revision_data_dir <- function() {
  # ## Spring 2026 raw/derived data root
  # Default to the repo-local revision data folder, but allow the same script to
  # be rerun on a cluster or external checkout by setting
  # CFWGS_SPRING2026_REVISION_DATA_DIR. Keeping this as a function avoids
  # hard-coding absolute user-specific paths throughout the final scripts.
  Sys.getenv(
    "CFWGS_SPRING2026_REVISION_DATA_DIR",
    unset = "Data_Spring_2026_Revisions"
  )
}

spring2026_revision_metadata_path <- function() {
  # ## Authoritative Spring 2026 metadata table
  # This CSV is the repo-shaped OICR revision table, not a filename-derived
  # reconstruction. It is the source of truth for patient IDs, sample IDs,
  # sample types, dates, timepoint labels, study labels, BAM names, and whether a
  # row is ready to append to the combined clinical metadata.
  Sys.getenv(
    "CFWGS_SPRING2026_REVISION_METADATA",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_repo_style_metadata.csv"
    )
  )
}

spring2026_revision_timepoint_override_path <- function() {
  # ## Manual timepoint overrides
  # Overrides are separated from the main metadata table so reviewer/manual
  # adjudications can be audited independently. They are matched by Sample_ID and
  # only replace timepoint_info after duplicate and unmatched IDs are rejected.
  Sys.getenv(
    "CFWGS_SPRING2026_TIMEPOINT_OVERRIDES",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_timepoint_overrides.csv"
    )
  )
}

spring2026_revision_primary_analysis_exclusion_path <- function() {
  # ## Primary-analysis exclusions
  # This optional table removes revision samples from automatic downstream
  # integration while preserving the reason in a separate auditable file. It
  # prevents ad hoc filtering rules from being hidden inside analysis scripts.
  Sys.getenv(
    "CFWGS_SPRING2026_PRIMARY_ANALYSIS_EXCLUSIONS",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_primary_analysis_exclusions.csv"
    )
  )
}

spring2026_revision_enabled <- function() {
  # ## Global revision toggle
  # The default is ON because these scripts now represent the revision-aware
  # manuscript workflow. Setting CFWGS_USE_SPRING2026_REVISION=0 restores the
  # historical inputs for debugging without editing code.
  identical(Sys.getenv("CFWGS_USE_SPRING2026_REVISION", unset = "1"), "1")
}

require_columns <- function(data, required, label) {
  # Fail early when a metadata/export contract changes. Silent NA-filled joins
  # are especially dangerous here because they can move revision rows between
  # cohorts, sample types, or timepoints without an obvious error.
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

xplus_charm_healthy_bams <- function() {
  # ## Exact XPlus CHARM healthy-control BAM allowlist
  # These TGL49 BAMs are healthy controls sequenced in the XPlus/CHARM context.
  # They must be identified by exact basename rather than by TGL49 prefix alone,
  # because other TGL49-like rows can be patient/dilution data and should not be
  # folded into the healthy-control reference distribution.
  c(
    "TGL49_0197_Cf_R_WG_HBC-01-0002-T2-P-DNA_NovaX.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0199_Cf_R_WG_HBC-01-0004-T2-P-DNA_NovaX.filter.deduped.recalibrated.0.3.bam",
    "TGL49_0255_Ct_n_WG.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0256_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0257_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0259_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0260_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0261_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0262_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0263_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0264_Pb_U_WG_HBC-01-0018-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0265_Pb_U_WG_HBC-01-0019-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.0.7.bam",
    "TGL49_0266_Pb_U_WG_HBC-01-0020-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0267_Cf_R_WG_HBC-01-0008-T1_NovaX.filter.deduped.recalibrated.0.55.bam",
    "TGL49_0267_Pb_U_WG_HBC-01-0021-T0-P-DNA-1.filter.deduped.recalibrated.0.7.bam",
    "TGL49_0268_Pb_U_WG_HBC-01-0022-T0-P-DNA-1.filter.deduped.recalibrated.bam",
    "TGL49_0302_Cf_R_WG_HBC-01-0023-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0303_Cf_R_WG_HBC-01-0024-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0304_Cf_R_WG_HBC-01-0025-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0305_Cf_R_WG_HBC-01-0026-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0306_Cf_R_WG_HBC-01-0027-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0307_Cf_R_WG_HBC-01-0028-T0-P-DNA.filter.deduped.recalibrated.bam"
  )
}

is_xplus_charm_healthy_bam <- function(x) {
  # Compare basenames so full paths, relative paths, and bare BAM filenames all
  # classify consistently.
  basename(as.character(x)) %in% xplus_charm_healthy_bams()
}

load_spring2026_revision_metadata <- function(required = FALSE) {
  # ## Load and validate revision metadata
  # Returns NULL when revision integration is disabled or optional and absent.
  # When a caller sets required=TRUE, absence becomes a hard failure because the
  # downstream analysis cannot be interpreted without this metadata authority.
  path <- spring2026_revision_metadata_path()
  if (!spring2026_revision_enabled()) return(NULL)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing Spring 2026 revision metadata: ", path, call. = FALSE)
    return(NULL)
  }
  metadata <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(
    metadata,
    c("Bam", "Patient", "Date_of_sample_collection", "Sample_type", "Timepoint",
      "Study", "Sample_ID", "timepoint_info", "integration_ready_for_combined_clinical_data"),
    "Spring 2026 revision metadata"
  )
  not_ready <- metadata$integration_ready_for_combined_clinical_data %in% FALSE
  if (any(not_ready, na.rm = TRUE)) {
    # Do not let pending-review rows leak into model/test-cohort outputs. The
    # handoff pipeline can carry unresolved rows, but final scripts only append
    # rows explicitly marked ready for combined clinical integration.
    stop(
      "Spring 2026 revision metadata has rows not ready for combined clinical integration: ",
      paste(metadata$Sample_ID[not_ready], collapse = ", "),
      call. = FALSE
    )
  }
  override_path <- spring2026_revision_timepoint_override_path()
  if (file.exists(override_path)) {
    # Apply manual timepoint corrections only after checking that the override
    # table is one row per Sample_ID and refers only to samples in the revision
    # metadata. This keeps the override file reviewable and prevents accidental
    # creation of new samples by typo.
    overrides <- utils::read.csv(override_path, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(
      overrides,
      c("Sample_ID", "timepoint_info"),
      "Spring 2026 revision timepoint override table"
    )
    duplicate_override_sample_id <- overrides$Sample_ID[
      !is.na(overrides$Sample_ID) & duplicated(overrides$Sample_ID)
    ]
    if (length(duplicate_override_sample_id)) {
      stop(
        "Spring 2026 revision timepoint override table has duplicate Sample_ID values: ",
        paste(unique(duplicate_override_sample_id), collapse = ", "),
        call. = FALSE
      )
    }
    unmatched_override_sample_id <- setdiff(overrides$Sample_ID, metadata$Sample_ID)
    if (length(unmatched_override_sample_id)) {
      stop(
        "Spring 2026 revision timepoint override table has Sample_ID values absent from revision metadata: ",
        paste(unmatched_override_sample_id, collapse = ", "),
        call. = FALSE
      )
    }
    idx <- match(overrides$Sample_ID, metadata$Sample_ID)
    metadata$timepoint_info_original_before_override <- metadata$timepoint_info
    metadata$timepoint_info[idx] <- overrides$timepoint_info
    metadata$timepoint_info_override_reason <- NA_character_
    if ("override_reason" %in% names(overrides)) {
      metadata$timepoint_info_override_reason[idx] <- overrides$override_reason
    }
  }
  exclusion_path <- spring2026_revision_primary_analysis_exclusion_path()
  if (file.exists(exclusion_path)) {
    # Apply sample-level exclusions after overrides. The exclusion reason remains
    # in the external CSV, while this loader returns only rows eligible for the
    # current primary-analysis workflow.
    exclusions <- utils::read.csv(exclusion_path, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(
      exclusions,
      c("Sample_ID", "exclusion_reason"),
      "Spring 2026 primary-analysis exclusion table"
    )
    duplicate_exclusion_sample_id <- exclusions$Sample_ID[
      !is.na(exclusions$Sample_ID) & duplicated(exclusions$Sample_ID)
    ]
    if (length(duplicate_exclusion_sample_id)) {
      stop(
        "Spring 2026 primary-analysis exclusion table has duplicate Sample_ID values: ",
        paste(unique(duplicate_exclusion_sample_id), collapse = ", "),
        call. = FALSE
      )
    }
    unmatched_exclusion_sample_id <- setdiff(exclusions$Sample_ID, metadata$Sample_ID)
    if (length(unmatched_exclusion_sample_id)) {
      stop(
        "Spring 2026 primary-analysis exclusion table has Sample_ID values absent from revision metadata: ",
        paste(unmatched_exclusion_sample_id, collapse = ", "),
        call. = FALSE
      )
    }
    metadata <- metadata[
      !(metadata$Sample_ID %in% exclusions$Sample_ID),
      ,
      drop = FALSE
    ]
  }
  metadata
}

augment_cohort_assignment_with_spring2026_revision <- function(cohort_df) {
  # ## Add revision patients to the cohort assignment table
  # Spring 2026 rows are treated as a test-cohort expansion. They should not
  # silently redefine the original training cohort, so new patients are appended
  # as Non-frontline while existing cohort labels are left unchanged.
  revision <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(revision)) return(cohort_df)

  require_columns(cohort_df, c("Patient", "Cohort"), "Cohort assignment table")

  revision_cohort <- revision %>%
    dplyr::filter(!is.na(.data$Patient), nzchar(.data$Patient)) %>%
    dplyr::distinct(Patient = as.character(.data$Patient)) %>%
    dplyr::mutate(Cohort = "Non-frontline")

  cohort_df %>%
    dplyr::mutate(Patient = as.character(.data$Patient)) %>%
    dplyr::bind_rows(
      revision_cohort %>%
        dplyr::filter(!.data$Patient %in% cohort_df$Patient)
    ) %>%
    dplyr::distinct(.data$Patient, .keep_all = TRUE)
}

coerce_revision_column_like_current <- function(current_col, revision_col, column_name) {
  # ## Type harmonization before bind_rows()
  # The revision CSV is read from disk and can arrive as character columns even
  # when the historical metadata column is Date, numeric, integer, or logical.
  # Coercing the revision column to the historical column's type preserves
  # downstream expectations and avoids accidental list/character promotion.
  if (grepl("date", column_name, ignore.case = TRUE)) {
    return(parse_date_safely(revision_col))
  }
  if (inherits(current_col, "Date")) {
    return(parse_date_safely(revision_col))
  }
  if (inherits(current_col, "POSIXct") || inherits(current_col, "POSIXt")) {
    return(as.POSIXct(revision_col, tz = attr(current_col, "tzone", exact = TRUE) %||% "UTC"))
  }
  if (is.character(current_col)) {
    return(as.character(revision_col))
  }
  if (is.numeric(current_col) && !is.logical(current_col)) {
    return(suppressWarnings(as.numeric(revision_col)))
  }
  if (is.integer(current_col)) {
    return(suppressWarnings(as.integer(revision_col)))
  }
  if (is.logical(current_col)) {
    if (all(is.na(revision_col) | revision_col %in% c(TRUE, FALSE, "TRUE", "FALSE", "true", "false", "1", "0"))) {
      return(as.logical(revision_col))
    }
    return(as.character(revision_col))
  }
  revision_col
}

parse_date_safely <- function(x) {
  # Accept the ISO dates emitted by the handoff pipeline plus common manual Excel
  # date formats. Unparseable non-empty values become NA rather than guessing.
  if (inherits(x, "Date")) return(x)
  x <- as.character(x)
  x[!nzchar(x)] <- NA_character_
  out <- suppressWarnings(tryCatch(
    as.Date(x),
    error = function(e) rep(as.Date(NA), length(x))
  ))
  unresolved <- is.na(out) & !is.na(x)
  if (any(unresolved)) {
    formats <- c("%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d", "%m-%d-%Y", "%d-%m-%Y")
    for (fmt in formats) {
      idx <- unresolved
      parsed <- suppressWarnings(as.Date(x[idx], format = fmt))
      idx_positions <- which(idx)
      out[idx_positions[!is.na(parsed)]] <- parsed[!is.na(parsed)]
      unresolved <- is.na(out) & !is.na(x)
      if (!any(unresolved)) break
    }
  }
  out
}

`%||%` <- function(x, y) {
  # Null-coalescing helper used for optional timezone attributes and similar
  # fields that may be absent.
  if (is.null(x) || length(x) == 0L) y else x
}

spring2026_normalize_patient_alias <- function(x) {
  # ## Alias normalization for overlap audits only
  # OICR/workbook rows may use MyP, MYP, IMG, hyphenated, underscored, or
  # zero-padded variants for the same patient. This function normalizes those
  # aliases to a patient-level audit key. It must not be used to rewrite
  # displayed Sample_ID, Bam, or manuscript labels.
  y <- toupper(trimws(as.character(x)))
  y[is.na(x)] <- NA_character_
  # MyP and IMG identifiers can name the same patient/sample in different
  # input systems. Normalize only patient-level overlap keys; do not rewrite
  # displayed Sample_ID/Bam values.
  y <- sub("^(MYP|IMG)[-_]?0*([0-9]+)$", "IMG-\\2", y)
  y
}

repair_historical_combined_clinical_metadata <- function(
    current,
    legacy_path = "combined_clinical_data_updated_Feb5_2025.csv") {
  # ## Restore fields accidentally blanked in the April metadata handoff
  # The April combined clinical table retained several historical SPORE/MyC rows
  # but lost sample IDs and collection dates needed to join clinical MRD results.
  # Fill only missing current fields from an older combined clinical metadata
  # source keyed by the same biological sample context. This is metadata repair,
  # not restoration from model/scoring outputs.
  if (!file.exists(legacy_path)) return(current)

  key_cols <- intersect(c("Patient", "Timepoint", "Sample_type", "timepoint_info"), names(current))
  fill_cols <- intersect(c("Date_of_sample_collection", "Sample_ID", "Study", "Bam"), names(current))
  if (!length(key_cols) || !length(fill_cols)) return(current)

  legacy <- readr::read_csv(legacy_path, show_col_types = FALSE)
  key_cols <- intersect(key_cols, names(legacy))
  fill_cols <- intersect(fill_cols, names(legacy))
  if (!length(key_cols) || !length(fill_cols)) return(current)

  first_non_missing <- function(x) {
    if (inherits(x, "Date")) {
      vals <- x[!is.na(x)]
      if (length(vals)) vals[[1]] else as.Date(NA)
    } else {
      vals <- x[!is.na(x) & nzchar(trimws(as.character(x)))]
      if (length(vals)) vals[[1]] else NA
    }
  }

  normalize_key_fields <- function(dat) {
    dat %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(key_cols), ~ trimws(as.character(.x))),
        dplyr::across(dplyr::any_of("Date_of_sample_collection"), parse_date_safely)
      )
  }

  current_repair_input <- normalize_key_fields(current)
  legacy_map <- normalize_key_fields(legacy) %>%
    dplyr::select(dplyr::all_of(c(key_cols, fill_cols))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(fill_cols), first_non_missing), .groups = "drop")

  joined <- current_repair_input %>%
    dplyr::left_join(legacy_map, by = key_cols, suffix = c("", ".legacy"))

  repair_matrix <- as.data.frame(
    lapply(fill_cols, function(col) {
      legacy_col <- paste0(col, ".legacy")
      if (!legacy_col %in% names(joined)) return(rep(FALSE, nrow(joined)))
      is.na(joined[[col]]) & !is.na(joined[[legacy_col]])
    })
  )
  names(repair_matrix) <- fill_cols

  repaired <- joined
  for (col in fill_cols) {
    legacy_col <- paste0(col, ".legacy")
    if (!legacy_col %in% names(repaired)) next
    repaired[[col]] <- dplyr::coalesce(repaired[[col]], repaired[[legacy_col]])
  }

  audit <- repaired %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(key_cols)),
      repaired_fields = apply(repair_matrix, 1, function(x) paste(names(x)[x], collapse = ";"))
    ) %>%
    dplyr::filter(nzchar(.data$repaired_fields)) %>%
    dplyr::distinct()

  if (nrow(audit) > 0L) {
    audit_dir <- file.path("Output_tables_2025", "clinical_support")
    if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(
      audit,
      file.path(audit_dir, "historical_combined_clinical_metadata_repair_audit.csv")
    )
  }

  repaired %>%
    dplyr::select(-dplyr::ends_with(".legacy"))
}

build_cfwgs_sample_identity_map <- function(
    canonical_metadata,
    legacy_paths = c(
      "combined_clinical_data_updated_April2025.csv",
      "combined_clinical_data_updated_Feb5_2025.csv"
    ),
    audit_path = file.path(
      "Output_tables_2025",
      "clinical_support",
      "cfwgs_sample_identity_map.csv"
    )) {
  # ## Central sample identity map for metadata joins
  # Revision rows may replace historical identifiers while older MRDetect files
  # still contain the historical BAM/VCF basenames. This map preserves those
  # source-level aliases and points them to the current canonical metadata row.
  # It is used only for joining raw inputs to metadata, not for final-table rescue.
  required <- c(
    "Bam", "Patient", "Date_of_sample_collection", "Sample_type",
    "Timepoint", "Study", "Sample_ID", "timepoint_info"
  )
  missing_required <- setdiff(required, names(canonical_metadata))
  if (length(missing_required) > 0L) {
    stop(
      "Cannot build cfWGS sample identity map; canonical metadata is missing: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  normalize_for_identity <- function(dat, source_label) {
    for (col in required) {
      if (!col %in% names(dat)) dat[[col]] <- NA
    }
    if (!"mutect2_pair_id" %in% names(dat)) dat$mutect2_pair_id <- NA_character_

    dat %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(required), as.character),
        Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection),
        VCF_clean_merge = sub("[.]filter.*", "", .data$Bam),
        identity_norm_patient = spring2026_normalize_patient_alias(.data$Patient),
        identity_primary_key = dplyr::if_else(
          !is.na(.data$identity_norm_patient) &
            !is.na(.data$Date_of_sample_collection) &
            !is.na(.data$Sample_type),
          paste(
            .data$identity_norm_patient,
            .data$Date_of_sample_collection,
            .data$Sample_type,
            sep = "|"
          ),
          NA_character_
        ),
        identity_timepoint_type_key = dplyr::if_else(
          !is.na(.data$identity_norm_patient) &
            !is.na(.data$Timepoint) &
            !is.na(.data$Sample_type),
          paste(
            .data$identity_norm_patient,
            .data$Timepoint,
            .data$Sample_type,
            sep = "|"
          ),
          NA_character_
        ),
        identity_source_table = source_label,
        identity_row_id = dplyr::row_number()
      )
  }

  canonical <- normalize_for_identity(canonical_metadata, "canonical")
  canonical_cols <- names(canonical)

  read_legacy_metadata <- function(path) {
    if (!file.exists(path)) return(NULL)
    dat <- readr::read_csv(path, show_col_types = FALSE)
    dat <- repair_historical_combined_clinical_metadata(dat)
    normalize_for_identity(dat, basename(path))
  }

  legacy <- dplyr::bind_rows(lapply(legacy_paths, read_legacy_metadata))
  for (col in setdiff(canonical_cols, names(legacy))) {
    legacy[[col]] <- NA
  }
  for (col in setdiff(names(legacy), canonical_cols)) {
    canonical[[col]] <- NA
  }
  legacy <- legacy[, union(canonical_cols, names(legacy)), drop = FALSE]
  canonical <- canonical[, union(canonical_cols, names(legacy)), drop = FALSE]
  canonical_cols <- names(canonical)

  canonical_primary_lookup <- canonical %>%
    dplyr::filter(!is.na(.data$identity_primary_key)) %>%
    dplyr::add_count(.data$identity_primary_key, name = "identity_key_n") %>%
    dplyr::filter(.data$identity_key_n == 1L) %>%
    dplyr::select(.data$identity_primary_key, canonical_identity_row_id = .data$identity_row_id)

  canonical_timepoint_lookup <- canonical %>%
    dplyr::filter(!is.na(.data$identity_timepoint_type_key)) %>%
    dplyr::add_count(.data$identity_timepoint_type_key, name = "identity_key_n") %>%
    dplyr::filter(.data$identity_key_n == 1L) %>%
    dplyr::select(
      .data$identity_timepoint_type_key,
      canonical_identity_row_id_timepoint = .data$identity_row_id
    )

  alias_links <- legacy %>%
    dplyr::select(
      alias_identity_row_id = .data$identity_row_id,
      alias_source_table = .data$identity_source_table,
      alias_Bam = .data$Bam,
      alias_VCF_clean_merge = .data$VCF_clean_merge,
      alias_Sample_ID = .data$Sample_ID,
      alias_Patient = .data$Patient,
      alias_Timepoint = .data$Timepoint,
      alias_Sample_type = .data$Sample_type,
      alias_timepoint_info = .data$timepoint_info,
      alias_Date_of_sample_collection = .data$Date_of_sample_collection,
      .data$identity_primary_key,
      .data$identity_timepoint_type_key
    ) %>%
    dplyr::left_join(canonical_primary_lookup, by = "identity_primary_key") %>%
    dplyr::left_join(canonical_timepoint_lookup, by = "identity_timepoint_type_key") %>%
    dplyr::mutate(
      canonical_identity_row_id = dplyr::coalesce(
        .data$canonical_identity_row_id,
        .data$canonical_identity_row_id_timepoint
      ),
      metadata_join_key_type = dplyr::case_when(
        !is.na(.data$canonical_identity_row_id) &
          !is.na(.data$identity_primary_key) ~ "patient_date_sample_type",
        !is.na(.data$canonical_identity_row_id) ~ "patient_timepoint_sample_type",
        TRUE ~ "legacy_self"
      )
    ) %>%
    dplyr::select(-.data$canonical_identity_row_id_timepoint)

  alias_map <- alias_links %>%
    dplyr::left_join(
      canonical %>%
        dplyr::select(
          canonical_identity_row_id = .data$identity_row_id,
          dplyr::all_of(canonical_cols)
        ),
      by = "canonical_identity_row_id"
    )

  legacy_self <- legacy %>%
    dplyr::filter(.data$identity_row_id %in% alias_links$alias_identity_row_id[is.na(alias_links$canonical_identity_row_id)])

  if (nrow(legacy_self) > 0L) {
    self_links <- alias_links %>%
      dplyr::filter(is.na(.data$canonical_identity_row_id)) %>%
      dplyr::select(
        .data$alias_identity_row_id, .data$alias_source_table, .data$alias_Bam,
        .data$alias_VCF_clean_merge, .data$alias_Sample_ID, .data$alias_Patient,
        .data$alias_Timepoint, .data$alias_Sample_type,
        .data$alias_timepoint_info, .data$alias_Date_of_sample_collection,
        .data$metadata_join_key_type
      )
    alias_map_self <- self_links %>%
      dplyr::left_join(
        legacy_self %>%
          dplyr::select(alias_identity_row_id = .data$identity_row_id, dplyr::all_of(canonical_cols)),
        by = "alias_identity_row_id"
      )
    alias_map <- dplyr::bind_rows(
      alias_map %>% dplyr::filter(!is.na(.data$canonical_identity_row_id)),
      alias_map_self
    )
  } else {
    alias_map <- alias_map %>% dplyr::filter(!is.na(.data$canonical_identity_row_id))
  }

  canonical_map <- canonical %>%
    dplyr::mutate(
      metadata_join_bam = .data$Bam,
      metadata_join_vcf_clean = .data$VCF_clean_merge,
      metadata_identity_source = "canonical",
      metadata_alias_sample_id = .data$Sample_ID,
      metadata_alias_patient = .data$Patient,
      metadata_join_key_type = "canonical",
      metadata_join_priority = 1L
    )

  alias_map <- alias_map %>%
    dplyr::mutate(
      metadata_join_bam = .data$alias_Bam,
      metadata_join_vcf_clean = .data$alias_VCF_clean_merge,
      metadata_identity_source = paste0("alias:", .data$alias_source_table),
      metadata_alias_sample_id = .data$alias_Sample_ID,
      metadata_alias_patient = .data$alias_Patient,
      metadata_join_priority = dplyr::if_else(
        .data$metadata_join_key_type == "legacy_self",
        3L,
        2L
      )
    )

  identity_map <- dplyr::bind_rows(canonical_map, alias_map) %>%
    dplyr::filter(
      (!is.na(.data$metadata_join_bam) & nzchar(.data$metadata_join_bam)) |
        (!is.na(.data$metadata_join_vcf_clean) & nzchar(.data$metadata_join_vcf_clean))
    ) %>%
    dplyr::arrange(.data$metadata_join_priority, .data$Patient, .data$Sample_ID, .data$metadata_join_bam) %>%
    dplyr::distinct(
      .data$metadata_join_bam,
      .data$metadata_join_vcf_clean,
      .data$Sample_ID,
      .keep_all = TRUE
    )

  audit_dir <- dirname(audit_path)
  if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(identity_map, audit_path)

  identity_map
}

read_combined_clinical_metadata_with_revision <- function(
    path = "combined_clinical_data_updated_April2025.csv",
    include_revision_extra = FALSE) {
  # ## Append revision metadata to the historical combined clinical table
  # This is the main integration helper used by final scripts. It keeps the
  # historical table schema by default, optionally retains revision-only columns,
  # validates duplicate identifiers, writes overlap audits, and prevents older
  # rows from coexisting with newer revision rows that represent the same
  # patient/date/sample-type or patient/timepoint/sample-type observation.
  if (!file.exists(path)) stop("Missing combined clinical metadata: ", path, call. = FALSE)
  current <- readr::read_csv(path, show_col_types = FALSE)
  current <- repair_historical_combined_clinical_metadata(current)
  revision <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(revision)) return(current)

  current_original <- current
  revision_original <- revision

  shared_cols <- if (isTRUE(include_revision_extra)) union(names(current), names(revision)) else names(current)
  extra_cols <- setdiff(names(revision), shared_cols)
  # Align schemas before binding. Missing columns are explicit NA rather than
  # being dropped silently; this matters for downstream joins that expect a
  # stable combined clinical table.
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  for (col in intersect(shared_cols, intersect(names(current), names(revision)))) {
    revision[[col]] <- coerce_revision_column_like_current(current[[col]], revision[[col]], col)
    if (grepl("date", col, ignore.case = TRUE) || inherits(revision[[col]], "Date")) {
      current[[col]] <- parse_date_safely(current[[col]])
    }
    if (is.character(revision[[col]]) && is.logical(current[[col]])) {
      current[[col]] <- as.character(current[[col]])
    }
  }
  revision_core <- revision[, shared_cols, drop = FALSE]
  current_core <- current[, shared_cols, drop = FALSE]

  duplicate_revision_sample_id <- revision_original$Sample_ID[
    !is.na(revision_original$Sample_ID) & duplicated(revision_original$Sample_ID)
  ]
  duplicate_revision_bam <- revision_original$Bam[
    !is.na(revision_original$Bam) & duplicated(revision_original$Bam)
  ]
  current_original <- current_original %>%
    dplyr::mutate(Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection))
  revision_original <- revision_original %>%
    dplyr::mutate(Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection))

  cross_sample_id <- intersect(
    # Exact Sample_ID/Bam overlaps are hard collisions and must be audited.
    unique(stats::na.omit(as.character(current_original$Sample_ID))),
    unique(stats::na.omit(as.character(revision_original$Sample_ID)))
  )
  cross_bam <- intersect(
    unique(stats::na.omit(as.character(current_original$Bam))),
    unique(stats::na.omit(as.character(revision_original$Bam)))
  )
  current_alias_keys <- current_original %>%
    # Alias/date/sample-type keys catch biologically same rows when the patient
    # label differs between IMG/MyP systems but collection date and compartment
    # agree.
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Date_of_sample_collection), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      alias_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Date_of_sample_collection,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$alias_key) %>%
    unique()
  revision_alias_keys <- revision_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Date_of_sample_collection), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      alias_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Date_of_sample_collection,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$alias_key) %>%
    unique()
  cross_alias_key <- intersect(current_alias_keys, revision_alias_keys)

  current_timepoint_type_keys <- current_original %>%
    # Patient/timepoint/sample-type keys catch same-context rows when dates are
    # missing or were revised between handoffs. They are deliberately used only
    # for overlap handling, not for deriving new labels.
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Timepoint), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      timepoint_type_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Timepoint,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$timepoint_type_key) %>%
    unique()
  revision_timepoint_type_keys <- revision_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Timepoint), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      timepoint_type_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Timepoint,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$timepoint_type_key) %>%
    unique()
  cross_timepoint_type_key <- intersect(current_timepoint_type_keys, revision_timepoint_type_keys)

  if (length(duplicate_revision_sample_id) || length(duplicate_revision_bam)) {
    # Internal duplicate revision IDs would make automatic binding ambiguous, so
    # stop before any downstream table can inherit the problem.
    stop(
      "Spring 2026 revision metadata has internal duplicate identifiers that need adjudication before automatic append.",
      "\nDuplicate revision Sample_ID: ", paste(unique(duplicate_revision_sample_id), collapse = ", "),
      "\nDuplicate revision Bam: ", paste(unique(duplicate_revision_bam), collapse = ", "),
      call. = FALSE
    )
  }

  if (length(cross_sample_id) || length(cross_bam) || length(cross_alias_key) || length(cross_timepoint_type_key)) {
    # Preserve evidence for every historical/revision overlap decision. The code
    # writes both sides of each overlap before removing the historical row from
    # the combined table, making the replacement auditable.
    audit_dir <- "Output_tables_2025"
    if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
    overlap_audit <- dplyr::bind_rows(
      current_original %>%
        dplyr::mutate(
          spring2026_alias_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Date_of_sample_collection,
            .data$Sample_type,
            sep = "|"
          ),
          spring2026_timepoint_type_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Timepoint,
            .data$Sample_type,
            sep = "|"
          )
        ) %>%
        dplyr::filter(
          Sample_ID %in% cross_sample_id |
            Bam %in% cross_bam |
            spring2026_alias_key %in% cross_alias_key |
            spring2026_timepoint_type_key %in% cross_timepoint_type_key
        ) %>%
        dplyr::mutate(metadata_source_for_overlap_audit = "current") %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character)),
      revision_original %>%
        dplyr::mutate(
          spring2026_alias_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Date_of_sample_collection,
            .data$Sample_type,
            sep = "|"
          ),
          spring2026_timepoint_type_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Timepoint,
            .data$Sample_type,
            sep = "|"
          )
        ) %>%
        dplyr::filter(
          Sample_ID %in% cross_sample_id |
            Bam %in% cross_bam |
            spring2026_alias_key %in% cross_alias_key |
            spring2026_timepoint_type_key %in% cross_timepoint_type_key
        ) %>%
        dplyr::mutate(metadata_source_for_overlap_audit = "spring2026_revision") %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    )
    readr::write_csv(
      overlap_audit,
      file.path(audit_dir, "spring2026_current_revision_identifier_overlap_audit.csv")
    )

    current_core <- current_core %>%
      # Prefer revision rows over older current rows when an overlap is detected,
      # because the Spring 2026 table is the current manually curated metadata
      # authority for those samples.
      dplyr::mutate(
        spring2026_alias_key = paste(
          spring2026_normalize_patient_alias(.data$Patient),
          .data$Date_of_sample_collection,
          .data$Sample_type,
          sep = "|"
        ),
        spring2026_timepoint_type_key = paste(
          spring2026_normalize_patient_alias(.data$Patient),
          .data$Timepoint,
          .data$Sample_type,
          sep = "|"
        )
      ) %>%
      dplyr::filter(
        !(Sample_ID %in% cross_sample_id),
        !(Bam %in% cross_bam),
        !(spring2026_alias_key %in% cross_alias_key),
        !(spring2026_timepoint_type_key %in% cross_timepoint_type_key)
      ) %>%
      dplyr::select(-spring2026_alias_key, -spring2026_timepoint_type_key)
  }

  combined <- dplyr::bind_rows(current_core, revision_core)

  if (length(extra_cols)) {
    # Store the revision-only column names as an attribute for scripts that need
    # to inspect whether extra handoff/provenance fields were present.
    attr(combined, "spring2026_revision_extra_columns") <- extra_cols
  }
  combined
}

spring2026_revision_maf_files <- function(sample_types) {
  # Return Spring 2026 MAFs whose mutect2_pair_id belongs to the requested sample
  # types. This uses metadata pair IDs rather than listing every MAF blindly, so
  # BM-only and blood-only mutation-processing scripts do not cross-contaminate.
  metadata <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(metadata)) return(character())
  maf_dir <- file.path(spring2026_revision_data_dir(), "MuTect2_All_Mafs")
  if (!dir.exists(maf_dir)) return(character())

  wanted_pairs <- unique(metadata$mutect2_pair_id[
    metadata$Sample_type %in% sample_types &
      !is.na(metadata$mutect2_pair_id) &
      nzchar(metadata$mutect2_pair_id)
  ])
  files <- list.files(maf_dir, pattern = "[.]maf$", full.names = TRUE)
  if (!length(files) || !length(wanted_pairs)) return(character())
  keep <- vapply(files, function(path) {
    any(startsWith(basename(path), paste0(wanted_pairs, ".")))
  }, logical(1))
  sort(files[keep])
}

spring2026_revision_files <- function(subdir, pattern) {
  # Small wrapper around list.files() that respects the revision root and returns
  # a deterministic sorted vector. Missing optional folders return character(0).
  root <- file.path(spring2026_revision_data_dir(), subdir)
  if (!dir.exists(root)) return(character())
  sort(list.files(root, pattern = pattern, full.names = TRUE, recursive = FALSE))
}

spring2026_pwgval_dilution_metadata_path <- function() {
  # Repo-local derived table for the PWGVAL/M4CHIP dilution samples. It is used
  # by both MRDetect dilution processing and downstream LOD model application.
  file.path(
    spring2026_revision_data_dir(),
    "derived_metadata",
    "pwgval_dilution_metadata.csv"
  )
}

spring2026_pwgval_dilution_tumor_fraction_path <- function() {
  # Optional cached tumor-fraction table for PWGVAL/M4CHIP dilution rows. If it
  # is absent, the helper rebuilds it from ichorCNA params when available.
  file.path(
    spring2026_revision_data_dir(),
    "derived_metadata",
    "pwgval_dilution_tumor_fraction.tsv"
  )
}

spring2026_pwgval_dilution_workbook_path <- function() {
  # The workbook is a fallback metadata source for dilution rows when the derived
  # CSV has not been materialized yet. Environment override supports cluster
  # layouts without changing the final script code.
  Sys.getenv(
    "CFWGS_SPRING2026_PWGVAL_DILUTION_WORKBOOK",
    unset = file.path(
      "Preparing shell scripts for cluster",
      "Dilution amounts for M4CHIP libs.xlsx"
    )
  )
}

build_spring2026_pwgval_dilution_metadata_from_workbook <- function(path = spring2026_pwgval_dilution_workbook_path()) {
  # ## Parse the PWGVAL/M4CHIP dilution workbook
  # Assumption: workbook Alias values begin with M4CHIP_#### and Group_ID encodes
  # the patient and replicate suffix, e.g. VA-09-00001-01. The constructed BAM
  # name must match the MRDetect queried BAM basename so joins are exact.
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "Package readxl is required to read Spring 2026 PWGVAL dilution workbook: ",
      path,
      call. = FALSE
    )
  }
  workbook <- readxl::read_excel(path, sheet = 1)
  names(workbook) <- gsub("[^A-Za-z0-9]+", "_", names(workbook))
  names(workbook) <- gsub("_$", "", names(workbook))
  require_columns(
    workbook,
    c("Name", "Alias", "Group_ID", "TF_Dilution", "Tissue_Origin", "Tissue_Type",
      "Library_Design", "Index_es", "i7_Index", "i5_Index", "External_Identifier"),
    "Spring 2026 PWGVAL dilution workbook"
  )

  group_id <- as.character(workbook$Group_ID)
  tf_dilution <- suppressWarnings(as.numeric(workbook$TF_Dilution))
  # The final token of Group_ID is the technical replicate identifier. The
  # leading VA-## token identifies the patient whose baseline mutation list must
  # match before the row is scored.
  replicate <- sub("^.*-", "", group_id)
  patient <- sub("^([A-Za-z]+-[0-9]+).*$", "\\1", group_id)
  m4chip_id <- sub("^(M4CHIP_[0-9]+).*$", "\\1", as.character(workbook$Alias))
  bam <- paste0(
    # Construct the pipeline-suite BAM basename from workbook fields rather than
    # relying on free-text filenames. This keeps the metadata reproducible from
    # the sequencing submission sheet.
    m4chip_id,
    "_",
    as.character(workbook$Tissue_Origin),
    "_",
    as.character(workbook$Tissue_Type),
    "_",
    as.character(workbook$Library_Design),
    "_",
    group_id,
    ".filter.deduped.recalibrated.bam"
  )

  workbook %>%
    dplyr::mutate(
      M4CHIP_ID = m4chip_id,
      Merge = paste0(group_id, "-P"),
      Patient = patient,
      replicate = replicate,
      Bam = bam,
      BAM = bam,
      LOD = tf_dilution,
      Sample_type = "Blood_plasma_cfDNA",
      Sample_type_Bam = "Blood_plasma_cfDNA",
      Timepoint = as.character(tf_dilution),
      timepoint_info = paste0("Dilution: ", as.character(tf_dilution)),
      Sample_ID = paste0(patient, "_PWGVAL_LOD_", as.character(tf_dilution), "_rep", replicate),
      dilution_series = "PWGVAL_M4CHIP",
      metadata_source = path
    ) %>%
    dplyr::rename(
      `TF_Dilution_pct` = "TF_Dilution",
      `Tissue Origin` = "Tissue_Origin",
      `Tissue Type` = "Tissue_Type",
      `Library Design` = "Library_Design",
      `Index(es)` = "Index_es",
      `i7 Index` = "i7_Index",
      `i5 Index` = "i5_Index",
      `External_Identifier` = "External_Identifier"
    )
}

load_spring2026_pwgval_dilution_metadata <- function(required = FALSE) {
  # Load cached PWGVAL dilution metadata when present; otherwise rebuild from the
  # workbook. This lets the pipeline run from a clean checkout if the workbook is
  # present, while still using the reviewed CSV when it has been generated.
  path <- spring2026_pwgval_dilution_metadata_path()
  if (!spring2026_revision_enabled()) return(NULL)
  if (!file.exists(path)) {
    metadata <- build_spring2026_pwgval_dilution_metadata_from_workbook()
    if (is.null(metadata)) {
      if (isTRUE(required)) {
        stop(
          "Missing Spring 2026 PWGVAL dilution metadata and workbook fallback: ",
          path,
          " / ",
          spring2026_pwgval_dilution_workbook_path(),
          call. = FALSE
        )
      }
      return(NULL)
    }
  } else {
    metadata <- readr::read_csv(path, show_col_types = FALSE)
  }
  require_columns(
    metadata,
    c("BAM", "Bam", "LOD", "Patient", "Sample_type", "Timepoint", "timepoint_info", "Sample_ID", "Merge"),
    "Spring 2026 PWGVAL dilution metadata"
  )
  metadata
}

build_spring2026_pwgval_dilution_tumor_fraction <- function() {
  # Join PWGVAL dilution metadata to ichorCNA parameter summaries. This table is
  # separate from the workbook dilution percentage: Tumor_fraction is the
  # observed ichorCNA estimate, while LOD/TF_Dilution_pct is the intended spike-in
  # dilution level.
  metadata <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
  params_path <- file.path(
    spring2026_revision_data_dir(),
    "Ichor_CNA_combined_params_summary_tumor_fraction_sex.tsv"
  )
  if (is.null(metadata) || !file.exists(params_path)) return(NULL)
  params <- readr::read_tsv(params_path, show_col_types = FALSE) %>%
    dplyr::filter(grepl("^M4CHIP_", .data$file)) %>%
    dplyr::mutate(Bam = sub("[.]params[.]txt$", "", .data$file))
  metadata %>%
    dplyr::mutate(Bam = sub("[.]bam$", "", .data$Bam)) %>%
    dplyr::inner_join(params, by = "Bam") %>%
    dplyr::transmute(
      Bam,
      Tumor_fraction = .data$tumor_fraction,
      Ploidy = .data$ploidy,
      Metadata_TF_Dilution_pct = .data$LOD,
      Patient = .data$Patient,
      Merge = .data$Merge,
      Group_ID = .data$Group_ID
    )
}

read_dilution_metadata_with_spring2026 <- function(path) {
  # Append Spring 2026 PWGVAL dilution rows to the historical dilution metadata.
  # Patient_Bam and Sample_type_Bam are filled from Patient/Sample_type when the
  # historical schema lacks those explicit BAM-level fields.
  if (!file.exists(path)) stop("Missing dilution metadata: ", path, call. = FALSE)
  current <- readr::read_csv(path, show_col_types = FALSE)
  revision <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
  if (is.null(revision)) return(current)

  shared_cols <- union(names(current), names(revision))
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  current <- current[, shared_cols, drop = FALSE]
  revision <- revision[, shared_cols, drop = FALSE]

  combined <- dplyr::bind_rows(current, revision)
  if (!"Patient_Bam" %in% names(combined)) combined$Patient_Bam <- NA_character_
  if (!"Sample_type_Bam" %in% names(combined)) combined$Sample_type_Bam <- NA_character_
  if (!"Patient" %in% names(combined)) combined$Patient <- NA_character_
  if (!"Sample_type" %in% names(combined)) combined$Sample_type <- NA_character_

  combined %>%
    dplyr::mutate(
      Patient_Bam = dplyr::coalesce(.data$Patient_Bam, .data$Patient),
      Sample_type_Bam = dplyr::coalesce(.data$Sample_type_Bam, .data$Sample_type)
    ) %>%
    dplyr::distinct(BAM, Bam, Merge, .keep_all = TRUE)
}

read_dilution_tumor_fraction_with_spring2026 <- function(path) {
  # Append Spring 2026 tumor-fraction rows to the historical dilution tumor-
  # fraction table. Distinct Bam keeps the first row per BAM after combining to
  # avoid duplicated LOD rows in downstream model application.
  if (!file.exists(path)) stop("Missing dilution tumor-fraction table: ", path, call. = FALSE)
  current <- readr::read_tsv(path, show_col_types = FALSE)
  revision_path <- spring2026_pwgval_dilution_tumor_fraction_path()
  if (!spring2026_revision_enabled()) return(current)
  revision <- if (file.exists(revision_path)) {
    readr::read_tsv(revision_path, show_col_types = FALSE)
  } else {
    build_spring2026_pwgval_dilution_tumor_fraction()
  }
  if (is.null(revision)) return(current)

  shared_cols <- union(names(current), names(revision))
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  current <- current[, shared_cols, drop = FALSE]
  revision <- revision[, shared_cols, drop = FALSE]

  dplyr::bind_rows(current, revision) %>%
    dplyr::distinct(Bam, .keep_all = TRUE)
}
