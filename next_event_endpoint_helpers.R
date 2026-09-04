################################################################################
## Shared next-event endpoint helpers
##
## Purpose:
##   For sample-level relapse/progression plots and survival summaries, select
##   the first known progression event after the sample date. If the first PFS
##   event happened before the sample, it is ignored and a later progression is
##   used when available; otherwise the row is censored at the latest known
##   follow-up/sample date after the sample.
##
## Workflow role:
##   These functions are sourced by 2_4, 2_4B, 4_1, and 4_1B. They contribute
##   the sample-relative progression/censor endpoints used in longitudinal and
##   survival figures but do not generate a figure or table directly.
##
## Main inputs:
##   * An RDS PFS table with Patient, Baseline_Date, Censor_date, and Relapsed.
##   * An RDS table containing all known Patient/Progression_date pairs.
##   * Optional patient follow-up tables, latest-date tables, and sample dates
##     that can provide later censoring dates.
##   * A sample-level data frame with patient and collection-date columns.
##
## Returned values:
##   load_next_event_endpoint_resources() returns the cleaned PFS table, all
##   progression dates, and censor-date candidates. add_next_event_endpoint()
##   appends the selected event/censor date, status, time from sample, source,
##   prior-event flags, and manual-override flag to the input rows.
##
## Key assumptions:
##   A progression up to event_grace_days before a sample is assigned to that
##   sample and its time is set to zero. Earlier progressions are ignored when a
##   later event or censor date is available. The two default day-zero overrides
##   are described at the add_next_event_endpoint() argument below.
##
## R packages used by these functions:
##   dplyr, lubridate, purrr, readr, tibble, and tidyr.
################################################################################

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))

  out <- suppressWarnings(as.Date(x))
  missing <- is.na(out) & !is.na(x) & nzchar(as.character(x))
  out[missing] <- suppressWarnings(lubridate::ymd(as.character(x[missing])))
  missing <- is.na(out) & !is.na(x) & nzchar(as.character(x))

  # 2026-07-29 audit fix: warn on order-dependent ambiguity.
  #
  # The mdy -> dmy fallback chain silently resolves ambiguous strings such as
  # "03/04/2021" as month/day, because mdy is tried first. Every clinical date
  # in the survival analyses flows through this function, so a silent
  # misinterpretation would shift an event or censor date by up to eleven months
  # without any diagnostic. Ambiguous values are now reported rather than
  # resolved quietly. Unambiguous slash dates (day > 12) are unaffected.
  if (any(missing)) {
    candidates <- as.character(x[missing])
    mdy_parsed <- suppressWarnings(lubridate::mdy(candidates))
    dmy_parsed <- suppressWarnings(lubridate::dmy(candidates))
    ambiguous <- !is.na(mdy_parsed) & !is.na(dmy_parsed) & mdy_parsed != dmy_parsed
    if (any(ambiguous)) {
      warning(
        "Ambiguous date string(s) resolved as month/day: ",
        paste(unique(candidates[ambiguous]), collapse = ", "),
        ". Supply an unambiguous format (ISO 8601) to remove this warning.",
        call. = FALSE
      )
    }
    out[missing] <- mdy_parsed
  }

  missing <- is.na(out) & !is.na(x) & nzchar(as.character(x))
  out[missing] <- suppressWarnings(lubridate::dmy(as.character(x[missing])))
  out
}

read_optional_followup_table <- function(rds_path = NULL, csv_path = NULL,
                                         fallback_col = "followup_end_date",
                                         source_label = "follow-up") {
  if (!is.null(rds_path) && file.exists(rds_path)) {
    tbl <- readRDS(rds_path)
  } else if (!is.null(csv_path) && file.exists(csv_path)) {
    tbl <- readr::read_csv(csv_path, show_col_types = FALSE)
  } else {
    return(tibble::tibble(Patient = character(), censor_candidate_date = as.Date(character()),
                          censor_candidate_source = character()))
  }

  if (!fallback_col %in% names(tbl)) {
    date_cols <- intersect(c("followup_end_date", "latest_date", "censor_date"), names(tbl))
    if (length(date_cols) == 0) {
      stop("No date column found in optional follow-up table: ", source_label, call. = FALSE)
    }
    fallback_col <- date_cols[[1]]
  }

  tbl %>%
    dplyr::transmute(
      Patient = as.character(.data$Patient),
      censor_candidate_date = as_date_safe(.data[[fallback_col]]),
      censor_candidate_source = source_label
    ) %>%
    dplyr::filter(!is.na(Patient), !is.na(censor_candidate_date))
}

load_next_event_endpoint_resources <- function(pfs_path,
                                               relapse_dates_path,
                                               followup_rds_path = NULL,
                                               followup_csv_path = NULL,
                                               latest_dates_paths = character(),
                                               sample_data = NULL,
                                               patient_col = "Patient",
                                               sample_date_col = "Date") {
  if (!file.exists(pfs_path)) {
    stop("Missing PFS endpoint table: ", pfs_path, call. = FALSE)
  }
  if (!file.exists(relapse_dates_path)) {
    stop("Missing full relapse/progression dates table: ", relapse_dates_path, call. = FALSE)
  }

  pfs <- readRDS(pfs_path) %>%
    dplyr::rename_with(tolower, dplyr::any_of(c("Baseline_Date", "Censor_date", "Relapsed"))) %>%
    dplyr::transmute(
      Patient = as.character(.data$Patient),
      baseline_date = as_date_safe(.data$baseline_date),
      first_pfs_date = as_date_safe(.data$censor_date),
      first_pfs_event = as.integer(.data$relapsed)
    )

  progression_dates <- readRDS(relapse_dates_path) %>%
    dplyr::transmute(
      Patient = as.character(.data$Patient),
      progression_date = as_date_safe(.data$Progression_date)
    ) %>%
    dplyr::filter(!is.na(Patient), !is.na(progression_date)) %>%
    dplyr::distinct()

  pfs_event_dates <- pfs %>%
    dplyr::filter(first_pfs_event == 1L, !is.na(first_pfs_date)) %>%
    dplyr::transmute(Patient, progression_date = first_pfs_date)

  progression_dates <- dplyr::bind_rows(progression_dates, pfs_event_dates) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Patient, progression_date)

  censor_candidates <- pfs %>%
    dplyr::transmute(
      Patient,
      censor_candidate_date = first_pfs_date,
      censor_candidate_source = dplyr::if_else(
        first_pfs_event == 1L,
        "first PFS progression date",
        "PFS censor date"
      )
    ) %>%
    dplyr::filter(!is.na(censor_candidate_date))

  followup_tbl <- read_optional_followup_table(
    rds_path = followup_rds_path,
    csv_path = followup_csv_path,
    fallback_col = "followup_end_date",
    source_label = "patient follow-up date"
  )

  latest_tbls <- purrr::map(latest_dates_paths, function(path) {
    if (!file.exists(path)) return(NULL)
    read_optional_followup_table(
      csv_path = path,
      fallback_col = "latest_date",
      source_label = paste0("latest date: ", basename(path))
    )
  }) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  sample_censor_tbl <- tibble::tibble()
  if (!is.null(sample_data)) {
    if (!all(c(patient_col, sample_date_col) %in% names(sample_data))) {
      stop("sample_data is missing patient/date columns for endpoint censor candidates.", call. = FALSE)
    }
    sample_censor_tbl <- sample_data %>%
      dplyr::transmute(
        Patient = as.character(.data[[patient_col]]),
        censor_candidate_date = as_date_safe(.data[[sample_date_col]]),
        censor_candidate_source = "latest sample date in analysis input"
      ) %>%
      dplyr::filter(!is.na(Patient), !is.na(censor_candidate_date)) %>%
      dplyr::group_by(Patient) %>%
      dplyr::summarise(
        censor_candidate_date = max(censor_candidate_date),
        censor_candidate_source = "latest sample date in analysis input",
        .groups = "drop"
      )
  }

  censor_candidates <- dplyr::bind_rows(
    censor_candidates,
    followup_tbl,
    latest_tbls,
    sample_censor_tbl
  ) %>%
    dplyr::filter(!is.na(Patient), !is.na(censor_candidate_date)) %>%
    dplyr::distinct()

  list(
    pfs = pfs,
    progression_dates = progression_dates,
    censor_candidates = censor_candidates
  )
}

add_next_event_endpoint <- function(data,
                                    endpoint_resources,
                                    sample_date_col = "sample_date",
                                    event_grace_days = 30L,
                                    patient_col = "Patient",
                                    # ------------------------------------------------------------------ #
                                    # Manual per-patient endpoint overrides.
                                    #
                                    # 2026-07-29 audit note. These two entries force the named
                                    # sample to `endpoint_status = 1L` at day 0, i.e. they assert
                                    # that the listed draw was taken at clinical progression.
                                    #
                                    # These are explicit manual endpoint assignments rather than
                                    # dates inferred by the general next-event rule. Pass
                                    # `relapse_sample_day0_overrides = NULL` to run the existing
                                    # sensitivity analysis without either assignment.
                                    #
                                    # Deliberately exposed as an argument (rather than hard-coded in
                                    # the body) so that the with/without comparison can be run
                                    # without editing this file.
                                    # ------------------------------------------------------------------ #
                                    relapse_sample_day0_overrides = tibble::tribble(
                                      ~Patient, ~relapse_sample_date,
                                      "CA-10", as.Date("2023-07-11"),
                                      "VA-07", as.Date("2021-01-01")
                                    )) {
  if (is.null(relapse_sample_day0_overrides)) {
    relapse_sample_day0_overrides <- tibble::tibble(
      Patient = character(), relapse_sample_date = as.Date(character())
    )
  }
  if (!all(c(patient_col, sample_date_col) %in% names(data))) {
    stop("Input data missing patient/sample-date columns for next-event endpoint.", call. = FALSE)
  }

  endpoint_tbl <- data %>%
    dplyr::mutate(.endpoint_row_id = dplyr::row_number()) %>%
    dplyr::transmute(
      .endpoint_row_id,
      Patient = as.character(.data[[patient_col]]),
      sample_date = as_date_safe(.data[[sample_date_col]])
    )

  event_candidates <- endpoint_tbl %>%
    dplyr::left_join(endpoint_resources$progression_dates, by = "Patient", relationship = "many-to-many") %>%
    dplyr::filter(
      !is.na(sample_date),
      !is.na(progression_date),
      progression_date >= sample_date - lubridate::days(event_grace_days)
    ) %>%
    dplyr::mutate(days_to_event_raw = as.numeric(progression_date - sample_date)) %>%
    dplyr::arrange(.endpoint_row_id, progression_date) %>%
    dplyr::group_by(.endpoint_row_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      .endpoint_row_id,
      next_progression_date = progression_date,
      next_progression_days_raw = days_to_event_raw,
      next_progression_days = pmax(days_to_event_raw, 0),
      endpoint_date = progression_date,
      endpoint_status = 1L,
      endpoint_type = "next_progression",
      endpoint_source = dplyr::if_else(
        days_to_event_raw < 0,
        paste0("progression within ", event_grace_days, "d pre-sample grace"),
        "progression date after sample"
      )
    )

  censor_candidates <- endpoint_tbl %>%
    dplyr::left_join(endpoint_resources$censor_candidates, by = "Patient", relationship = "many-to-many") %>%
    dplyr::filter(
      !is.na(sample_date),
      !is.na(censor_candidate_date),
      censor_candidate_date >= sample_date
    ) %>%
    dplyr::arrange(.endpoint_row_id, dplyr::desc(censor_candidate_date)) %>%
    dplyr::group_by(.endpoint_row_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      .endpoint_row_id,
      censor_endpoint_date = censor_candidate_date,
      censor_endpoint_days = as.numeric(censor_candidate_date - sample_date),
      censor_endpoint_source = censor_candidate_source
    )

  first_pfs_by_patient <- endpoint_resources$pfs %>%
    dplyr::select(Patient, baseline_date, first_pfs_date, first_pfs_event)

  ignored_prior_events <- endpoint_tbl %>%
    dplyr::left_join(endpoint_resources$progression_dates, by = "Patient", relationship = "many-to-many") %>%
    dplyr::filter(
      !is.na(sample_date),
      !is.na(progression_date),
      progression_date < sample_date - lubridate::days(event_grace_days)
    ) %>%
    dplyr::group_by(.endpoint_row_id) %>%
    dplyr::summarise(
      n_prior_progressions_before_sample = dplyr::n(),
      latest_prior_progression_date = max(progression_date),
      .groups = "drop"
    )

  row_endpoint <- endpoint_tbl %>%
    dplyr::left_join(first_pfs_by_patient, by = "Patient") %>%
    dplyr::left_join(event_candidates, by = ".endpoint_row_id") %>%
    dplyr::left_join(censor_candidates, by = ".endpoint_row_id") %>%
    dplyr::left_join(ignored_prior_events, by = ".endpoint_row_id") %>%
    dplyr::mutate(
      n_prior_progressions_before_sample = tidyr::replace_na(n_prior_progressions_before_sample, 0L),
      endpoint_date = dplyr::coalesce(endpoint_date, censor_endpoint_date),
      endpoint_status = tidyr::replace_na(endpoint_status, 0L),
      endpoint_type = dplyr::case_when(
        !is.na(next_progression_date) ~ endpoint_type,
        !is.na(censor_endpoint_date) ~ "censor",
        TRUE ~ "missing_endpoint"
      ),
      endpoint_source = dplyr::coalesce(endpoint_source, censor_endpoint_source),
      endpoint_days_from_sample = as.numeric(endpoint_date - sample_date),
      endpoint_days_from_sample = dplyr::if_else(
        endpoint_status == 1L & endpoint_days_from_sample < 0 &
          endpoint_days_from_sample >= -event_grace_days,
        0,
        endpoint_days_from_sample
      ),
      first_pfs_days_from_sample = as.numeric(first_pfs_date - sample_date),
      endpoint_uses_later_progression_after_first_pfs =
        endpoint_status == 1L &
        !is.na(first_pfs_date) &
        !is.na(next_progression_date) &
        next_progression_date > first_pfs_date,
      endpoint_ignores_prior_progression =
        n_prior_progressions_before_sample > 0
    ) %>%
    dplyr::mutate(
      force_relapse_sample_day0 = paste(Patient, sample_date) %in%
        paste(relapse_sample_day0_overrides$Patient,
              relapse_sample_day0_overrides$relapse_sample_date),
      endpoint_date = dplyr::if_else(
        force_relapse_sample_day0,
        sample_date,
        endpoint_date
      ),
      endpoint_status = dplyr::if_else(
        force_relapse_sample_day0,
        1L,
        endpoint_status
      ),
      endpoint_type = dplyr::if_else(
        force_relapse_sample_day0,
        "relapse_sample_day0",
        endpoint_type
      ),
      endpoint_source = dplyr::if_else(
        force_relapse_sample_day0,
        "relapse-labeled sample assigned event at day 0",
        endpoint_source
      ),
      endpoint_days_from_sample = as.numeric(endpoint_date - sample_date),
      endpoint_days_from_sample = dplyr::if_else(
        endpoint_status == 1L & endpoint_days_from_sample < 0 &
          endpoint_days_from_sample >= -event_grace_days,
        0,
        endpoint_days_from_sample
      ),
      endpoint_uses_later_progression_after_first_pfs = dplyr::if_else(
        force_relapse_sample_day0,
        FALSE,
        endpoint_uses_later_progression_after_first_pfs
      ),
      endpoint_ignores_prior_progression = dplyr::if_else(
        force_relapse_sample_day0,
        FALSE,
        endpoint_ignores_prior_progression
      )
    ) %>%
    dplyr::select(
      .endpoint_row_id,
      baseline_date,
      first_pfs_date,
      first_pfs_event,
      first_pfs_days_from_sample,
      next_progression_date,
      next_progression_days_raw,
      next_progression_days,
      endpoint_date,
      endpoint_status,
      endpoint_type,
      endpoint_source,
      endpoint_days_from_sample,
      force_relapse_sample_day0,
      n_prior_progressions_before_sample,
      latest_prior_progression_date,
      endpoint_uses_later_progression_after_first_pfs,
      endpoint_ignores_prior_progression
    )

  data %>%
    dplyr::mutate(.endpoint_row_id = dplyr::row_number()) %>%
    dplyr::left_join(row_endpoint, by = ".endpoint_row_id") %>%
    dplyr::select(-.endpoint_row_id)
}
