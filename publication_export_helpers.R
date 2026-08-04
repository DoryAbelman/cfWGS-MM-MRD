# Publication-facing export helpers
#
# Internal analysis objects retain their historical cohort codes because those
# codes are used in filters and frozen-model validation.  Final tables and
# source-data exports use the manuscript terminology requested for submission:
# Training and Testing.

publication_cohort_label <- function(x) {
  original_na <- is.na(x)
  value <- trimws(as.character(x))
  normalized <- tolower(gsub("[_ ]+", " ", value))
  normalized <- gsub("[[:space:]]+", " ", normalized)

  value[normalized %in% c(
    "frontline",
    "frontline cohort",
    "frontline induction-transplant",
    "training cohort"
  )] <- "Training"
  value[normalized %in% c(
    "non-frontline",
    "non frontline",
    "nonfrontline",
    "non-frontline cohort",
    "non frontline cohort",
    "nonfrontline cohort",
    "test cohort",
    "testing cohort"
  )] <- "Testing"
  value[normalized == "frontline omitted"] <- "Training omitted"
  value[original_na] <- NA_character_
  value
}

assert_publication_cohort_labels <- function(data, context = "publication export") {
  if (!is.data.frame(data)) {
    stop(context, " must be a data frame.", call. = FALSE)
  }

  stale_names <- names(data)[grepl("frontline", names(data), ignore.case = TRUE)]
  stale_values <- unique(unlist(lapply(data, function(column) {
    if (!is.character(column) && !is.factor(column)) return(character())
    values <- as.character(column)
    values[!is.na(values) & grepl("frontline", values, ignore.case = TRUE)]
  }), use.names = FALSE))

  if (length(stale_names) || length(stale_values)) {
    stop(
      context,
      " still contains internal cohort terminology after relabeling. Columns: ",
      if (length(stale_names)) paste(stale_names, collapse = ", ") else "none",
      "; values: ",
      if (length(stale_values)) paste(utils::head(stale_values, 10L), collapse = " | ") else "none",
      call. = FALSE
    )
  }
  invisible(data)
}

relabel_publication_cohorts <- function(data, context = "publication export") {
  if (!is.data.frame(data)) {
    stop(context, " must be a data frame.", call. = FALSE)
  }

  output <- data
  cohort_columns <- names(output)[vapply(
    output,
    function(column) {
      if (!is.character(column) && !is.factor(column)) return(FALSE)
      values <- as.character(column)
      any(!is.na(values) & grepl("frontline", values, ignore.case = TRUE))
    },
    logical(1)
  )]
  output[cohort_columns] <- lapply(output[cohort_columns], publication_cohort_label)

  names(output) <- gsub(
    "non[-_ ]?frontline",
    "testing",
    names(output),
    ignore.case = TRUE,
    perl = TRUE
  )
  names(output) <- gsub(
    "frontline",
    "training",
    names(output),
    ignore.case = TRUE,
    perl = TRUE
  )

  assert_publication_cohort_labels(output, context)
  output
}

relabel_publication_workbook_tables <- function(tables, context = "publication workbook") {
  if (!is.list(tables) || is.null(names(tables)) || any(!nzchar(names(tables)))) {
    stop(context, " must be a named list of data frames.", call. = FALSE)
  }
  output <- Map(
    function(data, sheet) relabel_publication_cohorts(
      data,
      paste0(context, " sheet ", sheet)
    ),
    tables,
    names(tables)
  )
  names(output) <- names(tables)
  output
}
