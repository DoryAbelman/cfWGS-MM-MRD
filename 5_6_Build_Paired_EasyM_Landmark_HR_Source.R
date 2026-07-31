#!/usr/bin/env Rscript

# Build the manuscript-facing paired cfWGS/EasyM landmark survival summary.
#
# Inputs are the current, dated patient-level source tables used for Extended
# Data Fig. 6G-H. Binary EasyM analyses use the isotype-specific Rapid Novor
# reference-threshold call already encoded in those tables. Continuous effects
# are reported per one standard-deviation increase.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(tibble)
})

project_root <- normalizePath(".", mustWork = TRUE)
source_dir <- file.path(project_root, "Output_tables_2025", "Source_data")
output_path <- file.path(
  source_dir,
  "Paired_EasyM_landmark_HR_source_data_2026-07-31.csv"
)

landmarks <- tribble(
  ~landmark, ~source_file, ~expected_n, ~expected_events,
  "One-year maintenance",
  "All_train_test_KM_source_data_1yr_maintenance_EasyM_reference_threshold_binary_2026-07-30.csv",
  18L, 8L,
  "Post-ASCT",
  "All_train_test_KM_source_data_post_transplant_EasyM_reference_threshold_binary_2026-07-30.csv",
  16L, 7L
)

predictors <- tribble(
  ~predictor, ~column, ~scale_continuous, ~definition,
  "BM-informed cfWGS call",
  "BM_zscore_only_detection_rate_call",
  FALSE,
  "Locked BM-informed cVAF-model call",
  "EasyM reference-threshold call",
  "EasyM_reference_threshold_binary",
  FALSE,
  paste(
    "Isotype-specific Rapid Novor reference threshold:",
    "negative if IgG <=1% or IgA/light-chain <=0.05% of baseline"
  ),
  "BM-informed cfWGS probability per 1 SD",
  "BM_zscore_only_detection_rate_prob",
  TRUE,
  "Locked BM-informed cVAF-model probability, standardized within landmark",
  "EasyM residual value per 1 SD",
  "EasyM_value",
  TRUE,
  "Residual monoclonal-protein percentage, standardized within landmark"
)

required_base <- c("Time_to_event", "Relapsed_Binary")

fit_one <- function(data, predictor_row, landmark, source_file) {
  column <- predictor_row$column[[1]]
  required <- c(required_base, column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", source_file, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  analysis <- data %>%
    select(all_of(required)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  if (nrow(analysis) < 5L || sum(analysis$Relapsed_Binary) < 2L) {
    stop(sprintf("Insufficient evaluable data for %s at %s", column, landmark), call. = FALSE)
  }
  if (dplyr::n_distinct(analysis[[column]]) < 2L) {
    stop(sprintf("No predictor variation for %s at %s", column, landmark), call. = FALSE)
  }

  model_column <- column
  if (isTRUE(predictor_row$scale_continuous[[1]])) {
    model_column <- paste0(column, "_z")
    analysis[[model_column]] <- as.numeric(scale(analysis[[column]]))
  }

  fit <- coxph(
    as.formula(sprintf("Surv(Time_to_event, Relapsed_Binary) ~ %s", model_column)),
    data = analysis,
    ties = "efron"
  )
  result <- summary(fit)

  tibble(
    Landmark = landmark,
    Predictor = predictor_row$predictor[[1]],
    Source_column = column,
    N = nrow(analysis),
    Events = sum(analysis$Relapsed_Binary),
    HR = unname(exp(coef(fit))[[1]]),
    CI_low = unname(result$conf.int[1, "lower .95"]),
    CI_high = unname(result$conf.int[1, "upper .95"]),
    P_value = unname(result$coefficients[1, "Pr(>|z|)"]),
    Definition = predictor_row$definition[[1]],
    Source_file = source_file
  )
}

results <- list()
for (landmark_index in seq_len(nrow(landmarks))) {
  landmark_row <- landmarks[landmark_index, ]
  source_path <- file.path(source_dir, landmark_row$source_file[[1]])
  if (!file.exists(source_path)) {
    stop(sprintf("Missing current EasyM landmark source: %s", source_path), call. = FALSE)
  }

  data <- read_csv(source_path, show_col_types = FALSE)
  if (
    nrow(data) != landmark_row$expected_n[[1]] ||
      sum(data$Relapsed_Binary, na.rm = TRUE) != landmark_row$expected_events[[1]]
  ) {
    stop(
      sprintf(
        "%s denominator drift: observed %d rows/%d events; expected %d/%d",
        landmark_row$landmark[[1]],
        nrow(data),
        sum(data$Relapsed_Binary, na.rm = TRUE),
        landmark_row$expected_n[[1]],
        landmark_row$expected_events[[1]]
      ),
      call. = FALSE
    )
  }

  for (predictor_index in seq_len(nrow(predictors))) {
    results[[length(results) + 1L]] <- fit_one(
      data,
      predictors[predictor_index, ],
      landmark_row$landmark[[1]],
      landmark_row$source_file[[1]]
    )
  }
}

result_table <- bind_rows(results)

expected_rounded <- tribble(
  ~Landmark, ~Predictor, ~HR_expected,
  "One-year maintenance", "BM-informed cfWGS call", 6.08,
  "One-year maintenance", "EasyM reference-threshold call", 4.63,
  "Post-ASCT", "BM-informed cfWGS call", 3.51,
  "Post-ASCT", "EasyM reference-threshold call", 1.29
)
verification <- result_table %>%
  inner_join(expected_rounded, by = c("Landmark", "Predictor")) %>%
  mutate(matches = round(HR, 2) == HR_expected)
if (nrow(verification) != nrow(expected_rounded) || any(!verification$matches)) {
  stop("Paired EasyM landmark hazard-ratio verification failed.", call. = FALSE)
}

write_csv(result_table, output_path, na = "NA")
message("WROTE ", output_path)
message("Rows: ", nrow(result_table))
