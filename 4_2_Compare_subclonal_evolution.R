###############################################################################
# 4_2_Compare_subclonal_evolution.R
#
# Purpose:
#   Copy-number-based assessment of emergent subclonal CNA events between
#   paired baseline and relapse cfDNA WGS samples (30-40x). Identifies CNAs
#   present at relapse but absent (or sub-threshold) at baseline, generating
#   per-patient longitudinal CNA plots, an emergent-event summary table, and
#   plain-language sentences for the manuscript Results section.
#
# Inputs:
#   - Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#       (integrated feature table; output of 2_0)
#   - cohort_assignment_table_updated.rds
#       (M4/SPORE/IMMAGINE cohort labels for patient filtering)
#
# Outputs:
#   - Subclonal_evolution_plots.pdf  (per-patient CNA track plots)
#   - Emergent_CNA_events.csv        (one row per emergent CNA event)
#   - Printed text: auto-generated manuscript sentences (stdout)
#
# Dependencies:
#   tidyverse, lubridate
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/4_2_Compare_subclonal_evolution.R
#
# Manuscript outputs created/updated:
#   - Extended Data Figure 10A-B: emergent CNA/subclonal evolution support
#     tables and frozen final ichorCNA figure components when available.
#
# Pipeline role:
#   The current final figure panels use genome-wide ichorCNA visualizations
#   assembled outside this R script. This script regenerates the in-repository
#   event tables that identify which CNAs are newly detected at relapse, and
#   places those provenance files beside the final figure artifacts.
#
# Author:    Dory Abelman
# Last update: August 2025
###############################################################################
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

## ---- 0. USER SETTINGS -------------------------------------------------------
in_rds        <- "Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds"
out_plot_pdf  <- "Subclonal_evolution_plots.pdf"
out_events_csv<- "Emergent_CNA_events.csv"
cohort_df <- readRDS("cohort_assignment_table_updated.rds")
outdir <- "Final Tables and Figures/"

## ---- 1. Final Tables and Figures/4E_performance_nested_folds_bm_validation_updated.png## ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# Shared manuscript-output helpers.
# Extended Data Figure 10A/B uses external ichorCNA component plots assembled
# outside this R script. This script contributes the in-repo support/provenance
# CSVs and preserves the frozen final PDF when it is present locally.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

## ---- 2. Load data -----------------------------------------------------------
all_df <- if (grepl("\\.rds$", in_rds, ignore.case = TRUE)) {
  readRDS(in_rds)
} else {
  read_csv(in_rds)          # assumes header row
}

## ---- 3. Keep cfDNA & harmonise flags ----------------------------------------
cfDNA_df <- all_df %>%
  filter(str_detect(Sample_type, regex("plasma", ignore_case = TRUE))) %>%
  filter(Patient %in% cohort_df$Patient) %>%
  mutate(
    # Baseline ≈ diagnosis / T0
    is_baseline = timepoint_info %in% c("Diagnosis", "Baseline"),
    # Relapse / progression
    is_relapse  = timepoint_info %in% c("Relapse", "Progression"),
    Sample_Date = as_date(Date_of_sample_collection)
  ) 

## ---- 4. Keep patients with baseline *and* relapse cfDNA ---------------------
pts_keep <- cfDNA_df %>%
  group_by(Patient) %>%
  summarise(has_bl = any(is_baseline),
            has_rl = any(is_relapse),
            .groups = "drop") %>%
  filter(has_bl & has_rl) %>%
  pull(Patient)

cfDNA_df <- cfDNA_df %>% filter(Patient %in% pts_keep)

## ---- 5. Simple FGA proxy + numeric CNA flags --------------------------------
event_cols <- c("del1p", "amp1q", "del13q", "del17p")

cfDNA_df <- cfDNA_df %>%
  mutate(across(all_of(event_cols), ~ as.numeric(.x))) %>%   # TRUE → 1, FALSE → 0
  rowwise() %>%
  mutate(FGA_proxy = sum(c_across(all_of(event_cols)), na.rm = TRUE) /
           length(event_cols)) %>%
  ungroup()

## ---- 6. Emergent CNA table --------------------------------------------------
# DEFINITION: An "emergent" CNA is one that is ABSENT (0) at baseline and
# PRESENT (1) at relapse.  The delta is computed per-patient as
#   relapse_value - baseline_value.
# A delta of +1 therefore captures 0->1 (newly acquired), whereas:
#   delta of  0 = no change (present at both or absent at both)
#   delta of -1 = lost at relapse (present at baseline, absent at relapse)
# Only delta == 1 events are retained as "Emergent" and written to CSV.
# NOTE: event_cols are binary (0/1) flags produced by 1_4_Process_CNA_Data.R.
emergent_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - .x[which(is_baseline)][1], .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(Event   = str_remove(Event, "^delta_"),
         Emergent = Delta == 1) %>%
  filter(Emergent) %>%
  select(Patient, Sample, Event) %>%
  arrange(Patient, Event)

edfig10_total_events_csv <- file.path(outdir, "Emergent_CNA_event_total.csv")
write_csv(emergent_tbl, edfig10_total_events_csv)


# Gains:
gains_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - first(.x), .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(
    Event    = str_remove(Event, "^delta_"),
    Emergent = Delta == 1
  ) %>%
  filter(Emergent) %>%
  select(Patient, Sample, Event)

edfig10_gain_events_csv <- file.path(outdir, "Emergent_CNA_event_gains.csv")
write_csv(gains_tbl, edfig10_gain_events_csv)

# Losses:
losses_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - first(.x), .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(
    Event = str_remove(Event, "^delta_"),
    Lost  = Delta == -1
  ) %>%
  filter(Lost) %>%
  select(Patient, Sample, Event)

edfig10_loss_events_csv <- file.path(outdir, "Emergent_CNA_loss_events.csv")
write_csv(losses_tbl, edfig10_loss_events_csv)

# MANUSCRIPT OUTPUT: Extended Data Figure 10A/B support data
# The current final figure panels are externally assembled from ichorCNA
# genome-wide CNA plots. The CSVs generated here are the in-repo provenance
# tables that identify emergent/gained/lost CNA events used to support those
# panels. They are copied beside the final figure artifacts so reviewers can
# trace the figure back to the tabular calls regenerated by this script.
ms_copy_artifact(
  source_path = edfig10_total_events_csv,
  artifact_id = "EDFIG10A_B_A",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10A support data: emergent CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

ms_copy_artifact(
  source_path = edfig10_gain_events_csv,
  artifact_id = "EDFIG10A_B_B",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10B support data: gained CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

ms_copy_artifact(
  source_path = edfig10_loss_events_csv,
  artifact_id = "EDFIG10A_B_B",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10B support data: lost CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

edfig10_final_pdf_candidates <- c(
  file.path("Manuscript_Exports", "02_extended_data_figures", "Extended_Data_Figure_10", "final_artifacts", "Extended_Data_Figure_10.pdf"),
  file.path("Figures_Exported", "Final_Feb2026", "Extended_Data_Figure_10.pdf"),
  file.path("reproducible_workflow", "outputs", "frozen", "extended_figures", "Extended_Data_Figure_10.pdf")
)
edfig10_final_pdf <- edfig10_final_pdf_candidates[file.exists(edfig10_final_pdf_candidates)][1]

if (!is.na(edfig10_final_pdf)) {
  ms_copy_artifact(
    source_path = edfig10_final_pdf,
    artifact_id = "EDFIG10A_B_A",
    role = "external_final_figure_pdf",
    description = "Extended Data Figure 10 final PDF: externally assembled ichorCNA component figure preserved as a frozen manuscript artifact.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )
  ms_copy_artifact(
    source_path = edfig10_final_pdf,
    artifact_id = "EDFIG10A_B_B",
    role = "external_final_figure_pdf",
    description = "Extended Data Figure 10 final PDF: externally assembled ichorCNA component figure preserved as a frozen manuscript artifact.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )
} else {
  warning(
    "No frozen Extended Data Figure 10 PDF was found. The in-repo support CSVs were generated, ",
    "but the final panel assembly depends on external ichorCNA component plots."
  )
}

## ---- 7. Per-patient plots (Tumour fraction + FGA proxy) ---------------------
plot_input <- cfDNA_df %>%
  filter(!is.na(Sample_Date), is.finite(as.numeric(Sample_Date)))

skipped_plot_rows <- cfDNA_df %>%
  filter(is.na(Sample_Date) | !is.finite(as.numeric(Sample_Date)))
if (nrow(skipped_plot_rows) > 0) {
  write_csv(skipped_plot_rows, file.path(outdir, "Subclonal_evolution_skipped_plot_rows_missing_dates.csv"))
}

plot_list <- plot_input %>%
  mutate(
    is_relapse = factor(is_relapse,
                        levels = c(FALSE, TRUE),
                        labels = c("Non-relapse", "Relapse"))
  ) %>%
  dplyr::group_nest(Patient) %>%        # <-- use dplyr::group_nest()
  mutate(
    plot = purrr::map2(
      data, Patient,
      ~ ggplot(.x, aes(x = Sample_Date)) +
        geom_line(aes(y = Tumor_Fraction),
                  colour = "steelblue", linewidth = 1) +
        geom_point(aes(y = Tumor_Fraction,
                       shape = is_relapse,
                       fill  = is_relapse),
                   colour = "black", size = 3) +
        geom_line(aes(y = FGA_proxy),
                  colour = "firebrick", linewidth = 1,
                  linetype = "dashed") +
        geom_point(aes(y = FGA_proxy),
                   colour = "firebrick", size = 3,
                   shape = 23, fill = "white") +
        scale_y_continuous(
          name     = "Tumour fraction (ichorCNA)",
          sec.axis = sec_axis(~., name = "FGA proxy (0–1)")
        ) +
        scale_shape_manual(values = c(21, 24)) +
        scale_fill_manual(values = c("grey80", "firebrick")) +
        labs(
          title = .y,                # patient_id now comes from .y
          x     = "Sample date"
        ) +
        theme_bw(base_size = 11) +
        theme(
          legend.position = "bottom",
          plot.title      = element_text(face = "bold")
        )
    )
  ) %>%
  filter(map_int(data, nrow) > 0)

pdf(out_plot_pdf, width = 7.5, height = 5.5)
if (nrow(plot_list) > 0) walk(plot_list$plot, print)
dev.off()


### Now get summary text
# ---- 8. Summary text for manuscript ----------------------------------------

# 1) Basic counts
df_keep        <- cfDNA_df
n_samples      <- nrow(df_keep)
patients       <- sort(unique(df_keep$Patient))
n_patients     <- length(patients)

# 2) Samples per patient
spp            <- df_keep %>%
  group_by(Patient) %>%
  summarise(n = dplyr::n(), .groups = "drop")
avg_spp        <- median(spp$n, na.rm = TRUE)
range_spp      <- range(spp$n, na.rm = TRUE)

# 3) Emergent‐CNA patients
emergent_pts   <- emergent_tbl %>% distinct(Patient) %>% pull(Patient)
n_emergent     <- length(emergent_pts)
pct_emergent   <- n_emergent / n_patients * 100

# 4) Days before progression for emergent events
#    (we look up Num_days_to_closest_relapse from cfDNA_df for those samples)
days_vec       <- cfDNA_df %>%
  filter(Sample %in% emergent_tbl$Sample) %>%
  pull(Num_days_to_closest_relapse)

mean_days      <- mean(days_vec, na.rm = TRUE)
range_days     <- range(days_vec, na.rm = TRUE)
iqr_days       <- IQR(days_vec, na.rm = TRUE)

# 5) Print formatted sentence
summary_text <- sprintf(
  "There were %d samples from %d patients with both baseline and progression cfDNA samples (median %.1f samples per patient, range %d–%d). Of these, %d/%d (%.1f%%) showed evidence of subclonal evolution via CNA changes, with emergent CNAs detected on average %.1f days before progression (range %d–%d; IQR %.1f days).",
  n_samples, n_patients,
  avg_spp, range_spp[1], range_spp[2],
  n_emergent, n_patients, pct_emergent,
  mean_days, range_days[1], range_days[2], iqr_days
)

cat(summary_text, "\n")

## Additional summaries

# 1) Add days-to-progression to emergent_tbl
emergent_tbl_days <- emergent_tbl %>%
  left_join(
    cfDNA_df %>% 
      select(Patient, Sample, Num_days_to_closest_relapse),
    by = c("Patient", "Sample")
  )

# 2) Overall summary
overall_days <- emergent_tbl_days$Num_days_to_closest_relapse
overall_summary <- tibble(
  n_events   = nrow(emergent_tbl_days),
  mean_days  = mean(overall_days, na.rm = TRUE),
  median_days= median(overall_days, na.rm = TRUE),
  range_days = paste(range(overall_days, na.rm = TRUE), collapse = "–"),
  iqr_days   = IQR(overall_days, na.rm = TRUE)
)

print(overall_summary)

# 3) Per–event type summary (optional)
per_event_summary <- emergent_tbl_days %>%
  group_by(Event) %>%
  summarise(
    n_events    = dplyr::n(),
    mean_days   = mean(Num_days_to_closest_relapse, na.rm = TRUE),
    median_days = median(Num_days_to_closest_relapse, na.rm = TRUE),
    range_days  = paste(range(Num_days_to_closest_relapse, na.rm = TRUE), collapse = "–"),
    iqr_days    = IQR(Num_days_to_closest_relapse, na.rm = TRUE),
    .groups     = "drop"
  )

print(per_event_summary)


### See second closest draw time

# 1) Identify patients of interest (e.g. those with emergent events)
patients_of_interest <- emergent_tbl %>%
  distinct(Patient) %>%
  pull(Patient)

# 2) For each patient, find the second‐closest cfDNA draw to progression
second_closest <- cfDNA_df %>%
  filter(Patient %in% patients_of_interest) %>%
  # keep only non‐relapse draws (so relapse itself with zero days doesn’t count)
  filter(Num_days_to_closest_relapse > 0) %>%
  group_by(Patient) %>%
  filter(timepoint_info != "Relapse") %>%
  arrange(Num_days_to_closest_relapse) %>%
  slice(1) %>%  # closest distance
  ungroup() %>%
  select(
    Patient,
    Sample_ID,
    Sample_Date       = Date_of_sample_collection,
    Days_before_relapse = Num_days_to_closest_relapse,
    Timepoint         = timepoint_info
  )

# Compute median, range, and IQR:
second_closest %>%
  summarise(
    median_days = median(Days_before_relapse, na.rm = TRUE),
    range_days  = paste(range(Days_before_relapse, na.rm = TRUE), collapse = "–"),
    iqr_days    = IQR(Days_before_relapse, na.rm = TRUE)
  ) %>%
  print()

# If you prefer separate base‐R values:
days <- second_closest$Days_before_relapse
cat("Median:", median(days, na.rm = TRUE), "\n")
cat("Range:", paste(range(days, na.rm = TRUE), collapse = "–"), "\n")
cat("IQR:", IQR(days, na.rm = TRUE), "\n")

print(second_closest)



### Now see how early tumor fraction rose 

# 1) Patients who progressed (same as before)
pts_prog <- emergent_tbl %>% distinct(Patient) %>% pull(Patient)

# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # nadir TF (lowest non‐relapse TF)
    tf_nadir    = min(Tumor_Fraction[!is_relapse], na.rm = TRUE),
    # days before progression when the nadir occurred
    days_nadir  = Num_days_to_closest_relapse[which.min(ifelse(is_relapse, Inf, Tumor_Fraction))],
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %d patients who showed novel CNAs at progression, median tumour fraction rose from %.1f%% (IQR %.1f–%.1f%%) at diagnosis to %.1f%% (IQR %.1f–%.1f%%) at relapse.  From each patient’s nadir (median %.1f%% [IQR %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (IQR %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days; IQR %.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_q[1], baseline_q[2],
  median(tf_summary$tf_relapse, na.rm = TRUE)  * 100, relapse_q[1], relapse_q[2],
  median(tf_summary$tf_nadir, na.rm = TRUE)    * 100, quantile(tf_summary$tf_nadir, .25, na.rm = TRUE) * 100, quantile(tf_summary$tf_nadir, .75, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise, na.rm = TRUE)     * 100, rise_q[1], rise_q[2],
  median(tf_summary$days_nadir, na.rm = TRUE), min(tf_summary$days_nadir, na.rm = TRUE), max(tf_summary$days_nadir, na.rm = TRUE), IQR(tf_summary$days_nadir, na.rm = TRUE)
)

cat(summary_sentence, "\n")


# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %d patients who showed novel CNAs at progression, median tumour fraction rose from %.1f%% (range %.1f–%.1f%%) at diagnosis to %.1f%% (range %.1f–%.1f%%) at relapse. From each patient’s nadir (median %.1f%% [range %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (range %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_range[1], baseline_range[2],
  median(tf_summary$tf_relapse,  na.rm = TRUE) * 100, relapse_range[1],  relapse_range[2],
  median(tf_summary$tf_nadir,    na.rm = TRUE) * 100, min(tf_summary$tf_nadir, na.rm = TRUE) * 100, max(tf_summary$tf_nadir, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise,     na.rm = TRUE) * 100, rise_range[1], rise_range[2],
  median(tf_summary$days_nadir,  na.rm = TRUE), days_range[1], days_range[2]
)

cat(summary_sentence, "\n")

write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_new_CNA_patients.csv"))


## Do now for just those who didn't show high risk CNAs 
# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(!Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # nadir TF (lowest non‐relapse TF)
    tf_nadir    = min(Tumor_Fraction[!is_relapse], na.rm = TRUE),
    # days before progression when the nadir occurred
    days_nadir  = Num_days_to_closest_relapse[which.min(ifelse(is_relapse, Inf, Tumor_Fraction))],
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %d patients who did not show novel CNAs at progression, median tumour fraction rose from %.1f%% (IQR %.1f–%.1f%%) at diagnosis to %.1f%% (IQR %.1f–%.1f%%) at relapse.  From each patient’s nadir (median %.1f%% [IQR %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (IQR %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days; IQR %.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_q[1], baseline_q[2],
  median(tf_summary$tf_relapse, na.rm = TRUE)  * 100, relapse_q[1], relapse_q[2],
  median(tf_summary$tf_nadir, na.rm = TRUE)    * 100, quantile(tf_summary$tf_nadir, .25, na.rm = TRUE) * 100, quantile(tf_summary$tf_nadir, .75, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise, na.rm = TRUE)     * 100, rise_q[1], rise_q[2],
  median(tf_summary$days_nadir, na.rm = TRUE), min(tf_summary$days_nadir, na.rm = TRUE), max(tf_summary$days_nadir, na.rm = TRUE), IQR(tf_summary$days_nadir, na.rm = TRUE)
)

cat(summary_sentence, "\n")


# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %d patients who did not show novel CNAs at progression, median tumour fraction rose from %.1f%% (range %.1f–%.1f%%) at diagnosis to %.1f%% (range %.1f–%.1f%%) at relapse. From each patient’s nadir (median %.1f%% [range %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (range %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_range[1], baseline_range[2],
  median(tf_summary$tf_relapse,  na.rm = TRUE) * 100, relapse_range[1],  relapse_range[2],
  median(tf_summary$tf_nadir,    na.rm = TRUE) * 100, min(tf_summary$tf_nadir, na.rm = TRUE) * 100, max(tf_summary$tf_nadir, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise,     na.rm = TRUE) * 100, rise_range[1], rise_range[2],
  median(tf_summary$days_nadir,  na.rm = TRUE), days_range[1], days_range[2]
)

cat(summary_sentence, "\n")

write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_new_CNA_patients.csv"))




### redo for all patients 
### Now see how early tumor fraction rose 

# 1) Patients who progressed (same as before)
pts_prog <- cfDNA_df %>% distinct(Patient) %>% pull(Patient)

# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # nadir TF (lowest non‐relapse TF)
    tf_nadir    = min(Tumor_Fraction[!is_relapse], na.rm = TRUE),
    # days before progression when the nadir occurred
    days_nadir  = Num_days_to_closest_relapse[which.min(ifelse(is_relapse, Inf, Tumor_Fraction))],
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %d patients with baseline and cfDNA samples, median tumour fraction rose from %.1f%% (IQR %.1f–%.1f%%) at diagnosis to %.1f%% (IQR %.1f–%.1f%%) at relapse. From each patient’s nadir (median %.1f%% [IQR %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (IQR %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days; IQR %.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_q[1], baseline_q[2],
  median(tf_summary$tf_relapse, na.rm = TRUE)  * 100, relapse_q[1], relapse_q[2],
  median(tf_summary$tf_nadir, na.rm = TRUE)    * 100,
  quantile(tf_summary$tf_nadir, .25, na.rm = TRUE) * 100,
  quantile(tf_summary$tf_nadir, .75, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise, na.rm = TRUE)     * 100, rise_q[1], rise_q[2],
  median(tf_summary$days_nadir, na.rm = TRUE), min(tf_summary$days_nadir, na.rm = TRUE),
  max(tf_summary$days_nadir, na.rm = TRUE), IQR(tf_summary$days_nadir, na.rm = TRUE)
)

cat(summary_sentence, "\n")

# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)

summary_sentence <- sprintf(
  "In the %.0f patients with baseline and progression cfDNA samples, median tumour fraction rose from %.1f%% (range %.1f–%.1f%%) at diagnosis to %.1f%% (range %.1f–%.1f%%) at relapse. From each patient’s nadir (median %.1f%% [range %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (range %.1f–%.1f%%), with the nadir detected a median of %.0f days before progression (range %.0f–%.0f days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline, na.rm = TRUE) * 100, baseline_range[1], baseline_range[2],
  median(tf_summary$tf_relapse,  na.rm = TRUE) * 100, relapse_range[1],  relapse_range[2],
  median(tf_summary$tf_nadir,    na.rm = TRUE) * 100,
  min(tf_summary$tf_nadir, na.rm = TRUE) * 100, max(tf_summary$tf_nadir, na.rm = TRUE) * 100,
  median(tf_summary$tf_rise,     na.rm = TRUE) * 100, rise_range[1], rise_range[2],
  median(tf_summary$days_nadir,  na.rm = TRUE), days_range[1], days_range[2]
)


cat(summary_sentence, "\n")

write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_all_CNA_patients.csv"))
