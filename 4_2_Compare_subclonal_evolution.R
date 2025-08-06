###############################################################################
# 4_2_Compare_subclonal_evolution.R
# Quick copy-number–based assessment of emergent subclones in longitudinal
# cfDNA WGS MRD samples (30×).  Produces per-patient plots + emergent-CNA table.
###############################################################################

## ---- 0. USER SETTINGS -------------------------------------------------------
in_rds        <- "Jan2025_exported_data/All_feature_data_May2025.rds"      # or "All_feature_data.csv"
out_plot_pdf  <- "Subclonal_evolution_plots.pdf"
out_events_csv<- "Emergent_CNA_events.csv"
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

## ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

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
  mutate(across(all_of(event_cols), ~ as.numeric(as.character(.x)))) %>%
  rowwise() %>%
  mutate(FGA_proxy = sum(c_across(all_of(event_cols)), na.rm = TRUE) /
           length(event_cols)) %>%
  ungroup()

## ---- 6. Emergent CNA table --------------------------------------------------
emergent_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - first(.x), .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(Event   = str_remove(Event, "^delta_"),
         Emergent = Delta == 1) %>%
  filter(Emergent) %>%
  select(Patient, Sample, Event) %>%
  arrange(Patient, Event)

write_csv(emergent_tbl, out_events_csv)

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

write_csv(gains_tbl, out_events_csv)

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

write_csv(losses_tbl, out_events_csv)

## ---- 7. Per-patient plots (Tumour fraction + FGA proxy) ---------------------
plot_list <- cfDNA_df %>%
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
  )

pdf(out_plot_pdf, width = 7.5, height = 5.5)
walk(plot_list$plot, print)
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
avg_spp        <- median(spp$n)
range_spp      <- range(spp$n)

# 3) Emergent‐CNA patients
emergent_pts   <- emergent_tbl %>% distinct(Patient) %>% pull(Patient)
n_emergent     <- length(emergent_pts)
pct_emergent   <- n_emergent / n_patients * 100

# 4) Days before progression for emergent events
#    (we look up Num_days_to_closest_relapse from cfDNA_df for those samples)
days_vec       <- cfDNA_df %>%
  filter(Sample %in% emergent_tbl$Sample) %>%
  pull(Num_days_to_closest_relapse)

mean_days      <- mean(days_vec)
range_days     <- range(days_vec)
iqr_days       <- IQR(days_vec)

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
  mean_days  = mean(overall_days),
  median_days= median(overall_days),
  range_days = paste(range(overall_days), collapse = "–"),
  iqr_days   = IQR(overall_days)
)

print(overall_summary)

# 3) Per–event type summary (optional)
per_event_summary <- emergent_tbl_days %>%
  group_by(Event) %>%
  summarise(
    n_events    = dplyr::n(),
    mean_days   = mean(Num_days_to_closest_relapse),
    median_days = median(Num_days_to_closest_relapse),
    range_days  = paste(range(Num_days_to_closest_relapse), collapse = "–"),
    iqr_days    = IQR(Num_days_to_closest_relapse),
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
    median_days = median(Days_before_relapse),
    range_days  = paste(range(Days_before_relapse), collapse = "–"),
    iqr_days    = IQR(Days_before_relapse)
  ) %>%
  print()

# If you prefer separate base‐R values:
days <- second_closest$Days_before_relapse
cat("Median:", median(days), "\n")
cat("Range:", paste(range(days), collapse = "–"), "\n")
cat("IQR:", IQR(days), "\n")

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
  "In the %d patients who progressed, median tumour fraction rose from %.1f%% (IQR %.1f–%.1f%%) at diagnosis to %.1f%% (IQR %.1f–%.1f%%) at relapse.  From each patient’s nadir (median %.1f%% [IQR %.1f–%.1f%%]), tumour fraction increased by a median of %.1f%% (IQR %.1f–%.1f%%), with the nadir detected a median of %d days before progression (range %d–%d days; IQR %d days).",
  nrow(tf_summary),
  median(tf_summary$tf_baseline) * 100, baseline_q[1], baseline_q[2],
  median(tf_summary$tf_relapse)  * 100, relapse_q[1], relapse_q[2],
  median(tf_summary$tf_nadir)    * 100, quantile(tf_summary$tf_nadir, .25) * 100, quantile(tf_summary$tf_nadir, .75) * 100,
  median(tf_summary$tf_rise)     * 100, rise_q[1], rise_q[2],
  median(tf_summary$days_nadir), min(tf_summary$days_nadir), max(tf_summary$days_nadir), IQR(tf_summary$days_nadir)
)

cat(summary_sentence, "\n")
