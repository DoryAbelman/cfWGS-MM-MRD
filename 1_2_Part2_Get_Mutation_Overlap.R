# =============================================================================
# Script: 1_2_Part2_Get_Mutation_Overlap.R
# =============================================================================
# ──────────────────────────────────────────────────────────────────────────────
# Title       : Mutation Overlap Analysis (BM vs. PB cfDNA)  
# Author      : Dory Abelman  
# Date        : 2025-06-30  
#  
# Description :  
#   1. Reads in bone-marrow and blood cfDNA MAFs (maftools::read.maf)  
#   2. Merges variant calls, creates per-mutation unique IDs  
#   3. Builds wide table for presence/absence by sample type  
#   4. For each patient:  
#        • draws a Venn diagram of BM vs. cfDNA calls  
#        • computes % overlap and plots barplot  
#  
# Usage       :  
#   Rscript 1_2_Part2_Get_Mutation_Overlap.R  
#   — or —  
#   source("1_2_Part2_Get_Mutation_Overlap.R")  
#  
# Dependencies:  
#   library(maftools)    # read.maf()  
#   library(dplyr)       # data wrangling  
#   library(tidyr)       # pivot_wider()  
#   library(ggplot2)     # barplot  
#   library(scales)      # optional formatting  
#   library(grid)        # grid graphics  
#   library(VennDiagram) # draw.pairwise.venn()  
#  
# Input files :  
#   - combined_maf_temp_bm_Jan2025.maf  
#   - combined_maf_temp_blood_Jan2025.maf  
#   - metada_df_mutation_comparison (data frame in workspace)  
#  
# Output files:  
#   - VennDiagram_<Patient>_May2024_updatedcolors.png  
#   - percent_overlap_barplot.png  
#  
# Notes       :  
#   • Ensure your working directory is set appropriately or use full paths  
#   • Make sure metada_df_mutation_comparison is loaded and correctly spelled  
# ──────────────────────────────────────────────────────────────────────────────


### Load libraries 
library(maftools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)
library(VennDiagram)
library(viridis)




# Read the MAF file using read.maf
df_bm <- readRDS("combined_maf_bm_dx.rds")
df_blood <- readRDS("combined_maf_blood_all_muts.rds")

## Remove rsids to match new filtering 
# filter out all rows with an RSID
df_blood <- df_blood %>%
  filter(is.na(dbSNP_RS) | dbSNP_RS == "" | dbSNP_RS == "." |
           !grepl("^rs", dbSNP_RS, ignore.case = TRUE))

df_bm <- df_bm %>%
  filter(is.na(dbSNP_RS) | dbSNP_RS == "" | dbSNP_RS == "." |
           !grepl("^rs", dbSNP_RS, ignore.case = TRUE))

# now combine if you like
combined_maf <- bind_rows(df_bm, df_blood)

rm(df_bm)
rm(df_blood)

## Venn diagram between mutations at diagnosis 
# Create unique identifiers for each mutation
combined_maf <- combined_maf %>%
  mutate(
    Unique_ID_patient = paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Patient, sep = ":"),
    Mutation_ID = paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")
  )

# Create a column to indicate presence (1) of the mutation
combined_maf <- combined_maf %>%
  mutate(Presence = 1)

## Re-join with info in case issues 
# Load in the patient info 
metada_df_mutation_comparison <- read_csv("combined_clinical_data_updated_April2025.csv")

# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Tumor_Sample_Barcode = Bam %>%
           # Remove _PG or _WG
           str_remove_all("_PG|_WG") %>%
           # Remove anything after ".filter", ".ded", or ".recalibrate"
           str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""))

# Add the Bam column to combined_maf
combined_maf <- combined_maf %>%
  mutate(Bam = paste0(Tumor_Sample_Barcode, ".filter.deduped.recalibrated.bam"))

# Modify the specific Tumor_Sample_Barcode in combined_maf with error
combined_maf <- combined_maf %>%
  mutate(Tumor_Sample_Barcode = ifelse(Tumor_Sample_Barcode == "TFRIM4_0189_Bm_P_ZC-02", 
                                       paste0(Tumor_Sample_Barcode, "-01-O-DNA"), 
                                       Tumor_Sample_Barcode))

combined_maf <- combined_maf %>%
  select(
    -Patient,
    -Date_of_sample_collection,
    -Sample_type,
    -Timepoint,
    -Study,
    -Sample_ID,
    -timepoint_info,
    -Relapsed,
    -Num_days_to_closest_relapse_absolute,
    -Num_days_to_closest_relapse_non_absolute,
    -Num_days_to_closest_relapse
  )

# Join with the metadata dataframe
combined_maf <- left_join(combined_maf %>% select(-Bam), metada_df_mutation_comparison, by = "Tumor_Sample_Barcode")

## First filter to just diagnosis timepoints
combined_maf <- combined_maf %>% filter(timepoint_info %in% c("Diagnosis", "Baseline"))

## Add cohort df 
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

# find patients in cohort_df not in combined_maf
missing_patients <- setdiff(combined_maf$Patient, combined_maf$Patient)

if (length(missing_patients) == 0) {
  message("✔ All patients in cohort_df are represented in combined_maf.")
} else {
  message("⚠️ These patients are missing from combined_maf: ", 
          paste(missing_patients, collapse = ", "))
}

# Pivot the dataframe to wide format
df_wide <- pivot_wider(
  combined_maf,
  id_cols = c("Hugo_Symbol", 'Patient', "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Unique_ID_patient"),
  names_from = Sample_type,
  values_from = c(VAF, t_depth, Presence),
  values_fn = {max}
)

## Filter for only pairs
patients_with_required_samples <- metada_df_mutation_comparison %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  group_by(Patient, Timepoint) %>%
  filter(any(Sample_type == "Blood_plasma_cfDNA") & any(Sample_type %in% c("BM_cells"))) %>%
  distinct(Patient, Timepoint)

# Join the filtered patients with the original df_wide to get the final filtered dataframe
df_wide <- df_wide %>%
  inner_join(patients_with_required_samples)


## Make plot
patients <- levels(as.factor(patients_with_required_samples$Patient))

## Now create a venn diagram for each patient 
# Loop over all patients
library(scales)
library(grid)
library(VennDiagram)

# Loop through each patient
for (i in 1:length(patients)) {
  
  # Subset the data for the current patient
  df_patient <- df_wide %>% filter(Patient == patients[i])
  
  # Get the mutations for each set
  bm_mutations <- df_patient %>% filter(Presence_BM_cells == 1) %>% pull(Unique_ID_patient)
  cf_mutations <- df_patient %>% filter(Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient)
  both_mutations <- df_patient %>% filter(Presence_BM_cells == 1 & Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient)
  
  # Create a new plot with a white background
  png(filename = paste0("VennDiagram_", patients[i], "_June2025_updatedcolors.png"), width = 2400, height = 2400, res = 300)
  
  # Clear the plot and set background to white
  grid.newpage()
  pushViewport(viewport(gp = gpar(fill = "white")))
  
  # Generate the Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(bm_mutations),
    area2 = length(cf_mutations),
    cross.area = length(intersect(bm_mutations, cf_mutations)),
    category = c("BM cells", "Blood cfDNA"),
    fill = c("#008080", "#FF6347"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("#008080", "#FF6347"),
    margin = 0.1
  )
  
  # Add title using grid graphics
  grid.text(paste0("Patient: ", patients[i]), x = 0.5, y = 0.95, gp = gpar(fontsize = 20))
  grid.draw(venn.plot)
  dev.off()
}


## Make barplot showing percent overlap 

# Create a dataframe to store the overlap percentages
overlap_data <- data.frame(Patient = character(), Percent_Overlap = numeric(), stringsAsFactors = FALSE)

# Loop through each patient
for (i in 1:length(patients)) {
  
  # Subset the data for the current patient
  df_patient <- df_wide %>% filter(Patient == patients[i])
  
  # Get the mutations for each set
  bm_mutations <- df_patient %>% filter(Presence_BM_cells == 1) %>% pull(Unique_ID_patient)
  cf_mutations <- df_patient %>% filter(Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient)
  
  # Calculate the number of overlapping mutations
  overlap_count <- length(intersect(bm_mutations, cf_mutations))
  
  # Calculate the total number of unique mutations
  total_unique_mutations <- length(unique(c(bm_mutations, cf_mutations)))
  
  # Calculate the percent overlap
  percent_overlap <- (overlap_count / total_unique_mutations) * 100
  
  # Store the results in the dataframe
  overlap_data <- rbind(overlap_data, data.frame(Patient = patients[i], Percent_Overlap = percent_overlap))
}

# Reorder the Patient factor levels based on Percent_Overlap in descending order
overlap_data <- overlap_data %>%
  arrange(desc(Percent_Overlap)) %>%
  mutate(Patient = factor(Patient, levels = unique(Patient)))

# Create the bar plot
overlap_plot <- ggplot(overlap_data, aes(x = Patient, y = Percent_Overlap)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Percent Overlap Between BM and PB cfDNA mutation calls",
    x = "Patient",
    y = "Percent Overlap"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

# Display the plot
print(overlap_plot)

# Save the plot
ggsave("percent_overlap_barplot_June2025.png", plot = overlap_plot, width = 6, height = 5, dpi = 500)


### Now get the overlap table and summarise stats 

# 1) Left-join the cohort assignments onto the overlap table
overlap_with_cohort <- overlap_data %>%
  left_join(cohort_df, by = "Patient")

# 2) How many patients have both baseline BM and cfDNA samples?
n_patients_overlap <- n_distinct(overlap_with_cohort$Patient)
message("Number of patients with both baseline BM and cfDNA: ", n_patients_overlap)

# 3) Overall summary of Percent_Overlap
overall_stats <- overlap_with_cohort %>%
  summarise(
    mean_overlap   = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap    = IQR(Percent_Overlap, na.rm = TRUE),
    min_overlap    = min(Percent_Overlap, na.rm = TRUE),
    max_overlap    = max(Percent_Overlap, na.rm = TRUE)
  )
print(overall_stats)

# 4) Summary by cohort (replace 'Cohort' with your actual cohort column name)
stats_by_cohort <- overlap_with_cohort %>%
  group_by(Cohort) %>%
  summarise(
    n               = dplyr::n(),
    mean_overlap    = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap  = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap     = IQR(Percent_Overlap, na.rm = TRUE),
    min_overlap     = min(Percent_Overlap, na.rm = TRUE),
    max_overlap     = max(Percent_Overlap, na.rm = TRUE)
  )
print(stats_by_cohort)


library(glue)
# Create descriptive sentences for each cohort
cohort_sentences <- stats_by_cohort %>%
  mutate(
    sentence = glue(
      "In the {Cohort} cohort, there were {n} matched baseline samples. ",
      "The percent overlap had a mean of {round(mean_overlap, 1)}% ",
      "(range {round(min_overlap, 1)}%–{round(max_overlap, 1)}%, ",
      "IQR {round(IQR_overlap, 1)}%)."
    )
  ) %>%
  pull(sentence)

# Print them
cat(cohort_sentences, sep = "\n")


### Now make figure 3C 
# 1. Sort overlap_data and add position for plotting
plot_df <- overlap_with_cohort %>%
  filter(Cohort == "Frontline") %>%
  arrange(Percent_Overlap) %>%
  mutate(
    Patient = factor(Patient, levels = Patient),      # keep order
    pos     = row_number()                            # y-position
  )

# 2. Calculate overall statistics
med  <- median(plot_df$Percent_Overlap, na.rm = TRUE)
iqrL <- quantile(plot_df$Percent_Overlap, 0.25, na.rm = TRUE)
iqrU <- quantile(plot_df$Percent_Overlap, 0.75, na.rm = TRUE)

# 3. Build lollipop plot
p_overlap <- ggplot(plot_df, aes(x = Percent_Overlap, y = Patient)) +
 # annotate(                         # IQR ribbon
 #   "rect", xmin = iqrL, xmax = iqrU,
 #   ymin = 0.5, ymax = nrow(plot_df) + 0.5,
 #   fill = "grey90"
 # ) +
  geom_vline(xintercept = med, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0, xend = Percent_Overlap, yend = Patient),
               colour = "grey65", size = 0.35) +
  geom_point(aes(colour = Percent_Overlap), size = 2) +
  scale_colour_viridis_c(
    option = "D", end = 0.9, name = "% overlap",
    guide  = guide_colourbar(barwidth = 0.4, barheight = 3)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Mutation overlap: BM ∩ cfDNA (%)",
    y = NULL,
    title    = "Patient-level overlap of baseline mutations (BM vs. cfDNA)",
    subtitle = glue::glue(
      "Median = {round(med,1)}%   |   IQR = {round(iqrL,1)}–{round(iqrU,1)}%"
    )
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title      = element_text(face = "bold", size = 9),
    plot.subtitle   = element_text(size = 7, margin = margin(b = 6)),
    axis.text.y     = element_text(size = 6),
    panel.grid.major.y = element_blank(),
    legend.position     = c(0.97, 0.05),     # inside bottom-right
    legend.justification = c(1, 0),
    legend.background    = element_rect(
      fill = scales::alpha("white", 0.7), colour = NA
    ),
    legend.key.width  = unit(0.3, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.title      = element_text(size = 6),
    legend.text       = element_text(size = 6)
  )

# 4. Save for Figure 3C
ggsave(
  filename = file.path(outdir, "Fig3C_mutation_overlap_lollipop.png"),
  plot     = p_overlap,
  width    = 4,   # a bit narrower
  height   = 5,
  dpi      = 600
)


# Export overlap_with_cohort table
write.csv(overlap_with_cohort, "mutation_overlap_with_cohort.csv", row.names = FALSE)
# or TSV:
# write.table(overlap_with_cohort, "overlap_with_cohort.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Export overall stats
write.csv(overall_stats, "overall_percent_overlap_stats.csv", row.names = FALSE)

# Export cohort-specific summary stats
write.csv(stats_by_cohort, "percent_overlap_stats_by_cohort.csv", row.names = FALSE)
