# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Define working directories
input_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/Others_stats"
output_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/Others_stats/Plots"
summary_file <- file.path(output_dir, "summary_statistics.txt")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory
setwd(input_dir)

# Open summary file for writing
summary_conn <- file(summary_file, open = "w")

# Variant-based statistics

# Variant quality
var_qual <- read_delim(file.path(input_dir, "output.lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

p <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Variant Quality Distribution")
ggsave(file.path(output_dir, "variant_quality_density.png"), plot = p)

summary_var_qual <- summary(var_qual$qual)
writeLines("### Variant Quality Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_qual), summary_conn)

# Variant mean depth
var_depth <- read_delim(file.path(input_dir, "output.ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

p <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  xlim(0, 100) +
  ggtitle("Mean Depth per Variant")
ggsave(file.path(output_dir, "variant_depth_density.png"), plot = p)

summary_var_depth <- summary(var_depth$mean_depth)
writeLines("\n### Variant Mean Depth Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_depth), summary_conn)

# Variant missingness
var_miss <- read_delim(file.path(input_dir, "output.lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

p <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Missingness per Variant")
ggsave(file.path(output_dir, "variant_missingness_density.png"), plot = p)

summary_var_miss <- summary(var_miss$fmiss)
writeLines("\n### Variant Missingness Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_miss), summary_conn)

# Minor allele frequency
var_freq <- read_delim(file.path(input_dir, "output.frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# Calculate minor allele frequency
var_freq <- var_freq %>%
  mutate(maf = pmin(as.numeric(a1), as.numeric(a2), na.rm = TRUE))

p <- ggplot(var_freq, aes(maf)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Minor Allele Frequency Distribution")
ggsave(file.path(output_dir, "maf_density.png"), plot = p)

summary_var_freq <- summary(var_freq$maf)
writeLines("\n### Minor Allele Frequency Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_freq), summary_conn)

# Individual-based statistics

# Mean depth per individual
ind_depth <- read_delim(file.path(input_dir, "output.idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

p <- ggplot(ind_depth, aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Mean Depth per Individual")
ggsave(file.path(output_dir, "individual_depth_histogram.png"), plot = p)

summary_ind_depth <- summary(ind_depth$depth)
writeLines("\n### Individual Mean Depth Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_depth), summary_conn)

# Missing data per individual
ind_miss <- read_delim(file.path(input_dir, "output.imiss"), delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

p <- ggplot(ind_miss, aes(fmiss)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Missing Data per Individual")
ggsave(file.path(output_dir, "individual_missing_data_histogram.png"), plot = p)

summary_ind_miss <- summary(ind_miss$fmiss)
writeLines("\n### Individual Missing Data Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_miss), summary_conn)

# Heterozygosity and inbreeding coefficient
ind_het <- read_delim(file.path(input_dir, "output.het"), delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

p <- ggplot(ind_het, aes(f)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Inbreeding Coefficient Distribution")
ggsave(file.path(output_dir, "inbreeding_coefficient_histogram.png"), plot = p)

summary_ind_het <- summary(ind_het$f)
writeLines("\n### Inbreeding Coefficient Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_het), summary_conn)

# Close summary file
close(summary_conn)

message("All plots and summary statistics have been saved to ", output_dir)

