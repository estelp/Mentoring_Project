## Install LEA package

install.packages("devtools")
devtools::install_github("bcm-uga/LEA")

# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses"
if (!dir.exists(lea_dir)) {
  dir.create(lea_dir, recursive = TRUE)
}

# Define the path for the input VCF file
vcf_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/vcf_filtered/filtered_snps.vcf"

# Check if the VCF file exists
if (!file.exists(vcf_file)) {
  stop(paste("Error: The specified VCF file does not exist:", vcf_file))
}

# Copy the VCF file to the working directory for LEA analysis
vcf_file_copy <- file.path(lea_dir, basename(vcf_file))
file.copy(vcf_file, vcf_file_copy, overwrite = TRUE)

# Read sample information from the VCF file to retain identity
vcf_data <- read.vcfR(vcf_file_copy)
sample_ids <- colnames(vcf_data@gt)[-1]  # Remove the chromosome information column

# Verify that the sample IDs were correctly retrieved
cat("Sample IDs extracted from the VCF file:\n")
cat(sample_ids, "\n")

# Convert the VCF file to GENO format while keeping sample identities
cat("Converting VCF file to GENO...\n")
output <- vcf2geno(vcf_file_copy)

# Rename and move the output files to the LEA directory
geno_file <- file.path(lea_dir, paste0(tools::file_path_sans_ext(basename(vcf_file)), ".geno"))
vcfsnp_file <- file.path(lea_dir, paste0(tools::file_path_sans_ext(basename(vcf_file)), ".vcfsnp"))
removed_file <- file.path(lea_dir, paste0(tools::file_path_sans_ext(basename(vcf_file)), ".removed"))

# Save the sample IDs in the GENO file
cat("Sample IDs saved in the GENO file.\n")
write.table(sample_ids, file = paste0(tools::file_path_sans_ext(geno_file), ".samples"), row.names = FALSE, col.names = FALSE)

# Display the output information
cat("Conversion completed. Files have been saved in:\n")
cat("- GENO file:", geno_file, "\n")
cat("- SNP information:", vcfsnp_file, "\n")
cat("- Removed lines:", removed_file, "\n")
cat("- Sample IDs saved in the file:", paste0(tools::file_path_sans_ext(geno_file), ".samples"), "\n")

# Load the GENO data after conversion
cat("Loading the GENO file...\n")
genotypes_matrix <- read.geno(geno_file)

# Check a preview of the genotype matrix
cat("Preview of the genotype matrix:\n")
head(genotypes_matrix)

# Convert the genotype matrix to numeric format, replacing 'NA' and empty cells with '9'
genotypes_matrix_numeric <- apply(genotypes_matrix, 2, function(x) {
  x[is.na(x) | x == ""] <- "9"  # Replace NA or empty with 9
  as.numeric(x)
})

# Verify that the conversion was done correctly
cat("Preview of the genotype matrix converted to numeric:\n")
head(genotypes_matrix_numeric)

# Save the genotype matrix in LFMM format
lfmm_file <- file.path(lea_dir, "genotypes.lfmm")
write.lfmm(genotypes_matrix_numeric, lfmm_file)

# Save the sample IDs in the LFMM file
cat("Sample IDs added to the LFMM file:\n")
write.table(sample_ids, paste0(tools::file_path_sans_ext(lfmm_file), "_samples.txt"), row.names = FALSE, col.names = FALSE)

# Optional: Additional analyses like genetic structure analysis with snmf
project <- snmf(file.path(lea_dir, "filtered_snps.geno"),
                K = 1:20, 
                entropy = TRUE, 
                repetitions = 20, 
                project = "new")

# Specify an output file for the PNG plot
cross_entropy_plot_file <- file.path(lea_dir, "cross_entropy_plot.png")

# Save the plot to a PNG file
png(cross_entropy_plot_file, width = 1800, height = 800)

# Plot the cross-entropy criterion without specifying xlab, ylab, and main (to avoid conflicts)
plot(project, 
     col = "blue", pch = 19, cex = 1.2)

# Close the PNG file to save the image
dev.off()

cat("The cross-entropy plot has been saved to:", cross_entropy_plot_file, "\n")

# Select the best run for K = 7 clusters
best_run <- which.min(cross.entropy(project, K = 7))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 7, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Define colors for clusters
cluster_colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange")

# Load sample IDs from the .samples file
samples_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses/filtered_snps.samples"
sample_ids <- read.table(samples_file, header = FALSE, stringsAsFactors = FALSE)[, 1]

# Create a larger output file for the figure
# Adjust the window size for the figure to be sufficiently large
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses/ancestry_matrix.png", width = 1800, height = 800)

# Create the ancestry matrix plot
barchart(project, K = 7, run = best_run,
         border = NA, space = 0,
         col = cluster_colors,
         xlab = "Individuals",       # Label for the x-axis
         ylab = "Proportions of ancestry",  # Label for the y-axis
         main = "Ancestry Matrix")

# Add the sample IDs to the x-axis
axis(1, at = 1:length(sample_ids), labels = sample_ids, las = 2, cex.axis = 0.5)  # Reduce the font size for the x-axis

# Close the PNG file to save the image
dev.off()

cat("The ancestry matrix has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses/ancestry_matrix.png\n")

