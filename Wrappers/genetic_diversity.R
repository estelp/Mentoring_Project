# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/2"
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
                K = 1:2, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 2 clusters
best_run <- which.min(cross.entropy(project, K = 2))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 2, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 2]  # Select cluster 2 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/2/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 2, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/2/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/3"
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
                K = 1:3, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 3 clusters
best_run <- which.min(cross.entropy(project, K = 3))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 3, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 3]  # Select cluster 3 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/3/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 3, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/3/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/4"
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
                K = 1:4, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 4 clusters
best_run <- which.min(cross.entropy(project, K = 4))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 4, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 4]  # Select cluster 4 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/4/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 4, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/4/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/5"
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
                K = 1:5, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 5 clusters
best_run <- which.min(cross.entropy(project, K = 5))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 5, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 5]  # Select cluster 5 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/5/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 5, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/5/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/6"
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
                K = 1:6, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 6 clusters
best_run <- which.min(cross.entropy(project, K = 6))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 6, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 6]  # Select cluster 6 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/6/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 6, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/6/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/7"
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
                K = 1:7, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 7]  # Select cluster 9 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/7/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 7, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/7/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/8"
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
                K = 1:8, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 8 clusters
best_run <- which.min(cross.entropy(project, K = 8))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 8, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 8]  # Select cluster 9 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/8/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 8, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/8/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")




# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/9"
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
                K = 1:9, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 9 clusters
best_run <- which.min(cross.entropy(project, K = 9))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 9, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 9]  # Select cluster 9 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/9/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 9, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red", "green"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/9/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red", "green")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")



# Load necessary libraries
library(LEA)
library(vcfR)

# Define the directory for LEA analysis
lea_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/10"
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
                K = 1:10, 
                entropy = TRUE, 
                repetitions = 10, 
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

# Select the best run for K = 10 clusters
best_run <- which.min(cross.entropy(project, K = 10))

# Extract the ancestry matrix (Q-matrix)
cat("Extracting the ancestry matrix...\n")
ancestry_matrix <- Q(project, K = 10, run = best_run)

# Add sample IDs to the ancestry matrix
ancestry_matrix_with_ids <- cbind(sample_ids, ancestry_matrix)

# Save the ancestry matrix with sample IDs
ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_with_ids.txt")
write.table(ancestry_matrix_with_ids, ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix with sample IDs saved:", ancestry_matrix_file, "\n")

# Sort individuals by degree of ancestry in a selected cluster
ancestry_cluster <- ancestry_matrix[, 10]  # Select cluster 10 for sorting
order_indices <- order(ancestry_cluster, decreasing = TRUE)  # Order by decreasing ancestry
sorted_sample_ids <- sample_ids[order_indices]
sorted_ancestry_matrix <- ancestry_matrix[order_indices, ]

# Save the sorted ancestry matrix
sorted_ancestry_matrix_with_ids <- cbind(sorted_sample_ids, sorted_ancestry_matrix)
sorted_ancestry_matrix_file <- file.path(lea_dir, "ancestry_matrix_sorted_with_ids.txt")
write.table(sorted_ancestry_matrix_with_ids, sorted_ancestry_matrix_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Ancestry matrix sorted by cluster proportions saved:", sorted_ancestry_matrix_file, "\n")

# Create the ancestry matrix plot (sorted)
png("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/10/ancestry_matrix_sorted.png", width = 1800, height = 800)

barchart(project, K = 10, run = best_run,
         border = NA, space = 0,
         col = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red", "green", "pink"),
         xlab = "Individuals",       
         ylab = "Proportions of ancestry",  
         main = "Ancestry Matrix (Sorted)")

# Add the sorted sample IDs to the x-axis
axis(1, at = 1:length(sorted_sample_ids), labels = sorted_sample_ids, las = 2, cex.axis = 0.5)

# Close the PNG file to save the image
dev.off()

cat("The sorted ancestry matrix plot has been saved to: /home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses_reprise/10/ancestry_matrix_sorted.png\n")

# Associate a color to each cluster
colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange", "red", "green", "pink")
cluster_colors <- colors[1:ncol(ancestry_matrix)]  # Assign a color to each cluster

# Identify the dominant cluster for each isolate
dominant_clusters <- apply(ancestry_matrix, 1, which.max)

# Associate isolates, dominant clusters, and colors
isolate_color_table <- data.frame(Sample_ID = sample_ids, 
                                  Dominant_Cluster = dominant_clusters, 
                                  Assigned_Color = cluster_colors[dominant_clusters])

# Save the table as a CSV file
isolate_color_file <- file.path(lea_dir, "isolate_colors.csv")
write.csv(isolate_color_table, isolate_color_file, row.names = FALSE)

cat("Tableau des isolats avec clusters et couleurs sauvegardé dans :", isolate_color_file, "\n")

