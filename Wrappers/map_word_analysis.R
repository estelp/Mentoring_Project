# Load necessary libraries
library(ggplot2)
library(ggforce)
library(sf)
library(dplyr)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyr)  # For the pivot_longer function
library(readxl)

# Define file paths
coords_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Topic/CIBIG_Coordonnates.xlsx"
ancestry_matrix <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses/ancestry_matrix_with_ids.txt"
output_map <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/RLEA_analyses/carte/ancestry_map.png"

# Read the Excel file containing coordinates
coords_data <- read_excel(coords_file, sheet = 1)  # sheet = 1 refers to the first sheet

# Preview the data
print(head(coords_data))

# Rename columns for longitude and latitude
coords_data_rename <- coords_data %>%
  rename(lon = Longitude, lat = Latitude) %>%
  select(lon, lat)  # Select only the necessary columns

# Ensure that coordinates are numeric (no text or other formats)
coords_data_rename$lon <- as.numeric(coords_data_rename$lon)
coords_data_rename$lat <- as.numeric(coords_data_rename$lat)

# Load the ancestry matrix
cat("Loading the ancestry matrix...\n")
ancestry_data <- read.table(ancestry_matrix, header = TRUE, stringsAsFactors = FALSE)

# Check the original column names
cat("Original column names in the ancestry matrix:\n")
print(colnames(ancestry_data))

# Replace the first column (sample IDs) with numeric values from 1 to 89
ancestry_data[, 1] <- as.character(1:nrow(ancestry_data))

# Verify the new column names
cat("Column names after modification:\n")
print(colnames(ancestry_data))

# Rename the columns of the ancestry matrix starting from the second column
colnames(ancestry_data)[2:ncol(ancestry_data)] <- paste0("Cluster", 1:(ncol(ancestry_data) - 1))

# Check the new column names
cat("Column names after renaming:\n")
print(colnames(ancestry_data))

# Combine the coordinates with the ancestry matrix
cat("Combining coordinates with the ancestry matrix...\n")
map_data <- cbind(coords_data_rename, ancestry_data)

# Verify the new column names
cat("Column names after combining:\n")
print(colnames(map_data))

# Reshape the data to have a "Cluster" column and a "Proportion" column
map_data_long <- map_data %>%
  pivot_longer(
    cols = starts_with("Cluster"),  # Select columns starting with "Cluster"
    names_to = "Cluster",           # New column for cluster names
    values_to = "Proportion"        # New column for proportions
  ) %>%
  mutate(
    Cluster = factor(Cluster, levels = paste0("Cluster", 1:7))  # Ensure a consistent order of clusters
  )

# Preview the reshaped data
cat("Preview of reshaped data:\n")
print(head(map_data_long))

# Load geographic data for the world map
cat("Loading geographic data...\n")
world <- st_as_sf(ne_countries(scale = "medium", returnclass = "sf"))

# Prepare data for circular charts
cat("Preparing data for circular charts...\n")
map_data_long <- map_data_long %>%
  group_by(sample_ids) %>%
  mutate(
    start_angle = cumsum(Proportion) * 2 * pi - Proportion * 2 * pi,  # Start angle
    end_angle = cumsum(Proportion) * 2 * pi                          # End angle
  )

# Ensure that 'lon' and 'lat' columns are numeric
map_data_long$lon <- as.numeric(map_data_long$lon)
map_data_long$lat <- as.numeric(map_data_long$lat)

# Check the data types after conversion
str(map_data_long)

# Create the map with circular diagrams
cat("Creating the map with circular diagrams...\n")

cluster_colors <- c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange")  # Colors for 7 clusters

map_plot <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray50") +  # Map background
  geom_arc_bar(
    data = map_data_long,
    aes(
      x0 = lon, y0 = lat,  # Central coordinates for the circle
      r0 = 0,              # Inner radius (full circle)
      r = 5,               # Outer radius (adjust based on point density)
      start = start_angle, # Start angle
      end = end_angle,     # End angle
      fill = Cluster       # Color based on the cluster
    ),
    alpha = 1  # Full opacity
  ) +
  scale_fill_manual(
    values = cluster_colors,  # Colors defined for the clusters
    name = "Cluster"          # Name for the legend
  ) +
  coord_sf(crs = 4326) +  # Project the map in "longlat" coordinate system
  
  # Customize axes
  scale_x_continuous(
    breaks = seq(-180, 180, by = 30),  # Set longitude tick marks
    labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))  # Add direction symbols
  ) +
  scale_y_continuous(
    breaks = seq(-90, 90, by = 30),  # Set latitude tick marks
    labels = function(x) paste0(abs(x), ifelse(x < 0, "°S", "°N"))  # Add direction symbols
  ) +
  labs(
    title = "Distribution of Isolates with Ancestry Proportions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +  # Minimalistic theme
  theme(
    legend.position = "right", # Legend position
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center the title
  )

# Display the map
print(map_plot)

# Save the map as a PNG image
cat("Saving the map as a PNG image...\n")
ggsave(output_map, plot = map_plot, width = 10, height = 7, dpi = 300)

cat("The map with ancestry proportions has been saved to:", output_map, "\n")

