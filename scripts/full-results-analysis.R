library(tidyverse)

presence_absence_table <- read_tsv("results/result_compound_presence_table_merged_with_metadata.tsv")

# list to store molecule information
results <- list(
  "hyodeoxycholic acid" = list(list(gnps_cluster = 340197, quality = "good", species = "2M+H")),
  "Indole-3-lactic acid" = list(list(gnps_cluster = 2627, quality = "questionable", species = "M+H")),
  "Tryptophan" = list(
    list(gnps_cluster = 2550, quality = "good", species = "M+H"),
    list(gnps_cluster = 162573, quality = "good", species = "2M+H"),
    list(gnps_cluster = 187158, quality = "good", species = "2M+Na")
  ),
  "Indole-3-propionic acid" = list(list(gnps_cluster = 1453, quality = "questionable", species = "M+H")),
  "4-hydroxyphenyllactic acid" = list(list(gnps_cluster = 371, quality = "questionable", species = "M+H-H20")),
  "Hippuric acid" = list(
    list(gnps_cluster = 659, quality = "good", species = "M+H"),
    list(gnps_cluster = 111243, quality = "good", species = "2M+H")
  ),
  "Serotonin" = list(list(gnps_cluster = 594, quality = "good", species = "M+H")),
  "D-phenyllactic acid" = list(list(gnps_cluster = 2813, quality = "good", species = "M+ACN+H")),
  "TMAO" = list(list(gnps_cluster = 131, quality = "questionable", species = "2M+H")),
  "Histamine-C8:0" = list(list(gnps_cluster = 8345, quality = "good", species = "M+H"))
)

# Convert the nested list to a dataframe
gnps_list <- list()
for (compound in names(results)) {
  entries <- results[[compound]]
  for (entry in entries) {
    gnps_list <- append(gnps_list, list(list(
      gnps_cluster_id = entry$gnps_cluster,
      compound = compound,
      quality = entry$quality,
      species = entry$species
    )))
  }
}

# Convert list to dataframe
gnps_df <- do.call(rbind, lapply(gnps_list, as.data.frame))

# Extract all cluster IDs from the gnps_df
cluster_ids <- as.character(gnps_df$gnps_cluster_id)

# Filter out samples where all molecules are absent (all 0s)
# First, identify the columns that contain cluster IDs
cluster_columns <- presence_absence_table %>% 
  select(all_of(cluster_ids))

# Filter to keep only rows where at least one molecule is present (sum > 0)
filtered_presence_absence_table <- presence_absence_table %>%
  filter(rowSums(select(., all_of(cluster_ids))) > 0)  %>% 
  filter(fermented == "yes")

# Count samples by food type
filtered_presence_absence_table %>% 
    group_by(sample_type_common) %>% 
    count()  %>% 
    arrange(desc(n))  %>% 
    print(n=61)

# Count occurrences of each molecule by food type
molecule_counts_by_food <- filtered_presence_absence_table %>%
  group_by(sample_type_common) %>%
  summarize(across(all_of(cluster_ids), sum), .groups = "drop") %>%
  # Keep only the sample_type_common column and numeric columns with sum > 0
  select(sample_type_common, all_of(cluster_ids[colSums(select(., all_of(cluster_ids))) > 0])) %>%
  arrange(desc(rowSums(across(-sample_type_common))))

# Print the counts
print(molecule_counts_by_food, n = 61)

# Update cluster_ids to only include those that remain in the dataframe
remaining_cluster_ids <- setdiff(names(molecule_counts_by_food), "sample_type_common")

# Create a more readable version with compound names instead of cluster IDs
# First, create a mapping from cluster ID to compound name
cluster_to_compound <- gnps_df %>%
  select(gnps_cluster_id, compound) %>%
  distinct()

# Function to rename columns using the mapping
rename_with_compounds <- function(df) {
  for (i in 1:nrow(cluster_to_compound)) {
    cluster_id <- as.character(cluster_to_compound$gnps_cluster_id[i])
    compound_name <- cluster_to_compound$compound[i]
    if (cluster_id %in% names(df)) {
      names(df)[names(df) == cluster_id] <- paste0(compound_name, " (", cluster_id, ")")
    }
  }
  return(df)
}

# Apply the renaming and print
molecule_counts_by_food_named <- rename_with_compounds(molecule_counts_by_food)
print(molecule_counts_by_food_named, n = 61)

# Create a heatmap of molecule presence across food types
molecule_counts_long <- molecule_counts_by_food %>%
  pivot_longer(cols = all_of(remaining_cluster_ids), 
               names_to = "cluster_id", 
               values_to = "count") %>%
  left_join(cluster_to_compound, by = c("cluster_id" = "gnps_cluster_id"))

# Plot the heatmap
ggplot(molecule_counts_long, aes(x = compound, y = sample_type_common, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Molecule Presence Across Fermented Food Types",
       x = "Compound",
       y = "Food Type",
       fill = "Count")

ggsave("molecule_heatmap.png", width = 12, height = 10)

