

# Remove all objects from the workspace
rm(list = ls())

# Load the data from an RDS file
data <- readRDS("mutiome_final.rds")

# Display the frequency table for phenotype and predict.id
table(data$phenotype)
table(data$predict.id)

library(dplyr)
library(ggplot2)

# Extract metadata from the Seurat object
metadata <- data@meta.data

# Calculate the count of each predict.id in each phenotype
count_data <- metadata %>% 
  group_by(phenotype, predict.id) %>% 
  summarise(count = n())

# Calculate the proportion of each predict.id in each phenotype
count_data <- count_data %>% 
  group_by(phenotype) %>% 
  mutate(proportion = count / sum(count))

# Plot a stacked bar chart
ggplot(count_data, aes(x = phenotype, y = proportion, fill = predict.id)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion", x = "Phenotype", fill = "Celltypes") +
  theme_minimal()

# Save the count data to a CSV file
write.csv(count_data, 'count_data.csv')
