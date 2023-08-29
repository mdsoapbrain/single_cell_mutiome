############################
##
## This script performs the following tasks:
## 1. Imports the required libraries and loads the data.
## 2. Calls peaks based on a combination of prediction IDs and phenotypes.
## 3. Saves the peak data to a CSV file and reads it back for further analysis.
## 4. Counts the number of peaks called in each group.
## 5. Filters data based on the "ROT" and "Control" phenotypes.
## 6. Removes duplicate entries based on specific columns.
## 7. Identifies common peaks between "ROT" and "Control" groups.
## 9. Counts the number of common peaks.
##
############################

############################
# Import the tidyverse package
library(tidyverse)

# Load the data
data <- readRDS("mutiome_final.rds")
DefaultAssay(data) <- "CTpeaks"

# Create a new column combining prediction ID and phenotype
data$combo_group <- paste0(data$predict.id, "_", data$phenotype)

# Call peaks
peaks <- CallPeaks(
  object = data,
  group.by = "combo_group",
  Idents = unique(data$combo_group),
  combine.peaks = TRUE,
  macs2.path = "/opt/anaconda3/envs/env01/bin/macs2",
  outdir = "/Volumes/Seagate/multiome_cp"
)

# Convert peaks to data frame and save as CSV
peaks_df <- as.data.frame(peaks)
write.csv(peaks_df, 'CTpeaks_annotated.csv')

# Read the saved CSV
dfc <- read_csv('CTpeaks_annotated.csv')

# Count the number of peaks called in each group
table(dfc$peak_called_in)

# Separate rows and columns for further analysis
peaks_df_expanded <- dfc %>%
  separate_rows(peak_called_in, sep = ",") %>%
  separate(peak_called_in, into = c("cluster", "phenotype"), sep = "_")

# Filter ROT phenotype and select relevant columns
ROT_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "ROT") %>%
  dplyr::select(seqnames, start, end, width, cluster)

# Remove duplicates
ROT_df_unique <- ROT_df[!duplicated(ROT_df[, c("seqnames", "start", "end")]), ]

# Filter control phenotype and select relevant columns
Control_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "control") %>%
  dplyr::select(seqnames, start, end, width, cluster)

# Remove duplicates
Control_df_unique <- Control_df[!duplicated(Control_df[, c("seqnames", "start", "end")])]

# Find common peaks between ROT and Control
Common_df <- peaks_df_expanded %>%
  dplyr::filter(paste(seqnames, start, end, sep = "_") %in% paste(ROT_df$seqnames, ROT_df$start, ROT_df$end, sep = "_") &
                  paste(seqnames, start, end, sep = "_") %in% paste(Control_df$seqnames, Control_df$start, Control_df$end, sep = "_")) %>%
  dplyr::select(seqnames, start, end, width, cluster)

# Remove duplicates
Common_df_unique <- Common_df[!duplicated(Common_df[, c("seqnames", "start", "end")]), ]

# Count the number of common peaks
table(Common_df)


####### Lift over from rn7 to hg38 #####

library(rtracklayer)

# Get a list of all BED files
bed_files <- list.files('./bed_file_rn7', pattern = "*.bed", full.names = TRUE)
output_path <- "./lifted_bed_files_hg/"
# Specify the path for the chain file
chain_path <- "rn7ToHg38.over.chain"
chain_path2 <- "hg38ToHg19.over.chain"
chain_data <- rtracklayer::import.chain(chain_path)
chain_data2 <- rtracklayer::import.chain(chain_path2)

# Perform conversion for each BED file
for (bed_file in bed_files) {
  # Read the BED file
  bed_gr <- rtracklayer::import(bed_file)
  
  # Use liftOver for conversion
  lifted_bed <- liftOver(bed_gr, chain_data)
  lifted_bed2 <- liftOver(lifted_bed, chain_data2)
  # Filter out empty GRanges objects
  filtered_bed <- lifted_bed2[sapply(lifted_bed2, function(x) length(x) > 0)]
  # Convert GRangesList to a single GRanges object
  flattened_bed <- unlist(filtered_bed)
  # Save the result as a new BED file
  output_file <- paste0(output_path, basename(bed_file))
  export(flattened_bed, output_file, format = "bed")
  
  # Print conversion statistics for each file
  cat(paste(basename(bed_file), ":", length(flattened_bed), "regions lifted\n"))
}

