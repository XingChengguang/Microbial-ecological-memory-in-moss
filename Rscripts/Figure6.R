
#Ternary plot
# Clear environment and set working directory
rm(list = ls())
setwd("input/Metagenomic/")

# Load necessary libraries
library(openxlsx)
library(magrittr)

# List of analyses to perform
analysis_types <- c("CAzy", "GO", "KO", "pfam")

# Loop through each analysis type
for (i in analysis_types) {
  
  # Load the relative abundance data
  asv <- read.table(paste0('../../../CoverM_results/Bacteria_', i, '_rpm.txt'), header = TRUE)
  meta <- read.xlsx("../../../metadata.xlsx", colNames = TRUE, rowNames = TRUE)
  
  # Subset the data to include only relevant metadata
  asv <- asv[, row.names(meta)]
  
  # Transpose and merge with metadata
  data <- as.data.frame(t(asv))
  all <- merge(meta, data, by = "row.names")
  row.names(all) <- all$Row.names
  
  # Remove unnecessary columns and calculate mean values for each group (CT, DTR, DTS)
  df <- all[, -c(1:3)]  # Keep only relevant columns
  grouped_df <- aggregate(. ~ all[, 4], data = df, FUN = mean)
  row.names(grouped_df) <- grouped_df$`all[, 4]`
  grouped_df <- grouped_df[, -1]
  grouped_df <- as.data.frame(t(grouped_df))
  grouped_df$meanScore <- rowMeans(grouped_df[, 1:3])
  
  # Combine with enrichment data
  enrich1 <- read.xlsx(paste0("Bacteria_Bacteria_", i, "_edge_dt1.xlsx"), colNames = TRUE, rowNames = TRUE)
  enrich7 <- read.xlsx(paste0("Bacteria_Bacteria_", i, "_edge_dt7.xlsx"), colNames = TRUE, rowNames = TRUE)
  
  enrich <- merge(enrich1, enrich7, by = "row.names", all = TRUE)
  enrich <- enrich[, c(1, 7, 13)]
  row.names(enrich) <- enrich$Row.names
  
  # Merge with grouped data
  all <- merge(enrich, grouped_df, by = "row.names")
  
  # Create enrichment categories
  all$a <- paste0(all$col.x, all$col.y)
  all$b <- ifelse(all$a == "EnrichedEnriched", "com",
                  ifelse(all$a == "EnrichedNon_sig", "DT1",
                         ifelse(all$a == "Non_sigEnriched", "DT7", "nosig")))
  
  all$c <- ifelse(all$a == "DepletedDepleted", "zcom",
                  ifelse(all$a == "DepletedNon_sig", "zDT1",
                         ifelse(all$a == "Non_sigDepleted", "zDT7", "nosig")))
  
  # Remove unnecessary columns and save the results
  all <- all[, -2]
  write.xlsx(all, paste0("./Ternary_plot/", i, "_tern_moss.xlsx"))
}

#heatmap

# Clear the workspace
rm(list = ls())

# Set working directory
setwd("./input/Metagenomic/C_N_P_S/Raw_reads/edgeRdtvsct/")

# Load required libraries
library(openxlsx)
library(pheatmap)

# Load data for DT1 and DT7
dt1 <- read.xlsx("Ccycle_edge_dt1.xlsx")  # Adjust filename accordingly
colnames(dt1)[1] <- "Gene"

dt7 <- read.xlsx("Ccycle_edge_dt7.xlsx")  # Adjust filename accordingly
colnames(dt7)[1] <- "Gene"

# Filter enriched genes
enrich1 <- dt1[dt1$col == "Enriched", ]
enrich7 <- dt7[dt7$col == "Enriched", ]

# Find common enriched genes
com <- intersect(enrich1$Gene, enrich7$Gene)

# Load raw expression data
data <- read.table("../Ccycle.raw.txt", sep = ' ', header = TRUE)  # Adjust file path accordingly
# data <- read.table("../Scycle.raw.txt", sep = '\t', header = TRUE)  # Uncomment and adjust for S cycle data

row.names(data) <- data$Gene
data <- data[, -1]

# Load metadata
group <- read.xlsx('../../../metadata.xlsx', colNames = TRUE, rowNames = TRUE)

# Ensure data columns are ordered by the metadata
data <- data[, row.names(group)]

# Filter genes based on enrichment
filter <- unique(c(enrich1$Gene, enrich7$Gene))
df <- data[filter, ]

# Load pathway annotation and filter it by the enriched genes
pathway <- read.xlsx("../../Cycs_Pathway_annatation.xlsx", sheet = 5)  # Adjust sheet number accordingly
path_df <- pathway[pathway$Gene %in% filter, ]

row.names(path_df) <- path_df$Gene
path_df <- path_df[order(path_df$Pathway), ]

# Reorder df to match the ordered pathway data
df <- df[path_df$Gene, ]

# Z-score transformation for heatmap
data4 <- t(scale(t(df)))

# Row annotation based on pathways
annotation_row <- path_df[, c("Gene", "Pathway")]

# Column annotation based on group metadata
annotation_col <- group[, c("Treatment", "Treatment_fine")]

# Generate heatmap
p <- pheatmap(data4,
              display_numbers = FALSE,
              cellwidth = 20, cellheight = 4,  # Adjust cell size
              cluster_row = FALSE, cluster_cols = FALSE,
              labels_col = group$Treatment_fine,
              show_colnames = TRUE,
              fontsize_row = 5, fontsize_col = 5,
              border_color = NA,  # Remove cell borders
              scale = "none",
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              color = colorRampPalette(c("#F4F2ED", "#DC1432"))(5))

# Save the heatmap
ggsave("C_enrich_heatmap.pdf", p, width = 6, height = 4, dpi = 300)

