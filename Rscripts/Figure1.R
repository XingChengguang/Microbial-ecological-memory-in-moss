# Clear environment
rm(list = ls())

# Set directories
input_dir <- "./data"
output_dir <- "./output"
setwd(output_dir)

# Load required libraries
library(openxlsx)
library(vegan)

# Read data

#bacteria
otu <- read.table("bacterial_asv.txt", header = TRUE, row.names = 1)
#fungi
otu <- read.table("fungi_asv.txt", header = TRUE, row.names = 1)
#metadata
metadata <- read.xlsx("metadata.xlsx", colNames = TRUE, rowNames = TRUE)

# Calculate relative abundance
relative_abundance <- otu / colSums(otu)

# PCoA analysis using Bray-Curtis distance
bray_dist <- vegdist(t(relative_abundance), method = "bray", binary = FALSE)
pcoa_result <- cmdscale(bray_dist, k = 3, eig = TRUE)

# Extract and merge PCoA results with metadata
pcoa_points <- as.data.frame(pcoa_result$points)
eig_percent <- round(pcoa_result$eig / sum(pcoa_result$eig) * 100, 1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
merged_data <- merge(pcoa_points, metadata, by = "row.names", all = TRUE)

# Plotting PCoA
library(ggplot2)

pcoa_plot <- ggplot(merged_data, aes(x = PCoA1, y = PCoA2, color = Type_fine)) +
  labs(x = paste0("PCoA 1 (", eig_percent[1], "%)"),
       y = paste0("PCoA 2 (", eig_percent[2], "%)")) +
  geom_point(size = 10, shape = 21) + 
  stat_ellipse(geom = "polygon", size = 0.6, level = 0.97, alpha = 0.03) +
  scale_color_manual(values = c("#469C76", "#6FB2E4", "#C6652A", "#E5BE60")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        panel.border = element_rect(color = "black", size = 1))

# Save plot
ggsave(filename = paste0(output_dir, "PCoA_plot.pdf"), plot = pcoa_plot, width = 6, height = 4, dpi = 300)

# Adonis analysis
adonis_result <- adonis2(t(relative_abundance) ~ Type_fine * Treatment_fine, 
                         data = metadata, method = "bray", permutations = 9999)

# Save Adonis result
write.xlsx(adonis_result, file = paste0(output_dir, "adonis_result.xlsx"), rowNames = TRUE)

# Anosim analysis
anosim_result <- anosim(bray_dist, metadata$Type, permutations = 999)
summary(anosim_result)
