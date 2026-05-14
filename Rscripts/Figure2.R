

#PCoA
# Clear environment
rm(list = ls())
# Set directories

input_dir <- "./data"
output_dir <- "./output"

# Load required libraries
library(openxlsx)
library(vegan)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

# Load data
#bacteria
otu <- read.table("bacterial_asv.txt", header = TRUE, row.names = 1)
#fungi
otu <- read.table("fungi_asv.txt", header = TRUE, row.names = 1)
#metadata
metadata <- read.xlsx("metadata.xlsx", colNames = TRUE, rowNames = TRUE)

# Calculate relative abundance
relative_abundance <- count_seq / colSums(count_seq)

# Filter metadata for specific conditions
metadata <- metadata[metadata$Type == "Moss" & metadata$Treatment_fine == "DT_S", ]
count_seq <- relative_abundance[, row.names(metadata)]
count_seq <- count_seq[rowSums(count_seq) > 0, ]

# PCoA analysis
bray_dist <- vegdist(t(count_seq), method = "bray", binary = FALSE)
pcoa_result <- cmdscale(bray_dist, k = 3, eig = TRUE)

# Extract and merge PCoA results with metadata
pcoa_points <- as.data.frame(pcoa_result$points)
eig_percent <- round(pcoa_result$eig / sum(pcoa_result$eig) * 100, 1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
merged_data <- merge(pcoa_points, metadata, by = "row.names", all = TRUE)

# Calculate centroids
centroid <- aggregate(cbind(PCoA1, PCoA2) ~ Group, data = merged_data, FUN = mean)
pcoa_results <- left_join(merged_data, centroid, by = "Group", suffix = c("", ".cen"))

# Plotting PCoA
colors <- brewer.pal(8, "Set2")

pcoa_plot <- ggplot(pcoa_results, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Group)) +
  geom_segment(aes(xend = PCoA1.cen, yend = PCoA2.cen, color = Group), show.legend = FALSE) +
  geom_label(data = centroid, aes(label = Group, fill = Group), size = 5, color = "white", show.legend = FALSE) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        panel.border = element_rect(color = "black", size = 1)) +
  labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
       y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""))

pcoa_plot

#ggsave(paste0(otudir, ".pdf"), pcoa_plot, width = 6, height = 4, dpi = 300)

# Adonis analysis
adonis_result <- adonis2(t(count_seq) ~ Group, data = metadata, method = "bray", permutations = 9999)
print(adonis_result)
# ANOSIM analysis
anosim_result <- anosim(bray_dist, metadata$Group, permutations = 999, distance = "bray")
print(anosim_result)

# MRPP analysis
mrpp_result <- mrpp(bray_dist, metadata$Group, permutations = 999, distance = "bray")
print(mrpp_result)


#Assembly
# Clear environment
rm(list = ls())

# Set directories

input_dir <- "./data"
output_dir <- "./output"

# Load required libraries
library(openxlsx)
library(ggpubr)
library(ggplot2)

# Load betaMNTD data
betaMNTD <- read.csv("./input_dir/S_MNTD.csv")
betaMNTD <- read.csv("./input_dir/Moss_MNTD.csv")
betaMNTD <- read.csv("./input_dir/Mat_MNTD.csv")

# Load sample group data
group <- read.xlsx("./metadata.xlsx", colNames = TRUE, rowNames = TRUE)

# Filter by sample type
group <- subset(group, Type %in% c("Moss", "Mat"))

# Extract sample pairs by treatment group
extract_sample_pairs <- function(treatment_group, data, group_data) {
  sample_names <- rownames(subset(group_data, Treatment_fine == treatment_group))
  subset(data, name1 %in% sample_names & name2 %in% sample_names)
}

betaMNTD_Control <- extract_sample_pairs("CT", betaMNTD, group)
betaMNTD_Treat_S <- extract_sample_pairs("DT_S", betaMNTD, group)
betaMNTD_Treat_R <- extract_sample_pairs("DT_R", betaMNTD, group)

# Calculate contribution rates
calculate_contribution_rate <- function(data, threshold) {
  nrow(data[abs(data$bNTI.wt) < threshold, ]) / nrow(data)
}

# Random factor contribution rate
random_control <- calculate_contribution_rate(betaMNTD_Control, 2)
random_treat_s <- calculate_contribution_rate(betaMNTD_Treat_S, 2)
random_treat_r <- calculate_contribution_rate(betaMNTD_Treat_R, 2)

# Deterministic factor contribution rate
deterministic_control <- 1 - random_control
deterministic_treat_s <- 1 - random_treat_s
deterministic_treat_r <- 1 - random_treat_r

# Combine data for plotting
betaMNTD_Control$group <- 'CT'
betaMNTD_Treat_S$group <- 'DT_S'
betaMNTD_Treat_R$group <- 'DT_R'
betaMNTD_group <- rbind(betaMNTD_Control, betaMNTD_Treat_S, betaMNTD_Treat_R)

# Plotting
plot_betaNTI <- function(data, output_file) {
  p <- ggplot(data, aes(group, bNTI.wt, color = group)) +
    stat_boxplot(geom = "errorbar", width = 0.2, size = 1) +
    geom_boxplot(width = 0.2, size = 1) +
    geom_violin(aes(fill = group, color = group), trim = TRUE, width = 1, alpha = 0.1) +
    geom_hline(aes(yintercept = 2), size = 0.6, linetype = "dashed", colour = "black") +
    geom_hline(aes(yintercept = -2), size = 0.6, linetype = "dashed", colour = "black") +
    labs(y = 'betaNTI') +
    theme_classic() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    scale_color_manual(values = c("#468B00", '#E18727', "#BD3C29")) +
    scale_y_continuous(limits = c(-7, 10))
  
  ggsave(paste0(otudir, output_file), p, width = 6, height = 4, dpi = 300)
}
