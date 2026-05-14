
#PBB:ct:dtr:dts
# Clear environment
rm(list = ls())

# Load required libraries
library(openxlsx)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(ggpmisc)
library(rfPermute)
library(ggridges)
library(ggpubr)
library(A3)
library(pheatmap)

# Set working directory
setwd("/Users/xingchengguang/Desktop/phd_data2/16s_V4/PPB/Version2/")

# Load and preprocess PPB data
ppb <- read.xlsx("ppb_moss_sum_new.xlsx", colNames = TRUE, rowNames = TRUE)
ppb <- ppb[, 1:48] %>% t() %>% as.data.frame()  # Transpose and convert to data frame

# Load metadata
meta <- read.xlsx("../metadata.xlsx", colNames = TRUE, rowNames = TRUE)
compart <- meta[meta$Type == "Moss", ]

# Combine data for plotting
df_plot <- cbind(compart, ppb)
df_plot <- df_plot[, c(4, 7:17)]
seq <- colnames(ppb[, 1:17])

# Convert data to long format
df_long <- df_plot %>%
  pivot_longer(cols = -c(Treatment_fine), names_to = "PPB", values_to = "Value")

# Set factor levels for PPB
df_long$PPB <- factor(df_long$PPB, levels = seq)

# Plotting
p2 <- ggplot(data = df_long, aes(PPB, Value, color = Treatment_fine)) +
  stat_boxplot(geom = "errorbar", width = 0.6, size = 1, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.6, size = 1, position = position_dodge(0.8), outlier.shape = NA) +
  geom_signif(comparisons = list(c("CT", "DT_R"), c("CT", "DT_S"), c("DT_R", "DT_S")),
              textsize = 2, test = t.test, step_increase = 0.1, map_signif_level = TRUE,
              margin_top = 0.05, position = position_dodge(0.8), size = 0.4) +
  geom_jitter(aes(PPB, Value), size = 3, alpha = 0.8, position = position_dodge(0.8)) +
  labs(y = "") +
  theme_classic() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"),
        axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("#00665C", "#0074b3", "#982b2b", "#512E03"))

p2


#
# Load necessary libraries
library(openxlsx)
library(magrittr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(ggpmisc)
library(ggpubr)
library(rstatix)
library(ggridges)
library(pheatmap)
library(rfPermute)
library(A3)

# Set working directory and output paths
set_working_dirs <- function() {
  output <- "/Users/xingchengguang/Desktop/phd_data2/Figures_V2/"
  setwd("/Users/xingchengguang/Desktop/phd_data2/16s_V4/Composition/")
  return(output)
}

output <- set_working_dirs()

# Function to load and preprocess data
load_data <- function() {
  genus <- read.xlsx("genus_relative.xlsx", colNames = T, rowNames = F)
  genus <- genus[,1:144] # Remove TP0
  ppb <- read.xlsx("/Users/xingchengguang/Desktop/phd_data2/16s_V4/PPB/PBB_Database.xlsx", colNames = T, rowNames = T)
  meta <- read.xlsx("../metadata.xlsx", colNames = T, rowNames = T)
  return(list(genus = genus, ppb = ppb, meta = meta))
}

data <- load_data()
genus <- data$genus
ppb <- data$ppb
meta <- data$meta

# Filter PPB data
filter_ppb <- function(genus, ppb) {
  my_data_ppb <- genus[grepl(paste(row.names(ppb), collapse = "|"), genus$Genus), "Genus"] %>% as.data.frame()
  my_data_ppb <- genus[genus$Genus %in% my_data_ppb$.,]
  row.names(my_data_ppb) <- my_data_ppb$Genus
  ppb <- my_data_ppb[,-1]
  return(ppb)
}

ppb <- filter_ppb(genus, ppb)

# Plotting Function
plot_boxplot <- function(df_plot, comparison, x_col, y_col, color_col, output_path, facet_col = NULL, facet_scale = "free") {
  p <- ggplot(data = df_plot, aes_string(x = x_col, y = y_col, color = color_col)) +
    stat_boxplot(geom = "errorbar", width = 0.3, size = 0.6) +
    geom_boxplot(width = 0.3, size = 0.6, geom = "errorbar") +
    geom_jitter(size = 3, shape = 1) +
    geom_signif(comparisons = comparison,
                textsize = 3,
                test = "wilcox.test",
                step_increase = 0.1,
                map_signif_level = T,
                margin_top = 0.05,
                position = "identity",
                inherit.aes = T,
                size = 0.4) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    scale_color_manual(values = c("#00665C","#BD7F2C","#8A4F09","#512E03" ))
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_col)), scales = facet_scale)
  }
  
  ggsave(output_path, p, width = 6, height = 4, dpi = 300)
}

# Prepare data for plotting
prepare_plot_data <- function(ppb, meta) {
  csum_ppb <- colSums(ppb) %>% as.data.frame()
  colnames(csum_ppb) <- "Sum"
  df_plot <- merge(meta, csum_ppb, by = "row.names")
  df_plot <- df_plot[df_plot$Type_fine != "OS",]
  df_plot$Type_fine <- factor(df_plot$Type_fine, levels = c("Moss", "Mat", "S"))
  return(df_plot)
}

df_plot <- prepare_plot_data(ppb, meta)

# Boxplot: Type Comparison
plot_boxplot(df_plot, list(c("Moss", "Mat"), c("Moss", "S"), c("Mat", "S")), 
             "Type_fine", "Sum", "Type_fine", paste0(output, "PPB_Type.pdf"))

# Boxplot: Treatment Comparison
plot_boxplot(df_plot, list(c("CT", "DT")), 
             "Treatment", "Sum", "Type_fine", paste0(output, "PPB_Treatment.pdf"), "Type_fine")

# Example of how to plot and save specific facets (e.g., for "Mat" or "Soil")
# You can add similar blocks for different facets or treatments as needed
moss <- meta[meta$Type == "Mat",]
df_plot_moss <- df_plot[rownames(moss),]
plot_boxplot(df_plot_moss, list(c("CT","DT_R"), c("CT","DT_S"), c("DT_R","DT_S")), 
             "Treatment_fine", "Sum", "Treatment_fine", paste0(output, "mat_Treatment.pdf"))

# Repeat similar steps for other specific plots or analyses as needed...

##TOP20 PBB

# Clear environment and load necessary libraries
rm(list = ls())
library("openxlsx")
library(magrittr)
library(reshape2)
library(ggplot2)
library(ggalluvial)

# Set working directory and output path
output <- "./output"
setwd("./input/")

# Load data
otu_relative <- read.xlsx("genus_relative.xlsx", colNames = TRUE, rowNames = TRUE)
otu_relative <- otu_relative[, 1:144]  # Remove TP0 columns

metadata <- read.xlsx("../metadata.xlsx", colNames = TRUE, rowNames = TRUE)
metadata <- metadata[!row.names(metadata) %in% c("TP0_1", "TP0_2", "TP0_3"),]

# Filter for Moss type
moss_metadata <- metadata[metadata$Type == "Moss", ]  # Filter for Moss samples
moss_otu_relative <- otu_relative[, row.names(moss_metadata)]

# Transpose and sort by mean relative abundance
data_pctg <- t(moss_otu_relative)
dat <- data_pctg[, order(colMeans(data_pctg), decreasing = TRUE)]

# Keep the top 25 genera
dat_top25 <- dat[, 1:25]
colnames(dat_top25) <- gsub("g__", "", colnames(dat_top25))  # Remove "g__" prefix
colnames(dat_top25)[5] <- "f__Beijerinckiaceae"  # Rename specific genus

# Calculate means across replicates for each time point
rows <- unique(sub("_\\d", "", row.names(dat_top25)))
mean_result <- apply(dat_top25, 2, function(x) {
  sapply(rows, function(y) {
    mean(x[grep(y, row.names(dat_top25))])
  })
})

# Merge with metadata
meta <- read.xlsx("meta_for_tp_line.xlsx", colNames = TRUE, rowNames = TRUE)
meta$TP <- paste0("TP", meta$TP)
meta_types <- meta[meta$Type == "Moss", ]
dat2 <- cbind(mean_result, meta_types)

# Convert to long format for plotting
dat_m <- melt(dat2)

# Plot the top 25 genera as an alluvial plot
p1 <- ggplot(data = dat_m, aes(x = TP, y = value, fill = variable, stratum = variable, alluvium = variable)) +
  geom_stratum(color = "black", width = 0.6, size = 0) +
  geom_flow(alpha = 0.5) +
  scale_fill_manual(values = c("#db6968","#4d97cd","#99cbeb","#459943",
                               "#fdc58f","#e8c559","#a3d393","#f8984e","#FFC0CB", "#FFA500", "#FFFF00", "#0000FF", "#800080",
                               "#FFD700", "#008000", "#FF0000",  "#808080",  "#800000", 
                               "#F0E68C", "#00CED1", "#FF1493", "#1E90FF", "#8B0000", "#FFDAB9", "#FF4500", "#8A2BE2",
                               "#008080", "#FA8072", "#4682B4", "#D2B48C", "#ADFF2F", "#4B0082", "#DC143C", "#2E8B57",
                               "#A0522D", "#6A5ACD", "#7FFF00")) +
  theme(axis.title.y = element_text(face = 'bold', color = 'black', size = 14),
        axis.title.x = element_text(face = 'bold', color = 'black', size = 10, vjust = -1.2),
        axis.text.y = element_text(face = 'bold', color = 'black', size = 10),
        axis.text.x = element_text(face = 'bold', color = 'black', size = 12, angle = 45, vjust = 0.5),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.key.height = unit(0.1, 'cm'),
        legend.text = element_text(face = 'bold', color = 'black', size = 5)) +
  labs(x = '', y = 'Relative Abundance', fill = NULL) +
  facet_grid(~Treatment_fine, drop = TRUE, scale = "free", space = "free_x")

# Save the plot
ggsave(paste0(output, "moss_top25_genues.pdf"), p1, width = 6, height = 4, dpi = 300)

