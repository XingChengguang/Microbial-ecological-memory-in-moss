# Clear environment
rm(list = ls())

# Load required libraries
library(WGCNA)
library(openxlsx)
library(graphlayouts)
library(tidyverse)
library(igraph)
library(ggraph)
library(oaqc)
library(ggforce)
library(concaveman)
library(psych)
library(reshape2)
library(ggsci)
library(ggsignif)
library(ggplot2)

# Set directories
dt_dir <- "./cross_domain_network/dt_fine/"
ct_dir <- "./cross_domain_network/ct_fine/"
setwd(ct_dir)

# Function to process data files
process_files <- function(directory, time_points) {
  df_list <- list()
  
  for (time_point in time_points) { 
    cor_df <- read.xlsx(paste0(time_point, "moss_net.xlsx"))
    
    # Identify bacteria or fungi
    cor_df$BF <- paste0(gsub("_.*", "", cor_df$Source), gsub("_.*", "", cor_df$Target))
    
    # Identify positive or negative correlation
    cor_df$PN <- ifelse(cor_df$cor > 0, "P", "N")
    
    # Combine correlation type and organism type
    cor_df$all <- paste0(cor_df$PN, "_", cor_df$BF)
    
    # Calculate frequencies and percentages
    freq <- table(cor_df$all)
    percentage <- prop.table(freq) * 100
    percentage_df <- as.data.frame(percentage)
    
    # Add time point to the data frame
    percentage_df$Time <- time_point
    df_list[[time_point]] <- percentage_df
  }
  
  return(do.call(rbind, df_list))
}

# Process CT and DT data
df_ct <- process_files(ct_dir, paste("TP", 1:8, sep = ""))
df_ct$Treat <- "CT"

df_dt <- process_files(dt_dir, paste("TP", 1:8, sep = ""))
df_dt$Treat <- rep(c("DTS", "DTR"), each = 6)

# Combine CT and DT data
df_all <- rbind(df_ct, df_dt)

# Generate comparison combinations
combinations <- combn(unique(df_all$Treat), 2, simplify = FALSE)

# Plot frequency distributions with significance testing
plot_frequency <- function(data, output_file) {
  p <- ggplot(data, aes(Treat, Freq, color = Treat)) +
    stat_boxplot(geom = "errorbar", width = 0.3, size = 0.3) +
    geom_boxplot(width = 0.3, size = 0.3) +
    geom_signif(comparisons = combinations,
                textsize = 2,
                test = "wilcox.test",
                step_increase = 0.1,
                map_signif_level = TRUE,
                margin_top = 0.05,
                position = "identity",
                inherit.aes = TRUE,
                size = 0.4) +
    geom_point(size = 1) +
    labs(x = "", y = "") +
    theme_classic() +
    scale_color_manual(values = c("#FC7F00", "#984EA3", "#1e90ff", "#1e90ff")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Var1, scales = "free")
  
  ggsave(output_file, p, width = 6, height = 4, dpi = 300)
}

plot_frequency(df_all, "../dt_fine/dts_P_N_B_F_sig.pdf")

# Process original data for CT and DT
df_ct_raw <- process_files(ct_dir, paste("TP", 1:8, sep = ""))
df_ct_raw$Treat <- "CT"

df_dt_raw <- process_files(dt_dir, paste("TP", 1:8, sep = ""))
df_dt_raw$Treat <- rep(c("DTS", "DTR"), each = 6)

df_all_raw <- rbind(df_ct_raw, df_dt_raw)

# Plot raw frequency distributions with significance testing
plot_frequency(df_all_raw, "../dt_fine/dts_P_N_B_F_sig_raw.pdf")

# Line plot for time series data
df_all_raw$Time <- gsub("TP", "", df_all_raw$Time)
df_all_raw$BF <- gsub("*._", "", df_all_raw$Var1)

line_plot <- function(data, output_file) {
  mypal <- pal_npg("nrc", alpha = 1)(6)
  
  p <- ggplot(data, aes(x = Time, y = Freq, color = Var1, shape = BF)) +
    geom_line(aes(as.numeric(Time), Freq), position = position_dodge(0.6), cex = 1) +
    scale_x_continuous(breaks = seq(0, 2, 7)) +
    facet_wrap(~ Treat, scales = "free") +
    theme_bw() +
    scale_color_manual(values = mypal)
  
  ggsave(output_file, p, width = 6, height = 4, dpi = 300)
}

line_plot(df_all_raw, "../dt_fine/dts_P_N_B_F_line.pdf")

# Box plots for different correlation types
box_plot_correlation <- function(data, treatment, output_file) {
  data <- subset(data, Treat == treatment)
  data$Var1 <- factor(data$Var1, levels = c("P_BB", "P_FB", "P_FF", "N_BB", "N_FB", "N_FF"))
  
  p <- ggplot(data, aes(Var1, Freq, color = Var1)) +
    stat_boxplot(geom = "errorbar", width = 0.3, size = 0.3) +
    geom_boxplot(width = 0.3, size = 0.3) +
    geom_point(size = 1) +
    labs(x = "", y = "") +
    theme_classic() +
    scale_color_manual(values = c("#FC7F00", "#FC7F00", "#FC7F00", "#984EA3", "#984EA3", "#984EA3")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, p, width = 6, height = 4, dpi = 300)
}

box_plot_correlation(df_ct_raw, "CT", "../ct_fine/P_N_B_F.pdf")
box_plot_correlation(df_dt_raw, "DTR", "../dt_fine/dtr_P_N_B_F.pdf")
box_plot_correlation(df_dt_raw, "DTS", "../dt_fine/dts_P_N_B_F.pdf")
