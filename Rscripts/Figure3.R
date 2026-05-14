#Clear environment
rm(list = ls())

# Set directories

input_dir <- "./data"
output_dir <- "./output"

# Load required libraries
library(openxlsx)
library(ggplot2)
library(ggpmisc)
library(ggalt)
library(ggsci)
library(ggpubr)
library(plyr)

# Load and preprocess data for resistant analysis
df_resistant <- read.xlsx("f_moss_bray_TPvsTP(n-1).xlsx", colNames = TRUE, sheet = 1)
df_resistant <- read.xlsx("B_mat_bray_TPvsTP(n-1).xlsx", colNames = TRUE, sheet = 1)
df_resistant$ID_N <- rep(1:8, each = 9)
df_resistant <- df_resistant[df_resistant$ID_N %in% c(1, 3, 5, 7), ]

# Plot for resistant analysis
plot_resistant <- ggplot(df_resistant, aes(x = ID_N, y = Bray_index)) +
  geom_point(size = 5, color = "#2b8cbe") +
  labs(x = "TP", y = "Bray-Curtis dissimilarity") +
  scale_x_continuous(breaks = seq(min(df_resistant$ID_N), max(df_resistant$ID_N), by = 2)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 30, face = "bold"),
        axis.title = element_text(size = 30, face = "bold"),
        panel.border = element_rect(color = "black", size = 1.5)) +
  scale_y_continuous(limits = c(0.2, 1)) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")), 
               formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top", color = "black")

plot_resistant
ggsave("../Figures_V2/f_mat_tpvstp(n-1)v2.pdf", plot_resistant, width = 6, height = 4, dpi = 300)

# Load and preprocess data for resilient analysis
df_resilient <- read.xlsx("f_moss_bray_TPvsTP.xlsx", colNames = TRUE, sheet = 1)
df_resilient <- read.xlsx("B_mat_bray_TPvsTP.xlsx", colNames = TRUE, sheet = 1)
df_resilient$ID_N <- rep(1:8, each = 9)
df_resilient <- df_resilient[df_resilient$ID_N %in% c(2, 4, 6, 8), ]

# Plot for resilient analysis
plot_resilient <- ggplot(df_resilient, aes(x = ID_N, y = Bray_index)) +
  geom_point(size = 5, color = "#2b8cbe") +
  labs(x = "TP", y = "Bray-Curtis dissimilarity") +
  scale_x_continuous(breaks = seq(min(df_resilient$ID_N), max(df_resilient$ID_N), by = 2)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 30, face = "bold"),
        axis.title = element_text(size = 30, face = "bold"),
        panel.border = element_rect(color = "black", size = 1.5)) +
  scale_y_continuous(limits = c(0.2, 1)) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")), 
               formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top", color = "black")

plot_resilient
ggsave("../Figures_V2/f_mat_tpvstpv2.pdf", plot_resilient, width = 6, height = 4, dpi = 300)
