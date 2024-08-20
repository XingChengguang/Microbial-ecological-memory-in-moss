
#SEM
# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(devtools)
library(plspm)
library(openxlsx)
library(tibble)
library(tidyr)
library(dplyr)
library(gtools)
library(readr)

# Set working directory
setwd("input/SEM/")

# Load and merge all files
files <- list.files()
df1 <- read.xlsx("../Treatment.xlsx")

for (i in files) {
  df <- read.xlsx(i)
  df1 <- merge(df1, df, by = "Sample")
}

sem <- df1[, -1]  # Remove the first column (Sample IDs)

# Define the path model
Treatment = c(0,0,0,0,0,0,0,0,0)
Metabolite = c(1,0,0,0,0,0,0,0,0)
PPB = c(1,1,0,0,0,0,0,0,0)
Fungi = c(1,1,0,0,0,0,0,0,0)
C = c(0,0,1,1,0,0,0,0,0)
M = c(0,0,1,1,0,0,0,0,0)
N = c(0,0,1,1,0,0,0,0,0)
P = c(0,0,1,1,0,0,0,0,0)
S = c(0,0,1,1,0,0,0,0,0)

foot_path <- rbind(Treatment, Metabolite, PPB, Fungi, C, M, N, P, S)
colnames(foot_path) <- rownames(foot_path)

# Visualize the model structure
innerplot(foot_path)

# Define the blocks corresponding to the model
foot_blocks <- list(1, 3, 6:15, 16:25, 26:34, 35:48, 49:58, 59:72, 73:78)

# Convert relevant columns to factors
sem$Treatment <- as.factor(sem$Treatment)
sem$Frequncy <- as.factor(sem$Frequncy)

# Perform PLS-PM
df_pls <- plspm(sem, foot_path, blocks = foot_blocks)

# Visualize loadings
plot(df_pls, what = "loadings", arr.width = 0.1, show.values = TRUE, lcol = 'gray')

# Visualize the path model
innerplot(df_pls, colpos = '#CA5023', colneg = '#457CC3', show.values = TRUE, lcol = 'gray20', box.lwd = 0)

# Summarize the PLS-PM results
summary(df_pls, rsquare = TRUE, standardized = TRUE, fit.measures = TRUE)

# Save the inner model results to a CSV file
a <- do.call(rbind.data.frame, df_pls$inner_model) %>%
  tibble::rownames_to_column() %>%
  tidyr::separate(col = rowname, into = c("Response", "Predictor"), sep = "[.]") %>%
  dplyr::filter(Predictor != "Intercept") %>%
  dplyr::mutate_if(is.numeric, ~ round(.x, 4)) %>%
  dplyr::mutate(" " = gtools::stars.pval(`Pr(>|t|)`))

print(a)

do.call(rbind.data.frame, df_pls$inner_model) %>%
  tibble::rownames_to_column() %>%
  tidyr::separate(col = rowname, into = c("Response", "Predictor"), sep = "[.]") %>%
  dplyr::filter(Predictor != "Intercept") %>%
  dplyr::relocate(from = Predictor, to = Response, weight = Estimate, p = `Pr(>|t|)`) %>%
  readr::write_csv("sem_table.csv")
