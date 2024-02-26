################################################################################
########################## Notes on provided code ##############################
################################################################################

# The code provide is intended for processing of Raman spectra. The code is broken
# into sections for each processing step to allow for trouble shooting.

################################################################################
############# Load libraries, set working directory, read file #################
################################################################################

library(tidyverse)
library(alkahest)
library(pracma)
library(ptw)
library(plotly)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load in the file saved from the txt-to-dataframe code
Raman_spectra <- read.csv("Raman_spectra.csv")

################################################################################
############# Interactive plot to identify outliers in data set ################
################################################################################

# Hover curser over lines that appear to be outliers

unique_groups <- unique(Raman_spectra$FileID)

p <- plot_ly(type = "scatter", mode = "lines", source = "source")

for (i in seq_along(unique_groups)) {
  group_value <- unique_groups[i]
  group_data <- subset(Raman_spectra, FileID == group_value)
  
  p <- add_trace(p, x = group_data$V1, y = group_data$V2, line = list(color = "gray", width = 2), text = group_value)
}

p <- layout(p, xaxis = list(title = "Wavenumber"), yaxis = list(title = "Intensity"), title = "Original unprocessed data")

p

################################################################################
######################### Savitzky-Golay smoothing  ############################
################################################################################

unique_groups <- unique(Raman_spectra$FileID)

smoothed_Raman_spectra <- data.frame()

for (group_value in unique_groups) {
  group_data <- subset(Raman_spectra, FileID == group_value)

  # m = window, p = degree of polynomial. Adjust as desired.
  smoothed_values <- smooth_savitzky(group_data$V1, group_data$V2, m = 11, p = 2)
  
  smoothed_Raman_spectra <- rbind(smoothed_Raman_spectra, 
                                data.frame(V1 = group_data$V1, 
                                           V2_smoothed = smoothed_values$y, 
                                           FileID = group_value))
}

head(smoothed_Raman_spectra)

################################################################################
########################## Baselining of spectra ###############################
################################################################################

# Only one baselineing method should be applied to spectra. Four are provided
# for user to select from:


##################### Asymmetric least squares (ALS) ###########################

unique_groups <- unique(smoothed_Raman_spectra$FileID)

als_baselined_Raman_spectra <- data.frame()
baseline_list <- list()

for (group_value in unique_groups) {
  group_data <- smoothed_Raman_spectra[smoothed_Raman_spectra$FileID == group_value, ]
  
  als_baseline <- asysm(group_data$V2_smoothed, lambda = 10^7, p = 0.01, maxit = 10)
  
  baseline_corrected <- group_data$V2_smoothed - als_baseline
  
  als_baselined_Raman_spectra <- rbind(als_baselined_Raman_spectra, 
                                      data.frame(V1 = group_data$V1,
                                                 V2_smoothed = group_data$V2_smoothed,
                                                 V2_baseline = baseline_corrected,
                                                 FileID = group_value))
  
  baseline_list[[group_value]] <- als_baseline
}

head(als_baselined_Raman_spectra)


############################### Polynomial #####################################

poly_baselined_Raman_spectra <- data.frame()
baseline_list <- list()

unique_groups <- unique(smoothed_Raman_spectra$FileID)

for (group_value in unique_groups) {
  group_data <- smoothed_Raman_spectra[smoothed_Raman_spectra$FileID == group_value, ]
  
  poly_baselined <- baseline_polynomial(group_data$V1, group_data$V2_smoothed, d = 2, tolerance = 0.02, stop = 1000) # Quadratic polynomial Weber et al 2021
  
  baseline_corrected <- group_data$V2_smoothed - poly_baselined$y
  
  poly_baselined_Raman_spectra <- rbind(poly_baselined_Raman_spectra, 
                                       data.frame(V1 = group_data$V1,
                                                  V2_smoothed = group_data$V2_smoothed,
                                                  V2_baseline = baseline_corrected,
                                                  FileID = group_value))
  
  baseline_list[[group_value]] <- poly_baselined
}

head(poly_baselined_Raman_spectra)


########################### Rolling ball (rb) ##################################

unique_groups <- unique(smoothed_Raman_spectra$FileID)

rb_baselined_Raman_spectra <- data.frame()
baseline_list <- list()

for (group_value in unique_groups) {
  group_data <- smoothed_Raman_spectra[smoothed_Raman_spectra$FileID == group_value, ]
  
  rb_baselined <- baseline_rollingball(group_data$V1, group_data$V2_smoothed, m = 201, s = 151)
  
  baseline_corrected <- signal_drift(group_data$V2_smoothed, lag = rb_baselined, subtract = TRUE)
  
  rb_baselined_Raman_spectra <- rbind(rb_baselined_Raman_spectra, 
                                     data.frame(V1 = group_data$V1,
                                                V2_smoothed = group_data$V2_smoothed,
                                                V2_baseline = baseline_corrected$y,
                                                FileID = group_value))
  
  baseline_list[[group_value]] <- rb_baselined
}

head(rb_baselined_Raman_spectra)


######## Statistics-sensitive non-linear iterative peak-clipping (SNIP) ########

unique_groups <- unique(smoothed_Raman_spectra$FileID)

SNIP_baselined_Raman_spectra <- data.frame()
baseline_list <- list()

for (group_value in unique_groups) {
  group_data <- smoothed_Raman_spectra[smoothed_Raman_spectra$FileID == group_value, ]
  
  SNIP_baselined <- baseline_snip(group_data$V1, group_data$V2_smoothed, LLS = FALSE, decreasing = FALSE, n = 100)
  
  baseline_corrected <- signal_drift(group_data$V2_smoothed, lag = SNIP_baselined, subtract = TRUE)
  
  SNIP_baselined_Raman_spectra <- rbind(SNIP_baselined_Raman_spectra, 
                                       data.frame(V1 = group_data$V1,
                                                  V2_smoothed = group_data$V2_smoothed,
                                                  V2_baseline = baseline_corrected$y,
                                                  FileID = group_value))
  
  baseline_list[[group_value]] <- SNIP_baselined
}

head(SNIP_baselined_Raman_spectra)


################################################################################
########################## Normalize trunc spectra #############################
################################################################################

normalize <- function(intensity) {
  (intensity - min(intensity)) / (sum(intensity) - min(intensity))
}

unique_groups <- unique(baselined_Raman_spectra$FileID)

normalized_Raman_spectra <- data.frame()

for (group_value in unique_groups) {
  group_data <- subset(baselined_Raman_spectra, FileID == group_value)
  
  normalized_values <- normalize(group_data$V2_baseline)
  
  normalized_Raman_spectra <- rbind(normalized_Raman_spectra, 
                                  data.frame(V1 = group_data$V1,
                                             V2_smoothed = group_data$V2_smoothed,
                                             V2_baseline = group_data$V2_baseline,
                                             V2_normalized = normalized_values,
                                             FileID = group_value))
}

head(normalized_Raman_spectra)


################################################################################
######################### Stack plots of processing ############################
################################################################################

par(mfrow = c(4, 1), mar = c(3, 3, 2, 1))

################################################################################

plot(Raman_spectra$V2 ~ Raman_spectra$V1, type = "n", xlab = "Wavenumber", ylab = "Intensity")

for (i in seq_along(unique_groups)) {
  group_value <- unique_groups[i]
  group_data <- subset(Raman_spectra, FileID == group_value)
  
  lines(group_data$V1, group_data$V2, col = "gray", lwd = 2)
}

title(main = "Original unprocessed data")

################################################################################

plot(smoothed_Raman_spectra$V2_smoothed ~ smoothed_Raman_spectra$V1, type = "n", xlab = "Wavenumber", ylab = "Intensity")

for (i in seq_along(unique_groups)) {
  group_value <- unique_groups[i]
  group_data <- subset(smoothed_Raman_spectra, FileID == group_value)
  
  lines(group_data$V1, group_data$V2_smoothed, col = "gray", lwd = 2)
}

title(main = "Smoothed data")

################################################################################

plot(baselined_Raman_spectra$V2_baseline ~ baselined_Raman_spectra$V1, type = "n", xlab = "Wavenumber", ylab = "Intensity")

for (i in seq_along(unique_groups)) {
  group_value <- unique_groups[i]
  group_data <- subset(baselined_Raman_spectra, FileID == group_value)
  
  lines(group_data$V1, group_data$V2_baseline, col = "gray", lwd = 2)
}

title(main = "Baselined data")

################################################################################

plot(normalized_Raman_spectra$V2_normalized ~ normalized_Raman_spectra$V1, type = "n", xlab = "Wavenumber", ylab = "Intensity")

for (i in seq_along(unique_groups)) {
  group_value <- unique_groups[i]
  group_data <- subset(normalized_Raman_spectra, FileID == group_value)
  
  lines(group_data$V1, group_data$V2_normalized, col = "gray", lwd = 2)
}

title(main = "Normalized data")


################################################################################
######################### Calculate area under curve ###########################
################################################################################

# Set wavenumber ranges

range_start_CD <- 2040
range_end_CD <- 2300
range_start_CH <- 2800
range_end_CH <- 3100

################################################################################

calculate_auc_CD <- function(group_data) {
  subset_data <- subset(group_data, V1 >= range_start_CD & V1 <= range_end_CD)
  
  baseline_auc <- trapz(c(range_start_CD, range_end_CD), c(subset_data$V2_normalized[1], subset_data$V2_normalized[nrow(subset_data)]))
  
  auc <- trapz(subset_data$V1, subset_data$V2_normalized)
  
  adjusted_auc <- auc - baseline_auc
  
  return(data.frame(AUC = adjusted_auc, Group = "CD"))
}

auc_CD_results <- normalized_Raman_spectra %>%
  group_by(FileID) %>%
  group_modify(~ calculate_auc_CD(.))

print(auc_CD_results)

################################################################################

calculate_auc_CH <- function(group_data) {
  subset_data <- subset(group_data, V1 >= range_start_CH & V1 <= range_end_CH)
  
  baseline_auc <- trapz(c(range_start_CH, range_end_CH), c(subset_data$V2_normalized[1], subset_data$V2_normalized[nrow(subset_data)]))
  
  auc <- trapz(subset_data$V1, subset_data$V2_normalized)
  
  adjusted_auc <- auc - baseline_auc
  
  return(data.frame(AUC = adjusted_auc, Group = "CH"))
}

auc_CH_results <- normalized_Raman_spectra %>%
  group_by(FileID) %>%
  group_modify(~ calculate_auc_CH(.))

print(auc_CH_results)

################################################################################

combined_auc_results <- bind_rows(auc_CD_results, auc_CH_results)

head(combined_auc_results)

Calculated_CD_results <- combined_auc_results %>%
  group_by(FileID) %>%
  mutate(Percentage_CD = ifelse((AUC[Group == "CD"] / sum(AUC)) * 100 < 0, 0, (AUC[Group == "CD"] / sum(AUC)) * 100)) %>%
  distinct(FileID, .keep_all = TRUE)

print(Calculated_CD_results)

write.csv(Calculated_CD_results, "Raman_CD_values.csv", row.names = FALSE)
