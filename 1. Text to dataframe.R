################################################################################
########################## Notes on provided code ##############################
################################################################################

# The code provide is intended for converting Raman spectra text files (.txt) to
# a dataframe for further processing. The start/end row allows for a specific 
# range of wavenumbers to be selected from the spectra, if desired. Otherwise
# list all the row number in the text file for the whole spectra.

################################################################################
########################## 0% txt files whole spectra ##########################
################################################################################

library(tidyverse)

setwd("C:\\Users\\location of folder with txt files")

read_specific_lines <- function(file, start_row, end_row) {
  data <- read.table(file, skip = start_row - 1, nrows = end_row - start_row + 1)
  return(data)
}

files <- list.files(pattern = "\\.txt$")

data_frames <- list()

# Define the specific range of rows you want to extract. Either all rows for the
# whole spectra or a truncated number of rows.
start_row <- 1
end_row <- 958

for (file in files) {
  selected_data <- read_specific_lines(file, start_row, end_row)
  selected_data$FileName <- file
  file_parts <- unlist(strsplit(gsub(".txt", "", file), "_"))
  for (i in 1:length(file_parts)) {
    col_name <- paste("Part", i, sep = "")
    selected_data[[col_name]] <- file_parts[i]
  }
  data_frames[[file]] <- selected_data
}

# Combine all data frames into a single data frame using dplyr's bind_rows
combined_whole_0pct_data <- dplyr::bind_rows(data_frames, .id = "FileID")

# Write the combined data to a CSV file
write.csv(combined_whole_0pct_data, "0pct_whole_spectra.csv", row.names = FALSE)

# Write the combined data to a rdata file
save(combined_whole_0pct_data, file="0pct_whole_spectra.rdata")


