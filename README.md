# Raman-spectra-processing
R code written for processing of Raman data using tidyverse (https://github.com/tidyverse), plotly (https://github.com/plotly), alkehest (https://github.com/tesselle/alkahest), pracma (https://github.com/cran/pracma), and pspline (https://github.com/cran/pspline) packages. This repository provides three R files for the analysis of Raman spectra:

1. Assuming Raman spectra has been saved as text files, R code 1 provides a function to convert all text files in a folder into a single dataframe that can be saved as a .Rdata or .csv file. This allows for all the unique spectra to be localized into one dataframe for processing. Additionally, two lines of code in the function can be edited to allow for the spectra to be truncated if desired.

2. Once the spectra have been converted into a dataframe, R code 2 provides several functions for processing of spectra:
  a. Plotly function that allows for outlier spectra, such as burnt cells, to be identified so the spectra can be removed from the dataframe.
  b. The Savitsky-Goley smoothing algorithm for removing noise from data.
  c. Several  algorithms for baseline removal of spectra.
  d. Scale-invariant sum normalization function to normalize all spectra to each other.
  e. Functions to calculate area under the curve for specific wavenumber ranges. The code is currently written for the comparision of the carbon-deuterium peak in the silent region of the spectra to the carbon-hydrogen peak.
  f. Calculaion of atom percent of deuterium (i.e., (CD /(CD + CH))*100).

3. If desired, R code 3 provides functions to calculate wavenumber ranges using the 2nd derivative.

![image](https://github.com/user-attachments/assets/11c375ee-17d9-4c31-bc3f-6312a094e832)
