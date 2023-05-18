# Mel_SpatialSig
# Melanoma spatial signatures
The R version is available on CRAN.
The codes cover all the in-house scripts of our data analysis pipeline generating spatial gene signatures for response to immunotherapy in Melanoma. 
# System Requirements
Hardware Requirements

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

Software Requirements

The CRAN packages should be compatible with Windows, Mac, and Linux operating systems. The packages have been tested on macOS: Monterey version 12.6.3

Before setting up the R packages, users should have RSTUDIO-2023.03.1-446 or higher, and several packages set up from CRAN.

# Installation Guide:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PackageName")
