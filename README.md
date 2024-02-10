# Mel_SpatialSig
Welcome to the Mel_SpatialSig repository, where we focus on leveraging R for the analysis of melanoma spatial signatures. This project is dedicated to harnessing spatial statistical analysis to better predict the outcome of immunotherapy treatment in advanced melanoma patients. Our goal is to contribute to advancements in melanoma research by providing tools for more precise diagnosis and the development of targeted treatment strategies.

# Project Overview
Purpose: We aim to improve the prediction of response or resistance to immunotherapies in melanoma patients. This goal is based on the hypothesis that current gene signatures predicting immunotherapy outcomes show only modest accuracy due to the lack of spatial information about cellular functions and molecular processes within tumors and their microenvironment. Experimental Design: we collected gene expression data spatially from three cellular compartments defined by CD68+macrophages, CD45+leukocytes and S100B+tumor cells in 55-immunotherapy-treated melanoma specimens using Digital Spatial Profiling-Whole Transcriptome Atlas (DSP-WTA). We developed a computational pipeline to discover compartment-specific gene signatures and determine if adding spatial information can improve patient stratification. Results: We achieved robust performance of compartment-specific signatures in predicting the outcome to ICI in the discovery cohort. Of the three signatures, S100B signature showed the best performance in the validation cohort (N=45). We also compared our compartment-specific signatures with published bulk signatures and found the S100B tumor spatial signature outperformed previous signatures. Within the 8-gene S100B signature, 5 genes (PSMB8, TAX1BP3, NOTCH3, LCP2, NQO1) with positive coefficients predict the response and 3 genes (KMT2C, OVCA2, MGRN1) with negative coefficients predict the resistance to treatment. Conclusion: We conclude that the spatially defined compartment signatures utilize tumor and TME-specific information, leading to more accurate prediction of treatment outcome, and thus merit prospective clinical assessment.

# Prerequisites
Before you begin, ensure you have the following requirements:
•	R (version 4.0.0 or newer)
•	Git for cloning the repository

# Installation Instructions
1.	Install R:
Download and install R from the Comprehensive R Archive Network (CRAN). 
install.packages(packages)

# How to Use 
To use the scripts, simply run them within an R environment. Each script is documented with comments explaining its purpose, input parameters, and output formats.

# Contributing 
Contributions to the Mel_SpatialSig project are welcome. To contribute, please fork the repository, create your feature branch, commit your changes, and push to the branch. Finally, open a pull request for review.

# Acknowledgments 
This project has been made possible thanks to the dedicated efforts of contributors. We appreciate all contributions that help improve the Mel_SpatialSig project.

# Contact
For more information or to inquire about collaboration opportunities, please contact thazin.aung@yale.edu or tznaung@gmail.com.
