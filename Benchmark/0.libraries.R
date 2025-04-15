

# List of required packages
required_packages <- c(
  "PatientProfiler", "tidyverse", "readxl", "ggeasy", "ggpubr", "fgsea", "ggplot2", 
  "cowplot", "writexl", "igraph", "stringr", "httr", "plyr", "pheatmap", "dendextend", 
  "RColorBrewer", "broom", "GGally", "devtools", "ggfortify", "ggrepel", "Rtsne", 
  "stats", "cluster", "gprofiler2", "MOFA2"
)

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {  # Check if the package is installed
    if (pkg %in% rownames(available.packages())) {
      install.packages(pkg)  # Install from CRAN if available
    } else {
      BiocManager::install(pkg, ask = FALSE, update = TRUE)  # Install from Bioconductor if needed
    }
  }
}

# Ensure that BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing packages
lapply(required_packages, install_if_missing)

# Load all packages
lapply(required_packages, library, character.only = TRUE)








