# List of required packages
packages <- c(
  "PatientProfiler", "tidyverse", "readxl", "ggeasy", "ggpubr", "broom", "fgsea", "ggplot2",
  "cowplot", "writexl", "igraph", "stringr", "httr", "plyr", "pheatmap", "purrr",
  "dendextend", "RColorBrewer", "broom", "GGally", "devtools", "ggfortify", 
  "ggrepel", "TCGAbiolinks", "dplyr", "readr", "survival", "survminer", 
  "gridExtra", "org.Hs.eg.db", "fs"
)

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to all packages
lapply(packages, install_if_missing)

# Load all packages
lapply(packages, library, character.only = TRUE)
