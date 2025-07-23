# List of required packages
packages <- c(
  "PatientProfiler", "tidyverse", "readxl", "ggeasy", "ggpubr", "fgsea", "ggplot2",
  "cowplot", "writexl", "igraph", "stringr", "httr", "plyr", "pheatmap",
  "dendextend", "RColorBrewer", "broom", "GGally", "devtools", "ggfortify",
  "ggrepel", "AnnotationDbi", "org.Hs.eg.db", "survival", "survminer", "gridExtra"
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
