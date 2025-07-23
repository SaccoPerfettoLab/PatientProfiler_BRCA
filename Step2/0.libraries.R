# List of required packages
packages <- c(
  "PatientProfiler", "broom", "cowplot", "dendextend", "devtools", "fgsea",
  "GGally", "ggeasy", "ggfortify", "ggplot2", "ggpubr", "ggrepel", "httr",
  "igraph", "plyr", "pheatmap", "readxl", "rstatix", "RColorBrewer",
  "stringr", "tidyverse", "writexl"
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
