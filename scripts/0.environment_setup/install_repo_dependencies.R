# Purpose: Install and validate project dependencies tracked with renv
# Run after: renv has been initialized for the project
# Inputs: none
# Outputs: installed packages in the project environment
# Must not change: package versions, installation order where intentional, and Bioconductor version handling
# Cleanup level: annotation and readability only, with no intended effect on installed package set

# -------------------------------------------------------------------------
# Configure Bioconductor behavior for this R version
# -------------------------------------------------------------------------

# When running under R 4.3.1, force Bioconductor 3.17 repositories.
# This helps avoid package resolution drift and GitHub/Bioconductor conflicts.
if (getRversion() == "4.3.1") {
  options(repos = BiocManager::repositories(version = "3.17"))
}

# Prevent remotes/renv from failing on warnings during installation
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# -------------------------------------------------------------------------
# Install core Bioconductor packages through renv for tracking
# -------------------------------------------------------------------------

# Ensure BiocManager itself is tracked by renv
renv::install("BiocManager")

# Main Bioconductor dependency used in the project
renv::install("DESeq2")

# -------------------------------------------------------------------------
# Install CRAN dependencies
# -------------------------------------------------------------------------

# NMF >= 0.23.0 is usually available from CRAN
renv::install("NMF")

# -------------------------------------------------------------------------
# Install Bioconductor packages that often need explicit version control
# -------------------------------------------------------------------------

# These packages have previously tended to update toward Bioconductor 3.18,
# so they are installed explicitly against version 3.17.
problem_packages <- c(
  "basilisk",
  "basilisk.utils",
  "ComplexHeatmap",
  "ensembldb",
  "OmnipathR",
  "treeio"
)

BiocManager::install("basilisk", version = "3.17", ask = FALSE, force = TRUE)
BiocManager::install("basilisk.utils", version = "3.17", ask = FALSE, force = TRUE)
BiocManager::install("treeio", version = "3.17", ask = FALSE, force = TRUE)
BiocManager::install("ComplexHeatmap", version = "3.17")
BiocManager::install("ensembldb", version = "3.17")
BiocManager::install("OmnipathR", version = "3.17")

# Annotation databases
BiocManager::install("org.Hs.eg.db", version = "3.17")
BiocManager::install("hgu133plus2.db", version = "3.17")

# Check whether installed Bioconductor packages are consistent
BiocManager::valid()

# -------------------------------------------------------------------------
# Install GitHub packages through renv for version locking
# -------------------------------------------------------------------------

# circlize and ComplexHeatmap
renv::install("jokergoo/circlize")
renv::install("jokergoo/ComplexHeatmap")

# NicheNet
renv::install("saeyslab/nichenetr")

# VISION
renv::install("YosefLab/VISION")

# CellChat. Installed last due to more complex dependencies
renv::install("sqjin/CellChat")

# LIANA. Left commented as in the original script
# renv::install("saezlab/liana")

# GPTCelltype
renv::install("Winnie09/GPTCelltype")

# GSVA
renv::install("rcastelo/GSVA@devel")

# -------------------------------------------------------------------------
# Final message
# -------------------------------------------------------------------------
message("✅ All GitHub and Bioconductor packages attempted. Review output for any issues.")