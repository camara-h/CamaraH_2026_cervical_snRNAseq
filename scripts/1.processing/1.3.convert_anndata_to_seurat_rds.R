# Purpose: Convert an AnnData (.h5ad) object to a Seurat-compatible RDS object
# Run after: required packages and Python anndata support are available
# Inputs: AnnData file path provided by command line or manually inside the script
# Outputs: Seurat object saved as .rds
# Must not change: conversion behavior, default output naming, and current manual batch-processing behavior
# Cleanup level: annotation and readability only, with no intended effect on outputs

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(optparse)
library(sceasy)
library(reticulate)

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option("--input", type = "character", help = "Path to the AnnData object."),
  make_option("--output", type = "character", help = "Output file name", default = NULL)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -------------------------------------------------------------------------
# Store command-line inputs
# -------------------------------------------------------------------------
path <- opt$input
h5_output <- opt$output

# -------------------------------------------------------------------------
# Manual file inspection and ad hoc conversion targets
# -------------------------------------------------------------------------
# Kept as in the original script because it may have been useful during
# interactive execution and testing.
list.files("~/Projects/Joy_2025/data/Zhang_2025/")

# -------------------------------------------------------------------------
# Convert one or more AnnData objects
# -------------------------------------------------------------------------
for (f in c(
  # here("output/1.processing/data/cleaned_objects/Zhang_bat_50pct_brown_50pct_white_adipo_100k.h5ad"),
  # here("output/1.processing/data/cleaned_objects/Zhang_bat_brown_white_adipo_100k.h5ad"),
  # "~/Projects/Joy_2025/data/Zhang_2025/GSE247719_PanSci_08_iWAT_adata_young_wt_male.h5ad",
  # "~/Projects/Joy_2025/data/Zhang_2025/GSE247719_PanSci_09_gWAT_adata_young_wt_male.h5ad",
  # "~/Projects/Joy_2025/data/Zhang_2025/GSE247719_PanSci_07_BAT_adata_young_wt_male.h5ad",
  "~/Projects/Joy_2025/data/Zhang_2025/GSE247719_PanSci_07_BAT_adata_young_wt.h5ad"
)) {
  path <- f
  
  dirpath <- dirname(path)
  filename_ext <- basename(path)
  filename_raw <- tools::file_path_sans_ext(filename_ext)
  
  print(filename_raw)
  
  # Ensure Python anndata is available before conversion
  reticulate::py_require("anndata")
  reticulate::import("anndata")
  
  # If running from command line and no explicit output was provided,
  # save the RDS next to the input file using the same base name.
  if (exists("opt")) {
    if (is.null(opt$output)) {
      h5_output <- paste(dirpath, "/", filename_raw, ".rds", sep = "")
    }
  } else {
    h5_output <- paste(dirpath, "/", filename_raw, ".rds", sep = "")
  }
  
  # Alternative output path logic kept commented exactly as in the original
  # dir.create(here("output/1.processing/data/cleaned_objects"), recursive = T)
  # h5_output <- here(paste0("output/1.processing/data/cleaned_objects/", filename_raw, ".rds"))
  
  sceasy::convertFormat(
    path,
    from = "anndata",
    to = "seurat",
    outFile = h5_output
  )
}