library(Seurat)
library(tidyverse)

# -----------------------------------------------------------------------------
# Create helper functions
# -----------------------------------------------------------------------------

standardize_ref_sample <- function(x) {
  x <- toupper(x)
  x <- sub("_GEX$", "", x)
  x <- sub("^NEVIS014_PP$", "NEVIS_014", x)
  x
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

seurat_ref_path <- "/home/camarah/storage/HCA_adipose/Objects/adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_modified_cleaned.rds"
output_path <- "../04_outputs/reference_gex"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
stopifnot(file.exists(seurat_ref_path))

# -----------------------------------------------------------------------------
# Load objects
# -----------------------------------------------------------------------------

seurat_ref <- readRDS(seurat_ref_path)


# -----------------------------------------------------------------------------
# Parse barcode and sample from the reference cell names
# -----------------------------------------------------------------------------

ref_cells <- seurat_ref@meta.data$cells

ref_barcode <- sub("_.*$", "", ref_cells)
ref_sample_raw <- sub("^[^_]+_", "", ref_cells)
ref_sample <- standardize_ref_sample(ref_sample_raw)


ref_canonical <- paste(ref_sample, ref_barcode, sep = "_")

# -----------------------------------------------------------------------------
# Store them in metadata
# -----------------------------------------------------------------------------
seurat_ref$barcode_raw <- ref_barcode
seurat_ref$sample_std <- ref_sample
seurat_ref$cell_canonical <- ref_canonical


# -----------------------------------------------------------------------------
# Extract, clean and save metadata as RDS
# -----------------------------------------------------------------------------
seurat_meta <- seurat_ref@meta.data
seurat_meta <- seurat_meta %>% select(starts_with("cell"), "clean_sample", "neck_region", "barcode_raw", "sample_std")

saveRDS(seurat_meta, file.path(output_path, "seurat_ref_meta.rds"))













