library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input paths
# -----------------------------------------------------------------------------

seurat_path <- here("data/9.multiome_data/multiome_seurat.rds")
seurat_ref_path <- here("data/adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_cleaned.rds")

stopifnot(file.exists(seurat_path))
stopifnot(file.exists(seurat_ref_path))

# -----------------------------------------------------------------------------
# Load objects
# -----------------------------------------------------------------------------

query_obj <- readRDS(seurat_path)
ref_obj <- readRDS(seurat_ref_path)

message("Loaded query object with ", ncol(query_obj), " cells.")
message("Available assays: ", paste(Assays(query_obj), collapse = ", "))
message("Available reductions: ", paste(Reductions(query_obj), collapse = ", "))

# -----------------------------------------------------------------------------
# Normalize peak annotation seqlevels to chr-prefixed style
# This helps keep annotation naming consistent with genome / peak conventions.
# -----------------------------------------------------------------------------

stopifnot("peaks" %in% Assays(query_obj))
DefaultAssay(query_obj) <- "peaks"
ann <- Annotation(query_obj)
seqlevels(ann) <- paste0("chr", gsub("^(chr)+", "", seqlevels(ann)))
Annotation(query_obj[["peaks"]]) <- ann

# -----------------------------------------------------------------------------
# Build reference metadata table for transfer
# -----------------------------------------------------------------------------

required_ref_cols <- c("clean_sample", "cells", "neck_region")
missing_ref_cols <- setdiff(required_ref_cols, colnames(ref_obj@meta.data))
if (length(missing_ref_cols) > 0) {
  stop("Reference object is missing required metadata columns: ",
       paste(missing_ref_cols, collapse = ", "))
}

ref_cell_labels <- ref_obj@meta.data %>%
  select(contains("cell_type"), clean_sample, cells, neck_region) %>%
  mutate(clean_sample = str_replace(clean_sample, "Nevis014_PP", "NEVIS_014"), # Reformat Nevis014_PP on reference
           barcode = paste(clean_sample, str_remove(cells, "_.*"), sep = "_")) %>%
  select(contains("cell_type"), barcode, neck_region)

# Match reference metadata to query object barcodes, using case-insensitive matching.
idx <- match(toupper(colnames(query_obj)), toupper(ref_cell_labels$barcode))

if (any(is.na(idx))) {
  print(sample(colnames(query_obj)[is.na(idx)],size = 10))
  stop("Some query barcodes were not found in the reference metadata.")
}

ref_cell_labels <- ref_cell_labels[idx, , drop = FALSE]

if (!all(toupper(ref_cell_labels$barcode) == toupper(colnames(query_obj)))) {
  stop("Reference metadata reordering failed. Barcodes are not aligned.")
}

# -----------------------------------------------------------------------------
# Append only new metadata columns to avoid silent overwriting or duplicate names
# -----------------------------------------------------------------------------

query_meta <- query_obj@meta.data
new_cols <- setdiff(colnames(ref_cell_labels), colnames(query_meta))
query_meta <- cbind(query_meta, ref_cell_labels[, new_cols, drop = FALSE])
query_obj@meta.data <- query_meta

# -----------------------------------------------------------------------------
# Output directory
# -----------------------------------------------------------------------------

out_dir_fig <- file.path(dirname(seurat_path), "output", "figures")
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Set identities for downstream plotting
# -----------------------------------------------------------------------------

stopifnot("cell_type_short" %in% colnames(query_obj@meta.data))
Idents(query_obj) <- "cell_type_short"
# -----------------------------------------------------------------------------
# Marker detection
# NOTE: This uses the current active identities already present in the object.
# If markers should be computed by transferred labels, set Idents() first.
# -----------------------------------------------------------------------------

# query_obj <- PrepSCTFindMarkers(query_obj)
# DefaultAssay(query_obj) <- "RNA"
# mrk_sct <- FindAllMarkers(query_obj, logfc.threshold = 1, min.pct = 0.2)

# -----------------------------------------------------------------------------
# UMAP plots
# -----------------------------------------------------------------------------

umap_reductions <- grep(
  x = Reductions(query_obj),
  pattern = "umap",
  ignore.case = TRUE,
  value = TRUE
)

stopifnot(length(umap_reductions) > 0)
stopifnot("sample" %in% colnames(query_obj@meta.data))

for (reduction in umap_reductions) {
  p <- DimPlot(query_obj, reduction = reduction, group.by = "sample")
  print(p)
  ggsave(
    filename = file.path(out_dir_fig, paste0("dimplot_by_sample_", reduction, ".pdf")),
    plot = p,
    width = 150,
    height = 100,
    units = "mm",
    bg = "white"
  )
}

for (reduction in umap_reductions) {
  p <- DimPlot(query_obj, reduction = reduction)
  print(p)
  ggsave(
    filename = file.path(out_dir_fig, paste0("dimplot_by_active_ident_", reduction, ".pdf")),
    plot = p,
    width = 150,
    height = 100,
    units = "mm",
    bg = "white"
  )
}

# -----------------------------------------------------------------------------
# QC plots
# -----------------------------------------------------------------------------

qc_features <- c("nCount_RNA", "nCount_peaks", "TSS.enrichment", "nucleosome_signal", "percent.mt")

p_qc_feature <- FeaturePlot(query_obj, qc_features)
p_qc_vln <- VlnPlot(query_obj, qc_features)

ggsave(
  filename = file.path(out_dir_fig, "featureplot_qc_metrics.pdf"),
  plot = p_qc_feature,
  width = 300,
  height = 150,
  units = "mm",
  bg = "white"
)

ggsave(
  filename = file.path(out_dir_fig, "vln_qc_metrics.pdf"),
  plot = p_qc_vln,
  width = 300,
  height = 150,
  units = "mm",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Coverage plots
# -----------------------------------------------------------------------------

stopifnot("peaks" %in% Assays(query_obj))
stopifnot("SCT" %in% Assays(query_obj))
stopifnot("cell_type_short" %in% colnames(query_obj@meta.data))

# ---- Plot coverage of canonical genes by cell_type_short ----
DefaultAssay(query_obj) <- "peaks"
Idents(query_obj) <- "cell_type_short"


features <- c(
  "ADIPOQ", "PPARGC1A", "PDGFRA", "DCN", "JAM2", "SHANK3", "PROX1",
  "STEAP4", "MYOCD", "CD163", "CPA3", "IL7R", "MS4A1", "CSF3R",
  "CDH19", "CNTNAP5", "PTH", "PRM1"
)

idents <- c("BAds", "WAds", "ASPC", "EndoCs", "SMCs", "Macro")

for (goi in features) {
  p <- CoveragePlot(
    query_obj,
    region = goi,
    features = goi,
    assay = "peaks",
    expression.assay = "RNA",
    idents = idents,
    peaks = FALSE) &
    scale_fill_manual(values = palette.use)
  
  ggsave(
    filename = file.path(out_dir_fig, paste0("coverage_", goi, ".pdf")),
    plot = p,
    width = 150,
    height = 150,
    units = "mm",
    bg = "white"
  )
}

# ---- Plot coverage of HOX genes in White Adipocytes from deep vs superficial ----
# save current identities
old_idents <- Idents(query_obj)

# Create a new identity. This is a workaround on defective subsetting of Seurat v5 for multiome
query_obj$cell_type_neck_region <- paste(query_obj$cell_type_short, query_obj$neck_region, sep = "_")

# Assing and select cells to plot
Idents(query_obj) <- "cell_type_neck_region"

idents_query <- paste0(c("WAds"), 
                       collapse = "|")

idents_to_plot <- grep(idents_query, unique(query_obj$cell_type_neck_region), value = T)
idents_to_plot <- grep(pattern = "Intermediate", idents_to_plot, value = T, invert = T)

# Select genes
features <- c("HOXA3", "HOXA4", "HOXA5", "HOXA6",
              "PAX3", "IRX2", "IRX5", "HOXD3", "HOXD4")

# Set palette
id_levels <- levels(Idents(query_obj))
neck_region_from_ident <- sub("^.*_", "", id_levels)
ident_colors <- setNames(neck.color[neck_region_from_ident], id_levels)
ident_colors

for (goi in features){
  
  scale_factor <- 2
  
# plot only the WAds-derived identity groups
  
  p <- CoveragePlot(
    query_obj,
    region = goi,
    features = goi,
    assay = "peaks",
    expression.assay = "SCT",
    peaks = TRUE,
    idents = idents_to_plot
  ) &
    scale_fill_manual(values = ident_colors) &
    theme(
      text = element_text(size = 5 * scale_factor),
      axis.title = element_text(size = 5 * scale_factor),
      axis.text = element_text(size = 5 * scale_factor),
      strip.text = element_text(size = 5 * scale_factor),
      legend.title = element_text(size = 5 * scale_factor),
      legend.text = element_text(size = 5 * scale_factor),
      axis.line = element_line(linewidth = 0.3 * scale_factor),
      axis.ticks = element_line(linewidth = 0.3 * scale_factor),
      # panel.border = element_rect(fill = NA, linewidth = 0.3 * scale_factor)
    )

ggsave(
  filename = file.path(OUTPUT_DIR, paste0("4G_coverage_", goi, "_", idents_query,"_by_neck.region.pdf")),
  plot = p,
  width = 65 * scale_factor,
  height = 50 * scale_factor,
  units = "mm",
  bg = "white"
)

}
# restore original identities afterward
Idents(query_obj) <- old_idents


