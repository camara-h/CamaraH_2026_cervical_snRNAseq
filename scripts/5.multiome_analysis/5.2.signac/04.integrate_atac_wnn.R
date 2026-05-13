#!/usr/bin/env Rscript

# Load necessary libraries
library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(optparse)
library(magrittr)
library(here)
library(purrr)
library(stringr)

# --- Command-line argument parsing ---
option_list <- list(
  make_option("--n_counts_cutoff", type = "integer", default = 750,
              help = "Minimum ATAC counts cutoff"),
  make_option("--TSS_cutoff", type = "numeric", default = 3,
              help = "TSS enrichment cutoff"),
  make_option("--mt_cutoff", type = "numeric", default = 2,
              help = "Mitochondrial percentage cutoff"),
  make_option("--UMI_cutoff", type = "integer", default = 500,
              help = "Minimum UMI cutoff"),
  make_option("--min_cell", type = "integer", default = 100,
              help = "Minimum cell threshold after QC to keep sample"),
  make_option("--kweight", type = "integer", default = 10,
              help = "k.weight parameter for ATAC integration"),
  make_option("--rna_dims", type = "integer", default = 20,
              help = "Number of Harmony dimensions to use for WNN"),
  make_option("--atac_dims", type = "integer", default = 30,
              help = "Number of integrated LSI dimensions to use for WNN"),
  make_option("--resolution", type = "double", default = 1,
              help = "Resolution for WNN clustering"),
  make_option("--data_dir", type = "character",
              default = "/home/camarah/storage/Tseng_Neck_Multiome_2026/03_cellranger/ref_2020/cellranger_arc/",
              help = "Folder with cellranger outputs"),
  make_option("--output_dir", type = "character",
              default = here("../04_outputs"),
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

n_counts_cutoff <- opt$n_counts_cutoff
TSS_cutoff <- opt$TSS_cutoff
mt_cutoff <- opt$mt_cutoff
UMI_cutoff <- opt$UMI_cutoff
min_cell_filter <- opt$min_cell
kweight <- opt$kweight
rna_dims <- opt$rna_dims
atac_dims <- opt$atac_dims
resolution <- opt$resolution
data_dir <- opt$data_dir
output_dir <- opt$output_dir

results_dir <- "/home/camarah/storage/Tseng_Neck_Multiome_2026/05_scripts/../04_outputs/750_atac_3_TSS_2_mt_500_UMI_gex_filtered"

filter <- paste(
  n_counts_cutoff, "_atac_",
  TSS_cutoff, "_TSS_",
  mt_cutoff, "_mt_",
  UMI_cutoff, "_UMI",
  sep = ""
)

tryCatch({
  
  # results_dir <- file.path(output_dir, filter)
  
  files <- list.files(
    results_dir,
    pattern = "_processed_seurat",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No Seurat RDS files found in '", results_dir, "' matching pattern '_processed_seurat'.")
  }
  
  message("Found ", length(files), " Seurat RDS files for integration.")
  
  filenames_only <- basename(files)
  sample_ids <- filenames_only %>%
    str_remove(pattern = "_processed_seurat.*")
  
  seurat_object_list <- files %>%
    set_names(sample_ids) %>%
    map(readRDS)
  
  seurat_object_list <- seurat_object_list[
    sapply(seurat_object_list, function(obj) !is.null(obj) && ncol(obj) > min_cell_filter)
  ]
  
  if (length(seurat_object_list) < 2) {
    stop("Fewer than 2 valid Seurat objects remain after filtering by min_cell.")
  }
  
  message("Proceeding with ", length(seurat_object_list), " objects.")
  
  # Make sure sample metadata is explicit and consistent
  for (i in seq_along(seurat_object_list)) {
    seurat_object_list[[i]]$sample <- names(seurat_object_list)[i]
  }
  
  # -------------------------------
  # 1. ATAC integration on per-sample objects
  # -------------------------------
  
  # Add sample names before Finding Integration Anchors
  for (nm in names(seurat_object_list)) {
    seurat_object_list[[nm]] <- RenameCells(
      seurat_object_list[[nm]],
      add.cell.id = nm
    )
  }
  
  for (i in seq_along(seurat_object_list)) {
    DefaultAssay(seurat_object_list[[i]]) <- "peaks"
  }
  
  # ------------------------------------
  # Optional: Find Anchors per pair to detect bad quality samples
  # ------------------------------------
  # pair_grid <- combn(sample_ids, 2, simplify = FALSE)
  # 
  # anchor_diag <- lapply(pair_grid, function(pair) {
  #   obj_pair <- seurat_object_list[pair]
  #   
  #   anchors_pair <- FindIntegrationAnchors(
  #     object.list = obj_pair,
  #     reduction = "rlsi",
  #     dims = 2:atac_dims
  #   )
  #   
  #   data.frame(
  #     sample1 = pair[1],
  #     sample2 = pair[2],
  #     n_anchors = nrow(anchors_pair@anchors)
  #   )
  # })
  # 
  # anchor_diag <- do.call(rbind, anchor_diag)
  # write.csv(anchor_diag, file.path(results_dir, "pairwise_anchor_counts.csv"), row.names = FALSE)
  # ----------------------------------
  # atac_anchors <- FindIntegrationAnchors(
  #   object.list = seurat_object_list,
  #   reduction = "rlsi",
  #   dims = 2:atac_dims
  # )

  # Merge with old reductions
  s.object <- merge(
    x = seurat_object_list[[1]],
    y = seurat_object_list[-1], 
    merge.dr = "TRUE"
  )
  
  # # Compute uncorrected merged-object LSI
  # DefaultAssay(s.object) <- "peaks"
  # s.object <- FindTopFeatures(s.object, min.cutoff = "q0")   # or your preferred cutoff
  # s.object <- RunTFIDF(s.object)
  # s.object <- RunSVD(s.object)
  
  # # Integrate ATAC embeddings onto merged object
  # atac_integrated <- IntegrateEmbeddings(
  #   anchorset = atac_anchors,
  #   reductions = s.object[["lsi"]],
  #   new.reduction.name = "integrated_lsi",
  #   dims.to.integrate = 2:atac_dims,
  #   k.weight = kweight
  # )
  
  # s.object[["integrated_lsi"]] <- atac_integrated[["integrated_lsi"]]
  
  # rm(atac_anchors, atac_integrated)
  
  # -------------------------------
  # 2. RNA integration with Harmony
  # -------------------------------
  # In your per-sample processing, SCT and PCA were already run.
  # We recompute PCA on the merged object so Harmony has a clean merged PCA input.
  DefaultAssay(s.object) <- "SCT"
  
  # Join layers if needed for Seurat v5 compatibility
  if (inherits(s.object[["SCT"]], "Assay5")) {
    s.object[["SCT"]] <- JoinLayers(s.object[["SCT"]])
  }
  
  # Harmony on sample
  s.object <- IntegrateLayers(
    object = s.object,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )
  
  message("Integrated PCA using Harmony")
  
  
  s.object <- FindNeighbors(s.object, reduction = "harmony", dims = 1:20)
  s.object <- FindClusters(s.object, resolution = .8, cluster.name = "harmony_clusters")
  
  s.object <- RunUMAP(s.object, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
  
  # Also save single-modality integrated UMAPs for QC
  s.object <- RunUMAP(
    s.object,
    reduction = "integrated_lsi",
    dims = 1:(atac_dims-1),
    reduction.name = "umap.integrated_lsi",
    reduction.key = "intLSIUMAP_"
  )
  
  s.object <- RunUMAP(
    s.object,
    reduction = "harmony",
    dims = 1:rna_dims,
    reduction.name = "umap.harmony",
    reduction.key = "harmonyUMAP_"
  )
  
  #Fix Annotation mistakes if carried
  DefaultAssay(s.object) <- "peaks"
  
  old <- seqlevels(Annotation(s.object))
  new <- gsub("^(chr)+", "chr", old)
  new[grepl(pattern = "M", x = new, ignore.case = FALSE)] <- "chrM"
  seqlevels(Annotation(s.object)) <- new
  
  DefaultAssay(s.object) <- "SCT"
  
  # -------------------
  # Infer RNA expression
  # -------------------
  DefaultAssay(s.object) <- "peaks"
  
  gene.activities <- GeneActivity(s.object)
  
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  s.object[['inferred.RNA']] <- CreateAssayObject(counts = gene.activities)
  DefaultAssay(s.object) <- "inferred.RNA"
  
  s.object <- NormalizeData(
    object = s.object,
    assay = 'inferred.RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(s.object$nCount_RNA)
  )
  
  
  saveRDS(s.object, file = file.path(results_dir, "integrated_seurat.rds"))
  
  message("Revised WNN integration script finished successfully.")
  
}, error = function(e) {
  message("Failed to integrate")
  message("Error message: ", e$message)
  # quit(save = "no", status = 1)
})

# basic checks per object
for (i in seq_along(seurat_object_list)) {
  cat("\n--- sample", i, "---\n")
  print(seurat_object_list[[i]])
  cat("assays:", Assays(seurat_object_list[[i]]), "\n")
  
  # for Seurat v5 layered assays
  if ("peaks" %in% Assays(seurat_object_list[[i]])) {
    cat("layers in peaks:\n")
    print(Layers(seurat_object_list[[i]][["peaks"]]))
    
    cat("counts layer exists? ")
    print("counts" %in% Layers(seurat_object_list[[i]][["peaks"]]))
    
    x <- try(LayerData(seurat_object_list[[i]][["peaks"]], layer = "counts"), silent = TRUE)
    if (!inherits(x, "try-error") && !is.null(x)) {
      cat("dim counts:", dim(x), "\n")
      cat("nonzero entries:", length(x@x), "\n")
    } else {
      cat("counts layer missing or unreadable\n")
    }
  }
}

sapply(seurat_object_list, ncol)

for (i in seq_along(seurat_object_list)) {
  cat("\n--- sample", i, "---\n")
  if ("peaks" %in% Assays(seurat_object_list[[i]])) {
    x <- try(LayerData(seurat_object_list[[i]][["peaks"]], layer = "counts"), silent = TRUE)
    if (!inherits(x, "try-error") && !is.null(x)) {
      cat("cells with >0 peak counts:", sum(Matrix::colSums(x) > 0), "\n")
      cat("median peak counts/cell:", median(Matrix::colSums(x)), "\n")
    }
  }
}


# basic summary table
sample_qc <- purrr::imap_dfr(seurat_object_list, function(obj, nm) {
  x <- LayerData(obj[["peaks"]], layer = "counts")
  tibble::tibble(
    sample = nm,
    n_cells = ncol(obj),
    n_peak_features = nrow(x),
    nnzero = length(x@x),
    nnzero_per_cell = length(x@x) / ncol(obj),
    median_counts_per_cell = median(Matrix::colSums(x)),
    median_peaks_per_cell = median(Matrix::colSums(x > 0))
  )
})

sample_qc