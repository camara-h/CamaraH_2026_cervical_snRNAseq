#!/usr/bin/env Rscript

# =========================================================
# Process per-sample multiome objects using:
#   1. a shared reduced peak universe
#   2. cached QC metrics
#   3. optional restriction to cells retained in the GEX reference
#
# This version is intentionally dedicated to the GEX-filtered workflow
# to keep logic simpler and easier to maintain.
# =========================================================

# ----------------------------
# Load required libraries
# ----------------------------
library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(optparse)
library(magrittr)
library(here)

# ----------------------------
# Parse command line arguments
# ----------------------------
option_list <- list(
  make_option("--n_counts_cutoff", type = "integer", default = 750,
              help = "Minimum ATAC counts cutoff [default %default]"),
  make_option("--TSS_cutoff", type = "numeric", default = 3,
              help = "Minimum TSS enrichment cutoff [default %default]"),
  make_option("--mt_cutoff", type = "numeric", default = 2,
              help = "Maximum mitochondrial percentage [default %default]"),
  make_option("--UMI_cutoff", type = "integer", default = 500,
              help = "Minimum RNA UMI cutoff [default %default]"),
  make_option("--min_cell", type = "integer", default = 100,
              help = "Minimum number of cells required to retain a sample [default %default]"),
  make_option("--data_dir", type = "character",
              default = "/home/camarah/storage/Tseng_Neck_Multiome_2026/03_cellranger/ref_2020/cellranger_arc/",
              help = "Folder containing Cell Ranger ARC outputs [default %default]"),
  make_option("--output_dir", type = "character",
              default = here("../04_outputs"),
              help = "Base output folder [default %default]"),
  make_option("--master_peak_path", type = "character",
              default = here("../04_outputs/master_peak_universe_gex_filtered/reduced_ranges_gex_filtered.rds"),
              help = "Path to shared reduced peak universe RDS [default %default]"),
  make_option("--qc_cache_path", type = "character",
              default = here("../04_outputs/qc_cache/qc_cache.rds"),
              help = "Path to cached per-cell QC metrics [default %default]"),
  make_option("--seurat_ref_meta_path", type = "character",
              default = here("../04_outputs/reference_gex/seurat_ref_meta.rds"),
              help = "Path to cleaned reference metadata RDS [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------
# Store parsed arguments
# ----------------------------
n_counts_cutoff <- opt$n_counts_cutoff
TSS_cutoff <- opt$TSS_cutoff
mt_cutoff <- opt$mt_cutoff
UMI_cutoff <- opt$UMI_cutoff
min_cell_filter <- opt$min_cell
data_dir <- opt$data_dir
output_dir <- opt$output_dir
master_peak_path <- opt$master_peak_path
qc_cache_path <- opt$qc_cache_path
seurat_ref_meta_path <- opt$seurat_ref_meta_path

output_suffix <- "_gex_filtered"

filter_label <- paste0(
  n_counts_cutoff, "_atac_",
  TSS_cutoff, "_TSS_",
  mt_cutoff, "_mt_",
  UMI_cutoff, "_UMI",
  output_suffix
)

results_dir <- file.path(output_dir, filter_label)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load gene annotation
# ----------------------------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

old_seqlevels <- seqlevels(annotation)
new_seqlevels <- paste0("chr", gsub("^(chr)+", "", old_seqlevels))
new_seqlevels[old_seqlevels == "MT"] <- "chrM"
seqlevels(annotation) <- new_seqlevels

message("Annotation loaded.")
message("Example seqlevels: ", paste(head(seqlevels(annotation)), collapse = ", "))

# ----------------------------
# Load shared resources
# ----------------------------
master_peak_universe <- readRDS(master_peak_path)
qc_cache <- readRDS(qc_cache_path)
ref_meta <- readRDS(seurat_ref_meta_path)

# Validate required reference metadata columns
required_ref_cols <- c("barcode_raw", "sample_std", "cell_canonical")
missing_ref_cols <- setdiff(required_ref_cols, colnames(ref_meta))

if (length(missing_ref_cols) > 0) {
  stop("Reference metadata is missing required columns: ",
       paste(missing_ref_cols, collapse = ", "))
}

# ----------------------------
# Initialize tracking objects
# ----------------------------
successful_samples <- character()
failed_samples <- character()

sample_timing <- data.frame(
  sample = character(),
  status = character(),
  elapsed_sec = numeric(),
  elapsed_min = numeric(),
  start_time = character(),
  end_time = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)

# ----------------------------
# Discover available samples
# ----------------------------
samples <- list.dirs(path = data_dir, full.names = FALSE, recursive = FALSE)

# ----------------------------
# Process each sample independently
# ----------------------------
for (sample.id in samples) {
  
  message("Processing sample: ", sample.id)
  
  start_time <- Sys.time()
  sample_status <- "failed"
  sample_error <- NA_character_
  
  tryCatch({
    
    # ----------------------------------------
    # 1. Read QC metrics for this sample
    # ----------------------------------------
    sample_qc <- qc_cache[qc_cache$sample == sample.id, ]
    
    if (nrow(sample_qc) == 0) {
      stop("No QC cache rows found for sample.")
    }
    
    qc_kept_cells <- sample_qc$raw_cell[
      sample_qc$nCount_peaks > n_counts_cutoff &
        sample_qc$TSS.enrichment >= TSS_cutoff &
        sample_qc$percent.mt < mt_cutoff &
        sample_qc$nCount_RNA > UMI_cutoff &
        sample_qc$nCount_peaks < 100000
    ]
    
    if (length(qc_kept_cells) == 0) {
      stop("No cells passed QC thresholds.")
    }
    
    if (length(qc_kept_cells) < min_cell_filter) {
      stop("Less than ", min_cell_filter, " cells remain after QC-passing. Discarting sample: ", sample.id)
    }
    
    message(length(qc_kept_cells), " cells passed QC in sample ", sample.id)
    
    # ----------------------------------------
    # 2. Read raw matrices and fragment path
    # ----------------------------------------
    h5_path <- file.path(data_dir, sample.id, "outs", "filtered_feature_bc_matrix.h5")
    fragpath <- file.path(data_dir, sample.id, "outs", "atac_fragments.tsv.gz")
    output_path <- file.path(results_dir, paste0(sample.id, "_processed_seurat", output_suffix, ".rds"))
    
    if (file.exists(output_path)) {
      stop("Output already exists: ", output_path)
    }
    
    counts <- Read10X_h5(h5_path)
    
    # ----------------------------------------
    # OPTIONAL: Read CellBender matrices and replace gene expression counts
    # ----------------------------------------
    # cellbender <- paste("./results/Cellbender_OutputFiles/", sample.id, "/raw_feature_bc_matrix_output_filtered.h5", sep = "")
    # cellbender_counts <- Read_CellBender_h5_Mat(cellbender)
    # cellbender_counts <- cellbender_counts[, colnames(cellbender_counts) %in% colnames(counts$`Gene Expression`)]
    # 
    # counts$`Gene Expression` <- cellbender_counts
    
    # ----------------------------------------
    # 3. Build RNA object
    # ----------------------------------------
    s.object <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA"
    )
    
    # ----------------------------------------
    # 4. Restrict to cells retained in GEX reference
    # ----------------------------------------
    ref_barcodes_this_sample <- ref_meta$barcode_raw[ref_meta$sample_std == sample.id]
    
    if (length(ref_barcodes_this_sample) == 0) {
      stop("No reference barcodes found for sample in reference metadata.")
    }
    
    # Final retained cells must:
    #   - exist in RNA matrix
    #   - exist in ATAC matrix
    #   - pass QC thresholds
    #   - be present in the reference GEX set
    common_cells <- Reduce(
      intersect,
      list(
        colnames(s.object),
        colnames(counts$Peaks),
        # qc_kept_cells, # Only included if you ALSO want to filter by ATAC QC
        ref_barcodes_this_sample
      )
    )
    
    if (length(common_cells) < min_cell_filter) {
      stop("Not enough cells remain after intersecting RNA, ATAC, QC-passing, and GEX-reference cells.")
    }
    
    message(length(common_cells), " cells retained after QC + GEX filtering for sample ", sample.id)
    
    s.object <- s.object[, common_cells]
    
    # ----------------------------------------
    # 5. Build fragment object and requantify
    #    against the shared peak universe
    # ----------------------------------------
    frags <- CreateFragmentObject(
      path = fragpath,
      cells = common_cells,
      validate.fragments = TRUE
    )
    
    peak_matrix <- FeatureMatrix(
      fragments = frags,
      features = master_peak_universe,
      cells = common_cells
    )
    
    s.object[["peaks"]] <- CreateChromatinAssay(
      counts = peak_matrix,
      sep = c(":", "-"),
      fragments = frags,
      annotation = annotation,
      min.cells = 0,
      min.features = 0
    )
    
    # ----------------------------------------
    # 6. Add cached QC metadata back to object
    # ----------------------------------------
    sample_qc_keep <- sample_qc[match(colnames(s.object), sample_qc$raw_cell), ]
    
    if (!all(sample_qc_keep$raw_cell == colnames(s.object))) {
      stop("QC metadata could not be aligned to retained cells.")
    }
    
    s.object$TSS.enrichment <- sample_qc_keep$TSS.enrichment
    s.object$nucleosome_signal <- sample_qc_keep$nucleosome_signal
    s.object$percent.mt <- sample_qc_keep$percent.mt
    s.object$sample <- sample.id
    s.object$qc_filter <- filter_label
    
    # ----------------------------------------
    # 7. RNA preprocessing
    # ----------------------------------------
    DefaultAssay(s.object) <- "RNA"
    
    s.object <- NormalizeData(
      object = s.object,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = median(s.object$nCount_RNA)
    )
    
    s.object <- SCTransform(s.object, verbose = FALSE)
    s.object <- RunPCA(s.object, npcs = 20, verbose = FALSE)
    s.object <- RunUMAP(
      s.object,
      reduction = "pca",
      dims = 1:20,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )
    
    # ----------------------------------------
    # 8. ATAC preprocessing
    # ----------------------------------------
    DefaultAssay(s.object) <- "peaks"
    
    s.object <- RunTFIDF(s.object)
    s.object <- FindTopFeatures(s.object, min.cutoff = 10)
    s.object <- RunSVD(s.object)
    s.object <- RunUMAP(
      s.object,
      reduction = "lsi",
      dims = 2:20,
      reduction.name = "umap.atac",
      reduction.key = "atacUMAP_"
    )
    
    # ----------------------------------------
    # 9. Save processed sample object
    # ----------------------------------------
    saveRDS(s.object, file = output_path)
    
    successful_samples <- c(successful_samples, sample.id)
    sample_status <- "success"
    sample_error <- NA_character_
    
    message("Successfully processed: ", sample.id)
    
  }, error = function(e) {
    message("Failed to process sample: ", sample.id)
    message("Error message: ", e$message)
    failed_samples <<- c(failed_samples, sample.id)
    sample_status <<- "failed"
    sample_error <<- e$message
    
  }, finally = {
    end_time <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    sample_timing <<- rbind(
      sample_timing,
      data.frame(
        sample = sample.id,
        status = sample_status,
        elapsed_sec = elapsed_sec,
        elapsed_min = elapsed_sec / 60,
        start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
        end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
        error_message = ifelse(is.na(sample_error), "", sample_error),
        stringsAsFactors = FALSE
      )
    )
    
    write.csv(
      sample_timing,
      file = file.path(results_dir, paste0("sample_processing_times", output_suffix, ".csv")),
      row.names = FALSE
    )
    
    message("Elapsed time for ", sample.id, ": ", round(elapsed_sec / 60, 2), " min")
  })
}

# ----------------------------
# Save run summary
# ----------------------------
writeLines(
  successful_samples,
  con = file.path(results_dir, paste0("successful_samples", output_suffix, ".txt"))
)

writeLines(
  failed_samples,
  con = file.path(results_dir, paste0("failed_samples", output_suffix, ".txt"))
)

message("Successfully processed ", length(successful_samples), " samples.")
if (length(failed_samples) > 0) {
  message("Failed samples: ", paste(failed_samples, collapse = ", "))
}