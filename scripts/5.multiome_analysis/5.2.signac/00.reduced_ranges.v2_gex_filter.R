#!/usr/bin/env Rscript

# =========================================================
# Build a shared reduced peak universe across samples.
# Optionally restrict to cells retained in the GEX reference.
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
  make_option("--TSS_cutoff", type = "integer", default = 2,
              help = "TSS enrichment cutoff [default %default]"),
  make_option("--mt_cutoff", type = "numeric", default = 0.1,
              help = "Mitochondrial percentage cutoff [default %default]"),
  make_option("--UMI_cutoff", type = "integer", default = 500,
              help = "Minimum UMI cutoff [default %default]"),
  make_option("--NSS_cutoff", type = "integer", default = 2,
              help = "Nucleosome signal cutoff [default %default]"),
  make_option("--filter_gex", type = "character", default = "TRUE",
              help = "Whether to restrict to cells present in the GEX reference. Use TRUE or FALSE [default %default]"),
  make_option("--data_dir", type = "character",
              default = "/home/camarah/storage/Tseng_Neck_Multiome_2026/03_cellranger/ref_2020/cellranger_arc/",
              help = "Folder containing Cell Ranger ARC outputs [default %default]"),
  make_option("--output_dir", type = "character",
              default = here("../04_outputs"),
              help = "Base output folder [default %default]"),
  make_option("--seurat_ref_path", type = "character",
              default = here("../04_outputs/reference_gex/seurat_ref_meta.rds"),
              help = "Path to reference GEX metadata RDS with cell_canonical column [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------
# Store parsed arguments
# ----------------------------
n_counts_cutoff <- opt$n_counts_cutoff
TSS_cutoff <- opt$TSS_cutoff
mt_cutoff <- opt$mt_cutoff
UMI_cutoff <- opt$UMI_cutoff
nucleosome_cutoff <- opt$NSS_cutoff
data_dir <- opt$data_dir
output_dir <- opt$output_dir
seurat_ref_path <- opt$seurat_ref_path

filter_gex <- toupper(opt$filter_gex) == "TRUE"
output_suffix <- if (filter_gex) "_gex_filtered" else ""

filter_label <- paste0(
  n_counts_cutoff, "_atac_",
  TSS_cutoff, "_TSS_",
  mt_cutoff, "_mt_",
  UMI_cutoff, "_UMI",
  output_suffix
)

# ----------------------------
# Load gene annotation
# ----------------------------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

old_seqlevels <- seqlevels(annotation)
new_seqlevels <- paste0("chr", gsub("^(chr)+", "", old_seqlevels))
new_seqlevels[grepl("^M(T)?$", old_seqlevels)] <- "chrM"
seqlevels(annotation) <- new_seqlevels

message("Annotation loaded.")
message("Example seqlevels: ", paste(head(seqlevels(annotation)), collapse = ", "))
message("Seqlevel style: ", paste(seqlevelsStyle(annotation), collapse = ", "))

# ----------------------------
# Load GEX reference metadata if requested
# ----------------------------
seurat_ref <- NULL
gex_cells <- NULL

if (filter_gex) {
  seurat_ref <- readRDS(seurat_ref_path)
  
  if (!"cell_canonical" %in% colnames(seurat_ref@meta.data)) {
    stop("The reference object does not contain a 'cell_canonical' column in meta.data.")
  }
  
  gex_cells <- unique(seurat_ref$cell_canonical)
  message("Loaded GEX reference with ", length(gex_cells), " canonical cell IDs.")
}

# ----------------------------
# Initialize containers
# ----------------------------
seurat_low_filter_files <- list()
successful_samples <- character()
failed_samples <- character()

samples <- list.dirs(path = data_dir, full.names = FALSE, recursive = FALSE)

# ----------------------------
# Process each sample
# ----------------------------
for (sample.id in samples) {
  
  message("Processing sample: ", sample.id)
  
  tryCatch({
    
    h5_path <- file.path(data_dir, sample.id, "outs", "filtered_feature_bc_matrix.h5")
    fragpath <- file.path(data_dir, sample.id, "outs", "atac_fragments.tsv.gz")
    
    counts <- Read10X_h5(h5_path)
    
    chrom_assay <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = fragpath,
      annotation = annotation,
      min.cells = 10,
      min.features = 200
    )
    
    s.object <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "peaks"
    )
    
    # Add sample prefix to cell names so they match canonical reference IDs
    s.object <- RenameCells(s.object, add.cell.id = sample.id)
    
    # Optionally restrict to cells retained in the GEX reference
    if (filter_gex) {
      cells_keep_this_sample <- intersect(colnames(s.object), gex_cells)
      s.object <- subset(s.object, cells = cells_keep_this_sample)
      
      if (ncol(s.object) == 0) {
        stop("No GEX-kept cells remain for sample.")
      }
    }
    
    DefaultAssay(s.object) <- "peaks"
    
    # Compute basic ATAC QC metrics for tracking
    s.object <- NucleosomeSignal(s.object)
    s.object <- TSSEnrichment(s.object)
    
    s.object$nucleosome_group <- ifelse(
      s.object$nucleosome_signal > nucleosome_cutoff,
      "NS+",
      "NS-"
    )
    
    s.object$sample <- sample.id
    
    seurat_low_filter_files[[sample.id]] <- s.object
    successful_samples <- c(successful_samples, sample.id)
    
    message("Successfully processed: ", sample.id, " with ", ncol(s.object), " cells.")
    
  }, error = function(e) {
    message("Failed to process sample: ", sample.id)
    message("Error message: ", e$message)
    failed_samples <<- c(failed_samples, sample.id)
  })
}

# ----------------------------
# Print run summary
# ----------------------------
message("Successfully processed ", length(successful_samples), " samples.")
if (length(successful_samples) > 0) {
  message("Successful samples: ", paste(successful_samples, collapse = ", "))
}

if (length(failed_samples) > 0) {
  message("Failed samples: ", paste(failed_samples, collapse = ", "))
}

if (length(seurat_low_filter_files) == 0) {
  stop("No valid samples were processed. Cannot build reduced ranges.")
}

# ----------------------------
# Build reduced peak universe
# ----------------------------
reduced_ranges <- lapply(seurat_low_filter_files, granges) %>%
  Reduce(c, .) %>%
  GenomicRanges::reduce()

peakwidths <- width(reduced_ranges)
reduced_ranges <- reduced_ranges[peakwidths < 10000 & peakwidths > 20]

message("Reduced ranges created with ", length(reduced_ranges), " peaks.")

# ----------------------------
# Save outputs
# ----------------------------
master_dir <- file.path(output_dir, paste0("master_peak_universe", output_suffix))
dir.create(master_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(reduced_ranges, file = file.path(master_dir, paste0("reduced_ranges", output_suffix, ".rds")))
writeLines(successful_samples, con = file.path(master_dir, paste0("successful_samples", output_suffix, ".txt")))
writeLines(failed_samples, con = file.path(master_dir, paste0("failed_samples", output_suffix, ".txt")))

message("Outputs saved to: ", master_dir)


