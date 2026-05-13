#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(optparse)
library(dplyr)
library(purrr)
library(here)

option_list <- list(
  make_option("--data_dir", type = "character",
              default = "/home/camarah/storage/Tseng_Neck_Multiome_2026/03_cellranger/ref_2020/cellranger_arc/",
              help = "Folder with cellranger outputs"),
  make_option("--output_dir", type = "character",
              default = here("../04_outputs"),
              help = "Output folder"),
  make_option("--qc_dirname", type = "character",
              default = "qc_cache",
              help = "Subfolder to store cached QC metrics")
)

opt <- parse_args(OptionParser(option_list = option_list))

data_dir <- opt$data_dir
output_dir <- opt$output_dir
qc_dirname <- opt$qc_dirname

qc_dir <- file.path(output_dir, qc_dirname)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

old <- seqlevels(annotation)
new <- paste0("chr", old)
new[old == "MT"] <- "chrM"
seqlevels(annotation) <- new

samples <- list.dirs(path = data_dir, full.names = FALSE, recursive = FALSE)

qc_list <- list()
successful_samples <- character()
failed_samples <- character()

for (sample.id in samples) {
  message("Computing QC for sample: ", sample.id)
  
  tryCatch({
    h5_path <- file.path(data_dir, sample.id, "outs/filtered_feature_bc_matrix.h5")
    fragpath <- file.path(data_dir, sample.id, "outs/atac_fragments.tsv.gz")
    
    counts <- Read10X_h5(h5_path)
    
    chrom_assay <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = fragpath,
      annotation = annotation,
      min.cells = 1,
      min.features = 1
    )
    
    s.object <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA"
    )
    
    common_cells <- intersect(colnames(s.object), colnames(chrom_assay))
    s.object <- s.object[, common_cells]
    chrom_assay <- subset(chrom_assay, cells = common_cells)
    
    s.object[["peaks"]] <- chrom_assay
    
    DefaultAssay(s.object) <- "peaks"
    s.object <- NucleosomeSignal(s.object)
    s.object <- TSSEnrichment(s.object, fast = FALSE)
    
    DefaultAssay(s.object) <- "RNA"
    s.object$percent.mt <- PercentageFeatureSet(s.object, pattern = "^MT-")
    s.object$sample <- sample.id
    s.object$cell <- paste(sample.id, colnames(s.object), sep = "_")
    
    qc_tbl <- s.object@meta.data %>%
      tibble::rownames_to_column("raw_cell") %>%
      mutate(
        sample = sample.id,
        cell = paste(sample, raw_cell, sep = "_")
      ) %>%
      select(
        cell, raw_cell, sample,
        nCount_RNA, nFeature_RNA,
        nCount_peaks, nFeature_peaks,
        percent.mt, TSS.enrichment, nucleosome_signal
      )
    
    qc_list[[sample.id]] <- qc_tbl
    successful_samples <- c(successful_samples, sample.id)
    
  }, error = function(e) {
    message("Failed QC for sample: ", sample.id)
    message("Error message: ", e$message)
    failed_samples <<- c(failed_samples, sample.id)
  })
}

qc_cache <- bind_rows(qc_list)

saveRDS(qc_cache, file = file.path(qc_dir, "qc_cache.rds"))
write.csv(qc_cache, file = file.path(qc_dir, "qc_cache.csv"), row.names = FALSE)

writeLines(successful_samples, con = file.path(qc_dir, "successful_samples.txt"))
writeLines(failed_samples, con = file.path(qc_dir, "failed_samples.txt"))

message("Saved QC cache to: ", file.path(qc_dir, "qc_cache.rds"))
message("Total cells in QC cache: ", nrow(qc_cache))