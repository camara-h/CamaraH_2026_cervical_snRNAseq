# Purpose: Run ssGSEA for a signature across bulk RNA-seq datasets and export
# per-sample signature scores together with CKMT2 expression
# Run after: cleaned bulk RNA-seq expression matrices, matching metadata files,
# and the gene signature CSV have been generated
# Inputs: bulk expression matrices, bulk metadata tables, signature gene list
# Outputs: one CSV per dataset containing ssGSEA score and CKMT2 expression
# Must not change: dataset matching logic, ssGSEA workflow, group ordering,
# and output file naming
# Cleanup level: readability, path harmonization, and removal of unused code only

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(GSVA)
library(GSEABase)
library(tidyverse)
library(here)

# -------------------------------------------------------------------------
# Define input and output paths
# -------------------------------------------------------------------------

# Directory containing cleaned bulk RNA-seq expression matrices
BULK_DIR <- here("output", "3.heat_signature", "3.public_bulkrnaseq", "data", "Cleaned_Datasets")

# Directory containing matching metadata tables for each bulk dataset
BULK_META_DIR <- here(
  "output", "3.heat_signature", "3.public_bulkrnaseq", "data", "Experimental_Design")

# Gene list in CSV format, one gene per row
GENE_LIST_CSV <- here(
  "output",
  "3.heat_signature",
  "data",
  "heat_signature.csv"
)

GENE_LIST_NAME <- "HEAT"

# Output directory for public bulk RNA-seq ssGSEA results
OUTPUT_DATA_DIR <- here(
  "output",
  "3.heat_signature",
  "3.public_bulkrnaseq",
  "data",
  "ssGSEA"
)

dir.create(OUTPUT_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Identify input files
# -------------------------------------------------------------------------

# Expression files
exp_files <- list.files(BULK_DIR, full.names = TRUE)
cleaned_files <- grep("_Cleaned", exp_files, value = TRUE)

# Metadata files
meta_files <- list.files(BULK_META_DIR, full.names = TRUE)

# Dataset names inferred from cleaned expression file names
dataset_names <- str_remove(basename(cleaned_files), "_Cleaned.csv")
dataset_names <- grep("Moushart", dataset_names, value = TRUE, invert = TRUE)


# -------------------------------------------------------------------------
# Run ssGSEA for each dataset
# -------------------------------------------------------------------------

# Save all per-dataset results into one supplemental workbook
supplemental_tables <- list()

for (dataset in dataset_names) {
  
  # -----------------------------------------------------------------------
  # Read expression matrix, signature genes, and metadata
  # -----------------------------------------------------------------------
  object <- read.csv(
    str_subset(cleaned_files, dataset)[1],
    header = TRUE,
    row.names = 1
  )
  
  gene_list <- read.csv(
    GENE_LIST_CSV,
    header = TRUE,
    col.names = c("Genes")
  )
  
  # Remove CMKT2 from signature for more reliable correlations
  gene_list <- gene_list[which(gene_list$Genes != "CKMT2"),]
  
  meta <- read.csv(str_subset(meta_files, dataset)[1])
  
  # -----------------------------------------------------------------------
  # Ensure sample groups follow the intended display order
  # -----------------------------------------------------------------------
  meta$group <- factor(
    meta$group,
    levels = unique(meta$group[order(meta$order)])
  )
  
  # -----------------------------------------------------------------------
  # Build the gene set collection and run ssGSEA
  # -----------------------------------------------------------------------
  heat_collection <- GeneSetCollection(
    GeneSet(
      gene_list,
      setName = GENE_LIST_NAME
    )
  )
  
  expr_set <- ExpressionSet(assayData = as.matrix(object))
  
  object.heat <- gsva(
    ssgseaParam(expr_set, heat_collection)
  )
  
  object.heat <- data.frame(t(object.heat@assayData$exprs)) %>%
    rownames_to_column(var = "sample") %>%
    left_join(meta, by = "sample")
  
  # -----------------------------------------------------------------------
  # Join ssGSEA scores with CKMT2 expression
  # -----------------------------------------------------------------------
  ckmt2vHEAT <- object.heat %>%
    left_join(
      object %>%
        dplyr::filter(rownames(object) %in% c("CKMT2")) %>%
        t() %>%
        data.frame() %>%
        rownames_to_column(var = "sample"),
      by = "sample"
    ) %>%
    rename(` ` = group)
  
  
  # -----------------------------
  #  Store data for supplemental table 
  # -----------------------------
  
  supplemental_tables[[dataset]] <- ckmt2vHEAT
  
  
  # -----------------------------------------------------------------------
  # Save results for the current dataset
  # -----------------------------------------------------------------------
  output_file <- paste0(dataset, "_", GENE_LIST_NAME, "vCKMT2.csv")
  output <- file.path(OUTPUT_DATA_DIR, output_file)
  
  write.csv(ckmt2vHEAT, output, row.names = FALSE)
  
  
  # -----------------------------
  #  Write summary 
  # -----------------------------
  if(exists("gene_overlap_summary")){
    gene_overlap_summary <- rbind(gene_overlap_summary,
                                  tibble::tibble(
                                    dataset = dataset,
                                    n_heat_genes = length(gene_list),
                                    n_heat_genes_detected = sum(gene_list %in% rownames(object)),
                                    missing_heat_genes = paste(setdiff(gene_list, rownames(object)), collapse = ";")
                                  )
    )
  }else{
    gene_overlap_summary <- tibble::tibble(
      dataset = dataset,
      n_heat_genes = length(gene_list),
      n_heat_genes_detected = sum(gene_list %in% rownames(object)),
      missing_heat_genes = paste(setdiff(gene_list, rownames(object)), collapse = ";")
    )
  }
}

# -----------------------------
#  Save summary 
# -----------------------------
write.csv(gene_overlap_summary, file.path(OUTPUT_DATA_DIR, "missing_heat_gsea_bulk.csv"), row.names = FALSE)

write.xlsx(
  supplemental_tables,
  file = file.path(OUTPUT_DATA_DIR, paste0(GENE_LIST_NAME, "_vCKMT2_summary.xlsx")),
  rowNames = FALSE,
  overwrite = TRUE
)
