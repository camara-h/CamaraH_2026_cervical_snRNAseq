# Purpose: Read public bulk RNA-seq datasets, standardize gene identifiers,
# and write cleaned expression matrices for downstream ssGSEA analysis
# Run after: raw public bulk datasets are available
# Inputs:
#   - single-file mode: one dataset via --input and --output
#   - manifest mode: a CSV with dataset, input, and output columns
# Outputs: cleaned expression matrix/matrices with gene symbols as row names
# Must not change: Ensembl-to-symbol mapping logic, duplicate removal rule,
# supported input formats, and final cleaned matrix structure
# Cleanup level: readability, annotation, and improved batch processing

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(biomaRt)
library(tidyverse)
library(optparse)
library(here)
library(openxlsx)

# -------------------------------------------------------------------------
# Define command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option("--input", type = "character", help = "Path to one RNA-seq dataset."),
  make_option("--output", type = "character", help = "Path to one cleaned dataset."),
  make_option(
    "--manifest",
    type = "character",
    help = "Path to CSV manifest with columns: dataset, input."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -------------------------------------------------------------------------
# Define default manifest path
# -------------------------------------------------------------------------
default_manifest <- here(
  "data",
  "5.published_datasets",
  "Public_BulkRNAseq",
  "public_bulk_manifest.csv"
)

manifest_path <- opt$manifest

# -----------------------------
#  Create OUTPUT path 
# -----------------------------

OUTPUT_DIR <- here("output", "3.heat_signature", "3.public_bulkrnaseq", "data", "Cleaned_Datasets")

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------
read_expression_file <- function(input) {
  file_ext <- tolower(tools::file_ext(input))
  
  df <- if (file_ext %in% c("tsv", "txt")) {
    read.delim(file = input, check.names = FALSE)
  } else if (file_ext == "csv") {
    read.csv(file = input, check.names = FALSE)
  } else if (file_ext %in% c("xlsx", "xls")) {
    openxlsx::read.xlsx(xlsxFile = input, rowNames = TRUE)
  } else {
    stop("Unsupported file format: ", input)
  }
  
  if (ncol(df) == 0) {
    stop("Input file has no columns.")
  }
  
  first_col <- df[[1]]
  
  if (!is.numeric(first_col)) {
    rownames(df) <- make.unique(as.character(first_col))
    df[[1]] <- NULL
  }
  
  df
}

clean_expression_matrix <- function(object, mart) {
  ensembl_ids <- rownames(object)
  is_ensembl <- all(grepl("^ENSG[0-9]+", ensembl_ids))

  if (is_ensembl) {
    gene_mapping <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = mart
    )

    colnames(gene_mapping) <- c("Ensembl_ID", "Gene_Symbol")

    object$Ensembl_ID <- rownames(object)
    object <- merge(
      gene_mapping,
      object,
      by.x = "Ensembl_ID",
      by.y = "Ensembl_ID",
      all.y = TRUE
    ) %>%
      dplyr::filter(Gene_Symbol != "") %>%
      dplyr::filter(!is.na(Gene_Symbol))

    object$Total_Expression <- rowSums(object[, -c(1, 2)])
    object <- object[order(-object$Total_Expression), ]
    object <- object[!duplicated(object$Gene_Symbol), ]

    object <- object[, -ncol(object)]
    rownames(object) <- object$Gene_Symbol
    object <- object[, -c(1, 2)]
  } else {
    warning("Row names do not appear to be Ensembl IDs. Assuming they are gene symbols.")
  }

  object
}

process_one_dataset <- function(input, output, mart) {
  object <- read_expression_file(input)
  object <- clean_expression_matrix(object, mart)

  out_dir <- dirname(output)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(object, output)
  message("Saved cleaned dataset to: ", output)
}

# -------------------------------------------------------------------------
# Connect to Ensembl once
# -------------------------------------------------------------------------
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Try the US East mirror
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "asia")

# -------------------------------------------------------------------------
# Run in single-file mode or manifest mode
# -------------------------------------------------------------------------
if (!is.null(opt$input) || !is.null(opt$output)) {
  if (is.null(opt$input) || is.null(opt$output)) {
    stop("Single-file mode requires both --input and --output.")
  }

  process_one_dataset(
    input = opt$input,
    output = opt$output,
    mart = mart
  )

} else {
  if (is.null(manifest_path)) {
    manifest_path <- default_manifest
  }

  if (!file.exists(manifest_path)) {
    stop("Manifest not found: ", manifest_path)
  }

  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

  required_cols <- c("dataset", "input")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  # seq_len(nrow(manifest))
  for (i in 5) {
    process_one_dataset(
      input = manifest$input[i],
      output = file.path(OUTPUT_DIR, paste0(manifest$dataset[i], "_Cleaned.csv")),
      mart = mart
    )
  }
}

# -----------------------------
#  Save summary 
# -----------------------------

cleaning_summary <- tibble::tibble(
  dataset = character(),
  n_genes_input = integer(),
  n_genes_after_symbol_mapping = integer(),
  n_genes_after_duplicate_removal = integer(),
  n_samples = integer(),
  input_file = character(),
  output_file = character()
)