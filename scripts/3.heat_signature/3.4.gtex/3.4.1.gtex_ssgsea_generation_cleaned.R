# Purpose: Run ssGSEA for the HEAT signature on GTEx adipose bulk RNA-seq
# datasets and export one score table per tissue
# Run after: the HEAT signature GMT has been generated
# Inputs: GTEx adipose GCT files, GTEx subject phenotype table, HEAT signature GMT
# Outputs: one *_ssGSEA.csv file per GTEx adipose tissue
# Must not change: ssGSEA workflow, tissue labels, and per-dataset output structure
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
INPUT_GTEX_SUB_GCT <- here(
  "data", "4.gtex_adipose", "gene_tpm_2022-06-06_v10_adipose_subcutaneous.gct"
  )
INPUT_GTEX_VIS_GCT <- here(
"data", "4.gtex_adipose", "gene_tpm_2022-06-06_v10_adipose_visceral_omentum.gct")

INPUT_GTEX_PHENO <- here(
  "data", "4.gtex_adipose", "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
)
INPUT_HEAT_GMT <- here(
  "output", "3.heat_signature", "data", "heat_signature.gmt"
  )

OUTPUT_DIR <- here("output", "3.heat_signature", "5.gtex", "data")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Build table of datasets to process
# -------------------------------------------------------------------------
argument_df <- tibble(
  input = c(INPUT_GTEX_SUB_GCT, INPUT_GTEX_VIS_GCT),
  tissue = c("subcutaneous", "visceral_omentum")
) %>%
  mutate(
    output = file.path(
      OUTPUT_DIR,
      paste0(tools::file_path_sans_ext(basename(input)), "_ssGSEA.csv")
    )
  )

# -------------------------------------------------------------------------
# Read the HEAT signature GMT and create a GeneSetCollection
# -------------------------------------------------------------------------
gmt_data <- readLines(INPUT_HEAT_GMT)

gmt_parsed <- lapply(gmt_data, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  list(
    name = fields[1],
    description = fields[2],
    genes = fields[-(1:2)]
  )
})

heat_signatures <- do.call(rbind, lapply(gmt_parsed, as.data.frame))

gene_set_df <- heat_signatures %>%
  filter(genes != "") %>%
  group_by(name) %>%
  summarise(gene_list = list(unique(genes)), .groups = "drop")

gene_sets <- lapply(seq_len(nrow(gene_set_df)), function(i) {
  GeneSet(gene_set_df$gene_list[[i]], setName = gene_set_df$name[i])
})

heat_collection <- GeneSetCollection(gene_sets)

# -------------------------------------------------------------------------
# Read GTEx subject phenotype annotations once
# -------------------------------------------------------------------------
annotation_df <- read.delim(INPUT_GTEX_PHENO) %>%
  mutate(subject = sub("(GTEX-)(.*)", "\\2", SUBJID))

# -------------------------------------------------------------------------
# Function to process one GTEx expression matrix
# -------------------------------------------------------------------------
process_gene_expression <- function(file_path, tissue_label, heat_collection, annotation_df) {
  # Load expression matrix
  expr_data <- read.delim(file_path, skip = 2, header = TRUE)

  # Compute total expression and keep the highest-expression row per gene symbol
  expr_data$Total_Expression <- rowSums(expr_data[, -c(1, 2)], na.rm = TRUE)
  expr_data <- expr_data[order(-expr_data$Total_Expression), ]
  expr_data <- expr_data[!duplicated(expr_data$Description), ]
  expr_data <- expr_data[!is.na(expr_data$Description), ]

  # Remove helper column, keep gene symbols as row names, and drop ID columns
  expr_data <- expr_data[, -ncol(expr_data)]
  rownames(expr_data) <- expr_data$Description
  expr_data <- expr_data[, -c(1, 2)]

  # Convert to ExpressionSet and run ssGSEA
  expr_set <- ExpressionSet(assayData = as.matrix(expr_data))
  gsva_result <- gsva(ssgseaParam(expr_set, heat_collection))

  # Convert to data frame and add subject-level metadata
  gsva_df <- data.frame(t(gsva_result@assayData$exprs)) %>%
    rownames_to_column(var = "SUBJID") %>%
    mutate(subject = sub("GTEX\\.([^.]+)\\..*", "\\1", SUBJID)) %>%
    left_join(annotation_df, by = "subject") %>%
    dplyr::select(-SUBJID.x, -SUBJID.y)

  gsva_df$tissue <- tissue_label
  gsva_df
}

# -------------------------------------------------------------------------
# Run ssGSEA for each GTEx adipose dataset
# -------------------------------------------------------------------------
for (i in seq_len(nrow(argument_df))) {
  input <- argument_df$input[i]
  tissue <- argument_df$tissue[i]
  output <- argument_df$output[i]

  results <- process_gene_expression(
    file_path = input,
    tissue_label = tissue,
    heat_collection = heat_collection,
    annotation_df = annotation_df
  )

  final_data <- results %>%
    mutate(SEX = ifelse(SEX == 1, "M", "F"))

  write_csv(final_data, output)
}
