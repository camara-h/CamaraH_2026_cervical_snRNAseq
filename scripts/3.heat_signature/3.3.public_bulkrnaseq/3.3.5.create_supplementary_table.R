# ============================================================
# Create one supplementary reproducibility workbook per dataset
# ssGSEA HEAT vs UCP1 only
# ============================================================

library(tidyverse)
library(openxlsx)
library(janitor)

# Set this to the parent folder containing:
# Cleaned_Datasets/
# Experimental_Design/
# ssGSEA/
base_dir <- here("output/3.heat_signature/3.public_bulkrnaseq/data/")

out_dir <- file.path(base_dir, "Supplementary_Reproducibility_ByDataset")
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

read_csv_clean <- function(file) {
  readr::read_csv(file, show_col_types = FALSE) %>%
    janitor::clean_names()
}

safe_sheet_name <- function(x) {
  x %>%
    stringr::str_replace_all("[\\[\\]\\*\\?/\\\\:]", "_") %>%
    stringr::str_sub(1, 31)
}

add_sheet <- function(wb, sheet_name, data) {
  sheet_name <- safe_sheet_name(sheet_name)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeDataTable(wb, sheet = sheet_name, x = data)
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(data), widths = "auto")
}

# Extract dataset names from files
# This assumes names like:
# Castella_Cleaned.csv
# Din_SampleGroups.csv
# Gavrila_HEATvUCP1.csv

get_dataset_from_cleaned <- function(x) {
  basename(x) %>%
    str_remove("\\.csv$") %>%
    str_remove("_Cleaned$") %>%
    str_remove(" - Clean$")
}

get_dataset_from_meta <- function(x) {
  basename(x) %>%
    str_remove("\\.csv$") %>%
    str_remove("_SampleGroups$")
}

get_dataset_from_ssgsea <- function(x) {
  basename(x) %>%
    str_remove("\\.csv$") %>%
    str_remove("_HEATvUCP1$")
}

# ------------------------------------------------------------
# Locate files
# ------------------------------------------------------------

cleaned_files <- list.files(
  file.path(base_dir, "Cleaned_Datasets"),
  pattern = "\\.csv$",
  full.names = TRUE
)

meta_files <- list.files(
  file.path(base_dir, "Experimental_Design"),
  pattern = "\\.csv$",
  full.names = TRUE
)

ssgsea_ucp1_files <- list.files(
  file.path(base_dir, "ssGSEA"),
  pattern = "_HEATvUCP1\\.csv$",
  full.names = TRUE
)

missing_file <- file.path(base_dir, "ssGSEA", "missing_heat_gsea_bulk.csv")

# Make lookup tables
cleaned_lookup <- tibble(
  dataset = get_dataset_from_cleaned(cleaned_files),
  cleaned_file = cleaned_files
)

meta_lookup <- tibble(
  dataset = get_dataset_from_meta(meta_files),
  meta_file = meta_files
)

ssgsea_lookup <- tibble(
  dataset = get_dataset_from_ssgsea(ssgsea_ucp1_files),
  ssgsea_file = ssgsea_ucp1_files
)

all_datasets <- union(cleaned_lookup$dataset, meta_lookup$dataset) %>%
  union(ssgsea_lookup$dataset) %>%
  sort()
all_datasets <- all_datasets[which(all_datasets != "Giroud")]

# Read missing genes table once
missing_heat <- read_csv_clean(missing_file)

# ------------------------------------------------------------
# Create one workbook per dataset
# ------------------------------------------------------------

# Start index to match manuscript Supp Table order
i <- 13
for (ds in all_datasets) {
  
  
  message("Creating workbook for: ", ds)
  
  wb <- openxlsx::createWorkbook()
  
  cleaned_file <- cleaned_lookup %>%
    filter(dataset == ds) %>%
    pull(cleaned_file)
  
  meta_file <- meta_lookup %>%
    filter(dataset == ds) %>%
    pull(meta_file)
  
  ssgsea_file <- ssgsea_lookup %>%
    filter(dataset == ds) %>%
    pull(ssgsea_file)
  
  readme <- tibble(
    sheet = c(
      "cleaned_dataset",
      "sample_groups",
      "ssGSEA_HEATvUCP1",
      "missing_heat_genes"
    ),
    description = c(
      "Cleaned expression or analysis-ready dataset used for this dataset.",
      "Sample group metadata / experimental design used for this dataset.",
      "ssGSEA HEAT vs UCP1 results for this dataset.",
      "Genes missing from HEAT/GSEA inputs. Filtered by dataset when dataset column is available; otherwise full missing-gene summary is shown."
    ),
    source_file = c(
      ifelse(length(cleaned_file) == 1, basename(cleaned_file), NA_character_),
      ifelse(length(meta_file) == 1, basename(meta_file), NA_character_),
      ifelse(length(ssgsea_file) == 1, basename(ssgsea_file), NA_character_),
      basename(missing_file)
    )
  )
  
  add_sheet(wb, "README", readme)
  
  if (length(cleaned_file) == 1) {
    cleaned_data <- read_csv_clean(cleaned_file)
    add_sheet(wb, "cleaned_dataset", cleaned_data)
  }
  
  if (length(meta_file) == 1) {
    meta_data <- read_csv_clean(meta_file)
    add_sheet(wb, "sample_groups", meta_data)
  }
  
  if (length(ssgsea_file) == 1) {
    ssgsea_data <- read_csv_clean(ssgsea_file)
    add_sheet(wb, "ssGSEA_HEATvUCP1", ssgsea_data)
  }
  
  # Try to filter missing gene table by dataset if there is a dataset-like column
  missing_for_ds <- missing_heat
  
  dataset_cols <- names(missing_heat)[names(missing_heat) %in% c("dataset", "study", "cohort", "source")]
  
  if (length(dataset_cols) > 0) {
    dataset_col <- dataset_cols[1]
    
    missing_for_ds <- missing_heat %>%
      filter(str_detect(
        string = as.character(.data[[dataset_col]]),
        pattern = regex(ds, ignore_case = TRUE)
      ))
  }
  
  add_sheet(wb, "missing_heat_genes", missing_for_ds)
  
  output_file <- file.path(
    out_dir,
    paste0("Supplementary_Table_",i, "_", ds, "_clean_bulk_ssGSEA_HEATvUCP1.xlsx")
  )
  
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Add index
  i <- i + 1
}