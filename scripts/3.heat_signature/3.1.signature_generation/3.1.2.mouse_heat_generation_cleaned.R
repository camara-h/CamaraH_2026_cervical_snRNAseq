# Purpose: Convert the human HEAT signature to mouse orthologues using gprofiler2
# Run after: the human HEAT signature has been generated and saved
# Inputs: HEAT signature RData file produced by the HEAT signature pipeline
# Outputs: mouse HEAT signature as GMT and CSV files
# Must not change: use of gprofiler2 orthologue mapping and the expected HEAT_signature object name
# Cleanup level: readability, annotation, and path harmonization only

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(gprofiler2)
library(readr)
library(dplyr)
library(here)

# -------------------------------------------------------------------------
# Define input and output paths
# -------------------------------------------------------------------------
INPUT_HEAT_RDATA <- here("output", "3.heat_signature", "data", "heat_signature.Rdata")

OUTPUT_MM_MAP <- here("output", "3.heat_signature", "data", "mouse_gorth_map.csv")
OUTPUT_MM_HEAT_GMT <- here("output", "3.heat_signature", "data", "mouse_heat.gmt")
OUTPUT_MM_HEAT_CSV <- here("output", "3.heat_signature", "data", "mouse_heat.csv")

# -------------------------------------------------------------------------
# Load the human HEAT signature
# -------------------------------------------------------------------------
load(INPUT_HEAT_RDATA)

# -------------------------------------------------------------------------
# Convert human genes to mouse orthologues
# -------------------------------------------------------------------------
mouse_map <- gorth(
  HEAT_signature,
  source_organism = "hsapiens",
  target_organism = "mmusculus",
  filter_na = TRUE
)
write_csv(mouse_map, OUTPUT_MM_MAP)
# -------------------------------------------------------------------------
# Save mouse HEAT signature in GMT format
# -------------------------------------------------------------------------
mouse_gmt <- paste(
  "Mouse_HEAT",
  "Hs_to_Mm_HEAT_converted_with_gprofiler2",
  paste(mouse_map$ortholog_name, collapse = "\t"),
  sep = "\t"
)

writeLines(mouse_gmt, OUTPUT_MM_HEAT_GMT)

# -------------------------------------------------------------------------
# Save mouse HEAT signature in CSV format
# -------------------------------------------------------------------------
mouse_heat_df <- mouse_map |>
  select(ortholog_name)

write_csv(mouse_heat_df, OUTPUT_MM_HEAT_CSV)

# -----------------------------
#  Save summary 
# -----------------------------
orthologue_summary <- tibble::tibble(
  metric = c(
    "n_human_heat_genes_input",
    "n_human_genes_with_mouse_orthologue",
    "n_mouse_orthologue_rows",
    "n_unique_mouse_orthologue_genes"
  ),
  value = c(
    length(HEAT_signature),
    dplyr::n_distinct(mouse_map$input),
    nrow(mouse_map),
    dplyr::n_distinct(mouse_map$ortholog_name)
  )
)

readr::write_csv(
  orthologue_summary,
  here("output", "3.heat_signature", "data", "mouse_heat_orthologue_summary.csv")
)

writeLines(
  capture.output(sessionInfo()),
  here("output", "3.heat_signature", "data", "sessionInfo_mouse_heat_generation.txt")
)
