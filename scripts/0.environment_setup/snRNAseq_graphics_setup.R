# Purpose: Define plotting themes, color palettes, and utility functions
# Run after: none. Source at the start of the project
# Inputs: none
# Outputs: objects and functions in the R environment
# Must not change: palette values, object names, and default theme settings
# Cleanup level: annotation and readability only, with no intended effect on outputs

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)

# -------------------------------------------------------------------------
# Versioning helper
# -------------------------------------------------------------------------
# Used for date-stamped outputs when needed downstream
today <- format(Sys.Date(), "%Y_%m_%d")

# -------------------------------------------------------------------------
# ggplot themes
# -------------------------------------------------------------------------

# Theme used for publication-style figures inspired by Nature Metabolism
theme_nature_metabolism <- function(base_size = 6, base_family = "Helvetica") {
  theme_pubr(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Base text
      text = element_text(
        size = base_size,
        family = base_family,
        color = "black"
      ),
      
      # Axes
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      
      # Titles
      plot.title = element_text(
        size = base_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        hjust = 0.5
      ),
      
      # Panel and grid
      panel.border = element_blank(),
      panel.grid = element_blank(),
      
      # Facets
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      
      # Legend
      legend.position = "right",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      
      # Plot margins
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

# Theme for plots where axes and related elements should be removed
theme_no_axis <- function() {
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank()
  )
}

# -------------------------------------------------------------------------
# Color palettes
# -------------------------------------------------------------------------

# Neck region palettes
neck.color_dark <- c(
  "Deep" = "#472713",
  "Intermediate" = "#f2ac0a",
  "Superficial" = "#ffd747"
)

neck.color_light <- c(
  "Deep" = "#7A5A46",
  "Intermediate" = "#f4bc3b",
  "Superficial" = "#ffdf6b"
)

neck.color <- c(
  "Deep" = "#bb2a35",
  "Intermediate" = "#157801",
  "Superficial" = "#01078f"
)

# Gene type palette
gene_type_palette <- c(
  "ADIPOQ" = "#ffe65a",
  "PPARG" = "#ffe65a",
  "PPARGC1A" = "#97600f",
  "PDGFRA" = "#629372",
  "DCN" = "#0098ae",
  "JAM2" = "#946b9f",
  "SHANK3" = "#b37f96",
  "PROX1" = "#946b9f",
  "STEAP4" = "#e0641c",
  "MYOCD" = "#fbb294",
  "CD163" = "#0461cf",
  "CPA3" = "#054fb9",
  "IL7R" = "#8babf1",
  "MS4A1" = "#b3c7f7",
  "CSF3R" = "#0073e6",
  "CDH19" = "#bc1523",
  "NRXN1" = "#d61f33",
  "PTH" = "#f02943",
  "PRM1" = "#fd5399"
)

# Broad cell group palette
group_palette <- c(
  "PreAdipocyte" = "#a4c76d",
  "ASPC" = "#8bbe9b",
  "Mesenchymal" = "#0098ae",
  "Endothelial_venous" = "#b37f96",
  "Endothelial_capillary" = "#dfa9c1",
  "Endothelial_arterial" = "#f7dae6",
  "Lymphatic_endothelial" = "#cc9fd7",
  "Endothelial" = "#dfa9c1",
  "Endothelial Cells" = "#dfa9c1",
  "T_Lymphocyte" = "#8babf1",
  "B_Lymphocyte" = "#b3c7f7",
  "Macrophage" = "#0461cf",
  "Mast" = "#054fb9",
  "Neutrophil" = "#0073e6",
  "Macrophage_ADIPOQ_pos" = "#001965",
  "PRM1_pos" = "#fd5399",
  "Pericyte" = "#e0641c",
  "Smooth_muscle" = "#fbb294",
  "Smooth_muscle_ADIPOQ_pos" = "#f87d27",
  "Neuron_associated_1" = "#d61f33",
  "Neuron_associated_2" = "#f02943",
  "Schwann" = "#f79399",
  "Unknown" = "#ede7f6",
  "White_adipocyte" = "#ffe65a",
  "Brown_adipocyte" = "#97600f"
)

# Detailed cell type palette
cell_type_palette <- c(
  "Mesenchymal" = "#0098ae",
  "ASPC_1" = "#629372",
  "ASPC_2" = "#76a986",
  "ASPC_3" = "#8bbe9b",
  "ASPC_4" = "#b7d8c1",
  "ASPC_5" = "#cce5d4",
  "PreAdipocyte" = "#a4c76d",
  "Endothelial_venous" = "#b37f96",
  "Endothelial_capillary" = "#dfa9c1",
  "Endothelial_arterial" = "#f7dae6",
  "Lymphatic_endothelial_1" = "#946b9f",
  "Lymphatic_endothelial_2" = "#cc9fd7",
  "Lymphatic_endothelial_3" = "#ecd7f1",
  "T_Lymphocyte" = "#8babf1",
  "B_Lymphocyte" = "#b3c7f7",
  "Macrophage" = "#0461cf",
  "Mast" = "#054fb9",
  "Neutrophil" = "#0073e6",
  "Macrophage_ADIPOQ_pos" = "#001965",
  "PRM1_pos" = "#fd5399",
  "Pericyte" = "#e0641c",
  "Smooth_muscle_1" = "#fbb294",
  "Smooth_muscle_2" = "#fdcbb8",
  "Smooth_muscle_ADIPOQ_pos" = "#f87d27",
  "Neuron_associated_1" = "#d61f33",
  "Neuron_associated_2" = "#f02943",
  "Schwann_myelinating" = "#bc1523",
  "Schwann_non_myelinating" = "#f79399",
  "Unknown_1" = "#d4cedc",
  "Unknown_2" = "#e0dae9",
  "Unknown_3" = "#ede7f6",
  "Unknown_4" = "#f4effa",
  "Unknown_5" = "#f7f2fd",
  "White_adipocyte_1" = "#ffe65a",
  "White_adipocyte_2" = "#fff0a6",
  "Brown_adipocyte" = "#97600f"
)

# Human-readable cell type labels
cell_type_palette2 <- c(
  "ADIPOQ+ Smooth Muscle Cells" = "#f87d27",
  "PRM1+ Cells" = "#fd5399",
  "B Lymphocytes" = "#b3c7f7",
  "Schwann Cells" = "#bc1523",
  "Mast Cells" = "#054fb9",
  "Brown Adipocytes" = "#97600f",
  "Mesenchymal Stem Cells" = "#0098ae",
  "Neutrophils" = "#0073e6",
  "Parathyroid Associated Cells" = "#f02943",
  "Lipid Associated Macrophages" = "#001965",
  "Pre-adipocytes" = "#a4c76d",
  "Pericytes" = "#e0641c",
  "T Lymphocytes" = "#8babf1",
  "Lymphatic Endothelial Cells" = "#946b9f",
  "Arterial Endothelial Cells" = "#f7dae6",
  "Capillary Endothelial Cells" = "#dfa9c1",
  "Venous Endothelial Cells" = "#b37f96",
  "Smooth Muscle Cells" = "#fbb294",
  "Neuron Associated Cells" = "#d61f33",
  "Macrophages" = "#0461cf",
  "Adipocyte Progenitor Cells" = "#629372",
  "White Adipocytes" = "#ffe65a"
)

# Short labels for figures with limited space
cell_type_short_pallette <- c(
  "PRM1+" = "#fd5399",
  "B Cells" = "#b3c7f7",
  "SchwannCs" = "#bc1523",
  "MastCs" = "#054fb9",
  "BAds" = "#97600f",
  "MSCs" = "#0098ae",
  "Neutrophils" = "#0073e6",
  "ParaThyCs" = "#f02943",
  "PreAds" = "#a4c76d",
  "Pericytes" = "#e0641c",
  "T Cells" = "#8babf1",
  "LECs" = "#946b9f",
  "EndoACs" = "#f7dae6",
  "EndoSCs" = "#dfa9c1",
  "EndoVCs" = "#b37f96",
  "EndoCs" = "#dfa9c1",
  "SMCs" = "#fbb294",
  "NeuroCs" = "#d61f33",
  "Macro" = "#0461cf",
  "ASPC" = "#629372",
  "WAds" = "#ffe65a"
)

# Combined palette object used across plots
palette.use <- c(
  group_palette,
  cell_type_palette,
  cell_type_palette2,
  cell_type_short_pallette
)

# Palette for Angueira dataset or related comparisons
anguiera_palette <- c(
  "Adipocytes" = "#ffe65a",
  "Immune Cells" = "#0461cf",
  "MFAP5+ Fibroblasts" = "#06ddfa",
  "PPARG+ Fibroblasts" = "#a4c76d",
  "SEMA5A+ Cells" = "#ff5b53",
  "Transitional Cells" = "#7b9292",
  "SMC" = "#fbb294",
  "Endothelial Cells" = "#dfa9c1",
  "PROX1+/NALCN+ Cells" = "#f74d97",
  "PPARG+ SMC Like Cells" = "#e0641c",
  "Mesothelial Cells" = "#502864",
  
  "Adipocytes_0" = "#fff4c1",
  "Adipocytes_1" = "#fff0a6",
  "Adipocytes_2" = "#ffe65a",
  "Adipocytes_3" = "#e8c846",
  "Adipocytes_4" = "#f6e39e",
  "Adipocytes_5" = "#97600f"
)

# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------

# Convert millimeters to inches, useful for figure sizing
mm_to_in <- function(mm) {
  mm / 25.4
}