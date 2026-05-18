# ---
# title: "5-S4-FigurePlotting"
# author: "Henrique"
# date: "2026-05-07"
# output: html_document
# ---

# # --- Setup ---
# ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

mm_to_in <- function(mm) mm / 25.4

# Load libraries
library(Seurat)
library(ggplot2)
library(VISION)
library(dplyr)
library(tidyr)
library(CellChat)
library(circlize) ### A tutorial on how to create your own circlize can be found here: https://r-graph-gallery.com/311-add-labels-to-hierarchical-edge-bundling.html
library(ComplexHeatmap)
library(stringr)
library(viridis)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(conflicted)
library(here)
conflicts_prefer(ggplot2::annotate)
conflicts_prefer(base::unname)
conflicts_prefer(tidygraph::filter)


# Define short cell names
cell_type_mapping <- c(
  "Brown Adipocytes" = "BAds",
  "White Adipocytes" = "WAds",
  "Lipid Associated Macrophages" = "PreAds",
  "ADIPOQ+ Smooth Muscle Cells" = "PreAds",
  "Pre-adipocytes" = "PreAds",
  "Adipocyte Progenitor Cells" = "ASPC",
  "Mesenchymal Stem Cells" = "MSCs",
  "Arterial Endothelial Cells" = "EndoCs",
  "Venous Endothelial Cells" = "EndoCs",
  "Capillary Endothelial Cells" = "EndoCs",
  "Lymphatic Endothelial Cells" = "LECs",
  "Pericytes" = "Pericytes",
  "Smooth Muscle Cells" = "SMCs",
  "Macrophages" = "Macro",
  "Mast Cells" = "MastCs",
  "T Lymphocytes" = "T Cells",
  "B Lymphocytes" = "B Cells",
  "Neutrophils" = "Neutrophils",
  "Schwann Cells" = "SchwannCs",
  "Neuron Associated Cells" = "NeuroCs",
  "Parathyroid Associated Cells" = "ParaThyCs",
  "PRM1+ Cells" = "PRM1+"
)

# OBS:  theme_nature_metabolism() and color palettes set on "scripts/0.environment_setup/snRNAseq_graphics_setup.R"
source(here("scripts/0.environment_setup/snRNAseq_graphics_setup.R"))

# --- Setup paths ---

## Inputs
SEURAT_INPUT <- here(
  "..",
  "NeckAdipose/ocean_code_data/cervical_at_gex_seurat.rds")

CELLCHAT_DEEP <- here(
  "..",
  "NeckAdipose/ocean_code_data/Figure_5_S4/cellchat_deep_pop.size.false_cell_type_short.rds")
CELLCHAT_SUPERFICIAL <- here(
  "..",
  "NeckAdipose/ocean_code_data/Figure_5_S4/cellchat_superficial_pop.size.false_cell_type_short.rds")
CELLCHAT_MERGED <- here(
  "..",
  "NeckAdipose/ocean_code_data/Figure_5_S4/cellchat_merged_pop.size.false_cell_type_short.rds")


MEBOCOST_COMM_CSV <- here(
  "..",
  "NeckAdipose/ocean_code_data/Figure_5_S4/communication_result_multi.csv")

# Outputs
OUTPUT_DIR <- here("..", "results", "Figure_5_S4")

dir.create(OUTPUT_DIR, showWarnings = F, recursive = T)
# ```

# # --- Load Seurat  CamaraH---
# ```{r}
s.object <- readRDS(SEURAT_INPUT)
# ```

# # --- CellChat ---
# ## --- Load the objects and create CC object
# ```{r}
# --- Load the cellchat objects
deep_cellchat <- readRDS(CELLCHAT_DEEP)
superficial_cellchat <- readRDS(CELLCHAT_SUPERFICIAL)
cellchat_merged <- readRDS(CELLCHAT_MERGED)
# ```

# ## --- Create object_list
# ```{r}
object_list <- list(Deep = deep_cellchat, Superficial = superficial_cellchat)
# ```


# ## --- Extract signaling information
# ```{r}
# Extract signaling information
df.sup <- as.data.frame(subsetCommunication(superficial_cellchat, thresh = 1))
df.deep <- as.data.frame(subsetCommunication(deep_cellchat, thresh = 1))

# Merge and filter significant
df.merge <- rbind(
  df.sup |> mutate(depot = "Superficial"),
  df.deep |> mutate(depot = "Deep")
) |>
  dplyr::filter(pval <= 0.05)

annotation_order <- df.merge |>
  count(annotation) |>
  arrange(-n) |>
  pull(annotation)
df.merge <- df.merge |> mutate(annotation = factor(annotation, levels = annotation_order))
# ```

# # --- Panel 5B. CellChat - Number of interactions per depot ---
# ## --- Create plot ---
# ```{r, fig.width= 1.968504, fig.height= 1.417323, warning=FALSE}
# Main Communication types
pd <- position_dodge(width = 0.9)

gg_comm_bar <- ggplot(df.merge, aes(x = annotation, fill = depot)) +
  geom_bar(position = pd) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count), group = depot),
    position = pd,
    vjust = -0.15,
    size = 5 / .pt
  ) +
  labs(
    x = "Signaling Type", y = str_wrap("Number of inferred interactions", 25), fill = "Depot",
    title = "Communication Types"
  ) +
  theme_nature_metabolism() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = neck.color) +
  scale_y_continuous(expand = expansion(c(0.05, 0.20)))
# ```

# ## --- Plotting and save ---
# ```{r, fig.width= 1.968504, fig.height= 1.224409, warning=FALSE}
mm_to_in(c(50, 31.1))
gg_comm_bar <- gg_comm_bar +
  # plot_annotation(tag_levels = list("d")) + theme(plot.tag = element_text(face = "bold", size = 6)) +
  # NoLegend() +
  theme(axis.title.x = element_blank())
gg_comm_bar

ggsave(file.path(OUTPUT_DIR, "5B_barplot_communication_types.pdf"), plot = gg_comm_bar, width = 50, height = 31.1, units = "mm", dpi = 600, bg = "white")
# ```

# # --- Panel S4C. CellChat - Number of interactions per depot ---
# ## --- Create plot ---
# ```{r}
# Main Communication types to Wads
gg_comm_to_wad <- ggplot(df.merge |> dplyr::filter(target == "WAds"), aes(x = annotation, fill = depot)) +
  geom_bar(position = pd) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count), group = depot),
    position = pd,
    vjust = -0.15,
    size = 5 / .pt
  ) +
  labs(x = "Signaling Type", y = str_wrap("Number of inferred interactions", 35), fill = "Depot", title = "Communication Types", subtitle = "White Adipocytes as Target") +
  theme_nature_metabolism() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = neck.color) +
  scale_y_continuous(expand = expansion(c(0.05, 0.15)))
# ```

# ## --- Plotting and saving ---
# ```{r, fig.width= 1.574803, fig.height= 1.771654, warning=FALSE}
mm_to_in(c(40, 45))
gg_comm_to_wad
ggsave(file.path(OUTPUT_DIR, "S4C_barplot_communication_types_to_wads.pdf"), plot = gg_comm_to_wad, width = 40, height = 45, units = "mm", dpi = 600, bg = "white")
# ```


# # --- Panel 5G. Communication flow analysis ---
# ## --- Set up helper function ---
# ```{r, fig.width=1.574803, fig.height=1.417323, warning=FALSE}
mm_to_in(c(40, 36))

# Helper function to plot only one signaling type consistently
customRankNet <- function(cellchat = cellchat_merge, ccht_comm_prob_df = df.merge, comm_type) {
  pathway.show <- ccht_comm_prob_df |>
    dplyr::filter(annotation == comm_type) |>
    pull(pathway_name)
  rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = "WAds", stacked = T, do.stat = TRUE, color.use = neck.color, signaling = pathway.show) + theme_nature_metabolism() + theme(legend.position = "top") +
    labs(title = comm_type, subtitle = "WAds as Receivers") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
# ```

# ## --- Plotting and save ---
# ```{r, fig.width=1.614173, fig.height=1.771654, warning=FALSE}
mm_to_in(c(41, 45))
customRankNet(cellchat_merge, df.merge, comm_type = "Secreted Signaling")
ggsave(file.path(OUTPUT_DIR, paste0("5G_communication_flow_to_wads_", "Secreted Signaling", ".pdf")), width = 41, height = 45, units = "mm", dpi = 600, bg = "white")
# ```
# # --- Panel S4E. Communication flow analysis ---
# ## --- Plotting and save ---
# ```{r, fig.width=1.574803, fig.height=1.574803, warning=FALSE}
mm_to_in(c(40, 40))
customRankNet(cellchat_merge, df.merge, comm_type = "ECM-Receptor")
ggsave(file.path(OUTPUT_DIR, paste0("S4E_communication_flow_to_wads_", "ECM-Receptor", ".pdf")), width = 40, height = 40, units = "mm", dpi = 600, bg = "white")
# ```

# # --- Panel S4G. Communication flow analysis ---
# ## --- Plotting and save ---

# ```{r, fig.width=1.574803, fig.height=1.574803, warning=FALSE}
customRankNet(cellchat_merge, df.merge, comm_type = "Cell-Cell Contact")
ggsave(file.path(OUTPUT_DIR, paste0("S4G_communication_flow_to_wads_", "Cell-Cell Contact", ".pdf")), width = 40, height = 40, units = "mm", dpi = 600, bg = "white")
# ```


# # --- Panels 5C,E and S4A. 2D communication plot ---
# ## --- Create plots ---
# ```{r, warning=FALSE}
# I have confirmed that the interaction strength is simply the sum of the communication probability within all the significant pathwys
# We can thus create our own plot

## Initialize empty list
gg_2d <- list()

# Start loop to plot for each signaling type
for (comm_type in unique(df.merge$annotation)) {
  df.filter <- df.merge |> dplyr::filter(annotation == comm_type)

  # Calculate strength
  ## Output signal strength
  out_strength_df <- df.filter |>
    group_by(source, depot) |>
    summarise(out_prob = sum(prob)) |>
    rename("cell" = "source")

  ## Input signal strength
  in_strength_df <- df.filter |>
    group_by(target, depot) |>
    summarise(in_prob = sum(prob)) |>
    rename("cell" = "target")

  # Communication count per cell (counts rows where cell appears as source or target).
  communication_count_df <- df.filter %>%
    mutate(target = ifelse(source == target, "", target)) %>% # Remove cell if autocrine to count autocrine just once
    select(depot, source, target) %>%
    pivot_longer(cols = c(source, target), values_to = "cell") %>%
    count(depot, cell, name = "comm_count")

  # Merge data frames
  strengt_df <- left_join(out_strength_df, in_strength_df, by = c("cell", "depot"))

  strengt_df <- strengt_df %>%
    left_join(communication_count_df, by = c("depot", "cell")) %>%
    mutate(
      comm_count = coalesce(comm_count, 0L),
      cell_label = cell
    )

  # Map names
  meta_lookup <- s.object@meta.data |>
    dplyr::select(cell_type, cell_type_short) |>
    dplyr::distinct(cell_type, .keep_all = TRUE) |>
    mutate(cell_type_short = as.character(cell_type_short))
  strengt_df <- dplyr::left_join(strengt_df, meta_lookup, by = c("cell" = "cell_type_short"))

  # Filter cells of interest
  if (comm_type == "Secreted Signaling") {
    strengt_df <- strengt_df |> mutate(cell_label = case_when(cell_type %in% c(
      "White Adipocytes",
      "Brown Adipocytes",
      "Lymphatic Endothelial Cells",
      "Schwann Cells"
    ) ~ cell, TRUE ~ NA))
  } else if (comm_type == "ECM-Receptor") {
    strengt_df <- strengt_df |> mutate(cell_label = case_when(cell_type %in% c(
      "White Adipocytes",
      "Brown Adipocytes",
      "Adipocyte Progenitor Cells",
      "Smooth Muscle Cells",
      "Pericytes",
      "Endothelial Cells",
      "Lymphatic Endothelial Cells",
      "Schwann Cells"
    ) ~ cell, TRUE ~ NA))
  } else if (comm_type == "Cell-Cell Contact") {
    strengt_df <- strengt_df |> mutate(cell_label = case_when(cell_type %in% c(
      "White Adipocytes",
      "B Lymphocytes",
      "Macrophages",
      "Mast Cells",
      "Neutrophils",
      "T Lymphocytes",
      "Endothelial Cells",
      "Lymphatic Endothelial Cells",
      "Schwann Cells"
    ) ~ cell, TRUE ~ NA))
  }


  sup_max <- strengt_df %>%
    ungroup() |>
    dplyr::filter(depot == "Superficial") %>%
    dplyr::summarise(
      x = max(out_prob, na.rm = TRUE) * 1.05,
      y = max(in_prob, na.rm = TRUE) * 1.05
    )

  # Plot 2D plot

  gg_2d[[comm_type]] <- ggplot(strengt_df, aes(x = out_prob, y = in_prob)) +
    geom_rect(
      data = sup_max,
      aes(xmin = -Inf, ymin = -Inf, xmax = x * 1.15, ymax = y * 1.15),
      inherit.aes = FALSE,
      fill = "grey90",
      alpha = 0.3,
      color = "black",
      linewidth = 0.2
    ) +
    geom_point(aes(fill = cell, size = comm_count), shape = 21, color = "black", alpha = 0.75) +
    geom_text_repel(
      aes(label = cell_label, color = cell),
      size = 5 / .pt,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      show.legend = FALSE, face = "bold"
    ) +
    facet_wrap(~depot) +
    scale_fill_manual(values = palette.use) +
    scale_color_manual(values = palette.use) +
    scale_size(range = c(1, 2))




  gg_2d[[comm_type]] <- gg_2d[[comm_type]] +
    theme_nature_metabolism() +
    NoLegend() +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.1, color = "grey90"),
      panel.grid.major.y = element_line(linewidth = 0.1, color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.spacing.x = unit(1.8, "mm"),
      axis.line = element_blank()
    ) +
    scale_x_continuous(
      expand = expansion(c(0.05, 0.15)),
      labels = scales::label_number(accuracy = 0.1, trim = TRUE)
    ) +
    scale_y_continuous(expand = expansion(c(0.05, 0.20))) +
    labs(
      title = comm_type,
      subtitle = "All to All",
      y = "Incoming Signal Strength",
      x = "Outgoing Signal Strength"
    )
}

# Rename for easier access
names(gg_2d) <- names(gg_2d) |>
  tolower() |>
  str_replace_all("[ -]", "_")
names(gg_2d)


# ```

# # --- Panels 5C. 2D communication plot - ECM-Receptor---
# ## --- Plotting and save ---
# ```{r, fig.width= 2.165354, fig.height= 1.299213, warning=FALSE}
mm_to_in(c(55, 33))

ggsave(file.path(OUTPUT_DIR, paste0("5C_2d_communication_plots_", "ecm_receptor", ".pdf")), plot = gg_2d[["ecm_receptor"]], width = 55, height = 33, units = "mm", dpi = 600, bg = "white")
# ```

# # --- Panels 5E. 2D communication plot - Secreted Signaling ---
# ## --- Plotting and save ---
# ```{r, fig.width= 2.165354, fig.height= 1.299213, warning=FALSE}
mm_to_in(c(55, 33))

ggsave(file.path(OUTPUT_DIR, paste0("5E_2d_communication_plots_", "secreted_signaling", ".pdf")), plot = gg_2d[["secreted_signaling"]], width = 55, height = 33, units = "mm", dpi = 600, bg = "white")
# ```

# # --- Panels S4A. 2D communication plot - Cell-Cell Contact ---
# ## --- Plotting and save ---
# ```{r, fig.width= 2.165354, fig.height= 1.692913, warning=FALSE}
mm_to_in(c(55, 43))

ggsave(file.path(OUTPUT_DIR, paste0("S4A_2d_communication_plots_", "cell_cell_contact", ".pdf")), plot = gg_2d[["cell_cell_contact"]], width = 55, height = 43, units = "mm", dpi = 600, bg = "white")
# ```

# # --- Panel S4D. CellChat Stacked communications ---
# ```{r, fig.width=1.771654, fig.height=4.429134}
mm_to_in(c(45, 112.5))
# Stacked flow ---
gg_ccht_ranknet <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = "WAds", stacked = F, do.stat = TRUE, color.use = neck.color) +
  theme_nature_metabolism() +
  theme(
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  labs(
    title = "CellChat\nInformation Flow",
    subtitle = "WAds as Receivers"
  )

gg_ccht_ranknet

# Save as png
ggsave(file.path(OUTPUT_DIR, "S4D_cellchat_ranknet.pdf"), gg_ccht_ranknet, width = 45, height = 112.5, units = "mm", bg = "white", dpi = 600)
# ```

# # --- Panels 5D,H and S4F,H. CellChat - Circle Plot of selected communications ---
# ## --- Define color pallete ---
# ```{r}
# Define the cells and filter the color pallete accordingly
cells <- levels(object_list[[1]]@idents)
color.use <- palette.use[cells]
color.use <- color.use[!is.na(color.use)]
color.use
# ```
# ## --- Plotting and save ---
# ```{r}
meta <- s.object@meta.data

celltype_map <- meta %>%
  dplyr::select(cell_type, cell_type_short) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  tibble::deframe()

# Helper function
rename_celltypes <- function(x, map) {
  x_clean <- stringr::str_remove(x, ".*~ ")
  dplyr::coalesce(unname(map[x_clean]), x_clean)
}

# Set up saving IDs
id_maps <- data.frame(
  pathway = c("COLLAGEN", "CD99", "BMP"),
  main_fig_id = c("D", NA, "H"),
  supp_fig_id = c("F", "H", NA)
)

# -----------------------------
#  Plot data
# -----------------------------

for (n in seq_len(length(id_maps))) {
  pathways.show <- id_maps$pathway[n]

  # Setting up a loop to get around issues when a single depot has significant communication
  weight.max <- 0
  for (i in seq_len(length(object_list))) {
    print(i)
    tryCatch(
      {
        weight <- getMaxWeight(object_list[i], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
        if (weight > weight.max) {
          weight.max <- weight
        }
      },
      error = function(e) {
        message(sprintf("Error at index %d: %s", i, e$message))
        # Continue to next iteration
      }
    )
  }

  # Plot in R
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:length(object_list)) {
    tryCatch(
      {
        netVisual_aggregate(object_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object_list)[i]), color.use = color.use)
      },
      error = function(e) {
        print(e$message)
      }
    )
  }

  if (pathways.show %in% c("COLLAGEN", "BMP")) {
  file_name <- paste0("5", id_maps$main_fig_id[n], "_cellchat_", tolower(pathways.show), "_circleplot.pdf")

  # Adjust dimension per plot
  if (pathways.show == "BMP") {
    plot.width <- 40
    plot.heigth <- 45
  } else if (pathways.show == "COLLAGEN") {
    plot.width <- 70
    plot.heigth <- 33
  } else {
    plot.width <- 50
    plot.heigth <- 40
  }

  # Save as pdf
  ## Note: This needs to be further edited using softwares like Inkscape or Illustrator
  pdf(file.path(OUTPUT_DIR, file_name), width = mm_to_in(plot.width * 3), height = mm_to_in(plot.heigth * 3))
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:length(object_list)) {
    tryCatch(
      {
        netVisual_aggregate(object_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object_list)[i]), color.use = color.use, )
      },
      error = function(e) {
        print(e$message)
      }
    )
  }

  dev.off()
  }
  
  if (pathways.show %in% c("COLLAGEN", "CD99")) {
    file_name <- paste0("S4", id_maps$supp_fig_id[n], "_cellchat_", tolower(pathways.show), "_WAd_as_target_circleplot.pdf")
    plot.width <- 50
    plot.heigth <- 40


    pdf(file.path(OUTPUT_DIR, file_name), width = mm_to_in(plot.width * 3), height = mm_to_in(plot.heigth * 3))
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in 1:length(object_list)) {
      tryCatch(
        {
          netVisual_aggregate(object_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object_list)[i]), color.use = color.use, targets.use = "WAds")
        },
        error = function(e) {
          print(e$message)
        }
      )
    }

    dev.off()
  }
}
# ```

# # --- Panel 5I. BMP Ligand-Receptor Violin Plots ---
# ## --- Filter dataset ---
# ```{r}
plot_list <- list()

Idents(s.object) <- "neck_region"
filtered_s.object <- subset(s.object, idents = c("Intermediate"), invert = TRUE)
Idents(filtered_s.object) <- "cell_type"
Idents(s.object) <- "cell_type"
# ```

# ## --- Plotting and save ---
# ```{r, fig.width=2.755906, fig.height=1.574803, warning=FALSE}
mm_to_in(c(70, 40))
gg_ccht_bmp_vln <- VlnPlot(filtered_s.object, features = c("BMP5", "BMP8B", "ACVR1", "ACVR2A", "BMPR2", "BMPR1A"), idents = c("Smooth Muscle Cells", "Schwann Cells", "White Adipocytes"), split.by = "neck_region", flip = T, stack = T, pt.size = 0) +
  scale_fill_manual(values = neck.color) +
  theme_nature_metabolism() +
  theme(axis.title = element_blank()) +
  theme(legend.position = "top") +
  scale_x_discrete(limits = c("Smooth Muscle Cells", "Schwann Cells", "White Adipocytes"))

gg_ccht_bmp_vln
# Save as png
ggsave(file.path(OUTPUT_DIR, "5I_cellchat_bmp_vln.pdf"), gg_ccht_bmp_vln, width = 70, height = 40, units = "mm", bg = "white", dpi = 600)
# ```
# ## --- Record N ---
# ```{r}
vln_plot_n <- gg_ccht_bmp_vln$data |> count(feature, ident, split)
vln_plot_n |>
  distinct(ident, split, n) |>
  arrange(n)
# ```


# # --- MEBOCOST ---
# ##  --- Load communication table ---
# ```{r}
# Load the communication table
mebocost_comm <- read.csv(MEBOCOST_COMM_CSV)

mebocost_comm %>% distinct(Condition, Sender)

# Filter for significant communications to White Adipocytes in Deep and Superficial
mebocost_comm_filtered <- mebocost_comm %>%
  rename("neck_region" = "Condition") %>%
  mutate(
    Sender = str_remove(Sender, "^.*~ "),
    Receiver = str_remove(Receiver, "^.*~ ")
  ) %>%
  filter(
    neck_region != "Intermediate",
    permutation_test_fdr <= 0.05,
    str_detect(Receiver, "White")
  )
# ```

# # --- Panel S4J. Mebocost - DotPlot of metabolites to White Adipocytes ---
# ## --- Create plots ---
# ```{r}
# Annotate sender region and link metabolite and sensor for plotting
graph.data <- mebocost_comm_filtered %>%
  mutate(
    Sender = paste0(Sender, " (", neck_region, ")"),
    Sender = str_replace(Sender, "Neuron Associated Cells", "Neuron-Like Cells"),
    Metabolite_Name = paste(Metabolite_Name, Sensor, sep = "->"),
    comm_score_log2 = log2(Norm_Commu_Score)
  ) |>
  dplyr::filter(str_detect(
    Metabolite_Name,
    paste0(c(
      "Dehydroepiandrosterone"
    ), collapse = "|")
  ))

# unique x-axis categories in order
sender_levels <- levels(factor(graph.data$Sender))

# Dataframe to add background boxes
box_df <- data.frame(
  Sender = sender_levels,
  x_center = seq_along(sender_levels),
  x_min = seq_along(sender_levels) - 0.5,
  x_max = seq_along(sender_levels) + 0.5,
  neck_region = graph.data$neck_region[match(sender_levels, graph.data$Sender)]
)

# Define the changing points by comparing left_shifted and right_shifted x_axist vector
cell_type_from_sender <- sub(" \\(.*", "", sender_levels)
change_points <- which(cell_type_from_sender[-1] != cell_type_from_sender[-length(cell_type_from_sender)]) + 0.5

# Create plot
gg_mebocost_wad_comm_dotplot <- ggplot(
  graph.data,
  aes(x = Sender, y = Metabolite_Name)
) +
  geom_vline(
    xintercept = change_points,
    linetype = "dashed",
    color = "black",
    linewidth = 0.3,
    alpha = 0.6
  ) +
  geom_point(
    aes(fill = comm_score_log2, size = comm_score_log2),
    shape = 21,
    color = "black"
  ) +
  scale_fill_viridis_c(
    option = "H",
    name = "Communication\nScore (log2)",
    guide = guide_legend(
      override.aes = list(shape = 21)
    )
  ) +
  scale_size_continuous(
    name = "Communication\nScore (log2)"
  ) +
  guides(
    fill = guide_legend(
      reverse = TRUE,
      override.aes = list(color = "black")
    ),
    size = guide_legend(
      reverse = TRUE
    )
  ) +
  theme_nature_metabolism() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "MEBOCOST Communication Score",
    subtitle = "White Adipocytes as Receivers",
    y = NULL
  )
# ```

# ## --- Plotting and save ---
# ```{r, fig.width=3.543307, fig.height=1.574803}
mm_to_in(c(90, 40))
gg_mebocost_wad_comm_dotplot <- gg_mebocost_wad_comm_dotplot +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )

gg_mebocost_wad_comm_dotplot

ggsave(filename = file.path(OUTPUT_DIR, "S4J_mebocost_white_adipocytes_comm_dotplot.pdf"), plot = gg_mebocost_wad_comm_dotplot, width = 90, height = 40, units = "mm", dpi = 300)
# ```
# # --- Panel 5K. Mebocost relative communication flow ---
# ## --- Create plot ---
# ```{r, fig.width=2.440945, fig.height= 4.429134}
# Calculate communication score by adding communication for each metabolite. by neck region
relative_comm_flow <- mebocost_comm_filtered %>%
  group_by(Metabolite_Name, neck_region) %>%
  summarise(sum_commu_score = sum(Commu_Score)) %>%
  group_by(Metabolite_Name) %>%
  mutate(relative_score = round(sum_commu_score / sum(sum_commu_score), 3))

# Order from most Deep enriched to most Superficial enriched
order <- relative_comm_flow %>%
  dplyr::filter(str_detect(neck_region, "Deep")) %>%
  arrange(relative_score) %>%
  pull(Metabolite_Name)
relative_comm_flow <- relative_comm_flow %>% mutate(Metabolite_Name = factor(Metabolite_Name, levels = order))

# Plot stacked proportions
gg_mebo_comm_flow_stacked <- ggplot(relative_comm_flow, aes(x = Metabolite_Name, y = relative_score)) +
  geom_col(aes(fill = neck_region)) +
  geom_hline(yintercept = .5, linetype = "dashed", color = "grey50") +
  theme_nature_metabolism() +
  scale_fill_manual(values = neck.color, name = NULL) +
  labs(title = "Metabolite Signaling", subtitle = "WAd as Receiver", y = "Relative information flow", x = NULL) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()
# ```

# ## --- Plotting and save ---
# ```{r, fig.width=2.165354, fig.height= 2.318898}
mm_to_in(c(55, 58.9))
gg_mebo_comm_flow_stacked <- gg_mebo_comm_flow_stacked

gg_mebo_comm_flow_stacked

ggsave(filename = file.path(OUTPUT_DIR, "5K_mebocost_wad_relative_information_flow_stacked.pdf"), plot = gg_mebo_comm_flow_stacked, width = 55, height = 58.9, units = "mm", dpi = 300)
# ```

# # --- Panel S4I. Mebocost relative communication flow (Raw) ---
# ## --- Create plot ---
# ```{r, fig.width=2.440945, fig.height= 4.429134}
# Order from strongest communication
order <- relative_comm_flow %>%
  group_by(Metabolite_Name) |>
  summarise(comm_strengt_total = sum(sum_commu_score)) |>
  arrange(comm_strengt_total) |>
  pull(Metabolite_Name)

relative_comm_flow <- relative_comm_flow %>% mutate(Metabolite_Name = factor(Metabolite_Name, levels = order))


# Plot dodged raw communication scores
gg_mebo_comm_flow_dodged <- ggplot(
  relative_comm_flow |> mutate(Metabolite_Name = str_wrap(Metabolite_Name, width = 80)),
  aes(x = Metabolite_Name, y = sum_commu_score)
) +
  geom_col(aes(fill = neck_region), position = position_dodge()) +
  theme_nature_metabolism() +
  scale_fill_manual(values = neck.color, name = NULL) +
  labs(
    title = "MEBOCOST\nInformation Flow",
    subtitle = "WAds as Receivers",
    y = "Information flow", x = NULL
  ) +
  scale_x_discrete(limits = str_wrap(order, width = 80)) +
  theme(legend.position = "top") +
  coord_flip(clip = "off")
# ```

# ## --- Plotting and save ---
# ```{r, fig.width=1.574803, fig.height= 4.724409}
mm_to_in(c(40, 120))
gg_mebo_comm_flow_dodged <- gg_mebo_comm_flow_dodged +
  theme(
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    axis.text.y = element_text(angle = 315, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(30, 1, 1, 1)
  )
gg_mebo_comm_flow_dodged

ggsave(filename = file.path(OUTPUT_DIR, "S4I_mebocost_wad_information_flow_dodged.pdf"), plot = gg_mebo_comm_flow_dodged, width = 40, height = 120, units = "mm", dpi = 300)
# ```

# # --- Panel 5L. Mebocost CirclePlot ---

# ## --- Setup helper function to position labels ---
# ```{r}
## Helper function
make_balanced_label_df <- function(lay, label_radius = 1.22, y_pad = 0.92) {
  max_y <- label_radius * y_pad

  label_df <- lay %>%
    dplyr::mutate(
      side = dplyr::if_else(x >= 0, "right", "left")
    ) %>%
    dplyr::group_by(side) %>%
    dplyr::arrange(dplyr::desc(y), .by_group = TRUE) %>%
    dplyr::mutate(
      y_lab = seq(from = max_y, to = -max_y, length.out = dplyr::n()),
      x_lab = dplyr::if_else(
        side == "right",
        sqrt(pmax(label_radius^2 - y_lab^2, 0)),
        -sqrt(pmax(label_radius^2 - y_lab^2, 0))
      ),
      hjust = dplyr::if_else(side == "right", 0, 1)
    ) %>%
    dplyr::ungroup()

  label_df
}
# ```

# ## --- Create plots ---
# ```{r}
meta <- s.object@meta.data

celltype_map <- meta %>%
  dplyr::select(cell_type, cell_type_short) %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  tibble::deframe()

rename_celltypes <- function(x, map) {
  x_clean <- stringr::str_remove(x, ".*~ ")
  dplyr::coalesce(unname(map[x_clean]), x_clean)
}

plot_list <- list()

for (metabolite_name in c("Dehydroepiandrosterone")) {
  ## Global max edge weight across depots for this metabolite
  max_edge <- mebocost_comm_filtered %>%
    dplyr::filter(
      stringr::str_detect(Metabolite_Name, metabolite_name),
      stringr::str_detect(neck_region, "Deep|Superficial"),
      stringr::str_detect(Receiver, "White Adipo")
    ) %>%
    dplyr::group_by(Sender, Receiver, Metabolite_Name, neck_region) %>%
    dplyr::summarise(edge_weight = sum(Commu_Score), .groups = "drop") %>%
    dplyr::pull(edge_weight) %>%
    max(na.rm = TRUE)

  for (depot in c("Deep", "Superficial")) {
    ## Edges for this panel
    edges <- mebocost_comm_filtered %>%
      dplyr::filter(
        stringr::str_detect(Metabolite_Name, metabolite_name),
        stringr::str_detect(neck_region, depot),
        stringr::str_detect(Receiver, "White Adipo")
      ) %>%
      dplyr::group_by(Sender, Receiver, Metabolite_Name) %>%
      dplyr::summarise(edge_weight = sum(Commu_Score), .groups = "drop") %>%
      dplyr::mutate(
        Sender = rename_celltypes(Sender, celltype_map),
        Receiver = rename_celltypes(Receiver, celltype_map)
      ) %>%
      dplyr::select(Sender, Receiver, edge_weight)

    edges <- mebocost_comm_filtered %>%
      dplyr::filter(
        str_detect(Metabolite_Name, metabolite_name),
        str_detect(neck_region, depot),
        str_detect(Receiver, "White Adipo")
      ) %>%
      group_by(Sender, Receiver, Metabolite_Name) %>%
      summarise(edge_weight = sum(Commu_Score), .groups = "drop") %>%
      mutate(
        Sender = rename_celltypes(Sender, celltype_map),
        Receiver = rename_celltypes(Receiver, celltype_map)
      ) %>%
      select(Sender, Receiver, edge_weight) %>%
      arrange(desc(edge_weight))

    ## All nodes, consistently renamed
    all_nodes <- unique(c(mebocost_comm$Sender, mebocost_comm$Receiver))
    all_nodes <- rename_celltypes(all_nodes, celltype_map)

    nodes_df <- tibble::tibble(name = unique(all_nodes))

    ## Graph
    g <- igraph::graph_from_data_frame(edges, vertices = nodes_df, directed = TRUE)
    g_tbl <- tidygraph::as_tbl_graph(g)

    node_from <- g_tbl %>%
      tibble::as_tibble() %>%
      dplyr::mutate(from = dplyr::row_number())

    new_g_tbl <- g_tbl %>%
      tidygraph::activate(edges) %>%
      dplyr::left_join(node_from, by = "from")

    ## Palette aligned to actual node order
    mebo.palette <- palette.use[nodes_df$name]

    ## Build layout once and extract coordinates
    lay <- ggraph::create_layout(new_g_tbl, layout = "circle")

    ## Fixed radial label placement
    label_df <- make_balanced_label_df(lay, label_radius = 1.22, y_pad = 0.999)

    ## Plot
    mebo_graph <- ggraph::ggraph(lay) +
      ggraph::geom_edge_arc(
        aes(color = name, width = edge_weight),
        strength = 0.1,
        show.legend = FALSE
      ) +
      ggraph::geom_edge_loop(
        aes(color = name, width = edge_weight),
        strength = 0.8,
        show.legend = FALSE
      ) +
      ggraph::geom_node_point(
        aes(fill = name),
        shape = 21,
        color = "black",
        size = 2
      ) +
      ggplot2::geom_text(
        data = label_df,
        mapping = ggplot2::aes(
          x = x_lab,
          y = y_lab,
          label = name,
          hjust = hjust
        ),
        size = 6 / .pt,
        inherit.aes = FALSE
      ) +
      ggraph::scale_edge_width(
        range = c(0.1, 1.5),
        limits = c(0, max_edge),
        name = "Communication\nScore"
      ) +
      ggraph::scale_edge_colour_manual(values = mebo.palette, guide = "none") +
      ggplot2::scale_fill_manual(values = mebo.palette, guide = "none") +
      coord_equal(
        xlim = c(-1.55, 1.55),
        ylim = c(-1.55, 1.55),
        clip = "off"
      ) + theme_void() + theme(
        plot.margin = margin(1, 12, 1, 19)
      )

    plot_list[[paste(metabolite_name, depot)]] <- mebo_graph
    
    # -----------------------------
    #  Save plot 
    # -----------------------------
    

    ggplot2::ggsave(
      filename = file.path(
        OUTPUT_DIR,
        paste0("5L_mebocost_", tolower(metabolite_name), "_", depot, ".pdf")
      ),
      plot = mebo_graph,
      width = 40,
      height = 40,
      units = "mm",
      dpi = 600
    )
  }
}
# ```

# ## --- Plotting ---
# ```{r, fig.width=1.574803, fig.height= 1.574803}
mm_to_in(40)
plot_list
# ```
