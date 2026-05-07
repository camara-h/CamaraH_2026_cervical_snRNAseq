

# --- Setup ---
library(Seurat)
library(CellChat)
library(tidyverse)

CELLCHAT_DEEP <- here("output/4.downstream_pipelines/cellchat/data/cellchat_deep_pop.size.false_cell_type_short.rds")
CELLCHAT_SUPERFICIAL <- here("output/4.downstream_pipelines/cellchat/data/cellchat_superficial_pop.size.false_cell_type_short.rds")
CELLCHAT_MERGED <- here("output/4.downstream_pipelines/cellchat/data/cellchat_merged_pop.size.false_cell_type_short.rds")


OUTPUT_DATA_DIR <- here("output", "4.downstream_pipelines", "cellchat", "data")

## --- Load the objects and create CC object
deep_cellchat <- readRDS(CELLCHAT_DEEP)
superficial_cellchat <- readRDS(CELLCHAT_SUPERFICIAL)
cellchat_merged <- readRDS(CELLCHAT_MERGED)


## --- Create object_list
object_list <- list(Deep = deep_cellchat, Superficial = superficial_cellchat)



## --- Extract signaling information

# Extract signaling information
df.sup <- as.data.frame(subsetCommunication(superficial_cellchat, thresh = 1))
df.deep <- as.data.frame(subsetCommunication(deep_cellchat, thresh = 1))

# Merge and filter significant
df.merge <- rbind(
  df.sup |> mutate(depot = "Superficial"),
  df.deep |> mutate(depot = "Deep")
)

write_csv(df.merge, file.path(OUTPUT_DATA_DIR, "cell_chat_communication_df.csv"))