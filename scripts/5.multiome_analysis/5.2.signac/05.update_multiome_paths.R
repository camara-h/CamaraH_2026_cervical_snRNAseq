library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(ggplot2)
library(tools)

# Function to update fragment paths in a Signac/Seurat object
update_fragment_paths <- function(seurat_obj, assay = "peaks", new_base_dir) {
  stopifnot(assay %in% names(seurat_obj@assays))
  
  frag_list <- seurat_obj[[assay]]@fragments
  
  if (length(frag_list) == 0) {
    message("No fragment objects found in assay: ", assay)
    return(seurat_obj)
  }
  
  old_paths <- vapply(frag_list, function(x) x@path, character(1))
  
  # Keep sample folder + outs + filename
  rel_paths <- vapply(old_paths, function(p) {
    parts <- strsplit(p, "/", fixed = TRUE)[[1]]
    paste(tail(parts, 3), collapse = "/")
  }, character(1))
  
  new_paths <- file.path(new_base_dir, rel_paths)
  
  # Update paths
  for (i in seq_along(frag_list)) {
    frag_list[[i]]@path <- new_paths[i]
  }
  
  seurat_obj[[assay]]@fragments <- frag_list
  
  # Print mapping for review
  path_map <- data.frame(
    old_path = old_paths,
    relative_path = rel_paths,
    new_path = new_paths,
    exists = file.exists(new_paths),
    stringsAsFactors = FALSE
  )
  
  print(path_map)
  
  missing_paths <- new_paths[!file.exists(new_paths)]
  if (length(missing_paths) > 0) {
    warning(
      "Some updated fragment paths do not exist yet:\n",
      paste(missing_paths, collapse = "\n")
    )
  }
  
  return(seurat_obj)
}

# -----------------------------------------------------------------------------
# Input paths
# -----------------------------------------------------------------------------

seurat_path <- here("data/multiome_seurat.rds")

stopifnot(file.exists(seurat_path))

# New local base directory where you unpacked the fragment bundle
new_base_dir <- here("data_onedrive/9.multiome_dataa/atac_fragments_bundle")



# -----------------------------------------------------------------------------
# Load objects
# -----------------------------------------------------------------------------

query_obj <- readRDS(seurat_path)

query_obj <- readRDS()

query_obj <- update_fragment_paths(
  seurat_obj = query_obj,
  assay = "peaks",
  new_base_dir = new_base_dir
)

vapply(query_obj[["peaks"]]@fragments, function(x) x@path, character(1))

saveRDS(query_obj, paste0(file_path_sans_ext(seurat_path), "_updated_fragment_paths.rds"))