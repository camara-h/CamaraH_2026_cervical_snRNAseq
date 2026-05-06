library(Seurat)
library(SeuratDisk)
library(optparse)

# --- Command-line argument parsing ---
option_list <- list(
  make_option("--input", type="character", help="Path to the Seurat object."),
  make_option("--output", type="character", help="Destination output folder ")
)

ensure_trailing_slash <- function(path) {
  if (grepl("/$", path)) path else paste0(path, "/")
}
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# --- Parse input paths ---
path = opt$input
output = opt$output
output <- ensure_trailing_slash(output)

dirpath = dirname(path)
filename_ext <- basename(path)  # "file.rds"
filename_raw <- tools::file_path_sans_ext(filename_ext) # "file"
print(filename_raw)

h5 = paste(output, filename_raw, ".h5Seurat", sep = "")

object = readRDS(path)
UpdateSeuratObject(object)
SaveH5Seurat(object, filename = h5)
Convert(h5, dest = "h5ad")

file.remove(h5)
