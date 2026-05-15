#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Stage Code Ocean input data by figure
# ============================================================
#
# Purpose:
#   Create a local upload-ready data folder for Code Ocean, organized
#   by figure, using hard links to avoid duplicating local disk usage.
#
# Output structure:
#   ocean_code_data/
#     figure_2_3_S2_S3/
#     Figure_4/
#     Figure_S1/
#     Figure_5_S4/
#     Figure_6_S5/
#
# Intended Code Ocean structure after upload:
#   /data/figure_2_3_S2_S3/
#   /data/Figure_4/
#   /data/Figure_S1/
#   /data/Figure_5_S4/
#   /data/Figure_6_S5/
#
# Notes:
#   - This script uses cp -al, which creates hard links.
#   - Hard links save local space and behave like regular files during upload.
#   - Do not edit files inside ocean_code_data, because hard-linked staged files
#     point to the same file content as the originals.
#   - Existing staged files/folders are skipped by default.
#   - If a source file changed and you want to restage it, delete the staged file
#     or delete the entire staged figure folder and rerun this script.
#
# Run from the project root:
#   bash scripts/stage_code_ocean_inputs_by_figure.sh
#
# ============================================================

ROOT_DEST="ocean_code_data"

FIG_2_3_S2_S3="$ROOT_DEST/figure_2_3_S2_S3"
FIG_4="$ROOT_DEST/Figure_4"
FIG_S1="$ROOT_DEST/Figure_S1"
FIG_5_S4="$ROOT_DEST/Figure_5_S4"
FIG_6_S5="$ROOT_DEST/Figure_6_S5"

SKIP_EXISTING=true

mkdir -p \
  "$FIG_2_3_S2_S3" \
  "$FIG_4" \
  "$FIG_S1" \
  "$FIG_5_S4" \
  "$FIG_6_S5"

stage_file_to() {
  dest="$1"
  src="$2"
  name="${3:-$(basename "$src")}"

  if [[ ! -f "$src" ]]; then
    echo "MISSING FILE: $src"
    exit 1
  fi

  mkdir -p "$dest"
  target="$dest/$name"

  if [[ -e "$target" && "$SKIP_EXISTING" == "true" ]]; then
    echo "SKIPPING existing file:"
    echo "  $target"
    return 0
  fi

  if [[ -e "$target" && "$SKIP_EXISTING" != "true" ]]; then
    rm -f "$target"
  fi

  echo "Staging file:"
  echo "  from: $src"
  echo "  to:   $target"

  cp -al "$src" "$target"
}

stage_dir_to() {
  dest="$1"
  src="$2"
  name="${3:-$(basename "$src")}"

  if [[ ! -d "$src" ]]; then
    echo "MISSING DIR: $src"
    exit 1
  fi

  mkdir -p "$dest"
  target="$dest/$name"

  if [[ -e "$target" && "$SKIP_EXISTING" == "true" ]]; then
    echo "SKIPPING existing directory:"
    echo "  $target"
    return 0
  fi

  if [[ -e "$target" && "$SKIP_EXISTING" != "true" ]]; then
    rm -rf "$target"
  fi

  echo "Staging directory:"
  echo "  from: $src"
  echo "  to:   $target"

  cp -al "$src" "$target"
}

make_manifest_for() {
  dest="$1"

  echo
  echo "Creating manifest for: $dest"
  find "$dest" -type f | sort > "$dest/input_manifest.txt"

  echo "Size:"
  du -sh "$dest"

  echo "Number of files:"
  find "$dest" -type f | wc -l
}

echo "Staging Code Ocean data inputs..."
echo "Root destination: $ROOT_DEST"
echo "Skip existing files/folders: $SKIP_EXISTING"
echo

# ============================================================
# Seurat GEX
# ============================================================

stage_file_to "$ROOT_DEST" \
  "data/adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_cleaned.rds"

#Rename  
mv "$ROOT_DEST/adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_cleaned.rds" "$ROOT_DEST/cervical_at_gex_seurat.rds"

# ============================================================
# Figure 2, Figure 3, Figure S2, Figure S3
# ============================================================

echo
echo "===== Staging Figure 2, Figure 3, Figure S2, Figure S3 ====="

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/data/heat_signature.Rdata"

stage_file_to "$FIG_2_3_S2_S3" \
  "data/6.internal_datasets/Tseng_Datasets/Ruidan Xue's Microarray data/20 arrays for hBAT hWAT tissues/Supplementary_Table_cervical_at_microarray.xlsx"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/data/Lynes_et_al_correlation.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "../Internal_datasets/Gavrila Lab Bulk RNAseq Data/logCPM_for_correlations.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/data/Salej_et_al_logCPM_correlation.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/data/genes_for_venn.Rdata"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/data/heat_signature.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/1.signature_generation/enrichr_results/enriched_heat.rds"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/4.downstream_pipelines/vision/data/heat_threshold/wad_heat_threshold.rds"

stage_dir_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/ssGSEA"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/GSE49795_Sondergaard_expression_matrix.rds"

stage_dir_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Cleaned_Datasets"

stage_dir_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Experimental_Design"

stage_file_to "$FIG_2_3_S2_S3" \
  "data/6.internal_datasets/Tseng_Datasets/Supplementary_Table_a41_fpkm.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "data/6.internal_datasets/Tseng_Datasets/Supplementary_Table_a38_tpm.xlsx"

stage_file_to "$FIG_2_3_S2_S3" \
  "data/5.published_datasets/Public_BulkRNAseq/Cero 2023/RNA-seq data set Cero et al.xlsx"

stage_dir_to "$FIG_2_3_S2_S3" \
  "output/4.downstream_pipelines/vision/data/vision_signature_scores"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/1.processing/data/cleaned_objects/AngueiraA_Dataset_adipo_recluster.rds"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/3.public_snrnaseq/data/Angueira_HEATvUCP1.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/4.downstream_pipelines/vision/data/vision_signature_scores/ShamsiF_Dataset_VISION_sigScore.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/4.downstream_pipelines/vision/data/cleaned_datasets_mouse/ShamsiF_Dataset.rds"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/4.downstream_pipelines/cibersort/data/deconvolution_cellular_proportions_table.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/5.gtex/data/HEAT_metadata_filtered.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/5.gtex/data/stratified_aov_statistics.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "output/3.heat_signature/5.gtex/data/significant_stratified_aov.csv"

stage_file_to "$FIG_2_3_S2_S3" \
  "local_data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDD.xlsx"

# ============================================================
# Figure 4
# ============================================================

echo
echo "===== Staging Figure 4 ====="

stage_file_to "$FIG_4" \
  "output/4.downstream_pipelines/de_analysis/data/deseq2_stats.csv"

stage_file_to "$FIG_4" \
  "output/4.downstream_pipelines/de_analysis/data/HallMark_GSEA.csv"

stage_file_to "$FIG_4" \
  "output/4.downstream_pipelines/de_analysis/data/GO_GSEA.csv"

stage_file_to "$FIG_4" \
  "output/4.downstream_pipelines/de_analysis/data/Reactome_GSEA.csv"

stage_file_to "$FIG_4" \
  "output/4.downstream_pipelines/de_analysis/data/heatmap_input.rds"

# ============================================================
# Figure S1
# ============================================================

echo
echo "===== Staging Figure S1 ====="

stage_file_to "$FIG_S1" \
  "../Internal_Datasets/Gavrila Lab Bulk RNAseq Data/TPM_values.csv"

stage_file_to "$FIG_S1" \
  "data/5.published_datasets/Public_BulkRNAseq/Din 2018/Din2018_countmatrix.txt"

stage_file_to "$FIG_S1" \
  "data/6.internal_datasets/Tseng_Datasets/Ruidan Xue's Microarray data/20 arrays for hBAT hWAT tissues/Supplementary_Table_cervical_at_microarray.xlsx"

stage_file_to "$FIG_S1" \
  "../Output/1 - Figure and Analysis/36.ReferenceMapping/aspc_combined.rds"

stage_file_to "$FIG_S1" \
  "../Output/1 - Figure and Analysis/36.ReferenceMapping/adipo_combined.rds"

stage_file_to "$FIG_S1" \
  "output/2.sobject_characterization/propeller/data/propeller_cellular_proportions_table.csv"

stage_file_to "$FIG_S1" \
  "output/2.sobject_characterization/propeller/data/proportions_stat_table.csv"

# ============================================================
# Figure 5 and Figure S4
# ============================================================

echo
echo "===== Staging Figure 5 and Figure S4 ====="

stage_file_to "$FIG_5_S4" \
  "output/4.downstream_pipelines/cellchat/data/cellchat_deep_pop.size.false_cell_type_short.rds"

stage_file_to "$FIG_5_S4" \
  "output/4.downstream_pipelines/cellchat/data/cellchat_superficial_pop.size.false_cell_type_short.rds"

stage_file_to "$FIG_5_S4" \
  "output/4.downstream_pipelines/cellchat/data/cellchat_merged_pop.size.false_cell_type_short.rds"

stage_file_to "$FIG_5_S4" \
  "output/4.downstream_pipelines/mebocost/data/communication_result_multi.csv"

# ============================================================
# Figure 6 and Figure S5
# ============================================================

echo
echo "===== Staging Figure 6 and Figure S5 ====="

stage_file_to "$FIG_6_S5" \
  "output/4.downstream_pipelines/de_analysis/data/deseq2_stats.csv"

stage_file_to "$FIG_6_S5" \
  "output/4.downstream_pipelines/vision/data/vision_objects/CamaraH_Dataset_VISION.rds"

stage_file_to "$FIG_6_S5" \
  "output/4.downstream_pipelines/de_analysis/data/CamaraH_pseudobulk_log_CPM_count_matrix_for_heatmap.rds"

stage_file_to "$FIG_6_S5" \
  "output/4.downstream_pipelines/de_analysis/data/CamaraH_metadata_for_heatmap.rds"

stage_file_to "$FIG_6_S5" \
  "output/4.downstream_pipelines/de_analysis/data/CamaraH_pseudobulk_count_matrix_for_ggplot.rds"

stage_dir_to "$FIG_6_S5" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Cleaned_Datasets"

stage_dir_to "$FIG_6_S5" \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Experimental_Design"

stage_file_to "$FIG_6_S5" \
  "data_onedrive/9.multiome_data/multiome_seurat.rds"

# ============================================================
# Manifests and final summary
# ============================================================

echo
echo "===== Creating manifests ====="

make_manifest_for "$FIG_2_3_S2_S3"
make_manifest_for "$FIG_4"
make_manifest_for "$FIG_S1"
make_manifest_for "$FIG_5_S4"
make_manifest_for "$FIG_6_S5"

echo
echo "All staging complete."
echo "Upload the contents of this folder to Code Ocean data:"
echo "$ROOT_DEST"
echo

echo "Created data folders:"
find "$ROOT_DEST" -maxdepth 2 -type d | sort

echo
echo "Combined staged size:"
du -sh "$ROOT_DEST"

echo
echo "Combined number of staged files:"
find "$ROOT_DEST" -type f | wc -l
