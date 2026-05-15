#!/usr/bin/env bash
set -euo pipefail

DEST="code_space_data/figure_2_3_S2_S3"
mkdir -p "$DEST"

stage_file() {
  src="$1"
  name="${2:-$(basename "$src")}"

  if [[ ! -f "$src" ]]; then
    echo "MISSING FILE: $src"
    exit 1
  fi

  echo "Staging file: $name"
  cp -a "$src" "$DEST/$name"
}

stage_dir() {
  src="$1"
  name="${2:-$(basename "$src")}"

  if [[ ! -d "$src" ]]; then
    echo "MISSING DIR: $src"
    exit 1
  fi

  echo "Staging directory: $name"
  rm -rf "$DEST/$name"
  cp -a "$src" "$DEST/$name"
}

echo "Staging Figure 2, Figure 3, Figure S2, Figure S3 inputs..."
echo "Destination: $DEST"
echo

# ----------------------------
# HEAT signature generation
# ----------------------------

stage_file \
  "output/3.heat_signature/data/heat_signature.Rdata"

stage_file \
  "data/6.internal_datasets/Tseng_Datasets/Ruidan Xue's Microarray data/20 arrays for hBAT hWAT tissues/Supplementary_Table_cervical_at_microarray.xlsx"

stage_file \
  "output/3.heat_signature/data/Lynes_et_al_correlation.csv"

stage_file \
  "../Internal_datasets/Gavrila Lab Bulk RNAseq Data/logCPM_for_correlations.csv"

stage_file \
  "output/3.heat_signature/data/Salej_et_al_logCPM_correlation.csv"

stage_file \
  "output/3.heat_signature/data/genes_for_venn.Rdata"

stage_file \
  "output/3.heat_signature/data/heat_signature.csv"


# ----------------------------
# HEAT Enrichr
# ----------------------------

stage_file \
  "output/3.heat_signature/1.signature_generation/enrichr_results/enriched_heat.rds"


# ----------------------------
# HEAT threshold
# ----------------------------

stage_file \
  "output/4.downstream_pipelines/vision/data/heat_threshold/wad_heat_threshold.rds"


# ----------------------------
# Bulk RNA-seq
# ----------------------------

stage_dir \
  "output/3.heat_signature/3.public_bulkrnaseq/data/ssGSEA"

stage_file \
  "output/3.heat_signature/3.public_bulkrnaseq/data/GSE49795_Sondergaard_expression_matrix.rds"

stage_dir \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Cleaned_Datasets"

stage_dir \
  "output/3.heat_signature/3.public_bulkrnaseq/data/Experimental_Design"

stage_file \
  "data/5.published_datasets/Public_BulkRNAseq/Cero 2023/RNA-seq data set Cero et al.xlsx"


# ----------------------------
# VISION
# ----------------------------

stage_dir \
  "output/4.downstream_pipelines/vision/data/vision_signature_scores"


# ----------------------------
# Angueira dataset
# ----------------------------

stage_file \
  "output/1.processing/data/cleaned_objects/AngueiraA_Dataset_adipo_recluster.rds"

stage_file \
  "output/3.heat_signature/3.public_snrnaseq/data/Angueira_HEATvUCP1.csv"


# ----------------------------
# Mouse datasets
# ----------------------------

stage_file \
  "output/4.downstream_pipelines/vision/data/vision_signature_scores/ShamsiF_Dataset_VISION_sigScore.csv"

stage_file \
  "output/4.downstream_pipelines/vision/data/cleaned_datasets_mouse/ShamsiF_Dataset.rds"


# ----------------------------
# CIBERSORTx
# ----------------------------

stage_file \
  "output/4.downstream_pipelines/cibersort/data/deconvolution_cellular_proportions_table.csv"


# ----------------------------
# GTEx
# ----------------------------

stage_file \
  "output/3.heat_signature/5.gtex/data/HEAT_metadata_filtered.csv"

stage_file \
  "output/3.heat_signature/5.gtex/data/stratified_aov_statistics.csv"

stage_file \
  "output/3.heat_signature/5.gtex/data/significant_stratified_aov.csv"

stage_file \
  "local_data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDD.xlsx"


# ----------------------------
# Manifest and summary
# ----------------------------

echo
echo "Creating manifest..."
find "$DEST" -type f | sort > "$DEST/input_manifest.txt"

echo
echo "Done staging inputs."
echo "Folder:"
echo "$DEST"
echo

echo "Size:"
du -sh "$DEST"
echo

echo "Files staged:"
find "$DEST" -maxdepth 2 -type f | sort