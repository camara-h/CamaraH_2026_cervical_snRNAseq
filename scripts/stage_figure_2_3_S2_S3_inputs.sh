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
  cp -al "$src" "$DEST/$name"
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
  cp -al "$src" "$DEST/$name"
}

# Files
stage_file "output/3.heat_signature/data/heat_signature.Rdata"
stage_file "data/6.internal_datasets/Tseng_Datasets/Ruidan Xue's Microarray data/20 arrays for hBAT hWAT tissues/Supplementary_Table_cervical_at_microarray.xlsx"
stage_file "output/3.heat_signature/data/Lynes_et_al_correlation.csv"
stage_file "../Internal_datasets/Gavrila Lab Bulk RNAseq Data/logCPM_for_correlations.csv"
stage_file "output/3.heat_signature/data/Salej_et_al_logCPM_correlation.csv"
stage_file "output/3.heat_signature/data/genes_for_venn.Rdata"
stage_file "output/3.heat_signature/data/heat_signature.csv"
stage_file "output/3.heat_signature/1.signature_generation/enrichr_results/enriched_heat.rds"
stage_file "output/4.downstream_pipelines/vision/data/heat_threshold/wad_heat_threshold.rds"

stage_file "output/3.heat_signature/3.public_bulkrnaseq/data/GSE49795_Sondergaard_expression_matrix.rds"
stage_file "data/6.internal_datasets/Tseng_Datasets/Supplementary_Table_a41_fpkm.csv"
stage_file "data/6.internal_datasets/Tseng_Datasets/Supplementary_Table_a38_tpm.xlsx"
stage_file "data/5.published_datasets/Public_BulkRNAseq/Cero 2023/RNA-seq data set Cero et al.xlsx"

# Directories
stage_dir "output/3.heat_signature/3.public_bulkrnaseq/data/ssGSEA"
stage_dir "output/3.heat_signature/3.public_bulkrnaseq/data/Cleaned_Datasets"
stage_dir "output/3.heat_signature/3.public_bulkrnaseq/data/Experimental_Design"

echo
echo "Done. Staged files are in:"
echo "$DEST"
echo

find "$DEST" -maxdepth 2 -print | sort
du -sh "$DEST"