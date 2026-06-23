library(Seurat)
library(tidyverse)
library(here)

s.object <- readRDS(here("data/adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_cleaned.rds"))

s.object@meta.data <- s.object@meta.data |> 
  select(nCount_log_10k_norm, nFeature_log_10k_norm, X_scvi_batch, X_scvi_labels, leiden_scVI, X_scvi_raw_norm_scaling, 
         sample, subject, neck_region, cell_type, cell_type_short, RNA_snn_res.0.8, RNA_snn_res.1.6, 
         Age, Gender, Race, Surgery, BMI, Fasting, Glucose, HbA1c, Leptin, Cpeptide, Insulin, Adiponectin, `HOMA-IR`, obesity)
  
write_rds(s.object, here("ocean_code_data/cervical_at_gex_seurat.rds"))
