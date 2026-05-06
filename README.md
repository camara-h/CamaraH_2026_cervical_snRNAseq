# Camara et al. 2026 - Analysis of Human Cervical Adipose Tissue single-nuclei RNA-seq + ATAC-seq
The repository is intended for transparency and reproducibility, not as a general-purpose software package.

Raw human snRNA-seq data are not included because of data-use/privacy restrictions.
Processed or public data access information is provided in the manuscript Data Availability statement.

## Scope

This repository is provided to document the analyses performed for the associated manuscript. 
It is not intended to be a general-purpose software package, and the code may require local path updates, external data access, or manual steps described in the README.

## What this repository contains

Analysis code for a human cervical adipose tissue single nucleus manuscript focused on depot specific cell states, brown adipocyte biology, thermogenic signature construction, and downstream cell cell communication and multimodal analyses.

This repository organizes the main analysis into a chronological pipeline:

1. environment setup
2. preprocessing and object conversion
3. Seurat object characterization
4. HEAT signature generation and validation
5. downstream pathway and communication analyses
6. multiome preprocessing and Signac analysis
7. manuscript figure preparation

The code is structured for reproducibility, but some steps still depend on data or external platforms that are not deposited directly in this repository, such as protected source datasets, CIBERSORTx runs, and large intermediate objects.

## Repository structure

```text
scripts/
  0.environment_setup/
  1.processing/
  2.sobject_characterization/
  3.heat_signature/
  4.downstream_pipelines/
  5.multiome_analysis/
  6.figure_preparation/
renv.lock
00.Tseng_neck_adipose.Rproj
```

### `scripts/0.environment_setup`
Project level setup scripts, including:
- R package installation
- plotting themes and shared color palettes used across manuscript figures

### `scripts/1.processing`
Early data preparation and conversion steps, including:
- Seurat object preprocessing and annotation
- conversion between Seurat and AnnData formats

### `scripts/2.sobject_characterization`
Characterization of the cleaned neck adipose Seurat object, including:
- composition analysis
- subtype reference mapping

### `scripts/3.heat_signature`
Generation and validation of the HEAT signature across multiple data sources:
- internal bulk and microarray correlations
- public snRNA-seq datasets
- public bulk RNA-seq datasets
- GTEx analyses

### `scripts/4.downstream_pipelines`
Downstream biological analyses, including:
- VISION signature scoring
- CIBERSORTx bootstrap reference generation and result extraction
- pseudobulk differential expression and GSEA
- CellChat
- MEBOCOST and COMPASS
- RNA velocity

### `scripts/5.multiome_analysis`
Multiome preprocessing and Signac workflows:
- Cell Ranger ARC style preprocessing
- ATAC and GEX processing
- WNN integration
- path update helpers for fragment files

### `scripts/6.figure_preparation`
Figure assembly scripts for the main manuscript and supplementary figures.

## Main biological workflow

At a high level, the analysis proceeds as follows.

### 1. Create and clean the main Seurat object
The main neck adipose object is processed and annotated in:

- `scripts/1.processing/1.1.seurat_object_preprocessing_cleaned.Rmd`

This script produces the cleaned Seurat object used by most downstream analyses.

### 2. Generate the HEAT signature
The HEAT signature is constructed by intersecting:

- genes positively correlated with UCP1 in internal deep neck bulk RNA-seq
- genes positively correlated with UCP1 in internal microarray data
- genes enriched in brown adipocytes in the snRNA-seq dataset

Main script:

- `scripts/3.heat_signature/3.1.signature_generation/3.1.1.heat_signature_generation_cleaned.Rmd`

### 3. Validate the signature in public datasets
The repository contains public dataset specific scripts for evaluating the HEAT signature in:
- published snRNA-seq datasets
- published bulk RNA-seq datasets
- GTEx

These scripts live under:

- `scripts/3.heat_signature/3.2.public_snrnaseq/`
- `scripts/3.heat_signature/3.3.public_bulkrnaseq/`
- `scripts/3.heat_signature/3.4.gtex/`

### 4. Run downstream pathway and communication analyses
Several complementary downstream analyses are included:

- VISION
- CIBERSORTx bootstrap deconvolution workflow
- pseudobulk differential expression
- CellChat
- MEBOCOST plus COMPASS
- RNA velocity

These live under:

- `scripts/4.downstream_pipelines/`

### 5. Generate manuscript figures
The figure scripts under `scripts/6.figure_preparation/` assemble the main and supplementary figures using outputs from earlier pipeline stages.

## Canonical inputs and outputs

The cleaned scripts are written assuming a project structure with top level `data/` and `output/` directories.

### Expected input directory
Large source objects and downloaded datasets are expected under:

```text
data/
```

Examples include:
- cleaned Seurat objects
- published dataset inputs
- metadata tables
- GMT signatures
- AnnData objects
- multiome fragment bundles

### Expected output directory
Intermediate and final analysis outputs are written under:

```text
output/
```

Subfolders generally follow the same stage numbering as the scripts.

## Reproducibility notes

### R environment
This project uses `renv` for package management.

Recommended setup:

```r
renv::restore()
```

Open the project from:

- `00.Tseng_neck_adipose.Rproj`

to ensure relative paths resolve correctly.

### Python environments
Some downstream workflows use Python specific environments rather than `renv`, especially:

- MEBOCOST and COMPASS
- RNA velocity with scVelo and velocyto

See:
- `scripts/4.downstream_pipelines/4.5.mebocost/environment.yml`
- `scripts/4.downstream_pipelines/4.5.mebocost/README_MEBOCOST_COMPASS.md`
- `scripts/4.downstream_pipelines/4.6.rna_velocity/README_run_velocyto_batch.md`

## Important manual or external steps

Some workflows are not fully end to end inside the repository and require manual or external execution.

### CIBERSORTx
The bootstrap reference script prepares sampled reference matrices, but the actual deconvolution step is performed externally in CIBERSORTx. Results are then downloaded and re-imported for aggregation.

Relevant scripts:
- `scripts/4.downstream_pipelines/4.2.cibersort/4.2.1.bootstrap_reference.Rmd`
- `scripts/4.downstream_pipelines/4.2.cibersort/4.2.2.extract_cibersortx_proportions_cleaned.Rmd`

### RNA velocity
The RNA velocity workflow depends on:
- per sample velocyto loom files
- a curated metadata table
- exported scVI embeddings
- a curated cell list

Relevant scripts:
- `scripts/4.downstream_pipelines/4.6.rna_velocity/4.6.1.run_velocyto_batch.py`
- `scripts/4.downstream_pipelines/4.6.rna_velocity/4.6.2.velocity_three_panels_cleaned.py`

### MEBOCOST and COMPASS
The MEBOCOST workflow exports an average expression matrix for COMPASS, runs COMPASS as a separate step, then resumes MEBOCOST using the COMPASS output folder.

Relevant files:
- `scripts/4.downstream_pipelines/4.5.mebocost/4.5.1.neck_multisample_mccc_analysis_cleaned.ipynb`
- `scripts/4.downstream_pipelines/4.5.mebocost/run_COMPASS_multi_cleaned.sh`

## Suggested execution order

A practical order for readers trying to understand or rerun the analysis is:

1. `scripts/0.environment_setup/`
2. `scripts/1.processing/`
3. `scripts/2.sobject_characterization/`
4. `scripts/3.heat_signature/`
5. `scripts/4.downstream_pipelines/`
6. `scripts/5.multiome_analysis/`
7. `scripts/6.figure_preparation/`

Within each folder, the numeric prefixes usually reflect the intended order.

## Public release notes

This repository contains analysis code and expects local or institutional access to large input data files. Not all datasets can be redistributed directly here.

Before reusing the code, users should expect to:
- recreate the local `data/` directory structure
- adjust any remaining machine specific paths if needed
- install the required R and Python environments
- supply external inputs for workflows that depend on protected or manually downloaded data

## Citation

If you use code or ideas from this repository, please cite the corresponding manuscript when available.

## License

The analysis code in this repository is released under the MIT License. See
`LICENSE` for details.

This license applies to the code, notebooks, scripts, and workflow files in this repository. It does not grant rights to redistribute protected human subject data, controlled-access datasets, third-party datasets, software, database
files, or external resources used by the analyses. 
Data access and reuse are governed by the original data providers, institutional approvals, and the manuscript Data Availability statement.


