# MEBOCOST and COMPASS setup

This project uses a Python workflow in which MEBOCOST is run first, COMPASS is run on the exported average expression matrix, and MEBOCOST is then resumed using the COMPASS output directory.

## What this environment is for

This environment is intended to support:

- running the multisample MEBOCOST notebook
- exporting `avg_exp_mat_multi.tsv`
- running the COMPASS shell script
- re-importing COMPASS outputs into the MEBOCOST notebook
- exporting `communication_result_multi.csv` for downstream R analysis

## Recommended environment setup

Create the conda environment from the template:

```bash
conda env create -f environment.yml
conda activate mebocost-compass
python -m ipykernel install --user --name mebocost-compass --display-name "Python (mebocost-compass)"
```

If the install from the YAML fails because of upstream package changes, install the core packages first and then install MEBOCOST and COMPASS separately:

```bash
conda create -n mebocost-compass python=3.10 -c conda-forge -c bioconda
conda activate mebocost-compass
conda install -c conda-forge -c bioconda numpy pandas scipy matplotlib seaborn scanpy anndata jupyterlab notebook ipykernel tqdm scikit-learn openpyxl
pip install git+https://github.com/kaifuchenlab/MEBOCOST.git
pip install git+https://github.com/YosefLab/Compass.git
```

## Expected project paths

These scripts assume a project layout similar to:

```text
NeckAdipose/
├── data/
│   └── adata_integrated_soupXoutput_with_infer_forR_temp_Dec18_cleaned.h5ad
├── scripts/
│   └── mebocost/
│       ├── mebocost.conf
│       └── run_COMPASS_multi.sh
├── output/
│   └── 4.downstream_pipelines/
│       ├── mebocost/
│       │   └── data/
│       └── compass/
│           └── data/
└── ...
```

## Files you should check before running

Update these paths in your notebook and shell script if needed:

- `PROJECT_ROOT`
- `INPUT_H5AD`
- `CONFIG_PATH`
- `MEBO_OUTPUT_DIR`
- `COMPASS_OUTPUT_DIR`
- `COMPASS_BIN`

## Workflow order

### 1. Run the MEBOCOST multisample notebook
This should generate at least:

- `NeckAdipose_commu_multi.pk`
- `avg_exp_mat_multi.tsv`

### 2. Run COMPASS
Run the shell script after `avg_exp_mat_multi.tsv` exists:

```bash
bash scripts/mebocost/run_COMPASS_multi.sh
```

This should generate a stable output directory such as:

- `output/4.downstream_pipelines/compass/data/compass_res_multi_Deep_Sup`

### 3. Resume the MEBOCOST notebook
Point `COMPASS_FOLDER` to the COMPASS output directory and rerun the downstream cells that:

- constrain communication with COMPASS flux results
- run differential communication analysis
- export `results_commu_multi.csv`
- export `communication_result_multi.csv`

## Notes on reproducibility

This project uses stable, non-dated file paths for the final pipeline so downstream scripts can rely on predictable locations.

If you want to preserve alternative runs, copy the output folders manually to an archive location rather than date-stamping the canonical pipeline outputs.

## Notes on `setup.py`

Do not include the upstream MEBOCOST `setup.py` in this project unless you are distributing a packaged fork of MEBOCOST itself. This analysis repo only needs:

- the environment definition
- the project-specific `mebocost.conf`
- the notebooks and shell scripts used in your workflow
