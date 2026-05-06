# RNA velocity loom generation with velocyto

This script runs `velocyto` across a batch of Cell Ranger sample folders and generates one `.loom` file per sample.

It replaces the older split scripts:

- `velocyto1.py`
- `velocyto2.py`
- `velocyto3.py`
- `velocyto4.py`
- `velocyto5.py`
- `velocyto6.py`

## What this script does

For each selected sample folder, the script:

1. finds the Cell Ranger BAM file
2. finds the barcode file
3. runs `velocyto run`
4. writes one output `.loom` file per sample

The script supports two barcode modes:

- `filtered`: uses `outs/filtered_feature_bc_matrix/barcodes.tsv`
- `raw`: uses `outs/raw_feature_bc_matrix/barcodes.tsv`

## Which mode is the main workflow?

The canonical workflow used in the downstream RNA velocity analysis is:

- `--barcode-mode filtered`

This is because the downstream velocity script reads loom files from:

- `.../results/alignment_new/velocyto_results`

The `raw` mode is kept only as an optional alternative or troubleshooting workflow.

## Required inputs

Each sample is expected to exist as a folder under the main results directory, for example:

```text
<results-dir>/<sample>/outs/possorted_genome_bam.bam
<results-dir>/<sample>/outs/filtered_feature_bc_matrix/barcodes.tsv
```

or, for raw mode:

```text
<results-dir>/<sample>/outs/raw_feature_bc_matrix/barcodes.tsv
```

The script also requires a reference GTF file.

## Main parameters

### `--results-dir`
Directory containing the per-sample Cell Ranger output folders.


### `--gtf-file`
Reference GTF file passed to `velocyto`.


### `--barcode-mode`
Choose one of:

- `filtered`
- `raw`

Default:

```text
filtered
```

### `--start`
Start index in the sorted sample list, inclusive.

### `--end`
End index in the sorted sample list, exclusive.

These allow the workload to be split across multiple runs.

### `--output-dir`
Optional manual override for the loom output directory.

If omitted, the script uses:

- filtered mode:
  `results_dir/velocyto_results`
- raw mode:
  `results_dir/velocyto_results/velocyto_unfiltered`

### `--velocyto-bin`
Path to the `velocyto` executable.

Default:

```text
velocyto
```

## Example usage

### Run the first 10 samples with filtered barcodes

```bash
python run_velocyto_batch.py --start 0 --end 10 --barcode-mode filtered
```

### Run the next 10 samples

```bash
python run_velocyto_batch.py --start 10 --end 20 --barcode-mode filtered
```

### Run the remaining samples

```bash
python run_velocyto_batch.py --start 20 --end 36 --barcode-mode filtered
```

### Optional raw-barcode run

```bash
python run_velocyto_batch.py --start 0 --end 10 --barcode-mode raw
```

## Expected outputs

In filtered mode, one `.loom` file per sample is written to:

```text
<results-dir>/velocyto_results/<sample>.loom
```

In raw mode, one `.loom` file per sample is written to:

```text
<results-dir>/velocyto_results/velocyto_unfiltered/<sample>.loom
```

## Notes for downstream analysis

The deposited downstream RNA velocity workflow uses the filtered-barcode loom files together with:

- a curated AnnData object
- metadata stored in the AnnData object
- dimensional reductions such as UMAP and possibly scVI embeddings

So for reproducing the main velocity figure, the filtered loom set is the important one to preserve.

## Environment notes

Before running this script, confirm that:

- `velocyto` is installed and available in your environment
- the reference GTF path is correct
- the Cell Ranger sample folders contain the expected BAM and barcode files

A simple test is:

```bash
velocyto --help
```

## Suggested deposition note

For repository deposition, this script can be described as:

> Batch script used to generate per-sample velocyto loom files from Cell Ranger outputs. The filtered-barcode mode was used for the deposited RNA velocity analysis; raw-barcode mode was retained only as an alternative workflow.
