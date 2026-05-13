CellRanger and CellRanger ARC Snakemake template

What this is
A small Snakemake workflow that runs CellRanger count for GEX only libraries and CellRanger ARC count for paired GEX plus ATAC libraries.
It is designed for HPC use and for clean reruns. Inputs are described in config.yaml, not hard coded in the Snakefile.

Folder layout you will get after runs
results/
  cellranger/
    SAMPLE_ID/
      outs/...
  cellranger_arc/
    ARC_SAMPLE_ID/
      outs/...

Files you edit
config.yaml
  Add your genome references and add one entry per sample.
  For ARC samples, define two libraries per sample, one Gene Expression and one Chromatin Accessibility.

How to run locally
snakemake -j 1 -p

How to run on a cluster
Use your existing Snakemake cluster setup. This template includes a simple example slurm profile folder.
If you already have a profile, just point Snakemake to it.

Notes
1. CellRanger creates a lot of intermediate files. Point run directories to a fast filesystem such as scratch.
2. Make sure your reference paths are readable by compute nodes.
3. Validate fastq paths and sample names with a quick file listing before you launch many samples.
