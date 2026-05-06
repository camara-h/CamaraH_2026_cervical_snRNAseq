#!/usr/bin/env python3
"""
Run velocyto on a batch of Cell Ranger sample folders.

This consolidated script replaces the older split scripts:
- velocyto1.py / velocyto2.py / velocyto3.py
- velocyto4.py / velocyto5.py / velocyto6.py

It supports:
- filtered barcodes output, which matches the downstream velocity analysis
- raw or unfiltered barcodes output for optional troubleshooting
- configurable sample ranges so jobs can still be split into batches

Example usage
-------------

# Main workflow, first 10 samples, filtered barcodes
python run_velocyto_batch.py --start 0 --end 10 --barcode-mode filtered

# Main workflow, samples 10 to 20
python run_velocyto_batch.py --start 10 --end 20 --barcode-mode filtered

# Optional raw-barcode workflow
python run_velocyto_batch.py --start 0 --end 10 --barcode-mode raw
"""

import argparse
import subprocess
from pathlib import Path


def find_file(filepath: Path):
    """Return filepath if it exists, otherwise try filepath.gz."""
    if filepath.exists():
        return filepath
    gz_path = Path(str(filepath) + ".gz")
    if gz_path.exists():
        return gz_path
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run velocyto over a selected range of Cell Ranger result folders."
    )

    parser.add_argument(
        "--results-dir",
        help="Directory containing per-sample Cell Ranger output folders.",
    )
    parser.add_argument(
        "--gtf-file",
        help="Reference GTF file used by velocyto.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Optional override for output directory. "
            "If omitted, this is derived from --barcode-mode."
        ),
    )
    parser.add_argument(
        "--barcode-mode",
        choices=["filtered", "raw"],
        default="filtered",
        help=(
            "Which barcode set to use. "
            "'filtered' writes to velocyto_results/ and matches the downstream "
            "velocity workflow. 'raw' writes to velocyto_results/velocyto_unfiltered/."
        ),
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Start index in the sorted sample list, inclusive.",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End index in the sorted sample list, exclusive.",
    )
    parser.add_argument(
        "--velocyto-bin",
        default="velocyto",
        help="Path to the velocyto executable.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    results_dir = Path(args.results_dir)
    gtf_file = Path(args.gtf_file)

    if args.output_dir is not None:
        output_dir = Path(args.output_dir)
    else:
        if args.barcode_mode == "filtered":
            output_dir = results_dir / "velocyto_results"
        else:
            output_dir = results_dir / "velocyto_results" / "velocyto_unfiltered"

    output_dir.mkdir(parents=True, exist_ok=True)

    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")
    if not gtf_file.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")

    samples = sorted([d.name for d in results_dir.iterdir() if d.is_dir()])
    selected_samples = samples[args.start:args.end]

    print(f"Barcode mode: {args.barcode_mode}")
    print(f"Results dir:  {results_dir}")
    print(f"Output dir:   {output_dir}")
    print(f"GTF file:     {gtf_file}")
    print(f"Sample range: {args.start}:{args.end}")
    print(f"N samples:    {len(selected_samples)}")

    if args.barcode_mode == "filtered":
        barcode_relpath = Path("outs/filtered_feature_bc_matrix/barcodes.tsv")
    else:
        barcode_relpath = Path("outs/raw_feature_bc_matrix/barcodes.tsv")

    for sample in selected_samples:
        print(f"\nProcessing: {sample}")

        loom_file = output_dir / f"{sample}.loom"
        if loom_file.exists():
            print(f"File {loom_file} already exists. Skipping...")
            continue

        sample_dir = results_dir / sample
        bam_path = sample_dir / "outs/possorted_genome_bam.bam"
        barcodes_path = sample_dir / barcode_relpath

        bam_file = find_file(bam_path)
        barcodes_file = find_file(barcodes_path)

        if bam_file is None or barcodes_file is None:
            print("Required BAM or barcode file not found. Skipping...")
            continue

        command = [
            args.velocyto_bin,
            "run",
            "-b",
            str(barcodes_file),
            "-o",
            str(output_dir),
            "--sampleid",
            sample,
            str(bam_file),
            str(gtf_file),
        ]

        try:
            subprocess.run(command, check=True)
        except Exception as e:
            print(f"Error while processing {sample}: {e}")


if __name__ == "__main__":
    main()
