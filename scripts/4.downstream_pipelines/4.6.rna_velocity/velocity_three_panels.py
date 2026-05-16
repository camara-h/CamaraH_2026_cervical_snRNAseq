#!/usr/bin/env python3
"""
Run RNA velocity adipocyte panels with scVelo.

This script is a standalone Python version of the velocity RMarkdown workflow.
It is intended to run locally first, then in Code Ocean with the same staged
input structure.

Expected input directory:
  velocyto_results/
    *.loom
    cell_names.txt
    scVI_embeddings.csv
    seurat_metadata.csv
    optional: cell_type_metadata.csv

Outputs:
  velocity_adipocytes_all.pdf
  velocity_adipocytes_Dsubject3_BAT.pdf
  velocity_adipocytes_Deep.pdf
  velocity_diagnostics.txt

Example local run from project root:
  conda activate scvelo_velocity
  python scripts/velocity_three_panels.py \
    --data-dir ocean_code_data/Figure_6_S5/velocyto_results \
    --out-dir output/6.figure_preparation/Figure_6_S5/velocity \
    --cache-h5ad output/6.figure_preparation/Figure_6_S5/velocity/filtered_velocity_input.h5ad \
    --use-cache

Example Code Ocean run:
  python /code/Figure_6_S5/velocity_three_panels.py \
    --data-dir /data/Figure_6_S5/velocyto_results \
    --out-dir /results/Figure_6_S5/velocity
"""

from __future__ import annotations

import argparse
import glob
import os
import random
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad
from pandas.core.arrays.categorical import CategoricalAccessor


# ---------------------------------------------------------------------
# Reproducibility and compatibility
# ---------------------------------------------------------------------

random.seed(1)
np.random.seed(1)

# Compatibility patch for newer NumPy versions used by some loom/anndata code.
if not hasattr(np, "string_"):
    np.string_ = np.bytes_
if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_

warnings.filterwarnings("ignore", category=DeprecationWarning, module="scvelo")
warnings.filterwarnings("ignore", category=FutureWarning, module="scvelo")


# ---------------------------------------------------------------------
# Plot settings
# ---------------------------------------------------------------------

ADIPOCYTE_PALETTE: Dict[str, str] = {
    "WAds": "#ffe65a",
    "BAds": "#97600f",
    "PreAds": "#a4c76d",
}


@contextmanager
def temporary_categorical_patch():
    """
    Patch used by rank_velocity_genes in some pandas versions.
    """
    original = CategoricalAccessor._delegate_property_set

    def patched(self, name, new_values):
        if name == "categories":
            self._parent = self._parent.rename_categories(new_values)
            return
        return original(self, name, new_values)

    CategoricalAccessor._delegate_property_set = patched
    try:
        yield
    finally:
        CategoricalAccessor._delegate_property_set = original


def log(msg: str, diagnostics_handle=None) -> None:
    print(msg, flush=True)
    if diagnostics_handle is not None:
        diagnostics_handle.write(str(msg) + "\n")
        diagnostics_handle.flush()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise NotADirectoryError(f"Missing {label}: {path}")


def save_intermediate_h5ad(adata_obj: ad.AnnData, cache_path: Optional[Path], diagnostics_handle=None) -> None:
    if cache_path is None:
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"Saving intermediate AnnData cache to: {cache_path}", diagnostics_handle)
    adata_obj.write_h5ad(cache_path)


def load_intermediate_h5ad(cache_path: Optional[Path], diagnostics_handle=None) -> Optional[ad.AnnData]:
    if cache_path is None:
        return None
    if not cache_path.is_file():
        log(f"Cache file requested but not found: {cache_path}", diagnostics_handle)
        return None
    log(f"Loading intermediate AnnData cache from: {cache_path}", diagnostics_handle)
    return sc.read_h5ad(cache_path)



def load_and_merge_looms(data_dir: Path, diagnostics_handle=None) -> ad.AnnData:
    """
    Load all top-level .loom files in data_dir and concatenate them.
    Cell names are transformed to match the Seurat-derived cell names.
    """
    loom_files = sorted(data_dir.glob("*.loom"))
    if not loom_files:
        raise FileNotFoundError(f"No .loom files found in {data_dir}")

    log(f"Number of loom files: {len(loom_files)}", diagnostics_handle)

    adatas = []
    keys = []

    for loom_path in loom_files:
        sample_name = loom_path.stem
        keys.append(sample_name)

        adata = sc.read_loom(str(loom_path), sparse=True)

        # Match the naming convention used in the original notebook.
        adata.obs_names = [f"{cell}_{sample_name}" for cell in adata.obs_names]

        adata.var["original_gene"] = adata.var_names
        adata.var_names_make_unique()

        log(f"Loaded {loom_path.name}: {adata.shape[0]} cells x {adata.shape[1]} genes", diagnostics_handle)
        adatas.append(adata)

    combined_adata = ad.concat(
        adatas,
        axis=0,
        join="outer",
        label="sample",
        keys=keys,
    )

    log(
        f"Combined loom object: {combined_adata.shape[0]} cells x {combined_adata.shape[1]} genes",
        diagnostics_handle,
    )

    # Original Rmd transformation:
    # name.split(":", 1)[1].replace("x_", "-1_")
    transformed_names = []
    for name in combined_adata.obs_names:
        if ":" in name:
            transformed = name.split(":", 1)[1]
        else:
            transformed = name
        transformed = transformed.replace("x_", "-1_")
        transformed_names.append(transformed)

    combined_adata.obs_names = transformed_names
    return combined_adata


def filter_to_curated_cells(
    combined_adata: ad.AnnData,
    cell_names_path: Path,
    diagnostics_handle=None,
) -> ad.AnnData:
    require_file(cell_names_path, "cell_names.txt")

    with cell_names_path.open("r") as handle:
        cell_names = [line.strip() for line in handle if line.strip()]

    filtered_adata = combined_adata[combined_adata.obs_names.isin(cell_names)].copy()

    log(
        f"After curated cell filtering: {filtered_adata.shape[0]} cells x {filtered_adata.shape[1]} genes",
        diagnostics_handle,
    )

    adata_cells = set(combined_adata.obs_names)
    cell_names_set = set(cell_names)
    common_cells = adata_cells & cell_names_set
    unique_to_adata = adata_cells - cell_names_set
    unique_to_cell_names = cell_names_set - adata_cells

    log(f"Common cells: {len(common_cells)}", diagnostics_handle)
    log(f"Unique to combined_adata: {len(unique_to_adata)}", diagnostics_handle)
    log(f"Unique to cell_names: {len(unique_to_cell_names)}", diagnostics_handle)

    if len(unique_to_cell_names) > 0:
        examples = sorted(unique_to_cell_names)[:20]
        log("Examples unique to cell_names:", diagnostics_handle)
        for example in examples:
            log(f"  {example}", diagnostics_handle)

    return filtered_adata


def add_scvi_embeddings(
    filtered_adata: ad.AnnData,
    scvi_embeddings_csv: Path,
    diagnostics_handle=None,
) -> ad.AnnData:
    require_file(scvi_embeddings_csv, "scVI_embeddings.csv")

    scvi_df = pd.read_csv(scvi_embeddings_csv, index_col=0)
    missing = filtered_adata.obs_names.difference(scvi_df.index)

    if len(missing) > 0:
        raise ValueError(
            f"scVI_embeddings.csv is missing {len(missing)} filtered cells. "
            f"First missing cells: {list(missing[:5])}"
        )

    scvi_df = scvi_df.loc[filtered_adata.obs_names]
    filtered_adata.obsm["X_scVI"] = scvi_df.values

    log(f"Added X_scVI embeddings: {filtered_adata.obsm['X_scVI'].shape}", diagnostics_handle)
    return filtered_adata


def add_seurat_metadata(
    filtered_adata: ad.AnnData,
    seurat_metadata_csv: Path,
    diagnostics_handle=None,
) -> ad.AnnData:
    require_file(seurat_metadata_csv, "seurat_metadata.csv")

    metadata = pd.read_csv(seurat_metadata_csv, index_col=0)
    missing = filtered_adata.obs_names.difference(metadata.index)

    if len(missing) > 0:
        raise ValueError(
            f"seurat_metadata.csv is missing {len(missing)} filtered cells. "
            f"First missing cells: {list(missing[:5])}"
        )

    filtered_metadata = metadata.loc[filtered_adata.obs_names].copy()
    filtered_adata.obs = filtered_metadata

    log(f"Added Seurat metadata with {filtered_adata.obs.shape[1]} columns", diagnostics_handle)
    return filtered_adata


def add_cell_type_short(
    filtered_adata: ad.AnnData,
    cell_type_metadata_csv: Optional[Path],
    diagnostics_handle=None,
) -> ad.AnnData:
    """
    Add cell_type_short.

    Preferred:
      cell_type_metadata.csv with columns barcode and cell_type_short.

    Fallback:
      Use cell_type_short already present in seurat_metadata.csv.
    """
    if cell_type_metadata_csv is not None and cell_type_metadata_csv.is_file():
        cell_type_meta = pd.read_csv(cell_type_metadata_csv, index_col="barcode")
        aligned = cell_type_meta.reindex(filtered_adata.obs_names)
        filtered_adata.obs["cell_type_short"] = aligned["cell_type_short"].astype("category")

        log(f"Added cell_type_short from {cell_type_metadata_csv}", diagnostics_handle)

        if "clean_sample" in cell_type_meta.columns:
            log(f"cell_type_metadata clean_sample count: {cell_type_meta['clean_sample'].nunique()}", diagnostics_handle)

    elif "cell_type_short" in filtered_adata.obs.columns:
        filtered_adata.obs["cell_type_short"] = filtered_adata.obs["cell_type_short"].astype("category")
        log("Using cell_type_short already present in seurat_metadata.csv", diagnostics_handle)

    else:
        raise FileNotFoundError(
            "Could not find cell_type_short. Provide cell_type_metadata.csv in the data directory "
            "or include a cell_type_short column in seurat_metadata.csv."
        )

    log("Unique values of cell_type_short:", diagnostics_handle)
    log(str(filtered_adata.obs["cell_type_short"].value_counts(dropna=False)), diagnostics_handle)
    log(
        f"Total NaN in cell_type_short: {int(filtered_adata.obs['cell_type_short'].isna().sum())}",
        diagnostics_handle,
    )

    return filtered_adata


def check_required_metadata_columns(filtered_adata: ad.AnnData) -> None:
    required = [
        "cell_type_short",
        "cell_type_group",
        "sample",
        "neck_region",
    ]
    missing = [col for col in required if col not in filtered_adata.obs.columns]
    if missing:
        raise KeyError(
            "Missing required metadata columns in filtered_adata.obs: "
            + ", ".join(missing)
            + "\nAvailable columns include: "
            + ", ".join(map(str, filtered_adata.obs.columns[:50]))
        )


def preprocess_velocity(
    filtered_adata: ad.AnnData,
    diagnostics_handle=None,
    preprocessing_mode: str = "auto",
) -> ad.AnnData:
    log("Running scVelo preprocessing on filtered AnnData...", diagnostics_handle)

    scv.pl.proportions(filtered_adata, show=False)

    if preprocessing_mode not in {"auto", "current", "legacy"}:
        raise ValueError("preprocessing_mode must be one of: auto, current, legacy")

    has_legacy = hasattr(scv.pp, "filter_genes_dispersion")

    if preprocessing_mode == "legacy" and not has_legacy:
        available = [name for name in dir(scv.pp) if "filter" in name or "normalize" in name]
        raise AttributeError(
            "Requested legacy preprocessing, but scv.pp.filter_genes_dispersion "
            "is not available in this scVelo version. Available relevant scv.pp functions: "
            + ", ".join(sorted(available))
        )

    if preprocessing_mode == "legacy" or (preprocessing_mode == "auto" and has_legacy):
        log("Using legacy preprocessing: filter_genes + normalize_per_cell + filter_genes_dispersion + log1p", diagnostics_handle)
        scv.pp.filter_genes(filtered_adata, min_shared_counts=20)
        scv.pp.normalize_per_cell(filtered_adata)
        # sc.pp.filter_genes_dispersion(filtered_adata, n_top_genes=10000)
        sc.pp.log1p(filtered_adata)
    else:
        log("Using compatibility preprocessing: filter_genes + normalize_per_cell + filter_genes_dispersion + scanpy.log1p", diagnostics_handle)
        
        scv.pp.filter_genes(filtered_adata, min_shared_counts=20)
        scv.pp.normalize_per_cell(filtered_adata)
        sc.pp.filter_genes_dispersion(
            filtered_adata,
            n_top_genes=10000,
            subset=True,
        )
        sc.pp.log1p(filtered_adata)

    sc.pp.neighbors(filtered_adata, n_neighbors=30, use_rep="X_scVI", n_pcs=None, random_state=0)
    sc.tl.umap(filtered_adata, random_state=0)

    scv.pp.moments(filtered_adata)
    scv.tl.velocity(filtered_adata)
    scv.tl.velocity_graph(filtered_adata)

    log("Finished full-object velocity preprocessing.", diagnostics_handle)

    counts = filtered_adata.obs["cell_type_group"].value_counts(dropna=False)
    log("cell_type_group counts:", diagnostics_handle)
    log(str(counts), diagnostics_handle)

    return filtered_adata


def plot_velocity_subset(
    adata_full: ad.AnnData,
    subset_mask,
    title: str,
    outfile: Path,
    color_key: str = "cell_type_short",
    palette: Dict[str, str] = ADIPOCYTE_PALETTE,
    n_neighbors: int = 15,
    diagnostics_handle=None,
) -> Optional[ad.AnnData]:
    """
    Create a subset, recompute neighbors/UMAP/moments/velocity graph,
    and save an RNA velocity plot.
    """
    sub = adata_full[subset_mask, :].copy()
    log(f"[{title}] initial subset: {sub.shape[0]} cells", diagnostics_handle)

    if sub.shape[0] < 10:
        log(f"[{title}] subset too small, skipping.", diagnostics_handle)
        return None

    if color_key in sub.obs.columns and palette is not None:
        keep = sub.obs[color_key].isin(list(palette.keys()))
        n_drop = int((~keep).sum())
        if n_drop > 0:
            log(f"[{title}] dropping {n_drop} cells with {color_key} outside palette", diagnostics_handle)
        sub = sub[keep, :].copy()

    log(f"[{title}] subset after palette filter: {sub.shape[0]} cells", diagnostics_handle)

    if sub.shape[0] < 10:
        log(f"[{title}] subset too small after palette filter, skipping.", diagnostics_handle)
        return None

    sub.obs[color_key] = sub.obs[color_key].astype("category")
    sub.obs[color_key] = sub.obs[color_key].cat.remove_unused_categories()

    uns_color_key = f"{color_key}_colors"
    if uns_color_key in sub.uns:
        del sub.uns[uns_color_key]

    n_neigh = min(n_neighbors, max(2, sub.shape[0] - 1))

    sc.pp.neighbors(
        sub,
        use_rep="X_scVI",
        n_neighbors=n_neigh,
        n_pcs=30,
        random_state=0,
    )

    try:
        sc.tl.leiden(
            sub,
            key_added="scVI_clusters",
            flavor="igraph",
            n_iterations=2,
            directed=False,
            random_state=0,
        )
    except TypeError:
        sc.tl.leiden(sub, key_added="scVI_clusters", random_state=0)

    sc.tl.umap(sub, random_state=0)

    scv.pp.moments(sub)
    scv.tl.velocity_graph(sub)

    with temporary_categorical_patch():
        try:
            scv.tl.rank_velocity_genes(sub, min_corr=0.3)
        except Exception as exc:
            log(f"[{title}] rank_velocity_genes failed, non-critical: {exc}", diagnostics_handle)

    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=600)

    cats_in_order = list(sub.obs[color_key].cat.categories)
    palette_list = [palette.get(c, "#cccccc") for c in cats_in_order]

    sc.pl.umap(
        sub,
        color=color_key,
        palette=palette_list,
        legend_loc="right margin",
        title=title,
        show=False,
        ax=ax,
    )

    scv.pl.velocity_embedding_grid(
        sub,
        basis="umap",
        color=color_key,
        palette=palette_list,
        ax=ax,
        arrow_length=5,
        show=False,
        density=1,
        linewidth=0.5,
        arrow_size=1,
        smooth=1.5,
        min_mass=10,
        arrow_color="black",
        alpha=0.6,
    )

    for cell_type in sub.obs[color_key].unique():
        mask = sub.obs[color_key] == cell_type
        umap_coords = sub.obsm["X_umap"][mask.values]
        if umap_coords.size == 0:
            continue

        x_centroid = np.median(umap_coords[:, 0])
        y_centroid = np.median(umap_coords[:, 1])

        ax.annotate(
            str(cell_type),
            (x_centroid, y_centroid),
            fontweight="bold",
            fontsize=10,
            color="black",
            ha="center",
            va="center",
            path_effects=[
                path_effects.withStroke(linewidth=3, foreground="white")
            ],
        )

    fig.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, format="pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    log(f"Saved to: {outfile}", diagnostics_handle)
    return sub


def write_subset_diagnostics(filtered_adata: ad.AnnData, mask_short, diagnostics_handle=None) -> None:
    adipo_obs = filtered_adata.obs.loc[mask_short]

    log(">>> Unique values of sample in WAds/PreAds/BAds cells:", diagnostics_handle)
    log(str(adipo_obs["sample"].value_counts(dropna=False)), diagnostics_handle)
    log("", diagnostics_handle)

    log(">>> Unique values of neck_region in WAds/PreAds/BAds cells:", diagnostics_handle)
    log(str(adipo_obs["neck_region"].value_counts(dropna=False)), diagnostics_handle)
    log("", diagnostics_handle)

    log(">>> Cross-tab sample x neck_region, adipocytes only:", diagnostics_handle)
    log(str(pd.crosstab(adipo_obs["sample"], adipo_obs["neck_region"], dropna=False)), diagnostics_handle)
    log("", diagnostics_handle)

    mask_sample = filtered_adata.obs["sample"] == "Dsubject3_BAT"
    mask_deep = filtered_adata.obs["neck_region"] == "Deep"

    mask_plot2 = mask_short & mask_sample
    mask_plot3 = mask_short & mask_deep

    n2 = int(mask_plot2.sum())
    n3 = int(mask_plot3.sum())
    n_overlap = int((mask_plot2 & mask_plot3).sum())
    n_only_2 = int((mask_plot2 & ~mask_plot3).sum())
    n_only_3 = int((~mask_plot2 & mask_plot3).sum())

    log(f"Plot 2 Dsubject3_BAT: {n2} cells", diagnostics_handle)
    log(f"Plot 3 Deep:          {n3} cells", diagnostics_handle)
    log(f"  Overlap:            {n_overlap}", diagnostics_handle)
    log(f"  Only Plot 2:         {n_only_2}", diagnostics_handle)
    log(f"  Only Plot 3:         {n_only_3}", diagnostics_handle)


def run_velocity(
    data_dir: Path,
    out_dir: Path,
    cell_type_metadata: Optional[Path],
    cache_h5ad: Optional[Path],
    use_cache: bool,
    preprocessing_mode: str,
) -> None:
    require_dir(data_dir, "velocity data directory")
    out_dir.mkdir(parents=True, exist_ok=True)

    diagnostics_path = out_dir / "velocity_diagnostics.txt"

    with diagnostics_path.open("w") as diag:
        log(f"Python: {sys.version}", diag)
        log(f"scanpy: {sc.__version__}", diag)
        log(f"scvelo: {scv.__version__}", diag)
        log(f"Data directory: {data_dir}", diag)
        log(f"Output directory: {out_dir}", diag)
        log("", diag)

        cell_names_path = data_dir / "cell_names.txt"
        scvi_embeddings_csv = data_dir / "scVI_embeddings.csv"
        seurat_metadata_csv = data_dir / "seurat_metadata.csv"

        if cell_type_metadata is None:
            candidate = data_dir / "cell_type_metadata.csv"
            cell_type_metadata = candidate if candidate.is_file() else None

        filtered_adata = None

        if use_cache:
            filtered_adata = load_intermediate_h5ad(cache_h5ad, diag)

        if filtered_adata is None:
            combined_adata = load_and_merge_looms(data_dir, diag)
            filtered_adata = filter_to_curated_cells(combined_adata, cell_names_path, diag)
            filtered_adata = add_scvi_embeddings(filtered_adata, scvi_embeddings_csv, diag)
            filtered_adata = add_seurat_metadata(filtered_adata, seurat_metadata_csv, diag)
            filtered_adata = add_cell_type_short(filtered_adata, cell_type_metadata, diag)
            save_intermediate_h5ad(filtered_adata, cache_h5ad, diag)

        check_required_metadata_columns(filtered_adata)

        filtered_adata = preprocess_velocity(
            filtered_adata,
            diagnostics_handle=diag,
            preprocessing_mode=preprocessing_mode,
        )

        selected_short = ["WAds", "PreAds", "BAds"]
        mask_short = filtered_adata.obs["cell_type_short"].isin(selected_short)

        write_subset_diagnostics(filtered_adata, mask_short, diag)

        plot_velocity_subset(
            filtered_adata,
            subset_mask=mask_short,
            title="RNA Velocity - Adipocytes (all samples)",
            outfile=out_dir / "velocity_adipocytes_all.pdf",
            diagnostics_handle=diag,
        )

        mask_sample = filtered_adata.obs["sample"] == "Dsubject3_BAT"
        mask_plot2 = mask_short & mask_sample

        log(f"Cells in Dsubject3_BAT subset: {int(mask_plot2.sum())}", diag)

        plot_velocity_subset(
            filtered_adata,
            subset_mask=mask_plot2,
            title="RNA Velocity - Adipocytes (Dsubject3_BAT)",
            outfile=out_dir / "velocity_adipocytes_Dsubject3_BAT.pdf",
            diagnostics_handle=diag,
        )

        mask_deep = filtered_adata.obs["neck_region"] == "Deep"
        mask_plot3 = mask_short & mask_deep

        log(f"Cells in Deep subset: {int(mask_plot3.sum())}", diag)

        plot_velocity_subset(
            filtered_adata,
            subset_mask=mask_plot3,
            title="RNA Velocity - Adipocytes (Deep neck region)",
            outfile=out_dir / "velocity_adipocytes_Deep.pdf",
            diagnostics_handle=diag,
        )

        # Final cross-tab diagnostic, if available.
        if "cell_type_res.1.6" in filtered_adata.obs.columns:
            log("", diag)
            log(">>> cell_type_short x cell_type_res.1.6 diagnostic:", diag)
            ct = pd.crosstab(
                filtered_adata.obs["cell_type_short"],
                filtered_adata.obs["cell_type_res.1.6"],
                dropna=False,
            )
            mask_rows = ct.index.isin(["WAds", "BAds", "PreAds"])
            mask_cols = (ct.loc[mask_rows].sum(axis=0) > 0)
            log(str(ct.loc[mask_rows, mask_cols]), diag)

    print(f"\nDone. Diagnostics written to: {diagnostics_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run scVelo velocity panels for Figure 6/S5."
    )

    parser.add_argument(
        "--data-dir",
        required=True,
        type=Path,
        help="Directory containing .loom files, cell_names.txt, scVI_embeddings.csv, and seurat_metadata.csv.",
    )

    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Directory where velocity PDFs and diagnostics will be written.",
    )

    parser.add_argument(
        "--cell-type-metadata",
        default=None,
        type=Path,
        help=(
            "Optional CSV with barcode and cell_type_short columns. "
            "If omitted, the script first looks for data-dir/cell_type_metadata.csv, "
            "then falls back to cell_type_short in seurat_metadata.csv."
        ),
    )

    parser.add_argument(
        "--cache-h5ad",
        default=None,
        type=Path,
        help=(
            "Optional path to save or read a filtered intermediate AnnData .h5ad file. "
            "This avoids reloading and concatenating all loom files after the first successful cache creation."
        ),
    )

    parser.add_argument(
        "--use-cache",
        action="store_true",
        help="If provided, use --cache-h5ad when it already exists.",
    )

    parser.add_argument(
        "--preprocessing-mode",
        default="auto",
        choices=["auto", "current", "legacy"],
        help=(
            "auto uses legacy preprocessing only if filter_genes_dispersion exists; "
            "current uses scv.pp.filter_and_normalize; "
            "legacy forces the old filter_genes_dispersion workflow."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_velocity(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        cell_type_metadata=args.cell_type_metadata,
        cache_h5ad=args.cache_h5ad,
        use_cache=args.use_cache,
        preprocessing_mode=args.preprocessing_mode,
    )


if __name__ == "__main__":
    main()
