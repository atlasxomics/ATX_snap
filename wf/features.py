import logging
import gc
import glob

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import numpy as np
import scanpy as sc
import scipy.sparse as sp

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def clean_adata(adata: anndata.AnnData) -> anndata.AnnData:
    obs = [
        "barcode", "n_genes_by_counts", "log1p_n_genes_by_counts", "total_counts",
        "log1p_total_counts", "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
        "pct_counts_in_top_200_genes", "pct_counts_in_top_500_genes", "total_counts_mt",
        "log1p_total_counts_mt", "pct_counts_mt"
    ]
    var = [
        "mt", "n_counts", "n_cells", "n_cells_by_counts", "mean_counts", "log1p_mean_counts",
        "pct_dropout_by_counts", "total_counts", "log1p_total_counts", "means",
        "dispersions", "dispersions_norm"
    ]

    rm_obs = [o for o in obs if o in adata.obs.keys()]
    rm_var = [v for v in var if v in adata.var.keys()]

    adata.obs.drop(rm_obs, axis=1, inplace=True)
    adata.var.drop(rm_var, axis=1, inplace=True)

    adata.varm.clear()
    adata.layers.clear()

    if adata.raw:
        adata.raw = None

    for uns in ["pca", "log1p"]:
        adata.uns.pop(uns, None)

    rm_obsm = [obsm for obsm in adata.obsm.keys()
               if obsm not in ["spatial", "X_umap"]]

    for obsm in rm_obsm:
        del adata.obsm[obsm]

    # Sparse to dense conversion
    if sp.issparse(adata.X):
        try:
            adata.X = adata.X.toarray()
        except Exception as e:
            logging.warning(f"Could not convert sparse .X to dense: {e}")

    try:
        adata.X = adata.X.astype(np.float16)
    except Exception as e:
        logging.warning(f"Cannot convert .X to float16: {e}")

    return adata


def clean_index_columns(*adatas: anndata.AnnData) -> None:
    """Remove _index columns from AnnData objects if they exist."""
    for adata in adatas:
        if hasattr(adata, 'raw') and adata.raw is not None:
            if "_index" in adata.raw.var.columns:
                adata.raw.var.drop(columns=['_index'], inplace=True)


def load_and_combine_data(suffix: str) -> anndata.AnnData:
    """Load and combine AnnData objects."""
    logging.info("Reading and combining gene AnnData...")
    files = glob.glob(f"*{suffix}.h5ad")
    adatas = [anndata.read_h5ad(file) for file in files]
    adata = sc.concat(adatas)

    # Clean up memory
    del adatas
    gc.collect()

    # Clean up index columns if they exist
    clean_index_columns(adata)

    return adata


def load_analysis_results(
    adata: anndata.AnnData,
    type: str,
    groups: List[str]
) -> None:
    """Load various analysis results into AnnData objects."""
    # Load differential analysis results
    if type == "gene":
        load_ranked_genes(adata)
    elif type == "motif":
        load_enriched_motifs(adata)

    # Load heatmap data
    load_heatmaps(adata, type)

    # Load volcano plots if condition analysis was performed
    if "condition" in groups:
        load_volcano_plots(adata, type)


def load_csv_files_to_uns(
    pattern: str,
    target_uns: Dict,
    dtype_spec: Optional[Dict] = None,
    index_col: Optional[int] = None,
    name_transform: Optional[callable] = None
) -> None:
    import pandas as pd
    """Generic function to load CSV files matching a pattern into uns
    dictionary."""
    try:
        files = glob.glob(pattern)
        for file in files:
            name = Path(file).stem
            if name_transform:
                name = name_transform(name, file)
            df = pd.read_csv(file, dtype=dtype_spec, index_col=index_col)
            target_uns[name] = df
    except Exception as e:
        logging.warning(f"Error loading files matching {pattern}: {e}")


def load_enriched_motifs(adata_motif: anndata.AnnData) -> None:
    """Load enriched motifs data."""
    logging.info("Adding enriched motifs...")
    load_csv_files_to_uns(
        "enrichedMotifs_*.csv", adata_motif.uns, dtype_spec={"group_name": str}
    )


def load_heatmaps(
    adata: anndata.AnnData,
    type: str
) -> None:
    """Load heatmap data for both genes and motifs."""
    logging.info("Adding heatmaps...")

    if type == "gene":
        # Gene heatmaps
        load_csv_files_to_uns(
            "genes_per_*_hm.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            index_col=0
        )

    elif type == "motif":
    # Motif heatmaps
        load_csv_files_to_uns(
            "motif_per_*_hm.csv",
            adata.uns,
            index_col=0
        )


def load_ranked_genes(adata_gene: anndata.AnnData) -> None:
    """Load ranked genes data."""
    logging.info("Adding ranked genes...")
    load_csv_files_to_uns(
        "ranked_genes_*.csv", adata_gene.uns, dtype_spec={"group_name": str}
    )


def load_volcano_plots(
    adata: anndata.AnnData, type: str
) -> None:
    """Load volcano plot data for genes and motifs."""
    if type == "gene":
        logging.info("Adding gene volcanos...")
        load_csv_files_to_uns(
            "volcanoMarkers_genes_*.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            name_transform=lambda name, file: volcano_name_transform(name, file, "genes")
        )

    elif type == "motif":
        logging.info("Adding motif volcanos...")
        load_csv_files_to_uns(
            "volcanoMarkers_motifs_*.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            name_transform=lambda name,
            file: volcano_name_transform(name, file, "motifs")
        )


def save_anndata_objects(
    adata: anndata.AnnData,
    suffix: str,
    base_dir: Path
) -> None:
    """Save full and reduced AnnData objects."""
    logging.info("Saving full adata...")
    # Save full objects
    adata.write(f"{base_dir}/combined{suffix}.h5ad")

    # Create and save reduced objects
    logging.info("Making reduced adata...")
    sm_adata = clean_adata(adata)

    logging.info("Saving reduced adata...")
    sm_adata.write(f"{base_dir}/combined_sm{suffix}.h5ad")


def transfer_auxiliary_data(
    adata: anndata.AnnData,
    data_paths: Dict[str, str],
    groups: List[str]
) -> None:
    """Transfer observation, UMAP, and spatial data to AnnData objects."""
    logging.info("Transferring auxiliary data...")

    # Transfer observation data
    transfer_obs_data(adata, data_paths['obs'], groups)

    # Transfer UMAP coordinates
    transfer_embedding_data(
        adata, data_paths['umap'], 'X_umap'
    )

    # Transfer spatial coordinates
    transfer_embedding_data(
        adata, data_paths['spatial'], 'spatial'
    )


def transfer_embedding_data(
    adata: anndata.AnnData,
    data_path: str,
    obsm_key: str
) -> None:
    """Transfer embedding data (UMAP/spatial) to AnnData objects."""
    import pandas as pd

    logging.info(f"Transferring {obsm_key}...")

    try:
        df = pd.read_csv(data_path, index_col=0)
        aligned_data_g = df.loc[adata_gene.obs_names].values
        aligned_data_m = df.loc[adata_motif.obs_names].values

        adata_gene.obsm[obsm_key] = aligned_data_g
        adata_motif.obsm[obsm_key] = aligned_data_m

    except (FileNotFoundError, KeyError) as e:
        logging.warning(f"Error loading {obsm_key} data: {e}")


def transfer_obs_data(
    adata: anndata.AnnData,
    obs_path: str,
    groups: List[str]
) -> None:
    """Transfer observation data to AnnData objects."""
    import pandas as pd

    logging.info("Transferring observation data...")
    try:
        obs = pd.read_csv(obs_path, index_col=0)
        if obs.empty:
            return

        obs_aligned_g = obs.reindex(adata_gene.obs.index)
        adata_gene.obs = obs_aligned_g

        obs_aligned_m = obs.reindex(adata_motif.obs.index)
        adata_motif.obs = obs_aligned_m

        # Ensure group columns are strings
        for group in groups:
            if group in adata.obs.columns and adata.obs[group].dtype != object:
                adata.obs[group] = adata.obs[group].astype(str)

    except FileNotFoundError:
        logging.warning(f"File not found: {obs_path}")


def volcano_name_transform(name: str, file: str, data_type: str) -> str:
    """Transform volcano plot file names."""
    file_name = Path(file).name
    treatment = file_name.replace(f"volcanoMarkers_{data_type}_", "").replace(".csv", "")
    return f"volcano_{treatment}"
