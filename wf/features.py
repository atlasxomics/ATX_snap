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


def load_and_combine_data() -> Tuple[anndata.AnnData, anndata.AnnData]:
    """Load and combine gene and motif AnnData objects."""
    logging.info("Reading and combining gene AnnData...")
    gene_files = glob.glob("*g_converted.h5ad")
    gene_adatas = [anndata.read_h5ad(file) for file in gene_files]
    adata_gene = sc.concat(gene_adatas)

    logging.info("Reading and combining motif AnnData...")
    motif_files = glob.glob("*m_converted.h5ad")
    motif_adatas = [anndata.read_h5ad(file) for file in motif_files]
    adata_motif = sc.concat(motif_adatas)

    # Clean up memory
    del gene_adatas, motif_adatas
    gc.collect()

    # Clean up index columns if they exist
    clean_index_columns(adata_gene, adata_motif)

    return adata_gene, adata_motif


def load_analysis_results(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    groups: List[str]
) -> None:
    """Load various analysis results into AnnData objects."""
    # Load differential analysis results
    load_ranked_genes(adata_gene)
    load_enriched_motifs(adata_motif)

    # Load heatmap data
    load_heatmaps(adata_gene, adata_motif)

    # Load volcano plots if condition analysis was performed
    if "condition" in groups:
        load_volcano_plots(adata_gene, adata_motif)


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
    adata_gene: anndata.AnnData, adata_motif: anndata.AnnData
) -> None:
    """Load heatmap data for both genes and motifs."""
    logging.info("Adding heatmaps...")

    # Gene heatmaps
    load_csv_files_to_uns(
        "genes_per_*_hm.csv",
        adata_gene.uns,
        dtype_spec={"cluster": str},
        index_col=0
    )

    # Motif heatmaps
    load_csv_files_to_uns(
        "motif_per_*_hm.csv",
        adata_motif.uns,
        index_col=0
    )


def load_ranked_genes(adata_gene: anndata.AnnData) -> None:
    """Load ranked genes data."""
    logging.info("Adding ranked genes...")
    load_csv_files_to_uns(
        "ranked_genes_*.csv", adata_gene.uns, dtype_spec={"group_name": str}
    )


def load_volcano_plots(
    adata_gene: anndata.AnnData, adata_motif: anndata.AnnData
) -> None:
    """Load volcano plot data for genes and motifs."""
    logging.info("Adding gene volcanos...")
    load_csv_files_to_uns(
        "volcanoMarkers_genes_*.csv",
        adata_gene.uns,
        dtype_spec={"cluster": str},
        name_transform=lambda name, file: volcano_name_transform(name, file, "genes")
    )

    logging.info("Adding motif volcanos...")
    load_csv_files_to_uns(
        "volcanoMarkers_motifs_*.csv",
        adata_motif.uns,
        dtype_spec={"cluster": str},
        name_transform=lambda name,
        file: volcano_name_transform(name, file, "motifs")
    )


def save_anndata_objects(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    base_dir: Path
) -> None:
    """Save full and reduced AnnData objects."""
    logging.info("Saving full adata...")

    # Save full objects
    adata_gene.write(f"{base_dir}/combined_ge.h5ad")
    adata_motif.write(f"{base_dir}/combined_motifs.h5ad")

    # Create and save reduced objects
    logging.info("Making reduced gene adata...")
    sm_adata_gene = clean_adata(adata_gene)
    sm_adata_motif = clean_adata(adata_motif)

    logging.info("Saving gene adata...")
    sm_adata_gene.write(f"{base_dir}/combined_sm_ge.h5ad")

    logging.info("Saving motif adata...")
    sm_adata_motif.write(f"{base_dir}/combined_sm_motifs.h5ad")


def transfer_auxiliary_data(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    data_paths: Dict[str, str],
    groups: List[str]
) -> None:
    """Transfer observation, UMAP, and spatial data to AnnData objects."""
    logging.info("Transferring auxiliary data...")

    # Transfer observation data
    transfer_obs_data(adata_gene, adata_motif, data_paths['obs'], groups)

    # Transfer UMAP coordinates
    transfer_embedding_data(
        adata_gene, adata_motif, data_paths['umap'], 'X_umap'
    )

    # Transfer spatial coordinates
    transfer_embedding_data(
        adata_gene, adata_motif, data_paths['spatial'], 'spatial'
    )


def transfer_embedding_data(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    data_path: str,
    obsm_key: str
) -> None:
    """Transfer embedding data (UMAP/spatial) to AnnData objects."""
    import pandas as pd

    logging.info(f"Transferring {obsm_key}...")

    try:
        df = pd.read_csv(data_path, index_col=0)
        aligned_data = df.loc[adata_gene.obs_names].values

        adata_gene.obsm[obsm_key] = aligned_data
        adata_motif.obsm[obsm_key] = aligned_data

    except (FileNotFoundError, KeyError) as e:
        logging.warning(f"Error loading {obsm_key} data: {e}")


def transfer_obs_data(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
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

        obs_aligned = obs.reindex(adata_gene.obs.index)
        adata_gene.obs = obs_aligned
        adata_motif.obs = obs_aligned

        # Ensure group columns are strings
        for group in groups:
            for adata in [adata_gene, adata_motif]:
                if group in adata.obs.columns and adata.obs[group].dtype != object:
                    adata.obs[group] = adata.obs[group].astype(str)

    except FileNotFoundError:
        logging.warning(f"File not found: {obs_path}")


def volcano_name_transform(name: str, file: str, data_type: str) -> str:
    """Transform volcano plot file names."""
    file_name = Path(file).name
    treatment = file_name.replace(f"volcanoMarkers_{data_type}_", "").replace(".csv", "")
    return f"volcano_{treatment}"
