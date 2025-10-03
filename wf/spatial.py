import anndata
import logging
from pathlib import Path
from typing import Optional

import wf.plotting as pl


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    adata.obsm["spatial"] = adata.obs[[y_key, x_key]].values

    return adata


def run_squidpy_analysis(
    adata: anndata.AnnData, figures_dir: Path, sample_key: Optional[str] = None
) -> anndata.AnnData:
    """Run Squidpy analysis and generate plots."""

    logging.info("Running squidpy...")
    adata = squidpy_analysis(adata, sample_key=sample_key)

    # Generate neighborhood plots
    logging.info("Making neighborhood plots...")
    group_dict = {"all": None}

    for group_name, group_value in group_dict.items():
        pl.plot_neighborhoods(
            adata, group_name, group_value, outdir=str(figures_dir)
        )

    return adata


def squidpy_analysis(
    adata: anndata.AnnData,
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, spatial_neighbors

    if not adata.obs[cluster_key].dtype.name == "category":
        adata.obs[cluster_key] = adata.obs["cluster"].astype("category")

    if sample_key:
        if not adata.obs[sample_key].dtype.name == "category":
            adata.obs[sample_key] = adata.obs[sample_key].astype("category")

    spatial_neighbors(
        adata, coord_type="grid", n_neighs=4, n_rings=1, library_key=sample_key
    )
    nhood_enrichment(
        adata, cluster_key=cluster_key, library_key=sample_key, seed=42
    )

    return adata
