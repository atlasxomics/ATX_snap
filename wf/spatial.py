import anndata
import logging
from pathlib import Path

import wf.plotting as pl


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    adata.obsm["spatial"] = adata.obs[[y_key, x_key]].values

    return adata


def run_squidpy_analysis(
    adata_gene: anndata.AnnData, figures_dir: Path
) -> anndata.AnnData:
    """Run Squidpy analysis and generate plots."""
    from squidpy.pl import ripley

    logging.info("Running squidpy...")
    adata_gene = squidpy_analysis(adata_gene)

    # Generate neighborhood plots
    logging.info("Making neighborhood plots...")
    group_dict = {"all": None}

    for group_name, group_value in group_dict.items():
        pl.plot_neighborhoods(
            adata_gene, group_name, group_value, outdir=str(figures_dir)
        )

    # Generate Ripley's plot
    logging.info("Running ripley...")
    ripley(adata_gene, cluster_key="cluster", mode="L", save="ripleys_L.pdf")

    return adata_gene


def squidpy_analysis(
    adata: anndata.AnnData, cluster_key: str = "cluster"
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, ripley, spatial_neighbors

    if not adata.obs["cluster"].dtype.name == "category":
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")

    spatial_neighbors(adata, coord_type="grid", n_neighs=4, n_rings=1)
    nhood_enrichment(adata, cluster_key=cluster_key)
    ripley(adata, cluster_key="cluster", mode="L", max_dist=500)

    return adata
