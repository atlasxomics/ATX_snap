import anndata
import squidpy as sq


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    adata.obsm["spatial"] = adata.obs[[y_key, x_key]].values

    return adata


def squidpy_analysis(
    adata: anndata.AnnData, cluster_key: str = "cluster"
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """

    if not adata.obs["cluster"].dtype.name == "category":
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")

    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=4, n_rings=1)
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    sq.gr.ripley(adata, cluster_key="cluster", mode="L", max_dist=500)

    return adata
