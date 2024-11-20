import anndata
import pandas as pd


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    xcor = adata.obs[x_key]
    ycor = adata.obs[y_key]

    spatial_df = pd.DataFrame({xcor.name: xcor, ycor.name: ycor})
    adata.obsm["spatial"] = spatial_df

    return adata
