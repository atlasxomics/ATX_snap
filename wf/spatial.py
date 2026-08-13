import anndata
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Iterable, Optional

import wf.plotting as pl


GROUPED_NHOOD_ENRICHMENT_KEY = "cluster_nhood_enrichment_by_group"
GROUPED_NHOOD_SCHEMA_VERSION = 1


def _ensure_obs_categorical(adata: anndata.AnnData, key: Optional[str]) -> None:
    if key is None or key not in adata.obs.columns:
        return

    if adata.obs[key].dtype.name == "category":
        adata.obs[key] = adata.obs[key].cat.remove_unused_categories()
    else:
        adata.obs[key] = adata.obs[key].astype("category")


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    # Backed SnapATAC2 exposes `.obs` as a PyDataFrameElem. Its Rust-backed
    # implementation does not support pandas-style multi-column slicing, so
    # materialize the small observation table before selecting coordinates.
    if hasattr(adata.obs, "to_pandas"):
        obs = adata.obs.to_pandas()
    else:
        obs = adata.obs
    spatial = obs[[y_key, x_key]].to_numpy(copy=True)

    # Negate row (y) so the coordinate follows Plotly convention (y increases upward).
    # squidpy's spatial_scatter inverts y internally, so static plots must reverse
    # this negation before calling squidpy (see plotting.plot_spatial).
    spatial[:, 1] *= -1
    adata.obsm["spatial"] = spatial

    return adata


def run_squidpy_analysis(
    adata: anndata.AnnData,
    figures_dir: Path,
    sample_key: Optional[str] = None,
    group_keys: Optional[Iterable[str]] = None,
) -> anndata.AnnData:
    """Run Squidpy analysis, including results consumed by backed Plots."""

    logging.info("Running squidpy...")
    adata = squidpy_analysis(adata, sample_key=sample_key)
    precompute_grouped_nhood_enrichment(
        adata,
        group_keys=group_keys or (),
        sample_key=sample_key,
    )

    # Generate neighborhood plots
    logging.info("Making neighborhood plots...")
    group_dict = {"all": None}

    for group_name, group_value in group_dict.items():
        pl.plot_neighborhoods(
            adata, group_name, group_value, outdir=str(figures_dir)
        )

    return adata


def _nhood_result_for_subgroup(
    adata: anndata.AnnData,
    mask: np.ndarray,
    cluster_key: str,
    sample_key: Optional[str],
    spatial_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute enrichment without copying the feature matrix from ``adata``."""
    obs_keys = [cluster_key]
    if sample_key is not None and sample_key != cluster_key:
        obs_keys.append(sample_key)

    subset = anndata.AnnData(obs=adata.obs.loc[mask, obs_keys].copy())
    subset.obsm[spatial_key] = np.asarray(adata.obsm[spatial_key])[mask].copy()
    _ensure_obs_categorical(subset, cluster_key)
    _ensure_obs_categorical(subset, sample_key)

    categories = subset.obs[cluster_key].cat.categories.astype(str).to_numpy()
    if subset.n_obs < 2:
        shape = (len(categories), len(categories))
        return categories, np.zeros(shape, dtype=float), np.zeros(shape, dtype=float)

    squidpy_analysis(
        subset,
        cluster_key=cluster_key,
        sample_key=sample_key,
        spatial_key=spatial_key,
    )
    result = subset.uns[f"{cluster_key}_nhood_enrichment"]
    return categories, np.asarray(result["zscore"]), np.asarray(result["count"])


def precompute_grouped_nhood_enrichment(
    adata: anndata.AnnData,
    group_keys: Iterable[str],
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None,
    spatial_key: Optional[str] = None,
) -> None:
    """Store subgroup neighborhood results in a backed-AnnData-friendly schema.

    The result contains only small enrichment matrices and their labels.  It is
    deliberately indexed with numeric HDF5-safe keys because user-provided group
    values can contain slashes and other characters that are unsafe as HDF5 paths.
    """
    if spatial_key is None:
        spatial_key = next(
            (key for key in ("spatial_offset", "spatial") if key in adata.obsm),
            None,
        )
    if spatial_key is None:
        raise KeyError(
            "Spatial coordinates were not found in `adata.obsm`; expected "
            "`spatial_offset` or `spatial`."
        )

    stored_groups = {}
    seen = set()
    for group_key in group_keys:
        if group_key in seen or group_key == cluster_key:
            continue
        seen.add(group_key)
        if group_key not in adata.obs:
            logging.warning(
                "Skipping neighborhood precomputation for missing obs key '%s'.",
                group_key,
            )
            continue

        values = pd.unique(adata.obs[group_key].dropna())
        stored_subgroups = {}
        for subgroup_index, group_value in enumerate(values):
            mask = (adata.obs[group_key] == group_value).to_numpy()
            logging.info(
                "Precomputing neighborhood enrichment for %s=%s (%d cells).",
                group_key,
                group_value,
                int(mask.sum()),
            )
            categories, zscore, count = _nhood_result_for_subgroup(
                adata,
                mask,
                cluster_key,
                sample_key,
                spatial_key,
            )
            stored_subgroups[str(subgroup_index)] = {
                "group_value": str(group_value),
                "cluster_categories": categories,
                "zscore": zscore,
                "count": count,
            }

        stored_groups[str(len(stored_groups))] = {
            "group_key": str(group_key),
            "subgroups": stored_subgroups,
        }

    adata.uns[GROUPED_NHOOD_ENRICHMENT_KEY] = {
        "schema_version": GROUPED_NHOOD_SCHEMA_VERSION,
        "cluster_key": cluster_key,
        "groups": stored_groups,
    }


def squidpy_analysis(
    adata: anndata.AnnData,
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None,
    spatial_key: Optional[str] = None,
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, spatial_neighbors

    _ensure_obs_categorical(adata, cluster_key)
    _ensure_obs_categorical(adata, sample_key)

    n_clusters = len(adata.obs[cluster_key].cat.categories)
    if spatial_key is None:
        spatial_key = next(
            (key for key in ("spatial_offset", "spatial") if key in adata.obsm),
            None,
        )
    if spatial_key is None:
        raise KeyError(
            "Spatial coordinates were not found in `adata.obsm`; expected "
            "`spatial_offset` or `spatial`."
        )

    spatial_neighbors(
        adata,
        spatial_key=spatial_key,
        coord_type="grid",
        n_neighs=4,
        n_rings=1,
        library_key=sample_key,
    )
    if n_clusters < 2:
        logging.warning(
            "Skipping Squidpy neighborhood enrichment because only "
            f"{n_clusters} cluster is present."
        )
        adata.uns[f"{cluster_key}_nhood_enrichment"] = {
            "zscore": np.zeros((n_clusters, n_clusters), dtype=float),
            "count": np.zeros((n_clusters, n_clusters), dtype=float),
            "skipped": "fewer_than_two_clusters",
        }
        return adata

    nhood_enrichment(
        adata, cluster_key=cluster_key, library_key=sample_key, seed=42
    )

    return adata


def run_spatial_autocorr(
    adata: anndata.AnnData,
    n_jobs: int = 4,
) -> pd.DataFrame:
    """Compute Moran's I per feature using squidpy. Returns DataFrame sorted by I descending."""
    import squidpy as sq

    sample_key = "sample" if "sample" in adata.obs.columns else None
    _ensure_obs_categorical(adata, sample_key)
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(
            adata, coord_type="grid", n_neighs=4, n_rings=1, library_key=sample_key
        )

    layer = None
    for candidate in ["log1p", "lognorm", "normalized"]:
        if candidate in adata.layers:
            layer = candidate
            break

    features_upper = pd.Index(adata.var_names.astype(str)).str.upper()
    keep = ~(
        features_upper.str.startswith("MT-")
        | features_upper.str.startswith("RPS")
        | features_upper.str.startswith("RPL")
        | features_upper.str.startswith("MTRNR")
    )
    if "highly_variable" in adata.var.columns:
        keep = keep & adata.var["highly_variable"].to_numpy()

    test_features = adata.var_names[keep].tolist()
    if not test_features:
        raise ValueError("No features remain after filtering for spatial autocorrelation.")

    logging.info(
        "Running spatial autocorrelation (Moran's I) on %d features (layer=%s).",
        len(test_features),
        layer or "X",
    )
    sq.gr.spatial_autocorr(
        adata, mode="moran", genes=test_features, layer=layer, n_perms=None, n_jobs=n_jobs
    )
    return adata.uns["moranI"].sort_values("I", ascending=False)
