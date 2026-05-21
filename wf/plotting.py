import anndata
import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc

from matplotlib.backends.backend_pdf import PdfPages
from typing import Callable, List, Optional

from wf.spatial import squidpy_analysis
from atx_common import filter_anndata


def _sanitize_nhood_enrichment(
    adata: anndata.AnnData,
    cluster_key: str = "cluster",
) -> None:
    key = f"{cluster_key}_nhood_enrichment"
    result = adata.uns.get(key)
    if not isinstance(result, dict) or "zscore" not in result:
        return

    zscore = result["zscore"]
    zscore_array = np.asarray(zscore, dtype=float)
    if np.isfinite(zscore_array).all():
        return

    logging.warning(
        "Neighborhood enrichment contains non-finite z-scores; replacing "
        "NaN/Inf with 0 for plotting."
    )
    sanitized = np.nan_to_num(zscore_array, nan=0.0, posinf=0.0, neginf=0.0)

    if isinstance(zscore, pd.DataFrame):
        result["zscore"] = pd.DataFrame(
            sanitized, index=zscore.index, columns=zscore.columns
        )
    else:
        result["zscore"] = sanitized


def _nhood_labels(adata: anndata.AnnData, cluster_key: str, n_labels: int) -> List[str]:
    if cluster_key in adata.obs:
        values = adata.obs[cluster_key]
        if hasattr(values, "cat"):
            labels = [str(label) for label in values.cat.categories]
        else:
            labels = sorted([str(label) for label in values.dropna().unique()])

        if len(labels) == n_labels:
            return labels

    return [str(i) for i in range(n_labels)]


def _fallback_nhood_heatmap(
    adata: anndata.AnnData,
    title: str,
    cluster_key: str = "cluster",
    cmap: str = "bwr",
    vmin: float = -50,
    vmax: float = 50,
):
    key = f"{cluster_key}_nhood_enrichment"
    result = adata.uns.get(key)
    if not isinstance(result, dict) or "zscore" not in result:
        raise ValueError(f"Missing neighborhood enrichment results at adata.uns['{key}']")

    zscore = np.nan_to_num(
        np.asarray(result["zscore"], dtype=float),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    labels = _nhood_labels(adata, cluster_key, zscore.shape[0])
    zscore_df = pd.DataFrame(zscore, index=labels, columns=labels)

    fig_size = max(6, min(18, 0.35 * len(labels)))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    sns.heatmap(
        zscore_df,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=0,
        square=True,
    )
    ax.set_title(title)
    ax.set_xlabel(cluster_key)
    ax.set_ylabel(cluster_key)
    return fig


def _plot_nhood_enrichment(
    adata: anndata.AnnData,
    title: str,
    cluster_key: str = "cluster",
    method: str = "single",
    cmap: str = "bwr",
    vmin: float = -50,
    vmax: float = 50,
):
    from squidpy.pl import nhood_enrichment

    result = adata.uns.get(f"{cluster_key}_nhood_enrichment")
    if isinstance(result, dict) and result.get("skipped") == "fewer_than_two_clusters":
        logging.warning(
            "Skipping neighborhood enrichment plot because fewer than two "
            "clusters are present."
        )
        return None

    _sanitize_nhood_enrichment(adata, cluster_key=cluster_key)
    try:
        return nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            method=method,
            title=title,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    except ValueError as e:
        if "finite values" not in str(e):
            raise

        logging.warning(
            "Squidpy neighborhood enrichment plot failed because the clustered "
            "distance matrix contained non-finite values. Falling back to an "
            "unclustered heatmap."
        )
        return _fallback_nhood_heatmap(
            adata,
            title=title,
            cluster_key=cluster_key,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )


def get_custom_group_order(groups):
    """Sorts groups numerically first, then alphabetically, with 'Union'
    at the end."""
    groups = [str(g) for g in groups]  # Ensure everything is a string
    numeric_groups = sorted([g for g in groups if g.isdigit()], key=int)
    non_numeric_groups = sorted(
        [g for g in groups if not g.isdigit() and g != "Union"]
    )
    return numeric_groups + non_numeric_groups + ["Union"]


def _get_page_saver(output_path: str) -> tuple[Callable, Callable, bool, List[str]]:
    """Return a page saver and closer for PDF or image outputs."""
    ext = os.path.splitext(output_path)[1].lower()
    page_idx = 1
    pdf = PdfPages(output_path) if ext == ".pdf" else None
    output_stem = os.path.splitext(output_path)[0] if ext else output_path
    output_ext = ext if ext else ".png"
    saved_paths: List[str] = []

    def save_page(fig):
        nonlocal page_idx
        if pdf is not None:
            pdf.savefig(fig)
            return

        image_path = f"{output_stem}_{page_idx:03d}{output_ext}"
        fig.savefig(image_path, dpi=200, bbox_inches="tight")
        saved_paths.append(image_path)
        page_idx += 1

    def close():
        if pdf is not None:
            pdf.close()

    return save_page, close, pdf is not None, saved_paths


def _write_html_gallery(
    output_path: str,
    title: str,
    image_paths: List[str],
    captions: Optional[List[str]] = None,
    html_output_path: Optional[str] = None,
) -> None:
    """Write a single HTML file that displays all saved image pages."""
    if len(image_paths) == 0:
        return

    html_path = (
        html_output_path
        if html_output_path is not None
        else f"{os.path.splitext(output_path)[0]}.html"
    )
    html_dir = os.path.dirname(html_path) or "."
    relative_paths = [os.path.relpath(path, start=html_dir) for path in image_paths]

    blocks = []
    for idx, rel_path in enumerate(relative_paths, start=1):
        caption = f"Page {idx}"
        if captions is not None and idx - 1 < len(captions):
            caption = captions[idx - 1]

        blocks.append(
            "<section class=\"page\">"
            f"<h2>{caption}</h2>"
            f"<img loading=\"lazy\" src=\"{rel_path}\" alt=\"{caption}\" />"
            "</section>"
        )

    body = "\n".join(blocks)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <style>
    body {{
      margin: 24px auto;
      max-width: 1600px;
      padding: 0 16px 32px;
      background: #f4f4f4;
      color: #1f1f1f;
      font-family: Arial, sans-serif;
    }}
    h1 {{
      margin: 0 0 8px;
    }}
    p {{
      margin: 0 0 20px;
      color: #404040;
    }}
    .page {{
      margin: 0 0 24px;
      background: #fff;
      border: 1px solid #ddd;
      border-radius: 8px;
      padding: 12px;
    }}
    .page h2 {{
      margin: 0 0 10px;
      font-size: 18px;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
    }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <p>{len(image_paths)} page{"s" if len(image_paths) != 1 else ""}</p>
  {body}
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def plot_neighborhoods(
    adata: anndata.AnnData,
    group: str,
    subgroups: Optional[List[str]],
    outdir: str = "figures"
):

    if group != "all" and subgroups:
        filtered_adatas = {}
        for sg in subgroups:
            filtered_adata = filter_anndata(adata, group, sg)
            squidpy_analysis(filtered_adata)
            filtered_adatas[sg] = filtered_adata

    plt.rcParams.update({'figure.autolayout': True})
    with PdfPages(f"{outdir}/{group}_neighborhoods.pdf") as pdf:

        if subgroups:
            for sg in subgroups:
                fig = _plot_nhood_enrichment(
                    filtered_adatas[sg],
                    title=f"{group} {sg}: Neighborhood enrichment",
                )
                if fig is None:
                    continue
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        elif group == "all":
            fig = _plot_nhood_enrichment(
                adata,
                title="All cells: Neighborhood enrichment",
            )
            if fig is None:
                return

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def plot_spatial(
    adata: anndata.AnnData,
    samples: List[str],
    color_by: str,
    output_path: str,
    pt_size: int = 75,
    html_output_path: Optional[str] = None,
) -> None:
    """Plot cells spatially, color by metadata stored in .obs. The function
    creates a plot for each run and saves paginated output, with four runs per
    page.
    """
    from squidpy.pl import spatial_scatter

    page_captions: List[str] = []
    save_page, close, is_pdf, image_paths = _get_page_saver(output_path)
    try:
        for i in range(0, len(samples), 4):

            sample_batch = samples[i:i + 4]
            fig, axs = plt.subplots(2, 2, figsize=(10, 10))
            axs = axs.flatten()

            for i, sample in enumerate(sample_batch):
                # squidpy inverts y internally; restore positive row values before plotting.
                sample_adata = adata[adata.obs["sample"] == sample].copy()
                sample_adata.obsm["spatial"][:, 1] *= -1

                spatial_scatter(
                    sample_adata,
                    color=color_by,
                    size=pt_size,
                    shape=None,
                    library_id=sample,
                    ax=axs[i],
                    title=f"{sample}: {color_by}"
                )
                axs[i].axis("off")

            # Ensure empty plots are not displayed
            for j in range(len(sample_batch), 4):
                axs[j].axis("off")

            plt.tight_layout()

            page_captions.append("Samples: " + ", ".join(sample_batch))
            save_page(fig)
            plt.close(fig)
    finally:
        close()
    if not is_pdf:
        _write_html_gallery(
            output_path,
            title=f"Spatial Plots: {color_by}",
            image_paths=image_paths,
            captions=page_captions,
            html_output_path=html_output_path,
        )


def plot_spatial_qc(
    adata: anndata.AnnData,
    samples: List[str],
    qc_metrics: List[str],
    output_path: str,
    pt_size: int = 25,
    html_output_path: Optional[str] = None,
):
    """Generates a grid of spatial scatter plots for each sample and QC metric,
    saving them into paginated output. Each row corresponds to a sample and
    each column to a QC metric.
    """
    from squidpy.pl import spatial_scatter

    rows_per_page = 3
    cols_per_page = len(qc_metrics)

    page_captions: List[str] = []
    save_page, close, is_pdf, image_paths = _get_page_saver(output_path)
    try:
        for i in range(0, len(samples), rows_per_page):

            sample_batch = samples[i:i + rows_per_page]

            # Create a figure for the current page
            fig, axs = plt.subplots(
                len(sample_batch),
                cols_per_page,
                figsize=(cols_per_page * 5, len(sample_batch) * 5)
            )

            # If  one sample, make axs a list
            if len(sample_batch) == 1:
                axs = [axs]

            for row_idx, sample in enumerate(sample_batch):
                # squidpy inverts y internally; restore positive row values before plotting.
                sample_adata = adata[adata.obs['sample'] == sample].copy()
                sample_adata.obsm["spatial"][:, 1] *= -1

                for col_idx, qc_metric in enumerate(qc_metrics):

                    ax = axs[row_idx][col_idx]
                    spatial_scatter(
                        sample_adata,
                        color=qc_metric,
                        size=pt_size,
                        shape=None,
                        ax=ax,
                        library_id=sample,
                        title=f"{sample} : {qc_metric}",
                        colorbar=False
                    )
                    cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.7)

            plt.tight_layout()
            page_captions.append("Samples: " + ", ".join(sample_batch))
            save_page(fig)
            plt.close(fig)
    finally:
        close()
    if not is_pdf:
        _write_html_gallery(
            output_path,
            title="Spatial QC Plots",
            image_paths=image_paths,
            captions=page_captions,
            html_output_path=html_output_path,
        )


def plot_stacked_peaks(
    data: pd.DataFrame,
    group_by: str,
    color_by: str,
    group: str,
    output_path: str,
):
    """
    Plots a stacked bar chart of peak counts, grouped and colored by specified
    columns.
    Parameters:
    - data (pd.DataFrame): Input DataFrame containing peak data.
    - group_by (str): Column to use for x-axis grouping (e.g., "group").
    - color_by (str): Column to color the bars by (e.g., "peakType").
    - group (str): Group name for the title.
    - output_path (str): File path to save the plot as a PDF.
    Returns:
    - None (saves the plot to the specified path).
    """

    # Ensure required columns exist
    if group_by not in data.columns or color_by not in data.columns:
        logging.error(
            f"Missing required columns: {group_by} or {color_by} not found in data."
        )
        return

    # Ensure 'group' is treated as a string
    data["group"] = data["group"].astype(str)

    # Sort the groups numerically first, then alphabetically, with "Union" at the end
    custom_order = get_custom_group_order(data["group"].unique())

    # Convert 'group' column to categorical with specified order
    data["group"] = pd.Categorical(
        data["group"], categories=custom_order, ordered=True
    )

    # Create plot
    plt.figure(figsize=(8, 5))
    sns.histplot(
        data=data,
        x=group_by,
        hue=color_by,
        multiple="stack",
        discrete=True,
        alpha=0.5,
        edgecolor="gray",
        linewidth=0.5,
        shrink=0.9  # Adjust bar width to add spacing
    )

    plt.xlabel("Group")
    plt.ylabel("Peak Count")
    plt.title(f"Peak Counts by {group} and PeakType")

    plt.xticks(rotation=45)

    # Save figure
    plt.savefig(output_path, format="pdf", bbox_inches="tight")

    plt.close()


def plot_umaps(
    adata: anndata.AnnData,
    groups: List[str],
    output_path: str,
    html_output_path: Optional[str] = None,
) -> None:
    """Create a figure with UMAPs colored categorical metadata.
    """

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    axs = axs.flatten()

    for i in range(len(groups)):
        group = groups[i]
        sc.pl.umap(
            adata,
            s=10,
            color=group,
            ax=axs[i],
            show=False,
            title=f"UMAP: colored by {group}"
        )

    # Ensure empty plots are not displayed
    for j in range(len(axs)):
        axs[j].axis("off")

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    if html_output_path is not None and output_path.lower().endswith(".png"):
        _write_html_gallery(
            output_path,
            title="UMAP Plots",
            image_paths=[output_path],
            captions=["UMAP overview"],
            html_output_path=html_output_path,
        )
