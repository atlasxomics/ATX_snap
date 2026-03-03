import anndata
import logging
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scanpy as sc

from pathlib import Path
from typing import List, Optional

from wf.spatial import squidpy_analysis
from wf.utils import filter_anndata


def get_custom_group_order(groups):
    """Sorts groups numerically first, then alphabetically, with 'Union'
    at the end."""
    groups = [str(g) for g in groups]  # Ensure everything is a string
    numeric_groups = sorted([g for g in groups if g.isdigit()], key=int)
    non_numeric_groups = sorted(
        [g for g in groups if not g.isdigit() and g != "Union"]
    )
    return numeric_groups + non_numeric_groups + ["Union"]


def _page_output_path(output_path: str, page_num: int, total_pages: int) -> str:
    path = Path(output_path)
    if total_pages <= 1:
        return str(path)
    return str(path.with_name(f"{path.stem}_page_{page_num:02d}{path.suffix}"))


def _save_figure(
    fig: plt.Figure, output_path: str, page_num: int = 1, total_pages: int = 1
) -> None:
    fig.savefig(
        _page_output_path(output_path, page_num, total_pages),
        dpi=200,
        bbox_inches="tight"
    )


def plot_neighborhoods(
    adata: anndata.AnnData,
    group: str,
    subgroups: Optional[List[str]],
    outdir: str = "figures"
):
    from squidpy.pl import nhood_enrichment

    if group != "all" and subgroups:
        filtered_adatas = {}
        for sg in subgroups:
            filtered_adata = filter_anndata(adata, group, sg)
            squidpy_analysis(filtered_adata)
            filtered_adatas[sg] = filtered_adata

    plt.rcParams.update({'figure.autolayout': True})

    if subgroups:
        for sg in subgroups:
            fig = nhood_enrichment(
                filtered_adatas[sg],
                cluster_key="cluster",
                method="single",
                title=f"{group} {sg}: Neighborhood enrichment",
                cmap="bwr",
                vmin=-50,
                vmax=50,
            )
            _save_figure(fig, f"{outdir}/{group}_{sg}_neighborhoods.png")
            plt.close(fig)

    elif group == "all":
        fig = nhood_enrichment(
            adata,
            cluster_key="cluster",
            method="single",
            title="All cells: Neighborhood enrichment",
            cmap="bwr",
            vmin=-50,
            vmax=50,
        )

        _save_figure(fig, f"{outdir}/{group}_neighborhoods.png")
        plt.close(fig)


def plot_spatial(
    adata: anndata.AnnData,
    samples: List[str],
    color_by: str,
    output_path: str,
    pt_size: int = 75
) -> None:
    """Plot cells spatially, color by metadata stored in .obs.
    Creates up to four runs per page and saves as PNG.
    """
    from squidpy.pl import spatial_scatter

    total_pages = max(1, (len(samples) + 3) // 4)
    page_num = 0
    for start in range(0, len(samples), 4):
        page_num += 1
        sample_batch = samples[start:start + 4]
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        axs = axs.flatten()

        for idx, sample in enumerate(sample_batch):

            spatial_scatter(
                adata[adata.obs["sample"] == sample],
                color=color_by,
                size=pt_size,
                shape=None,
                library_id=sample,
                ax=axs[idx],
                title=f"{sample}: {color_by}"
            )
            axs[idx].axis("off")

        # Ensure empty plots are not displayed
        for j in range(len(sample_batch), 4):
            axs[j].axis("off")

        plt.tight_layout()

        _save_figure(fig, output_path, page_num=page_num, total_pages=total_pages)
        plt.close(fig)


def plot_spatial_qc(
    adata: anndata.AnnData,
    samples: List[str],
    qc_metrics: List[str],
    output_path: str,
    pt_size: int = 25
):
    """Generate a grid of spatial scatter plots for each sample and QC metric.
    Each row corresponds to a sample and each column to a QC metric.
    """
    from squidpy.pl import spatial_scatter

    rows_per_page = 3
    cols_per_page = len(qc_metrics)
    total_pages = max(1, (len(samples) + rows_per_page - 1) // rows_per_page)

    page_num = 0
    for i in range(0, len(samples), rows_per_page):
        page_num += 1

        sample_batch = samples[i:i + rows_per_page]

        # Create a figure for the current page
        fig, axs = plt.subplots(
            len(sample_batch),
            cols_per_page,
            figsize=(cols_per_page * 5, len(sample_batch) * 5)
        )

        # If one sample, make axs a list
        if len(sample_batch) == 1:
            axs = [axs]

        for row_idx, sample in enumerate(sample_batch):
            for col_idx, qc_metric in enumerate(qc_metrics):

                ax = axs[row_idx][col_idx]
                spatial_scatter(
                    adata[adata.obs['sample'] == sample],
                    color=qc_metric,
                    size=pt_size,
                    shape=None,
                    ax=ax,
                    library_id=sample,
                    title=f"{sample} : {qc_metric}",
                    colorbar=False
                )
                fig.colorbar(ax.collections[0], ax=ax, shrink=0.7)

        plt.tight_layout()
        _save_figure(fig, output_path, page_num=page_num, total_pages=total_pages)
        plt.close(fig)


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
    - output_path (str): File path to save the plot.
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
    plt.savefig(output_path, dpi=200, bbox_inches="tight")

    plt.close()


def plot_umaps(
    adata: anndata.AnnData, groups: List[str], output_path: str
) -> None:
    """Create a figure with UMAPs colored categorical metadata.
    """

    _, axs = plt.subplots(2, 2, figsize=(10, 10))
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

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
