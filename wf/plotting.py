import anndata
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from matplotlib.backends.backend_pdf import PdfPages
from typing import List


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
            title=f"UMAP: colored by {groups[0]}"
        )

    # Ensure empty plots are not displayed
    for j in range(len(axs)):
        axs[j].axis("off")

    plt.tight_layout()

    plt.savefig(output_path)


def plot_spatial(
    adata: anndata.AnnData,
    samples: List[str],
    color_by: str,
    output_path: str,
    pt_size: int = 75
) -> None:
    """Plot cells spatially, color by metadata stored in .obs. The function
    creates a plot for each run and saves to a .pdf, with four runs per page.
    """

    with PdfPages(output_path) as pdf:
        for i in range(0, len(samples), 4):

            sample_batch = samples[i:i + 4]
            fig, axs = plt.subplots(2, 2, figsize=(10, 10))
            axs = axs.flatten()

            for i, sample in enumerate(sample_batch):

                sq.pl.spatial_scatter(
                    adata[adata.obs["sample"] == sample],
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

        pdf.savefig(fig)
        plt.close(fig)


def plot_spatial_qc(
    adata: anndata.AnnData,
    samples: List[str],
    qc_metrics: List[str],
    output_path: str,
    pt_size: int = 25
):
    """Generates a grid of spatial scatter plots for each sample and QC metric,
    saving them into a PDF.  Each row corresponds to a sample and each column
    to a QC metric.
    """

    rows_per_page = 3
    cols_per_page = len(qc_metrics)

    with PdfPages(output_path) as pdf:
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
                for col_idx, qc_metric in enumerate(qc_metrics):

                    ax = axs[row_idx][col_idx]
                    sq.pl.spatial_scatter(
                        adata[adata.obs['sample'] == sample],
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
            pdf.savefig(fig)
            plt.close(fig)
