import anndata
import matplotlib.pyplot as plt
import plotly.io as pio
import scanpy as sc
import snapatac2 as snap
import squidpy as sq

from matplotlib.backends.backend_pdf import PdfPages
from plotly.subplots import make_subplots
from typing import List


def plot_umaps(
    adata: anndata.AnnData, groups: List[str], output_path: str
) -> None:
    """Create a figure with UMAPs colored by categorical metadata using
    snapatac2.pl.umap.
    """
    # Create a 2x2 grid for subplots
    rows = 2
    cols = 2
    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=[f"UMAP: colored by {group}" for group in groups]
    )

    # Iterate over groups and generate UMAP plots
    for idx, group in enumerate(groups):
        row = idx // cols + 1
        col = idx % cols + 1

        umap_plot = snap.pl.umap(
            adata,
            color=group,
            show=False,
            interactive=False
        )

        # Add the trace from the UMAP plot to the subplot
        for trace in umap_plot.data:
            fig.add_trace(trace, row=row, col=col)

    # Update the layout
    fig.update_layout(
        height=800,
        width=800,
        title="UMAPs Colored by Groups",
        showlegend=False,
        plot_bgcolor="rgba(0,0,0,0)"
    )

    # Save the figure
    pio.write_image(fig, output_path)


def plot_spatial(
    adata: anndata.AnnData,
    samples: List[str],
    color_by: str,
    output_path: str,
    pt_size: int = 75
) -> None:
    """
    Plot spatial projections, color by metadata stored in .obs.
    Creates a plot for each batch of samples and saves to a PDF using Plotly.
    """
    # Create a Plotly PDF writer
    pdf_pages = []

    for i in range(0, len(samples), 4):
        sample_batch = samples[i:i + 4]

        # Create a 2x2 grid layout for the samples in the batch
        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=[f"{sample}: {color_by}" for sample in sample_batch]
        )

        for j, sample in enumerate(sample_batch):
            # Subset the data for the current sample
            obs_indices = (adata.obs["sample"] == sample).to_numpy().nonzero()[0]
            adata_sub = adata.subset(obs_indices=obs_indices, out=sample)[0]

            # Generate the spatial plot using snapatac2
            umap_fig = snap.pl.umap(
                adata_sub,
                color=color_by,
                use_rep="spatial",
                marker_size=pt_size,
                show=False  # Ensure it returns a Plotly figure
            )

            # Add each trace from the UMAP Plotly figure to the subplot
            row, col = divmod(j, 2)
            for trace in umap_fig.data:
                fig.add_trace(trace, row=row + 1, col=col + 1)

        # Adjust layout for better aesthetics
        fig.update_layout(
            height=800,
            width=800,
            showlegend=False,
            title=f"Spatial Plots: {color_by}"
        )

        # Save the current figure to the list of PDF pages
        pdf_pages.append(fig)

    # Save all figures to a PDF file
    with open(output_path, "wb") as f:
        pio.write_image(pdf_pages, file=f)


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
