from typing import List

from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

from wf.task import motif_task, snap_task
from wf.utils import Run, Genome


metadata = LatchMetadata(
    display_name="atx_snap",
    author=LatchAuthor(
        name="James McGann",
        email="jamesm@atlasxomics.com",
        github="github.com/atlasxomics"
    ),
    repository="https://github.com/atlasxomics/ATX_snap",
    license="MIT",
    parameters={
        "runs": LatchParameter(
            display_name="runs",
            description="List of runs to be analyzed; each run must contain a \
                run_id and fragments.tsv file; optional: condition, tissue \
                position file for filtering on/off tissue. Note that multiple \
                Conditions must be separted by '_' (i.e., Female-control).",
            batch_table_column=True
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "resolution": LatchParameter(
            display_name="clustering resolution",
            description="Clustering resolution for Leiden algorithm; higher \
                values result in more clusters.",
            batch_table_column=True
        ),
        "leiden_iters": LatchParameter(
            display_name="leiden iterations",
            description="Number of iterations for the leiden algorithm to \
                perform when assigning labels to clusters. 'Positive values \
                above 2 define the total number of  iterations to perform, -1 \
                has the algorithm run until it  reaches its  optimal \
                clustering.' - SnapATAC2 docs.",
            batch_table_column=True,
            hidden=True
        ),
        "min_cluster_size": LatchParameter(
            display_name="minimum cells per cluster",
            description="Minimum number of cells in a cluster.",
            batch_table_column=True,
            hidden=True
        ),
        "min_tss": LatchParameter(
            display_name="minimum TSS",
            description="The minimum numeric transcription start site (TSS) \
                enrichment score required for a cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="The minimum number of mapped fragments required  \
                per cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the \
                tile matrix.",
            batch_table_column=True,
            hidden=True
        ),
        "clustering_iters": LatchParameter(
            display_name="clustering iterations",
            description="Iterations performed when selecting variable \
                features for tile matrix. 'If greater than 1, this function \
                will perform iterative clustering and feature selection based \
                on variable features found using previous clustering results. \
                This is similar to the procedure implemented in ArchR... \
                Default value is 1, which means no iterative clustering is \
                 performed.'- SnapATAC2 docs",
            batch_table_column=True
        ),
        "n_features": LatchParameter(
            display_name="number of features",
            description="Number of features to be selected as 'most \
                accessible' in tile matrix.",
            batch_table_column=True
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of output directory in snap_outs/",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
    }
)


@workflow(metadata)
def snap_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    resolution: float = 1.0,
    leiden_iters: int = -1,
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
    tile_size: int = 5000,
    n_features: int = 25000,
    clustering_iters: int = 1
) -> LatchDir:

    results = snap_task(
        runs=runs,
        genome=genome,
        project_name=project_name,
        resolution=resolution,
        leiden_iters=leiden_iters,
        min_cluster_size=min_cluster_size,
        min_tss=min_tss,
        min_frags=min_frags,
        tile_size=tile_size,
        n_features=n_features,
        clustering_iters=clustering_iters
    )

    cluster_peaks, motifs = motif_task(
        input_dir=results,
        genome=genome,
        project_name=project_name
    )

    return results


if __name__ == "__main__":
    import anndata

    motif_task(
        cluster_peaks=anndata.read_h5ad("cluster_peaks.h5ad"),
        genome=Genome.hg38,
        project_name="latch_dev"
    )
