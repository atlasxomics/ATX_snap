from typing import List

from atx_common import Genome
from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    LatchRule,
)

from wf.task import (
    combine_gene_h5ads_task,
    gene_stats_task,
    genes_task,
    make_adata,
    motifs_task,
    registry_task,
)
from wf.utils import Run


metadata = LatchMetadata(
    display_name="atx_snap",
    author=LatchAuthor(
        name="James McGann",
        email="jamesm@atlasxomics.com",
        github="github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/ATX_snap",
    license="MIT",
    parameters={
        "runs": LatchParameter(
            display_name="runs",
            description="List of runs to be analyzed; each run must contain a "
            "run_id and fragments.tsv file; optional: condition, alternative "
            "sample name. Spaces in condition labels are normalized to '_' for "
            "downstream ArchR grouping.",
            batch_table_column=True,
            samplesheet=True,
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of subfolder in output directory.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the "
            "tile matrix.",
            batch_table_column=True,
        ),
        "n_features": LatchParameter(
            display_name="number of features",
            description="Number of features to be selected as most accessible "
            "in the tile matrix.",
            batch_table_column=True,
        ),
        "n_comps": LatchParameter(
            display_name="number of components",
            description="Number of components/dimensions to keep during "
            "dimensionality reduction with `snap.tl.spectral`.",
            batch_table_column=True,
        ),
        "resolution": LatchParameter(
            display_name="clustering resolution",
            description="Clustering resolution for Leiden algorithm; higher "
            "values result in more clusters.",
            batch_table_column=True,
        ),
        "clustering_iters": LatchParameter(
            display_name="clustering iterations",
            description="Iterations performed when selecting variable features "
            "for tile matrix.",
            batch_table_column=True,
        ),
        "gene_stats_threads": LatchParameter(
            display_name="gene statistics threads",
            description="Number of ArchR threads to use while running gene "
            "differential statistics. Values above 50 are clamped to 50.",
            batch_table_column=True,
        ),
        "leiden_iters": LatchParameter(
            display_name="leiden iterations",
            description="Number of iterations for the leiden algorithm.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_cluster_size": LatchParameter(
            display_name="minimum cells per cluster",
            description="Minimum number of cells in a cluster.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_tss": LatchParameter(
            display_name="minimum TSS",
            description="Minimum transcription start site enrichment score "
            "required for a cell to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="Minimum number of mapped fragments required per cell "
            "to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "output_dir": LatchParameter(
            display_name="output directory",
            description="Folder in Latch Data to save outputs; defaults to "
            "`snap_outs`. Outputs are saved in a subfolder named with the "
            "project name.",
            batch_table_column=True,
            hidden=True,
        ),
    },
)


@workflow(metadata)
def snap_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    tile_size: int = 5000,
    n_features: int = 25000,
    n_comps: int = 30,
    resolution: float = 1.0,
    clustering_iters: int = 1,
    gene_stats_threads: int = 25,
    leiden_iters: int = -1,
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
    output_dir: LatchDir = LatchDir("latch:///snap_outs/"),
) -> LatchDir:
    """
    SnapATAC2 and ArchR analysis for spatial ATAC runs.
    """

    results, _groups = make_adata(
        runs=runs,
        genome=genome,
        project_name=project_name,
        resolution=resolution,
        leiden_iters=leiden_iters,
        n_comps=n_comps,
        min_cluster_size=min_cluster_size,
        min_tss=min_tss,
        min_frags=min_frags,
        tile_size=tile_size,
        n_features=n_features,
        clustering_iters=clustering_iters,
        output_dir=output_dir,
    )

    gene_results = genes_task(
        runs=runs,
        results_dir=results,
        project_name=project_name,
        genome=genome,
    )

    results_ge = combine_gene_h5ads_task(
        runs=runs,
        results_dir=results,
        gene_results_dir=gene_results,
        project_name=project_name,
    )

    results_motifs = motifs_task(
        runs=runs,
        results_dir=results_ge,
        project_name=project_name,
        genome=genome,
    )

    results_with_gene_stats = gene_stats_task(
        runs=runs,
        results_dir=results_motifs,
        project_name=project_name,
        gene_stats_threads=gene_stats_threads,
    )

    uploaded_results = registry_task(runs=runs, results=results_with_gene_stats)

    return uploaded_results
