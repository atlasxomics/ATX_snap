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
    complete_results_task,
    gene_stats_task,
    gene_project_task,
    genes_task,
    make_adata,
    make_anndata_dataset_task,
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
        "include_y_chromosome": LatchParameter(
            display_name="include y chromosome",
            description="Include Y chromosome regions in SnapATAC2 and ArchR "
            "matrices/peak calling. Defaults to excluding Y chromosome regions.",
            batch_table_column=True,
            hidden=True,
        ),
        "output_dir": LatchParameter(
            display_name="output directory",
            description="Folder in Latch Data to save outputs; defaults to "
            "`atac_analysis_snap`. Outputs are saved in a subfolder named "
            "with the project name.",
            batch_table_column=True,
            hidden=True,
        ),
    },
)


resume_from_gene_export_metadata = LatchMetadata(
    display_name="atx_snap_resume_from_gene_export",
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
            description="The same run metadata used for the original workflow run.",
            batch_table_column=True,
            samplesheet=True,
        ),
        "results_dir": LatchParameter(
            display_name="checkpointed results directory",
            description=(
                "The results directory containing both the make_adata outputs "
                "and the ArchR project uploaded by gene_project_task."
            ),
            batch_table_column=True,
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="The reference genome used for the original workflow run.",
            batch_table_column=True,
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="The project name used for the original workflow run.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
        "include_y_chromosome": LatchParameter(
            display_name="include y chromosome",
            description="Use the same setting as the original workflow run.",
            batch_table_column=True,
        ),
    },
)


# @workflow(metadata)
# def snap_workflow(
#     runs: List[Run],
#     genome: Genome,
#     project_name: str,
#     tile_size: int = 5000,
#     n_features: int = 25000,
#     n_comps: int = 30,
#     resolution: float = 1.0,
#     clustering_iters: int = 1,
#     leiden_iters: int = -1,
#     min_cluster_size: int = 20,
#     min_tss: float = 2.0,
#     min_frags: int = 10,
#     include_y_chromosome: bool = False,
#     output_dir: LatchDir = LatchDir("latch:///atac_analysis_snap/"),
# ) -> LatchDir:
#     """
#     SnapATAC2 and ArchR analysis for spatial ATAC runs.
#     """

#     anndata_dataset = make_anndata_dataset_task(
#         runs=runs,
#         genome=genome,
#         project_name=project_name,
#         min_tss=min_tss,
#         min_frags=min_frags,
#         include_y_chromosome=include_y_chromosome,
#         tile_size=tile_size,
#         output_dir=output_dir,
#     )

#     results, _groups = make_adata(
#         runs=runs,
#         anndata_dataset=anndata_dataset,
#         genome=genome,
#         project_name=project_name,
#         resolution=resolution,
#         leiden_iters=leiden_iters,
#         n_comps=n_comps,
#         min_cluster_size=min_cluster_size,
#         min_tss=min_tss,
#         min_frags=min_frags,
#         include_y_chromosome=include_y_chromosome,
#         tile_size=tile_size,
#         n_features=n_features,
#         clustering_iters=clustering_iters,
#         output_dir=output_dir,
#     )

#     gene_project = gene_project_task(
#         runs=runs,
#         results_dir=results,
#         project_name=project_name,
#         genome=genome,
#         include_y_chromosome=include_y_chromosome,
#     )

#     gene_results = genes_task(
#         runs=runs,
#         results_dir=results,
#         gene_project_dir=gene_project,
#         project_name=project_name,
#         genome=genome,
#         include_y_chromosome=include_y_chromosome,
#     )

#     results_ge = combine_gene_h5ads_task(
#         runs=runs,
#         results_dir=results,
#         gene_results_dir=gene_results,
#         project_name=project_name,
#     )

#     results_motifs = motifs_task(
#         runs=runs,
#         results_dir=results,
#         gene_results_dir=gene_results,
#         project_name=project_name,
#         genome=genome,
#         include_y_chromosome=include_y_chromosome,
#     )

#     results_with_gene_stats = gene_stats_task(
#         runs=runs,
#         gene_results_dir=gene_results,
#         gene_expression_results_dir=results_ge,
#         results_root=results,
#         project_name=project_name,
#     )

#     final_results = complete_results_task(
#         base_results_dir=results,
#         gene_results_dir=gene_results,
#         gene_expression_results_dir=results_ge,
#         gene_stats_results_dir=results_with_gene_stats,
#         motif_results_dir=results_motifs,
#     )

#     uploaded_results = registry_task(runs=runs, results=final_results)

#     return uploaded_results


@workflow(resume_from_gene_export_metadata)
def snap_workflow(
    runs: List[Run],
    results_dir: LatchDir,
    genome: Genome,
    project_name: str,
    include_y_chromosome: bool = False,
) -> LatchDir:
    """Resume ATX Snap analysis from chunked gene export.

    Uses the checkpointed ArchR project already stored in the results directory
    and runs gene export and all remaining tasks without rebuilding the project.
    """

    gene_results = genes_task(
        runs=runs,
        results_dir=results_dir,
        gene_project_dir=results_dir,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
    )

    results_ge = combine_gene_h5ads_task(
        runs=runs,
        results_dir=results_dir,
        gene_results_dir=gene_results,
        project_name=project_name,
    )

    results_motifs = motifs_task(
        runs=runs,
        results_dir=results_dir,
        gene_results_dir=gene_results,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
    )

    results_with_gene_stats = gene_stats_task(
        runs=runs,
        gene_results_dir=gene_results,
        gene_expression_results_dir=results_ge,
        results_root=results_dir,
        project_name=project_name,
    )

    final_results = complete_results_task(
        base_results_dir=results_dir,
        gene_results_dir=gene_results,
        gene_expression_results_dir=results_ge,
        gene_stats_results_dir=results_with_gene_stats,
        motif_results_dir=results_motifs,
    )

    return registry_task(runs=runs, results=final_results)
