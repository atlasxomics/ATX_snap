from typing import List

from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

# from wf.task import registry_task, snap_task
from wf.task import (
    combine_gene_h5ads_task,
    gene_stats_task,
    genes_task,
    motifs_task,
    patch_gene_stats_task,
    registry_task,
)
from atx_common import Genome
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
            description="List of runs to be analyzed; each run must contain a \
                run_id and fragments.tsv file; optional: condition, \
                alternative sample name. Spaces in condition labels are \
                normalized to '_' for downstream ArchR grouping.",
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
            description="Project name used for the existing make_adata output.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
        "results_dir": LatchParameter(
            display_name="make_adata results directory",
            description="Directory produced by make_adata containing the \
                required tables inputs (`obs.csv`, `X_umap.csv`, \
                `spatial.csv`, `spectral.csv`). Gene and motif outputs will \
                be written back into this same directory.",
            batch_table_column=True,
        ),
        "gene_artifacts_dir": LatchParameter(
            display_name="gene artifacts directory",
            description="Directory containing saved gene-stage artifacts from \
                a failed run, including the ArchRProject, `*_SeuratObj.rds`, \
                and any existing `*_g_converted.h5ad` files. The gene task \
                resumes by converting missing h5ad files from Seurat objects.",
            batch_table_column=True,
        ),
    },
)


gene_stats_metadata = LatchMetadata(
    display_name="atx_snap_gene_stats_patch",
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
            description="Runs from the completed ATX_snap launch. Sample names \
                are used to recreate sample-name gene statistic outputs.",
            batch_table_column=True,
            samplesheet=True,
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Project name used for the completed ATX_snap output.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
        "results_dir": LatchParameter(
            display_name="completed ATX_snap results directory",
            description="Completed results directory containing \
                `<project>_ArchRProject` and `combined_sm_ge.h5ad`.",
            batch_table_column=True,
        ),
        "gene_stats_threads": LatchParameter(
            display_name="gene statistics threads",
            description="Number of ArchR threads to use while running the \
                post-hoc gene differential statistics. Values above 50 are \
                clamped to 50 by the task.",
            batch_table_column=True,
        ),
    },
)


# @workflow(metadata)
# def snap_workflow(
#     runs: List[Run],
#     genome: Genome,
#     project_name: str,
#     results_dir: LatchDir,
#     gene_artifacts_dir: LatchDir,
# ) -> None:
#     """
#     Temporary entrypoint that skips `make_adata` and resumes from an existing
#     results directory produced by that task. The gene stage is split so
#     per-run Seurat and h5ad artifacts are uploaded before the combined gene
#     h5ad is built.
#     """

#     gene_results = genes_task(
#         runs=runs,
#         results_dir=results_dir,
#         gene_artifacts_dir=gene_artifacts_dir,
#         project_name=project_name,
#         genome=genome,
#     )

#     results_ge = combine_gene_h5ads_task(
#         runs=runs,
#         results_dir=results_dir,
#         gene_results_dir=gene_results,
#         project_name=project_name,
#     )

#     outdir_motifs = motifs_task(
#         runs=runs,
#         results_dir=results_ge,
#         project_name=project_name,
#         genome=genome,
#     )

#     uploaded_results = registry_task(runs=runs, results=outdir_motifs)

#     return uploaded_results


@workflow(gene_stats_metadata)
def snap_workflow(
    runs: List[Run],
    project_name: str,
    results_dir: LatchDir,
    gene_stats_threads: int = 1,
) -> LatchDir:
    """
    Post-hoc workflow for running the skipped gene differential statistics and
    patching the real outputs into `combined_sm_ge.h5ad`.
    """

    gene_stats = gene_stats_task(
        runs=runs,
        results_dir=results_dir,
        project_name=project_name,
        gene_stats_threads=gene_stats_threads,
    )

    patched_results = patch_gene_stats_task(
        runs=runs,
        results_dir=results_dir,
        gene_stats_dir=gene_stats,
        project_name=project_name,
    )

    return patched_results


if __name__ == "__main__":

    from latch.types import LatchDir, LatchFile

    # genes_task(
    #     runs=[Run(
    #         run_id="demo",
    #         fragments_file=LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
    #         spatial_dir=LatchDir("latch://atx-illumina.mount/Images_spatial/D1033/spatial/"),
    #         condition="control",
    #     )],
    #     outdir="latch://13502.account/snap_outs/demo_001716",
    #     genome=Genome("hg38"),
    #     project_name="develop_archrGenes",
    # )
    motifs_task(
        runs=[Run(
            run_id="demo",
            fragments_file=LatchFile("latch://13502.account/atx-archive-latch-1762285824.0858006/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
            spatial_dir=LatchDir("latch://atx-illumina.mount/Images_spatial/D1033/spatial/"),
            condition="demo",
        )],
        results_dir=LatchDir("latch://13502.account/Processed_Data/demo_002307_snap"),
        project_name="demo_002307_snap",
        genome=Genome("hg38"),
    )
