from typing import List

from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

# from wf.task import registry_task, snap_task
from wf.task import (
    genes_task,
    motifs_task,
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
    },
)


@workflow(metadata)
def snap_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    results_dir: LatchDir,
) -> None:
    """
    Temporary entrypoint that skips `make_adata` and resumes from an existing
    results directory produced by that task.
    """

    results_ge = genes_task(
        runs=runs,
        results_dir=results_dir,
        project_name=project_name,
        genome=genome,
    )

    outdir_motifs = motifs_task(
        runs=runs,
        results_dir=results_ge,
        project_name=project_name,
        genome=genome,
    )

    uploaded_results = registry_task(runs=runs, results=outdir_motifs)

    return uploaded_results


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
