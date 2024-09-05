from typing import List

from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

from wf.task import snap_task
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
                        run_id and fragments.tsv file; optional: condition, \
                        tissue position file for filtering on/off tissue.  \
                        Note that multiple Conditions must be separted by '_' \
                        (i.e., Female-control-old).",
            batch_table_column=True,
            samplesheet=True
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "min_tss": LatchParameter(
            display_name="minimum TSS",
            description="The minimum numeric transcription start site (TSS) \
                        enrichment score required for a cell to pass \
                        filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="The minimum number of mapped ATAC-seq fragments \
                        required per cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the \
                        TileMatrix.",
            batch_table_column=True,
            hidden=True
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
    min_tss: float = 2.0,
    min_frags: int = 0,
    tile_size: int = 5000
) -> LatchDir:

    return snap_task(
        runs=runs,
        genome=genome,
        min_tss=min_tss,
        min_frags=min_frags,
        tile_size=tile_size,
        project_name=project_name
    )
