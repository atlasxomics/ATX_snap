from atx_common import Genome
from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter, LatchRule

from wf.peak_recovery import restore_peaks_task


metadata = LatchMetadata(
    display_name="ArchR peak recovery",
    author=LatchAuthor(name="AtlasXomics Inc.", email="jamesm@atlasxomics.com",
                       github="https://github.com/atlasxomics"),
    repository="https://github.com/atlasxomics/ATX_snap",
    license="MIT",
    parameters={
        "archr_project": LatchParameter(
            display_name="existing ArchRProject",
            description="Select the project folder containing Save-ArchR-Project.rds "
            "and ArrowFiles, not the parent results folder.",
        ),
        "genome": LatchParameter(display_name="genome", description="Genome used by the original workflow."),
        "project_name": LatchParameter(
            display_name="output project name",
            rules=[LatchRule(regex=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
                             message="Use letters, numbers, dots, underscores or hyphens; start with a letter or number.")],
        ),
        "output_dir": LatchParameter(display_name="output directory"),
        "group_by": LatchParameter(
            display_name="peak grouping override",
            description="Leave empty to select the last group from the original workflow: "
            "last condition column, Sample, or Clusters. Override if a later group failed in the original run.",
        ),
        "include_y_chromosome": LatchParameter(
            display_name="include Y chromosome", hidden=True,
            description="Use the same setting as the original analysis.",
        ),
    },
)


@workflow(metadata)
def snap_workflow(
    archr_project: LatchDir,
    genome: Genome,
    project_name: str,
    output_dir: LatchDir = LatchDir("latch:///archr_peak_recovery/"),
    group_by: str = "",
    include_y_chromosome: bool = False,
) -> LatchDir:
    """Restore ArchR peak calling results.

    Reuse an existing ArchRProject's Arrow files to rebuild its active peak set,
    PeakMatrix and motif annotations, saving the repaired project separately.
    """
    return restore_peaks_task(
        archr_project=archr_project,
        genome=genome,
        project_name=project_name,
        output_dir=output_dir,
        group_by=group_by,
        include_y_chromosome=include_y_chromosome,
    )
