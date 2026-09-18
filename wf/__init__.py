from atx_common import Genome
from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter, LatchRule

from wf.motif_recovery import restore_motifs_task


metadata = LatchMetadata(
    display_name="ArchR motif recovery",
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
    },
)


@workflow(metadata)
def snap_workflow(
    archr_project: LatchDir,
    genome: Genome,
    project_name: str,
    output_dir: LatchDir = LatchDir("latch:///archr_motif_recovery/"),
) -> LatchDir:
    """Restore ArchR motif deviations.

    Preserve existing genes and peaks, rebuild motif annotations and deviations,
    and save a validated copy of the repaired project.
    """
    return restore_motifs_task(
        archr_project=archr_project,
        genome=genome,
        project_name=project_name,
        output_dir=output_dir,
    )
