import logging
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import mkdtemp

from atx_common import Genome
from latch.resources.tasks import custom_task
from latch.types import LatchDir

logging.basicConfig(format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO)


@custom_task(cpu=8, memory=64, storage_gib=2000, retries=0)
def restore_peaks_task(
    archr_project: LatchDir,
    genome: Genome,
    project_name: str,
    output_dir: LatchDir,
    group_by: str,
    include_y_chromosome: bool,
) -> LatchDir:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", project_name):
        raise ValueError("Invalid output project name.")
    destination = f"{output_dir.remote_path.rstrip('/')}/{project_name}"
    source_remote = archr_project.remote_path.rstrip("/")
    # Never upload into the source project or one of its ancestors/children.
    if (destination == source_remote or destination.startswith(source_remote + "/")
            or source_remote.startswith(destination + "/")):
        raise ValueError("Choose an output location separate from the input project.")
    source = Path(archr_project.local_path)
    if not (source / "Save-ArchR-Project.rds").is_file() or not (source / "ArrowFiles").is_dir():
        raise ValueError("Input must be an ArchRProject folder containing Save-ArchR-Project.rds and ArrowFiles.")

    results = Path(mkdtemp(prefix="archr_peak_recovery_"))
    project = results / f"{project_name}_ArchRProject"
    # Work exclusively on a copy: ArchR writes PeakMatrix into the Arrow files.
    shutil.copytree(source, project)
    subprocess.run([
        "Rscript", "/root/wf/R/restore_peaks.R", str(project), genome.value,
        group_by, str(include_y_chromosome).lower(), "8",
    ], check=True)
    if not (project / "Save-ArchR-Project.rds").is_file():
        raise FileNotFoundError("The recovered ArchRProject was not saved.")
    logging.info("Uploading recovered project to %s", destination)
    return LatchDir(str(results), destination)
