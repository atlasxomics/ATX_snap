import logging

from dataclasses import dataclass
from typing import Dict

from latch.types import LatchDir, LatchFile

# Shared AtlasXomics utilities — re-exported for `import wf.utils as utils`
from atx_common import (  # noqa: F401
    Genome,
    copy_peak_files,
    create_output_directories,
    filter_anndata,
    get_blacklist_path,
    get_channels,
    get_genome_fasta,
    get_groups,
    get_LatchFile,
    hg38_chrsizes,
    mm10_chrsizes,
    move_files_to_directory,
    organize_outputs,
    pt_sizes,
    ref_dict,
    rnor6_chrsizes,
    sanitize_condition,
)


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


# ---------------------------------------------------------------------------
# Workflow-specific types and helpers
# ---------------------------------------------------------------------------
@dataclass
class Run:
    run_id: str
    sample_name: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"


def get_data_paths(outdir: LatchDir) -> Dict[str, str]:
    """Get paths to required data files."""
    base_path = outdir.remote_path
    return {
        "obs": LatchFile(f"{base_path}/tables/obs.csv").local_path,
        "spatial": LatchFile(f"{base_path}/tables/spatial.csv").local_path,
        "umap": LatchFile(f"{base_path}/tables/X_umap.csv").local_path,
        "spectral": LatchFile(f"{base_path}/tables/spectral.csv")
    }
