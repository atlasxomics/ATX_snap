from dataclasses import dataclass
from enum import Enum

import snapatac2 as snap

from latch.types import LatchFile, LatchDir


class Genome(Enum):
    mm10 = "mm10"
    hg38 = "hg38"
    rnor6 = "rnor6"


@dataclass
class Run:
    run_id: str
    fragments_file: LatchFile
    condition: str = "None"
    spatial_dir: LatchDir = LatchDir(
        "latch:///spatials/demo/spatial/"
    )
    positions_file: LatchFile = LatchFile(
        "latch:///spatials/demo/spatial/tissue_positions_list.csv"
    )


genome_dict = {
    "mm10": snap.genome.mm10,
    "hg38": snap.genome.hg38
}
