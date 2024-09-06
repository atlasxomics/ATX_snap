from dataclasses import dataclass
from enum import Enum

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


def get_genome_fasta(genome: str) -> LatchFile:
    """Download reference genome fasta files from latch-public"""

    fasta_paths = {
        "mm10": "s3://latch-public/test-data/13502/GRCm38_genome.fa",
        "hg38":  "s3://latch-public/test-data/13502/GRCh38_genome.fa",
        "rnor6": "s3://latch-public/test-data/13502/Rnor6_genome.fa"
    }

    return LatchFile(fasta_paths[genome])
