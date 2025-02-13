import anndata
import json
import logging

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List

from latch.types import LatchFile, LatchDir


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


# Map DBiT channels to plot point sizes for various spatial plots
pt_sizes = {
    50: {"dim": 75, "qc": 25},
    96: {"dim": 10, "qc": 5},
    210: {"dim": 0.25, "qc": 0.25},
    220: {"dim": 0.25, "qc": 0.25}
}

mm10_chrsizes = {
    "chr1": 195471971, "chr2": 182113224, "chr3": 160039680, "chr4": 156508116,
    "chr5": 151834684, "chr6": 149736546, "chr7": 145441459, "chr8": 129401213,
    "chr9": 124595110, "chr10": 130694993, "chr11": 122082543,
    "chr12": 120129022, "chr13": 120421639, "chr14": 124902244,
    "chr15": 104043685, "chr16": 98207768, "chr17": 94987271,
    "chr18": 90702639, "chr19": 61431566, "chrX": 171031299,
    "chrY": 91744698, "chrM": 16299
}

hg38_chrsizes = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622,
    "chr12": 133275309, "chr13": 114364328, "chr14": 107043718,
    "chr15": 101991189, "chr16": 90338345, "chr17": 83257441,
    "chr18": 80373285, "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569
}

ref_dict = {
    "mm10": [
        mm10_chrsizes,
        "/root/references/gencode_vM25_GRCm38.gff3.gz",
        "/root/references/mm10_genes.csv",
        "/root/references/mm10_exons.csv",
        "/root/references/mm10_promoters.csv",
    ],
    "hg38": [
        hg38_chrsizes,
        "/root/references/gencode_v41_GRCh38.gff3.gz",
        "/root/references/hg38_genes.csv",
        "/root/references/hg38_exons.csv",
        "/root/references/hg38_promoters.csv",
    ]
}


class Genome(Enum):
    mm10 = "mm10"
    hg38 = "hg38"
    rnor6 = "rnor6"


@dataclass
class Run:
    run_id: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"


def filter_anndata(
    adata: anndata.AnnData, group: str, subgroup: List[str]
) -> anndata.AnnData:
    return adata[adata.obs[group] == subgroup]


def get_channels(run: Run):
    spatial_dir = run.spatial_dir.local_path
    metadata_json = f"{spatial_dir}/metadata.json"

    with open(metadata_json, "r") as f:
        metadata = json.load(f)
        channels = metadata["numChannels"]

    return channels


def get_genome_fasta(genome: str) -> LatchFile:
    """Download reference genome fasta files from latch-public"""

    fasta_paths = {
        "mm10": "s3://latch-public/test-data/13502/GRCm38_genome.fa",
        "hg38":  "s3://latch-public/test-data/13502/GRCh38_genome.fa",
        "rnor6": "s3://latch-public/test-data/13502/Rnor6_genome.fa"
    }

    return LatchFile(fasta_paths[genome])


def get_groups(runs: List[Run]):
    """Set 'groups' list for differential analysis"""

    samples = [run.run_id for run in runs]
    conditions = list({run.condition for run in runs})

    groups = ["cluster"]
    if len(samples) > 1:
        groups.append("sample")
    if len(conditions) > 1:
        groups.append("condition")

    return groups


def get_LatchFile(directory: LatchDir, file_name: str) -> LatchFile:

    try:
        files = [file for file in directory.iterdir()
                 if isinstance(file, LatchFile) and
                 Path(file.path).name == file_name]

        if len(files) == 1:
            return files[0]
        elif len(files) == 0:
            raise FileNotFoundError(
                f"No file {file_name} found in {directory.remote_path}"
            )
        elif len(files) > 1:
            raise FileNotFoundError(
                f"Multiple files {file_name} found in {directory.remote_path}"
            )

    except Exception as e:
        logging.error(f"Failed to find file '{file_name}'; error {e}")
        return None
