import anndata
import glob
import json
import logging
import os
import subprocess

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List

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

rnor6_chrsizes = {
    "chr1": 282763074, "chr2": 266435125, "chr4": 184226339, "chr3": 177699992,
    "chr5": 173707219, "chrX": 159970021, "chr6": 147991367, "chr7": 145729302,
    "chr8": 133307652, "chr9": 122095297, "chr14": 115493446, "chr13": 114033958,
    "chr10": 112626471, "chr15": 111246239, "chr17": 90843779, "chr16": 90668790,
    "chr11": 90463843, "chr18": 88201929, "chr19": 62275575, "chr20": 56205956,
    "chr12": 52716770, "chrY": 3310458, "chrM": 16313,
}

ref_dict = {
    "mm10": [
        mm10_chrsizes, "/root/references/gencode_vM25_GRCm38.gff3.gz",
    ],
    "hg38": [
        hg38_chrsizes, "/root/references/gencode_v41_GRCh38.gff3.gz",
    ],
    "rnor6": [
        rnor6_chrsizes, "/root/references/rn6_liftoff_mm10_RefSeq.gff3.gz"
    ]
}


class Genome(Enum):
    mm10 = "mm10"
    hg38 = "hg38"
    rnor6 = "rnor6"


@dataclass
class Run:
    run_id: str
    sample_name: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"


def copy_peak_files(project_name: str, dirs: Dict[str, Path]) -> None:

    table_dir = f"/root/{project_name}/{project_name}_ArchRProject/PeakCalls/"
    plot_dir = f"/root/{project_name}/{project_name}_ArchRProject/Plots/"
    if os.path.exists(table_dir):
        peak_csvs = glob.glob(f"{table_dir}/*.csv")
        if not peak_csvs:
            logging.warning(f"No peaks csv files found in {table_dir}")
        else:
            subprocess.run(["cp"] + peak_csvs + [str(dirs['tables'])])
    else:
        logging.warning(f"No {table_dir} found")

    if os.path.exists(plot_dir):
        peaks_pdfs = glob.glob(f"{plot_dir}/*.pdf")
        if not peaks_pdfs:
            logging.warning(f"No peaks csv files found in {plot_dir}")
        else:
            subprocess.run(["cp"] + peaks_pdfs + [str(dirs['figures'])])
    else:
        logging.warning(f"No {plot_dir} found")


def create_output_directories(project_name: str) -> Dict[str, Path]:
    """Create necessary output directories and return paths."""
    base_dir = Path(f"/root/{project_name}")
    figures_dir = base_dir / "figures"
    tables_dir = base_dir / "tables"

    for directory in [base_dir, figures_dir, tables_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    return {
        'base': base_dir,
        'figures': figures_dir,
        'tables': tables_dir
    }


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


def get_data_paths(outdir: LatchDir) -> Dict[str, str]:
    """Get paths to required data files."""
    base_path = outdir.remote_path
    return {
        "obs": LatchFile(f"{base_path}/tables/obs.csv").local_path,
        "spatial": LatchFile(f"{base_path}/tables/spatial.csv").local_path,
        "umap": LatchFile(f"{base_path}/tables/X_umap.csv").local_path,
        "spectral": LatchFile(f"{base_path}/tables/spectral.csv")
    }


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


def move_files_to_directory(patterns: List[str], target_dir: Path) -> None:
    """Move files matching patterns to target directory."""
    files_to_move = []
    for pattern in patterns:
        files_to_move.extend(glob.glob(pattern))

    if files_to_move:
        subprocess.run(['mv'] + files_to_move + [str(target_dir)])


def organize_outputs(project_name: str, dirs: Dict[str, Path]) -> None:
    """Move output files to appropriate directories."""
    logging.info("Moving outputs to output directory...")

    # Move main project files
    project_patterns = [f'{project_name}_*', '*.rds', '*.h5ad']
    move_files_to_directory(project_patterns, dirs['base'])

    # Move tables
    csv_files = glob.glob('*.csv')
    if csv_files:
        subprocess.run(['mv'] + csv_files + [str(dirs['tables'])])

    # Move figures (excluding Rplots.pdf)
    figures = [fig for fig in glob.glob('*.pdf') if fig != 'Rplots.pdf']
    if figures:
        subprocess.run(['mv'] + figures + [str(dirs['figures'])])
