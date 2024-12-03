import anndata
import json
import numpy as np
import snapatac2 as snap

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional

from latch.types import LatchFile, LatchDir


# Map DBiT channels to plot point sizes for various spatial plots
pt_sizes = {
    50: {"dim": 75, "qc": 25},
    96: {"dim": 10, "qc": 5},
    210: {"dim": 5, "qc": 0.5},
    220: {"dim": 5, "qc": 0.5}
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
    positions_file: LatchFile
    condition: str = "None"


def copy_adata(
    adata: anndata.AnnData,
    groups: List[str],
    obs: Optional[List[str]] = ["n_fragment", "tsse", "log10_frags"],
    obsm: Optional[List[str]] = ["spatial", "X_umap"]
) -> anndata.AnnData:
    """From SnapATAC2 backend, make a lightweight AnnData copy for plotting.
    """
    if "sample" not in groups:
        groups.append("sample")
    new_adata = anndata.AnnData()

    for group in groups:
        new_adata.obs[group] = adata.obs[group]

    new_adata.obs_names = adata.obs_names

    for ob in obs:
        new_adata.obs[ob] = adata.obs[ob]

    for ob in obsm:
        new_adata.obsm[ob] = adata.obsm[ob]
        if type(new_adata.obsm[ob]) is not np.ndarray:
            new_adata.obsm[ob] = new_adata.obsm[ob].to_numpy()

    return new_adata


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


def refresh_adata(adata: anndata.AnnData, file_name: str) -> anndata.AnnData:
    """Running with snapATAC2 backend results in .h5ad files frequently being
    closed, necessitating that they be regularly reopened with r+ permissions
    in order to be modified.  Here, we ensure the object is closed and then
    reopen with r+.
    """
    adata = adata.close()
    adata = snap.read(f"{file_name}.h5ad", "r+")
    return adata
