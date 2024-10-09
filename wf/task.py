import anndata
import logging
import numpy as np
import os
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap
import subprocess

from typing import List, Tuple

from latch import message
from latch.resources.tasks import custom_task, small_gpu_task
from latch.types import LatchDir, LatchFile

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp
import wf.utils as utils


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@small_gpu_task()
def snap_task(
    runs: List[utils.Run],
    genome: utils.Genome,
    resolution: float,
    leiden_iters: int,
    min_cluster_size: int,
    min_tss: float,
    min_frags: int,
    tile_size: int,
    n_features: int,
    clustering_iters: int,
    project_name: str
) -> LatchDir:

    print('we made it')
    import rapids_singlecell as rsc
    
    print(dir(rsc))


@custom_task(cpu=62, memory=975, storage_gib=4949)
def motif_task(
    input_dir: LatchDir, genome: utils.Genome, project_name: str
) -> Tuple[LatchFile, LatchFile]:
    """Get Anndata object with motifs matrix from cluster peak matrix.  We
    seperated into a seperate task because of the high memory requirements.
    """

    print('motif issh')



if __name__ == "__main__":

    logging.info("Plotting SnapATAC peak heatmap...")
    # Perform SnapATAC marker peaks and heatmap

    anndata_peak = anndata.read_h5ad("cluster_peaks.h5ad")
    group = "cluster"
    genome = "hg38"

    marker_peaks = snap.tl.marker_regions(
        anndata_peak, groupby=group, pvalue=0.05
    )
    snap.pl.regions(
        anndata_peak,
        groupby=group,
        peaks=marker_peaks,
        interactive=False,
        out_file="snap_peak_heatmap.pdf"
    )

    logging.info("Plotting SnapATAC motif heatmap...")
    motifs = snap.tl.motif_enrichment(
        motifs=snap.datasets.cis_bp(unique=True),
        regions=anndata_peak,
        genome_fasta=(
            snap.genome.mm10 if genome == "mm10" else snap.genome.hg38
        )
    )
    snap.pl.motif_enrichment(
        motifs,
        max_fdr=0.0001,
        height=1600,
        interactive=False,
        out_file="motaf_enrichment.pdf"
    )
