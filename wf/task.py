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
from latch.resources.tasks import custom_task
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


@custom_task(cpu=32, memory=128, storage_gib=4949)
def motif_task(
    input_dir: LatchDir,
    runs: List[utils.Run],
    genome: utils.Genome,
    project_name: str
) -> Tuple[LatchFile, LatchDir]:
    """Get Anndata object with motifs matrix from cluster peak matrix.  We
    seperated into a seperate task because of the high memory requirements.
    """

    logging.info("Downloading data from previous step...")
    anndata_path = f"{input_dir.local_path}/cluster_peaks.h5ad"
    cluster_peaks = anndata.read_h5ad(anndata_path)

    groups = utils.get_groups(runs)

    out_dir = "/root/motifs"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    logging.info("Downloading reference genome for motifs...")
    genome = genome.value  # Convert to str
    fasta = utils.get_genome_fasta(genome)

    logging.info("Preparing peak matrix for motifs...")
    cluster_peaks = ft.get_motifs(cluster_peaks, fasta.local_path)
    cluster_peaks.write("cluster_peaks.h5ad")

    # Have to convert X to float64 for pc.compute_deviations
    cluster_peaks.X = cluster_peaks.X.astype(np.float64)

    logging.info("Computing motif deviation matrix...")
    adata_motif = pc.compute_deviations(
        cluster_peaks, n_jobs=1, chunk_size=10000
    )

    # Copy over cell data
    adata_motif.obs = cluster_peaks.obs
    adata_motif.obsm = cluster_peaks.obsm

    ft.rank_features(
        adata_motif, groups=groups, feature_type="motifs", save=out_dir
    )

    # Plot heatmap for motifs
    sc.pl.rank_genes_groups_matrixplot(
        adata_motif,
        n_genes=5,
        groupby="cluster",
        values_to_plot="scores",
        key="cluster_motifs",
        min_logfoldchange=0.1,
        save="motifs"
    )

    # Move scanpy plots
    subprocess.run([f"mv /root/figures/* {figures_dir}"], shell=True)

    adata_motif.write(f"{out_dir}/combined_motifs.h5ad")

    logging.info("Uploading motif data to Latch ...")
    return (
        LatchFile(
            "cluster_peaks.h5ad",
            f"latch:///snap_outs/{project_name}/cluster_peaks.h5ad"
        ),
        LatchDir(
            out_dir,
            f"latch:///snap_outs/{project_name}/motifs"
        )
    )


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
