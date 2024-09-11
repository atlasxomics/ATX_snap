import logging
import numpy as np
import os
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap

from typing import List

from latch.resources.tasks import custom_task
from latch.types import LatchDir

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp

from wf.utils import Genome, Run, get_genome_fasta


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@custom_task(cpu=62, memory=384, storage_gib=4949)
def snap_task(
    runs: List[Run],
    genome: Genome,
    resolution: float,
    iterations: int,
    min_cluster_size: int,
    min_tss: float,
    min_frags: int,
    tile_size: int,
    project_name: str
) -> LatchDir:

    samples = [run.run_id for run in runs]
    conditions = list({run.condition for run in runs})

    # Set 'groups' list for differential analysis
    groups = ["cluster"]
    if len(samples) > 1:
        groups.append("sample")
    if len(conditions) > 1:
        groups.append("condition")
    logging.info(f"Comparing features amoung groups {groups}.")

    qc_metrics = ["n_fragment", "log10_frags", "tsse"]

    genome = genome.value  # Convert to str

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    # Preprocessing -----------------------------------------------------------
    logging.info("Creating AnnData objects...")
    adatas = pp.make_anndatas(runs, genome, min_frags=min_frags)
    adatas = pp.filter_adatas(adatas, min_tss=min_tss)

    logging.info("Adding tile matrix to objects...")
    adatas = pp.add_tilematrix(adatas, tile_size=tile_size, n_features=25000)

    logging.info("Combining objects...")
    adata = pp.combine_anndata(adatas, samples, filename="combined")
    snap.pp.select_features(adata, n_features=25000)

    logging.info("Performing dimensionality reduction...")
    adata = pp.add_clusters(adata, resolution, iterations, min_cluster_size)
    adata = sp.add_spatial(adata)  # Add spatial coordinates to tixels

    # Plotting --
    pl.plot_umaps(adata, groups, f"{figures_dir}/umap.pdf")
    pl.plot_spatial(
        adata, samples, "cluster", f"{figures_dir}/spatial_dim.pdf"
    )
    pl.plot_spatial_qc(
        adata, samples, qc_metrics, f"{figures_dir}/spatial_qc.pdf"
    )

    # Genes ------------------------------------------------------------------
    logging.info("Making gene matrix...")
    adata_gene = ft.make_geneadata(adata, genome)

    for group in groups:

        logging.info(f"Finding marker genes for {group}s...")
        sc.tl.rank_genes_groups(
            adata_gene,
            groupby=group,
            method="t-test",
            key_added=f"{group}_genes"
        )

        # Write marker genes to csv
        sc.get.rank_genes_groups_df(
            adata_gene,
            group=None,
            key=f"{group}_genes",
            pval_cutoff=0.05,
            log2fc_min=0.1
        ).to_csv(f"{out_dir}/marker_genes_per_{group}.csv", index=False)

    adata_gene.write(f"{out_dir}/combined_ge.h5ad")

    # Peaks ------------------------------------------------------------------
    peak_mats = {}
    for group in groups:

        logging.info(f"Calling peaks for {group}s...")
        snap.tl.macs3(
            adata,
            groupby=group,
            shift=-75,
            extsize=150,
            qvalue=0.1,
            key_added=f"{group}_peaks"
        )

        logging.info("Making peak matrix AnnData...")
        anndata_peak = ft.make_peakmatrix(
            adata, genome, f"{group}_peaks", log_norm=True
        )

        peak_mats[group] = anndata_peak

        logging.info("Finded marker peaks ...")
        sc.tl.rank_genes_groups(
            peak_mats[group], groupby=group, method="wilcoxon"
        )

        anndata_peak.write(f"{out_dir}/{group}_peaks.h5ad")  # Save AnnData

        sc.get.rank_genes_groups_df(  # Save as csv
            peak_mats[group], group=None, pval_cutoff=0.05, log2fc_min=0.1
        ).to_csv(f"{out_dir}/marker_peaks_per_{group}.csv", index=False)

    adata.write(f"{out_dir}/combined.h5ad")

    # Motifs -----------------------------------------------------------------
    logging.info("Downloading reference genome for motifs...")
    fasta = get_genome_fasta(genome)

    logging.info("Preparing peak matrix for motifs...")
    cluster_peaks = peak_mats["cluster"]  # Only motifs for cluster peaks
    cluster_peaks = ft.get_motifs(cluster_peaks, fasta.local_path)
    cluster_peaks.write(f"{out_dir}/cluster_peaks.h5ad")  # Save with motifs

    # Have to convert X to float64 for pc.compute_deviations
    cluster_peaks.X = cluster_peaks.X.astype(np.float64)

    logging.info("Computing motif deviation matrix...")
    adata_motif = pc.compute_deviations(cluster_peaks, n_jobs=90)
    adata_motif.write(f"{out_dir}/combined_motifs.h5ad")

    # Fin --------------------------------------------------------------------

    return LatchDir(
        out_dir, f"latch:///snap_outs/{project_name}"
    )
