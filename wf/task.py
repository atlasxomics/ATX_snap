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
from latch.resources.tasks import custom_task, small_gpu_task, large_gpu_task
from latch.types import LatchDir, LatchFile

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp
import wf.utils as utils

import rapids_singlecell as rsc
logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@large_gpu_task()
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

    samples = [run.run_id for run in runs]
    conditions = list({run.condition for run in runs})

    # Get channels for specifying plot point size, use max for now...
    channels = max({utils.get_channels(run) for run in runs})

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

    if min_frags == 0:
        logging.warning("Minimum fragments set to 0.")
        message(
            typ="warning",
            data={"title": "min_frags", "body": "Minimum fragments set to 0."}
        )

    # Preprocessing -----------------------------------------------------------
    logging.info("Creating AnnData objects...")
    adatas = pp.make_anndatas(runs, genome, min_frags=min_frags)
    adatas = pp.filter_adatas(adatas, min_tss=min_tss)

    logging.info("Adding tile matrix to objects...")
    snap.pp.add_tile_matrix(adatas, bin_size=tile_size)

    if len(samples) > 1:
        logging.info("Combining objects...")
        adata = pp.combine_anndata(adatas, samples, filename="combined")
    else:
        adata = adatas[0]

    logging.info(
        f"Selecting features with {n_features} features and \
        {clustering_iters} clustering iteration(s)"
    )
    snap.pp.select_features(
        adata, n_features=n_features, max_iter=clustering_iters
    )

    logging.info("Performing dimensionality reduction...")
    adata = pp.add_clusters(adata, resolution, leiden_iters, min_cluster_size)
    adata = sp.add_spatial(adata)  # Add spatial coordinates to tixels

    # Plotting --
    pl.plot_umaps(adata, groups, f"{figures_dir}/umap.pdf")
    pl.plot_spatial(
        adata,
        samples,
        "cluster",
        f"{figures_dir}/spatial_dim.pdf",
        pt_size=utils.pt_sizes[channels]["dim"]
    )
    pl.plot_spatial_qc(
        adata,
        samples,
        qc_metrics,
        f"{figures_dir}/spatial_qc.pdf",
        pt_size=utils.pt_sizes[channels]["qc"]
    )

    # Genes ------------------------------------------------------------------
    logging.info("Making gene matrix...")
    adata_gene = ft.make_geneadata(rsc, adata, genome)
    adata_gene.obs.to_csv("gene_metadata.csv")
    rsc.get.anndata_to_GPU(adata_gene)
    for group in groups:

        print(group)
        logging.info(f"Finding marker genes for {group}s...")
        rsc.tl.rank_genes_groups_logreg(
            adata_gene,
            groupby=group,
            method="t-test",
            use_raw=False
        )

        rsc.get.anndata_to_CPU(adata_gene)
        # Write marker genes to csv
        sc.get.rank_genes_groups_df(
            adata_gene,
            group=None,
            key=f"rank_genes_groups",
            pval_cutoff=0.05,
            log2fc_min=0.1
        ).to_csv(f"{out_dir}/marker_genes_per_{group}.csv", index=False)

    # # Plot heatmap for genes
    sc.pl.rank_genes_groups_matrixplot(
        adata_gene,
        n_genes=5,
        groupby="cluster",
        values_to_plot="scores",
        key="rank_genes_groups",
        min_logfoldchange=0.1,
        save="genes"
    )

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
            rsc, adata, genome, f"{group}_peaks", log_norm=True
        )

        peak_mats[group] = anndata_peak
        rsc.get.anndata_to_GPU(anndata_peak)
        
        logging.info("Finded marker peaks ...")
        rsc.tl.rank_genes_groups_logreg(
            peak_mats[group], groupby=group, method="wilcoxon"
        )

        rsc.get.anndata_to_CPU(anndata_peak)
        logging.info("Writing peak matrix ...")
        anndata_peak.write(f"{out_dir}/{group}_peaks.h5ad")  # Save AnnData

        logging.info("Writing marker peaks to .csv ...")
        sc.get.rank_genes_groups_df(  # Save as csv
            peak_mats[group], group=None, pval_cutoff=0.05, log2fc_min=0.1
        ).to_csv(f"{out_dir}/marker_peaks_per_{group}.csv", index=False)

    logging.info("Writing combined anndata with peaks ...")
    adata.write(f"{out_dir}/combined.h5ad")

    # Move scanpy plots
    subprocess.run([f"mv /root/figures/* {figures_dir}"], shell=True)

    logging.info("Uploading data to Latch ...")
    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")


@large_gpu_task()
def motif_task(
    input_dir: LatchDir, genome: utils.Genome, project_name: str
) -> Tuple[LatchFile, LatchFile]:
    """Get Anndata object with motifs matrix from cluster peak matrix.  We
    seperated into a seperate task because of the high memory requirements.
    """

    logging.info("Downloading data from previous step...")
    anndata_path = f"{input_dir.local_path}/cluster_peaks.h5ad"
    cluster_peaks = anndata.read_h5ad(anndata_path)

    logging.info("Downloading reference genome for motifs...")
    genome = genome.value  # Convert to str
    fasta = utils.get_genome_fasta(genome)

    logging.info("Preparing peak matrix for motifs...")
    cluster_peaks = ft.get_motifs(cluster_peaks, fasta.local_path)
    cluster_peaks.write("cluster_peaks.h5ad")

    # Have to convert X to float64 for pc.compute_deviations
    cluster_peaks.X = cluster_peaks.X.astype(np.float64)

    logging.info("Computing motif deviation matrix...")
    adata_motif = pc.compute_deviations(cluster_peaks, n_jobs=90)

    # Copy over cell data
    adata_motif.obs = cluster_peaks.obs
    adata_motif.obsm = cluster_peaks.obsm

    adata_motif.write("combined_motifs.h5ad")

    return (
        LatchFile(
            "cluster_peaks.h5ad",
            f"latch:///snap_outs/{project_name}/cluster_peaks.h5ad"
        ),
        LatchFile(
            "combined_motifs.h5ad",
            f"latch:///snap_outs/{project_name}/combined_motifs.h5ad"
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
