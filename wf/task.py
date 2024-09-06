import logging
import os
import pandas as pd
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap

from typing import List
from pyjaspar import jaspardb

from latch.resources.tasks import large_task
from latch.types import LatchDir

import wf.preprocessing as pp
import wf.spatial as sp
import wf.features as ft

from wf.utils import Genome, Run


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@large_task
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
    groups = ["cluster", "sample", "condition"]
    genome = genome.value  # Convert to str

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    # Preprocessing ------
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

    # Gene matrix ------
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

    # Peak calling ------
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

        logging.info
        sc.tl.rank_genes_groups(
            peak_mats[group], groupby=group, method="wilcoxon"
        )

        anndata_peak.write(f"{out_dir}/group_peaks.h5ad")  # Save AnnData
        sc.get.rank_genes_groups_df(  # Save as csv
            peak_mats[group], group=None, pval_cutoff=0.05, log2fc_min=0.1
        ).to_csv(f"{out_dir}/marker_peaks_per_{group}.csv", index=False)

    # Fin ------
    adata.write(f"{out_dir}/combined.h5ad")

    return LatchDir(
        out_dir, f"latch:///snap_outs/{project_name}"
    )
