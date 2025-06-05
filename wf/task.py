import glob
import logging
import os
import pickle
import subprocess
from typing import List

import anndata
import gc
import numpy as np
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap
from latch import message
from latch.registry.table import Table
from latch.resources.tasks import custom_task, large_gpu_task, small_task
from latch.types import LatchDir, LatchFile

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp
import wf.utils as utils
from wf.peaks import call_peaks_macs3_gpu


@custom_task(cpu=62, memory=512, storage_gib=1000)
def make_adata(
    runs: List[utils.Run],
    genome: utils.Genome,
    project_name: str,
    resolution: float,
    leiden_iters: int,
    n_comps: int,
    min_cluster_size: int,
    min_tss: float,
    min_frags: int,
    tile_size: int,
    n_features: int,
    clustering_iters: int,
) -> tuple[LatchDir, List[str]]:
    import pandas as pd

    samples = [run.run_id for run in runs]

    # Get channels for specifying plot point size, use max for now...
    channels = max({utils.get_channels(run) for run in runs})

    groups = utils.get_groups(runs)
    logging.info(f"Comparing features amoung groups {groups}.")

    qc_metrics = ["n_fragment", "log10_frags", "tsse"]

    genome = genome.value  # Convert to str

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{out_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)

    # Save input parameters to csv
    parameters = [
        ["project name", project_name],
        ["genome", genome],
        ["tile size", tile_size],
        ["number features", n_features],
        ["leiden iterations", leiden_iters],
        ["number of components", n_comps],
        ["minimum cluster size", min_cluster_size],
        ["minimum TSS", min_tss],
        ["minimum fragments", min_frags],
        ["clustering_resolution", resolution],
        ["clustering iterations", clustering_iters],
    ]

    pd.DataFrame(parameters, columns=["parameter", "value"]).to_csv(
        f"{tables_dir}/metadata.csv", index=False
    )

    if min_frags == 0:
        logging.warning("Minimum fragments set to 0.")
        message(
            typ="warning",
            data={"title": "min_frags", "body": "Minimum fragments set to 0."},
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
    snap.pp.select_features(adata, n_features=n_features, max_iter=clustering_iters)

    logging.info("Performing dimensionality reduction...")
    adata = pp.add_clusters(adata, resolution, n_comps, leiden_iters, min_cluster_size)

    adata = sp.add_spatial(adata)  # Add spatial coordinates to tixels

    logging.info("Creating coverages for groups...")
    for group in groups:
        coverage_dir = f"{out_dir}/{group}_coverages"
        os.makedirs(coverage_dir, exist_ok=True)
        snap.ex.export_coverage(
            adata, groupby=group, suffix=f"{group}.bedgraph.gz"
        )
        bgs = glob.glob("*.bedgraph.gz")
        subprocess.run(["mv"] + bgs + [coverage_dir])
    logging.info("Finished coverages for groups...")

    pl.plot_umaps(adata, groups, f"{figures_dir}/umap.pdf")
    pl.plot_spatial(
        adata,
        samples,
        "cluster",
        f"{figures_dir}/spatial_dim.pdf",
        pt_size=utils.pt_sizes[channels]["dim"],
    )
    pl.plot_spatial_qc(
        adata,
        samples,
        qc_metrics,
        f"{figures_dir}/spatial_qc.pdf",
        pt_size=utils.pt_sizes[channels]["qc"],
    )

    subprocess.run([f"mv /root/figures/* {figures_dir}"], shell=True)

    # Save critical data for gene matrix
    adata.obs.to_csv(f"{tables_dir}/obs.csv", index=True)

    umap_df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
    umap_df.to_csv(f"{tables_dir}/X_umap.csv")

    spatial_df = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names)
    spatial_df.to_csv(f"{tables_dir}/spatial.csv")

    adata.write(f"{out_dir}/combined.h5ad")

    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}"), groups


@custom_task(cpu=50, memory=975, storage_gib=2000)
def make_adata_gene(
    runs: List[utils.Run],
    outdir: LatchDir,
    project_name: str,
    groups: List[str],
    genome: utils.Genome,
) -> LatchDir:
    import pandas as pd
    from squidpy.pl import ripley

    # Read in data tables
    obs_path = LatchFile(f"{outdir.remote_path}/tables/obs.csv").local_path
    spatial_path = LatchFile(f"{outdir.remote_path}/tables/spatial.csv").local_path
    umap_path = LatchFile(f"{outdir.remote_path}/tables/X_umap.csv").local_path

    # adata = anndata.read_h5ad(data_path.local_path)
    genome = genome.value

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{out_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)

    logging.info("Making gene matrix...")

    # run subprocess R script to make .h5ad file
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_objs.R',
        project_name,
        genome,
        obs_path,
    ]

    runs = [
        (
            f'{run.run_id},'
            f'{run.fragments_file.local_path},'
            f'{run.condition},'
            f'{utils.get_LatchFile(run.spatial_dir, "tissue_positions_list.csv").local_path},'
            f'{run.spatial_dir.local_path},'
        )
        for run in runs
    ]
    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd, check=True)
    logging.info("Reading and combining gene AnnData...")
    h5ads_g = glob.glob("*g_converted.h5ad")
    adatas_g = [anndata.read_h5ad(h5ad) for h5ad in h5ads_g]
    adata_gene = sc.concat(adatas_g)

    logging.info("Reading and combining motif AnnData...")
    h5ads_m = glob.glob("*m_converted.h5ad")
    adatas_m = [anndata.read_h5ad(h5ad) for h5ad in h5ads_m]
    adata_motif = sc.concat(adatas_m)

    del adatas
    gc.collect()

    if "_index" in adata_gene.raw.var:  # Do this for some stupid reason
        adata_gene.raw.var.drop(columns=['_index'], inplace=True)
    if "_index" in adata_motif.raw.var:  # Do this for some stupid reason
        adata_motif.raw.var.drop(columns=['_index'], inplace=True)

    logging.info("Transferring auxiliary data...")
    try:
        obs = pd.read_csv(obs_path, index_col=0)
    except FileNotFoundError:
        logging.warning(f"File not found: {obs_path}")
        obs = None
    if obs is not None and not obs.empty:
        obs_aligned = obs.reindex(adata_gene.obs.index)
        adata_gene.obs = obs_aligned
        adata_motif.obs = obs_aligned
        for group in groups:
            if adata_gene.obs[group].dtype != object:  # Ensure groups are str
                adata_gene.obs[group] = adata_gene.obs[group].astype(str)
            if adata_motif.obs[group].dtype != object:  # Ensure groups are str
                adata_motif.obs[group] = adata_motif.obs[group].astype(str)

    logging.info("Transferring umap...")
    # Convert to DataFrame and include cell names
    umap_df = pd.read_csv(umap_path, index_col=0)
    umap_aligned = umap_df.loc[adata_gene.obs_names].values
    adata_gene.obsm["X_umap"] = umap_aligned
    adata_motif.obsm["X_umap"] = umap_aligned

    logging.info("Transferring spatial...")
    spatial_df = pd.read_csv(spatial_path, index_col=0)
    spatial_aligned = spatial_df.loc[adata_gene.obs_names].values
    adata_gene.obsm["spatial"] = spatial_aligned
    adata_motif.obsm["spatial"] = spatial_aligned

    logging.info("Running squidpy...")
    # Neighbrohood enrichment plot, Ripley's plot
    adata_gene = sp.squidpy_analysis(adata_gene)

    group_dict = dict()
    group_dict["all"] = None

    logging.info("Making neighborhood plots...")
    for group in group_dict.keys():
        pl.plot_neighborhoods(adata_gene, group, group_dict[group], outdir=figures_dir)

    logging.info("Running ripley...")
    ripley(adata_gene, cluster_key="cluster", mode="L", save="ripleys_L.pdf")

    # Add add DA gene tables
    logging.info("Adding ranked genes...")
    try:
        marker_files = glob.glob("ranked_genes_*.csv")
        for file in marker_files:
            name = file.split("/")[-1]
            name = name.replace(".csv", "")
            df = pd.read_csv(file, dtype={"group_name": str})
            adata_gene.uns[name] = df
    except Exception as e:
        logging.warning(f"Error {e} loading marker genes files.")

    # Add add DA motif tables
    try:
        marker_files = glob.glob("enrichedMotifs_*.csv")
        for file in marker_files:
            name = file.split("/")[-1]
            name = name.replace(".csv", "")
            df = pd.read_csv(file, dtype={"group_name": str})
            adata_motif.uns[name] = df
    except Exception as e:
        logging.warning(f"Error {e} loading marker genes files.")

    # Add archr genes heatmap
    logging.info("Adding heatmaps...")
    try:
        hm_files = glob.glob("genes_per_*_hm.csv")
        for file in hm_files:
            name = file.split("/")[-1]
            name = name.replace(".csv", "")
            df = pd.read_csv(file, dtype={"cluster": str}, index_col=0)
            adata_gene.uns[name] = df
    except Exception as e:
        logging.warning(f"Error {e} loading heatmap files.")

    # Add archr motifs heatmap
    try:
        hm_files = glob.glob("motif_per_*_hm.csv")
        for file in hm_files:
            name = file.split("/")[-1]
            name = name.replace(".csv", "")
            df = pd.read_csv(file, index_col=0)
            adata_motif.uns[name] = df
    except Exception as e:
        logging.warning(f"Error {e} loading heatmap files.")

    # Add archr gene volcano plots
    logging.info("Adding gene volcanos...")
    if "condition" in groups:
        try:
            volcano_files = glob.glob("volcanoMarkers_genes_*.csv")
            for file in volcano_files:
                name = file.split("/")[-1]
                treatment = name.replace("volcanoMarkers_genes_", "").replace(".csv", "")
                df = pd.read_csv(file, dtype={"cluster": str})
                adata_gene.uns[f"volcano_{treatment}"] = df
        except Exception as e:
            logging.warning(f"Error {e} loading volcano files.")

    logging.info("Adding motif volcanos...")
    # Add archr motif volcano plots
    if "condition" in groups:
        try:
            volcano_files = glob.glob("volcanoMarkers_motifs_*.csv") 
            for file in volcano_files:
                name = file.split("/")[-1]
                treatment = name.replace("volcanoMarkers_motifs_", "").replace(".csv", "")
                df = pd.read_csv(file, dtype={"cluster": str})
                adata_motif.uns[f"volcano_{treatment}"] = df
        except Exception as e:
            logging.warning(f"Error {e} loading volcano files.")

    logging.info("Moving outputs to ourdir...")
    # Move data files to subfolder
    project_dirs = glob.glob(f'{project_name}_*')
    seurat_objs = glob.glob('*.rds')
    h5_files = glob.glob('*.h5ad')
    _mv_cmd = (['mv'] + project_dirs + seurat_objs + h5_files + [out_dir])
    subprocess.run(_mv_cmd)

    # Move tables into subfolder
    csv_tables = glob.glob('*.csv')
    if len(csv_tables) > 0:
        _mv_tables_cmd = ['mv'] + csv_tables + [str(tables_dir)]
        subprocess.run(_mv_tables_cmd)

    # Move figures into subfolder
    figures = [fig for fig in glob.glob('*.pdf') if fig != 'Rplots.pdf']
    if len(figures) > 0:
        _mv_figures_cmd = ['mv'] + figures + [str(figures_dir)]
        subprocess.run(_mv_figures_cmd)

    logging.info("Saving full adata...")
    adata_gene.write(f"{out_dir}/combined_ge.h5ad")

    logging.info("Making reduced gene adata...")
    # Reduce size of anndata object for Plots
    sm_adata = ft.clean_adata(adata_gene)
    sm_adata_m = ft.clean_adata(adata_motif)

    logging.info("Saving gene adata...")
    adata_gene.write(f"{out_dir}/combined_ge.h5ad")
    sm_adata.write(f"{out_dir}/combined_sm_ge.h5ad")

    logging.info("Saving motif adata...")
    adata_motif.write(f"{out_dir}/combined_motifs.h5ad")
    sm_adata_m.write(f"{out_dir}/combined_sm_motifs.h5ad")

    logging.info("Uploading data to Latch...")
    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")


@large_gpu_task
def call_peaks(
    outdir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    groups: List[str],
) -> LatchDir:
    import pandas as pd

    genome = genome.value

    data_path = LatchFile(f"{outdir.remote_path}/combined.h5ad")
    adata = anndata.read_h5ad(data_path.local_path)

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{out_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)

    peak_mats = {}
    for group in groups:
        logging.info(f"Calling peaks for {group}s...")
        adata = call_peaks_macs3_gpu(
            adata=adata,
            groupby_key=group,
            d_treat=150,
            d_ctrl=10000,
            max_gap=30,
            peak_amp=150,
            q_thresh=0.1,
        )

        # Annotate all peaks for stacked bargraph
        logging.info("Annotating and plotting all peaks...")

        # Load reference features from CSV files
        feat_paths = utils.ref_dict[genome][2:5]
        feats = [pd.read_csv(feat) for feat in feat_paths]

        all_peaks = ft.make_plotting_peaks(
            adata, f"{group}_peaks", genome, feats
        )
        pl.plot_stacked_peaks(
            all_peaks,
            "group",
            "peakType",
            group,
            f"{figures_dir}/{group}Peaks_stackedBar.pdf"
        )

        logging.info("Making peak matrix AnnData...")
        anndata_peak = ft.make_peakmatrix(
            adata, genome, f"{group}_peaks", log_norm=False
        )

        logging.info("Writing peak matrix...")
        anndata_peak.write(f"{out_dir}/{group}_peaks.h5ad")  # Save AnnData

        peak_mats[group] = anndata_peak

        logging.info("Writing combined anndata with peaks info...")
        peaks = list(peak_mats["cluster"].var_names)
        snap.metrics.frip(adata, {"cluster_peaks": peaks})

        adata.write(f"{out_dir}/combined.h5ad")

        logging.info("Writing all peaks to disk...")
        peaks_path = f"{out_dir}/peaks.pkl"
        with open(peaks_path, "wb") as f:
            pickle.dump(peak_mats, f)

    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")


@custom_task(cpu=124, memory=192, storage_gib=1000)
def rank_peaks(
    outdir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    groups: List[str],
) -> LatchDir:
    import pandas as pd

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{out_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)

    genome = genome.value
    feats = [pd.read_csv(feat) for feat in utils.ref_dict[genome][2:5]]

    for group in groups:
        data_path = LatchFile(f"{outdir.remote_path}/{group}_peaks.h5ad")
        anndata_peak = anndata.read_h5ad(data_path.local_path)

        logging.info(f"Finding marker peaks {group}...")

        peaks_df = ft.rank_differential_peaks(anndata_peak, group)
        anndata_peak.uns["rank_peaks_groups"] = peaks_df

        logging.info("Writing marker peaks to .csv ...")
        peaks_df = peaks_df[peaks_df["adjusted p-value"] <= 0.05]
        peaks_df = ft.reformat_peak_df(peaks_df, "names", group_col="group")
        peaks_df = ft.annotate_peaks(peaks_df, feats)
        peaks_df.to_csv(f"{tables_dir}/marker_peaks_per_{group}.csv", index=False)

        logging.info("Writing ranked peak matrix...")
        anndata_peak.write(f"{out_dir}/{group}_peaks.h5ad")  # Save AnnData

    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")


@custom_task(cpu=62, memory=512, storage_gib=1000)
def motifs_task(
    outdir: LatchDir,
    project_name: str,
    groups: List[str],
    genome: utils.Genome,
) -> LatchDir:

    genome = genome.value
    data_path = LatchFile(f"{outdir.remote_path}/combined.h5ad")
    adata = anndata.read_h5ad(data_path.local_path)
    peaks_path = LatchFile(f"{outdir.remote_path}/peaks.pkl")
    with open(peaks_path.local_path, "rb") as f:
        peak_mats = pickle.load(f)
    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    figures_dir = f"{out_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{out_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)

    # Calculate the medians for each sample, create a DataFrame
    grouped = adata.obs.groupby("sample")

    medians_df = grouped.agg(
        {"n_fragment": "median", "tsse": "median", "cluster_peaks": "median"}
    ).reset_index()
    # Rename columns to match the desired output
    medians_df.rename(
        columns={"sample": "run_id", "tsse": "tss", "cluster_peaks": "frip"},
        inplace=True,
    )
    medians_df.to_csv(f"{tables_dir}/medians.csv", index=False)

    # Motifs ------------------------------------------------------------------
    # Get Anndata object with motifs matrix from cluster peak matrix.
    cluster_peaks = peak_mats["cluster"]

    logging.info("Preparing peak matrix for motifs...")
    fasta = utils.get_genome_fasta(genome)

    cluster_peaks = ft.get_motifs(cluster_peaks, fasta.local_path)
    cluster_peaks.write(f"{out_dir}/cluster_peaks.h5ad")

    # Have to convert X to float64 for pc.compute_deviations
    cluster_peaks.X = cluster_peaks.X.astype(np.float64)

    logging.info("Computing motif deviation matrix...")
    adata_motif = pc.compute_deviations(cluster_peaks, n_jobs=1, chunk_size=10000)

    # Copy over cell data
    adata_motif.obs = cluster_peaks.obs
    adata_motif.obsm = cluster_peaks.obsm

    ft.rank_features(adata_motif, groups=groups, feature_type="motifs", save=tables_dir)

    # Plot heatmap for motifs
    sc.pl.rank_genes_groups_matrixplot(
        adata_motif,
        n_genes=5,
        groupby="cluster",
        values_to_plot="scores",
        key="cluster_motifs",
        min_logfoldchange=0.1,
        save="motifs",
    )

    sm_adata = ft.clean_adata(adata_motif)

    adata_motif.write(f"{out_dir}/combined_motifs.h5ad")
    sm_adata.write(f"{out_dir}/combined_sm_motifs.h5ad")

    # Upload data -----------------------------------------------------------

    logging.info("Uploading data to Latch ...")

    # Move scanpy plots
    subprocess.run([f"mv /root/figures/* {figures_dir}"], shell=True)
    return LatchDir(out_dir, f"latch:///snap_outs/{project_name}")


@small_task(cache=True)
def registry_task(runs: List[utils.Run], results: LatchDir) -> LatchDir:
    tbl = Table(id="761")

    logging.info("Uploading results to Runs Table in Registry...")

    for run in runs:
        logging.info(f"Adding {run.run_id} results to Registry...")

        with tbl.update() as updater:
            updater.upsert_record(
                name=run.run_id,
                fragments_file=run.fragments_file,
                spatial_directory=run.spatial_dir,
                condition=run.condition,
                atx_snap_outs=results,
            )

    logging.info("Done uploading to Registry.")
    return results
