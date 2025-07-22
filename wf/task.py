import glob
import logging
import os
import subprocess
from typing import List

import snapatac2 as snap
from latch import message
from latch.registry.table import Table
from latch.resources.tasks import custom_task, small_task
from latch.types import LatchDir

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp
import wf.utils as utils


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
            adata,
            groupby=group,
            suffix=f"{group}.bw",
            bin_size=10,
            output_format="bigwig",
        )
        bws = glob.glob("*.bw")
        subprocess.run(["mv"] + bws + [coverage_dir])
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
def genes_task(
    runs: List[utils.Run],
    outdir: LatchDir,
    project_name: str,
    groups: List[str],
    genome: utils.Genome,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(outdir)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)

    logging.info("Running ArchR analysis...")
    # run subprocess R script to make .h5ad file
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_genes.R',
        project_name,
        genome,
        data_paths['obs'],
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

    # # Load and combine data
    # adata_gene = ft.load_and_combine_data("g_converted")

    # # Transfer auxiliary data to combined AnnData
    # ft.transfer_auxiliary_data(adata_gene, data_paths, groups)

    # # Run spatial analysis
    # adata_gene = sp.run_squidpy_analysis(adata_gene, dirs["figures"])

    # # Load differential analysis results
    # ft.load_analysis_results(adata_gene, "gene", groups)

    # # Organize outputs
    # utils.organize_outputs(project_name, dirs)

    # # Save AnnData
    # ft.save_anndata_objects(adata_gene, "_ge", dirs['base'])

    logging.info("Uploading data to Latch...")
    return LatchDir(str(dirs['base']), f"latch:///snap_outs/{project_name}")


@custom_task(cpu=50, memory=975, storage_gib=2000)
def motifs_task(
    runs: List[utils.Run],
    outdir: LatchDir,
    project_name: str,
    groups: List[str],
    genome: utils.Genome,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(outdir)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)

    # Download ArchRProject
    archrproj_path = LatchDir(
        f"{outdir.remote_path}/{project_name}_ArchRProject"
    ).local_path

    logging.info("Running ArchR analysis...")
    # run subprocess R script to make .h5ad file
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_motifs.R',
        project_name,
        genome,
        data_paths['obs'],
        archrproj_path,
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

    # Load and combine data
    adata_motif = ft.load_and_combine_data("m_converted")

    # Transfer auxiliary data to combined AnnData
    ft.transfer_auxiliary_data(adata_motif, data_paths, groups)

    # Load differential analysis results
    ft.load_analysis_results(adata_motif, "motif", groups)

    # Organize outputs
    utils.organize_outputs(project_name, dirs)

    # Save AnnData
    ft.save_anndata_objects(adata_motif, "_motifs", dirs['base'])

    logging.info("Uploading data to Latch...")
    return LatchDir(str(dirs['base']), f"latch:///snap_outs/{project_name}")


@small_task(cache=True)
def registry_task(runs: List[utils.Run], results: LatchDir) -> LatchDir:
    try:
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
    except Exception as e:
        logging.warning(f"Unexpected {e=}, {type(e)=}")
        return results
