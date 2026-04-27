import json
import logging
import os
import shutil
import subprocess

from pathlib import Path
from typing import Dict, List

import snapatac2 as snap

from latch import message
from latch.registry.table import Table
from latch.resources.tasks import custom_task, small_task
from latch.types import LatchDir
from latch.types.plots import (
    PlotsArtifactBindings, PlotsArtifactTemplate, PlotsArtifact, Widget
)

import wf.features as ft
import wf.plotting as pl
import wf.preprocessing as pp
import wf.spatial as sp
import wf.utils as utils


def _dirs_for_base(base_dir: Path) -> Dict[str, Path]:
    figures_dir = base_dir / "figures"
    tables_dir = base_dir / "tables"

    for directory in [base_dir, figures_dir, tables_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    return {
        "base": base_dir,
        "figures": figures_dir,
        "tables": tables_dir,
    }


def _copy_required_input_tables(data_paths: Dict[str, str], tables_dir: Path) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)

    for local_path in data_paths.values():
        source = Path(local_path)
        if not source.exists():
            logging.warning(f"Input table does not exist and was not copied: {source}")
            continue

        destination = tables_dir / source.name
        try:
            if source.resolve() == destination.resolve():
                continue
        except FileNotFoundError:
            pass

        shutil.copy2(source, destination)


def _copy_directory_contents(source_dir: Path, destination_dir: Path) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)

    if not source_dir.exists():
        logging.warning(f"Directory does not exist and was not copied: {source_dir}")
        return

    for source in source_dir.iterdir():
        destination = destination_dir / source.name
        if source.is_dir():
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)


def _gene_stats_run_args(runs: List[utils.Run]) -> List[str]:
    return [
        f"{run.run_id},{run.sample_name}"
        for run in runs
    ]


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
    output_dir: LatchDir
) -> tuple[LatchDir, List[str]]:
    import pandas as pd

    output_dir = output_dir.remote_path + "/" + project_name

    samples = [run.run_id for run in runs]

    # Get channels for specifying plot point size, use max for now...
    channels = max({utils.get_channels(run) for run in runs})

    groups = utils.get_groups(runs)
    logging.info(f"Comparing features amoung groups {groups}.")

    qc_metrics = ["n_fragment", "log10_frags", "tsse"]

    genome = genome.value  # Convert to str
    blacklist = utils.get_blacklist_path(genome)
    if blacklist is None:
        logging.warning("Proceeding without blacklist filtering.")
    else:
        logging.info(f"Using blacklist: {blacklist}")

    result_dir = f"/root/{project_name}"
    os.makedirs(result_dir, exist_ok=True)

    figures_dir = f"{result_dir}/figures"
    os.makedirs(figures_dir, exist_ok=True)

    tables_dir = f"{result_dir}/tables"
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
        f"{tables_dir}/input_parameters.csv", index=False
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
    snap.pp.select_features(
        adata,
        n_features=n_features,
        max_iter=clustering_iters,
        blacklist=blacklist,
    )

    logging.info("Performing dimensionality reduction...")
    adata, spectral_key = pp.add_clusters(
        adata, resolution, n_comps, leiden_iters, min_cluster_size
    )

    adata = sp.add_spatial(adata)  # Add spatial coordinates to tixels

    logging.info("Creating coverages for groups...")
    coverage_groups = groups if "sample" in groups else groups + ["sample"]
    for group in coverage_groups:
        coverage_dir = Path(result_dir) / f"{group}_coverages"
        coverage_dir.mkdir(parents=True, exist_ok=True)

        snap.ex.export_coverage(
            adata,
            groupby=group,
            suffix=f"{group}.bw",
            bin_size=10,
            output_format="bigwig",
            out_dir=coverage_dir,
            blacklist=blacklist,
        )

    # Optionally duplicate sample coverages with user-provided sample names.
    sample_name_map = {}
    for run in runs:
        sample_name = str(run.sample_name).strip() if run.sample_name is not None else ""
        if sample_name and sample_name.lower() != "none":
            safe_sample_name = sample_name.replace("/", "_")
            if safe_sample_name != sample_name:
                logging.warning(
                    f"Replacing '/' with '_' in sample_name '{sample_name}'."
                )
            sample_name_map[run.run_id] = safe_sample_name

    if sample_name_map:
        sample_cov_dir = Path(result_dir) / "sample_coverages"
        for source_bw in sample_cov_dir.glob("*.bw"):
            matched_run = None
            for run_id in sample_name_map:
                if source_bw.name.startswith(f"{run_id}."):
                    matched_run = run_id
                    break

            if matched_run is None:
                continue

            target_name = (
                f"{sample_name_map[matched_run]}"
                f"{source_bw.name[len(matched_run):]}"
            )
            target_bw = sample_cov_dir / target_name
            if target_bw.exists():
                logging.warning(
                    f"Skipping sample_name coverage copy; {target_name} already exists."
                )
                continue

            shutil.copy2(source_bw, target_bw)

    logging.info("Finished coverages for groups...")

    pl.plot_umaps(
        adata,
        groups,
        f"{figures_dir}/umap.png",
        html_output_path=f"{result_dir}/umap.html",
    )
    pl.plot_spatial(
        adata,
        samples,
        "cluster",
        f"{figures_dir}/spatial_dim.png",
        pt_size=utils.pt_sizes[channels]["dim"],
        html_output_path=f"{result_dir}/spatial_dim.html",
    )
    pl.plot_spatial_qc(
        adata,
        samples,
        qc_metrics,
        f"{figures_dir}/spatial_qc.png",
        pt_size=utils.pt_sizes[channels]["qc"],
        html_output_path=f"{result_dir}/spatial_qc.html",
    )

    subprocess.run([f"mv /root/figures/* {figures_dir}"], shell=True)

    # Save critical data for gene matrix
    adata.obs.to_csv(f"{tables_dir}/obs.csv", index=True)

    umap_df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
    umap_df.to_csv(f"{tables_dir}/X_umap.csv")

    spatial_df = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names)
    spatial_df.to_csv(f"{tables_dir}/spatial.csv")

    spectral_df = pd.DataFrame(adata.obsm[spectral_key], index=adata.obs_names)
    spectral_df.to_csv(f"{tables_dir}/spectral.csv")

    adata.write(f"{result_dir}/combined.h5ad")

    return LatchDir(result_dir, output_dir), groups


@custom_task(cpu=30, memory=700, storage_gib=2000)
def genes_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(results_dir)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)
    _copy_required_input_tables(data_paths, dirs["tables"])

    logging.info("Running ArchR analysis...")
    # run subprocess R script to make .h5ad file
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_genes.R',
        project_name,
        genome,
        data_paths['obs'],
        data_paths['spectral'],
    ]

    position_files = {}
    missing_positions = []
    for run in runs:
        position_file = utils.get_LatchFile(
            run.spatial_dir, "tissue_positions_list.csv"
        )
        if position_file is None:
            missing_positions.append(f"{run.run_id} ({run.spatial_dir.remote_path})")
            continue
        position_files[run.run_id] = position_file

    if missing_positions:
        missing_str = ", ".join(missing_positions)
        raise FileNotFoundError(
            "Unable to resolve 'tissue_positions_list.csv' for one or more runs: "
            f"{missing_str}. Ensure each spatial directory contains exactly one "
            "'tissue_positions_list.csv' file."
        )

    runs = [
        (
            f'{run.run_id},'
            f'{run.fragments_file.local_path},'
            f'{utils.sanitize_condition(run.condition)},'
            f'{position_files[run.run_id].local_path},'
            f'{run.spatial_dir.local_path},'
            f'{run.sample_name}'
        )
        for run in runs
    ]
    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd, check=True)

    # Stage ArchRProject, Seurat objects, per-run h5ads, and R-side tables.
    utils.organize_outputs(project_name, dirs)

    logging.info("Uploading gene-stage artifacts to Latch...")
    return LatchDir(str(dirs['base']), results_dir.remote_path)


@custom_task(cpu=16, memory=512, storage_gib=2000)
def combine_gene_h5ads_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    gene_results_dir: LatchDir,
    project_name: str,
) -> LatchDir:

    # Read in data tables from the original make_adata output.
    data_paths = utils.get_data_paths(results_dir)
    groups = utils.get_groups(runs)

    # Continue in the uploaded gene artifact directory so a combine failure
    # does not discard the already-materialized Seurat and per-run h5ad files.
    dirs = _dirs_for_base(Path(gene_results_dir.local_path))
    _copy_required_input_tables(data_paths, dirs["tables"])

    # Load and combine data
    adata_gene = ft.load_and_combine_data(
        "g_converted",
        input_dir=dirs["base"],
        temp_dir=dirs["base"],
    )

    # Transfer auxiliary data to combined AnnData
    ft.transfer_auxiliary_data(adata_gene, data_paths, groups)

    # Run spatial analysis
    sample_key = "sample" if "sample" in groups else None
    adata_gene = sp.run_squidpy_analysis(
        adata_gene, dirs["figures"], sample_key
    )

    # Load differential analysis results
    ft.load_analysis_results(
        adata_gene,
        "gene",
        groups,
        input_dir=dirs["tables"],
    )

    # Save AnnData
    ft.save_anndata_objects(adata_gene, "_ge", dirs["base"])

    logging.info("Uploading data to Latch...")
    return LatchDir(str(dirs['base']), results_dir.remote_path)


@custom_task(cpu=26, memory=700, storage_gib=2000)
def gene_stats_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    project_name: str,
    gene_stats_threads: int,
) -> LatchDir:

    import anndata

    local_results = Path(results_dir.local_path)
    dirs = _dirs_for_base(local_results)
    archrproj_path = local_results / f"{project_name}_ArchRProject"
    if not archrproj_path.exists():
        raise FileNotFoundError(
            f"Could not find ArchRProject at {archrproj_path}. "
            "Run this task after the gene and motif artifact tasks."
        )

    marker_threads = max(1, min(int(gene_stats_threads), 25))
    if marker_threads != gene_stats_threads:
        logging.warning(
            f"Clamping gene_stats_threads from {gene_stats_threads} to "
            f"{marker_threads}; this task reserves at most 50 CPUs."
        )

    output_dir = Path(f"/root/{project_name}_gene_stats")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(
        f"Running full gene differential statistics with {marker_threads} thread(s)."
    )
    _gene_stats_cmd = [
        "Rscript",
        "/root/wf/R/archr_gene_stats.R",
        project_name,
        str(archrproj_path),
        str(output_dir),
        str(marker_threads),
    ]
    _gene_stats_cmd.extend(_gene_stats_run_args(runs))
    subprocess.run(_gene_stats_cmd, check=True)

    _copy_directory_contents(output_dir / "tables", dirs["tables"])
    _copy_directory_contents(output_dir / "figures", dirs["figures"])

    adata_path = local_results / "combined_sm_ge.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(
            f"Could not find combined_sm_ge.h5ad at {adata_path}."
        )

    logging.info(f"Patching gene statistics into {adata_path}...")
    adata_gene_sm = anndata.read_h5ad(adata_path)
    ft.load_analysis_results(
        adata_gene_sm,
        "gene",
        utils.get_groups(runs),
        input_dir=dirs["tables"],
    )

    ft._sanitize_uns_for_h5ad(adata_gene_sm.uns)
    adata_gene_sm.write(adata_path)

    logging.info("Uploading results with gene differential statistics to Latch...")
    return LatchDir(str(local_results), results_dir.remote_path)


@custom_task(cpu=50, memory=512, storage_gib=1000)
def motifs_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(results_dir)
    groups = utils.get_groups(runs)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)

    # Download ArchRProject
    archrproj_path = f"{results_dir.remote_path}/{project_name}_ArchRProject"
    archrproj_path = LatchDir(archrproj_path).local_path

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

    position_files = {}
    missing_positions = []
    for run in runs:
        position_file = utils.get_LatchFile(
            run.spatial_dir, "tissue_positions_list.csv"
        )
        if position_file is None:
            missing_positions.append(f"{run.run_id} ({run.spatial_dir.remote_path})")
            continue
        position_files[run.run_id] = position_file

    if missing_positions:
        missing_str = ", ".join(missing_positions)
        raise FileNotFoundError(
            "Unable to resolve 'tissue_positions_list.csv' for one or more runs: "
            f"{missing_str}. Ensure each spatial directory contains exactly one "
            "'tissue_positions_list.csv' file."
        )

    runs = [
        (
            f'{run.run_id},'
            f'{run.fragments_file.local_path},'
            f'{utils.sanitize_condition(run.condition)},'
            f'{position_files[run.run_id].local_path},'
            f'{run.spatial_dir.local_path},'
            f'{run.sample_name}'
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
    utils.organize_outputs(project_name, dirs, exclude_pattern="*_hm.csv")

    # Save AnnData
    ft.save_anndata_objects(adata_motif, "_motifs", dirs['base'])

    logging.info("Copying ArchR peak files to top directory...")
    utils.copy_peak_files(project_name, dirs)

    logging.info("Making Plots Artifact...")
    artifact = PlotsArtifact(
        bindings=PlotsArtifactBindings(
            plot_templates=[
                PlotsArtifactTemplate(
                    template_id="760",
                    widgets=[
                        Widget(
                            transform_id="401983",
                            key="data_path",
                            value=results_dir.remote_path
                        ),
                        Widget(
                            transform_id="401994",
                            key="coverages_genome",
                            value=genome
                        )
                    ],
                )
            ]
        )
    )

    artifact_dict = artifact.asdict()

    artifacts_dir = dirs["base"] / "Launch_Plots"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    with open(artifacts_dir / "artifact.json", "w") as f:
        json.dump(artifact_dict, f, indent=2)

    logging.info("Uploading data to Latch...")
    return LatchDir(str(dirs['base']), results_dir.remote_path)


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
