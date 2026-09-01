import fnmatch
import glob
import json
import logging
import os
import shutil
import subprocess

from pathlib import Path
from typing import Dict, List, Set, Optional

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


# Object files are organized into dedicated subfolders under the project
# directory: Seurat .rds objects in seurat_objects/, AnnData .h5ad objects in
# anndata/.
SEURAT_SUBDIR = "seurat_objects"
ANNDATA_SUBDIR = "anndata"


def _seurat_dir(base_dir: Path) -> Path:
    d = Path(base_dir) / SEURAT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def _anndata_dir(base_dir: Path) -> Path:
    d = Path(base_dir) / ANNDATA_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def _rebuild_gene_h5ad_from_rds(rds_path: Path, h5ad_path: Path) -> None:
    """Recreate a legacy gene H5AD without rerunning ArchR imputation."""
    import tempfile

    rds_path = Path(rds_path)
    h5ad_path = Path(h5ad_path)
    if not rds_path.is_file():
        raise FileNotFoundError(
            f"Cannot rebuild {h5ad_path.name}; Seurat RDS is missing: {rds_path}"
        )

    with tempfile.TemporaryDirectory(
        prefix="gene_h5ad_rebuild_",
        dir=h5ad_path.parent,
    ) as temporary_dir:
        rebuilt_path = Path(temporary_dir) / h5ad_path.name
        subprocess.run(
            [
                "Rscript",
                "/root/wf/R/rebuild_gene_h5ad.R",
                str(rds_path),
                str(rebuilt_path),
            ],
            check=True,
        )
        if not rebuilt_path.is_file():
            raise FileNotFoundError(
                f"RDS conversion did not create expected H5AD: {rebuilt_path}"
            )
        os.replace(rebuilt_path, h5ad_path)


def _ensure_dense_h5ad_x(
    h5ad_path: Path,
    recovery_rds_path: Optional[Path] = None,
) -> None:
    """Ensure a per-sample gene H5AD stores X as a dense HDF5 dataset."""
    import anndata
    import h5py
    import scipy.sparse as sparse

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.is_file():
        if recovery_rds_path is not None:
            logging.warning(
                f"Rebuilding missing H5AD from Seurat RDS: {h5ad_path.name}"
            )
            _rebuild_gene_h5ad_from_rds(recovery_rds_path, h5ad_path)
            _ensure_dense_h5ad_x(h5ad_path)
            return
        raise FileNotFoundError(f"Expected gene H5AD was not created: {h5ad_path}")

    try:
        with h5py.File(h5ad_path, "r") as handle:
            x_node = handle.get("X")
            if x_node is None:
                raise ValueError(f"H5AD has no X matrix: {h5ad_path}")
            if isinstance(x_node, h5py.Dataset):
                logging.info(f"Verified dense H5AD X: {h5ad_path.name}")
                return
    except (OSError, ValueError) as error:
        if recovery_rds_path is None:
            raise
        logging.warning(
            f"Rebuilding unreadable H5AD from Seurat RDS: {h5ad_path.name}: "
            f"{error}"
        )
        _rebuild_gene_h5ad_from_rds(recovery_rds_path, h5ad_path)
        _ensure_dense_h5ad_x(h5ad_path)
        return

    if recovery_rds_path is not None:
        logging.warning(
            f"Rebuilding legacy sparse H5AD from Seurat RDS: {h5ad_path.name}"
        )
        _rebuild_gene_h5ad_from_rds(recovery_rds_path, h5ad_path)
        _ensure_dense_h5ad_x(h5ad_path)
        return

    logging.info(f"Converting sparse H5AD X to dense: {h5ad_path.name}")
    adata = anndata.read_h5ad(h5ad_path)
    if sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    # Older SeuratDisk files can contain reserved `_index` columns and a raw
    # copy of the same gene assay. The downstream workflow does not use raw and
    # already removes it before final output, so discard that duplicate here.
    ft.clean_index_columns(adata)
    if adata.raw is not None:
        logging.info(f"Dropping duplicate raw matrix from {h5ad_path.name}")
        adata.raw = None

    temporary_path = h5ad_path.with_name(
        f".{h5ad_path.stem}.dense.tmp.h5ad"
    )
    try:
        adata.write_h5ad(temporary_path, compression="gzip")
        os.replace(temporary_path, h5ad_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()

    with h5py.File(h5ad_path, "r") as handle:
        if not isinstance(handle["X"], h5py.Dataset):
            raise ValueError(f"Unable to store H5AD X densely: {h5ad_path}")


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


def _move_directory_delta(
    source_dir: Path,
    destination_dir: Path,
    skip_relative_paths: Set[Path],
) -> None:
    """Stage generated files without retaining a second local copy."""
    destination_dir.mkdir(parents=True, exist_ok=True)

    if not source_dir.exists():
        logging.warning(f"Directory does not exist and was not staged: {source_dir}")
        return

    for source in source_dir.rglob("*"):
        if not source.is_file():
            continue

        relative_path = source.relative_to(source_dir)
        if relative_path in skip_relative_paths:
            continue

        destination = destination_dir / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))


def _copy_relative_files(
    source_dir: Path,
    destination_dir: Path,
    relative_paths: List[Path],
) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)

    for relative_path in relative_paths:
        source = source_dir / relative_path
        if not source.exists():
            logging.warning(f"Expected stage output does not exist: {source}")
            continue

        destination = destination_dir / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source, destination)
            continue

        shutil.copy2(source, destination)


def _copy_archr_reports(
    archr_project_dir: Path,
    dirs: Dict[str, Path],
) -> None:
    """Promote small ArchR reports before their checkpoint is cleaned up."""
    archr_project_dir = Path(archr_project_dir)
    report_types = {
        "*.csv": dirs["tables"],
        "*.pdf": dirs["figures"],
    }

    for pattern, destination_dir in report_types.items():
        destination_dir.mkdir(parents=True, exist_ok=True)
        for source in sorted(archr_project_dir.rglob(pattern)):
            if not source.is_file():
                continue

            destination = destination_dir / source.name
            if destination.exists():
                # Preserve both reports if different ArchR subdirectories use
                # the same basename, while keeping normal filenames unchanged.
                relative_parent = source.parent.relative_to(archr_project_dir)
                parent_prefix = "__".join(relative_parent.parts)
                if parent_prefix:
                    destination = destination_dir / f"{parent_prefix}__{source.name}"

            shutil.copy2(source, destination)
            logging.info(f"Preserved ArchR report: {destination}")


def _fresh_stage_dir(project_name: str, stage_name: str) -> Path:
    stage_dir = Path(f"/root/{project_name}_{stage_name}_stage")
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    return stage_dir


def _gene_stats_run_args(runs: List[utils.Run]) -> List[str]:
    return [
        f"{run.run_id},{run.sample_name}"
        for run in runs
    ]


def _gene_project_run_args(runs: List[utils.Run]) -> List[str]:
    """Build ArchR run arguments without downloading unused spatial inputs."""
    return [
        (
            f"{run.run_id},"
            f"{run.fragments_file.local_path},"
            f"{utils.sanitize_condition(run.condition)},"
            f"unused,unused,"
            f"{run.sample_name}"
        )
        for run in runs
    ]


def _gene_export_run_args(runs: List[utils.Run]) -> List[str]:
    """Build run arguments needed to construct per-sample spatial objects."""
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

    return [
        (
            f"{run.run_id},unused,"
            f"{utils.sanitize_condition(run.condition)},"
            f"{position_files[run.run_id].local_path},"
            f"{run.spatial_dir.local_path},"
            f"{run.sample_name}"
        )
        for run in runs
    ]


def _unique_glob_matches(
    patterns: List[str],
    exclude_pattern: Optional[str] = None,
) -> List[str]:
    matches = []
    seen = set()

    for pattern in patterns:
        for match in glob.glob(pattern):
            if exclude_pattern and fnmatch.fnmatch(match, exclude_pattern):
                continue
            if match in seen:
                continue

            seen.add(match)
            matches.append(match)

    return matches


def _move_files_to_directory(files: List[str], target_dir: Path) -> None:
    if files:
        subprocess.run(["mv"] + files + [str(target_dir)], check=True)


def _organize_outputs(
    project_name: str,
    dirs: Dict[str, Path],
    exclude_pattern: Optional[str] = None,
) -> None:
    logging.info("Moving outputs to output directory...")

    project_dirs = _unique_glob_matches([f"{project_name}_*"])
    _move_files_to_directory(project_dirs, dirs["base"])

    # Route object files into dedicated subfolders
    rds_files = _unique_glob_matches(["*.rds"])
    if rds_files:
        _move_files_to_directory(rds_files, _seurat_dir(dirs["base"]))

    h5ad_files = _unique_glob_matches(["*.h5ad"])
    if h5ad_files:
        _move_files_to_directory(h5ad_files, _anndata_dir(dirs["base"]))

    csv_files = _unique_glob_matches(["*.csv"], exclude_pattern=exclude_pattern)
    _move_files_to_directory(csv_files, dirs["tables"])

    figures = [
        figure for figure in _unique_glob_matches(["*.pdf"])
        if figure != "Rplots.pdf"
    ]
    _move_files_to_directory(figures, dirs["figures"])


@custom_task(cpu=32, memory=512, storage_gib=1000)
def make_anndata_dataset_task(
    runs: List[utils.Run],
    genome: utils.Genome,
    project_name: str,
    min_tss: float,
    min_frags: int,
    include_y_chromosome: bool,
    tile_size: int,
    output_dir: LatchDir,
) -> LatchDir:
    """Create and upload the backed multi-sample dataset as its own stage."""
    genome_name = genome.value
    stage_dir = _fresh_stage_dir(project_name, "anndata_dataset")

    logging.info("Creating AnnData objects...")
    adatas = pp.make_anndatas(runs, genome_name, min_frags=min_frags)
    adatas = pp.filter_adatas(adatas, min_tss=min_tss)

    logging.info("Adding tile matrix to objects...")
    excluded_chroms = ["chrM", "M"]
    if not include_y_chromosome:
        excluded_chroms.extend(["chrY", "Y"])
    snap.pp.add_tile_matrix(
        adatas,
        bin_size=tile_size,
        exclude_chroms=excluded_chroms,
    )

    samples = [run.run_id for run in runs]
    logging.info("Persisting combined AnnDataSet before materialization...")
    pp.persist_anndata_dataset(adatas, samples, stage_dir)

    remote_path = (
        f"{output_dir.remote_path.rstrip('/')}/{project_name}"
        "/checkpoints/anndata_dataset"
    )
    return LatchDir(str(stage_dir), remote_path)


@custom_task(cpu=48, memory=512, storage_gib=2000)
def make_adata(
    runs: List[utils.Run],
    anndata_dataset: LatchDir,
    genome: utils.Genome,
    project_name: str,
    resolution: float,
    leiden_iters: int,
    n_comps: int,
    min_cluster_size: int,
    min_tss: float,
    min_frags: int,
    include_y_chromosome: bool,
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
        ["include y chromosome", include_y_chromosome],
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

    # Materialization happens in a separate task from construction of the
    # AnnDataSet, so the per-sample in-memory objects have already been freed.
    # Keep the combined object in memory for compatibility with the existing
    # analysis code; the task boundary still reduces the peak memory footprint.
    combined_path = _anndata_dir(Path(result_dir)) / "combined.h5ad"
    adata = pp.materialize_anndata_dataset(
        anndata_dataset.local_path,
    )

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

    ft.add_spatial_offset(adata)
    adata.write(combined_path)

    return LatchDir(result_dir, output_dir), groups


@custom_task(cpu=30, memory=150, storage_gib=2000)
def gene_project_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    include_y_chromosome: bool,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(results_dir)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)
    _copy_required_input_tables(data_paths, dirs["tables"])

    logging.info("Building checkpointed ArchR gene project...")
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_genes.R',
        project_name,
        genome,
        data_paths['obs'],
        data_paths['spectral'],
        str(include_y_chromosome).lower(),
        'prepare',
    ]
    _archr_cmd.extend(_gene_project_run_args(runs))
    subprocess.run(_archr_cmd, check=True)

    # Persist the project and its on-disk imputation weights before any dense
    # GeneScoreMatrix extraction begins.
    _organize_outputs(project_name, dirs)

    checkpoint_file = (
        dirs["base"]
        / f"{project_name}_ArchRProject"
        / "Save-ArchR-Project.rds"
    )
    if not checkpoint_file.is_file():
        raise FileNotFoundError(
            f"ArchR gene-project checkpoint was not created at {checkpoint_file}."
        )

    delta_dir = _fresh_stage_dir(project_name, "gene_project_delta")
    _move_directory_delta(
        dirs["base"],
        delta_dir,
        {
            Path("tables/obs.csv"),
            Path("tables/spatial.csv"),
            Path("tables/X_umap.csv"),
            Path("tables/spectral.csv"),
        },
    )

    logging.info("Uploading the ArchR gene-project checkpoint to Latch...")
    return LatchDir(str(delta_dir), results_dir.remote_path)


@custom_task(cpu=24, memory=112, storage_gib=2000)
def genes_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    gene_project_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    include_y_chromosome: bool,
) -> LatchDir:

    data_paths = utils.get_data_paths(results_dir)
    genome = genome.value

    dirs = utils.create_output_directories(project_name)
    _copy_required_input_tables(data_paths, dirs["tables"])

    logging.info("Exporting checkpointed ArchR gene scores in bounded chunks...")
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_genes.R',
        project_name,
        genome,
        data_paths['obs'],
        data_paths['spectral'],
        str(include_y_chromosome).lower(),
        'export',
        gene_project_dir.local_path,
    ]
    _archr_cmd.extend(_gene_export_run_args(runs))
    subprocess.run(_archr_cmd, check=True)

    # Seurat v5 stores assay counts sparsely. Convert each exported H5AD back
    # to a dense X dataset and verify it before uploading the gene artifacts.
    for run in runs:
        _ensure_dense_h5ad_x(Path(f"{run.run_id}_g_converted.h5ad"))

    # Stage Seurat objects, per-run h5ads, and R-side tables. The ArchRProject
    # has already been persisted by gene_project_task at the same remote root.
    _organize_outputs(project_name, dirs)

    missing_outputs = []
    seurat_dir = _seurat_dir(dirs["base"])
    anndata_dir = _anndata_dir(dirs["base"])
    for run in runs:
        rds_path = seurat_dir / f"{run.run_id}_SeuratObj.rds"
        h5ad_path = anndata_dir / f"{run.run_id}_g_converted.h5ad"
        if not rds_path.is_file():
            missing_outputs.append(str(rds_path))
        if not h5ad_path.is_file():
            missing_outputs.append(str(h5ad_path))

    if missing_outputs:
        raise FileNotFoundError(
            "Gene export completed without all expected per-run outputs: "
            + ", ".join(missing_outputs)
        )

    delta_dir = _fresh_stage_dir(project_name, "genes_delta")
    _move_directory_delta(
        dirs["base"],
        delta_dir,
        {
            Path("tables/obs.csv"),
            Path("tables/spatial.csv"),
            Path("tables/X_umap.csv"),
            Path("tables/spectral.csv"),
        },
    )

    logging.info("Uploading gene-stage artifacts to Latch...")
    return LatchDir(str(delta_dir), results_dir.remote_path)


@custom_task(cpu=6, memory=300, storage_gib=2000)
def combine_gene_h5ads_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    gene_results_dir: LatchDir,
    project_name: str,
) -> LatchDir:
    # Continue in the uploaded gene artifact directory so a combine failure
    # does not discard the already-materialized Seurat and per-run h5ad files.
    dirs = _dirs_for_base(Path(gene_results_dir.local_path))

    # Load and combine data (per-run h5ads live in the anndata/ subfolder)
    anndata_dir = _anndata_dir(dirs["base"])
    seurat_dir = _seurat_dir(dirs["base"])
    for run in runs:
        _ensure_dense_h5ad_x(
            anndata_dir / f"{run.run_id}_g_converted.h5ad",
            recovery_rds_path=seurat_dir / f"{run.run_id}_SeuratObj.rds",
        )

    adata_gene = ft.load_and_combine_data(
        "g_converted",
        input_dir=anndata_dir,
        temp_dir=anndata_dir,
    )

    # Persist the expensive repair/combination result before any Squidpy work.
    # The final gene-expression outputs are already float32, so cast at this
    # checkpoint to halve its size without changing downstream precision.
    ft._cast_X_dtype(adata_gene, "float32", "combined gene checkpoint")
    checkpoint_dir = _fresh_stage_dir(project_name, "gene_combined")
    checkpoint_path = (
        _anndata_dir(checkpoint_dir) / "combined_g_pre_squidpy.h5ad"
    )
    logging.info(f"Saving dense pre-Squidpy checkpoint to {checkpoint_path}...")
    adata_gene.write_h5ad(checkpoint_path)
    if not checkpoint_path.is_file():
        raise FileNotFoundError(
            f"Dense pre-Squidpy checkpoint was not created: {checkpoint_path}"
        )

    logging.info("Uploading dense combined gene checkpoint...")
    checkpoint_remote_path = (
        f"{results_dir.remote_path.rstrip('/')}/checkpoints/gene_combined"
    )
    return LatchDir(str(checkpoint_dir), checkpoint_remote_path)


@custom_task(cpu=16, memory=600, storage_gib=4000)
def gene_spatial_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    gene_results_dir: LatchDir,
    gene_combined_dir: LatchDir,
    project_name: str,
) -> LatchDir:
    import anndata

    checkpoint_path = (
        Path(gene_combined_dir.local_path)
        / ANNDATA_SUBDIR
        / "combined_g_pre_squidpy.h5ad"
    )
    if not checkpoint_path.is_file():
        raise FileNotFoundError(
            f"Could not find dense pre-Squidpy checkpoint at {checkpoint_path}."
        )

    logging.info(f"Loading dense pre-Squidpy checkpoint from {checkpoint_path}...")
    adata_gene = anndata.read_h5ad(checkpoint_path)
    groups = utils.get_groups(runs)
    data_paths = utils.get_data_paths(results_dir)
    dirs = _dirs_for_base(_fresh_stage_dir(project_name, "gene_spatial_work"))
    anndata_dir = _anndata_dir(dirs["base"])

    # Cell-ID alignment and auxiliary transfer are intentionally downstream of
    # the durable dense checkpoint so failures here never repeat combination.
    ft.transfer_auxiliary_data(adata_gene, data_paths, groups)
    missing_obsm = [
        key for key in ("X_umap", "spatial") if key not in adata_gene.obsm
    ]
    if missing_obsm:
        raise KeyError(
            "Dense gene checkpoint is missing required auxiliary coordinates: "
            + ", ".join(missing_obsm)
        )
    if not adata_gene.obs_names.is_unique:
        raise ValueError(
            "Gene observation names remain non-unique after auxiliary alignment."
        )

    # Run spatial analysis
    sample_key = "sample" if "sample" in adata_gene.obs.columns else None
    adata_gene = sp.run_squidpy_analysis(
        adata_gene,
        dirs["figures"],
        sample_key,
        group_keys=groups,
    )

    # Spatially variable genes
    samples = [run.run_id for run in runs]
    try:
        svg_df = sp.run_spatial_autocorr(adata_gene, n_jobs=4)
        svg_df.to_csv(dirs["tables"] / "svg_genes.csv")
        pl.plot_svg_spatial(
            adata_gene,
            svg_df,
            samples,
            str(dirs["figures"] / "svg_spatial_genes.png"),
            modality="Genes",
            top_n=10,
            html_output_path=str(dirs["figures"] / "svg_spatial_genes.html"),
        )
    except Exception as e:
        warning = f"Spatial autocorrelation (genes) failed: {e}"
        logging.warning(warning)
        message(typ="warning", data={"title": "SVG genes failed", "body": warning})

    # Load differential analysis results
    ft.load_analysis_results(
        adata_gene,
        "gene",
        groups,
        input_dir=Path(gene_results_dir.local_path) / "tables",
    )
    # Save AnnData (combined objects go into the anndata/ subfolder)
    ft.save_anndata_objects(
        adata_gene,
        "_ge",
        anndata_dir,
        full_x_dtype="float32",
    )

    delta_dir = _fresh_stage_dir(project_name, "gene_expression_delta")

    svg_table_paths = [
        p.relative_to(dirs["base"]) for p in dirs["tables"].glob("svg_*.csv")
    ]
    svg_figure_paths = [
        p.relative_to(dirs["base"]) for p in dirs["figures"].glob("svg_spatial_genes*")
    ]

    _copy_relative_files(
        dirs["base"],
        delta_dir,
        [
            Path(ANNDATA_SUBDIR) / "combined_ge.h5ad",
            Path(ANNDATA_SUBDIR) / "combined_sm_ge.h5ad",
            Path("figures/all_neighborhoods.pdf"),
        ] + svg_table_paths + svg_figure_paths,
    )

    logging.info("Uploading gene-expression stage artifacts to Latch...")
    return LatchDir(str(delta_dir), results_dir.remote_path)


@custom_task(cpu=26, memory=500, storage_gib=4000)
def gene_stats_task(
    runs: List[utils.Run],
    gene_results_dir: LatchDir,
    gene_expression_results_dir: LatchDir,
    results_root: LatchDir,
    project_name: str,
) -> LatchDir:

    import anndata

    local_gene_results = Path(gene_results_dir.local_path)
    local_gene_expression = Path(gene_expression_results_dir.local_path)
    archrproj_path = local_gene_results / f"{project_name}_ArchRProject"
    if not archrproj_path.exists():
        raise FileNotFoundError(
            f"Could not find ArchRProject at {archrproj_path}. "
            "Run this task after the gene artifact task."
        )

    source_adata_path = (
        local_gene_expression / ANNDATA_SUBDIR / "combined_sm_ge.h5ad"
    )
    if not source_adata_path.exists():
        raise FileNotFoundError(
            f"Could not find combined_sm_ge.h5ad at {source_adata_path}."
        )

    delta_dir = _fresh_stage_dir(project_name, "gene_stats_delta")
    dirs = _dirs_for_base(delta_dir)

    marker_threads = 25

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

    adata_path = _anndata_dir(delta_dir) / "combined_sm_ge.h5ad"
    shutil.copy2(source_adata_path, adata_path)

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
    return LatchDir(str(delta_dir), results_root.remote_path)


@custom_task(cpu=8, memory=32, storage_gib=1000)
def motif_coverages_task(
    gene_results_dir: LatchDir,
    project_name: str,
) -> LatchDir:
    """Generate and persist the cluster group-coverage checkpoint."""
    archrproj_remote_path = (
        f"{gene_results_dir.remote_path.rstrip('/')}"
        f"/{project_name}_ArchRProject"
    )
    archrproj_path = Path(LatchDir(archrproj_remote_path).local_path)
    if not archrproj_path.exists():
        raise FileNotFoundError(
            f"Could not find ArchRProject at {archrproj_path}."
        )

    checkpoint_dir = _fresh_stage_dir(project_name, "motif_coverages")
    checkpoint_project = checkpoint_dir / "ArchRProject"
    subprocess.run(
        [
            "Rscript",
            "/root/wf/R/archr_motif_coverages.R",
            str(archrproj_path),
            str(checkpoint_project),
            "8",
        ],
        check=True,
    )
    checkpoint_rds = checkpoint_project / "Save-ArchR-Project.rds"
    if not checkpoint_project.is_dir() or not checkpoint_rds.is_file():
        raise FileNotFoundError(
            "ArchR coverage checkpoint is incomplete; expected project file "
            f"at {checkpoint_rds}."
        )

    checkpoint_remote_path = (
        f"{gene_results_dir.remote_path.rstrip('/')}"
        "/checkpoints/motifs/coverages"
    )
    logging.info("Uploading motif group-coverage checkpoint...")
    return LatchDir(str(checkpoint_dir), checkpoint_remote_path)


@custom_task(cpu=8, memory=64, storage_gib=1000)
def motif_peaks_task(
    motif_coverages_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    include_y_chromosome: bool,
) -> LatchDir:
    """Create peaks and motif annotations from the coverage checkpoint."""
    coverage_project = Path(motif_coverages_dir.local_path) / "ArchRProject"
    if not coverage_project.is_dir():
        raise FileNotFoundError(
            "Could not find the coverage-checkpoint ArchRProject at "
            f"{coverage_project}."
        )

    genome_sizes = {
        "hg38": 3.3e9,
        "mm10": 3.0e9,
        "mm39": 2.7e9,
        "rnor6": 2.9e9,
    }
    genome_name = genome.value
    if genome_name not in genome_sizes:
        raise ValueError(f"No genome size configured for genome: {genome_name}")

    checkpoint_dir = _fresh_stage_dir(project_name, "motif_peaks")
    checkpoint_project = checkpoint_dir / "ArchRProject"
    subprocess.run(
        [
            "Rscript",
            "/root/wf/R/archr_motif_peaks.R",
            str(coverage_project),
            str(checkpoint_project),
            genome_name,
            str(include_y_chromosome).lower(),
            str(genome_sizes[genome_name]),
            "8",
        ],
        check=True,
    )
    checkpoint_rds = checkpoint_project / "Save-ArchR-Project.rds"
    if not checkpoint_project.is_dir() or not checkpoint_rds.is_file():
        raise FileNotFoundError(
            "Annotated peak checkpoint is incomplete; expected project file "
            f"at {checkpoint_rds}."
        )

    coverages_remote_path = motif_coverages_dir.remote_path.rstrip("/")
    checkpoint_root = coverages_remote_path.rsplit("/", 1)[0]
    checkpoint_remote_path = f"{checkpoint_root}/peaks"
    logging.info("Uploading annotated motif peak checkpoint...")
    return LatchDir(str(checkpoint_dir), checkpoint_remote_path)


@custom_task(cpu=8, memory=264, storage_gib=1000)
def motifs_task(
    runs: List[utils.Run],
    results_dir: LatchDir,
    motif_peaks_dir: LatchDir,
    project_name: str,
    genome: utils.Genome,
    include_y_chromosome: bool,
) -> LatchDir:

    # Read in data tables
    data_paths = utils.get_data_paths(results_dir)
    groups = utils.get_groups(runs)

    genome = genome.value

    # Create output dirs
    dirs = utils.create_output_directories(project_name)

    # Load the durable ArchR checkpoint containing cluster peaks and motif
    # annotations. Coverage generation and peak calling run in prior tasks.
    archrproj_path = Path(motif_peaks_dir.local_path) / "ArchRProject"
    if not archrproj_path.is_dir():
        raise FileNotFoundError(
            f"Could not find annotated motif ArchRProject at {archrproj_path}."
        )

    logging.info("Running ArchR analysis...")
    # run subprocess R script to make .h5ad file
    _archr_cmd = [
        'Rscript',
        '/root/wf/R/archr_motifs.R',
        project_name,
        genome,
        data_paths['obs'],
        str(archrproj_path),
        str(include_y_chromosome).lower(),
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

    run_args = [
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
    _archr_cmd.extend(run_args)
    subprocess.run(_archr_cmd, check=True)

    # Load and combine data
    adata_motif = ft.load_and_combine_data("m_converted")

    # Transfer auxiliary data to combined AnnData
    ft.transfer_auxiliary_data(adata_motif, data_paths, groups)

    # Spatially variable motifs
    samples = [run.run_id for run in runs]
    try:
        svg_df_motifs = sp.run_spatial_autocorr(adata_motif, n_jobs=4)
        svg_df_motifs.to_csv(dirs["tables"] / "svg_motifs.csv")
        pl.plot_svg_spatial(
            adata_motif,
            svg_df_motifs,
            samples,
            str(dirs["figures"] / "svg_spatial_motifs.png"),
            modality="Motifs",
            top_n=10,
            html_output_path=str(dirs["figures"] / "svg_spatial_motifs.html"),
        )
    except Exception as e:
        warning = f"Spatial autocorrelation (motifs) failed: {e}"
        logging.warning(warning)
        message(typ="warning", data={"title": "SVG motifs failed", "body": warning})

    # Load differential analysis results
    ft.load_analysis_results(adata_motif, "motif", groups)

    # Organize outputs
    _organize_outputs(project_name, dirs, exclude_pattern="*_hm.csv")

    # Peak calling runs in a checkpoint task. Preserve its small summary tables
    # and figures in the final output before successful-workflow cleanup removes
    # the large checkpoint project.
    _copy_archr_reports(archrproj_path, dirs)

    # Save AnnData (combined objects go into the anndata/ subfolder)
    ft.save_anndata_objects(adata_motif, "_motifs", _anndata_dir(dirs['base']))

    logging.info("Copying ArchR peak files to top directory...")
    utils.copy_peak_files(project_name, dirs)

    logging.info("Making Plots Artifact...")
    artifact = PlotsArtifact(
        bindings=PlotsArtifactBindings(
            plot_templates=[
                PlotsArtifactTemplate(
                    template_id="1090",
                    widgets=[
                        Widget(
                            transform_id="433333",
                            key="data_path",
                            value=results_dir.remote_path
                        ),
                        Widget(
                            transform_id="433323",
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
def complete_results_task(
    base_results_dir: LatchDir,
    gene_results_dir: LatchDir,
    gene_expression_results_dir: LatchDir,
    gene_stats_results_dir: LatchDir,
    motif_results_dir: LatchDir,
) -> LatchDir:
    _ = (
        gene_results_dir,
        gene_expression_results_dir,
        gene_stats_results_dir,
        motif_results_dir,
    )
    return base_results_dir


@small_task(cache=False)
def cleanup_checkpoints_task(results: LatchDir) -> LatchDir:
    """Remove durable intermediates after every result branch succeeds."""
    from latch.ldata.path import LPath

    results_remote_path = results.remote_path
    if results_remote_path is None:
        raise ValueError("Cannot clean checkpoints without a remote results path.")

    results_remote_path = results_remote_path.rstrip("/")
    if not results_remote_path.startswith("latch://"):
        raise ValueError(
            "Checkpoint cleanup only supports Latch Data paths; got "
            f"{results_remote_path}."
        )

    checkpoint_remote_path = f"{results_remote_path}/checkpoints"
    if checkpoint_remote_path in {"latch://checkpoints", "latch:///checkpoints"}:
        raise ValueError(
            f"Refusing to clean unsafe checkpoint path: {checkpoint_remote_path}"
        )

    checkpoint_path = LPath(checkpoint_remote_path)
    if checkpoint_path.exists():
        logging.info(
            "Removing completed workflow checkpoints from "
            f"{checkpoint_remote_path}..."
        )
        checkpoint_path.rmr()
        logging.info("Checkpoint cleanup complete.")
    else:
        logging.info(
            f"No workflow checkpoints found at {checkpoint_remote_path}."
        )

    # Reuse the existing remote directory rather than uploading a local copy.
    return LatchDir(results_remote_path)


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
