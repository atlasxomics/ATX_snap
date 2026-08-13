import json
import logging
import math

from typing import List, Tuple, Union
from pathlib import Path

import anndata
import numpy as np
import snapatac2 as snap
from scipy.sparse import vstack

from wf.utils import Run, get_LatchFile, ref_dict, sanitize_condition

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def add_clusters(
    adata: anndata.AnnData,
    resolution: float,
    n_comps: int,
    leiden_iters: int,
    min_cluster_size: int,
) -> Tuple[anndata.AnnData, str]:
    """Perform dimensionality reduction, batch correction, umap, clustering."""

    # First reduce to n_comps demensions
    snap.tl.spectral(adata, n_comps=n_comps, features="selected")

    try:
        n_runs = len(adata.obs["sample"].unique())
    except KeyError as e:
        raise KeyError(
            f"Exception {e}: Please add metadata to combined AnnData."
        )

    if n_runs > 1:
        logging.info("Performing batch correction with Harmony...")
        snap.pp.harmony(adata, batch="sample", max_iter_harmony=20)
        rep = "X_spectral_harmony"
    else:
        rep = "X_spectral"

    # Add umap, nearest neightbors, clusters to .obs
    snap.tl.umap(adata, use_rep=rep)
    snap.pp.knn(adata, use_rep=rep)
    snap.tl.leiden(
        adata,
        resolution=resolution,
        n_iterations=leiden_iters,
        min_cluster_size=min_cluster_size,
        key_added="cluster",
    )

    return adata, rep


def add_metadata(
    run: Run, adata: anndata.AnnData, positions_file: Union[Path, str]
) -> anndata.AnnData:
    """Add metadata and spatial info .obs of AnnData."""
    import pandas as pd

    # Read in tissue_positions file from spatial/
    positions = pd.read_csv(positions_file, header=None)
    positions.columns = ["barcode", "on_off", "row", "col", "xcor", "ycor"]

    # Match barcodes in adata/fragments_file
    positions["barcode"] = positions["barcode"] + "-1"

    # Merge fragments file with Anndata.obs
    adata.obs["barcode"] = adata.obs.index
    adata.obs = adata.obs.merge(positions, on="barcode", how="left")

    # Set run_id, condition, sample_name
    adata.obs["sample"] = run.run_id
    adata.obs["condition"] = sanitize_condition(run.condition)
    adata.obs["sample_name"] = run.sample_name

    # Ensure obs_names unique
    adata.obs_names = [
        run_id + "#" + bc
        for run_id, bc in zip(adata.obs["sample"], adata.obs["barcode"])
    ]

    return adata


def _validate_anndata_inputs(
    adatas: List[anndata.AnnData], names: List[str]
) -> None:
    """Validate inputs before constructing a multi-sample AnnDataSet."""
    import pandas as pd

    if len(adatas) == 0:
        raise ValueError("No AnnData objects provided for combination.")
    if len(adatas) != len(names):
        raise ValueError(
            f"Received {len(adatas)} AnnData objects but {len(names)} names."
        )

    # Guard against duplicate feature names within an AnnData object.
    for i, adata in enumerate(adatas):
        idx = pd.Index(adata.var_names)
        if idx.has_duplicates:
            dupes = idx[idx.duplicated()].unique().tolist()
            preview = ", ".join(map(str, dupes[:5]))
            raise ValueError(
                "Duplicate var_names detected in input AnnData "
                f"(index {i}). Example duplicates: {preview}"
            )

    # Guard against feature mismatches across AnnData objects.
    ref_var_names = np.asarray(adatas[0].var_names)
    for i, adata in enumerate(adatas[1:], start=1):
        cur_var_names = np.asarray(adata.var_names)
        if (
            cur_var_names.shape[0] != ref_var_names.shape[0]
            or not np.array_equal(cur_var_names, ref_var_names)
        ):
            raise ValueError(
                "Inconsistent var_names across AnnData inputs. "
                f"Input 0 has {ref_var_names.shape[0]} features, "
                f"input {i} has {cur_var_names.shape[0]}. "
                "All inputs must have identical feature coordinates/order."
            )


def persist_anndata_dataset(
    adatas: List[anndata.AnnData],
    names: List[str],
    output_dir: Union[Path, str],
) -> Path:
    """Persist a multi-sample AnnDataSet without materializing its combined X.

    The returned directory is self-contained so it can be uploaded as a task
    output and reopened by a later task.  Separating these operations prevents
    the in-memory per-sample objects and the materialized combined object from
    existing at the same time.
    """
    import pandas as pd

    _validate_anndata_inputs(adatas, names)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    component_paths = [
        output_dir / f"sample_{index:04d}_backend.h5ad"
        for index in range(len(adatas))
    ]

    # Input AnnData must be backend for AnnDataSet, not in-memory.
    logging.info("Converting AnnData objects to backend...")
    adatas_be = [
        convert_tobackend(adata, path)
        for adata, path in zip(adatas, component_paths)
    ]

    logging.info("Creating AnnDataSet...")
    dataset_path = output_dir / "combined.h5ads"
    adataset = snap.AnnDataSet(
        adatas=[(name, adata) for name, adata in zip(names, adatas_be)],
        filename=dataset_path,
    )
    logging.info(f"AnnDataSet created with shape {adataset.shape}")

    # We have seen the dataset lose var_names, ensure them here.
    if len(adataset.var_names) == 0:
        logging.warning(
            "AnnDataSet lost var_names; restoring from input genomic feature names."
        )
        adataset.var_names = list(adatas[0].var_names)

    combined_obs = pd.concat([adata.obs for adata in adatas])

    # Ensure obs_names unique
    combined_obs.index = [
        run_id + "#" + bc
        for run_id, bc in zip(
            combined_obs["sample"], combined_obs["barcode"]
        )
    ]
    adataset.obs = combined_obs
    adataset.obs_names = combined_obs.index

    obs_path = output_dir / "combined_obs.pkl"
    combined_obs.to_pickle(obs_path)

    manifest = {
        "dataset": dataset_path.name,
        "obs": obs_path.name,
        "components": [path.name for path in component_paths],
    }
    with (output_dir / "manifest.json").open("w") as stream:
        json.dump(manifest, stream)

    adataset.close()
    for adata_backend in adatas_be:
        adata_backend.close()

    return dataset_path


def materialize_anndata_dataset(
    dataset_dir: Union[Path, str],
) -> anndata.AnnData:
    """Load a persisted AnnDataSet as an in-memory combined AnnData."""
    import pandas as pd

    dataset_dir = Path(dataset_dir)

    with (dataset_dir / "manifest.json").open() as stream:
        manifest = json.load(stream)

    dataset_path = dataset_dir / manifest["dataset"]
    logging.info(f"Opening persisted AnnDataSet at {dataset_path}...")
    adataset = snap.read_dataset(
        dataset_path,
        adata_files_update=dataset_dir,
    )

    logging.info("Materializing combined AnnData in memory...")
    combined_adata = adataset.to_adata()
    combined_obs = pd.read_pickle(dataset_dir / manifest["obs"])
    combined_adata.obs = combined_obs
    combined_adata.obs_names = combined_obs.index

    # AnnDataSet does not inherit .obsm; add manually :/
    component_adatas = [
        snap.read(dataset_dir / component, backed="r")
        for component in manifest["components"]
    ]
    frags = vstack(
        [adata.obsm["fragment_paired"] for adata in component_adatas],
        format="csr",
    )
    combined_adata.obsm["fragment_paired"] = frags

    for component_adata in component_adatas:
        component_adata.close()
    adataset.close()

    return combined_adata


def convert_tobackend(
    adata: anndata.AnnData, filename: Union[Path, str]
) -> anndata.AnnData:
    """Create a new backend AnnData object; necessary for creating AnnDataSet;
    saves each AnnData object to disk as .h5ad.
    """

    adata_backend = snap.AnnData(
        filename=filename,
        X=adata.X,
        obs=adata.obs,
        var=adata.var,
        uns=adata.uns,
        obsm=dict(adata.obsm),
    )
    adata_backend.obs_names = adata.obs_names

    return adata_backend


def filter_adatas(
    adatas: List[anndata.AnnData], min_tss: float = 2.0
) -> List[anndata.AnnData]:
    """Filter AnnData by on/off tissue tixels, TSS enrichment, max frag counts."""

    # Filter 'off tissue' tixels
    try:
        adatas = [adata[adata.obs["on_off"] == 1] for adata in adatas]
    except KeyError as e:
        logging.warning(f"Exception {e}: no positions data found in AnnData.obs")

    # Filter cells by tss, max_counts
    snap.pp.filter_cells(adatas, min_tsse=min_tss, min_counts=None, max_counts=1e7)

    return adatas


def make_anndatas(
    runs: List[Run], genome: str, min_frags: int
) -> List[anndata.AnnData]:
    """Basic preprocessing for snapATAC2 analysis; converts fragement_tsv.gz
    files into list of _in memory_ AnnData objects. QCs, metadata and spatial
    data are stored in AnnData.obs.
    """

    # As 'in_memory' so we can add metadata to .obs
    adatas = snap.pp.import_fragments(
        [run.fragments_file.local_path for run in runs],
        chrom_sizes=ref_dict[genome][0],
        min_num_fragments=min_frags,
        sorted_by_barcode=False,
    )

    position_files = {}
    missing_positions = []
    for run in runs:
        position_file = get_LatchFile(run.spatial_dir, "tissue_positions_list.csv")
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

    # Add run_id, condition, spatial info to .obs, TSS enrichment
    adatas = [add_metadata(run, adata, position_files[run.run_id].local_path)
              for run, adata in zip(runs, adatas)]

    # Add addtional QCs
    snap.metrics.tsse(adatas, ref_dict[genome][1])

    for adata in adatas:
        if min_frags == 0:  # Convert 0 to NA
            logging.warning("Converting 0's to NA in .obs['n_fragment']")
            adata.obs["n_fragment"] = adata.obs["n_fragment"].apply(
                lambda x: np.nan if x <= 0 else x
            )

        adata.obs["log10_frags"] = adata.obs["n_fragment"].apply(math.log10)

    return adatas
