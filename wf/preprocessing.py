import anndata
import logging
import math
import numpy as np
import pandas as pd
import snapatac2 as snap

from scipy.sparse import vstack
from typing import List

from wf.utils import Genome, Run, refresh_adata


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


def add_clusters(
    adata: anndata.AnnData,
    resolution: float,
    leiden_iters: int,
    min_cluster_size: int
) -> anndata.AnnData:
    """Perform dimensionality reduction, batch correction, umap, clustering.
    """

    # Dimensionality reduction
    snap.tl.spectral(adata)

    try:
        n_runs = len(adata.obs["sample"].unique())
    except KeyError as e:
        logging.warning(
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
        key_added="cluster"
    )

    return adata


def add_metadata(run: Run, adata: anndata.AnnData) -> anndata.AnnData:
    """Add metadata and spatial info .obs of AnnData.
    """

    # Read in tissue_positions file from spatial/
    positions = pd.read_csv(run.positions_file, header=None)
    positions.columns = ["barcode", "on_off", "row", "col", "xcor", "ycor"]

    # Match barcodes in adata/fragments_file
    positions["barcode"] = positions["barcode"] + "-1"
    positions.set_index("barcode", inplace=True)

    # Merge fragments file with Anndata.obs
    adata.obs["barcode"] = adata.obs_names

    obs_names_series = pd.Series(adata.obs_names, index=adata.obs_names)
    for col in positions.columns:
        adata.obs[col] = obs_names_series.map(positions[col])

    # Set run_id, condition
    adata.obs["sample"] = [run.run_id] * len(adata.obs_names)
    adata.obs["condition"] = [run.condition] * len(adata.obs_names)

    # Ensure obs_names unique
    adata.obs_names = [
        run_id + "#" + bc for
        run_id, bc in zip(adata.obs["sample"], adata.obs["barcode"])
    ]

    return adata


def combine_anndata(
    adatas: List[anndata.AnnData],
    names: List[str],
    filename: str = "combined"
) -> anndata.AnnData:
    """Combines a list of AnnData objects into a combined AnnData object.
    Converts in-memory AnnData to backend, saving to disk as .h5ad.
    Combines as AnnDataSet (object written to disk as .h5ad), then back to
    AnnData.
    """

    # obs/obsm not inherited automatically, need to be set at beginning
    obs = adatas[0].to_memory().obs.columns
    frags = vstack([adata.obsm["fragment_paired"] for adata in adatas])

    logging.info("Creating AnnDataSet...")
    adataset = snap.AnnDataSet(
        adatas=[(name, adata) for name, adata in zip(names, adatas)],
        filename=f"{filename}.h5ad"
    )
    logging.info(f"AnnDataSet created with shape {adataset.shape}")

    # We have seen the dataset lose var_names, ensure them here.
    if len(adataset.var_names) == 0:
        adataset.var_names = [str(i) for i in range(len(adatas[0].var_names))]

    # AnnDataSet does inherit .obs; add manually :/
    for ob in obs:
        adataset.obs[ob] = adataset.adatas.obs[ob]

    # Ensure obs_names unique
    adataset.obs_names = [
        run_id + "#" + bc for
        run_id, bc in zip(adataset.obs["sample"], adataset.obs["barcode"])
    ]

    adataset.obsm["fragment_paired"] = frags

    return adataset


def filter_adatas(
    adatas: List[anndata.AnnData], min_tss: float = 2.0
) -> List[anndata.AnnData]:
    """Filter AnnData by on/off tissue tixels, TSS enrichment, max frag counts.
    """

    # Filter 'off tissue' tixels
    print("filtering")
    for adata in adatas:
        if "on_off" in adata.obs:
            try:
                obs_indices = (adata.obs["on_off"] == 1).to_numpy().nonzero()[0]
                adata.subset(obs_indices=obs_indices)
            except Exception as e:
                logging.warning(f"Exception {e}: Error while subsetting AnnData.")
        else:
            logging.warning("No 'on_off' data found in AnnData.obs.")

    # Filter cells by tss, max_counts
    snap.pp.filter_cells(
        adatas, min_tsse=min_tss, min_counts=None, max_counts=1e7
    )

    return adatas


def make_anndatas(
    runs: List[Run],
    genome: Genome,
    min_frags: int
) -> List[anndata.AnnData]:
    """Basic preprocessing for snapATAC2 analysis; converts fragement_tsv.gz
    files into list of _in memory_ AnnData objects. QCs, metadata and spatial
    data are stored in AnnData.obs.
    """

    # Can't use a dict because of flyte
    genome_ref = snap.genome.mm10 if genome == "mm10" else snap.genome.hg38

    adatas = snap.pp.import_data(
        [run.fragments_file.local_path for run in runs],
        chrom_sizes=genome_ref,
        min_num_fragments=min_frags,
        sorted_by_barcode=False,
        file=[f"{run.run_id}.h5ad" for run in runs]
    )

    # Read back into memory to ensure write access, it's dumb
    adatas = [
        refresh_adata(adata, run.run_id) for adata, run in zip(adatas, runs)
    ]

    # Add run_id, condition, spatial info to .obs, TSS enrichment
    adatas = [add_metadata(run, adata) for run, adata in zip(runs, adatas)]

    # Add addtional QCs
    snap.metrics.tsse(adatas, genome_ref)

    adatas = [
        refresh_adata(adata, run.run_id) for adata, run in zip(adatas, runs)
    ]

    for adata in adatas:

        if min_frags == 0:  # Convert 0 to NA
            logging.warning("Converting 0's to NA in .obs['n_fragment']")
            adata.obs["n_fragment"] = pd.Series(adata.obs["n_fragment"]).apply(
                lambda x: np.nan if x <= 0 else x
            )

        adata.obs["log10_frags"] = pd.Series(adata.obs["n_fragment"]).apply(math.log10)

    return adatas
