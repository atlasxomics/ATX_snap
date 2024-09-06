import anndata
import logging
import math
import pandas as pd
import snapatac2 as snap

from scipy.sparse import vstack
from typing import List

from wf.utils import Genome, Run


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


def add_clusters(adata: anndata.AnnData) -> anndata.AnnData:
    """Perform dimensionality reduction, batch correction, umap, clustering.
    """

    # Dimensionality reduction
    logging.info("Performing dimensionality reduction...")
    snap.tl.spectral(adata)

    try:
        n_runs = len(adata.obs["run_id"].unique())
    except KeyError as e:
        logging.warn(
            f"Exception {e}: Please add metadata to combined AnnData."
        )

    if n_runs > 1:
        logging.info("Performing batch correction with Harmony...")
        snap.pp.harmony(adata, batch="run_id", max_iter_harmony=20)
        rep = "X_spectral_harmony"
    else:
        rep = "X_spectral"

    # Add umap, nearest neightbors, clusters to .obs
    snap.tl.umap(adata, use_rep=rep)
    snap.pp.knn(adata, use_rep=rep)
    snap.tl.leiden(adata)

    return adata


def add_metadata(run: Run, adata: anndata.AnnData) -> anndata.AnnData:
    """Add metadata and spatial info .obs of AnnData.
    """

    # Read in tissue_positions file from spatial/
    positions = pd.read_csv(run.positions_file, header=None)
    positions.columns = ["barcode", "on_off", "row", "col", "xcor", "ycor"]

    # Match barcodes in adata/fragments_file
    positions["barcode"] = positions["barcode"] + "-1"

    # Merge fragments file with Anndata.obs
    adata.obs["barcode"] = adata.obs.index
    adata.obs = adata.obs.merge(positions, on="barcode", how="left")

    # Set run_id, codition
    adata.obs["run_id"] = run.run_id
    adata.obs["condition"] = run.condition

    # Ensure obs_names unique
    adata.obs_names = [
        run_id + "#" + bc for
        run_id, bc in zip(adata.obs["run_id"], adata.obs["barcode"])
    ]

    return adata


def add_tilematrix(
    adatas: List[anndata.AnnData],
    tile_size: int = 5000,
    n_features: int = 25000
) -> List[anndata.AnnData]:
    """Add tile matrix to AnnData.X, and select variable features.
    """

    snap.pp.add_tile_matrix(adatas, bin_size=tile_size)
    snap.pp.select_features(adatas, n_features=n_features)

    return adatas


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

    # Input AnnData must be backend for AnnDataSet, not in-memory.
    logging.info("Converting AnnData objects to backend...")
    adatas_be = [
        convert_tobackend(adata, f"{name}") for
        adata, name in zip(adatas, names)
    ]

    # Convert to AnnDataSet; this is what tutorial does, can we just combine
    # in-memory AnnData objects and bypass this read/write crap?
    logging.info("Creating AnnDataSet...")
    adataset = snap.AnnDataSet(
        adatas=[(name, adata) for name, adata in zip(names, adatas_be)],
        filename=f"{filename}.h5ad"
    )

    # Convert back to AnnData so we can add metadata :/
    combined_adata = adataset.to_adata()

    # AnnDataSet does inherit .obs; add manually :/
    combined_adata.obs = pd.concat([adata.obs for adata in adatas])

    # Ensure obs_names unique
    combined_adata.obs_names = [
        run_id + "#" + bc for
        run_id, bc in
        zip(combined_adata.obs["run_id"], combined_adata.obs["barcode"])
    ]

    # AnnDataSet does not inherit .obsm; add manually :/
    frags = vstack([adata.obsm["fragment_paired"] for adata in adatas])
    combined_adata.obsm["fragment_paired"] = frags

    return combined_adata


def convert_tobackend(
    adata: anndata.AnnData, filename: str
) -> anndata.AnnData:
    """Create a new backend AnnData object; necessary for creating AnnDataSet;
    saves each AnnData object to disk as .h5ad.
    """

    adata_backend = snap.AnnData(
        filename=f"{filename}_backend.h5ad",
        X=adata.X,
        obs=adata.obs,
        var=adata.var,
        uns=adata.uns,
        obsm=dict(adata.obsm)
    )
    adata_backend.obs_names = adata.obs_names

    return adata_backend


def filter_adatas(
    adatas: List[anndata.AnnData], min_tss: float = 2.0
) -> List[anndata.AnnData]:
    """Filter AnnData by on/off tissue tixels, TSS enrichment, max frag counts.
    """

    # Filter 'off tissue' tixels
    try:
        adatas = [adata[adata.obs["on_off"] == 1] for adata in adatas]
    except KeyError as e:
        logging.warn(f"Exception {e}: no positions data found in AnnData.obs")

    # Filter cells by tss, max_counts
    snap.pp.filter_cells(adatas, min_tsse=min_tss, max_counts=1e7)

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

    # As 'in_memory' so we can add metadata to .obs
    adatas = snap.pp.import_data(
        [run.fragments_file.local_path for run in runs],
        chrom_sizes=genome_ref,
        min_num_fragments=min_frags,
        sorted_by_barcode=False
    )

    # Add run_id, condition, spatial info to .obs, TSS enrichment
    adatas = [add_metadata(run, adata) for run, adata in zip(runs, adatas)]

    # Add addtional QCs
    snap.metrics.tsse(adatas, genome_ref)
    for adata in adatas:
        adata.obs["log10_frags"] = adata.obs["n_fragment"].apply(math.log10)

    return adatas
