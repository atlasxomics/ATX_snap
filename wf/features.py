import anndata
import logging
import numpy as np
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap

from pyjaspar import jaspardb
from typing import List, Optional


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


def get_motifs(
    adata: anndata.AnnData, fasta_path: str, release: str = "JASPAR2024"
) -> anndata.AnnData:
    """With the python implementation of chromVAR, pychromvar, map motifs to
    peaks; return the AnnData with peak sequences (uns.['peak_seq']), sequence
    gc (.var['gc_bias']), background peaks (.varm['bg_peaks']), and motifs
    (varm['motif_match']).
    """

    jdb_obj = jaspardb(release=release)
    motifs = jdb_obj.fetch_motifs(
        collection="CORE",
        tax_group=["vertebrates"]
    )

    pc.add_peak_seq(adata, genome_file=fasta_path, delimiter=":|-")
    pc.add_gc_bias(adata)
    pc.get_bg_peaks(adata)

    pc.match_motif(adata, motifs=motifs)

    return adata


def make_peakmatrix(
    adata: anndata.AnnData,
    genome: str,
    key: str,
    log_norm: bool = True,
    obs: Optional[List[str]] = ["n_fragment", "tsse", "log10_frags"],
    obsm: Optional[List[str]] = ["spatial", "X_umap"]
) -> anndata.AnnData:
    """Given an AnnData object with macs2 peak calls stored in .uns[key],
    returns a new AnnData object with X a peak count matrix.
    """

    peaks = adata.uns[key]

    if not isinstance(peaks, dict):  # Convert to dict for merge_peaks()
        peaks = {"0": adata.uns[key]}

    # Can't use a dict because of flyte
    genome_ref = snap.genome.mm10 if genome == "mm10" else snap.genome.hg38
    merged_peaks = snap.tl.merge_peaks(peaks, genome_ref)

    adata_p = snap.pp.make_peak_matrix(adata, use_rep=merged_peaks["Peaks"])

    # Copy over cell data
    for ob in obs:
        adata_p.obs[ob] = adata.obs[ob]
    for ob in obsm:
        adata_p.obsm[ob] = adata.obsm[ob]
        if type(adata_p.obsm[ob]) is not np.ndarray:
            adata_p.obsm[ob] = adata_p.obsm[ob].to_numpy()

    if log_norm:
        sc.pp.log1p(adata_p)

    return adata_p


def make_geneadata(
    adata: anndata.AnnData,
    genome: str,
    min_counts: int = 1,
    min_cells: int = 1,
) -> anndata.AnnData:
    """Create an AnnData object where X is a Gene Expression Matrix; .obs is
    inherited from input AnnData; filter genes with low cells, counts.
    Parameters recapitulate ArchR defaults.
    """

    # Can't use a dict because of flyte
    genome_ref = snap.genome.mm10 if genome == "mm10" else snap.genome.hg38

    # New AnnData, parameters to match ArchR
    logging.info("Creating gene matrix...")
    adata_ge = snap.pp.make_gene_matrix(
        adata,
        gene_anno=genome_ref,
        upstream=5000,  # ArchR default
        downstream=0,  # ArchR default
        include_gene_body=True  # Use genebody, not TSS, per ArchR
    )

    # Copy adata .obsm
    for obsm in ["X_umap", "X_spectral_harmony", "spatial"]:
        try:
            adata_ge.obsm[obsm] = adata.obsm[obsm]
        except Exception as e:
            logging.warning(
                f"Exception {e}: no annotation {obsm} found for observations."
            )

    # Remove genes with no cells, counts; per sc, one metric per call...
    logging.info("Removing mitochondtrial genes, filtering...")

    print(f"Pre-filtering shape: {adata_ge.shape}")
    adata_ge.var["mt"] = adata_ge.var_names.str.startswith("MT-")
    adata_ge = adata_ge[:, ~adata_ge.var["mt"]].copy()
    print(f"post-filtering shape: {adata_ge.shape}")

    sc.pp.filter_genes(adata_ge, min_counts=min_counts)
    sc.pp.filter_genes(adata_ge, min_cells=min_cells)

    logging.info("Normalizing matrix and computing log...")
    sc.pp.normalize_total(adata_ge)
    sc.pp.log1p(adata_ge)

    logging.info("Batch correction with MAGIC...")
    sc.external.pp.magic(adata_ge, solver="approximate")

    sc.pp.calculate_qc_metrics(
        adata_ge, qc_vars="mt", inplace=True, log1p=True
    )

    return adata_ge


def make_motifmatrix(
    adata: anndata.AnnData, n_jobs: int = -1
) -> anndata.AnnData:
    """Return a AnnData object with X as a motif deviation matrix.
    """

    if adata.X.dtype != 'float64':
        adata.X = adata.X.astype(np.float64)

    return pc.compute_deviations(adata, n_jobs=n_jobs)


def rank_features(
    adata: anndata.AnnData,
    groups: List[str],
    feature_type: str,
    save: Optional[str],
    use_raw: bool = False,
    pval_cutoff: float = 0.05
):
    """For each metadata cell grouping provided, add gene ranking information;
    if 'save' is a string, a csv of the rank data is saved to a directory
    specified by the string.
    """

    for group in groups:

        logging.info(f"Finding marker genes for {group}s...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=group,
            method="t-test",
            key_added=f"{group}_{feature_type}",
            use_raw=use_raw
        )

        # Write marker genes to csv
        sc.get.rank_genes_groups_df(
            adata,
            group=None,
            key=f"{group}_{feature_type}",
            pval_cutoff=pval_cutoff,
        ).to_csv(f"{save}/marker_{feature_type}_per_{group}.csv", index=False)
