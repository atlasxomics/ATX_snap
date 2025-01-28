import logging
from typing import List, Optional

import anndata
import numpy as np
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap
from pyjaspar import jaspardb

from wf.utils import ref_dict

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
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
    motifs = jdb_obj.fetch_motifs(collection="CORE", tax_group=["vertebrates"])

    pc.add_peak_seq(adata, genome_file=fasta_path, delimiter=":|-")
    pc.add_gc_bias(adata)
    pc.get_bg_peaks(adata)

    pc.match_motif(adata, motifs=motifs)

    return adata


def make_peakmatrix(
    adata: anndata.AnnData, genome: str, key: str, log_norm: bool = True
) -> anndata.AnnData:
    """Given an AnnData object with macs2 peak calls stored in .uns[key],
    returns a new AnnData object with X a peak count matrix.
    """
    peaks = adata.uns[key]

    if not isinstance(peaks, dict):  # Convert to dict for merge_peaks()
        peaks = {"0": adata.uns[key]}

    # Can't use a dict because of flyte
    merged_peaks = snap.tl.merge_peaks(peaks, ref_dict[genome][0])

    adata_p = snap.pp.make_peak_matrix(adata, use_rep=merged_peaks["Peaks"])

    # Copy over cell data
    adata_p.obs = adata.obs
    adata_p.obsm = adata.obsm

    if log_norm:
        sc.pp.log1p(adata_p)

    return adata_p


def make_geneadata(
    adata: anndata.AnnData,
    genome: str,
    min_counts: int = 1,
    min_cells: int = 1,
) -> anndata.AnnData:
    """Create an AnnData object where X is a Gene Expression Matrix (GEM); .obs
    is inherited from input AnnData; filter genes with low cells, counts.
    Parameters recapitulate ArchR defaults.  snap.pp.make_gene_matrix() first
    creates a GEM of putative raw counts.  This matrix is normalized, log
    transformed, and imputed with MAGIC.  In the returned AnnData object, .X is
    the imputated matrix, .raw.X is the log transformed, normalized matrix; raw
    counts are stored in .layers["raw_counts"].
    """

    # New AnnData, parameters to match ArchR
    logging.info("Creating gene matrix...")
    adata_ge = snap.pp.make_gene_matrix(
        adata,
        gene_anno=ref_dict[genome][1],
        upstream=5000,  # ArchR default
        downstream=0,  # ArchR default
        include_gene_body=True,  # Use genebody, not TSS, per ArchR
    )

    # Save raw counts in layers
    adata_ge.layers["raw_counts"] = adata_ge.X.copy()

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

    sc.pp.calculate_qc_metrics(adata_ge, qc_vars="mt", inplace=True, log1p=True)
    print("Done.... MAGIC+QC")
    return adata_ge


def make_motifmatrix(adata: anndata.AnnData, n_jobs: int = -1) -> anndata.AnnData:
    """Return a AnnData object with X as a motif deviation matrix."""
    if adata.X.dtype != "float64":
        adata.X = adata.X.astype(np.float64)

    return pc.compute_deviations(adata, n_jobs=n_jobs)


def rank_features(
    adata: anndata.AnnData,
    groups: List[str],
    feature_type: str,
    save: Optional[str],
    use_raw: bool = False,
    pval_cutoff: float = 0.05,
    logfoldchange_cutoff: float = 0.1,
):
    """For each metadata cell grouping provided, add gene ranking information;
    if 'save' is a string, a csv of the rank data is saved to a directory
    specified by the string.
    """
    import pandas as pd

    for group in groups:
        logging.info(f"Finding marker genes for {group}s...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=group,
            method="t-test",
            key_added=f"{group}_{feature_type}",
            use_raw=use_raw,
        )
        df = sc.get.rank_genes_groups_df(
            adata, group=None, key=f"{group}_{feature_type}"
        )

        # Write ranked features to csv
        df.to_csv(f"{save}/ranked_{feature_type}_per_{group}.csv", index=False)

        # Filter and summarize marker features
        df = df[df["pvals_adj"] <= pval_cutoff]
        df = df[abs(df["logfoldchanges"]) > logfoldchange_cutoff]

        counts = df.groupby("group").agg("count")["names"]
        counts["total"] = counts.sum()
        counts = pd.DataFrame(counts).reset_index()

        counts.rename(
            columns={
                "group": f"{group}",
                "names": f"number differential {feature_type}",
            },
            inplace=True,
        )
        counts.to_csv(
            f"{save}/differential_{feature_type}_per_{group}_counts.csv", index=False
        )
