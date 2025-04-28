import logging
from typing import List, Optional

import anndata
import numpy as np
import pandas as pd
import pychromvar as pc
import scanpy as sc
import snapatac2 as snap

from pybedtools import BedTool
from pyjaspar import jaspardb

from wf.utils import ref_dict

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def annotate_peaks(
    peaks_df: pd.DataFrame,
    features: List[pd.DataFrame]
) -> pd.DataFrame:

    peaks_df = reformat_peak_df(peaks_df)
    peaks = make_peak_bed(peaks_df)
    genes, exons, promoters = [BedTool.from_dataframe(feature)
                               for feature in features]

    nearest_genes = peaks.closest(genes.sort(), d=True, t="first")

    nearest_df = nearest_genes.to_dataframe(
        names=["peak_chrom", "peak_start", "peak_end", "peak_id",
               "gene_chrom", "gene_start", "gene_end", "gene_name",
               "gene_score", "gene_strand", "distance"]
    )
    nearest_df.index = nearest_df["peak_id"]

    # Add distance and nearest gene info to original dataframe
    peaks_df["distToGeneStart"] = nearest_df["distance"]
    peaks_df["nearestGene"] = nearest_df["gene_name"]

    # Find overlaps with different features
    promoter_overlaps = peaks.intersect(promoters, u=True).to_dataframe()
    gene_overlaps = peaks.intersect(genes, u=True).to_dataframe()
    exon_overlaps = peaks.intersect(exons, u=True).to_dataframe()

    # Initialize all peaks as Distal
    peaks_df["peakType"] = "Distal"

    # Create boolean masks for overlaps
    peak_ids = peaks_df.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
    )

    # Apply ArchR's logic for peak type classification
    peaks_in_genes = peak_ids.isin(gene_overlaps.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
    ))
    peaks_in_exons = peak_ids.isin(exon_overlaps.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
    ))
    peaks_in_promoters = peak_ids.isin(promoter_overlaps.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
    ))

    # Set peak types following ArchR's logic
    peaks_df.loc[peaks_in_genes & peaks_in_exons, "peakType"] = "Exonic"
    peaks_df.loc[peaks_in_genes & ~peaks_in_exons, "peakType"] = "Intronic"
    peaks_df.loc[peaks_in_promoters, "peakType"] = "Promoter"

    peaks_df = peaks_df.drop("id", axis=1)

    return peaks_df


def clean_adata(adata: anndata.AnnData) -> anndata.AnnData:

    obs = ["barcode", "n_genes_by_counts", "log1p_n_genes_by_counts", "total_counts", "log1p_total_counts", "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes", "pct_counts_in_top_200_genes", "pct_counts_in_top_500_genes", "total_counts_mt", "log1p_total_counts_mt", "pct_counts_mt"]
    var = ["mt", "n_counts", "n_cells", "n_cells_by_counts", "mean_counts", "log1p_mean_counts", "pct_dropout_by_counts", "total_counts", "log1p_total_counts"]

    rm_obs = [o for o in obs if o in adata.obs.keys()]
    rm_var = [v for v in var if v in adata.var.keys()]

    adata.obs.drop(rm_obs, axis=1, inplace=True)
    adata.var.drop(rm_var, axis=1, inplace=True)

    del adata.varm
    del adata.layers

    if adata.raw:
        del adata.raw

    for uns in ["pca", "log1p"]:
        if uns in adata.uns.keys():
            del adata.uns[uns]

    rm_obsm = [obsm for obsm in adata.obsm.keys()
               if obsm not in ["spatial", "X_umap"]]

    for obsm in rm_obsm:
        del adata.obsm[obsm]

    try:
        adata.X = adata.X.astype(np.float16)
    except:
        logging.warning("Cannot convert .X to float")

    return adata


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

    # Copy data for neighborhood enrichment plots
        try:
            adata_ge.obsp["spatial_connectivities"] = adata.obsp["spatial_connectivities"]
            adata_ge.uns["cluster_nhood_enrichment"] = adata.uns["cluster_nhood_enrichment"]
        except Exception as e:
            logging.warning(
                f"Exception {e}: Could not copy neighborhood enrichment data."
            )

    # Remove genes with no cells, counts; per sc, one metric per call...
    logging.info("Removing mitochondtrial genes, filtering...")

    print(f"Pre-filtering shape: {adata_ge.shape}")
    adata_ge.var["mt"] = adata_ge.var_names.str.startswith(("MT-", "Mt-", "mt-"))
    adata_ge = adata_ge[:, ~adata_ge.var["mt"]].copy()
    print(f"post-filtering shape: {adata_ge.shape}")

    sc.pp.filter_genes(adata_ge, min_counts=min_counts)
    sc.pp.filter_genes(adata_ge, min_cells=min_cells)

    logging.info("Normalizing matrix and computing log...")
    sc.pp.normalize_total(adata_ge)
    sc.pp.log1p(adata_ge)

    logging.info("Regressing out n_fragments...")
    regressed_ge = sc.pp.regress_out(adata_ge, keys=["n_fragment"], n_jobs=-1)

    logging.info("Batch correction with MAGIC...")
    sc.external.pp.magic(adata_ge, solver="approximate", n_jobs=-1)

    return adata_ge


def make_motifmatrix(adata: anndata.AnnData, n_jobs: int = -1) -> anndata.AnnData:
    """Return a AnnData object with X as a motif deviation matrix."""
    if adata.X.dtype != "float64":
        adata.X = adata.X.astype(np.float64)

    return pc.compute_deviations(adata, n_jobs=n_jobs)


def make_peak_bed(peaks_df: pd.DataFrame, sort: bool = True) -> BedTool:

    columns = ["chrom", "start", "end", "id"]

    try:
        peaks_df[columns]
    except KeyError as e:
        logging.warning(f"Error {e}: Expected column names missing.")

    bed = BedTool.from_dataframe(peaks_df[columns])

    return bed.sort() if sort else bed


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


def reformat_peak_df(peaks_df: pd.DataFrame) -> pd.DataFrame:
    """ Expects a DataFrame from sc.rank_genes_groups_df with column 'names,
    containing coordinates in the format 'chr1:631135-631636'; returns a
    DataFrame reformatted for for compatibility with BedTools.
    """

    try:
        peaks_df[["names", "group"]]
    except KeyError as e:
        logging.warning(f"Error {e}: Expected column names missing.")

    peaks_df[["chrom", "range"]] = peaks_df["names"].str.split(":", expand=True)
    peaks_df[["start", "end"]] = peaks_df["range"].str.split("-", expand=True)

    # Make a unique identifier for mapping
    peaks_df["group"] = peaks_df["group"].astype(str)
    peaks_df["id"] = peaks_df["group"].str.cat(peaks_df["names"], sep=":")
    peaks_df.index = peaks_df["id"]

    return peaks_df
