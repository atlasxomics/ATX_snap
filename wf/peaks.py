from os import makedirs
from shutil import rmtree

import anndata
import gmacs as gm
import numpy as np
import polars as pl
import snapatac2._snapatac2 as _snapatac2


def call_peaks_macs3_gpu(
    adata: anndata.AnnData,
    groupby_key: str,
    d_treat=150,
    d_ctrl: int = 10000,
    max_gap: int = 30,
    peak_amp: int = 150,
    q_thresh: float = 0.1,
) -> anndata.AnnData:
    """
    Wrapper around gmacs function, peak caller routine implemented for identifying peaks,
    as a potential faster alternative for MACS3.
    """
    genome_length = adata.uns["reference_sequences"]["reference_seq_length"].sum()
    groupby = list(adata.obs[groupby_key])
    tmpdir = "tmp"
    makedirs(tmpdir, exist_ok=True)
    fragments = _snapatac2.export_tags(adata, tmpdir, group_by=groupby)

    peaks_dict = {}
    for f in fragments:
        print("cluster....", f)
        for cluster in fragments[f]:
            d, num_reads = gm.load_bedfile(cluster)
            peaks = gm.gmacs(
                intervals_by_chrom=d,
                num_reads=num_reads,
                d_treat=d_treat,
                d_ctrl=d_ctrl,
                genome_length=genome_length,
                max_gap=max_gap,
                peak_amp=peak_amp,
                q_thresh=q_thresh,
            )
            peaks["score"] = peaks["score"].astype(np.uint16)
            peaks = pl.from_pandas(peaks)
            peaks = peaks.with_columns(
                [
                    pl.col("chrom").cast(pl.String),
                    pl.col("start").cast(pl.UInt64),
                    pl.col("end").cast(pl.UInt64),
                    pl.col("name").cast(pl.String),
                    pl.col("score").cast(pl.UInt16),
                    pl.col("strand").cast(pl.String),
                    pl.col("signal_value").cast(pl.Float64),
                    pl.col("p_value").cast(pl.Float64),
                    pl.col("q_value").cast(pl.Float64),
                    pl.col("peak").cast(pl.UInt64),
                ]
            )
            peaks_dict[f] = peaks
            print(f, num_reads, len(peaks))
    rmtree(tmpdir)
    if adata.isbacked:
        adata.uns[f"{groupby_key}_peaks"] = peaks_dict
    else:
        adata.uns[f"{groupby_key}_peaks"] = {
            k: v.to_pandas() for k, v in peaks_dict.items()
        }
    return adata
