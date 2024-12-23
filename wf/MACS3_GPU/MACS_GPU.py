import csv
import datetime as dt
from collections import defaultdict

import cupy as cp
import numpy as np
import pandas as pd
import polars as ps
import zstandard as zstd
from cupyx.scipy.ndimage import maximum_filter1d
from cupyx.scipy.special import gammaincc, pdtr, pdtrc


def Window_Max(vector, window):
    """
    Function to return the maximum value in a sliding window.
    inputs:
        - vector: cupy vector of the coverage along the chromosome
        - window: size of the sliding window
    output:
        - result: cupy vector containing the rolling max
    """
    result = maximum_filter1d(vector, size=window, mode="constant", cval=0)
    return result


def Window_Mean(vector, window):
    """
    Function to return the average in a sliding window.
    inputs:
        - vector: cupy vector of the coverage along the chromosome
        - window: size of the sliding window
    output:
        - result: cupy vector containing the rolling mean
    """
    rolling_max = (
        pd.Series(cp.asnumpy(vector)).rolling(window=window + 1, center=True).mean()
    )
    result = rolling_max.fillna(0).to_numpy()
    return cp.asarray(result)


def compute_q_values(p_values):
    """
    Function to compute the FDR corrected qvalues.
    inputs:
        - p_values: cupy vector of p-values at every point along the chromosome
    outputs:
        - q_values: cupy vector of q-values at every point along the chromosome
    """
    log_p = cp.log10(p_values)
    unique_p_values, inverse_indices, counts = cp.unique(
        log_p, return_inverse=True, return_counts=True
    )
    unique_p_values = cp.power(10, unique_p_values)
    n = cp.sum(counts)
    sorted_indices = cp.argsort(unique_p_values)
    sorted_unique_p_values = unique_p_values[sorted_indices]
    ranks = cp.concatenate((cp.asarray([1]), counts[sorted_indices][:-1]))
    ranks = cp.cumsum(ranks)
    q_values_unique = (sorted_unique_p_values * n) / ranks
    q_values_unique = q_values_unique[::-1]
    q_values_unique = np.minimum.accumulate(cp.asnumpy(q_values_unique))[::-1]
    q_values_unique = cp.asarray(q_values_unique)
    q_values_unique_unsorted = cp.empty_like(q_values_unique)
    q_values_unique_unsorted[sorted_indices] = q_values_unique
    original_order_q_values = q_values_unique_unsorted[inverse_indices]
    del q_values_unique
    return original_order_q_values


def Load_BedFile_df(zst_file):
    with open(zst_file, "rb") as f:
        dctx = zstd.ZstdDecompressor()
        with dctx.stream_reader(f) as reader:
            df = ps.read_csv(
                reader,
                separator="\t",
                new_columns=["chr", "start", "end", "strand", "qual"],
            )
    num_reads = len(df)
    result_df = df.group_by("chr").agg(
        [ps.col("start").implode(), ps.col("end").implode()]
    )

    # Convert to desired dictionary format
    result = {
        row["chr"]: {"start": row["start"], "end": row["end"]}
        for row in result_df.to_dicts()
    }
    return result, num_reads


def Pileup(starts, ends, chr_length, offset):
    """
    Function to compute the pileup at every point along the chromosome.
    Pileup is the total number of reads aligning to the chromosome.
    inputs:
        - p_values: cupy vector of p-values at every point along the chromosome
    outputs:
        - q_values: cupy vector of q-values at every point along the chromosome
    """
    start_poss = np.concatenate((starts - offset, ends - offset))
    end_poss = np.concatenate((starts + offset, ends + offset))
    cov_vec = cp.zeros(chr_length, dtype=np.int32)
    start_poss = np.clip(start_poss, 0, chr_length - 1)
    end_poss = np.clip(end_poss, 0, chr_length - 1)
    cp.add.at(cov_vec, start_poss, 1)
    cp.add.at(cov_vec, end_poss, -1)
    cov_vec = cp.cumsum(cov_vec)

    return cov_vec


def compute_peaks(q_vals, treat, ctrl, thresh=0.1):
    """
    Function to filter peaks based q-values. Return the indices where q-value < thresh
    inputs:
        - q_vals: cupy vector of q-values
        - thresh: Threshold parameter (default 0.1)
    outputs:
        - peaks: cupy vector of indices where q-values < thresh
    """
    peaks = cp.where((q_vals >= thresh))[0]
    return peaks


def merge_consecutive(arr, max_gap=30):
    """
    Function to merge overlapping and consecutive peak calls.
    Returns a matrix with two columns encoding the start and ends of peaks.
    inputs:
        - arr: cupy vector of coordinates returned by compute_peaks.
        - max_gap: parameter to merge peaks within "max_gap". Default [30]
    outputs:
        - M: a matrix with two columns encoding the start and ends of peaks
    """
    boundaries = cp.where(cp.diff(arr) >= max_gap)[0] + 1
    split_indices = cp.concatenate(
        (cp.asarray([0]), boundaries, cp.asarray([len(arr)]))
    )
    starts = arr[split_indices[:-1]]
    ends = arr[split_indices[1:] - 1]
    M = cp.stack([starts, ends], axis=1)
    del split_indices, boundaries, starts, ends
    return M


def Filter_Peaks(peaks, peak_amplitude=150):
    """
    Function to remove short peaks. Removes peaks smaller than 150 bases.
    inputs:
        - peaks: cupy matrix of coordinates returned by merge_consecutive.
        - peak_amplitude: integer representing the peak length
    outputs:
        - merged_filt: a matrix with two columns encoding the start and ends of peaks
    """
    peak_lengths = peaks[:, 1] - peaks[:, 0]
    peak_ind = cp.where(peak_lengths >= peak_amplitude)[0]
    merged_filt = peaks[peak_ind]
    return merged_filt


def calculate_peak_summits(peaks, signal):
    """
    Function to compute peak summits. Peak summit is the highest point along the peak.
    The highest point in the peak corresponds to the point with the lowest q-value
    inputs:
        - peaks: a nx2 numpy matrix encoding the peak coordinates
        - signal: a vector based on which peaks are called. Here the signal corresponds to the q-values
    outputs:
        - min_values: q-values corresponding to peaks
        - arg_min_indices: indices of the peak summits along the chromosomes.
    """
    starts = peaks[:, 0]
    ends = peaks[:, 1]
    max_values = np.zeros(len(starts))
    arg_max_indices = np.zeros(len(starts), dtype=np.int64)
    for i in range(0, len(starts)):
        am = starts[i] + np.argmin(signal[starts[i] : ends[i]])
        arg_max_indices[i] = am
    return arg_max_indices


def extract_values(indexes, vec):
    return vec[indexes]


def FDR(unique_p_values, unique_p_counts):
    """
    Compute the FDR corrected values using the Benjamini-Hochberg procedure.

    Parameters:
    - unique_p_values (cp.ndarray): Vector of p-values (assumed unique).
    - unique_p_counts (cp.ndarray): Vector of counts corresponding to the p-values.

    Returns:
    - pq_table (dict): Mapping of p-values to their q-values.
    """
    # Step 1: Sort unique p-values and counts
    unique_p_values = -1 * unique_p_values  # Negate p-values for sorting
    sorted_indices = cp.argsort(unique_p_values)[::-1]
    sorted_unique_p_values = unique_p_values[sorted_indices]
    sorted_counts = unique_p_counts[sorted_indices]
    cumulative_k = cp.cumsum(sorted_counts) - sorted_counts + 1
    sorted_unique_p_values = cp.where(
        sorted_unique_p_values == -cp.inf, -cp.inf, sorted_unique_p_values
    )
    total_counts = cp.sum(unique_p_counts)

    f = cp.log10(total_counts)
    q_values = cp.asnumpy(sorted_unique_p_values + (cp.log10(cumulative_k) - f))

    q_values_np = np.zeros(len(q_values))
    preq_q = float("inf")
    for i in range(len(q_values)):
        q = max(min(preq_q, q_values[i]), 0)
        preq_q = q
        q_values_np[i] = q

    p_values_original = -1 * sorted_unique_p_values  # Convert back to original p-values
    pq_table = dict(zip(cp.asnumpy(p_values_original), q_values_np))

    return pq_table


def replace_with_dict(array, mapping):
    """
    Function to replace every element in the array with value in the dictionary. Utility function to replace p_values with q_values.
    inputs:
        - array: a cupy vector
        - mapping: dictionary to replace the elements of the array
    outputs:
        - result: a vector whose values are replaced with the values from mapping.
    """
    keys = cp.array(list(mapping.keys()))
    values = cp.array(list(mapping.values()))

    sorted_indices = cp.argsort(keys)
    sorted_keys = keys[sorted_indices]
    sorted_values = values[sorted_indices]

    idx = cp.searchsorted(sorted_keys, array)
    valid = (idx < len(sorted_keys)) & (sorted_keys[idx] == array)

    result = cp.where(valid, sorted_values[idx], array)

    return result


def make_PQ_tables(d, num_reads, d_treat=150, d_ctrl=10000, genome_length=3088286401):
    """
    This is a method to compute the pq-table. Across all chromosomes, we estimate significance values. Upon computing the frequencies of the
    p-values, we compute calculate their ranks, and correct the p-values with Benjamini Hochberg correction procedure. The function returns a
    dictionary of p-values as keys, and their correspoding q-values as it values. Since the alignment is done across all chromosomes, in one
    go, we compute do this procedure first, to prevent biasing towards any one chromosome.

    inputs:
        - d : dictionary of coordinates of starts and ends for every chromosome.
              d = {"chr1":{"start":[], "end":[]}, "chr2":{"start":[], "end":[]}}
        - num_reads: number of reads
        - d: distance for extending alignments to compute peaks
        - d_ctrl: distance for extending alignments for building the null hypothesis
        - genome_length: mappable genome length

    outputs:
        - pq_table: A dictionary of q-value for every p-value
            {-14.533:-11:566, ...}
        Note the p and q-values are log10 transformed.
    """
    unique_p_values = cp.asarray([])
    unique_p_counts = cp.asarray([])
    for c in d:
        starts = np.array(d[c]["start"][0])
        ends = np.array(d[c]["end"][0])

        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()

        chrom_length = cp.max(ends)

        scale = d_treat / d_ctrl
        lambda_bg = 2 * d_treat * num_reads / genome_length

        pileup_treat = Pileup(starts, ends, chrom_length, int(d_treat / 2))
        pileup_ctrl = Pileup(starts, ends, chrom_length, int(d_ctrl / 2))
        pileup_ctrl = cp.maximum(pileup_ctrl * scale, lambda_bg)
        p_values = compute_poisson_cdfs(pileup_treat, pileup_ctrl)
        p_values, p_counts = cp.unique(p_values, return_counts=True)
        all_values = cp.concatenate((unique_p_values, p_values))
        all_counts = cp.concatenate((unique_p_counts, p_counts))
        merged_values, inverse_indices = cp.unique(all_values, return_inverse=True)
        merged_counts = cp.zeros_like(merged_values, dtype=all_counts.dtype)
        cp.add.at(merged_counts, inverse_indices, all_counts)

        unique_p_values = merged_values
        unique_p_counts = merged_counts

        del (
            pileup_ctrl,
            pileup_treat,
            p_values,
            merged_values,
            inverse_indices,
            merged_counts,
            all_values,
            all_counts,
        )

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        cp._default_memory_pool.free_all_blocks()
    pq_table = FDR(unique_p_values, unique_p_counts)
    return pq_table


def compute_poisson_cdfs(observations, lambdas):
    p_vals = pdtrc(observations, lambdas)
    return cp.log10(p_vals)


def call_peaks(
    starts,
    ends,
    pq_table,
    num_reads,
    q_thresh=0.1,
    d=150,
    d_ctrl=10000,
    genome_length=3088286401,
    max_gap=30,
    peak_amp=150,
):
    """
    Wrapper for the peak calling routines. This function computes the pileups from the data and the null hypothesis(pileup_ctrl).
    We compute the p_values, q values, identify and filter peaks.
    inputs:
        - starts: a list of start coordinates for the read alignments
        - ends: a list of end coordinates for the read alignments
        - num_reads: number of reads
        - q_thresh: threshold for q-values to identify peaks
        - d: distance for extending alignments to compute peaks
        - d_ctrl: distance for extending alignments for building the null hypothesis
        - genome_length: mappable genome length
        - max_gap: gaps for extending peaks
        - peak_amp: Amplitude for filtering peaks.
    outputs:
        - peaks_summit_q_vals: q-values of peak summits
        - peak_summits_args: peak summits
        - merged_peaks: Matrix of peak coordinates
    """
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()

    chrom_length = cp.max(ends)

    scale = d / d_ctrl
    lambda_bg = 2 * d * num_reads / genome_length
    q_thresh = -cp.log10(q_thresh)
    peak_amp = peak_amp - 1

    pileup_treat = Pileup(starts, ends, chrom_length, int(d / 2))
    pileup_ctrl = Pileup(starts, ends, chrom_length, int(d_ctrl / 2))
    pileup_ctrl = cp.maximum(pileup_ctrl * scale, lambda_bg)
    print(dt.datetime.now(), "Pileup")

    p_values = compute_poisson_cdfs(pileup_treat, pileup_ctrl)
    print(dt.datetime.now(), "Computed P Values")

    q_values = replace_with_dict(p_values, pq_table)

    p_values[p_values == -cp.inf] = -1000
    q_values[q_values == cp.inf] = 1000
    print(dt.datetime.now(), "Computed q Values")

    peaks = compute_peaks(q_values, pileup_treat, pileup_ctrl, q_thresh)
    print(dt.datetime.now(), "Calculated Peaks")

    if len(peaks) == 0:
        return pd.DataFrame()

    m = merge_consecutive(peaks, max_gap)
    filtered_peaks = Filter_Peaks(m, peak_amp)
    if len(filtered_peaks) == 0:
        return pd.DataFrame()

    print(
        dt.datetime.now(),
        "Merged Peaks",
        len(filtered_peaks),
        cp.max(filtered_peaks[:, 1] - filtered_peaks[:, 0]),
    )

    Peaks = cp.asnumpy(peaks)
    merged_peaks = cp.asnumpy(filtered_peaks)
    peak_summits_args = calculate_peak_summits(merged_peaks, cp.asnumpy(p_values))
    q_summit = cp.asnumpy(extract_values(peak_summits_args, q_values))
    p_summit = cp.asnumpy(extract_values(peak_summits_args, p_values))
    treat_summit = cp.asnumpy(extract_values(peak_summits_args, pileup_treat))
    ctrl_summit = cp.asnumpy(extract_values(peak_summits_args, pileup_ctrl))

    df_op = pd.DataFrame(
        data={
            "start": merged_peaks[:, 0],
            "end": merged_peaks[:, 1],
            "peak": peak_summits_args - merged_peaks[:, 0],
            "signal_value": (treat_summit + 1) / (ctrl_summit + 1),
            "p_value": -1 * p_summit,
            "q_value": q_summit,
            "pileup": treat_summit,
        }
    )
    df_op["name"] = "."
    df_op["score"] = np.minimum(np.array(df_op["q_value"] * 10, dtype=np.int64), 1000)
    df_op["strand"] = "."

    print(dt.datetime.now(), "Calulcated peaks summits")

    del q_values, p_values, pileup_treat, m, peaks, pileup_ctrl, filtered_peaks

    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()
    cp._default_memory_pool.free_all_blocks()

    return df_op


def MACS_GPU(
    d,
    num_reads,
    q_thresh=0.1,
    d_treat=150,
    d_ctrl=10000,
    genome_length=3088286401,
    max_gap=30,
    peak_amp=150,
):
    """
    This is a reimplementation of MACS3 used for calling peaks on ChIP and ATAC Seq data to work on high throughput
    sequencing runs using GPUs. Currently the functionalities are limited to narrow-peak calling.
    Inputs:
        - d: Dictionary of alignment coordinates.
             d = {
                    "chr1":{
                        start:[1000, 1500, ....],
                        end: [1200, 1770,...]
                        },
                    "chr2":{
                    ....
                    },
                    .....
                }
        - num_reads: total number of alignments
        - q_thresh: therehold for q-values
        - d_treat: analogous to the "extsize" parameter in MACS3
        - d_ctrl: window length for ctrl
        - genome_length: length of the genome
        - max_gap: gaps for extending peaks
        - peak_amp: Amplitude for filtering peaks.

    Output:
        - Peaks: A dataframe containing summaries about the peaks.
        ┌───────┬───────────┬───────────┬──────┬───┬──────────────┬───────────┬───────────┬──────┐
        │ chrom ┆ start     ┆ end       ┆ name ┆ … ┆ signal_value ┆ p_value   ┆ q_value   ┆ peak │
        │ ---   ┆ ---       ┆ ---       ┆ ---  ┆   ┆ ---          ┆ ---       ┆ ---       ┆ ---  │
        │ str   ┆ u64       ┆ u64       ┆ str  ┆   ┆ f64          ┆ f64       ┆ f64       ┆ u64  │
        ╞═══════╪═══════════╪═══════════╪══════╪═══╪══════════════╪═══════════╪═══════════╪══════╡
        │ chr1  ┆ 10192     ┆ 10576     ┆ .    ┆ … ┆ 4.428408     ┆ 6.571632  ┆ 3.621155  ┆ 54   │
        │ chr1  ┆ 183382    ┆ 183531    ┆ .    ┆ … ┆ 3.389831     ┆ 4.42147   ┆ 1.819646  ┆ 51   │
        │ chr1  ┆ 183676    ┆ 183825    ┆ .    ┆ … ┆ 3.389831     ┆ 4.42147   ┆ 1.819646  ┆ 34   │
        │ chr1  ┆ 267980    ┆ 268129    ┆ .    ┆ … ┆ 2.657045     ┆ 3.487542  ┆ 1.043066  ┆ 0    │
        │ chr1  ┆ 629867    ┆ 630036    ┆ .    ┆ … ┆ 10.326087    ┆ 18.87     ┆ 14.647516 ┆ 50   │
        │ …     ┆ …         ┆ …         ┆ …    ┆ … ┆ …            ┆ …         ┆ …         ┆ …    │
        │ chrX  ┆ 155843550 ┆ 155843699 ┆ .    ┆ … ┆ 3.542727     ┆ 4.981606  ┆ 2.255885  ┆ 48   │
        │ chrX  ┆ 155880726 ┆ 155880905 ┆ .    ┆ … ┆ 4.733728     ┆ 6.16004   ┆ 3.276946  ┆ 51   │
        │ chrX  ┆ 155881107 ┆ 155881564 ┆ .    ┆ … ┆ 7.831325     ┆ 12.405941 ┆ 8.816095  ┆ 133  │
        │ chrX  ┆ 155995895 ┆ 155996044 ┆ .    ┆ … ┆ 3.542727     ┆ 4.981606  ┆ 2.255885  ┆ 46   │
        │ chrY  ┆ 56837273  ┆ 56837422  ┆ .    ┆ … ┆ 4.615385     ┆ 6.105977  ┆ 3.229534  ┆ 107  │
        └───────┴───────────┴───────────┴──────┴───┴──────────────┴───────────┴───────────┴──────┘
    """
    start = dt.datetime.now()
    pq_table = make_PQ_tables(
        d, num_reads, genome_length=genome_length, d_treat=d_treat, d_ctrl=d_ctrl
    )
    end = dt.datetime.now()
    print(f"Calculated PQ Table....., {(end-start).total_seconds()/60.0}, minutes...")

    wf_start = dt.datetime.now()
    peaks = pd.DataFrame()

    for chr in d:
        start = dt.datetime.now()
        starts = np.array(d[chr]["start"][0])
        ends = np.array(d[chr]["end"][0])
        print(start, chr)
        df_chr = call_peaks(
            starts,
            ends,
            pq_table,
            num_reads,
            max_gap=max_gap,
            q_thresh=q_thresh,
            peak_amp=peak_amp,
            genome_length=genome_length,
            d=d_treat,
            d_ctrl=d_ctrl,
        )
        df_chr["chrom"] = chr
        peaks = pd.concat([peaks, df_chr], axis=0)
        end = dt.datetime.now()
        print(end, (end - start).total_seconds(), "Done\n")
    peaks = peaks[
        [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signal_value",
            "p_value",
            "q_value",
            "peak",
        ]
    ]
    wf_end = dt.datetime.now()
    print(f"Total time elapsed..., {(wf_end-wf_start).total_seconds()/60.0}, minutes")

    return peaks
