import logging
import gc
import glob

from pathlib import Path
from typing import Dict, List, Optional

import anndata
from anndata.experimental import concat_on_disk
import numpy as np
import scipy.sparse as sp

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def clean_adata(adata: anndata.AnnData) -> anndata.AnnData:
    obs = [
        "barcode", "n_genes_by_counts", "log1p_n_genes_by_counts", "total_counts",
        "log1p_total_counts", "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
        "pct_counts_in_top_200_genes", "pct_counts_in_top_500_genes", "total_counts_mt",
        "log1p_total_counts_mt", "pct_counts_mt"
    ]
    var = [
        "mt", "n_counts", "n_cells", "n_cells_by_counts", "mean_counts", "log1p_mean_counts",
        "pct_dropout_by_counts", "total_counts", "log1p_total_counts", "means",
        "dispersions", "dispersions_norm"
    ]

    rm_obs = [o for o in obs if o in adata.obs.keys()]
    rm_var = [v for v in var if v in adata.var.keys()]

    adata.obs.drop(rm_obs, axis=1, inplace=True)
    adata.var.drop(rm_var, axis=1, inplace=True)

    adata.varm.clear()
    adata.layers.clear()

    if adata.raw:
        adata.raw = None

    for uns in ["pca", "log1p"]:
        adata.uns.pop(uns, None)

    rm_obsm = [obsm for obsm in adata.obsm.keys()
               if obsm not in ["spatial", "X_umap"]]

    for obsm in rm_obsm:
        del adata.obsm[obsm]

    # Sparse to dense conversion
    if sp.issparse(adata.X):
        try:
            adata.X = adata.X.toarray()
        except Exception as e:
            logging.warning(f"Could not convert sparse .X to dense: {e}")

    try:
        adata.X = adata.X.astype(np.float16)
    except Exception as e:
        logging.warning(f"Cannot convert .X to float16: {e}")

    return adata


def _cast_X_dtype(
    adata: anndata.AnnData,
    dtype: str,
    label: str,
) -> None:
    try:
        adata.X = adata.X.astype(dtype)
    except Exception as e:
        logging.warning(f"Cannot convert {label} .X to {dtype}: {e}")


def _sanitize_dataframe_for_h5ad(df) -> None:
    """Ensure object columns can be written as H5AD string arrays."""
    import pandas as pd

    obj_cols = df.select_dtypes(include=["object"]).columns
    if len(obj_cols) == 0:
        return

    for col in obj_cols:
        series = df[col]
        non_null = series.dropna()
        if non_null.empty:
            df[col] = series.astype(str)
            continue

        # If all non-null values are numeric, keep as numeric.
        if non_null.map(
            lambda x: isinstance(x, (int, float, np.integer, np.floating))
        ).all():
            df[col] = pd.to_numeric(series, errors="coerce")
            continue

        # If all non-null values are already strings, keep as strings.
        if non_null.map(lambda x: isinstance(x, str)).all():
            df[col] = series.astype(str)
            continue

        # Fallback: coerce to string to avoid h5py conversion errors.
        df[col] = series.astype(str)


def _sanitize_uns_for_h5ad(uns: Dict) -> None:
    """Recursively sanitize .uns contents for H5AD writing."""
    import pandas as pd

    for _, value in list(uns.items()):
        if isinstance(value, pd.DataFrame):
            _sanitize_dataframe_for_h5ad(value)
        elif isinstance(value, dict):
            _sanitize_uns_for_h5ad(value)
        elif isinstance(value, (list, tuple)):
            for item in value:
                if isinstance(item, pd.DataFrame):
                    _sanitize_dataframe_for_h5ad(item)
                elif isinstance(item, dict):
                    _sanitize_uns_for_h5ad(item)


def clean_index_columns(*adatas: anndata.AnnData) -> None:
    """Remove _index columns from AnnData objects if they exist."""
    for adata in adatas:
        if hasattr(adata, 'raw') and adata.raw is not None:
            if "_index" in adata.raw.var.columns:
                adata.raw.var.drop(columns=['_index'], inplace=True)


def _ensure_anndata_root_encoding(path: Path) -> None:
    """Add missing AnnData encoding metadata for older H5AD writers."""
    import h5py

    with h5py.File(path, "r") as handle:
        encoding_type = handle.attrs.get("encoding-type")
        encoding_version = handle.attrs.get("encoding-version")

    if isinstance(encoding_type, bytes):
        encoding_type = encoding_type.decode()
    if isinstance(encoding_version, bytes):
        encoding_version = encoding_version.decode()

    if encoding_type != "anndata":
        try:
            backed = anndata.read_h5ad(path, backed="r")
            backed.file.close()
        except Exception as e:
            raise ValueError(
                f"{path} is not a readable AnnData h5ad file. "
                f"Root encoding-type={encoding_type!r}."
            ) from e

    with h5py.File(path, "r+") as handle:
        if encoding_type != "anndata":
            logging.warning(
                "Adding missing AnnData root encoding attrs to h5ad written by an "
                f"older converter: {path}"
            )
            handle.attrs["encoding-type"] = "anndata"
            handle.attrs["encoding-version"] = encoding_version or "0.1.0"

        for group_name in ["layers", "obsm", "varm", "obsp", "varp", "uns"]:
            if group_name not in handle:
                logging.warning(
                    f"Adding missing empty AnnData group '{group_name}' to {path}"
                )
                group = handle.create_group(group_name)
            else:
                group = handle[group_name]

            if group.attrs.get("encoding-type") is None:
                group.attrs["encoding-type"] = "dict"
                group.attrs["encoding-version"] = "0.1.0"


def load_and_combine_data(
    suffix: str,
    input_dir: Path = Path("."),
    temp_dir: Optional[Path] = None,
) -> anndata.AnnData:
    """Load and combine AnnData objects."""
    logging.info("Combining AnnData objects on disk...")
    input_dir = Path(input_dir)
    temp_dir = Path(temp_dir) if temp_dir is not None else input_dir
    temp_dir.mkdir(parents=True, exist_ok=True)
    combined_path = temp_dir / f"combined_{suffix}.tmp.h5ad"
    if combined_path.exists():
        combined_path.unlink()

    files = sorted(
        str(path)
        for path in input_dir.glob(f"*{suffix}.h5ad")
        if path.resolve() != combined_path.resolve()
    )
    if len(files) == 0:
        raise FileNotFoundError(
            f"No AnnData files found matching '*{suffix}.h5ad' in {input_dir}."
        )

    logging.info(f"Found {len(files)} AnnData file(s) to combine.")
    for file in files:
        _ensure_anndata_root_encoding(Path(file))

    # Match anndata/scanpy concat defaults while avoiding a list of full
    # in-memory AnnData objects plus a second combined copy.
    concat_on_disk(
        files,
        combined_path,
        axis=0,
        join="inner",
        merge=None,
        uns_merge=None,
        label=None,
        keys=None,
        index_unique=None,
        pairwise=False,
    )

    adata = anndata.read_h5ad(combined_path)
    combined_path.unlink(missing_ok=True)
    gc.collect()

    # Clean up index columns if they exist
    clean_index_columns(adata)

    return adata


def _is_valid_sample_name(value: object) -> bool:
    if value is None:
        return False

    sample_name = str(value).strip()
    return sample_name != "" and sample_name.lower() != "none"


def _get_sample_name_map(adata: anndata.AnnData) -> Dict[str, str]:
    if "sample" not in adata.obs.columns or "sample_name" not in adata.obs.columns:
        return {}

    sample_df = adata.obs[["sample", "sample_name"]].copy()
    if sample_df.empty:
        return {}

    sample_name_map = {}
    sample_df["sample"] = sample_df["sample"].astype(str)

    for sample, group in sample_df.groupby("sample", sort=False):
        valid_names = [
            str(name).strip()
            for name in group["sample_name"]
            if _is_valid_sample_name(name)
        ]

        if len(valid_names) == 0:
            continue

        unique_names = list(dict.fromkeys(valid_names))
        if len(unique_names) > 1:
            logging.warning(
                f"Multiple sample_name values found for sample '{sample}'; "
                f"using '{unique_names[0]}'."
            )

        sample_name_map[sample] = unique_names[0]

    return sample_name_map


def _remap_sample_labels_in_df(df, sample_name_map: Dict[str, str]):
    remapped = df.copy()

    remapped.columns = [
        sample_name_map.get(str(col), col) for col in remapped.columns
    ]
    remapped.index = [
        sample_name_map.get(str(idx), idx) for idx in remapped.index
    ]

    obj_cols = remapped.select_dtypes(include=["object"]).columns
    for col in obj_cols:
        remapped[col] = remapped[col].map(
            lambda x: sample_name_map.get(x, x) if isinstance(x, str) else x
        )

    return remapped


def _sample_name_key(key: str) -> str:
    if "sample_name" in key:
        return key

    parts = key.split("_")
    changed = False
    for i, part in enumerate(parts):
        if part == "sample":
            parts[i] = "sample_name"
            changed = True

    if changed:
        return "_".join(parts)

    return key


def _add_sample_name_results(adata: anndata.AnnData) -> None:
    import pandas as pd

    sample_name_map = _get_sample_name_map(adata)
    if len(sample_name_map) == 0:
        return

    for key, value in list(adata.uns.items()):
        if "sample" not in key or "sample_name" in key:
            continue

        mapped_key = _sample_name_key(key)
        if mapped_key == key or mapped_key in adata.uns:
            continue

        if isinstance(value, pd.DataFrame):
            adata.uns[mapped_key] = _remap_sample_labels_in_df(
                value, sample_name_map
            )
        else:
            adata.uns[mapped_key] = value


def load_analysis_results(
    adata: anndata.AnnData,
    type: str,
    groups: List[str],
    input_dir: Path = Path("."),
) -> None:
    """Load various analysis results into AnnData objects."""
    # Load differential analysis results
    if type == "gene":
        load_ranked_genes(adata, input_dir=input_dir)
    elif type == "motif":
        load_enriched_motifs(adata, input_dir=input_dir)

    # Load heatmap data
    load_heatmaps(adata, type, input_dir=input_dir)

    # Load volcano plots if condition analysis was performed
    if "condition" in groups:
        load_volcano_plots(adata, type, input_dir=input_dir)

    _add_sample_name_results(adata)


def load_csv_files_to_uns(
    pattern: str,
    target_uns: Dict,
    dtype_spec: Optional[Dict] = None,
    index_col: Optional[int] = None,
    name_transform: Optional[callable] = None,
    input_dir: Path = Path("."),
) -> None:
    import pandas as pd
    """Generic function to load CSV files matching a pattern into uns
    dictionary."""
    try:
        files = glob.glob(str(Path(input_dir) / pattern))
        for file in files:
            name = Path(file).stem
            if name_transform:
                name = name_transform(name, file)
            df = pd.read_csv(file, dtype=dtype_spec, index_col=index_col)
            _sanitize_dataframe_for_h5ad(df)
            target_uns[name] = df
    except Exception as e:
        logging.warning(f"Error loading files matching {pattern}: {e}")


def load_enriched_motifs(
    adata_motif: anndata.AnnData,
    input_dir: Path = Path("."),
) -> None:
    """Load enriched motifs data."""
    logging.info("Adding enriched motifs...")
    load_csv_files_to_uns(
        "enrichedMotifs_*.csv",
        adata_motif.uns,
        dtype_spec={"group_name": str},
        input_dir=input_dir,
    )


def load_heatmaps(
    adata: anndata.AnnData,
    type: str,
    input_dir: Path = Path("."),
) -> None:
    """Load heatmap data for both genes and motifs."""
    logging.info("Adding heatmaps...")

    if type == "gene":
        # Gene heatmaps
        load_csv_files_to_uns(
            "genes_per_*_hm.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            index_col=0,
            input_dir=input_dir,
        )

    elif type == "motif":
        # Motif heatmaps
        load_csv_files_to_uns(
            "motif_per_*_hm.csv",
            adata.uns,
            index_col=0,
            input_dir=input_dir,
        )


def load_ranked_genes(
    adata_gene: anndata.AnnData,
    input_dir: Path = Path("."),
) -> None:
    """Load ranked genes data."""
    logging.info("Adding ranked genes...")
    load_csv_files_to_uns(
        "ranked_genes_*.csv",
        adata_gene.uns,
        dtype_spec={"group_name": str},
        input_dir=input_dir,
    )


def load_volcano_plots(
    adata: anndata.AnnData,
    type: str,
    input_dir: Path = Path("."),
) -> None:
    """Load volcano plot data for genes and motifs."""
    if type == "gene":
        logging.info("Adding gene volcanos...")
        load_csv_files_to_uns(
            "volcanoMarkers_genes_*.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            name_transform=lambda name, file: volcano_name_transform(name, file, "genes"),
            input_dir=input_dir,
        )

    elif type == "motif":
        logging.info("Adding motif volcanos...")
        load_csv_files_to_uns(
            "volcanoMarkers_motifs_*.csv",
            adata.uns,
            dtype_spec={"cluster": str},
            name_transform=lambda name,
            file: volcano_name_transform(name, file, "motifs"),
            input_dir=input_dir,
        )


def save_anndata_objects(
    adata: anndata.AnnData,
    suffix: str,
    base_dir: Path,
    full_x_dtype: Optional[str] = None,
) -> None:
    """Save full and reduced AnnData objects."""
    logging.info("Saving full adata...")
    if full_x_dtype is not None:
        _cast_X_dtype(adata, full_x_dtype, "full adata")
    _sanitize_uns_for_h5ad(adata.uns)
    # Save full objects
    adata.write(f"{base_dir}/combined{suffix}.h5ad")

    # Create and save reduced objects
    logging.info("Making reduced adata...")
    sm_adata = clean_adata(adata)

    logging.info("Saving reduced adata...")
    _sanitize_uns_for_h5ad(sm_adata.uns)
    sm_adata.write(f"{base_dir}/combined_sm{suffix}.h5ad")


def transfer_auxiliary_data(
    adata: anndata.AnnData,
    data_paths: Dict[str, str],
    groups: List[str]
) -> None:
    """Transfer observation, UMAP, and spatial data to AnnData objects."""
    logging.info("Transferring auxiliary data...")

    # Transfer observation data
    transfer_obs_data(adata, data_paths['obs'], groups)

    # Transfer UMAP coordinates
    transfer_embedding_data(
        adata, data_paths['umap'], 'X_umap'
    )

    # Transfer spatial coordinates
    transfer_embedding_data(
        adata, data_paths['spatial'], 'spatial'
    )


def transfer_embedding_data(
    adata: anndata.AnnData,
    data_path: str,
    obsm_key: str
) -> None:
    """Transfer embedding data (UMAP/spatial) to AnnData objects."""
    import pandas as pd

    logging.info(f"Transferring {obsm_key}...")

    try:
        df = pd.read_csv(data_path, index_col=0)
        aligned_data = df.loc[adata.obs_names].values

        adata.obsm[obsm_key] = aligned_data

    except (FileNotFoundError, KeyError) as e:
        logging.warning(f"Error loading {obsm_key} data: {e}")


def transfer_obs_data(
    adata: anndata.AnnData,
    obs_path: str,
    groups: List[str]
) -> None:
    """Transfer observation data to AnnData objects."""
    import pandas as pd

    logging.info("Transferring observation data...")
    try:
        obs = pd.read_csv(obs_path, index_col=0)
        if obs.empty:
            return

        obs_aligned = obs.reindex(adata.obs.index)
        adata.obs = obs_aligned

        # Ensure group columns are strings
        for group in groups:
            if group in adata.obs.columns and adata.obs[group].dtype != object:
                adata.obs[group] = adata.obs[group].astype(str)

    except FileNotFoundError:
        logging.warning(f"File not found: {obs_path}")


def volcano_name_transform(name: str, file: str, data_type: str) -> str:
    """Transform volcano plot file names."""
    file_name = Path(file).name
    treatment = file_name.replace(f"volcanoMarkers_{data_type}_", "").replace(".csv", "")
    return f"volcano_{treatment}"
