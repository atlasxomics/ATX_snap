"""Regression tests for duplicated labels in motif result tables.

Run with: python -m unittest discover -s tests -p 'test_features_h5ad.py'
"""

import importlib.util
from pathlib import Path
import tempfile
import unittest

import anndata
import numpy as np
import pandas as pd


# Load the helpers without importing the Latch workflow entrypoint.
spec = importlib.util.spec_from_file_location(
    "features", Path(__file__).resolve().parents[1] / "wf" / "features.py"
)
features = importlib.util.module_from_spec(spec)
spec.loader.exec_module(features)


class H5ADSanitizationTests(unittest.TestCase):
    def test_reserved_sample_name_full_and_reduced_round_trip(self):
        frame = pd.DataFrame(
            {"run1": [1.0, 2.0], "run2": [3.0, 4.0], "_index-1": [5.0, 6.0]},
            index=["motif1", "motif2"],
        )
        mapped = features._remap_sample_labels_in_df(
            frame, {"run1": "_index", "run2": "_index"}
        )
        adata = anndata.AnnData(X=np.zeros((2, 1)))
        adata.uns["motifs"] = {"results": mapped}
        expected = frame.copy()
        expected.columns = ["_index-3", "_index-2", "_index-1"]
        with tempfile.TemporaryDirectory() as directory:
            features.save_anndata_objects(adata, "_motifs", Path(directory))
            for name in ("combined_motifs.h5ad", "combined_sm_motifs.h5ad"):
                restored = anndata.read_h5ad(Path(directory) / name)
                pd.testing.assert_frame_equal(
                    restored.uns["motifs"]["results"], expected
                )
        features._sanitize_dataframe_for_h5ad(mapped)
        pd.testing.assert_frame_equal(mapped, expected)

    def test_reserved_column_rename_avoids_named_index(self):
        frame = pd.DataFrame(
            {"_index": ["A", "B"], "_index-1": [10, 20]},
            index=pd.Index(["motif1", "motif2"], name="_index-2"),
        )
        features._sanitize_dataframe_for_h5ad(frame)
        self.assertEqual(frame.columns.tolist(), ["_index-3", "_index-1"])
        self.assertEqual(frame["_index-3"].tolist(), ["A", "B"])
        self.assertEqual(frame.index.name, "_index-2")
        adata = anndata.AnnData(X=np.zeros((2, 1)))
        adata.uns["results"] = frame
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "motifs.h5ad"
            adata.write_h5ad(path)
            restored = anndata.read_h5ad(path)
        pd.testing.assert_frame_equal(restored.uns["results"], frame)

    def test_duplicate_object_columns_preserve_values_and_types(self):
        frame = pd.DataFrame(
            [[1, "A", 10], [None, "B", 20]],
            columns=["sample", "sample", "sample-1"],
            dtype=object,
        )
        features._sanitize_dataframe_for_h5ad(frame)
        self.assertEqual(frame.columns.tolist(), ["sample", "sample-2", "sample-1"])
        self.assertEqual(frame["sample"].iloc[0], 1)
        self.assertTrue(pd.isna(frame["sample"].iloc[1]))
        self.assertEqual(frame["sample-2"].tolist(), ["A", "B"])
        self.assertEqual(frame["sample-1"].tolist(), [10, 20])
        self.assertTrue(pd.api.types.is_numeric_dtype(frame["sample"]))
        expected = frame.copy()
        features._sanitize_dataframe_for_h5ad(frame)
        pd.testing.assert_frame_equal(frame, expected)

    def test_shared_sample_names_preserve_separate_runs(self):
        frame = pd.DataFrame(
            [["run1", "run2", "unmapped"]],
            columns=["run1", "run2", "shared-1"],
        )
        original = frame.copy()
        mapped = features._remap_sample_labels_in_df(
            frame, {"run1": "shared", "run2": "shared"}
        )
        self.assertEqual(mapped.columns.tolist(), ["shared", "shared-2", "shared-1"])
        self.assertEqual(mapped.iloc[0].tolist(), ["shared", "shared", "unmapped"])
        pd.testing.assert_frame_equal(frame, original)

    def test_nested_uns_round_trip_with_numeric_and_object_duplicates(self):
        frame = pd.DataFrame({
            "first": [1.0, 2.0],
            "second": [3.0, 4.0],
            "label": ["A", "B"],
            "other_label": ["C", "D"],
        }, index=["motif1", "motif2"])
        frame.columns = ["sample", "sample", "label", "label"]
        adata = anndata.AnnData(X=np.zeros((2, 1)))
        adata.uns["motifs"] = {"results": frame}
        features._sanitize_uns_for_h5ad(adata.uns)
        self.assertEqual(frame.columns.tolist(), ["sample", "sample-1", "label", "label-1"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "motifs.h5ad"
            adata.write_h5ad(path)
            restored = anndata.read_h5ad(path)
        pd.testing.assert_frame_equal(restored.uns["motifs"]["results"], frame)

    def test_unique_columns_and_empty_duplicate_columns(self):
        frame = pd.DataFrame({"score": [1.0], "name": ["motif"]})
        original = frame.copy()
        features._sanitize_dataframe_for_h5ad(frame)
        pd.testing.assert_frame_equal(frame, original)
        empty = pd.DataFrame(columns=["sample", "sample"])
        features._sanitize_dataframe_for_h5ad(empty)
        self.assertEqual(empty.columns.tolist(), ["sample", "sample-1"])


if __name__ == "__main__":
    unittest.main()
