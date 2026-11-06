"""Exercise input integrity and recovery without loading the Latch entrypoint."""

import ast
import logging
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
from typing import Optional
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from test_features_h5ad import features


# These helpers only need the standard library and features, not Flyte/Snap.
task_path = Path(__file__).resolve().parents[1] / "wf" / "task.py"
tree = ast.parse(task_path.read_text())
tree.body = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name in {"_ensure_dense_h5ad_x", "_rebuild_gene_h5ad_from_rds"}
]
namespace = dict(
    Path=Path, Optional=Optional, logging=logging, os=os,
    subprocess=subprocess, ft=features,
)
exec(compile(tree, str(task_path), "exec"), namespace)
task = SimpleNamespace(**namespace)


def write_fixture(path, corrupt_dataset=None):
    with h5py.File(path, "w") as handle:
        for name in ("X", "obsm/spatial"):
            dataset = handle.create_dataset(
                name, data=np.arange(64).reshape(16, 4),
                chunks=(4, 4), fletcher32=True,
            )
            if name == corrupt_dataset:
                offset = dataset.id.get_chunk_info_by_coord((12, 0)).byte_offset
        handle.create_dataset("uns/scalar", data=1)
        handle.create_dataset("uns/empty", shape=(0, 4))
        handle.create_dataset("uns/null", dtype="f")
    if corrupt_dataset:
        # Corrupt a late chunk while leaving the file and X metadata readable.
        with open(path, "r+b") as stream:
            stream.seek(offset)
            first_byte = stream.read(1)
            stream.seek(offset)
            stream.write(bytes([first_byte[0] ^ 0xFF]))


class GeneH5ADRecoveryTests(unittest.TestCase):
    def test_reads_every_element_with_bounded_buffers(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "valid.h5ad"
            write_fixture(path)
            reads = []
            original = h5py.Dataset.__getitem__

            def read(dataset, selection):
                result = original(dataset, selection)
                if dataset.name == "/X":
                    reads.append(result)
                return result

            with patch.object(h5py.Dataset, "__getitem__", read):
                features.validate_h5ad_readable(path, max_block_bytes=16)
            self.assertTrue(all(array.nbytes <= 16 for array in reads))
            self.assertEqual(sum(array.size for array in reads), 64)
            np.testing.assert_array_equal(
                np.sort(np.concatenate([array.ravel() for array in reads])),
                np.arange(64),
            )

    def test_late_corruption_names_file_and_dataset(self):
        for dataset in ("X", "obsm/spatial"):
            with self.subTest(dataset=dataset), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "broken.h5ad"
                write_fixture(path, corrupt_dataset=dataset)
                with h5py.File(path, "r") as handle:
                    self.assertIsInstance(handle["X"], h5py.Dataset)
                with self.assertRaises(OSError) as raised:
                    task._ensure_dense_h5ad_x(path)
                self.assertIn(str(path), str(raised.exception))
                self.assertIn(f"dataset /{dataset}", str(raised.exception))

    def test_recovers_from_rds_and_validates_before_replacement(self):
        for replacement_corrupt in (False, True):
            with self.subTest(corrupt=replacement_corrupt), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "sample.h5ad"
                rds = Path(directory) / "sample.rds"
                rds.touch()
                write_fixture(path, corrupt_dataset="X")
                original = path.read_bytes()

                def rebuild(command, check):
                    self.assertTrue(check)
                    self.assertEqual(command[2], str(rds))
                    write_fixture(
                        command[3], "X" if replacement_corrupt else None
                    )

                with patch.object(subprocess, "run", side_effect=rebuild) as run:
                    if replacement_corrupt:
                        with self.assertRaises(OSError):
                            task._ensure_dense_h5ad_x(path, rds)
                        self.assertEqual(path.read_bytes(), original)
                    else:
                        task._ensure_dense_h5ad_x(path, rds)
                        features.validate_h5ad_readable(path)
                    run.assert_called_once()
                self.assertFalse(list(Path(directory).glob("gene_h5ad_rebuild_*")))


if __name__ == "__main__":
    unittest.main()
