import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path


# Load the standalone integrity helper without importing the workflow SDK.
spec = importlib.util.spec_from_file_location(
    "project_artifacts", Path(__file__).parents[1] / "wf/project_artifacts.py"
)
artifacts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(artifacts)


class UploadIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.stage = Path(self.tmp.name) / "stage"
        (self.stage / "ArrowFiles").mkdir(parents=True)
        (self.stage / "Save-ArchR-Project.rds").write_bytes(b"completed metadata")
        (self.stage / "ArrowFiles/sample.arrow").write_bytes(b"genes peaks motifs")
        self.digest = artifacts.write_project_manifest(self.stage)
        self.uploaded = Path(self.tmp.name) / "uploaded"
        shutil.copytree(self.stage, self.uploaded)

    def test_download_matches_staged_project(self):
        artifacts.verify_project_manifest(self.uploaded, self.digest)

    def test_stale_arrow_overwrite_is_rejected(self):
        (self.uploaded / "ArrowFiles/sample.arrow").write_bytes(b"genes only")
        with self.assertRaisesRegex(ValueError, "sample.arrow"):
            artifacts.verify_project_manifest(self.uploaded, self.digest)

    def test_missing_arrow_is_rejected(self):
        (self.uploaded / "ArrowFiles/sample.arrow").unlink()
        with self.assertRaisesRegex(ValueError, "sample.arrow"):
            artifacts.verify_project_manifest(self.uploaded, self.digest)

    def test_entire_stale_project_and_manifest_are_rejected(self):
        (self.uploaded / "ArrowFiles/sample.arrow").write_bytes(b"older project")
        artifacts.write_project_manifest(self.uploaded)
        with self.assertRaisesRegex(ValueError, "manifest differs"):
            artifacts.verify_project_manifest(self.uploaded, self.digest)


if __name__ == "__main__":
    unittest.main()
