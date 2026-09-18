"""Integrity checks for the final ArchRProject upload."""

import hashlib
import json
from pathlib import Path


MANIFEST = "archr_upload_manifest.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_project_manifest(project: Path) -> str:
    files = {
        str(path.relative_to(project)): _sha256(path)
        for path in sorted(project.rglob("*"))
        if path.is_file() and path != project / MANIFEST
    }
    if "Save-ArchR-Project.rds" not in files or not any(
        name.startswith("ArrowFiles/") for name in files
    ):
        raise ValueError("Cannot publish an incomplete ArchRProject.")
    (project / MANIFEST).write_text(json.dumps(files, indent=2) + "\n")
    return _sha256(project / MANIFEST)


def verify_project_manifest(project: Path, expected_manifest_hash: str) -> None:
    if _sha256(project / MANIFEST) != expected_manifest_hash:
        raise ValueError("Uploaded project manifest differs from the publisher's manifest.")
    files = json.loads((project / MANIFEST).read_text())
    if "Save-ArchR-Project.rds" not in files or not any(
        name.startswith("ArrowFiles/") for name in files
    ):
        raise ValueError("Uploaded manifest does not describe an ArchRProject.")
    for name, expected in files.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Invalid path in project manifest: {name}")
        path = project / relative
        if not path.is_file() or _sha256(path) != expected:
            raise ValueError(f"Uploaded ArchRProject file is missing or changed: {name}")
