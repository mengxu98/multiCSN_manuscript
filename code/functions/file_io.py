import json
from pathlib import Path
import shutil
import sys
import tempfile

from runtime import sha256


def directory_hashes(root, exclude=()):
    """Fingerprint regular files; never follow links outside an artifact tree."""
    root = Path(root)
    hashes = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if relative.parts[0] in exclude:
            continue
        if path.is_symlink():
            raise ValueError(f"Artifact contains a symbolic link: {path}")
        if path.is_file():
            hashes[str(relative)] = sha256(path)
    return hashes


def publish_directory(candidate, target, backup_root, preserve=()):
    """Back up and verify the old tree, then replace it with rollback on failure."""
    candidate, target = Path(candidate), Path(target)
    if candidate.is_symlink() or target.is_symlink():
        raise ValueError("Publication requires real directories")
    if not candidate.is_dir() or candidate.stat().st_dev != target.parent.stat().st_dev:
        raise ValueError("Prepare the candidate on the target filesystem")
    if any((candidate / name).exists() for name in preserve):
        raise ValueError("Candidate contains a directory reserved for preservation")
    if not target.exists():
        candidate.rename(target)
        return None

    original = directory_hashes(target, exclude=preserve)
    backup_root = Path(backup_root)
    backup_root.mkdir(parents=True, exist_ok=True)
    backup = Path(tempfile.mkdtemp(prefix="previous-", dir=backup_root))
    previous = backup / "contents"
    previous.mkdir()
    for entry in target.iterdir():
        if entry.name in preserve:
            continue
        if entry.is_dir():
            shutil.copytree(entry, previous / entry.name)
        else:
            shutil.copy2(entry, previous / entry.name)
    if directory_hashes(previous) != original or directory_hashes(target, exclude=preserve) != original:
        raise RuntimeError(f"Backup verification failed: {backup}")
    (backup / "sha256.json").write_text(json.dumps(original, indent=2) + "\n")

    rollback = Path(tempfile.mkdtemp(prefix=f".{target.name}-rollback-", dir=target.parent))
    rollback.rmdir()
    moved = []
    try:
        target.rename(rollback)
        for name in preserve:
            if (rollback / name).exists():
                (rollback / name).rename(candidate / name)
                moved.append(name)
        candidate.rename(target)
    except BaseException:
        for name in reversed(moved):
            (candidate / name).rename(rollback / name)
        if rollback.exists():
            rollback.rename(target)
        raise
    try:
        shutil.rmtree(rollback)
    except OSError as error:
        print(f"Published {target}; old temporary tree retained at {rollback}: {error}", file=sys.stderr)
    return backup


def write_artifact_manifest(root):
    import pandas as pd

    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if path.name == "artifact_manifest.tsv":
            continue
        rows.append({
            "file": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    pd.DataFrame(rows).to_csv(root / "artifact_manifest.tsv", sep="\t", index=False)
