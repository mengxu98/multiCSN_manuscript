from pathlib import Path
import shutil
import sys
import tempfile



def publish_directory(candidate, target, preserve=()):
    """Replace a directory with rollback on failure and no retained old copy."""
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
