#!/usr/bin/env python3
"""Check that the complete set of figure PDFs was refreshed in this run."""
import json
from pathlib import Path
REPO = Path(__file__).absolute().parents[2]
CHECKPOINT = REPO / "logs/plotting/render_started.json"
PANELS = {"fig1": "", "fig2": "abcd", "fig3": "abcdefghij", "fig4": "abcdefghij", "fig5": "abcdefghij",
          **{f"figs{n}": "" for n in range(1, 7)}, "figs7": "abcde", "figs8": "",
          "figs9": "ab", "figs10": "abcde", "figs11": "ab", "figs12": ""}

def verify_figures():
    if not CHECKPOINT.is_file():
        raise RuntimeError("Missing render start record; run bash 06_plotting.sh.")
    started = json.loads(CHECKPOINT.read_text())["started_ns"]
    expected = [REPO / f"figures/{name}{panel}.pdf" for name, panels in PANELS.items()
                for panel in ("", *panels)]
    assert len(expected) == 65
    missing = [str(path) for path in expected if not path.is_file()]
    stale = [str(path) for path in expected if path.is_file() and path.stat().st_mtime_ns <= started]
    if missing or stale:
        raise RuntimeError(f"Missing rendered panels: {missing}; panels not refreshed: {stale}")
    CHECKPOINT.unlink()
    print(f"PASS: Rendered {len(expected)} expected PDF assets.")

if __name__ == "__main__":
    verify_figures()
