from __future__ import annotations

import os
import sys
from pathlib import Path


def _resolve_dsrnaseeker_root() -> Path:
    """Resolve the dsRNASeeker source tree without hard-coding private HPC paths.

    Local/HPC development may override with DSRNASEEKER_ROOT. Public deployments
    default to the checked-out repository root: <repo>/webapp/app/bridge.py -> <repo>.
    """
    here = Path(__file__).resolve()
    repo_root = here.parents[2]
    return Path(os.environ.get("DSRNASEEKER_ROOT", str(repo_root))).resolve()


def load_dsrnaseeker():
    root = _resolve_dsrnaseeker_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    try:
        from modules.priority import add_priority_columns
        from modules.supervised_benchmark import fixed_grouped_oof
    except Exception as exc:
        raise RuntimeError(
            "Could not import dsRNASeeker from the configured source tree. "
            "For local development, set DSRNASEEKER_ROOT to the dsRNASeeker repository."
        ) from exc

    return root, add_priority_columns, fixed_grouped_oof
