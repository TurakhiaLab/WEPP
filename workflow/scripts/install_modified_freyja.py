#!/usr/bin/env python3
"""Install WEPP's modified Freyja (src/Freyja) into the active conda env once.

The marker and lock files live inside $CONDA_PREFIX, so every project that
shares the same conda env via `--conda-prefix` sees the same state:
the first run installs, subsequent runs from any project are no-ops.

Lock uses fcntl.flock, which works identically on Linux and macOS.

Usage: install_modified_freyja.py <path-to-src/Freyja>
"""
import fcntl, os, subprocess, sys

conda_prefix = os.environ.get("CONDA_PREFIX") or sys.exit(
    "CONDA_PREFIX is empty; run snakemake with --use-conda"
)
freyja_dir = sys.argv[1]
marker = os.path.join(conda_prefix, ".wepp_freyja_installed")
lock   = os.path.join(conda_prefix, ".wepp_freyja_install.lock")

with open(lock, "w") as fh:
    fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
    if not os.path.exists(marker):
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--no-deps", "-e", freyja_dir]
        )
        open(marker, "w").close()
