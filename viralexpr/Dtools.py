"""
Small utility functions for viralexpr pipeline.
"""

import os
import shutil


def mkdirs(dirs):
    """Create directories if not exist."""
    if isinstance(dirs, (str, os.PathLike)):
        dirs = [dirs]
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def rmdirs(dirs):
    """Remove directories if exist."""
    if isinstance(dirs, (str, os.PathLike)):
        dirs = [dirs]
    for d in dirs:
        if os.path.isdir(d):
            shutil.rmtree(d)


def mk_ckp(ckp_file):
    """Make checkpoint file. Do not overwrite if already exists."""
    if not os.path.exists(ckp_file):
        with open(ckp_file, 'w') as f:
            pass


def check_ckp(ckp_file, overwrite=False):
    """Check if checkpoint exists. Returns True if can skip."""
    if overwrite:
        return False
    return os.path.exists(ckp_file)


def test_s(cmd):
    """Test if a command exists in PATH."""
    return shutil.which(cmd) is not None


def which(cmd):
    """Return full path of command or None."""
    return shutil.which(cmd)