"""Put the repository root on ``sys.path`` so scripts and tests in this directory can import the
V1 modules that stay at the repository root (``config``, ``cascade``, ``pileup_counts`` ...).

Importing this module has no other effect.  Run scripts from the repository root, e.g.
``python3 DIR/script.py``: relative data/results paths ("data/...", "results/...") are resolved
against the current directory, exactly as before the files were moved.
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_root = _here
while not os.path.exists(os.path.join(_root, "binomial_baseline.py")):
    _parent = os.path.dirname(_root)
    if _parent == _root:
        raise RuntimeError("repository root (binomial_baseline.py) not found above " + _here)
    _root = _parent
if _root not in sys.path:
    sys.path.insert(0, _root)
