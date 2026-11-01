"""pytest hook: make the repository-root V1 modules importable (same logic as _bootstrap.py)."""
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
