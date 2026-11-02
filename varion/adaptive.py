"""R3.1 depth-adaptive configuration choice. Pure functions: no I/O except loading the thresholds file.

Two indel configurations exist:
    A    joint multi-allelic genotyping
    A+B  A plus a repeat-specific error model estimated from the BAM itself
The thresholds that decide between them were chosen on development data (see docs/R3.1_THRESHOLD_SELECTION.md)
and live in thresholds.json together with their provenance. A threshold of null means "the data did not
support switching": auto mode then uses A everywhere.
"""
import json
from dataclasses import dataclass, asdict
from pathlib import Path

MODES = ("auto", "A", "A+B")
THRESHOLDS_FILE = Path(__file__).with_name("thresholds.json")


@dataclass
class Thresholds:
    global_min_depth: float | None   # contig-level median depth at/above which B is estimated; None = never in auto
    local_min_depth: int | None      # per-locus read depth at/above which a locus uses B; None = no gate
    status: str = "unknown"
    provenance: str = ""

    @classmethod
    def load(cls, path=None):
        d = json.loads(Path(path or THRESHOLDS_FILE).read_text())
        return cls(d.get("global_min_depth"), d.get("local_min_depth"), d.get("status", "unknown"), d.get("provenance", ""))


@dataclass
class Decision:
    config: str                 # "A" or "A+B"
    use_repeat_model: bool
    str_min_depth: int          # per-locus gate passed to the engine (0 = every locus)
    reason: str
    source: str                 # "manual" or "auto"
    global_depth: float | None

    def as_dict(self):
        return asdict(self)


def choose_config(mode, global_depth, thresholds):
    """Select the configuration for one contig.

    mode: "auto" | "A" | "A+B" (manual override). global_depth: median window depth of the contig, or None.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    if mode == "A":
        return Decision("A", False, 0, "manual override: A", "manual", global_depth)
    if mode == "A+B":
        return Decision("A+B", True, 0, "manual override: A+B at every locus (no depth gate)", "manual", global_depth)
    # auto
    if global_depth is None:
        return Decision("A", False, 0, "auto: depth could not be measured; using the cheaper, safer A", "auto", None)
    g = thresholds.global_min_depth
    if g is None:
        return Decision("A", False, 0, f"auto: no depth was found at which B reliably helps ({thresholds.status}); using A",
                        "auto", global_depth)
    if global_depth < g:
        return Decision("A", False, 0, f"auto: contig depth {global_depth:.1f}x < {g:g}x switching threshold; using A",
                        "auto", global_depth)
    gate = thresholds.local_min_depth or 0
    return Decision("A+B", True, int(gate),
                    f"auto: contig depth {global_depth:.1f}x >= {g:g}x; using A+B, repeat model applied at loci with >= {gate} reads",
                    "auto", global_depth)
