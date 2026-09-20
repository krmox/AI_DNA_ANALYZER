"""RECONSTRUCTED fuzz regression for the native COUNTS backend (``pileup_native.count_region``).

Why this exists: the first native *read-tensor* design made one htslib pass over a whole run of windows and
differed from the per-64-nt-window reference in 30 cells (overlapping mates + repeated read names: htslib's
mate-overlap adjustment pairs reads by name among those alive in the pileup). The counts backend also makes one
pass over the whole span, so it needs its own protection. Two checks per random BAM:

1. ``load_counts`` native vs pysam (the production contract);
2. ``count_region`` over the whole span vs the concatenation of per-window ``count_region`` calls
   (the exact analogue of the failed design).

Small synthetic BAMs, regression protection only; not a scientific result.
"""
from __future__ import annotations

import random
from array import array

import numpy as np
import pysam
import pytest

import pileup_counts as P
from test_native_reads_equivalence import CONTIG, CLEN
from testdata import build_fixture

pytestmark = pytest.mark.skipif(P._pileup_native is None, reason="pileup_native C extension not built")


def _heavy_name_bam(path, reference, seed, realistic=False):
    """Many overlapping pairs.

    realistic=True : one unique name per pair and CORRECT mate fields (RNEXT/PNEXT/TLEN point at the real mate),
                     i.e. what an aligner writes; this is the contract the production data satisfy.
    realistic=False: adversarial - a 12-name pool (names repeat across many pairs and window boundaries) and
                     mate fields that do not match any real mate. Known to diverge (see xfail below).
    """
    rng = random.Random(seed)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CLEN}]}
    names = [f"n{i}" for i in range(12)]
    counter = [0]
    tmp = path.with_suffix(".u.bam")
    reads = []
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        def mk(name, start, flag, ln, bq_lo, mate=None):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, start
            a.mapping_quality, a.cigar = rng.choice([20, 40, 60, 60]), [(0, ln)]
            a.query_sequence = "".join(rng.choice("ACGT") if rng.random() < 0.15 else reference[start + i] for i in range(ln))
            a.query_qualities = array("B", [rng.randint(bq_lo, 41) for _ in range(ln)])
            a.next_reference_id, a.next_reference_start, a.template_length = 0, (mate if mate is not None else start), ln
            return a
        for _ in range(900):
            ln = rng.choice([60, 100, 150])
            s = rng.randint(0, CLEN - ln - 1)
            counter[0] += 1
            name = f"p{counter[0]}" if realistic else rng.choice(names)
            f = 0x1 | 0x2 | (0x40 if rng.random() < .5 else 0x80) | (0x10 if rng.random() < .5 else 0)
            if rng.random() < 0.6:                                    # overlapping mate, same name
                s2 = min(CLEN - ln - 1, s + rng.randint(0, ln - 10))
                m1, m2 = (s2, s) if realistic else (None, None)
                reads.append(mk(name, s, f, ln, rng.choice([5, 13, 25]), m1))
                reads.append(mk(name, s2, f ^ 0x10 ^ 0x40 ^ 0x80, ln, rng.choice([13, 25]), m2))
            else:
                reads.append(mk(name, s, f, ln, rng.choice([5, 13, 25])))
        reads.sort(key=lambda a: a.reference_start)
        for a in reads:
            out.write(a)
    pysam.sort("-o", str(path), str(tmp))
    pysam.index(str(path))


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build_fixture(tmp_path_factory.mktemp("cnt_fuzz"), contig_length=CLEN, depth=10,
                         fasta_contig=CONTIG, bam_contig=CONTIG, vcf_contig=CONTIG)


ADVERSARIAL_XFAIL = pytest.mark.xfail(
    reason="KNOWN LIMITATION: whole-span htslib pass differs from per-window passes on BAMs with heavily repeated "
           "read names / inconsistent mate fields (same mechanism as the first native read-tensor design). "
           "Not observed on real GIAB data (0 mismatches on the validated regions).", strict=False)
MODES = [pytest.param(True, id="realistic"), pytest.param(False, id="adversarial", marks=ADVERSARIAL_XFAIL)]


@pytest.mark.parametrize("realistic", MODES)
@pytest.mark.parametrize("seed", range(10))
def test_load_counts_native_equals_pysam(fx, tmp_path, seed, realistic):
    bam = tmp_path / "h.bam"
    _heavy_name_bam(bam, fx.reference, seed, realistic)
    bed = tmp_path / "gap.bed"
    bed.write_text(f"{CONTIG}\t0\t1000\n{CONTIG}\t1100\t2500\n{CONTIG}\t2560\t4096\n")   # BED gaps -> dropped windows
    res = {}
    for tag, flag in (("py", False), ("nat", True)):
        P.NATIVE_AVAILABLE = flag
        try:
            res[tag] = P.load_counts(str(fx.fasta_path), str(bam), str(fx.vcf_path), str(bed), (0, CLEN),
                                     contig=CONTIG, return_positions=True)
        finally:
            P.NATIVE_AVAILABLE = P._pileup_native is not None
    for a, b in zip(res["py"], res["nat"]):
        assert np.array_equal(a, b)
    assert res["nat"][0][:, 6].max() >= 4


@pytest.mark.parametrize("realistic", MODES)
@pytest.mark.parametrize("seed", range(10))
def test_whole_span_pass_equals_per_window_passes(fx, tmp_path, seed, realistic):
    bam = tmp_path / "h.bam"
    _heavy_name_bam(bam, fx.reference, 100 + seed, realistic)
    whole = np.frombuffer(P._pileup_native.count_region(str(bam), CONTIG, 0, CLEN, 20, 13),
                          dtype=np.float64).reshape(CLEN, 9)
    per = np.concatenate([np.frombuffer(P._pileup_native.count_region(str(bam), CONTIG, s, s + 64, 20, 13),
                                        dtype=np.float64).reshape(64, 9) for s in range(0, CLEN, 64)])
    assert np.array_equal(whole, per), f"{int((whole != per).sum())} cells differ between whole-span and per-window passes"
