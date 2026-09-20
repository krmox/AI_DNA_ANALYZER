"""RECONSTRUCTED extraction-integrity regression tests (coordinate bug C1 and related invariants).

NOT the original ``test_hg005_extraction_integrity.py`` (lost with a wiped temporary worktree). The original
HG005 cache is not available, so these use synthetic fixtures (testdata.py) and assert invariants only; no
expected outputs are copied from any historical result.

Invariants: positions are true genomic coordinates (never ``chunk_start + row_index`` after a dropped window);
reference_index equals the FASTA base at every position; truth-SNP REF equals the FASTA base; shapes; determinism;
native == pysam.
"""
from __future__ import annotations

import numpy as np
import pysam
import pytest

import pileup_counts as P
from config import NUCLEOTIDES
from testdata import build_fixture

CONTIG = "chr21"
TOKENS = ("N",) + NUCLEOTIDES


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build_fixture(tmp_path_factory.mktemp("integrity"), contig_length=4096, depth=15,
                         fasta_contig=CONTIG, bam_contig=CONTIG, vcf_contig=CONTIG)


@pytest.fixture()
def gap_bed(tmp_path):
    bed = tmp_path / "gap.bed"
    bed.write_text(f"{CONTIG}\t0\t300\n{CONTIG}\t400\t1500\n{CONTIG}\t1600\t4096\n")   # drops windows 256-319, 384-447, 1472-1599
    return bed


def _run(fx, bed, native, region=(0, 4096)):
    P.NATIVE_AVAILABLE = native and P._pileup_native is not None
    try:
        return P.load_counts(str(fx.fasta_path), str(fx.bam_path), str(fx.vcf_path), str(bed), region,
                             contig=CONTIG, return_positions=True)
    finally:
        P.NATIVE_AVAILABLE = P._pileup_native is not None


@pytest.mark.parametrize("native", [False, True])
class TestIntegrity:
    def test_positions_are_true_coordinates_not_row_index(self, fx, gap_bed, native):
        counts, labels, pos = _run(fx, gap_bed, native)
        assert pos.shape == (len(counts),) == labels.shape
        assert (np.diff(pos) > 0).all()
        assert (pos % 64 == np.arange(len(pos)) % 64).all(), "positions must stay window-aligned"
        naive = pos[0] + np.arange(len(pos))
        assert not np.array_equal(pos, naive), "fixture must contain dropped windows, else this test proves nothing"
        for lo, hi in [(300 // 64 * 64, 320), (384, 448)]:          # dropped windows are absent
            assert not ((pos >= lo) & (pos < hi)).any()

    def test_reference_index_matches_fasta_at_every_position(self, fx, gap_bed, native):
        counts, _, pos = _run(fx, gap_bed, native)
        fa = pysam.FastaFile(str(fx.fasta_path))
        ref = fa.fetch(CONTIG, 0, 4096).upper()
        want = np.array([TOKENS.index(ref[p]) if ref[p] in TOKENS else 0 for p in pos], dtype=np.float64)
        assert np.array_equal(counts[:, 9], want)

    def test_truth_snp_ref_matches_fasta_and_labels_land_on_it(self, fx, gap_bed, native):
        counts, labels, pos = _run(fx, gap_bed, native)
        snps = [v for v in fx.variants if v.kind == "snp"]
        assert snps
        hit = 0
        for v in snps:
            assert fx.reference[v.position].upper() == v.ref_allele.upper()
            idx = np.searchsorted(pos, v.position)
            if idx < len(pos) and pos[idx] == v.position:
                assert labels[idx] != 0, "truth SNP locus must carry a non-Normal label"
                hit += 1
        assert hit > 0
        assert (labels != 0).sum() >= hit

    def test_dimensions_and_dtypes(self, fx, gap_bed, native):
        counts, labels, pos = _run(fx, gap_bed, native)
        assert counts.shape[1] == P.COUNT_DIM == 10 and counts.dtype == np.float64
        assert labels.dtype == np.int64 and pos.dtype == np.int64
        assert len(counts) % 64 == 0

    def test_deterministic(self, fx, gap_bed, native):
        a, b = _run(fx, gap_bed, native), _run(fx, gap_bed, native)
        assert all(np.array_equal(x, y) for x, y in zip(a, b))

    def test_all_windows_dropped_returns_empty_shaped_arrays(self, fx, tmp_path, native):
        bed = tmp_path / "tiny.bed"
        bed.write_text(f"{CONTIG}\t10\t20\n")
        counts, labels, pos = _run(fx, bed, native)
        assert counts.shape == (0, P.COUNT_DIM) and labels.shape == (0,) and pos.shape == (0,)


@pytest.mark.skipif(P._pileup_native is None, reason="pileup_native not built")
def test_native_equals_python_with_dropped_windows(fx, gap_bed):
    a, b = _run(fx, gap_bed, False), _run(fx, gap_bed, True)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
