"""Regression suite for the reads_native C/htslib acceleration of load_reads().

Every test compares the read-level tensor and labels from the native path against the original
pure-Python ``window_reads`` path, byte for byte, through the public ``load_reads`` entry point.
Skipped (not failed) when the extension is not built. A randomised "fuzz" BAM exercises every
admission/encoding branch of ``ReadLevelPileupProvider._read_row``: duplicates, QC-fail,
secondary, supplementary, low MAPQ, orphans, overlapping mates, tied read names, insertions,
deletions, ref-skips, soft clips, missing (0xff) qualities, missing NM tags, differently-typed
integer NM tags, low base qualities, N bases, both strands, and depth > MAX_READS
(hash-ordered truncation).
"""
from __future__ import annotations

import hashlib
import random
from array import array
from pathlib import Path

import numpy as np
import pysam
import pytest

import read_level_pileup as R
from testdata import build_fixture

pytestmark = pytest.mark.skipif(
    not R.NATIVE_READS_AVAILABLE, reason="reads_native C extension not built")

CONTIG, CLEN = "chr21", 4096


def _both(fasta, bam, vcf, bed, start, end, max_reads=48):
    out = {}
    for tag, flag in (("py", False), ("native", True)):
        R.NATIVE_READS_AVAILABLE = flag
        try:
            out[tag] = R.load_reads(str(fasta), str(bam), str(vcf), str(bed) if bed else None,
                                    (start, end), contig=CONTIG, max_reads=max_reads)
        finally:
            R.NATIVE_READS_AVAILABLE = R._reads_native is not None
    return out["py"], out["native"]


def _assert_identical(py, nat):
    (pr, pl), (nr, nl) = py, nat
    assert pr.shape == nr.shape and pr.dtype == nr.dtype == np.uint8
    assert np.array_equal(pr, nr), f"{int((pr != nr).sum())} tensor cells differ"
    assert np.array_equal(pl, nl)


def _fuzz_bam(path: Path, reference: str, seed: int, nm_float: bool = False, no_seq: bool = False,
              unique_names: bool = False) -> None:
    rng = random.Random(seed)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CLEN}]}
    unsorted = path.with_suffix(".unsorted.bam")
    names = [f"q{i}" for i in range(400)]
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        def mk(name, start, flag, mapq, ops, mate=None):
            a = pysam.AlignedSegment(out.header)
            seq, qref = [], start
            cig = []
            for op, ln in ops:
                cig.append((op, ln))
                if op in (0, 1, 4):
                    for _ in range(ln):
                        if op == 0:
                            b = reference[qref] if qref < CLEN else "A"
                            qref += 1
                        else:
                            b = rng.choice("ACGT")
                        r = rng.random()
                        if r < 0.03:
                            b = rng.choice("ACGT")
                        elif r < 0.04:
                            b = "N"
                        seq.append(b)
                elif op in (2, 3):
                    qref += ln
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, start
            a.mapping_quality = mapq
            a.cigar = cig
            a.query_sequence = "".join(seq)
            if rng.random() < 0.05:
                a.query_qualities = None            # 0xff -> quality 0
            else:
                lo = rng.choice([2, 10, 12, 13, 14, 30])
                a.query_qualities = array("B", [rng.randint(lo, 41) if rng.random() > 0.05 else rng.randint(0, 12)
                                                for _ in seq])
            t = rng.random()
            if nm_float and t < 0.3:
                a.set_tag("NM", float(rng.randint(0, 4)), "f")
            elif t < 0.7:
                a.set_tag("NM", rng.randint(0, 6), rng.choice(["C", "S", "I", "c", "s", "i"]))
            if mate:
                a.next_reference_id, a.next_reference_start = 0, mate[0]
                a.template_length = mate[1]
            return a

        reads = []
        for i in range(1500):
            rl = rng.choice([50, 100, 150, 250])
            start = rng.randint(0, CLEN - rl - 1)
            if rng.random() < 0.3:                       # dense cluster -> depth > 48
                start = rng.randint(1000, 1100)
            ops = [(0, rl)]
            r = rng.random()
            if r < 0.12 and rl > 60:
                ops = [(0, 20), (2, rng.randint(1, 5)), (0, rl - 20)]
            elif r < 0.22 and rl > 60:
                ops = [(0, 25), (1, rng.randint(1, 4)), (0, rl - 25)]
            elif r < 0.27 and rl > 60:
                ops = [(0, 30), (3, rng.randint(5, 40)), (0, rl - 30)]
            elif r < 0.37:
                ops = [(4, rng.randint(2, 10)), (0, rl - 10)]
            flag = 0x1 | 0x2 | (0x40 if rng.random() < .5 else 0x80)
            if rng.random() < 0.5:
                flag |= 0x10
            x = rng.random()
            if x < 0.03: flag |= 0x400
            elif x < 0.05: flag |= 0x200
            elif x < 0.07: flag |= 0x100
            elif x < 0.09: flag |= 0x800
            elif x < 0.12: flag &= ~0x2                  # orphan (paired, not proper)
            elif x < 0.14: flag &= ~0x1 & ~0x2           # unpaired
            mapq = rng.choice([0, 19, 20, 20, 30, 60, 60, 60])
            name = f"u{i}" if unique_names else rng.choice(names)   # name reuse -> tied hash keys
            reads.append(mk(name, start, flag, mapq, ops, mate=(max(0, start - 30), 200)))
            if rng.random() < 0.3:                       # overlapping mate, same qname
                s2 = min(CLEN - 101, start + rng.randint(0, 60))
                reads.append(mk(name, s2, flag ^ 0x10 ^ 0x40 ^ 0x80, mapq, [(0, 100)], mate=(start, -200)))
        if no_seq:
            a = mk("noseq", 500, 0, 60, [(0, 50)])
            a.query_sequence = None
            reads.append(a)
        reads.sort(key=lambda a: a.reference_start)
        for a in reads:
            out.write(a)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build_fixture(tmp_path_factory.mktemp("reads_eq"), contig_length=CLEN, depth=10,
                         fasta_contig=CONTIG, bam_contig=CONTIG, vcf_contig=CONTIG)


class TestFuzz:
    @pytest.mark.parametrize("seed", range(8))
    def test_fuzz_bam(self, fx, tmp_path, seed):
        bam = tmp_path / f"fuzz{seed}.bam"
        _fuzz_bam(bam, fx.reference, seed)
        py, nat = _both(fx.fasta_path, bam, fx.vcf_path, fx.bed_path, 0, CLEN)
        _assert_identical(py, nat)
        reads = py[0]
        assert int(((reads[..., 6] & 1) > 0).sum(axis=1).max()) == 48, "fixture must exceed MAX_READS somewhere"
        assert int((reads[..., 6] & R.FLAG_GAP).astype(bool).sum()) > 0
        assert int((reads[..., 6] & R.FLAG_INSERTION).astype(bool).sum()) > 0
        assert int((reads[..., 6] & R.FLAG_DELETION_NEXT).astype(bool).sum()) > 0

    def test_small_max_reads(self, fx, tmp_path):
        bam = tmp_path / "fz.bam"
        _fuzz_bam(bam, fx.reference, 99)
        py, nat = _both(fx.fasta_path, bam, fx.vcf_path, fx.bed_path, 0, CLEN, max_reads=7)
        _assert_identical(py, nat)

    def test_bed_gaps_and_partial_regions(self, fx, tmp_path):
        bam = tmp_path / "fz.bam"
        _fuzz_bam(bam, fx.reference, 5)
        bed = tmp_path / "gap.bed"
        bed.write_text(f"{CONTIG}\t0\t640\n{CONTIG}\t704\t2048\n{CONTIG}\t2100\t3300\n")
        for (s, e) in [(0, 1536), (100, 3000), (1000, 1200), (2200, 3400)]:
            py, nat = _both(fx.fasta_path, bam, fx.vcf_path, bed, s, e)
            _assert_identical(py, nat)

    @pytest.mark.parametrize("threads", [2, 4])
    def test_thread_count_does_not_change_bytes(self, fx, tmp_path, monkeypatch, threads):
        bam = tmp_path / "fz.bam"
        _fuzz_bam(bam, fx.reference, 11)
        monkeypatch.setenv("AI_DNA_ANALYZER_READ_THREADS", str(threads))
        py, nat = _both(fx.fasta_path, bam, fx.vcf_path, fx.bed_path, 0, CLEN)
        _assert_identical(py, nat)


class TestFallback:
    def test_float_nm_falls_back_to_python_and_stays_identical(self, fx, tmp_path):
        import reads_native
        bam = tmp_path / "f.bam"
        _fuzz_bam(bam, fx.reference, 3, nm_float=True)
        with pytest.raises(RuntimeError):
            reads_native.read_windows(str(bam), CONTIG, np.asarray([0, 64, 128, 192], dtype=np.int64),
                                      64, 20, 13, 48, 1)
        py, nat = _both(fx.fasta_path, bam, fx.vcf_path, fx.bed_path, 0, CLEN)
        _assert_identical(py, nat)

    def test_missing_sequence_read_defers_to_python_behaviour(self, fx, tmp_path):
        # A primary read with SEQ='*' makes the original Python path raise TypeError;
        # the native path must not silently "fix" it.
        bam = tmp_path / "ns.bam"
        _fuzz_bam(bam, fx.reference, 4, no_seq=True)
        for flag in (False, True):
            R.NATIVE_READS_AVAILABLE = flag
            try:
                with pytest.raises(TypeError):
                    R.load_reads(str(fx.fasta_path), str(bam), str(fx.vcf_path), str(fx.bed_path),
                                 (0, CLEN), contig=CONTIG)
            finally:
                R.NATIVE_READS_AVAILABLE = R._reads_native is not None


class TestHash:
    def test_blake2b_key_matches_hashlib(self):
        import reads_native
        rng = random.Random(1)
        for length in list(range(0, 260)) + [500, 1000]:
            for _ in range(5):
                b = bytes(rng.randrange(256) for _ in range(length))
                assert reads_native.hash_name(b) == int.from_bytes(
                    hashlib.blake2b(b, digest_size=8).digest(), "big")


import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[0])  # project root (path-independent)
REAL = Path(_ROOT + "/data/giab_hg002_chr21_12Mb")


@pytest.mark.skipif(not (REAL / "hg002_chr21_32_44M_15x.bam").exists(), reason="real HG002 slice absent")
class TestRealData:
    @pytest.mark.parametrize("start", [32_000_000, 36_500_000])
    def test_real_hg002_15x(self, start):
        args = (_ROOT + "/data/reference/chr21_full.fa",
                str(REAL / "hg002_chr21_32_44M_15x.bam"), str(REAL / "hg002_chr21_32_44M.vcf.gz"),
                str(REAL / "hg002_chr21_32_44M_highconf.bed"))
        out = {}
        for tag, flag in (("py", False), ("native", True)):
            R.NATIVE_READS_AVAILABLE = flag
            try:
                out[tag] = R.load_reads(*args, (start, start + 20_000), contig="chr21")
            finally:
                R.NATIVE_READS_AVAILABLE = R._reads_native is not None
        _assert_identical(out["py"], out["native"])
