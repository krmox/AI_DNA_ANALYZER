"""Regression tests for generate_mock_bam.py's variant-splicing logic.

Pins the bug fixed here: Insertion/Deletion VCF records used to be read but
never spliced into any read, so every CIGAR was a flat all-match block and
the pileup carried zero evidence for those labels. These tests fail loudly
if that regresses — a read spanning a planted insertion/deletion must show
a real 'I'/'D' CIGAR op of the correct length and the correct query-length
delta versus the reference.
"""

from __future__ import annotations

import random

import pysam
import pytest

from generate_mock_bam import (
    CIGAR_DELETION,
    CIGAR_INSERTION,
    CIGAR_MATCH,
    build_read_window,
    classify_variant,
    generate_reads,
)

# Deliberately not a repeating pattern: positions 5/6 are fixed to 'A'/'T' so
# every test below can assert exact anchor/deleted bases without ambiguity.
REFERENCE = list("GGGGG" + "AT" + "GGGG" * 8)  # 39 bases, index 5 = 'A', index 6 = 'T'


class TestClassifyVariant:
    def test_snp(self) -> None:
        assert classify_variant("A", "T") == "snp"

    def test_insertion(self) -> None:
        assert classify_variant("A", "AT") == "insertion"

    def test_deletion(self) -> None:
        assert classify_variant("AT", "A") == "deletion"

    def test_mnp_is_left_unclassified_matching_original_behaviour(self) -> None:
        # The original generator only ever handled len(ref) == len(alt) == 1;
        # any other equal-length pair (an MNP) was silently ignored. Keep
        # that exact behaviour rather than newly supporting it.
        assert classify_variant("AC", "GT") is None


class TestBuildReadWindowSnp:
    def test_alt_read_gets_substitution(self) -> None:
        variants = {5: ("A", "T")}
        assert REFERENCE[5] == "A"
        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=True)
        assert query[5] == "T"
        assert cigar == [(CIGAR_MATCH, 10)]

    def test_ref_read_is_untouched(self) -> None:
        variants = {5: ("A", "T")}
        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=False)
        assert query == "".join(REFERENCE[0:10])
        assert cigar == [(CIGAR_MATCH, 10)]


class TestBuildReadWindowInsertion:
    def test_planted_insertion_produces_real_i_op(self) -> None:
        # A -> ATT: 1-base anchor + 2 inserted bases with no reference coordinate.
        variants = {5: ("A", "ATT")}
        assert REFERENCE[5] == "A"

        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=True)

        # Anchor + preceding bases (6M), the 2 inserted bases (2I), then the
        # remaining 4 reference positions in the window (4M).
        assert cigar == [(CIGAR_MATCH, 6), (CIGAR_INSERTION, 2), (CIGAR_MATCH, 4)]
        # Reference-consuming ops (M) must still sum to the window width.
        ref_consumed = sum(length for op, length in cigar if op in (CIGAR_MATCH, CIGAR_DELETION))
        assert ref_consumed == 10
        # Query gained exactly the inserted length versus a plain reference copy.
        assert len(query) == 10 + 2
        assert query[5] == "A"
        assert query[6:8] == "TT"

    def test_non_alt_read_has_no_insertion(self) -> None:
        variants = {5: ("A", "ATT")}
        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=False)
        assert cigar == [(CIGAR_MATCH, 10)]
        assert len(query) == 10

    def test_insertion_at_window_boundary_is_skipped_not_corrupted(self) -> None:
        # ref_consumed (len(ref)=1) still fits, but a longer ref requiring
        # more positions than remain must fall back to a plain copy instead
        # of emitting a truncated / malformed CIGAR.
        variants = {9: ("AC", "ACGG")}  # len(ref)=2 needs positions 9-10; window ends at 10
        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=True)
        assert cigar == [(CIGAR_MATCH, 10)]
        assert len(query) == 10


class TestBuildReadWindowDeletion:
    def test_planted_deletion_produces_real_d_op(self) -> None:
        # AT -> A: 1-base anchor kept, 1 reference base actually deleted.
        variants = {5: ("AT", "A")}
        assert REFERENCE[5:7] == ["A", "T"]

        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=True)

        assert cigar == [(CIGAR_MATCH, 6), (CIGAR_DELETION, 1), (CIGAR_MATCH, 3)]
        ref_consumed = sum(length for op, length in cigar if op in (CIGAR_MATCH, CIGAR_DELETION))
        assert ref_consumed == 10
        # Query lost exactly the deleted length versus a plain reference copy.
        assert len(query) == 10 - 1

    def test_multi_base_deletion_length_matches_ref_alt_delta(self) -> None:
        # ACGT -> A: anchor kept, 3 bases actually deleted.
        variants = {5: ("ACGT", "A")}
        query, cigar = build_read_window(0, 15, REFERENCE, variants, carries_alt=True)

        deletion_ops = [length for op, length in cigar if op == CIGAR_DELETION]
        assert deletion_ops == [3]
        assert len(query) == 15 - 3

    def test_non_alt_read_has_no_deletion(self) -> None:
        variants = {5: ("AT", "A")}
        query, cigar = build_read_window(0, 10, REFERENCE, variants, carries_alt=False)
        assert cigar == [(CIGAR_MATCH, 10)]
        assert len(query) == 10


class TestGenerateReadsIntegration:
    """End-to-end: real pysam BAM I/O over a tiny synthetic contig."""

    @pytest.fixture()
    def small_reference(self) -> list[str]:
        random.seed(0)
        return [random.choice("ACGT") for _ in range(2000)]

    def test_regenerated_bam_has_indel_cigar_ops_matching_vcf_deltas(
        self, tmp_path, small_reference
    ) -> None:
        variants = {
            500: ("A", "ATTT"),  # insertion, +3
            1000: ("ACGT", "A"),  # deletion, -3
        }
        rng = random.Random(42)
        reads = generate_reads(small_reference, variants, read_len=150, coverage=40, rng=rng)
        assert reads, "expected at least one simulated read"

        fasta_path = tmp_path / "ref.fa"
        with open(fasta_path, "w") as handle:
            handle.write(">contig1\n")
            handle.write("".join(small_reference) + "\n")
        pysam.faidx(str(fasta_path))

        bam_path = tmp_path / "reads.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": len(small_reference), "SN": "contig1"}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
            for index, (start, query_sequence, cigar_ops) in enumerate(reads):
                alignment = pysam.AlignedSegment()
                alignment.query_name = f"read_{index}"
                alignment.query_sequence = query_sequence
                alignment.flag = 0
                alignment.reference_id = 0
                alignment.reference_start = start
                alignment.mapping_quality = 60
                alignment.cigar = cigar_ops
                alignment.query_qualities = [30] * len(query_sequence)
                out_bam.write(alignment)

        sorted_path = tmp_path / "reads.sorted.bam"
        pysam.sort("-o", str(sorted_path), str(bam_path))
        pysam.index(str(sorted_path))

        bam = pysam.AlignmentFile(str(sorted_path), "rb")
        insertion_lengths = set()
        deletion_lengths = set()
        for read in bam.fetch("contig1"):
            if not read.cigartuples:
                continue
            for op, length in read.cigartuples:
                if op == CIGAR_INSERTION:
                    insertion_lengths.add(length)
                elif op == CIGAR_DELETION:
                    deletion_lengths.add(length)
        bam.close()

        # The planted insertion is +3 bases (ATTT vs A) and the planted
        # deletion is -3 bases (ACGT vs A); with 40x coverage over a
        # 2000bp contig and a seeded RNG we expect to observe both.
        assert 3 in insertion_lengths, "no 3bp insertion CIGAR observed for the planted insertion"
        assert 3 in deletion_lengths, "no 3bp deletion CIGAR observed for the planted deletion"
