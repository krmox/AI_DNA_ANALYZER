"""Regression tests for generate_mock_vcf.py's REF-allele selection.

Pins the bug fixed here: REF used to be drawn at random from ACGT instead
of read from the reference FASTA, so it matched the true reference base
only ~25% of the time. For SNPs that made ALT land *on* the real reference
base at chance rate, producing labelled loci whose reads are identical to
non-variant reads — phantom variants with no evidence. These tests fail
loudly if REF ever stops being read from the FASTA.
"""

from __future__ import annotations

import random

import pysam
import pytest

from generate_mock_bam import classify_variant
from generate_mock_vcf import make_record

BASES = "ACGT"

SNP_DRAW = 0.1
INSERTION_DRAW = 0.85
DELETION_DRAW = 0.95


@pytest.fixture
def reference() -> str:
    """A fixture contig with no repeating period, so a wrong REF can't match by luck."""
    rng = random.Random(1234)
    return "".join(rng.choice(BASES) for _ in range(2000))


def _all_draws() -> list[float]:
    return [SNP_DRAW, INSERTION_DRAW, DELETION_DRAW]


class TestRefComesFromFasta:
    def test_ref_matches_fasta_at_every_position_for_every_variant_class(
        self, reference: str
    ) -> None:
        rng = random.Random(7)
        checked = 0
        for pos0 in range(0, len(reference) - 1, 7):
            for draw in _all_draws():
                record = make_record(reference, pos0, draw, rng)
                assert record is not None, "pure-ACGT fixture must never be skipped"
                ref, _ = record
                assert ref == reference[pos0 : pos0 + len(ref)]
                checked += 1
        # The old bug matched at ~25%; a partial pass is still a failure.
        assert checked > 500

    def test_snp_alt_is_never_the_reference_base(self, reference: str) -> None:
        rng = random.Random(11)
        for pos0 in range(len(reference)):
            ref, alt = make_record(reference, pos0, SNP_DRAW, rng)
            assert classify_variant(ref, alt) == "snp"
            assert alt != reference[pos0], "ALT on the reference base is a phantom variant"

    def test_insertion_anchor_is_the_reference_base(self, reference: str) -> None:
        rng = random.Random(12)
        for pos0 in range(0, len(reference), 13):
            ref, alt = make_record(reference, pos0, INSERTION_DRAW, rng)
            assert classify_variant(ref, alt) == "insertion"
            assert ref == reference[pos0]
            assert alt.startswith(ref)

    def test_deletion_ref_spans_two_real_reference_bases(self, reference: str) -> None:
        rng = random.Random(13)
        for pos0 in range(0, len(reference) - 1, 13):
            ref, alt = make_record(reference, pos0, DELETION_DRAW, rng)
            assert classify_variant(ref, alt) == "deletion"
            assert ref == reference[pos0 : pos0 + 2]
            assert alt == ref[0]


class TestNonAcgtPositionsAreSkipped:
    """A non-ACGT REF cannot survive the BAM generator's ACGT sanitisation,
    so such sites must be dropped rather than emitted as evidence-free labels."""

    def test_n_base_is_skipped_for_every_variant_class(self) -> None:
        rng = random.Random(3)
        sequence = "ACNGT"
        for draw in _all_draws():
            assert make_record(sequence, 2, draw, rng) is None

    def test_deletion_spanning_an_n_is_skipped(self) -> None:
        assert make_record("ACNGT", 1, DELETION_DRAW, random.Random(4)) is None

    def test_deletion_at_the_last_base_is_skipped(self) -> None:
        assert make_record("ACGT", 3, DELETION_DRAW, random.Random(5)) is None


class TestGeneratedVcfFixture:
    """End-to-end: the checked-in mock VCF's REF must match the checked-in FASTA
    at 100% of records — the exact figure the phantom-variant diagnosis measured."""

    def test_every_record_ref_matches_the_fasta(self) -> None:
        fasta_path = "data/giab_chr21/chr21_slice.fa"
        vcf_path = "data/giab_chr21/chr21_slice.vcf.gz"
        try:
            fasta = pysam.FastaFile(fasta_path)
        except (OSError, ValueError):
            pytest.skip(f"{fasta_path} not present in this checkout")

        with fasta:
            sequence = fasta.fetch(fasta.references[0])
            mismatches = []
            phantom_snps = 0
            total = 0
            with pysam.VariantFile(vcf_path) as vcf:
                for record in vcf.fetch():
                    total += 1
                    pos0 = record.pos - 1  # VCF POS is 1-based
                    if sequence[pos0 : pos0 + len(record.ref)].upper() != record.ref:
                        mismatches.append(record.pos)
                    if (
                        classify_variant(record.ref, record.alts[0]) == "snp"
                        and record.alts[0] == sequence[pos0].upper()
                    ):
                        phantom_snps += 1

        assert total > 0
        assert not mismatches, f"{len(mismatches)}/{total} REF alleles disagree with the FASTA"
        assert phantom_snps == 0, f"{phantom_snps}/{total} SNPs place ALT on the reference base"
