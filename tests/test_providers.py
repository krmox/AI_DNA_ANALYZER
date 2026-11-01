"""Regression tests for the real-data provider, run against pysam fixtures.

These exercise genuine file I/O: a real ``.fai``, a real BAM index, a real
tabix index, real CIGAR strings, and a deliberate ``chr`` prefix mismatch
between the VCF and the other two files.

What they establish is that the *parser* is correct. They say nothing about
whether the model performs on GIAB — see ``testdata.py`` for the gap between
a fixture and real chr21.
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import pytest
import torch

pysam = pytest.importorskip("pysam", reason="pysam is required for provider tests")

from bimamba_variant_caller.config import (  # noqa: E402
    GAP_ID,
    LABEL_DELETION,
    LABEL_INSERTION,
    LABEL_NORMAL,
    LABEL_SNP,
)
from bimamba_variant_caller.loss import weights_from_counts  # noqa: E402
from bimamba_variant_caller.providers import (  # noqa: E402
    GenomicWindow,
    GiabAlignmentProvider,
    ProviderDataset,
    load_bed_regions,
    summarise_label_distribution,
)
from bimamba_variant_caller.testdata import build_fixture  # noqa: E402


@pytest.fixture(scope="module")
def fixture(tmp_path_factory: pytest.TempPathFactory):
    """Build one shared fixture; generation is slow relative to the tests."""
    directory = tmp_path_factory.mktemp("giab")
    return build_fixture(directory, contig_length=3000, depth=20, variant_spacing=150)


@pytest.fixture()
def provider(fixture) -> GiabAlignmentProvider:
    """An opened provider over the fixture's variant-bearing region."""
    instance = GiabAlignmentProvider(
        fasta_path=fixture.fasta_path,
        bam_path=fixture.bam_path,
        vcf_path=fixture.vcf_path,
        seq_len=64,
        contig="chr21",
        high_confidence_bed=fixture.bed_path,
        region=(200, 2800),
        min_depth=3,
    )
    instance.open()
    yield instance
    instance.close()


class TestGenomicWindow:
    """Interval invariants."""

    def test_length_is_half_open(self) -> None:
        assert len(GenomicWindow("chr21", 100, 164)) == 64

    @pytest.mark.parametrize("bounds", [(-1, 10), (50, 50), (60, 50)])
    def test_rejects_invalid_intervals(self, bounds: Tuple[int, int]) -> None:
        with pytest.raises(ValueError):
            GenomicWindow("chr21", *bounds)


class TestContigNormalisation:
    """The fixture uses '21' in the VCF and 'chr21' elsewhere on purpose."""

    def test_resolves_prefixed_and_bare_spellings(self) -> None:
        resolve = GiabAlignmentProvider._normalise_contig
        assert resolve("chr21", ["chr21", "chr22"]) == "chr21"
        assert resolve("chr21", ["21", "22"]) == "21"
        assert resolve("21", ["chr21"]) == "chr21"

    def test_unknown_contig_raises_loudly(self) -> None:
        """A silent mismatch yields an empty dataset and a degenerate model."""
        with pytest.raises(ValueError, match="not found"):
            GiabAlignmentProvider._normalise_contig("chrX", ["chr21", "chr22"])

    def test_provider_reconciles_across_files(self, provider: GiabAlignmentProvider) -> None:
        assert provider._fasta_contig == "chr21"
        assert provider._bam_contig == "chr21"
        assert provider._vcf_contig == "21"


class TestReferenceFetch:
    """FASTA extraction and coordinate alignment."""

    def test_reference_matches_the_source_sequence(self, provider, fixture) -> None:
        window = GenomicWindow("chr21", 500, 564)
        assert provider._fetch_reference(window) == fixture.reference[500:564]

    def test_reference_is_uppercased(self, provider) -> None:
        sequence = provider._fetch_reference(GenomicWindow("chr21", 300, 364))
        assert sequence == sequence.upper()

    def test_reference_is_padded_past_the_contig_end(self, provider) -> None:
        sequence = provider._fetch_reference(GenomicWindow("chr21", 2960, 3024))
        assert len(sequence) == 64
        assert sequence.endswith("N")


class TestLabelFetch:
    """VCF parsing, including the 1-based to 0-based conversion."""

    def test_planted_variants_are_labelled_at_the_right_offset(
        self, provider, fixture
    ) -> None:
        """Labels must land on the correct 0-based coordinate.

        SNPs and insertions are marked on the record position itself;
        deletions start one base later, because the VCF anchor is not
        deleted.
        """
        expected = {"snp": LABEL_SNP, "insertion": LABEL_INSERTION, "deletion": LABEL_DELETION}

        checked = 0
        for variant in fixture.variants:
            if not 200 <= variant.position < 2700:
                continue
            window = GenomicWindow("chr21", variant.position - 10, variant.position + 54)
            labels = provider._fetch_labels(window)
            checked += 1

            marked_at = 11 if variant.kind == "deletion" else 10
            assert labels[marked_at] == expected[variant.kind], (
                f"{variant.kind} at {variant.position} landed on the wrong offset"
            )
            if variant.kind == "deletion":
                assert labels[10] == LABEL_NORMAL, "the deletion anchor is not deleted"
        assert checked > 0

    def test_variant_free_windows_are_all_normal(self, provider, fixture) -> None:
        occupied = {variant.position for variant in fixture.variants}
        for start in range(210, 2600, 64):
            window = GenomicWindow("chr21", start, start + 64)
            if any(start <= position < start + 70 for position in occupied):
                continue
            assert set(provider._fetch_labels(window)) == {LABEL_NORMAL}

    def test_deletions_skip_the_anchor_base(self) -> None:
        """Regression: the anchor in 'ACGT -> A' is present, not deleted.

        Labelling from the record position marked a base that was never
        deleted and missed the last one that was, which pinned the anchor
        column's VAF at 0 and diluted the deletion class.
        """
        label, offset, span = GiabAlignmentProvider._classify_allele("ACGT", "A")
        assert label == LABEL_DELETION
        assert offset == 1, "deletion labels must start after the anchor"
        assert span == 3

    def test_insertions_mark_only_the_anchor(self) -> None:
        label, offset, span = GiabAlignmentProvider._classify_allele("A", "ACGT")
        assert label == LABEL_INSERTION
        assert (offset, span) == (0, 1)

    def test_mnp_decomposes_across_its_span(self) -> None:
        label, offset, span = GiabAlignmentProvider._classify_allele("AC", "GT")
        assert label == LABEL_SNP
        assert (offset, span) == (0, 2)

    def test_snp_covers_exactly_one_base(self) -> None:
        assert GiabAlignmentProvider._classify_allele("A", "G") == (LABEL_SNP, 0, 1)

    def test_empty_alleles_are_ignored(self) -> None:
        assert GiabAlignmentProvider._classify_allele("", "A") == (None, 0, 0)


class TestPileupConsensus:
    """The consensus, VAF and depth logic that makes real data usable."""

    def test_snp_loci_always_differ_from_the_reference(self, provider) -> None:
        """Establishes the pileup is calling variants, not echoing reference."""
        checked = 0
        for window in provider:
            snp = window.labels == LABEL_SNP
            if not bool(snp.any()):
                continue
            checked += 1
            assert bool((window.input_ids[snp] != window.reference_ids[snp]).all())
        assert checked > 0, "expected at least one SNP in the fixture region"

    def test_background_loci_agree_with_the_reference(self, provider) -> None:
        """Sequencing errors must be filtered out by the VAF threshold.

        At 1% error and 20x depth an error reaches VAF ~0.01, far below the
        0.15 threshold, so the consensus should stay on the reference base.
        """
        mismatched = 0
        total = 0
        for window in provider:
            normal = window.labels == LABEL_NORMAL
            total += int(normal.sum())
            mismatched += int((window.input_ids[normal] != window.reference_ids[normal]).sum())
        assert total > 0
        assert mismatched / total < 0.05

    def test_vaf_separates_variant_from_background(self, provider) -> None:
        """Heterozygous sites sit near 0.5; background near 0."""
        variant_vaf: list[float] = []
        background_vaf: list[float] = []
        for window in provider:
            assert window.vaf is not None
            snp = window.labels == LABEL_SNP
            normal = window.labels == LABEL_NORMAL
            variant_vaf += window.vaf[snp].tolist()
            background_vaf += window.vaf[normal].tolist()

        assert variant_vaf and background_vaf
        assert sum(variant_vaf) / len(variant_vaf) > 0.3
        assert sum(background_vaf) / len(background_vaf) < 0.1

    def test_deletion_gaps_produce_non_zero_vaf(self, provider) -> None:
        """Regression: gaps were excluded from VAF, pinning deletions at 0.

        A read spanning a locus with a deletion positively asserts the base
        is absent, so it is variant evidence rather than missing data.
        """
        deletion_vaf: list[float] = []
        for window in provider:
            assert window.vaf is not None
            deletion = window.labels == LABEL_DELETION
            if bool(deletion.any()):
                deletion_vaf += window.vaf[deletion].tolist()

        assert deletion_vaf, "expected deletions in the fixture region"
        mean = sum(deletion_vaf) / len(deletion_vaf)
        assert mean > 0.2, f"deletion VAF collapsed to {mean:.3f}"
        assert max(deletion_vaf) > 0.0

    def test_all_variant_classes_share_a_vaf_scale(self, provider) -> None:
        """One threshold must serve substitutions and indels alike."""
        by_class: dict[int, list[float]] = {
            LABEL_SNP: [],
            LABEL_INSERTION: [],
            LABEL_DELETION: [],
            LABEL_NORMAL: [],
        }
        for window in provider:
            assert window.vaf is not None
            for label, value in zip(window.labels.tolist(), window.vaf.tolist()):
                by_class[int(label)].append(value)

        background = by_class[LABEL_NORMAL]
        assert sum(background) / len(background) < 0.1

        for label in (LABEL_SNP, LABEL_INSERTION, LABEL_DELETION):
            values = by_class[label]
            assert values, f"no columns for class {label}"
            assert sum(values) / len(values) > 0.2

    def test_deletion_support_property(self) -> None:
        from bimamba_variant_caller.providers import LocusPileup

        locus = LocusPileup(
            depth=10,
            base_counts={"A": 6, "-": 4},
            consensus_base="A",
            reference_base="A",
            vaf=0.4,
            quality=35.0,
            insertion_sequences=[],
        )
        assert locus.deletion_support == pytest.approx(0.4)

        empty = LocusPileup(0, {}, "N", "A", 0.0, 0.0, [])
        assert empty.deletion_support == 0.0

    def test_depth_is_populated_and_positive(self, provider) -> None:
        """Regression: depth was inert in the single-read representation."""
        for window in provider:
            assert float(window.depth.max()) > 0.0
            assert float(window.depth.mean()) > 1.0

    def test_low_mapq_reads_are_excluded(self, fixture) -> None:
        """MAPQ filtering is the main defence against repeat-region calls."""
        strict = GiabAlignmentProvider(
            fixture.fasta_path,
            fixture.bam_path,
            fixture.vcf_path,
            seq_len=64,
            contig="chr21",
            region=(200, 600),
            min_mapping_quality=99,  # above the fixture's MAPQ 60
            min_depth=0,
        )
        with strict:
            assert all(float(window.depth.max()) == 0.0 for window in strict)

    def test_min_depth_produces_n_not_a_guess(self, fixture) -> None:
        """Below the depth floor the caller must abstain, not invent a base."""
        shallow = GiabAlignmentProvider(
            fixture.fasta_path,
            fixture.bam_path,
            fixture.vcf_path,
            seq_len=64,
            contig="chr21",
            region=(200, 600),
            min_depth=10_000,
        )
        with shallow:
            window = shallow[0]
            assert bool((window.input_ids == 0).all())


class TestGapPadding:
    """Column alignment must survive indels."""

    def test_all_channels_share_one_width(self, provider) -> None:
        for window in provider:
            assert window.reference_ids.shape == (64,)
            assert window.input_ids.shape == (64,)
            assert window.labels.shape == (64,)
            assert window.base_quality.shape == (64,)
            assert window.depth.shape == (64,)
            assert window.vaf is not None and window.vaf.shape == (64,)

    def test_insertion_columns_are_gapped_in_the_reference(self, provider) -> None:
        """Regression: without this the downstream frame shifts."""
        found = False
        for window in provider:
            gaps = window.reference_ids == GAP_ID
            if bool(gaps.any()):
                found = True
                # A reference gap means an inserted base, never another gap.
                assert bool((window.input_ids[gaps] != GAP_ID).all())
        assert found, "expected at least one insertion column in the fixture"

    def test_gap_columns_carry_no_depth_signal_confusion(self, provider) -> None:
        for window in provider:
            gaps = window.input_ids == GAP_ID
            if not bool(gaps.any()):
                continue
            assert bool((window.base_quality[gaps] == 0.0).all())


class TestProviderPlumbing:
    """Windowing, BED filtering and the Dataset adapter."""

    def test_windows_tile_the_region_by_stride(self, provider) -> None:
        windows = [window.window for window in provider]
        assert windows[0] is not None
        for earlier, later in zip(windows, windows[1:]):
            assert later.start - earlier.start == provider.stride

    def test_max_windows_caps_the_scan(self, fixture) -> None:
        capped = GiabAlignmentProvider(
            fixture.fasta_path,
            fixture.bam_path,
            fixture.vcf_path,
            seq_len=64,
            contig="chr21",
            region=(200, 2800),
            max_windows=3,
        )
        with capped:
            assert len(capped) == 3

    def test_bed_filtering_excludes_outside_windows(self, fixture, tmp_path: Path) -> None:
        narrow = tmp_path / "narrow.bed"
        narrow.write_text("chr21\t400\t800\n")

        restricted = GiabAlignmentProvider(
            fixture.fasta_path,
            fixture.bam_path,
            fixture.vcf_path,
            seq_len=64,
            contig="chr21",
            region=(200, 2800),
            high_confidence_bed=narrow,
        )
        with restricted:
            for window in restricted:
                assert window.window is not None
                assert 400 <= window.window.start and window.window.end <= 800

    def test_out_of_range_index_raises(self, provider) -> None:
        with pytest.raises(IndexError):
            provider[len(provider)]

    def test_dataset_adapter_yields_model_ready_samples(self, provider) -> None:
        dataset = ProviderDataset(provider)
        sample = dataset[0]
        for key in ("reference_ids", "input_ids", "labels", "base_quality", "depth"):
            assert sample[key].shape == (64,)

    def test_rejects_invalid_thresholds(self, fixture) -> None:
        with pytest.raises(ValueError, match="vaf_threshold"):
            GiabAlignmentProvider(
                fixture.fasta_path, fixture.bam_path, fixture.vcf_path, vaf_threshold=1.5
            )
        with pytest.raises(ValueError, match="stride"):
            GiabAlignmentProvider(
                fixture.fasta_path, fixture.bam_path, fixture.vcf_path, stride=0
            )


class TestBedLoading:
    """BED is already 0-based half-open; no conversion may be applied."""

    def test_intervals_are_read_verbatim(self, tmp_path: Path) -> None:
        bed = tmp_path / "regions.bed"
        bed.write_text("chr21\t100\t200\nchr21\t300\t400\n")
        assert load_bed_regions(bed, "chr21") == [(100, 200), (300, 400)]

    def test_overlapping_intervals_are_merged(self, tmp_path: Path) -> None:
        bed = tmp_path / "overlap.bed"
        bed.write_text("chr21\t100\t250\nchr21\t200\t400\n")
        assert load_bed_regions(bed, "chr21") == [(100, 400)]

    def test_contig_prefix_is_ignored_when_matching(self, tmp_path: Path) -> None:
        bed = tmp_path / "bare.bed"
        bed.write_text("21\t100\t200\n")
        assert load_bed_regions(bed, "chr21") == [(100, 200)]

    def test_other_contigs_and_headers_are_skipped(self, tmp_path: Path) -> None:
        bed = tmp_path / "mixed.bed"
        bed.write_text("track name=x\n#comment\nchr22\t0\t500\nchr21\t100\t200\n")
        assert load_bed_regions(bed, "chr21") == [(100, 200)]


class TestImbalanceWeighting:
    """Weight calibration at real genomic density."""

    def test_sqrt_weighting_tames_the_ratio(self) -> None:
        """At 1:1000, raw inverse frequency is what caused false positives."""
        counts = torch.tensor([1_000_000.0, 700.0, 150.0, 150.0])
        raw = weights_from_counts(counts, alpha=1.0, max_weight_cap=1e9)
        smoothed = weights_from_counts(counts, alpha=0.5, max_weight_cap=1e9)
        assert float(raw.max() / raw.min()) > 1000
        assert float(smoothed.max() / smoothed.min()) < 200

    def test_cap_binds_on_unnormalised_weights(self) -> None:
        """The spec's formula has no normalisation, so the cap is active."""
        counts = torch.tensor([1_000_000.0, 700.0, 150.0, 150.0])
        capped = weights_from_counts(counts, alpha=0.5, max_weight_cap=10.0, normalise=False)
        assert float(capped.max()) == pytest.approx(10.0)

    def test_alpha_zero_disables_weighting(self) -> None:
        counts = torch.tensor([1_000_000.0, 700.0, 150.0, 150.0])
        weights = weights_from_counts(counts, alpha=0.0, max_weight_cap=10.0)
        assert torch.allclose(weights, torch.ones(4), atol=1e-5)

    def test_absent_class_does_not_divide_by_zero(self) -> None:
        counts = torch.tensor([1000.0, 0.0, 10.0, 10.0])
        assert torch.isfinite(weights_from_counts(counts)).all()

    def test_rejects_invalid_arguments(self) -> None:
        with pytest.raises(ValueError, match="alpha"):
            weights_from_counts(torch.ones(4), alpha=-1.0)
        with pytest.raises(ValueError, match="max_weight_cap"):
            weights_from_counts(torch.ones(4), max_weight_cap=0.0)


class TestDiagnosticSummary:
    """The counting utility behind the diagnostic script."""

    def test_counts_by_class_name(self) -> None:
        counts = summarise_label_distribution([0, 0, 1, 2, 3, 3])
        assert counts == {"Normal": 2, "SNP": 1, "Insertion": 1, "Deletion": 2}

    def test_absent_classes_report_zero(self) -> None:
        assert summarise_label_distribution([0, 0])["SNP"] == 0
