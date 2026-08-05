"""Unit tests for the BiMamba variant caller package.

Run with::

    python -m pytest bimamba_variant_caller/tests -q

The tests that matter most are the dataset alignment ones: every past bug in
this project came from the three channels drifting out of lockstep or from
an event being attributed to the wrong index.
"""

from __future__ import annotations

import pytest
import torch

from bimamba_variant_caller.config import (
    GAP_ID,
    LABEL_DELETION,
    LABEL_INSERTION,
    LABEL_NORMAL,
    LABEL_SNP,
    DataConfig,
    ExperimentConfig,
    ModelConfig,
    TrainingConfig,
)
from bimamba_variant_caller.dataset import (
    DNATokenizer,
    SyntheticVariantDataset,
    build_dataloaders,
)
from bimamba_variant_caller.loss import (
    FocalLoss,
    compute_calibrated_class_weights,
    estimate_class_counts,
)
from bimamba_variant_caller.model import BiMambaBlock, VariantCaller
from bimamba_variant_caller.train import (
    EarlyStopping,
    Trainer,
    build_scheduler,
    compute_metrics,
    validate,
)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


class TestConfig:
    """Configuration dataclasses reject invalid experiments at construction."""

    def test_defaults_are_valid(self) -> None:
        config = ExperimentConfig()
        assert config.model.vocab_size == 6
        assert config.model.num_classes == 4
        assert config.data.val_seed is not None, "validation split must be deterministic"

    @pytest.mark.parametrize("rate", [-0.1, 1.5])
    def test_rejects_out_of_range_mutation_rate(self, rate: float) -> None:
        with pytest.raises(ValueError, match="mutation_rate"):
            DataConfig(mutation_rate=rate)

    def test_rejects_unknown_fusion(self) -> None:
        with pytest.raises(ValueError, match="fusion"):
            ModelConfig(fusion="attention")  # type: ignore[arg-type]

    def test_rejects_warmup_longer_than_training(self) -> None:
        with pytest.raises(ValueError, match="warmup_epochs"):
            TrainingConfig(num_epochs=3, warmup_epochs=5)

    def test_rejects_label_smoothing_of_one(self) -> None:
        with pytest.raises(ValueError, match="label_smoothing"):
            TrainingConfig(label_smoothing=1.0)


# ---------------------------------------------------------------------------
# Tokenizer
# ---------------------------------------------------------------------------


class TestDNATokenizer:
    """Encoding and decoding round-trip over the five-letter alphabet."""

    def test_encode_maps_alphabet_to_expected_ids(self) -> None:
        tokenizer = DNATokenizer()
        assert tokenizer.encode("NACGT").tolist() == [0, 1, 2, 3, 4]

    def test_round_trip_preserves_sequence(self) -> None:
        tokenizer = DNATokenizer()
        sequence = "ACGTNNACGT"
        assert tokenizer.decode(tokenizer.encode(sequence)) == sequence

    def test_unknown_characters_degrade_to_n(self) -> None:
        tokenizer = DNATokenizer()
        assert tokenizer.encode("AXG").tolist() == [1, 0, 3]

    def test_decode_accepts_plain_sequence(self) -> None:
        assert DNATokenizer().decode([1, 2, 3, 4]) == "ACGT"


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


class TestSyntheticVariantDataset:
    """The generator must keep all three channels index-aligned."""

    def test_shapes_are_exactly_seq_len(self) -> None:
        dataset = SyntheticVariantDataset(20, seq_len=32, mutation_rate=0.3, seed=1)
        for index in range(len(dataset)):
            sample = dataset[index]
            assert sample["input_ids"].shape == (32,)
            assert sample["reference_ids"].shape == (32,)
            assert sample["labels"].shape == (32,)

    def test_insertions_do_not_change_emitted_length(self) -> None:
        """Insertions emit two tokens, so truncation must hold the length."""
        dataset = SyntheticVariantDataset(30, seq_len=16, mutation_rate=0.9, seed=2)
        for index in range(len(dataset)):
            assert dataset[index]["input_ids"].numel() == 16

    def test_labels_stay_in_valid_range(self) -> None:
        dataset = SyntheticVariantDataset(20, seq_len=32, mutation_rate=0.4, seed=3)
        for index in range(len(dataset)):
            labels = dataset[index]["labels"]
            assert int(labels.min()) >= 0
            assert int(labels.max()) <= 3

    def test_zero_mutation_rate_yields_identical_channels(self) -> None:
        dataset = SyntheticVariantDataset(10, seq_len=32, mutation_rate=0.0, seed=4, noise_rate=0.0)
        for index in range(len(dataset)):
            sample = dataset[index]
            assert torch.equal(sample["input_ids"], sample["reference_ids"])
            assert torch.equal(sample["labels"], torch.zeros(32, dtype=torch.long))

    def test_normal_columns_match_reference_absent_noise(self) -> None:
        """Without frame shifts or noise, a Normal label implies a match.

        Under gap padding this now holds even with indels present, which is
        the whole point of the alignment upgrade.
        """
        dataset = SyntheticVariantDataset(
            40, seq_len=48, mutation_rate=0.25, seed=5, noise_rate=0.0
        )
        for index in range(len(dataset)):
            sample = dataset[index]
            normal = sample["labels"] == LABEL_NORMAL
            observed = sample["input_ids"][normal]
            reference = sample["reference_ids"][normal]
            # Padding is (N, N) and matches too, so equality must hold throughout.
            assert torch.equal(observed, reference)

    def test_snp_columns_mismatch_reference(self) -> None:
        """SNPs exclude the reference base, so the column always mismatches."""
        dataset = SyntheticVariantDataset(80, seq_len=64, mutation_rate=0.1, seed=6, noise_rate=0.0)
        checked = 0
        for index in range(len(dataset)):
            sample = dataset[index]
            snp = sample["labels"] == LABEL_SNP
            if not bool(snp.any()):
                continue
            checked += 1
            assert bool((sample["input_ids"][snp] != sample["reference_ids"][snp]).all())
        assert checked > 0, "expected at least one SNP to check"

    def test_deletion_labels_are_emitted(self) -> None:
        """Deletions defer their label to the next emitted token.

        They no longer imply an N in the observed channel: the gap is a
        frame shift, and only a trailing-edge deletion emits an explicit N.
        """
        dataset = SyntheticVariantDataset(40, seq_len=48, mutation_rate=0.3, seed=7)
        assert any(
            bool((dataset[index]["labels"] == LABEL_DELETION).any()) for index in range(len(dataset))
        )

    def test_insertion_zones_are_fully_labelled(self) -> None:
        """Every base of the inserted zone carries the Insertion label."""
        dataset = SyntheticVariantDataset(
            40, seq_len=64, mutation_rate=0.2, seed=8, max_insertion_len=6
        )
        seen_multi_base_zone = False
        for index in range(len(dataset)):
            labels = dataset[index]["labels"].tolist()
            run = 0
            for label in labels + [LABEL_NORMAL]:
                if label == LABEL_INSERTION:
                    run += 1
                    continue
                if run > 1:
                    seen_multi_base_zone = True
                run = 0
        assert seen_multi_base_zone, "expected at least one multi-base insertion zone"

    def test_insertion_columns_are_gapped_in_the_reference(self) -> None:
        """Gap-padded alignment: an insertion column has GAP in reference."""
        dataset = SyntheticVariantDataset(
            40, seq_len=64, mutation_rate=0.25, seed=21, max_insertion_len=6
        )
        checked = 0
        for index in range(len(dataset)):
            sample = dataset[index]
            insertion = sample["labels"] == LABEL_INSERTION
            if not bool(insertion.any()):
                continue
            checked += 1
            assert bool((sample["reference_ids"][insertion] == GAP_ID).all())
            assert bool((sample["input_ids"][insertion] != GAP_ID).all())
        assert checked > 0

    def test_deletion_columns_are_gapped_in_the_observed(self) -> None:
        """Gap-padded alignment: a deletion column has GAP in observed."""
        dataset = SyntheticVariantDataset(
            40, seq_len=64, mutation_rate=0.25, seed=23, max_insertion_len=6
        )
        checked = 0
        for index in range(len(dataset)):
            sample = dataset[index]
            deletion = sample["labels"] == LABEL_DELETION
            if not bool(deletion.any()):
                continue
            checked += 1
            assert bool((sample["input_ids"][deletion] == GAP_ID).all())
            assert bool((sample["reference_ids"][deletion] != GAP_ID).all())
        assert checked > 0

    def test_no_frame_shift_remains(self) -> None:
        """The regression this upgrade exists to fix.

        Under gap padding, a Normal column must agree with its reference
        counterpart unless a sequencing error hit it. Previously ~42% of
        Normal columns mismatched because indels shifted the frame.
        """
        dataset = SyntheticVariantDataset(
            60, seq_len=64, mutation_rate=0.2, seed=22, max_insertion_len=6, noise_rate=0.0
        )
        for index in range(len(dataset)):
            sample = dataset[index]
            normal = sample["labels"] == LABEL_NORMAL
            assert torch.equal(sample["input_ids"][normal], sample["reference_ids"][normal])

    def test_channels_are_column_aligned(self) -> None:
        """All five channels must describe the same columns."""
        dataset = SyntheticVariantDataset(20, seq_len=48, mutation_rate=0.2, seed=24)
        for index in range(len(dataset)):
            sample = dataset[index]
            for key in ("reference_ids", "input_ids", "labels", "base_quality", "depth"):
                assert sample[key].shape == (48,)

    def test_noise_is_never_labelled(self) -> None:
        """Sequencing errors must not appear in the label array.

        With indels disabled, every mismatch on a Normal position is by
        construction a sequencing error, and must still read as Normal.
        """
        noisy = SyntheticVariantDataset(
            40, seq_len=64, mutation_rate=0.0, seed=31, noise_rate=0.05
        )
        saw_error = False
        for index in range(len(noisy)):
            sample = noisy[index]
            assert torch.equal(sample["labels"], torch.zeros(64, dtype=torch.long))
            if bool((sample["input_ids"] != sample["reference_ids"]).any()):
                saw_error = True
        assert saw_error, "expected the noise model to introduce some errors"

    def test_zero_noise_rate_leaves_sequence_untouched(self) -> None:
        clean = SyntheticVariantDataset(10, seq_len=48, mutation_rate=0.0, seed=32, noise_rate=0.0)
        for index in range(len(clean)):
            sample = clean[index]
            assert torch.equal(sample["input_ids"], sample["reference_ids"])

    def test_noise_rate_is_approximately_honoured(self) -> None:
        rate = 0.05
        dataset = SyntheticVariantDataset(
            200, seq_len=64, mutation_rate=0.0, seed=33, noise_rate=rate
        )
        mismatches = sum(
            int((dataset[i]["input_ids"] != dataset[i]["reference_ids"]).sum()) for i in range(200)
        )
        observed = mismatches / (200 * 64)
        assert rate * 0.7 < observed < rate * 1.3

    def test_rejects_invalid_noise_and_insertion_bounds(self) -> None:
        with pytest.raises(ValueError, match="noise_rate"):
            SyntheticVariantDataset(4, seq_len=8, mutation_rate=0.1, noise_rate=1.5)
        with pytest.raises(ValueError, match="max_insertion_len"):
            SyntheticVariantDataset(4, seq_len=8, mutation_rate=0.1, max_insertion_len=0)

    def test_errors_get_lower_quality_than_correct_bases(self) -> None:
        """The planted signal that makes noise separable from a true SNP.

        Without this correlation the quality channel is decoration: a
        sequencing error and a SNP are identical in the base channels.
        """
        dataset = SyntheticVariantDataset(
            120, seq_len=64, mutation_rate=0.0, seed=41, noise_rate=0.15
        )
        error_quality: List[float] = []
        clean_quality: List[float] = []

        for index in range(len(dataset)):
            sample = dataset[index]
            is_error = sample["input_ids"] != sample["reference_ids"]
            real = sample["input_ids"] != 0
            error_quality += sample["base_quality"][is_error & real].tolist()
            clean_quality += sample["base_quality"][~is_error & real].tolist()

        assert error_quality and clean_quality
        mean_error = sum(error_quality) / len(error_quality)
        mean_clean = sum(clean_quality) / len(clean_quality)
        assert mean_error < mean_clean - 10.0

    def test_gap_columns_carry_zero_quality_and_depth(self) -> None:
        """A gap has no base call, so it has no confidence or coverage."""
        dataset = SyntheticVariantDataset(30, seq_len=64, mutation_rate=0.3, seed=42)
        for index in range(len(dataset)):
            sample = dataset[index]
            gaps = sample["input_ids"] == GAP_ID
            if not bool(gaps.any()):
                continue
            assert bool((sample["base_quality"][gaps] == 0.0).all())
            assert bool((sample["depth"][gaps] == 0.0).all())

    def test_quality_stays_within_phred_bounds(self) -> None:
        dataset = SyntheticVariantDataset(
            30, seq_len=64, mutation_rate=0.1, seed=43, noise_rate=0.05, max_phred=60
        )
        for index in range(len(dataset)):
            quality = dataset[index]["base_quality"]
            assert float(quality.min()) >= 0.0
            assert float(quality.max()) <= 60.0

    def test_depth_is_positive_on_real_base_calls(self) -> None:
        dataset = SyntheticVariantDataset(20, seq_len=64, mutation_rate=0.1, seed=44)
        for index in range(len(dataset)):
            sample = dataset[index]
            real = (sample["input_ids"] != 0) & (sample["input_ids"] != GAP_ID)
            assert bool((sample["depth"][real] > 0).all())

    def test_seeded_dataset_is_stable_across_repeated_passes(self) -> None:
        """Regression: val_loss is only comparable if the split never moves.

        A single shared RNG is reproducible across runs but advances on every
        access, so a second pass returned different sequences and every epoch
        silently validated on fresh data.
        """
        dataset = SyntheticVariantDataset(6, seq_len=32, mutation_rate=0.2, seed=123)
        first_pass = [dataset[index]["input_ids"].clone() for index in range(6)]
        second_pass = [dataset[index]["input_ids"].clone() for index in range(6)]
        assert all(torch.equal(a, b) for a, b in zip(first_pass, second_pass))

    def test_seeded_access_is_order_independent(self) -> None:
        """Shuffled traversal must not change what a given index returns."""
        dataset = SyntheticVariantDataset(6, seq_len=32, mutation_rate=0.2, seed=123)
        sequential = [dataset[index]["labels"].clone() for index in range(6)]
        reverse = {index: dataset[index]["labels"].clone() for index in reversed(range(6))}
        assert all(torch.equal(sequential[index], reverse[index]) for index in range(6))

    def test_distinct_indices_yield_distinct_samples(self) -> None:
        """Per-index seeding must not collapse the dataset into one sample."""
        dataset = SyntheticVariantDataset(16, seq_len=32, mutation_rate=0.2, seed=123)
        unique = {tuple(dataset[index]["input_ids"].tolist()) for index in range(16)}
        assert len(unique) == 16

    def test_out_of_range_index_raises(self) -> None:
        dataset = SyntheticVariantDataset(4, seq_len=16, mutation_rate=0.1, seed=10)
        with pytest.raises(IndexError):
            dataset[4]

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"num_samples": 0, "seq_len": 8, "mutation_rate": 0.1},
            {"num_samples": 4, "seq_len": 0, "mutation_rate": 0.1},
            {"num_samples": 4, "seq_len": 8, "mutation_rate": 2.0},
        ],
    )
    def test_rejects_invalid_construction(self, kwargs: dict) -> None:
        with pytest.raises(ValueError):
            SyntheticVariantDataset(**kwargs)

    def test_build_dataloaders_returns_batched_tensors(self) -> None:
        config = DataConfig(train_samples=8, val_samples=8, seq_len=16, batch_size=4)
        train_loader, val_loader, _, _ = build_dataloaders(config)
        batch = next(iter(train_loader))
        assert batch["input_ids"].shape == (4, 16)
        assert batch["reference_ids"].shape == (4, 16)
        assert batch["labels"].shape == (4, 16)
        assert len(val_loader) == 2


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class TestModel:
    """Shape contracts, bidirectionality and input validation."""

    @pytest.fixture()
    def small_config(self) -> ModelConfig:
        return ModelConfig(d_model=16, num_layers=2, dropout=0.0)

    @staticmethod
    def _inputs(batch: int, length: int) -> tuple:
        """Build a full set of four input channels."""
        return (
            torch.randint(0, 6, (batch, length)),
            torch.randint(0, 6, (batch, length)),
            torch.rand(batch, length) * 60.0,
            torch.rand(batch, length) * 50.0,
        )

    def test_forward_returns_per_token_logits(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config)
        assert model(*self._inputs(3, 20)).shape == (3, 20, 4)

    def test_quality_channel_changes_the_output(self, small_config: ModelConfig) -> None:
        """If quality were ignored, noise and SNPs stay inseparable."""
        model = VariantCaller(small_config).eval()
        observed, reference, _, depth = self._inputs(1, 16)
        low = torch.full((1, 16), 5.0)
        high = torch.full((1, 16), 40.0)
        with torch.no_grad():
            assert not torch.allclose(
                model(observed, reference, low, depth),
                model(observed, reference, high, depth),
            )

    def test_depth_channel_changes_the_output(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config).eval()
        observed, reference, quality, _ = self._inputs(1, 16)
        with torch.no_grad():
            assert not torch.allclose(
                model(observed, reference, quality, torch.full((1, 16), 2.0)),
                model(observed, reference, quality, torch.full((1, 16), 90.0)),
            )

    def test_disabled_channels_are_not_constructed(self) -> None:
        """Ablation must remove the parameters, not merely ignore them."""
        model = VariantCaller(ModelConfig(use_base_quality=False, use_depth=False))
        assert model.input_fusion.quality_embedding is None
        assert model.input_fusion.depth_proj is None

        full = VariantCaller(ModelConfig(use_base_quality=True, use_depth=True))
        assert full.count_parameters() > model.count_parameters()

    def test_missing_required_channel_raises(self) -> None:
        model = VariantCaller(ModelConfig(d_model=16, num_layers=1, use_base_quality=True))
        with pytest.raises(ValueError, match="base_quality is required"):
            model(torch.randint(0, 6, (1, 8)), torch.randint(0, 6, (1, 8)))

    def test_auxiliary_channel_shape_mismatch_raises(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config)
        with pytest.raises(ValueError, match="base_quality"):
            model(
                torch.randint(0, 6, (1, 8)),
                torch.randint(0, 6, (1, 8)),
                torch.rand(1, 12),
                torch.rand(1, 8),
            )

    @pytest.mark.parametrize("fusion", ["concat", "sum"])
    def test_both_fusion_modes_produce_same_shape(self, fusion: str) -> None:
        model = VariantCaller(ModelConfig(d_model=16, num_layers=2, fusion=fusion))  # type: ignore[arg-type]
        assert model(*self._inputs(2, 12)).shape == (2, 12, 4)

    def test_sum_fusion_has_no_projection_layer(self) -> None:
        assert VariantCaller(ModelConfig(fusion="sum")).input_fusion.fusion_proj is None
        assert VariantCaller(ModelConfig(fusion="concat")).input_fusion.fusion_proj is not None

    def test_reference_channel_changes_the_output(self, small_config: ModelConfig) -> None:
        """If the reference were ignored the task would be unsolvable."""
        model = VariantCaller(small_config).eval()
        observed, _, quality, depth = self._inputs(1, 16)
        with torch.no_grad():
            first = model(observed, observed, quality, depth)
            second = model(observed, torch.roll(observed, shifts=3, dims=1), quality, depth)
        assert not torch.allclose(first, second)

    def test_mismatched_shapes_raise(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config)
        with pytest.raises(ValueError, match="identical shapes"):
            model(torch.randint(0, 6, (2, 10)), torch.randint(0, 6, (2, 12)))

    def test_non_2d_input_raises(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config)
        with pytest.raises(ValueError, match="2-D"):
            model(torch.randint(0, 6, (10,)), torch.randint(0, 6, (10,)))

    def test_block_is_genuinely_bidirectional(self) -> None:
        """A later token must influence an earlier token's representation.

        The perturbation must vary across features: every core begins with a
        LayerNorm over the feature dimension, which subtracts the per-position
        mean and therefore erases a uniform shift entirely.
        """
        torch.manual_seed(0)
        block = BiMambaBlock(d_model=8, conv_kernel_size=3, dropout=0.0).eval()
        base = torch.randn(1, 10, 8)

        perturbed = base.clone()
        perturbed[0, 9, :] = torch.linspace(-4.0, 4.0, 8)  # final position only

        with torch.no_grad():
            baseline = block(base)
            changed = block(perturbed)

        # Position 8 sits to the left of the perturbation, so a causal-only
        # block would leave it untouched.
        assert not torch.allclose(baseline[0, 8], changed[0, 8], atol=1e-5)

        # Position 0 is beyond the kernel's reach in either direction, which
        # confirms the effect above came from the backward pass and not from
        # some global mixing.
        assert torch.allclose(baseline[0, 0], changed[0, 0], atol=1e-6)

    def test_gradients_reach_both_embedding_tables(self, small_config: ModelConfig) -> None:
        model = VariantCaller(small_config)
        model(*self._inputs(2, 12)).sum().backward()
        fusion = model.input_fusion
        assert fusion.observed_embedding.weight.grad is not None
        assert fusion.reference_embedding.weight.grad is not None
        assert fusion.quality_embedding is not None
        assert fusion.quality_embedding.weight.grad is not None
        assert fusion.depth_proj is not None
        assert fusion.depth_proj.weight.grad is not None

    def test_parameter_count_is_positive(self, small_config: ModelConfig) -> None:
        assert VariantCaller(small_config).count_parameters() > 0


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------


class TestFocalLoss:
    """Focal modulation, alpha weighting and reduction behaviour."""

    def test_gamma_zero_matches_cross_entropy(self) -> None:
        logits = torch.randn(32, 4)
        targets = torch.randint(0, 4, (32,))
        focal = FocalLoss(gamma=0.0)(logits, targets)
        expected = torch.nn.functional.cross_entropy(logits, targets)
        assert torch.allclose(focal, expected, atol=1e-6)

    def test_confident_predictions_are_down_weighted(self) -> None:
        """The whole point of gamma: easy tokens contribute less."""
        confident = torch.tensor([[10.0, 0.0, 0.0, 0.0]])
        target = torch.tensor([0])
        focal = FocalLoss(gamma=2.0, reduction="none")(confident, target)
        plain = FocalLoss(gamma=0.0, reduction="none")(confident, target)
        assert float(focal.item()) < float(plain.item())

    def test_alpha_scales_the_targeted_class(self) -> None:
        logits = torch.randn(16, 4)
        targets = torch.zeros(16, dtype=torch.long)
        unweighted = FocalLoss(gamma=0.0, reduction="none")(logits, targets)
        weighted = FocalLoss(
            alpha=torch.tensor([2.0, 1.0, 1.0, 1.0]), gamma=0.0, reduction="none"
        )(logits, targets)
        assert torch.allclose(weighted, 2.0 * unweighted, atol=1e-6)

    def test_loss_is_non_negative(self) -> None:
        loss = FocalLoss(gamma=2.0)(torch.randn(64, 4), torch.randint(0, 4, (64,)))
        assert float(loss.item()) >= 0.0

    @pytest.mark.parametrize("reduction,expected", [("none", (16,)), ("mean", ()), ("sum", ())])
    def test_reduction_shapes(self, reduction: str, expected: tuple) -> None:
        loss = FocalLoss(gamma=1.0, reduction=reduction)(
            torch.randn(16, 4), torch.randint(0, 4, (16,))
        )
        assert loss.shape == expected

    def test_label_smoothing_raises_loss_on_perfect_prediction(self) -> None:
        logits = torch.tensor([[20.0, 0.0, 0.0, 0.0]])
        target = torch.tensor([0])
        without = FocalLoss(gamma=0.0, label_smoothing=0.0)(logits, target)
        with_smoothing = FocalLoss(gamma=0.0, label_smoothing=0.1)(logits, target)
        assert float(with_smoothing.item()) > float(without.item())

    def test_rejects_invalid_arguments(self) -> None:
        with pytest.raises(ValueError, match="gamma"):
            FocalLoss(gamma=-1.0)
        with pytest.raises(ValueError, match="reduction"):
            FocalLoss(reduction="median")

    def test_rejects_shape_mismatch(self) -> None:
        with pytest.raises(ValueError, match="token count"):
            FocalLoss()(torch.randn(8, 4), torch.randint(0, 4, (5,)))


class TestClassWeights:
    """Weight calibration must correct imbalance without over-correcting."""

    def test_counts_sum_to_total_tokens(self) -> None:
        dataset = SyntheticVariantDataset(10, seq_len=32, mutation_rate=0.1, seed=11)
        counts = estimate_class_counts(dataset, num_classes=4, sample_size=10)
        assert counts.shape == (4,)
        assert float(counts.sum().item()) == pytest.approx(10 * 32)

    def test_rare_classes_outweigh_the_dominant_one(self) -> None:
        dataset = SyntheticVariantDataset(50, seq_len=64, mutation_rate=0.05, seed=12)
        weights = compute_calibrated_class_weights(dataset, 4, sample_size=50)
        assert float(weights[0]) < 1.0, "Normal is dominant and must be down-weighted"
        assert all(float(weights[label]) > float(weights[0]) for label in (1, 2, 3))

    def test_cap_is_respected(self) -> None:
        dataset = SyntheticVariantDataset(50, seq_len=64, mutation_rate=0.01, seed=13)
        weights = compute_calibrated_class_weights(dataset, 4, sample_size=50, max_weight_cap=2.0)
        assert float(weights.max().item()) <= 2.0

    def test_sqrt_smoothing_beats_raw_inverse_frequency(self) -> None:
        """Guards the specific regression that caused mass false positives."""
        dataset = SyntheticVariantDataset(50, seq_len=64, mutation_rate=0.05, seed=14)
        counts = estimate_class_counts(dataset, 4, sample_size=50)
        raw_ratio = float((counts.sum() / (4 * counts)).max() / (counts.sum() / (4 * counts)).min())

        weights = compute_calibrated_class_weights(dataset, 4, sample_size=50)
        calibrated_ratio = float(weights.max() / weights.min())
        assert calibrated_ratio < raw_ratio


# ---------------------------------------------------------------------------
# Metrics, early stopping, scheduling
# ---------------------------------------------------------------------------


class TestMetrics:
    """F1 computation and the distinction between the two macro averages."""

    def test_perfect_prediction_scores_one(self) -> None:
        targets = torch.tensor([0, 1, 2, 3, 0, 1])
        metrics = compute_metrics(targets.clone(), targets)
        assert metrics["macro_f1"] == pytest.approx(1.0)
        assert metrics["f1_mutations_macro"] == pytest.approx(1.0)

    def test_all_normal_collapse_is_caught_by_mutation_f1(self) -> None:
        """macro_f1 stays misleadingly high; the mutation average does not."""
        targets = torch.tensor([0] * 95 + [1, 1, 2, 2, 3])
        predictions = torch.zeros_like(targets)
        metrics = compute_metrics(predictions, targets)
        assert metrics["f1_mutations_macro"] == pytest.approx(0.0)
        assert metrics["f1_class_0"] > 0.9

    def test_precision_and_recall_are_computed_per_class(self) -> None:
        targets = torch.tensor([0, 0, 1, 1])
        predictions = torch.tensor([0, 1, 1, 1])
        metrics = compute_metrics(predictions, targets)
        assert metrics["precision_class_1"] == pytest.approx(2 / 3)
        assert metrics["recall_class_1"] == pytest.approx(1.0)

    def test_absent_class_scores_zero_without_dividing_by_zero(self) -> None:
        targets = torch.tensor([0, 0, 0])
        metrics = compute_metrics(targets.clone(), targets)
        assert metrics["f1_class_3"] == pytest.approx(0.0)

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="differ in length"):
            compute_metrics(torch.tensor([0, 1]), torch.tensor([0, 1, 2]))


class TestEarlyStopping:
    """Stopping logic and best-checkpoint retention."""

    @pytest.fixture()
    def model(self) -> torch.nn.Module:
        return torch.nn.Linear(4, 4)

    def test_stops_after_patience_without_improvement(self, model: torch.nn.Module) -> None:
        stopper = EarlyStopping(patience=2, mode="min")
        stopper.step(1.0, model, 1)
        assert not stopper.should_stop
        stopper.step(1.5, model, 2)
        assert not stopper.should_stop
        stopper.step(1.4, model, 3)
        assert stopper.should_stop

    def test_improvement_resets_the_counter(self, model: torch.nn.Module) -> None:
        stopper = EarlyStopping(patience=2, mode="min")
        stopper.step(1.0, model, 1)
        stopper.step(1.2, model, 2)
        assert stopper.counter == 1
        stopper.step(0.5, model, 3)
        assert stopper.counter == 0
        assert stopper.best_epoch == 3

    def test_max_mode_tracks_increasing_scores(self, model: torch.nn.Module) -> None:
        stopper = EarlyStopping(patience=2, mode="max")
        assert stopper.step(0.5, model, 1)
        assert stopper.step(0.7, model, 2)
        assert not stopper.step(0.6, model, 3)

    def test_restores_best_weights(self) -> None:
        model = torch.nn.Linear(2, 2)
        stopper = EarlyStopping(patience=3, mode="min")

        with torch.no_grad():
            model.weight.fill_(1.0)
        stopper.step(0.1, model, 1)

        with torch.no_grad():
            model.weight.fill_(9.0)
        stopper.step(5.0, model, 2)

        stopper.restore_best(model)
        assert torch.allclose(model.weight, torch.ones_like(model.weight))

    def test_rejects_invalid_arguments(self) -> None:
        with pytest.raises(ValueError, match="patience"):
            EarlyStopping(patience=0)
        with pytest.raises(ValueError, match="mode"):
            EarlyStopping(mode="maximise")  # type: ignore[arg-type]


class TestScheduler:
    """Warmup then cosine annealing, stepped once per epoch."""

    def test_learning_rate_rises_during_warmup(self) -> None:
        model = torch.nn.Linear(4, 4)
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
        scheduler = build_scheduler(optimizer, warmup_epochs=3, total_epochs=10)

        first = optimizer.param_groups[0]["lr"]
        scheduler.step()
        second = optimizer.param_groups[0]["lr"]
        assert second > first
        assert first == pytest.approx(1e-4, rel=1e-3)

    def test_zero_warmup_is_pure_cosine(self) -> None:
        optimizer = torch.optim.AdamW(torch.nn.Linear(4, 4).parameters(), lr=1e-3)
        scheduler = build_scheduler(optimizer, warmup_epochs=0, total_epochs=10)
        assert optimizer.param_groups[0]["lr"] == pytest.approx(1e-3)
        scheduler.step()
        assert optimizer.param_groups[0]["lr"] < 1e-3

    def test_rejects_warmup_longer_than_training(self) -> None:
        optimizer = torch.optim.AdamW(torch.nn.Linear(4, 4).parameters(), lr=1e-3)
        with pytest.raises(ValueError, match="warmup_epochs"):
            build_scheduler(optimizer, warmup_epochs=10, total_epochs=5)


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------


class TestTrainingIntegration:
    """End-to-end smoke tests on a deliberately tiny configuration."""

    @pytest.fixture()
    def tiny_config(self) -> ExperimentConfig:
        return ExperimentConfig(
            data=DataConfig(train_samples=64, val_samples=32, seq_len=32, batch_size=8),
            model=ModelConfig(d_model=32, num_layers=2, dropout=0.0),
            training=TrainingConfig(num_epochs=4, warmup_epochs=1, weight_sample_size=32),
        )

    def test_fit_records_history_and_reduces_loss(self, tiny_config: ExperimentConfig) -> None:
        """Loss must fall over training.

        Measured across a handful of epochs rather than two: on the noisy
        generator a single epoch-to-epoch step is not reliably monotonic.
        """
        trainer = Trainer(tiny_config, verbose=False)
        history = trainer.fit()
        assert len(history) == tiny_config.training.num_epochs
        assert history[-1].train_loss < history[0].train_loss

    def test_evaluate_returns_full_metric_set(self, tiny_config: ExperimentConfig) -> None:
        trainer = Trainer(tiny_config, verbose=False)
        trainer.fit()
        val_loss, metrics = trainer.evaluate()
        assert val_loss >= 0.0
        assert "f1_mutations_macro" in metrics
        for class_id in range(4):
            assert f"f1_class_{class_id}" in metrics

    def test_validate_is_deterministic_in_eval_mode(self, tiny_config: ExperimentConfig) -> None:
        """Dropout must be inactive, or val_loss would be noise."""
        trainer = Trainer(tiny_config, verbose=False)
        first_loss, first_metrics = trainer.evaluate()
        second_loss, second_metrics = trainer.evaluate()
        assert first_loss == pytest.approx(second_loss)
        assert first_metrics["macro_f1"] == pytest.approx(second_metrics["macro_f1"])

    def test_epoch_report_renders(self, tiny_config: ExperimentConfig) -> None:
        trainer = Trainer(tiny_config, verbose=False)
        line = trainer.fit()[0].format_line()
        assert "macro F1" in line
        assert "F1[Insertion]" in line
