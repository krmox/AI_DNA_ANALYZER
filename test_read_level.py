"""Tests for the read-level representation and model.

The properties that have to hold before any GPU time is spent:

1. the aggregate counts fall out of the read-level tensor *exactly* -- if they
   do not, the two representations describe different reads and no comparison
   between them is meaningful;
2. extraction and read ordering are deterministic and independent of BAM
   traversal order, so the network cannot learn traversal order;
3. the padding mask is respected -- padded rows must not influence the locus
   embedding, or "depth" leaks in through the aggregation;
4. the aggregation is permutation-invariant, so read order carries no signal;
5. no feature encodes reference-discordance or any truth-derived quantity.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from read_level_pileup import (
    FEATURE_DIM,
    FEATURE_NAMES,
    FLAG_COUNTED,
    FLAG_GAP,
    FLAG_INSERTION,
    FLAG_NEAR_END,
    FLAG_VALID,
    MAX_READS,
    READ_DIM,
    build_features,
    read_mask,
    read_sort_key,
    reconstruct_counts,
)

SEQ_LEN = 8


def make_locus(rows: list[list[int]], max_reads: int = MAX_READS) -> np.ndarray:
    """Pad a list of read rows into one ``[R, READ_DIM]`` locus matrix."""
    matrix = np.zeros((max_reads, READ_DIM), dtype=np.uint8)
    matrix[:, 0] = 255
    for index, row in enumerate(rows):
        matrix[index] = row
    return matrix


def row(code: int, bq: int = 35, mq: int = 60, strand: int = 0, position: int = 128,
        end_distance: int = 200, flags: int = FLAG_VALID, mismatch: int = 0) -> list[int]:
    """One compact read row with sensible defaults."""
    return [code, bq, mq, strand, position, end_distance, flags, mismatch]


@pytest.fixture
def synthetic() -> np.ndarray:
    """Three loci exercising counted bases, low quality, gaps and padding."""
    counted = FLAG_VALID | FLAG_COUNTED
    locus_a = make_locus([
        row(0, bq=35, mq=60, flags=counted),
        row(0, bq=30, mq=55, flags=counted),
        row(2, bq=40, mq=60, flags=counted),
        row(1, bq=5, mq=60, flags=FLAG_VALID),                    # fails BQ filter
        row(4, bq=0, mq=45, flags=FLAG_VALID | FLAG_GAP),         # deletion
        row(3, bq=38, mq=60, flags=counted | FLAG_INSERTION),     # insertion after
    ])
    locus_b = make_locus([row(1, bq=20, mq=25, flags=counted) for _ in range(5)])
    locus_c = make_locus([])                                       # zero depth
    return np.stack([locus_a, locus_b, locus_c])


# --- 1. exact reconstruction of the aggregate representation -----------------


def test_reconstruction_matches_hand_computed_counts(synthetic):
    counts = reconstruct_counts(synthetic, np.array([1, 2, 1]))
    # Locus A: A=2, C=0 (the C failed the quality filter), G=1, T=1.
    assert counts[0, 0] == 2 and counts[0, 1] == 0
    assert counts[0, 2] == 1 and counts[0, 3] == 1
    assert counts[0, 4] == 1                       # one gap read
    assert counts[0, 5] == 1                       # one insertion-carrying read
    assert counts[0, 6] == 6                       # depth counts every admitted read
    assert counts[0, 7] == 35 + 30 + 40 + 38       # quality_sum over counted bases only
    assert counts[0, 8] == 60 + 55 + 60 + 60 + 45 + 60   # mapping_sum over all reads
    assert counts[1, 1] == 5 and counts[1, 6] == 5
    assert counts[2].tolist()[:9] == [0] * 9       # empty locus contributes nothing


def test_low_quality_bases_excluded_from_counts_but_kept_as_rows(synthetic):
    """The read the binomial caller never sees is still present in the tensor."""
    counts = reconstruct_counts(synthetic, np.array([1, 2, 1]))
    assert counts[0, 1] == 0                       # not counted
    assert counts[0, 6] == 6                       # but included in depth
    features = build_features(synthetic)
    assert features[0, 3, FEATURE_NAMES.index("low_quality_base")] == 1.0


def test_reference_index_is_passed_through_untouched(synthetic):
    counts = reconstruct_counts(synthetic, np.array([1, 2, 3]))
    assert counts[:, 9].tolist() == [1.0, 2.0, 3.0]


def test_padding_contributes_nothing_to_counts():
    empty = make_locus([])[None]
    counts = reconstruct_counts(empty, np.array([1]))
    assert counts[0, :9].sum() == 0


# --- 2. determinism and ordering --------------------------------------------


def test_read_sort_key_is_stable_across_processes():
    """A salted hash would make extraction irreproducible between runs."""
    assert read_sort_key("read_001") == read_sort_key("read_001")
    assert read_sort_key("read_001") != read_sort_key("read_002")
    # Value pinned so a change in hashing is caught rather than silently
    # reshuffling every cached tensor.
    assert read_sort_key("read_001") == 15545452545807633358


def test_read_sort_key_is_uncorrelated_with_read_index():
    keys = [read_sort_key(f"read_{i:05d}") for i in range(500)]
    ranks = np.argsort(np.argsort(keys))
    correlation = np.corrcoef(ranks, np.arange(500))[0, 1]
    assert abs(correlation) < 0.15


# --- 3. feature construction and masking ------------------------------------


def test_features_have_expected_shape_and_no_nans(synthetic):
    features = build_features(synthetic)
    assert features.shape == (3, MAX_READS, FEATURE_DIM)
    assert np.all(np.isfinite(features))


def test_padded_rows_are_all_zero(synthetic):
    features = build_features(synthetic)
    mask = read_mask(synthetic)
    assert np.all(features[~mask] == 0.0)
    assert mask[0].sum() == 6 and mask[1].sum() == 5 and mask[2].sum() == 0


def test_base_one_hot_is_exclusive_and_only_for_real_bases(synthetic):
    features = build_features(synthetic)
    mask = read_mask(synthetic)
    one_hot = features[..., 0:5]
    assert np.all(one_hot.sum(-1)[mask] <= 1.0)
    assert np.all(np.isin(one_hot, (0.0, 1.0)))


def test_feature_ranges_are_bounded(synthetic):
    features = build_features(synthetic)
    assert features.min() >= 0.0
    assert features.max() <= 4.3            # base quality 255/60 is the ceiling


def test_ablation_zeroes_only_its_group(synthetic):
    full = build_features(synthetic)
    dropped = build_features(synthetic, drop_groups=("strand",))
    strand = FEATURE_NAMES.index("strand")
    assert np.all(dropped[..., strand] == 0.0)
    others = [i for i in range(FEATURE_DIM) if i != strand]
    assert np.array_equal(full[..., others], dropped[..., others])


def test_ablation_preserves_tensor_shape(synthetic):
    """Ablations must not change capacity, only information."""
    for group in ("quality", "strand", "position", "indel"):
        assert build_features(synthetic, drop_groups=(group,)).shape == \
            build_features(synthetic).shape


# --- 4. no forbidden information --------------------------------------------


def test_no_feature_encodes_reference_discordance():
    """Two loci with identical reads must give identical features regardless of
    what the reference base is -- the extractor never sees it."""
    reads = make_locus([row(0, flags=FLAG_VALID | FLAG_COUNTED)])[None]
    a = reconstruct_counts(reads, np.array([1]))       # reference A
    b = reconstruct_counts(reads, np.array([2]))       # reference C
    assert np.array_equal(a[:, :9], b[:, :9])
    assert np.array_equal(build_features(reads), build_features(reads))


def test_feature_names_contain_no_truth_derived_quantity():
    forbidden = ("label", "truth", "vcf", "alt", "variant", "snp", "genotype",
                 "is_ref", "differs")
    for name in FEATURE_NAMES:
        assert not any(token in name.lower() for token in forbidden)


# --- 5. model: masking, permutation invariance, parameter count -------------


def test_read_encoder_is_permutation_invariant():
    from model_read_level import ReadLevelMamba
    torch.manual_seed(0)
    model = ReadLevelMamba().eval()
    features = torch.randn(2, SEQ_LEN, MAX_READS, FEATURE_DIM)
    mask = torch.zeros(2, SEQ_LEN, MAX_READS, dtype=torch.bool)
    mask[..., :6] = True

    reference = torch.eye(5)[torch.randint(0, 5, (2, SEQ_LEN))]
    with torch.no_grad():
        base = model(read_features=features, read_valid=mask, reference_onehot=reference)
        permutation = torch.randperm(MAX_READS)
        permuted = model(read_features=features[:, :, permutation],
                         read_valid=mask[:, :, permutation], reference_onehot=reference)
    assert torch.allclose(base, permuted, atol=1e-5)


def test_padded_reads_do_not_affect_the_locus_embedding():
    """Changing the contents of masked rows must change nothing."""
    from model_read_level import ReadLevelMamba
    torch.manual_seed(0)
    model = ReadLevelMamba().eval()
    features = torch.randn(2, SEQ_LEN, MAX_READS, FEATURE_DIM)
    mask = torch.zeros(2, SEQ_LEN, MAX_READS, dtype=torch.bool)
    mask[..., :4] = True

    reference = torch.eye(5)[torch.randint(0, 5, (2, SEQ_LEN))]
    polluted = features.clone()
    polluted[:, :, 4:] = 999.0
    with torch.no_grad():
        assert torch.allclose(
            model(read_features=features, read_valid=mask, reference_onehot=reference),
            model(read_features=polluted, read_valid=mask, reference_onehot=reference),
            atol=1e-6)


def test_zero_depth_locus_produces_finite_output():
    """An all-padding locus must not divide by zero."""
    from model_read_level import ReadLevelMamba
    torch.manual_seed(0)
    model = ReadLevelMamba().eval()
    features = torch.zeros(1, SEQ_LEN, MAX_READS, FEATURE_DIM)
    mask = torch.zeros(1, SEQ_LEN, MAX_READS, dtype=torch.bool)
    reference = torch.eye(5)[torch.randint(0, 5, (1, SEQ_LEN))]
    with torch.no_grad():
        out = model(read_features=features, read_valid=mask, reference_onehot=reference)
    assert torch.all(torch.isfinite(out))


def test_parameter_count_is_comparable_to_the_aggregate_model():
    from model_raw_pileup import RawPileupMamba
    from model_read_level import ReadLevelMamba
    aggregate = sum(p.numel() for p in RawPileupMamba().parameters())
    read_level = sum(p.numel() for p in ReadLevelMamba().parameters())
    assert read_level < 1.25 * aggregate, (read_level, aggregate)


def test_output_shape_is_per_locus_four_class():
    from model_read_level import ReadLevelMamba
    model = ReadLevelMamba().eval()
    features = torch.randn(3, SEQ_LEN, MAX_READS, FEATURE_DIM)
    mask = torch.ones(3, SEQ_LEN, MAX_READS, dtype=torch.bool)
    reference = torch.eye(5)[torch.randint(0, 5, (3, SEQ_LEN))]
    with torch.no_grad():
        out = model(read_features=features, read_valid=mask, reference_onehot=reference)
    assert out.shape == (3, SEQ_LEN, 4)


def test_torch_feature_builder_matches_numpy_exactly(synthetic):
    """The GPU path and the validated numpy path must not diverge."""
    from read_level_pileup import build_features_torch, resolve_drop_indices
    reference = build_features(synthetic)
    produced = build_features_torch(torch.tensor(synthetic)).numpy()
    assert np.array_equal(reference, produced)

    for group in ("quality", "strand", "position", "indel"):
        reference = build_features(synthetic, drop_groups=(group,))
        produced = build_features_torch(torch.tensor(synthetic),
                                        resolve_drop_indices((group,))).numpy()
        assert np.array_equal(reference, produced), group


def test_torch_feature_builder_matches_numpy_on_real_cache():
    """Same check against real extracted data, not just synthetic rows."""
    from read_level_pileup import build_features_torch
    blob = np.load("data/giab_hg002_15x_readlevel/test_15x_reads.npz", allow_pickle=True)
    sample = blob["reads"][:2000]
    assert np.array_equal(build_features(sample),
                          build_features_torch(torch.tensor(sample)).numpy())


# --- 6. controls preserve exactly what they claim to preserve ---------------


def test_decouple_preserves_aggregate_counts_exactly(synthetic):
    """CONTROL 6 must leave every quantity the binomial caller uses untouched."""
    from read_level_controls import decouple_base_from_attributes
    reference = np.array([1, 2, 1])
    transformed = decouple_base_from_attributes(synthetic, seed=7)
    assert np.array_equal(reconstruct_counts(synthetic, reference),
                          reconstruct_counts(transformed, reference))


def test_decouple_preserves_attribute_marginals(synthetic):
    from read_level_controls import decouple_base_from_attributes
    transformed = decouple_base_from_attributes(synthetic, seed=7)
    mask = read_mask(synthetic)
    for column in (1, 2, 3, 4, 5):
        assert sorted(synthetic[0, :, column][mask[0]]) == \
            sorted(transformed[0, :, column][mask[0]])


def test_decouple_actually_changes_the_association():
    """A locus with a clean quality/base association must be disrupted."""
    from read_level_controls import decouple_base_from_attributes
    counted = FLAG_VALID | FLAG_COUNTED
    locus = make_locus([row(0, bq=40, flags=counted) for _ in range(8)] +
                       [row(1, bq=14, flags=counted) for _ in range(8)])[None]
    transformed = decouple_base_from_attributes(locus, seed=3)
    high = transformed[0, :16][transformed[0, :16, 0] == 0][:, 1]
    assert not np.all(high == 40)


def test_permute_reads_preserves_counts_and_multiset(synthetic):
    from read_level_controls import permute_reads
    reference = np.array([1, 2, 1])
    transformed = permute_reads(synthetic, seed=5)
    assert np.array_equal(reconstruct_counts(synthetic, reference),
                          reconstruct_counts(transformed, reference))
    for locus in range(synthetic.shape[0]):
        assert sorted(map(tuple, synthetic[locus].tolist())) == \
            sorted(map(tuple, transformed[locus].tolist()))


def test_aggregate_only_preserves_counts_and_removes_variation(synthetic):
    from read_level_controls import aggregate_only
    transformed = aggregate_only(synthetic)
    reference = np.array([1, 2, 1])
    original_counts = reconstruct_counts(synthetic, reference)
    new_counts = reconstruct_counts(transformed, reference)
    # Allele counts, depth, gaps and insertions must be identical.
    assert np.array_equal(original_counts[:, [0, 1, 2, 3, 4, 5, 6]],
                          new_counts[:, [0, 1, 2, 3, 4, 5, 6]])
    # Per-read variation in quality/strand/position is gone.
    mask = read_mask(transformed)
    for column in (1, 2, 3, 4, 5):
        values = transformed[0, :, column][mask[0]]
        assert len(set(values.tolist())) == 1, column


def test_aggregate_only_preserves_mean_quality_to_rounding(synthetic):
    from read_level_controls import aggregate_only
    mask = read_mask(synthetic)
    before = synthetic[0, :, 2][mask[0]].mean()
    after = aggregate_only(synthetic)[0, :, 2][mask[0]].mean()
    assert abs(before - after) <= 0.5


def test_controls_are_deterministic(synthetic):
    from read_level_controls import decouple_base_from_attributes, permute_reads
    assert np.array_equal(decouple_base_from_attributes(synthetic, seed=11),
                          decouple_base_from_attributes(synthetic, seed=11))
    assert np.array_equal(permute_reads(synthetic, seed=11),
                          permute_reads(synthetic, seed=11))


def test_reference_base_changes_the_prediction():
    """The locus reference must actually reach the output: identical reads with
    a different reference base are a different biological situation."""
    from model_read_level import ReadLevelMamba
    torch.manual_seed(0)
    model = ReadLevelMamba().eval()
    features = torch.randn(1, SEQ_LEN, MAX_READS, FEATURE_DIM)
    mask = torch.ones(1, SEQ_LEN, MAX_READS, dtype=torch.bool)
    reference_a = torch.eye(5)[torch.full((1, SEQ_LEN), 1)]
    reference_c = torch.eye(5)[torch.full((1, SEQ_LEN), 2)]
    with torch.no_grad():
        a = model(read_features=features, read_valid=mask, reference_onehot=reference_a)
        c = model(read_features=features, read_valid=mask, reference_onehot=reference_c)
    assert not torch.allclose(a, c, atol=1e-4)
