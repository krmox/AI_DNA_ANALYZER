// Spec §5: per-feature distances d_k, combined distance D_ij, similarity S_ij.
#pragma once

#include <vector>

#include "nindels/pattern.hpp"

namespace nindels {

struct SimilarityParams {
    // w_k (spec §5.2). Default 1 for every feature: ARBITRARY, the spec says
    // they must be learned; no learning objective is given.
    double w_type = 1, w_length = 1, w_motif = 1, w_context = 1, w_hp = 1, w_repeat = 1;
    // Scale parameters (spec §5.1 s_L; s_h and s_r are the same construction for
    // homopolymer length and repeat copies). Set from training data with
    // estimate_scales().
    double s_L = 1, s_h = 1, s_r = 1;

    void validate() const;  // w_k >= 0, sum(w_k) > 0, s > 0
};

int levenshtein(std::string_view a, std::string_view b);

// LevenshteinDistance(a,b) / max(len a, len b); 0 when both are empty.
double d_sequence(std::string_view a, std::string_view b);

struct FeatureDistances {
    double type, length, motif, context, hp, repeat;  // each in [0,1]
};

FeatureDistances feature_distances(const Pattern& i, const Pattern& j, const SimilarityParams& p);

// D_ij = sum(w_k d_k) / sum(w_k)
double pattern_distance(const Pattern& i, const Pattern& j, const SimilarityParams& p);

// S_ij = 1 - D_ij. Not a probability (spec §16 rule 6).
inline double pattern_similarity(const Pattern& i, const Pattern& j, const SimilarityParams& p) {
    return 1.0 - pattern_distance(i, j, p);
}

// s_L, s_h, s_r := 95th percentile of the raw |log(1+a) - log(1+b)| over
// n_pairs deterministic random pairs of training patterns (pairs with a zero
// difference are excluded). The 95th percentile is an ARBITRARY choice.
SimilarityParams estimate_scales(const std::vector<Pattern>& training, SimilarityParams base,
                                 size_t n_pairs = 200000, uint64_t seed = 1);

}  // namespace nindels
