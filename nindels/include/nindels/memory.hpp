// Spec §6, §7, §7.1, §14: compact pattern memory, indexed retrieval and the
// accumulated positive/negative statistics.
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "nindels/similarity.hpp"

namespace nindels {

// One independently labelled observation used to build the memory. The label
// must come from an external truth set, never from nindels' own calls (spec §6).
struct LabeledPattern {
    Pattern pattern;
    bool positive = false;
    int32_t tid = -1;
    int64_t pos = 0;
};

struct MemoryParams {
    double theta = 0.15;              // ARBITRARY: leader-clustering radius in D units
    int top_k = 32;                   // ARBITRARY: classes retrieved per query
    int reservoir = 64;               // ARBITRARY: members kept per class for the medoid
    int max_leaders_per_bucket = 4096;// ARBITRARY: bound on classes per index bucket
    int locus_merge_bp = 50;          // ARBITRARY: events in the same 50-bp bin are one locus
};

// M_c = (mu_c, n_c_positive, n_c_negative, metadata_c)
struct PatternClass {
    Pattern medoid;           // mu_c
    uint32_t n_pos = 0;       // n_c_positive
    uint32_t n_neg = 0;       // n_c_negative
    uint32_t loci = 0;        // distinct training loci behind the class (for N_eff_loci)
};

// Result of one memory lookup (spec §7, §13 step 9).
struct ContextStats {
    int retrieved = 0;
    double n_positive = 0;   // N_positive_i = sum S_ic n_c_positive
    double n_negative = 0;   // N_negative_i = sum S_ic n_c_negative
    double p_context = 0;    // (N_pos + alpha) / (N_pos + N_neg + alpha + beta)
    Opt n_eff_kish;          // spec §7.1 formula; missing when nothing was retrieved
    Opt n_eff_loci;          // same formula with distinct loci instead of observations
    Opt s_max;               // highest similarity to a retrieved class
    std::vector<std::pair<uint32_t, double>> closest;  // (class id, S) best first
};

class PatternMemory {
public:
    static PatternMemory build(std::vector<LabeledPattern> items, const SimilarityParams& sim,
                               const MemoryParams& mp, std::string metadata);

    ContextStats query(const Pattern& p) const;

    void save(const std::string& path) const;
    static PatternMemory load(const std::string& path);

    const std::vector<PatternClass>& classes() const { return classes_; }
    double alpha() const { return alpha_; }
    double beta() const { return beta_; }
    const SimilarityParams& similarity_params() const { return sim_; }
    const std::string& metadata() const { return metadata_; }

    static uint32_t bucket_key(IndelType t, int length_bin, int hp_bin, int period);
    static int length_bin(int length);  // {1,2,3-5,6-10,11-20,21-50,>50} -> 0..6
    static int hp_bin(int hp);          // {<=3,4-7,8-12,>12} -> 0..3

private:
    void rebuild_index();
    void estimate_alpha_beta();

    SimilarityParams sim_;
    MemoryParams mp_;
    std::string metadata_;  // metadata_c, shared: platform, reference, training set, ...
    double alpha_ = 1, beta_ = 1;
    std::vector<PatternClass> classes_;
    std::unordered_map<uint32_t, std::vector<uint32_t>> index_;
};

}  // namespace nindels
