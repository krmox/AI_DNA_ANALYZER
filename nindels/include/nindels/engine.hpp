// Spec §13: PROCESS_INDEL — the complete per-candidate algorithm.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "nindels/evidence.hpp"
#include "nindels/likelihood.hpp"
#include "nindels/memory.hpp"
#include "nindels/model.hpp"
#include "nindels/pattern.hpp"
#include "nindels/reassess.hpp"

namespace nindels {

enum class Decision { Pass, LowProb, Uncertain, NoModel };
const char* decision_name(Decision d);

struct EngineConfig {
    PatternParams pattern;
    EvidenceParams evidence;
    LikelihoodParams likelihood;
    GenotypePriors priors;
    double count_gl_error = 0.01;  // ARBITRARY indel error rate for the count model
    ComplexityParams complexity;
    ReassessParams reassess;
    NucleotideParams nucleotide;
    bool force_reassess = false;   // every candidate gets step 13 (used to build training tables)
    int min_depth = 4;             // ARBITRARY: fewer informative reads -> UNCERTAIN
    int ood_limit = 2;             // ARBITRARY: >= 2 features outside training range -> UNCERTAIN
    double tau_call = 0.5;         // must be replaced by the value chosen on validation data
};

// Temporary candidate storage (spec §14): everything needed for one indel.
struct CandidateInput {
    Indel indel;
    std::vector<ReadObservation> obs;  // one per informative read
    std::vector<HmmRead> reads;        // Ref/Alt/Other reads for the PairHMM
};

struct CandidateResult {
    Indel indel;
    Pattern pattern;                    // step 2
    Evidence evidence;                  // steps 3-5
    std::optional<ContextStats> context;// steps 6-9 (absent without a memory)
    Complexity complexity;              // step 10
    GenotypeLik count_genotype;         // diploid genotype from read counts
    double p_evidence = 0;              // non-ref posterior of the count model
    double p_base = 0;                  // stage-1 probability used for U
    double U = 0;
    Opt D;                              // |P_context - P_evidence|; missing without memory
    double R = 0, tau = 0;              // steps 11-12
    bool reassessed = false;            // step 13
    int hmm_reads = 0;
    std::optional<GenotypeLik> hmm_genotype;
    Opt delta_ell_original;             // spec §8 as written
    Opt consistency;
    NucleotideResult nucleotide;
    FeatureMap features;                // step 14
    Opt p_final;                        // step 15
    int ood = 0;
    Decision decision = Decision::NoModel;  // step 16
    int gt = 0;                         // 0 = 0/0, 1 = 0/1, 2 = 1/1
};

CandidateResult process_indel(const CandidateInput& in, const RefSeq& ref, const PatternMemory* memory,
                              const LogisticModel* base_model, const LogisticModel* final_model,
                              const EngineConfig& cfg);

// Every feature the engine can produce, in a fixed order. A model selects a
// subset by name; the variants (with/without memory, original vs diploid
// likelihood, Kish vs distinct-locus N_eff) are different subsets.
FeatureMap make_features(const CandidateResult& r);
const std::vector<std::string>& all_feature_names();
// Names of the features produced by step 13 (missing when not reassessed).
const std::vector<std::string>& reassessment_feature_names();

}  // namespace nindels
