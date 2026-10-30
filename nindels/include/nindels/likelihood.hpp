// Spec §8: local read-level verification.
//   P(R_r | H) from a 3-state PairHMM (base qualities, indel errors, alignment
//   uncertainty), H_ref / H_alt construction, the original
//   Delta_ell = ell(H_alt) - ell(H_ref), and the diploid genotype likelihoods
//   GL(0/0), GL(0/1), GL(1/1) required by the author's correction.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "nindels/evidence.hpp"

namespace nindels {

struct PairHmmParams {
    int gap_open_q = 45;  // GATK default, not re-fitted
    int gap_ext_q = 10;   // GATK default, not re-fitted
};

// ln P(read | haplotype). The read must align end to end; it may start
// anywhere on the haplotype (free start, as in GATK's PairHMM).
// With band >= 0, only columns within +-band of the diagonal "read base i at
// hap column i + diag" are computed (falls back to the full matrix if the band
// leaves the haplotype).
double pairhmm_loglik(std::string_view read, const uint8_t* quals, std::string_view hap,
                      const PairHmmParams& p = {}, int64_t diag = 0, int64_t band = -1);

// A read prepared for likelihood computation: aligned part only (soft clips
// removed), with its reference span and the allele its CIGAR implies.
struct HmmRead {
    std::string bases;
    std::vector<uint8_t> quals;
    int64_t ref_start = 0, ref_end = 0;  // [start, end) on the reference
    uint8_t mapq = 0;
    ReadAllele cigar_allele = ReadAllele::Ref;
    uint64_t name_hash = 0;
    std::vector<int16_t> window_bq;  // base quality at event-centred positions, -1 if not aligned
};

// H_ref and H_alt over the reference window [start, end).
struct Haplotypes {
    std::string ref, alt;
    int64_t start = 0;
    Indel indel;
    // Index in H_alt of reference position g (g outside a deletion).
    int64_t alt_index(int64_t g) const;
};

Haplotypes build_haplotypes(const RefSeq& ref, const Indel& x, int64_t start, int64_t end);

struct ReadLik {
    double lp_ref = 0, lp_alt = 0;  // ln P(R|H_ref), ln P(R|H_alt), uncapped
    uint8_t mapq = 0;
    ReadAllele cigar_allele = ReadAllele::Ref;
    bool best_alt() const { return lp_alt > lp_ref; }
};

struct LikelihoodParams {
    PairHmmParams hmm;
    int max_reads = 250;   // ARBITRARY cost cap; reads chosen by smallest name hash
    int hap_pad = 20;      // ARBITRARY: haplotype bases beyond each read end
    int band = 20;         // ARBITRARY: DP band half-width is band + indel length (-1 = full DP)
};

// P(R|H_ref) and P(R|H_alt) for up to max_reads reads.
std::vector<ReadLik> read_likelihoods(const RefSeq& ref, const Indel& x, const std::vector<HmmRead>& reads,
                                      const LikelihoodParams& p, std::vector<size_t>* used = nullptr);

// Spec §8 as originally written: sum ln P(R|H_alt) - sum ln P(R|H_ref).
// Kept as an intentionally limited baseline: it tests hom-alt against hom-ref.
double delta_ell_original(const std::vector<ReadLik>& liks);

struct GenotypePriors {
    // ARBITRARY: 0/0 : 0/1 : 1/1 = 1 : 1e-4 : 5e-5 (ln scale below).
    std::array<double, 3> ln_prior = {0.0, -9.210340371976184, -9.903487552536127};
};

struct GenotypeLik {
    std::array<double, 3> gl{};         // ln P(reads | g), g = 0/0, 0/1, 1/1
    std::array<double, 3> posterior{};  // normalised with the priors
    int gt = 0;                         // argmax posterior
    double qual = 0;   // ln[ P(reads|0/1) + P(reads|1/1) ] - ln P(reads|0/0): non-ref vs ref LR
    double gq = 0;     // phred of 1 - posterior[gt], capped at 99
};

// Diploid model: GL(g) = sum_r ln[ (2-g)/2 P(R_r|H_ref) + g/2 P(R_r|H_alt) ].
// With mapq_cap, each read's two likelihoods are mixed with a floor so that a
// read can shift the comparison by at most MAPQ/10 * ln(10).
GenotypeLik diploid_gl(const std::vector<ReadLik>& liks, bool mapq_cap, const GenotypePriors& pr);

// Binomial diploid model on counts only (no HMM): P(alt read | g) =
// err, 0.5, 1-err. Gives every candidate a genotype.
GenotypeLik count_gl(int n_alt, int n_ref, double err, const GenotypePriors& pr);

// Fraction of Ref/Alt reads whose CIGAR allele agrees with their HMM-best haplotype.
Opt hmm_consistency(const std::vector<ReadLik>& liks);

}  // namespace nindels
