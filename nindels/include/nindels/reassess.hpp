// Spec §9 (aggressiveness A, priority R_i, tau_review), §10 (region
// complexity C_r) and §11 (nucleotide-level evidence Q_k).
#pragma once

#include <cmath>
#include <vector>

#include "nindels/likelihood.hpp"
#include "nindels/pattern.hpp"

namespace nindels {

// ------------------------------------------------------------- §10 complexity

struct ComplexityParams {
    // v_1..v_4: equal, ARBITRARY (spec gives no values). Must be >= 0, sum 1.
    double v_hp = 0.25, v_repeat = 0.25, v_lowcomplexity = 0.25, v_alignment = 0.25;
    double hp_norm = 15;    // ARBITRARY: H_r = min(1, hp_length / 15)
    double span_norm = 50;  // ARBITRARY: R_r = min(1, repeat_span / 50)
    void validate() const;
};

struct Complexity {
    double H = 0, R = 0, L = 0;  // normalized homopolymer, repeat, low-complexity (1 - entropy/2)
    Opt G;                        // normalized alignment instability (from reads)
    double C = 0;                 // C_r; if G is missing the other weights are renormalized
};

Complexity region_complexity(const Pattern& p, Opt instability, const ComplexityParams& cp);

// ---------------------------------------------------------- §9 reassessment

struct ReassessParams {
    double a_U = 1.0 / 3, a_C = 1.0 / 3, a_D = 1.0 / 3;  // ARBITRARY equal weights, sum 1
    double tau_strict = 0.6;  // ARBITRARY until set experimentally (spec §9.2)
    double tau_loose = 0.2;   // ARBITRARY until set experimentally (spec §9.2)
    double A = 0.5;           // aggressiveness, 0..1
    void validate() const;
};

// U_i = 1 - |2 P_base - 1|. P_base must be a calibrated probability.
inline double uncertainty(double p_base) { return 1.0 - std::abs(2.0 * p_base - 1.0); }

// tau_review(A) = tau_strict - A (tau_strict - tau_loose)
inline double tau_review(const ReassessParams& p) { return p.tau_strict - p.A * (p.tau_strict - p.tau_loose); }

// R_i = a_U U + a_C C + a_D D
inline double reassess_priority(double U, double C, double D, const ReassessParams& p) {
    return p.a_U * U + p.a_C * C + p.a_D * D;
}

// ------------------------------------------------------ §11 nucleotide level

struct NucleotideParams {
    // w_1..w_4 for B, M, A, S: equal, ARBITRARY. Must be >= 0.
    double w_bq = 0.25, w_mq = 0.25, w_align = 0.25, w_support = 0.25;
    double bq_norm = 40, mq_norm = 60;  // ARBITRARY normalisers to [0,1]
    int half_window = 10;               // positions event-10 .. event+10
};

struct NucleotideResult {
    std::vector<Opt> q;  // Q_k per window position; missing where no read is aligned
    Opt q_min;
    int argmin_offset = 0;  // position of q_min relative to the event
};

// Q_k = w1 B_k + w2 M_k + w3 A_k + w4 S_k, evaluated against the complete
// haplotypes: A_k and S_k use each read's HMM-best haplotype (spec §11 note),
// never an independent per-base decision.
//   B_k: mean base quality at k / bq_norm
//   M_k: mean MAPQ of reads aligned at k / mq_norm
//   A_k: fraction of reads aligned at k whose CIGAR allele agrees with their HMM-best haplotype
//   S_k: 1 - |fraction of reads at k whose HMM-best haplotype is H_alt - f|
// reads[i] pairs with liks[i].
NucleotideResult nucleotide_scores(const std::vector<const HmmRead*>& reads, const std::vector<ReadLik>& liks,
                                   Opt f, const NucleotideParams& p);

}  // namespace nindels
