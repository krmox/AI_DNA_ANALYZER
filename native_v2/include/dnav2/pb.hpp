// Poisson-binomial LLR: frozen V1 semantics of
//   quality_error_model.poisson_binomial_llr(evidence, k, alt_index=alt_index)
// (the variant that reproduces the cached 300x pb_llr bit-for-bit; the min(k, n_counted) variant
//  is NOT the V1 contract).
//
// Two implementations, required to be bit-identical to each other and to the NumPy oracle:
//   pb_reference  - transparent 1:1 port of the Python loops (full 49-wide DP, every slot).
//   pb_native_v2  - same arithmetic, exact work reduction:
//        (a) unused slots are skipped        (V1: `dp = where(active, updated, dp)`)
//        (b) DP truncated to columns 0..k    (column j depends only on columns j and j-1, and
//                                             only column k is read out)
//        (c) column 0 is a running sum       (logaddexp(x, -inf) == x exactly)
//        (d) log(q), log1p(-q) come from a 256-entry table indexed by Phred (q is a pure
//            function of Q, so the table holds the very same doubles V1 computes per read)
// Requires -ffp-contract=off.
#pragma once
#include <cmath>
#include <cstdint>
#include "binomial.hpp"
#include "numerics.hpp"

namespace dnav2 {

constexpr double EPSILON_FLOOR = 1e-4;   // quality_error_model.EPSILON_FLOOR
constexpr double EPSILON_CEILING = 0.25; // quality_error_model.EPSILON_CEILING

// One retained read slot of the V1 read tensor, reduced to what PB reads
struct PbSlot {
    uint8_t base_code;  // 0..3 A,C,G,T (4 gap, 5 other, 255 pad)
    uint8_t bq;         // Phred
    uint8_t counted;    // FLAG_COUNTED (A/C/G/T base with BQ >= 13)
};

// error probability exactly as V1: clip(10 ** (-Q / 10), floor, ceiling)
inline double phred_p(uint8_t bq) {
    return clip(std::pow(10.0, -static_cast<double>(bq) / 10.0), EPSILON_FLOOR, EPSILON_CEILING);
}
// hypothesis success probabilities, expressions in V1's evaluation order
inline double q_h0(double p) { return p / 3.0; }
inline double q_het(double p) { return 0.5 * (1.0 - p) + 0.5 * (p / 3.0); }
inline double q_hom(double p) { return 1.0 - p + p / 3.0; }

inline double pb_llr_from(double h0, double het, double hom, int n_counted) {
    double alt = het > hom ? het : hom;  // np.maximum
    double llr = alt - h0;
    if (n_counted <= 0) llr = 0.0;
    if (!std::isfinite(llr)) llr = 0.0;  // nan_to_num(nan=0, posinf=0, neginf=0)
    return llr;
}

// ------------------------------------------------------------------ reference ------------
inline double log_pb_pmf_reference(const double q_in[MAX_READS], const bool active[MAX_READS], int k) {
    constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
    double dp[MAX_READS + 1], nd[MAX_READS + 1];
    for (int j = 0; j <= MAX_READS; j++) dp[j] = NEG_INF;
    dp[0] = 0.0;
    for (int s = 0; s < MAX_READS; s++) {
        double q = clip(active[s] ? q_in[s] : 0.0, 1e-12, 1 - 1e-12);
        double lq = std::log(q), lnq = std::log1p(-q);
        for (int j = 0; j <= MAX_READS; j++) {
            double stay = dp[j] + lnq;
            double step = j >= 1 ? dp[j - 1] + lq : NEG_INF;
            double upd = logaddexp(stay, step);
            nd[j] = active[s] ? upd : dp[j];
        }
        for (int j = 0; j <= MAX_READS; j++) dp[j] = nd[j];
    }
    return dp[k];
}

// slots: the MAX_READS tensor rows in tensor (hash) order; alt_index 0..3
inline double pb_reference(const PbSlot slots[MAX_READS], int alt_index) {
    double p[MAX_READS];
    bool active[MAX_READS];
    int n_counted = 0, k = 0;
    for (int s = 0; s < MAX_READS; s++) {
        p[s] = phred_p(slots[s].bq);
        active[s] = slots[s].counted != 0;
        if (active[s]) { n_counted++; if (slots[s].base_code == alt_index) k++; }
    }
    double q0[MAX_READS], qh[MAX_READS], qo[MAX_READS];
    for (int s = 0; s < MAX_READS; s++) { q0[s] = q_h0(p[s]); qh[s] = q_het(p[s]); qo[s] = q_hom(p[s]); }
    double h0 = log_pb_pmf_reference(q0, active, k);
    double het = log_pb_pmf_reference(qh, active, k);
    double hom = log_pb_pmf_reference(qo, active, k);
    return pb_llr_from(h0, het, hom, n_counted);
}

// ------------------------------------------------------------------ native v2 ------------
struct PbTables {
    double lq[3][256], lnq[3][256];   // [hypothesis][Phred]
    PbTables() {
        for (int q = 0; q < 256; q++) {
            double p = phred_p(static_cast<uint8_t>(q));
            double qq[3] = {q_h0(p), q_het(p), q_hom(p)};
            for (int h = 0; h < 3; h++) {
                double c = clip(qq[h], 1e-12, 1 - 1e-12);
                lq[h][q] = std::log(c);
                lnq[h][q] = std::log1p(-c);
            }
        }
    }
};
inline const PbTables& pb_tables() { static const PbTables t; return t; }

// Compact locus: only the counted reads, in tensor order.  Returns pb LLR.
inline double pb_native_v2(const uint8_t* code, const uint8_t* bq, int n_counted, int alt_index,
                           const PbTables& T = pb_tables()) {
    int k = 0;
    for (int i = 0; i < n_counted; i++) k += (code[i] == alt_index);
    double res[3];
    for (int h = 0; h < 3; h++) {
        double dp[MAX_READS + 1];
        dp[0] = 0.0;
        for (int j = 1; j <= k; j++) dp[j] = -std::numeric_limits<double>::infinity();
        const double* lq = T.lq[h];
        const double* lnq = T.lnq[h];
        for (int i = 0; i < n_counted; i++) {
            const double a = lq[bq[i]], b = lnq[bq[i]];
            int top = (i + 1 < k) ? i + 1 : k;
            for (int j = top; j >= 1; j--) dp[j] = logaddexp(dp[j] + b, dp[j - 1] + a);
            dp[0] += b;
        }
        res[h] = dp[k];
    }
    return pb_llr_from(res[0], res[1], res[2], n_counted);
}

}  // namespace dnav2
