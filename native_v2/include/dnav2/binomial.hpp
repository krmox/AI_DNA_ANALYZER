// Frozen V1 Binomial caller (binomial_baseline.BinomialVariantCaller.score_counts with the fixed
// error_rate=0.01, include_homozygous=True) for ONE locus, plus the frozen router/threshold
// constants.  Semantics are copied, not redesigned.
#pragma once
#include <cmath>
#include <cstdint>
#include <vector>
#include "numerics.hpp"

namespace dnav2 {

// ---- FROZEN CONSTANTS (never tuned) ---------------------------------------------------------
constexpr double FROZEN_BINOMIAL_THRESHOLD = 7.0;
constexpr double FROZEN_PB_THRESHOLD = 10.5;
constexpr double FROZEN_ROUTER_CUTOFF = 5.411872376933351;
constexpr double FROZEN_EPSILON = 0.01;      // Binomial fixed error rate AND Method C eps
constexpr int MAX_READS = 48;                // frozen read-level cap

// Result of the Binomial scoring of one locus
struct BinomialScore {
    double n;        // usable A/C/G/T observations (float, as in V1)
    double k;        // observations supporting the candidate ALT
    int alt_index;   // 0..3 = A,C,G,T
    double llr;
};

// counts4 = count_A..count_T; ref_index: 0=N, 1..4 = A,C,G,T (pileup_counts convention)
inline BinomialScore binomial_score(const double counts4[4], int ref_index) {
    static const double eps = FROZEN_EPSILON;
    static const double p0 = clip(eps / 3.0, 1e-12, 1 - 1e-12);
    static const double p_het = clip(0.5 * (1.0 - eps) + 0.5 * (eps / 3.0), 1e-12, 1 - 1e-12);
    static const double p_hom = clip(1.0 - eps + eps / 3.0, 1e-12, 1 - 1e-12);

    double n = ((counts4[0] + counts4[1]) + counts4[2]) + counts4[3];
    // candidate ALT: masked argmax (first maximum wins, as numpy.argmax)
    double masked[4] = {counts4[0], counts4[1], counts4[2], counts4[3]};
    if (ref_index >= 1 && ref_index <= 4) masked[ref_index - 1] = -1.0;
    int alt = 0;
    double mx = masked[0];
    for (int i = 1; i < 4; i++)
        if (masked[i] > mx) { mx = masked[i]; alt = i; }
    double k = counts4[alt];
    if (mx < 0) k = 0.0;

    // scipy.stats.binom.logpmf: combiln + xlogy(k,p) + xlog1py(n-k,-p)
    auto logpmf = [&](double p) {
        double combiln = gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1));
        return (combiln + xlogy(k, p)) + xlog1py(n - k, -p);
    };
    double h0 = logpmf(p0), het = logpmf(p_het), hom = logpmf(p_hom);
    double alt_ll = het > hom ? het : hom;   // np.maximum (no NaN reachable)
    double llr = alt_ll - h0;
    if (n <= 0) llr = 0.0;
    if (!std::isfinite(llr)) llr = 0.0;      // nan_to_num(nan=0, posinf=0, neginf=0)
    return {n, k, alt, llr};
}

inline bool is_routed(double binomial_llr) {
    return std::fabs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF;
}

}  // namespace dnav2
