// Frozen Method C genotyping (experimental/chr20_validation optimized_method_c, eps = 0.01),
// one record at a time.  Copied semantics; requires -ffp-contract=off.
#pragma once
#include <cmath>
#include <cstdint>
#include "binomial.hpp"
#include "numerics.hpp"

namespace dnav2 {

constexpr double GQ_CAP = 99.0;

struct MethodCResult {
    int gt;        // 0 = "0/1"(forced from 0/0), 1 = "0/1", 2 = "1/1" reported class index into GT_STR
    int gq;        // int(round(float(round(gq, 1))))
    bool forced;   // best class was 0/0 and was forced to 0/1
};
inline const char* gt_string(int gt) { return gt == 2 ? "1/1" : "0/1"; }

inline MethodCResult method_c(double k, double n) {
    static const double EPS = FROZEN_EPSILON;
    static const double p[3] = {clip(0.0 * (1 - EPS) + 1.0 * EPS, 1e-12, 1 - 1e-12),
                                clip(0.5 * (1 - EPS) + 0.5 * EPS, 1e-12, 1 - 1e-12),
                                clip(1.0 * (1 - EPS) + 0.0 * EPS, 1e-12, 1 - 1e-12)};
    static const double lp[3] = {std::log(p[0]), std::log(p[1]), std::log(p[2])};
    static const double l1p[3] = {std::log(1 - p[0]), std::log(1 - p[1]), std::log(1 - p[2])};
    static const double LN10 = std::log(10.0);

    double logp[3];
    if (n == 0) {
        logp[0] = logp[1] = logp[2] = 0.0;
    } else {
        double coeff = (gammaln(n + 1) - gammaln(k + 1)) - gammaln(n - k + 1);
        for (int c = 0; c < 3; c++) logp[c] = (coeff + k * lp[c]) + (n - k) * l1p[c];
    }
    // argsort(-logp): stable for 3 elements (numpy insertion sort); first two entries
    int best = 0;
    for (int c = 1; c < 3; c++) if (logp[c] > logp[best]) best = c;
    int second = -1;
    for (int c = 0; c < 3; c++) {
        if (c == best) continue;
        if (second < 0 || logp[c] > logp[second]) second = c;
    }
    double gq_raw = std::fmin(10.0 * (logp[best] - logp[second]) / LN10, GQ_CAP);
    gq_raw = std::fmax(gq_raw, 0.0);
    bool theta0 = best == 0;
    double gq_capped = theta0 ? std::fmin(gq_raw, 5.0) : gq_raw;
    double gq_out = std::rint(gq_capped * 10.0) / 10.0;        // np.round(x, 1)
    int gq_int = static_cast<int>(std::nearbyint(gq_out));      // python round(): half-even
    return {best == 2 ? 2 : 1, gq_int, theta0};
}

}  // namespace dnav2
