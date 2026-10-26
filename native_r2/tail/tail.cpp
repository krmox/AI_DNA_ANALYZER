#include "tail.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace dnav2 {

std::string ood_reason(const Calibration& c, const OodEnvelope& e) {
    if (!c.valid) return "calibration_failed";
    if (!e.set) return "";
    char b[160];
    auto out = [&](const char* what, double v, double lo, double hi) {
        std::snprintf(b, sizeof b, "%s=%.4g_outside_[%.4g,%.4g]", what, v, lo, hi);
        return std::string(b);
    };
    if (c.eps_hat < e.eps_lo || c.eps_hat > e.eps_hi) return out("eps", c.eps_hat, e.eps_lo, e.eps_hi);
    if (c.mq_p90 < e.mq_lo || c.mq_p90 > e.mq_hi) return out("mapq_p90", c.mq_p90, e.mq_lo, e.mq_hi);
    if (c.bq_p50 < e.bq_lo || c.bq_p50 > e.bq_hi) return out("bq_p50", c.bq_p50, e.bq_lo, e.bq_hi);
    if (c.depth_med < e.depth_lo || c.depth_med > e.depth_hi) return out("depth_med", c.depth_med, e.depth_lo, e.depth_hi);
    return "";
}

double hist_quantile(const std::vector<uint64_t>& h, double q) {
    uint64_t total = 0;
    for (uint64_t v : h) total += v;
    if (total == 0) return 0.0;
    const double target = q * static_cast<double>(total);
    uint64_t acc = 0;
    for (size_t i = 0; i < h.size(); i++) {
        acc += h[i];
        if (static_cast<double>(acc) >= target) return static_cast<double>(i);
    }
    return static_cast<double>(h.size() - 1);
}

void pattern_features(const std::vector<CandPoint>& c, int64_t span_lo, int64_t span_hi, const Calibration& cal,
                      std::vector<PatternOut>& out) {
    const size_t n = c.size();
    out.assign(n, PatternOut{0.0, 0.0, 0.0});
    if (n == 0) return;
    // expected counts floor at 0.5 so that sparse expectations do not explode the ratio
    const double rate = cal.cand_per_bp;
    constexpr int64_t R1 = 100, R2 = 5000;
    constexpr double SHRINK_PRIOR = 2.0, SHRINK_VAF = 0.2;
    size_t lo1 = 0, hi1 = 0, lo2 = 0, hi2 = 0;   // [lo,hi) index ranges within +-R of candidate j
    std::vector<double> lowv_prefix(n + 1, 0.0);
    for (size_t j = 0; j < n; j++) lowv_prefix[j + 1] = lowv_prefix[j] + (c[j].vaf < LOWVAF_CUT ? 1.0 : 0.0);
    for (size_t j = 0; j < n; j++) {
        const int64_t p = c[j].pos;
        while (lo1 < n && c[lo1].pos < p - R1) lo1++;
        while (hi1 < n && c[hi1].pos <= p + R1) hi1++;
        while (lo2 < n && c[lo2].pos < p - R2) lo2++;
        while (hi2 < n && c[hi2].pos <= p + R2) hi2++;
        const double n1 = static_cast<double>(hi1 - lo1) - 1.0;   // exclude self
        const double n2 = static_cast<double>(hi2 - lo2) - 1.0;
        const int64_t w1 = std::min(p + R1 + 1, span_hi) - std::max(p - R1, span_lo) - 1;
        const int64_t w2 = std::min(p + R2 + 1, span_hi) - std::max(p - R2, span_lo) - 1;
        out[j].dens100 = n1 / std::max(0.5, rate * static_cast<double>(std::max<int64_t>(w1, 1)));
        out[j].dens5k = n2 / std::max(0.5, rate * static_cast<double>(std::max<int64_t>(w2, 1)));
        double low = lowv_prefix[hi2] - lowv_prefix[lo2] - (c[j].vaf < LOWVAF_CUT ? 1.0 : 0.0);
        out[j].lowvaf5k = (low + SHRINK_VAF * SHRINK_PRIOR) / (n2 + SHRINK_PRIOR);
    }
}

}  // namespace dnav2
