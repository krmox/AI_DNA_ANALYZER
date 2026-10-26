// Varion R2 "tail": per-sample self-calibration, out-of-distribution guard and the three-scale pattern layer.
// Pure functions where possible (no I/O, no globals) so they can be unit-tested and run under sanitizers.
// Nothing here reads truth labels.
#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace dnav2 {

// Sample statistics measured on a deterministic ~1 % sample of the contig (see calibrate.cpp).
struct Calibration {
    bool valid = false;
    int64_t n_cols = 0;          // sampled columns that entered the statistics
    int64_t n_bases = 0;         // bases that entered the error-rate estimate
    double depth_mean = 0, depth_med = 0;   // filtered depth (MAPQ>=20), per sampled column
    double mq_p50 = 0, mq_p90 = 0;          // read MAPQ quantiles (all reads, before filters)
    double bq_p10 = 0, bq_p50 = 0, bq_p90 = 0;
    double eps_hat = 0.01;       // non-reference base rate in columns without a variant (clamped)
    double cand_per_bp = 0;      // candidate columns (ALT count >= 3, VAF >= 0.05) per sampled base
    double lowvaf_base = 0;      // fraction of candidates with VAF < 0.3
    int read_len_med = 0;
};

constexpr double EPS_CLAMP_LO = 0.002, EPS_CLAMP_HI = 0.03;
constexpr double LOWVAF_CUT = 0.3;

// Envelope of calibration values seen in training; a sample outside it triggers the conservative fallback caller.
struct OodEnvelope {
    bool set = false;
    double eps_lo = 0, eps_hi = 1, mq_lo = 0, mq_hi = 1e9, bq_lo = 0, bq_hi = 1e9, depth_lo = 0, depth_hi = 1e9;
};
// Empty string when inside the envelope (or no envelope / invalid calibration -> reported as a reason).
std::string ood_reason(const Calibration& c, const OodEnvelope& e);

struct CandPoint { int64_t pos; double vaf; };   // candidates of one task, sorted by pos
struct PatternOut { double dens100, dens5k, lowvaf5k; };

// Pattern layer: candidate density at +-100 bp and +-5 kb relative to the contig-wide expectation (cal.cand_per_bp), and the
// shrunk fraction of low-VAF candidates within +-5 kb (excluding the candidate itself). Windows are clipped to [span_lo, span_hi).
void pattern_features(const std::vector<CandPoint>& cands, int64_t span_lo, int64_t span_hi, const Calibration& cal,
                      std::vector<PatternOut>& out);

// Quantile of an integer histogram (q in [0,1]); 0 for an empty histogram.
double hist_quantile(const std::vector<uint64_t>& h, double q);

}  // namespace dnav2
