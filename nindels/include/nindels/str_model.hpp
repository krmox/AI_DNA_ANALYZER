// Experiment B: repeat-aware indel error model.
//
// eps(p, u) = probability that a read shows a whole-unit slippage in a short
// tandem repeat of period p and u units. It is estimated from the BAM being
// called, without truth labels: binomial-mixture EM over sampled STR loci
// (hom-ref: Bin(n, eps), het: Bin(n, 1/2), hom-alt: Bin(n, 1 - eps)).
// The PairHMM then uses a position-specific gap-open probability inside each
// STR, calibrated so that the HMM's one-unit-slip probability matches eps.
#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "nindels/likelihood.hpp"

namespace nindels {

// A maximal tandem run [start, end) of period p on a sequence.
struct StrRun {
    int64_t start = 0, end = 0;
    int period = 0;
    int units() const { return static_cast<int>((end - start) / period); }
};

// Minimum units for a run to count as an STR: homopolymers >= 4 bp, other periods >= 3 units.
// ARBITRARY (shorter runs keep the default gap-open).
int str_min_units(int period);

// Maximal runs of period 1..6, keeping for each position the longest run that
// covers it (ties: smaller period). Only runs with >= str_min_units units.
std::vector<StrRun> find_str_runs(std::string_view seq);

struct StrLocusCount {
    int period = 0, units = 0;
    int n = 0;  // spanning fragments
    int s = 0;  // fragments with a whole-unit indel inside the run
};

class StrErrorModel {
public:
    static constexpr int kMaxPeriod = 6;
    static constexpr int kMaxUnits = 40;  // units >= 40 share one stratum

    // EM per (p, u) stratum; strata with < min_loci loci inherit from (p, u-1).
    void fit(const std::vector<StrLocusCount>& counts, int min_loci = 30);
    // Gap-open per stratum such that P(read with one-unit deletion | H) / P(clean read | H)
    // equals eps / (2 (1 - eps)) on a synthetic STR (insertions get the same open probability).
    void calibrate(const PairHmmParams& hmm);

    // Gap-open probability per haplotype position (index j = gap opening at/after hap base j).
    std::vector<double> gap_open_profile(std::string_view hap, double base_gop) const;

    double eps(int p, int u) const { return eps_[p][clamp_u(u)]; }
    double gop(int p, int u) const { return gop_[p][clamp_u(u)]; }
    int loci(int p, int u) const { return loci_[p][clamp_u(u)]; }

    void save(const std::string& path) const;
    static StrErrorModel load(const std::string& path);

private:
    static int clamp_u(int u) { return u < 0 ? 0 : (u > kMaxUnits ? kMaxUnits : u); }
    double eps_[kMaxPeriod + 1][kMaxUnits + 1] = {};   // < 0: unknown
    double gop_[kMaxPeriod + 1][kMaxUnits + 1] = {};   // < 0: use the default
    int loci_[kMaxPeriod + 1][kMaxUnits + 1] = {};
    bool fitted_ = false;
};

}  // namespace nindels
