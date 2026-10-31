#include "nindels/reassess.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nindels {

static void check_simplex(std::initializer_list<double> w, const char* what) {
    double s = 0;
    for (double v : w) {
        if (!(v >= 0)) throw std::invalid_argument(std::string(what) + ": weights must be >= 0");
        s += v;
    }
    if (std::fabs(s - 1.0) > 1e-9) throw std::invalid_argument(std::string(what) + ": weights must sum to 1");
}

void ComplexityParams::validate() const {
    check_simplex({v_hp, v_repeat, v_lowcomplexity, v_alignment}, "complexity");
}

void ReassessParams::validate() const {
    check_simplex({a_U, a_C, a_D}, "reassessment");
    if (!(A >= 0 && A <= 1)) throw std::invalid_argument("A must be in [0,1]");
    if (!(tau_loose <= tau_strict)) throw std::invalid_argument("tau_loose must be <= tau_strict");
}

Complexity region_complexity(const Pattern& p, Opt instability, const ComplexityParams& cp) {
    Complexity c;
    c.H = std::min(1.0, p.hp_length / cp.hp_norm);
    c.R = std::min(1.0, p.repeat_span / cp.span_norm);
    c.L = std::clamp(1.0 - p.entropy / 2.0, 0.0, 1.0);
    c.G = instability;
    double num = cp.v_hp * c.H + cp.v_repeat * c.R + cp.v_lowcomplexity * c.L;
    double den = cp.v_hp + cp.v_repeat + cp.v_lowcomplexity;
    if (c.G) { num += cp.v_alignment * *c.G; den += cp.v_alignment; }
    c.C = den > 0 ? num / den : 0.0;  // explicit renormalisation when G is missing
    return c;
}

NucleotideResult nucleotide_scores(const std::vector<const HmmRead*>& reads, const std::vector<ReadLik>& liks,
                                   Opt f, const NucleotideParams& p) {
    const int W = 2 * p.half_window + 1;
    NucleotideResult out;
    out.q.assign(static_cast<size_t>(W), std::nullopt);
    for (int k = 0; k < W; ++k) {
        double bq = 0, mq = 0, agree = 0, alt_best = 0;
        int n = 0, n_decided = 0;
        for (size_t i = 0; i < reads.size(); ++i) {
            const HmmRead& r = *reads[i];
            if (static_cast<int>(r.window_bq.size()) != W || r.window_bq[k] < 0) continue;
            ++n;
            bq += r.window_bq[k];
            mq += r.mapq;
            alt_best += liks[i].best_alt();
            if (r.cigar_allele != ReadAllele::Other) {
                ++n_decided;
                agree += (r.cigar_allele == ReadAllele::Alt) == liks[i].best_alt();
            }
        }
        if (n == 0) continue;
        double num = p.w_bq * std::min(1.0, bq / n / p.bq_norm) + p.w_mq * std::min(1.0, mq / n / p.mq_norm);
        double den = p.w_bq + p.w_mq;
        if (n_decided > 0) { num += p.w_align * agree / n_decided; den += p.w_align; }
        if (f) { num += p.w_support * (1.0 - std::fabs(alt_best / n - *f)); den += p.w_support; }
        out.q[static_cast<size_t>(k)] = den > 0 ? num / den : 0.0;
    }
    for (int k = 0; k < W; ++k) {
        const Opt& v = out.q[static_cast<size_t>(k)];
        if (v && (!out.q_min || *v < *out.q_min)) { out.q_min = v; out.argmin_offset = k - p.half_window; }
    }
    return out;
}

}  // namespace nindels
