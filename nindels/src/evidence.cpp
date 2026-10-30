#include "nindels/evidence.hpp"

#include <cmath>

namespace nindels {

static Opt ratio(double num, double den) {
    if (den <= 0) return std::nullopt;
    return num / den;
}

Evidence summarize_evidence(const std::vector<ReadObservation>& obs, const EvidenceParams& p) {
    Evidence e;
    double alt_fwd = 0, alt_rev = 0, ref_fwd = 0, ref_rev = 0;
    double rp_alt = 0, rp_ref = 0, clip_alt = 0, clip_all = 0, mm_alt = 0, split = 0;
    double mq_alt = 0, mq0 = 0, mqlow = 0, bq_alt = 0, unstable = 0;

    for (const ReadObservation& r : obs) {
        const bool clipped = r.soft_clip >= p.clip_threshold;
        const bool low = r.mapq < p.low_mapq;
        ++e.n_total;
        clip_all += clipped;
        split += r.split;
        mq0 += r.mapq == 0;
        mqlow += low;
        unstable += (low || clipped);
        switch (r.allele) {
            case ReadAllele::Alt:
                ++e.n_alt;
                (r.reverse ? alt_rev : alt_fwd) += 1;
                rp_alt += r.rel_pos;
                clip_alt += clipped;
                mm_alt += r.mismatch_rate;
                mq_alt += r.mapq;
                bq_alt += r.flank_bq;
                break;
            case ReadAllele::Ref:
                ++e.n_ref;
                (r.reverse ? ref_rev : ref_fwd) += 1;
                rp_ref += r.rel_pos;
                break;
            case ReadAllele::Other:
                ++e.n_other;
                break;
        }
    }

    e.f = ratio(e.n_alt, e.n_alt + e.n_ref);
    e.l_depth = std::log1p(static_cast<double>(e.n_alt + e.n_ref));

    e.alt_fwd_frac = ratio(alt_fwd, e.n_alt);
    if (e.n_alt > 0 && e.n_ref > 0)
        e.strand_bias = std::fabs(std::log(((alt_fwd + 0.5) * (ref_rev + 0.5)) /
                                           ((alt_rev + 0.5) * (ref_fwd + 0.5))));

    e.readpos_alt = ratio(rp_alt, e.n_alt);
    if (e.n_alt > 0 && e.n_ref > 0) e.readpos_diff = *e.readpos_alt - rp_ref / e.n_ref;

    e.clip_alt = ratio(clip_alt, e.n_alt);
    e.clip_all = ratio(clip_all, e.n_total);
    e.mismatch_alt = ratio(mm_alt, e.n_alt);
    e.other_frac = ratio(e.n_other, e.n_total);
    e.split_frac = ratio(split, e.n_total);
    e.mapq_alt = ratio(mq_alt, e.n_alt);
    e.mq0_frac = ratio(mq0, e.n_total);
    e.mqlow_frac = ratio(mqlow, e.n_total);
    e.bq_alt = ratio(bq_alt, e.n_alt);
    e.instability = ratio(unstable, e.n_total);
    return e;
}

}  // namespace nindels
