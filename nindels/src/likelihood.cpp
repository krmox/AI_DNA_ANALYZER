#include "nindels/likelihood.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace nindels {

namespace {

constexpr double kLn10 = 2.302585092994046;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double phred_to_p(int q) { return std::pow(10.0, -q / 10.0); }

double logsumexp(double a, double b) {
    if (a == kNegInf) return b;
    if (b == kNegInf) return a;
    double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
}

struct Transitions {
    double mm, gm, mi, ii;  // match->match, gap->match, match->gap, gap->gap
};

Transitions transitions(const PairHmmParams& p) {
    double go = phred_to_p(p.gap_open_q), ge = phred_to_p(p.gap_ext_q);
    return {1.0 - 2.0 * go, 1.0 - ge, go, ge};
}

// Base-quality error probabilities, precomputed (pow() in the inner loop was
// the dominant cost in the first profile).
struct QualTable {
    double match[64], mismatch[64];
    QualTable() {
        for (int q = 0; q < 64; ++q) {
            double e = phred_to_p(std::clamp(q, 2, 60));
            match[q] = 1.0 - e;
            mismatch[q] = e / 3.0;
        }
    }
};
const QualTable kQual;

double prior(char r, char h, uint8_t q) {
    const int qq = std::min<int>(q, 63);
    return (r == h || r == 'N' || h == 'N') ? kQual.match[qq] : kQual.mismatch[qq];
}

// Linear-space recursion with a large initial scale (as GATK does); returns
// ln P or NaN on underflow. With band >= 0 only hap columns within
// [i + diag - band, i + diag + band] are computed for read base i (1-based).
double pairhmm_scaled(std::string_view r, const uint8_t* q, std::string_view h, const Transitions& t,
                      int64_t diag, int64_t band) {
    const int64_t n = static_cast<int64_t>(r.size()), m = static_cast<int64_t>(h.size());
    const double kScale = std::ldexp(1.0, 1020);
    std::vector<double> M0(static_cast<size_t>(m + 2), 0), I0(static_cast<size_t>(m + 2), 0),
        D0(static_cast<size_t>(m + 2), kScale / static_cast<double>(m));
    std::vector<double> M1(static_cast<size_t>(m + 2), 0), I1(static_cast<size_t>(m + 2), 0),
        D1(static_cast<size_t>(m + 2), 0);
    D0[0] = D0[static_cast<size_t>(m + 1)] = 0;
    int64_t lo = 1, hi = m;
    for (int64_t i = 1; i <= n; ++i) {
        if (band >= 0) {
            lo = std::max<int64_t>(1, i + diag - band);
            hi = std::min<int64_t>(m, i + diag + band);
            if (lo > hi) return std::numeric_limits<double>::quiet_NaN();
        }
        const double pm = kQual.match[std::min<int>(q[i - 1], 63)], px = kQual.mismatch[std::min<int>(q[i - 1], 63)];
        const char rb = r[static_cast<size_t>(i - 1)];
        M1[static_cast<size_t>(lo - 1)] = I1[static_cast<size_t>(lo - 1)] = D1[static_cast<size_t>(lo - 1)] = 0;
        for (int64_t j = lo; j <= hi; ++j) {
            const size_t J = static_cast<size_t>(j);
            const char hb = h[J - 1];
            const double e = (rb == hb || rb == 'N' || hb == 'N') ? pm : px;
            M1[J] = e * (t.mm * M0[J - 1] + t.gm * (I0[J - 1] + D0[J - 1]));
            I1[J] = t.mi * M0[J] + t.ii * I0[J];
            D1[J] = t.mi * M1[J - 1] + t.ii * D1[J - 1];
        }
        if (hi + 1 <= m + 1) M1[static_cast<size_t>(hi + 1)] = I1[static_cast<size_t>(hi + 1)] = D1[static_cast<size_t>(hi + 1)] = 0;
        std::swap(M0, M1);
        std::swap(I0, I1);
        std::swap(D0, D1);
    }
    double s = 0;
    for (int64_t j = lo; j <= hi; ++j) s += M0[static_cast<size_t>(j)] + I0[static_cast<size_t>(j)];
    if (!(s > 0) || !std::isfinite(s)) return std::numeric_limits<double>::quiet_NaN();
    return std::log(s) - std::log(kScale);
}

// Same recursion in log space; slow path used only when the scaled one underflows.
double pairhmm_log(std::string_view r, const uint8_t* q, std::string_view h, const Transitions& t) {
    const size_t n = r.size(), m = h.size();
    const double lmm = std::log(t.mm), lgm = std::log(t.gm), lmi = std::log(t.mi), lii = std::log(t.ii);
    std::vector<double> M0(m + 1, kNegInf), I0(m + 1, kNegInf), D0(m + 1, -std::log(static_cast<double>(m)));
    std::vector<double> M1(m + 1), I1(m + 1), D1(m + 1);
    for (size_t i = 1; i <= n; ++i) {
        M1[0] = I1[0] = D1[0] = kNegInf;
        for (size_t j = 1; j <= m; ++j) {
            double from = logsumexp(lmm + M0[j - 1], lgm + logsumexp(I0[j - 1], D0[j - 1]));
            M1[j] = std::log(prior(r[i - 1], h[j - 1], q[i - 1])) + from;
            I1[j] = logsumexp(lmi + M0[j], lii + I0[j]);
            D1[j] = logsumexp(lmi + M1[j - 1], lii + D1[j - 1]);
        }
        std::swap(M0, M1);
        std::swap(I0, I1);
        std::swap(D0, D1);
    }
    double s = kNegInf;
    for (size_t j = 1; j <= m; ++j) s = logsumexp(s, logsumexp(M0[j], I0[j]));
    return s;
}

}  // namespace

double pairhmm_loglik(std::string_view read, const uint8_t* quals, std::string_view hap, const PairHmmParams& p,
                      int64_t diag, int64_t band) {
    if (read.empty()) return 0.0;
    if (hap.empty()) return kNegInf;
    Transitions t = transitions(p);
    double v = pairhmm_scaled(read, quals, hap, t, diag, band);
    if (std::isnan(v) && band >= 0) v = pairhmm_scaled(read, quals, hap, t, 0, -1);
    if (std::isnan(v)) v = pairhmm_log(read, quals, hap, t);
    return v;
}

int64_t Haplotypes::alt_index(int64_t g) const {
    const int L = indel.length();
    if (g < indel.pos) return g - start;
    if (indel.type == IndelType::Insertion) return g - start + L;
    if (g < indel.pos + L) return indel.pos - start;  // inside the deleted block
    return g - start - L;
}

Haplotypes build_haplotypes(const RefSeq& ref, const Indel& x, int64_t start, int64_t end) {
    start = std::max(start, ref.offset);
    end = std::min(end, ref.end());
    Haplotypes h;
    h.start = start;
    h.indel = x;
    h.ref = ref.sub(start, end - start);
    const int64_t cut = x.pos - start;
    h.alt = h.ref.substr(0, static_cast<size_t>(std::max<int64_t>(0, cut)));
    if (x.type == IndelType::Insertion) {
        h.alt += x.seq;
        h.alt += h.ref.substr(static_cast<size_t>(cut));
    } else {
        const int64_t after = cut + x.length();
        if (after < static_cast<int64_t>(h.ref.size())) h.alt += h.ref.substr(static_cast<size_t>(after));
    }
    return h;
}

std::vector<ReadLik> read_likelihoods(const RefSeq& ref, const Indel& x, const std::vector<HmmRead>& reads,
                                      const LikelihoodParams& p, std::vector<size_t>* used) {
    std::vector<size_t> order(reads.size());
    std::iota(order.begin(), order.end(), 0);
    if (static_cast<int>(order.size()) > p.max_reads) {
        std::nth_element(order.begin(), order.begin() + p.max_reads, order.end(),
                         [&](size_t a, size_t b) { return reads[a].name_hash < reads[b].name_hash; });
        order.resize(static_cast<size_t>(p.max_reads));
        std::sort(order.begin(), order.end());
    }
    std::vector<ReadLik> out;
    if (order.empty()) return out;

    int64_t lo = x.pos, hi = x.pos + x.ref_span();
    for (size_t k : order) {
        lo = std::min(lo, reads[k].ref_start);
        hi = std::max(hi, reads[k].ref_end);
    }
    const int pad = p.hap_pad + x.length();
    Haplotypes H = build_haplotypes(ref, x, lo - pad, hi + pad);
    const int64_t nref = static_cast<int64_t>(H.ref.size()), nalt = static_cast<int64_t>(H.alt.size());

    out.reserve(order.size());
    for (size_t k : order) {
        const HmmRead& r = reads[k];
        const int64_t g0 = r.ref_start - pad, g1 = r.ref_end + pad;
        int64_t a0 = std::clamp<int64_t>(g0 - H.start, 0, nref), a1 = std::clamp<int64_t>(g1 - H.start, 0, nref);
        int64_t b0 = std::clamp<int64_t>(H.alt_index(g0), 0, nalt), b1 = std::clamp<int64_t>(H.alt_index(g1), 0, nalt);
        std::string_view href(H.ref.data() + a0, static_cast<size_t>(a1 - a0));
        std::string_view halt(H.alt.data() + b0, static_cast<size_t>(b1 - b0));
        ReadLik L;
        // Banded around the read's mapped diagonal: read base 1 sits at slice
        // column (ref_start - slice start) + 1 on H_ref, and likewise on H_alt.
        const int64_t band = p.band < 0 ? -1 : p.band + x.length();
        const int64_t diag_ref = (r.ref_start - H.start) - a0;
        const int64_t diag_alt = std::clamp<int64_t>(H.alt_index(r.ref_start), 0, nalt) - b0;
        L.lp_ref = pairhmm_loglik(r.bases, r.quals.data(), href, p.hmm, diag_ref, band);
        L.lp_alt = pairhmm_loglik(r.bases, r.quals.data(), halt, p.hmm, diag_alt, band);
        L.mapq = r.mapq;
        L.cigar_allele = r.cigar_allele;
        out.push_back(L);
    }
    if (used) *used = std::move(order);
    return out;
}

double delta_ell_original(const std::vector<ReadLik>& liks) {
    double a = 0, r = 0;
    for (const ReadLik& l : liks) { a += l.lp_alt; r += l.lp_ref; }
    return a - r;
}

static void finish(GenotypeLik& g, const GenotypePriors& pr) {
    double lp[3], z = kNegInf;
    for (int k = 0; k < 3; ++k) { lp[k] = g.gl[k] + pr.ln_prior[k]; z = logsumexp(z, lp[k]); }
    for (int k = 0; k < 3; ++k) g.posterior[k] = std::exp(lp[k] - z);
    g.gt = static_cast<int>(std::max_element(g.posterior.begin(), g.posterior.end()) - g.posterior.begin());
    double err = std::max(1e-10, 1.0 - g.posterior[g.gt]);
    g.gq = std::min(99.0, -10.0 * std::log10(err));
    g.qual = logsumexp(g.gl[1], g.gl[2]) - g.gl[0];
}

GenotypeLik diploid_gl(const std::vector<ReadLik>& liks, bool mapq_cap, const GenotypePriors& pr) {
    GenotypeLik g;
    for (const ReadLik& l : liks) {
        double r = l.lp_ref, a = l.lp_alt;
        if (mapq_cap) {
            // ln[(1-eps) P + eps * max(P_ref, P_alt)], eps = 10^(-MAPQ/10)
            double eps = std::min(1.0, phred_to_p(l.mapq));
            double floor = std::max(r, a) + std::log(eps);
            double keep = eps < 1.0 ? std::log1p(-eps) : kNegInf;
            r = logsumexp(keep + r, floor);
            a = logsumexp(keep + a, floor);
        }
        g.gl[0] += r;
        g.gl[1] += logsumexp(r, a) - std::log(2.0);
        g.gl[2] += a;
    }
    finish(g, pr);
    return g;
}

GenotypeLik count_gl(int n_alt, int n_ref, double err, const GenotypePriors& pr) {
    GenotypeLik g;
    const double p[3] = {err, 0.5, 1.0 - err};
    for (int k = 0; k < 3; ++k) g.gl[k] = n_alt * std::log(p[k]) + n_ref * std::log1p(-p[k]);
    finish(g, pr);
    return g;
}

Opt hmm_consistency(const std::vector<ReadLik>& liks) {
    int n = 0, agree = 0;
    for (const ReadLik& l : liks) {
        if (l.cigar_allele == ReadAllele::Other || l.lp_alt == l.lp_ref) continue;
        ++n;
        agree += (l.cigar_allele == ReadAllele::Alt) == l.best_alt();
    }
    if (n == 0) return std::nullopt;
    return static_cast<double>(agree) / n;
}

}  // namespace nindels
