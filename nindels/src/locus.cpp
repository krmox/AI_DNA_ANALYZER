#include "nindels/locus.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <unordered_set>

#include "nindels/str_model.hpp"
#include "nindels/table_io.hpp"

namespace nindels {

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double lse(double a, double b) {
    if (a == kNegInf) return b;
    if (b == kNegInf) return a;
    double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
}

int n_alt_reads(const CandidateInput& c) {
    int n = 0;
    for (const auto& o : c.obs) n += o.allele == ReadAllele::Alt;
    return n;
}

}  // namespace

std::vector<std::vector<size_t>> cluster_candidates(const RefSeq& ref, const std::vector<const CandidateInput*>& in) {
    std::vector<std::vector<size_t>> groups;
    int64_t group_end = std::numeric_limits<int64_t>::min();
    for (size_t i = 0; i < in.size(); ++i) {
        const Indel& x = in[i]->indel;
        const Pattern pat = extract_pattern(ref, x);
        const int64_t a = x.pos - 1;
        const int64_t b = x.pos + std::max({x.ref_span(), pat.repeat_span, pat.hp_length}) + 1;
        if (groups.empty() || a > group_end) {
            groups.push_back({i});
            group_end = b;
        } else {
            groups.back().push_back(i);
            group_end = std::max(group_end, b);
        }
    }
    return groups;
}

LocusCall genotype_locus(const RefSeq& ref, const std::vector<const CandidateInput*>& group_in,
                         const LikelihoodParams& lp, const LocusParams& p) {
    // Alleles: the K best-supported candidates.
    std::vector<const CandidateInput*> group = group_in;
    std::stable_sort(group.begin(), group.end(),
                     [](const CandidateInput* x, const CandidateInput* y) { return n_alt_reads(*x) > n_alt_reads(*y); });
    if (static_cast<int>(group.size()) > p.max_alleles) group.resize(static_cast<size_t>(p.max_alleles));
    LocusCall call;
    for (const auto* c : group) call.alleles.push_back(c->indel);
    const int K = static_cast<int>(call.alleles.size());

    // Reads: union over candidates (one entry per read), capped by name hash.
    std::vector<const HmmRead*> reads;
    std::unordered_set<uint64_t> seen;
    for (const auto* c : group_in)
        for (const auto& r : c->reads) {
            uint64_t key = mix64(r.name_hash ^ static_cast<uint64_t>(r.ref_start) * 0x9E3779B97F4A7C15ULL);
            if (seen.insert(key).second) reads.push_back(&r);
        }
    std::sort(reads.begin(), reads.end(), [](const HmmRead* a, const HmmRead* b) {
        return a->name_hash != b->name_hash ? a->name_hash < b->name_hash : a->ref_start < b->ref_start;
    });
    if (static_cast<int>(reads.size()) > lp.max_reads) reads.resize(static_cast<size_t>(lp.max_reads));
    call.depth = static_cast<int>(reads.size());
    call.ad.assign(static_cast<size_t>(K + 1), 0);
    if (reads.empty()) return call;

    // Haplotypes over a common reference window.
    int max_len = 0;
    int64_t lo = std::numeric_limits<int64_t>::max(), hi = std::numeric_limits<int64_t>::min();
    for (const auto& x : call.alleles) {
        max_len = std::max(max_len, x.length());
        lo = std::min(lo, x.pos);
        hi = std::max(hi, x.pos + x.ref_span());
    }
    for (const auto* r : reads) { lo = std::min(lo, r->ref_start); hi = std::max(hi, r->ref_end); }
    const int pad = lp.hap_pad + max_len;
    std::vector<Haplotypes> H;  // H[0].ref is the reference haplotype; H[a-1].alt is allele a
    for (const auto& x : call.alleles) H.push_back(build_haplotypes(ref, x, lo - pad, hi + pad));
    std::vector<const std::string*> hap(static_cast<size_t>(K + 1));
    hap[0] = &H[0].ref;
    for (int a = 1; a <= K; ++a) hap[static_cast<size_t>(a)] = &H[static_cast<size_t>(a - 1)].alt;
    const int64_t ws = H[0].start;
    std::vector<std::vector<double>> gop(static_cast<size_t>(K + 1));
    if (lp.str) {
        const double base = std::pow(10.0, -lp.hmm.gap_open_q / 10.0);
        for (int a = 0; a <= K; ++a) gop[static_cast<size_t>(a)] = lp.str->gap_open_profile(*hap[static_cast<size_t>(a)], base);
    }
    // reference coordinate -> index on haplotype a
    auto index_on = [&](int a, int64_t g) -> int64_t {
        return a == 0 ? g - ws : H[static_cast<size_t>(a - 1)].alt_index(g);
    };
    const int64_t band = lp.band < 0 ? -1 : lp.band + max_len;

    // l[r][a] = ln P(read r | H_a), then the MAPQ mixture (mapping uncertainty).
    std::vector<std::vector<double>> l(reads.size(), std::vector<double>(static_cast<size_t>(K + 1)));
    for (size_t r = 0; r < reads.size(); ++r) {
        const HmmRead& rd = *reads[r];
        double best = kNegInf;
        for (int a = 0; a <= K; ++a) {
            const std::string& h = *hap[static_cast<size_t>(a)];
            const int64_t n = static_cast<int64_t>(h.size());
            const int64_t b0 = std::clamp<int64_t>(index_on(a, rd.ref_start - pad), 0, n);
            const int64_t b1 = std::clamp<int64_t>(index_on(a, rd.ref_end + pad), 0, n);
            const int64_t diag = std::clamp<int64_t>(index_on(a, rd.ref_start), 0, n) - b0;
            std::string_view slice(h.data() + b0, static_cast<size_t>(std::max<int64_t>(0, b1 - b0)));
            double v = pairhmm_loglik(rd.bases, rd.quals.data(), slice, lp.hmm, diag, band,
                                      lp.str ? gop[static_cast<size_t>(a)].data() + b0 : nullptr);
            l[r][static_cast<size_t>(a)] = v;
            best = std::max(best, v);
        }
        const double eps = std::min(1.0, std::pow(10.0, -rd.mapq / 10.0));
        const double keep = eps < 1.0 ? std::log1p(-eps) : kNegInf, flo = best + std::log(eps);
        int arg = -1, ties = 0;
        double top = kNegInf;
        for (int a = 0; a <= K; ++a) {
            double& v = l[r][static_cast<size_t>(a)];
            if (v > top + 1e-9) { top = v; arg = a; ties = 0; }
            else if (std::fabs(v - top) <= 1e-9) ++ties;
            v = lse(keep + v, flo);
        }
        if (arg >= 0 && ties == 0) ++call.ad[static_cast<size_t>(arg)];
    }

    // Genotype likelihoods and posteriors over unordered pairs.
    struct G { int a, b; double lp; };
    std::vector<G> gts;
    const double ln2 = std::log(2.0);
    for (int a = 0; a <= K; ++a)
        for (int b = a; b <= K; ++b) {
            double gl = 0;
            for (size_t r = 0; r < reads.size(); ++r)
                gl += a == b ? l[r][static_cast<size_t>(a)]
                             : lse(l[r][static_cast<size_t>(a)], l[r][static_cast<size_t>(b)]) - ln2;
            double prior = (a == 0 && b == 0) ? 0.0
                         : (a == 0) ? std::log(p.theta)
                         : (a == b) ? std::log(p.theta / 2)
                                    : std::log(p.theta * p.hetalt_factor);
            gts.push_back({a, b, gl + prior});
        }
    double z = kNegInf;
    for (const G& g : gts) z = lse(z, g.lp);
    size_t best = 0, best_nonref = 1;
    for (size_t i = 0; i < gts.size(); ++i) {
        if (gts[i].lp > gts[best].lp) best = i;
        if (i > 0 && gts[i].lp > gts[best_nonref].lp) best_nonref = i;
    }
    const double p00 = std::exp(gts[0].lp - z);
    call.qual = std::min(9999.0, -10.0 * std::log10(std::max(p00, 1e-300)));
    const G& chosen = best == 0 ? gts[best_nonref] : gts[best];
    call.g1 = chosen.a;
    call.g2 = chosen.b;
    call.gq = std::min(99.0, -10.0 * std::log10(std::max(1e-10, 1.0 - std::exp(chosen.lp - z))));
    call.pass = best != 0 && call.qual >= p.min_qual;
    return call;
}

void write_locus_vcf(std::ostream& o, const std::string& chrom, const RefSeq& ref, const LocusCall& c) {
    if (c.alleles.empty() || (c.g1 == 0 && c.g2 == 0)) return;
    std::vector<int> used;
    for (int a : {c.g1, c.g2})
        if (a > 0 && std::find(used.begin(), used.end(), a) == used.end()) used.push_back(a);
    int64_t anchor = std::numeric_limits<int64_t>::max(), end = 0;
    for (int a : used) {
        const Indel& x = c.alleles[static_cast<size_t>(a - 1)];
        anchor = std::min(anchor, x.pos - 1);
        end = std::max(end, x.pos + x.ref_span());
    }
    const std::string REF = ref.sub(anchor, end - anchor);
    std::vector<std::string> alts;
    for (int a : used) {
        const Indel& x = c.alleles[static_cast<size_t>(a - 1)];
        std::string s = ref.sub(anchor, x.pos - anchor);
        if (x.type == IndelType::Insertion) s += x.seq;
        s += ref.sub(x.pos + x.ref_span(), end - (x.pos + x.ref_span()));
        alts.push_back(s);
    }
    auto idx = [&](int a) {
        if (a == 0) return 0;
        return static_cast<int>(std::find(used.begin(), used.end(), a) - used.begin()) + 1;
    };
    o << chrom << '\t' << anchor + 1 << "\t.\t" << REF << '\t';
    for (size_t i = 0; i < alts.size(); ++i) o << (i ? "," : "") << alts[i];
    o << '\t' << format_opt(std::round(c.qual * 100) / 100) << '\t' << (c.pass ? "PASS" : "LowProb") << '\t'
      << "K=" << c.alleles.size() << ";DP=" << c.depth << "\tGT:GQ:AD\t" << idx(c.g1) << '/' << idx(c.g2) << ':'
      << static_cast<int>(std::round(c.gq)) << ':' << c.ad[0];
    for (int a : used) o << ',' << c.ad[static_cast<size_t>(a)];
    o << '\n';
}

}  // namespace nindels
