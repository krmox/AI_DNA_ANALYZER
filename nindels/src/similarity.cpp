#include "nindels/similarity.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nindels {

void SimilarityParams::validate() const {
    const double w[] = {w_type, w_length, w_motif, w_context, w_hp, w_repeat};
    double s = 0;
    for (double v : w) {
        if (!(v >= 0)) throw std::invalid_argument("similarity weight must be >= 0");
        s += v;
    }
    if (!(s > 0)) throw std::invalid_argument("sum of similarity weights must be > 0");
    if (!(s_L > 0 && s_h > 0 && s_r > 0)) throw std::invalid_argument("scale parameters must be > 0");
}

int levenshtein(std::string_view a, std::string_view b) {
    if (a.size() < b.size()) std::swap(a, b);
    std::vector<int> prev(b.size() + 1), cur(b.size() + 1);
    for (size_t j = 0; j <= b.size(); ++j) prev[j] = static_cast<int>(j);
    for (size_t i = 1; i <= a.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (size_t j = 1; j <= b.size(); ++j) {
            int sub = prev[j - 1] + (a[i - 1] != b[j - 1]);
            cur[j] = std::min({sub, prev[j] + 1, cur[j - 1] + 1});
        }
        std::swap(prev, cur);
    }
    return prev[b.size()];
}

double d_sequence(std::string_view a, std::string_view b) {
    size_t m = std::max(a.size(), b.size());
    if (m == 0) return 0.0;
    return static_cast<double>(levenshtein(a, b)) / static_cast<double>(m);
}

static double log_dist(double a, double b, double scale) {
    return std::min(1.0, std::fabs(std::log1p(a) - std::log1p(b)) / scale);
}

FeatureDistances feature_distances(const Pattern& i, const Pattern& j, const SimilarityParams& p) {
    FeatureDistances d;
    d.type = i.type == j.type ? 0.0 : 1.0;
    d.length = log_dist(i.length, j.length, p.s_L);
    d.motif = d_sequence(i.motif, j.motif);
    d.context = d_sequence(i.context, j.context);
    d.hp = log_dist(i.hp_length, j.hp_length, p.s_h);
    // ASSUMPTION: half unit mismatch, half copy-number difference.
    d.repeat = 0.5 * (i.period != j.period ? 1.0 : 0.0) + 0.5 * log_dist(i.copies, j.copies, p.s_r);
    return d;
}

double pattern_distance(const Pattern& i, const Pattern& j, const SimilarityParams& p) {
    FeatureDistances d = feature_distances(i, j, p);
    double num = p.w_type * d.type + p.w_length * d.length + p.w_motif * d.motif +
                 p.w_context * d.context + p.w_hp * d.hp + p.w_repeat * d.repeat;
    double den = p.w_type + p.w_length + p.w_motif + p.w_context + p.w_hp + p.w_repeat;
    return num / den;
}

static double percentile(std::vector<double>& v, double q) {
    if (v.empty()) return 1.0;
    size_t k = static_cast<size_t>(q * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<long>(k), v.end());
    return v[k];
}

SimilarityParams estimate_scales(const std::vector<Pattern>& tr, SimilarityParams base, size_t n_pairs,
                                 uint64_t seed) {
    if (tr.size() < 2) return base;
    std::vector<double> dl, dh, dr;
    uint64_t s = seed;
    for (size_t k = 0; k < n_pairs; ++k) {
        s = mix64(s);
        const Pattern& a = tr[s % tr.size()];
        s = mix64(s);
        const Pattern& b = tr[s % tr.size()];
        auto push = [](std::vector<double>& v, double x, double y) {
            double d = std::fabs(std::log1p(x) - std::log1p(y));
            if (d > 0) v.push_back(d);
        };
        push(dl, a.length, b.length);
        push(dh, a.hp_length, b.hp_length);
        push(dr, a.copies, b.copies);
    }
    base.s_L = std::max(1e-6, percentile(dl, 0.95));
    base.s_h = std::max(1e-6, percentile(dh, 0.95));
    base.s_r = std::max(1e-6, percentile(dr, 0.95));
    return base;
}

}  // namespace nindels
