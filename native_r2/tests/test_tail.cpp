// Unit tests for the R2 tail: histogram quantiles, OOD envelope, pattern layer (vs brute force + fuzz), eps-aware binomial, model tail gate.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "dnav2/binomial.hpp"
#include "dnav2/model.hpp"
#include "calibrate.hpp"
#include "tail.hpp"

using namespace dnav2;
static int failures = 0, checks = 0;
#define CHECK(cond) do { checks++; if (!(cond)) { failures++; std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #cond); } } while (0)
static bool near(double a, double b, double tol = 1e-12) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b)); }

static void test_hist() {
    std::vector<uint64_t> h(10, 0);
    CHECK(hist_quantile(h, 0.5) == 0.0);
    h[7] = 5;
    CHECK(hist_quantile(h, 0.0) == 7.0 || hist_quantile(h, 0.0) == 0.0);   // q=0 -> first bin reaching 0 mass
    CHECK(hist_quantile(h, 0.5) == 7.0);
    CHECK(hist_quantile(h, 1.0) == 7.0);
    h[2] = 5;
    CHECK(hist_quantile(h, 0.5) == 2.0);
    CHECK(hist_quantile(h, 0.51) == 7.0);
}

static void test_ood() {
    Calibration c; OodEnvelope e;
    CHECK(!ood_reason(c, e).empty());                 // invalid calibration is always a reason
    c.valid = true; c.eps_hat = 0.01; c.mq_p90 = 70; c.bq_p50 = 37; c.depth_med = 30;
    CHECK(ood_reason(c, e).empty());                  // no envelope -> accept
    e.set = true; e.eps_lo = 0.002; e.eps_hi = 0.03; e.mq_lo = 50; e.mq_hi = 80; e.bq_lo = 20; e.bq_hi = 45; e.depth_lo = 4; e.depth_hi = 600;
    CHECK(ood_reason(c, e).empty());
    Calibration d = c; d.eps_hat = 0.05; CHECK(ood_reason(d, e).find("eps") == 0);
    d = c; d.mq_p90 = 40; CHECK(ood_reason(d, e).find("mapq_p90") == 0);
    d = c; d.bq_p50 = 10; CHECK(ood_reason(d, e).find("bq_p50") == 0);
    d = c; d.depth_med = 2; CHECK(ood_reason(d, e).find("depth_med") == 0);
    d = c; d.depth_med = 600; CHECK(ood_reason(d, e).empty());   // bounds inclusive
}

static void brute(const std::vector<CandPoint>& c, int64_t lo_span, int64_t hi_span, const Calibration& cal, std::vector<PatternOut>& out) {
    out.assign(c.size(), PatternOut{0, 0, 0});
    for (size_t j = 0; j < c.size(); j++) {
        double n1 = 0, n2 = 0, low = 0;
        for (size_t k = 0; k < c.size(); k++) {
            if (k == j) continue;
            const int64_t d = std::llabs(c[k].pos - c[j].pos);
            if (d <= 100) n1++;
            if (d <= 5000) { n2++; if (c[k].vaf < LOWVAF_CUT) low++; }
        }
        const int64_t p = c[j].pos;
        const int64_t w1 = std::min(p + 101, hi_span) - std::max(p - 100, lo_span) - 1;
        const int64_t w2 = std::min(p + 5001, hi_span) - std::max(p - 5000, lo_span) - 1;
        out[j].dens100 = n1 / std::max(0.5, cal.cand_per_bp * double(std::max<int64_t>(w1, 1)));
        out[j].dens5k = n2 / std::max(0.5, cal.cand_per_bp * double(std::max<int64_t>(w2, 1)));
        out[j].lowvaf5k = (low + 0.4) / (n2 + 2.0);
    }
}

static void test_pattern() {
    Calibration cal; cal.valid = true; cal.cand_per_bp = 0.001;
    std::vector<PatternOut> o;
    pattern_features({}, 0, 1000, cal, o);
    CHECK(o.empty());
    pattern_features({{500, 0.5}}, 0, 1000, cal, o);
    CHECK(o.size() == 1 && o[0].dens100 == 0.0 && o[0].dens5k == 0.0);
    CHECK(near(o[0].lowvaf5k, 0.4 / 2.0));            // shrunk prior only
    // fuzz against the brute-force definition
    std::mt19937_64 rng(20261027);
    for (int trial = 0; trial < 200; trial++) {
        const int n = static_cast<int>(rng() % 60);
        const int64_t span_lo = static_cast<int64_t>(rng() % 1000), span_hi = span_lo + 20000 + static_cast<int64_t>(rng() % 100000);
        std::vector<CandPoint> c;
        int64_t p = span_lo;
        for (int i = 0; i < n; i++) {
            p += 1 + static_cast<int64_t>(rng() % 3000);
            if (p >= span_hi) break;
            c.push_back({p, double(rng() % 1000) / 1000.0});
        }
        cal.cand_per_bp = double(rng() % 100) / 20000.0;   // includes 0 (floor of the expectation must prevent division by zero)
        std::vector<PatternOut> a, b;
        pattern_features(c, span_lo, span_hi, cal, a);
        brute(c, span_lo, span_hi, cal, b);
        CHECK(a.size() == b.size());
        for (size_t j = 0; j < a.size() && j < b.size(); j++) {
            CHECK(near(a[j].dens100, b[j].dens100)); CHECK(near(a[j].dens5k, b[j].dens5k)); CHECK(near(a[j].lowvaf5k, b[j].lowvaf5k));
            CHECK(std::isfinite(a[j].dens100) && std::isfinite(a[j].dens5k) && std::isfinite(a[j].lowvaf5k));
        }
    }
}

static void test_binomial_eps() {
    // eps-aware scoring at the frozen eps is bit-identical to the frozen function
    std::mt19937_64 rng(7);
    for (int i = 0; i < 20000; i++) {
        double c4[4] = {double(rng() % 40), double(rng() % 40), double(rng() % 40), double(rng() % 40)};
        const int ref = static_cast<int>(rng() % 5);
        const BinomialScore a = binomial_score(c4, ref), b = binomial_score_eps(c4, ref, FROZEN_EPSILON);
        CHECK(a.llr == b.llr && a.k == b.k && a.n == b.n && a.alt_index == b.alt_index);
    }
    double c4[4] = {0, 0, 0, 0};
    CHECK(binomial_score_eps(c4, 1, 0.005).llr == 0.0);   // empty column
    double d4[4] = {20, 0, 10, 0};
    CHECK(binomial_score_eps(d4, 1, 0.002).llr > binomial_score_eps(d4, 1, 0.03).llr);   // lower error rate -> more evidence
}

static void test_model_tail() {
    const char* path = "/tmp/dnav2_test_tail.model";
    {
        std::ofstream f(path);
        f << "name T\nkind logistic\nkmin 3\nintercept 0.5\nzthr 0.0\nfeat vaf 0.4 0.3 2.0\n"
             "tail_feat depth_rel 1.0 0.2 -1.5\ngate_depth_rel 1.3\nood_eps 0.002 0.03\n";
    }
    Model m;
    CHECK(m.load(path).empty());
    CHECK(m.ood.set && m.ood.eps_hi == 0.03);
    double f[NFEAT];
    for (int j = 0; j < NFEAT; j++) f[j] = 0.0;
    f[5] = 0.5;                                        // vaf
    f[FI_DEPTH_REL] = 1.0;
    const double base = m.logit(f, false);
    CHECK(near(base, 0.5 + 2.0 * (0.5 - 0.4) / 0.3));
    CHECK(!m.tail_gate(f));
    f[FI_DEPTH_REL] = 1.5;
    CHECK(m.tail_gate(f));
    const double with_tail = m.logit(f, false);
    CHECK(near(with_tail, base + (-1.5) * ((1.5 - 1.0) / 0.2)));
    m.tail_strength = 0.0;
    CHECK(near(m.logit(f, false), base));              // strength 0 disables the tail
    m.tail_strength = 1.0;
    f[FI_DEPTH_REL] = 1.0; f[FI_FORCED] = 1.0;         // researcher-forced region opens the gate
    CHECK(m.tail_gate(f));
    // feat2: cap before log1p
    { std::ofstream g(path); g << "kind logistic\nintercept 0\nzthr 0\nfeat2 dens5k 1.0 2.0 1.5 1 4\ntail_feat2 depth_rel 0 1 0.5 0 2\ngate_depth_rel 1.0\n"; }
    Model m2; CHECK(m2.load(path).empty());
    double g[NFEAT]; for (int j = 0; j < NFEAT; j++) g[j] = 0.0;
    g[FI_DENS5K] = 100.0; g[FI_DEPTH_REL] = 5.0;
    CHECK(near(m2.logit(g, false), 1.5 * ((std::log1p(4.0) - 1.0) / 2.0) + 0.5 * ((2.0 - 0.0) / 1.0)));
    { std::ofstream g2(path); g2 << "kind logistic\nfeat2 vaf 0 0 1 0 inf\n"; }
    Model m3; CHECK(!m3.load(path).empty());               // zero scale rejected
    // unknown keys and unknown tail features must be rejected
    { std::ofstream g(path); g << "kind logistic\ntail_feat nonexistent 0 1 1\n"; }
    Model bad; CHECK(!bad.load(path).empty());
    std::remove(path);
}

static void test_calibrate_errors() {
    Calibration c;
    std::vector<Window> none;
    CHECK(!calibrate("/nonexistent.bam", "/nonexistent.fa", "chr1", none, nullptr, 10, 10, 2, c).empty());
    CHECK(!c.valid);
    std::vector<Window> w = {{0}, {64}};
    CHECK(!calibrate("/nonexistent.bam", "/nonexistent.fa", "chr1", w, nullptr, 10, 10, 2, c).empty());
}

int main() {
    test_hist(); test_ood(); test_pattern(); test_binomial_eps(); test_model_tail(); test_calibrate_errors();
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
