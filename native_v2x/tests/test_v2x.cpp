// V2.x unit tests: (1) pb_llr_n == pb_native_v2 == pb_reference (bit-exact) on random reads; (2) PB variants (caps, block
// consensus) on hand-checkable inputs; (3) model file parse + logit/veto decisions against hard-coded Python values.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include "dnav2/features.hpp"
#include "dnav2/model.hpp"
#include "dnav2/pb.hpp"
using namespace dnav2;
static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { std::printf("FAIL: %s\n", msg); fails++; } } while (0)
static uint64_t bits(double d) { uint64_t u; std::memcpy(&u, &d, 8); return u; }
int main() {
    std::mt19937_64 rng(12345);
    std::vector<double> dp;
    int n_pb = 0;
    for (int it = 0; it < 20000; it++) {
        int n = 1 + rng() % 48;
        uint8_t code[48], bq[48]; PbSlot slots[48];
        for (int i = 0; i < 48; i++) slots[i] = {255, 0, 0};
        int alt = rng() % 4;
        double vaf = (rng() % 100) / 100.0;
        for (int i = 0; i < n; i++) {
            code[i] = (rng() % 100) < vaf * 100 ? alt : (alt + 1) % 4;
            bq[i] = 13 + rng() % 30;
            slots[i] = {code[i], bq[i], 1};
        }
        double a = pb_native_v2(code, bq, n, alt), b = pb_llr_n(code, bq, n, alt, dp), c = pb_reference(slots, alt);
        CHECK(bits(a) == bits(b), "pb_llr_n != pb_native_v2");
        CHECK(bits(a) == bits(c), "pb_native_v2 != pb_reference");
        n_pb++;
    }
    std::printf("pb_llr_n vs pb_native_v2 vs pb_reference: %d random loci checked\n", n_pb);
    // PbStore variants: 130 admitted reads, all counted, ALT every 5th read; block structure known
    PbStore st; st.alt = 1; st.n_adm = 130;
    for (int i = 0; i < 130; i++) { st.code.push_back(i % 5 == 0 ? 1 : 0); st.bq.push_back(30); }
    st.cnt_at.resize(131); for (int i = 0; i <= 130; i++) st.cnt_at[i] = i;
    double p48 = pb_store_eval(st, PBV_CAP48, dp), pall = pb_store_eval(st, PBV_ALL, dp);
    double direct48 = pb_llr_n(st.code.data(), st.bq.data(), 48, 1, dp);
    CHECK(bits(p48) == bits(direct48), "CAP48 != first 48 reads");
    double b0 = pb_llr_n(st.code.data(), st.bq.data(), 48, 1, dp), b1 = pb_llr_n(st.code.data() + 48, st.bq.data() + 48, 48, 1, dp);
    double mn = pb_store_eval(st, PBV_BLK_MIN, dp), mx = pb_store_eval(st, PBV_BLK_MAX, dp), md = pb_store_eval(st, PBV_BLK_MED, dp);
    CHECK(bits(mn) == bits(std::min(b0, b1)) && bits(mx) == bits(std::max(b0, b1)), "block min/max (2 full blocks; 34-read tail ignored)");
    CHECK(bits(md) == bits(0.5 * (std::min(b0, b1) + std::max(b0, b1))), "block median of 2");
    CHECK(pall > p48, "PB(all 130 reads) must exceed PB(48) at the same VAF");
    std::printf("PB variants: cap48=%.4f all=%.4f blk[min,med,max]=%.4f %.4f %.4f\n", p48, pall, mn, md, mx);
    // model file
    { std::ofstream f("/tmp/v2x_test.model");
      f << "name t\nkind logistic\nkmin 3\nintercept -0.3\nzthr 0.0\nfeat depth 5.0 0.5 1.5\nfeat vaf 0.2 0.1 -2.0\n"; }
    Model m; std::string e = m.load("/tmp/v2x_test.model");
    CHECK(e.empty(), "model load");
    double fv[NFEAT] = {0}; fv[0] = 200.0; fv[5] = 0.35;
    double z = m.logit(fv, false);
    const double expect = -0.3 + 1.5 * ((std::log1p(200.0) - 5.0) / 0.5) + (-2.0) * ((0.35 - 0.2) / 0.1);
    CHECK(std::fabs(z - expect) < 1e-12, "logit value");
    int64_t ev = 0;
    CHECK(m.decide(fv, false, 20.0, nullptr, dp, &ev) == (z >= 0.0), "logistic decision");
    { std::ofstream f("/tmp/v2x_test2.model"); f << "name s\nkind stump\nveto vaf le 0.4\n"; }
    Model s; CHECK(s.load("/tmp/v2x_test2.model").empty(), "stump load");
    CHECK(!s.decide(fv, false, 20.0, nullptr, dp, &ev), "veto vaf<=0.4 must veto vaf 0.35");
    fv[5] = 0.5; CHECK(s.decide(fv, false, 20.0, nullptr, dp, &ev), "no veto at vaf 0.5 and Binomial>=7");
    CHECK(!s.decide(fv, false, 6.9, nullptr, dp, &ev), "Binomial<7 not called");
    std::printf("logit=%.15g expected=%.15g\n", z, expect);
    std::printf(fails ? "V2X UNIT TESTS FAILED (%d)\n" : "V2X UNIT TESTS PASSED\n", fails);
    return fails ? 1 : 0;
}
