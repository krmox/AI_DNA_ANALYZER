// C ABI over the kernels, used only by the test-suite (ctypes) and the oracle comparisons.
#include <cstdint>
#include <cstring>
#include "dnav2/binomial.hpp"
#include "dnav2/methodc.hpp"
#include "dnav2/pb.hpp"

extern "C" {
void v2_log1p_cephes(const double* x, double* out, int64_t n) {
    for (int64_t i = 0; i < n; i++) out[i] = dnav2::log1p_cephes(x[i]);
}
void v2_gammaln(const double* x, double* out, int64_t n) {
    for (int64_t i = 0; i < n; i++) out[i] = dnav2::gammaln(x[i]);
}
// counts: [n,10] float64 (pileup_counts layout).  out_*: [n]
void v2_binomial(const double* counts, int64_t n, double* llr, int32_t* alt, double* k, double* nn) {
    for (int64_t i = 0; i < n; i++) {
        auto s = dnav2::binomial_score(counts + i * 10, static_cast<int>(counts[i * 10 + 9]));
        llr[i] = s.llr; alt[i] = s.alt_index; k[i] = s.k; nn[i] = s.n;
    }
}
// reads: [n,48,8] uint8 V1 tensor rows.  mode 0 = reference, 1 = native v2 (compact)
void v2_pb(const uint8_t* reads, const int32_t* alt, int64_t n, int mode, double* out) {
    for (int64_t i = 0; i < n; i++) {
        const uint8_t* r = reads + i * 48 * 8;
        if (mode == 0) {
            dnav2::PbSlot s[48];
            for (int j = 0; j < 48; j++) s[j] = {r[j * 8], r[j * 8 + 1], (uint8_t)((r[j * 8 + 6] >> 1) & 1)};
            out[i] = dnav2::pb_reference(s, alt[i]);
        } else {
            uint8_t code[48], bq[48]; int m = 0;
            for (int j = 0; j < 48; j++)
                if ((r[j * 8 + 6] >> 1) & 1) { code[m] = r[j * 8]; bq[m] = r[j * 8 + 1]; m++; }
            out[i] = dnav2::pb_native_v2(code, bq, m, alt[i]);
        }
    }
}
// Method C: outputs gt index (1='0/1', 2='1/1'), gq int, forced flag
void v2_methodc(const double* k, const double* n, int64_t cnt, int32_t* gt, int32_t* gq, int32_t* forced) {
    for (int64_t i = 0; i < cnt; i++) {
        auto r = dnav2::method_c(k[i], n[i]);
        gt[i] = r.gt; gq[i] = r.gq; forced[i] = r.forced;
    }
}
}
