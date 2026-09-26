// Numerical primitives that reproduce, bit for bit, the NumPy/SciPy functions the frozen V1
// pipeline calls.  Everything here must be compiled with -ffp-contract=off (no FMA fusion) or
// the results stop being identical to V1.  See tests/test_kernels.py for the bit-equality tests.
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>

namespace dnav2 {

// ---- scipy.special.gammaln == cephes lgam (xsf::cephes::lgam) -------------------------------
namespace detail {
inline double polevl(double x, const double* c, int n) {
    double a = c[0];
    for (int i = 1; i <= n; i++) a = a * x + c[i];
    return a;
}
inline double p1evl(double x, const double* c, int n) {
    double a = x + c[0];
    for (int i = 1; i < n; i++) a = a * x + c[i];
    return a;
}
constexpr double LGAM_A[] = {8.11614167470508450300E-4, -5.95061904284301438324E-4,
                             7.93650340457716943945E-4, -2.77777777730099687205E-3,
                             8.33333333333331927722E-2};
constexpr double LGAM_B[] = {-1.37825152569120859100E3, -3.88016315134637840924E4,
                             -3.31612992738871184744E5, -1.16237097492762307383E6,
                             -1.72173700820839662146E6, -8.53555664245765465627E5};
constexpr double LGAM_C[] = {-3.51815701436523470549E2, -1.70642106651881159223E4,
                             -2.20528590553854454839E5, -1.13933444367982507207E6,
                             -2.53252307177582951285E6, -2.01889141433532773231E6};
constexpr double LS2PI = 0.91893853320467274178;
constexpr double MAXLGM = 2.556348e305;
}  // namespace detail

// Only the x >= 1 domain is needed (arguments are count+1); anything else returns NaN so a
// misuse cannot silently produce a plausible number.
inline double gammaln(double x) {
    using namespace detail;
    if (!(x >= 1.0) || !std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
    if (x < 13.0) {
        double z = 1.0, p = 0.0, u = x;
        while (u >= 3.0) { p -= 1.0; u = x + p; z *= u; }
        while (u < 2.0) { z /= u; p += 1.0; u = x + p; }
        if (u == 2.0) return std::log(z);
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, LGAM_B, 5) / p1evl(x, LGAM_C, 6);
        return std::log(z) + p;
    }
    if (x > MAXLGM) return std::numeric_limits<double>::infinity();
    double q = (x - 0.5) * std::log(x) - x + LS2PI;
    if (x > 1.0e8) return q;
    double p = 1.0 / (x * x);
    if (x >= 1000.0)
        q += ((7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p +
              0.0833333333333333333333) / x;
    else
        q += polevl(p, LGAM_A, 4) / x;
    return q;
}

// scipy.special.log1p == cephes log1p (NOT glibc log1p: they differ in the last bit for many
// arguments, e.g. log1p(-p_het) used by the Binomial arm).  NumPy's np.log1p (used by the PB
// arm) is glibc's, so PB keeps std::log1p.
inline double log1p_cephes(double x) {
    static const double LP[] = {4.5270000862445199635215E-5, 4.9854102823193375972212E-1,
                                6.5787325942061044846969E0, 2.9911919328553073277375E1,
                                6.0949667980987787057556E1, 5.7112963590585538103336E1,
                                2.0039553499201281259648E1};
    static const double LQ[] = {1.5062909083469192043167E1, 8.3047565967967209469434E1,
                                2.2176239823732856465394E2, 3.0909872225312059774938E2,
                                2.1642788614495947685003E2, 6.0118660497603843919306E1};
    double z = 1.0 + x;
    if (z < 0.70710678118654752440 || z > 1.41421356237309504880) return std::log(z);
    z = x * x;
    z = -0.5 * z + x * (z * detail::polevl(x, LP, 6) / detail::p1evl(x, LQ, 6));
    return x + z;
}

// scipy.special.xlogy / xlog1py for finite arguments
inline double xlogy(double x, double y) { return x == 0.0 ? 0.0 : x * std::log(y); }
inline double xlog1py(double x, double y) { return x == 0.0 ? 0.0 : x * log1p_cephes(y); }

// numpy.logaddexp (npy_logaddexp)
inline double logaddexp(double x, double y) {
    if (x == y) return x + 0.693147180559945309417232121458176568;
    const double tmp = x - y;
    if (tmp > 0) return x + std::log1p(std::exp(-tmp));
    if (tmp <= 0) return y + std::log1p(std::exp(tmp));
    return tmp;  // NaN
}

// numpy.clip on a scalar
inline double clip(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

}  // namespace dnav2
