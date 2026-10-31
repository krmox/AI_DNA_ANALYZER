#include "nindels/str_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nindels {

int str_min_units(int period) { return period == 1 ? 4 : 3; }

std::vector<StrRun> find_str_runs(std::string_view seq) {
    const int64_t n = static_cast<int64_t>(seq.size());
    std::vector<StrRun> all;
    for (int p = 1; p <= StrErrorModel::kMaxPeriod; ++p) {
        int64_t j = p;
        while (j < n) {
            if (seq[j] != seq[j - p] || seq[j] == 'N') { ++j; continue; }
            int64_t k = j;
            while (k < n && seq[k] == seq[k - p] && seq[k] != 'N') ++k;
            StrRun r{j - p, k, p};
            // a period-p run must not be a shorter-period run in disguise
            bool reducible = false;
            for (int q = 1; q < p && !reducible; ++q)
                if (p % q == 0) {
                    bool all_eq = true;
                    for (int64_t t = r.start + q; t < r.start + p && all_eq; ++t) all_eq = seq[t] == seq[t - q];
                    reducible = all_eq;
                }
            if (!reducible && r.units() >= str_min_units(p)) all.push_back(r);
            j = k + 1;
        }
    }
    // keep, per position, the longest covering run (ties: smaller period)
    std::sort(all.begin(), all.end(), [](const StrRun& a, const StrRun& b) {
        int64_t la = a.end - a.start, lb = b.end - b.start;
        return la != lb ? la > lb : a.period < b.period;
    });
    std::vector<char> taken(static_cast<size_t>(n), 0);
    std::vector<StrRun> out;
    for (const StrRun& r : all) {
        bool overlap = false;
        for (int64_t t = r.start; t < r.end && !overlap; ++t) overlap = taken[static_cast<size_t>(t)];
        if (overlap) continue;
        for (int64_t t = r.start; t < r.end; ++t) taken[static_cast<size_t>(t)] = 1;
        out.push_back(r);
    }
    std::sort(out.begin(), out.end(), [](const StrRun& a, const StrRun& b) { return a.start < b.start; });
    return out;
}

namespace {

// Three-component binomial mixture EM; returns eps.
double fit_eps(const std::vector<std::pair<int, int>>& ns) {
    double eps = 0.01, pi[3] = {0.9, 0.07, 0.03};
    for (int it = 0; it < 200; ++it) {
        double num = 0, den = 0, w[3] = {0, 0, 0};
        const double p[3] = {eps, 0.5, 1 - eps};
        for (auto [n, s] : ns) {
            double lp[3], mx = -1e300;
            for (int z = 0; z < 3; ++z) {
                lp[z] = std::log(pi[z]) + s * std::log(p[z]) + (n - s) * std::log1p(-p[z]);
                mx = std::max(mx, lp[z]);
            }
            double r[3], tot = 0;
            for (int z = 0; z < 3; ++z) tot += (r[z] = std::exp(lp[z] - mx));
            for (int z = 0; z < 3; ++z) { r[z] /= tot; w[z] += r[z]; }
            num += r[0] * s + r[2] * (n - s);
            den += (r[0] + r[2]) * n;
        }
        const double total = w[0] + w[1] + w[2];
        for (int z = 0; z < 3; ++z) pi[z] = std::max(1e-6, w[z] / total);
        double next = den > 0 ? std::clamp(num / den, 1e-6, 0.3) : eps;
        if (std::fabs(next - eps) < 1e-9) { eps = next; break; }
        eps = next;
    }
    return eps;
}

std::string synthetic_unit(int p) {
    static const char* units[] = {"", "A", "AC", "AAC", "AAAC", "AAAAC", "AAAAAC"};
    return units[p];
}

}  // namespace

void StrErrorModel::fit(const std::vector<StrLocusCount>& counts, int min_loci) {
    std::vector<std::pair<int, int>> by[kMaxPeriod + 1][kMaxUnits + 1];
    for (const auto& c : counts)
        if (c.period >= 1 && c.period <= kMaxPeriod && c.n >= 5) by[c.period][clamp_u(c.units)].emplace_back(c.n, c.s);
    for (int p = 1; p <= kMaxPeriod; ++p) {
        double last = -1;
        for (int u = 0; u <= kMaxUnits; ++u) {
            loci_[p][u] = static_cast<int>(by[p][u].size());
            if (u < str_min_units(p)) { eps_[p][u] = -1; continue; }
            if (loci_[p][u] >= min_loci) last = fit_eps(by[p][u]);
            eps_[p][u] = last;  // inherits from (p, u-1) when sparse; -1 if nothing below
        }
    }
    fitted_ = true;
}

void StrErrorModel::calibrate(const PairHmmParams& hmm) {
    const double base = std::pow(10.0, -hmm.gap_open_q / 10.0);
    const std::string left = "GATCCTAGGCATTGACCGTA", right = "TGCAGGTCATGCATTGCAGC";
    for (int p = 1; p <= kMaxPeriod; ++p)
        for (int u = 0; u <= kMaxUnits; ++u) {
            gop_[p][u] = -1;
            if (eps_[p][u] <= 0) continue;
            std::string unit = synthetic_unit(p), rep;
            for (int k = 0; k < std::max(u, 2); ++k) rep += unit;
            const std::string hap = left + rep + right;
            const std::string clean = hap, slip = left + rep.substr(static_cast<size_t>(p)) + right;
            std::vector<uint8_t> q(hap.size(), 30);
            const double target = eps_[p][u] / (2.0 * (1.0 - eps_[p][u]));
            double g = base;
            for (int it = 0; it < 8; ++it) {
                std::vector<double> prof(hap.size(), base);
                for (size_t j = left.size(); j <= left.size() + rep.size() && j < prof.size(); ++j) prof[j] = g;
                double lc = pairhmm_loglik(clean, q.data(), hap, hmm, 0, -1, prof.data());
                double ls = pairhmm_loglik(slip, q.data(), hap, hmm, 0, -1, prof.data());
                double R = std::exp(ls - lc);
                double next = std::clamp(g * target / std::max(R, 1e-300), base, 0.1);  // 0.1 cap: ARBITRARY
                if (std::fabs(next - g) < 1e-3 * g) { g = next; break; }
                g = next;
            }
            gop_[p][u] = g;
        }
}

std::vector<double> StrErrorModel::gap_open_profile(std::string_view hap, double base_gop) const {
    std::vector<double> prof(hap.size(), base_gop);
    for (const StrRun& r : find_str_runs(hap)) {
        const double g = gop_[r.period][clamp_u(r.units())];
        if (g <= 0) continue;
        // gaps may open at any base of the run and right after its last base
        for (int64_t j = r.start; j <= r.end && j < static_cast<int64_t>(prof.size()); ++j)
            prof[static_cast<size_t>(j)] = std::max(prof[static_cast<size_t>(j)], g);
    }
    return prof;
}

void StrErrorModel::save(const std::string& path) const {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(10);
    o << "#period\tunits\tloci\teps\tgap_open\n";
    for (int p = 1; p <= kMaxPeriod; ++p)
        for (int u = str_min_units(p); u <= kMaxUnits; ++u)
            o << p << '\t' << u << '\t' << loci_[p][u] << '\t' << eps_[p][u] << '\t' << gop_[p][u] << '\n';
}

StrErrorModel StrErrorModel::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    StrErrorModel m;
    for (auto& row : m.eps_) std::fill(std::begin(row), std::end(row), -1.0);
    for (auto& row : m.gop_) std::fill(std::begin(row), std::end(row), -1.0);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int p, u, n;
        double e, g;
        if (!(ss >> p >> u >> n >> e >> g)) throw std::runtime_error("bad STR model line: " + line);
        if (p < 1 || p > kMaxPeriod || u < 0 || u > kMaxUnits) continue;
        m.loci_[p][u] = n;
        m.eps_[p][u] = e;
        m.gop_[p][u] = g;
    }
    m.fitted_ = true;
    return m;
}

}  // namespace nindels
