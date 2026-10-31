#include "nindels/model.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nindels {

namespace {

double sigmoid(double z) {
    return z >= 0 ? 1.0 / (1.0 + std::exp(-z)) : std::exp(z) / (1.0 + std::exp(z));
}

// Solve A x = b for symmetric positive definite A (Cholesky). A is n x n row-major.
std::vector<double> solve_spd(std::vector<double> A, std::vector<double> b, size_t n) {
    for (size_t j = 0; j < n; ++j) {
        double d = A[j * n + j];
        for (size_t k = 0; k < j; ++k) d -= A[j * n + k] * A[j * n + k];
        if (d <= 0) throw std::runtime_error("logistic fit: Hessian not positive definite");
        d = std::sqrt(d);
        A[j * n + j] = d;
        for (size_t i = j + 1; i < n; ++i) {
            double s = A[i * n + j];
            for (size_t k = 0; k < j; ++k) s -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / d;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        double s = b[i];
        for (size_t k = 0; k < i; ++k) s -= A[i * n + k] * b[k];
        b[i] = s / A[i * n + i];
    }
    for (size_t i = n; i-- > 0;) {
        double s = b[i];
        for (size_t k = i + 1; k < n; ++k) s -= A[k * n + i] * b[k];
        b[i] = s / A[i * n + i];
    }
    return b;
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return 0;
    size_t k = static_cast<size_t>(q * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<long>(k), v.end());
    return v[k];
}

}  // namespace

std::vector<double> LogisticModel::design(const std::vector<Opt>& x) const {
    // [1, standardized values..., indicators...]
    std::vector<double> d;
    d.reserve(1 + names_.size() * 2);
    d.push_back(1.0);
    for (size_t j = 0; j < names_.size(); ++j) d.push_back(x[j] ? (*x[j] - mean_[j]) / sd_[j] : 0.0);
    for (size_t j = 0; j < names_.size(); ++j)
        if (has_indicator_[j]) d.push_back(x[j] ? 0.0 : 1.0);
    return d;
}

LogisticModel LogisticModel::fit(const std::vector<std::string>& names, const std::vector<std::vector<Opt>>& rows,
                                 const std::vector<int>& labels, double lambda, int max_iter) {
    if (rows.size() != labels.size() || rows.empty()) throw std::invalid_argument("logistic fit: bad input");
    const size_t F = names.size();
    LogisticModel m;
    m.names_ = names;
    m.lambda_ = lambda;
    m.has_indicator_.assign(F, 0);
    m.mean_.assign(F, 0);
    m.sd_.assign(F, 1);
    m.lo_.assign(F, 0);
    m.hi_.assign(F, 0);
    for (size_t j = 0; j < F; ++j) {
        std::vector<double> v;
        for (const auto& r : rows) {
            if (r[j]) v.push_back(*r[j]);
            else m.has_indicator_[j] = 1;
        }
        if (v.empty()) throw std::invalid_argument("feature '" + names[j] + "' is missing in every training row");
        double mu = 0, var = 0;
        for (double x : v) mu += x;
        mu /= static_cast<double>(v.size());
        for (double x : v) var += (x - mu) * (x - mu);
        var /= std::max<double>(1.0, static_cast<double>(v.size() - 1));
        m.mean_[j] = mu;
        m.sd_[j] = var > 1e-12 ? std::sqrt(var) : 1.0;
        m.lo_[j] = quantile(v, 0.005);
        m.hi_[j] = quantile(v, 0.995);
    }

    std::vector<std::vector<double>> X;
    X.reserve(rows.size());
    for (const auto& r : rows) X.push_back(m.design(r));
    const size_t P = X.front().size();
    std::vector<double> beta(P, 0.0);

    for (int it = 0; it < max_iter; ++it) {
        std::vector<double> H(P * P, 0.0), g(P, 0.0);
        for (size_t i = 0; i < X.size(); ++i) {
            const auto& x = X[i];
            double z = 0;
            for (size_t a = 0; a < P; ++a) z += beta[a] * x[a];
            double p = sigmoid(z), w = std::max(p * (1 - p), 1e-12), resid = labels[i] - p;
            for (size_t a = 0; a < P; ++a) {
                if (x[a] == 0.0) continue;
                g[a] += resid * x[a];
                for (size_t b = 0; b <= a; ++b) H[a * P + b] += w * x[a] * x[b];
            }
        }
        for (size_t a = 0; a < P; ++a)
            for (size_t b = 0; b < a; ++b) H[b * P + a] = H[a * P + b];
        for (size_t a = 1; a < P; ++a) {  // no penalty on the intercept
            H[a * P + a] += lambda;
            g[a] -= lambda * beta[a];
        }
        H[0] += 1e-9;
        std::vector<double> step = solve_spd(H, g, P);
        double mx = 0;
        for (size_t a = 0; a < P; ++a) { beta[a] += step[a]; mx = std::max(mx, std::fabs(step[a])); }
        if (mx < 1e-8) break;
    }
    m.beta0_ = beta[0];
    m.beta_.assign(beta.begin() + 1, beta.end());
    return m;
}

double LogisticModel::predict(const std::vector<Opt>& x) const {
    std::vector<double> d = design(x);
    double z = beta0_;
    for (size_t a = 1; a < d.size(); ++a) z += beta_[a - 1] * d[a];
    return sigmoid(z);
}

std::vector<Opt> LogisticModel::select(const FeatureMap& f) const {
    std::vector<Opt> x(names_.size());
    for (size_t j = 0; j < names_.size(); ++j) {
        auto it = std::find_if(f.begin(), f.end(), [&](const NamedFeature& nf) { return nf.name == names_[j]; });
        if (it == f.end()) throw std::invalid_argument("model needs feature '" + names_[j] + "' which was not computed");
        x[j] = it->value;
    }
    return x;
}

double LogisticModel::predict(const FeatureMap& f) const { return predict(select(f)); }

int LogisticModel::out_of_range(const std::vector<Opt>& x) const {
    int n = 0;
    for (size_t j = 0; j < names_.size(); ++j) {
        if (!x[j]) n += has_indicator_[j] ? 0 : 1;
        else if (*x[j] < lo_[j] || *x[j] > hi_[j]) ++n;
    }
    return n;
}

int LogisticModel::out_of_range(const FeatureMap& f) const { return out_of_range(select(f)); }

void LogisticModel::save(const std::string& path) const {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write model " + path);
    o.precision(17);
    o << "#nindels-logistic\tv1\n#lambda\t" << lambda_ << "\n#intercept\t" << beta0_ << "\n";
    o << "#name\tmean\tsd\tlo\thi\tbeta\thas_indicator\tbeta_missing\n";
    size_t ind = names_.size();
    for (size_t j = 0; j < names_.size(); ++j) {
        o << names_[j] << '\t' << mean_[j] << '\t' << sd_[j] << '\t' << lo_[j] << '\t' << hi_[j] << '\t' << beta_[j]
          << '\t' << has_indicator_[j] << '\t' << (has_indicator_[j] ? beta_[ind++] : 0.0) << '\n';
    }
}

LogisticModel LogisticModel::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read model " + path);
    LogisticModel m;
    std::vector<double> bval, bind;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        if (line[0] == '#') {
            std::string tag;
            ss >> tag;
            if (tag == "#lambda") ss >> m.lambda_;
            else if (tag == "#intercept") ss >> m.beta0_;
            continue;
        }
        std::string name;
        double mu, sd, lo, hi, b, bm;
        int hi_ind;
        ss >> name >> mu >> sd >> lo >> hi >> b >> hi_ind >> bm;
        if (!ss) throw std::runtime_error("malformed model line: " + line);
        m.names_.push_back(name);
        m.mean_.push_back(mu);
        m.sd_.push_back(sd);
        m.lo_.push_back(lo);
        m.hi_.push_back(hi);
        m.has_indicator_.push_back(hi_ind);
        bval.push_back(b);
        if (hi_ind) bind.push_back(bm);
    }
    m.beta_ = bval;
    m.beta_.insert(m.beta_.end(), bind.begin(), bind.end());
    return m;
}

}  // namespace nindels
