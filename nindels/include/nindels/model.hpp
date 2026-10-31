// Spec §12: z_i = beta_0 + sum beta_j X_ij, P_final = 1 / (1 + exp(-z_i)).
// Logistic regression fitted by IRLS (Newton) with an L2 penalty on
// standardized features. Missing values are handled explicitly: a feature that
// was missing in training gets a 0/1 missingness indicator, and its value is
// replaced by the training mean only together with that indicator.
#pragma once

#include <string>
#include <vector>

#include "nindels/types.hpp"

namespace nindels {

struct NamedFeature {
    std::string name;
    Opt value;
};
using FeatureMap = std::vector<NamedFeature>;

class LogisticModel {
public:
    // rows[i][j] is feature names[j] of example i; labels are 0/1.
    static LogisticModel fit(const std::vector<std::string>& names, const std::vector<std::vector<Opt>>& rows,
                             const std::vector<int>& labels, double lambda = 1e-2, int max_iter = 100);

    double predict(const std::vector<Opt>& x) const;  // P_final
    double predict(const FeatureMap& f) const;        // features looked up by name
    // Number of features outside the training [0.5%, 99.5%] range, plus
    // features that are missing although they were never missing in training.
    int out_of_range(const std::vector<Opt>& x) const;
    int out_of_range(const FeatureMap& f) const;
    std::vector<Opt> select(const FeatureMap& f) const;  // throws if a model feature is absent

    void save(const std::string& path) const;
    static LogisticModel load(const std::string& path);

    const std::vector<std::string>& names() const { return names_; }
    const std::vector<double>& coefficients() const { return beta_; }
    double intercept() const { return beta0_; }

private:
    std::vector<double> design(const std::vector<Opt>& x) const;

    std::vector<std::string> names_;
    std::vector<int> has_indicator_;
    std::vector<double> mean_, sd_, lo_, hi_;
    std::vector<double> beta_;  // value coefficients, then indicator coefficients
    double beta0_ = 0;
    double lambda_ = 0;
};

}  // namespace nindels
