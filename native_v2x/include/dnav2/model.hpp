// V2.x decision layer for candidate loci.  Parameters are read from a plain-text model file fitted OFFLINE on the
// development split only (analysis/export_model.py); the engine never learns anything.
//   kind logistic : z = intercept + sum_j coef_j * (x'_j - mean_j)/scale_j   (x' = log1p(x) for count-like features)
//                   call = z >= zthr ; optional PB exception band [zlo,zhi]: PB variant decides (frozen 10.5)
//   kind stump    : V2 decision (routed -> PB variant >= 10.5, else Binomial >= 7.0) minus hard vetoes
// Frozen constants (7.0 / router cutoff / 10.5) are used exactly where V2 uses them and never modified.
#pragma once
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include "binomial.hpp"
#include "features.hpp"

namespace dnav2 {

struct Model {
    std::string name = "none", kind = "none";
    int kmin = 3;
    int routed_variant = PBV_CAP48;
    int band_variant = -1;
    double zthr = 0, zlo = 0, zhi = 0, intercept = 0;
    struct Feat { int idx; bool log1p; double mean, scale, coef; };   // idx == -1: routed flag
    std::vector<Feat> feats;
    struct Veto { int idx; bool ge; double thr; };
    std::vector<Veto> vetos;

    bool active() const { return kind != "none"; }

    static bool is_log1p(const std::string& n) {
        return n == "depth" || n == "n" || n == "k" || n == "nref" || n == "k2" || n == "alt_lowbq";
    }
    static int feat_index(const std::string& n) {
        if (n == "routed_f") return -1;
        for (int i = 0; i < NFEAT; i++) if (n == FEAT_NAMES[i]) return i;
        return -2;
    }

    std::string load(const std::string& path) {
        std::ifstream in(path);
        if (!in) return "cannot open model " + path;
        std::string key;
        while (in >> key) {
            if (key == "name") in >> name;
            else if (key == "kind") in >> kind;
            else if (key == "kmin") in >> kmin;
            else if (key == "routed_variant") in >> routed_variant;
            else if (key == "band_variant") in >> band_variant;
            else if (key == "zthr") in >> zthr;
            else if (key == "zlo") in >> zlo;
            else if (key == "zhi") in >> zhi;
            else if (key == "intercept") in >> intercept;
            else if (key == "feat") {
                std::string n; Feat f{};
                in >> n >> f.mean >> f.scale >> f.coef;
                f.idx = feat_index(n);
                if (f.idx == -2) return "unknown feature in model: " + n;
                f.log1p = is_log1p(n);
                feats.push_back(f);
            } else if (key == "veto") {
                std::string n, dir; Veto v{};
                in >> n >> dir >> v.thr;
                v.idx = feat_index(n);
                if (v.idx < 0) return "unknown veto feature: " + n;
                v.ge = dir == "ge";
                vetos.push_back(v);
            } else return "unknown model key: " + key;
        }
        if (kind != "logistic" && kind != "stump") return "model kind must be logistic|stump";
        return "";
    }

    double logit(const double f[NFEAT], bool routed) const {
        double z = intercept;
        for (const Feat& t : feats) {
            double x = t.idx < 0 ? (routed ? 1.0 : 0.0) : f[t.idx];
            if (t.log1p) x = std::log1p(x);
            z += t.coef * ((x - t.mean) / t.scale);
        }
        return z;
    }

    // Decision for one candidate.  `store` may be null (PB input not kept: locus outside the PB band).
    // *pb_evals counts PB evaluations actually performed.
    bool decide(const double f[NFEAT], bool routed, double bllr, const PbStore* store, std::vector<double>& dp,
                int64_t* pb_evals) const {
        if (kind == "logistic") {
            const double z = logit(f, routed);
            bool call = z >= zthr;
            if (band_variant >= 0 && store && z >= zlo && z <= zhi) {
                call = pb_store_eval(*store, band_variant, dp) >= FROZEN_PB_THRESHOLD;
                (*pb_evals)++;
            }
            return call;
        }
        bool call;
        if (routed) {
            call = store ? pb_store_eval(*store, routed_variant, dp) >= FROZEN_PB_THRESHOLD : false;
            (*pb_evals)++;
        } else call = bllr >= FROZEN_BINOMIAL_THRESHOLD;
        for (const Veto& v : vetos) if (call && (v.ge ? f[v.idx] >= v.thr : f[v.idx] <= v.thr)) call = false;
        return call;
    }
};

}  // namespace dnav2
