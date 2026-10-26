// V2.x decision layer for candidate loci.  Parameters are read from a plain-text model file fitted OFFLINE on the
// development split only (analysis/export_model.py); the engine never learns anything.
//   kind logistic : z = intercept + sum_j coef_j * (x'_j - mean_j)/scale_j   (x' = log1p(x) for count-like features)
//                   call = z >= zthr ; optional PB exception band [zlo,zhi]: PB variant decides (frozen 10.5)
//   kind stump    : V2 decision (routed -> PB variant >= 10.5, else Binomial >= 7.0) minus hard vetoes
// Frozen constants (7.0 / router cutoff / 10.5) are used exactly where V2 uses them and never modified.
#pragma once
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include "binomial.hpp"
#include "features.hpp"
#include "tail.hpp"

namespace dnav2 {

struct Model {
    std::string name = "none", kind = "none";
    int kmin = 3;
    int routed_variant = PBV_CAP48;
    int band_variant = -1;
    double zthr = 0, zlo = 0, zhi = 0, intercept = 0;
    struct Feat { int idx; bool log1p; double mean, scale, coef; double cap = 1e300; };   // idx == -1: routed flag; cap: upper clip applied before log1p
    std::vector<Feat> feats;
    struct Veto { int idx; bool ge; double thr; };
    std::vector<Veto> vetos;
    // R2 tail: second-stage linear term applied only where the difficulty gate fires (see logit()).
    std::vector<Feat> tail_feats;
    double tail_bias = 0.0;               // constant added inside the tail term (only where the gate fires)
    double gate_depth_rel = 1e9, gate_lowvaf5k = 1e9, gate_dens5k = 1e9;
    double tail_strength = 1.0;           // runtime multiplier (researcher parameter); 0 disables the tail
    OodEnvelope ood;                      // training envelope for the out-of-distribution guard

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
            } else if (key == "feat2" || key == "tail_feat2") {
                // explicit transform: name mean scale coef log1p(0|1) cap   (cap "inf" = none)
                std::string n, capstr; Feat f{}; int lg = 0;
                in >> n >> f.mean >> f.scale >> f.coef >> lg >> capstr;
                f.idx = feat_index(n);
                if (f.idx == -2 || (f.idx < 0 && key == "tail_feat2")) return "unknown feature in model: " + n;
                f.log1p = lg != 0;
                f.cap = capstr == "inf" ? 1e300 : std::atof(capstr.c_str());
                if (!(f.scale > 0.0)) return "non-positive scale for feature: " + n;
                (key == "feat2" ? feats : tail_feats).push_back(f);
            } else if (key == "tail_bias") in >> tail_bias;
            else if (key == "tail_feat") {
                std::string n; Feat f{};
                in >> n >> f.mean >> f.scale >> f.coef;
                f.idx = feat_index(n);
                if (f.idx < 0) return "unknown tail feature in model: " + n;
                f.log1p = is_log1p(n);
                tail_feats.push_back(f);
            } else if (key == "gate_depth_rel") in >> gate_depth_rel;
            else if (key == "gate_lowvaf5k") in >> gate_lowvaf5k;
            else if (key == "gate_dens5k") in >> gate_dens5k;
            else if (key == "ood_eps") { in >> ood.eps_lo >> ood.eps_hi; ood.set = true; }
            else if (key == "ood_mq") { in >> ood.mq_lo >> ood.mq_hi; ood.set = true; }
            else if (key == "ood_bq") { in >> ood.bq_lo >> ood.bq_hi; ood.set = true; }
            else if (key == "ood_depth") { in >> ood.depth_lo >> ood.depth_hi; ood.set = true; }
            else if (key == "veto") {
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

    static double transform(const Feat& t, double x) {
        if (x > t.cap) x = t.cap;
        return t.log1p ? std::log1p(x) : x;
    }

    // Difficulty gate of the tail: researcher-forced region, collapsed-copy signature (depth well above the sample's own level),
    // dense candidate neighbourhood, or a noisy (low-VAF-rich) neighbourhood.
    bool tail_gate(const double f[NFEAT]) const {
        return f[FI_FORCED] > 0.5 || f[FI_DEPTH_REL] > gate_depth_rel || f[FI_DENS5K] > gate_dens5k || f[FI_LOWVAF5K] > gate_lowvaf5k;
    }

    double logit(const double f[NFEAT], bool routed) const {
        double z = intercept;
        for (const Feat& t : feats) {
            double x = t.idx < 0 ? (routed ? 1.0 : 0.0) : f[t.idx];
            x = transform(t, x);
            z += t.coef * ((x - t.mean) / t.scale);
        }
        if (!tail_feats.empty() && tail_strength != 0.0 && tail_gate(f)) {
            double adj = tail_bias;
            for (const Feat& t : tail_feats) {
                const double x = transform(t, f[t.idx]);
                adj += t.coef * ((x - t.mean) / t.scale);
            }
            z += tail_strength * adj;
        }
        return z;
    }

    // Decision for one candidate.  `store` may be null (PB input not kept: locus outside the PB band).
    // *pb_evals counts PB evaluations actually performed.
    bool decide(const double f[NFEAT], bool routed, double bllr, const PbStore* store, std::vector<double>& dp,
                int64_t* pb_evals, double* score = nullptr) const {
        if (kind == "logistic") {
            const double z = logit(f, routed);
            if (score) *score = z - zthr;
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
