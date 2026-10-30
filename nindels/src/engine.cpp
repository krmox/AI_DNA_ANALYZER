#include "nindels/engine.hpp"

#include <cmath>

namespace nindels {

const char* decision_name(Decision d) {
    switch (d) {
        case Decision::Pass: return "PASS";
        case Decision::LowProb: return "LowProb";
        case Decision::Uncertain: return "UNCERTAIN";
        case Decision::NoModel: return "NoModel";
    }
    return "?";
}

const std::vector<std::string>& all_feature_names() {
    static const std::vector<std::string> names = {
        // read evidence E_i (spec §2.1-§4)
        "f", "l_depth", "alt_fwd_frac", "strand_bias", "readpos_alt", "readpos_diff", "clip_alt", "clip_all",
        "mismatch_alt", "other_frac", "split_frac", "mapq_alt", "mq0_frac", "mqlow_frac", "bq_alt",
        // diploid genotype from counts
        "count_qual",
        // pattern memory (spec §7, §7.1)
        "p_context", "n_eff_kish", "n_eff_loci", "s_max",
        // region complexity (spec §10)
        "c_r",
        // step 13 reassessment (spec §8, §11)
        "hmm_qual", "hmm_p_hom", "delta_ell_orig", "hmm_consistency", "q_min"};
    return names;
}

const std::vector<std::string>& reassessment_feature_names() {
    static const std::vector<std::string> names = {"hmm_qual", "hmm_p_hom", "delta_ell_orig", "hmm_consistency",
                                                   "q_min"};
    return names;
}

FeatureMap make_features(const CandidateResult& r) {
    const Evidence& e = r.evidence;
    FeatureMap f;
    f.reserve(all_feature_names().size());
    auto add = [&](const char* n, Opt v) { f.push_back({n, v}); };
    add("f", e.f);
    add("l_depth", e.l_depth);
    add("alt_fwd_frac", e.alt_fwd_frac);
    add("strand_bias", e.strand_bias);
    add("readpos_alt", e.readpos_alt);
    add("readpos_diff", e.readpos_diff);
    add("clip_alt", e.clip_alt);
    add("clip_all", e.clip_all);
    add("mismatch_alt", e.mismatch_alt);
    add("other_frac", e.other_frac);
    add("split_frac", e.split_frac);
    add("mapq_alt", e.mapq_alt);
    add("mq0_frac", e.mq0_frac);
    add("mqlow_frac", e.mqlow_frac);
    add("bq_alt", e.bq_alt);
    add("count_qual", r.count_genotype.qual);
    if (r.context) {
        add("p_context", r.context->p_context);
        add("n_eff_kish", r.context->n_eff_kish);
        add("n_eff_loci", r.context->n_eff_loci);
        add("s_max", r.context->s_max);
    } else {
        add("p_context", std::nullopt);
        add("n_eff_kish", std::nullopt);
        add("n_eff_loci", std::nullopt);
        add("s_max", std::nullopt);
    }
    add("c_r", r.complexity.C);
    if (r.hmm_genotype) {
        add("hmm_qual", r.hmm_genotype->qual);
        add("hmm_p_hom", r.hmm_genotype->posterior[2]);
    } else {
        add("hmm_qual", std::nullopt);
        add("hmm_p_hom", std::nullopt);
    }
    add("delta_ell_orig", r.delta_ell_original);
    add("hmm_consistency", r.consistency);
    add("q_min", r.nucleotide.q_min);
    return f;
}

static FeatureMap without_reassessment(FeatureMap f) {
    for (auto& nf : f)
        for (const auto& n : reassessment_feature_names())
            if (nf.name == n) nf.value = std::nullopt;
    return f;
}

CandidateResult process_indel(const CandidateInput& in, const RefSeq& ref, const PatternMemory* memory,
                              const LogisticModel* base_model, const LogisticModel* final_model,
                              const EngineConfig& cfg) {
    CandidateResult r;
    // 1. normalized representation (idempotent if already normalized)
    r.indel = left_normalize(ref, in.indel);
    // 2. pattern features P_i
    r.pattern = extract_pattern(ref, r.indel, cfg.pattern);
    // 3-5. read evidence E_i, f_i, L_depth
    r.evidence = summarize_evidence(in.obs, cfg.evidence);
    r.count_genotype = count_gl(r.evidence.n_alt, r.evidence.n_ref, cfg.count_gl_error, cfg.priors);
    r.p_evidence = 1.0 - r.count_genotype.posterior[0];
    // 6-9. memory retrieval, similarity, aggregated statistics
    if (memory) r.context = memory->query(r.pattern);
    // 10. region complexity
    r.complexity = region_complexity(r.pattern, r.evidence.instability, cfg.complexity);

    // 11. reassessment priority R_i. P_base is the stage-1 model when one is
    // supplied (fitted and calibrated on training/validation data); otherwise
    // the count-model posterior, which is a model probability but NOT a
    // validated calibrated one.
    FeatureMap stage1 = without_reassessment(make_features(r));
    r.p_base = base_model ? base_model->predict(stage1) : r.p_evidence;
    r.U = uncertainty(r.p_base);
    const ReassessParams& rp = cfg.reassess;
    if (r.context) {
        r.D = std::fabs(r.context->p_context - r.p_evidence);
        r.R = reassess_priority(r.U, r.complexity.C, *r.D, rp);
    } else {
        // D_i unavailable: renormalize over the remaining terms (explicit, not D = 0).
        double den = rp.a_U + rp.a_C;
        r.R = den > 0 ? (rp.a_U * r.U + rp.a_C * r.complexity.C) / den : 0.0;
    }
    // 12. tau_review(A). A only decides who is reassessed; it is not a feature.
    r.tau = tau_review(rp);

    // 13. local verification on H_ref / H_alt + nucleotide-level inspection
    r.reassessed = cfg.force_reassess || r.R >= r.tau;
    if (r.reassessed && !in.reads.empty()) {
        std::vector<size_t> used;
        std::vector<ReadLik> liks = read_likelihoods(ref, r.indel, in.reads, cfg.likelihood, &used);
        r.hmm_reads = static_cast<int>(liks.size());
        r.hmm_genotype = diploid_gl(liks, /*mapq_cap=*/true, cfg.priors);
        r.delta_ell_original = delta_ell_original(liks);
        r.consistency = hmm_consistency(liks);
        std::vector<const HmmRead*> rr;
        rr.reserve(used.size());
        for (size_t k : used) rr.push_back(&in.reads[k]);
        r.nucleotide = nucleotide_scores(rr, liks, r.evidence.f, cfg.nucleotide);
    }

    // 14. final feature vector
    r.features = make_features(r);

    // 15-16. P_final and decision
    if (final_model) {
        r.p_final = final_model->predict(r.features);
        r.ood = final_model->out_of_range(r.features);
    }
    if (!r.evidence.f || r.evidence.n_total < cfg.min_depth) r.decision = Decision::Uncertain;
    else if (!r.p_final) r.decision = Decision::NoModel;
    else if (r.ood >= cfg.ood_limit) r.decision = Decision::Uncertain;
    else r.decision = *r.p_final >= cfg.tau_call ? Decision::Pass : Decision::LowProb;

    // Genotype: HMM diploid model when reassessed, else the count model. A
    // called variant is never reported as 0/0.
    const GenotypeLik& g = r.hmm_genotype ? *r.hmm_genotype : r.count_genotype;
    r.gt = g.gt;
    if (r.gt == 0) r.gt = g.posterior[2] > g.posterior[1] ? 2 : 1;
    return r;
}

}  // namespace nindels
