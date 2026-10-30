// Unit tests for the core formulas on synthetic inputs (no data needed).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "nindels/engine.hpp"
#include "nindels/locus.hpp"
#include "nindels/str_model.hpp"
#include "nindels/table_io.hpp"

using namespace nindels;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond)                                                          \
    do {                                                                     \
        if (cond) ++g_pass;                                                  \
        else { ++g_fail; std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); } \
    } while (0)
#define NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

static RefSeq make_ref(const std::string& s, int64_t off = 1000) { return RefSeq{s, off}; }

static void test_normalize() {
    //                 0123456789
    RefSeq ref = make_ref("GGCAAAAAATTCACACAGG");
    // deleting the last A of the run normalizes to the run start
    Indel d{0, 1000 + 8, IndelType::Deletion, "A"};
    Indel n = left_normalize(ref, d);
    CHECK(n.pos == 1003 && n.seq == "A");
    // inserting "CA" at the end of CACACA moves to its start
    Indel ins{0, 1000 + 17, IndelType::Insertion, "CA"};
    Indel m = left_normalize(ref, ins);
    CHECK(m.pos == 1011);
    CHECK(m.seq == "TC" || m.seq == "CA" || m.seq == "AC");
    // insertion of a base unrelated to its left neighbour stays put
    Indel z{0, 1000 + 2, IndelType::Insertion, "T"};
    CHECK(left_normalize(ref, z).pos == 1002);
}

static void test_pattern() {
    RefSeq ref = make_ref("GGCAAAAAATTCACACAGGTTTT");
    Pattern p = extract_pattern(ref, Indel{0, 1003, IndelType::Deletion, "A"});
    CHECK(p.hp_length == 6 && p.hp_base == 'A');
    CHECK(p.period == 1 && p.copies == 6 && p.repeat_span == 6);
    Pattern q = extract_pattern(ref, Indel{0, 1011, IndelType::Insertion, "CA"});
    CHECK(q.period == 2 && q.copies == 3);
    NEAR(shannon_entropy("AAAA"), 0.0, 1e-12);
    NEAR(shannon_entropy("ACGT"), 2.0, 1e-12);
}

static void test_similarity() {
    CHECK(levenshtein("kitten", "sitting") == 3);
    NEAR(d_sequence("", ""), 0.0, 0);
    NEAR(d_sequence("AC", ""), 1.0, 0);
    Pattern a, b;
    a.type = IndelType::Insertion;
    a.length = 1; a.motif = "A"; a.context = "ACGTACGT"; a.hp_length = 5;
    b = a;
    SimilarityParams sp;
    NEAR(pattern_distance(a, b, sp), 0.0, 1e-12);
    NEAR(pattern_similarity(a, b, sp), 1.0, 1e-12);
    b.type = IndelType::Deletion;
    FeatureDistances d = feature_distances(a, b, sp);
    NEAR(d.type, 1.0, 0);
    NEAR(pattern_distance(a, b, sp), 1.0 / 6.0, 1e-12);  // only d_type differs, 6 equal weights
    // d_length formula: min(1, |ln(1+1) - ln(1+3)| / s_L)
    b = a; b.length = 3;
    sp.s_L = 2.0;
    NEAR(feature_distances(a, b, sp).length, std::fabs(std::log(2.0) - std::log(4.0)) / 2.0, 1e-12);
    bool threw = false;
    try { SimilarityParams bad; bad.w_type = -1; bad.validate(); } catch (...) { threw = true; }
    CHECK(threw);
}

static void test_memory() {
    std::vector<LabeledPattern> items;
    Pattern hp;
    hp.type = IndelType::Deletion; hp.length = 1; hp.motif = "A"; hp.context = "CCCCCAAAAAAAAAGGGGG";
    hp.hp_length = 10; hp.hp_base = 'A'; hp.period = 1; hp.copies = 10; hp.repeat_span = 10;
    Pattern uniq;
    uniq.type = IndelType::Insertion; uniq.length = 4; uniq.motif = "ACGT"; uniq.context = "GATTACAGATTACAGATT";
    uniq.hp_length = 1; uniq.period = 4; uniq.copies = 1; uniq.repeat_span = 4;
    for (int i = 0; i < 40; ++i) items.push_back({hp, i % 4 == 0, 0, 1000 + i * 1000});    // 10 pos, 30 neg
    for (int i = 0; i < 20; ++i) items.push_back({uniq, i % 10 != 0, 0, 900000 + i * 1000}); // 18 pos, 2 neg
    PatternMemory m = PatternMemory::build(items, SimilarityParams{}, MemoryParams{}, "test");
    CHECK(m.classes().size() == 2);
    ContextStats c = m.query(hp);
    CHECK(c.retrieved >= 1);
    NEAR(*c.s_max, 1.0, 1e-12);
    // P_context = (N_pos + alpha) / (N_pos + N_neg + alpha + beta)
    NEAR(c.p_context, (c.n_positive + m.alpha()) / (c.n_positive + c.n_negative + m.alpha() + m.beta()), 1e-12);
    CHECK(c.p_context < 0.5);
    CHECK(m.query(uniq).p_context > 0.5);
    // N_eff (Kish) with a single class equals its observation count, and the
    // distinct-locus version equals its number of distinct loci.
    PatternMemory one = PatternMemory::build({items.begin(), items.begin() + 40}, SimilarityParams{}, MemoryParams{}, "");
    ContextStats c1 = one.query(hp);
    NEAR(*c1.n_eff_kish, 40.0, 1e-9);
    NEAR(*c1.n_eff_loci, 40.0, 1e-9);
    // Same locus observed 40 times: Kish still says 40, distinct-locus says 1.
    std::vector<LabeledPattern> dup(40, LabeledPattern{hp, true, 0, 5000});
    ContextStats cd = PatternMemory::build(dup, SimilarityParams{}, MemoryParams{}, "").query(hp);
    NEAR(*cd.n_eff_kish, 40.0, 1e-9);
    NEAR(*cd.n_eff_loci, 1.0, 1e-9);
    // round trip
    m.save("/tmp/nindels_test_memory.tsv");
    PatternMemory r = PatternMemory::load("/tmp/nindels_test_memory.tsv");
    NEAR(r.query(hp).p_context, c.p_context, 1e-12);
}

static void test_evidence() {
    std::vector<ReadObservation> obs;
    for (int i = 0; i < 30; ++i) obs.push_back({ReadAllele::Alt, i % 2 == 0, 60, 0.3f, 0, 35, 0.01f, false});
    for (int i = 0; i < 70; ++i) obs.push_back({ReadAllele::Ref, i % 2 == 0, 60, 0.3f, 0, 35, 0.0f, false});
    Evidence e = summarize_evidence(obs);
    NEAR(*e.f, 0.30, 1e-12);  // spec §3 worked example
    NEAR(e.l_depth, std::log(101.0), 1e-12);
    NEAR(*e.alt_fwd_frac, 0.5, 1e-12);
    CHECK(*e.strand_bias < 0.05);
    Evidence empty = summarize_evidence({});
    CHECK(!empty.f.has_value());  // denominator 0 -> unavailable
}

static std::vector<uint8_t> quals(size_t n, uint8_t q = 30) { return std::vector<uint8_t>(n, q); }

static void test_pairhmm() {
    std::string hap = "TTGACCTAGGCATTAGCAGGCTACGATCGGAT";
    std::string read = hap.substr(5, 20);
    auto q = quals(read.size());
    double ll = pairhmm_loglik(read, q.data(), hap);
    std::string other = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    CHECK(ll > pairhmm_loglik(read, q.data(), other) + 50);
    CHECK(ll < 0);
}

// Builds reads from a haplotype sequence for the likelihood tests.
static HmmRead read_from(const std::string& hapseq, int64_t ref_start, int64_t ref_end, size_t from, size_t len,
                         ReadAllele a, uint64_t h) {
    HmmRead r;
    r.bases = hapseq.substr(from, len);
    r.quals = quals(len);
    r.ref_start = ref_start;
    r.ref_end = ref_end;
    r.mapq = 60;
    r.cigar_allele = a;
    r.name_hash = h;
    r.window_bq.assign(21, 30);
    return r;
}

static void test_genotypes() {
    std::string left = "GATCCTAGGCATTGACCGTATGCAGGTCA", right = "TGCATTGCAGCTTACGATCAGGCTAAC";
    std::string refs = left + "CTG" + right;   // deletion of CTG
    RefSeq ref{refs, 0};
    Indel del{0, static_cast<int64_t>(left.size()), IndelType::Deletion, "CTG"};
    Haplotypes H = build_haplotypes(ref, del, 0, static_cast<int64_t>(refs.size()));
    CHECK(H.alt == left + right);

    std::vector<HmmRead> het, hom;
    for (int i = 0; i < 10; ++i) {
        het.push_back(read_from(H.ref, 5, 55, 5, 50, ReadAllele::Ref, 2 * i));
        het.push_back(read_from(H.alt, 5, 55, 5, 47, ReadAllele::Alt, 2 * i + 1));
        hom.push_back(read_from(H.alt, 5, 55, 5, 47, ReadAllele::Alt, i));
    }
    LikelihoodParams lp;
    GenotypePriors pr;
    auto lh = read_likelihoods(ref, del, het, lp);
    GenotypeLik g = diploid_gl(lh, true, pr);
    CHECK(g.gt == 1);                       // heterozygous reads -> 0/1
    CHECK(g.gl[1] > g.gl[0] && g.gl[1] > g.gl[2]);
    // The original Delta_ell is a hom-alt vs hom-ref test. For a balanced het
    // the ref reads penalise H_alt as much as the alt reads penalise H_ref, so it
    // is near zero even though the variant is real (the reason for the correction).
    double dl_het = delta_ell_original(lh);
    CHECK(std::fabs(dl_het) < 0.25 * std::fabs(g.gl[1] - g.gl[0]) + 5.0);
    auto lhh = read_likelihoods(ref, del, hom, lp);
    CHECK(diploid_gl(lhh, true, pr).gt == 2);
    CHECK(delta_ell_original(lhh) > 20);
    NEAR(*hmm_consistency(lh), 1.0, 1e-12);
    // Banded DP must agree with the full DP on these reads.
    LikelihoodParams full = lp;
    full.band = -1;
    auto lf = read_likelihoods(ref, del, het, full);
    for (size_t i = 0; i < lf.size(); ++i) {
        NEAR(lh[i].lp_ref, lf[i].lp_ref, 1e-6);
        NEAR(lh[i].lp_alt, lf[i].lp_alt, 1e-6);
    }

    GenotypeLik cg = count_gl(15, 15, 0.01, pr);
    CHECK(cg.gt == 1);
    CHECK(count_gl(0, 30, 0.01, pr).gt == 0);
    CHECK(count_gl(30, 0, 0.01, pr).gt == 2);

    // MAPQ 0 reads cannot move the genotype.
    for (auto& l : lhh) l.mapq = 0;
    GenotypeLik z = diploid_gl(lhh, true, pr);
    NEAR(z.gl[0], z.gl[2], 1e-9);
}

static void test_reassess_and_complexity() {
    ReassessParams rp;
    rp.A = 0; NEAR(tau_review(rp), rp.tau_strict, 1e-12);
    rp.A = 1; NEAR(tau_review(rp), rp.tau_loose, 1e-12);
    NEAR(uncertainty(0.5), 1.0, 1e-12);
    NEAR(uncertainty(0.0), 0.0, 1e-12);
    NEAR(uncertainty(1.0), 0.0, 1e-12);
    NEAR(reassess_priority(1, 0, 0, rp), 1.0 / 3, 1e-12);
    bool threw = false;
    try { ReassessParams b; b.a_U = 0.9; b.validate(); } catch (...) { threw = true; }
    CHECK(threw);

    Pattern p;
    p.hp_length = 15; p.repeat_span = 50; p.entropy = 0;
    ComplexityParams cp;
    Complexity c = region_complexity(p, 1.0, cp);
    NEAR(c.C, 1.0, 1e-12);
    Complexity c2 = region_complexity(p, std::nullopt, cp);  // G missing: renormalised, not zero
    NEAR(c2.C, 1.0, 1e-12);
}

static void test_logistic() {
    // y depends on x1 only; x2 is noise and sometimes missing.
    std::vector<std::vector<Opt>> X;
    std::vector<int> y;
    uint64_t s = 7;
    for (int i = 0; i < 4000; ++i) {
        s = mix64(s);
        double x1 = static_cast<double>(s % 10000) / 10000.0 * 6 - 3;
        s = mix64(s);
        double u = static_cast<double>(s % 10000) / 10000.0;
        double p = 1 / (1 + std::exp(-(0.5 + 2 * x1)));
        s = mix64(s);
        Opt x2 = (s % 3 == 0) ? Opt{} : Opt{static_cast<double>(s % 100) / 100.0};
        X.push_back({x1, x2});
        y.push_back(u < p);
    }
    LogisticModel m = LogisticModel::fit({"x1", "x2"}, X, y, 1e-6);
    // Standardized coefficient = 2 * sd(x1) ~ 2 * sqrt(3) for uniform(-3, 3)
    NEAR(m.coefficients()[0], 2 * std::sqrt(3.0), 0.35);
    CHECK(m.predict(std::vector<Opt>{2.0, std::nullopt}) > 0.95);
    CHECK(m.predict(std::vector<Opt>{-2.0, 0.5}) < 0.05);
    m.save("/tmp/nindels_test_model.tsv");
    LogisticModel r = LogisticModel::load("/tmp/nindels_test_model.tsv");
    NEAR(r.predict(std::vector<Opt>{0.3, std::nullopt}), m.predict(std::vector<Opt>{0.3, std::nullopt}), 1e-12);
}

static void test_engine() {
    std::string left = "GATCCTAGGCATTGACCGTATGCAGGTCA", right = "TGCATTGCAGCTTACGATCAGGCTAAC";
    std::string refs = left + "CTG" + right;
    RefSeq ref{refs, 0};
    Indel del{0, static_cast<int64_t>(left.size()), IndelType::Deletion, "CTG"};
    Haplotypes H = build_haplotypes(ref, del, 0, static_cast<int64_t>(refs.size()));
    CandidateInput in;
    in.indel = del;
    for (int i = 0; i < 12; ++i) {
        bool alt = i % 2;
        in.obs.push_back({alt ? ReadAllele::Alt : ReadAllele::Ref, i % 4 < 2, 60, 0.3f, 0, 35, 0.0f, false});
        in.reads.push_back(alt ? read_from(H.alt, 5, 55, 5, 47, ReadAllele::Alt, i)
                               : read_from(H.ref, 5, 55, 5, 50, ReadAllele::Ref, i));
    }
    EngineConfig cfg;
    cfg.force_reassess = true;
    CandidateResult r = process_indel(in, ref, nullptr, nullptr, nullptr, cfg);
    CHECK(r.reassessed && r.hmm_genotype.has_value());
    CHECK(r.decision == Decision::Pass);  // untrained genotype rule
    CHECK(r.gt == 1);
    CHECK(r.features.size() == all_feature_names().size());
    for (size_t i = 0; i < r.features.size(); ++i) CHECK(r.features[i].name == all_feature_names()[i]);
    // A changes who is reassessed, never the features of a reassessed candidate.
    EngineConfig cfg0 = cfg;
    cfg0.force_reassess = false;
    cfg0.reassess.A = 0;
    cfg0.reassess.tau_strict = cfg0.reassess.tau_loose = 2.0;  // nobody reassessed
    CandidateResult r0 = process_indel(in, ref, nullptr, nullptr, nullptr, cfg0);
    CHECK(!r0.reassessed && !r0.hmm_genotype.has_value());
    CHECK(*r0.features[0].value == *r.features[0].value);
}

static void test_str_model() {
    //                      0123456789012345678901234
    auto runs = find_str_runs("GCTAAAAAAGTCACACACATGG");
    CHECK(runs.size() == 2);
    CHECK(runs[0].period == 1 && runs[0].start == 3 && runs[0].end == 9);
    CHECK(runs[1].period == 2 && runs[1].start == 11 && runs[1].units() == 4);  // CACACACA
    CHECK(find_str_runs("GCTAAAGC").empty());  // homopolymer of 3: below minimum

    // EM recovers eps from a mixture of hom-ref (eps = 0.05), het and hom-alt loci.
    std::vector<StrLocusCount> c;
    uint64_t s = 3;
    for (int k = 0; k < 2000; ++k) {
        s = mix64(s);
        int z = (s % 100) < 85 ? 0 : ((s % 100) < 95 ? 1 : 2);
        double pr = z == 0 ? 0.05 : (z == 1 ? 0.5 : 0.95);
        int n = 30, x = 0;
        for (int r = 0; r < n; ++r) { s = mix64(s); x += (static_cast<double>(s % 1000000) / 1e6) < pr; }
        c.push_back({1, 10, n, x});
    }
    StrErrorModel m;
    m.fit(c);
    NEAR(m.eps(1, 10), 0.05, 0.005);
    CHECK(m.eps(1, 11) == m.eps(1, 10));  // sparse stratum inherits from u-1
    m.calibrate(PairHmmParams{});
    CHECK(m.gop(1, 10) > 1e-3);
    // calibrated ratio: P(1-unit deletion read) / P(clean read) ~= eps / (2 (1 - eps))
    std::string hap = "GATCCTAGGCATTGACCGTA" + std::string(10, 'A') + "TGCAGGTCATGCATTGCAGC";
    std::string slip = "GATCCTAGGCATTGACCGTA" + std::string(9, 'A') + "TGCAGGTCATGCATTGCAGC";
    auto prof = m.gap_open_profile(hap, std::pow(10.0, -4.5));
    std::vector<uint8_t> q(hap.size(), 30);
    double R = std::exp(pairhmm_loglik(slip, q.data(), hap, {}, 0, -1, prof.data()) -
                        pairhmm_loglik(hap, q.data(), hap, {}, 0, -1, prof.data()));
    NEAR(R, 0.05 / (2 * 0.95), 0.05 / (2 * 0.95) * 0.05);
    m.save("/tmp/nindels_test_str.tsv");
    NEAR(StrErrorModel::load("/tmp/nindels_test_str.tsv").gop(1, 10), m.gop(1, 10), 1e-9);
}

static void test_locus() {
    std::string left = "GATCCTAGGCATTGACCGTATGCAGGTCAG", right = "CTGCATTGCAGCTTACGATCAGGCTAACG";
    std::string refs = left + std::string(10, 'A') + right;
    RefSeq ref{refs, 0};
    const int64_t hp = static_cast<int64_t>(left.size());
    Indel d1{0, hp, IndelType::Deletion, "A"}, d2{0, hp, IndelType::Deletion, "AA"};
    std::string h1 = left + std::string(9, 'A') + right, h2 = left + std::string(8, 'A') + right;
    auto mk = [&](const std::string& hs, ReadAllele a, uint64_t id) {
        return read_from(hs, 2, static_cast<int64_t>(refs.size()) - 2, 2, hs.size() - 4, a, id);
    };
    // 1/2: half the reads -A, half -AA, no reference reads.
    CandidateInput c1, c2;
    c1.indel = d1;
    c2.indel = d2;
    for (int i = 0; i < 10; ++i) {
        c1.reads.push_back(mk(h1, ReadAllele::Alt, 100 + i));
        c1.obs.push_back({ReadAllele::Alt, false, 60, 0.3f, 0, 30, 0, false});
        c2.reads.push_back(mk(h2, ReadAllele::Alt, 200 + i));
        c2.obs.push_back({ReadAllele::Alt, false, 60, 0.3f, 0, 30, 0, false});
    }
    LocusCall c = genotype_locus(ref, {&c1, &c2}, LikelihoodParams{}, LocusParams{});
    CHECK(c.pass);
    CHECK(c.g1 != 0 && c.g2 != 0 && c.g1 != c.g2);  // 1/2, not two homozygous calls
    CHECK(c.ad[0] == 0 && c.ad[1] == 10 && c.ad[2] == 10);
    // the biallelic model (baseline) calls -A as 1/1 on the same reads
    std::vector<HmmRead> all = c1.reads;
    for (auto r : c2.reads) { r.cigar_allele = ReadAllele::Other; all.push_back(r); }
    auto lk = read_likelihoods(ref, d1, all, LikelihoodParams{});
    CHECK(diploid_gl(lk, true, GenotypePriors{}).gt == 2);

    // 0/1 with reference reads
    CandidateInput h;
    h.indel = d1;
    for (int i = 0; i < 10; ++i) {
        h.reads.push_back(mk(refs, ReadAllele::Ref, 300 + i));
        h.reads.push_back(mk(h1, ReadAllele::Alt, 400 + i));
    }
    LocusCall ch = genotype_locus(ref, {&h}, LikelihoodParams{}, LocusParams{});
    CHECK(ch.pass && ch.g1 == 0 && ch.g2 == 1);
    // 6 slipped reads out of 30 in a 10-bp homopolymer. Constant Q45: each
    // slipped read is worth ~+7.4 nats for 0/1, which beats the prior; with the
    // estimated eps ~ 2% it is worth ~+3.8 and the locus stays hom-ref.
    CandidateInput n;
    n.indel = d1;
    for (int i = 0; i < 24; ++i) n.reads.push_back(mk(refs, ReadAllele::Ref, 500 + i));
    for (int i = 0; i < 6; ++i) n.reads.push_back(mk(h1, ReadAllele::Alt, 600 + i));
    std::vector<StrLocusCount> cnt;
    for (int k = 0; k < 100; ++k) cnt.push_back({1, 10, 50, k % 10 == 0 ? 2 : 1});
    StrErrorModel sm;
    sm.fit(cnt);
    sm.calibrate(PairHmmParams{});
    LikelihoodParams with_str;
    with_str.str = &sm;
    LocusCall cn0 = genotype_locus(ref, {&n}, LikelihoodParams{}, LocusParams{});
    LocusCall cn1 = genotype_locus(ref, {&n}, with_str, LocusParams{});
    CHECK(cn0.pass);   // constant Q45 gap-open: called 0/1 (the baseline failure)
    CHECK(!cn1.pass);  // data-estimated slippage: not called
    CHECK(cn1.qual < cn0.qual);
    std::printf("  locus test: 6/30 slipped reads QUAL constant-Q45 %.1f vs STR model %.1f (eps %.4f)\n", cn0.qual,
                cn1.qual, sm.eps(1, 10));
}

static void test_table_io() {
    std::string line = "1\tNA\t-2.5";
    auto v = split_tab(line);
    CHECK(v.size() == 3);
    CHECK(!parse_opt("NA").has_value());
    NEAR(*parse_opt("-2.5"), -2.5, 0);
    CHECK(format_opt(std::nullopt) == "NA");
}

// Synthetic labelled feature table for the CLI smoke test: true heterozygous
// deletions with balanced, clean evidence vs. errors in a long homopolymer
// with low, strand-biased support.
static void write_synthetic_table(const std::string& path) {
    std::string left = "GATCCTAGGCATTGACCGTATGCAGGTCA", right = "TGCATTGCAGCTTACGATCAGGCTAAC";
    std::ofstream o(path);
    std::vector<std::pair<CandidateResult, int>> rows;
    uint64_t s = 11;
    for (int i = 0; i < 400; ++i) {
        const bool truth = i % 2 == 0;
        std::string mid = truth ? "CTG" : "AAAAAAAAAAAA";
        std::string refs = left + mid + right;
        RefSeq ref{refs, 100000LL * i};
        Indel x{0, ref.offset + static_cast<int64_t>(left.size()), IndelType::Deletion, truth ? "CTG" : "A"};
        Haplotypes H = build_haplotypes(ref, x, ref.offset, ref.end());
        CandidateInput in;
        in.indel = x;
        s = mix64(s);
        const int n_alt = truth ? 10 + static_cast<int>(s % 6) : 2 + static_cast<int>(s % 3);
        const int n_ref = truth ? 10 + static_cast<int>((s >> 8) % 6) : 25;
        for (int k = 0; k < n_alt + n_ref; ++k) {
            bool alt = k < n_alt;
            bool rev = truth ? k % 2 : (alt ? true : k % 2);
            in.obs.push_back({alt ? ReadAllele::Alt : ReadAllele::Ref, rev, 60, alt && !truth ? 0.05f : 0.3f, 0,
                              truth ? 35.0f : 22.0f, truth ? 0.0f : 0.03f, false});
            const std::string& hs = alt ? H.alt : H.ref;
            in.reads.push_back(read_from(hs, ref.offset + 3, ref.end() - 3, 3, hs.size() - 6,
                                         alt ? ReadAllele::Alt : ReadAllele::Ref, static_cast<uint64_t>(k)));
        }
        EngineConfig cfg;
        cfg.force_reassess = true;
        rows.emplace_back(process_indel(in, ref, nullptr, nullptr, nullptr, cfg), truth ? 1 : 0);
        if (i == 0) write_table_header(o, {"truth_indels_in_bed=" + std::to_string(200)});
        write_table_row(o, "chrS", ref, rows.back().first, rows.back().second);
    }
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--write-table") {
        write_synthetic_table(argv[2]);
        return 0;
    }
    test_normalize();
    test_pattern();
    test_similarity();
    test_memory();
    test_evidence();
    test_pairhmm();
    test_genotypes();
    test_reassess_and_complexity();
    test_logistic();
    test_engine();
    test_table_io();
    test_str_model();
    test_locus();
    std::printf("%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
