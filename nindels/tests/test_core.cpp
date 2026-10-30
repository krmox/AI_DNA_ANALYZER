// Unit tests for the core formulas on synthetic inputs (no data needed).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "nindels/engine.hpp"
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
    std::printf("%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
