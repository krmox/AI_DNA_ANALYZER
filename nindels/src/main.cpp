// nindels-call: command line entry point.
//
//   call          BAM -> candidates -> Pattern Engine -> VCF and/or feature table
//   build-memory  labelled feature tables -> pattern memory (spec §6)
//   fit           labelled feature table -> logistic model (spec §12)
//   select-tau    scored validation table -> tau_call maximising F1
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "nindels/engine.hpp"
#include "nindels/table_io.hpp"
#include "nindels/vcf_out.hpp"

#ifdef NINDELS_HAVE_HTSLIB
#include "nindels/bam_scan.hpp"
#include "nindels/truth.hpp"
#endif

using namespace nindels;

namespace {

struct Args {
    std::map<std::string, std::vector<std::string>> kv;
    std::set<std::string> flags;

    static Args parse(int argc, char** argv, int from, const std::set<std::string>& flag_names) {
        Args a;
        for (int i = from; i < argc; ++i) {
            std::string k = argv[i];
            if (k.rfind("--", 0) != 0) throw std::invalid_argument("unexpected argument: " + k);
            k = k.substr(2);
            if (flag_names.count(k)) { a.flags.insert(k); continue; }
            if (i + 1 >= argc) throw std::invalid_argument("--" + k + " needs a value");
            a.kv[k].push_back(argv[++i]);
        }
        return a;
    }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
    std::string get(const std::string& k, const std::string& def = "") const {
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second.back();
    }
    std::string need(const std::string& k) const {
        if (!has(k)) throw std::invalid_argument("missing --" + k);
        return get(k);
    }
    double num(const std::string& k, double def) const { return has(k) ? std::stod(get(k)) : def; }
    const std::vector<std::string>& all(const std::string& k) const {
        static const std::vector<std::string> none;
        auto it = kv.find(k);
        return it == kv.end() ? none : it->second;
    }
};

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string x;
    while (std::getline(ss, x, ','))
        if (!x.empty()) out.push_back(x);
    return out;
}

// Parity of a 2-Mb block (default). Used for two-fold pattern memory on the
// training chromosome, so a training row never sees its own label in memory.
std::string parity_of(int64_t pos, int64_t block) { return (pos / block) % 2 == 0 ? "even" : "odd"; }

// ------------------------------------------------------------------ presets

// Feature subsets = the benchmarked variants. "spec-x" is the spec §12 X_i verbatim.
const std::vector<std::string> kEvidence = {"f", "l_depth", "alt_fwd_frac", "strand_bias", "readpos_alt",
                                            "readpos_diff", "clip_alt", "clip_all", "mismatch_alt",
                                            "other_frac", "split_frac", "mapq_alt", "mq0_frac", "mqlow_frac",
                                            "bq_alt", "count_qual"};

std::vector<std::string> preset(const std::string& name) {
    auto cat = [](std::initializer_list<std::vector<std::string>> parts) {
        std::vector<std::string> v;
        for (const auto& p : parts) v.insert(v.end(), p.begin(), p.end());
        return v;
    };
    const std::vector<std::string> mem_k = {"p_context", "n_eff_kish", "s_max"};
    const std::vector<std::string> mem_l = {"p_context", "n_eff_loci", "s_max"};
    const std::vector<std::string> dip = {"hmm_qual", "hmm_p_hom", "hmm_consistency", "q_min"};
    const std::vector<std::string> orig = {"delta_ell_orig", "hmm_consistency", "q_min"};
    if (name == "full") return cat({kEvidence, mem_k, {"c_r"}, dip});
    if (name == "full-neff-loci") return cat({kEvidence, mem_l, {"c_r"}, dip});
    if (name == "no-memory") return cat({kEvidence, {"c_r"}, dip});
    if (name == "dl-orig") return cat({kEvidence, mem_k, {"c_r"}, orig});
    if (name == "stage1") return cat({kEvidence, mem_k, {"c_r"}});
    if (name == "stage1-no-memory") return cat({kEvidence, {"c_r"}});
    if (name == "evidence-only") return kEvidence;
    if (name == "spec-x") return {"f", "l_depth", "p_context", "n_eff_kish", "s_max", "delta_ell_orig", "c_r"};
    throw std::invalid_argument("unknown preset " + name);
}

// ---------------------------------------------------------------- build-memory

int cmd_build_memory(const Args& a) {
    MemoryParams mp;
    mp.theta = a.num("theta", mp.theta);
    mp.top_k = static_cast<int>(a.num("top-k", mp.top_k));
    const std::string parity = a.get("parity", "all");
    const int64_t block = static_cast<int64_t>(a.num("block", 2e6));
    std::vector<LabeledPattern> items;
    std::vector<Pattern> pats;
    int skipped = 0;
    for (const std::string& path : a.all("table")) {
        Table t = read_table(path);
        const int cl = t.col("label"), cp = t.col("pos0"), cc = t.col("chrom");
        std::map<std::string, int32_t> chrom_ids;
        for (const auto& row : t.rows) {
            int label = std::stoi(row[cl]);
            int64_t pos = std::stoll(row[cp]);
            if (label < 0) { ++skipped; continue; }
            if (parity != "all" && parity_of(pos, block) != parity) continue;
            auto [it, _] = chrom_ids.emplace(row[cc], static_cast<int32_t>(chrom_ids.size()));
            items.push_back({pattern_from_row(t, row), label == 1, it->second, pos});
            pats.push_back(items.back().pattern);
        }
    }
    if (items.empty()) throw std::runtime_error("no labelled rows");
    SimilarityParams sp = estimate_scales(pats, SimilarityParams{});
    PatternMemory m = PatternMemory::build(std::move(items), sp, mp, a.get("metadata", "."));
    m.save(a.need("out"));
    size_t npos = 0, nneg = 0;
    for (const auto& c : m.classes()) { npos += c.n_pos; nneg += c.n_neg; }
    std::fprintf(stderr, "memory: %zu classes from %zu pos / %zu neg (%d unlabelled skipped); alpha=%.3f beta=%.3f "
                 "s_L=%.3f s_h=%.3f s_r=%.3f\n", m.classes().size(), npos, nneg, skipped, m.alpha(), m.beta(),
                 sp.s_L, sp.s_h, sp.s_r);
    return 0;
}

// ------------------------------------------------------------------------ fit

int cmd_fit(const Args& a) {
    std::vector<std::string> names = a.has("features") ? split_csv(a.get("features")) : preset(a.get("preset", "full"));
    const double mask_frac = a.num("mask-reassess", 0.5);
    const uint64_t seed = static_cast<uint64_t>(a.num("seed", 1));
    std::vector<std::vector<Opt>> rows;
    std::vector<int> labels;
    for (const std::string& path : a.all("table")) {
        Table t = read_table(path);
        std::vector<int> idx;
        for (const auto& n : names) idx.push_back(t.col(n));
        const int cl = t.col("label"), cp = t.col("pos0");
        for (const auto& row : t.rows) {
            int label = std::stoi(row[cl]);
            if (label < 0) continue;
            std::vector<Opt> x;
            for (int i : idx) x.push_back(parse_opt(row[static_cast<size_t>(i)]));
            // Missingness augmentation: a deterministic fraction of rows lose
            // their step-13 features, so the model learns the not-reassessed case.
            uint64_t h = mix64(seed ^ static_cast<uint64_t>(std::stoll(row[cp])));
            if (static_cast<double>(h % 1000000) / 1e6 < mask_frac)
                for (size_t j = 0; j < names.size(); ++j)
                    for (const auto& rn : reassessment_feature_names())
                        if (names[j] == rn) x[j] = std::nullopt;
            rows.push_back(std::move(x));
            labels.push_back(label);
        }
    }
    LogisticModel m = LogisticModel::fit(names, rows, labels, a.num("lambda", 1e-2));
    m.save(a.need("out"));
    int npos = static_cast<int>(std::count(labels.begin(), labels.end(), 1));
    std::fprintf(stderr, "fit: %zu rows (%d pos), %zu features, intercept %.4f\n", rows.size(), npos, names.size(),
                 m.intercept());
    for (size_t j = 0; j < names.size(); ++j)
        std::fprintf(stderr, "  %-16s beta=%+.4f\n", names[j].c_str(), m.coefficients()[j]);
    return 0;
}

// ----------------------------------------------------------------- select-tau

// F1 over labelled candidates. Recall denominator = truth indels in the BED
// (from the table header), so truth indels that were never candidates count
// as false negatives. Candidate-level matching is exact-allele; final numbers
// must come from vcfeval.
int cmd_select_tau(const Args& a) {
    Table t = read_table(a.need("table"));
    int n_truth = -1;
    for (const auto& c : t.comments)
        if (c.rfind("truth_indels_in_bed=", 0) == 0) n_truth = std::stoi(c.substr(20));
    if (n_truth <= 0) throw std::runtime_error("table lacks ##truth_indels_in_bed");
    const int cl = t.col("label"), cp = t.col("p_final"), cd = t.col("decision");
    std::vector<std::pair<double, int>> scored;
    for (const auto& row : t.rows) {
        int label = std::stoi(row[cl]);
        Opt p = parse_opt(row[cp]);
        if (label < 0 || !p || row[cd] == "UNCERTAIN") continue;
        scored.emplace_back(*p, label);
    }
    if (scored.empty()) throw std::runtime_error("no labelled rows with p_final (run call with --model)");
    double best_f1 = -1, best_tau = 0.5;
    std::printf("tau\tTP\tFP\tFN\tprecision\trecall\tF1\n");
    for (int k = 1; k < 100; ++k) {
        double tau = k / 100.0;
        int tp = 0, fp = 0;
        for (auto [p, l] : scored)
            if (p >= tau) (l ? tp : fp) += 1;
        int fn = n_truth - tp;
        double prec = tp + fp ? static_cast<double>(tp) / (tp + fp) : 0, rec = static_cast<double>(tp) / n_truth;
        double f1 = prec + rec > 0 ? 2 * prec * rec / (prec + rec) : 0;
        std::printf("%.2f\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\n", tau, tp, fp, fn, prec, rec, f1);
        if (f1 > best_f1) { best_f1 = f1; best_tau = tau; }
    }
    std::fprintf(stderr, "best tau_call=%.2f F1=%.5f (candidate-level, exact allele)\n", best_tau, best_f1);
    return 0;
}

#ifdef NINDELS_HAVE_HTSLIB

struct Region {
    std::string chrom;
    int64_t start = 0, end = -1;  // 0-based half-open; end -1 = contig end
};

Region parse_region(const std::string& s) {
    Region r;
    size_t c = s.find(':');
    r.chrom = s.substr(0, c);
    if (c != std::string::npos) {
        std::string rest = s.substr(c + 1);
        rest.erase(std::remove(rest.begin(), rest.end(), ','), rest.end());
        size_t d = rest.find('-');
        r.start = std::stoll(rest.substr(0, d)) - 1;
        if (d != std::string::npos) r.end = std::stoll(rest.substr(d + 1));
    }
    return r;
}

int cmd_call(const Args& a, const std::string& cmdline) {
    const std::string bam = a.need("bam"), fasta = a.need("ref");
    Region reg = parse_region(a.need("region"));
    const int threads = std::max(1, static_cast<int>(a.num("threads", 1)));

    EngineConfig cfg;
    cfg.reassess.A = a.num("A", cfg.reassess.A);
    cfg.reassess.tau_strict = a.num("tau-strict", cfg.reassess.tau_strict);
    cfg.reassess.tau_loose = a.num("tau-loose", cfg.reassess.tau_loose);
    cfg.reassess.validate();
    cfg.complexity.validate();
    cfg.tau_call = a.num("tau-call", cfg.tau_call);
    cfg.rule_min_qual = a.num("min-qual", cfg.rule_min_qual);
    cfg.force_reassess = a.flags.count("force-reassess") > 0;
    cfg.likelihood.max_reads = static_cast<int>(a.num("max-hmm-reads", cfg.likelihood.max_reads));
    ScanParams sp;
    sp.min_alt = static_cast<int>(a.num("min-alt", sp.min_alt));
    sp.min_alt_frac = a.num("min-alt-frac", sp.min_alt_frac);
    sp.chunk = static_cast<int64_t>(a.num("chunk", static_cast<double>(sp.chunk)));
    const std::string parity = a.get("parity", "all");
    const int64_t block = static_cast<int64_t>(a.num("block", 2e6));

    std::unique_ptr<PatternMemory> memory;
    std::unique_ptr<LogisticModel> base_model, final_model;
    if (a.has("memory")) memory = std::make_unique<PatternMemory>(PatternMemory::load(a.get("memory")));
    if (a.has("base-model")) base_model = std::make_unique<LogisticModel>(LogisticModel::load(a.get("base-model")));
    if (a.has("model")) final_model = std::make_unique<LogisticModel>(LogisticModel::load(a.get("model")));

    BamScanner probe(bam, fasta);
    const int32_t tid = probe.tid(reg.chrom);
    if (tid < 0) throw std::runtime_error("contig not in BAM: " + reg.chrom);
    int64_t contig_len = 0;
    for (auto& [n, l] : probe.contigs())
        if (n == reg.chrom) contig_len = l;
    if (reg.end < 0 || reg.end > contig_len) reg.end = contig_len;

    const int64_t ref_pad = 20000;
    RefSeq ref = load_reference(fasta, reg.chrom, std::max<int64_t>(0, reg.start - ref_pad),
                                std::min(contig_len, reg.end + ref_pad));

    std::unique_ptr<BedRegions> bed;
    std::unique_ptr<TruthSet> truth;
    if (a.has("bed")) bed = std::make_unique<BedRegions>(load_bed(a.get("bed"), reg.chrom, reg.start, reg.end));
    if (a.has("truth"))
        truth = std::make_unique<TruthSet>(load_truth(a.get("truth"), reg.chrom, reg.start, reg.end, tid, ref, bed.get()));
    std::vector<Indel> external;
    if (a.has("candidates-vcf")) external = load_vcf_indels(a.get("candidates-vcf"), reg.chrom, reg.start, reg.end, tid, ref);

    std::ofstream vcf, table;
    if (a.has("out-vcf")) {
        vcf.open(a.get("out-vcf"));
        VcfHeaderInfo h;
        h.sample = a.get("sample", probe.sample_name());
        h.contigs = probe.contigs();
        h.extra.push_back("nindels_command=" + cmdline);
        write_vcf_header(vcf, h);
    }
    if (a.has("out-table")) {
        table.open(a.get("out-table"));
        std::vector<std::string> comments = {"nindels_command=" + cmdline};
        if (truth) {
            comments.push_back("truth_indels_in_bed=" + std::to_string(truth->n_in_bed));
            comments.push_back("truth_complex_alleles_skipped=" + std::to_string(truth->n_complex));
        }
        write_table_header(table, comments);
    }

    std::vector<std::pair<int64_t, int64_t>> chunks;
    for (int64_t s = reg.start; s < reg.end; s += sp.chunk) chunks.emplace_back(s, std::min(reg.end, s + sp.chunk));

    size_t n_cand = 0, n_pass = 0, n_reassessed = 0;
    std::atomic<long long> t_scan_us{0}, t_engine_us{0};
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto us = [](auto d) { return std::chrono::duration_cast<std::chrono::microseconds>(d).count(); };
    std::vector<std::unique_ptr<BamScanner>> scanners;
    for (int t = 0; t < threads; ++t) scanners.push_back(std::make_unique<BamScanner>(bam, fasta));

    // Batches of `threads` chunks: scan them in parallel, then score every
    // candidate of the batch with a shared work index (per-candidate load
    // balancing), then write in genomic order.
    for (size_t c0 = 0; c0 < chunks.size(); c0 += static_cast<size_t>(threads)) {
        const size_t c1 = std::min(chunks.size(), c0 + static_cast<size_t>(threads));
        std::vector<std::vector<CandidateInput>> scanned(c1 - c0);
        std::vector<std::string> errors(static_cast<size_t>(threads));
        {
            std::vector<std::thread> pool;
            for (size_t c = c0; c < c1; ++c)
                pool.emplace_back([&, c] {
                    try {
                        auto t0 = now();
                        scanned[c - c0] = scanners[c - c0]->scan(ref, tid, chunks[c].first, chunks[c].second, sp,
                                                                 a.has("candidates-vcf") ? &external : nullptr,
                                                                 cfg.nucleotide.half_window);
                        t_scan_us += us(now() - t0);
                    } catch (const std::exception& e) {
                        errors[c - c0] = e.what();
                    }
                });
            for (auto& th : pool) th.join();
        }
        for (const auto& e : errors)
            if (!e.empty()) throw std::runtime_error(e);
        std::vector<const CandidateInput*> todo;
        for (auto& v : scanned)
            for (auto& in : v)
                if (parity == "all" || parity_of(in.indel.pos, block) == parity) todo.push_back(&in);
        std::vector<CandidateResult> results(todo.size());
        std::atomic<size_t> next{0};
        {
            std::vector<std::thread> pool;
            for (int t = 0; t < threads; ++t)
                pool.emplace_back([&, t] {
                    try {
                        auto t0 = now();
                        for (size_t k; (k = next.fetch_add(1)) < todo.size();)
                            results[k] = process_indel(*todo[k], ref, memory.get(), base_model.get(),
                                                       final_model.get(), cfg);
                        t_engine_us += us(now() - t0);
                    } catch (const std::exception& e) {
                        errors[static_cast<size_t>(t)] = e.what();
                    }
                });
            for (auto& th : pool) th.join();
        }
        for (const auto& e : errors)
            if (!e.empty()) throw std::runtime_error(e);
        for (const auto& r : results) {
            ++n_cand;
            n_pass += r.decision == Decision::Pass;
            n_reassessed += r.reassessed;
            if (vcf.is_open()) write_vcf_record(vcf, reg.chrom, ref, r);
            if (table.is_open())
                write_table_row(table, reg.chrom, ref, r, truth ? truth_label(*truth, bed.get(), r.indel) : -1);
        }
    }
    std::fprintf(stderr, "call: %s:%lld-%lld candidates=%zu reassessed=%zu PASS=%zu\n", reg.chrom.c_str(),
                 static_cast<long long>(reg.start + 1), static_cast<long long>(reg.end), n_cand, n_reassessed, n_pass);
    std::fprintf(stderr, "thread-seconds: scan+evidence %.1f, engine (pattern+HMM+decision) %.1f\n",
                 t_scan_us.load() / 1e6, t_engine_us.load() / 1e6);
    if (truth)
        std::fprintf(stderr, "truth: %d simple indels in BED/region, %d complex alleles skipped\n", truth->n_in_bed,
                     truth->n_complex);
    return 0;
}
#endif

void usage() {
    std::fprintf(stderr,
                 "nindels-call <command> [options]\n\n"
                 "  call          --bam B --ref FA --region chr[:s-e] [--memory M] [--model F] [--base-model S]\n"
                 "                [--A 0.5] [--tau-strict 0.6] [--tau-loose 0.2] [--tau-call 0.5] [--min-qual 20]\n"
                 "                [--out-vcf V] [--out-table T] [--truth VCF --bed BED] [--threads N]\n"
                 "                [--candidates-vcf V (shadow mode)] [--force-reassess] [--parity even|odd]\n"
                 "                [--min-alt 2] [--min-alt-frac 0.03] [--max-hmm-reads 250]\n"
                 "  build-memory  --table T [--table T2 ...] --out M [--parity even|odd|all] [--theta 0.15] [--top-k 32]\n"
                 "  fit           --table T [--table ...] --out F (--preset P | --features a,b,c) [--lambda 0.01]\n"
                 "                [--mask-reassess 0.5]   presets: full full-neff-loci no-memory dl-orig stage1\n"
                 "                stage1-no-memory evidence-only spec-x\n"
                 "  select-tau    --table T   (scored, labelled validation table)\n");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];
    std::string cmdline;
    for (int i = 0; i < argc; ++i) cmdline += (i ? " " : "") + std::string(argv[i]);
    try {
        Args a = Args::parse(argc, argv, 2, {"force-reassess"});
        if (cmd == "build-memory") return cmd_build_memory(a);
        if (cmd == "fit") return cmd_fit(a);
        if (cmd == "select-tau") return cmd_select_tau(a);
#ifdef NINDELS_HAVE_HTSLIB
        if (cmd == "call") return cmd_call(a, cmdline);
#else
        if (cmd == "call") { std::fprintf(stderr, "built without htslib: 'call' unavailable\n"); return 2; }
#endif
        usage();
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
