#include "dnav2/engine.hpp"
#include <sys/resource.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include "dnav2/binomial.hpp"
#include "dnav2/extract.hpp"
#include "dnav2/features.hpp"
#include "dnav2/frame.hpp"
#include "dnav2/methodc.hpp"
#include "dnav2/model.hpp"
#include "dnav2/pb.hpp"

#ifdef DNAV2_PROFILE
#include <time.h>
#define PROF_T(v) uint64_t v = ([]{ timespec t; clock_gettime(CLOCK_MONOTONIC_RAW, &t); return uint64_t(t.tv_sec) * 1000000000ull + t.tv_nsec; })()
#define PROF_ADD(field, v) (dnav2::prof().field += (v))
#else
#define PROF_T(v) (void)0
#define PROF_ADD(field, v) (void)0
#endif

namespace dnav2 {
namespace {

using Clock = std::chrono::steady_clock;
inline double secs(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

const char* REF_TOKENS = "NACGT";

struct Task { size_t w0, nw; };

struct TaskResult {
    bool done = false;
    std::string vcf_cascade, vcf_pbonly, routed_tsv, dump_tsv, feat_tsv;
    int64_t n_cand = 0, n_pb_model = 0, n_model_calls = 0, n_model_changed = 0;
    std::vector<double> binom_llr, pb_llr;    // only when dump_llr_prefix
    int64_t n_loci = 0, n_frame = 0, n_routed = 0, n_casc = 0, n_pbo = 0, n_casc_forced = 0, n_pbo_forced = 0, n_pb_eval = 0;
    double t_count = 0, t_ref = 0, t_bin = 0, t_rows = 0, t_pb = 0, t_gt = 0;
    std::string error;
};

struct Shared {
    const Config* cfg;
    const std::vector<Window>* windows;
    const std::vector<uint8_t>* labels;
    std::vector<Task> tasks;
    std::vector<TaskResult> results;
    std::atomic<size_t> next{0};
    std::mutex mu;
    size_t next_flush = 0;
    // outputs
    FILE *f_casc = nullptr, *f_pbo = nullptr, *f_routed = nullptr, *f_dump = nullptr, *f_feat = nullptr, *f_bin = nullptr, *f_pbb = nullptr;
    const Model* model = nullptr;
    const std::vector<int64_t>* loci;   // sorted force/dump loci (may be empty)
    Summary sum;
    std::string error;
};

std::string hexf(double v) { char b[64]; std::snprintf(b, sizeof b, "%a", v); return b; }

struct RecordAcc { double k, n; int64_t pos0, depth; char ref, alt; };

void genotype_lines(const std::vector<RecordAcc>& recs, const Config& cfg, std::string& out, int64_t& forced) {
    for (const auto& r : recs) {
        MethodCResult m = method_c(r.k, r.n);
        if (m.forced) forced++;
        char line[256];
        std::snprintf(line, sizeof line, "%s\t%lld\t.\t%c\t%c\t.\tPASS\tDP=%lld\tGT:GQ\t%s:%d\n", cfg.contig.c_str(),
                      static_cast<long long>(r.pos0 + 1), r.ref, r.alt, static_cast<long long>(r.depth),
                      gt_string(m.gt), m.gq);
        out += line;
    }
}

struct PbCtx {
    int64_t ts;
    const std::vector<uint8_t>* alt;
    std::vector<double>* pb;
    std::vector<int32_t>* n_counted;
    std::vector<int32_t>* k_tensor;
    std::vector<uint8_t>* evaluated;
    std::unordered_map<int64_t, std::vector<uint8_t>>* rows_keep;  // dump loci -> 48x8 bytes
    const std::vector<int64_t>* loci;
    bool keep_rows;
    int64_t n_eval = 0;
};

void pb_emit(void* vctx, int64_t pos, const RowEntry* kept, int n_kept, int) {
    auto* c = static_cast<PbCtx*>(vctx);
    const size_t li = static_cast<size_t>(pos - c->ts);
    uint8_t code[MAX_READS], bq[MAX_READS];
    int m = 0;
    for (int i = 0; i < n_kept; i++)
        if (kept[i].row[6] & 2) { code[m] = kept[i].row[0]; bq[m] = kept[i].row[1]; m++; }
    const int alt = (*c->alt)[li];
    int k = 0;
    for (int i = 0; i < m; i++) k += (code[i] == alt);
    (*c->pb)[li] = pb_native_v2(code, bq, m, alt);
    (*c->n_counted)[li] = m;
    (*c->k_tensor)[li] = k;
    (*c->evaluated)[li] = 1;
    c->n_eval++;
    if (c->keep_rows && std::binary_search(c->loci->begin(), c->loci->end(), pos)) {
        std::vector<uint8_t> rows(MAX_READS * 8, 0);
        for (int j = 0; j < MAX_READS; j++) rows[j * 8] = 255;
        for (int j = 0; j < n_kept; j++) std::memcpy(&rows[j * 8], kept[j].row, 8);
        (*c->rows_keep)[pos] = std::move(rows);
    }
}

constexpr double FEAT_PB_LO = -5.0, FEAT_PB_HI = 60.0;   // PB@m dumped only for loci where PB could plausibly matter

struct CandRaw { int64_t i = 0; ReadStats rs; PbStore store; bool has_pb = false; };

// local-context stats of locus i from the task-local per-column arrays (counts only; no truth, no frame)
CtxStats compute_ctx(int64_t i, int64_t n, const std::vector<CountRow>& cnt, const std::vector<double>& kk,
                     const std::vector<double>& nn, const std::vector<uint8_t>& ref, int alt) {
    CtxStats c;
    auto cand = [&](int64_t j) { return kk[j] >= 3 && nn[j] > 0 && kk[j] / nn[j] >= 0.05; };
    auto nonref = [&](int64_t j) { const int r = ref[j]; return r >= 1 ? nn[j] - cnt[j].v[r - 1] : 0.0; };
    const int64_t lo25 = std::max<int64_t>(0, i - 25), hi25 = std::min<int64_t>(n - 1, i + 25);
    double snr = 0, sn = 0, sd = 0; int nd = 0;
    for (int64_t j = lo25; j <= hi25; j++) {
        if (j == i) continue;
        const int64_t d = j > i ? j - i : i - j;
        snr += nonref(j); sn += nn[j]; sd += cnt[j].v[6]; nd++;
        if (cand(j)) {
            c.nc25++;
            if (d <= 10) c.nc10++;
            if (d <= 5) c.nc5++;
            if (d <= 2) c.nc2++;
            if (d < c.nearest) c.nearest = static_cast<int>(d);
            const double v = kk[j] / nn[j];
            if (v > c.neigh_max_vaf25) c.neigh_max_vaf25 = v;
        }
    }
    c.lmr25 = sn > 0 ? snr / sn : 0.0;
    c.depth_ratio25 = (nd > 0 && sd > 0) ? cnt[i].v[6] / (sd / nd) : 1.0;
    double g = 0, dd = 0;
    for (int64_t j = std::max<int64_t>(0, i - 10); j <= std::min<int64_t>(n - 1, i + 10); j++) {
        if (j == i) continue;
        g += cnt[j].v[4] + cnt[j].v[5]; dd += cnt[j].v[6];
    }
    c.lgap10 = dd > 0 ? g / dd : 0.0;
    const int r = ref[i];
    int64_t a = i, b = i;
    while (a > 0 && ref[a - 1] == r && i - a < 30) a--;
    while (b + 1 < n && ref[b + 1] == r && b - i < 30) b++;
    c.hp_len = static_cast<int>(b - a + 1);
    c.left_is_alt = (i > 0 && ref[i - 1] == alt + 1);
    c.right_is_alt = (i + 1 < n && ref[i + 1] == alt + 1);
    for (int64_t j = std::max<int64_t>(0, i - 10); j <= std::min<int64_t>(n - 3, i + 8); j++)
        c.dinuc10 += (ref[j] != 0 && ref[j] == ref[j + 2]);
    c.is_ts = r >= 1 && (((r - 1) ^ alt) == 2);
    c.ref_cg = (r == 2 || r == 3);
    c.cpg = (r == 2 && i + 1 < n && ref[i + 1] == 3) || (r == 3 && i > 0 && ref[i - 1] == 2);
    return c;
}

std::string process_task(Shared& S, const Task& t, BamHandle& bam, Fai* fai, std::vector<RowEntry>& scratch, TaskResult& R) {
    const Config& cfg = *S.cfg;
    const auto& W = *S.windows;
    const int64_t ts = W[t.w0].start;
    const int64_t n = static_cast<int64_t>(t.nw) * SEQ_LEN;
    const int64_t te = ts + n;
    Filters filt;
    Clock::time_point a, b;

    std::vector<CountRow> counts(static_cast<size_t>(n));
    std::memset(counts.data(), 0, counts.size() * sizeof(CountRow));
    std::vector<uint8_t> refidx(static_cast<size_t>(n));
    a = Clock::now();
    // windows inside a task are contiguous
    if (auto e = fetch_ref_idx(fai, cfg.contig, ts, te, refidx.data()); !e.empty()) return e;
    b = Clock::now(); R.t_ref += secs(a, b);

    std::vector<double> bllr(n, 0.0), kk(n, 0.0), nn(n, 0.0);
    std::vector<uint8_t> alt(n, 0), frame(n, 0), routed(n, 0), scored(n, 0);
    std::vector<double> pb(n, 0.0);
    std::vector<int32_t> ncnt(n, 0), ktens(n, 0);
    std::vector<uint8_t> evaluated(n, 0);
    std::unordered_map<int64_t, std::vector<uint8_t>> rows_keep;
    const bool dump = !cfg.dump_loci.empty();
    std::vector<CandRaw> cands_out;

    auto score_locus = [&](int64_t i) {
        const CountRow& c = counts[i];
        double c4[4] = {double(c.v[0]), double(c.v[1]), double(c.v[2]), double(c.v[3])};
        BinomialScore s = binomial_score(c4, refidx[i]);
        bllr[i] = s.llr; alt[i] = static_cast<uint8_t>(s.alt_index); kk[i] = s.k; nn[i] = s.n;
        frame[i] = (*S.labels)[static_cast<size_t>(ts + i)] <= LABEL_SNP;
        routed[i] = frame[i] && is_routed(s.llr);
        scored[i] = 1;
    };
    // force/dump loci bitmask (per locus)
    std::vector<uint64_t> force(t.nw, 0);
    if (!S.loci->empty()) {
        auto lo = std::lower_bound(S.loci->begin(), S.loci->end(), ts);
        for (; lo != S.loci->end() && *lo < te; ++lo) force[(*lo - ts) / SEQ_LEN] |= 1ULL << ((*lo - ts) % SEQ_LEN);
    }
    auto is_forced = [&](int64_t i) { return (force[i / SEQ_LEN] >> (i % SEQ_LEN)) & 1ULL; };

    PbCtx ctx{ts, &alt, &pb, &ncnt, &ktens, &evaluated, &rows_keep, S.loci, dump};

    if (!cfg.exact_windows) {
        // ---- FUSED: one streaming pileup; counts, Binomial, router and (where needed) PB per column ----
        struct Fused {
            const Config* cfg; int64_t ts; const Filters* filt;
            std::vector<CountRow>* counts; std::function<void(int64_t)>* score; std::function<bool(int64_t)>* forced;
            std::vector<uint8_t>* routed; std::vector<RowEntry>* scratch; PbCtx* pb; const char* err = nullptr;
            std::vector<KeyEntry> ke; RowEntry sel[MAX_READS]; bool lazy_hash = false;
            double t_pb = 0;
            std::vector<CandRaw>* cands = nullptr; PbScratch* pbs = nullptr; int kmin = 3; bool feat_pb = false;
            const std::vector<uint8_t>* frame = nullptr; const std::vector<double>* kk = nullptr; const std::vector<uint8_t>* alt = nullptr;
            const std::vector<uint8_t>* refidx = nullptr; const std::vector<double>* bllr = nullptr; bool model_on = false;
        } F{&cfg, ts, &filt, &counts, nullptr, nullptr, &routed, &scratch, &ctx};
        std::function<void(int64_t)> f_score = score_locus;
        std::function<bool(int64_t)> f_forced = is_forced;
        F.score = &f_score; F.forced = &f_forced;
        PbScratch pbs;
        if (!cfg.dump_features.empty() || S.model) {
            F.model_on = S.model != nullptr;
            F.cands = &cands_out; F.pbs = &pbs; F.kmin = cfg.feat_kmin; F.feat_pb = cfg.feat_pb;
            F.frame = &frame; F.kk = &kk; F.alt = &alt; F.refidx = &refidx; F.bllr = &bllr;
        }
        F.lazy_hash = cfg.mode == Mode::CASCADE;   // few columns need rows: hash on demand instead of per read
        a = Clock::now();
        auto col = [](void* v, int64_t pos, const bam_pileup1_t* pil, int n_plp) {
            Fused& f = *static_cast<Fused*>(v);
            if (f.err) return;
            const int64_t i = pos - f.ts;
            PROF_T(t0);
            (*f.counts)[i] = count_column(pil, n_plp, *f.filt);
            PROF_T(t1); PROF_ADD(t_count_body, t1 - t0);
            (*f.score)(i);
            PROF_T(t2); PROF_ADD(t_build, t2 - t1);
            if (f.cands && (*f.frame)[i] && (*f.kk)[i] >= f.kmin) {
                CandRaw cr;
                cr.i = i;
                const int rb = (*f.refidx)[i] - 1;
                read_stats(pil, n_plp, *f.filt, (*f.alt)[i], rb, 32, cr.rs);
                const double bl = (*f.bllr)[i];
                if (f.feat_pb && bl >= FEAT_PB_LO && bl <= FEAT_PB_HI) {
                    cr.has_pb = pb_store_build(pil, n_plp, *f.filt, (*f.alt)[i], f.lazy_hash, *f.pbs, cr.store);
                }
                f.cands->push_back(cr);
            }
            const bool need = f.cfg->mode == Mode::PB_ALL || ((*f.routed)[i] && !(f.model_on && (*f.kk)[i] >= f.kmin)) || (*f.forced)(i);
            if (!need) return;
            int n = 0, keep;
            const RowEntry* kept;
            if (f.cfg->reference_rows) {
                n = build_rows(pil, n_plp, *f.filt, f.pb->keep_rows, *f.scratch, &f.err);
                PROF_T(t3); PROF_ADD(t_hash, t3 - t2);          // "t_hash" slot reused: row build
                if (f.err || n <= 0) return;
                keep = select_kept(*f.scratch, n, MAX_READS);
                kept = f.scratch->data();
            } else {
                keep = select_rows_lazy(pil, n_plp, *f.filt, f.pb->keep_rows, f.lazy_hash, MAX_READS, f.ke, f.sel, &n, &f.err);
                if (f.err || keep <= 0) return;
                kept = f.sel;
            }
            PROF_T(t4); PROF_ADD(t_select, t4 - t2);
            pb_emit(f.pb, pos, kept, keep, n);
            PROF_T(t5); PROF_ADD(t_emit, t5 - t4);
            PROF_ADD(n_loci_emitted, 1);
        };
        PROF_T(pass0);
        if (const char* e = pileup_pass(bam, ts, te, col, &F, !F.lazy_hash)) return e;
        if (F.err) return F.err;
        PROF_T(pass1); PROF_ADD(t_mplp, pass1 - pass0);   // whole pass incl. callbacks; callbacks subtracted in report
        b = Clock::now(); R.t_count += secs(a, b);   // fused: one pass; not split by stage
        for (int64_t i = 0; i < n; i++) if (!scored[i]) score_locus(i);   // uncovered columns: zero counts
    } else {
        // ---- V1-literal: separate counts pass, then a FRESH pileup per 64-bp window ----
        a = Clock::now();
        if (const char* e = count_pass(bam, ts, te, filt, counts.data())) return e;
        b = Clock::now(); R.t_count += secs(a, b);
        a = Clock::now();
        for (int64_t i = 0; i < n; i++) score_locus(i);
        b = Clock::now(); R.t_bin += secs(a, b);
        std::vector<uint64_t> want(t.nw, 0);
        if (cfg.mode == Mode::PB_ALL) std::fill(want.begin(), want.end(), ~0ULL);
        else for (int64_t i = 0; i < n; i++) if (routed[i]) want[i / SEQ_LEN] |= 1ULL << (i % SEQ_LEN);
        for (size_t w = 0; w < t.nw; w++) want[w] |= force[w];
        a = Clock::now();
        for (size_t w = 0; w < t.nw; w++) {
            if (!want[w]) continue;
            if (const char* e = window_pass(bam, ts + int64_t(w) * SEQ_LEN, ts + int64_t(w + 1) * SEQ_LEN, filt, MAX_READS,
                                            dump, want[w], scratch, pb_emit, &ctx)) return e;
        }
        b = Clock::now(); R.t_rows += secs(a, b);
    }
    R.n_loci = n;
    for (int64_t i = 0; i < n; i++) { R.n_frame += frame[i]; R.n_routed += routed[i]; }
    R.n_pb_eval = ctx.n_eval;

    std::vector<uint8_t> xdec(n, 0), xcall(n, 0);
    if (!cfg.dump_features.empty() || S.model) {
        char buf[64];
        std::vector<double> dpv;
        for (const CandRaw& cr : cands_out) {
            const int64_t i = cr.i;
            const CountRow& c = counts[i];
            CtxStats cx = compute_ctx(i, n, counts, kk, nn, refidx, alt[i]);
            double fv[NFEAT];
            derive_features(cr.rs, static_cast<int>(c.v[6]), static_cast<int>(nn[i]), static_cast<int>(kk[i]),
                            cr.rs.n_ref, cr.rs.n_other, bllr[i], cx, static_cast<int>(c.v[4]), static_cast<int>(c.v[5]), fv);
            if (S.model) {
                const bool call = S.model->decide(fv, routed[i], bllr[i], cr.has_pb ? &cr.store : nullptr, dpv, &R.n_pb_model);
                xdec[i] = 1; xcall[i] = call;
            }
            if (cfg.dump_features.empty()) continue;
            std::snprintf(buf, sizeof buf, "%lld\t%d\t%d\t%d\t%d\t%d", (long long)(ts + i), refidx[i], alt[i],
                          int(routed[i]), cr.has_pb ? cr.store.n_adm : -1, int((*S.labels)[ts + i]));
            R.feat_tsv += buf;
            for (int j = 0; j < NFEAT; j++) { std::snprintf(buf, sizeof buf, "\t%.17g", fv[j]); R.feat_tsv += buf; }
            for (int j = 0; j < NPB + 3 + 1; j++) {
                double v = 0;
                if (!cr.has_pb) { R.feat_tsv += "\tnan"; continue; }
                if (j < NPB) v = pb_store_eval(cr.store, j, dpv);
                else if (j < NPB + 3) v = pb_store_eval(cr.store, PBV_BLK_MED + (j - NPB), dpv);   // med, min, max
                else { int nb = 0; for (int b0 = 0; b0 + 48 <= cr.store.n_adm || b0 == 0; b0 += 48) { nb++; if (b0 + 48 >= cr.store.n_adm) break; } v = cr.store.n_adm == 0 ? 0.0 : nb; }
                std::snprintf(buf, sizeof buf, "\t%.17g", v); R.feat_tsv += buf;
            }
            R.feat_tsv += '\n';
        }
        R.n_cand = static_cast<int64_t>(cands_out.size());
    }

    // ---- calls + genotype ----------------------------------------------------------------
    a = Clock::now();
    std::vector<RecordAcc> rc, rp;
    auto make_record = [&](int64_t i, std::vector<RecordAcc>& v) {
        const CountRow& c = counts[i];
        int ri = refidx[i];
        if (ri == 0) return;                          // ref == N -> skipped (V1 extract_records)
        // Method C ALT: most-supported non-reference base, first maximum in A,C,G,T order, must be > 0
        double bc[4] = {double(c.v[0]), double(c.v[1]), double(c.v[2]), double(c.v[3])};
        int best = -1; double bv = 0;
        for (int j = 0; j < 4; j++) {
            if (j == ri - 1) continue;
            if (best < 0 || bc[j] > bv) { best = j; bv = bc[j]; }
        }
        if (best < 0 || bv <= 0) return;
        double k = bv, ref_n = bc[ri - 1];
        v.push_back({k, k + ref_n, ts + i, static_cast<int64_t>(c.v[6]), REF_TOKENS[ri], REF_TOKENS[best + 1]});
    };
    std::vector<uint8_t> ccall(n, 0), pcall(n, 0);
    for (int64_t i = 0; i < n; i++) {
        if (!frame[i]) continue;
        bool bcall = bllr[i] >= FROZEN_BINOMIAL_THRESHOLD;
        bool pc = evaluated[i] ? (pb[i] >= FROZEN_PB_THRESHOLD) : false;
        if (cfg.mode == Mode::PB_ALL) pcall[i] = pc;
        ccall[i] = routed[i] ? pc : bcall;
        if (xdec[i]) ccall[i] = xcall[i];
        if (ccall[i]) make_record(i, rc);
        if (cfg.mode == Mode::PB_ALL && pcall[i]) make_record(i, rp);
    }
    if (S.model) {
        R.n_pb_eval += R.n_pb_model;
        for (int64_t i = 0; i < n; i++) if (xdec[i]) { R.n_model_calls += xcall[i]; R.n_model_changed += (xcall[i] != (routed[i] ? (evaluated[i] && pb[i] >= FROZEN_PB_THRESHOLD) : (bllr[i] >= FROZEN_BINOMIAL_THRESHOLD))); }
    }
    genotype_lines(rc, cfg, R.vcf_cascade, R.n_casc_forced);
    R.n_casc = static_cast<int64_t>(rc.size());
    if (cfg.mode == Mode::PB_ALL) {
        genotype_lines(rp, cfg, R.vcf_pbonly, R.n_pbo_forced);
        R.n_pbo = static_cast<int64_t>(rp.size());
    }
    b = Clock::now(); R.t_gt += secs(a, b);

    // ---- auxiliary outputs -------------------------------------------------------------------
    for (int64_t i = 0; i < n; i++) {
        if (!routed[i]) continue;
        char line[512];
        const CountRow& c = counts[i];
        std::snprintf(line, sizeof line, "%lld\t%d\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%d\t%d\t%d\n", (long long)(ts + i),
                      refidx[i], c.v[0], c.v[1], c.v[2], c.v[3], c.v[6], hexf(bllr[i]).c_str(), hexf(pb[i]).c_str(),
                      int(bllr[i] >= FROZEN_BINOMIAL_THRESHOLD), int(pb[i] >= FROZEN_PB_THRESHOLD), int(ccall[i]));
        R.routed_tsv += line;
    }
    if (dump) {
        for (auto lo = std::lower_bound(S.loci->begin(), S.loci->end(), ts); lo != S.loci->end() && *lo < te; ++lo) {
            const int64_t i = *lo - ts;
            const CountRow& c = counts[i];
            std::ostringstream o;
            o << *lo << '\t' << int(refidx[i]);
            for (int j = 0; j < 9; j++) o << '\t' << c.v[j];
            o << '\t' << hexf(bllr[i]) << '\t' << int(alt[i]) << '\t' << static_cast<long long>(kk[i]) << '\t'
              << static_cast<long long>(nn[i]) << '\t' << int((*S.labels)[*lo]) << '\t' << int(frame[i]) << '\t'
              << int(routed[i]) << '\t' << int(bllr[i] >= FROZEN_BINOMIAL_THRESHOLD) << '\t' << int(evaluated[i]) << '\t'
              << hexf(pb[i]) << '\t' << ncnt[i] << '\t' << ktens[i] << '\t' << int(pb[i] >= FROZEN_PB_THRESHOLD) << '\t'
              << int(ccall[i]) << '\t';
            auto it = rows_keep.find(*lo);
            if (it == rows_keep.end()) o << '-';
            else { static const char* H = "0123456789abcdef"; for (uint8_t v : it->second) o << H[v >> 4] << H[v & 15]; }
            o << '\n';
            R.dump_tsv += o.str();
        }
    }
    if (!cfg.dump_llr_prefix.empty()) { R.binom_llr = bllr; R.pb_llr = pb; }
    return "";
}

void flush_ready(Shared& S) {   // caller holds S.mu
    while (S.next_flush < S.results.size() && S.results[S.next_flush].done) {
        TaskResult& r = S.results[S.next_flush];
        if (S.f_casc) std::fwrite(r.vcf_cascade.data(), 1, r.vcf_cascade.size(), S.f_casc);
        if (S.f_pbo) std::fwrite(r.vcf_pbonly.data(), 1, r.vcf_pbonly.size(), S.f_pbo);
        if (S.f_routed) std::fwrite(r.routed_tsv.data(), 1, r.routed_tsv.size(), S.f_routed);
        if (S.f_feat) std::fwrite(r.feat_tsv.data(), 1, r.feat_tsv.size(), S.f_feat);
        if (S.f_dump) std::fwrite(r.dump_tsv.data(), 1, r.dump_tsv.size(), S.f_dump);
        if (S.f_bin) std::fwrite(r.binom_llr.data(), 8, r.binom_llr.size(), S.f_bin);
        if (S.f_pbb) std::fwrite(r.pb_llr.data(), 8, r.pb_llr.size(), S.f_pbb);
        Summary& s = S.sum;
        s.n_loci += r.n_loci; s.n_frame += r.n_frame; s.n_routed += r.n_routed;
        s.n_cascade_records += r.n_casc; s.n_pbonly_records += r.n_pbo;
        s.n_cand += r.n_cand; s.n_model_calls += r.n_model_calls; s.n_model_changed += r.n_model_changed; s.n_pb_model += r.n_pb_model;
        s.n_cascade_forced += r.n_casc_forced; s.n_pbonly_forced += r.n_pbo_forced; s.n_pb_evaluated += r.n_pb_eval;
        s.t_count_pass += r.t_count; s.t_ref_labels += r.t_ref; s.t_binomial += r.t_bin;
        s.t_window_rows += r.t_rows; s.t_pb += r.t_pb; s.t_genotype += r.t_gt;
        // release the big buffers
        TaskResult empty; empty.done = true;
        r = std::move(empty);
        S.next_flush++;
    }
}

const char* HEADER_TMPL =
    "##fileformat=VCFv4.2\n"
    "##source=ai_dna_analyzer_%s_%s_C\n"
    "##contig=<ID=%s,length=%lld>\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth (passing-filter reads)\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (frozen Method C: binomial genotype posterior argmax over {0/0,0/1,1/1}, eps=0.01 fixed)\">\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality (10*log10(P_best/P_second), capped at 99)\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n";

}  // namespace

std::string run(const Config& cfg, Summary& out) {
    auto t_start = Clock::now();
    Shared S;
    S.cfg = &cfg;
    std::vector<Window> windows;
    std::vector<uint8_t> labels;
    std::string e;
    if (!(e = build_windows(cfg.fasta, cfg.bed, cfg.contig, cfg.contig_len, windows)).empty()) return e;
    if (!(e = load_labels(cfg.truth, cfg.contig, cfg.contig_len, labels)).empty()) return e;

    if (cfg.skip_windows > 0 || cfg.limit_windows >= 0) {
        size_t a = std::min<size_t>(windows.size(), static_cast<size_t>(cfg.skip_windows));
        size_t b = cfg.limit_windows < 0 ? windows.size() : std::min<size_t>(windows.size(), a + cfg.limit_windows);
        windows = std::vector<Window>(windows.begin() + a, windows.begin() + b);
    }
    std::vector<int64_t> loci;
    if (!cfg.loci_file.empty()) {
        std::ifstream in(cfg.loci_file);
        int64_t v;
        while (in >> v) loci.push_back(v);
        std::sort(loci.begin(), loci.end());
        loci.erase(std::unique(loci.begin(), loci.end()), loci.end());
        // restrict to windows containing these loci
        std::vector<Window> sub;
        for (const auto& w : windows) {
            auto lo = std::lower_bound(loci.begin(), loci.end(), w.start);
            if (lo != loci.end() && *lo < w.start + SEQ_LEN) sub.push_back(w);
        }
        windows.swap(sub);
    }
    if (cfg.setup_only) {
        FILE* fw = std::fopen((cfg.out_dir + "/windows_" + cfg.tag + ".i64").c_str(), "wb");
        FILE* fl = std::fopen((cfg.out_dir + "/labels_" + cfg.tag + ".u8").c_str(), "wb");
        if (!fw || !fl) return "cannot write setup dumps";
        for (const auto& w : windows) {
            int64_t st = w.start;
            std::fwrite(&st, 8, 1, fw);
            std::fwrite(&labels[static_cast<size_t>(st)], 1, SEQ_LEN, fl);
        }
        std::fclose(fw); std::fclose(fl);
        out.n_windows = static_cast<int64_t>(windows.size());
        out.n_loci = out.n_windows * SEQ_LEN;
        out.t_setup = secs(t_start, Clock::now());
        return "";
    }
    Model model;
    if (!cfg.model_file.empty()) {
        if (!(e = model.load(cfg.model_file)).empty()) return e;
        if (cfg.mode == Mode::PB_ALL) return "--model is cascade-mode only";
        S.model = &model;
    }
    S.windows = &windows; S.labels = &labels; S.loci = &loci;
    S.sum.n_windows = static_cast<int64_t>(windows.size());

    // tasks: runs of consecutive windows, at most task_windows each
    for (size_t i = 0; i < windows.size();) {
        size_t j = i + 1;
        while (j < windows.size() && j - i < static_cast<size_t>(cfg.task_windows) &&
               windows[j].start == windows[j - 1].start + SEQ_LEN) j++;
        S.tasks.push_back({i, j - i});
        i = j;
    }
    S.results.resize(S.tasks.size());

    auto open_out = [&](const std::string& name, FILE*& f) -> bool {
        f = std::fopen((cfg.out_dir + "/" + name).c_str(), "wb");
        return f != nullptr;
    };
    if (cfg.write_vcf) {
        if (!open_out("ai_cascade_" + cfg.tag + ".vcf", S.f_casc)) return "cannot write cascade VCF";
        std::fprintf(S.f_casc, HEADER_TMPL, cfg.tag.c_str(), "cascade", cfg.contig.c_str(), (long long)cfg.contig_len, cfg.sample.c_str());
        if (cfg.mode == Mode::PB_ALL) {
            if (!open_out("ai_pb_only_" + cfg.tag + ".vcf", S.f_pbo)) return "cannot write pb-only VCF";
            std::fprintf(S.f_pbo, HEADER_TMPL, cfg.tag.c_str(), "pb_only", cfg.contig.c_str(), (long long)cfg.contig_len, cfg.sample.c_str());
        }
    }
    if (!open_out("routed_loci_" + cfg.tag + ".tsv", S.f_routed)) return "cannot write routed table";
    std::fprintf(S.f_routed, "pos0\tref_idx\tcA\tcC\tcG\tcT\tdepth\tbinomial_llr_hex\tpb_llr_hex\tbinomial_call\tpb_call\tcascade_call\n");
    if (!cfg.dump_loci.empty()) {
        S.f_dump = std::fopen(cfg.dump_loci.c_str(), "wb");
        if (!S.f_dump) return "cannot write dump";
        std::fprintf(S.f_dump, "pos0\tref_idx\tc0\tc1\tc2\tc3\tc4\tc5\tc6\tc7\tc8\tbinomial_llr_hex\talt\tk\tn\tlabel\tframe\trouted\t"
                               "binomial_call\tpb_evaluated\tpb_llr_hex\tn_counted_tensor\tk_tensor\tpb_call\tcascade_call\trows_hex\n");
    }
    if (!cfg.dump_features.empty()) {
        S.f_feat = std::fopen(cfg.dump_features.c_str(), "wb");
        if (!S.f_feat) return "cannot write feature dump";
        std::fprintf(S.f_feat, "pos0\tref_idx\talt\trouted\tn_adm\tlabel");
        for (int j = 0; j < NFEAT; j++) std::fprintf(S.f_feat, "\t%s", FEAT_NAMES[j]);
        for (int j = 0; j < NPB; j++) std::fprintf(S.f_feat, "\tpb%d", PB_CAPS[j]);
        std::fprintf(S.f_feat, "\tpbblk_med\tpbblk_min\tpbblk_max\tpbblk_n\n");
    }
    if (!cfg.dump_llr_prefix.empty()) {
        S.f_bin = std::fopen((cfg.dump_llr_prefix + ".binom.f64").c_str(), "wb");
        S.f_pbb = std::fopen((cfg.dump_llr_prefix + ".pb.f64").c_str(), "wb");
        if (!S.f_bin || !S.f_pbb) return "cannot write llr dumps";
    }
    out = S.sum;
    out.t_setup = secs(t_start, Clock::now());

    const int nthreads = std::max(1, cfg.threads);
    std::vector<std::thread> pool;
    for (int ti = 0; ti < nthreads; ti++) {
        pool.emplace_back([&]() {
            BamHandle bam;
            std::string err = bam.open(cfg.bam, cfg.contig, cfg.decomp_threads);
            Fai* fai = fai_open(cfg.fasta);
            if (err.empty() && !fai) err = "cannot open FASTA in worker";
            std::vector<RowEntry> scratch(4096);
            for (;;) {
                size_t ix = S.next.fetch_add(1);
                if (ix >= S.tasks.size()) break;
                TaskResult r;
                if (err.empty()) err = process_task(S, S.tasks[ix], bam, fai, scratch, r);
                r.done = true;
                std::lock_guard<std::mutex> g(S.mu);
                if (!err.empty() && S.error.empty()) S.error = err;
                S.results[ix] = std::move(r);
                flush_ready(S);
                if (!err.empty()) break;
            }
            fai_close(fai);
        });
    }
    for (auto& th : pool) th.join();
    for (FILE* f : {S.f_casc, S.f_pbo, S.f_routed, S.f_dump, S.f_feat, S.f_bin, S.f_pbb}) if (f) std::fclose(f);
    if (!S.error.empty()) return S.error;

    out = S.sum;
    out.threads = nthreads;
    out.t_setup = secs(t_start, t_start) + 0;   // filled below
    out.t_wall = secs(t_start, Clock::now());
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    out.t_cpu_user = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6;
    out.t_cpu_sys = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1e6;
    out.peak_rss_kb = ru.ru_maxrss;
    return "";
}

std::string summary_json(const Config& c, const Summary& s) {
    char b[2048];
    std::snprintf(b, sizeof b,
        "{\n  \"tag\": \"%s\", \"mode\": \"%s\", \"threads\": %d, \"decomp_threads\": %d, \"task_windows\": %d,\n"
        "  \"n_windows\": %lld, \"n_loci\": %lld, \"n_frame_scored\": %lld, \"n_routed_pb\": %lld,\n"
        "  \"n_pb_evaluated\": %lld, \"n_candidates\": %lld, \"n_pb_model_evals\": %lld, \"n_model_calls\": %lld, \"n_model_changed_vs_v2\": %lld, \"model\": \"%s\",\n"
        "  \"cascade\": {\"n_calls\": %lld, \"n_forced_00_to_01\": %lld},\n"
        "  \"pb_only\": {\"n_calls\": %lld, \"n_forced_00_to_01\": %lld},\n"
        "  \"cpu_seconds_by_stage\": {\"count_pass\": %.3f, \"ref_fetch\": %.3f, \"binomial_router\": %.3f, "
        "\"window_pass_and_pb\": %.3f, \"genotype_vcf_text\": %.3f},\n"
        "  \"wall_seconds\": %.3f, \"cpu_user_seconds\": %.3f, \"cpu_sys_seconds\": %.3f, \"peak_rss_kb\": %lld\n}\n",
        c.tag.c_str(), c.mode == Mode::CASCADE ? "cascade" : "pb_all", s.threads, c.decomp_threads, c.task_windows,
        (long long)s.n_windows, (long long)s.n_loci, (long long)s.n_frame, (long long)s.n_routed, (long long)s.n_pb_evaluated, (long long)s.n_cand, (long long)s.n_pb_model, (long long)s.n_model_calls, (long long)s.n_model_changed, c.model_file.c_str(),
        (long long)s.n_cascade_records, (long long)s.n_cascade_forced, (long long)s.n_pbonly_records, (long long)s.n_pbonly_forced,
        s.t_count_pass, s.t_ref_labels, s.t_binomial, s.t_window_rows, s.t_genotype, s.t_wall, s.t_cpu_user, s.t_cpu_sys,
        (long long)s.peak_rss_kb);
    return b;
}

}  // namespace dnav2
