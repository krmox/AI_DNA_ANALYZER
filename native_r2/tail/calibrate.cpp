#include "calibrate.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>
#include "dnav2/binomial.hpp"
#include "dnav2/extract.hpp"

namespace dnav2 {
namespace {

constexpr int DEPTH_BINS = 2048;
constexpr int SUBSAMPLE_READ_COLS = 8;   // read-level histograms (MAPQ/BQ) are filled on every 8th column only
constexpr int MIN_COL_DEPTH_FOR_EPS = 10;
constexpr double EPS_COL_MAX_VAF = 0.1;

struct Acc {
    const uint8_t* ref = nullptr;          // reference index of [start, end)
    int64_t start = 0;
    const std::vector<std::pair<int64_t, int64_t>>* clip = nullptr;
    Filters filt;
    std::vector<uint64_t> depth_hist = std::vector<uint64_t>(DEPTH_BINS, 0);
    std::vector<uint64_t> mq_hist = std::vector<uint64_t>(256, 0);
    std::vector<uint64_t> bq_hist = std::vector<uint64_t>(256, 0);
    std::vector<uint64_t> rl_hist = std::vector<uint64_t>(1024, 0);
    int64_t depth_sum = 0;
    int64_t cols = 0, eps_bases = 0, eps_nonref = 0, cand = 0, cand_low = 0;
    int64_t col_counter = 0;
};

void on_column(void* v, int64_t pos, const bam_pileup1_t* pil, int n_plp) {
    Acc& a = *static_cast<Acc*>(v);
    if (a.clip && !in_intervals(*a.clip, pos)) return;
    const CountRow c = count_column(pil, n_plp, a.filt);
    const int depth = static_cast<int>(c.v[6]);
    a.cols++;
    a.depth_sum += depth;
    a.depth_hist[static_cast<size_t>(std::min(depth, DEPTH_BINS - 1))]++;
    if ((a.col_counter++ % SUBSAMPLE_READ_COLS) == 0) {
        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t* p = pil + i;
            a.mq_hist[p->b->core.qual]++;
            if (!p->is_del && !p->is_refskip) {
                const uint8_t* q = bam_get_qual(p->b);
                if (q[0] != 0xff) a.bq_hist[q[p->qpos]]++;
            }
            const int rl = p->b->core.l_qseq;
            if (rl > 0) a.rl_hist[static_cast<size_t>(std::min(rl, 1023))]++;
        }
    }
    const int ri = a.ref[pos - a.start];
    if (ri < 1 || ri > 4) return;
    const double n = static_cast<double>(c.v[0]) + c.v[1] + c.v[2] + c.v[3];
    if (n < MIN_COL_DEPTH_FOR_EPS) return;
    const double nonref = n - c.v[ri - 1];
    double c4[4] = {double(c.v[0]), double(c.v[1]), double(c.v[2]), double(c.v[3])};
    const BinomialScore s = binomial_score(c4, ri);
    if (s.k >= 3 && s.k / s.n >= 0.05) {
        a.cand++;
        if (s.k / s.n < LOWVAF_CUT) a.cand_low++;
    }
    if (nonref / n < EPS_COL_MAX_VAF) { a.eps_bases += static_cast<int64_t>(n); a.eps_nonref += static_cast<int64_t>(nonref); }
}

}  // namespace

namespace {
void merge_into(Acc& dst, const Acc& src) {
    for (size_t i = 0; i < dst.depth_hist.size(); i++) dst.depth_hist[i] += src.depth_hist[i];
    for (size_t i = 0; i < dst.mq_hist.size(); i++) dst.mq_hist[i] += src.mq_hist[i];
    for (size_t i = 0; i < dst.bq_hist.size(); i++) dst.bq_hist[i] += src.bq_hist[i];
    for (size_t i = 0; i < dst.rl_hist.size(); i++) dst.rl_hist[i] += src.rl_hist[i];
    dst.depth_sum += src.depth_sum; dst.cols += src.cols; dst.eps_bases += src.eps_bases; dst.eps_nonref += src.eps_nonref;
    dst.cand += src.cand; dst.cand_low += src.cand_low;
}
}  // namespace

std::string calibrate(const std::string& bam, const std::string& fasta, const std::string& contig,
                      const std::vector<Window>& windows, const std::vector<std::pair<int64_t, int64_t>>* bed_clip,
                      int blocks, int block_windows, int threads, Calibration& out) {
    out = Calibration{};
    if (windows.empty() || blocks <= 0 || block_windows <= 0) return "calibration: nothing to sample";
    const size_t nw = windows.size();
    const size_t span = std::min<size_t>(static_cast<size_t>(block_windows), nw);
    const size_t nblocks = std::min<size_t>(static_cast<size_t>(blocks), nw / span ? nw / span : 1);
    // Every block's statistics are integer sums, so the merged result does not depend on the thread count or scheduling.
    const size_t nthreads = std::max<size_t>(1, std::min<size_t>(static_cast<size_t>(std::max(threads, 1)), nblocks));
    std::vector<Acc> accs(nthreads);
    std::vector<std::string> errs(nthreads);
    std::atomic<size_t> next{0};
    auto worker = [&](size_t t) {
        BamHandle bh;
        if (std::string e = bh.open(bam, contig); !e.empty()) { errs[t] = "calibration: " + e; return; }
        Fai* fai = fai_open(fasta);
        if (!fai) { errs[t] = "calibration: cannot open FASTA"; return; }
        Acc& a = accs[t];
        a.clip = bed_clip;
        std::vector<uint8_t> ref;
        for (size_t b = next.fetch_add(1); b < nblocks && errs[t].empty(); b = next.fetch_add(1)) {
            const size_t w0 = nblocks == 1 ? 0 : b * ((nw - span) / (nblocks - 1));
            size_t w1 = w0;
            while (w1 + 1 < nw && w1 + 1 - w0 < span && windows[w1 + 1].start == windows[w1].start + SEQ_LEN) w1++;
            const int64_t s = windows[w0].start, e = windows[w1].start + SEQ_LEN;
            ref.assign(static_cast<size_t>(e - s), 0);
            if (std::string r = fetch_ref_idx(fai, contig, s, e, ref.data()); !r.empty()) { errs[t] = "calibration: " + r; break; }
            a.ref = ref.data(); a.start = s;
            if (const char* r = pileup_pass(bh, s, e, on_column, &a, false)) errs[t] = std::string("calibration: ") + r;
        }
        fai_close(fai);
    };
    std::vector<std::thread> pool;
    for (size_t t = 1; t < nthreads; t++) pool.emplace_back(worker, t);
    worker(0);
    for (auto& th : pool) th.join();
    for (const std::string& e : errs) if (!e.empty()) return e;
    Acc a;
    for (const Acc& x : accs) merge_into(a, x);
    if (a.cols < 1000) return "calibration: too few covered columns";
    out.n_cols = a.cols;
    out.n_bases = a.eps_bases;
    out.depth_mean = static_cast<double>(a.depth_sum) / static_cast<double>(a.cols);
    out.depth_med = hist_quantile(a.depth_hist, 0.5);
    out.mq_p50 = hist_quantile(a.mq_hist, 0.5);
    out.mq_p90 = hist_quantile(a.mq_hist, 0.9);
    out.bq_p10 = hist_quantile(a.bq_hist, 0.1);
    out.bq_p50 = hist_quantile(a.bq_hist, 0.5);
    out.bq_p90 = hist_quantile(a.bq_hist, 0.9);
    out.read_len_med = static_cast<int>(hist_quantile(a.rl_hist, 0.5));
    if (a.eps_bases > 0)
        out.eps_hat = std::min(EPS_CLAMP_HI, std::max(EPS_CLAMP_LO, static_cast<double>(a.eps_nonref) / static_cast<double>(a.eps_bases)));
    out.cand_per_bp = static_cast<double>(a.cand) / static_cast<double>(a.cols);
    out.lowvaf_base = a.cand > 0 ? static_cast<double>(a.cand_low) / static_cast<double>(a.cand) : 0.0;
    out.valid = out.depth_med > 0 && out.mq_p90 > 0 && out.bq_p50 > 0;
    return out.valid ? "" : "calibration: degenerate statistics";
}

}  // namespace dnav2
