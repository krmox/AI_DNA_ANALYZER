#include "dnav2/extract.hpp"
#include <cmath>
#include <cstdlib>

#ifdef DNAV2_PROFILE
#include <time.h>
#include <atomic>
#include <cstdio>
namespace dnav2 {
static std::atomic<uint64_t> G[13];
static inline uint64_t now_ns() { timespec t; clock_gettime(CLOCK_MONOTONIC_RAW, &t); return uint64_t(t.tv_sec) * 1000000000ull + t.tv_nsec; }
struct ProfTL { ProfCounters c; ~ProfTL() { uint64_t* p = reinterpret_cast<uint64_t*>(&c); for (int i = 0; i < 13; i++) G[i] += p[i]; } };
ProfCounters& prof() { static thread_local ProfTL t; return t.c; }
void prof_dump(const char* path) {
    const char* names[13] = {"t_mplp_ns","t_build_ns","t_select_ns","t_emit_ns","t_hash_ns","n_hash","n_loci_emitted","n_reads_admitted","n_reads_seen","t_count_mplp_ns","t_count_body_ns","n_count_loci","n_count_reads"};
    FILE* f = std::fopen(path, "w");
    for (int i = 0; i < 13; i++) std::fprintf(f, "%s\t%llu\n", names[i], (unsigned long long)G[i].load());
    std::fclose(f);
}
}
#define PROF_NOW() dnav2::now_ns()
#else
#define PROF_NOW() 0
#endif

namespace dnav2 {

std::string BamHandle::open(const std::string& path, const std::string& contig, int n_dec) {
    fp = sam_open(path.c_str(), "r");
    if (!fp) return "cannot open BAM " + path;
    hts_set_cache_size(fp, 8 * 1024 * 1024);
    if (n_dec > 0) hts_set_threads(fp, n_dec);
    hdr = sam_hdr_read(fp);
    if (!hdr) return "cannot read BAM header";
    idx = sam_index_load(fp, path.c_str());
    if (!idx) return "cannot load BAM index";
    tid = sam_hdr_name2tid(hdr, contig.c_str());
    if (tid < 0) return "contig not in BAM: " + contig;
    return "";
}
BamHandle::~BamHandle() {
    if (idx) hts_idx_destroy(idx);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp) sam_close(fp);
}

namespace {
struct IterState { samFile* fp; hts_itr_t* iter; };

// identical read callback to V1 (pysam stepper="samtools" subset: default flag filter + orphans)
int read_cb(void* data, bam1_t* b) {
    auto* st = static_cast<IterState*>(data);
    int ret;
    for (;;) {
        ret = sam_itr_next(st->fp, st->iter, b);
        if (ret < 0) break;
        uint16_t flag = b->core.flag;
        if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if ((flag & BAM_FPAIRED) && !(flag & BAM_FPROPER_PAIR)) continue;
        break;
    }
    return ret;
}
int plp_construct(void*, const bam1_t* b, bam_pileup_cd* cd) {
    const char* qn = bam_get_qname(b);
#ifdef DNAV2_PROFILE
    uint64_t t0 = PROF_NOW();
#endif
    cd->i = static_cast<int64_t>(blake2b_key64(reinterpret_cast<const uint8_t*>(qn), std::strlen(qn)));
#ifdef DNAV2_PROFILE
    prof().t_hash += PROF_NOW() - t0; prof().n_hash++;
#endif
    return 0;
}
inline int clip255(long v) { return v > 255 ? 255 : static_cast<int>(v); }
}  // namespace

CountRow count_column(const bam_pileup1_t* pil, int n_plp, const Filters& f) {
    uint32_t c[4] = {0, 0, 0, 0};
    uint32_t gap = 0, ins = 0, depth = 0, qsum = 0, msum = 0;
    for (int i = 0; i < n_plp; i++) {
        const bam_pileup1_t* p = pil + i;
        bam1_t* bb = p->b;
        int mapq = bb->core.qual;
        if (mapq < f.min_mapq) continue;
        uint16_t flag = bb->core.flag;
        if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        depth++;
        msum += static_cast<uint32_t>(mapq);
        if (p->is_del || p->is_refskip) { gap++; continue; }
        int qpos = p->qpos;
        uint8_t* qual = bam_get_qual(bb);
        int bq = (qual[0] == 0xff) ? 0 : static_cast<int>(qual[qpos]);
        if (bq < f.min_baseq) continue;
        int col;
        switch (seq_nt16_str[bam_seqi(bam_get_seq(bb), qpos)]) {
            case 'A': col = 0; break; case 'C': col = 1; break;
            case 'G': col = 2; break; case 'T': col = 3; break; default: col = -1;
        }
        if (col < 0) continue;
        c[col]++;
        qsum += static_cast<uint32_t>(bq);
        if (p->indel > 0) ins++;
    }
    CountRow r;
    r.v[0] = c[0]; r.v[1] = c[1]; r.v[2] = c[2]; r.v[3] = c[3];
    r.v[4] = gap; r.v[5] = ins; r.v[6] = depth; r.v[7] = qsum; r.v[8] = msum;
    return r;
}

const char* count_pass(BamHandle& b, int64_t start, int64_t end, const Filters& f, CountRow* out) {
    hts_itr_t* iter = sam_itr_queryi(b.idx, b.tid, static_cast<hts_pos_t>(start), static_cast<hts_pos_t>(end));
    if (!iter) return "cannot create iterator";
    IterState st{b.fp, iter};
    void* data_arr[1] = {&st};
    bam_mplp_t mplp = bam_mplp_init(1, read_cb, data_arr);
    bam_mplp_set_maxcnt(mplp, 8000);
    bam_mplp_init_overlaps(mplp);
    int tid_out, pos, n_plp_arr[1];
    const bam_pileup1_t* pil_arr[1];
    while (bam_mplp_auto(mplp, &tid_out, &pos, n_plp_arr, pil_arr) > 0) {
        if (pos < start || pos >= end) continue;
        out[pos - start] = count_column(pil_arr[0], n_plp_arr[0], f);
    }
    bam_mplp_destroy(mplp);
    hts_itr_destroy(iter);
    return nullptr;
}

const char* pileup_pass(BamHandle& b, int64_t start, int64_t end, ColumnFn fn, void* ctx, bool hash_at_construct) {
    hts_itr_t* iter = sam_itr_queryi(b.idx, b.tid, static_cast<hts_pos_t>(start), static_cast<hts_pos_t>(end));
    if (!iter) return "cannot create iterator";
    IterState st{b.fp, iter};
    void* data_arr[1] = {&st};
    bam_mplp_t mplp = bam_mplp_init(1, read_cb, data_arr);
    bam_mplp_set_maxcnt(mplp, 8000);
    bam_mplp_init_overlaps(mplp);
    if (hash_at_construct) bam_mplp_constructor(mplp, plp_construct);
    int tid_out, pos, n_plp_arr[1];
    const bam_pileup1_t* pil_arr[1];
    while (bam_mplp_auto(mplp, &tid_out, &pos, n_plp_arr, pil_arr) > 0) {
        if (pos < start || pos >= end) continue;
        fn(ctx, pos, pil_arr[0], n_plp_arr[0]);
    }
    bam_mplp_destroy(mplp);
    hts_itr_destroy(iter);
    return nullptr;
}

int build_rows(const bam_pileup1_t* pil, int n_plp, const Filters& f, bool full_rows,
               std::vector<RowEntry>& ents, const char** err) {
    if (static_cast<size_t>(n_plp) > ents.size()) ents.resize(static_cast<size_t>(n_plp) * 2);
    int n = 0;
    for (int i = 0; i < n_plp; i++) {
        const bam_pileup1_t* p = pil + i;
        bam1_t* bb = p->b;
        int mapq = bb->core.qual;
        if (mapq < f.min_mapq) continue;
        uint16_t flag = bb->core.flag;
        if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        int length = bb->core.l_qseq;
        if (length <= 0) { *err = "l_qseq==0 (V1 cannot reproduce)"; return -1; }
        RowEntry* e = &ents[n];
        e->key = static_cast<uint64_t>(p->cd.i);
        e->idx = static_cast<uint32_t>(n);
        uint8_t* r = e->row;
        int flags = 1;
        if (full_rows) {
            if (p->indel > 0) flags |= 8; else if (p->indel < 0) flags |= 16;
        }
        long nm = 0;
        int md = 0;
        if (full_rows) {
            uint8_t* aux = bam_aux_get(bb, "NM");
            if (aux) {
                char t = static_cast<char>(*aux);
                if (t == 'c' || t == 'C' || t == 's' || t == 'S' || t == 'i' || t == 'I')
                    nm = static_cast<long>(bam_aux2i(aux));
                else { *err = "non-integer NM (V1 cannot reproduce)"; return -1; }
            }
            md = clip255(static_cast<long>(std::nearbyint(2550.0 * static_cast<double>(nm) /
                                                          static_cast<double>(length > 1 ? length : 1))));
        }
        int strand = (flag & BAM_FREVERSE) ? 1 : 0;
        if (p->is_del || p->is_refskip) {
            r[0] = 4; r[1] = 0; r[2] = static_cast<uint8_t>(mapq); r[3] = static_cast<uint8_t>(strand);
            r[4] = 0; r[5] = 0; r[6] = static_cast<uint8_t>(flags | 4); r[7] = static_cast<uint8_t>(md);
        } else {
            int qpos = p->qpos;
            uint8_t* qual = bam_get_qual(bb);
            int q = (qual[0] == 0xff) ? 0 : static_cast<int>(qual[qpos]);
            int code = 5;
            switch (seq_nt16_str[bam_seqi(bam_get_seq(bb), qpos)]) {
                case 'A': code = 0; break; case 'C': code = 1; break;
                case 'G': code = 2; break; case 'T': code = 3; break; default: break;
            }
            if (code <= 3 && q >= f.min_baseq) flags |= 2;
            r[0] = static_cast<uint8_t>(code); r[1] = static_cast<uint8_t>(q);
            r[2] = static_cast<uint8_t>(mapq); r[3] = static_cast<uint8_t>(strand);
            if (full_rows) {
                int rest = length - 1 - qpos; if (rest < 0) rest = 0;
                int end_dist = qpos < rest ? qpos : rest;
                if (end_dist < 5) flags |= 32;
                double half = static_cast<double>(length) / 2.0; if (half < 1.0) half = 1.0;
                r[4] = static_cast<uint8_t>(clip255(static_cast<long>(std::nearbyint(
                    255.0 * static_cast<double>(qpos) / static_cast<double>(length > 1 ? length : 1)))));
                r[5] = static_cast<uint8_t>(clip255(static_cast<long>(std::nearbyint(
                    255.0 * static_cast<double>(end_dist) / half))));
            } else { r[4] = r[5] = 0; }
            r[6] = static_cast<uint8_t>(flags); r[7] = static_cast<uint8_t>(md);
        }
        n++;
    }
    return n;
}

// one row exactly as build_rows builds it (same code path semantics)
static inline bool make_row(const bam_pileup1_t* p, const Filters& f, bool full_rows, uint8_t* r, const char** err) {
    bam1_t* bb = p->b;
    const int mapq = bb->core.qual;
    const int length = bb->core.l_qseq;
    int flags = 1;
    if (full_rows) {
        if (p->indel > 0) flags |= 8; else if (p->indel < 0) flags |= 16;
    }
    long nm = 0;
    int md = 0;
    if (full_rows) {
        uint8_t* aux = bam_aux_get(bb, "NM");
        if (aux) {
            char t = static_cast<char>(*aux);
            if (t == 'c' || t == 'C' || t == 's' || t == 'S' || t == 'i' || t == 'I')
                nm = static_cast<long>(bam_aux2i(aux));
            else { *err = "non-integer NM (V1 cannot reproduce)"; return false; }
        }
        md = clip255(static_cast<long>(std::nearbyint(2550.0 * static_cast<double>(nm) /
                                                      static_cast<double>(length > 1 ? length : 1))));
    }
    const int strand = (bb->core.flag & BAM_FREVERSE) ? 1 : 0;
    if (p->is_del || p->is_refskip) {
        r[0] = 4; r[1] = 0; r[2] = static_cast<uint8_t>(mapq); r[3] = static_cast<uint8_t>(strand);
        r[4] = 0; r[5] = 0; r[6] = static_cast<uint8_t>(flags | 4); r[7] = static_cast<uint8_t>(md);
        return true;
    }
    const int qpos = p->qpos;
    uint8_t* qual = bam_get_qual(bb);
    const int q = (qual[0] == 0xff) ? 0 : static_cast<int>(qual[qpos]);
    int code = 5;
    switch (seq_nt16_str[bam_seqi(bam_get_seq(bb), qpos)]) {
        case 'A': code = 0; break; case 'C': code = 1; break;
        case 'G': code = 2; break; case 'T': code = 3; break; default: break;
    }
    if (code <= 3 && q >= f.min_baseq) flags |= 2;
    r[0] = static_cast<uint8_t>(code); r[1] = static_cast<uint8_t>(q);
    r[2] = static_cast<uint8_t>(mapq); r[3] = static_cast<uint8_t>(strand);
    if (full_rows) {
        int rest = length - 1 - qpos; if (rest < 0) rest = 0;
        int end_dist = qpos < rest ? qpos : rest;
        if (end_dist < 5) flags |= 32;
        double half = static_cast<double>(length) / 2.0; if (half < 1.0) half = 1.0;
        r[4] = static_cast<uint8_t>(clip255(static_cast<long>(std::nearbyint(
            255.0 * static_cast<double>(qpos) / static_cast<double>(length > 1 ? length : 1)))));
        r[5] = static_cast<uint8_t>(clip255(static_cast<long>(std::nearbyint(
            255.0 * static_cast<double>(end_dist) / half))));
    } else { r[4] = r[5] = 0; }
    r[6] = static_cast<uint8_t>(flags); r[7] = static_cast<uint8_t>(md);
    return true;
}

static bool key_less(const KeyEntry& x, const KeyEntry& y) {
    return x.key != y.key ? x.key < y.key : x.idx < y.idx;
}

int select_rows_lazy(const bam_pileup1_t* pil, int n_plp, const Filters& f, bool full_rows, bool lazy_hash,
                     int max_reads, std::vector<KeyEntry>& ke, RowEntry* out, int* n_admitted, const char** err) {
    if (static_cast<size_t>(n_plp) > ke.size()) ke.resize(static_cast<size_t>(n_plp) * 2);
    int n = 0;
    for (int i = 0; i < n_plp; i++) {
        const bam_pileup1_t* p = pil + i;
        bam1_t* bb = p->b;
        if (bb->core.qual < f.min_mapq) continue;
        const uint16_t flag = bb->core.flag;
        if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if (bb->core.l_qseq <= 0) { *err = "l_qseq==0 (V1 cannot reproduce)"; return -1; }
        uint64_t key;
        if (lazy_hash) {
            const char* qn = bam_get_qname(bb);
            key = blake2b_key64(reinterpret_cast<const uint8_t*>(qn), std::strlen(qn));
        } else key = static_cast<uint64_t>(p->cd.i);
        ke[n] = {key, static_cast<uint32_t>(n), static_cast<uint32_t>(i)};
        n++;
    }
    *n_admitted = n;
    if (n == 0) return 0;
    const int keep = n < max_reads ? n : max_reads;
    if (n > keep) std::nth_element(ke.begin(), ke.begin() + keep, ke.begin() + n, key_less);
    std::sort(ke.begin(), ke.begin() + keep, key_less);
    for (int j = 0; j < keep; j++) {
        out[j].key = ke[j].key; out[j].idx = ke[j].idx;
        if (!make_row(pil + ke[j].slot, f, full_rows, out[j].row, err)) return -1;
    }
    if (full_rows && n > keep) {   // V1 refuses non-integer NM on ANY admitted read of the column
        for (int j = keep; j < n; j++) {
            uint8_t tmp[8];
            if (!make_row(pil + ke[j].slot, f, true, tmp, err)) return -1;
        }
    }
    return keep;
}

static bool entry_less(const RowEntry& x, const RowEntry& y) {
    return x.key != y.key ? x.key < y.key : x.idx < y.idx;
}

int select_kept(std::vector<RowEntry>& ents, int n, int max_reads) {
    int keep = n < max_reads ? n : max_reads;
    // V1: stable qsort of ALL admitted reads, keep the smallest `max_reads`.  (key, idx) is a strict
    // total order, so selecting the smallest `keep` and sorting only those yields the same sequence.
    if (n > keep) std::nth_element(ents.begin(), ents.begin() + keep, ents.begin() + n, entry_less);
    std::sort(ents.begin(), ents.begin() + keep, entry_less);
    return keep;
}

const char* window_pass(BamHandle& b, int64_t wstart, int64_t wend, const Filters& f, int max_reads,
                        bool full_rows, uint64_t want, std::vector<RowEntry>& ents,
                        EmitFn emit, void* ctx) {
    hts_itr_t* iter = sam_itr_queryi(b.idx, b.tid, static_cast<hts_pos_t>(wstart), static_cast<hts_pos_t>(wend));
    if (!iter) return "cannot create iterator";
    IterState st{b.fp, iter};
    void* data_arr[1] = {&st};
    bam_mplp_t mplp = bam_mplp_init(1, read_cb, data_arr);
    bam_mplp_set_maxcnt(mplp, 8000);
    bam_mplp_init_overlaps(mplp);
    bam_mplp_constructor(mplp, plp_construct);
    int tid_out, pos, n_plp_arr[1];
    const bam_pileup1_t* pil_arr[1];
    const char* err = nullptr;
    while (!err && bam_mplp_auto(mplp, &tid_out, &pos, n_plp_arr, pil_arr) > 0) {
        if (pos < wstart || pos >= wend) continue;
        if (!((want >> (pos - wstart)) & 1ULL)) continue;
        int n = build_rows(pil_arr[0], n_plp_arr[0], f, full_rows, ents, &err);
        if (err || n <= 0) continue;
        int keep = select_kept(ents, n, max_reads);
        emit(ctx, pos, ents.data(), keep, n);
    }
    bam_mplp_destroy(mplp);
    hts_itr_destroy(iter);
    return err;
}

}  // namespace dnav2
