#include "dnav2/features.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include "dnav2/pb.hpp"

namespace dnav2 {

const char* const FEAT_NAMES[NFEAT] = {
    "depth", "n", "k", "nref", "k2", "vaf", "bllr",
    "alt_bq_mean", "ref_bq_mean", "alt_bq_hi_frac", "alt_lowbq", "alt_mq_mean", "ref_mq_mean", "alt_mq_lt40_frac",
    "alt_fwd_frac", "ref_fwd_frac", "strand_diff", "alt_min_strand",
    "alt_end_mean", "ref_end_mean", "alt_near_end_frac",
    "alt_nm_mean", "ref_nm_mean", "alt_clip_frac", "ref_clip_frac", "alt_uniq_frac", "gap_frac", "ins_frac", "del_next_frac",
    "nc2", "nc5", "nc10", "nc25", "nearest", "lmr25", "lgap10", "depth_ratio25", "neigh_max_vaf25",
    "hp_len", "left_is_alt", "right_is_alt", "dinuc10", "is_ts", "ref_cg", "cpg",
    "bq_diff", "mq_diff", "lowmq_frac"};

namespace {
inline int base_code(const bam1_t* bb, int qpos) {
    switch (seq_nt16_str[bam_seqi(bam_get_seq(bb), qpos)]) {
        case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1;
    }
}
inline bool clipped(const bam1_t* bb) {
    const uint32_t* cg = bam_get_cigar(bb);
    const int nc = bb->core.n_cigar;
    if (nc == 0) return false;
    return bam_cigar_op(cg[0]) == BAM_CSOFT_CLIP || bam_cigar_op(cg[nc - 1]) == BAM_CSOFT_CLIP;
}
inline int nm_of(const bam1_t* bb) {
    uint8_t* aux = bam_aux_get(const_cast<bam1_t*>(bb), "NM");
    if (!aux) return 0;
    char t = static_cast<char>(*aux);
    if (t == 'c' || t == 'C' || t == 's' || t == 'S' || t == 'i' || t == 'I') return static_cast<int>(bam_aux2i(aux));
    return 0;
}
inline bool key_less(const KeyEntry& x, const KeyEntry& y) { return x.key != y.key ? x.key < y.key : x.idx < y.idx; }
}  // namespace

void read_stats(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, int ref, int ref_sample, ReadStats& s) {
    s = ReadStats();
    s.n_plp = n_plp;
    std::vector<int64_t> starts;
    for (int i = 0; i < n_plp; i++) {
        const bam_pileup1_t* p = pil + i;
        const bam1_t* bb = p->b;
        const int mapq = bb->core.qual;
        if (mapq < f.min_mapq) continue;
        const uint16_t flag = bb->core.flag;
        if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        s.n_admitted++;
        if (p->indel < 0) s.del_next++;
        if (p->is_del || p->is_refskip) continue;
        const int qpos = p->qpos;
        const uint8_t* qual = bam_get_qual(bb);
        const int bq = (qual[0] == 0xff) ? 0 : static_cast<int>(qual[qpos]);
        const int code = base_code(bb, qpos);
        if (code < 0) continue;
        if (bq < f.min_baseq) { if (code == alt) s.alt_lowbq++; continue; }
        const bool rev = (flag & BAM_FREVERSE) != 0;
        const int rest = bb->core.l_qseq - 1 - qpos;
        const int edist = qpos < rest ? qpos : (rest < 0 ? 0 : rest);
        if (code == alt) {
            s.n_alt++;
            (rev ? s.alt_rev : s.alt_fwd)++;
            s.alt_bq += bq; s.alt_mq += mapq;
            s.alt_bq_hi += bq >= 20; s.alt_mq_lt40 += mapq < 40;
            s.alt_end += edist; s.alt_near_end += edist < 10;
            s.alt_nm += nm_of(bb); s.alt_clip += clipped(bb);
            starts.push_back((static_cast<int64_t>(bb->core.pos) << 1) | (rev ? 1 : 0));
        } else if (code == ref) {
            s.n_ref++;
            (rev ? s.ref_rev : s.ref_fwd)++;
            s.ref_bq += bq; s.ref_mq += mapq;
            if (s.n_ref_samp < ref_sample) {   // pricier per-read features on a bounded sample only
                s.n_ref_samp++;
                s.ref_end += edist; s.ref_nm += nm_of(bb); s.ref_clip += clipped(bb);
            }
        } else s.n_other++;
    }
    std::sort(starts.begin(), starts.end());
    s.alt_uniq = static_cast<int>(std::unique(starts.begin(), starts.end()) - starts.begin());
}

double pb_llr_n(const uint8_t* code, const uint8_t* bq, int n_counted, int alt_index, std::vector<double>& dp) {
    const PbTables& T = pb_tables();
    int k = 0;
    for (int i = 0; i < n_counted; i++) k += (code[i] == alt_index);
    if (dp.size() < static_cast<size_t>(k) + 1) dp.resize(static_cast<size_t>(k) + 1);
    double res[3];
    for (int h = 0; h < 3; h++) {
        dp[0] = 0.0;
        for (int j = 1; j <= k; j++) dp[j] = -std::numeric_limits<double>::infinity();
        const double* lq = T.lq[h];
        const double* lnq = T.lnq[h];
        for (int i = 0; i < n_counted; i++) {
            const double a = lq[bq[i]], b = lnq[bq[i]];
            int top = (i + 1 < k) ? i + 1 : k;
            for (int j = top; j >= 1; j--) dp[j] = logaddexp(dp[j] + b, dp[j - 1] + a);
            dp[0] += b;
        }
        res[h] = dp[k];
    }
    return pb_llr_from(res[0], res[1], res[2], n_counted);
}

bool pb_store_build(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, bool lazy_hash, PbScratch& sc, PbStore& st) {
    auto& ke = sc.ke;
    ke.clear();
    for (int i = 0; i < n_plp; i++) {
        const bam_pileup1_t* p = pil + i;
        const bam1_t* bb = p->b;
        if (bb->core.qual < f.min_mapq) continue;
        const uint16_t flag = bb->core.flag;
        if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if (bb->core.l_qseq <= 0) return false;
        uint64_t key;
        if (lazy_hash) {
            const char* qn = bam_get_qname(bb);
            key = blake2b_key64(reinterpret_cast<const uint8_t*>(qn), std::strlen(qn));
        } else key = static_cast<uint64_t>(p->cd.i);
        ke.push_back({key, static_cast<uint32_t>(ke.size()), static_cast<uint32_t>(i)});
    }
    const int n = static_cast<int>(ke.size());
    st.n_adm = n; st.alt = alt;
    st.code.clear(); st.bq.clear(); st.cnt_at.assign(n + 1, 0);
    if (n == 0) return true;
    std::sort(ke.begin(), ke.end(), key_less);
    st.code.reserve(n); st.bq.reserve(n);
    for (int j = 0; j < n; j++) {
        const bam_pileup1_t* p = pil + ke[j].slot;
        if (!(p->is_del || p->is_refskip)) {
            const bam1_t* bb = p->b;
            const uint8_t* qual = bam_get_qual(bb);
            const int q = (qual[0] == 0xff) ? 0 : static_cast<int>(qual[p->qpos]);
            const int code = base_code(bb, p->qpos);
            if (code >= 0 && q >= f.min_baseq) { st.code.push_back(static_cast<uint8_t>(code)); st.bq.push_back(static_cast<uint8_t>(q)); }
        }
        st.cnt_at[j + 1] = static_cast<int>(st.code.size());
    }
    return true;
}

double pb_store_eval(const PbStore& st, int variant, std::vector<double>& dp) {
    const int n = st.n_adm;
    if (n == 0) return 0.0;
    if (variant <= PBV_ALL) {
        const int cap = PB_CAPS[variant];
        const int m = cap <= 0 || cap > n ? n : cap;
        return pb_llr_n(st.code.data(), st.bq.data(), st.cnt_at[m], st.alt, dp);
    }
    // block consensus: PB (frozen semantics) on each disjoint block of 48 admitted reads in hash order
    std::vector<double> bl;
    for (int b0 = 0; b0 + 48 <= n || b0 == 0; b0 += 48) {
        const int b1 = b0 + 48 < n ? b0 + 48 : n;
        bl.push_back(pb_llr_n(st.code.data() + st.cnt_at[b0], st.bq.data() + st.cnt_at[b0], st.cnt_at[b1] - st.cnt_at[b0], st.alt, dp));
        if (b1 == n) break;
    }
    std::sort(bl.begin(), bl.end());
    const size_t nb = bl.size();
    if (variant == PBV_BLK_MIN) return bl.front();
    if (variant == PBV_BLK_MAX) return bl.back();
    return (nb % 2) ? bl[nb / 2] : 0.5 * (bl[nb / 2 - 1] + bl[nb / 2]);
}

bool pb_variants(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, bool lazy_hash, PbScratch& sc,
                 double out[NPB], double blk[NBLK], int* n_adm) {
    PbStore st;
    if (!pb_store_build(pil, n_plp, f, alt, lazy_hash, sc, st)) return false;
    *n_adm = st.n_adm;
    for (int c = 0; c < NPB; c++) out[c] = pb_store_eval(st, c, sc.dp);
    blk[0] = pb_store_eval(st, PBV_BLK_MED, sc.dp);
    blk[1] = pb_store_eval(st, PBV_BLK_MIN, sc.dp);
    blk[2] = pb_store_eval(st, PBV_BLK_MAX, sc.dp);
    int nb = 0;
    for (int b0 = 0; b0 + 48 <= st.n_adm || b0 == 0; b0 += 48) { nb++; if (b0 + 48 >= st.n_adm) break; }
    blk[3] = st.n_adm == 0 ? 0.0 : nb;
    return true;
}

void derive_features(const ReadStats& s, int depth, int n, int k, int nref, int k2, double bllr, const CtxStats& c,
                     int gap, int ins, double f[NFEAT]) {
    auto div = [](double a, double b) { return b > 0 ? a / b : 0.0; };
    const double ka = s.n_alt, kr = s.n_ref;
    f[0] = depth; f[1] = n; f[2] = k; f[3] = nref; f[4] = k2; f[5] = div(k, n); f[6] = bllr;
    f[7] = div(s.alt_bq, ka); f[8] = div(s.ref_bq, kr);
    f[9] = div(s.alt_bq_hi, ka); f[10] = s.alt_lowbq;
    f[11] = div(s.alt_mq, ka); f[12] = div(s.ref_mq, kr); f[13] = div(s.alt_mq_lt40, ka);
    f[14] = div(s.alt_fwd, ka); f[15] = div(s.ref_fwd, kr);
    f[16] = std::fabs(f[14] - f[15]); f[17] = std::min(s.alt_fwd, s.alt_rev);
    f[18] = div(s.alt_end, ka); f[19] = div(s.ref_end, s.n_ref_samp); f[20] = div(s.alt_near_end, ka);
    f[21] = div(s.alt_nm, ka); f[22] = div(s.ref_nm, s.n_ref_samp);
    f[23] = div(s.alt_clip, ka); f[24] = div(s.ref_clip, s.n_ref_samp);
    f[25] = div(s.alt_uniq, ka);
    f[26] = div(gap, depth); f[27] = div(ins, depth); f[28] = div(s.del_next, depth);
    f[29] = c.nc2; f[30] = c.nc5; f[31] = c.nc10; f[32] = c.nc25; f[33] = c.nearest;
    f[34] = c.lmr25; f[35] = c.lgap10; f[36] = c.depth_ratio25; f[37] = c.neigh_max_vaf25;
    f[38] = c.hp_len; f[39] = c.left_is_alt; f[40] = c.right_is_alt; f[41] = c.dinuc10;
    f[42] = c.is_ts; f[43] = c.ref_cg; f[44] = c.cpg;
    f[45] = f[8] - f[7]; f[46] = f[12] - f[11];
    f[47] = div(s.n_plp - s.n_admitted, s.n_plp);
}

}  // namespace dnav2
