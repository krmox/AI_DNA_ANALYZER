// V2.x cheap evidence layer (Mechanism A), local-context layer (Mechanism D) and PB@m helper (Mechanism B).
// Everything is computed from data the existing single streaming pileup already touches; nothing here
// changes Binomial / router / PB / Method C semantics.  Candidate loci only (k >= K_MIN): the cost is
// proportional to the (tiny) candidate set, not to the 56 M scored loci.
#pragma once
#include <cstdint>
#include <vector>
#include "extract.hpp"

namespace dnav2 {

constexpr int NFEAT = 48;
extern const char* const FEAT_NAMES[NFEAT];

// Raw per-read accumulators of one candidate column (integers only: deterministic, order independent).
struct ReadStats {
    int n_alt = 0, n_ref = 0, n_other = 0;
    int alt_fwd = 0, alt_rev = 0, ref_fwd = 0, ref_rev = 0;
    int64_t alt_bq = 0, ref_bq = 0, alt_mq = 0, ref_mq = 0;
    int alt_bq_hi = 0, alt_lowbq = 0, alt_mq_lt40 = 0;
    int64_t alt_end = 0, ref_end = 0;
    int alt_near_end = 0;
    int n_ref_samp = 0;
    int64_t alt_nm = 0, ref_nm = 0;
    int alt_clip = 0, ref_clip = 0;
    int alt_uniq = 0;
    int del_next = 0;
    int n_admitted = 0, n_plp = 0;
};

// One pass over the pileup column for the candidate ALT base (0..3) and REF base (0..3, -1 if N).
// `ref_sample` bounds the number of REF reads inspected for the (aux/cigar) features that cost more per read.
void read_stats(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, int ref, int ref_sample, ReadStats& s);

// PB LLR (frozen V1 PB semantics, alt_index recount) on the first `m` admitted reads in V1 hash order
// (m<=0 -> all admitted reads).  Bit-identical to pb_native_v2 when m == 48.
constexpr int NPB = 6;
constexpr int NBLK = 4;   // block-consensus outputs: median, min, max of PB over disjoint 48-read blocks, number of blocks
constexpr int PB_CAPS[NPB] = {24, 48, 64, 96, 128, 0};
struct PbScratch { std::vector<KeyEntry> ke; std::vector<double> dp; };
// Compact copy of a column's PB input: counted reads in V1 hash order + prefix counts (cnt_at[m] = counted reads among
// the first m admitted).  ~2 bytes/read; kept only for the few loci whose PB might be needed, evaluated lazily.
struct PbStore { std::vector<uint8_t> code, bq; std::vector<int> cnt_at; int n_adm = 0; int alt = 0; };
enum PbVariant { PBV_CAP24, PBV_CAP48, PBV_CAP64, PBV_CAP96, PBV_CAP128, PBV_ALL, PBV_BLK_MED, PBV_BLK_MIN, PBV_BLK_MAX, PBV_N };
bool pb_store_build(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, bool lazy_hash, PbScratch& sc, PbStore& st);
double pb_store_eval(const PbStore& st, int variant, std::vector<double>& dp);
// returns false on error (l_qseq==0); out[NPB]; *n_adm = admitted reads
bool pb_variants(const bam_pileup1_t* pil, int n_plp, const Filters& f, int alt, bool lazy_hash, PbScratch& sc,
                 double out[NPB], double blk[NBLK], int* n_adm);

// PB LLR for a compact list of counted reads (any length)
double pb_llr_n(const uint8_t* code, const uint8_t* bq, int n_counted, int alt_index, std::vector<double>& dp);

// Derived features: fills f[NFEAT] from the raw stats; context fields are filled separately by the engine
// (they need neighbouring columns) and passed in `ctx`.
struct CtxStats {
    int nc2 = 0, nc5 = 0, nc10 = 0, nc25 = 0, nearest = 26;
    double lmr25 = 0, lgap10 = 0, depth_ratio25 = 1, neigh_max_vaf25 = 0;
    int hp_len = 1, left_is_alt = 0, right_is_alt = 0, dinuc10 = 0, is_ts = 0, ref_cg = 0, cpg = 0;
};
void derive_features(const ReadStats& s, int depth, int n, int k, int nref, int k2, double bllr,
                     const CtxStats& c, int gap, int ins, double f[NFEAT]);

}  // namespace dnav2
