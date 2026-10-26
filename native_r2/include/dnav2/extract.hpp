// BAM extraction for the V2 engine.  Two passes, each replicating one V1 pass exactly:
//   count_pass  == native/pileup_native.c  (one streaming pileup over a run of windows)
//   window_pass == native/reads_native.c   (one FRESH pileup per 64-bp window; V1 read-tensor
//                                           semantics incl. per-window mate-overlap handling)
// but without materialising any tensor: window_pass hands each locus' <=48 retained reads to a
// callback and forgets them.
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "blake2b.hpp"

namespace dnav2 {

#ifdef DNAV2_PROFILE
struct ProfCounters {   // nanoseconds / counts, summed over threads at thread exit
    uint64_t t_mplp = 0, t_build = 0, t_select = 0, t_emit = 0, t_hash = 0, n_hash = 0, n_loci_emitted = 0,
             n_reads_admitted = 0, n_reads_seen = 0, t_count_mplp = 0, t_count_body = 0, n_count_loci = 0, n_count_reads = 0;
};
ProfCounters& prof();   // thread-local accumulator
void prof_dump(const char* path);
#endif

struct Filters { int min_mapq = 20; int min_baseq = 13; };

struct CountRow { uint32_t v[9]; };  // A C G T gap ins depth qsum msum  (V1 columns 0..8)

struct BamHandle {
    samFile* fp = nullptr; sam_hdr_t* hdr = nullptr; hts_idx_t* idx = nullptr; int tid = -1;
    std::string open(const std::string& path, const std::string& contig, int n_decomp_threads = 0);
    ~BamHandle();
};

// Streaming pileup over [start,end): fills out[end-start] (must be zero-initialised).
// Returns nullptr on success or an error string.
const char* count_pass(BamHandle& b, int64_t start, int64_t end, const Filters& f, CountRow* out);

struct RowEntry { uint64_t key; uint32_t idx; uint8_t row[8]; };

// One V1 window pass over [wstart,wend).  For every locus with >=1 admitted read, calls
//   emit(pos, kept, n_kept, n_admitted)
// where kept = the min(48, n_admitted) retained reads sorted by (blake2b key, arrival index),
// exactly V1's tensor order.  full_rows=false fills only base_code, base_quality and the flags
// byte (all PB needs); true reproduces all 8 V1 columns.  `want` (64-bit mask over window
// offsets, bit i = locus wstart+i) restricts which loci are emitted; the pileup itself always
// runs over the whole window (that is what V1 does).
using EmitFn = void (*)(void* ctx, int64_t pos, const RowEntry* kept, int n_kept, int n_admitted);
const char* window_pass(BamHandle& b, int64_t wstart, int64_t wend, const Filters& f, int max_reads,
                        bool full_rows, uint64_t want, std::vector<RowEntry>& scratch,
                        EmitFn emit, void* ctx);

// ---- column-level helpers shared by both paths ------------------------------------------------
// Counts of one pileup column, replicating pileup_native.c's inner loop.
CountRow count_column(const bam_pileup1_t* pil, int n_plp, const Filters& f);

// Admitted reads of one column as V1 tensor rows (unsorted, arrival order in ents[0..n)).  Returns n,
// or -1 with *err set when V1 could not reproduce the row (l_qseq==0 / non-integer NM).
int build_rows(const bam_pileup1_t* pil, int n_plp, const Filters& f, bool full_rows,
               std::vector<RowEntry>& ents, const char** err);

// V1 ordering: keep the `max_reads` smallest (blake2b key, arrival idx); leaves them sorted in ents[0..keep).
int select_kept(std::vector<RowEntry>& ents, int n, int max_reads);

// Lazy variant (same result as build_rows + select_kept, less work): admitted reads are first reduced to
// (key, arrival idx, pileup slot); only the `max_reads` selected reads get their row bytes built.
// `lazy_hash`: the pileup was run WITHOUT the read-construct hash hook, compute BLAKE2b from the qname here
// (only worthwhile when few columns need rows, e.g. cascade).  Fills out[0..keep) and returns keep;
// *n_admitted = admitted reads; returns -1 + *err as build_rows.
struct KeyEntry { uint64_t key; uint32_t idx; uint32_t slot; };
int select_rows_lazy(const bam_pileup1_t* pil, int n_plp, const Filters& f, bool full_rows, bool lazy_hash,
                     int max_reads, std::vector<KeyEntry>& ke, RowEntry* out, int* n_admitted, const char** err);

// Fused streaming pass: ONE mplp over [start,end) (read-name hash attached at read construction);
// fn(ctx, pos, pil, n_plp) is called for every pileup column inside [start,end).
using ColumnFn = void (*)(void* ctx, int64_t pos, const bam_pileup1_t* pil, int n_plp);
const char* pileup_pass(BamHandle& b, int64_t start, int64_t end, ColumnFn fn, void* ctx, bool hash_at_construct = true);

}  // namespace dnav2
