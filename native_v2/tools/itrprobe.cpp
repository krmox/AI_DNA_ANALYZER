// How many BAM records does a fresh 64-bp iterator query actually decode?  (tests the 16-kb linear-index hypothesis)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "htslib/sam.h"
int main(int argc, char** argv) {
    samFile* fp = sam_open(argv[1], "r"); hts_set_cache_size(fp, 8 << 20);
    sam_hdr_t* h = sam_hdr_read(fp); hts_idx_t* idx = sam_index_load(fp, argv[1]);
    int tid = sam_hdr_name2tid(h, "chr20"); long start = atol(argv[2]); int nw = atoi(argv[3]);
    bam1_t* b = bam_init1();
    double t_itr = 0; long returned = 0;
    auto T0 = std::chrono::steady_clock::now();
    for (int w = 0; w < nw; w++) {
        auto t0 = std::chrono::steady_clock::now();
        hts_itr_t* it = sam_itr_queryi(idx, tid, start + 64L * w, start + 64L * (w + 1));
        while (sam_itr_next(fp, it, b) >= 0) returned++;
        hts_itr_destroy(it);
        t_itr += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    // linear-index-driven decode volume, measured by re-scanning from the same file offset is not exposed;
    // report timing and returned reads per window instead
    std::printf("windows=%d returned_reads=%ld (%.1f/window) itr_seconds=%.3f (%.1f us/window)\n", nw, returned,
                double(returned) / nw, t_itr, 1e6 * t_itr / nw);
    return 0;
}
