// Floor measurement: BGZF inflate + BAM record decode of a region, nothing else.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "htslib/sam.h"
int main(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "bamscan BAM contig start end [decomp_threads]\n"); return 2; }
    samFile* fp = sam_open(argv[1], "r");
    if (argc > 5) hts_set_threads(fp, std::atoi(argv[5]));
    sam_hdr_t* h = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, argv[1]);
    int tid = sam_hdr_name2tid(h, argv[2]);
    auto t0 = std::chrono::steady_clock::now();
    hts_itr_t* it = sam_itr_queryi(idx, tid, atoll(argv[3]), atoll(argv[4]));
    bam1_t* b = bam_init1();
    long n = 0, bases = 0;
    while (sam_itr_next(fp, it, b) >= 0) { n++; bases += b->core.l_qseq; }
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("reads=%ld bases=%ld seconds=%.3f reads_per_s=%.0f\n", n, bases, s, n / s);
    return 0;
}
