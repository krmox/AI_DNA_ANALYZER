// Evaluation frame of the frozen V1 protocol: which 64-bp windows exist (BED containment + N
// filter, providers.GiabAlignmentProvider._build_windows) and which loci are scored (truth label
// in {0 = normal, 1 = SNP}, providers._fetch_labels).  The truth VCF is used ONLY to reproduce
// V1's scored-locus frame; it never influences a call, a score or the router.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace dnav2 {

constexpr int LABEL_NORMAL = 0, LABEL_SNP = 1, LABEL_INSERTION = 2, LABEL_DELETION = 3;
constexpr int SEQ_LEN = 64;   // V1 window width == stride

struct Window { int64_t start; };  // [start, start+64)

struct Frame {
    std::string contig;
    int64_t contig_len = 0;
    std::vector<Window> windows;        // ordered
    std::vector<uint8_t> labels;        // per contig position, truth-derived (frame filter only)
    // reference index per position of [start,end): 0=N/other, 1..4 = A,C,G,T (uppercased)
};

// Load BED intervals for `contig` (merged like providers.load_bed_regions).
std::vector<std::pair<int64_t, int64_t>> load_bed(const std::string& path, const std::string& contig);

// Build the window list exactly as V1 (region (0, contig_len), stride 64).
// Errors are returned as a non-empty string.
std::string build_windows(const std::string& fasta, const std::string& bed, const std::string& contig,
                          int64_t contig_len, std::vector<Window>& out);

// Truth-derived per-position label array (LABEL_*), V1 _fetch_labels semantics applied globally.
std::string load_labels(const std::string& truth_vcf, const std::string& contig, int64_t contig_len,
                        std::vector<uint8_t>& labels);

// Reference index for [start,end) (0=N/other, 1..4=ACGT).  Not thread-safe per Fai handle.
struct Fai;
Fai* fai_open(const std::string& fasta);
void fai_close(Fai*);
std::string fetch_ref_idx(Fai*, const std::string& contig, int64_t start, int64_t end, uint8_t* out,
                          int* n_count = nullptr);

}  // namespace dnav2
