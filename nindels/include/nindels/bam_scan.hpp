// Independent candidate discovery from a BAM (htslib), and construction of
// the per-candidate read evidence used by the engine.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "nindels/engine.hpp"
#include "nindels/str_model.hpp"

namespace nindels {

struct ScanParams {
    int min_alt = 2;              // ARBITRARY: reads with the same normalized indel to keep a candidate
    double min_alt_frac = 0.03;   // ARBITRARY: minimum f to keep a candidate
    int64_t chunk = 250000;       // bases per work unit
    int64_t read_pad = 1000;      // reads fetched beyond the chunk
    int span_pad = 3;             // ARBITRARY: a Ref read must cover the event (and its repeat) +-3 bp
    int mismatch_flank = 20;      // ARBITRARY: window for the alt-read mismatch rate
    int bq_flank = 5;             // ARBITRARY: window for the flanking base quality
    // Experiment C: fragment-level evidence. Collapse unmarked duplicates
    // (same start, strand, mate start, insert size, read1/2) and keep one mate
    // per fragment at each candidate (the one with higher flanking BQ).
    bool fragments = false;
    // Which mate represents the fragment: false = higher flanking BQ (C),
    // true = event farther from the read end, i.e. the more reliable alignment (C2).
    bool fragment_pick_by_end = false;
};

// Loads reference [start, end) of a contig with faidx.
RefSeq load_reference(const std::string& fasta, const std::string& chrom, int64_t start, int64_t end);

class BamScanner {
public:
    BamScanner(const std::string& bam_path, const std::string& ref_fasta_for_cram);
    ~BamScanner();
    BamScanner(const BamScanner&) = delete;
    BamScanner& operator=(const BamScanner&) = delete;

    int tid(const std::string& chrom) const;
    std::vector<std::pair<std::string, int64_t>> contigs() const;
    std::string sample_name() const;  // first @RG SM, or "SAMPLE"

    // Candidates whose normalized position lies in [start, end). If external is
    // given (shadow mode), those indels are evaluated instead of discovered ones.
    std::vector<CandidateInput> scan(const RefSeq& ref, int32_t tid, int64_t start, int64_t end,
                                     const ScanParams& p, const std::vector<Indel>* external = nullptr,
                                     int window_half = 10);

    // Experiment B estimation: fragments spanning run +-3 bp (MAPQ >= min_mapq,
    // one read per fragment) and how many carry a whole-unit indel inside the run.
    StrLocusCount count_str_locus(const StrRun& run, int32_t tid, int min_mapq = 20);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace nindels
