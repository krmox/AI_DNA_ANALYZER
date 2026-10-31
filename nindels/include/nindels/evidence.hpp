// Spec §2.1, §3, §4: evidence vector E_i, support fraction f_i and L_depth.
// Every feature computed here is passed to the final model (no dead features).
#pragma once

#include <vector>

#include "nindels/types.hpp"

namespace nindels {

enum class ReadAllele : uint8_t { Ref, Alt, Other };

// What one informative read says about one candidate.
struct ReadObservation {
    ReadAllele allele = ReadAllele::Ref;
    bool reverse = false;
    uint8_t mapq = 0;
    float rel_pos = 0;        // min(dist to read start, dist to read end) / read length, 0..0.5
    int soft_clip = 0;        // total soft-clipped bases of the read
    float flank_bq = 0;       // mean base quality around the event in this read
    float mismatch_rate = 0;  // mismatches / aligned bases near the event (event excluded)
    bool split = false;       // read has a supplementary alignment (SA tag): mapping conflict
};

struct EvidenceParams {
    int clip_threshold = 10;  // ARBITRARY: soft clip >= 10 bp counts as "clipped"
    int low_mapq = 20;        // ARBITRARY: MAPQ < 20 counts as low
};

struct Evidence {
    int n_alt = 0, n_ref = 0, n_other = 0;
    int n_total = 0;      // informative reads, all alleles (local depth)
    Opt f;                // n_alt / (n_alt + n_ref); missing if the denominator is 0
    double l_depth = 0;   // ln(1 + n_alt + n_ref)
    // b_i: strand
    Opt alt_fwd_frac;     // alt reads on the forward strand / n_alt
    Opt strand_bias;      // |log odds ratio| of alt/ref x fwd/rev, +0.5 pseudo-counts
    // read position
    Opt readpos_alt;      // mean rel_pos of alt reads
    Opt readpos_diff;     // readpos_alt - mean rel_pos of ref reads
    // a_i: alignment
    Opt clip_alt;         // fraction of alt reads with soft clip >= threshold
    Opt clip_all;         // same over all informative reads
    Opt mismatch_alt;     // mean mismatch rate of alt reads near the event
    Opt other_frac;       // n_other / n_total: competing indel alleles at the locus
    Opt split_frac;       // fraction of informative reads with a supplementary alignment
    // q_i: quality
    Opt mapq_alt;         // mean MAPQ of alt reads
    Opt mq0_frac;         // fraction of informative reads with MAPQ 0
    Opt mqlow_frac;       // fraction with MAPQ < low_mapq
    Opt bq_alt;           // mean flanking base quality of alt reads
    // G_r input for region complexity (spec §10)
    Opt instability;      // fraction of informative reads with low MAPQ or clip >= threshold
};

Evidence summarize_evidence(const std::vector<ReadObservation>& obs, const EvidenceParams& p = {});

}  // namespace nindels
