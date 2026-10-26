// Per-sample self-calibration pre-pass (R2). Deterministic: evenly spaced blocks of consecutive windows, one thread.
#pragma once
#include <string>
#include <utility>
#include <vector>
#include "dnav2/frame.hpp"
#include "tail.hpp"

namespace dnav2 {

// Samples `blocks` runs of `block_windows` consecutive 64-bp windows spread evenly over `windows` (only windows inside the BED clip
// intervals count) and fills `out`. Returns "" on success; on failure out.valid stays false and the reason is returned.
std::string calibrate(const std::string& bam, const std::string& fasta, const std::string& contig,
                      const std::vector<Window>& windows, const std::vector<std::pair<int64_t, int64_t>>* bed_clip,
                      int blocks, int block_windows, int threads, Calibration& out);

}  // namespace dnav2
