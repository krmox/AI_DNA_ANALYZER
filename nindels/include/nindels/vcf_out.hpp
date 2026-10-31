// Spec §13 step 17: VCF output carrying the decision, the confidence, the
// closest patterns, the pattern statistics and the reassessment information.
#pragma once

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "nindels/engine.hpp"

namespace nindels {

struct VcfHeaderInfo {
    std::string sample = "SAMPLE";
    std::vector<std::pair<std::string, int64_t>> contigs;
    std::vector<std::string> extra;  // extra ## lines (command, parameters, file hashes)
};

void write_vcf_header(std::ostream& o, const VcfHeaderInfo& h);
void write_vcf_record(std::ostream& o, const std::string& chrom, const RefSeq& ref, const CandidateResult& r);

}  // namespace nindels
