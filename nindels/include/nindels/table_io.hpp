// Tab-separated candidate tables: one row per candidate with its pattern, its
// truth label (if known) and every feature. Used to build the memory, fit the
// models and inspect results.
#pragma once

#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

#include "nindels/engine.hpp"

namespace nindels {

std::vector<std::string> split_tab(std::string_view line);
Opt parse_opt(std::string_view s);  // "NA" -> missing
std::string format_opt(Opt v);

// VCF REF/ALT with the anchor base before the event, and the 1-based POS.
struct VcfAlleles {
    int64_t pos1;
    std::string ref, alt;
};
VcfAlleles vcf_alleles(const RefSeq& ref, const Indel& x);

// Label codes: 1 true (in truth set), 0 false (inside the confident regions,
// not in the truth set), -1 unlabelled (outside the confident regions or no truth given).
void write_table_header(std::ostream& o, const std::vector<std::string>& comments);
void write_table_row(std::ostream& o, const std::string& chrom, const RefSeq& ref, const CandidateResult& r,
                     int label);

struct Table {
    std::vector<std::string> comments;  // lines starting with "##"
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;
    int col(const std::string& name) const;  // throws if absent
};
Table read_table(const std::string& path);

Pattern pattern_from_row(const Table& t, const std::vector<std::string>& row);

}  // namespace nindels
