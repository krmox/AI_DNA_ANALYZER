// Independent validation labels (spec §6: only externally validated labels
// feed the pattern statistics): truth VCF + confident-region BED.
#pragma once

#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "nindels/types.hpp"

namespace nindels {

struct BedRegions {
    std::vector<std::pair<int64_t, int64_t>> iv;  // sorted, 0-based half-open
    // True if the VCF record of x (anchor base .. last deleted base) lies in one interval.
    bool contains(const Indel& x) const;
};

BedRegions load_bed(const std::string& path, const std::string& chrom, int64_t start, int64_t end);

struct TruthSet {
    std::unordered_set<Indel, IndelHash> indels;  // normalized simple indels
    int n_in_bed = 0;     // truth indels inside the BED and region (recall denominator)
    int n_complex = 0;    // non-simple indel alleles that cannot be matched exactly
};

TruthSet load_truth(const std::string& vcf, const std::string& chrom, int64_t start, int64_t end,
                    int32_t tid, const RefSeq& ref, const BedRegions* bed);

// 1 = in truth, 0 = inside the BED but not in truth, -1 = unlabelled.
int truth_label(const TruthSet& t, const BedRegions* bed, const Indel& x);

// Simple-indel conversion of one VCF allele pair; false for SNPs / MNPs / complex alleles.
bool allele_to_indel(int64_t pos0, std::string ref, std::string alt, int32_t tid, Indel& out, bool& complex);

// Indels from any VCF (shadow mode: candidates from another caller).
std::vector<Indel> load_vcf_indels(const std::string& vcf, const std::string& chrom, int64_t start, int64_t end,
                                   int32_t tid, const RefSeq& ref);

}  // namespace nindels
