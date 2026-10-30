// Experiment A: joint multi-allelic genotyping of an indel locus.
//
// Candidates whose footprints (event + its repeat) overlap form one locus with
// alleles 1..K and reference allele 0. Every read is scored against every
// haplotype with the PairHMM; genotypes are unordered allele pairs, so a read
// from a third allele is no longer forced to support "the" alternative.
#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "nindels/engine.hpp"

namespace nindels {

struct LocusParams {
    int max_alleles = 5;          // ARBITRARY: alleles kept per locus (by supporting reads)
    double theta = 1e-4;          // ARBITRARY prior: P(0/a) = theta, P(a/a) = theta/2
    double hetalt_factor = 0.1;   // ARBITRARY prior: P(a/b) = theta * hetalt_factor
    double min_qual = 20;         // phred posterior of non-reference (same value as the baseline rule)
};

struct LocusCall {
    std::vector<Indel> alleles;    // index a-1 = allele a
    int g1 = 0, g2 = 0;            // called genotype (allele indices, 0 = reference)
    double qual = 0;               // -10 log10 P(0/0 | reads)
    double gq = 0;                 // -10 log10 (1 - P(called genotype | reads)), capped at 99
    std::vector<int> ad;           // reads whose best haplotype is allele 0..K
    int depth = 0;
    bool pass = false;
};

// Groups of indices into `inputs` (sorted by position) forming one locus each.
std::vector<std::vector<size_t>> cluster_candidates(const RefSeq& ref, const std::vector<const CandidateInput*>& inputs);

LocusCall genotype_locus(const RefSeq& ref, const std::vector<const CandidateInput*>& group,
                         const LikelihoodParams& lp, const LocusParams& p);

void write_locus_vcf(std::ostream& o, const std::string& chrom, const RefSeq& ref, const LocusCall& c);

}  // namespace nindels
