// Spec §2 and §13 steps 1-2: normalization and the pattern vector
// P_i = (t_i, L_i, m_i, c_i, h_i, r_i).
#pragma once

#include "nindels/types.hpp"

namespace nindels {

// Shift an indel to its leftmost equivalent position in the reference.
Indel left_normalize(const RefSeq& ref, Indel x);

struct Pattern {
    IndelType type = IndelType::Deletion;  // t_i
    int length = 0;                        // L_i
    std::string motif;                     // m_i
    std::string context;                   // c_i: flank left + flank right of the event
    int hp_length = 0;                     // h_i: longest homopolymer run touching the event
    char hp_base = 'N';
    int period = 0;       // r_i: smallest repeat unit (1..6) that tiles m_i, 0 if none
    int copies = 0;       // r_i: consecutive copies of that unit in the reference at the site
    int repeat_span = 0;  // r_i: period * copies
    double entropy = 0;   // r_i: Shannon entropy of c_i in bits (0..2), local complexity
};

struct PatternParams {
    int context_flank = 20;  // ASSUMPTION (arbitrary): bases of reference on each side
    int max_period = 6;      // ASSUMPTION (arbitrary): longest tandem unit considered
};

Pattern extract_pattern(const RefSeq& ref, const Indel& x, const PatternParams& p = {});

// Shannon entropy (bits) of the A/C/G/T composition of s; N ignored.
double shannon_entropy(std::string_view s);

}  // namespace nindels
