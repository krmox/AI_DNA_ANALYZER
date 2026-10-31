#include "nindels/pattern.hpp"

#include <cmath>

namespace nindels {

Indel left_normalize(const RefSeq& ref, Indel x) {
    const int L = x.length();
    if (L == 0) return x;
    if (x.type == IndelType::Deletion) {
        // Deleted block [pos, pos+L) can slide left while ref[pos-1] == ref[pos+L-1].
        while (ref.contains(x.pos - 1)) {
            char a = ref.at(x.pos - 1);
            if (a == 'N' || a != ref.at(x.pos + L - 1)) break;
            --x.pos;
        }
        x.seq = ref.sub(x.pos, L);
    } else {
        // Insertion of s before pos equals insertion of ref[pos-1] + s[0..L-1)
        // before pos-1 whenever ref[pos-1] == s[L-1].
        while (ref.contains(x.pos - 1)) {
            char a = ref.at(x.pos - 1);
            if (a == 'N' || a != x.seq.back()) break;
            x.seq.pop_back();
            x.seq.insert(x.seq.begin(), a);
            --x.pos;
        }
    }
    return x;
}

double shannon_entropy(std::string_view s) {
    double n[4] = {0, 0, 0, 0};
    double tot = 0;
    for (char c : s) {
        int k = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
        if (k >= 0) { n[k] += 1; tot += 1; }
    }
    if (tot == 0) return 0.0;
    double h = 0;
    for (double v : n)
        if (v > 0) h -= (v / tot) * std::log2(v / tot);
    return h;
}

static int run_right(const RefSeq& ref, int64_t pos) {
    char b = ref.at(pos);
    if (b == 'N') return 0;
    int n = 0;
    while (ref.at(pos + n) == b) ++n;
    return n;
}

static int run_left(const RefSeq& ref, int64_t pos) {
    char b = ref.at(pos);
    if (b == 'N') return 0;
    int n = 0;
    while (ref.at(pos - n) == b) ++n;
    return n;
}

Pattern extract_pattern(const RefSeq& ref, const Indel& x, const PatternParams& pp) {
    Pattern p;
    p.type = x.type;
    p.length = x.length();
    p.motif = x.seq;

    const int F = pp.context_flank;
    p.context = ref.sub(x.pos - F, F) + ref.sub(x.pos + x.ref_span(), F);
    p.entropy = shannon_entropy(p.context);

    // h_i: homopolymer run in the reference that starts at the (left-normalized)
    // event or ends immediately before it.
    int r = run_right(ref, x.pos), l = run_left(ref, x.pos - 1);
    if (r >= l) { p.hp_length = r; p.hp_base = ref.at(x.pos); }
    else        { p.hp_length = l; p.hp_base = ref.at(x.pos - 1); }

    // r_i: smallest unit that tiles the motif exactly.
    for (int per = 1; per <= pp.max_period && per <= p.length; ++per) {
        if (p.length % per) continue;
        bool ok = true;
        for (int i = per; i < p.length && ok; ++i) ok = p.motif[i] == p.motif[i % per];
        if (ok) { p.period = per; break; }
    }
    if (p.period > 0) {
        const std::string unit = p.motif.substr(0, p.period);
        int c = 0;
        while (ref.sub(x.pos + static_cast<int64_t>(c) * p.period, p.period) == unit) ++c;
        p.copies = c;
        p.repeat_span = c * p.period;
    }
    return p;
}

}  // namespace nindels
