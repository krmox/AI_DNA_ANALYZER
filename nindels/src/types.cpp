#include "nindels/types.hpp"

#include <algorithm>
#include <tuple>

namespace nindels {

uint64_t hash_string(std::string_view s, uint64_t seed) {
    uint64_t h = mix64(seed ^ 0x51ED270B27A1F3C9ULL);
    for (unsigned char c : s) h = mix64(h ^ c);
    return h;
}

size_t IndelHash::operator()(const Indel& x) const noexcept {
    uint64_t h = mix64(static_cast<uint64_t>(x.tid) * 31 + static_cast<uint64_t>(x.pos));
    h = mix64(h ^ static_cast<uint64_t>(x.type));
    return static_cast<size_t>(hash_string(x.seq, h));
}

bool operator<(const Indel& a, const Indel& b) {
    return std::tie(a.tid, a.pos, a.type, a.seq) < std::tie(b.tid, b.pos, b.type, b.seq);
}

std::string RefSeq::sub(int64_t pos, int64_t len) const {
    int64_t b = std::max(pos, offset);
    int64_t e = std::min(pos + len, end());
    if (e <= b) return {};
    return seq.substr(static_cast<size_t>(b - offset), static_cast<size_t>(e - b));
}

}  // namespace nindels
