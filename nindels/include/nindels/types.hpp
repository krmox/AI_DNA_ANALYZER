// Core value types shared by every stage of the Indel Pattern Engine.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace nindels {

// A feature value that may be unavailable (spec §3, §12: missing values are
// never silently replaced by zero).
using Opt = std::optional<double>;

enum class IndelType : uint8_t { Insertion = 0, Deletion = 1 };

// A left-normalized indel. pos is 0-based:
//   deletion  -> first deleted reference base,
//   insertion -> reference base the sequence is inserted in front of.
// seq is the inserted or deleted bases (m_i in spec §2).
struct Indel {
    int32_t tid = -1;
    int64_t pos = 0;
    IndelType type = IndelType::Deletion;
    std::string seq;

    int length() const { return static_cast<int>(seq.size()); }
    // Reference bases consumed by the event (0 for insertions).
    int ref_span() const { return type == IndelType::Deletion ? length() : 0; }
    bool operator==(const Indel&) const = default;
};

struct IndelHash {
    size_t operator()(const Indel& x) const noexcept;
};

bool operator<(const Indel& a, const Indel& b);

// Reference sequence (upper case) covering [offset, offset + seq.size()).
struct RefSeq {
    std::string seq;
    int64_t offset = 0;

    char at(int64_t pos) const {
        int64_t i = pos - offset;
        return (i >= 0 && i < static_cast<int64_t>(seq.size())) ? seq[i] : 'N';
    }
    bool contains(int64_t pos) const {
        return pos >= offset && pos < offset + static_cast<int64_t>(seq.size());
    }
    // Bases [pos, pos+len) clipped to the stored range.
    std::string sub(int64_t pos, int64_t len) const;
    int64_t end() const { return offset + static_cast<int64_t>(seq.size()); }
};

// Deterministic 64-bit mixing (splitmix64), used for every pseudo-random choice
// so that runs are reproducible (spec §16 rule 14).
inline uint64_t mix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
uint64_t hash_string(std::string_view s, uint64_t seed = 0);

}  // namespace nindels
