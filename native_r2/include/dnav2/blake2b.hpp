// BLAKE2b(digest_size=8) read-name key, big-endian: the frozen V1 read ordering key.
// Algorithm text taken from native/reads_native.c (V1), which is unit-tested against
// hashlib.blake2b; V2 re-tests it in tests/test_blake2b.py.
#pragma once
#include <cstdint>
#include <cstring>
#include <cstddef>
namespace dnav2 {
/* ---- BLAKE2b (RFC 7693), unkeyed, single-call ------------------------- */
inline constexpr uint64_t B2_IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL,
    0xa54ff53a5f1d36f1ULL, 0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
    0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL};
inline constexpr uint8_t B2_SIGMA[12][16] = {
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, {14,10,4,8,9,15,13,6,1,12,0,2,11,7,5,3},
    {11,8,12,0,5,2,15,13,10,14,3,6,7,1,9,4}, {7,9,3,1,13,12,11,14,2,6,5,10,4,0,15,8},
    {9,0,5,7,2,4,10,15,14,1,11,12,6,8,3,13}, {2,12,6,10,0,11,8,3,4,13,7,5,15,14,1,9},
    {12,5,1,15,14,13,4,10,0,7,6,3,9,2,8,11}, {13,11,7,14,12,1,3,9,5,0,15,4,8,6,2,10},
    {6,15,14,9,11,3,0,8,12,2,13,7,1,4,10,5}, {10,2,8,4,7,6,1,5,15,11,9,14,3,12,13,0},
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, {14,10,4,8,9,15,13,6,1,12,0,2,11,7,5,3}};
#define ROTR64(x, n) (((x) >> (n)) | ((x) << (64 - (n))))
#define B2G(a, b, c, d, x, y) do { \
    v[a] = v[a] + v[b] + (x); v[d] = ROTR64(v[d] ^ v[a], 32); \
    v[c] = v[c] + v[d];       v[b] = ROTR64(v[b] ^ v[c], 24); \
    v[a] = v[a] + v[b] + (y); v[d] = ROTR64(v[d] ^ v[a], 16); \
    v[c] = v[c] + v[d];       v[b] = ROTR64(v[b] ^ v[c], 63); } while (0)

inline void b2_compress(uint64_t h[8], const uint8_t blk[128], uint64_t t, int last) {
    uint64_t v[16], m[16];
    for (int i = 0; i < 16; i++) {
        uint64_t w = 0;
        for (int j = 7; j >= 0; j--) w = (w << 8) | blk[i * 8 + j];
        m[i] = w;
    }
    for (int i = 0; i < 8; i++) { v[i] = h[i]; v[i + 8] = B2_IV[i]; }
    v[12] ^= t;
    if (last) v[14] = ~v[14];
    for (int r = 0; r < 12; r++) {
        const uint8_t *s = B2_SIGMA[r];
        B2G(0, 4, 8, 12, m[s[0]], m[s[1]]);   B2G(1, 5, 9, 13, m[s[2]], m[s[3]]);
        B2G(2, 6, 10, 14, m[s[4]], m[s[5]]);  B2G(3, 7, 11, 15, m[s[6]], m[s[7]]);
        B2G(0, 5, 10, 15, m[s[8]], m[s[9]]);  B2G(1, 6, 11, 12, m[s[10]], m[s[11]]);
        B2G(2, 7, 8, 13, m[s[12]], m[s[13]]); B2G(3, 4, 9, 14, m[s[14]], m[s[15]]);
    }
    for (int i = 0; i < 8; i++) h[i] ^= v[i] ^ v[i + 8];
}

/* blake2b(digest_size=8) of msg, bytes -> big-endian uint64 (== Python
 * int.from_bytes(digest, "big")). */
inline uint64_t blake2b_key64(const uint8_t *msg, size_t len) {
    uint64_t h[8];
    for (int i = 0; i < 8; i++) h[i] = B2_IV[i];
    h[0] ^= 0x01010000ULL ^ 8ULL;
    uint8_t blk[128];
    size_t off = 0;
    while (len - off > 128) {
        memcpy(blk, msg + off, 128);
        off += 128;
        b2_compress(h, blk, (uint64_t)off, 0);
    }
    size_t rem = len - off;
    memset(blk, 0, 128);
    if (rem) memcpy(blk, msg + off, rem);
    b2_compress(h, blk, (uint64_t)len, 1);
    uint8_t out[8];
    for (int i = 0; i < 8; i++) out[i] = (uint8_t)(h[0] >> (8 * i));
    uint64_t key = 0;
    for (int i = 0; i < 8; i++) key = (key << 8) | out[i];
    return key;
}


}  // namespace dnav2
