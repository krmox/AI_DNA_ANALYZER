// Stage 2 / H1: BED-edge handling of build_windows on a synthetic contig (no external data).
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unistd.h>
#include "dnav2/frame.hpp"
using namespace dnav2;
static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { std::printf("FAIL: %s\n", msg); fails++; } } while (0)

static std::vector<int64_t> starts(const std::vector<Window>& w) {
    std::vector<int64_t> s;
    for (const auto& x : w) s.push_back(x.start);
    return s;
}

int main() {
    char dir[] = "/tmp/test_frame_XXXXXX";
    if (!mkdtemp(dir)) return 2;
    const std::string fasta = std::string(dir) + "/ref.fa", bed = std::string(dir) + "/t.bed";
    const int64_t len = 512;
    std::ofstream(fasta) << ">chr1\n" << std::string(len, 'A') << "\n";
    std::ofstream(bed) << "chr1\t100\t300\nchr1\t330\t340\n";

    std::vector<Window> strict, edge;
    CHECK(build_windows(fasta, bed, "chr1", len, strict).empty(), "strict build");
    CHECK(build_windows(fasta, bed, "chr1", len, edge, true).empty(), "edge build");
    CHECK((starts(strict) == std::vector<int64_t>{128, 192}), "strict keeps only windows fully inside the BED");
    CHECK((starts(edge) == std::vector<int64_t>{64, 128, 192, 256, 320}), "edge keeps every window touching the BED");

    // The defect: a BED-interior locus at 110 lies in window [64,128), which strict mode drops.
    auto scanned = [](const std::vector<Window>& w, int64_t pos) {
        for (const auto& x : w) if (pos >= x.start && pos < x.start + SEQ_LEN) return true;
        return false;
    };
    CHECK(!scanned(strict, 110), "strict mode loses locus 110 (documented H1 defect)");
    CHECK(scanned(edge, 110), "edge mode scans locus 110");
    CHECK(scanned(edge, 335), "edge mode scans locus 335 of the short second interval");

    const auto iv = load_bed(bed, "chr1");
    CHECK(!in_intervals(iv, 99) && in_intervals(iv, 100) && in_intervals(iv, 299) && !in_intervals(iv, 300), "interval 1 edges");
    CHECK(!in_intervals(iv, 329) && in_intervals(iv, 330) && in_intervals(iv, 339) && !in_intervals(iv, 340), "interval 2 edges");

    std::printf(fails ? "FRAME TESTS FAILED (%d)\n" : "FRAME TESTS PASSED\n", fails);
    std::system((std::string("rm -rf ") + dir).c_str());
    return fails ? 1 : 0;
}
