#include "dnav2/frame.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"

namespace dnav2 {

struct Fai { faidx_t* f; };
Fai* fai_open(const std::string& fasta) {
    faidx_t* f = fai_load(fasta.c_str());
    return f ? new Fai{f} : nullptr;
}
void fai_close(Fai* f) { if (f) { fai_destroy(f->f); delete f; } }

std::string fetch_ref_idx(Fai* fai, const std::string& contig, int64_t start, int64_t end, uint8_t* out,
                          int* n_count) {
    hts_pos_t len = 0;
    char* s = faidx_fetch_seq64(fai->f, contig.c_str(), start, end - 1, &len);
    if (!s || len != end - start) { free(s); return "faidx fetch failed"; }
    int nn = 0;
    for (int64_t i = 0; i < len; i++) {
        char c = static_cast<char>(s[i] >= 'a' && s[i] <= 'z' ? s[i] - 32 : s[i]);   // .upper()
        if (c == 'N') nn++;
        out[i] = c == 'A' ? 1 : c == 'C' ? 2 : c == 'G' ? 3 : c == 'T' ? 4 : 0;
    }
    free(s);
    if (n_count) *n_count = nn;
    return "";
}

std::vector<std::pair<int64_t, int64_t>> load_bed(const std::string& path, const std::string& contig) {
    std::ifstream in(path);
    std::vector<std::pair<int64_t, int64_t>> iv;
    std::string line;
    const std::string bare = contig.rfind("chr", 0) == 0 ? contig.substr(3) : contig;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
        std::istringstream ss(line);
        std::string c; int64_t s, e;
        if (!(ss >> c >> s >> e)) continue;
        std::string b = c.rfind("chr", 0) == 0 ? c.substr(3) : c;
        if (b != bare) continue;
        iv.emplace_back(s, e);
    }
    std::sort(iv.begin(), iv.end());
    std::vector<std::pair<int64_t, int64_t>> merged;
    for (auto& p : iv) {
        if (!merged.empty() && p.first <= merged.back().second)
            merged.back().second = std::max(merged.back().second, p.second);
        else merged.push_back(p);
    }
    return merged;
}

std::string build_windows(const std::string& fasta, const std::string& bed, const std::string& contig,
                          int64_t contig_len, std::vector<Window>& out) {
    auto iv = load_bed(bed, contig);
    if (iv.empty()) return "no BED intervals for contig";
    Fai* fai = fai_open(fasta);
    if (!fai) return "cannot open FASTA (.fai required): " + fasta;
    std::vector<uint8_t> buf(SEQ_LEN);
    size_t k = 0;   // merged intervals are disjoint & sorted: containment == in the interval that
                    // starts at or before the window (V1's prefix-max test)
    for (int64_t start = 0; start <= contig_len - SEQ_LEN; start += SEQ_LEN) {
        int64_t end = start + SEQ_LEN;
        while (k + 1 < iv.size() && iv[k + 1].first <= start) k++;
        if (!(iv[k].first <= start && iv[k].second >= end)) continue;
        int nn = 0;
        std::string e = fetch_ref_idx(fai, contig, start, end, buf.data(), &nn);
        if (!e.empty()) { fai_close(fai); return e; }
        if (nn > SEQ_LEN / 2) continue;
        out.push_back({start});
    }
    fai_close(fai);
    return "";
}

static void split(const char* s, char d, std::vector<std::string>& out) {
    out.clear();
    std::string cur;
    for (; *s; s++) { if (*s == d) { out.push_back(cur); cur.clear(); } else cur.push_back(*s); }
    out.push_back(cur);
}

std::string load_labels(const std::string& truth_vcf, const std::string& contig, int64_t contig_len,
                        std::vector<uint8_t>& labels) {
    labels.assign(static_cast<size_t>(contig_len), LABEL_NORMAL);
    htsFile* fp = hts_open(truth_vcf.c_str(), "r");
    if (!fp) return "cannot open truth VCF";
    tbx_t* tbx = tbx_index_load(truth_vcf.c_str());
    if (!tbx) { hts_close(fp); return "cannot load truth VCF tabix index"; }
    hts_itr_t* it = tbx_itr_querys(tbx, contig.c_str());
    if (!it) { tbx_destroy(tbx); hts_close(fp); return "contig not in truth VCF"; }
    kstring_t ks = {0, 0, nullptr};
    std::vector<std::string> f, alts, fmt, smp, al;
    while (tbx_itr_next(fp, tbx, it, &ks) >= 0) {
        split(ks.s, '\t', f);
        if (f.size() != 10) { free(ks.s); return "truth VCF record without exactly one sample"; }
        const int64_t start = std::atoll(f[1].c_str()) - 1;
        const std::string& ref = f[3];
        if (f[4] == ".") continue;                       // record.alts is None
        // homozygous-reference sample -> skipped (providers._is_homozygous_reference)
        split(f[8].c_str(), ':', fmt);
        split(f[9].c_str(), ':', smp);
        int gi = -1;
        for (size_t i = 0; i < fmt.size(); i++) if (fmt[i] == "GT") { gi = static_cast<int>(i); break; }
        if (gi >= 0 && static_cast<size_t>(gi) < smp.size()) {
            std::string g = smp[gi];
            std::replace(g.begin(), g.end(), '|', '/');
            split(g.c_str(), '/', al);
            bool hom_ref = true;
            for (auto& a : al) if (!(a == "0" || a == ".")) hom_ref = false;
            if (hom_ref) continue;
        }
        split(f[4].c_str(), ',', alts);
        for (const auto& alt : alts) {
            if (alt.empty() || alt[0] == '<') continue;
            if (ref.empty()) continue;
            int label, off, span;
            if (ref.size() == alt.size()) { label = LABEL_SNP; off = 0; span = static_cast<int>(ref.size()); }
            else if (alt.size() > ref.size()) { label = LABEL_INSERTION; off = 0; span = 1; }
            else { label = LABEL_DELETION; off = 1; span = std::max<int>(1, static_cast<int>(ref.size() - alt.size())); }
            for (int o = off; o < off + span; o++) {
                int64_t p = start + o;
                if (p >= 0 && p < contig_len) labels[static_cast<size_t>(p)] = static_cast<uint8_t>(label);
            }
        }
    }
    free(ks.s);
    tbx_itr_destroy(it);
    tbx_destroy(tbx);
    hts_close(fp);
    return "";
}

}  // namespace dnav2
