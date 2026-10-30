#include "nindels/truth.hpp"

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "nindels/pattern.hpp"

namespace nindels {

bool BedRegions::contains(const Indel& x) const {
    const int64_t a = x.pos - 1, b = x.pos + x.ref_span();  // [a, b)
    auto it = std::upper_bound(iv.begin(), iv.end(), std::make_pair(a, INT64_MAX));
    if (it == iv.begin()) return false;
    --it;
    return a >= it->first && b <= it->second;
}

BedRegions load_bed(const std::string& path, const std::string& chrom, int64_t start, int64_t end) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read BED " + path);
    BedRegions r;
    std::string line, c;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0) continue;
        std::istringstream ss(line);
        int64_t s, e;
        if (!(ss >> c >> s >> e) || c != chrom) continue;
        s = std::max(s, start);
        e = std::min(e, end);
        if (e > s) r.iv.emplace_back(s, e);
    }
    std::sort(r.iv.begin(), r.iv.end());
    return r;
}

bool allele_to_indel(int64_t pos0, std::string ref, std::string alt, int32_t tid, Indel& out, bool& complex) {
    complex = false;
    if (ref.size() == alt.size()) return false;
    if (alt.find_first_not_of("ACGTacgt") != std::string::npos || ref.find_first_not_of("ACGTNacgtn") != std::string::npos)
        return false;  // symbolic / breakend / missing alleles
    std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
    std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);
    while (ref.size() > 1 && alt.size() > 1 && ref.back() == alt.back()) { ref.pop_back(); alt.pop_back(); }
    while (ref.size() > 1 && alt.size() > 1 && ref[0] == alt[0]) { ref.erase(0, 1); alt.erase(0, 1); ++pos0; }
    if (ref.size() == 1 && alt.size() > 1 && ref[0] == alt[0]) {
        out = Indel{tid, pos0 + 1, IndelType::Insertion, alt.substr(1)};
        return true;
    }
    if (alt.size() == 1 && ref.size() > 1 && ref[0] == alt[0]) {
        out = Indel{tid, pos0 + 1, IndelType::Deletion, ref.substr(1)};
        return true;
    }
    complex = true;
    return false;
}

namespace {

template <class F>
void for_each_record(const std::string& vcf, const std::string& chrom, int64_t start, int64_t end, F&& f) {
    bcf_srs_t* sr = bcf_sr_init();
    std::string region = chrom + ":" + std::to_string(start + 1) + "-" + std::to_string(end);
    if (bcf_sr_set_regions(sr, region.c_str(), 0) < 0 || !bcf_sr_add_reader(sr, vcf.c_str())) {
        bcf_sr_destroy(sr);
        throw std::runtime_error("cannot open indexed VCF " + vcf);
    }
    bcf_hdr_t* hdr = bcf_sr_get_header(sr, 0);
    while (bcf_sr_next_line(sr)) f(hdr, bcf_sr_get_line(sr, 0));
    bcf_sr_destroy(sr);
}

}  // namespace

TruthSet load_truth(const std::string& vcf, const std::string& chrom, int64_t start, int64_t end, int32_t tid,
                    const RefSeq& ref, const BedRegions* bed) {
    TruthSet t;
    int32_t* gt = nullptr;
    int ngt = 0;
    for_each_record(vcf, chrom, start, end, [&](bcf_hdr_t* hdr, bcf1_t* rec) {
        bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT);
        std::vector<int> used;
        if (bcf_hdr_nsamples(hdr) > 0 && bcf_get_genotypes(hdr, rec, &gt, &ngt) > 0) {
            for (int i = 0; i < std::min(ngt, 2); ++i)
                if (!bcf_gt_is_missing(gt[i]) && bcf_gt_allele(gt[i]) > 0) used.push_back(bcf_gt_allele(gt[i]));
        } else {
            for (int a = 1; a < rec->n_allele; ++a) used.push_back(a);
        }
        std::sort(used.begin(), used.end());
        used.erase(std::unique(used.begin(), used.end()), used.end());
        for (int a : used) {
            Indel x;
            bool cx = false;
            if (!allele_to_indel(rec->pos, rec->d.allele[0], rec->d.allele[a], tid, x, cx)) {
                if (cx && (!bed || bed->contains(Indel{tid, rec->pos + 1, IndelType::Insertion, "N"}))) ++t.n_complex;
                continue;
            }
            x = left_normalize(ref, x);
            if (x.pos < start || x.pos >= end) continue;
            if (t.indels.insert(x).second && (!bed || bed->contains(x))) ++t.n_in_bed;
        }
    });
    free(gt);
    return t;
}

int truth_label(const TruthSet& t, const BedRegions* bed, const Indel& x) {
    if (bed && !bed->contains(x)) return -1;
    return t.indels.count(x) ? 1 : 0;
}

std::vector<Indel> load_vcf_indels(const std::string& vcf, const std::string& chrom, int64_t start, int64_t end,
                                   int32_t tid, const RefSeq& ref) {
    std::vector<Indel> out;
    for_each_record(vcf, chrom, start, end, [&](bcf_hdr_t*, bcf1_t* rec) {
        bcf_unpack(rec, BCF_UN_STR);
        for (int a = 1; a < rec->n_allele; ++a) {
            Indel x;
            bool cx;
            if (allele_to_indel(rec->pos, rec->d.allele[0], rec->d.allele[a], tid, x, cx))
                out.push_back(left_normalize(ref, x));
        }
    });
    return out;
}

}  // namespace nindels
