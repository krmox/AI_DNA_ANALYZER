#include "nindels/bam_scan.hpp"

#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_map>

namespace nindels {

RefSeq load_reference(const std::string& fasta, const std::string& chrom, int64_t start, int64_t end) {
    faidx_t* fai = fai_load(fasta.c_str());
    if (!fai) throw std::runtime_error("cannot load fasta index for " + fasta);
    hts_pos_t len = 0;
    char* s = faidx_fetch_seq64(fai, chrom.c_str(), start, end - 1, &len);
    fai_destroy(fai);
    if (!s || len <= 0) {
        free(s);
        throw std::runtime_error("cannot fetch " + chrom + " from " + fasta);
    }
    RefSeq r;
    r.offset = start;
    r.seq.assign(s, static_cast<size_t>(len));
    free(s);
    for (char& c : r.seq) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
    }
    return r;
}

namespace {

struct Event {
    Indel indel;    // normalized
    int qpos = 0;   // query position of the event in the read (before normalization)
};

struct AlnRead {
    int64_t start = 0, end = 0;
    uint8_t mapq = 0;
    bool rev = false, split = false;
    int clip_l = 0, clip_r = 0, qlen = 0;
    uint64_t name_hash = 0;
    uint64_t dup_key = 0;
    std::string seq;
    std::vector<uint8_t> qual;
    std::vector<uint32_t> cigar;
    std::vector<Event> events;
};

// Query position aligned to reference position g, or -1 (deleted / not covered).
int qpos_at(const AlnRead& r, int64_t g) {
    int64_t rpos = r.start;
    int q = 0;
    for (uint32_t c : r.cigar) {
        const int op = bam_cigar_op(c), len = static_cast<int>(bam_cigar_oplen(c));
        const int t = bam_cigar_type(op);  // bit 1: consumes query, bit 2: consumes reference
        if ((t & 1) && (t & 2)) {
            if (g >= rpos && g < rpos + len) return q + static_cast<int>(g - rpos);
            rpos += len;
            q += len;
        } else if (t & 2) {
            if (g >= rpos && g < rpos + len) return -1;
            rpos += len;
        } else if (t & 1) {
            q += len;
        }
        if (rpos > g) break;
    }
    return -1;
}

// Mismatches / aligned bases over reference [a, b).
float mismatch_rate(const AlnRead& r, const RefSeq& ref, int64_t a, int64_t b) {
    int64_t rpos = r.start;
    int q = 0, n = 0, mm = 0;
    for (uint32_t c : r.cigar) {
        const int op = bam_cigar_op(c), len = static_cast<int>(bam_cigar_oplen(c));
        const int t = bam_cigar_type(op);
        if ((t & 1) && (t & 2)) {
            for (int i = 0; i < len; ++i) {
                int64_t g = rpos + i;
                if (g < a || g >= b) continue;
                char rb = ref.at(g), qb = r.seq[static_cast<size_t>(q + i)];
                if (rb == 'N' || qb == 'N') continue;
                ++n;
                mm += rb != qb;
            }
            rpos += len;
            q += len;
        } else if (t & 2) {
            rpos += len;
        } else if (t & 1) {
            q += len;
        }
        if (rpos >= b) break;
    }
    return n ? static_cast<float>(mm) / static_cast<float>(n) : 0.0f;
}

void load_read(const bam1_t* b, const RefSeq& ref, int32_t tid, AlnRead& r) {
    r.start = b->core.pos;
    r.end = bam_endpos(b);
    r.mapq = b->core.qual;
    r.rev = bam_is_rev(b);
    r.split = bam_aux_get(b, "SA") != nullptr;
    r.qlen = b->core.l_qseq;
    r.name_hash = hash_string(bam_get_qname(b));
    r.dup_key = mix64(mix64(mix64(static_cast<uint64_t>(b->core.pos) * 2 + (bam_is_rev(b) ? 1 : 0)) ^
                            static_cast<uint64_t>(b->core.mpos + 1)) ^
                      (static_cast<uint64_t>(b->core.isize) * 4 + ((b->core.flag & BAM_FREAD1) ? 1 : 0) +
                       ((b->core.flag & BAM_FREAD2) ? 2 : 0)));
    const uint8_t* s = bam_get_seq(b);
    const uint8_t* q = bam_get_qual(b);
    r.seq.resize(static_cast<size_t>(r.qlen));
    r.qual.assign(q, q + r.qlen);
    for (int i = 0; i < r.qlen; ++i) r.seq[static_cast<size_t>(i)] = seq_nt16_str[bam_seqi(s, i)];
    const uint32_t* cig = bam_get_cigar(b);
    r.cigar.assign(cig, cig + b->core.n_cigar);
    r.events.clear();
    r.clip_l = r.clip_r = 0;

    int64_t rpos = r.start;
    int qp = 0;
    bool seen_aligned = false;
    for (uint32_t c : r.cigar) {
        const int op = bam_cigar_op(c), len = static_cast<int>(bam_cigar_oplen(c));
        switch (op) {
            case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                rpos += len; qp += len; seen_aligned = true; break;
            case BAM_CINS: {
                std::string ins = r.seq.substr(static_cast<size_t>(qp), static_cast<size_t>(len));
                if (ins.find('N') == std::string::npos)
                    r.events.push_back({left_normalize(ref, Indel{tid, rpos, IndelType::Insertion, ins}), qp});
                qp += len;
                break;
            }
            case BAM_CDEL: {
                std::string del = ref.sub(rpos, len);
                if (static_cast<int>(del.size()) == len && del.find('N') == std::string::npos)
                    r.events.push_back({left_normalize(ref, Indel{tid, rpos, IndelType::Deletion, del}), qp});
                rpos += len;
                break;
            }
            case BAM_CREF_SKIP: rpos += len; break;
            case BAM_CSOFT_CLIP:
                (seen_aligned ? r.clip_r : r.clip_l) += len;
                qp += len;
                break;
            default: break;
        }
    }
}

}  // namespace

struct BamScanner::Impl {
    samFile* fp = nullptr;
    sam_hdr_t* hdr = nullptr;
    hts_idx_t* idx = nullptr;
};

BamScanner::BamScanner(const std::string& bam_path, const std::string& ref_fasta) : impl_(new Impl) {
    impl_->fp = sam_open(bam_path.c_str(), "r");
    if (!impl_->fp) throw std::runtime_error("cannot open " + bam_path);
    if (!ref_fasta.empty()) hts_set_fai_filename(impl_->fp, ref_fasta.c_str());
    impl_->hdr = sam_hdr_read(impl_->fp);
    impl_->idx = sam_index_load(impl_->fp, bam_path.c_str());
    if (!impl_->hdr || !impl_->idx) throw std::runtime_error("cannot read header/index of " + bam_path);
}

BamScanner::~BamScanner() {
    if (impl_->idx) hts_idx_destroy(impl_->idx);
    if (impl_->hdr) sam_hdr_destroy(impl_->hdr);
    if (impl_->fp) sam_close(impl_->fp);
}

int BamScanner::tid(const std::string& chrom) const { return sam_hdr_name2tid(impl_->hdr, chrom.c_str()); }

std::vector<std::pair<std::string, int64_t>> BamScanner::contigs() const {
    std::vector<std::pair<std::string, int64_t>> out;
    for (int i = 0; i < sam_hdr_nref(impl_->hdr); ++i)
        out.emplace_back(sam_hdr_tid2name(impl_->hdr, i), sam_hdr_tid2len(impl_->hdr, i));
    return out;
}

std::string BamScanner::sample_name() const {
    kstring_t ks = {0, 0, nullptr};
    std::string out = "SAMPLE";
    if (sam_hdr_find_tag_pos(impl_->hdr, "RG", 0, "SM", &ks) == 0) out.assign(ks.s, ks.l);
    free(ks.s);
    return out;
}

StrLocusCount BamScanner::count_str_locus(const StrRun& run, int32_t tid, int min_mapq) {
    StrLocusCount c;
    c.period = run.period;
    c.units = run.units();
    const int64_t a = run.start - 3, e = run.end + 3;
    hts_itr_t* itr = sam_itr_queryi(impl_->idx, tid, a, e);
    if (!itr) return c;
    bam1_t* b = bam_init1();
    std::unordered_map<uint64_t, bool> frag;  // name hash -> slipped (first mate seen wins)
    while (sam_itr_next(impl_->fp, itr, b) >= 0) {
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) continue;
        if (b->core.qual < min_mapq || b->core.pos > a || bam_endpos(b) < e) continue;
        const uint64_t h = hash_string(bam_get_qname(b));
        if (frag.count(h)) continue;
        bool slipped = false;
        int64_t rpos = b->core.pos;
        const uint32_t* cig = bam_get_cigar(b);
        for (uint32_t k = 0; k < b->core.n_cigar; ++k) {
            const int op = bam_cigar_op(cig[k]), len = static_cast<int>(bam_cigar_oplen(cig[k]));
            if ((op == BAM_CINS || op == BAM_CDEL) && rpos >= run.start - 1 && rpos <= run.end + 1 &&
                len % run.period == 0)
                slipped = true;
            if (bam_cigar_type(op) & 2) rpos += len;
        }
        frag.emplace(h, slipped);
    }
    bam_destroy1(b);
    hts_itr_destroy(itr);
    c.n = static_cast<int>(frag.size());
    for (auto& [h, sl] : frag) c.s += sl;
    return c;
}

std::vector<CandidateInput> BamScanner::scan(const RefSeq& ref, int32_t tid, int64_t start, int64_t end,
                                             const ScanParams& p, const std::vector<Indel>* external,
                                             int window_half) {
    // 1. Load reads overlapping the chunk (+ padding).
    std::vector<AlnRead> reads;
    hts_itr_t* itr = sam_itr_queryi(impl_->idx, tid, std::max<int64_t>(0, start - p.read_pad), end + p.read_pad);
    if (!itr) throw std::runtime_error("BAM iterator failed");
    bam1_t* b = bam_init1();
    int64_t max_span = 0;
    while (sam_itr_next(impl_->fp, itr, b) >= 0) {
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) continue;
        if (b->core.n_cigar == 0 || b->core.l_qseq == 0) continue;
        reads.emplace_back();
        load_read(b, ref, tid, reads.back());
        max_span = std::max(max_span, reads.back().end - reads.back().start);
    }
    bam_destroy1(b);
    hts_itr_destroy(itr);

    if (p.fragments) {
        // Unmarked duplicates: among reads with the same duplicate key keep the
        // smallest name hash. Both mates of a duplicate pair make the same
        // choice, so whole fragments are kept or dropped together.
        std::unordered_map<uint64_t, uint64_t> best;
        for (const AlnRead& r : reads) {
            auto [it, fresh] = best.emplace(r.dup_key, r.name_hash);
            if (!fresh) it->second = std::min(it->second, r.name_hash);
        }
        std::vector<AlnRead> kept;
        kept.reserve(reads.size());
        for (AlnRead& r : reads)
            if (best[r.dup_key] == r.name_hash) kept.push_back(std::move(r));
        reads.swap(kept);
    }

    // 2. Candidate alleles: normalized events with >= min_alt supporting reads.
    std::vector<Indel> cands;
    if (external) {
        for (const Indel& x : *external)
            if (x.pos >= start && x.pos < end) cands.push_back(left_normalize(ref, x));
    } else {
        std::unordered_map<Indel, int, IndelHash> count;
        for (const AlnRead& r : reads)
            for (const Event& e : r.events)
                if (e.indel.pos >= start && e.indel.pos < end) ++count[e.indel];
        for (auto& [x, n] : count)
            if (n >= p.min_alt) cands.push_back(x);
    }
    std::sort(cands.begin(), cands.end());
    cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

    // 3. Evidence for each candidate from the reads overlapping it.
    std::vector<int64_t> starts(reads.size());
    for (size_t i = 0; i < reads.size(); ++i) starts[i] = reads[i].start;

    std::vector<CandidateInput> out;
    out.reserve(cands.size());
    const int W = 2 * window_half + 1;
    for (const Indel& x : cands) {
        Pattern pat = extract_pattern(ref, x);
        const int64_t rs = x.pos - p.span_pad;
        const int64_t re = x.pos + std::max({x.ref_span(), pat.repeat_span, pat.hp_length}) + p.span_pad;

        CandidateInput ci;
        ci.indel = x;
        size_t hi = static_cast<size_t>(std::lower_bound(starts.begin(), starts.end(), re) - starts.begin());
        for (size_t k = hi; k-- > 0;) {
            const AlnRead& r = reads[k];
            if (r.start < rs - max_span) break;
            if (r.end <= rs) continue;

            ReadAllele allele;
            int qpos = -1;
            bool overlapping_other = false, has_alt = false;
            for (const Event& e : r.events) {
                if (e.indel == x) { has_alt = true; qpos = e.qpos; break; }
                if (e.indel.pos < re && e.indel.pos + std::max(1, e.indel.ref_span()) > rs) overlapping_other = true;
            }
            if (has_alt) allele = ReadAllele::Alt;
            else if (overlapping_other) allele = ReadAllele::Other;
            else if (r.start <= rs && r.end >= re) allele = ReadAllele::Ref;
            else continue;  // read does not span the event: uninformative

            if (qpos < 0) qpos = qpos_at(r, x.pos);
            if (qpos < 0) qpos = r.qlen / 2;

            ReadObservation o;
            o.allele = allele;
            o.reverse = r.rev;
            o.mapq = r.mapq;
            o.rel_pos = static_cast<float>(std::min(qpos, r.qlen - qpos)) / static_cast<float>(std::max(1, r.qlen));
            o.soft_clip = r.clip_l + r.clip_r;
            int q0 = std::max(0, qpos - p.bq_flank);
            int q1 = std::min(r.qlen, qpos + p.bq_flank + (x.type == IndelType::Insertion ? x.length() : 0));
            double sbq = 0;
            for (int i = q0; i < q1; ++i) sbq += r.qual[static_cast<size_t>(i)];
            o.flank_bq = q1 > q0 ? static_cast<float>(sbq / (q1 - q0)) : 0.0f;
            o.mismatch_rate = mismatch_rate(r, ref, x.pos - p.mismatch_flank, x.pos + x.ref_span() + p.mismatch_flank);
            o.split = r.split;
            ci.obs.push_back(o);

            HmmRead h;
            const size_t a = static_cast<size_t>(r.clip_l), n = static_cast<size_t>(r.qlen - r.clip_l - r.clip_r);
            h.bases = r.seq.substr(a, n);
            h.quals.assign(r.qual.begin() + static_cast<long>(a), r.qual.begin() + static_cast<long>(a + n));
            h.ref_start = r.start;
            h.ref_end = r.end;
            h.mapq = r.mapq;
            h.cigar_allele = allele;
            h.name_hash = r.name_hash;
            h.window_bq.assign(static_cast<size_t>(W), -1);
            for (int k2 = 0; k2 < W; ++k2) {
                int qq = qpos_at(r, x.pos - window_half + k2);
                if (qq >= 0) h.window_bq[static_cast<size_t>(k2)] = r.qual[static_cast<size_t>(qq)];
            }
            ci.reads.push_back(std::move(h));
        }
        // Evidence was collected newest-first; restore read order for determinism.
        std::reverse(ci.obs.begin(), ci.obs.end());
        std::reverse(ci.reads.begin(), ci.reads.end());

        if (p.fragments) {
            // One observation per fragment: overlapping mates share the
            // molecule (and its PCR errors). Keep the mate with the higher
            // flanking base quality; ties keep the first in read order.
            std::unordered_map<uint64_t, size_t> first;
            std::vector<char> drop(ci.obs.size(), 0);
            for (size_t k2 = 0; k2 < ci.obs.size(); ++k2) {
                auto [it, fresh] = first.emplace(ci.reads[k2].name_hash, k2);
                if (fresh) continue;
                size_t j = it->second;
                const bool better = p.fragment_pick_by_end
                    ? (ci.obs[k2].rel_pos != ci.obs[j].rel_pos ? ci.obs[k2].rel_pos > ci.obs[j].rel_pos
                                                               : ci.obs[k2].flank_bq > ci.obs[j].flank_bq)
                    : ci.obs[k2].flank_bq > ci.obs[j].flank_bq;
                if (better) { drop[j] = 1; it->second = k2; }
                else drop[k2] = 1;
            }
            size_t w = 0;
            for (size_t k2 = 0; k2 < ci.obs.size(); ++k2)
                if (!drop[k2]) {
                    ci.obs[w] = ci.obs[k2];
                    ci.reads[w] = std::move(ci.reads[k2]);
                    ++w;
                }
            ci.obs.resize(w);
            ci.reads.resize(w);
        }

        if (!external) {
            int n_alt = 0, n_ref = 0;
            for (const auto& o : ci.obs) { n_alt += o.allele == ReadAllele::Alt; n_ref += o.allele == ReadAllele::Ref; }
            if (n_alt < p.min_alt) continue;
            if (static_cast<double>(n_alt) / (n_alt + n_ref) < p.min_alt_frac) continue;
        }
        out.push_back(std::move(ci));
    }
    return out;
}

}  // namespace nindels
