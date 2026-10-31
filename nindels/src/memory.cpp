#include "nindels/memory.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_set>

namespace nindels {

int PatternMemory::length_bin(int L) {
    if (L <= 1) return 0;
    if (L == 2) return 1;
    if (L <= 5) return 2;
    if (L <= 10) return 3;
    if (L <= 20) return 4;
    if (L <= 50) return 5;
    return 6;
}

int PatternMemory::hp_bin(int h) {
    if (h <= 3) return 0;
    if (h <= 7) return 1;
    if (h <= 12) return 2;
    return 3;
}

uint32_t PatternMemory::bucket_key(IndelType t, int lb, int hb, int period) {
    return (static_cast<uint32_t>(t) << 12) | (static_cast<uint32_t>(lb) << 8) |
           (static_cast<uint32_t>(hb) << 4) | static_cast<uint32_t>(period);
}

static uint32_t key_of(const Pattern& p) {
    return PatternMemory::bucket_key(p.type, PatternMemory::length_bin(p.length),
                                     PatternMemory::hp_bin(p.hp_length), std::min(p.period, 15));
}

PatternMemory PatternMemory::build(std::vector<LabeledPattern> items, const SimilarityParams& sim,
                                   const MemoryParams& mp, std::string metadata) {
    sim.validate();
    PatternMemory m;
    m.sim_ = sim;
    m.mp_ = mp;
    m.metadata_ = std::move(metadata);

    // Deterministic order: genomic position, then allele.
    std::sort(items.begin(), items.end(), [](const LabeledPattern& a, const LabeledPattern& b) {
        return std::tie(a.tid, a.pos, a.pattern.type, a.pattern.motif) <
               std::tie(b.tid, b.pos, b.pattern.type, b.pattern.motif);
    });

    struct Building {
        Pattern leader;
        std::vector<std::pair<uint64_t, Pattern>> reservoir;  // smallest hashes kept
        std::unordered_set<uint64_t> loci;
        uint32_t n_pos = 0, n_neg = 0;
    };
    std::vector<Building> work;
    std::unordered_map<uint32_t, std::vector<uint32_t>> leaders;

    for (const LabeledPattern& it : items) {
        const uint32_t key = key_of(it.pattern);
        auto& bucket = leaders[key];
        double best = std::numeric_limits<double>::infinity();
        uint32_t best_id = 0;
        for (uint32_t id : bucket) {
            double d = pattern_distance(it.pattern, work[id].leader, sim);
            if (d < best) { best = d; best_id = id; }
        }
        uint32_t cid;
        if (best < mp.theta || (!bucket.empty() && static_cast<int>(bucket.size()) >= mp.max_leaders_per_bucket)) {
            cid = best_id;
        } else {
            cid = static_cast<uint32_t>(work.size());
            work.push_back({it.pattern, {}, {}, 0, 0});
            bucket.push_back(cid);
        }
        Building& b = work[cid];
        (it.positive ? b.n_pos : b.n_neg) += 1;
        b.loci.insert(mix64(static_cast<uint64_t>(it.tid) << 40 ^
                            static_cast<uint64_t>(it.pos / mp.locus_merge_bp)));
        uint64_t h = hash_string(it.pattern.motif, mix64(static_cast<uint64_t>(it.tid) * 1000003ULL +
                                                         static_cast<uint64_t>(it.pos)));
        if (static_cast<int>(b.reservoir.size()) < mp.reservoir) {
            b.reservoir.emplace_back(h, it.pattern);
        } else {
            auto worst = std::max_element(b.reservoir.begin(), b.reservoir.end(),
                                          [](auto& x, auto& y) { return x.first < y.first; });
            if (h < worst->first) *worst = {h, it.pattern};
        }
    }

    m.classes_.reserve(work.size());
    for (Building& b : work) {
        // Medoid: reservoir member with the smallest summed distance to the others.
        std::sort(b.reservoir.begin(), b.reservoir.end(), [](auto& x, auto& y) { return x.first < y.first; });
        size_t best = 0;
        double best_sum = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < b.reservoir.size(); ++i) {
            double s = 0;
            for (size_t j = 0; j < b.reservoir.size(); ++j)
                if (i != j) s += pattern_distance(b.reservoir[i].second, b.reservoir[j].second, sim);
            if (s < best_sum) { best_sum = s; best = i; }
        }
        PatternClass c;
        c.medoid = b.reservoir.empty() ? b.leader : b.reservoir[best].second;
        c.n_pos = b.n_pos;
        c.n_neg = b.n_neg;
        c.loci = static_cast<uint32_t>(b.loci.size());
        m.classes_.push_back(std::move(c));
    }
    m.rebuild_index();
    m.estimate_alpha_beta();
    return m;
}

void PatternMemory::rebuild_index() {
    index_.clear();
    for (uint32_t i = 0; i < classes_.size(); ++i) index_[key_of(classes_[i].medoid)].push_back(i);
}

// alpha, beta by beta-binomial method of moments over classes with >= 5
// observations (spec §7 "selected using training data"). Falls back to
// alpha = beta = 1 (uniform prior) when the moments are not usable.
void PatternMemory::estimate_alpha_beta() {
    std::vector<double> p;
    for (const auto& c : classes_) {
        uint32_t n = c.n_pos + c.n_neg;
        if (n >= 5) p.push_back(static_cast<double>(c.n_pos) / n);
    }
    alpha_ = beta_ = 1.0;
    if (p.size() < 10) return;
    double mean = 0, var = 0;
    for (double x : p) mean += x;
    mean /= static_cast<double>(p.size());
    for (double x : p) var += (x - mean) * (x - mean);
    var /= static_cast<double>(p.size() - 1);
    if (mean <= 0 || mean >= 1 || var <= 0) return;
    double s = mean * (1 - mean) / var - 1;
    s = std::clamp(s, 1.0, 1000.0);  // ARBITRARY clamp to keep the prior proper and bounded
    alpha_ = mean * s;
    beta_ = (1 - mean) * s;
}

ContextStats PatternMemory::query(const Pattern& p) const {
    ContextStats out;
    const int lb = length_bin(p.length), hb = hp_bin(p.hp_length);
    std::vector<std::pair<uint32_t, double>> cand;
    // Indexed retrieval: same type, neighbouring length / homopolymer bins, any period.
    for (int l = std::max(0, lb - 1); l <= std::min(6, lb + 1); ++l)
        for (int h = std::max(0, hb - 1); h <= std::min(3, hb + 1); ++h)
            for (int per = 0; per <= 15; ++per) {
                auto it = index_.find(bucket_key(p.type, l, h, per));
                if (it == index_.end()) continue;
                for (uint32_t id : it->second)
                    cand.emplace_back(id, pattern_similarity(p, classes_[id].medoid, sim_));
            }
    auto better = [](const auto& a, const auto& b) {
        return a.second != b.second ? a.second > b.second : a.first < b.first;
    };
    size_t k = std::min(cand.size(), static_cast<size_t>(mp_.top_k));
    std::partial_sort(cand.begin(), cand.begin() + static_cast<long>(k), cand.end(), better);
    cand.resize(k);

    double W = 0, W2 = 0, Wl = 0, Wl2 = 0;
    for (auto [id, S] : cand) {
        const PatternClass& c = classes_[id];
        const double n = c.n_pos + c.n_neg;
        out.n_positive += S * c.n_pos;
        out.n_negative += S * c.n_neg;
        W += S * n;
        W2 += S * S * n;
        Wl += S * c.loci;
        Wl2 += S * S * c.loci;
    }
    out.retrieved = static_cast<int>(k);
    out.p_context = (out.n_positive + alpha_) / (out.n_positive + out.n_negative + alpha_ + beta_);
    if (W2 > 0) out.n_eff_kish = W * W / W2;
    if (Wl2 > 0) out.n_eff_loci = Wl * Wl / Wl2;
    if (k > 0) out.s_max = cand.front().second;
    out.closest = std::move(cand);
    return out;
}

// ---------------------------------------------------------------- persistence

static std::string or_dot(const std::string& s) { return s.empty() ? "." : s; }
static std::string from_dot(const std::string& s) { return s == "." ? std::string() : s; }

void PatternMemory::save(const std::string& path) const {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write memory file " + path);
    o.precision(17);
    o << "#nindels-memory\tv1\n";
    o << "#metadata\t" << or_dot(metadata_) << "\n";
    o << "#similarity\t" << sim_.w_type << ' ' << sim_.w_length << ' ' << sim_.w_motif << ' '
      << sim_.w_context << ' ' << sim_.w_hp << ' ' << sim_.w_repeat << ' ' << sim_.s_L << ' ' << sim_.s_h
      << ' ' << sim_.s_r << "\n";
    o << "#memory\t" << mp_.theta << ' ' << mp_.top_k << ' ' << mp_.reservoir << ' '
      << mp_.max_leaders_per_bucket << ' ' << mp_.locus_merge_bp << "\n";
    o << "#prior\t" << alpha_ << ' ' << beta_ << "\n";
    o << "#type\tlength\tmotif\tcontext\thp_length\thp_base\tperiod\tcopies\trepeat_span\tentropy\tn_pos\tn_neg\tloci\n";
    for (const auto& c : classes_) {
        const Pattern& p = c.medoid;
        o << (p.type == IndelType::Insertion ? "INS" : "DEL") << '\t' << p.length << '\t' << or_dot(p.motif)
          << '\t' << or_dot(p.context) << '\t' << p.hp_length << '\t' << p.hp_base << '\t' << p.period << '\t'
          << p.copies << '\t' << p.repeat_span << '\t' << p.entropy << '\t' << c.n_pos << '\t' << c.n_neg
          << '\t' << c.loci << '\n';
    }
}

PatternMemory PatternMemory::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read memory file " + path);
    PatternMemory m;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        if (line[0] == '#') {
            std::string tag;
            ss >> tag;
            if (tag == "#metadata") { std::string rest; std::getline(ss >> std::ws, rest); m.metadata_ = from_dot(rest); }
            else if (tag == "#similarity")
                ss >> m.sim_.w_type >> m.sim_.w_length >> m.sim_.w_motif >> m.sim_.w_context >> m.sim_.w_hp >>
                    m.sim_.w_repeat >> m.sim_.s_L >> m.sim_.s_h >> m.sim_.s_r;
            else if (tag == "#memory")
                ss >> m.mp_.theta >> m.mp_.top_k >> m.mp_.reservoir >> m.mp_.max_leaders_per_bucket >>
                    m.mp_.locus_merge_bp;
            else if (tag == "#prior") ss >> m.alpha_ >> m.beta_;
            continue;
        }
        PatternClass c;
        std::string type, motif, ctx;
        ss >> type >> c.medoid.length >> motif >> ctx >> c.medoid.hp_length >> c.medoid.hp_base >>
            c.medoid.period >> c.medoid.copies >> c.medoid.repeat_span >> c.medoid.entropy >> c.n_pos >>
            c.n_neg >> c.loci;
        if (!ss) throw std::runtime_error("malformed memory line: " + line);
        c.medoid.type = type == "INS" ? IndelType::Insertion : IndelType::Deletion;
        c.medoid.motif = from_dot(motif);
        c.medoid.context = from_dot(ctx);
        m.classes_.push_back(std::move(c));
    }
    m.sim_.validate();
    m.rebuild_index();
    return m;
}

}  // namespace nindels
