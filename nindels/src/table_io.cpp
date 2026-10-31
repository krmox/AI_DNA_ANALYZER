#include "nindels/table_io.hpp"

#include <charconv>
#include <fstream>
#include <ostream>
#include <stdexcept>

namespace nindels {

std::vector<std::string> split_tab(std::string_view line) {
    std::vector<std::string> out;
    size_t b = 0;
    while (true) {
        size_t e = line.find('\t', b);
        out.emplace_back(line.substr(b, e == std::string_view::npos ? std::string_view::npos : e - b));
        if (e == std::string_view::npos) break;
        b = e + 1;
    }
    return out;
}

Opt parse_opt(std::string_view s) {
    if (s.empty() || s == "NA") return std::nullopt;
    return std::stod(std::string(s));
}

std::string format_opt(Opt v) {
    if (!v) return "NA";
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof buf, *v);
    return std::string(buf, res.ptr);
}

VcfAlleles vcf_alleles(const RefSeq& ref, const Indel& x) {
    VcfAlleles a;
    const char anchor = ref.at(x.pos - 1);
    a.pos1 = x.pos;  // 0-based anchor index x.pos-1 == 1-based x.pos
    if (x.type == IndelType::Deletion) {
        a.ref = std::string(1, anchor) + x.seq;
        a.alt = std::string(1, anchor);
    } else {
        a.ref = std::string(1, anchor);
        a.alt = std::string(1, anchor) + x.seq;
    }
    return a;
}

static const char* kFixedColumns[] = {"chrom",  "pos",     "ref",    "alt",         "type",    "length",
                                      "motif",  "context", "hp_length", "hp_base",   "period",  "copies",
                                      "repeat_span", "entropy", "label", "n_alt", "n_ref", "n_other", "n_total"};
static const char* kTailColumns[] = {"p_evidence", "p_base", "U", "D", "R", "tau", "reassessed", "hmm_reads",
                                     "q_argmin", "gt", "gq", "p_final", "ood", "decision", "pos0"};

void write_table_header(std::ostream& o, const std::vector<std::string>& comments) {
    for (const auto& c : comments) o << "##" << c << '\n';
    bool first = true;
    auto col = [&](const std::string& n) { o << (first ? "" : "\t") << n; first = false; };
    for (const char* c : kFixedColumns) col(c);
    for (const auto& n : all_feature_names()) col(n);
    for (const char* c : kTailColumns) col(c);
    o << '\n';
}

void write_table_row(std::ostream& o, const std::string& chrom, const RefSeq& ref, const CandidateResult& r,
                     int label) {
    const VcfAlleles a = vcf_alleles(ref, r.indel);
    const Pattern& p = r.pattern;
    const Evidence& e = r.evidence;
    o << chrom << '\t' << a.pos1 << '\t' << a.ref << '\t' << a.alt << '\t'
      << (p.type == IndelType::Insertion ? "INS" : "DEL") << '\t' << p.length << '\t' << p.motif << '\t'
      << (p.context.empty() ? "." : p.context) << '\t' << p.hp_length << '\t' << p.hp_base << '\t' << p.period
      << '\t' << p.copies << '\t' << p.repeat_span << '\t' << format_opt(p.entropy) << '\t' << label << '\t'
      << e.n_alt << '\t' << e.n_ref << '\t' << e.n_other << '\t' << e.n_total;
    for (const auto& f : r.features) o << '\t' << format_opt(f.value);
    const GenotypeLik& g = r.hmm_genotype ? *r.hmm_genotype : r.count_genotype;
    o << '\t' << format_opt(r.p_evidence) << '\t' << format_opt(r.p_base) << '\t' << format_opt(r.U) << '\t'
      << format_opt(r.D) << '\t' << format_opt(r.R) << '\t' << format_opt(r.tau) << '\t' << int(r.reassessed)
      << '\t' << r.hmm_reads << '\t' << r.nucleotide.argmin_offset << '\t' << r.gt << '\t' << format_opt(g.gq)
      << '\t' << format_opt(r.p_final) << '\t' << r.ood << '\t' << decision_name(r.decision) << '\t'
      << r.indel.pos << '\n';
}

int Table::col(const std::string& name) const {
    for (size_t i = 0; i < columns.size(); ++i)
        if (columns[i] == name) return static_cast<int>(i);
    throw std::invalid_argument("table has no column '" + name + "'");
}

Table read_table(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read table " + path);
    Table t;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line.rfind("##", 0) == 0) { t.comments.push_back(line.substr(2)); continue; }
        if (t.columns.empty()) { t.columns = split_tab(line); continue; }
        t.rows.push_back(split_tab(line));
        if (t.rows.back().size() != t.columns.size())
            throw std::runtime_error("table row has wrong number of columns: " + line.substr(0, 80));
    }
    return t;
}

Pattern pattern_from_row(const Table& t, const std::vector<std::string>& row) {
    Pattern p;
    p.type = row[t.col("type")] == "INS" ? IndelType::Insertion : IndelType::Deletion;
    p.length = std::stoi(row[t.col("length")]);
    p.motif = row[t.col("motif")];
    p.context = row[t.col("context")] == "." ? "" : row[t.col("context")];
    p.hp_length = std::stoi(row[t.col("hp_length")]);
    p.hp_base = row[t.col("hp_base")].empty() ? 'N' : row[t.col("hp_base")][0];
    p.period = std::stoi(row[t.col("period")]);
    p.copies = std::stoi(row[t.col("copies")]);
    p.repeat_span = std::stoi(row[t.col("repeat_span")]);
    p.entropy = std::stod(row[t.col("entropy")]);
    return p;
}

}  // namespace nindels
