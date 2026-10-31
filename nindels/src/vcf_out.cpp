#include "nindels/vcf_out.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>

#include "nindels/table_io.hpp"

namespace nindels {

void write_vcf_header(std::ostream& o, const VcfHeaderInfo& h) {
    o << "##fileformat=VCFv4.2\n";
    o << "##source=nindels-call\n";
    for (const auto& [name, len] : h.contigs) o << "##contig=<ID=" << name << ",length=" << len << ">\n";
    o << "##FILTER=<ID=PASS,Description=\"P_final >= tau_call\">\n"
         "##FILTER=<ID=LowProb,Description=\"P_final < tau_call\">\n"
         "##FILTER=<ID=UNCERTAIN,Description=\"Insufficient evidence or outside the model's validated range\">\n"
         "##FILTER=<ID=NoModel,Description=\"No final model supplied; features only\">\n"
         "##INFO=<ID=P,Number=1,Type=Float,Description=\"P_final (spec 12)\">\n"
         "##INFO=<ID=F,Number=1,Type=Float,Description=\"f = n_alt/(n_alt+n_ref)\">\n"
         "##INFO=<ID=NALT,Number=1,Type=Integer,Description=\"Reads supporting the indel\">\n"
         "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\"Reads supporting the reference\">\n"
         "##INFO=<ID=NOTH,Number=1,Type=Integer,Description=\"Reads carrying a different indel here\">\n"
         "##INFO=<ID=PCTX,Number=1,Type=Float,Description=\"P_context from pattern memory\">\n"
         "##INFO=<ID=NEFFK,Number=1,Type=Float,Description=\"N_effective, spec 7.1 (Kish)\">\n"
         "##INFO=<ID=NEFFL,Number=1,Type=Float,Description=\"N_effective over distinct training loci\">\n"
         "##INFO=<ID=SMAX,Number=1,Type=Float,Description=\"Highest similarity to a stored class\">\n"
         "##INFO=<ID=K,Number=1,Type=Integer,Description=\"Alleles genotyped jointly at the locus (--genotyper multi)\">\n"
         "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads used for joint genotyping\">\n"
         "##INFO=<ID=CLS,Number=.,Type=String,Description=\"Closest classes id:S (up to 3)\">\n"
         "##INFO=<ID=CR,Number=1,Type=Float,Description=\"Region complexity C_r\">\n"
         "##INFO=<ID=RPRI,Number=1,Type=Float,Description=\"Reassessment priority R_i\">\n"
         "##INFO=<ID=TAUR,Number=1,Type=Float,Description=\"tau_review(A)\">\n"
         "##INFO=<ID=REASSESSED,Number=0,Type=Flag,Description=\"Local verification was run\">\n"
         "##INFO=<ID=DLO,Number=1,Type=Float,Description=\"Original Delta_ell = ell(H_alt) - ell(H_ref)\">\n"
         "##INFO=<ID=HQ,Number=1,Type=Float,Description=\"Diploid ln LR non-ref vs ref (PairHMM)\">\n"
         "##INFO=<ID=QMIN,Number=1,Type=Float,Description=\"Minimum nucleotide-level Q_k\">\n"
         "##INFO=<ID=QMINPOS,Number=1,Type=Integer,Description=\"Offset of QMIN from the event\">\n"
         "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
         "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
         "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth per allele\">\n"
         "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Informative reads\">\n";
    for (const auto& e : h.extra) o << "##" << e << '\n';
    o << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << h.sample << '\n';
}

static void info(std::ostream& o, bool& first, const char* key, Opt v) {
    if (!v) return;
    o << (first ? "" : ";") << key << '=' << format_opt(*v);
    first = false;
}

void write_vcf_record(std::ostream& o, const std::string& chrom, const RefSeq& ref, const CandidateResult& r) {
    const VcfAlleles a = vcf_alleles(ref, r.indel);
    const GenotypeLik& g = r.hmm_genotype ? *r.hmm_genotype : r.count_genotype;
    // QUAL: phred-scaled 1 - P_final when a model exists, else the genotype-model LR.
    double qual = r.p_final ? std::min(99.0, -10.0 * std::log10(std::max(1e-10, 1.0 - *r.p_final)))
                            : std::max(0.0, g.qual * 10.0 / std::log(10.0));
    o << chrom << '\t' << a.pos1 << "\t.\t" << a.ref << '\t' << a.alt << '\t' << format_opt(std::round(qual * 100) / 100)
      << '\t' << decision_name(r.decision) << '\t';
    bool first = true;
    info(o, first, "P", r.p_final);
    info(o, first, "F", r.evidence.f);
    info(o, first, "NALT", r.evidence.n_alt);
    info(o, first, "NREF", r.evidence.n_ref);
    info(o, first, "NOTH", r.evidence.n_other);
    if (r.context) {
        info(o, first, "PCTX", r.context->p_context);
        info(o, first, "NEFFK", r.context->n_eff_kish);
        info(o, first, "NEFFL", r.context->n_eff_loci);
        info(o, first, "SMAX", r.context->s_max);
        if (!r.context->closest.empty()) {
            o << (first ? "" : ";") << "CLS=";
            first = false;
            for (size_t i = 0; i < std::min<size_t>(3, r.context->closest.size()); ++i)
                o << (i ? "," : "") << r.context->closest[i].first << ':' << format_opt(r.context->closest[i].second);
        }
    }
    info(o, first, "CR", r.complexity.C);
    info(o, first, "RPRI", r.R);
    info(o, first, "TAUR", r.tau);
    if (r.reassessed) { o << (first ? "" : ";") << "REASSESSED"; first = false; }
    info(o, first, "DLO", r.delta_ell_original);
    if (r.hmm_genotype) info(o, first, "HQ", r.hmm_genotype->qual);
    info(o, first, "QMIN", r.nucleotide.q_min);
    if (r.nucleotide.q_min) info(o, first, "QMINPOS", r.nucleotide.argmin_offset);
    if (first) o << '.';
    const char* gt = r.gt == 2 ? "1/1" : "0/1";
    o << "\tGT:GQ:AD:DP\t" << gt << ':' << static_cast<int>(std::round(g.gq)) << ':' << r.evidence.n_ref << ','
      << r.evidence.n_alt << ':' << r.evidence.n_total << '\n';
}

}  // namespace nindels
