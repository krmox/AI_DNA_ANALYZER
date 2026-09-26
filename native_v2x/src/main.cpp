#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include "dnav2/engine.hpp"
#include "dnav2/extract.hpp"

static void usage() {
    std::fprintf(stderr,
        "dnav2 run --bam B --ref FASTA --truth TRUTH.vcf.gz --bed BED --out-dir D [options]\n"
        "  --contig chr20 --contig-len 64444167 --sample HG002 --tag chr20_300x\n"
        "  --mode cascade|pb_all   (cascade: PB only on routed loci; pb_all: PB everywhere + PB-only VCF)\n"
        "  --threads N  --decomp-threads N  --task-windows N\n"
        "  --loci-file F --dump-loci OUT.tsv   (gate: restrict to windows of these 0-based loci, dump per-locus fields)\n"
        "  --exact-windows   (V1-literal: counts pass + fresh pileup per window; default is one fused streaming pileup)\n  --skip-windows M --limit-windows N   (benchmark slice of the window list)\n  --dump-llr-prefix P   (write P.binom.f64 / P.pb.f64 for every scored locus)\n");
}

int main(int argc, char** argv) {
    if (argc < 2 || std::strcmp(argv[1], "run") != 0) { usage(); return 2; }
    dnav2::Config c;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        auto val = [&]() -> std::string { if (i + 1 >= argc) { usage(); std::exit(2); } return argv[++i]; };
        if (a == "--bam") c.bam = val(); else if (a == "--ref") c.fasta = val();
        else if (a == "--truth") c.truth = val(); else if (a == "--bed") c.bed = val();
        else if (a == "--out-dir") c.out_dir = val(); else if (a == "--contig") c.contig = val();
        else if (a == "--contig-len") c.contig_len = std::atoll(val().c_str());
        else if (a == "--sample") c.sample = val(); else if (a == "--tag") c.tag = val();
        else if (a == "--threads") c.threads = std::atoi(val().c_str());
        else if (a == "--decomp-threads") c.decomp_threads = std::atoi(val().c_str());
        else if (a == "--task-windows") c.task_windows = std::atoi(val().c_str());
        else if (a == "--mode") { std::string m = val(); c.mode = m == "pb_all" ? dnav2::Mode::PB_ALL : dnav2::Mode::CASCADE; }
        else if (a == "--loci-file") c.loci_file = val(); else if (a == "--dump-loci") c.dump_loci = val();
        else if (a == "--dump-llr-prefix") c.dump_llr_prefix = val();
        else if (a == "--model") c.model_file = val();
        else if (a == "--dump-features") c.dump_features = val();
        else if (a == "--feat-kmin") c.feat_kmin = std::atoi(val().c_str());
        else if (a == "--no-feat-pb") c.feat_pb = false;
        else if (a == "--setup-only") c.setup_only = true;
        else if (a == "--exact-windows") c.exact_windows = true;
        else if (a == "--reference-rows") c.reference_rows = true;
        else if (a == "--skip-windows") c.skip_windows = std::atoll(val().c_str());
        else if (a == "--limit-windows") c.limit_windows = std::atoll(val().c_str());
        else { usage(); return 2; }
    }
    if (c.bam.empty() || c.fasta.empty() || c.truth.empty() || c.bed.empty() || c.out_dir.empty()) { usage(); return 2; }
    dnav2::Summary s;
    std::string err = dnav2::run(c, s);
    if (!err.empty()) { std::fprintf(stderr, "ERROR: %s\n", err.c_str()); return 1; }
    if (c.setup_only) { std::printf("setup: windows=%lld loci=%lld setup_s=%.2f\n", (long long)s.n_windows, (long long)s.n_loci, s.t_setup); return 0; }
#ifdef DNAV2_PROFILE
    dnav2::prof_dump((c.out_dir + "/profile_" + c.tag + ".tsv").c_str());
#endif
    std::string js = dnav2::summary_json(c, s);
    std::ofstream(c.out_dir + "/summary_" + c.tag + ".json") << js;
    std::fputs(js.c_str(), stdout);
    return 0;
}
