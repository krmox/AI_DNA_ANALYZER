// V2 streaming engine: BAM -> counts -> Binomial -> router -> (PB on routed only | PB everywhere)
// -> Method C -> VCF.  No read tensor is ever materialised: PB consumes each locus' <=48 reads
// straight from the window pass.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace dnav2 {

enum class Mode { CASCADE, PB_ALL };

struct Config {
    std::string bam, fasta, truth, bed, contig = "chr20", sample = "HG002", tag = "chr20_300x", out_dir;
    int64_t contig_len = 64444167;
    int threads = 1;
    int decomp_threads = 0;          // htslib BGZF decompression threads per worker
    int task_windows = 2048;
    Mode mode = Mode::CASCADE;
    std::string loci_file;           // optional: 0-based positions; restrict to their windows + force PB there
    std::string dump_loci;           // per-locus TSV for --loci-file loci (gate)
    std::string dump_llr_prefix;     // full float64 arrays of binomial / pb LLR (all loci, window order)
    std::string model_file;          // V2.x decision model (see model.hpp); empty = V2 control decisions
    std::string dump_features;       // V2.x: per-candidate feature table (TSV)
    int feat_kmin = 3;               // candidate = frame locus with ALT count >= feat_kmin
    bool feat_pb = true;             // also dump PB@{24,48,64,96,128,all}
    bool write_vcf = true;
    int64_t skip_windows = 0, limit_windows = -1;   // benchmark slice of the window list
    bool reference_rows = false;     // use build_rows+select_kept (straight port) instead of the lazy selector
    bool exact_windows = false;      // V1-literal path: separate counts pass + FRESH pileup per 64-bp window
    bool setup_only = false;         // build windows + labels, dump them, stop (frame equivalence test)
    bool quiet = false;
};

struct Summary {
    int64_t n_windows = 0, n_loci = 0, n_frame = 0, n_routed = 0;
    int64_t n_cascade_records = 0, n_pbonly_records = 0, n_cascade_forced = 0, n_pbonly_forced = 0;
    int64_t n_pb_evaluated = 0, n_cand = 0, n_model_calls = 0, n_model_changed = 0, n_pb_model = 0;
    // CPU seconds summed over workers, per stage
    double t_count_pass = 0, t_ref_labels = 0, t_binomial = 0, t_window_rows = 0, t_pb = 0, t_genotype = 0;
    double t_setup = 0, t_wall = 0, t_cpu_user = 0, t_cpu_sys = 0;
    int64_t peak_rss_kb = 0;
    int threads = 1;
};

// returns "" on success
std::string run(const Config& cfg, Summary& out);
std::string summary_json(const Config& cfg, const Summary& s);

}  // namespace dnav2
