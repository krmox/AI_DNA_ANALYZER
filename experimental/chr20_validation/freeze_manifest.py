"""Freeze record for the chr20 holdout validation (written BEFORE any chr20 call is made).
usage: python3 freeze_manifest.py OUT.json [--analysis FILE ...]"""
import hashlib, json, os, subprocess, sys, platform
import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
W = _ROOT
M = _ROOT
os.chdir(W); sys.path.insert(0, W)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()

MODEL_FILES = ["cascade.py", "experimental/genotype_layer/method_c_regression.py", "quality_error_model.py",
               "binomial_baseline.py", "robustness_benchmark.py", "bench_v12_stage1.py", "bench_v13_robustness.py",
               "bench_v14_crosschrom.py", "bench_v19_unseen.py", "config.py", "locus_evidence.py", "train_raw_pileup.py",
               "native/pileup_native.c"]
OPTIMISED_FILES = ["read_level_pileup.py", "providers.py", "pileup_counts.py", "extract_bench_v12.py",
                   "native/reads_native.c", "native/build.sh"]
out = {"purpose": "chr20 holdout validation freeze (rebuilt after /tmp wipe on 2026-09-20)",
       "created_local": subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip()}
out["git_base_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
out["git_branch"] = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True).stdout.strip()
diff = subprocess.run(["git", "diff", "--", "*.py", "*.sh", "*.c"], capture_output=True, text=True).stdout
out["working_tree_diff_sha256_tracked_py_sh_c"] = hashlib.sha256(diff.encode()).hexdigest()
out["model_files_sha256"] = {f: sha(f) for f in MODEL_FILES if os.path.exists(f)}
out["optimised_files_sha256"] = {f: sha(f) for f in OPTIMISED_FILES if os.path.exists(f)}
import cascade
out["frozen_constants_as_imported"] = {
    "FROZEN_BINOMIAL_THRESHOLD": cascade.FROZEN_BINOMIAL_THRESHOLD, "FROZEN_PB_THRESHOLD": cascade.FROZEN_PB_THRESHOLD,
    "FROZEN_ROUTER_CUTOFF": cascade.FROZEN_ROUTER_CUTOFF,
    "match_expected": (cascade.FROZEN_BINOMIAL_THRESHOLD, cascade.FROZEN_PB_THRESHOLD, cascade.FROZEN_ROUTER_CUTOFF) == (7.0, 10.5, 5.411872376933351)}
out["method_c"] = {"eps": 0.01, "gq_cap": 99.0, "zero_zero_forced_to_0/1_gq_cap": 5.0}
D = f"{M}/data/giab_hg002_chr20"
data = {"chr20_full.bam": f"{D}/chr20_full.bam", "chr20_15x.bam": f"{D}/chr20_15x.bam",
        "chr20.vcf.gz (truth)": f"{D}/chr20.vcf.gz", "chr20_highconf.bed": f"{D}/chr20_highconf.bed",
        "reference chr20_full.fa": f"{M}/data/reference/chr20_full.fa"}
out["data_sha256"] = {k: sha(v) for k, v in data.items() if os.path.exists(v)}
out["data_sources"] = {
    "bam": "GIAB HG002_NA24385_son NIST_Illumina_2x250bps novoalign HG002.GRCh38.2x250.bam (region chr20)",
    "truth_vcf": "GIAB NISTv4.2.1 GRCh38 HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz (chr20)",
    "bed": "GIAB NISTv4.2.1 GRCh38 HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed (chr20 rows)",
    "downsample": "samtools view -s 42.2114 (seed 42; 15 / 70.9603 native mean depth over the chr20 high-confidence BED); achieved 14.99x"}
out["environment"] = {"python": platform.python_version(), "platform": platform.platform(),
    "cpu": subprocess.run("lscpu | grep 'Model name'", shell=True, capture_output=True, text=True).stdout.strip(),
    "samtools": subprocess.run("samtools --version | head -1", shell=True, capture_output=True, text=True).stdout.strip(),
    "hap_py_image": "quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0 @sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0"}
for mod in ("numpy", "scipy", "pysam", "psutil"):
    out["environment"][mod] = __import__(mod).__version__
if "--analysis" in sys.argv:
    out["analysis_scripts_sha256"] = {f: sha(f) for f in sys.argv[sys.argv.index("--analysis") + 1:] if os.path.exists(f)}
json.dump(out, open(sys.argv[1], "w"), indent=2)
print(json.dumps(out["frozen_constants_as_imported"]))
