"""Cumulative stage decomposition of baseline load_reads on identical windows (unprofiled wall + CPU). usage: [START END]"""
import time, resource, sys
from common import *
import numpy as np
from read_level_pileup import ReadLevelPileupProvider, read_sort_key
import read_level_pileup as R
use_original_window_scan(); R.NATIVE_READS_AVAILABLE = False
d = HG005; region = tuple(int(x) for x in sys.argv[1:3]) if len(sys.argv) > 2 else (1_000_000, 1_060_000)
prov = ReadLevelPileupProvider(fasta_path=d["fasta"], bam_path=d["bam"], vcf_path=d["vcf"], contig=d["contig"], region=region, seq_len=64, high_confidence_bed=d["bed"], max_reads=48)
def cpu(): r = resource.getrusage(resource.RUSAGE_SELF); return r.ru_utime + r.ru_stime
res = {}
with prov:
    t = time.perf_counter(); c = cpu(); windows = prov._build_windows(); res["build_windows"] = (time.perf_counter() - t, cpu() - c)
    bam, contig = prov._bam, prov._bam_contig; res["n_windows"] = len(windows); res["n_loci"] = 64 * len(windows)
    def run(name, fn):
        t = time.perf_counter(); c = cpu(); n = sum(fn(w) for w in windows); res[name] = (time.perf_counter() - t, cpu() - c, n)
    def cols(w): return bam.pileup(contig, w.start, w.end, truncate=True, min_base_quality=0, stepper="samtools")
    def it(w): return sum(1 for _ in cols(w))
    def it_pl(w):
        n = 0
        for col in cols(w):
            for r in col.pileups: n += 1
        return n
    def it_aln(w):
        n = 0
        for col in cols(w):
            for r in col.pileups:
                a = r.alignment; _ = (a.mapping_quality, a.is_duplicate, a.is_qcfail, a.is_unmapped, a.is_secondary, a.is_supplementary, a.is_reverse); n += 1
        return n
    def it_tag(w):
        n = 0
        for col in cols(w):
            for r in col.pileups:
                a = r.alignment; _ = (a.mapping_quality, a.is_duplicate, a.is_qcfail, a.is_unmapped, a.is_secondary, a.is_supplementary, a.is_reverse)
                _ = r.indel; _ = a.query_length
                try: _ = int(a.get_tag("NM"))
                except (KeyError, ValueError): pass
                n += 1
        return n
    def it_seq(w):
        n = 0
        for col in cols(w):
            for r in col.pileups:
                a = r.alignment; _ = (a.mapping_quality, a.is_duplicate, a.is_qcfail, a.is_unmapped, a.is_secondary, a.is_supplementary, a.is_reverse)
                _ = r.indel; _ = a.query_length
                try: _ = int(a.get_tag("NM"))
                except (KeyError, ValueError): pass
                if r.is_del or r.is_refskip: n += 1; continue
                p = r.query_position
                _ = int(a.query_qualities[p]) if a.query_qualities else 0
                _ = a.query_sequence[p].upper(); n += 1
        return n
    def it_hash(w):
        n = 0
        for col in cols(w):
            for r in col.pileups: read_sort_key(r.alignment.query_name); n += 1
        return n
    for nm, f in [("1_pysam_columns_only", it), ("2_+PileupRead_objects", it_pl), ("3_+alignment_and_flags", it_aln),
                  ("4_+NM_tag_indel_len", it_tag), ("5_+qpos_qual_seq", it_seq), ("6_hash_only(name+blake2b) [alone, incl. PileupRead objects]", it_hash)]:
        run(nm, f)
    t = time.perf_counter(); c = cpu()
    for w in windows: prov.window_reads(w)
    res["7_full_window_reads"] = (time.perf_counter() - t, cpu() - c, 0)
t = time.perf_counter(); c = cpu(); R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"]); res["8_load_reads_e2e"] = (time.perf_counter() - t, cpu() - c, 0)
for k, v in res.items(): print(k, v)
