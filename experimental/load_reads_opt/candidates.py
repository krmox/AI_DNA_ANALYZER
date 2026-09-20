"""load_reads candidates; each returns (reads, labels) with the contract of read_level_pileup.load_reads."""
import numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import read_level_pileup as R
from read_level_pileup import ReadLevelPileupProvider, READ_DIM, _PAD_CODE

def baseline(d, region):   # C0: original pure-Python path + original O(#intervals) BED scan
    from common import use_original_window_scan
    use_original_window_scan(); R.NATIVE_READS_AVAILABLE = False
    return R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"])

class _MemoProvider(ReadLevelPileupProvider):   # C1: pure-Python micro-optimisation
    def _locus_reads(self, column):
        rows, keys, memo = [], [], self._memo
        for read in column.pileups:
            row = self._read_row(read)
            if row is not None:
                rows.append(row); name = read.alignment.query_name
                k = memo.get(name)
                if k is None: k = memo[name] = R.read_sort_key(name)
                keys.append(k)
        matrix = self._empty_locus()
        if not rows: return matrix
        order = np.argsort(np.asarray(keys, dtype=np.uint64), kind="stable")
        sel = order[: self.max_reads]
        matrix[: sel.size] = np.asarray(rows, dtype=np.uint8)[sel]
        return matrix

def microopt(d, region):   # one pileup pass per contiguous run of windows + per-name hash memo (NOT window-exact; timing only)
    from common import use_original_window_scan
    use_original_window_scan(); R.NATIVE_READS_AVAILABLE = False
    prov = _MemoProvider(fasta_path=d["fasta"], bam_path=d["bam"], vcf_path=d["vcf"], contig=d["contig"], region=region,
                         seq_len=64, high_confidence_bed=d["bed"], max_reads=48)
    prov._memo = {}; reads, labels = [], []
    with prov:
        windows = prov._build_windows(); prov._windows = windows; i = 0
        while i < len(windows):
            j = i
            while j + 1 < len(windows) and windows[j + 1].start == windows[j].end: j += 1
            s, e = windows[i].start, windows[j].end
            block = np.zeros((e - s, prov.max_reads, READ_DIM), dtype=np.uint8); block[:, :, 0] = _PAD_CODE
            for col in prov._bam.pileup(prov._bam_contig, s, e, truncate=True, min_base_quality=0, stepper="samtools"):
                block[col.reference_pos - s] = prov._locus_reads(col)
            for w in windows[i:j + 1]:
                reads.append(block[w.start - s: w.end - s]); labels.append(np.asarray(prov._fetch_labels(w), dtype=np.int64))
            i = j + 1; prov._memo.clear()
    return np.concatenate(reads), np.concatenate(labels)

def native(d, region):
    R.NATIVE_READS_AVAILABLE = True
    return R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"])

def native_origscan(d, region):   # native C + ORIGINAL BED scan: isolates the bisect fix
    from common import use_original_window_scan
    use_original_window_scan(); R.NATIVE_READS_AVAILABLE = True
    return R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"])

def chunks(region, chunk_bp):
    s, e = region
    return [(a, min(a + chunk_bp, e)) for a in range(s, e, chunk_bp)]
def run_serial(fn, d, region, chunk_bp): return [fn(d, c) for c in chunks(region, chunk_bp)]
def run_threads(fn, d, region, chunk_bp, workers):
    with ThreadPoolExecutor(workers) as ex: return list(ex.map(lambda c: fn(d, c), chunks(region, chunk_bp)))
def _proc_task(args):
    fn_name, d, c = args
    import candidates
    return getattr(candidates, fn_name)(d, c)
def run_procs(fn_name, d, region, chunk_bp, workers):
    import multiprocessing as mp
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as ex:
        return list(ex.map(_proc_task, [(fn_name, d, c) for c in chunks(region, chunk_bp)]))
