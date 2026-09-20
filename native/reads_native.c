/* reads_native: htslib construction of the read-level tensor that
 * read_level_pileup.ReadLevelPileupProvider.window_reads() builds in Python.
 *
 * Equivalence contract (must match read_level_pileup.py exactly):
 *   - pileup: one pass PER WINDOW, == pysam pileup(truncate=True, min_base_quality=0,
 *     stepper="samtools") defaults == mplp with the same read callback as pileup_native.c
 *     (FUNMAP|FSECONDARY|FQCFAIL|FDUP filtered; ignore_orphans), ignore_overlaps
 *     on, max depth 8000. A fresh iterator + mplp per window reproduces pysam's
 *     per-window calls exactly, including htslib's mate-overlap quality adjustment,
 *     whose pairing depends on which reads are alive in the pileup (a longer pass
 *     can differ when read names repeat).
 *   - admission (_read_row): mapq >= min_mapq; not dup/qcfail/unmapped;
 *     not secondary/supplementary.
 *   - row = [base_code, base_quality, mapq, strand, read_position,
 *            end_distance, flags, mismatch_density] (uint8), same arithmetic
 *     order (double) and Python round() == round-half-even (nearbyint).
 *   - per-locus order: stable sort by 64-bit BLAKE2b(digest_size=8, read name)
 *     interpreted big-endian; keep the max_reads smallest.
 *   - unused rows: base_code = 255, all other columns 0.
 * Anything this code cannot reproduce bit-for-bit (l_qseq == 0, non-integer NM
 * tag) raises RuntimeError so the Python wrapper can fall back to the
 * pure-Python path instead of silently diverging.
 */
#include <Python.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "htslib/sam.h"
#include "htslib/hts.h"

/* ---- BLAKE2b (RFC 7693), unkeyed, single-call ------------------------- */
static const uint64_t B2_IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL,
    0xa54ff53a5f1d36f1ULL, 0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
    0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL};
static const uint8_t B2_SIGMA[12][16] = {
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, {14,10,4,8,9,15,13,6,1,12,0,2,11,7,5,3},
    {11,8,12,0,5,2,15,13,10,14,3,6,7,1,9,4}, {7,9,3,1,13,12,11,14,2,6,5,10,4,0,15,8},
    {9,0,5,7,2,4,10,15,14,1,11,12,6,8,3,13}, {2,12,6,10,0,11,8,3,4,13,7,5,15,14,1,9},
    {12,5,1,15,14,13,4,10,0,7,6,3,9,2,8,11}, {13,11,7,14,12,1,3,9,5,0,15,4,8,6,2,10},
    {6,15,14,9,11,3,0,8,12,2,13,7,1,4,10,5}, {10,2,8,4,7,6,1,5,15,11,9,14,3,12,13,0},
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, {14,10,4,8,9,15,13,6,1,12,0,2,11,7,5,3}};
#define ROTR64(x, n) (((x) >> (n)) | ((x) << (64 - (n))))
#define B2G(a, b, c, d, x, y) do { \
    v[a] = v[a] + v[b] + (x); v[d] = ROTR64(v[d] ^ v[a], 32); \
    v[c] = v[c] + v[d];       v[b] = ROTR64(v[b] ^ v[c], 24); \
    v[a] = v[a] + v[b] + (y); v[d] = ROTR64(v[d] ^ v[a], 16); \
    v[c] = v[c] + v[d];       v[b] = ROTR64(v[b] ^ v[c], 63); } while (0)

static void b2_compress(uint64_t h[8], const uint8_t blk[128], uint64_t t, int last) {
    uint64_t v[16], m[16];
    for (int i = 0; i < 16; i++) {
        uint64_t w = 0;
        for (int j = 7; j >= 0; j--) w = (w << 8) | blk[i * 8 + j];
        m[i] = w;
    }
    for (int i = 0; i < 8; i++) { v[i] = h[i]; v[i + 8] = B2_IV[i]; }
    v[12] ^= t;
    if (last) v[14] = ~v[14];
    for (int r = 0; r < 12; r++) {
        const uint8_t *s = B2_SIGMA[r];
        B2G(0, 4, 8, 12, m[s[0]], m[s[1]]);   B2G(1, 5, 9, 13, m[s[2]], m[s[3]]);
        B2G(2, 6, 10, 14, m[s[4]], m[s[5]]);  B2G(3, 7, 11, 15, m[s[6]], m[s[7]]);
        B2G(0, 5, 10, 15, m[s[8]], m[s[9]]);  B2G(1, 6, 11, 12, m[s[10]], m[s[11]]);
        B2G(2, 7, 8, 13, m[s[12]], m[s[13]]); B2G(3, 4, 9, 14, m[s[14]], m[s[15]]);
    }
    for (int i = 0; i < 8; i++) h[i] ^= v[i] ^ v[i + 8];
}

/* blake2b(digest_size=8) of msg, bytes -> big-endian uint64 (== Python
 * int.from_bytes(digest, "big")). */
static uint64_t blake2b_key64(const uint8_t *msg, size_t len) {
    uint64_t h[8];
    for (int i = 0; i < 8; i++) h[i] = B2_IV[i];
    h[0] ^= 0x01010000ULL ^ 8ULL;
    uint8_t blk[128];
    size_t off = 0;
    while (len - off > 128) {
        memcpy(blk, msg + off, 128);
        off += 128;
        b2_compress(h, blk, (uint64_t)off, 0);
    }
    size_t rem = len - off;
    memset(blk, 0, 128);
    if (rem) memcpy(blk, msg + off, rem);
    b2_compress(h, blk, (uint64_t)len, 1);
    uint8_t out[8];
    for (int i = 0; i < 8; i++) out[i] = (uint8_t)(h[0] >> (8 * i));
    uint64_t key = 0;
    for (int i = 0; i < 8; i++) key = (key << 8) | out[i];
    return key;
}

/* ---- pileup plumbing (identical read callback to pileup_native.c) ----- */
typedef struct { samFile *fp; hts_itr_t *iter; } read_iter_state;

static int read_bam_cb(void *data, bam1_t *b) {
    read_iter_state *st = (read_iter_state *)data;
    int ret;
    for (;;) {
        ret = sam_itr_next(st->fp, st->iter, b);
        if (ret < 0) break;
        uint16_t flag = b->core.flag;
        if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if ((flag & BAM_FPAIRED) && !(flag & BAM_FPROPER_PAIR)) continue;
        break;
    }
    return ret;
}

/* Constructor: runs once per read when it enters the pileup; caches the read
 * name hash so it is computed once per read, not once per covered locus. */
static int plp_construct(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    (void)data;
    const char *qn = bam_get_qname(b);
    cd->i = (int64_t)blake2b_key64((const uint8_t *)qn, strlen(qn));
    return 0;
}

typedef struct { uint64_t key; uint32_t idx; uint8_t row[8]; } entry_t;

static int entry_cmp(const void *a, const void *b) {
    const entry_t *x = (const entry_t *)a, *y = (const entry_t *)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    return x->idx < y->idx ? -1 : (x->idx > y->idx);
}

static inline int clip255(long v) { return v > 255 ? 255 : (int)v; }

/* One pileup pass over [wstart, wend), writing rows for the (wend - wstart)
 * loci into out (already pad-initialised). */
static const char *window_pass(samFile *fp, hts_idx_t *idx, int tid, long long wstart, long long wend,
                               int min_mapq, int min_baseq, int max_reads, uint8_t *out,
                               entry_t **pents, size_t *pcap) {
    size_t row_bytes = (size_t)max_reads * 8;
    hts_itr_t *iter = sam_itr_queryi(idx, tid, (hts_pos_t)wstart, (hts_pos_t)wend);
    if (!iter) return "cannot create iterator";
    read_iter_state st = {fp, iter};
    void *data_arr[1] = {&st};
    bam_mplp_t mplp = bam_mplp_init(1, read_bam_cb, data_arr);
    bam_mplp_set_maxcnt(mplp, 8000);
    bam_mplp_init_overlaps(mplp);
    bam_mplp_constructor(mplp, plp_construct);

    int tid_out, pos, n_plp_arr[1];
    const bam_pileup1_t *pil_arr[1];
    const char *err = NULL;
    while (!err && bam_mplp_auto(mplp, &tid_out, &pos, n_plp_arr, pil_arr) > 0) {
        if (pos < wstart || pos >= wend) continue;
        const bam_pileup1_t *pil = pil_arr[0];
        int n_plp = n_plp_arr[0];
        if ((size_t)n_plp > *pcap) {
            size_t ncap = (size_t)n_plp * 2;
            entry_t *tmp = (entry_t *)realloc(*pents, sizeof(entry_t) * ncap);
            if (!tmp) { err = "out of memory"; break; }
            *pents = tmp; *pcap = ncap;
        }
        entry_t *ents = *pents;
        int n = 0;
        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t *p = pil + i;
            bam1_t *b = p->b;
            int mapq = b->core.qual;
            if (mapq < min_mapq) continue;
            uint16_t flag = b->core.flag;
            if (flag & (BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)) continue;
            if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

            int length = b->core.l_qseq;
            if (length <= 0) { err = "l_qseq==0"; break; }

            int flags = 1;
            if (p->indel > 0) flags |= 8; else if (p->indel < 0) flags |= 16;
            int strand = (flag & BAM_FREVERSE) ? 1 : 0;

            long nm = 0;
            uint8_t *aux = bam_aux_get(b, "NM");
            if (aux) {
                char t = (char)*aux;
                if (t == 'c' || t == 'C' || t == 's' || t == 'S' || t == 'i' || t == 'I')
                    nm = (long)bam_aux2i(aux);
                else { err = "non-integer NM"; break; }
            }
            int md = clip255((long)nearbyint(2550.0 * (double)nm / (double)(length > 1 ? length : 1)));

            entry_t *e = &ents[n];
            e->key = (uint64_t)p->cd.i;
            e->idx = (uint32_t)n;
            uint8_t *r = e->row;
            if (p->is_del || p->is_refskip) {
                r[0] = 4; r[1] = 0; r[2] = (uint8_t)mapq; r[3] = (uint8_t)strand;
                r[4] = 0; r[5] = 0; r[6] = (uint8_t)(flags | 4); r[7] = (uint8_t)md;
            } else {
                int qpos = p->qpos;
                uint8_t *qual = bam_get_qual(b);
                int q = (qual[0] == 0xff) ? 0 : (int)qual[qpos];
                int code = 5;
                switch (seq_nt16_str[bam_seqi(bam_get_seq(b), qpos)]) {
                    case 'A': code = 0; break; case 'C': code = 1; break;
                    case 'G': code = 2; break; case 'T': code = 3; break;
                    default: break;
                }
                if (code <= 3 && q >= min_baseq) flags |= 2;
                int rest = length - 1 - qpos; if (rest < 0) rest = 0;
                int end_dist = qpos < rest ? qpos : rest;
                if (end_dist < 5) flags |= 32;
                double half = (double)length / 2.0; if (half < 1.0) half = 1.0;
                r[0] = (uint8_t)code; r[1] = (uint8_t)q; r[2] = (uint8_t)mapq; r[3] = (uint8_t)strand;
                r[4] = (uint8_t)clip255((long)nearbyint(255.0 * (double)qpos / (double)(length > 1 ? length : 1)));
                r[5] = (uint8_t)clip255((long)nearbyint(255.0 * (double)end_dist / half));
                r[6] = (uint8_t)flags; r[7] = (uint8_t)md;
            }
            n++;
        }
        if (err || n == 0) continue;
        qsort(ents, (size_t)n, sizeof(entry_t), entry_cmp);
        int keep = n < max_reads ? n : max_reads;
        uint8_t *dst = out + (size_t)(pos - wstart) * row_bytes;
        for (int k = 0; k < keep; k++) memcpy(dst + (size_t)k * 8, ents[k].row, 8);
    }
    bam_mplp_destroy(mplp);
    hts_itr_destroy(iter);
    return err;
}

typedef struct {
    const char *bam_path, *contig;
    const int64_t *ws;
    Py_ssize_t n_windows;
    int seq_len, min_mapq, min_baseq, max_reads;
    uint8_t *out;
    size_t win_bytes;
    volatile Py_ssize_t *next;          /* shared work counter */
    pthread_mutex_t *mu;
    const char *err;                    /* per-thread */
    int io_error;                       /* 1 = cannot open/index/contig (IOError) */
} worker_t;

#define WORK_BLOCK 8

static void *worker_main(void *arg) {
    worker_t *w = (worker_t *)arg;
    samFile *fp = sam_open(w->bam_path, "r");
    if (!fp) { w->err = "cannot open BAM"; w->io_error = 1; return NULL; }
    hts_set_cache_size(fp, 8 * 1024 * 1024);   /* adjacent windows re-read the same BGZF blocks */
    sam_hdr_t *hdr = sam_hdr_read(fp);
    hts_idx_t *idx = hdr ? sam_index_load(fp, w->bam_path) : NULL;
    int tid = (hdr && idx) ? sam_hdr_name2tid(hdr, w->contig) : -1;
    if (tid < 0) { w->err = "cannot open/index BAM or find contig"; w->io_error = 1; }
    entry_t *ents = (entry_t *)malloc(sizeof(entry_t) * 8192);
    size_t cap = 8192;
    while (!w->err) {
        pthread_mutex_lock(w->mu);
        Py_ssize_t lo = *w->next; *w->next += WORK_BLOCK;
        pthread_mutex_unlock(w->mu);
        if (lo >= w->n_windows) break;
        Py_ssize_t hi = lo + WORK_BLOCK < w->n_windows ? lo + WORK_BLOCK : w->n_windows;
        for (Py_ssize_t i = lo; i < hi && !w->err; i++)
            w->err = window_pass(fp, idx, tid, w->ws[i], w->ws[i] + w->seq_len, w->min_mapq,
                                 w->min_baseq, w->max_reads, w->out + (size_t)i * w->win_bytes, &ents, &cap);
    }
    free(ents);
    if (idx) hts_idx_destroy(idx);
    if (hdr) sam_hdr_destroy(hdr);
    sam_close(fp);
    return NULL;
}

/* read_windows(bam, contig, starts: int64 buffer, seq_len, min_mapq, min_baseq, max_reads, n_threads)
 *   -> bytes (n_windows * seq_len, max_reads, 8)
 * Windows are independent passes, so any n_threads yields identical bytes. */
static PyObject *read_windows(PyObject *self, PyObject *args) {
    const char *bam_path, *contig;
    Py_buffer starts;
    int seq_len, min_mapq, min_baseq, max_reads, n_threads;
    if (!PyArg_ParseTuple(args, "ssy*iiiii", &bam_path, &contig, &starts, &seq_len,
                          &min_mapq, &min_baseq, &max_reads, &n_threads))
        return NULL;
    Py_ssize_t n_windows = starts.len / (Py_ssize_t)sizeof(int64_t);
    if (seq_len < 1 || max_reads < 1 || n_windows < 1 || n_threads < 1 || n_threads > 64) {
        PyBuffer_Release(&starts);
        return PyErr_Format(PyExc_ValueError, "bad seq_len/max_reads/window list/n_threads");
    }
    size_t row_bytes = (size_t)max_reads * 8;
    size_t win_bytes = (size_t)seq_len * row_bytes;
    size_t total = (size_t)n_windows * win_bytes;
    PyObject *result = PyBytes_FromStringAndSize(NULL, (Py_ssize_t)total);
    if (!result) { PyBuffer_Release(&starts); return NULL; }
    uint8_t *out = (uint8_t *)PyBytes_AS_STRING(result);
    memset(out, 0, total);
    for (size_t i = 0; i < (size_t)n_windows * (size_t)seq_len * (size_t)max_reads; i++) out[i * 8] = 255;

    volatile Py_ssize_t next = 0;
    pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
    worker_t *ws_ = (worker_t *)calloc((size_t)n_threads, sizeof(worker_t));
    pthread_t *th = (pthread_t *)calloc((size_t)n_threads, sizeof(pthread_t));
    int *started = (int *)calloc((size_t)n_threads, sizeof(int));
    for (int t = 0; t < n_threads; t++) {
        worker_t *w = &ws_[t];
        w->bam_path = bam_path; w->contig = contig; w->ws = (const int64_t *)starts.buf;
        w->n_windows = n_windows; w->seq_len = seq_len; w->min_mapq = min_mapq;
        w->min_baseq = min_baseq; w->max_reads = max_reads; w->out = out;
        w->win_bytes = win_bytes; w->next = &next; w->mu = &mu;
    }
    Py_BEGIN_ALLOW_THREADS
    for (int t = 1; t < n_threads; t++)
        started[t] = (pthread_create(&th[t], NULL, worker_main, &ws_[t]) == 0);
    worker_main(&ws_[0]);
    for (int t = 1; t < n_threads; t++) if (started[t]) pthread_join(th[t], NULL);
    Py_END_ALLOW_THREADS

    const char *err = NULL; int io_error = 0;
    for (int t = 0; t < n_threads; t++) if (ws_[t].err && !err) { err = ws_[t].err; io_error = ws_[t].io_error; }
    free(ws_); free(th); free(started);
    PyBuffer_Release(&starts);
    if (err) {
        Py_DECREF(result);
        if (io_error) return PyErr_Format(PyExc_IOError, "%s: %s %s", err, bam_path, contig);
        return PyErr_Format(PyExc_RuntimeError, "reads_native cannot reproduce Python path: %s", err);
    }
    return result;
}

static PyObject *hash_name(PyObject *self, PyObject *args) {
    const char *s; Py_ssize_t n;
    if (!PyArg_ParseTuple(args, "y#", &s, &n)) return NULL;
    return PyLong_FromUnsignedLongLong(blake2b_key64((const uint8_t *)s, (size_t)n));
}

static PyMethodDef Methods[] = {
    {"read_windows", read_windows, METH_VARARGS,
     "read_windows(bam, contig, starts_int64_buffer, seq_len, min_mapq, min_baseq, max_reads, n_threads) -> bytes\n"
     "uint8 buffer (n_windows*seq_len, max_reads, 8): the read-level tensor, one pileup pass per window."},
    {"hash_name", hash_name, METH_VARARGS, "blake2b-8 big-endian key of a bytes name (test hook)"},
    {NULL, NULL, 0, NULL}};
static struct PyModuleDef mod = {PyModuleDef_HEAD_INIT, "reads_native", NULL, -1, Methods};
PyMODINIT_FUNC PyInit_reads_native(void) { return PyModule_Create(&mod); }
